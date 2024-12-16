#!/usr/bin/env python
# coding: utf-8
import os
import shutil
import math
import numpy as np
import copy
import matplotlib.pyplot as plt




atmNToSymb=\
    {1:'H',2:'He',3:'Li',4:'Be',5:'B',6:'C',7:'N',8:'O',9:'F',10:'Ne',11:'Na',12:'Mg',13:'Al',14:'Si',15:'P',16:'S',
     17:'Cl',18:'Ar',19:'K',20:'Ca',21:'Sc',22:'Ti',23:'V',24:'Cr',25:'Mn',26:'Fe',27:'Co',28:'Ni',29:'Cu',30:'Zn',
     31:'Ga',32:'Ge',33:'As',34:'Se',35:'Br',36:'Kr',37:'Rb',38:'Sr',39:'Y',40:'Zr',41:'Nb',42:'Mo',43:'Tc',44:'Ru',
     45:'Rh',46:'Pd',47:'Ag',48:'Cd',49:'In',50:'Sn',51:'Sb',52:'Te',53:'I',54:'Xe',55:'Cs',56:'Ba',57:'La',58:'Ce',
     59:'Pr',60:'Nd',61:'Pm',62:'Sm',63:'Eu',64:'Gd',65:'Tb',66:'Dy',67:'Ho',68:'Er',69:'Tm',70:'Yb',71:'Lu',72:'Hf',
     73:'Ta',74:'W',75:'Re',76:'Os',77:'Ir',78:'Pt',79:'Au',80:'Hg',81:'Tl',82:'Pb',83:'Bi',84:'Po',85:'At',86:'Rn',
     87:'Fr',88:'Ra',89:'Ac',90:'Th',91:'Pa',92:'U',93:'Np',94:'Pu',95:'Am',96:'Cm',97:'Bk',98:'Cf',99:'Es',100:'Fm',
     101:'Md',102:'No',103:'Lr',104:'Rf',105:'Db',106:'Sg',107:'Bh',108:'Hs',109:'Mt',110:'Ds',111:'Rg',112:'Cn',
     113:'Nh',114:'Fl',115:'Mc',116:'Lv',117:'Ts',118:'Og'}

atmSymbToN=\
    {'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,
     'Cl':17,'Ar':18,'K':19,'Ca':20,'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,
     'Ga':31,'Ge':32,'As':33,'Se':34,'Br':35,'Kr':36,'Rb':37,'Sr':38,'Y':39,'Zr':40,'Nb':41,'Mo':42,'Tc':43,'Ru':44,
     'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,'Sn':50,'Sb':51,'Te':52,'I':53,'Xe':54,'Cs':55,'Ba':56,'La':57,'Ce':58,
     'Pr':59,'Nd':60,'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,'Dy':66,'Ho':67,'Er':68,'Tm':69,'Yb':70,'Lu':71,'Hf':72,
     'Ta':73,'W':74,'Re':75,'Os':76,'Ir':77,'Pt':78,'Au':79,'Hg':80,'Tl':81,'Pb':82,'Bi':83,'Po':84,'At':85,'Rn':86,
     'Fr':87,'Ra':88,'Ac':89,'Th':90,'Pa':91,'U':92,'Np':93,'Pu':94,'Am':95,'Cm':96,'Bk':97,'Cf':98,'Es':99,'Fm':100,
     'Md':101,'No':102,'Lr':103,'Rf':104,'Db':105,'Sg':106,'Bh':107,'Hs':108,'Mt':109,'Ds':110,'Rg':111,'Cn':112,
     'Nh':113,'Fl':114,'Mc':115,'Lv':116,'Ts':117,'Og':118}

sqm_head = """Run semi-empirical minimization
 &qmmm
    qm_theory='DFTB3', grms_tol=0.0005,
 scfconv=1.d-10, ndiis_attempts=700, qmcharge=0, maxcyc=0,
 /\n"""

#sqm_head = \
#"""Run semi-empirical minimization
# &qmmm
#    qm_theory='DFTB3', grms_tol=0.0005,
# scfconv=1.d-10, ndiis_attempts=700,   qmcharge=0, dftb_slko_path=\'PATH\', maxcyc=0,
# /\n"""


def get_gaussian_orig(mu: float, sigma: float, x: float) -> float:
    return (1 / (np.sqrt(2 * np.pi * (sigma ** 2)))) * np.exp(- (x - mu) ** 2 / (2 * sigma ** 2))


def get_gaussian(w: float, x: float, x1: float, x2: float, sigma: float = None) -> float:
    if abs(x1 - x2) <= 0.005:
        return 0.0

    # Determino los límites del intervalo
    lb = min(x1, x2)
    hb = max(x1, x2)

    if not sigma:
        mu = (hb + lb) / 2
        sigma = (mu - lb)*0.3

    return (w * (np.exp((x - lb) * (hb - x) / ((sigma) ** 2))))


def get_bounded_harmonic(x: float, x1: float, x2: float, top_val: float = None) -> float:
    lb = 0
    hb = 0
    if x1 > x2:
        lb = x2
        hb = x1
    elif x1 < x2:
        lb = x1
        hb = x2
    elif abs(x1 - x2) <= 0.05:
        return 0.0

    if x > hb or x < lb:
        return 0.0

    m = 2 * math.pi / (hb - lb)
    n = -((hb + 3 * lb) * math.pi / (2 * (hb - lb)))
    return 0.5 * top_val * (math.sin(m * x + n) + 1.0)


class Atom:

    def __init__(self, elem: str = 'X', coord: list = []):
        self.__elem = elem
        self.__coord = np.array(coord)
        if elem in atmSymbToN:
            self.__atmN = atmSymbToN[elem]
        else:
            self.__atmN = 0

    def get_element(self) -> str:
        return self.__elem

    def get_coord(self) -> list:
        return self.__coord

    def get_atomic_number(self) -> int:
        return self.__atmN

    def set_element(self, elem: str) -> None:
        if elem in atmSymbToN:
            self.__elem = elem
            self.__atmN = atmSymbToN[elem]
        else:
            print("The input element: " + str(elem) + "  is not a valid input")

    def get_distance(self, atm) -> float:
        x, y, z = atm.get_coord()
        return math.sqrt(pow(self.__coord[0] - x, 2) + pow(self.__coord[1] - y, 2) + pow(self.__coord[2] - z, 2))

    def match_elem(self, inp) -> bool:
        if self.__elem == inp:
            return True
        else:
            return False


class XYZHandler:
    """Clase para leer un fichero XYZ y obtener:
    1) las distribuciones de distancias para un par de elementos
    2) dada una función de términos repulsivos en splines para un par de elementos, obtener el total de esta energía por
     cada frame
    """

    def __init__(self, xyz_path: str = None):
        self.path = xyz_path
        self.frames = {}
        self.elemToIndex = {}
        self.indexToElem = {}
        self.read_xyz(file=xyz_path)
        self.verify_frames()
        self.pdb_names = {}
        self.set_pdb_names()

    def set_pdb_names(self) -> None:
        for e in self.elemToIndex:
            for i in range(0, len(self.elemToIndex[e])):
                self.pdb_names[self.elemToIndex[e][i]] = e + str(i)

    def get_pdb_name(self, indx: int) -> str:
        return self.pdb_names[indx]

    def verify_frames(self):
        for frame in self.frames:
            for id_ in self.frames[frame]:
                if self.frames[frame][id_].get_element() != self.indexToElem[id_]:
                    print("Error! el símbolo químico de un átomo no coincide")
                    print("Frame: " + str(frame) + " Atom index: " + str(id_))

    def read_xyz(self, file: str = None) -> None:

        if not file:
            return

        xyz_file = open(file, 'r')
        lines = xyz_file.readlines()
        xyz_file.close()
        numb = int(lines[0].strip())
        frame = {}
        findex = 0
        aindex = 1
        read_struct = True

        jump = False

        for i in range(2, len(lines)):
            line = lines[i]

            if len(line.split()) == 4:
                if read_struct:
                    self.indexToElem[aindex] = line.split()[0]
                    if not line.split()[0] in self.elemToIndex:
                        self.elemToIndex[line.split()[0]] = []

                    self.elemToIndex[line.split()[0]].append(aindex)

                atom = Atom(elem=line.split()[0],
                            coord=[float(line.split()[1]), float(line.split()[2]), float(line.split()[3])])

                frame[aindex] = atom
                aindex += 1

            else:
                if not jump:
                    findex += 1
                    if findex == 1:
                        read_struct = False
                    if findex > 0:
                        self.frames[findex] = copy.deepcopy(frame)

                    aindex = 1
                    frame = {}
                    jump = True

                else:
                    jump = False

        if len(frame) != 0:
            findex += 1
            self.frames[findex] = copy.deepcopy(frame)


    def get_distribution(self, elem1: str, elem2: str, frame: int) -> list:
        dlist = []
        if elem1 == elem2:
            list1 = sorted(self.elemToIndex[elem1])
            for a in range(0, len(list1)):
                for b in range(a+1, len(list1)):
                    dlist.append(self.frames[frame][list1[a]].get_distance(self.frames[frame][list1[b]]))

        else:
            list1 = sorted(self.elemToIndex[elem1])
            list2 = sorted(self.elemToIndex[elem2])
            for i in list1:
                for j in list2:
                    dlist.append(self.frames[frame][i].get_distance(self.frames[frame][j]))
        return dlist

    def write_xyz(self, file: str = None, append=False, frames: list = []) -> None:
        if not file:
            return

        xyz_file = None
        if append:
            xyz_file = open(file, 'a')
        else:
            xyz_file = open(file, 'w')

        if len(frames) == 0:
            frames = sorted(self.frames.keys())

        for frm in frames:
            line = ''
            if frm in self.frames:
                natm = len(self.frames[frm])
                xyz_file.write(str(natm)+'\n\n')
                for a_id in sorted(self.frames[frm].keys()):
                    line = self.frames[frm][a_id].get_element() + '        ' +                            "{:.6f}".format(self.frames[frm][a_id].get_coord()[0]) + '       ' +                           "{:.6f}".format(self.frames[frm][a_id].get_coord()[1]) + '       ' +                            "{:.6f}".format(self.frames[frm][a_id].get_coord()[2]) + '\n'
                    xyz_file.write(line)
        xyz_file.close()

    def write_sqm_in(self, base_path: str = '', frames: list = [], slko_path: str = 'slko') -> None:
        if len(frames) == 0:
            frames = sorted(self.frames.keys())

        for frm in frames:
            sqm_file = open(base_path+str(frm)+'.in', 'w')
            line = ''
            if frm in self.frames:
                sqm_file.write(sqm_head)
                for a_id in sorted(self.frames[frm].keys()):
                    line = '   ' + str(self.frames[frm][a_id].get_atomic_number()) + '   ' + self.pdb_names[a_id] \
                           + '        ' + "{:.3f}".format(self.frames[frm][a_id].get_coord()[0]) + '       ' \
                           + "{:.3f}".format(self.frames[frm][a_id].get_coord()[1]) + '       ' \
                           + "{:.3f}".format(self.frames[frm][a_id].get_coord()[2]) + '\n'

                    sqm_file.write(line)
                sqm_file.close()

class PotPlotRecipe():
    def __init__(self):
        self.ylabel = 'Potential function (EV)'
        self.xlabel = 'Distance (angstroms)'
        self.title = 'Potential function vs distance'


class DFTBParam:

    def __init__(self, ener: str='EV'):
        self.HtoKJ: float = 2625.5
        self.HtoKCal: float = 627.5
        self.HtoEV: float = 27.211

        self.plotter = PotPlotRecipe()
        self.plotter.xlabel = 'Distance (angstroms)'

        if ener == 'EV':
            self.__efactor = self.HtoEV
            self.plotter.ylabel = 'Potential Function (EV)'
            self.__ener = ener

        elif ener == 'KJ':
            self.__efactor = self.HtoKJ
            self.plotter.ylabel = 'Potential Function (KJ)'
            self.__ener = ener

        elif ener == 'KCal':
            self.__efactor = self.HtoKCal
            self.plotter.ylabel = 'Potential Function (KCal)'
            self.__ener = ener

        else:
            self.__efactor = 1
            self.plotter.ylabel = 'Potential Function (Hartree)'
            self.__ener = 'EV'

        self.plotter.title = 'Potential function vs distance'

        self.x = np.array([])
        self.y = np.array([])
        self.mod_y = np.array([])
        self.__elect_lines: list[str] = []
        self.__spline_lines: list[str] = []
        self.__doc_lines: list[str] = []
        self.__spl_dict: dict = {}
        self.__nInt: int = 0
        self.__cutoff: float = -1.0
        self.__a1: float = 0
        self.__a2: float = 0
        self.__a3: float = 0
        self.__cuts: list[str] = []

    def set_potent_curve(self, set_mod_y: bool = False) -> None:
        self.x, self.y = self.compute_spline_potential()
        if set_mod_y:
            self.mod_y = copy.deepcopy(self.y)

    def set_elect_lines(self, elect_lines: list) -> None:
        self.__elect_lines = elect_lines

    def get_elect_lines(self) -> list:
        return self.__elect_lines

    def set_doc_lines(self, doc_lines: list) -> None:
        self.__doc_lines = doc_lines

    def get_doc_lines(self) -> list:
        return self.__doc_lines

    def __set_splines(self, spline_lines: list) -> None:
        self.__spline_lines = spline_lines
        self.__nInt = int(spline_lines[1].split()[0])
        self.__cutoff = float(spline_lines[1].split()[1])
        self.__a1 = float(spline_lines[2].split()[0])
        self.__a2 = float(spline_lines[2].split()[1])
        self.__a3 = float(spline_lines[2].split()[2])

        spl_dict: dict = {}
        cuts: list = []

        for i in range(1, self.__nInt + 1):
            start = spline_lines[2 + i].split()[0]
            spl_dict[start] = [float(x) for x in spline_lines[2 + i].split()[2:]]
            cuts.append(start)

        self.__spl_dict = spl_dict
        self.__cuts = cuts

    def set_splines(self, nInt: int, cutoff: float,
                    a1: float, a2: float, a3: float,
                    cuts: list, spl_dict: dict) -> None:
        if len(spl_dict) != nInt:
            print("Warning! the number of intervals does not match the nInt value, please, introduce"
                  "the spline's info again")
            return

        elif len(spl_dict) != len(cuts):
            print('Warning! the length of the knots\' list does not match the length of the '
                  'dictionary with the splines\' coefficients, please review the data before trying '
                  'to modify the splines')
            return

        for key in cuts:
            if not key in spl_dict:
                print("Warning! the key: {} was not found in the provided dictionary".format(key))
                return

        self.__nInt = nInt
        self.__cutoff = cutoff
        self.__a1 = a1
        self.__a2 = a2
        self.__a3 = a3
        self.__cuts = cuts
        self.__spl_dict = spl_dict

        spline_lines = ["Spline\n", '{}  {}\n'.format(nInt, cutoff),
                        '{}  {}  {}\n'.format('{:.15e}'.format(a1),
                                              '{:.15e}'.format(a2),
                                              '{:.15e}'.format(a3))]
        for i in range(0, len(cuts) - 1):
            tmp = '{}  {}'.format(cuts[i], cuts[i + 1])
            for coef in spl_dict[cuts[i]]:
                tmp = tmp + '  {}'.format('{:.15e}'.format(coef))
            spline_lines.append(tmp + '\n')

        tmp = '{}  {}'.format(cuts[len(cuts) - 1], '{:.5f}'.format(cutoff))
        for coef in spl_dict[cuts[len(cuts) - 1]]:
            tmp = tmp + '  {}'.format('{:.15e}'.format(coef))

        spline_lines.append(tmp + '\n')
        self.__spline_lines = spline_lines

    def set_spline_lines(self, spline_lines: list) -> None:
        self.__spline_lines = spline_lines

    def get_spline_lines(self) -> list:
        return self.__spline_lines

    def set_spl_dict(self, spl: dict) -> None:
        self.__spl_dict = spl

    def get_spl_dict(self) -> dict:
        return self.__spl_dict

    def set_nInt(self, nInt: int) -> None:
        self.__nInt = nInt

    def get_nInt(self) -> int:
        return self.__nInt

    def set_cutoff(self, cutoff: float) -> None:
        self.__cutoff = cutoff

    def get_cutoff(self) -> float:
        return self.__cutoff

    def set_a1(self, a1: float) -> None:
        self.__a1 = a1

    def get_a1(self) -> float:
        return self.__a1

    def set_a2(self, a2: float) -> None:
        self.__a2 = a2

    def get_a2(self) -> float:
        return self.__a2

    def set_a3(self, a3: float) -> None:
        self.__a3 = a3

    def get_a3(self) -> float:
        return self.__a3

    def set_cuts(self, cuts: list) -> None:
        self.__cuts = cuts

    def get_cuts(self) -> list:
        return self.__cuts

    def read_params(self, file_path) -> None:
        pf = open(file_path, 'r')
        lines = pf.readlines()
        pf.close()

        elect_lines: list[str] = []
        spline_lines: list[str] = []
        doc_lines: list[str] = []
        el = True
        sp = False
        dc = False

        for line in lines:
            if not sp:
                if "Spline" == line.strip():
                    el = False
                    sp = True

            else:
                if "<Documentation>" in line:
                    sp = False
                    dc = True

            if el:
                elect_lines.append(line)

            elif sp:
                spline_lines.append(line)

            elif dc:
                doc_lines.append(line)

        self.__elect_lines = elect_lines
        self.__set_splines(spline_lines)
        self.__doc_lines = doc_lines
        self.set_potent_curve(set_mod_y=True)

    def write_params(self, param_path: str):
        param = open(param_path, 'w')
        param.writelines(self.__elect_lines)
        param.writelines(self.__spline_lines)
        param.writelines(self.__doc_lines)
        param.close()

    def plot_spline_potential(self):
        x, y = self.compute_spline_potential()
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.set_title(self.plotter.title)
        ax.set_xlabel(self.plotter.xlabel)
        ax.set_ylabel(self.plotter.ylabel)
        plt.plot(x, y, label="potential", linestyle='dashed', linewidth=1, marker='o',
                          markerfacecolor='blue', markersize=6)
        plt.grid()
        plt.tight_layout()
        plt.minorticks_on()
        plt.legend()
        plt.show()

    def compute_spline_potential(self, step: float = 0.05):
        initial = float(self.__cuts[0])
        final = self.__cutoff
        if step == 0:
            step = (final - initial) / 100
        dist = 0
        dist_list = []
        val_list = []
        while dist <= final:
            dist = dist + step
            dist_list.append(dist)
            val_list.append(self.get_potential(dist))
        return np.array(dist_list), np.array(val_list)

    def get_potential(self, dist: float) -> float:
        initial = float(self.__cuts[0])
        final = self.__cutoff
        if dist < initial:
            return (math.exp(-self.__a1 * dist + self.__a2) + self.__a3) * self.__efactor

        elif dist > final:
            return 0.0

        else:
            val = 0
            for i in range(0, len(self.__cuts) - 1):
                if float(self.__cuts[i]) <= dist < float(self.__cuts[i + 1]):
                    for c in range(0, len(self.__spl_dict[self.__cuts[i]])):
                        coef = self.__spl_dict[self.__cuts[i]][c]
                        val = val + coef * ((dist - float(self.__cuts[i])) ** c)

            if val == 0:
                if dist >= float(self.__cuts[-1]):
                    for c in range(0, len(self.__spl_dict[self.__cuts[-1]])):
                        coef = self.__spl_dict[self.__cuts[-1]][c]
                        val = val + coef * ((dist - float(self.__cuts[-1])) ** c)

            return val * self.__efactor

    def get_potential_plus_gaussian(self, dist: float, x1: float, x2: float, w: float = 0.1,
                                    sigma: float = None) -> float:
        return self.get_potential(dist) + get_gaussian(w=w, x=dist, x1=x1, x2=x2, sigma=sigma) * self.__efactor

    def compute_spline_potential_plus_gaussian(self, x1: float, x2: float, w: float = 0.1, sigma: float = None):
        initial = float(self.__cuts[0])
        final = self.__cutoff
        step = (final - initial) / 100
        dist = 0
        dist_list = []
        val_list = []
        while dist <= final:
            dist = dist + step
            dist_list.append(dist)
            val_list.append(self.get_potential_plus_gaussian(dist=dist, x1=x1, x2=x2, w=w, sigma=sigma))
        return dist_list, val_list

    def plot_spline_potential_plus_gaussian(self, x1: float, x2: float, w: float = 0.1, sigma: float = None):
        x, y = self.compute_spline_potential_plus_gaussian(x1=x1, x2=x2, w=w, sigma=sigma)
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.set_title(self.plotter.title)
        ax.set_xlabel(self.plotter.xlabel)
        ax.set_ylabel(self.plotter.ylabel)
        plt.plot(x, y, label="potential+gaussian", linestyle='dashed', linewidth=1, marker='*',
                          markerfacecolor='red', markersize=6)
        plt.grid()
        plt.tight_layout()
        plt.minorticks_on()
        plt.legend()
        plt.show()

    def sum_gaussian(self, x1: float, x2: float, set_curve: bool = False, w: float = 0.0, sigma: float = None) -> None:
        if set_curve:
            self.set_potent_curve(set_mod_y=True)

        for i in range(0, len(self.x)):
            val = self.x[i]
            self.mod_y[i] = self.mod_y[i] + get_gaussian(x=val, x1=x1, x2=x2, w=w, sigma=sigma) * self.__efactor

    def sum_bounded_gaussian(self, x1: float, x2: float, set_curve: bool = False, top_val: float = None) -> None:
        if not top_val:
            return

        sign = int(top_val/abs(top_val))
        w = sign*self.__efactor*1.00000e-5
        print('w: {}, denom: {}'.format(w, top_val/(w*self.__efactor)))
        sigma = abs(x2-x1)/math.sqrt(4*math.log(top_val/(w*self.__efactor)))
        print('Sigma: {}, w: {}'.format(sigma, w))
        self.sum_gaussian(x1=x1, x2=x2, set_curve=set_curve, w=w, sigma=sigma)

    def sum_bounded_harmonic(self, x1: float, x2: float, set_curve: bool = False, top_val: float = None)\
            -> None:
        if not top_val:
            return

        if set_curve:
            self.set_potent_curve(set_mod_y=True)

        for i in range(0, len(self.x)):
            val = self.x[i]
            self.mod_y[i] = self.mod_y[i] + get_bounded_harmonic(x=val, x1=x1, x2=x2, top_val=top_val)

    def plot_potent_curve(self, xlimits=None, ylimits=None):
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.set_title(self.plotter.title)
        ax.set_xlabel(self.plotter.xlabel)
        ax.set_ylabel(self.plotter.ylabel)

        plt.plot(self.x, self.y, label="potential", linestyle='dashed', linewidth=1, marker='o',
                          markerfacecolor='blue', markersize=6)
        plt.plot(self.x, self.mod_y, label="potential+gaussian", linestyle='dashed', linewidth=1, marker='*',
                          markerfacecolor='red', markersize=6)
        plt.grid()
        plt.tight_layout()
        plt.minorticks_on()
        if xlimits:
            plt.xlim(xlimits[0], xlimits[1])
        if ylimits:
            plt.ylim(ylimits[0], ylimits[1])

        plt.legend()
        plt.show()

    def compute_full_potential(self, dist_list: list) -> float:
        pot = 0.0
        for i in dist_list:
            pot = pot + self.get_potential(i)

        return pot

    def compute_full_potential_in_range(self, dist_list: list, lb: float, hb: float) -> float:
        pot = 0.0
        for i in dist_list:
            if (hb > i) and (lb <= i):
                pot = pot + self.get_potential(i)

        return pot

    def fit_splines(self, index: int = 0) -> None:
        if index >= (len(self.x)-3) or index <= 3:
            print("The provided index is out of permitted range")
            return

        from scipy.optimize import curve_fit

        def exp_funct(x: float, a1: float, a2: float, a3: float) -> float:
            y = np.exp(-a1*x + a2) + a3
            return y

        x1 = self.x[:index+1]
        y1 = self.mod_y[:index+1]/self.__efactor
        fit = curve_fit(exp_funct, x1, y1)
        #print('Setting exponential curve: \nORIGINAL: a1={:.15e}, a2={:.15e}, a3={:.15e}'.
        #      format(self.__a1, self.__a2, self.__a3) + '\n')
        self.set_a1(fit[0][0])
        self.set_a2(fit[0][1])
        self.set_a3(fit[0][2])
        #print('Setting exponential curve: \nORIGINAL: a1={:.15e}, a2={:.15e}, a3={:.15e}'.
        #      format(self.__a1, self.__a2, self.__a3) + '\n')

        der1 = -self.__a1*math.exp(-self.__a1 * self.x[index+1] + self.__a2)

        der2 = self.__spl_dict[self.__cuts[-2]][1]

        x2 = self.x[index:-2]
        y2 = self.mod_y[index:-2]/self.__efactor

        new_cuts = ['{:.5f}'.format(i) for i in self.x[index:-2]]

        fifth_spl = [self.__cuts[-1], self.__spl_dict[self.__cuts[-1]]]
        new_spline_dict = {}

        from scipy.interpolate import CubicSpline
        cs = CubicSpline(x2, y2, bc_type=((1, der1), (1, der2)))

        for i in range(0, len(cs.c[0])):
            new_spline_dict[new_cuts[i]] = [cs.c[3][i], cs.c[2][i], cs.c[1][i], cs.c[0][i]]
            #print('X: ' + new_cuts[i] + ' COEF: c0={:.15e}, c1={:.15e}, c2={:.15e}, c3={:.15e}'.
            #      format(cs.c[3][i], cs.c[2][i], cs.c[1][i], cs.c[0][i]) + '\n')

        new_spline_dict[fifth_spl[0]] = fifth_spl[1]

        self.set_splines(nInt=len(new_spline_dict), spl_dict=new_spline_dict, cutoff=self.__cutoff, cuts=new_cuts,
                         a1=fit[0][0], a2=fit[0][1], a3=fit[0][2])


class ParamHandler:

    def __init__(self, orig_param_file: str = None, traj_file: str = None):

        self.orig_param = DFTBParam(ener='KCal')
        self.mod_param = DFTBParam(ener='KCal')

        if orig_param_file:
            self.orig_param.read_params(orig_param_file)
            self.mod_param.read_params(orig_param_file)
            self.orig_param_file = orig_param_file
        else:
            self.orig_param_file = None

        self.trajectory = XYZHandler(traj_file)

    def sum_bounded_gaussian(self, x1: float = None, x2: float = None, set_curve: bool = True, top_val: int = 10,
                             index: int = 20):
        self.mod_param.sum_bounded_gaussian(x1=x1, x2=x2, set_curve=set_curve, top_val=top_val)
        self.mod_param.fit_splines(index=index)

    def sum_bounded_harmonic(self, x1: float, x2: float, set_curve: bool = False, top_val: float = None,
                             index: int = 20):
        self.mod_param.sum_bounded_harmonic(x1=x1, x2=x2, set_curve=set_curve, top_val=top_val)
        self.mod_param.fit_splines(index=index)

    def reset_mod_params(self):
        if self.orig_param_file:
            self.mod_param = DFTBParam(ener='KCal')
            self.mod_param.read_params(self.orig_param_file)

    def modify_params_to(self, param_file: str = None, output_param_file: str = None):
        aux_param = DFTBParam(ener='KCal')
        aux_param.read_params(param_file)
        aux_param.set_splines(nInt=self.mod_param.get_nInt(), cutoff=self.mod_param.get_cutoff(),
                              a1=self.mod_param.get_a1(), a2=self.mod_param.get_a2(), a3=self.mod_param.get_a3(),
                              cuts=self.mod_param.get_cuts(), spl_dict=self.mod_param.get_spl_dict())
        aux_param.write_params(param_path=output_param_file)

    def get_potential_from_dist(self, elem1='C', elem2='H'):
        dat = []
        for i in range(1, len(self.trajectory.frames)+1):
            dat.append(self.orig_param.compute_full_potential(self.trajectory.get_distribution(elem1=elem1, elem2=elem2,
                                                                                               frame=i)))
        return np.arange(1, len(dat)+1), np.array(dat)

    def get_mod_potential_from_dist(self, elem1='C', elem2='H'):
        dat = []
        for i in range(1, len(self.trajectory.frames) + 1):
            dat.append(self.mod_param.compute_full_potential(self.trajectory.get_distribution(elem1=elem1, elem2=elem2,
                                                                                               frame=i)))
        return np.arange(1, len(dat)+1), np.array(dat)

    def run_sqm(self, sqm_path: str = None, amber_home: str = None, run_path: str = None):
        if amber_home:
            os.environ["AMBERHOME"] = amber_home

        import shutil
        if run_path:
            if os.path.exists(run_path):
                shutil.rmtree(run_path)
                os.makedirs(run_path)
            else:
                os.makedirs(run_path)
        else:
            run_path = './run_path'
            if os.path.exists(run_path):
                shutil.rmtree(run_path)
                os.makedirs(run_path)
            else:
                os.makedirs(run_path)

        self.trajectory.write_sqm_in(base_path=run_path + '/chunks', slko_path='')
        cmds = []
        for item in os.listdir(run_path):
            if item.endswith('.in'):
                cmd = [sqm_path, '-O', '-i', run_path + '/' + item, '-o',
                       run_path + '/' + item.split('.')[0] + '.dat']
                cmds.append(cmd)

        outputs_list = [run_path + "/chunks" + str(i) + ".dat" for i in range(1, len(cmds)+1)]
        from subprocess import Popen, PIPE
        nCPUs = os.cpu_count()*2
        procs = []
        for i in range(len(cmds)):
            procs.append(Popen(cmds[i], stdout=PIPE, stderr=PIPE))
            if i % (nCPUs - 1) == 0:
                for proc in procs:
                    proc.wait()
                procs = []

            elif i == (len(cmds) - 1):
                for proc in procs:
                    proc.wait()

        energ = scf_line_extract(outputs_list)
        x, y = zip(*energ)
        y1 = [y[i] - min(y) for i in range(len(y))]
        return x, y1


def scf_line_extract(outputs_list: list) -> list:
    temp1 = []
    for item in outputs_list:
        # print(item)
        with open(item, 'r') as temp2:
            for line in temp2.readlines():
                # print('Line: ' + line)
                if 'Total SCF energy' in line:

                    temp1.append((int(item.split('/')[-1].split('s')[1].split('.')[0]), float(line.split()[4])))
                else:
                    pass
    return temp1