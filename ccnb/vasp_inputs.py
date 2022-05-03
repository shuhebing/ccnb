from pymatgen.io.vasp.inputs import Incar, Kpoints, Potcar
from pymatgen.core import Structure
from functools import reduce
import sys, os
from pymatgen.core.sites import Site

from ccnb.vasp_parse import parse_potcar
from ccnb.path_poscar import path_poscar


class vasp_inputs(object):

    def __init__(self, struc):
        self.incar_params = {
            'pp': 'pbe',
            'istart': 0,
            'icharg': 2,
            'lwave': False,
            'lcharg': False,
            'lreal': 'auto',
            'prec': 'accurate',
            'gga': 'pe',
            'encut': 400,  #!!!!
            'ediff': 1e-4,  #!!!!
            'ediffg': -0.05,  #!!!!
            'nelm': 60,  #!!!!
            'nelmin': 2,  #!!!!
            'nsw': 1000,  #!!!!
            'isif': 2,  #!!!!
            'sigma': 0.2,  #!!!!
            'npar': 2,
            'ibrion': 2,
            'ismear': -5
        }
        self.neb_incar = {
            'spring': -5,
            'lclimb': True,
            'iopt': 1,  #!!!!
            'potim': 0,
            'ichain': 0,
            'ibrion': 3,  #!!!!
            'ismear': 0,
            'algo': 'fast'
        }
        self.vasp_path = 'vasp'
        self.struc = struc
        self.kpts = ()
        self.symbols = []
        self.set_kpts()
        self.set_symbols()
        self.dir = ''

    def set_incar(self, dict):
        for k in dict:
            self.incar_params[k] = dict[k]

    def set_lsf(self, vasp_path):
        self.vasp_path = vasp_path

    def set_kpts(self):
        self.kpts = ((int(round(30 / self.struc.lattice.a)),
                      int(round(30 / self.struc.lattice.b)),
                      int(round(30 / self.struc.lattice.c))), )

    def set_symbols(self):
        func = lambda x, y: x if y in x else x + [y]
        species = reduce(func, [
            [],
        ] + self.struc.species)
        for i in range(len(species)):
            self.symbols.append(str(species[i]))

    def potcar(self):
        try:
            potc = Potcar(symbols=self.symbols)  #functional=u'PBE' default
            potc.write_file(str(self.dir) + '/POTCAR')
        except OSError as errorinfo:
            print(str(errorinfo))
            return 0

    def incar(self, neb):
        file = open(str(self.dir) + '/POTCAR')
        p = parse_potcar(file.read())
        self.set_incar({'encut': 1.5 * p.get_max_enmax()})
        file.close()

        if neb == True:
            self.incar_params.update(self.neb_incar)
        incar = Incar(self.incar_params)
        incar.write_file(str(self.dir) + '/INCAR')

    def kpoints(self):
        kp = Kpoints(comment=u'Automatic mesh',
                     kpts=self.kpts,
                     style=Kpoints.supported_modes.Monkhorst)
        kp.write_file(str(self.dir) + '/KPOINTS')

    def lsf(self, neb):
        n = 16
        if neb:
            n = 16 * self.incar_params['images']
        filename = os.path.split(self.dir)[1]
        file = open(str(self.dir) + '/vasp.lsf', 'wb')
        file.write('#!/bin/sh\n'.encode('ascii'))
        file.write('#BSUB -cwd.\n'.encode('ascii'))
        file.write('#BSUB -e %J.err\n'.encode('ascii'))
        file.write('#BSUB -o %J.out\n'.encode('ascii'))
        file.write('#BSUB -q priority\n'.encode('ascii'))
        file.write(('#BSUB -n ' + str(n) + '\n').encode('ascii'))
        file.write(('#BSUB -J ' + str(filename) + '\n').encode('ascii'))
        file.write('#BSUB -x\n'.encode('ascii'))
        file.write(
            'ncpus=`cat $LSB_DJOB_HOSTFILE | wc -l ` \n'.encode('ascii'))
        file.write(('mpirun  -hostfile $LSB_DJOB_HOSTFILE  -np  ${ncpus} ' +
                    self.vasp_path + '\n').encode('ascii'))
        file.close()

    def inputs(self, dir, neb=False):
        self.dir = dir
        self.potcar()
        self.incar(neb)
        self.kpoints()
        self.lsf(neb)
        path_poscar(dir + '/POSCAR1', dir + '/POSCAR2', dir + '/path.cif')


#path includes the start and end sites
def interpolate(dir, struc1, struc2, path):
    nimages = len(path) - 1
    images = struc1.interpolate(struc2, nimages, True)
    dir = os.dirname()

    i = 0
    for struc in images:
        struc.translate_sites(0,
                              path[i] - struc.sites[0].frac_coords,
                              frac_coords=True,
                              to_unit_cell=True)
        num = ('%02d' % i)
        if not os.path.exists(dir + '/' + num):
            os.mkdir(dir + '/' + num)
        struc.to(filename=dir + '/' + num + '/POSCAR')
        i = i + 1


# struc = Structure.from_file('test/POSCAR')
# vi = vasp_inputs('Li',struc)
# vi.set_incar({'nelect':128})
# vi.inputs('test',False)
# vi.path_poscar([1,1,1])
