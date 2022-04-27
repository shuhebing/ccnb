#author:cst
#
#input files:
#1.POSCAR(parent relaxed structure)
#2.POSCAR1 POSCAR2 (relaxed start and end structure)
#
#parameter:
#1.nelect for neb calculation
#2.migration_ion, default value 'Li'
#
#out file:
#neb packages for vasp:include POSCAR1 POSCAR2 (start and end structure nedd to relax)

import numpy as np
import sys, os
from pymatgen.core import Structure

from ccnb.zip_paths import zip_path
from ccnb.vasp_inputs import vasp_inputs


class neb_packages(object):
    #dir: the poscar dir
    #filename: structure name,such as icsd_5
    def init(self, filename, ion='Li'):
        self.nelect = 128  ######notice this
        self.migration_ion = ion
        #self.paths=[] #all paths from start to end
        self.filename = filename
        self.dir = os.path.dirname(filename)
        struc = Structure.from_file(self.dir + '/POSCAR')
        self.vi = vasp_inputs(struc)
        self.vi.set_incar({'nelect': self.nelect})

    def deal_path_periodicity(self, path):
        new_path = []
        num = 0
        #print('path:',path)
        for coord in path:
            coord = np.mod(coord, 1)
            #if(num%2==0):
            new_path.append(coord)
            #num=num+1
        #print('new_path:',new_path)
        return new_path

    def from_file(self, path_file):
        self.paths = np.load(path_file)

    def from_list(self, paths):
        self.paths = paths

    def data_parse(self):
        if not os.path.exists(self.dir + '/paths'):
            os.mkdir(self.dir + '/paths')
        num = 0
        for path in self.paths:
            # 处理坐标周期性问题
            p = self.deal_path_periodicity(path)
            dir = self.dir + '/paths/path_' + str(num)
            # dir=self.paths_dir+'/path_'+str(num)
            if self.judge_site(p[0], p[-1], dir):
                continue
            struc = Structure.from_file(self.dir + '/POSCAR')
            struc.to(filename=dir + '/POSCAR')
            for i in range(1, len(p) - 1):
                struc.insert(0, 'He', p[i])
            struc.to(fmt='cif', filename=dir + '/path.cif')
            struc1 = Structure.from_file(dir + '/POSCAR_base')
            struc1.insert(
                0, self.migration_ion,
                path[0])  # at the first site insert the interval point
            struc1.to(filename=dir + "/POSCAR1")
            struc2 = Structure.from_file(dir + '/POSCAR_base')
            struc2.insert(
                0, self.migration_ion,
                path[-1])  # at the first site insert the interval point
            struc2.to(filename=dir + "/POSCAR2")
            self.vi.set_incar({'images': len(p) - 2})

            self.vi.inputs(dir, True)
            num = num + 1
        zip_path(self.dir + '/paths', self.filename + '_neb_paths.zip')

    # product the base poscar for insert site
    def judge_site(self, p1, p2, dir):
        start = None
        end = None
        # print('start:%s', p1)
        # print('end:%s', p2)
        struc = Structure.from_file(self.dir + '/POSCAR')
        struc.insert(0, self.migration_ion, p2)
        struc.insert(0, self.migration_ion, p1)
        dis = struc.get_distance(0, 1)
        # print(dir, 'dis:', dis)
        # if (dis > 6): #if distance greater than 5, then can not participate in neb calculation
        # print('dis too long!!!!!!!!')
        # return 1
        for i in range(2, len(struc.sites)):
            # print('%s:%s', i,struc.sites[i].frac_coords)
            dis1 = struc.get_distance(i, 0)
            dis2 = struc.get_distance(i, 1)
            # print('dis:',dis1,dis2)
            if dis1 < 0.5:
                start = i
                # print('dis1:',dis1)
            if dis2 < 0.5:
                end = i
                # print('dis2:',dis2)
        # print('start:%s', start)
        # print('end:%s',end)
        if start:
            struc.sites.pop(start)
            if end:
                if end > start:
                    struc.sites.pop(end - 1)
                elif end < start:
                    struc.sites.pop(end)
        elif end:
            struc.sites.pop(end)
        else:
            print('no site')
            return 1
        struc.sites.pop(0)
        struc.sites.pop(0)
        if not os.path.exists(dir):
            os.mkdir(dir)
        struc.to(filename=dir + '/POSCAR_base')
