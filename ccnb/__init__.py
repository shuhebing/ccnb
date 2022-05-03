import os
from pathlib import Path
from struct import unpack
from monty.io import zopen
import numpy as np
import pymatgen
from pymatgen.core import Structure
from pymatgen.core import SETTINGS
from cavd.local_environment import CifParser_new
from ccnb.bvse_cavd import MigrationNetwork
from ccnb.mergecluster import Void, Channel
from ccnb.neb_packages import neb_packages
from ccnb.bvse import bv_calculation
from ccnb.cavd_channel import cal_channel_cavd
from rich.progress import Progress


class Progressbar:

    def __init__(self):
        self.progress = Progress()
        self.task = None

    def set_task(self, workstr: str, total):
        self.task = self.progress.add_task(workstr, total)

    def updatetask(self, i: int):
        self.progress.update(self.task, advance=0.9)
        print('bvse compute have finished {0:.2f}%'.format(i))


def get_channel_cavd(filename,
                     migrant,
                     ntol=0.02,
                     rad_flag=True,
                     lower=0,
                     upper=10.0,
                     rad_dict=None):
    """
    calculate interstitial network and Voronoi network by CAVD
    :param filename: CIF filename
    :param migrant: mobile ion
    :param lower: the lower threshold  of ion migration
                  to select the suitable interstitial network
    :param upper: the upper threshold  of ion migration
                  to select the suitable interstitial network
    :return: RT value, NET file for saving
             interstitial network and Voronoi network
    """
    conn_val = cal_channel_cavd(filename,
                                migrant,
                                ntol=ntol,
                                rad_flag=rad_flag,
                                lower=lower,
                                upper=upper,
                                rad_dict=rad_dict)
    return conn_val


def get_bvse(filename_cif,
             moveion='Li',
             valenceofmoveion=1,
             resolution=0.1,
             progress=None):
    """
    calculate BVSE landscape
    :param filename_cif: cif filename
    :param moveion:  move ion
    :param valenceofmoveion: valence of move ion
    :param resolution: resolution for calculating BVSE landsacpe
    :return: migration energy barrier, numpy Binary.
             Pgrid file for saving and visualization BVSE landscape
    """
    prgbar = Progressbar()
    prgbar.set_task('bvse computing progress...', total=100)
    barrier = bv_calculation(filename_cif,
                             moveion=moveion,
                             valenceofmoveion=valenceofmoveion,
                             resolution=resolution,
                             progress=prgbar.updatetask)
    return barrier


def load_struc(filename_cif):
    """
    read structure from CIF file
    :param filename_cif: cif filename
    :return: Structure object defined by pymatgen
    """
    with zopen(filename_cif, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    structure = parser.get_structures(primitive=False)[0]
    return structure


def load_bvse_from_npy(filename_bvse):
    """
    Read BVSE potential from numpy binary file
    :param filename_bvse: numpy binary filename
    :return: Three dimensional array
    """
    energy_landscape = np.load(filename_bvse)
    return energy_landscape


def load_bvse_from_pgrid(filename_bvse):
    '''
     Read BVSE potential from pgrid file
    :param filename_bvse: potential filename
    :return: three dimensional array
    '''
    binfile = open(filename_bvse, 'rb')
    data = binfile.read(112)
    data1 = binfile.read(12)
    num = unpack('3i', data1)
    total_num = str(num[0] * num[1] * num[2]) + 'f'
    data2 = binfile.read(28)
    data_num = binfile.read()
    num_energy = unpack(total_num, data_num)
    energy_landscape = np.array(num_energy).reshape((num[0], num[1], num[2]))
    return energy_landscape


def load_voids_channels_from_file(filename_cavd):
    """
    Read interstices and channel segments from net file calculated by cavd
    :param filename_cavd: NET filename
    :return: dict, interstices and channel segments
    """
    voids_dict = {}
    channels_dict = {}
    flag_p = 0
    flag_n = 0
    file = open(filename_cavd, 'r')
    for line in file.readlines():
        if 'Interstitial' in line:
            flag_p = 1
            flag_n = 0
            continue
        if 'Connection' in line:
            flag_p = 0
            flag_n = 1
            continue
        if flag_p == 1:
            line = line.split()
            if len(line) > 3:
                void = Void()
                void.id = int(line[0])
                void.label = int(line[1])
                void.coord = [
                    np.float64(line[2]),
                    np.float64(line[3]),
                    np.float64(line[4])
                ]
                void.radii = np.float64(line[5])
                voids_dict[void.id] = void
        if flag_n == 1:
            line = line.split()
            if len(line) > 4:
                channel = Channel()
                channel.start = int(line[0])
                channel.end = int(line[1])
                channel.phase = [int(line[2]), int(line[3]), int(line[4])]
                channel.coord = [
                    np.float64(line[5]),
                    np.float64(line[6]),
                    np.float64(line[7])
                ]
                channel.radii = np.float64(line[8])
                channels_dict[(channel.start, channel.end)] = channel
    return voids_dict, channels_dict


def get_non_equivalent_paths_between_latticesite(filename_CIF,
                                                 filename_BVSE,
                                                 filename_CAVD,
                                                 energythreshold,
                                                 moveion='Li'):
    """
    get nonequivalent transport pathways between lattice sites
    :param filename_CIF:  CIF filename
    :param filename_BVSE: BVSE landscape numpy binary
    :param filename_CAVD: net filename
    :param energythreshold: energythreshold is used to delete path segments 
           in the transport network
    that are higher than this threshold
    :param moveion: move ion
    :return: the fractional coordinates of each image on the MEP
    """
    voids, channels = load_voids_channels_from_file(filename_CAVD)
    struc = load_struc(filename_CIF)
    energy = load_bvse_from_npy(filename_BVSE)
    mn = MigrationNetwork(struc,
                          energy,
                          voids,
                          channels,
                          filename_CIF,
                          moveion=moveion,
                          ismergecluster=True,
                          energythreshold=energythreshold,
                          iscalnonequalchannels=False)
    mn.cal_nonequl_paths()
    mn.save_data(filename_CIF)
    mn.showenergy(filename_CIF)
    return mn
    # return mn.paths_position，return mn.paths_position


def get_migration_networks_voids(filename_CIF,
                                 filename_BVSE,
                                 filename_CAVD,
                                 energythreshold,
                                 moveion='Li',
                                 mergecluster=True,
                                 clusterradii=0.75):
    """
    获得迁移网络中的间隙点
    """
    voids, channels = load_voids_channels_from_file(filename_CAVD)
    struc = Structure.from_file(filename_CIF)
    energy = load_bvse_from_npy(filename_BVSE)
    mn = MigrationNetwork(struc,
                          energy,
                          voids,
                          channels,
                          filename_CIF,
                          moveion=moveion,
                          ismergecluster=mergecluster,
                          energythreshold=energythreshold,
                          iscalnonequalchannels=False,
                          clusterradii=clusterradii)
    for void in mn._voids.values():
        struc.insert(0, 'He', void.coord)
    struc.to(fmt='cif', filename=filename_CIF + '_voids.cif')
    return mn


def get_non_equivalent_paths_between_voids(filename_CIF,
                                           filename_BVSE,
                                           filename_CAVD,
                                           energythreshold=None,
                                           moveion='Li'):
    voids, channels = load_voids_channels_from_file(filename_CAVD)
    struc = load_struc(filename_CIF)
    energy = load_bvse_from_npy(filename_BVSE)
    mn = MigrationNetwork(struc,
                          energy,
                          voids,
                          channels,
                          filename_CIF,
                          moveion=moveion,
                          ismergecluster=True,
                          energythreshold=energythreshold,
                          iscalnonequalchannels=True)
    return mn.cal_nonequal_mep_between_voids()


def configure_neb_packet(filename_CIF, mep, moveion="Li"):
    struc1 = load_struc(filename_CIF)
    path_cif_file = Path(filename_CIF).parent.joinpath('POSCAR')

    struc1.to(fmt='POSCAR', filename=str(path_cif_file))
    file_prefix = filename_CIF.split('.')[0]
    #SETTINGS['PMG_VASP_PSP_DIR'] = os.path.abspath("..") + '\ccnb\psp_resources'
    n = neb_packages()
    n.init(file_prefix, moveion)
    n.from_list(mep)
    n.data_parse()


def all_cal(filename_cif,
            moveion='Li',
            valenceofmoveion=1,
            energythreshold=None):
    bv_calculation(filename_cif,
                   moveion=moveion,
                   valenceofmoveion=valenceofmoveion,
                   resolution=0.1)
    cal_channel_cavd(filename_cif,
                     moveion,
                     ntol=0.02,
                     rad_flag=True,
                     lower=0,
                     upper=10.0,
                     rad_dict=None)
    filename_BVSE = filename_cif.split('.')[0] + '.npy'
    filename_CAVD = filename_cif.split('.')[0] + '_origin.net'
    meps = get_non_equivalent_paths_between_latticesite(filename_cif,
                                                        filename_BVSE,
                                                        filename_CAVD,
                                                        energythreshold,
                                                        moveion=moveion)
    configure_neb_packet(filename_cif, meps, moveion=moveion)
