import os,re
from struct import unpack
from monty.io import zopen
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen import SETTINGS
from cavd.local_environment import  CifParser_new
from ccnb.bvse_cavd import MigrationNetwork
from ccnb.mergecluster import Void, Channel
from ccnb.neb_packages import neb_packages
from ccnb.bvse import bv_calculation
from ccnb.cavd_channel import cal_channel_cavd


def load_struc(filename_cif):
    """
    从cif文件中读取结构，调用的是pymatgen库中Structure类方法
    :param filename_cif: cif文件名
    :return: pymatgen定义的Structure对象
    """
    with zopen(filename_cif, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    structure = parser.get_structures(primitive=False)[0]
    return structure


def load_bvse_from_npy(filename_bvse):
    """
    从numpy二进制文件中读取数据
    :param filename_bvse: 保存BVSE能量的文件
    :return: 三维数组
    """
    energy_landscape = np.load(filename_bvse)
    return energy_landscape


def load_bvse_from_pgrid(filename_bvse):
    '''
    从grd文件中读取数据
    :param filename_bvse:
    :return:
    '''
    binfile = open(filename_bvse, 'rb')
    data = binfile.read(112)
    data1 = binfile.read(12)
    num = unpack('3i', data1)
    total_num = str(num[0] * num[1] * num[2]) + 'f'
    data2 = binfile.read(28)
    data_num = binfile.read()
    num_energy = unpack(total_num,data_num)
    energy_landscape = np.array(num_energy ).reshape((num[0],num[1],num[2]))
    return energy_landscape


def load_voids_channels_from_file(filename_cavd):
    """
    从CAVD计算出的NET文件中读取间隙、导通信息
    :param filename_cavd: 要读取的文件名
    :return: 返回voids 和channels 两个字典
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
                void.coord = [np.float64(line[2]), np.float64(line[3]), np.float64(line[4])]
                void.radii = np.float64(line[5])
                voids_dict[void.id] = void
        if flag_n == 1:
            line = line.split()
            if len(line) > 4:
                channel = Channel()
                channel.start = int(line[0])
                channel.end = int(line[1])
                channel.phase = [int(line[2]), int(line[3]), int(line[4])]
                channel.coord = [np.float64(line[5]), np.float64(line[6]), np.float64(line[7])]
                channel.radii = np.float64(line[8])
                channels_dict[(channel.start, channel.end)] = channel
    return voids_dict, channels_dict


def get_channel_cavd(filename, migrant, ntol=0.02, rad_flag=True, lower=0, upper=10.0, rad_dict=None):
    conn_val = cal_channel_cavd(filename, migrant, ntol=ntol, rad_flag=rad_flag,
                                lower=lower, upper=upper, rad_dict=rad_dict)
    return conn_val


def get_bvse(filename_cif, moveion='Li', valenceofmoveion=1, resolution=0.1):
    ea = bv_calculation(filename_cif, moveion=moveion, valenceofmoveion=valenceofmoveion, resolution=resolution)
    return ea


def get_non_equivalent_paths_between_latticesite(filename_CIF, filename_BVSE, filename_CAVD, energythreshold, moveion='Li'):
    voids, channels = load_voids_channels_from_file(filename_CAVD)
    struc = load_struc(filename_CIF)
    energy = load_bvse_from_npy(filename_BVSE)
    mn = MigrationNetwork(struc, energy, voids, channels,filename_CIF, moveion=moveion, ismergecluster=True,
                          energythreshold=energythreshold, iscalnonequalchannels=False)

    mn.cal_nonequl_paths()
    return mn.paths_position


def get_non_equivalent_paths_between_voids(filename_CIF, filename_BVSE, filename_CAVD, energythreshold=None, moveion='Li'):
    voids, channels = load_voids_channels_from_file(filename_CAVD)
    struc = load_struc(filename_CIF)
    energy = load_bvse_from_npy(filename_BVSE)
    mn = MigrationNetwork(struc, energy, voids, channels,filename_CIF, moveion=moveion, ismergecluster=True,
                          energythreshold=energythreshold, iscalnonequalchannels=True)
    return mn.cal_nonequalpath_between_voids()


def configure_neb_packet(filename_CIF, mep, moveion="Li"):
    struc = load_struc(filename_CIF)
    struc.to(fmt='POSCAR', filename=os.path.join(filename_CIF.split('/')[0], 'POSCAR'))
    file_prefix = filename_CIF.split('.')[0]
    SETTINGS['PMG_VASP_PSP_DIR'] = 'ccnb/psp_resources'
    n = neb_packages()
    n.init(file_prefix, moveion)
    n.from_list(mep)
    n.data_parse()


def all_cal(filename_cif,moveion='Li',valenceofmoveion=1,energythreshold=None):
    bv_calculation(filename_cif, moveion=moveion, valenceofmoveion=valenceofmoveion, resolution=0.2)
    cal_channel_cavd(filename_cif, moveion, ntol=0.02, rad_flag=True, lower=0, upper=10.0, rad_dict=None)
    filename_BVSE = filename_cif.split('.')[0]+'.npy'
    filename_CAVD = filename_cif.split('.')[0] + '_origin.net'
    file_prefix = filename_cif.split('.')[0]
    meps = get_non_equivalent_paths_between_latticesite(filename_cif, filename_BVSE, filename_CAVD,
                                                   energythreshold, moveion=moveion)
    configure_neb_packet(filename_cif, meps, moveion=moveion)