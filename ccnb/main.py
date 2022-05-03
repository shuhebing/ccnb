"""


"""
import argparse
from pathlib import Path
from ccnb import cal_channel_cavd
from ccnb import get_bvse
from ccnb import get_migration_networks_voids


def get_voids(filename_CIF,
              energythreshold,
              moveion,
              mergecluster=True,
              clusterradii=0.75):
    filename_CIF = Path(filename_CIF)
    filename_CAVD = str(filename_CIF.parent.joinpath(
        filename_CIF.stem)) + ".net"
    filename_BVSE = str(filename_CIF.parent.joinpath(
        filename_CIF.stem)) + "_BVSE.npy"
    if not Path(filename_CAVD).exists():
        print('{} file no exist,\n please cavd calculate first'.format(
            filename_CAVD))
        return None
    if not Path(filename_BVSE).exists():
        print('{} file no exist,\n please bvse calculate first'.format(
            filename_BVSE))
        return None
    mn = get_migration_networks_voids(str(filename_CIF),
                                      filename_BVSE,
                                      filename_CAVD,
                                      energythreshold,
                                      moveion=moveion,
                                      mergecluster=mergecluster,
                                      clusterradii=clusterradii)
    return mn


def main():
    parser = argparse.ArgumentParser(description='ccnb 电解质材料离子输运特征计算分析程序')
    parser.add_argument('-f',
                        '--structure_file',
                        dest='struct_file',
                        type=str,
                        required=True,
                        help='结构文件名(cif格式)')
    parser.add_argument('-i',
                        '--move_ion',
                        dest='ion',
                        type=str,
                        default='Li',
                        required=True,
                        help='迁移离子类型')
    parser.add_argument('-c',
                        '--calculation_type',
                        dest='cal_type',
                        type=str,
                        default='bvse',
                        choices=['bvse', 'cavd', 'find_voids'],
                        required=True,
                        help='计算类型选择')
    bvse = parser.add_argument_group('bvse')
    bvse.add_argument('-v',
                      '--valence',
                      type=float,
                      default=1,
                      dest='valence',
                      help='迁移离子化合价')
    bvse.add_argument('-r',
                      '--resolution',
                      type=float,
                      default=0.1,
                      dest='resolution',
                      help='计算区域分辨率（单位埃）')
    cavd = parser.add_argument_group('cavd')
    cavd.add_argument('--ntol',
                      dest='ntol',
                      type=float,
                      default=0.02,
                      help='计算容限')
    cavd.add_argument('--radius_flag',
                      dest='rad_flag',
                      type=bool,
                      default=True,
                      help='几何分析时是否考虑离子半径')
    cavd.add_argument('-l',
                      '--lower',
                      dest='lower',
                      type=float,
                      default=0.5,
                      help='通道大小下限值(单位埃)')
    cavd.add_argument('-u',
                      '--upper',
                      dest='upper',
                      type=float,
                      default=1.0,
                      help='通道大小上限值(单位埃)')

    find_voids = parser.add_argument_group('find_voids')
    find_voids.add_argument('--energy',
                            dest='energy',
                            type=float,
                            default=0.5,
                            help='能量筛选上限')
    find_voids.add_argument('--cluster',
                            dest='cluster',
                            type=bool,
                            default=True,
                            help='是否进行合并团簇间隙点')
    find_voids.add_argument('-cr',
                            '--cluster_radii',
                            dest='cluster_radii',
                            type=float,
                            default=0.8,
                            help='合并团簇范围半径')
    args = parser.parse_args()
    if args.cal_type == 'cavd':
        RT = cal_channel_cavd(args.struct_file,
                              migrant=args.ion,
                              ntol=args.ntol,
                              rad_flag=args.rad_flag,
                              lower=args.lower,
                              upper=args.upper,
                              rad_dict=None)
        print(RT)
    if args.cal_type == 'bvse':
        barrier = get_bvse(args.struct_file,
                           moveion=args.ion,
                           valenceofmoveion=args.valence,
                           resolution=args.resolution)
        print(barrier)
    if args.cal_type == 'find_voids':
        mn = get_voids(filename_CIF=args.struct_file,
                       energythreshold=args.energy,
                       moveion=args.ion,
                       mergecluster=args.cluster,
                       clusterradii=args.cluster_radii)
        if mn:
            print(mn)
        else:
            print('Calculation failed!!!')


if __name__ == "__main__":
    main()
