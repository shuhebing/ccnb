"""


"""
import argparse
from pathlib import Path
from ccnb import cal_channel_cavd
from ccnb import get_bvse
def main():
    parser = argparse.ArgumentParser(description='ccnb 电解质材料离子输运特征计算分析程序')
    parser.add_argument('-f', '--structure_file',dest='struct_file', type=str,
                        required=True, help='结构文件名(cif格式)')
    parser.add_argument('-mi', '--move_ion', dest='ion',
                        type=str, default='Li', required=True, help='迁移离子类型')
    parser.add_argument('-c', '--calculation_type', dest='cal_type', type=str,
                        default='bvse', choices=['bvse', 'cavd', 'all'], required=True, help='计算类型选择')
    bvse=parser.add_argument_group('bvse')
    bvse.add_argument('-v','--valence',type=float,default=1,dest='valence',help='迁移离子化合价')
    bvse.add_argument('-r','--resolution',type=float,default=0.1,dest='resolution',help='计算区域分辨率（单位埃）')
    cavd=parser.add_argument_group('cavd')
    cavd.add_argument('--ntol',dest='ntol',type=float,default=0.02,help='计算容限')
    cavd.add_argument('--radius_flag',dest='rad_flag',type=bool,default=True,help='几何分析时是否考虑离子半径')
    cavd.add_argument('-l','--lower',dest='lower',type=float,default=0.5,help='通道大小下限值(单位埃)')
    cavd.add_argument('-u','--upper',dest='upper',type=float,default=1.0,help='通道大小上限值(单位埃)')
    args = parser.parse_args()
    if args.cal_type == 'cavd':
        RT =cal_channel_cavd(args.struct_file, migrant=args.ion, ntol=args.ntol, rad_flag=args.rad_flag, lower=args.lower, upper=args.upper, rad_dict=None)
        print(RT)
    if args.cal_type == 'bvse':
        barrier = get_bvse(args.struct_file, moveion=args.ion, valenceofmoveion=args.valence, resolution=args.resolution)
        print(barrier)
if __name__ == "__main__":
    main()
    