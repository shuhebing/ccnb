from ccnb import load_struc
from ccnb import get_non_equivalent_paths_between_latticesite
from ccnb import cal_channel_cavd
from ccnb import get_bvse


if __name__ == "__main__":
    filename_CAVD = "LGPS/LGPS.net"
    filename_BVSE = "LGPS/LGPS.npy"
    filename_CIF = "LGPS/LGPS.cif"
    struc = load_struc(filename_CIF)
    #all_cal(filename_CIF, moveion='Li', valenceofmoveion=1, energythreshold=1.4)
    c=cal_channel_cavd(filename_CIF, 'Li', ntol=0.02, rad_flag=True, lower=0, upper=10.0, rad_dict=None)
    #ea = get_bvse(filename_CIF, moveion='Mg', valenceofmoveion=2, resolution=0.1)
    #non_equivalent_paths_between_latticesite(filename_CIF, filename_BVSE, filename_CAVD, energythreshold=1.4, moveion='Li')



