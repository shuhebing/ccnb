from pathlib import Path
from ccnb import load_struc
from ccnb import get_non_equivalent_paths_between_latticesite, get_non_equivalent_paths_between_voids
from ccnb import cal_channel_cavd
from ccnb import get_bvse
from ccnb import all_cal, configure_neb_packet


if __name__ == "__main__":
    filename_CIF = Path('D:\\temp\\zhouzy\\LiTaPO-hb.cif')
    filename_CAVD = str(filename_CIF.parent.joinpath(filename_CIF.stem))+".net"
    filename_BVSE = str(filename_CIF.parent.joinpath(filename_CIF.stem))+"_BVSE.npy"
    filename_CIF = str(filename_CIF)
    # Read structure from CIF file
    #struc = load_struc(filename_CIF)
    #print(len(struc.sites))
    # calculate interstital  network by CAVD
    #RT =cal_channel_cavd(filename_CIF, 'Li', ntol=0.02, rad_flag=True, lower=0.5, upper=1.2, rad_dict=None)
    #print(RT)
    # calculate bond valence site example
    #barrier = get_bvse(filename_CIF, moveion='Li', valenceofmoveion=1, resolution=0.1)
    #print(barrier)

    # energythreshold is used to delete path segments in the transport network that are higher than this threshold
    mep = get_non_equivalent_paths_between_latticesite(filename_CIF, filename_BVSE, filename_CAVD, energythreshold=-1.0, moveion='Li')
    print(mep)
    configure_neb_packet(filename_CIF, mep)
   # all_cal(filename_CIF, moveion='Li', valenceofmoveion=1, energythreshold=1.4)
