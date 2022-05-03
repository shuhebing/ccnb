
from cavd.local_environment import CifParser_new
from ccnb.bvse_cavd import MigrationNetwork
from ccnb.mergecluster import Void, Channel
from ccnb.neb_packages import neb_packages
from ccnb.bvse import bv_calculation
from ccnb import (load_struc,load_bvse_from_npy,
          load_voids_channels_from_file,get_bvse)
from ccnb.cavd_channel import cal_channel_cavd
from ccnb.MigrationNetworkLatticeSites import MigrationNetworkLatticeSites


if __name__ == "__main__":
    filename_CIF = 'tzw_2021032612573832.cif'
    # calculate interstital  network by CAVD
    #RT =cal_channel_cavd(filename_CIF, 'Li', ntol=0.02, rad_flag=True, lower=0.4, upper=10.0, rad_dict=None)
    #print(RT)
    # calculate bond valence site example
    barrier = get_bvse(filename_CIF, moveion='Li', valenceofmoveion=1, resolution=0.1)
    print(barrier)
    #filename_CAVD = "icsd_221\\icsd_221.net"
    #filename_BVSE = "icsd_221\\icsd_221_BVSE.npy"
    #voids, channels = load_voids_channels_from_file(filename_CAVD)
    #struc = load_struc(filename_CIF)
    #energy = load_bvse_from_npy(filename_BVSE)
    #print(voids)
    #mn = MigrationNetworkLatticeSites(struc, energy, voids, channels, filename_CIF, moveion="Li")
    #mn.cal_paths()
    #mn.save_data(filename_CIF)



