import os, re
from monty.io import zopen
import numpy as np
from cavd.netio import *
from cavd.channel import Channel
from cavd.netstorage import AtomNetwork, connection_values_list
from cavd.local_environment import CifParser_new, LocalEnvirCom
from cavd.get_Symmetry import get_labeled_vornet
from pymatgen.core. structure import Structure

def cal_channel_cavd(filename,
                     migrant,
                     ntol=0.02,
                     rad_flag=True,
                     lower=0,
                     upper=10.0,
                     rad_dict=None):
    with zopen(filename, "rt") as f:
        input_string = f.read()
    parser = CifParser_new.from_string(input_string)
    stru = parser.get_structures(primitive=False)[0]
    #stru = Structure.from_file(filename)
    #species = [str(sp).replace("Specie ", "") for sp in stru.species]
    #elements = [re.sub('[^a-zA-Z]', '', sp) for sp in species]
    elements = stru.symbol_set
    #symm_number, symm_sybol = parser.get_symme()
    sitesym = parser.get_sym_opt()
    if migrant not in elements:
        raise ValueError(
            "The input migrant ion not in the input structure! Please check it."
        )
    effec_radii, migrant_radius, migrant_alpha, nei_dises, coordination_list = LocalEnvirCom(
        stru, migrant)
    radii = {}
    if rad_flag:
        if rad_dict:
            radii = rad_dict
        else:
            radii = effec_radii
    print(radii)
    atmnet = AtomNetwork.read_from_RemoveMigrantCif(filename, migrant, radii,
                                                    rad_flag)
    vornet, edge_centers, fcs, faces = atmnet.perform_voronoi_decomposition(
        True, ntol)
    # add_fcs_vornet = vornet.add_facecenters(faces)
    prefixname = filename.replace(".cif", "")
    symprec = 0.01
    sym_vornet, voids = get_labeled_vornet(vornet, sitesym, symprec)
    writeNETFile(prefixname + "_origin.net", atmnet, sym_vornet)
    channels = Channel.findChannels2(vornet, atmnet, lower, upper,
                                     prefixname + ".net")
    Channel.writeToVESTA(channels, atmnet, prefixname)
    conn_val = connection_values_list(prefixname + ".resex", sym_vornet)
    return conn_val
