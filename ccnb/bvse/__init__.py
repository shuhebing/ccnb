from pickle import NONE
import ase.io
import os
import ccnb.bvse.Structure as Structure
import ccnb.bvse.BVAnalysis as BVAnalysis


def bv_calculation(filename,
                   moveion='Li',
                   valenceofmoveion=1,
                   resolution=0.1,
                   progress=None):
    atoms = ase.io.read(filename, store_tags=True)
    struc = Structure.Structure()
    struc.GetAseStructure(atoms)
    bvs = BVAnalysis.BVAnalysis()
    bvs.SetStructure(struc)
    bvs.SetMoveIon(moveion)
    bvs.ValenceOfMoveIon = valenceofmoveion
    bvs.SetLengthResolution(resolution)
    bvs.CaluBVSE(progress)
    bv_data = bvs.get_data()
    bvs.SaveBVSEData(os.path.splitext(filename)[0])
    bvs.SaveData(
        os.path.splitext(filename)[0] + '.bvse.grd', bv_data['BVSE'], 'BVSE')
    return bvs.Ea['BVSE']
