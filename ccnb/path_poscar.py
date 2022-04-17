from pymatgen.core import Structure
import os, sys

# from vasp_inputs import interpolate


#product the neb packages for one path
#posc1:file of POSCAR1POSCAR1 posc2:file of POSCAR3 posc_path:file of POSCAR_path
def path_poscar(posc1, posc2, posc_path):
    struc1 = Structure.from_file(posc1)
    struc2 = Structure.from_file(posc2)
    s = Structure.from_file(posc_path)

    path = []
    path.append(struc1.sites[0].frac_coords)
    for site in s.sites:
        if site.specie.symbol == 'He':
            path.append(site.frac_coords)
    path.append(struc2.sites[0].frac_coords)

    nimages = len(path) - 1
    images = struc1.interpolate(struc2, nimages, True)
    dir = os.path.dirname(posc_path)

    i = 0
    for struc in images:
        struc.translate_sites(0,
                              path[i] - struc.sites[0].frac_coords,
                              frac_coords=True,
                              to_unit_cell=True)
        num = ('%02d' % i)
        if not os.path.exists(dir + '/' + num):
            os.mkdir(dir + '/' + num)
        struc.to(filename=dir + '/' + num + '/POSCAR')
        i = i + 1

    # interpolate(dir, struc1, struc2, path)


#path_poscar(sys.argv[0], sys.argv[1], sys.argv[2])
