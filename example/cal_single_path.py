from pymatgen.core.structure import Structure
from pymatgen.core.sites import PeriodicSite
from scipy.interpolate import RegularGridInterpolator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import numpy as np
import math
import copy
import numpy.linalg as la
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline as spline


def pathfind(start, end, V, n_images=21, dr=None, h=0.1, k=0.17, min_iter=100, max_iter=10000, max_tol=5e-6):
    if not dr:
        dr = np.array([1.0 / (V.shape[0] - 1), 1.0 / (V.shape[1] - 1), 1.0 / (V.shape[2] - 1)])
    else:
        dr = np.array(dr, dtype=float)
    keff = k * dr * n_images
    h0 = h
    g1 = np.linspace(0, 1, n_images)
    s0 = start
    s1 = end
    s = np.array([g * (s1 - s0) for g in g1]) + s0
    ds = s - np.roll(s, 1, axis=0)
    ds[0] = (ds[0] - ds[0])
    ls = np.cumsum(la.norm(ds, axis=1))
    ls = ls / ls[-1]
    fi = interp1d(ls, s, axis=0)
    s = fi(g1)
    # Evaluate initial distances (for elastic equilibrium)
    ds0_plus = s - np.roll(s, 1, axis=0)
    ds0_minus = s - np.roll(s, -1, axis=0)
    ds0_plus[0] = (ds0_plus[0] - ds0_plus[0])
    ds0_minus[-1] = (ds0_minus[-1] - ds0_minus[-1])
    dV = np.gradient(V)
    for step in range(0, max_iter):
        if step > min_iter:
            h = h0 * np.exp(-2.0 * (step - min_iter) / max_iter)
        else:
            h = h0
        # Calculate forces acting on string
        d = V.shape
        s0 = s
        edV = np.array([[dV[0][int(pt[0]) % d[0]][int(pt[1]) % d[1]][int(pt[2]) % d[2]] / dr[0],

                         dV[1][int(pt[0]) % d[0]][int(pt[1]) % d[1]][int(pt[2]) % d[2]] / dr[1],

                         dV[2][int(pt[0]) % d[0]][int(pt[1]) % d[1]][int(pt[2]) % d[2]] / dr[2]] for pt in s])

        # Update according to force due to potential and string elasticity
        ds_plus = s - np.roll(s, 1, axis=0)
        ds_minus = s - np.roll(s, -1, axis=0)
        ds_plus[0] = (ds_plus[0] - ds_plus[0])
        ds_minus[-1] = (ds_minus[-1] - ds_minus[-1])
        Fpot = edV
        Fel = keff * (la.norm(ds_plus) - la.norm(ds0_plus)) * (ds_plus / la.norm(ds_plus))
        Fel += keff * (la.norm(ds_minus) - la.norm(ds0_minus)) * (ds_minus / la.norm(ds_minus))
        s = s - h * (Fpot + Fel)
        s[0] = s0[0]
        s[-1] = s0[-1]
        ds = s - np.roll(s, 1, axis=0)
        ds[0] = (ds[0] - ds[0])
        ls = np.cumsum(la.norm(ds, axis=1))
        ls = ls / ls[-1]
        fi = interp1d(ls, s, axis=0)
        s = fi(g1)
        tol = la.norm((s - s0) * dr) / n_images / h
        if (tol > 1e9):
            s = [[0, 0, 0]]
            break
        if (step > min_iter and tol < max_tol):
            break
    return s


def get_dis(struc,p1, p2):
    temp_site1 = PeriodicSite('Ar', p1, struc.lattice)
    temp_site2 = PeriodicSite('Ar', p2, struc.lattice)
    dis = temp_site1.distance(temp_site2)
    return dis


def calpointenergy(energy, point):
    x = np.linspace(0, energy.shape[0] - 1, energy.shape[0])
    y = np.linspace(0, energy.shape[1] - 1, energy.shape[1])
    z = np.linspace(0, energy.shape[2] - 1, energy.shape[2])
    interpolating_function = RegularGridInterpolator((x, y, z), energy)
    return interpolating_function(np.array(point))


def pathenergy(energy,struc, p):
    energy_path = []
    for point_of_p in p:
        point_temp = [0.0,0.0,0.0]
        for i in range(len(point_temp)):
            if point_of_p[i] < 0:
                 point_temp[i] = (point_of_p[i]+1)*(energy.shape[i]-1)
            elif point_of_p[i] <= 1:
                point_temp[i] = point_of_p[i] * (energy.shape[i]-1)
            else:
                point_temp[i] = (point_of_p[i]-1) * (energy.shape[i] - 1)
        energy_point_temp = calpointenergy(energy, point_temp)
        energy_path.append(energy_point_temp)
    dis_path = [0.0]
    tol = 0
    for i in range(len(p)-1):
        dist = get_dis(struc,p[i], p[i+1])
        tol += dist
        dis_path.append(tol)
    for i in range(len(dis_path)):
        dis_path[i] = (dis_path[i]/dis_path[-1]) *100
    migs = np.zeros(shape=(len(dis_path), 2))
    for i in range(len(energy_path)):
        migs[i][1] = energy_path[i]
        migs[i][0] = dis_path[i]
    return migs


def calpath(endpoints,struc, energy,count_images = 11):
    paths = []
    for i in range(len(endpoints)-1):
        start_f = endpoints[i]
        end_f = endpoints[i+1]
        start = np.array([start_f[0] * (energy.shape[0] - 1), start_f[1] * (energy.shape[1] - 1),
                          start_f[2] * (energy.shape[2] - 1)])
        end = np.array(
            [end_f[0] * (energy.shape[0] - 1), end_f[1] * (energy.shape[1] - 1), end_f[2] * (energy.shape[2] - 1)])
        #count_images = max(5, int(get_dis(struc,start_f, end_f) // 0.2))
        p= pathfind(start, end, energy, n_images=count_images,
                     dr=[struc.lattice.a / (energy.shape[0] ),
                         struc.lattice.b / (energy.shape[1] ),
                         struc.lattice.c / (energy.shape[2] )], h=0.01, k=0.17, min_iter=100,
                     max_iter=50000, max_tol=5e-6)
        print(p)
        for p1 in p:
            p1[0] = round(p1[0] / (energy.shape[0] - 1), 5)
            p1[1] = round(p1[1] / (energy.shape[1] - 1), 5)
            p1[2] = round(p1[2] / (energy.shape[2] - 1), 5)
        if i == 0:
            paths = copy.deepcopy(p)
        else:
            paths = np.append(paths, p[1:], axis=0)
    return paths


def paths_2_cif(paths, p_e, filename):
    with open(filename.split(".")[0] + 'path_energy_message.txt', 'w') as f:
        for i in range(len(paths)):
            f.write(str(p_e[i][0]) + '  ' + str(p_e[i][1]) + '\n ')
        for i in range(len(paths)):
            f.write(str(paths[i][0]) + '  ' + str(paths[i][1]) + ' ' + str(paths[i][1]) + '\n ')
    struc1 = Structure.from_file(filename)
    #struc.to(fmt='POSCAR', filename=filename.split(".")[0] + '_POSCAR')
    #struc1 = Structure.from_file(filename.split(".")[0] + '_POSCAR')
    for j in range(1, len(paths) - 1):
        struc1.insert(0, 'He', paths[j])
    struc1.to(fmt="CIF", filename=filename.split(".")[0] + 'path')


def showpathsenergy(p_e):
    xcoords = []
    ycoords = []
    for j in range(len(p_e)):
        xcoords.append(p_e[j][0])
        ycoords.append(p_e[j][1])

    x_new = np.linspace(min(xcoords), max(xcoords), 100)
    y_smooth = spline(xcoords, ycoords)(x_new)
    plt.figure(figsize=(6, 4))
    plt.scatter(xcoords, ycoords, color='k', marker='o')
    plt.plot(x_new, y_smooth,linewidth=2, color ='k')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlabel("Migration path(%) ", fontsize=20)
    plt.ylabel("Energy (eV)", fontsize=20)
    plt.show()


if __name__ == "__main__":
    filename_BVSE = 'D:\\temp\\zhouzy\\JF1-new-2_BVSE.npy'
    filename_CIF = 'D:\\temp\\zhouzy\\JF1-new-2-hb.cif'
    struc = Structure.from_file(filename_CIF)
    energy = np.load(filename_BVSE)
    energy = energy - np.amin(energy)
    #Li5->Li4->Li1->L2->Li3
    endpoints = [[0.754 , 0.353,  0.361],
                 [0.746 , 0.147 , 0.639]]
    p =calpath(endpoints, struc,energy,count_images=10)
    #count_images 是每段路径插入点的个数
    print(p)  #输出坐标
    p_e = pathenergy(energy, struc, p)
    print(p_e)  #输出能量
    showpathsenergy(p_e)


