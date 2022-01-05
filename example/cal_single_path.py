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
    # Set parameters
    if not dr:
        dr = np.array([1.0 / (V.shape[0] - 1), 1.0 / (V.shape[1] - 1), 1.0 / (V.shape[2] - 1)])
    else:
        dr = np.array(dr, dtype=float)
    keff = k * dr * n_images
    h0 = h
    # 初使化 string
    g1 = np.linspace(0, 1, n_images)  # 创建等差数列
    s0 = start  # 初使结构的迁移离子的坐标
    s1 = end  # 最终结构的迁移离子的坐标
    s = np.array([g * (s1 - s0) for g in g1]) + s0  # s是一个3*n的矩阵
    ds = s - np.roll(s, 1, axis=0)  # np.roll函数的意思是将s，沿着axis的方向，滚动1个长度，相当于把最后一个元素放到第一个
    ds[0] = (ds[0] - ds[0])  # 把第一个元素置为（0，0，0）
    ls = np.cumsum(la.norm(ds, axis=1))  # norm求ds的范数  cumsum计算轴向元素累加和，返回由中间结果组成的数组
    ls = ls / ls[-1]  # 归一化
    fi = interp1d(ls, s, axis=0)  # 插值
    s = fi(g1)
    # Evaluate initial distances (for elastic equilibrium)
    ds0_plus = s - np.roll(s, 1, axis=0)  # 正向
    ds0_minus = s - np.roll(s, -1, axis=0)  # 负向
    ds0_plus[0] = (ds0_plus[0] - ds0_plus[0])
    ds0_minus[-1] = (ds0_minus[-1] - ds0_minus[-1])
    dV = np.gradient(V)  # 计算梯度

    # Evolve string
    for step in range(0, max_iter):
        if step > min_iter:
            h = h0 * np.exp(-2.0 * (step - min_iter) / max_iter)  # 逐步衰减步长以减少震荡
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
        # 每次迭代保持两端值保持不变
        s[0] = s0[0]
        s[-1] = s0[-1]
        # Reparametrize string            #更新参数
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
            # print("Converged at step {}".format(step))
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
    # energy = energy - np.amin(energy)
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
        #dis_path[i] = "%.2f%%" % (dis_path[i] * 100)
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
    plt.figure(figsize=(6, 4))  # 创建绘图对象
    plt.scatter(xcoords, ycoords, color='k', marker='o')
    plt.plot(x_new, y_smooth,linewidth=2, color ='k')
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    '''
    poly = np.polyfit(xcoords, ycoords, deg=24)  # 最小二乘法多项式拟合
    z = np.polyval(poly, xcoords)
    plt.figure(figsize=(8, 6))  # 创建绘图对象
    plt.plot(xcoords, z, linewidth=2, color="black")
    plt.scatter(xcoords, ycoords, color='k', marker='o')'''
    # plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1], fontsize=35)
    #plt.yticks([0.00,0.05,0.10,0.15,0.20,0.25,0.30, 0.35], fontsize=35)
    plt.xlabel("Migration path(%) ", fontsize=20)  # X轴标签
    plt.ylabel("Energy (eV)", fontsize=20)  # Y轴标签
    plt.show()


if __name__ == "__main__":
    
    filename_BVSE = 'LGPS\\LGPS.npy'
    filename_CIF = 'LGPS\\LGPS.cif'
    # filename_cif是文件名，filename_npy是BVSE程序计算出的一个后缀名为npy的能量文件，endpoints是计算的坐标,可以是只给起始点坐标，或者指定路径中间的点
    struc = Structure.from_file(filename_CIF)
    energy = np.load(filename_BVSE)
    energy = energy - np.amin(energy)
    #Li5->Li4->Li1->L2->Li3
    endpoints = [[0.75873 , 0.78301,  0.38703],[0.50071 , 0.50837 , 0.25337 ],
                 [ 0.25019 , 0.28335 , 0.11941],
                 [0.23286,  0.25683 , 0.44804],
                 [0.25439,  0.23225 , 0.76262]]
    p =calpath(endpoints, struc,energy,count_images=11)
    #count_images 是每段路径插入点的个数
    print(p)  #输出坐标
    p_e = pathenergy(energy, struc, p)
    print(p_e)  #输出能量
    showpathsenergy(p_e)


