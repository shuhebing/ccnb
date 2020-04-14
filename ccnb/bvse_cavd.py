import operator
import copy
import math
import numpy as np
import numpy.linalg as la
import networkx as nx
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.sites import PeriodicSite
import ccnb.mergecluster as mc
from ccnb.mergecluster import Void, Channel


class MigrationPath(object):
    """
    用来计算两个间隙点之间的近似最小能量路径，以及路径上每点的能量、瓶颈点的能量、迁移能垒
    """
    def __init__(self, endpoints, phases, energy, struc, countimages=11):
        """
        :param endpoints: 迁移路径端点坐标，类型为两维列表
        :param phase: 路径的周期位移矢量，类型为两维列表
        :param energy: BVSE势函数，类型为三维数组
        :param struc: 结构，类型为pymatgen structure对象
        """
        self._energy = energy
        self._struc = struc
        self._endpoints = endpoints
        self._phase = phases
        self._countimages = countimages
        self._path = []
        self._path_energy = []
        self._max_energy = None
        self._max_position = None
        self._min_energy = None
        self._barrier = None
        self.cal_path()
        self.cal_energy()

    @property
    def path(self):
        """
        :return: 返回路径上各点的坐标
        """
        return self._path

    @property
    def path_energy(self):
        """
        :return: 返回路径上各点的能量，类型为两维列表，第一列表示各点的距离，第二列表示各点的能量
        """
        return self._path_energy

    @property
    def energy_max(self):
        """
        :return: 返回瓶颈点的能量
        """
        return self._max_energy

    @property
    def position_max(self):
        """
        :return: 返回瓶颈点的能量
        """
        return self._max_position

    @property
    def energy_min(self):
        """
        :return: 返回路径上各点的坐标
        """
        return self._min_energy

    @property
    def barrier(self):
        """
        :return: 返回迁移能垒
        """
        return self._barrier

    def cal_path(self):
        """
        计算起始点之间的最小能量路径，
        :return: 返回值为两维数组，其中每一点为路径上某一点的坐标
        """
        for i in range(len(self._endpoints)-1):
            start_f = [j for j in self._endpoints[i]]
            for j in range(len(start_f)):
                if start_f[j] > 1:
                    start_f[j] -= 1
                if start_f[j] < 0.0:
                    start_f[j] += 1
            end_f = [j for j in self._endpoints[i+1]]
            for j in range(len(end_f)):
                if end_f[j] > 1:
                    end_f[j] -= 1
                if end_f[j] < 0.0:
                    end_f[j] += 1
            if operator.eq(self._phase[i], [0, 0, 0]):
                start = np.array([start_f[0]*(self._energy.shape[0] - 1),
                                  start_f[1] * (self._energy.shape[1] - 1),
                                  start_f[2] * (self._energy.shape[2] - 1)])
                end = np.array([end_f[0] * (self._energy.shape[0] - 1),
                                end_f[1] * (self._energy.shape[1] - 1),
                                end_f[2] * (self._energy.shape[2] - 1)])
                self._countimages = max(5, int(self.get_dis(self._endpoints[i], self._endpoints[i+1])//0.1))
                temp_path = self.pathfind(start, end, self._energy, n_images=self._countimages,
                                          dr=[self._struc.lattice.a / (self._energy.shape[0] - 1),
                                              self._struc.lattice.b / (self._energy.shape[1] - 1),
                                              self._struc.lattice.c / (self._energy.shape[2] - 1)],
                                          h=0.002, k=0.17, min_iter=100,max_iter=10000, max_tol=5e-6)
                if len(temp_path) > 4:
                    for p1 in temp_path:
                        p1[0] = round(p1[0] / (self._energy.shape[0] - 1), 4)
                        p1[1] = round(p1[1] / (self._energy.shape[1] - 1), 4)
                        p1[2] = round(p1[2] / (self._energy.shape[2] - 1), 4)
                if i == 0:
                    self._path = copy.deepcopy(temp_path)
                else:
                    self._path = np.append(self._path, temp_path[1:], axis=0)
            else:
                expan_energy2 = np.concatenate((self._energy, self._energy, self._energy), axis=0)
                expan_energy3 = np.concatenate((expan_energy2, expan_energy2, expan_energy2), axis=1)
                expan_energy = np.concatenate((expan_energy3, expan_energy3, expan_energy3), axis=2)
                start =np.array([start_f[0] * (self._energy.shape[0] - 1) + self._energy.shape[0],
                                 start_f[1] * (self._energy.shape[1] - 1) + self._energy.shape[1],
                                 start_f[2] * (self._energy.shape[2] - 1) + self._energy.shape[2]])
                end = np.array([end_f[0] * (self._energy.shape[0] - 1) + (1+self._phase[i][0])*self._energy.shape[0],
                                end_f[1] * (self._energy.shape[1] - 1) + (1+self._phase[i][1])*self._energy.shape[1],
                                end_f[2] * (self._energy.shape[2] - 1) + (1+self._phase[i][2])*self._energy.shape[2]])
                self._countimages = max(5, int(self.get_dis(self._endpoints[i], self._endpoints[i + 1]) // 0.2))
                temp_path = self.pathfind(start, end, expan_energy, n_images=self._countimages,
                                          dr=[self._struc.lattice.a / (self._energy.shape[0] - 1),
                                              self._struc.lattice.b / (self._energy.shape[1] - 1),
                                              self._struc.lattice.c / (self._energy.shape[2] - 1)],
                                          h=0.02, k=0.17, min_iter=100, max_iter=10000, max_tol=5e-6)
                if len(temp_path) > 4:
                    for p1 in temp_path:
                        p1[0] = round((p1[0]-self._energy.shape[0]) / (self._energy.shape[0]-1), 4)
                        p1[1] = round((p1[1]-self._energy.shape[1]) / (self._energy.shape[1]-1), 4)
                        p1[2] = round((p1[2]-self._energy.shape[2]) / (self._energy.shape[2]-1), 4)
                if i == 0:
                    self._path = copy.deepcopy(temp_path)
                else:
                    self._path = np.append(self._path, temp_path[1:], axis=0)

    @staticmethod
    def pathfind(start, end, V, n_images=21, dr=None, h=0.001, k=0.17, min_iter=100, max_iter=10000, max_tol=5e-6):
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
        # print(s)
        # Evaluate initial distances (for elastic equilibrium)
        ds0_plus = s - np.roll(s, 1, axis=0)  # 正向
        ds0_minus = s - np.roll(s, -1, axis=0)  # 负向
        ds0_plus[0] = (ds0_plus[0] - ds0_plus[0])
        ds0_minus[-1] = (ds0_minus[-1] - ds0_minus[-1])
        dV = np.gradient(V)  # 计算梯度

        # Evolve string
        for step in range(0, max_iter):
            # if step > min_iter:
            # h = h0 * np.exp(-2.0 * (step - min_iter) / max_iter)  # 逐步衰减步长以减少震荡
            # else:
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
            # if (step % 100 == 0):
            # print ("Step {} - ds = {}".format(step, tol))
        return s

    def cal_point_energy(self, point):
        """
        计算结构中任意一点的能量
        :param point: 一维列表，表示该点的分数坐标
        :return: float型，该点的能量
        """

        x = np.linspace(0, 1, self._energy.shape[0])
        y = np.linspace(0, 1, self._energy.shape[1])
        z = np.linspace(0, 1, self._energy.shape[2])
        interpolating_function = RegularGridInterpolator((x, y, z), self._energy)
        point = [math.modf(i)[0] for i in point]
        for i in range(len(point)):
            if point[i] < 0.0:
                point[i] += 1
        return interpolating_function(np.array(point))[0]

    def get_dis(self, p1, p2):
        """
        计算结构中任意两点之间的距离，该距离考虑周期性
        :param p1: 坐标，一维数组
        :param p2: 坐标，一维数组
        :return: float，距离
        """
        temp_site1 = PeriodicSite('Ar', p1, self._struc.lattice)
        temp_site2 = PeriodicSite('Ar', p2, self._struc.lattice)
        dis = temp_site1.distance(temp_site2)
        return dis

    def cal_energy(self):
        """
        计算路径上点的能量，以及瓶颈处的能量和迁移能垒
        """
        energy_path = []
        for point in self._path:
            point_energy = self.cal_point_energy(point)
            energy_path.append(point_energy)
        dis_path = [0.0]
        tol = 0
        for i in range(len(self._path) - 1):
            dist = self.get_dis(self._path[i], self._path[i + 1])
            tol += dist
            dis_path.append(tol)
        self._path_energy = np.zeros(shape=(len(dis_path), 2))
        for i in range(len(energy_path)):
            self._path_energy[i][1] = energy_path[i]
            self._path_energy[i][0] = dis_path[i]
        self._max_energy = max(energy_path)
        self._min_energy = min(energy_path)
        self._max_position = self._path[energy_path.index(self._max_energy)]
        self._barrier = self._max_energy - self._min_energy

    def path_to_poscar(self, filename):
        """
        可视化路径
        :param filename: 结构的cif文件名
        """
        struc = Structure.from_str(open(filename).read(), fmt="cif")
        struc.to(fmt='POSCAR', filename=filename.split(".")[0] + 'POSCAR')
        struc1 = Structure.from_file(filename.split(".")[0] + 'POSCAR')
        for j in range(1, len(self._path) - 1):
            struc1.insert(0, 'He', self._path[j])
        struc1.to(filename=filename.split(".")[0] + 'pathPOSCAR')


class MigrationNetwork(object):
    def __init__(self, struc, energy, voids_dict, channels_dict, filename_CIF, moveion='Li',ismergecluster=False,
                 energythreshold=None, iscalnonequalchannels = True):
        self._moveion = moveion
        self._voids = {}
        self._channels = {}
        self._energy = energy
        self._struc = struc
        self._nonequalchannels = {}
        self._mignet = None
        self._nonequl_paths = {}
        self._iscalnonequalchannels = iscalnonequalchannels
        self._energythreshold = energythreshold
        if ismergecluster:
            mergevoid = mc.MergeCluster(voids_dict, channels_dict, self._struc, self._energy, filename_CIF,
                                        clusterradii=0.5)
            for void in mergevoid.mergedvoids:
                self._voids[void.id] = void
            for channel in mergevoid.mergedchannels:
                self._channels[(channel.start, channel.end)] = channel
        else:
            self._voids = voids_dict
            self._channels = channels_dict
        self.init_voids_channels()
        self.select_voids_channels_by_energy()

    @property
    def paths_position(self):
        return [path['points_path'] for key, path in self._nonequl_paths.items()]

    @property
    def paths_energy(self):
        return [path['energys_path'] for key, path in self._nonequl_paths.items()]


    def init_voids_channels(self, radii_threadhold=0.35):
        small_voids = [void_id for void_id, void in self._voids.items() if void.radii < radii_threadhold]
        self._voids = {void_id: void for void_id, void in self._voids.items() if void.radii >= radii_threadhold}
        self._channels = {channel_id: channel for channel_id, channel in self._channels.items()
                          if channel.start not in small_voids and channel.end not in small_voids}
        self._channels = {channel_id: channel for channel_id, channel in self._channels.items()
                          if channel.radii >= radii_threadhold}

        for void_id, void in self._voids.items():
            void.energy = self.cal_point_energy(void.coord)
        for channel_id, channel in self._channels.items():
            channel.dist = self.get_dis(self._voids[channel.start].coord, self._voids[channel.end].coord)
        if self._iscalnonequalchannels:
            i = 0
            for channel_id, channel in self._channels.items():
                if self._voids[channel.start].label < self._voids[channel.end].label:
                    key = (self._voids[channel.start].label, self._voids[channel.end].label,
                           round(channel.dist, 0))
                else:
                    key = (self._voids[channel.end].label, self._voids[channel.start].label,
                           round(channel.dist, 0))
                if key not in self._nonequalchannels.keys():
                    mp = MigrationPath([self._voids[channel.start].coord, self._voids[channel.end].coord],
                                       [channel.phase], self._energy, self._struc)
                    channel.energy = mp.energy_max
                    channel.label = i
                    self._nonequalchannels[key] = {"label": i, "energy": channel.energy}
                    i += 1
                else:
                    channel.label = self._nonequalchannels[key]["label"]
                    channel.energy = self._nonequalchannels[key]["energy"]
        else:
            i = 0
            for channel_id, channel in self._channels.items():
                if self._voids[channel.start].label < self._voids[channel.end].label:
                    key = (self._voids[channel.start].label, self._voids[channel.end].label,
                           round(channel.dist, 0))
                else:
                    key = (self._voids[channel.end].label, self._voids[channel.start].label,
                           round(channel.dist, 0))
                if key not in self._nonequalchannels.keys():
                    channel.label = i
                    self._nonequalchannels[key] = {"label": i}
                    i += 1
                else:
                    channel.label = self._nonequalchannels[key]["label"]

    def select_voids_channels_by_energy(self):
        if self._energythreshold:
            if self._iscalnonequalchannels:
                highenergy_voids = [void_id for void_id, void in self._voids.items() if
                                    void.energy >= self._energythreshold]
                self._voids = {void_id: void for void_id, void in self._voids.items() if
                               void.energy < self._energythreshold}
                self._channels = {channel_id: channel for channel_id, channel in self._channels.items()
                                  if channel.start not in highenergy_voids and channel.end not in highenergy_voids}
                self._channels = {channel_id: channel for channel_id, channel in self._channels.items()
                                  if channel.energy < self._energythreshold}
            else:
                self._voids = {void_id: void for void_id, void in self._voids.items() if
                               void.energy < self._energythreshold}

    def cal_nonequalpath_between_voids(self):
        if self._iscalnonequalchannels:
            nonequalpath = []
            print(self._voids.keys())
            for channelid,channel in self._nonequalchannels.items():
                print(channelid)
                if channelid[0] in self._voids.keys() and channelid[1] in self._voids.keys():
                    print(self._voids[channelid[0]].energy)
                    path = {'start':channelid[0],'end':channelid[1],'startenergy': self._voids[channelid[0]].energy,
                            'endenergy': self._voids[channelid[1]].energy,'bnenergy':channel['energy']}
                    nonequalpath.append(path)
            return nonequalpath
        else:
            raise("please set iscalnonequalchannels parameter to True!")

    def cal_nonequl_paths(self):
        self._mignet = nx.DiGraph()
        for id_void, void in self._voids.items():
            self._mignet.add_node(void.id, label=void.label, coord=void.coord, energy=void.energy)
        for channel_id, channel in self._channels.items():  # 添加边
            # channel_barrier = max(0, channel.energy-self._voids[channel.start].energy)
            self._mignet.add_edge(channel.start, channel.end, label=channel.label, phase=channel.phase)
        self.cal_noneqpaths_cavd()
        self.cal_noneqpaths_bvse()

    def cal_point_energy(self, point):
        # 计算能量场中一点的能量
        x = np.linspace(0, 1, self._energy.shape[0])
        y = np.linspace(0, 1, self._energy.shape[1])
        z = np.linspace(0, 1, self._energy.shape[2])
        interpolating_function = RegularGridInterpolator((x, y, z), self._energy)
        point = [math.modf(i)[0] for i in point]
        for i in range(len(point)):
            if point[i] < 0.0:
                point[i] += 1
        return interpolating_function(np.array(point))[0]

    def get_dis(self, p1, p2):
        temp_site1 = PeriodicSite('Ar', p1, self._struc.lattice)
        temp_site2 = PeriodicSite('Ar', p2, self._struc.lattice)
        dis = temp_site1.distance(temp_site2)
        return dis

    def cal_noneqpaths_cavd(self):
        """
        使用cavd计算出的间隙网络计算出非等价迁移路径
        """
        labels_positions = []
        positions_moveion = []
        positions_pair = {}
        a = SpacegroupAnalyzer(self._struc, symprec=0.1)
        symm_structure = a.get_symmetrized_structure()
        for i, sites in enumerate(symm_structure.equivalent_sites):
            for site in sites:
                if site.specie.symbol == self._moveion:
                    temp_position = {}
                    temp_position['fracoord'] = site.frac_coords
                    temp_position['atom_site_label'] = site._atom_site_label
                    mini_dis = 1000
                    for void in self._voids.values():
                        dis = self.get_dis(temp_position['fracoord'], void.coord)
                        if dis < mini_dis:
                            mini_dis = dis
                            temp_position['label'] = void.label
                            temp_position['id_CAVD'] = void.id
                            if mini_dis < 0.1:
                                break
                    if temp_position['label'] not in labels_positions:
                        labels_positions.append(temp_position['label'])
                    positions_moveion.append(temp_position)
        for i in range(len(positions_moveion) - 1):
            for j in range(i + 1, len(positions_moveion)):
                if positions_moveion[i]['label'] <= positions_moveion[j]['label']:
                    key = (positions_moveion[i]['atom_site_label'], positions_moveion[j]['atom_site_label'])
                else:
                    key = (positions_moveion[j]['atom_site_label'], positions_moveion[i]['atom_site_label'])
                if positions_moveion[i]['id_CAVD'] < positions_moveion[j]['id_CAVD']:
                    value = [positions_moveion[i]['id_CAVD'], positions_moveion[j]['id_CAVD']]
                else:
                    value = [positions_moveion[j]['id_CAVD'], positions_moveion[i]['id_CAVD']]
                if key not in positions_pair.keys():
                    positions_pair[key] = [value]
                else:
                    positions_pair[key].append(value)
        for pair_label, position_pairs in positions_pair.items():
            for position_pair in position_pairs:
                try:
                    #voids_path = nx.shortest_path(self._mignet, source=position_pair[0],
                                                  # target=position_pair[1], weight='weight')
                    voids_path = list(nx.shortest_path(self._mignet,
                                                          source=position_pair[0], target=position_pair[1]))
                except nx.NetworkXNoPath:
                    voids_path = None
                if 0 < len(voids_path) < 7:
                    temppath = {'phase_path': [0, 0, 0], 'phases_path_segment':[], 'points_path':[],
                                'energys_path':[], 'pair_label': pair_label}
                    labels_path = []
                    for i in range(len(voids_path)-1):
                        # 计算路径片段的label
                        labels_path.append(self._mignet[voids_path[i]][voids_path[i+1]]['label'])
                        temppath['phases_path_segment'].append(self._mignet[voids_path[i]][voids_path[i+1]]['phase'])
                        temppath['phase_path'][0] += self._mignet[voids_path[i]][voids_path[i+1]]['phase'][0]
                        temppath['phase_path'][1] += self._mignet[voids_path[i]][voids_path[i+1]]['phase'][1]
                        temppath['phase_path'][2] += self._mignet[voids_path[i]][voids_path[i+1]]['phase'][2]
                    tag0 = True  # tag0用来判断是否含晶格位点
                    voidlabels_path = [self._mignet.node[voids_path[0]]["label"]]
                    for i in range(1, len(voids_path)-1):
                        lab = self._mignet.node[voids_path[i]]["label"]
                        if lab not in labels_positions:
                            voidlabels_path.append(lab)
                        else:
                            tag0 = False
                            break
                    voidlabels_path.append(self._mignet.node[voids_path[-1]]["label"])
                    temppath['voidlabels_path'] = voidlabels_path
                    temppath['voidid_path'] = voids_path
                    temppath['voidcoord_path'] = [self._voids[id].coord for id in voids_path]
                    temppath['barrier_path'] = None
                    if tag0:
                        print(labels_path)
                        if labels_path[0] < labels_path[-1]:
                            key_path = tuple(labels_path)
                        else:
                            key_path = tuple(labels_path[::-1])
                        if key_path not in self._nonequl_paths.keys():
                            self._nonequl_paths[key_path] = temppath
                        else:
                            if operator.eq(temppath['phase_path'], [0, 0, 0]) and \
                               not (operator.eq(self._nonequl_paths[key_path]['phase_path'], [0, 0, 0])):
                                self._nonequl_paths[tuple(labels_path)] = temppath

    def cal_noneqpaths_bvse(self):
        """
            计算非等价迁移路径上的各点坐标和能量
        """
        for path in self._nonequl_paths.values():
            coords_path = path['voidcoord_path']
            phases_path = path['phases_path_segment']
            mp = MigrationPath(coords_path, phases_path, self._energy, self._struc)
            path['points_path'] = mp.path
            path['energys_path'] = mp.path_energy
            path['barrier_path'] = mp.barrier
            print(path['pair_label'][0], path['pair_label'][1], round(path['barrier_path'],3))
            '''struc = Structure.from_str(open(self._filename_cif).read(), fmt="cif")
            struc.to(fmt='POSCAR', filename=self._filename_cif.split(".")[0] + '_POSCAR')
            struc1 = Structure.from_file(self._filename_cif.split(".")[0] + '_POSCAR')
            for j in range(1, len(path['points_path']) - 1):
                struc1.insert(0, 'H', path['points_path'][j])
            struc1.to(filename=self._filename_cif.split(".")[0] + str(i) + '_POSCARpath')
            print(path['voidlabels_path'], path['barrier_path'])
            i += 1'''

    def save_data(self,filename):
        """
        保存数据
        """
        with open(filename.split(".")[0] + '_nonequalpaths.txt', 'w') as f:
            i = 0
            for nonqul_id, path in self._nonequl_paths.items():
                if len(path['points_path']) > 0:
                    f.write('nonequalpath id:' + str(i) + '\n')
                    i += 1
                    points = path['points_path']
                    energys = path['energys_path']
                    for j in range(len(points)):
                        f.write(str(j) + "\t" + str(points[j][0]) + '  ' + str(points[j][1])+ '  '+ str(points[j][2])
                                + '\t'+ str(energys[j][0]) + '  '+ str(energys[j][1]) + '\n ')

    def showenergy(self, filename):
        for nonqul_id, path in self._nonequl_paths.items():
            if len(path['points_path']) > 0:
                xcoords = [path['energys_path'][j][0] for j in range(len(path['energys_path']))]
                ycoords = [path['energys_path'][j][1] for j in range(len(path['energys_path']))]
                ycoords = [y - min(ycoords) for y in ycoords]
                poly = np.polyfit(xcoords, ycoords, deg=7)  # 最小二乘法多项式拟合
                z = np.polyval(poly, xcoords)
                plt.figure(figsize=(8, 6))
                plt.plot(xcoords, z, linewidth=2, color="black")
                plt.xticks(fontsize=16)
                plt.yticks(fontsize=16)
                plt.scatter(xcoords, z, color='k', marker='o')
                plt.xlabel("Reaction Coordinate ", fontsize=18)  # X轴标签
                plt.ylabel("Energy", fontsize=18)  # Y轴标签
                save_path = filename.split(".")[0]+'-'.join([str(i) for i in nonqul_id]) + '.png'
                plt.savefig(save_path)
                plt.show()


def paths_to_poscar(filename_cavd, filename_cif, filename_bvse,energythreshold=None):
    voids_list = []
    channels_list = []
    flag_p = 0
    flag_n = 0
    file = open(filename_cavd, 'r')
    for line in file.readlines():
        if 'Interstitial' in line:
            flag_p = 1
            flag_n = 0
            continue
        if 'Connection' in line:
            flag_p = 0
            flag_n = 1
            continue
        if flag_p == 1:
            line = line.split()
            if len(line) > 4:
                void = Void()
                void.id = int(line[0])
                void.label = int(line[1])
                void.coord = [np.float64(line[2]), np.float64(line[3]), np.float64(line[4])]
                void.radii = np.float64(line[5])
                voids_list.append(void)
        if flag_n == 1:
            line = line.split()
            if len(line) > 5:
                channel = Channel()
                channel.start = int(line[0])
                channel.end = int(line[1])
                channel.phase = [int(line[2]), int(line[3]), int(line[4])]
                channel.coord = [np.float64(line[5]), np.float64(line[6]), np.float64(line[7])]
                channel.radii = np.float64(line[8])
                channels_list.append(channel)
    voids = {}
    channels = {}
    for void in voids_list:
        voids[void.id] = void
    for channel in channels_list:
        if channel.start < channel.end:
            print(channel.start, channel.end)
            channels[(channel.start, channel.end)] = channel
    struc = Structure.from_str(open(filename_cif).read(), fmt="cif")
    struc.to(fmt='POSCAR', filename=filename_cif.split(".")[0] + '_POSCAR')
    struc1 = Structure.from_file(filename_cif.split(".")[0] + '_POSCAR')
    path_endpoints = []
    energy = np.load(filename_bvse)
    for channel_id, channel in channels.items():
        mp = MigrationPath([voids[channel.start].coord, voids[channel.end].coord],
                           [channel.phase], energy, struc)
        path = mp.path
        energy_max = mp.energy_max
        if energy_max > energythreshold:
            print(channel.start, channel.end)
        else:
            path_endpoints.append([channel.start, channel.end])
            for j in range(1, len(path) - 1):
                struc1.insert(0, 'Ne', path[j])
    for channel in path_endpoints:
        struc1.insert(0, 'He', voids[channel[0]].coord)
        struc1.insert(0, 'He', voids[channel[1]].coord)
    struc1.to(filename=filename_cif.split(".")[0] + '_POSCAR')

#filename_CAVD = "example/LGPS/LGPS.net"
#filename_BVSE = "example/LGPS/LGPS.npy"
filename_CIF = "example/LGPS/LGPS.cif"
#paths_to_poscar(filename_CAVD, filename_CIF,filename_BVSE)


