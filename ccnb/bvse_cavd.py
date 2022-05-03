import operator
import copy
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
    Calculate the  minimum energy path between endpoints, as well as the energy of each image on
    the path, the energy of bottleneck site and the migration energy barrier
    """

    def __init__(self, endpoints, phases, energy, struc, countimages=11):
        """
        :param endpoints: two-dimensional array,migration path endpoint coordinates
        :param phase: two-dimensional array.Periodic displacement vector of path
        :param energy: three-dimensional array.Bvse potential
        :param struc: pymatgen structure object
        """
        self._energy = energy
        self._struc = struc
        self._endpoints = endpoints
        self._phase = phases  # 例如[[0,0,0],[0,-1,0]]
        self._countimages = countimages
        # 路径上每点的坐标
        self._path = []
        # 路径上每点的能量
        self._path_energy = []
        # 最高点的能量
        self._max_energy = None
        self._max_position = None  # 最高点的位置
        self._min_energy = None  # 最低点的能量
        self._barrier = None  # 最高点的能量-最低点的能量
        self._energy_function = None  # 能量计算插值函数
        self.__set_energy_function()
        self.cal_path()  # 计算路径
        self.cal_energy()  # 计算路径上每点的能量

    @property
    def path(self):
        """
        Return the coordinates of each image on the MEP
        """
        return self._path

    @property
    def path_energy(self):
        """
        Return the energy of each image on the MEP. The type is a two-dimensional array.
        The first column represents the distance of each image, and the second column represents
        the energy of each image
        """
        return self._path_energy

    @property
    def energy_max(self):
        """
        :return: BVSE energy of bottleneck site
        """
        return self._max_energy

    @property
    def position_max(self):
        """
        :return: the coordinate of  bottleneck site
        """
        return self._max_position

    @property
    def energy_min(self):
        """
        :return:  the BVSE energy of the image with the lowest BVSE energy value
        """
        return self._min_energy

    @property
    def barrier(self):
        """
        :return: migration barrier
        """
        return self._barrier

    def cal_path(self):
        """
        Calculate the MEP between endpoints
        :return: Two dimensional array，every row is the coordinate of each image on the MEP
        """
        energy_shape = np.array(self._energy.shape)
        for i in range(len(self._endpoints) - 1):
            start_f = [j for j in self._endpoints[i]]
            start_f = np.mod(start_f, 1)
            end_f = [j for j in self._endpoints[i + 1]]
            end_f = np.mod(end_f, 1)
            if operator.eq(self._phase[i], [0, 0, 0]):
                start = start_f * (energy_shape - 1)
                end = end_f * (energy_shape - 1)
                self._countimages = max(
                    4,
                    int(
                        self.get_dis(self._endpoints[i],
                                     self._endpoints[i + 1]) // 0.2))
                temp_path = self.pathfind(
                    start,
                    end,
                    self._energy,
                    n_images=self._countimages,
                    dr=[
                        self._struc.lattice.a / (self._energy.shape[0] - 1),
                        self._struc.lattice.b / (self._energy.shape[1] - 1),
                        self._struc.lattice.c / (self._energy.shape[2] - 1)
                    ],
                    h=0.002,
                    k=0.17,
                    min_iter=100,
                    max_iter=10000,
                    max_tol=5e-6)
                if len(temp_path) > 3:
                    temp_path = np.round(temp_path / (energy_shape - 1), 4)
                    temp_path = np.mod(temp_path, 1.0)
                if len(self._path):
                    self._path = np.append(self._path, temp_path[1:], axis=0)
                else:
                    self._path = copy.deepcopy(temp_path)
            else:
                expan_energy2 = np.concatenate(
                    (self._energy, self._energy, self._energy), axis=0)
                expan_energy3 = np.concatenate(
                    (expan_energy2, expan_energy2, expan_energy2), axis=1)
                expan_energy = np.concatenate(
                    (expan_energy3, expan_energy3, expan_energy3), axis=2)
                start = np.array([
                    start_f[0] * (self._energy.shape[0] - 1) +
                    self._energy.shape[0], start_f[1] *
                    (self._energy.shape[1] - 1) + self._energy.shape[1],
                    start_f[2] * (self._energy.shape[2] - 1) +
                    self._energy.shape[2]
                ])
                end = end_f * (energy_shape - 1) + (
                    1 + np.array(self._phase[i])) * energy_shape
                self._countimages = max(
                    4,
                    int(
                        self.get_dis(self._endpoints[i],
                                     self._endpoints[i + 1]) // 0.2))
                temp_path = self.pathfind(
                    start,
                    end,
                    expan_energy,
                    n_images=self._countimages,
                    dr=[
                        self._struc.lattice.a / (self._energy.shape[0] - 1),
                        self._struc.lattice.b / (self._energy.shape[1] - 1),
                        self._struc.lattice.c / (self._energy.shape[2] - 1)
                    ],
                    h=0.02,
                    k=0.17,
                    min_iter=100,
                    max_iter=10000,
                    max_tol=5e-6)
                if len(temp_path) > 3:
                    temp_path = np.round(
                        (temp_path - energy_shape) / (energy_shape - 1), 4)
                    temp_path = np.mod(temp_path, 1.0)
                if len(self._path):
                    self._path = np.append(self._path, temp_path[1:], axis=0)
                else:
                    self._path = copy.deepcopy(temp_path)

    @staticmethod
    def pathfind(start,
                 end,
                 V,
                 n_images=21,
                 dr=None,
                 h=0.001,
                 k=0.17,
                 min_iter=100,
                 max_iter=10000,
                 max_tol=5e-6):
        # created by pymatgen
        if not dr:
            dr = np.array([
                1.0 / (V.shape[0] - 1), 1.0 / (V.shape[1] - 1),
                1.0 / (V.shape[2] - 1)
            ])
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

        # Evlve string
        for step in range(0, max_iter):
            # if step > min_iter:
            # h = h0 * np.exp(-2.0 * (step - min_iter) / max_iter)  # 逐步衰减步长以减少震荡
            # else:
            h = h0
            # Calculate forces acting on string
            d = V.shape
            s0 = s
            edV = np.array([[
                dV[0][int(pt[0]) % d[0]][int(pt[1]) % d[1]][int(pt[2]) % d[2]]
                / dr[0],
                dV[1][int(pt[0]) % d[0]][int(pt[1]) % d[1]][int(pt[2]) % d[2]]
                / dr[1],
                dV[2][int(pt[0]) % d[0]][int(pt[1]) % d[1]][int(pt[2]) % d[2]]
                / dr[2]
            ] for pt in s])
            # Update according to force due to potential and string elasticity
            ds_plus = s - np.roll(s, 1, axis=0)
            ds_minus = s - np.roll(s, -1, axis=0)
            ds_plus[0] = (ds_plus[0] - ds_plus[0])
            ds_minus[-1] = (ds_minus[-1] - ds_minus[-1])
            Fpot = edV
            Fel = keff * (la.norm(ds_plus) -
                          la.norm(ds0_plus)) * (ds_plus / la.norm(ds_plus))
            Fel += keff * (la.norm(ds_minus) -
                           la.norm(ds0_minus)) * (ds_minus / la.norm(ds_minus))
            s = s - h * (Fpot + Fel)
            # Keep the coordinates of the endpoints unchanged in each iteration
            s[0] = s0[0]
            s[-1] = s0[-1]
            # Reparametrize string
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

    def __set_energy_function(self):
        x = np.linspace(0, 1, self._energy.shape[0])
        y = np.linspace(0, 1, self._energy.shape[1])
        z = np.linspace(0, 1, self._energy.shape[2])
        self._energy_function = RegularGridInterpolator((x, y, z),
                                                        self._energy)

    def cal_point_energy(self, point):
        """
        Calculate the BVSE energy according to the discrete grid coordinates of the site point
        :param point: array，the discrete grid coordinates of the site point
        :return: float，BVSE energy
        """
        return self._energy_function(np.array(np.mod(point, 1.0)))[0]

    def get_dis(self, p1, p2):
        """
        Calculate the distance between p1 and p2 in the structure, which considers periodicity
        :param p1: One-dimensional array such as [0.11,0.22,0.33]
        :param p2:，One-dimensional
        :return: float，the distance between p1 and p2
        """
        temp_site1 = PeriodicSite('Ar', p1, self._struc.lattice)
        temp_site2 = PeriodicSite('Ar', p2, self._struc.lattice)
        dis = temp_site1.distance(temp_site2)
        return dis

    def cal_energy(self):
        """
        Calculate the BVSE energy of points on the MEP, the BVSE energy at the bottleneck and energy barriers
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
        Visual MEP by Poscar file
        :param filename: cif filename
        """
        struc = Structure.from_str(open(filename).read(), fmt="cif")
        struc.to(fmt='POSCAR', filename=filename.split(".")[0] + 'POSCAR')
        struc1 = Structure.from_file(filename.split(".")[0] + 'POSCAR')
        for j in range(1, len(self._path) - 1):
            struc1.insert(0, 'He', self._path[j])
        struc1.to(filename=filename.split(".")[0] + 'pathPOSCAR')


class MigrationNetwork(object):

    def __init__(self,
                 struc,
                 energy,
                 voids_dict,
                 channels_dict,
                 filename_CIF,
                 moveion='Li',
                 ismergecluster=False,
                 energythreshold=None,
                 iscalnonequalchannels=True,
                 clusterradii=0.75,
                 void_radii_threshold=[0.5, 1.0]):
        """
        Calculate the transport network,analysis transport pathways
        """
        self._moveion = moveion  # 迁移离子
        self._voids = {}  # CAVD计算的所有间隙点
        self._channels = {}  # CAVD计算的所有通道段
        self._energy = energy  # BVSE三维势场
        self._struc = struc  # Pymatgen读取的结构
        self._nonequalchannels = {}  # 以间隙为起始点的所有非等价路径
        self._mignet = None  # 间隙网络图 networkx
        self._nonequl_paths = {}  # 以晶格位为起始点的所有非等价路径
        # true、false 用于判断是否计算以间隙为起始点的所有非等价路径
        self._iscalnonequalchannels = iscalnonequalchannels
        # 用于筛选间隙和通道段的阈值
        if energythreshold:
            self._energythreshold = energythreshold + np.min(energy)
        if ismergecluster:  # 选择是否合并间隙簇
            mergevoid = mc.MergeCluster(voids_dict,
                                        channels_dict,
                                        self._struc,
                                        self._energy,
                                        filename_CIF,
                                        clusterradii=clusterradii)  # 调用合并间隙簇类
            mergevoid.to_net(filename_CIF)
            for void in mergevoid.mergedvoids:
                self._voids[void.id] = void
            for channel in mergevoid.mergedchannels:
                self._channels[(channel.start, channel.end)] = channel
        else:
            self._voids = voids_dict
            self._channels = channels_dict
        self.init_voids_channels()

    @property
    def paths_position(self):  # 返回以晶格位为起始点的所有非等价路径上每一点的坐标
        return [
            path['points_path'] for key, path in self._nonequl_paths.items()
        ]

    @property
    def paths_energy(self):  #返回以晶格位为起始点的所有非等价路径上每一点的能量
        return [
            path['energys_path'] for key, path in self._nonequl_paths.items()
        ]

    def init_voids_channels(self):
        for void_id, void in self._voids.items():
            void.energy = self.cal_point_energy(void.coord)  # 计算每一个间隙点的能量
        for channel_id, channel in self._channels.items():  # 计算每一条通道段的长度
            channel.dist = self.get_dis(self._voids[channel.start].coord,
                                        self._voids[channel.end].coord)
        if self._iscalnonequalchannels:  # 是否对每一条间隙之间的通道段计算BVSE MEP
            self.cal_nonequal_mep_between_voids()
        else:
            i = 0
            for channel_id, channel in self._channels.items(
            ):  # 计算每一条通道段的label
                if self._voids[channel.start].label < self._voids[
                        channel.end].label:
                    key = (self._voids[channel.start].label,
                           self._voids[channel.end].label,
                           round(channel.dist, 0))  # 起始点label 和长度
                else:
                    key = (self._voids[channel.end].label,
                           self._voids[channel.start].label,
                           round(channel.dist, 0))
                if key not in self._nonequalchannels.keys():
                    channel.label = i
                    self._nonequalchannels[key] = {"label": i}
                    i += 1
                else:
                    channel.label = self._nonequalchannels[key]["label"]
        if self._energythreshold:
            highenergy_voids = [
                void_id for void_id, void in self._voids.items()
                if void.energy >= self._energythreshold
            ]
            self._voids = {
                void_id: void
                for void_id, void in self._voids.items()
                if void.energy < self._energythreshold
            }
            if self._iscalnonequalchannels:
                self._channels = {
                    channel_id: channel
                    for channel_id, channel in self._channels.items()
                    if channel.start not in highenergy_voids
                    and channel.end not in highenergy_voids
                }
                self._channels = {
                    channel_id: channel
                    for channel_id, channel in self._channels.items()
                    if channel.energy < self._energythreshold
                }
            else:
                self._channels = {
                    channel_id: channel
                    for channel_id, channel in self._channels.items()
                    if channel.start not in highenergy_voids
                    and channel.end not in highenergy_voids
                }

    def cal_nonequal_mep_between_voids(self):
        # 对每一条间隙之间的通道段计算BVSE MEP
        i = 0
        for channel_id, channel in self._channels.items():
            if self._voids[channel.start].label < self._voids[
                    channel.end].label:
                key = (self._voids[channel.start].label,
                       self._voids[channel.end].label, round(channel.dist, 0))
            else:
                key = (self._voids[channel.end].label,
                       self._voids[channel.start].label,
                       round(channel.dist, 0))
            if key not in self._nonequalchannels.keys():
                mp = MigrationPath([
                    self._voids[channel.start].coord,
                    self._voids[channel.end].coord
                ], [channel.phase], self._energy, self._struc)
                channel.energy = mp.energy_max
                channel.label = i
                self._nonequalchannels[key] = {
                    "start": channel.start,
                    "end": channel.end,
                    "label": i,
                    "energy": channel.energy
                }
                i += 1
            else:
                channel.label = self._nonequalchannels[key]["label"]
                channel.energy = self._nonequalchannels[key]["energy"]
        nonequalpath = []
        for channelid, channel in self._nonequalchannels.items():
            path = {
                'start': channel["start"],
                'end': channel["end"],
                'startenergy': self._voids[channel["start"]].energy,
                'endenergy': self._voids[channel["end"]].energy,
                'bnenergy': channel['energy']
            }
            nonequalpath.append(path)
            print(path)
        return nonequalpath

    def cal_nonequl_paths(self):
        self._mignet = nx.DiGraph()
        for id_void, void in self._voids.items():
            self._mignet.add_node(void.id,
                                  label=void.label,
                                  coord=void.coord,
                                  energy=void.energy)
        for channel_id, channel in self._channels.items():
            self._mignet.add_edge(channel.start,
                                  channel.end,
                                  label=channel.label,
                                  phase=channel.phase)
        self.cal_noneqpaths_cavd()
        self.cal_noneqpaths_bvse()

    def cal_point_energy(self, point):
        # Calculate the bvse energy  according to the discrete grid coordinates of the site
        x = np.linspace(0, 1, self._energy.shape[0])
        y = np.linspace(0, 1, self._energy.shape[1])
        z = np.linspace(0, 1, self._energy.shape[2])
        interpolating_function = RegularGridInterpolator((x, y, z),
                                                         self._energy)
        point = np.mod(point, 1.0)
        return interpolating_function(np.array(point))[0]

    def get_dis(self, p1, p2):
        temp_site1 = PeriodicSite('Ar', p1, self._struc.lattice)
        temp_site2 = PeriodicSite('Ar', p2, self._struc.lattice)
        dis = temp_site1.distance(temp_site2)
        return dis

    def cal_noneqpaths_cavd(self):
        """
        calculate the label of each transport pathway
        between mobile ion lattice site
        """
        # 保存所有晶格位的类型  比如 Li1,Li2，Li3
        labels_positions = []
        # 保存所有晶格位在间隙网络中的编号，label，坐标
        positions_moveion = []
        # 字典，对所有的晶格位两两配对 key是晶格位的类型
        # 比如 （Li1,Li2） value是在间隙网络中的编号 比如（5，10）
        positions_pair = {}
        a = SpacegroupAnalyzer(self._struc, symprec=0.1)
        symm_structure = a.get_symmetrized_structure()
        # 这个循环计算positions_moveion
        # 晶格位在间隙网络中的编号，label，坐标
        for i, sites in enumerate(symm_structure.equivalent_sites):
            for site in sites:
                if site.specie.symbol == self._moveion:
                    temp_position = {}
                    temp_position['fracoord'] = site.frac_coords
                    temp_position['atom_site_label'] = site._atom_site_label
                    mini_dis = 1000
                    for void in self._voids.values():
                        dis = self.get_dis(temp_position['fracoord'],
                                           void.coord)
                        if dis < mini_dis:
                            mini_dis = dis
                            temp_position['label'] = void.label
                            temp_position['id_CAVD'] = void.id
                            if mini_dis < 0.05:
                                break
                    if temp_position['label'] not in labels_positions:
                        labels_positions.append(temp_position['label'])
                    positions_moveion.append(temp_position)
        # 这个循环计算positions_pair
        for i in range(len(positions_moveion) - 1):
            for j in range(i + 1, len(positions_moveion)):
                if positions_moveion[i]['label'] <= positions_moveion[j][
                        'label']:
                    key = (positions_moveion[i]['atom_site_label'],
                           positions_moveion[j]['atom_site_label'])
                else:
                    key = (positions_moveion[j]['atom_site_label'],
                           positions_moveion[i]['atom_site_label'])
                if positions_moveion[i]['id_CAVD'] < positions_moveion[j][
                        'id_CAVD']:
                    value = [
                        positions_moveion[i]['id_CAVD'],
                        positions_moveion[j]['id_CAVD']
                    ]
                else:
                    value = [
                        positions_moveion[j]['id_CAVD'],
                        positions_moveion[i]['id_CAVD']
                    ]
                if key not in positions_pair.keys():
                    positions_pair[key] = [value]
                else:
                    positions_pair[key].append(value)
        for pair_label, position_pairs in positions_pair.items():
            for position_pair in position_pairs:
                try:
                    voids_path = list(
                        nx.shortest_path(self._mignet,
                                         source=position_pair[0],
                                         target=position_pair[1]))
                # 计算在间隙网络中的最短路径，
                # 比如Li-Li2 结果为[1,5,7,9,2] 每一个数字表示间隙网络中的编号
                except nx.NetworkXNoPath:
                    # 如果没有路径即为空，不考虑
                    voids_path = []
                if 0 < len(voids_path) < 7:
                    temppath = {
                        'phase_path': [0, 0, 0],
                        'phases_path_segment': [],
                        'points_path': [],
                        'energys_path': [],
                        'pair_label': pair_label
                    }
                    labels_path = []  # 在该晶格位之间的每一通道段的label
                    for i in range(len(voids_path) - 1):
                        # phases_path_segment表示每一段的label，phase_path表示整条路径的label
                        labels_path.append(self._mignet[voids_path[i]][
                            voids_path[i + 1]]['label'])
                        temppath['phases_path_segment'].append(self._mignet[
                            voids_path[i]][voids_path[i + 1]]['phase'])
                        temppath['phase_path'][0] += self._mignet[
                            voids_path[i]][voids_path[i + 1]]['phase'][0]
                        temppath['phase_path'][1] += self._mignet[
                            voids_path[i]][voids_path[i + 1]]['phase'][1]
                        temppath['phase_path'][2] += self._mignet[
                            voids_path[i]][voids_path[i + 1]]['phase'][2]
                    tag0 = True  # 用于判断这条路径上每一间隙点是否为晶格位
                    voidlabels_path = [
                        self._mignet.node[voids_path[0]]["label"]
                    ]
                    for i in range(1, len(voids_path) - 1):
                        lab = self._mignet.node[voids_path[i]]["label"]
                        if lab not in labels_positions:
                            voidlabels_path.append(lab)
                        else:
                            # 如果是就不计算
                            tag0 = False
                            break
                    voidlabels_path.append(
                        self._mignet.node[voids_path[-1]]["label"])
                    # 路径上每一间隙点的label
                    temppath['voidlabels_path'] = voidlabels_path
                    # 路径上每一间隙点的编号
                    temppath['voidid_path'] = voids_path
                    # 路径上每一间隙点的坐标
                    temppath['voidcoord_path'] = [
                        self._voids[id].coord for id in voids_path
                    ]
                    temppath['barrier_path'] = None
                    if tag0 and len(labels_path) > 0:
                        # 对起始点排序，编号小的在前，
                        if labels_path[0] < labels_path[-1]:
                            key_path = tuple(labels_path)
                        else:
                            key_path = tuple(labels_path[::-1])
                        if key_path not in self._nonequl_paths.keys(
                        ):  # 判断是否在非等价路径集合里
                            self._nonequl_paths[key_path] = temppath
                        else:
                            # 和他同一label的路径
                            if operator.eq(temppath['phase_path'],
                                           [0, 0, 0]) and not (operator.eq(
                                               self._nonequl_paths[key_path]
                                               ['phase_path'], [0, 0, 0])):
                                # 如果这条路径有和他同一label的路径，但是这条和他同一label的路径是跨包的，
                                # 这条路径不跨包，就用它代替，尽可能保证选出的路径是不跨包的
                                self._nonequl_paths[tuple(
                                    labels_path)] = temppath

    def cal_noneqpaths_bvse(self):
        """
            Calculate the mep under BVSE potential
        """
        for path in self._nonequl_paths.values():
            coords_path = path['voidcoord_path']
            phases_path = path['phases_path_segment']
            mp = MigrationPath(coords_path, phases_path, self._energy,
                               self._struc)
            path['points_path'] = mp.path
            path['energys_path'] = mp.path_energy
            path['barrier_path'] = mp.barrier

    def save_data(self, filename):
        """
        Save nonequivalent transport pathways to ***_nonequalpaths.txt
        """
        with open(filename.split(".")[0] + '_nonequalpaths.txt', 'w') as f:
            i = 0
            for nonqul_id, path in self._nonequl_paths.items():
                if len(path['points_path']) > 0:
                    f.write('nonequalpath id:' + str(i) + '\n')
                    i += 1
                    f.write(path['pair_label'][0] + "  ->  " +
                            path['pair_label'][1] + "  energy barrier  " +
                            str(round(path['barrier_path'], 3)) + '\n ')
                    points = path['points_path']
                    energys = path['energys_path']
                    for j in range(len(points)):
                        f.write(
                            str(j) + "\t" + str(points[j][0]) + '  ' +
                            str(points[j][1]) + '  ' + str(points[j][2]) +
                            '\t' + str(energys[j][0]) + '  ' +
                            str(energys[j][1]) + '\n ')

    def showenergy(self, filename):
        index = 0
        for nonqul_id, path in self._nonequl_paths.items():
            if len(path['points_path']) > 0:
                xcoords = [
                    path['energys_path'][j][0]
                    for j in range(len(path['energys_path']))
                ]
                ycoords = [
                    path['energys_path'][j][1]
                    for j in range(len(path['energys_path']))
                ]
                # ycoords = [y - min(ycoords) for y in ycoords]
                poly = np.polyfit(xcoords, ycoords, deg=9)
                z = np.polyval(poly, xcoords)
                plt.figure(figsize=(8, 6))
                plt.plot(xcoords, z, linewidth=2, color="black")
                plt.xticks(fontsize=16)
                plt.yticks(fontsize=16)
                plt.scatter(xcoords, ycoords, color='k', marker='o')
                plt.xlabel("Reaction Coordinate ", fontsize=18)
                plt.ylabel("Energy", fontsize=18)
                save_path = filename.split(".")[0]+"_"+str(index) \
                            + '_' + path['pair_label'][0] \
                            + '-' + path['pair_label'][1] + '_'\
                            + '-'.join([str(i) for i in nonqul_id]) \
                            + '.svg'
                index = index + 1
                plt.savefig(save_path)
                plt.show()


def paths_to_poscar(filename_cavd,
                    filename_cif,
                    filename_bvse,
                    energythreshold=None):
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
                void.coord = [
                    np.float64(line[2]),
                    np.float64(line[3]),
                    np.float64(line[4])
                ]
                void.radii = np.float64(line[5])
                voids_list.append(void)
        if flag_n == 1:
            line = line.split()
            if len(line) > 5:
                channel = Channel()
                channel.start = int(line[0])
                channel.end = int(line[1])
                channel.phase = [int(line[2]), int(line[3]), int(line[4])]
                channel.coord = [
                    np.float64(line[5]),
                    np.float64(line[6]),
                    np.float64(line[7])
                ]
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
        mp = MigrationPath(
            [voids[channel.start].coord, voids[channel.end].coord],
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
