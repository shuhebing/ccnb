import operator
import math
import numpy as np
import networkx as nx
from scipy.interpolate import RegularGridInterpolator
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.sites import PeriodicSite
import ccnb.mergecluster as mc
from ccnb.bvse_cavd import MigrationPath


class MigrationNetworkLatticeSites(object):

    def __init__(self,
                 struc,
                 energy,
                 voids_dict,
                 channels_dict,
                 filename_CIF,
                 moveion='Li',
                 ismergecluster=False):
        self._moveion = moveion
        self._voids = {}
        self._channels = {}
        self._energy = energy
        self._struc = struc
        self._mignet = None
        self._sites = {}
        self._nonequalchannels = {}
        self._nonequl_paths = {}
        self.all_paths = {}
        self._id_barrier_phase = {}
        if ismergecluster:
            mergevoid = mc.MergeCluster(voids_dict,
                                        channels_dict,
                                        self._struc,
                                        self._energy,
                                        filename_CIF,
                                        clusterradii=0.5)
            mergevoid.to_net(filename_CIF)
            for void in mergevoid.mergedvoids:
                self._voids[void.id] = void
            for channel in mergevoid.mergedchannels:
                self._channels[(channel.start, channel.end)] = channel
        else:
            self._voids = voids_dict
            self._channels = channels_dict
        self.init_voids_channels()

    def init_voids_channels(self):
        for void_id, void in self._voids.items():
            void.energy = self.cal_point_energy(void.coord)
        for channel_id, channel in self._channels.items():
            channel.dist = self.get_dis(self._voids[channel.start].coord,
                                        self._voids[channel.end].coord)
        i = 0
        for channel_id, channel in self._channels.items():
            if self._voids[channel.start].label < self._voids[
                    channel.end].label:
                key = (self._voids[channel.start].label,
                       self._voids[channel.end].label, round(channel.dist, 1))
            else:
                key = (self._voids[channel.end].label,
                       self._voids[channel.start].label,
                       round(channel.dist, 1))
            if key not in self._nonequalchannels.keys():
                channel.label = i
                self._nonequalchannels[key] = {"label": i}
                i += 1
            else:
                channel.label = self._nonequalchannels[key]["label"]

    def cal_paths(self):
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
        self.cal_paths_cavd()
        self.cal_noneqpaths_bvse()
        for key_path, path in self.all_paths.items():
            phase2 = [0, 0, 0]
            for i in range(3):
                if path['phase_path'][i] != 0:
                    phase2[i] = -1 * path['phase_path'][i]
            if path['voidlabels_path'][0] < path['voidlabels_path'][-1]:
                self._id_barrier_phase[(path['voidid_path'][0], path['voidid_path'][-1])] = \
                    [self._nonequl_paths[path['key_path']]['barrier_begin_to_end'], path['phase_path']]
                self._id_barrier_phase[(path['voidid_path'][-1], path['voidid_path'][0])] = \
                    [self._nonequl_paths[path['key_path']]['barrier_end_to_begin'], phase2]
            else:
                self._id_barrier_phase[(path['voidid_path'][0], path['voidid_path'][-1])] = \
                    [self._nonequl_paths[path['key_path']]['barrier_end_to_begin'], path['phase_path']]
                self._id_barrier_phase[(path['voidid_path'][-1], path['voidid_path'][0])] = \
                    [self._nonequl_paths[path['key_path']]['barrier_begin_to_end'], phase2]
        for site in self._sites.values():
            site["neighbor"] = []
            site['barriers'] = []
            site['phases'] = []
            for endpoints, barrier_phase in self._id_barrier_phase.items():
                if endpoints[0] == site['id_CAVD']:
                    site["neighbor"].append(endpoints[1])
                    site['barriers'].append(barrier_phase[0])
                    site['phases'].append(barrier_phase[1])

    def fac2cart(self, coord):
        """
        分数坐标转换成笛卡尔坐标
        """
        return np.dot(coord, self._struc.lattice.matrix)

    def cart2fac(self, coord):
        """
        笛卡尔坐标转换成分数坐标
        """
        return np.dot(coord, np.linalg.inv(self._struc.lattice.matrix))

    def cal_point_energy(self, point):
        # Calculate the bvse energy  according to the discrete grid coordinates of the site
        x = np.linspace(0, 1, self._energy.shape[0])
        y = np.linspace(0, 1, self._energy.shape[1])
        z = np.linspace(0, 1, self._energy.shape[2])
        interpolating_function = RegularGridInterpolator((x, y, z),
                                                         self._energy)
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

    def cal_paths_cavd(self):
        """
        calculate the label of each transport pathway between mobile ion lattice site
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
                    temp_position['atom_site_label'] = site._atom_site_label
                    mini_dis = 1000
                    for void in self._voids.values():
                        dis = self.get_dis(site.frac_coords, void.coord)
                        if dis < mini_dis:
                            mini_dis = dis
                            temp_position['label'] = void.label
                            temp_position['id_CAVD'] = void.id
                            if mini_dis < 0.01:
                                break
                    temp_position['fracoord'] = self._voids[
                        temp_position['id_CAVD']].coord
                    temp_position['cartcoord'] = self.fac2cart(
                        temp_position['fracoord'])
                    if temp_position['label'] not in labels_positions:
                        labels_positions.append(temp_position['label'])
                    positions_moveion.append(temp_position)
                    print(temp_position)
                    self._sites[temp_position['id_CAVD']] = temp_position
        for i in range(len(positions_moveion) - 1):
            for j in range(i + 1, len(positions_moveion)):
                key = (positions_moveion[j]['atom_site_label'],
                       positions_moveion[i]['atom_site_label'])
                if positions_moveion[i]['label'] < positions_moveion[j][
                        'label']:
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
                except nx.NetworkXNoPath:
                    voids_path = []
                if 0 < len(voids_path) < 7:
                    temppath = {
                        'phase_path': [0, 0, 0],
                        'phases_path_segment': [],
                        'points_path': [],
                        'energys_path': [],
                        'pair_label': pair_label
                    }
                    labels_path = []
                    for i in range(len(voids_path) - 1):
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
                    tag0 = True
                    voidlabels_path = [
                        self._mignet.node[voids_path[0]]["label"]
                    ]
                    for i in range(1, len(voids_path) - 1):
                        lab = self._mignet.node[voids_path[i]]["label"]
                        if lab not in labels_positions:
                            voidlabels_path.append(lab)
                        else:
                            tag0 = False
                            break
                    voidlabels_path.append(
                        self._mignet.node[voids_path[-1]]["label"])
                    temppath['voidlabels_path'] = voidlabels_path
                    temppath['voidid_path'] = voids_path
                    temppath['voidcoord_path'] = [
                        self._voids[id].coord for id in voids_path
                    ]
                    temppath['barrier_path'] = None
                    if tag0:
                        if labels_path[0] < labels_path[-1]:
                            key_path = tuple(labels_path)
                        else:
                            key_path = tuple(labels_path[::-1])
                        temppath['key_path'] = key_path
                        self.all_paths[(voids_path[0],
                                        voids_path[-1])] = temppath
                        if key_path not in self._nonequl_paths.keys():
                            self._nonequl_paths[key_path] = temppath
                        else:
                            if operator.eq(temppath['phase_path'], [0, 0, 0]) and \
                               not (operator.eq(self._nonequl_paths[key_path]['phase_path'], [0, 0, 0])):
                                self._nonequl_paths[tuple(
                                    labels_path)] = temppath

    def cal_noneqpaths_bvse(self):
        for path in self._nonequl_paths.values():
            coords_path = path['voidcoord_path']
            phases_path = path['phases_path_segment']
            mp = MigrationPath(coords_path, phases_path, self._energy,
                               self._struc)
            path['points_path'] = mp.path
            path['energys_path'] = mp.path_energy
            energys = []
            for i in range(len(path['energys_path'])):
                energys.append(path['energys_path'][i][1])
            i_max = 0
            for i in range(len(energys)):
                if energys[i] == mp.energy_max:
                    i_max = i
                    break
            if path['voidlabels_path'][0] < path['voidlabels_path'][0]:
                path['barrier_begin_to_end'] = mp.energy_max - min(
                    energys[:i_max])
                path['barrier_end_to_begin'] = mp.energy_max - min(
                    energys[i_max:])
            else:
                path['barrier_begin_to_end'] = mp.energy_max - min(
                    energys[i_max:])
                path['barrier_end_to_begin'] = mp.energy_max - min(
                    energys[:i_max])

    def save_data(self, filename):
        """
        Save nonequivalent transport pathways to ***_nonequalpaths.txt
        """
        with open(filename.split(".")[0] + '_mc.txt', 'w') as f:
            f.write('Lattice_basis_vector ')
            f.write('\n')
            f.write(str(self._struc.lattice))
            f.write('\n\n')
            for site in self._sites.values():
                print(site)
                f.write('Site_number ' + str(site['id_CAVD']) + '\n')
                f.write('Cartesian_coordinates ' + str(site['cartcoord'][0]) +
                        ' ' + str(site['cartcoord'][1]) + ' ' +
                        str(site['cartcoord'][2]) + '\n')
                f.write('Normalized_coordinates ' + str(site['fracoord'][0]) +
                        ' ' + str(site['fracoord'][1]) + ' ' +
                        str(site['fracoord'][2]) + '\n')
                f.write('Neighbor_sites')
                for nb in site['neighbor']:
                    f.write(' ' + str(nb))
                f.write('\n')
                f.write('Neighbor_sites_phase')
                for ph in site['phases']:
                    f.write(' ' + str(ph[0]) + ' ' + str(ph[1]) + ' ' +
                            str(ph[2]) + ',')
                f.write('\n')
                f.write('Barrier')
                for ba in site['barriers']:
                    f.write(' ' + str(ba))
                f.write('\n')
                f.write('Site_types ' + str(site['atom_site_label']) + '\n')
                f.write('\n')
