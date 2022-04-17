from monty.io import zopen
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import cKDTree
import networkx as nx
import ase.spacegroup as spg
from pymatgen.core.sites import PeriodicSite
from cavd.local_environment import CifParser_new


class Void(object):

    def __init__(self):
        self.id = None
        self.label = None
        self.coord = None
        self.radii = None
        self.energy = None


class Channel(object):

    def __init__(self):
        self.start = None
        self.end = None
        self.phase = None
        self.coord = None
        self.radii = None
        self.dist = None
        self.label = None
        self.energy = None


class MergeCluster(object):

    def __init__(self,
                 voids_dict,
                 channels_dict,
                 structure,
                 bvse,
                 filename_cif,
                 clusterradii=0.5):
        self._struc = structure
        self._energy = bvse
        self._voids = voids_dict
        self._channels = channels_dict
        self._mergedvoids = []
        self._mergedchannels = []
        self._clusterradii = clusterradii
        self._clusters = []
        self.init_voids_channels()
        self.find_clusters()
        self.handle_voidsandchannels()
        self.cal_void_label(filename_cif)

    @property
    def mergedvoids(self):
        """
        Return to all merged voids
        """
        return self._mergedvoids

    @property
    def mergedchannels(self):
        """
        Return to all merged channel segments
        """
        return self._mergedchannels

    def get_absolute_dis(self, p1, p2):
        """
        Calculate the distance between two sites without considering periodicity
        :param p1: fractional coordinates of site p1，such as[0.5，0.5，0.5]
        :param p2: fractional coordinates of site p2
        :return: float, the distance between two sites
        """
        coord1 = np.array(self.fac2cart(p1))
        coord2 = np.array(self.fac2cart(p2))
        diff = np.subtract(coord1, coord2)
        return np.linalg.norm(diff)

    def get_period_dis(self, p1, p2):
        """
        Calculate the distance between two sites with considering periodicity
        :param p1: fractional coordinates of site p1，such as[0.5，0.5，0.5]
        :param p2: fractional coordinates of site p2
        :return:  the distance between two sites
        """
        temp_site1 = PeriodicSite('Ar', p1, self._struc.lattice)
        temp_site2 = PeriodicSite('Ar', p2, self._struc.lattice)
        dis = temp_site1.distance(temp_site2)
        return dis

    def fac2cart(self, coord):
        """
        Converting fractional coordinates to Cartesian coordinates
        """
        return np.dot(coord, self._struc.lattice.matrix)

    def cart2fac(self, coord):
        """
        Converting Cartesian coordinates to fractional coordinates
        """
        return np.dot(coord, np.linalg.inv(self._struc.lattice.matrix))

    def cal_point_energy(self, point):
        """
        Calculate the bond valence site energy of a site in crystal structure
        :param point: fractional coordinates of a site
        :return: float,the bond valence site energy of a site
        """
        x = np.linspace(0, 1, self._energy.shape[0])
        y = np.linspace(0, 1, self._energy.shape[1])
        z = np.linspace(0, 1, self._energy.shape[2])
        interpolating_function = RegularGridInterpolator((x, y, z),
                                                         self._energy)
        for i in range(len(point)):
            if point[i] < 0.0:
                point[i] += 1
            if point[i] > 1.0:
                point[i] -= 1
        return interpolating_function(np.array(point))[0]

    def init_voids_channels(self, radii_threadhold=0.5):
        """
        Initialize interstices and channel segments,delete too small interstices and channel segments
        :param radii_threadhold: float, filter threshold
        """
        for void_id, void in self._voids.items():
            void.energy = self.cal_point_energy(void.coord)
        small_voids = [
            void_id for void_id, void in self._voids.items()
            if void.radii < radii_threadhold
        ]
        self._voids = {
            void_id: void
            for void_id, void in self._voids.items()
            if void.radii >= radii_threadhold
        }
        self._channels = {
            channel_id: channel
            for channel_id, channel in self._channels.items() if
            channel.start not in small_voids and channel.end not in small_voids
        }
        self._channels = {
            channel_id: channel
            for channel_id, channel in self._channels.items()
            if channel.radii >= radii_threadhold
        }

    def find_clusters(self):
        """
        find each clusters in interstital network
        """
        coords = []
        pair_voids = []
        voids_list = [void for void_id, void in self._voids.items()]
        for void in voids_list:
            coords.append(self.fac2cart(void.coord))
        coord_tree = cKDTree(coords)
        all_pair_voids = [
            i for i in coord_tree.query_pairs(r=self._clusterradii)
        ]
        for i in all_pair_voids:
            pair_voids.append([voids_list[i[0]].id, voids_list[i[1]].id])
        if len(pair_voids) > 0:
            graph_clusters = nx.Graph()
            for e in pair_voids:
                graph_clusters.add_edge(e[0], e[1])
            queue_clusters = []
            for sc in nx.connected_component_subgraphs(graph_clusters):
                queue_clusters.append(list(sc.nodes))
            while queue_clusters:
                subv_in = queue_clusters.pop(0)
                subv_out = []
                temp_coord = [0.0, 0.0, 0.0]
                for subv in subv_in:
                    temp_coord[0] += self.fac2cart(self._voids[subv].coord)[0]
                    temp_coord[1] += self.fac2cart(self._voids[subv].coord)[1]
                    temp_coord[2] += self.fac2cart(self._voids[subv].coord)[2]
                centre_coord = self.cart2fac([
                    temp_coord[0] / len(subv_in), temp_coord[1] / len(subv_in),
                    temp_coord[2] / len(subv_in)
                ])
                for i in range(len(subv_in)):
                    if self.get_period_dis(
                            self._voids[subv_in[i]].coord,
                            centre_coord) > self._clusterradii + 0.2:
                        subv_out.append(subv_in[i])
                for subv in subv_out:
                    subv_in.remove(subv)
                self._clusters.append(subv_in)
                if len(subv_out) > 1:
                    pair_subvout = []
                    for i in range(len(subv_out)):
                        for j in range(i + 1, len(subv_out)):
                            if self.get_period_dis(
                                    self._voids[subv_out[i]].coord,
                                    self._voids[subv_out[j]].coord
                            ) < self._clusterradii:
                                pair_subvout.append([subv_out[i], subv_out[j]])
                    if len(pair_subvout) > 0:
                        graph_subvout = nx.Graph()
                        for e in pair_subvout:
                            graph_subvout.add_edge(e[0], e[1])
                        for sc in nx.connected_component_subgraphs(
                                graph_subvout):
                            queue_clusters.append(list(sc.nodes))
        print("Number of clusters", len(self._clusters))
        print(self._clusters)

    def handle_voidsandchannels(self):
        """
        handle clusters
        """
        mignet = nx.Graph()
        for void_id, void in self._voids.items():
            mignet.add_node(void.id,
                            label=void.label,
                            coord=void.coord,
                            radii=void.radii)
        for e_id, e in self._channels.items():
            if e.start < e.end:
                phase1 = e.phase
                phase2 = [0, 0, 0]
                for i in range(3):
                    if phase1[i] != 0:
                        phase2[i] = -1 * phase1[i]
                mignet.add_edge(e.start,
                                e.end,
                                phase1=phase1,
                                phase2=phase2,
                                coord1=e.coord,
                                radii1=e.radii,
                                dist1=e.dist)
        if len(self._clusters) > 0:
            for i in range(len(self._clusters)):
                minenergy = 20000
                centervoid_id = self._clusters[i][0]
                for void_id in self._clusters[i]:
                    if minenergy > self._voids[void_id].energy:
                        centervoid_id = void_id
                        minenergy = self._voids[void_id].energy
                center_void = Void()
                center_void.id = centervoid_id
                center_void.label = self._voids[centervoid_id].label
                center_void.coord = self._voids[centervoid_id].coord
                center_void.radii = self._voids[centervoid_id].radii
                center_void.energy = self._voids[centervoid_id].energy
                tempedges = []
                nearvoids = []
                for nearvoid in list(mignet.adj[centervoid_id].keys()):
                    if nearvoid not in self._clusters[i]:
                        nearvoids.append(nearvoid)
                        if centervoid_id < nearvoid:
                            start = centervoid_id
                            end = nearvoid
                        else:
                            end = centervoid_id
                            start = nearvoid
                        tempedges.append({
                            "from":
                            start,
                            "to":
                            end,
                            "phase1":
                            mignet[centervoid_id][nearvoid]['phase1'],
                            "phase2":
                            mignet[centervoid_id][nearvoid]['phase2'],
                            "coord1":
                            mignet[centervoid_id][nearvoid]['coord1'],
                            "radii1":
                            mignet[centervoid_id][nearvoid]['radii1'],
                            "dist1":
                            mignet[centervoid_id][nearvoid]['dist1']
                        })
                for id in self._clusters[i]:
                    if id != center_void.id:
                        for nearvoid in list(mignet.adj[id].keys()):
                            if nearvoid not in self._clusters[
                                    i] and nearvoid not in nearvoids:
                                if centervoid_id < nearvoid:
                                    start = centervoid_id
                                    end = nearvoid
                                else:
                                    end = centervoid_id
                                    start = nearvoid
                                if id < nearvoid:
                                    ph_cen_nearvoid = mignet[id][nearvoid][
                                        'phase1']
                                    ph_nearvoid_cen = mignet[id][nearvoid][
                                        'phase2']
                                else:
                                    ph_cen_nearvoid = mignet[id][nearvoid][
                                        'phase2']
                                    ph_nearvoid_cen = mignet[id][nearvoid][
                                        'phase1']
                                if centervoid_id < nearvoid:
                                    ph1 = ph_cen_nearvoid
                                    ph2 = ph_nearvoid_cen
                                else:
                                    ph2 = ph_cen_nearvoid
                                    ph1 = ph_nearvoid_cen
                                tempedges.append({
                                    "from":
                                    start,
                                    "to":
                                    end,
                                    "phase1":
                                    ph1,
                                    "phase2":
                                    ph2,
                                    "coord1":
                                    mignet[id][nearvoid]['coord1'],
                                    "radii1":
                                    mignet[id][nearvoid]['radii1'],
                                    "dist1":
                                    mignet[id][nearvoid]['dist1']
                                })
                for void in self._clusters[i]:
                    mignet.remove_node(void)
                mignet.add_node(center_void.id,
                                label=center_void.label,
                                coord=center_void.coord,
                                radii=center_void.radii)
                for e in tempedges:
                    mignet.add_edge(e['from'],
                                    e['to'],
                                    phase1=e['phase1'],
                                    phase2=e['phase2'],
                                    coord1=e['coord1'],
                                    radii1=e['radii1'],
                                    dist1=e['dist1'])
        for nd in mignet.nodes():
            tempvoid = Void()
            tempvoid.id = nd
            tempvoid.label = mignet.nodes[nd]['label']
            tempvoid.coord = mignet.nodes[nd]['coord']
            tempvoid.radii = mignet.nodes[nd]['radii']
            self._mergedvoids.append(tempvoid)
        for edge in mignet.edges():
            tempchannel1 = Channel()
            tempchannel2 = Channel()
            tempchannel1.start = edge[0]
            tempchannel1.end = edge[1]
            tempchannel2.end = edge[0]
            tempchannel2.start = edge[1]
            if edge[0] < edge[1]:
                tempchannel1.phase = mignet[edge[0]][edge[1]]["phase1"]
                tempchannel2.phase = mignet[edge[0]][edge[1]]["phase2"]
            else:
                tempchannel1.phase = mignet[edge[0]][edge[1]]["phase2"]
                tempchannel2.phase = mignet[edge[0]][edge[1]]["phase1"]
            tempchannel1.coord = mignet[edge[0]][edge[1]]["coord1"]
            tempchannel2.coord = mignet[edge[0]][edge[1]]["coord1"]
            tempchannel1.radii = mignet[edge[0]][edge[1]]["radii1"]
            tempchannel2.radii = mignet[edge[0]][edge[1]]["radii1"]
            tempchannel1.dist = mignet[edge[0]][edge[1]]["dist1"]
            tempchannel2.dist = mignet[edge[0]][edge[1]]["dist1"]
            self._mergedchannels.append(tempchannel1)
            self._mergedchannels.append(tempchannel2)

    @staticmethod
    def tag_sites(sitesym, scaled_positions, symprec=1e-5):
        scaled = np.around(np.array(scaled_positions, ndmin=2), 8)
        scaled %= 1.0
        scaled %= 1.0
        np.set_printoptions(suppress=True)
        tags = -np.ones((len(scaled), ), dtype=int)
        tagdis = 100 * np.ones((len(scaled), ), dtype=float)
        rot, trans = spg.spacegroup.parse_sitesym(sitesym)
        siteskdTree = cKDTree(scaled)
        for i in range(len(scaled)):
            if tags[i] == -1:
                curpos = scaled[i]
                sympos = np.dot(rot, curpos) + trans
                sympos %= 1.0
                sympos %= 1.0
                sympos = np.unique(np.around(sympos, 8), axis=0)
                min_dis, min_ids = siteskdTree.query(sympos, k=1)
                select = min_dis < symprec
                select_ids = min_ids[select]
                tags[select_ids] = i
                tagdis[select_ids] = min_dis[select]
        return tags, tagdis

    def cal_void_label(self, filename_cif, symprec=0.1):
        """
        Classification of interstices based on symmetry
        """
        with zopen(filename_cif, "rt") as f:
            input_string = f.read()
        parser = CifParser_new.from_string(input_string)
        sitesym = parser.get_sym_opt()
        voids_positions = []
        for void in self.mergedvoids:
            voids_positions.append(void.coord)
        tags, tagdis = self.tag_sites(sitesym, voids_positions, symprec)
        for i in range(len(tags)):
            self.mergedvoids[i].label = tags[i]

    def to_net(self, filename):
        """
        Save the merged interstices and channel fragments to a net file
        :param filename: output will be written to a file
        """
        with open(filename.split(".")[0] + '_mergecluster.net', 'w') as f:
            f.write('Interstitial table:\n')
            for void in self.mergedvoids:
                f.write(
                    str(void.id) + "\t" + str(void.label) + "\t " +
                    str(void.coord[0]) + " " + str(void.coord[1]) + " " +
                    str(void.coord[2]) + "\t " + str(void.radii) + "\n")
            f.write('Connection table:\n')
            for channel in self.mergedchannels:
                f.write(
                    str(channel.start) + "\t " + str(channel.end) + "\t " +
                    str(int(channel.phase[0])) + " " +
                    str(int(channel.phase[1])) + " " +
                    str(int(channel.phase[2])) + "\t " +
                    str(channel.coord[0]) + " " + str(channel.coord[1]) + " " +
                    str(channel.coord[2]) + "\t " + str(channel.radii) + "\n")

    def get_clusters(self):
        """
        Find all interstice clusters in a periodic unit cell
        :return: list, all interstice clusters
        """
        clusters = []
        for cluster in self._clusters:
            cluster_temp = []
            for void_id in cluster:
                cluster_temp.append({
                    "void_id": void_id,
                    "void_coord": self._voids[void_id].coord
                })
            clusters.append(cluster_temp)
        return clusters
