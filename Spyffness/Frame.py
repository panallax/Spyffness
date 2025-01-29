import numpy as np
import networkx as nx
from scipy import sparse

from Spyffness.Beam import Beam
from Spyffness.Node import Node
from Spyffness.Material import Material


class Frame():

    def __init__(self):
        self.Nodes = {} # {ID: Node}
        self.Beams = {} # {ID: Beam}


    def addMaterial(self, form_config= True, E= None, Iy= None, Iz= None, G= None, J= None, A= None):
        if form_config:
            self.Material = Material(True)
        else:
            self.Material = Material(False, E, Iy, Iz, G, J, A)

    def addNode(self, ID, coords):
        if ID:
            self.Nodes[ID] = Node(ID, coords)
        else:
            self.Nodes[len(self.Nodes)] = Node(len(self.Nodes), coords)

    def addBeam(self, ID, node1, node2, A = None):
        if ID:
            self.Beams[ID] = Beam(ID, self.Nodes[node1], self.Nodes[node2], self.Material, A)
        else:
            self.Beams[len(self.Beams)] = Beam(len(self.Beams), self.Nodes[node1], self.Nodes[node2], self.Material, A)

    def newSections(self, A):
        for i, beam in enumerate(self.Beams):
            beam.A = A[i]
    
    def fromGraph(self, G, A = None):
        """
        The nodes are ordered on the z-axis. Member nodes are also ordered according to the z-coordinate.
        """
        sorted_nodes = sorted(G.nodes, key= lambda k: k[2])

        for i, node in enumerate(sorted_nodes):
            self.addNode(i, np.array(node))
        

        sorted_edges = [sorted([sorted(x, key= lambda k: k[2]) for x in G.edges()], key= lambda k: k[0][2])]
        for i, edge in enumerate(sorted_edges):
            node1= self.__getNodeIDfromCoords(edge[0])
            node2= self.__getNodeIDfromCoords(edge[1])
            self.addBeam(i, node1, node2, A)

    def __getNodeIDfromCoords(self, coords):
        for id, node in self.Nodes.items():
            if np.array_equal(node.coords, coords):
                return id 
            
    
    def K(self):
        n = len(self.Nodes)
        k = np.zeros((6*n, 6*n))
        for beam in self.Beams.values():
            kglob = beam.Kglob()
            i = beam.node1.ID
            j = beam.node2.ID
            k[6*i:6*i+6, 6*i:6*i+6] += kglob[:6,:6]
            k[6*j:6*j+6, 6*j:6*j+6] += kglob[6:,6:]
            k[6*i:6*i+6, 6*j:6*j+6] += kglob[:6,6:]
            k[6*j:6*j+6, 6*i:6*i+6] += kglob[6:,:6]

        return k
    
    def setLoads(self, node, load_vect):
        self.Nodes[node].Loads = load_vect
    
    def setSupport(self, node, support_vect):
        self.Nodes[node].support = support_vect
    
    def fixAllBottomNodes(self):
        bottom_coord = self.__findBottomNodesCoord()
        for node in self.Nodes.values():
            if node.coords[2] == bottom_coord:
                node.fixedSupport()

    def pinAllBottomNodes(self):
        for node in self.Nodes.values():
            if node.coords[2] == 0:
                node.pinnedSupport()

    def setCompresionLoad(self, load):
        top_coord = self.__findTopNodesCoord()
        for node in self.Nodes.values():
            if node.coords[2] == top_coord:
                node.Loads[2] = -load

    def __getSupports(self):
        supports = np.zeros(6*len(self.Nodes))
        for i, node in self.Nodes.items():
            supports[6*i:6*i+6] = node.support

        if np.count_nonzero(supports) == 0:
            raise ValueError("No supports have been set")
        
        return supports
    
    def __getLoads(self):
        loads = np.zeros(6*len(self.Nodes))
        for i, node in self.Nodes.items():
            loads[6*i:6*i+6] = node.Loads
        
        if np.count_nonzero(loads) == 0:
            raise ValueError("No loads have been set")
        
        return loads
    
    def __getDisplacements(self):
        displacements = np.zeros(6*len(self.Nodes))
        for i, node in self.Nodes.items():
            displacements[6*i:6*i+6] = node.Displacements

        return displacements
    
    def __getReactions(self):
        reactions = np.zeros(6*len(self.Nodes))
        for i, node in self.Nodes.items():
            reactions[6*i:6*i+6] = node.Reactions

        return reactions
    
    def solve(self):
        k = self.K()
        # print(sparse.csr_matrix(k))
        supports = self.__getSupports()
        loads = self.__getLoads()

        displacements = self.__getDisplacements()
        reactions = self.__getReactions()
        free = np.where(supports == 1)[0]
        fixed = np.where(supports == 0)[0]
        kfree = k[free,:][:,free]
        kfixed = k[fixed,:][:,fixed]

        loads = loads[free]
        displacements[free] = np.linalg.solve(kfree, loads)
        print(displacements)

        reactions[fixed] = np.dot(kfixed, displacements[fixed])
        print(reactions)
        self.__setDisplacements(displacements)
        self.__setReactions(reactions)

    def __setDisplacements(self, displacements):
        for i, node in self.Nodes.items():
            node.Displacements = displacements[6*i:6*i+6]

    
    def __setReactions(self, reactions):
        for i, node in self.Nodes.items():
            node.Reactions = reactions[6*i:6*i+6]

    def __findBottomNodesCoord(self):
        min_coord = 0
        for node in self.Nodes.values():
            if node.coords[2] < min_coord:
                min_coord = node.coords[2]

        return min_coord
    
    def __findTopNodesCoord(self):
        max_coord = 0
        for node in self.Nodes.values():
            if node.coords[2] > max_coord:
                max_coord = node.coords[2]

        return max_coord