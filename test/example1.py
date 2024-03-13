import numpy as np
import json 
import networkx as nx

from Spyffness import Frame

struct = Frame()

struct.addNode(0, [0,0,0])
struct.addNode(1, [0,0,3])
struct.addNode(2, [0,3,0])

struct.addBeam(0, 0, 1)
struct.addBeam(1, 1, 2)
struct.addBeam(2, 2, 0)

struct.addMaterial(False, E= 210e9, Iy= 0.0001, Iz= 0.0001, G= 80e9, J= 0.0001, A= 0.01)

struct.fixAllBottomNodes()
struct.setComprensionLoad(1000)

struct.solve()
