import numpy as np
import json 
import networkx as nx

from Spyffness.Frame import Frame

struct = Frame()

struct.addNode(0, [0,0,0])
struct.addNode(1, [0,0,3])
struct.addNode(2, [0,3,0])

struct.addMaterial(False, E= 30000, Iy= 100, Iz= 100, G= 10000, J= 50, A= 10)

struct.addBeam(0, 0, 1)
struct.addBeam(1, 1, 2)
struct.addBeam(2, 0, 2)

# print(struct.K())
struct.fixAllBottomNodes()
struct.setCompresionLoad(1000)

struct.solve()

print(struct.Beams[0].Floc())
print(struct.Beams[1].Floc())
print(struct.Beams[2].Floc())