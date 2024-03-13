import numpy as np

class Node():


    def __init__(self, ID, coords):

        self.ID = ID
        self.coords = np.array(coords)

        self.Loads = np.zeros(6)
        self.Displacements = np.zeros(6)
        self.Reactions = np.zeros(6)

        self.support = np.ones(6)

    def distance(self, other):
        return np.linalg.norm(self.coords - other.coords)
    
    def fixedSupport(self):
        self.support = np.zeros(6)
    
    def pinnedSupport(self, thetax = 1, thetay = 1, thetaz = 1):
        self.support = np.array([0,0,0,thetax,thetay,thetaz])
    
