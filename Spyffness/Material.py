import config

class Material():

    def __init__(self, from_config = True, E= None, Iy= None, Iz= None, G= None, J= None, A= None):

        if from_config:
            self.E = config.E
            self.Iy = config.Iy
            self.Iz = config.Iz
            self.G = config.G
            self.J = config.J
            self.A = config.A

        else:
            self.E = E
            self.Iy = Iy
            self.Iz = Iz
            self.G = G
            self.J = J
            self.A = A
