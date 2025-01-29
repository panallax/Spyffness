import numpy as np

class Beam():

    def __init__(self, ID, node1, node2, material, A= None):

        self.ID = ID
        self.node1 = node1
        self.node2 = node2

        self.E = material.E
        self.Iy = material.Iy
        self.Iz = material.Iz
        self.G = material.G
        self.J = material.J

        if not A:
            self.A = material.A

        self.Ta = self.T()
    
    def L(self):
        return self.node1.distance(self.node2)
    
    def Kloc(self):
        """
        Local stiffness matrix
        """
        A = self.A
        E = self.E
        Iz = self.Iz
        Iy = self.Iy
        G = self.G
        J = self.J
        L = self.L()

        k = np.array([[A*E/L,  0,             0,             0,      0,            0,            -A*E/L, 0,             0,             0,      0,            0],
                    [0,      12*E*Iz/L**3,  0,             0,      0,            6*E*Iz/L**2,  0,      -12*E*Iz/L**3, 0,             0,      0,            6*E*Iz/L**2],
                    [0,      0,             12*E*Iy/L**3,  0,      -6*E*Iy/L**2, 0,            0,      0,             -12*E*Iy/L**3, 0,      -6*E*Iy/L**2, 0],
                    [0,      0,             0,             G*J/L,  0,            0,            0,      0,             0,             -G*J/L, 0,            0],
                    [0,      0,             -6*E*Iy/L**2,  0,      4*E*Iy/L,     0,            0,      0,             6*E*Iy/L**2,   0,      2*E*Iy/L,     0],
                    [0,      6*E*Iz/L**2,   0,             0,      0,            4*E*Iz/L,     0,      -6*E*Iz/L**2,  0,             0,      0,            2*E*Iz/L],
                    [-A*E/L, 0,             0,             0,      0,            0,            A*E/L,  0,             0,             0,      0,            0],
                    [0,      -12*E*Iz/L**3, 0,             0,      0,            -6*E*Iz/L**2, 0,      12*E*Iz/L**3,  0,             0,      0,            -6*E*Iz/L**2],
                    [0,      0,             -12*E*Iy/L**3, 0,      6*E*Iy/L**2,  0,            0,      0,             12*E*Iy/L**3,  0,      6*E*Iy/L**2,  0],
                    [0,      0,             0,             -G*J/L, 0,            0,            0,      0,             0,             G*J/L,  0,            0],
                    [0,      0,             -6*E*Iy/L**2,  0,      2*E*Iy/L,     0,            0,      0,             6*E*Iy/L**2,   0,      4*E*Iy/L,     0],
                    [0,      6*E*Iz/L**2,   0,             0,      0,            2*E*Iz/L,     0,      -6*E*Iz/L**2,  0,             0,      0,            4*E*Iz/L]])

        return k
    
 
    def T(self):
        """
        Transformation matrix
        """
        x1, y1, z1 = self.node1.coords
        x2, y2, z2 = self.node2.coords
        L = self.L()

        x = [(x2-x1)/L, (y2-y1)/L, (z2-z1)/L]
        if np.isclose(x1, x2) and np.isclose(z1, z2):
            if y2 > y1:
                y = [-1, 0, 0]
                z = [0, 0, 1]
            else:
                y = [1, 0, 0]
                z = [0, 0, 1]

        elif np.isclose(y1, y2):
            y = [0, 1, 0]
            z = np.cross(x, y)
            z = np.divide(z, np.sqrt(z[0]**2 + z[1]**2 + z[2]**2))

        else:
            proj = [x2-x1, 0, z2-z1]

            if y2 > y1:
                z = np.cross(proj, x)
            else:
                z = np.cross(x, proj)

            z = np.divide(z, np.sqrt(z[0]**2 + z[1]**2 + z[2]**2))
            y = np.cross(z, x)
            y = np.divide(y, np.sqrt(y[0]**2 + y[1]**2 + y[2]**2))

        dirCos = np.array([x, y, z])

        T = np.zeros((12, 12))
        T[0:3, 0:3] = dirCos
        T[3:6, 3:6] = dirCos
        T[6:9, 6:9] = dirCos
        T[9:12, 9:12] = dirCos

        return T
    
    def Kglob(self):
        """
        Global stiffness matrix
        """
        return np.matmul(np.matmul(np.linalg.inv(self.T()), self.Kloc()), self.T())

    
    def nodalLoadVector(self):
        """
        Nodal load vector in global coordinates
        """
        return np.matmul(np.matmul(self.Kloc(), np.linalg.inv(self.Ta)), self.D())
    
    def D(self):
        """
        Displacement vector in global coordinates
        """
        return np.append(self.node1.Displacements, self.node2.Displacements)
    
    def Dloc(self):
        """
        Displacement vector in local coordinates
        """
        return np.matmul(np.linalg.inv(self.Ta), self.D())
    
    def Floc(self):
        """
        Force vector in local coordinates
        """
        return np.matmul(self.Kloc(), self.Dloc())