import matplotlib.pyplot as plt


class Visualization():

    def __init__(self, model):
        self.Beams = model.Beams
        self.Nodes = model.Nodes


    def plot(self):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for beam in self.Beams.values():
            x = [beam.node1.coords[0], beam.node2.coords[0]]
            y = [beam.node1.coords[1], beam.node2.coords[1]]
            z = [beam.node1.coords[2], beam.node2.coords[2]]
            ax.plot(x, y, z)
        plt.show()

    def plotDisplacements(self, scale = 1):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for beam in self.Beams.values():
            x = [beam.node1.coords[0] + scale*beam.node1.Displacements[0], beam.node2.coords[0] + scale*beam.node2.Displacements[0]]
            y = [beam.node1.coords[1] + scale*beam.node1.Displacements[1], beam.node2.coords[1] + scale*beam.node2.Displacements[1]]
            z = [beam.node1.coords[2] + scale*beam.node1.Displacements[2], beam.node2.coords[2] + scale*beam.node2.Displacements[2]]
            ax.plot(x, y, z)
        plt.show()

    def plotReactions(self, scale = 1):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for beam in self.Beams.values():
            x = [beam.node1.coords[0] - scale*beam.node1.Reactions[0], beam.node2.coords[0] - scale*beam.node2.Reactions[0]]
            y = [beam.node1.coords[1] - scale*beam.node1.Reactions[1], beam.node2.coords[1] - scale*beam.node2.Reactions[1]]
            z = [beam.node1.coords[2] - scale*beam.node1.Reactions[2], beam.node2.coords[2] - scale*beam.node2.Reactions[2]]
            ax.plot(x, y, z)
        plt.show()

    def plotLoads(self, scale = 1):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for beam in self.Beams.values():
            x = [beam.node1.coords[0] + scale*beam.node1.Loads[0], beam.node2.coords[0] + scale*beam.node2.Loads[0]]
            y = [beam.node1.coords[1] + scale*beam.node1.Loads[1], beam.node2.coords[1] + scale*beam.node2.Loads[1]]
            z = [beam.node1.coords[2] + scale*beam.node1.Loads[2], beam.node2.coords[2] + scale*beam.node2.Loads[2]]
            ax.plot(x, y, z)
        plt.show()

    def plotSupports(self, scale = 1):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for beam in self.Beams.values():
            x = [beam.node1.coords[0] + scale*beam.node1.support[0], beam.node2.coords[0] + scale*beam.node2.support[0]]
            y = [beam.node1.coords[1] + scale*beam.node1.support[1], beam.node2.coords[1] + scale*beam.node2.support[1]]
            z = [beam.node1.coords[2] + scale*beam.node1.support[2], beam.node2.coords[2] + scale*beam.node2.support[2]]
            ax.plot(x, y, z)
        plt.show()

    