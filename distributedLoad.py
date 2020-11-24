import matplotlib.pyplot as plt
import numpy as np
from aircraftProperties import AircraftProperties
from scipy.integrate import quad
from scipy.interpolate import interp1d

class distributedLoad:
    semispan = AircraftProperties.Planform["span"]/2

    integrationLimit = 50
    integrationFidelity = 20

    def __init__(self, limits):
        self.limits = limits
        self.distributedLoads = None
        self.pointLoads = []
        self.distributedMoments = None

        self.pointMoments = []

        self.shearForce = None

    def getShearForce(self, fidelity=integrationFidelity, limit=integrationLimit):
        if self.shearForce is not None:
            return self.shearForce

        yList = np.linspace(self.limits[0], self.limits[1], fidelity, endpoint=True)


        shear = []

        if self.distributedLoads is None:
            for yVal in yList:
                for force in self.pointLoads:
                    if force[0] > yVal:
                        shear.append(force[1])
                    else:
                        shear.append(0)
        else:
            for yVal in yList:
                i = quad(self.distributedLoads, yVal, self.limits[1], limit=limit)
                shearAtPoint = i[0]
                for force in self.pointLoads:
                    if force[0] > yVal:
                        shearAtPoint += force[1]
                shear.append(shearAtPoint)

        shearFunc = interp1d(yList, shear)

        self.shearForce = shearFunc
        return shearFunc

    def getMomentDistribution(self, fidelity=integrationFidelity, limit=integrationLimit):
        if self.shearForce is None:
            return self.distributedMoments
        if self.distributedMoments is not None:
            return self.distributedMoments

        yList = np.linspace(self.limits[0], self.limits[1], fidelity, endpoint=True)

        moment = []

        for yVal in yList:
            i = quad(self.shearForce, yVal, self.limits[1], limit=limit)
            momentAtPoint = -i[0]
            for moment in self.pointMoments:
                if moment[0] > yVal:
                    momentAtPoint += moment[1]
            moment.append(momentAtPoint)

        momentFunc = interp1d(yList, moment)

        self.distributedMoments = momentFunc

        return momentFunc

    def addPointLoad(self, magnitude, point, name=""):
        if magnitude != 0:
            self.pointLoads.append([point, magnitude, name])
        return

    def addDistibutedLoad(self, load):
        if self.distributedLoads is not None:
            self.distributedLoads = lambda y: load(y) + self.distributedLoads(y)
        else:
            self.distributedLoads = load

    def addMomentDistribution(self, moment):
        if self.distributedMoments is not None:
            self.distributedMoments = lambda y: moment(y) + self.distributedLoads(y)
        else:
            self.distributedMoments = moment

    def drawLoads(self, fidelity=50):
        yList = np.linspace(0, self.semispan, fidelity, endpoint=True)

        plt.axhline(y=0, color='k')
        if self.distributedLoads is not None:
            plt.plot(yList, [self.distributedLoads(i) for i in yList])
        plt.xlabel('semi-span [m]')
        plt.ylabel('Force [N]')


        #draw point forces
        if self.pointLoads != []:
            origin = [[], []]
            vector = [[], []]
            for load in self.pointLoads:
                origin[0].append(load[0])
                vector[0].append(0)
                if load[1] < 0:
                    origin[1].append(abs(load[1]))
                    vector[1].append(load[1])
                else:
                    origin[1].append(0)
                    vector[1].append(load[1])

                plt.annotate(load[2], [load[0], abs(load[1])])
            plt.quiver(origin[0], origin[1], vector[0], vector[1], color="red", angles='xy', scale_units='xy', scale=1, width=0.003)
            plt.xlim(self.limits)


        plt.show()

    def drawShear(self, fidelity=50):
        yList = np.linspace(0, self.semispan, fidelity, endpoint=True)
        plt.axhline(y=0, color='k')
        if self.shearForce is not None:
            plt.plot(yList, [self.shearForce(i) for i in yList])
        plt.title("Shear Force")
        plt.xlabel('semi-span [m]')
        plt.ylabel('Shear Force [N]')
        plt.xlim(self.limits)
        plt.show()

    def drawMoment(self, fidelity=50):
        yList = np.linspace(0, self.semispan, fidelity, endpoint=True)
        plt.axhline(y=0, color='k')
        if self.distributedMoments is not None:
            plt.plot(yList, [self.distributedMoments(i) for i in yList])
        plt.title("Moment")
        plt.xlabel('semi-span [m]')
        plt.ylabel('Moment [Nm]')
        plt.xlim(self.limits)
        plt.show()
