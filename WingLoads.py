import matplotlib.pyplot as plt
import numpy as np
from aircraftProperties import AircraftProperties
from scipy.integrate import quad
from scipy.interpolate import interp1d
from mpl_toolkits.mplot3d import Axes3D

class WingLoads:
    # wing distributed load with origin at wingbox centroid at wingroot, y spanwise along centroid (body to tip), z normal (upwards), x tangential (aftwards)
    __semispan = AircraftProperties.Planform["span"]/2
    __indeces= {
        "x": 0,
        "y": 1,
        "z": 2,
    }
    __externalForcesNames = [r"External Forces x", r"External Forces y", r"External Forces z"]
    __externalMomentNames = [r"External Moments x", r"External Torque", r"External Moments z"]
    __internalMomentNames = [r"$M_x$", r"Torque", r"$M_z$"]
    __internalForceNames = [r"$V_x$", r"Invalid", r"$V_z$"]

    integrationLimit = 50
    integrationFidelity = 20

    def __init__(self, pointForces, distributedForces, distributedMoments):
        self.limits = [0, self.__semispan]

        self.pointForces = pointForces

        self.distributedForces = [lambda y: 0, lambda y: 0, lambda y: 0]
        for i in range(3):
            if distributedForces[i] is not None:
                self.distributedForces[i] = distributedForces[i]

        self.distributedMoments = [lambda y: 0, lambda y: 0, lambda y: 0]
        for i in range(3):
            if distributedMoments[i] is not None:
                self.distributedMoments[i] = distributedMoments[i]

        #pos, magnitude
        self.perAxisPointForces = self.decomposePointLoads(self.pointForces)
        self.perAxisPointMoments = self.getPointMoments(self.pointForces)

        self.__normalForce = None
        self.__internalShear = [None, lambda y: 0, None] # transverse shear at y axis is 0
        self.__internalMoments = [None, None, None]

    def decomposePointLoads(self, pointForces):
        loads = [[], [], []]

        for force in pointForces:
            for axis in range(3):
                if abs(force.forceVector[axis]) > 0.1: # same as for point moments
                    loads[axis].append([force.forceVector[axis], force.originPoint[1], force.identifier])

        return loads

    def getPointMoments(self, pointForces):
        moments = [[], [], []]

        for force in pointForces:
            y = force.originPoint[1]
            pointMoments = force.getMomentsAroundPoint([0, y, 0])
            for i in range(3):
                if abs(pointMoments[i]) > 0.1: # since rotating a point force often transforms 0 values to smth like e-16 (smth similar to floating point inacuracy happening), therefore moments smaller than 0.1 are discarded since their contribtuions are extremely small anyway (0.1Nm<<100KNm)
                    moments[i].append([pointMoments[i], y, force.identifier])

        return moments

    def getNormalForce(self, fidelity=integrationFidelity, limit=integrationLimit):
        if self.__normalForce is not None:
            return self.__normalForce

        yList = np.linspace(self.limits[0], self.limits[1], fidelity, endpoint=True)

        normalForce = []

        for yVal in yList:
            i = quad(self.distributedForces[1], yVal, self.limits[1], limit=limit)
            normalAtPoint = i[0]
            for force in self.perAxisPointForces[1]:
                if force[1] > yVal:
                    normalAtPoint += force[0]
            normalForce.append(normalAtPoint)

        normalFunc = interp1d(yList, normalForce)

        self.__normalForce = normalFunc
        return normalFunc

    def getShearForce(self, axis, fidelity=integrationFidelity, limit=integrationLimit):
        if self.__internalShear[axis] is not None:
            return self.__internalShear[axis]

        yList = np.linspace(self.limits[0], self.limits[1], fidelity, endpoint=True)

        internalForce = []

        for yVal in yList:
            i = quad(self.distributedForces[axis], yVal, self.limits[1], limit=limit)
            shearAtPoint = i[0]
            for force in self.perAxisPointForces[axis]:
                if force[1] > yVal:
                    shearAtPoint += force[0]
            internalForce.append(shearAtPoint)

        shearFunc = interp1d(yList, internalForce)

        self.__internalShear[axis] = shearFunc
        return shearFunc

    def getInternalMoment(self, axis, fidelity=integrationFidelity, limit=integrationLimit):
        if self.__internalMoments[axis] is not None:
            return self.__internalMoments[axis]

        bendingAxis = [2, 1, 0]

        if self.__internalShear[bendingAxis[axis]] is None:
            self.getShearForce(bendingAxis[axis])

        yList = np.linspace(self.limits[0], self.limits[1], fidelity, endpoint=True)

        internalMoment = []

        for yVal in yList:
            # contribtuion from shear force
            i = quad(self.__internalShear[bendingAxis[axis]], yVal, self.limits[1], limit=limit)
            momentAtPoint = -i[0]
            # contribution from distributed moment
            j = quad(self.distributedMoments[axis], yVal, self.limits[1], limit=limit)
            momentAtPoint -= j[0]

            #contribution from point moments created by point forces
            for moment in self.perAxisPointMoments[axis]:
                if moment[1] > yVal:
                    momentAtPoint -= moment[0]

            internalMoment.append(momentAtPoint)

        if axis == 1:
            momentFunc = interp1d(yList, [-i for i in internalMoment])
        else:
            momentFunc = interp1d(yList, internalMoment)

        self.__internalMoments[axis] = momentFunc
        return momentFunc

    def drawExternalForces(self, axis, fidelity=50):
        yList = np.linspace(0, self.__semispan, fidelity, endpoint=True)
        fig, plots = plt.subplots(1)
        plots.axhline(y=0, color='k')
        plots.set_xlim(self.limits)
        plots.plot(yList, [self.distributedForces[axis](i) for i in yList])
        plots.set_title(self.__externalForcesNames[axis])
        plots.set(xlabel='semi-span [m]', ylabel='force [N]')

        origin = [[], []]
        vector = [[], []]
        for load in self.perAxisPointForces[axis]:
            origin[0].append(load[1])
            vector[0].append(0)
            if load[0] < 0:
                origin[1].append(abs(load[0]))
                vector[1].append(load[0])
            else:
                origin[1].append(0)
                vector[1].append(load[0])

            plots.annotate(load[2], [load[1], abs(load[0])])
            s = [0 for n in range(len(origin[0]))]
            plots.scatter(load[1], abs(load[0]), s=s)
        plots.quiver(origin[0], origin[1], vector[0], vector[1], color="red", angles='xy', scale_units='xy', scale=1, width=0.003)
        plt.show()

    def drawExternalMoments(self, axis, fidelity=50):
        yList = np.linspace(0, self.__semispan, fidelity, endpoint=True)
        fig, plots = plt.subplots(1)
        plots.axhline(y=0, color='k')
        plots.set_xlim(self.limits)
        plots.plot(yList, [self.distributedMoments[axis](i) for i in yList])
        plots.set_title(self.__externalMomentNames[axis])
        plots.set(xlabel='semi-span [m]', ylabel='moment [Nm]')

        origin = [[], []]
        vector = [[], []]
        for load in self.perAxisPointMoments[axis]:
            plots.annotate(("Moment " + load[2] + ": " + str(round(load[0], 1)) + " [Nm]"), [load[1], 0], color="red")
            plots.scatter(load[1], 0, color="red")
        plots.quiver(origin[0], origin[1], vector[0], vector[1], color="red", angles='xy', scale_units='xy', scale=1, width=0.003)
        plt.show()

    def drawShearForce(self, axis, fidelity=50):
        internalForce = self.getShearForce(axis)
        yList = np.linspace(0, self.__semispan, fidelity, endpoint=True)
        fig, plots = plt.subplots(1)
        plots.axhline(y=0, color='k')
        plots.set_xlim(self.limits)
        plots.plot(yList, [internalForce(i) for i in yList])
        plots.set_title(self.__internalForceNames[axis])
        plots.set(xlabel='semi-span [m]', ylabel='[N]')

        origin = [[], []]
        vector = [[], []]

        plots.quiver(origin[0], origin[1], vector[0], vector[1], color="red", angles='xy', scale_units='xy', scale=1, width=0.003)
        plt.show()

    def drawInternalMoment(self, axis, fidelity=50):
        internalMoment = self.getInternalMoment(axis)
        yList = np.linspace(0, self.__semispan, fidelity, endpoint=True)
        fig, plots = plt.subplots(1)
        plots.axhline(y=0, color='k')
        plots.set_xlim(self.limits)
        plots.plot(yList, [internalMoment(i) for i in yList])
        plots.set_title(self.__internalMomentNames[axis])
        plots.set(xlabel='semi-span [m]', ylabel='[Nm]')

        origin = [[], []]
        vector = [[], []]

        plots.quiver(origin[0], origin[1], vector[0], vector[1], color="red", angles='xy', scale_units='xy', scale=1, width=0.003)
        plt.show()

    def drawNormalForce(self, fidelity=50):
        internalForce = self.getNormalForce()
        yList = np.linspace(0, self.__semispan, fidelity, endpoint=True)
        fig, plots = plt.subplots(1)
        plots.axhline(y=0, color='k')
        plots.set_xlim(self.limits)
        plots.plot(yList, [internalForce(i) for i in yList])
        plots.set_title("Normal Force")
        plots.set(xlabel='semi-span [m]', ylabel='[N]')

        origin = [[], []]
        vector = [[], []]

        plots.quiver(origin[0], origin[1], vector[0], vector[1], color="red", angles='xy', scale_units='xy', scale=1, width=0.003)
        plt.show()



class PointForce:
    def __init__(self, originPoint, forceVector, identifier=""):
        self.originPoint = originPoint
        self.forceVector = forceVector
        self.identifier = identifier

    #get moments around point
    def getMomentsAroundPoint(self, point):
        r = [self.originPoint[i] - point[i] for i in range(3)]
        #rxF
        return np.cross(r, self.forceVector)

    def getMagnitude(self):
        magnitudeSquared = 0

        for i in self.forceVector:
            magnitudeSquared += i**2

        return magnitudeSquared**0.5

    #rotate by angle around axis through origin https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    def getRotatedVectorAroundAxis(self, theta, uvecAxis):
        cos = np.cos(theta)
        sin = np.sin(theta)
        rotationMatrix = np.array([ [cos+uvecAxis[0]**2*(1-cos), uvecAxis[0]*uvecAxis[1]*(1-cos)-uvecAxis[2]*sin, uvecAxis[0]*uvecAxis[2]*(1-cos)+uvecAxis[1]*sin],
                                    [uvecAxis[1]*uvecAxis[0]*(1-cos)+uvecAxis[2]*sin, cos+uvecAxis[1]**2*(1-cos), uvecAxis[1]*uvecAxis[2]*(1-cos)-uvecAxis[0]*sin],
                                    [uvecAxis[2]*uvecAxis[0]*(1-cos)-uvecAxis[1]*sin, uvecAxis[2]*uvecAxis[1]*(1-cos)+uvecAxis[0]*sin, cos+uvecAxis[2]**2*(1-cos)]])

        return PointForce(self.originPoint.copy(), rotationMatrix.dot(self.forceVector), self.identifier)

    def draw(self):
        mag = self.getMagnitude()
        vect = [(i/mag)*5 for i in self.forceVector]
        figure = plt.figure()
        plot = figure.add_subplot(111, projection='3d')
        plot.quiver(*self.originPoint, *vect, color="red")

        plot.plot([-2.5, -2.5, 2.5, 2.5, -2.5], [0, 29, 29, 0, 0], [0, 0, 0, 0, 0])# draw rough planform
        plot.quiver(-2.5, 14.5, 0, -5, 0, 0, color="blue")

        plot.set_title(self.identifier)

        plot.set_xlim(-15, 15)
        plot.set_ylim(-0.5, 30)
        plot.set_zlim(-15, 15)
        plt.show()
