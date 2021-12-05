from scipy.constants import epsilon_0, pi
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator




def calculateZ(chargeLocation, observationLocation, delta_s): #calculate potential contribution of a charged patch to an observation point on the center of a patch.

    if np.array_equal(chargeLocation, observationLocation): #if charged patch is on/near the observation patch, approximate
        Z = 1 / (2 * epsilon_0 * np.sqrt(pi * delta_s))
        return Z
    else: #use potential formula
        distanceRmRn = np.linalg.norm(chargeLocation - observationLocation) #find magnitude of distance vector
        Z = 1 / (4 * pi * epsilon_0 * distanceRmRn)
        return Z


def getPatchVectorLocations(w, d, subdivisions, numberOfPatches): #get an array of all N 3d position vectors of the centers of the patches
    delta_w = w / subdivisions #get the increment between centers of patch
    patchLocations = np.zeros((numberOfPatches, 3)) #center vector locations of all N patches
    numberOfPatchesPerPlate = int(numberOfPatches / 2)

    for i in range(numberOfPatchesPerPlate): #bottom plate where z is 0
        yPatchIndex = i % subdivisions
        xPatchIndex = np.floor(i / subdivisions)

        x = 0.5 * delta_w + xPatchIndex * delta_w
        y = 0.5 * delta_w + yPatchIndex * delta_w
        z = 0

        locationVector = np.array([x, y, z])
        patchLocations[i] = locationVector

    for i in range(numberOfPatchesPerPlate): #top plate where z is d

        yPatchIndex = i % subdivisions
        xPatchIndex = np.floor(i / subdivisions)

        x = 0.5 * delta_w + xPatchIndex * delta_w
        y = 0.5 * delta_w + yPatchIndex * delta_w
        z = d

        locationVector = np.array([x, y, z])
        patchLocations[i + numberOfPatchesPerPlate] = locationVector

    return patchLocations

        





#parameters

V_0 = 1 #positive plate is V_0 / 2
w = 10 #width of plates in mm
d = 3 #gap between plates in mm
subdivisions = 40 #how many discrete intervals along the sides of the plates. Each plate will be subdivisions x subdivisions in terms of patches

numberOfPatches = np.square(subdivisions) * 2 #calculate the total number of patches
delta_s = (2 / numberOfPatches) * np.square(w) #get delta_s discrete area

patchVectorLocations = getPatchVectorLocations(w, d, subdivisions, numberOfPatches)

Z = np.zeros((numberOfPatches, numberOfPatches)) #Z is N x N

for m in range(Z.shape[0]):

    observationLocation = patchVectorLocations[m] #every row of Z is a fixed observation point

    for n in range(Z.shape[1]): #vary where the charged plate is that contributes the potential

        Z_mn = calculateZ(patchVectorLocations[n], observationLocation, delta_s)
        Z[m][n] = Z_mn


V = np.zeros(numberOfPatches) #V is N x 1

for i in range(numberOfPatches): 
    #construct V such that the first half elements are the voltages of the bottom plate's patches, -V_0 / 2
    #second half elements are voltages of top plate's patches, V_0 / 2
    if np.floor(i / (numberOfPatches / 2)) == 0:
        V[i] = - V_0 / 2
    else:
        V[i] = V_0 / 2


Z_inv = np.linalg.inv(Z)

Q = np.matmul(Z_inv, V) #Q charge vector is found with matrix algebra

numberOfPatchesPerPlate = int(numberOfPatches / 2)

#post process the Q vector for plotting distribution of one plate

x = np.zeros(numberOfPatchesPerPlate)
y = np.zeros(numberOfPatchesPerPlate)
surfaceChargeDensity = np.zeros(numberOfPatchesPerPlate)
totalChargePerPlate = 0

for i in range(numberOfPatchesPerPlate):
    locationVector = patchVectorLocations[i + numberOfPatchesPerPlate] #offsetting by N/2 since we want the top plate's positive charges
    x[i] = locationVector[0]
    y[i] = locationVector[1] #extract the x,y coordinate of the top plate patches
    surfaceChargeDensity[i] = Q[i + numberOfPatchesPerPlate] / delta_s #extract the discrete charge densities of top plate patches by dividing by patch area
    totalChargePerPlate = totalChargePerPlate + Q[i + numberOfPatchesPerPlate] #sum up all discrete charges to find plate's total charge


calculatedCapacitance = totalChargePerPlate / V_0 #with the summed up charges of one plate, divide by V_0 to extract capacitance


#plot the top plate's charge distribution

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
fig.suptitle('Top plate surface charge distribution, no. of patches per plate: {}, width: {}mm, gap: {}mm, capacitance: {}mF'.format(numberOfPatchesPerPlate, w, d, calculatedCapacitance), fontsize=12)
surf = ax.plot_trisurf(x, y, surfaceChargeDensity, cmap=cm.viridis)
ax.set_xlabel('X (mm)')
ax.set_ylabel('Y (mm)')
ax.set_zlabel('\n\n\n\ncharge density\n(mC / mm^2)')
plt.show()











#vary d and observe calculated capacitance


d_plot = np.arange(0.5, 40, 0.1) #distances from 1 to 40mm
w = 10 #fixed 1cm width
V_0 = 1
c_calc = [] #calculated capacitances using discretized method
c_formula = []#using approximation formula
subdivisions = 5
numberOfPatches = np.square(subdivisions) * 2 
delta_s = (2 / numberOfPatches) * np.square(w)

for d in d_plot: #vary d and calculate capacitance using same discrete method

    print(d)

    patchVectorLocations = getPatchVectorLocations(w, d, subdivisions, numberOfPatches)

    Z = np.zeros((numberOfPatches, numberOfPatches))

    for m in range(Z.shape[0]):

        observationLocation = patchVectorLocations[m] 

        for n in range(Z.shape[1]):

            Z_mn = calculateZ(patchVectorLocations[n], observationLocation, delta_s)
            Z[m][n] = Z_mn


    V = np.zeros(numberOfPatches)

    for i in range(numberOfPatches):
        if np.floor(i / (numberOfPatches / 2)) == 0:
            V[i] = - V_0 / 2
        else:
            V[i] = V_0 / 2

    Z_inv = np.linalg.inv(Z)

    Q = np.matmul(Z_inv, V)

    numberOfPatchesPerPlate = int(numberOfPatches / 2)

    x = np.zeros(numberOfPatchesPerPlate)
    y = np.zeros(numberOfPatchesPerPlate)
    surfaceChargeDensity = np.zeros(numberOfPatchesPerPlate)
    totalChargePerPlate = 0

    for i in range(numberOfPatchesPerPlate):
        locationVector = patchVectorLocations[i + numberOfPatchesPerPlate]
        x[i] = locationVector[0]
        y[i] = locationVector[1]
        surfaceChargeDensity[i] = Q[i + numberOfPatchesPerPlate] / delta_s
        totalChargePerPlate = totalChargePerPlate + Q[i + numberOfPatchesPerPlate]

    capacitance_calculated = totalChargePerPlate / V_0

    c_calc.append(capacitance_calculated)

    capacitance_formula = epsilon_0 * (np.square(w * w) / d)

    c_formula.append(capacitance_formula)

    
plt.plot(d_plot, c_calc, label = "Capcitance calculated using discrete method")
plt.plot(d_plot, c_formula, label = "Capcitance calculated using approximation formula")
plt.xlabel('Gap distance (mm)')
plt.ylabel('Calculated capacitance (mF)')
plt.legend()
plt.show()





