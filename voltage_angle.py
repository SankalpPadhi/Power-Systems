# Computation of node voltages, totalrealpowerloss and total reactive power loss in the given distribution network
# Importing modules
from matplotlib import pyplot as plt
from matplotlib import style
style.use('ggplot')
import numpy as np
import pandas as pd

# Defining input matrices
D1 = np.array([[1,2,3,4,5,2,7,4],[2,3,4,5,6,7,8,9]])
D2 = np.array([[2,3,4,5,6,7,8,9],[3,4,5,6,9,0,0,0],[4,5,6,9,0,0,0,0],[5,6,0,0,0,0,0,0],[6,0,0,0,0,0,0,0],[7,8,0,0,0,0,0,0],[8,0,0,0,0,0,0,0],[9,0,0,0,0,0,0,0]])

# Defining parameters
NB = 9
LN = 8
maxIter = 2
tol = 0.00001
MVABase = 100

# Initialising vectors
V = np.ones(NB, dtype=complex)
N = np.array([8,5,4,2,1,2,1,1])
Il = np.zeros(LN, dtype=complex)
Ib = np.zeros(LN, dtype=complex)
S = np.array([0.004+0.003j,0.005+0.004j,0.002+0.004j,0.004+0.000j,0.004+0.000j,0.001+0.000j,0.001+0.002j,0.001+0.000j])
Z = np.array([0.2 + 0.3721j, 0.3388+ 0.1j, 0+0.75196j, 0.1+0.4j,0.4+0.3j,0.5+0.6j,0.35+0.4j,0.4+0.2j])

# Computation of node voltages
for t in range(maxIter):

    def idx(x):
        return x-2
    def idx1(x):
        return x-1
    def loadCurrent(S,V,Il):
        for i in range(LN):
            Il[i] = np.conj(S[i])/np.conj(V[i+1])
        return Il
    def branchCurrent (Il,Ib):
        for i in range(LN):
            for j in range(N[i]):
                Ib[i] += Il[idx(D2[i,j])]
        return Ib
    def nodeVoltage(Ib,V):
        for i in range(LN):
            V[idx1(D1[1,i])] =  V[idx1(D1[0,i])] - Ib[i]*Z[i]
        return V

    Il = loadCurrent(S,V,Il)
    Ib = branchCurrent(Il,Ib)
    V = nodeVoltage(Ib,V)

# Computation of losses
complexSourcePower = V[0]*np.conj(Ib[0])
realSourcePower = np.real(complexSourcePower)
reactiveSourcePower = np.imag(complexSourcePower)
totalRealPowerLoss = (realSourcePower - np.sum(np.real(S)))*MVABase*1000
totalReactivePowerLoss = (reactiveSourcePower - np.sum(np.imag(S)))*MVABase*1000

# Defining output vectors
loadCurrentVec = pd.DataFrame(Il)
branchCurrentVec = pd.DataFrame(Ib)
branchCurrentVecMag = np.abs(Ib)
nodeVoltageVec = pd.DataFrame(V)
nodeVoltageVecMag = np.abs(V)
volMag = pd.DataFrame(nodeVoltageVecMag)
branchCurrentMag = pd.DataFrame(branchCurrentVecMag)

# Printing results
print(loadCurrentVec)
print(branchCurrentVec)
print(nodeVoltageVec)
print(totalRealPowerLoss)
print(totalReactivePowerLoss)
print(volMag)
print(branchCurrentMag)

# Plotting results
slNo = np.arange(NB)
slNo1 = np.arange(LN)
plt.plot(slNo,nodeVoltageVecMag, label = "node voltage magnitude")
plt.plot(slNo1,branchCurrentVecMag, label = "branch current magnitude")
plt.title("Variation of node voltage magnitude against node number")
plt.xlabel("Node number")
plt.ylabel("Magnitude of node voltage")
plt.legend()
plt.show()
