import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import style
style.use('ggplot')

# Analysis of balanced radial distribution network
# Defining input matrices
D1 = np.array([[1,2,3,4,5,2,7,4],[2,3,4,5,6,7,8,9]])
D2 = np.array([[2,3,4,5,6,7,8,9],[3,4,5,6,9,0,0,0],[4,5,6,9,0,0,0,0],[5,6,0,0,0,0,0,0],[6,0,0,0,0,0,0,0],[7,8,0,0,0,0,0,0],[8,0,0,0,0,0,0,0],[9,0,0,0,0,0,0,0]])
E1 = np.array([[2,3,4,5,6,7,8],[3,4,5,8,0,0,0],[4,5,8,0,0,0,0],[5,0,0,0,0,0,0],[0,0,0,0,0,0,0],[7,0,0,0,0,0,0],[0,0,0,0,0,0,0],[0,0,0,0,0,0,0]])

# Defining input vectors
N = np.array([8,5,4,2,1,2,1,1])
B = np.array([7,4,3,1,0,1,0,0])
S = np.array([0.004+0.003j,0.005+0.004j,0.002+0.004j,0.004+0.000j,0.004+0.000j,0.001+0.000j,0.001+0.002j,0.001+0.000j])
Z = np.array([0.2 + 0.3721j, 0.3388+ 0.1j, 0+0.75196j, 0.1+0.4j,0.4+0.3j,0.5+0.6j,0.35+0.4j,0.4+0.2j])

# Defining real and imaginary parts
r = np.real(Z)
x = np.imag(Z)
PL = np.real(S)
QL = np.imag(S)
QL[0] = QL[0] -0.0017

# Defining parameters
slNO = np.array([1,2,3,4,5,6,7,8,9])
NB = 9
LN = 8
maxIter = 3
tol = 0.000000001

# Initialsing node voltages and losses
V = np.ones(NB)
Vold = np.ones(NB)
Ploss = np.zeros(LN)
Qloss = np.zeros(LN)
P = np.zeros(NB)
Q = np.zeros(NB)

# Initialsing functional parameters
g = np.zeros(LN)
h = np.zeros(LN)
A = np.zeros(LN)
D = np.zeros(LN)

# Calculation of node voltages and branch losses
def idx(x):
    return x-2

def idx1(x):
    return x-1

def totalPower(m,n,o):
    for i in range(LN):
        for j in range(N[i]):
            g[i] += o[idx(D2[i,j])]
        for j in range(B[i]):
            h[i] += n[idx(E1[i,j])]
        m[i+1] = g[i] + h[i]
    return m

P = totalPower(P,Ploss,PL)
g = np.zeros(LN)
h = np.zeros(LN)
Q = totalPower(Q,Qloss,QL)
g = np.zeros(LN)
h = np.zeros(LN)

for j in range(maxIter):
    for i in range(LN):
    
        A[i] = P[idx1(D1[1,i])]*r[i] + Q[idx1(D1[1,i])]*x[i] - 0.5*V[idx1(D1[0,i])]**2
        D[i] = np.sqrt(A[i]**2 - ((r[i]**2 + x[i]**2)*(P[idx1(D1[1,i])]**2 + Q[idx1(D1[1,i])]**2)))
        V[idx1(D1[1,i])] = np.sqrt(D[i] - A[i])
        Ploss[i] = (r[i]*(P[idx1(D1[1,i])]**2 + Q[idx1(D1[1,i])]**2))/V[idx1(D1[1,i])]**2
        Qloss[i] = (x[i]*(P[idx1(D1[1,i])]**2 + Q[idx1(D1[1,i])]**2))/V[idx1(D1[1,i])]**2

    P = totalPower(P,Ploss,PL)
    g = np.zeros(LN)
    h = np.zeros(LN)
    Q = totalPower(Q,Qloss,QL)
    g = np.zeros(LN)
    h = np.zeros(LN)
    
# Defining results
totalRealPowerLoss = np.sum(Ploss)
totalReactivePowerLoss = np.sum(Qloss)
volMag = pd.DataFrame(V)
pLoss = pd.DataFrame(Ploss)
qLoss = pd.DataFrame(Qloss)
p = pd.DataFrame(P)
q = pd.DataFrame(Q)

# Analysis of voltage stability 
# Initialising b,c and SI vectors
b = np.zeros(LN)
c = np.zeros(LN)
SI = np.zeros(NB)
SImin = 0
crtclnode = 0

# Calculation of voltage stability index
SI[0] = 50 # just for the sake of calculation
for i in range(LN):
    b[i] = V[idx1(D1[0,i])]**2 -2*P[idx1(D1[1,i])]*r[i] -2*Q[idx1(D1[1,i])]*x[i]
    c[i] = (r[i]**2 + x[i]**2)*(P[idx1(D1[1,i])]**2 + Q[idx1(D1[1,i])]**2)
    SI[i+1] = b[i]**2 -4*c[i]
    SImin = np.min(SI)

stbleIndex = pd.DataFrame(SI[1:9])
crtclnode = np.argmin(SI)

# Printing results
print(volMag)
print(pLoss)
print(qLoss)
print(p)
print(q)
print(totalRealPowerLoss)
print(totalReactivePowerLoss)
print(stbleIndex)
print(SImin)
print(crtclnode)

# Plotting the variation of voltage stability index and node voltage
plt.plot(slNO,V,label = 'node voltage')
plt.plot(slNO[1:9],SI[1:9], label ='voltage stability index')
plt.title("Variation of voltage magnitude and voltage stability index")
plt.ylabel("Voltage magnitude and voltage stability index")
plt.xlabel("slNo")
plt.legend()
plt.show()
