# Finding optimal value of shunt capacitor and optimal position of shunt capacitor placement
# Importing modules
from matplotlib import pyplot as plt
from matplotlib import style
style.use('ggplot')
import numpy as np
import pandas as pd

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

# Defining ratings of shunt capacitor
QLsh = np.array([-0.0005,-0.0006,-0.0007,-0.0008,-0.0009,-0.001,-0.0011,-0.0012,-0.0013
,-0.0014,-0.0015,-0.0016,-0.0017,-0.0018,-0.0019,-0.002,-0.0021,-0.0022,-0.0023,-0.0024
,-0.0025,-0.0026,-0.0027,-0.0028,-0.0029,-0.003])

# Defining parameters
slNO = np.array([1,2,3,4,5,6,7,8,9])
NB = 9
LN = 8
maxIter = 4
tol = 0.000000001
w = len(QLsh)

# Initialsing node voltages and losses
V = np.ones(NB)
Vold = np.ones(NB)
Ploss = np.zeros(LN)
Qloss = np.zeros(LN)
P = np.zeros(NB)
Q = np.zeros(NB)

# Initialsing real power loss vector
rPLoss = np.zeros(w)
newVar = np.zeros(w)

# Initialsing functional parameters
g = np.zeros(LN)
h = np.zeros(LN)
A = np.zeros(LN)
D = np.zeros(LN)

# Calculation of optimal value of shunt capacitor and optimal position of it's placement
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

for k in range(w):

    QL[2] = QL[2] + QLsh[k]
        
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
    
    totalRealPowerLoss = np.sum(Ploss)
    rPLoss[k] = totalRealPowerLoss

# Printing results  
rPLossVec = pd.DataFrame(rPLoss)
print(rPLossVec)
print(np.argmin(rPLoss))
print(np.min(rPLoss))

# Printing results
slNo = np.arange(w)
plt.plot(slNo,rPLoss)
plt.title("Variation of total real power loss against ratings of shunt capacitor")
plt.xlabel("Rating of shunt capacitor in p.u KVAr")
plt.ylabel("Total real power loss in p.u MW")
plt.show()
