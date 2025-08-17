import numpy as np
import math

M = np.array([[1, 0], [0, 4]])              #Input Mass matrix
K = np.array([[16, 0], [0, 1]])             #Input Stiffness matrix
#---------Input damping coefficients----------#
alpha_d = 0 
beta_d = 0
C = alpha_d*M+beta_d*K
ft = 0                                      #Input initial time
f_0 = np.array([0.016*math.cos(ft), 1])     #Input Force vector
N = 6                                       #Input Number of iterations

#---------Input Newmark Beta constants----------#
gamma = 0.5 
beta = 0.25

#---------Check critical time step---------#
[om,ve]=np.linalg.eig(np.linalg.inv(M)*K)
omega = math.sqrt(max(om))
T_min = 2/omega
print("t_critical:\n", T_min)

t = 0.1                                      #Input time step

q = np.zeros((N, 2))
v = np.zeros((N, 2))
a = np.zeros((N, 2))


q[0, :] = np.array([0, 0])  
v[0, :] = np.array([0, 0])  
a[0, :] = np.linalg.inv(M).dot(f_0 - C.dot(v[0,:])-K.dot(q[0, :]))  


P1 = np.linalg.inv(M / (beta * t**2) + C*gamma / (beta* t) + K)
P2 = (1 - 2 * beta) * M / (2 * beta) + (gamma-2*beta) * t * C/ (2 * beta)
P3 = M / (beta * t) + (gamma - beta) * C / beta
P4 = M / (beta * t**2) + C*gamma / (beta* t)


for i in range(1, N):
    ft=i*t
    fo = np.array([0.016*math.cos(ft), 1])          #Input Force vector
    q[i, :] = P1.dot(fo + P2.dot(a[i-1, :]) + P3.dot(v[i-1, :]) + P4.dot(q[i-1, :]))
  
    v[i, :] = ((2 * beta - gamma) * t * a[i-1, :] / (2 * beta) + 
               (beta - gamma) * v[i-1, :] / beta - 
               gamma * q[i-1, :] / (beta * t) + 
               gamma * q[i, :] / (beta * t))
    

    a[i, :] = ((2 * beta - 1) * a[i-1, :] / (2 * beta) - 
               v[i-1, :] / (beta * t) - 
               q[i-1, :] / (beta * t**2) + 
               q[i, :] / (beta * t**2))


print("Displacements:\n", q)
print("Velocities:\n", v)
print("Accelerations:\n", a)
