import numpy as np 
import matplotlib.pyplot as plt

tmax = 31449600
dt = 86400
N = int(tmax/dt)
r0 = np.array([149.6e9,0.0])
v0 = np.array([0,29208.36397])
G = 6.67348e-11
M = 1.989e30

def acceleration(r):
    return (-G*M*(r))/((np.sqrt(((r[0])**2)+((r[1])**2)))**3)

def verlet(r0, v0, dt, N):
    v = np.zeros((N,2))
    r = np.zeros((N,2))   
    r[0] = r0
    v[0] = v0
    
    for i in range(0,N-1):
        r[i+1] = r[i] + v[i]*dt + (1/2)*acceleration(r[i])*(dt**2)
        v[i+1] = v[i] + (1/2)*(acceleration(r[i])+acceleration(r[i+1]))*dt
        
    return r, v 

position, velocity = verlet(r0, v0, dt, N)

plt.figure(figsize=(5,5))
x_values = position[:, 0]
y_values = position[:, 1]
plt.plot(x_values,y_values, "b-")
plt.scatter(0,0,color="yellow",s=500)
plt.title("Earth's Orbit around the sun")
plt.show