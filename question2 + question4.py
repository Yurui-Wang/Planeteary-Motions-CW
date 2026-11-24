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

#Question 4
tmax = 1e8
dt = 86400
N = int(tmax/dt)

v01 = np.array([0,20000])
v02 = np.array([0,42000])
v02reverse = np.array([0,-42000])
v03 = np.array([0,32000])
v04 = np.array([0,34000])

r0 = np.array([149.6e9,0.0])

v01ver = verlet(r0, v01, dt, N)
v02r = np.concatenate([verlet(r0, v02, dt, N)[0],np.flip(verlet(r0, v02reverse, dt, N)[0], 0)])
v03ver = verlet(r0, v03, dt, N)
v04ver = verlet(r0, v04, dt, N)

plt.figure()
plt.plot(v01ver[0][:,0],v01ver[0][:,1], color="steelblue",label="v0 = 32000m/s ")
plt.plot(v02r[:,0],v02r[:,1], color="palevioletred",label="v0 = 42000m/s ")
plt.plot(v03ver[0][:,0],v03ver[0][:,1], color="plum",label="v0 = 20000m/s ")
plt.plot(v04ver[0][:,0],v04ver[0][:,1], color="lightsalmon",label="v0 = 37000m/s ")
plt.scatter(0,0,color="yellow",s=500)
y_min = -3e11
y_max = 3e11
x_min = -3e11
x_max = 3e11
plt.ylim([y_min,y_max])
plt.xlim([x_min,x_max])
plt.title("Earth's Orbit around the sun")
plt.legend()
plt.show

def peri(v0r, N):
    return min([np.sqrt(((v0r[i][0])**2)+((v0r[i][1])**2)) for i in range(N-1)])

def aph(v0r, N):
    return max([np.sqrt(((v0r[i][0])**2)+((v0r[i][1])**2)) for i in range(N-1)])

def a(v0r, N):
    return (aph(v0r, N) + peri(v0r, N))/2

def b(v0r, N):
    return np.sqrt(aph(v0r, N) * peri(v0r,N))

def e(v0r, N):
    return (aph(v0r, N) - peri(v0r, N))/(aph(v0r, N) + peri(v0r, N))

for i in [v01ver, v03ver, v04ver]:
    print("For v0 = {}:".format(i[1][0][1]))
    print("Perihelion = {}".format(peri(i[0],N)))
    print("Aphelion = {}".format(aph(i[0],N)))
    print()
    print("a = {}".format(a(i[0],N)))
    print("b = {}".format(b(i[0],N)))
    print("e = {}".format(e(i[0],N)))
    print()
