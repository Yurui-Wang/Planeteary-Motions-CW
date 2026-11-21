import numpy as np 
import matplotlib.pyplot as plt

tmax = 31449600
dt = 86400
N = int(tmax/dt)
r = np.zeros((N,2))
r0 = np.array([149.6e9,0.0])
v = np.zeros((N,2))
v0 = np.array([0,29208.36397])
r[0] = r0
v[0] = v0


def a(rk):
    return (-(6.67348e-11)*(1.989e30)*(rk))/((np.sqrt(((rk[0])**2)+((rk[1])**2)))**3)

for i in range(0,N-1):
    r[i+1] = r[i] + v[i]*dt + (1/2)*a(r[i])*(dt**2)
    v[i+1] = v[i] + (1/2)*(a(r[i])+a(r[i+1]))*dt
    r[i+1] = r[i] + v[i]*dt + (1/2)*a(r[i])*(dt**2)
    v[i+1] = v[i] + (1/2)*(a(r[i])+a(r[i+1]))*dt

plt.figure()
x_values = [r[i][0] for i in range(0,N)]
y_values = [r[i][1] for i in range(0,N)]
plt.plot(x_values,y_values, "b-")
plt.scatter(0,0,color="yellow",s=500)
plt.show