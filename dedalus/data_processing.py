import numpy as np 
import csv 
import matplotlib.pyplot as plt 
from scipy import interpolate

with open("init.csv") as f:
    file = csv.reader(f, delimiter=',')
    x = []
    u = []
    for row in file:
        x.append(float(row[0]))
        u.append(float(row[1]))

inds = np.array(x).argsort()
x = sorted(x)
u = np.array(u)[inds]
plt.plot(x,u)
tck = interpolate.splrep(x,u)
xnew = np.linspace(-40*np.pi,40*np.pi,500)
# unew = interpolate.splev(xnew, tck)
unew = np.interp(xnew, x, u)
# print(unew)
np.savetxt("u_init.txt",unew, delimiter=',')
plt.plot(xnew,unew)
plt.show()