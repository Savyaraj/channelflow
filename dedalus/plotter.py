"""
Script for plotting flowfields and bifurcation diagrams

"""

import numpy as np
import time
import matplotlib.pyplot as plt

def plot_array():
    arr = np.load("./flowfield.npy")
    plt.plot(np.linspace(0,2*2*np.pi/1.25,len(arr)), arr, 'k-')
    plt.xlabel("$x/L_c$",fontsize=15)
    plt.ylabel("$u(x)$",fontsize=15)
    plt.title("Solution Profile, $r=-0.01$",fontsize=15)
    plt.savefig("periodic_soln_0.01.png")
    plt.show()

def plot_bifurcation(filename):
    file = np.genfromtxt(filename,dtype=float)
    print(file.shape)
    N_arr = file[:,1]
    mu_arr = file[:,0]
    plt.plot(mu_arr, N_arr*0.25/0.147, 'k-')
    plt.xlabel("$r$",fontsize=15)
    plt.ylabel("$||u||_2$",fontsize=15)
    plt.title("Bifurcation diagram - periodic branch",fontsize=15)
    file = np.genfromtxt('./MuE_other_side.asc',dtype=float)
    N_arr = file[1:57,1]
    mu_arr = file[1:57,0]
    plt.plot(mu_arr, N_arr*0.25/0.147, 'k-')
    plt.savefig("bifurcation.png")
    plt.show()

if __name__=="__main__":
    plot_bifurcation('./MuE.asc')
    # plot_array()