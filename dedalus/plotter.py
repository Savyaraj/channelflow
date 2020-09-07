"""
Script for plotting flowfields and bifurcation diagrams

"""

import numpy as np
import time
import matplotlib.pyplot as plt

def plot_array():
    arr = np.load("./flowfield.npy")
    plt.plot(np.linspace(-40*np.pi,40*np.pi,len(arr)), arr, 'k-')
    plt.xlabel("$x/L_c$",fontsize=15)
    plt.ylabel("$u(x)$",fontsize=15)
    plt.title("Solution Profile, $r=-0.002$",fontsize=15)
    plt.savefig("localized_solution.png")
    plt.savefig("localized_solution.pdf")
    plt.show()

def plot_bifurcation(filenames):

    for filename in filenames:
        file = np.genfromtxt(filename,dtype=float)
        print(file.shape)
        N_arr = file[:-1,1]
        mu_arr = file[:-1,0]
        plt.plot(mu_arr, N_arr, 'k-')
        plt.xlabel("$r$",fontsize=15)
        plt.ylabel("$||u||_2$",fontsize=15)
        plt.title("Bifurcation diagram - localized solutions",fontsize=15)
    # file = np.genfromtxt('./MuE_other_side.asc',dtype=float)
    # N_arr = file[1:57,1]
    # mu_arr = file[1:57,0]
    # plt.plot(mu_arr, N_arr*0.56/0.34, 'k-')
    plt.savefig("bifurcation.png")
    plt.show()

if __name__=="__main__":
    # plot_bifurcation(['./MuE_1.asc','./MuE_2.asc'])
    plot_bifurcation(['./MuE.asc'])
    # plot_array()