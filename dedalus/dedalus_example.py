"""
1D Korteweg-de Vries / Burgers equation

This script should be ran serially (because it is 1D), and creates a space-time
plot of the computed solution.

"""

import numpy as np
import time
import matplotlib.pyplot as plt

from dedalus import public as de
from dedalus.extras.plot_tools import quad_mesh, pad_limits

import logging
logger = logging.getLogger(__name__)

class DedalusPy:
    
    def __init__(self):
        self.domain_setup()
        self.problem_setup()
        self.build_solver()
        self.init_problem()
        print("-----------dedalus problem setup complete--------------\n")

        self.sum = 0

    def domain_setup(self):
        # Bases and domain
        self.x_basis = de.Fourier('x', 1024, interval=(0, 8), dealias=3/2)
        self.domain = de.Domain([self.x_basis], np.float64)

    def problem_setup(self):
        # Problem
        self.problem = de.IVP(self.domain, variables=['u', 'ux', 'uxx'])
        self.problem.parameters['a'] = 2e-4
        self.problem.parameters['b'] = 1e-4
        self.problem.add_equation("dt(u) - a*dx(ux) - b*dx(uxx) = -u*ux")
        self.problem.add_equation("ux - dx(u) = 0")
        self.problem.add_equation("uxx - dx(ux) = 0")

    def build_solver(self):
        # Build solver
        self.solver = self.problem.build_solver(de.timesteppers.SBDF2)
        self.solver.stop_wall_time = 60
        self.solver.stop_iteration = 5000

    def init_problem(self):
        # Initial conditions
        self.x = self.domain.grid(0)
        self.u = self.solver.state['u']
        self.ux = self.solver.state['ux']
        self.uxx = self.solver.state['uxx']

        self.n = 20
        self.u['g'] = np.log(1 + np.cosh(self.n)**2/np.cosh(self.n*self.x)**2) / (2*self.n)
        self.u.differentiate(0, out=self.ux)
        self.ux.differentiate(0, out=self.uxx)

        # Store data for final plot
        self.u.set_scales(1)
        self.u_list = [np.copy(self.u['g'])]
        self.t_list = [self.solver.sim_time]

    def advance(self, u_init = None):
        # Main loop
        dt = 2e-3

        print("Reading input array in Python \n \n Array length: "+str(len(u_init))+"\n")
        for i in range(min(len(u_init),len(self.u['g']))):
                self.u['g'][i] = u_init[i]

        self.solver.step(dt)
        print("Length of output array in python is:"+str(len(self.u['g']))+"\n 2-Norm of the array is:"+str(np.linalg.norm(self.u['g'],2))+"\n")
        return(self.u['g'])

    def add(self,i):
        self.sum = self.sum + i

    def printSum(self):
        print(self.sum)
    

    

    

    

    

# # try:
# #     logger.info('Starting loop')
# #     start_time = time.time()
# #     while solver.proceed:
# #         solver.step(dt)
# #         if solver.iteration % 20 == 0:
# #             u.set_scales(1)
# #             u_list.append(np.copy(u['g']))
# #             t_list.append(solver.sim_time)
# #         if solver.iteration % 100 == 0:
# #             logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
# # except:
# #     logger.error('Exception raised, triggering end of main loop.')
# #     raise
# # finally:
# #     end_time = time.time()
# #     logger.info('Iterations: %i' %solver.iteration)
# #     logger.info('Sim end time: %f' %solver.sim_time)
# #     logger.info('Run time: %.2f sec' %(end_time-start_time))
# #     logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))

# # # Create space-time plot
# # u_array = np.array(u_list)
# # t_array = np.array(t_list)
# # xmesh, ymesh = quad_mesh(x=x, y=t_array)
# # plt.figure()
# # plt.pcolormesh(xmesh, ymesh, u_array, cmap='RdBu_r')
# # plt.axis(pad_limits(xmesh, ymesh))
# # plt.colorbar()
# # plt.xlabel('x')
# # plt.ylabel('t')
# # plt.title('KdV-Burgers, (a,b)=(%g,%g)' %(problem.parameters['a'], problem.parameters['b']))
# # plt.savefig('kdv_burgers.png')

