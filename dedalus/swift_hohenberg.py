"""
1D Swift-Hohenberg equation 

borrowed from Burke and Knobloch (2006)

"""

import numpy as np
import time
import matplotlib.pyplot as plt
import math 
from mpi4py import MPI
from dedalus import public as de
from dedalus.extras.plot_tools import quad_mesh, pad_limits

import itertools
import logging
logger = logging.getLogger(__name__)

class DedalusPy:
    
    def __init__(self):
        self.domain_setup()
        self.problem_setup(-0.002)
        self.build_solver()
        self.init_problem()
        print("-----------dedalus problem setup complete--------------\n")

    def domain_setup(self):
        # Bases and domain
        self.Nx = 100
        self.Lc = 2*np.pi/0.5   # (2*pi/qc)
        self.Lx = 20*self.Lc
        self.Ny = 4
        self.Nz = 4
        self.Nd = 3
        self.scale = 3/2
        self.x_basis = de.Fourier('x', self.Nx, interval=(0, self.Lx), dealias=self.scale)
        self.domain = de.Domain([self.x_basis], grid_dtype=np.float64)

    def problem_setup(self, mu):
        # Problem
        self.problem = de.IVP(self.domain, variables=['u', 'ux', 'uxx'])
        self.problem.parameters['r'] = mu
        self.problem.parameters['qc'] = 0.5
        self.problem.parameters['v'] = 0.41
        self.problem.parameters['g'] = 1
        # self.problem.parameters['u0'] = np.ones(self.Nx)
        # self.problem.add_equation("dt(u) - a*dx(ux) - b*dx(uxx) = -u*ux")
        self.problem.add_equation("dt(u) -(r-qc*qc*qc*qc)*u + 2*qc*qc*uxx + dx(dx(uxx)) = v*u*u - g*u*u*u")
        self.problem.add_equation("ux - dx(u) = 0")
        self.problem.add_equation("uxx - dx(ux) = 0")

    def build_solver(self):
        # Build solver
        self.solver = self.problem.build_solver(de.timesteppers.RK443)
        self.solver.stop_wall_time = np.inf
        self.solver.stop_iteration = 5000

    def init_problem(self):
        # Initial conditions
        self.x = self.domain.grid(0)
        self.u = self.solver.state['u']
        self.ux = self.solver.state['ux']
        self.uxx = self.solver.state['uxx']

        self.n = 20
        r = self.problem.parameters['r']
        qc = self.problem.parameters['qc']
        gamma3 = 8.356
        # self.u['g'] = 0.1*np.ones(len(self.u['g']))
        # self.u['g'] = 1*np.cos(2*np.pi*(self.x-self.Lx/2)/self.Lc) #0.1*np.random.rand(len(self.u['g'])) #.random.rand(self.Nx) #0.1*np.sin(2*self.x)  #
        self.u['g'] = 2*math.sqrt(-2*r/gamma3)*(1/np.cosh((self.x-self.Lx/2)*math.sqrt(-r)/(2*qc)))*(0.15 + np.cos(qc*(self.x-self.Lx/2)))
        self.u.differentiate(0, out=self.ux)
        self.ux.differentiate(0, out=self.uxx)

        # Store data for final plot
        self.u.set_scales(1)
        self.u_list = [np.copy(self.u['g'])]
        self.ux_list = [np.copy(self.ux['g'])]
        self.t_list = [self.solver.sim_time]

    def advance(self, T, u_init = None):
        # Main loop
        dt = 2e-2 #self.CFL.compute_dt()
        init_time = self.solver.sim_time
        self.u.set_scales(1)
        if u_init is not None:
            if len(u_init) == self.Nx:
                print("Reading input array in Python of length: "+str(len(u_init)))
                self.u['g'][:] = u_init[:]
            # if len(u_init) == self.Nx*self.Ny*self.Nz*self.scale*self.Nd: 
            #     print("Reading input array in Python of length: "+str(len(u_init)))
            #     for i, nx, ny, nz in itertools.product(range(self.Nd),range(self.Nx),range(self.Ny),range(self.Nz)):
            #         if i == 0: self.u['g'][nx] = u_init[i*(self.Ny+self.Nx+self.Nz)+ny*(self.Nx+self.Nz)+nx*self.Nz+nz]
                    # elif i == 1: self.v[nx,ny,nz] = u_init[i*(self.Ny+self.Nx+self.Nz)+ny*(self.Nx+self.Nz)+nx*self.Nz+nz]
                    # elif i == 2: self.w[nx,ny,nz] = u_init[i*(self.Ny+self.Nx+self.Nz)+ny*(self.Nx+self.Nz)+nx*self.Nz+nz]
        # T = 100*dt
        # print("2-Norm before timestepping is: "+str(np.linalg.norm(self.u['g'],2)))
        # print("u before timestepping is: "+str(self.u['g'][0]))
        print("Integrating in dedalus, T = "+str(T))
        u_temp = np.copy(self.u['g'])
        self.u.set_scales(self.scale)
        while init_time + T > self.solver.sim_time:
            self.solver.step(dt)
        self.u.set_scales(1)
        print("G: "+str(np.linalg.norm(self.u['g']-u_temp,1)/T))
        # print("u after timestepping is: "+str(self.u['g'][0]))
        # print("Length of output array in python is:"+str(len(self.u['g']))+"\n)
        print("2-Norm of u:"+str(np.linalg.norm(self.u['g'],2))+'\n')
        return(self.u['g'])

    def updateMu(self, mu):
        self.domain_setup()
        self.problem_setup(mu)
        self.build_solver()
        self.init_problem()
        print("updating mu to: "+str(mu))
        return 

    def save_array(self, arr):
        np.save("./flowfield", arr)
        return 

if __name__=="__main__":
    sh = DedalusPy() 
    try:
        logger.info('Starting loop')
        start_time = time.time()
        # sh.problem.parameters['r'] = 0.1
        while sh.solver.proceed:
            dt = 2e-2
            sh.solver.step(dt)
            if sh.solver.iteration % 20 == 0:
                sh.u.set_scales(1)
                sh.u_list.append(np.copy(sh.u['g']))
                sh.ux.set_scales(1)
                sh.ux_list.append(np.copy(sh.ux['g']))
                sh.t_list.append(sh.solver.sim_time)
            if sh.solver.iteration % 1000 == 0:
                logger.info('Iteration: %i, Time: %e, dt: %e' %(sh.solver.iteration, sh.solver.sim_time, dt))
    except:
        logger.error('Exception raised, triggering end of main loop.')
        raise
    finally:
        end_time = time.time()
        logger.info('Iterations: %i' %sh.solver.iteration)
        logger.info('Sim end time: %f' %sh.solver.sim_time)
        logger.info('Run time: %.2f sec' %(end_time-start_time))
        logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*sh.domain.dist.comm_cart.size))

    # Create space-time plot
    u_array = np.array(sh.u_list)
    ux_array = np.array(sh.ux_list)
    print(str(math.sqrt(2)*np.linalg.norm(u_array[-1],2))+'\n')
    print(u_array[-1])
    print(np.linalg.norm(u_array[-1]-u_array[0],1)/1000)
    t_array = np.array(sh.t_list)
    xmesh, ymesh = quad_mesh(x=sh.x, y=t_array)
    plt.figure()
    plt.pcolormesh(xmesh, ymesh, u_array, cmap='RdBu_r')
    plt.axis(pad_limits(xmesh, ymesh))
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('t')
    plt.title('Swift-Honenberg, r=%g' %(sh.problem.parameters['r']))
    # plt.savefig('Swift_Honenberg.pdf',format='pdf')
    plt.savefig('Swift_Honenberg.png')
    plt.show()
    # sh.save_array(u_array[-1], "./")

