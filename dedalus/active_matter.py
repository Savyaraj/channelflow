"""
Continuum model of bacterial turbulence 

borrowed from 'Meso-scale turbulence in living fluids, Wensink et al. (2012)'

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
        self.dt = 3e-2
        self.N = 128
        self.T0 = -0.045
        self.T2 = 9.1125e-5
        self.a = self.T0**2/(4*self.T2) - 0.002 
        self.k_max = np.sqrt(-self.T0/(2*self.T2))
        self.timestepper = de.timesteppers.RK443
        self.Niter = 10000
        self.domain_setup()
        self.problem_setup(self.a)
        self.build_solver(self.Niter)
        self.init_problem()
        print("-----------dedalus problem setup complete--------------\n")

    def domain_setup(self):
        # Bases and domain
        self.Nx = self.N
        self.Ny = self.N 
        self.Lc = 2*np.pi/self.k_max #np.pi 
        self.Lx = 10*self.Lc 
        self.Ly = 10*self.Lc 
        self.Nz = 6
        self.Nd = 3
        self.scale = 3/2
        self.x_basis = de.Fourier('x', self.Nx, interval=(0, self.Lx), dealias=self.scale)
        self.y_basis = de.Fourier('y', self.Ny, interval=(0, self.Ly), dealias=self.scale)
        self.domain = de.Domain([self.x_basis, self.y_basis], grid_dtype=np.float64)

    def problem_setup(self, mu = -1):
        # Problem

        #variables: flowfield (u,v)(x,y), pressure p(x,y), vorticity w(x,y) 

        self.problem = de.IVP(self.domain, variables=['u','v','w','p']) 
        self.problem.parameters['a'] = mu
        self.problem.parameters['b'] = 0.5
        self.problem.parameters['Tau_0'] = self.T0
        self.problem.parameters['Tau_2'] = self.T2
        self.problem.parameters['S'] = -2.5   

        #Incompressibility 
        self.problem.add_equation("dx(u) + dy(v) = 0", condition="(nx!=0) or (ny!=0)")
        self.problem.add_equation("p = 0",condition="(nx==0) and (ny==0)")    

        #dynamics (-S)*dx(u+v) 
        self.problem.add_equation("dt(u) + dx(p) + a*u - Tau_0*(d(u,x=2) + d(u,y=2)) + Tau_2*(d(u,x=4) + d(u,y=4)) = \
                                   -(1-S)*(u*dx(u) + v*dy(u)) + (-S/2)*dx(u*u) - b*(u*u+v*v)*u")
        self.problem.add_equation("dt(v) + dy(p) + a*v - Tau_0*(d(v,x=2) + d(v,y=2)) + Tau_2*(d(v,x=4) + d(v,y=4)) = \
                                   -(1-S)*(u*dx(v) + v*dy(v)) + (-S/2)*dy(v*v) - b*(u*u+v*v)*v") 
        # self.problem.add_equation("dt(u) + dx(p) + a*u - Tau_0*d(u,x=2) + Tau_2*d(u,x=4) = 0")
        # self.problem.add_equation("dt(v) + dy(p) + a*v - Tau_0*d(v,y=2) + Tau_2*d(v,y=4) = 0")

        #vorticity
        self.problem.add_equation("w - dx(v) + dy(u) = 0")

    def build_solver(self, n_iter = 1000):
        # Build solver
        self.solver = self.problem.build_solver(self.timestepper)
        self.solver.stop_wall_time = np.inf
        self.solver.stop_iteration = n_iter 
        # self.solver.stop_sim_time = 180

    def init_problem(self):
        # Initial conditions
        self.x = self.domain.grid(0)
        self.y = self.domain.grid(1)
        self.u = self.solver.state['u']
        self.v = self.solver.state['v']
        self.w = self.solver.state['w']
        self.p = self.solver.state['p']
        self.u.set_scales(1)
        self.v.set_scales(1)
        #set initial flowfield at random perturbation from queiscent state
        # self.u['g'] = 0.1*np.random.rand(np.shape(self.p['g'])[0],np.shape(self.p['g'])[1])
        # k_max = 15
        self.u['g'] = 1e-2*np.cos(self.k_max*self.y)*np.ones(self.u['g'].shape)
        # self.v['g'] = -1*np.sin(k_max*self.x)*self.y*np.ones(self.v['g'].shape)
        # self.v['g'] = 0.1*np.random.rand(np.shape(self.p['g'])[0],np.shape(self.p['g'])[1])
        
        # Store data for final plot
        self.u_list = [np.copy(self.u['g'])]
        self.v_list = [np.copy(self.v['g'])]
        self.w_list = [np.copy(self.w['g'])]
        self.t_list = [self.solver.sim_time]

    def advance(self, T, u_init = None):
        # Main loop
        init_time = self.solver.sim_time
        self.u.set_scales(1)
        self.v.set_scales(1)
        u_chflow = np.zeros([self.Nx,self.Nz,self.Ny,self.Nd])

        if u_init is not None:
            if len(u_init) == self.Nx*self.Ny*self.Nz*self.Nd: 
                print("Reading input array in Python of length: "+str(len(u_init)))
                ind = 0
                for i in range(self.Nd):
                    for ny in range(self.Nz):
                        for nx in range(self.Nx):
                            for nz in range(self.Ny):
                                u_chflow[nx,ny,nz,i] = u_init[ind]
                                ind = ind + 1

        self.u['g'] = u_chflow[:,2,:,0]
        self.v['g'] = u_chflow[:,2,:,2]

        print("Integrating in dedalus, T = "+str(T))

        print("u before:"+str(np.linalg.norm(self.u['g'][0]))+'\n')#str(np.linalg.norm(self.u['g'],2))+'\n')
        
        while init_time + T > self.solver.sim_time:
            self.solver.step(self.dt)
        self.u.set_scales(1)

        print("u after:"+str(np.linalg.norm(self.u['g'][0]))+'\n')#str(np.linalg.norm(self.u['g'],2))+'\n')
        # print("G: "+str(np.linalg.norm(self.u['g']-u_temp,1)/T))

        u_out = np.zeros(len(u_init))
        ind = 0
        for i in range(self.Nd):
            for ny in range(self.Nz):
                for nx in range(self.Nx):
                    for nz in range(self.Ny):
                        if ny == 2 and i == 0: u_out[ind] = self.u['g'][nx,nz]
                        elif ny == 2 and i == 2: u_out[ind] = self.v['g'][nx,nz]
                        ind = ind + 1
        print("G: "+str(np.linalg.norm(u_out-u_init,2)/T))
        return(u_out)


    # def advance(self, T, u_init = None):
    #     self.rank = MPI.COMM_WORLD.Get_rank()
    #     self.size = int(self.Nx/MPI.COMM_WORLD.Get_size())
    #     print("This is process no. "+str(MPI.COMM_WORLD.Get_rank()))
    #     # Main loop
    #     init_time = self.solver.sim_time
    #     self.u.set_scales(1)
    #     self.v.set_scales(1)
    #     u_chflow = np.zeros([self.Nd,self.Nz,self.Nx,self.Ny])
        
    #     if u_init is not None:
    #         if len(u_init) == self.Nx*self.Ny*self.Nz*self.Nd: 
    #             print("Reading input array in Python of length: "+str(len(u_init)))
    #             ind = 0
    #             for i in range(self.Nd):
    #                 for ny in range(self.Nz):
    #                     for nx in range(self.Nx):
    #                         for nz in range(self.Ny):
    #                             u_chflow[i,ny,nx,nz] = u_init[ind]
    #                             ind = ind + 1
        
    #     print(np.shape(u_chflow[0,2,:,self.rank*self.size:(self.rank+1)*self.size]))
        
    #     self.u['g'] = u_chflow[0,2,:,self.rank*self.size:(self.rank+1)*self.size]
    #     self.v['g'] = u_chflow[2,2,:,self.rank*self.size:(self.rank+1)*self.size]

    #     print("Integrating in dedalus, T = "+str(T))

    #     u_temp = np.copy(self.u['g'])
    #     print("2-Norm of u before:"+str(np.linalg.norm(self.u['g'],2))+'\n')
        
    #     while init_time + T > self.solver.sim_time:
    #         self.solver.step(self.dt)
    #     self.u.set_scales(1)

    #     print("2-Norm of u after:"+str(np.linalg.norm(self.u['g'],2))+'\n')
    #     print("u after:"+str(self.u['g'][0])+'\n')
    #     # print("G: "+str(np.linalg.norm(self.u['g']-u_temp,1)/T))

    #     u_out = np.zeros(len(u_init))
    #     ind = 0
    #     for i in range(self.Nd):
    #         for ny in range(int(self.Nz/self.size)):
    #             for nx in range(self.Nx):
    #                 for nz in range(self.Ny):
    #                     if ny == 2 and i == 0: u_out[ind] = self.u['g'][nx,nz+self.rank*self.size]
    #                     elif ny == 2 and i == 2: u_out[ind] = self.v['g'][nx,nz+self.rank*self.size]
    #                     ind = ind + 1
    #     print("G: "+str(np.linalg.norm(u_out-u_init,2)/T))
    #     return(u_out)

    def updateMu(self, mu):
        self.domain_setup()
        self.problem_setup(mu)
        self.build_solver()
        self.init_problem()
        print("updating mu to: "+str(mu))
        return 

    def save_array(self, arr):
        u_chflow = np.zeros([self.Nx,self.Nz,self.Ny,self.Nd])

        ind = 0
        for i in range(self.Nd):
            for ny in range(self.Nz):
                for nx in range(self.Nx):
                    for nz in range(self.Ny):
                        u_chflow[nx,ny,nz,i] = arr[ind]
                        ind = ind + 1
        np.save("./flowfield", u_chflow)
        return 

# def main():
if __name__=="__main__":
    # arr = np.ones(10)
    # print(arr)
    # return np.zeros(len(arr))
    
    Dd = DedalusPy()
    try:
        logger.info('Starting loop')
        start_time = time.time()
        # sh.problem.parameters['r'] = 0.1
        while Dd.solver.proceed:
            Dd.solver.step(Dd.dt)
            if Dd.solver.iteration % 20 == 0:
                Dd.u.set_scales(1)
                Dd.u_list.append(np.copy(Dd.u['g']))
                Dd.v.set_scales(1)
                Dd.v_list.append(np.copy(Dd.v['g']))
                Dd.w.set_scales(1)
                Dd.w_list.append(np.copy(Dd.w['g']))
                Dd.t_list.append(Dd.solver.sim_time)
            if Dd.solver.iteration % 100 == 0:
                logger.info('Iteration: %i, Time: %e, dt: %e, 2-norm: %e' %(Dd.solver.iteration, Dd.solver.sim_time, Dd.dt,np.linalg.norm(Dd.u['g'][0],2)))
                # plt.imshow(Dd.u['g'], cmap='hot', interpolation='nearest') # see previous comments
                # plt.pause(0.001)
    except:
        logger.error('Exception raised, triggering end of main loop.')
        raise
    finally:
        end_time = time.time()
        logger.info('Iterations: %i' %Dd.solver.iteration)
        logger.info('Sim end time: %f' %Dd.solver.sim_time)
        logger.info('Run time: %.2f sec' %(end_time-start_time))
        logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*Dd.domain.dist.comm_cart.size))
        

    # Create space-time plot
    u_array = np.array(Dd.u_list)
    v_array = np.array(Dd.v_list) 
    w_array = np.array(Dd.w_list) 
    t_array = np.array(Dd.t_list)

    plt.figure()
    # xmesh, ymesh = quad_mesh(x=Dd.x[:,0], y=Dd.y[:,0])
    # xmesh, ymesh = np.meshgrid(Dd.x,Dd.y)
    # plt.pcolormesh(xmesh, ymesh, w_array[-1], cmap='RdBu_r')
    # plt.axis(pad_limits(xmesh, ymesh))
    # plt.colorbar()
    plt.imshow(u_array[-1])
    plt.xlabel('x')
    plt.ylabel('y')

    #  -------------edited for travelling wave ----------------
    plt.title('Flowfield , alpha=%g' %(Dd.problem.parameters['a']))   
    # plt.savefig('Vorticity.png')
    plt.show()

# if __name__=="__main__":
#     main()