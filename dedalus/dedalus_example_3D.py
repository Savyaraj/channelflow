"""
Dedalus script for 3D Rayleigh-Benard convection.

"""

import numpy as np
from mpi4py import MPI
import time
import itertools

from dedalus import public as de
from dedalus.extras import flow_tools

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
        # Parameters
        Lx, Ly, Lz = (25., 25., 1.)
        self.Nd = 3
        self.Nx = 8
        self.Ny = 4
        self.Nz = 8
        self.scale = 1
        epsilon = 0.8
        self.Pr = 1.0

        Ra_crit = 1707.762  # No-slip top & bottom
        self.Ra = Ra_crit * (1 + epsilon)
        self.F = self.Ra*self.Pr
        # Create bases and self.domain
        self.x_basis = de.Fourier('x', self.Nx, interval=(0, Lx), dealias=self.scale)
        self.y_basis = de.Fourier('y', self.Ny, interval=(0, Ly), dealias=self.scale)
        self.z_basis = de.Chebyshev('z', self.Nz, interval=(-Lz/2, Lz/2), dealias=self.scale)
        self.domain = de.Domain([self.x_basis, self.y_basis, self.z_basis], grid_dtype=np.float64)

    def problem_setup(self):
        # 2D Boussinesq hydrodynamics
        self.problem = de.IVP(self.domain, variables=['p','b','u','v','w','bz','uz','vz','wz'], time='t')
        self.problem.meta['p','b','w','bz','wz']['x','y']['parity'] = 1
        self.problem.meta['v','vz']['x']['parity'] = 1
        self.problem.meta['v','vz']['y']['parity'] = -1
        self.problem.meta['u','uz']['x']['parity'] = -1
        self.problem.meta['u','uz']['y']['parity'] = 1
        self.problem.parameters['P'] = 1
        self.problem.parameters['R'] = self.Pr
        self.problem.parameters['F'] = self.F

        self.problem.add_equation("dx(u) + dy(v) + wz = 0")
        self.problem.add_equation("dt(b) - P*(dx(dx(b)) + dy(dy(b)) + dz(bz))             = - u*dx(b) - v*dy(b) - w*bz")
        self.problem.add_equation("dt(u) - R*(dx(dx(u)) + dy(dy(u)) + dz(uz)) + dx(p)     = - u*dx(u) - v*dy(u) - w*uz")
        self.problem.add_equation("dt(v) - R*(dx(dx(v)) + dy(dy(v)) + dz(vz)) + dy(p)     = - u*dx(v) - v*dy(v) - w*vz")
        self.problem.add_equation("dt(w) - R*(dx(dx(w)) + dy(dy(w)) + dz(wz)) + dz(p) - b = - u*dx(w) - v*dy(w) - w*wz")
        self.problem.add_equation("bz - dz(b) = 0")
        self.problem.add_equation("uz - dz(u) = 0")
        self.problem.add_equation("vz - dz(v) = 0")
        self.problem.add_equation("wz - dz(w) = 0")

        self.problem.add_bc("left(b) = -left(F*z)")
        self.problem.add_bc("left(u) = 0")
        self.problem.add_bc("left(v) = 0")
        self.problem.add_bc("left(w) = 0")
        self.problem.add_bc("right(b) = -right(F*z)")
        self.problem.add_bc("right(u) = 0")
        self.problem.add_bc("right(v) = 0")
        self.problem.add_bc("right(w) = 0", condition="(nx != 0) or (ny != 0)")
        self.problem.add_bc("integ_z(p) = 0", condition="(nx == 0) and (ny == 0)")

    def build_solver(self):
        # Build solver
        self.solver = self.problem.build_solver(de.timesteppers.MCNAB2)
        logger.info('Solver built')
        # Integration parameters
        self.solver.stop_sim_time = 100
        self.solver.stop_wall_time = 60 * 60.
        self.solver.stop_iteration = np.inf

    def init_problem(self):
        # Initial conditions
        self.u = self.solver.state['u']['g']
        self.v = self.solver.state['v']['g']
        self.w = self.solver.state['w']['g']
        z = self.domain.grid(2)
        b = self.solver.state['b']
        bz = self.solver.state['bz']

        # Random perturbations, initialized globally for same results in parallel
        gshape = self.domain.dist.grid_layout.global_shape(scales=1)
        slices = self.domain.dist.grid_layout.slices(scales=1)
        rand = np.random.RandomState(seed=23)
        noise = rand.standard_normal(gshape)[slices]

        # Linear background + perturbations damped at walls
        zb, zt = self.z_basis.interval
        pert =  1e-3 * noise * (zt - z) * (z - zb)
        b['g'] = -self.F*(z - pert)
        b.differentiate('z', out=bz)

        # Flow properties
        self.CFL = flow_tools.CFL(self.solver, initial_dt=1e-4, cadence=5, safety=0.5, max_change=1.5, min_change=0.5, max_dt=0.05)
        self.CFL.add_velocities(('u','v','w'))
        self.flow = flow_tools.GlobalFlowProperty(self.solver, cadence=10)
        self.flow.add_property("sqrt(u*u + v*v + w*w) / R", name='Re')

    def advance(self, T, u_init = None):
        # Main loop

        dt = self.CFL.compute_dt()
        init_time = self.solver.sim_time
        if u_init is not None:
            print("Reading input array in Python \n \n Array length: "+str(len(u_init))+"\n")
            if len(u_init) == self.Nx*self.Ny*self.Nz*self.scale*self.Nd: 
                for i, nx, ny, nz in itertools.product(range(self.Nd),range(self.Nx),range(self.Ny),range(self.Nz)):
                    if i == 0: self.u[nx,ny,nz] = u_init[i*(self.Ny+self.Nx+self.Nz)+ny*(self.Nx+self.Nz)+nx*self.Nz+nz]
                    elif i == 1: self.v[nx,ny,nz] = u_init[i*(self.Ny+self.Nx+self.Nz)+ny*(self.Nx+self.Nz)+nx*self.Nz+nz]
                    elif i == 2: self.w[nx,ny,nz] = u_init[i*(self.Ny+self.Nx+self.Nz)+ny*(self.Nx+self.Nz)+nx*self.Nz+nz]
        
        while init_time + T > self.solver.sim_time:
            print(self.solver.sim_time)
            self.solver.step(dt)
        # print("Length of output array in python is:"+str(len(self.u['g']))+"\n 2-Norm of the array is:"+str(np.linalg.norm(self.u['g'],2))+"\n")
        return(np.stack((self.solver.state['u']['g'],self.solver.state['v']['g'],self.solver.state['w']['g']),axis=3))

    def add(self,i):
        self.sum = self.sum + i

    def printSum(self):
        print(self.sum)
    

if __name__=="__main__":
    depy = DedalusPy()
    print(depy.advance(np.zeros(768)).shape)
# # Main loop
# end_init_time = time.time()
# logger.info('Initialization time: %f' %(end_init_time-start_init_time))
# try:
#     logger.info('Starting loop')
#     start_run_time = time.time()
#     while solver.ok:
#         dt = CFL.compute_dt()
#         solver.step(dt)
#         if (solver.iteration-1) % 100 == 0:
#             logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
#             logger.info('Max Re = %f' %flow.max('Re'))
# except:
#     logger.error('Exception raised, triggering end of main loop.')
#     raise
# finally:
#     end_run_time = time.time()
#     logger.info('Iterations: %i' %solver.iteration)
#     logger.info('Sim end time: %f' %solver.sim_time)
#     logger.info('Run time: %.2f sec' %(end_run_time-start_run_time))
#     logger.info('Run time: %f cpu-hr' %((end_run_time-start_run_time)/60/60*self.domain.dist.comm_cart.size))

