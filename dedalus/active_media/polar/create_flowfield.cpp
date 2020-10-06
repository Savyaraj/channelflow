/**
 * Script to verify Dynamical Systems Interface for dedalus
 */

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include<channelflow/dedalusdsi.h>
#include "channelflow/flowfield.h"
#include "nsolver/nsolver.h"
#include <Eigen/Dense>
#include "channelflow/diffops.h"

#include <fstream>
#include <Python.h>
#include "numpy/arrayobject.h"
#include "numpy/numpyconfig.h"

using namespace chflow;
using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]) {
    cfMPI_Init(&argc, &argv);
    {
        ArgList args(argc, argv, "find an invariant solution using Newton-Krylov-hookstep algorithm");

        /** Choose the Newton algorithm to be used. Currently, two options are available: simple Newton without any
         * trust region optimization, and Newton with Hookstep (default). For the simple Newton, you can choose either a
         * full-space algorithm to solve the Newton equations (-solver "eigen") or between the two iterative algorithms
         * GMRES and BiCGStab. Newton-Hookstep requires GMRES. Note that the available parameters depend on your choice
         * of the algorithm.
         */

        unique_ptr<Newton> N;
        NewtonSearchFlags searchflags(args);
        searchflags.save(searchflags.outdir);
        N = unique_ptr<Newton>(new NewtonAlgorithm(searchflags));

        DNSFlags dnsflags(args, searchflags.laurette);
        TimeStep dt(dnsflags);

        bool Rxsearch, Rzsearch, Tsearch;
        Rxsearch = searchflags.xrelative;
        Rzsearch = searchflags.zrelative;
        Tsearch = searchflags.solntype == PeriodicOrbit ? true : false;

        const bool Tnormalize = (Tsearch || searchflags.laurette) ? false : true;

        /** Read in remaining arguments */

        args.section("Program options");
        const string sigmastr =
            args.getstr("-sigma", "--sigma", "", "file containing sigma of sigma f^T(u) - u = 0 (default == identity)");
        const Real unormalize = args.getreal("-un", "--unormalize", 0.0, "lower bound in energy for search");
        const int nproc0 =
            args.getint("-np0", "--nproc0", 0, "number of MPI-processes for transpose/number of parallel ffts");
        const int nproc1 = args.getint("-np1", "--nproc1", 0, "number of MPI-processes for one fft");
        const bool msinit =
            args.getflag("-MSinit", "--MSinitials", "read different files as the initial guesses for different shoots");
        const string uname = args.getstr(1, "<flowfield>", "initial guess for the solution");

        args.check();
        args.save();
        WriteProcessInfo(argc, argv);
        dnsflags.save();
        cout << dnsflags << endl;

        CfMPI* cfmpi = &CfMPI::getInstance(nproc0, nproc1);

        FlowField u1(uname, cfmpi);

        FlowField fu(uname, cfmpi);

        FieldSymmetry sigma;
        if (sigmastr.length() != 0)
            sigma = FieldSymmetry(sigmastr);

        /** Construct the dynamical-systems interface object depending on the given parameters. Current options are
         * either standard (f(u) via forward time integration) or Laurette (f(u) via Laurettes method)
         */
        unique_ptr<dedalusDSI> dsi;
        dsi = unique_ptr<dedalusDSI>(new dedalusDSI(dnsflags, sigma, 0, dt, Tsearch, Rxsearch, Rzsearch, Tnormalize, unormalize,
                                          u1, N->getLogstream()));

        Real L = 2*pi; 
        Real beta = 0.5;
        Real alpha = -0.5;
        int Nxz = 128;
        FlowField u(Nxz, 4, Nxz, 3, L, L, 0, 1, cfmpi);
        int Nx = u.Nx(), Ny = u.Ny(), Nz = u.Nz(), Nd = u.Nd();

        u.makePhysical();
        for (int i=0; i<Nd; ++i){
            for (int ny=0; ny<Ny; ++ny){
                for (int nx=0; nx<Nx; ++nx){
                    for (int nz=0; nz<Nz; ++nz){
                        Real x = nx*L/Nx - L/2;
                        if (i==0 && ny==2){u(nx, ny, nz, i) = 1.0;}//0.4870017730666727;}// + 0.01*}
                        // if (i==0 && ny==2){u(nx, ny, nz, i) = sqrt(-alpha/beta);}//0.4870017730666727;}// + 0.01*(float)rand()/RAND_MAX;}
                        // else if (i==2 && ny==2){u(nx, ny, nz, i) = -x;}//(float)rand()/RAND_MAX;}//0.4870017730666727;}// + 0.01*(float)rand()/RAND_MAX;}
                        // if (i==2 && ny==2){u(nx, ny, nz, i) = 2*sqrt(-2*r/gamma3)*(1/cosh(x*sqrt(-r)/(2*qc)))*(0.15+cos(qc*x))+0.005*(float)rand()/RAND_MAX;}//0.4870017730666727;}// + 0.01*(float)rand()/RAND_MAX;}
                        
                        else{u(nx, ny, nz, i) = 0.0;}
                        // if (i==0 && ny==2 && nx==5 && nz==5){cout<<"u_create:"<<u(nx,ny,nz,i)<<endl;}
                        // if(i==0 && ny==2 && nz==0){cout<<u(nx,ny,nz,i)<<endl;} 
                    }
                }
            }
        }
        u.makeSpectral();
        u.writeNetCDF("./sample.nc");

        VectorXd x;
        field2vector(u, x);
        VectorXd y = dsi->eval(x);

        // u.makeSpectral();
        // FlowField u2(Nxz, 4, Nxz, 3, L, L, 0, 1, cfmpi);
        // VectorXd x;
        // field2vector(u, x);
        // vector2field(x, u2);
        // u.makePhysical();
        // u2.makePhysical();
        // // int Nx = u1.Nx(), Ny = u1.Ny(), Nz = u1.Nz(), Nd = u1.Nd();
        // for (int i=0; i<Nd; ++i){
        //     for (int ny=0; ny<Ny; ++ny){
        //         for (int nx=0; nx<Nx; ++nx){
        //             for (int nz=0; nz<Nz; ++nz){
        //                 if((u(nx,ny,nz,i) - u2(nx,ny,nz,i))*(u(nx,ny,nz,i) - u2(nx,ny,nz,i))>1e-6){
        //                     cout<<u(nx,ny,nz,i)<<'\t'<<u2(nx,ny,nz,i)<<endl;
        //                 }
        //                 // if(u2(nx,ny,nz,i)>1e-10){cout<<u2(nx,ny,nz,i)<<endl;}
        //             }
        //         }
        //     }
        // }
        
    }
}
