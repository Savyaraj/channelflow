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

        Real qc = 0.5;
        Real Lc = 2*pi/qc; //2*pi/qc
        Real Lx = 20*Lc; 
        Real gamma3 = 8.356;
        Real r = -0.002;

        FlowField u(100, 4, 4, 3, Lx, 1, 0, 1, cfmpi);
        int Nx = u.Nx(), Ny = u.Ny(), Nz = u.Nz(), Nd = u.Nd();
        
        // Real u_init[Nx];
        // ifstream file("u_init.txt");
        // if(file.is_open()){
        //     for (int i=0; i<Nx; ++i){
        //         file>>u_init[i];
        //     }
        // }

        u.makePhysical();
        for (int i=0; i<Nd; ++i){
            for (int ny=0; ny<Ny; ++ny){
                for (int nx=0; nx<Nx; ++nx){
                    for (int nz=0; nz<Nz; ++nz){
                        Real x = nx*Lx/Nx - Lx/2;
                        // if (i==2 && ny==2){u(nx, ny, nz, i) =  0.1 + 0.2*cos(2*pi*x/Lc);}//0.4870017730666727;}// + 0.01*(float)rand()/RAND_MAX;}
                        // if (i==2 && ny==2){u(nx, ny, nz, i) = 2*sqrt(-2*r/gamma3)*(1/cosh(x*sqrt(-r)/(2*qc)))*(1+0.01*(float)rand()/RAND_MAX)*(0.15+cos(qc*x+pi));}//0.4870017730666727;}// + 0.01*(float)rand()/RAND_MAX;}
                        if (i==2 && ny==2){u(nx, ny, nz, i) = (1+0.0*rand()/RAND_MAX)*sin(2*pi*x/(5*Lc));}
                        // if (i==2 && ny==2){u(nx, ny, nz, i) = u_init[nx];}
                        // if (i==2 && ny==2){u(nx, ny, nz, i) = 0.34473189;}// + 0.01*(float)rand()/RAND_MAX;}
                        else{u(nx, ny, nz, i) = 0.0;}
                        if(i==2 && ny==2 && nz==0){cout<<u(nx,ny,nz,i)<<endl;} 
                    }
                }
            }
        }
        u.makeSpectral();
        u.writeNetCDF("./sample.nc");

        // dsi->save_array(u);


        // u.makeSpectral();
        // for (int nx = 0 ; nx < u.Mx(); nx++){
        //     // u.cmplx(nx,0,0,0) = 0.6278+0.0*I;
        //     cout<<u.cmplx(nx,1,0,0)<<endl;
        // } 
        // fixdivnoslip(u);
        // u.makeSpectral();
        // VectorXd x;
        // field2vector(u, x);
        // for (int i = 0; i< x.size(); i++){
        //     if(x(i) != 0){cout<<i<<"\t"<<x(i)<<endl;}
        // }
        // cout<<L2Norm2(x)<<endl;
        // vector2field(x, u);
        // u.makeSpectral();
        // FlowField u_temp(100, 6, 6, 3, 1, 1, 0, 1, cfmpi);
        // u_temp.makePhysical();
        // for (int nx = 0; nx<u_temp.Nx(); nx++){
        //     u_temp(nx,0,0,0) = 0.6278; 
        // }
        // vector2field(x,u_temp);
        // u_temp.makePhysical();
        // for (int nx = 0 ; nx < u.Mx(); nx++){
        //     cout<<u.cmplx(nx,1,0,0)<<endl;
        // }
        // u.makePhysical();
        // cout<<u.Lx()<<u.Ly()<<u.Lz();
        // FlowField uopt("./ubest",NULL);
        // uopt.makePhysical();
        // cout<<"u is:"<<uopt(0,0,0,0)<<endl;
        // VectorXd x; 
        // field2vector(u, x);

        // Real norm = L2Norm2(x);
        // cout<<"2-norm of the equilibrium solution is: "<<norm<<endl;
        // for (int i=0; i<10;i++){
        //     cout<<x[i]<<endl;
        // }

        // FlowField u1(1024, 10, 10, 3, 1, 1, 0, 1, cfmpi);

        // Real T=0.001;

        // u.makePhysical();
        // for (int nx = 0; nx <u.Nx(); ++nx){
        //     cout<<u(nx,1,0,0)<<endl;
        // }

        // cout<<u.Nx()<<"\t"<<u.Ny()<<"\t"<<u.Nz()<<endl;
        // u.makeSpectral();
        // for (int i=0; i<u.Nd(); ++i){
        //     for (int ny=0; ny<u.My(); ++ny){
        //         for (int nx=0; nx<u.Mx(); ++nx){
        //             for (int nz=0; nz<u.Mz(); ++nz){
        //                 if (Re(u.cmplx(nx,ny,nz,i))!=0 || Im(u.cmplx(nx,ny,nz,i))!=0){
        //                     cout<<"mx: "<<nx<<" my: "<<ny<<" mz: "<<nz<<" i: "<<i<<" u: "<<u.cmplx(nx,ny,nz,i)<<endl;
        //                     }
        //             }
        //         }
        //     }
        // }

        // cout<<"After conversion: "<<endl;
        // VectorXd x;

        // field2vector(u, x);
        // vector2field(x, u1);
        // for (int i=0; i<u.Nd(); ++i){
        //     for (int ny=0; ny<u.My(); ++ny){
        //         for (int nx=0; nx<u.Mx(); ++nx){
        //             for (int nz=0; nz<u.Mz(); ++nz){
        //                 if (Re(u.cmplx(nx,ny,nz,i))!=0 || Im(u.cmplx(nx,ny,nz,i))!=0){
        //                     cout<<"mx: "<<nx<<" my: "<<ny<<" mz: "<<nz<<" i: "<<i<<" u: "<<u.cmplx(nx,ny,nz,i)<<endl;
        //                     }
        //             }
        //         }
        //     }
        // }
        // for (int nx = 0; nx <u.Mx(); ++nx){
        //     cout<<u.cmplx(nx,1,0,0)<<endl;
        // }
        // u.makePhysical();
        // for (int nx = 0; nx <u.Nx(); ++nx){
        //     cout<<u(nx,2,0,0)<<endl;
        // }

        // VectorXd y = dsi->eval(x);

        // VectorXd z = dsi->eval(y);
        // cout<<"Output vector size is: "<<y.size()<<endl;

        // u1.makeSpectral();
        // FlowField u2(1000, 4, 4, 3, Lx, 1, 0, 1, cfmpi);
        // VectorXd x;
        // field2vector(u1, x);
        // vector2field(x, u2);
        // u1.makePhysical();
        // u2.makePhysical();
        // // int Nx = u1.Nx(), Ny = u1.Ny(), Nz = u1.Nz(), Nd = u1.Nd();
        // for (int i=0; i<Nd; ++i){
        //     for (int ny=0; ny<Ny; ++ny){
        //         for (int nx=0; nx<Nx; ++nx){
        //             for (int nz=0; nz<Nz; ++nz){
        //                 if((u1(nx,ny,nz,i) - u2(nx,ny,nz,i))*(u1(nx,ny,nz,i) - u2(nx,ny,nz,i))>1e-10){
        //                     cout<<u1(nx,ny,nz,i)<<'\t'<<u2(nx,ny,nz,i)<<endl;
        //                 }
        //                 // if(u2(nx,ny,nz,i)>1e-10){cout<<u2(nx,ny,nz,i)<<endl;}
        //             }
        //         }
        //     }
        // }
        
    }
}
