/**
 * Channelflow Dynamical System Interface
 *
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */

#include "dedalusdsi.h"
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define NPE_PY_ARRAY_OBJECT PyArrayObject

#include "cfbasics/mathdefs.h"
#include "channelflow/diffops.h"
#include "channelflow/poissonsolver.h"

#include <Python.h>
#include "numpy/arrayobject.h"
#include "channelflow/flowfield.h"
#include <Eigen/Dense>
#include "numpy/numpyconfig.h"

using namespace std;
using namespace Eigen;

namespace chflow {

dedalusDSI::dedalusDSI() {}

dedalusDSI::dedalusDSI(FieldSymmetry sigma, PoincareCondition* h, Real T,bool Tsearch, bool xrelative,
             bool zrelative, bool Tnormalize, Real Unormalize, const FlowField& u, ostream* os)
        :DSI(os),
        cfmpi_(u.cfmpi()),
        sigma_(sigma),
        h_(h),
        Tsearch_(Tsearch),
        xrelative_(xrelative),
        zrelative_(zrelative),
        //-----------------TODO------------------------
        Tinit_(T),
        // Tinit_(dnsflags_.T),
        axinit_(sigma.ax()),
        azinit_(sigma.az()),
        Tnormalize_(Tnormalize),
        Unormalize_(Unormalize),
        Nx_(u.Nx()),
        Ny_(u.Ny()),
        Nz_(u.Nz()),
        Nd_(u.Nd()),
        Lx_(u.Lx()),
        Lz_(u.Lz()),
        ya_(u.a()),
        yb_(u.b()),
        uunk_(u.taskid() == 0 ? xrelative_ + zrelative_ + Tsearch_ : 0),
        de_(init_dedalus()),
        muName_("None"),
        mu_(0){}

VectorXd dedalusDSI::eval(const VectorXd& x) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    Real T;
    extractVector(x, u, sigma_, T);
    FlowField Gu(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    // T = 0.2;
    G_dedalus(u, T, h_, sigma_, Gu, Tnormalize_, Unormalize_, de_, *os_);
    VectorXd Gx(VectorXd::Zero(x.rows()));
    field2vector(Gu, Gx);  // This does not change the size of Gx and automatically leaves the last entries zero
    return Gx;
}

VectorXd dedalusDSI::eval(const VectorXd& x0, const VectorXd& x1, bool symopt) {
    FlowField u0(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField u1(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    Real T0, T1;
    FieldSymmetry sigma0, sigma1;
    extractVector(x0, u0, sigma0, T0);
    extractVector(x1, u1, sigma1, T1);

    FlowField Gu(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);

    f(u0, T0, Gu, *os_, de_);
    Real funorm = L2Norm3d(Gu);
    if (symopt)
        Gu *= sigma0;
    Gu -= u1;

    if (Tnormalize_)
        Gu *= 1.0 / T0;
    if (Unormalize_ != 0.0)
        Gu *= 1. / sqrt(abs(funorm * (Unormalize_ - funorm)));

    VectorXd Gx(VectorXd::Zero(x0.rows()));

    field2vector(Gu, Gx);  // This does not change the size of Gx and automatically leaves the last entries zero

    return Gx;
}


//----------------------TODO----------------------------
void dedalusDSI::save(const VectorXd& x, const string filebase, const string outdir, const bool fieldsonly) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    extractVector(x, u, sigma, T);
    u.save(outdir + "u" + filebase);
    // dnsflags_.T = T;

    if (!fieldsonly) {
        string fs = fieldstats(u);
        if (u.taskid() == 0) {
            if (xrelative_ || zrelative_ || !sigma.isIdentity())
                sigma.save(outdir + "sigma" + filebase);
            if (Tsearch_)
                chflow::save(T, outdir + "T" + filebase);
            ofstream fout((outdir + "fieldconverge.asc").c_str(), ios::app);
            long pos = fout.tellp();
            if (pos == 0)
                fout << fieldstatsheader() << endl;
            fout << fs << endl;
            fout.close();
            // dnsflags_.save(outdir);
            // save_sp(T,outdir);
        }
    }
}

void dedalusDSI::saveEigenvec(const VectorXd& ev, const string label, const string outdir) {
    FlowField ef(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(ev, ef);
    ef *= 1.0 / L2Norm(ef);
    ef.save(outdir + "ef" + label);
}

void dedalusDSI::saveEigenvec(const VectorXd& evA, const VectorXd& evB, const string label1, const string label2,
                         const string outdir) {
    FlowField efA(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField efB(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(evA, efA);
    vector2field(evB, efB);
    Real c = 1.0 / sqrt(L2Norm2(efA) + L2Norm2(efB));
    efA *= c;
    efB *= c;
    efA.save(outdir + "ef" + label1);
    efB.save(outdir + "ef" + label2);
}

string dedalusDSI::stats(const VectorXd& x) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(x, u);
    return fieldstats_t(u, mu_);
}

pair<string, string> dedalusDSI::stats_minmax(const VectorXd& x) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField Gu(u);
    FieldSymmetry sigma;
    Real T;
    extractVector(x, u, sigma, T);

    vector<Real> stats = fieldstats_vector(u);
    vector<Real> minstats(stats);
    vector<Real> maxstats(stats);

    // quick hack to avoid new interface or creating simple f()
    // TimeStep dt = TimeStep(dnsflags_.dt, 0, 1, 1, 0, 0, false);
    // int fcount = 0;
    // PoincareCondition* h = 0;
    // Real CFL = 0.0;
    ostream muted_os(0);
    Real timep = T / 100.0;

    *os_ << "Using flag -orbOut: Calculate minmax-statistics of periodic orbit." << endl;
    for (int t = 0; t < 100; t++) {
        f(u, timep, Gu, muted_os, de_);
        stats = fieldstats_vector(u);
        for (uint i = 0; i < stats.size(); i++) {
            minstats[i] = (minstats[i] < stats[i]) ? minstats[i] : stats[i];
            maxstats[i] = (maxstats[i] > stats[i]) ? maxstats[i] : stats[i];
        }
        u = Gu;
    }
    // Return string
    stringstream smin;
    stringstream smax;
    smin << setw(8) << mu_;
    smax << setw(8) << mu_;
    for (uint i = 0; i < stats.size(); i++) {
        smin << setw(14) << minstats[i];
        smax << setw(14) << maxstats[i];
    }

    pair<string, string> minmax;
    minmax = make_pair(smin.str(), smax.str());
    return minmax;
}

string dedalusDSI::statsHeader() { return fieldstatsheader_t(muName_); }

void dedalusDSI::updateMu(Real mu) {
    DSI::updateMu(mu);
    mu_ = mu;    
    PyObject *method =  PyUnicode_FromString("updateMu");
    PyObject* mu_py = Py_BuildValue("d", mu);
    PyObject* muName_py = Py_BuildValue("s", muName_.c_str());
	PyObject_CallMethodObjArgs(de_, method, mu_py, muName_py, NULL);
}

void dedalusDSI::chooseMu(string muName) { muName_ = muName; }

string dedalusDSI::printMu() { return muName_; }

//-------------------TODO-----------------------------------
void dedalusDSI::saveParameters(string searchdir) { return; }

/// after finding new solution fix phases
void dedalusDSI::phaseShift(VectorXd& x) {
    FlowField unew(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FieldSymmetry sigma;
    Real T;
    extractVector(x, unew, sigma, T);
    // vector2field (x,unew);
    const int phasehackcoord = 0;  // Those values were fixed in continuesoln anyway
    const parity phasehackparity = Odd;
    const Real phasehackguess = 0.0;

    if (zphasehack_) {
        FieldSymmetry tau = zfixphasehack(unew, phasehackguess, phasehackcoord, phasehackparity);
        cout << "fixing z phase of potential solution with phase shift tau == " << tau << endl;
        unew *= tau;
    }
    if (xphasehack_) {
        FieldSymmetry tau = xfixphasehack(unew, phasehackguess, phasehackcoord, phasehackparity);
        cout << "fixing x phase of potential solution with phase shift tau == " << tau << endl;
        unew *= tau;
    }
    if (uUbasehack_) {

// ---------------TODO----------------------        
        // cout
        //     << "fixing u+Ubase decomposition so that <du/dy> = 0 at walls (i.e. Ubase balances mean pressure gradient))"
        //     << endl;
        // Real ubulk = Re(unew.profile(0, 0, 0)).mean();
        // if (abs(ubulk) < 1e-15)
        //     ubulk = 0.0;

        // ChebyCoeff Ubase =
        //     laminarProfile(dnsflags_.nu, dnsflags_.constraint, dnsflags_.dPdx, dnsflags_.Ubulk - ubulk, dnsflags_.Vsuck,
        //                    unew.a(), unew.b(), dnsflags_.ulowerwall, dnsflags_.uupperwall, unew.Ny());

        // fixuUbasehack(unew, Ubase);
    }
    makeVector(unew, sigma, T, x);
}

void dedalusDSI::phaseShift(MatrixXd& y) {
    FlowField unew(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    VectorXd yvec;
    FieldSymmetry sigma;
    Real T;

    const int phasehackcoord = 0;  // Those values were fixed in continuesoln anyway
    const parity phasehackparity = Odd;
    const Real phasehackguess = 0.0;

    FieldSymmetry taux(0.0, 0.0);
    FieldSymmetry tauz(0.0, 0.0);

    extractVector(y.col(0), unew, sigma, T);

    if (xphasehack_) {
        taux = xfixphasehack(unew, phasehackguess, phasehackcoord, phasehackparity);
        cout << "fixing x phase of potential solution with phase shift tau == " << taux << endl;
    }
    if (zphasehack_) {
        tauz = zfixphasehack(unew, phasehackguess, phasehackcoord, phasehackparity);
        cout << "fixing z phase of potential solution with phase shift tau == " << tauz << endl;
    }

    for (int i = 0; i < y.cols(); i++) {
        extractVector(y.col(i), unew, sigma, T);
        unew *= taux;
        unew *= tauz;
        makeVector(unew, sigma, T, yvec);
        y.col(i) = yvec;
    }
}

Real dedalusDSI::DSIL2Norm(const VectorXd& x) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(x, u);
    return L2Norm(u);
}

void dedalusDSI::makeVector(const FlowField& u, const FieldSymmetry& sigma, const Real T, VectorXd& x) {
    int taskid = u.taskid();
    int uunk = field2vector_size(u);                         // # of variables for u unknonwn
    const int Tunk = (Tsearch_ && taskid == 0) ? uunk : -1;  // index for T unknown
    const int xunk = (xrelative_ && taskid == 0) ? uunk + Tsearch_ : -1;
    const int zunk = (zrelative_ && taskid == 0) ? uunk + Tsearch_ + xrelative_ : -1;
    int Nunk = (taskid == 0) ? uunk + Tsearch_ + xrelative_ + zrelative_ : uunk;

    //   VectorXd x(Nunk);
    if (x.rows() < Nunk)
        x.resize(Nunk);
    field2vector(u, x);

    if (taskid == 0) {
        if (Tsearch_)
            x(Tunk) = T;
        if (xrelative_)
            x(xunk) = sigma.ax();
        if (zrelative_)
            x(zunk) = sigma.az();
    }
}

void dedalusDSI::extractVector(const VectorXd& x, FlowField& u, FieldSymmetry& sigma, Real& T) {
    int uunk = field2vector_size(u);  // Important for arclength method

    vector2field(x, u);
    const int Tunk = uunk + Tsearch_ - 1;
    const int xunk = uunk + Tsearch_ + xrelative_ - 1;
    const int zunk = uunk + Tsearch_ + xrelative_ + zrelative_ - 1;
    Real ax, az;

#ifdef HAVE_MPI
    if (u.taskid() == 0) {
#else
    {
#endif
        T = Tsearch_ ? x(Tunk) : Tinit_;
        ax = xrelative_ ? x(xunk) : axinit_;
        az = zrelative_ ? x(zunk) : azinit_;
    }

#ifdef HAVE_MPI
    MPI_Bcast(&T, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&az, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif
    sigma = FieldSymmetry(sigma_.sx(), sigma_.sy(), sigma_.sz(), ax, az, sigma_.s());
}

Real dedalusDSI::extractT(const VectorXd& x) {
    Real Tvec;
    FieldSymmetry sigma;
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    extractVector(x, u, sigma, Tvec);
    return Tvec;
}

Real dedalusDSI::extractXshift(const VectorXd& x) {
    Real Tvec;
    FieldSymmetry sigma;
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    extractVector(x, u, sigma, Tvec);
    return sigma.ax();
}

Real dedalusDSI::extractZshift(const VectorXd& x) {
    Real Tvec;
    FieldSymmetry sigma;
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    extractVector(x, u, sigma, Tvec);
    return sigma.az();
}

VectorXd dedalusDSI::xdiff(const VectorXd& a) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(a, u);
    VectorXd dadx(a.size());
    dadx.setZero();
    u = chflow::xdiff(u);
    field2vector(u, dadx);
    dadx *= 1. / L2Norm(dadx);
    return dadx;
}

VectorXd dedalusDSI::zdiff(const VectorXd& a) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    vector2field(a, u);
    VectorXd dadx(a.size());
    dadx.setZero();
    u = chflow::zdiff(u);
    field2vector(u, dadx);
    dadx *= 1. / L2Norm(dadx);
    return dadx;
}

VectorXd dedalusDSI::tdiff(const VectorXd& a, Real epsDt) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    FlowField edudtf(u);
    vector2field(a, u);
    f(u, epsDt, edudtf, *os_, de_);
    edudtf -= u;
    VectorXd dadt(a.size());
    field2vector(edudtf, dadt);
    dadt *= 1. / L2Norm(dadt);
    return dadt;
}

Real dedalusDSI::observable(VectorXd& x) {
    printout("computing mean dissipation", false);
    FlowField uarg(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    Real T;
    extractVector(x, uarg, sigma_, T);

    // just the current dissipation, was :
    // obs = meandissipation (u, T, dnsflags_, dt_, solntype);
    FlowField u(uarg);

    ChebyCoeff Ubase(u.Ny(), u.a(), u.b());

    //--------------------TODO----------------------
    // if (dnsflags_.baseflow == LinearBase)
    //     Ubase[1] = 1;
    // else if (dnsflags_.baseflow == ParabolicBase) {
    //     Ubase[0] = 0.5;
    //     Ubase[2] = -0.5;
    // } else if (dnsflags_.baseflow == LaminarBase) {
    //     Real ubulk = Re(uarg.profile(0, 0, 0)).mean();
    //     if (abs(ubulk) < 1e-15)
    //         ubulk = 0.0;
    //     Ubase =
    //         laminarProfile(dnsflags_.nu, dnsflags_.constraint, dnsflags_.dPdx, dnsflags_.Ubulk - ubulk, dnsflags_.Vsuck,
    //                        uarg.a(), uarg.b(), dnsflags_.ulowerwall, dnsflags_.uupperwall, uarg.Ny());
    // }
    printout(" of EQB/TW...", false);
    u += Ubase;
    printout("done");
    return dissipation(u);
}

Real dedalusDSI::tph_observable(VectorXd& x) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    Real T;
    extractVector(x, u, sigma_, T);
    return Ecf(u);
}

void dedalusDSI::save_array(FlowField& u){

    u.makePhysical();
	float* Array = new float[u.Nx()];
	for (int nx=0; nx<u.Nx(); ++nx){
        Array[nx] = u(nx, 2, 0, 2);
    }

    PyObject *py_array;
	float *ptr = Array;
	npy_intp dims[1] = { u.Nx() };
	py_array = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT, ptr);
	PyObject *method =  PyUnicode_FromString("save_array");
	PyObject_CallMethodObjArgs(de_, method, py_array, NULL);
    u.makeSpectral();

}

PyObject* init_dedalus(void){

    // setenv("PYTHONPATH", ".", 1);
	//Building python module and creating DedalusPy object
	cout<<"Creating DedalusPy object...\n";
	PyObject *module_name, *module, *dict, *python_class, *object;

	Py_Initialize();
    import_array();
    PyRun_SimpleString("from dedalus import public as de"); 
    cout<<"Python interpreter initialized...\n";
	module_name = PyUnicode_FromString("active_matter"); //dedalus_example_3D
    cout<<"Python module created...\n";
	module = PyImport_Import(module_name);
    if (module==NULL){
        cout<<"Failed to import!!!\n";
    }
    else
    {
        cout<<"Python module imported...\n";   
    }
	Py_DECREF(module_name);
	dict = PyModule_GetDict(module);
	Py_DECREF(module);
	python_class = PyDict_GetItemString(dict, "DedalusPy");
	Py_DECREF(dict);
	object = PyObject_CallObject(python_class, nullptr);
	Py_DECREF(python_class); 
    return object;    
}

void advance_dedalus(FlowField& u, Real T, PyObject* de){

    int Nx = u.Nx(), Ny = u.Ny(), Nz = u.Nz(), Nd = u.Nd();
    u.makePhysical(); 
	// Create an array to pass in dedalus
	int len_u = Nx*Ny*Nz*Nd;
    // int len_u = Nx;
	float* Array = new float[len_u];
    int ind = 0;
	for (int i=0; i<Nd; ++i){
		for (int ny=0; ny<Ny; ++ny){
			for (int nx=0; nx<Nx; ++nx){
				for (int nz=0; nz<Nz; ++nz){
                    Array[ind] = u(nx,ny,nz,i);
                    ind++;
				}
			}
		}
	}

    PyObject *py_array;
	float *ptr = Array;
	npy_intp dims[1] = { len_u };
	py_array = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT, ptr);
    PyObject* T_py = Py_BuildValue("d", T);
	PyObject *method =  PyUnicode_FromString("advance");

	PyObject *u_out = PyObject_CallMethodObjArgs(de, method, T_py, py_array, NULL);
    npy_intp u_c_ind[1]{0}; 

	double* u_c = reinterpret_cast<double*>(PyArray_GetPtr((PyArrayObject*)u_out, u_c_ind));
	//Convert this array back to Flowfield format
    ind = 0;
	for (int i=0; i<Nd; ++i){
		for (int ny=0; ny<Ny; ++ny){
			for (int nx=0; nx<Nx; ++nx){
				for (int nz=0; nz<Nz; ++nz){
                    u(nx,ny,nz,i) = u_c[ind];
                    ind++;
				}
			}
		}
	}

    u.makeSpectral(); 

    FlowField u_temp = u;
    fixdivnoslip(u_temp);
    if(L2Norm2_3d(u_temp-u)>1e-10){cout<<"Warning: flowfield is no longer divergence free: "<<L2Norm2_3d(u_temp-u)<<endl;}
    return; 
}

// return f^{N dt}(u) = time-(N dt) DNS integration of u
void f(const FlowField& u, Real T, FlowField& f_u, ostream& os, PyObject* de) {
    os << "f(u, N, dt, f_u, flags, dt) : " << flush;
    // vector<FlowField> fields(2);
    // fields[0] = u;
    // DNS dns(fields, flags);
    // dns.advance(fields, N);
    // f_u = fields[0];
    FlowField u_temp = u;
    // PyObject* de = init_dedalus();
    advance_dedalus(u_temp, T, de);
    f_u = u_temp;
    return;
}

// G(x) = G(u,sigma) = (sigma f^T(u) - u) for orbits
void G_dedalus(const FlowField& u, Real& T, PoincareCondition* h, const FieldSymmetry& sigma, FlowField& Gu,
       bool Tnormalize, Real Unormalize, PyObject* de,
       ostream& os) {
    // f(u, T, h, Gu, flags, dt, fcount, CFL, os);
    f(u, T, Gu, os, de); 
    Real funorm = L2Norm3d(Gu);
    Gu *= sigma;
    Gu -= u;
    if (Tnormalize)
        Gu *= 1.0 / T;
    if (Unormalize != 0.0) {
        Gu *= 1. / sqrt(abs(funorm * (Unormalize - funorm)));
    }
}


}  // namespace chflow
