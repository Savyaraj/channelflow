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

dedalusDSI::dedalusDSI(DNSFlags& dnsflags, FieldSymmetry sigma, PoincareCondition* h, TimeStep dt, bool Tsearch, bool xrelative,
             bool zrelative, bool Tnormalize, Real Unormalize, const FlowField& u, ostream* os)
        :cfDSI(dnsflags, sigma, h, dt, Tsearch, xrelative, zrelative, Tnormalize, Unormalize, u, os),
        de_(init_dedalus()){}

VectorXd dedalusDSI::eval(const VectorXd& x) {
    FlowField u(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    Real T;
    extractVector(x, u, sigma_, T);
    FlowField Gu(Nx_, Ny_, Nz_, Nd_, Lx_, Lz_, ya_, yb_, cfmpi_);
    G_dedalus(u, T, h_, sigma_, Gu, dnsflags_, dt_, Tnormalize_, Unormalize_, fcount_, CFL_,de_, *os_);
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

    f(u0, T0, Gu, dnsflags_, *os_, de_);
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

void dedalusDSI::updateMu(Real mu){
    cfDSI::updateMu(mu);
    
    PyObject *method =  PyUnicode_FromString("updateMu");
    PyObject* mu_py = Py_BuildValue("d", mu);
	PyObject_CallMethodObjArgs(de_, method, mu_py, NULL);
    return;
}

void dedalusDSI::save_array(FlowField& u){

    u.makePhysical();
    int len_u = u.Nx()*u.Ny()*u.Nz()*u.Nd();
	float* Array = new float[len_u];
    int ind = 0;
	for (int i=0; i<u.Nd(); ++i){
		for (int ny=0; ny<u.Ny(); ++ny){
			for (int nx=0; nx<u.Nx(); ++nx){
				for (int nz=0; nz<u.Nz(); ++nz){
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
void f(const FlowField& u, Real T, FlowField& f_u, const DNSFlags& flags_, ostream& os, PyObject* de) {
    os << "f(u, N, dt, f_u, flags, dt) : " << flush;
    DNSFlags flags(flags_);
    FlowField u_temp = u;
    advance_dedalus(u_temp, T, de);
    f_u = u_temp;
    return;
}

// G(x) = G(u,sigma) = (sigma f^T(u) - u) for orbits
void G_dedalus(const FlowField& u, Real& T, PoincareCondition* h, const FieldSymmetry& sigma, FlowField& Gu,
       const DNSFlags& flags, const TimeStep& dt, bool Tnormalize, Real Unormalize, int& fcount, Real& CFL, PyObject* de,
       ostream& os) {
    f(u, T, Gu, flags, os, de); 
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
