#include <stdio.h>
#include <iostream>
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "numpy/arrayobject.h"
#include "channelflow/flowfield.h"
#include <Eigen/Dense>
#include "numpy/numpyconfig.h"
// #ifdef NPY_1_7_API_VERSION
// #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
// #define NPE_PY_ARRAY_OBJECT PyArrayObject
// #else
// #define NPE_PY_ARRAY_OBJECT PyObject
// #endif

int main() {

	setenv("PYTHONPATH", ".", 1);

	//Building python module and creating DedalusPy object
	std::cout<<"Creating DedalusPy object...\n";
	PyObject *module_name, *module, *dict, *python_class, *object;

	Py_Initialize();
	import_array()
	module_name = PyUnicode_FromString("swift_hohenberg");
	module = PyImport_Import(module_name);
	Py_DECREF(module_name);
	dict = PyModule_GetDict(module);
	Py_DECREF(module);
	python_class = PyDict_GetItemString(dict, "DedalusPy");
	Py_DECREF(dict);
	object = PyObject_CallObject(python_class, nullptr);
	Py_DECREF(python_class);

	//Creating flowfield object
	std::cout<<"Creating flowfield...\n";
	const int Nx = 1024, Ny = 1, Nz = 1, Nd = 1;
    const double Lx = 1.0, Lz = 1.0;
    const double a = 0.0, b = 1.0;
	chflow::FlowField u_init(Nx, Ny, Nz, Nd, Lx, Lz, a, b, NULL); //, chflow::Physical, chflow::Physical);
    // chflow::FlowField u_init("newfield",NULL); 
    std::cout<<"Flowfield created!"<<std::endl;
	u_init.makePhysical();

	//Converting flowfield to vector object for further manipulation
	// std::cout<<"Converting flowfield into vector...\n";
	// Eigen::VectorXd v(3);
	// chflow::field2vector(u_init, v);
	// std::cout<<"Size of converted vector is: "<<v.size()<<std::endl;
	// Eigen::VectorXd v(Nx,Ny,Nz,Nd);

	// Create an array to pass in dedalus
	int len_u = Nx*Ny*Nz*Nd;
	float Array [len_u];
	for (int i=0; i<Nd; ++i){
		for (int ny=0; ny<Ny; ++ny){
			for (int nx=0; nx<Nx; ++nx){
				for (int nz=0; nz<Nz; ++nz){
					Array[i*(Ny+Nx+Nz)+ny*(Nx+Nz)+nx*Nz+nz] = u_init(nx, ny, nz, i);
				}

			}
		}

	}
	
	PyObject *py_array;
	float *ptr = Array;
	npy_intp dims[1] = { len_u };
	// double *ptr = v.data();
	// npy_intp dims[1] = { v.size() };
	py_array = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT, ptr);
	double T = 10;
	PyObject* T_py = Py_BuildValue("d", T);
	//Test the instance method "advance" for timestepping and print output
	// PyObject *u = PyObject_CallMethod(object, "advance","()");
	PyObject *method =  PyUnicode_FromString("advance");
	PyObject *u = PyObject_CallMethodObjArgs(object, method, T_py, py_array, NULL);

	const char* dt_type = Py_TYPE(u)->tp_name;
	std::cout<<"The data type of output is: "<<dt_type<<std::endl;

	int dim = PyArray_Size(u);
	std::cout<<"The size of output array is: "<<dim<<std::endl;

	npy_intp u_c_ind[4]{0,0,0,0}; 
	double* u_c = reinterpret_cast<double*>(PyArray_GetPtr((PyArrayObject*)u, u_c_ind));

    // double* u_c = reinterpret_cast<double*>(PyArray_DATA(contig));
	for (int i=0; i<10; i++){
		std::cout<<u_c[i,0,0,0]<<std::endl;
	}

	//Convert this array back to Flowfield format
	std::cout<<"Converting output array to flowfield...\n \n";
	// Eigen::VectorXd v_out(dim);
	// for (int i=0; i<10; i++){
	// 	v_out(i) = u_c[i];
	// }
	chflow::FlowField u_out(Nx, Ny, Nz, Nd, Lx, Lz, a, b, NULL, chflow::Physical, chflow::Physical);
	// chflow::vector2field(v_out, u_out);

	for (int i=0; i<Nd; ++i){
		for (int ny=0; ny<Ny; ++ny){
			for (int nx=0; nx<Nx; ++nx){
				for (int nz=0; nz<Nz; ++nz){
					u_out(nx, ny, nz, i) = u_c[nx,ny,nz,i];
				}

			}
		}

	}

	u_out.makeSpectral();

	std::cout<<"Size of output array is "<<field2vector_size(u_out)<<std::endl;

	std::cout<<"Converted!! \n \n";

	std::getchar();
	Py_Finalize();

	return (0);
}