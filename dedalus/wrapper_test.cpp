#include <stdio.h>
#include <iostream>
#include <Python.h>
#include "numpy/arrayobject.h"

#include "numpy/numpyconfig.h"
#ifdef NPY_1_7_API_VERSION
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define NPE_PY_ARRAY_OBJECT PyArrayObject
#else
#define NPE_PY_ARRAY_OBJECT PyObject
#endif

int main() {

	setenv("PYTHONPATH", ".", 1);

	PyObject *module_name, *module, *dict, *python_class, *object;

	// Initializes the Python interpreter
	Py_Initialize();
	import_array()
	module_name = PyUnicode_FromString(
		"dedalus_example");

	// Load the module object
	module = PyImport_Import(module_name);
	if (module == nullptr) {
		PyErr_Print();
		std::cerr << "Fails to import the module.\n";
	return 1;
	}
	else{
		std::cout<<"Module import successful.\n";
	}
	
	Py_DECREF(module_name);

	// dict is a borrowed reference.
	dict = PyModule_GetDict(module);
	if (dict == nullptr) {
		PyErr_Print();
		std::cerr << "Fails to get the dictionary.\n";
		return 1;
	}
	else{
		std::cout<<"Dictionary retrieval successful.\n";
	}
	
	Py_DECREF(module);

	// Builds the name of a callable class
	python_class = PyDict_GetItemString(dict, "DedalusPy");
	if (python_class == nullptr) {
		PyErr_Print();
		std::cerr << "Fails to get the Python class.\n";
		return 1;
	}
	else {
		std::cout<<"Python class build successful.\n";
	}
	
	Py_DECREF(dict);

	// Creates an instance of the class
	if (PyCallable_Check(python_class)) {
		object = PyObject_CallObject(python_class, nullptr);
		Py_DECREF(python_class);
		std::cout<<"Class object generation successful.\n";
	} 
	else {
		std::cout << "Cannot instantiate the Python class" << std::endl;
		Py_DECREF(python_class);
	return 1;
	}

	//Test the instance method "advance" for timestepping and print output
	PyObject *u = PyObject_CallMethod(object, "advance","()"); 

	const char* dt_type = Py_TYPE(u)->tp_name;
	std::cout<<"The data type of output is: "<<dt_type<<std::endl;

	int dim = PyArray_Size(u);
	std::cout<<"The size of output array is: "<<dim<<std::endl;

    double* u_c = reinterpret_cast<double*>(PyArray_DATA(u));
	for (int i=0; i<10; i++){
		std::cout<<u_c[i]<<std::endl;
	}

	std::getchar();
	Py_Finalize();

	return (0);
}