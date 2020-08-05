#include <Python.h>
#include <iostream>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <vector>
#include <cstdlib> // setenv

namespace np = boost::python::numpy;
namespace python = boost::python;
using namespace std;

int main()
{
  // Allow Python to load modules from the current directory.
  setenv("PYTHONPATH", ".", 1);
  // Initialize Python.
  Py_Initialize();
//   np::initialize();
  
  try
  {
    // >>> import dedalus_example
    python::object my_python_class_module = python::import("dedalus_example");

    // >>> create instance of the DedalusPy class
    python::object deda = my_python_class_module.attr("DedalusPy")();

    // >>> call advance method from dedalus 
    python::object u = deda.attr("advance")();
    
    cout<<"Length of array in C++ is: "<<python::len(u)<<endl;
    std::string object_classname = boost::python::extract<std::string>(u.attr("__class__").attr("__name__"));
    std::cout<<"Datatype of the array in C++ is: "<<object_classname<<std::endl;

    cout<<"The first 5 elements of the output are: "<<endl;
    for (int i=0; i<5; i++){
		std::cout<<python::extract<double>(u[i])<<std::endl;
	}
  }
  catch (const python::error_already_set&)
  {
    PyErr_Print();
    return 1;
  }

  // Do not call Py_Finalize() with Boost.Python.
}

