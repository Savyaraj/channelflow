/**
 * Dynamical System Interface (DSI) for communicating with dedalus 
 *
 * This file is a part of channelflow version 2.0, https://channelflow.ch .
 * License is GNU GPL version 2 or later: ./LICENSE
 */// #ifdef NPY_1_7_API_VERSION
// #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#ifndef DEDALUSDSI_H
#define DEDALUSDSI_H
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define NPE_PY_ARRAY_OBJECT PyArrayObject

#include <sys/stat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "cfbasics/cfvector.h"
#include "channelflow/cfmpi.h"
#include "channelflow/chebyshev.h"
#include "channelflow/dns.h"
#include "channelflow/flowfield.h"
#include "channelflow/periodicfunc.h"
#include "channelflow/symmetry.h"
#include "channelflow/tausolver.h"
#include "channelflow/utilfuncs.h"
#include "nsolver/nsolver.h"
#include "channelflow/cfdsi.h"

#include <Python.h>
#include "numpy/arrayobject.h"
#include "channelflow/flowfield.h"
#include <Eigen/Dense>
#include "numpy/numpyconfig.h"


namespace chflow {

class dedalusDSI : public cfDSI {
   public:
   
    /** \brief default constructor */
    dedalusDSI();
    virtual ~dedalusDSI() {}

    /** \brief Initialize dedalusDSI */
    dedalusDSI(DNSFlags& dnsflags, FieldSymmetry sigma, PoincareCondition* h, TimeStep dt, bool Tsearch, bool xrelative,
          bool zrelative, bool Tnormalize, Real Unormalize, const FlowField& u, std::ostream* os = &std::cout);
    Eigen::VectorXd eval(const Eigen::VectorXd& x);
    Eigen::VectorXd eval(const Eigen::VectorXd& x0, const Eigen::VectorXd& x1, bool symopt);

    void updateMu(Real mu); 
    
    void save_array(FlowField& u);
    
    private:
      PyObject* de_;
};

    PyObject* init_dedalus(void);
    
    void advance_dedalus(FlowField& u, Real T, PyObject* de); 

    void f(const FlowField& u, Real T, FlowField& f_u, const DNSFlags& flags_, std::ostream& os, PyObject* de);

    void G_dedalus(const FlowField& u, Real& T, PoincareCondition* h, const FieldSymmetry& sigma, FlowField& Gu,
      const DNSFlags& flags, const TimeStep& dt, bool Tnormalize, Real Unormalize, int& fcount, Real& CFL, PyObject*de,
      std::ostream& os = std::cout);

}  // namespace chflow

#endif
