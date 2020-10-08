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

#include <Python.h>
#include "numpy/arrayobject.h"
#include "channelflow/flowfield.h"
#include <Eigen/Dense>
#include "numpy/numpyconfig.h"


namespace chflow {

enum class deParameter {
    T,
    Re,
    P,
    Ub,
    Uw,
    ReP,
    Theta,
    ThArc,
    ThLx,
    ThLz,
    Lx,
    Lz,
    Aspect,
    Diag,
    Lt,
    Vs,
    ReVs,
    H,
    HVs,
    Rot
};

// converts the string from "fieldstats" in diffops to a vector of Reals
std::vector<Real> fieldstats_vector(const FlowField& u);

class dedalusDSI : public DSI {
   public:
   
    /** \brief default constructor */
    dedalusDSI();
    virtual ~dedalusDSI() {}

    /** \brief Initialize dedalusDSI */
    dedalusDSI(DNSFlags& dnsflags, FieldSymmetry sigma, PoincareCondition* h, TimeStep dt, bool Tsearch, bool xrelative,
          bool zrelative, bool Tnormalize, Real Unormalize, const FlowField& u, std::ostream* os = &std::cout);
    Eigen::VectorXd eval(const Eigen::VectorXd& x);
    Eigen::VectorXd eval(const Eigen::VectorXd& x0, const Eigen::VectorXd& x1, bool symopt);
    void save(const Eigen::VectorXd& x, const std::string filebase, const std::string outdir = "./",
              const bool fieldsonly = false) override;
    void saveEigenvec(const Eigen::VectorXd& x, const std::string label, const std::string outdir) override;
    void saveEigenvec(const Eigen::VectorXd& x1, const Eigen::VectorXd& x2, const std::string label1,
                      const std::string label2, const std::string outdir) override;

    Real DSIL2Norm(const Eigen::VectorXd& x) override;
    std::string stats(const Eigen::VectorXd& x) override;
    std::pair<std::string, std::string> stats_minmax(const Eigen::VectorXd& x) override;
    std::string statsHeader() override;
    void makeVector(const FlowField& u, const FieldSymmetry& sigma, const Real T, Eigen::VectorXd& x);
    void extractVector(const Eigen::VectorXd& x, FlowField& u, FieldSymmetry& sigma, Real& T);
    void toVector(const std::vector<FlowField>& u, const FieldSymmetry& sigma, const Real T, Eigen::VectorXd& x){};

    /// \name Compute derivatives of FlowField corresponding to this vector
    Eigen::VectorXd xdiff(const Eigen::VectorXd& a) override;
    Eigen::VectorXd zdiff(const Eigen::VectorXd& a) override;
    Eigen::VectorXd tdiff(const Eigen::VectorXd& a, Real epsDt) override;

    /// \name Handle continuation parameter
    void updateMu(Real mu) override;
    void chooseMu(std::string muName);
    void chooseMu(deParameter mu);
    std::string printMu() override;  // document
    void saveParameters(std::string searchdir) override;
    deParameter s2cPar(std::string muName);
    std::string cPar2s(deParameter cPar);
    void phaseShift(Eigen::VectorXd& x) override;
    void phaseShift(Eigen::MatrixXd& y) override;
    inline void setPhaseShifts(bool xphasehack, bool zphasehack, bool uUbasehack);
    Real observable(Eigen::VectorXd& x) override;

    Real tph_observable(Eigen::VectorXd& x) override;
    Real extractT(const Eigen::VectorXd& x) override;
    Real extractXshift(const Eigen::VectorXd& x) override;
    Real extractZshift(const Eigen::VectorXd& x) override;

    Real getCFL() const { return CFL_; };
    bool XrelSearch() const override { return xrelative_; };
    bool ZrelSearch() const override { return zrelative_; };
    bool Tsearch() const override { return Tsearch_; };

    void save_array(FlowField& u);

   protected:
    DNSFlags dnsflags_;
    CfMPI* cfmpi_;
    FieldSymmetry sigma_;
    PoincareCondition* h_;
    TimeStep dt_;
    bool Tsearch_;
    bool xrelative_;
    bool zrelative_;
    Real Tinit_;
    Real axinit_;
    Real azinit_;
    bool Tnormalize_;
    Real Unormalize_;
    int fcount_;
    int Nx_;
    int Ny_;
    int Nz_;
    int Nd_;
    Real Lx_;
    Real Lz_;
    Real ya_;
    Real yb_;
    Real CFL_;
    int uunk_;
    bool xphasehack_;
    bool zphasehack_;
    bool uUbasehack_;

    deParameter cPar_ = deParameter::T;
    
    private:
      PyObject* de_;
};

inline void dedalusDSI::setPhaseShifts(bool xphasehack, bool zphasehack, bool uUbasehack) {
    xphasehack_ = xphasehack;
    zphasehack_ = zphasehack;
    uUbasehack_ = uUbasehack;
}

PyObject* init_dedalus(void);

void advance_dedalus(FlowField& u, Real T, PyObject* de); 

void f(const FlowField& u, Real T, FlowField& f_u, const DNSFlags& flags_, std::ostream& os, PyObject* de);

void G_dedalus(const FlowField& u, Real& T, PoincareCondition* h, const FieldSymmetry& sigma, FlowField& Gu,
const DNSFlags& flags, const TimeStep& dt, bool Tnormalize, Real Unormalize, int& fcount, Real& CFL, PyObject*de,
std::ostream& os = std::cout);

}  // namespace chflow

#endif
