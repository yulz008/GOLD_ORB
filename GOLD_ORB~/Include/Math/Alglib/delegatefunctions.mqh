//+------------------------------------------------------------------+
//|                                            delegatefunctions.mqh |
//|            Copyright 2003-2012 Sergey Bochkanov (ALGLIB project) |
//|                   Copyright 2012-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
//| Implementation of ALGLIB library in MetaQuotes Language 5        |
//|                                                                  |
//| The features of the library include:                             |
//| - Linear algebra (direct algorithms, EVD, SVD)                   |
//| - Solving systems of linear and non-linear equations             |
//| - Interpolation                                                  |
//| - Optimization                                                   |
//| - FFT (Fast Fourier Transform)                                   |
//| - Numerical integration                                          |
//| - Linear and nonlinear least-squares fitting                     |
//| - Ordinary differential equations                                |
//| - Computation of special functions                               |
//| - Descriptive statistics and hypothesis testing                  |
//| - Data analysis - classification, regression                     |
//| - Implementing linear algebra algorithms, interpolation, etc.    |
//|   in high-precision arithmetic (using MPFR)                      |
//|                                                                  |
//| This file is free software; you can redistribute it and/or       |
//| modify it under the terms of the GNU General Public License as   |
//| published by the Free Software Foundation (www.fsf.org); either  |
//| version 2 of the License, or (at your option) any later version. |
//|                                                                  |
//| This program is distributed in the hope that it will be useful,  |
//| but WITHOUT ANY WARRANTY; without even the implied warranty of   |
//| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the     |
//| GNU General Public License for more details.                     |
//+------------------------------------------------------------------+
#include <Object.mqh>
#include "matrix.mqh"
//+------------------------------------------------------------------+
//| Calculates f(arg), stores result to func                         |
//+------------------------------------------------------------------+
class CNDimensional_Func
  {
public:
                     CNDimensional_Func(void);
                    ~CNDimensional_Func(void);

   virtual void      Func(double &x[],double &func,CObject &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNDimensional_Func::CNDimensional_Func(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNDimensional_Func::~CNDimensional_Func(void)
  {

  }
//+------------------------------------------------------------------+
//| Empty function body                                              |
//+------------------------------------------------------------------+
void CNDimensional_Func::Func(double &x[],double &func,CObject &obj)
  {

  }
//+------------------------------------------------------------------+
//| calculates func = f(arg), grad[i] = df(arg)/d(arg[i])            |
//+------------------------------------------------------------------+
class CNDimensional_Grad
  {
public:
   //--- constructor, destructor
                     CNDimensional_Grad(void);
                    ~CNDimensional_Grad(void);
   //--- virtual method
   virtual void      Grad(double &x[],double &func,double &grad[],CObject &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNDimensional_Grad::CNDimensional_Grad(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNDimensional_Grad::~CNDimensional_Grad(void)
  {

  }
//+------------------------------------------------------------------+
//| Empty function body                                              |
//+------------------------------------------------------------------+
void CNDimensional_Grad::Grad(double &x[],double &func,double &grad[],
                              CObject &obj)
  {

  }
//+------------------------------------------------------------------+
//| Calculates func = f(arg), grad[i] = df(arg)/d(arg[i]),           |
//| hess[i,j] = d2f(arg)/(d(arg[i])*d(arg[j]))                       |
//+------------------------------------------------------------------+
class CNDimensional_Hess
  {
public:
   //--- constructor, destructor
                     CNDimensional_Hess(void);
                    ~CNDimensional_Hess(void);
   //--- virtual method
   virtual void      Hess(double &x[],double &func,double &grad[],CMatrixDouble &hess,CObject &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNDimensional_Hess::CNDimensional_Hess(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNDimensional_Hess::~CNDimensional_Hess(void)
  {

  }
//+------------------------------------------------------------------+
//| Empty function body                                              |
//+------------------------------------------------------------------+
void CNDimensional_Hess::Hess(double &x[],double &func,double &grad[],
                              CMatrixDouble &hess,CObject &obj)
  {

  }
//+------------------------------------------------------------------+
//| Calculates vector function f(arg), stores result to fi           |
//+------------------------------------------------------------------+
class CNDimensional_FVec
  {
public:
   //--- constructor, destructor
                     CNDimensional_FVec(void);
                    ~CNDimensional_FVec(void);
   //--- virtual method
   virtual void      FVec(double &x[],double &fi[],CObject &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNDimensional_FVec::CNDimensional_FVec(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNDimensional_FVec::~CNDimensional_FVec(void)
  {

  }
//+------------------------------------------------------------------+
//| Empty function body                                              |
//+------------------------------------------------------------------+
void CNDimensional_FVec::FVec(double &x[],double &fi[],CObject &obj)
  {

  }
//+------------------------------------------------------------------+
//| Calculates f[i] = fi(arg), jac[i,j] = df[i](arg)/d(arg[j])       |
//+------------------------------------------------------------------+
class CNDimensional_Jac
  {
public:
   //--- constructor, destructor
                     CNDimensional_Jac(void);
                    ~CNDimensional_Jac(void);
   //--- virtual method
   virtual void      Jac(double &x[],double &fi[],CMatrixDouble &jac,
                         CObject &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNDimensional_Jac::CNDimensional_Jac(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNDimensional_Jac::~CNDimensional_Jac(void)
  {

  }
//+------------------------------------------------------------------+
//| Empty function body                                              |
//+------------------------------------------------------------------+
void CNDimensional_Jac::Jac(double &x[],double &fi[],CMatrixDouble &jac,
                            CObject &obj)
  {

  }
//+------------------------------------------------------------------+
//| Calculates f(p,q), stores result to func                         |
//+------------------------------------------------------------------+
class CNDimensional_PFunc
  {
public:
   //--- constructor, destructor
                     CNDimensional_PFunc(void);
                    ~CNDimensional_PFunc(void);
   //--- virtual method
   virtual void      PFunc(double &c[],double &x[],double &func,CObject &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNDimensional_PFunc::CNDimensional_PFunc(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNDimensional_PFunc::~CNDimensional_PFunc(void)
  {

  }
//+------------------------------------------------------------------+
//| Empty function body                                              |
//+------------------------------------------------------------------+
void CNDimensional_PFunc::PFunc(double &c[],double &x[],double &func,
                                CObject &obj)
  {

  }
//+------------------------------------------------------------------+
//| Calculates func = f(p,q), grad[i] = df(p,q)/d(p[i])              |
//+------------------------------------------------------------------+
class CNDimensional_PGrad
  {
public:
   //--- constructor, destructor
                     CNDimensional_PGrad(void);
                    ~CNDimensional_PGrad(void);
   //--- virtual method
   virtual void      PGrad(double &c[],double &x[],double &func,double &grad[],CObject &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNDimensional_PGrad::CNDimensional_PGrad(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNDimensional_PGrad::~CNDimensional_PGrad(void)
  {

  }
//+------------------------------------------------------------------+
//| Empty function body                                              |
//+------------------------------------------------------------------+
void CNDimensional_PGrad::PGrad(double &c[],double &x[],double &func,
                                double &grad[],CObject &obj)
  {

  }
//+------------------------------------------------------------------+
//| Calculates func = f(p,q), grad[i] = df(p,q)/d(p[i]),             |
//| hess[i,j] = d2f(p,q)/(d(p[i])*d(p[j]))                           |
//+------------------------------------------------------------------+
class CNDimensional_PHess
  {
public:
   //--- constructor, destructor
                     CNDimensional_PHess(void);
                    ~CNDimensional_PHess(void);
   //--- virtual method
   virtual void      PHess(double &c[],double &x[],double &func,double &grad[],CMatrixDouble &hess,CObject &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNDimensional_PHess::CNDimensional_PHess(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNDimensional_PHess::~CNDimensional_PHess(void)
  {

  }
//+------------------------------------------------------------------+
//| Empty function body                                              |
//+------------------------------------------------------------------+
void CNDimensional_PHess::PHess(double &c[],double &x[],double &func,
                                double &grad[],CMatrixDouble &hess,
                                CObject &obj)
  {

  }
//+------------------------------------------------------------------+
//| Callbacks for ODE solvers: calculates dy/dx for given y[] and x  |
//+------------------------------------------------------------------+
class CNDimensional_ODE_RP
  {
public:
   //--- constructor, destructor
                     CNDimensional_ODE_RP(void);
                    ~CNDimensional_ODE_RP(void);
   //--- virtual method
   virtual void      ODE_RP(double &y[],double x,double &dy[],CObject &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNDimensional_ODE_RP::CNDimensional_ODE_RP(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNDimensional_ODE_RP::~CNDimensional_ODE_RP(void)
  {

  }
//+------------------------------------------------------------------+
//| Empty function body                                              |
//+------------------------------------------------------------------+
void CNDimensional_ODE_RP::ODE_RP(double &y[],double x,double &dy[],
                                  CObject &obj)
  {

  }
//+------------------------------------------------------------------+
//| Callbacks for integrators: calculates f(x) for given x           |
//| (additional parameters xminusa and bminusx contain x-a and b-x)  |
//+------------------------------------------------------------------+
class CIntegrator1_Func
  {
public:
   //--- constructor, destructor
                     CIntegrator1_Func(void);
                    ~CIntegrator1_Func(void);
   //--- virtual method
   virtual void      Int_Func(double x,double xminusa,double bminusx,double &y,CObject &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CIntegrator1_Func::CIntegrator1_Func(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CIntegrator1_Func::~CIntegrator1_Func(void)
  {

  }
//+------------------------------------------------------------------+
//| Empty function body                                              |
//+------------------------------------------------------------------+
void CIntegrator1_Func::Int_Func(double x,double xminusa,double bminusx,
                                 double &y,CObject &obj)
  {

  }
//+------------------------------------------------------------------+
//| Callbacks for progress reports: reports current position of      |
//| optimization algo                                                |
//+------------------------------------------------------------------+
class CNDimensional_Rep
  {
public:
   //--- constructor, destructor
                     CNDimensional_Rep(void);
                    ~CNDimensional_Rep(void);
   //--- virtual method
   virtual void      Rep(double &arg[],double func,CObject &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNDimensional_Rep::CNDimensional_Rep(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNDimensional_Rep::~CNDimensional_Rep(void)
  {

  }
//+------------------------------------------------------------------+
//| Empty function body                                              |
//+------------------------------------------------------------------+
void CNDimensional_Rep::Rep(double &arg[],double func,CObject &obj)
  {

  }
//+------------------------------------------------------------------+
