//+------------------------------------------------------------------+
//|                                                 optimization.mqh |
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
#include "matrix.mqh"
#include "ap.mqh"
#include "alglibinternal.mqh"
#include "linalg.mqh"
//+------------------------------------------------------------------+
//| This object stores state of the nonlinear CG optimizer.          |
//| You should use ALGLIB functions to work with this object.        |
//+------------------------------------------------------------------+
class CMinCGState
  {
public:
   //--- variables
   int               m_n;
   double            m_epsg;
   double            m_epsf;
   double            m_epsx;
   int               m_maxits;
   double            m_stpmax;
   double            m_suggestedstep;
   bool              m_xrep;
   bool              m_drep;
   int               m_cgtype;
   int               m_prectype;
   int               m_vcnt;
   double            m_diffstep;
   int               m_nfev;
   int               m_mcstage;
   int               m_k;
   double            m_fold;
   double            m_stp;
   double            m_curstpmax;
   double            m_laststep;
   double            m_lastscaledstep;
   int               m_mcinfo;
   bool              m_innerresetneeded;
   bool              m_terminationneeded;
   double            m_trimthreshold;
   int               m_rstimer;
   double            m_f;
   bool              m_needf;
   bool              m_needfg;
   bool              m_xupdated;
   bool              m_algpowerup;
   bool              m_lsstart;
   bool              m_lsend;
   RCommState        m_rstate;
   int               m_repiterationscount;
   int               m_repnfev;
   int               m_repterminationtype;
   int               m_debugrestartscount;
   CLinMinState      m_lstate;
   double            m_fbase;
   double            m_fm2;
   double            m_fm1;
   double            m_fp1;
   double            m_fp2;
   double            m_betahs;
   double            m_betady;
   //--- arrays
   double            m_xk[];
   double            m_dk[];
   double            m_xn[];
   double            m_dn[];
   double            m_d[];
   double            m_x[];
   double            m_yk[];
   double            m_s[];
   double            m_g[];
   double            m_diagh[];
   double            m_diaghl2[];
   double            m_work0[];
   double            m_work1[];
   //--- matrix
   CMatrixDouble     m_vcorr;
   //--- constructor, destructor
                     CMinCGState(void);
                    ~CMinCGState(void);
   //--- copy
   void              Copy(CMinCGState &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinCGState::CMinCGState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinCGState::~CMinCGState(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMinCGState::Copy(CMinCGState &obj)
  {
//--- copy variables
   m_n=obj.m_n;
   m_epsg=obj.m_epsg;
   m_epsf=obj.m_epsf;
   m_epsx=obj.m_epsx;
   m_maxits=obj.m_maxits;
   m_stpmax=obj.m_stpmax;
   m_suggestedstep=obj.m_suggestedstep;
   m_xrep=obj.m_xrep;
   m_drep=obj.m_drep;
   m_cgtype=obj.m_cgtype;
   m_prectype=obj.m_prectype;
   m_vcnt=obj.m_vcnt;
   m_diffstep=obj.m_diffstep;
   m_nfev=obj.m_nfev;
   m_mcstage=obj.m_mcstage;
   m_k=obj.m_k;
   m_fold=obj.m_fold;
   m_stp=obj.m_stp;
   m_curstpmax=obj.m_curstpmax;
   m_laststep=obj.m_laststep;
   m_lastscaledstep=obj.m_lastscaledstep;
   m_mcinfo=obj.m_mcinfo;
   m_innerresetneeded=obj.m_innerresetneeded;
   m_terminationneeded=obj.m_terminationneeded;
   m_trimthreshold=obj.m_trimthreshold;
   m_rstimer=obj.m_rstimer;
   m_f=obj.m_f;
   m_needf=obj.m_needf;
   m_needfg=obj.m_needfg;
   m_xupdated=obj.m_xupdated;
   m_algpowerup=obj.m_algpowerup;
   m_lsstart=obj.m_lsstart;
   m_lsend=obj.m_lsend;
   m_repiterationscount=obj.m_repiterationscount;
   m_repnfev=obj.m_repnfev;
   m_repterminationtype=obj.m_repterminationtype;
   m_debugrestartscount=obj.m_debugrestartscount;
   m_fbase=obj.m_fbase;
   m_fm2=obj.m_fm2;
   m_fm1=obj.m_fm1;
   m_fp1=obj.m_fp1;
   m_fp2=obj.m_fp2;
   m_betahs=obj.m_betahs;
   m_betady=obj.m_betady;
   m_rstate.Copy(obj.m_rstate);
   m_lstate.Copy(obj.m_lstate);
//--- copy arrays
   ArrayCopy(m_xk,obj.m_xk);
   ArrayCopy(m_dk,obj.m_dk);
   ArrayCopy(m_xn,obj.m_xn);
   ArrayCopy(m_dn,obj.m_dn);
   ArrayCopy(m_d,obj.m_d);
   ArrayCopy(m_x,obj.m_x);
   ArrayCopy(m_yk,obj.m_yk);
   ArrayCopy(m_s,obj.m_s);
   ArrayCopy(m_g,obj.m_g);
   ArrayCopy(m_diagh,obj.m_diagh);
   ArrayCopy(m_diaghl2,obj.m_diaghl2);
   ArrayCopy(m_work0,obj.m_work0);
   ArrayCopy(m_work1,obj.m_work1);
//--- matrix
   m_vcorr=obj.m_vcorr;
  }
//+------------------------------------------------------------------+
//| This object stores state of the nonlinear CG optimizer.          |
//| You should use ALGLIB functions to work with this object.        |
//+------------------------------------------------------------------+
class CMinCGStateShell
  {
private:
   CMinCGState       m_innerobj;

public:
   //--- constructors, destructor
                     CMinCGStateShell(void);
                     CMinCGStateShell(CMinCGState &obj);
                    ~CMinCGStateShell(void);
   //--- methods
   bool              GetNeedF(void);
   void              SetNeedF(const bool b);
   bool              GetNeedFG(void);
   void              SetNeedFG(const bool b);
   bool              GetXUpdated(void);
   void              SetXUpdated(const bool b);
   double            GetF(void);
   void              SetF(const double d);
   CMinCGState      *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinCGStateShell::CMinCGStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMinCGStateShell::CMinCGStateShell(CMinCGState &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinCGStateShell::~CMinCGStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needf                          |
//+------------------------------------------------------------------+
bool CMinCGStateShell::GetNeedF(void)
  {
//--- return result
   return(m_innerobj.m_needf);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needf                         |
//+------------------------------------------------------------------+
void CMinCGStateShell::SetNeedF(const bool b)
  {
//--- change value
   m_innerobj.m_needf=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needfg                         |
//+------------------------------------------------------------------+
bool CMinCGStateShell::GetNeedFG(void)
  {
//--- return result
   return(m_innerobj.m_needfg);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needfg                        |
//+------------------------------------------------------------------+
void CMinCGStateShell::SetNeedFG(const bool b)
  {
//--- change value
   m_innerobj.m_needfg=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable xupdated                       |
//+------------------------------------------------------------------+
bool CMinCGStateShell::GetXUpdated(void)
  {
//--- return result
   return(m_innerobj.m_xupdated);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable xupdated                      |
//+------------------------------------------------------------------+
void CMinCGStateShell::SetXUpdated(const bool b)
  {
//--- change value
   m_innerobj.m_xupdated=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable f                              |
//+------------------------------------------------------------------+
double CMinCGStateShell::GetF(void)
  {
//--- return result
   return(m_innerobj.m_f);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable f                             |
//+------------------------------------------------------------------+
void CMinCGStateShell::SetF(const double d)
  {
//--- change value
   m_innerobj.m_f=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMinCGState *CMinCGStateShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CMinCG                                       |
//+------------------------------------------------------------------+
class CMinCGReport
  {
public:
   int               m_iterationscount;
   int               m_nfev;
   int               m_terminationtype;
   //--- constructor, destructor
                     CMinCGReport(void);
                    ~CMinCGReport(void);
   //--- copy
   void              Copy(CMinCGReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinCGReport::CMinCGReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinCGReport::~CMinCGReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMinCGReport::Copy(CMinCGReport &obj)
  {
//--- copy variables
   m_iterationscount=obj.m_iterationscount;
   m_nfev=obj.m_nfev;
   m_terminationtype=obj.m_terminationtype;
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CMinCGReport                     |
//+------------------------------------------------------------------+
class CMinCGReportShell
  {
private:
   CMinCGReport      m_innerobj;
public:
   //--- constructors, destructor
                     CMinCGReportShell(void);
                     CMinCGReportShell(CMinCGReport &obj);
                    ~CMinCGReportShell(void);
   //--- methods
   int               GetIterationsCount(void);
   void              SetIterationsCount(const int i);
   int               GetNFev(void);
   void              SetNFev(const int i);
   int               GetTerminationType(void);
   void              SetTerminationType(const int i);
   CMinCGReport     *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinCGReportShell::CMinCGReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMinCGReportShell::CMinCGReportShell(CMinCGReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinCGReportShell::~CMinCGReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable iterationscount                |
//+------------------------------------------------------------------+
int CMinCGReportShell::GetIterationsCount(void)
  {
//--- return result
   return(m_innerobj.m_iterationscount);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable iterationscount               |
//+------------------------------------------------------------------+
void CMinCGReportShell::SetIterationsCount(const int i)
  {
//--- change value
   m_innerobj.m_iterationscount=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable nfev                           |
//+------------------------------------------------------------------+
int CMinCGReportShell::GetNFev(void)
  {
//--- return result
   return(m_innerobj.m_nfev);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable nfev                          |
//+------------------------------------------------------------------+
void CMinCGReportShell::SetNFev(const int i)
  {
//--- change value
   m_innerobj.m_nfev=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable terminationtype                |
//+------------------------------------------------------------------+
int CMinCGReportShell::GetTerminationType(void)
  {
//--- return result
   return(m_innerobj.m_terminationtype);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable terminationtype               |
//+------------------------------------------------------------------+
void CMinCGReportShell::SetTerminationType(const int i)
  {
//--- change value
   m_innerobj.m_terminationtype=i;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMinCGReport *CMinCGReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Conjugate gradient optimizer                                     |
//+------------------------------------------------------------------+
class CMinCG
  {
private:
   //--- private methods
   static void       ClearRequestFields(CMinCGState &state);
   static void       PreconditionedMultiply(CMinCGState &state,double &x[],double &work0[],double &work1[]);
   static double     PreconditionedMultiply2(CMinCGState &state,double &x[],double &y[],double &work0[],double &work1[]);
   static void       MinCGInitInternal(const int n,const double diffstep,CMinCGState &state);
   //--- auxiliary functions for MinCGIteration
   static void       Func_lbl_rcomm(CMinCGState &state,int n,int i,double betak,double v,double vv);
   static bool       Func_lbl_18(CMinCGState &state,int &n,int &i,double &betak,double &v,double &vv);
   static bool       Func_lbl_19(CMinCGState &state,int &n,int &i,double &betak,double &v,double &vv);
   static bool       Func_lbl_22(CMinCGState &state,int &n,int &i,double &betak,double &v,double &vv);
   static bool       Func_lbl_24(CMinCGState &state,int &n,int &i,double &betak,double &v,double &vv);
   static bool       Func_lbl_26(CMinCGState &state,int &n,int &i,double &betak,double &v,double &vv);
   static bool       Func_lbl_28(CMinCGState &state,int &n,int &i,double &betak,double &v,double &vv);
   static bool       Func_lbl_30(CMinCGState &state,int &n,int &i,double &betak,double &v,double &vv);
   static bool       Func_lbl_31(CMinCGState &state,int &n,int &i,double &betak,double &v,double &vv);
   static bool       Func_lbl_33(CMinCGState &state,int &n,int &i,double &betak,double &v,double &vv);
   static bool       Func_lbl_34(CMinCGState &state,int &n,int &i,double &betak,double &v,double &vv);
   static bool       Func_lbl_37(CMinCGState &state,int &n,int &i,double &betak,double &v,double &vv);
   static bool       Func_lbl_39(CMinCGState &state,int &n,int &i,double &betak,double &v,double &vv);
public:
   //--- class constants
   static const int  m_rscountdownlen;
   static const double m_gtol;
   //--- constructor, destructor
                     CMinCG(void);
                    ~CMinCG(void);
   //--- public methods
   static void       MinCGCreate(const int n,double &x[],CMinCGState &state);
   static void       MinCGCreateF(const int n,double &x[],const double diffstep,CMinCGState &state);
   static void       MinCGSetCond(CMinCGState &state,const double epsg,const double epsf,double epsx,const int maxits);
   static void       MinCGSetScale(CMinCGState &state,double &s[]);
   static void       MinCGSetXRep(CMinCGState &state,const bool needxrep);
   static void       MinCGSetDRep(CMinCGState &state,const bool needdrep);
   static void       MinCGSetCGType(CMinCGState &state,int cgtype);
   static void       MinCGSetStpMax(CMinCGState &state,const double stpmax);
   static void       MinCGSuggestStep(CMinCGState &state,const double stp);
   static void       MinCGSetPrecDefault(CMinCGState &state);
   static void       MinCGSetPrecDiag(CMinCGState &state,double &d[]);
   static void       MinCGSetPrecScale(CMinCGState &state);
   static void       MinCGResults(CMinCGState &state,double &x[],CMinCGReport &rep);
   static void       MinCGResultsBuf(CMinCGState &state,double &x[],CMinCGReport &rep);
   static void       MinCGRestartFrom(CMinCGState &state,double &x[]);
   static void       MinCGSetPrecDiagFast(CMinCGState &state,double &d[]);
   static void       MinCGSetPrecLowRankFast(CMinCGState &state,double &d1[],double &c[],CMatrixDouble &v,const int vcnt);
   static void       MinCGSetPrecVarPart(CMinCGState &state,double &d2[]);
   static bool       MinCGIteration(CMinCGState &state);
  };
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const int    CMinCG::m_rscountdownlen=10;
const double CMinCG::m_gtol=0.3;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinCG::CMinCG(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinCG::~CMinCG(void)
  {

  }
//+------------------------------------------------------------------+
//|         NONLINEAR CONJUGATE GRADIENT METHOD                      |
//| DESCRIPTION:                                                     |
//| The subroutine minimizes function F(x) of N arguments by using   |
//| one of the nonlinear conjugate gradient methods.                 |
//| These CG methods are globally convergent (even on non-convex     |
//| functions) as long as grad(f) is Lipschitz continuous in a some  |
//| neighborhood of the L = { x : f(x)<=f(x0) }.                     |
//| REQUIREMENTS:                                                    |
//| Algorithm will request following information during its          |
//| operation:                                                       |
//| * function value F and its gradient G (simultaneously) at given  |
//|   point X                                                        |
//| USAGE:                                                           |
//| 1. User initializes algorithm state with MinCGCreate() call      |
//| 2. User tunes solver parameters with MinCGSetCond(),             |
//|    MinCGSetStpMax() and other functions                          |
//| 3. User calls MinCGOptimize() function which takes algorithm     |
//|    state and pointer (delegate, etc.) to callback function which |
//|    calculates F/G.                                               |
//| 4. User calls MinCGResults() to get solution                     |
//| 5. Optionally, user may call MinCGRestartFrom() to solve another |
//|    problem with same N but another starting point and/or another |
//|    function. MinCGRestartFrom() allows to reuse already          |
//|    initialized structure.                                        |
//| INPUT PARAMETERS:                                                |
//|     N       -   problem dimension, N>0:                          |
//|                 * if given, only leading N elements of X are used|
//|                 * if not given, automatically determined from    |
//|                   size of X                                      |
//|     X       -   starting point, array[0..N-1].                   |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//+------------------------------------------------------------------+
static void CMinCG::MinCGCreate(const int n,double &x[],CMinCGState &state)
  {
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N too small!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- function call
   MinCGInitInternal(n,0.0,state);
//--- function call
   MinCGRestartFrom(state,x);
  }
//+------------------------------------------------------------------+
//| The subroutine is finite difference variant of MinCGCreate().    |
//| It uses finite differences in order to differentiate target      |
//| function.                                                        |
//| Description below contains information which is specific to this |
//| function only. We recommend to read comments on MinCGCreate() in |
//| order to get more information about creation of CG optimizer.    |
//| INPUT PARAMETERS:                                                |
//|     N       -   problem dimension, N>0:                          |
//|                 * if given, only leading N elements of X are     |
//|                   used                                           |
//|                 * if not given, automatically determined from    |
//|                   size of X                                      |
//|     X       -   starting point, array[0..N-1].                   |
//|     DiffStep-   differentiation step, >0                         |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//| NOTES:                                                           |
//| 1. algorithm uses 4-point central formula for differentiation.   |
//| 2. differentiation step along I-th axis is equal to              |
//|    DiffStep*S[I] where S[] is scaling vector which can be set by |
//|    MinCGSetScale() call.                                         |
//| 3. we recommend you to use moderate values of differentiation    |
//|    step. Too large step will result in too large truncation      |
//|    errors, while too small step will result in too large         |
//|    numerical errors. 1.0E-6 can be good value to start with.     |
//| 4. Numerical differentiation is very inefficient - one gradient  |
//|    calculation needs 4*N function evaluations. This function will|
//|    work for any N - either small (1...10), moderate (10...100) or|
//|    large  (100...). However, performance penalty will be too     |
//|    severe for any N's except for small ones.                     |
//|    We should also say that code which relies on numerical        |
//|    differentiation is less robust and precise. L-BFGS  needs     |
//|    exact gradient values. Imprecise gradient may slow down       |
//|    convergence, especially on highly nonlinear problems.         |
//|    Thus we recommend to use this function for fast prototyping   |
//|    on small- dimensional problems only, and to implement         |
//|    analytical gradient as soon as possible.                      |
//+------------------------------------------------------------------+
static void CMinCG::MinCGCreateF(const int n,double &x[],const double diffstep,
                                 CMinCGState &state)
  {
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N too small!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(diffstep),__FUNCTION__+": DiffStep is infinite or NaN!"))
      return;
//--- check
   if(!CAp::Assert(diffstep>0.0,__FUNCTION__+": DiffStep is non-positive!"))
      return;
//--- function call
   MinCGInitInternal(n,diffstep,state);
//--- function call
   MinCGRestartFrom(state,x);
  }
//+------------------------------------------------------------------+
//| This function sets stopping conditions for CG optimization       |
//| algorithm.                                                       |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     EpsG    -   >=0                                              |
//|                 The subroutine finishes its work if the condition|
//|                 |v|<EpsG is satisfied, where:                    |
//|                 * |.| means Euclidian norm                       |
//|                 * v - scaled gradient vector, v[i]=g[i]*s[i]     |
//|                 * g - gradient                                   |
//|                 * s - scaling coefficients set by MinCGSetScale()|
//|     EpsF    -   >=0                                              |
//|                 The subroutine finishes its work if on k+1-th    |
//|                 iteration the condition |F(k+1)-F(k)| <=         |
//|                 <= EpsF*max{|F(k)|,|F(k+1)|,1} is satisfied.     |
//|     EpsX    -   >=0                                              |
//|                 The subroutine finishes its work if on k+1-th    |
//|                 iteration the condition |v|<=EpsX is fulfilled,  |
//|                 where:                                           |
//|                 * |.| means Euclidian norm                       |
//|                 * v - scaled step vector, v[i]=dx[i]/s[i]        |
//|                 * dx - ste pvector, dx=X(k+1)-X(k)               |
//|                 * s - scaling coefficients set by MinCGSetScale()|
//|     MaxIts  -   maximum number of iterations. If MaxIts=0, the   |
//|                 number of iterations is unlimited.               |
//| Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will|
//| lead to automatic stopping criterion selection (small EpsX).     |
//+------------------------------------------------------------------+
static void CMinCG::MinCGSetCond(CMinCGState &state,const double epsg,
                                 const double epsf,double epsx,const int maxits)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsg),__FUNCTION__+": EpsG is not finite number!"))
      return;
//--- check
   if(!CAp::Assert(epsg>=0.0,__FUNCTION__+": negative EpsG!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsf),__FUNCTION__+": EpsF is not finite number!"))
      return;
//--- check
   if(!CAp::Assert(epsf>=0.0,__FUNCTION__+": negative EpsF!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsx),__FUNCTION__+": EpsX is not finite number!"))
      return;
//--- check
   if(!CAp::Assert(epsx>=0.0,__FUNCTION__+": negative EpsX!"))
      return;
//--- check
   if(!CAp::Assert(maxits>=0,__FUNCTION__+": negative MaxIts!"))
      return;
//--- check
   if(epsg==0.0 && epsf==0.0 && epsx==0.0 && maxits==0)
      epsx=1.0E-6;
//--- change values
   state.m_epsg=epsg;
   state.m_epsf=epsf;
   state.m_epsx=epsx;
   state.m_maxits=maxits;
  }
//+------------------------------------------------------------------+
//| This function sets scaling coefficients for CG optimizer.        |
//| ALGLIB optimizers use scaling matrices to test stopping          |
//| conditions (step size and gradient are scaled before comparison  |
//| with tolerances). Scale of the I-th variable is a translation    |
//| invariant measure of:                                            |
//| a) "how large" the variable is                                   |
//| b) how large the step should be to make significant changes in   |
//|    the function                                                  |
//| Scaling is also used by finite difference variant of CG          |
//| optimizer - step along I-th axis is equal to DiffStep*S[I].      |
//| In most optimizers (and in the CG too) scaling is NOT a form of  |
//| preconditioning. It just affects stopping conditions. You should |
//| set preconditioner by separate call to one of the                |
//| MinCGSetPrec...() functions.                                     |
//| There is special preconditioning mode, however, which uses       |
//| scaling coefficients to form diagonal preconditioning matrix.    |
//| You can turn this mode on, if you want. But you should understand|
//| that scaling is not the same thing as preconditioning - these are|
//| two different, although related forms of tuning solver.          |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure stores algorithm state                 |
//|     S       -   array[N], non-zero scaling coefficients          |
//|                 S[i] may be negative, sign doesn't matter.       |
//+------------------------------------------------------------------+
static void CMinCG::MinCGSetScale(CMinCGState &state,double &s[])
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(CAp::Len(s)>=state.m_n,__FUNCTION__+": Length(S)<N"))
      return;
   for(i=0;i<=state.m_n-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(s[i]),__FUNCTION__+": S contains infinite or NAN elements"))
         return;
      //--- check
      if(!CAp::Assert(s[i]!=0.0,__FUNCTION__+": S contains zero elements"))
         return;
      state.m_s[i]=MathAbs(s[i]);
     }
  }
//+------------------------------------------------------------------+
//| This function turns on/off reporting.                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     NeedXRep-   whether iteration reports are needed or not      |
//| If NeedXRep is True, algorithm will call rep() callback function |
//| if it is provided to MinCGOptimize().                            |
//+------------------------------------------------------------------+
static void CMinCG::MinCGSetXRep(CMinCGState &state,const bool needxrep)
  {
//--- change value
   state.m_xrep=needxrep;
  }
//+------------------------------------------------------------------+
//| This function turns on/off line search reports.                  |
//| These reports are described in more details in developer-only    | 
//| comments on MinCGState object.                                   |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     NeedDRep-   whether line search reports are needed or not    |
//| This function is intended for private use only. Turning it on    |
//| artificially may cause program failure.                          |
//+------------------------------------------------------------------+
static void CMinCG::MinCGSetDRep(CMinCGState &state,const bool needdrep)
  {
//--- change value
   state.m_drep=needdrep;
  }
//+------------------------------------------------------------------+
//| This function sets CG algorithm.                                 |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     CGType  -   algorithm type:                                  |
//|                 * -1    automatic selection of the best          |
//|                         algorithm                                |
//|                 * 0     DY (Dai and Yuan) algorithm              |
//|                 * 1     Hybrid DY-HS algorithm                   |
//+------------------------------------------------------------------+
static void CMinCG::MinCGSetCGType(CMinCGState &state,int cgtype)
  {
//--- check
   if(!CAp::Assert(cgtype>=-1 && cgtype<=1,__FUNCTION__+": incorrect CGType!"))
      return;
//--- check
   if(cgtype==-1)
      cgtype=1;
//--- change value
   state.m_cgtype=cgtype;
  }
//+------------------------------------------------------------------+
//| This function sets maximum step length                           |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     StpMax  -   maximum step length, >=0. Set StpMax to 0.0, if  |
//|                 you don't want to limit step length.             |
//| Use this subroutine when you optimize target function which      |
//| contains exp() or other fast growing functions, and optimization |
//| algorithm makes too large steps which leads to overflow. This    |
//| function allows us to reject steps that are too large (and       |
//| therefore expose us to the possible overflow) without actually   |
//| calculating function value at the x+stp*d.                       |
//+------------------------------------------------------------------+
static void CMinCG::MinCGSetStpMax(CMinCGState &state,const double stpmax)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(stpmax),__FUNCTION__+": StpMax is not finite!"))
      return;
//--- check
   if(!CAp::Assert(stpmax>=0.0,__FUNCTION__+": StpMax<0!"))
      return;
//--- change value
   state.m_stpmax=stpmax;
  }
//+------------------------------------------------------------------+
//| This function allows to suggest initial step length to the CG    |
//| algorithm.                                                       |
//| Suggested step length is used as starting point for the line     |
//| search. It can be useful when you have badly scaled problem, i.e.|
//| when ||grad|| (which is used as initial estimate for the first   |
//| step) is many orders of magnitude different from the desired     |
//| step.                                                            |
//| Line search may fail on such problems without good estimate of   |
//| initial step length. Imagine, for example, problem with          |
//| ||grad||=10^50 and desired step equal to 0.1 Line search         |
//| function will use 10^50 as initial step, then it will decrease   |
//| step length by 2 (up to 20 attempts) and will get 10^44, which is|
//| still too large.                                                 |
//| This function allows us to tell than line search should be       |
//| started from some moderate step length, like 1.0, so algorithm   |
//| will be able to detect desired step length in a several searches.|
//| Default behavior (when no step is suggested) is to use           |
//| preconditioner, if it is available, to generate initial estimate |
//| of step length.                                                  |
//| This function influences only first iteration of algorithm. It   |
//| should be called between MinCGCreate/MinCGRestartFrom() call and |
//| MinCGOptimize call. Suggested step is ignored if you have        |
//| preconditioner.                                                  |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure used to store algorithm state.         |
//|     Stp     -   initial estimate of the step length.             |
//|                 Can be zero (no estimate).                       |
//+------------------------------------------------------------------+
static void CMinCG::MinCGSuggestStep(CMinCGState &state,const double stp)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(stp),__FUNCTION__+": Stp is infinite or NAN"))
      return;
//--- check
   if(!CAp::Assert(stp>=0.0,__FUNCTION__+": Stp<0"))
      return;
//--- change value
   state.m_suggestedstep=stp;
  }
//+------------------------------------------------------------------+
//| Modification of the preconditioner: preconditioning is turned    |
//| off.                                                             |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//| NOTE: you can change preconditioner "on the fly", during         |
//| algorithm iterations.                                            |
//+------------------------------------------------------------------+
static void CMinCG::MinCGSetPrecDefault(CMinCGState &state)
  {
//--- change values
   state.m_prectype=0;
   state.m_innerresetneeded=true;
  }
//+------------------------------------------------------------------+
//| Modification of the preconditioner: diagonal of approximate      |
//| Hessian is used.                                                 |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     D       -   diagonal of the approximate Hessian,             |
//|                 array[0..N-1], (if larger, only leading N        |
//|                 elements are used).                              |
//| NOTE: you can change preconditioner "on the fly", during         |
//| algorithm iterations.                                            |
//| NOTE 2: D[i] should be positive. Exception will be thrown        |
//| otherwise.                                                       |
//| NOTE 3: you should pass diagonal of approximate Hessian - NOT    |
//| ITS INVERSE.                                                     |
//+------------------------------------------------------------------+
static void CMinCG::MinCGSetPrecDiag(CMinCGState &state,double &d[])
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(CAp::Len(d)>=state.m_n,__FUNCTION__+": D is too short"))
      return;
   for(i=0;i<=state.m_n-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(d[i]),__FUNCTION__+": D contains infinite or NAN elements"))
         return;
      //--- check
      if(!CAp::Assert((double)(d[i])>0.0,__FUNCTION__+": D contains non-positive elements"))
         return;
     }
//--- function call
   MinCGSetPrecDiagFast(state,d);
  }
//+------------------------------------------------------------------+
//| Modification of the preconditioner: scale-based diagonal         |
//| preconditioning.                                                 |
//| This preconditioning mode can be useful when you don't have      |
//| approximate diagonal of Hessian, but you know that your variables|
//| are badly scaled (for example, one variable is in [1,10], and    |
//| another in [1000,100000]), and most part of the ill-conditioning |
//| comes from different scales of vars.                             |
//| In this case simple scale-based  preconditioner,                 |
//| with H[i] = 1/(s[i]^2), can greatly improve convergence.         |
//| IMPRTANT: you should set scale of your variables with            |
//| MinCGSetScale() call (before or after MinCGSetPrecScale() call). |
//| Without knowledge of the scale of your variables scale-based     |
//| preconditioner will be just unit matrix.                         |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//| NOTE: you can change preconditioner "on the fly", during         |
//| algorithm iterations.                                            |
//+------------------------------------------------------------------+
static void CMinCG::MinCGSetPrecScale(CMinCGState &state)
  {
//--- change values
   state.m_prectype=3;
   state.m_innerresetneeded=true;
  }
//+------------------------------------------------------------------+
//| Conjugate gradient results                                       |
//| INPUT PARAMETERS:                                                |
//|     State   -   algorithm state                                  |
//| OUTPUT PARAMETERS:                                               |
//|     X       -   array[0..N-1], solution                          |
//|     Rep     -   optimization report:                             |
//|                 * Rep.TerminationType completetion code:         |
//|                     *  1    relative function improvement is no  |
//|                             more than EpsF.                      |
//|                     *  2    relative step is no more than EpsX.  |
//|                     *  4    gradient norm is no more than EpsG   |
//|                     *  5    MaxIts steps was taken               |
//|                     *  7    stopping conditions are too          |
//|                             stringent, further improvement is    |
//|                             impossible, we return best X found   |
//|                             so far                               |
//|                     *  8    terminated by user                   |
//|                 * Rep.IterationsCount contains iterations count  |
//|                 * NFEV countains number of function calculations |
//+------------------------------------------------------------------+
static void CMinCG::MinCGResults(CMinCGState &state,double &x[],CMinCGReport &rep)
  {
//--- reset memory
   ArrayResizeAL(x,0);
//--- function call
   MinCGResultsBuf(state,x,rep);
  }
//+------------------------------------------------------------------+
//| Conjugate gradient results                                       |
//| Buffered implementation of MinCGResults(), which uses            |
//| pre-allocated buffer to store X[]. If buffer size is too small,  |
//| it resizes buffer.It is intended to be used in the inner cycles  |
//| of performance critical algorithms where array reallocation      |
//| penalty is too large to be ignored.                              |
//+------------------------------------------------------------------+
static void CMinCG::MinCGResultsBuf(CMinCGState &state,double &x[],CMinCGReport &rep)
  {
//--- create a variable
   int i_=0;
//--- check
   if(CAp::Len(x)<state.m_n)
      ArrayResizeAL(x,state.m_n);
//--- copy
   for(i_=0;i_<=state.m_n-1;i_++)
      x[i_]=state.m_xn[i_];
//--- change values
   rep.m_iterationscount=state.m_repiterationscount;
   rep.m_nfev=state.m_repnfev;
   rep.m_terminationtype=state.m_repterminationtype;
  }
//+------------------------------------------------------------------+
//| This  subroutine  restarts  CG  algorithm from new point. All    |
//| optimization parameters are left unchanged.                      |
//| This function allows to solve multiple optimization problems     |
//| (which must have same number of dimensions) without object       |
//| reallocation penalty.                                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure used to store algorithm state.         |
//|     X       -   new starting point.                              |
//+------------------------------------------------------------------+
static void CMinCG::MinCGRestartFrom(CMinCGState &state,double &x[])
  {
//--- create a variable
   int i_=0;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=state.m_n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,state.m_n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- copy
   for(i_=0;i_<=state.m_n-1;i_++)
      state.m_x[i_]=x[i_];
//--- function call
   MinCGSuggestStep(state,0.0);
//--- allocation
   ArrayResizeAL(state.m_rstate.ia,2);
   ArrayResizeAL(state.m_rstate.ra,3);
//--- change value
   state.m_rstate.stage=-1;
//--- function call
   ClearRequestFields(state);
  }
//+------------------------------------------------------------------+
//| Faster version of MinCGSetPrecDiag(), for time-critical parts of |
//| code, without safety checks.                                     |
//+------------------------------------------------------------------+
static void CMinCG::MinCGSetPrecDiagFast(CMinCGState &state,double &d[])
  {
//--- create a variable
   int i=0;
//--- function call
   CApServ::RVectorSetLengthAtLeast(state.m_diagh,state.m_n);
//--- function call
   CApServ::RVectorSetLengthAtLeast(state.m_diaghl2,state.m_n);
//--- change values
   state.m_prectype=2;
   state.m_vcnt=0;
   state.m_innerresetneeded=true;
//--- copy
   for(i=0;i<=state.m_n-1;i++)
     {
      state.m_diagh[i]=d[i];
      state.m_diaghl2[i]=0.0;
     }
  }
//+------------------------------------------------------------------+
//| This function sets low-rank preconditioner for Hessian matrix    |
//| H=D+V'*C*V, where:                                               |
//| * H is a Hessian matrix, which is approximated by D/V/C          |
//| * D=D1+D2 is a diagonal matrix, which includes two positive      |
//|   definite terms:                                                |
//|   * constant term D1 (is not updated or infrequently updated)    |
//|   * variable term D2 (can be cheaply updated from iteration to   |
//|     iteration)                                                   |
//| * V is a low-rank correction                                     |
//| * C is a diagonal factor of low-rank correction                  |
//| Preconditioner P is calculated using approximate Woodburry       |
//| formula:                                                         |
//|     P  = D^(-1) - D^(-1)*V'*(C^(-1)+V*D1^(-1)*V')^(-1)*V*D^(-1)  |
//|        = D^(-1) - D^(-1)*VC'*VC*D^(-1),                          |
//| where                                                            |
//|     VC = sqrt(B)*V                                               |
//|     B  = (C^(-1)+V*D1^(-1)*V')^(-1)                              |
//| Note that B is calculated using constant term (D1) only, which   |
//| allows us to update D2 without recalculation of B or VC. Such    |
//| preconditioner is exact when D2 is zero. When D2 is non-zero, it |
//| is only approximation, but very good and cheap one.              |
//| This function accepts D1, V, C.                                  |
//| D2 is set to zero by default.                                    |
//| Cost of this update is O(N*VCnt*VCnt), but D2 can be updated in  |
//| just O(N) by MinCGSetPrecVarPart.                                |
//+------------------------------------------------------------------+
static void CMinCG::MinCGSetPrecLowRankFast(CMinCGState &state,double &d1[],
                                            double &c[],CMatrixDouble &v,
                                            const int vcnt)
  {
//--- create variables
   int    i=0;
   int    i_=0;
   int    j=0;
   int    k=0;
   int    n=0;
   double t=0;
//--- create matrix
   CMatrixDouble b;
//--- check
   if(vcnt==0)
     {
      //--- function call
      MinCGSetPrecDiagFast(state,d1);
      return;
     }
//--- initialization
   n=state.m_n;
   b.Resize(vcnt,vcnt);
//--- function call
   CApServ::RVectorSetLengthAtLeast(state.m_diagh,n);
//--- function call
   CApServ::RVectorSetLengthAtLeast(state.m_diaghl2,n);
//--- function call
   CApServ::RMatrixSetLengthAtLeast(state.m_vcorr,vcnt,n);
   state.m_prectype=2;
   state.m_vcnt=vcnt;
   state.m_innerresetneeded=true;
//--- copy
   for(i=0;i<=n-1;i++)
     {
      state.m_diagh[i]=d1[i];
      state.m_diaghl2[i]=0.0;
     }
//--- calculation
   for(i=0;i<=vcnt-1;i++)
     {
      for(j=i;j<=vcnt-1;j++)
        {
         t=0;
         for(k=0;k<=n-1;k++)
            t=t+v[i][k]*v[j][k]/d1[k];
         b[i].Set(j,t);
        }
      b[i].Set(i,b[i][i]+1.0/c[i]);
     }
//--- check
   if(!CTrFac::SPDMatrixCholeskyRec(b,0,vcnt,true,state.m_work0))
     {
      state.m_vcnt=0;
      return;
     }
//--- calculation
   for(i=0;i<=vcnt-1;i++)
     {
      for(i_=0;i_<=n-1;i_++)
         state.m_vcorr[i].Set(i_,v[i][i_]);
      //--- change values
      for(j=0;j<=i-1;j++)
        {
         t=b[j][i];
         for(i_=0;i_<=n-1;i_++)
            state.m_vcorr[i].Set(i_,state.m_vcorr[i][i_]-t*state.m_vcorr[j][i_]);
        }
      t=1/b[i][i];
      //--- change values
      for(i_=0;i_<=n-1;i_++)
         state.m_vcorr[i].Set(i_,t*state.m_vcorr[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| This function updates variable part (diagonal matrix D2)         |
//| of low-rank preconditioner.                                      |
//| This update is very cheap and takes just O(N) time.              |
//| It has no effect with default preconditioner.                    |
//+------------------------------------------------------------------+
static void CMinCG::MinCGSetPrecVarPart(CMinCGState &state,double &d2[])
  {
//--- create variables
   int i=0;
   int n=0;
//--- initialization
   n=state.m_n;
//--- copy
   for(i=0;i<=n-1;i++)
      state.m_diaghl2[i]=d2[i];
  }
//+------------------------------------------------------------------+
//| Clears request fileds (to be sure that we don't forgot to clear  |
//| something)                                                       |
//+------------------------------------------------------------------+
static void CMinCG::ClearRequestFields(CMinCGState &state)
  {
//--- change values
   state.m_needf=false;
   state.m_needfg=false;
   state.m_xupdated=false;
   state.m_lsstart=false;
   state.m_lsend=false;
   state.m_algpowerup=false;
  }
//+------------------------------------------------------------------+
//| This function calculates preconditioned product H^(-1)*x and     |
//| stores result back into X. Work0[] and Work1[] are used as       |
//| temporaries (size must be at least N; this function doesn't      |
//| allocate arrays).                                                |
//+------------------------------------------------------------------+
static void CMinCG::PreconditionedMultiply(CMinCGState &state,double &x[],
                                           double &work0[],double &work1[])
  {
//--- create variables
   int    i=0;
   int    n=0;
   int    vcnt=0;
   double v=0;
   int    i_=0;
//--- initialization
   n=state.m_n;
   vcnt=state.m_vcnt;
//--- check
   if(state.m_prectype==0)
      return;
//--- check
   if(state.m_prectype==3)
     {
      for(i=0;i<=n-1;i++)
         x[i]=x[i]*state.m_s[i]*state.m_s[i];
      //--- exit the function
      return;
     }
//--- check
   if(!CAp::Assert(state.m_prectype==2,__FUNCTION__+": internal error (unexpected PrecType)"))
      return;
//--- handle part common for VCnt=0 and VCnt<>0
   for(i=0;i<=n-1;i++)
      x[i]=x[i]/(state.m_diagh[i]+state.m_diaghl2[i]);
//--- if VCnt>0
   if(vcnt>0)
     {
      //--- calculation work0
      for(i=0;i<=vcnt-1;i++)
        {
         v=0.0;
         for(i_=0;i_<=n-1;i_++)
            v+=state.m_vcorr[i][i_]*x[i_];
         work0[i]=v;
        }
      //--- calculation work1
      for(i=0;i<=n-1;i++)
         work1[i]=0;
      for(i=0;i<=vcnt-1;i++)
        {
         v=work0[i];
         for(i_=0;i_<=n-1;i_++)
            state.m_work1[i_]=state.m_work1[i_]+v*state.m_vcorr[i][i_];
        }
      //--- change x
      for(i=0;i<=n-1;i++)
         x[i]=x[i]-state.m_work1[i]/(state.m_diagh[i]+state.m_diaghl2[i]);
     }
  }
//+------------------------------------------------------------------+
//| This function calculates preconditioned product x'*H^(-1)*y.     |
//| Work0[] and Work1[] are used as temporaries (size must be at     |
//| least N; this function doesn't allocate arrays).                 |
//+------------------------------------------------------------------+
static double CMinCG::PreconditionedMultiply2(CMinCGState &state,double &x[],
                                              double &y[],double &work0[],
                                              double &work1[])
  {
//--- create variables
   double result=0;
   int    i=0;
   int    n=0;
   int    vcnt=0;
   double v0=0;
   double v1=0;
   int    i_=0;
//--- initialization
   n=state.m_n;
   vcnt=state.m_vcnt;
//--- no preconditioning
   if(state.m_prectype==0)
     {
      v0=0.0;
      for(i_=0;i_<=n-1;i_++)
         v0+=x[i_]*y[i_];
      //--- return result
      return(v0);
     }
//--- check
   if(state.m_prectype==3)
     {
      result=0;
      for(i=0;i<=n-1;i++)
         result=result+x[i]*state.m_s[i]*state.m_s[i]*y[i];
      //--- return result
      return(result);
     }
//--- check
   if(!CAp::Assert(state.m_prectype==2,__FUNCTION__+": internal error (unexpected PrecType)"))
      return(EMPTY_VALUE);
//--- low rank preconditioning
   result=0.0;
   for(i=0;i<=n-1;i++)
      result=result+x[i]*y[i]/(state.m_diagh[i]+state.m_diaghl2[i]);
//--- check
   if(vcnt>0)
     {
      //--- prepare arrays
      for(i=0;i<=n-1;i++)
        {
         work0[i]=x[i]/(state.m_diagh[i]+state.m_diaghl2[i]);
         work1[i]=y[i]/(state.m_diagh[i]+state.m_diaghl2[i]);
        }
      for(i=0;i<=vcnt-1;i++)
        {
         //--- calculation
         v0=0.0;
         for(i_=0;i_<=n-1;i_++)
            v0+=work0[i_]*state.m_vcorr[i][i_];
         v1=0.0;
         for(i_=0;i_<=n-1;i_++)
            v1+=work1[i_]*state.m_vcorr[i][i_];
         //--- get result
         result=result-v0*v1;
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Internal initialization subroutine                               |
//+------------------------------------------------------------------+
static void CMinCG::MinCGInitInternal(const int n,const double diffstep,
                                      CMinCGState &state)
  {
//--- create a variable
   int i=0;
//--- initialization
   state.m_n=n;
   state.m_diffstep=diffstep;
//--- function call
   MinCGSetCond(state,0,0,0,0);
//--- function call
   MinCGSetXRep(state,false);
//--- function call
   MinCGSetDRep(state,false);
//--- function call
   MinCGSetStpMax(state,0);
//--- function call
   MinCGSetCGType(state,-1);
//--- function call
   MinCGSetPrecDefault(state);
//--- allocation
   ArrayResizeAL(state.m_xk,n);
   ArrayResizeAL(state.m_dk,n);
   ArrayResizeAL(state.m_xn,n);
   ArrayResizeAL(state.m_dn,n);
   ArrayResizeAL(state.m_x,n);
   ArrayResizeAL(state.m_d,n);
   ArrayResizeAL(state.m_g,n);
   ArrayResizeAL(state.m_work0,n);
   ArrayResizeAL(state.m_work1,n);
   ArrayResizeAL(state.m_yk,n);
   ArrayResizeAL(state.m_s,n);
//--- copy
   for(i=0;i<=n-1;i++)
      state.m_s[i]=1.0;
  }
//+------------------------------------------------------------------+
//| NOTES:                                                           |
//| 1. This function has two different implementations: one which    |
//|    uses exact (analytical) user-supplied  gradient, and one which| 
//|    uses function value only and numerically differentiates       |
//|    function in order to obtain gradient.                         |
//|    Depending on the specific function used to create optimizer   |
//|    object (either MinCGCreate() for analytical gradient or       |
//|    MinCGCreateF() for numerical differentiation) you should      |
//|    choose appropriate variant of MinCGOptimize() - one which     |
//|    accepts function AND gradient or one which accepts function   |
//|    ONLY.                                                         |
//|    Be careful to choose variant of MinCGOptimize() which         |
//|    corresponds to your optimization scheme! Table below lists    |
//|    different combinations of callback (function/gradient) passed |
//|    to MinCGOptimize() and specific function used to create       |
//|    optimizer.                                                    |
//|                   |         USER PASSED TO MinCGOptimize()       |
//|    CREATED WITH   |  function only   |  function and gradient    |
//|    ------------------------------------------------------------  |
//|    MinCGCreateF() |     work                FAIL                 |
//|    MinCGCreate()  |     FAIL                work                 |
//|    Here "FAIL" denotes inappropriate combinations of optimizer   |
//|    creation function and MinCGOptimize() version. Attemps to use |
//|    such combination (for example, to create optimizer with       |
//|    MinCGCreateF() and to pass gradient information to            |
//|    MinCGOptimize()) will lead to exception being thrown. Either  |
//|    you did not pass gradient when it WAS needed or you passed    |
//|    gradient when it was NOT needed.                              |
//+------------------------------------------------------------------+
static bool CMinCG::MinCGIteration(CMinCGState &state)
  {
//--- create variables
   int    n=0;
   int    i=0;
   double betak=0;
   double v=0;
   double vv=0;
   int    i_=0;
//--- This code initializes locals by:
//--- * random values determined during code
//---   generation - on first subroutine call
//--- * values from previous call - on subsequent calls
   if(state.m_rstate.stage>=0)
     {
      //--- initialization
      n=state.m_rstate.ia[0];
      i=state.m_rstate.ia[1];
      betak=state.m_rstate.ra[0];
      v=state.m_rstate.ra[1];
      vv=state.m_rstate.ra[2];
     }
   else
     {
      //--- initialization
      n=-983;
      i=-989;
      betak=-834;
      v=900;
      vv=-287;
     }
//--- check
   if(state.m_rstate.stage==0)
     {
      //--- change value
      state.m_needfg=false;
      //--- function call, return result
      return(Func_lbl_18(state,n,i,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==1)
     {
      //--- change values
      state.m_fbase=state.m_f;
      i=0;
      //--- function call, return result
      return(Func_lbl_19(state,n,i,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==2)
     {
      //--- change values
      state.m_fm2=state.m_f;
      state.m_x[i]=v-0.5*state.m_diffstep*state.m_s[i];
      state.m_rstate.stage=3;
      //--- Saving state
      Func_lbl_rcomm(state,n,i,betak,v,vv);
      return(true);
     }
//--- check
   if(state.m_rstate.stage==3)
     {
      //--- change values
      state.m_fm1=state.m_f;
      state.m_x[i]=v+0.5*state.m_diffstep*state.m_s[i];
      state.m_rstate.stage=4;
      //--- Saving state
      Func_lbl_rcomm(state,n,i,betak,v,vv);
      return(true);
     }
//--- check
   if(state.m_rstate.stage==4)
     {
      //--- change values
      state.m_fp1=state.m_f;
      state.m_x[i]=v+state.m_diffstep*state.m_s[i];
      state.m_rstate.stage=5;
      //--- Saving state
      Func_lbl_rcomm(state,n,i,betak,v,vv);
      return(true);
     }
//--- check
   if(state.m_rstate.stage==5)
     {
      //--- change values
      state.m_fp2=state.m_f;
      state.m_x[i]=v;
      state.m_g[i]=(8*(state.m_fp1-state.m_fm1)-(state.m_fp2-state.m_fm2))/(6*state.m_diffstep*state.m_s[i]);
      i=i+1;
      //--- function call, return result
      return(Func_lbl_19(state,n,i,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==6)
     {
      //--- change value
      state.m_algpowerup=false;
      //--- function call, return result
      return(Func_lbl_22(state,n,i,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==7)
     {
      //--- change value
      state.m_xupdated=false;
      //--- function call, return result
      return(Func_lbl_24(state,n,i,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==8)
     {
      //--- change value
      state.m_lsstart=false;
      //--- function call, return result
      return(Func_lbl_28(state,n,i,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==9)
     {
      //--- change value
      state.m_needfg=false;
      //--- function call, return result
      return(Func_lbl_33(state,n,i,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==10)
     {
      //--- change values
      state.m_fbase=state.m_f;
      i=0;
      //--- function call, return result
      return(Func_lbl_34(state,n,i,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==11)
     {
      //--- change values
      state.m_fm2=state.m_f;
      state.m_x[i]=v-0.5*state.m_diffstep*state.m_s[i];
      state.m_rstate.stage=12;
      //--- Saving state
      Func_lbl_rcomm(state,n,i,betak,v,vv);
      return(true);
     }
//--- check
   if(state.m_rstate.stage==12)
     {
      //--- change values
      state.m_fm1=state.m_f;
      state.m_x[i]=v+0.5*state.m_diffstep*state.m_s[i];
      state.m_rstate.stage=13;
      //--- Saving state
      Func_lbl_rcomm(state,n,i,betak,v,vv);
      return(true);
     }
//--- check
   if(state.m_rstate.stage==13)
     {
      //--- change values
      state.m_fp1=state.m_f;
      state.m_x[i]=v+state.m_diffstep*state.m_s[i];
      state.m_rstate.stage=14;
      //--- Saving state
      Func_lbl_rcomm(state,n,i,betak,v,vv);
      return(true);
     }
//--- check
   if(state.m_rstate.stage==14)
     {
      //--- change values
      state.m_fp2=state.m_f;
      state.m_x[i]=v;
      state.m_g[i]=(8*(state.m_fp1-state.m_fm1)-(state.m_fp2-state.m_fm2))/(6*state.m_diffstep*state.m_s[i]);
      i=i+1;
      //--- function call, return result
      return(Func_lbl_34(state,n,i,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==15)
     {
      //--- change value
      state.m_lsend=false;
      //--- function call, return result
      return(Func_lbl_37(state,n,i,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==16)
     {
      //--- change value
      state.m_xupdated=false;
      //--- function call, return result
      return(Func_lbl_39(state,n,i,betak,v,vv));
     }
//--- Routine body
//--- Prepare
   n=state.m_n;
   state.m_repterminationtype=0;
   state.m_repiterationscount=0;
   state.m_repnfev=0;
   state.m_debugrestartscount=0;
//--- Preparations continue:
//--- * set XK
//--- * calculate F/G
//--- * set DK to -G
//--- * powerup algo (it may change preconditioner)
//--- * apply preconditioner to DK
//--- * report update of X
//--- * check stopping conditions for G
   for(i_=0;i_<=n-1;i_++)
      state.m_xk[i_]=state.m_x[i_];
//--- change value
   state.m_terminationneeded=false;
//--- function call
   ClearRequestFields(state);
//--- check
   if(state.m_diffstep!=0.0)
     {
      state.m_needf=true;
      state.m_rstate.stage=1;
      //--- Saving state
      Func_lbl_rcomm(state,n,i,betak,v,vv);
      //--- return result
      return(true);
     }
//--- change values
   state.m_needfg=true;
   state.m_rstate.stage=0;
//--- Saving state
   Func_lbl_rcomm(state,n,i,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinCGIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static void CMinCG::Func_lbl_rcomm(CMinCGState &state,int n,int i,
                                   double betak,double v,double vv)
  {
//--- save
   state.m_rstate.ia[0]=n;
   state.m_rstate.ia[1]=i;
   state.m_rstate.ra[0]=betak;
   state.m_rstate.ra[1]=v;
   state.m_rstate.ra[2]=vv;
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinCGIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinCG::Func_lbl_18(CMinCGState &state,int &n,int &i,
                                double &betak,double &v,double &vv)
  {
//--- check
   if(!state.m_drep)
      return(Func_lbl_22(state,n,i,betak,v,vv));
//--- Report algorithm powerup (if needed)
   ClearRequestFields(state);
//--- change values
   state.m_algpowerup=true;
   state.m_rstate.stage=6;
//--- Saving state
   Func_lbl_rcomm(state,n,i,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinCGIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinCG::Func_lbl_19(CMinCGState &state,int &n,int &i,
                                double &betak,double &v,double &vv)
  {
//--- check
   if(i>n-1)
     {
      state.m_f=state.m_fbase;
      state.m_needf=false;
      //--- function call, return result
      return(Func_lbl_18(state,n,i,betak,v,vv));
     }
//--- change values
   v=state.m_x[i];
   state.m_x[i]=v-state.m_diffstep*state.m_s[i];
   state.m_rstate.stage=2;
//--- Saving state
   Func_lbl_rcomm(state,n,i,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinCGIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinCG::Func_lbl_22(CMinCGState &state,int &n,int &i,
                                double &betak,double &v,double &vv)
  {
//--- function call
   COptServ::TrimPrepare(state.m_f,state.m_trimthreshold);
   for(int i_=0;i_<=n-1;i_++)
      state.m_dk[i_]=-state.m_g[i_];
//--- function call
   PreconditionedMultiply(state,state.m_dk,state.m_work0,state.m_work1);
//--- check
   if(!state.m_xrep)
      return(Func_lbl_24(state,n,i,betak,v,vv));
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=7;
//--- Saving state
   Func_lbl_rcomm(state,n,i,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinCGIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinCG::Func_lbl_24(CMinCGState &state,int &n,int &i,
                                double &betak,double &v,double &vv)
  {
//--- check
   if(state.m_terminationneeded)
     {
      for(int i_=0;i_<=n-1;i_++)
         state.m_xn[i_]=state.m_xk[i_];
      state.m_repterminationtype=8;
      //--- return result
      return(false);
     }
//--- change value
   v=0;
   for(i=0;i<=n-1;i++)
      v=v+CMath::Sqr(state.m_g[i]*state.m_s[i]);
//--- check
   if(MathSqrt(v)<=state.m_epsg)
     {
      for(int i_=0;i_<=n-1;i_++)
         state.m_xn[i_]=state.m_xk[i_];
      state.m_repterminationtype=4;
      //--- return result
      return(false);
     }
//--- change values
   state.m_repnfev=1;
   state.m_k=0;
   state.m_fold=state.m_f;
//--- Choose initial step.
//--- Apply preconditioner,if we have something other than default.
   if(state.m_prectype==2||state.m_prectype==3)
     {
      //--- because we use preconditioner,step length must be equal
      //--- to the norm of DK
      v=0.0;
      for(int i_=0;i_<=n-1;i_++)
         v+=state.m_dk[i_]*state.m_dk[i_];
      state.m_laststep=MathSqrt(v);
     }
   else
     {
      //--- No preconditioner is used,we try to use suggested step
      if(state.m_suggestedstep>0.0)
         state.m_laststep=state.m_suggestedstep;
      else
        {
         //--- change value
         v=0.0;
         for(int i_=0;i_<=n-1;i_++)
            v+=state.m_g[i_]*state.m_g[i_];
         v=MathSqrt(v);
         //--- check
         if(state.m_stpmax==0.0)
            state.m_laststep=MathMin(1.0/v,1);
         else
            state.m_laststep=MathMin(1.0/v,state.m_stpmax);
        }
     }
//--- Main cycle
   state.m_rstimer=m_rscountdownlen;
//--- function call, return result
   return(Func_lbl_26(state,n,i,betak,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinCGIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinCG::Func_lbl_26(CMinCGState &state,int &n,int &i,
                                double &betak,double &v,double &vv)
  {
//--- * clear reset flag
//--- * clear termination flag
//--- * store G[k] for later calculation of Y[k]
//--- * prepare starting point and direction and step length for line search
   state.m_innerresetneeded=false;
   state.m_terminationneeded=false;
   for(int i_=0;i_<=n-1;i_++)
      state.m_yk[i_]=-state.m_g[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_d[i_]=state.m_dk[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xk[i_];
//--- change values
   state.m_mcstage=0;
   state.m_stp=1.0;
//--- function call
   CLinMin::LinMinNormalized(state.m_d,state.m_stp,n);
//--- check
   if(state.m_laststep!=0.0)
      state.m_stp=state.m_laststep;
   state.m_curstpmax=state.m_stpmax;
//--- Report beginning of line search (if needed)
//--- Terminate algorithm,if user request was detected
   if(!state.m_drep)
      return(Func_lbl_28(state,n,i,betak,v,vv));
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_lsstart=true;
   state.m_rstate.stage=8;
//--- Saving state
   Func_lbl_rcomm(state,n,i,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinCGIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinCG::Func_lbl_28(CMinCGState &state,int &n,int &i,
                                double &betak,double &v,double &vv)
  {
//--- check
   if(state.m_terminationneeded)
     {
      for(int i_=0;i_<=n-1;i_++)
         state.m_xn[i_]=state.m_x[i_];
      state.m_repterminationtype=8;
      //--- return result
      return(false);
     }
//--- Minimization along D
   CLinMin::MCSrch(n,state.m_x,state.m_f,state.m_g,state.m_d,state.m_stp,state.m_curstpmax,m_gtol,state.m_mcinfo,state.m_nfev,state.m_work0,state.m_lstate,state.m_mcstage);
//--- function call, return result
   return(Func_lbl_30(state,n,i,betak,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinCGIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinCG::Func_lbl_30(CMinCGState &state,int &n,int &i,
                                double &betak,double &v,double &vv)
  {
//--- check
   if(state.m_mcstage==0)
      return(Func_lbl_31(state,n,i,betak,v,vv));
//--- Calculate function/gradient using either
//--- analytical gradient supplied by user
//--- or finite difference approximation.
//--- "Trim" function in order to handle near-singularity points.
   ClearRequestFields(state);
//--- check
   if((double)(state.m_diffstep)!=0.0)
     {
      state.m_needf=true;
      state.m_rstate.stage=10;
      //--- Saving state
      Func_lbl_rcomm(state,n,i,betak,v,vv);
      //--- return result
      return(true);
     }
//--- change values
   state.m_needfg=true;
   state.m_rstate.stage=9;
//--- Saving state
   Func_lbl_rcomm(state,n,i,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinCGIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinCG::Func_lbl_31(CMinCGState &state,int &n,int &i,
                                double &betak,double &v,double &vv)
  {
//--- * report end of line search
//--- * store current point to XN
//--- * report iteration
//--- * terminate algorithm if user request was detected
   if(!state.m_drep)
      return(Func_lbl_37(state,n,i,betak,v,vv));
//--- Report end of line search (if needed)
   ClearRequestFields(state);
//--- change values
   state.m_lsend=true;
   state.m_rstate.stage=15;
//--- Saving state
   Func_lbl_rcomm(state,n,i,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinCGIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinCG::Func_lbl_33(CMinCGState &state,int &n,int &i,
                                double &betak,double &v,double &vv)
  {
//--- function call
   COptServ::TrimFunction(state.m_f,state.m_g,n,state.m_trimthreshold);
//--- Call MCSRCH again
   CLinMin::MCSrch(n,state.m_x,state.m_f,state.m_g,state.m_d,state.m_stp,state.m_curstpmax,m_gtol,state.m_mcinfo,state.m_nfev,state.m_work0,state.m_lstate,state.m_mcstage);
//--- function call, return result
   return(Func_lbl_30(state,n,i,betak,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinCGIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinCG::Func_lbl_34(CMinCGState &state,int &n,int &i,
                                double &betak,double &v,double &vv)
  {
//--- check
   if(i>n-1)
     {
      state.m_f=state.m_fbase;
      state.m_needf=false;
      //--- function call, return result
      return(Func_lbl_33(state,n,i,betak,v,vv));
     }
//--- change values
   v=state.m_x[i];
   state.m_x[i]=v-state.m_diffstep*state.m_s[i];
   state.m_rstate.stage=11;
//--- Saving state
   Func_lbl_rcomm(state,n,i,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinCGIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinCG::Func_lbl_37(CMinCGState &state,int &n,int &i,
                                double &betak,double &v,double &vv)
  {
//--- copy
   for(int i_=0;i_<=n-1;i_++)
      state.m_xn[i_]=state.m_x[i_];
//--- check
   if(!state.m_xrep)
      return(Func_lbl_39(state,n,i,betak,v,vv));
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=16;
//--- Saving state
   Func_lbl_rcomm(state,n,i,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinCGIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinCG::Func_lbl_39(CMinCGState &state,int &n,int &i,
                                double &betak,double &v,double &vv)
  {
//--- check
   if(state.m_terminationneeded)
     {
      for(int i_=0;i_<=n-1;i_++)
         state.m_xn[i_]=state.m_x[i_];
      state.m_repterminationtype=8;
      //--- return result
      return(false);
     }
//--- Line search is finished.
//--- * calculate BetaK
//--- * calculate DN
//--- * update timers
//--- * calculate step length
   if(state.m_mcinfo==1 && !state.m_innerresetneeded)
     {
      //--- Standard Wolfe conditions hold
      //--- Calculate Y[K] and D[K]'*Y[K]
      for(int i_=0;i_<=n-1;i_++)
         state.m_yk[i_]=state.m_yk[i_]+state.m_g[i_];
      //--- change value
      vv=0.0;
      for(int i_=0;i_<=n-1;i_++)
         vv+=state.m_yk[i_]*state.m_dk[i_];
      //--- Calculate BetaK according to DY formula
      v=PreconditionedMultiply2(state,state.m_g,state.m_g,state.m_work0,state.m_work1);
      state.m_betady=v/vv;
      //--- Calculate BetaK according to HS formula
      v=PreconditionedMultiply2(state,state.m_g,state.m_yk,state.m_work0,state.m_work1);
      state.m_betahs=v/vv;
      //--- Choose BetaK
      if(state.m_cgtype==0)
         betak=state.m_betady;
      //--- check
      if(state.m_cgtype==1)
         betak=MathMax(0,MathMin(state.m_betady,state.m_betahs));
     }
   else
     {
      //--- Something is wrong (may be function is too wild or too flat)
      //--- or we just have to restart algo.
      //--- We'll set BetaK=0,which will restart CG algorithm.
      //--- We can stop later (during normal checks) if stopping conditions are met.
      betak=0;
      state.m_debugrestartscount=state.m_debugrestartscount+1;
     }
//--- check
   if(state.m_repiterationscount>0&&state.m_repiterationscount%(3+n)==0)
     {
      //--- clear Beta every N iterations
      betak=0;
     }
//--- check
   if(state.m_mcinfo==1||state.m_mcinfo==5)
      state.m_rstimer=m_rscountdownlen;
   else
      state.m_rstimer=state.m_rstimer-1;
   for(int i_=0;i_<=n-1;i_++)
      state.m_dn[i_]=-state.m_g[i_];
//--- function call
   PreconditionedMultiply(state,state.m_dn,state.m_work0,state.m_work1);
   for(int i_=0;i_<=n-1;i_++)
      state.m_dn[i_]=state.m_dn[i_]+betak*state.m_dk[i_];
//--- change values
   state.m_laststep=0;
   state.m_lastscaledstep=0.0;
   for(i=0;i<=n-1;i++)
     {
      state.m_laststep=state.m_laststep+CMath::Sqr(state.m_d[i]);
      state.m_lastscaledstep=state.m_lastscaledstep+CMath::Sqr(state.m_d[i]/state.m_s[i]);
     }
//--- change values
   state.m_laststep=state.m_stp*MathSqrt(state.m_laststep);
   state.m_lastscaledstep=state.m_stp*MathSqrt(state.m_lastscaledstep);
//--- Update information.
//--- Check stopping conditions.
   state.m_repnfev=state.m_repnfev+state.m_nfev;
   state.m_repiterationscount=state.m_repiterationscount+1;
//--- check
   if(state.m_repiterationscount>=state.m_maxits&&state.m_maxits>0)
     {
      //--- Too many iterations
      state.m_repterminationtype=5;
      //--- return result
      return(false);
     }
//--- change value
   v=0;
   for(i=0;i<=n-1;i++)
      v=v+CMath::Sqr(state.m_g[i]*state.m_s[i]);
//--- check
   if(MathSqrt(v)<=state.m_epsg)
     {
      //--- Gradient is small enough
      state.m_repterminationtype=4;
      //--- return result
      return(false);
     }
//--- check
   if(!state.m_innerresetneeded)
     {
      //--- These conditions are checked only when no inner reset was requested by user
      if(state.m_fold-state.m_f<=state.m_epsf*MathMax(MathAbs(state.m_fold),MathMax(MathAbs(state.m_f),1.0)))
        {
         //--- F(k+1)-F(k) is small enough
         state.m_repterminationtype=1;
         //--- return result
         return(false);
        }
      //--- check
      if(state.m_lastscaledstep<=state.m_epsx)
        {
         //--- X(k+1)-X(k) is small enough
         state.m_repterminationtype=2;
         //--- return result
         return(false);
        }
     }
//--- check
   if(state.m_rstimer<=0)
     {
      //--- Too many subsequent restarts
      state.m_repterminationtype=7;
      //--- return result
      return(false);
     }
//--- Shift Xk/Dk,update other information
   for(int i_=0;i_<=n-1;i_++)
      state.m_xk[i_]=state.m_xn[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_dk[i_]=state.m_dn[i_];
//--- change values
   state.m_fold=state.m_f;
   state.m_k=state.m_k+1;
//--- function call, return result
   return(Func_lbl_26(state,n,i,betak,v,vv));
  }
//+------------------------------------------------------------------+
//| This object stores nonlinear optimizer state.                    |
//| You should use functions provided by MinBLEIC subpackage to work |
//| with this object                                                 |
//+------------------------------------------------------------------+
class CMinBLEICState
  {
public:
   //--- variables
   int               m_nmain;
   int               m_nslack;
   double            m_innerepsg;
   double            m_innerepsf;
   double            m_innerepsx;
   double            m_outerepsx;
   double            m_outerepsi;
   int               m_maxits;
   bool              m_xrep;
   double            m_stpmax;
   double            m_diffstep;
   int               m_prectype;
   double            m_f;
   bool              m_needf;
   bool              m_needfg;
   bool              m_xupdated;
   RCommState        m_rstate;
   int               m_repinneriterationscount;
   int               m_repouteriterationscount;
   int               m_repnfev;
   int               m_repterminationtype;
   double            m_repdebugeqerr;
   double            m_repdebugfs;
   double            m_repdebugff;
   double            m_repdebugdx;
   int               m_itsleft;
   double            m_trimthreshold;
   int               m_cecnt;
   int               m_cedim;
   double            m_v0;
   double            m_v1;
   double            m_v2;
   double            m_t;
   double            m_errfeas;
   double            m_gnorm;
   double            m_mpgnorm;
   double            m_mba;
   int               m_variabletofreeze;
   double            m_valuetofreeze;
   double            m_fbase;
   double            m_fm2;
   double            m_fm1;
   double            m_fp1;
   double            m_fp2;
   double            m_xm1;
   double            m_xp1;
   CMinCGState       m_cgstate;
   CMinCGReport      m_cgrep;
   int               m_optdim;
   //--- arrays
   double            m_diaghoriginal[];
   double            m_diagh[];
   double            m_x[];
   double            m_g[];
   double            m_xcur[];
   double            m_xprev[];
   double            m_xstart[];
   double            m_xend[];
   double            m_lastg[];
   int               m_ct[];
   double            m_xe[];
   bool              m_hasbndl[];
   bool              m_hasbndu[];
   double            m_bndloriginal[];
   double            m_bnduoriginal[];
   double            m_bndleffective[];
   double            m_bndueffective[];
   bool              m_activeconstraints[];
   double            m_constrainedvalues[];
   double            m_transforms[];
   double            m_seffective[];
   double            m_soriginal[];
   double            m_w[];
   double            m_tmp0[];
   double            m_tmp1[];
   double            m_tmp2[];
   double            m_r[];
   //--- matrix
   CMatrixDouble     m_ceoriginal;
   CMatrixDouble     m_ceeffective;
   CMatrixDouble     m_cecurrent;
   CMatrixDouble     m_lmmatrix;
   //--- constructor, destructor
                     CMinBLEICState(void);
                    ~CMinBLEICState(void);
   //--- copy
   void              Copy(CMinBLEICState &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinBLEICState::CMinBLEICState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinBLEICState::~CMinBLEICState(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMinBLEICState::Copy(CMinBLEICState &obj)
  {
//--- copy variables
   m_nmain=obj.m_nmain;
   m_nslack=obj.m_nslack;
   m_innerepsg=obj.m_innerepsg;
   m_innerepsf=obj.m_innerepsf;
   m_innerepsx=obj.m_innerepsx;
   m_outerepsx=obj.m_outerepsx;
   m_outerepsi=obj.m_outerepsi;
   m_maxits=obj.m_maxits;
   m_xrep=obj.m_xrep;
   m_stpmax=obj.m_stpmax;
   m_diffstep=obj.m_diffstep;
   m_prectype=obj.m_prectype;
   m_f=obj.m_f;
   m_needf=obj.m_needf;
   m_needfg=obj.m_needfg;
   m_xupdated=obj.m_xupdated;
   m_repinneriterationscount=obj.m_repinneriterationscount;
   m_repouteriterationscount=obj.m_repouteriterationscount;
   m_repnfev=obj.m_repnfev;
   m_repterminationtype=obj.m_repterminationtype;
   m_repdebugeqerr=obj.m_repdebugeqerr;
   m_repdebugfs=obj.m_repdebugfs;
   m_repdebugff=obj.m_repdebugff;
   m_repdebugdx=obj.m_repdebugdx;
   m_itsleft=obj.m_itsleft;
   m_trimthreshold=obj.m_trimthreshold;
   m_cecnt=obj.m_cecnt;
   m_cedim=obj.m_cedim;
   m_v0=obj.m_v0;
   m_v1=obj.m_v1;
   m_v2=obj.m_v2;
   m_t=obj.m_t;
   m_errfeas=obj.m_errfeas;
   m_gnorm=obj.m_gnorm;
   m_mpgnorm=obj.m_mpgnorm;
   m_mba=obj.m_mba;
   m_variabletofreeze=obj.m_variabletofreeze;
   m_valuetofreeze=obj.m_valuetofreeze;
   m_fbase=obj.m_fbase;
   m_fm2=obj.m_fm2;
   m_fm1=obj.m_fm1;
   m_fp1=obj.m_fp1;
   m_fp2=obj.m_fp2;
   m_xm1=obj.m_xm1;
   m_xp1=obj.m_xp1;
   m_optdim=obj.m_optdim;
   m_rstate.Copy(obj.m_rstate);
   m_cgstate.Copy(obj.m_cgstate);
   m_cgrep.Copy(obj.m_cgrep);
//--- copy arrays
   ArrayCopy(m_diaghoriginal,obj.m_diaghoriginal);
   ArrayCopy(m_diagh,obj.m_diagh);
   ArrayCopy(m_x,obj.m_x);
   ArrayCopy(m_g,obj.m_g);
   ArrayCopy(m_xcur,obj.m_xcur);
   ArrayCopy(m_xprev,obj.m_xprev);
   ArrayCopy(m_xstart,obj.m_xstart);
   ArrayCopy(m_xend,obj.m_xend);
   ArrayCopy(m_lastg,obj.m_lastg);
   ArrayCopy(m_ct,obj.m_ct);
   ArrayCopy(m_xe,obj.m_xe);
   ArrayCopy(m_hasbndl,obj.m_hasbndl);
   ArrayCopy(m_hasbndu,obj.m_hasbndu);
   ArrayCopy(m_bndloriginal,obj.m_bndloriginal);
   ArrayCopy(m_bnduoriginal,obj.m_bnduoriginal);
   ArrayCopy(m_bndleffective,obj.m_bndleffective);
   ArrayCopy(m_bndueffective,obj.m_bndueffective);
   ArrayCopy(m_activeconstraints,obj.m_activeconstraints);
   ArrayCopy(m_constrainedvalues,obj.m_constrainedvalues);
   ArrayCopy(m_transforms,obj.m_transforms);
   ArrayCopy(m_seffective,obj.m_seffective);
   ArrayCopy(m_soriginal,obj.m_soriginal);
   ArrayCopy(m_w,obj.m_w);
   ArrayCopy(m_tmp0,obj.m_tmp0);
   ArrayCopy(m_tmp1,obj.m_tmp1);
   ArrayCopy(m_tmp2,obj.m_tmp2);
   ArrayCopy(m_r,obj.m_r);
//--- copy matrix
   m_ceoriginal=obj.m_ceoriginal;
   m_ceeffective=obj.m_ceeffective;
   m_cecurrent=obj.m_cecurrent;
   m_lmmatrix=obj.m_lmmatrix;
  }
//+------------------------------------------------------------------+
//| This object stores nonlinear optimizer state.                    |
//| You should use functions provided by MinBLEIC subpackage to work |
//| with this object                                                 |
//+------------------------------------------------------------------+
class CMinBLEICStateShell
  {
private:
   CMinBLEICState    m_innerobj;
public:
   //--- constructors, destructor
                     CMinBLEICStateShell(void);
                     CMinBLEICStateShell(CMinBLEICState &obj);
                    ~CMinBLEICStateShell(void);
   //--- methods
   bool              GetNeedF(void);
   void              SetNeedF(const bool b);
   bool              GetNeedFG(void);
   void              SetNeedFG(const bool b);
   bool              GetXUpdated(void);
   void              SetXUpdated(const bool b);
   double            GetF(void);
   void              SetF(const double d);
   CMinBLEICState   *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinBLEICStateShell::CMinBLEICStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMinBLEICStateShell::CMinBLEICStateShell(CMinBLEICState &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinBLEICStateShell::~CMinBLEICStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needf                          |
//+------------------------------------------------------------------+
bool CMinBLEICStateShell::GetNeedF(void)
  {
//--- return result
   return(m_innerobj.m_needf);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needf                         |
//+------------------------------------------------------------------+
void CMinBLEICStateShell::SetNeedF(const bool b)
  {
//--- change value
   m_innerobj.m_needf=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needfg                         |
//+------------------------------------------------------------------+
bool CMinBLEICStateShell::GetNeedFG(void)
  {
//--- return result
   return(m_innerobj.m_needfg);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needfg                        |
//+------------------------------------------------------------------+
void CMinBLEICStateShell::SetNeedFG(const bool b)
  {
//--- change value
   m_innerobj.m_needfg=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable xupdated                       |
//+------------------------------------------------------------------+
bool CMinBLEICStateShell::GetXUpdated(void)
  {
//--- return result
   return(m_innerobj.m_xupdated);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable xupdated                      |
//+------------------------------------------------------------------+
void CMinBLEICStateShell::SetXUpdated(const bool b)
  {
//--- change value
   m_innerobj.m_xupdated=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable f                              |
//+------------------------------------------------------------------+
double CMinBLEICStateShell::GetF(void)
  {
//--- return result
   return(m_innerobj.m_f);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable f                             |
//+------------------------------------------------------------------+
void CMinBLEICStateShell::SetF(const double d)
  {
//--- change value
   m_innerobj.m_f=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMinBLEICState *CMinBLEICStateShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| This structure stores optimization report:                       |
//| * InnerIterationsCount      number of inner iterations           |
//| * OuterIterationsCount      number of outer iterations           |
//| * NFEV                      number of gradient evaluations       |
//| * TerminationType           termination type (see below)         |
//| TERMINATION CODES                                                |
//| TerminationType field contains completion code,which can be:     |
//|   -10   unsupported combination of algorithm settings:           |
//|         1) StpMax is set to non-zero value,                      |
//|         AND 2) non-default preconditioner is used.               |
//|         You can't use both features at the same moment,          |
//|         so you have to choose one of them (and to turn           |
//|         off another one).                                        |
//|   -3    inconsistent constraints. Feasible point is              |
//|         either nonexistent or too hard to find. Try to           |
//|         restart optimizer with better initial                    |
//|         approximation                                            |
//|    4    conditions on constraints are fulfilled                  |
//|         with error less than or equal to EpsC                    |
//|    5    MaxIts steps was taken                                   |
//|    7    stopping conditions are too stringent,                   |
//|         further improvement is impossible,                       |
//|         X contains best point found so far.                      |
//| ADDITIONAL FIELDS                                                |
//| There are additional fields which can be used for debugging:     |
//| * DebugEqErr                error in the equality constraints    |
//|                             (2-norm)                             |
//| * DebugFS                   f,calculated at projection of initial|
//|                             point to the feasible set            |
//| * DebugFF                   f,calculated at the final point      |
//| * DebugDX                   |X_start-X_final|                    |
//+------------------------------------------------------------------+
class CMinBLEICReport
  {
public:
   //--- variables
   int               m_inneriterationscount;
   int               m_outeriterationscount;
   int               m_nfev;
   int               m_terminationtype;
   double            m_debugeqerr;
   double            m_debugfs;
   double            m_debugff;
   double            m_debugdx;
   //--- constructor, destructor
                     CMinBLEICReport(void);
                    ~CMinBLEICReport(void);
   //--- copy
   void              Copy(CMinBLEICReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinBLEICReport::CMinBLEICReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinBLEICReport::~CMinBLEICReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMinBLEICReport::Copy(CMinBLEICReport &obj)
  {
//--- copy variables
   m_inneriterationscount=obj.m_inneriterationscount;
   m_outeriterationscount=obj.m_outeriterationscount;
   m_nfev=obj.m_nfev;
   m_terminationtype=obj.m_terminationtype;
   m_debugeqerr=obj.m_debugeqerr;
   m_debugfs=obj.m_debugfs;
   m_debugff=obj.m_debugff;
   m_debugdx=obj.m_debugdx;
  }
//+------------------------------------------------------------------+
//| This structure stores optimization report:                       |
//| * InnerIterationsCount      number of inner iterations           |
//| * OuterIterationsCount      number of outer iterations           |
//| * NFEV                      number of gradient evaluations       |
//| * TerminationType           termination type (see below)         |
//| TERMINATION CODES                                                |
//| TerminationType field contains completion code,which can be:     |
//|   -10   unsupported combination of algorithm settings:           |
//|         1) StpMax is set to non-zero value,                      |
//|         AND 2) non-default preconditioner is used.               |
//|         You can't use both features at the same moment,          |
//|         so you have to choose one of them (and to turn           |
//|         off another one).                                        |
//|   -3    inconsistent constraints. Feasible point is              |
//|         either nonexistent or too hard to find. Try to           |
//|         restart optimizer with better initial                    |
//|         approximation                                            |
//|    4    conditions on constraints are fulfilled                  |
//|         with error less than or equal to EpsC                    |
//|    5    MaxIts steps was taken                                   |
//|    7    stopping conditions are too stringent,                   |
//|         further improvement is impossible,                       |
//|         X contains best point found so far.                      |
//| ADDITIONAL FIELDS                                                |
//| There are additional fields which can be used for debugging:     |
//| * DebugEqErr                error in the equality constraints    |
//|                             (2-norm)                             |
//| * DebugFS                   f,calculated at projection of initial|
//|                             point to the feasible set            |
//| * DebugFF                   f,calculated at the final point      |
//| * DebugDX                   |X_start-X_final|                    |
//+------------------------------------------------------------------+
class CMinBLEICReportShell
  {
private:
   CMinBLEICReport   m_innerobj;
public:
   //--- constructors, destructor
                     CMinBLEICReportShell(void);
                     CMinBLEICReportShell(CMinBLEICReport &obj);
                    ~CMinBLEICReportShell(void);
   //--- methods
   int               GetInnerIterationsCount(void);
   void              SetInnerIterationsCount(const int i);
   int               GetOuterIterationsCount(void);
   void              SetOuterIterationsCount(const int i);
   int               GetNFev(void);
   void              SetNFev(const int i);
   int               GetTerminationType(void);
   void              SetTerminationType(const int i);
   double            GetDebugEqErr(void);
   void              SetDebugEqErr(const double d);
   double            GetDebugFS(void);
   void              SetDebugFS(const double d);
   double            GetDebugFF(void);
   void              SetDebugFF(const double d);
   double            GetDebugDX(void);
   void              SetDebugDX(const double d);
   CMinBLEICReport *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinBLEICReportShell::CMinBLEICReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMinBLEICReportShell::CMinBLEICReportShell(CMinBLEICReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinBLEICReportShell::~CMinBLEICReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable inneriterationscount           |
//+------------------------------------------------------------------+
int CMinBLEICReportShell::GetInnerIterationsCount(void)
  {
//--- return result
   return(m_innerobj.m_inneriterationscount);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable inneriterationscount          |
//+------------------------------------------------------------------+
void CMinBLEICReportShell::SetInnerIterationsCount(const int i)
  {
//--- change value
   m_innerobj.m_inneriterationscount=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable outeriterationscount           |
//+------------------------------------------------------------------+
int CMinBLEICReportShell::GetOuterIterationsCount(void)
  {
//--- return result
   return(m_innerobj.m_outeriterationscount);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable outeriterationscount          |
//+------------------------------------------------------------------+
void CMinBLEICReportShell::SetOuterIterationsCount(const int i)
  {
//--- change value
   m_innerobj.m_outeriterationscount=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable nfev                           |
//+------------------------------------------------------------------+
int CMinBLEICReportShell::GetNFev(void)
  {
//--- return result
   return(m_innerobj.m_nfev);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable nfev                          |
//+------------------------------------------------------------------+
void CMinBLEICReportShell::SetNFev(const int i)
  {
//--- change value
   m_innerobj.m_nfev=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable terminationtype                |
//+------------------------------------------------------------------+
int CMinBLEICReportShell::GetTerminationType(void)
  {
//--- return result
   return(m_innerobj.m_terminationtype);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable terminationtype               |
//+------------------------------------------------------------------+
void CMinBLEICReportShell::SetTerminationType(const int i)
  {
//--- change value
   m_innerobj.m_terminationtype=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable debugeqerr                     |
//+------------------------------------------------------------------+
double CMinBLEICReportShell::GetDebugEqErr(void)
  {
//--- return result
   return(m_innerobj.m_debugeqerr);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable debugeqerr                    |
//+------------------------------------------------------------------+
void CMinBLEICReportShell::SetDebugEqErr(const double d)
  {
//--- change value
   m_innerobj.m_debugeqerr=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable debugfs                        |
//+------------------------------------------------------------------+
double CMinBLEICReportShell::GetDebugFS(void)
  {
//--- return result
   return(m_innerobj.m_debugfs);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable debugfs                       |
//+------------------------------------------------------------------+
void CMinBLEICReportShell::SetDebugFS(const double d)
  {
//--- change value
   m_innerobj.m_debugfs=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable debugff                        |
//+------------------------------------------------------------------+
double CMinBLEICReportShell::GetDebugFF(void)
  {
//--- return result
   return(m_innerobj.m_debugff);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable debugff                       |
//+------------------------------------------------------------------+
void CMinBLEICReportShell::SetDebugFF(const double d)
  {
//--- change value
   m_innerobj.m_debugff=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable debugdx                        |
//+------------------------------------------------------------------+
double CMinBLEICReportShell::GetDebugDX(void)
  {
//--- return result
   return(m_innerobj.m_debugdx);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable debugdx                       |
//+------------------------------------------------------------------+
void CMinBLEICReportShell::SetDebugDX(const double d)
  {
//--- change value
   m_innerobj.m_debugdx=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMinBLEICReport *CMinBLEICReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Bound constrained optimization with additional linear equality   |
//| and inequality constraints                                       |
//+------------------------------------------------------------------+
class CMinBLEIC
  {
private:
   //--- private methods
   static void       ClearRequestFields(CMinBLEICState &state);
   static void       UnscalePoint(CMinBLEICState &state,double &xscaled[],double &xunscaled[]);
   static void       ProjectPointAndUnscale(CMinBLEICState &state,double &xscaled[],double &xunscaled[],double &rscaled[],double &rnorm2);
   static void       ScaleGradientAndExpand(CMinBLEICState &state,double &gunscaled[],double &gscaled[]);
   static void       ModifyTargetFunction(CMinBLEICState &state,double &x[],double &r[],const double rnorm2,double &f,double &g[],double &gnorm,double &mpgnorm);
   static bool       AdditionalCheckForConstraints(CMinBLEICState &state,double &x[]);
   static void       RebuildCEXE(CMinBLEICState &state);
   static void       MakeGradientProjection(CMinBLEICState &state,double &pg[]);
   static bool       PrepareConstraintMatrix(CMinBLEICState &state,double &x[],double &g[],double &px[],double &pg[]);
   static void       MinBLEICInitInternal(const int n,double &x[],const double diffstep,CMinBLEICState &state);
   //--- auxiliary functions for MinBLEICIteration
   static void       Func_lbl_rcomm(CMinBLEICState &state,int nmain,int nslack,int m,int i,int j,bool b,double v,double vv);
   static bool       Func_lbl_14(CMinBLEICState &state,int &nmain,int &nslack,int &m,int &i,int &j,bool &b,double &v,double &vv);
   static bool       Func_lbl_15(CMinBLEICState &state,int &nmain,int &nslack,int &m,int &i,int &j,bool &b,double &v,double &vv);
   static bool       Func_lbl_16(CMinBLEICState &state,int &nmain,int &nslack,int &m,int &i,int &j,bool &b,double &v,double &vv);
   static bool       Func_lbl_17(CMinBLEICState &state,int &nmain,int &nslack,int &m,int &i,int &j,bool &b,double &v,double &vv);
   static bool       Func_lbl_18(CMinBLEICState &state,int &nmain,int &nslack,int &m,int &i,int &j,bool &b,double &v,double &vv);
   static bool       Func_lbl_19(CMinBLEICState &state,int &nmain,int &nslack,int &m,int &i,int &j,bool &b,double &v,double &vv);
   static bool       Func_lbl_22(CMinBLEICState &state,int &nmain,int &nslack,int &m,int &i,int &j,bool &b,double &v,double &vv);
   static bool       Func_lbl_23(CMinBLEICState &state,int &nmain,int &nslack,int &m,int &i,int &j,bool &b,double &v,double &vv);
   static bool       Func_lbl_31(CMinBLEICState &state,int &nmain,int &nslack,int &m,int &i,int &j,bool &b,double &v,double &vv);
public:
   //--- class constants
   static const double m_svdtol;
   static const double m_maxouterits;
   //--- constructor, destructor
                     CMinBLEIC(void);
                    ~CMinBLEIC(void);
   //--- public methods
   static void       MinBLEICCreate(const int n,double &x[],CMinBLEICState &state);
   static void       MinBLEICCreateF(const int n,double &x[],const double diffstep,CMinBLEICState &state);
   static void       MinBLEICSetBC(CMinBLEICState &state,double &bndl[],double &bndu[]);
   static void       MinBLEICSetLC(CMinBLEICState &state,CMatrixDouble &c,int &ct[],const int k);
   static void       MinBLEICSetInnerCond(CMinBLEICState &state,const double epsg,const double epsf,const double epsx);
   static void       MinBLEICSetOuterCond(CMinBLEICState &state,const double epsx,const double epsi);
   static void       MinBLEICSetScale(CMinBLEICState &state,double &s[]);
   static void       MinBLEICSetPrecDefault(CMinBLEICState &state);
   static void       MinBLEICSetPrecDiag(CMinBLEICState &state,double &d[]);
   static void       MinBLEICSetPrecScale(CMinBLEICState &state);
   static void       MinBLEICSetMaxIts(CMinBLEICState &state,const int maxits);
   static void       MinBLEICSetXRep(CMinBLEICState &state,const bool needxrep);
   static void       MinBLEICSetStpMax(CMinBLEICState &state,const double stpmax);
   static void       MinBLEICResults(CMinBLEICState &state,double &x[],CMinBLEICReport &rep);
   static void       MinBLEICResultsBuf(CMinBLEICState &state,double &x[],CMinBLEICReport &rep);
   static void       MinBLEICRestartFrom(CMinBLEICState &state,double &x[]);
   static bool       MinBLEICIteration(CMinBLEICState &state);
  };
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const double CMinBLEIC::m_svdtol=100;
const double CMinBLEIC::m_maxouterits=20;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinBLEIC::CMinBLEIC(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinBLEIC::~CMinBLEIC(void)
  {

  }
//+------------------------------------------------------------------+
//|                      BOUND CONSTRAINED OPTIMIZATION              |
//|        WITH ADDITIONAL LINEAR EQUALITY AND INEQUALITY CONSTRAINTS|
//| DESCRIPTION:                                                     |
//| The subroutine minimizes function F(x) of N arguments subject to |
//| any combination of:                                              |
//| * bound constraints                                              |
//| * linear inequality constraints                                  |
//| * linear equality constraints                                    |
//| REQUIREMENTS:                                                    |
//| * user must provide function value and gradient                  |
//| * starting point X0 must be feasible or                          |
//|   not too far away from the feasible set                         |
//| * grad(f) must be Lipschitz continuous on a level set:           |
//|   L = { x : f(x)<=f(x0) }                                        |
//| * function must be defined everywhere on the feasible set F      |
//| USAGE:                                                           |
//| Constrained optimization if far more complex than the            |
//| unconstrained one. Here we give very brief outline of the BLEIC  |
//| optimizer. We strongly recommend you to read examples in the     |
//| ALGLIB Reference Manual and to read ALGLIB User Guide on         |
//| optimization, which is available at                              |
//| http://www.alglib.net/optimization/                              |
//| 1. User initializes algorithm state with MinBLEICCreate() call   |
//| 2. USer adds boundary and/or linear constraints by calling       |
//|    MinBLEICSetBC() and MinBLEICSetLC() functions.                |
//| 3. User sets stopping conditions for underlying unconstrained    |
//|    solver with MinBLEICSetInnerCond() call.                      |
//|    This function controls accuracy of underlying optimization    |
//|    algorithm.                                                    |
//| 4. User sets stopping conditions for outer iteration by calling  |
//|    MinBLEICSetOuterCond() function.                              |
//|    This function controls handling of boundary and inequality    |
//|    constraints.                                                  |
//| 5. Additionally, user may set limit on number of internal        |
//|    iterations by MinBLEICSetMaxIts() call.                       |
//|    This function allows to prevent algorithm from looping        |
//|    forever.                                                      |
//| 6. User calls MinBLEICOptimize() function which takes algorithm  |
//|    state and pointer (delegate, etc.) to callback function       |
//|    which calculates F/G.                                         |
//| 7. User calls MinBLEICResults() to get solution                  |
//| 8. Optionally user may call MinBLEICRestartFrom() to solve       |
//|    another problem with same N but another starting point.       |
//|    MinBLEICRestartFrom() allows to reuse already initialized     |
//|    structure.                                                    |
//| INPUT PARAMETERS:                                                |
//|     N       -   problem dimension, N>0:                          |
//|                 * if given, only leading N elements of X are     |
//|                   used                                           |
//|                 * if not given, automatically determined from    |
//|                   size ofX                                       |
//|     X       -   starting point, array[N]:                        |
//|                 * it is better to set X to a feasible point      |
//|                 * but X can be infeasible, in which case         |
//|                   algorithm will try to find feasible point      |
//|                   first, using X as initial approximation.       |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure stores algorithm state                 |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICCreate(const int n,double &x[],CMinBLEICState &state)
  {
//--- create matrix
   CMatrixDouble c;
//--- create array
   int ct[];
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- function call
   MinBLEICInitInternal(n,x,0.0,state);
  }
//+------------------------------------------------------------------+
//| The subroutine is finite difference variant of MinBLEICCreate(). |
//| It uses finite differences in order to differentiate target      |
//| function.                                                        |
//| Description below contains information which is specific to this |
//| function only. We recommend to read comments on MinBLEICCreate() |
//| in order to get more information about creation of BLEIC         |
//| optimizer.                                                       |
//| INPUT PARAMETERS:                                                |
//|     N       -   problem dimension, N>0:                          |
//|                 * if given, only leading N elements of X are used|
//|                 * if not given, automatically determined from    |
//|                   size of X                                      |
//|     X       -   starting point, array[0..N-1].                   |
//|     DiffStep-   differentiation step, >0                         |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//| NOTES:                                                           |
//| 1. algorithm uses 4-point central formula for differentiation.   |
//| 2. differentiation step along I-th axis is equal to DiffStep*S[I]|
//|    where S[] is scaling vector which can be set by               |
//|    MinBLEICSetScale() call.                                      |
//| 3. we recommend you to use moderate values of differentiation    |
//|    step. Too large step will result in too large truncation      |
//|    errors, while too small step will result in too large         |
//|    numerical errors. 1.0E-6 can be good value to start with.     |
//| 4. Numerical differentiation is very inefficient - one  gradient |
//|    calculation needs 4*N function evaluations. This function will|
//|    work for any N - either small (1...10), moderate (10...100) or|
//|    large (100...). However, performance penalty will be too      |
//|    severe for any N's except for small ones.                     |
//|    We should also say that code which relies on numerical        |
//|    differentiation is less robust and precise. CG needs exact    |
//|    gradient values. Imprecise gradient may slow down convergence,|
//|    especially on highly nonlinear problems.                      |
//|    Thus we recommend to use this function for fast prototyping on|
//|    small - dimensional problems only, and to implement analytical|
//|    gradient as soon as possible.                                 |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICCreateF(const int n,double &x[],
                                       const double diffstep,
                                       CMinBLEICState &state)
  {
//--- create matrix
   CMatrixDouble c;
//--- create array
   int ct[];
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(diffstep),__FUNCTION__+": DiffStep is infinite or NaN!"))
      return;
//--- check
   if(!CAp::Assert(diffstep>0.0,__FUNCTION__+": DiffStep is non-positive!"))
      return;
//--- function call
   MinBLEICInitInternal(n,x,diffstep,state);
  }
//+------------------------------------------------------------------+
//| This function sets boundary constraints for BLEIC optimizer.     |
//| Boundary constraints are inactive by default (after initial      |
//| creation). They are preserved after algorithm restart with       |
//| MinBLEICRestartFrom().                                           |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure stores algorithm state                 |
//|     BndL    -   lower bounds, array[N].                          |
//|                 If some (all) variables are unbounded, you may   |
//|                 specify very small number or -INF.               |
//|     BndU    -   upper bounds, array[N].                          |
//|                 If some (all) variables are unbounded, you may   |
//|                 specify very large number or +INF.               |
//| NOTE 1: it is possible to specify BndL[i]=BndU[i]. In this case  |
//| I-th variable will be "frozen" at X[i]=BndL[i]=BndU[i].          |
//| NOTE 2: this solver has following useful properties:             |
//| * bound constraints are always satisfied exactly                 |
//| * function is evaluated only INSIDE area specified by bound      |
//|   constraints, even when numerical differentiation is used       |
//|   (algorithm adjusts nodes according to boundary constraints)    |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICSetBC(CMinBLEICState &state,double &bndl[],
                                     double &bndu[])
  {
//--- create variables
   int i=0;
   int n=0;
//--- initialization
   n=state.m_nmain;
//--- check
   if(!CAp::Assert(CAp::Len(bndl)>=n,__FUNCTION__+": Length(BndL)<N"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(bndu)>=n,__FUNCTION__+": Length(BndU)<N"))
      return;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(bndl[i]) || CInfOrNaN::IsNegativeInfinity(bndl[i]),__FUNCTION__+": BndL contains NAN or +INF"))
         return;
      //--- check
      if(!CAp::Assert(CMath::IsFinite(bndu[i]) || CInfOrNaN::IsPositiveInfinity(bndu[i]),__FUNCTION__+": BndU contains NAN or -INF"))
         return;
      //--- change values
      state.m_bndloriginal[i]=bndl[i];
      state.m_hasbndl[i]=CMath::IsFinite(bndl[i]);
      state.m_bnduoriginal[i]=bndu[i];
      state.m_hasbndu[i]=CMath::IsFinite(bndu[i]);
     }
  }
//+------------------------------------------------------------------+
//| This function sets linear constraints for BLEIC optimizer.       |
//| Linear constraints are inactive by default (after initial        |
//| creation). They are preserved after algorithm restart with       |
//| MinBLEICRestartFrom().                                           |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure previously allocated with              |
//|                MinBLEICCreate call.                              |
//|     C       -   linear constraints, array[K,N+1].                |
//|                 Each row of C represents one constraint, either  |
//|                 equality or inequality (see below):              |
//|                 * first N elements correspond to coefficients,   |
//|                 * last element corresponds to the right part.    |
//|                 All elements of C (including right part) must be |
//|                 finite.                                          |
//|     CT      -   type of constraints, array[K]:                   |
//|                 * if CT[i]>0, then I-th constraint is            |
//|                   C[i,*]*x >= C[i,n+1]                           |
//|                 * if CT[i]=0, then I-th constraint is            |
//|                   C[i,*]*x  = C[i,n+1]                           |
//|                 * if CT[i]<0, then I-th constraint is            |
//|                   C[i,*]*x <= C[i,n+1]                           |
//|     K       -   number of equality/inequality constraints, K>=0: |
//|                 * if given, only leading K elements of C/CT are  |
//|                   used                                           |
//|                 * if not given, automatically determined from    |
//|                   sizes of C/CT                                  |
//| NOTE 1: linear (non-bound) constraints are satisfied only        |
//| approximately:                                                   |
//| * there always exists some minor violation (about Epsilon in     |
//|   magnitude) due to rounding errors                              |
//| * numerical differentiation, if used, may lead to function       |
//|   evaluations outside of the feasible area, because algorithm    |
//|   does NOT change numerical differentiation formula according to |
//|   linear constraints.                                            |
//| If you want constraints to be satisfied exactly, try to          |
//| reformulate your problem in such manner that all constraints will|
//| become boundary ones (this kind of constraints is always         |
//| satisfied exactly, both in the final solution and in all         |
//| intermediate points).                                            |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICSetLC(CMinBLEICState &state,CMatrixDouble &c,
                                     int &ct[],const int k)
  {
//--- create variables
   int nmain=0;
   int i=0;
   int i_=0;
//--- initialization
   nmain=state.m_nmain;
//--- First,check for errors in the inputs
   if(!CAp::Assert(k>=0,__FUNCTION__+": K<0"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(c)>=nmain+1||k==0,__FUNCTION__+": Cols(C)<N+1"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(c)>=k,__FUNCTION__+": Rows(C)<K"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(ct)>=k,__FUNCTION__+": Length(CT)<K"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(c,k,nmain+1),__FUNCTION__+": C contains infinite or NaN values!"))
      return;
//--- Determine number of constraints,
//--- allocate space and copy
   state.m_cecnt=k;
//--- function call
   CApServ::RMatrixSetLengthAtLeast(state.m_ceoriginal,state.m_cecnt,nmain+1);
//--- function call
   CApServ::IVectorSetLengthAtLeast(state.m_ct,state.m_cecnt);
//--- calculation
   for(i=0;i<=k-1;i++)
     {
      state.m_ct[i]=ct[i];
      for(i_=0;i_<=nmain;i_++)
         state.m_ceoriginal[i].Set(i_,c[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| This function sets stopping conditions for the underlying        |
//| nonlinear CG optimizer. It controls overall accuracy of solution.|
//| These conditions should be strict enough in order for algorithm  |
//| to converge.                                                     |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     EpsG    -   >=0                                              |
//|                 The subroutine finishes its work if the condition|
//|                 |v|<EpsG is satisfied, where:                    |
//|                 * |.| means Euclidian norm                       |
//|                 * v - scaled gradient vector, v[i]=g[i]*s[i]     |
//|                 * g - gradient                                   |
//|                 * s - scaling coefficients set by                |
//|                       MinBLEICSetScale()                         |
//|     EpsF    -   >=0                                              |
//|                 The subroutine finishes its work if on k+1-th    |
//|                 iteration the condition |F(k+1)-F(k)| <=         |
//|                 <= EpsF*max{|F(k)|,|F(k+1)|,1} is satisfied.     |
//|     EpsX    -   >=0                                              |
//|                 The subroutine finishes its work if on k+1-th    |
//|                 iteration the condition |v|<=EpsX is fulfilled,  |
//|                 where:                                           |
//|                 * |.| means Euclidian norm                       |
//|                 * v - scaled step vector, v[i]=dx[i]/s[i]        |
//|                 * dx - ste pvector, dx=X(k+1)-X(k)               |
//|                 * s - scaling coefficients set by                |
//|                       MinBLEICSetScale()                         |
//| Passing EpsG=0, EpsF=0 and EpsX=0 (simultaneously) will lead to  |
//| automatic stopping criterion selection.                          |
//| These conditions are used to terminate inner iterations. However,|
//| you need to tune termination conditions for outer iterations too.|
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICSetInnerCond(CMinBLEICState &state,const double epsg,
                                            const double epsf,const double epsx)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsg),__FUNCTION__+": EpsG is not finite number"))
      return;
//--- check
   if(!CAp::Assert(epsg>=0.0,__FUNCTION__+": negative EpsG"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsf),__FUNCTION__+": EpsF is not finite number"))
      return;
//--- check
   if(!CAp::Assert(epsf>=0.0,__FUNCTION__+": negative EpsF"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsx),__FUNCTION__+": EpsX is not finite number"))
      return;
//--- check
   if(!CAp::Assert(epsx>=0.0,__FUNCTION__+": negative EpsX"))
      return;
//--- change values
   state.m_innerepsg=epsg;
   state.m_innerepsf=epsf;
   state.m_innerepsx=epsx;
  }
//+------------------------------------------------------------------+
//| This function sets stopping conditions for outer iteration of    |
//| BLEIC algo.                                                      |
//| These conditions control accuracy of constraint handling and     |
//| amount of infeasibility allowed in the solution.                 |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     EpsX    -   >0, stopping condition on outer iteration step   |
//|                     length                                       |
//|     EpsI    -   >0, stopping condition on infeasibility          |
//| Both EpsX and EpsI must be non-zero.                             |
//| MEANING OF EpsX                                                  |
//| EpsX is a stopping condition for outer iterations. Algorithm will|
//| stop when solution of the current modified subproblem will be    |
//| within EpsX (using 2-norm) of the previous solution.             |
//| MEANING OF EpsI                                                  |
//| EpsI controls feasibility properties - algorithm won't stop until|
//| all inequality constraints will be satisfied with error (distance|
//| from current point to the feasible area) at most EpsI.           |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICSetOuterCond(CMinBLEICState &state,const double epsx,
                                            const double epsi)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsx),__FUNCTION__+": EpsX is not finite number"))
      return;
//--- check
   if(!CAp::Assert(epsx>0.0,__FUNCTION__+": non-positive EpsX"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsi),__FUNCTION__+": EpsI is not finite number"))
      return;
//--- check
   if(!CAp::Assert((double)(epsi)>0.0,__FUNCTION__+": non-positive EpsI"))
      return;
//--- change values
   state.m_outerepsx=epsx;
   state.m_outerepsi=epsi;
  }
//+------------------------------------------------------------------+
//| This function sets scaling coefficients for BLEIC optimizer.     |
//| ALGLIB optimizers use scaling matrices to test stopping          |
//| conditions (step size and gradient are scaled before comparison  |
//| with tolerances). Scale of the I-th variable is a translation    |
//| invariant measure of:                                            |
//| a) "how large" the variable is                                   |
//| b) how large the step should be to make significant changes in   |
//| the function                                                     |
//| Scaling is also used by finite difference variant of the         |
//| optimizer - step along I-th axis is equal to DiffStep*S[I].      |
//| In most optimizers (and in the BLEIC too) scaling is NOT a form  |
//| of preconditioning. It just affects stopping conditions. You     |
//| should set preconditioner by separate call to one of the         |
//| MinBLEICSetPrec...() functions.                                  |
//| There is a special preconditioning mode, however, which uses     |
//| scaling coefficients to form diagonal preconditioning matrix.    |
//| You can turn this mode on, if you want. But you should understand|
//| that scaling is not the same thing as preconditioning - these are|
//| two different, although related forms of tuning solver.          |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure stores algorithm state                 |
//|     S       -   array[N], non-zero scaling coefficients          |
//|                 S[i] may be negative, sign doesn't matter.       |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICSetScale(CMinBLEICState &state,double &s[])
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(CAp::Len(s)>=state.m_nmain,__FUNCTION__+": Length(S)<N"))
      return;
   for(i=0;i<=state.m_nmain-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(s[i]),__FUNCTION__+": S contains infinite or NAN elements"))
         return;
      //--- check
      if(!CAp::Assert(s[i]!=0.0,__FUNCTION__+": S contains zero elements"))
         return;
      //--- change values
      state.m_soriginal[i]=MathAbs(s[i]);
     }
  }
//+------------------------------------------------------------------+
//| Modification of the preconditioner: preconditioning is turned    |
//| off.                                                             |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICSetPrecDefault(CMinBLEICState &state)
  {
//--- change value
   state.m_prectype=0;
  }
//+------------------------------------------------------------------+
//| Modification of the preconditioner: diagonal of approximate      |
//| Hessian is used.                                                 |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     D       -   diagonal of the approximate Hessian,             |
//|                 array[0..N-1], (if larger, only leading N        |
//|                 elements are used).                              |
//| NOTE 1: D[i] should be positive. Exception will be thrown        |
//|         otherwise.                                               |
//| NOTE 2: you should pass diagonal of approximate Hessian - NOT    |
//|         ITS INVERSE.                                             |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICSetPrecDiag(CMinBLEICState &state,double &d[])
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(CAp::Len(d)>=state.m_nmain,__FUNCTION__+": D is too short"))
      return;
   for(i=0;i<=state.m_nmain-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(d[i]),__FUNCTION__+": D contains infinite or NAN elements"))
         return;
      //--- check
      if(!CAp::Assert((double)(d[i])>0.0,__FUNCTION__+": D contains non-positive elements"))
         return;
     }
//--- function call
   CApServ::RVectorSetLengthAtLeast(state.m_diaghoriginal,state.m_nmain);
   state.m_prectype=2;
//--- copy
   for(i=0;i<=state.m_nmain-1;i++)
      state.m_diaghoriginal[i]=d[i];
  }
//+------------------------------------------------------------------+
//| Modification of the preconditioner: scale-based diagonal         |
//| preconditioning.                                                 |
//| This preconditioning mode can be useful when you don't have      |
//| approximate diagonal of Hessian, but you know that your variables|
//| are badly scaled (for example, one variable is in [1,10], and    |
//| another in [1000,100000]), and most part of the ill-conditioning |
//| comes from different scales of vars.                             |
//| In this case simple scale-based preconditioner, with H[i] =      |
//| = 1/(s[i]^2), can greatly improve convergence.                   |
//| IMPRTANT: you should set scale of your variables with            |
//| MinBLEICSetScale() call (before or after MinBLEICSetPrecScale()  |
//| call). Without knowledge of the scale of your variables          |
//| scale-based preconditioner will be just unit matrix.             |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICSetPrecScale(CMinBLEICState &state)
  {
//--- change value
   state.m_prectype=3;
  }
//+------------------------------------------------------------------+
//| This function allows to stop algorithm after specified number of |
//| inner iterations.                                                |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     MaxIts  -   maximum number of inner iterations.              |
//|                 If MaxIts=0, the number of iterations is         |
//|                 unlimited.                                       |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICSetMaxIts(CMinBLEICState &state,
                                         const int maxits)
  {
//--- check
   if(!CAp::Assert(maxits>=0,__FUNCTION__+": negative MaxIts!"))
      return;
//--- change value
   state.m_maxits=maxits;
  }
//+------------------------------------------------------------------+
//| This function turns on/off reporting.                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     NeedXRep-   whether iteration reports are needed or not      |
//| If NeedXRep is True, algorithm will call rep() callback function |
//| if it is provided to MinBLEICOptimize().                         |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICSetXRep(CMinBLEICState &state,
                                       const bool needxrep)
  {
//--- change value
   state.m_xrep=needxrep;
  }
//+------------------------------------------------------------------+
//| This function sets maximum step length                           |
//| IMPORTANT: this feature is hard to combine with preconditioning. |
//| You can't set upper limit on step length, when you solve         |
//| optimization problem with linear (non-boundary) constraints AND  |
//| preconditioner turned on.                                        |
//| When non-boundary constraints are present, you have to either a) |
//| use preconditioner, or b) use upper limit on step length. YOU    |
//| CAN'T USE BOTH! In this case algorithm will terminate with       |
//| appropriate error code.                                          |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     StpMax  -   maximum step length, >=0. Set StpMax to 0.0, if  |
//|                 you don't want to limit step length.             |
//| Use this subroutine when you optimize target function which      |
//| contains exp() or other fast growing functions, and optimization |
//| algorithm makes too large steps which lead to overflow. This     |
//| function allows us to reject steps that are too large (and       |
//| therefore expose us to the possible overflow) without actually   |
//| calculating function value at the x+stp*d.                       |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICSetStpMax(CMinBLEICState &state,
                                         const double stpmax)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(stpmax),__FUNCTION__+": StpMax is not finite!"))
      return;
//--- check
   if(!CAp::Assert(stpmax>=0.0,__FUNCTION__+": StpMax<0!"))
      return;
//--- change value
   state.m_stpmax=stpmax;
  }
//+------------------------------------------------------------------+
//| BLEIC results                                                    |
//| INPUT PARAMETERS:                                                |
//|     State   -   algorithm state                                  |
//| OUTPUT PARAMETERS:                                               |
//|     X       -   array[0..N-1], solution                          |
//|     Rep     -   optimization report. You should check Rep.       |
//|                 TerminationType in order to distinguish          |
//|                 successful termination from unsuccessful one.    |
//|                 More information about fields of this structure  |
//|                 can be found in the comments on MinBLEICReport   |
//|                 datatype.                                        |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICResults(CMinBLEICState &state,double &x[],
                                       CMinBLEICReport &rep)
  {
//--- reset memory
   ArrayResizeAL(x,0);
//--- function call
   MinBLEICResultsBuf(state,x,rep);
  }
//+------------------------------------------------------------------+
//| BLEIC results                                                    |
//| Buffered implementation of MinBLEICResults() which uses          |
//| pre-allocated buffer to store X[]. If buffer size is too small,  |
//| it resizes buffer. It is intended to be used in the inner cycles |
//| of performance critical algorithms where array reallocation      |
//| penalty is too large to be ignored.                              |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICResultsBuf(CMinBLEICState &state,double &x[],
                                          CMinBLEICReport &rep)
  {
//--- create variables
   int i=0;
   int i_=0;
//--- check
   if(CAp::Len(x)<state.m_nmain)
      ArrayResizeAL(x,state.m_nmain);
//--- change values
   rep.m_inneriterationscount=state.m_repinneriterationscount;
   rep.m_outeriterationscount=state.m_repouteriterationscount;
   rep.m_nfev=state.m_repnfev;
   rep.m_terminationtype=state.m_repterminationtype;
//--- check
   if(state.m_repterminationtype>0)
     {
      for(i_=0;i_<=state.m_nmain-1;i_++)
         x[i_]=state.m_xend[i_];
     }
   else
     {
      for(i=0;i<=state.m_nmain-1;i++)
         x[i]=CInfOrNaN::NaN();
     }
//--- change values
   rep.m_debugeqerr=state.m_repdebugeqerr;
   rep.m_debugfs=state.m_repdebugfs;
   rep.m_debugff=state.m_repdebugff;
   rep.m_debugdx=state.m_repdebugdx;
  }
//+------------------------------------------------------------------+
//| This subroutine restarts algorithm from new point.               |
//| All optimization parameters (including constraints) are left     |
//| unchanged.                                                       |
//| This function allows to solve multiple optimization problems     |
//| (which must have same number of dimensions) without object       |
//| reallocation penalty.                                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure previously allocated with              |
//|                 MinBLEICCreate call.                             |
//|     X       -   new starting point.                              |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICRestartFrom(CMinBLEICState &state,double &x[])
  {
//--- create variables
   int n=0;
   int i_=0;
//--- initialization
   n=state.m_nmain;
//--- First,check for errors in the inputs
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- Set XC
   for(i_=0;i_<=n-1;i_++)
      state.m_xstart[i_]=x[i_];
//--- prepare RComm facilities
   ArrayResizeAL(state.m_rstate.ia,5);
   ArrayResizeAL(state.m_rstate.ba,1);
   ArrayResizeAL(state.m_rstate.ra,2);
   state.m_rstate.stage=-1;
//--- function call
   ClearRequestFields(state);
  }
//+------------------------------------------------------------------+
//| Clears request fileds (to be sure that we don't forget to clear  |
//| something)                                                       |
//+------------------------------------------------------------------+
static void CMinBLEIC::ClearRequestFields(CMinBLEICState &state)
  {
//--- change values
   state.m_needf=false;
   state.m_needfg=false;
   state.m_xupdated=false;
  }
//+------------------------------------------------------------------+
//| This functions "unscales" point, i.e. it makes transformation    |
//| from scaled variables to unscaled ones. Only leading NMain       |
//| variables are copied from XUnscaled to XScaled.                  |
//+------------------------------------------------------------------+
static void CMinBLEIC::UnscalePoint(CMinBLEICState &state,double &xscaled[],
                                    double &xunscaled[])
  {
//--- create variables
   int    i=0;
   double v=0;
//--- calculation
   for(i=0;i<=state.m_nmain-1;i++)
     {
      v=xscaled[i]*state.m_transforms[i];
      //--- check
      if(state.m_hasbndl[i])
        {
         //--- check
         if(v<state.m_bndloriginal[i])
            v=state.m_bndloriginal[i];
        }
      //--- check
      if(state.m_hasbndu[i])
        {
         //--- check
         if(v>state.m_bnduoriginal[i])
            v=state.m_bnduoriginal[i];
        }
      xunscaled[i]=v;
     }
  }
//+------------------------------------------------------------------+
//| This function:                                                   |
//| 1. makes projection of XScaled into equality constrained subspace|
//|    (X is modified in-place)                                      |
//| 2. stores residual from the projection into R                    |
//| 3. unscales projected XScaled and stores result into XUnscaled   |
//|    with additional enforcement                                   |
//| It calculates set of additional values which are used later for  |
//| modification of the target function F.                           |
//| INPUT PARAMETERS:                                                |
//|     State   -   optimizer state (we use its fields to get        |
//|                 information about constraints)                   |
//|     X       -   vector being projected                           |
//|     R       -   preallocated buffer, used to store residual from |
//|                 projection                                       |
//| OUTPUT PARAMETERS:                                               |
//|     X       -   projection of input X                            |
//|     R       -   residual                                         |
//|     RNorm   -   residual norm squared, used later to modify      |
//|                 target function                                  |
//+------------------------------------------------------------------+
static void CMinBLEIC::ProjectPointAndUnscale(CMinBLEICState &state,
                                              double &xscaled[],double &xunscaled[],
                                              double &rscaled[],double &rnorm2)
  {
//--- create variables
   double v=0;
   int    i=0;
   int    nmain=0;
   int    nslack=0;
   int    i_=0;
//--- initialization
   rnorm2=0;
   nmain=state.m_nmain;
   nslack=state.m_nslack;
//--- * subtract XE from XScaled
//--- * project XScaled
//--- * calculate norm of deviation from null space,store it in RNorm2
//--- * calculate residual from projection,store it in R
//--- * add XE to XScaled
//--- * unscale variables
   for(i_=0;i_<=nmain+nslack-1;i_++)
      xscaled[i_]=xscaled[i_]-state.m_xe[i_];
   rnorm2=0;
   for(i=0;i<=nmain+nslack-1;i++)
      rscaled[i]=0;
//--- change values
   for(i=0;i<=nmain+nslack-1;i++)
     {
      //--- check
      if(state.m_activeconstraints[i])
        {
         v=xscaled[i];
         xscaled[i]=0;
         rscaled[i]=rscaled[i]+v;
         rnorm2=rnorm2+CMath::Sqr(v);
        }
     }
//--- calculation
   for(i=0;i<=state.m_cecnt-1;i++)
     {
      v=0.0;
      for(i_=0;i_<=nmain+nslack-1;i_++)
         v+=xscaled[i_]*state.m_cecurrent[i][i_];
      for(i_=0;i_<=nmain+nslack-1;i_++)
         xscaled[i_]=xscaled[i_]-v*state.m_cecurrent[i][i_];
      for(i_=0;i_<=nmain+nslack-1;i_++)
         rscaled[i_]=rscaled[i_]+v*state.m_cecurrent[i][i_];
      rnorm2=rnorm2+CMath::Sqr(v);
     }
   for(i_=0;i_<=nmain+nslack-1;i_++)
      xscaled[i_]=xscaled[i_]+state.m_xe[i_];
//--- function call
   UnscalePoint(state,xscaled,xunscaled);
  }
//+------------------------------------------------------------------+
//| This function scales and copies NMain elements of GUnscaled into |
//| GScaled. Other NSlack components of GScaled are set to zero.     |
//+------------------------------------------------------------------+
static void CMinBLEIC::ScaleGradientAndExpand(CMinBLEICState &state,
                                              double &gunscaled[],
                                              double &gscaled[])
  {
//--- create a variable
   int i=0;
//--- change values
   for(i=0;i<=state.m_nmain-1;i++)
      gscaled[i]=gunscaled[i]*state.m_transforms[i];
   for(i=0;i<=state.m_nslack-1;i++)
      gscaled[state.m_nmain+i]=0;
  }
//+------------------------------------------------------------------+
//| This subroutine applies modifications to the target function     |
//| given by its value F and gradient G at the projected point X     |
//| which lies in the equality constrained subspace.                 |
//| Following modifications are applied:                             |
//| * modified barrier functions to handle inequality constraints    |
//|   (both F and G are modified)                                    |
//| * projection of gradient into equality constrained subspace      |
//|   (only G is modified)                                           |
//| * quadratic penalty for deviations from equality constrained     |
//|   subspace (both F and G are modified)                           |
//| It also calculates gradient norm (three different norms for three|
//| different types of gradient), feasibility and complementary      |
//| slackness errors.                                                |
//| INPUT PARAMETERS:                                                |
//|     State   -   optimizer state (we use its fields to get        |
//|                 information about constraints)                   |
//|     X       -   point (projected into equality constrained       |
//|                 subspace)                                        |
//|     R       -   residual from projection                         |
//|     RNorm2  -   residual norm squared                            |
//|     F       -   function value at X                              |
//|     G       -   function gradient at X                           |
//| OUTPUT PARAMETERS:                                               |
//|     F       -   modified function value at X                     |
//|     G       -   modified function gradient at X                  |
//|     GNorm   -   2-norm of unmodified G                           |
//|     MPGNorm -   2-norm of modified G                             |
//|     MBA     -   minimum argument of barrier functions.           |
//|                 If X is strictly feasible, it is greater than    |
//|                 zero.                                            |
//|                 If X lies on a boundary, it is zero.             |
//|                 It is negative for infeasible X.                 |
//|     FIErr   -   2-norm of feasibility error with respect to      |
//|                 inequality/bound constraints                     |
//|     CSErr   -   2-norm of complementarity slackness error        |
//+------------------------------------------------------------------+
static void CMinBLEIC::ModifyTargetFunction(CMinBLEICState &state,double &x[],
                                            double &r[],const double rnorm2,
                                            double &f,double &g[],
                                            double &gnorm,double &mpgnorm)
  {
//--- create variables
   double v=0;
   int    i=0;
   int    nmain=0;
   int    nslack=0;
   bool   hasconstraints;
   int    i_=0;
//--- initialization
   gnorm=0;
   mpgnorm=0;
   nmain=state.m_nmain;
   nslack=state.m_nslack;
   hasconstraints=false;
//--- GNorm
   v=0.0;
   for(i_=0;i_<=nmain+nslack-1;i_++)
      v+=g[i_]*g[i_];
   gnorm=MathSqrt(v);
//--- Process equality constraints:
//--- * modify F to handle penalty term for equality constraints
//--- * project gradient on null space of equality constraints
//--- * add penalty term for equality constraints to gradient
   f=f+rnorm2;
   for(i=0;i<=nmain+nslack-1;i++)
     {
      //--- check
      if(state.m_activeconstraints[i])
         g[i]=0;
     }
   for(i=0;i<=state.m_cecnt-1;i++)
     {
      v=0.0;
      //--- change values
      for(i_=0;i_<=nmain+nslack-1;i_++)
         v+=g[i_]*state.m_cecurrent[i][i_];
      for(i_=0;i_<=nmain+nslack-1;i_++)
         g[i_]=g[i_]-v*state.m_cecurrent[i][i_];
     }
   for(i_=0;i_<=nmain+nslack-1;i_++)
      g[i_]=g[i_]+2*r[i_];
//--- MPGNorm
   v=0.0;
   for(i_=0;i_<=nmain+nslack-1;i_++)
      v+=g[i_]*g[i_];
   mpgnorm=MathSqrt(v);
  }
//+------------------------------------------------------------------+
//| This function makes additional check for constraints which can be|
//| activated.                                                       |
//| We try activate constraints one by one, but it is possible that  |
//| several constraints should be activated during one iteration. It |
//| this case only one of them (probably last) will be activated.    |
//| This function will fix it - it will pass through constraints and |
//| activate those which are at the boundary or beyond it.           |
//| It will return True, if at least one constraint was activated by |
//| this function.                                                   |
//+------------------------------------------------------------------+
static bool CMinBLEIC::AdditionalCheckForConstraints(CMinBLEICState &state,
                                                     double &x[])
  {
//--- create variables
   bool result;
   int  i=0;
   int  nmain=0;
   int  nslack=0;
//--- initialization
   result=false;
   nmain=state.m_nmain;
   nslack=state.m_nslack;
//--- calculation
   for(i=0;i<=nmain-1;i++)
     {
      //--- check
      if(!state.m_activeconstraints[i])
        {
         //--- check
         if(state.m_hasbndl[i])
           {
            //--- check
            if(x[i]<=state.m_bndleffective[i])
              {
               state.m_activeconstraints[i]=true;
               state.m_constrainedvalues[i]=state.m_bndleffective[i];
               result=true;
              }
           }
         //--- check
         if(state.m_hasbndu[i])
           {
            //--- check
            if(x[i]>=state.m_bndueffective[i])
              {
               state.m_activeconstraints[i]=true;
               state.m_constrainedvalues[i]=state.m_bndueffective[i];
               result=true;
              }
           }
        }
     }
   for(i=0;i<=nslack-1;i++)
     {
      //--- check
      if(!state.m_activeconstraints[nmain+i])
        {
         //--- check
         if(x[nmain+i]<=0.0)
           {
            state.m_activeconstraints[nmain+i]=true;
            state.m_constrainedvalues[nmain+i]=0;
            result=true;
           }
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| This function rebuilds CECurrent and XE according to current set |
//| of active bound constraints.                                     |
//+------------------------------------------------------------------+
static void CMinBLEIC::RebuildCEXE(CMinBLEICState &state)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    nmain=0;
   int    nslack=0;
   double v=0;
   int    i_=0;
//--- initialization
   nmain=state.m_nmain;
   nslack=state.m_nslack;
//--- function call
   CAblas::RMatrixCopy(state.m_cecnt,nmain+nslack+1,state.m_ceeffective,0,0,state.m_cecurrent,0,0);
//--- calculation
   for(i=0;i<=state.m_cecnt-1;i++)
     {
      //--- "Subtract" active bound constraints from I-th linear constraint
      for(j=0;j<=nmain+nslack-1;j++)
        {
         //--- check
         if(state.m_activeconstraints[j])
           {
            state.m_cecurrent[i].Set(nmain+nslack,state.m_cecurrent[i][nmain+nslack]-state.m_cecurrent[i][j]*state.m_constrainedvalues[j]);
            state.m_cecurrent[i].Set(j,0.0);
           }
        }
      //--- Reorthogonalize I-th constraint with respect to previous ones
      //--- NOTE: we also update right part,which is CECurrent[...,NMain+NSlack].
      for(k=0;k<=i-1;k++)
        {
         v=0.0;
         for(i_=0;i_<=nmain+nslack-1;i_++)
            v+=state.m_cecurrent[k][i_]*state.m_cecurrent[i][i_];
         for(i_=0;i_<=nmain+nslack;i_++)
            state.m_cecurrent[i].Set(i_,state.m_cecurrent[i][i_]-v*state.m_cecurrent[k][i_]);
        }
      //--- Calculate norm of I-th row of CECurrent. Fill by zeros,if it is
      //--- too small. Normalize otherwise.
      //--- NOTE: we also scale last column of CECurrent (right part)
      v=0.0;
      for(i_=0;i_<=nmain+nslack-1;i_++)
         v+=state.m_cecurrent[i][i_]*state.m_cecurrent[i][i_];
      v=MathSqrt(v);
      //--- check
      if(v>10000*CMath::m_machineepsilon)
        {
         v=1/v;
         for(i_=0;i_<=nmain+nslack;i_++)
            state.m_cecurrent[i].Set(i_,v*state.m_cecurrent[i][i_]);
        }
      else
        {
         for(j=0;j<=nmain+nslack;j++)
            state.m_cecurrent[i].Set(j,0);
        }
     }
//--- change values
   for(j=0;j<=nmain+nslack-1;j++)
      state.m_xe[j]=0;
   for(i=0;i<=nmain+nslack-1;i++)
     {
      //--- check
      if(state.m_activeconstraints[i])
         state.m_xe[i]=state.m_xe[i]+state.m_constrainedvalues[i];
     }
//--- change values
   for(i=0;i<=state.m_cecnt-1;i++)
     {
      v=state.m_cecurrent[i][nmain+nslack];
      for(i_=0;i_<=nmain+nslack-1;i_++)
         state.m_xe[i_]=state.m_xe[i_]+v*state.m_cecurrent[i][i_];
     }
  }
//+------------------------------------------------------------------+
//| This function projects gradient onto equality constrained        |
//| subspace                                                         |
//+------------------------------------------------------------------+
static void CMinBLEIC::MakeGradientProjection(CMinBLEICState &state,
                                              double &pg[])
  {
//--- create variables
   int    i=0;
   int    nmain=0;
   int    nslack=0;
   double v=0;
   int    i_=0;
//--- initialization
   nmain=state.m_nmain;
   nslack=state.m_nslack;
   for(i=0;i<=nmain+nslack-1;i++)
     {
      //--- check
      if(state.m_activeconstraints[i])
         pg[i]=0;
     }
//--- calculation
   for(i=0;i<=state.m_cecnt-1;i++)
     {
      v=0.0;
      for(i_=0;i_<=nmain+nslack-1;i_++)
         v+=pg[i_]*state.m_cecurrent[i][i_];
      for(i_=0;i_<=nmain+nslack-1;i_++)
         pg[i_]=pg[i_]-v*state.m_cecurrent[i][i_];
     }
  }
//+------------------------------------------------------------------+
//| This function prepares equality constrained subproblem:          |
//| 1. X is used to activate constraints (if there are ones which are|
//|    still inactive, but should be activated).                     |
//| 2. constraints matrix CEOrt is copied to CECurrent and modified  | 
//|    according to the list of active bound constraints             |
//|    (corresponding elements are filled by zeros and               |
//|    reorthogonalized).                                            |
//| 3. XE - least squares solution of equality constraints - is      |
//|    recalculated                                                  |
//| 4. X is copied to PX and projected onto equality constrained     |
//|    subspace                                                      |
//| 5. inactive constraints are checked against PX - if there is at  |
//|    least one which should be activated, we activate it and move  |
//|    back to (2)                                                   |
//| 6. as result, PX is feasible with respect to bound constraints - |
//|    step (5) guarantees it. But PX can be infeasible with respect |
//|    to equality ones, because step (2) is done without checks for |
//|    consistency. As the final step, we check that PX is feasible. |
//|    If not, we return False. True is returned otherwise.          |
//| If this algorithm returned True, then:                           |
//| * X is not changed                                               |
//| * PX contains projection of X onto constrained subspace          |
//| * G is not changed                                               |
//| * PG contains projection of G onto constrained subspace          |
//| * PX is feasible with respect to all constraints                 |
//| * all constraints which are active at PX, are activated          |
//+------------------------------------------------------------------+
static bool CMinBLEIC::PrepareConstraintMatrix(CMinBLEICState &state,double &x[],
                                               double &g[],double &px[],double &pg[])
  {
//--- create variables
   int    i=0;
   int    nmain=0;
   int    nslack=0;
   double v=0;
   double ferr=0;
   int    i_=0;
//--- initialization
   nmain=state.m_nmain;
   nslack=state.m_nslack;
//--- Step 1
   AdditionalCheckForConstraints(state,x);
//--- Steps 2-5
   do
     {
      //--- Steps 2-3
      RebuildCEXE(state);
      //--- Step 4
      //--- Calculate PX,PG
      for(i_=0;i_<=nmain+nslack-1;i_++)
         px[i_]=x[i_];
      for(i_=0;i_<=nmain+nslack-1;i_++)
         px[i_]=px[i_]-state.m_xe[i_];
      for(i_=0;i_<=nmain+nslack-1;i_++)
         pg[i_]=g[i_];
      for(i=0;i<=nmain+nslack-1;i++)
        {
         //--- check
         if(state.m_activeconstraints[i])
           {
            px[i]=0;
            pg[i]=0;
           }
        }
      //--- calculation
      for(i=0;i<=state.m_cecnt-1;i++)
        {
         //--- change values
         v=0.0;
         for(i_=0;i_<=nmain+nslack-1;i_++)
            v+=px[i_]*state.m_cecurrent[i][i_];
         for(i_=0;i_<=nmain+nslack-1;i_++)
            px[i_]=px[i_]-v*state.m_cecurrent[i][i_];
         //--- change values
         v=0.0;
         for(i_=0;i_<=nmain+nslack-1;i_++)
            v+=pg[i_]*state.m_cecurrent[i][i_];
         for(i_=0;i_<=nmain+nslack-1;i_++)
            pg[i_]=pg[i_]-v*state.m_cecurrent[i][i_];
        }
      for(i_=0;i_<=nmain+nslack-1;i_++)
         px[i_]=px[i_]+state.m_xe[i_];
      //--- Step 5 (loop condition below)
     }
   while(AdditionalCheckForConstraints(state,px));
//--- Step 6
   ferr=0;
   for(i=0;i<=state.m_cecnt-1;i++)
     {
      v=0.0;
      for(i_=0;i_<=nmain+nslack-1;i_++)
         v+=px[i_]*state.m_ceeffective[i][i_];
      //--- change values
      v=v-state.m_ceeffective[i][nmain+nslack];
      ferr=MathMax(ferr,MathAbs(v));
     }
//--- check
   if(ferr<=state.m_outerepsi)
      return(true);
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Internal initialization subroutine                               |
//+------------------------------------------------------------------+
static void CMinBLEIC::MinBLEICInitInternal(const int n,double &x[],
                                            const double diffstep,
                                            CMinBLEICState &state)
  {
//--- create a variable
   int i=0;
//--- create matrix
   CMatrixDouble c;
//--- create array
   int ct[];
//--- initialization
   state.m_nmain=n;
   state.m_optdim=0;
   state.m_diffstep=diffstep;
//--- allocation
   ArrayResizeAL(state.m_bndloriginal,n);
   ArrayResizeAL(state.m_bndleffective,n);
   ArrayResizeAL(state.m_hasbndl,n);
   ArrayResizeAL(state.m_bnduoriginal,n);
   ArrayResizeAL(state.m_bndueffective,n);
   ArrayResizeAL(state.m_hasbndu,n);
   ArrayResizeAL(state.m_xstart,n);
   ArrayResizeAL(state.m_soriginal,n);
   ArrayResizeAL(state.m_x,n);
   ArrayResizeAL(state.m_g,n);
   for(i=0;i<=n-1;i++)
     {
      state.m_bndloriginal[i]=CInfOrNaN::NegativeInfinity();
      state.m_hasbndl[i]=false;
      state.m_bnduoriginal[i]=CInfOrNaN::PositiveInfinity();
      state.m_hasbndu[i]=false;
      state.m_soriginal[i]=1.0;
     }
//--- function call
   MinBLEICSetLC(state,c,ct,0);
//--- function call
   MinBLEICSetInnerCond(state,0.0,0.0,0.0);
//--- function call
   MinBLEICSetOuterCond(state,1.0E-6,1.0E-6);
//--- function call
   MinBLEICSetMaxIts(state,0);
//--- function call
   MinBLEICSetXRep(state,false);
//--- function call
   MinBLEICSetStpMax(state,0.0);
//--- function call
   MinBLEICSetPrecDefault(state);
//--- function call
   MinBLEICRestartFrom(state,x);
  }
//+------------------------------------------------------------------+
//| NOTES:                                                           |
//| 1. This function has two different implementations: one which    |
//|    uses exact (analytical) user-supplied gradient, and one which |
//|    uses function value only and numerically differentiates       |
//|    function in order to obtain gradient.                         |
//|    Depending on the specific function used to create optimizer   |
//|    object (either MinBLEICCreate() for analytical gradient or    |
//|    MinBLEICCreateF() for numerical differentiation) you should   |
//|    choose appropriate variant of MinBLEICOptimize() - one which  |
//|    accepts function AND gradient or one which accepts function   |
//|    ONLY.                                                         |
//|    Be careful to choose variant of MinBLEICOptimize() which      |
//|    corresponds to your optimization scheme! Table below lists    |
//|    different combinations of callback (function/gradient) passed |
//|    to MinBLEICOptimize() and specific function used to create    |
//|    optimizer.                                                    |
//|                      |         USER PASSED TO MinBLEICOptimize() |
//|    CREATED WITH      |  function only   |  function and gradient |
//|    ------------------------------------------------------------  |
//|    MinBLEICCreateF() |     work                FAIL              |
//|    MinBLEICCreate()  |     FAIL                work              |
//|    Here "FAIL" denotes inappropriate combinations of optimizer   |
//|    creation function and MinBLEICOptimize() version. Attemps to  |
//|    use such combination (for example, to create optimizer with   |
//|    MinBLEICCreateF() and to pass gradient information to         |
//|    MinCGOptimize()) will lead to exception being thrown. Either  |
//|    you did not pass gradient when it WAS needed or you passed    |
//|    gradient when it was NOT needed.                              |
//+------------------------------------------------------------------+
static bool CMinBLEIC::MinBLEICIteration(CMinBLEICState &state)
  {
//--- create variables
   int    nmain=0;
   int    nslack=0;
   int    m=0;
   int    i=0;
   int    j=0;
   double v=0;
   double vv=0;
   bool   b;
   int    i_=0;
//--- This code initializes locals by:
//--- * random values determined during code
//---   generation - on first subroutine call
//--- * values from previous call - on subsequent calls
   if(state.m_rstate.stage>=0)
     {
      //--- initialization
      nmain=state.m_rstate.ia[0];
      nslack=state.m_rstate.ia[1];
      m=state.m_rstate.ia[2];
      i=state.m_rstate.ia[3];
      j=state.m_rstate.ia[4];
      b=state.m_rstate.ba[0];
      v=state.m_rstate.ra[0];
      vv=state.m_rstate.ra[1];
     }
   else
     {
      //--- initialization
      nmain=-983;
      nslack=-989;
      m=-834;
      i=900;
      j=-287;
      b=false;
      v=214;
      vv=-338;
     }
//--- check
   if(state.m_rstate.stage==0)
     {
      //--- change value
      state.m_needfg=false;
      //--- function call, return result
      return(Func_lbl_14(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- check
   if(state.m_rstate.stage==1)
     {
      //--- change value
      state.m_needf=false;
      //--- function call, return result
      return(Func_lbl_14(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- check
   if(state.m_rstate.stage==2)
     {
      //--- change value
      state.m_needfg=false;
      //--- function call, return result
      return(Func_lbl_22(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- check
   if(state.m_rstate.stage==3)
     {
      //--- change values
      state.m_fbase=state.m_f;
      i=0;
      //--- function call, return result
      return(Func_lbl_23(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- check
   if(state.m_rstate.stage==4)
     {
      //--- change values
      state.m_fm2=state.m_f;
      state.m_x[i]=v-0.5*state.m_diffstep*state.m_soriginal[i];
      state.m_rstate.stage=5;
      //--- Saving state
      Func_lbl_rcomm(state,nmain,nslack,m,i,j,b,v,vv);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==5)
     {
      //--- change values
      state.m_fm1=state.m_f;
      state.m_x[i]=v+0.5*state.m_diffstep*state.m_soriginal[i];
      state.m_rstate.stage=6;
      //--- Saving state
      Func_lbl_rcomm(state,nmain,nslack,m,i,j,b,v,vv);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==6)
     {
      //--- change values
      state.m_fp1=state.m_f;
      state.m_x[i]=v+state.m_diffstep*state.m_soriginal[i];
      state.m_rstate.stage=7;
      //--- Saving state
      Func_lbl_rcomm(state,nmain,nslack,m,i,j,b,v,vv);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==7)
     {
      //--- change values
      state.m_fp2=state.m_f;
      state.m_g[i]=(8*(state.m_fp1-state.m_fm1)-(state.m_fp2-state.m_fm2))/(6*state.m_diffstep*state.m_soriginal[i]);
      state.m_x[i]=v;
      i=i+1;
      //--- function call, return result
      return(Func_lbl_23(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- check
   if(state.m_rstate.stage==8)
     {
      //--- change values
      state.m_fm1=state.m_f;
      state.m_xp1=MathMin(v+state.m_diffstep*state.m_soriginal[i],state.m_bnduoriginal[i]);
      state.m_x[i]=state.m_xp1;
      state.m_rstate.stage=9;
      //--- Saving state
      Func_lbl_rcomm(state,nmain,nslack,m,i,j,b,v,vv);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==9)
     {
      //--- change values
      state.m_fp1=state.m_f;
      state.m_g[i]=(state.m_fp1-state.m_fm1)/(state.m_xp1-state.m_xm1);
      state.m_x[i]=v;
      i=i+1;
      //--- function call, return result
      return(Func_lbl_23(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- check
   if(state.m_rstate.stage==10)
     {
      //--- change value
      state.m_xupdated=false;
      //--- function call, return result
      return(Func_lbl_17(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- check
   if(state.m_rstate.stage==11)
     {
      //--- change value
      state.m_needfg=false;
      //--- function call, return result
      return(Func_lbl_31(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- check
   if(state.m_rstate.stage==12)
     {
      //--- change value
      state.m_needf=false;
      //--- function call, return result
      return(Func_lbl_31(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- Routine body
//--- Prepare:
//--- * calculate number of slack variables
//--- * initialize locals
//--- * initialize debug fields
//--- * make quick check
   nmain=state.m_nmain;
   nslack=0;
   for(i=0;i<=state.m_cecnt-1;i++)
     {
      //--- check
      if(state.m_ct[i]!=0)
         nslack=nslack+1;
     }
//--- change values
   state.m_nslack=nslack;
   state.m_repterminationtype=0;
   state.m_repinneriterationscount=0;
   state.m_repouteriterationscount=0;
   state.m_repnfev=0;
   state.m_repdebugeqerr=0.0;
   state.m_repdebugfs=CInfOrNaN::NaN();
   state.m_repdebugff=CInfOrNaN::NaN();
   state.m_repdebugdx=CInfOrNaN::NaN();
//--- check
   if(state.m_stpmax!=0.0 && state.m_prectype!=0)
     {
      state.m_repterminationtype=-10;
      //--- return result
      return(false);
     }
//--- allocate
   CApServ::RVectorSetLengthAtLeast(state.m_r,nmain+nslack);
   CApServ::RVectorSetLengthAtLeast(state.m_diagh,nmain+nslack);
   CApServ::RVectorSetLengthAtLeast(state.m_tmp0,nmain+nslack);
   CApServ::RVectorSetLengthAtLeast(state.m_tmp1,nmain+nslack);
   CApServ::RVectorSetLengthAtLeast(state.m_tmp2,nmain+nslack);
   CApServ::RMatrixSetLengthAtLeast(state.m_cecurrent,state.m_cecnt,nmain+nslack+1);
   CApServ::BVectorSetLengthAtLeast(state.m_activeconstraints,nmain+nslack);
   CApServ::RVectorSetLengthAtLeast(state.m_constrainedvalues,nmain+nslack);
   CApServ::RVectorSetLengthAtLeast(state.m_lastg,nmain+nslack);
   CApServ::RVectorSetLengthAtLeast(state.m_xe,nmain+nslack);
   CApServ::RVectorSetLengthAtLeast(state.m_xcur,nmain+nslack);
   CApServ::RVectorSetLengthAtLeast(state.m_xprev,nmain+nslack);
   CApServ::RVectorSetLengthAtLeast(state.m_xend,nmain);
//--- Create/restart optimizer.
//--- State.OptDim is used to determine current state of optimizer.
   if(state.m_optdim!=nmain+nslack)
     {
      for(i=0;i<=nmain+nslack-1;i++)
         state.m_tmp1[i]=0.0;
      //--- function call
      CMinCG::MinCGCreate(nmain+nslack,state.m_tmp1,state.m_cgstate);
      state.m_optdim=nmain+nslack;
     }
//--- Prepare transformation.
//--- MinBLEIC's handling of preconditioner matrix is somewhat unusual -
//--- instead of incorporating it into algorithm and making implicit
//--- scaling (as most optimizers do) BLEIC optimizer uses explicit
//--- scaling - it solves problem in the scaled parameters space S,
//--- making transition between scaled (S) and unscaled (X) variables
//--- every time we ask for function value.
//--- Following fields are calculated here:
//--- * TransformS         X[i]=TransformS[i]*S[i],array[NMain]
//--- * SEffective         "effective" scale of the variables after
//---                      transformation,array[NMain+NSlack]
   CApServ::RVectorSetLengthAtLeast(state.m_transforms,nmain);
   for(i=0;i<=nmain-1;i++)
     {
      //--- check
      if(state.m_prectype==2)
        {
         state.m_transforms[i]=1/MathSqrt(state.m_diaghoriginal[i]);
         continue;
        }
      //--- check
      if(state.m_prectype==3)
        {
         state.m_transforms[i]=state.m_soriginal[i];
         continue;
        }
      state.m_transforms[i]=1;
     }
//--- function call
   CApServ::RVectorSetLengthAtLeast(state.m_seffective,nmain+nslack);
   for(i=0;i<=nmain-1;i++)
      state.m_seffective[i]=state.m_soriginal[i]/state.m_transforms[i];
   for(i=0;i<=nslack-1;i++)
      state.m_seffective[nmain+i]=1;
//--- function call
   CMinCG::MinCGSetScale(state.m_cgstate,state.m_seffective);
//--- Pre-process constraints
//--- * check consistency of bound constraints
//--- * add slack vars,convert problem to the bound/equality
//---   constrained one
//--- We calculate here:
//--- * BndLEffective - lower bounds after transformation of variables (see above)
//--- * BndUEffective - upper bounds after transformation of variables (see above)
//--- * CEEffective - matrix of equality constraints for transformed variables
   for(i=0;i<=nmain-1;i++)
     {
      //--- check
      if(state.m_hasbndl[i])
         state.m_bndleffective[i]=state.m_bndloriginal[i]/state.m_transforms[i];
      //--- check
      if(state.m_hasbndu[i])
         state.m_bndueffective[i]=state.m_bnduoriginal[i]/state.m_transforms[i];
     }
   for(i=0;i<=nmain-1;i++)
     {
      //--- check
      if(state.m_hasbndl[i] && state.m_hasbndu[i])
        {
         //--- check
         if(state.m_bndleffective[i]>state.m_bndueffective[i])
           {
            state.m_repterminationtype=-3;
            //--- return result
            return(false);
           }
        }
     }
//--- function call
   CApServ::RMatrixSetLengthAtLeast(state.m_ceeffective,state.m_cecnt,nmain+nslack+1);
//--- change value
   m=0;
//--- calculation
   for(i=0;i<=state.m_cecnt-1;i++)
     {
      //--- NOTE: when we add slack variable,we use V=max(abs(CE[i,...])) as
      //--- coefficient before it in order to make linear equations better
      //--- conditioned.
      v=0;
      for(j=0;j<=nmain-1;j++)
        {
         state.m_ceeffective[i].Set(j,state.m_ceoriginal[i][j]*state.m_transforms[j]);
         v=MathMax(v,MathAbs(state.m_ceeffective[i][j]));
        }
      //--- check
      if(v==0.0)
         v=1;
      for(j=0;j<=nslack-1;j++)
         state.m_ceeffective[i].Set(nmain+j,0.0);
      state.m_ceeffective[i].Set(nmain+nslack,state.m_ceoriginal[i][nmain]);
      //--- check
      if(state.m_ct[i]<0)
        {
         state.m_ceeffective[i].Set(nmain+m,v);
         m=m+1;
        }
      //--- check
      if(state.m_ct[i]>0)
        {
         state.m_ceeffective[i].Set(nmain+m,-v);
         m=m+1;
        }
     }
//--- Find feasible point.
//--- 0. Convert from unscaled values (as stored in XStart) to scaled
//---    ones
//--- 1. calculate values of slack variables such that starting
//---    point satisfies inequality constraints (after conversion to
//---    equality ones) as much as possible.
//--- 2. use PrepareConstraintMatrix() function,which forces X
//---    to be strictly feasible.
   for(i=0;i<=nmain-1;i++)
      state.m_tmp0[i]=state.m_xstart[i]/state.m_transforms[i];
//--- change value
   m=0;
//--- calculation
   for(i=0;i<=state.m_cecnt-1;i++)
     {
      v=0.0;
      for(i_=0;i_<=nmain-1;i_++)
         v+=state.m_ceeffective[i][i_]*state.m_tmp0[i_];
      //--- check
      if(state.m_ct[i]<0)
        {
         state.m_tmp0[nmain+m]=state.m_ceeffective[i][nmain+nslack]-v;
         m=m+1;
        }
      //--- check
      if(state.m_ct[i]>0)
        {
         state.m_tmp0[nmain+m]=v-state.m_ceeffective[i][nmain+nslack];
         m=m+1;
        }
     }
//--- change values
   for(i=0;i<=nmain+nslack-1;i++)
      state.m_tmp1[i]=0;
   for(i=0;i<=nmain+nslack-1;i++)
      state.m_activeconstraints[i]=false;
//--- function call
   b=PrepareConstraintMatrix(state,state.m_tmp0,state.m_tmp1,state.m_xcur,state.m_tmp2);
   state.m_repdebugeqerr=0.0;
//--- calculation
   for(i=0;i<=state.m_cecnt-1;i++)
     {
      v=0.0;
      for(i_=0;i_<=nmain+nslack-1;i_++)
         v+=state.m_ceeffective[i][i_]*state.m_xcur[i_];
      state.m_repdebugeqerr=state.m_repdebugeqerr+CMath::Sqr(v-state.m_ceeffective[i][nmain+nslack]);
     }
   state.m_repdebugeqerr=MathSqrt(state.m_repdebugeqerr);
//--- check
   if(!b)
     {
      state.m_repterminationtype=-3;
      //--- return result
      return(false);
     }
//--- Initialize RepDebugFS with function value at initial point
   UnscalePoint(state,state.m_xcur,state.m_x);
   ClearRequestFields(state);
//--- check
   if(state.m_diffstep!=0.0)
     {
      state.m_needf=true;
      state.m_rstate.stage=1;
      //--- Saving state
      Func_lbl_rcomm(state,nmain,nslack,m,i,j,b,v,vv);
      //--- return result
      return(true);
     }
//--- change values
   state.m_needfg=true;
   state.m_rstate.stage=0;
//--- Saving state
   Func_lbl_rcomm(state,nmain,nslack,m,i,j,b,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinBLEICIteration. Is a product to get rid|
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static void CMinBLEIC::Func_lbl_rcomm(CMinBLEICState &state,int nmain,
                                      int nslack,int m,int i,int j,
                                      bool b,double v,double vv)
  {
//--- save
   state.m_rstate.ia[0]=nmain;
   state.m_rstate.ia[1]=nslack;
   state.m_rstate.ia[2]=m;
   state.m_rstate.ia[3]=i;
   state.m_rstate.ia[4]=j;
   state.m_rstate.ba[0]=b;
   state.m_rstate.ra[0]=v;
   state.m_rstate.ra[1]=vv;
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinBLEICIteration. Is a product to get rid|
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinBLEIC::Func_lbl_14(CMinBLEICState &state,int &nmain,
                                   int &nslack,int &m,int &i,int &j,
                                   bool &b,double &v,double &vv)
  {
//--- function call
   COptServ::TrimPrepare(state.m_f,state.m_trimthreshold);
   state.m_repnfev=state.m_repnfev+1;
   state.m_repdebugfs=state.m_f;
//--- Outer cycle
   state.m_itsleft=state.m_maxits;
//--- copy
   for(int i_=0;i_<=nmain+nslack-1;i_++)
      state.m_xprev[i_]=state.m_xcur[i_];
//--- function call, return result
   return(Func_lbl_15(state,nmain,nslack,m,i,j,b,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinBLEICIteration. Is a product to get rid|
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinBLEIC::Func_lbl_15(CMinBLEICState &state,int &nmain,
                                   int &nslack,int &m,int &i,int &j,
                                   bool &b,double &v,double &vv)
  {
//--- check
   if(!CAp::Assert(state.m_prectype==0 || state.m_stpmax==0.0,"MinBLEIC: internal error (-10)"))
      return(false);
//--- Inner cycle: CG with projections and penalty functions
   for(int i_=0;i_<=nmain+nslack-1;i_++)
      state.m_tmp0[i_]=state.m_xcur[i_];
   for(i=0;i<=nmain+nslack-1;i++)
     {
      state.m_tmp1[i]=0;
      state.m_activeconstraints[i]=false;
     }
//--- check
   if(!PrepareConstraintMatrix(state,state.m_tmp0,state.m_tmp1,state.m_xcur,state.m_tmp2))
     {
      state.m_repterminationtype=-3;
      //--- return result
      return(false);
     }
   for(i=0;i<=nmain+nslack-1;i++)
      state.m_activeconstraints[i]=false;
   RebuildCEXE(state);
//--- function call
   CMinCG::MinCGRestartFrom(state.m_cgstate,state.m_xcur);
//--- function call
   CMinCG::MinCGSetCond(state.m_cgstate,state.m_innerepsg,state.m_innerepsf,state.m_innerepsx,state.m_itsleft);
//--- function call
   CMinCG::MinCGSetXRep(state.m_cgstate,state.m_xrep);
//--- function call
   CMinCG::MinCGSetDRep(state.m_cgstate,true);
//--- function call
   CMinCG::MinCGSetStpMax(state.m_cgstate,state.m_stpmax);
//--- function call, return result
   return(Func_lbl_17(state,nmain,nslack,m,i,j,b,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinBLEICIteration. Is a product to get rid|
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinBLEIC::Func_lbl_16(CMinBLEICState &state,int &nmain,
                                   int &nslack,int &m,int &i,int &j,
                                   bool &b,double &v,double &vv)
  {
//--- We've stopped,fill debug information
   state.m_repdebugeqerr=0.0;
   for(i=0;i<=state.m_cecnt-1;i++)
     {
      //--- change value
      v=0.0;
      for(int i_=0;i_<=nmain+nslack-1;i_++)
         v+=state.m_ceeffective[i][i_]*state.m_xcur[i_];
      state.m_repdebugeqerr=state.m_repdebugeqerr+CMath::Sqr(v-state.m_ceeffective[i][nmain+nslack]);
     }
//--- change values
   state.m_repdebugeqerr=MathSqrt(state.m_repdebugeqerr);
   state.m_repdebugdx=0;
   for(i=0;i<=nmain-1;i++)
      state.m_repdebugdx=state.m_repdebugdx+CMath::Sqr(state.m_xcur[i]-state.m_xstart[i]);
   state.m_repdebugdx=MathSqrt(state.m_repdebugdx);
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinBLEICIteration. Is a product to get rid|
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinBLEIC::Func_lbl_17(CMinBLEICState &state,int &nmain,
                                   int &nslack,int &m,int &i,int &j,
                                   bool &b,double &v,double &vv)
  {
//--- check
   if(!CMinCG::MinCGIteration(state.m_cgstate))
      return(Func_lbl_18(state,nmain,nslack,m,i,j,b,v,vv));
//--- process different requests/reports of inner optimizer
   if(state.m_cgstate.m_algpowerup)
     {
      for(i=0;i<=nmain+nslack-1;i++)
         state.m_activeconstraints[i]=false;
      //--- cycle
      do
        {
         //--- function call
         RebuildCEXE(state);
         for(int i_=0;i_<=nmain+nslack-1;i_++)
            state.m_tmp1[i_]=state.m_cgstate.m_g[i_];
         //--- function call
         MakeGradientProjection(state,state.m_tmp1);
         b=false;
         for(i=0;i<=nmain-1;i++)
           {
            //--- check
            if(!state.m_activeconstraints[i])
              {
               //--- check
               if(state.m_hasbndl[i])
                 {
                  //--- check
                  if(state.m_cgstate.m_x[i]==state.m_bndleffective[i] && state.m_tmp1[i]>=0.0)
                    {
                     //--- change values
                     state.m_activeconstraints[i]=true;
                     state.m_constrainedvalues[i]=state.m_bndleffective[i];
                     b=true;
                    }
                 }
               //--- check
               if(state.m_hasbndu[i])
                 {
                  //--- check
                  if(state.m_cgstate.m_x[i]==state.m_bndueffective[i] && state.m_tmp1[i]<=0.0)
                    {
                     //--- change values
                     state.m_activeconstraints[i]=true;
                     state.m_constrainedvalues[i]=state.m_bndueffective[i];
                     b=true;
                    }
                 }
              }
           }
         for(i=0;i<=nslack-1;i++)
           {
            //--- check
            if(!state.m_activeconstraints[nmain+i])
              {
               //--- check
               if(state.m_cgstate.m_x[nmain+i]==0.0 && state.m_tmp1[nmain+i]>=0.0)
                 {
                  //--- change values
                  state.m_activeconstraints[nmain+i]=true;
                  state.m_constrainedvalues[nmain+i]=0;
                  b=true;
                 }
              }
           }
        }
      while(b);
      //--- copy
      for(int i_=0;i_<=nmain+nslack-1;i_++)
         state.m_cgstate.m_g[i_]=state.m_tmp1[i_];
      //--- function call, return result
      return(Func_lbl_17(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- check
   if(state.m_cgstate.m_lsstart)
     {
      //--- Beginning of the line search: set upper limit on step size
      //--- to prevent algo from leaving feasible area.
      state.m_variabletofreeze=-1;
      //--- check
      if((double)(state.m_cgstate.m_curstpmax)==0.0)
         state.m_cgstate.m_curstpmax=1.0E50;
      for(i=0;i<=nmain-1;i++)
        {
         //--- check
         if(state.m_hasbndl[i] && state.m_cgstate.m_d[i]<0.0)
           {
            //--- change values
            v=state.m_cgstate.m_curstpmax;
            vv=state.m_cgstate.m_x[i]-state.m_bndleffective[i];
            //--- check
            if(vv<0.0)
               vv=0;
            state.m_cgstate.m_curstpmax=CApServ::SafeMinPosRV(vv,-state.m_cgstate.m_d[i],state.m_cgstate.m_curstpmax);
            //--- check
            if(state.m_cgstate.m_curstpmax<v)
              {
               state.m_variabletofreeze=i;
               state.m_valuetofreeze=state.m_bndleffective[i];
              }
           }
         //--- check
         if(state.m_hasbndu[i] && state.m_cgstate.m_d[i]>0.0)
           {
            //--- change values
            v=state.m_cgstate.m_curstpmax;
            vv=state.m_bndueffective[i]-state.m_cgstate.m_x[i];
            //--- check
            if(vv<0.0)
               vv=0;
            state.m_cgstate.m_curstpmax=CApServ::SafeMinPosRV(vv,state.m_cgstate.m_d[i],state.m_cgstate.m_curstpmax);
            //--- check
            if(state.m_cgstate.m_curstpmax<v)
              {
               state.m_variabletofreeze=i;
               state.m_valuetofreeze=state.m_bndueffective[i];
              }
           }
        }
      for(i=0;i<=nslack-1;i++)
        {
         //--- check
         if(state.m_cgstate.m_d[nmain+i]<0.0)
           {
            //--- change values
            v=state.m_cgstate.m_curstpmax;
            vv=state.m_cgstate.m_x[nmain+i];
            //--- check
            if(vv<0.0)
               vv=0;
            state.m_cgstate.m_curstpmax=CApServ::SafeMinPosRV(vv,-state.m_cgstate.m_d[nmain+i],state.m_cgstate.m_curstpmax);
            //--- check
            if(state.m_cgstate.m_curstpmax<v)
              {
               state.m_variabletofreeze=nmain+i;
               state.m_valuetofreeze=0;
              }
           }
        }
      //--- check
      if(state.m_cgstate.m_curstpmax==0.0)
        {
         //--- change values
         state.m_activeconstraints[state.m_variabletofreeze]=true;
         state.m_constrainedvalues[state.m_variabletofreeze]=state.m_valuetofreeze;
         state.m_cgstate.m_x[state.m_variabletofreeze]=state.m_valuetofreeze;
         state.m_cgstate.m_terminationneeded=true;
        }
      //--- function call, return result
      return(Func_lbl_17(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- check
   if(state.m_cgstate.m_lsend)
     {
      //--- Line search just finished.
      //--- Maybe we should activate some constraints?
      b=state.m_cgstate.m_stp>=state.m_cgstate.m_curstpmax && state.m_variabletofreeze>=0;
      //--- check
      if(b)
        {
         state.m_activeconstraints[state.m_variabletofreeze]=true;
         state.m_constrainedvalues[state.m_variabletofreeze]=state.m_valuetofreeze;
        }
      //--- Additional activation of constraints
      b=b || AdditionalCheckForConstraints(state,state.m_cgstate.m_x);
      //--- If at least one constraint was activated we have to rebuild constraint matrices
      if(b)
        {
         //--- copy
         for(int i_=0;i_<=nmain+nslack-1;i_++)
            state.m_tmp0[i_]=state.m_cgstate.m_x[i_];
         for(int i_=0;i_<=nmain+nslack-1;i_++)
            state.m_tmp1[i_]=state.m_lastg[i_];
         //--- check
         if(!PrepareConstraintMatrix(state,state.m_tmp0,state.m_tmp1,state.m_cgstate.m_x,state.m_cgstate.m_g))
           {
            state.m_repterminationtype=-3;
            //--- return result
            return(false);
           }
         state.m_cgstate.m_innerresetneeded=true;
        }
      //--- function call, return result
      return(Func_lbl_17(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- check
   if(!state.m_cgstate.m_needfg)
      return(Func_lbl_19(state,nmain,nslack,m,i,j,b,v,vv));
//--- copy
   for(int i_=0;i_<=nmain+nslack-1;i_++)
      state.m_tmp1[i_]=state.m_cgstate.m_x[i_];
//--- function call
   ProjectPointAndUnscale(state,state.m_tmp1,state.m_x,state.m_r,vv);
//--- function call
   ClearRequestFields(state);
//--- check
   if(state.m_diffstep!=0.0)
     {
      //--- change values
      state.m_needf=true;
      state.m_rstate.stage=3;
      //--- Saving state
      Func_lbl_rcomm(state,nmain,nslack,m,i,j,b,v,vv);
      //--- return result
      return(true);
     }
//--- change values
   state.m_needfg=true;
   state.m_rstate.stage=2;
//--- Saving state
   Func_lbl_rcomm(state,nmain,nslack,m,i,j,b,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinBLEICIteration. Is a product to get rid|
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinBLEIC::Func_lbl_18(CMinBLEICState &state,int &nmain,
                                   int &nslack,int &m,int &i,int &j,
                                   bool &b,double &v,double &vv)
  {
//--- function call
   CMinCG::MinCGResults(state.m_cgstate,state.m_xcur,state.m_cgrep);
//--- function call
   UnscalePoint(state,state.m_xcur,state.m_xend);
//--- change values
   state.m_repinneriterationscount=state.m_repinneriterationscount+state.m_cgrep.m_iterationscount;
   state.m_repouteriterationscount=state.m_repouteriterationscount+1;
   state.m_repnfev=state.m_repnfev+state.m_cgrep.m_nfev;
//--- Update RepDebugFF with function value at current point
   UnscalePoint(state,state.m_xcur,state.m_x);
   ClearRequestFields(state);
//--- check
   if(state.m_diffstep!=0.0)
     {
      //--- change values
      state.m_needf=true;
      state.m_rstate.stage=12;
      //--- Saving state
      Func_lbl_rcomm(state,nmain,nslack,m,i,j,b,v,vv);
      //--- return result
      return(true);
     }
//--- change values
   state.m_needfg=true;
   state.m_rstate.stage=11;
//--- Saving state
   Func_lbl_rcomm(state,nmain,nslack,m,i,j,b,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinBLEICIteration. Is a product to get rid|
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinBLEIC::Func_lbl_19(CMinBLEICState &state,int &nmain,
                                   int &nslack,int &m,int &i,int &j,
                                   bool &b,double &v,double &vv)
  {
//--- check
   if(!state.m_cgstate.m_xupdated)
      return(Func_lbl_17(state,nmain,nslack,m,i,j,b,v,vv));
//--- Report
   UnscalePoint(state,state.m_cgstate.m_x,state.m_x);
   state.m_f=state.m_cgstate.m_f;
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=10;
//--- Saving state
   Func_lbl_rcomm(state,nmain,nslack,m,i,j,b,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinBLEICIteration. Is a product to get rid|
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinBLEIC::Func_lbl_22(CMinBLEICState &state,int &nmain,
                                   int &nslack,int &m,int &i,int &j,
                                   bool &b,double &v,double &vv)
  {
//--- check
   if(state.m_f<state.m_trimthreshold)
     {
      //--- normal processing
      state.m_cgstate.m_f=state.m_f;
      //--- function call
      ScaleGradientAndExpand(state,state.m_g,state.m_cgstate.m_g);
      //--- copy
      for(int i_=0;i_<=nmain+nslack-1;i_++)
         state.m_lastg[i_]=state.m_cgstate.m_g[i_];
      //--- function call
      ModifyTargetFunction(state,state.m_tmp1,state.m_r,vv,state.m_cgstate.m_f,state.m_cgstate.m_g,state.m_gnorm,state.m_mpgnorm);
     }
   else
     {
      //--- function value is too high,trim it
      state.m_cgstate.m_f=state.m_trimthreshold;
      for(i=0;i<=nmain+nslack-1;i++)
         state.m_cgstate.m_g[i]=0.0;
     }
//--- function call, return result
   return(Func_lbl_17(state,nmain,nslack,m,i,j,b,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinBLEICIteration. Is a product to get rid|
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinBLEIC::Func_lbl_23(CMinBLEICState &state,int &nmain,
                                   int &nslack,int &m,int &i,int &j,
                                   bool &b,double &v,double &vv)
  {
//--- check
   if(i>nmain-1)
     {
      //--- change values
      state.m_f=state.m_fbase;
      state.m_needf=false;
      //--- function call, return result
      return(Func_lbl_22(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- change values
   v=state.m_x[i];
   b=false;
//--- check
   if(state.m_hasbndl[i])
      b=b || v-state.m_diffstep*state.m_soriginal[i]<state.m_bndloriginal[i];
//--- check
   if(state.m_hasbndu[i])
      b=b || v+state.m_diffstep*state.m_soriginal[i]>state.m_bnduoriginal[i];
//--- check
   if(b)
     {
      //--- change values
      state.m_xm1=MathMax(v-state.m_diffstep*state.m_soriginal[i],state.m_bndloriginal[i]);
      state.m_x[i]=state.m_xm1;
      state.m_rstate.stage=8;
      //--- Saving state
      Func_lbl_rcomm(state,nmain,nslack,m,i,j,b,v,vv);
      //--- return result
      return(true);
     }
//--- change values
   state.m_x[i]=v-state.m_diffstep*state.m_soriginal[i];
   state.m_rstate.stage=4;
//--- Saving state
   Func_lbl_rcomm(state,nmain,nslack,m,i,j,b,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinBLEICIteration. Is a product to get rid|
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinBLEIC::Func_lbl_31(CMinBLEICState &state,int &nmain,
                                   int &nslack,int &m,int &i,int &j,
                                   bool &b,double &v,double &vv)
  {
//--- change values
   state.m_repnfev=state.m_repnfev+1;
   state.m_repdebugff=state.m_f;
//--- Check for stopping:
//--- * "normal",outer step size is small enough,infeasibility is within bounds
//--- * "inconsistent", if Lagrange multipliers increased beyond threshold given by MaxLagrangeMul
//--- * "too stringent",in other cases
   v=0;
   for(i=0;i<=nmain-1;i++)
      v=v+CMath::Sqr((state.m_xcur[i]-state.m_xprev[i])/state.m_seffective[i]);
   v=MathSqrt(v);
//--- check
   if(v<=(double)(state.m_outerepsx))
     {
      state.m_repterminationtype=4;
      //--- function call, return result
      return(Func_lbl_16(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- check
   if(state.m_maxits>0)
     {
      state.m_itsleft=state.m_itsleft-state.m_cgrep.m_iterationscount;
      //--- check
      if(state.m_itsleft<=0)
        {
         state.m_repterminationtype=5;
         //--- function call, return result
         return(Func_lbl_16(state,nmain,nslack,m,i,j,b,v,vv));
        }
     }
//--- check
   if(state.m_repouteriterationscount>=m_maxouterits)
     {
      state.m_repterminationtype=5;
      //--- function call, return result
      return(Func_lbl_16(state,nmain,nslack,m,i,j,b,v,vv));
     }
//--- Next iteration
   for(int i_=0;i_<=nmain+nslack-1;i_++)
      state.m_xprev[i_]=state.m_xcur[i_];
//--- function call, return result
   return(Func_lbl_15(state,nmain,nslack,m,i,j,b,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CMinLBFGS                                    |
//+------------------------------------------------------------------+
class CMinLBFGSState
  {
public:
   //--- variables
   int               m_n;
   int               m_m;
   double            m_epsg;
   double            m_epsf;
   double            m_epsx;
   int               m_maxits;
   bool              m_xrep;
   double            m_stpmax;
   double            m_diffstep;
   int               m_nfev;
   int               m_mcstage;
   int               m_k;
   int               m_q;
   int               m_p;
   double            m_stp;
   double            m_fold;
   double            m_trimthreshold;
   int               m_prectype;
   double            m_gammak;
   double            m_fbase;
   double            m_fm2;
   double            m_fm1;
   double            m_fp1;
   double            m_fp2;
   double            m_f;
   bool              m_needf;
   bool              m_needfg;
   bool              m_xupdated;
   RCommState        m_rstate;
   int               m_repiterationscount;
   int               m_repnfev;
   int               m_repterminationtype;
   CLinMinState      m_lstate;
   //--- arrays
   double            m_s[];
   double            m_rho[];
   double            m_theta[];
   double            m_d[];
   double            m_work[];
   double            m_diagh[];
   double            m_autobuf[];
   double            m_x[];
   double            m_g[];
   //--- matrix
   CMatrixDouble     m_yk;
   CMatrixDouble     m_sk;
   CMatrixDouble     m_denseh;
   //--- constructor, destructor
                     CMinLBFGSState(void);
                    ~CMinLBFGSState(void);
   //--- copy
   void              Copy(CMinLBFGSState &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinLBFGSState::CMinLBFGSState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinLBFGSState::~CMinLBFGSState(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMinLBFGSState::Copy(CMinLBFGSState &obj)
  {
//--- copy variables
   m_n=obj.m_n;
   m_m=obj.m_m;
   m_epsg=obj.m_epsg;
   m_epsf=obj.m_epsf;
   m_epsx=obj.m_epsx;
   m_maxits=obj.m_maxits;
   m_xrep=obj.m_xrep;
   m_stpmax=obj.m_stpmax;
   m_diffstep=obj.m_diffstep;
   m_nfev=obj.m_nfev;
   m_mcstage=obj.m_mcstage;
   m_k=obj.m_k;
   m_q=obj.m_q;
   m_p=obj.m_p;
   m_stp=obj.m_stp;
   m_fold=obj.m_fold;
   m_trimthreshold=obj.m_trimthreshold;
   m_prectype=obj.m_prectype;
   m_gammak=obj.m_gammak;
   m_fbase=obj.m_fbase;
   m_fm2=obj.m_fm2;
   m_fm1=obj.m_fm1;
   m_fp1=obj.m_fp1;
   m_fp2=obj.m_fp2;
   m_f=obj.m_f;
   m_needf=obj.m_needf;
   m_needfg=obj.m_needfg;
   m_xupdated=obj.m_xupdated;
   m_repiterationscount=obj.m_repiterationscount;
   m_repnfev=obj.m_repnfev;
   m_repterminationtype=obj.m_repterminationtype;
   m_rstate.Copy(obj.m_rstate);
   m_lstate.Copy(obj.m_lstate);
//--- copy arrays
   ArrayCopy(m_s,obj.m_s);
   ArrayCopy(m_rho,obj.m_rho);
   ArrayCopy(m_theta,obj.m_theta);
   ArrayCopy(m_d,obj.m_d);
   ArrayCopy(m_work,obj.m_work);
   ArrayCopy(m_diagh,obj.m_diagh);
   ArrayCopy(m_autobuf,obj.m_autobuf);
   ArrayCopy(m_x,obj.m_x);
   ArrayCopy(m_g,obj.m_g);
//--- copy matrix
   m_yk=obj.m_yk;
   m_sk=obj.m_sk;
   m_denseh=obj.m_denseh;
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CMinLBFGSState                   |
//+------------------------------------------------------------------+
class CMinLBFGSStateShell
  {
private:
   CMinLBFGSState    m_innerobj;
public:
   //--- constructors, destructor
                     CMinLBFGSStateShell(void);
                     CMinLBFGSStateShell(CMinLBFGSState &obj);
                    ~CMinLBFGSStateShell(void);
   //--- methods
   bool              GetNeedF(void);
   void              SetNeedF(const bool b);
   bool              GetNeedFG(void);
   void              SetNeedFG(const bool b);
   bool              GetXUpdated(void);
   void              SetXUpdated(const bool b);
   double            GetF(void);
   void              SetF(const double d);
   CMinLBFGSState   *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinLBFGSStateShell::CMinLBFGSStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMinLBFGSStateShell::CMinLBFGSStateShell(CMinLBFGSState &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinLBFGSStateShell::~CMinLBFGSStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needf                          |
//+------------------------------------------------------------------+
bool CMinLBFGSStateShell::GetNeedF(void)
  {
//--- return result
   return(m_innerobj.m_needf);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needf                         |
//+------------------------------------------------------------------+
void CMinLBFGSStateShell::SetNeedF(const bool b)
  {
//--- change value
   m_innerobj.m_needf=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needfg                         |
//+------------------------------------------------------------------+
bool CMinLBFGSStateShell::GetNeedFG(void)
  {
//--- return result
   return(m_innerobj.m_needfg);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needfg                        |
//+------------------------------------------------------------------+
void CMinLBFGSStateShell::SetNeedFG(const bool b)
  {
//--- change value
   m_innerobj.m_needfg=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable xupdated                       |
//+------------------------------------------------------------------+
bool CMinLBFGSStateShell::GetXUpdated(void)
  {
//--- return result
   return(m_innerobj.m_xupdated);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable xupdated                      |
//+------------------------------------------------------------------+
void CMinLBFGSStateShell::SetXUpdated(const bool b)
  {
//--- change value
   m_innerobj.m_xupdated=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable f                              |
//+------------------------------------------------------------------+
double CMinLBFGSStateShell::GetF(void)
  {
//--- return result
   return(m_innerobj.m_f);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable f                             |
//+------------------------------------------------------------------+
void CMinLBFGSStateShell::SetF(const double d)
  {
//--- change value
   m_innerobj.m_f=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMinLBFGSState *CMinLBFGSStateShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CMinLFBFGS                                   |
//+------------------------------------------------------------------+
class CMinLBFGSReport
  {
public:
   //--- variables
   int               m_iterationscount;
   int               m_nfev;
   int               m_terminationtype;
   //--- constructor, destructor
                     CMinLBFGSReport(void);
                    ~CMinLBFGSReport(void);
   //--- copy
   void              Copy(CMinLBFGSReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinLBFGSReport::CMinLBFGSReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinLBFGSReport::~CMinLBFGSReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMinLBFGSReport::Copy(CMinLBFGSReport &obj)
  {
//--- copy variables
   m_iterationscount=obj.m_iterationscount;
   m_nfev=obj.m_nfev;
   m_terminationtype=obj.m_terminationtype;
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CMinLBFGSReport                  |
//+------------------------------------------------------------------+
class CMinLBFGSReportShell
  {
private:
   CMinLBFGSReport   m_innerobj;
public:
   //--- constructors, destructor
                     CMinLBFGSReportShell(void);
                     CMinLBFGSReportShell(CMinLBFGSReport &obj);
                    ~CMinLBFGSReportShell(void);
   //--- methods
   int               GetIterationsCount(void);
   void              SetIterationsCount(const int i);
   int               GetNFev(void);
   void              SetNFev(const int i);
   int               GetTerminationType(void);
   void              SetTerminationType(const int i);
   CMinLBFGSReport *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinLBFGSReportShell::CMinLBFGSReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMinLBFGSReportShell::CMinLBFGSReportShell(CMinLBFGSReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinLBFGSReportShell::~CMinLBFGSReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable iterationscount                |
//+------------------------------------------------------------------+
int CMinLBFGSReportShell::GetIterationsCount(void)
  {
//--- return result
   return(m_innerobj.m_iterationscount);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable iterationscount               |
//+------------------------------------------------------------------+
void CMinLBFGSReportShell::SetIterationsCount(const int i)
  {
//--- change value
   m_innerobj.m_iterationscount=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable nfev                           |
//+------------------------------------------------------------------+
int CMinLBFGSReportShell::GetNFev(void)
  {
//--- return result
   return(m_innerobj.m_nfev);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable nfev                          |
//+------------------------------------------------------------------+
void CMinLBFGSReportShell::SetNFev(const int i)
  {
//--- change value
   m_innerobj.m_nfev=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable terminationtype                |
//+------------------------------------------------------------------+
int CMinLBFGSReportShell::GetTerminationType(void)
  {
//--- return result
   return(m_innerobj.m_terminationtype);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable terminationtype               |
//+------------------------------------------------------------------+
void CMinLBFGSReportShell::SetTerminationType(const int i)
  {
//--- change value
   m_innerobj.m_terminationtype=i;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMinLBFGSReport *CMinLBFGSReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Limited memory BFGS method for large scale optimization          |
//+------------------------------------------------------------------+
class CMinLBFGS
  {
private:
   //--- private method
   static void       ClearRequestFields(CMinLBFGSState &state);
   //--- auxiliary functions for MinLBFGSiteration
   static void       Func_lbl_rcomm(CMinLBFGSState &state,int n,int m,int i,int j,int ic,int mcinfo,double v,double vv);
   static bool       Func_lbl_16(CMinLBFGSState &state,int &n,int &m,int &i,int &j,int &ic,int &mcinfo,double &v,double &vv);
   static bool       Func_lbl_19(CMinLBFGSState &state,int &n,int &m,int &i,int &j,int &ic,int &mcinfo,double &v,double &vv);
   static bool       Func_lbl_21(CMinLBFGSState &state,int &n,int &m,int &i,int &j,int &ic,int &mcinfo,double &v,double &vv);
   static bool       Func_lbl_23(CMinLBFGSState &state,int &n,int &m,int &i,int &j,int &ic,int &mcinfo,double &v,double &vv);
   static bool       Func_lbl_27(CMinLBFGSState &state,int &n,int &m,int &i,int &j,int &ic,int &mcinfo,double &v,double &vv);
   static bool       Func_lbl_30(CMinLBFGSState &state,int &n,int &m,int &i,int &j,int &ic,int &mcinfo,double &v,double &vv);
public:
   //--- constant
   static const double m_gtol;
   //--- constructor, destructor
                     CMinLBFGS(void);
                    ~CMinLBFGS(void);
   //--- public methods
   static void       MinLBFGSCreate(const int n,const int m,double &x[],CMinLBFGSState &state);
   static void       MinLBFGSCreateF(const int n,const int m,double &x[],const double diffstep,CMinLBFGSState &state);
   static void       MinLBFGSSetCond(CMinLBFGSState &state,const double epsg,const double epsf,double epsx,const int maxits);
   static void       MinLBFGSSetXRep(CMinLBFGSState &state,const bool needxrep);
   static void       MinLBFGSSetStpMax(CMinLBFGSState &state,const double stpmax);
   static void       MinLBFGSSetScale(CMinLBFGSState &state,double &s[]);
   static void       MinLBFGSCreateX(const int n,const int m,double &x[],int flags,const double diffstep,CMinLBFGSState &state);
   static void       MinLBFGSSetPrecDefault(CMinLBFGSState &state);
   static void       MinLBFGSSetPrecCholesky(CMinLBFGSState &state,CMatrixDouble &p,const bool isupper);
   static void       MinLBFGSSetPrecDiag(CMinLBFGSState &state,double &d[]);
   static void       MinLBFGSSetPrecScale(CMinLBFGSState &state);
   static void       MinLBFGSResults(CMinLBFGSState &state,double &x[],CMinLBFGSReport &rep);
   static void       MinLBFGSresultsbuf(CMinLBFGSState &state,double &x[],CMinLBFGSReport &rep);
   static void       MinLBFGSRestartFrom(CMinLBFGSState &state,double &x[]);
   static bool       MinLBFGSIteration(CMinLBFGSState &state);
  };
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const double CMinLBFGS::m_gtol=0.4;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinLBFGS::CMinLBFGS(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinLBFGS::~CMinLBFGS(void)
  {

  }
//+------------------------------------------------------------------+
//|         LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION  |
//| DESCRIPTION:                                                     |
//| The subroutine minimizes function F(x) of N arguments by using a |
//| quasi - Newton method (LBFGS scheme) which is optimized to use a |
//| minimum  amount of memory.                                       |
//| The subroutine generates the approximation of an inverse Hessian |
//| matrix by using information about the last M steps of the        |
//| algorithm (instead of N). It lessens a required amount of memory |
//| from a value of order N^2 to a value of order 2*N*M.             |
//| REQUIREMENTS:                                                    |
//| Algorithm will request following information during its          |
//| operation:                                                       |
//| * function value F and its gradient G (simultaneously) at given  |
//| point X                                                          |
//| USAGE:                                                           |
//| 1. User initializes algorithm state with MinLBFGSCreate() call   |
//| 2. User tunes solver parameters with MinLBFGSSetCond()           |
//|    MinLBFGSSetStpMax() and other functions                       |
//| 3. User calls MinLBFGSOptimize() function which takes algorithm  |
//|    state and pointer (delegate, etc.) to callback function which |
//|    calculates F/G.                                               |
//| 4. User calls MinLBFGSResults() to get solution                  |
//| 5. Optionally user may call MinLBFGSRestartFrom() to solve       |
//|    another problem with same N/M but another starting point      |
//|    and/or another function. MinLBFGSRestartFrom() allows to reuse|
//|    already initialized structure.                                |
//| INPUT PARAMETERS:                                                |
//|     N       -   problem dimension. N>0                           |
//|     M       -   number of corrections in the BFGS scheme of      |
//|                 Hessian approximation update. Recommended value: |
//|                 3<=M<=7. The smaller value causes worse          |
//|                 convergence, the bigger will not cause a         |
//|                 considerably better convergence, but will cause  |
//|                 a fall in  the performance. M<=N.                |
//|     X       -   initial solution approximation, array[0..N-1].   |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//| NOTES:                                                           |
//| 1. you may tune stopping conditions with MinLBFGSSetCond()       |
//|    function                                                      |
//| 2. if target function contains exp() or other fast growing       |
//|    functions, and optimization algorithm makes too large steps   |
//|    which leads to overflow, use MinLBFGSSetStpMax() function to  |
//|    bound algorithm's steps. However, L-BFGS rarely needs such a  |
//|    tuning.                                                       |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSCreate(const int n,const int m,double &x[],
                                      CMinLBFGSState &state)
  {
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1"))
      return;
//--- check
   if(!CAp::Assert(m<=n,__FUNCTION__+": M>N"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- function call
   MinLBFGSCreateX(n,m,x,0,0.0,state);
  }
//+------------------------------------------------------------------+
//| The subroutine is finite difference variant of MinLBFGSCreate(). |
//| It uses finite differences in order to differentiate target      |
//| function.                                                        |
//| Description below contains information which is specific to this |
//| function only. We recommend to read comments on MinLBFGSCreate() |
//| in order to get more information about creation of LBFGS         |
//| optimizer.                                                       |
//| INPUT PARAMETERS:                                                |
//|     N       -   problem dimension, N>0:                          |
//|                 * if given, only leading N elements of X are used|
//|                 * if not given, automatically determined from    |
//|                   size of X                                      |
//|     M       -   number of corrections in the BFGS scheme of      |
//|                 Hessian approximation update. Recommended value: |
//|                 3<=M<=7. The smaller value causes worse          |
//|                 convergence, the bigger will not cause a         |
//|                 considerably better convergence, but will cause a|
//|                 fall in  the performance. M<=N.                  |
//|     X       -   starting point, array[0..N-1].                   |
//|     DiffStep-   differentiation step, >0                         |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//| NOTES:                                                           |
//| 1. algorithm uses 4-point central formula for differentiation.   |
//| 2. differentiation step along I-th axis is equal to DiffStep*S[I]|
//|    where S[] is scaling vector which can be set by               |
//|    MinLBFGSSetScale() call.                                      |
//| 3. we recommend you to use moderate values of differentiation    |
//|    step. Too large step will result in too large truncation      |
//|    errors, while too small step will result in too large         |
//|    numerical errors. 1.0E-6 can be good value to start with.     |
//| 4. Numerical differentiation is very inefficient - one gradient  |
//|    calculation needs 4*N function evaluations. This function will|
//|    work for any N - either small (1...10), moderate (10...100) or|
//|    large (100...). However, performance penalty will be too      |
//|    severe for any N's except for small ones.                     |
//|    We should also say that code which relies on numerical        |
//|    differentiation is less robust and precise. LBFGS needs exact |
//|    gradient values. Imprecise gradient may slow down convergence,|
//|    especially on highly nonlinear problems.                      |
//|    Thus we recommend to use this function for fast prototyping on|
//|    small- dimensional problems only, and to implement analytical |
//|    gradient as soon as possible.                                 |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSCreateF(const int n,const int m,double &x[],
                                       const double diffstep,CMinLBFGSState &state)
  {
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N too small!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1"))
      return;
//--- check
   if(!CAp::Assert(m<=n,__FUNCTION__+": M>N"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(diffstep),__FUNCTION__+": DiffStep is infinite or NaN!"))
      return;
//--- check
   if(!CAp::Assert(diffstep>0.0,__FUNCTION__+": DiffStep is non-positive!"))
      return;
//--- function call
   MinLBFGSCreateX(n,m,x,0,diffstep,state);
  }
//+------------------------------------------------------------------+
//| This function sets stopping conditions for L-BFGS optimization   |
//| algorithm.                                                       |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     EpsG    -   >=0                                              |
//|                 The subroutine finishes its work if the condition|
//|                 |v|<EpsG is satisfied, where:                    |
//|                 * |.| means Euclidian norm                       |
//|                 * v - scaled gradient vector, v[i]=g[i]*s[i]     |
//|                 * g - gradient                                   |
//|                 * s - scaling coefficients set by                |
//|                       MinLBFGSSetScale()                         |
//|     EpsF    -   >=0                                              |
//|                 The subroutine finishes its work if on k+1-th    |
//|                 iteration the condition |F(k+1)-F(k)| <=         |
//|                 <= EpsF*max{|F(k)|,|F(k+1)|,1} is satisfied.     |
//|     EpsX    -   >=0                                              |
//|                 The subroutine finishes its work if on k+1-th    |
//|                 iteration the condition |v|<=EpsX is fulfilled,  |
//|                 where:                                           |
//|                 * |.| means Euclidian norm                       |
//|                 * v - scaled step vector, v[i]=dx[i]/s[i]        |
//|                 * dx - ste pvector, dx=X(k+1)-X(k)               |
//|                 * s - scaling coefficients set by                |
//|                   MinLBFGSSetScale()                             |
//|     MaxIts  -   maximum number of iterations. If MaxIts=0, the   |
//|                 number of iterations is unlimited.               |
//| Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will|
//| lead to automatic stopping criterion selection (small EpsX).     |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSSetCond(CMinLBFGSState &state,const double epsg,
                                       const double epsf,double epsx,
                                       const int maxits)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsg),__FUNCTION__+": EpsG is not finite number!"))
      return;
//--- check
   if(!CAp::Assert(epsg>=0.0,__FUNCTION__+": negative EpsG!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsf),__FUNCTION__+": EpsF is not finite number!"))
      return;
//--- check
   if(!CAp::Assert(epsf>=0.0,__FUNCTION__+": negative EpsF!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsx),__FUNCTION__+": EpsX is not finite number!"))
      return;
//--- check
   if(!CAp::Assert(epsx>=0.0,__FUNCTION__+": negative EpsX!"))
      return;
//--- check
   if(!CAp::Assert(maxits>=0,__FUNCTION__+": negative MaxIts!"))
      return;
//--- check
   if(((epsg==0.0 && epsf==0.0) && epsx==0.0) && maxits==0)
      epsx=1.0E-6;
//--- change values
   state.m_epsg=epsg;
   state.m_epsf=epsf;
   state.m_epsx=epsx;
   state.m_maxits=maxits;
  }
//+------------------------------------------------------------------+
//| This function turns on/off reporting.                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     NeedXRep-   whether iteration reports are needed or not      |
//| If NeedXRep is True, algorithm will call rep() callback function |
//| if it is provided to MinLBFGSOptimize().                         |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSSetXRep(CMinLBFGSState &state,const bool needxrep)
  {
//--- change value
   state.m_xrep=needxrep;
  }
//+------------------------------------------------------------------+
//| This function sets maximum step length                           |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     StpMax  -   maximum step length, >=0. Set StpMax to 0.0      |
//|                 (default), if you don't want to limit step       |
//|                 length.                                          |
//| Use this subroutine when you optimize target function which      |
//| contains exp() or other fast growing functions, and optimization |
//| algorithm makes too large steps which leads to overflow. This    |
//| function allows us to reject steps that are too large (and       |
//| therefore expose us to the possible overflow) without actually   |
//| calculating function value at the x+stp*d.                       |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSSetStpMax(CMinLBFGSState &state,
                                         const double stpmax)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(stpmax),__FUNCTION__+": StpMax is not finite!"))
      return;
//--- check
   if(!CAp::Assert(stpmax>=0.0,__FUNCTION__+": StpMax<0!"))
      return;
//--- change value
   state.m_stpmax=stpmax;
  }
//+------------------------------------------------------------------+
//| This function sets scaling coefficients for LBFGS optimizer.     |
//| ALGLIB optimizers use scaling matrices to test stopping          |
//| conditions (step size and gradient are scaled before comparison  |
//| with tolerances). Scale of the I-th variable is a translation    |
//| invariant measure of:                                            |
//| a) "how large" the variable is                                   |
//| b) how large the step should be to make significant changes in   |
//|    the function                                                  |
//| Scaling is also used by finite difference variant of the         |
//| optimizer - step along I-th axis is equal to DiffStep*S[I].      |
//| In most optimizers (and in the LBFGS too) scaling is NOT a form  |
//| of preconditioning. It just affects stopping conditions. You     |
//| should set preconditioner by separate call to one of the         |
//| MinLBFGSSetPrec...() functions.                                  |
//| There is special preconditioning mode, however, which uses       |
//| scaling coefficients to form diagonal preconditioning matrix.    |
//| You can turn this mode on, if you want. But  you should          |
//| understand that scaling is not the same thing as                 |
//| preconditioning - these are two different, although related      |
//| forms of tuning solver.                                          |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure stores algorithm state                 |
//|     S       -   array[N], non-zero scaling coefficients          |
//|                 S[i] may be negative, sign doesn't matter.       |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSSetScale(CMinLBFGSState &state,double &s[])
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(CAp::Len(s)>=state.m_n,__FUNCTION__+": Length(S)<N"))
      return;
   for(i=0;i<=state.m_n-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(s[i]),__FUNCTION__+": S contains infinite or NAN elements"))
         return;
      //--- check
      if(!CAp::Assert(s[i]!=0.0,__FUNCTION__+": S contains zero elements"))
         return;
      state.m_s[i]=MathAbs(s[i]);
     }
  }
//+------------------------------------------------------------------+
//| Extended subroutine for internal use only.                       |
//| Accepts additional parameters:                                   |
//|     Flags - additional settings:                                 |
//|             * Flags = 0     means no additional settings         |
//|             * Flags = 1     "do not allocate memory". used when  |
//|                             solving a many subsequent tasks with |
//|                             same N/M values. First call MUST be  |
//|                             without this flag bit set, subsequent|
//|                             calls of MinLBFGS with same          |
//|                             MinLBFGSState structure can set Flags|
//|                             to 1.                                |
//|     DiffStep - numerical differentiation step                    |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSCreateX(const int n,const int m,double &x[],
                                       int flags,const double diffstep,
                                       CMinLBFGSState &state)
  {
//--- create variables
   bool allocatemem;
   int  i=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N too small!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M too small!"))
      return;
//--- check
   if(!CAp::Assert(m<=n,__FUNCTION__+": M too large!"))
      return;
//--- Initialize
   state.m_diffstep=diffstep;
   state.m_n=n;
   state.m_m=m;
   allocatemem=flags%2==0;
   flags=flags/2;
//--- check
   if(allocatemem)
     {
      //--- allocation
      ArrayResizeAL(state.m_rho,m);
      ArrayResizeAL(state.m_theta,m);
      state.m_yk.Resize(m,n);
      state.m_sk.Resize(m,n);
      ArrayResizeAL(state.m_d,n);
      ArrayResizeAL(state.m_x,n);
      ArrayResizeAL(state.m_s,n);
      ArrayResizeAL(state.m_g,n);
      ArrayResizeAL(state.m_work,n);
     }
//--- function call
   MinLBFGSSetCond(state,0,0,0,0);
//--- function call
   MinLBFGSSetXRep(state,false);
//--- function call
   MinLBFGSSetStpMax(state,0);
//--- function call
   MinLBFGSRestartFrom(state,x);
//--- change values
   for(i=0;i<=n-1;i++)
      state.m_s[i]=1.0;
   state.m_prectype=0;
  }
//+------------------------------------------------------------------+
//| Modification of the preconditioner: default preconditioner       |
//| (simple scaling, same for all elements of X) is used.            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//| NOTE: you can change preconditioner "on the fly", during         |
//| algorithm iterations.                                            |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSSetPrecDefault(CMinLBFGSState &state)
  {
//--- change value
   state.m_prectype=0;
  }
//+------------------------------------------------------------------+
//| Modification of the preconditioner: Cholesky factorization of    |
//| approximate Hessian is used.                                     |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     P       -   triangular preconditioner, Cholesky factorization|
//|                 of the approximate Hessian. array[0..N-1,0..N-1],|
//|                 (if larger, only leading N elements are used).   |
//|     IsUpper -   whether upper or lower triangle of P is given    |
//|                 (other triangle is not referenced)               |
//| After call to this function preconditioner is changed to P (P is |
//| copied into the internal buffer).                                |
//| NOTE: you can change preconditioner "on the fly", during         |
//| algorithm iterations.                                            |
//| NOTE 2: P should be nonsingular. Exception will be thrown        |
//| otherwise.                                                       |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSSetPrecCholesky(CMinLBFGSState &state,
                                               CMatrixDouble &p,
                                               const bool isupper)
  {
//--- create variables
   int    i=0;
   double mx=0;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteRTrMatrix(p,state.m_n,isupper),__FUNCTION__+": P contains infinite or NAN values!"))
      return;
//--- initialization
   mx=0;
   for(i=0;i<=state.m_n-1;i++)
      mx=MathMax(mx,MathAbs(p[i][i]));
//--- check
   if(!CAp::Assert((double)(mx)>0.0,__FUNCTION__+": P is strictly singular!"))
      return;
//--- check
   if(CAp::Rows(state.m_denseh)<state.m_n||CAp::Cols(state.m_denseh)<state.m_n)
      state.m_denseh.Resize(state.m_n,state.m_n);
//--- initialization
   state.m_prectype=1;
//--- check
   if(isupper)
      CAblas::RMatrixCopy(state.m_n,state.m_n,p,0,0,state.m_denseh,0,0);
   else
      CAblas::RMatrixTranspose(state.m_n,state.m_n,p,0,0,state.m_denseh,0,0);
  }
//+------------------------------------------------------------------+
//| Modification of the preconditioner: diagonal of approximate      |
//| Hessian is used.                                                 |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     D       -   diagonal of the approximate Hessian,             |
//|                 array[0..N-1], (if larger, only leading N        |
//|                 elements are used).                              |
//| NOTE: you can change preconditioner "on the fly", during         |
//| algorithm iterations.                                            |
//| NOTE 2: D[i] should be positive. Exception will be thrown        |
//| otherwise.                                                       |
//| NOTE 3: you should pass diagonal of approximate Hessian - NOT    |
//| ITS INVERSE.                                                     |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSSetPrecDiag(CMinLBFGSState &state,double &d[])
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(CAp::Len(d)>=state.m_n,__FUNCTION__+": D is too short"))
      return;
   for(i=0;i<=state.m_n-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(d[i]),__FUNCTION__+": D contains infinite or NAN elements"))
         return;
      //--- check
      if(!CAp::Assert((double)(d[i])>0.0,__FUNCTION__+": D contains non-positive elements"))
         return;
     }
//--- function call
   CApServ::RVectorSetLengthAtLeast(state.m_diagh,state.m_n);
//--- change values
   state.m_prectype=2;
   for(i=0;i<=state.m_n-1;i++)
      state.m_diagh[i]=d[i];
  }
//+------------------------------------------------------------------+
//| Modification of the preconditioner: scale-based diagonal         |
//| preconditioning.                                                 |
//| This preconditioning mode can be useful when you don't have      |
//| approximate diagonal of Hessian, but you know that your variables|
//| are badly scaled (for example, one variable is in [1,10], and    |
//| another in [1000,100000]), and most part of the ill-conditioning |
//| comes from different scales of vars.                             |
//| In this case simple scale-based preconditioner, with H[i] =      |
//| = 1/(s[i]^2), can greatly improve convergence.                   |
//| IMPRTANT: you should set scale of your variables with            |
//| MinLBFGSSetScale() call (before or after MinLBFGSSetPrecScale()  |
//| call). Without knowledge of the scale of your variables          |
//| scale-based preconditioner will be just unit matrix.             |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSSetPrecScale(CMinLBFGSState &state)
  {
//--- change values
   state.m_prectype=3;
  }
//+------------------------------------------------------------------+
//| L-BFGS algorithm results                                         |
//| INPUT PARAMETERS:                                                |
//|     State   -   algorithm state                                  |
//| OUTPUT PARAMETERS:                                               |
//|     X       -   array[0..N-1], solution                          |
//|     Rep     -   optimization report:                             |
//|                 * Rep.TerminationType completetion code:         |
//|                     * -2    rounding errors prevent further      |
//|                             improvement. X contains best point   |
//|                             found.                               |
//|                     * -1    incorrect parameters were specified  |
//|                     *  1    relative function improvement is no  |
//|                             more than EpsF.                      |
//|                     *  2    relative step is no more than EpsX.  |
//|                     *  4    gradient norm is no more than EpsG   |
//|                     *  5    MaxIts steps was taken               |
//|                     *  7    stopping conditions are too          |
//|                             stringent, further improvement is    |
//|                             impossible                           |
//|                 * Rep.IterationsCount contains iterations count  |
//|                 * NFEV countains number of function calculations |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSResults(CMinLBFGSState &state,double &x[],
                                       CMinLBFGSReport &rep)
  {
//--- reset memory
   ArrayResizeAL(x,0);
//--- function call
   MinLBFGSresultsbuf(state,x,rep);
  }
//+------------------------------------------------------------------+
//| L-BFGS algorithm results                                         |
//| Buffered implementation of MinLBFGSResults which uses            |
//| pre-allocated buffer to store X[]. If buffer size is too small,  |
//| it resizes buffer. It is intended to be used in the inner cycles |
//| of performance critical algorithms where array reallocation      |
//| penalty is too large to be ignored.                              |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSresultsbuf(CMinLBFGSState &state,double &x[],
                                          CMinLBFGSReport &rep)
  {
//--- create a variable
   int i_=0;
//--- check
   if(CAp::Len(x)<state.m_n)
      ArrayResizeAL(x,state.m_n);
//--- copy
   for(i_=0;i_<=state.m_n-1;i_++)
      x[i_]=state.m_x[i_];
//--- change values
   rep.m_iterationscount=state.m_repiterationscount;
   rep.m_nfev=state.m_repnfev;
   rep.m_terminationtype=state.m_repterminationtype;
  }
//+------------------------------------------------------------------+
//| This  subroutine restarts LBFGS algorithm from new point. All    |
//| optimization parameters are left unchanged.                      |
//| This function allows to solve multiple optimization problems     |
//| (which must have same number of dimensions) without object       |
//| reallocation penalty.                                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure used to store algorithm state          |
//|     X       -   new starting point.                              |
//+------------------------------------------------------------------+
static void CMinLBFGS::MinLBFGSRestartFrom(CMinLBFGSState &state,double &x[])
  {
//--- create a variable
   int i_=0;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=state.m_n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,state.m_n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- copy
   for(i_=0;i_<=state.m_n-1;i_++)
      state.m_x[i_]=x[i_];
//--- allocation
   ArrayResizeAL(state.m_rstate.ia,6);
   ArrayResizeAL(state.m_rstate.ra,2);
   state.m_rstate.stage=-1;
//--- function call
   ClearRequestFields(state);
  }
//+------------------------------------------------------------------+
//| Clears request fileds (to be sure that we don't forgot to clear  |
//| something)                                                       |
//+------------------------------------------------------------------+
static void CMinLBFGS::ClearRequestFields(CMinLBFGSState &state)
  {
//--- change values
   state.m_needf=false;
   state.m_needfg=false;
   state.m_xupdated=false;
  }
//+------------------------------------------------------------------+
//| NOTES:                                                           |
//| 1. This function has two different implementations: one which    |
//|    uses exact (analytical) user-supplied gradient, and one which |
//|    uses function value only and numerically differentiates       |
//|    function in order to obtain gradient.                         |
//|    Depending on the specific function used to create optimizer   |
//|    object (either MinLBFGSCreate() for analytical gradient or    |
//|    MinLBFGSCreateF() for numerical differentiation) you should   |
//|    choose appropriate variant of MinLBFGSOptimize() - one which  |
//|    accepts function AND gradient or one which accepts function   |
//|    ONLY.                                                         |
//|    Be careful to choose variant of MinLBFGSOptimize() which      |
//|    corresponds to your optimization scheme! Table below lists    |
//|    different combinations of callback (function/gradient) passed |
//|    to MinLBFGSOptimize() and specific function used to create    |
//|    optimizer.                                                    |
//|                      |         USER PASSED TO MinLBFGSOptimize() |
//|    CREATED WITH      |  function only   |  function and gradient |
//|    ------------------------------------------------------------  |
//|    MinLBFGSCreateF() |     work                FAIL              |
//|    MinLBFGSCreate()  |     FAIL                work              |
//|    Here "FAIL" denotes inappropriate combinations of optimizer   |
//|    creation function and MinLBFGSOptimize() version. Attemps to  |
//|    use such combination (for example, to create optimizer with   |
//|    MinLBFGSCreateF() and to pass gradient information to         |
//|    MinCGOptimize()) will lead to exception being thrown. Either  |
//|    you did not pass gradient when it WAS needed or you passed    |
//|    gradient when it was NOT needed.                              |
//+------------------------------------------------------------------+
static bool CMinLBFGS::MinLBFGSIteration(CMinLBFGSState &state)
  {
//--- create variables
   int    n=0;
   int    m=0;
   int    i=0;
   int    j=0;
   int    ic=0;
   int    mcinfo=0;
   double v=0;
   double vv=0;
   int    i_=0;
//--- Reverse communication preparations
//--- I know it looks ugly, but it works the same way
//--- anywhere from C++ to Python.
//--- This code initializes locals by:
//--- * random values determined during code
//---   generation-on first subroutine call
//--- * values from previous call-on subsequent calls
   if(state.m_rstate.stage>=0)
     {
      //--- initialization
      n=state.m_rstate.ia[0];
      m=state.m_rstate.ia[1];
      i=state.m_rstate.ia[2];
      j=state.m_rstate.ia[3];
      ic=state.m_rstate.ia[4];
      mcinfo=state.m_rstate.ia[5];
      v=state.m_rstate.ra[0];
      vv=state.m_rstate.ra[1];
     }
   else
     {
      //--- initialization
      n=-983;
      m=-989;
      i=-834;
      j=900;
      ic=-287;
      mcinfo=364;
      v=214;
      vv=-338;
     }
//--- check
   if(state.m_rstate.stage==0)
     {
      //--- change value
      state.m_needfg=false;
      //--- function call
      COptServ::TrimPrepare(state.m_f,state.m_trimthreshold);
      //--- check
      if(!state.m_xrep)
         return(Func_lbl_19(state,n,m,i,j,ic,mcinfo,v,vv));
      //--- function call
      ClearRequestFields(state);
      //--- change values
      state.m_xupdated=true;
      state.m_rstate.stage=6;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==1)
     {
      //--- change values
      state.m_fbase=state.m_f;
      i=0;
      //--- function call, return result
      return(Func_lbl_16(state,n,m,i,j,ic,mcinfo,v,vv));
     }
//--- check
   if(state.m_rstate.stage==2)
     {
      //--- change values
      state.m_fm2=state.m_f;
      state.m_x[i]=v-0.5*state.m_diffstep*state.m_s[i];
      state.m_rstate.stage=3;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==3)
     {
      //--- change values
      state.m_fm1=state.m_f;
      state.m_x[i]=v+0.5*state.m_diffstep*state.m_s[i];
      state.m_rstate.stage=4;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==4)
     {
      //--- change values
      state.m_fp1=state.m_f;
      state.m_x[i]=v+state.m_diffstep*state.m_s[i];
      state.m_rstate.stage=5;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==5)
     {
      //--- change values
      state.m_fp2=state.m_f;
      state.m_x[i]=v;
      state.m_g[i]=(8*(state.m_fp1-state.m_fm1)-(state.m_fp2-state.m_fm2))/(6*state.m_diffstep*state.m_s[i]);
      i=i+1;
      //--- function call, return result
      return(Func_lbl_16(state,n,m,i,j,ic,mcinfo,v,vv));
     }
//--- check
   if(state.m_rstate.stage==6)
     {
      //--- change value
      state.m_xupdated=false;
      //--- function call, return result
      return(Func_lbl_19(state,n,m,i,j,ic,mcinfo,v,vv));
     }
//--- check
   if(state.m_rstate.stage==7)
     {
      //--- change value
      state.m_needfg=false;
      //--- function call
      COptServ::TrimFunction(state.m_f,state.m_g,n,state.m_trimthreshold);
      //--- function call
      CLinMin::MCSrch(n,state.m_x,state.m_f,state.m_g,state.m_d,state.m_stp,state.m_stpmax,m_gtol,mcinfo,state.m_nfev,state.m_work,state.m_lstate,state.m_mcstage);
      //--- function call, return result
      return(Func_lbl_23(state,n,m,i,j,ic,mcinfo,v,vv));
     }
//--- check
   if(state.m_rstate.stage==8)
     {
      //--- change values
      state.m_fbase=state.m_f;
      i=0;
      //--- function call, return result
      return(Func_lbl_27(state,n,m,i,j,ic,mcinfo,v,vv));
     }
//--- check
   if(state.m_rstate.stage==9)
     {
      //--- change values
      state.m_fm2=state.m_f;
      state.m_x[i]=v-0.5*state.m_diffstep*state.m_s[i];
      state.m_rstate.stage=10;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==10)
     {
      //--- change values
      state.m_fm1=state.m_f;
      state.m_x[i]=v+0.5*state.m_diffstep*state.m_s[i];
      state.m_rstate.stage=11;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==11)
     {
      //--- change values
      state.m_fp1=state.m_f;
      state.m_x[i]=v+state.m_diffstep*state.m_s[i];
      state.m_rstate.stage=12;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==12)
     {
      //--- change values
      state.m_fp2=state.m_f;
      state.m_x[i]=v;
      state.m_g[i]=(8*(state.m_fp1-state.m_fm1)-(state.m_fp2-state.m_fm2))/(6*state.m_diffstep*state.m_s[i]);
      i=i+1;
      //--- function call, return result
      return(Func_lbl_27(state,n,m,i,j,ic,mcinfo,v,vv));
     }
//--- check
   if(state.m_rstate.stage==13)
     {
      //--- change value
      state.m_xupdated=false;
      //--- function call, return result
      return(Func_lbl_30(state,n,m,i,j,ic,mcinfo,v,vv));
     }
//--- Routine body
//--- Unload frequently used variables from State structure
//--- (just for typing convinience)
   n=state.m_n;
   m=state.m_m;
   state.m_repterminationtype=0;
   state.m_repiterationscount=0;
   state.m_repnfev=0;
//--- Calculate F/G at the initial point
   ClearRequestFields(state);
//--- check
   if(state.m_diffstep!=0.0)
     {
      //--- change values
      state.m_needf=true;
      state.m_rstate.stage=1;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
      //--- return result
      return(true);
     }
//--- change values
   state.m_needfg=true;
   state.m_rstate.stage=0;
//--- Saving state
   Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLBFGSIteration. Is a product to get    |
//| rid of the operator unconditional jump goto.                     |
//+------------------------------------------------------------------+
static void CMinLBFGS::Func_lbl_rcomm(CMinLBFGSState &state,int n,int m,
                                      int i,int j,int ic,int mcinfo,
                                      double v,double vv)
  {
//--- save
   state.m_rstate.ia[0]=n;
   state.m_rstate.ia[1]=m;
   state.m_rstate.ia[2]=i;
   state.m_rstate.ia[3]=j;
   state.m_rstate.ia[4]=ic;
   state.m_rstate.ia[5]=mcinfo;
   state.m_rstate.ra[0]=v;
   state.m_rstate.ra[1]=vv;
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLBFGSIteration. Is a product to get    |
//| rid of the operator unconditional jump goto.                     |
//+------------------------------------------------------------------+
static bool CMinLBFGS::Func_lbl_16(CMinLBFGSState &state,int &n,int &m,
                                   int &i,int &j,int &ic,int &mcinfo,
                                   double &v,double &vv)
  {
//--- check
   if(i>n-1)
     {
      //--- change values
      state.m_f=state.m_fbase;
      state.m_needf=false;
      //--- function call
      COptServ::TrimPrepare(state.m_f,state.m_trimthreshold);
      //--- check
      if(!state.m_xrep)
         return(Func_lbl_19(state,n,m,i,j,ic,mcinfo,v,vv));
      //--- function call
      ClearRequestFields(state);
      //--- change values
      state.m_xupdated=true;
      state.m_rstate.stage=6;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
      //--- return result
      return(true);
     }
//--- change values
   v=state.m_x[i];
   state.m_x[i]=v-state.m_diffstep*state.m_s[i];
   state.m_rstate.stage=2;
//--- Saving state
   Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLBFGSIteration. Is a product to get    |
//| rid of the operator unconditional jump goto.                     |
//+------------------------------------------------------------------+
static bool CMinLBFGS::Func_lbl_19(CMinLBFGSState &state,int &n,int &m,
                                   int &i,int &j,int &ic,int &mcinfo,
                                   double &v,double &vv)
  {
//--- change values
   state.m_repnfev=1;
   state.m_fold=state.m_f;
//--- calculation
   v=0;
   for(i=0;i<=n-1;i++)
      v=v+CMath::Sqr(state.m_g[i]*state.m_s[i]);
//--- check
   if(MathSqrt(v)<=state.m_epsg)
     {
      state.m_repterminationtype=4;
      //--- return result
      return(false);
     }
//--- Choose initial step and direction.
//--- Apply preconditioner, if we have something other than default.
   for(int i_=0;i_<=n-1;i_++)
      state.m_d[i_]=-state.m_g[i_];
//--- check
   if(state.m_prectype==0)
     {
      //--- Default preconditioner is used, but we can't use it before iterations will start
      v=0.0;
      for(int i_=0;i_<=n-1;i_++)
         v+=state.m_g[i_]*state.m_g[i_];
      v=MathSqrt(v);
      //--- check
      if(state.m_stpmax==0.0)
         state.m_stp=MathMin(1.0/v,1);
      else
         state.m_stp=MathMin(1.0/v,state.m_stpmax);
     }
//--- check
   if(state.m_prectype==1)
     {
      //--- Cholesky preconditioner is used
      CFbls::FblsCholeskySolve(state.m_denseh,1.0,n,true,state.m_d,state.m_autobuf);
      state.m_stp=1;
     }
//--- check
   if(state.m_prectype==2)
     {
      //--- diagonal approximation is used
      for(i=0;i<=n-1;i++)
         state.m_d[i]=state.m_d[i]/state.m_diagh[i];
      state.m_stp=1;
     }
//--- check
   if(state.m_prectype==3)
     {
      //--- scale-based preconditioner is used
      for(i=0;i<=n-1;i++)
         state.m_d[i]=state.m_d[i]*state.m_s[i]*state.m_s[i];
      state.m_stp=1;
     }
//--- Main cycle
   state.m_k=0;
//--- function call, return result
   return(Func_lbl_21(state,n,m,i,j,ic,mcinfo,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLBFGSIteration. Is a product to get    |
//| rid of the operator unconditional jump goto.                     |
//+------------------------------------------------------------------+
static bool CMinLBFGS::Func_lbl_21(CMinLBFGSState &state,int &n,int &m,
                                   int &i,int &j,int &ic,int &mcinfo,
                                   double &v,double &vv)
  {
//--- Main cycle: prepare to 1-D line search
   state.m_p=state.m_k % m;
   state.m_q=MathMin(state.m_k, m-1);
//--- Store X[k], G[k]
   for(int i_=0;i_<=n-1;i_++)
      state.m_sk[state.m_p].Set(i_,-state.m_x[i_]);
   for(int i_=0;i_<=n-1;i_++)
      state.m_yk[state.m_p].Set(i_,-state.m_g[i_]);
//--- Minimize F(x+alpha*d)
//--- Calculate S[k], Y[k]
   state.m_mcstage=0;
//--- check
   if(state.m_k!=0)
      state.m_stp=1.0;
//--- function call
   CLinMin::LinMinNormalized(state.m_d,state.m_stp,n);
//--- function call
   CLinMin::MCSrch(n,state.m_x,state.m_f,state.m_g,state.m_d,state.m_stp,state.m_stpmax,m_gtol,mcinfo,state.m_nfev,state.m_work,state.m_lstate,state.m_mcstage);
//--- function call, return result
   return(Func_lbl_23(state,n,m,i,j,ic,mcinfo,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLBFGSIteration. Is a product to get    |
//| rid of the operator unconditional jump goto.                     |
//+------------------------------------------------------------------+
static bool CMinLBFGS::Func_lbl_23(CMinLBFGSState &state,int &n,int &m,
                                   int &i,int &j,int &ic,int &mcinfo,
                                   double &v,double &vv)
  {
//--- check
   if(state.m_mcstage==0)
     {
      //--- check
      if(!state.m_xrep)
         return(Func_lbl_30(state,n,m,i,j,ic,mcinfo,v,vv));
      //--- report
      ClearRequestFields(state);
      //--- change values
      state.m_xupdated=true;
      state.m_rstate.stage=13;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
      //--- return result
      return(true);
     }
//--- function call
   ClearRequestFields(state);
//--- check
   if((double)(state.m_diffstep)!=0.0)
     {
      //--- change values
      state.m_needf=true;
      state.m_rstate.stage=8;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
      //--- return result
      return(true);
     }
//--- change values
   state.m_needfg=true;
   state.m_rstate.stage=7;
//--- Saving state
   Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLBFGSIteration. Is a product to get    |
//| rid of the operator unconditional jump goto.                     |
//+------------------------------------------------------------------+
static bool CMinLBFGS::Func_lbl_27(CMinLBFGSState &state,int &n,int &m,
                                   int &i,int &j,int &ic,int &mcinfo,
                                   double &v,double &vv)
  {
//--- check
   if(i>n-1)
     {
      //--- change values
      state.m_f=state.m_fbase;
      state.m_needf=false;
      //--- function call
      COptServ::TrimFunction(state.m_f,state.m_g,n,state.m_trimthreshold);
      //--- function call
      CLinMin::MCSrch(n,state.m_x,state.m_f,state.m_g,state.m_d,state.m_stp,state.m_stpmax,m_gtol,mcinfo,state.m_nfev,state.m_work,state.m_lstate,state.m_mcstage);
      //--- function call, return result
      return(Func_lbl_23(state,n,m,i,j,ic,mcinfo,v,vv));
     }
//--- change values
   v=state.m_x[i];
   state.m_x[i]=v-state.m_diffstep*state.m_s[i];
   state.m_rstate.stage=9;
//--- Saving state
   Func_lbl_rcomm(state,n,m,i,j,ic,mcinfo,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLBFGSIteration. Is a product to get    |
//| rid of the operator unconditional jump goto.                     |
//+------------------------------------------------------------------+
static bool CMinLBFGS::Func_lbl_30(CMinLBFGSState &state,int &n,int &m,
                                   int &i,int &j,int &ic,int &mcinfo,
                                   double &v,double &vv)
  {
//--- change values
   state.m_repnfev=state.m_repnfev+state.m_nfev;
   state.m_repiterationscount=state.m_repiterationscount+1;
//--- calculation
   for(int i_=0;i_<=n-1;i_++)
      state.m_sk[state.m_p].Set(i_,state.m_sk[state.m_p][i_]+state.m_x[i_]);
   for(int i_=0;i_<=n-1;i_++)
      state.m_yk[state.m_p].Set(i_,state.m_yk[state.m_p][i_]+state.m_g[i_]);
//--- Stopping conditions
   if(state.m_repiterationscount>=state.m_maxits&&state.m_maxits>0)
     {
      //--- Too many iterations
      state.m_repterminationtype=5;
      //--- return result
      return(false);
     }
//--- change value
   v=0;
   for(i=0;i<=n-1;i++)
      v=v+CMath::Sqr(state.m_g[i]*state.m_s[i]);
//--- check
   if(MathSqrt(v)<=state.m_epsg)
     {
      //--- Gradient is small enough
      state.m_repterminationtype=4;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_fold-state.m_f<=state.m_epsf*MathMax(MathAbs(state.m_fold),MathMax(MathAbs(state.m_f),1.0)))
     {
      //--- F(k+1)-F(k) is small enough
      state.m_repterminationtype=1;
      //--- return result
      return(false);
     }
//--- change value
   v=0;
   for(i=0;i<=n-1;i++)
      v=v+CMath::Sqr(state.m_sk[state.m_p][i]/state.m_s[i]);
//--- check
   if((double)(MathSqrt(v))<=(double)(state.m_epsx))
     {
      //--- X(k+1)-X(k) is small enough
      state.m_repterminationtype=2;
      //--- return result
      return(false);
     }
//--- If Wolfe conditions are satisfied, we can update
//--- limited memory model.
//--- However, if conditions are not satisfied (NFEV limit is met,
//--- function is too wild, ...), we'll skip L-BFGS update
   if(mcinfo!=1)
     {
      //--- Skip update.
      //--- In such cases we'll initialize search direction by
      //--- antigradient vector, because it  leads to more
      //--- transparent code with less number of special cases
      state.m_fold=state.m_f;
      for(int i_=0;i_<=n-1;i_++)
         state.m_d[i_]=-state.m_g[i_];
     }
   else
     {
      //--- Calculate Rho[k], GammaK
      v=0.0;
      for(int i_=0;i_<=n-1;i_++)
         v+=state.m_yk[state.m_p][i_]*state.m_sk[state.m_p][i_];
      //--- change value
      vv=0.0;
      for(int i_=0;i_<=n-1;i_++)
         vv+=state.m_yk[state.m_p][i_]*state.m_yk[state.m_p][i_];
      //--- check
      if(v==0.0 || vv==0.0)
        {
         //--- Rounding errors make further iterations impossible.
         state.m_repterminationtype=-2;
         //--- return result
         return(false);
        }
      //--- change values
      state.m_rho[state.m_p]=1/v;
      state.m_gammak=v/vv;
      //---  Calculate d(k+1)=-H(k+1)*g(k+1)
      //---  for I:=K downto K-Q do
      //---      V=s(i)^T * work(iteration:I)
      //---      theta(i)=V
      //---      work(iteration:I+1)=work(iteration:I)-V*Rho(i)*y(i)
      //---  work(last iteration)=H0*work(last iteration)-preconditioner
      //---  for I:=K-Q to K do
      //---      V=y(i)^T*work(iteration:I)
      //---      work(iteration:I+1)=work(iteration:I) +(-V+theta(i))*Rho(i)*s(i)
      //---  NOW WORK CONTAINS d(k+1)
      for(int i_=0;i_<=n-1;i_++)
         state.m_work[i_]=state.m_g[i_];
      for(i=state.m_k;i>=state.m_k-state.m_q;i--)
        {
         ic=i%m;
         v=0.0;
         for(int i_=0;i_<=n-1;i_++)
            v+=state.m_sk[ic][i_]*state.m_work[i_];
         //--- change values
         state.m_theta[ic]=v;
         vv=v*state.m_rho[ic];
         for(int i_=0;i_<=n-1;i_++)
            state.m_work[i_]=state.m_work[i_]-vv*state.m_yk[ic][i_];
        }
      //--- check
      if(state.m_prectype==0)
        {
         //--- Simple preconditioner is used
         v=state.m_gammak;
         for(int i_=0;i_<=n-1;i_++)
            state.m_work[i_]=v*state.m_work[i_];
        }
      //--- check
      if(state.m_prectype==1)
        {
         //--- Cholesky preconditioner is used
         CFbls::FblsCholeskySolve(state.m_denseh,1,n,true,state.m_work,state.m_autobuf);
        }
      //--- check
      if(state.m_prectype==2)
        {
         //--- diagonal approximation is used
         for(i=0;i<=n-1;i++)
           {
            state.m_work[i]=state.m_work[i]/state.m_diagh[i];
           }
        }
      //--- check
      if(state.m_prectype==3)
        {
         //--- scale-based preconditioner is used
         for(i=0;i<=n-1;i++)
            state.m_work[i]=state.m_work[i]*state.m_s[i]*state.m_s[i];
        }
      //--- calculation
      for(i=state.m_k-state.m_q;i<=state.m_k;i++)
        {
         ic=i%m;
         v=0.0;
         for(int i_=0;i_<=n-1;i_++)
            v+=state.m_yk[ic][i_]*state.m_work[i_];
         //--- change value
         vv=state.m_rho[ic]*(-v+state.m_theta[ic]);
         for(int i_=0;i_<=n-1;i_++)
           {
            state.m_work[i_]=state.m_work[i_]+vv*state.m_sk[ic][i_];
           }
        }
      for(int i_=0;i_<=n-1;i_++)
         state.m_d[i_]=-state.m_work[i_];
      //--- Next step
      state.m_fold=state.m_f;
      state.m_k=state.m_k+1;
     }
//--- function call, return result
   return(Func_lbl_21(state,n,m,i,j,ic,mcinfo,v,vv));
  }
//+------------------------------------------------------------------+
//| This object stores nonlinear optimizer state.                    |
//| You should use functions provided by MinQP subpackage to work    |
//| with this object                                                 |
//+------------------------------------------------------------------+
class CMinQPState
  {
public:
   //--- variables
   int               m_n;
   int               m_algokind;
   int               m_akind;
   bool              m_havex;
   double            m_constterm;
   int               m_repinneriterationscount;
   int               m_repouteriterationscount;
   int               m_repncholesky;
   int               m_repnmv;
   int               m_repterminationtype;
   CApBuff           m_buf;
   //--- arrays
   double            m_diaga[];
   double            m_b[];
   double            m_bndl[];
   double            m_bndu[];
   bool              m_havebndl[];
   bool              m_havebndu[];
   double            m_xorigin[];
   double            m_startx[];
   double            m_xc[];
   double            m_gc[];
   int               m_activeconstraints[];
   int               m_prevactiveconstraints[];
   double            m_workbndl[];
   double            m_workbndu[];
   double            m_tmp0[];
   double            m_tmp1[];
   int               m_itmp0[];
   int               m_p2[];
   double            m_bufb[];
   double            m_bufx[];
   //--- matrix
   CMatrixDouble     m_densea;
   CMatrixDouble     m_bufa;
   //--- constructor, destructor
                     CMinQPState(void);
                    ~CMinQPState(void);
   //--- copy
   void              Copy(CMinQPState &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinQPState::CMinQPState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinQPState::~CMinQPState(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMinQPState::Copy(CMinQPState &obj)
  {
//--- copy variables
   m_n=obj.m_n;
   m_algokind=obj.m_algokind;
   m_akind=obj.m_akind;
   m_havex=obj.m_havex;
   m_constterm=obj.m_constterm;
   m_repinneriterationscount=obj.m_repinneriterationscount;
   m_repouteriterationscount=obj.m_repouteriterationscount;
   m_repncholesky=obj.m_repncholesky;
   m_repnmv=obj.m_repnmv;
   m_repterminationtype=obj.m_repterminationtype;
   m_buf.Copy(obj.m_buf);
//--- copy arrays
   ArrayCopy(m_diaga,obj.m_diaga);
   ArrayCopy(m_b,obj.m_b);
   ArrayCopy(m_bndl,obj.m_bndl);
   ArrayCopy(m_bndu,obj.m_bndu);
   ArrayCopy(m_havebndl,obj.m_havebndl);
   ArrayCopy(m_havebndu,obj.m_havebndu);
   ArrayCopy(m_xorigin,obj.m_xorigin);
   ArrayCopy(m_startx,obj.m_startx);
   ArrayCopy(m_xc,obj.m_xc);
   ArrayCopy(m_gc,obj.m_gc);
   ArrayCopy(m_activeconstraints,obj.m_activeconstraints);
   ArrayCopy(m_prevactiveconstraints,obj.m_prevactiveconstraints);
   ArrayCopy(m_workbndl,obj.m_workbndl);
   ArrayCopy(m_workbndu,obj.m_workbndu);
   ArrayCopy(m_tmp0,obj.m_tmp0);
   ArrayCopy(m_tmp1,obj.m_tmp1);
   ArrayCopy(m_itmp0,obj.m_itmp0);
   ArrayCopy(m_p2,obj.m_p2);
   ArrayCopy(m_bufb,obj.m_bufb);
   ArrayCopy(m_bufx,obj.m_bufx);
//--- copy matrix
   m_densea=obj.m_densea;
   m_bufa=obj.m_bufa;
  }
//+------------------------------------------------------------------+
//| This object stores nonlinear optimizer state.                    |
//| You should use functions provided by MinQP subpackage to work    |
//| with this object                                                 |
//+------------------------------------------------------------------+
class CMinQPStateShell
  {
private:
   CMinQPState       m_innerobj;
public:
   //--- constructors, destructor
                     CMinQPStateShell(void);
                     CMinQPStateShell(CMinQPState &obj);
                    ~CMinQPStateShell(void);
   //--- method
   CMinQPState      *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinQPStateShell::CMinQPStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMinQPStateShell::CMinQPStateShell(CMinQPState &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinQPStateShell::~CMinQPStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMinQPState *CMinQPStateShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| This structure stores optimization report:                       |
//| * InnerIterationsCount      number of inner iterations           |
//| * OuterIterationsCount      number of outer iterations           |
//| * NCholesky                 number of Cholesky decomposition     |
//| * NMV                       number of matrix-vector products     |
//|                             (only products calculated as part of |
//|                             iterative process are counted)       |
//| * TerminationType           completion code (see below)          |
//| Completion codes:                                                |
//| * -5    inappropriate solver was used:                           |
//|         * Cholesky solver for semidefinite or indefinite problems|
//|         * Cholesky solver for problems with non-boundary         |
//|           constraints                                            |
//| * -3    inconsistent constraints (or, maybe, feasible point is   |
//|         too hard to find). If you are sure that constraints are  |
//|         feasible, try to restart optimizer with better initial   |
//|         approximation.                                           |
//| *  4    successful completion                                    |
//| *  5    MaxIts steps was taken                                   |
//| *  7    stopping conditions are too stringent,                   |
//|         further improvement is impossible,                       |
//|         X contains best point found so far.                      |
//+------------------------------------------------------------------+
class CMinQPReport
  {
public:
   //--- variables
   int               m_inneriterationscount;
   int               m_outeriterationscount;
   int               m_nmv;
   int               m_ncholesky;
   int               m_terminationtype;
   //--- constructor, destructor
                     CMinQPReport(void);
                    ~CMinQPReport(void);
   //--- copy
   void              Copy(CMinQPReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinQPReport::CMinQPReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinQPReport::~CMinQPReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMinQPReport::Copy(CMinQPReport &obj)
  {
//--- copy variables
   m_inneriterationscount=obj.m_inneriterationscount;
   m_outeriterationscount=obj.m_outeriterationscount;
   m_nmv=obj.m_nmv;
   m_ncholesky=obj.m_ncholesky;
   m_terminationtype=obj.m_terminationtype;
  }
//+------------------------------------------------------------------+
//| This structure stores optimization report:                       |
//| * InnerIterationsCount      number of inner iterations           |
//| * OuterIterationsCount      number of outer iterations           |
//| * NCholesky                 number of Cholesky decomposition     |
//| * NMV                       number of matrix-vector products     |
//|                             (only products calculated as part of |
//|                             iterative process are counted)       |
//| * TerminationType           completion code (see below)          |
//| Completion codes:                                                |
//| * -5    inappropriate solver was used:                           |
//|         * Cholesky solver for semidefinite or indefinite problems|
//|         * Cholesky solver for problems with non-boundary         |
//|           constraints                                            |
//| * -3    inconsistent constraints (or, maybe, feasible point is   |
//|         too hard to find). If you are sure that constraints are  |
//|         feasible, try to restart optimizer with better initial   |
//|         approximation.                                           |
//| *  4    successful completion                                    |
//| *  5    MaxIts steps was taken                                   |
//| *  7    stopping conditions are too stringent,                   |
//|         further improvement is impossible,                       |
//|         X contains best point found so far.                      |
//+------------------------------------------------------------------+
class CMinQPReportShell
  {
private:
   CMinQPReport      m_innerobj;
public:
   //--- constructors, destructor
                     CMinQPReportShell(void);
                     CMinQPReportShell(CMinQPReport &obj);
                    ~CMinQPReportShell(void);
   //--- methods
   int               GetInnerIterationsCount(void);
   void              SetInnerIterationsCount(const int i);
   int               GetOuterIterationsCount(void);
   void              SetOuterIterationsCount(const int i);
   int               GetNMV(void);
   void              SetNMV(const int i);
   int               GetNCholesky(void);
   void              SetNCholesky(const int i);
   int               GetTerminationType(void);
   void              SetTerminationType(const int i);
   CMinQPReport     *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinQPReportShell::CMinQPReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMinQPReportShell::CMinQPReportShell(CMinQPReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinQPReportShell::~CMinQPReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable inneriterationscount           |
//+------------------------------------------------------------------+
int CMinQPReportShell::GetInnerIterationsCount(void)
  {
//--- return result
   return(m_innerobj.m_inneriterationscount);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable inneriterationscount          |
//+------------------------------------------------------------------+
void CMinQPReportShell::SetInnerIterationsCount(const int i)
  {
//--- change value
   m_innerobj.m_inneriterationscount=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable outeriterationscount           |
//+------------------------------------------------------------------+
int CMinQPReportShell::GetOuterIterationsCount(void)
  {
//--- return result
   return(m_innerobj.m_outeriterationscount);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable outeriterationscount          |
//+------------------------------------------------------------------+
void CMinQPReportShell::SetOuterIterationsCount(const int i)
  {
//--- change value
   m_innerobj.m_outeriterationscount=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable nmv                            |
//+------------------------------------------------------------------+
int CMinQPReportShell::GetNMV(void)
  {
//--- return result
   return(m_innerobj.m_nmv);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable nmv                           |
//+------------------------------------------------------------------+
void CMinQPReportShell::SetNMV(const int i)
  {
//--- change value
   m_innerobj.m_nmv=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable ncholesky                      |
//+------------------------------------------------------------------+
int CMinQPReportShell::GetNCholesky(void)
  {
//--- return result
   return(m_innerobj.m_ncholesky);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable ncholesky                     |
//+------------------------------------------------------------------+
void CMinQPReportShell::SetNCholesky(const int i)
  {
//--- change value
   m_innerobj.m_ncholesky=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable terminationtype                |
//+------------------------------------------------------------------+
int CMinQPReportShell::GetTerminationType(void)
  {
//--- return result
   return(m_innerobj.m_terminationtype);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable terminationtype               |
//+------------------------------------------------------------------+
void CMinQPReportShell::SetTerminationType(const int i)
  {
//--- change value
   m_innerobj.m_terminationtype=i;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMinQPReport *CMinQPReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Constrained quadratic programming                                |
//+------------------------------------------------------------------+
class CMinQP
  {
private:
   //--- private methods
   static void       MinQPGrad(CMinQPState &state);
   static double     MinQPXTAX(CMinQPState &state,double &x[]);
public:
   //--- constructor, destructor
                     CMinQP(void);
                    ~CMinQP(void);
   //--- public methods
   static void       MinQPCreate(const int n,CMinQPState &state);
   static void       MinQPSetLinearTerm(CMinQPState &state,double &b[]);
   static void       MinQPSetQuadraticTerm(CMinQPState &state,CMatrixDouble &a,const bool isupper);
   static void       MinQPSetStartingPoint(CMinQPState &state,double &x[]);
   static void       MinQPSetOrigin(CMinQPState &state,double &xorigin[]);
   static void       MinQPSetAlgoCholesky(CMinQPState &state);
   static void       MinQPSetBC(CMinQPState &state,double &bndl[],double &bndu[]);
   static void       MinQPOptimize(CMinQPState &state);
   static void       MinQPResults(CMinQPState &state,double &x[],CMinQPReport &rep);
   static void       MinQPResultsBuf(CMinQPState &state,double &x[],CMinQPReport &rep);
   static void       MinQPSetLinearTermFast(CMinQPState &state,double &b[]);
   static void       MinQPSetQuadraticTermFast(CMinQPState &state,CMatrixDouble &a,const bool isupper,const double s);
   static void       MinQPRewriteDiagonal(CMinQPState &state,double &s[]);
   static void       MinQPSetStartingPointFast(CMinQPState &state,double &x[]);
   static void       MinQPSetOriginFast(CMinQPState &state,double &xorigin[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinQP::CMinQP(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinQP::~CMinQP(void)
  {

  }
//+------------------------------------------------------------------+
//|                     CONSTRAINED QUADRATIC PROGRAMMING            |
//| The subroutine creates QP optimizer. After initial creation, it  |
//| contains default optimization problem with zero quadratic and    |
//| linear terms and no constraints. You should set quadratic/linear |
//| terms with calls to functions provided by MinQP subpackage.      |
//| INPUT PARAMETERS:                                                |
//|     N       -   problem size                                     |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   optimizer with zero quadratic/linear terms       |
//|                 and no constraints                               |
//+------------------------------------------------------------------+
static void CMinQP::MinQPCreate(const int n,CMinQPState &state)
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1"))
      return;
//--- initialize QP solver
   state.m_n=n;
   state.m_akind=-1;
   state.m_repterminationtype=0;
//--- allocation
   ArrayResizeAL(state.m_b,n);
   ArrayResizeAL(state.m_bndl,n);
   ArrayResizeAL(state.m_bndu,n);
   ArrayResizeAL(state.m_workbndl,n);
   ArrayResizeAL(state.m_workbndu,n);
   ArrayResizeAL(state.m_havebndl,n);
   ArrayResizeAL(state.m_havebndu,n);
   ArrayResizeAL(state.m_startx,n);
   ArrayResizeAL(state.m_xorigin,n);
   ArrayResizeAL(state.m_xc,n);
   ArrayResizeAL(state.m_gc,n);
//--- initialization
   for(i=0;i<=n-1;i++)
     {
      state.m_b[i]=0.0;
      state.m_workbndl[i]=CInfOrNaN::NegativeInfinity();
      state.m_workbndu[i]=CInfOrNaN::PositiveInfinity();
      state.m_havebndl[i]=false;
      state.m_havebndu[i]=false;
      state.m_startx[i]=0.0;
      state.m_xorigin[i]=0.0;
     }
   state.m_havex=false;
//--- function call
   MinQPSetAlgoCholesky(state);
  }
//+------------------------------------------------------------------+
//| This function sets linear term for QP solver.                    |
//| By default, linear term is zero.                                 |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     B       -   linear term, array[N].                           |
//+------------------------------------------------------------------+
static void CMinQP::MinQPSetLinearTerm(CMinQPState &state,double &b[])
  {
//--- create a variable
   int n=0;
//--- initialization
   n=state.m_n;
//--- check
   if(!CAp::Assert(CAp::Len(b)>=n,__FUNCTION__+": Length(B)<N"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(b,n),__FUNCTION__+": B contains infinite or NaN elements"))
      return;
//--- function call
   MinQPSetLinearTermFast(state,b);
  }
//+------------------------------------------------------------------+
//| This function sets quadratic term for QP solver.                 |
//| By default quadratic term is zero.                               |
//| IMPORTANT: this solver minimizes following  function:            |
//|     f(x) = 0.5*x'*A*x + b'*x.                                    |
//| Note that quadratic term has 0.5 before it. So if you want to    |
//| minimize                                                         |
//|     f(x) = x^2 + x                                               |
//| you should rewrite your problem as follows:                      |
//|     f(x) = 0.5*(2*x^2) + x                                       |
//| and your matrix A will be equal to [[2.0]], not to [[1.0]]       |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     A       -   matrix, array[N,N]                               |
//|     IsUpper -   (optional) storage type:                         |
//|                 * if True, symmetric matrix A is given by its    |
//|                   upper triangle, and the lower triangle isn?t   |
//|                   used                                           |
//|                 * if False, symmetric matrix A is given by its   |
//|                   lower triangle, and the upper triangle isn?t   |
//|                   used                                           |
//|                 * if not given, both lower and upper triangles   |
//|                   must be filled.                                |
//+------------------------------------------------------------------+
static void CMinQP::MinQPSetQuadraticTerm(CMinQPState &state,CMatrixDouble &a,
                                          const bool isupper)
  {
//--- create a variable
   int n=0;
//--- initialization
   n=state.m_n;
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": Rows(A)<N"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": Cols(A)<N"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteRTrMatrix(a,n,isupper),__FUNCTION__+": A contains infinite or NaN elements"))
      return;
//--- function call
   MinQPSetQuadraticTermFast(state,a,isupper,0.0);
  }
//+------------------------------------------------------------------+
//| This function sets starting point for QP solver. It is useful to |
//| have good initial approximation to the solution, because it will |
//| increase speed of convergence and identification of active       |
//| constraints.                                                     |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     X       -   starting point, array[N].                        |
//+------------------------------------------------------------------+
static void CMinQP::MinQPSetStartingPoint(CMinQPState &state,double &x[])
  {
//--- create a variable
   int n=0;
//--- initialization
   n=state.m_n;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(B)<N"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN elements"))
      return;
//--- function call
   MinQPSetStartingPointFast(state,x);
  }
//+------------------------------------------------------------------+
//| This  function sets origin for QP solver. By default, following  |
//| QP program is solved:                                            |
//|     min(0.5*x'*A*x+b'*x)                                         |
//| This function allows to solve different problem:                 |
//|     min(0.5*(x-x_origin)'*A*(x-x_origin)+b'*(x-x_origin))        |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     XOrigin -   origin, array[N].                                |
//+------------------------------------------------------------------+
static void CMinQP::MinQPSetOrigin(CMinQPState &state,double &xorigin[])
  {
//--- create a variable
   int n=0;
//--- initialization
   n=state.m_n;
//--- check
   if(!CAp::Assert(CAp::Len(xorigin)>=n,__FUNCTION__+": Length(B)<N"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(xorigin,n),__FUNCTION__+": B contains infinite or NaN elements"))
      return;
//--- function call
   MinQPSetOriginFast(state,xorigin);
  }
//+------------------------------------------------------------------+
//| This function tells solver to use Cholesky-based algorithm.      |
//| Cholesky-based algorithm can be used when:                       |
//| * problem is convex                                              |
//| * there is no constraints or only boundary constraints are       |
//|   present                                                        |
//| This algorithm has O(N^3) complexity for unconstrained problem   |
//| and is up to several times slower on bound constrained problems  |
//| (these additional iterations are needed to identify active       |
//| constraints).                                                    |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//+------------------------------------------------------------------+
static void CMinQP::MinQPSetAlgoCholesky(CMinQPState &state)
  {
//--- change value
   state.m_algokind=1;
  }
//+------------------------------------------------------------------+
//| This function sets boundary constraints for QP solver            |
//| Boundary constraints are inactive by default (after initial      |
//| creation). After being set, they are preserved until explicitly  |
//| turned off with another SetBC() call.                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure stores algorithm state                 |
//|     BndL    -   lower bounds, array[N].                          |
//|                 If some (all) variables are unbounded, you may   |
//|                 specify very small number or -INF (latter is     |
//|                 recommended because it will allow solver to use  |
//|                 better algorithm).                               |
//|     BndU    -   upper bounds, array[N].                          |
//|                 If some (all) variables are unbounded, you may   |
//|                 specify very large number or +INF (latter is     |
//|                 recommended because it will allow solver to use  |
//|                 better algorithm).                               |
//| NOTE: it is possible to specify BndL[i]=BndU[i]. In this case    |
//| I-th variable will be "frozen" at X[i]=BndL[i]=BndU[i].          |
//+------------------------------------------------------------------+
static void CMinQP::MinQPSetBC(CMinQPState &state,double &bndl[],double &bndu[])
  {
//--- create variables
   int i=0;
   int n=0;
//--- initialization
   n=state.m_n;
//--- check
   if(!CAp::Assert(CAp::Len(bndl)>=n,__FUNCTION__+": Length(BndL)<N"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(bndu)>=n,__FUNCTION__+": Length(BndU)<N"))
      return;
//--- change values
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(bndl[i]) || CInfOrNaN::IsNegativeInfinity(bndl[i]),__FUNCTION__+": BndL contains NAN or +INF"))
         return;
      //--- check
      if(!CAp::Assert(CMath::IsFinite(bndu[i]) || CInfOrNaN::IsPositiveInfinity(bndu[i]),__FUNCTION__+": BndU contains NAN or -INF"))
         return;
      //--- change values
      state.m_bndl[i]=bndl[i];
      state.m_havebndl[i]=CMath::IsFinite(bndl[i]);
      state.m_bndu[i]=bndu[i];
      state.m_havebndu[i]=CMath::IsFinite(bndu[i]);
     }
  }
//+------------------------------------------------------------------+
//| This function solves quadratic programming problem.              |
//| You should call it after setting solver options with             |
//| MinQPSet...() calls.                                             |
//| INPUT PARAMETERS:                                                |
//|     State   -   algorithm state                                  |
//| You should use MinQPResults() function to access results after   |
//| calls to this function.                                          |
//+------------------------------------------------------------------+
static void CMinQP::MinQPOptimize(CMinQPState &state)
  {
//--- create variables
   int    n=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    nbc=0;
   int    nlc=0;
   int    nactive=0;
   int    nfree=0;
   double f=0;
   double fprev=0;
   double v=0;
   bool   b;
   int    i_=0;
//--- initialization
   n=state.m_n;
   state.m_repterminationtype=-5;
   state.m_repinneriterationscount=0;
   state.m_repouteriterationscount=0;
   state.m_repncholesky=0;
   state.m_repnmv=0;
//--- check correctness of constraints
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(state.m_havebndl[i] && state.m_havebndu[i])
        {
         //--- check
         if(state.m_bndl[i]>state.m_bndu[i])
           {
            state.m_repterminationtype=-3;
            return;
           }
        }
     }
//--- count number of bound and linear constraints
   nbc=0;
   nlc=0;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(state.m_havebndl[i])
         nbc=nbc+1;
      //--- check
      if(state.m_havebndu[i])
         nbc=nbc+1;
     }
//--- Our formulation of quadratic problem includes origin point,
//--- i.m_e. we have F(x-x_origin) which is minimized subject to
//--- constraints on x,instead of having simply F(x).
//--- Here we make transition from non-zero origin to zero one.
//--- In order to make such transition we have to:
//--- 1. subtract x_origin from x_start
//--- 2. modify constraints
//--- 3. solve problem
//--- 4. add x_origin to solution
//--- There is alternate solution - to modify quadratic function
//--- by expansion of multipliers containing (x-x_origin),but
//--- we prefer to modify constraints,because it is a) more precise
//--- and b) easier to to.
//--- Parts (1)-(2) are done here. After this block is over,
//--- we have:
//--- * XC,which stores shifted XStart (if we don't have XStart,
//---   value of XC will be ignored later)
//--- * WorkBndL,WorkBndU,which store modified boundary constraints.
   for(i=0;i<=n-1;i++)
     {
      state.m_xc[i]=state.m_startx[i]-state.m_xorigin[i];
      //--- check
      if(state.m_havebndl[i])
         state.m_workbndl[i]=state.m_bndl[i]-state.m_xorigin[i];
      //--- check
      if(state.m_havebndu[i])
         state.m_workbndu[i]=state.m_bndu[i]-state.m_xorigin[i];
     }
//--- modify starting point XC according to boundary constraints
   if(state.m_havex)
     {
      //--- We have starting point in XC,so we just have to bound it
      for(i=0;i<=n-1;i++)
        {
         //--- check
         if(state.m_havebndl[i])
           {
            //--- check
            if(state.m_xc[i]<state.m_workbndl[i])
               state.m_xc[i]=state.m_workbndl[i];
           }
         //--- check
         if(state.m_havebndu[i])
           {
            //--- check
            if(state.m_xc[i]>state.m_workbndu[i])
               state.m_xc[i]=state.m_workbndu[i];
           }
        }
     }
   else
     {
      //--- We don't have starting point,so we deduce it from
      //--- constraints (if they are present).
      //--- NOTE: XC contains some meaningless values from previous block
      //--- which are ignored by code below.
      for(i=0;i<=n-1;i++)
        {
         //--- check
         if(state.m_havebndl[i] && state.m_havebndu[i])
           {
            state.m_xc[i]=0.5*(state.m_workbndl[i]+state.m_workbndu[i]);
            //--- check
            if(state.m_xc[i]<state.m_workbndl[i])
               state.m_xc[i]=state.m_workbndl[i];
            //--- check
            if(state.m_xc[i]>state.m_workbndu[i])
               state.m_xc[i]=state.m_workbndu[i];
            //--- continue iteration
            continue;
           }
         //--- check
         if(state.m_havebndl[i])
           {
            state.m_xc[i]=state.m_workbndl[i];
            continue;
           }
         //--- check
         if(state.m_havebndu[i])
           {
            state.m_xc[i]=state.m_workbndu[i];
            continue;
           }
         state.m_xc[i]=0;
        }
     }
//--- Select algo
   if(state.m_algokind==1 && state.m_akind==0)
     {
      //--- Cholesky-based algorithm for dense bound constrained problems.
      //--- This algorithm exists in two variants:
      //--- * unconstrained one,which can solve problem using only one NxN
      //---   double matrix
      //--- * bound constrained one,which needs two NxN matrices
      //--- We will try to solve problem using unconstrained algorithm,
      //--- and will use bound constrained version only when constraints
      //--- are actually present
      if(nbc==0 && nlc==0)
        {
         //--- "Simple" unconstrained version
         CApServ::RVectorSetLengthAtLeast(state.m_tmp0,n);
         //--- function call
         CApServ::RVectorSetLengthAtLeast(state.m_bufb,n);
         //--- calculation
         state.m_densea[0].Set(0,state.m_diaga[0]);
         for(k=1;k<=n-1;k++)
           {
            for(i_=0;i_<=k-1;i_++)
               state.m_densea[i_].Set(k,state.m_densea[k][i_]);
            state.m_densea[k].Set(k,state.m_diaga[k]);
           }
         //--- change values
         for(i_=0;i_<=n-1;i_++)
            state.m_bufb[i_]=state.m_b[i_];
         state.m_repncholesky=1;
         //--- check
         if(!CTrFac::SPDMatrixCholeskyRec(state.m_densea,0,n,true,state.m_tmp0))
           {
            state.m_repterminationtype=-5;
            return;
           }
         //--- function call
         CFbls::FblsCholeskySolve(state.m_densea,1.0,n,true,state.m_bufb,state.m_tmp0);
         for(i_=0;i_<=n-1;i_++)
            state.m_xc[i_]=-state.m_bufb[i_];
         for(i_=0;i_<=n-1;i_++)
            state.m_xc[i_]=state.m_xc[i_]+state.m_xorigin[i_];
         //--- change values
         state.m_repouteriterationscount=1;
         state.m_repterminationtype=4;
         //--- exit the function
         return;
        }
      //--- General bound constrained algo
      CApServ::RMatrixSetLengthAtLeast(state.m_bufa,n,n);
      //--- function call
      CApServ::RVectorSetLengthAtLeast(state.m_bufb,n);
      //--- function call
      CApServ::RVectorSetLengthAtLeast(state.m_bufx,n);
      //--- function call
      CApServ::IVectorSetLengthAtLeast(state.m_activeconstraints,n);
      //--- function call
      CApServ::IVectorSetLengthAtLeast(state.m_prevactiveconstraints,n);
      //--- function call
      CApServ::RVectorSetLengthAtLeast(state.m_tmp0,n);
      //--- Prepare constraints vectors:
      //--- * ActiveConstraints - constraints active at current step
      //--- * PrevActiveConstraints - constraints which were active at previous step
      //--- Elements of constraints vectors can be:
      //--- *  0 - inactive
      //--- *  1 - active
      //--- * -1 - undefined (used to initialize PrevActiveConstraints before first iteration)
      for(i=0;i<=n-1;i++)
         state.m_prevactiveconstraints[i]=-1;
      //--- Main cycle
      fprev=CMath::m_maxrealnumber;
      while(true)
        {
         //--- * calculate gradient at XC
         //--- * determine active constraints
         //--- * break if there is no free variables or
         //---   there were no changes in the list of active constraints
         MinQPGrad(state);
         nactive=0;
         for(i=0;i<=n-1;i++)
           {
            state.m_activeconstraints[i]=0;
            //--- check
            if(state.m_havebndl[i])
              {
               //--- check
               if(state.m_xc[i]<=state.m_workbndl[i] && state.m_gc[i]>=0.0)
                  state.m_activeconstraints[i]=1;
              }
            //--- check
            if(state.m_havebndu[i])
              {
               //--- check
               if(state.m_xc[i]>=state.m_workbndu[i] && state.m_gc[i]<=0.0)
                  state.m_activeconstraints[i]=1;
              }
            //--- check
            if(state.m_havebndl[i] && state.m_havebndu[i])
              {
               //--- check
               if(state.m_workbndl[i]==state.m_workbndu[i])
                  state.m_activeconstraints[i]=1;
              }
            //--- check
            if(state.m_activeconstraints[i]>0)
               nactive=nactive+1;
           }
         nfree=n-nactive;
         //--- check
         if(nfree==0)
            break;
         b=false;
         for(i=0;i<=n-1;i++)
           {
            //--- check
            if(state.m_activeconstraints[i]!=state.m_prevactiveconstraints[i])
               b=true;
           }
         //--- check
         if(!b)
            break;
         //--- * copy A,B and X to buffer
         //--- * rearrange BufA,BufB and BufX,in such way that active variables come first,
         //---   inactive are moved to the tail. We use sorting subroutine
         //---   to solve this problem.
         state.m_bufa[0].Set(0,state.m_diaga[0]);
         for(k=1;k<=n-1;k++)
           {
            for(i_=0;i_<=k-1;i_++)
               state.m_bufa[k].Set(i_,state.m_densea[k][i_]);
            for(i_=0;i_<=k-1;i_++)
               state.m_bufa[i_].Set(k,state.m_densea[k][i_]);
            state.m_bufa[k].Set(k,state.m_diaga[k]);
           }
         //--- change values
         for(i_=0;i_<=n-1;i_++)
            state.m_bufb[i_]=state.m_b[i_];
         for(i_=0;i_<=n-1;i_++)
            state.m_bufx[i_]=state.m_xc[i_];
         for(i=0;i<=n-1;i++)
            state.m_tmp0[i]=state.m_activeconstraints[i];
         //--- function call
         CTSort::TagSortBuf(state.m_tmp0,n,state.m_itmp0,state.m_p2,state.m_buf);
         for(k=0;k<=n-1;k++)
           {
            //--- check
            if(state.m_p2[k]!=k)
              {
               //--- swap
               v=state.m_bufb[k];
               state.m_bufb[k]=state.m_bufb[state.m_p2[k]];
               state.m_bufb[state.m_p2[k]]=v;
               v=state.m_bufx[k];
               state.m_bufx[k]=state.m_bufx[state.m_p2[k]];
               state.m_bufx[state.m_p2[k]]=v;
              }
           }
         for(i=0;i<=n-1;i++)
           {
            for(i_=0;i_<=n-1;i_++)
               state.m_tmp0[i_]=state.m_bufa[i][i_];
            for(k=0;k<=n-1;k++)
              {
               //--- check
               if(state.m_p2[k]!=k)
                 {
                  //--- swap
                  v=state.m_tmp0[k];
                  state.m_tmp0[k]=state.m_tmp0[state.m_p2[k]];
                  state.m_tmp0[state.m_p2[k]]=v;
                 }
              }
            for(i_=0;i_<=n-1;i_++)
               state.m_bufa[i].Set(i_,state.m_tmp0[i_]);
           }
         for(i=0;i<=n-1;i++)
           {
            //--- check
            if(state.m_p2[i]!=i)
              {
               for(i_=0;i_<=n-1;i_++)
                  state.m_tmp0[i_]=state.m_bufa[i][i_];
               for(i_=0;i_<=n-1;i_++)
                  state.m_bufa[i].Set(i_,state.m_bufa[state.m_p2[i]][i_]);
               for(i_=0;i_<=n-1;i_++)
                  state.m_bufa[state.m_p2[i]].Set(i_,state.m_tmp0[i_]);
              }
           }
         //--- Now we have A and B in BufA and BufB,variables are rearranged
         //--- into two groups: Xf - free variables,Xc - active (fixed) variables,
         //--- and our quadratic problem can be written as
         //---                         ( Af  Ac  )   ( Xf )                 ( Xf )
         //--- F(X)=0.5* ( Xf' Xc' ) * (         ) * (    ) + ( Bf' Bc' ) * (    )
         //---                         ( Ac' Acc )   ( Xc )                 ( Xc )
         //--- we want to convert to the optimization with respect to Xf,
         //--- treating Xc as constant term. After expansion of expression above
         //--- we get
         //--- F(Xf)=0.5*Xf'*Af*Xf + (Bf+Ac*Xc)'*Xf + 0.5*Xc'*Acc*Xc
         //--- We will update BufB using this expression and calculate
         //--- constant term.
         CAblas::RMatrixMVect(nfree,nactive,state.m_bufa,0,nfree,0,state.m_bufx,nfree,state.m_tmp0,0);
         for(i_=0;i_<=nfree-1;i_++)
            state.m_bufb[i_]=state.m_bufb[i_]+state.m_tmp0[i_];
         state.m_constterm=0.0;
         for(i=nfree;i<=n-1;i++)
           {
            state.m_constterm=state.m_constterm+0.5*state.m_bufx[i]*state.m_bufa[i][i]*state.m_bufx[i];
            for(j=i+1;j<=n-1;j++)
               state.m_constterm=state.m_constterm+state.m_bufx[i]*state.m_bufa[i][j]*state.m_bufx[j];
           }
         //--- Now we are ready to minimize F(Xf)...
         state.m_repncholesky=state.m_repncholesky+1;
         //--- check
         if(!CTrFac::SPDMatrixCholeskyRec(state.m_bufa,0,nfree,true,state.m_tmp0))
           {
            state.m_repterminationtype=-5;
            return;
           }
         //--- function call
         CFbls::FblsCholeskySolve(state.m_bufa,1.0,nfree,true,state.m_bufb,state.m_tmp0);
         for(i_=0;i_<=nfree-1;i_++)
            state.m_bufx[i_]=-state.m_bufb[i_];
         //--- ...m_and to copy results back to XC.
         //--- It is done in several steps:
         //--- * original order of variables is restored
         //--- * result is copied back to XC
         //--- * XC is bounded with respect to bound constraints
         for(k=n-1;k>=0;k--)
           {
            //--- check
            if(state.m_p2[k]!=k)
              {
               v=state.m_bufx[k];
               state.m_bufx[k]=state.m_bufx[state.m_p2[k]];
               state.m_bufx[state.m_p2[k]]=v;
              }
           }
         for(i_=0;i_<=n-1;i_++)
            state.m_xc[i_]=state.m_bufx[i_];
         for(i=0;i<=n-1;i++)
           {
            //--- check
            if(state.m_havebndl[i])
              {
               //--- check
               if(state.m_xc[i]<state.m_workbndl[i])
                  state.m_xc[i]=state.m_workbndl[i];
              }
            //--- check
            if(state.m_havebndu[i])
              {
               //--- check
               if(state.m_xc[i]>state.m_workbndu[i])
                  state.m_xc[i]=state.m_workbndu[i];
              }
           }
         //--- Calculate F,compare it with FPrev.
         //--- Break if F>=FPrev
         //--- (sometimes possible at extremum due to numerical noise).
         f=0.0;
         for(i_=0;i_<=n-1;i_++)
            f+=state.m_b[i_]*state.m_xc[i_];
         f=f+MinQPXTAX(state,state.m_xc);
         //--- check
         if(f>=fprev)
            break;
         fprev=f;
         //--- Update PrevActiveConstraints
         for(i=0;i<=n-1;i++)
            state.m_prevactiveconstraints[i]=state.m_activeconstraints[i];
         //--- Update report-related fields
         state.m_repouteriterationscount=state.m_repouteriterationscount+1;
        }
      //--- change values
      state.m_repterminationtype=4;
      for(i_=0;i_<=n-1;i_++)
         state.m_xc[i_]=state.m_xc[i_]+state.m_xorigin[i_];
      //--- exit the function
      return;
     }
  }
//+------------------------------------------------------------------+
//| QP solver results                                                |
//| INPUT PARAMETERS:                                                |
//|     State   -   algorithm state                                  |
//| OUTPUT PARAMETERS:                                               |
//|     X       -   array[0..N-1], solution                          |
//|     Rep     -   optimization report. You should check Rep.       |
//|                 TerminationType, which contains completion code, |
//|                 and you may check another fields which contain   |
//|                 another information about algorithm functioning. |
//+------------------------------------------------------------------+
static void CMinQP::MinQPResults(CMinQPState &state,double &x[],CMinQPReport &rep)
  {
//--- reset memory
   ArrayResizeAL(x,0);
//--- function call
   MinQPResultsBuf(state,x,rep);
  }
//+------------------------------------------------------------------+
//| QP results                                                       |
//| Buffered implementation of MinQPResults() which uses             |
//| pre-allocated  buffer to store X[]. If buffer size is too small, |
//| it resizes buffer. It is intended to be used in the inner cycles |
//| of performance critical algorithms where array reallocation      |
//| penalty is too large to be ignored.                              |
//+------------------------------------------------------------------+
static void CMinQP::MinQPResultsBuf(CMinQPState &state,double &x[],
                                    CMinQPReport &rep)
  {
//--- create a variable
   int i_=0;
//--- check
   if(CAp::Len(x)<state.m_n)
      ArrayResizeAL(x,state.m_n);
//--- copy
   for(i_=0;i_<=state.m_n-1;i_++)
      x[i_]=state.m_xc[i_];
//--- change values
   rep.m_inneriterationscount=state.m_repinneriterationscount;
   rep.m_outeriterationscount=state.m_repouteriterationscount;
   rep.m_nmv=state.m_repnmv;
   rep.m_terminationtype=state.m_repterminationtype;
  }
//+------------------------------------------------------------------+
//| Fast version of MinQPSetLinearTerm(), which doesn't check its    |
//| arguments. For internal use only.                                |
//+------------------------------------------------------------------+
static void CMinQP::MinQPSetLinearTermFast(CMinQPState &state,double &b[])
  {
//--- create variables
   int n=0;
   int i_=0;
//--- initialization
   n=state.m_n;
   for(i_=0;i_<=n-1;i_++)
      state.m_b[i_]=b[i_];
  }
//+------------------------------------------------------------------+
//| Fast version of MinQPSetQuadraticTerm(), which doesn't check its |
//| arguments.                                                       |
//| It accepts additional parameter - shift S, which allows to       |
//| "shift" matrix A by adding s*I to A. S must be positive (although|
//| it is not checked).                                              |
//| For internal use only.                                           |
//+------------------------------------------------------------------+
static void CMinQP::MinQPSetQuadraticTermFast(CMinQPState &state,CMatrixDouble &a,
                                              const bool isupper,const double s)
  {
//--- create variables
   int k=0;
   int n=0;
   int i_=0;
//--- We store off-diagonal part of A in the lower triangle of DenseA.
//--- Diagonal elements of A are stored in the DiagA.
//--- Diagonal of DenseA and uppper triangle are used as temporaries.
//--- Why such complex storage? Because it:
//--- 1. allows us to easily recover from exceptions (lower triangle
//---    is unmodified during execution as well as DiagA,and on entry
//---    we will always find unmodified matrix)
//--- 2. allows us to make Cholesky decomposition in the upper triangle
//---    of DenseA or to do other SPD-related operations.
   n=state.m_n;
   state.m_akind=0;
//--- function call
   CApServ::RMatrixSetLengthAtLeast(state.m_densea,n,n);
//--- function call
   CApServ::RVectorSetLengthAtLeast(state.m_diaga,n);
//--- check
   if(isupper)
     {
      for(k=0;k<=n-2;k++)
        {
         //--- change values
         state.m_diaga[k]=a[k][k]+s;
         for(i_=k+1;i_<=n-1;i_++)
            state.m_densea[i_].Set(k,a[k][i_]);
        }
      state.m_diaga[n-1]=a[n-1][n-1]+s;
     }
   else
     {
      state.m_diaga[0]=a[0][0]+s;
      for(k=1;k<=n-1;k++)
        {
         //--- change values
         for(i_=0;i_<=k-1;i_++)
            state.m_densea[k].Set(i_,a[k][i_]);
         state.m_diaga[k]=a[k][k]+s;
        }
     }
  }
//+------------------------------------------------------------------+
//| Interna lfunction which allows to rewrite diagonal of quadratic  |
//| term. For internal use only.                                     |
//| This function can be used only when you have dense A and already |
//| made MinQPSetQuadraticTerm(Fast) call.                           |
//+------------------------------------------------------------------+
static void CMinQP::MinQPRewriteDiagonal(CMinQPState &state,double &s[])
  {
//--- create variables
   int k=0;
   int n=0;
//--- check
   if(!CAp::Assert(state.m_akind==0,__FUNCTION__+": internal error (AKind<>0)"))
      return;
//--- initialization
   n=state.m_n;
   for(k=0;k<=n-1;k++)
      state.m_diaga[k]=s[k];
  }
//+------------------------------------------------------------------+
//| Fast version of MinQPSetStartingPoint(), which doesn't check its |
//| arguments. For internal use only.                                |
//+------------------------------------------------------------------+
static void CMinQP::MinQPSetStartingPointFast(CMinQPState &state,double &x[])
  {
//--- create variables
   int n=0;
   int i_=0;
//--- initialization
   n=state.m_n;
   for(i_=0;i_<=n-1;i_++)
      state.m_startx[i_]=x[i_];
   state.m_havex=true;
  }
//+------------------------------------------------------------------+
//| Fast version of MinQPSetOrigin(), which doesn't check its        |
//| arguments. For internal use only.                                |
//+------------------------------------------------------------------+
static void CMinQP::MinQPSetOriginFast(CMinQPState &state,double &xorigin[])
  {
//--- create variables
   int n=0;
   int i_=0;
//--- initialization
   n=state.m_n;
   for(i_=0;i_<=n-1;i_++)
      state.m_xorigin[i_]=xorigin[i_];
  }
//+------------------------------------------------------------------+
//| This function calculates gradient of quadratic function at XC and|
//| stores it in the GC.                                             |
//+------------------------------------------------------------------+
static void CMinQP::MinQPGrad(CMinQPState &state)
  {
//--- create variables
   int    n=0;
   int    i=0;
   double v=0;
   int    i_=0;
//--- initialization
   n=state.m_n;
//--- check
   if(!CAp::Assert(state.m_akind==-1||state.m_akind==0,__FUNCTION__+": internal error"))
      return;
//--- zero A
   if(state.m_akind==-1)
     {
      for(i_=0;i_<=n-1;i_++)
         state.m_gc[i_]=state.m_b[i_];
      //--- exit the function
      return;
     }
//--- dense A
   if(state.m_akind==0)
     {
      for(i_=0;i_<=n-1;i_++)
         state.m_gc[i_]=state.m_b[i_];
      state.m_gc[0]=state.m_gc[0]+state.m_diaga[0]*state.m_xc[0];
      //--- calculation
      for(i=1;i<=n-1;i++)
        {
         v=0.0;
         for(i_=0;i_<=i-1;i_++)
            v+=state.m_densea[i][i_]*state.m_xc[i_];
         state.m_gc[i]=state.m_gc[i]+v+state.m_diaga[i]*state.m_xc[i];
         v=state.m_xc[i];
         //--- change values
         for(i_=0;i_<=i-1;i_++)
            state.m_gc[i_]=state.m_gc[i_]+v*state.m_densea[i][i_];
        }
      //--- exit the function
      return;
     }
  }
//+------------------------------------------------------------------+
//| This function calculates x'*A*x for given X.                     |
//+------------------------------------------------------------------+
static double CMinQP::MinQPXTAX(CMinQPState &state,double &x[])
  {
//--- create variables
   double result=0;
   int    n=0;
   int    i=0;
   int    j=0;
//--- initialization
   n=state.m_n;
//--- check
   if(!CAp::Assert(state.m_akind==-1 || state.m_akind==0,__FUNCTION__+": internal error"))
      return(EMPTY_VALUE);
   result=0;
//--- zero A
   if(state.m_akind==-1)
      return(0.0);
//--- dense A
   if(state.m_akind==0)
     {
      result=0;
      for(i=0;i<=n-1;i++)
        {
         for(j=0;j<=i-1;j++)
            result=result+state.m_densea[i][j]*x[i]*x[j];
         //--- get result
         result=result+0.5*state.m_diaga[i]*CMath::Sqr(x[i]);
        }
      //--- return result
      return(result);
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Levenberg-Marquardt optimizer.                                   |
//| This structure should be created using one of the                |
//| MinLMCreate() functions. You should not access its fields        |
//| directly; use ALGLIB functions to work with it.                  |
//+------------------------------------------------------------------+
class CMinLMState
  {
public:
   //--- variables
   int               m_n;
   int               m_m;
   double            m_diffstep;
   double            m_epsg;
   double            m_epsf;
   double            m_epsx;
   int               m_maxits;
   bool              m_xrep;
   double            m_stpmax;
   int               m_maxmodelage;
   bool              m_makeadditers;
   double            m_f;
   bool              m_needf;
   bool              m_needfg;
   bool              m_needfgh;
   bool              m_needfij;
   bool              m_needfi;
   bool              m_xupdated;
   int               m_algomode;
   bool              m_hasf;
   bool              m_hasfi;
   bool              m_hasg;
   double            m_fbase;
   double            m_lambdav;
   double            m_nu;
   int               m_modelage;
   bool              m_deltaxready;
   bool              m_deltafready;
   int               m_repiterationscount;
   int               m_repterminationtype;
   int               m_repnfunc;
   int               m_repnjac;
   int               m_repngrad;
   int               m_repnhess;
   int               m_repncholesky;
   RCommState        m_rstate;
   double            m_actualdecrease;
   double            m_predicteddecrease;
   double            m_xm1;
   double            m_xp1;
   CMinLBFGSState    m_internalstate;
   CMinLBFGSReport   m_internalrep;
   CMinQPState       m_qpstate;
   CMinQPReport      m_qprep;
   //--- arrays
   double            m_x[];
   double            m_fi[];
   double            m_g[];
   double            m_xbase[];
   double            m_fibase[];
   double            m_gbase[];
   double            m_bndl[];
   double            m_bndu[];
   bool              m_havebndl[];
   bool              m_havebndu[];
   double            m_s[];
   double            m_xdir[];
   double            m_deltax[];
   double            m_deltaf[];
   double            m_choleskybuf[];
   double            m_tmp0[];
   double            m_fm1[];
   double            m_fp1[];
   //--- matrix
   CMatrixDouble     m_j;
   CMatrixDouble     m_h;
   CMatrixDouble     m_quadraticmodel;
   //--- constructor, destructor
                     CMinLMState(void);
                    ~CMinLMState(void);
   //--- copy
   void              Copy(CMinLMState &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinLMState::CMinLMState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinLMState::~CMinLMState(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMinLMState::Copy(CMinLMState &obj)
  {
//--- copy variables
   m_n=obj.m_n;
   m_m=obj.m_m;
   m_diffstep=obj.m_diffstep;
   m_epsg=obj.m_epsg;
   m_epsf=obj.m_epsf;
   m_epsx=obj.m_epsx;
   m_maxits=obj.m_maxits;
   m_xrep=obj.m_xrep;
   m_stpmax=obj.m_stpmax;
   m_maxmodelage=obj.m_maxmodelage;
   m_makeadditers=obj.m_makeadditers;
   m_f=obj.m_f;
   m_needf=obj.m_needf;
   m_needfg=obj.m_needfg;
   m_needfgh=obj.m_needfgh;
   m_needfij=obj.m_needfij;
   m_needfi=obj.m_needfi;
   m_xupdated=obj.m_xupdated;
   m_algomode=obj.m_algomode;
   m_hasf=obj.m_hasf;
   m_hasfi=obj.m_hasfi;
   m_hasg=obj.m_hasg;
   m_fbase=obj.m_fbase;
   m_lambdav=obj.m_lambdav;
   m_nu=obj.m_nu;
   m_modelage=obj.m_modelage;
   m_deltaxready=obj.m_deltaxready;
   m_deltafready=obj.m_deltafready;
   m_repiterationscount=obj.m_repiterationscount;
   m_repterminationtype=obj.m_repterminationtype;
   m_repnfunc=obj.m_repnfunc;
   m_repnjac=obj.m_repnjac;
   m_repngrad=obj.m_repngrad;
   m_repnhess=obj.m_repnhess;
   m_repncholesky=obj.m_repncholesky;
   m_actualdecrease=obj.m_actualdecrease;
   m_predicteddecrease=obj.m_predicteddecrease;
   m_xm1=obj.m_xm1;
   m_xp1=obj.m_xp1;
   m_rstate.Copy(obj.m_rstate);
   m_internalstate.Copy(obj.m_internalstate);
   m_internalrep.Copy(obj.m_internalrep);
   m_qpstate.Copy(obj.m_qpstate);
   m_qprep.Copy(obj.m_qprep);
//--- copy arrays
   ArrayCopy(m_x,obj.m_x);
   ArrayCopy(m_fi,obj.m_fi);
   ArrayCopy(m_g,obj.m_g);
   ArrayCopy(m_xbase,obj.m_xbase);
   ArrayCopy(m_fibase,obj.m_fibase);
   ArrayCopy(m_gbase,obj.m_gbase);
   ArrayCopy(m_bndl,obj.m_bndl);
   ArrayCopy(m_bndu,obj.m_bndu);
   ArrayCopy(m_havebndl,obj.m_havebndl);
   ArrayCopy(m_havebndu,obj.m_havebndu);
   ArrayCopy(m_s,obj.m_s);
   ArrayCopy(m_xdir,obj.m_xdir);
   ArrayCopy(m_deltax,obj.m_deltax);
   ArrayCopy(m_deltaf,obj.m_deltaf);
   ArrayCopy(m_choleskybuf,obj.m_choleskybuf);
   ArrayCopy(m_tmp0,obj.m_tmp0);
   ArrayCopy(m_fm1,obj.m_fm1);
   ArrayCopy(m_fp1,obj.m_fp1);
//--- copy matrix
   m_j=obj.m_j;
   m_h=obj.m_h;
   m_quadraticmodel=obj.m_quadraticmodel;
  }
//+------------------------------------------------------------------+
//| Levenberg-Marquardt optimizer.                                   |
//| This structure should be created using one of the                |
//| MinLMCreate() functions. You should not access its fields        |
//| directly; use ALGLIB functions to work with it.                  |
//+------------------------------------------------------------------+
class CMinLMStateShell
  {
private:
   CMinLMState       m_innerobj;
public:
   //--- constructors, destructor
                     CMinLMStateShell(void);
                     CMinLMStateShell(CMinLMState &obj);
                    ~CMinLMStateShell(void);
   //--- methods
   bool              GetNeedF(void);
   void              SetNeedF(const bool b);
   bool              GetNeedFG(void);
   void              SetNeedFG(const bool b);
   bool              GetNeedFGH(void);
   void              SetNeedFGH(const bool b);
   bool              GetNeedFI(void);
   void              SetNeedFI(const bool b);
   bool              GetNeedFIJ(void);
   void              SetNeedFIJ(const bool b);
   bool              GetXUpdated(void);
   void              SetXUpdated(const bool b);
   double            GetF(void);
   void              SetF(const double d);
   CMinLMState      *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinLMStateShell::CMinLMStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMinLMStateShell::CMinLMStateShell(CMinLMState &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinLMStateShell::~CMinLMStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needf                          |
//+------------------------------------------------------------------+
bool CMinLMStateShell::GetNeedF(void)
  {
//--- return result
   return(m_innerobj.m_needf);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needf                         |
//+------------------------------------------------------------------+
void CMinLMStateShell::SetNeedF(const bool b)
  {
//--- change value
   m_innerobj.m_needf=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needfg                         |
//+------------------------------------------------------------------+
bool CMinLMStateShell::GetNeedFG(void)
  {
//--- return result
   return(m_innerobj.m_needfg);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needfg                        |
//+------------------------------------------------------------------+
void CMinLMStateShell::SetNeedFG(const bool b)
  {
//--- change value
   m_innerobj.m_needfg=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needfgh                        |
//+------------------------------------------------------------------+
bool CMinLMStateShell::GetNeedFGH(void)
  {
//--- return result
   return(m_innerobj.m_needfgh);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needfgh                       |
//+------------------------------------------------------------------+
void CMinLMStateShell::SetNeedFGH(const bool b)
  {
//--- change value
   m_innerobj.m_needfgh=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needfi                         |
//+------------------------------------------------------------------+
bool CMinLMStateShell::GetNeedFI(void)
  {
//--- return result
   return(m_innerobj.m_needfi);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needfi                        |
//+------------------------------------------------------------------+
void CMinLMStateShell::SetNeedFI(const bool b)
  {
//--- change value
   m_innerobj.m_needfi=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needfij                        |
//+------------------------------------------------------------------+
bool CMinLMStateShell::GetNeedFIJ(void)
  {
//--- return result
   return(m_innerobj.m_needfij);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needfij                       |
//+------------------------------------------------------------------+
void CMinLMStateShell::SetNeedFIJ(const bool b)
  {
//--- change value
   m_innerobj.m_needfij=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable xupdated                       |
//+------------------------------------------------------------------+
bool CMinLMStateShell::GetXUpdated(void)
  {
//--- return result
   return(m_innerobj.m_xupdated);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable xupdated                      |
//+------------------------------------------------------------------+
void CMinLMStateShell::SetXUpdated(const bool b)
  {
//--- change value
   m_innerobj.m_xupdated=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable f                              |
//+------------------------------------------------------------------+
double CMinLMStateShell::GetF(void)
  {
//--- return result
   return(m_innerobj.m_f);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable f                             |
//+------------------------------------------------------------------+
void CMinLMStateShell::SetF(const double d)
  {
//--- change value
   m_innerobj.m_f=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMinLMState *CMinLMStateShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Optimization report, filled by MinLMResults() function           |
//| FIELDS:                                                          |
//| * TerminationType, completetion code:                            |
//|     * -9    derivative correctness check failed;                 |
//|             see Rep.WrongNum, Rep.WrongI, Rep.WrongJ for         |
//|             more information.                                    |
//|     *  1    relative function improvement is no more than        |
//|             EpsF.                                                |
//|     *  2    relative step is no more than EpsX.                  |
//|     *  4    gradient is no more than EpsG.                       |
//|     *  5    MaxIts steps was taken                               |
//|     *  7    stopping conditions are too stringent,               |
//|             further improvement is impossible                    |
//| * IterationsCount, contains iterations count                     |
//| * NFunc, number of function calculations                         |
//| * NJac, number of Jacobi matrix calculations                     |
//| * NGrad, number of gradient calculations                         |
//| * NHess, number of Hessian calculations                          |
//| * NCholesky, number of Cholesky decomposition calculations       |
//+------------------------------------------------------------------+
class CMinLMReport
  {
public:
   //--- variables
   int               m_iterationscount;
   int               m_terminationtype;
   int               m_nfunc;
   int               m_njac;
   int               m_ngrad;
   int               m_nhess;
   int               m_ncholesky;
   //--- constructor, destructor
                     CMinLMReport(void);
                    ~CMinLMReport(void);
   //--- copy
   void              Copy(CMinLMReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinLMReport::CMinLMReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinLMReport::~CMinLMReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMinLMReport::Copy(CMinLMReport &obj)
  {
//--- copy variables
   m_iterationscount=obj.m_iterationscount;
   m_terminationtype=obj.m_terminationtype;
   m_nfunc=obj.m_nfunc;
   m_njac=obj.m_njac;
   m_ngrad=obj.m_ngrad;
   m_nhess=obj.m_nhess;
   m_ncholesky=obj.m_ncholesky;
  }
//+------------------------------------------------------------------+
//| Optimization report, filled by MinLMResults() function           |
//| FIELDS:                                                          |
//| * TerminationType, completetion code:                            |
//|     * -9    derivative correctness check failed;                 |
//|             see Rep.WrongNum, Rep.WrongI, Rep.WrongJ for         |
//|             more information.                                    |
//|     *  1    relative function improvement is no more than        |
//|             EpsF.                                                |
//|     *  2    relative step is no more than EpsX.                  |
//|     *  4    gradient is no more than EpsG.                       |
//|     *  5    MaxIts steps was taken                               |
//|     *  7    stopping conditions are too stringent,               |
//|             further improvement is impossible                    |
//| * IterationsCount, contains iterations count                     |
//| * NFunc, number of function calculations                         |
//| * NJac, number of Jacobi matrix calculations                     |
//| * NGrad, number of gradient calculations                         |
//| * NHess, number of Hessian calculations                          |
//| * NCholesky, number of Cholesky decomposition calculations       |
//+------------------------------------------------------------------+
class CMinLMReportShell
  {
private:
   CMinLMReport      m_innerobj;
public:
   //--- constructors, destructor
                     CMinLMReportShell(void);
                     CMinLMReportShell(CMinLMReport &obj);
                    ~CMinLMReportShell(void);
   //--- methods
   int               GetIterationsCount(void);
   void              SetIterationsCount(const int i);
   int               GetTerminationType(void);
   void              SetTerminationType(const int i);
   int               GetNFunc(void);
   void              SetNFunc(const int i);
   int               GetNJAC(void);
   void              SetNJAC(const int i);
   int               GetNGrad(void);
   void              SetNGrad(const int i);
   int               GetNHess(void);
   void              SetNHess(const int i);
   int               GetNCholesky(void);
   void              SetNCholesky(const int i);
   CMinLMReport     *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinLMReportShell::CMinLMReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMinLMReportShell::CMinLMReportShell(CMinLMReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinLMReportShell::~CMinLMReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable iterationscount                |
//+------------------------------------------------------------------+
int CMinLMReportShell::GetIterationsCount(void)
  {
//--- return result
   return(m_innerobj.m_iterationscount);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable iterationscount               |
//+------------------------------------------------------------------+
void CMinLMReportShell::SetIterationsCount(const int i)
  {
//--- change value
   m_innerobj.m_iterationscount=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable terminationtype                |
//+------------------------------------------------------------------+
int CMinLMReportShell::GetTerminationType(void)
  {
//--- return result
   return(m_innerobj.m_terminationtype);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable terminationtype               |
//+------------------------------------------------------------------+
void CMinLMReportShell::SetTerminationType(const int i)
  {
//--- change value
   m_innerobj.m_terminationtype=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable nfunc                          |
//+------------------------------------------------------------------+
int CMinLMReportShell::GetNFunc(void)
  {
//--- return result
   return(m_innerobj.m_nfunc);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable nfunc                         |
//+------------------------------------------------------------------+
void CMinLMReportShell::SetNFunc(const int i)
  {
//--- change value
   m_innerobj.m_nfunc=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable njac                           |
//+------------------------------------------------------------------+
int CMinLMReportShell::GetNJAC(void)
  {
//--- return result
   return(m_innerobj.m_njac);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable njac                          |
//+------------------------------------------------------------------+
void CMinLMReportShell::SetNJAC(const int i)
  {
//--- change value
   m_innerobj.m_njac=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable ngrad                          |
//+------------------------------------------------------------------+
int CMinLMReportShell::GetNGrad(void)
  {
//--- return result
   return(m_innerobj.m_ngrad);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable ngrad                         |
//+------------------------------------------------------------------+
void CMinLMReportShell::SetNGrad(const int i)
  {
//--- change value
   m_innerobj.m_ngrad=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable nhess                          |
//+------------------------------------------------------------------+
int CMinLMReportShell::GetNHess(void)
  {
//--- return result
   return(m_innerobj.m_nhess);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable nhess                         |
//+------------------------------------------------------------------+
void CMinLMReportShell::SetNHess(const int i)
  {
//--- change value
   m_innerobj.m_nhess=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable ncholesky                      |
//+------------------------------------------------------------------+
int CMinLMReportShell::GetNCholesky(void)
  {
//--- return result
   return(m_innerobj.m_ncholesky);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable ncholesky                     |
//+------------------------------------------------------------------+
void CMinLMReportShell::SetNCholesky(const int i)
  {
//--- change value
   m_innerobj.m_ncholesky=i;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMinLMReport *CMinLMReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Levenberg-Marquardt method                                       |
//+------------------------------------------------------------------+
class CMinLM
  {
private:
   //--- private methods
   static void       LMPRepare(const int n,const int m,bool havegrad,CMinLMState &state);
   static void       ClearRequestFields(CMinLMState &state);
   static bool       IncreaseLambda(double &lambdav,double &nu);
   static void       DecreaseLambda(double &lambdav,double &nu);
   static double     BoundedScaledAntigradNorm(CMinLMState &state,double &x[],double &g[]);
   //--- auxiliary functions for MinLMIteration
   static void       Func_lbl_rcomm(CMinLMState &state,int n,int m,int iflag,int i,int k,bool bflag,double v,double s,double t);
   static bool       Func_lbl_16(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_19(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_20(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_21(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_22(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_24(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_25(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_28(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_31(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_39(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_40(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_41(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_48(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_49(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
   static bool       Func_lbl_55(CMinLMState &state,int &n,int &m,int &iflag,int &i,int &k,bool &bflag,double &v,double &s,double &t);
public:
   //--- class constants
   static const int  m_lmmodefj;
   static const int  m_lmmodefgj;
   static const int  m_lmmodefgh;
   static const int  m_lmflagnopreLBFGS;
   static const int  m_lmflagnointLBFGS;
   static const int  m_lmpreLBFGSm;
   static const int  m_lmintLBFGSits;
   static const int  m_lbfgsnorealloc;
   static const double m_lambdaup;
   static const double m_lambdadown;
   static const double m_suspiciousnu;
   static const int  m_smallmodelage;
   static const int  m_additers;
   //--- constructor, destructor
                     CMinLM(void);
                    ~CMinLM(void);
   //--- public methods
   static void       MinLMCreateVJ(const int n,const int m,double &x[],CMinLMState &state);
   static void       MinLMCreateV(const int n,const int m,double &x[],const double diffstep,CMinLMState &state);
   static void       MinLMCreateFGH(const int n,double &x[],CMinLMState &state);
   static void       MinLMSetCond(CMinLMState &state,const double epsg,const double epsf,double epsx,const int maxits);
   static void       MinLMSetXRep(CMinLMState &state,const bool needxrep);
   static void       MinLMSetStpMax(CMinLMState &state,const double stpmax);
   static void       MinLMSetScale(CMinLMState &state,double &s[]);
   static void       MinLMSetBC(CMinLMState &state,double &bndl[],double &bndu[]);
   static void       MinLMSetAccType(CMinLMState &state,int acctype);
   static void       MinLMResults(CMinLMState &state,double &x[],CMinLMReport &rep);
   static void       MinLMResultsBuf(CMinLMState &state,double &x[],CMinLMReport &rep);
   static void       MinLMRestartFrom(CMinLMState &state,double &x[]);
   static void       MinLMCreateVGJ(const int n,const int m,double &x[],CMinLMState &state);
   static void       MinLMCreateFGJ(const int n,const int m,double &x[],CMinLMState &state);
   static void       MinLMCreateFJ(const int n,const int m,double &x[],CMinLMState &state);
   static bool       MinLMIteration(CMinLMState &state);
  };
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const int    CMinLM::m_lmmodefj=0;
const int    CMinLM::m_lmmodefgj=1;
const int    CMinLM::m_lmmodefgh=2;
const int    CMinLM::m_lmflagnopreLBFGS=1;
const int    CMinLM::m_lmflagnointLBFGS=2;
const int    CMinLM::m_lmpreLBFGSm=5;
const int    CMinLM::m_lmintLBFGSits=5;
const int    CMinLM::m_lbfgsnorealloc=1;
const double CMinLM::m_lambdaup=2.0;
const double CMinLM::m_lambdadown=0.33;
const double CMinLM::m_suspiciousnu=16;
const int    CMinLM::m_smallmodelage=3;
const int    CMinLM::m_additers=5;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinLM::CMinLM(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinLM::~CMinLM(void)
  {

  }
//+------------------------------------------------------------------+
//|                 IMPROVED LEVENBERG-MARQUARDT METHOD FOR          |
//|                  NON-LINEAR LEAST SQUARES OPTIMIZATION           |
//| DESCRIPTION:                                                     |
//| This function is used to find minimum of function which is       |
//| represented as sum of squares:                                   |
//|     F(x) = f[0]^2(x[0],...,x[n-1]) + ... +                       |
//|     + f[m-1]^2(x[0],...,x[n-1])                                  |
//| using value of function vector f[] and Jacobian of f[].          |
//| REQUIREMENTS:                                                    |
//| This algorithm will request following information during its     |
//| operation:                                                       |
//| * function vector f[] at given point X                           |
//| * function vector f[] and Jacobian of f[] (simultaneously) at    |
//| given point                                                      |
//| There are several overloaded versions of MinLMOptimize()         |
//| function which correspond to different LM-like optimization      |
//| algorithms provided by this unit. You should choose version which|
//| accepts fvec() and jac() callbacks. First one is used to         |
//| calculate f[] at given point, second one calculates f[] and      |
//| Jacobian df[i]/dx[j].                                            |
//| You can try to initialize MinLMState structure with VJ function  |
//| and then use incorrect version  of  MinLMOptimize() (for example,|
//| version which works with general form function and does not      |
//| provide Jacobian), but it will lead to exception being thrown    |
//| after first attempt to calculate Jacobian.                       |
//| USAGE:                                                           |
//| 1. User initializes algorithm state with MinLMCreateVJ() call    |
//| 2. User tunes solver parameters with MinLMSetCond(),             |
//|    MinLMSetStpMax() and other functions                          |
//| 3. User calls MinLMOptimize() function which takes algorithm     |
//|    state and callback functions.                                 |
//| 4. User calls MinLMResults() to get solution                     |
//| 5. Optionally, user may call MinLMRestartFrom() to solve another |
//|    problem with same N/M but another starting point and/or       |
//|    another function. MinLMRestartFrom() allows to reuse already  |
//|    initialized structure.                                        |
//| INPUT PARAMETERS:                                                |
//|     N       -   dimension, N>1                                   |
//|                 * if given, only leading N elements of X are     |
//|                   used                                           |
//|                 * if not given, automatically determined from    |
//|                   size of X                                      |
//|     M       -   number of functions f[i]                         |
//|     X       -   initial solution, array[0..N-1]                  |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//| NOTES:                                                           |
//| 1. you may tune stopping conditions with MinLMSetCond() function |
//| 2. if target function contains exp() or other fast growing       |
//|    functions, and optimization algorithm makes too large steps   |
//|    which leads to overflow, use MinLMSetStpMax() function to     |
//|    bound algorithm's steps.                                      |
//+------------------------------------------------------------------+
static void CMinLM::MinLMCreateVJ(const int n,const int m,double &x[],
                                  CMinLMState &state)
  {
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- initialize,check parameters
   state.m_n=n;
   state.m_m=m;
   state.m_algomode=1;
   state.m_hasf=false;
   state.m_hasfi=true;
   state.m_hasg=false;
//--- second stage of initialization
   LMPRepare(n,m,false,state);
   MinLMSetAccType(state,0);
   MinLMSetCond(state,0,0,0,0);
   MinLMSetXRep(state,false);
   MinLMSetStpMax(state,0);
   MinLMRestartFrom(state,x);
  }
//+------------------------------------------------------------------+
//|                 IMPROVED LEVENBERG-MARQUARDT METHOD FOR          |
//|                  NON-LINEAR LEAST SQUARES OPTIMIZATION           |
//| DESCRIPTION:                                                     |
//| This function is used to find minimum of function which is       |
//| represented as sum of squares:                                   |
//|     F(x) = f[0]^2(x[0],...,x[n-1]) + ... +                       |
//|     + f[m-1]^2(x[0],...,x[n-1])                                  |
//| using value of function vector f[] only. Finite differences are  |
//| used to calculate Jacobian.                                      |
//| REQUIREMENTS:                                                    |
//| This algorithm will request following information during its     |
//| operation:                                                       |
//| * function vector f[] at given point X                           |
//| There are several overloaded versions of MinLMOptimize() function|
//| which correspond to different LM-like optimization algorithms    |
//| provided by this unit. You should choose version which accepts   |
//| fvec() callback.                                                 |
//| You can try to initialize MinLMState structure with VJ function  |
//| and then use incorrect version of MinLMOptimize() (for example,  |
//| version which works with general form function and does not      |
//| accept function vector), but it will lead to exception being     |
//| thrown after first attempt to calculate Jacobian.                |
//| USAGE:                                                           |
//| 1. User initializes algorithm state with MinLMCreateV() call     |
//| 2. User tunes solver parameters with MinLMSetCond(),             |
//|    MinLMSetStpMax() and other functions                          |
//| 3. User calls MinLMOptimize() function which takes algorithm     |
//|    state and callback functions.                                 |
//| 4. User calls MinLMResults() to get solution                     |
//| 5. Optionally, user may call MinLMRestartFrom() to solve another |
//|    problem with same N/M but another starting point and/or       |
//|    another function. MinLMRestartFrom() allows to reuse already  |
//|    initialized structure.                                        |
//| INPUT PARAMETERS:                                                |
//|     N       -   dimension, N>1                                   |
//|                 * if given, only leading N elements of X are     |
//|                   used                                           |
//|                 * if not given, automatically determined from    |
//|                   size of X                                      |
//|     M       -   number of functions f[i]                         |
//|     X       -   initial solution, array[0..N-1]                  |
//|     DiffStep-   differentiation step, >0                         |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//| See also MinLMIteration, MinLMResults.                           |
//| NOTES:                                                           |
//| 1. you may tune stopping conditions with MinLMSetCond() function |
//| 2. if target function contains exp() or other fast growing       |
//|    functions, and optimization algorithm makes too large steps   |
//|    which leads to overflow, use MinLMSetStpMax() function to     |
//|    bound algorithm's steps.                                      |
//+------------------------------------------------------------------+
static void CMinLM::MinLMCreateV(const int n,const int m,double &x[],
                                 const double diffstep,CMinLMState &state)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(diffstep),__FUNCTION__+": DiffStep is not finite!"))
      return;
//--- check
   if(!CAp::Assert(diffstep>0.0,__FUNCTION__+": DiffStep<=0!"))
      return;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- initialize
   state.m_n=n;
   state.m_m=m;
   state.m_algomode=0;
   state.m_hasf=false;
   state.m_hasfi=true;
   state.m_hasg=false;
   state.m_diffstep=diffstep;
//--- second stage of initialization
   LMPRepare(n,m,false,state);
   MinLMSetAccType(state,1);
   MinLMSetCond(state,0,0,0,0);
   MinLMSetXRep(state,false);
   MinLMSetStpMax(state,0);
   MinLMRestartFrom(state,x);
  }
//+------------------------------------------------------------------+
//|     LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION  |
//| DESCRIPTION:                                                     |
//| This function is used to find minimum of general form (not       |
//| "sum-of-squares") function                                       |
//|     F = F(x[0], ..., x[n-1])                                     |
//| using its gradient and Hessian. Levenberg-Marquardt modification |
//| with L-BFGS pre-optimization and internal pre-conditioned L-BFGS |
//| optimization after each Levenberg-Marquardt step is used.        |
//| REQUIREMENTS:                                                    |
//| This algorithm will request following information during its     |
//| operation:                                                       |
//| * function value F at given point X                              |
//| * F and gradient G (simultaneously) at given point X             |
//| * F, G and Hessian H (simultaneously) at given point X           |
//| There are several overloaded versions of  MinLMOptimize()        |
//| function which correspond to different LM-like optimization      |
//| algorithms provided by this unit. You should choose version which|
//| accepts func(), grad() and hess() function pointers. First       |
//| pointer is used to calculate F at given point, second one        |
//| calculates F(x) and grad F(x), third one calculates F(x), grad   |
//| F(x), hess F(x).                                                 |
//| You can try to initialize MinLMState structure with FGH-function |
//| and then use incorrect version of MinLMOptimize() (for example,  |
//| version which does not provide Hessian matrix), but it will lead |
//| to exception being thrown after first attempt to calculate       |
//| Hessian.                                                         |
//| USAGE:                                                           |
//| 1. User initializes algorithm state with MinLMCreateFGH() call   |
//| 2. User tunes solver parameters with MinLMSetCond(),             |
//|    MinLMSetStpMax() and other functions                          |
//| 3. User calls MinLMOptimize() function which takes algorithm     |
//|    state and pointers (delegates, etc.) to callback functions.   |
//| 4. User calls MinLMResults() to get solution                     |
//| 5. Optionally, user may call MinLMRestartFrom() to solve another |
//|    problem with same N but another starting point and/or another |
//|    function. MinLMRestartFrom() allows to reuse already          |
//|    initialized structure.                                        |
//| INPUT PARAMETERS:                                                |
//|     N       -   dimension, N>1                                   |
//|                 * if given, only leading N elements of X are     |
//|                   used                                           |
//|                 * if not given, automatically determined from    |
//|                   size of X                                      |
//|     X       -   initial solution, array[0..N-1]                  |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//| NOTES:                                                           |
//| 1. you may tune stopping conditions with MinLMSetCond() function |
//| 2. if target function contains exp() or other fast growing       |
//|    functions, and optimization algorithm makes too large steps   |
//|    which leads to overflow, use MinLMSetStpMax() function to     |
//|    bound algorithm's steps.                                      |
//+------------------------------------------------------------------+
static void CMinLM::MinLMCreateFGH(const int n,double &x[],CMinLMState &state)
  {
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- initialize
   state.m_n=n;
   state.m_m=0;
   state.m_algomode=2;
   state.m_hasf=true;
   state.m_hasfi=false;
   state.m_hasg=true;
//--- init2
   LMPRepare(n,0,true,state);
   MinLMSetAccType(state,2);
   MinLMSetCond(state,0,0,0,0);
   MinLMSetXRep(state,false);
   MinLMSetStpMax(state,0);
   MinLMRestartFrom(state,x);
  }
//+------------------------------------------------------------------+
//| This function sets stopping conditions for Levenberg-Marquardt   |
//| optimization algorithm.                                          |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     EpsG    -   >=0                                              |
//|                 The  subroutine  finishes its work if the        |
//|                 condition |v|<EpsG is satisfied, where:          |
//|                 * |.| means Euclidian norm                       |
//|                 * v - scaled gradient vector, v[i]=g[i]*s[i]     |
//|                 * g - gradient                                   |
//|                 * s - scaling coefficients set by MinLMSetScale()|
//|     EpsF    -   >=0                                              |
//|                 The  subroutine  finishes  its work if on k+1-th |
//|                 iteration the condition |F(k+1)-F(k)| <=         |
//|                 <= EpsF*max{|F(k)|,|F(k+1)|,1} is satisfied.     |
//|     EpsX    -   >=0                                              |
//|                 The subroutine finishes its work if on k+1-th    |
//|                 iteration the condition |v|<=EpsX is fulfilled,  |
//|                 where:                                           |
//|                 * |.| means Euclidian norm                       |
//|                 * v - scaled step vector, v[i]=dx[i]/s[i]        |
//|                 * dx - ste pvector, dx=X(k+1)-X(k)               |
//|                 * s - scaling coefficients set by MinLMSetScale()|
//|     MaxIts  -   maximum number of iterations. If MaxIts=0, the   |
//|                 number of iterations is unlimited. Only          |
//|                 Levenberg-Marquardt iterations are counted       |
//|                 (L-BFGS/CG iterations are NOT counted because    |
//|                 their cost is very low compared to that of LM).  |
//| Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will|
//| lead to automatic stopping criterion selection (small EpsX).     |
//+------------------------------------------------------------------+
static void CMinLM::MinLMSetCond(CMinLMState &state,const double epsg,
                                 const double epsf,double epsx,
                                 const int maxits)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsg),__FUNCTION__+": EpsG is not finite number!"))
      return;
//--- check
   if(!CAp::Assert(epsg>=0.0,__FUNCTION__+": negative EpsG!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsf),__FUNCTION__+": EpsF is not finite number!"))
      return;
//--- check
   if(!CAp::Assert(epsf>=0.0,__FUNCTION__+": negative EpsF!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsx),__FUNCTION__+": EpsX is not finite number!"))
      return;
//--- check
   if(!CAp::Assert(epsx>=0.0,__FUNCTION__+": negative EpsX!"))
      return;
//--- check
   if(!CAp::Assert(maxits>=0,__FUNCTION__+": negative MaxIts!"))
      return;
//--- check
   if(((epsg==0.0 && epsf==0.0) && epsx==0.0) && maxits==0)
      epsx=1.0E-6;
//--- change values
   state.m_epsg=epsg;
   state.m_epsf=epsf;
   state.m_epsx=epsx;
   state.m_maxits=maxits;
  }
//+------------------------------------------------------------------+
//| This function turns on/off reporting.                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     NeedXRep-   whether iteration reports are needed or not      |
//| If NeedXRep is True, algorithm will call rep() callback function |
//| if it is provided to MinLMOptimize(). Both Levenberg-Marquardt   |
//| and internal L-BFGS iterations are reported.                     |
//+------------------------------------------------------------------+
static void CMinLM::MinLMSetXRep(CMinLMState &state,const bool needxrep)
  {
//--- change value
   state.m_xrep=needxrep;
  }
//+------------------------------------------------------------------+
//| This function sets maximum step length                           |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     StpMax  -   maximum step length, >=0. Set StpMax to 0.0, if  |
//|                 you don't want to limit step length.             |
//| Use this subroutine when you optimize target function which      |
//| contains exp() or other fast growing functions, and optimization |
//| algorithm makes too large steps which leads to overflow. This    |
//| function allows us to reject steps that are too large (and       |
//| therefore expose us to the possible overflow) without actually   |
//| calculating function value at the x+stp*d.                       |
//| NOTE: non-zero StpMax leads to moderate performance degradation  |
//| because intermediate step of preconditioned L-BFGS optimization  |
//| is incompatible with limits on step size.                        |
//+------------------------------------------------------------------+
static void CMinLM::MinLMSetStpMax(CMinLMState &state,const double stpmax)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(stpmax),__FUNCTION__+": StpMax is not finite!"))
      return;
//--- check
   if(!CAp::Assert(stpmax>=0.0,__FUNCTION__+": StpMax<0!"))
      return;
//--- change value
   state.m_stpmax=stpmax;
  }
//+------------------------------------------------------------------+
//| This function sets scaling coefficients for LM optimizer.        |
//| ALGLIB optimizers use scaling matrices to test stopping          |
//| conditions (step size and gradient are scaled before comparison  |
//| with tolerances). Scale of the I-th variable is a translation    |
//| invariant measure of:                                            |
//| a) "how large" the variable is                                   |
//| b) how large the step should be to make significant changes in   |
//| the function                                                     |
//| Generally, scale is NOT considered to be a form of               |
//| preconditioner. But LM optimizer is unique in that it uses       |
//| scaling matrix both in the stopping condition tests and as       |
//| Marquardt damping factor.                                        |
//| Proper scaling is very important for the algorithm performance.  |
//| It is less important for the quality of results, but still has   |
//| some influence (it is easier to converge when variables are      |
//| properly scaled, so premature stopping is possible when very     |
//| badly scalled variables are combined with relaxed stopping       |
//| conditions).                                                     |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure stores algorithm state                 |
//|     S       -   array[N], non-zero scaling coefficients          |
//|                 S[i] may be negative, sign doesn't matter.       |
//+------------------------------------------------------------------+
static void CMinLM::MinLMSetScale(CMinLMState &state,double &s[])
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(CAp::Len(s)>=state.m_n,__FUNCTION__+": Length(S)<N"))
      return;
   for(i=0;i<=state.m_n-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(s[i]),__FUNCTION__+": S contains infinite or NAN elements"))
         return;
      //--- check
      if(!CAp::Assert(s[i]!=0.0,__FUNCTION__+": S contains zero elements"))
         return;
      //--- change value
      state.m_s[i]=MathAbs(s[i]);
     }
  }
//+------------------------------------------------------------------+
//| This function sets boundary constraints for LM optimizer         |
//| Boundary constraints are inactive by default (after initial      |
//| creation). They are preserved until explicitly turned off with   |
//| another SetBC() call.                                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure stores algorithm state                 |
//|     BndL    -   lower bounds, array[N].                          |
//|                 If some (all) variables are unbounded, you may   |
//|                 specify very small number or -INF (latter is     |
//|                 recommended because it will allow solver to use  |
//|                 better algorithm).                               |
//|     BndU    -   upper bounds, array[N].                          |
//|                 If some (all) variables are unbounded, you may   |
//|                 specify very large number or +INF (latter is     |
//|                 recommended because it will allow solver to use  |
//|                 better algorithm).                               |
//| NOTE 1: it is possible to specify BndL[i]=BndU[i]. In this case  |
//| I-th variable will be "frozen" at X[i]=BndL[i]=BndU[i].          |
//| NOTE 2: this solver has following useful properties:             |
//| * bound constraints are always satisfied exactly                 |
//| * function is evaluated only INSIDE area specified by bound      |
//|   constraints or at its boundary                                 |
//+------------------------------------------------------------------+
static void CMinLM::MinLMSetBC(CMinLMState &state,double &bndl[],double &bndu[])
  {
//--- create variables
   int i=0;
   int n=0;
//--- initialization
   n=state.m_n;
//--- check
   if(!CAp::Assert(CAp::Len(bndl)>=n,__FUNCTION__+": Length(BndL)<N"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(bndu)>=n,__FUNCTION__+": Length(BndU)<N"))
      return;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(bndl[i]) || CInfOrNaN::IsNegativeInfinity(bndl[i]),"MinLMSetBC: BndL contains NAN or +INF"))
         return;
      //--- check
      if(!CAp::Assert(CMath::IsFinite(bndu[i]) || CInfOrNaN::IsPositiveInfinity(bndu[i]),"MinLMSetBC: BndU contains NAN or -INF"))
         return;
      //--- change values
      state.m_bndl[i]=bndl[i];
      state.m_havebndl[i]=CMath::IsFinite(bndl[i]);
      state.m_bndu[i]=bndu[i];
      state.m_havebndu[i]=CMath::IsFinite(bndu[i]);
     }
  }
//+------------------------------------------------------------------+
//| This function is used to change acceleration settings            |
//| You can choose between three acceleration strategies:            |
//| * AccType=0, no acceleration.                                    |
//| * AccType=1, secant updates are used to update quadratic model   |
//|   after each iteration. After fixed number of iterations (or     |
//|   after model breakdown) we recalculate quadratic model using    |
//|   analytic Jacobian or finite differences. Number of secant-based|
//|   iterations depends on optimization settings: about 3           |
//|   iterations - when we have analytic Jacobian, up to 2*N         |
//|   iterations - when we use finite differences to calculate       |
//|   Jacobian.                                                      |
//| AccType=1 is recommended when Jacobian calculation cost is       |
//| prohibitive high (several Mx1 function vector calculations       |
//| followed by several NxN Cholesky factorizations are faster than  |
//| calculation of one M*N  Jacobian). It should also be used when we|
//| have no Jacobian, because finite difference approximation takes  |
//| too much time to compute.                                        |
//| Table below list optimization protocols (XYZ protocol corresponds|
//| to MinLMCreateXYZ) and acceleration types they support (and use  |
//| by default).                                                     |
//| ACCELERATION TYPES SUPPORTED BY OPTIMIZATION PROTOCOLS:          |
//| protocol    0   1   comment                                      |
//| V           +   +                                                |
//| VJ          +   +                                                |
//| FGH         +                                                    |
//| DAFAULT VALUES:                                                  |
//| protocol    0   1   comment                                      |
//| V               x   without acceleration it is so slooooooooow   |
//| VJ          x                                                    |
//| FGH         x                                                    |
//| NOTE: this  function should be called before optimization.       |
//| Attempt to call it during algorithm iterations may result in     |
//| unexpected behavior.                                             |
//| NOTE: attempt to call this function with unsupported             |
//| protocol/acceleration combination will result in exception being |
//| thrown.                                                          |
//+------------------------------------------------------------------+
static void CMinLM::MinLMSetAccType(CMinLMState &state,int acctype)
  {
//--- check
   if(!CAp::Assert((acctype==0||acctype==1)||acctype==2,__FUNCTION__+": incorrect AccType!"))
      return;
//--- check
   if(acctype==2)
      acctype=0;
//--- check
   if(acctype==0)
     {
      state.m_maxmodelage=0;
      state.m_makeadditers=false;
      //--- exit the function
      return;
     }
//--- check
   if(acctype==1)
     {
      //--- check
      if(!CAp::Assert(state.m_hasfi,__FUNCTION__+": AccType=1 is incompatible with current protocol!"))
         return;
      //--- check
      if(state.m_algomode==0)
         state.m_maxmodelage=2*state.m_n;
      else
         state.m_maxmodelage=m_smallmodelage;
      state.m_makeadditers=false;
      //--- exit the function
      return;
     }
  }
//+------------------------------------------------------------------+
//| Levenberg-Marquardt algorithm results                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   algorithm state                                  |
//| OUTPUT PARAMETERS:                                               |
//|     X       -   array[0..N-1], solution                          |
//|     Rep     -   optimization report;                             |
//|                 see comments for this structure for more info.   |
//+------------------------------------------------------------------+
static void CMinLM::MinLMResults(CMinLMState &state,double &x[],
                                 CMinLMReport &rep)
  {
//--- reset memory
   ArrayResizeAL(x,0);
//--- function call
   MinLMResultsBuf(state,x,rep);
  }
//+------------------------------------------------------------------+
//| Levenberg-Marquardt algorithm results                            |
//| Buffered implementation of MinLMResults(), which uses            |
//| pre-allocated buffer to store X[]. If buffer size is too small,  |
//| it resizes buffer. It is intended to be used in the inner cycles |
//| of performance critical algorithms where array reallocation      |
//| penalty is too large to be ignored.                              |
//+------------------------------------------------------------------+
static void CMinLM::MinLMResultsBuf(CMinLMState &state,double &x[],
                                    CMinLMReport &rep)
  {
//--- create a variable
   int i_=0;
//--- check
   if(CAp::Len(x)<state.m_n)
      ArrayResizeAL(x,state.m_n);
//--- copy
   for(i_=0;i_<=state.m_n-1;i_++)
      x[i_]=state.m_x[i_];
//--- change values
   rep.m_iterationscount=state.m_repiterationscount;
   rep.m_terminationtype=state.m_repterminationtype;
   rep.m_nfunc=state.m_repnfunc;
   rep.m_njac=state.m_repnjac;
   rep.m_ngrad=state.m_repngrad;
   rep.m_nhess=state.m_repnhess;
   rep.m_ncholesky=state.m_repncholesky;
  }
//+------------------------------------------------------------------+
//| This subroutine restarts LM algorithm from new point. All        |
//| optimization parameters are left unchanged.                      |
//| This function allows to solve multiple optimization problems     |
//| (which must have same number of dimensions) without object       |
//| reallocation penalty.                                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure used for reverse communication         |
//|                 previously allocated with MinLMCreateXXX call.   |
//|     X       -   new starting point.                              |
//+------------------------------------------------------------------+
static void CMinLM::MinLMRestartFrom(CMinLMState &state,double &x[])
  {
//--- create a variable
   int i_=0;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=state.m_n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,state.m_n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- copy
   for(i_=0;i_<=state.m_n-1;i_++)
      state.m_xbase[i_]=x[i_];
//--- allocation
   ArrayResizeAL(state.m_rstate.ia,5);
   ArrayResizeAL(state.m_rstate.ba,1);
   ArrayResizeAL(state.m_rstate.ra,3);
   state.m_rstate.stage=-1;
//--- function call
   ClearRequestFields(state);
  }
//+------------------------------------------------------------------+
//| This is obsolete function.                                       |
//| Since ALGLIB 3.3 it is equivalent to MinLMCreateVJ().            |
//+------------------------------------------------------------------+
static void CMinLM::MinLMCreateVGJ(const int n,const int m,double &x[],
                                   CMinLMState &state)
  {
//--- function call
   MinLMCreateVJ(n,m,x,state);
  }
//+------------------------------------------------------------------+
//| This is obsolete function.                                       |
//| Since ALGLIB 3.3 it is equivalent to MinLMCreateFJ().            |
//+------------------------------------------------------------------+
static void CMinLM::MinLMCreateFGJ(const int n,const int m,double &x[],
                                   CMinLMState &state)
  {
//--- function call
   MinLMCreateFJ(n,m,x,state);
  }
//+------------------------------------------------------------------+
//| This function is considered obsolete since ALGLIB 3.1.0 and is   |
//| present for backward compatibility only. We recommend to use     |
//| MinLMCreateVJ, which provides similar, but more consistent and   |
//| feature-rich interface.                                          |
//+------------------------------------------------------------------+
static void CMinLM::MinLMCreateFJ(const int n,const int m,double &x[],
                                  CMinLMState &state)
  {
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- initialize
   state.m_n=n;
   state.m_m=m;
   state.m_algomode=1;
   state.m_hasf=true;
   state.m_hasfi=false;
   state.m_hasg=false;
//--- init 2
   LMPRepare(n,m,true,state);
//--- function call
   MinLMSetAccType(state,0);
//--- function call
   MinLMSetCond(state,0,0,0,0);
//--- function call
   MinLMSetXRep(state,false);
//--- function call
   MinLMSetStpMax(state,0);
//--- function call
   MinLMRestartFrom(state,x);
  }
//+------------------------------------------------------------------+
//| Prepare internal structures (except for RComm).                  |
//| Note: M must be zero for FGH mode, non-zero for V/VJ/FJ/FGJ mode.|
//+------------------------------------------------------------------+
static void CMinLM::LMPRepare(const int n,const int m,bool havegrad,
                              CMinLMState &state)
  {
//--- create a variable
   int i=0;
//--- check
   if(n<=0 || m<0)
      return;
//--- check
   if(havegrad)
      ArrayResizeAL(state.m_g,n);
//--- check
   if(m!=0)
     {
      //--- allocation
      state.m_j.Resize(m,n);
      ArrayResizeAL(state.m_fi,m);
      ArrayResizeAL(state.m_fibase,m);
      ArrayResizeAL(state.m_deltaf,m);
      ArrayResizeAL(state.m_fm1,m);
      ArrayResizeAL(state.m_fp1,m);
     }
   else
      state.m_h.Resize(n,n);
//--- allocation
   ArrayResizeAL(state.m_x,n);
   ArrayResizeAL(state.m_deltax,n);
   state.m_quadraticmodel.Resize(n,n);
   ArrayResizeAL(state.m_xbase,n);
   ArrayResizeAL(state.m_gbase,n);
   ArrayResizeAL(state.m_xdir,n);
   ArrayResizeAL(state.m_tmp0,n);
//--- prepare internal L-BFGS
   for(i=0;i<=n-1;i++)
      state.m_x[i]=0;
//--- function call
   CMinLBFGS::MinLBFGSCreate(n,MathMin(m_additers,n),state.m_x,state.m_internalstate);
//--- function call
   CMinLBFGS::MinLBFGSSetCond(state.m_internalstate,0.0,0.0,0.0,MathMin(m_additers,n));
//--- Prepare internal QP solver
   CMinQP::MinQPCreate(n,state.m_qpstate);
//--- function call
   CMinQP::MinQPSetAlgoCholesky(state.m_qpstate);
//--- Prepare boundary constraints
   ArrayResizeAL(state.m_bndl,n);
   ArrayResizeAL(state.m_bndu,n);
   ArrayResizeAL(state.m_havebndl,n);
   ArrayResizeAL(state.m_havebndu,n);
   for(i=0;i<=n-1;i++)
     {
      //--- change values
      state.m_bndl[i]=CInfOrNaN::NegativeInfinity();
      state.m_havebndl[i]=false;
      state.m_bndu[i]=CInfOrNaN::PositiveInfinity();
      state.m_havebndu[i]=false;
     }
//--- Prepare scaling matrix
   ArrayResizeAL(state.m_s,n);
   for(i=0;i<=n-1;i++)
      state.m_s[i]=1.0;
  }
//+------------------------------------------------------------------+
//| Clears request fileds (to be sure that we don't forgot to clear  |
//| something)                                                       |
//+------------------------------------------------------------------+
static void CMinLM::ClearRequestFields(CMinLMState &state)
  {
//--- change values
   state.m_needf=false;
   state.m_needfg=false;
   state.m_needfgh=false;
   state.m_needfij=false;
   state.m_needfi=false;
   state.m_xupdated=false;
  }
//+------------------------------------------------------------------+
//| Increases lambda, returns False when there is a danger of        |
//| overflow                                                         |
//+------------------------------------------------------------------+
static bool CMinLM::IncreaseLambda(double &lambdav,double &nu)
  {
//--- create variables
   bool   result;
   double lnlambda=0;
   double lnnu=0;
   double lnlambdaup=0;
   double lnmax=0;
//--- initialization
   result=false;
   lnlambda=MathLog(lambdav);
   lnlambdaup=MathLog(m_lambdaup);
   lnnu=MathLog(nu);
   lnmax=MathLog(CMath::m_maxrealnumber);
//--- check
   if(lnlambda+lnlambdaup+lnnu>0.25*lnmax)
     {
      //--- return result
      return(result);
     }
//--- check
   if(lnnu+MathLog(2)>lnmax)
      return(result);
//--- change values
   lambdav=lambdav*m_lambdaup*nu;
   nu=nu*2;
   result=true;
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Decreases lambda, but leaves it unchanged when there is danger of|
//| underflow.                                                       |
//+------------------------------------------------------------------+
static void CMinLM::DecreaseLambda(double &lambdav,double &nu)
  {
//--- initialization
   nu=1;
//--- check
   if(MathLog(lambdav)+MathLog(m_lambdadown)<MathLog(CMath::m_minrealnumber))
      lambdav=CMath::m_minrealnumber;
   else
      lambdav=lambdav*m_lambdadown;
  }
//+------------------------------------------------------------------+
//| Returns norm of bounded scaled anti-gradient.                    |
//| Bounded antigradient is a vector obtained from anti-gradient by  |
//| zeroing components which point outwards:                         |
//|     result = norm(v)                                             |
//|     v[i]=0     if ((-g[i]<0)and(x[i]=bndl[i])) or                |
//|                   ((-g[i]>0)and(x[i]=bndu[i]))                   |
//|     v[i]=-g[i]*s[i] otherwise, where s[i] is a scale for I-th    |
//|     variable                                                     |
//| This function may be used to check a stopping criterion.         |
//+------------------------------------------------------------------+
static double CMinLM::BoundedScaledAntigradNorm(CMinLMState &state,
                                                double &x[],double &g[])
  {
//--- create variables
   double result=0;
   int    n=0;
   int    i=0;
   double v=0;
//--- initialization
   result=0;
   n=state.m_n;
   for(i=0;i<=n-1;i++)
     {
      v=-(g[i]*state.m_s[i]);
      //--- check
      if(state.m_havebndl[i])
        {
         //--- check
         if(x[i]<=state.m_bndl[i] && (double)(-g[i])<0.0)
            v=0;
        }
      //--- check
      if(state.m_havebndu[i])
        {
         //--- check
         if(x[i]>=state.m_bndu[i] && (double)(-g[i])>0.0)
            v=0;
        }
      result=result+CMath::Sqr(v);
     }
//--- return result
   return(MathSqrt(result));
  }
//+------------------------------------------------------------------+
//| NOTES:                                                           |
//| 1. Depending on function used to create state structure, this    |
//|    algorithm may accept Jacobian and/or Hessian and/or gradient. |
//|    According to the said above, there ase several versions of    |
//|    this function, which accept different sets of callbacks.      |
//|    This flexibility opens way to subtle errors - you may create  |
//|    state with MinLMCreateFGH() (optimization using Hessian), but |
//|    call function which does not accept Hessian. So when          |
//|    algorithm will request Hessian, there will be no callback to  |
//|    call. In this case exception will be thrown.                  |
//|    Be careful to avoid such errors because there is no way to    |
//|    find them at compile time - you can see them at runtime only. |
//+------------------------------------------------------------------+
static bool CMinLM::MinLMIteration(CMinLMState &state)
  {
//--- create variables
   int    n=0;
   int    m=0;
   bool   bflag;
   int    iflag=0;
   double v=0;
   double s=0;
   double t=0;
   int    i=0;
   int    k=0;
   int    i_=0;
//--- This code initializes locals by:
//--- * random values determined during code
//---   generation - on first subroutine call
//--- * values from previous call - on subsequent calls
   if(state.m_rstate.stage>=0)
     {
      //--- initialization
      n=state.m_rstate.ia[0];
      m=state.m_rstate.ia[1];
      iflag=state.m_rstate.ia[2];
      i=state.m_rstate.ia[3];
      k=state.m_rstate.ia[4];
      bflag=state.m_rstate.ba[0];
      v=state.m_rstate.ra[0];
      s=state.m_rstate.ra[1];
      t=state.m_rstate.ra[2];
     }
   else
     {
      //--- initialization
      n=-983;
      m=-989;
      iflag=-834;
      i=900;
      k=-287;
      bflag=false;
      v=214;
      s=-338;
      t=-686;
     }
//--- check
   if(state.m_rstate.stage==0)
     {
      //--- change value
      state.m_needf=false;
      //--- function call, return result
      return(Func_lbl_19(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- check
   if(state.m_rstate.stage==1)
     {
      //--- change values
      state.m_needfi=false;
      v=0.0;
      for(i_=0;i_<=m-1;i_++)
         v+=state.m_fi[i_]*state.m_fi[i_];
      state.m_f=v;
      //--- function call, return result
      return(Func_lbl_19(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- check
   if(state.m_rstate.stage==2)
     {
      //--- change value
      state.m_xupdated=false;
      //--- function call, return result
      return(Func_lbl_16(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- check
   if(state.m_rstate.stage==3)
     {
      //--- change value
      state.m_repnfunc=state.m_repnfunc+1;
      //--- copy
      for(i_=0;i_<=m-1;i_++)
         state.m_fm1[i_]=state.m_fi[i_];
      for(i_=0;i_<=n-1;i_++)
         state.m_x[i_]=state.m_xbase[i_];
      state.m_x[k]=state.m_x[k]+state.m_s[k]*state.m_diffstep;
      //--- check
      if(state.m_havebndl[k])
         state.m_x[k]=MathMax(state.m_x[k],state.m_bndl[k]);
      //--- check
      if(state.m_havebndu[k])
         state.m_x[k]=MathMin(state.m_x[k],state.m_bndu[k]);
      state.m_xp1=state.m_x[k];
      //--- function call
      ClearRequestFields(state);
      state.m_needfi=true;
      state.m_rstate.stage=4;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==4)
     {
      //--- change value
      state.m_repnfunc=state.m_repnfunc+1;
      //--- copy
      for(i_=0;i_<=m-1;i_++)
         state.m_fp1[i_]=state.m_fi[i_];
      v=state.m_xp1-state.m_xm1;
      //--- check
      if(v!=0.0)
        {
         v=1/v;
         for(i_=0;i_<=m-1;i_++)
            state.m_j[i_].Set(k,v*state.m_fp1[i_]);
         for(i_=0;i_<=m-1;i_++)
            state.m_j[i_].Set(k,state.m_j[i_][k]-v*state.m_fm1[i_]);
        }
      else
        {
         for(i=0;i<=m-1;i++)
            state.m_j[i].Set(k,0);
        }
      k=k+1;
      //--- function call, return result
      return(Func_lbl_28(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- check
   if(state.m_rstate.stage==5)
     {
      //--- change values
      state.m_needfi=false;
      state.m_repnfunc=state.m_repnfunc+1;
      state.m_repnjac=state.m_repnjac+1;
      //--- New model
      state.m_modelage=0;
      //--- function call, return result
      return(Func_lbl_25(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- check
   if(state.m_rstate.stage==6)
     {
      //--- change values
      state.m_needfij=false;
      state.m_repnfunc=state.m_repnfunc+1;
      state.m_repnjac=state.m_repnjac+1;
      //--- New model
      state.m_modelage=0;
      //--- function call, return result
      return(Func_lbl_25(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- check
   if(state.m_rstate.stage==7)
     {
      //--- change values
      state.m_needfgh=false;
      state.m_repnfunc=state.m_repnfunc+1;
      state.m_repngrad=state.m_repngrad+1;
      state.m_repnhess=state.m_repnhess+1;
      //--- function call
      CAblas::RMatrixCopy(n,n,state.m_h,0,0,state.m_quadraticmodel,0,0);
      for(i_=0;i_<=n-1;i_++)
         state.m_gbase[i_]=state.m_g[i_];
      state.m_fbase=state.m_f;
      //--- set control variables
      bflag=true;
      state.m_modelage=0;
      //--- function call, return result
      return(Func_lbl_31(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- check
   if(state.m_rstate.stage==8)
     {
      //--- change value
      state.m_xupdated=false;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_rstate.stage==9)
     {
      //--- change value
      state.m_xupdated=false;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_rstate.stage==10)
     {
      //--- change values
      state.m_needfi=false;
      v=0.0;
      for(i_=0;i_<=m-1;i_++)
         v+=state.m_fi[i_]*state.m_fi[i_];
      state.m_f=v;
      //--- copy
      for(i_=0;i_<=m-1;i_++)
         state.m_deltaf[i_]=state.m_fi[i_];
      for(i_=0;i_<=m-1;i_++)
         state.m_deltaf[i_]=state.m_deltaf[i_]-state.m_fibase[i_];
      state.m_deltafready=true;
      //--- function call, return result
      return(Func_lbl_48(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- check
   if(state.m_rstate.stage==11)
     {
      //--- change value
      state.m_needf=false;
      //--- function call, return result
      return(Func_lbl_48(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- check
   if(state.m_rstate.stage==12)
     {
      //--- change value
      state.m_xupdated=false;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_rstate.stage==13)
     {
      //--- change value
      state.m_xupdated=false;
      //--- function call, return result
      return(Func_lbl_55(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- check
   if(state.m_rstate.stage==14)
     {
      //--- change value
      state.m_xupdated=false;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_rstate.stage==15)
     {
      //--- change value
      state.m_xupdated=false;
      //--- return result
      return(false);
     }
//--- Routine body
//--- prepare
   n=state.m_n;
   m=state.m_m;
   state.m_repiterationscount=0;
   state.m_repterminationtype=0;
   state.m_repnfunc=0;
   state.m_repnjac=0;
   state.m_repngrad=0;
   state.m_repnhess=0;
   state.m_repncholesky=0;
//--- check consistency of constraints
//--- set constraints
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(state.m_havebndl[i] && state.m_havebndu[i])
        {
         //--- check
         if(state.m_bndl[i]>state.m_bndu[i])
           {
            state.m_repterminationtype=-3;
            //--- return result
            return(false);
           }
        }
     }
//--- function call
   CMinQP::MinQPSetBC(state.m_qpstate,state.m_bndl,state.m_bndu);
//--- Initial report of current point
//--- Note 1: we rewrite State.X twice because
//--- user may accidentally change it after first call.
//--- Note 2: we set NeedF or NeedFI depending on what
//--- information about function we have.
   if(!state.m_xrep)
      return(Func_lbl_16(state,n,m,iflag,i,k,bflag,v,s,t));
   for(i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
//--- function call
   ClearRequestFields(state);
//--- check
   if(!state.m_hasf)
     {
      //--- check
      if(!CAp::Assert(state.m_hasfi,"MinLM: internal error 2!"))
         return(false);
      //--- change values
      state.m_needfi=true;
      state.m_rstate.stage=1;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
      //--- return result
      return(true);
     }
//--- change values
   state.m_needf=true;
   state.m_rstate.stage=0;
//--- Saving state
   Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static void CMinLM::Func_lbl_rcomm(CMinLMState &state,int n,int m,int iflag,
                                   int i,int k,bool bflag,double v,
                                   double s,double t)
  {
//--- save
   state.m_rstate.ia[0]=n;
   state.m_rstate.ia[1]=m;
   state.m_rstate.ia[2]=iflag;
   state.m_rstate.ia[3]=i;
   state.m_rstate.ia[4]=k;
   state.m_rstate.ba[0]=bflag;
   state.m_rstate.ra[0]=v;
   state.m_rstate.ra[1]=s;
   state.m_rstate.ra[2]=t;
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_16(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
//--- Prepare control variables
   state.m_nu=1;
   state.m_lambdav=-CMath::m_maxrealnumber;
   state.m_modelage=state.m_maxmodelage+1;
   state.m_deltaxready=false;
   state.m_deltafready=false;
//--- Main cycle.
//--- We move through it until either:
//--- * one of the stopping conditions is met
//--- * we decide that stopping conditions are too stringent
//---   and break from cycle
   return(Func_lbl_20(state,n,m,iflag,i,k,bflag,v,s,t));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_19(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
//--- change value
   state.m_repnfunc=state.m_repnfunc+1;
//--- copy
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=2;
//--- Saving state
   Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_20(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
//--- First,we have to prepare quadratic model for our function.
//--- We use BFlag to ensure that model is prepared;
//--- if it is false at the end of this block,something went wrong.
//--- We may either calculate brand new model or update old one.
//--- Before this block we have:
//--- * State.XBase            - current position.
//--- * State.DeltaX           - if DeltaXReady is True
//--- * State.DeltaF           - if DeltaFReady is True
//--- After this block is over,we will have:
//--- * State.XBase            - base point (unchanged)
//--- * State.FBase            - F(XBase)
//--- * State.GBase            - linear term
//--- * State.QuadraticModel   - quadratic term
//--- * State.LambdaV          - current estimate for lambda
//--- We also clear DeltaXReady/DeltaFReady flags
//--- after initialization is done.
   bflag=false;
   if(!(state.m_algomode==0 || state.m_algomode==1))
      return(Func_lbl_22(state,n,m,iflag,i,k,bflag,v,s,t));
//--- Calculate f[] and Jacobian
   if(!(state.m_modelage>state.m_maxmodelage||!(state.m_deltaxready  &state.m_deltafready)))
      return(Func_lbl_24(state,n,m,iflag,i,k,bflag,v,s,t));
//--- Refresh model (using either finite differences or analytic Jacobian)
   if(state.m_algomode!=0)
     {
      //--- Obtain f[] and Jacobian
      for(int i_=0;i_<=n-1;i_++)
         state.m_x[i_]=state.m_xbase[i_];
      //--- function call
      ClearRequestFields(state);
      //--- change values
      state.m_needfij=true;
      state.m_rstate.stage=6;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
      //--- return result
      return(true);
     }
//--- Optimization using F values only.
//--- Use finite differences to estimate Jacobian.
   if(!CAp::Assert(state.m_hasfi,"MinLMIteration: internal error when estimating Jacobian (no f[])"))
      return(false);
   k=0;
//--- function call, return result
   return(Func_lbl_28(state,n,m,iflag,i,k,bflag,v,s,t));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_21(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
//--- Lambda is too large,we have to break iterations.
   state.m_repterminationtype=7;
   if(!state.m_xrep)
      return(false);
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
   state.m_f=state.m_fbase;
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=15;
//--- Saving state
   Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_22(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
//--- check
   if(state.m_algomode!=2)
      return(Func_lbl_31(state,n,m,iflag,i,k,bflag,v,s,t));
//--- check
   if(!CAp::Assert(!state.m_hasfi,"MinLMIteration: internal error (HasFI is True in Hessian-based mode)"))
      return(false);
//--- Obtain F,G,H
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_needfgh=true;
   state.m_rstate.stage=7;
//--- Saving state
   Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_24(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
//--- State.J contains Jacobian or its current approximation;
//--- refresh it using secant updates:
//--- f(x0+dx)=f(x0) + J*dx,
//--- J_new=J_old + u*h'
//--- h=x_new-x_old
//--- u=(f_new - f_old - J_old*h)/(h'h)
//--- We can explicitly generate h and u,but it is
//--- preferential to do in-place calculations. Only
//--- I-th row of J_old is needed to calculate u[I],
//--- so we can update J row by row in one pass.
//--- NOTE: we expect that State.XBase contains new point,
//--- State.FBase contains old point,State.DeltaX and
//--- State.DeltaY contain updates from last step.
   if(!CAp::Assert(state.m_deltaxready && state.m_deltafready,"MinLMIteration: uninitialized DeltaX/DeltaF"))
      return(false);
   t=0.0;
   for(int i_=0;i_<=n-1;i_++)
      t+=state.m_deltax[i_]*state.m_deltax[i_];
//--- check
   if(!CAp::Assert(t!=0.0,"MinLM: internal error (T=0)"))
      return(false);
   for(i=0;i<=m-1;i++)
     {
      //--- change value
      v=0.0;
      for(int i_=0;i_<=n-1;i_++)
         v+=state.m_j[i][i_]*state.m_deltax[i_];
      v=(state.m_deltaf[i]-v)/t;
      for(int i_=0;i_<=n-1;i_++)
         state.m_j[i].Set(i_,state.m_j[i][i_]+v*state.m_deltax[i_]);
     }
   for(int i_=0;i_<=m-1;i_++)
      state.m_fi[i_]=state.m_fibase[i_];
   for(int i_=0;i_<=m-1;i_++)
      state.m_fi[i_]=state.m_fi[i_]+state.m_deltaf[i_];
//--- Increase model age
   state.m_modelage=state.m_modelage+1;
//--- function call, return result
   return(Func_lbl_25(state,n,m,iflag,i,k,bflag,v,s,t));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_25(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
//--- Generate quadratic model:
//---     f(xbase+dx)=
//---                =(f0 + J*dx)'(f0 + J*dx)
//---                =f0^2 + dx'J'f0 + f0*J*dx + dx'J'J*dx
//---                =f0^2 + 2*f0*J*dx + dx'J'J*dx
//--- Note that we calculate 2*(J'J) instead of J'J because
//--- our quadratic model is based on Tailor decomposition,
//--- i.m_e. it has 0.5 before quadratic term.
   CAblas::RMatrixGemm(n,n,m,2.0,state.m_j,0,0,1,state.m_j,0,0,0,0.0,state.m_quadraticmodel,0,0);
   CAblas::RMatrixMVect(n,m,state.m_j,0,0,1,state.m_fi,0,state.m_gbase,0);
   for(int i_=0;i_<=n-1;i_++)
      state.m_gbase[i_]=2*state.m_gbase[i_];
//--- change value
   v=0.0;
   for(int i_=0;i_<=m-1;i_++)
      v+=state.m_fi[i_]*state.m_fi[i_];
   state.m_fbase=v;
//--- copy
   for(int i_=0;i_<=m-1;i_++)
      state.m_fibase[i_]=state.m_fi[i_];
//--- set control variables
   bflag=true;
//--- function call, return result
   return(Func_lbl_22(state,n,m,iflag,i,k,bflag,v,s,t));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_28(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
//--- check
   if(k>n-1)
     {
      //--- Calculate F(XBase)
      for(int i_=0;i_<=n-1;i_++)
         state.m_x[i_]=state.m_xbase[i_];
      //--- function call
      ClearRequestFields(state);
      //--- change values
      state.m_needfi=true;
      state.m_rstate.stage=5;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
      //--- return result
      return(true);
     }
//--- We guard X[k] from leaving [BndL,BndU].
//--- In case BndL=BndU,we assume that derivative in this direction is zero.
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
   state.m_x[k]=state.m_x[k]-state.m_s[k]*state.m_diffstep;
//--- check
   if(state.m_havebndl[k])
      state.m_x[k]=MathMax(state.m_x[k],state.m_bndl[k]);
//--- check
   if(state.m_havebndu[k])
      state.m_x[k]=MathMin(state.m_x[k],state.m_bndu[k]);
   state.m_xm1=state.m_x[k];
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_needfi=true;
   state.m_rstate.stage=3;
//--- Saving state
   Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_31(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
//--- check
   if(!CAp::Assert(bflag,"MinLM: internal integrity check failed!"))
      return(false);
//--- change values
   state.m_deltaxready=false;
   state.m_deltafready=false;
//--- If Lambda is not initialized,initialize it using quadratic model
   if(state.m_lambdav<0.0)
     {
      state.m_lambdav=0;
      for(i=0;i<=n-1;i++)
         state.m_lambdav=MathMax(state.m_lambdav,MathAbs(state.m_quadraticmodel[i][i])*CMath::Sqr(state.m_s[i]));
      state.m_lambdav=0.001*state.m_lambdav;
      //--- check
      if(state.m_lambdav==0.0)
         state.m_lambdav=1;
     }
//--- Test stopping conditions for function gradient
   if(BoundedScaledAntigradNorm(state,state.m_xbase,state.m_gbase)>state.m_epsg)
     {
      //--- Find value of Levenberg-Marquardt damping parameter which:
      //--- * leads to positive definite damped model
      //--- * within bounds specified by StpMax
      //--- * generates step which decreases function value
      //--- After this block IFlag is set to:
      //--- * -3,if constraints are infeasible
      //--- * -2,if model update is needed (either Lambda growth is too large
      //---       or step is too short,but we can't rely on model and stop iterations)
      //--- * -1,if model is fresh,Lambda have grown too large,termination is needed
      //--- *  0,if everything is OK,continue iterations
      //--- State.Nu can have any value on enter,but after exit it is set to 1.0
      iflag=-99;
      //--- function call, return result
      return(Func_lbl_39(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- check
   if(state.m_modelage!=0)
     {
      //--- Model is not fresh,we should refresh it and test
      //--- conditions once more
      state.m_modelage=state.m_maxmodelage+1;
      //--- function call, return result
      return(Func_lbl_20(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- Model is fresh,we can rely on it and terminate algorithm
   state.m_repterminationtype=4;
//--- check
   if(!state.m_xrep)
      return(false);
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
   state.m_f=state.m_fbase;
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=8;
//--- Saving state
   Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_39(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
//--- Do we need model update?
   if(state.m_modelage>0 && state.m_nu>=m_suspiciousnu)
     {
      iflag=-2;
      //--- function call, return result
      return(Func_lbl_40(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- Setup quadratic solver and solve quadratic programming problem.
//--- After problem is solved we'll try to bound step by StpMax
//--- (Lambda will be increased if step size is too large).
//--- We use BFlag variable to indicate that we have to increase Lambda.
//--- If it is False,we will try to increase Lambda and move to new iteration.
   bflag=true;
//--- function call
   CMinQP::MinQPSetStartingPointFast(state.m_qpstate,state.m_xbase);
//--- function call
   CMinQP::MinQPSetOriginFast(state.m_qpstate,state.m_xbase);
//--- function call
   CMinQP::MinQPSetLinearTermFast(state.m_qpstate,state.m_gbase);
//--- function call
   CMinQP::MinQPSetQuadraticTermFast(state.m_qpstate,state.m_quadraticmodel,true,0.0);
   for(i=0;i<=n-1;i++)
      state.m_tmp0[i]=state.m_quadraticmodel[i][i]+state.m_lambdav/CMath::Sqr(state.m_s[i]);
//--- function call
   CMinQP::MinQPRewriteDiagonal(state.m_qpstate,state.m_tmp0);
//--- function call
   CMinQP::MinQPOptimize(state.m_qpstate);
//--- function call
   CMinQP::MinQPResultsBuf(state.m_qpstate,state.m_xdir,state.m_qprep);
//--- check
   if(state.m_qprep.m_terminationtype>0)
     {
      //--- successful solution of QP problem
      for(int i_=0;i_<=n-1;i_++)
         state.m_xdir[i_]=state.m_xdir[i_]-state.m_xbase[i_];
      v=0.0;
      for(int i_=0;i_<=n-1;i_++)
         v+=state.m_xdir[i_]*state.m_xdir[i_];
      //--- check
      if(CMath::IsFinite(v))
        {
         v=MathSqrt(v);
         //--- check
         if((state.m_stpmax>0.0)&&(v>state.m_stpmax))
            bflag=false;
        }
      else
         bflag=false;
     }
   else
     {
      //--- Either problem is non-convex (increase LambdaV) or constraints are inconsistent
      if(!CAp::Assert(state.m_qprep.m_terminationtype==-3 || state.m_qprep.m_terminationtype==-5,"MinLM: unexpected completion code from QP solver"))
         return(false);
      //--- check
      if(state.m_qprep.m_terminationtype==-3)
        {
         iflag=-3;
         //--- function call, return result
         return(Func_lbl_40(state,n,m,iflag,i,k,bflag,v,s,t));
        }
      bflag=false;
     }
//--- check
   if(!bflag)
     {
      //--- Solution failed:
      //--- try to increase lambda to make matrix positive definite and continue.
      if(!IncreaseLambda(state.m_lambdav,state.m_nu))
        {
         iflag=-1;
         //--- function call, return result
         return(Func_lbl_40(state,n,m,iflag,i,k,bflag,v,s,t));
        }
      //--- function call, return result
      return(Func_lbl_39(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- Step in State.XDir and it is bounded by StpMax.
//--- We should check stopping conditions on step size here.
//--- DeltaX,which is used for secant updates,is initialized here.
//--- This code is a bit tricky because sometimes XDir<>0,but
//--- it is so small that XDir+XBase==XBase (in finite precision
//--- arithmetics). So we set DeltaX to XBase,then
//--- add XDir,and then subtract XBase to get exact value of
//--- DeltaX.
//--- Step length is estimated using DeltaX.
//--- NOTE: stopping conditions are tested
//--- for fresh models only (ModelAge=0)
   for(int i_=0;i_<=n-1;i_++)
      state.m_deltax[i_]=state.m_xbase[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_deltax[i_]=state.m_deltax[i_]+state.m_xdir[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_deltax[i_]=state.m_deltax[i_]-state.m_xbase[i_];
   state.m_deltaxready=true;
//--- change value
   v=0.0;
   for(i=0;i<=n-1;i++)
      v=v+CMath::Sqr(state.m_deltax[i]/state.m_s[i]);
   v=MathSqrt(v);
//--- check
   if(v>state.m_epsx)
      return(Func_lbl_41(state,n,m,iflag,i,k,bflag,v,s,t));
//--- check
   if(state.m_modelage!=0)
     {
      //--- Step is suspiciously short,but model is not fresh
      //--- and we can't rely on it.
      iflag=-2;
      //--- function call, return result
      return(Func_lbl_40(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- Step is too short,model is fresh and we can rely on it.
//--- Terminating.
   state.m_repterminationtype=2;
   if(!state.m_xrep)
      return(false);
//--- copy
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
   state.m_f=state.m_fbase;
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=9;
//--- Saving state
   Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_40(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
   state.m_nu=1;
//--- check
   if(!CAp::Assert(iflag>=-3 && iflag<=0,"MinLM: internal integrity check failed!"))
      return(false);
//--- check
   if(iflag==-3)
     {
      state.m_repterminationtype=-3;
      //--- return result
      return(false);
     }
//--- check
   if(iflag==-2)
     {
      state.m_modelage=state.m_maxmodelage+1;
      //--- function call, return result
      return(Func_lbl_20(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- check
   if(iflag==-1)
      return(Func_lbl_21(state,n,m,iflag,i,k,bflag,v,s,t));
//--- Levenberg-Marquardt step is ready.
//--- Compare predicted vs. actual decrease and decide what to do with lambda.
//--- NOTE: we expect that State.DeltaX contains direction of step,
//--- State.F contains function value at new point.
   if(!CAp::Assert(state.m_deltaxready,"MinLM: deltaX is not ready"))
      return(false);
   t=0;
   for(i=0;i<=n-1;i++)
     {
      //--- change values
      v=0.0;
      for(int i_=0;i_<=n-1;i_++)
         v+=state.m_quadraticmodel[i][i_]*state.m_deltax[i_];
      t=t+state.m_deltax[i]*state.m_gbase[i]+0.5*state.m_deltax[i]*v;
     }
//--- change values
   state.m_predicteddecrease=-t;
   state.m_actualdecrease=-(state.m_f-state.m_fbase);
//--- check
   if(state.m_predicteddecrease<=0.0)
      return(Func_lbl_21(state,n,m,iflag,i,k,bflag,v,s,t));
   v=state.m_actualdecrease/state.m_predicteddecrease;
//--- check
   if(v>=0.1)
      return(Func_lbl_49(state,n,m,iflag,i,k,bflag,v,s,t));
//--- check
   if(IncreaseLambda(state.m_lambdav,state.m_nu))
      return(Func_lbl_49(state,n,m,iflag,i,k,bflag,v,s,t));
//--- Lambda is too large,we have to break iterations.
   state.m_repterminationtype=7;
//--- check
   if(!state.m_xrep)
      return(false);
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
   state.m_f=state.m_fbase;
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=12;
//--- Saving state
   Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_41(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
//--- Let's evaluate new step:
//--- a) if we have Fi vector,we evaluate it using rcomm,and
//---    then we manually calculate State.F as sum of squares of Fi[]
//--- b) if we have F value,we just evaluate it through rcomm interface
//--- We prefer (a) because we may need Fi vector for additional
//--- iterations
   if(!CAp::Assert(state.m_hasfi|state.m_hasf,"MinLM: internal error 2!"))
      return(false);
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_x[i_]+state.m_xdir[i_];
//--- function call
   ClearRequestFields(state);
//--- check
   if(!state.m_hasfi)
     {
      //--- change values
      state.m_needf=true;
      state.m_rstate.stage=11;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
      //--- return result
      return(true);
     }
//--- change values
   state.m_needfi=true;
   state.m_rstate.stage=10;
//--- Saving state
   Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_48(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
   state.m_repnfunc=state.m_repnfunc+1;
//--- check
   if(state.m_f>=state.m_fbase)
     {
      //--- Increase lambda and continue
      if(!IncreaseLambda(state.m_lambdav,state.m_nu))
        {
         iflag=-1;
         //--- function call, return result
         return(Func_lbl_40(state,n,m,iflag,i,k,bflag,v,s,t));
        }
      //--- function call, return result
      return(Func_lbl_39(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- We've found our step!
   iflag=0;
//--- function call, return result
   return(Func_lbl_40(state,n,m,iflag,i,k,bflag,v,s,t));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_49(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
//--- check
   if(v>0.5)
      DecreaseLambda(state.m_lambdav,state.m_nu);
//--- Accept step,report it and
//--- test stopping conditions on iterations count and function decrease.
//--- NOTE: we expect that State.DeltaX contains direction of step,
//--- State.F contains function value at new point.
//--- NOTE2: we should update XBase ONLY. In the beginning of the next
//--- iteration we expect that State.FIBase is NOT updated and
//--- contains old value of a function vector.
   for(int i_=0;i_<=n-1;i_++)
      state.m_xbase[i_]=state.m_xbase[i_]+state.m_deltax[i_];
//--- check
   if(!state.m_xrep)
      return(Func_lbl_55(state,n,m,iflag,i,k,bflag,v,s,t));
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=13;
//--- Saving state
   Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinLMIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CMinLM::Func_lbl_55(CMinLMState &state,int &n,int &m,int &iflag,
                                int &i,int &k,bool &bflag,double &v,
                                double &s,double &t)
  {
   state.m_repiterationscount=state.m_repiterationscount+1;
//--- check
   if(state.m_repiterationscount>=state.m_maxits&&state.m_maxits>0)
      state.m_repterminationtype=5;
//--- check
   if(state.m_modelage==0)
     {
      //--- check
      if(MathAbs(state.m_f-state.m_fbase)<=state.m_epsf*MathMax(1,MathMax(MathAbs(state.m_f),MathAbs(state.m_fbase))))
         state.m_repterminationtype=1;
     }
//--- check
   if(state.m_repterminationtype<=0)
     {
      state.m_modelage=state.m_modelage+1;
      //--- function call, return result
      return(Func_lbl_20(state,n,m,iflag,i,k,bflag,v,s,t));
     }
//--- check
   if(!state.m_xrep)
      return(false);
//--- Report: XBase contains new point,F contains function value at new point
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=14;
//--- Saving state
   Func_lbl_rcomm(state,n,m,iflag,i,k,bflag,v,s,t);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CMinComp                                     |
//+------------------------------------------------------------------+
class CMinASAState
  {
public:
   //--- variables
   int               m_n;
   double            m_epsg;
   double            m_epsf;
   double            m_epsx;
   int               m_maxits;
   bool              m_xrep;
   double            m_stpmax;
   int               m_cgtype;
   int               m_k;
   int               m_nfev;
   int               m_mcstage;
   int               m_curalgo;
   int               m_acount;
   double            m_mu;
   double            m_finit;
   double            m_dginit;
   double            m_fold;
   double            m_stp;
   double            m_laststep;
   double            m_f;
   bool              m_needfg;
   bool              m_xupdated;
   RCommState        m_rstate;
   int               m_repiterationscount;
   int               m_repnfev;
   int               m_repterminationtype;
   int               m_debugrestartscount;
   CLinMinState      m_lstate;
   double            m_betahs;
   double            m_betady;
   //--- arrays
   double            m_bndl[];
   double            m_bndu[];
   double            m_ak[];
   double            m_xk[];
   double            m_dk[];
   double            m_an[];
   double            m_xn[];
   double            m_dn[];
   double            m_d[];
   double            m_work[];
   double            m_yk[];
   double            m_gc[];
   double            m_x[];
   double            m_g[];
   //--- constructor, destructor
                     CMinASAState(void);
                    ~CMinASAState(void);
   //--- copy
   void              Copy(CMinASAState &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinASAState::CMinASAState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinASAState::~CMinASAState(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMinASAState::Copy(CMinASAState &obj)
  {
//--- copy variables
   m_n=obj.m_n;
   m_epsg=obj.m_epsg;
   m_epsf=obj.m_epsf;
   m_epsx=obj.m_epsx;
   m_maxits=obj.m_maxits;
   m_xrep=obj.m_xrep;
   m_stpmax=obj.m_stpmax;
   m_cgtype=obj.m_cgtype;
   m_k=obj.m_k;
   m_nfev=obj.m_nfev;
   m_mcstage=obj.m_mcstage;
   m_curalgo=obj.m_curalgo;
   m_acount=obj.m_acount;
   m_mu=obj.m_mu;
   m_finit=obj.m_finit;
   m_dginit=obj.m_dginit;
   m_fold=obj.m_fold;
   m_stp=obj.m_stp;
   m_laststep=obj.m_laststep;
   m_f=obj.m_f;
   m_needfg=obj.m_needfg;
   m_xupdated=obj.m_xupdated;
   m_repiterationscount=obj.m_repiterationscount;
   m_repnfev=obj.m_repnfev;
   m_repterminationtype=obj.m_repterminationtype;
   m_debugrestartscount=obj.m_debugrestartscount;
   m_betahs=obj.m_betahs;
   m_betady=obj.m_betady;
   m_rstate.Copy(obj.m_rstate);
   m_lstate.Copy(obj.m_lstate);
//--- copy arrays
   ArrayCopy(m_bndl,obj.m_bndl);
   ArrayCopy(m_bndu,obj.m_bndu);
   ArrayCopy(m_ak,obj.m_ak);
   ArrayCopy(m_xk,obj.m_xk);
   ArrayCopy(m_dk,obj.m_dk);
   ArrayCopy(m_an,obj.m_an);
   ArrayCopy(m_xn,obj.m_xn);
   ArrayCopy(m_dn,obj.m_dn);
   ArrayCopy(m_d,obj.m_d);
   ArrayCopy(m_work,obj.m_work);
   ArrayCopy(m_yk,obj.m_yk);
   ArrayCopy(m_gc,obj.m_gc);
   ArrayCopy(m_x,obj.m_x);
   ArrayCopy(m_g,obj.m_g);
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CMinASAState                     |
//+------------------------------------------------------------------+
class CMinASAStateShell
  {
private:
   CMinASAState      m_innerobj;
public:
   //--- constructors, destructor
                     CMinASAStateShell(void);
                     CMinASAStateShell(CMinASAState &obj);
                    ~CMinASAStateShell(void);
   //--- methods
   bool              GetNeedFG(void);
   void              SetNeedFG(const bool b);
   bool              GetXUpdated(void);
   void              SetXUpdated(const bool b);
   double            GetF(void);
   void              SetF(const double d);
   CMinASAState     *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinASAStateShell::CMinASAStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMinASAStateShell::CMinASAStateShell(CMinASAState &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinASAStateShell::~CMinASAStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needfg                         |
//+------------------------------------------------------------------+
bool CMinASAStateShell::GetNeedFG(void)
  {
//--- return result
   return(m_innerobj.m_needfg);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needfg                        |
//+------------------------------------------------------------------+
void CMinASAStateShell::SetNeedFG(const bool b)
  {
//--- change value
   m_innerobj.m_needfg=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable xupdated                       |
//+------------------------------------------------------------------+
bool CMinASAStateShell::GetXUpdated(void)
  {
//--- return result
   return(m_innerobj.m_xupdated);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable xupdated                      |
//+------------------------------------------------------------------+
void CMinASAStateShell::SetXUpdated(const bool b)
  {
//--- change value
   m_innerobj.m_xupdated=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable f                              |
//+------------------------------------------------------------------+
double CMinASAStateShell::GetF(void)
  {
//--- return result
   return(m_innerobj.m_f);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable f                             |
//+------------------------------------------------------------------+
void CMinASAStateShell::SetF(const double d)
  {
//--- change value
   m_innerobj.m_f=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMinASAState *CMinASAStateShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CMinComp                                     |
//+------------------------------------------------------------------+
class CMinASAReport
  {
public:
   //--- variables
   int               m_iterationscount;
   int               m_nfev;
   int               m_terminationtype;
   int               m_activeconstraints;
   //--- constructor, destructor
                     CMinASAReport(void);
                    ~CMinASAReport(void);
   //--- copy
   void              Copy(CMinASAReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinASAReport::CMinASAReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinASAReport::~CMinASAReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMinASAReport::Copy(CMinASAReport &obj)
  {
//--- copy variables
   m_iterationscount=obj.m_iterationscount;
   m_nfev=obj.m_nfev;
   m_terminationtype=obj.m_terminationtype;
   m_activeconstraints=obj.m_activeconstraints;
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CMinASAReport                    |
//+------------------------------------------------------------------+
class CMinASAReportShell
  {
private:
   CMinASAReport     m_innerobj;
public:
   //--- constructors, destructor
                     CMinASAReportShell(void);
                     CMinASAReportShell(CMinASAReport &obj);
                    ~CMinASAReportShell(void);
   //--- methods
   int               GetIterationsCount(void);
   void              SetIterationsCount(const int i);
   int               GetNFev(void);
   void              SetNFev(const int i);
   int               GetTerminationType(void);
   void              SetTerminationType(const int i);
   int               GetActiveConstraints(void);
   void              SetActiveConstraints(const int i);
   CMinASAReport    *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinASAReportShell::CMinASAReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMinASAReportShell::CMinASAReportShell(CMinASAReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinASAReportShell::~CMinASAReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable iterationscount                |
//+------------------------------------------------------------------+
int CMinASAReportShell::GetIterationsCount(void)
  {
//--- return result
   return(m_innerobj.m_iterationscount);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable iterationscount               |
//+------------------------------------------------------------------+
void CMinASAReportShell::SetIterationsCount(const int i)
  {
//--- change value
   m_innerobj.m_iterationscount=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable nfev                           |
//+------------------------------------------------------------------+
int CMinASAReportShell::GetNFev(void)
  {
//--- return result
   return(m_innerobj.m_nfev);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable nfev                          |
//+------------------------------------------------------------------+
void CMinASAReportShell::SetNFev(const int i)
  {
//--- change value
   m_innerobj.m_nfev=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable terminationtype                |
//+------------------------------------------------------------------+
int CMinASAReportShell::GetTerminationType(void)
  {
//--- return result
   return(m_innerobj.m_terminationtype);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable terminationtype               |
//+------------------------------------------------------------------+
void CMinASAReportShell::SetTerminationType(const int i)
  {
//--- change value
   m_innerobj.m_terminationtype=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable activeconstraints              |
//+------------------------------------------------------------------+
int CMinASAReportShell::GetActiveConstraints(void)
  {
//--- return result
   return(m_innerobj.m_activeconstraints);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable activeconstraints             |
//+------------------------------------------------------------------+
void CMinASAReportShell::SetActiveConstraints(const int i)
  {
//--- change value
   m_innerobj.m_activeconstraints=i;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMinASAReport *CMinASAReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Backward compatibility functions                                 |
//+------------------------------------------------------------------+
class CMinComp
  {
private:
   //--- private methods
   static double     ASABoundedAntigradNorm(CMinASAState &state);
   static double     ASAGINorm(CMinASAState &state);
   static double     ASAD1Norm(CMinASAState &state);
   static bool       ASAUIsEmpty(CMinASAState &state);
   static void       ClearRequestFields(CMinASAState &state);
   //--- auxiliary functions for MinASAIteration
   static void       Func_lbl_rcomm(CMinASAState &state,int n,int i,int mcinfo,int diffcnt,bool b,bool stepfound,double betak,double v,double vv);
   static bool       Func_lbl_15(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_17(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_19(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_21(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_24(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_26(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_27(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_29(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_31(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_35(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_39(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_43(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_49(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_51(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_52(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_53(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_55(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_59(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_63(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
   static bool       Func_lbl_65(CMinASAState &state,int &n,int &i,int &mcinfo,int &diffcnt,bool &b,bool &stepfound,double &betak,double &v,double &vv);
public:
   //--- class constants
   static const int  m_n1;
   static const int  m_n2;
   static const double m_stpmin;
   static const double m_gtol;
   static const double m_gpaftol;
   static const double m_gpadecay;
   static const double m_asarho;
   //--- constructor, destructor
                     CMinComp(void);
                    ~CMinComp(void);
   //--- public methods
   static void       MinLBFGSSetDefaultPreconditioner(CMinLBFGSState &state);
   static void       MinLBFGSSetCholeskyPreconditioner(CMinLBFGSState &state,CMatrixDouble &p,const bool isupper);
   static void       MinBLEICSetBarrierWidth(CMinBLEICState &state,const double mu);
   static void       MinBLEICSetBarrierDecay(CMinBLEICState &state,const double mudecay);
   static void       MinASACreate(const int n,double &x[],double &bndl[],double &bndu[],CMinASAState &state);
   static void       MinASASetCond(CMinASAState &state,const double epsg,const double epsf,double epsx,const int maxits);
   static void       MinASASetXRep(CMinASAState &state,const bool needxrep);
   static void       MinASASetAlgorithm(CMinASAState &state,int algotype);
   static void       MinASASetStpMax(CMinASAState &state,const double stpmax);
   static void       MinASAResults(CMinASAState &state,double &x[],CMinASAReport &rep);
   static void       MinASAResultsBuf(CMinASAState &state,double &x[],CMinASAReport &rep);
   static void       MinASARestartFrom(CMinASAState &state,double &x[],double &bndl[],double &bndu[]);
   static bool       MinASAIteration(CMinASAState &state);
  };
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const int    CMinComp::m_n1=2;
const int    CMinComp::m_n2=2;
const double CMinComp::m_stpmin=1.0E-300;
const double CMinComp::m_gtol=0.3;
const double CMinComp::m_gpaftol=0.0001;
const double CMinComp::m_gpadecay=0.5;
const double CMinComp::m_asarho=0.5;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMinComp::CMinComp(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMinComp::~CMinComp(void)
  {

  }
//+------------------------------------------------------------------+
//| Obsolete function, use MinLBFGSSetPrecDefault() instead.         |
//+------------------------------------------------------------------+
static void CMinComp::MinLBFGSSetDefaultPreconditioner(CMinLBFGSState &state)
  {
//--- function call
   CMinLBFGS::MinLBFGSSetPrecDefault(state);
  }
//+------------------------------------------------------------------+
//| Obsolete function, use MinLBFGSSetCholeskyPreconditioner()       |
//| instead.                                                         |
//+------------------------------------------------------------------+
static void CMinComp::MinLBFGSSetCholeskyPreconditioner(CMinLBFGSState &state,
                                                        CMatrixDouble &p,
                                                        const bool isupper)
  {
//--- function call
   CMinLBFGS::MinLBFGSSetPrecCholesky(state,p,isupper);
  }
//+------------------------------------------------------------------+
//| This is obsolete function which was used by previous version of  |
//| the  BLEIC optimizer. It does nothing in the current version of  |
//| BLEIC.                                                           |
//+------------------------------------------------------------------+
static void CMinComp::MinBLEICSetBarrierWidth(CMinBLEICState &state,
                                              const double mu)
  {

  }
//+------------------------------------------------------------------+
//| This is obsolete function which was used by previous version of  |
//| the  BLEIC optimizer. It does nothing in the current version of  |
//| BLEIC.                                                           |
//+------------------------------------------------------------------+
static void CMinComp::MinBLEICSetBarrierDecay(CMinBLEICState &state,
                                              const double mudecay)
  {

  }
//+------------------------------------------------------------------+
//| Obsolete optimization algorithm.                                 |
//| Was replaced by MinBLEIC subpackage.                             |
//+------------------------------------------------------------------+
static void CMinComp::MinASACreate(const int n,double &x[],double &bndl[],
                                   double &bndu[],CMinASAState &state)
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N too small!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(bndl)>=n,__FUNCTION__+": Length(BndL)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(bndl,n),__FUNCTION__+": BndL contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(bndu)>=n,__FUNCTION__+": Length(BndU)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(bndu,n),__FUNCTION__+": BndU contains infinite or NaN values!"))
      return;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(!CAp::Assert((double)(bndl[i])<=(double)(bndu[i]),__FUNCTION__+": inconsistent bounds!"))
         return;
      //--- check
      if(!CAp::Assert((double)(bndl[i])<=x[i],__FUNCTION__+": infeasible X!"))
         return;
      //--- check
      if(!CAp::Assert(x[i]<=(double)(bndu[i]),__FUNCTION__+": infeasible X!"))
         return;
     }
//--- Initialize
   state.m_n=n;
   MinASASetCond(state,0,0,0,0);
   MinASASetXRep(state,false);
   MinASASetStpMax(state,0);
   MinASASetAlgorithm(state,-1);
//--- allocation
   ArrayResizeAL(state.m_bndl,n);
   ArrayResizeAL(state.m_bndu,n);
   ArrayResizeAL(state.m_ak,n);
   ArrayResizeAL(state.m_xk,n);
   ArrayResizeAL(state.m_dk,n);
   ArrayResizeAL(state.m_an,n);
   ArrayResizeAL(state.m_xn,n);
   ArrayResizeAL(state.m_dn,n);
   ArrayResizeAL(state.m_x,n);
   ArrayResizeAL(state.m_d,n);
   ArrayResizeAL(state.m_g,n);
   ArrayResizeAL(state.m_gc,n);
   ArrayResizeAL(state.m_work,n);
   ArrayResizeAL(state.m_yk,n);
//--- function call
   MinASARestartFrom(state,x,bndl,bndu);
  }
//+------------------------------------------------------------------+
//| Obsolete optimization algorithm.                                 |
//| Was replaced by MinBLEIC subpackage.                             |
//+------------------------------------------------------------------+
static void CMinComp::MinASASetCond(CMinASAState &state,const double epsg,
                                    const double epsf,double epsx,
                                    const int maxits)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsg),__FUNCTION__+": EpsG is not finite number!"))
      return;
//--- check
   if(!CAp::Assert(epsg>=0.0,__FUNCTION__+": negative EpsG!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsf),__FUNCTION__+": EpsF is not finite number!"))
      return;
//--- check
   if(!CAp::Assert(epsf>=0.0,__FUNCTION__+": negative EpsF!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsx),__FUNCTION__+": EpsX is not finite number!"))
      return;
//--- check
   if(!CAp::Assert(epsx>=0.0,__FUNCTION__+": negative EpsX!"))
      return;
//--- check
   if(!CAp::Assert(maxits>=0,__FUNCTION__+": negative MaxIts!"))
      return;
//--- check
   if(((epsg==0.0 && epsf==0.0) && epsx==0.0) && maxits==0)
      epsx=1.0E-6;
//--- change values
   state.m_epsg=epsg;
   state.m_epsf=epsf;
   state.m_epsx=epsx;
   state.m_maxits=maxits;
  }
//+------------------------------------------------------------------+
//| Obsolete optimization algorithm.                                 |
//| Was replaced by MinBLEIC subpackage.                             |
//+------------------------------------------------------------------+
static void CMinComp::MinASASetXRep(CMinASAState &state,const bool needxrep)
  {
//--- change value
   state.m_xrep=needxrep;
  }
//+------------------------------------------------------------------+
//| Obsolete optimization algorithm.                                 |
//| Was replaced by MinBLEIC subpackage.                             |
//+------------------------------------------------------------------+
static void CMinComp::MinASASetAlgorithm(CMinASAState &state,int algotype)
  {
//--- check
   if(!CAp::Assert(algotype>=-1 && algotype<=1,__FUNCTION__+": incorrect AlgoType!"))
      return;
//--- check
   if(algotype==-1)
      algotype=1;
//--- change value
   state.m_cgtype=algotype;
  }
//+------------------------------------------------------------------+
//| Obsolete optimization algorithm.                                 |
//| Was replaced by MinBLEIC subpackage.                             |
//+------------------------------------------------------------------+
static void CMinComp::MinASASetStpMax(CMinASAState &state,const double stpmax)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(stpmax),__FUNCTION__+": StpMax is not finite!"))
      return;
//--- check
   if(!CAp::Assert(stpmax>=0.0,__FUNCTION__+": StpMax<0!"))
      return;
//--- change value
   state.m_stpmax=stpmax;
  }
//+------------------------------------------------------------------+
//| Obsolete optimization algorithm.                                 |
//| Was replaced by MinBLEIC subpackage.                             |
//+------------------------------------------------------------------+
static void CMinComp::MinASAResults(CMinASAState &state,double &x[],CMinASAReport &rep)
  {
//--- reset memory
   ArrayResizeAL(x,0);
//--- function call
   MinASAResultsBuf(state,x,rep);
  }
//+------------------------------------------------------------------+
//| Obsolete optimization algorithm.                                 |
//| Was replaced by MinBLEIC subpackage.                             |
//+------------------------------------------------------------------+
static void CMinComp::MinASAResultsBuf(CMinASAState &state,double &x[],
                                       CMinASAReport &rep)
  {
//--- create variables
   int i=0;
   int i_=0;
//--- check
   if(CAp::Len(x)<state.m_n)
      ArrayResizeAL(x,state.m_n);
//--- copy
   for(i_=0;i_<=state.m_n-1;i_++)
      x[i_]=state.m_x[i_];
//--- change values
   rep.m_iterationscount=state.m_repiterationscount;
   rep.m_nfev=state.m_repnfev;
   rep.m_terminationtype=state.m_repterminationtype;
   rep.m_activeconstraints=0;
   for(i=0;i<=state.m_n-1;i++)
     {
      //--- check
      if(state.m_ak[i]==0.0)
         rep.m_activeconstraints=rep.m_activeconstraints+1;
     }
  }
//+------------------------------------------------------------------+
//| Obsolete optimization algorithm.                                 |
//| Was replaced by MinBLEIC subpackage.                             |
//+------------------------------------------------------------------+
static void CMinComp::MinASARestartFrom(CMinASAState &state,double &x[],double &bndl[],double &bndu[])
  {
//--- create a variable
   int i_=0;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=state.m_n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,state.m_n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(bndl)>=state.m_n,__FUNCTION__+": Length(BndL)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(bndl,state.m_n),__FUNCTION__+": BndL contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(bndu)>=state.m_n,__FUNCTION__+": Length(BndU)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(bndu,state.m_n),__FUNCTION__+": BndU contains infinite or NaN values!"))
      return;
//--- copy
   for(i_=0;i_<=state.m_n-1;i_++)
      state.m_x[i_]=x[i_];
   for(i_=0;i_<=state.m_n-1;i_++)
      state.m_bndl[i_]=bndl[i_];
   for(i_=0;i_<=state.m_n-1;i_++)
      state.m_bndu[i_]=bndu[i_];
   state.m_laststep=0;
//--- allocation
   ArrayResizeAL(state.m_rstate.ia,4);
   ArrayResizeAL(state.m_rstate.ba,2);
   ArrayResizeAL(state.m_rstate.ra,3);
   state.m_rstate.stage=-1;
//--- function call
   ClearRequestFields(state);
  }
//+------------------------------------------------------------------+
//| Returns norm of bounded anti-gradient.                           |
//| Bounded antigradient is a vector obtained from anti-gradient by  |
//| zeroing components which point outwards:                         |
//|     result = norm(v)                                             |
//|     v[i]=0     if ((-g[i]<0)and(x[i]=bndl[i])) or                |
//|                   ((-g[i]>0)and(x[i]=bndu[i]))                   |
//|     v[i]=-g[i] otherwise                                         |
//| This function may be used to check a stopping criterion.         |
//+------------------------------------------------------------------+
static double CMinComp::ASABoundedAntigradNorm(CMinASAState &state)
  {
//--- create variables
   double result=0;
   int    i=0;
   double v=0;
//--- initialization
   result=0;
   for(i=0;i<=state.m_n-1;i++)
     {
      v=-state.m_g[i];
      //--- check
      if(state.m_x[i]==state.m_bndl[i] && -state.m_g[i]<0.0)
         v=0;
      //--- check
      if(state.m_x[i]==state.m_bndu[i] && -state.m_g[i]>0.0)
         v=0;
      result=result+CMath::Sqr(v);
     }
//--- return result
   return(MathSqrt(result));
  }
//+------------------------------------------------------------------+
//| Returns norm of GI(x).                                           |
//| GI(x) is a gradient vector whose components associated with      |
//| active constraints are zeroed. It differs from bounded           |
//| anti-gradient because components of GI(x) are zeroed             |
//| independently of sign(g[i]), and anti-gradient's components are  |
//| zeroed with respect to both constraint and sign.                 |
//+------------------------------------------------------------------+
static double CMinComp::ASAGINorm(CMinASAState &state)
  {
//--- create variables
   double result=0;
   int    i=0;
//--- initialization
   result=0;
   for(i=0;i<=state.m_n-1;i++)
     {
      //--- check
      if(state.m_x[i]!=state.m_bndl[i] && state.m_x[i]!=state.m_bndu[i])
         result=result+CMath::Sqr(state.m_g[i]);
     }
//--- return result
   return(MathSqrt(result));
  }
//+------------------------------------------------------------------+
//| Returns norm(D1(State.X))                                        |
//| For a meaning of D1 see 'NEW ACTIVE SET ALGORITHM FOR BOX        |
//| CONSTRAINED OPTIMIZATION' by WILLIAM W. HAGER AND HONGCHAO ZHANG.|
//+------------------------------------------------------------------+
static double CMinComp::ASAD1Norm(CMinASAState &state)
  {
//--- create variables
   double result=0;
   int    i=0;
//--- initialization
   result=0;
   for(i=0;i<=state.m_n-1;i++)
      result=result+CMath::Sqr(CApServ::BoundVal(state.m_x[i]-state.m_g[i],state.m_bndl[i],state.m_bndu[i])-state.m_x[i]);
//--- return result
   return(MathSqrt(result));
  }
//+------------------------------------------------------------------+
//| Returns True, if U set is empty.                                 |
//| * State.X is used as point,                                      |
//| * State.G - as gradient,                                         |
//| * D is calculated within function (because State.D may have      |
//|   different meaning depending on current optimization algorithm) |
//| For a meaning of U see 'NEW ACTIVE SET ALGORITHM FOR BOX         |
//| CONSTRAINED OPTIMIZATION' by WILLIAM W. HAGER AND HONGCHAO ZHANG.|
//+------------------------------------------------------------------+
static bool CMinComp::ASAUIsEmpty(CMinASAState &state)
  {
//--- create variables
   int    i=0;
   double d=0;
   double d2=0;
   double d32=0;
//--- initialization
   d=ASAD1Norm(state);
   d2=MathSqrt(d);
   d32=d*d2;
   for(i=0;i<=state.m_n-1;i++)
     {
      //--- check
      if(MathAbs(state.m_g[i])>=d2 && MathMin(state.m_x[i]-state.m_bndl[i],state.m_bndu[i]-state.m_x[i])>=d32)
         return(false);
     }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Clears request fileds (to be sure that we don't forgot to clear  |
//| something)                                                       |
//+------------------------------------------------------------------+
static void CMinComp::ClearRequestFields(CMinASAState &state)
  {
//--- change values
   state.m_needfg=false;
   state.m_xupdated=false;
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
static bool CMinComp::MinASAIteration(CMinASAState &state)
  {
//--- create variables
   int    n=0;
   int    i=0;
   double betak=0;
   double v=0;
   double vv=0;
   int    mcinfo=0;
   bool   b;
   bool   stepfound;
   int    diffcnt=0;
   int    i_=0;
//--- This code initializes locals by:
//--- * random values determined during code
//---   generation - on first subroutine call
//--- * values from previous call - on subsequent calls
   if(state.m_rstate.stage>=0)
     {
      //--- initialization
      n=state.m_rstate.ia[0];
      i=state.m_rstate.ia[1];
      mcinfo=state.m_rstate.ia[2];
      diffcnt=state.m_rstate.ia[3];
      b=state.m_rstate.ba[0];
      stepfound=state.m_rstate.ba[1];
      betak=state.m_rstate.ra[0];
      v=state.m_rstate.ra[1];
      vv=state.m_rstate.ra[2];
     }
   else
     {
      //--- initialization
      n=-983;
      i=-989;
      mcinfo=-834;
      diffcnt=900;
      b=true;
      stepfound=false;
      betak=214;
      v=-338;
      vv=-686;
     }
//--- check
   if(state.m_rstate.stage==0)
     {
      //--- change value
      state.m_needfg=false;
      //--- check
      if(!state.m_xrep)
         return(Func_lbl_15(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
      //--- progress report
      ClearRequestFields(state);
      state.m_xupdated=true;
      state.m_rstate.stage=1;
      //--- Saving state
      Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==1)
     {
      //--- change value
      state.m_xupdated=false;
      //--- function call, return result
      return(Func_lbl_15(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==2)
     {
      //--- change values
      state.m_needfg=false;
      state.m_repnfev=state.m_repnfev+1;
      stepfound=state.m_f<=state.m_finit+m_gpaftol*state.m_dginit;
      //--- function call, return result
      return(Func_lbl_24(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==3)
     {
      //--- change values
      state.m_needfg=false;
      state.m_repnfev=state.m_repnfev+1;
      //--- check
      if(state.m_stp<=m_stpmin)
        {
         for(i_=0;i_<=n-1;i_++)
            state.m_xn[i_]=state.m_x[i_];
         //--- function call, return result
         return(Func_lbl_26(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
        }
      //--- check
      if(state.m_f<=state.m_finit+state.m_stp*m_gpaftol*state.m_dginit)
        {
         //--- copy
         for(i_=0;i_<=n-1;i_++)
            state.m_xn[i_]=state.m_x[i_];
         //--- function call, return result
         return(Func_lbl_26(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
        }
      //--- change value
      state.m_stp=state.m_stp*m_gpadecay;
      //--- function call, return result
      return(Func_lbl_27(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==4)
     {
      //--- change value
      state.m_xupdated=false;
      //--- function call, return result
      return(Func_lbl_29(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==5)
     {
      //--- change value
      state.m_xupdated=false;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_rstate.stage==6)
     {
      //--- change value
      state.m_xupdated=false;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_rstate.stage==7)
     {
      //--- change value
      state.m_xupdated=false;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_rstate.stage==8)
     {
      //--- change value
      state.m_xupdated=false;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_rstate.stage==9)
     {
      //--- change value
      state.m_needfg=false;
      //--- postprocess data: zero components of G corresponding to
      //--- the active constraints
      for(i=0;i<=n-1;i++)
        {
         //--- check
         if(state.m_x[i]==state.m_bndl[i] || state.m_x[i]==state.m_bndu[i])
            state.m_gc[i]=0;
         else
            state.m_gc[i]=state.m_g[i];
        }
      CLinMin::MCSrch(n,state.m_xn,state.m_f,state.m_gc,state.m_d,state.m_stp,state.m_stpmax,m_gtol,mcinfo,state.m_nfev,state.m_work,state.m_lstate,state.m_mcstage);
      //--- function call, return result
      return(Func_lbl_51(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==10)
     {
      //--- change value
      state.m_xupdated=false;
      //--- function call, return result
      return(Func_lbl_53(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
     }
//--- check
   if(state.m_rstate.stage==11)
     {
      //--- change value
      state.m_xupdated=false;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_rstate.stage==12)
     {
      //--- change value
      state.m_xupdated=false;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_rstate.stage==13)
     {
      //--- change value
      state.m_xupdated=false;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_rstate.stage==14)
     {
      //--- change value
      state.m_xupdated=false;
      //--- return result
      return(false);
     }
//--- Routine body
//--- Prepare
   n=state.m_n;
   state.m_repterminationtype=0;
   state.m_repiterationscount=0;
   state.m_repnfev=0;
   state.m_debugrestartscount=0;
   state.m_cgtype=1;
//--- copy
   for(i_=0;i_<=n-1;i_++)
      state.m_xk[i_]=state.m_x[i_];
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(state.m_xk[i]==state.m_bndl[i] || state.m_xk[i]==state.m_bndu[i])
         state.m_ak[i]=0;
      else
         state.m_ak[i]=1;
     }
//--- change values
   state.m_mu=0.1;
   state.m_curalgo=0;
//--- Calculate F/G,initialize algorithm
   ClearRequestFields(state);
//--- change values
   state.m_needfg=true;
   state.m_rstate.stage=0;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static void CMinComp::Func_lbl_rcomm(CMinASAState &state,int n,int i,
                                     int mcinfo,int diffcnt,bool b,
                                     bool stepfound,double betak,
                                     double v,double vv)
  {
//--- save
   state.m_rstate.ia[0]=n;
   state.m_rstate.ia[1]=i;
   state.m_rstate.ia[2]=mcinfo;
   state.m_rstate.ia[3]=diffcnt;
   state.m_rstate.ba[0]=b;
   state.m_rstate.ba[1]=stepfound;
   state.m_rstate.ra[0]=betak;
   state.m_rstate.ra[1]=v;
   state.m_rstate.ra[2]=vv;
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_15(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- check
   if(ASABoundedAntigradNorm(state)<=state.m_epsg)
     {
      state.m_repterminationtype=4;
      //--- return result
      return(false);
     }
   state.m_repnfev=state.m_repnfev+1;
//--- Main cycle
//--- At the beginning of new iteration:
//--- * CurAlgo stores current algorithm selector
//--- * State.XK,State.F and State.G store current X/F/G
//--- * State.AK stores current set of active constraints
   return(Func_lbl_17(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_17(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- GPA algorithm
   if(state.m_curalgo!=0)
      return(Func_lbl_19(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- change values
   state.m_k=0;
   state.m_acount=0;
//--- function call, return result
   return(Func_lbl_21(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_19(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- CG algorithm
   if(state.m_curalgo!=1)
      return(Func_lbl_17(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- first,check that there are non-active constraints.
//--- move to GPA algorithm,if all constraints are active
   b=true;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(state.m_ak[i]!=0.0)
        {
         b=false;
         break;
        }
     }
//--- check
   if(b)
     {
      state.m_curalgo=0;
      //--- function call, return result
      return(Func_lbl_17(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
     }
//--- CG iterations
   state.m_fold=state.m_f;
   for(int i_=0;i_<=n-1;i_++)
      state.m_xk[i_]=state.m_x[i_];
   for(i=0;i<=n-1;i++)
     {
      //--- change values
      state.m_dk[i]=-(state.m_g[i]*state.m_ak[i]);
      state.m_gc[i]=state.m_g[i]*state.m_ak[i];
     }
//--- function call, return result
   return(Func_lbl_49(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_21(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- Determine Dk=proj(xk - gk)-xk
   for(i=0;i<=n-1;i++)
      state.m_d[i]=CApServ::BoundVal(state.m_xk[i]-state.m_g[i],state.m_bndl[i],state.m_bndu[i])-state.m_xk[i];
//--- Armijo line search.
//--- * exact search with alpha=1 is tried first,
//---   'exact' means that we evaluate f() EXACTLY at
//---   bound(x-g,bndl,bndu),without intermediate floating
//---   point operations.
//--- * alpha<1 are tried if explicit search wasn't successful
//--- Result is placed into XN.
//--- Two types of search are needed because we can't
//--- just use second type with alpha=1 because in finite
//--- precision arithmetics (x1-x0)+x0 may differ from x1.
//--- So while x1 is correctly bounded (it lie EXACTLY on
//--- boundary,if it is active),(x1-x0)+x0 may be
//--- not bounded.
   v=0.0;
   for(int i_=0;i_<=n-1;i_++)
      v+=state.m_d[i_]*state.m_g[i_];
//--- change values
   state.m_dginit=v;
   state.m_finit=state.m_f;
//--- check
   if(!(ASAD1Norm(state)<=state.m_stpmax || state.m_stpmax==0.0))
     {
      stepfound=false;
      //--- function call, return result
      return(Func_lbl_24(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
     }
//--- Try alpha=1 step first
   for(i=0;i<=n-1;i++)
      state.m_x[i]=CApServ::BoundVal(state.m_xk[i]-state.m_g[i],state.m_bndl[i],state.m_bndu[i]);
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_needfg=true;
   state.m_rstate.stage=2;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_24(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- check
   if(!stepfound)
     {
      //--- alpha=1 is too large,try smaller values
      state.m_stp=1;
      //--- function call
      CLinMin::LinMinNormalized(state.m_d,state.m_stp,n);
      //--- change values
      state.m_dginit=state.m_dginit/state.m_stp;
      state.m_stp=m_gpadecay*state.m_stp;
      //--- check
      if(state.m_stpmax>0.0)
         state.m_stp=MathMin(state.m_stp,state.m_stpmax);
      //--- function call, return result
      return(Func_lbl_27(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
     }
//--- we are at the boundary(ies)
   for(int i_=0;i_<=n-1;i_++)
      state.m_xn[i_]=state.m_x[i_];
   state.m_stp=1;
//--- function call, return result
   return(Func_lbl_26(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_26(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
   state.m_repiterationscount=state.m_repiterationscount+1;
//--- check
   if(!state.m_xrep)
      return(Func_lbl_29(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- progress report
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=4;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_27(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
   v=state.m_stp;
//--- copy
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xk[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_x[i_]+v*state.m_d[i_];
   ClearRequestFields(state);
//--- change values
   state.m_needfg=true;
   state.m_rstate.stage=3;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_29(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- Calculate new set of active constraints.
//--- Reset counter if active set was changed.
//--- Prepare for the new iteration
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(state.m_xn[i]==state.m_bndl[i] || state.m_xn[i]==state.m_bndu[i])
         state.m_an[i]=0;
      else
         state.m_an[i]=1;
     }
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(state.m_ak[i]!=state.m_an[i])
        {
         state.m_acount=-1;
         break;
        }
     }
   state.m_acount=state.m_acount+1;
//--- copy
   for(int i_=0;i_<=n-1;i_++)
      state.m_xk[i_]=state.m_xn[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_ak[i_]=state.m_an[i_];
//--- Stopping conditions
   if(!(state.m_repiterationscount>=state.m_maxits&&state.m_maxits>0))
      return(Func_lbl_31(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- Too many iterations
   state.m_repterminationtype=5;
//--- check
   if(!state.m_xrep)
      return(false);
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=5;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_31(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- check
   if(ASABoundedAntigradNorm(state)>state.m_epsg)
      return(Func_lbl_35(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- Gradient is small enough
   state.m_repterminationtype=4;
//--- check
   if(!state.m_xrep)
      return(false);
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=6;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_35(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- change value
   v=0.0;
   for(int i_=0;i_<=n-1;i_++)
      v+=state.m_d[i_]*state.m_d[i_];
//--- check
   if(MathSqrt(v)*state.m_stp>state.m_epsx)
      return(Func_lbl_39(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- Step size is too small,no further improvement is
//--- possible
   state.m_repterminationtype=2;
//--- check
   if(!state.m_xrep)
      return(false);
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=7;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_39(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- check
   if(state.m_finit-state.m_f>state.m_epsf*MathMax(MathAbs(state.m_finit),MathMax(MathAbs(state.m_f),1.0)))
      return(Func_lbl_43(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- F(k+1)-F(k) is small enough
   state.m_repterminationtype=1;
//--- check
   if(!state.m_xrep)
      return(false);
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=8;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_43(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- Decide - should we switch algorithm or not
   if(ASAUIsEmpty(state))
     {
      //--- check
      if(ASAGINorm(state)>=state.m_mu*ASAD1Norm(state))
        {
         state.m_curalgo=1;
         //--- function call, return result
         return(Func_lbl_19(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
        }
      else
         state.m_mu=state.m_mu*m_asarho;
     }
   else
     {
      //--- check
      if(state.m_acount==m_n1)
        {
         //--- check
         if(ASAGINorm(state)>=state.m_mu*ASAD1Norm(state))
           {
            state.m_curalgo=1;
            //--- function call, return result
            return(Func_lbl_19(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
           }
        }
     }
//--- Next iteration
   state.m_k=state.m_k+1;
//--- function call, return result
   return(Func_lbl_21(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_49(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- Store G[k] for later calculation of Y[k]
   for(i=0;i<=n-1;i++)
      state.m_yk[i]=-state.m_gc[i];
//--- Make a CG step in direction given by DK[]:
//--- * calculate step. Step projection into feasible set
//---   is used. It has several benefits: a) step may be
//---   found with usual line search,b) multiple constraints
//---   may be activated with one step,c) activated constraints
//---   are detected in a natural way - just compare x[i] with
//---   bounds
//--- * update active set,set B to True,if there
//---   were changes in the set.
   for(int i_=0;i_<=n-1;i_++)
      state.m_d[i_]=state.m_dk[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_xn[i_]=state.m_xk[i_];
//--- change values
   state.m_mcstage=0;
   state.m_stp=1;
//--- function call
   CLinMin::LinMinNormalized(state.m_d,state.m_stp,n);
//--- check
   if(state.m_laststep!=0.0)
      state.m_stp=state.m_laststep;
//--- function call
   CLinMin::MCSrch(n,state.m_xn,state.m_f,state.m_gc,state.m_d,state.m_stp,state.m_stpmax,m_gtol,mcinfo,state.m_nfev,state.m_work,state.m_lstate,state.m_mcstage);
//--- function call, return result
   return(Func_lbl_51(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_51(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- check
   if(state.m_mcstage==0)
      return(Func_lbl_52(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- preprocess data: bound State.XN so it belongs to the
//--- feasible set and store it in the State.X
   for(i=0;i<=n-1;i++)
      state.m_x[i]=CApServ::BoundVal(state.m_xn[i],state.m_bndl[i],state.m_bndu[i]);
//--- RComm
   ClearRequestFields(state);
//--- change values
   state.m_needfg=true;
   state.m_rstate.stage=9;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_52(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
   diffcnt=0;
   for(i=0;i<=n-1;i++)
     {
      //--- XN contains unprojected result,project it,
      //--- save copy to X (will be used for progress reporting)
      state.m_xn[i]=CApServ::BoundVal(state.m_xn[i],state.m_bndl[i],state.m_bndu[i]);
      //--- update active set
      if(state.m_xn[i]==state.m_bndl[i] || state.m_xn[i]==state.m_bndu[i])
         state.m_an[i]=0;
      else
         state.m_an[i]=1;
      //--- check
      if(state.m_an[i]!=state.m_ak[i])
         diffcnt=diffcnt+1;
      state.m_ak[i]=state.m_an[i];
     }
   for(int i_=0;i_<=n-1;i_++)
      state.m_xk[i_]=state.m_xn[i_];
//--- change values
   state.m_repnfev=state.m_repnfev+state.m_nfev;
   state.m_repiterationscount=state.m_repiterationscount+1;
//--- check
   if(!state.m_xrep)
      return(Func_lbl_53(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- progress report
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=10;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_53(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- Update info about step length
   v=0.0;
   for(int i_=0;i_<=n-1;i_++)
      v+=state.m_d[i_]*state.m_d[i_];
   state.m_laststep=MathSqrt(v)*state.m_stp;
//--- Check stopping conditions.
   if(ASABoundedAntigradNorm(state)>state.m_epsg)
      return(Func_lbl_55(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- Gradient is small enough
   state.m_repterminationtype=4;
//--- check
   if(!state.m_xrep)
      return(false);
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=11;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_55(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- check
   if(!(state.m_repiterationscount>=state.m_maxits && state.m_maxits>0))
      return(Func_lbl_59(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- Too many iterations
   state.m_repterminationtype=5;
//--- check
   if(!state.m_xrep)
      return(false);
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=12;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_59(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- check
   if(!(ASAGINorm(state)>=state.m_mu*ASAD1Norm(state) && diffcnt==0))
      return(Func_lbl_63(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- These conditions (EpsF/EpsX) are explicitly or implicitly
//--- related to the current step size and influenced
//--- by changes in the active constraints.
//--- For these reasons they are checked only when we don't
//--- want to 'unstick' at the end of the iteration and there
//--- were no changes in the active set.
//--- NOTE: consition |G|>=Mu*|D1| must be exactly opposite
//--- to the condition used to switch back to GPA. At least
//--- one inequality must be strict,otherwise infinite cycle
//--- may occur when |G|=Mu*|D1| (we DON'T test stopping
//--- conditions and we DON'T switch to GPA,so we cycle
//--- indefinitely).
   if(state.m_fold-state.m_f>state.m_epsf*MathMax(MathAbs(state.m_fold),MathMax(MathAbs(state.m_f),1.0)))
      return(Func_lbl_65(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- F(k+1)-F(k) is small enough
   state.m_repterminationtype=1;
   if(!state.m_xrep)
      return(false);
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=13;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_63(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- Check conditions for switching
   if(ASAGINorm(state)<state.m_mu*ASAD1Norm(state))
     {
      state.m_curalgo=0;
      //--- function call, return result
      return(Func_lbl_17(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
     }
//--- check
   if(diffcnt>0)
     {
      //--- check
      if(ASAUIsEmpty(state) || diffcnt>=m_n2)
         state.m_curalgo=1;
      else
         state.m_curalgo=0;
      //--- function call, return result
      return(Func_lbl_17(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
     }
//--- Calculate D(k+1)
//--- Line search may result in:
//--- * maximum feasible step being taken (already processed)
//--- * point satisfying Wolfe conditions
//--- * some kind of error (CG is restarted by assigning 0.0 to Beta)
   if(mcinfo==1)
     {
      //--- Standard Wolfe conditions are satisfied:
      //--- * calculate Y[K] and BetaK
      for(int i_=0;i_<=n-1;i_++)
         state.m_yk[i_]=state.m_yk[i_]+state.m_gc[i_];
      //--- change value
      vv=0.0;
      for(int i_=0;i_<=n-1;i_++)
         vv+=state.m_yk[i_]*state.m_dk[i_];
      //--- change value
      v=0.0;
      for(int i_=0;i_<=n-1;i_++)
         v+=state.m_gc[i_]*state.m_gc[i_];
      state.m_betady=v/vv;
      //--- change value
      v=0.0;
      for(int i_=0;i_<=n-1;i_++)
         v+=state.m_gc[i_]*state.m_yk[i_];
      state.m_betahs=v/vv;
      //--- check
      if(state.m_cgtype==0)
         betak=state.m_betady;
      //--- check
      if(state.m_cgtype==1)
         betak=MathMax(0,MathMin(state.m_betady,state.m_betahs));
     }
   else
     {
      //--- Something is wrong (may be function is too wild or too flat).
      //--- We'll set BetaK=0,which will restart CG algorithm.
      //--- We can stop later (during normal checks) if stopping conditions are met.
      betak=0;
      state.m_debugrestartscount=state.m_debugrestartscount+1;
     }
//--- change values
   for(int i_=0;i_<=n-1;i_++)
      state.m_dn[i_]=-state.m_gc[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_dn[i_]=state.m_dn[i_]+betak*state.m_dk[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_dk[i_]=state.m_dn[i_];
//--- update other information
   state.m_fold=state.m_f;
   state.m_k=state.m_k+1;
//--- function call, return result
   return(Func_lbl_49(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for MinASAIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CMinComp::Func_lbl_65(CMinASAState &state,int &n,int &i,
                                  int &mcinfo,int &diffcnt,bool &b,
                                  bool &stepfound,double &betak,
                                  double &v,double &vv)
  {
//--- check
   if(state.m_laststep>state.m_epsx)
      return(Func_lbl_63(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv));
//--- X(k+1)-X(k) is small enough
   state.m_repterminationtype=2;
//--- check
   if(!state.m_xrep)
      return(false);
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=14;
//--- Saving state
   Func_lbl_rcomm(state,n,i,mcinfo,diffcnt,b,stepfound,betak,v,vv);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
