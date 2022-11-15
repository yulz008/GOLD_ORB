//+------------------------------------------------------------------+
//|                                                diffequations.mqh |
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
#include "alglibinternal.mqh"
//+------------------------------------------------------------------+
//| Auxiliary class for CODESolver                                   |
//+------------------------------------------------------------------+
class CODESolverState
  {
public:
   int               m_n;
   int               m_m;
   double            m_xscale;
   double            m_h;
   double            m_eps;
   bool              m_fraceps;
   int               m_repterminationtype;
   int               m_repnfev;
   int               m_solvertype;
   bool              m_needdy;
   double            m_x;
   RCommState        m_rstate;
   //--- arrays
   double            m_yc[];
   double            m_escale[];
   double            m_xg[];
   double            m_y[];
   double            m_dy[];
   double            m_yn[];
   double            m_yns[];
   double            m_rka[];
   double            m_rkc[];
   double            m_rkcs[];
   //--- matrices
   CMatrixDouble     m_ytbl;
   CMatrixDouble     m_rkb;
   CMatrixDouble     m_rkk;

public:
                     CODESolverState(void);
                    ~CODESolverState(void);

   void              Copy(CODESolverState &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CODESolverState::CODESolverState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CODESolverState::~CODESolverState(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CODESolverState::Copy(CODESolverState &obj)
  {
//--- copy variables
   m_n=obj.m_n;
   m_m=obj.m_m;
   m_xscale=obj.m_xscale;
   m_h=obj.m_h;
   m_eps=obj.m_eps;
   m_fraceps=obj.m_fraceps;
   m_repterminationtype=obj.m_repterminationtype;
   m_repnfev=obj.m_repnfev;
   m_solvertype=obj.m_solvertype;
   m_needdy=obj.m_needdy;
   m_x=obj.m_x;
   m_rstate.Copy(obj.m_rstate);
//--- copy arrays
   ArrayCopy(m_yc,obj.m_yc);
   ArrayCopy(m_escale,obj.m_escale);
   ArrayCopy(m_xg,obj.m_xg);
   ArrayCopy(m_y,obj.m_y);
   ArrayCopy(m_dy,obj.m_dy);
   ArrayCopy(m_yn,obj.m_yn);
   ArrayCopy(m_yns,obj.m_yns);
   ArrayCopy(m_rka,obj.m_rka);
   ArrayCopy(m_rkc,obj.m_rkc);
   ArrayCopy(m_rkcs,obj.m_rkcs);
//--- copy matrices
   m_ytbl=obj.m_ytbl;
   m_rkb=obj.m_rkb;
   m_rkk=obj.m_rkk;
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CODESolverState                  |
//+------------------------------------------------------------------+
class CODESolverStateShell
  {
private:
   CODESolverState   m_innerobj;
public:
   //--- constructors, destructor
                     CODESolverStateShell(void);
                     CODESolverStateShell(CODESolverState &obj);
                    ~CODESolverStateShell(void);
   //--- methods
   bool              GetNeedDY(void);
   void              SetNeedDY(const bool b);
   double            GetX(void);
   void              SetX(const double d);
   CODESolverState *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CODESolverStateShell::CODESolverStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CODESolverStateShell::CODESolverStateShell(CODESolverState &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CODESolverStateShell::~CODESolverStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needdy                         |
//+------------------------------------------------------------------+
bool CODESolverStateShell::GetNeedDY(void)
  {
//--- return result
   return(m_innerobj.m_needdy);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needdy                        |
//+------------------------------------------------------------------+
void CODESolverStateShell::SetNeedDY(const bool b)
  {
//--- change value
   m_innerobj.m_needdy=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable x                              |
//+------------------------------------------------------------------+
double CODESolverStateShell::GetX(void)
  {
//--- return result
   return(m_innerobj.m_x);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable x                             |
//+------------------------------------------------------------------+
void CODESolverStateShell::SetX(const double d)
  {
//--- change value
   m_innerobj.m_x=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CODESolverState *CODESolverStateShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CODESolver                                   |
//+------------------------------------------------------------------+
class CODESolverReport
  {
public:
   //--- class variables
   int               m_nfev;
   int               m_terminationtype;
   //--- constructor, destructor
                     CODESolverReport(void);
                    ~CODESolverReport(void);
   //--- copy
   void              Copy(CODESolverReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CODESolverReport::CODESolverReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CODESolverReport::~CODESolverReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CODESolverReport::Copy(CODESolverReport &obj)
  {
//--- copy variables
   m_nfev=obj.m_nfev;
   m_terminationtype=obj.m_terminationtype;
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CODESolverReport                 |
//+------------------------------------------------------------------+
class CODESolverReportShell
  {
private:
   CODESolverReport  m_innerobj;
public:
   //--- constructor, destructor
                     CODESolverReportShell(void);
                     CODESolverReportShell(CODESolverReport &obj);
                    ~CODESolverReportShell(void);
   //--- methods
   int               GetNFev(void);
   void              SetNFev(const int i);
   int               GetTerminationType(void);
   void              SetTerminationType(const int i);
   CODESolverReport *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CODESolverReportShell::CODESolverReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CODESolverReportShell::CODESolverReportShell(CODESolverReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CODESolverReportShell::~CODESolverReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable nfev                           |
//+------------------------------------------------------------------+
int CODESolverReportShell::GetNFev(void)
  {
//--- return result
   return(m_innerobj.m_nfev);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable nfev                          |
//+------------------------------------------------------------------+
void CODESolverReportShell::SetNFev(const int i)
  {
//--- change value
   m_innerobj.m_nfev=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable terminationtype                |
//+------------------------------------------------------------------+
int CODESolverReportShell::GetTerminationType(void)
  {
//--- return result
   return(m_innerobj.m_terminationtype);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable terminationtype               |
//+------------------------------------------------------------------+
void CODESolverReportShell::SetTerminationType(const int i)
  {
//--- change value
   m_innerobj.m_terminationtype=i;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CODESolverReport *CODESolverReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Solution of ordinary differential equations                      |
//+------------------------------------------------------------------+
class CODESolver
  {
private:
   //--- private method
   static void       ODESolverInit(int solvertype,double &y[],const int n,double &x[],const int m,const double eps,double h,CODESolverState &state);
   //--- auxiliary functions for ODESolverIteration
   static void       Func_lbl_rcomm(CODESolverState &state,int n,int m,int i,int j,int k,int klimit,bool gridpoint,double xc,double v,double h,double h2,double err,double maxgrowpow);
   static bool       Func_lbl_6(CODESolverState &state,int &n,int &m,int &i,int &j,int &k,int &klimit,bool &gridpoint,double &xc,double &v,double &h,double &h2,double &err,double &maxgrowpow);
   static bool       Func_lbl_8(CODESolverState &state,int &n,int &m,int &i,int &j,int &k,int &klimit,bool &gridpoint,double &xc,double &v,double &h,double &h2,double &err,double &maxgrowpow);
   static bool       Func_lbl_10(CODESolverState &state,int &n,int &m,int &i,int &j,int &k,int &klimit,bool &gridpoint,double &xc,double &v,double &h,double &h2,double &err,double &maxgrowpow);
public:
   //--- class constants
   static const double m_odesolvermaxgrow;
   static const double m_odesolvermaxshrink;
   //--- constructor, destructor
                     CODESolver(void);
                    ~CODESolver(void);
   //--- public methods
   static void       ODESolverRKCK(double &y[],const int n,double &x[],const int m,const double eps,const double h,CODESolverState &state);
   static void       ODESolverResults(CODESolverState &state,int &m,double &xtbl[],CMatrixDouble &ytbl,CODESolverReport &rep);
   static bool       ODESolverIteration(CODESolverState &state);
  };
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const double CODESolver::m_odesolvermaxgrow=3.0;
const double CODESolver::m_odesolvermaxshrink=10.0;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CODESolver::CODESolver(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CODESolver::~CODESolver(void)
  {

  }
//+------------------------------------------------------------------+
//| Cash-Karp adaptive ODE solver.                                   |
//| This subroutine solves ODE  Y'=f(Y,x) with initial conditions    |
//| Y(xs)=Ys (here Y may be single variable or vector of N variables)|
//| INPUT PARAMETERS:                                                |
//|     Y       -   initial conditions, array[0..N-1].               |
//|                 contains values of Y[] at X[0]                   |
//|     N       -   system size                                      |
//|     X       -   points at which Y should be tabulated,           |
//|                 array[0..M-1] integrations starts at X[0], ends  |
//|                 at X[M-1], intermediate values at X[i] are       |
//|                 returned too.                                    |
//|                 SHOULD BE ORDERED BY ASCENDING OR BY DESCENDING!!|
//|     M       -   number of intermediate points + first point +    |
//|                 last point:                                      |
//|                 * M>2 means that you need both Y(X[M-1]) and M-2 |
//|                   values at intermediate points                  |
//|                 * M=2 means that you want just to integrate from |
//|                   X[0] to X[1] and don't interested in           |
//|                   intermediate values.                           |
//|                 * M=1 means that you don't want to integrate :)  |
//|                   it is degenerate case, but it will be handled  |
//|                   correctly.                                     |
//|                 * M<1 means error                                |
//|     Eps     -   tolerance (absolute/relative error on each step  |
//|                 will be less than Eps). When passing:            |
//|                 * Eps>0, it means desired ABSOLUTE error         |
//|                 * Eps<0, it means desired RELATIVE error.        |
//|                   Relative errors are calculated with respect to |
//|                   maximum values of Y seen so far. Be careful to |
//|                   use this criterion when starting from Y[] that |
//|                   are close to zero.                             |
//|     H       -   initial step lenth, it will be adjusted          |
//|                 automatically after the first step. If H=0, step |
//|                 will be selected automatically (usualy it will   |
//|                 be equal to 0.001 of min(x[i]-x[j])).            |
//| OUTPUT PARAMETERS                                                |
//|     State   -   structure which stores algorithm state between   |
//|                 subsequent calls of OdeSolverIteration. Used     |
//|                 for reverse communication. This structure should |
//|                 be passed to the OdeSolverIteration subroutine.  |
//| SEE ALSO                                                         |
//|     AutoGKSmoothW, AutoGKSingular, AutoGKIteration, AutoGKResults|
//+------------------------------------------------------------------+
static void CODESolver::ODESolverRKCK(double &y[],const int n,double &x[],
                                      const int m,const double eps,
                                      const double h,CODESolverState &state)
  {
//--- check
   if(!CAp::Assert(n>=1,"ODESolverRKCK: N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,"ODESolverRKCK: M<1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,"ODESolverRKCK: Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=m,"ODESolverRKCK: Length(X)<M!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),"ODESolverRKCK: Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,m),"ODESolverRKCK: Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(eps),"ODESolverRKCK: Eps is not finite!"))
      return;
//--- check
   if(!CAp::Assert(eps!=0.0,"ODESolverRKCK: Eps is zero!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(h),"ODESolverRKCK: H is not finite!"))
      return;
//--- function call
   ODESolverInit(0,y,n,x,m,eps,h,state);
  }
//+------------------------------------------------------------------+
//| ODE solver results                                               |
//| Called after OdeSolverIteration returned False.                  |
//| INPUT PARAMETERS:                                                |
//|     State   -   algorithm state (used by OdeSolverIteration).    |
//| OUTPUT PARAMETERS:                                               |
//|     M       -   number of tabulated values, M>=1                 |
//|     XTbl    -   array[0..M-1], values of X                       |
//|     YTbl    -   array[0..M-1,0..N-1], values of Y in X[i]        |
//|     Rep     -   solver report:                                   |
//|                 * Rep.TerminationType completetion code:         |
//|                     * -2    X is not ordered  by                 |
//|                             ascending/descending or there are    |
//|                             non-distinct X[], i.e.  X[i]=X[i+1]  |
//|                     * -1    incorrect parameters were specified  |
//|                     *  1    task has been solved                 |
//|                 * Rep.NFEV contains number of function           |
//|                   calculations                                   |
//+------------------------------------------------------------------+
static void CODESolver::ODESolverResults(CODESolverState &state,int &m,
                                         double &xtbl[],CMatrixDouble &ytbl,
                                         CODESolverReport &rep)
  {
//--- create variables
   double v=0;
   int    i=0;
   int    i_=0;
//--- initialization
   m=0;
   rep.m_terminationtype=state.m_repterminationtype;
//--- check
   if(rep.m_terminationtype>0)
     {
      //--- change values
      m=state.m_m;
      rep.m_nfev=state.m_repnfev;
      //--- allocation
      ArrayResizeAL(xtbl,state.m_m);
      v=state.m_xscale;
      //--- calculation
      for(i_=0;i_<=state.m_m-1;i_++)
         xtbl[i_]=v*state.m_xg[i_];
      //--- allocation
      ytbl.Resize(state.m_m,state.m_n);
      for(i=0;i<=state.m_m-1;i++)
        {
         for(i_=0;i_<=state.m_n-1;i_++)
            ytbl[i].Set(i_,state.m_ytbl[i][i_]);
        }
     }
   else
      rep.m_nfev=0;
  }
//+------------------------------------------------------------------+
//| Internal initialization subroutine                               |
//+------------------------------------------------------------------+
static void CODESolver::ODESolverInit(int solvertype,double &y[],const int n,
                                      double &x[],const int m,const double eps,
                                      double h,CODESolverState &state)
  {
//--- create variables
   int    i=0;
   double v=0;
   int    i_=0;
//--- Prepare RComm
   ArrayResizeAL(state.m_rstate.ia,6);
   ArrayResizeAL(state.m_rstate.ba,1);
   ArrayResizeAL(state.m_rstate.ra,6);
   state.m_rstate.stage=-1;
   state.m_needdy=false;
//--- check parameters.
   if((n<=0 || m<1) || eps==0.0)
     {
      state.m_repterminationtype=-1;
      return;
     }
//--- check
   if(h<0.0)
      h=-h;
//--- quick exit if necessary.
//--- after this block we assume that M>1
   if(m==1)
     {
      //--- change values
      state.m_repnfev=0;
      state.m_repterminationtype=1;
      state.m_ytbl.Resize(1,n);
      for(i_=0;i_<=n-1;i_++)
         state.m_ytbl[0].Set(i_,y[i_]);
      //--- allocation
      ArrayResizeAL(state.m_xg,m);
      for(i_=0;i_<=m-1;i_++)
         state.m_xg[i_]=x[i_];
      //--- exit the function
      return;
     }
//--- check again: correct order of X[]
   if(x[1]==x[0])
     {
      state.m_repterminationtype=-2;
      return;
     }
   for(i=1;i<=m-1;i++)
     {
      //--- check
      if((x[1]>x[0] && x[i]<=x[i-1]) || (x[1]<x[0] && x[i]>=x[i-1]))
        {
         state.m_repterminationtype=-2;
         return;
        }
     }
//--- auto-select H if necessary
   if(h==0.0)
     {
      v=MathAbs(x[1]-x[0]);
      for(i=2;i<=m-1;i++)
         v=MathMin(v,MathAbs(x[i]-x[i-1]));
      h=0.001*v;
     }
//--- store parameters
   state.m_n=n;
   state.m_m=m;
   state.m_h=h;
   state.m_eps=MathAbs(eps);
   state.m_fraceps=eps<0.0;
//--- allocation
   ArrayResizeAL(state.m_xg,m);
   for(i_=0;i_<=m-1;i_++)
      state.m_xg[i_]=x[i_];
//--- check
   if(x[1]>x[0])
      state.m_xscale=1;
   else
     {
      state.m_xscale=-1;
      for(i_=0;i_<=m-1;i_++)
         state.m_xg[i_]=-1*state.m_xg[i_];
     }
//--- allocation
   ArrayResizeAL(state.m_yc,n);
   for(i_=0;i_<=n-1;i_++)
      state.m_yc[i_]=y[i_];
//--- change values
   state.m_solvertype=solvertype;
   state.m_repterminationtype=0;
//--- Allocate arrays
   ArrayResizeAL(state.m_y,n);
   ArrayResizeAL(state.m_dy,n);
  }
//+------------------------------------------------------------------+
//| Iterative method                                                 |
//+------------------------------------------------------------------+
static bool CODESolver::ODESolverIteration(CODESolverState &state)
  {
//--- create variables
   int    n=0;
   int    m=0;
   int    i=0;
   int    j=0;
   int    k=0;
   double xc=0;
   double v=0;
   double h=0;
   double h2=0;
   bool   gridpoint;
   double err=0;
   double maxgrowpow=0;
   int    klimit=0;
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
      i=state.m_rstate.ia[2];
      j=state.m_rstate.ia[3];
      k=state.m_rstate.ia[4];
      klimit=state.m_rstate.ia[5];
      gridpoint=state.m_rstate.ba[0];
      xc=state.m_rstate.ra[0];
      v=state.m_rstate.ra[1];
      h=state.m_rstate.ra[2];
      h2=state.m_rstate.ra[3];
      err=state.m_rstate.ra[4];
      maxgrowpow=state.m_rstate.ra[5];
     }
   else
     {
      //--- initialization
      n=-983;
      m=-989;
      i=-834;
      j=900;
      k=-287;
      klimit=364;
      gridpoint=false;
      xc=-338;
      v=-686;
      h=912;
      h2=585;
      err=497;
      maxgrowpow=-271;
     }
//--- check
   if(state.m_rstate.stage==0)
     {
      //--- change values
      state.m_needdy=false;
      state.m_repnfev=state.m_repnfev+1;
      v=h*state.m_xscale;
      for(i_=0;i_<=n-1;i_++)
         state.m_rkk[k].Set(i_,v*state.m_dy[i_]);
      //--- update YN/YNS
      v=state.m_rkc[k];
      for(i_=0;i_<=n-1;i_++)
         state.m_yn[i_]=state.m_yn[i_]+v*state.m_rkk[k][i_];
      v=state.m_rkcs[k];
      for(i_=0;i_<=n-1;i_++)
         state.m_yns[i_]=state.m_yns[i_]+v*state.m_rkk[k][i_];
      k=k+1;
      return(Func_lbl_8(state,n,m,i,j,k,klimit,gridpoint,xc,v,h,h2,err,maxgrowpow));
     }
//--- Routine body
//--- prepare
   if(state.m_repterminationtype!=0)
      return(false);
//--- change values
   n=state.m_n;
   m=state.m_m;
   h=state.m_h;
   maxgrowpow=MathPow(m_odesolvermaxgrow,5);
   state.m_repnfev=0;
//--- some preliminary checks for internal errors
//--- after this we assume that H>0 and M>1
   if(!CAp::Assert(state.m_h>0.0,"ODESolver: internal error"))
      return(false);
//--- check
   if(!CAp::Assert(m>1,"ODESolverIteration: internal error"))
      return(false);
//--- choose solver
   if(state.m_solvertype!=0)
      return(false);
//--- Cask-Karp solver
//--- Prepare coefficients table.
//--- Check it for errors
   ArrayResizeAL(state.m_rka,6);
//--- calculation
   state.m_rka[0]=0;
   state.m_rka[1]=1.0/5.0;
   state.m_rka[2]=3.0/10.0;
   state.m_rka[3]=3.0/5.0;
   state.m_rka[4]=1;
   state.m_rka[5]=7.0/8.0;
   state.m_rkb.Resize(6,5);
   state.m_rkb[1].Set(0,1.0/5.0);
   state.m_rkb[2].Set(0,3.0/40.0);
   state.m_rkb[2].Set(1,9.0/40.0);
   state.m_rkb[3].Set(0,3.0/10.0);
   state.m_rkb[3].Set(1,-(9.0/10.0));
   state.m_rkb[3].Set(2,6.0/5.0);
   state.m_rkb[4].Set(0,-(11.0/54.0));
   state.m_rkb[4].Set(1,5.0/2.0);
   state.m_rkb[4].Set(2,-(70.0/27.0));
   state.m_rkb[4].Set(3,35.0/27.0);
   state.m_rkb[5].Set(0,1631.0/55296.0);
   state.m_rkb[5].Set(1,175.0/512.0);
   state.m_rkb[5].Set(2,575.0/13824.0);
   state.m_rkb[5].Set(3,44275.0/110592.0);
   state.m_rkb[5].Set(4,253.0/4096.0);
//--- allocation
   ArrayResizeAL(state.m_rkc,6);
//--- calculation
   state.m_rkc[0]=37.0/378.0;
   state.m_rkc[1]=0;
   state.m_rkc[2]=250.0/621.0;
   state.m_rkc[3]=125.0/594.0;
   state.m_rkc[4]=0;
   state.m_rkc[5]=512.0/1771.0;
//--- allocation
   ArrayResizeAL(state.m_rkcs,6);
//--- calculation
   state.m_rkcs[0]=2825.0/27648.0;
   state.m_rkcs[1]=0;
   state.m_rkcs[2]=18575.0/48384.0;
   state.m_rkcs[3]=13525.0/55296.0;
   state.m_rkcs[4]=277.0/14336.0;
   state.m_rkcs[5]=1.0/4.0;
   state.m_rkk.Resize(6,n);
//--- Main cycle consists of two iterations:
//--- * outer where we travel from X[i-1] to X[i]
//--- * inner where we travel inside [X[i-1],X[i]]
   state.m_ytbl.Resize(m,n);
   ArrayResizeAL(state.m_escale,n);
   ArrayResizeAL(state.m_yn,n);
   ArrayResizeAL(state.m_yns,n);
//--- change value
   xc=state.m_xg[0];
   for(i_=0;i_<=n-1;i_++)
      state.m_ytbl[0].Set(i_,state.m_yc[i_]);
   for(j=0;j<=n-1;j++)
      state.m_escale[j]=0;
   i=1;
//--- check
   if(i>m-1)
     {
      state.m_repterminationtype=1;
      //--- return result
      return(false);
     }
//--- begin inner iteration
   return(Func_lbl_6(state,n,m,i,j,k,klimit,gridpoint,xc,v,h,h2,err,maxgrowpow));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for ODESolverIteration. Is a product to get   |
//| rid of the operator unconditional jump goto                      |
//+------------------------------------------------------------------+
static void CODESolver::Func_lbl_rcomm(CODESolverState &state,int n,int m,
                                       int i,int j,int k,int klimit,
                                       bool gridpoint,double xc,double v,
                                       double h,double h2,double err,
                                       double maxgrowpow)
  {
//--- save
   state.m_rstate.ia[0]=n;
   state.m_rstate.ia[1]=m;
   state.m_rstate.ia[2]=i;
   state.m_rstate.ia[3]=j;
   state.m_rstate.ia[4]=k;
   state.m_rstate.ia[5]=klimit;
   state.m_rstate.ba[0]=gridpoint;
   state.m_rstate.ra[0]=xc;
   state.m_rstate.ra[1]=v;
   state.m_rstate.ra[2]=h;
   state.m_rstate.ra[3]=h2;
   state.m_rstate.ra[4]=err;
   state.m_rstate.ra[5]=maxgrowpow;
  }
//+------------------------------------------------------------------+
//| Auxiliary function for ODESolverIteration. Is a product to get   |
//| rid of the operator unconditional jump goto                      |
//+------------------------------------------------------------------+
static bool CODESolver::Func_lbl_6(CODESolverState &state,int &n,int &m,
                                   int &i,int &j,int &k,int &klimit,
                                   bool &gridpoint,double &xc,double &v,
                                   double &h,double &h2,double &err,
                                   double &maxgrowpow)
  {
//--- truncate step if needed (beyond right boundary).
//--- determine should we store X or not
   if(xc+h>=state.m_xg[i])
     {
      h=state.m_xg[i]-xc;
      gridpoint=true;
     }
   else
      gridpoint=false;
//--- Update error scale maximums
//--- These maximums are initialized by zeros,
//--- then updated every iterations.
   for(j=0;j<=n-1;j++)
      state.m_escale[j]=MathMax(state.m_escale[j],MathAbs(state.m_yc[j]));
//--- make one step:
//--- 1. calculate all info needed to do step
//--- 2. update errors scale maximums using values/derivatives
//---    obtained during (1)
//--- Take into account that we use scaling of X to reduce task
//--- to the form where x[0] < x[1] < ... < x[n-1]. So X is
//--- replaced by x=xscale*t,and dy/dx=f(y,x) is replaced
//--- by dy/dt=xscale*f(y,xscale*t).
   for(int i_=0;i_<=n-1;i_++)
      state.m_yn[i_]=state.m_yc[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_yns[i_]=state.m_yc[i_];
   k=0;
//--- function call, return result
   return(Func_lbl_8(state,n,m,i,j,k,klimit,gridpoint,xc,v,h,h2,err,maxgrowpow));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for ODESolverIteration. Is a product to get   |
//| rid of the operator unconditional jump goto                      |
//+------------------------------------------------------------------+
static bool CODESolver::Func_lbl_8(CODESolverState &state,int &n,int &m,
                                   int &i,int &j,int &k,int &klimit,
                                   bool &gridpoint,double &xc,double &v,
                                   double &h,double &h2,double &err,
                                   double &maxgrowpow)
  {
//--- check
   if(k>5)
      return(Func_lbl_10(state,n,m,i,j,k,klimit,gridpoint,xc,v,h,h2,err,maxgrowpow));
//--- prepare data for the next update of YN/YNS
   state.m_x=state.m_xscale*(xc+state.m_rka[k]*h);
   for(int i_=0;i_<=n-1;i_++)
      state.m_y[i_]=state.m_yc[i_];
//--- calculation
   for(j=0;j<=k-1;j++)
     {
      v=state.m_rkb[k][j];
      for(int i_=0;i_<=n-1;i_++)
         state.m_y[i_]=state.m_y[i_]+v*state.m_rkk[j][i_];
     }
   state.m_needdy=true;
   state.m_rstate.stage=0;
//--- Saving state
   Func_lbl_rcomm(state,n,m,i,j,k,klimit,gridpoint,xc,v,h,h2,err,maxgrowpow);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for ODESolverIteration. Is a product to get   |
//| rid of the operator unconditional jump goto                      |
//+------------------------------------------------------------------+
static bool CODESolver::Func_lbl_10(CODESolverState &state,int &n,int &m,
                                    int &i,int &j,int &k,int &klimit,
                                    bool &gridpoint,double &xc,double &v,
                                    double &h,double &h2,double &err,
                                    double &maxgrowpow)
  {
//--- estimate error
   err=0;
   for(j=0;j<=n-1;j++)
     {
      //--- check
      if(!state.m_fraceps)
        {
         //--- absolute error is estimated
         err=MathMax(err,MathAbs(state.m_yn[j]-state.m_yns[j]));
        }
      else
        {
         //--- Relative error is estimated
         v=state.m_escale[j];
         //--- check
         if(v==0.0)
            v=1;
         err=MathMax(err,MathAbs(state.m_yn[j]-state.m_yns[j])/v);
        }
     }
//--- calculate new step,restart if necessary
   if(maxgrowpow*err<=state.m_eps)
      h2=m_odesolvermaxgrow*h;
   else
      h2=h*MathPow(state.m_eps/err,0.2);
//--- check
   if(h2<h/m_odesolvermaxshrink)
      h2=h/m_odesolvermaxshrink;
//--- check
   if(err>state.m_eps)
     {
      h=h2;
      //--- begin inner iteration
      return(Func_lbl_6(state,n,m,i,j,k,klimit,gridpoint,xc,v,h,h2,err,maxgrowpow));
     }
//--- advance position
   xc=xc+h;
   for(int i_=0;i_<=n-1;i_++)
      state.m_yc[i_]=state.m_yn[i_];
//--- update H
   h=h2;
//--- break on grid point
   if(gridpoint)
     {
      //--- save result
      for(int i_=0;i_<=n-1;i_++)
         state.m_ytbl[i].Set(i_,state.m_yc[i_]);
      i=i+1;
      //--- check
      if(i>m-1)
        {
         state.m_repterminationtype=1;
         return(false);
        }
      //--- begin inner iteration
      return(Func_lbl_6(state,n,m,i,j,k,klimit,gridpoint,xc,v,h,h2,err,maxgrowpow));
     }
//--- begin inner iteration
   return(Func_lbl_6(state,n,m,i,j,k,klimit,gridpoint,xc,v,h,h2,err,maxgrowpow));
  }
//+------------------------------------------------------------------+
