//+------------------------------------------------------------------+
//|                                                      solvers.mqh |
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
//| Auxiliary class for CDenseSolver                                 |
//+------------------------------------------------------------------+
class CDenseSolverReport
  {
public:
   double            m_r1;
   double            m_rinf;

                     CDenseSolverReport(void);
                    ~CDenseSolverReport(void);

   void              Copy(CDenseSolverReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDenseSolverReport::CDenseSolverReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CDenseSolverReport::~CDenseSolverReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CDenseSolverReport::Copy(CDenseSolverReport &obj)
  {
//--- copy variables
   m_r1=obj.m_r1;
   m_rinf=obj.m_rinf;
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CDenseSolverReport               |
//+------------------------------------------------------------------+
class CDenseSolverReportShell
  {
private:
   CDenseSolverReport m_innerobj;
public:
   //--- constructors, destructor
                     CDenseSolverReportShell(void);
                     CDenseSolverReportShell(CDenseSolverReport &obj);
                    ~CDenseSolverReportShell(void);
   //--- methods
   double            GetR1(void);
   void              SetR1(const double d);
   double            GetRInf(void);
   void              SetRInf(const double d);
   CDenseSolverReport *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDenseSolverReportShell::CDenseSolverReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CDenseSolverReportShell::CDenseSolverReportShell(CDenseSolverReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CDenseSolverReportShell::~CDenseSolverReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable r1                             |
//+------------------------------------------------------------------+
double CDenseSolverReportShell::GetR1(void)
  {
//--- return result
   return(m_innerobj.m_r1);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable r1                            |
//+------------------------------------------------------------------+
void CDenseSolverReportShell::SetR1(const double d)
  {
//--- change value
   m_innerobj.m_r1=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable rinf                           |
//+------------------------------------------------------------------+
double CDenseSolverReportShell::GetRInf(void)
  {
//--- return result
   return(m_innerobj.m_rinf);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable rinf                          |
//+------------------------------------------------------------------+
void CDenseSolverReportShell::SetRInf(const double d)
  {
//--- change value
   m_innerobj.m_rinf=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CDenseSolverReport *CDenseSolverReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CDenseSolver                                 |
//+------------------------------------------------------------------+
class CDenseSolverLSReport
  {
public:
   double            m_r2;
   CMatrixDouble     m_cx;
   int               m_n;
   int               m_k;
   //--- constructor, destructor
                     CDenseSolverLSReport(void);
                    ~CDenseSolverLSReport(void);
   //--- copy
   void              Copy(CDenseSolverLSReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDenseSolverLSReport::CDenseSolverLSReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CDenseSolverLSReport::~CDenseSolverLSReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CDenseSolverLSReport::Copy(CDenseSolverLSReport &obj)
  {
//--- copy variables
   m_r2=obj.m_r2;
   m_n=obj.m_n;
   m_k=obj.m_k;
//--- copy matrix
   m_cx=obj.m_cx;
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CDenseSolverLSReport             |
//+------------------------------------------------------------------+
class CDenseSolverLSReportShell
  {
private:
   CDenseSolverLSReport m_innerobj;
public:
   //--- constructors, destructor
                     CDenseSolverLSReportShell(void);
                     CDenseSolverLSReportShell(CDenseSolverLSReport &obj);
                    ~CDenseSolverLSReportShell(void);
   //--- methods
   double            GetR2(void);
   void              SetR2(const double d);
   int               GetN(void);
   void              SetN(const int i);
   int               GetK(void);
   void              SetK(const int i);
   CDenseSolverLSReport *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDenseSolverLSReportShell::CDenseSolverLSReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CDenseSolverLSReportShell::CDenseSolverLSReportShell(CDenseSolverLSReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CDenseSolverLSReportShell::~CDenseSolverLSReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable r2                             |
//+------------------------------------------------------------------+
double CDenseSolverLSReportShell::GetR2(void)
  {
//--- return result
   return(m_innerobj.m_r2);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable r2                            |
//+------------------------------------------------------------------+
void CDenseSolverLSReportShell::SetR2(const double d)
  {
//--- change value
   m_innerobj.m_r2=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable n                              |
//+------------------------------------------------------------------+
int CDenseSolverLSReportShell::GetN(void)
  {
//--- return result
   return(m_innerobj.m_n);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable n                             |
//+------------------------------------------------------------------+
void CDenseSolverLSReportShell::SetN(const int i)
  {
//--- change value
   m_innerobj.m_n=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable k                              |
//+------------------------------------------------------------------+
int CDenseSolverLSReportShell::GetK(void)
  {
//--- return result
   return(m_innerobj.m_k);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable k                             |
//+------------------------------------------------------------------+
void CDenseSolverLSReportShell::SetK(const int i)
  {
//--- change value
   m_innerobj.m_k=i;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CDenseSolverLSReport *CDenseSolverLSReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Dense solver                                                     |
//+------------------------------------------------------------------+
class CDenseSolver
  {
private:
   //--- private methods
   static void       RMatrixLUSolveInternal(CMatrixDouble &lua,int &p[],const double scalea,const int n,CMatrixDouble &a,const bool havea,CMatrixDouble &b,const int m,int &info,CDenseSolverReport &rep,CMatrixDouble &x);
   static void       SPDMatrixCholeskySolveInternal(CMatrixDouble &cha,const double sqrtscalea,const int n,const bool isupper,CMatrixDouble &a,const bool havea,CMatrixDouble &b,const int m,int &info,CDenseSolverReport &rep,CMatrixDouble &x);
   static void       CMatrixLUSolveInternal(CMatrixComplex &lua,int &p[],const double scalea,const int n,CMatrixComplex &a,const bool havea,CMatrixComplex &b,const int m,int &info,CDenseSolverReport &rep,CMatrixComplex &x);
   static void       HPDMatrixCholeskySolveInternal(CMatrixComplex &cha,const double sqrtscalea,const int n,const bool isupper,CMatrixComplex &a,const bool havea,CMatrixComplex &b,const int m,int &info,CDenseSolverReport &rep,CMatrixComplex &x);
   static int        CDenseSolverRFSMax(const int n,const double r1,const double rinf);
   static int        CDenseSolverRFSMaxV2(const int n,const double r2);
   static void       RBasicLUSolve(CMatrixDouble &lua,int &p[],const double scalea,const int n,double &xb[],double &tmp[]);
   static void       SPDBasicCholeskySolve(CMatrixDouble &cha,const double sqrtscalea,const int n,const bool isupper,double &xb[],double &tmp[]);
   static void       CBasicLUSolve(CMatrixComplex &lua,int &p[],const double scalea,const int n,al_complex &xb[],al_complex &tmp[]);
   static void       HPDBasicCholeskySolve(CMatrixComplex &cha,const double sqrtscalea,const int n,const bool isupper,al_complex &xb[],al_complex &tmp[]);
public:
   //--- constructor, destructor
                     CDenseSolver(void);
                    ~CDenseSolver(void);
   //--- public methods
   static void       RMatrixSolve(CMatrixDouble &a,const int n,double &b[],int &info,CDenseSolverReport &rep,double &x[]);
   static void       RMatrixSolveM(CMatrixDouble &a,const int n,CMatrixDouble &b,const int m,const bool rfs,int &info,CDenseSolverReport &rep,CMatrixDouble &x);
   static void       RMatrixLUSolve(CMatrixDouble &lua,int &p[],const int n,double &b[],int &info,CDenseSolverReport &rep,double &x[]);
   static void       RMatrixLUSolveM(CMatrixDouble &lua,int &p[],const int n,CMatrixDouble &b,const int m,int &info,CDenseSolverReport &rep,CMatrixDouble &x);
   static void       RMatrixMixedSolve(CMatrixDouble &a,CMatrixDouble &lua,int &p[],const int n,double &b[],int &info,CDenseSolverReport &rep,double &x[]);
   static void       RMatrixMixedSolveM(CMatrixDouble &a,CMatrixDouble &lua,int &p[],const int n,CMatrixDouble &b,const int m,int &info,CDenseSolverReport &rep,CMatrixDouble &x);
   static void       CMatrixSolveM(CMatrixComplex &a,const int n,CMatrixComplex &b,const int m,const bool rfs,int &info,CDenseSolverReport &rep,CMatrixComplex &x);
   static void       CMatrixSolve(CMatrixComplex &a,const int n,al_complex &b[],int &info,CDenseSolverReport &rep,al_complex &x[]);
   static void       CMatrixLUSolveM(CMatrixComplex &lua,int &p[],const int n,CMatrixComplex &b,const int m,int &info,CDenseSolverReport &rep,CMatrixComplex &x);
   static void       CMatrixLUSolve(CMatrixComplex &lua,int &p[],const int n,al_complex &b[],int &info,CDenseSolverReport &rep,al_complex &x[]);
   static void       CMatrixMixedSolveM(CMatrixComplex &a,CMatrixComplex &lua,int &p[],const int n,CMatrixComplex &b,const int m,int &info,CDenseSolverReport &rep,CMatrixComplex &x);
   static void       CMatrixMixedSolve(CMatrixComplex &a,CMatrixComplex &lua,int &p[],const int n,al_complex &b[],int &info,CDenseSolverReport &rep,al_complex &x[]);
   static void       SPDMatrixSolveM(CMatrixDouble &a,const int n,const bool isupper,CMatrixDouble &b,const int m,int &info,CDenseSolverReport &rep,CMatrixDouble &x);
   static void       SPDMatrixSolve(CMatrixDouble &a,const int n,const bool isupper,double &b[],int &info,CDenseSolverReport &rep,double &x[]);
   static void       SPDMatrixCholeskySolveM(CMatrixDouble &cha,const int n,const bool isupper,CMatrixDouble &b,const int m,int &info,CDenseSolverReport &rep,CMatrixDouble &x);
   static void       SPDMatrixCholeskySolve(CMatrixDouble &cha,const int n,const bool isupper,double &b[],int &info,CDenseSolverReport &rep,double &x[]);
   static void       HPDMatrixSolveM(CMatrixComplex &a,const int n,const bool isupper,CMatrixComplex &b,const int m,int &info,CDenseSolverReport &rep,CMatrixComplex &x);
   static void       HPDMatrixSolve(CMatrixComplex &a,const int n,const bool isupper,al_complex &b[],int &info,CDenseSolverReport &rep,al_complex &x[]);
   static void       HPDMatrixCholeskySolveM(CMatrixComplex &cha,const int n,const bool isupper,CMatrixComplex &b,const int m,int &info,CDenseSolverReport &rep,CMatrixComplex &x);
   static void       HPDMatrixCholeskySolve(CMatrixComplex &cha,const int n,const bool isupper,al_complex &b[],int &info,CDenseSolverReport &rep,al_complex &x[]);
   static void       RMatrixSolveLS(CMatrixDouble &a,const int nrows,const int ncols,double &b[],double threshold,int &info,CDenseSolverLSReport &rep,double &x[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDenseSolver::CDenseSolver(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CDenseSolver::~CDenseSolver(void)
  {

  }
//+------------------------------------------------------------------+
//| Dense solver.                                                    |
//| This subroutine solves a system A*x=b, where A is NxN            |
//| non-denegerate real matrix, x and b are vectors.                 |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * condition number estimation                                    |
//| * iterative refinement                                           |
//| * O(N^3) complexity                                              |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1,0..N-1], system matrix              |
//|     N       -   size of A                                        |
//|     B       -   array[0..N-1], right part                        |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   return code:                                     |
//|                 * -3    A is singular, or VERY close to singular.|
//|                         X is filled by zeros in such cases.      |
//|                 * -1    N<=0 was passed                          |
//|                 *  1    task is solved (but matrix A may be      |
//|                         ill-conditioned, check R1/RInf parameters|
//|                         for condition numbers).                  |
//|     Rep     -   solver report, see below for more info           |
//|     X       -   array[0..N-1], it contains:                      |
//|                 * solution of A*x=b if A is non-singular         |
//|                   (well-conditioned or ill-conditioned, but not  |
//|                   very close to singular)                        |
//|                 * zeros, if A is singular or VERY close to       |
//|                   singular (in this case Info=-3).               |
//| SOLVER REPORT                                                    |
//| Subroutine sets following fields of the Rep structure:           |
//| * R1        reciprocal of condition number: 1/cond(A), 1-norm.   |
//| * RInf      reciprocal of condition number: 1/cond(A), inf-norm. |
//+------------------------------------------------------------------+
static void CDenseSolver::RMatrixSolve(CMatrixDouble &a,const int n,double &b[],
                                       int &info,CDenseSolverReport &rep,
                                       double &x[])
  {
//--- create a variable
   int i_=0;
//--- create matrix
   CMatrixDouble bm;
   CMatrixDouble xm;
//--- initialization
   info=0;
//--- check
   if(n<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   bm.Resize(n,1);
//--- filling
   for(i_=0;i_<=n-1;i_++)
      bm[i_].Set(0,b[i_]);
//--- function call
   RMatrixSolveM(a,n,bm,1,true,info,rep,xm);
//--- allocation
   ArrayResizeAL(x,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      x[i_]=xm[i_][0];
  }
//+------------------------------------------------------------------+
//| Dense solver.                                                    |
//| Similar to RMatrixSolve() but solves task with multiple right    |
//| parts (where b and x are NxM matrices).                          |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * condition number estimation                                    |
//| * optional iterative refinement                                  |
//| * O(N^3+M*N^2) complexity                                        |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1,0..N-1], system matrix              |
//|     N       -   size of A                                        |
//|     B       -   array[0..N-1,0..M-1], right part                 |
//|     M       -   right part size                                  |
//|     RFS     -   iterative refinement switch:                     |
//|                 * True - refinement is used.                     |
//|                   Less performance, more precision.              |
//|                 * False - refinement is not used.                |
//|                   More performance, less precision.              |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve                          |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::RMatrixSolveM(CMatrixDouble &a,const int n,
                                        CMatrixDouble &b,const int m,
                                        const bool rfs,int &info,
                                        CDenseSolverReport &rep,
                                        CMatrixDouble &x)
  {
//--- create variables
   double scalea=0;
   int i=0;
   int j=0;
   int i_=0;
//--- create matrix
   CMatrixDouble da;
   CMatrixDouble emptya;
//--- create array
   int p[];
//--- initialization
   info=0;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   da.Resize(n,n);
//--- 1. scale matrix,max(|A[i][j]|)
//--- 2. factorize scaled matrix
//--- 3. solve
   scalea=0;
   for(i=0;i<=n-1;i++)
     {
      for(j=0;j<=n-1;j++)
         scalea=MathMax(scalea,MathAbs(a[i][j]));
     }
//--- check
   if(scalea==0.0)
      scalea=1;
//--- change values
   scalea=1/scalea;
   for(i=0;i<=n-1;i++)
     {
      for(i_=0;i_<=n-1;i_++)
         da[i].Set(i_,a[i][i_]);
     }
//--- function call
   CTrFac::RMatrixLU(da,n,n,p);
//--- check
   if(rfs)
      RMatrixLUSolveInternal(da,p,scalea,n,a,true,b,m,info,rep,x);
   else
      RMatrixLUSolveInternal(da,p,scalea,n,emptya,false,b,m,info,rep,x);
  }
//+------------------------------------------------------------------+
//| Dense solver.                                                    |
//| This subroutine solves a system A*X=B, where A is NxN            |
//| non-denegerate real matrix given by its LU decomposition, X and  |
//| B are NxM real matrices.                                         |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * O(N^2) complexity                                              |
//| * condition number estimation                                    |
//| No iterative refinement is provided because exact form of        |
//| original matrix is not known to subroutine. Use RMatrixSolve or  |
//| RMatrixMixedSolve if you need iterative refinement.              |
//| INPUT PARAMETERS                                                 |
//|     LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU|
//|                 result                                           |
//|     P       -   array[0..N-1], pivots array, RMatrixLU result    |
//|     N       -   size of A                                        |
//|     B       -   array[0..N-1], right part                        |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve                          |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::RMatrixLUSolve(CMatrixDouble &lua,int &p[],
                                         const int n,double &b[],
                                         int &info,CDenseSolverReport &rep,
                                         double &x[])
  {
//--- create matrix
   CMatrixDouble bm;
   CMatrixDouble xm;
//--- create a variable
   int i_=0;
//--- initialization
   info=0;
//--- check
   if(n<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   bm.Resize(n,1);
//--- filling
   for(i_=0;i_<=n-1;i_++)
      bm[i_].Set(0,b[i_]);
//--- function call
   RMatrixLUSolveM(lua,p,n,bm,1,info,rep,xm);
//--- allocation
   ArrayResizeAL(x,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      x[i_]=xm[i_][0];
  }
//+------------------------------------------------------------------+
//| Dense solver.                                                    |
//| Similar to RMatrixLUSolve() but solves task with multiple right  |
//| parts (where b and x are NxM matrices).                          |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * O(M*N^2) complexity                                            |
//| * condition number estimation                                    |
//| No iterative refinement is provided because exact form of        |
//| original matrix is not known to subroutine. Use RMatrixSolve or  |
//| RMatrixMixedSolve if you need iterative refinement.              |
//| INPUT PARAMETERS                                                 |
//|     LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU|
//|                 result                                           |
//|     P       -   array[0..N-1], pivots array, RMatrixLU result    |
//|     N       -   size of A                                        |
//|     B       -   array[0..N-1,0..M-1], right part                 |
//|     M       -   right part size                                  |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve                          |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::RMatrixLUSolveM(CMatrixDouble &lua,int &p[],
                                          const int n,CMatrixDouble &b,
                                          const int m,int &info,
                                          CDenseSolverReport &rep,
                                          CMatrixDouble &x)
  {
//--- create matrix
   CMatrixDouble emptya;
//--- create variables
   int    i=0;
   int    j=0;
   double scalea=0;
//--- initialization
   info=0;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
//--- 1. scale matrix,max(|U[i][j]|)
//---    we assume that LU is in its normal form,i.e. |L[i][j]|<=1
//--- 2. solve
   scalea=0;
   for(i=0;i<=n-1;i++)
     {
      for(j=i;j<=n-1;j++)
         scalea=MathMax(scalea,MathAbs(lua[i][j]));
     }
//--- check
   if(scalea==0.0)
      scalea=1;
//--- change values
   scalea=1/scalea;
//--- function call
   RMatrixLUSolveInternal(lua,p,scalea,n,emptya,false,b,m,info,rep,x);
  }
//+------------------------------------------------------------------+
//| Dense solver.                                                    |
//| This subroutine solves a system A*x=b, where BOTH ORIGINAL A AND |
//| ITS LU DECOMPOSITION ARE KNOWN. You can use it if for some       |
//| reasons you have both A and its LU decomposition.                |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * condition number estimation                                    |
//| * iterative refinement                                           |
//| * O(N^2) complexity                                              |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1,0..N-1], system matrix              |
//|     LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU|
//|                 result                                           |
//|     P       -   array[0..N-1], pivots array, RMatrixLU result    |
//|     N       -   size of A                                        |
//|     B       -   array[0..N-1], right part                        |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolveM                         |
//|     Rep     -   same as in RMatrixSolveM                         |
//|     X       -   same as in RMatrixSolveM                         |
//+------------------------------------------------------------------+
static void CDenseSolver::RMatrixMixedSolve(CMatrixDouble &a,CMatrixDouble &lua,
                                            int &p[],const int n,double &b[],
                                            int &info,CDenseSolverReport &rep,
                                            double &x[])
  {
//--- create a variable
   int i_=0;
//--- create matrix
   CMatrixDouble bm;
   CMatrixDouble xm;
//--- initialization
   info=0;
//--- check
   if(n<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   bm.Resize(n,1);
//--- filling
   for(i_=0;i_<=n-1;i_++)
      bm[i_].Set(0,b[i_]);
//--- function call
   RMatrixMixedSolveM(a,lua,p,n,bm,1,info,rep,xm);
//--- allocation
   ArrayResizeAL(x,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      x[i_]=xm[i_][0];
  }
//+------------------------------------------------------------------+
//| Dense solver.                                                    |
//| Similar to RMatrixMixedSolve() but solves task with multiple     |
//| right parts (where b and x are NxM matrices).                    |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * condition number estimation                                    |
//| * iterative refinement                                           |
//| * O(M*N^2) complexity                                            |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1,0..N-1], system matrix              |
//|     LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU|
//|                 result                                           |
//|     P       -   array[0..N-1], pivots array, RMatrixLU result    |
//|     N       -   size of A                                        |
//|     B       -   array[0..N-1,0..M-1], right part                 |
//|     M       -   right part size                                  |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolveM                         |
//|     Rep     -   same as in RMatrixSolveM                         |
//|     X       -   same as in RMatrixSolveM                         |
//+------------------------------------------------------------------+
static void CDenseSolver::RMatrixMixedSolveM(CMatrixDouble &a,CMatrixDouble &lua,
                                             int &p[],const int n,CMatrixDouble &b,
                                             const int m,int &info,
                                             CDenseSolverReport &rep,
                                             CMatrixDouble &x)
  {
//--- create variables
   double scalea=0;
   int    i=0;
   int    j=0;
//--- initialization
   info=0;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
//--- 1. scale matrix,max(|A[i][j]|)
//--- 2. factorize scaled matrix
//--- 3. solve
   scalea=0;
   for(i=0;i<=n-1;i++)
     {
      for(j=0;j<=n-1;j++)
         scalea=MathMax(scalea,MathAbs(a[i][j]));
     }
//--- check
   if(scalea==0.0)
      scalea=1;
//--- change values
   scalea=1/scalea;
//--- function call
   RMatrixLUSolveInternal(lua,p,scalea,n,a,true,b,m,info,rep,x);
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixSolveM(), but for complex matrices. |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * condition number estimation                                    |
//| * iterative refinement                                           |
//| * O(N^3+M*N^2) complexity                                        |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1,0..N-1], system matrix              |
//|     N       -   size of A                                        |
//|     B       -   array[0..N-1,0..M-1], right part                 |
//|     M       -   right part size                                  |
//|     RFS     -   iterative refinement switch:                     |
//|                 * True - refinement is used.                     |
//|                   Less performance, more precision.              |
//|                 * False - refinement is not used.                |
//|                   More performance, less precision.              |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve                          |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::CMatrixSolveM(CMatrixComplex &a,const int n,
                                        CMatrixComplex &b,const int m,
                                        const bool rfs,int &info,
                                        CDenseSolverReport &rep,
                                        CMatrixComplex &x)
  {
//--- create variables
   double scalea=0;
   int    i=0;
   int    j=0;
   int    i_=0;
//--- create array
   int p[];
//--- create matrix
   CMatrixComplex da;
   CMatrixComplex emptya;
//--- initialization
   info=0;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   da.Resize(n,n);
//--- 1. scale matrix,max(|A[i][j]|)
//--- 2. factorize scaled matrix
//--- 3. solve
   scalea=0;
   for(i=0;i<=n-1;i++)
     {
      for(j=0;j<=n-1;j++)
         scalea=MathMax(scalea,CMath::AbsComplex(a[i][j]));
     }
//--- check
   if(scalea==0.0)
      scalea=1;
//--- change values
   scalea=1/scalea;
   for(i=0;i<=n-1;i++)
     {
      for(i_=0;i_<=n-1;i_++)
         da[i].Set(i_,a[i][i_]);
     }
//--- function call
   CTrFac::CMatrixLU(da,n,n,p);
//--- check
   if(rfs)
      CMatrixLUSolveInternal(da,p,scalea,n,a,true,b,m,info,rep,x);
   else
      CMatrixLUSolveInternal(da,p,scalea,n,emptya,false,b,m,info,rep,x);
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixSolve(), but for complex matrices.  |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * condition number estimation                                    |
//| * iterative refinement                                           |
//| * O(N^3) complexity                                              |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1,0..N-1], system matrix              |
//|     N       -   size of A                                        |
//|     B       -   array[0..N-1], right part                        |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve                          |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::CMatrixSolve(CMatrixComplex &a,const int n,
                                       al_complex &b[],int &info,
                                       CDenseSolverReport &rep,al_complex &x[])
  {
//--- create a variable
   int i_=0;
//--- create matrix
   CMatrixComplex bm;
   CMatrixComplex xm;
//--- initialization
   info=0;
//--- check
   if(n<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   bm.Resize(n,1);
//--- filling
   for(i_=0;i_<=n-1;i_++)
      bm[i_].Set(0,b[i_]);
//--- function call
   CMatrixSolveM(a,n,bm,1,true,info,rep,xm);
//--- allocation
   ArrayResizeAL(x,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      x[i_]=xm[i_][0];
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixLUSolveM(), but for complex         |
//| matrices.                                                        |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * O(M*N^2) complexity                                            |
//| * condition number estimation                                    |
//| No iterative refinement is provided because exact form of        |
//| original matrix is not known to subroutine. Use CMatrixSolve or  |
//| CMatrixMixedSolve if you need iterative refinement.              |
//| INPUT PARAMETERS                                                 |
//|     LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU|
//|                 result                                           |
//|     P       -   array[0..N-1], pivots array, RMatrixLU result    |
//|     N       -   size of A                                        |
//|     B       -   array[0..N-1,0..M-1], right part                 |
//|     M       -   right part size                                  |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve                          |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::CMatrixLUSolveM(CMatrixComplex &lua,int &p[],
                                          const int n,CMatrixComplex &b,
                                          const int m,int &info,
                                          CDenseSolverReport &rep,
                                          CMatrixComplex &x)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double scalea=0;
//--- create matrix
   CMatrixComplex emptya;
//--- initialization
   info=0;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
//--- 1. scale matrix,max(|U[i][j]|)
//---    we assume that LU is in its normal form,i.e. |L[i][j]|<=1
//--- 2. solve
   scalea=0;
   for(i=0;i<=n-1;i++)
     {
      for(j=i;j<=n-1;j++)
         scalea=MathMax(scalea,CMath::AbsComplex(lua[i][j]));
     }
//--- check
   if(scalea==0.0)
      scalea=1;
//--- change values
   scalea=1/scalea;
//--- function call
   CMatrixLUSolveInternal(lua,p,scalea,n,emptya,false,b,m,info,rep,x);
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixLUSolve(), but for complex matrices.|
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * O(N^2) complexity                                              |
//| * condition number estimation                                    |
//| No iterative refinement is provided because exact form of        |
//| original matrix is not known to subroutine. Use CMatrixSolve or  |
//| CMatrixMixedSolve if you need iterative refinement.              |
//| INPUT PARAMETERS                                                 |
//|     LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU|
//|                 result                                           |
//|     P       -   array[0..N-1], pivots array, CMatrixLU result    |
//|     N       -   size of A                                        |
//|     B       -   array[0..N-1], right part                        |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve                          |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::CMatrixLUSolve(CMatrixComplex &lua,int &p[],
                                         const int n,al_complex &b[],int &info,
                                         CDenseSolverReport &rep,al_complex &x[])
  {
//--- create matrix
   CMatrixComplex bm;
   CMatrixComplex xm;
//--- create a variable
   int i_=0;
//--- initialization
   info=0;
//--- check
   if(n<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   bm.Resize(n,1);
//--- filling
   for(i_=0;i_<=n-1;i_++)
      bm[i_].Set(0,b[i_]);
//--- function call
   CMatrixLUSolveM(lua,p,n,bm,1,info,rep,xm);
//--- allocation
   ArrayResizeAL(x,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      x[i_]=xm[i_][0];
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixMixedSolveM(), but for complex      |
//| matrices.                                                        |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * condition number estimation                                    |
//| * iterative refinement                                           |
//| * O(M*N^2) complexity                                            |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1,0..N-1], system matrix              |
//|     LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU|
//|                 result                                           |
//|     P       -   array[0..N-1], pivots array, CMatrixLU result    |
//|     N       -   size of A                                        |
//|     B       -   array[0..N-1,0..M-1], right part                 |
//|     M       -   right part size                                  |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolveM                         |
//|     Rep     -   same as in RMatrixSolveM                         |
//|     X       -   same as in RMatrixSolveM                         |
//+------------------------------------------------------------------+
static void CDenseSolver::CMatrixMixedSolveM(CMatrixComplex &a,CMatrixComplex &lua,
                                             int &p[],const int n,CMatrixComplex &b,
                                             const int m,int &info,CDenseSolverReport &rep,
                                             CMatrixComplex &x)
  {
//--- create variables
   double scalea=0;
   int    i=0;
   int    j=0;
//--- initialization
   info=0;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
//--- 1. scale matrix,max(|A[i][j]|)
//--- 2. factorize scaled matrix
//--- 3. solve
   scalea=0;
   for(i=0;i<=n-1;i++)
     {
      for(j=0;j<=n-1;j++)
         scalea=MathMax(scalea,CMath::AbsComplex(a[i][j]));
     }
//--- check
   if(scalea==0.0)
      scalea=1;
//--- change values
   scalea=1/scalea;
//--- function call
   CMatrixLUSolveInternal(lua,p,scalea,n,a,true,b,m,info,rep,x);
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixMixedSolve(), but for complex       |
//| matrices.                                                        |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * condition number estimation                                    |
//| * iterative refinement                                           |
//| * O(N^2) complexity                                              |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1,0..N-1], system matrix              |
//|     LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU|
//|                 result                                           |
//|     P       -   array[0..N-1], pivots array, CMatrixLU result    |
//|     N       -   size of A                                        |
//|     B       -   array[0..N-1], right part                        |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolveM                         |
//|     Rep     -   same as in RMatrixSolveM                         |
//|     X       -   same as in RMatrixSolveM                         |
//+------------------------------------------------------------------+
static void CDenseSolver::CMatrixMixedSolve(CMatrixComplex &a,CMatrixComplex &lua,
                                            int &p[],const int n,al_complex &b[],
                                            int &info,CDenseSolverReport &rep,
                                            al_complex &x[])
  {
//--- create matrix
   CMatrixComplex bm;
   CMatrixComplex xm;
//--- create a variable
   int i_=0;
//--- initialization
   info=0;
//--- check
   if(n<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   bm.Resize(n,1);
//--- filling
   for(i_=0;i_<=n-1;i_++)
      bm[i_].Set(0,b[i_]);
//--- function call
   CMatrixMixedSolveM(a,lua,p,n,bm,1,info,rep,xm);
//--- allocation
   ArrayResizeAL(x,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      x[i_]=xm[i_][0];
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixSolveM(), but for symmetric positive|
//| definite matrices.                                               |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * condition number estimation                                    |
//| * O(N^3+M*N^2) complexity                                        |
//| * matrix is represented by its upper or lower triangle           |
//| No iterative refinement is provided because such partial         |
//| representation of matrix does not allow efficient calculation of |
//| extra-precise matrix-vector products for large matrices. Use     |
//| RMatrixSolve or RMatrixMixedSolve if you need iterative          |
//| refinement.                                                      |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1,0..N-1], system matrix              |
//|     N       -   size of A                                        |
//|     IsUpper -   what half of A is provided                       |
//|     B       -   array[0..N-1,0..M-1], right part                 |
//|     M       -   right part size                                  |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve.                         |
//|                 Returns -3 for non-SPD matrices.                 |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::SPDMatrixSolveM(CMatrixDouble &a,const int n,
                                          const bool isupper,CMatrixDouble &b,
                                          const int m,int &info,
                                          CDenseSolverReport &rep,
                                          CMatrixDouble &x)
  {
//--- create variables
   double sqrtscalea=0;
   int    i=0;
   int    j=0;
   int    j1=0;
   int    j2=0;
   int    i_=0;
//--- create matrix
   CMatrixDouble da;
//--- initialization
   info=0;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   da.Resize(n,n);
//--- 1. scale matrix,max(|A[i][j]|)
//--- 2. factorize scaled matrix
//--- 3. solve
   sqrtscalea=0;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(isupper)
        {
         j1=i;
         j2=n-1;
        }
      else
        {
         j1=0;
         j2=i;
        }
      //--- calculation
      for(j=j1;j<=j2;j++)
         sqrtscalea=MathMax(sqrtscalea,MathAbs(a[i][j]));
     }
//--- check
   if(sqrtscalea==0.0)
      sqrtscalea=1;
//--- change values
   sqrtscalea=1/sqrtscalea;
   sqrtscalea=MathSqrt(sqrtscalea);
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(isupper)
        {
         j1=i;
         j2=n-1;
        }
      else
        {
         j1=0;
         j2=i;
        }
      //--- calculation
      for(i_=j1;i_<=j2;i_++)
         da[i].Set(i_,a[i][i_]);
     }
//--- check
   if(!CTrFac::SPDMatrixCholesky(da,n,isupper))
     {
      //--- allocation
      x.Resize(n,m);
      for(i=0;i<=n-1;i++)
        {
         for(j=0;j<=m-1;j++)
            x[i].Set(j,0);
        }
      //--- change values
      rep.m_r1=0;
      rep.m_rinf=0;
      info=-3;
      //--- exit the function
      return;
     }
   info=1;
//--- function call
   SPDMatrixCholeskySolveInternal(da,sqrtscalea,n,isupper,a,true,b,m,info,rep,x);
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixSolve(), but for SPD matrices.      |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * condition number estimation                                    |
//| * O(N^3) complexity                                              |
//| * matrix is represented by its upper or lower triangle           |
//| No iterative refinement is provided because such partial         |
//| representation of matrix does not allow efficient calculation of |
//| extra-precise  matrix-vector products for large matrices. Use    |
//| RMatrixSolve or RMatrixMixedSolve if you need iterative          |
//| refinement.                                                      |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1,0..N-1], system matrix              |
//|     N       -   size of A                                        |
//|     IsUpper -   what half of A is provided                       |
//|     B       -   array[0..N-1], right part                        |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve                          |
//|                 Returns -3 for non-SPD matrices.                 |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::SPDMatrixSolve(CMatrixDouble &a,const int n,
                                         const bool isupper,double &b[],
                                         int &info,CDenseSolverReport &rep,
                                         double &x[])
  {
//--- create a variable
   int i_=0;
//--- create matrix
   CMatrixDouble bm;
   CMatrixDouble xm;
//--- initialization
   info=0;
//--- check
   if(n<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   bm.Resize(n,1);
//--- filling
   for(i_=0;i_<=n-1;i_++)
      bm[i_].Set(0,b[i_]);
//--- function call
   SPDMatrixSolveM(a,n,isupper,bm,1,info,rep,xm);
//--- allocation
   ArrayResizeAL(x,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      x[i_]=xm[i_][0];
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixLUSolveM(), but for SPD matrices    |
//| represented by their Cholesky decomposition.                     |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * O(M*N^2) complexity                                            |
//| * condition number estimation                                    |
//| * matrix is represented by its upper or lower triangle           |
//| No iterative refinement is provided because such partial         |
//| representation of matrix does not allow efficient calculation of |
//| extra-precise  matrix-vector products for large matrices. Use    |
//| RMatrixSolve or RMatrixMixedSolve  if  you need iterative        |
//| refinement.                                                      |
//| INPUT PARAMETERS                                                 |
//|     CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,    |
//|                 SPDMatrixCholesky result                         |
//|     N       -   size of CHA                                      |
//|     IsUpper -   what half of CHA is provided                     |
//|     B       -   array[0..N-1,0..M-1], right part                 |
//|     M       -   right part size                                  |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve                          |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::SPDMatrixCholeskySolveM(CMatrixDouble &cha,const int n,
                                                  const bool isupper,CMatrixDouble &b,
                                                  const int m,int &info,
                                                  CDenseSolverReport &rep,
                                                  CMatrixDouble &x)
  {
//--- create variables
   double sqrtscalea=0;
   int    i=0;
   int    j=0;
   int    j1=0;
   int    j2=0;
//--- create matrix
   CMatrixDouble emptya;
//--- initialization
   info=0;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
//--- 1. scale matrix,max(|U[i][j]|)
//--- 2. factorize scaled matrix
//--- 3. solve
   sqrtscalea=0;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(isupper)
        {
         j1=i;
         j2=n-1;
        }
      else
        {
         j1=0;
         j2=i;
        }
      //--- calculation
      for(j=j1;j<=j2;j++)
         sqrtscalea=MathMax(sqrtscalea,MathAbs(cha[i][j]));
     }
//--- check
   if(sqrtscalea==0.0)
      sqrtscalea=1;
//--- change values
   sqrtscalea=1/sqrtscalea;
//--- function call
   SPDMatrixCholeskySolveInternal(cha,sqrtscalea,n,isupper,emptya,false,b,m,info,rep,x);
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixLUSolve(), but for  SPD matrices    |
//| represented by their Cholesky decomposition.                     |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * O(N^2) complexity                                              |
//| * condition number estimation                                    |
//| * matrix is represented by its upper or lower triangle           |
//| No iterative refinement is provided because such partial         |
//| representation of matrix does not allow efficient calculation of |
//| extra-precise  matrix-vector products for large matrices. Use    |
//| RMatrixSolve or RMatrixMixedSolve  if  you need iterative        |
//| refinement.                                                      |
//| INPUT PARAMETERS                                                 |
//|     CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,    |
//|                 SPDMatrixCholesky result                         |
//|     N       -   size of A                                        |
//|     IsUpper -   what half of CHA is provided                     |
//|     B       -   array[0..N-1], right part                        |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve                          |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::SPDMatrixCholeskySolve(CMatrixDouble &cha,const int n,
                                                 const bool isupper,double &b[],
                                                 int &info,CDenseSolverReport &rep,
                                                 double &x[])
  {
//--- create a variable
   int i_=0;
//--- create matrix
   CMatrixDouble bm;
   CMatrixDouble xm;
//--- initialization
   info=0;
//--- check
   if(n<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   bm.Resize(n,1);
//--- filling
   for(i_=0;i_<=n-1;i_++)
      bm[i_].Set(0,b[i_]);
//--- function call
   SPDMatrixCholeskySolveM(cha,n,isupper,bm,1,info,rep,xm);
//--- allocation
   ArrayResizeAL(x,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      x[i_]=xm[i_][0];
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixSolveM(), but for Hermitian positive|
//| definite matrices.                                               |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * condition number estimation                                    |
//| * O(N^3+M*N^2) complexity                                        |
//| * matrix is represented by its upper or lower triangle           |
//| No iterative refinement is provided because such partial         |
//| representation of matrix does not allow efficient calculation of |
//| extra-precise  matrix-vector products for large matrices. Use    |
//| RMatrixSolve or RMatrixMixedSolve if you need iterative          |
//| refinement.                                                      |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1,0..N-1], system matrix              |
//|     N       -   size of A                                        |
//|     IsUpper -   what half of A is provided                       |
//|     B       -   array[0..N-1,0..M-1], right part                 |
//|     M       -   right part size                                  |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve.                         |
//|                 Returns -3 for non-HPD matrices.                 |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::HPDMatrixSolveM(CMatrixComplex &a,const int n,
                                          const bool isupper,CMatrixComplex &b,
                                          const int m,int &info,
                                          CDenseSolverReport &rep,CMatrixComplex &x)
  {
//--- create variables
   double sqrtscalea=0;
   int    i=0;
   int    j=0;
   int    j1=0;
   int    j2=0;
   int    i_=0;
//--- create matrix
   CMatrixComplex da;
//--- initialization
   info=0;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   da.Resize(n,n);
//--- 1. scale matrix,max(|A[i][j]|)
//--- 2. factorize scaled matrix
//--- 3. solve
   sqrtscalea=0;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(isupper)
        {
         j1=i;
         j2=n-1;
        }
      else
        {
         j1=0;
         j2=i;
        }
      //--- calculation
      for(j=j1;j<=j2;j++)
         sqrtscalea=MathMax(sqrtscalea,CMath::AbsComplex(a[i][j]));
     }
//--- check
   if(sqrtscalea==0.0)
      sqrtscalea=1;
//--- change values
   sqrtscalea=1/sqrtscalea;
   sqrtscalea=MathSqrt(sqrtscalea);
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(isupper)
        {
         j1=i;
         j2=n-1;
        }
      else
        {
         j1=0;
         j2=i;
        }
      //--- calculation
      for(i_=j1;i_<=j2;i_++)
         da[i].Set(i_,a[i][i_]);
     }
//--- check
   if(!CTrFac::HPDMatrixCholesky(da,n,isupper))
     {
      //--- allocation
      x.Resize(n,m);
      for(i=0;i<=n-1;i++)
        {
         for(j=0;j<=m-1;j++)
            x[i].Set(j,0);
        }
      //--- change values
      rep.m_r1=0;
      rep.m_rinf=0;
      info=-3;
      //--- exit the function
      return;
     }
   info=1;
//--- function call
   HPDMatrixCholeskySolveInternal(da,sqrtscalea,n,isupper,a,true,b,m,info,rep,x);
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixSolve(), but for Hermitian positive |
//| definite matrices.                                               |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * condition number estimation                                    |
//| * O(N^3) complexity                                              |
//| * matrix is represented by its upper or lower triangle           |
//| No iterative refinement is provided because such partial         |
//| representation of matrix does not allow efficient calculation of |
//| extra-precise matrix-vector products for large matrices. Use     |
//| RMatrixSolve or RMatrixMixedSolve if you need iterative          |
//| refinement.                                                      |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1,0..N-1], system matrix              |
//|     N       -   size of A                                        |
//|     IsUpper -   what half of A is provided                       |
//|     B       -   array[0..N-1], right part                        |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve                          |
//|                 Returns -3 for non-HPD matrices.                 |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::HPDMatrixSolve(CMatrixComplex &a,const int n,
                                         const bool isupper,al_complex &b[],
                                         int &info,CDenseSolverReport &rep,
                                         al_complex &x[])
  {
//--- create a variable
   int i_=0;
//--- create matrix
   CMatrixComplex bm;
   CMatrixComplex xm;
//--- initialization
   info=0;
//--- check
   if(n<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   bm.Resize(n,1);
//--- filling
   for(i_=0;i_<=n-1;i_++)
      bm[i_].Set(0,b[i_]);
//--- function call
   HPDMatrixSolveM(a,n,isupper,bm,1,info,rep,xm);
//--- allocation
   ArrayResizeAL(x,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      x[i_]=xm[i_][0];
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixLUSolveM(), but for HPD matrices    |
//| represented by their Cholesky decomposition.                     |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * O(M*N^2) complexity                                            |
//| * condition number estimation                                    |
//| * matrix is represented by its upper or lower triangle           |
//| No iterative refinement is provided because such partial         |
//| representation of matrix does not allow efficient calculation of |
//| extra-precise matrix-vector products for large matrices. Use     |
//| RMatrixSolve or RMatrixMixedSolve if you need iterative          |
//| refinement.                                                      |
//| INPUT PARAMETERS                                                 |
//|     CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,    |
//|                 HPDMatrixCholesky result                         |
//|     N       -   size of CHA                                      |
//|     IsUpper -   what half of CHA is provided                     |
//|     B       -   array[0..N-1,0..M-1], right part                 |
//|     M       -   right part size                                  |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve                          |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::HPDMatrixCholeskySolveM(CMatrixComplex &cha,
                                                  const int n,const bool isupper,
                                                  CMatrixComplex &b,const int m,
                                                  int &info,CDenseSolverReport &rep,
                                                  CMatrixComplex &x)
  {
//--- create variables
   double sqrtscalea=0;
   int    i=0;
   int    j=0;
   int    j1=0;
   int    j2=0;
//--- create matrix
   CMatrixComplex emptya;
//--- initialization
   info=0;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
//--- 1. scale matrix,max(|U[i][j]|)
//--- 2. factorize scaled matrix
//--- 3. solve
   sqrtscalea=0;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(isupper)
        {
         j1=i;
         j2=n-1;
        }
      else
        {
         j1=0;
         j2=i;
        }
      //--- calculation
      for(j=j1;j<=j2;j++)
        {
         sqrtscalea=MathMax(sqrtscalea,CMath::AbsComplex(cha[i][j]));
        }
     }
//--- check
   if(sqrtscalea==0.0)
     {
      sqrtscalea=1;
     }
//--- change values
   sqrtscalea=1/sqrtscalea;
//--- function call
   HPDMatrixCholeskySolveInternal(cha,sqrtscalea,n,isupper,emptya,false,b,m,info,rep,x);
  }
//+------------------------------------------------------------------+
//| Dense solver. Same as RMatrixLUSolve(), but for  HPD matrices    |
//| represented by their Cholesky decomposition.                     |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * O(N^2) complexity                                              |
//| * condition number estimation                                    |
//| * matrix is represented by its upper or lower triangle           |
//| No iterative refinement is provided because such partial         |
//| representation of matrix does not allow efficient calculation of |
//| extra-precise  matrix-vector products for large matrices. Use    |
//| RMatrixSolve or RMatrixMixedSolve if you need iterative          |
//| refinement.                                                      |
//| INPUT PARAMETERS                                                 |
//|     CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,    |
//|                 SPDMatrixCholesky result                         |
//|     N       -   size of A                                        |
//|     IsUpper -   what half of CHA is provided                     |
//|     B       -   array[0..N-1], right part                        |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   same as in RMatrixSolve                          |
//|     Rep     -   same as in RMatrixSolve                          |
//|     X       -   same as in RMatrixSolve                          |
//+------------------------------------------------------------------+
static void CDenseSolver::HPDMatrixCholeskySolve(CMatrixComplex &cha,
                                                 const int n,const bool isupper,
                                                 al_complex &b[],int &info,
                                                 CDenseSolverReport &rep,
                                                 al_complex &x[])
  {
//--- create a variable
   int i_=0;
//--- create matrix
   CMatrixComplex bm;
   CMatrixComplex xm;
//--- initialization
   info=0;
//--- check
   if(n<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   bm.Resize(n,1);
//--- filling
   for(i_=0;i_<=n-1;i_++)
      bm[i_].Set(0,b[i_]);
//--- function call
   HPDMatrixCholeskySolveM(cha,n,isupper,bm,1,info,rep,xm);
//--- allocation
   ArrayResizeAL(x,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      x[i_]=xm[i_][0];
  }
//+------------------------------------------------------------------+
//| Dense solver.                                                    |
//| This subroutine finds solution of the linear system A*X=B with   |
//| non-square, possibly degenerate A. System is solved in the least |
//| squares sense, and general least squares solution  X = X0 + CX*y |
//| which  minimizes |A*X-B| is returned. If A is non-degenerate,    |
//| solution in the usual sense is returned                          |
//| Algorithm features:                                              |
//| * automatic detection of degenerate cases                        |
//| * iterative refinement                                           |
//| * O(N^3) complexity                                              |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..NRows-1,0..NCols-1], system matrix      |
//|     NRows   -   vertical size of A                               |
//|     NCols   -   horizontal size of A                             |
//|     B       -   array[0..NCols-1], right part                    |
//|     Threshold-  a number in [0,1]. Singular values beyond        |
//|                 Threshold are considered  zero.  Set it to 0.0,  |
//|                 if you don't understand what it means, so the    |
//|                 solver will choose good value on its own.        |
//| OUTPUT PARAMETERS                                                |
//|     Info    -   return code:                                     |
//|                 * -4    SVD subroutine failed                    |
//|                 * -1    if NRows<=0 or NCols<=0 or Threshold<0   |
//|                         was passed                               |
//|                 *  1    if task is solved                        |
//|     Rep     -   solver report, see below for more info           |
//|     X       -   array[0..N-1,0..M-1], it contains:               |
//|                 * solution of A*X=B if A is non-singular         |
//|                   (well-conditioned or ill-conditioned, but not  |
//|                   very close to singular)                        |
//|                 * zeros, if A is singular or VERY close to       |
//|                   singular (in this case Info=-3).               |
//| SOLVER REPORT                                                    |
//| Subroutine sets following fields of the Rep structure:           |
//| * R2        reciprocal of condition number: 1/cond(A), 2-norm.   |
//| * N         = NCols                                              |
//| * K         dim(Null(A))                                         |
//| * CX        array[0..N-1,0..K-1], kernel of A.                   |
//|             Columns of CX store such vectors that A*CX[i]=0.     |
//+------------------------------------------------------------------+
static void CDenseSolver::RMatrixSolveLS(CMatrixDouble &a,const int nrows,
                                         const int ncols,double &b[],
                                         double threshold,int &info,
                                         CDenseSolverLSReport &rep,
                                         double &x[])
  {
//--- create matrix
   CMatrixDouble u;
   CMatrixDouble vt;
//--- create arrays
   double sv[];
   double rp[];
   double utb[];
   double sutb[];
   double tmp[];
   double ta[];
   double tx[];
   double buf[];
   double w[];
//--- create variables
   int    i=0;
   int    j=0;
   int    nsv=0;
   int    kernelidx=0;
   double v=0;
   double verr=0;
   bool   svdfailed;
   bool   zeroa;
   int    rfs=0;
   int    nrfs=0;
   bool   terminatenexttime;
   bool   smallerr;
   int    i_=0;
//--- initialization
   info=0;
//--- check
   if((nrows<=0 || ncols<=0) || threshold<0.0)
     {
      info=-1;
      return;
     }
//--- check
   if(threshold==0.0)
      threshold=1000*CMath::m_machineepsilon;
//--- Factorize A first
   svdfailed=!CSingValueDecompose::RMatrixSVD(a,nrows,ncols,1,2,2,sv,u,vt);
//--- check
   if(sv[0]==0.0)
      zeroa=true;
   else
      zeroa=false;
//--- check
   if(svdfailed || zeroa)
     {
      //--- check
      if(svdfailed)
         info=-4;
      else
         info=1;
      //--- allocation
      ArrayResizeAL(x,ncols);
      for(i=0;i<=ncols-1;i++)
         x[i]=0;
      //--- change values
      rep.m_n=ncols;
      rep.m_k=ncols;
      rep.m_cx.Resize(ncols,ncols);
      for(i=0;i<=ncols-1;i++)
        {
         for(j=0;j<=ncols-1;j++)
           {
            //--- check
            if(i==j)
               rep.m_cx[i].Set(j,1);
            else
               rep.m_cx[i].Set(j,0);
           }
        }
      rep.m_r2=0;
      //--- exit the function
      return;
     }
   nsv=MathMin(ncols,nrows);
//--- check
   if(nsv==ncols)
      rep.m_r2=sv[nsv-1]/sv[0];
   else
      rep.m_r2=0;
//--- change values
   rep.m_n=ncols;
   info=1;
//--- Iterative refinement of xc combined with solution:
//--- 1. xc=0
//--- 2. calculate r=bc-A*xc using extra-precise dot product
//--- 3. solve A*y=r
//--- 4. update x:=x+r
//--- 5. goto 2
//--- This cycle is executed until one of two things happens:
//--- 1. maximum number of iterations reached
//--- 2. last iteration decreased error to the lower limit
   ArrayResizeAL(utb,nsv);
   ArrayResizeAL(sutb,nsv);
   ArrayResizeAL(x,ncols);
   ArrayResizeAL(tmp,ncols);
   ArrayResizeAL(ta,ncols+1);
   ArrayResizeAL(tx,ncols+1);
   ArrayResizeAL(buf,ncols+1);
//--- initialization
   for(i=0;i<=ncols-1;i++)
      x[i]=0;
   kernelidx=nsv;
   for(i=0;i<=nsv-1;i++)
     {
      //--- check
      if(sv[i]<=threshold*sv[0])
        {
         kernelidx=i;
         break;
        }
     }
//--- change values
   rep.m_k=ncols-kernelidx;
   nrfs=CDenseSolverRFSMaxV2(ncols,rep.m_r2);
   terminatenexttime=false;
//--- allocation
   ArrayResizeAL(rp,nrows);
   for(rfs=0;rfs<=nrfs;rfs++)
     {
      //--- check
      if(terminatenexttime)
         break;
      //--- calculate right part
      if(rfs==0)
        {
         for(i_=0;i_<=nrows-1;i_++)
            rp[i_]=b[i_];
        }
      else
        {
         smallerr=true;
         for(i=0;i<=nrows-1;i++)
           {
            //--- copy
            for(i_=0;i_<=ncols-1;i_++)
               ta[i_]=a[i][i_];
            ta[ncols]=-1;
            //--- copy
            for(i_=0;i_<=ncols-1;i_++)
               tx[i_]=x[i_];
            tx[ncols]=b[i];
            //--- function call
            CXblas::XDot(ta,tx,ncols+1,buf,v,verr);
            rp[i]=-v;
            smallerr=smallerr && MathAbs(v)<4*verr;
           }
         //--- check
         if(smallerr)
            terminatenexttime=true;
        }
      //--- solve A*dx=rp
      for(i=0;i<=ncols-1;i++)
         tmp[i]=0;
      for(i=0;i<=nsv-1;i++)
         utb[i]=0;
      //--- change values
      for(i=0;i<=nrows-1;i++)
        {
         v=rp[i];
         for(i_=0;i_<=nsv-1;i_++)
            utb[i_]=utb[i_]+v*u[i][i_];
        }
      for(i=0;i<=nsv-1;i++)
        {
         //--- check
         if(i<kernelidx)
            sutb[i]=utb[i]/sv[i];
         else
            sutb[i]=0;
        }
      //--- change values
      for(i=0;i<=nsv-1;i++)
        {
         v=sutb[i];
         for(i_=0;i_<=ncols-1;i_++)
            tmp[i_]=tmp[i_]+v*vt[i][i_];
        }
      //--- update x:  x:=x+dx
      for(i_=0;i_<=ncols-1;i_++)
         x[i_]=x[i_]+tmp[i_];
     }
//--- fill CX
   if(rep.m_k>0)
     {
      //--- allocation
      rep.m_cx.Resize(ncols,rep.m_k);
      for(i=0;i<=rep.m_k-1;i++)
        {
         for(i_=0;i_<=ncols-1;i_++)
            rep.m_cx[i_].Set(i,vt[kernelidx+i][i_]);
        }
     }
  }
//+------------------------------------------------------------------+
//| Internal LU solver                                               |
//+------------------------------------------------------------------+
static void CDenseSolver::RMatrixLUSolveInternal(CMatrixDouble &lua,int &p[],
                                                 const double scalea,const int n,
                                                 CMatrixDouble &a,const bool havea,
                                                 CMatrixDouble &b,const int m,
                                                 int &info,CDenseSolverReport &rep,
                                                 CMatrixDouble &x)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    rfs=0;
   int    nrfs=0;
   double v=0;
   double verr=0;
   double mxb=0;
   double scaleright=0;
   bool   smallerr;
   bool   terminatenexttime;
   int    i_=0;
//--- create arrays
   double xc[];
   double y[];
   double bc[];
   double xa[];
   double xb[];
   double tx[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(scalea>0.0))
      return;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(p[i]>n-1 || p[i]<i)
        {
         info=-1;
         return;
        }
     }
//--- allocation
   x.Resize(n,m);
   ArrayResizeAL(y,n);
   ArrayResizeAL(xc,n);
   ArrayResizeAL(bc,n);
   ArrayResizeAL(tx,n+1);
   ArrayResizeAL(xa,n+1);
   ArrayResizeAL(xb,n+1);
//--- estimate condition number,test for near singularity
   rep.m_r1=CRCond::RMatrixLURCond1(lua,n);
   rep.m_rinf=CRCond::RMatrixLURCondInf(lua,n);
//--- check
   if(rep.m_r1<CRCond::RCondThreshold() || rep.m_rinf<CRCond::RCondThreshold())
     {
      for(i=0;i<=n-1;i++)
        {
         for(j=0;j<=m-1;j++)
            x[i].Set(j,0);
        }
      //--- change values
      rep.m_r1=0;
      rep.m_rinf=0;
      info=-3;
      //--- exit the function
      return;
     }
   info=1;
//--- solve
   for(k=0;k<=m-1;k++)
     {
      //--- copy B to contiguous storage
      for(i_=0;i_<=n-1;i_++)
         bc[i_]=b[i_][k];
      //--- Scale right part:
      //--- * MX stores max(|Bi|)
      //--- * ScaleRight stores actual scaling applied to B when solving systems
      //---   it is chosen to make |scaleRight*b| close to 1.
      mxb=0;
      for(i=0;i<=n-1;i++)
         mxb=MathMax(mxb,MathAbs(bc[i]));
      //--- check
      if(mxb==0.0)
         mxb=1;
      //--- change values
      scaleright=1/mxb;
      //--- First,non-iterative part of solution process.
      //--- We use separate code for this task because
      //--- XDot is quite slow and we want to save time.
      for(i_=0;i_<=n-1;i_++)
         xc[i_]=bc[i_]*scaleright;
      //--- function call
      RBasicLUSolve(lua,p,scalea,n,xc,tx);
      //--- Iterative refinement of xc:
      //--- * calculate r=bc-A*xc using extra-precise dot product
      //--- * solve A*y=r
      //--- * update x:=x+r
      //--- This cycle is executed until one of two things happens:
      //--- 1. maximum number of iterations reached
      //--- 2. last iteration decreased error to the lower limit
      if(havea)
        {
         //--- calculation
         nrfs=CDenseSolverRFSMax(n,rep.m_r1,rep.m_rinf);
         terminatenexttime=false;
         for(rfs=0;rfs<=nrfs-1;rfs++)
           {
            //--- check
            if(terminatenexttime)
               break;
            //--- generate right part
            smallerr=true;
            for(i_=0;i_<=n-1;i_++)
               xb[i_]=xc[i_];
            for(i=0;i<=n-1;i++)
              {
               //--- change values
               for(i_=0;i_<=n-1;i_++)
                  xa[i_]=a[i][i_]*scalea;
               xa[n]=-1;
               xb[n]=bc[i]*scaleright;
               //--- function call
               CXblas::XDot(xa,xb,n+1,tx,v,verr);
               y[i]=-v;
               smallerr=smallerr && MathAbs(v)<4*verr;
              }
            //--- check
            if(smallerr)
               terminatenexttime=true;
            //--- solve and update
            RBasicLUSolve(lua,p,scalea,n,y,tx);
            for(i_=0;i_<=n-1;i_++)
               xc[i_]=xc[i_]+y[i_];
           }
        }
      //--- Store xc.
      //--- Post-scale result.
      v=scalea*mxb;
      for(i_=0;i_<=n-1;i_++)
         x[i_].Set(k,xc[i_]*v);
     }
  }
//+------------------------------------------------------------------+
//| Internal Cholesky solver                                         |
//+------------------------------------------------------------------+
static void CDenseSolver::SPDMatrixCholeskySolveInternal(CMatrixDouble &cha,
                                                         const double sqrtscalea,
                                                         const int n,const bool isupper,
                                                         CMatrixDouble &a,
                                                         const bool havea,
                                                         CMatrixDouble &b,
                                                         const int m,int &info,
                                                         CDenseSolverReport &rep,
                                                         CMatrixDouble &x)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   double v=0;
   double mxb=0;
   double scaleright=0;
   int    i_=0;
//--- create arrays
   double xc[];
   double y[];
   double bc[];
   double xa[];
   double xb[];
   double tx[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(sqrtscalea>0.0))
      return;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   x.Resize(n,m);
   ArrayResizeAL(y,n);
   ArrayResizeAL(xc,n);
   ArrayResizeAL(bc,n);
   ArrayResizeAL(tx,n+1);
   ArrayResizeAL(xa,n+1);
   ArrayResizeAL(xb,n+1);
//--- estimate condition number,test for near singularity
   rep.m_r1=CRCond::SPDMatrixCholeskyRCond(cha,n,isupper);
   rep.m_rinf=rep.m_r1;
//--- check
   if(rep.m_r1<CRCond::RCondThreshold())
     {
      for(i=0;i<=n-1;i++)
        {
         for(j=0;j<=m-1;j++)
            x[i].Set(j,0);
        }
      //--- change values
      rep.m_r1=0;
      rep.m_rinf=0;
      info=-3;
      //--- exit the function
      return;
     }
   info=1;
//--- solve
   for(k=0;k<=m-1;k++)
     {
      //--- copy B to contiguous storage
      for(i_=0;i_<=n-1;i_++)
         bc[i_]=b[i_][k];
      //--- Scale right part:
      //--- * MX stores max(|Bi|)
      //--- * ScaleRight stores actual scaling applied to B when solving systems
      //---   it is chosen to make |scaleRight*b| close to 1.
      mxb=0;
      for(i=0;i<=n-1;i++)
         mxb=MathMax(mxb,MathAbs(bc[i]));
      //--- check
      if(mxb==0.0)
         mxb=1;
      scaleright=1/mxb;
      //--- First,non-iterative part of solution process.
      //--- We use separate code for this task because
      //--- XDot is quite slow and we want to save time.
      for(i_=0;i_<=n-1;i_++)
         xc[i_]=bc[i_]*scaleright;
      //--- function call
      SPDBasicCholeskySolve(cha,sqrtscalea,n,isupper,xc,tx);
      //--- Store xc.
      //--- Post-scale result.
      v=CMath::Sqr(sqrtscalea)*mxb;
      for(i_=0;i_<=n-1;i_++)
         x[i_].Set(k,xc[i_]*v);
     }
  }
//+------------------------------------------------------------------+
//| Internal LU solver                                               |
//+------------------------------------------------------------------+
static void CDenseSolver::CMatrixLUSolveInternal(CMatrixComplex &lua,int &p[],
                                                 const double scalea,const int n,
                                                 CMatrixComplex &a,const bool havea,
                                                 CMatrixComplex &b,const int m,
                                                 int &info,CDenseSolverReport &rep,
                                                 CMatrixComplex &x)
  {
//--- create variables
   int        i=0;
   int        j=0;
   int        k=0;
   int        rfs=0;
   int        nrfs=0;
   al_complex v=0;
   double     verr=0;
   double     mxb=0;
   double     scaleright=0;
   bool       smallerr;
   bool       terminatenexttime;
   int        i_=0;
//--- create arrays
   al_complex xc[];
   al_complex y[];
   al_complex bc[];
   al_complex xa[];
   al_complex xb[];
   al_complex tx[];
   double     tmpbuf[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(scalea>0.0))
      return;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(p[i]>n-1 || p[i]<i)
        {
         info=-1;
         return;
        }
     }
//--- allocation
   x.Resize(n,m);
   ArrayResizeAL(y,n);
   ArrayResizeAL(xc,n);
   ArrayResizeAL(bc,n);
   ArrayResizeAL(tx,n);
   ArrayResizeAL(xa,n+1);
   ArrayResizeAL(xb,n+1);
   ArrayResizeAL(tmpbuf,2*n+2);
//--- estimate condition number,test for near singularity
   rep.m_r1=CRCond::CMatrixLURCond1(lua,n);
   rep.m_rinf=CRCond::CMatrixLURCondInf(lua,n);
//--- check
   if(rep.m_r1<CRCond::RCondThreshold() || rep.m_rinf<CRCond::RCondThreshold())
     {
      for(i=0;i<=n-1;i++)
        {
         for(j=0;j<=m-1;j++)
            x[i].Set(j,0);
        }
      //--- change values
      rep.m_r1=0;
      rep.m_rinf=0;
      info=-3;
      //--- exit the function
      return;
     }
   info=1;
//--- solve
   for(k=0;k<=m-1;k++)
     {
      //--- copy B to contiguous storage
      for(i_=0;i_<=n-1;i_++)
         bc[i_]=b[i_][k];
      //--- Scale right part:
      //--- * MX stores max(|Bi|)
      //--- * ScaleRight stores actual scaling applied to B when solving systems
      //---   it is chosen to make |scaleRight*b| close to 1.
      mxb=0;
      for(i=0;i<=n-1;i++)
         mxb=MathMax(mxb,CMath::AbsComplex(bc[i]));
      //--- check
      if(mxb==0.0)
         mxb=1;
      scaleright=1/mxb;
      //--- First,non-iterative part of solution process.
      //--- We use separate code for this task because
      //--- XDot is quite slow and we want to save time.
      for(i_=0;i_<=n-1;i_++)
         xc[i_]=bc[i_]*scaleright;
      //--- function call
      CBasicLUSolve(lua,p,scalea,n,xc,tx);
      //--- Iterative refinement of xc:
      //--- * calculate r=bc-A*xc using extra-precise dot product
      //--- * solve A*y=r
      //--- * update x:=x+r
      //--- This cycle is executed until one of two things happens:
      //--- 1. maximum number of iterations reached
      //--- 2. last iteration decreased error to the lower limit
      if(havea)
        {
         //--- calculation
         nrfs=CDenseSolverRFSMax(n,rep.m_r1,rep.m_rinf);
         terminatenexttime=false;
         for(rfs=0;rfs<=nrfs-1;rfs++)
           {
            //--- check
            if(terminatenexttime)
               break;
            //--- generate right part
            smallerr=true;
            //--- copy
            for(i_=0;i_<=n-1;i_++)
               xb[i_]=xc[i_];
            for(i=0;i<=n-1;i++)
              {
               //--- change values
               for(i_=0;i_<=n-1;i_++)
                  xa[i_]=a[i][i_]*scalea;
               xa[n]=-1;
               xb[n]=bc[i]*scaleright;
               //--- function call
               CXblas::XCDot(xa,xb,n+1,tmpbuf,v,verr);
               y[i]=-v;
               smallerr=smallerr && CMath::AbsComplex(v)<4*verr;
              }
            //--- check
            if(smallerr)
               terminatenexttime=true;
            //--- solve and update
            CBasicLUSolve(lua,p,scalea,n,y,tx);
            for(i_=0;i_<=n-1;i_++)
               xc[i_]=xc[i_]+y[i_];
           }
        }
      //--- Store xc.
      //--- Post-scale result.
      v=scalea*mxb;
      for(i_=0;i_<=n-1;i_++)
         x[i_].Set(k,xc[i_]*v);
     }
  }
//+------------------------------------------------------------------+
//| Internal Cholesky solver                                         |
//+------------------------------------------------------------------+
static void CDenseSolver::HPDMatrixCholeskySolveInternal(CMatrixComplex &cha,
                                                         const double sqrtscalea,
                                                         const int n,
                                                         const bool isupper,
                                                         CMatrixComplex &a,
                                                         const bool havea,
                                                         CMatrixComplex &b,
                                                         const int m,int &info,
                                                         CDenseSolverReport &rep,
                                                         CMatrixComplex &x)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   double v=0;
   double mxb=0;
   double scaleright=0;
   int    i_=0;
//--- create arrays
   al_complex xc[];
   al_complex y[];
   al_complex bc[];
   al_complex xa[];
   al_complex xb[];
   al_complex tx[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(sqrtscalea>0.0))
      return;
//--- prepare: check inputs,allocate space...
   if(n<=0 || m<=0)
     {
      info=-1;
      return;
     }
//--- allocation
   x.Resize(n,m);
   ArrayResizeAL(y,n);
   ArrayResizeAL(xc,n);
   ArrayResizeAL(bc,n);
   ArrayResizeAL(tx,n+1);
   ArrayResizeAL(xa,n+1);
   ArrayResizeAL(xb,n+1);
//--- estimate condition number,test for near singularity
   rep.m_r1=CRCond::HPDMatrixCholeskyRCond(cha,n,isupper);
   rep.m_rinf=rep.m_r1;
//--- check
   if(rep.m_r1<CRCond::RCondThreshold())
     {
      for(i=0;i<=n-1;i++)
        {
         for(j=0;j<=m-1;j++)
            x[i].Set(j,0);
        }
      //--- change values
      rep.m_r1=0;
      rep.m_rinf=0;
      info=-3;
      //--- exit the function
      return;
     }
   info=1;
//--- solve
   for(k=0;k<=m-1;k++)
     {
      //--- copy B to contiguous storage
      for(i_=0;i_<=n-1;i_++)
         bc[i_]=b[i_][k];
      //--- Scale right part:
      //--- * MX stores max(|Bi|)
      //--- * ScaleRight stores actual scaling applied to B when solving systems
      //---   it is chosen to make |scaleRight*b| close to 1.
      mxb=0;
      for(i=0;i<=n-1;i++)
         mxb=MathMax(mxb,CMath::AbsComplex(bc[i]));
      //--- check
      if(mxb==0.0)
         mxb=1;
      scaleright=1/mxb;
      //--- First,non-iterative part of solution process.
      //--- We use separate code for this task because
      //--- XDot is quite slow and we want to save time.
      for(i_=0;i_<=n-1;i_++)
         xc[i_]=bc[i_]*scaleright;
      //--- function call
      HPDBasicCholeskySolve(cha,sqrtscalea,n,isupper,xc,tx);
      //--- Store xc.
      //--- Post-scale result.
      v=CMath::Sqr(sqrtscalea)*mxb;
      for(i_=0;i_<=n-1;i_++)
         x[i_].Set(k,xc[i_]*v);
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine.                                             |
//| Returns maximum count of RFS iterations as function of:          |
//| 1. machine epsilon                                               |
//| 2. task size.                                                    |
//| 3. condition number                                              |
//+------------------------------------------------------------------+
static int CDenseSolver::CDenseSolverRFSMax(const int n,const double r1,
                                            const double rinf)
  {
//--- return result
   return(5);
  }
//+------------------------------------------------------------------+
//| Internal subroutine.                                             |
//| Returns maximum count of RFS iterations as function of:          |
//| 1. machine epsilon                                               |
//| 2. task size.                                                    |
//| 3. norm-2 condition number                                       |
//+------------------------------------------------------------------+
static int CDenseSolver::CDenseSolverRFSMaxV2(const int n,const double r2)
  {
//--- return result
   return(CDenseSolverRFSMax(n,0,0));
  }
//+------------------------------------------------------------------+
//| Basic LU solver for ScaleA*PLU*x = y.                            |
//| This subroutine assumes that:                                    |
//| * L is well-scaled, and it is U which needs scaling by ScaleA.   |
//| * A=PLU is well-conditioned, so no zero divisions or overflow may|
//|   occur                                                          |
//+------------------------------------------------------------------+
static void CDenseSolver::RBasicLUSolve(CMatrixDouble &lua,int &p[],
                                        const double scalea,const int n,
                                        double &xb[],double &tmp[])
  {
//--- create variables
   int    i=0;
   double v=0;
   int    i_=0;
//--- swap
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(p[i]!=i)
        {
         v=xb[i];
         xb[i]=xb[p[i]];
         xb[p[i]]=v;
        }
     }
//--- calculation
   for(i=1;i<=n-1;i++)
     {
      v=0.0;
      //--- change value
      for(i_=0;i_<=i-1;i_++)
         v+=lua[i][i_]*xb[i_];
      //--- shift
      xb[i]=xb[i]-v;
     }
//--- change value
   xb[n-1]=xb[n-1]/(lua[n-1][n-1]*scalea);
   for(i=n-2;i>=0;i--)
     {
      //--- calculation
      for(i_=i+1;i_<=n-1;i_++)
         tmp[i_]=lua[i][i_]*scalea;
      v=0.0;
      //--- change value
      for(i_=i+1;i_<=n-1;i_++)
         v+=tmp[i_]*xb[i_];
      //--- get result
      xb[i]=(xb[i]-v)/(lua[i][i]*scalea);
     }
  }
//+------------------------------------------------------------------+
//| Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.             |
//| This subroutine assumes that:                                    |
//| * A*ScaleA is well scaled                                        |
//| * A is well-conditioned, so no zero divisions or overflow may    |
//| occur                                                            |
//+------------------------------------------------------------------+
static void CDenseSolver::SPDBasicCholeskySolve(CMatrixDouble &cha,
                                                const double sqrtscalea,
                                                const int n,const bool isupper,
                                                double &xb[],double &tmp[])
  {
//--- create variables
   int    i=0;
   double v=0;
   int    i_=0;
//--- A=L*L' or A=U'*U
   if(isupper)
     {
      //--- Solve U'*y=b first.
      for(i=0;i<=n-1;i++)
        {
         xb[i]=xb[i]/(sqrtscalea*cha[i][i]);
         //--- check
         if(i<n-1)
           {
            v=xb[i];
            //--- calculation
            for(i_=i+1;i_<=n-1;i_++)
               tmp[i_]=sqrtscalea*cha[i][i_];
            for(i_=i+1;i_<=n-1;i_++)
               xb[i_]=xb[i_]-v*tmp[i_];
           }
        }
      //--- Solve U*x=y then.
      for(i=n-1;i>=0;i--)
        {
         //--- check
         if(i<n-1)
           {
            for(i_=i+1;i_<=n-1;i_++)
               tmp[i_]=sqrtscalea*cha[i][i_];
            v=0.0;
            //--- change value
            for(i_=i+1;i_<=n-1;i_++)
               v+=tmp[i_]*xb[i_];
            //--- shift
            xb[i]=xb[i]-v;
           }
         xb[i]=xb[i]/(sqrtscalea*cha[i][i]);
        }
     }
   else
     {
      //--- Solve L*y=b first
      for(i=0;i<=n-1;i++)
        {
         //--- check
         if(i>0)
           {
            for(i_=0;i_<=i-1;i_++)
               tmp[i_]=sqrtscalea*cha[i][i_];
            //--- change value
            v=0.0;
            for(i_=0;i_<=i-1;i_++)
               v+=tmp[i_]*xb[i_];
            //--- shift
            xb[i]=xb[i]-v;
           }
         xb[i]=xb[i]/(sqrtscalea*cha[i][i]);
        }
      //--- Solve L'*x=y then.
      for(i=n-1;i>=0;i--)
        {
         xb[i]=xb[i]/(sqrtscalea*cha[i][i]);
         //--- check
         if(i>0)
           {
            v=xb[i];
            //--- calculation
            for(i_=0;i_<=i-1;i_++)
               tmp[i_]=sqrtscalea*cha[i][i_];
            for(i_=0;i_<=i-1;i_++)
               xb[i_]=xb[i_]-v*tmp[i_];
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Basic LU solver for ScaleA*PLU*x = y.                            |
//| This subroutine assumes that:                                    |
//| * L is well-scaled, and it is U which needs scaling by ScaleA.   |
//| * A=PLU is well-conditioned, so no zero divisions or overflow may|
//|   occur                                                          |
//+------------------------------------------------------------------+
static void CDenseSolver::CBasicLUSolve(CMatrixComplex &lua,int &p[],
                                        const double scalea,const int n,
                                        al_complex &xb[],al_complex &tmp[])
  {
//--- create variables
   int        i=0;
   al_complex v=0;
   int        i_=0;
//--- swap
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(p[i]!=i)
        {
         v=xb[i];
         xb[i]=xb[p[i]];
         xb[p[i]]=v;
        }
     }
   for(i=1;i<=n-1;i++)
     {
      v=0.0;
      //--- calculation
      for(i_=0;i_<=i-1;i_++)
         v+=lua[i][i_]*xb[i_];
      //--- shift
      xb[i]=xb[i]-v;
     }
//--- change values
   xb[n-1]=xb[n-1]/(lua[n-1][n-1]*scalea);
   for(i=n-2;i>=0;i--)
     {
      for(i_=i+1;i_<=n-1;i_++)
         tmp[i_]=lua[i][i_]*scalea;
      //--- calculation
      v=0.0;
      for(i_=i+1;i_<=n-1;i_++)
         v+=tmp[i_]*xb[i_];
      //--- get result
      xb[i]=(xb[i]-v)/(lua[i][i]*scalea);
     }
  }
//+------------------------------------------------------------------+
//| Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.             |
//| This subroutine assumes that:                                    |
//| * A*ScaleA is well scaled                                        |
//| * A is well-conditioned, so no zero divisions or overflow may    |
//|   occur                                                          |
//+------------------------------------------------------------------+
static void CDenseSolver::HPDBasicCholeskySolve(CMatrixComplex &cha,
                                                const double sqrtscalea,
                                                const int n,const bool isupper,
                                                al_complex &xb[],al_complex &tmp[])
  {
//--- create variables
   int        i=0;
   al_complex v=0;
   int        i_=0;
//--- A=L*L' or A=U'*U
   if(isupper)
     {
      //--- Solve U'*y=b first.
      for(i=0;i<=n-1;i++)
        {
         xb[i]=xb[i]/(CMath::Conj(cha[i][i])*sqrtscalea);
         //--- check
         if(i<n-1)
           {
            v=xb[i];
            //--- calculation
            for(i_=i+1;i_<=n-1;i_++)
               tmp[i_]=CMath::Conj(cha[i][i_])*sqrtscalea;
            for(i_=i+1;i_<=n-1;i_++)
               xb[i_]=xb[i_]-v*tmp[i_];
           }
        }
      //--- Solve U*x=y then.
      for(i=n-1;i>=0;i--)
        {
         //--- check
         if(i<n-1)
           {
            for(i_=i+1;i_<=n-1;i_++)
               tmp[i_]=cha[i][i_]*sqrtscalea;
            //--- change value
            v=0.0;
            for(i_=i+1;i_<=n-1;i_++)
               v+=tmp[i_]*xb[i_];
            //--- shift
            xb[i]=xb[i]-v;
           }
         xb[i]=xb[i]/(cha[i][i]*sqrtscalea);
        }
     }
   else
     {
      //--- Solve L*y=b first
      for(i=0;i<=n-1;i++)
        {
         //--- check
         if(i>0)
           {
            for(i_=0;i_<=i-1;i_++)
               tmp[i_]=cha[i][i_]*sqrtscalea;
            //--- change value
            v=0.0;
            for(i_=0;i_<=i-1;i_++)
               v+=tmp[i_]*xb[i_];
            //--- shift
            xb[i]=xb[i]-v;
           }
         xb[i]=xb[i]/(cha[i][i]*sqrtscalea);
        }
      //--- Solve L'*x=y then.
      for(i=n-1;i>=0;i--)
        {
         xb[i]=xb[i]/(CMath::Conj(cha[i][i])*sqrtscalea);
         //--- check
         if(i>0)
           {
            v=xb[i];
            //--- calculation
            for(i_=0;i_<=i-1;i_++)
               tmp[i_]=CMath::Conj(cha[i][i_])*sqrtscalea;
            for(i_=0;i_<=i-1;i_++)
               xb[i_]=xb[i_]-v*tmp[i_];
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CNlEq                                        |
//+------------------------------------------------------------------+
class CNlEqState
  {
public:
   //--- variables
   int               m_n;
   int               m_m;
   double            m_epsf;
   int               m_maxits;
   bool              m_xrep;
   double            m_stpmax;
   double            m_f;
   bool              m_needf;
   bool              m_needfij;
   bool              m_xupdated;
   RCommState        m_rstate;
   int               m_repiterationscount;
   int               m_repnfunc;
   int               m_repnjac;
   int               m_repterminationtype;
   double            m_fbase;
   double            m_fprev;
   //--- arrays
   double            m_x[];
   double            m_fi[];
   double            m_xbase[];
   double            m_candstep[];
   double            m_rightpart[];
   double            m_cgbuf[];
   //--- matrix
   CMatrixDouble     m_j;
   //--- constructor, destructor
                     CNlEqState(void);
                    ~CNlEqState(void);
   //--- copy
   void              Copy(CNlEqState &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNlEqState::CNlEqState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNlEqState::~CNlEqState(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CNlEqState::Copy(CNlEqState &obj)
  {
//--- copy variables
   m_n=obj.m_n;
   m_m=obj.m_m;
   m_epsf=obj.m_epsf;
   m_maxits=obj.m_maxits;
   m_xrep=obj.m_xrep;
   m_stpmax=obj.m_stpmax;
   m_f=obj.m_f;
   m_needf=obj.m_needf;
   m_needfij=obj.m_needfij;
   m_xupdated=obj.m_xupdated;
   m_repiterationscount=obj.m_repiterationscount;
   m_repnfunc=obj.m_repnfunc;
   m_repnjac=obj.m_repnjac;
   m_repterminationtype=obj.m_repterminationtype;
   m_fbase=obj.m_fbase;
   m_fprev=obj.m_fprev;
   m_rstate.Copy(obj.m_rstate);
//--- copy arrays
   ArrayCopy(m_x,obj.m_x);
   ArrayCopy(m_fi,obj.m_fi);
   ArrayCopy(m_xbase,obj.m_xbase);
   ArrayCopy(m_candstep,obj.m_candstep);
   ArrayCopy(m_rightpart,obj.m_rightpart);
   ArrayCopy(m_cgbuf,obj.m_cgbuf);
//--- copy matrix
   m_j=obj.m_j;
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CNlEqState                       |
//+------------------------------------------------------------------+
class CNlEqStateShell
  {
private:
   CNlEqState        m_innerobj;
public:
   //--- constructors, destructor
                     CNlEqStateShell(void);
                     CNlEqStateShell(CNlEqState &obj);
                    ~CNlEqStateShell(void);
   //--- methods
   bool              GetNeedF(void);
   void              SetNeedF(const bool b);
   bool              GetNeedFIJ(void);
   void              SetNeedFIJ(const bool b);
   bool              GetXUpdated(void);
   void              SetXUpdated(const bool b);
   double            GetF(void);
   void              SetF(const double d);
   CNlEqState       *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNlEqStateShell::CNlEqStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CNlEqStateShell::CNlEqStateShell(CNlEqState &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNlEqStateShell::~CNlEqStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needf                          |
//+------------------------------------------------------------------+
bool CNlEqStateShell::GetNeedF(void)
  {
//--- return result
   return(m_innerobj.m_needf);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needf                         |
//+------------------------------------------------------------------+
void CNlEqStateShell::SetNeedF(const bool b)
  {
//--- change value
   m_innerobj.m_needf=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needfij                        |
//+------------------------------------------------------------------+
bool CNlEqStateShell::GetNeedFIJ(void)
  {
//--- return result
   return(m_innerobj.m_needfij);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needfij                       |
//+------------------------------------------------------------------+
void CNlEqStateShell::SetNeedFIJ(const bool b)
  {
//--- change value
   m_innerobj.m_needfij=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable xupdated                       |
//+------------------------------------------------------------------+
bool CNlEqStateShell::GetXUpdated(void)
  {
//--- return result
   return(m_innerobj.m_xupdated);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable xupdated                      |
//+------------------------------------------------------------------+
void CNlEqStateShell::SetXUpdated(const bool b)
  {
//--- change value
   m_innerobj.m_xupdated=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable f                              |
//+------------------------------------------------------------------+
double CNlEqStateShell::GetF(void)
  {
//--- return result
   return(m_innerobj.m_f);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable f                             |
//+------------------------------------------------------------------+
void CNlEqStateShell::SetF(const double d)
  {
//--- change value
   m_innerobj.m_f=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CNlEqState *CNlEqStateShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CNlEq                                        |
//+------------------------------------------------------------------+
class CNlEqReport
  {
public:
   //--- variables
   int               m_iterationscount;
   int               m_nfunc;
   int               m_njac;
   int               m_terminationtype;
   //--- constructor, destructor
                     CNlEqReport(void);
                    ~CNlEqReport(void);
   //--- copy
   void              Copy(CNlEqReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNlEqReport::CNlEqReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNlEqReport::~CNlEqReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CNlEqReport::Copy(CNlEqReport &obj)
  {
//--- copy variables
   m_iterationscount=obj.m_iterationscount;
   m_nfunc=obj.m_nfunc;
   m_njac=obj.m_njac;
   m_terminationtype=obj.m_terminationtype;
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CNlEqReport                      |
//+------------------------------------------------------------------+
class CNlEqReportShell
  {
private:
   CNlEqReport       m_innerobj;
public:
   //--- constructors, destructor
                     CNlEqReportShell(void);
                     CNlEqReportShell(CNlEqReport &obj);
                    ~CNlEqReportShell(void);
   //--- methods
   int               GetIterationsCount(void);
   void              SetIterationsCount(const int i);
   int               GetNFunc(void);
   void              SetNFunc(const int i);
   int               GetNJac(void);
   void              SetNJac(const int i);
   int               GetTerminationType(void);
   void              SetTerminationType(const int i);
   CNlEqReport      *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNlEqReportShell::CNlEqReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CNlEqReportShell::CNlEqReportShell(CNlEqReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNlEqReportShell::~CNlEqReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable iterationscount                |
//+------------------------------------------------------------------+
int CNlEqReportShell::GetIterationsCount(void)
  {
//--- return result
   return(m_innerobj.m_iterationscount);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable iterationscount               |
//+------------------------------------------------------------------+
void CNlEqReportShell::SetIterationsCount(const int i)
  {
//--- change value
   m_innerobj.m_iterationscount=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable nfunc                          |
//+------------------------------------------------------------------+
int CNlEqReportShell::GetNFunc(void)
  {
//--- return result
   return(m_innerobj.m_nfunc);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable nfunc                         |
//+------------------------------------------------------------------+
void CNlEqReportShell::SetNFunc(const int i)
  {
//--- change value
   m_innerobj.m_nfunc=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable njac                           |
//+------------------------------------------------------------------+
int CNlEqReportShell::GetNJac(void)
  {
//--- return result
   return(m_innerobj.m_njac);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable njac                          |
//+------------------------------------------------------------------+
void CNlEqReportShell::SetNJac(const int i)
  {
//--- change value
   m_innerobj.m_njac=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable terminationtype                |
//+------------------------------------------------------------------+
int CNlEqReportShell::GetTerminationType(void)
  {
//--- return result
   return(m_innerobj.m_terminationtype);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable terminationtype               |
//+------------------------------------------------------------------+
void CNlEqReportShell::SetTerminationType(const int i)
  {
//--- change value
   m_innerobj.m_terminationtype=i;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CNlEqReport *CNlEqReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Solving systems of nonlinear equations                           |
//+------------------------------------------------------------------+
class CNlEq
  {
private:
   //--- private methods
   static void       ClearRequestFields(CNlEqState &state);
   static bool       IncreaseLambda(double &lambdav,double &nu,const double lambdaup);
   static void       DecreaseLambda(double &lambdav,double &nu,const double lambdadown);
   //--- auxiliary functions forNlEqiteration
   static void       Func_lbl_rcomm(CNlEqState &state,const int n,const int m,const int i,const bool b,const double lambdaup,const double lambdadown,const double lambdav,const double rho,const double mu,const double stepnorm);
   static void       Func_lbl_7(CNlEqState &state,const int n);
   static bool       Func_lbl_5(CNlEqState &state,double &lambdaup,double &lambdadown,double &lambdav,double &rho);
   static bool       Func_lbl_11(CNlEqState &state,const double stepnorm);
   static int        Func_lbl_10(CNlEqState &state,const int n,const int m,const int i,const bool b,const double lambdaup,const double lambdadown,const double lambdav,const double rho,const double mu,const double stepnorm);
   static int        Func_lbl_9(CNlEqState &state,int &n,int &m,int &i,bool &b,const double lambdaup,const double lambdadown,double &lambdav,const double rho,const double mu,double &stepnorm);
public:
   //--- constant
   static const int  m_armijomaxfev;
   //--- constructor, destructor
                     CNlEq(void);
                    ~CNlEq(void);
   //--- public methods
   static void       NlEqCreateLM(const int n,const int m,double &x[],CNlEqState &state);
   static void       NlEqSetCond(CNlEqState &state,double epsf,const int maxits);
   static void       NlEqSetXRep(CNlEqState &state,const bool needxrep);
   static void       NlEqSetStpMax(CNlEqState &state,const double stpmax);
   static void       NlEqResults(CNlEqState &state,double &x[],CNlEqReport &rep);
   static void       NlEqResultsBuf(CNlEqState &state,double &x[],CNlEqReport &rep);
   static void       NlEqRestartFrom(CNlEqState &state,double &x[]);
   static bool       NlEqIteration(CNlEqState &state);
  };
//+------------------------------------------------------------------+
//| Initialize constant                                              |
//+------------------------------------------------------------------+
const int CNlEq::m_armijomaxfev=20;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNlEq::CNlEq(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNlEq::~CNlEq(void)
  {

  }
//+------------------------------------------------------------------+
//|                 LEVENBERG-MARQUARDT-LIKE NONLINEAR SOLVER        |
//| DESCRIPTION:                                                     |
//| This algorithm solves system of nonlinear equations              |
//|     F[0](x[0], ..., x[n-1])   = 0                                |
//|     F[1](x[0], ..., x[n-1])   = 0                                |
//|     ...                                                          |
//|     F[M-1](x[0], ..., x[n-1]) = 0                                |
//| with M/N do not necessarily coincide. Algorithm converges        |
//| quadratically under following conditions:                        |
//|     * the solution set XS is nonempty                            |
//|     * for some xs in XS there exist such neighbourhood N(xs)     |
//|       that:                                                      |
//|       * vector function F(x) and its Jacobian J(x) are           |
//|         continuously differentiable on N                         |
//|       * ||F(x)|| provides local error bound on N, i.e. there     |
//|         exists such c1, that ||F(x)||>c1*distance(x,XS)          |
//| Note that these conditions are much more weaker than usual       |
//| non-singularity conditions. For example, algorithm will converge |
//| for any affine function F (whether its Jacobian singular or not).|
//| REQUIREMENTS:                                                    |
//| Algorithm will request following information during its          |
//| operation:                                                       |
//| * function vector F[] and Jacobian matrix at given point X       |
//| * value of merit function f(x)=F[0]^2(x)+...+F[M-1]^2(x) at given|
//| point X                                                          |
//| USAGE:                                                           |
//| 1. User initializes algorithm state with NLEQCreateLM() call     |
//| 2. User tunes solver parameters with NLEQSetCond(),              |
//|    NLEQSetStpMax() and other functions                           |
//| 3. User calls NLEQSolve() function which takes algorithm state   |
//|    and pointers (delegates, etc.) to callback functions which    |
//|    calculate merit function value and Jacobian.                  |
//| 4. User calls NLEQResults() to get solution                      |
//| 5. Optionally, user may call NLEQRestartFrom() to solve another  |
//|    problem with same parameters (N/M) but another starting point |
//|    and/or another function vector. NLEQRestartFrom() allows to   |
//|    reuse already initialized structure.                          |
//| INPUT PARAMETERS:                                                |
//|     N       -   space dimension, N>1:                            |
//|                 * if provided, only leading N elements of X are  |
//|                   used                                           |
//|                 * if not provided, determined automatically from |
//|                   size of X                                      |
//|     M       -   system size                                      |
//|     X       -   starting point                                   |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//| NOTES:                                                           |
//| 1. you may tune stopping conditions with NLEQSetCond() function  |
//| 2. if target function contains exp() or other fast growing       |
//|    functions, and optimization algorithm makes too large steps   |
//|    which leads to overflow, use NLEQSetStpMax() function to bound|
//|    algorithm's steps.                                            |
//| 3. this algorithm is a slightly modified implementation of the   |
//|    method described in 'Levenberg-Marquardt method for           |
//|    constrained nonlinear equations with strong local convergence |
//|    properties' by Christian Kanzow Nobuo Yamashita and Masao     |
//|    Fukushima and further developed in  'On the convergence of a  |
//|    New Levenberg-Marquardt Method' by Jin-yan Fan and Ya-Xiang   |
//|    Yuan.                                                         |
//+------------------------------------------------------------------+
static void CNlEq::NlEqCreateLM(const int n,const int m,double &x[],
                                CNlEqState &state)
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
//--- Initialize
   state.m_n=n;
   state.m_m=m;
   NlEqSetCond(state,0,0);
   NlEqSetXRep(state,false);
   NlEqSetStpMax(state,0);
//--- allocation
   ArrayResizeAL(state.m_x,n);
   ArrayResizeAL(state.m_xbase,n);
   state.m_j.Resize(m,n);
   ArrayResizeAL(state.m_fi,m);
   ArrayResizeAL(state.m_rightpart,n);
   ArrayResizeAL(state.m_candstep,n);
//--- function call
   NlEqRestartFrom(state,x);
  }
//+------------------------------------------------------------------+
//| This function sets stopping conditions for the nonlinear solver  |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     EpsF    -   >=0                                              |
//|                 The subroutine finishes  its work if on k+1-th   |
//|                 iteration the condition ||F||<=EpsF is satisfied |
//|     MaxIts  -   maximum number of iterations. If MaxIts=0, the   |
//|                 number of iterations is unlimited.               |
//| Passing EpsF=0 and MaxIts=0 simultaneously will lead to          |
//| automatic stopping criterion selection (small EpsF).             |
//| NOTES:                                                           |
//+------------------------------------------------------------------+
static void CNlEq::NlEqSetCond(CNlEqState &state,double epsf,const int maxits)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsf),__FUNCTION__+": EpsF is not finite number!"))
      return;
//--- check
   if(!CAp::Assert(epsf>=0.0,__FUNCTION__+": negative EpsF!"))
      return;
//--- check
   if(!CAp::Assert(maxits>=0,__FUNCTION__+": negative MaxIts!"))
      return;
//--- check
   if(epsf==0.0 && maxits==0)
      epsf=1.0E-6;
//--- change values
   state.m_epsf=epsf;
   state.m_maxits=maxits;
  }
//+------------------------------------------------------------------+
//| This function turns on/off reporting.                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     NeedXRep-   whether iteration reports are needed or not      |
//| If NeedXRep is True, algorithm will call rep() callback function |
//| if it is provided to NLEQSolve().                                |
//+------------------------------------------------------------------+
static void CNlEq::NlEqSetXRep(CNlEqState &state,const bool needxrep)
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
//| Use this subroutine when target function contains exp() or other |
//| fast growing functions, and algorithm makes too large steps which|
//| lead to overflow. This function allows us to reject steps that   |
//| are too large (and therefore expose us to the possible overflow) |
//| without actually calculating function value at the x+stp*d.      |
//+------------------------------------------------------------------+
static void CNlEq::NlEqSetStpMax(CNlEqState &state,const double stpmax)
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
//| NLEQ solver results                                              |
//| INPUT PARAMETERS:                                                |
//|     State   -   algorithm state.                                 |
//| OUTPUT PARAMETERS:                                               |
//|     X       -   array[0..N-1], solution                          |
//|     Rep     -   optimization report:                             |
//|                 * Rep.TerminationType completetion code:         |
//|                     * -4    ERROR: algorithm has converged to the|
//|                             stationary point Xf which is local   |
//|                             minimum of f=F[0]^2+...+F[m-1]^2,    |
//|                             but is not solution of nonlinear     |
//|                             system.                              |
//|                     *  1    sqrt(f)<=EpsF.                       |
//|                     *  5    MaxIts steps was taken               |
//|                     *  7    stopping conditions are too          |
//|                             stringent, further improvement is    |
//|                             impossible                           |
//|                 * Rep.IterationsCount contains iterations count  |
//|                 * NFEV countains number of function calculations |
//|                 * ActiveConstraints contains number of active    |
//|                   constraints                                    |
//+------------------------------------------------------------------+
static void CNlEq::NlEqResults(CNlEqState &state,double &x[],CNlEqReport &rep)
  {
   ArrayResizeAL(x,0);
//--- function call
   NlEqResultsBuf(state,x,rep);
  }
//+------------------------------------------------------------------+
//| NLEQ solver results                                              |
//| Buffered implementation of NLEQResults(), which uses             |
//| pre-allocated buffer to store X[]. If buffer size is too small,  |
//| it resizes buffer. It is intended to be used in the inner cycles |
//| of performance critical algorithms where array reallocation      |
//| penalty is too large to be ignored.                              |
//+------------------------------------------------------------------+
static void CNlEq::NlEqResultsBuf(CNlEqState &state,double &x[],CNlEqReport &rep)
  {
//--- create a variable
   int i_=0;
//--- check
   if(CAp::Len(x)<state.m_n)
      ArrayResizeAL(x,state.m_n);
//--- copy
   for(i_=0;i_<=state.m_n-1;i_++)
      x[i_]=state.m_xbase[i_];
//--- change values
   rep.m_iterationscount=state.m_repiterationscount;
   rep.m_nfunc=state.m_repnfunc;
   rep.m_njac=state.m_repnjac;
   rep.m_terminationtype=state.m_repterminationtype;
  }
//+------------------------------------------------------------------+
//| This subroutine restarts CG algorithm from new point. All        |
//| optimization parameters are left unchanged.                      |
//| This function allows to solve multiple optimization problems     |
//| (which must have same number of dimensions) without object       |
//| reallocation penalty.                                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure used for reverse communication         |
//|                 previously allocated with MinCGCreate call.      |
//|     X       -   new starting point.                              |
//|     BndL    -   new lower bounds                                 |
//|     BndU    -   new upper bounds                                 |
//+------------------------------------------------------------------+
static void CNlEq::NlEqRestartFrom(CNlEqState &state,double &x[])
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
   ArrayResizeAL(state.m_rstate.ia,3);
   ArrayResizeAL(state.m_rstate.ba,1);
   ArrayResizeAL(state.m_rstate.ra,6);
   state.m_rstate.stage=-1;
//--- function call
   ClearRequestFields(state);
  }
//+------------------------------------------------------------------+
//| Clears request fileds (to be sure that we don't forgot to clear  |
//| something)                                                       |
//+------------------------------------------------------------------+
static void CNlEq::ClearRequestFields(CNlEqState &state)
  {
//--- change values
   state.m_needf=false;
   state.m_needfij=false;
   state.m_xupdated=false;
  }
//+------------------------------------------------------------------+
//| Increases lambda,returns False when there is a danger of         |
//| overflow                                                         |
//+------------------------------------------------------------------+
static bool CNlEq::IncreaseLambda(double &lambdav,double &nu,const double lambdaup)
  {
//--- create variables
   double lnlambda=0;
   double lnnu=0;
   double lnlambdaup=0;
   double lnmax=0;
//--- initialization
   lnlambda=MathLog(lambdav);
   lnlambdaup=MathLog(lambdaup);
   lnnu=MathLog(nu);
   lnmax=0.5*MathLog(CMath::m_maxrealnumber);
//--- check
   if(lnlambda+lnlambdaup+lnnu>lnmax)
      return(false);
//--- check
   if(lnnu+MathLog(2)>lnmax)
      return(false);
//--- change values
   lambdav=lambdav*lambdaup*nu;
   nu=nu*2;
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Decreases lambda, but leaves it unchanged when there is danger of|
//| underflow.                                                       |
//+------------------------------------------------------------------+
static void CNlEq::DecreaseLambda(double &lambdav,double &nu,const double lambdadown)
  {
//--- initialization
   nu=1;
//--- check
   if(MathLog(lambdav)+MathLog(lambdadown)<MathLog(CMath::m_minrealnumber))
      lambdav=CMath::m_minrealnumber;
   else
      lambdav=lambdav*lambdadown;
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
static bool CNlEq::NlEqIteration(CNlEqState &state)
  {
//--- create variables
   int    n=0;
   int    m=0;
   int    i=0;
   double lambdaup=0;
   double lambdadown=0;
   double lambdav=0;
   double rho=0;
   double mu=0;
   double stepnorm=0;
   bool   b;
   int    i_=0;
   int    temp;
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
      b=state.m_rstate.ba[0];
      lambdaup=state.m_rstate.ra[0];
      lambdadown=state.m_rstate.ra[1];
      lambdav=state.m_rstate.ra[2];
      rho=state.m_rstate.ra[3];
      mu=state.m_rstate.ra[4];
      stepnorm=state.m_rstate.ra[5];
     }
   else
     {
      //--- initialization
      n=-983;
      m=-989;
      i=-834;
      b=false;
      lambdaup=-287;
      lambdadown=364;
      lambdav=214;
      rho=-338;
      mu=-686;
      stepnorm=912;
     }
//--- check
   if(state.m_rstate.stage==0)
     {
      //--- change values
      state.m_needf=false;
      state.m_repnfunc=state.m_repnfunc+1;
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         state.m_xbase[i_]=state.m_x[i_];
      //--- change values
      state.m_fbase=state.m_f;
      state.m_fprev=CMath::m_maxrealnumber;
      //--- check
      if(!state.m_xrep)
        {
         //--- check
         if(!Func_lbl_5(state,lambdaup,lambdadown,lambdav,rho))
            return(false);
         //--- function call
         Func_lbl_7(state,n);
         //--- Saving state
         Func_lbl_rcomm(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
         //--- return result
         return(true);
        }
      //--- progress report
      ClearRequestFields(state);
      state.m_xupdated=true;
      state.m_rstate.stage=1;
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==1)
     {
      //--- change value
      state.m_xupdated=false;
      //--- check
      if(!Func_lbl_5(state,lambdaup,lambdadown,lambdav,rho))
         return(false);
      //--- function call
      Func_lbl_7(state,n);
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==2)
     {
      //--- change values
      state.m_needfij=false;
      state.m_repnfunc=state.m_repnfunc+1;
      state.m_repnjac=state.m_repnjac+1;
      //--- function call
      CAblas::RMatrixMVect(n,m,state.m_j,0,0,1,state.m_fi,0,state.m_rightpart,0);
      for(i_=0;i_<=n-1;i_++)
         state.m_rightpart[i_]=-1*state.m_rightpart[i_];
      //--- Inner cycle: find good lambda
      temp=Func_lbl_9(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
      //--- check
      if(temp==-1)
         return(false);
      //--- check
      if(temp==1)
         return(true);
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==3)
     {
      //--- change values
      state.m_needf=false;
      state.m_repnfunc=state.m_repnfunc+1;
      //--- check
      if(state.m_f<state.m_fbase)
        {
         //--- function value decreased,move on
         DecreaseLambda(lambdav,rho,lambdadown);
         temp=Func_lbl_10(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
         //--- check
         if(temp==-1)
            return(false);
         //--- check
         if(temp==1)
            return(true);
         //--- Saving state
         Func_lbl_rcomm(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
         //--- return result
         return(true);
        }
      //--- check
      if(!IncreaseLambda(lambdav,rho,lambdaup))
        {
         //--- Lambda is too large (near overflow),force zero step and break
         stepnorm=0;
         for(i_=0;i_<=n-1;i_++)
            state.m_x[i_]=state.m_xbase[i_];
         state.m_f=state.m_fbase;
         temp=Func_lbl_10(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
         //--- check
         if(temp==-1)
            return(false);
         //--- check
         if(temp==1)
            return(true);
         //--- Saving state
         Func_lbl_rcomm(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
         //--- return result
         return(true);
        }
      temp=Func_lbl_9(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
      //--- check
      if(temp==-1)
         return(false);
      //--- check
      if(temp==1)
         return(true);
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
      //--- return result
      return(true);
     }
//--- check
   if(state.m_rstate.stage==4)
     {
      //--- change value
      state.m_xupdated=false;
      //--- check
      if(!Func_lbl_11(state,stepnorm))
         return(false);
      //--- Now,iteration is finally over
      Func_lbl_7(state,n);
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
      //--- return result
      return(true);
     }
//--- Routine body
//--- Prepare
   n=state.m_n;
   m=state.m_m;
   state.m_repterminationtype=0;
   state.m_repiterationscount=0;
   state.m_repnfunc=0;
   state.m_repnjac=0;
//--- Calculate F/G,initialize algorithm
   ClearRequestFields(state);
   state.m_needf=true;
   state.m_rstate.stage=0;
//--- Saving state
   Func_lbl_rcomm(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for NlEqiteration. Is a product to get rid of |
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static void CNlEq::Func_lbl_rcomm(CNlEqState &state,const int n,const int m,
                                  const int i,const bool b,const double lambdaup,
                                  const double lambdadown,const double lambdav,
                                  const double rho,const double mu,const double stepnorm)
  {
//--- save
   state.m_rstate.ia[0]=n;
   state.m_rstate.ia[1]=m;
   state.m_rstate.ia[2]=i;
   state.m_rstate.ba[0]=b;
   state.m_rstate.ra[0]=lambdaup;
   state.m_rstate.ra[1]=lambdadown;
   state.m_rstate.ra[2]=lambdav;
   state.m_rstate.ra[3]=rho;
   state.m_rstate.ra[4]=mu;
   state.m_rstate.ra[5]=stepnorm;
  }
//+------------------------------------------------------------------+
//| Auxiliary function for NlEqiteration. Is a product to get rid of |
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static void CNlEq::Func_lbl_7(CNlEqState &state,const int n)
  {
//--- Get Jacobian;
//--- before we get to this point we already have State.XBase filled
//--- with current point and State.FBase filled with function value
//--- at XBase
   ClearRequestFields(state);
   state.m_needfij=true;
//--- copy
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
   state.m_rstate.stage=2;
  }
//+------------------------------------------------------------------+
//| Auxiliary function for NlEqiteration. Is a product to get rid of |
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CNlEq::Func_lbl_5(CNlEqState &state,double &lambdaup,
                              double &lambdadown,double &lambdav,
                              double &rho)
  {
//--- check
   if(state.m_f<=CMath::Sqr(state.m_epsf))
     {
      state.m_repterminationtype=1;
      return(false);
     }
//--- change values
   lambdaup=10;
   lambdadown=0.3;
   lambdav=0.001;
   rho=1;
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for NlEqiteration. Is a product to get rid of |
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CNlEq::Func_lbl_11(CNlEqState &state,const double stepnorm)
  {
//--- Test stopping conditions on F,step (zero/non-zero) and MaxIts;
//--- If one of the conditions is met,RepTerminationType is changed.
   if(MathSqrt(state.m_f)<=state.m_epsf)
      state.m_repterminationtype=1;
//--- check
   if(stepnorm==0.0 && state.m_repterminationtype==0)
      state.m_repterminationtype=-4;
//--- check
   if(state.m_repiterationscount>=state.m_maxits && state.m_maxits>0)
      state.m_repterminationtype=5;
//--- check
   if(state.m_repterminationtype!=0)
      return(false);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for NlEqiteration. Is a product to get rid of |
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static int CNlEq::Func_lbl_10(CNlEqState &state,const int n,const int m,
                              const int i,const bool b,const double lambdaup,
                              const double lambdadown,const double lambdav,
                              const double rho,const double mu,
                              const double stepnorm)
  {
//--- Accept step:
//--- * new position
//--- * new function value
   state.m_fbase=state.m_f;
   for(int i_=0;i_<=n-1;i_++)
      state.m_xbase[i_]=state.m_xbase[i_]+stepnorm*state.m_candstep[i_];
   state.m_repiterationscount=state.m_repiterationscount+1;
//--- Report new iteration
   if(!state.m_xrep)
     {
      //--- check
      if(!Func_lbl_11(state,stepnorm))
         return(-1);
      //--- Now,iteration is finally over
      Func_lbl_7(state,n);
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
      //--- return result
      return(1);
     }
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_f=state.m_fbase;
//--- copy
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
   state.m_rstate.stage=4;
//--- return result
   return(0);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for NlEqiteration. Is a product to get rid of |
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static int CNlEq::Func_lbl_9(CNlEqState &state,int &n,int &m,int &i,bool &b,
                             const double lambdaup,const double lambdadown,
                             double &lambdav,const double rho,const double mu,
                             double &stepnorm)
  {
//--- Solve (J^T*J + (Lambda+Mu)*I)*y=J^T*F
//--- to get step d=-y where:
//--- * Mu=||F|| - is damping parameter for nonlinear system
//--- * Lambda   - is additional Levenberg-Marquardt parameter
//---              for better convergence when far away from minimum
   for(i=0;i<=n-1;i++)
      state.m_candstep[i]=0;
//--- function call
   CFbls::FblsSolveCGx(state.m_j,m,n,lambdav,state.m_rightpart,state.m_candstep,state.m_cgbuf);
//--- Normalize step (it must be no more than StpMax)
   stepnorm=0;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(state.m_candstep[i]!=0.0)
        {
         stepnorm=1;
         break;
        }
     }
   CLinMin::LinMinNormalized(state.m_candstep,stepnorm,n);
//--- check
   if(state.m_stpmax!=0.0)
      stepnorm=MathMin(stepnorm,state.m_stpmax);
//--- Test new step - is it good enough?
//--- * if not,Lambda is increased and we try again.
//--- * if step is good,we decrease Lambda and move on.
//--- We can break this cycle on two occasions:
//--- * step is so small that x+step==x (in floating point arithmetics)
//--- * lambda is so large
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_x[i_]+stepnorm*state.m_candstep[i_];
   b=true;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(state.m_x[i]!=state.m_xbase[i])
        {
         b=false;
         break;
        }
     }
//--- check
   if(b)
     {
      //--- Step is too small,force zero step and break
      stepnorm=0;
      for(int i_=0;i_<=n-1;i_++)
         state.m_x[i_]=state.m_xbase[i_];
      state.m_f=state.m_fbase;
      //--- function call
      int temp=Func_lbl_10(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
      //--- check
      if(temp!=0)
         return(temp);
      //--- Saving state
      Func_lbl_rcomm(state,n,m,i,b,lambdaup,lambdadown,lambdav,rho,mu,stepnorm);
      //--- return result
      return(1);
     }
//--- function call
   ClearRequestFields(state);
//--- change values
   state.m_needf=true;
   state.m_rstate.stage=3;
//--- return result
   return(0);
  }
//+------------------------------------------------------------------+
