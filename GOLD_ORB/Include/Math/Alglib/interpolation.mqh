//+------------------------------------------------------------------+
//|                                                interpolation.mqh |
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
#include "alglibmisc.mqh"
#include "optimization.mqh"
#include "solvers.mqh"
#include "integration.mqh"
//+------------------------------------------------------------------+
//| IDW interpolant.                                                 |
//+------------------------------------------------------------------+
class CIDWInterpolant
  {
public:
   int               m_n;
   int               m_nx;
   int               m_d;
   double            m_r;
   int               m_nw;
   CKDTree           m_tree;
   int               m_modeltype;
   int               m_debugsolverfailures;
   double            m_debugworstrcond;
   double            m_debugbestrcond;
   //--- arrays
   double            m_xbuf[];
   int               m_tbuf[];
   double            m_rbuf[];
   //--- matrices
   CMatrixDouble     m_q;
   CMatrixDouble     m_xybuf;
   
public:   
                     CIDWInterpolant(void);
                    ~CIDWInterpolant(void);

   void              Copy(CIDWInterpolant &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CIDWInterpolant::CIDWInterpolant(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CIDWInterpolant::~CIDWInterpolant(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CIDWInterpolant::Copy(CIDWInterpolant &obj)
  {
//--- copy variables
   m_n=obj.m_n;
   m_nx=obj.m_nx;
   m_d=obj.m_d;
   m_r=obj.m_r;
   m_nw=obj.m_nw;
   m_modeltype=obj.m_modeltype;
   m_debugsolverfailures=obj.m_debugsolverfailures;
   m_debugworstrcond=obj.m_debugworstrcond;
   m_debugbestrcond=obj.m_debugbestrcond;
   m_tree.Copy(obj.m_tree);
//--- copy arrays
   ArrayCopy(m_xbuf,obj.m_xbuf);
   ArrayCopy(m_tbuf,obj.m_tbuf);
   ArrayCopy(m_rbuf,obj.m_rbuf);
//--- copy matrices
   m_q=obj.m_q;
   m_xybuf=obj.m_xybuf;
  }
//+------------------------------------------------------------------+
//| IDW interpolant.                                                 |
//+------------------------------------------------------------------+
class CIDWInterpolantShell
  {
private:
   CIDWInterpolant   m_innerobj;
public:
   //--- constructors, destructor
                     CIDWInterpolantShell(void);
                     CIDWInterpolantShell(CIDWInterpolant &obj);
                    ~CIDWInterpolantShell(void);
   //--- method
   CIDWInterpolant *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CIDWInterpolantShell::CIDWInterpolantShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CIDWInterpolantShell::CIDWInterpolantShell(CIDWInterpolant &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CIDWInterpolantShell::~CIDWInterpolantShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CIDWInterpolant *CIDWInterpolantShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Inverse distance weighting interpolation                         |
//+------------------------------------------------------------------+
class CIDWInt
  {
private:
   //--- private methods
   static double     IDWCalcQ(CIDWInterpolant &z,double &x[],const int k);
   static void       IDWInit1(const int n,const int nx,const int d,int nq,int nw,CIDWInterpolant &z);
   static void       IDWInternalSolver(double &y[],double &w[],CMatrixDouble &fmatrix,double &temp[],const int n,const int m,int &info,double &x[],double &taskrcond);
public:
   //--- class constants
   static const double m_idwqfactor;
   static const int  m_idwkmin;
   //--- constructor, destructor
                     CIDWInt(void);
                    ~CIDWInt(void);
   //--- public methods
   static double     IDWCalc(CIDWInterpolant &z,double &x[]);
   static void       IDWBuildModifiedShepard(CMatrixDouble &xy,const int n,const int nx,const int d,int nq,int nw,CIDWInterpolant &z);
   static void       IDWBuildModifiedShepardR(CMatrixDouble &xy,const int n,const int nx,const double r,CIDWInterpolant &z);
   static void       IDWBuildNoisy(CMatrixDouble &xy,const int n,const int nx,const int d,int nq,int nw,CIDWInterpolant &z);
  };
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const double CIDWInt::m_idwqfactor=1.5;
const int    CIDWInt::m_idwkmin=5;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CIDWInt::CIDWInt(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CIDWInt::~CIDWInt(void)
  {

  }
//+------------------------------------------------------------------+
//| IDW interpolation                                                |
//| INPUT PARAMETERS:                                                |
//|     Z   -   IDW interpolant built with one of model building     |
//|             subroutines.                                         |
//|     X   -   array[0..NX-1], interpolation point                  |
//| Result:                                                          |
//|     IDW interpolant Z(X)                                         |
//+------------------------------------------------------------------+
static double CIDWInt::IDWCalc(CIDWInterpolant &z,double &x[])
  {
//--- create variables
   double result=0;
   int    nx=0;
   int    i=0;
   int    k=0;
   double r=0;
   double s=0;
   double w=0;
   double v1=0;
   double v2=0;
   double d0=0;
   double di=0;
//--- these initializers are not really necessary,
//--- but without them compiler complains about uninitialized locals
   k=0;
//--- Query
   if(z.m_modeltype==0)
     {
      //--- NQ/NW-based model
      nx=z.m_nx;
      k=CNearestNeighbor::KDTreeQueryKNN(z.m_tree,x,z.m_nw,true);
      //--- function call
      CNearestNeighbor::KDTreeQueryResultsDistances(z.m_tree,z.m_rbuf);
      //--- function call
      CNearestNeighbor::KDTreeQueryResultsTags(z.m_tree,z.m_tbuf);
     }
   if(z.m_modeltype==1)
     {
      //--- R-based model
      nx=z.m_nx;
      k=CNearestNeighbor::KDTreeQueryRNN(z.m_tree,x,z.m_r,true);
      //--- function call
      CNearestNeighbor::KDTreeQueryResultsDistances(z.m_tree,z.m_rbuf);
      //--- function call
      CNearestNeighbor::KDTreeQueryResultsTags(z.m_tree,z.m_tbuf);
      if(k<m_idwkmin)
        {
         //--- we need at least IDWKMin points
         k=CNearestNeighbor::KDTreeQueryKNN(z.m_tree,x,m_idwkmin,true);
         //--- function call
         CNearestNeighbor::KDTreeQueryResultsDistances(z.m_tree,z.m_rbuf);
         //--- function call
         CNearestNeighbor::KDTreeQueryResultsTags(z.m_tree,z.m_tbuf);
        }
     }
//--- initialize weights for linear/quadratic members calculation.
//--- NOTE 1: weights are calculated using NORMALIZED modified
//--- Shepard's formula. Original formula gives w(i)=sqr((R-di)/(R*di)),
//--- where di is i-th distance,R is max(di). Modified formula have
//--- following form:
//---     w_mod(i)=1,if di=d0
//---     w_mod(i)=w(i)/w(0),if di<>d0
//--- NOTE 2: self-match is USED for this query
//--- NOTE 3: last point almost always gain zero weight,but it MUST
//--- be used for fitting because sometimes it will gain NON-ZERO
//--- weight - for example,when all distances are equal.
   r=z.m_rbuf[k-1];
   d0=z.m_rbuf[0];
   result=0;
   s=0;
   for(i=0;i<=k-1;i++)
     {
      di=z.m_rbuf[i];
      //--- check
      if(di==d0)
        {
         //--- distance is equal to shortest,set it 1.0
         //--- without explicitly calculating (which would give
         //--- us same result,but 'll expose us to the risk of
         //--- division by zero).
         w=1;
        }
      else
        {
         //--- use normalized formula
         v1=(r-di)/(r-d0);
         v2=d0/di;
         w=CMath::Sqr(v1*v2);
        }
      //--- change result
      result=result+w*IDWCalcQ(z,x,z.m_tbuf[i]);
      s=s+w;
     }
//--- return result
   return(result/s);
  }
//+------------------------------------------------------------------+
//| IDW interpolant using modified Shepard method for uniform point  |
//| distributions.                                                   |
//| INPUT PARAMETERS:                                                |
//|     XY  -   X and Y values, array[0..N-1,0..NX].                 |
//|             First NX columns contain X-values, last column       |
//|             contain Y-values.                                    |
//|     N   -   number of nodes, N>0.                                |
//|     NX  -   space dimension, NX>=1.                              |
//|     D   -   nodal function type, either:                         |
//|             * 0     constant  model. Just for demonstration only,|
//|                     worst model ever.                            |
//|             * 1     linear model, least squares fitting. Simpe   |
//|                     model for datasets too small for quadratic   |
//|                     models                                       |
//|             * 2     quadratic model, least squares fitting.      |
//|                     Best model available (if your dataset is     |
//|                     large enough).                               |
//|             * -1    "fast" linear model, use with caution!!! It  |
//|                     is significantly faster than linear/quadratic|
//|                     and better than constant model. But it is    |
//|                     less robust (especially in the presence of   |
//|                     noise).                                      |
//|     NQ  -   number of points used to calculate  nodal  functions |
//|             (ignored for constant models). NQ should be LARGER   |
//|             than:                                                |
//|             * max(1.5*(1+NX),2^NX+1) for linear model,           |
//|             * max(3/4*(NX+2)*(NX+1),2^NX+1) for quadratic model. |
//|             Values less than this threshold will be silently     |
//|             increased.                                           |
//|     NW  -   number of points used to calculate weights and to    |
//|             interpolate. Required: >=2^NX+1, values less than    |
//|             this  threshold will be silently increased.          |
//|             Recommended value: about 2*NQ                        |
//| OUTPUT PARAMETERS:                                               |
//|     Z   -   IDW interpolant.                                     |
//| NOTES:                                                           |
//|   * best results are obtained with quadratic models, worst - with|
//|     constant models                                              |
//|   * when N is large, NQ and NW must be significantly smaller than|
//|     N  both to obtain optimal performance and to obtain optimal  |
//|     accuracy. In 2 or 3-dimensional tasks NQ=15 and NW=25 are    |
//|     good values to start with.                                   |
//|   * NQ and NW may be greater than N. In such cases they will be  |
//|     automatically decreased.                                     |
//|   * this subroutine is always succeeds (as long as correct       |
//|     parameters are passed).                                      |
//|   * see 'Multivariate Interpolation of Large Sets of Scattered   |
//|     Data' by Robert J. Renka for more information on this        |
//|     algorithm.                                                   |
//|   * this subroutine assumes that point distribution is uniform at|
//|     the small scales. If it isn't  -  for example, points are    |
//|     concentrated along "lines", but "lines" distribution is      |
//|     uniform at the larger scale - then you should use            |
//|     IDWBuildModifiedShepardR()                                   |
//+------------------------------------------------------------------+
static void CIDWInt::IDWBuildModifiedShepard(CMatrixDouble &xy,const int n,
                                             const int nx,const int d,
                                             int nq,int nw,CIDWInterpolant &z)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    j2=0;
   int    j3=0;
   double v=0;
   double r=0;
   double s=0;
   double d0=0;
   double di=0;
   double v1=0;
   double v2=0;
   int    nc=0;
   int    offs=0;
   int    info=0;
   double taskrcond=0;
   int    i_=0;
//--- create arrays
   double x[];
   double qrbuf[];
   double y[];
   double w[];
   double qsol[];
   double temp[];
   int    tags[];
//--- create matrix
   CMatrixDouble qxybuf;
   CMatrixDouble fmatrix;
//--- these initializers are not really necessary,
//--- but without them compiler complains about uninitialized locals
   nc=0;
//--- assertions
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(nx>=1,__FUNCTION__+": NX<1!"))
      return;
//--- check
   if(!CAp::Assert(d>=-1 && d<=2,__FUNCTION__+": D<>-1 and D<>0 and D<>1 and D<>2!"))
      return;
//--- Correct parameters if needed
   if(d==1)
     {
      nq=MathMax(nq,(int)MathCeil(m_idwqfactor*(1+nx))+1);
      nq=MathMax(nq,(int)MathRound(MathPow(2,nx))+1);
     }
//--- check
   if(d==2)
     {
      nq=MathMax(nq,(int)MathCeil(m_idwqfactor*(nx+2)*(nx+1)/2)+1);
      nq=MathMax(nq,(int)MathRound(MathPow(2,nx))+1);
     }
//--- change values
   nw=MathMax(nw,(int)MathRound(MathPow(2,nx))+1);
   nq=MathMin(nq,n);
   nw=MathMin(nw,n);
//--- primary initialization of Z
   IDWInit1(n,nx,d,nq,nw,z);
   z.m_modeltype=0;
//--- Create KD-tree
   ArrayResizeAL(tags,n);
   for(i=0;i<=n-1;i++)
      tags[i]=i;
//--- function call
   CNearestNeighbor::KDTreeBuildTagged(xy,tags,n,nx,1,2,z.m_tree);
//--- build nodal functions
   ArrayResizeAL(temp,nq+1);
   ArrayResizeAL(x,nx);
   ArrayResizeAL(qrbuf,nq);
   qxybuf.Resize(nq,nx+1);
//--- check
   if(d==-1)
      ArrayResizeAL(w,nq);
//--- check
   if(d==1)
     {
      //--- allocation
      ArrayResizeAL(y,nq);
      ArrayResizeAL(w,nq);
      ArrayResizeAL(qsol,nx);
      //--- NX for linear members,
      //--- 1 for temporary storage
      fmatrix.Resize(nq,nx+1);
     }
//--- check
   if(d==2)
     {
      //--- allocation
      ArrayResizeAL(y,nq);
      ArrayResizeAL(w,nq);
      ArrayResizeAL(qsol,nx+(int)MathRound(nx*(nx+1)*0.5));
      //--- NX for linear members,
      //--- Round(NX*(NX+1)*0.5) for quadratic model,
      //--- 1 for temporary storage
      fmatrix.Resize(nq,nx+(int)MathRound(nx*(nx+1)*0.5)+1);
     }
   for(i=0;i<=n-1;i++)
     {
      //--- Initialize center and function value.
      //--- If D=0 it is all what we need
      for(i_=0;i_<=nx;i_++)
         z.m_q[i].Set(i_,xy[i][i_]);
      //--- check
      if(d==0)
         continue;
      //--- calculate weights for linear/quadratic members calculation.
      //--- NOTE 1: weights are calculated using NORMALIZED modified
      //--- Shepard's formula. Original formula is w(i)=sqr((R-di)/(R*di)),
      //--- where di is i-th distance,R is max(di). Modified formula have
      //--- following form:
      //---     w_mod(i)=1,if di=d0
      //---     w_mod(i)=w(i)/w(0),if di<>d0
      //--- NOTE 2: self-match is NOT used for this query
      //--- NOTE 3: last point almost always gain zero weight,but it MUST
      //--- be used for fitting because sometimes it will gain NON-ZERO
      //--- weight - for example,when all distances are equal.
      for(i_=0;i_<=nx-1;i_++)
         x[i_]=xy[i][i_];
      k=CNearestNeighbor::KDTreeQueryKNN(z.m_tree,x,nq,false);
      //--- function call
      CNearestNeighbor::KDTreeQueryResultsXY(z.m_tree,qxybuf);
      //--- function call
      CNearestNeighbor::KDTreeQueryResultsDistances(z.m_tree,qrbuf);
      r=qrbuf[k-1];
      d0=qrbuf[0];
      //--- calculation
      for(j=0;j<=k-1;j++)
        {
         di=qrbuf[j];
         //--- check
         if(di==d0)
           {
            //--- distance is equal to shortest,set it 1.0
            //--- without explicitly calculating (which would give
            //--- us same result,but 'll expose us to the risk of
            //--- division by zero).
            w[j]=1;
           }
         else
           {
            //--- use normalized formula
            v1=(r-di)/(r-d0);
            v2=d0/di;
            w[j]=CMath::Sqr(v1*v2);
           }
        }
      //--- calculate linear/quadratic members
      if(d==-1)
        {
         //--- "Fast" linear nodal function calculated using
         //--- inverse distance weighting
         for(j=0;j<=nx-1;j++)
            x[j]=0;
         s=0;
         //--- calculation
         for(j=0;j<=k-1;j++)
           {
            //--- calculate J-th inverse distance weighted gradient:
            //---     grad_k=(y_j-y_k)*(x_j-x_k)/sqr(norm(x_j-x_k))
            //---     grad=sum(wk*grad_k)/sum(w_k)
            v=0;
            for(j2=0;j2<=nx-1;j2++)
               v=v+CMath::Sqr(qxybuf[j][j2]-xy[i][j2]);
            //--- Although x_j<>x_k,sqr(norm(x_j-x_k)) may be zero due to
            //--- underflow. If it is,we assume than J-th gradient is zero
            //--- (i.m_e. don't add anything)
            if(v!=0.0)
              {
               for(j2=0;j2<=nx-1;j2++)
                  x[j2]=x[j2]+w[j]*(qxybuf[j][nx]-xy[i][nx])*(qxybuf[j][j2]-xy[i][j2])/v;
              }
            s=s+w[j];
           }
         for(j=0;j<=nx-1;j++)
            z.m_q[i].Set(nx+1+j,x[j]/s);
        }
      else
        {
         //--- Least squares models: build
         if(d==1)
           {
            //--- Linear nodal function calculated using
            //--- least squares fitting to its neighbors
            for(j=0;j<=k-1;j++)
              {
               for(j2=0;j2<=nx-1;j2++)
                  fmatrix[j].Set(j2,qxybuf[j][j2]-xy[i][j2]);
               y[j]=qxybuf[j][nx]-xy[i][nx];
              }
            nc=nx;
           }
         //--- check
         if(d==2)
           {
            //--- Quadratic nodal function calculated using
            //--- least squares fitting to its neighbors
            for(j=0;j<=k-1;j++)
              {
               offs=0;
               for(j2=0;j2<=nx-1;j2++)
                 {
                  fmatrix[j].Set(offs,qxybuf[j][j2]-xy[i][j2]);
                  offs=offs+1;
                 }
               //--- calculation
               for(j2=0;j2<=nx-1;j2++)
                 {
                  for(j3=j2;j3<=nx-1;j3++)
                    {
                     fmatrix[j].Set(offs,(qxybuf[j][j2]-xy[i][j2])*(qxybuf[j][j3]-xy[i][j3]));
                     offs=offs+1;
                    }
                 }
               y[j]=qxybuf[j][nx]-xy[i][nx];
              }
            nc=nx+(int)MathRound(nx*(nx+1)*0.5);
           }
         //--- function call
         IDWInternalSolver(y,w,fmatrix,temp,k,nc,info,qsol,taskrcond);
         //--- Least squares models: copy results
         if(info>0)
           {
            //--- LLS task is solved,copy results
            z.m_debugworstrcond=MathMin(z.m_debugworstrcond,taskrcond);
            z.m_debugbestrcond=MathMax(z.m_debugbestrcond,taskrcond);
            for(j=0;j<=nc-1;j++)
               z.m_q[i].Set(nx+1+j,qsol[j]);
           }
         else
           {
            //--- Solver failure,very strange,but we will use
            //--- zero values to handle it.
            z.m_debugsolverfailures=z.m_debugsolverfailures+1;
            for(j=0;j<=nc-1;j++)
               z.m_q[i].Set(nx+1+j,0);
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| IDW interpolant using modified Shepard method for non-uniform    |
//| datasets.                                                        |
//| This type of model uses constant nodal functions and interpolates|
//| using all nodes which are closer than user-specified radius R. It|
//| may be used when points distribution is non-uniform at the small |
//| scale, but it is at the distances as large as R.                 |
//| INPUT PARAMETERS:                                                |
//|     XY  -   X and Y values, array[0..N-1,0..NX].                 |
//|             First NX columns contain X-values, last column       |
//|             contain Y-values.                                    |
//|     N   -   number of nodes, N>0.                                |
//|     NX  -   space dimension, NX>=1.                              |
//|     R   -   radius, R>0                                          |
//| OUTPUT PARAMETERS:                                               |
//|     Z   -   IDW interpolant.                                     |
//| NOTES:                                                           |
//| * if there is less than IDWKMin points within R-ball, algorithm  |
//|   selects IDWKMin closest ones, so that continuity properties of |
//|   interpolant are preserved even far from points.                |
//+------------------------------------------------------------------+
static void CIDWInt::IDWBuildModifiedShepardR(CMatrixDouble &xy,const int n,
                                              const int nx,const double r,
                                              CIDWInterpolant &z)
  {
//--- create variables
   int i=0;
   int i_=0;
//--- create array
   int tags[];
//--- assertions
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(nx>=1,__FUNCTION__+": NX<1!"))
      return;
//--- check
   if(!CAp::Assert(r>0.0,__FUNCTION__+": R<=0!"))
      return;
//--- primary initialization of Z
   IDWInit1(n,nx,0,0,n,z);
   z.m_modeltype=1;
   z.m_r=r;
//--- Create KD-tree
   ArrayResizeAL(tags,n);
   for(i=0;i<=n-1;i++)
      tags[i]=i;
//--- function call
   CNearestNeighbor::KDTreeBuildTagged(xy,tags,n,nx,1,2,z.m_tree);
//--- build nodal functions
   for(i=0;i<=n-1;i++)
     {
      for(i_=0;i_<=nx;i_++)
         z.m_q[i].Set(i_,xy[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| IDW model for noisy data.                                        |
//| This subroutine may be used to handle noisy data, i.e. data with |
//| noise in OUTPUT values. It differs from IDWBuildModifiedShepard()|
//| in the following aspects:                                        |
//| * nodal functions are not constrained to pass through nodes:     |
//|   Qi(xi)<>yi, i.e. we have fitting instead of interpolation.     |
//| * weights which are used during least squares fitting stage are  |
//|   all equal to 1.0 (independently of distance)                   |
//| * "fast"-linear or constant nodal functions are not supported    |
//|   (either not robust enough or too rigid)                        |
//| This problem require far more complex tuning than interpolation  |
//| problems.                                                        |
//| Below you can find some recommendations regarding this problem:  |
//| * focus on tuning NQ; it controls noise reduction. As for NW, you|
//|   can just make it equal to 2*NQ.                                |
//| * you can use cross-validation to determine optimal NQ.          |
//| * optimal NQ is a result of complex tradeoff between noise level |
//|   (more noise = larger NQ required) and underlying function      |
//|   complexity (given fixed N, larger NQ means smoothing of compex |
//|   features in the data). For example, NQ=N will reduce noise to  |
//|   the minimum level possible, but you will end up with just      |
//|   constant/linear/quadratic (depending on D) least squares       |
//|   model for the whole dataset.                                   |
//| INPUT PARAMETERS:                                                |
//|     XY  -   X and Y values, array[0..N-1,0..NX].                 |
//|             First NX columns contain X-values, last column       |
//|             contain Y-values.                                    |
//|     N   -   number of nodes, N>0.                                |
//|     NX  -   space dimension, NX>=1.                              |
//|     D   -   nodal function degree, either:                       |
//|             * 1     linear model, least squares fitting. Simpe   |
//|                     model for datasets too small for quadratic   |
//|                     models (or for very noisy problems).         |
//|             * 2     quadratic model, least squares fitting. Best |
//|                     model available (if your dataset is large    |
//|                     enough).                                     |
//|     NQ  -   number of points used to calculate nodal functions.  |
//|             NQ should be significantly larger than 1.5 times the |
//|             number of coefficients in a nodal function to        |
//|             overcome effects of noise:                           |
//|             * larger than 1.5*(1+NX) for linear model,           |
//|             * larger than 3/4*(NX+2)*(NX+1) for quadratic model. |
//|             Values less than this threshold will be silently     |
//|             increased.                                           |
//|     NW  -   number of points used to calculate weights and to    |
//|             interpolate. Required: >=2^NX+1, values less than    |
//|             this threshold will be silently increased.           |
//|             Recommended value: about 2*NQ or larger              |
//| OUTPUT PARAMETERS:                                               |
//|     Z   -   IDW inte rpolant.                                    |
//| NOTES:                                                           |
//|   * best results are obtained with quadratic models, linear      |
//|     models are not recommended to use unless you are pretty sure |
//|     that it is what you want                                     |
//|   * this subroutine is always succeeds (as long as correct       |
//|     parameters are passed).                                      |
//|   * see 'Multivariate Interpolation of Large Sets of Scattered   |
//|     Data' by Robert J. Renka for more information on this        |
//|     algorithm.                                                   |
//+------------------------------------------------------------------+
static void CIDWInt::IDWBuildNoisy(CMatrixDouble &xy,const int n,const int nx,
                                   const int d,int nq,int nw,CIDWInterpolant &z)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    j2=0;
   int    j3=0;
   double v=0;
   int    nc=0;
   int    offs=0;
   double taskrcond=0;
   int    info=0;
   int    i_=0;
//--- create arrays
   double x[];
   double qrbuf[];
   double y[];
   double w[];
   double qsol[];
   int    tags[];
   double temp[];
//--- create matrix
   CMatrixDouble qxybuf;
   CMatrixDouble fmatrix;
//--- these initializers are not really necessary,
//--- but without them compiler complains about uninitialized locals
   nc=0;
//--- assertions
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(nx>=1,__FUNCTION__+": NX<1!"))
      return;
//--- check
   if(!CAp::Assert(d>=1 && d<=2,__FUNCTION__+": D<>1 and D<>2!"))
      return;
//--- Correct parameters if needed
   if(d==1)
      nq=MathMax(nq,(int)MathCeil(m_idwqfactor*(1+nx))+1);
//--- check
   if(d==2)
      nq=MathMax(nq,(int)MathCeil(m_idwqfactor*(nx+2)*(nx+1)/2)+1);
//--- change values
   nw=MathMax(nw,(int)MathRound(MathPow(2,nx))+1);
   nq=MathMin(nq,n);
   nw=MathMin(nw,n);
//--- primary initialization of Z
   IDWInit1(n,nx,d,nq,nw,z);
   z.m_modeltype=0;
//--- Create KD-tree
   ArrayResizeAL(tags,n);
   for(i=0;i<=n-1;i++)
      tags[i]=i;
//--- function call
   CNearestNeighbor::KDTreeBuildTagged(xy,tags,n,nx,1,2,z.m_tree);
//--- build nodal functions
//--- (special algorithm for noisy data is used)
   ArrayResizeAL(temp,nq+1);
   ArrayResizeAL(x,nx);
   ArrayResizeAL(qrbuf,nq);
   qxybuf.Resize(nq,nx+1);
//--- check
   if(d==1)
     {
      //--- allocation
      ArrayResizeAL(y,nq);
      ArrayResizeAL(w,nq);
      ArrayResizeAL(qsol,1+nx);
      //--- 1 for constant member,
      //--- NX for linear members,
      //--- 1 for temporary storage
      fmatrix.Resize(nq,1+nx+1);
     }
//--- check
   if(d==2)
     {
      //--- allocation
      ArrayResizeAL(y,nq);
      ArrayResizeAL(w,nq);
      ArrayResizeAL(qsol,1+nx+(int)MathRound(nx*(nx+1)*0.5));
      //--- 1 for constant member,
      //--- NX for linear members,
      //--- Round(NX*(NX+1)*0.5) for quadratic model,
      //--- 1 for temporary storage
      fmatrix.Resize(nq,1+nx+(int)MathRound(nx*(nx+1)*0.5)+1);
     }
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      //--- Initialize center.
      for(i_=0;i_<=nx-1;i_++)
         z.m_q[i].Set(i_,xy[i][i_]);
      //--- Calculate linear/quadratic members
      //--- using least squares fit
      //--- NOTE 1: all weight are equal to 1.0
      //--- NOTE 2: self-match is USED for this query
      for(i_=0;i_<=nx-1;i_++)
         x[i_]=xy[i][i_];
      k=CNearestNeighbor::KDTreeQueryKNN(z.m_tree,x,nq,true);
      //--- function call
      CNearestNeighbor::KDTreeQueryResultsXY(z.m_tree,qxybuf);
      //--- function call
      CNearestNeighbor::KDTreeQueryResultsDistances(z.m_tree,qrbuf);
      //--- check
      if(d==1)
        {
         //--- Linear nodal function calculated using
         //--- least squares fitting to its neighbors
         for(j=0;j<=k-1;j++)
           {
            fmatrix[j].Set(0,1.0);
            for(j2=0;j2<=nx-1;j2++)
               fmatrix[j].Set(1+j2,qxybuf[j][j2]-xy[i][j2]);
            //--- change values
            y[j]=qxybuf[j][nx];
            w[j]=1;
           }
         nc=1+nx;
        }
      //--- check
      if(d==2)
        {
         //--- Quadratic nodal function calculated using
         //--- least squares fitting to its neighbors
         for(j=0;j<=k-1;j++)
           {
            fmatrix[j].Set(0,1);
            offs=1;
            for(j2=0;j2<=nx-1;j2++)
              {
               fmatrix[j].Set(offs,qxybuf[j][j2]-xy[i][j2]);
               offs=offs+1;
              }
            //--- calculation
            for(j2=0;j2<=nx-1;j2++)
              {
               for(j3=j2;j3<=nx-1;j3++)
                 {
                  fmatrix[j].Set(offs,(qxybuf[j][j2]-xy[i][j2])*(qxybuf[j][j3]-xy[i][j3]));
                  offs=offs+1;
                 }
              }
            //--- change values
            y[j]=qxybuf[j][nx];
            w[j]=1;
           }
         nc=1+nx+(int)MathRound(nx*(nx+1)*0.5);
        }
      //--- function call
      IDWInternalSolver(y,w,fmatrix,temp,k,nc,info,qsol,taskrcond);
      //--- Least squares models: copy results
      if(info>0)
        {
         //--- LLS task is solved,copy results
         z.m_debugworstrcond=MathMin(z.m_debugworstrcond,taskrcond);
         z.m_debugbestrcond=MathMax(z.m_debugbestrcond,taskrcond);
         for(j=0;j<=nc-1;j++)
            z.m_q[i].Set(nx+j,qsol[j]);
        }
      else
        {
         //--- Solver failure,very strange,but we will use
         //--- zero values to handle it.
         z.m_debugsolverfailures=z.m_debugsolverfailures+1;
         v=0;
         for(j=0;j<=k-1;j++)
            v=v+qxybuf[j][nx];
         z.m_q[i].Set(nx,v/k);
         for(j=0;j<=nc-2;j++)
            z.m_q[i].Set(nx+1+j,0);
        }
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine: K-th nodal function calculation             |
//+------------------------------------------------------------------+
static double CIDWInt::IDWCalcQ(CIDWInterpolant &z,double &x[],const int k)
  {
//--- create variables
   double result=0;
   int    nx=0;
   int    i=0;
   int    j=0;
   int    offs=0;
//--- initialization
   nx=z.m_nx;
//--- constant member
   result=z.m_q[k][nx];
//--- linear members
   if(z.m_d>=1)
     {
      for(i=0;i<=nx-1;i++)
         result=result+z.m_q[k][nx+1+i]*(x[i]-z.m_q[k][i]);
     }
//--- quadratic members
   if(z.m_d>=2)
     {
      offs=nx+1+nx;
      for(i=0;i<=nx-1;i++)
        {
         for(j=i;j<=nx-1;j++)
           {
            result=result+z.m_q[k][offs]*(x[i]-z.m_q[k][i])*(x[j]-z.m_q[k][j]);
            offs=offs+1;
           }
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Initialization of internal structures.                           |
//| It assumes correctness of all parameters.                        |
//+------------------------------------------------------------------+
static void CIDWInt::IDWInit1(const int n,const int nx,const int d,
                              int nq,int nw,CIDWInterpolant &z)
  {
//--- initialization
   z.m_debugsolverfailures=0;
   z.m_debugworstrcond=1.0;
   z.m_debugbestrcond=0;
   z.m_n=n;
   z.m_nx=nx;
   z.m_d=0;
//--- check
   if(d==1)
      z.m_d=1;
//--- check
   if(d==2)
      z.m_d=2;
//--- check
   if(d==-1)
      z.m_d=1;
   z.m_nw=nw;
//--- check
   if(d==-1)
      z.m_q.Resize(n,2*nx+1);
//--- check
   if(d==0)
      z.m_q.Resize(n,nx+1);
//--- check
   if(d==1)
      z.m_q.Resize(n,2*nx+1);
//--- check
   if(d==2)
      z.m_q.Resize(n,nx+1+nx+(int)MathRound(nx*(nx+1)*0.5));
//--- allocation
   ArrayResizeAL(z.m_tbuf,nw);
   ArrayResizeAL(z.m_rbuf,nw);
   z.m_xybuf.Resize(nw,nx+1);
   ArrayResizeAL(z.m_xbuf,nx);
  }
//+------------------------------------------------------------------+
//| Linear least squares solver for small tasks.                     |
//| Works faster than standard ALGLIB solver in non-degenerate       |
//| cases (due to absense of internal allocations and optimized      |
//| row/colums). In  degenerate cases it calls standard solver, which|
//| results in small performance penalty associated with preliminary |
//| steps.                                                           |
//| INPUT PARAMETERS:                                                |
//|     Y           array[0..N-1]                                    |
//|     W           array[0..N-1]                                    |
//|     FMatrix     array[0..N-1,0..M], have additional column for   |
//|                 temporary values                                 |
//|     Temp        array[0..N]                                      |
//+------------------------------------------------------------------+
static void CIDWInt::IDWInternalSolver(double &y[],double &w[],CMatrixDouble &fmatrix,
                                       double &temp[],const int n,const int m,
                                       int &info,double &x[],double &taskrcond)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   double tau=0;
   int    i_=0;
   int    i1_=0;
//--- create array
   double b[];
//--- object of class
   CDenseSolverLSReport srep;
//--- set up info
   info=1;
//--- prepare matrix
   for(i=0;i<=n-1;i++)
     {
      fmatrix[i].Set(m,y[i]);
      v=w[i];
      for(i_=0;i_<=m;i_++)
         fmatrix[i].Set(i_,v*fmatrix[i][i_]);
     }
//--- use either fast algorithm or general algorithm
   if(m<=n)
     {
      //--- QR decomposition
      //--- We assume that M<=N (we would have called LSFit() otherwise)
      for(i=0;i<=m-1;i++)
        {
         //--- check
         if(i<n-1)
           {
            i1_=i-1;
            for(i_=1;i_<=n-i;i_++)
               temp[i_]=fmatrix[i_+i1_][i];
            //--- function call
            CReflections::GenerateReflection(temp,n-i,tau);
            //--- change values
            fmatrix[i].Set(i,temp[1]);
            temp[1]=1;
            //--- calculation
            for(j=i+1;j<=m;j++)
              {
               i1_=1-i;
               v=0.0;
               for(i_=i;i_<=n-1;i_++)
                  v+=fmatrix[i_][j]*temp[i_+i1_];
               //--- change values
               v=tau*v;
               i1_=1-i;
               for(i_=i;i_<=n-1;i_++)
                  fmatrix[i_].Set(j,fmatrix[i_][j]-v*temp[i_+i1_]);
              }
           }
        }
      //--- Check condition number
      taskrcond=CRCond::RMatrixTrRCondInf(fmatrix,m,true,false);
      //--- use either fast algorithm for non-degenerate cases
      //--- or slow algorithm for degenerate cases
      if(taskrcond>10000*n*CMath::m_machineepsilon)
        {
         //--- solve triangular system R*x=FMatrix[0:M-1,M]
         //--- using fast algorithm,then exit
         x[m-1]=fmatrix[m-1][m]/fmatrix[m-1][m-1];
         for(i=m-2;i>=0;i--)
           {
            v=0.0;
            for(i_=i+1;i_<=m-1;i_++)
               v+=fmatrix[i][i_]*x[i_];
            x[i]=(fmatrix[i][m]-v)/fmatrix[i][i];
           }
        }
      else
        {
         //--- use more general algorithm
         ArrayResizeAL(b,m);
         for(i=0;i<=m-1;i++)
           {
            for(j=0;j<=i-1;j++)
               fmatrix[i].Set(j,0.0);
            b[i]=fmatrix[i][m];
           }
         //--- function call
         CDenseSolver::RMatrixSolveLS(fmatrix,m,m,b,10000*CMath::m_machineepsilon,info,srep,x);
        }
     }
   else
     {
      //--- use more general algorithm
      ArrayResizeAL(b,n);
      for(i=0;i<=n-1;i++)
         b[i]=fmatrix[i][m];
      //--- function call
      CDenseSolver::RMatrixSolveLS(fmatrix,n,m,b,10000*CMath::m_machineepsilon,info,srep,x);
      taskrcond=srep.m_r2;
     }
  }
//+------------------------------------------------------------------+
//| Barycentric interpolant.                                         |
//+------------------------------------------------------------------+
class CBarycentricInterpolant
  {
public:
   //--- variables
   int               m_n;
   double            m_sy;
   //--- arrays
   double            m_x[];
   double            m_y[];
   double            m_w[];
   //--- constructor, destructor
                     CBarycentricInterpolant(void);
                    ~CBarycentricInterpolant(void);
   //--- copy
   void              Copy(CBarycentricInterpolant &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CBarycentricInterpolant::CBarycentricInterpolant(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CBarycentricInterpolant::~CBarycentricInterpolant(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CBarycentricInterpolant::Copy(CBarycentricInterpolant &obj)
  {
//--- copy variables
   m_n=obj.m_n;
   m_sy=obj.m_sy;
//--- copy arrays
   ArrayCopy(m_x,obj.m_x);
   ArrayCopy(m_y,obj.m_y);
   ArrayCopy(m_w,obj.m_w);
  }
//+------------------------------------------------------------------+
//| Barycentric interpolant.                                         |
//+------------------------------------------------------------------+
class CBarycentricInterpolantShell
  {
private:
   CBarycentricInterpolant m_innerobj;
public:
   //--- constructor, destructor
                     CBarycentricInterpolantShell(void);
                     CBarycentricInterpolantShell(CBarycentricInterpolant &obj);
                    ~CBarycentricInterpolantShell(void);
   //--- method
   CBarycentricInterpolant *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CBarycentricInterpolantShell::CBarycentricInterpolantShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CBarycentricInterpolantShell::CBarycentricInterpolantShell(CBarycentricInterpolant &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CBarycentricInterpolantShell::~CBarycentricInterpolantShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CBarycentricInterpolant *CBarycentricInterpolantShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Rational interpolation                                           |
//+------------------------------------------------------------------+
class CRatInt
  {
private:
   //--- private method
   static void       BarycentricNormalize(CBarycentricInterpolant &b);
public:
   //--- constructor, destructor
                     CRatInt(void);
                    ~CRatInt(void);
   //--- public methods
   static double     BarycentricCalc(CBarycentricInterpolant &b,const double t);
   static void       BarycentricDiff1(CBarycentricInterpolant &b,double t,double &f,double &df);
   static void       BarycentricDiff2(CBarycentricInterpolant &b,const double t,double &f,double &df,double &d2f);
   static void       BarycentricLinTransX(CBarycentricInterpolant &b,const double ca,const double cb);
   static void       BarycentricLinTransY(CBarycentricInterpolant &b,const double ca,const double cb);
   static void       BarycentricUnpack(CBarycentricInterpolant &b,int &n,double &x[],double &y[],double &w[]);
   static void       BarycentricBuildXYW(double &x[],double &y[],double &w[],const int n,CBarycentricInterpolant &b);
   static void       BarycentricBuildFloaterHormann(double &x[],double &y[],const int n,int d,CBarycentricInterpolant &b);
   static void       BarycentricCopy(CBarycentricInterpolant &b,CBarycentricInterpolant &b2);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CRatInt::CRatInt(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CRatInt::~CRatInt(void)
  {

  }
//+------------------------------------------------------------------+
//| Rational interpolation using barycentric formula                 |
//| F(t)=SUM(i=0,n-1,w[i]*f[i]/(t-x[i])) / SUM(i=0,n-1,w[i]/(t-x[i]))|
//| Input parameters:                                                |
//|     B   -   barycentric interpolant built with one of model      |
//|             building subroutines.                                |
//|     T   -   interpolation point                                  |
//| Result:                                                          |
//|     barycentric interpolant F(t)                                 |
//+------------------------------------------------------------------+
static double CRatInt::BarycentricCalc(CBarycentricInterpolant &b,const double t)
  {
//--- create variables
   double s1=0;
   double s2=0;
   double s=0;
   double v=0;
   int    i=0;
//--- check
   if(!CAp::Assert(!CInfOrNaN::IsInfinity(t),__FUNCTION__+": infinite T!"))
      return(EMPTY_VALUE);
//--- special case: NaN
   if(CInfOrNaN::IsNaN(t))
      return(CInfOrNaN::NaN());
//--- special case: N=1
   if(b.m_n==1)
      return(b.m_sy*b.m_y[0]);
//--- Here we assume that task is normalized,i.m_e.:
//--- 1. abs(Y[i])<=1
//--- 2. abs(W[i])<=1
//--- 3. X[] is ordered
   s=MathAbs(t-b.m_x[0]);
   for(i=0;i<=b.m_n-1;i++)
     {
      v=b.m_x[i];
      //--- check
      if(v==(double)(t))
         return(b.m_sy*b.m_y[i]);
      v=MathAbs(t-v);
      //--- check
      if(v<s)
         s=v;
     }
//--- change values
   s1=0;
   s2=0;
//--- calculation
   for(i=0;i<=b.m_n-1;i++)
     {
      v=s/(t-b.m_x[i]);
      v=v*b.m_w[i];
      s1=s1+v*b.m_y[i];
      s2=s2+v;
     }
//--- return result
   return(b.m_sy*s1/s2);
  }
//+------------------------------------------------------------------+
//| Differentiation of barycentric interpolant: first derivative.    |
//| Algorithm used in this subroutine is very robust and should not  |
//| fail until provided with values too close to MaxRealNumber       |
//| (usually  MaxRealNumber/N or greater will overflow).             |
//| INPUT PARAMETERS:                                                |
//|     B   -   barycentric interpolant built with one of model      |
//|             building subroutines.                                |
//|     T   -   interpolation point                                  |
//| OUTPUT PARAMETERS:                                               |
//|     F   -   barycentric interpolant at T                         |
//|     DF  -   first derivative                                     |
//+------------------------------------------------------------------+
static void CRatInt::BarycentricDiff1(CBarycentricInterpolant &b,double t,
                                      double &f,double &df)
  {
//--- create variables
   double v=0;
   double vv=0;
   int    i=0;
   int    k=0;
   double n0=0;
   double n1=0;
   double d0=0;
   double d1=0;
   double s0=0;
   double s1=0;
   double xk=0;
   double xi=0;
   double xmin=0;
   double xmax=0;
   double xscale1=0;
   double xoffs1=0;
   double xscale2=0;
   double xoffs2=0;
   double xprev=0;
//--- initialization
   f=0;
   df=0;
//--- check
   if(!CAp::Assert(!CInfOrNaN::IsInfinity(t),__FUNCTION__+": infinite T!"))
      return;
//--- special case: NaN
   if(CInfOrNaN::IsNaN(t))
     {
      //--- change values
      f=CInfOrNaN::NaN();
      df=CInfOrNaN::NaN();
      //--- exit the function
      return;
     }
//--- special case: N=1
   if(b.m_n==1)
     {
      //--- change values
      f=b.m_sy*b.m_y[0];
      df=0;
      //--- exit the function
      return;
     }
//--- check
   if(b.m_sy==0.0)
     {
      //--- change values
      f=0;
      df=0;
      //--- exit the function
      return;
     }
//--- check
   if(!CAp::Assert(b.m_sy>0.0,__FUNCTION__+": internal error"))
      return;
//--- We assume than N>1 and B.SY>0. Find:
//--- 1. pivot point (X[i] closest to T)
//--- 2. width of interval containing X[i]
   v=MathAbs(b.m_x[0]-t);
   k=0;
   xmin=b.m_x[0];
   xmax=b.m_x[0];
//--- calculation
   for(i=1;i<=b.m_n-1;i++)
     {
      vv=b.m_x[i];
      //--- check
      if(MathAbs(vv-t)<v)
        {
         v=MathAbs(vv-t);
         k=i;
        }
      //--- change values
      xmin=MathMin(xmin,vv);
      xmax=MathMax(xmax,vv);
     }
//--- pivot point found,calculate dNumerator and dDenominator
   xscale1=1/(xmax-xmin);
   xoffs1=-(xmin/(xmax-xmin))+1;
   xscale2=2;
   xoffs2=-3;
   t=t*xscale1+xoffs1;
   t=t*xscale2+xoffs2;
   xk=b.m_x[k];
   xk=xk*xscale1+xoffs1;
   xk=xk*xscale2+xoffs2;
   v=t-xk;
   n0=0;
   n1=0;
   d0=0;
   d1=0;
   xprev=-2;
//--- calculation
   for(i=0;i<=b.m_n-1;i++)
     {
      //--- change values
      xi=b.m_x[i];
      xi=xi*xscale1+xoffs1;
      xi=xi*xscale2+xoffs2;
      //--- check
      if(!CAp::Assert(xi>xprev,__FUNCTION__+": points are too close!"))
         return;
      xprev=xi;
      //--- check
      if(i!=k)
        {
         vv=CMath::Sqr(t-xi);
         s0=(t-xk)/(t-xi);
         s1=(xk-xi)/vv;
        }
      else
        {
         s0=1;
         s1=0;
        }
      //--- change values
      vv=b.m_w[i]*b.m_y[i];
      n0=n0+s0*vv;
      n1=n1+s1*vv;
      vv=b.m_w[i];
      d0=d0+s0*vv;
      d1=d1+s1*vv;
     }
//--- change values
   f=b.m_sy*n0/d0;
   df=(n1*d0-n0*d1)/CMath::Sqr(d0);
//--- check
   if(df!=0.0)
      df=MathSign(df)*MathExp(MathLog(MathAbs(df))+MathLog(b.m_sy)+MathLog(xscale1)+MathLog(xscale2));
  }
//+------------------------------------------------------------------+
//| Differentiation of barycentric interpolant: first/second         |
//| derivatives.                                                     |
//| INPUT PARAMETERS:                                                |
//|     B   -   barycentric interpolant built with one of model      |
//|             building subroutines.                                |
//|     T   -   interpolation point                                  |
//| OUTPUT PARAMETERS:                                               |
//|     F   -   barycentric interpolant at T                         |
//|     DF  -   first derivative                                     |
//|     D2F -   second derivative                                    |
//| NOTE: this algorithm may fail due to overflow/underflor if used  |
//| on data whose values are close to MaxRealNumber or MinRealNumber.|
//| Use more robust BarycentricDiff1() subroutine in such cases.     |
//+------------------------------------------------------------------+
static void CRatInt::BarycentricDiff2(CBarycentricInterpolant &b,const double t,
                                      double &f,double &df,double &d2f)
  {
//--- create variables
   double v=0;
   double vv=0;
   int    i=0;
   int    k=0;
   double n0=0;
   double n1=0;
   double n2=0;
   double d0=0;
   double d1=0;
   double d2=0;
   double s0=0;
   double s1=0;
   double s2=0;
   double xk=0;
   double xi=0;
//--- initialization
   f=0;
   df=0;
   d2f=0;
//--- check
   if(!CAp::Assert(!CInfOrNaN::IsInfinity(t),__FUNCTION__+": infinite T!"))
      return;
//--- special case: NaN
   if(CInfOrNaN::IsNaN(t))
     {
      //--- change values
      f=CInfOrNaN::NaN();
      df=CInfOrNaN::NaN();
      d2f=CInfOrNaN::NaN();
      //--- exit the function
      return;
     }
//--- special case: N=1
   if(b.m_n==1)
     {
      //--- change values
      f=b.m_sy*b.m_y[0];
      df=0;
      d2f=0;
      //--- exit the function
      return;
     }
//--- check
   if(b.m_sy==0.0)
     {
      //--- change values
      f=0;
      df=0;
      d2f=0;
      //--- exit the function
      return;
     }
//--- We assume than N>1 and B.SY>0. Find:
//--- 1. pivot point (X[i] closest to T)
//--- 2. width of interval containing X[i]
   if(!CAp::Assert(b.m_sy>0.0,__FUNCTION__+": internal error"))
      return;
//--- change values
   f=0;
   df=0;
   d2f=0;
   v=MathAbs(b.m_x[0]-t);
   k=0;
   for(i=1;i<=b.m_n-1;i++)
     {
      vv=b.m_x[i];
      //--- check
      if(MathAbs(vv-t)<v)
        {
         v=MathAbs(vv-t);
         k=i;
        }
     }
//--- pivot point found, calculate dNumerator and dDenominator
   xk=b.m_x[k];
   v=t-xk;
   n0=0;
   n1=0;
   n2=0;
   d0=0;
   d1=0;
   d2=0;
//--- calculation
   for(i=0;i<=b.m_n-1;i++)
     {
      //--- check
      if(i!=k)
        {
         xi=b.m_x[i];
         vv=CMath::Sqr(t-xi);
         s0=(t-xk)/(t-xi);
         s1=(xk-xi)/vv;
         s2=-(2*(xk-xi)/(vv*(t-xi)));
        }
      else
        {
         s0=1;
         s1=0;
         s2=0;
        }
      //--- change values
      vv=b.m_w[i]*b.m_y[i];
      n0=n0+s0*vv;
      n1=n1+s1*vv;
      n2=n2+s2*vv;
      vv=b.m_w[i];
      d0=d0+s0*vv;
      d1=d1+s1*vv;
      d2=d2+s2*vv;
     }
//--- change values
   f=b.m_sy*n0/d0;
   df=b.m_sy*(n1*d0-n0*d1)/CMath::Sqr(d0);
   d2f=b.m_sy*((n2*d0-n0*d2)*CMath::Sqr(d0)-(n1*d0-n0*d1)*2*d0*d1)/CMath::Sqr(CMath::Sqr(d0));
  }
//+------------------------------------------------------------------+
//| This subroutine performs linear transformation of the argument.  |
//| INPUT PARAMETERS:                                                |
//|     B       -   rational interpolant in barycentric form         |
//|     CA, CB  -   transformation coefficients: x = CA*t + CB       |
//| OUTPUT PARAMETERS:                                               |
//|     B       -   transformed interpolant with X replaced by T     |
//+------------------------------------------------------------------+
static void CRatInt::BarycentricLinTransX(CBarycentricInterpolant &b,
                                          const double ca,const double cb)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
//--- special case,replace by constant F(CB)
   if(ca==0.0)
     {
      b.m_sy=BarycentricCalc(b,cb);
      v=1;
      for(i=0;i<=b.m_n-1;i++)
        {
         b.m_y[i]=1;
         b.m_w[i]=v;
         v=-v;
        }
      //--- exit the function
      return;
     }
//--- general case: CA<>0
   for(i=0;i<=b.m_n-1;i++)
      b.m_x[i]=(b.m_x[i]-cb)/ca;
//--- check
   if(ca<0.0)
     {
      for(i=0;i<=b.m_n-1;i++)
        {
         //--- check
         if(i<b.m_n-1-i)
           {
            //--- change values
            j=b.m_n-1-i;
            v=b.m_x[i];
            b.m_x[i]=b.m_x[j];
            b.m_x[j]=v;
            v=b.m_y[i];
            b.m_y[i]=b.m_y[j];
            b.m_y[j]=v;
            v=b.m_w[i];
            b.m_w[i]=b.m_w[j];
            b.m_w[j]=v;
           }
         else
            break;
        }
     }
  }
//+------------------------------------------------------------------+
//| This subroutine performs linear transformation of the barycentric|
//| interpolant.                                                     |
//| INPUT PARAMETERS:                                                |
//|     B       -   rational interpolant in barycentric form         |
//|     CA, CB  -   transformation coefficients: B2(x) = CA*B(x) + CB|
//| OUTPUT PARAMETERS:                                               |
//|     B       -   transformed interpolant                          |
//+------------------------------------------------------------------+
static void CRatInt::BarycentricLinTransY(CBarycentricInterpolant &b,
                                          const double ca,const double cb)
  {
//--- create variables
   int    i=0;
   double v=0;
   int    i_=0;
//--- calculation
   for(i=0;i<=b.m_n-1;i++)
      b.m_y[i]=ca*b.m_sy*b.m_y[i]+cb;
//--- change value
   b.m_sy=0;
   for(i=0;i<=b.m_n-1;i++)
      b.m_sy=MathMax(b.m_sy,MathAbs(b.m_y[i]));
//--- check
   if(b.m_sy>0.0)
     {
      v=1/b.m_sy;
      //--- calculation
      for(i_=0;i_<=b.m_n-1;i_++)
         b.m_y[i_]=v*b.m_y[i_];
     }
  }
//+------------------------------------------------------------------+
//| Extracts X/Y/W arrays from rational interpolant                  |
//| INPUT PARAMETERS:                                                |
//|     B   -   barycentric interpolant                              |
//| OUTPUT PARAMETERS:                                               |
//|     N   -   nodes count, N>0                                     |
//|     X   -   interpolation nodes, array[0..N-1]                   |
//|     F   -   function values, array[0..N-1]                       |
//|     W   -   barycentric weights, array[0..N-1]                   |
//+------------------------------------------------------------------+
static void CRatInt::BarycentricUnpack(CBarycentricInterpolant &b,int &n,
                                       double &x[],double &y[],double &w[])
  {
//--- create variables
   double v=0;
   int    i_=0;
//--- initialization
   n=b.m_n;
//--- allocation
   ArrayResizeAL(x,n);
   ArrayResizeAL(y,n);
   ArrayResizeAL(w,n);
//--- initialization
   v=b.m_sy;
//--- copy
   for(i_=0;i_<=n-1;i_++)
      x[i_]=b.m_x[i_];
   for(i_=0;i_<=n-1;i_++)
      y[i_]=v*b.m_y[i_];
   for(i_=0;i_<=n-1;i_++)
      w[i_]=b.m_w[i_];
  }
//+------------------------------------------------------------------+
//| Rational interpolant from X/Y/W arrays                           |
//| F(t)=SUM(i=0,n-1,w[i]*f[i]/(t-x[i])) / SUM(i=0,n-1,w[i]/(t-x[i]))|
//| INPUT PARAMETERS:                                                |
//|     X   -   interpolation nodes, array[0..N-1]                   |
//|     F   -   function values, array[0..N-1]                       |
//|     W   -   barycentric weights, array[0..N-1]                   |
//|     N   -   nodes count, N>0                                     |
//| OUTPUT PARAMETERS:                                               |
//|     B   -   barycentric interpolant built from (X, Y, W)         |
//+------------------------------------------------------------------+
static void CRatInt::BarycentricBuildXYW(double &x[],double &y[],double &w[],
                                         const int n,CBarycentricInterpolant &b)
  {
//--- create a variable
   int i_=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- fill X/Y/W
   ArrayResizeAL(b.m_x,n);
   ArrayResizeAL(b.m_y,n);
   ArrayResizeAL(b.m_w,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      b.m_x[i_]=x[i_];
   for(i_=0;i_<=n-1;i_++)
      b.m_y[i_]=y[i_];
   for(i_=0;i_<=n-1;i_++)
      b.m_w[i_]=w[i_];
   b.m_n=n;
//--- Normalize
   BarycentricNormalize(b);
  }
//+------------------------------------------------------------------+
//| Rational interpolant without poles                               |
//| The subroutine constructs the rational interpolating function    |
//| without real poles (see 'Barycentric rational interpolation with |
//| no poles and high rates of approximation', Michael S. Floater.   |
//| and Kai Hormann, for more information on this subject).          |
//| Input parameters:                                                |
//|     X   -   interpolation nodes, array[0..N-1].                  |
//|     Y   -   function values, array[0..N-1].                      |
//|     N   -   number of nodes, N>0.                                |
//|     D   -   order of the interpolation scheme, 0 <= D <= N-1.    |
//|             D<0 will cause an error.                             |
//|             D>=N it will be replaced with D=N-1.                 |
//|             if you don't know what D to choose, use small value  |
//|             about 3-5.                                           |
//| Output parameters:                                               |
//|     B   -   barycentric interpolant.                             |
//| Note:                                                            |
//|     this algorithm always succeeds and calculates the weights    |
//|     with close to machine precision.                             |
//+------------------------------------------------------------------+
static void CRatInt::BarycentricBuildFloaterHormann(double &x[],double &y[],
                                                    const int n,int d,
                                                    CBarycentricInterpolant &b)
  {
//--- create variables
   double s0=0;
   double s=0;
   double v=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    i_=0;
//--- create arrays
   int    perm[];
   double wtemp[];
   double sortrbuf[];
   double sortrbuf2[];
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(d>=0,__FUNCTION__+": incorrect D!"))
      return;
//--- Prepare
   if(d>n-1)
      d=n-1;
   b.m_n=n;
//--- special case: N=1
   if(n==1)
     {
      //--- allocation
      ArrayResizeAL(b.m_x,n);
      ArrayResizeAL(b.m_y,n);
      ArrayResizeAL(b.m_w,n);
      //--- change values
      b.m_x[0]=x[0];
      b.m_y[0]=y[0];
      b.m_w[0]=1;
      //--- function call
      BarycentricNormalize(b);
      //--- exit the function
      return;
     }
//--- Fill X/Y
   ArrayResizeAL(b.m_x,n);
   ArrayResizeAL(b.m_y,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      b.m_x[i_]=x[i_];
   for(i_=0;i_<=n-1;i_++)
      b.m_y[i_]=y[i_];
//--- function call
   CTSort::TagSortFastR(b.m_x,b.m_y,sortrbuf,sortrbuf2,n);
//--- Calculate Wk
   ArrayResizeAL(b.m_w,n);
   s0=1;
   for(k=1;k<=d;k++)
      s0=-s0;
//--- calculation
   for(k=0;k<=n-1;k++)
     {
      //--- Wk
      s=0;
      for(i=(int)(MathMax(k-d,0));i<=MathMin(k,n-1-d);i++)
        {
         v=1;
         for(j=i;j<=i+d;j++)
           {
            //--- check
            if(j!=k)
               v=v/MathAbs(b.m_x[k]-b.m_x[j]);
           }
         s=s+v;
        }
      b.m_w[k]=s0*s;
      //--- Next S0
      s0=-s0;
     }
//--- Normalize
   BarycentricNormalize(b);
  }
//+------------------------------------------------------------------+
//| Copying of the barycentric interpolant (for internal use only)   |
//| INPUT PARAMETERS:                                                |
//|     B   -   barycentric interpolant                              |
//| OUTPUT PARAMETERS:                                               |
//|     B2  -   copy(B1)                                             |
//+------------------------------------------------------------------+
static void CRatInt::BarycentricCopy(CBarycentricInterpolant &b,
                                     CBarycentricInterpolant &b2)
  {
//--- create a variable
   int i_=0;
//--- change values
   b2.m_n=b.m_n;
   b2.m_sy=b.m_sy;
//--- allocation
   ArrayResizeAL(b2.m_x,b2.m_n);
   ArrayResizeAL(b2.m_y,b2.m_n);
   ArrayResizeAL(b2.m_w,b2.m_n);
//--- copy
   for(i_=0;i_<=b2.m_n-1;i_++)
      b2.m_x[i_]=b.m_x[i_];
   for(i_=0;i_<=b2.m_n-1;i_++)
      b2.m_y[i_]=b.m_y[i_];
   for(i_=0;i_<=b2.m_n-1;i_++)
      b2.m_w[i_]=b.m_w[i_];
  }
//+------------------------------------------------------------------+
//| Normalization of barycentric interpolant:                        |
//| * B.N, B.X, B.Y and B.W are initialized                          |
//| * B.SY is NOT initialized                                        |
//| * Y[] is normalized, scaling coefficient is stored in B.SY       |
//| * W[] is normalized, no scaling coefficient is stored            |
//| * X[] is sorted                                                  |
//| Internal subroutine.                                             |
//+------------------------------------------------------------------+
static void CRatInt::BarycentricNormalize(CBarycentricInterpolant &b)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    j2=0;
   double v=0;
   int    i_=0;
//--- create arrays
   int p1[];
   int p2[];
//--- Normalize task: |Y|<=1,|W|<=1,sort X[]
   b.m_sy=0;
   for(i=0;i<=b.m_n-1;i++)
      b.m_sy=MathMax(b.m_sy,MathAbs(b.m_y[i]));
//--- check
   if(b.m_sy>0.0 && MathAbs(b.m_sy-1)>10*CMath::m_machineepsilon)
     {
      v=1/b.m_sy;
      for(i_=0;i_<=b.m_n-1;i_++)
         b.m_y[i_]=v*b.m_y[i_];
     }
//--- change value
   v=0;
   for(i=0;i<=b.m_n-1;i++)
      v=MathMax(v,MathAbs(b.m_w[i]));
//--- check
   if(v>0.0 && MathAbs(v-1)>10*CMath::m_machineepsilon)
     {
      v=1/v;
      for(i_=0;i_<=b.m_n-1;i_++)
         b.m_w[i_]=v*b.m_w[i_];
     }
   for(i=0;i<=b.m_n-2;i++)
     {
      //--- check
      if(b.m_x[i+1]<b.m_x[i])
        {
         //--- function call
         CTSort::TagSort(b.m_x,b.m_n,p1,p2);
         //--- calculation
         for(j=0;j<=b.m_n-1;j++)
           {
            j2=p2[j];
            v=b.m_y[j];
            b.m_y[j]=b.m_y[j2];
            b.m_y[j2]=v;
            v=b.m_w[j];
            b.m_w[j]=b.m_w[j2];
            b.m_w[j2]=v;
           }
         break;
        }
     }
  }
//+------------------------------------------------------------------+
//| Polynomial interpolant                                           |
//+------------------------------------------------------------------+
class CPolInt
  {
public:
   //--- constructor, destructor
                     CPolInt(void);
                    ~CPolInt(void);
   //--- methods
   static void       PolynomialBar2Cheb(CBarycentricInterpolant &p,const double a,const double b,double &t[]);
   static void       PolynomialCheb2Bar(double &t[],const int n,const double a,const double b,CBarycentricInterpolant &p);
   static void       PolynomialBar2Pow(CBarycentricInterpolant &p,const double c,const double s,double &a[]);
   static void       PolynomialPow2Bar(double &a[],const int n,const double c,const double s,CBarycentricInterpolant &p);
   static void       PolynomialBuild(double &cx[],double &cy[],const int n,CBarycentricInterpolant &p);
   static void       PolynomialBuildEqDist(const double a,const double b,double &y[],const int n,CBarycentricInterpolant &p);
   static void       PolynomialBuildCheb1(const double a,const double b,double &y[],const int n,CBarycentricInterpolant &p);
   static void       PolynomialBuildCheb2(const double a,const double b,double &y[],const int n,CBarycentricInterpolant &p);
   static double     PolynomialCalcEqDist(const double a,const double b,double &f[],const int n,const double t);
   static double     PolynomialCalcCheb1(const double a,const double b,double &f[],const int n,double t);
   static double     PolynomialCalcCheb2(const double a,const double b,double &f[],const int n,double t);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CPolInt::CPolInt(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CPolInt::~CPolInt(void)
  {

  }
//+------------------------------------------------------------------+
//| Conversion from barycentric representation to Chebyshev basis.   |
//| This function has O(N^2) complexity.                             |
//| INPUT PARAMETERS:                                                |
//|     P   -   polynomial in barycentric form                       |
//|     A,B -   base interval for Chebyshev polynomials (see below)  |
//|             A<>B                                                 |
//| OUTPUT PARAMETERS                                                |
//|     T   -   coefficients of Chebyshev representation;            |
//|             P(x) = sum { T[i]*Ti(2*(x-A)/(B-A)-1), i=0..N-1 },   |
//|             where Ti - I-th Chebyshev polynomial.                |
//| NOTES:                                                           |
//|     barycentric interpolant passed as P may be either polynomial |
//|     obtained from polynomial interpolation/ fitting or rational  |
//|     function which is NOT polynomial. We can't distinguish       |
//|     between these two cases, and this algorithm just tries to    |
//|     work assuming that P IS a polynomial. If not, algorithm will |
//|     return results, but they won't have any meaning.             |
//+------------------------------------------------------------------+
static void CPolInt::PolynomialBar2Cheb(CBarycentricInterpolant &p,
                                        const double a,const double b,
                                        double &t[])
  {
//--- create variables
   int    i=0;
   int    k=0;
   double v=0;
   int    i_=0;
//--- create arrays
   double vp[];
   double vx[];
   double tk[];
   double tk1[];
//--- check
   if(!CAp::Assert(CMath::IsFinite(a),__FUNCTION__+": A is not finite!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(b),__FUNCTION__+": B is not finite!"))
      return;
//--- check
   if(!CAp::Assert(a!=b,__FUNCTION__+": A=B!"))
      return;
//--- check
   if(!CAp::Assert(p.m_n>0,__FUNCTION__+": P is not correctly initialized barycentric interpolant!"))
      return;
//--- Calculate function values on a Chebyshev grid
   ArrayResizeAL(vp,p.m_n);
   ArrayResizeAL(vx,p.m_n);
   for(i=0;i<=p.m_n-1;i++)
     {
      vx[i]=MathCos(M_PI*(i+0.5)/p.m_n);
      vp[i]=CRatInt::BarycentricCalc(p,0.5*(vx[i]+1)*(b-a)+a);
     }
//--- T[0]
   ArrayResizeAL(t,p.m_n);
   v=0;
   for(i=0;i<=p.m_n-1;i++)
      v=v+vp[i];
   t[0]=v/p.m_n;
//--- other T's.
//--- NOTES:
//--- 1. TK stores T{k} on VX,TK1 stores T{k-1} on VX
//--- 2. we can do same calculations with fast DCT,but it
//---    * adds dependencies
//---    * still leaves us with O(N^2) algorithm because
//---      preparation of function values is O(N^2) process
   if(p.m_n>1)
     {
      //--- allocation
      ArrayResizeAL(tk,p.m_n);
      ArrayResizeAL(tk1,p.m_n);
      for(i=0;i<=p.m_n-1;i++)
        {
         tk[i]=vx[i];
         tk1[i]=1;
        }
      //--- calculation
      for(k=1;k<=p.m_n-1;k++)
        {
         //--- calculate discrete product of function vector and TK
         v=0.0;
         for(i_=0;i_<=p.m_n-1;i_++)
            v+=tk[i_]*vp[i_];
         t[k]=v/(0.5*p.m_n);
         //--- Update TK and TK1
         for(i=0;i<=p.m_n-1;i++)
           {
            v=2*vx[i]*tk[i]-tk1[i];
            tk1[i]=tk[i];
            tk[i]=v;
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Conversion from Chebyshev basis to barycentric representation.   |
//| This function has O(N^2) complexity.                             |
//| INPUT PARAMETERS:                                                |
//|     T   -   coefficients of Chebyshev representation;            |
//|             P(x) = sum { T[i]*Ti(2*(x-A)/(B-A)-1), i=0..N },     |
//|             where Ti - I-th Chebyshev polynomial.                |
//|     N   -   number of coefficients:                              |
//|             * if given, only leading N elements of T are used    |
//|             * if not given, automatically determined from size   |
//|               of T                                               |
//|     A,B -   base interval for Chebyshev polynomials (see above)  |
//|             A<B                                                  |
//| OUTPUT PARAMETERS                                                |
//|     P   -   polynomial in barycentric form                       |
//+------------------------------------------------------------------+
static void CPolInt::PolynomialCheb2Bar(double &t[],const int n,const double a,
                                        const double b,CBarycentricInterpolant &p)
  {
//--- create variables
   int    i=0;
   int    k=0;
   double tk=0;
   double tk1=0;
   double vx=0;
   double vy=0;
   double v=0;
//--- create array
   double y[];
//--- check
   if(!CAp::Assert(CMath::IsFinite(a),__FUNCTION__+": A is not finite!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(b),__FUNCTION__+": B is not finite!"))
      return;
//--- check
   if(!CAp::Assert(a!=b,__FUNCTION__+": A=B!"))
      return;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(t)>=n,__FUNCTION__+": Length(T)<N"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(t,n),__FUNCTION__+": T[] contains INF or NAN"))
      return;
//--- Calculate function values on a Chebyshev grid spanning [-1,+1]
   ArrayResizeAL(y,n);
   for(i=0;i<=n-1;i++)
     {
      //--- Calculate value on a grid spanning [-1,+1]
      vx=MathCos(M_PI*(i+0.5)/n);
      vy=t[0];
      tk1=1;
      tk=vx;
      //--- change values
      for(k=1;k<=n-1;k++)
        {
         vy=vy+t[k]*tk;
         v=2*vx*tk-tk1;
         tk1=tk;
         tk=v;
        }
      y[i]=vy;
     }
//--- Build barycentric interpolant,map grid from [-1,+1] to [A,B]
   PolynomialBuildCheb1(a,b,y,n,p);
  }
//+------------------------------------------------------------------+
//| Conversion from barycentric representation to power basis.       |
//| This function has O(N^2) complexity.                             |
//| INPUT PARAMETERS:                                                |
//|     P   -   polynomial in barycentric form                       |
//|     C   -   offset (see below); 0.0 is used as default value.    |
//|     S   -   scale (see below); 1.0 is used as default value.     |
//|             S<>0.                                                |
//| OUTPUT PARAMETERS                                                |
//|     A   -   coefficients,                                        |
//|             P(x) = sum { A[i]*((X-C)/S)^i, i=0..N-1 }            |
//|     N   -   number of coefficients (polynomial degree plus 1)    |
//| NOTES:                                                           |
//| 1.  this function accepts offset and scale, which can be set to  |
//|     improve numerical properties of polynomial. For example, if  |
//|     P was obtained as result of interpolation on [-1,+1], you can|
//|     set C=0 and S=1 and represent P as sum of 1, x, x^2, x^3 and |
//|     so on. In most cases you it is exactly what you need.        |
//|     However, if your interpolation model was built on [999,1001],|
//|     you will see significant growth of numerical errors when     |
//|     using {1, x, x^2, x^3} as basis. Representing P as sum of 1, |
//|     (x-1000), (x-1000)^2, (x-1000)^3 will be better option. Such |
//|     representation can be  obtained  by  using 1000.0 as offset  |
//|     C and 1.0 as scale S.                                        |
//| 2.  power basis is ill-conditioned and tricks described above    |
//|     can't solve this problem completely. This function  will     |
//|     return coefficients in any case, but for N>8 they will become|
//|     unreliable. However, N's less than 5 are pretty safe.        |
//| 3.  barycentric interpolant passed as P may be either polynomial |
//|     obtained from polynomial interpolation/ fitting or rational  |
//|     function which is NOT polynomial. We can't distinguish       |
//|     between these two cases, and this algorithm just tries to    |
//|     work assuming that P IS a polynomial. If not, algorithm will |
//|     return results, but they won't have any meaning.             |
//+------------------------------------------------------------------+
static void CPolInt::PolynomialBar2Pow(CBarycentricInterpolant &p,
                                       const double c,const double s,
                                       double &a[])
  {
//--- create variables
   int    i=0;
   int    k=0;
   double e=0;
   double d=0;
   double v=0;
   int    i_=0;
//--- create arrays
   double vp[];
   double vx[];
   double tk[];
   double tk1[];
   double t[];
//--- check
   if(!CAp::Assert(CMath::IsFinite(c),__FUNCTION__+": C is not finite!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(s),__FUNCTION__+": S is not finite!"))
      return;
//--- check
   if(!CAp::Assert(s!=0.0,__FUNCTION__+": S=0!"))
      return;
//--- check
   if(!CAp::Assert(p.m_n>0,__FUNCTION__+": P is not correctly initialized barycentric interpolant!"))
      return;
//--- Calculate function values on a Chebyshev grid
   ArrayResizeAL(vp,p.m_n);
   ArrayResizeAL(vx,p.m_n);
   for(i=0;i<=p.m_n-1;i++)
     {
      vx[i]=MathCos(M_PI*(i+0.5)/p.m_n);
      vp[i]=CRatInt::BarycentricCalc(p,s*vx[i]+c);
     }
//--- T[0]
   ArrayResizeAL(t,p.m_n);
   v=0;
   for(i=0;i<=p.m_n-1;i++)
      v=v+vp[i];
   t[0]=v/p.m_n;
//--- other T's.
//--- NOTES:
//--- 1. TK stores T{k} on VX,TK1 stores T{k-1} on VX
//--- 2. we can do same calculations with fast DCT,but it
//---    * adds dependencies
//---    * still leaves us with O(N^2) algorithm because
//---      preparation of function values is O(N^2) process
   if(p.m_n>1)
     {
      //--- allocation
      ArrayResizeAL(tk,p.m_n);
      ArrayResizeAL(tk1,p.m_n);
      for(i=0;i<=p.m_n-1;i++)
        {
         tk[i]=vx[i];
         tk1[i]=1;
        }
      //--- calculation
      for(k=1;k<=p.m_n-1;k++)
        {
         //--- calculate discrete product of function vector and TK
         v=0.0;
         for(i_=0;i_<=p.m_n-1;i_++)
            v+=tk[i_]*vp[i_];
         t[k]=v/(0.5*p.m_n);
         //--- Update TK and TK1
         for(i=0;i<=p.m_n-1;i++)
           {
            v=2*vx[i]*tk[i]-tk1[i];
            tk1[i]=tk[i];
            tk[i]=v;
           }
        }
     }
//--- Convert from Chebyshev basis to power basis
   ArrayResizeAL(a,p.m_n);
   for(i=0;i<=p.m_n-1;i++)
      a[i]=0;
   d=0;
//--- calculation
   for(i=0;i<=p.m_n-1;i++)
     {
      for(k=i;k<=p.m_n-1;k++)
        {
         e=a[k];
         a[k]=0;
         //--- check
         if(i<=1 && k==i)
            a[k]=1;
         else
           {
            //--- check
            if(i!=0)
               a[k]=2*d;
            //--- check
            if(k>i+1)
               a[k]=a[k]-a[k-2];
           }
         d=e;
        }
      //--- change values
      d=a[i];
      e=0;
      k=i;
      //--- cycle
      while(k<=p.m_n-1)
        {
         e=e+a[k]*t[k];
         k=k+2;
        }
      a[i]=e;
     }
  }
//+------------------------------------------------------------------+
//| Conversion from power basis to barycentric representation.       |
//| This function has O(N^2) complexity.                             |
//| INPUT PARAMETERS:                                                |
//|     A   -   coefficients, P(x)=sum { A[i]*((X-C)/S)^i, i=0..N-1 }|
//|     N   -   number of coefficients (polynomial degree plus 1)    |
//|             * if given, only leading N elements of A are used    |
//|             * if not given, automatically determined from size   |
//|               of A                                               |
//|     C   -   offset (see below); 0.0 is used as default value.    |
//|     S   -   scale (see below); 1.0 is used as default value.     |
//|             S<>0.                                                |
//| OUTPUT PARAMETERS                                                |
//|     P   -   polynomial in barycentric form                       |
//| NOTES:                                                           |
//| 1.  this function accepts offset and scale, which can be set to  |
//|     improve numerical properties of polynomial. For example, if  |
//|     you interpolate on [-1,+1], you can set C=0 and S=1 and      |
//|     convert from sum of 1, x, x^2, x^3 and so on. In most cases  |
//|     you it is exactly what you need.                             |
//|     However, if your interpolation model was built on [999,1001],|
//|     you will see significant growth of numerical errors when     |
//|     using {1, x, x^2, x^3} as input basis. Converting from sum   |
//|     of 1, (x-1000), (x-1000)^2, (x-1000)^3 will be better option |
//|     (you have to specify 1000.0 as offset C and 1.0 as scale S). |
//| 2.  power basis is ill-conditioned and tricks described above    |
//|     can't solve this problem completely. This function will      |
//|     return barycentric model in any case, but for N>8 accuracy   |
//|     well degrade. However, N's less than 5 are pretty safe.      |
//+------------------------------------------------------------------+
static void CPolInt::PolynomialPow2Bar(double &a[],const int n,const double c,
                                       const double s,CBarycentricInterpolant &p)
  {
//--- create variables
   int    i=0;
   int    k=0;
   double vx=0;
   double vy=0;
   double px=0;
//--- create array
   double y[];
//--- check
   if(!CAp::Assert(CMath::IsFinite(c),__FUNCTION__+": C is not finite!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(s),__FUNCTION__+": S is not finite!"))
      return;
//--- check
   if(!CAp::Assert(s!=0.0,__FUNCTION__+": S is zero!"))
      return;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(a)>=n,__FUNCTION__+": Length(A)<N"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(a,n),__FUNCTION__+": A[] contains INF or NAN"))
      return;
//--- Calculate function values on a Chebyshev grid spanning [-1,+1]
   ArrayResizeAL(y,n);
   for(i=0;i<=n-1;i++)
     {
      //--- Calculate value on a grid spanning [-1,+1]
      vx=MathCos(M_PI*(i+0.5)/n);
      vy=a[0];
      px=vx;
      //--- calculation
      for(k=1;k<=n-1;k++)
        {
         vy=vy+px*a[k];
         px=px*vx;
        }
      y[i]=vy;
     }
//--- Build barycentric interpolant,map grid from [-1,+1] to [A,B]
   PolynomialBuildCheb1(c-s,c+s,y,n,p);
  }
//+------------------------------------------------------------------+
//| Lagrange intepolant: generation of the model on the general grid.|
//| This function has O(N^2) complexity.                             |
//| INPUT PARAMETERS:                                                |
//|     X   -   abscissas, array[0..N-1]                             |
//|     Y   -   function values, array[0..N-1]                       |
//|     N   -   number of points, N>=1                               |
//| OUTPUT PARAMETERS                                                |
//|     P   -   barycentric model which represents Lagrange          |
//|             interpolant (see ratint unit info and                |
//|             BarycentricCalc() description for more information). |
//+------------------------------------------------------------------+
static void CPolInt::PolynomialBuild(double &cx[],double &cy[],const int n,
                                     CBarycentricInterpolant &p)
  {
//--- create variables
   int    j=0;
   int    k=0;
   double b=0;
   double a=0;
   double v=0;
   double mx=0;
   int    i_=0;
//--- create arrays
   double w[];
   double sortrbuf[];
   double sortrbuf2[];
   double x[];
   double y[];
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- function call
   CTSort::TagSortFastR(x,y,sortrbuf,sortrbuf2,n);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(x,n),__FUNCTION__+": at least two consequent points are too close!"))
      return;
//--- calculate W[j]
//--- multi-pass algorithm is used to avoid overflow
   ArrayResizeAL(w,n);
   a=x[0];
   b=x[0];
   for(j=0;j<=n-1;j++)
     {
      w[j]=1;
      a=MathMin(a,x[j]);
      b=MathMax(b,x[j]);
     }
//--- calculation
   for(k=0;k<=n-1;k++)
     {
      //--- W[K] is used instead of 0.0 because
      //--- cycle on J does not touch K-th element
      //--- and we MUST get maximum from ALL elements
      mx=MathAbs(w[k]);
      for(j=0;j<=n-1;j++)
        {
         //--- check
         if(j!=k)
           {
            v=(b-a)/(x[j]-x[k]);
            w[j]=w[j]*v;
            mx=MathMax(mx,MathAbs(w[j]));
           }
        }
      //--- check
      if(k%5==0)
        {
         //--- every 5-th run we renormalize W[]
         v=1/mx;
         for(i_=0;i_<=n-1;i_++)
            w[i_]=v*w[i_];
        }
     }
//--- function call
   CRatInt::BarycentricBuildXYW(x,y,w,n,p);
  }
//+------------------------------------------------------------------+
//| Lagrange intepolant: generation of the model on equidistant grid.|
//| This function has O(N) complexity.                               |
//| INPUT PARAMETERS:                                                |
//|     A   -   left boundary of [A,B]                               |
//|     B   -   right boundary of [A,B]                              |
//|     Y   -   function values at the nodes, array[0..N-1]          |
//|     N   -   number of points, N>=1                               |
//|             for N=1 a constant model is constructed.             |
//| OUTPUT PARAMETERS                                                |
//|     P   -   barycentric model which represents Lagrange          |
//|             interpolant (see ratint unit info and                |
//|             BarycentricCalc() description for more information). |
//+------------------------------------------------------------------+
static void CPolInt::PolynomialBuildEqDist(const double a,const double b,
                                           double &y[],const int n,
                                           CBarycentricInterpolant &p)
  {
//--- create variables
   int    i=0;
   double v=0;
//--- create arrays
   double w[];
   double x[];
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(a),__FUNCTION__+": A is infinite or NaN!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(b),__FUNCTION__+": B is infinite or NaN!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(b!=a,__FUNCTION__+": B=A!"))
      return;
//--- check
   if(!CAp::Assert((double)(a+(b-a)/n)!=a,__FUNCTION__+": B is too close to A!"))
      return;
//--- Special case: N=1
   if(n==1)
     {
      //--- allocation
      ArrayResizeAL(x,1);
      ArrayResizeAL(w,1);
      x[0]=0.5*(b+a);
      w[0]=1;
      //--- function call
      CRatInt::BarycentricBuildXYW(x,y,w,1,p);
      //--- exit the function
      return;
     }
//--- general case
   ArrayResizeAL(x,n);
   ArrayResizeAL(w,n);
   v=1;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      w[i]=v;
      x[i]=a+(b-a)*i/(n-1);
      v=-(v*(n-1-i));
      v=v/(i+1);
     }
//--- function call
   CRatInt::BarycentricBuildXYW(x,y,w,n,p);
  }
//+------------------------------------------------------------------+
//| Lagrange intepolant on Chebyshev grid (first kind).              |
//| This function has O(N) complexity.                               |
//| INPUT PARAMETERS:                                                |
//|     A   -   left boundary of [A,B]                               |
//|     B   -   right boundary of [A,B]                              |
//|     Y   -   function values at the nodes, array[0..N-1],         |
//|             Y[I] = Y(0.5*(B+A) + 0.5*(B-A)*Cos(PI*(2*i+1)/(2*n)))|
//|     N   -   number of points, N>=1                               |
//|             for N=1 a constant model is constructed.             |
//| OUTPUT PARAMETERS                                                |
//|     P   -   barycentric model which represents Lagrange          |
//|             interpolant (see ratint unit info and                |
//|             BarycentricCalc() description for more information). |
//+------------------------------------------------------------------+
static void CPolInt::PolynomialBuildCheb1(const double a,const double b,
                                          double &y[],const int n,
                                          CBarycentricInterpolant &p)
  {
//--- create variables
   int    i=0;
   double v=0;
   double t=0;
//--- create arrays
   double w[];
   double x[];
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(a),__FUNCTION__+": A is infinite or NaN!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(b),__FUNCTION__+": B is infinite or NaN!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(b!=a,__FUNCTION__+": B=A!"))
      return;
//--- Special case: N=1
   if(n==1)
     {
      //--- allocation
      ArrayResizeAL(x,1);
      ArrayResizeAL(w,1);
      x[0]=0.5*(b+a);
      w[0]=1;
      //--- function call
      CRatInt::BarycentricBuildXYW(x,y,w,1,p);
      //--- exit the function
      return;
     }
//--- general case
   ArrayResizeAL(x,n);
   ArrayResizeAL(w,n);
   v=1;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      t=MathTan(0.5*M_PI*(2*i+1)/(2*n));
      w[i]=2*v*t/(1+CMath::Sqr(t));
      x[i]=0.5*(b+a)+0.5*(b-a)*(1-CMath::Sqr(t))/(1+CMath::Sqr(t));
      v=-v;
     }
//--- function call
   CRatInt::BarycentricBuildXYW(x,y,w,n,p);
  }
//+------------------------------------------------------------------+
//| Lagrange intepolant on Chebyshev grid (second kind).             |
//| This function has O(N) complexity.                               |
//| INPUT PARAMETERS:                                                |
//|     A   -   left boundary of [A,B]                               |
//|     B   -   right boundary of [A,B]                              |
//|     Y   -   function values at the nodes, array[0..N-1],         |
//|             Y[I] = Y(0.5*(B+A) + 0.5*(B-A)*Cos(PI*i/(n-1)))      |
//|     N   -   number of points, N>=1                               |
//|             for N=1 a constant model is constructed.             |
//| OUTPUT PARAMETERS                                                |
//|     P   -   barycentric model which represents Lagrange          |
//|             interpolant (see ratint unit info and                |
//|             BarycentricCalc() description for more information). |
//+------------------------------------------------------------------+
static void CPolInt::PolynomialBuildCheb2(const double a,const double b,
                                          double &y[],const int n,
                                          CBarycentricInterpolant &p)
  {
//--- create variables
   int    i=0;
   double v=0;
//--- create arrays
   double w[];
   double x[];
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(a),__FUNCTION__+": A is infinite or NaN!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(b),__FUNCTION__+": B is infinite or NaN!"))
      return;
//--- check
   if(!CAp::Assert(b!=a,__FUNCTION__+": B=A!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- Special case: N=1
   if(n==1)
     {
      //--- allocation
      ArrayResizeAL(x,1);
      ArrayResizeAL(w,1);
      x[0]=0.5*(b+a);
      w[0]=1;
      //--- function call
      CRatInt::BarycentricBuildXYW(x,y,w,1,p);
      //--- exit the function
      return;
     }
//--- general case
   ArrayResizeAL(x,n);
   ArrayResizeAL(w,n);
   v=1;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(i==0 || i==n-1)
         w[i]=v*0.5;
      else
         w[i]=v;
      x[i]=0.5*(b+a)+0.5*(b-a)*MathCos(M_PI*i/(n-1));
      v=-v;
     }
//--- function call
   CRatInt::BarycentricBuildXYW(x,y,w,n,p);
  }
//+------------------------------------------------------------------+
//| Fast equidistant polynomial interpolation function with O(N)     |
//| complexity                                                       |
//| INPUT PARAMETERS:                                                |
//|     A   -   left boundary of [A,B]                               |
//|     B   -   right boundary of [A,B]                              |
//|     F   -   function values, array[0..N-1]                       |
//|     N   -   number of points on equidistant grid, N>=1           |
//|             for N=1 a constant model is constructed.             |
//|     T   -   position where P(x) is calculated                    |
//| RESULT                                                           |
//|     value of the Lagrange interpolant at T                       |
//| IMPORTANT                                                        |
//|     this function provides fast interface which is not           |
//|     overflow-safe nor it is very precise.                        |
//|     the best option is to use PolynomialBuildEqDist() or         |
//|     BarycentricCalc() subroutines unless you are pretty sure that|
//|     your data will not result in overflow.                       |
//+------------------------------------------------------------------+
static double CPolInt::PolynomialCalcEqDist(const double a,const double b,
                                            double &f[],const int n,const double t)
  {
//--- create variables
   double s1=0;
   double s2=0;
   double v=0;
   double threshold=0;
   double s=0;
   double h=0;
   int    i=0;
   int    j=0;
   double w=0;
   double x=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Len(f)>=n,__FUNCTION__+": Length(F)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CMath::IsFinite(a),__FUNCTION__+": A is infinite or NaN!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CMath::IsFinite(b),__FUNCTION__+": B is infinite or NaN!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(f,n),__FUNCTION__+": F contains infinite or NaN values!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(b!=a,__FUNCTION__+": B=A!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(!CInfOrNaN::IsInfinity(t),__FUNCTION__+": T is infinite!"))
      return(EMPTY_VALUE);
//--- Special case: T is NAN
   if(CInfOrNaN::IsNaN(t))
      return(CInfOrNaN::NaN());
//--- Special case: N=1
   if(n==1)
      return(f[0]);
//--- First,decide: should we use "safe" formula (guarded
//--- against overflow) or fast one?
   threshold=MathSqrt(CMath::m_minrealnumber);
   j=0;
   s=t-a;
//--- calculation
   for(i=1;i<=n-1;i++)
     {
      x=a+(double)i/(double)(n-1)*(b-a);
      //--- check
      if(MathAbs(t-x)<MathAbs(s))
        {
         s=t-x;
         j=i;
        }
     }
//--- check
   if(s==0.0)
      return(f[j]);
//--- check
   if(MathAbs(s)>threshold)
     {
      //--- use fast formula
      j=-1;
      s=1.0;
     }
//--- Calculate using safe or fast barycentric formula
   s1=0;
   s2=0;
   w=1.0;
   h=(b-a)/(n-1);
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(i!=j)
        {
         v=s*w/(t-(a+i*h));
         s1=s1+v*f[i];
         s2=s2+v;
        }
      else
        {
         v=w;
         s1=s1+v*f[i];
         s2=s2+v;
        }
      //--- change values
      w=-(w*(n-1-i));
      w=w/(i+1);
     }
//--- return result
   return(s1/s2);
  }
//+------------------------------------------------------------------+
//| Fast polynomial interpolation function on Chebyshev points (first|
//| kind) with O(N) complexity.                                      |
//| INPUT PARAMETERS:                                                |
//|     A   -   left boundary of [A,B]                               |
//|     B   -   right boundary of [A,B]                              |
//|     F   -   function values, array[0..N-1]                       |
//|     N   -   number of points on Chebyshev grid (first kind),     |
//|             X[i] = 0.5*(B+A) + 0.5*(B-A)*Cos(PI*(2*i+1)/(2*n))   |
//|             for N=1 a constant model is constructed.             |
//|     T   -   position where P(x) is calculated                    |
//| RESULT                                                           |
//|     value of the Lagrange interpolant at T                       |
//| IMPORTANT                                                        |
//|     this function provides fast interface which is not           |
//|     overflow-safe nor it is very precise                         |
//|     the best option is to use  PolIntBuildCheb1() or             |
//|     BarycentricCalc() subroutines unless you are pretty sure that|
//|     your data will not result in overflow.                       |
//+------------------------------------------------------------------+
static double CPolInt::PolynomialCalcCheb1(const double a,const double b,
                                           double &f[],const int n,double t)
  {
//--- create variables
   double s1=0;
   double s2=0;
   double v=0;
   double threshold=0;
   double s=0;
   int    i=0;
   int    j=0;
   double a0=0;
   double delta=0;
   double alpha=0;
   double beta=0;
   double ca=0;
   double sa=0;
   double tempc=0;
   double temps=0;
   double x=0;
   double w=0;
   double p1=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Len(f)>=n,__FUNCTION__+": Length(F)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CMath::IsFinite(a),__FUNCTION__+": A is infinite or NaN!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CMath::IsFinite(b),__FUNCTION__+": B is infinite or NaN!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(f,n),__FUNCTION__+": F contains infinite or NaN values!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(b!=a,__FUNCTION__+": B=A!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(!CInfOrNaN::IsInfinity(t),__FUNCTION__+": T is infinite!"))
      return(EMPTY_VALUE);
//--- Special case: T is NAN
   if(CInfOrNaN::IsNaN(t))
      return(CInfOrNaN::NaN());
//--- Special case: N=1
   if(n==1)
      return(f[0]);
//--- Prepare information for the recurrence formula
//--- used to calculate sin(pi*(2j+1)/(2n+2)) and
//--- cos(pi*(2j+1)/(2n+2)):
//--- A0=pi/(2n+2)
//--- Delta=pi/(n+1)
//--- Alpha=2 sin^2 (Delta/2)
//--- Beta=sin(Delta)
//--- so that sin(..)=sin(A0+j*delta) and cos(..)=cos(A0+j*delta).
//--- Then we use
//--- sin(x+delta)=sin(x) - (alpha*sin(x) - beta*cos(x))
//--- cos(x+delta)=cos(x) - (alpha*cos(x) - beta*sin(x))
//--- to repeatedly calculate sin(..) and cos(..).
   threshold=MathSqrt(CMath::m_minrealnumber);
   t=(t-0.5*(a+b))/(0.5*(b-a));
   a0=M_PI/(2*(n-1)+2);
   delta=2*M_PI/(2*(n-1)+2);
   alpha=2*CMath::Sqr(MathSin(delta/2));
   beta=MathSin(delta);
//--- First, decide: should we use "safe" formula (guarded
//--- against overflow) or fast one?
   ca=MathCos(a0);
   sa=MathSin(a0);
   j=0;
   x=ca;
   s=t-x;
//--- calculation
   for(i=1;i<=n-1;i++)
     {
      //--- Next X[i]
      temps=sa-(alpha*sa-beta*ca);
      tempc=ca-(alpha*ca+beta*sa);
      sa=temps;
      ca=tempc;
      x=ca;
      //--- Use X[i]
      if(MathAbs(t-x)<MathAbs(s))
        {
         s=t-x;
         j=i;
        }
     }
//--- check
   if(s==0.0)
      return(f[j]);
//--- check
   if(MathAbs(s)>threshold)
     {
      //--- use fast formula
      j=-1;
      s=1.0;
     }
//--- Calculate using safe or fast barycentric formula
   s1=0;
   s2=0;
   ca=MathCos(a0);
   sa=MathSin(a0);
   p1=1.0;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      //--- Calculate X[i],W[i]
      x=ca;
      w=p1*sa;
      //--- Proceed
      if(i!=j)
        {
         v=s*w/(t-x);
         s1=s1+v*f[i];
         s2=s2+v;
        }
      else
        {
         v=w;
         s1=s1+v*f[i];
         s2=s2+v;
        }
      //--- Next CA,SA,P1
      temps=sa-(alpha*sa-beta*ca);
      tempc=ca-(alpha*ca+beta*sa);
      sa=temps;
      ca=tempc;
      p1=-p1;
     }
//--- return result
   return(s1/s2);
  }
//+------------------------------------------------------------------+
//| Fast polynomial interpolation function on Chebyshev points       |
//| (second kind) with O(N) complexity.                              |
//| INPUT PARAMETERS:                                                |
//|     A   -   left boundary of [A,B]                               |
//|     B   -   right boundary of [A,B]                              |
//|     F   -   function values, array[0..N-1]                       |
//|     N   -   number of points on Chebyshev grid (second kind),    |
//|             X[i] = 0.5*(B+A) + 0.5*(B-A)*Cos(PI*i/(n-1))         |
//|             for N=1 a constant model is constructed.             |
//|     T   -   position where P(x) is calculated                    |
//| RESULT                                                           |
//|     value of the Lagrange interpolant at T                       |
//| IMPORTANT                                                        |
//|     this function provides fast interface which is not           |
//|     overflow-safe nor it is very precise.                        |
//|     the best option is to use PolIntBuildCheb2() or              |
//|     BarycentricCalc() subroutines unless you are pretty sure that|
//|     your data will not result in overflow.                       |
//+------------------------------------------------------------------+
static double CPolInt::PolynomialCalcCheb2(const double a,const double b,
                                           double &f[],const int n,double t)
  {
//--- create variables
   double s1=0;
   double s2=0;
   double v=0;
   double threshold=0;
   double s=0;
   int    i=0;
   int    j=0;
   double a0=0;
   double delta=0;
   double alpha=0;
   double beta=0;
   double ca=0;
   double sa=0;
   double tempc=0;
   double temps=0;
   double x=0;
   double w=0;
   double p1=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Len(f)>=n,__FUNCTION__+": Length(F)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CMath::IsFinite(a),__FUNCTION__+": A is infinite or NaN!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CMath::IsFinite(b),__FUNCTION__+": B is infinite or NaN!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(b!=a,__FUNCTION__+": B=A!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(f,n),__FUNCTION__+": F contains infinite or NaN values!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(!CInfOrNaN::IsInfinity(t),__FUNCTION__+": T is infinite!"))
      return(EMPTY_VALUE);
//--- Special case: T is NAN
   if(CInfOrNaN::IsNaN(t))
      return(CInfOrNaN::NaN());
//--- Special case: N=1
   if(n==1)
      return(f[0]);
//--- Prepare information for the recurrence formula
//--- used to calculate sin(pi*i/n) and
//--- cos(pi*i/n):
//--- A0=0
//--- Delta=pi/n
//--- Alpha=2 sin^2 (Delta/2)
//--- Beta=sin(Delta)
//--- so that sin(..)=sin(A0+j*delta) and cos(..)=cos(A0+j*delta).
//--- Then we use
//--- sin(x+delta)=sin(x) - (alpha*sin(x) - beta*cos(x))
//--- cos(x+delta)=cos(x) - (alpha*cos(x) - beta*sin(x))
//--- to repeatedly calculate sin(..) and cos(..).
   threshold=MathSqrt(CMath::m_minrealnumber);
   t=(t-0.5*(a+b))/(0.5*(b-a));
   a0=0.0;
   delta=M_PI/(n-1);
   alpha=2*CMath::Sqr(MathSin(delta/2));
   beta=MathSin(delta);
//--- First,decide: should we use "safe" formula (guarded
//--- against overflow) or fast one?
   ca=MathCos(a0);
   sa=MathSin(a0);
   j=0;
   x=ca;
   s=t-x;
//--- calculation
   for(i=1;i<=n-1;i++)
     {
      //--- Next X[i]
      temps=sa-(alpha*sa-beta*ca);
      tempc=ca-(alpha*ca+beta*sa);
      sa=temps;
      ca=tempc;
      x=ca;
      //--- Use X[i]
      if(MathAbs(t-x)<MathAbs(s))
        {
         s=t-x;
         j=i;
        }
     }
//--- check
   if(s==0.0)
      return(f[j]);
//--- check
   if(MathAbs(s)>threshold)
     {
      //--- use fast formula
      j=-1;
      s=1.0;
     }
//--- Calculate using safe or fast barycentric formula
   s1=0;
   s2=0;
   ca=MathCos(a0);
   sa=MathSin(a0);
   p1=1.0;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      //--- Calculate X[i],W[i]
      x=ca;
      //--- check
      if(i==0 || i==n-1)
         w=0.5*p1;
      else
         w=1.0*p1;
      //--- Proceed
      if(i!=j)
        {
         v=s*w/(t-x);
         s1=s1+v*f[i];
         s2=s2+v;
        }
      else
        {
         v=w;
         s1=s1+v*f[i];
         s2=s2+v;
        }
      //--- Next CA,SA,P1
      temps=sa-(alpha*sa-beta*ca);
      tempc=ca-(alpha*ca+beta*sa);
      sa=temps;
      ca=tempc;
      p1=-p1;
     }
//--- return result
   return(s1/s2);
  }
//+------------------------------------------------------------------+
//| 1-dimensional spline inteprolant                                 |
//+------------------------------------------------------------------+
class CSpline1DInterpolant
  {
public:
   //--- variables
   bool              m_periodic;
   int               m_n;
   int               m_k;
   //--- arrays
   double            m_x[];
   double            m_c[];
   //--- constructor, destructor
                     CSpline1DInterpolant(void);
                    ~CSpline1DInterpolant(void);
   //--- copy
   void              Copy(CSpline1DInterpolant &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSpline1DInterpolant::CSpline1DInterpolant(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSpline1DInterpolant::~CSpline1DInterpolant(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CSpline1DInterpolant::Copy(CSpline1DInterpolant &obj)
  {
//--- copy variables
   m_periodic=obj.m_periodic;
   m_n=obj.m_n;
   m_k=obj.m_k;
//--- copy arrays
   ArrayCopy(m_x,obj.m_x);
   ArrayCopy(m_c,obj.m_c);
  }
//+------------------------------------------------------------------+
//| 1-dimensional spline inteprolant                                 |
//+------------------------------------------------------------------+
class CSpline1DInterpolantShell
  {
private:
   CSpline1DInterpolant m_innerobj;
public:
   //--- constructors, destructor
                     CSpline1DInterpolantShell(void);
                     CSpline1DInterpolantShell(CSpline1DInterpolant &obj);
                    ~CSpline1DInterpolantShell(void);
   //--- method
   CSpline1DInterpolant *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSpline1DInterpolantShell::CSpline1DInterpolantShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CSpline1DInterpolantShell::CSpline1DInterpolantShell(CSpline1DInterpolant &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSpline1DInterpolantShell::~CSpline1DInterpolantShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CSpline1DInterpolant *CSpline1DInterpolantShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| 1-dimensional spline interpolation                               |
//+------------------------------------------------------------------+
class CSpline1D
  {
private:
   //--- private methods
   static void       Spline1DGridDiffCubicInternal(double &x[],double &y[],const int n,const int boundltype,const double boundl,const int boundrtype,const double boundr,double &d[],double &a1[],double &a2[],double &a3[],double &b[],double &dt[]);
   static void       HeapSortPoints(double &x[],double &y[],const int n);
   static void       HeapSortPPoints(double &x[],double &y[],int &p[],const int n);
   static void       SolveTridiagonal(double &a[],double &cb[],double &c[],double &cd[],const int n,double &x[]);
   static void       SolveCyclicTridiagonal(double &a[],double &cb[],double &c[],double &d[],const int n,double &x[]);
   static double     DiffThreePoint(double t,const double x0,const double f0,double x1,const double f1,double x2,const double f2);
public:
   //--- constructor, destructor
                     CSpline1D(void);
                    ~CSpline1D(void);
   //--- public methods
   static void       Spline1DBuildLinear(double &cx[],double &cy[],const int n,CSpline1DInterpolant &c);
   static void       Spline1DBuildCubic(double &cx[],double &cy[],const int n,const int boundltype,const double boundl,const int boundrtype,const double boundr,CSpline1DInterpolant &c);
   static void       Spline1DGridDiffCubic(double &cx[],double &cy[],const int n,const int boundltype,const double boundl,const int boundrtype,const double boundr,double &d[]);
   static void       Spline1DGridDiff2Cubic(double &cx[],double &cy[],const int n,const int boundltype,const double boundl,const int boundrtype,const double boundr,double &d1[],double &d2[]);
   static void       Spline1DConvCubic(double &cx[],double &cy[],const int n,const int boundltype,const double boundl,const int boundrtype,const double boundr,double &cx2[],const int n2,double &y2[]);
   static void       Spline1DConvDiffCubic(double &cx[],double &cy[],const int n,const int boundltype,const double boundl,const int boundrtype,const double boundr,double &cx2[],const int n2,double &y2[],double &d2[]);
   static void       Spline1DConvDiff2Cubic(double &cx[],double &cy[],const int n,const int boundltype,const double boundl,const int boundrtype,const double boundr,double &cx2[],const int n2,double &y2[],double &d2[],double &dd2[]);
   static void       Spline1DBuildCatmullRom(double &cx[],double &cy[],const int n,const int boundtype,const double tension,CSpline1DInterpolant &c);
   static void       Spline1DBuildHermite(double &cx[],double &cy[],double &cd[],const int n,CSpline1DInterpolant &c);
   static void       Spline1DBuildAkima(double &cx[],double &cy[],const int n,CSpline1DInterpolant &c);
   static double     Spline1DCalc(CSpline1DInterpolant &c,double x);
   static void       Spline1DDiff(CSpline1DInterpolant &c,double x,double &s,double &ds,double &d2s);
   static void       Spline1DCopy(CSpline1DInterpolant &c,CSpline1DInterpolant &cc);
   static void       Spline1DUnpack(CSpline1DInterpolant &c,int &n,CMatrixDouble &tbl);
   static void       Spline1DLinTransX(CSpline1DInterpolant &c,const double a,const double b);
   static void       Spline1DLinTransY(CSpline1DInterpolant &c,const double a,const double b);
   static double     Spline1DIntegrate(CSpline1DInterpolant &c,double x);
   static void       Spline1DConvDiffInternal(double &xold[],double &yold[],double &dold[],const int n,double &x2[],const int n2,double &y[],const bool needy,double &d1[],const bool needd1,double &d2[],const bool needd2);
   static void       HeapSortDPoints(double &x[],double &y[],double &d[],const int n);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSpline1D::CSpline1D(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSpline1D::~CSpline1D(void)
  {

  }
//+------------------------------------------------------------------+
//| This subroutine builds linear spline interpolant                 |
//| INPUT PARAMETERS:                                                |
//|     X   -   spline nodes, array[0..N-1]                          |
//|     Y   -   function values, array[0..N-1]                       |
//|     N   -   points count (optional):                             |
//|             * N>=2                                               |
//|             * if given, only first N points are used to build    |
//|               spline                                             |
//|             * if not given, automatically detected from X/Y      |
//|               sizes (len(X) must be equal to len(Y))             |
//| OUTPUT PARAMETERS:                                               |
//|     C   -   spline interpolant                                   |
//| ORDER OF POINTS                                                  |
//| Subroutine automatically sorts points, so caller may pass        |
//| unsorted array.                                                  |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DBuildLinear(double &cx[],double &cy[],
                                           const int n,CSpline1DInterpolant &c)
  {
//--- create a variable
   int i=0;
//--- create arrays
   double x[];
   double y[];
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
//--- check
   if(!CAp::Assert(n>1,__FUNCTION__+": N<2!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check and sort points
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- function call
   HeapSortPoints(x,y,n);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(x,n),__FUNCTION__+": at least two consequent points are too close!"))
      return;
//--- Build
   c.m_periodic=false;
   c.m_n=n;
   c.m_k=3;
//--- allocation
   ArrayResizeAL(c.m_x,n);
   ArrayResizeAL(c.m_c,4*(n-1));
//--- copy
   for(i=0;i<=n-1;i++)
      c.m_x[i]=x[i];
//--- calculation
   for(i=0;i<=n-2;i++)
     {
      c.m_c[4*i+0]=y[i];
      c.m_c[4*i+1]=(y[i+1]-y[i])/(x[i+1]-x[i]);
      c.m_c[4*i+2]=0;
      c.m_c[4*i+3]=0;
     }
  }
//+------------------------------------------------------------------+
//| This subroutine builds cubic spline interpolant.                 |
//| INPUT PARAMETERS:                                                |
//|     X           -   spline nodes, array[0..N-1].                 |
//|     Y           -   function values, array[0..N-1].              |
//| OPTIONAL PARAMETERS:                                             |
//|     N           -   points count:                                |
//|                     * N>=2                                       |
//|                     * if given, only first N points are used to  |
//|                       build spline                               |
//|                     * if not given, automatically detected from  |
//|                       X/Y sizes (len(X) must be equal to len(Y)) |
//|     BoundLType  -   boundary condition type for the left boundary|
//|     BoundL      -   left boundary condition (first or second     |
//|                     derivative, depending on the BoundLType)     |
//|     BoundRType  -   boundary condition type for the right        |
//|                     boundary                                     |
//|     BoundR      -   right boundary condition (first or second    |
//|                     derivative, depending on the BoundRType)     |
//| OUTPUT PARAMETERS:                                               |
//|     C           -   spline interpolant                           |
//| ORDER OF POINTS                                                  |
//| Subroutine automatically sorts points, so caller may pass        |
//| unsorted array.                                                  |
//| SETTING BOUNDARY VALUES:                                         |
//| The BoundLType/BoundRType parameters can have the following      |
//| values:                                                          |
//|     * -1, which corresonds to the periodic (cyclic) boundary     |
//|           conditions. In this case:                              |
//|           * both BoundLType and BoundRType must be equal to -1.  |
//|           * BoundL/BoundR are ignored                            |
//|           * Y[last] is ignored (it is assumed to be equal to     |
//|             Y[first]).                                           |
//|     *  0, which  corresponds to the parabolically terminated     |
//|           spline (BoundL and/or BoundR are ignored).             |
//|     *  1, which corresponds to the first derivative boundary     |
//|           condition                                              |
//|     *  2, which corresponds to the second derivative boundary    |
//|           condition                                              |
//|     *  by default, BoundType=0 is used                           |
//| PROBLEMS WITH PERIODIC BOUNDARY CONDITIONS:                      |
//| Problems with periodic boundary conditions have                  |
//| Y[first_point]=Y[last_point]. However, this subroutine doesn't   |
//| require you to specify equal values for the first and last       |
//| points - it automatically forces them to be equal by copying     |
//| Y[first_point] (corresponds to the leftmost, minimal X[]) to     |
//| Y[last_point]. However it is recommended to pass consistent      |
//| values of Y[], i.e. to make Y[first_point]=Y[last_point].        |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DBuildCubic(double &cx[],double &cy[],
                                          const int n,const int boundltype,
                                          const double boundl,const int boundrtype,
                                          const double boundr,CSpline1DInterpolant &c)
  {
//--- create a variable
   int ylen=0;
//--- create arrays
   double a1[];
   double a2[];
   double a3[];
   double b[];
   double dt[];
   double d[];
   int    p[];
   double x[];
   double y[];
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
//--- check correctness of boundary conditions
   if(!CAp::Assert(((boundltype==-1 || boundltype==0) || boundltype==1) || boundltype==2,__FUNCTION__+": incorrect BoundLType!"))
      return;
//--- check
   if(!CAp::Assert(((boundrtype==-1 || boundrtype==0) || boundrtype==1) || boundrtype==2,__FUNCTION__+": incorrect BoundRType!"))
      return;
//--- check
   if(!CAp::Assert((boundrtype==-1 && boundltype==-1) || (boundrtype!=-1 && boundltype!=-1),__FUNCTION__+": incorrect BoundLType/BoundRType!"))
      return;
//--- check
   if(boundltype==1 || boundltype==2)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(boundl),__FUNCTION__+": BoundL is infinite or NAN!"))
         return;
     }
//--- check
   if(boundrtype==1 || boundrtype==2)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(boundr),__FUNCTION__+": BoundR is infinite or NAN!"))
         return;
     }
//--- check lengths of arguments
   if(!CAp::Assert(n>=2,__FUNCTION__+": N<2!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check and sort points
   ylen=n;
//--- check
   if(boundltype==-1)
      ylen=n-1;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,ylen),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- function call
   HeapSortPPoints(x,y,p,n);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(x,n),__FUNCTION__+": at least two consequent points are too close!"))
      return;
//--- Now we've checked and preordered everything,
//--- so we can call internal function to calculate derivatives,
//--- and then build Hermite spline using these derivatives
   Spline1DGridDiffCubicInternal(x,y,n,boundltype,boundl,boundrtype,boundr,d,a1,a2,a3,b,dt);
   Spline1DBuildHermite(x,y,d,n,c);
//--- check
   if(boundltype==-1 || boundrtype==-1)
      c.m_periodic=1;
   else
      c.m_periodic=0;
  }
//+------------------------------------------------------------------+
//| This function solves following problem: given table y[] of       |
//| function values at nodes x[], it calculates and returns table of |
//| function derivatives  d[] (calculated at the same nodes x[]).    |
//| This function yields same result as Spline1DBuildCubic() call    |
//| followed  by sequence of Spline1DDiff() calls, but it can be     |
//| several times faster when called for ordered X[] and X2[].       |
//| INPUT PARAMETERS:                                                |
//|     X           -   spline nodes                                 |
//|     Y           -   function values                              |
//| OPTIONAL PARAMETERS:                                             |
//|     N           -   points count:                                |
//|                     * N>=2                                       |
//|                     * if given, only first N points are used     |
//|                     * if not given, automatically detected from  |
//|                       X/Y sizes (len(X) must be equal to len(Y)) |
//|     BoundLType  -   boundary condition type for the left boundary|
//|     BoundL      -   left boundary condition (first or second     |
//|                     derivative, depending on the BoundLType)     |
//|     BoundRType  -   boundary condition type for the right        |
//|                     boundary                                     |
//|     BoundR      -   right boundary condition (first or second    |
//|                     derivative, depending on the BoundRType)     |
//| OUTPUT PARAMETERS:                                               |
//|     D           -   derivative values at X[]                     |
//| ORDER OF POINTS                                                  |
//| Subroutine automatically sorts points, so caller may pass        |
//| unsorted array. Derivative values are correctly reordered on     |
//| return, so D[I] is always equal to S'(X[I]) independently of     |
//| points order.                                                    |
//| SETTING BOUNDARY VALUES:                                         |
//| The BoundLType/BoundRType parameters can have the following      |
//| values:                                                          |
//|     * -1, which corresonds to the periodic (cyclic) boundary     |
//|           conditions. In this case:                              |
//|           * both BoundLType and BoundRType must be equal to -1.  |
//|           * BoundL/BoundR are ignored                            |
//|           * Y[last] is ignored (it is assumed to be equal to     |
//|           Y[first]).                                             |
//|     *  0, which corresponds to the parabolically terminated      |
//|           spline (BoundL and/or BoundR are ignored).             |
//|     *  1, which corresponds to the first derivative boundary     |
//|           condition                                              |
//|     *  2, which corresponds to the second derivative boundary    |
//|           condition                                              |
//|     *  by default, BoundType=0 is used                           |
//| PROBLEMS WITH PERIODIC BOUNDARY CONDITIONS:                      |
//| Problems with periodic boundary conditions have                  |
//| Y[first_point]=Y[last_point]. However, this subroutine doesn't   |
//| require you to specify equal values for the first and last       |
//| points - it automatically forces them to be equal by copying     |
//| Y[first_point] (corresponds to the leftmost, minimal X[]) to     |
//| Y[last_point]. However it is recommended to pass consistent      |
//| values of Y[], i.e. to make Y[first_point]=Y[last_point].        |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DGridDiffCubic(double &cx[],double &cy[],
                                             const int n,const int boundltype,
                                             const double boundl,const int boundrtype,
                                             const double boundr,double &d[])
  {
//--- create variables
   int i=0;
   int ylen=0;
   int i_=0;
//--- create arrays
   double a1[];
   double a2[];
   double a3[];
   double b[];
   double dt[];
   int    p[];
   double x[];
   double y[];
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
//--- check correctness of boundary conditions
   if(!CAp::Assert(((boundltype==-1 || boundltype==0) || boundltype==1) || boundltype==2,__FUNCTION__+": incorrect BoundLType!"))
      return;
//--- check
   if(!CAp::Assert(((boundrtype==-1 || boundrtype==0) || boundrtype==1) || boundrtype==2,__FUNCTION__+": incorrect BoundRType!"))
      return;
//--- check
   if(!CAp::Assert((boundrtype==-1 && boundltype==-1) || (boundrtype!=-1 && boundltype!=-1),__FUNCTION__+": incorrect BoundLType/BoundRType!"))
      return;
//--- check
   if(boundltype==1 || boundltype==2)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(boundl),__FUNCTION__+": BoundL is infinite or NAN!"))
         return;
     }
//--- check
   if(boundrtype==1 || boundrtype==2)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(boundr),__FUNCTION__+": BoundR is infinite or NAN!"))
         return;
     }
//--- check lengths of arguments
   if(!CAp::Assert(n>=2,__FUNCTION__+": N<2!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check and sort points
   ylen=n;
//--- check
   if(boundltype==-1)
      ylen=n-1;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,ylen),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- function call
   HeapSortPPoints(x,y,p,n);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(x,n),__FUNCTION__+": at least two consequent points are too close!"))
      return;
//--- Now we've checked and preordered everything,
//--- so we can call internal function.
   Spline1DGridDiffCubicInternal(x,y,n,boundltype,boundl,boundrtype,boundr,d,a1,a2,a3,b,dt);
//--- Remember that HeapSortPPoints() call?
//--- Now we have to reorder them back.
   if(CAp::Len(dt)<n)
      ArrayResizeAL(dt,n);
//--- copy
   for(i=0;i<=n-1;i++)
      dt[p[i]]=d[i];
   for(i_=0;i_<=n-1;i_++)
      d[i_]=dt[i_];
  }
//+------------------------------------------------------------------+
//| This function solves following problem: given table y[] of       |
//| function values at nodes x[], it calculates and returns tables of|
//| first and second function derivatives d1[] and d2[] (calculated  |
//| at the same nodes x[]).                                          |
//| This function yields same result as Spline1DBuildCubic() call    |
//| followed  by sequence of Spline1DDiff() calls, but it can be     |
//| several times faster when called for ordered X[] and X2[].       |
//| INPUT PARAMETERS:                                                |
//|     X           -   spline nodes                                 |
//|     Y           -   function values                              |
//| OPTIONAL PARAMETERS:                                             |
//|     N           -   points count:                                |
//|                     * N>=2                                       |
//|                     * if given, only first N points are used     |
//|                     * if not given, automatically detected from  |
//|                       X/Y sizes (len(X) must be equal to len(Y)) |
//|     BoundLType  -   boundary condition type for the left boundary|
//|     BoundL      -   left boundary condition (first or second     |
//|                     derivative, depending on the BoundLType)     |
//|     BoundRType  -   boundary condition type for the right        |
//|                     boundary                                     |
//|     BoundR      -   right boundary condition (first or second    |
//|                     derivative, depending on the BoundRType)     |
//| OUTPUT PARAMETERS:                                               |
//|     D1          -   S' values at X[]                             |
//|     D2          -   S'' values at X[]                            |
//| ORDER OF POINTS                                                  |
//| Subroutine automatically sorts points, so caller may pass        |
//| unsorted array. Derivative values are correctly reordered on     |
//| return, so D[I] is always equal to S'(X[I]) independently of     |
//| points order.                                                    |
//| SETTING BOUNDARY VALUES:                                         |
//| The BoundLType/BoundRType parameters can have the following      |
//| values:                                                          |
//|     * -1, which corresonds to the periodic (cyclic) boundary     |
//|           conditions. In this case:                              |
//|           * both BoundLType and BoundRType must be equal to -1.  |
//|           * BoundL/BoundR are ignored                            |
//|           * Y[last] is ignored (it is assumed to be equal to     |
//|           Y[first]).                                             |
//|     *  0, which corresponds to the parabolically terminated      |
//|           spline (BoundL and/or BoundR are ignored).             |
//|     *  1, which corresponds to the first derivative boundary     |
//|           condition                                              |
//|     *  2, which corresponds to the second derivative boundary    |
//|           condition                                              |
//|     *  by default, BoundType=0 is used                           |
//| PROBLEMS WITH PERIODIC BOUNDARY CONDITIONS:                      |
//| Problems with periodic boundary conditions have                  |
//| Y[first_point]=Y[last_point].                                    |
//| However, this subroutine doesn't require you to specify equal    |
//| values for the first and last points - it automatically forces   |
//| them to be equal by copying Y[first_point] (corresponds to the   |
//| leftmost, minimal X[]) to Y[last_point]. However it is           |
//| recommended to pass consistent values of Y[], i.e. to make       |
//| Y[first_point]=Y[last_point].                                    |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DGridDiff2Cubic(double &cx[],double &cy[],
                                              const int n,const int boundltype,
                                              const double boundl,const int boundrtype,
                                              const double boundr,double &d1[],double &d2[])
  {
//--- create variables
   int    i=0;
   int    ylen=0;
   double delta=0;
   double delta2=0;
   double delta3=0;
   double s0=0;
   double s1=0;
   double s2=0;
   double s3=0;
   int    i_=0;
//--- create arrays
   double a1[];
   double a2[];
   double a3[];
   double b[];
   double dt[];
   int    p[];
   double x[];
   double y[];
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
//--- check correctness of boundary conditions
   if(!CAp::Assert(((boundltype==-1 || boundltype==0) || boundltype==1) || boundltype==2,__FUNCTION__+": incorrect BoundLType!"))
      return;
//--- check
   if(!CAp::Assert(((boundrtype==-1 || boundrtype==0) || boundrtype==1) || boundrtype==2,__FUNCTION__+": incorrect BoundRType!"))
      return;
//--- check
   if(!CAp::Assert((boundrtype==-1 && boundltype==-1) || (boundrtype!=-1 && boundltype!=-1),__FUNCTION__+": incorrect BoundLType/BoundRType!"))
      return;
//--- check
   if(boundltype==1 || boundltype==2)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(boundl),__FUNCTION__+": BoundL is infinite or NAN!"))
         return;
     }
//--- check
   if(boundrtype==1 || boundrtype==2)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(boundr),__FUNCTION__+": BoundR is infinite or NAN!"))
         return;
     }
//--- check lengths of arguments
   if(!CAp::Assert(n>=2,__FUNCTION__+": N<2!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check and sort points
   ylen=n;
//--- check
   if(boundltype==-1)
      ylen=n-1;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,ylen),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- function call
   HeapSortPPoints(x,y,p,n);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(x,n),__FUNCTION__+": at least two consequent points are too close!"))
      return;
//--- Now we've checked and preordered everything,
//--- so we can call internal function.
//--- After this call we will calculate second derivatives
//--- (manually,by converting to the power basis)
   Spline1DGridDiffCubicInternal(x,y,n,boundltype,boundl,boundrtype,boundr,d1,a1,a2,a3,b,dt);
//--- allocation
   ArrayResizeAL(d2,n);
   delta=0;
   s2=0;
   s3=0;
//--- calculation
   for(i=0;i<=n-2;i++)
     {
      //--- We convert from Hermite basis to the power basis.
      //--- Si is coefficient before x^i.
      //--- Inside this cycle we need just S2,
      //--- because we calculate S'' exactly at spline node,
      //--- (only x^2 matters at x=0),but after iterations
      //--- will be over,we will need other coefficients
      //--- to calculate spline value at the last node.
      delta=x[i+1]-x[i];
      delta2=CMath::Sqr(delta);
      delta3=delta*delta2;
      s0=y[i];
      s1=d1[i];
      s2=(3*(y[i+1]-y[i])-2*d1[i]*delta-d1[i+1]*delta)/delta2;
      s3=(2*(y[i]-y[i+1])+d1[i]*delta+d1[i+1]*delta)/delta3;
      d2[i]=2*s2;
     }
   d2[n-1]=2*s2+6*s3*delta;
//--- Remember that HeapSortPPoints() call?
//--- Now we have to reorder them back.
   if(CAp::Len(dt)<n)
      ArrayResizeAL(dt,n);
//--- copy
   for(i=0;i<=n-1;i++)
      dt[p[i]]=d1[i];
   for(i_=0;i_<=n-1;i_++)
      d1[i_]=dt[i_];
   for(i=0;i<=n-1;i++)
      dt[p[i]]=d2[i];
   for(i_=0;i_<=n-1;i_++)
      d2[i_]=dt[i_];
  }
//+------------------------------------------------------------------+
//| This function solves following problem: given table y[] of       |
//| function values at old nodes x[]  and new nodes  x2[], it        |
//| calculates and returns table of function values y2[] (calculated |
//| at x2[]).                                                        |
//| This function yields same result as Spline1DBuildCubic() call    |
//| followed  by sequence of Spline1DDiff() calls, but it can be     |
//| several times faster when called for ordered X[] and X2[].       |
//| INPUT PARAMETERS:                                                |
//|     X           -   old spline nodes                             |
//|     Y           -   function values                              |
//|     X2           -  new spline nodes                             |
//| OPTIONAL PARAMETERS:                                             |
//|     N           -   points count:                                |
//|                     * N>=2                                       |
//|                     * if given, only first N points from X/Y are |
//|                       used                                       |
//|                     * if not given, automatically detected from  |
//|                       X/Y sizes (len(X) must be equal to len(Y)) |
//|     BoundLType  -   boundary condition type for the left boundary|
//|     BoundL      -   left boundary condition (first or second     |
//|                     derivative, depending on the BoundLType)     |
//|     BoundRType  -   boundary condition type for the right        |
//|                     boundary                                     |
//|     BoundR      -   right boundary condition (first or second    |
//|                     derivative, depending on the BoundRType)     |
//|     N2          -   new points count:                            |
//|                     * N2>=2                                      |
//|                     * if given, only first N2 points from X2 are |
//|                       used                                       |
//|                     * if not given, automatically detected from  |
//|                       X2 size                                    |
//| OUTPUT PARAMETERS:                                               |
//|     F2          -   function values at X2[]                      |
//| ORDER OF POINTS                                                  |
//| Subroutine automatically sorts points, so caller  may pass       |
//| unsorted array. Function  values  are correctly reordered on     |
//| return, so F2[I] is always equal to S(X2[I]) independently of    |
//| points order.                                                    |
//| SETTING BOUNDARY VALUES:                                         |
//| The BoundLType/BoundRType parameters can have the following      |
//| values:                                                          |
//|     * -1, which corresonds to the periodic (cyclic) boundary     |
//|           conditions. In this case:                              |
//|           * both BoundLType and BoundRType must be equal to -1.  |
//|           * BoundL/BoundR are ignored                            |
//|           * Y[last] is ignored (it is assumed to be equal to     |
//|           Y[first]).                                             |
//|     *  0, which corresponds to the parabolically terminated      |
//|           spline (BoundL and/or BoundR are ignored).             |
//|     *  1, which corresponds to the first derivative boundary     |
//|           condition                                              |
//|     *  2, which corresponds to the second derivative boundary    |
//|           condition                                              |
//|     *  by default, BoundType=0 is used                           |
//| PROBLEMS WITH PERIODIC BOUNDARY CONDITIONS:                      |
//| Problems with periodic boundary conditions have                  |
//| Y[first_point]=Y[last_point]. However, this subroutine doesn't   |
//| require you to specify equal values for the first and last       |
//| points - it automatically forces them to be equal by copying     |
//| Y[first_point] (corresponds to the leftmost, minimal X[]) to     |
//| Y[last_point]. However it is recommended to pass consistent      |
//| values of Y[], i.e. to make Y[first_point]=Y[last_point].        |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DConvCubic(double &cx[],double &cy[],
                                         const int n,const int boundltype,
                                         const double boundl,const int boundrtype,
                                         const double boundr,double &cx2[],
                                         const int n2,double &y2[])
  {
//--- create variables
   int    i=0;
   int    ylen=0;
   double t=0;
   double t2=0;
   int    i_=0;
//--- create arrays
   double a1[];
   double a2[];
   double a3[];
   double b[];
   double d[];
   double dt[];
   double d1[];
   double d2[];
   int    p[];
   int    p2[];
   double x[];
   double y[];
   double x2[];
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
   ArrayCopy(x2,cx2);
//--- check correctness of boundary conditions
   if(!CAp::Assert(((boundltype==-1 || boundltype==0) || boundltype==1) || boundltype==2,__FUNCTION__+": incorrect BoundLType!"))
      return;
//--- check
   if(!CAp::Assert(((boundrtype==-1 || boundrtype==0) || boundrtype==1) || boundrtype==2,__FUNCTION__+": incorrect BoundRType!"))
      return;
//--- check
   if(!CAp::Assert((boundrtype==-1 && boundltype==-1) || (boundrtype!=-1 && boundltype!=-1),__FUNCTION__+": incorrect BoundLType/BoundRType!"))
      return;
//--- check
   if(boundltype==1 || boundltype==2)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(boundl),__FUNCTION__+": BoundL is infinite or NAN!"))
         return;
     }
//--- check
   if(boundrtype==1 || boundrtype==2)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(boundr),__FUNCTION__+": BoundR is infinite or NAN!"))
         return;
     }
//--- check lengths of arguments
   if(!CAp::Assert(n>=2,__FUNCTION__+": N<2!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(n2>=2,__FUNCTION__+": N2<2!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x2)>=n2,__FUNCTION__+": Length(X2)<N2!"))
      return;
//--- check and sort X/Y
   ylen=n;
//--- check
   if(boundltype==-1)
      ylen=n-1;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,ylen),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x2,n2),__FUNCTION__+": X2 contains infinite or NAN values!"))
      return;
//--- function call
   HeapSortPPoints(x,y,p,n);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(x,n),__FUNCTION__+": at least two consequent points are too close!"))
      return;
//--- set up DT (we will need it below)
   ArrayResizeAL(dt,MathMax(n,n2));
//--- sort X2:
//--- * use fake array DT because HeapSortPPoints() needs both integer AND real arrays
//--- * if we have periodic problem,wrap points
//--- * sort them,store permutation at P2
   if(boundrtype==-1 && boundltype==-1)
     {
      for(i=0;i<=n2-1;i++)
        {
         t=x2[i];
         CApServ::ApPeriodicMap(t,x[0],x[n-1],t2);
         x2[i]=t;
        }
     }
//--- function call
   HeapSortPPoints(x2,dt,p2,n2);
//--- Now we've checked and preordered everything,so we:
//--- * call internal GridDiff() function to get Hermite form of spline
//--- * convert using internal Conv() function
//--- * convert Y2 back to original order
   Spline1DGridDiffCubicInternal(x,y,n,boundltype,boundl,boundrtype,boundr,d,a1,a2,a3,b,dt);
   Spline1DConvDiffInternal(x,y,d,n,x2,n2,y2,true,d1,false,d2,false);
//--- check
   if(!CAp::Assert(CAp::Len(dt)>=n2,__FUNCTION__+": internal error!"))
      return;
//--- copy
   for(i=0;i<=n2-1;i++)
      dt[p2[i]]=y2[i];
   for(i_=0;i_<=n2-1;i_++)
      y2[i_]=dt[i_];
  }
//+------------------------------------------------------------------+
//| This function solves following problem: given table y[] of       |
//| function values at old nodes x[] and new nodes x2[], it          |
//| calculates and returns table of function values y2[] and         |
//| derivatives d2[] (calculated at x2[]).                           |
//| This function yields same result as Spline1DBuildCubic() call    |
//| followed by sequence of Spline1DDiff() calls, but it can be      |
//| several times faster when called for ordered X[] and X2[].       |
//| INPUT PARAMETERS:                                                |
//|     X           -   old spline nodes                             |
//|     Y           -   function values                              |
//|     X2           -  new spline nodes                             |
//| OPTIONAL PARAMETERS:                                             |
//|     N           -   points count:                                |
//|                     * N>=2                                       |
//|                     * if given, only first N points from X/Y are |
//|                       used                                       |
//|                     * if not given, automatically detected from  |
//|                       X/Y sizes (len(X) must be equal to len(Y)) |
//|     BoundLType  -   boundary condition type for the left boundary|
//|     BoundL      -   left boundary condition (first or second     |
//|                     derivative, depending on the BoundLType)     |
//|     BoundRType  -   boundary condition type for the right        |
//|                     boundary                                     |
//|     BoundR      -   right boundary condition (first or second    |
//|                     derivative, depending on the BoundRType)     |
//|     N2          -   new points count:                            |
//|                     * N2>=2                                      |
//|                     * if given, only first N2 points from X2 are |
//|                       used                                       |
//|                     * if not given, automatically detected from  |
//|                       X2 size                                    |
//| OUTPUT PARAMETERS:                                               |
//|     F2          -   function values at X2[]                      |
//|     D2          -   first derivatives at X2[]                    |
//| ORDER OF POINTS                                                  |
//| Subroutine automatically sorts points, so caller may pass        |
//| unsorted array. Function  values  are correctly reordered on     |
//| return, so F2[I] is always equal to S(X2[I]) independently of    |
//| points order.                                                    |
//| SETTING BOUNDARY VALUES:                                         |
//| The BoundLType/BoundRType parameters can have the following      |
//| values:                                                          |
//|     * -1, which corresonds to the periodic (cyclic) boundary     |
//|           conditions. In this case:                              |
//|           * both BoundLType and BoundRType must be equal to -1.  |
//|           * BoundL/BoundR are ignored                            |
//|           * Y[last] is ignored (it is assumed to be equal to     |
//|             Y[first]).                                           |
//|     *  0, which corresponds to the parabolically terminated      |
//|           spline (BoundL and/or BoundR are ignored).             |
//|     *  1, which corresponds to the first derivative boundary     |
//|           condition                                              |
//|     *  2, which corresponds to the second derivative boundary    |
//|           condition                                              |
//|     *  by default, BoundType=0 is used                           |
//| PROBLEMS WITH PERIODIC BOUNDARY CONDITIONS:                      |
//| Problems with periodic boundary conditions have                  |
//| Y[first_point]=Y[last_point]. However, this subroutine doesn't   |
//| require you to specify equal values for the first and last       |
//| points - it automatically forces them to be equal by copying     |
//| Y[first_point] (corresponds to the leftmost, minimal X[]) to     |
//| Y[last_point]. However it is recommended to pass consistent      |
//| values of Y[], i.e. to make Y[first_point]=Y[last_point].        |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DConvDiffCubic(double &cx[],double &cy[],
                                             const int n,const int boundltype,
                                             const double boundl,const int boundrtype,
                                             const double boundr,double &cx2[],
                                             const int n2,double &y2[],double &d2[])
  {
//--- create variables
   int    i=0;
   int    ylen=0;
   double t=0;
   double t2=0;
   int    i_=0;
//--- create arrays
   double a1[];
   double a2[];
   double a3[];
   double b[];
   double d[];
   double dt[];
   double rt1[];
   int    p[];
   int    p2[];
   double x[];
   double y[];
   double x2[];
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
   ArrayCopy(x2,cx2);
//--- check correctness of boundary conditions
   if(!CAp::Assert(((boundltype==-1 || boundltype==0) || boundltype==1) || boundltype==2,__FUNCTION__+": incorrect BoundLType!"))
      return;
//--- check
   if(!CAp::Assert(((boundrtype==-1 || boundrtype==0) || boundrtype==1) || boundrtype==2,__FUNCTION__+": incorrect BoundRType!"))
      return;
//--- check
   if(!CAp::Assert((boundrtype==-1 && boundltype==-1) || (boundrtype!=-1 && boundltype!=-1),__FUNCTION__+": incorrect BoundLType/BoundRType!"))
      return;
//--- check
   if(boundltype==1 || boundltype==2)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(boundl),__FUNCTION__+": BoundL is infinite or NAN!"))
         return;
     }
//--- check
   if(boundrtype==1 || boundrtype==2)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(boundr),__FUNCTION__+": BoundR is infinite or NAN!"))
         return;
     }
//--- check lengths of arguments
   if(!CAp::Assert(n>=2,__FUNCTION__+": N<2!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(n2>=2,__FUNCTION__+"Spline1DConvDiffCubic: N2<2!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x2)>=n2,__FUNCTION__+": Length(X2)<N2!"))
      return;
//--- check and sort X/Y
   ylen=n;
//--- check
   if(boundltype==-1)
      ylen=n-1;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,ylen),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x2,n2),__FUNCTION__+": X2 contains infinite or NAN values!"))
      return;
//--- function call
   HeapSortPPoints(x,y,p,n);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(x,n),__FUNCTION__+": at least two consequent points are too close!"))
      return;
//--- set up DT (we will need it below)
   ArrayResizeAL(dt,MathMax(n,n2));
//--- sort X2:
//--- * use fake array DT because HeapSortPPoints() needs both integer AND real arrays
//--- * if we have periodic problem,wrap points
//--- * sort them,store permutation at P2
   if(boundrtype==-1 && boundltype==-1)
     {
      for(i=0;i<=n2-1;i++)
        {
         t=x2[i];
         CApServ::ApPeriodicMap(t,x[0],x[n-1],t2);
         x2[i]=t;
        }
     }
//--- function call
   HeapSortPPoints(x2,dt,p2,n2);
//--- Now we've checked and preordered everything,so we:
//--- * call internal GridDiff() function to get Hermite form of spline
//--- * convert using internal Conv() function
//--- * convert Y2 back to original order
   Spline1DGridDiffCubicInternal(x,y,n,boundltype,boundl,boundrtype,boundr,d,a1,a2,a3,b,dt);
   Spline1DConvDiffInternal(x,y,d,n,x2,n2,y2,true,d2,true,rt1,false);
//--- check
   if(!CAp::Assert(CAp::Len(dt)>=n2,__FUNCTION__+": internal error!"))
      return;
//--- copy
   for(i=0;i<=n2-1;i++)
      dt[p2[i]]=y2[i];
   for(i_=0;i_<=n2-1;i_++)
      y2[i_]=dt[i_];
   for(i=0;i<=n2-1;i++)
      dt[p2[i]]=d2[i];
   for(i_=0;i_<=n2-1;i_++)
      d2[i_]=dt[i_];
  }
//+------------------------------------------------------------------+
//| This function solves following problem: given table y[] of       |
//| function values at old nodes x[] and new nodes  x2[], it         |
//| calculates and returns table of function values y2[], first and  |
//| second derivatives d2[] and dd2[] (calculated at x2[]).          |
//| This function yields same result as Spline1DBuildCubic() call    |
//| followed  by sequence of Spline1DDiff() calls, but it can be     |
//| several times faster when called for ordered X[] and X2[].       |
//| INPUT PARAMETERS:                                                |
//|     X           -   old spline nodes                             |
//|     Y           -   function values                              |
//|     X2           -  new spline nodes                             |
//| OPTIONAL PARAMETERS:                                             |
//|     N           -   points count:                                |
//|                     * N>=2                                       |
//|                     * if given, only first N points from X/Y are |
//|                       used                                       |
//|                     * if not given, automatically detected from  |
//|                       X/Y sizes (len(X) must be equal to len(Y)) |
//|     BoundLType  -   boundary condition type for the left boundary|
//|     BoundL      -   left boundary condition (first or second     |
//|                     derivative, depending on the BoundLType)     |
//|     BoundRType  -   boundary condition type for the right        |
//|                     boundary                                     |
//|     BoundR      -   right boundary condition (first or second    |
//|                     derivative, depending on the BoundRType)     |
//|     N2          -   new points count:                            |
//|                     * N2>=2                                      |
//|                     * if given, only first N2 points from X2 are |
//|                       used                                       |
//|                     * if not given, automatically detected from  |
//|                       X2 size                                    |
//| OUTPUT PARAMETERS:                                               |
//|     F2          -   function values at X2[]                      |
//|     D2          -   first derivatives at X2[]                    |
//|     DD2         -   second derivatives at X2[]                   |
//| ORDER OF POINTS                                                  |
//| Subroutine automatically sorts points, so caller may pass        |
//| unsorted array. Function values are correctly reordered on       |
//| return, so F2[I] is always equal to S(X2[I]) independently of    |
//| points order.                                                    |
//| SETTING BOUNDARY VALUES:                                         |
//| The BoundLType/BoundRType parameters can have the following      |
//| values:                                                          |
//|     * -1, which corresonds to the periodic (cyclic) boundary     |
//|           conditions. In this case:                              |
//|           * both BoundLType and BoundRType must be equal to -1.  |
//|           * BoundL/BoundR are ignored                            |
//|           * Y[last] is ignored (it is assumed to be equal to     |
//|             Y[first]).                                           |
//|     *  0, which corresponds to the parabolically terminated      |
//|           spline (BoundL and/or BoundR are ignored).             |
//|     *  1, which corresponds to the first derivative boundary     |
//|           condition                                              |
//|     *  2, which corresponds to the second derivative boundary    |
//|           condition                                              |
//|     *  by default, BoundType=0 is used                           |
//| PROBLEMS WITH PERIODIC BOUNDARY CONDITIONS:                      |
//| Problems with periodic boundary conditions have                  |
//| Y[first_point]=Y[last_point]. However, this subroutine doesn't   |
//| require you to specify equal values for the first and last       |
//| points - it automatically forces them to be equal by copying     |
//| Y[first_point] (corresponds to the leftmost, minimal X[]) to     |
//| Y[last_point]. However it is recommended to pass consistent      |
//| values of Y[], i.e. to make Y[first_point]=Y[last_point].        |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DConvDiff2Cubic(double &cx[],double &cy[],
                                              const int n,const int boundltype,
                                              const double boundl,
                                              const int boundrtype,
                                              const double boundr,double &cx2[],
                                              const int n2,double &y2[],
                                              double &d2[],double &dd2[])
  {
//--- create variables
   int    i=0;
   int    ylen=0;
   double t=0;
   double t2=0;
   int    i_=0;
//--- create arrays
   double a1[];
   double a2[];
   double a3[];
   double b[];
   double d[];
   double dt[];
   int    p[];
   int    p2[];
   double x[];
   double y[];
   double x2[];
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
   ArrayCopy(x2,cx2);
//--- check correctness of boundary conditions
   if(!CAp::Assert(((boundltype==-1 || boundltype==0) || boundltype==1) || boundltype==2,__FUNCTION__+": incorrect BoundLType!"))
      return;
//--- check
   if(!CAp::Assert(((boundrtype==-1 || boundrtype==0) || boundrtype==1) || boundrtype==2,__FUNCTION__+": incorrect BoundRType!"))
      return;
//--- check
   if(!CAp::Assert((boundrtype==-1 && boundltype==-1) || (boundrtype!=-1 && boundltype!=-1),__FUNCTION__+": incorrect BoundLType/BoundRType!"))
      return;
//--- check
   if(boundltype==1 || boundltype==2)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(boundl),__FUNCTION__+": BoundL is infinite or NAN!"))
         return;
     }
//--- check
   if(boundrtype==1 || boundrtype==2)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(boundr),__FUNCTION__+": BoundR is infinite or NAN!"))
         return;
     }
//--- check lengths of arguments
   if(!CAp::Assert(n>=2,__FUNCTION__+": N<2!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(n2>=2,__FUNCTION__+": N2<2!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x2)>=n2,__FUNCTION__+": Length(X2)<N2!"))
      return;
//--- check and sort X/Y
   ylen=n;
//--- check
   if(boundltype==-1)
      ylen=n-1;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,ylen),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x2,n2),__FUNCTION__+": X2 contains infinite or NAN values!"))
      return;
//--- function call
   HeapSortPPoints(x,y,p,n);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(x,n),__FUNCTION__+": at least two consequent points are too close!"))
      return;
//--- set up DT (we will need it below)
   ArrayResizeAL(dt,MathMax(n,n2));
//--- sort X2:
//--- * use fake array DT because HeapSortPPoints() needs both integer AND real arrays
//--- * if we have periodic problem,wrap points
//--- * sort them,store permutation at P2
   if(boundrtype==-1 && boundltype==-1)
     {
      for(i=0;i<=n2-1;i++)
        {
         t=x2[i];
         CApServ::ApPeriodicMap(t,x[0],x[n-1],t2);
         x2[i]=t;
        }
     }
//--- function call
   HeapSortPPoints(x2,dt,p2,n2);
//--- Now we've checked and preordered everything,so we:
//--- * call internal GridDiff() function to get Hermite form of spline
//--- * convert using internal Conv() function
//--- * convert Y2 back to original order
   Spline1DGridDiffCubicInternal(x,y,n,boundltype,boundl,boundrtype,boundr,d,a1,a2,a3,b,dt);
   Spline1DConvDiffInternal(x,y,d,n,x2,n2,y2,true,d2,true,dd2,true);
//--- check
   if(!CAp::Assert(CAp::Len(dt)>=n2,__FUNCTION__+": internal error!"))
      return;
//--- copy
   for(i=0;i<=n2-1;i++)
      dt[p2[i]]=y2[i];
   for(i_=0;i_<=n2-1;i_++)
      y2[i_]=dt[i_];
   for(i=0;i<=n2-1;i++)
      dt[p2[i]]=d2[i];
   for(i_=0;i_<=n2-1;i_++)
      d2[i_]=dt[i_];
   for(i=0;i<=n2-1;i++)
      dt[p2[i]]=dd2[i];
   for(i_=0;i_<=n2-1;i_++)
      dd2[i_]=dt[i_];
  }
//+------------------------------------------------------------------+
//| This subroutine builds Catmull-Rom spline interpolant.           |
//| INPUT PARAMETERS:                                                |
//|     X           -   spline nodes, array[0..N-1].                 |
//|     Y           -   function values, array[0..N-1].              |
//| OPTIONAL PARAMETERS:                                             |
//|     N           -   points count:                                |
//|                     * N>=2                                       |
//|                     * if given, only first N points are used to  |
//|                       build spline                               |
//|                     * if not given, automatically detected from  |
//|                       X/Y sizes (len(X) must be equal to len(Y)) |
//|     BoundType   -   boundary condition type:                     |
//|                     * -1 for periodic boundary condition         |
//|                     *  0 for parabolically terminated spline     |
//|                        (default)                                 |
//|     Tension     -   tension parameter:                           |
//|                     * tension=0   corresponds to classic         |
//|                       Catmull-Rom spline (default)               |
//|                     * 0<tension<1 corresponds to more general    |
//|                       form - cardinal spline                     |
//| OUTPUT PARAMETERS:                                               |
//|     C           -   spline interpolant                           |
//| ORDER OF POINTS                                                  |
//| Subroutine automatically sorts points, so caller may pass        |
//| unsorted array.                                                  |
//| PROBLEMS WITH PERIODIC BOUNDARY CONDITIONS:                      |
//| Problems with periodic boundary conditions have                  |
//| Y[first_point]=Y[last_point]. However, this subroutine doesn't   |
//| require you to specify equal values for the first and last       |
//| points - it automatically forces them to be equal by copying     |
//| Y[first_point] (corresponds to the leftmost, minimal X[]) to     |
//| Y[last_point]. However it is recommended to pass consistent      |
//| values of Y[], i.e. to make Y[first_point]=Y[last_point].        |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DBuildCatmullRom(double &cx[],double &cy[],
                                               const int n,const int boundtype,
                                               const double tension,
                                               CSpline1DInterpolant &c)
  {
//--- create a variable
   int i=0;
//--- create arrays
   double d[];
   double x[];
   double y[];
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
//--- check
   if(!CAp::Assert(n>=2,__FUNCTION__+": N<2!"))
      return;
//--- check
   if(!CAp::Assert(boundtype==-1 || boundtype==0,__FUNCTION__+": incorrect BoundType!"))
      return;
//--- check
   if(!CAp::Assert((double)(tension)>=0.0,__FUNCTION__+": Tension<0!"))
      return;
//--- check
   if(!CAp::Assert((double)(tension)<=(double)(1),__FUNCTION__+": Tension>1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check and sort points
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- function call
   HeapSortPoints(x,y,n);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(x,n),__FUNCTION__+": at least two consequent points are too close!"))
      return;
//--- Special cases:
//--- * N=2,parabolic terminated boundary condition on both ends
//--- * N=2,periodic boundary condition
   if(n==2 && boundtype==0)
     {
      //--- Just linear spline
      Spline1DBuildLinear(x,y,n,c);
      return;
     }
   if(n==2 && boundtype==-1)
     {
      //--- Same as cubic spline with periodic conditions
      Spline1DBuildCubic(x,y,n,-1,0.0,-1,0.0,c);
      return;
     }
//--- Periodic or non-periodic boundary conditions
   if(boundtype==-1)
     {
      //--- Periodic boundary conditions
      y[n-1]=y[0];
      //--- allocation
      ArrayResizeAL(d,n);
      d[0]=(y[1]-y[n-2])/(2*(x[1]-x[0]+x[n-1]-x[n-2]));
      for(i=1;i<=n-2;i++)
         d[i]=(1-tension)*(y[i+1]-y[i-1])/(x[i+1]-x[i-1]);
      d[n-1]=d[0];
      //--- Now problem is reduced to the cubic Hermite spline
      Spline1DBuildHermite(x,y,d,n,c);
      c.m_periodic=true;
     }
   else
     {
      //--- Non-periodic boundary conditions
      ArrayResizeAL(d,n);
      for(i=1;i<=n-2;i++)
         d[i]=(1-tension)*(y[i+1]-y[i-1])/(x[i+1]-x[i-1]);
      d[0]=2*(y[1]-y[0])/(x[1]-x[0])-d[1];
      d[n-1]=2*(y[n-1]-y[n-2])/(x[n-1]-x[n-2])-d[n-2];
      //--- Now problem is reduced to the cubic Hermite spline
      Spline1DBuildHermite(x,y,d,n,c);
     }
  }
//+------------------------------------------------------------------+
//| This subroutine builds Hermite spline interpolant.               |
//| INPUT PARAMETERS:                                                |
//|     X           -   spline nodes, array[0..N-1]                  |
//|     Y           -   function values, array[0..N-1]               |
//|     D           -   derivatives, array[0..N-1]                   |
//|     N           -   points count (optional):                     |
//|                     * N>=2                                       |
//|                     * if given, only first N points are used to  |
//|                       build spline                               |
//|                     * if not given, automatically detected from  |
//|                       X/Y sizes (len(X) must be equal to len(Y)) |
//| OUTPUT PARAMETERS:                                               |
//|     C           -   spline interpolant.                          |
//| ORDER OF POINTS                                                  |
//| Subroutine automatically sorts points, so caller may pass        |
//| unsorted array.                                                  |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DBuildHermite(double &cx[],double &cy[],
                                            double &cd[],const int n,
                                            CSpline1DInterpolant &c)
  {
//--- create variables
   int    i=0;
   double delta=0;
   double delta2=0;
   double delta3=0;
//--- create arrays
   double x[];
   double y[];
   double d[];
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
   ArrayCopy(d,cd);
//--- check
   if(!CAp::Assert(n>=2,__FUNCTION__+": N<2!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(d)>=n,__FUNCTION__+": Length(D)<N!"))
      return;
//--- check and sort points
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(d,n),__FUNCTION__+": D contains infinite or NAN values!"))
      return;
   HeapSortDPoints(x,y,d,n);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(x,n),__FUNCTION__+": at least two consequent points are too close!"))
      return;
//--- Build
   ArrayResizeAL(c.m_x,n);
   ArrayResizeAL(c.m_c,4*(n-1));
//--- change values
   c.m_periodic=false;
   c.m_k=3;
   c.m_n=n;
//--- copy
   for(i=0;i<=n-1;i++)
      c.m_x[i]=x[i];
//--- calculation
   for(i=0;i<=n-2;i++)
     {
      delta=x[i+1]-x[i];
      delta2=CMath::Sqr(delta);
      delta3=delta*delta2;
      c.m_c[4*i+0]=y[i];
      c.m_c[4*i+1]=d[i];
      c.m_c[4*i+2]=(3*(y[i+1]-y[i])-2*d[i]*delta-d[i+1]*delta)/delta2;
      c.m_c[4*i+3]=(2*(y[i]-y[i+1])+d[i]*delta+d[i+1]*delta)/delta3;
     }
  }
//+------------------------------------------------------------------+
//| This subroutine builds Akima spline interpolant                  |
//| INPUT PARAMETERS:                                                |
//|     X           -   spline nodes, array[0..N-1]                  |
//|     Y           -   function values, array[0..N-1]               |
//|     N           -   points count (optional):                     |
//|                     * N>=5                                       |
//|                     * if given, only first N points are used to  |
//|                       build spline                               |
//|                     * if not given, automatically detected from  |
//|                       X/Y sizes (len(X) must be equal to len(Y)) |
//| OUTPUT PARAMETERS:                                               |
//|     C           -   spline interpolant                           |
//| ORDER OF POINTS                                                  |
//| Subroutine automatically sorts points, so caller may pass        |
//| unsorted array.                                                  |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DBuildAkima(double &cx[],double &cy[],
                                          const int n,CSpline1DInterpolant &c)
  {
//--- create a variable
   int i=0;
//--- create arrays
   double d[];
   double w[];
   double diff[];
   double x[];
   double y[];
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
//--- check
   if(!CAp::Assert(n>=5,__FUNCTION__+": N<5!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check and sort points
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- function call
   HeapSortPoints(x,y,n);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(x,n),__FUNCTION__+": at least two consequent points are too close!"))
      return;
//--- Prepare W (weights),Diff (divided differences)
   ArrayResizeAL(w,n-1);
   ArrayResizeAL(diff,n-1);
   for(i=0;i<=n-2;i++)
      diff[i]=(y[i+1]-y[i])/(x[i+1]-x[i]);
   for(i=1;i<=n-2;i++)
      w[i]=MathAbs(diff[i]-diff[i-1]);
//--- Prepare Hermite interpolation scheme
   ArrayResizeAL(d,n);
   for(i=2;i<=n-3;i++)
     {
      //--- check
      if(MathAbs(w[i-1])+MathAbs(w[i+1])!=0.0)
         d[i]=(w[i+1]*diff[i-1]+w[i-1]*diff[i])/(w[i+1]+w[i-1]);
      else
         d[i]=((x[i+1]-x[i])*diff[i-1]+(x[i]-x[i-1])*diff[i])/(x[i+1]-x[i-1]);
     }
//--- change values
   d[0]=DiffThreePoint(x[0],x[0],y[0],x[1],y[1],x[2],y[2]);
   d[1]=DiffThreePoint(x[1],x[0],y[0],x[1],y[1],x[2],y[2]);
   d[n-2]=DiffThreePoint(x[n-2],x[n-3],y[n-3],x[n-2],y[n-2],x[n-1],y[n-1]);
   d[n-1]=DiffThreePoint(x[n-1],x[n-3],y[n-3],x[n-2],y[n-2],x[n-1],y[n-1]);
//--- Build Akima spline using Hermite interpolation scheme
   Spline1DBuildHermite(x,y,d,n,c);
  }
//+------------------------------------------------------------------+
//| This subroutine calculates the value of the spline at the given  |
//| point X.                                                         |
//| INPUT PARAMETERS:                                                |
//|     C   -   spline interpolant                                   |
//|     X   -   point                                                |
//| Result:                                                          |
//|     S(x)                                                         |
//+------------------------------------------------------------------+
static double CSpline1D::Spline1DCalc(CSpline1DInterpolant &c,double x)
  {
//--- create variables
   int    l=0;
   int    r=0;
   int    m=0;
   double t=0;
//--- check
   if(!CAp::Assert(c.m_k==3,__FUNCTION__+": internal error"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(!CInfOrNaN::IsInfinity(x),__FUNCTION__+": infinite X!"))
      return(EMPTY_VALUE);
//--- special case: NaN
   if(CInfOrNaN::IsNaN(x))
      return(CInfOrNaN::NaN());
//--- correct if periodic
   if(c.m_periodic)
      CApServ::ApPeriodicMap(x,c.m_x[0],c.m_x[c.m_n-1],t);
//--- Binary search in the [ x[0],...,x[n-2] ] (x[n-1] is not included)
   l=0;
   r=c.m_n-2+1;
   while(l!=r-1)
     {
      m=(l+r)/2;
      //--- check
      if(c.m_x[m]>=x)
         r=m;
      else
         l=m;
     }
//--- Interpolation
   x=x-c.m_x[l];
   m=4*l;
//--- return result
   return(c.m_c[m]+x*(c.m_c[m+1]+x*(c.m_c[m+2]+x*c.m_c[m+3])));
  }
//+------------------------------------------------------------------+
//| This subroutine differentiates the spline.                       |
//| INPUT PARAMETERS:                                                |
//|     C   -   spline interpolant.                                  |
//|     X   -   point                                                |
//| Result:                                                          |
//|     S   -   S(x)                                                 |
//|     DS  -   S'(x)                                                |
//|     D2S -   S''(x)                                               |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DDiff(CSpline1DInterpolant &c,double x,
                                    double &s,double &ds,double &d2s)
  {
//--- create variables
   int    l=0;
   int    r=0;
   int    m=0;
   double t=0;
//--- initialization
   s=0;
   ds=0;
   d2s=0;
//--- check
   if(!CAp::Assert(c.m_k==3,__FUNCTION__+": internal error"))
      return;
//--- check
   if(!CAp::Assert(!CInfOrNaN::IsInfinity(x),__FUNCTION__+": infinite X!"))
      return;
//--- special case: NaN
   if(CInfOrNaN::IsNaN(x))
     {
      //--- change values
      s=CInfOrNaN::NaN();
      ds=CInfOrNaN::NaN();
      d2s=CInfOrNaN::NaN();
      //--- exit the function
      return;
     }
//--- correct if periodic
   if(c.m_periodic)
      CApServ::ApPeriodicMap(x,c.m_x[0],c.m_x[c.m_n-1],t);
//--- Binary search
   l=0;
   r=c.m_n-2+1;
   while(l!=r-1)
     {
      m=(l+r)/2;
      //--- check
      if(c.m_x[m]>=x)
         r=m;
      else
         l=m;
     }
//--- Differentiation
   x=x-c.m_x[l];
   m=4*l;
   s=c.m_c[m]+x*(c.m_c[m+1]+x*(c.m_c[m+2]+x*c.m_c[m+3]));
   ds=c.m_c[m+1]+2*x*c.m_c[m+2]+3*CMath::Sqr(x)*c.m_c[m+3];
   d2s=2*c.m_c[m+2]+6*x*c.m_c[m+3];
  }
//+------------------------------------------------------------------+
//| This subroutine makes the copy of the spline.                    |
//| INPUT PARAMETERS:                                                |
//|     C   -   spline interpolant.                                  |
//| Result:                                                          |
//|     CC  -   spline copy                                          |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DCopy(CSpline1DInterpolant &c,CSpline1DInterpolant &cc)
  {
//--- create a variable
   int i_=0;
//--- change values
   cc.m_periodic=c.m_periodic;
   cc.m_n=c.m_n;
   cc.m_k=c.m_k;
//--- allocation
   ArrayResizeAL(cc.m_x,cc.m_n);
//--- copy
   for(i_=0;i_<=cc.m_n-1;i_++)
      cc.m_x[i_]=c.m_x[i_];
//--- allocation
   ArrayResizeAL(cc.m_c,(cc.m_k+1)*(cc.m_n-1));
//--- copy
   for(i_=0;i_<=(cc.m_k+1)*(cc.m_n-1)-1;i_++)
      cc.m_c[i_]=c.m_c[i_];
  }
//+------------------------------------------------------------------+
//| This subroutine unpacks the spline into the coefficients table.  |
//| INPUT PARAMETERS:                                                |
//|     C   -   spline interpolant.                                  |
//|     X   -   point                                                |
//| Result:                                                          |
//|     Tbl -   coefficients table, unpacked format, array[0..N-2,   |
//|             0..5].                                               |
//|             For I = 0...N-2:                                     |
//|                 Tbl[I,0] = X[i]                                  |
//|                 Tbl[I,1] = X[i+1]                                |
//|                 Tbl[I,2] = C0                                    |
//|                 Tbl[I,3] = C1                                    |
//|                 Tbl[I,4] = C2                                    |
//|                 Tbl[I,5] = C3                                    |
//|             On [x[i], x[i+1]] spline is equals to:               |
//|                 S(x) = C0 + C1*t + C2*t^2 + C3*t^3               |
//|                 t = x-x[i]                                       |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DUnpack(CSpline1DInterpolant &c,int &n,
                                      CMatrixDouble &tbl)
  {
//--- create variables
   int i=0;
   int j=0;
//--- allocation
   tbl.Resize(c.m_n-2+1,2+c.m_k+1);
//--- initialization
   n=c.m_n;
//--- Fill
   for(i=0;i<=n-2;i++)
     {
      tbl[i].Set(0,c.m_x[i]);
      tbl[i].Set(1,c.m_x[i+1]);
      for(j=0;j<=c.m_k;j++)
         tbl[i].Set(2+j,c.m_c[(c.m_k+1)*i+j]);
     }
  }
//+------------------------------------------------------------------+
//| This subroutine performs linear transformation of the spline     |
//| argument.                                                        |
//| INPUT PARAMETERS:                                                |
//|     C   -   spline interpolant.                                  |
//|     A, B-   transformation coefficients: x = A*t + B             |
//| Result:                                                          |
//|     C   -   transformed spline                                   |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DLinTransX(CSpline1DInterpolant &c,const double a,
                                         const double b)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    n=0;
   double v=0;
   double dv=0;
   double d2v=0;
//--- create arrays
   double x[];
   double y[];
   double d[];
//--- initialization
   n=c.m_n;
//--- Special case: A=0
   if(a==0.0)
     {
      v=Spline1DCalc(c,b);
      for(i=0;i<=n-2;i++)
        {
         c.m_c[(c.m_k+1)*i]=v;
         for(j=1;j<=c.m_k;j++)
            c.m_c[(c.m_k+1)*i+j]=0;
        }
      //--- exit the function
      return;
     }
//--- General case: A<>0.
//--- Unpack,X,Y,dY/dX.
//--- Scale and pack again.
   if(!CAp::Assert(c.m_k==3,__FUNCTION__+": internal error"))
      return;
//--- allocation
   ArrayResizeAL(x,n);
   ArrayResizeAL(y,n);
   ArrayResizeAL(d,n);
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      x[i]=c.m_x[i];
      Spline1DDiff(c,x[i],v,dv,d2v);
      x[i]=(x[i]-b)/a;
      y[i]=v;
      d[i]=a*dv;
     }
//--- function call
   Spline1DBuildHermite(x,y,d,n,c);
  }
//+------------------------------------------------------------------+
//| This subroutine performs linear transformation of the spline.    |
//| INPUT PARAMETERS:                                                |
//|     C   -   spline interpolant.                                  |
//|     A,B-   transformation coefficients: S2(x)=A*S(x) + B         |
//| Result:                                                          |
//|     C   -   transformed spline                                   |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DLinTransY(CSpline1DInterpolant &c,const double a,
                                         const double b)
  {
//--- create variables
   int i=0;
   int j=0;
   int n=0;
//--- initialization
   n=c.m_n;
//--- calculation
   for(i=0;i<=n-2;i++)
     {
      c.m_c[(c.m_k+1)*i]=a*c.m_c[(c.m_k+1)*i]+b;
      for(j=1;j<=c.m_k;j++)
         c.m_c[(c.m_k+1)*i+j]=a*c.m_c[(c.m_k+1)*i+j];
     }
  }
//+------------------------------------------------------------------+
//| This subroutine integrates the spline.                           |
//| INPUT PARAMETERS:                                                |
//|     C   -   spline interpolant.                                  |
//|     X   -   right bound of the integration interval [a, x],      |
//|             here 'a' denotes min(x[])                            |
//| Result:                                                          |
//|     integral(S(t)dt,a,x)                                         |
//+------------------------------------------------------------------+
static double CSpline1D::Spline1DIntegrate(CSpline1DInterpolant &c,double x)
  {
//--- create variables
   double result=0;
   int    n=0;
   int    i=0;
   int    j=0;
   int    l=0;
   int    r=0;
   int    m=0;
   double w=0;
   double v=0;
   double t=0;
   double intab=0;
   double additionalterm=0;
//--- initialization
   n=c.m_n;
//--- Periodic splines require special treatment. We make
//--- following transformation:
//---     integral(S(t)dt,A,X)=integral(S(t)dt,A,Z)+AdditionalTerm
//--- here X may lie outside of [A,B],Z lies strictly in [A,B],
//--- AdditionalTerm is equals to integral(S(t)dt,A,B) times some
//--- integer number (may be zero).
   if(c.m_periodic && (x<c.m_x[0] || x>c.m_x[c.m_n-1]))
     {
      //--- compute integral(S(x)dx,A,B)
      intab=0;
      for(i=0;i<=c.m_n-2;i++)
        {
         w=c.m_x[i+1]-c.m_x[i];
         m=(c.m_k+1)*i;
         intab=intab+c.m_c[m]*w;
         v=w;
         for(j=1;j<=c.m_k;j++)
           {
            v=v*w;
            intab=intab+c.m_c[m+j]*v/(j+1);
           }
        }
      //--- map X into [A,B]
      CApServ::ApPeriodicMap(x,c.m_x[0],c.m_x[c.m_n-1],t);
      additionalterm=t*intab;
     }
   else
      additionalterm=0;
//--- Binary search in the [ x[0],...,x[n-2] ] (x[n-1] is not included)
   l=0;
   r=n-2+1;
   while(l!=r-1)
     {
      m=(l+r)/2;
      //--- check
      if(c.m_x[m]>=x)
         r=m;
      else
         l=m;
     }
//--- Integration
   result=0;
   for(i=0;i<=l-1;i++)
     {
      w=c.m_x[i+1]-c.m_x[i];
      m=(c.m_k+1)*i;
      result=result+c.m_c[m]*w;
      v=w;
      //--- calculation
      for(j=1;j<=c.m_k;j++)
        {
         v=v*w;
         result=result+c.m_c[m+j]*v/(j+1);
        }
     }
//--- change values
   w=x-c.m_x[l];
   m=(c.m_k+1)*l;
   v=w;
   result=result+c.m_c[m]*w;
//--- calculation
   for(j=1;j<=c.m_k;j++)
     {
      v=v*w;
      result=result+c.m_c[m+j]*v/(j+1);
     }
//--- return result
   return(result+additionalterm);
  }
//+------------------------------------------------------------------+
//| Internal version of Spline1DConvDiff                             |
//| Converts from Hermite spline given by grid XOld to new grid X2   |
//| INPUT PARAMETERS:                                                |
//|     XOld    -   old grid                                         |
//|     YOld    -   values at old grid                               |
//|     DOld    -   first derivative at old grid                     |
//|     N       -   grid size                                        |
//|     X2      -   new grid                                         |
//|     N2      -   new grid size                                    |
//|     Y       -   possibly preallocated output array               |
//|                 (reallocate if too small)                        |
//|     NeedY   -   do we need Y?                                    |
//|     D1      -   possibly preallocated output array               |
//|                 (reallocate if too small)                        |
//|     NeedD1  -   do we need D1?                                   |
//|     D2      -   possibly preallocated output array               |
//|                 (reallocate if too small)                        |
//|     NeedD2  -   do we need D1?                                   |
//| OUTPUT ARRAYS:                                                   |
//|     Y       -   values, if needed                                |
//|     D1      -   first derivative, if needed                      |
//|     D2      -   second derivative, if needed                     |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DConvDiffInternal(double &xold[],double &yold[],
                                                double &dold[],const int n,
                                                double &x2[],const int n2,
                                                double &y[],const bool needy,
                                                double &d1[],const bool needd1,
                                                double &d2[],const bool needd2)
  {
//--- create variables
   int    intervalindex=0;
   int    pointindex=0;
   bool   havetoadvance;
   double c0=0;
   double c1=0;
   double c2=0;
   double c3=0;
   double a=0;
   double b=0;
   double w=0;
   double w2=0;
   double w3=0;
   double fa=0;
   double fb=0;
   double da=0;
   double db=0;
   double t=0;
//--- Prepare space
   if(needy && CAp::Len(y)<n2)
      ArrayResizeAL(y,n2);
//--- check
   if(needd1 && CAp::Len(d1)<n2)
      ArrayResizeAL(d1,n2);
//--- check
   if(needd2 && CAp::Len(d2)<n2)
      ArrayResizeAL(d2,n2);
//--- These assignments aren't actually needed
//--- (variables are initialized in the loop below),
//--- but without them compiler will complain about uninitialized locals
   c0=0;
   c1=0;
   c2=0;
   c3=0;
   a=0;
   b=0;
//--- Cycle
   intervalindex=-1;
   pointindex=0;
//--- calculation
   while(true)
     {
      //--- are we ready to exit?
      if(pointindex>=n2)
         break;
      t=x2[pointindex];
      //--- do we need to advance interval?
      havetoadvance=false;
      //--- check
      if(intervalindex==-1)
         havetoadvance=true;
      else
        {
         //--- check
         if(intervalindex<n-2)
            havetoadvance=t>=b;
        }
      //--- check
      if(havetoadvance)
        {
         //--- change values
         intervalindex=intervalindex+1;
         a=xold[intervalindex];
         b=xold[intervalindex+1];
         w=b-a;
         w2=w*w;
         w3=w*w2;
         fa=yold[intervalindex];
         fb=yold[intervalindex+1];
         da=dold[intervalindex];
         db=dold[intervalindex+1];
         c0=fa;
         c1=da;
         c2=(3*(fb-fa)-2*da*w-db*w)/w2;
         c3=(2*(fa-fb)+da*w+db*w)/w3;
         continue;
        }
      //--- Calculate spline and its derivatives using power basis
      t=t-a;
      if(needy)
         y[pointindex]=c0+t*(c1+t*(c2+t*c3));
      //--- check
      if(needd1)
         d1[pointindex]=c1+2*t*c2+3*t*t*c3;
      //--- check
      if(needd2)
         d2[pointindex]=2*c2+6*t*c3;
      //--- change value
      pointindex=pointindex+1;
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine. Heap sort.                                  |
//+------------------------------------------------------------------+
static void CSpline1D::HeapSortDPoints(double &x[],double &y[],double &d[],
                                       const int n)
  {
//--- create variables
   int i=0;
   int i_=0;
//--- create arrays
   double rbuf[];
   int    ibuf[];
   double rbuf2[];
   int    ibuf2[];
//--- allocation
   ArrayResizeAL(ibuf,n);
   ArrayResizeAL(rbuf,n);
   for(i=0;i<=n-1;i++)
      ibuf[i]=i;
//--- function call
   CTSort::TagSortFastI(x,ibuf,rbuf2,ibuf2,n);
//--- copy
   for(i=0;i<=n-1;i++)
      rbuf[i]=y[ibuf[i]];
   for(i_=0;i_<=n-1;i_++)
      y[i_]=rbuf[i_];
   for(i=0;i<=n-1;i++)
      rbuf[i]=d[ibuf[i]];
   for(i_=0;i_<=n-1;i_++)
      d[i_]=rbuf[i_];
  }
//+------------------------------------------------------------------+
//| Internal version of Spline1DGridDiffCubic.                       |
//| Accepts pre-ordered X/Y, temporary arrays (which may be          |
//| preallocated, if you want to save time, or not) and output array |
//| (which may be preallocated too).                                 |
//| Y is passed as var-parameter because we may need to force last   |
//| element  to be equal to the first one (if periodic boundary      |
//| conditions are specified).                                       |
//+------------------------------------------------------------------+
static void CSpline1D::Spline1DGridDiffCubicInternal(double &x[],double &y[],
                                                     const int n,const int boundltype,
                                                     const double boundl,
                                                     const int boundrtype,
                                                     const double boundr,
                                                     double &d[],double &a1[],
                                                     double &a2[],double &a3[],
                                                     double &b[],double &dt[])
  {
//--- create variables
   int i=0;
   int i_=0;
//--- allocate arrays
   if(CAp::Len(d)<n)
      ArrayResizeAL(d,n);
//--- check
   if(CAp::Len(a1)<n)
      ArrayResizeAL(a1,n);
//--- check
   if(CAp::Len(a2)<n)
      ArrayResizeAL(a2,n);
//--- check
   if(CAp::Len(a3)<n)
      ArrayResizeAL(a3,n);
//--- check
   if(CAp::Len(b)<n)
      ArrayResizeAL(b,n);
//--- check
   if(CAp::Len(dt)<n)
      ArrayResizeAL(dt,n);
//--- Special cases:
//--- * N=2,parabolic terminated boundary condition on both ends
//--- * N=2,periodic boundary condition
   if((n==2 && boundltype==0) && boundrtype==0)
     {
      //--- change values
      d[0]=(y[1]-y[0])/(x[1]-x[0]);
      d[1]=d[0];
      //--- exit the function
      return;
     }
//--- check
   if((n==2 && boundltype==-1) && boundrtype==-1)
     {
      //--- change values
      d[0]=0;
      d[1]=0;
      //--- exit the function
      return;
     }
//--- Periodic and non-periodic boundary conditions are
//--- two separate classes
   if(boundrtype==-1 && boundltype==-1)
     {
      //--- Periodic boundary conditions
      y[n-1]=y[0];
      //--- Boundary conditions at N-1 points
      //--- (one point less because last point is the same as first point).
      a1[0]=x[1]-x[0];
      a2[0]=2*(x[1]-x[0]+x[n-1]-x[n-2]);
      a3[0]=x[n-1]-x[n-2];
      b[0]=3*(y[n-1]-y[n-2])/(x[n-1]-x[n-2])*(x[1]-x[0])+3*(y[1]-y[0])/(x[1]-x[0])*(x[n-1]-x[n-2]);
      //--- calculation
      for(i=1;i<=n-2;i++)
        {
         //--- Altough last point is [N-2],we use X[N-1] and Y[N-1]
         //--- (because of periodicity)
         a1[i]=x[i+1]-x[i];
         a2[i]=2*(x[i+1]-x[i-1]);
         a3[i]=x[i]-x[i-1];
         b[i]=3*(y[i]-y[i-1])/(x[i]-x[i-1])*(x[i+1]-x[i])+3*(y[i+1]-y[i])/(x[i+1]-x[i])*(x[i]-x[i-1]);
        }
      //--- Solve,add last point (with index N-1)
      SolveCyclicTridiagonal(a1,a2,a3,b,n-1,dt);
      for(i_=0;i_<=n-2;i_++)
         d[i_]=dt[i_];
      d[n-1]=d[0];
     }
   else
     {
      //--- Non-periodic boundary condition.
      //--- Left boundary conditions.
      if(boundltype==0)
        {
         //--- change values
         a1[0]=0;
         a2[0]=1;
         a3[0]=1;
         b[0]=2*(y[1]-y[0])/(x[1]-x[0]);
        }
      //--- check
      if(boundltype==1)
        {
         //--- change values
         a1[0]=0;
         a2[0]=1;
         a3[0]=0;
         b[0]=boundl;
        }
      //--- check
      if(boundltype==2)
        {
         //--- change values
         a1[0]=0;
         a2[0]=2;
         a3[0]=1;
         b[0]=3*(y[1]-y[0])/(x[1]-x[0])-0.5*boundl*(x[1]-x[0]);
        }
      //--- Central conditions
      for(i=1;i<=n-2;i++)
        {
         a1[i]=x[i+1]-x[i];
         a2[i]=2*(x[i+1]-x[i-1]);
         a3[i]=x[i]-x[i-1];
         b[i]=3*(y[i]-y[i-1])/(x[i]-x[i-1])*(x[i+1]-x[i])+3*(y[i+1]-y[i])/(x[i+1]-x[i])*(x[i]-x[i-1]);
        }
      //--- Right boundary conditions
      if(boundrtype==0)
        {
         //--- change values
         a1[n-1]=1;
         a2[n-1]=1;
         a3[n-1]=0;
         b[n-1]=2*(y[n-1]-y[n-2])/(x[n-1]-x[n-2]);
        }
      //--- check
      if(boundrtype==1)
        {
         //--- change values
         a1[n-1]=0;
         a2[n-1]=1;
         a3[n-1]=0;
         b[n-1]=boundr;
        }
      //--- check
      if(boundrtype==2)
        {
         //--- change values
         a1[n-1]=1;
         a2[n-1]=2;
         a3[n-1]=0;
         b[n-1]=3*(y[n-1]-y[n-2])/(x[n-1]-x[n-2])+0.5*boundr*(x[n-1]-x[n-2]);
        }
      //--- Solve
      SolveTridiagonal(a1,a2,a3,b,n,d);
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine. Heap sort.                                  |
//+------------------------------------------------------------------+
static void CSpline1D::HeapSortPoints(double &x[],double &y[],const int n)
  {
//--- create arrays
   double bufx[];
   double bufy[];
//--- function call
   CTSort::TagSortFastR(x,y,bufx,bufy,n);
  }
//+------------------------------------------------------------------+
//| Internal subroutine. Heap sort.                                  |
//| Accepts:                                                         |
//|     X, Y    -   points                                           |
//|     P       -   empty or preallocated array                      |
//| Returns:                                                         |
//|     X, Y    -   sorted by X                                      |
//|     P       -   array of permutations; I-th position of output   |
//| arrays X/Y contains(X[P[I]],Y[P[I]])                             |
//+------------------------------------------------------------------+
static void CSpline1D::HeapSortPPoints(double &x[],double &y[],int &p[],
                                       const int n)
  {
//--- create variables
   int i=0;
   int i_=0;
//--- create arrays
   double rbuf[];
   int    ibuf[];
//--- check
   if(CAp::Len(p)<n)
      ArrayResizeAL(p,n);
//--- allocation
   ArrayResizeAL(rbuf,n);
//--- initialization
   for(i=0;i<=n-1;i++)
      p[i]=i;
//--- function call
   CTSort::TagSortFastI(x,p,rbuf,ibuf,n);
//--- copy
   for(i=0;i<=n-1;i++)
      rbuf[i]=y[p[i]];
   for(i_=0;i_<=n-1;i_++)
      y[i_]=rbuf[i_];
  }
//+------------------------------------------------------------------+
//| Internal subroutine. Tridiagonal solver. Solves                  |
//| ( B[0] C[0]                      )                               |
//| ( A[1] B[1] C[1]                 )                               |
//| (      A[2] B[2] C[2]            )                               |
//| (            ..........          ) * X=D                         |
//| (            ..........          )                               |
//| (           A[N-2] B[N-2] C[N-2] )                               |
//| (                  A[N-1] B[N-1] )                               |
//+------------------------------------------------------------------+
static void CSpline1D::SolveTridiagonal(double &a[],double &cb[],double &c[],
                                        double &cd[],const int n,double &x[])
  {
//--- create variables
   int    k=0;
   double t=0;
//--- create arrays
   double d[];
   double b[];
//--- copy arrays
   ArrayCopy(d,cd);
   ArrayCopy(b,cb);
//--- check
   if(CAp::Len(x)<n)
      ArrayResizeAL(x,n);
//--- calculation
   for(k=1;k<=n-1;k++)
     {
      t=a[k]/b[k-1];
      b[k]=b[k]-t*c[k-1];
      d[k]=d[k]-t*d[k-1];
     }
   x[n-1]=d[n-1]/b[n-1];
   for(k=n-2;k>=0;k--)
      x[k]=(d[k]-c[k]*x[k+1])/b[k];
  }
//+------------------------------------------------------------------+
//| Internal subroutine. Cyclic tridiagonal solver. Solves           |
//| ( B[0] C[0]                 A[0] )                               |
//| ( A[1] B[1] C[1]                 )                               |
//| (      A[2] B[2] C[2]            )                               |
//| (            ..........          ) * X=D                         |
//| (            ..........          )                               |
//| (           A[N-2] B[N-2] C[N-2] )                               |
//| ( C[N-1]           A[N-1] B[N-1] )                               |
//+------------------------------------------------------------------+
static void CSpline1D::SolveCyclicTridiagonal(double &a[],double &cb[],
                                              double &c[],double &d[],
                                              const int n,double &x[])
  {
//--- create variables
   int    k=0;
   double alpha=0;
   double beta=0;
   double gamma=0;
//--- create arrays
   double y[];
   double z[];
   double u[];
   double b[];
//--- copy array
   ArrayCopy(b,cb);
//--- check
   if(CAp::Len(x)<n)
      ArrayResizeAL(x,n);
//--- change values
   beta=a[0];
   alpha=c[n-1];
   gamma=-b[0];
   b[0]=2*b[0];
   b[n-1]=b[n-1]-alpha*beta/gamma;
//--- allocation
   ArrayResizeAL(u,n);
//--- initialization
   for(k=0;k<=n-1;k++)
      u[k]=0;
   u[0]=gamma;
   u[n-1]=alpha;
//--- function call
   SolveTridiagonal(a,b,c,d,n,y);
//--- function call
   SolveTridiagonal(a,b,c,u,n,z);
//--- calculation
   for(k=0;k<=n-1;k++)
      x[k]=y[k]-(y[0]+beta/gamma*y[n-1])/(1+z[0]+beta/gamma*z[n-1])*z[k];
  }
//+------------------------------------------------------------------+
//| Internal subroutine. Three-point differentiation                 |
//+------------------------------------------------------------------+
static double CSpline1D::DiffThreePoint(double t,const double x0,const double f0,
                                        double x1,const double f1,double x2,
                                        const double f2)
  {
//--- create variables
   double a=0;
   double b=0;
//--- change values
   t=t-x0;
   x1=x1-x0;
   x2=x2-x0;
   a=(f2-f0-x2/x1*(f1-f0))/(CMath::Sqr(x2)-x1*x2);
   b=(f1-f0-a*CMath::Sqr(x1))/x1;
//--- return result
   return(2*a*t+b);
  }
//+------------------------------------------------------------------+
//| Polynomial fitting report:                                       |
//|     TaskRCond       reciprocal of task's condition number        |
//|     RMSError        RMS error                                    |
//|     AvgError        average error                                |
//|     AvgRelError     average relative error (for non-zero Y[I])   |
//|     MaxError        maximum error                                |
//+------------------------------------------------------------------+
class CPolynomialFitReport
  {
public:
   //--- variables
   double            m_taskrcond;
   double            m_rmserror;
   double            m_avgerror;
   double            m_avgrelerror;
   double            m_maxerror;
   //--- constructor, destructor
                     CPolynomialFitReport(void);
                    ~CPolynomialFitReport(void);
   //--- copy
   void              Copy(CPolynomialFitReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CPolynomialFitReport::CPolynomialFitReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CPolynomialFitReport::~CPolynomialFitReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CPolynomialFitReport::Copy(CPolynomialFitReport &obj)
  {
//--- copy variables
   m_taskrcond=obj.m_taskrcond;
   m_rmserror=obj.m_rmserror;
   m_avgerror=obj.m_avgerror;
   m_avgrelerror=obj.m_avgrelerror;
   m_maxerror=obj.m_maxerror;
  }
//+------------------------------------------------------------------+
//| Polynomial fitting report:                                       |
//|     TaskRCond       reciprocal of task's condition number        |
//|     RMSError        RMS error                                    |
//|     AvgError        average error                                |
//|     AvgRelError     average relative error (for non-zero Y[I])   |
//|     MaxError        maximum error                                |
//+------------------------------------------------------------------+
class CPolynomialFitReportShell
  {
private:
   CPolynomialFitReport m_innerobj;
public:
   //--- constructors, destructor
                     CPolynomialFitReportShell(void);
                     CPolynomialFitReportShell(CPolynomialFitReport &obj);
                    ~CPolynomialFitReportShell(void);
   //--- methods
   double            GetTaskRCond(void);
   void              SetTaskRCond(const double d);
   double            GetRMSError(void);
   void              SetRMSError(const double d);
   double            GetAvgError(void);
   void              SetAvgError(const double d);
   double            GetAvgRelError(void);
   void              SetAvgRelError(const double d);
   double            GetMaxError(void);
   void              SetMaxError(const double d);
   CPolynomialFitReport *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CPolynomialFitReportShell::CPolynomialFitReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CPolynomialFitReportShell::CPolynomialFitReportShell(CPolynomialFitReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CPolynomialFitReportShell::~CPolynomialFitReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable taskrcond                      |
//+------------------------------------------------------------------+
double CPolynomialFitReportShell::GetTaskRCond(void)
  {
//--- return result
   return(m_innerobj.m_taskrcond);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable taskrcond                     |
//+------------------------------------------------------------------+
void CPolynomialFitReportShell::SetTaskRCond(const double d)
  {
//--- change value
   m_innerobj.m_taskrcond=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable rmserror                       |
//+------------------------------------------------------------------+
double CPolynomialFitReportShell::GetRMSError(void)
  {
//--- return result
   return(m_innerobj.m_rmserror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable rmserror                      |
//+------------------------------------------------------------------+
void CPolynomialFitReportShell::SetRMSError(const double d)
  {
//--- change value
   m_innerobj.m_rmserror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgerror                       |
//+------------------------------------------------------------------+
double CPolynomialFitReportShell::GetAvgError(void)
  {
//--- return result
   return(m_innerobj.m_avgerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgerror                      |
//+------------------------------------------------------------------+
void CPolynomialFitReportShell::SetAvgError(const double d)
  {
//--- change value
   m_innerobj.m_avgerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgrelerror                    |
//+------------------------------------------------------------------+
double CPolynomialFitReportShell::GetAvgRelError(void)
  {
//--- return result
   return(m_innerobj.m_avgrelerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgrelerror                   |
//+------------------------------------------------------------------+
void CPolynomialFitReportShell::SetAvgRelError(const double d)
  {
//--- change value
   m_innerobj.m_avgrelerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable maxerror                       |
//+------------------------------------------------------------------+
double CPolynomialFitReportShell::GetMaxError(void)
  {
//--- return result
   return(m_innerobj.m_maxerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable maxerror                      |
//+------------------------------------------------------------------+
void CPolynomialFitReportShell::SetMaxError(const double d)
  {
//--- change value
   m_innerobj.m_maxerror=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CPolynomialFitReport *CPolynomialFitReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Barycentric fitting report:                                      |
//|     RMSError        RMS error                                    |
//|     AvgError        average error                                |
//|     AvgRelError     average relative error (for non-zero Y[I])   |
//|     MaxError        maximum error                                |
//|     TaskRCond       reciprocal of task's condition number        |
//+------------------------------------------------------------------+
class CBarycentricFitReport
  {
public:
   //--- variables
   double            m_taskrcond;
   int               m_dbest;
   double            m_rmserror;
   double            m_avgerror;
   double            m_avgrelerror;
   double            m_maxerror;
   //--- constructor, destructor
                     CBarycentricFitReport(void);
                    ~CBarycentricFitReport(void);
   //--- copy
   void              Copy(CBarycentricFitReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CBarycentricFitReport::CBarycentricFitReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CBarycentricFitReport::~CBarycentricFitReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CBarycentricFitReport::Copy(CBarycentricFitReport &obj)
  {
//--- copy variables
   m_taskrcond=obj.m_taskrcond;
   m_dbest=obj.m_dbest;
   m_rmserror=obj.m_rmserror;
   m_avgerror=obj.m_avgerror;
   m_avgrelerror=obj.m_avgrelerror;
   m_maxerror=obj.m_maxerror;
  }
//+------------------------------------------------------------------+
//| Barycentric fitting report:                                      |
//|     RMSError        RMS error                                    |
//|     AvgError        average error                                |
//|     AvgRelError     average relative error (for non-zero Y[I])   |
//|     MaxError        maximum error                                |
//|     TaskRCond       reciprocal of task's condition number        |
//+------------------------------------------------------------------+
class CBarycentricFitReportShell
  {
private:
   CBarycentricFitReport m_innerobj;
public:
   //--- constructors, destructor
                     CBarycentricFitReportShell(void);
                     CBarycentricFitReportShell(CBarycentricFitReport &obj);
                    ~CBarycentricFitReportShell(void);
   //--- methods
   double            GetTaskRCond(void);
   void              SetTaskRCond(const double d);
   int               GetDBest(void);
   void              SetDBest(const int i);
   double            GetRMSError(void);
   void              SetRMSError(const double d);
   double            GetAvgError(void);
   void              SetAvgError(const double d);
   double            GetAvgRelError(void);
   void              SetAvgRelError(const double d);
   double            GetMaxError(void);
   void              SetMaxError(const double d);
   CBarycentricFitReport *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CBarycentricFitReportShell::CBarycentricFitReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CBarycentricFitReportShell::CBarycentricFitReportShell(CBarycentricFitReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CBarycentricFitReportShell::~CBarycentricFitReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable taskrcond                      |
//+------------------------------------------------------------------+
double CBarycentricFitReportShell::GetTaskRCond(void)
  {
//--- return result
   return(m_innerobj.m_taskrcond);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable taskrcond                     |
//+------------------------------------------------------------------+
void CBarycentricFitReportShell::SetTaskRCond(const double d)
  {
//--- change value
   m_innerobj.m_taskrcond=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable dbest                          |
//+------------------------------------------------------------------+
int CBarycentricFitReportShell::GetDBest(void)
  {
//--- return result
   return(m_innerobj.m_dbest);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable dbest                         |
//+------------------------------------------------------------------+
void CBarycentricFitReportShell::SetDBest(const int i)
  {
//--- change value
   m_innerobj.m_dbest=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable rmserror                       |
//+------------------------------------------------------------------+
double CBarycentricFitReportShell::GetRMSError(void)
  {
//--- return result
   return(m_innerobj.m_rmserror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable rmserror                      |
//+------------------------------------------------------------------+
void CBarycentricFitReportShell::SetRMSError(const double d)
  {
//--- change value
   m_innerobj.m_rmserror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgerror                       |
//+------------------------------------------------------------------+
double CBarycentricFitReportShell::GetAvgError(void)
  {
//--- return result
   return(m_innerobj.m_avgerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgerror                      |
//+------------------------------------------------------------------+
void CBarycentricFitReportShell::SetAvgError(const double d)
  {
//--- change value
   m_innerobj.m_avgerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgrelerror                    |
//+------------------------------------------------------------------+
double CBarycentricFitReportShell::GetAvgRelError(void)
  {
//--- return result
   return(m_innerobj.m_avgrelerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgrelerror                   |
//+------------------------------------------------------------------+
void CBarycentricFitReportShell::SetAvgRelError(const double d)
  {
//--- change value
   m_innerobj.m_avgrelerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable maxerror                       |
//+------------------------------------------------------------------+
double CBarycentricFitReportShell::GetMaxError(void)
  {
//--- return result
   return(m_innerobj.m_maxerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable maxerror                      |
//+------------------------------------------------------------------+
void CBarycentricFitReportShell::SetMaxError(const double d)
  {
//--- change value
   m_innerobj.m_maxerror=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CBarycentricFitReport *CBarycentricFitReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Spline fitting report:                                           |
//|     RMSError        RMS error                                    |
//|     AvgError        average error                                |
//|     AvgRelError     average relative error (for non-zero Y[I])   |
//|     MaxError        maximum error                                |
//| Fields below are filled by obsolete functions (Spline1DFitCubic, |
//| Spline1DFitHermite). Modern fitting functions do NOT fill these  | 
//| fields:                                                          |
//|     TaskRCond       reciprocal of task's condition number        |
//+------------------------------------------------------------------+
class CSpline1DFitReport
  {
public:
   //--- variables
   double            m_taskrcond;
   double            m_rmserror;
   double            m_avgerror;
   double            m_avgrelerror;
   double            m_maxerror;
   //--- constructor, destructor
                     CSpline1DFitReport(void);
                    ~CSpline1DFitReport(void);
   //--- copy
   void              Copy(CSpline1DFitReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSpline1DFitReport::CSpline1DFitReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSpline1DFitReport::~CSpline1DFitReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CSpline1DFitReport::Copy(CSpline1DFitReport &obj)
  {
//--- copy variables
   m_taskrcond=obj.m_taskrcond;
   m_rmserror=obj.m_rmserror;
   m_avgerror=obj.m_avgerror;
   m_avgrelerror=obj.m_avgrelerror;
   m_maxerror=obj.m_maxerror;
  }
//+------------------------------------------------------------------+
//| Spline fitting report:                                           |
//|     RMSError        RMS error                                    |
//|     AvgError        average error                                |
//|     AvgRelError     average relative error (for non-zero Y[I])   |
//|     MaxError        maximum error                                |
//| Fields below are filled by obsolete functions (Spline1DFitCubic, |
//| Spline1DFitHermite). Modern fitting functions do NOT fill these  | 
//| fields:                                                          |
//|     TaskRCond       reciprocal of task's condition number        |
//+------------------------------------------------------------------+
class CSpline1DFitReportShell
  {
private:
   CSpline1DFitReport m_innerobj;
public:
   //--- constructors, destructor
                     CSpline1DFitReportShell(void);
                     CSpline1DFitReportShell(CSpline1DFitReport &obj);
                    ~CSpline1DFitReportShell(void);
   //--- methods
   double            GetTaskRCond(void);
   void              SetTaskRCond(const double d);
   double            GetRMSError(void);
   void              SetRMSError(const double d);
   double            GetAvgError(void);
   void              SetAvgError(const double d);
   double            GetAvgRelError(void);
   void              SetAvgRelError(const double d);
   double            GetMaxError(void);
   void              SetMaxError(const double d);
   CSpline1DFitReport *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSpline1DFitReportShell::CSpline1DFitReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CSpline1DFitReportShell::CSpline1DFitReportShell(CSpline1DFitReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSpline1DFitReportShell::~CSpline1DFitReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable taskrcond                      |
//+------------------------------------------------------------------+
double CSpline1DFitReportShell::GetTaskRCond(void)
  {
//--- return result
   return(m_innerobj.m_taskrcond);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable taskrcond                     |
//+------------------------------------------------------------------+
void CSpline1DFitReportShell::SetTaskRCond(const double d)
  {
//--- change value
   m_innerobj.m_taskrcond=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable rmserror                       |
//+------------------------------------------------------------------+
double CSpline1DFitReportShell::GetRMSError(void)
  {
//--- return result
   return(m_innerobj.m_rmserror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable rmserror                      |
//+------------------------------------------------------------------+
void CSpline1DFitReportShell::SetRMSError(const double d)
  {
//--- change value
   m_innerobj.m_rmserror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgerror                       |
//+------------------------------------------------------------------+
double CSpline1DFitReportShell::GetAvgError(void)
  {
//--- return result
   return(m_innerobj.m_avgerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgerror                      |
//+------------------------------------------------------------------+
void CSpline1DFitReportShell::SetAvgError(const double d)
  {
//--- change value
   m_innerobj.m_avgerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgrelerror                    |
//+------------------------------------------------------------------+
double CSpline1DFitReportShell::GetAvgRelError(void)
  {
//--- return result
   return(m_innerobj.m_avgrelerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgrelerror                   |
//+------------------------------------------------------------------+
void CSpline1DFitReportShell::SetAvgRelError(const double d)
  {
//--- change value
   m_innerobj.m_avgrelerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable maxerror                       |
//+------------------------------------------------------------------+
double CSpline1DFitReportShell::GetMaxError(void)
  {
//--- return result
   return(m_innerobj.m_maxerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable maxerror                      |
//+------------------------------------------------------------------+
void CSpline1DFitReportShell::SetMaxError(const double d)
  {
//--- change value
   m_innerobj.m_maxerror=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CSpline1DFitReport *CSpline1DFitReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Least squares fitting report:                                    |
//|     TaskRCond       reciprocal of task's condition number        |
//|     IterationsCount number of internal iterations                |
//|     RMSError        RMS error                                    |
//|     AvgError        average error                                |
//|     AvgRelError     average relative error (for non-zero Y[I])   |
//|     MaxError        maximum error                                |
//|     WRMSError       weighted RMS error                           |
//+------------------------------------------------------------------+
class CLSFitReport
  {
public:
   //--- variables
   double            m_taskrcond;
   int               m_iterationscount;
   double            m_rmserror;
   double            m_avgerror;
   double            m_avgrelerror;
   double            m_maxerror;
   double            m_wrmserror;
   //--- constructor, destructor
                     CLSFitReport(void);
                    ~CLSFitReport(void);
   //--- copy
   void              Copy(CLSFitReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLSFitReport::CLSFitReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLSFitReport::~CLSFitReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CLSFitReport::Copy(CLSFitReport &obj)
  {
//--- copy variables
   m_taskrcond=obj.m_taskrcond;
   m_iterationscount=obj.m_iterationscount;
   m_rmserror=obj.m_rmserror;
   m_avgerror=obj.m_avgerror;
   m_avgrelerror=obj.m_avgrelerror;
   m_maxerror=obj.m_maxerror;
   m_wrmserror=obj.m_wrmserror;
  }
//+------------------------------------------------------------------+
//| Least squares fitting report:                                    |
//|     TaskRCond       reciprocal of task's condition number        |
//|     IterationsCount number of internal iterations                |
//|     RMSError        RMS error                                    |
//|     AvgError        average error                                |
//|     AvgRelError     average relative error (for non-zero Y[I])   |
//|     MaxError        maximum error                                |
//|     WRMSError       weighted RMS error                           |
//+------------------------------------------------------------------+
class CLSFitReportShell
  {
private:
   CLSFitReport      m_innerobj;
public:
   //--- constructors, destructor
                     CLSFitReportShell(void);
                     CLSFitReportShell(CLSFitReport &obj);
                    ~CLSFitReportShell(void);
   //--- methods
   double            GetTaskRCond(void);
   void              SetTaskRCond(const double d);
   int               GetIterationsCount(void);
   void              SetIterationsCount(const int i);
   double            GetRMSError(void);
   void              SetRMSError(const double d);
   double            GetAvgError(void);
   void              SetAvgError(const double d);
   double            GetAvgRelError(void);
   void              SetAvgRelError(const double d);
   double            GetMaxError(void);
   void              SetMaxError(const double d);
   double            GetWRMSError(void);
   void              SetWRMSError(const double d);
   CLSFitReport     *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLSFitReportShell::CLSFitReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CLSFitReportShell::CLSFitReportShell(CLSFitReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLSFitReportShell::~CLSFitReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable taskrcond                      |
//+------------------------------------------------------------------+
double CLSFitReportShell::GetTaskRCond(void)
  {
//--- return result
   return(m_innerobj.m_taskrcond);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable taskrcond                     |
//+------------------------------------------------------------------+
void CLSFitReportShell::SetTaskRCond(const double d)
  {
//--- change value
   m_innerobj.m_taskrcond=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable iterationscount                |
//+------------------------------------------------------------------+
int CLSFitReportShell::GetIterationsCount(void)
  {
//--- return result
   return(m_innerobj.m_iterationscount);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable iterationscount               |
//+------------------------------------------------------------------+
void CLSFitReportShell::SetIterationsCount(const int i)
  {
//--- change value
   m_innerobj.m_iterationscount=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable rmserror                       |
//+------------------------------------------------------------------+
double CLSFitReportShell::GetRMSError(void)
  {
//--- return result
   return(m_innerobj.m_rmserror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable rmserror                      |
//+------------------------------------------------------------------+
void CLSFitReportShell::SetRMSError(const double d)
  {
//--- change value
   m_innerobj.m_rmserror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgerror                       |
//+------------------------------------------------------------------+
double CLSFitReportShell::GetAvgError(void)
  {
//--- return result
   return(m_innerobj.m_avgerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgerror                      |
//+------------------------------------------------------------------+
void CLSFitReportShell::SetAvgError(const double d)
  {
//--- change value
   m_innerobj.m_avgerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgrelerror                    |
//+------------------------------------------------------------------+
double CLSFitReportShell::GetAvgRelError(void)
  {
//--- return result
   return(m_innerobj.m_avgrelerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgrelerror                   |
//+------------------------------------------------------------------+
void CLSFitReportShell::SetAvgRelError(const double d)
  {
//--- change value
   m_innerobj.m_avgrelerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable maxerror                       |
//+------------------------------------------------------------------+
double CLSFitReportShell::GetMaxError(void)
  {
//--- return result
   return(m_innerobj.m_maxerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable maxerror                      |
//+------------------------------------------------------------------+
void CLSFitReportShell::SetMaxError(const double d)
  {
//--- change value
   m_innerobj.m_maxerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable wrmserror                      |
//+------------------------------------------------------------------+
double CLSFitReportShell::GetWRMSError(void)
  {
//--- return result
   return(m_innerobj.m_wrmserror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable wrmserror                     |
//+------------------------------------------------------------------+
void CLSFitReportShell::SetWRMSError(const double d)
  {
//--- change value
   m_innerobj.m_wrmserror=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CLSFitReport *CLSFitReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Nonlinear fitter.                                                |
//| You should use ALGLIB functions to work with fitter.             |
//| Never try to access its fields directly!                         |
//+------------------------------------------------------------------+
class CLSFitState
  {
public:
   //--- variables
   int               m_optalgo;
   int               m_m;
   int               m_k;
   double            m_f;
   double            m_epsf;
   double            m_epsx;
   int               m_maxits;
   double            m_stpmax;
   bool              m_xrep;
   int               m_npoints;
   int               m_nweights;
   int               m_wkind;
   int               m_wits;
   bool              m_xupdated;
   bool              m_needf;
   bool              m_needfg;
   bool              m_needfgh;
   int               m_pointindex;
   int               m_repiterationscount;
   int               m_repterminationtype;
   double            m_reprmserror;
   double            m_repavgerror;
   double            m_repavgrelerror;
   double            m_repmaxerror;
   double            m_repwrmserror;
   CMinLMState       m_optstate;
   CMinLMReport      m_optrep;
   int               m_prevnpt;
   int               m_prevalgo;
   RCommState        m_rstate;
   //--- arrays
   double            m_s[];
   double            m_bndl[];
   double            m_bndu[];
   double            m_tasky[];
   double            m_w[];
   double            m_x[];
   double            m_c[];
   double            m_g[];
   //--- matrices
   CMatrixDouble     m_taskx;
   CMatrixDouble     m_h;
   //--- constructor, destructor
                     CLSFitState(void);
                    ~CLSFitState(void);
   //--- copy
   void              Copy(CLSFitState &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLSFitState::CLSFitState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLSFitState::~CLSFitState(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CLSFitState::Copy(CLSFitState &obj)
  {
//--- copy variables
   m_optalgo=obj.m_optalgo;
   m_m=obj.m_m;
   m_k=obj.m_k;
   m_f=obj.m_f;
   m_epsf=obj.m_epsf;
   m_epsx=obj.m_epsx;
   m_maxits=obj.m_maxits;
   m_stpmax=obj.m_stpmax;
   m_xrep=obj.m_xrep;
   m_npoints=obj.m_npoints;
   m_nweights=obj.m_nweights;
   m_wkind=obj.m_wkind;
   m_wits=obj.m_wits;
   m_xupdated=obj.m_xupdated;
   m_needf=obj.m_needf;
   m_needfg=obj.m_needfg;
   m_needfgh=obj.m_needfgh;
   m_pointindex=obj.m_pointindex;
   m_repiterationscount=obj.m_repiterationscount;
   m_repterminationtype=obj.m_repterminationtype;
   m_reprmserror=obj.m_reprmserror;
   m_repavgerror=obj.m_repavgerror;
   m_repavgrelerror=obj.m_repavgrelerror;
   m_repmaxerror=obj.m_repmaxerror;
   m_repwrmserror=obj.m_repwrmserror;
   m_prevnpt=obj.m_prevnpt;
   m_prevalgo=obj.m_prevalgo;
   m_optstate.Copy(obj.m_optstate);
   m_optrep.Copy(obj.m_optrep);
   m_rstate.Copy(obj.m_rstate);
//--- copy arrays
   ArrayCopy(m_s,obj.m_s);
   ArrayCopy(m_bndl,obj.m_bndl);
   ArrayCopy(m_bndu,obj.m_bndu);
   ArrayCopy(m_tasky,obj.m_tasky);
   ArrayCopy(m_w,obj.m_w);
   ArrayCopy(m_x,obj.m_x);
   ArrayCopy(m_c,obj.m_c);
   ArrayCopy(m_g,obj.m_g);
//--- copy matrices
   m_taskx=obj.m_taskx;
   m_h=obj.m_h;
  }
//+------------------------------------------------------------------+
//| Nonlinear fitter.                                                |
//| You should use ALGLIB functions to work with fitter.             |
//| Never try to access its fields directly!                         |
//+------------------------------------------------------------------+
class CLSFitStateShell
  {
private:
   CLSFitState       m_innerobj;
public:
   //--- constructors, destructor
                     CLSFitStateShell(void);
                     CLSFitStateShell(CLSFitState &obj);
                    ~CLSFitStateShell(void);
   //--- methods
   bool              GetNeedF(void);
   void              SetNeedF(const bool b);
   bool              GetNeedFG(void);
   void              SetNeedFG(const bool b);
   bool              GetNeedFGH(void);
   void              SetNeedFGH(const bool b);
   bool              GetXUpdated(void);
   void              SetXUpdated(const bool b);
   double            GetF(void);
   void              SetF(const double d);
   CLSFitState      *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLSFitStateShell::CLSFitStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CLSFitStateShell::CLSFitStateShell(CLSFitState &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLSFitStateShell::~CLSFitStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needf                          |
//+------------------------------------------------------------------+
bool CLSFitStateShell::GetNeedF(void)
  {
//--- return result
   return(m_innerobj.m_needf);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needf                         |
//+------------------------------------------------------------------+
void CLSFitStateShell::SetNeedF(const bool b)
  {
//--- change value
   m_innerobj.m_needf=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needfg                         |
//+------------------------------------------------------------------+
bool CLSFitStateShell::GetNeedFG(void)
  {
//--- return result
   return(m_innerobj.m_needfg);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needfg                        |
//+------------------------------------------------------------------+
void CLSFitStateShell::SetNeedFG(const bool b)
  {
//--- change value
   m_innerobj.m_needfg=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable needfgh                        |
//+------------------------------------------------------------------+
bool CLSFitStateShell::GetNeedFGH(void)
  {
//--- return result
   return(m_innerobj.m_needfgh);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable needfgh                       |
//+------------------------------------------------------------------+
void CLSFitStateShell::SetNeedFGH(const bool b)
  {
//--- change value
   m_innerobj.m_needfgh=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable xupdated                       |
//+------------------------------------------------------------------+
bool CLSFitStateShell::GetXUpdated(void)
  {
//--- return result
   return(m_innerobj.m_xupdated);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable xupdated                      |
//+------------------------------------------------------------------+
void CLSFitStateShell::SetXUpdated(const bool b)
  {
//--- change value
   m_innerobj.m_xupdated=b;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable f                              |
//+------------------------------------------------------------------+
double CLSFitStateShell::GetF(void)
  {
//--- return result
   return(m_innerobj.m_f);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable f                             |
//+------------------------------------------------------------------+
void CLSFitStateShell::SetF(const double d)
  {
//--- change value
   m_innerobj.m_f=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CLSFitState *CLSFitStateShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Least squares fitting                                            |
//+------------------------------------------------------------------+
class CLSFit
  {
private:
   //--- private methods
   static void       Spline1DFitInternal(const int st,double &cx[],double &cy[],double &cw[],const int n,double &cxc[],double &cyc[],int &dc[],const int k,const int m,int &info,CSpline1DInterpolant &s,CSpline1DFitReport &rep);
   static void       LSFitLinearInternal(double &y[],double &w[],CMatrixDouble &fmatrix,const int n,const int m,int &info,double &c[],CLSFitReport &rep);
   static void       LSFitClearRequestFields(CLSFitState &state);
   static void       BarycentricCalcBasis(CBarycentricInterpolant &b,const double t,double &y[]);
   static void       InternalChebyshevFit(double &x[],double &y[],double &w[],const int n,double &cxc[],double &cyc[],int &dc[],const int k,const int m,int &info,double &c[],CLSFitReport &rep);
   static void       BarycentricFitWCFixedD(double &cx[],double &cy[],double &cw[],const int n,double &cxc[],double &cyc[],int &dc[],const int k,const int m,const int d,int &info,CBarycentricInterpolant &b,CBarycentricFitReport &rep);
   //--- auxiliary functions for LSFitIteration
   static void       Func_lbl_rcomm(CLSFitState &state,int n,int m,int k,int i,int j,double v,double vv,double relcnt);
   static bool       Func_lbl_7(CLSFitState &state,int &n,int &m,int &k,int &i,int &j,double &v,double &vv,double &relcnt);
   static bool       Func_lbl_11(CLSFitState &state,int &n,int &m,int &k,int &i,int &j,double &v,double &vv,double &relcnt);
   static bool       Func_lbl_14(CLSFitState &state,int &n,int &m,int &k,int &i,int &j,double &v,double &vv,double &relcnt);
   static bool       Func_lbl_16(CLSFitState &state,int &n,int &m,int &k,int &i,int &j,double &v,double &vv,double &relcnt);
   static bool       Func_lbl_21(CLSFitState &state,int &n,int &m,int &k,int &i,int &j,double &v,double &vv,double &relcnt);
   static bool       Func_lbl_24(CLSFitState &state,int &n,int &m,int &k,int &i,int &j,double &v,double &vv,double &relcnt);
   static bool       Func_lbl_26(CLSFitState &state,int &n,int &m,int &k,int &i,int &j,double &v,double &vv,double &relcnt);
   static bool       Func_lbl_29(CLSFitState &state,int &n,int &m,int &k,int &i,int &j,double &v,double &vv,double &relcnt);
   static bool       Func_lbl_31(CLSFitState &state,int &n,int &m,int &k,int &i,int &j,double &v,double &vv,double &relcnt);
   static bool       Func_lbl_38(CLSFitState &state,int &n,int &m,int &k,int &i,int &j,double &v,double &vv,double &relcnt);
public:
   //--- class constant
   static const int  m_rfsmax;
   //--- constructor, destructor
                     CLSFit(void);
                    ~CLSFit(void);
   //--- public methods
   static void       PolynomialFit(double &x[],double &y[],const int n,const int m,int &info,CBarycentricInterpolant &p,CPolynomialFitReport &rep);
   static void       PolynomialFitWC(double &cx[],double &cy[],double &cw[],const int n,double &cxc[],double &cyc[],int &dc[],const int k,const int m,int &info,CBarycentricInterpolant &p,CPolynomialFitReport &rep);
   static void       BarycentricFitFloaterHormannWC(double &x[],double &y[],double &w[],const int n,double &xc[],double &yc[],int &dc[],const int k,const int m,int &info,CBarycentricInterpolant &b,CBarycentricFitReport &rep);
   static void       BarycentricFitFloaterHormann(double &x[],double &y[],const int n,const int m,int &info,CBarycentricInterpolant &b,CBarycentricFitReport &rep);
   static void       Spline1DFitPenalized(double &cx[],double &cy[],const int n,const int m,const double rho,int &info,CSpline1DInterpolant &s,CSpline1DFitReport &rep);
   static void       Spline1DFitPenalizedW(double &cx[],double &cy[],double &cw[],const int n,const int m,double rho,int &info,CSpline1DInterpolant &s,CSpline1DFitReport &rep);
   static void       Spline1DFitCubicWC(double &x[],double &y[],double &w[],const int n,double &xc[],double &yc[],int &dc[],const int k,const int m,int &info,CSpline1DInterpolant &s,CSpline1DFitReport &rep);
   static void       Spline1DFitHermiteWC(double &x[],double &y[],double &w[],const int n,double &xc[],double &yc[],int &dc[],const int k,const int m,int &info,CSpline1DInterpolant &s,CSpline1DFitReport &rep);
   static void       Spline1DFitCubic(double &x[],double &y[],const int n,const int m,int &info,CSpline1DInterpolant &s,CSpline1DFitReport &rep);
   static void       Spline1DFitHermite(double &x[],double &y[],const int n,const int m,int &info,CSpline1DInterpolant &s,CSpline1DFitReport &rep);
   static void       LSFitLinearW(double &y[],double &w[],CMatrixDouble &fmatrix,const int n,const int m,int &info,double &c[],CLSFitReport &rep);
   static void       LSFitLinearWC(double &cy[],double &w[],CMatrixDouble &fmatrix,CMatrixDouble &ccmatrix,const int n,const int m,const int k,int &info,double &c[],CLSFitReport &rep);
   static void       LSFitLinear(double &y[],CMatrixDouble &fmatrix,const int n,const int m,int &info,double &c[],CLSFitReport &rep);
   static void       LSFitLinearC(double &cy[],CMatrixDouble &fmatrix,CMatrixDouble &cmatrix,const int n,const int m,const int k,int &info,double &c[],CLSFitReport &rep);
   static void       LSFitCreateWF(CMatrixDouble &x,double &y[],double &w[],double &c[],const int n,const int m,const int k,const double diffstep,CLSFitState &state);
   static void       LSFitCreateF(CMatrixDouble &x,double &y[],double &c[],const int n,const int m,const int k,const double diffstep,CLSFitState &state);
   static void       LSFitCreateWFG(CMatrixDouble &x,double &y[],double &w[],double &c[],const int n,const int m,const int k,bool cheapfg,CLSFitState &state);
   static void       LSFitCreateFG(CMatrixDouble &x,double &y[],double &c[],const int n,const int m,const int k,const bool cheapfg,CLSFitState &state);
   static void       LSFitCreateWFGH(CMatrixDouble &x,double &y[],double &w[],double &c[],const int n,const int m,const int k,CLSFitState &state);
   static void       LSFitCreateFGH(CMatrixDouble &x,double &y[],double &c[],const int n,const int m,const int k,CLSFitState &state);
   static void       LSFitSetCond(CLSFitState &state,const double epsf,const double epsx,const int maxits);
   static void       LSFitSetStpMax(CLSFitState &state,const double stpmax);
   static void       LSFitSetXRep(CLSFitState &state,const bool needxrep);
   static void       LSFitSetScale(CLSFitState &state,double &s[]);
   static void       LSFitSetBC(CLSFitState &state,double &bndl[],double &bndu[]);
   static void       LSFitResults(CLSFitState &state,int &info,double &c[],CLSFitReport &rep);
   static void       LSFitScaleXY(double &x[],double &y[],double &w[],const int n,double &xc[],double &yc[],int &dc[],const int k,double &xa,double &xb,double &sa,double &sb,double &xoriginal[],double &yoriginal[]);
   static bool       LSFitIteration(CLSFitState &state);
  };
//+------------------------------------------------------------------+
//| Initialize constant                                              |
//+------------------------------------------------------------------+
const int CLSFit::m_rfsmax=10;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLSFit::CLSFit(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLSFit::~CLSFit(void)
  {

  }
//+------------------------------------------------------------------+
//| Fitting by polynomials in barycentric form. This function        |
//| provides simple unterface for unconstrained unweighted fitting.  |
//| See PolynomialFitWC() if you need constrained fitting.           |
//| Task is linear, so linear least squares solver is used.          |
//| Complexity of this computational scheme is O(N*M^2), mostly      |
//| dominated by least squares solver                                |
//| SEE ALSO:                                                        |
//|     PolynomialFitWC()                                            |
//| INPUT PARAMETERS:                                                |
//|     X   -   points, array[0..N-1].                               |
//|     Y   -   function values, array[0..N-1].                      |
//|     N   -   number of points, N>0                                |
//|             * if given, only leading N elements of X/Y are used  |
//|             * if not given, automatically determined from sizes  |
//|               of X/Y                                             |
//|     M   -   number of basis functions (= polynomial_degree + 1), |
//|             M>=1                                                 |
//| OUTPUT PARAMETERS:                                               |
//|     Info-   same format as in LSFitLinearW() subroutine:         |
//|             * Info>0    task is solved                           |
//|             * Info<=0   an error occured:                        |
//|                         -4 means inconvergence of internal SVD   |
//|     P   -   interpolant in barycentric form.                     |
//|     Rep -   report, same format as in LSFitLinearW() subroutine. |
//|             Following fields are set:                            |
//|             * RMSError      rms error on the (X,Y).              |
//|             * AvgError      average error on the (X,Y).          |
//|             * AvgRelError   average relative error on the        |
//|                             non-zero Y                           |
//|             * MaxError      maximum error                        |
//|                             NON-WEIGHTED ERRORS ARE CALCULATED   |
//| NOTES:                                                           |
//|     you can convert P from barycentric form to the power or      |
//|     Chebyshev basis with PolynomialBar2Pow() or                  |
//|     PolynomialBar2Cheb() functions from POLINT subpackage.       |
//+------------------------------------------------------------------+
static void CLSFit::PolynomialFit(double &x[],double &y[],const int n,
                                  const int m,int &info,
                                  CBarycentricInterpolant &p,
                                  CPolynomialFitReport &rep)
  {
//--- create a variable
   int i=0;
//--- create arrays
   double w[];
   double xc[];
   double yc[];
   int    dc[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(m>0,__FUNCTION__+": M<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- allocation
   ArrayResizeAL(w,n);
//--- initialization
   for(i=0;i<=n-1;i++)
      w[i]=1;
//--- function call
   PolynomialFitWC(x,y,w,n,xc,yc,dc,0,m,info,p,rep);
  }
//+------------------------------------------------------------------+
//| Weighted  fitting by polynomials in barycentric form, with       |
//| constraints  on function values or first derivatives.            |
//| Small regularizing term is used when solving constrained tasks   |
//| (to improve stability).                                          |
//| Task is linear, so linear least squares solver is used.          |
//| Complexity of this computational scheme is O(N*M^2), mostly      |
//| dominated by least squares solver                                |
//| SEE ALSO:                                                        |
//|     PolynomialFit()                                              |
//| INPUT PARAMETERS:                                                |
//|     X   -   points, array[0..N-1].                               |
//|     Y   -   function values, array[0..N-1].                      |
//|     W   -   weights, array[0..N-1]                               |
//|             Each summand in square sum of approximation          |
//|             deviations from given values is multiplied by the    |
//|             square of corresponding weight. Fill it by 1's if you|
//|             don't want to solve weighted task.                   |
//|     N   -   number of points, N>0.                               |
//|             * if given, only leading N elements of X/Y/W are used|
//|             * if not given, automatically determined from sizes  |
//|               of X/Y/W                                           |
//|     XC  -   points where polynomial values/derivatives are       |
//|             constrained, array[0..K-1].                          |
//|     YC  -   values of constraints, array[0..K-1]                 |
//|     DC  -   array[0..K-1], types of constraints:                 |
//|             * DC[i]=0   means that P(XC[i])=YC[i]                |
//|             * DC[i]=1   means that P'(XC[i])=YC[i]               |
//|             SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS   |
//|     K   -   number of constraints, 0<=K<M.                       |
//|             K=0 means no constraints (XC/YC/DC are not used in   |
//|             such cases)                                          |
//|     M   -   number of basis functions (= polynomial_degree + 1), |
//|             M>=1                                                 |
//| OUTPUT PARAMETERS:                                               |
//|     Info-   same format as in LSFitLinearW() subroutine:         |
//|             * Info>0    task is solved                           |
//|             * Info<=0   an error occured:                        |
//|                         -4 means inconvergence of internal SVD   |
//|                         -3 means inconsistent constraints        |
//|     P   -   interpolant in barycentric form.                     |
//|     Rep -   report, same format as in LSFitLinearW() subroutine. |
//|             Following fields are set:                            |
//|             * RMSError      rms error on the (X,Y).              |
//|             * AvgError      average error on the (X,Y).          |
//|             * AvgRelError   average relative error on the        |
//|                             non-zero Y                           |
//|             * MaxError      maximum error                        |
//|                             NON-WEIGHTED ERRORS ARE CALCULATED   |
//| IMPORTANT:                                                       |
//|     this subroitine doesn't calculate task's condition number    |
//|     for K<>0.                                                    |
//| NOTES:                                                           |
//|     you can convert P from barycentric form to the power or      |
//|     Chebyshev basis with PolynomialBar2Pow() or                  |
//|     PolynomialBar2Cheb() functions from POLINT subpackage.       |
//| SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:                 |
//| Setting constraints can lead to undesired results, like          |
//| ill-conditioned behavior, or inconsistency being detected.       |
//| From the other side, it allows us to improve quality of the fit. |
//| Here we summarize our experience with constrained regression     |
//| splines:                                                         |
//| * even simple constraints can be inconsistent, see Wikipedia     |
//|   article on this subject:                                       |
//|   http://en.wikipedia.org/wiki/Birkhoff_interpolation            |
//| * the greater is M (given fixed constraints), the more chances   |
//|   that constraints will be consistent                            |
//| * in the general case, consistency of constraints is NOT         |
//|   GUARANTEED.                                                    |
//| * in the one special cases, however, we can guarantee            |
//|   consistency. This case  is:  M>1  and constraints on the       |
//|   function values (NOT DERIVATIVES)                              |
//| Our final recommendation is to use constraints WHEN AND ONLY when|
//| you can't solve your task without them. Anything beyond special  |
//| cases given above is not guaranteed and may result in            |
//| inconsistency.                                                   |
//+------------------------------------------------------------------+
static void CLSFit::PolynomialFitWC(double &cx[],double &cy[],double &cw[],
                                    const int n,double &cxc[],double &cyc[],
                                    int &dc[],const int k,const int m,
                                    int &info,CBarycentricInterpolant &p,
                                    CPolynomialFitReport &rep)
  {
//--- create variables
   double xa=0;
   double xb=0;
   double sa=0;
   double sb=0;
   int    i=0;
   int    j=0;
   double u=0;
   double v=0;
   double s=0;
   int    relcnt=0;
//--- create arrays
   double xoriginal[];
   double yoriginal[];
   double y2[];
   double w2[];
   double tmp[];
   double tmp2[];
   double bx[];
   double by[];
   double bw[];
   double x[];
   double y[];
   double w[];
   double xc[];
   double yc[];
//--- object of class
   CLSFitReport lrep;
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
   ArrayCopy(w,cw);
   ArrayCopy(xc,cxc);
   ArrayCopy(yc,cyc);
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(m>0,__FUNCTION__+": M<=0!"))
      return;
//--- check
   if(!CAp::Assert(k>=0,__FUNCTION__+": K<0!"))
      return;
//--- check
   if(!CAp::Assert(k<m,__FUNCTION__+": K>=M!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(w)>=n,__FUNCTION__+": Length(W)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(xc)>=k,__FUNCTION__+": Length(XC)<K!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(yc)>=k,__FUNCTION__+": Length(YC)<K!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(dc)>=k,__FUNCTION__+": Length(DC)<K!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(w,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(xc,k),__FUNCTION__+": XC contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(yc,k),__FUNCTION__+": YC contains infinite or NaN values!"))
      return;
   for(i=0;i<=k-1;i++)
     {
      //--- check
      if(!CAp::Assert(dc[i]==0 || dc[i]==1,__FUNCTION__+": one of DC[] is not 0 or 1!"))
         return;
     }
//--- Scale X,Y,XC,YC.
//--- Solve scaled problem using internal Chebyshev fitting function.
   LSFitScaleXY(x,y,w,n,xc,yc,dc,k,xa,xb,sa,sb,xoriginal,yoriginal);
   InternalChebyshevFit(x,y,w,n,xc,yc,dc,k,m,info,tmp,lrep);
//--- check
   if(info<0)
      return;
//--- Generate barycentric model and scale it
//--- * BX,BY store barycentric model nodes
//--- * FMatrix is reused (remember - it is at least MxM,what we need)
//--- Model intialization is done in O(M^2). In principle,it can be
//--- done in O(M*log(M)),but before it we solved task with O(N*M^2)
//--- complexity,so it is only a small amount of total time spent.
   ArrayResizeAL(bx,m);
   ArrayResizeAL(by,m);
   ArrayResizeAL(bw,m);
   ArrayResizeAL(tmp2,m);
   s=1;
//--- calculation
   for(i=0;i<=m-1;i++)
     {
      //--- check
      if(m!=1)
         u=MathCos(M_PI*i/(m-1));
      else
         u=0;
      v=0;
      for(j=0;j<=m-1;j++)
        {
         //--- check
         if(j==0)
            tmp2[j]=1;
         else
           {
            //--- check
            if(j==1)
               tmp2[j]=u;
            else
               tmp2[j]=2*u*tmp2[j-1]-tmp2[j-2];
           }
         v=v+tmp[j]*tmp2[j];
        }
      //--- change values
      bx[i]=u;
      by[i]=v;
      bw[i]=s;
      //--- check
      if(i==0 || i==m-1)
         bw[i]=0.5*bw[i];
      s=-s;
     }
//--- function call
   CRatInt::BarycentricBuildXYW(bx,by,bw,m,p);
//--- function call
   CRatInt::BarycentricLinTransX(p,2/(xb-xa),-((xa+xb)/(xb-xa)));
//--- function call
   CRatInt::BarycentricLinTransY(p,sb-sa,sa);
//--- Scale absolute errors obtained from LSFitLinearW.
//--- Relative error should be calculated separately
//--- (because of shifting/scaling of the task)
   rep.m_taskrcond=lrep.m_taskrcond;
   rep.m_rmserror=lrep.m_rmserror*(sb-sa);
   rep.m_avgerror=lrep.m_avgerror*(sb-sa);
   rep.m_maxerror=lrep.m_maxerror*(sb-sa);
   rep.m_avgrelerror=0;
   relcnt=0;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(yoriginal[i]!=0.0)
        {
         rep.m_avgrelerror=rep.m_avgrelerror+MathAbs(CRatInt::BarycentricCalc(p,xoriginal[i])-yoriginal[i])/MathAbs(yoriginal[i]);
         relcnt=relcnt+1;
        }
     }
//--- check
   if(relcnt!=0)
      rep.m_avgrelerror=rep.m_avgrelerror/relcnt;
  }
//+------------------------------------------------------------------+
//| Weghted rational least squares fitting using Floater-Hormann     |
//| rational functions with optimal D chosen from [0,9], with        |
//| constraints and individual weights.                              |
//| Equidistant grid with M node on [min(x),max(x)] is used to build |
//| basis functions. Different values of D are tried, optimal D      |
//| (least WEIGHTED root mean square error) is chosen. Task is       |
//| linear, so linear least squares solver is used. Complexity of    |
//| this computational scheme is O(N*M^2) (mostly dominated by the   |
//| least squares solver).                                           |
//| SEE ALSO                                                         |
//| * BarycentricFitFloaterHormann(), "lightweight" fitting without  |
//|   invididual weights and constraints.                            |
//| INPUT PARAMETERS:                                                |
//|     X   -   points, array[0..N-1].                               |
//|     Y   -   function values, array[0..N-1].                      |
//|     W   -   weights, array[0..N-1]                               |
//|             Each summand in square sum of approximation          |
//|             deviations from given values is multiplied by the    |
//|             square of corresponding weight. Fill it by 1's if    |
//|             you don't want to solve weighted task.               |
//|     N   -   number of points, N>0.                               |
//|     XC  -   points where function values/derivatives are         |
//|             constrained, array[0..K-1].                          |
//|     YC  -   values of constraints, array[0..K-1]                 |
//|     DC  -   array[0..K-1], types of constraints:                 |
//|             * DC[i]=0   means that S(XC[i])=YC[i]                |
//|             * DC[i]=1   means that S'(XC[i])=YC[i]               |
//|             SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS   |
//|     K   -   number of constraints, 0<=K<M.                       |
//|             K=0 means no constraints (XC/YC/DC are not used in   |
//|             such cases)                                          |
//|     M   -   number of basis functions ( = number_of_nodes),      |
//|             M>=2.                                                |
//| OUTPUT PARAMETERS:                                               |
//|     Info-   same format as in LSFitLinearWC() subroutine.        |
//|             * Info>0    task is solved                           |
//|             * Info<=0   an error occured:                        |
//|                         -4 means inconvergence of internal SVD   |
//|                         -3 means inconsistent constraints        |
//|                         -1 means another errors in parameters    |
//|                            passed (N<=0, for example)            |
//|     B   -   barycentric interpolant.                             |
//|     Rep -   report, same format as in LSFitLinearWC() subroutine.|
//|             Following fields are set:                            |
//|             * DBest         best value of the D parameter        |
//|             * RMSError      rms error on the (X,Y).              |
//|             * AvgError      average error on the (X,Y).          |
//|             * AvgRelError   average relative error on the        |
//|                             non-zero Y                           |
//|             * MaxError      maximum error                        |
//|                             NON-WEIGHTED ERRORS ARE CALCULATED   |
//| IMPORTANT:                                                       |
//|     this subroutine doesn't calculate task's condition number    |
//|     for K<>0.                                                    |
//| SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:                 |
//| Setting constraints can lead to undesired results, like          |
//| ill-conditioned behavior, or inconsistency being detected. From  |
//| the other side, it allows us to improve quality of the fit. Here |
//| we summarize our experience with constrained barycentric         |
//| interpolants:                                                    |
//| * excessive constraints can be inconsistent. Floater-Hormann     |
//|   basis functions aren't as flexible as splines (although they   |
//|   are very smooth).                                              |
//| * the more evenly constraints are spread across [min(x),max(x)], |
//|   the more chances that they will be consistent                  |
//| * the greater is M (given  fixed  constraints), the more chances |
//|   that constraints will be consistent                            |
//| * in the general case, consistency of constraints IS NOT         |
//|   GUARANTEED.                                                    |
//| * in the several special cases, however, we CAN guarantee        |
//|   consistency.                                                   |
//| * one of this cases is constraints on the function VALUES at the |
//|   interval boundaries. Note that consustency of the constraints  |
//|   on the function DERIVATIVES is NOT guaranteed (you can use in  |
//|   such cases cubic splines which are more flexible).             |
//| * another special case is ONE constraint on the function value   |
//|   (OR, but not AND, derivative) anywhere in the interval         |
//| Our final recommendation is to use constraints WHEN AND ONLY     |
//| WHEN you can't solve your task without them. Anything beyond     |
//| special cases given above is not guaranteed and may result in    |
//| inconsistency.                                                   |
//+------------------------------------------------------------------+
static void CLSFit::BarycentricFitFloaterHormannWC(double &x[],double &y[],
                                                   double &w[],const int n,
                                                   double &xc[],double &yc[],
                                                   int &dc[],const int k,
                                                   const int m,int &info,
                                                   CBarycentricInterpolant &b,
                                                   CBarycentricFitReport &rep)
  {
//--- create variables
   int    d=0;
   int    i=0;
   double wrmscur=0;
   double wrmsbest=0;
   int    locinfo=0;
//--- objects of classes
   CBarycentricInterpolant locb;
   CBarycentricFitReport locrep;
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(m>0,__FUNCTION__+": M<=0!"))
      return;
//--- check
   if(!CAp::Assert(k>=0,__FUNCTION__+": K<0!"))
      return;
//--- check
   if(!CAp::Assert(k<m,__FUNCTION__+": K>=M!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(w)>=n,__FUNCTION__+": Length(W)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(xc)>=k,__FUNCTION__+": Length(XC)<K!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(yc)>=k,__FUNCTION__+": Length(YC)<K!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(dc)>=k,__FUNCTION__+": Length(DC)<K!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(w,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(xc,k),__FUNCTION__+": XC contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(yc,k),__FUNCTION__+": YC contains infinite or NaN values!"))
      return;
   for(i=0;i<=k-1;i++)
     {
      //--- check
      if(!CAp::Assert(dc[i]==0 || dc[i]==1,__FUNCTION__+": one of DC[] is not 0 or 1!"))
         return;
     }
//--- Find optimal D
//--- Info is -3 by default (degenerate constraints).
//--- If LocInfo will always be equal to -3,Info will remain equal to -3.
//--- If at least once LocInfo will be -4,Info will be -4.
   wrmsbest=CMath::m_maxrealnumber;
   rep.m_dbest=-1;
   info=-3;
//--- calculation
   for(d=0;d<=MathMin(9,n-1);d++)
     {
      //--- function call
      BarycentricFitWCFixedD(x,y,w,n,xc,yc,dc,k,m,d,locinfo,locb,locrep);
      //--- check
      if(!CAp::Assert((locinfo==-4 || locinfo==-3) || locinfo>0,__FUNCTION__+": unexpected result from BarycentricFitWCFixedD!"))
         return;
      //--- check
      if(locinfo>0)
        {
         //--- Calculate weghted RMS
         wrmscur=0;
         for(i=0;i<=n-1;i++)
            wrmscur=wrmscur+CMath::Sqr(w[i]*(y[i]-CRatInt::BarycentricCalc(locb,x[i])));
         wrmscur=MathSqrt(wrmscur/n);
         //--- check
         if(wrmscur<wrmsbest || rep.m_dbest<0)
           {
            //--- function call
            CRatInt::BarycentricCopy(locb,b);
            //--- change values
            rep.m_dbest=d;
            info=1;
            rep.m_rmserror=locrep.m_rmserror;
            rep.m_avgerror=locrep.m_avgerror;
            rep.m_avgrelerror=locrep.m_avgrelerror;
            rep.m_maxerror=locrep.m_maxerror;
            rep.m_taskrcond=locrep.m_taskrcond;
            wrmsbest=wrmscur;
           }
        }
      else
        {
         //--- check
         if(locinfo!=-3 && info<0)
            info=locinfo;
        }
     }
  }
//+------------------------------------------------------------------+
//| Rational least squares fitting using Floater-Hormann rational    |
//| functions with optimal D chosen from [0,9].                      |
//| Equidistant grid with M node on [min(x),max(x)] is used to build |
//| basis functions. Different values of D are tried, optimal D      |
//| (least root mean square error) is chosen.  Task is linear, so    |
//| linear least squares solver is used. Complexity of this          |
//| computational scheme is O(N*M^2) (mostly dominated by the least  |
//| squares solver).                                                 |
//| INPUT PARAMETERS:                                                |
//|     X   -   points, array[0..N-1].                               |
//|     Y   -   function values, array[0..N-1].                      |
//|     N   -   number of points, N>0.                               |
//|     M   -   number of basis functions ( = number_of_nodes), M>=2.|
//| OUTPUT PARAMETERS:                                               |
//|     Info-   same format as in LSFitLinearWC() subroutine.        |
//|             * Info>0    task is solved                           |
//|             * Info<=0   an error occured:                        |
//|                         -4 means inconvergence of internal SVD   |
//|                         -3 means inconsistent constraints        |
//|     B   -   barycentric interpolant.                             |
//|     Rep -   report, same format as in LSFitLinearWC() subroutine.|
//|             Following fields are set:                            |
//|             * DBest         best value of the D parameter        |
//|             * RMSError      rms error on the (X,Y).              |
//|             * AvgError      average error on the (X,Y).          |
//|             * AvgRelError   average relative error on the        |
//|                             non-zero Y                           |
//|             * MaxError      maximum error                        |
//|                             NON-WEIGHTED ERRORS ARE CALCULATED   |
//+------------------------------------------------------------------+
static void CLSFit::BarycentricFitFloaterHormann(double &x[],double &y[],
                                                 const int n,const int m,
                                                 int &info,CBarycentricInterpolant &b,
                                                 CBarycentricFitReport &rep)
  {
//--- create arrays
   double w[];
   double xc[];
   double yc[];
   int    dc[];
//--- create a variable
   int i=0;
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(m>0,__FUNCTION__+": M<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- allocation
   ArrayResizeAL(w,n);
//--- initialization
   for(i=0;i<=n-1;i++)
      w[i]=1;
//--- function call
   BarycentricFitFloaterHormannWC(x,y,w,n,xc,yc,dc,0,m,info,b,rep);
  }
//+------------------------------------------------------------------+
//| Rational least squares fitting using Floater-Hormann rational    |
//| functions with optimal D chosen from [0,9].                      |
//| Equidistant grid with M node on [min(x),max(x)] is used to build |
//| basis functions. Different values of D are tried, optimal D      |
//| (least root mean square error) is chosen. Task is linear, so     |
//| linear least squares solver is used. Complexity of this          |
//| computational scheme is O(N*M^2) (mostly dominated by the least  |
//| squares solver).                                                 |
//| INPUT PARAMETERS:                                                |
//|     X   -   points, array[0..N-1].                               |
//|     Y   -   function values, array[0..N-1].                      |
//|     N   -   number of points, N>0.                               |
//|     M   -   number of basis functions ( = number_of_nodes), M>=2.|
//| OUTPUT PARAMETERS:                                               |
//|     Info-   same format as in LSFitLinearWC() subroutine.        |
//|             * Info>0    task is solved                           |
//|             * Info<=0   an error occured:                        |
//|                         -4 means inconvergence of internal SVD   |
//|                         -3 means inconsistent constraints        |
//|     B   -   barycentric interpolant.                             |
//|     Rep -   report, same format as in LSFitLinearWC() subroutine.|
//|             Following fields are set:                            |
//|             * DBest         best value of the D parameter        |
//|             * RMSError      rms error on the (X,Y).              |
//|             * AvgError      average error on the (X,Y).          |
//|             * AvgRelError   average relative error on the        |
//|                             non-zero Y                           |
//|             * MaxError      maximum error                        |
//|                             NON-WEIGHTED ERRORS ARE CALCULATED   |
//+------------------------------------------------------------------+
static void CLSFit::Spline1DFitPenalized(double &cx[],double &cy[],
                                         const int n,const int m,
                                         const double rho,int &info,
                                         CSpline1DInterpolant &s,
                                         CSpline1DFitReport &rep)
  {
//--- create a variable
   int i=0;
//--- create arrays
   double w[];
   double x[];
   double y[];
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=4,__FUNCTION__+": M<4!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(rho),__FUNCTION__+": Rho is infinite!"))
      return;
//--- allocation
   ArrayResizeAL(w,n);
//--- initialization
   for(i=0;i<=n-1;i++)
      w[i]=1;
//--- function call
   Spline1DFitPenalizedW(x,y,w,n,m,rho,info,s,rep);
  }
//+------------------------------------------------------------------+
//| Weighted fitting by penalized cubic spline.                      |
//| Equidistant grid with M nodes on [min(x,xc),max(x,xc)] is used to|
//| build basis functions. Basis functions are cubic splines with    |
//| natural boundary conditions. Problem is regularized by adding    |
//| non-linearity penalty to the usual least squares penalty         |
//| function:                                                        |
//|     S(x) = arg min { LS + P }, where                             |
//|     LS   = SUM { w[i]^2*(y[i] - S(x[i]))^2 } - least squares     |
//|            penalty                                               |
//|     P    = C*10^rho*integral{ S''(x)^2*dx } - non-linearity      |
//|            penalty                                               |
//|     rho  - tunable constant given by user                        |
//|     C    - automatically determined scale parameter,             |
//|            makes penalty invariant with respect to scaling of X, |
//|            Y, W.                                                 |
//| INPUT PARAMETERS:                                                |
//|     X   -   points, array[0..N-1].                               |
//|     Y   -   function values, array[0..N-1].                      |
//|     W   -   weights, array[0..N-1]                               |
//|             Each summand in square sum of approximation          |
//|             deviations from given values is multiplied by the    |
//|             square of corresponding weight. Fill it by 1's if    |
//|             you don't want to solve weighted problem.            |
//|     N   -   number of points (optional):                         |
//|             * N>0                                                |
//|             * if given, only first N elements of X/Y/W are       |
//|               processed                                          |
//|             * if not given, automatically determined from X/Y/W  |
//|               sizes                                              |
//|     M   -   number of basis functions ( = number_of_nodes), M>=4.|
//|     Rho -   regularization constant passed by user. It penalizes |
//|             nonlinearity in the regression spline. It is         |
//|             logarithmically scaled, i.e. actual value of         |
//|             regularization constant is calculated as 10^Rho. It  |
//|             is automatically scaled so that:                     |
//|             * Rho=2.0 corresponds to moderate amount of          |
//|               nonlinearity                                       |
//|             * generally, it should be somewhere in the           |
//|               [-8.0,+8.0]                                        |
//|             If you do not want to penalize nonlineary,           |
//|             pass small Rho. Values as low as -15 should work.    |
//| OUTPUT PARAMETERS:                                               |
//|     Info-   same format as in LSFitLinearWC() subroutine.        |
//|             * Info>0    task is solved                           |
//|             * Info<=0   an error occured:                        |
//|                         -4 means inconvergence of internal SVD   |
//|                            or Cholesky decomposition; problem    |
//|                            may be too ill-conditioned (very      |
//|                            rare)                                 |
//|     S   -   spline interpolant.                                  |
//|     Rep -   Following fields are set:                            |
//|             * RMSError      rms error on the (X,Y).              |
//|             * AvgError      average error on the (X,Y).          |
//|             * AvgRelError   average relative error on the        |
//|                             non-zero Y                           |
//|             * MaxError      maximum error                        |
//|                             NON-WEIGHTED ERRORS ARE CALCULATED   |
//| IMPORTANT:                                                       |
//|     this subroitine doesn't calculate task's condition number    |
//|     for K<>0.                                                    |
//| NOTE 1: additional nodes are added to the spline outside of the  |
//| fitting interval to force linearity when x<min(x,xc) or          |
//| x>max(x,xc). It is done for consistency - we penalize            |
//| non-linearity at [min(x,xc),max(x,xc)], so it is natural to      |
//| force linearity outside of this interval.                        |
//| NOTE 2: function automatically sorts points, so caller may pass  |
//| unsorted array.                                                  |
//+------------------------------------------------------------------+
static void CLSFit::Spline1DFitPenalizedW(double &cx[],double &cy[],
                                          double &cw[],const int n,
                                          const int m,double rho,
                                          int &info,CSpline1DInterpolant &s,
                                          CSpline1DFitReport &rep)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    b=0;
   double v=0;
   double relcnt=0;
   double xa=0;
   double xb=0;
   double sa=0;
   double sb=0;
   double pdecay=0;
   double tdecay=0;
   double fdmax=0;
   double admax=0;
   double fa=0;
   double ga=0;
   double fb=0;
   double gb=0;
   double lambdav=0;
   int    i_=0;
   int    i1_=0;
//--- create arrays
   double xoriginal[];
   double yoriginal[];
   double fcolumn[];
   double y2[];
   double w2[];
   double xc[];
   double yc[];
   int    dc[];
   double bx[];
   double by[];
   double bd1[];
   double bd2[];
   double tx[];
   double ty[];
   double td[];
   double rightpart[];
   double c[];
   double tmp0[];
   double x[];
   double y[];
   double w[];
//--- create matrix
   CMatrixDouble fmatrix;
   CMatrixDouble amatrix;
   CMatrixDouble d2matrix;
   CMatrixDouble nmatrix;
//--- objects of classes
   CSpline1DInterpolant bs;
   CFblsLinCgState      cgstate;
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
   ArrayCopy(w,cw);
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=4,__FUNCTION__+": M<4!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(w)>=n,__FUNCTION__+": Length(W)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(w,n),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(rho),__FUNCTION__+": Rho is infinite!"))
      return;
//--- Prepare LambdaV
   v=-(MathLog(CMath::m_machineepsilon)/MathLog(10));
//--- check
   if(rho<-v)
      rho=-v;
//--- check
   if(rho>v)
      rho=v;
   lambdav=MathPow(10,rho);
//--- Sort X,Y,W
   CSpline1D::HeapSortDPoints(x,y,w,n);
//--- Scale X,Y,XC,YC
   LSFitScaleXY(x,y,w,n,xc,yc,dc,0,xa,xb,sa,sb,xoriginal,yoriginal);
//--- Allocate space
   fmatrix.Resize(n,m);
   amatrix.Resize(m,m);
   d2matrix.Resize(m,m);
   ArrayResizeAL(bx,m);
   ArrayResizeAL(by,m);
   ArrayResizeAL(fcolumn,n);
   nmatrix.Resize(m,m);
   ArrayResizeAL(rightpart,m);
   ArrayResizeAL(tmp0,MathMax(m,n));
   ArrayResizeAL(c,m);
//--- Fill:
//--- * FMatrix by values of basis functions
//--- * TmpAMatrix by second derivatives of I-th function at J-th point
//--- * CMatrix by constraints
   fdmax=0;
   for(b=0;b<=m-1;b++)
     {
      //--- Prepare I-th basis function
      for(j=0;j<=m-1;j++)
        {
         bx[j]=(double)(2*j)/(double)(m-1)-1;
         by[j]=0;
        }
      by[b]=1;
      //--- function call
      CSpline1D::Spline1DGridDiff2Cubic(bx,by,m,2,0.0,2,0.0,bd1,bd2);
      //--- function call
      CSpline1D::Spline1DBuildCubic(bx,by,m,2,0.0,2,0.0,bs);
      //--- Calculate B-th column of FMatrix
      //--- Update FDMax (maximum column norm)
      CSpline1D::Spline1DConvCubic(bx,by,m,2,0.0,2,0.0,x,n,fcolumn);
      for(i_=0;i_<=n-1;i_++)
         fmatrix[i_].Set(b,fcolumn[i_]);
      v=0;
      for(i=0;i<=n-1;i++)
         v=v+CMath::Sqr(w[i]*fcolumn[i]);
      fdmax=MathMax(fdmax,v);
      //--- Fill temporary with second derivatives of basis function
      for(i_=0;i_<=m-1;i_++)
         d2matrix[b].Set(i_,bd2[i_]);
     }
//--- * calculate penalty matrix A
//--- * calculate max of diagonal elements of A
//--- * calculate PDecay - coefficient before penalty matrix
   for(i=0;i<=m-1;i++)
     {
      for(j=i;j<=m-1;j++)
        {
         //--- calculate integral(B_i''*B_j'') where B_i and B_j are
         //--- i-th and j-th basis splines.
         //--- B_i and B_j are piecewise linear functions.
         v=0;
         for(b=0;b<=m-2;b++)
           {
            //--- change values
            fa=d2matrix[i][b];
            fb=d2matrix[i][b+1];
            ga=d2matrix[j][b];
            gb=d2matrix[j][b+1];
            v=v+(bx[b+1]-bx[b])*(fa*ga+(fa*(gb-ga)+ga*(fb-fa))/2+(fb-fa)*(gb-ga)/3);
           }
         amatrix[i].Set(j,v);
         amatrix[j].Set(i,v);
        }
     }
//--- change values
   admax=0;
   for(i=0;i<=m-1;i++)
      admax=MathMax(admax,MathAbs(amatrix[i][i]));
   pdecay=lambdav*fdmax/admax;
//--- Calculate TDecay for Tikhonov regularization
   tdecay=fdmax*(1+pdecay)*10*CMath::m_machineepsilon;
//--- Prepare system
//--- NOTE: FMatrix is spoiled during this process
   for(i=0;i<=n-1;i++)
     {
      v=w[i];
      for(i_=0;i_<=m-1;i_++)
         fmatrix[i].Set(i_,v*fmatrix[i][i_]);
     }
//--- function call
   CAblas::RMatrixGemm(m,m,n,1.0,fmatrix,0,0,1,fmatrix,0,0,0,0.0,nmatrix,0,0);
   for(i=0;i<=m-1;i++)
     {
      for(j=0;j<=m-1;j++)
         nmatrix[i].Set(j,nmatrix[i][j]+pdecay*amatrix[i][j]);
     }
//--- calculation
   for(i=0;i<=m-1;i++)
      nmatrix[i].Set(i,nmatrix[i][i]+tdecay);
   for(i=0;i<=m-1;i++)
      rightpart[i]=0;
//--- change values
   for(i=0;i<=n-1;i++)
     {
      v=y[i]*w[i];
      for(i_=0;i_<=m-1;i_++)
         rightpart[i_]=rightpart[i_]+v*fmatrix[i][i_];
     }
//--- Solve system
   if(!CTrFac::SPDMatrixCholesky(nmatrix,m,true))
     {
      info=-4;
      return;
     }
//--- function call
   CFbls::FblsCholeskySolve(nmatrix,1.0,m,true,rightpart,tmp0);
//--- copy
   for(i_=0;i_<=m-1;i_++)
      c[i_]=rightpart[i_];
//--- add nodes to force linearity outside of the fitting interval
   CSpline1D::Spline1DGridDiffCubic(bx,c,m,2,0.0,2,0.0,bd1);
//--- allocation
   ArrayResizeAL(tx,m+2);
   ArrayResizeAL(ty,m+2);
   ArrayResizeAL(td,m+2);
//--- copy
   i1_=-1;
   for(i_=1;i_<=m;i_++)
      tx[i_]=bx[i_+i1_];
   i1_=-1;
   for(i_=1;i_<=m;i_++)
      ty[i_]=rightpart[i_+i1_];
   i1_=-1;
   for(i_=1;i_<=m;i_++)
      td[i_]=bd1[i_+i1_];
//--- change values
   tx[0]=tx[1]-(tx[2]-tx[1]);
   ty[0]=ty[1]-td[1]*(tx[2]-tx[1]);
   td[0]=td[1];
   tx[m+1]=tx[m]+(tx[m]-tx[m-1]);
   ty[m+1]=ty[m]+td[m]*(tx[m]-tx[m-1]);
   td[m+1]=td[m];
//--- function call
   CSpline1D::Spline1DBuildHermite(tx,ty,td,m+2,s);
//--- function call
   CSpline1D::Spline1DLinTransX(s,2/(xb-xa),-((xa+xb)/(xb-xa)));
//--- function call
   CSpline1D::Spline1DLinTransY(s,sb-sa,sa);
//--- change value
   info=1;
//--- Fill report
   rep.m_rmserror=0;
   rep.m_avgerror=0;
   rep.m_avgrelerror=0;
   rep.m_maxerror=0;
   relcnt=0;
//--- function call
   CSpline1D::Spline1DConvCubic(bx,rightpart,m,2,0.0,2,0.0,x,n,fcolumn);
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      //--- change values
      v=(sb-sa)*fcolumn[i]+sa;
      rep.m_rmserror=rep.m_rmserror+CMath::Sqr(v-yoriginal[i]);
      rep.m_avgerror=rep.m_avgerror+MathAbs(v-yoriginal[i]);
      //--- check
      if(yoriginal[i]!=0.0)
        {
         rep.m_avgrelerror=rep.m_avgrelerror+MathAbs(v-yoriginal[i])/MathAbs(yoriginal[i]);
         relcnt=relcnt+1;
        }
      rep.m_maxerror=MathMax(rep.m_maxerror,MathAbs(v-yoriginal[i]));
     }
//--- change values
   rep.m_rmserror=MathSqrt(rep.m_rmserror/n);
   rep.m_avgerror=rep.m_avgerror/n;
//--- check
   if(relcnt!=0.0)
      rep.m_avgrelerror=rep.m_avgrelerror/relcnt;
  }
//+------------------------------------------------------------------+
//| Weighted fitting by cubic spline, with constraints on function   |
//| values or derivatives.                                           |
//| Equidistant grid with M-2 nodes on [min(x,xc),max(x,xc)] is used |
//| to build basis functions. Basis functions are cubic splines with |
//| continuous second derivatives and non-fixed first derivatives at |
//| interval ends. Small regularizing term is used when solving      |
//| constrained tasks (to improve stability).                        |
//| Task is linear, so linear least squares solver is used.          |
//| Complexity of this computational scheme is O(N*M^2), mostly      |
//| dominated by least squares solver                                |
//| SEE ALSO                                                         |
//|     Spline1DFitHermiteWC()  -   fitting by Hermite splines (more |
//|                                 flexible, less smooth)           |
//|     Spline1DFitCubic()      -   "lightweight" fitting by cubic   |
//|                                 splines, without invididual      |
//|                                 weights and constraints          |
//| INPUT PARAMETERS:                                                |
//|     X   -   points, array[0..N-1].                               |
//|     Y   -   function values, array[0..N-1].                      |
//|     W   -   weights, array[0..N-1]                               |
//|             Each summand in square sum of approximation          |
//|             deviations from given values is multiplied by the    |
//|             square of corresponding weight. Fill it by 1's if you|
//|             don't want to solve weighted task.                   |
//|     N   -   number of points (optional):                         |
//|             * N>0                                                |
//|             * if given, only first N elements of X/Y/W are       |
//|               processed                                          |
//|             * if not given, automatically determined from X/Y/W  |
//|               sizes                                              |
//|     XC  -   points where spline values/derivatives are           |
//|             constrained, array[0..K-1].                          |
//|     YC  -   values of constraints, array[0..K-1]                 |
//|     DC  -   array[0..K-1], types of constraints:                 |
//|             * DC[i]=0   means that S(XC[i])=YC[i]                |
//|             * DC[i]=1   means that S'(XC[i])=YC[i]               |
//|             SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS   |
//|     K   -   number of constraints (optional):                    |
//|             * 0<=K<M.                                            |
//|             * K=0 means no constraints (XC/YC/DC are not used)   |
//|             * if given, only first K elements of XC/YC/DC are    |
//|               used                                               |
//|             * if not given, automatically determined from        |
//|               XC/YC/DC                                           |
//|     M   -   number of basis functions ( = number_of_nodes+2),    |
//|             M>=4.                                                |
//| OUTPUT PARAMETERS:                                               |
//|     Info-   same format as in LSFitLinearWC() subroutine.        |
//|             * Info>0    task is solved                           |
//|             * Info<=0   an error occured:                        |
//|                         -4 means inconvergence of internal SVD   |
//|                         -3 means inconsistent constraints        |
//|     S   -   spline interpolant.                                  |
//|     Rep -   report, same format as in LSFitLinearWC() subroutine.|
//|             Following fields are set:                            |
//|             * RMSError      rms error on the (X,Y).              |
//|             * AvgError      average error on the (X,Y).          |
//|             * AvgRelError   average relative error on the        |
//|                             non-zero Y                           |
//|             * MaxError      maximum error                        |
//|                             NON-WEIGHTED ERRORS ARE CALCULATED   |
//| IMPORTANT:                                                       |
//|     this subroitine doesn't calculate task's condition number    |
//|     for K<>0.                                                    |
//| ORDER OF POINTS                                                  |
//| Subroutine automatically sorts points, so caller may pass        |
//| unsorted array.                                                  |
//| SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:                 |
//| Setting constraints can lead  to undesired  results, like        |
//| ill-conditioned behavior, or inconsistency being detected. From  |
//| the other side, it allows us to improve quality of the fit.      |
//| Here we summarize our experience with constrained regression     |
//| splines:                                                         |
//| * excessive constraints can be inconsistent. Splines are         |
//|   piecewise cubic functions, and it is easy to create an         |
//|   example, where large number of constraints  concentrated in    |
//|   small area will result in inconsistency. Just because spline   |
//|   is not flexible enough to satisfy all of them. And same        |
//|   constraints spread across the [min(x),max(x)] will be          |
//|   perfectly consistent.                                          |
//| * the more evenly constraints are spread across [min(x),max(x)], |
//|   the more chances that they will be consistent                  |
//| * the greater is M (given fixed constraints), the more chances   |
//|   that constraints will be consistent                            |
//| * in the general case, consistency of constraints IS NOT         |
//|   GUARANTEED.                                                    |
//| * in the several special cases, however, we CAN guarantee        |
//|   consistency.                                                   |
//| * one of this cases is constraints on the function values        |
//|   AND/OR its derivatives at the interval boundaries.             |
//| * another special case is ONE constraint on the function value   |
//|   (OR, but not AND, derivative) anywhere in the interval         |
//| Our final recommendation is to use constraints WHEN AND ONLY WHEN|
//| you can't solve your task without them. Anything beyond special  |
//| cases given above is not guaranteed and may result in            |
//| inconsistency.                                                   |
//+------------------------------------------------------------------+
static void CLSFit::Spline1DFitCubicWC(double &x[],double &y[],double &w[],
                                       const int n,double &xc[],double &yc[],
                                       int &dc[],const int k,const int m,
                                       int &info,CSpline1DInterpolant &s,
                                       CSpline1DFitReport &rep)
  {
//--- create a variable
   int i=0;
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=4,__FUNCTION__+": M<4!"))
      return;
//--- check
   if(!CAp::Assert(k>=0,__FUNCTION__+": K<0!"))
      return;
//--- check
   if(!CAp::Assert(k<m,__FUNCTION__+": K>=M!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(w)>=n,__FUNCTION__+": Length(W)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(xc)>=k,__FUNCTION__+": Length(XC)<K!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(yc)>=k,__FUNCTION__+": Length(YC)<K!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(dc)>=k,__FUNCTION__+": Length(DC)<K!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(w,n),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(xc,k),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(yc,k),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
   for(i=0;i<=k-1;i++)
     {
      //--- check
      if(!CAp::Assert(dc[i]==0 || dc[i]==1,__FUNCTION__+": DC[i] is neither 0 or 1!"))
         return;
     }
//--- function call
   Spline1DFitInternal(0,x,y,w,n,xc,yc,dc,k,m,info,s,rep);
  }
//+------------------------------------------------------------------+
//| Weighted fitting by Hermite spline, with constraints on function |
//| values or first derivatives.                                     |
//| Equidistant grid with M nodes on [min(x,xc),max(x,xc)] is used to|
//| build basis functions. Basis functions are Hermite splines. Small|
//| regularizing term is used when solving constrained tasks (to     |
//| improve stability).                                              |
//| Task is linear, so linear least squares solver is used.          |
//| Complexity of this computational scheme is O(N*M^2), mostly      |
//| dominated by least squares solver                                |
//| SEE ALSO                                                         |
//|     Spline1DFitCubicWC()    -   fitting by Cubic splines (less   |
//|                                 flexible, more smooth)           |
//|     Spline1DFitHermite()    -   "lightweight" Hermite fitting,   |
//|                                 without invididual weights and   |
//|                                 constraints                      |
//| INPUT PARAMETERS:                                                |
//|     X   -   points, array[0..N-1].                               |
//|     Y   -   function values, array[0..N-1].                      |
//|     W   -   weights, array[0..N-1]                               |
//|             Each summand in square sum of approximation          |
//|             deviations from given values is multiplied by the    |
//|             square of corresponding weight. Fill it by 1's if    |
//|             you don't want to solve weighted task.               |
//|     N   -   number of points (optional):                         |
//|             * N>0                                                |
//|             * if given, only first N elements of X/Y/W are       |
//|               processed                                          |
//|             * if not given, automatically determined from X/Y/W  |
//|               sizes                                              |
//|     XC  -   points where spline values/derivatives are           |
//|             constrained, array[0..K-1].                          |
//|     YC  -   values of constraints, array[0..K-1]                 |
//|     DC  -   array[0..K-1], types of constraints:                 |
//|             * DC[i]=0   means that S(XC[i])=YC[i]                |
//|             * DC[i]=1   means that S'(XC[i])=YC[i]               |
//|             SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS   |
//|     K   -   number of constraints (optional):                    |
//|             * 0<=K<M.                                            |
//|             * K=0 means no constraints (XC/YC/DC are not used)   |
//|             * if given, only first K elements of XC/YC/DC are    |
//|               used                                               |
//|             * if not given, automatically determined from        |
//|               XC/YC/DC                                           |
//|     M   -   number of basis functions (= 2 * number of nodes),   |
//|             M>=4,                                                |
//|             M IS EVEN!                                           |
//| OUTPUT PARAMETERS:                                               |
//|     Info-   same format as in LSFitLinearW() subroutine:         |
//|             * Info>0    task is solved                           |
//|             * Info<=0   an error occured:                        |
//|                         -4 means inconvergence of internal SVD   |
//|                         -3 means inconsistent constraints        |
//|                         -2 means odd M was passed (which is not  |
//|                            supported)                            |
//|                         -1 means another errors in parameters    |
//|                            passed (N<=0, for example)            |
//|     S   -   spline interpolant.                                  |
//|     Rep -   report, same format as in LSFitLinearW() subroutine. |
//|             Following fields are set:                            |
//|             * RMSError      rms error on the (X,Y).              |
//|             * AvgError      average error on the (X,Y).          |
//|             * AvgRelError   average relative error on the        |
//|                             non-zero Y                           |
//|             * MaxError      maximum error                        |
//|                             NON-WEIGHTED ERRORS ARE CALCULATED   |
//| IMPORTANT:                                                       |
//|     this subroitine doesn't calculate task's condition number    |
//|     for K<>0.                                                    |
//| IMPORTANT:                                                       |
//|     this subroitine supports only even M's                       |
//| ORDER OF POINTS                                                  |
//| ubroutine automatically sorts points, so caller may pass         |
//| unsorted array.                                                  |
//| SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:                 |
//| Setting constraints can lead to undesired results, like          |
//| ill-conditioned behavior, or inconsistency being detected. From  |
//| the other side, it allows us to improve quality of the fit. Here |
//| we summarize our experience  with constrained regression splines:|
//| * excessive constraints can be inconsistent. Splines are         |
//|   piecewise cubic functions, and it is easy to create an example,|
//|   where large number of constraints concentrated in small area   |
//|   will result in inconsistency. Just because spline is not       |
//|   flexible enough to satisfy all of them. And same constraints   |
//|   spread across the [min(x),max(x)] will be perfectly consistent.|
//| * the more evenly constraints are spread across [min(x),max(x)], |
//|   the more chances that they will be consistent                  |
//| * the greater is M (given  fixed  constraints), the more chances |
//|   that constraints will be consistent                            |
//| * in the general case, consistency of constraints is NOT         |
//|   GUARANTEED.                                                    |
//| * in the several special cases, however, we can guarantee        |
//|   consistency.                                                   |
//| * one of this cases is M>=4 and constraints on the function      |
//|   value (AND/OR its derivative) at the interval boundaries.      |
//| * another special case is M>=4 and ONE constraint on the         |
//|   function value (OR, BUT NOT AND, derivative) anywhere in       |
//|   [min(x),max(x)]                                                |
//| Our final recommendation is to use constraints WHEN AND ONLY when|
//| you can't solve your task without them. Anything beyond  special |
//| cases given above is not guaranteed and may result in            |
//| inconsistency.                                                   |
//+------------------------------------------------------------------+
static void CLSFit::Spline1DFitHermiteWC(double &x[],double &y[],double &w[],
                                         const int n,double &xc[],double &yc[],
                                         int &dc[],const int k,const int m,
                                         int &info,CSpline1DInterpolant &s,
                                         CSpline1DFitReport &rep)
  {
//--- create a variable
   int i=0;
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=4,__FUNCTION__+": M<4!"))
      return;
//--- check
   if(!CAp::Assert(m%2==0,__FUNCTION__+": M is odd!"))
      return;
//--- check
   if(!CAp::Assert(k>=0,__FUNCTION__+": K<0!"))
      return;
//--- check
   if(!CAp::Assert(k<m,__FUNCTION__+": K>=M!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(w)>=n,__FUNCTION__+": Length(W)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(xc)>=k,__FUNCTION__+": Length(XC)<K!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(yc)>=k,__FUNCTION__+": Length(YC)<K!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(dc)>=k,__FUNCTION__+": Length(DC)<K!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(w,n),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(xc,k),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(yc,k),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
   for(i=0;i<=k-1;i++)
     {
      //--- check
      if(!CAp::Assert(dc[i]==0 || dc[i]==1,__FUNCTION__+": DC[i] is neither 0 or 1!"))
         return;
     }
//--- function call
   Spline1DFitInternal(1,x,y,w,n,xc,yc,dc,k,m,info,s,rep);
  }
//+------------------------------------------------------------------+
//| Least squares fitting by cubic spline.                           |
//| This subroutine is "lightweight" alternative for more complex    |
//| and feature - rich Spline1DFitCubicWC(). See Spline1DFitCubicWC()|
//| for more information about subroutine parameters (we don't       |
//| duplicate it here because of length)                             |
//+------------------------------------------------------------------+
static void CLSFit::Spline1DFitCubic(double &x[],double &y[],const int n,
                                     const int m,int &info,
                                     CSpline1DInterpolant &s,
                                     CSpline1DFitReport &rep)
  {
//--- create a variable
   int i=0;
//--- create arrays
   double w[];
   double xc[];
   double yc[];
   int    dc[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=4,__FUNCTION__+": M<4!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- allocation
   ArrayResizeAL(w,n);
//--- initialization
   for(i=0;i<=n-1;i++)
      w[i]=1;
//--- function call
   Spline1DFitCubicWC(x,y,w,n,xc,yc,dc,0,m,info,s,rep);
  }
//+------------------------------------------------------------------+
//| Least squares fitting by Hermite spline.                         |
//| This subroutine is "lightweight" alternative for more complex    |
//| and feature - rich Spline1DFitHermiteWC(). See                   |
//| Spline1DFitHermiteWC() description for more information about    |
//| subroutine parameters (we don't duplicate it here because of     |
//| length).                                                         |
//+------------------------------------------------------------------+
static void CLSFit::Spline1DFitHermite(double &x[],double &y[],const int n,
                                       const int m,int &info,
                                       CSpline1DInterpolant &s,
                                       CSpline1DFitReport &rep)
  {
//--- create a variable
   int i=0;
//--- create arrays
   double w[];
   double xc[];
   double yc[];
   int    dc[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=4,__FUNCTION__+": M<4!"))
      return;
//--- check
   if(!CAp::Assert(m%2==0,__FUNCTION__+": M is odd!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(x)>=n,__FUNCTION__+": Length(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": Length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(x,n),__FUNCTION__+": X contains infinite or NAN values!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NAN values!"))
      return;
//--- allocation
   ArrayResizeAL(w,n);
//--- initialization
   for(i=0;i<=n-1;i++)
      w[i]=1;
//--- function call
   Spline1DFitHermiteWC(x,y,w,n,xc,yc,dc,0,m,info,s,rep);
  }
//+------------------------------------------------------------------+
//| Weighted linear least squares fitting.                           |
//| QR decomposition is used to reduce task to MxM, then triangular  |
//| solver or SVD-based solver is used depending on condition number |
//| of the system. It allows to maximize speed and retain decent     |
//| accuracy.                                                        |
//| INPUT PARAMETERS:                                                |
//|     Y       -   array[0..N-1] Function values in N points.       |
//|     W       -   array[0..N-1] Weights corresponding to function  |
//|                 values. Each summand in square sum of            |
//|                 approximation deviations from given values is    |
//|                 multiplied by the square of corresponding weight.|
//|     FMatrix -   a table of basis functions values,               |
//|                 array[0..N-1, 0..M-1]. FMatrix[I, J] - value of  |
//|                 J-th basis function in I-th point.               |
//|     N       -   number of points used. N>=1.                     |
//|     M       -   number of basis functions, M>=1.                 |
//| OUTPUT PARAMETERS:                                               |
//|     Info    -   error code:                                      |
//|                 * -4    internal SVD decomposition subroutine    |
//|                         failed (very rare and for degenerate     |
//|                         systems only)                            |
//|                 * -1    incorrect N/M were specified             |
//|                 *  1    task is solved                           |
//|     C       -   decomposition coefficients, array[0..M-1]        |
//|     Rep     -   fitting report. Following fields are set:        |
//|                 * Rep.TaskRCond     reciprocal of condition      |
//|                                     number                       |
//|                 * RMSError          rms error on the (X,Y).      |
//|                 * AvgError          average error on the (X,Y).  |
//|                 * AvgRelError       average relative error on the|
//|                                     non-zero Y                   |
//|                 * MaxError          maximum error                |
//|                                     NON-WEIGHTED ERRORS ARE      |
//|                                     CALCULATED                   |
//+------------------------------------------------------------------+
static void CLSFit::LSFitLinearW(double &y[],double &w[],CMatrixDouble &fmatrix,
                                 const int n,const int m,int &info,
                                 double &c[],CLSFitReport &rep)
  {
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(w)>=n,__FUNCTION__+": length(W)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(w,n),__FUNCTION__+": W contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(fmatrix)>=n,__FUNCTION__+": rows(FMatrix)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(fmatrix)>=m,__FUNCTION__+": cols(FMatrix)<M!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(fmatrix,n,m),__FUNCTION__+": FMatrix contains infinite or NaN values!"))
      return;
//--- function call
   LSFitLinearInternal(y,w,fmatrix,n,m,info,c,rep);
  }
//+------------------------------------------------------------------+
//| Weighted constained linear least squares fitting.                |
//| This is variation of LSFitLinearW(), which searchs for           |
//| min|A*x=b| given that K additional constaints C*x=bc are         |
//| satisfied. It reduces original task to modified one: min|B*y-d|  |
//| WITHOUT constraints, then LSFitLinearW() is called.              |
//| INPUT PARAMETERS:                                                |
//|     Y       -   array[0..N-1] Function values in  N  points.     |
//|     W       -   array[0..N-1] Weights corresponding to function  |
//|                 values. Each summand in square sum of            |
//|                 approximation deviations from given values is    |
//|                 multiplied by the square of corresponding        |
//|                 weight.                                          |
//|     FMatrix -   a table of basis functions values,               |
//|                 array[0..N-1, 0..M-1]. FMatrix[I,J] - value of   |
//|                 J-th basis function in I-th point.               |
//|     CMatrix -   a table of constaints, array[0..K-1,0..M].       |
//|                 I-th row of CMatrix corresponds to I-th linear   |
//|                 constraint: CMatrix[I,0]*C[0] + ... +            |
//|                 + CMatrix[I,M-1]*C[M-1] = CMatrix[I,M]           |
//|     N       -   number of points used. N>=1.                     |
//|     M       -   number of basis functions, M>=1.                 |
//|     K       -   number of constraints, 0 <= K < M                |
//|                 K=0 corresponds to absence of constraints.       |
//| OUTPUT PARAMETERS:                                               |
//|     Info    -   error code:                                      |
//|                 * -4    internal SVD decomposition subroutine    |
//|                         failed (very rare and for degenerate     |
//|                         systems only)                            |
//|                 * -3    either too many constraints (M or more), |
//|                         degenerate constraints (some constraints |
//|                         are repetead twice) or inconsistent      |
//|                         constraints were specified.              |
//|                 *  1    task is solved                           |
//|     C       -   decomposition coefficients, array[0..M-1]        |
//|     Rep     -   fitting report. Following fields are set:        |
//|                 * RMSError          rms error on the (X,Y).      |
//|                 * AvgError          average error on the (X,Y).  |
//|                 * AvgRelError       average relative error on the|
//|                                     non-zero Y                   |
//|                 * MaxError          maximum error                |
//|                                     NON-WEIGHTED ERRORS ARE      |
//|                                     CALCULATED                   |
//| IMPORTANT:                                                       |
//|     this subroitine doesn't calculate task's condition number    |
//|     for K<>0.                                                    |
//+------------------------------------------------------------------+
static void CLSFit::LSFitLinearWC(double &cy[],double &w[],CMatrixDouble &fmatrix,
                                  CMatrixDouble &ccmatrix,const int n,
                                  const int m,const int k,int &info,
                                  double &c[],CLSFitReport &rep)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   int    i_=0;
//--- create arrays
   double tau[];
   double tmp[];
   double c0[];
   double y[];
//--- create matrix
   CMatrixDouble q;
   CMatrixDouble f2;
   CMatrixDouble cmatrix;
//--- copy array
   ArrayCopy(y,cy);
//--- copy matrix
   cmatrix=ccmatrix;
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1!"))
      return;
//--- check
   if(!CAp::Assert(k>=0,__FUNCTION__+": K<0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(w)>=n,__FUNCTION__+": length(W)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(w,n),__FUNCTION__+": W contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(fmatrix)>=n,__FUNCTION__+": rows(FMatrix)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(fmatrix)>=m,__FUNCTION__+": cols(FMatrix)<M!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(fmatrix,n,m),__FUNCTION__+": FMatrix contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(cmatrix)>=k,__FUNCTION__+": rows(CMatrix)<K!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(cmatrix)>=m+1 || k==0,__FUNCTION__+": cols(CMatrix)<M+1!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(cmatrix,k,m+1),__FUNCTION__+": CMatrix contains infinite or NaN values!"))
      return;
//--- check
   if(k>=m)
     {
      info=-3;
      return;
     }
//--- Solve
   if(k==0)
     {
      //--- no constraints
      LSFitLinearInternal(y,w,fmatrix,n,m,info,c,rep);
     }
   else
     {
      //--- First,find general form solution of constraints system:
      //--- * factorize C=L*Q
      //--- * unpack Q
      //--- * fill upper part of C with zeros (for RCond)
      //--- We got C=C0+Q2'*y where Q2 is lower M-K rows of Q.
      COrtFac::RMatrixLQ(cmatrix,k,m,tau);
      COrtFac::RMatrixLQUnpackQ(cmatrix,k,m,tau,m,q);
      for(i=0;i<=k-1;i++)
        {
         for(j=i+1;j<=m-1;j++)
            cmatrix[i].Set(j,0.0);
        }
      //--- check
      if(CRCond::RMatrixLURCondInf(cmatrix,k)<1000*CMath::m_machineepsilon)
        {
         info=-3;
         return;
        }
      //--- allocation
      ArrayResizeAL(tmp,k);
      //--- calculation
      for(i=0;i<=k-1;i++)
        {
         //--- check
         if(i>0)
           {
            v=0.0;
            for(i_=0;i_<=i-1;i_++)
               v+=cmatrix[i][i_]*tmp[i_];
           }
         else
            v=0;
         //--- change values
         tmp[i]=(cmatrix[i][m]-v)/cmatrix[i][i];
        }
      //--- allocation
      ArrayResizeAL(c0,m);
      //--- calculation
      for(i=0;i<=m-1;i++)
         c0[i]=0;
      for(i=0;i<=k-1;i++)
        {
         v=tmp[i];
         for(i_=0;i_<=m-1;i_++)
            c0[i_]=c0[i_]+v*q[i][i_];
        }
      //--- Second,prepare modified matrix F2=F*Q2' and solve modified task
      ArrayResizeAL(tmp,MathMax(n,m)+1);
      f2.Resize(n,m-k);
      //--- function call
      CBlas::MatrixVectorMultiply(fmatrix,0,n-1,0,m-1,false,c0,0,m-1,-1.0,y,0,n-1,1.0);
      //--- function call
      CBlas::MatrixMatrixMultiply(fmatrix,0,n-1,0,m-1,false,q,k,m-1,0,m-1,true,1.0,f2,0,n-1,0,m-k-1,0.0,tmp);
      //--- function call
      LSFitLinearInternal(y,w,f2,n,m-k,info,tmp,rep);
      rep.m_taskrcond=-1;
      //--- check
      if(info<=0)
         return;
      //--- then,convert back to original answer: C=C0 + Q2'*Y0
      ArrayResizeAL(c,m);
      for(i_=0;i_<=m-1;i_++)
         c[i_]=c0[i_];
      //--- function call
      CBlas::MatrixVectorMultiply(q,k,m-1,0,m-1,true,tmp,0,m-k-1,1.0,c,0,m-1,1.0);
     }
  }
//+------------------------------------------------------------------+
//| Linear least squares fitting.                                    |
//| QR decomposition is used to reduce task to MxM, then triangular  |
//| solver or SVD-based solver is used depending on condition number |
//| of the system. It allows to maximize speed and retain decent     |
//| accuracy.                                                        |
//| INPUT PARAMETERS:                                                |
//|     Y       -   array[0..N-1] Function values in  N  points.     |
//|     FMatrix -   a table of basis functions values,               |
//|                 array[0..N-1, 0..M-1].                           |
//|                 FMatrix[I, J] - value of J-th basis function in  |
//|                 I-th point.                                      |
//|     N       -   number of points used. N>=1.                     |
//|     M       -   number of basis functions, M>=1.                 |
//| OUTPUT PARAMETERS:                                               |
//|     Info    -   error code:                                      |
//|                 * -4    internal SVD decomposition subroutine    |
//|                         failed (very rare and for degenerate     |
//|                         systems only)                            |
//|                 *  1    task is solved                           |
//|     C       -   decomposition coefficients, array[0..M-1]        |
//|     Rep     -   fitting report. Following fields are set:        |
//|                 * Rep.TaskRCond     reciprocal of condition      |
//|                                     number                       |
//|                 * RMSError          rms error on the (X,Y).      |
//|                 * AvgError          average error on the (X,Y).  |
//|                 * AvgRelError       average relative error on the|
//|                                     non-zero Y                   |
//|                 * MaxError          maximum error                |
//|                                     NON-WEIGHTED ERRORS ARE      |
//|                                     CALCULATED                   |
//+------------------------------------------------------------------+
static void CLSFit::LSFitLinear(double &y[],CMatrixDouble &fmatrix,
                                const int n,const int m,int &info,
                                double &c[],CLSFitReport &rep)
  {
//--- create a variable
   int i=0;
//--- create array
   double w[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(fmatrix)>=n,__FUNCTION__+": rows(FMatrix)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(fmatrix)>=m,__FUNCTION__+": cols(FMatrix)<M!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(fmatrix,n,m),__FUNCTION__+": FMatrix contains infinite or NaN values!"))
      return;
//--- allocation
   ArrayResizeAL(w,n);
//--- initialization
   for(i=0;i<=n-1;i++)
      w[i]=1;
//--- function call
   LSFitLinearInternal(y,w,fmatrix,n,m,info,c,rep);
  }
//+------------------------------------------------------------------+
//| Constained linear least squares fitting.                         |
//| This is variation of LSFitLinear(), which searchs for min|A*x=b| |
//| given that K additional constaints C*x=bc are satisfied. It      |
//| reduces original task to modified one: min|B*y-d| WITHOUT        |
//| constraints, then LSFitLinear() is called.                       |
//| INPUT PARAMETERS:                                                |
//|     Y       -   array[0..N-1] Function values in N points.       |
//|     FMatrix -   a table of basis functions values,               |
//|                 array[0..N-1, 0..M-1]. FMatrix[I,J] - value of   |
//|                 J-th basis function in I-th point.               |
//|     CMatrix -   a table of constaints, array[0..K-1,0..M].       |
//|                 I-th row of CMatrix corresponds to I-th linear   |
//|                 constraint: CMatrix[I,0]*C[0] + ... +            |
//|                 + CMatrix[I,M-1]*C[M-1] = CMatrix[I,M]           |
//|     N       -   number of points used. N>=1.                     |
//|     M       -   number of basis functions, M>=1.                 |
//|     K       -   number of constraints, 0 <= K < M                |
//|                 K=0 corresponds to absence of constraints.       |
//| OUTPUT PARAMETERS:                                               |
//|     Info    -   error code:                                      |
//|                 * -4    internal SVD decomposition subroutine    |
//|                         failed (very rare and for degenerate     |
//|                         systems only)                            |
//|                 * -3    either too many constraints (M or more), |
//|                         degenerate constraints (some constraints |
//|                         are repetead twice) or inconsistent      |
//|                         constraints were specified.              |
//|                 *  1    task is solved                           |
//|     C       -   decomposition coefficients, array[0..M-1]        |
//|     Rep     -   fitting report. Following fields are set:        |
//|                 * RMSError          rms error on the (X,Y).      |
//|                 * AvgError          average error on the (X,Y).  |
//|                 * AvgRelError       average relative error on the|
//|                                     non-zero Y                   |
//|                 * MaxError          maximum error                |
//|                                     NON-WEIGHTED ERRORS ARE      |
//|                                     CALCULATED                   |
//| IMPORTANT:                                                       |
//|     this subroitine doesn't calculate task's condition number    |
//|     for K<>0.                                                    |
//+------------------------------------------------------------------+
static void CLSFit::LSFitLinearC(double &cy[],CMatrixDouble &fmatrix,
                                 CMatrixDouble &cmatrix,const int n,
                                 const int m,const int k,int &info,
                                 double &c[],CLSFitReport &rep)
  {
//--- create a variable
   int i=0;
//--- create arrays
   double w[];
   double y[];
//--- copy array
   ArrayCopy(y,cy);
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1!"))
      return;
//--- check
   if(!CAp::Assert(k>=0,__FUNCTION__+": K<0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(fmatrix)>=n,__FUNCTION__+": rows(FMatrix)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(fmatrix)>=m,__FUNCTION__+": cols(FMatrix)<M!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(fmatrix,n,m),__FUNCTION__+": FMatrix contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(cmatrix)>=k,__FUNCTION__+": rows(CMatrix)<K!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(cmatrix)>=m+1 || k==0,__FUNCTION__+": cols(CMatrix)<M+1!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(cmatrix,k,m+1),__FUNCTION__+": CMatrix contains infinite or NaN values!"))
      return;
//--- allocation
   ArrayResizeAL(w,n);
//--- initialization
   for(i=0;i<=n-1;i++)
      w[i]=1;
//--- function call
   LSFitLinearWC(y,w,fmatrix,cmatrix,n,m,k,info,c,rep);
  }
//+------------------------------------------------------------------+
//| Weighted nonlinear least squares fitting using function values   |
//| only.                                                            |
//| Combination of numerical differentiation and secant updates is   |
//| used to obtain function Jacobian.                                |
//| Nonlinear task min(F(c)) is solved, where                        |
//|     F(c) = (w[0]*(f(c,x[0])-y[0]))^2 + ... +                     |
//|     + (w[n-1]*(f(c,x[n-1])-y[n-1]))^2,                           |
//|     * N is a number of points,                                   |
//|     * M is a dimension of a space points belong to,              |
//|     * K is a dimension of a space of parameters being fitted,    |
//|     * w is an N-dimensional vector of weight coefficients,       |
//|     * x is a set of N points, each of them is an M-dimensional   |
//|       vector,                                                    |
//|     * c is a K-dimensional vector of parameters being fitted     |
//| This subroutine uses only f(c,x[i]).                             |
//| INPUT PARAMETERS:                                                |
//|     X       -   array[0..N-1,0..M-1], points (one row = one      |
//|                 point)                                           |
//|     Y       -   array[0..N-1], function values.                  |
//|     W       -   weights, array[0..N-1]                           |
//|     C       -   array[0..K-1], initial approximation to the      |
//|                 solution,                                        |
//|     N       -   number of points, N>1                            |
//|     M       -   dimension of space                               |
//|     K       -   number of parameters being fitted                |
//|     DiffStep-   numerical differentiation step;                  |
//|                 should not be very small or large;               |
//|                 large = loss of accuracy                         |
//|                 small = growth of round-off errors               |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//+------------------------------------------------------------------+
static void CLSFit::LSFitCreateWF(CMatrixDouble &x,double &y[],double &w[],
                                  double &c[],const int n,const int m,
                                  const int k,const double diffstep,
                                  CLSFitState &state)
  {
//--- create variables
   int i=0;
   int i_=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1!"))
      return;
//--- check
   if(!CAp::Assert(k>=1,__FUNCTION__+": K<1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(c)>=k,__FUNCTION__+": length(C)<K!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(c,k),__FUNCTION__+": C contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(w)>=n,__FUNCTION__+": length(W)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(w,n),__FUNCTION__+": W contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(x)>=n,__FUNCTION__+": rows(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(x)>=m,__FUNCTION__+": cols(X)<M!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(x,n,m),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(diffstep),__FUNCTION__+": DiffStep is not finite!"))
      return;
//--- check
   if(!CAp::Assert((double)(diffstep)>0.0,__FUNCTION__+": DiffStep<=0!"))
      return;
//--- initialization
   state.m_npoints=n;
   state.m_nweights=n;
   state.m_wkind=1;
   state.m_m=m;
   state.m_k=k;
//--- function call
   LSFitSetCond(state,0.0,0.0,0);
//--- function call
   LSFitSetStpMax(state,0.0);
//--- function call
   LSFitSetXRep(state,false);
//--- allocation
   state.m_taskx.Resize(n,m);
   ArrayResizeAL(state.m_tasky,n);
   ArrayResizeAL(state.m_w,n);
   ArrayResizeAL(state.m_c,k);
   ArrayResizeAL(state.m_x,m);
//--- copy
   for(i_=0;i_<=k-1;i_++)
      state.m_c[i_]=c[i_];
   for(i_=0;i_<=n-1;i_++)
      state.m_w[i_]=w[i_];
   for(i=0;i<=n-1;i++)
     {
      for(i_=0;i_<=m-1;i_++)
         state.m_taskx[i].Set(i_,x[i][i_]);
      state.m_tasky[i]=y[i];
     }
//--- allocation
   ArrayResizeAL(state.m_s,k);
   ArrayResizeAL(state.m_bndl,k);
   ArrayResizeAL(state.m_bndu,k);
//--- change values
   for(i=0;i<=k-1;i++)
     {
      state.m_s[i]=1.0;
      state.m_bndl[i]=CInfOrNaN::NegativeInfinity();
      state.m_bndu[i]=CInfOrNaN::PositiveInfinity();
     }
//--- change values
   state.m_optalgo=0;
   state.m_prevnpt=-1;
   state.m_prevalgo=-1;
//--- function call
   CMinLM::MinLMCreateV(k,n,state.m_c,diffstep,state.m_optstate);
//--- function call
   LSFitClearRequestFields(state);
//--- allocation
   ArrayResizeAL(state.m_rstate.ia,5);
   ArrayResizeAL(state.m_rstate.ra,3);
   state.m_rstate.stage=-1;
  }
//+------------------------------------------------------------------+
//| Nonlinear least squares fitting using function values only.      |
//| Combination of numerical differentiation and secant updates is   |
//| used to obtain function Jacobian.                                |
//| Nonlinear task min(F(c)) is solved, where                        |
//|     F(c) = (f(c,x[0])-y[0])^2 + ... + (f(c,x[n-1])-y[n-1])^2,    |
//|     * N is a number of points,                                   |
//|     * M is a dimension of a space points belong to,              |
//|     * K is a dimension of a space of parameters being fitted,    |
//|     * w is an N-dimensional vector of weight coefficients,       |
//|     * x is a set of N points, each of them is an M-dimensional   |
//|       vector,                                                    |
//|     * c is a K-dimensional vector of parameters being fitted     |
//| This subroutine uses only f(c,x[i]).                             |
//| INPUT PARAMETERS:                                                |
//|     X       -   array[0..N-1,0..M-1], points (one row = one      |
//|                 point)                                           |
//|     Y       -   array[0..N-1], function values.                  |
//|     C       -   array[0..K-1], initial approximation to the      |
//|                 solution,                                        |
//|     N       -   number of points, N>1                            |
//|     M       -   dimension of space                               |
//|     K       -   number of parameters being fitted                |
//|     DiffStep-   numerical differentiation step;                  |
//|                 should not be very small or large;               |
//|                 large = loss of accuracy                         |
//|                 small = growth of round-off errors               |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//+------------------------------------------------------------------+
static void CLSFit::LSFitCreateF(CMatrixDouble &x,double &y[],double &c[],
                                 const int n,const int m,const int k,
                                 const double diffstep,CLSFitState &state)
  {
//--- create variables
   int i=0;
   int i_=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1!"))
      return;
//--- check
   if(!CAp::Assert(k>=1,__FUNCTION__+": K<1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(c)>=k,__FUNCTION__+": length(C)<K!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(c,k),__FUNCTION__+": C contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(x)>=n,__FUNCTION__+": rows(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(x)>=m,__FUNCTION__+": cols(X)<M!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(x,n,m),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(x)>=n,__FUNCTION__+": rows(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(x)>=m,__FUNCTION__+": cols(X)<M!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(x,n,m),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(diffstep),__FUNCTION__+": DiffStep is not finite!"))
      return;
//--- check
   if(!CAp::Assert((double)(diffstep)>0.0,__FUNCTION__+": DiffStep<=0!"))
      return;
//--- initialization
   state.m_npoints=n;
   state.m_wkind=0;
   state.m_m=m;
   state.m_k=k;
//--- function call
   LSFitSetCond(state,0.0,0.0,0);
//--- function call
   LSFitSetStpMax(state,0.0);
//--- function call
   LSFitSetXRep(state,false);
//--- allocation
   state.m_taskx.Resize(n,m);
   ArrayResizeAL(state.m_tasky,n);
   ArrayResizeAL(state.m_c,k);
   ArrayResizeAL(state.m_x,m);
//--- copy
   for(i_=0;i_<=k-1;i_++)
      state.m_c[i_]=c[i_];
   for(i=0;i<=n-1;i++)
     {
      for(i_=0;i_<=m-1;i_++)
         state.m_taskx[i].Set(i_,x[i][i_]);
      state.m_tasky[i]=y[i];
     }
//--- allocation
   ArrayResizeAL(state.m_s,k);
   ArrayResizeAL(state.m_bndl,k);
   ArrayResizeAL(state.m_bndu,k);
//--- change values
   for(i=0;i<=k-1;i++)
     {
      state.m_s[i]=1.0;
      state.m_bndl[i]=CInfOrNaN::NegativeInfinity();
      state.m_bndu[i]=CInfOrNaN::PositiveInfinity();
     }
//--- change values
   state.m_optalgo=0;
   state.m_prevnpt=-1;
   state.m_prevalgo=-1;
//--- function call
   CMinLM::MinLMCreateV(k,n,state.m_c,diffstep,state.m_optstate);
//--- function call
   LSFitClearRequestFields(state);
//--- allocation
   ArrayResizeAL(state.m_rstate.ia,5);
   ArrayResizeAL(state.m_rstate.ra,3);
   state.m_rstate.stage=-1;
  }
//+------------------------------------------------------------------+
//| Weighted nonlinear least squares fitting using gradient only.    |
//| Nonlinear task min(F(c)) is solved, where                        |
//|     F(c) = (w[0]*(f(c,x[0])-y[0]))^2 + ... +                     |
//|     + (w[n-1]*(f(c,x[n-1])-y[n-1]))^2,                           |
//|     * N is a number of points,                                   |
//|     * M is a dimension of a space points belong to,              |
//|     * K is a dimension of a space of parameters being fitted,    |
//|     * w is an N-dimensional vector of weight coefficients,       |
//|     * x is a set of N points, each of them is an M-dimensional   |
//|       vector,                                                    |
//|     * c is a K-dimensional vector of parameters being fitted     |
//| This subroutine uses only f(c,x[i]) and its gradient.            |
//| INPUT PARAMETERS:                                                |
//|     X       -   array[0..N-1,0..M-1], points (one row = one      |
//|                 point)                                           |
//|     Y       -   array[0..N-1], function values.                  |
//|     W       -   weights, array[0..N-1]                           |
//|     C       -   array[0..K-1], initial approximation to the      |
//|                 solution,                                        |
//|     N       -   number of points, N>1                            |
//|     M       -   dimension of space                               |
//|     K       -   number of parameters being fitted                |
//|     CheapFG -   boolean flag, which is:                          |
//|                 * True if both function and gradient calculation |
//|                        complexity are less than O(M^2). An       |
//|                        improved algorithm can be used which      |
//|                        corresponds to FGJ scheme from MINLM unit.|
//|                 * False otherwise.                               |
//|                        Standard Jacibian-bases                   |
//|                        Levenberg-Marquardt algo will be used (FJ |
//|                        scheme).                                  |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//| See also:                                                        |
//|     LSFitResults                                                 |
//|     LSFitCreateFG (fitting without weights)                      |
//|     LSFitCreateWFGH (fitting using Hessian)                      |
//|     LSFitCreateFGH (fitting using Hessian, without weights)      |
//+------------------------------------------------------------------+
static void CLSFit::LSFitCreateWFG(CMatrixDouble &x,double &y[],double &w[],
                                   double &c[],const int n,const int m,
                                   const int k,bool cheapfg,CLSFitState &state)
  {
//--- create variables
   int i=0;
   int i_=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1!"))
      return;
//--- check
   if(!CAp::Assert(k>=1,__FUNCTION__+": K<1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(c)>=k,__FUNCTION__+": length(C)<K!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(c,k),__FUNCTION__+": C contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(w)>=n,__FUNCTION__+": length(W)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(w,n),__FUNCTION__+": W contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(x)>=n,__FUNCTION__+": rows(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(x)>=m,__FUNCTION__+": cols(X)<M!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(x,n,m),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- initialization
   state.m_npoints=n;
   state.m_nweights=n;
   state.m_wkind=1;
   state.m_m=m;
   state.m_k=k;
//--- function call
   LSFitSetCond(state,0.0,0.0,0);
//--- function call
   LSFitSetStpMax(state,0.0);
//--- function call
   LSFitSetXRep(state,false);
//--- allocation
   state.m_taskx.Resize(n,m);
   ArrayResizeAL(state.m_tasky,n);
   ArrayResizeAL(state.m_w,n);
   ArrayResizeAL(state.m_c,k);
   ArrayResizeAL(state.m_x,m);
   ArrayResizeAL(state.m_g,k);
//--- copy
   for(i_=0;i_<=k-1;i_++)
      state.m_c[i_]=c[i_];
   for(i_=0;i_<=n-1;i_++)
      state.m_w[i_]=w[i_];
   for(i=0;i<=n-1;i++)
     {
      for(i_=0;i_<=m-1;i_++)
         state.m_taskx[i].Set(i_,x[i][i_]);
      state.m_tasky[i]=y[i];
     }
//--- allocation
   ArrayResizeAL(state.m_s,k);
   ArrayResizeAL(state.m_bndl,k);
   ArrayResizeAL(state.m_bndu,k);
//--- change values
   for(i=0;i<=k-1;i++)
     {
      state.m_s[i]=1.0;
      state.m_bndl[i]=CInfOrNaN::NegativeInfinity();
      state.m_bndu[i]=CInfOrNaN::PositiveInfinity();
     }
//--- change values
   state.m_optalgo=1;
   state.m_prevnpt=-1;
   state.m_prevalgo=-1;
//--- check
   if(cheapfg)
      CMinLM::MinLMCreateVGJ(k,n,state.m_c,state.m_optstate);
   else
      CMinLM::MinLMCreateVJ(k,n,state.m_c,state.m_optstate);
//--- function call
   LSFitClearRequestFields(state);
//--- allocation
   ArrayResizeAL(state.m_rstate.ia,5);
   ArrayResizeAL(state.m_rstate.ra,3);
   state.m_rstate.stage=-1;
  }
//+------------------------------------------------------------------+
//| Nonlinear least squares fitting using gradient only, without     |
//| individual weights.                                              |
//| Nonlinear task min(F(c)) is solved, where                        |
//|     F(c) = ((f(c,x[0])-y[0]))^2 + ... + ((f(c,x[n-1])-y[n-1]))^2,|
//|     * N is a number of points,                                   |
//|     * M is a dimension of a space points belong to,              |
//|     * K is a dimension of a space of parameters being fitted,    |
//|     * x is a set of N points, each of them is an M-dimensional   |
//|       vector,                                                    |
//|     * c is a K-dimensional vector of parameters being fitted     |
//| This subroutine uses only f(c,x[i]) and its gradient.            |
//| INPUT PARAMETERS:                                                |
//|     X       -   array[0..N-1,0..M-1], points (one row = one      |
//|                 point)                                           |
//|     Y       -   array[0..N-1], function values.                  |
//|     C       -   array[0..K-1], initial approximation to the      |
//|                 solution,                                        |
//|     N       -   number of points, N>1                            |
//|     M       -   dimension of space                               |
//|     K       -   number of parameters being fitted                |
//|     CheapFG -   boolean flag, which is:                          |
//|                 * True  if both function and gradient calculation|
//|                         complexity are less than O(M^2). An      |
//|                         improved algorithm can be used which     |
//|                         corresponds to FGJ scheme from MINLM     |
//|                         unit.                                    |
//|                 * False otherwise.                               |
//|                         Standard Jacibian-bases                  |
//|                         Levenberg-Marquardt algo will be used    |
//|                         (FJ scheme).                             |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//+------------------------------------------------------------------+
static void CLSFit::LSFitCreateFG(CMatrixDouble &x,double &y[],double &c[],
                                  const int n,const int m,const int k,
                                  const bool cheapfg,CLSFitState &state)
  {
//--- create variables
   int i=0;
   int i_=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1!"))
      return;
//--- check
   if(!CAp::Assert(k>=1,__FUNCTION__+": K<1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(c)>=k,__FUNCTION__+": length(C)<K!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(c,k),__FUNCTION__+": C contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(x)>=n,__FUNCTION__+": rows(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(x)>=m,__FUNCTION__+": cols(X)<M!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(x,n,m),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(x)>=n,__FUNCTION__+": rows(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(x)>=m,__FUNCTION__+": cols(X)<M!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(x,n,m),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- initialization
   state.m_npoints=n;
   state.m_wkind=0;
   state.m_m=m;
   state.m_k=k;
   LSFitSetCond(state,0.0,0.0,0);
   LSFitSetStpMax(state,0.0);
   LSFitSetXRep(state,false);
//--- allocation
   state.m_taskx.Resize(n,m);
   ArrayResizeAL(state.m_tasky,n);
   ArrayResizeAL(state.m_c,k);
   ArrayResizeAL(state.m_x,m);
   ArrayResizeAL(state.m_g,k);
//--- copy
   for(i_=0;i_<=k-1;i_++)
      state.m_c[i_]=c[i_];
   for(i=0;i<=n-1;i++)
     {
      for(i_=0;i_<=m-1;i_++)
         state.m_taskx[i].Set(i_,x[i][i_]);
      state.m_tasky[i]=y[i];
     }
//--- allocation
   ArrayResizeAL(state.m_s,k);
   ArrayResizeAL(state.m_bndl,k);
   ArrayResizeAL(state.m_bndu,k);
//--- change values
   for(i=0;i<=k-1;i++)
     {
      state.m_s[i]=1.0;
      state.m_bndl[i]=CInfOrNaN::NegativeInfinity();
      state.m_bndu[i]=CInfOrNaN::PositiveInfinity();
     }
//--- change values
   state.m_optalgo=1;
   state.m_prevnpt=-1;
   state.m_prevalgo=-1;
//--- check
   if(cheapfg)
      CMinLM::MinLMCreateVGJ(k,n,state.m_c,state.m_optstate);
   else
      CMinLM::MinLMCreateVJ(k,n,state.m_c,state.m_optstate);
//--- function call
   LSFitClearRequestFields(state);
//--- allocation
   ArrayResizeAL(state.m_rstate.ia,5);
   ArrayResizeAL(state.m_rstate.ra,3);
   state.m_rstate.stage=-1;
  }
//+------------------------------------------------------------------+
//| Weighted nonlinear least squares fitting using gradient/Hessian. |
//| Nonlinear task min(F(c)) is solved, where                        |
//|     F(c) = (w[0]*(f(c,x[0])-y[0]))^2 + ... +                     |
//|     (w[n-1]*(f(c,x[n-1])-y[n-1]))^2,                             |
//|     * N is a number of points,                                   |
//|     * M is a dimension of a space points belong to,              |
//|     * K is a dimension of a space of parameters being fitted,    |
//|     * w is an N-dimensional vector of weight coefficients,       |
//|     * x is a set of N points, each of them is an M-dimensional   |
//|     vector,                                                      |
//|     * c is a K-dimensional vector of parameters being fitted     |
//| This subroutine uses f(c,x[i]), its gradient and its Hessian.    |
//| INPUT PARAMETERS:                                                |
//|     X       -   array[0..N-1,0..M-1], points (one row = one      |
//|                 point)                                           |
//|     Y       -   array[0..N-1], function values.                  |
//|     W       -   weights, array[0..N-1]                           |
//|     C       -   array[0..K-1], initial approximation to the      |
//|                 solution,                                        |
//|     N       -   number of points, N>1                            |
//|     M       -   dimension of space                               |
//|     K       -   number of parameters being fitted                |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//+------------------------------------------------------------------+
static void CLSFit::LSFitCreateWFGH(CMatrixDouble &x,double &y[],double &w[],
                                    double &c[],const int n,const int m,
                                    const int k,CLSFitState &state)
  {
//--- create variables
   int i=0;
   int i_=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1!"))
      return;
//--- check
   if(!CAp::Assert(k>=1,__FUNCTION__+": K<1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(c)>=k,__FUNCTION__+": length(C)<K!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(c,k),__FUNCTION__+": C contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(w)>=n,__FUNCTION__+": length(W)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(w,n),__FUNCTION__+": W contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(x)>=n,__FUNCTION__+": rows(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(x)>=m,__FUNCTION__+": cols(X)<M!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(x,n,m),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- initialization
   state.m_npoints=n;
   state.m_nweights=n;
   state.m_wkind=1;
   state.m_m=m;
   state.m_k=k;
//--- function call
   LSFitSetCond(state,0.0,0.0,0);
//--- function call
   LSFitSetStpMax(state,0.0);
//--- function call
   LSFitSetXRep(state,false);
//--- allocation
   state.m_taskx.Resize(n,m);
   ArrayResizeAL(state.m_tasky,n);
   ArrayResizeAL(state.m_w,n);
   ArrayResizeAL(state.m_c,k);
   state.m_h.Resize(k,k);
   ArrayResizeAL(state.m_x,m);
   ArrayResizeAL(state.m_g,k);
//--- copy
   for(i_=0;i_<=k-1;i_++)
      state.m_c[i_]=c[i_];
   for(i_=0;i_<=n-1;i_++)
      state.m_w[i_]=w[i_];
   for(i=0;i<=n-1;i++)
     {
      for(i_=0;i_<=m-1;i_++)
         state.m_taskx[i].Set(i_,x[i][i_]);
      state.m_tasky[i]=y[i];
     }
//--- allocation
   ArrayResizeAL(state.m_s,k);
   ArrayResizeAL(state.m_bndl,k);
   ArrayResizeAL(state.m_bndu,k);
//--- change values
   for(i=0;i<=k-1;i++)
     {
      state.m_s[i]=1.0;
      state.m_bndl[i]=CInfOrNaN::NegativeInfinity();
      state.m_bndu[i]=CInfOrNaN::PositiveInfinity();
     }
//--- change values
   state.m_optalgo=2;
   state.m_prevnpt=-1;
   state.m_prevalgo=-1;
//--- function call
   CMinLM::MinLMCreateFGH(k,state.m_c,state.m_optstate);
//--- function call
   LSFitClearRequestFields(state);
//--- allocation
   ArrayResizeAL(state.m_rstate.ia,5);
   ArrayResizeAL(state.m_rstate.ra,3);
   state.m_rstate.stage=-1;
  }
//+------------------------------------------------------------------+
//| Nonlinear least squares fitting using gradient/Hessian, without  |
//| individial weights.                                              |
//| Nonlinear task min(F(c)) is solved, where                        |
//|     F(c) = ((f(c,x[0])-y[0]))^2 + ... +                          |
//|     ((f(c,x[n-1])-y[n-1]))^2,                                    |
//|     * N is a number of points,                                   |
//|     * M is a dimension of a space points belong to,              |
//|     * K is a dimension of a space of parameters being fitted,    |
//|     * x is a set of N points, each of them is an M-dimensional   |
//|       vector,                                                    |
//|     * c is a K-dimensional vector of parameters being fitted     |
//| This subroutine uses f(c,x[i]), its gradient and its Hessian.    |
//| INPUT PARAMETERS:                                                |
//|     X       -   array[0..N-1,0..M-1], points (one row = one      |
//|                 point)                                           |
//|     Y       -   array[0..N-1], function values.                  |
//|     C       -   array[0..K-1], initial approximation to the      |
//|                 solution,                                        |
//|     N       -   number of points, N>1                            |
//|     M       -   dimension of space                               |
//|     K       -   number of parameters being fitted                |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which stores algorithm state           |
//+------------------------------------------------------------------+
static void CLSFit::LSFitCreateFGH(CMatrixDouble &x,double &y[],double &c[],
                                   const int n,const int m,const int k,
                                   CLSFitState &state)
  {
//--- create variables
   int i=0;
   int i_=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- check
   if(!CAp::Assert(m>=1,__FUNCTION__+": M<1!"))
      return;
//--- check
   if(!CAp::Assert(k>=1,__FUNCTION__+": K<1!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(c)>=k,__FUNCTION__+": length(C)<K!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(c,k),__FUNCTION__+": C contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(y)>=n,__FUNCTION__+": length(Y)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(y,n),__FUNCTION__+": Y contains infinite or NaN values!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(x)>=n,__FUNCTION__+": rows(X)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(x)>=m,__FUNCTION__+": cols(X)<M!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(x,n,m),__FUNCTION__+": X contains infinite or NaN values!"))
      return;
//--- initialization
   state.m_npoints=n;
   state.m_wkind=0;
   state.m_m=m;
   state.m_k=k;
//--- function call
   LSFitSetCond(state,0.0,0.0,0);
//--- function call
   LSFitSetStpMax(state,0.0);
//--- function call
   LSFitSetXRep(state,false);
//--- allocation
   state.m_taskx.Resize(n,m);
   ArrayResizeAL(state.m_tasky,n);
   ArrayResizeAL(state.m_c,k);
   state.m_h.Resize(k,k);
   ArrayResizeAL(state.m_x,m);
   ArrayResizeAL(state.m_g,k);
//--- copy
   for(i_=0;i_<=k-1;i_++)
      state.m_c[i_]=c[i_];
   for(i=0;i<=n-1;i++)
     {
      for(i_=0;i_<=m-1;i_++)
         state.m_taskx[i].Set(i_,x[i][i_]);
      state.m_tasky[i]=y[i];
     }
//--- allocation
   ArrayResizeAL(state.m_s,k);
   ArrayResizeAL(state.m_bndl,k);
   ArrayResizeAL(state.m_bndu,k);
//--- change values
   for(i=0;i<=k-1;i++)
     {
      state.m_s[i]=1.0;
      state.m_bndl[i]=CInfOrNaN::NegativeInfinity();
      state.m_bndu[i]=CInfOrNaN::PositiveInfinity();
     }
//--- change values
   state.m_optalgo=2;
   state.m_prevnpt=-1;
   state.m_prevalgo=-1;
//--- function call
   CMinLM::MinLMCreateFGH(k,state.m_c,state.m_optstate);
//--- function call
   LSFitClearRequestFields(state);
//--- allocation
   ArrayResizeAL(state.m_rstate.ia,5);
   ArrayResizeAL(state.m_rstate.ra,3);
   state.m_rstate.stage=-1;
  }
//+------------------------------------------------------------------+
//| Stopping conditions for nonlinear least squares fitting.         |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     EpsF    -   stopping criterion. Algorithm stops if           |
//|                 |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}   |
//|     EpsX    -   >=0                                              |
//|                 The subroutine finishes its work if on k+1-th    |
//|                 iteration the condition |v|<=EpsX is fulfilled,  |
//|                 where:                                           |
//|                 * |.| means Euclidian norm                       |
//|                 * v - scaled step vector, v[i]=dx[i]/s[i]        |
//|                 * dx - ste pvector, dx=X(k+1)-X(k)               |
//|                 * s - scaling coefficients set by LSFitSetScale()|
//|     MaxIts  -   maximum number of iterations. If MaxIts=0, the   |
//|                 number of iterations is unlimited. Only          |
//|                 Levenberg-Marquardt iterations are counted       |
//|                 (L-BFGS/CG iterations are NOT counted because    |
//|                 their cost is very low compared to that of LM).  |
//| NOTE                                                             |
//| Passing EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to|
//| automatic stopping criterion selection (according to the scheme  |
//| used by MINLM unit).                                             |
//+------------------------------------------------------------------+
static void CLSFit::LSFitSetCond(CLSFitState &state,const double epsf,
                                 const double epsx,const int maxits)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsf),__FUNCTION__+": EpsF is not finite!"))
      return;
//--- check
   if(!CAp::Assert((double)(epsf)>=0.0,__FUNCTION__+": negative EpsF!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(epsx),__FUNCTION__+": EpsX is not finite!"))
      return;
//--- check
   if(!CAp::Assert((double)(epsx)>=0.0,__FUNCTION__+": negative EpsX!"))
      return;
//--- check
   if(!CAp::Assert(maxits>=0,__FUNCTION__+": negative MaxIts!"))
      return;
//--- change values
   state.m_epsf=epsf;
   state.m_epsx=epsx;
   state.m_maxits=maxits;
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
static void CLSFit::LSFitSetStpMax(CLSFitState &state,const double stpmax)
  {
//--- check
   if(!CAp::Assert(stpmax>=0.0,__FUNCTION__+": StpMax<0!"))
      return;
//--- change value
   state.m_stpmax=stpmax;
  }
//+------------------------------------------------------------------+
//| This function turns on/off reporting.                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure which stores algorithm state           |
//|     NeedXRep-   whether iteration reports are needed or not      |
//| When reports are needed, State.C (current parameters) and State. |
//| F (current value of fitting function) are reported.              |
//+------------------------------------------------------------------+
static void CLSFit::LSFitSetXRep(CLSFitState &state,const bool needxrep)
  {
//--- change value
   state.m_xrep=needxrep;
  }
//+------------------------------------------------------------------+
//| This function sets scaling coefficients for underlying optimizer.|
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
static void CLSFit::LSFitSetScale(CLSFitState &state,double &s[])
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(CAp::Len(s)>=state.m_k,__FUNCTION__+": Length(S)<K"))
      return;
   for(i=0;i<=state.m_k-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(s[i]),__FUNCTION__+": S contains infinite or NAN elements"))
         return;
      //--- check
      if(!CAp::Assert((double)(s[i])!=0.0,__FUNCTION__+": S contains infinite or NAN elements"))
         return;
      //--- change values
      state.m_s[i]=s[i];
     }
  }
//+------------------------------------------------------------------+
//| This function sets boundary constraints for underlying optimizer |
//| Boundary constraints are inactive by default (after initial      |
//| creation). They are preserved until explicitly turned off with   |
//| another SetBC() call.                                            |
//| INPUT PARAMETERS:                                                |
//|     State   -   structure stores algorithm state                 |
//|     BndL    -   lower bounds, array[K].                          |
//|                 If some (all) variables are unbounded, you may   |
//|                 specify very small number or -INF (latter is     |
//|                 recommended because it will allow solver to use  |
//|                 better algorithm).                               |
//|     BndU    -   upper bounds, array[K].                          |
//|                 If some (all) variables are unbounded, you may   |
//|                 specify very large number or +INF (latter is     |
//|                 recommended because it will allow solver to use  |
//|                 better algorithm).                               |
//| NOTE 1: it is possible to specify BndL[i]=BndU[i]. In this case  |
//| I-th variable will be "frozen" at X[i]=BndL[i]=BndU[i].          |
//| NOTE 2: unlike other constrained optimization algorithms, this   |
//| solver has following useful properties:                          |
//| * bound constraints are always satisfied exactly                 |
//| * function is evaluated only INSIDE area specified by bound      |
//|   constraints                                                    |
//+------------------------------------------------------------------+
static void CLSFit::LSFitSetBC(CLSFitState &state,double &bndl[],
                               double &bndu[])
  {
//--- create variables
   int i=0;
   int k=0;
//--- initialization
   k=state.m_k;
//--- check
   if(!CAp::Assert(CAp::Len(bndl)>=k,__FUNCTION__+": Length(BndL)<K"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(bndu)>=k,__FUNCTION__+": Length(BndU)<K"))
      return;
   for(i=0;i<=k-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(bndl[i]) || CInfOrNaN::IsNegativeInfinity(bndl[i]),__FUNCTION__+": BndL contains NAN or +INF"))
         return;
      //--- check
      if(!CAp::Assert(CMath::IsFinite(bndu[i]) || CInfOrNaN::IsPositiveInfinity(bndu[i]),__FUNCTION__+": BndU contains NAN or -INF"))
         return;
      //--- check
      if(CMath::IsFinite(bndl[i]) && CMath::IsFinite(bndu[i]))
        {
         //--- check
         if(!CAp::Assert(bndl[i]<=bndu[i],__FUNCTION__+": BndL[i]>BndU[i]"))
            return;
        }
      //--- change values
      state.m_bndl[i]=bndl[i];
      state.m_bndu[i]=bndu[i];
     }
  }
//+------------------------------------------------------------------+
//| Nonlinear least squares fitting results.                         |
//| Called after return from LSFitFit().                             |
//| INPUT PARAMETERS:                                                |
//|     State   -   algorithm state                                  |
//| OUTPUT PARAMETERS:                                               |
//|     Info    -   completetion code:                               |
//|                     *  1    relative function improvement is no  |
//|                             more than EpsF.                      |
//|                     *  2    relative step is no more than EpsX.  |
//|                     *  4    gradient norm is no more than EpsG   |
//|                     *  5    MaxIts steps was taken               |
//|                     *  7    stopping conditions are too          |
//|                             stringent, further improvement is    |
//|                             impossible                           |
//|     C       -   array[0..K-1], solution                          |
//|     Rep     -   optimization report. Following fields are set:   |
//|                 * Rep.TerminationType completetion code:         |
//|                 * RMSError          rms error on the (X,Y).      |
//|                 * AvgError          average error on the (X,Y).  |
//|                 * AvgRelError       average relative error on the|
//|                                     non-zero Y                   |
//|                 * MaxError          maximum error                |
//|                                     NON-WEIGHTED ERRORS ARE      |
//|                                     CALCULATED                   |
//|                 * WRMSError         weighted rms error on the    |
//|                                     (X,Y).                       |
//+------------------------------------------------------------------+
static void CLSFit::LSFitResults(CLSFitState &state,int &info,double &c[],
                                 CLSFitReport &rep)
  {
//--- create a variable
   int i_=0;
//--- initialization
   info=state.m_repterminationtype;
//--- check
   if(info>0)
     {
      //--- allocation
      ArrayResizeAL(c,state.m_k);
      for(i_=0;i_<=state.m_k-1;i_++)
         c[i_]=state.m_c[i_];
      //--- change values
      rep.m_rmserror=state.m_reprmserror;
      rep.m_wrmserror=state.m_repwrmserror;
      rep.m_avgerror=state.m_repavgerror;
      rep.m_avgrelerror=state.m_repavgrelerror;
      rep.m_maxerror=state.m_repmaxerror;
      rep.m_iterationscount=state.m_repiterationscount;
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine: automatic scaling for LLS tasks.            |
//| NEVER CALL IT DIRECTLY!                                          |
//| Maps abscissas to [-1,1], standartizes ordinates and             |
//| correspondingly scales constraints. It also scales weights so    |
//| that max(W[i])=1                                                 |
//| Transformations performed:                                       |
//| * X, XC         [XA,XB] => [-1,+1]                               |
//|                 transformation makes min(X)=-1, max(X)=+1        |
//| * Y             [SA,SB] => [0,1]                                 |
//|                 transformation makes mean(Y)=0, stddev(Y)=1      |
//| * YC            transformed accordingly to SA, SB, DC[I]         |
//+------------------------------------------------------------------+
static void CLSFit::LSFitScaleXY(double &x[],double &y[],double &w[],
                                 const int n,double &xc[],double &yc[],
                                 int &dc[],const int k,double &xa,
                                 double &xb,double &sa,double &sb,
                                 double &xoriginal[],double &yoriginal[])
  {
//--- create variables
   double xmin=0;
   double xmax=0;
   int    i=0;
   double mx=0;
   int    i_=0;
//--- initialization
   xa=0;
   xb=0;
   sa=0;
   sb=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": incorrect N"))
      return;
//--- check
   if(!CAp::Assert(k>=0,__FUNCTION__+": incorrect K"))
      return;
//--- Calculate xmin/xmax.
//--- Force xmin<>xmax.
   xmin=x[0];
   xmax=x[0];
   for(i=1;i<=n-1;i++)
     {
      xmin=MathMin(xmin,x[i]);
      xmax=MathMax(xmax,x[i]);
     }
   for(i=0;i<=k-1;i++)
     {
      xmin=MathMin(xmin,xc[i]);
      xmax=MathMax(xmax,xc[i]);
     }
//--- check
   if(xmin==xmax)
     {
      //--- check
      if(xmin==0.0)
        {
         xmin=-1;
         xmax=1;
        }
      else
        {
         //--- check
         if(xmin>0.0)
            xmin=0.5*xmin;
         else
            xmax=0.5*xmax;
        }
     }
//--- Transform abscissas: map [XA,XB] to [0,1]
//--- Store old X[] in XOriginal[] (it will be used
//--- to calculate relative error).
   ArrayResizeAL(xoriginal,n);
   for(i_=0;i_<=n-1;i_++)
      xoriginal[i_]=x[i_];
//--- change values
   xa=xmin;
   xb=xmax;
   for(i=0;i<=n-1;i++)
      x[i]=2*(x[i]-0.5*(xa+xb))/(xb-xa);
//--- calculation
   for(i=0;i<=k-1;i++)
     {
      //--- check
      if(!CAp::Assert(dc[i]>=0,__FUNCTION__+": internal error!"))
         return;
      xc[i]=2*(xc[i]-0.5*(xa+xb))/(xb-xa);
      yc[i]=yc[i]*MathPow(0.5*(xb-xa),dc[i]);
     }
//--- Transform function values: map [SA,SB] to [0,1]
//--- SA=mean(Y),
//--- SB=SA+stddev(Y).
//--- Store old Y[] in YOriginal[] (it will be used
//--- to calculate relative error).
   ArrayResizeAL(yoriginal,n);
   for(i_=0;i_<=n-1;i_++)
      yoriginal[i_]=y[i_];
   sa=0;
   for(i=0;i<=n-1;i++)
      sa=sa+y[i];
   sa=sa/n;
//--- change value
   sb=0;
   for(i=0;i<=n-1;i++)
      sb=sb+CMath::Sqr(y[i]-sa);
   sb=MathSqrt(sb/n)+sa;
//--- check
   if(sb==sa)
      sb=2*sa;
//--- check
   if(sb==sa)
      sb=sa+1;
   for(i=0;i<=n-1;i++)
      y[i]=(y[i]-sa)/(sb-sa);
   for(i=0;i<=k-1;i++)
     {
      //--- check
      if(dc[i]==0)
         yc[i]=(yc[i]-sa)/(sb-sa);
      else
         yc[i]=yc[i]/(sb-sa);
     }
//--- Scale weights
   mx=0;
   for(i=0;i<=n-1;i++)
      mx=MathMax(mx,MathAbs(w[i]));
//--- check
   if(mx!=0.0)
     {
      for(i=0;i<=n-1;i++)
         w[i]=w[i]/mx;
     }
  }
//+------------------------------------------------------------------+
//| Internal spline fitting subroutine                               |
//+------------------------------------------------------------------+
static void CLSFit::Spline1DFitInternal(const int st,double &cx[],double &cy[],
                                        double &cw[],const int n,double &cxc[],
                                        double &cyc[],int &dc[],const int k,
                                        const int m,int &info,
                                        CSpline1DInterpolant &s,
                                        CSpline1DFitReport &rep)
  {
//--- create variables
   double v0=0;
   double v1=0;
   double v2=0;
   double mx=0;
   int    i=0;
   int    j=0;
   int    relcnt=0;
   double xa=0;
   double xb=0;
   double sa=0;
   double sb=0;
   double bl=0;
   double br=0;
   double decay=0;
   int    i_=0;
//--- create arrays
   double y2[];
   double w2[];
   double sx[];
   double sy[];
   double sd[];
   double tmp[];
   double xoriginal[];
   double yoriginal[];
   double x[];
   double y[];
   double w[];
   double xc[];
   double yc[];
//--- create matrix
   CMatrixDouble fmatrix;
   CMatrixDouble cmatrix;
//--- objects of classes
   CLSFitReport         lrep;
   CSpline1DInterpolant s2;
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
   ArrayCopy(w,cw);
   ArrayCopy(xc,cxc);
   ArrayCopy(yc,cyc);
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(st==0 || st==1,__FUNCTION__+": internal error!"))
      return;
//--- check
   if(st==0 && m<4)
     {
      info=-1;
      return;
     }
//--- check
   if(st==1 && m<4)
     {
      info=-1;
      return;
     }
//--- check
   if((n<1 || k<0) || k>=m)
     {
      info=-1;
      return;
     }
   for(i=0;i<=k-1;i++)
     {
      info=0;
      //--- check
      if(dc[i]<0)
         info=-1;
      //--- check
      if(dc[i]>1)
         info=-1;
      //--- check
      if(info<0)
         return;
     }
//--- check
   if(st==1 && m%2!=0)
     {
      //--- Hermite fitter must have even number of basis functions
      info=-2;
      return;
     }
//--- weight decay for correct handling of task which becomes
//--- degenerate after constraints are applied
   decay=10000*CMath::m_machineepsilon;
//--- Scale X,Y,XC,YC
   LSFitScaleXY(x,y,w,n,xc,yc,dc,k,xa,xb,sa,sb,xoriginal,yoriginal);
//--- allocate space,initialize:
//--- * SX     -   grid for basis functions
//--- * SY     -   values of basis functions at grid points
//--- * FMatrix-   values of basis functions at X[]
//--- * CMatrix-   values (derivatives) of basis functions at XC[]
   ArrayResizeAL(y2,n+m);
   ArrayResizeAL(w2,n+m);
   fmatrix.Resize(n+m,m);
//--- check
   if(k>0)
      cmatrix.Resize(k,m+1);
//--- check
   if(st==0)
     {
      //--- allocate space for cubic spline
      ArrayResizeAL(sx,m-2);
      ArrayResizeAL(sy,m-2);
      for(j=0;j<=m-2-1;j++)
         sx[j]=(double)(2*j)/(double)(m-2-1)-1;
     }
//--- check
   if(st==1)
     {
      //--- allocate space for Hermite spline
      ArrayResizeAL(sx,m/2);
      ArrayResizeAL(sy,m/2);
      ArrayResizeAL(sd,m/2);
      for(j=0;j<=m/2-1;j++)
         sx[j]=(double)(2*j)/(double)(m/2-1)-1;
     }
//--- Prepare design and constraints matrices:
//--- * fill constraints matrix
//--- * fill first N rows of design matrix with values
//--- * fill next M rows of design matrix with regularizing term
//--- * append M zeros to Y
//--- * append M elements,mean(abs(W)) each,to W
   for(j=0;j<=m-1;j++)
     {
      //--- prepare Jth basis function
      if(st==0)
        {
         //--- cubic spline basis
         for(i=0;i<=m-2-1;i++)
            sy[i]=0;
         bl=0;
         br=0;
         //--- check
         if(j<m-2)
            sy[j]=1;
         //--- check
         if(j==m-2)
            bl=1;
         //--- check
         if(j==m-1)
            br=1;
         //--- function call
         CSpline1D::Spline1DBuildCubic(sx,sy,m-2,1,bl,1,br,s2);
        }
      //--- check
      if(st==1)
        {
         //--- Hermite basis
         for(i=0;i<=m/2-1;i++)
           {
            sy[i]=0;
            sd[i]=0;
           }
         //--- check
         if(j%2==0)
            sy[j/2]=1;
         else
            sd[j/2]=1;
         //--- function call
         CSpline1D::Spline1DBuildHermite(sx,sy,sd,m/2,s2);
        }
      //--- values at X[],XC[]
      for(i=0;i<=n-1;i++)
         fmatrix[i].Set(j,CSpline1D::Spline1DCalc(s2,x[i]));
      for(i=0;i<=k-1;i++)
        {
         //--- check
         if(!CAp::Assert(dc[i]>=0 && dc[i]<=2,__FUNCTION__+": internal error!"))
            return;
         //--- function call
         CSpline1D::Spline1DDiff(s2,xc[i],v0,v1,v2);
         //--- check
         if(dc[i]==0)
            cmatrix[i].Set(j,v0);
         //--- check
         if(dc[i]==1)
            cmatrix[i].Set(j,v1);
         //--- check
         if(dc[i]==2)
            cmatrix[i].Set(j,v2);
        }
     }
//--- calculation
   for(i=0;i<=k-1;i++)
      cmatrix[i].Set(m,yc[i]);
   for(i=0;i<=m-1;i++)
     {
      for(j=0;j<=m-1;j++)
        {
         //--- check
         if(i==j)
            fmatrix[n+i].Set(j,decay);
         else
            fmatrix[n+i].Set(j,0);
        }
     }
//--- allocation
   ArrayResizeAL(y2,n+m);
   ArrayResizeAL(w2,n+m);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      y2[i_]=y[i_];
   for(i_=0;i_<=n-1;i_++)
      w2[i_]=w[i_];
//--- change value
   mx=0;
   for(i=0;i<=n-1;i++)
      mx=mx+MathAbs(w[i]);
   mx=mx/n;
   for(i=0;i<=m-1;i++)
     {
      y2[n+i]=0;
      w2[n+i]=mx;
     }
//--- Solve constrained task
   if(k>0)
     {
      //--- solve using regularization
      LSFitLinearWC(y2,w2,fmatrix,cmatrix,n+m,m,k,info,tmp,lrep);
     }
   else
     {
      //--- no constraints,no regularization needed
      LSFitLinearWC(y,w,fmatrix,cmatrix,n,m,k,info,tmp,lrep);
     }
//--- check
   if(info<0)
      return;
//--- Generate spline and scale it
   if(st==0)
     {
      //--- cubic spline basis
      for(i_=0;i_<=m-2-1;i_++)
         sy[i_]=tmp[i_];
      //--- function call
      CSpline1D::Spline1DBuildCubic(sx,sy,m-2,1,tmp[m-2],1,tmp[m-1],s);
     }
//--- check
   if(st==1)
     {
      //--- Hermite basis
      for(i=0;i<=m/2-1;i++)
        {
         sy[i]=tmp[2*i];
         sd[i]=tmp[2*i+1];
        }
      //--- function call
      CSpline1D::Spline1DBuildHermite(sx,sy,sd,m/2,s);
     }
//--- function call
   CSpline1D::Spline1DLinTransX(s,2/(xb-xa),-((xa+xb)/(xb-xa)));
//--- function call
   CSpline1D::Spline1DLinTransY(s,sb-sa,sa);
//--- Scale absolute errors obtained from LSFitLinearW.
//--- Relative error should be calculated separately
//--- (because of shifting/scaling of the task)
   rep.m_taskrcond=lrep.m_taskrcond;
   rep.m_rmserror=lrep.m_rmserror*(sb-sa);
   rep.m_avgerror=lrep.m_avgerror*(sb-sa);
   rep.m_maxerror=lrep.m_maxerror*(sb-sa);
   rep.m_avgrelerror=0;
   relcnt=0;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(yoriginal[i]!=0.0)
        {
         rep.m_avgrelerror=rep.m_avgrelerror+MathAbs(CSpline1D::Spline1DCalc(s,xoriginal[i])-yoriginal[i])/MathAbs(yoriginal[i]);
         relcnt=relcnt+1;
        }
     }
//--- check
   if(relcnt!=0)
      rep.m_avgrelerror=rep.m_avgrelerror/relcnt;
  }
//+------------------------------------------------------------------+
//| Internal fitting subroutine                                      |
//+------------------------------------------------------------------+
static void CLSFit::LSFitLinearInternal(double &y[],double &w[],
                                        CMatrixDouble &fmatrix,const int n,
                                        const int m,int &info,double &c[],
                                        CLSFitReport &rep)
  {
//--- create variables
   double threshold=0;
   int    i=0;
   int    j=0;
   double v=0;
   int    relcnt=0;
   int    i_=0;
//--- create arrays
   double b[];
   double wmod[];
   double tau[];
   double sv[];
   double tmp[];
   double utb[];
   double sutb[];
//--- create matrix
   CMatrixDouble ft;
   CMatrixDouble q;
   CMatrixDouble l;
   CMatrixDouble r;
   CMatrixDouble u;
   CMatrixDouble vt;
//--- initialization
   info=0;
//--- check
   if(n<1 || m<1)
     {
      info=-1;
      return;
     }
//--- initialization
   info=1;
   threshold=MathSqrt(CMath::m_machineepsilon);
//--- Degenerate case,needs special handling
   if(n<m)
     {
      //--- Create design matrix.
      ft.Resize(n,m);
      ArrayResizeAL(b,n);
      ArrayResizeAL(wmod,n);
      for(j=0;j<=n-1;j++)
        {
         v=w[j];
         for(i_=0;i_<=m-1;i_++)
            ft[j].Set(i_,v*fmatrix[j][i_]);
         //--- change values
         b[j]=w[j]*y[j];
         wmod[j]=1;
        }
      //--- LQ decomposition and reduction to M=N
      ArrayResizeAL(c,m);
      for(i=0;i<=m-1;i++)
         c[i]=0;
      rep.m_taskrcond=0;
      //--- function call
      COrtFac::RMatrixLQ(ft,n,m,tau);
      //--- function call
      COrtFac::RMatrixLQUnpackQ(ft,n,m,tau,n,q);
      //--- function call
      COrtFac::RMatrixLQUnpackL(ft,n,m,l);
      //--- function call
      LSFitLinearInternal(b,wmod,l,n,n,info,tmp,rep);
      //--- check
      if(info<=0)
         return;
      //--- calculation
      for(i=0;i<=n-1;i++)
        {
         v=tmp[i];
         for(i_=0;i_<=m-1;i_++)
            c[i_]=c[i_]+v*q[i][i_];
        }
      //--- exit the function
      return;
     }
//--- N>=M. Generate design matrix and reduce to N=M using
//--- QR decomposition.
   ft.Resize(n,m);
   ArrayResizeAL(b,n);
   for(j=0;j<=n-1;j++)
     {
      v=w[j];
      for(i_=0;i_<=m-1;i_++)
         ft[j].Set(i_,v*fmatrix[j][i_]);
      b[j]=w[j]*y[j];
     }
//--- function call
   COrtFac::RMatrixQR(ft,n,m,tau);
//--- function call
   COrtFac::RMatrixQRUnpackQ(ft,n,m,tau,m,q);
//--- function call
   COrtFac::RMatrixQRUnpackR(ft,n,m,r);
//--- allocation
   ArrayResizeAL(tmp,m);
   for(i=0;i<=m-1;i++)
      tmp[i]=0;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      v=b[i];
      for(i_=0;i_<=m-1;i_++)
         tmp[i_]=tmp[i_]+v*q[i][i_];
     }
//--- allocation
   ArrayResizeAL(b,m);
//--- copy
   for(i_=0;i_<=m-1;i_++)
      b[i_]=tmp[i_];
//--- R contains reduced MxM design upper triangular matrix,
//--- B contains reduced Mx1 right part.
//--- Determine system condition number and decide
//--- should we use triangular solver (faster) or
//--- SVD-based solver (more stable).
//--- We can use LU-based RCond estimator for this task.
   rep.m_taskrcond=CRCond::RMatrixLURCondInf(r,m);
//--- check
   if(rep.m_taskrcond>threshold)
     {
      //--- use QR-based solver
      ArrayResizeAL(c,m);
      c[m-1]=b[m-1]/r[m-1][m-1];
      //--- calculation
      for(i=m-2;i>=0;i--)
        {
         v=0.0;
         for(i_=i+1;i_<=m-1;i_++)
            v+=r[i][i_]*c[i_];
         c[i]=(b[i]-v)/r[i][i];
        }
     }
   else
     {
      //--- use SVD-based solver
      if(!CSingValueDecompose::RMatrixSVD(r,m,m,1,1,2,sv,u,vt))
        {
         info=-4;
         return;
        }
      //--- allocation
      ArrayResizeAL(utb,m);
      ArrayResizeAL(sutb,m);
      for(i=0;i<=m-1;i++)
         utb[i]=0;
      //--- calculation
      for(i=0;i<=m-1;i++)
        {
         v=b[i];
         for(i_=0;i_<=m-1;i_++)
            utb[i_]=utb[i_]+v*u[i][i_];
        }
      //--- check
      if(sv[0]>0.0)
        {
         rep.m_taskrcond=sv[m-1]/sv[0];
         for(i=0;i<=m-1;i++)
           {
            //--- check
            if(sv[i]>threshold*sv[0])
               sutb[i]=utb[i]/sv[i];
            else
               sutb[i]=0;
           }
        }
      else
        {
         //--- change values
         rep.m_taskrcond=0;
         for(i=0;i<=m-1;i++)
            sutb[i]=0;
        }
      //--- allocation
      ArrayResizeAL(c,m);
      for(i=0;i<=m-1;i++)
         c[i]=0;
      //--- calculation
      for(i=0;i<=m-1;i++)
        {
         v=sutb[i];
         for(i_=0;i_<=m-1;i_++)
            c[i_]=c[i_]+v*vt[i][i_];
        }
     }
//--- calculate errors
   rep.m_rmserror=0;
   rep.m_avgerror=0;
   rep.m_avgrelerror=0;
   rep.m_maxerror=0;
   relcnt=0;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      v=0.0;
      for(i_=0;i_<=m-1;i_++)
         v+=fmatrix[i][i_]*c[i_];
      //--- change values
      rep.m_rmserror=rep.m_rmserror+CMath::Sqr(v-y[i]);
      rep.m_avgerror=rep.m_avgerror+MathAbs(v-y[i]);
      //--- check
      if(y[i]!=0.0)
        {
         rep.m_avgrelerror=rep.m_avgrelerror+MathAbs(v-y[i])/MathAbs(y[i]);
         relcnt=relcnt+1;
        }
      rep.m_maxerror=MathMax(rep.m_maxerror,MathAbs(v-y[i]));
     }
//--- change values
   rep.m_rmserror=MathSqrt(rep.m_rmserror/n);
   rep.m_avgerror=rep.m_avgerror/n;
//--- check
   if(relcnt!=0)
      rep.m_avgrelerror=rep.m_avgrelerror/relcnt;
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CLSFit::LSFitClearRequestFields(CLSFitState &state)
  {
//--- change values
   state.m_needf=false;
   state.m_needfg=false;
   state.m_needfgh=false;
   state.m_xupdated=false;
  }
//+------------------------------------------------------------------+
//| Internal subroutine, calculates barycentric basis functions.     |
//| Used for efficient simultaneous calculation of N basis functions.|
//+------------------------------------------------------------------+
static void CLSFit::BarycentricCalcBasis(CBarycentricInterpolant &b,
                                         const double t,double &y[])
  {
//--- create variables
   double s2=0;
   double s=0;
   double v=0;
   int    i=0;
   int    j=0;
   int    i_=0;
//--- special case: N=1
   if(b.m_n==1)
     {
      y[0]=1;
      return;
     }
//--- Here we assume that task is normalized,i.m_e.:
//--- 1. abs(Y[i])<=1
//--- 2. abs(W[i])<=1
//--- 3. X[] is ordered
//--- First,we decide: should we use "safe" formula (guarded
//--- against overflow) or fast one?
   s=MathAbs(t-b.m_x[0]);
   for(i=0;i<=b.m_n-1;i++)
     {
      v=b.m_x[i];
      //--- check
      if(v==t)
        {
         for(j=0;j<=b.m_n-1;j++)
            y[j]=0;
         y[i]=1;
         //--- exit the function
         return;
        }
      //--- change value
      v=MathAbs(t-v);
      //--- check
      if(v<s)
         s=v;
     }
   s2=0;
//--- calculation
   for(i=0;i<=b.m_n-1;i++)
     {
      v=s/(t-b.m_x[i]);
      v=v*b.m_w[i];
      y[i]=v;
      s2=s2+v;
     }
//--- change values
   v=1/s2;
   for(i_=0;i_<=b.m_n-1;i_++)
      y[i_]=v*y[i_];
  }
//+------------------------------------------------------------------+
//| This is internal function for Chebyshev fitting.                 |
//| It assumes that input data are normalized:                       |
//| * X/XC belong to [-1,+1],                                        |
//| * mean(Y)=0, stddev(Y)=1.                                        |
//| It does not checks inputs for errors.                            |
//| This function is used to fit general (shifted) Chebyshev models, | 
//| power basis models or barycentric models.                        |
//| INPUT PARAMETERS:                                                |
//|     X   -   points, array[0..N-1].                               |
//|     Y   -   function values, array[0..N-1].                      |
//|     W   -   weights, array[0..N-1]                               |
//|     N   -   number of points, N>0.                               |
//|     XC  -   points where polynomial values/derivatives are       |
//|             constrained, array[0..K-1].                          |
//|     YC  -   values of constraints, array[0..K-1]                 |
//|     DC  -   array[0..K-1], types of constraints:                 |
//|             * DC[i]=0   means that P(XC[i])=YC[i]                |
//|             * DC[i]=1   means that P'(XC[i])=YC[i]               |
//|     K   -   number of constraints, 0<=K<M.                       |
//|             K=0 means no constraints (XC/YC/DC are not used in   |
//|             such cases)                                          |
//|     M   -   number of basis functions (= polynomial_degree + 1), |
//|             M>=1                                                 |
//| OUTPUT PARAMETERS:                                               |
//|     Info-   same format as in LSFitLinearW() subroutine:         |
//|             * Info>0    task is solved                           |
//|             * Info<=0   an error occured:                        |
//|                         -4 means inconvergence of internal SVD   |
//|                         -3 means inconsistent constraints        |
//|     C   -   interpolant in Chebyshev form; [-1,+1] is used as    |
//|             base interval                                        |
//|     Rep -   report, same format as in LSFitLinearW() subroutine. |
//|             Following fields are set:                            |
//|             * RMSError      rms error on the (X,Y).              |
//|             * AvgError      average error on the (X,Y).          |
//|             * AvgRelError   average relative error on the        |
//|                             non-zero Y                           |
//|             * MaxError      maximum error                        |
//|                             NON-WEIGHTED ERRORS ARE CALCULATED   |
//| IMPORTANT:                                                       |
//|     this subroitine doesn't calculate task's condition number for|
//|     K<>0.                                                        |
//+------------------------------------------------------------------+
static void CLSFit::InternalChebyshevFit(double &x[],double &y[],double &w[],
                                         const int n,double &cxc[],double &cyc[],
                                         int &dc[],const int k,const int m,
                                         int &info,double &c[],CLSFitReport &rep)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double mx=0;
   double decay=0;
   int    i_=0;
//--- create arrays
   double y2[];
   double w2[];
   double tmp[];
   double tmp2[];
   double tmpdiff[];
   double bx[];
   double by[];
   double bw[];
   double xc[];
   double yc[];
//--- create matrix
   CMatrixDouble fmatrix;
   CMatrixDouble cmatrix;
//--- copy arrays
   ArrayCopy(xc,cxc);
   ArrayCopy(yc,cyc);
//--- initialization
   info=0;
//--- weight decay for correct handling of task which becomes
//--- degenerate after constraints are applied
   decay=10000*CMath::m_machineepsilon;
//--- allocate space,initialize/fill:
//--- * FMatrix-   values of basis functions at X[]
//--- * CMatrix-   values (derivatives) of basis functions at XC[]
//--- * fill constraints matrix
//--- * fill first N rows of design matrix with values
//--- * fill next M rows of design matrix with regularizing term
//--- * append M zeros to Y
//--- * append M elements,mean(abs(W)) each,to W
   ArrayResizeAL(y2,n+m);
   ArrayResizeAL(w2,n+m);
   ArrayResizeAL(tmp,m);
   ArrayResizeAL(tmpdiff,m);
   fmatrix.Resize(n+m,m);
//--- check
   if(k>0)
      cmatrix.Resize(k,m+1);
//--- Fill design matrix,Y2,W2:
//--- * first N rows with basis functions for original points
//--- * next M rows with decay terms
   for(i=0;i<=n-1;i++)
     {
      //--- prepare Ith row
      //--- use Tmp for calculations to avoid multidimensional arrays overhead
      for(j=0;j<=m-1;j++)
        {
         //--- check
         if(j==0)
            tmp[j]=1;
         else
           {
            //--- check
            if(j==1)
               tmp[j]=x[i];
            else
               tmp[j]=2*x[i]*tmp[j-1]-tmp[j-2];
           }
        }
      //--- copy
      for(i_=0;i_<=m-1;i_++)
         fmatrix[i].Set(i_,tmp[i_]);
     }
//--- calculation
   for(i=0;i<=m-1;i++)
     {
      for(j=0;j<=m-1;j++)
        {
         //--- check
         if(i==j)
            fmatrix[n+i].Set(j,decay);
         else
            fmatrix[n+i].Set(j,0);
        }
     }
//--- copy
   for(i_=0;i_<=n-1;i_++)
      y2[i_]=y[i_];
   for(i_=0;i_<=n-1;i_++)
      w2[i_]=w[i_];
//--- change value
   mx=0;
   for(i=0;i<=n-1;i++)
      mx=mx+MathAbs(w[i]);
   mx=mx/n;
//--- change values
   for(i=0;i<=m-1;i++)
     {
      y2[n+i]=0;
      w2[n+i]=mx;
     }
//--- fill constraints matrix
   for(i=0;i<=k-1;i++)
     {
      //--- prepare Ith row
      //--- use Tmp for basis function values,
      //--- TmpDiff for basos function derivatives
      for(j=0;j<=m-1;j++)
        {
         //--- check
         if(j==0)
           {
            tmp[j]=1;
            tmpdiff[j]=0;
           }
         else
           {
            //--- check
            if(j==1)
              {
               tmp[j]=xc[i];
               tmpdiff[j]=1;
              }
            else
              {
               tmp[j]=2*xc[i]*tmp[j-1]-tmp[j-2];
               tmpdiff[j]=2*(tmp[j-1]+xc[i]*tmpdiff[j-1])-tmpdiff[j-2];
              }
           }
        }
      //--- check
      if(dc[i]==0)
        {
         for(i_=0;i_<=m-1;i_++)
            cmatrix[i].Set(i_,tmp[i_]);
        }
      //--- check
      if(dc[i]==1)
        {
         for(i_=0;i_<=m-1;i_++)
            cmatrix[i].Set(i_,tmpdiff[i_]);
        }
      cmatrix[i].Set(m,yc[i]);
     }
//--- Solve constrained task
   if(k>0)
     {
      //--- solve using regularization
      LSFitLinearWC(y2,w2,fmatrix,cmatrix,n+m,m,k,info,c,rep);
     }
   else
     {
      //--- no constraints,no regularization needed
      LSFitLinearWC(y,w,fmatrix,cmatrix,n,m,0,info,c,rep);
     }
//--- check
   if(info<0)
      return;
  }
//+------------------------------------------------------------------+
//| Internal Floater-Hormann fitting subroutine for fixed D          |
//+------------------------------------------------------------------+
static void CLSFit::BarycentricFitWCFixedD(double &cx[],double &cy[],
                                           double &cw[],const int n,
                                           double &cxc[],double &cyc[],
                                           int &dc[],const int k,const int m,
                                           const int d,int &info,
                                           CBarycentricInterpolant &b,
                                           CBarycentricFitReport &rep)
  {
//--- create variables
   double v0=0;
   double v1=0;
   double mx=0;
   int    i=0;
   int    j=0;
   int    relcnt=0;
   double xa=0;
   double xb=0;
   double sa=0;
   double sb=0;
   double decay=0;
   int    i_=0;
//--- create arrays
   double y2[];
   double w2[];
   double sx[];
   double sy[];
   double sbf[];
   double xoriginal[];
   double yoriginal[];
   double tmp[];
   double x[];
   double y[];
   double w[];
   double xc[];
   double yc[];
//--- create matrix
   CMatrixDouble fmatrix;
   CMatrixDouble cmatrix;
//--- objects of classes
   CLSFitReport            lrep;
   CBarycentricInterpolant b2;
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
   ArrayCopy(w,cw);
   ArrayCopy(xc,cxc);
   ArrayCopy(yc,cyc);
//--- initialization
   info=0;
//--- check
   if(((n<1 || m<2) || k<0) || k>=m)
     {
      info=-1;
      return;
     }
   for(i=0;i<=k-1;i++)
     {
      info=0;
      //--- check
      if(dc[i]<0)
         info=-1;
      //--- check
      if(dc[i]>1)
         info=-1;
      //--- check
      if(info<0)
         return;
     }
//--- weight decay for correct handling of task which becomes
//--- degenerate after constraints are applied
   decay=10000*CMath::m_machineepsilon;
//--- Scale X,Y,XC,YC
   LSFitScaleXY(x,y,w,n,xc,yc,dc,k,xa,xb,sa,sb,xoriginal,yoriginal);
//--- allocate space,initialize:
//--- * FMatrix-   values of basis functions at X[]
//--- * CMatrix-   values (derivatives) of basis functions at XC[]
   ArrayResizeAL(y2,n+m);
   ArrayResizeAL(w2,n+m);
   fmatrix.Resize(n+m,m);
//--- check
   if(k>0)
      cmatrix.Resize(k,m+1);
//--- allocation
   ArrayResizeAL(y2,n+m);
   ArrayResizeAL(w2,n+m);
//--- Prepare design and constraints matrices:
//--- * fill constraints matrix
//--- * fill first N rows of design matrix with values
//--- * fill next M rows of design matrix with regularizing term
//--- * append M zeros to Y
//--- * append M elements,mean(abs(W)) each,to W
   ArrayResizeAL(sx,m);
   ArrayResizeAL(sy,m);
   ArrayResizeAL(sbf,m);
   for(j=0;j<=m-1;j++)
      sx[j]=(double)(2*j)/(double)(m-1)-1;
   for(i=0;i<=m-1;i++)
      sy[i]=1;
//--- function call
   CRatInt::BarycentricBuildFloaterHormann(sx,sy,m,d,b2);
//--- change value
   mx=0;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      //--- function call
      BarycentricCalcBasis(b2,x[i],sbf);
      for(i_=0;i_<=m-1;i_++)
         fmatrix[i].Set(i_,sbf[i_]);
      //--- change values
      y2[i]=y[i];
      w2[i]=w[i];
      mx=mx+MathAbs(w[i])/n;
     }
//--- calculation
   for(i=0;i<=m-1;i++)
     {
      for(j=0;j<=m-1;j++)
        {
         //--- check
         if(i==j)
            fmatrix[n+i].Set(j,decay);
         else
            fmatrix[n+i].Set(j,0);
        }
      //--- change values
      y2[n+i]=0;
      w2[n+i]=mx;
     }
//--- check
   if(k>0)
     {
      for(j=0;j<=m-1;j++)
        {
         for(i=0;i<=m-1;i++)
            sy[i]=0;
         sy[j]=1;
         //--- function call
         CRatInt::BarycentricBuildFloaterHormann(sx,sy,m,d,b2);
         //--- calculation
         for(i=0;i<=k-1;i++)
           {
            //--- check
            if(!CAp::Assert(dc[i]>=0 && dc[i]<=1,__FUNCTION__+": internal error!"))
               return;
            //--- function call
            CRatInt::BarycentricDiff1(b2,xc[i],v0,v1);
            //--- check
            if(dc[i]==0)
               cmatrix[i].Set(j,v0);
            //--- check
            if(dc[i]==1)
               cmatrix[i].Set(j,v1);
           }
        }
      for(i=0;i<=k-1;i++)
         cmatrix[i].Set(m,yc[i]);
     }
//--- Solve constrained task
   if(k>0)
     {
      //--- solve using regularization
      LSFitLinearWC(y2,w2,fmatrix,cmatrix,n+m,m,k,info,tmp,lrep);
     }
   else
     {
      //--- no constraints,no regularization needed
      LSFitLinearWC(y,w,fmatrix,cmatrix,n,m,k,info,tmp,lrep);
     }
//--- check
   if(info<0)
      return;
//--- Generate interpolant and scale it
   for(i_=0;i_<=m-1;i_++)
      sy[i_]=tmp[i_];
//--- function call
   CRatInt::BarycentricBuildFloaterHormann(sx,sy,m,d,b);
//--- function call
   CRatInt::BarycentricLinTransX(b,2/(xb-xa),-((xa+xb)/(xb-xa)));
//--- function call
   CRatInt::BarycentricLinTransY(b,sb-sa,sa);
//--- Scale absolute errors obtained from LSFitLinearW.
//--- Relative error should be calculated separately
//--- (because of shifting/scaling of the task)
   rep.m_taskrcond=lrep.m_taskrcond;
   rep.m_rmserror=lrep.m_rmserror*(sb-sa);
   rep.m_avgerror=lrep.m_avgerror*(sb-sa);
   rep.m_maxerror=lrep.m_maxerror*(sb-sa);
   rep.m_avgrelerror=0;
   relcnt=0;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(yoriginal[i]!=0.0)
        {
         rep.m_avgrelerror=rep.m_avgrelerror+MathAbs(CRatInt::BarycentricCalc(b,xoriginal[i])-yoriginal[i])/MathAbs(yoriginal[i]);
         relcnt=relcnt+1;
        }
     }
//--- check
   if(relcnt!=0)
      rep.m_avgrelerror=rep.m_avgrelerror/relcnt;
  }
//+------------------------------------------------------------------+
//| NOTES:                                                           |
//| 1. this algorithm is somewhat unusual because it works with      |
//|    parameterized function f(C,X), where X is a function argument |
//|    (we have many points which are characterized by different     |
//|    argument values), and C is a parameter to fit.                |
//|    For example, if we want to do linear fit by f(c0,c1,x) =      |
//|    = c0*x+c1, then x will be argument, and {c0,c1} will be       |
//|    parameters.                                                   |
//|    It is important to understand that this algorithm finds       |
//|    minimum in the space of function PARAMETERS (not arguments),  | 
//|    so it needs derivatives of f() with respect to C, not X.      |
//|    In the example above it will need f=c0*x+c1 and               |
//|    {df/dc0,df/dc1} = {x,1} instead of {df/dx} = {c0}.            |
//| 2. Callback functions accept C as the first parameter, and X as  |
//|    the second                                                    |
//| 3. If state was created with LSFitCreateFG(), algorithm needs    |
//|    just function and its gradient, but if state was created with |
//|    LSFitCreateFGH(), algorithm will need function, gradient and  |
//|    Hessian.                                                      |
//|    According to the said above, there ase several versions of    |
//|    this function, which accept different sets of callbacks.      |
//|    This flexibility opens way to subtle errors - you may create  |
//|    state with LSFitCreateFGH() (optimization using Hessian), but |
//|    call function which does not accept Hessian. So when algorithm|
//|    will request Hessian, there will be no callback to call. In   |
//|    this case exception will be thrown.                           |
//|    Be careful to avoid such errors because there is no way to    |
//|    find them at compile time - you can see them at runtime only. |
//+------------------------------------------------------------------+
static bool CLSFit::LSFitIteration(CLSFitState &state)
  {
//--- create variables
   int    n=0;
   int    m=0;
   int    k=0;
   int    i=0;
   int    j=0;
   double v=0;
   double vv=0;
   double relcnt=0;
   int    i_=0;
//--- Reverse communication preparations
//--- I know it looks ugly,but it works the same way
//--- anywhere from C++ to Python.
//
//--- This code initializes locals by:
//--- * random values determined during code
//---   generation - on first subroutine call
//--- * values from previous call - on subsequent calls
   if(state.m_rstate.stage>=0)
     {
      //--- initialization
      n=state.m_rstate.ia[0];
      m=state.m_rstate.ia[1];
      k=state.m_rstate.ia[2];
      i=state.m_rstate.ia[3];
      j=state.m_rstate.ia[4];
      v=state.m_rstate.ra[0];
      vv=state.m_rstate.ra[1];
      relcnt=state.m_rstate.ra[2];
     }
   else
     {
      //--- initialization
      n=-983;
      m=-989;
      k=-834;
      i=900;
      j=-287;
      v=364;
      vv=214;
      relcnt=-338;
     }
//--- check
   if(state.m_rstate.stage==0)
     {
      //--- change value
      state.m_needf=false;
      //--- check
      if(state.m_wkind==1)
         vv=state.m_w[i];
      else
         vv=1.0;
      //--- change values
      state.m_optstate.m_fi[i]=vv*(state.m_f-state.m_tasky[i]);
      i=i+1;
      //--- function call, return result
      return(Func_lbl_11(state,n,m,k,i,j,v,vv,relcnt));
     }
//--- check
   if(state.m_rstate.stage==1)
     {
      //--- change value
      state.m_needf=false;
      //--- check
      if(state.m_wkind==1)
         vv=state.m_w[i];
      else
         vv=1.0;
      //--- change values
      state.m_optstate.m_f=state.m_optstate.m_f+CMath::Sqr(vv*(state.m_f-state.m_tasky[i]));
      i=i+1;
      //--- function call, return result
      return(Func_lbl_16(state,n,m,k,i,j,v,vv,relcnt));
     }
//--- check
   if(state.m_rstate.stage==2)
     {
      //--- change value
      state.m_needfg=false;
      //--- check
      if(state.m_wkind==1)
         vv=state.m_w[i];
      else
         vv=1.0;
      //--- change values
      state.m_optstate.m_f=state.m_optstate.m_f+CMath::Sqr(vv*(state.m_f-state.m_tasky[i]));
      v=CMath::Sqr(vv)*2*(state.m_f-state.m_tasky[i]);
      for(i_=0;i_<=k-1;i_++)
         state.m_optstate.m_g[i_]=state.m_optstate.m_g[i_]+v*state.m_g[i_];
      i=i+1;
      //--- function call, return result
      return(Func_lbl_21(state,n,m,k,i,j,v,vv,relcnt));
     }
//--- check
   if(state.m_rstate.stage==3)
     {
      //--- change value
      state.m_needfg=false;
      //--- check
      if(state.m_wkind==1)
         vv=state.m_w[i];
      else
         vv=1.0;
      //--- change values
      state.m_optstate.m_fi[i]=vv*(state.m_f-state.m_tasky[i]);
      for(i_=0;i_<=k-1;i_++)
         state.m_optstate.m_j[i].Set(i_,vv*state.m_g[i_]);
      i=i+1;
      //--- function call, return result
      return(Func_lbl_26(state,n,m,k,i,j,v,vv,relcnt));
     }
//--- check
   if(state.m_rstate.stage==4)
     {
      //--- change value
      state.m_needfgh=false;
      //--- check
      if(state.m_wkind==1)
         vv=state.m_w[i];
      else
         vv=1.0;
      //--- change values
      state.m_optstate.m_f=state.m_optstate.m_f+CMath::Sqr(vv*(state.m_f-state.m_tasky[i]));
      v=CMath::Sqr(vv)*2*(state.m_f-state.m_tasky[i]);
      for(i_=0;i_<=k-1;i_++)
         state.m_optstate.m_g[i_]=state.m_optstate.m_g[i_]+v*state.m_g[i_];
      //--- calculation
      for(j=0;j<=k-1;j++)
        {
         v=2*CMath::Sqr(vv)*state.m_g[j];
         for(i_=0;i_<=k-1;i_++)
            state.m_optstate.m_h[j].Set(i_,state.m_optstate.m_h[j][i_]+v*state.m_g[i_]);
         v=2*CMath::Sqr(vv)*(state.m_f-state.m_tasky[i]);
         for(i_=0;i_<=k-1;i_++)
            state.m_optstate.m_h[j].Set(i_,state.m_optstate.m_h[j][i_]+v*state.m_h[j][i_]);
        }
      i=i+1;
      //--- function call, return result
      return(Func_lbl_31(state,n,m,k,i,j,v,vv,relcnt));
     }
//--- check
   if(state.m_rstate.stage==5)
     {
      //--- change value
      state.m_xupdated=false;
      //--- function call, return result
      return(Func_lbl_7(state,n,m,k,i,j,v,vv,relcnt));
     }
//--- check
   if(state.m_rstate.stage==6)
     {
      //--- change values
      state.m_needf=false;
      v=state.m_f;
      //--- check
      if(state.m_wkind==1)
         vv=state.m_w[i];
      else
         vv=1.0;
      //--- change values
      state.m_reprmserror=state.m_reprmserror+CMath::Sqr(v-state.m_tasky[i]);
      state.m_repwrmserror=state.m_repwrmserror+CMath::Sqr(vv*(v-state.m_tasky[i]));
      state.m_repavgerror=state.m_repavgerror+MathAbs(v-state.m_tasky[i]);
      //--- check
      if(state.m_tasky[i]!=0.0)
        {
         state.m_repavgrelerror=state.m_repavgrelerror+MathAbs(v-state.m_tasky[i])/MathAbs(state.m_tasky[i]);
         relcnt=relcnt+1;
        }
      //--- change values
      state.m_repmaxerror=MathMax(state.m_repmaxerror,MathAbs(v-state.m_tasky[i]));
      i=i+1;
      //--- function call, return result
      return(Func_lbl_38(state,n,m,k,i,j,v,vv,relcnt));
     }
//--- Routine body
//--- init
   if(state.m_wkind==1)
     {
      //--- check
      if(!CAp::Assert(state.m_npoints==state.m_nweights,__FUNCTION__+": number of points is not equal to the number of weights"))
         return(false);
     }
//--- change values
   n=state.m_npoints;
   m=state.m_m;
   k=state.m_k;
//--- function call
   CMinLM::MinLMSetCond(state.m_optstate,0.0,state.m_epsf,state.m_epsx,state.m_maxits);
//--- function call
   CMinLM::MinLMSetStpMax(state.m_optstate,state.m_stpmax);
//--- function call
   CMinLM::MinLMSetXRep(state.m_optstate,state.m_xrep);
//--- function call
   CMinLM::MinLMSetScale(state.m_optstate,state.m_s);
//--- function call
   CMinLM::MinLMSetBC(state.m_optstate,state.m_bndl,state.m_bndu);
//--- Optimize
   return(Func_lbl_7(state,n,m,k,i,j,v,vv,relcnt));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for LSFitIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static void CLSFit::Func_lbl_rcomm(CLSFitState &state,int n,int m,int k,int i,
                                   int j,double v,double vv,double relcnt)
  {
//--- Saving state
   state.m_rstate.ia[0]=n;
   state.m_rstate.ia[1]=m;
   state.m_rstate.ia[2]=k;
   state.m_rstate.ia[3]=i;
   state.m_rstate.ia[4]=j;
   state.m_rstate.ra[0]=v;
   state.m_rstate.ra[1]=vv;
   state.m_rstate.ra[2]=relcnt;
  }
//+------------------------------------------------------------------+
//| Auxiliary function for LSFitIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CLSFit::Func_lbl_7(CLSFitState &state,int &n,int &m,int &k,int &i,
                               int &j,double &v,double &vv,double &relcnt)
  {
//--- check
   if(!CMinLM::MinLMIteration(state.m_optstate))
     {
      //--- function call
      CMinLM::MinLMResults(state.m_optstate,state.m_c,state.m_optrep);
      //--- change values
      state.m_repterminationtype=state.m_optrep.m_terminationtype;
      state.m_repiterationscount=state.m_optrep.m_iterationscount;
      //--- calculate errors
      if(state.m_repterminationtype<=0)
         return(false);
      //--- change values
      state.m_reprmserror=0;
      state.m_repwrmserror=0;
      state.m_repavgerror=0;
      state.m_repavgrelerror=0;
      state.m_repmaxerror=0;
      relcnt=0;
      i=0;
      //--- function call, return result
      return(Func_lbl_38(state,n,m,k,i,j,v,vv,relcnt));
     }
//--- check
   if(!state.m_optstate.m_needfi)
     {
      //--- check
      if(!state.m_optstate.m_needf)
         return(Func_lbl_14(state,n,m,k,i,j,v,vv,relcnt));
      //--- calculate F=sum (wi*(f(xi,c)-yi))^2
      state.m_optstate.m_f=0;
      i=0;
      //--- function call, return result
      return(Func_lbl_16(state,n,m,k,i,j,v,vv,relcnt));
     }
//--- calculate f[]=wi*(f(xi,c)-yi)
   i=0;
//--- function call, return result
   return(Func_lbl_11(state,n,m,k,i,j,v,vv,relcnt));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for LSFitIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CLSFit::Func_lbl_11(CLSFitState &state,int &n,int &m,int &k,int &i,
                                int &j,double &v,double &vv,double &relcnt)
  {
//--- check
   if(i>n-1)
      return(Func_lbl_7(state,n,m,k,i,j,v,vv,relcnt));
//--- copy
   for(int i_=0;i_<=k-1;i_++)
      state.m_c[i_]=state.m_optstate.m_x[i_];
   for(int i_=0;i_<=m-1;i_++)
      state.m_x[i_]=state.m_taskx[i][i_];
   state.m_pointindex=i;
//--- function call
   LSFitClearRequestFields(state);
//--- change values
   state.m_needf=true;
   state.m_rstate.stage=0;
//--- Saving state
   Func_lbl_rcomm(state,n,m,k,i,j,v,vv,relcnt);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for LSFitIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CLSFit::Func_lbl_14(CLSFitState &state,int &n,int &m,int &k,int &i,
                                int &j,double &v,double &vv,double &relcnt)
  {
//--- check
   if(!state.m_optstate.m_needfg)
     {
      //--- check
      if(!state.m_optstate.m_needfij)
         return(Func_lbl_24(state,n,m,k,i,j,v,vv,relcnt));
      //--- calculate Fi/jac(Fi)
      i=0;
      //--- function call, return result
      return(Func_lbl_26(state,n,m,k,i,j,v,vv,relcnt));
     }
//--- calculate F/gradF
   state.m_optstate.m_f=0;
   for(i=0;i<=k-1;i++)
      state.m_optstate.m_g[i]=0;
   i=0;
//--- function call, return result
   return(Func_lbl_21(state,n,m,k,i,j,v,vv,relcnt));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for LSFitIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CLSFit::Func_lbl_16(CLSFitState &state,int &n,int &m,int &k,int &i,
                                int &j,double &v,double &vv,double &relcnt)
  {
//--- check
   if(i>n-1)
      return(Func_lbl_7(state,n,m,k,i,j,v,vv,relcnt));
//--- copy
   for(int i_=0;i_<=k-1;i_++)
      state.m_c[i_]=state.m_optstate.m_x[i_];
   for(int i_=0;i_<=m-1;i_++)
      state.m_x[i_]=state.m_taskx[i][i_];
   state.m_pointindex=i;
//--- function call
   LSFitClearRequestFields(state);
//--- change values
   state.m_needf=true;
   state.m_rstate.stage=1;
//--- Saving state
   Func_lbl_rcomm(state,n,m,k,i,j,v,vv,relcnt);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for LSFitIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CLSFit::Func_lbl_21(CLSFitState &state,int &n,int &m,int &k,int &i,
                                int &j,double &v,double &vv,double &relcnt)
  {
//--- check
   if(i>n-1)
      return(Func_lbl_7(state,n,m,k,i,j,v,vv,relcnt));
//--- copy
   for(int i_=0;i_<=k-1;i_++)
      state.m_c[i_]=state.m_optstate.m_x[i_];
   for(int i_=0;i_<=m-1;i_++)
      state.m_x[i_]=state.m_taskx[i][i_];
   state.m_pointindex=i;
//--- function call
   LSFitClearRequestFields(state);
//--- change values
   state.m_needfg=true;
   state.m_rstate.stage=2;
//--- Saving state
   Func_lbl_rcomm(state,n,m,k,i,j,v,vv,relcnt);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for LSFitIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CLSFit::Func_lbl_24(CLSFitState &state,int &n,int &m,int &k,int &i,
                                int &j,double &v,double &vv,double &relcnt)
  {
//--- check
   if(!state.m_optstate.m_needfgh)
      return(Func_lbl_29(state,n,m,k,i,j,v,vv,relcnt));
//--- calculate F/grad(F)/hess(F)
   state.m_optstate.m_f=0;
   for(i=0;i<=k-1;i++)
      state.m_optstate.m_g[i]=0;
   for(i=0;i<=k-1;i++)
     {
      for(j=0;j<=k-1;j++)
         state.m_optstate.m_h[i].Set(j,0);
     }
//--- change value
   i=0;
//--- function call, return result
   return(Func_lbl_31(state,n,m,k,i,j,v,vv,relcnt));
  }
//+------------------------------------------------------------------+
//| Auxiliary function for LSFitIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CLSFit::Func_lbl_26(CLSFitState &state,int &n,int &m,int &k,int &i,
                                int &j,double &v,double &vv,double &relcnt)
  {
//--- check
   if(i>n-1)
      return(Func_lbl_7(state,n,m,k,i,j,v,vv,relcnt));
//--- copy
   for(int i_=0;i_<=k-1;i_++)
      state.m_c[i_]=state.m_optstate.m_x[i_];
   for(int i_=0;i_<=m-1;i_++)
      state.m_x[i_]=state.m_taskx[i][i_];
   state.m_pointindex=i;
//--- function call
   LSFitClearRequestFields(state);
//--- change values
   state.m_needfg=true;
   state.m_rstate.stage=3;
//--- Saving state
   Func_lbl_rcomm(state,n,m,k,i,j,v,vv,relcnt);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for LSFitIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CLSFit::Func_lbl_29(CLSFitState &state,int &n,int &m,int &k,int &i,
                                int &j,double &v,double &vv,double &relcnt)
  {
//--- check
   if(!state.m_optstate.m_xupdated)
      return(Func_lbl_7(state,n,m,k,i,j,v,vv,relcnt));
//--- Report new iteration
   for(int i_=0;i_<=k-1;i_++)
      state.m_c[i_]=state.m_optstate.m_x[i_];
   state.m_f=state.m_optstate.m_f;
//--- function call
   LSFitClearRequestFields(state);
//--- change values
   state.m_xupdated=true;
   state.m_rstate.stage=5;
//--- Saving state
   Func_lbl_rcomm(state,n,m,k,i,j,v,vv,relcnt);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for LSFitIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CLSFit::Func_lbl_31(CLSFitState &state,int &n,int &m,int &k,int &i,
                                int &j,double &v,double &vv,double &relcnt)
  {
//--- check
   if(i>n-1)
      return(Func_lbl_7(state,n,m,k,i,j,v,vv,relcnt));
//--- copy
   for(int i_=0;i_<=k-1;i_++)
      state.m_c[i_]=state.m_optstate.m_x[i_];
   for(int i_=0;i_<=m-1;i_++)
      state.m_x[i_]=state.m_taskx[i][i_];
   state.m_pointindex=i;
//--- function call
   LSFitClearRequestFields(state);
//--- change values
   state.m_needfgh=true;
   state.m_rstate.stage=4;
//--- Saving state
   Func_lbl_rcomm(state,n,m,k,i,j,v,vv,relcnt);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for LSFitIteration. Is a product to get rid of|
//| the operator unconditional jump goto.                            |
//+------------------------------------------------------------------+
static bool CLSFit::Func_lbl_38(CLSFitState &state,int &n,int &m,int &k,int &i,
                                int &j,double &v,double &vv,double &relcnt)
  {
//--- check
   if(i>n-1)
     {
      //--- change values
      state.m_reprmserror=MathSqrt(state.m_reprmserror/n);
      state.m_repwrmserror=MathSqrt(state.m_repwrmserror/n);
      state.m_repavgerror=state.m_repavgerror/n;
      //--- check
      if(relcnt!=0.0)
         state.m_repavgrelerror=state.m_repavgrelerror/relcnt;
      //--- return result
      return(false);
     }
//--- copy
   for(int i_=0;i_<=k-1;i_++)
      state.m_c[i_]=state.m_c[i_];
   for(int i_=0;i_<=m-1;i_++)
      state.m_x[i_]=state.m_taskx[i][i_];
   state.m_pointindex=i;
//--- function call
   LSFitClearRequestFields(state);
//--- change values
   state.m_needf=true;
   state.m_rstate.stage=6;
//--- Saving state
   Func_lbl_rcomm(state,n,m,k,i,j,v,vv,relcnt);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Parametric spline inteprolant: 2-dimensional curve.              |
//| You should not try to access its members directly - use          |
//| PSpline2XXXXXXXX() functions instead.                            |
//+------------------------------------------------------------------+
class CPSpline2Interpolant
  {
public:
   //--- variables
   int               m_n;
   bool              m_periodic;
   CSpline1DInterpolant m_x;
   CSpline1DInterpolant m_y;
   //--- array
   double            m_p[];
   //--- constructor, destructor
                     CPSpline2Interpolant(void);
                    ~CPSpline2Interpolant(void);
   //--- copy
   void              Copy(CPSpline2Interpolant &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CPSpline2Interpolant::CPSpline2Interpolant(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CPSpline2Interpolant::~CPSpline2Interpolant(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CPSpline2Interpolant::Copy(CPSpline2Interpolant &obj)
  {
//--- copy variables
   m_n=obj.m_n;
   m_periodic=obj.m_periodic;
   m_x.Copy(obj.m_x);
   m_y.Copy(obj.m_y);
//--- copy array
   ArrayCopy(m_p,obj.m_p);
  }
//+------------------------------------------------------------------+
//| Parametric spline inteprolant: 2-dimensional curve.              |
//| You should not try to access its members directly - use          |
//| PSpline2XXXXXXXX() functions instead.                            |
//+------------------------------------------------------------------+
class CPSpline2InterpolantShell
  {
private:
   CPSpline2Interpolant m_innerobj;
public:
   //--- constructors, destructor
                     CPSpline2InterpolantShell(void);
                     CPSpline2InterpolantShell(CPSpline2Interpolant &obj);
                    ~CPSpline2InterpolantShell(void);
   //--- method
   CPSpline2Interpolant *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CPSpline2InterpolantShell::CPSpline2InterpolantShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CPSpline2InterpolantShell::CPSpline2InterpolantShell(CPSpline2Interpolant &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CPSpline2InterpolantShell::~CPSpline2InterpolantShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CPSpline2Interpolant *CPSpline2InterpolantShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Parametric spline inteprolant: 3-dimensional curve.              |
//| You should not try to access its members directly - use          |
//| PSpline3XXXXXXXX() functions instead.                            |
//+------------------------------------------------------------------+
class CPSpline3Interpolant
  {
public:
   //--- variables
   int               m_n;
   bool              m_periodic;
   CSpline1DInterpolant m_x;
   CSpline1DInterpolant m_y;
   CSpline1DInterpolant m_z;
   //--- array
   double            m_p[];
   //--- constructor, destructor
                     CPSpline3Interpolant(void);
                    ~CPSpline3Interpolant(void);
   //--- copy
   void              Copy(CPSpline3Interpolant &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CPSpline3Interpolant::CPSpline3Interpolant(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CPSpline3Interpolant::~CPSpline3Interpolant(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CPSpline3Interpolant::Copy(CPSpline3Interpolant &obj)
  {
//--- copy variables
   m_n=obj.m_n;
   m_periodic=obj.m_periodic;
   m_x.Copy(obj.m_x);
   m_y.Copy(obj.m_y);
   m_z.Copy(obj.m_z);
//--- copy array
   ArrayCopy(m_p,obj.m_p);
  }
//+------------------------------------------------------------------+
//| Parametric spline inteprolant: 3-dimensional curve.              |
//| You should not try to access its members directly - use          |
//| PSpline3XXXXXXXX() functions instead.                            |
//+------------------------------------------------------------------+
class CPSpline3InterpolantShell
  {
private:
   CPSpline3Interpolant m_innerobj;
public:
   //--- constructors, destructor
                     CPSpline3InterpolantShell(void);
                     CPSpline3InterpolantShell(CPSpline3Interpolant &obj);
                    ~CPSpline3InterpolantShell(void);
   //--- method
   CPSpline3Interpolant *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CPSpline3InterpolantShell::CPSpline3InterpolantShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CPSpline3InterpolantShell::CPSpline3InterpolantShell(CPSpline3Interpolant &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CPSpline3InterpolantShell::~CPSpline3InterpolantShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CPSpline3Interpolant *CPSpline3InterpolantShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Parametric spline                                                |
//+------------------------------------------------------------------+
class CPSpline
  {
private:
   //--- private methods
   static void       PSpline2Par(CMatrixDouble &xy,const int n,const int pt,double &p[]);
   static void       PSpline3Par(CMatrixDouble &xy,const int n,const int pt,double &p[]);
public:
   //--- constructor, destructor
                     CPSpline(void);
                    ~CPSpline(void);
   //--- public methods
   static void       PSpline2Build(CMatrixDouble &cxy,const int n,const int st,const int pt,CPSpline2Interpolant &p);
   static void       PSpline3Build(CMatrixDouble &cxy,const int n,const int st,const int pt,CPSpline3Interpolant &p);
   static void       PSpline2BuildPeriodic(CMatrixDouble &cxy,const int n,const int st,const int pt,CPSpline2Interpolant &p);
   static void       PSpline3BuildPeriodic(CMatrixDouble &cxy,const int n,const int st,const int pt,CPSpline3Interpolant &p);
   static void       PSpline2ParameterValues(CPSpline2Interpolant &p,int &n,double &t[]);
   static void       PSpline3ParameterValues(CPSpline3Interpolant &p,int &n,double &t[]);
   static void       PSpline2Calc(CPSpline2Interpolant &p,double t,double &x,double &y);
   static void       PSpline3Calc(CPSpline3Interpolant &p,double t,double &x,double &y,double &z);
   static void       PSpline2Tangent(CPSpline2Interpolant &p,double t,double &x,double &y);
   static void       PSpline3Tangent(CPSpline3Interpolant &p,double t,double &x,double &y,double &z);
   static void       PSpline2Diff(CPSpline2Interpolant &p,double t,double &x,double &dx,double &y,double &dy);
   static void       PSpline3Diff(CPSpline3Interpolant &p,double t,double &x,double &dx,double &y,double &dy,double &z,double &dz);
   static void       PSpline2Diff2(CPSpline2Interpolant &p,double t,double &x,double &dx,double &d2x,double &y,double &dy,double &d2y);
   static void       PSpline3Diff2(CPSpline3Interpolant &p,double t,double &x,double &dx,double &d2x,double &y,double &dy,double &d2y,double &z,double &dz,double &d2z);
   static double     PSpline2ArcLength(CPSpline2Interpolant &p,const double a,const double b);
   static double     PSpline3ArcLength(CPSpline3Interpolant &p,const double a,const double b);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CPSpline::CPSpline(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CPSpline::~CPSpline(void)
  {

  }
//+------------------------------------------------------------------+
//| This function builds non-periodic 2-dimensional parametric       |
//| spline which starts at (X[0],Y[0]) and ends at (X[N-1],Y[N-1]).  |
//| INPUT PARAMETERS:                                                |
//|     XY  -   points, array[0..N-1,0..1].                          |
//|             XY[I,0:1] corresponds to the Ith point.              |
//|             Order of points is important!                        |
//|     N   -   points count, N>=5 for Akima splines, N>=2 for other |
//|             types of splines.                                    |
//|     ST  -   spline type:                                         |
//|             * 0     Akima spline                                 |
//|             * 1     parabolically terminated Catmull-Rom spline  |
//|                     (Tension=0)                                  |
//|             * 2     parabolically terminated cubic spline        |
//|     PT  -   parameterization type:                               |
//|             * 0     uniform                                      |
//|             * 1     chord length                                 |
//|             * 2     centripetal                                  |
//| OUTPUT PARAMETERS:                                               |
//|     P   -   parametric spline interpolant                        |
//| NOTES:                                                           |
//| * this function assumes that there all consequent points are     |
//|   distinct. I.e. (x0,y0)<>(x1,y1), (x1,y1)<>(x2,y2),             |
//|   (x2,y2)<>(x3,y3) and so on. However, non-consequent points may |
//|   coincide, i.e. we can have (x0,y0) = (x2,y2).                  |
//+------------------------------------------------------------------+
static void CPSpline::PSpline2Build(CMatrixDouble &cxy,const int n,
                                    const int st,const int pt,
                                    CPSpline2Interpolant &p)
  {
//--- create a variable
   int i_=0;
//--- create array
   double tmp[];
//--- copy matrix
   CMatrixDouble xy;
   xy=cxy;
//--- check
   if(!CAp::Assert(st>=0 && st<=2,__FUNCTION__+": incorrect spline type!"))
      return;
//--- check
   if(!CAp::Assert(pt>=0 && pt<=2,__FUNCTION__+": incorrect parameterization type!"))
      return;
//--- check
   if(st==0)
     {
      //--- check
      if(!CAp::Assert(n>=5,__FUNCTION__+": N<5 (minimum value for Akima splines)!"))
         return;
     }
   else
     {
      //--- check
      if(!CAp::Assert(n>=2,__FUNCTION__+": N<2!"))
         return;
     }
//--- Prepare
   p.m_n=n;
   p.m_periodic=false;
//--- allocation
   ArrayResizeAL(tmp,n);
//--- Build parameterization,check that all parameters are distinct
   PSpline2Par(xy,n,pt,p.m_p);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(p.m_p,n),__FUNCTION__+": consequent points are too close!"))
      return;
//--- Build splines
   if(st==0)
     {
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][0];
      //--- function call
      CSpline1D::Spline1DBuildAkima(p.m_p,tmp,n,p.m_x);
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][1];
      //--- function call
      CSpline1D::Spline1DBuildAkima(p.m_p,tmp,n,p.m_y);
     }
//--- check
   if(st==1)
     {
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][0];
      //--- function call
      CSpline1D::Spline1DBuildCatmullRom(p.m_p,tmp,n,0,0.0,p.m_x);
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][1];
      //--- function call
      CSpline1D::Spline1DBuildCatmullRom(p.m_p,tmp,n,0,0.0,p.m_y);
     }
//--- check
   if(st==2)
     {
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][0];
      //--- function call
      CSpline1D::Spline1DBuildCubic(p.m_p,tmp,n,0,0.0,0,0.0,p.m_x);
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][1];
      //--- function call
      CSpline1D::Spline1DBuildCubic(p.m_p,tmp,n,0,0.0,0,0.0,p.m_y);
     }
  }
//+------------------------------------------------------------------+
//| This function builds non-periodic 3-dimensional parametric spline|
//| which starts at (X[0],Y[0],Z[0]) and ends at                     |
//| (X[N-1],Y[N-1],Z[N-1]).                                          |
//| Same as PSpline2Build() function, but for 3D, so we won't        |
//| duplicate its description here.                                  |
//+------------------------------------------------------------------+
static void CPSpline::PSpline3Build(CMatrixDouble &cxy,const int n,
                                    const int st,const int pt,
                                    CPSpline3Interpolant &p)
  {
//--- create a variable
   int i_=0;
//--- create array
   double tmp[];
//--- copy matrix
   CMatrixDouble xy;
   xy=cxy;
//--- check
   if(!CAp::Assert(st>=0 && st<=2,__FUNCTION__+": incorrect spline type!"))
      return;
//--- check
   if(!CAp::Assert(pt>=0 && pt<=2,__FUNCTION__+": incorrect parameterization type!"))
      return;
//--- check
   if(st==0)
     {
      //--- check
      if(!CAp::Assert(n>=5,__FUNCTION__+": N<5 (minimum value for Akima splines)!"))
         return;
     }
   else
     {
      //--- check
      if(!CAp::Assert(n>=2,__FUNCTION__+"PSpline3Build: N<2!"))
         return;
     }
//--- Prepare
   p.m_n=n;
   p.m_periodic=false;
//--- allocation
   ArrayResizeAL(tmp,n);
//--- Build parameterization,check that all parameters are distinct
   PSpline3Par(xy,n,pt,p.m_p);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(p.m_p,n),__FUNCTION__+": consequent points are too close!"))
      return;
//--- Build splines
   if(st==0)
     {
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][0];
      //--- function call
      CSpline1D::Spline1DBuildAkima(p.m_p,tmp,n,p.m_x);
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][1];
      //--- function call
      CSpline1D::Spline1DBuildAkima(p.m_p,tmp,n,p.m_y);
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][2];
      //--- function call
      CSpline1D::Spline1DBuildAkima(p.m_p,tmp,n,p.m_z);
     }
//--- check
   if(st==1)
     {
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][0];
      //--- function call
      CSpline1D::Spline1DBuildCatmullRom(p.m_p,tmp,n,0,0.0,p.m_x);
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][1];
      //--- function call
      CSpline1D::Spline1DBuildCatmullRom(p.m_p,tmp,n,0,0.0,p.m_y);
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][2];
      //--- function call
      CSpline1D::Spline1DBuildCatmullRom(p.m_p,tmp,n,0,0.0,p.m_z);
     }
//--- check
   if(st==2)
     {
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][0];
      //--- function call
      CSpline1D::Spline1DBuildCubic(p.m_p,tmp,n,0,0.0,0,0.0,p.m_x);
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][1];
      //--- function call
      CSpline1D::Spline1DBuildCubic(p.m_p,tmp,n,0,0.0,0,0.0,p.m_y);
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         tmp[i_]=xy[i_][2];
      //--- function call
      CSpline1D::Spline1DBuildCubic(p.m_p,tmp,n,0,0.0,0,0.0,p.m_z);
     }
  }
//+------------------------------------------------------------------+
//| This function builds periodic 2-dimensional parametric spline    |
//| which starts at (X[0],Y[0]), goes through all points to          |
//| (X[N-1],Y[N-1]) and then back to (X[0],Y[0]).                    |
//| INPUT PARAMETERS:                                                |
//|     XY  -   points, array[0..N-1,0..1].                          |
//|             XY[I,0:1] corresponds to the Ith point.              |
//|             XY[N-1,0:1] must be different from XY[0,0:1].        |
//|             Order of points is important!                        |
//|     N   -   points count, N>=3 for other types of splines.       |
//|     ST  -   spline type:                                         |
//|             * 1     Catmull-Rom spline (Tension=0) with cyclic   |
//|                     boundary conditions                          |
//|             * 2     cubic spline with cyclic boundary conditions |
//|     PT  -   parameterization type:                               |
//|             * 0     uniform                                      |
//|             * 1     chord length                                 |
//|             * 2     centripetal                                  |
//| OUTPUT PARAMETERS:                                               |
//|     P   -   parametric spline interpolant                        |
//| NOTES:                                                           |
//| * this function assumes that there all consequent points are     |
//|   distinct. I.e. (x0,y0)<>(x1,y1), (x1,y1)<>(x2,y2),             |
//|   (x2,y2)<>(x3,y3) and so on. However, non-consequent points may |
//|   coincide, i.e. we can  have (x0,y0) = (x2,y2).                 |
//| * last point of sequence is NOT equal to the first point. You    |
//|   shouldn't make curve "explicitly periodic" by making them      |
//|   equal.                                                         |
//+------------------------------------------------------------------+
static void CPSpline::PSpline2BuildPeriodic(CMatrixDouble &cxy,const int n,
                                            const int st,const int pt,
                                            CPSpline2Interpolant &p)
  {
//--- create a variable
   int i_=0;
//--- create array
   double tmp[];
//--- create matrix
   CMatrixDouble xyp;
   CMatrixDouble xy;
//--- copy matrix
   xy=cxy;
//--- check
   if(!CAp::Assert(st>=1 && st<=2,__FUNCTION__+": incorrect spline type!"))
      return;
//--- check
   if(!CAp::Assert(pt>=0 && pt<=2,__FUNCTION__+": incorrect parameterization type!"))
      return;
//--- check
   if(!CAp::Assert(n>=3,__FUNCTION__+": N<3!"))
      return;
//--- Prepare
   p.m_n=n;
   p.m_periodic=true;
//--- allocation
   ArrayResizeAL(tmp,n+1);
   xyp.Resize(n+1,2);
//--- change values
   for(i_=0;i_<=n-1;i_++)
      xyp[i_].Set(0,xy[i_][0]);
   for(i_=0;i_<=n-1;i_++)
      xyp[i_].Set(1,xy[i_][1]);
   for(i_=0;i_<=1;i_++)
      xyp[n].Set(i_,xy[0][i_]);
//--- Build parameterization,check that all parameters are distinct
   PSpline2Par(xyp,n+1,pt,p.m_p);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(p.m_p,n+1),__FUNCTION__+": consequent (or first and last) points are too close!"))
      return;
//--- Build splines
   if(st==1)
     {
      //--- copy
      for(i_=0;i_<=n;i_++)
         tmp[i_]=xyp[i_][0];
      //--- function call
      CSpline1D::Spline1DBuildCatmullRom(p.m_p,tmp,n+1,-1,0.0,p.m_x);
      //--- copy
      for(i_=0;i_<=n;i_++)
         tmp[i_]=xyp[i_][1];
      //--- function call
      CSpline1D::Spline1DBuildCatmullRom(p.m_p,tmp,n+1,-1,0.0,p.m_y);
     }
//--- check
   if(st==2)
     {
      //--- copy
      for(i_=0;i_<=n;i_++)
         tmp[i_]=xyp[i_][0];
      //--- function call
      CSpline1D::Spline1DBuildCubic(p.m_p,tmp,n+1,-1,0.0,-1,0.0,p.m_x);
      //--- copy
      for(i_=0;i_<=n;i_++)
         tmp[i_]=xyp[i_][1];
      //--- function call
      CSpline1D::Spline1DBuildCubic(p.m_p,tmp,n+1,-1,0.0,-1,0.0,p.m_y);
     }
  }
//+------------------------------------------------------------------+
//| This function builds periodic 3-dimensional parametric spline    |
//| which starts at (X[0],Y[0],Z[0]), goes through all points to     |
//| (X[N-1],Y[N-1],Z[N-1]) and then back to (X[0],Y[0],Z[0]).        |
//| Same as PSpline2Build() function, but for 3D, so we won't        |
//| duplicate its description here.                                  |
//+------------------------------------------------------------------+
static void CPSpline::PSpline3BuildPeriodic(CMatrixDouble &cxy,const int n,
                                            const int st,const int pt,
                                            CPSpline3Interpolant &p)
  {
//--- create a variable
   int i_=0;
//--- create array
   double tmp[];
//--- create matrix
   CMatrixDouble xyp;
   CMatrixDouble xy;
//--- copy matrix
   xy=cxy;
//--- check
   if(!CAp::Assert(st>=1 && st<=2,__FUNCTION__+": incorrect spline type!"))
      return;
//--- check
   if(!CAp::Assert(pt>=0 && pt<=2,__FUNCTION__+": incorrect parameterization type!"))
      return;
//--- check
   if(!CAp::Assert(n>=3,__FUNCTION__+": N<3!"))
      return;
//--- Prepare
   p.m_n=n;
   p.m_periodic=true;
//--- allocation
   ArrayResizeAL(tmp,n+1);
   xyp.Resize(n+1,3);
//--- change values
   for(i_=0;i_<=n-1;i_++)
      xyp[i_].Set(0,xy[i_][0]);
   for(i_=0;i_<=n-1;i_++)
      xyp[i_].Set(1,xy[i_][1]);
   for(i_=0;i_<=n-1;i_++)
      xyp[i_].Set(2,xy[i_][2]);
   for(i_=0;i_<=2;i_++)
      xyp[n].Set(i_,xy[0][i_]);
//--- Build parameterization,check that all parameters are distinct
   PSpline3Par(xyp,n+1,pt,p.m_p);
//--- check
   if(!CAp::Assert(CApServ::AreDistinct(p.m_p,n+1),__FUNCTION__+": consequent (or first and last) points are too close!"))
      return;
//--- Build splines
   if(st==1)
     {
      //--- copy
      for(i_=0;i_<=n;i_++)
         tmp[i_]=xyp[i_][0];
      //--- function call
      CSpline1D::Spline1DBuildCatmullRom(p.m_p,tmp,n+1,-1,0.0,p.m_x);
      //--- copy
      for(i_=0;i_<=n;i_++)
         tmp[i_]=xyp[i_][1];
      //--- function call
      CSpline1D::Spline1DBuildCatmullRom(p.m_p,tmp,n+1,-1,0.0,p.m_y);
      //--- copy
      for(i_=0;i_<=n;i_++)
         tmp[i_]=xyp[i_][2];
      //--- function call
      CSpline1D::Spline1DBuildCatmullRom(p.m_p,tmp,n+1,-1,0.0,p.m_z);
     }
//--- check
   if(st==2)
     {
      //--- copy
      for(i_=0;i_<=n;i_++)
         tmp[i_]=xyp[i_][0];
      //--- function call
      CSpline1D::Spline1DBuildCubic(p.m_p,tmp,n+1,-1,0.0,-1,0.0,p.m_x);
      //--- copy
      for(i_=0;i_<=n;i_++)
         tmp[i_]=xyp[i_][1];
      //--- function call
      CSpline1D::Spline1DBuildCubic(p.m_p,tmp,n+1,-1,0.0,-1,0.0,p.m_y);
      //--- copy
      for(i_=0;i_<=n;i_++)
         tmp[i_]=xyp[i_][2];
      //--- function call
      CSpline1D::Spline1DBuildCubic(p.m_p,tmp,n+1,-1,0.0,-1,0.0,p.m_z);
     }
  }
//+------------------------------------------------------------------+
//| This function returns vector of parameter values correspoding to |
//| points.                                                          |
//| I.e. for P created from (X[0],Y[0])...(X[N-1],Y[N-1]) and        |
//| U=TValues(P) we have                                             |
//|     (X[0],Y[0]) = PSpline2Calc(P,U[0]),                          |
//|     (X[1],Y[1]) = PSpline2Calc(P,U[1]),                          |
//|     (X[2],Y[2]) = PSpline2Calc(P,U[2]),                          |
//|     ...                                                          |
//| INPUT PARAMETERS:                                                |
//|     P   -   parametric spline interpolant                        |
//| OUTPUT PARAMETERS:                                               |
//|     N   -   array size                                           |
//|     T   -   array[0..N-1]                                        |
//| NOTES:                                                           |
//| * for non-periodic splines U[0]=0, U[0]<U[1]<...<U[N-1], U[N-1]=1|
//| * for periodic splines     U[0]=0, U[0]<U[1]<...<U[N-1], U[N-1]<1|
//+------------------------------------------------------------------+
static void CPSpline::PSpline2ParameterValues(CPSpline2Interpolant &p,
                                              int &n,double &t[])
  {
//--- create a variable
   int i_=0;
//--- initialization
   n=0;
//--- check
   if(!CAp::Assert(p.m_n>=2,__FUNCTION__+": internal error!"))
      return;
//--- initialization
   n=p.m_n;
//--- allocation
   ArrayResizeAL(t,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      t[i_]=p.m_p[i_];
   t[0]=0;
//--- check
   if(!p.m_periodic)
      t[n-1]=1;
  }
//+------------------------------------------------------------------+
//| This function returns vector of parameter values correspoding to |
//| points.                                                          |
//| Same as PSpline2ParameterValues(), but for 3D.                   |
//+------------------------------------------------------------------+
static void CPSpline::PSpline3ParameterValues(CPSpline3Interpolant &p,
                                              int &n,double &t[])
  {
//--- create a variable
   int i_=0;
//--- initialization
   n=0;
//--- check
   if(!CAp::Assert(p.m_n>=2,__FUNCTION__+": internal error!"))
      return;
//--- initialization
   n=p.m_n;
//--- allocation
   ArrayResizeAL(t,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      t[i_]=p.m_p[i_];
   t[0]=0;
//--- check
   if(!p.m_periodic)
      t[n-1]=1;
  }
//+------------------------------------------------------------------+
//| This function calculates the value of the parametric spline for a| 
//| given value of parameter T                                       |
//| INPUT PARAMETERS:                                                |
//|     P   -   parametric spline interpolant                        |
//|     T   -   point:                                               |
//|             * T in [0,1] corresponds to interval spanned by      |
//|               points                                             |
//|             * for non-periodic splines T<0 (or T>1) correspond to|
//|               parts of the curve before the first (after the     |
//|               last) point                                        |
//|             * for periodic splines T<0 (or T>1) are projected    |
//|               into [0,1] by making T=T-floor(T).                 |
//| OUTPUT PARAMETERS:                                               |
//|     X   -   X-position                                           |
//|     Y   -   Y-position                                           |
//+------------------------------------------------------------------+
static void CPSpline::PSpline2Calc(CPSpline2Interpolant &p,double t,
                                   double &x,double &y)
  {
//--- initialization
   x=0;
   y=0;
//--- check
   if(p.m_periodic)
      t=t-(int)MathFloor(t);
//--- function call
   x=CSpline1D::Spline1DCalc(p.m_x,t);
//--- function call
   y=CSpline1D::Spline1DCalc(p.m_y,t);
  }
//+------------------------------------------------------------------+
//| This function calculates the value of the parametric spline for a|
//| given value of parameter T.                                      |
//| INPUT PARAMETERS:                                                |
//|     P   -   parametric spline interpolant                        |
//|     T   -   point:                                               |
//|             * T in [0,1] corresponds to interval spanned by      |
//|               points                                             |
//|             * for non-periodic splines T<0 (or T>1) correspond   |
//|               to parts of the curve before the first (after the  |
//|               last) point                                        |
//|             * for periodic splines T<0 (or T>1) are projected    |
//|               into [0,1] by making T=T-floor(T).                 |
//| OUTPUT PARAMETERS:                                               |
//|     X   -   X-position                                           |
//|     Y   -   Y-position                                           |
//|     Z   -   Z-position                                           |
//+------------------------------------------------------------------+
static void CPSpline::PSpline3Calc(CPSpline3Interpolant &p,double t,
                                   double &x,double &y,double &z)
  {
//--- initialization
   x=0;
   y=0;
   z=0;
//--- check
   if(p.m_periodic)
      t=t-(int)MathFloor(t);
//--- function call
   x=CSpline1D::Spline1DCalc(p.m_x,t);
//--- function call
   y=CSpline1D::Spline1DCalc(p.m_y,t);
//--- function call
   z=CSpline1D::Spline1DCalc(p.m_z,t);
  }
//+------------------------------------------------------------------+
//| This function calculates tangent vector for a given value of     |
//| parameter T                                                      |
//| INPUT PARAMETERS:                                                |
//|     P   -   parametric spline interpolant                        |
//|     T   -   point:                                               |
//|             * T in [0,1] corresponds to interval spanned by      |
//|               points                                             |
//|             * for non-periodic splines T<0 (or T>1) correspond to|
//|               parts of the curve before the first (after the     |
//|               last) point                                        |
//|             * for periodic splines T<0 (or T>1) are projected    |
//|               into [0,1] by making T=T-floor(T).                 |
//| OUTPUT PARAMETERS:                                               |
//|     X    -   X-component of tangent vector (normalized)          |
//|     Y    -   Y-component of tangent vector (normalized)          |
//| NOTE:                                                            |
//|     X^2+Y^2 is either 1 (for non-zero tangent vector) or 0.      |
//+------------------------------------------------------------------+
static void CPSpline::PSpline2Tangent(CPSpline2Interpolant &p,double t,
                                      double &x,double &y)
  {
//--- create variables
   double v=0;
   double v0=0;
   double v1=0;
//--- initialization
   x=0;
   y=0;
//--- check
   if(p.m_periodic)
      t=t-(int)MathFloor(t);
//--- function call
   PSpline2Diff(p,t,v0,x,v1,y);
//--- check
   if(x!=0.0 || y!=0.0)
     {
      //--- this code is a bit more complex than X^2+Y^2 to avoid
      //--- overflow for large values of X and Y.
      v=CApServ::SafePythag2(x,y);
      x=x/v;
      y=y/v;
     }
  }
//+------------------------------------------------------------------+
//| This function calculates tangent vector for a given value of     |
//| parameter T                                                      |
//| INPUT PARAMETERS:                                                |
//|     P   -   parametric spline interpolant                        |
//|     T   -   point:                                               |
//|             * T in [0,1] corresponds to interval spanned by      |
//|               points                                             |
//|             * for non-periodic splines T<0 (or T>1) correspond to|
//|               parts of the curve before the first (after the     |
//|               last) point                                        |
//|             * for periodic splines T<0 (or T>1) are projected    |
//|               into [0,1] by making T=T-floor(T).                 |
//| OUTPUT PARAMETERS:                                               |
//|     X    -   X-component of tangent vector (normalized)          |
//|     Y    -   Y-component of tangent vector (normalized)          |
//|     Z    -   Z-component of tangent vector (normalized)          |
//| NOTE:                                                            |
//|     X^2+Y^2+Z^2 is either 1 (for non-zero tangent vector) or 0.  |
//+------------------------------------------------------------------+
static void CPSpline::PSpline3Tangent(CPSpline3Interpolant &p,double t,
                                      double &x,double &y,double &z)
  {
//--- create variables
   double v=0;
   double v0=0;
   double v1=0;
   double v2=0;
//--- initialization
   x=0;
   y=0;
   z=0;
//--- check
   if(p.m_periodic)
      t=t-(int)MathFloor(t);
//--- function call
   PSpline3Diff(p,t,v0,x,v1,y,v2,z);
//--- check
   if((x!=0.0 || y!=0.0) || z!=0.0)
     {
      //--- function call
      v=CApServ::SafePythag3(x,y,z);
      //--- change values
      x=x/v;
      y=y/v;
      z=z/v;
     }
  }
//+------------------------------------------------------------------+
//| This function calculates derivative, i.e. it returns             |
//| (dX/dT,dY/dT).                                                   |
//| INPUT PARAMETERS:                                                |
//|     P   -   parametric spline interpolant                        |
//|     T   -   point:                                               |
//|             * T in [0,1] corresponds to interval spanned by      |
//|               points                                             |
//|             * for non-periodic splines T<0 (or T>1) correspond to|
//|               parts of the curve before the first (after the     |
//|               last) point                                        |
//|             * for periodic splines T<0 (or T>1) are projected    |
//|               into [0,1] by making T=T-floor(T).                 |
//| OUTPUT PARAMETERS:                                               |
//|     X   -   X-value                                              |
//|     DX  -   X-derivative                                         |
//|     Y   -   Y-value                                              |
//|     DY  -   Y-derivative                                         |
//+------------------------------------------------------------------+
static void CPSpline::PSpline2Diff(CPSpline2Interpolant &p,double t,
                                   double &x,double &dx,double &y,
                                   double &dy)
  {
//--- create a variable
   double d2s=0;
//--- change values
   x=0;
   dx=0;
   y=0;
   dy=0;
//--- check
   if(p.m_periodic)
      t=t-(int)MathFloor(t);
//--- function call
   CSpline1D::Spline1DDiff(p.m_x,t,x,dx,d2s);
//--- function call
   CSpline1D::Spline1DDiff(p.m_y,t,y,dy,d2s);
  }
//+------------------------------------------------------------------+
//| This function calculates derivative, i.e. it returns             |
//| (dX/dT,dY/dT,dZ/dT).                                             |
//| INPUT PARAMETERS:                                                |
//|     P   -   parametric spline interpolant                        |
//|     T   -   point:                                               |
//|             * T in [0,1] corresponds to interval spanned by      |
//|               points                                             |
//|             * for non-periodic splines T<0 (or T>1) correspond to|
//|               parts of the curve before the first (after the     |
//|               last) point                                        |
//|             * for periodic splines T<0 (or T>1) are projected    |
//|               into [0,1] by making T=T-floor(T).                 |
//| OUTPUT PARAMETERS:                                               |
//|     X   -   X-value                                              |
//|     DX  -   X-derivative                                         |
//|     Y   -   Y-value                                              |
//|     DY  -   Y-derivative                                         |
//|     Z   -   Z-value                                              |
//|     DZ  -   Z-derivative                                         |
//+------------------------------------------------------------------+
static void CPSpline::PSpline3Diff(CPSpline3Interpolant &p,double t,
                                   double &x,double &dx,double &y,
                                   double &dy,double &z,double &dz)
  {
//--- create a variable
   double d2s=0;
//--- initialization
   x=0;
   dx=0;
   y=0;
   dy=0;
   z=0;
   dz=0;
//--- check
   if(p.m_periodic)
      t=t-(int)MathFloor(t);
//--- function call
   CSpline1D::Spline1DDiff(p.m_x,t,x,dx,d2s);
//--- function call
   CSpline1D::Spline1DDiff(p.m_y,t,y,dy,d2s);
//--- function call
   CSpline1D::Spline1DDiff(p.m_z,t,z,dz,d2s);
  }
//+------------------------------------------------------------------+
//| This function calculates first and second derivative with respect|
//| to T.                                                            |
//| INPUT PARAMETERS:                                                |
//|     P   -   parametric spline interpolant                        |
//|     T   -   point:                                               |
//|             * T in [0,1] corresponds to interval spanned by      |
//|               points                                             |
//|             * for non-periodic splines T<0 (or T>1) correspond to|
//|               parts of the curve before the first (after the     |
//|               last) point                                        |
//|             * for periodic splines T<0 (or T>1) are projected    |
//|               into [0,1] by making T=T-floor(T).                 |
//| OUTPUT PARAMETERS:                                               |
//|     X   -   X-value                                              |
//|     DX  -   derivative                                           |
//|     D2X -   second derivative                                    |
//|     Y   -   Y-value                                              |
//|     DY  -   derivative                                           |
//|     D2Y -   second derivative                                    |
//+------------------------------------------------------------------+
static void CPSpline::PSpline2Diff2(CPSpline2Interpolant &p,double t,
                                    double &x,double &dx,double &d2x,
                                    double &y,double &dy,double &d2y)
  {
//--- initialization
   x=0;
   dx=0;
   d2x=0;
   y=0;
   dy=0;
   d2y=0;
//--- check
   if(p.m_periodic)
      t=t-(int)MathFloor(t);
//--- function call
   CSpline1D::Spline1DDiff(p.m_x,t,x,dx,d2x);
//--- function call
   CSpline1D::Spline1DDiff(p.m_y,t,y,dy,d2y);
  }
//+------------------------------------------------------------------+
//| This function calculates first and second derivative with respect|
//| to T.                                                            |
//| INPUT PARAMETERS:                                                |
//|     P   -   parametric spline interpolant                        |
//|     T   -   point:                                               |
//|             * T in [0,1] corresponds to interval spanned by      |
//|               points                                             |
//|             * for non-periodic splines T<0 (or T>1) correspond to|
//|               parts of the curve before the first (after the     |
//|               last) point                                        |
//|             * for periodic splines T<0 (or T>1) are projected    |
//|               into [0,1] by making T=T-floor(T).                 |
//| OUTPUT PARAMETERS:                                               |
//|     X   -   X-value                                              |
//|     DX  -   derivative                                           |
//|     D2X -   second derivative                                    |
//|     Y   -   Y-value                                              |
//|     DY  -   derivative                                           |
//|     D2Y -   second derivative                                    |
//|     Z   -   Z-value                                              |
//|     DZ  -   derivative                                           |
//|     D2Z -   second derivative                                    |
//+------------------------------------------------------------------+
static void CPSpline::PSpline3Diff2(CPSpline3Interpolant &p,double t,
                                    double &x,double &dx,double &d2x,
                                    double &y,double &dy,double &d2y,
                                    double &z,double &dz,double &d2z)
  {
//--- initialization
   x=0;
   dx=0;
   d2x=0;
   y=0;
   dy=0;
   d2y=0;
   z=0;
   dz=0;
   d2z=0;
//--- check
   if(p.m_periodic)
      t=t-(int)MathFloor(t);
//--- function call
   CSpline1D::Spline1DDiff(p.m_x,t,x,dx,d2x);
//--- function call
   CSpline1D::Spline1DDiff(p.m_y,t,y,dy,d2y);
//--- function call
   CSpline1D::Spline1DDiff(p.m_z,t,z,dz,d2z);
  }
//+------------------------------------------------------------------+
//| This function calculates arc length, i.e. length of curve between|
//| t=a and t=b.                                                     |
//| INPUT PARAMETERS:                                                |
//|     P   -   parametric spline interpolant                        |
//|     A,B -   parameter values corresponding to arc ends:          |
//|             * B>A will result in positive length returned        |
//|             * B<A will result in negative length returned        |
//| RESULT:                                                          |
//|     length of arc starting at T=A and ending at T=B.             |
//+------------------------------------------------------------------+
static double CPSpline::PSpline2ArcLength(CPSpline2Interpolant &p,const double a,
                                          const double b)
  {
//--- create variables
   double result=0;
   double sx=0;
   double dsx=0;
   double d2sx=0;
   double sy=0;
   double dsy=0;
   double d2sy=0;
//--- objects of classes
   CAutoGKState  state;
   CAutoGKReport rep;
//--- function call
   CAutoGK::AutoGKSmooth(a,b,state);
//--- cycle
   while(CAutoGK::AutoGKIteration(state))
     {
      CSpline1D::Spline1DDiff(p.m_x,state.m_x,sx,dsx,d2sx);
      CSpline1D::Spline1DDiff(p.m_y,state.m_x,sy,dsy,d2sy);
      state.m_f=CApServ::SafePythag2(dsx,dsy);
     }
//--- function call
   CAutoGK::AutoGKResults(state,result,rep);
//--- check
   if(!CAp::Assert(rep.m_terminationtype>0,__FUNCTION__+": internal error!"))
      return(EMPTY_VALUE);
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| This function calculates arc length, i.e. length of curve between|
//| t=a and t=b.                                                     |
//| INPUT PARAMETERS:                                                |
//|     P   -   parametric spline interpolant                        |
//|     A,B -   parameter values corresponding to arc ends:          |
//|             * B>A will result in positive length returned        |
//|             * B<A will result in negative length returned        |
//| RESULT:                                                          |
//|     length of arc starting at T=A and ending at T=B.             |
//+------------------------------------------------------------------+
static double CPSpline::PSpline3ArcLength(CPSpline3Interpolant &p,const double a,
                                          const double b)
  {
//--- create variables
   double result=0;
   double sx=0;
   double dsx=0;
   double d2sx=0;
   double sy=0;
   double dsy=0;
   double d2sy=0;
   double sz=0;
   double dsz=0;
   double d2sz=0;
//--- objects of classes
   CAutoGKState  state;
   CAutoGKReport rep;
//--- function call
   CAutoGK::AutoGKSmooth(a,b,state);
//--- cycle
   while(CAutoGK::AutoGKIteration(state))
     {
      CSpline1D::Spline1DDiff(p.m_x,state.m_x,sx,dsx,d2sx);
      CSpline1D::Spline1DDiff(p.m_y,state.m_x,sy,dsy,d2sy);
      CSpline1D::Spline1DDiff(p.m_z,state.m_x,sz,dsz,d2sz);
      state.m_f=CApServ::SafePythag3(dsx,dsy,dsz);
     }
//--- function call
   CAutoGK::AutoGKResults(state,result,rep);
//--- check
   if(!CAp::Assert(rep.m_terminationtype>0,__FUNCTION__+": internal error!"))
      return(EMPTY_VALUE);
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Builds non-periodic parameterization for 2-dimensional spline    |
//+------------------------------------------------------------------+
static void CPSpline::PSpline2Par(CMatrixDouble &xy,const int n,const int pt,
                                  double &p[])
  {
//--- create variables
   double v=0;
   int    i=0;
   int    i_=0;
//--- check
   if(!CAp::Assert(pt>=0 && pt<=2,__FUNCTION__+": internal error!"))
      return;
//--- Build parameterization:
//--- * fill by non-normalized values
//--- * normalize them so we have P[0]=0,P[N-1]=1.
   ArrayResizeAL(p,n);
//--- check
   if(pt==0)
     {
      for(i=0;i<=n-1;i++)
         p[i]=i;
     }
//--- check
   if(pt==1)
     {
      p[0]=0;
      //--- calculation
      for(i=1;i<=n-1;i++)
         p[i]=p[i-1]+CApServ::SafePythag2(xy[i][0]-xy[i-1][0],xy[i][1]-xy[i-1][1]);
     }
//--- check
   if(pt==2)
     {
      p[0]=0;
      //--- calculation
      for(i=1;i<=n-1;i++)
         p[i]=p[i-1]+MathSqrt(CApServ::SafePythag2(xy[i][0]-xy[i-1][0],xy[i][1]-xy[i-1][1]));
     }
//--- change value
   v=1/p[n-1];
//--- calculation
   for(i_=0;i_<=n-1;i_++)
      p[i_]=v*p[i_];
  }
//+------------------------------------------------------------------+
//| Builds non-periodic parameterization for 3-dimensional spline    |
//+------------------------------------------------------------------+
static void CPSpline::PSpline3Par(CMatrixDouble &xy,const int n,const int pt,
                                  double &p[])
  {
//--- create variables
   double v=0;
   int    i=0;
   int    i_=0;
//--- check
   if(!CAp::Assert(pt>=0 && pt<=2,__FUNCTION__+": internal error!"))
      return;
//--- Build parameterization:
//--- * fill by non-normalized values
//--- * normalize them so we have P[0]=0,P[N-1]=1.
   ArrayResizeAL(p,n);
//--- check
   if(pt==0)
     {
      for(i=0;i<=n-1;i++)
         p[i]=i;
     }
//--- check
   if(pt==1)
     {
      p[0]=0;
      //--- calculation
      for(i=1;i<=n-1;i++)
         p[i]=p[i-1]+CApServ::SafePythag3(xy[i][0]-xy[i-1][0],xy[i][1]-xy[i-1][1],xy[i][2]-xy[i-1][2]);
     }
//--- check
   if(pt==2)
     {
      p[0]=0;
      //--- calculation
      for(i=1;i<=n-1;i++)
         p[i]=p[i-1]+MathSqrt(CApServ::SafePythag3(xy[i][0]-xy[i-1][0],xy[i][1]-xy[i-1][1],xy[i][2]-xy[i-1][2]));
     }
//--- change value
   v=1/p[n-1];
//--- calculation
   for(i_=0;i_<=n-1;i_++)
      p[i_]=v*p[i_];
  }
//+------------------------------------------------------------------+
//| 2-dimensional spline inteprolant                                 |
//+------------------------------------------------------------------+
class CSpline2DInterpolant
  {
public:
   //--- variable
   int               m_k;
   //--- array
   double            m_c[];
   //--- constructor, destructor
                     CSpline2DInterpolant(void);
                    ~CSpline2DInterpolant(void);
   //--- copy
   void              Copy(CSpline2DInterpolant &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSpline2DInterpolant::CSpline2DInterpolant(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSpline2DInterpolant::~CSpline2DInterpolant(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CSpline2DInterpolant::Copy(CSpline2DInterpolant &obj)
  {
//--- copy variable
   m_k=obj.m_k;
//--- copy array
   ArrayCopy(m_c,obj.m_c);
  }
//+------------------------------------------------------------------+
//| 2-dimensional spline inteprolant                                 |
//+------------------------------------------------------------------+
class CSpline2DInterpolantShell
  {
private:
   CSpline2DInterpolant m_innerobj;
public:
   //--- constructors, destructor
                     CSpline2DInterpolantShell(void);
                     CSpline2DInterpolantShell(CSpline2DInterpolant &obj);
                    ~CSpline2DInterpolantShell(void);
   //--- method
   CSpline2DInterpolant *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSpline2DInterpolantShell::CSpline2DInterpolantShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CSpline2DInterpolantShell::CSpline2DInterpolantShell(CSpline2DInterpolant &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSpline2DInterpolantShell::~CSpline2DInterpolantShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CSpline2DInterpolant *CSpline2DInterpolantShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| 2-dimensional spline interpolation                               |
//+------------------------------------------------------------------+
class CSpline2D
  {
private:
   //--- private method
   static void       BicubicCalcDerivatives(CMatrixDouble &a,double &x[],double &y[],const int m,const int n,CMatrixDouble &dx,CMatrixDouble &dy,CMatrixDouble &dxy);
public:
   //--- constructor, destructor
                     CSpline2D(void);
                    ~CSpline2D(void);
   //--- public methods
   static void       Spline2DBuildBilinear(double &cx[],double &cy[],CMatrixDouble &cf,const int m,const int n,CSpline2DInterpolant &c);
   static void       Spline2DBuildBicubic(double &cx[],double &cy[],CMatrixDouble &cf,const int m,const int n,CSpline2DInterpolant &c);
   static double     Spline2DCalc(CSpline2DInterpolant &c,const double x,const double y);
   static void       Spline2DDiff(CSpline2DInterpolant &c,const double x,const double y,double &f,double &fx,double &fy,double &fxy);
   static void       Spline2DUnpack(CSpline2DInterpolant &c,int &m,int &n,CMatrixDouble &tbl);
   static void       Spline2DLinTransXY(CSpline2DInterpolant &c,double ax,double bx,double ay,double by);
   static void       Spline2DLinTransF(CSpline2DInterpolant &c,const double a,const double b);
   static void       Spline2DCopy(CSpline2DInterpolant &c,CSpline2DInterpolant &cc);
   static void       Spline2DResampleBicubic(CMatrixDouble &a,const int oldheight,const int oldwidth,CMatrixDouble &b,const int newheight,const int newwidth);
   static void       Spline2DResampleBilinear(CMatrixDouble &a,const int oldheight,const int oldwidth,CMatrixDouble &b,const int newheight,const int newwidth);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSpline2D::CSpline2D(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSpline2D::~CSpline2D(void)
  {

  }
//+------------------------------------------------------------------+
//| This subroutine builds bilinear spline coefficients table.       |
//| Input parameters:                                                |
//|     X   -   spline abscissas, array[0..N-1]                      |
//|     Y   -   spline ordinates, array[0..M-1]                      |
//|     F   -   function values, array[0..M-1,0..N-1]                |
//|     M,N -   grid size, M>=2, N>=2                                |
//| Output parameters:                                               |
//|     C   -   spline interpolant                                   |
//+------------------------------------------------------------------+
static void CSpline2D::Spline2DBuildBilinear(double &cx[],double &cy[],
                                             CMatrixDouble &cf,const int m,
                                             const int n,CSpline2DInterpolant &c)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    tblsize=0;
   int    shift=0;
   double t=0;
//--- create arrays
   double x[];
   double y[];
//--- create matrix
   CMatrixDouble dx;
   CMatrixDouble dy;
   CMatrixDouble dxy;
   CMatrixDouble f;
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
//--- copy matrix
   f=cf;
//--- check
   if(!CAp::Assert(n>=2 && m>=2,__FUNCTION__+": N<2 or M<2!"))
      return;
//--- Sort points
   for(j=0;j<=n-1;j++)
     {
      k=j;
      for(i=j+1;i<=n-1;i++)
        {
         //--- check
         if(x[i]<x[k])
            k=i;
        }
      //--- check
      if(k!=j)
        {
         for(i=0;i<=m-1;i++)
           {
            //--- swap
            t=f[i][j];
            f[i].Set(j,f[i][k]);
            f[i].Set(k,t);
           }
         //--- swap
         t=x[j];
         x[j]=x[k];
         x[k]=t;
        }
     }
//--- calculation
   for(i=0;i<=m-1;i++)
     {
      k=i;
      for(j=i+1;j<=m-1;j++)
        {
         //--- check
         if(y[j]<y[k])
            k=j;
        }
      //--- check
      if(k!=i)
        {
         for(j=0;j<=n-1;j++)
           {
            //--- swap
            t=f[i][j];
            f[i].Set(j,f[k][j]);
            f[k].Set(j,t);
           }
         //--- swap
         t=y[i];
         y[i]=y[k];
         y[k]=t;
        }
     }
//--- Fill C:
//---  C[0]            -   length(C)
//---  C[1]            -   type(C):
//---                      -1=bilinear interpolant
//---                      -3=general cubic spline
//---                           (see BuildBicubicSpline)
//---  C[2]:
//---      N (x count)
//---  C[3]:
//---      M (y count)
//---  C[4]...C[4+N-1]:
//---      x[i],i=0...N-1
//---  C[4+N]...C[4+N+M-1]:
//---      y[i],i=0...M-1
//---  C[4+N+M]...C[4+N+M+(N*M-1)]:
//---      f(i,j) table. f(0,0),f(0,1),f(0,2) and so on...
   c.m_k=1;
   tblsize=4+n+m+n*m;
//--- allocation
   ArrayResizeAL(c.m_c,tblsize);
//--- change values
   c.m_c[0]=tblsize;
   c.m_c[1]=-1;
   c.m_c[2]=n;
   c.m_c[3]=m;
//--- copy
   for(i=0;i<=n-1;i++)
      c.m_c[4+i]=x[i];
   for(i=0;i<=m-1;i++)
      c.m_c[4+n+i]=y[i];
//--- calculation
   for(i=0;i<=m-1;i++)
     {
      for(j=0;j<=n-1;j++)
        {
         shift=i*n+j;
         c.m_c[4+n+m+shift]=f[i][j];
        }
     }
  }
//+------------------------------------------------------------------+
//| This subroutine builds bicubic spline coefficients table.        |
//| Input parameters:                                                |
//|     X   -   spline abscissas, array[0..N-1]                      |
//|     Y   -   spline ordinates, array[0..M-1]                      |
//|     F   -   function values, array[0..M-1,0..N-1]                |
//|     M,N -   grid size, M>=2, N>=2                                |
//| Output parameters:                                               |
//|     C   -   spline interpolant                                   |
//+------------------------------------------------------------------+
static void CSpline2D::Spline2DBuildBicubic(double &cx[],double &cy[],
                                            CMatrixDouble &cf,const int m,
                                            const int n,CSpline2DInterpolant &c)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    tblsize=0;
   int    shift=0;
   double t=0;
//--- create arrays
   double x[];
   double y[];
//--- create matrix
   CMatrixDouble dx;
   CMatrixDouble dy;
   CMatrixDouble dxy;
   CMatrixDouble f;
//--- copy arrays
   ArrayCopy(x,cx);
   ArrayCopy(y,cy);
//--- copy matrix
   f=cf;
//--- check
   if(!CAp::Assert(n>=2 && m>=2,__FUNCTION__+": N<2 or M<2!"))
      return;
//--- Sort points
   for(j=0;j<=n-1;j++)
     {
      k=j;
      for(i=j+1;i<=n-1;i++)
        {
         //--- check
         if(x[i]<x[k])
            k=i;
        }
      //--- check
      if(k!=j)
        {
         for(i=0;i<=m-1;i++)
           {
            //--- swap
            t=f[i][j];
            f[i].Set(j,f[i][k]);
            f[i].Set(k,t);
           }
         //--- swap
         t=x[j];
         x[j]=x[k];
         x[k]=t;
        }
     }
//--- calculation
   for(i=0;i<=m-1;i++)
     {
      k=i;
      for(j=i+1;j<=m-1;j++)
        {
         //--- check
         if(y[j]<y[k])
            k=j;
        }
      //--- check
      if(k!=i)
        {
         for(j=0;j<=n-1;j++)
           {
            //--- swap
            t=f[i][j];
            f[i].Set(j,f[k][j]);
            f[k].Set(j,t);
           }
         //--- swap
         t=y[i];
         y[i]=y[k];
         y[k]=t;
        }
     }
//--- Fill C:
//---  C[0]            -   length(C)
//---  C[1]            -   type(C):
//---                      -1=bilinear interpolant
//---                           (see BuildBilinearInterpolant)
//---                      -3=general cubic spline
//---  C[2]:
//---      N (x count)
//---  C[3]:
//---      M (y count)
//---  C[4]...C[4+N-1]:
//---      x[i],i=0...N-1
//---  C[4+N]...C[4+N+M-1]:
//---      y[i],i=0...M-1
//---  C[4+N+M]...C[4+N+M+(N*M-1)]:
//---      f(i,j) table. f(0,0),f(0,1),f(0,2) and so on...
//---  C[4+N+M+N*M]...C[4+N+M+(2*N*M-1)]:
//---      df(i,j)/dx table.
//---  C[4+N+M+2*N*M]...C[4+N+M+(3*N*M-1)]:
//---      df(i,j)/dy table.
//---  C[4+N+M+3*N*M]...C[4+N+M+(4*N*M-1)]:
//---      d2f(i,j)/dxdy table.
   c.m_k=3;
   tblsize=4+n+m+4*n*m;
//--- allocation
   ArrayResizeAL(c.m_c,tblsize);
//--- change values
   c.m_c[0]=tblsize;
   c.m_c[1]=-3;
   c.m_c[2]=n;
   c.m_c[3]=m;
//--- copy
   for(i=0;i<=n-1;i++)
      c.m_c[4+i]=x[i];
   for(i=0;i<=m-1;i++)
      c.m_c[4+n+i]=y[i];
//--- function call
   BicubicCalcDerivatives(f,x,y,m,n,dx,dy,dxy);
//--- change values
   for(i=0;i<=m-1;i++)
     {
      for(j=0;j<=n-1;j++)
        {
         shift=i*n+j;
         c.m_c[4+n+m+shift]=f[i][j];
         c.m_c[4+n+m+n*m+shift]=dx[i][j];
         c.m_c[4+n+m+2*n*m+shift]=dy[i][j];
         c.m_c[4+n+m+3*n*m+shift]=dxy[i][j];
        }
     }
  }
//+------------------------------------------------------------------+
//| This subroutine calculates the value of the bilinear or bicubic  |
//| spline at the given point X.                                     |
//| Input parameters:                                                |
//|     C   -   coefficients table.                                  |
//|             Built by BuildBilinearSpline or BuildBicubicSpline.  |
//|     X, Y-   point                                                |
//| Result:                                                          |
//|     S(x,y)                                                       |
//+------------------------------------------------------------------+
static double CSpline2D::Spline2DCalc(CSpline2DInterpolant &c,const double x,
                                      const double y)
  {
//--- create variables
   double v=0;
   double vx=0;
   double vy=0;
   double vxy=0;
//--- function call
   Spline2DDiff(c,x,y,v,vx,vy,vxy);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| This subroutine calculates the value of the bilinear or bicubic  |
//| spline at the given point X and its derivatives.                 |
//| Input parameters:                                                |
//|     C   -   spline interpolant.                                  |
//|     X, Y-   point                                                |
//| Output parameters:                                               |
//|     F   -   S(x,y)                                               |
//|     FX  -   dS(x,y)/dX                                           |
//|     FY  -   dS(x,y)/dY                                           |
//|     FXY -   d2S(x,y)/dXdY                                        |
//+------------------------------------------------------------------+
static void CSpline2D::Spline2DDiff(CSpline2DInterpolant &c,const double x,
                                    const double y,double &f,double &fx,
                                    double &fy,double &fxy)
  {
//--- create variables
   int    n=0;
   int    m=0;
   double t=0;
   double dt=0;
   double u=0;
   double du=0;
   int    ix=0;
   int    iy=0;
   int    l=0;
   int    r=0;
   int    h=0;
   int    shift1=0;
   int    s1=0;
   int    s2=0;
   int    s3=0;
   int    s4=0;
   int    sf=0;
   int    sfx=0;
   int    sfy=0;
   int    sfxy=0;
   double y1=0;
   double y2=0;
   double y3=0;
   double y4=0;
   double v=0;
   double t0=0;
   double t1=0;
   double t2=0;
   double t3=0;
   double u0=0;
   double u1=0;
   double u2=0;
   double u3=0;
//--- initialization
   f=0;
   fx=0;
   fy=0;
   fxy=0;
//--- check
   if(!CAp::Assert((int)MathRound(c.m_c[1])==-1 || (int)MathRound(c.m_c[1])==-3,__FUNCTION__+": incorrect C!"))
      return;
//--- initialization
   n=(int)MathRound(c.m_c[2]);
   m=(int)MathRound(c.m_c[3]);
//--- Binary search in the [ x[0],...,x[n-2] ] (x[n-1] is not included)
   l=4;
   r=4+n-2+1;
   while(l!=r-1)
     {
      h=(l+r)/2;
      //--- check
      if(c.m_c[h]>=x)
         r=h;
      else
         l=h;
     }
//--- change values
   t=(x-c.m_c[l])/(c.m_c[l+1]-c.m_c[l]);
   dt=1.0/(c.m_c[l+1]-c.m_c[l]);
   ix=l-4;
//--- Binary search in the [ y[0],...,y[m-2] ] (y[m-1] is not included)
   l=4+n;
   r=4+n+(m-2)+1;
   while(l!=r-1)
     {
      h=(l+r)/2;
      //--- check
      if(c.m_c[h]>=y)
         r=h;
      else
         l=h;
     }
//--- change values
   u=(y-c.m_c[l])/(c.m_c[l+1]-c.m_c[l]);
   du=1.0/(c.m_c[l+1]-c.m_c[l]);
   iy=l-(4+n);
//--- Prepare F,dF/dX,dF/dY,d2F/dXdY
   f=0;
   fx=0;
   fy=0;
   fxy=0;
//--- Bilinear interpolation
   if((int)MathRound(c.m_c[1])==-1)
     {
      //--- calculation
      shift1=4+n+m;
      y1=c.m_c[shift1+n*iy+ix];
      y2=c.m_c[shift1+n*iy+(ix+1)];
      y3=c.m_c[shift1+n*(iy+1)+(ix+1)];
      y4=c.m_c[shift1+n*(iy+1)+ix];
      f=(1-t)*(1-u)*y1+t*(1-u)*y2+t*u*y3+(1-t)*u*y4;
      fx=(-((1-u)*y1)+(1-u)*y2+u*y3-u*y4)*dt;
      fy=(-((1-t)*y1)-t*y2+t*y3+(1-t)*y4)*du;
      fxy=(y1-y2+y3-y4)*du*dt;
      //--- exit the function
      return;
     }
//--- Bicubic interpolation
   if((int)MathRound(c.m_c[1])==-3)
     {
      //--- Prepare info
      t0=1;
      t1=t;
      t2=CMath::Sqr(t);
      t3=t*t2;
      u0=1;
      u1=u;
      u2=CMath::Sqr(u);
      u3=u*u2;
      sf=4+n+m;
      sfx=4+n+m+n*m;
      sfy=4+n+m+2*n*m;
      sfxy=4+n+m+3*n*m;
      s1=n*iy+ix;
      s2=n*iy+(ix+1);
      s3=n*(iy+1)+(ix+1);
      s4=n*(iy+1)+ix;
      //--- Calculate
      v=1*c.m_c[sf+s1];
      f=f+v*t0*u0;
      v=1*c.m_c[sfy+s1]/du;
      f=f+v*t0*u1;
      fy=fy+1*v*t0*u0*du;
      v=-(3*c.m_c[sf+s1])+3*c.m_c[sf+s4]-2*c.m_c[sfy+s1]/du-1*c.m_c[sfy+s4]/du;
      f=f+v*t0*u2;
      fy=fy+2*v*t0*u1*du;
      v=2*c.m_c[sf+s1]-2*c.m_c[sf+s4]+1*c.m_c[sfy+s1]/du+1*c.m_c[sfy+s4]/du;
      f=f+v*t0*u3;
      fy=fy+3*v*t0*u2*du;
      v=1*c.m_c[sfx+s1]/dt;
      f=f+v*t1*u0;
      fx=fx+1*v*t0*u0*dt;
      v=1*c.m_c[sfxy+s1]/(dt*du);
      f=f+v*t1*u1;
      fx=fx+1*v*t0*u1*dt;
      fy=fy+1*v*t1*u0*du;
      fxy=fxy+1*v*t0*u0*dt*du;
      v=-(3*c.m_c[sfx+s1]/dt)+3*c.m_c[sfx+s4]/dt-2*c.m_c[sfxy+s1]/(dt*du)-1*c.m_c[sfxy+s4]/(dt*du);
      f=f+v*t1*u2;
      fx=fx+1*v*t0*u2*dt;
      fy=fy+2*v*t1*u1*du;
      fxy=fxy+2*v*t0*u1*dt*du;
      v=2*c.m_c[sfx+s1]/dt-2*c.m_c[sfx+s4]/dt+1*c.m_c[sfxy+s1]/(dt*du)+1*c.m_c[sfxy+s4]/(dt*du);
      f=f+v*t1*u3;
      fx=fx+1*v*t0*u3*dt;
      fy=fy+3*v*t1*u2*du;
      fxy=fxy+3*v*t0*u2*dt*du;
      v=-(3*c.m_c[sf+s1])+3*c.m_c[sf+s2]-2*c.m_c[sfx+s1]/dt-1*c.m_c[sfx+s2]/dt;
      f=f+v*t2*u0;
      fx=fx+2*v*t1*u0*dt;
      v=-(3*c.m_c[sfy+s1]/du)+3*c.m_c[sfy+s2]/du-2*c.m_c[sfxy+s1]/(dt*du)-1*c.m_c[sfxy+s2]/(dt*du);
      f=f+v*t2*u1;
      fx=fx+2*v*t1*u1*dt;
      fy=fy+1*v*t2*u0*du;
      fxy=fxy+2*v*t1*u0*dt*du;
      v=9*c.m_c[sf+s1]-9*c.m_c[sf+s2]+9*c.m_c[sf+s3]-9*c.m_c[sf+s4]+6*c.m_c[sfx+s1]/dt+3*c.m_c[sfx+s2]/dt-3*c.m_c[sfx+s3]/dt-6*c.m_c[sfx+s4]/dt+6*c.m_c[sfy+s1]/du-6*c.m_c[sfy+s2]/du-3*c.m_c[sfy+s3]/du+3*c.m_c[sfy+s4]/du+4*c.m_c[sfxy+s1]/(dt*du)+2*c.m_c[sfxy+s2]/(dt*du)+1*c.m_c[sfxy+s3]/(dt*du)+2*c.m_c[sfxy+s4]/(dt*du);
      f=f+v*t2*u2;
      fx=fx+2*v*t1*u2*dt;
      fy=fy+2*v*t2*u1*du;
      fxy=fxy+4*v*t1*u1*dt*du;
      v=-(6*c.m_c[sf+s1])+6*c.m_c[sf+s2]-6*c.m_c[sf+s3]+6*c.m_c[sf+s4]-4*c.m_c[sfx+s1]/dt-2*c.m_c[sfx+s2]/dt+2*c.m_c[sfx+s3]/dt+4*c.m_c[sfx+s4]/dt-3*c.m_c[sfy+s1]/du+3*c.m_c[sfy+s2]/du+3*c.m_c[sfy+s3]/du-3*c.m_c[sfy+s4]/du-2*c.m_c[sfxy+s1]/(dt*du)-1*c.m_c[sfxy+s2]/(dt*du)-1*c.m_c[sfxy+s3]/(dt*du)-2*c.m_c[sfxy+s4]/(dt*du);
      f=f+v*t2*u3;
      fx=fx+2*v*t1*u3*dt;
      fy=fy+3*v*t2*u2*du;
      fxy=fxy+6*v*t1*u2*dt*du;
      v=2*c.m_c[sf+s1]-2*c.m_c[sf+s2]+1*c.m_c[sfx+s1]/dt+1*c.m_c[sfx+s2]/dt;
      f=f+v*t3*u0;
      fx=fx+3*v*t2*u0*dt;
      v=2*c.m_c[sfy+s1]/du-2*c.m_c[sfy+s2]/du+1*c.m_c[sfxy+s1]/(dt*du)+1*c.m_c[sfxy+s2]/(dt*du);
      f=f+v*t3*u1;
      fx=fx+3*v*t2*u1*dt;
      fy=fy+1*v*t3*u0*du;
      fxy=fxy+3*v*t2*u0*dt*du;
      v=-(6*c.m_c[sf+s1])+6*c.m_c[sf+s2]-6*c.m_c[sf+s3]+6*c.m_c[sf+s4]-3*c.m_c[sfx+s1]/dt-3*c.m_c[sfx+s2]/dt+3*c.m_c[sfx+s3]/dt+3*c.m_c[sfx+s4]/dt-4*c.m_c[sfy+s1]/du+4*c.m_c[sfy+s2]/du+2*c.m_c[sfy+s3]/du-2*c.m_c[sfy+s4]/du-2*c.m_c[sfxy+s1]/(dt*du)-2*c.m_c[sfxy+s2]/(dt*du)-1*c.m_c[sfxy+s3]/(dt*du)-1*c.m_c[sfxy+s4]/(dt*du);
      f=f+v*t3*u2;
      fx=fx+3*v*t2*u2*dt;
      fy=fy+2*v*t3*u1*du;
      fxy=fxy+6*v*t2*u1*dt*du;
      v=4*c.m_c[sf+s1]-4*c.m_c[sf+s2]+4*c.m_c[sf+s3]-4*c.m_c[sf+s4]+2*c.m_c[sfx+s1]/dt+2*c.m_c[sfx+s2]/dt-2*c.m_c[sfx+s3]/dt-2*c.m_c[sfx+s4]/dt+2*c.m_c[sfy+s1]/du-2*c.m_c[sfy+s2]/du-2*c.m_c[sfy+s3]/du+2*c.m_c[sfy+s4]/du+1*c.m_c[sfxy+s1]/(dt*du)+1*c.m_c[sfxy+s2]/(dt*du)+1*c.m_c[sfxy+s3]/(dt*du)+1*c.m_c[sfxy+s4]/(dt*du);
      f=f+v*t3*u3;
      fx=fx+3*v*t2*u3*dt;
      fy=fy+3*v*t3*u2*du;
      fxy=fxy+9*v*t2*u2*dt*du;
      //--- exit the function
      return;
     }
  }
//+------------------------------------------------------------------+
//| This subroutine unpacks two-dimensional spline into the          |
//| coefficients table                                               |
//| Input parameters:                                                |
//|     C   -   spline interpolant.                                  |
//| Result:                                                          |
//|     M, N-   grid size (x-axis and y-axis)                        |
//|     Tbl -   coefficients table, unpacked format,                 |
//|             [0..(N-1)*(M-1)-1, 0..19].                           |
//|             For I = 0...M-2, J=0..N-2:                           |
//|                 K =  I*(N-1)+J                                   |
//|                 Tbl[K,0] = X[j]                                  |
//|                 Tbl[K,1] = X[j+1]                                |
//|                 Tbl[K,2] = Y[i]                                  |
//|                 Tbl[K,3] = Y[i+1]                                |
//|                 Tbl[K,4] = C00                                   |
//|                 Tbl[K,5] = C01                                   |
//|                 Tbl[K,6] = C02                                   |
//|                 Tbl[K,7] = C03                                   |
//|                 Tbl[K,8] = C10                                   |
//|                 Tbl[K,9] = C11                                   |
//|                 ...                                              |
//|                 Tbl[K,19] = C33                                  |
//|             On each grid square spline is equals to:             |
//|                 S(x) = SUM(c[i,j]*(x^i)*(y^j), i=0..3, j=0..3)   |
//|                 t = x-x[j]                                       |
//|                 u = y-y[i]                                       |
//+------------------------------------------------------------------+
static void CSpline2D::Spline2DUnpack(CSpline2DInterpolant &c,int &m,
                                      int &n,CMatrixDouble &tbl)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    ci=0;
   int    cj=0;
   int    k=0;
   int    p=0;
   int    shift=0;
   int    s1=0;
   int    s2=0;
   int    s3=0;
   int    s4=0;
   int    sf=0;
   int    sfx=0;
   int    sfy=0;
   int    sfxy=0;
   double y1=0;
   double y2=0;
   double y3=0;
   double y4=0;
   double dt=0;
   double du=0;
//--- initialization
   m=0;
   n=0;
//--- check
   if(!CAp::Assert((int)MathRound(c.m_c[1])==-3 || (int)MathRound(c.m_c[1])==-1,__FUNCTION__+": incorrect C!"))
      return;
//--- initialization
   n=(int)MathRound(c.m_c[2]);
   m=(int)MathRound(c.m_c[3]);
//--- allocation
   tbl.Resize((n-1)*(m-1),20);
//--- Fill
   for(i=0;i<=m-2;i++)
     {
      for(j=0;j<=n-2;j++)
        {
         //--- calculation
         p=i*(n-1)+j;
         tbl[p].Set(0,c.m_c[4+j]);
         tbl[p].Set(1,c.m_c[4+j+1]);
         tbl[p].Set(2,c.m_c[4+n+i]);
         tbl[p].Set(3,c.m_c[4+n+i+1]);
         dt=1/(tbl[p][1]-tbl[p][0]);
         du=1/(tbl[p][3]-tbl[p][2]);
         //--- Bilinear interpolation
         if((int)MathRound(c.m_c[1])==-1)
           {
            for(k=4;k<=19;k++)
               tbl[p].Set(k,0);
            //--- calculation
            shift=4+n+m;
            y1=c.m_c[shift+n*i+j];
            y2=c.m_c[shift+n*i+(j+1)];
            y3=c.m_c[shift+n*(i+1)+(j+1)];
            y4=c.m_c[shift+n*(i+1)+j];
            tbl[p].Set(4,y1);
            tbl[p].Set(4+1*4+0,y2-y1);
            tbl[p].Set(4+0*4+1,y4-y1);
            tbl[p].Set(4+1*4+1,y3-y2-y4+y1);
           }
         //--- Bicubic interpolation
         if((int)MathRound(c.m_c[1])==-3)
           {
            //--- change values
            sf=4+n+m;
            sfx=4+n+m+n*m;
            sfy=4+n+m+2*n*m;
            sfxy=4+n+m+3*n*m;
            s1=n*i+j;
            s2=n*i+(j+1);
            s3=n*(i+1)+(j+1);
            s4=n*(i+1)+j;
            //--- change values
            tbl[p].Set(4+0*4+0,1*c.m_c[sf+s1]);
            tbl[p].Set(4+0*4+1,1*c.m_c[sfy+s1]/du);
            tbl[p].Set(4+0*4+2,-(3*c.m_c[sf+s1])+3*c.m_c[sf+s4]-2*c.m_c[sfy+s1]/du-1*c.m_c[sfy+s4]/du);
            tbl[p].Set(4+0*4+3,2*c.m_c[sf+s1]-2*c.m_c[sf+s4]+1*c.m_c[sfy+s1]/du+1*c.m_c[sfy+s4]/du);
            tbl[p].Set(4+1*4+0,1*c.m_c[sfx+s1]/dt);
            tbl[p].Set(4+1*4+1,1*c.m_c[sfxy+s1]/(dt*du));
            tbl[p].Set(4+1*4+2,-(3*c.m_c[sfx+s1]/dt)+3*c.m_c[sfx+s4]/dt-2*c.m_c[sfxy+s1]/(dt*du)-1*c.m_c[sfxy+s4]/(dt*du));
            tbl[p].Set(4+1*4+3,2*c.m_c[sfx+s1]/dt-2*c.m_c[sfx+s4]/dt+1*c.m_c[sfxy+s1]/(dt*du)+1*c.m_c[sfxy+s4]/(dt*du));
            tbl[p].Set(4+2*4+0,-(3*c.m_c[sf+s1])+3*c.m_c[sf+s2]-2*c.m_c[sfx+s1]/dt-1*c.m_c[sfx+s2]/dt);
            tbl[p].Set(4+2*4+1,-(3*c.m_c[sfy+s1]/du)+3*c.m_c[sfy+s2]/du-2*c.m_c[sfxy+s1]/(dt*du)-1*c.m_c[sfxy+s2]/(dt*du));
            tbl[p].Set(4+2*4+2,9*c.m_c[sf+s1]-9*c.m_c[sf+s2]+9*c.m_c[sf+s3]-9*c.m_c[sf+s4]+6*c.m_c[sfx+s1]/dt+3*c.m_c[sfx+s2]/dt-3*c.m_c[sfx+s3]/dt-6*c.m_c[sfx+s4]/dt+6*c.m_c[sfy+s1]/du-6*c.m_c[sfy+s2]/du-3*c.m_c[sfy+s3]/du+3*c.m_c[sfy+s4]/du+4*c.m_c[sfxy+s1]/(dt*du)+2*c.m_c[sfxy+s2]/(dt*du)+1*c.m_c[sfxy+s3]/(dt*du)+2*c.m_c[sfxy+s4]/(dt*du));
            tbl[p].Set(4+2*4+3,-(6*c.m_c[sf+s1])+6*c.m_c[sf+s2]-6*c.m_c[sf+s3]+6*c.m_c[sf+s4]-4*c.m_c[sfx+s1]/dt-2*c.m_c[sfx+s2]/dt+2*c.m_c[sfx+s3]/dt+4*c.m_c[sfx+s4]/dt-3*c.m_c[sfy+s1]/du+3*c.m_c[sfy+s2]/du+3*c.m_c[sfy+s3]/du-3*c.m_c[sfy+s4]/du-2*c.m_c[sfxy+s1]/(dt*du)-1*c.m_c[sfxy+s2]/(dt*du)-1*c.m_c[sfxy+s3]/(dt*du)-2*c.m_c[sfxy+s4]/(dt*du));
            tbl[p].Set(4+3*4+0,2*c.m_c[sf+s1]-2*c.m_c[sf+s2]+1*c.m_c[sfx+s1]/dt+1*c.m_c[sfx+s2]/dt);
            tbl[p].Set(4+3*4+1,2*c.m_c[sfy+s1]/du-2*c.m_c[sfy+s2]/du+1*c.m_c[sfxy+s1]/(dt*du)+1*c.m_c[sfxy+s2]/(dt*du));
            tbl[p].Set(4+3*4+2,-(6*c.m_c[sf+s1])+6*c.m_c[sf+s2]-6*c.m_c[sf+s3]+6*c.m_c[sf+s4]-3*c.m_c[sfx+s1]/dt-3*c.m_c[sfx+s2]/dt+3*c.m_c[sfx+s3]/dt+3*c.m_c[sfx+s4]/dt-4*c.m_c[sfy+s1]/du+4*c.m_c[sfy+s2]/du+2*c.m_c[sfy+s3]/du-2*c.m_c[sfy+s4]/du-2*c.m_c[sfxy+s1]/(dt*du)-2*c.m_c[sfxy+s2]/(dt*du)-1*c.m_c[sfxy+s3]/(dt*du)-1*c.m_c[sfxy+s4]/(dt*du));
            tbl[p].Set(4+3*4+3,4*c.m_c[sf+s1]-4*c.m_c[sf+s2]+4*c.m_c[sf+s3]-4*c.m_c[sf+s4]+2*c.m_c[sfx+s1]/dt+2*c.m_c[sfx+s2]/dt-2*c.m_c[sfx+s3]/dt-2*c.m_c[sfx+s4]/dt+2*c.m_c[sfy+s1]/du-2*c.m_c[sfy+s2]/du-2*c.m_c[sfy+s3]/du+2*c.m_c[sfy+s4]/du+1*c.m_c[sfxy+s1]/(dt*du)+1*c.m_c[sfxy+s2]/(dt*du)+1*c.m_c[sfxy+s3]/(dt*du)+1*c.m_c[sfxy+s4]/(dt*du));
           }
         //--- Rescale Cij
         for(ci=0;ci<=3;ci++)
           {
            for(cj=0;cj<=3;cj++)
               tbl[p].Set(4+ci*4+cj,tbl[p][4+ci*4+cj]*MathPow(dt,ci)*MathPow(du,cj));
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| This subroutine performs linear transformation of the spline     |
//| argument.                                                        |
//| Input parameters:                                                |
//|     C       -   spline interpolant                               |
//|     AX, BX  -   transformation coefficients: x = A*t + B         |
//|     AY, BY  -   transformation coefficients: y = A*u + B         |
//| Result:                                                          |
//|     C   -   transformed spline                                   |
//+------------------------------------------------------------------+
static void CSpline2D::Spline2DLinTransXY(CSpline2DInterpolant &c,double ax,
                                          double bx,double ay,double by)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    n=0;
   int    m=0;
   double v=0;
   int    typec=0;
//--- create arrays
   double x[];
   double y[];
//--- create matrix
   CMatrixDouble f;
//--- initialization
   typec=(int)MathRound(c.m_c[1]);
//--- check
   if(!CAp::Assert(typec==-3 || typec==-1,__FUNCTION__+": incorrect C!"))
      return;
//--- initialization
   n=(int)MathRound(c.m_c[2]);
   m=(int)MathRound(c.m_c[3]);
//--- allocation
   ArrayResizeAL(x,n);
   ArrayResizeAL(y,m);
   f.Resize(m,n);
//--- copy
   for(j=0;j<=n-1;j++)
      x[j]=c.m_c[4+j];
   for(i=0;i<=m-1;i++)
      y[i]=c.m_c[4+n+i];
   for(i=0;i<=m-1;i++)
     {
      for(j=0;j<=n-1;j++)
         f[i].Set(j,c.m_c[4+n+m+i*n+j]);
     }
//--- Special case: AX=0 or AY=0
   if(ax==0.0)
     {
      //--- change values
      for(i=0;i<=m-1;i++)
        {
         v=Spline2DCalc(c,bx,y[i]);
         for(j=0;j<=n-1;j++)
            f[i].Set(j,v);
        }
      //--- check
      if(typec==-3)
         Spline2DBuildBicubic(x,y,f,m,n,c);
      //--- check
      if(typec==-1)
         Spline2DBuildBilinear(x,y,f,m,n,c);
      //--- change values
      ax=1;
      bx=0;
     }
//--- check
   if(ay==0.0)
     {
      //--- change values
      for(j=0;j<=n-1;j++)
        {
         v=Spline2DCalc(c,x[j],by);
         for(i=0;i<=m-1;i++)
            f[i].Set(j,v);
        }
      //--- check
      if(typec==-3)
         Spline2DBuildBicubic(x,y,f,m,n,c);
      //--- check
      if(typec==-1)
         Spline2DBuildBilinear(x,y,f,m,n,c);
      //--- change values
      ay=1;
      by=0;
     }
//--- General case: AX<>0,AY<>0
//--- Unpack,scale and pack again.
   for(j=0;j<=n-1;j++)
      x[j]=(x[j]-bx)/ax;
   for(i=0;i<=m-1;i++)
      y[i]=(y[i]-by)/ay;
//--- check
   if(typec==-3)
      Spline2DBuildBicubic(x,y,f,m,n,c);
//--- check
   if(typec==-1)
      Spline2DBuildBilinear(x,y,f,m,n,c);
  }
//+------------------------------------------------------------------+
//| This subroutine performs linear transformation of the spline.    |
//| Input parameters:                                                |
//|     C   -   spline interpolant.                                  |
//|     A, B-   transformation coefficients: S2(x,y) = A*S(x,y) + B  |
//| Output parameters:                                               |
//|     C   -   transformed spline                                   |
//+------------------------------------------------------------------+
static void CSpline2D::Spline2DLinTransF(CSpline2DInterpolant &c,const double a,
                                         const double b)
  {
//--- create variables
   int i=0;
   int j=0;
   int n=0;
   int m=0;
   int typec=0;
//--- create arrays
   double x[];
   double y[];
//--- create matrix
   CMatrixDouble f;
//--- initialization
   typec=(int)MathRound(c.m_c[1]);
//--- check
   if(!CAp::Assert(typec==-3 || typec==-1,__FUNCTION__+": incorrect C!"))
      return;
//--- initialization
   n=(int)MathRound(c.m_c[2]);
   m=(int)MathRound(c.m_c[3]);
//--- allocation
   ArrayResizeAL(x,n);
   ArrayResizeAL(y,m);
   f.Resize(m,n);
//--- copy
   for(j=0;j<=n-1;j++)
      x[j]=c.m_c[4+j];
   for(i=0;i<=m-1;i++)
      y[i]=c.m_c[4+n+i];
   for(i=0;i<=m-1;i++)
     {
      for(j=0;j<=n-1;j++)
         f[i].Set(j,a*c.m_c[4+n+m+i*n+j]+b);
     }
//--- check
   if(typec==-3)
      Spline2DBuildBicubic(x,y,f,m,n,c);
//--- check
   if(typec==-1)
      Spline2DBuildBilinear(x,y,f,m,n,c);
  }
//+------------------------------------------------------------------+
//| This subroutine makes the copy of the spline model.              |
//| Input parameters:                                                |
//|     C   -   spline interpolant                                   |
//| Output parameters:                                               |
//|     CC  -   spline copy                                          |
//+------------------------------------------------------------------+
static void CSpline2D::Spline2DCopy(CSpline2DInterpolant &c,CSpline2DInterpolant &cc)
  {
//--- create variables
   int n=0;
   int i_=0;
//--- check
   if(!CAp::Assert(c.m_k==1 || c.m_k==3,__FUNCTION__+": incorrect C!"))
      return;
//--- change values
   cc.m_k=c.m_k;
   n=(int)MathRound(c.m_c[0]);
//--- allocation
   ArrayResizeAL(cc.m_c,n);
//--- copy
   for(i_=0;i_<=n-1;i_++)
      cc.m_c[i_]=c.m_c[i_];
  }
//+------------------------------------------------------------------+
//| Bicubic spline resampling                                        |
//| Input parameters:                                                |
//|     A           -   function values at the old grid,             |
//|                     array[0..OldHeight-1, 0..OldWidth-1]         |
//|     OldHeight   -   old grid height, OldHeight>1                 |
//|     OldWidth    -   old grid width, OldWidth>1                   |
//|     NewHeight   -   new grid height, NewHeight>1                 |
//|     NewWidth    -   new grid width, NewWidth>1                   |
//| Output parameters:                                               |
//|     B           -   function values at the new grid,             |
//|                     array[0..NewHeight-1, 0..NewWidth-1]         |
//+------------------------------------------------------------------+
static void CSpline2D::Spline2DResampleBicubic(CMatrixDouble &a,const int oldheight,
                                               const int oldwidth,CMatrixDouble &b,
                                               const int newheight,const int newwidth)
  {
//--- create variables
   int i=0;
   int j=0;
   int mw=0;
   int mh=0;
//--- create arrays
   double x[];
   double y[];
//--- create matrix
   CMatrixDouble buf;
//--- object of class
   CSpline1DInterpolant c;
//--- check
   if(!CAp::Assert(oldwidth>1 && oldheight>1,__FUNCTION__+": width/height less than 1"))
      return;
//--- check
   if(!CAp::Assert(newwidth>1 && newheight>1,__FUNCTION__+": width/height less than 1"))
      return;
//--- Prepare
   mw=MathMax(oldwidth,newwidth);
   mh=MathMax(oldheight,newheight);
//--- allocation
   b.Resize(newheight,newwidth);
   buf.Resize(oldheight,newwidth);
   ArrayResizeAL(x,MathMax(mw,mh));
   ArrayResizeAL(y,MathMax(mw,mh));
//--- Horizontal interpolation
   for(i=0;i<=oldheight-1;i++)
     {
      //--- Fill X,Y
      for(j=0;j<=oldwidth-1;j++)
        {
         x[j]=(double)j/(double)(oldwidth-1);
         y[j]=a[i][j];
        }
      //--- Interpolate and place result into temporary matrix
      CSpline1D::Spline1DBuildCubic(x,y,oldwidth,0,0.0,0,0.0,c);
      for(j=0;j<=newwidth-1;j++)
         buf[i].Set(j,CSpline1D::Spline1DCalc(c,(double)j/(double)(newwidth-1)));
     }
//--- Vertical interpolation
   for(j=0;j<=newwidth-1;j++)
     {
      //--- Fill X,Y
      for(i=0;i<=oldheight-1;i++)
        {
         x[i]=(double)i/(double)(oldheight-1);
         y[i]=buf[i][j];
        }
      //--- Interpolate and place result into B
      CSpline1D::Spline1DBuildCubic(x,y,oldheight,0,0.0,0,0.0,c);
      for(i=0;i<=newheight-1;i++)
         b[i].Set(j,CSpline1D::Spline1DCalc(c,(double)i/(double)(newheight-1)));
     }
  }
//+------------------------------------------------------------------+
//| Bilinear spline resampling                                       |
//| Input parameters:                                                |
//|     A           -   function values at the old grid,             |
//|                     array[0..OldHeight-1, 0..OldWidth-1]         |
//|     OldHeight   -   old grid height, OldHeight>1                 |
//|     OldWidth    -   old grid width, OldWidth>1                   |
//|     NewHeight   -   new grid height, NewHeight>1                 |
//|     NewWidth    -   new grid width, NewWidth>1                   |
//| Output parameters:                                               |
//|     B           -   function values at the new grid,             |
//|                     array[0..NewHeight-1, 0..NewWidth-1]         |
//+------------------------------------------------------------------+
static void CSpline2D::Spline2DResampleBilinear(CMatrixDouble &a,const int oldheight,
                                                const int oldwidth,CMatrixDouble &b,
                                                const int newheight,const int newwidth)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    l=0;
   int    c=0;
   double t=0;
   double u=0;
//--- allocation
   b.Resize(newheight,newwidth);
   for(i=0;i<=newheight-1;i++)
     {
      for(j=0;j<=newwidth-1;j++)
        {
         //--- calculation
         l=i*(oldheight-1)/(newheight-1);
         //--- check
         if(l==oldheight-1)
            l=oldheight-2;
         //--- calculation
         u=(double)i/(double)(newheight-1)*(oldheight-1)-l;
         c=j*(oldwidth-1)/(newwidth-1);
         //--- check
         if(c==oldwidth-1)
            c=oldwidth-2;
         //--- calculation
         t=(double)(j*(oldwidth-1))/(double)(newwidth-1)-c;
         b[i].Set(j,(1-t)*(1-u)*a[l][c]+t*(1-u)*a[l][c+1]+t*u*a[l+1][c+1]+(1-t)*u*a[l+1][c]);
        }
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine.                                             |
//| Calculation of the first derivatives and the cross-derivative.   |
//+------------------------------------------------------------------+
static void CSpline2D::BicubicCalcDerivatives(CMatrixDouble &a,double &x[],
                                              double &y[],const int m,
                                              const int n,CMatrixDouble &dx,
                                              CMatrixDouble &dy,CMatrixDouble &dxy)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double s=0;
   double ds=0;
   double d2s=0;
//--- create arrays
   double xt[];
   double ft[];
//--- object of class
   CSpline1DInterpolant c;
//--- allocation
   dx.Resize(m,n);
   dy.Resize(m,n);
   dxy.Resize(m,n);
//--- dF/dX
   ArrayResizeAL(xt,n);
   ArrayResizeAL(ft,n);
   for(i=0;i<=m-1;i++)
     {
      for(j=0;j<=n-1;j++)
        {
         xt[j]=x[j];
         ft[j]=a[i][j];
        }
      //--- function call
      CSpline1D::Spline1DBuildCubic(xt,ft,n,0,0.0,0,0.0,c);
      for(j=0;j<=n-1;j++)
        {
         //--- function call
         CSpline1D::Spline1DDiff(c,x[j],s,ds,d2s);
         dx[i].Set(j,ds);
        }
     }
//--- dF/dY
   ArrayResizeAL(xt,m);
   ArrayResizeAL(ft,m);
   for(j=0;j<=n-1;j++)
     {
      for(i=0;i<=m-1;i++)
        {
         xt[i]=y[i];
         ft[i]=a[i][j];
        }
      //--- function call
      CSpline1D::Spline1DBuildCubic(xt,ft,m,0,0.0,0,0.0,c);
      for(i=0;i<=m-1;i++)
        {
         //--- function call
         CSpline1D::Spline1DDiff(c,y[i],s,ds,d2s);
         dy[i].Set(j,ds);
        }
     }
//--- d2F/dXdY
   ArrayResizeAL(xt,n);
   ArrayResizeAL(ft,n);
   for(i=0;i<=m-1;i++)
     {
      for(j=0;j<=n-1;j++)
        {
         xt[j]=x[j];
         ft[j]=dy[i][j];
        }
      //--- function call
      CSpline1D::Spline1DBuildCubic(xt,ft,n,0,0.0,0,0.0,c);
      for(j=0;j<=n-1;j++)
        {
         //--- function call
         CSpline1D::Spline1DDiff(c,x[j],s,ds,d2s);
         dxy[i].Set(j,ds);
        }
     }
  }
//+------------------------------------------------------------------+
