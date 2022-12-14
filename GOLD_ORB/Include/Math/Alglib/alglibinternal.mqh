//+------------------------------------------------------------------+
//|                                               alglibinternal.mqh |
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
//| published by the Free Software Foundation (www.fsf.org);either   |
//| version 2 of the License, or (at your option) any later version. |
//|                                                                  |
//| This program is distributed in the hope that it will be useful,  |
//| but WITHOUT ANY WARRANTY;without even the implied warranty of    |
//| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the     |
//| GNU General Public License for more details.                     |
//+------------------------------------------------------------------+
#include "ap.mqh"
//+------------------------------------------------------------------+
//| Class stores serialized codes                                    |
//+------------------------------------------------------------------+
class CSCodes
  {
public:
                     CSCodes(void);
                    ~CSCodes(void);

   static int        GetRDFSerializationCode(void)    { return(1); }
   static int        GetKDTreeSerializationCode(void) { return(2); }
   static int        GetMLPSerializationCode(void)    { return(3); }
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSCodes::CSCodes(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSCodes::~CSCodes(void)
  {

  }
//+------------------------------------------------------------------+
//| Buffers for internal functions which need buffers:               |
//| * check for size of the buffer you want to use.                  |
//| * if buffer is too small, resize it; leave unchanged, if it is   |
//| larger than needed.                                              |
//| * use it.                                                        |
//| We can pass this structure to multiple functions;  after first   |
//| run through functions buffer sizes will be finally determined,   |
//| and  on  a next run no allocation will be required.              |
//+------------------------------------------------------------------+
class CApBuff
  {
public:
   //--- arrays
   int               m_ia0[];
   int               m_ia1[];
   int               m_ia2[];
   int               m_ia3[];
   double            m_ra0[];
   double            m_ra1[];
   double            m_ra2[];
   double            m_ra3[];
   //--- constructor, destructor
                     CApBuff(void);
                    ~CApBuff(void);
   //--- copy
   void              Copy(CApBuff &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CApBuff::CApBuff(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CApBuff::~CApBuff(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CApBuff::Copy(CApBuff &obj)
  {
//--- copy arrays
   ArrayCopy(m_ia0,obj.m_ia0);
   ArrayCopy(m_ia1,obj.m_ia1);
   ArrayCopy(m_ia2,obj.m_ia2);
   ArrayCopy(m_ia3,obj.m_ia3);
   ArrayCopy(m_ra0,obj.m_ra0);
   ArrayCopy(m_ra1,obj.m_ra1);
   ArrayCopy(m_ra2,obj.m_ra2);
   ArrayCopy(m_ra3,obj.m_ra3);
  }
//+------------------------------------------------------------------+
//| Basic functions                                                  |
//+------------------------------------------------------------------+
class CApServ
  {
public:
                     CApServ(void);
                    ~CApServ(void);
   //--- generate interpolation
   static void       TaskGenInt1D(const double a,const double b,const int n,double &x[],double &y[]);
   static void       TaskGenInt1DEquidist(const double a,const double b,const int n,double &x[],double &y[]);
   static void       TaskGenInt1DCheb1(const double a,const double b,const int n,double &x[],double &y[]);
   static void       TaskGenInt1DCheb2(const double a,const double b,const int n,double &x[],double &y[]);
   //--- distinct
   static bool       AreDistinct(double &x[],const int n);
   //--- resize arrays
   static void       BVectorSetLengthAtLeast(bool &x[],const int n);
   static void       IVectorSetLengthAtLeast(int &x[],const int n);
   static void       RVectorSetLengthAtLeast(double &x[],const int n);
   static void       RMatrixSetLengthAtLeast(CMatrixDouble &x,const int m,const int n);
   //--- resize matrix
   static void       RMatrixResize(CMatrixDouble &x,const int m,const int n);
   //--- check to infinity
   static bool       IsFiniteVector(const double &x[],const int n);
   static bool       IsFiniteComplexVector(al_complex &z[],const int n);
   static bool       IsFiniteMatrix(const CMatrixDouble &x,const int m,const int n);
   static bool       IsFiniteComplexMatrix(CMatrixComplex &x,const int m,const int n);
   static bool       IsFiniteRTrMatrix(CMatrixDouble &x,const int n,const bool isupper);
   static bool       IsFiniteCTrMatrix(CMatrixComplex &x,const int n,const bool isupper);
   static bool       IsFiniteOrNaNMatrix(CMatrixDouble &x,const int m,const int n);
   //--- safe methods
   static double     SafePythag2(const double x,const double y);
   static double     SafePythag3(double x,double y,double z);
   static int        SafeRDiv(double x,double y,double &r);
   static double     SafeMinPosRV(const double x,const double y,const double v);
   static void       ApPeriodicMap(double &x,const double a,const double b,double &k);
   static double     BoundVal(const double x,const double b1,const double b2);
   //--- serialization/unserialization
   static void       AllocComplex(CSerializer &s,al_complex &v);
   static void       SerializeComplex(CSerializer &s,al_complex &v);
   static al_complex UnserializeComplex(CSerializer &s);
   static void       AllocRealArray(CSerializer &s,double &v[],int n);
   static void       SerializeRealArray(CSerializer &s,double &v[],int n);
   static void       UnserializeRealArray(CSerializer &s,double &v[]);
   static void       AllocIntegerArray(CSerializer &s,int &v[],int n);
   static void       SerializeIntegerArray(CSerializer &s,int &v[],int n);
   static void       UnserializeIntegerArray(CSerializer &s,int &v[]);
   static void       AllocRealMatrix(CSerializer &s,CMatrixDouble &v,int n0,int n1);
   static void       SerializeRealMatrix(CSerializer &s,CMatrixDouble &v,int n0,int n1);
   static void       UnserializeRealMatrix(CSerializer &s,CMatrixDouble &v);
   //--- copy
   static void       CopyIntegerArray(int &src[],int &dst[]);
   static void       CopyRealArray(double &src[],double &dst[]);
   static void       CopyRealMatrix(CMatrixDouble &src,CMatrixDouble &dst);
   //--- check array
   static int        RecSearch(int &a[],const int nrec,const int nheader,int i0,int i1,int &b[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CApServ::CApServ(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CApServ::~CApServ(void)
  {

  }
//+------------------------------------------------------------------+
//| This  function  generates  1-dimensional  general  interpolation |
//| task with moderate Lipshitz constant (close to 1.0)              |
//| If N=1 then suborutine generates only one point at the middle    |
//| of [A,B]                                                         |
//+------------------------------------------------------------------+
static void CApServ::TaskGenInt1D(const double a,const double b,const int n,
                                  double &x[],double &y[])
  {
//--- create variables
   int    i=0;
   double h=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- allocation
   ArrayResizeAL(x,n);
   ArrayResizeAL(y,n);
//--- check
   if(n>1)
     {
      //--- change values
      x[0]=a;
      y[0]=2*CMath::RandomReal()-1;
      h=(b-a)/(n-1);
      for(i=1;i<=n-1;i++)
        {
         //--- check
         if(i!=n-1)
            x[i]=a+(i+0.2*(2*CMath::RandomReal()-1))*h;
         else
            x[i]=b;
         y[i]=y[i-1]+(2*CMath::RandomReal()-1)*(x[i]-x[i-1]);
        }
     }
   else
     {
      //--- change values
      x[0]=0.5*(a+b);
      y[0]=2*CMath::RandomReal()-1;
     }
  }
//+------------------------------------------------------------------+
//| This function generates  1-dimensional equidistant interpolation |
//| task withmoderate Lipshitz constant(close to 1.0)                |
//| If N=1 then suborutine generates only one point at the middle    |
//| of[A,B]                                                          |
//+------------------------------------------------------------------+
static void CApServ::TaskGenInt1DEquidist(const double a,const double b,
                                          const int n,double &x[],double &y[])
  {
//--- create variables
   int    i=0;
   double h=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- allocation
   ArrayResizeAL(x,n);
   ArrayResizeAL(y,n);
//--- check
   if(n>1)
     {
      //--- change values
      x[0]=a;
      y[0]=2*CMath::RandomReal()-1;
      h=(b-a)/(n-1);
      for(i=1;i<=n-1;i++)
        {
         x[i]=a+i*h;
         y[i]=y[i-1]+(2*CMath::RandomReal()-1)*h;
        }
     }
   else
     {
      //--- change values
      x[0]=0.5*(a+b);
      y[0]=2*CMath::RandomReal()-1;
     }
  }
//+------------------------------------------------------------------+
//| This function generates  1-dimensional Chebyshev-1 interpolation |
//| task with moderate Lipshitz constant(close to 1.0)               |
//| If N=1 then suborutine generates only one point at the middle    |
//| of[A,B]                                                          |
//+------------------------------------------------------------------+
static void CApServ::TaskGenInt1DCheb1(const double a,const double b,
                                       const int n,double &x[],double &y[])
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- allocation
   ArrayResizeAL(x,n);
   ArrayResizeAL(y,n);
//--- check
   if(n>1)
     {
      for(i=0;i<=n-1;i++)
        {
         x[i]=0.5*(b+a)+0.5*(b-a)*MathCos(M_PI*(2*i+1)/(2*n));
         //--- check
         if(i==0)
            y[i]=2*CMath::RandomReal()-1;
         else
            y[i]=y[i-1]+(2*CMath::RandomReal()-1)*(x[i]-x[i-1]);
        }
     }
   else
     {
      //--- change values
      x[0]=0.5*(a+b);
      y[0]=2*CMath::RandomReal()-1;
     }
  }
//+------------------------------------------------------------------+
//| This function generates  1-dimensional Chebyshev-2 interpolation |
//| task with moderate Lipshitz constant(close to 1.0)               |
//| If N=1 then suborutine generates only one point at the middle    |
//| of[A,B]                                                          |
//+------------------------------------------------------------------+
static void CApServ::TaskGenInt1DCheb2(const double a,const double b,
                                       const int n,double &x[],double &y[])
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- allocation
   ArrayResizeAL(x,n);
   ArrayResizeAL(y,n);
//--- check
   if(n>1)
     {
      for(i=0;i<=n-1;i++)
        {
         x[i]=0.5*(b+a)+0.5*(b-a)*MathCos(M_PI*i/(n-1));
         //--- check
         if(i==0)
            y[i]=2*CMath::RandomReal()-1;
         else
            y[i]=y[i-1]+(2*CMath::RandomReal()-1)*(x[i]-x[i-1]);
        }
     }
   else
     {
      //--- change values
      x[0]=0.5*(a+b);
      y[0]=2*CMath::RandomReal()-1;
     }
  }
//+------------------------------------------------------------------+
//| This function checks that all values from X[] are distinct.      |
//| It does more than just usual floating point comparison:          |
//| * first, it calculates max(X) and min(X)                         |
//| * second, it maps X[] from [min,max] to [1,2]                    |
//| * only at this stage actual comparison is done                   |
//| The meaning of such check is to ensure that all values are       |
//| "distinct enough" and will not cause interpolation subroutine    |
//| to fail.                                                         |
//|  NOTE:                                                           |
//|     X[] must be sorted by ascending (subroutine ASSERT's it)     |
//+------------------------------------------------------------------+
static bool CApServ::AreDistinct(double &x[],const int n)
  {
//--- create variables
   double a=0;
   double b=0;
   int    i=0;
   bool   nonsorted;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": internal error (N<1)"))
      return(false);
//--- check
   if(n==1)
     {
      //--- everything is alright, it is up to caller to decide whether it
      //--- can interpolate something with just one point
      return(true);
     }
//--- initialization
   a=x[0];
   b=x[0];
   nonsorted=false;
   for(i=1;i<=n-1;i++)
     {
      a=MathMin(a,x[i]);
      b=MathMax(b,x[i]);
      nonsorted=nonsorted || x[i-1]>=x[i];
     }
//--- check
   if(!CAp::Assert(!nonsorted,__FUNCTION__+": internal error (not sorted)"))
      return(false);
   for(i=1;i<=n-1;i++)
     {
      //--- check
      if((x[i]-a)/(b-a)+1==(x[i-1]-a)/(b-a)+1)
         return(false);
     }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| If Length(X)<N, resizes X                                        |
//+------------------------------------------------------------------+
static void CApServ::BVectorSetLengthAtLeast(bool &x[],const int n)
  {
//--- check
   if(CAp::Len(x)<n)
      ArrayResizeAL(x,n);
  }
//+------------------------------------------------------------------+
//| If Length(X)<N, resizes X                                        |
//+------------------------------------------------------------------+
static void CApServ::IVectorSetLengthAtLeast(int &x[],const int n)
  {
//--- check
   if(CAp::Len(x)<n)
      ArrayResizeAL(x,n);
  }
//+------------------------------------------------------------------+
//| If Length(X)<N , resizes X                                       |
//+------------------------------------------------------------------+
static void CApServ::RVectorSetLengthAtLeast(double &x[],const int n)
  {
//--- check
   if(CAp::Len(x)<n)
      ArrayResizeAL(x,n);
  }
//+------------------------------------------------------------------+
//| If Cols(X)<N or Rows(X)<M, resizes X                             |
//+------------------------------------------------------------------+
static void CApServ::RMatrixSetLengthAtLeast(CMatrixDouble &x,const int m,
                                             const int n)
  {
//--- check
   if(CAp::Rows(x)<m || CAp::Cols(x)<n)
      x.Resize(m,n);
  }
//+------------------------------------------------------------------+
//| Resizes X and:                                                   |
//| * preserves old contents of X                                    |
//| * fills new elements by zeros                                    |
//+------------------------------------------------------------------+
static void CApServ::RMatrixResize(CMatrixDouble &x,const int m,const int n)
  {
//--- create variables
   int i=0;
   int j=0;
   int m2=0;
   int n2=0;
//--- create matrix
   CMatrixDouble oldx;
//--- initialization
   m2=CAp::Rows(x);
   n2=CAp::Cols(x);
//--- function call
   CAp::Swap(x,oldx);
//--- resize
   x.Resize(m,n);
//--- filling
   for(i=0;i<=m-1;i++)
     {
      for(j=0;j<=n-1;j++)
        {
         //--- check
         if(i<m2 && j<n2)
            x[i].Set(j,oldx[i][j]);
         else
            x[i].Set(j,0.0);
        }
     }
  }
//+------------------------------------------------------------------+
//| This function checks that all values from X[] are finite         |
//+------------------------------------------------------------------+
static bool CApServ::IsFiniteComplexVector(al_complex &z[],const int n)
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(n>=0,__FUNCTION__+": internal error (N<0)"))
      return(false);
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(!CMath::IsFinite(z[i].re) || !CMath::IsFinite(z[i].im))
         return(false);
     }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| This function checks that all values from X[] are finite         |
//+------------------------------------------------------------------+
static bool CApServ::IsFiniteVector(const double &x[],const int n)
  {
//--- check
   if(!CAp::Assert(n>=0,__FUNCTION__+": the error variable"))
      return(false);
//--- is finite?
   for(int i=0;i<n;i++)
      if(!CMath::IsFinite(x[i]))
         return(false);
//--- is finite
   return(true);
  }
//+------------------------------------------------------------------+
//| This function checks that all values from X[0..M-1,0..N-1]       |
//| are finite                                                       |
//+------------------------------------------------------------------+
static bool CApServ::IsFiniteMatrix(const CMatrixDouble &x,const int m,
                                    const int n)
  {
//--- check
   if(!CAp::Assert(n>=0,__FUNCTION__+": the error variable"))
      return(false);
//--- check
   if(!CAp::Assert(m>=0,__FUNCTION__+": the error variable"))
      return(false);
//--- is finite?
   for(int i=0;i<m;i++)
      for(int j=0;j<n;j++)
         //--- check
         if(!CMath::IsFinite(x[i][j]))
            return(false);
//--- is finite
   return(true);
  }
//+------------------------------------------------------------------+
//| This function checks that all values from X[0..M-1,0..N-1]       |
//| are finite                                                       |
//+------------------------------------------------------------------+
static bool CApServ::IsFiniteComplexMatrix(CMatrixComplex &x,const int m,
                                           const int n)
  {
//--- create variables
   int  i=0;
   int  j=0;
//--- check
   if(!CAp::Assert(n>=0,__FUNCTION__+": internal error (N<0)"))
      return(false);
//--- check
   if(!CAp::Assert(m>=0,__FUNCTION__+": internal error (M<0)"))
      return(false);
//--- is finite?
   for(i=0;i<m;i++)
     {
      for(j=0;j<n;j++)
         //--- check
         if(!CMath::IsFinite(x[i][j].re) || !CMath::IsFinite(x[i][j].im))
            return(false);
     }
//--- is finite
   return(true);
  }
//+------------------------------------------------------------------+
//| This function checks that all values from upper/lower triangle of|
//| X[0..N-1,0..N-1] are finite                                      |
//+------------------------------------------------------------------+
static bool CApServ::IsFiniteRTrMatrix(CMatrixDouble &x,const int n,
                                       const bool isupper)
  {
//--- create variables
   int i=0;
   int j1=0;
   int j2=0;
   int j=0;
//--- check
   if(!CAp::Assert(n>=0,__FUNCTION__+": internal error (N<0)"))
      return(false);
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
      for(j=j1;j<=j2;j++)
        {
         //--- check
         if(!CMath::IsFinite(x[i][j]))
            return(false);
        }
     }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| This function checks that all values from upper/lower triangle of|
//| X[0..N-1,0..N-1] are finite                                      |
//+------------------------------------------------------------------+
static bool CApServ::IsFiniteCTrMatrix(CMatrixComplex &x,const int n,
                                       const bool isupper)
  {
//--- create variables
   int i=0;
   int j1=0;
   int j2=0;
   int j=0;
//--- check
   if(!CAp::Assert(n>=0,__FUNCTION__+": internal error (N<0)"))
      return(false);
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
      for(j=j1; j<=j2; j++)
        {
         //--- check
         if(!CMath::IsFinite(x[i][j].re) || !CMath::IsFinite(x[i][j].im))
            return(false);
        }
     }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| This function checks that all values from X[0..M-1,0..N-1] are   |
//| finite or NaN's.                                                 |
//+------------------------------------------------------------------+
static bool CApServ::IsFiniteOrNaNMatrix(CMatrixDouble &x,const int m,
                                         const int n)
  {
//--- create variables
   int i=0;
   int j=0;
//--- check
   if(!CAp::Assert(n>=0,__FUNCTION__+": internal error (N<0)"))
      return(false);
//--- check
   if(!CAp::Assert(m>=0,__FUNCTION__+": internal error (M<0)"))
      return(false);
   for(i=0;i<=m-1;i++)
     {
      for(j=0;j<=n-1;j++)
        {
         //--- check
         if(!(CMath::IsFinite(x[i][j]) || CInfOrNaN::IsNaN(x[i][j])))
            return(false);
        }
     }
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Safe sqrt(x^2+y^2)                                               |
//+------------------------------------------------------------------+
static double CApServ::SafePythag2(const double x,const double y)
  {
//--- create variables
   double result=0;
   double w=0;
   double xabs=0;
   double yabs=0;
   double z=0;
//--- initialization
   xabs=MathAbs(x);
   yabs=MathAbs(y);
   w=MathMax(xabs,yabs);
   z=MathMin(xabs,yabs);
//--- check
   if(z==0.0)
      result=w;
   else
      result=w*MathSqrt(1+CMath::Sqr(z/w));
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Safe sqrt(x^2+y^2)                                               |
//+------------------------------------------------------------------+
static double CApServ::SafePythag3(double x,double y,double z)
  {
//--- create a variable
   double w=0;
//--- initialization
   w=MathMax(MathAbs(x),MathMax(MathAbs(y),MathAbs(z)));
//--- check
   if(w==0.0)
      return(0);
//--- change values
   x=x/w;
   y=y/w;
   z=z/w;
//--- return result
   return(w*MathSqrt(CMath::Sqr(x)+CMath::Sqr(y)+CMath::Sqr(z)));
  }
//+------------------------------------------------------------------+
//| Safe division.                                                   |
//| This function attempts to calculate R=X/Y without overflow.      |
//| It returns:                                                      |
//| * +1, if abs(X/Y)>=MaxRealNumber or undefined - overflow-like    |
//|       situation (no overlfow is generated, R is either NAN,      |
//|       PosINF, NegINF)                                            |
//| *  0, if MinRealNumber<abs(X/Y)<MaxRealNumber or X=0, Y<>0       |
//|       (R contains result, may be zero)                           |
//| * -1, if 0<abs(X/Y)<MinRealNumber - underflow-like situation     |
//|       (R contains zero; it corresponds to underflow)             |
//| No overflow is generated in any case.                            |
//+------------------------------------------------------------------+
static int CApServ::SafeRDiv(double x,double y,double &r)
  {
//--- create variables
   int result=0;
//--- initialization
   r=0;
//--- Two special cases:
//--- * Y=0
//--- * X=0 and Y<>0
   if(y==0.0)
     {
      result=1;
      //--- check
      if(x==0.0)
         r=CInfOrNaN::NaN();
      //--- check
      if(x>0.0)
         r=CInfOrNaN::PositiveInfinity();
      //--- check
      if(x<0.0)
         r=CInfOrNaN::NegativeInfinity();
      //--- return result
      return(result);
     }
//--- check
   if(x==0.0)
     {
      r=0;
      result=0;
      //--- return result
      return(result);
     }
//--- make Y>0
   if(y<0.0)
     {
      x=-x;
      y=-y;
     }
//--- check
   if(y>=1.0)
     {
      r=x/y;
      //--- check
      if(MathAbs(r)<=CMath::m_minrealnumber)
        {
         result=-1;
         r=0;
        }
      else
         result=0;
     }
   else
     {
      //--- check
      if(MathAbs(x)>=CMath::m_maxrealnumber*y)
        {
         //--- check
         if(x>0.0)
            r=CInfOrNaN::PositiveInfinity();
         else
            r=CInfOrNaN::NegativeInfinity();
         result=1;
        }
      else
        {
         r=x/y;
         result=0;
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| This function calculates "safe" min(X/Y,V) for positive finite X,|
//| Y, V. No overflow is generated in any case.                      |
//+------------------------------------------------------------------+
static double CApServ::SafeMinPosRV(const double x,const double y,const double v)
  {
//--- create variables
   double result=0;
   double r=0;
//--- check
   if(y>=1.0)
     {
      //--- Y>=1, we can safely divide by Y
      r=x/y;
      result=v;
      //--- check
      if(v>r)
         result=r;
      else
         result=v;
     }
   else
     {
      //--- Y<1, we can safely multiply by Y
      if(x<v*y)
         result=x/y;
      else
         result=v;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| This function makes periodic mapping of X to [A,B].              |
//| It accepts X, A, B (A>B). It returns T which lies in  [A,B] and  |
//| integer K, such that X = T + K*(B-A).                            |
//| NOTES:                                                           |
//| * K is represented as real value, although actually it is integer|
//| * T is guaranteed to be in [A,B]                                 |
//| * T replaces X                                                   |
//+------------------------------------------------------------------+
static void CApServ::ApPeriodicMap(double &x,const double a,const double b,
                                   double &k)
  {
//--- initialization
   k=0;
//--- check
   if(!CAp::Assert(a<b,__FUNCTION__+": internal error!"))
      return;
//--- initialization
   k=(int)MathFloor((x-a)/(b-a));
   x=x-k*(b-a);
//--- change values
   while(x<a)
     {
      x=x+(b-a);
      k=k-1;
     }
//--- change values
   while(x>b)
     {
      x=x-(b-a);
      k=k+1;
     }
//--- change values
   x=MathMax(x,a);
   x=MathMin(x,b);
  }
//+------------------------------------------------------------------+
//| 'bounds' value: maps X to [B1,B2]                                |
//+------------------------------------------------------------------+
static double CApServ::BoundVal(const double x,const double b1,const double b2)
  {
//--- check
   if(x<=b1)
      return(b1);
//--- check
   if(x>=b2)
      return(b2);
//--- return result
   return(x);
  }
//+------------------------------------------------------------------+
//| Allocation of serializer: complex value                          |
//+------------------------------------------------------------------+
static void CApServ::AllocComplex(CSerializer &s,al_complex &v)
  {
//--- entry
   s.Alloc_Entry();
   s.Alloc_Entry();
  }
//+------------------------------------------------------------------+
//| Serialization: complex value                                     |
//+------------------------------------------------------------------+
static void CApServ::SerializeComplex(CSerializer &s,al_complex &v)
  {
//--- serialization
   s.Serialize_Double(v.re);
   s.Serialize_Double(v.im);
  }
//+------------------------------------------------------------------+
//| Unserialization: complex value                                   |
//+------------------------------------------------------------------+
static al_complex CApServ::UnserializeComplex(CSerializer &s)
  {
//--- create a variable
   al_complex result;
//--- unserialization
   result.re=s.Unserialize_Double();
   result.im=s.Unserialize_Double();
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Allocation of serializer: real array                             |
//+------------------------------------------------------------------+
static void CApServ::AllocRealArray(CSerializer &s,double &v[],int n)
  {
//--- create a variable
   int i=0;
//--- check
   if(n<0)
      n=CAp::Len(v);
//--- entry
   s.Alloc_Entry();
   for(i=0;i<=n-1;i++)
      s.Alloc_Entry();
  }
//+------------------------------------------------------------------+
//| Serialization: real array                                        |
//+------------------------------------------------------------------+
static void CApServ::SerializeRealArray(CSerializer &s,double &v[],int n)
  {
//--- create a variable
   int i=0;
//--- check
   if(n<0)
      n=CAp::Len(v);
//--- serialization
   s.Serialize_Int(n);
   for(i=0;i<=n-1;i++)
      s.Serialize_Double(v[i]);
  }
//+------------------------------------------------------------------+
//| Unserialization: real array                                      |
//+------------------------------------------------------------------+
static void CApServ::UnserializeRealArray(CSerializer &s,double &v[])
  {
//--- create variables
   int    n=0;
   int    i=0;
   double t=0;
//--- unserialization
   n=s.Unserialize_Int();
//--- check
   if(n==0)
      return;
//--- allocation
   ArrayResizeAL(v,n);
//--- unserialization
   for(i=0;i<=n-1;i++)
     {
      t=s.Unserialize_Double();
      v[i]=t;
     }
  }
//+------------------------------------------------------------------+
//| Allocation of serializer: Integer array                          |
//+------------------------------------------------------------------+
static void CApServ::AllocIntegerArray(CSerializer &s,int &v[],int n)
  {
//--- create a variable
   int i=0;
//--- check
   if(n<0)
      n=CAp::Len(v);
//--- entry
   s.Alloc_Entry();
   for(i=0;i<=n-1;i++)
      s.Alloc_Entry();
  }
//+------------------------------------------------------------------+
//| Serialization: Integer array                                     |
//+------------------------------------------------------------------+
static void CApServ::SerializeIntegerArray(CSerializer &s,int &v[],int n)
  {
//--- create a variable
   int i=0;
//--- check
   if(n<0)
      n=CAp::Len(v);
//--- serialization
   s.Serialize_Int(n);
   for(i=0;i<=n-1;i++)
      s.Serialize_Int(v[i]);
  }
//+------------------------------------------------------------------+
//| Unserialization: Integer array                                   |
//+------------------------------------------------------------------+
static void CApServ::UnserializeIntegerArray(CSerializer &s,int &v[])
  {
//--- create variables
   int n=0;
   int i=0;
   int t=0;
//--- unserialization
   n=s.Unserialize_Int();
//--- check
   if(n==0)
      return;
//--- allocation
   ArrayResizeAL(v,n);
   for(i=0;i<=n-1;i++)
     {
      t=s.Unserialize_Int();
      v[i]=t;
     }
  }
//+------------------------------------------------------------------+
//| Allocation of serializer: real matrix                            |
//+------------------------------------------------------------------+
static void CApServ::AllocRealMatrix(CSerializer &s,CMatrixDouble &v,int n0,int n1)
  {
//--- create variables
   int i=0;
   int j=0;
//--- check
   if(n0<0)
      n0=CAp::Rows(v);
//--- check
   if(n1<0)
      n1=CAp::Cols(v);
//--- entry
   s.Alloc_Entry();
   s.Alloc_Entry();
   for(i=0;i<=n0-1;i++)
     {
      for(j=0;j<=n1-1;j++)
         s.Alloc_Entry();
     }
  }
//+------------------------------------------------------------------+
//| Serialization: real matrix                                       |
//+------------------------------------------------------------------+
static void CApServ::SerializeRealMatrix(CSerializer &s,CMatrixDouble &v,int n0,int n1)
  {
//--- create variables
   int i=0;
   int j=0;
//--- check
   if(n0<0)
      n0=CAp::Rows(v);
//--- check
   if(n1<0)
      n1=CAp::Cols(v);
//--- serialization
   s.Serialize_Int(n0);
//--- serialization
   s.Serialize_Int(n1);
   for(i=0;i<=n0-1;i++)
      for(j=0;j<=n1-1;j++)
         s.Serialize_Double(v[i][j]);
  }
//+------------------------------------------------------------------+
//| Unserialization: real matrix                                     |
//+------------------------------------------------------------------+
static void CApServ::UnserializeRealMatrix(CSerializer &s,CMatrixDouble &v)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    n0=0;
   int    n1=0;
   double t=0;
//--- unserialization
   n0=s.Unserialize_Int();
   n1=s.Unserialize_Int();
//--- check
   if(n0==0 || n1==0)
      return;
//--- resize
   v.Resize(n0,n1);
//--- unserialization
   for(i=0;i<=n0-1;i++)
      for(j=0;j<=n1-1;j++)
        {
         t=s.Unserialize_Double();
         v[i].Set(j,t);
        }
  }
//+------------------------------------------------------------------+
//| Copy integer array                                               |
//+------------------------------------------------------------------+
static void CApServ::CopyIntegerArray(int &src[],int &dst[])
  {
   int i=0;
//--- check
   if(CAp::Len(src)>0)
     {
      //--- allocation
      ArrayResizeAL(dst,CAp::Len(src));
      //--- copy
      for(i=0;i<=CAp::Len(src)-1;i++)
         dst[i]=src[i];
     }
  }
//+------------------------------------------------------------------+
//| Copy real array                                                  |
//+------------------------------------------------------------------+
static void CApServ::CopyRealArray(double &src[],double &dst[])
  {
   int i=0;
//--- check
   if(CAp::Len(src)>0)
     {
      ArrayResizeAL(dst,CAp::Len(src));
      for(i=0;i<=CAp::Len(src)-1;i++)
         dst[i]=src[i];
     }
  }
//+------------------------------------------------------------------+
//| Copy real matrix                                                 |
//+------------------------------------------------------------------+
static void CApServ::CopyRealMatrix(CMatrixDouble &src,CMatrixDouble &dst)
  {
//--- create variables
   int i=0;
   int j=0;
//--- check
   if(CAp::Rows(src)>0 && CAp::Cols(src)>0)
     {
      dst.Resize(CAp::Rows(src),CAp::Cols(src));
      //--- copy
      for(i=0;i<=CAp::Rows(src)-1;i++)
         for(j=0;j<=CAp::Cols(src)-1;j++)
            dst[i].Set(j,src[i][j]);
     }
  }
//+------------------------------------------------------------------+
//| This function searches integer array. Elements in this array are |
//| actually records, each NRec elements wide. Each record has unique|
//| header - NHeader integer values, which identify it. Records are  |
//| lexicographically sorted by header.                              |
//| Records are identified by their index, not offset                |
//| (offset = NRec*index).                                           |
//| This function searches A (records with indices [I0,I1)) for a    |
//| record with header B. It returns index of this record            |
//| (not offset!), or -1 on failure.                                 |
//+------------------------------------------------------------------+
static int CApServ::RecSearch(int &a[],const int nrec,const int nheader,int i0,int i1,int &b[])
  {
//--- create variables
   int mididx=0;
   int cflag=0;
   int k=0;
   int offs=0;
//--- cycle
   while(true)
     {
      //--- check
      if(i0>=i1)
         break;
      //--- change values
      mididx=(i0+i1)/2;
      offs=nrec*mididx;
      cflag=0;
      for(k=0;k<=nheader-1;k++)
        {
         //--- check
         if(a[offs+k]<b[k])
           {
            cflag=-1;
            break;
           }
         //--- check
         if(a[offs+k]>b[k])
           {
            cflag=1;
            break;
           }
        }
      //--- check
      if(cflag==0)
        {
         return(mididx);
        }
      //--- check
      if(cflag<0)
         i0=mididx+1;
      else
         i1=mididx;
     }
//--- return result
   return(-1);
  }
//+------------------------------------------------------------------+
//| Tag Sort                                                         |
//+------------------------------------------------------------------+
class CTSort
  {
private:
   //--- private methods
   static void       TagSortFastIRec(double &a[],int &b[],double &bufa[],int &bufb[],const int i1,const int i2);
   static void       TagSortFastRRec(double &a[],double &b[],double &bufa[],double &bufb[],const int i1,const int i2);
   static void       TagSortFastRec(double &a[],double &bufa[],const int i1,const int i2);
public:
                     CTSort(void);
                    ~CTSort(void);
   //--- public methods
   static void       TagSort(double &a[],const int n,int &p1[],int &p2[]);
   static void       TagSortBuf(double &a[],const int n,int &p1[],int &p2[],CApBuff &buf);
   static void       TagSortFastI(double &a[],int &b[],double &bufa[],int &bufb[],const int n);
   static void       TagSortFastR(double &a[],double &b[],double &bufa[],double &bufb[],const int n);
   static void       TagSortFast(double &a[],double &bufa[],const int n);
   static void       TagHeapPushI(double &a[],int &b[],int &n,const double va,const int vb);
   static void       TagHeapReplaceTopI(double &a[],int &b[],const int n,const double va,const int vb);
   static void       TagHeapPopI(double &a[],int &b[],int &n);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CTSort::CTSort(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CTSort::~CTSort(void)
  {

  }
//+------------------------------------------------------------------+
//| This function sorts array of real keys by ascending.             |
//| Its results are:                                                 |
//| * sorted array A                                                 |
//| * permutation tables P1, P2                                      |
//| Algorithm outputs permutation tables using two formats:          |
//| * as usual permutation of [0..N-1]. If P1[i]=j, then sorted A[i] |
//|   contains value which was moved there from J-th position.       |
//| * as a sequence of pairwise permutations. Sorted A[] may be      |
//|    obtained byswaping A[i] and A[P2[i]] for all i from 0 to N-1. |
//| INPUT PARAMETERS:                                                |
//|     A       -   unsorted array                                   |
//|     N       -   array size                                       |
//| OUPUT PARAMETERS:                                                |
//|     A       -   sorted array                                     |
//|     P1, P2  -   permutation tables, array[N]                     |
//| NOTES:                                                           |
//|     this function assumes that A[] is finite; it doesn't checks  |
//|     that condition. All other conditions (size of input arrays,  |
//|     etc.) are not checked too.                                   |
//+------------------------------------------------------------------+
static void CTSort::TagSort(double &a[],const int n,int &p1[],int &p2[])
  {
//--- create a variable
   CApBuff buf;
//--- function call
   TagSortBuf(a,n,p1,p2,buf);
  }
//+------------------------------------------------------------------+
//| Buffered variant of TagSort, which accepts preallocated output   |
//| arrays as well as special structure for buffered allocations. If |
//| arrays are too short, they are reallocated. If they are large    |
//| enough, no memoryallocation is done.                             |
//| It is intended to be used in the performance-critical parts of   |
//| code, where additional allocations can lead to severe performance|
//| degradation                                                      |
//+------------------------------------------------------------------+
static void CTSort::TagSortBuf(double &a[],const int n,int &p1[],int &p2[],
                               CApBuff &buf)
  {
//--- create variables
   int i=0;
   int lv=0;
   int lp=0;
   int rv=0;
   int rp=0;
//--- Special cases
   if(n<=0)
      return;
//--- check
   if(n==1)
     {
      //--- function call
      CApServ::IVectorSetLengthAtLeast(p1,1);
      //--- function call
      CApServ::IVectorSetLengthAtLeast(p2,1);
      p1[0]=0;
      p2[0]=0;
      //--- exit the function
      return;
     }
//--- General case, N>1: prepare permutations table P1
   CApServ::IVectorSetLengthAtLeast(p1,n);
   for(i=0;i<=n-1;i++)
      p1[i]=i;
//--- General case, N>1: sort, update P1
   CApServ::RVectorSetLengthAtLeast(buf.m_ra0,n);
//--- function call
   CApServ::IVectorSetLengthAtLeast(buf.m_ia0,n);
   TagSortFastI(a,p1,buf.m_ra0,buf.m_ia0,n);
//--- General case, N>1: fill permutations table P2
//--- To fill P2 we maintain two arrays:
//--- * PV (Buf.IA0), Position(Value). PV[i] contains position of I-th key at the moment
//--- * VP (Buf.IA1), Value(Position). VP[i] contains key which has position I at the moment
//--- At each step we making permutation of two items:
//--- Left,which is given by position/value pair LP/LV
//--- and Right,which is given by RP/RV
//--- and updating PV[] and VP[] correspondingly.
   CApServ::IVectorSetLengthAtLeast(buf.m_ia0,n);
//--- function call
   CApServ::IVectorSetLengthAtLeast(buf.m_ia1,n);
//--- function call
   CApServ::IVectorSetLengthAtLeast(p2,n);
   for(i=0;i<=n-1;i++)
     {
      buf.m_ia0[i]=i;
      buf.m_ia1[i]=i;
     }
   for(i=0;i<=n-1;i++)
     {
      //--- calculate LP, LV, RP, RV
      lp=i;
      lv=buf.m_ia1[lp];
      rv=p1[i];
      rp=buf.m_ia0[rv];
      //--- Fill P2
      p2[i]=rp;
      //--- update PV and VP
      buf.m_ia1[lp]=rv;
      buf.m_ia1[rp]=lv;
      buf.m_ia0[lv]=rp;
      buf.m_ia0[rv]=lp;
     }
  }
//+------------------------------------------------------------------+
//| Same as TagSort, but optimized for real keys and integer labels. |
//| A is sorted, and same permutations are applied to B.             |
//| NOTES:                                                           |
//| 1.  this function assumes that A[] is finite; it doesn't checks  |
//|     that condition. All other conditions (size of input arrays,  |
//|     etc.) are not checked too.                                   |
//| 2.  this function uses two buffers, BufA and BufB, each is N     |
//|     elements large. They may be preallocated (which will save    |
//|     some time) or not, in which case function will automatically | 
//|     allocate memory.                                             |
//+------------------------------------------------------------------+
static void CTSort::TagSortFastI(double &a[],int &b[],double &bufa[],
                                 int &bufb[],const int n)
  {
//--- create variables
   int    i=0;
   int    j=0;
   bool   isascending;
   bool   isdescending;
   double tmpr=0;
   int    tmpi=0;
//--- Special case
   if(n<=1)
      return;
//--- Test for already sorted set
   isascending=true;
   isdescending=true;
   for(i=1;i<=n-1;i++)
     {
      isascending=isascending && a[i]>=a[i-1];
      isdescending=isdescending && a[i]<=a[i-1];
     }
//--- check
   if(isascending)
      return;
//--- check
   if(isdescending)
     {
      for(i=0;i<=n-1;i++)
        {
         j=n-1-i;
         //--- check
         if(j<=i)
            break;
         //--- swap
         tmpr=a[i];
         a[i]=a[j];
         a[j]=tmpr;
         tmpi=b[i];
         b[i]=b[j];
         b[j]=tmpi;
        }
      //--- exit the function
      return;
     }
//--- General case
   if(CAp::Len(bufa)<n)
      ArrayResizeAL(bufa,n);
//--- check
   if(CAp::Len(bufb)<n)
      ArrayResizeAL(bufb,n);
//--- function call
   TagSortFastIRec(a,b,bufa,bufb,0,n-1);
  }
//+------------------------------------------------------------------+
//| Same as TagSort, but optimized for real keys and real labels.    |
//| A is sorted, and same permutations are applied to B.             |
//| NOTES:                                                           |
//| 1.  this function assumes that A[] is finite; it doesn't checks  | 
//|     etc.) are not that condition. All other conditions (size of  |
//|     input arrays, checked too.                                   |
//| 2.  this function uses two buffers, BufA and BufB, each is N     |
//|     elements large. They may be preallocated (which will save    |
//|     some time) or not, in whichcase function will automatically  |
//|     allocate memory.                                             |
//+------------------------------------------------------------------+
static void CTSort::TagSortFastR(double &a[],double &b[],double &bufa[],
                                 double &bufb[],const int n)
  {
//--- create variables
   int    i=0;
   int    j=0;
   bool   isascending;
   bool   isdescending;
   double tmpr=0;
//--- Special case
   if(n<=1)
      return;
//--- Test for already sorted set
   isascending=true;
   isdescending=true;
   for(i=1;i<=n-1;i++)
     {
      isascending=isascending && a[i]>=a[i-1];
      isdescending=isdescending && a[i]<=a[i-1];
     }
//--- check
   if(isascending)
      return;
//--- check
   if(isdescending)
     {
      for(i=0;i<=n-1;i++)
        {
         j=n-1-i;
         //--- check
         if(j<=i)
            break;
         //--- swap
         tmpr=a[i];
         a[i]=a[j];
         a[j]=tmpr;
         tmpr=b[i];
         b[i]=b[j];
         b[j]=tmpr;
        }
      //--- exit the function
      return;
     }
//--- General case
   if(CAp::Len(bufa)<n)
      ArrayResizeAL(bufa,n);
//--- check
   if(CAp::Len(bufb)<n)
      ArrayResizeAL(bufb,n);
//--- function call
   TagSortFastRRec(a,b,bufa,bufb,0,n-1);
  }
//+------------------------------------------------------------------+
//| Same as TagSort, but optimized for real keys without labels.     |
//| A is sorted, and that's all.                                     |
//| NOTES:                                                           |
//| 1.  this function assumes that A[] is finite; it doesn't checks  |
//|     that condition. All other conditions (size of input arrays,  |
//|     etc.) are not checked too.                                   |
//| 2.  this function uses buffer, BufA, which is N elements large.  |
//|     It may be preallocated (which will save some time) or not,   |
//|     in which casefunction will automatically allocate memory.    |
//+------------------------------------------------------------------+
static void CTSort::TagSortFast(double &a[],double &bufa[],const int n)
  {
//--- Special case
   if(n<=1)
      return;
//--- create variables
   int    i;
   int    j;
   bool   isAsCending;
   bool   isDesCending;
   double tmpr;
//--- Test for already sorted set
   isAsCending=true;
   isDesCending=true;
   for(i=1;i<n;i++)
     {
      isAsCending=isAsCending && a[i]>=a[i-1];
      isDesCending=isDesCending && a[i]<=a[i-1];
     }
//--- check
   if(isAsCending)
      return;
//--- check
   if(isDesCending)
     {
      for(i=0;i<n;i++)
        {
         j=n-1-i;
         if(j<=i)
            break;
         tmpr=a[i];
         a[i]=a[j];
         a[j]=tmpr;
        }
      //--- exit the function
      return;
     }
//--- General case
   if(CAp::Len(bufa)<n)
      ArrayResizeAL(bufa,n);
//--- function call
   TagSortFastRec(a,bufa,0,n-1);
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Heap operations: adds element to the heap                        |
//| PARAMETERS:                                                      |
//|     A       -   heap itself, must be at least array[0..N]        |
//|     B       -   array of integer tags, which are updated         |
//|                 according to permutations in the heap            |
//|     N       -   size of the heap (without new element).          |
//|                 updated on output                                |
//|     VA      -   value of the element being added                 |
//|     VB      -   value of the tag                                 |
//+------------------------------------------------------------------+
static void CTSort::TagHeapPushI(double &a[],int &b[],int &n,const double va,
                                 const int vb)
  {
//--- create variables
   int    j=0;
   int    k=0;
   double v=0;
//--- check
   if(n<0)
      return;
//--- N=0 is a special case
   if(n==0)
     {
      a[0]=va;
      b[0]=vb;
      n=n+1;
      //--- exit the function
      return;
     }
//--- add current point to the heap
//--- (add to the bottom, then move up)
//--- we don't write point to the heap
//--- until its final position is determined
//--- (it allow us to reduce number of array access operations)
   j=n;
   n=n+1;
   while(j>0)
     {
      k=(j-1)/2;
      v=a[k];
      //--- check
      if(v<va)
        {
         //--- swap with higher element
         a[j]=v;
         b[j]=b[k];
         j=k;
        }
      else
        {
         //--- element in its place. terminate.
         break;
        }
     }
//--- change values
   a[j]=va;
   b[j]=vb;
  }
//+------------------------------------------------------------------+
//| Heap operations: replaces top element with new element           |
//| (which is moved down)                                            |
//| PARAMETERS:                                                      |
//|     A       -   heap itself, must be at least array[0..N-1]      |
//|     B       -   array of integer tags, which are updated         |
//|                 according to permutations in the heap            |
//|     N       -   size of the heap                                 |
//|     VA      -   value of the element which replaces top element  |
//|     VB      -   value of the tag                                 |
//+------------------------------------------------------------------+
static void CTSort::TagHeapReplaceTopI(double &a[],int &b[],const int n,
                                       const double va,const int vb)
  {
//--- create variables
   int    j=0;
   int    k1=0;
   int    k2=0;
   double v=0;
   double v1=0;
   double v2=0;
//--- check
   if(n<1)
      return;
//--- N=1 is a special case
   if(n==1)
     {
      a[0]=va;
      b[0]=vb;
      //--- exit the function
      return;
     }
//--- move down through heap:
//--- * J  -   current element
//--- * K1 -   first child (always exists)
//--- * K2 -   second child (may not exists)
//--- we don't write point to the heap
//--- until its final position is determined
//--- (it allow us to reduce number of array access operations)
   j=0;
   k1=1;
   k2=2;
   while(k1<n)
     {
      //--- check
      if(k2>=n)
        {
         //--- only one child.
         //--- swap and terminate (because this child
         //--- have no siblings due to heap structure)
         v=a[k1];
         //--- check
         if(v>va)
           {
            a[j]=v;
            b[j]=b[k1];
            j=k1;
           }
         break;
        }
      else
        {
         //--- two childs
         v1=a[k1];
         v2=a[k2];
         //--- check
         if(v1>v2)
           {
            //--- check
            if(va<v1)
              {
               a[j]=v1;
               b[j]=b[k1];
               j=k1;
              }
            else
               break;
           }
         else
           {
            //--- check
            if(va<v2)
              {
               a[j]=v2;
               b[j]=b[k2];
               j=k2;
              }
            else
               break;
           }
         //--- change values
         k1=2*j+1;
         k2=2*j+2;
        }
     }
//--- change values
   a[j]=va;
   b[j]=vb;
  }
//+------------------------------------------------------------------+
//| Heap operations: pops top element from the heap                  |
//| PARAMETERS:                                                      |
//|     A       -   heap itself, must be at least array[0..N-1]      |
//|     B       -   array of integer tags, which are updated         |
//|                 according to permutations in the heap            |
//|     N       -   size of the heap, N>=1                           |
//| On output top element is moved to A[N-1], B[N-1], heap is        |
//| reordered, N is decreased by 1.                                  |
//+------------------------------------------------------------------+
static void CTSort::TagHeapPopI(double &a[],int &b[],int &n)
  {
//--- create variables
   double va=0;
   int    vb=0;
//--- check
   if(n<1)
      return;
//--- N=1 is a special case
   if(n==1)
     {
      n=0;
      return;
     }
//--- swap top element and last element,
//--- then reorder heap
   va=a[n-1];
   vb=b[n-1];
   a[n-1]=a[0];
   b[n-1]=b[0];
   n=n-1;
//--- function call
   TagHeapReplaceTopI(a,b,n,va,vb);
  }
//+------------------------------------------------------------------+
//| Internal TagSortFastI: sorts A[I1...I2] (both bounds are         |
//| included), applies same permutations to B.                       |
//+------------------------------------------------------------------+
static void CTSort::TagSortFastIRec(double &a[],int &b[],double &bufa[],
                                    int &bufb[],const int i1,const int i2)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    cntless=0;
   int    cnteq=0;
   int    cntgreater=0;
   double tmpr=0;
   int    tmpi=0;
   double v0=0;
   double v1=0;
   double v2=0;
   double vp=0;
//--- Fast exit
   if(i2<=i1)
      return;
//--- Non-recursive sort for small arrays
   if(i2-i1<=16)
     {
      for(j=i1+1;j<=i2;j++)
        {
         //--- Search elements [I1..J-1] for place to insert Jth element.
         //--- This code stops immediately if we can leave A[J] at J-th position
         //--- (all elements have same value of A[J] larger than any of them)
         tmpr=a[j];
         tmpi=j;
         for(k=j-1;k>=i1;k--)
           {
            //--- check
            if(a[k]<=tmpr)
               break;
            tmpi=k;
           }
         k=tmpi;
         //--- Insert Jth element into Kth position
         if(k!=j)
           {
            //--- change values
            tmpr=a[j];
            tmpi=b[j];
            for(i=j-1;i>=k;i--)
              {
               a[i+1]=a[i];
               b[i+1]=b[i];
              }
            a[k]=tmpr;
            b[k]=tmpi;
           }
        }
      //--- exit the function
      return;
     }
//--- Quicksort: choose pivot
//--- Here we assume that I2-I1>=2
   v0=a[i1];
   v1=a[i1+(i2-i1)/2];
   v2=a[i2];
//--- check
   if(v0>v1)
     {
      tmpr=v1;
      v1=v0;
      v0=tmpr;
     }
//--- check
   if(v1>v2)
     {
      tmpr=v2;
      v2=v1;
      v1=tmpr;
     }
//--- check
   if(v0>v1)
     {
      tmpr=v1;
      v1=v0;
      v0=tmpr;
     }
   vp=v1;
//--- now pass through A/B and:
//--- * move elements that are LESS than VP to the left of A/B
//--- * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
//--- * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
//--- * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
//--- * move elements from the left of BufA/BufB to the end of A/B
   cntless=0;
   cnteq=0;
   cntgreater=0;
   for(i=i1;i<=i2;i++)
     {
      v0=a[i];
      //--- check
      if(v0<vp)
        {
         //--- LESS
         k=i1+cntless;
         //--- check
         if(i!=k)
           {
            a[k]=v0;
            b[k]=b[i];
           }
         cntless=cntless+1;
         continue;
        }
      //--- check
      if(v0==vp)
        {
         //--- EQUAL
         k=i2-cnteq;
         bufa[k]=v0;
         bufb[k]=b[i];
         cnteq=cnteq+1;
         continue;
        }
      //--- GREATER
      k=i1+cntgreater;
      bufa[k]=v0;
      bufb[k]=b[i];
      cntgreater=cntgreater+1;
     }
//--- change values
   for(i=0;i<=cnteq-1;i++)
     {
      j=i1+cntless+cnteq-1-i;
      k=i2+i-(cnteq-1);
      a[j]=bufa[k];
      b[j]=bufb[k];
     }
//--- change values
   for(i=0;i<=cntgreater-1;i++)
     {
      j=i1+cntless+cnteq+i;
      k=i1+i;
      a[j]=bufa[k];
      b[j]=bufb[k];
     }
//--- Sort left and right parts of the array (ignoring middle part)
   TagSortFastIRec(a,b,bufa,bufb,i1,i1+cntless-1);
//--- function call
   TagSortFastIRec(a,b,bufa,bufb,i1+cntless+cnteq,i2);
  }
//+------------------------------------------------------------------+
//| Internal TagSortFastR: sorts A[I1...I2] (both bounds are         |
//| included), applies same permutations to B.                       |
//+------------------------------------------------------------------+
static void CTSort::TagSortFastRRec(double &a[],double &b[],double &bufa[],
                                    double &bufb[],const int i1,const int i2)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   double tmpr=0;
   double tmpr2=0;
   int    tmpi=0;
   int    cntless=0;
   int    cnteq=0;
   int    cntgreater=0;
   double v0=0;
   double v1=0;
   double v2=0;
   double vp=0;
//--- Fast exit
   if(i2<=i1)
      return;
//--- Non-recursive sort for small arrays
   if(i2-i1<=16)
     {
      for(j=i1+1;j<=i2;j++)
        {
         //--- Search elements [I1..J-1] for place to insert Jth element.
         //--- This code stops immediatly if we can leave A[J] at J-th position
         //--- (all elements have same value of A[J] larger than any of them)
         tmpr=a[j];
         tmpi=j;
         for(k=j-1;k>=i1;k--)
           {
            //--- check
            if(a[k]<=tmpr)
               break;
            tmpi=k;
           }
         k=tmpi;
         //--- Insert Jth element into Kth position
         if(k!=j)
           {
            //--- change values
            tmpr=a[j];
            tmpr2=b[j];
            for(i=j-1;i>=k;i--)
              {
               a[i+1]=a[i];
               b[i+1]=b[i];
              }
            a[k]=tmpr;
            b[k]=tmpr2;
           }
        }
      //--- exit the function
      return;
     }
//--- Quicksort: choose pivot
//--- Here we assume that I2-I1>=16
   v0=a[i1];
   v1=a[i1+(i2-i1)/2];
   v2=a[i2];
//--- check
   if(v0>v1)
     {
      tmpr=v1;
      v1=v0;
      v0=tmpr;
     }
//--- check
   if(v1>v2)
     {
      tmpr=v2;
      v2=v1;
      v1=tmpr;
     }
//--- check
   if(v0>v1)
     {
      tmpr=v1;
      v1=v0;
      v0=tmpr;
     }
   vp=v1;
//--- now pass through A/B and:
//--- * move elements that are LESS than VP to the left of A/B
//--- * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
//--- * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
//--- * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
//--- * move elements from the left of BufA/BufB to the end of A/B
   cntless=0;
   cnteq=0;
   cntgreater=0;
   for(i=i1;i<=i2;i++)
     {
      v0=a[i];
      //--- check
      if(v0<vp)
        {
         //--- LESS
         k=i1+cntless;
         //--- check
         if(i!=k)
           {
            a[k]=v0;
            b[k]=b[i];
           }
         cntless=cntless+1;
         continue;
        }
      //--- check
      if(v0==vp)
        {
         //--- EQUAL
         k=i2-cnteq;
         bufa[k]=v0;
         bufb[k]=b[i];
         cnteq=cnteq+1;
         continue;
        }
      //--- GREATER
      k=i1+cntgreater;
      bufa[k]=v0;
      bufb[k]=b[i];
      cntgreater=cntgreater+1;
     }
//--- change values
   for(i=0;i<=cnteq-1;i++)
     {
      j=i1+cntless+cnteq-1-i;
      k=i2+i-(cnteq-1);
      a[j]=bufa[k];
      b[j]=bufb[k];
     }
//--- change values
   for(i=0;i<=cntgreater-1;i++)
     {
      j=i1+cntless+cnteq+i;
      k=i1+i;
      a[j]=bufa[k];
      b[j]=bufb[k];
     }
//--- Sort left and right parts of the array (ignoring middle part)
   TagSortFastRRec(a,b,bufa,bufb,i1,i1+cntless-1);
//--- function call
   TagSortFastRRec(a,b,bufa,bufb,i1+cntless+cnteq,i2);
  }
//+------------------------------------------------------------------+
//| Internal TagSortFastI: sorts A[I1...I2] (both bounds are         |
//| included), applies same permutations to B.                       |
//+------------------------------------------------------------------+
static void CTSort::TagSortFastRec(double &a[],double &bufa[],
                                   const int i1,const int i2)
  {
//--- check
   if(i2<=i1)
      return;
//--- create variables
   int    cntLess;
   int    cntEq;
   int    cntGreat;
   int    i;
   int    j;
   int    k;
   double tmpr;
   int    tmpi;
   double v0;
   double v1;
   double v2;
   double vp;
//--- Non-recursive sort for small arrays
   if(i2-i1<=16)
     {
      for(j=i1+1;j<=i2;j++)
        {
         //--- Search elements [I1..J-1] for place to insert Jth element.
         //--- This code stops immediatly if we can leave A[J] at J-th position
         //--- (all elements have same value of A[J] larger than any of them)
         tmpr=a[j];
         tmpi=j;
         for(k=j-1;k>=i1;k--)
           {
            //--- check
            if(a[k]<=tmpr)
               break;
            tmpi=k;
           }
         k=tmpi;
         //--- Insert Jth element into Kth position
         if(k!=j)
           {
            tmpr=a[j];
            for(i=j-1;i>=k;i--)
               a[i+1]=a[i];
            a[k]=tmpr;
           }
        }
      return;
     }
//--- Quicksort: choose pivot
//--- Here we assume that I2-I1>=16
   v0=a[i1];
   v1=a[i1+(i2-i1)/2];
   v2=a[i2];
//--- check
   if(v0>v1)
     {
      tmpr=v1;
      v1=v0;
      v0=tmpr;
     }
//--- check
   if(v1>v2)
     {
      tmpr=v2;
      v2=v1;
      v1=tmpr;
     }
//--- check
   if(v0>v1)
     {
      tmpr=v1;
      v1=v0;
      v0=tmpr;
     }
   vp=v1;
//--- now pass through A/B and:
//--- * move elements that are LESS than VP to the left of A/B
//--- * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
//--- * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
//--- * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
//--- * move elements from the left of BufA/BufB to the end of A/B
   cntLess=0;
   cntEq=0;
   cntGreat=0;
   for(i=i1;i<=i2;i++)
     {
      v0=a[i];
      //--- check
      if(v0<vp)
        {
         //--- LESS
         k=i1+cntLess;
         if(i!=k)
            a[k]=v0;
         cntLess++;
         continue;
        }
      //--- check
      if(v0==vp)
        {
         //--- EQUAL
         k=i2-cntEq;
         bufa[k]=v0;
         cntEq++;
         continue;
        }
      //--- GREATER
      k=i1+cntGreat;
      bufa[k]=v0;
      cntGreat++;
     }
//--- change values
   for(i=0;i<cntEq;i++)
     {
      j=i1+cntLess+cntEq-1-i;
      k=i2+i-(cntEq-1);
      a[j]=bufa[k];
     }
   for(i=0;i<cntGreat;i++)
     {
      j=i1+cntLess+cntEq+i;
      k=i1+i;
      a[j]=bufa[k];
     }
//--- Sort left and right parts of the array (ignoring middle part)
   TagSortFastRec(a,bufa,i1,i1+cntLess-1);
   TagSortFastRec(a,bufa,i1+cntLess+cntEq,i2);
   return;
  }
//+------------------------------------------------------------------+
//| Calculation of basic statistical properties                      |
//+------------------------------------------------------------------+
class CBasicStatOps
  {
public:
   //--- constructor, destructor
                     CBasicStatOps(void);
                    ~CBasicStatOps(void);
   //--- method
   static void       RankX(double &x[],const int n,CApBuff &buf);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CBasicStatOps::CBasicStatOps(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CBasicStatOps::~CBasicStatOps(void)
  {

  }
//+------------------------------------------------------------------+
//| Internal ranking subroutine                                      |
//+------------------------------------------------------------------+
static void CBasicStatOps::RankX(double &x[],const int n,CApBuff &buf)
  {
//--- Prepare
   if(n<1)
      return;
//--- check
   if(n==1)
     {
      x[0]=1;
      return;
     }
//--- create variables
   int    i;
   int    j;
   int    k;
   int    t;
   double tmp;
   int    tmpi;
//--- check
   if(CAp::Len(buf.m_ra1)<n)
      ArrayResizeAL(buf.m_ra1,n);
//--- check
   if(CAp::Len(buf.m_ia1)<n)
      ArrayResizeAL(buf.m_ia1,n);
//--- copy
   for(i=0;i<n;i++)
     {
      buf.m_ra1[i]=x[i];
      buf.m_ia1[i]=i;
     }
//--- sort {R, C}
   if(n!=1)
     {
      i=2;
      //--- cycle
      do
        {
         t=i;
         //--- cycle
         while(t!=1)
           {
            k=t/2;
            //--- check
            if(buf.m_ra1[k-1]>=buf.m_ra1[t-1])
               t=1;
            else
              {
               //--- swap
               tmp=buf.m_ra1[k-1];
               buf.m_ra1[k-1]=buf.m_ra1[t-1];
               buf.m_ra1[t-1]=tmp;
               tmpi=buf.m_ia1[k-1];
               buf.m_ia1[k-1]=buf.m_ia1[t-1];
               buf.m_ia1[t-1]=tmpi;
               //--- set value
               t=k;
              }
           }
         //--- next iteration
         i=i+1;
        }
      while(i<=n);
      //--- set value
      i=n-1;
      //--- cycle
      do
        {
         //--- swap
         tmp=buf.m_ra1[i];
         buf.m_ra1[i]=buf.m_ra1[0];
         buf.m_ra1[0]=tmp;
         tmpi=buf.m_ia1[i];
         buf.m_ia1[i]=buf.m_ia1[0];
         buf.m_ia1[0]=tmpi;
         //--- set value
         t=1;
         while(t!=0)
           {
            k=2*t;
            if(k>i)
               t=0;
            else
              {
               //--- check
               if(k<i)
                  if(buf.m_ra1[k]>buf.m_ra1[k-1])
                     k++;
               //--- check
               if(buf.m_ra1[t-1]>=buf.m_ra1[k-1])
                  t=0;
               else
                 {
                  //--- swap
                  tmp=buf.m_ra1[k-1];
                  buf.m_ra1[k-1]=buf.m_ra1[t-1];
                  buf.m_ra1[t-1]=tmp;
                  tmpi=buf.m_ia1[k-1];
                  buf.m_ia1[k-1]=buf.m_ia1[t-1];
                  buf.m_ia1[t-1]=tmpi;
                  //--- set value
                  t=k;
                 }
              }
           }
         i=i-1;
        }
      while(i>=1);
     }
//--- compute tied ranks
   i=0;
   while(i<n)
     {
      j=i+1;
      while(j<n)
        {
         //--- check
         if(buf.m_ra1[j]!=buf.m_ra1[i])
            break;
         j=j+1;
        }
      for(k=i;k<j;k++)
         buf.m_ra1[k]=1+(double)(i+j-1)/2.0;
      //--- set value
      i=j;
     }
//--- back to x
   for(i=0;i<n;i++)
      x[buf.m_ia1[i]]=buf.m_ra1[i];
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Class includes templates for future functions                    |
//+------------------------------------------------------------------+
class CAblasF
  {
public:
   //--- constructor, destructor
                     CAblasF(void);
                    ~CAblasF(void);
   //--- return false
   static bool       CMatrixRank1F(const int m,const int n,CMatrixComplex &a,int ia,int ja,al_complex &u[],int iu,al_complex &v[],int iv);
   static bool       RMatrixRank1F(const int m,const int n,CMatrixDouble &a,int ia,int ja,double &u[],int iu,double &v[],int iv);
   static bool       CMatrixMVF(const int m,const int n,CMatrixComplex &a,int ia,int ja,int opa,al_complex &x[],int ix,al_complex &y[],int iy);
   static bool       RMatrixMVF(const int m,const int n,CMatrixDouble &a,int ia,int ja,int opa,double &x[],int ix,double &y[],int iy);
   static bool       CMatrixRightTrsMF(const int m,const int n,CMatrixComplex &a,const int i1,int j1,const bool isupper,bool isunit,int optype,CMatrixComplex &x,int i2,int j2);
   static bool       CMatrixLeftTrsMF(const int m,const int n,CMatrixComplex &a,const int i1,int j1,const bool isupper,bool isunit,int optype,CMatrixComplex &x,int i2,int j2);
   static bool       RMatrixRightTrsMF(const int m,const int n,CMatrixDouble &a,const int i1,int j1,const bool isupper,bool isunit,int optype,CMatrixDouble &x,int i2,int j2);
   static bool       RMatrixLeftTrsMF(const int m,const int n,CMatrixDouble &a,const int i1,int j1,const bool isupper,bool isunit,int optype,CMatrixDouble &x,int i2,int j2);
   static bool       CMatrixSyrkF(const int n,int k,double alpha,CMatrixComplex &a,int ia,int ja,int optypea,double beta,CMatrixComplex &c,int ic,int jc,bool isupper);
   static bool       RMatrixSyrkF(const int n,int k,double alpha,CMatrixDouble &a,int ia,int ja,int optypea,double beta,CMatrixDouble &c,int ic,int jc,bool isupper);
   static bool       RMatrixGemmF(const int m,const int n,int k,double alpha,CMatrixDouble &a,int ia,int ja,int optypea,CMatrixDouble &b,int ib,int jb,int optypeb,double beta,CMatrixDouble &c,int ic,int jc);
   static bool       CMatrixGemmF(const int m,const int n,int k,al_complex &alpha,CMatrixComplex &a,int ia,int ja,int optypea,CMatrixComplex &b,int ib,int jb,int optypeb,al_complex &beta,CMatrixComplex &c,int ic,int jc);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CAblasF::CAblasF(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CAblasF::~CAblasF(void)
  {

  }
//+------------------------------------------------------------------+
//| Fast kernel                                                      |
//+------------------------------------------------------------------+
static bool CAblasF::CMatrixRank1F(const int m,const int n,CMatrixComplex &a,
                                   int ia,int ja,al_complex &u[],int iu,
                                   al_complex &v[],int iv)
  {
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Fast kernel                                                      |
//+------------------------------------------------------------------+
static bool CAblasF::RMatrixRank1F(const int m,const int n,CMatrixDouble &a,
                                   int ia,int ja,double &u[],int iu,
                                   double &v[],int iv)
  {
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Fast kernel                                                      |
//+------------------------------------------------------------------+
static bool CAblasF::CMatrixMVF(const int m,const int n,CMatrixComplex &a,
                                int ia,int ja,int opa,al_complex &x[],int ix,
                                al_complex &y[],int iy)
  {
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Fast kernel                                                      |
//+------------------------------------------------------------------+
static bool CAblasF::RMatrixMVF(const int m,const int n,CMatrixDouble &a,
                                int ia,int ja,int opa,double &x[],int ix,
                                double &y[],int iy)
  {
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Fast kernel                                                      |
//+------------------------------------------------------------------+
static bool CAblasF::CMatrixRightTrsMF(const int m,const int n,CMatrixComplex &a,
                                       const int i1,int j1,const bool isupper,
                                       bool isunit,int optype,CMatrixComplex &x,
                                       int i2,int j2)
  {
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Fast kernel                                                      |
//+------------------------------------------------------------------+
static bool CAblasF::CMatrixLeftTrsMF(const int m,const int n,CMatrixComplex &a,
                                      const int i1,int j1,const bool isupper,
                                      bool isunit,int optype,CMatrixComplex &x,
                                      int i2,int j2)
  {
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Fast kernel                                                      |
//+------------------------------------------------------------------+
static bool CAblasF::RMatrixRightTrsMF(const int m,const int n,CMatrixDouble &a,
                                       const int i1,int j1,const bool isupper,
                                       bool isunit,int optype,CMatrixDouble &x,
                                       int i2,int j2)
  {
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Fast kernel                                                      |
//+------------------------------------------------------------------+
static bool CAblasF::RMatrixLeftTrsMF(const int m,const int n,CMatrixDouble &a,
                                      const int i1,int j1,const bool isupper,
                                      bool isunit,int optype,CMatrixDouble &x,
                                      int i2,int j2)
  {
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Fast kernel                                                      |
//+------------------------------------------------------------------+
static bool CAblasF::CMatrixSyrkF(const int n,int k,double alpha,CMatrixComplex &a,
                                  int ia,int ja,int optypea,double beta,
                                  CMatrixComplex &c,int ic,int jc,bool isupper)
  {
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Fast kernel                                                      |
//+------------------------------------------------------------------+
static bool CAblasF::RMatrixSyrkF(const int n,int k,double alpha,CMatrixDouble &a,
                                  int ia,int ja,int optypea,double beta,
                                  CMatrixDouble &c,int ic,int jc,bool isupper)
  {
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Fast kernel                                                      |
//+------------------------------------------------------------------+
static bool CAblasF::RMatrixGemmF(const int m,const int n,int k,double alpha,
                                  CMatrixDouble &a,int ia,int ja,int optypea,
                                  CMatrixDouble &b,int ib,int jb,int optypeb,
                                  double beta,CMatrixDouble &c,int ic,int jc)
  {
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Fast kernel                                                      |
//+------------------------------------------------------------------+
static bool CAblasF::CMatrixGemmF(const int m,const int n,int k,al_complex &alpha,
                                  CMatrixComplex &a,int ia,int ja,int optypea,
                                  CMatrixComplex &b,int ib,int jb,int optypeb,
                                  al_complex &beta,CMatrixComplex &c,int ic,int jc)
  {
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Computing matrix-vector and matrix-matrix                        |
//+------------------------------------------------------------------+
class CBlas
  {
public:
   //--- constructor, destructor
                     CBlas(void);
                    ~CBlas(void);
   //--- methods
   static double     VectorNorm2(double &x[],const int i1,const int i2);
   static int        VectorIdxAbsMax(double &x[],const int i1,const int i2);
   static int        ColumnIdxAbsMax(CMatrixDouble &x,const int i1,const int i2,const int j);
   static int        RowIdxAbsMax(CMatrixDouble &x,const int j1,const int j2,const int i);
   static double     UpperHessenberg1Norm(CMatrixDouble &a,const int i1,const int i2,const int j1,const int j2,double &work[]);
   static void       CopyMatrix(CMatrixDouble &a,const int is1,const int is2,const int js1,const int js2,CMatrixDouble &b,const int id1,const int id2,const int jd1,const int jd2);
   static void       InplaceTranspose(CMatrixDouble &a,const int i1,const int i2,const int j1,const int j2,double &work[]);
   static void       CopyAndTranspose(CMatrixDouble &a,const int is1,const int is2,const int js1,const int js2,CMatrixDouble &b,const int id1,const int id2,const int jd1,const int jd2);
   static void       MatrixVectorMultiply(CMatrixDouble &a,const int i1,const int i2,const int j1,const int j2,const bool trans,double &x[],const int ix1,const int ix2,const double alpha,double &y[],const int iy1,const int iy2,const double beta);
   static double     PyThag2(double x,double y);
   static void       MatrixMatrixMultiply(CMatrixDouble &a,const int ai1,const int ai2,const int aj1,const int aj2,const bool transa,CMatrixDouble &b,const int bi1,const int bi2,const int bj1,const int bj2,const bool transb,const double alpha,CMatrixDouble &c,const int ci1,const int ci2,const int cj1,const int cj2,const double beta,double &work[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CBlas::CBlas(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CBlas::~CBlas(void)
  {

  }
//+------------------------------------------------------------------+
//| Vector norm                                                      |
//+------------------------------------------------------------------+
static double CBlas::VectorNorm2(double &x[],const int i1,const int i2)
  {
//--- create variables
   int    n=i2-i1+1;
   int    ix=0;
   double absxi=0;
   double scl=0;
   double ssq=1;
//--- check
   if(n<1)
      return(0);
//--- check
   if(n==1)
      return(MathAbs(x[i1]));
//--- norm
   for(ix=i1;ix<=i2;ix++)
     {
      //--- check
      if(x[ix]!=0.0)
        {
         absxi=MathAbs(x[ix]);
         //--- check
         if(scl<absxi)
           {
            ssq=1+ssq*CMath::Sqr(scl/absxi);
            scl=absxi;
           }
         else
            ssq=ssq+CMath::Sqr(absxi/scl);
        }
     }
//--- return result
   return(scl*MathSqrt(ssq));
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static int CBlas::VectorIdxAbsMax(double &x[],const int i1,const int i2)
  {
//--- create variables
   int    i=0;
   double a=0;
   int    result=i1;
//--- calculation
   a=MathAbs(x[result]);
   for(i=i1+1;i<=i2;i++)
     {
      //--- check
      if(MathAbs(x[i])>MathAbs(x[result]))
         result=i;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static int CBlas::ColumnIdxAbsMax(CMatrixDouble &x,const int i1,const int i2,const int j)
  {
//--- create variables
   int result=i1;
   int i=0;
   double a=0;
//--- calculation
   a=MathAbs(x[result][j]);
   for(i=i1+1;i<=i2;i++)
     {
      //--- check
      if(MathAbs(x[i][j])>MathAbs(x[result][j]))
         result=i;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static int CBlas::RowIdxAbsMax(CMatrixDouble &x,const int j1,const int j2,const int i)
  {
//--- create variables
   int    result=j1;
   int    j=0;
   double a=0;
//--- calculation
   a=MathAbs(x[i][result]);
   for(j=j1+1;j<=j2;j++)
     {
      //--- check
      if(MathAbs(x[i][j])>MathAbs(x[i][result]))
         result=j;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Upper Hessenberg norm                                            |
//+------------------------------------------------------------------+
static double CBlas::UpperHessenberg1Norm(CMatrixDouble &a,const int i1,
                                          const int i2,const int j1,
                                          const int j2,double &work[])
  {
//--- create variables
   double result=0;
   int    i=0;
   int    j=0;
//--- check
   if(!CAp::Assert(i2-i1==j2-j1,__FUNCTION__+": I2-I1!=J2-J1!"))
      return(EMPTY_VALUE);
   for(j=j1;j<=j2;j++)
      work[j]=0;
   for(i=i1;i<=i2;i++)
     {
      for(j=MathMax(j1,j1+i-i1-1);j<=j2;j++)
         work[j]=work[j]+MathAbs(a[i][j]);
     }
//--- get result
   result=0;
   for(j=j1;j<=j2;j++)
      result=MathMax(result,work[j]);
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Copy matrix                                                      |
//+------------------------------------------------------------------+
static void CBlas::CopyMatrix(CMatrixDouble &a,const int is1,const int is2,
                              const int js1,const int js2,CMatrixDouble &b,
                              const int id1,const int id2,const int jd1,const int jd2)
  {
//--- create variables
   int isrc=0;
   int idst=0;
   int i_=0;
   int i1_=0;
//--- check
   if(is1>is2 || js1>js2)
      return;
//--- check
   if(!CAp::Assert(is2-is1==id2-id1,__FUNCTION__+": different sizes!"))
      return;
//--- check
   if(!CAp::Assert(js2-js1==jd2-jd1,__FUNCTION__+": different sizes!"))
      return;
//--- copy
   for(isrc=is1;isrc<=is2;isrc++)
     {
      idst=isrc-is1+id1;
      i1_=js1-jd1;
      for(i_=jd1;i_<=jd2;i_++)
         b[idst].Set(i_,a[isrc][i_+i1_]);
     }
  }
//+------------------------------------------------------------------+
//| Matrix transpose                                                 |
//+------------------------------------------------------------------+
static void CBlas::InplaceTranspose(CMatrixDouble &a,const int i1,const int i2,
                                    const int j1,const int j2,double &work[])
  {
//--- create variables
   int i=0;
   int j=0;
   int ips=0;
   int jps=0;
   int l=0;
   int i_=0;
   int i1_=0;
//--- check
   if(i1>i2 || j1>j2)
      return;
//--- check
   if(!CAp::Assert(i1-i2==j1-j2,__FUNCTION__+": incorrect array size!"))
      return;
   for(i=i1;i<=i2-1;i++)
     {
      //--- change values
      j=j1+i-i1;
      ips=i+1;
      jps=j1+ips-i1;
      l=i2-i;
      i1_=ips-1;
      //--- transpose
      for(i_=1;i_<=l;i_++)
         work[i_]=a[i_+i1_][j];
      i1_=jps-ips;
      for(i_=ips;i_<=i2;i_++)
         a[i_].Set(j,a[i][i_+i1_]);
      i1_=1-jps;
      for(i_=jps;i_<=j2;i_++)
         a[i].Set(i_,work[i_+i1_]);
     }
  }
//+------------------------------------------------------------------+
//| Copy and transpose matrix                                        |
//+------------------------------------------------------------------+
static void CBlas::CopyAndTranspose(CMatrixDouble &a,const int is1,const int is2,
                                    const int js1,const int js2,CMatrixDouble &b,
                                    const int id1,const int id2,const int jd1,const int jd2)
  {
//--- create variables
   int isrc=0;
   int jdst=0;
   int i_=0;
   int i1_=0;
//--- check
   if(is1>is2 || js1>js2)
      return;
//--- check
   if(!CAp::Assert(is2-is1==jd2-jd1,__FUNCTION__+": different sizes!"))
      return;
//--- check
   if(!CAp::Assert(js2-js1==id2-id1,__FUNCTION__+": different sizes!"))
      return;
//--- copy and transpose
   for(isrc=is1;isrc<=is2;isrc++)
     {
      jdst=isrc-is1+jd1;
      i1_=js1-id1;
      for(i_=id1;i_<=id2;i_++)
         b[i_].Set(jdst,a[isrc][i_+i1_]);
     }
  }
//+------------------------------------------------------------------+
//| Matrix vector multiply                                           |
//+------------------------------------------------------------------+
static void CBlas::MatrixVectorMultiply(CMatrixDouble &a,const int i1,const int i2,
                                        const int j1,const int j2,const bool trans,
                                        double &x[],const int ix1,const int ix2,
                                        const double alpha,double &y[],const int iy1,
                                        const int iy2,const double beta)
  {
//--- create variables
   int    i=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(!trans)
     {
      //--- y := alpha*A*x + beta*y;
      if(i1>i2 || j1>j2)
         return;
      //--- check
      if(!CAp::Assert(j2-j1==ix2-ix1,__FUNCTION__+": A and X dont match!"))
         return;
      //--- check
      if(!CAp::Assert(i2-i1==iy2-iy1,__FUNCTION__+": A and Y dont match!"))
         return;
      //--- beta*y
      if(beta==0.0)
        {
         for(i=iy1;i<=iy2;i++)
            y[i]=0;
        }
      else
        {
         for(i_=iy1;i_<=iy2;i_++)
            y[i_]=beta*y[i_];
        }
      //--- alpha*A*x
      for(i=i1;i<=i2;i++)
        {
         i1_=ix1-j1;
         v=0.0;
         for(i_=j1;i_<=j2;i_++)
            v+=a[i][i_]*x[i_+i1_];
         y[iy1+i-i1]=y[iy1+i-i1]+alpha*v;
        }
     }
   else
     {
      //--- y := alpha*A'*x + beta*y;
      if(i1>i2 || j1>j2)
         return;
      //--- check
      if(!CAp::Assert(i2-i1==ix2-ix1,"MatrixVectorMultiply: A and X dont match!"))
         return;
      //--- check
      if(!CAp::Assert(j2-j1==iy2-iy1,"MatrixVectorMultiply: A and Y dont match!"))
         return;
      //--- beta*y
      if(beta==0.0)
        {
         for(i=iy1;i<=iy2;i++)
            y[i]=0;
        }
      else
        {
         for(i_=iy1;i_<=iy2;i_++)
            y[i_]=beta*y[i_];
        }
      //--- alpha*A'*x
      for(i=i1;i<=i2;i++)
        {
         v=alpha*x[ix1+i-i1];
         i1_=j1-iy1;
         for(i_=iy1;i_<=iy2;i_++)
            y[i_]=y[i_]+v*a[i][i_+i1_];
        }
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static double CBlas::PyThag2(double x,double y)
  {
//--- create variables
   double result=0;
   double w=0;
   double xabs=0;
   double yabs=0;
   double z=0;
//--- initialization
   xabs=MathAbs(x);
   yabs=MathAbs(y);
   w=MathMax(xabs,yabs);
   z=MathMin(xabs,yabs);
//--- check
   if(z==0.0)
      result=w;
   else
      result=w*MathSqrt(1+CMath::Sqr(z/w));
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Matrix matrix multiply                                           |
//+------------------------------------------------------------------+
static void CBlas::MatrixMatrixMultiply(CMatrixDouble &a,const int ai1,const int ai2,
                                        const int aj1,const int aj2,const bool transa,
                                        CMatrixDouble &b,const int bi1,const int bi2,
                                        const int bj1,const int bj2,const bool transb,
                                        const double alpha,CMatrixDouble &c,const int ci1,
                                        const int ci2,const int cj1,const int cj2,
                                        const double beta,double &work[])
  {
//--- create variables
   int    arows=0;
   int    acols=0;
   int    brows=0;
   int    bcols=0;
   int    crows=0;
   int    ccols=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    l=0;
   int    r=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- Setup
   if(!transa)
     {
      arows=ai2-ai1+1;
      acols=aj2-aj1+1;
     }
   else
     {
      arows=aj2-aj1+1;
      acols=ai2-ai1+1;
     }
//--- check
   if(!transb)
     {
      brows=bi2-bi1+1;
      bcols=bj2-bj1+1;
     }
   else
     {
      brows=bj2-bj1+1;
      bcols=bi2-bi1+1;
     }
//--- check
   if(!CAp::Assert(acols==brows,__FUNCTION__+": incorrect matrix sizes!"))
      return;
//--- check
   if(arows<=0 || acols<=0 || brows<=0 || bcols<=0)
      return;
   crows=arows;
   ccols=bcols;
//--- Test WORK
   i=MathMax(arows,acols);
   i=MathMax(brows,i);
   i=MathMax(i,bcols);
   work[1]=0;
   work[i]=0;
//--- Prepare C
   if(beta==0.0)
     {
      for(i=ci1;i<=ci2;i++)
         for(j=cj1;j<=cj2;j++)
            c[i].Set(j,0);
     }
   else
     {
      for(i=ci1;i<=ci2;i++)
         for(i_=cj1;i_<=cj2;i_++)
            c[i].Set(i_,beta*c[i][i_]);
     }
//--- A*B
   if(!transa && !transb)
     {
      for(l=ai1;l<=ai2;l++)
        {
         for(r=bi1;r<=bi2;r++)
           {
            //--- change values
            v=alpha*a[l][aj1+r-bi1];
            k=ci1+l-ai1;
            i1_=bj1-cj1;
            for(i_=cj1;i_<=cj2;i_++)
               c[k].Set(i_,c[k][i_]+v*b[r][i_+i1_]);
           }
        }
      //--- exit the function
      return;
     }
//--- A*B'
   if(!transa && transb)
     {
      //--- check
      if(arows*acols<brows*bcols)
        {
         for(r=bi1;r<=bi2;r++)
           {
            for(l=ai1;l<=ai2;l++)
              {
               //--- change values
               i1_=bj1-aj1;
               v=0.0;
               for(i_=aj1;i_<=aj2;i_++)
                  v+=a[l][i_]*b[r][i_+i1_];
               c[ci1+l-ai1].Set(cj1+r-bi1,c[ci1+l-ai1][cj1+r-bi1]+alpha*v);
              }
           }
         //--- exit the function
         return;
        }
      else
        {
         for(l=ai1;l<=ai2;l++)
           {
            for(r=bi1;r<=bi2;r++)
              {
               //--- change values
               i1_=bj1-aj1;
               v=0.0;
               for(i_=aj1;i_<=aj2;i_++)
                  v+=a[l][i_]*b[r][i_+i1_];
               c[ci1+l-ai1].Set(cj1+r-bi1,c[ci1+l-ai1][cj1+r-bi1]+alpha*v);
              }
           }
         //--- exit the function
         return;
        }
     }
//--- A'*B
   if(transa && !transb)
     {
      for(l=aj1;l<=aj2;l++)
        {
         for(r=bi1;r<=bi2;r++)
           {
            //--- change values
            v=alpha*a[ai1+r-bi1][l];
            k=ci1+l-aj1;
            i1_=bj1-cj1;
            for(i_=cj1;i_<=cj2;i_++)
               c[k].Set(i_,c[k][i_]+v*b[r][i_+i1_]);
           }
        }
      //--- exit the function
      return;
     }
//--- A'*B'
   if(transa && transb)
     {
      //--- check
      if(arows*acols<brows*bcols)
        {
         for(r=bi1;r<=bi2;r++)
           {
            k=cj1+r-bi1;
            for(i=1;i<=crows;i++)
               work[i]=0.0;
            for(l=ai1;l<=ai2;l++)
              {
               //--- change values
               v=alpha*b[r][bj1+l-ai1];
               i1_=aj1-1;
               for(i_=1;i_<=crows;i_++)
                  work[i_]=work[i_]+v*a[l][i_+i1_];
              }
            i1_=1-ci1;
            for(i_=ci1;i_<=ci2;i_++)
               c[i_].Set(k,c[i_][k]+work[i_+i1_]);
           }
         //--- exit the function
         return;
        }
      else
        {
         for(l=aj1;l<=aj2;l++)
           {
            k=ai2-ai1+1;
            i1_=ai1-1;
            for(i_=1;i_<=k;i_++)
               work[i_]=a[i_+i1_][l];
            for(r=bi1;r<=bi2;r++)
              {
               //--- change values
               i1_=bj1-1;
               v=0.0;
               for(i_=1;i_<=k;i_++)
                  v+=work[i_]*b[r][i_+i1_];
               c[ci1+l-aj1].Set(cj1+r-bi1,c[ci1+l-aj1][cj1+r-bi1]+alpha*v);
              }
           }
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Work with the Hermitian matrix                                   |
//+------------------------------------------------------------------+
class CHblas
  {
public:
   //--- constructor, destructor
                     CHblas(void);
                    ~CHblas(void);
   //--- methods
   static void       HermitianMatrixVectorMultiply(CMatrixComplex &a,const bool isupper,const int i1,const int i2,al_complex &x[],al_complex &alpha,al_complex &y[]);
   static void       HermitianRank2Update(CMatrixComplex &a,const bool isupper,const int i1,const int i2,al_complex &x[],al_complex &y[],al_complex &t[],al_complex &alpha);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CHblas::CHblas(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CHblas::~CHblas(void)
  {

  }
//+------------------------------------------------------------------+
//| Multiply                                                         |
//+------------------------------------------------------------------+
static void CHblas::HermitianMatrixVectorMultiply(CMatrixComplex &a,const bool isupper,
                                                  const int i1,const int i2,al_complex &x[],
                                                  al_complex &alpha,al_complex &y[])
  {
//--- create variables
   int        i=0;
   int        ba1=0;
   int        ba2=0;
   int        by1=0;
   int        by2=0;
   int        bx1=0;
   int        bx2=0;
   int        n=i2-i1+1;
   int        i_=0;
   int        i1_=0;
   al_complex v=0;
//--- check
   if(n<=0)
      return;
//--- Let A = L + D + U, where
//---  L is strictly lower triangular (main diagonal is zero)
//---  D is diagonal
//---  U is strictly upper triangular (main diagonal is zero)
//--- A*x = L*x + D*x + U*x
//--- Calculate D*x first
   for(i=i1;i<=i2;i++)
      y[i-i1+1]=a[i][i]*x[i-i1+1];
//--- Add L*x + U*x
   if(isupper)
     {
      for(i=i1;i<=i2-1;i++)
        {
         //--- Add L*x to the result
         v=x[i-i1+1];
         by1=i-i1+2;
         by2=n;
         ba1=i+1;
         ba2=i2;
         i1_=ba1-by1;
         for(i_=by1;i_<=by2;i_++)
            y[i_]=y[i_]+v*CMath::Conj(a[i][i_+i1_]);
         //--- Add U*x to the result
         bx1=i-i1+2;
         bx2=n;
         ba1=i+1;
         ba2=i2;
         i1_=ba1-bx1;
         v=0.0;
         for(i_=bx1;i_<=bx2;i_++)
            v+=x[i_]*a[i][i_+i1_];
         y[i-i1+1]=y[i-i1+1]+v;
        }
     }
   else
     {
      for(i=i1+1;i<=i2;i++)
        {
         //--- Add L*x to the result
         bx1=1;
         bx2=i-i1;
         ba1=i1;
         ba2=i-1;
         i1_=ba1-bx1;
         v=0.0;
         for(i_=bx1;i_<=bx2;i_++)
            v+=x[i_]*a[i][i_+i1_];
         y[i-i1+1]=y[i-i1+1]+v;
         //--- change parameters
         v=x[i-i1+1];
         by1=1;
         by2=i-i1;
         ba1=i1;
         ba2=i-1;
         i1_=ba1-by1;
         //--- Add U*x to the result
         for(i_=by1;i_<=by2;i_++)
            y[i_]=y[i_]+v*CMath::Conj(a[i][i_+i1_]);
        }
     }
//--- get result
   for(i_=1;i_<=n;i_++)
      y[i_]=alpha*y[i_];
  }
//+------------------------------------------------------------------+
//| Update matrix                                                    |
//+------------------------------------------------------------------+
static void CHblas::HermitianRank2Update(CMatrixComplex &a,const bool isupper,
                                         const int i1,const int i2,al_complex &x[],
                                         al_complex &y[],al_complex &t[],al_complex &alpha)
  {
//--- create variables
   int        i=0;
   int        tp1=0;
   int        tp2=0;
   al_complex v=0;
   int        i_=0;
   int        i1_=0;
//--- check
   if(isupper)
     {
      for(i=i1;i<=i2;i++)
        {
         //--- change values
         tp1=i+1-i1;
         tp2=i2-i1+1;
         v=alpha*x[i+1-i1];
         //--- change t
         for(i_=tp1;i_<=tp2;i_++)
            t[i_]=v*CMath::Conj(y[i_]);
         v=CMath::Conj(alpha)*y[i+1-i1];
         //--- change t
         for(i_=tp1;i_<=tp2;i_++)
            t[i_]=t[i_]+v*CMath::Conj(x[i_]);
         i1_=tp1-i;
         //--- change a
         for(i_=i;i_<=i2;i_++)
            a[i].Set(i_,a[i][i_]+t[i_+i1_]);
        }
     }
   else
     {
      for(i=i1;i<=i2;i++)
        {
         //--- change values
         tp1=1;
         tp2=i+1-i1;
         v=alpha*x[i+1-i1];
         //--- change t
         for(i_=tp1;i_<=tp2;i_++)
            t[i_]=v*CMath::Conj(y[i_]);
         v=CMath::Conj(alpha)*y[i+1-i1];
         //--- change t
         for(i_=tp1;i_<=tp2;i_++)
            t[i_]=t[i_]+v*CMath::Conj(x[i_]);
         i1_=tp1-i1;
         //--- change a
         for(i_=i1;i_<=i;i_++)
            a[i].Set(i_,a[i][i_]+t[i_+i1_]);
        }
     }
  }
//+------------------------------------------------------------------+
//| Reflections                                                      |
//+------------------------------------------------------------------+
class CReflections
  {
public:
   //--- constructor, destructor
                     CReflections(void);
                    ~CReflections(void);
   //--- methods
   static void       GenerateReflection(double &x[],const int n,double &tau);
   static void       ApplyReflectionFromTheLeft(CMatrixDouble &c,const double tau,const double &v[],const int m1,const int m2,const int n1,const int n2,double &work[]);
   static void       ApplyReflectionFromTheRight(CMatrixDouble &c,const double tau,const double &v[],const int m1,const int m2,const int n1,const int n2,double &work[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CReflections::CReflections(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CReflections::~CReflections(void)
  {

  }
//+------------------------------------------------------------------+
//| Generation of an elementary reflection transformation            |
//| The subroutine generates elementary reflection H of order N, so  |
//| that, for a given X, the following equality holds true:          |
//|     ( X(1) )   ( Beta )                                          |
//| H * (  ..  ) = (  0   )                                          |
//|     ( X(n) )   (  0   )                                          |
//| where                                                            |
//|               ( V(1) )                                           |
//| H = 1 - Tau * (  ..  ) * ( V(1), ..., V(n) )                     |
//|               ( V(n) )                                           |
//| where the first component of vector V equals 1.                  |
//| Input parameters:                                                |
//|     X   -   vector. Array whose index ranges within [1..N].      |
//|     N   -   reflection order.                                    |
//| Output parameters:                                               |
//|     X   -   components from 2 to N are replaced with vector V.   |
//|             The first component is replaced with parameter Beta. |
//|     Tau -   scalar value Tau. If X is a null vector, Tau equals  |
//|             0, otherwise 1 <= Tau <= 2.                          |
//| This subroutine is the modification of the DLARFG subroutines    |
//| from the LAPACK library.                                         |
//| MODIFICATIONS:                                                   |
//|     24.12.2005 sign(Alpha) was replaced with an analogous to the |
//|     Fortran SIGN code.                                           |
//|   -- LAPACK auxiliary routine (version 3.0) --                   |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      September 30, 1994                                          |
//+------------------------------------------------------------------+
static void CReflections::GenerateReflection(double &x[],const int n,
                                             double &tau)
  {
//--- check
   if(n<=1)
     {
      tau=0;
      //--- exit the function
      return;
     }
//--- create variables
   int    j=0;
   double alpha=0;
   double xnorm=0;
   double v=0;
   double beta=0;
   double mx=0;
   double s=0;
   int    i_=0;
//--- Scale if needed (to avoid overflow/underflow during intermediate
//--- calculations).
   for(j=1;j<=n;j++)
      mx=MathMax(MathAbs(x[j]),mx);
   s=1;
//--- check
   if(mx!=0.0)
     {
      //--- check
      if(mx<=CMath::m_minrealnumber/CMath::m_machineepsilon)
        {
         //--- change parameters
         s=CMath::m_minrealnumber/CMath::m_machineepsilon;
         v=1/s;
         //--- change x
         for(i_=1;i_<=n;i_++)
            x[i_]=v*x[i_];
         mx=mx*v;
        }
      else
        {
         //--- check
         if(mx>=CMath::m_maxrealnumber*CMath::m_machineepsilon)
           {
            //--- change parameters
            s=CMath::m_maxrealnumber*CMath::m_machineepsilon;
            v=1/s;
            //--- change x
            for(i_=1;i_<=n;i_++)
               x[i_]=v*x[i_];
            mx=mx*v;
           }
        }
     }
//--- XNORM = DNRM2( N-1, X, INCX )
   alpha=x[1];
   xnorm=0;
//--- check
   if(mx!=0.0)
     {
      for(j=2;j<=n;j++)
         xnorm=xnorm+CMath::Sqr(x[j]/mx);
      xnorm=MathSqrt(xnorm)*mx;
     }
//--- check
   if(xnorm==0.0)
     {
      //--- H = I
      tau=0;
      x[1]=x[1]*s;
      //--- exit the function
      return;
     }
//--- general case
   mx=MathMax(MathAbs(alpha),MathAbs(xnorm));
   beta=-(mx*MathSqrt(CMath::Sqr(alpha/mx)+CMath::Sqr(xnorm/mx)));
//--- check
   if(alpha<0.0)
      beta=-beta;
//--- change parameters
   tau=(beta-alpha)/beta;
   v=1/(alpha-beta);
//--- change x
   for(i_=2;i_<=n;i_++)
      x[i_]=v*x[i_];
   x[1]=beta;
//--- Scale back outputs
   x[1]=x[1]*s;
  }
//+------------------------------------------------------------------+
//| Application of an elementary reflection to a rectangular matrix  |
//| of size MxN                                                      |
//| The algorithm pre-multiplies the matrix by an elementary         |
//| reflection transformation which is given by column V and scalar  |
//| Tau (see the description of the GenerateReflection procedure).   |
//| Not the whole matrix but only a part of it is transformed (rows  |
//| from M1 to M2, columns from N1 to N2). Only the elements of this |
//| submatrix are changed.                                           |
//| Input parameters:                                                |
//|     C       -   matrix to be transformed.                        |
//|     Tau     -   scalar defining the transformation.              |
//|     V       -   column defining the transformation.              |
//|                 Array whose index ranges within [1..M2-M1+1].    |
//|     M1, M2  -   range of rows to be transformed.                 |
//|     N1, N2  -   range of columns to be transformed.              |
//|     WORK    -   working array whose indexes goes from N1 to N2.  |
//| Output parameters:                                               |
//|     C       -   the result of multiplying the input matrix C by  |
//|                 the transformation matrix which is given by Tau  | 
//|                 and V. If N1>N2 or M1>M2, C is not modified.     |
//|   -- LAPACK auxiliary routine (version 3.0) --                   |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      September 30, 1994                                          |
//+------------------------------------------------------------------+
static void CReflections::ApplyReflectionFromTheLeft(CMatrixDouble &c,const double tau,
                                                     const double &v[],const int m1,
                                                     const int m2,const int n1,
                                                     const int n2,double &work[])
  {
//--- check
   if(tau==0.0 || n1>n2 || m1>m2)
      return;
//--- create variables
   double t=0;
   int    i=0;
   int    vm=0;
   int    i_=0;
//--- w := C' * v
   vm=m2-m1+1;
   for(i=n1;i<=n2;i++)
      work[i]=0;
   for(i=m1;i<=m2;i++)
     {
      t=v[i+1-m1];
      //--- change array
      for(i_=n1;i_<=n2;i_++)
         work[i_]+=t*c[i][i_];
     }
//--- C := C - tau * v * w'
   for(i=m1;i<=m2;i++)
     {
      t=v[i-m1+1]*tau;
      for(i_=n1;i_<=n2;i_++)
         c[i].Set(i_,c[i][i_]-t*work[i_]);
     }
  }
//+------------------------------------------------------------------+
//| Application of an elementary reflection to a rectangular matrix  |
//| of size MxN                                                      |
//| The algorithm post-multiplies the matrix by an elementary        |
//| reflection transformation which is given by column V and scalar  |
//| Tau (see the description of the GenerateReflection procedure).   |
//| Not the whole matrix but only a part of it is transformed (rows  |
//| from M1 to M2, columns from N1 to N2). Only the elements of this |
//| submatrix are changed.                                           |
//| Input parameters:                                                |
//|     C       -   matrix to be transformed.                        |
//|     Tau     -   scalar defining the transformation.              |
//|     V       -   column defining the transformation.              |
//|                 Array whose index ranges within [1..N2-N1+1].    |
//|     M1, M2  -   range of rows to be transformed.                 |
//|     N1, N2  -   range of columns to be transformed.              |
//|     WORK    -   working array whose indexes goes from M1 to M2.  |
//| Output parameters:                                               |
//|     C       -   the result of multiplying the input matrix C by  |
//|                 the transformation matrix which is given by Tau  |
//|                 and V. If N1>N2 or M1>M2, C is not modified.     |
//|   -- LAPACK auxiliary routine (version 3.0) --                   |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      September 30, 1994                                          |
//+------------------------------------------------------------------+
static void CReflections::ApplyReflectionFromTheRight(CMatrixDouble &c,const double tau,
                                                      const double &v[],const int m1,
                                                      const int m2,const int n1,
                                                      const int n2,double &work[])
  {
//--- check
   if(tau==0.0 || n1>n2 || m1>m2)
      return;
//--- create variables
   double t=0;
   int    i=0;
   int    vm=n2-n1+1;
   int    i_=0;
   int    i1_=0;
//--- change matrix
   for(i=m1;i<=m2;i++)
     {
      i1_=1-n1;
      t=0.0;
      //--- calculation parameters
      for(i_=n1;i_<=n2;i_++)
         t+=c[i][i_]*v[i_+i1_];
      t=t*tau;
      i1_=1-n1;
      for(i_=n1;i_<=n2;i_++)
         c[i].Set(i_,c[i][i_]-t*v[i_+i1_]);
     }
  }
//+------------------------------------------------------------------+
//| Complex reflections                                              |
//+------------------------------------------------------------------+
class CComplexReflections
  {
public:
   //--- constructor, destructor
   void              CComplexReflections(void);
   void             ~CComplexReflections(void);
   //--- methods
   static void       ComplexGenerateReflection(al_complex &x[],const int n,al_complex &tau);
   static void       ComplexApplyReflectionFromTheLeft(CMatrixComplex &c,al_complex &tau,al_complex &v[],const int m1,const int m2,const int n1,const int n2,al_complex &work[]);
   static void       ComplexApplyReflectionFromTheRight(CMatrixComplex &c,al_complex &tau,al_complex &v[],const int m1,const int m2,const int n1,const int n2,al_complex &work[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CComplexReflections::CComplexReflections(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CComplexReflections::~CComplexReflections(void)
  {

  }
//+------------------------------------------------------------------+
//| Generation of an elementary complex reflection transformation    |
//| The subroutine generates elementary complex reflection H of      |
//| order N, so that, for a given X, the following equality holds    |
//| true:                                                            |
//|      ( X(1) )   ( Beta )                                         |
//| H' * (  ..  ) = (  0   ),   H'*H = I,   Beta is a real number    |
//|      ( X(n) )   (  0   )                                         |
//| where                                                            |
//|               ( V(1) )                                           |
//| H = 1 - Tau * (  ..  ) * ( conj(V(1)), ..., conj(V(n)) )         |
//|               ( V(n) )                                           |
//| where the first component of vector V equals 1.                  |
//| Input parameters:                                                |
//|     X   -   vector. Array with elements [1..N].                  |
//|     N   -   reflection order.                                    |
//| Output parameters:                                               |
//|     X   -   components from 2 to N are replaced by vector V.     |
//|             The first component is replaced with parameter Beta. |
//|     Tau -   scalar value Tau.                                    |
//| This subroutine is the modification of CLARFG subroutines  from  |
//| the LAPACK library. It has similar functionality except for the  |
//| fact that it  doesn?t handle errors when intermediate results    |
//| cause an overflow.                                               |
//|   -- LAPACK auxiliary routine (version 3.0) --                   |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      September 30, 1994                                          |
//+------------------------------------------------------------------+
static void CComplexReflections::ComplexGenerateReflection(al_complex &x[],
                                                           const int n,
                                                           al_complex &tau)
  {
//--- check
   if(n<=0)
     {
      tau=0;
      //--- exit the function
      return;
     }
//--- create variables
   int        j=0;
   al_complex alpha=0;
   double     alphi=0;
   double     alphr=0;
   double     beta=0;
   double     xnorm=0;
   double     mx=0;
   al_complex t=0;
   double     s=1;
   al_complex v=0;
   int        i_=0;
   al_complex One(1,0);
//--- Scale if needed (to avoid overflow/underflow during intermediate
//--- calculations).
   for(j=1;j<=n;j++)
      mx=MathMax(CMath::AbsComplex(x[j]),mx);
//--- check
   if(mx!=0)
     {
      //--- check
      if(mx<1)
        {
         s=MathSqrt(CMath::m_minrealnumber);
         v=1/s;
         //--- change x
         for(i_=1;i_<=n;i_++)
            x[i_]=v*x[i_];
        }
      else
        {
         s=MathSqrt(CMath::m_maxrealnumber);
         v=1/s;
         //--- change x
         for(i_=1;i_<=n;i_++)
            x[i_]=v*x[i_];
        }
     }
//--- calculate
   alpha=x[1];
   mx=0;
   for(j=2;j<=n;j++)
      mx=MathMax(CMath::AbsComplex(x[j]),mx);
   xnorm=0;
//--- check
   if(mx!=0)
     {
      for(j=2;j<=n;j++)
        {
         t=x[j]/mx;
         xnorm=xnorm+(t*CMath::Conj(t)).re;
        }
      xnorm=MathSqrt(xnorm)*mx;
     }
//--- change parameters
   alphr=alpha.re;
   alphi=alpha.im;
//--- check
   if((xnorm==0) && (alphi==0))
     {
      tau=0;
      x[1]=x[1]*s;
      //--- exit the function
      return;
     }
//--- change parameters
   mx=MathMax(MathAbs(alphr),MathAbs(alphi));
   mx=MathMax(mx,MathAbs(xnorm));
   beta=-(mx*MathSqrt(CMath::Sqr(alphr/mx)+CMath::Sqr(alphi/mx)+CMath::Sqr(xnorm/mx)));
//--- check
   if(alphr<0)
      beta=-beta;
//--- change parameters
   tau.re=(beta-alphr)/beta;
   tau.im=-(alphi/beta);
   alpha=One/(alpha-beta);
//--- check
   if(n>1)
     {
      //--- change x
      for(i_=2;i_<=n;i_++)
         x[i_]=alpha*x[i_];
     }
   alpha=beta;
   x[1]=alpha;
//--- Scale back
   x[1]=x[1]*s;
  }
//+------------------------------------------------------------------+
//| Application of an elementary reflection to a rectangular matrix  |
//| of size MxN                                                      |
//| The  algorithm  pre-multiplies  the  matrix  by  an  elementary  |
//| reflection transformation  which  is  given  by  column  V  and  | 
//| scalar  Tau (see the description of the GenerateReflection). Not | 
//| the whole matrix  but  only  a part of it is transformed (rows   | 
//| from M1 to M2, columns from N1 to N2). Only the elements of this | 
//| submatrix are changed.                                           |
//| Note: the matrix is multiplied by H, not by H'.   If  it  is     | 
//| required  to multiply the matrix by H', it is necessary to pass  |
//| Conj(Tau) instead of Tau.                                        |
//| Input parameters:                                                |
//|     C       -   matrix to be transformed.                        |
//|     Tau     -   scalar defining transformation.                  |
//|     V       -   column defining transformation.                  |
//|                 Array whose index ranges within [1..M2-M1+1]     |
//|     M1, M2  -   range of rows to be transformed.                 |
//|     N1, N2  -   range of columns to be transformed.              |
//|     WORK    -   working array whose index goes from N1 to N2.    |
//| Output parameters:                                               |
//|     C       -   the result of multiplying the input matrix C by  |
//|                 the transformation matrix which is given by Tau  |
//|                 and V. If N1>N2 or M1>M2, C is not modified.     |
//|   -- LAPACK auxiliary routine (version 3.0) --                   |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      September 30, 1994                                          |
//+------------------------------------------------------------------+
static void CComplexReflections::ComplexApplyReflectionFromTheLeft(CMatrixComplex &c,
                                                                   al_complex &tau,
                                                                   al_complex &v[],
                                                                   const int m1,
                                                                   const int m2,
                                                                   const int n1,
                                                                   const int n2,
                                                                   al_complex &work[])
  {
   al_complex Zero(0,0);
//--- check
   if(tau==Zero || n1>n2 || m1>m2)
      return;
//--- create variables
   al_complex t=0;
   int        i=0;
   int        vm=0;
   int        i_=0;
//--- w := C^T * conj(v)
   vm=m2-m1+1;
   for(i=n1;i<=n2;i++)
      work[i]=0;
   for(i=m1;i<=m2;i++)
     {
      t=CMath::Conj(v[i+1-m1]);
      for(i_=n1;i_<=n2;i_++)
         work[i_]=work[i_]+t*c[i][i_];
     }
//--- C := C - tau * v * w^T
   for(i=m1;i<=m2;i++)
     {
      t=v[i-m1+1]*tau;
      for(i_=n1;i_<=n2;i_++)
         c[i].Set(i_,c[i][i_]-t*work[i_]);
     }
  }
//+------------------------------------------------------------------+
//| Application of an elementary reflection to a rectangular matrix  |
//| of size MxN                                                      |
//| The  algorithm  post-multiplies  the  matrix  by  an elementary  |
//| reflection transformation  which  is  given  by  column  V  and  |
//| scalar  Tau (see the description  of  the  GenerateReflection).  | 
//| Not the whole matrix but only a part  of  it  is  transformed    |
//| (rows from M1 to M2, columns from N1 to N2). Only the elements   |
//| of this submatrix are changed.                                   |
//| Input parameters:                                                |
//|     C       -   matrix to be transformed.                        |
//|     Tau     -   scalar defining transformation.                  |
//|     V       -   column defining transformation.                  |
//|                 Array whose index ranges within [1..N2-N1+1]     |
//|     M1, M2  -   range of rows to be transformed.                 |
//|     N1, N2  -   range of columns to be transformed.              |
//|     WORK    -   working array whose index goes from M1 to M2.    |
//| Output parameters:                                               |
//|     C       -   the result of multiplying the input matrix C by  |
//|                 the transformation matrix which is given by Tau  |
//|                 and V. If N1>N2 or M1>M2, C is not modified.     |
//|   -- LAPACK auxiliary routine (version 3.0) --                   |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      September 30, 1994                                          |
//+------------------------------------------------------------------+
static void CComplexReflections::ComplexApplyReflectionFromTheRight(CMatrixComplex &c,
                                                                    al_complex &tau,
                                                                    al_complex &v[],
                                                                    const int m1,
                                                                    const int m2,
                                                                    const int n1,
                                                                    const int n2,
                                                                    al_complex &work[])
  {
   al_complex Zero(0,0);
//--- check
   if(tau==Zero || n1>n2 || m1>m2)
      return;
//--- create variables
   al_complex t=0;
   int        i=0;
   int        vm=0;
   int        i_=0;
   int        i1_=0;
//--- w := C * v
   vm=n2-n1+1;
   for(i=m1;i<=m2;i++)
     {
      i1_=1-n1;
      t=0.0;
      //--- change values
      for(i_=n1;i_<=n2;i_++)
         t+=c[i][i_]*v[i_+i1_];
      work[i]=t;
     }
//--- C := C - w * conj(v^T)
   for(i_=1;i_<=vm;i_++)
      v[i_]=CMath::Conj(v[i_]);
//--- get result
   for(i=m1;i<=m2;i++)
     {
      t=work[i]*tau;
      i1_=1-n1;
      for(i_=n1;i_<=n2;i_++)
         c[i].Set(i_,c[i][i_]-t*v[i_+i1_]);
     }
   for(i_=1;i_<=vm;i_++)
      v[i_]=CMath::Conj(v[i_]);
  }
//+------------------------------------------------------------------+
//| Work with the symmetric matrix                                   |
//+------------------------------------------------------------------+
class CSblas
  {
public:
   //--- constructor, destructor
                     CSblas(void);
                    ~CSblas(void);
   //--- methods
   static void       SymmetricMatrixVectorMultiply(const CMatrixDouble &a,const bool isupper,const int i1,const int i2,const double &x[],const double alpha,double &y[]);
   static void       SymmetricRank2Update(CMatrixDouble &a,const bool isupper,const int i1,const int i2,const double &x[],const double &y[],double &t[],const double alpha);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSblas::CSblas(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSblas::~CSblas(void)
  {

  }
//+------------------------------------------------------------------+
//| Multiply                                                         |
//+------------------------------------------------------------------+
static void CSblas::SymmetricMatrixVectorMultiply(const CMatrixDouble &a,
                                                  const bool isupper,
                                                  const int i1,const int i2,
                                                  const double &x[],
                                                  const double alpha,
                                                  double &y[])
  {
//--- create variables
   int    i=0;
   int    ba1=0;
   int    ba2=0;
   int    by1=0;
   int    by2=0;
   int    bx1=0;
   int    bx2=0;
   int    n=i2-i1+1;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(n<=0)
      return;
//--- Let A = L + D + U, where
//---  L is strictly lower triangular (main diagonal is zero)
//---  D is diagonal
//---  U is strictly upper triangular (main diagonal is zero)
//--- A*x = L*x + D*x + U*x
//--- Calculate D*x first
   for(i=i1;i<=i2;i++)
      y[i-i1+1]=a[i][i]*x[i-i1+1];
//--- Add L*x + U*x
   if(isupper)
     {
      for(i=i1;i<i2;i++)
        {
         //--- Add L*x to the result
         v=x[i-i1+1];
         by1=i-i1+2;
         by2=n;
         ba1=i+1;
         ba2=i2;
         i1_=ba1-by1;
         for(i_=by1;i_<=by2;i_++)
            y[i_]=y[i_]+v*a[i][i_+i1_];
         //--- Add U*x to the result
         bx1=i-i1+2;
         bx2=n;
         ba1=i+1;
         ba2=i2;
         i1_=ba1-bx1;
         v=0.0;
         for(i_=bx1;i_<=bx2;i_++)
            v+=x[i_]*a[i][i_+i1_];
         y[i-i1+1]=y[i-i1+1]+v;
        }
     }
   else
     {
      for(i=i1+1;i<=i2;i++)
        {
         //--- Add L*x to the result
         bx1=1;
         bx2=i-i1;
         ba1=i1;
         ba2=i-1;
         i1_=(ba1)-(bx1);
         v=0.0;
         for(i_=bx1;i_<=bx2;i_++)
            v+=x[i_]*a[i][i_+i1_];
         y[i-i1+1]=y[i-i1+1]+v;
         //--- Add U*x to the result
         v=x[i-i1+1];
         by1=1;
         by2=i-i1;
         ba1=i1;
         ba2=i-1;
         i1_=ba1-by1;
         for(i_=by1;i_<=by2;i_++)
            y[i_]=y[i_]+v*a[i][i_+i1_];
        }
     }
//--- get result
   for(i_=1;i_<=n;i_++)
      y[i_]=alpha*y[i_];
  }
//+------------------------------------------------------------------+
//| Update matrix                                                    |
//+------------------------------------------------------------------+
static void CSblas::SymmetricRank2Update(CMatrixDouble &a,const bool isupper,
                                         const int i1,const int i2,
                                         const double &x[],const double &y[],
                                         double &t[],const double alpha)
  {
//--- create variables
   int    i=0;
   int    tp1=0;
   int    tp2=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(isupper)
     {
      for(i=i1;i<=i2;i++)
        {
         //--- change values
         tp1=i+1-i1;
         tp2=i2-i1+1;
         v=x[i+1-i1];
         //--- change t
         for(i_=tp1;i_<=tp2;i_++)
            t[i_]=v*y[i_];
         v=y[i+1-i1];
         //--- change t
         for(i_=tp1;i_<=tp2;i_++)
            t[i_]=t[i_]+v*x[i_];
         //--- change t
         for(i_=tp1;i_<=tp2;i_++)
            t[i_]=alpha*t[i_];
         i1_=tp1-i;
         //--- change a
         for(i_=i;i_<=i2;i_++)
            a[i].Set(i_,a[i][i_]+t[i_+i1_]);
        }
     }
   else
     {
      for(i=i1;i<=i2;i++)
        {
         //--- change values
         tp1=1;
         tp2=i+1-i1;
         v=x[i+1-i1];
         //--- change t
         for(i_=tp1;i_<=tp2;i_++)
            t[i_]=v*y[i_];
         v=y[i+1-i1];
         //--- change t
         for(i_=tp1;i_<=tp2;i_++)
            t[i_]=t[i_]+v*x[i_];
         //--- change t
         for(i_=tp1;i_<=tp2;i_++)
            t[i_]=alpha*t[i_];
         i1_=tp1-i1;
         //--- change a
         for(i_=i1;i_<=i;i_++)
            a[i].Set(i_,a[i][i_]+t[i_+i1_]);
        }
     }
  }
//+------------------------------------------------------------------+
//| Rotations                                                        |
//+------------------------------------------------------------------+
class CRotations
  {
public:
   //--- constructor, destructor
                     CRotations(void);
                    ~CRotations(void);
   //--- methods
   static void       ApplyRotationsFromTheLeft(const bool isforward,const int m1,const int m2,const int n1,const int n2,double &c[],double &s[],CMatrixDouble &a,double &work[]);
   static void       ApplyRotationsFromTheRight(const bool isforward,const int m1,const int m2,const int n1,const int n2,double &c[],double &s[],CMatrixDouble &a,double &work[]);
   static void       GenerateRotation(const double f,const double g,double &cs,double &sn,double &r);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CRotations::CRotations(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CRotations::~CRotations(void)
  {

  }
//+------------------------------------------------------------------+
//| Application of a sequence of  elementary rotations to a matrix   |
//| The algorithm pre-multiplies the matrix by a sequence of rotation|
//| transformations which is given by arrays C and S. Depending on   |
//| the value of the IsForward parameter either 1 and 2, 3 and 4 and |
//| so on (if IsForward=true) rows are rotated, or the rows N and    |
//| N-1, N-2 and N-3 and so on, are rotated.                         |
//| Not the whole matrix but only a part of it is transformed (rows  |
//| from M1 to M2, columns from N1 to N2). Only the elements of this |
//| submatrix are changed.                                           |
//| Input parameters:                                                |
//|     IsForward   -   the sequence of the rotation application.    |
//|     M1,M2       -   the range of rows to be transformed.         |
//|     N1, N2      -   the range of columns to be transformed.      |
//|     C,S         -   transformation coefficients.                 |
//|                     Array whose index ranges within [1..M2-M1].  |
//|     A           -   processed matrix.                            |
//|     WORK        -   working array whose index ranges within      |
//|                     [N1..N2].                                    |
//| Output parameters:                                               |
//|     A           -   transformed matrix.                          |
//| Utility subroutine.                                              |
//+------------------------------------------------------------------+
static void CRotations::ApplyRotationsFromTheLeft(const bool isforward,
                                                  const int m1,const int m2,
                                                  const int n1,const int n2,
                                                  double &c[],double &s[],
                                                  CMatrixDouble &a,double &work[])
  {
//--- create variables
   int    j=0;
   int    jp1=0;
   double ctemp=0;
   double stemp=0;
   double temp=0;
   int    i_=0;
//--- check
   if(m1>m2 || n1>n2)
      return;
//--- Form  P * A
   if(isforward)
     {
      //--- check
      if(n1!=n2)
        {
         //--- Common case: N1<>N2
         for(j=m1;j<m2;j++)
           {
            ctemp=c[j-m1+1];
            stemp=s[j-m1+1];
            //--- check
            if(ctemp!=1.0 || stemp!=0.0)
              {
               jp1=j+1;
               //--- prepare array
               for(i_=n1;i_<=n2;i_++)
                  work[i_]=ctemp*a[jp1][i_];
               for(i_=n1;i_<=n2;i_++)
                  work[i_]=work[i_]-stemp*a[j][i_];
               //--- get result
               for(i_=n1;i_<=n2;i_++)
                  a[j].Set(i_,ctemp*a[j][i_]);
               for(i_=n1;i_<=n2;i_++)
                  a[j].Set(i_,a[j][i_]+stemp*a[jp1][i_]);
               for(i_=n1;i_<=n2;i_++)
                  a[jp1].Set(i_,work[i_]);
              }
           }
        }
      else
        {
         //--- Special case: N1=N2
         for(j=m1;j<m2;j++)
           {
            ctemp=c[j-m1+1];
            stemp=s[j-m1+1];
            //--- check
            if(ctemp!=1.0 || stemp!=0.0)
              {
               temp=a[j+1][n1];
               //--- get result
               a[j+1].Set(n1,ctemp*temp-stemp*a[j][n1]);
               a[j].Set(n1,stemp*temp+ctemp*a[j][n1]);
              }
           }
        }
     }
   else
     {
      if(n1!=n2)
        {
         //--- Common case: N1<>N2
         for(j=m2-1;j>=m1;j--)
           {
            ctemp=c[j-m1+1];
            stemp=s[j-m1+1];
            //--- check
            if(ctemp!=1.0 || stemp!=0.0)
              {
               jp1=j+1;
               //--- prepare array
               for(i_=n1;i_<=n2;i_++)
                  work[i_]=ctemp*a[jp1][i_];
               for(i_=n1;i_<=n2;i_++)
                  work[i_]=work[i_]-stemp*a[j][i_];
               for(i_=n1;i_<=n2;i_++)
                  a[j].Set(i_,ctemp*a[j][i_]);
               //--- get result
               for(i_=n1;i_<=n2;i_++)
                  a[j].Set(i_,a[j][i_]+stemp*a[jp1][i_]);
               for(i_=n1;i_<=n2;i_++)
                  a[jp1].Set(i_,work[i_]);
              }
           }
        }
      else
        {
         //--- Special case: N1=N2
         for(j=m2-1;j>=m1;j--)
           {
            ctemp=c[j-m1+1];
            stemp=s[j-m1+1];
            //--- check
            if(ctemp!=1.0 || stemp!=0.0)
              {
               temp=a[j+1][n1];
               //--- get result
               a[j+1].Set(n1,ctemp*temp-stemp*a[j][n1]);
               a[j].Set(n1,stemp*temp+ctemp*a[j][n1]);
              }
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Application of a sequence of  elementary rotations to a matrix   |
//| The algorithm post-multiplies the matrix by a sequence of        |
//| rotation transformations which is given by arrays C and S.       |
//| Depending on the value of the IsForward parameter either 1 and 2,|
//| 3 and 4 and so on (if IsForward=true) rows are rotated, or the   |
//| rows N and N-1, N-2 and N-3 and so on are rotated.               |
//| Not the whole matrix but only a part of it is transformed (rows  | 
//| from M1 to M2, columns from N1 to N2). Only the elements of this |
//| submatrix are changed.                                           |
//| Input parameters:                                                |
//|     IsForward   -   the sequence of the rotation application.    |
//|     M1,M2       -   the range of rows to be transformed.         |
//|     N1, N2      -   the range of columns to be transformed.      |
//|     C,S         -   transformation coefficients.                 |
//|                     Array whose index ranges within [1..N2-N1].  |
//|     A           -   processed matrix.                            |
//|     WORK        -   working array whose index ranges within      |
//|                     [M1..M2].                                    |
//| Output parameters:                                               |
//|     A           -   transformed matrix.                          |
//| Utility subroutine.                                              |
//+------------------------------------------------------------------+
static void CRotations::ApplyRotationsFromTheRight(const bool isforward,
                                                   const int m1,const int m2,
                                                   const int n1,const int n2,
                                                   double &c[],double &s[],
                                                   CMatrixDouble &a,
                                                   double &work[])
  {
//--- create variables
   int    j=0;
   int    jp1=0;
   double ctemp=0;
   double stemp=0;
   double temp=0;
   int    i_=0;
//--- Form A * P'
   if(isforward)
     {
      //--- check
      if(m1!=m2)
        {
         //--- Common case: M1<>M2
         for(j=n1;j<=n2-1;j++)
           {
            ctemp=c[j-n1+1];
            stemp=s[j-n1+1];
            //--- check
            if(ctemp!=1.0 || stemp!=0.0)
              {
               jp1=j+1;
               //--- prepare array
               for(i_=m1;i_<=m2;i_++)
                  work[i_]=ctemp*a[i_][jp1];
               for(i_=m1;i_<=m2;i_++)
                  work[i_]=work[i_]-stemp*a[i_][j];
               //--- get result
               for(i_=m1;i_<=m2;i_++)
                  a[i_].Set(j,ctemp*a[i_][j]);
               for(i_=m1;i_<=m2;i_++)
                  a[i_].Set(j,a[i_][j]+stemp*a[i_][jp1]);
               for(i_=m1;i_<=m2;i_++)
                  a[i_].Set(jp1,work[i_]);
              }
           }
        }
      else
        {
         //--- Special case: M1=M2
         for(j=n1;j<=n2-1;j++)
           {
            ctemp=c[j-n1+1];
            stemp=s[j-n1+1];
            //--- check
            if(ctemp!=1.0 || stemp!=0.0)
              {
               temp=a[m1][j+1];
               //--- get result
               a[m1].Set(j+1,ctemp*temp-stemp*a[m1][j]);
               a[m1].Set(j,stemp*temp+ctemp*a[m1][j]);
              }
           }
        }
     }
   else
     {
      //--- check
      if(m1!=m2)
        {
         //--- Common case: M1<>M2
         for(j=n2-1;j>=n1;j--)
           {
            ctemp=c[j-n1+1];
            stemp=s[j-n1+1];
            //--- check
            if(ctemp!=1.0 || stemp!=0.0)
              {
               jp1=j+1;
               //--- prepare array
               for(i_=m1;i_<=m2;i_++)
                  work[i_]=ctemp*a[i_][jp1];
               for(i_=m1;i_<=m2;i_++)
                  work[i_]=work[i_]-stemp*a[i_][j];
               //--- get result
               for(i_=m1;i_<=m2;i_++)
                  a[i_].Set(j,ctemp*a[i_][j]);
               for(i_=m1;i_<=m2;i_++)
                  a[i_].Set(j,a[i_][j]+stemp*a[i_][jp1]);
               for(i_=m1;i_<=m2;i_++)
                  a[i_].Set(jp1,work[i_]);
              }
           }
        }
      else
        {
         //--- Special case: M1=M2
         for(j=n2-1;j>=n1;j--)
           {
            ctemp=c[j-n1+1];
            stemp=s[j-n1+1];
            //--- check
            if(ctemp!=1.0 || stemp!=0.0)
              {
               temp=a[m1][j+1];
               //--- get result
               a[m1].Set(j+1,ctemp*temp-stemp*a[m1][j]);
               a[m1].Set(j,stemp*temp+ctemp*a[m1][j]);
              }
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| The subroutine generates the elementary rotation, so that:       |
//| [  CS  SN  ]  .  [ F ]  =  [ R ]                                 |
//| [ -SN  CS  ]     [ G ]     [ 0 ]                                 |
//| CS**2 + SN**2 = 1                                                |
//+------------------------------------------------------------------+
static void CRotations::GenerateRotation(const double f,const double g,
                                         double &cs,double &sn,double &r)
  {
//--- create variables
   double f1=0;
   double g1=0;
//--- check
   if(g==0.0)
     {
      //--- get result
      cs=1;
      sn=0;
      r=f;
     }
   else
     {
      //--- check
      if(f==0.0)
        {
         //--- get result
         cs=0;
         sn=1;
         r=g;
        }
      else
        {
         f1=f;
         g1=g;
         //--- check
         if(MathAbs(f1)>MathAbs(g1))
            r=MathAbs(f1)*MathSqrt(1+CMath::Sqr(g1/f1));
         else
            r=MathAbs(g1)*MathSqrt(1+CMath::Sqr(f1/g1));
         cs=f1/r;
         sn=g1/r;
         //--- check
         if(MathAbs(f)>MathAbs(g) && cs<0.0)
           {
            //--- get result
            cs=-cs;
            sn=-sn;
            r=-r;
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Schur decomposition of Hessenberg matrix                         |
//+------------------------------------------------------------------+
class CHsSchur
  {
private:
   //--- private methods
   static void       InternalAuxSchur(const bool wantt,const bool wantz,const int n,const int ilo,const int ihi,CMatrixDouble &h,double &wr[],double &wi[],const int iloz,const int ihiz,CMatrixDouble &z,double &work[],double &workv3[],double &workc1[],double &works1[],int &info);
   static void       Aux2x2Schur(double &a,double &b,double &c,double &d,double &rt1r,double &rt1i,double &rt2r,double &rt2i,double &cs,double &sn);
   static double     ExtSchurSign(const double a,const double b);
   static int        ExtSchurSignToone(const double b);
public:
   //--- constructor, destructor
                     CHsSchur(void);
                    ~CHsSchur(void);
   //--- public methods
   static bool       UpperHessenbergSchurDecomposition(CMatrixDouble &h,const int n,CMatrixDouble &s);
   static void       InternalSchurDecomposition(CMatrixDouble &h,const int n,const int tneeded,const int zneeded,double &wr[],double &wi[],CMatrixDouble &z,int &info);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CHsSchur::CHsSchur(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CHsSchur::~CHsSchur(void)
  {

  }
//+------------------------------------------------------------------+
//| Subroutine performing the Schur decomposition of a matrix in     |
//| upper Hessenberg form using the QR algorithm with multiple       |
//| shifts.                                                          |
//| The  source matrix  H is  represented as S'*H*S = T, where H -   |
//| matrix in upper Hessenberg form, S - orthogonal matrix (Schur    |
//| vectors), T - upper quasi-triangular matrix (with blocks of sizes| 
//|  1x1  and  2x2  on  the main diagonal).                          |
//| Input parameters:                                                |
//|     H   -   matrix to be decomposed.                             |
//|             Array whose indexes range within [1..N, 1..N].       |
//|     N   -   size of H, N>=0.                                     |
//| Output parameters:                                               |
//|     H   ?   contains the matrix T.                               |
//|             Array whose indexes range within [1..N, 1..N].       |
//|             All elements below the blocks on the main diagonal   |
//|             are equal to 0.                                      |
//|     S   -   contains Schur vectors.                              |
//|             Array whose indexes range within [1..N, 1..N].       |
//| Note 1:                                                          |
//|     The block structure of matrix T could be easily recognized:  |
//|     since all the elements below the blocks are zeros, the       |
//|     elements a[i+1,i] which are equal to 0 show the block border.|
//| Note 2:                                                          |
//|     the algorithm  performance  depends  on  the  value  of  the |
//|     internal parameter NS of InternalSchurDecomposition          |
//|     subroutine which  defines the number of shifts in the QR     |
//|     algorithm (analog of the block width in block matrix         |
//|     algorithms in linear algebra). If you require  maximum       |
//|     performance on your machine, it is recommended to            |
//|     adjust thisparameter manually.                               |
//| Result:                                                          |
//|     True, if the algorithm has converged and the parameters H and| 
//|         S contain the result.                                    |
//|     False, if the algorithm has not converged.                   |
//| Algorithm implemented on the basis of subroutine DHSEQR          |
//| (LAPACK 3.0 library).                                            |
//+------------------------------------------------------------------+
static bool CHsSchur::UpperHessenbergSchurDecomposition(CMatrixDouble &h,
                                                        const int n,
                                                        CMatrixDouble &s)
  {
//--- create variables
   bool result;
   int  info=0;
//--- create arrays
   double wi[];
   double wr[];
//--- function call
   InternalSchurDecomposition(h,n,1,2,wr,wi,s,info);
   result=info==0;
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CHsSchur::InternalSchurDecomposition(CMatrixDouble &h,const int n,
                                                 const int tneeded,
                                                 const int zneeded,
                                                 double &wr[],double &wi[],
                                                 CMatrixDouble &z,int &info)
  {
//--- create variables
   int    i=0;
   int    i1=0;
   int    i2=0;
   int    ierr=0;
   int    ii=0;
   int    itemp=0;
   int    itn=0;
   int    its=0;
   int    j=0;
   int    k=0;
   int    l=0;
   int    maxb=0;
   int    nr=0;
   int    ns=0;
   int    nv=0;
   double absw=0;
   double ovfl=0;
   double smlnum=0;
   double tau=0;
   double temp=0;
   double tst1=0;
   double ulp=0;
   double unfl=0;
   bool   initz;
   bool   wantt;
   bool   wantz;
   double cnst=0;
   bool   failflag;
   int    p1=0;
   int    p2=0;
   double vt=0;
   int    i_=0;
   int    i1_=0;
//--- create arrays
   double v[];
   double vv[];
   double work[];
   double workc1[];
   double works1[];
   double workv3[];
   double tmpwr[];
   double tmpwi[];
   CMatrixDouble s;
//--- initialization
   info=0;
//--- Set the order of the multi-shift QR algorithm to be used.
//--- If you want to tune algorithm, change this values
   ns=12;
   maxb=50;
//--- Now 2 < NS <= MAXB < NH.
   maxb=(int)MathMax(3,maxb);
   ns=MathMin(maxb,ns);
//--- Initialize
   cnst=1.5;
//--- allocation
   ArrayResizeAL(work,(int)MathMax(n,1)+1);
   s.Resize(ns+1,ns+1);
   ArrayResizeAL(v,ns+2);
   ArrayResizeAL(vv,ns+2);
   ArrayResizeAL(wr,(int)MathMax(n,1)+1);
   ArrayResizeAL(wi,(int)MathMax(n,1)+1);
   ArrayResizeAL(workc1,2);
   ArrayResizeAL(works1,2);
   ArrayResizeAL(workv3,4);
   ArrayResizeAL(tmpwr,(int)MathMax(n,1)+1);
   ArrayResizeAL(tmpwi,(int)MathMax(n,1)+1);
//--- check
   if(!CAp::Assert(n>=0,"InternalSchurDecomposition: incorrect N!"))
      return;
//--- check
   if(!CAp::Assert(tneeded==0 || tneeded==1,"InternalSchurDecomposition: incorrect TNeeded!"))
      return;
//--- check
   if(!CAp::Assert((zneeded==0 || zneeded==1) || zneeded==2,"InternalSchurDecomposition: incorrect ZNeeded!"))
      return;
//--- initialization
   wantt=tneeded==1;
   initz=zneeded==2;
   wantz=zneeded!=0;
   info=0;
//--- Initialize Z, if necessary
   if(initz)
     {
      z.Resize(n+1,n+1);
      for(i=1;i<=n;i++)
        {
         for(j=1;j<=n;j++)
           {
            //--- check
            if(i==j)
               z[i].Set(j,1);
            else
               z[i].Set(j,0);
           }
        }
     }
//--- Quick return if possible
   if(n==0)
      return;
//--- check
   if(n==1)
     {
      wr[1]=h[1][1];
      wi[1]=0;
      //--- exit the function
      return;
     }
//--- Set rows and columns 1 to N to zero below the first
//--- subdiagonal.
   for(j=1;j<=n-2;j++)
      for(i=j+2;i<=n;i++)
         h[i].Set(j,0);
//--- Test if N is sufficiently small
   if((ns<=2 || ns>n) || maxb>=n)
     {
      //--- Use the standard double-shift algorithm
      InternalAuxSchur(wantt,wantz,n,1,n,h,wr,wi,1,n,z,work,workv3,workc1,works1,info);
      //--- fill entries under diagonal blocks of T with zeros
      if(wantt)
        {
         j=1;
         while(j<=n)
           {
            //--- check
            if(wi[j]==0.0)
              {
               for(i=j+1;i<=n;i++)
                  h[i].Set(j,0);
               j=j+1;
              }
            else
              {
               for(i=j+2;i<=n;i++)
                 {
                  h[i].Set(j,0);
                  h[i].Set(j+1,0);
                 }
               j=j+2;
              }
           }
        }
      //--- exit the function
      return;
     }
//--- change values
   unfl=CMath::m_minrealnumber;
   ovfl=1/unfl;
   ulp=2*CMath::m_machineepsilon;
   smlnum=unfl*(n/ulp);
//--- I1 and I2 are the indices of the first row and last column of H
//--- to which transformations must be applied. If eigenvalues only are
//--- being computed, I1 and I2 are set inside the main loop.
   i1=1;
   i2=n;
//--- ITN is the total number of multiple-shift QR iterations allowed.
   itn=30*n;
//--- main block
   i=n;
   while(true)
     {
      l=1;
      //--- check
      if(i<1)
        {
         //--- fill entries under diagonal blocks of T with zeros
         if(wantt)
           {
            j=1;
            while(j<=n)
              {
               //--- check
               if(wi[j]==0.0)
                 {
                  for(i=j+1;i<=n;i++)
                     h[i].Set(j,0);
                  j=j+1;
                 }
               else
                 {
                  for(i=j+2;i<=n;i++)
                    {
                     h[i].Set(j,0);
                     h[i].Set(j+1,0);
                    }
                  j=j+2;
                 }
              }
           }
         //--- Exit
         return;
        }
      //--- Perform multiple-shift QR iterations on rows and columns ILO to I
      //--- until a submatrix of order at most MAXB splits off at the bottom
      //--- because a subdiagonal element has become negligible.
      failflag=true;
      for(its=0;its<=itn;its++)
        {
         //--- Look for a single small subdiagonal element.
         for(k=i;k>=l+1;k--)
           {
            tst1=MathAbs(h[k-1][k-1])+MathAbs(h[k][k]);
            //--- check
            if(tst1==0.0)
               tst1=CBlas::UpperHessenberg1Norm(h,l,i,l,i,work);
            //--- check
            if(MathAbs(h[k][k-1])<=MathMax(ulp*tst1,smlnum))
               break;
           }
         l=k;
         //--- check
         if(l>1)
            //--- H(L,L-1) is negligible.
            h[l].Set(l-1,0);
         //--- Exit from loop if a submatrix of order <= MAXB has split off.
         if(l>=i-maxb+1)
           {
            failflag=false;
            break;
           }
         //--- Now the active submatrix is in rows and columns L to I. If
         //--- eigenvalues only are being computed, only the active submatrix
         //--- need be transformed.
         if(its==20 || its==30)
           {
            //--- Exceptional shifts.
            for(ii=i-ns+1;ii<=i;ii++)
              {
               wr[ii]=cnst*(MathAbs(h[ii][ii-1])+MathAbs(h[ii][ii]));
               wi[ii]=0;
              }
           }
         else
           {
            //--- Use eigenvalues of trailing submatrix of order NS as shifts.
            CBlas::CopyMatrix(h,i-ns+1,i,i-ns+1,i,s,1,ns,1,ns);
            //--- function call
            InternalAuxSchur(false,false,ns,1,ns,s,tmpwr,tmpwi,1,ns,z,work,workv3,workc1,works1,ierr);
            for(p1=1;p1<=ns;p1++)
              {
               wr[i-ns+p1]=tmpwr[p1];
               wi[i-ns+p1]=tmpwi[p1];
              }
            //--- check
            if(ierr>0)
              {
               //--- If DLAHQR failed to compute all NS eigenvalues, use the
               //--- unconverged diagonal elements as the remaining shifts.
               for(ii=1;ii<=ierr;ii++)
                 {
                  wr[i-ns+ii]=s[ii][ii];
                  wi[i-ns+ii]=0;
                 }
              }
           }
         //--- Form the first column of (G-w(1)) (G-w(2)) . . . (G-w(ns))
         //--- where G is the Hessenberg submatrix H(L:I,L:I) and w is
         //--- the vector of shifts (stored in WR and WI). The result is
         //--- stored in the local array V.
         v[1]=1;
         for(ii=2;ii<=ns+1;ii++)
            v[ii]=0;
         nv=1;
         for(j=i-ns+1;j<=i;j++)
           {
            //--- check
            if(wi[j]>=0.0)
              {
               //--- check
               if(wi[j]==0.0)
                 {
                  //--- real shift
                  p1=nv+1;
                  for(i_=1;i_<=p1;i_++)
                     vv[i_]=v[i_];
                  //--- function call
                  CBlas::MatrixVectorMultiply(h,l,l+nv,l,l+nv-1,false,vv,1,nv,1.0,v,1,nv+1,-wr[j]);
                  nv=nv+1;
                 }
               else
                 {
                  //--- check
                  if(wi[j]>0.0)
                    {
                     //--- complex conjugate pair of shifts
                     p1=nv+1;
                     for(i_=1;i_<=p1;i_++)
                        vv[i_]=v[i_];
                     //--- function call
                     CBlas::MatrixVectorMultiply(h,l,l+nv,l,l+nv-1,false,v,1,nv,1.0,vv,1,nv+1,-(2*wr[j]));
                     //--- function call
                     itemp=CBlas::VectorIdxAbsMax(vv,1,nv+1);
                     temp=1/MathMax(MathAbs(vv[itemp]),smlnum);
                     p1=nv+1;
                     for(i_=1;i_<=p1;i_++)
                        vv[i_]=temp*vv[i_];
                     //--- function call
                     absw=CBlas::PyThag2(wr[j],wi[j]);
                     temp=temp*absw*absw;
                     //--- function call
                     CBlas::MatrixVectorMultiply(h,l,l+nv+1,l,l+nv,false,vv,1,nv+1,1.0,v,1,nv+2,temp);
                     nv=nv+2;
                    }
                 }
               //--- Scale V(1:NV) so that max(abs(V(i))) = 1. If V is zero,
               //--- reset it to the unit vector.
               itemp=CBlas::VectorIdxAbsMax(v,1,nv);
               temp=MathAbs(v[itemp]);
               //--- check
               if(temp==0.0)
                 {
                  v[1]=1;
                  for(ii=2;ii<=nv;ii++)
                     v[ii]=0;
                 }
               else
                 {
                  temp=MathMax(temp,smlnum);
                  vt=1/temp;
                  for(i_=1;i_<=nv;i_++)
                     v[i_]=vt*v[i_];
                 }
              }
           }
         //--- Multiple-shift QR step
         for(k=l;k<=i-1;k++)
           {
            //--- The first iteration of this loop determines a reflection G
            //--- from the vector V and applies it from left and right to H,
            //--- thus creating a nonzero bulge below the subdiagonal.
            //--- Each subsequent iteration determines a reflection G to
            //--- restore the Hessenberg form in the (K-1)th column, and thus
            //--- chases the bulge one step toward the bottom of the active
            //--- submatrix. NR is the order of G.
            nr=MathMin(ns+1,i-k+1);
            //--- check
            if(k>l)
              {
               //--- change values
               p1=k-1;
               p2=k+nr-1;
               i1_=k-1;
               for(i_=1;i_<=nr;i_++)
                  v[i_]=h[i_+i1_][p1];
              }
            //--- function call
            CReflections::GenerateReflection(v,nr,tau);
            //--- check
            if(k>l)
              {
               h[k].Set(k-1,v[1]);
               for(ii=k+1;ii<=i;ii++)
                  h[ii].Set(k-1,0);
              }
            v[1]=1;
            //--- Apply G from the left to transform the rows of the matrix in
            //--- columns K to I2.
            CReflections::ApplyReflectionFromTheLeft(h,tau,v,k,k+nr-1,k,i2,work);
            //--- Apply G from the right to transform the columns of the
            //--- matrix in rows I1 to min(K+NR,I).
            CReflections::ApplyReflectionFromTheRight(h,tau,v,i1,MathMin(k+nr,i),k,k+nr-1,work);
            //--- check
            if(wantz)
               //--- Accumulate transformations in the matrix Z
               CReflections::ApplyReflectionFromTheRight(z,tau,v,1,n,k,k+nr-1,work);
           }
        }
      //--- Failure to converge in remaining number of iterations
      if(failflag)
        {
         info=i;
         //--- exit the function
         return;
        }
      //--- A submatrix of order <= MAXB in rows and columns L to I has split
      //--- off. Use the double-shift QR algorithm to handle it.
      InternalAuxSchur(wantt,wantz,n,l,i,h,wr,wi,1,n,z,work,workv3,workc1,works1,info);
      //--- check
      if(info>0)
         return;
      //--- Decrement number of remaining iterations, and return to start of
      //--- the main loop with a new value of I.
      itn=itn-its;
      i=l-1;
     }
  }
//+------------------------------------------------------------------+
//| The auxiliary function                                          
//+------------------------------------------------------------------+
static void CHsSchur::InternalAuxSchur(const bool wantt,const bool wantz,
                                       const int n,const int ilo,const int ihi,
                                       CMatrixDouble &h,double &wr[],double &wi[],
                                       const int iloz,const int ihiz,CMatrixDouble &z,
                                       double &work[],double &workv3[],
                                       double &workc1[],double &works1[],int &info)
  {
//--- create variables
   int    i=0;
   int    i1=0;
   int    i2=0;
   int    itn=0;
   int    its=0;
   int    j=0;
   int    k=0;
   int    l=0;
   int    m=0;
   int    nh=0;
   int    nr=0;
   int    nz=0;
   double ave=0;
   double cs=0;
   double disc=0;
   double h00=0;
   double h10=0;
   double h11=0;
   double h12=0;
   double h21=0;
   double h22=0;
   double h33=0;
   double h33s=0;
   double h43h34=0;
   double h44=0;
   double h44s=0;
   double ovfl=0;
   double s=0;
   double smlnum=0;
   double sn=0;
   double sum=0;
   double t1=0;
   double t2=0;
   double t3=0;
   double tst1=0;
   double unfl=0;
   double v1=0;
   double v2=0;
   double v3=0;
   bool   failflag;
   double dat1=0;
   double dat2=0;
   int    p1=0;
   double him1im1=0;
   double him1i=0;
   double hiim1=0;
   double hii=0;
   double wrim1=0;
   double wri=0;
   double wiim1=0;
   double wii=0;
   double ulp=0;
//--- initialization
   info=0;
   dat1=0.75;
   dat2=-0.4375;
   ulp=CMath::m_machineepsilon;
//--- Quick return if possible
   if(n==0)
      return;
//--- check
   if(ilo==ihi)
     {
      wr[ilo]=h[ilo][ilo];
      wi[ilo]=0;
      //--- exit the function
      return;
     }
//--- initialization
   nh=ihi-ilo+1;
   nz=ihiz-iloz+1;
//--- Set machine-dependent constants for the stopping criterion.
//--- If norm(H) <= sqrt(OVFL), overflow should not occur.
   unfl=CMath::m_minrealnumber;
   ovfl=1/unfl;
   smlnum=unfl*(nh/ulp);
//--- I1 and I2 are the indices of the first row and last column of H
//--- to which transformations must be applied. If eigenvalues only are
//--- being computed, I1 and I2 are set inside the main loop.
   i1=1;
   i2=n;
//--- ITN is the total number of QR iterations allowed.
   itn=30*nh;
//--- The main loop begins here. I is the loop index and decreases from
//--- IHI to ILO in steps of 1 or 2. Each iteration of the loop works
//--- with the active submatrix in rows and columns L to I.
//--- Eigenvalues I+1 to IHI have already converged. Either L = ILO or
//--- H(L,L-1) is negligible so that the matrix splits.
   i=ihi;
   while(true)
     {
      l=ilo;
      //--- check
      if(i<ilo)
         return;
      //--- Perform QR iterations on rows and columns ILO to I until a
      //--- submatrix of order 1 or 2 splits off at the bottom because a
      //--- subdiagonal element has become negligible.
      failflag=true;
      for(its=0;its<=itn;its++)
        {
         //--- Look for a single small subdiagonal element.
         for(k=i;k>=l+1;k--)
           {
            tst1=MathAbs(h[k-1][k-1])+MathAbs(h[k][k]);
            //--- check
            if(tst1==0.0)
               tst1=CBlas::UpperHessenberg1Norm(h,l,i,l,i,work);
            //--- check
            if(MathAbs(h[k][k-1])<=MathMax(ulp*tst1,smlnum))
               break;
           }
         l=k;
         //--- check
         if(l>ilo)
            //--- H(L,L-1) is negligible
            h[l].Set(l-1,0);
         //--- Exit from loop if a submatrix of order 1 or 2 has split off.
         if(l>=i-1)
           {
            failflag=false;
            break;
           }
         //--- Now the active submatrix is in rows and columns L to I. If
         //--- eigenvalues only are being computed, only the active submatrix
         //--- need be transformed.
         if(its==10 || its==20)
           {
            //--- Exceptional shift.
            s=MathAbs(h[i][i-1])+MathAbs(h[i-1][i-2]);
            h44=dat1*s+h[i][i];
            h33=h44;
            h43h34=dat2*s*s;
           }
         else
           {
            //--- Prepare to use Francis' double shift
            //--- (i.e. 2nd degree generalized Rayleigh quotient)
            h44=h[i][i];
            h33=h[i-1][i-1];
            h43h34=h[i][i-1]*h[i-1][i];
            s=h[i-1][i-2]*h[i-1][i-2];
            disc=(h33-h44)*0.5;
            disc=disc*disc+h43h34;
            //--- check
            if(disc>0.0)
              {
               //--- Real roots: use Wilkinson's shift twice
               disc=MathSqrt(disc);
               ave=0.5*(h33+h44);
               //--- check
               if(MathAbs(h33)-MathAbs(h44)>0.0)
                 {
                  h33=h33*h44-h43h34;
                  h44=h33/(ExtSchurSign(disc,ave)+ave);
                 }
               else
                  h44=ExtSchurSign(disc,ave)+ave;
               h33=h44;
               h43h34=0;
              }
           }
         //--- Look for two consecutive small subdiagonal elements.
         for(m=i-2;m>=l;m--)
           {
            //--- Determine the effect of starting the double-shift QR
            //--- iteration at row M, and see if this would make H(M,M-1)
            //--- negligible.
            h11=h[m][m];
            h22=h[m+1][m+1];
            h21=h[m+1][m];
            h12=h[m][m+1];
            h44s=h44-h11;
            h33s=h33-h11;
            v1=(h33s*h44s-h43h34)/h21+h12;
            v2=h22-h11-h33s-h44s;
            v3=h[m+2][m+1];
            s=MathAbs(v1)+MathAbs(v2)+MathAbs(v3);
            v1=v1/s;
            v2=v2/s;
            v3=v3/s;
            workv3[1]=v1;
            workv3[2]=v2;
            workv3[3]=v3;
            //--- check
            if(m==l)
               break;
            h00=h[m-1][m-1];
            h10=h[m][m-1];
            tst1=MathAbs(v1)*(MathAbs(h00)+MathAbs(h11)+MathAbs(h22));
            //--- check
            if(MathAbs(h10)*(MathAbs(v2)+MathAbs(v3))<=ulp*tst1)
               break;
           }
         //--- Double-shift QR step
         for(k=m;k<=i-1;k++)
           {
            //--- The first iteration of this loop determines a reflection G
            //--- from the vector V and applies it from left and right to H,
            //--- thus creating a nonzero bulge below the subdiagonal.
            //--- Each subsequent iteration determines a reflection G to
            //--- restore the Hessenberg form in the (K-1)th column, and thus
            //--- chases the bulge one step toward the bottom of the active
            //--- submatrix. NR is the order of G.
            nr=(int)MathMin(3,i-k+1);
            //--- check
            if(k>m)
              {
               for(p1=1;p1<=nr;p1++)
                  workv3[p1]=h[k+p1-1][k-1];
              }
            //--- function call
            CReflections::GenerateReflection(workv3,nr,t1);
            //--- check
            if(k>m)
              {
               h[k].Set(k-1,workv3[1]);
               h[k+1].Set(k-1,0);
               //--- check
               if(k<i-1)
                  h[k+2].Set(k-1,0);
              }
            else
              {
               //--- check
               if(m>l)
                  h[k].Set(k-1,-h[k][k-1]);
              }
            v2=workv3[2];
            t2=t1*v2;
            //--- check
            if(nr==3)
              {
               v3=workv3[3];
               t3=t1*v3;
               //--- Apply G from the left to transform the rows of the matrix
               //--- in columns K to I2.
               for(j=k;j<=i2;j++)
                 {
                  sum=h[k][j]+v2*h[k+1][j]+v3*h[k+2][j];
                  h[k].Set(j,h[k][j]-sum*t1);
                  h[k+1].Set(j,h[k+1][j]-sum*t2);
                  h[k+2].Set(j,h[k+2][j]-sum*t3);
                 }
               //--- Apply G from the right to transform the columns of the
               //--- matrix in rows I1 to min(K+3,I).
               for(j=i1;j<=MathMin(k+3,i);j++)
                 {
                  sum=h[j][k]+v2*h[j][k+1]+v3*h[j][k+2];
                  h[j].Set(k,h[j][k]-sum*t1);
                  h[j].Set(k+1,h[j][k+1]-sum*t2);
                  h[j].Set(k+2,h[j][k+2]-sum*t3);
                 }
               //--- check
               if(wantz)
                 {
                  //--- Accumulate transformations in the matrix Z
                  for(j=iloz;j<=ihiz;j++)
                    {
                     sum=z[j][k]+v2*z[j][k+1]+v3*z[j][k+2];
                     z[j].Set(k,z[j][k]-sum*t1);
                     z[j].Set(k+1,z[j][k+1]-sum*t2);
                     z[j].Set(k+2,z[j][k+2]-sum*t3);
                    }
                 }
              }
            else
              {
               //--- check
               if(nr==2)
                 {
                  //--- Apply G from the left to transform the rows of the matrix
                  //--- in columns K to I2.
                  for(j=k;j<=i2;j++)
                    {
                     sum=h[k][j]+v2*h[k+1][j];
                     h[k].Set(j,h[k][j]-sum*t1);
                     h[k+1].Set(j,h[k+1][j]-sum*t2);
                    }
                  //--- Apply G from the right to transform the columns of the
                  //--- matrix in rows I1 to min(K+3,I).
                  for(j=i1;j<=i;j++)
                    {
                     sum=h[j][k]+v2*h[j][k+1];
                     h[j].Set(k,h[j][k]-sum*t1);
                     h[j].Set(k+1,h[j][k+1]-sum*t2);
                    }
                  //--- check
                  if(wantz)
                    {
                     //--- Accumulate transformations in the matrix Z
                     for(j=iloz;j<=ihiz;j++)
                       {
                        sum=z[j][k]+v2*z[j][k+1];
                        z[j].Set(k,z[j][k]-sum*t1);
                        z[j].Set(k+1,z[j][k+1]-sum*t2);
                       }
                    }
                 }
              }
           }
        }
      //--- check
      if(failflag)
        {
         //--- Failure to converge in remaining number of iterations
         info=i;
         //--- exit the function
         return;
        }
      //--- check
      if(l==i)
        {
         //--- H(I,I-1) is negligible: one eigenvalue has converged.
         wr[i]=h[i][i];
         wi[i]=0;
        }
      else
        {
         //--- check
         if(l==i-1)
           {
            //--- H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
            //--- Transform the 2-by-2 submatrix to standard Schur form,
            //--- and compute and store the eigenvalues.
            him1im1=h[i-1][i-1];
            him1i=h[i-1][i];
            hiim1=h[i][i-1];
            hii=h[i][i];
            //--- function call
            Aux2x2Schur(him1im1,him1i,hiim1,hii,wrim1,wiim1,wri,wii,cs,sn);
            //--- change values
            wr[i-1]=wrim1;
            wi[i-1]=wiim1;
            wr[i]=wri;
            wi[i]=wii;
            h[i-1].Set(i-1,him1im1);
            h[i-1].Set(i,him1i);
            h[i].Set(i-1,hiim1);
            h[i].Set(i,hii);
            //--- check
            if(wantt)
              {
               //--- Apply the transformation to the rest of H.
               if(i2>i)
                 {
                  workc1[1]=cs;
                  works1[1]=sn;
                  //--- function call
                  CRotations::ApplyRotationsFromTheLeft(true,i-1,i,i+1,i2,workc1,works1,h,work);
                 }
               workc1[1]=cs;
               works1[1]=sn;
               //--- function call
               CRotations::ApplyRotationsFromTheRight(true,i1,i-2,i-1,i,workc1,works1,h,work);
              }
            //--- check
            if(wantz)
              {
               //--- Apply the transformation to Z.
               workc1[1]=cs;
               works1[1]=sn;
               //--- function call
               CRotations::ApplyRotationsFromTheRight(true,iloz,iloz+nz-1,i-1,i,workc1,works1,z,work);
              }
           }
        }
      //--- Decrement number of remaining iterations, and return to start of
      //--- the main loop with new value of I.
      itn=itn-its;
      i=l-1;
     }
  }
//+------------------------------------------------------------------+
//| 2x2                                                              |
//+------------------------------------------------------------------+
static void CHsSchur::Aux2x2Schur(double &a,double &b,double &c,double &d,
                                  double &rt1r,double &rt1i,double &rt2r,
                                  double &rt2i,double &cs,double &sn)
  {
//--- create variables
   double multpl=0;
   double aa=0;
   double bb=0;
   double bcmax=0;
   double bcmis=0;
   double cc=0;
   double cs1=0;
   double dd=0;
   double eps=0;
   double p=0;
   double sab=0;
   double sac=0;
   double scl=0;
   double sigma=0;
   double sn1=0;
   double tau=0;
   double temp=0;
   double z=0;
//--- initialization
   rt1r=0;
   rt1i=0;
   rt2r=0;
   rt2i=0;
   cs=0;
   sn=0;
   multpl=4.0;
   eps=CMath::m_machineepsilon;
//--- check
   if(c==0.0)
     {
      cs=1;
      sn=0;
     }
   else
     {
      //--- check
      if(b==0.0)
        {
         //--- Swap rows and columns
         cs=0;
         sn=1;
         temp=d;
         d=a;
         a=temp;
         b=-c;
         c=0;
        }
      else
        {
         //--- check
         if(a-d==0.0 && ExtSchurSignToone(b)!=ExtSchurSignToone(c))
           {
            cs=1;
            sn=0;
           }
         else
           {
            //--- change values
            temp=a-d;
            p=0.5*temp;
            bcmax=MathMax(MathAbs(b),MathAbs(c));
            bcmis=MathMin(MathAbs(b),MathAbs(c))*ExtSchurSignToone(b)*ExtSchurSignToone(c);
            scl=MathMax(MathAbs(p),bcmax);
            z=p/scl*p+bcmax/scl*bcmis;
            //--- If Z is of the order of the machine accuracy, postpone the
            //--- decision on the nature of eigenvalues
            if(z>=multpl*eps)
              {
               //--- Real eigenvalues. Compute A and D.
               z=p+ExtSchurSign(MathSqrt(scl)*MathSqrt(z),p);
               a=d+z;
               d=d-bcmax/z*bcmis;
               //--- Compute B and the rotation matrix
               tau=CBlas::PyThag2(c,z);
               cs=z/tau;
               sn=c/tau;
               b=b-c;
               c=0;
              }
            else
              {
               //--- Complex eigenvalues, or real (almost) equal eigenvalues.
               //--- Make diagonal elements equal.
               sigma=b+c;
               tau=CBlas::PyThag2(sigma,temp);
               cs=MathSqrt(0.5*(1+MathAbs(sigma)/tau));
               sn=-(p/(tau*cs)*ExtSchurSign(1,sigma));
               //--- Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
               //---         [ CC  DD ]   [ C  D ] [ SN  CS ]
               aa=a*cs+b*sn;
               bb=-(a*sn)+b*cs;
               cc=c*cs+d*sn;
               dd=-(c*sn)+d*cs;
               //--- Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
               //---         [ C  D ]   [-SN  CS ] [ CC  DD ]
               a=aa*cs+cc*sn;
               b=bb*cs+dd*sn;
               c=-(aa*sn)+cc*cs;
               d=-(bb*sn)+dd*cs;
               temp=0.5*(a+d);
               a=temp;
               d=temp;
               //--- check
               if(c!=0.0)
                 {
                  //--- check
                  if(b!=0.0)
                    {
                     //--- check
                     if(ExtSchurSignToone(b)==ExtSchurSignToone(c))
                       {
                        //--- Real eigenvalues: reduce to upper triangular form
                        sab=MathSqrt(MathAbs(b));
                        sac=MathSqrt(MathAbs(c));
                        //--- function call
                        p=ExtSchurSign(sab*sac,c);
                        tau=1/MathSqrt(MathAbs(b+c));
                        a=temp+p;
                        d=temp-p;
                        b=b-c;
                        c=0;
                        cs1=sab*tau;
                        sn1=sac*tau;
                        temp=cs*cs1-sn*sn1;
                        sn=cs*sn1+sn*cs1;
                        cs=temp;
                       }
                    }
                  else
                    {
                     //--- change values
                     b=-c;
                     c=0;
                     temp=cs;
                     cs=-sn;
                     sn=temp;
                    }
                 }
              }
           }
        }
     }
//--- Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
   rt1r=a;
   rt2r=d;
//--- check
   if(c==0.0)
     {
      rt1i=0;
      rt2i=0;
     }
   else
     {
      rt1i=MathSqrt(MathAbs(b))*MathSqrt(MathAbs(c));
      rt2i=-rt1i;
     }
  }
//+------------------------------------------------------------------+
//| Schur sign                                                       |
//+------------------------------------------------------------------+
static double CHsSchur::ExtSchurSign(const double a,const double b)
  {
//--- create variables
   double result=0;
//--- check
   if(b>=0.0)
      result=MathAbs(a);
   else
      result=-MathAbs(a);
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Schur sign (1 or -1)                                             |
//+------------------------------------------------------------------+
static int CHsSchur::ExtSchurSignToone(const double b)
  {
//--- create variables
   int result=0;
//--- check
   if(b>=0.0)
      result=1;
   else
      result=-1;
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Safe solutions for tridiagonal linear matrix                     |
//+------------------------------------------------------------------+
class CTrLinSolve
  {
public:
   //--- constructor, destructor
                     CTrLinSolve(void);
                    ~CTrLinSolve(void);
   //--- methods
   static void       RMatrixTrSafeSolve(CMatrixDouble &a,const int n,double &x[],double &s,const bool isupper,const bool istrans,const bool isunit);
   static void       SafeSolveTriangular(CMatrixDouble &a,const int n,double &x[],double &s,const bool isupper,const bool istrans,const bool isunit,const bool normin,double &cnorm[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CTrLinSolve::CTrLinSolve(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CTrLinSolve::~CTrLinSolve(void)
  {

  }
//+------------------------------------------------------------------+
//| Utility subroutine performing the "safe" solution of system of   |
//| linear equations with triangular coefficient matrices.           |
//| The subroutine uses scaling and solves the scaled system A*x=s*b | 
//| (where  s is  a  scalar  value)  instead  of  A*x=b,  choosing   |
//| s  so  that x can be represented by a floating-point number. The | 
//| closer the system  gets  to  a  singular, the less s is. If the  |
//| system is singular, s=0 and x contains the non-trivial solution  |
//| of equation A*x=0.                                               |
//| The feature of an algorithm is that it could not cause an        |
//| overflow  or  a division by zero regardless of the matrix used   |
//| as the input.                                                    |
//| The algorithm can solve systems of equations with  upper/lower   | 
//| triangular matrices,  with/without unit diagonal, and systems of | 
//| type A*x=b or A'*x=b (where A' is a transposed matrix A).        |
//| Input parameters:                                                |
//|     A       -   system matrix. Array whose indexes range within  |
//|                 [0..N-1, 0..N-1].                                |
//|     N       -   size of matrix A.                                |
//|     X       -   right-hand member of a system.                   |
//|                 Array whose index ranges within [0..N-1].        |
//|     IsUpper -   matrix type. If it is True, the system matrix is | 
//|                 the upper triangular and is located in  the      |
//|                 corresponding  part  of matrix A.                |
//|     Trans   -   problem type. If it is True, the problem to be   | 
//|                 solved  is A'*x=b, otherwise it is A*x=b.        |
//|     Isunit  -   matrix type. If it is True, the system matrix has|  
//|                 a  unit diagonal (the elements on the main       |
//|                 diagonal are  not  used in the calculation       |
//|                 process), otherwise the matrix is considered to  |
//|                 be a general triangular matrix.                  |
//| Output parameters:                                               |
//|     X       -   solution. Array whose index ranges within        |
//|                 [0..N-1].                                        |
//|     S       -   scaling factor.                                  |
//|   -- LAPACK auxiliary routine (version 3.0) --                   |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      June 30, 1992                                               |
//+------------------------------------------------------------------+
static void CTrLinSolve::RMatrixTrSafeSolve(CMatrixDouble &a,const int n,
                                            double &x[],double &s,
                                            const bool isupper,
                                            const bool istrans,
                                            const bool isunit)
  {
//--- create variables
   bool normin;
   int  i=0;
   int  i_=0;
   int  i1_=0;
//--- create arrays
   double cnorm[];
   double x1[];
//--- create matrix
   CMatrixDouble a1;
//--- initialization
   s=0;
//--- From 0-based to 1-based
   normin=false;
//--- allocation
   a1.Resize(n+1,n+1);
   ArrayResizeAL(x1,n+1);
   for(i=1;i<=n;i++)
     {
      i1_=-1;
      for(i_=1;i_<=n;i_++)
         a1[i].Set(i_,a[i-1][i_+i1_]);
     }
//--- change values
   i1_=-1;
   for(i_=1;i_<=n;i_++)
      x1[i_]=x[i_+i1_];
//--- Solve 1-based
   SafeSolveTriangular(a1,n,x1,s,isupper,istrans,isunit,normin,cnorm);
//--- From 1-based to 0-based
   i1_=1;
   for(i_=0;i_<=n-1;i_++)
      x[i_]=x1[i_+i1_];
  }
//+------------------------------------------------------------------+
//| Obsolete 1-based subroutine.                                     |
//| See RMatrixTRSafeSolve for 0-based replacement.                  |
//+------------------------------------------------------------------+
static void CTrLinSolve::SafeSolveTriangular(CMatrixDouble &a,const int n,
                                             double &x[],double &s,
                                             const bool isupper,
                                             const bool istrans,
                                             const bool isunit,
                                             const bool normin,
                                             double &cnorm[])
  {
//--- create variables
   int    i=0;
   int    imax=0;
   int    j=0;
   int    jfirst=0;
   int    jinc=0;
   int    jlast=0;
   int    jm1=0;
   int    jp1=0;
   int    ip1=0;
   int    im1=0;
   int    k=0;
   int    flg=0;
   double v=0;
   double vd=0;
   double bignum=0;
   double grow=0;
   double rec=0;
   double smlnum=0;
   double sumj=0;
   double tjj=0;
   double tjjs=0;
   double tmax=0;
   double tscal=0;
   double uscal=0;
   double xbnd=0;
   double xj=0;
   double xmax=0;
   bool   notran;
   bool   upper;
   bool   nounit;
   int    i_=0;
//--- initialization
   s=0;
   upper=isupper;
   notran=!istrans;
   nounit=!isunit;
//--- these initializers are not really necessary,
//--- but without them compiler complains about uninitialized locals
   tjjs=0;
//--- Quick return if possible
   if(n==0)
      return;
//--- Determine machine dependent parameters to control overflow.
   smlnum=CMath::m_minrealnumber/(CMath::m_machineepsilon*2);
   bignum=1/smlnum;
   s=1;
//--- check
   if(!normin)
     {
      ArrayResizeAL(cnorm,n+1);
      //--- Compute the 1-norm of each column,not including the diagonal.
      if(upper)
        {
         //--- A is upper triangular.
         for(j=1;j<=n;j++)
           {
            v=0;
            for(k=1;k<=j-1;k++)
               v=v+MathAbs(a[k][j]);
            cnorm[j]=v;
           }
        }
      else
        {
         //--- A is lower triangular.
         for(j=1;j<=n-1;j++)
           {
            v=0;
            for(k=j+1;k<=n;k++)
               v=v+MathAbs(a[k][j]);
            cnorm[j]=v;
           }
         cnorm[n]=0;
        }
     }
//--- Scale the column norms by TSCAL if the maximum element in CNORM is
//--- greater than BIGNUM.
   imax=1;
   for(k=2;k<=n;k++)
     {
      //--- check
      if(cnorm[k]>cnorm[imax])
         imax=k;
     }
   tmax=cnorm[imax];
//--- check
   if(tmax<=bignum)
      tscal=1;
   else
     {
      tscal=1/(smlnum*tmax);
      for(i_=1;i_<=n;i_++)
         cnorm[i_]=tscal*cnorm[i_];
     }
//--- Compute a bound on the computed solution vector to see if the
//--- Level 2 BLAS routine DTRSV can be used.
   j=1;
   for(k=2;k<=n;k++)
     {
      //--- check
      if(MathAbs(x[k])>MathAbs(x[j]))
         j=k;
     }
//--- change values
   xmax=MathAbs(x[j]);
   xbnd=xmax;
//--- check
   if(notran)
     {
      //--- Compute the growth in A * x=b.
      if(upper)
        {
         jfirst=n;
         jlast=1;
         jinc=-1;
        }
      else
        {
         jfirst=1;
         jlast=n;
         jinc=1;
        }
      //--- check
      if(tscal!=1.0)
         grow=0;
      else
        {
         //--- check
         if(nounit)
           {
            //--- A is non-unit triangular.
            //--- Compute GROW=1/G(j) and XBND=1/M(j).
            //--- Initially,G(0)=max{x(i),i=1,...,n}.
            grow=1/MathMax(xbnd,smlnum);
            xbnd=grow;
            j=jfirst;
            while((jinc>0 && j<=jlast) || (jinc<0 && j>=jlast))
              {
               //--- Exit the loop if the growth factor is too small.
               if(grow<=smlnum)
                  break;
               //--- M(j)=G(j-1) / abs(A(j,j))
               tjj=MathAbs(a[j][j]);
               xbnd=MathMin(xbnd,MathMin(1,tjj)*grow);
               //--- check
               if(tjj+cnorm[j]>=smlnum)
                 {
                  //--- G(j)=G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
                  grow=grow*(tjj/(tjj+cnorm[j]));
                 }
               else
                 {
                  //--- G(j) could overflow,set GROW to 0.
                  grow=0;
                 }
               //--- check
               if(j==jlast)
                  grow=xbnd;
               j=j+jinc;
              }
           }
         else
           {
            //--- A is unit triangular.
            //--- Compute GROW=1/G(j), where G(0)=max{x(i), i=1,...,n}.
            grow=MathMin(1,1/MathMax(xbnd,smlnum));
            j=jfirst;
            while((jinc>0 && j<=jlast) || (jinc<0 && j>=jlast))
              {
               //--- Exit the loop if the growth factor is too small.
               if(grow<=smlnum)
                  break;
               //--- G(j) = G(j-1)*( 1 + CNORM(j) )
               grow=grow*(1/(1+cnorm[j]));
               j=j+jinc;
              }
           }
        }
     }
   else
     {
      //--- Compute the growth in A' * x = b.
      if(upper)
        {
         jfirst=1;
         jlast=n;
         jinc=1;
        }
      else
        {
         jfirst=n;
         jlast=1;
         jinc=-1;
        }
      //--- check
      if(tscal!=1.0)
         grow=0;
      else
        {
         //--- check
         if(nounit)
           {
            //--- A is non-unit triangular.
            //--- Compute GROW=1/G(j) and XBND=1/M(j).
            //--- Initially, M(0)=max{x(i), i=1,...,n}.
            grow=1/MathMax(xbnd,smlnum);
            xbnd=grow;
            j=jfirst;
            while((jinc>0 && j<=jlast) || (jinc<0 && j>=jlast))
              {
               //--- Exit the loop if the growth factor is too small.
               if(grow<=smlnum)
                  break;
               //--- G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
               xj=1+cnorm[j];
               grow=MathMin(grow,xbnd/xj);
               //--- M(j)=M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
               tjj=MathAbs(a[j][j]);
               //--- check
               if(xj>tjj)
                  xbnd=xbnd*(tjj/xj);
               //--- check
               if(j==jlast)
                  grow=MathMin(grow,xbnd);
               j=j+jinc;
              }
           }
         else
           {
            //--- A is unit triangular.
            //--- Compute GROW=1/G(j), where G(0)=max{x(i), i=1,...,n}.
            grow=MathMin(1,1/MathMax(xbnd,smlnum));
            j=jfirst;
            while((jinc>0 && j<=jlast) || (jinc<0 && j>=jlast))
              {
               //--- Exit the loop if the growth factor is too small.
               if(grow<=smlnum)
                  break;
               //--- G(j)=( 1 + CNORM(j) )*G(j-1)
               xj=1+cnorm[j];
               grow=grow/xj;
               j=j+jinc;
              }
           }
        }
     }
   if(grow*tscal>smlnum)
     {
      //--- Use the Level 2 BLAS solve if the reciprocal of the bound on
      //--- elements of X is not too small.
      if((upper && notran) || (!upper && !notran))
        {
         //--- check
         if(nounit)
            vd=a[n][n];
         else
            vd=1;
         x[n]=x[n]/vd;
         for(i=n-1;i>=1;i--)
           {
            ip1=i+1;
            //--- check
            if(upper)
              {
               v=0.0;
               for(i_=ip1;i_<=n;i_++)
                  v+=a[i][i_]*x[i_];
              }
            else
              {
               v=0.0;
               for(i_=ip1;i_<=n;i_++)
                  v+=a[i_][i]*x[i_];
              }
            //--- check
            if(nounit)
               vd=a[i][i];
            else
               vd=1;
            x[i]=(x[i]-v)/vd;
           }
        }
      else
        {
         //--- check
         if(nounit)
            vd=a[1][1];
         else
            vd=1;
         x[1]=x[1]/vd;
         for(i=2;i<=n;i++)
           {
            im1=i-1;
            //--- check
            if(upper)
              {
               v=0.0;
               for(i_=1;i_<=im1;i_++)
                  v+=a[i_][i]*x[i_];
              }
            else
              {
               v=0.0;
               for(i_=1;i_<=im1;i_++)
                  v+=a[i][i_]*x[i_];
              }
            //--- check
            if(nounit)
               vd=a[i][i];
            else
               vd=1;
            x[i]=(x[i]-v)/vd;
           }
        }
     }
   else
     {
      //--- Use a Level 1 BLAS solve, scaling intermediate results.
      if(xmax>bignum)
        {
         //--- Scale X so that its components are less than or equal to
         //--- BIGNUM in absolute value.
         s=bignum/xmax;
         for(i_=1;i_<=n;i_++)
            x[i_]=s*x[i_];
         xmax=bignum;
        }
      //--- check
      if(notran)
        {
         //--- Solve A * x = b
         j=jfirst;
         while((jinc>0 && j<=jlast) || (jinc<0 && j>=jlast))
           {
            //--- Compute x(j)=b(j) / A(j,j), scaling x if necessary.
            xj=MathAbs(x[j]);
            flg=0;
            //--- check
            if(nounit)
               tjjs=a[j][j]*tscal;
            else
              {
               tjjs=tscal;
               //--- check
               if(tscal==1.0)
                  flg=100;
              }
            //--- check
            if(flg!=100)
              {
               tjj=MathAbs(tjjs);
               if(tjj>smlnum)
                 {
                  //--- abs(A(j,j)) > SMLNUM:
                  if(tjj<1.0)
                    {
                     //--- check
                     if(xj>(double)(tjj*bignum))
                       {
                        //--- Scale x by 1/b(j).
                        rec=1/xj;
                        for(i_=1;i_<=n;i_++)
                          {
                           x[i_]=rec*x[i_];
                          }
                        s=s*rec;
                        xmax=xmax*rec;
                       }
                    }
                  x[j]=x[j]/tjjs;
                  xj=MathAbs(x[j]);
                 }
               else
                 {
                  //--- check
                  if(tjj>0.0)
                    {
                     //--- 0 < abs(A(j,j)) <=SMLNUM:
                     if(xj>(double)(tjj*bignum))
                       {
                        //--- Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
                        //--- to avoid overflow when dividing by A(j,j).
                        rec=tjj*bignum/xj;
                        //--- check
                        if(cnorm[j]>1.0)
                          {
                           //--- Scale by 1/CNORM(j) to avoid overflow when
                           //--- multiplying x(j) times column j.
                           rec=rec/cnorm[j];
                          }
                        for(i_=1;i_<=n;i_++)
                           x[i_]=rec*x[i_];
                        s=s*rec;
                        xmax=xmax*rec;
                       }
                     x[j]=x[j]/tjjs;
                     xj=MathAbs(x[j]);
                    }
                  else
                    {
                     //--- A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                     //--- scale = 0, and compute a solution to A*x = 0.
                     for(i=1;i<=n;i++)
                        x[i]=0;
                     //--- change values
                     x[j]=1;
                     xj=1;
                     s=0;
                     xmax=0;
                    }
                 }
              }
            //--- Scale x if necessary to avoid overflow when adding a
            //--- multiple of column j of A.
            if(xj>1.0)
              {
               rec=1/xj;
               //--- check
               if(cnorm[j]>(bignum-xmax)*rec)
                 {
                  //--- Scale x by 1/(2*abs(x(j))).
                  rec=rec*0.5;
                  for(i_=1;i_<=n;i_++)
                     x[i_]=rec*x[i_];
                  s=s*rec;
                 }
              }
            else
              {
               //--- check
               if(xj*cnorm[j]>bignum-xmax)
                 {
                  //--- Scale x by 1/2.
                  for(i_=1;i_<=n;i_++)
                     x[i_]=0.5*x[i_];
                  s=s*0.5;
                 }
              }
            //--- check
            if(upper)
              {
               //--- check
               if(j>1)
                 {
                  //--- Compute the update
                  //--- x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
                  v=x[j]*tscal;
                  jm1=j-1;
                  //--- change x
                  for(i_=1;i_<=jm1;i_++)
                     x[i_]=x[i_]-v*a[i_][j];
                  i=1;
                  for(k=2;k<=j-1;k++)
                    {
                     //--- check
                     if(MathAbs(x[k])>MathAbs(x[i]))
                        i=k;
                    }
                  xmax=MathAbs(x[i]);
                 }
              }
            else
              {
               //--- check
               if(j<n)
                 {
                  //--- Compute the update
                  //--- x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
                  jp1=j+1;
                  v=x[j]*tscal;
                  //--- change x
                  for(i_=jp1;i_<=n;i_++)
                     x[i_]=x[i_]-v*a[i_][j];
                  i=j+1;
                  for(k=j+2;k<=n;k++)
                    {
                     //--- check
                     if(MathAbs(x[k])>MathAbs(x[i]))
                        i=k;
                    }
                  xmax=MathAbs(x[i]);
                 }
              }
            j=j+jinc;
           }
        }
      else
        {
         //--- Solve A' * x = b
         j=jfirst;
         while((jinc>0 && j<=jlast) || (jinc<0 && j>=jlast))
           {
            //--- Compute x(j) = b(j) - sum A(k,j)*x(k).
            //--- k<>j
            xj=MathAbs(x[j]);
            uscal=tscal;
            rec=1/MathMax(xmax,1);
            //--- check
            if(cnorm[j]>(bignum-xj)*rec)
              {
               //--- If x(j) could overflow,scale x by 1/(2*XMAX).
               rec=rec*0.5;
               //--- check
               if(nounit)
                  tjjs=a[j][j]*tscal;
               else
                  tjjs=tscal;
               tjj=MathAbs(tjjs);
               //--- check
               if(tjj>1.0)
                 {
                  //--- Divide by A(j,j) when scaling x if A(j,j) > 1.
                  rec=MathMin(1,rec*tjj);
                  uscal=uscal/tjjs;
                 }
               //--- check
               if(rec<1.0)
                 {
                  for(i_=1;i_<=n;i_++)
                     x[i_]=rec*x[i_];
                  s=s*rec;
                  xmax=xmax*rec;
                 }
              }
            sumj=0;
            //--- check
            if(uscal==1.0)
              {
               //--- If the scaling needed for A in the dot product is 1,
               //--- call DDOT to perform the dot product.
               if(upper)
                 {
                  //--- check
                  if(j>1)
                    {
                     jm1=j-1;
                     sumj=0.0;
                     for(i_=1;i_<=jm1;i_++)
                        sumj+=a[i_][j]*x[i_];
                    }
                  else
                     sumj=0;
                 }
               else
                 {
                  //--- check
                  if(j<n)
                    {
                     jp1=j+1;
                     sumj=0.0;
                     for(i_=jp1;i_<=n;i_++)
                        sumj+=a[i_][j]*x[i_];
                    }
                 }
              }
            else
              {
               //--- Otherwise, use in-line code for the dot product.
               if(upper)
                 {
                  for(i=1;i<=j-1;i++)
                    {
                     v=a[i][j]*uscal;
                     sumj=sumj+v*x[i];
                    }
                 }
               else
                 {
                  //--- check
                  if(j<n)
                    {
                     for(i=j+1;i<=n;i++)
                       {
                        v=a[i][j]*uscal;
                        sumj=sumj+v*x[i];
                       }
                    }
                 }
              }
            if(uscal==tscal)
              {
               //--- Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
               //--- was not used to scale the dotproduct.
               x[j]=x[j]-sumj;
               xj=MathAbs(x[j]);
               flg=0;
               //--- check
               if(nounit)
                  tjjs=a[j][j]*tscal;
               else
                 {
                  tjjs=tscal;
                  //--- check
                  if(tscal==1.0)
                     flg=150;
                 }
               //--- Compute x(j) = x(j) / A(j,j), scaling if necessary.
               if(flg!=150)
                 {
                  tjj=MathAbs(tjjs);
                  //--- check
                  if(tjj>smlnum)
                    {
                     //--- abs(A(j,j)) > SMLNUM:
                     if(tjj<1.0)
                       {
                        //--- check
                        if(xj>(double)(tjj*bignum))
                          {
                           //--- Scale X by 1/abs(x(j)).
                           rec=1/xj;
                           for(i_=1;i_<=n;i_++)
                              x[i_]=rec*x[i_];
                           s=s*rec;
                           xmax=xmax*rec;
                          }
                       }
                     x[j]=x[j]/tjjs;
                    }
                  else
                    {
                     //--- check
                     if(tjj>0.0)
                       {
                        //--- 0 < abs(A(j,j)) <=SMLNUM:
                        if(xj>(double)(tjj*bignum))
                          {
                           //--- Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
                           rec=tjj*bignum/xj;
                           for(i_=1;i_<=n;i_++)
                              x[i_]=rec*x[i_];
                           s=s*rec;
                           xmax=xmax*rec;
                          }
                        x[j]=x[j]/tjjs;
                       }
                     else
                       {
                        //--- A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                        //--- scale = 0, and compute a solution to A'*x = 0.
                        for(i=1;i<=n;i++)
                           x[i]=0;
                        x[j]=1;
                        s=0;
                        xmax=0;
                       }
                    }
                 }
              }
            else
              {
               //--- Compute x(j) := x(j) / A(j,j)  - sumj if the dot
               //--- product has already been divided by 1/A(j,j).
               x[j]=x[j]/tjjs-sumj;
              }
            xmax=MathMax(xmax,MathAbs(x[j]));
            j=j+jinc;
           }
        }
      s=s/tscal;
     }
//--- Scale the column norms by 1/TSCAL for return.
   if(tscal!=1.0)
     {
      v=1/tscal;
      for(i_=1;i_<=n;i_++)
         cnorm[i_]=v*cnorm[i_];
     }
  }
//+------------------------------------------------------------------+
//| Safe solvers                                                     |
//+------------------------------------------------------------------+
class CSafeSolve
  {
private:
   //--- private method
   static bool       CBasicSolveAndUpdate(al_complex &alpha,al_complex &beta,const double lnmax,const double bnorm,const double maxgrowth,double &xnorm,al_complex &x);
public:
   //--- constructor, destructor
                     CSafeSolve(void);
                    ~CSafeSolve(void);
   //--- public methods
   static bool       RMatrixScaledTrSafeSolve(CMatrixDouble &a,const double sa,const int n,double &x[],const bool isupper,const int trans,const bool isunit,const double maxgrowth);
   static bool       CMatrixScaledTrSafeSolve(CMatrixComplex &a,const double sa,const int n,al_complex &x[],const bool isupper,const int trans,const bool isunit,const double maxgrowth);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSafeSolve::CSafeSolve(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSafeSolve::~CSafeSolve(void)
  {

  }
//+------------------------------------------------------------------+
//| Real implementation of CMatrixScaledTRSafeSolve                  |
//+------------------------------------------------------------------+
static bool CSafeSolve::RMatrixScaledTrSafeSolve(CMatrixDouble &a,const double sa,
                                                 const int n,double &x[],
                                                 const bool isupper,const int trans,
                                                 const bool isunit,const double maxgrowth)
  {
//--- create variables
   bool       result;
   double     lnmax=0;
   double     nrmb=0;
   double     nrmx=0;
   double     vr=0;
   al_complex alpha=0;
   al_complex beta=0;
   al_complex cx=0;
   int        i_=0;
   int        i=0;
//--- create array
   double tmp[];
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return(false);
//--- check
   if(!CAp::Assert(trans==0 || trans==1,__FUNCTION__+": incorrect Trans!"))
      return(false);
//--- initialization
   result=true;
   lnmax=MathLog(CMath::m_maxrealnumber);
//--- Quick return if possible
   if(n<=0)
      return(result);
//--- Load norms: right part and X
   nrmb=0;
   for(i=0;i<n;i++)
      nrmb=MathMax(nrmb,MathAbs(x[i]));
   nrmx=0;
//--- Solve
   ArrayResizeAL(tmp,n);
   result=true;
//--- check
   if(isupper && trans==0)
     {
      //--- U*x = b
      for(i=n-1;i>=0;i--)
        {
         //--- Task is reduced to alpha*x[i] = beta
         if(isunit)
            alpha=sa;
         else
            alpha=a[i][i]*sa;
         //--- check
         if(i<n-1)
           {
            for(i_=i+1;i_<n;i_++)
               tmp[i_]=sa*a[i][i_];
            //--- calculation
            vr=0.0;
            for(i_=i+1;i_<n;i_++)
               vr+=tmp[i_]*x[i_];
            beta=x[i]-vr;
           }
         else
            beta=x[i];
         //--- solve alpha*x[i] = beta
         result=CBasicSolveAndUpdate(alpha,beta,lnmax,nrmb,maxgrowth,nrmx,cx);
         //--- check
         if(!result)
            return(result);
         //--- change values
         x[i]=cx.re;
        }
      //--- return result
      return(result);
     }
//--- check
   if(!isupper && trans==0)
     {
      //--- L*x = b
      for(i=0;i<n;i++)
        {
         //--- Task is reduced to alpha*x[i] = beta
         if(isunit)
            alpha=sa;
         else
            alpha=a[i][i]*sa;
         //--- check
         if(i>0)
           {
            for(i_=0;i_<i;i_++)
               tmp[i_]=sa*a[i][i_];
            //--- calculation
            vr=0.0;
            for(i_=0;i_<i;i_++)
               vr+=tmp[i_]*x[i_];
            beta=x[i]-vr;
           }
         else
            beta=x[i];
         //--- solve alpha*x[i] = beta
         result=CBasicSolveAndUpdate(alpha,beta,lnmax,nrmb,maxgrowth,nrmx,cx);
         //--- check
         if(!result)
            return(result);
         //--- change values
         x[i]=cx.re;
        }
      //--- return result
      return(result);
     }
//--- check
   if(isupper && trans==1)
     {
      //--- U^T*x = b
      for(i=0;i<n;i++)
        {
         //--- Task is reduced to alpha*x[i] = beta
         if(isunit)
            alpha=sa;
         else
            alpha=a[i][i]*sa;
         beta=x[i];
         //--- solve alpha*x[i] = beta
         result=CBasicSolveAndUpdate(alpha,beta,lnmax,nrmb,maxgrowth,nrmx,cx);
         //--- check
         if(!result)
            return(result);
         //--- change values
         x[i]=cx.re;
         //--- update the rest of right part
         if(i<n-1)
           {
            vr=cx.re;
            for(i_=i+1;i_<n;i_++)
               tmp[i_]=sa*a[i][i_];
            for(i_=i+1;i_<n;i_++)
               x[i_]=x[i_]-vr*tmp[i_];
           }
        }
      //--- return result
      return(result);
     }
//--- check
   if(!isupper && trans==1)
     {
      //--- L^T*x = b
      for(i=n-1;i>=0;i--)
        {
         //--- Task is reduced to alpha*x[i] = beta
         if(isunit)
            alpha=sa;
         else
            alpha=a[i][i]*sa;
         beta=x[i];
         //--- solve alpha*x[i] = beta
         result=CBasicSolveAndUpdate(alpha,beta,lnmax,nrmb,maxgrowth,nrmx,cx);
         //--- check
         if(!result)
            return(result);
         //--- change values
         x[i]=cx.re;
         //--- update the rest of right part
         if(i>0)
           {
            vr=cx.re;
            for(i_=0;i_<i;i_++)
               tmp[i_]=sa*a[i][i_];
            for(i_=0;i_<i;i_++)
               x[i_]=x[i_]-vr*tmp[i_];
           }
        }
      //--- return result
      return(result);
     }
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Internal subroutine for safe solution of                         |
//|     SA*op(A)=b                                                   |
//| where A is NxN upper/lower triangular/unitriangular matrix, op(A)|
//| is either identity transform, transposition or Hermitian         |
//| transposition, SA is a scaling factor such that max(|SA*A[i,j]|) |
//| is close to 1.0 in magnutude.                                    |
//| This subroutine limits relative growth of solution (in inf-norm) |
//| by MaxGrowth, returning False if growth exceeds MaxGrowth.       |
//| Degenerate or near-degenerate matrices are handled correctly     |
//| (False is returned) as long as MaxGrowth is significantly less   |
//| than MaxRealNumber/norm(b).                                      |
//+------------------------------------------------------------------+
static bool CSafeSolve::CMatrixScaledTrSafeSolve(CMatrixComplex &a,const double sa,
                                                 const int n,al_complex &x[],
                                                 const bool isupper,const int trans,
                                                 const bool isunit,const double maxgrowth)
  {
//--- create variables
   al_complex Csa;
   bool       result;
   double     lnmax=0;
   double     nrmb=0;
   double     nrmx=0;
   al_complex alpha=0;
   al_complex beta=0;
   al_complex vc=0;
   int        i=0;
   int        i_=0;
//--- create array
   al_complex tmp[];
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return(false);
//--- check
   if(!CAp::Assert((trans==0 || trans==1) || trans==2,__FUNCTION__+": incorrect Trans!"))
      return(false);
//--- initialization
   result=true;
   lnmax=MathLog(CMath::m_maxrealnumber);
//--- Quick return if possible
   if(n<=0)
      return(result);
//--- Load norms: right part and X
   nrmb=0;
   for(i=0;i<n;i++)
      nrmb=MathMax(nrmb,CMath::AbsComplex(x[i]));
   nrmx=0;
//--- Solve
   ArrayResizeAL(tmp,n);
   result=true;
//--- check
   if(isupper && trans==0)
     {
      //--- U*x = b
      for(i=n-1;i>=0;i--)
        {
         //--- Task is reduced to alpha*x[i] = beta
         if(isunit)
            alpha=sa;
         else
            alpha=a[i][i]*sa;
         //--- check
         if(i<n-1)
           {
            for(i_=i+1;i_<n;i_++)
              {
               Csa=sa;
               tmp[i_]=Csa*a[i][i_];
              }
            //--- calculation
            vc=0.0;
            for(i_=i+1;i_<n;i_++)
               vc+=tmp[i_]*x[i_];
            beta=x[i]-vc;
           }
         else
            beta=x[i];
         //--- solve alpha*x[i] = beta
         result=CBasicSolveAndUpdate(alpha,beta,lnmax,nrmb,maxgrowth,nrmx,vc);
         //--- check
         if(!result)
            return(result);
         x[i]=vc;
        }
      //--- return result
      return(result);
     }
//--- check
   if(!isupper && trans==0)
     {
      //--- L*x = b
      for(i=0;i<n;i++)
        {
         //--- Task is reduced to alpha*x[i] = beta
         if(isunit)
            alpha=sa;
         else
            alpha=a[i][i]*sa;
         //--- check
         if(i>0)
           {
            for(i_=0;i_<i;i_++)
              {
               Csa=sa;
               tmp[i_]=Csa*a[i][i_];
              }
            //--- calculation
            vc=0.0;
            for(i_=0;i_<i;i_++)
               vc+=tmp[i_]*x[i_];
            beta=x[i]-vc;
           }
         else
            beta=x[i];
         //--- solve alpha*x[i] = beta
         result=CBasicSolveAndUpdate(alpha,beta,lnmax,nrmb,maxgrowth,nrmx,vc);
         //--- check
         if(!result)
            return(result);
         x[i]=vc;
        }
      //--- return result
      return(result);
     }
//--- check
   if(isupper && trans==1)
     {
      //--- U^T*x = b
      for(i=0;i<n;i++)
        {
         //--- Task is reduced to alpha*x[i] = beta
         if(isunit)
            alpha=sa;
         else
            alpha=a[i][i]*sa;
         beta=x[i];
         //--- solve alpha*x[i] = beta
         result=CBasicSolveAndUpdate(alpha,beta,lnmax,nrmb,maxgrowth,nrmx,vc);
         //--- check
         if(!result)
            return(result);
         x[i]=vc;
         //--- update the rest of right part
         if(i<n-1)
           {
            for(i_=i+1;i_<n;i_++)
              {
               Csa=sa;
               tmp[i_]=Csa*a[i][i_];
              }
            for(i_=i+1;i_<n;i_++)
              {
               x[i_]=x[i_]-vc*tmp[i_];
              }
           }
        }
      //--- return result
      return(result);
     }
//--- check
   if(!isupper && trans==1)
     {
      //--- L^T*x = b
      for(i=n-1;i>=0;i--)
        {
         //--- Task is reduced to alpha*x[i] = beta
         if(isunit)
            alpha=sa;
         else
            alpha=a[i][i]*sa;
         beta=x[i];
         //--- solve alpha*x[i] = beta
         result=CBasicSolveAndUpdate(alpha,beta,lnmax,nrmb,maxgrowth,nrmx,vc);
         //--- check
         if(!result)
            return(result);
         x[i]=vc;
         //--- update the rest of right part
         if(i>0)
           {
            for(i_=0;i_<i;i_++)
              {
               Csa=sa;
               tmp[i_]=Csa*a[i][i_];
              }
            for(i_=0;i_<i;i_++)
              {
               x[i_]=x[i_]-vc*tmp[i_];
              }
           }
        }
      //--- return result
      return(result);
     }
//--- check
   if(isupper && trans==2)
     {
      //--- U^H*x=b
      for(i=0;i<n;i++)
        {
         //--- Task is reduced to alpha*x[i] = beta
         if(isunit)
            alpha=sa;
         else
            alpha=CMath::Conj(a[i][i])*sa;
         beta=x[i];
         //--- solve alpha*x[i] = beta
         result=CBasicSolveAndUpdate(alpha,beta,lnmax,nrmb,maxgrowth,nrmx,vc);
         //--- check
         if(!result)
            return(result);
         x[i]=vc;
         //--- update the rest of right part
         if(i<n-1)
           {
            for(i_=i+1;i_<n;i_++)
              {
               Csa=sa;
               tmp[i_]=Csa*CMath::Conj(a[i][i_]);
              }
            for(i_=i+1;i_<n;i_++)
              {
               x[i_]=x[i_]-vc*tmp[i_];
              }
           }
        }
      //--- return result
      return(result);
     }
//--- check
   if(!isupper && trans==2)
     {
      //--- L^T*x = b
      for(i=n-1;i>=0;i--)
        {
         //--- Task is reduced to alpha*x[i] = beta
         if(isunit)
            alpha=sa;
         else
            alpha=CMath::Conj(a[i][i])*sa;
         beta=x[i];
         //--- solve alpha*x[i] = beta
         result=CBasicSolveAndUpdate(alpha,beta,lnmax,nrmb,maxgrowth,nrmx,vc);
         //--- check
         if(!result)
            return(result);
         x[i]=vc;
         //--- update the rest of right part
         if(i>0)
           {
            for(i_=0;i_<i;i_++)
              {
               Csa=sa;
               tmp[i_]=Csa*CMath::Conj(a[i][i_]);
              }
            for(i_=0;i_<i;i_++)
              {
               x[i_]=x[i_]-vc*tmp[i_];
              }
           }
        }
      //--- return result
      return(result);
     }
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| complex basic solver-updater for reduced linear system           |
//|     alpha*x[i] = beta                                            |
//| solves this equation and updates it in overlfow-safe manner      |
//| (keeping track of relative growth of solution).                  |
//| Parameters:                                                      |
//|     Alpha   -   alpha                                            |
//|     Beta    -   beta                                             |
//|     LnMax   -   precomputed Ln(MaxRealNumber)                    |
//|     BNorm   -   inf-norm of b (right part of original system)    |
//|     MaxGrowth-  maximum growth of norm(x) relative to norm(b)    |
//|     XNorm   -   inf-norm of other components of X (which are     |
//|                 already processed) it is updated by              |
//|                 CBasicSolveAndUpdate.                            |
//|     X       -   solution                                         |
//+------------------------------------------------------------------+
static bool CSafeSolve::CBasicSolveAndUpdate(al_complex &alpha,al_complex &beta,
                                             const double lnmax,const double bnorm,
                                             const double maxgrowth,
                                             double &xnorm,al_complex &x)
  {
//--- create variables
   bool   result;
   double v=0;
//--- initialization
   x=0;
   result=false;
//--- check
   if(alpha==0)
      return(result);
//--- check
   if(beta!=0)
     {
      //--- alpha*x[i]=beta
      v=MathLog(CMath::AbsComplex(beta))-MathLog(CMath::AbsComplex(alpha));
      //--- check
      if(v>lnmax)
         return(result);
      x=beta/alpha;
     }
   else
     {
      //--- alpha*x[i]=0
      x=0;
     }
//--- update NrmX, test growth limit
   xnorm=MathMax(xnorm,CMath::AbsComplex(x));
//--- check
   if(xnorm>maxgrowth*bnorm)
      return(result);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Dot-product                                                      |
//+------------------------------------------------------------------+
class CXblas
  {
private:
   //--- private methods
   static void       XSum(double &w[],const double mx,const int n,double &r,double &rerr);
   static double     XFastPow(const double r,const int n);
public:
   //--- constructor, destructor
                     CXblas(void);
                    ~CXblas(void);
   //--- public methods
   static void       XDot(double &a[],double &b[],const int n,double &temp[],double &r,double &rerr);
   static void       XCDot(al_complex &a[],al_complex &b[],const int n,double &temp[],al_complex &r,double &rerr);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CXblas::CXblas(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CXblas::~CXblas(void)
  {

  }
//+------------------------------------------------------------------+
//| More precise dot-product. Absolute error of subroutine result is |
//| about 1 ulp of max(MX,V), where:                                 |
//|     MX = max( |a[i]*b[i]| )                                      |
//|     V  = |(a,b)|                                                 |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1], vector 1                          |
//|     B       -   array[0..N-1], vector 2                          |
//|     N       -   vectors length, N<2^29.                          |
//|     Temp    -   array[0..N-1], pre-allocated temporary storage   |
//| OUTPUT PARAMETERS                                                |
//|     R       -   (A,B)                                            |
//|     RErr    -   estimate of error. This estimate accounts for    |
//|                 both errors during  calculation of (A,B) and     |
//|                errors introduced by rounding of A and B to fit in|
//| double (about 1 ulp).                                            |
//+------------------------------------------------------------------+
static void CXblas::XDot(double &a[],double &b[],const int n,double &temp[],
                         double &r,double &rerr)
  {
//--- create variables
   int    i=0;
   double mx=0;
   double v=0;
//--- initialization
   r=0;
   rerr=0;
//--- special cases:
//--- * N=0
   if(n==0)
     {
      r=0;
      rerr=0;
      //--- exit the function
      return;
     }
   mx=0;
//--- calculations
   for(i=0;i<=n-1;i++)
     {
      v=a[i]*b[i];
      temp[i]=v;
      mx=MathMax(mx,MathAbs(v));
     }
//--- check
   if(mx==0.0)
     {
      r=0;
      rerr=0;
      //--- exit the function
      return;
     }
//--- function call
   XSum(temp,mx,n,r,rerr);
  }
//+------------------------------------------------------------------+
//| More precise complex dot-product. Absolute error of subroutine   |
//| result is about 1 ulp of max(MX,V), where:                       |
//|     MX = max( |a[i]*b[i]| )                                      |
//|     V  = |(a,b)|                                                 |
//| INPUT PARAMETERS                                                 |
//|     A       -   array[0..N-1], vector 1                          |
//|     B       -   array[0..N-1], vector 2                          |
//|     N       -   vectors length, N<2^29.                          |
//|     Temp    -   array[0..2*N-1], pre-allocated temporary storage |
//| OUTPUT PARAMETERS                                                |
//|     R       -   (A,B)                                            |
//|     RErr    -   estimate of error. This estimate accounts for    |
//|                 both errors during  calculation of (A,B) and     |
//|                 errors introduced by rounding of A and B to fit  |
//|                 in double (about 1 ulp).                         |
//+------------------------------------------------------------------+
static void CXblas::XCDot(al_complex &a[],al_complex &b[],const int n,double &temp[],
                          al_complex &r,double &rerr)
  {
//--- create variables
   int    i=0;
   double mx=0;
   double v=0;
   double rerrx=0;
   double rerry=0;
//--- initialization
   r=0;
   rerr=0;
//--- special cases:
//--- * N=0
   if(n==0)
     {
      r=0;
      rerr=0;
      //--- exit the function
      return;
     }
//--- calculate real part
   mx=0;
   for(i=0;i<=n-1;i++)
     {
      //--- change values
      v=a[i].re*b[i].re;
      temp[2*i+0]=v;
      mx=MathMax(mx,MathAbs(v));
      v=-(a[i].im*b[i].im);
      temp[2*i+1]=v;
      mx=MathMax(mx,MathAbs(v));
     }
//--- check
   if(mx==0.0)
     {
      r.re=0;
      rerrx=0;
     }
   else
      XSum(temp,mx,2*n,r.re,rerrx);
//--- calculate imaginary part
   mx=0;
   for(i=0;i<=n-1;i++)
     {
      //--- change values
      v=a[i].re*b[i].im;
      temp[2*i+0]=v;
      mx=MathMax(mx,MathAbs(v));
      v=a[i].im*b[i].re;
      temp[2*i+1]=v;
      mx=MathMax(mx,MathAbs(v));
     }
//--- check
   if(mx==0.0)
     {
      r.im=0;
      rerry=0;
     }
   else
      XSum(temp,mx,2*n,r.im,rerry);
//--- total error
   if(rerrx==0.0 && rerry==0.0)
      rerr=0;
   else
      rerr=MathMax(rerrx,rerry)*MathSqrt(1+CMath::Sqr(MathMin(rerrx,rerry)/MathMax(rerrx,rerry)));
  }
//+------------------------------------------------------------------+
//| Internal subroutine for extra-precise calculation of SUM(w[i]).  |
//| INPUT PARAMETERS:                                                |
//|     W   -   array[0..N-1], values to be added                    |
//|             W is modified during calculations.                   |
//|     MX  -   max(W[i])                                            |
//|     N   -   array size                                           |
//| OUTPUT PARAMETERS:                                               |
//|     R   -   SUM(w[i])                                            |
//|     RErr-   error estimate for R                                 |
//+------------------------------------------------------------------+
static void CXblas::XSum(double &w[],const double mx,const int n,double &r,
                         double &rerr)
  {
//--- create variables
   int    i=0;
   int    k=0;
   int    ks=0;
   double v=0;
   double s=0;
   double ln2=0;
   double chunk=0;
   double invchunk=0;
   bool   allzeros;
   int    i_=0;
//--- initialization
   r=0;
   rerr=0;
//--- special cases:
//--- * N=0
//--- * N is too large to use integer arithmetics
   if(n==0)
     {
      r=0;
      rerr=0;
      //--- exit the function
      return;
     }
//--- check
   if(mx==0.0)
     {
      r=0;
      rerr=0;
      //--- exit the function
      return;
     }
//--- check
   if(!CAp::Assert(n<536870912,__FUNCTION__+": N is too large!"))
      return;
//--- Prepare
   ln2=MathLog(2);
   rerr=mx*CMath::m_machineepsilon;
//--- 1. find S such that 0.5<=S*MX<1
//--- 2. multiply W by S, so task is normalized in some sense
//--- 3. S:=1/S so we can obtain original vector multiplying by S
   k=(int)MathRound(MathLog(mx)/ln2);
   s=XFastPow(2,-k);
//--- change s
   while(s*mx>=1.0)
      s=0.5*s;
   while(s*mx<0.5)
      s=2*s;
   for(i_=0;i_<=n-1;i_++)
      w[i_]=s*w[i_];
   s=1/s;
//--- find Chunk=2^M such that N*Chunk<2^29
//--- we have chosen upper limit (2^29) with enough space left
//--- to tolerate possible problems with rounding and N's close
//--- to the limit, so we don't want to be very strict here.
   k=(int)(MathLog((double)536870912/(double)n)/ln2);
   chunk=XFastPow(2,k);
//--- check
   if(chunk<2.0)
      chunk=2;
   invchunk=1/chunk;
//--- calculate result
   r=0;
   for(i_=0;i_<=n-1;i_++)
      w[i_]=chunk*w[i_];
//--- cycle
   while(true)
     {
      //--- change values
      s=s*invchunk;
      allzeros=true;
      ks=0;
      for(i=0;i<=n-1;i++)
        {
         v=w[i];
         k=(int)(v);
         //--- check
         if(v!=k)
            allzeros=false;
         w[i]=chunk*(v-k);
         ks=ks+k;
        }
      r=r+s*ks;
      v=MathAbs(r);
      //--- check
      if(allzeros || s*n+mx==mx)
         break;
     }
//--- correct error
   rerr=MathMax(rerr,MathAbs(r)*CMath::m_machineepsilon);
  }
//+------------------------------------------------------------------+
//| Fast Pow                                                         |
//+------------------------------------------------------------------+
static double CXblas::XFastPow(const double r,const int n)
  {
//--- create a variable
   double result=0;
//--- check
   if(n>0)
     {
      //--- check
      if(n%2==0)
         result=CMath::Sqr(XFastPow(r,n/2));
      else
         result=r*XFastPow(r,n-1);
      //--- return result
      return(result);
     }
//--- check
   if(n==0)
      result=1;
//--- check
   if(n<0)
      result=XFastPow(1/r,-n);
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CLinMin                                      |
//+------------------------------------------------------------------+
class CLinMinState
  {
public:
   //--- class variables
   bool              m_brackt;
   bool              m_stage1;
   int               m_infoc;
   double            m_dg;
   double            m_dgm;
   double            m_dginit;
   double            m_dgtest;
   double            m_dgx;
   double            m_dgxm;
   double            m_dgy;
   double            m_dgym;
   double            m_finit;
   double            m_ftest1;
   double            m_fm;
   double            m_fx;
   double            m_fxm;
   double            m_fy;
   double            m_fym;
   double            m_stx;
   double            m_sty;
   double            m_stmin;
   double            m_stmax;
   double            m_width;
   double            m_width1;
   double            m_xtrapf;
   //--- constructor, destructor
                     CLinMinState(void);
                    ~CLinMinState(void);
   //--- create a copy
   void              Copy(CLinMinState &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLinMinState::CLinMinState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLinMinState::~CLinMinState(void)
  {

  }
//+------------------------------------------------------------------+
//| Create a copy                                                    |
//+------------------------------------------------------------------+
void CLinMinState::Copy(CLinMinState &obj)
  {
//--- copy variables
   m_brackt=obj.m_brackt;
   m_stage1=obj.m_stage1;
   m_infoc=obj.m_infoc;
   m_dg=obj.m_dg;
   m_dgm=obj.m_dgm;
   m_dginit=obj.m_dginit;
   m_dgtest=obj.m_dgtest;
   m_dgx=obj.m_dgx;
   m_dgxm=obj.m_dgxm;
   m_dgy=obj.m_dgy;
   m_dgym=obj.m_dgym;
   m_finit=obj.m_finit;
   m_ftest1=obj.m_ftest1;
   m_fm=obj.m_fm;
   m_fx=obj.m_fx;
   m_fxm=obj.m_fxm;
   m_fy=obj.m_fy;
   m_fym=obj.m_fym;
   m_stx=obj.m_stx;
   m_sty=obj.m_sty;
   m_stmin=obj.m_stmin;
   m_stmax=obj.m_stmax;
   m_width=obj.m_width;
   m_width1=obj.m_width1;
   m_xtrapf=obj.m_xtrapf;
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CLinMin                                      |
//+------------------------------------------------------------------+
class CArmijoState
  {
public:
   bool              m_needf;
   double            m_x[];
   double            m_f;
   int               m_n;
   double            m_xbase[];
   double            m_s[];
   double            m_stplen;
   double            m_fcur;
   double            m_stpmax;
   int               m_fmax;
   int               m_nfev;
   int               m_info;
   RCommState        m_rstate;
   //--- constructor, destructor
                     CArmijoState(void);
                    ~CArmijoState(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CArmijoState::CArmijoState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CArmijoState::~CArmijoState(void)
  {

  }
//+------------------------------------------------------------------+
//| Minimization of linear forms                                     |
//+------------------------------------------------------------------+
class CLinMin
  {
private:
   //--- private methods
   static void       MCStep(double &stx,double &fx,double &dx,double &sty,double &fy,double &dy,double &stp,double fp,double dp,bool &m_brackt,double stmin,double stmax,int &info);
   //--- auxiliary functions for ArmijoIteration
   static void       Func_lbl_rcomm(CArmijoState &state,int n,double v);
   static bool       Func_lbl_6(CArmijoState &state,int &n,double &v);
   static bool       Func_lbl_10(CArmijoState &state,int &n,double &v);
public:
   //--- class constants
   static const double m_ftol;
   static const double m_xtol;
   static const int  m_maxfev;
   static const double m_stpmin;
   static const double m_defstpmax;
   static const double m_armijofactor;
   //--- constructor, destructor
                     CLinMin(void);
                    ~CLinMin(void);
   //--- public methods
   static void       LinMinNormalized(double &d[],double &stp,const int n);
   static void       MCSrch(const int n,double &x[],double &f,double &g[],double &s[],double &stp,double stpmax,double gtol,int &info,int &nfev,double &wa[],CLinMinState &state,int &stage);
   static void       ArmijoCreate(const int n,double &x[],const double f,double &s[],const double stp,const double stpmax,const int ffmax,CArmijoState &state);
   static void       ArmijoResults(CArmijoState &state,int &info,double &stp,double &f);
   static bool       ArmijoIteration(CArmijoState &state);
  };
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const double CLinMin::m_ftol=0.001;
const double CLinMin::m_xtol=100*CMath::m_machineepsilon;
const int    CLinMin::m_maxfev=20;
const double CLinMin::m_stpmin=1.0E-50;
const double CLinMin::m_defstpmax=1.0E+50;
const double CLinMin::m_armijofactor=1.3;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLinMin::CLinMin(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLinMin::~CLinMin(void)
  {

  }
//+------------------------------------------------------------------+
//| Normalizes direction/step pair: makes |D|=1,scales Stp.          |
//| If |D|=0,it returns,leavind D/Stp unchanged.                     |
//+------------------------------------------------------------------+
static void CLinMin::LinMinNormalized(double &d[],double &stp,const int n)
  {
//--- create variables
   double mx=0;
   double s=0;
   int    i=0;
   int    i_=0;
//--- first, scale D to avoid underflow/overflow durng squaring
   mx=0;
   for(i=0;i<=n-1;i++)
      mx=MathMax(mx,MathAbs(d[i]));
//--- check
   if(mx==0.0)
      return;
   s=1/mx;
   for(i_=0;i_<=n-1;i_++)
      d[i_]=s*d[i_];
   stp=stp/s;
//--- normalize D
   s=0.0;
   for(i_=0;i_<=n-1;i_++)
      s+=d[i_]*d[i_];
   s=1/MathSqrt(s);
   for(i_=0;i_<=n-1;i_++)
      d[i_]=s*d[i_];
   stp=stp/s;
  }
//+------------------------------------------------------------------+
//| The purpose of MCSrch is to find a step which satisfies a        |
//| sufficient decrease condition and a curvature condition.         |
//| At each stage the subroutine updates an interval of uncertainty  |
//| with endpoints stx and sty. The interval of uncertainty is       |
//| initially chosen so that it contains a minimizer of the modified |
//| function                                                         |
//|     f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).                   |
//| If a step is obtained for which the modified function has a      |
//| nonpositive function value and nonnegative derivative, then the  |
//| interval of uncertainty is chosen so that it contains a minimizer| 
//| of f(x+stp*s).                                                   |
//| The  algorithm  is  designed to find a step which satisfies the  |
//| sufficient decrease condition                                    |
//|     f(x+stp*s) .le. f(x) + ftol*stp*(gradf(x)'s),                |
//| and the curvature condition                                      |
//|     abs(gradf(x+stp*s)'s)) .le. gtol*abs(gradf(x)'s).            |
//| If  ftol is less than gtol and if, for example, the function is  |
//| bounded below, then there is always a step which satisfies both  |
//| conditions. If no step can be found which satisfies both         |
//| conditions, then the algorithm usually stops when rounding errors| 
//| prevent further progress. In this case stp only satisfies the    |
//| sufficient decrease condition.                                   |
//| :::::::::::::important notes:::::::::::::                        |
//| note 1:                                                          |
//| This routine  guarantees that it will stop at the last point     |
//| where function value was calculated. It won't make several       |
//| additional function evaluations after finding good point. So if  |
//| you store function evaluations requested by this routine, you can|
//| be sure that last one is the point where we've stopped.          |
//| NOTE 2:                                                          |
//| when 0<StpMax<StpMin, algorithm will terminate with INFO=5 and   |
//| Stp=0.0                                                          |
//| :::::::::::::::::::::::::::::::::::::::::                        |
//| Parameters descriprion                                           |
//| Stage is zero on first call, zero on final exit                  |
//| N is a positive integer input variable set to the number of      |
//| variables.                                                       |
//| X is  an  array  of  length n. on input it must contain the base |
//| point for the line search. on output it contains x+stp*s.        |
//| F is  a  variable. on input it must contain the value of f at x. | 
//| on output it contains the value of f at x + stp*s.               |
//| G is an array of length n. on input it must contain the gradient |
//| of f at x. on output it contains the gradient of f at x + stp*s. |
//| S is an input array of length n which specifies the search       |
//| direction.                                                       |
//| Stp  is  a nonnegative variable. on input stp contains an initial| 
//| estimate of a satisfactory step. on output stp contains the final|
//| estimate.                                                        |
//| Ftol and gtol are nonnegative input variables. termination occurs|
//| when the sufficient decrease condition and the directional       |
//| derivative condition are satisfied.                              |
//| Xtol is a nonnegative input variable. termination occurs when the|
//| relative width of the interval of uncertainty is at most xtol.   |
//| Stpmin and stpmax are nonnegative input variables which specify  |
//| lower and upper bounds for the step.                             |
//| Maxfev is a positive integer input variable. termination occurs  |
//| when the number of calls to fcn is at least maxfev by the end of |
//| an iteration.                                                    |
//| Info is an integer output variable set as follows:               |
//|     info = 0  improper input parameters.                         |
//|     info = 1  the sufficient decrease condition and the          |
//|               directional derivative condition hold.             |
//|     info = 2  relative width of the interval of uncertainty      |
//|              is at most xtol.                                    |
//|     info = 3  number of calls to fcn has reached maxfev.         |
//|     info = 4  the step is at the lower bound stpmin.             |
//|     info = 5  the step is at the upper bound stpmax.             |
//|     info = 6  rounding errors prevent further progress.          |
//|               there may not be a step which satisfies the        |
//|               sufficient decrease and curvature conditions.      |
//|               tolerances may be too small.                       |
//| Nfev is an integer output variable set to the number of calls to |
//| fcn.                                                             |
//| wa is a work array of length n.                                  |
//| argonne national laboratory. minpack project. june 1983          |
//| Jorge J. More', David J. Thuente                                 |
//+------------------------------------------------------------------+
static void CLinMin::MCSrch(const int n,double &x[],double &f,double &g[],
                            double &s[],double &stp,double stpmax,double gtol,
                            int &info,int &nfev,double &wa[],
                            CLinMinState &state,int &stage)
  {
//--- create variables
   double v=0;
   double p5=0;
   double p66=0;
   double zero=0;
   int    i_=0;
//--- init
   p5=0.5;
   p66=0.66;
   state.m_xtrapf=4.0;
   zero=0;
//--- check
   if(stpmax==0.0)
      stpmax=m_defstpmax;
//--- check
   if(stp<m_stpmin)
      stp=m_stpmin;
//--- check
   if(stp>stpmax)
      stp=stpmax;
//--- Main cycle
   while(true)
     {
      //--- check
      if(stage==0)
        {
         //--- NEXT
         stage=2;
         continue;
        }
      //--- check
      if(stage==2)
        {
         state.m_infoc=1;
         info=0;
         //--- check the input parameters for errors.
         if(stpmax<m_stpmin && stpmax>0.0)
           {
            info=5;
            stp=0.0;
            //--- exit the function
            return;
           }
         //--- check
         if(n<=0 || stp<=0.0 || m_ftol<0.0 || gtol<zero || m_xtol<zero || m_stpmin<zero || stpmax<m_stpmin || m_maxfev<=0)
           {
            stage=0;
            return;
           }
         //--- compute the initial gradient in the search direction
         //--- and check that s is a descent direction.
         v=0.0;
         for(i_=0;i_<=n-1;i_++)
            v+=g[i_]*s[i_];
         state.m_dginit=v;
         //--- check
         if(state.m_dginit>=0.0)
           {
            stage=0;
            return;
           }
         //--- initialize local variables.
         state.m_brackt=false;
         state.m_stage1=true;
         nfev=0;
         state.m_finit=f;
         state.m_dgtest=m_ftol*state.m_dginit;
         state.m_width=stpmax-m_stpmin;
         state.m_width1=state.m_width/p5;
         for(i_=0;i_<=n-1;i_++)
            wa[i_]=x[i_];
         //--- the variables stx,fx,dgx contain the values of the step,
         //--- function,and directional derivative at the best step.
         //--- the variables sty,fy,dgy contain the value of the step,
         //--- function,and derivative at the other endpoint of
         //--- the interval of uncertainty.
         //--- the variables stp,f,dg contain the values of the step,
         //--- function,and derivative at the current step.
         state.m_stx=0;
         state.m_fx=state.m_finit;
         state.m_dgx=state.m_dginit;
         state.m_sty=0;
         state.m_fy=state.m_finit;
         state.m_dgy=state.m_dginit;
         //--- next
         stage=3;
         continue;
        }
      //--- check
      if(stage==3)
        {
         //--- start of iteration.
         //--- set the minimum and maximum steps to correspond
         //--- to the present interval of uncertainty.
         if(state.m_brackt)
           {
            //--- check
            if(state.m_stx<state.m_sty)
              {
               state.m_stmin=state.m_stx;
               state.m_stmax=state.m_sty;
              }
            else
              {
               state.m_stmin=state.m_sty;
               state.m_stmax=state.m_stx;
              }
           }
         else
           {
            state.m_stmin=state.m_stx;
            state.m_stmax=stp+state.m_xtrapf*(stp-state.m_stx);
           }
         //--- force the step to be within the bounds stpmax and stpmin.
         if(stp>stpmax)
            stp=stpmax;
         //--- check
         if(stp<m_stpmin)
            stp=m_stpmin;
         //--- if an unusual termination is to occur then let
         //--- stp be the lowest point obtained so far.
         if((state.m_brackt && (stp<=state.m_stmin || stp>=state.m_stmax)) || nfev>=m_maxfev-1 || state.m_infoc==0 || 
            (state.m_brackt && state.m_stmax-state.m_stmin<=m_xtol*state.m_stmax))
           {
            stp=state.m_stx;
           }
         //--- evaluate the function and gradient at stp
         //--- and compute the directional derivative.
         for(i_=0;i_<=n-1;i_++)
            x[i_]=wa[i_];
         for(i_=0;i_<=n-1;i_++)
            x[i_]=x[i_]+stp*s[i_];
         //--- NEXT
         stage=4;
         return;
        }
      //--- check
      if(stage==4)
        {
         info=0;
         nfev=nfev+1;
         v=0.0;
         for(i_=0;i_<=n-1;i_++)
            v+=g[i_]*s[i_];
         state.m_dg=v;
         state.m_ftest1=state.m_finit+stp*state.m_dgtest;
         //--- test for convergence.
         if((state.m_brackt && (stp<=state.m_stmin || stp>=state.m_stmax)) || state.m_infoc==0)
            info=6;
         //--- check
         if((stp==stpmax && f<=state.m_ftest1) && state.m_dg<=state.m_dgtest)
            info=5;
         //--- check
         if(stp==m_stpmin && (f>state.m_ftest1 || state.m_dg>=state.m_dgtest))
            info=4;
         //--- check
         if(nfev>=m_maxfev)
            info=3;
         //--- check
         if(state.m_brackt && state.m_stmax-state.m_stmin<=m_xtol*state.m_stmax)
            info=2;
         //--- check
         if(f<=state.m_ftest1 && MathAbs(state.m_dg)<=-(gtol*state.m_dginit))
            info=1;
         //--- check for termination.
         if(info!=0)
           {
            stage=0;
            return;
           }
         //--- in the first stage we seek a step for which the modified
         //--- function has a nonpositive value and nonnegative derivative.
         if((state.m_stage1 && f<=state.m_ftest1) && state.m_dg>=MathMin(m_ftol,gtol)*state.m_dginit)
            state.m_stage1=false;
         //--- a modified function is used to predict the step only if
         //--- we have not obtained a step for which the modified
         //--- function has a nonpositive function value and nonnegative
         //--- derivative,and if a lower function value has been
         //--- obtained but the decrease is not sufficient.
         if((state.m_stage1 && f<=state.m_fx) && f>state.m_ftest1)
           {
            //--- define the modified function and derivative values.
            state.m_fm=f-stp*state.m_dgtest;
            state.m_fxm=state.m_fx-state.m_stx*state.m_dgtest;
            state.m_fym=state.m_fy-state.m_sty*state.m_dgtest;
            state.m_dgm=state.m_dg-state.m_dgtest;
            state.m_dgxm=state.m_dgx-state.m_dgtest;
            state.m_dgym=state.m_dgy-state.m_dgtest;
            //--- call cstep to update the interval of uncertainty
            //--- and to compute the new step.
            MCStep(state.m_stx,state.m_fxm,state.m_dgxm,state.m_sty,state.m_fym,state.m_dgym,stp,state.m_fm,state.m_dgm,state.m_brackt,state.m_stmin,state.m_stmax,state.m_infoc);
            //--- reset the function and gradient values for f.
            state.m_fx=state.m_fxm+state.m_stx*state.m_dgtest;
            state.m_fy=state.m_fym+state.m_sty*state.m_dgtest;
            state.m_dgx=state.m_dgxm+state.m_dgtest;
            state.m_dgy=state.m_dgym+state.m_dgtest;
           }
         else
           {
            //--- call mcstep to update the interval of uncertainty
            //--- and to compute the new step.
            MCStep(state.m_stx,state.m_fx,state.m_dgx,state.m_sty,state.m_fy,state.m_dgy,stp,f,state.m_dg,state.m_brackt,state.m_stmin,state.m_stmax,state.m_infoc);
           }
         //--- force a sufficient decrease in the size of the
         //--- interval of uncertainty.
         if(state.m_brackt)
           {
            //--- check
            if(MathAbs(state.m_sty-state.m_stx)>=p66*state.m_width1)
               stp=state.m_stx+p5*(state.m_sty-state.m_stx);
            state.m_width1=state.m_width;
            state.m_width=MathAbs(state.m_sty-state.m_stx);
           }
         //--- next.
         stage=3;
         continue;
        }
     }
  }
//+------------------------------------------------------------------+
//| These functions perform Armijo line search using at most FMAX    |
//| function evaluations. It doesn't enforce some kind of            |
//| "sufficient decrease" criterion - it just tries different Armijo |
//| steps and returns optimum found so far.                          |
//| Optimization is done using F-rcomm interface:                    |
//| * ArmijoCreate initializes State structure                       |
//|   (reusing previously allocated buffers)                         |
//| * ArmijoIteration is subsequently called                         |
//| * ArmijoResults returns results                                  |
//| INPUT PARAMETERS:                                                |
//|     N       -   problem size                                     |
//|     X       -   array[N], starting point                         |
//|     F       -   F(X+S*STP)                                       |
//|     S       -   step direction, S>0                              |
//|     STP     -   step length                                      |
//|     STPMAX  -   maximum value for STP or zero (if no limit is    |
//|                 imposed)                                         |
//|     FMAX    -   maximum number of function evaluations           |
//|     State   -   optimization state                               |
//+------------------------------------------------------------------+
static void CLinMin::ArmijoCreate(const int n,double &x[],const double f,
                                  double &s[],const double stp,const double stpmax,
                                  const int ffmax,CArmijoState &state)
  {
//--- create a variable
   int i_=0;
//--- check
   if(CAp::Len(state.m_x)<n)
      ArrayResizeAL(state.m_x,n);
//--- check
   if(CAp::Len(state.m_xbase)<n)
      ArrayResizeAL(state.m_xbase,n);
//--- check
   if(CAp::Len(state.m_s)<n)
      ArrayResizeAL(state.m_s,n);
//--- copy
   state.m_stpmax=stpmax;
   state.m_fmax=ffmax;
   state.m_stplen=stp;
   state.m_fcur=f;
   state.m_n=n;
   for(i_=0;i_<=n-1;i_++)
      state.m_xbase[i_]=x[i_];
   for(i_=0;i_<=n-1;i_++)
      state.m_s[i_]=s[i_];
//--- allocation
   ArrayResizeAL(state.m_rstate.ia,1);
   ArrayResizeAL(state.m_rstate.ra,1);
   state.m_rstate.stage=-1;
  }
//+------------------------------------------------------------------+
//| Results of Armijo search                                         |
//| OUTPUT PARAMETERS:                                               |
//|     INFO    -   on output it is set to one of the return codes:  |
//|                 * 0     improper input params                    |
//|                 * 1     optimum step is found with at most FMAX  |
//|                         evaluations                              |
//|                 * 3     FMAX evaluations were used,              |
//|                         X contains optimum found so far          |
//|                 * 4     step is at lower bound STPMIN            |
//|                 * 5     step is at upper bound                   |
//|     STP     -   step length (in case of failure it is still      |
//|                 returned)                                        |
//|     F       -   function value (in case of failure it is still   |
//|                 returned)                                        |
//+------------------------------------------------------------------+
static void CLinMin::ArmijoResults(CArmijoState &state,int &info,
                                   double &stp,double &f)
  {
//--- change values
   info=state.m_info;
   stp=state.m_stplen;
   f=state.m_fcur;
  }
//+------------------------------------------------------------------+
//| Internal subroutine for MCSrch                                   |
//+------------------------------------------------------------------+
static void CLinMin::MCStep(double &stx,double &fx,double &dx,double &sty,
                            double &fy,double &dy,double &stp,double fp,
                            double dp,bool &m_brackt,double stmin,
                            double stmax,int &info)
  {
//--- create variables
   bool   bound;
   double gamma=0;
   double p=0;
   double q=0;
   double r=0;
   double s=0;
   double sgnd=0;
   double stpc=0;
   double stpf=0;
   double stpq=0;
   double theta=0;
//--- initialization
   info=0;
//--- check the input parameters for errors.
   if(((m_brackt && (stp<=MathMin(stx,sty) || stp>=MathMax(stx,sty))) || dx*(stp-stx)>=0.0) || stmax<stmin)
      return;
//--- determine if the derivatives have opposite sign.
   sgnd=dp*(dx/MathAbs(dx));
//--- first case. a higher function value.
//--- the minimum is bracketed. if the cubic step is closer
//--- to stx than the quadratic step,the cubic step is taken,
//--- else the average of the cubic and quadratic steps is taken.
   if(fp>fx)
     {
      //--- initialization
      info=1;
      bound=true;
      theta=3*(fx-fp)/(stp-stx)+dx+dp;
      s=MathMax(MathAbs(theta),MathMax(MathAbs(dx),MathAbs(dp)));
      gamma=s*MathSqrt(CMath::Sqr(theta/s)-dx/s*(dp/s));
      //--- check
      if(stp<stx)
         gamma=-gamma;
      //--- initialization
      p=gamma-dx+theta;
      q=gamma-dx+gamma+dp;
      r=p/q;
      stpc=stx+r*(stp-stx);
      stpq=stx+dx/((fx-fp)/(stp-stx)+dx)/2*(stp-stx);
      //--- check
      if(MathAbs(stpc-stx)<MathAbs(stpq-stx))
         stpf=stpc;
      else
         stpf=stpc+(stpq-stpc)/2;
      m_brackt=true;
     }
   else
     {
      //--- check
      if(sgnd<0.0)
        {
         //--- second case. a lower function value and derivatives of
         //--- opposite sign. the minimum is bracketed. if the cubic
         //--- step is closer to stx than the quadratic (secant) step,
         //--- the cubic step is taken, else the quadratic step is taken.
         info=2;
         bound=false;
         theta=3*(fx-fp)/(stp-stx)+dx+dp;
         s=MathMax(MathAbs(theta),MathMax(MathAbs(dx),MathAbs(dp)));
         gamma=s*MathSqrt(CMath::Sqr(theta/s)-dx/s*(dp/s));
         //--- check
         if(stp>stx)
            gamma=-gamma;
         //--- initialization
         p=gamma-dp+theta;
         q=gamma-dp+gamma+dx;
         r=p/q;
         stpc=stp+r*(stx-stp);
         stpq=stp+dp/(dp-dx)*(stx-stp);
         //--- check
         if(MathAbs(stpc-stp)>MathAbs(stpq-stp))
            stpf=stpc;
         else
            stpf=stpq;
         m_brackt=true;
        }
      else
        {
         //--- check
         if(MathAbs(dp)<MathAbs(dx))
           {
            //--- third case. a lower function value,derivatives of the
            //--- same sign, and the magnitude of the derivative decreases.
            //--- the cubic step is only used if the cubic tends to infinity
            //--- in the direction of the step or if the minimum of the cubic
            //--- is beyond stp. otherwise the cubic step is defined to be
            //--- either stpmin or stpmax. the quadratic (secant) step is also
            //--- computed and if the minimum is bracketed then the the step
            //--- closest to stx is taken, else the step farthest away is taken.
            info=3;
            bound=true;
            theta=3*(fx-fp)/(stp-stx)+dx+dp;
            s=MathMax(MathAbs(theta),MathMax(MathAbs(dx),MathAbs(dp)));
            //--- the case gamma=0 only arises if the cubic does not tend
            //--- to infinity in the direction of the step.
            gamma=s*MathSqrt(MathMax(0,CMath::Sqr(theta/s)-dx/s*(dp/s)));
            //--- check
            if(stp>stx)
               gamma=-gamma;
            p=gamma-dp+theta;
            q=gamma+(dx-dp)+gamma;
            r=p/q;
            //--- check
            if(r<0.0 && (double)(gamma)!=0.0)
               stpc=stp+r*(stx-stp);
            else
              {
               //--- check
               if(stp>stx)
                  stpc=stmax;
               else
                  stpc=stmin;
              }
            stpq=stp+dp/(dp-dx)*(stx-stp);
            //--- check
            if(m_brackt)
              {
               //--- check
               if(MathAbs(stp-stpc)<MathAbs(stp-stpq))
                  stpf=stpc;
               else
                  stpf=stpq;
              }
            else
              {
               //--- check
               if(MathAbs(stp-stpc)>MathAbs(stp-stpq))
                  stpf=stpc;
               else
                  stpf=stpq;
              }
           }
         else
           {
            //--- fourth case. a lower function value,derivatives of the
            //--- same sign, and the magnitude of the derivative does
            //--- not decrease. if the minimum is not bracketed, the step
            //--- is either stpmin or stpmax, else the cubic step is taken.
            info=4;
            bound=false;
            //--- check
            if(m_brackt)
              {
               theta=3*(fp-fy)/(sty-stp)+dy+dp;
               s=MathMax(MathAbs(theta),MathMax(MathAbs(dy),MathAbs(dp)));
               gamma=s*MathSqrt(CMath::Sqr(theta/s)-dy/s*(dp/s));
               //--- check
               if(stp>sty)
                  gamma=-gamma;
               //--- initialization
               p=gamma-dp+theta;
               q=gamma-dp+gamma+dy;
               r=p/q;
               stpc=stp+r*(sty-stp);
               stpf=stpc;
              }
            else
              {
               //--- check
               if(stp>stx)
                  stpf=stmax;
               else
                  stpf=stmin;
              }
           }
        }
     }
//--- update the interval of uncertainty. this update does not
//--- depend on the new step or the case analysis above.
   if(fp>fx)
     {
      //--- set value
      sty=stp;
      fy=fp;
      dy=dp;
     }
   else
     {
      //--- check
      if(sgnd<0.0)
        {
         //--- set value
         sty=stx;
         fy=fx;
         dy=dx;
        }
      //--- set value
      stx=stp;
      fx=fp;
      dx=dp;
     }
//--- compute the new step and safeguard it.
   stpf=MathMin(stmax,stpf);
   stpf=MathMax(stmin,stpf);
   stp=stpf;
//--- check
   if(m_brackt && bound)
     {
      //--- check
      if(sty>stx)
         stp=MathMin(stx+0.66*(sty-stx),stp);
      else
         stp=MathMax(stx+0.66*(sty-stx),stp);
     }
  }
//+------------------------------------------------------------------+
//| This is rcomm-based search function                              |
//+------------------------------------------------------------------+
static bool CLinMin::ArmijoIteration(CArmijoState &state)
  {
//--- create variables
   double v=0;
   int    n=0;
   int    i_=0;
//--- This code initializes locals by:
//--- * random values determined during code
//---   generation - on first subroutine call
//--- * values from previous call - on subsequent calls
   if(state.m_rstate.stage>=0)
     {
      //--- initialization
      n=state.m_rstate.ia[0];
      v=state.m_rstate.ra[0];
     }
   else
     {
      //--- initialization
      n=-983;
      v=-989;
     }
//--- check
   if(state.m_rstate.stage==0)
     {
      state.m_nfev=state.m_nfev+1;
      //--- check
      if(state.m_f>=state.m_fcur)
        {
         //--- Decrease length
         v=state.m_stplen/m_armijofactor;
         //--- copy
         for(i_=0;i_<=n-1;i_++)
            state.m_x[i_]=state.m_xbase[i_];
         for(i_=0;i_<=n-1;i_++)
            state.m_x[i_]=state.m_x[i_]+v*state.m_s[i_];
         state.m_rstate.stage=2;
         //--- Saving state
         Func_lbl_rcomm(state,n,v);
         //--- return result
         return(true);
        }
      //--- change values
      state.m_stplen=v;
      state.m_fcur=state.m_f;
      //--- function call, return result
      return(Func_lbl_6(state,n,v));
     }
//--- check
   if(state.m_rstate.stage==1)
     {
      state.m_nfev=state.m_nfev+1;
      //--- make decision
      if(state.m_f<state.m_fcur)
        {
         //--- change values
         state.m_stplen=v;
         state.m_fcur=state.m_f;
        }
      else
        {
         state.m_info=1;
         //--- return result
         return(false);
        }
      //--- function call, return result
      return(Func_lbl_6(state,n,v));
     }
//--- check
   if(state.m_rstate.stage==2)
     {
      state.m_nfev=state.m_nfev+1;
      //--- check
      if(state.m_f>=state.m_fcur)
        {
         //--- Nothing to be done
         state.m_info=1;
         //--- return result
         return(false);
        }
      //--- change values
      state.m_stplen=state.m_stplen/m_armijofactor;
      state.m_fcur=state.m_f;
      //--- function call, return result
      return(Func_lbl_10(state,n,v));
     }
//--- check
   if(state.m_rstate.stage==3)
     {
      state.m_nfev=state.m_nfev+1;
      //--- make decision
      if(state.m_f<state.m_fcur)
        {
         //--- change values
         state.m_stplen=state.m_stplen/m_armijofactor;
         state.m_fcur=state.m_f;
        }
      else
        {
         state.m_info=1;
         //--- return result
         return(false);
        }
      return(Func_lbl_10(state,n,v));
     }
//--- Routine body
   if((state.m_stplen<=0.0 || state.m_stpmax<0.0) || state.m_fmax<2)
     {
      state.m_info=0;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_stplen<=m_stpmin)
     {
      state.m_info=4;
      //--- return result
      return(false);
     }
//--- change values
   n=state.m_n;
   state.m_nfev=0;
//--- We always need F
   state.m_needf=true;
//--- Bound StpLen
   if(state.m_stplen>state.m_stpmax && state.m_stpmax!=0.0)
      state.m_stplen=state.m_stpmax;
//--- Increase length
   v=state.m_stplen*m_armijofactor;
//--- check
   if(v>state.m_stpmax && state.m_stpmax!=0.0)
      v=state.m_stpmax;
//--- copy
   for(i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
   for(i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_x[i_]+v*state.m_s[i_];
   state.m_rstate.stage=0;
//--- Saving state
   Func_lbl_rcomm(state,n,v);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for ArmijoIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static void CLinMin::Func_lbl_rcomm(CArmijoState &state,int n,double v)
  {
//--- save
   state.m_rstate.ia[0]=n;
   state.m_rstate.ra[0]=v;
  }
//+------------------------------------------------------------------+
//| Auxiliary function for ArmijoIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CLinMin::Func_lbl_6(CArmijoState &state,int &n,double &v)
  {
//--- test stopping conditions
   if(state.m_nfev>=state.m_fmax)
     {
      state.m_info=3;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_stplen>=state.m_stpmax)
     {
      state.m_info=5;
      //--- return result
      return(false);
     }
//--- evaluate F
   v=state.m_stplen*m_armijofactor;
//--- check
   if(v>state.m_stpmax && state.m_stpmax!=0.0)
      v=state.m_stpmax;
//--- copy
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_x[i_]+v*state.m_s[i_];
   state.m_rstate.stage=1;
//--- Saving state
   Func_lbl_rcomm(state,n,v);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for ArmijoIteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CLinMin::Func_lbl_10(CArmijoState &state,int &n,double &v)
  {
//--- test stopping conditions
   if(state.m_nfev>=state.m_fmax)
     {
      state.m_info=3;
      //--- return result
      return(false);
     }
//--- check
   if(state.m_stplen<=m_stpmin)
     {
      state.m_info=4;
      //--- return result
      return(false);
     }
//--- evaluate F
   v=state.m_stplen/m_armijofactor;
//--- copy
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_xbase[i_];
   for(int i_=0;i_<=n-1;i_++)
      state.m_x[i_]=state.m_x[i_]+v*state.m_s[i_];
   state.m_rstate.stage=3;
//--- Saving state
   Func_lbl_rcomm(state,n,v);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary functions for other classes                            |
//+------------------------------------------------------------------+
class COptServ
  {
public:
   //--- constructor, destructor
                     COptServ(void);
                    ~COptServ(void);
   //--- methods
   static void       TrimPrepare(const double f,double &threshold);
   static void       TrimFunction(double &f,double &g[],const int n,const double threshold);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
COptServ::COptServ(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
COptServ::~COptServ(void)
  {

  }
//+------------------------------------------------------------------+
//| This subroutine is used to prepare threshold value which will be |
//| used for trimming of the target function (see comments on        |
//| TrimFunction() for more information).                            |
//| This function accepts only one parameter: function value at the  |
//| starting point. It returns threshold which will be used for      |
//| trimming.                                                        |
//+------------------------------------------------------------------+
static void COptServ::TrimPrepare(const double f,double &threshold)
  {
//--- calculation
   threshold=10*(MathAbs(f)+1);
  }
//+------------------------------------------------------------------+
//| This subroutine is used to "trim" target function, i.e. to do    |
//| following transformation:                                        |
//|                    { {F,G}          if F<Threshold               |
//|     {F_tr, G_tr} = {                                             |
//|                    { {Threshold, 0} if F>=Threshold              |
//| Such transformation allows us to solve problems with             |
//| singularities by redefining function in such way that it becomes |
//| bounded from above.                                              |
//+------------------------------------------------------------------+
static void COptServ::TrimFunction(double &f,double &g[],const int n,
                                   const double threshold)
  {
//--- create a variable
   int i=0;
//--- check
   if(f>=threshold)
     {
      f=threshold;
      for(i=0;i<=n-1;i++)
         g[i]=0.0;
     }
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CFtBase                                      |
//+------------------------------------------------------------------+
class CFtPlan
  {
public:
   //--- arrays
   int               m_plan[];
   double            m_precomputed[];
   double            m_tmpbuf[];
   double            m_stackbuf[];
   //--- constructor, destructor
                     CFtPlan(void);
                    ~CFtPlan(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CFtPlan::CFtPlan(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CFtPlan::~CFtPlan(void)
  {

  }
//+------------------------------------------------------------------+
//| Generation FFT, FHT plans                                        |
//+------------------------------------------------------------------+
class CFtBase
  {
private:
   //--- private methods
   static void       FtBaseGeneratePlanRec(const int n,int tasktype,CFtPlan &plan,int &plansize,int &precomputedsize,int &planarraysize,int &tmpmemsize,int &stackmemsize,int stackptr);
   static void       FtBasePrecomputePlanRec(CFtPlan &plan,const int entryoffset,const int stackptr);
   static void       FFtTwCalc(double &a[],const int aoffset,const int n1,const int n2);
   static void       InternalComplexLinTranspose(double &a[],const int m,const int n,const int astart,double &buf[]);
   static void       InternalRealLinTranspose(double &a[],const int m,const int n,const int astart,double &buf[]);
   static void       FFtICLTRec(double &a[],const int astart,const int astride,double &b[],const int bstart,const int bstride,const int m,const int n);
   static void       FFtIRLTRec(double &a[],const int astart,const int astride,double &b[],const int bstart,const int bstride,const int m,const int n);
   static void       FtBaseFindSmoothRec(const int n,const int seed,const int leastfactor,int &best);
   static void       FFtArrayResize(int &a[],int &asize,const int newasize);
   static void       ReFFHt(double &a[],const int n,const int offs);
public:
   //--- class constants
   static const int  m_ftbaseplanentrysize;
   static const int  m_ftbasecffttask;
   static const int  m_ftbaserfhttask;
   static const int  m_ftbaserffttask;
   static const int  m_fftcooleytukeyplan;
   static const int  m_fftbluesteinplan;
   static const int  m_fftcodeletplan;
   static const int  m_fhtcooleytukeyplan;
   static const int  m_fhtcodeletplan;
   static const int  m_fftrealcooleytukeyplan;
   static const int  m_fftemptyplan;
   static const int  m_fhtn2plan;
   static const int  m_ftbaseupdatetw;
   static const int  m_ftbasecodeletrecommended;
   static const double m_ftbaseinefficiencyfactor;
   static const int  m_ftbasemaxsmoothfactor;
   //--- constructor, destructor
                     CFtBase(void);
                    ~CFtBase(void);
   //--- public methods
   static void       FtBaseGenerateComplexFFtPlan(const int n,CFtPlan &plan);
   static void       FtBaseGenerateRealFFtPlan(const int n,CFtPlan &plan);
   static void       FtBaseGenerateRealFHtPlan(const int n,CFtPlan &plan);
   static void       FtBaseExecutePlan(double &a[],const int aoffset,const int n,CFtPlan &plan);
   static void       FtBaseExecutePlanRec(double &a[],const int aoffset,CFtPlan &plan,const int entryoffset,const int stackptr);
   static void       FtBaseFactorize(const int n,const int tasktype,int &n1,int &n2);
   static bool       FtBaseIsSmooth(int n);
   static int        FtBaseFindSmooth(const int n);
   static int        FtBaseFindSmoothEven(const int n);
   static double     FtBaseGetFlopEstimate(const int n);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CFtBase::CFtBase(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CFtBase::~CFtBase(void)
  {

  }
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const int    CFtBase::m_ftbaseplanentrysize=8;
const int    CFtBase::m_ftbasecffttask=0;
const int    CFtBase::m_ftbaserfhttask=1;
const int    CFtBase::m_ftbaserffttask=2;
const int    CFtBase::m_fftcooleytukeyplan=0;
const int    CFtBase::m_fftbluesteinplan=1;
const int    CFtBase::m_fftcodeletplan=2;
const int    CFtBase::m_fhtcooleytukeyplan=3;
const int    CFtBase::m_fhtcodeletplan=4;
const int    CFtBase::m_fftrealcooleytukeyplan=5;
const int    CFtBase::m_fftemptyplan=6;
const int    CFtBase::m_fhtn2plan=999;
const int    CFtBase::m_ftbaseupdatetw=4;
const int    CFtBase::m_ftbasecodeletrecommended=5;
const double CFtBase::m_ftbaseinefficiencyfactor=1.3;
const int    CFtBase::m_ftbasemaxsmoothfactor=5;
//+------------------------------------------------------------------+
//| This subroutine generates FFT plan - a decomposition of a        |
//| N-length FFT to the more simpler operations. Plan consists of the|
//| root entry and the child entries.                                |
//| Subroutine parameters:                                           |
//|     N               task size                                    |
//| Output parameters:                                               |
//|     Plan            plan                                         |
//+------------------------------------------------------------------+
static void CFtBase::FtBaseGenerateComplexFFtPlan(const int n,CFtPlan &plan)
  {
//--- create variables
   int planarraysize=0;
   int plansize=0;
   int precomputedsize=0;
   int tmpmemsize=0;
   int stackmemsize=0;
   int stackptr=0;
//--- initialization
   planarraysize=1;
   plansize=0;
   precomputedsize=0;
   stackmemsize=0;
   stackptr=0;
   tmpmemsize=2*n;
//--- allocation
   ArrayResizeAL(plan.m_plan,planarraysize);
//--- function call
   FtBaseGeneratePlanRec(n,m_ftbasecffttask,plan,plansize,precomputedsize,planarraysize,tmpmemsize,stackmemsize,stackptr);
//--- check
   if(!CAp::Assert(stackptr==0,__FUNCTION__+": stack ptr!"))
      return;
//--- allocation
   ArrayResizeAL(plan.m_stackbuf,(int)MathMax(stackmemsize,1));
//--- allocation
   ArrayResizeAL(plan.m_tmpbuf,(int)MathMax(tmpmemsize,1));
//--- allocation
   ArrayResizeAL(plan.m_precomputed,(int)MathMax(precomputedsize,1));
   stackptr=0;
//--- function call
   FtBasePrecomputePlanRec(plan,0,stackptr);
//--- check
   if(!CAp::Assert(stackptr==0,__FUNCTION__+": stack ptr!"))
      return;
  }
//+------------------------------------------------------------------+
//| Generates real FFT plan                                          |
//+------------------------------------------------------------------+
static void CFtBase::FtBaseGenerateRealFFtPlan(const int n,CFtPlan &plan)
  {
//--- create variables
   int planarraysize=0;
   int plansize=0;
   int precomputedsize=0;
   int tmpmemsize=0;
   int stackmemsize=0;
   int stackptr=0;
//--- initialization
   planarraysize=1;
   plansize=0;
   precomputedsize=0;
   stackmemsize=0;
   stackptr=0;
   tmpmemsize=2*n;
//--- allocation
   ArrayResizeAL(plan.m_plan,planarraysize);
//--- function call
   FtBaseGeneratePlanRec(n,m_ftbaserffttask,plan,plansize,precomputedsize,planarraysize,tmpmemsize,stackmemsize,stackptr);
//--- check
   if(!CAp::Assert(stackptr==0,__FUNCTION__+": stack ptr!"))
      return;
//--- allocation
   ArrayResizeAL(plan.m_stackbuf,(int)MathMax(stackmemsize,1));
//--- allocation
   ArrayResizeAL(plan.m_tmpbuf,(int)MathMax(tmpmemsize,1));
//--- allocation
   ArrayResizeAL(plan.m_precomputed,(int)MathMax(precomputedsize,1));
   stackptr=0;
//--- function call
   FtBasePrecomputePlanRec(plan,0,stackptr);
//--- check
   if(!CAp::Assert(stackptr==0,__FUNCTION__+": stack ptr!"))
      return;
  }
//+------------------------------------------------------------------+
//| Generates real FHT plan                                          |
//+------------------------------------------------------------------+
static void CFtBase::FtBaseGenerateRealFHtPlan(const int n,CFtPlan &plan)
  {
//--- create variables
   int planarraysize=0;
   int plansize=0;
   int precomputedsize=0;
   int tmpmemsize=0;
   int stackmemsize=0;
   int stackptr=0;
//--- initialization
   planarraysize=1;
   plansize=0;
   precomputedsize=0;
   stackmemsize=0;
   stackptr=0;
   tmpmemsize=n;
//--- allocation
   ArrayResizeAL(plan.m_plan,planarraysize);
//--- function call
   FtBaseGeneratePlanRec(n,m_ftbaserfhttask,plan,plansize,precomputedsize,planarraysize,tmpmemsize,stackmemsize,stackptr);
//--- check
   if(!CAp::Assert(stackptr==0,__FUNCTION__+": stack ptr!"))
      return;
//--- allocation
   ArrayResizeAL(plan.m_stackbuf,(int)MathMax(stackmemsize,1));
//--- allocation
   ArrayResizeAL(plan.m_tmpbuf,(int)MathMax(tmpmemsize,1));
//--- allocation
   ArrayResizeAL(plan.m_precomputed,(int)MathMax(precomputedsize,1));
   stackptr=0;
//--- function call
   FtBasePrecomputePlanRec(plan,0,stackptr);
//--- check
   if(!CAp::Assert(stackptr==0,__FUNCTION__+": stack ptr!"))
      return;
  }
//+------------------------------------------------------------------+
//| This subroutine executes FFT/FHT plan.                           |
//| If Plan is a:                                                    |
//| * complex FFT plan  -   sizeof(A)=2*N,                           |
//|                         A contains interleaved real/imaginary    |
//|                         values                                   |
//| * real FFT plan     -   sizeof(A)=2*N,                           |
//|                         A contains real values interleaved with  |
//|                         zeros                                    |
//| * real FHT plan     -   sizeof(A)=2*N,                           |
//|                         A contains real values interleaved with  |
//|                         zeros                                    |
//+------------------------------------------------------------------+
static void CFtBase::FtBaseExecutePlan(double &a[],const int aoffset,
                                       const int n,CFtPlan &plan)
  {
//--- create a variable
   int stackptr=0;
//--- function call
   FtBaseExecutePlanRec(a,aoffset,plan,0,stackptr);
  }
//+------------------------------------------------------------------+
//| Recurrent subroutine for the FTBaseExecutePlan                   |
//| Parameters:                                                      |
//|     A           FFT'ed array                                     |
//|     AOffset     offset of the FFT'ed part (distance is measured  |
//|                 in doubles)                                      |
//+------------------------------------------------------------------+
static void CFtBase::FtBaseExecutePlanRec(double &a[],const int aoffset,
                                          CFtPlan &plan,const int entryoffset,
                                          const int stackptr)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    n1=0;
   int    n2=0;
   int    n=0;
   int    m=0;
   int    offs=0;
   int    offs1=0;
   int    offs2=0;
   int    offsa=0;
   int    offsb=0;
   int    offsp=0;
   double hk=0;
   double hnk=0;
   double x=0;
   double y=0;
   double bx=0;
   double by=0;
   double emptyarray[];
   double a0x=0;
   double a0y=0;
   double a1x=0;
   double a1y=0;
   double a2x=0;
   double a2y=0;
   double a3x=0;
   double a3y=0;
   double v0=0;
   double v1=0;
   double v2=0;
   double v3=0;
   double t1x=0;
   double t1y=0;
   double t2x=0;
   double t2y=0;
   double t3x=0;
   double t3y=0;
   double t4x=0;
   double t4y=0;
   double t5x=0;
   double t5y=0;
   double m1x=0;
   double m1y=0;
   double m2x=0;
   double m2y=0;
   double m3x=0;
   double m3y=0;
   double m4x=0;
   double m4y=0;
   double m5x=0;
   double m5y=0;
   double s1x=0;
   double s1y=0;
   double s2x=0;
   double s2y=0;
   double s3x=0;
   double s3y=0;
   double s4x=0;
   double s4y=0;
   double s5x=0;
   double s5y=0;
   double c1=0;
   double c2=0;
   double c3=0;
   double c4=0;
   double c5=0;
   double tmp[];
   int    i_=0;
   int    i1_=0;
//--- check
   if(plan.m_plan[entryoffset+3]==m_fftemptyplan)
      return;
//--- check
   if(plan.m_plan[entryoffset+3]==m_fftcooleytukeyplan)
     {
      //--- Cooley-Tukey plan
      //--- * transposition
      //--- * row-wise FFT
      //--- * twiddle factors:
      //---   - TwBase is a basis twiddle factor for I=1,J=1
      //---   - TwRow is a twiddle factor for a second element in a row (J=1)
      //---   - Tw is a twiddle factor for a current element
      //--- * transposition again
      //--- * row-wise FFT again
      n1=plan.m_plan[entryoffset+1];
      n2=plan.m_plan[entryoffset+2];
      //--- function call
      InternalComplexLinTranspose(a,n1,n2,aoffset,plan.m_tmpbuf);
      for(i=0;i<=n2-1;i++)
        {
         //--- function call
         FtBaseExecutePlanRec(a,aoffset+i*n1*2,plan,plan.m_plan[entryoffset+5],stackptr);
        }
      //--- function call
      FFtTwCalc(a,aoffset,n1,n2);
      //--- function call
      InternalComplexLinTranspose(a,n2,n1,aoffset,plan.m_tmpbuf);
      for(i=0;i<=n1-1;i++)
        {
         //--- function call
         FtBaseExecutePlanRec(a,aoffset+i*n2*2,plan,plan.m_plan[entryoffset+6],stackptr);
        }
      //--- function call
      InternalComplexLinTranspose(a,n1,n2,aoffset,plan.m_tmpbuf);
      //--- exit the function
      return;
     }
//--- check
   if(plan.m_plan[entryoffset+3]==m_fftrealcooleytukeyplan)
     {
      //--- Cooley-Tukey plan
      //--- * transposition
      //--- * row-wise FFT
      //--- * twiddle factors:
      //---   - TwBase is a basis twiddle factor for I=1,J=1
      //---   - TwRow is a twiddle factor for a second element in a row (J=1)
      //---   - Tw is a twiddle factor for a current element
      //--- * transposition again
      //--- * row-wise FFT again
      n1=plan.m_plan[entryoffset+1];
      n2=plan.m_plan[entryoffset+2];
      //--- function call
      InternalComplexLinTranspose(a,n2,n1,aoffset,plan.m_tmpbuf);
      for(i=0;i<=n1/2-1;i++)
        {
         //--- pack two adjacent smaller real FFT's together,
         //--- make one complex FFT,
         //--- unpack result
         offs=aoffset+2*i*n2*2;
         for(k=0;k<=n2-1;k++)
            a[offs+2*k+1]=a[offs+2*n2+2*k+0];
         //--- function call
         FtBaseExecutePlanRec(a,offs,plan,plan.m_plan[entryoffset+6],stackptr);
         //--- change values
         plan.m_tmpbuf[0]=a[offs+0];
         plan.m_tmpbuf[1]=0;
         plan.m_tmpbuf[2*n2+0]=a[offs+1];
         plan.m_tmpbuf[2*n2+1]=0;
         for(k=1;k<=n2-1;k++)
           {
            //--- change values
            offs1=2*k;
            offs2=2*n2+2*k;
            hk=a[offs+2*k+0];
            hnk=a[offs+2*(n2-k)+0];
            plan.m_tmpbuf[offs1+0]=0.5*(hk+hnk);
            plan.m_tmpbuf[offs2+1]=-(0.5*(hk-hnk));
            hk=a[offs+2*k+1];
            hnk=a[offs+2*(n2-k)+1];
            plan.m_tmpbuf[offs2+0]=0.5*(hk+hnk);
            plan.m_tmpbuf[offs1+1]=0.5*(hk-hnk);
           }
         i1_=-offs;
         for(i_=offs;i_<=offs+2*n2*2-1;i_++)
            a[i_]=plan.m_tmpbuf[i_+i1_];
        }
      //--- check
      if(n1%2!=0)
        {
         //--- function call
         FtBaseExecutePlanRec(a,aoffset+(n1-1)*n2*2,plan,plan.m_plan[entryoffset+6],stackptr);
        }
      //--- function call
      FFtTwCalc(a,aoffset,n2,n1);
      //--- function call
      InternalComplexLinTranspose(a,n1,n2,aoffset,plan.m_tmpbuf);
      for(i=0;i<=n2-1;i++)
        {
         //--- function call
         FtBaseExecutePlanRec(a,aoffset+i*n1*2,plan,plan.m_plan[entryoffset+5],stackptr);
        }
      //--- function call
      InternalComplexLinTranspose(a,n2,n1,aoffset,plan.m_tmpbuf);
      //--- exit the function
      return;
     }
//--- check
   if(plan.m_plan[entryoffset+3]==m_fhtcooleytukeyplan)
     {
      //--- Cooley-Tukey FHT plan:
      //--- * transpose                    \
      //--- * smaller FHT's                |
      //--- * pre-process                  |
      //--- * multiply by twiddle factors  | corresponds to multiplication by H1
      //--- * post-process                 |
      //--- * transpose again              /
      //--- * multiply by H2 (smaller FHT's)
      //--- * final transposition
      //--- For more details see Vitezslav Vesely,"Fast algorithms
      //--- of Fourier and Hartley transform and their implementation in MATLAB",
      //--- page 31.
      n1=plan.m_plan[entryoffset+1];
      n2=plan.m_plan[entryoffset+2];
      n=n1*n2;
      //--- function call
      InternalRealLinTranspose(a,n1,n2,aoffset,plan.m_tmpbuf);
      for(i=0;i<=n2-1;i++)
        {
         //--- function call
         FtBaseExecutePlanRec(a,aoffset+i*n1,plan,plan.m_plan[entryoffset+5],stackptr);
        }
      for(i=0;i<=n2-1;i++)
        {
         for(j=0;j<=n1-1;j++)
           {
            //--- change values
            offsa=aoffset+i*n1;
            hk=a[offsa+j];
            hnk=a[offsa+(n1-j)%n1];
            offs=2*(i*n1+j);
            plan.m_tmpbuf[offs+0]=-(0.5*(hnk-hk));
            plan.m_tmpbuf[offs+1]=0.5*(hk+hnk);
           }
        }
      //--- function call
      FFtTwCalc(plan.m_tmpbuf,0,n1,n2);
      for(j=0;j<=n1-1;j++)
         a[aoffset+j]=plan.m_tmpbuf[2*j+0]+plan.m_tmpbuf[2*j+1];
      //--- check
      if(n2%2==0)
        {
         offs=2*(n2/2)*n1;
         offsa=aoffset+n2/2*n1;
         for(j=0;j<=n1-1;j++)
            a[offsa+j]=plan.m_tmpbuf[offs+2*j+0]+plan.m_tmpbuf[offs+2*j+1];
        }
      for(i=1;i<=(n2+1)/2-1;i++)
        {
         //--- change values
         offs=2*i*n1;
         offs2=2*(n2-i)*n1;
         offsa=aoffset+i*n1;
         for(j=0;j<=n1-1;j++)
            a[offsa+j]=plan.m_tmpbuf[offs+2*j+1]+plan.m_tmpbuf[offs2+2*j+0];
         offsa=aoffset+(n2-i)*n1;
         for(j=0;j<=n1-1;j++)
            a[offsa+j]=plan.m_tmpbuf[offs+2*j+0]+plan.m_tmpbuf[offs2+2*j+1];
        }
      //--- function call
      InternalRealLinTranspose(a,n2,n1,aoffset,plan.m_tmpbuf);
      for(i=0;i<=n1-1;i++)
        {
         //--- function call
         FtBaseExecutePlanRec(a,aoffset+i*n2,plan,plan.m_plan[entryoffset+6],stackptr);
        }
      //--- function call
      InternalRealLinTranspose(a,n1,n2,aoffset,plan.m_tmpbuf);
      //--- exit the function
      return;
     }
//--- check
   if(plan.m_plan[entryoffset+3]==m_fhtn2plan)
     {
      //--- Cooley-Tukey FHT plan
      n1=plan.m_plan[entryoffset+1];
      n2=plan.m_plan[entryoffset+2];
      n=n1*n2;
      ReFFHt(a,n,aoffset);
      //--- exit the function
      return;
     }
//--- check
   if(plan.m_plan[entryoffset+3]==m_fftcodeletplan)
     {
      n1=plan.m_plan[entryoffset+1];
      n2=plan.m_plan[entryoffset+2];
      n=n1*n2;
      //--- check
      if(n==2)
        {
         //--- change values
         a0x=a[aoffset+0];
         a0y=a[aoffset+1];
         a1x=a[aoffset+2];
         a1y=a[aoffset+3];
         v0=a0x+a1x;
         v1=a0y+a1y;
         v2=a0x-a1x;
         v3=a0y-a1y;
         a[aoffset+0]=v0;
         a[aoffset+1]=v1;
         a[aoffset+2]=v2;
         a[aoffset+3]=v3;
         //--- exit the function
         return;
        }
      //--- check
      if(n==3)
        {
         //--- change values
         offs=plan.m_plan[entryoffset+7];
         c1=plan.m_precomputed[offs+0];
         c2=plan.m_precomputed[offs+1];
         a0x=a[aoffset+0];
         a0y=a[aoffset+1];
         a1x=a[aoffset+2];
         a1y=a[aoffset+3];
         a2x=a[aoffset+4];
         a2y=a[aoffset+5];
         t1x=a1x+a2x;
         t1y=a1y+a2y;
         a0x=a0x+t1x;
         a0y=a0y+t1y;
         m1x=c1*t1x;
         m1y=c1*t1y;
         m2x=c2*(a1y-a2y);
         m2y=c2*(a2x-a1x);
         s1x=a0x+m1x;
         s1y=a0y+m1y;
         a1x=s1x+m2x;
         a1y=s1y+m2y;
         a2x=s1x-m2x;
         a2y=s1y-m2y;
         a[aoffset+0]=a0x;
         a[aoffset+1]=a0y;
         a[aoffset+2]=a1x;
         a[aoffset+3]=a1y;
         a[aoffset+4]=a2x;
         a[aoffset+5]=a2y;
         //--- exit the function
         return;
        }
      //--- check
      if(n==4)
        {
         //--- change values
         a0x=a[aoffset+0];
         a0y=a[aoffset+1];
         a1x=a[aoffset+2];
         a1y=a[aoffset+3];
         a2x=a[aoffset+4];
         a2y=a[aoffset+5];
         a3x=a[aoffset+6];
         a3y=a[aoffset+7];
         t1x=a0x+a2x;
         t1y=a0y+a2y;
         t2x=a1x+a3x;
         t2y=a1y+a3y;
         m2x=a0x-a2x;
         m2y=a0y-a2y;
         m3x=a1y-a3y;
         m3y=a3x-a1x;
         a[aoffset+0]=t1x+t2x;
         a[aoffset+1]=t1y+t2y;
         a[aoffset+4]=t1x-t2x;
         a[aoffset+5]=t1y-t2y;
         a[aoffset+2]=m2x+m3x;
         a[aoffset+3]=m2y+m3y;
         a[aoffset+6]=m2x-m3x;
         a[aoffset+7]=m2y-m3y;
         //--- exit the function
         return;
        }
      //--- check
      if(n==5)
        {
         //--- change values
         offs=plan.m_plan[entryoffset+7];
         c1=plan.m_precomputed[offs+0];
         c2=plan.m_precomputed[offs+1];
         c3=plan.m_precomputed[offs+2];
         c4=plan.m_precomputed[offs+3];
         c5=plan.m_precomputed[offs+4];
         t1x=a[aoffset+2]+a[aoffset+8];
         t1y=a[aoffset+3]+a[aoffset+9];
         t2x=a[aoffset+4]+a[aoffset+6];
         t2y=a[aoffset+5]+a[aoffset+7];
         t3x=a[aoffset+2]-a[aoffset+8];
         t3y=a[aoffset+3]-a[aoffset+9];
         t4x=a[aoffset+6]-a[aoffset+4];
         t4y=a[aoffset+7]-a[aoffset+5];
         t5x=t1x+t2x;
         t5y=t1y+t2y;
         a[aoffset+0]=a[aoffset+0]+t5x;
         a[aoffset+1]=a[aoffset+1]+t5y;
         m1x=c1*t5x;
         m1y=c1*t5y;
         m2x=c2*(t1x-t2x);
         m2y=c2*(t1y-t2y);
         m3x=-(c3*(t3y+t4y));
         m3y=c3*(t3x+t4x);
         m4x=-(c4*t4y);
         m4y=c4*t4x;
         m5x=-(c5*t3y);
         m5y=c5*t3x;
         s3x=m3x-m4x;
         s3y=m3y-m4y;
         s5x=m3x+m5x;
         s5y=m3y+m5y;
         s1x=a[aoffset+0]+m1x;
         s1y=a[aoffset+1]+m1y;
         s2x=s1x+m2x;
         s2y=s1y+m2y;
         s4x=s1x-m2x;
         s4y=s1y-m2y;
         a[aoffset+2]=s2x+s3x;
         a[aoffset+3]=s2y+s3y;
         a[aoffset+4]=s4x+s5x;
         a[aoffset+5]=s4y+s5y;
         a[aoffset+6]=s4x-s5x;
         a[aoffset+7]=s4y-s5y;
         a[aoffset+8]=s2x-s3x;
         a[aoffset+9]=s2y-s3y;
         //--- exit the function
         return;
        }
     }
//--- check
   if(plan.m_plan[entryoffset+3]==m_fhtcodeletplan)
     {
      //--- change values
      n1=plan.m_plan[entryoffset+1];
      n2=plan.m_plan[entryoffset+2];
      n=n1*n2;
      //--- check
      if(n==2)
        {
         //--- change values
         a0x=a[aoffset+0];
         a1x=a[aoffset+1];
         a[aoffset+0]=a0x+a1x;
         a[aoffset+1]=a0x-a1x;
         //--- exit the function
         return;
        }
      //--- check
      if(n==3)
        {
         //--- change values
         offs=plan.m_plan[entryoffset+7];
         c1=plan.m_precomputed[offs+0];
         c2=plan.m_precomputed[offs+1];
         a0x=a[aoffset+0];
         a1x=a[aoffset+1];
         a2x=a[aoffset+2];
         t1x=a1x+a2x;
         a0x=a0x+t1x;
         m1x=c1*t1x;
         m2y=c2*(a2x-a1x);
         s1x=a0x+m1x;
         a[aoffset+0]=a0x;
         a[aoffset+1]=s1x-m2y;
         a[aoffset+2]=s1x+m2y;
         //--- exit the function
         return;
        }
      //--- check
      if(n==4)
        {
         //--- change values
         a0x=a[aoffset+0];
         a1x=a[aoffset+1];
         a2x=a[aoffset+2];
         a3x=a[aoffset+3];
         t1x=a0x+a2x;
         t2x=a1x+a3x;
         m2x=a0x-a2x;
         m3y=a3x-a1x;
         a[aoffset+0]=t1x+t2x;
         a[aoffset+1]=m2x-m3y;
         a[aoffset+2]=t1x-t2x;
         a[aoffset+3]=m2x+m3y;
         //--- exit the function
         return;
        }
      //--- check
      if(n==5)
        {
         //--- change values
         offs=plan.m_plan[entryoffset+7];
         c1=plan.m_precomputed[offs+0];
         c2=plan.m_precomputed[offs+1];
         c3=plan.m_precomputed[offs+2];
         c4=plan.m_precomputed[offs+3];
         c5=plan.m_precomputed[offs+4];
         t1x=a[aoffset+1]+a[aoffset+4];
         t2x=a[aoffset+2]+a[aoffset+3];
         t3x=a[aoffset+1]-a[aoffset+4];
         t4x=a[aoffset+3]-a[aoffset+2];
         t5x=t1x+t2x;
         v0=a[aoffset+0]+t5x;
         a[aoffset+0]=v0;
         m2x=c2*(t1x-t2x);
         m3y=c3*(t3x+t4x);
         s3y=m3y-c4*t4x;
         s5y=m3y+c5*t3x;
         s1x=v0+c1*t5x;
         s2x=s1x+m2x;
         s4x=s1x-m2x;
         a[aoffset+1]=s2x-s3y;
         a[aoffset+2]=s4x-s5y;
         a[aoffset+3]=s4x+s5y;
         a[aoffset+4]=s2x+s3y;
         //--- exit the function
         return;
        }
     }
//--- check
   if(plan.m_plan[entryoffset+3]==m_fftbluesteinplan)
     {
      //--- Bluestein plan:
      //--- 1. multiply by precomputed coefficients
      //--- 2. make convolution: forward FFT,multiplication by precomputed FFT
      //---    and backward FFT. backward FFT is represented as
      //---        invfft(x)=fft(x')'/M
      //---    for performance reasons reduction of inverse FFT to
      //---    forward FFT is merged with multiplication of FFT components
      //---    and last stage of Bluestein's transformation.
      //--- 3. post-multiplication by Bluestein factors
      n=plan.m_plan[entryoffset+1];
      m=plan.m_plan[entryoffset+4];
      offs=plan.m_plan[entryoffset+7];
      for(i=stackptr+2*n;i<=stackptr+2*m-1;i++)
         plan.m_stackbuf[i]=0;
      //--- change values
      offsp=offs+2*m;
      offsa=aoffset;
      offsb=stackptr;
      for(i=0;i<=n-1;i++)
        {
         //--- change values
         bx=plan.m_precomputed[offsp+0];
         by=plan.m_precomputed[offsp+1];
         x=a[offsa+0];
         y=a[offsa+1];
         plan.m_stackbuf[offsb+0]=x*bx-y*-by;
         plan.m_stackbuf[offsb+1]=x*-by+y*bx;
         offsp=offsp+2;
         offsa=offsa+2;
         offsb=offsb+2;
        }
      //--- function call
      FtBaseExecutePlanRec(plan.m_stackbuf,stackptr,plan,plan.m_plan[entryoffset+5],stackptr+2*2*m);
      offsb=stackptr;
      offsp=offs;
      for(i=0;i<=m-1;i++)
        {
         //--- change values
         x=plan.m_stackbuf[offsb+0];
         y=plan.m_stackbuf[offsb+1];
         bx=plan.m_precomputed[offsp+0];
         by=plan.m_precomputed[offsp+1];
         plan.m_stackbuf[offsb+0]=x*bx-y*by;
         plan.m_stackbuf[offsb+1]=-(x*by+y*bx);
         offsb=offsb+2;
         offsp=offsp+2;
        }
      //--- function call
      FtBaseExecutePlanRec(plan.m_stackbuf,stackptr,plan,plan.m_plan[entryoffset+5],stackptr+2*2*m);
      offsb=stackptr;
      offsp=offs+2*m;
      offsa=aoffset;
      for(i=0;i<=n-1;i++)
        {
         //--- change values
         x=plan.m_stackbuf[offsb+0]/m;
         y=-(plan.m_stackbuf[offsb+1]/m);
         bx=plan.m_precomputed[offsp+0];
         by=plan.m_precomputed[offsp+1];
         a[offsa+0]=x*bx-y*-by;
         a[offsa+1]=x*-by+y*bx;
         offsp=offsp+2;
         offsa=offsa+2;
         offsb=offsb+2;
        }
      //--- exit the function
      return;
     }
  }
//+------------------------------------------------------------------+
//| Returns good factorization N=N1*N2.                              |
//| Usually N1<=N2 (but not always - small N's may be exception).    |
//| if N1<>1 then N2<>1.                                             |
//| Factorization is chosen depending on task type and codelets we   |
//| have.                                                            |
//+------------------------------------------------------------------+
static void CFtBase::FtBaseFactorize(const int n,const int tasktype,
                                     int &n1,int &n2)
  {
//--- create a variable
   int j=0;
//--- initialization
   n1=0;
   n2=0;
//--- try to find good codelet
   if(n1*n2!=n)
     {
      for(j=m_ftbasecodeletrecommended;j>=2;j--)
        {
         //--- check
         if(n%j==0)
           {
            n1=j;
            n2=n/j;
            break;
           }
        }
     }
//--- try to factorize N
   if(n1*n2!=n)
     {
      for(j=m_ftbasecodeletrecommended+1;j<=n-1;j++)
        {
         //--- check
         if(n%j==0)
           {
            n1=j;
            n2=n/j;
            break;
           }
        }
     }
//--- looks like N is prime :(
   if(n1*n2!=n)
     {
      n1=1;
      n2=n;
     }
//--- normalize
   if(n2==1 && n1!=1)
     {
      n2=n1;
      n1=1;
     }
  }
//+------------------------------------------------------------------+
//| Is number smooth?                                                |
//+------------------------------------------------------------------+
static bool CFtBase::FtBaseIsSmooth(int n)
  {
//--- create a variable
   int i=0;
//--- change n
   for(i=2;i<=m_ftbasemaxsmoothfactor;i++)
     {
      while(n%i==0)
         n=n/i;
     }
//--- check
   if(n==1)
      return(true);
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Returns smallest smooth (divisible only by 2, 3, 5) number that  |
//| is greater than or equal to max(N,2)                             |
//+------------------------------------------------------------------+
static int CFtBase::FtBaseFindSmooth(const int n)
  {
//--- create a variable
   int best=0;
//--- calculation
   best=2;
   while(best<n)
      best=2*best;
//--- function call
   FtBaseFindSmoothRec(n,1,2,best);
//--- return result
   return(best);
  }
//+------------------------------------------------------------------+
//| Returns smallest smooth (divisible only by 2, 3, 5) even number  |
//| that is greater than or equal to max(N,2)                        |
//+------------------------------------------------------------------+
static int CFtBase::FtBaseFindSmoothEven(const int n)
  {
//--- create a variable
   int best=0;
//--- calculation
   best=2;
   while(best<n)
      best=2*best;
//--- function call
   FtBaseFindSmoothRec(n,2,2,best);
//--- return result
   return(best);
  }
//+------------------------------------------------------------------+
//| Returns estimate of FLOP count for the FFT.                      |
//| It is only an estimate based on operations count for the PERFECT |
//| FFT and relative inefficiency of the algorithm actually used.    |
//| N should be power of 2, estimates are badly wrong for            |
//| non-power-of-2 N's.                                              |
//+------------------------------------------------------------------+
static double CFtBase::FtBaseGetFlopEstimate(const int n)
  {
//--- return result
   return(m_ftbaseinefficiencyfactor*(4*n*MathLog(n)/MathLog(2)-6*n+8));
  }
//+------------------------------------------------------------------+
//| Recurrent subroutine for the FFTGeneratePlan:                    |
//| PARAMETERS:                                                      |
//|     N                   plan size                                |
//|     IsReal              whether input is real or not.            |
//|                         subroutine MUST NOT ignore this flag     |
//|                         because real inputs comes with           |
//|                         non-initialized imaginary parts, so      |
//|                         ignoring this flag will result in        |
//|                         corrupted output                         |
//|     HalfOut             whether full output or only half of it   |
//|                         from 0 to floor(N/2) is needed. This flag|
//|                         may be ignored if doing so will simplify |
//|                         calculations                             |
//|     Plan                plan array                               |
//|     PlanSize            size of used part (in integers)          |
//|     PrecomputedSize     size of precomputed array allocated yet  |
//|     PlanArraySize       plan array size (actual)                 |
//|     TmpMemSize          temporary memory required size           |
//|     BluesteinMemSize    temporary memory required size           |
//+------------------------------------------------------------------+
static void CFtBase::FtBaseGeneratePlanRec(const int n,int tasktype,
                                           CFtPlan &plan,int &plansize,
                                           int &precomputedsize,
                                           int &planarraysize,
                                           int &tmpmemsize,
                                           int &stackmemsize,int stackptr)
  {
//--- create variables
   int k=0;
   int m=0;
   int n1=0;
   int n2=0;
   int esize=0;
   int entryoffset=0;
//--- prepare
   if(plansize+m_ftbaseplanentrysize>planarraysize)
      FFtArrayResize(plan.m_plan,planarraysize,8*planarraysize);
   entryoffset=plansize;
   esize=m_ftbaseplanentrysize;
   plansize=plansize+esize;
//--- if N=1,generate empty plan and exit
   if(n==1)
     {
      //--- change values
      plan.m_plan[entryoffset+0]=esize;
      plan.m_plan[entryoffset+1]=-1;
      plan.m_plan[entryoffset+2]=-1;
      plan.m_plan[entryoffset+3]=m_fftemptyplan;
      plan.m_plan[entryoffset+4]=-1;
      plan.m_plan[entryoffset+5]=-1;
      plan.m_plan[entryoffset+6]=-1;
      plan.m_plan[entryoffset+7]=-1;
      //--- exit the function
      return;
     }
//--- generate plans
   FtBaseFactorize(n,tasktype,n1,n2);
//--- check
   if(tasktype==m_ftbasecffttask || tasktype==m_ftbaserffttask)
     {
      //--- complex FFT plans
      if(n1!=1)
        {
         //--- Cooley-Tukey plan (real or complex)
         //--- Note that child plans are COMPLEX
         //--- (whether plan itself is complex or not).
         tmpmemsize=MathMax(tmpmemsize,2*n1*n2);
         plan.m_plan[entryoffset+0]=esize;
         plan.m_plan[entryoffset+1]=n1;
         plan.m_plan[entryoffset+2]=n2;
         //--- check
         if(tasktype==m_ftbasecffttask)
            plan.m_plan[entryoffset+3]=m_fftcooleytukeyplan;
         else
            plan.m_plan[entryoffset+3]=m_fftrealcooleytukeyplan;
         plan.m_plan[entryoffset+4]=0;
         plan.m_plan[entryoffset+5]=plansize;
         //--- function call
         FtBaseGeneratePlanRec(n1,m_ftbasecffttask,plan,plansize,precomputedsize,planarraysize,tmpmemsize,stackmemsize,stackptr);
         plan.m_plan[entryoffset+6]=plansize;
         //--- function call
         FtBaseGeneratePlanRec(n2,m_ftbasecffttask,plan,plansize,precomputedsize,planarraysize,tmpmemsize,stackmemsize,stackptr);
         plan.m_plan[entryoffset+7]=-1;
         //--- exit the function
         return;
        }
      else
        {
         //--- check
         if(((n==2 || n==3) || n==4) || n==5)
           {
            //--- hard-coded plan
            plan.m_plan[entryoffset+0]=esize;
            plan.m_plan[entryoffset+1]=n1;
            plan.m_plan[entryoffset+2]=n2;
            plan.m_plan[entryoffset+3]=m_fftcodeletplan;
            plan.m_plan[entryoffset+4]=0;
            plan.m_plan[entryoffset+5]=-1;
            plan.m_plan[entryoffset+6]=-1;
            plan.m_plan[entryoffset+7]=precomputedsize;
            //--- check
            if(n==3)
               precomputedsize=precomputedsize+2;
            //--- check
            if(n==5)
               precomputedsize=precomputedsize+5;
            //--- exit the function
            return;
           }
         else
           {
            //--- Bluestein's plan
            //--- Select such M that M>=2*N-1,M is composite,and M's
            //--- factors are 2,3,5
            k=2*n2-1;
            m=FtBaseFindSmooth(k);
            tmpmemsize=MathMax(tmpmemsize,2*m);
            plan.m_plan[entryoffset+0]=esize;
            plan.m_plan[entryoffset+1]=n2;
            plan.m_plan[entryoffset+2]=-1;
            plan.m_plan[entryoffset+3]=m_fftbluesteinplan;
            plan.m_plan[entryoffset+4]=m;
            plan.m_plan[entryoffset+5]=plansize;
            stackptr=stackptr+2*2*m;
            stackmemsize=MathMax(stackmemsize,stackptr);
            //--- function call
            FtBaseGeneratePlanRec(m,m_ftbasecffttask,plan,plansize,precomputedsize,planarraysize,tmpmemsize,stackmemsize,stackptr);
            stackptr=stackptr-2*2*m;
            plan.m_plan[entryoffset+6]=-1;
            plan.m_plan[entryoffset+7]=precomputedsize;
            precomputedsize=precomputedsize+2*m+2*n;
            //--- exit the function
            return;
           }
        }
     }
//--- check
   if(tasktype==m_ftbaserfhttask)
     {
      //--- real FHT plans
      if(n1!=1)
        {
         //--- Cooley-Tukey plan
         tmpmemsize=MathMax(tmpmemsize,2*n1*n2);
         plan.m_plan[entryoffset+0]=esize;
         plan.m_plan[entryoffset+1]=n1;
         plan.m_plan[entryoffset+2]=n2;
         plan.m_plan[entryoffset+3]=m_fhtcooleytukeyplan;
         plan.m_plan[entryoffset+4]=0;
         plan.m_plan[entryoffset+5]=plansize;
         //--- function call
         FtBaseGeneratePlanRec(n1,tasktype,plan,plansize,precomputedsize,planarraysize,tmpmemsize,stackmemsize,stackptr);
         plan.m_plan[entryoffset+6]=plansize;
         //--- function call
         FtBaseGeneratePlanRec(n2,tasktype,plan,plansize,precomputedsize,planarraysize,tmpmemsize,stackmemsize,stackptr);
         plan.m_plan[entryoffset+7]=-1;
         //--- exit the function
         return;
        }
      else
        {
         //--- N2 plan
         plan.m_plan[entryoffset+0]=esize;
         plan.m_plan[entryoffset+1]=n1;
         plan.m_plan[entryoffset+2]=n2;
         plan.m_plan[entryoffset+3]=m_fhtn2plan;
         plan.m_plan[entryoffset+4]=0;
         plan.m_plan[entryoffset+5]=-1;
         plan.m_plan[entryoffset+6]=-1;
         plan.m_plan[entryoffset+7]=-1;
         //--- check
         if(((n==2 || n==3) || n==4) || n==5)
           {
            //--- hard-coded plan
            plan.m_plan[entryoffset+0]=esize;
            plan.m_plan[entryoffset+1]=n1;
            plan.m_plan[entryoffset+2]=n2;
            plan.m_plan[entryoffset+3]=m_fhtcodeletplan;
            plan.m_plan[entryoffset+4]=0;
            plan.m_plan[entryoffset+5]=-1;
            plan.m_plan[entryoffset+6]=-1;
            plan.m_plan[entryoffset+7]=precomputedsize;
            //--- check
            if(n==3)
               precomputedsize=precomputedsize+2;
            //--- check
            if(n==5)
               precomputedsize=precomputedsize+5;
            //--- exit the function
            return;
           }
         //--- exit the function
         return;
        }
     }
  }
//+------------------------------------------------------------------+
//| Recurrent subroutine for precomputing FFT plans                  |
//+------------------------------------------------------------------+
static void CFtBase::FtBasePrecomputePlanRec(CFtPlan &plan,const int entryoffset,
                                             const int stackptr)
  {
//--- create variables
   int    i=0;
   int    n1=0;
   int    n2=0;
   int    n=0;
   int    m=0;
   int    offs=0;
   double v=0;
   double emptyarray[];
   double bx=0;
   double by=0;
//--- check
   if((plan.m_plan[entryoffset+3]==m_fftcooleytukeyplan || plan.m_plan[entryoffset+3]==m_fftrealcooleytukeyplan) || plan.m_plan[entryoffset+3]==m_fhtcooleytukeyplan)
     {
      //--- function call
      FtBasePrecomputePlanRec(plan,plan.m_plan[entryoffset+5],stackptr);
      //--- function call
      FtBasePrecomputePlanRec(plan,plan.m_plan[entryoffset+6],stackptr);
      //--- exit the function
      return;
     }
//--- check
   if(plan.m_plan[entryoffset+3]==m_fftcodeletplan || plan.m_plan[entryoffset+3]==m_fhtcodeletplan)
     {
      n1=plan.m_plan[entryoffset+1];
      n2=plan.m_plan[entryoffset+2];
      n=n1*n2;
      //--- check
      if(n==3)
        {
         offs=plan.m_plan[entryoffset+7];
         plan.m_precomputed[offs+0]=MathCos(2*M_PI/3)-1;
         plan.m_precomputed[offs+1]=MathSin(2*M_PI/3);
         //--- exit the function
         return;
        }
      //--- check
      if(n==5)
        {
         offs=plan.m_plan[entryoffset+7];
         v=2*M_PI/5;
         plan.m_precomputed[offs+0]=(MathCos(v)+MathCos(2*v))/2-1;
         plan.m_precomputed[offs+1]=(MathCos(v)-MathCos(2*v))/2;
         plan.m_precomputed[offs+2]=-MathSin(v);
         plan.m_precomputed[offs+3]=-(MathSin(v)+MathSin(2*v));
         plan.m_precomputed[offs+4]=MathSin(v)-MathSin(2*v);
         //--- exit the function
         return;
        }
     }
//--- check
   if(plan.m_plan[entryoffset+3]==m_fftbluesteinplan)
     {
      //--- function call
      FtBasePrecomputePlanRec(plan,plan.m_plan[entryoffset+5],stackptr);
      n=plan.m_plan[entryoffset+1];
      m=plan.m_plan[entryoffset+4];
      offs=plan.m_plan[entryoffset+7];
      for(i=0;i<=2*m-1;i++)
         plan.m_precomputed[offs+i]=0;
      //--- change values
      for(i=0;i<=n-1;i++)
        {
         bx=MathCos(M_PI*CMath::Sqr(i)/n);
         by=MathSin(M_PI*CMath::Sqr(i)/n);
         plan.m_precomputed[offs+2*i+0]=bx;
         plan.m_precomputed[offs+2*i+1]=by;
         plan.m_precomputed[offs+2*m+2*i+0]=bx;
         plan.m_precomputed[offs+2*m+2*i+1]=by;
         //--- check
         if(i>0)
           {
            plan.m_precomputed[offs+2*(m-i)+0]=bx;
            plan.m_precomputed[offs+2*(m-i)+1]=by;
           }
        }
      //--- function call
      FtBaseExecutePlanRec(plan.m_precomputed,offs,plan,plan.m_plan[entryoffset+5],stackptr);
      //--- exit the function
      return;
     }
  }
//+------------------------------------------------------------------+
//| Twiddle factors calculation                                      |
//+------------------------------------------------------------------+
static void CFtBase::FFtTwCalc(double &a[],const int aoffset,const int n1,
                               const int n2)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    n=0;
   int    idx=0;
   int    offs=0;
   double x=0;
   double y=0;
   double twxm1=0;
   double twy=0;
   double twbasexm1=0;
   double twbasey=0;
   double twrowxm1=0;
   double twrowy=0;
   double tmpx=0;
   double tmpy=0;
   double v=0;
//--- initialization
   n=n1*n2;
   v=-(2*M_PI/n);
   twbasexm1=-(2*CMath::Sqr(MathSin(0.5*v)));
   twbasey=MathSin(v);
   twrowxm1=0;
   twrowy=0;
//--- calculation
   for(i=0;i<=n2-1;i++)
     {
      twxm1=0;
      twy=0;
      for(j=0;j<=n1-1;j++)
        {
         //--- change values
         idx=i*n1+j;
         offs=aoffset+2*idx;
         x=a[offs+0];
         y=a[offs+1];
         tmpx=x*twxm1-y*twy;
         tmpy=x*twy+y*twxm1;
         a[offs+0]=x+tmpx;
         a[offs+1]=y+tmpy;
         //--- update Tw: Tw(new)=Tw(old)*TwRow
         if(j<n1-1)
           {
            //--- check
            if(j%m_ftbaseupdatetw==0)
              {
               v=-(2*M_PI*i*(j+1)/n);
               twxm1=-(2*CMath::Sqr(MathSin(0.5*v)));
               twy=MathSin(v);
              }
            else
              {
               tmpx=twrowxm1+twxm1*twrowxm1-twy*twrowy;
               tmpy=twrowy+twxm1*twrowy+twy*twrowxm1;
               twxm1=twxm1+tmpx;
               twy=twy+tmpy;
              }
           }
        }
      //--- update TwRow: TwRow(new)=TwRow(old)*TwBase
      if(i<n2-1)
        {
         //--- check
         if(j%m_ftbaseupdatetw==0)
           {
            v=-(2*M_PI*(i+1)/n);
            twrowxm1=-(2*CMath::Sqr(MathSin(0.5*v)));
            twrowy=MathSin(v);
           }
         else
           {
            tmpx=twbasexm1+twrowxm1*twbasexm1-twrowy*twbasey;
            tmpy=twbasey+twrowxm1*twbasey+twrowy*twbasexm1;
            twrowxm1=twrowxm1+tmpx;
            twrowy=twrowy+tmpy;
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Linear transpose: transpose complex matrix stored in             |
//| 1-dimensional array                                              |
//+------------------------------------------------------------------+
static void CFtBase::InternalComplexLinTranspose(double &a[],const int m,
                                                 const int n,const int astart,
                                                 double &buf[])
  {
//--- create variables
   int i_=0;
   int i1_=0;
//--- function call
   FFtICLTRec(a,astart,n,buf,0,m,m,n);
   i1_=-astart;
//--- calculation
   for(i_=astart;i_<=astart+2*m*n-1;i_++)
      a[i_]=buf[i_+i1_];
  }
//+------------------------------------------------------------------+
//| Linear transpose: transpose real matrix stored in 1-dimensional  |
//| array                                                            |
//+------------------------------------------------------------------+
static void CFtBase::InternalRealLinTranspose(double &a[],const int m,
                                              const int n,const int astart,
                                              double &buf[])
  {
//--- create variables
   int i_=0;
   int i1_=0;
//--- function call
   FFtIRLTRec(a,astart,n,buf,0,m,m,n);
   i1_=-astart;
//--- calculation
   for(i_=astart;i_<=astart+m*n-1;i_++)
      a[i_]=buf[i_+i1_];
  }
//+------------------------------------------------------------------+
//| Recurrent subroutine for a InternalComplexLinTranspose           |
//| Write A^T to B, where:                                           |
//| * A is m*n complex matrix stored in array A as pairs of          |
//|   real/image values, beginning from AStart position, with AStride|
//|   stride                                                         |
//| * B is n*m complex matrix stored in array B as pairs of          |
//|   real/image values, beginning from BStart position, with BStride|
//|   stride                                                         |
//| stride is measured in complex numbers, i.e. in real/image pairs. |
//+------------------------------------------------------------------+
static void CFtBase::FFtICLTRec(double &a[],const int astart,const int astride,
                                double &b[],const int bstart,const int bstride,
                                const int m,const int n)
  {
//--- create variables
   int i=0;
   int j=0;
   int idx1=0;
   int idx2=0;
   int m2=0;
   int m1=0;
   int n1=0;
//--- check
   if(m==0 || n==0)
      return;
//--- check
   if(MathMax(m,n)<=8)
     {
      m2=2*bstride;
      for(i=0;i<=m-1;i++)
        {
         //--- calculation
         idx1=bstart+2*i;
         idx2=astart+2*i*astride;
         for(j=0;j<=n-1;j++)
           {
            b[idx1+0]=a[idx2+0];
            b[idx1+1]=a[idx2+1];
            idx1=idx1+m2;
            idx2=idx2+2;
           }
        }
      //--- exit the function
      return;
     }
//--- check
   if(n>m)
     {
      //--- New partition:
      //--- "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
      //---                                  ( B2 )
      n1=n/2;
      //--- check
      if(n-n1>=8 && n1%8!=0)
         n1=n1+(8-n1%8);
      //--- check
      if(!CAp::Assert(n-n1>0))
         return;
      //--- function call
      FFtICLTRec(a,astart,astride,b,bstart,bstride,m,n1);
      //--- function call
      FFtICLTRec(a,astart+2*n1,astride,b,bstart+2*n1*bstride,bstride,m,n-n1);
     }
   else
     {
      //--- New partition:
      //--- "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
      //---                     ( A2 )
      m1=m/2;
      //--- check
      if(m-m1>=8 && m1%8!=0)
         m1=m1+(8-m1%8);
      //--- check
      if(!CAp::Assert(m-m1>0))
         return;
      //--- function call
      FFtICLTRec(a,astart,astride,b,bstart,bstride,m1,n);
      //--- function call
      FFtICLTRec(a,astart+2*m1*astride,astride,b,bstart+2*m1,bstride,m-m1,n);
     }
  }
//+------------------------------------------------------------------+
//| Recurrent subroutine for a InternalRealLinTranspose              |
//+------------------------------------------------------------------+
static void CFtBase::FFtIRLTRec(double &a[],const int astart,const int astride,
                                double &b[],const int bstart,const int bstride,
                                const int m,const int n)
  {
//--- create variables
   int i=0;
   int j=0;
   int idx1=0;
   int idx2=0;
   int m1=0;
   int n1=0;
//--- check
   if(m==0 || n==0)
      return;
//--- check
   if(MathMax(m,n)<=8)
     {
      for(i=0;i<=m-1;i++)
        {
         //--- calculation
         idx1=bstart+i;
         idx2=astart+i*astride;
         for(j=0;j<=n-1;j++)
           {
            b[idx1]=a[idx2];
            idx1=idx1+bstride;
            idx2=idx2+1;
           }
        }
      //--- exit the function
      return;
     }
//--- check
   if(n>m)
     {
      //--- New partition:
      //--- "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
      //---                                  ( B2 )
      n1=n/2;
      //--- check
      if(n-n1>=8 && n1%8!=0)
         n1=n1+(8-n1%8);
      //--- check
      if(!CAp::Assert(n-n1>0))
         return;
      //--- function call
      FFtIRLTRec(a,astart,astride,b,bstart,bstride,m,n1);
      //--- function call
      FFtIRLTRec(a,astart+n1,astride,b,bstart+n1*bstride,bstride,m,n-n1);
     }
   else
     {
      //--- New partition:
      //--- "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
      //---                     ( A2 )
      m1=m/2;
      //--- check
      if(m-m1>=8 && m1%8!=0)
         m1=m1+(8-m1%8);
      //--- check
      if(!CAp::Assert(m-m1>0))
         return;
      //--- function call
      FFtIRLTRec(a,astart,astride,b,bstart,bstride,m1,n);
      //--- function call
      FFtIRLTRec(a,astart+m1*astride,astride,b,bstart+m1,bstride,m-m1,n);
     }
  }
//+------------------------------------------------------------------+
//| recurrent subroutine for FFTFindSmoothRec                        |
//+------------------------------------------------------------------+
static void CFtBase::FtBaseFindSmoothRec(const int n,const int seed,
                                         const int leastfactor,int &best)
  {
//--- check
   if(!CAp::Assert(m_ftbasemaxsmoothfactor<=5,__FUNCTION__+": internal error!"))
      return;
//--- check
   if(seed>=n)
     {
      best=MathMin(best,seed);
      return;
     }
//--- check
   if(leastfactor<=2)
     {
      //--- function call
      FtBaseFindSmoothRec(n,seed*2,2,best);
     }
//--- check
   if(leastfactor<=3)
     {
      //--- function call
      FtBaseFindSmoothRec(n,seed*3,3,best);
     }
//--- check
   if(leastfactor<=5)
     {
      //--- function call
      FtBaseFindSmoothRec(n,seed*5,5,best);
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine: array resize                                |
//+------------------------------------------------------------------+
static void CFtBase::FFtArrayResize(int &a[],int &asize,const int newasize)
  {
//--- create variables
   int tmp[];
   int i=0;
//--- allocation
   ArrayResizeAL(tmp,asize);
   for(i=0;i<=asize-1;i++)
      tmp[i]=a[i];
//--- allocation
   ArrayResizeAL(a,newasize);
   for(i=0;i<=asize-1;i++)
      a[i]=tmp[i];
//--- get result
   asize=newasize;
  }
//+------------------------------------------------------------------+
//| Reference FHT stub                                               |
//+------------------------------------------------------------------+
static void CFtBase::ReFFHt(double &a[],const int n,const int offs)
  {
//--- create array
   double buf[];
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- allocation
   ArrayResizeAL(buf,n);
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      v=0;
      for(j=0;j<=n-1;j++)
         v=v+a[offs+j]*(MathCos(2*M_PI*i*j/n)+MathSin(2*M_PI*i*j/n));
      buf[i]=v;
     }
//--- get result
   for(i=0;i<=n-1;i++)
      a[offs+i]=buf[i];
  }
//+------------------------------------------------------------------+
//| Auxiliary class for calculation mathematical functions           |
//+------------------------------------------------------------------+
class CNearUnitYUnit
  {
public:
   //--- constructor, destructor
                     CNearUnitYUnit(void);
                    ~CNearUnitYUnit(void);
   //--- methods
   static double     NULog1p(const double x);
   static double     NUExp1m(const double x);
   static double     NUCos1m(const double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNearUnitYUnit::CNearUnitYUnit(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNearUnitYUnit::~CNearUnitYUnit(void)
  {

  }
//+------------------------------------------------------------------+
//| Log                                                              |
//+------------------------------------------------------------------+
static double CNearUnitYUnit::NULog1p(const double x)
  {
//--- create variables
   double z=1.0+x;
   double lp=0;
   double lq=0;
//--- check
   if(z<0.70710678118654752440 || z>1.41421356237309504880)
      return(MathLog(z));
//--- calculation result
   z=x*x;
   lp=4.5270000862445199635215E-5;
   lp=lp*x+4.9854102823193375972212E-1;
   lp=lp*x+6.5787325942061044846969E0;
   lp=lp*x+2.9911919328553073277375E1;
   lp=lp*x+6.0949667980987787057556E1;
   lp=lp*x+5.7112963590585538103336E1;
   lp=lp*x+2.0039553499201281259648E1;
   lq=1.0000000000000000000000E0;
   lq=lq*x+1.5062909083469192043167E1;
   lq=lq*x+8.3047565967967209469434E1;
   lq=lq*x+2.2176239823732856465394E2;
   lq=lq*x+3.0909872225312059774938E2;
   lq=lq*x+2.1642788614495947685003E2;
   lq=lq*x+6.0118660497603843919306E1;
   z=-(0.5*z)+x*(z*lp/lq);
//--- return result
   return(x+z);
  }
//+------------------------------------------------------------------+
//| Exp                                                              |
//+------------------------------------------------------------------+
static double CNearUnitYUnit::NUExp1m(const double x)
  {
//--- create variables
   double r;
   double xx;
   double ep;
   double eq;
//--- check
   if(x<-0.5 || x>0.5)
      return(MathExp(x)-1.0);
//--- calculation result
   xx=x*x;
   ep=1.2617719307481059087798E-4;
   ep=ep*xx+3.0299440770744196129956E-2;
   ep=ep*xx+9.9999999999999999991025E-1;
   eq=3.0019850513866445504159E-6;
   eq=eq*xx+2.5244834034968410419224E-3;
   eq=eq*xx+2.2726554820815502876593E-1;
   eq=eq*xx+2.0000000000000000000897E0;
   r=x*ep;
   r=r/(eq-r);
//--- return result
   return(r+r);
  }
//+------------------------------------------------------------------+
//| Cos                                                              |
//+------------------------------------------------------------------+
static double CNearUnitYUnit::NUCos1m(const double x)
  {
//--- create variables
   double xx;
   double c;
//--- check
   if(x<-0.25*M_PI || x>0.25*M_PI)
      return(MathCos(x)-1);
//--- get result
   xx=x*x;
   c=4.7377507964246204691685E-14;
   c=c*xx-1.1470284843425359765671E-11;
   c=c*xx+2.0876754287081521758361E-9;
   c=c*xx-2.7557319214999787979814E-7;
   c=c*xx+2.4801587301570552304991E-5;
   c=c*xx-1.3888888888888872993737E-3;
   c=c*xx+4.1666666666666666609054E-2;
//--- return result
   return(-(0.5*xx)+xx*xx*c);
  }
//+------------------------------------------------------------------+
