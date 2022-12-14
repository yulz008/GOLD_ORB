//+------------------------------------------------------------------+
//|                                                       linalg.mqh |
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
//| but WITHOUT ANY WARRANTY;without even the implied warranty of    |
//| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the     |
//| GNU General Public License for more details.                     |
//+------------------------------------------------------------------+
#include "alglibinternal.mqh"
#include "alglibmisc.mqh"
//+------------------------------------------------------------------+
//| Work with matrix forms                                           |
//+------------------------------------------------------------------+
class CAblas
  {
private:
   //--- split lenght
   static void       AblasInternalSplitLength(const int n,const int nb,int &n1,int &n2);
   //--- real numbers
   static void       RMatrixSyrk2(const int n,const int k,const double alpha,const CMatrixDouble &a,const int ia,const int ja,const int optypea,const double beta,CMatrixDouble &c,const int ic,const int jc,const bool isUpper);
   static void       RMatrixGemmK(const int m,const int n,const int k,const double alpha,const CMatrixDouble &a,const int ia,const int ja,const int optypea,const CMatrixDouble &b,const int ib,const int jb,const int optypeb,const double beta,CMatrixDouble &c,const int ic,const int jc);
   static void       RMatrixRightTrsM2(const int m,const int n,CMatrixDouble &a,const int i1,const int j1,const bool isUpper,const bool isUnit,const int optype,CMatrixDouble &x,const int i2,const int j2);
   static void       RMatrixLeftTrsM2(const int m,const int n,CMatrixDouble &a,const int i1,const int j1,const bool isUpper,const bool isUnit,const int optype,CMatrixDouble &x,const int i2,const int j2);
   //--- complex numbers
   static void       CMatrixSyrk2(const int n,const int k,const double alpha,const CMatrixComplex &a,const int ia,const int ja,const int optypea,const double beta,CMatrixComplex &c,const int ic,const int jc,const bool isUpper);
   static void       CMatrixGemmk(const int m,const int n,const int k,al_complex &alpha,const CMatrixComplex &a,const int ia,const int ja,const int optypea,const CMatrixComplex &b,const int ib,const int jb,const int optypeb,al_complex &beta,CMatrixComplex &c,const int ic,const int jc);
   static void       CMatrixRightTrsM2(const int m,const int n,CMatrixComplex &a,const int i1,const int j1,const bool isUpper,const bool isUnit,const int optype,CMatrixComplex &x,const int i2,const int j2);
   static void       CMatrixLeftTrsM2(const int m,const int n,CMatrixComplex &a,const int i1,const int j1,const bool isUpper,const bool isUnit,const int optype,CMatrixComplex &x,const int i2,const int j2);
public:
                     CAblas(void);
                    ~CAblas(void);
   //--- size
   static int        AblasBlockSize(void)        { return(32);}
   static int        AblasMicroBlockSize(void)   { return(8); }
   static int        AblasComplexBlockSize(void) { return(24);}
   //--- split lenght
   static void       AblasSplitLength(const CMatrixDouble &a,const int n,int &n1,int &n2);
   static void       AblasComplexSplitLength(const CMatrixComplex &a,const int n,int &n1,int &n2);
   //--- real numbers
   static void       RMatrixSyrk(const int n,const int k,const double alpha,const CMatrixDouble &a,const int ia,const int ja,const int optypea,const double beta,CMatrixDouble &c,const int ic,const int jc,const bool isUpper);
   static void       RMatrixGemm(const int m,const int n,const int k,const double alpha,const CMatrixDouble &a,const int ia,const int ja,const int optypea,const CMatrixDouble &b,const int ib,const int jb,const int optypeb,const double beta,CMatrixDouble &c,const int ic,const int jc);
   static void       RMatrixTranspose(const int m,const int n,const CMatrixDouble &a,const int ia,const int ja,CMatrixDouble &b,const int ib,const int jb);
   static void       RMatrixCopy(const int m,const int n,const CMatrixDouble &a,const int ia,const int ja,CMatrixDouble &b,const int ib,const int jb);
   static void       RMatrixRank1(const int m,const int n,CMatrixDouble &a,const int ia,const int ja,const double &u[],const int iu,const double &v[],const int iv);
   static void       RMatrixMVect(const int m,const int n,const CMatrixDouble &a,const int ia,const int ja,const int opa,const double &x[],const int ix,double &y[],const int iy);
   static void       RMatrixRightTrsM(const int m,const int n,CMatrixDouble &a,const int i1,const int j1,const bool isUpper,const bool isUnit,const int optype,CMatrixDouble &x,const int i2,const int j2);
   static void       RMatrixLeftTrsM(const int m,const int n,CMatrixDouble &a,const int i1,const int j1,const bool isUpper,const bool isUnit,const int optype,CMatrixDouble &x,const int i2,const int j2);
   //--- complex numbers
   static void       CMatrixSyrk(const int n,const int k,const double alpha,CMatrixComplex &a,const int ia,const int ja,const int optypea,const double beta,CMatrixComplex &c,const int ic,const int jc,const bool isUpper);
   static void       CMatrixGemm(const int m,const int n,const int k,al_complex &alpha,CMatrixComplex &a,const int ia,const int ja,const int optypea,CMatrixComplex &b,const int ib,const int jb,const int optypeb,al_complex &beta,CMatrixComplex &c,const int ic,const int jc);
   static void       CMatrixTranspose(const int m,const int n,const CMatrixComplex &a,const int ia,const int ja,CMatrixComplex &b,const int ib,const int jb);
   static void       CMatrixCopy(const int m,const int n,const CMatrixComplex &a,const int ia,const int ja,CMatrixComplex &b,const int ib,const int jb);
   static void       CMatrixRank1(const int m,const int n,CMatrixComplex &a,const int ia,const int ja,const al_complex &u[],const int iu,const al_complex &v[],const int iv);
   static void       CMatrixMVect(const int m,const int n,const CMatrixComplex &a,const int ia,const int ja,const int opa,const al_complex &x[],const int ix,al_complex &y[],const int iy);
   static void       CMatrixRightTrsM(const int m,const int n,CMatrixComplex &a,const int i1,const int j1,const bool isUpper,const bool isUnit,const int optype,CMatrixComplex &x,const int i2,const int j2);
   static void       CMatrixLeftTrsM(const int m,const int n,CMatrixComplex &a,const int i1,const int j1,const bool isUpper,const bool isUnit,const int optype,CMatrixComplex &x,const int i2,const int j2);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CAblas::CAblas(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CAblas::~CAblas(void)
  {

  }
//+------------------------------------------------------------------+
//| Same as CMatrixSYRK, but for real matrices                       |
//| OpType may be only 0 or 1.                                       |
//+------------------------------------------------------------------+
static void CAblas::RMatrixSyrk(const int n,const int k,const double alpha,
                                const CMatrixDouble &a,const int ia,const int ja,
                                const int optypea,const double beta,CMatrixDouble &c,
                                const int ic,const int jc,const bool isUpper)
  {
//--- create variables
   int s1=0;
   int s2=0;
   int bs=0;
//--- calculation size
   bs=AblasBlockSize();
//--- check
   if(n<=bs && k<=bs)
     {
      RMatrixSyrk2(n,k,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
      return;
     }
//--- check
   if(k>=n)
     {
      //--- Split K
      AblasSplitLength(a,k,s1,s2);
      //--- check
      if(optypea==0)
        {
         RMatrixSyrk(n,s1,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
         RMatrixSyrk(n,s2,alpha,a,ia,ja+s1,optypea,1.0,c,ic,jc,isUpper);
        }
      else
        {
         RMatrixSyrk(n,s1,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
         RMatrixSyrk(n,s2,alpha,a,ia+s1,ja,optypea,1.0,c,ic,jc,isUpper);
        }
     }
   else
     {
      //--- Split N
      AblasSplitLength(a,n,s1,s2);
      //--- check
      if(optypea==0 && isUpper)
        {
         RMatrixSyrk(s1,k,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
         RMatrixGemm(s1,s2,k,alpha,a,ia,ja,0,a,ia+s1,ja,1,beta,c,ic,jc+s1);
         RMatrixSyrk(s2,k,alpha,a,ia+s1,ja,optypea,beta,c,ic+s1,jc+s1,isUpper);
         //--- exit the function
         return;
        }
      //--- check
      if(optypea==0 && !isUpper)
        {
         RMatrixSyrk(s1,k,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
         RMatrixGemm(s2,s1,k,alpha,a,ia+s1,ja,0,a,ia,ja,1,beta,c,ic+s1,jc);
         RMatrixSyrk(s2,k,alpha,a,ia+s1,ja,optypea,beta,c,ic+s1,jc+s1,isUpper);
         //--- exit the function
         return;
        }
      //--- check
      if(optypea!=0 && isUpper)
        {
         RMatrixSyrk(s1,k,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
         RMatrixGemm(s1,s2,k,alpha,a,ia,ja,1,a,ia,ja+s1,0,beta,c,ic,jc+s1);
         RMatrixSyrk(s2,k,alpha,a,ia,ja+s1,optypea,beta,c,ic+s1,jc+s1,isUpper);
         //--- exit the function
         return;
        }
      //--- check
      if(optypea!=0 && !isUpper)
        {
         RMatrixSyrk(s1,k,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
         RMatrixGemm(s2,s1,k,alpha,a,ia,ja+s1,1,a,ia,ja,0,beta,c,ic+s1,jc);
         RMatrixSyrk(s2,k,alpha,a,ia,ja+s1,optypea,beta,c,ic+s1,jc+s1,isUpper);
         //--- exit the function
         return;
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Splits matrix length in two parts, left part should match ABLAS  |
//| block size                                                       |
//| INPUT PARAMETERS                                                 |
//|     A   -   real matrix, is passed to ensure that we didn't split|
//|             complex matrix using real splitting subroutine.      |
//|             matrix itself is not changed.                        |
//|     N   -   length, N>0                                          |
//| OUTPUT PARAMETERS                                                |
//|     N1  -   length                                               |
//|     N2  -   length                                               |
//| N1+N2=N, N1>=N2, N2 may be zero                                  |
//+------------------------------------------------------------------+
static void CAblas::AblasSplitLength(const CMatrixDouble &a,const int n,
                                     int &n1,int &n2)
  {
//--- initialization
   n1=0;
   n2=0;
//--- check
   if(n>AblasBlockSize())
      AblasInternalSplitLength(n,AblasBlockSize(),n1,n2);
   else
      AblasInternalSplitLength(n,AblasMicroBlockSize(),n1,n2);
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Complex ABLASSplitLength                                         |
//+------------------------------------------------------------------+
static void CAblas::AblasComplexSplitLength(const CMatrixComplex &a,const int n,
                                            int &n1,int &n2)
  {
//--- check
   if(n>AblasComplexBlockSize())
      AblasInternalSplitLength(n,AblasComplexBlockSize(),n1,n2);
   else
      AblasInternalSplitLength(n,AblasMicroBlockSize(),n1,n2);
  }
//+------------------------------------------------------------------+
//| Same as CMatrixGEMM, but for real numbers.                       |
//| OpType may be only 0 or 1.                                       |
//+------------------------------------------------------------------+
static void CAblas::RMatrixGemm(const int m,const int n,const int k,const double alpha,
                                const CMatrixDouble &a,const int ia,const int ja,
                                const int optypea,const CMatrixDouble &b,
                                const int ib,const int jb,const int optypeb,
                                const double beta,CMatrixDouble &c,const int ic,
                                const int jc)
  {
//--- create variables
   int s1=0;
   int s2=0;
   int bs;
//--- calculation size
   bs=AblasBlockSize();
//--- check
   if(m<=bs && n<=bs && k<=bs)
     {
      RMatrixGemmK(m,n,k,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
      return;
     }
//--- check
   if(m>=n && m>=k)
     {
      //--- A*B = (A1 A2)^T*B
      AblasSplitLength(a,m,s1,s2);
      //--- check
      if(optypea==0)
        {
         RMatrixGemm(s1,n,k,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         RMatrixGemm(s2,n,k,alpha,a,ia+s1,ja,optypea,b,ib,jb,optypeb,beta,c,ic+s1,jc);
        }
      else
        {
         RMatrixGemm(s1,n,k,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         RMatrixGemm(s2,n,k,alpha,a,ia,ja+s1,optypea,b,ib,jb,optypeb,beta,c,ic+s1,jc);
        }
      return;
     }
//--- check
   if(n>=m && n>=k)
     {
      //--- A*B = A*(B1 B2)
      AblasSplitLength(a,n,s1,s2);
      //--- check
      if(optypeb==0)
        {
         RMatrixGemm(m,s1,k,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         RMatrixGemm(m,s2,k,alpha,a,ia,ja,optypea,b,ib,jb+s1,optypeb,beta,c,ic,jc+s1);
        }
      else
        {
         RMatrixGemm(m,s1,k,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         RMatrixGemm(m,s2,k,alpha,a,ia,ja,optypea,b,ib+s1,jb,optypeb,beta,c,ic,jc+s1);
        }
      return;
     }
//--- check
   if(k>=m && k>=n)
     {
      //--- A*B = (A1 A2)*(B1 B2)^T
      AblasSplitLength(a,k,s1,s2);
      //--- check
      if(optypea==0 && optypeb==0)
        {
         RMatrixGemm(m,n,s1,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         RMatrixGemm(m,n,s2,alpha,a,ia,ja+s1,optypea,b,ib+s1,jb,optypeb,1.0,c,ic,jc);
        }
      //--- check
      if(optypea==0 && optypeb!=0)
        {
         RMatrixGemm(m,n,s1,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         RMatrixGemm(m,n,s2,alpha,a,ia,ja+s1,optypea,b,ib,jb+s1,optypeb,1.0,c,ic,jc);
        }
      //--- check
      if(optypea!=0 && optypeb==0)
        {
         RMatrixGemm(m,n,s1,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         RMatrixGemm(m,n,s2,alpha,a,ia+s1,ja,optypea,b,ib+s1,jb,optypeb,1.0,c,ic,jc);
        }
      //--- check
      if(optypea!=0 && optypeb!=0)
        {
         RMatrixGemm(m,n,s1,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         RMatrixGemm(m,n,s2,alpha,a,ia+s1,ja,optypea,b,ib,jb+s1,optypeb,1.0,c,ic,jc);
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Level 2 subrotuine                                               |
//+------------------------------------------------------------------+
static void CAblas::RMatrixSyrk2(const int n,const int k,const double alpha,
                                 const CMatrixDouble &a,const int ia,const int ja,
                                 const int optypea,const double beta,CMatrixDouble &c,
                                 const int ic,const int jc,const bool isUpper)
  {
//--- check
   if((alpha==0.0 || k==0.0) && beta==1.0)
      return;
//--- create variables
   int    i=0;
   int    j=0;
   int    j1=0;
   int    j2=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(optypea==0)
     {
      //--- C=alpha*A*A^H+beta*C
      for(i=0;i<n;i++)
        {
         if(isUpper)
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
            if(alpha!=0 && k>0)
              {
               v=0.0;
               for(i_=ja;i_<ja+k;i_++)
                  v+=a[ia+i][i_]*a[ia+j][i_];
              }
            else
               v=0;
            //--- check
            if(beta==0)
               c[ic+i].Set(jc+j,alpha*v);
            else
               c[ic+i].Set(jc+j,beta*c[ic+i][jc+j]+alpha*v);
           }
        }
      return;
     }
   else
     {
      //--- C=alpha*A^H*A+beta*C
      for(i=0;i<n;i++)
        {
         //--- check
         if(isUpper)
           {
            j1=i;
            j2=n-1;
           }
         else
           {
            j1=0;
            j2=i;
           }
         //--- check
         if(beta==0)
           {
            for(j=j1;j<=j2;j++)
               c[ic+i].Set(jc+j,0);
           }
         else
           {
            for(i_=jc+j1;i_<=jc+j2;i_++)
               c[ic+i].Set(i_,beta*c[ic+i][i_]);
           }
        }
      for(i=0;i<k;i++)
        {
         for(j=0;j<n;j++)
           {
            //--- check
            if(isUpper)
              {
               j1=j;
               j2=n-1;
              }
            else
              {
               j1=0;
               j2=j;
              }
            //--- change values
            v=alpha*a[ia+i][ja+j];
            i1_=(ja+j1)-(jc+j1);
            for(i_=jc+j1;i_<=jc+j2;i_++)
               c[ic+j].Set(i_,c[ic+j][i_]+v*a[ia+i][i_+i1_]);
           }
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Complex ABLASSplitLength                                         |
//+------------------------------------------------------------------+
static void CAblas::AblasInternalSplitLength(const int n,const int nb,
                                             int &n1,int &n2)
  {
//--- initialization
   int r=0;
   n1=0;
   n2=0;
//--- check
   if(n<=nb)
     {
      //--- Block size, no further splitting
      n1=n;
      n2=0;
     }
   else
     {
      //--- Greater than block size
      if(n%nb!=0)
        {
         //--- Split remainder
         n2=n%nb;
         n1=n-n2;
        }
      else
        {
         //--- Split on block boundaries
         n2=n/2;
         n1=n-n2;
         //--- check
         if(n1%nb==0)
            return;
         r=nb-n1%nb;
         n1+=r;
         n2-=r;
        }
     }
  }
//+------------------------------------------------------------------+
//| GEMM kernel                                                      |
//+------------------------------------------------------------------+
static void CAblas::RMatrixGemmK(const int m,const int n,const int k,const double alpha,
                                 const CMatrixDouble &a,const int ia,const int ja,
                                 const int optypea,const CMatrixDouble &b,const int ib,
                                 const int jb,const int optypeb,const double beta,
                                 CMatrixDouble &c,const int ic,const int jc)
  {
//--- check
   if(m*n==0)
      return;
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- if K=0, then C=Beta*C
   if(k==0)
     {
      //--- check
      if(beta!=1)
        {
         //--- check
         if(beta!=0)
           {
            for(i=0;i<m;i++)
               for(j=0;j<n;j++)
                  c[ic+i].Set(jc+j,beta*c[ic+i][jc+j]);
           }
         else
           {
            for(i=0;i<m;i++)
               for(j=0;j<n;j++)
                  c[ic+i].Set(jc+j,0);
           }
        }
      return;
     }
//--- check
   if(optypea==0 && optypeb!=0)
     {
      //--- a*b'
      for(i=0;i<m;i++)
        {
         for(j=0;j<n;j++)
           {
            //--- check
            if(k==0 || alpha==0)
               v=0;
            else
              {
               i1_=jb-ja;
               v=0.0;
               for(i_=ja;i_<ja+k;i_++)
                  v+=a[ia+i][i_]*b[ib+j][i_+i1_];
              }
            //--- check
            if(beta==0)
               c[ic+i].Set(jc+j,alpha*v);
            else
               c[ic+i].Set(jc+j,beta*c[ic+i][jc+j]+alpha*v);
           }
        }
      return;
     }
//--- check
   if(optypea==0 && optypeb==0)
     {
      //--- a*b
      for(i=0;i<m;i++)
        {
         //--- check
         if(beta!=0)
           {
            for(i_=jc;i_<jc+n;i_++)
               c[ic+i].Set(i_,beta*c[ic+i][i_]);
           }
         else
           {
            for(j=0;j<n;j++)
               c[ic+i].Set(jc+j,0);
           }
         //--- check
         if(alpha!=0)
           {
            for(j=0;j<k;j++)
              {
               v=alpha*a[ia+i][ja+j];
               i1_=jb-jc;
               for(i_=jc;i_<jc+n;i_++)
                  c[ic+i].Set(i_,c[ic+i][i_]+v*b[ib+j][i_+i1_]);
              }
           }
        }
      return;
     }
//--- check
   if(optypea!=0 && optypeb!=0)
     {
      //--- a`*b`
      for(i=0;i<m;i++)
        {
         for(j=0;j<n;j++)
           {
            //--- check
            if(alpha==0)
               v=0;
            else
              {
               i1_=(jb)-(ia);
               v=0.0;
               for(i_=ia;i_<ia+k;i_++)
                  v+=a[i_][ja+i]*b[ib+j][i_+i1_];
              }
            //--- check
            if(beta==0)
               c[ic+i].Set(jc+j,alpha*v);
            else
               c[ic+i].Set(jc+j,beta*c[ic+i][jc+j]+alpha*v);
           }
        }
      return;
     }
//--- check
   if(optypea!=0 && optypeb==0)
     {
      //--- a`*b
      if(beta==0)
        {
         for(i=0;i<m;i++)
            for(j=0;j<n;j++)
               c[ic+i].Set(jc+j,0);
        }
      else
        {
         for(i=0;i<m;i++)
            for(i_=jc;i_<jc+n;i_++)
               c[ic+i].Set(i_,beta*c[ic+i][i_]);
        }
      //--- check
      if(alpha!=0)
        {
         for(j=0;j<k;j++)
            for(i=0;i<m;i++)
              {
               v=alpha*a[ia+j][ja+i];
               i1_=jb-jc;
               for(i_=jc;i_<jc+n;i_++)
                  c[ic+i].Set(i_,c[ic+i][i_]+v*b[ib+j][i_+i1_]);
              }
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Cache-oblivous complex "copy-and-transpose"                      |
//| Input parameters:                                                |
//|     M   -   number of rows                                       |
//|     N   -   number of columns                                    |
//|     A   -   source matrix, MxN submatrix is copied and transposed|
//|     IA  -   submatrix offset (row index)                         |
//|     JA  -   submatrix offset (column index)                      |
//|     A   -   destination matrix                                   |
//|     IB  -   submatrix offset (row index)                         |
//|     JB  -   submatrix offset (column index)                      |
//+------------------------------------------------------------------+
static void CAblas::CMatrixTranspose(const int m,const int n,const CMatrixComplex &a,
                                     const int ia,const int ja,CMatrixComplex &b,
                                     const int ib,const int jb)
  {
//--- create variables
   int i=0;
   int s1=0;
   int s2=0;
   int i_=0;
   int i1_=0;
//--- check
   if(m<=2*AblasComplexBlockSize() && n<=2*AblasComplexBlockSize())
     {
      //--- base case
      for(i=0;i<m;i++)
        {
         i1_=ja-ib;
         for(i_=ib;i_<ib+n;i_++)
            b[i_].Set(jb+i,a[ia+i][i_+i1_]);
        }
     }
   else
     {
      //--- Cache-oblivious recursion
      if(m>n)
        {
         //--- split
         AblasComplexSplitLength(a,m,s1,s2);
         //--- function call
         CMatrixTranspose(s1,n,a,ia,ja,b,ib,jb);
         CMatrixTranspose(s2,n,a,ia+s1,ja,b,ib,jb+s1);
        }
      else
        {
         //--- split
         AblasComplexSplitLength(a,n,s1,s2);
         //--- function call
         CMatrixTranspose(m,s1,a,ia,ja,b,ib,jb);
         CMatrixTranspose(m,s2,a,ia,ja+s1,b,ib+s1,jb);
        }
     }
  }
//+------------------------------------------------------------------+
//| Cache-oblivous real "copy-and-transpose"                         |
//| Input parameters:                                                |
//|     M   -   number of rows                                       |
//|     N   -   number of columns                                    |
//|     A   -   source matrix, MxN submatrix is copied and transposed|
//|     IA  -   submatrix offset (row index)                         |
//|     JA  -   submatrix offset (column index)                      |
//|     A   -   destination matrix                                   |
//|     IB  -   submatrix offset (row index)                         |
//|     JB  -   submatrix offset (column index)                      |
//+------------------------------------------------------------------+
static void CAblas::RMatrixTranspose(const int m,const int n,const CMatrixDouble &a,
                                     const int ia,const int ja,CMatrixDouble &b,
                                     const int ib,const int jb)
  {
//--- create variables
   int i=0;
   int s1=0;
   int s2=0;
   int i_=0;
   int i1_=0;
//--- check
   if(m<=2*AblasBlockSize() && n<=2*AblasBlockSize())
     {
      //--- base case
      for(i=0;i<m;i++)
        {
         i1_=ja-ib;
         for(i_=ib;i_<ib+n;i_++)
            b[i_].Set(jb+i,a[ia+i][i_+i1_]);
        }
     }
   else
     {
      //--- Cache-oblivious recursion
      if(m>n)
        {
         //--- split
         AblasSplitLength(a,m,s1,s2);
         //--- function call
         RMatrixTranspose(s1,n,a,ia,ja,b,ib,jb);
         RMatrixTranspose(s2,n,a,ia+s1,ja,b,ib,jb+s1);
        }
      else
        {
         //--- split
         AblasSplitLength(a,n,s1,s2);
         //--- function call
         RMatrixTranspose(m,s1,a,ia,ja,b,ib,jb);
         RMatrixTranspose(m,s2,a,ia,ja+s1,b,ib+s1,jb);
        }
     }
  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//| Input parameters:                                                |
//|     M   -   number of rows                                       |
//|     N   -   number of columns                                    |
//|     A   -   source matrix, MxN submatrix is copied and transposed|
//|     IA  -   submatrix offset (row index)                         |
//|     JA  -   submatrix offset (column index)                      |
//|     B   -   destination matrix                                   |
//|     IB  -   submatrix offset (row index)                         |
//|     JB  -   submatrix offset (column index)                      |
//+------------------------------------------------------------------+
static void CAblas::CMatrixCopy(const int m,const int n,const CMatrixComplex &a,
                                const int ia,const int ja,CMatrixComplex &b,
                                const int ib,const int jb)
  {
//--- create variables
   int i=0;
   int i_=0;
   int i1_=0;
//--- copy
   for(i=0;i<m;i++)
     {
      i1_=ja-jb;
      for(i_=jb;i_<jb+n;i_++)
         b[ib+i].Set(i_,a[ia+i][i_+i1_]);
     }
  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//| Input parameters:                                                |
//|     M   -   number of rows                                       |
//|     N   -   number of columns                                    |
//|     A   -   source matrix, MxN submatrix is copied and transposed|
//|     IA  -   submatrix offset (row index)                         |
//|     JA  -   submatrix offset (column index)                      |
//|     B   -   destination matrix                                   |
//|     IB  -   submatrix offset (row index)                         |
//|     JB  -   submatrix offset (column index)                      |
//+------------------------------------------------------------------+
static void CAblas::RMatrixCopy(const int m,const int n,const CMatrixDouble &a,
                                const int ia,const int ja,CMatrixDouble &b,
                                const int ib,const int jb)
  {
//--- create variables
   int i=0;
   int i_=0;
   int i1_=0;
//--- copy
   for(i=0;i<m;i++)
     {
      i1_=ja-jb;
      for(i_=jb;i_<jb+n;i_++)
         b[ib+i].Set(i_,a[ia+i][i_+i1_]);
     }
  }
//+------------------------------------------------------------------+
//| Rank-1 correction: A := A + u*v'                                 |
//| INPUT PARAMETERS:                                                |
//|     M   -   number of rows                                       |
//|     N   -   number of columns                                    |
//|     A   -   target matrix, MxN submatrix is updated              |
//|     IA  -   submatrix offset (row index)                         |
//|     JA  -   submatrix offset (column index)                      |
//|     U   -   vector #1                                            |
//|     IU  -   subvector offset                                     |
//|     V   -   vector #2                                            |
//|     IV  -   subvector offset                                     |
//+------------------------------------------------------------------+
static void CAblas::CMatrixRank1(const int m,const int n,CMatrixComplex &a,
                                 const int ia,const int ja,const al_complex &u[],
                                 const int iu,const al_complex &v[],const int iv)
  {
//--- check
   if(m==0 || n==0)
      return;
//--- create variables
   int        i=0;
   al_complex s=0;
   int        i_=0;
   int        i1_=0;
//--- correction
   for(i=0;i<m;i++)
     {
      s=u[iu+i];
      i1_=iv-ja;
      for(i_=ja;i_<ja+n;i_++)
         a[ia+i].Set(i_,a[ia+i][i_]+s*v[i_+i1_]);
     }
  }
//+------------------------------------------------------------------+
//| Rank-1 correction: A := A + u*v'                                 |
//| INPUT PARAMETERS:                                                |
//|     M   -   number of rows                                       |
//|     N   -   number of columns                                    |
//|     A   -   target matrix, MxN submatrix is updated              |
//|     IA  -   submatrix offset (row index)                         |
//|     JA  -   submatrix offset (column index)                      |
//|     U   -   vector #1                                            |
//|     IU  -   subvector offset                                     |
//|     V   -   vector #2                                            |
//|     IV  -   subvector offset                                     |
//+------------------------------------------------------------------+
static void CAblas::RMatrixRank1(const int m,const int n,CMatrixDouble &a,
                                 const int ia,const int ja,const double &u[],
                                 const int iu,const double &v[],const int iv)
  {
//--- check
   if(m==0 || n==0)
      return;
//--- create variables
   int    i=0;
   double s=0;
   int    i_=0;
   int    i1_=0;
//--- correction
   for(i=0;i<m;i++)
     {
      s=u[iu+i];
      i1_=iv-ja;
      for(i_=ja;i_<ja+n;i_++)
         a[ia+i].Set(i_,a[ia+i][i_]+s*v[i_+i1_]);
     }
  }
//+------------------------------------------------------------------+
//| Matrix-vector product: y := op(A)*x                              |
//| INPUT PARAMETERS:                                                |
//|     M   -   number of rows of op(A)                              |
//|             M>=0                                                 |
//|     N   -   number of columns of op(A)                           |
//|             N>=0                                                 |
//|     A   -   target matrix                                        |
//|     IA  -   submatrix offset (row index)                         |
//|     JA  -   submatrix offset (column index)                      |
//|     OpA -   operation type:                                      |
//|             * OpA=0     =>  op(A) = A                            |
//|             * OpA=1     =>  op(A) = A^T                          |
//|             * OpA=2     =>  op(A) = A^H                          |
//|     X   -   input vector                                         |
//|     IX  -   subvector offset                                     |
//|     IY  -   subvector offset                                     |
//| OUTPUT PARAMETERS:                                               |
//|     Y   -   vector which stores result                           |
//| if M=0, then subroutine does nothing.                            |
//| if N=0, Y is filled by zeros.                                    |
//+------------------------------------------------------------------+
static void CAblas::CMatrixMVect(const int m,const int n,const CMatrixComplex &a,
                                 const int ia,const int ja,const int opa,
                                 const al_complex &x[],const int ix,al_complex &y[],
                                 const int iy)
  {
//--- create variables
   int        i=0;
   al_complex v=0;
   int        i_=0;
   int        i1_=0;
//--- check
   if(m==0)
      return;
//--- check
   if(n==0)
     {
      for(i=0;i<m;i++)
         y[iy+i]=0;
      //--- exit the function
      return;
     }
//--- check
   if(opa==0)
     {
      //--- y = A*x
      for(i=0;i<m;i++)
        {
         i1_=ix-ja;
         v=0.0;
         for(i_=ja;i_<ja+n;i_++)
            v+=a[ia+i][i_]*x[i_+i1_];
         //--- get y
         y[iy+i]=v;
        }
      //--- exit the function
      return;
     }
//--- check
   if(opa==1)
     {
      //--- y = A^T*x
      for(i=0;i<m;i++)
         y[iy+i]=0;
      for(i=0;i<n;i++)
        {
         v=x[ix+i];
         i1_=ja-iy;
         for(i_=iy;i_<=iy+m-1;i_++)
            y[i_]+=v*a[ia+i][i_+i1_];
        }
      //--- exit the function
      return;
     }
//--- check
   if(opa==2)
     {
      //--- y = A^H*x
      for(i=0;i<m;i++)
         y[iy+i]=0;
      for(i=0;i<n;i++)
        {
         v=x[ix+i];
         i1_=ja-iy;
         for(i_=iy;i_<=iy+m-1;i_++)
            y[i_]+=v*CMath::Conj(a[ia+i][i_+i1_]);
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Matrix-vector product: y := op(A)*x                              |
//| INPUT PARAMETERS:                                                |
//|     M   -   number of rows of op(A)                              |
//|     N   -   number of columns of op(A)                           |
//|     A   -   target matrix                                        |
//|     IA  -   submatrix offset (row index)                         |
//|     JA  -   submatrix offset (column index)                      |
//|     OpA -   operation type:                                      |
//|             * OpA=0     =>  op(A) = A                            |
//|             * OpA=1     =>  op(A) = A^T                          |
//|     X   -   input vector                                         |
//|     IX  -   subvector offset                                     |
//|     IY  -   subvector offset                                     |
//| OUTPUT PARAMETERS:                                               |
//|     Y   -   vector which stores result                           |
//| if M=0, then subroutine does nothing.                            |
//| if N=0, Y is filled by zeros.                                    |
//+------------------------------------------------------------------+
static void CAblas::RMatrixMVect(const int m,const int n,const CMatrixDouble &a,
                                 const int ia,const int ja,const int opa,
                                 const double &x[],const int ix,double &y[],
                                 const int iy)
  {
//--- create variables
   int    i=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(m==0)
      return;
//--- check
   if(n==0)
     {
      for(i=0;i<m;i++)
         y[iy+i]=0;
      //--- exit the function
      return;
     }
//--- check
   if(opa==0)
     {
      //--- y = A*x
      for(i=0;i<m;i++)
        {
         i1_=ix-ja;
         v=0.0;
         for(i_=ja;i_<ja+n;i_++)
            v+=a[ia+i][i_]*x[i_+i1_];
         y[iy+i]=v;
        }
      //--- exit the function
      return;
     }
//--- check
   if(opa==1)
     {
      //--- y = A^T*x
      for(i=0;i<m;i++)
         y[iy+i]=0;
      for(i=0;i<n;i++)
        {
         v=x[ix+i];
         i1_=ja-iy;
         for(i_=iy;i_<iy+m;i_++)
            y[i_]+=v*a[ia+i][i_+i1_];
        }
      //--- exit the function
      return;
     }
  }
//+------------------------------------------------------------------+
//| This subroutine calculates X*op(A^-1) where:                     |
//| * X is MxN general matrix                                        |
//| * A is NxN upper/lower triangular/unitriangular matrix           |
//| * "op" may be identity transformation, transposition, conjugate  |
//|   transposition                                                  |
//| Multiplication result replaces X.                                |
//| Cache-oblivious algorithm is used.                               |
//| INPUT PARAMETERS                                                 |
//|     N   -   matrix size, N>=0                                    |
//|     M   -   matrix size, N>=0                                    |
//|     A       -   matrix, actial matrix is stored in               |
//|                 A[I1:I1+N-1,J1:J1+N-1]                           |
//|     I1      -   submatrix offset                                 |
//|     J1      -   submatrix offset                                 |
//|     IsUpper -   whether matrix is upper triangular               |
//|     IsUnit  -   whether matrix is unitriangular                  |
//|     OpType  -   transformation type:                             |
//|                 * 0 - no transformation                          |
//|                 * 1 - transposition                              |
//|                 * 2 - conjugate transposition                    |
//|     C   -   matrix, actial matrix is stored in                   |
//|             C[I2:I2+M-1,J2:J2+N-1]                               |
//|     I2  -   submatrix offset                                     |
//|     J2  -   submatrix offset                                     |
//+------------------------------------------------------------------+
static void CAblas::CMatrixRightTrsM(const int m,const int n,CMatrixComplex &a,
                                     const int i1,const int j1,const bool isUpper,
                                     const bool isUnit,const int optype,
                                     CMatrixComplex &x,const int i2,const int j2)
  {
//--- create variables
   al_complex Alpha(-1,0);
   al_complex Beta(1,0);
   int        s1=0;
   int        s2=0;
   int        bs=AblasComplexBlockSize();
//--- check
   if(m<=bs && n<=bs)
     {
      //--- basic algorithm
      CMatrixRightTrsM2(m,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
      //--- exit the function
      return;
     }
//--- check
   if(m>=n)
     {
      //--- Split X: X*A = (X1 X2)^T*A
      AblasComplexSplitLength(a,m,s1,s2);
      CMatrixRightTrsM(s1,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
      CMatrixRightTrsM(s2,n,a,i1,j1,isUpper,isUnit,optype,x,i2+s1,j2);
     }
   else
     {
      //--- Split A:
      //---               (A1  A12)
      //--- X*op(A) = X*op(       )
      //---               (     A2)
      //---
      //--- Different variants depending on
      //--- IsUpper/OpType combinations
      AblasComplexSplitLength(a,n,s1,s2);
      //--- check
      if(isUpper && optype==0)
        {
         //---                  (A1  A12)-1
         //--- X*A^-1 = (X1 X2)*(       )
         //---                  (     A2)
         CMatrixRightTrsM(m,s1,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         CMatrixGemm(m,s2,s1,Alpha,x,i2,j2,0,a,i1,j1+s1,0,Beta,x,i2,j2+s1);
         CMatrixRightTrsM(m,s2,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2,j2+s1);
         //--- exit the function
         return;
        }
      //--- check
      if(isUpper && optype!=0)
        {
         //---                  (A1'     )-1
         //--- X*A^-1 = (X1 X2)*(        )
         //---                  (A12' A2')
         CMatrixRightTrsM(m,s2,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2,j2+s1);
         CMatrixGemm(m,s1,s2,Alpha,x,i2,j2+s1,0,a,i1,j1+s1,optype,Beta,x,i2,j2);
         CMatrixRightTrsM(m,s1,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         //--- exit the function
         return;
        }
      //--- check
      if(!isUpper && optype==0)
        {
         //---                  (A1     )-1
         //--- X*A^-1 = (X1 X2)*(       )
         //---                  (A21  A2)
         CMatrixRightTrsM(m,s2,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2,j2+s1);
         CMatrixGemm(m,s1,s2,Alpha,x,i2,j2+s1,0,a,i1+s1,j1,0,Beta,x,i2,j2);
         CMatrixRightTrsM(m,s1,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         //--- exit the function
         return;
        }
      //--- check
      if(!isUpper && optype!=0)
        {
         //---                  (A1' A21')-1
         //--- X*A^-1 = (X1 X2)*(        )
         //---                  (     A2')
         CMatrixRightTrsM(m,s1,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         CMatrixGemm(m,s2,s1,Alpha,x,i2,j2,0,a,i1+s1,j1,optype,Beta,x,i2,j2+s1);
         CMatrixRightTrsM(m,s2,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2,j2+s1);
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| This subroutine calculates op(A^-1)*X where:                     |
//| * X is MxN general matrix                                        |
//| * A is MxM upper/lower triangular/unitriangular matrix           |
//| * "op" may be identity transformation, transposition, conjugate  |
//|   transposition                                                  |
//| Multiplication result replaces X.                                |
//| Cache-oblivious algorithm is used.                               |
//| INPUT PARAMETERS                                                 |
//|     N   -   matrix size, N>=0                                    |
//|     M   -   matrix size, N>=0                                    |
//|     A       -   matrix, actial matrix is stored in               |
//|                 A[I1:I1+M-1,J1:J1+M-1]                           |
//|     I1      -   submatrix offset                                 |
//|     J1      -   submatrix offset                                 |
//|     IsUpper -   whether matrix is upper triangular               |
//|     IsUnit  -   whether matrix is unitriangular                  |
//|     OpType  -   transformation type:                             |
//|                 * 0 - no transformation                          |
//|                 * 1 - transposition                              |
//|                 * 2 - conjugate transposition                    |
//|     C   -   matrix, actial matrix is stored in                   |
//|             C[I2:I2+M-1,J2:J2+N-1]                               |
//|     I2  -   submatrix offset                                     |
//|     J2  -   submatrix offset                                     |
//+------------------------------------------------------------------+
static void CAblas::CMatrixLeftTrsM(const int m,const int n,CMatrixComplex &a,
                                    const int i1,const int j1,const bool isUpper,
                                    const bool isUnit,const int optype,
                                    CMatrixComplex &x,const int i2,const int j2)
  {
//--- create variables
   al_complex Alpha(-1,0);
   al_complex Beta(1,0);
   int        s1=0;
   int        s2=0;
   int        bs=AblasComplexBlockSize();
//--- check
   if(m<=bs && n<=bs)
     {
      //--- basic algorithm
      CMatrixLeftTrsM2(m,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
      //--- exit the function
      return;
     }
//--- check
   if(n>=m)
     {
      //--- Split X: op(A)^-1*X = op(A)^-1*(X1 X2)
      AblasComplexSplitLength(x,n,s1,s2);
      CMatrixLeftTrsM(m,s1,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
      CMatrixLeftTrsM(m,s2,a,i1,j1,isUpper,isUnit,optype,x,i2,j2+s1);
     }
   else
     {
      //--- Split A
      AblasComplexSplitLength(a,m,s1,s2);
      //--- check
      if(isUpper && optype==0)
        {
         //---           (A1  A12)-1  ( X1 )
         //--- A^-1*X* = (       )   *(    )
         //---           (     A2)    ( X2 )
         CMatrixLeftTrsM(s2,n,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2+s1,j2);
         CMatrixGemm(s1,n,s2,Alpha,a,i1,j1+s1,0,x,i2+s1,j2,0,Beta,x,i2,j2);
         CMatrixLeftTrsM(s1,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         //--- exit the function
         return;
        }
      //--- check
      if(isUpper && optype!=0)
        {
         //---          (A1'     )-1 ( X1 )
         //--- A^-1*X = (        )  *(    )
         //---          (A12' A2')   ( X2 )
         CMatrixLeftTrsM(s1,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         CMatrixGemm(s2,n,s1,Alpha,a,i1,j1+s1,optype,x,i2,j2,0,Beta,x,i2+s1,j2);
         CMatrixLeftTrsM(s2,n,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2+s1,j2);
         //--- exit the function
         return;
        }
      //--- check
      if(!isUpper && optype==0)
        {
         //---          (A1     )-1 ( X1 )
         //--- A^-1*X = (       )  *(    )
         //---          (A21  A2)   ( X2 )
         CMatrixLeftTrsM(s1,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         CMatrixGemm(s2,n,s1,Alpha,a,i1+s1,j1,0,x,i2,j2,0,Beta,x,i2+s1,j2);
         CMatrixLeftTrsM(s2,n,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2+s1,j2);
         //--- exit the function
         return;
        }
      //--- check
      if(!isUpper && optype!=0)
        {
         //---          (A1' A21')-1 ( X1 )
         //--- A^-1*X = (        )  *(    )
         //---          (     A2')   ( X2 )
         CMatrixLeftTrsM(s2,n,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2+s1,j2);
         CMatrixGemm(s1,n,s2,Alpha,a,i1+s1,j1,optype,x,i2+s1,j2,0,Beta,x,i2,j2);
         CMatrixLeftTrsM(s1,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Same as CMatrixRightTRSM, but for real matrices                  |
//| OpType may be only 0 or 1.                                       |
//+------------------------------------------------------------------+
static void CAblas::RMatrixRightTrsM(const int m,const int n,CMatrixDouble &a,
                                     const int i1,const int j1,const bool isUpper,
                                     const bool isUnit,const int optype,
                                     CMatrixDouble &x,const int i2,const int j2)
  {
//--- create variables
   int s1=0;
   int s2=0;
   int bs=AblasBlockSize();
//--- check
   if(m<=bs && n<=bs)
     {
      //--- basic algorithm
      RMatrixRightTrsM2(m,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
      //--- exit the function
      return;
     }
//--- check
   if(m>=n)
     {
      //--- Split X: X*A = (X1 X2)^T*A
      AblasSplitLength(a,m,s1,s2);
      RMatrixRightTrsM(s1,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
      RMatrixRightTrsM(s2,n,a,i1,j1,isUpper,isUnit,optype,x,i2+s1,j2);
     }
   else
     {
      //--- Split A:
      //---               (A1  A12)
      //--- X*op(A) = X*op(       )
      //---               (     A2)
      //--- Different variants depending on
      //--- IsUpper/OpType combinations
      AblasSplitLength(a,n,s1,s2);
      //--- check
      if(isUpper && optype==0)
        {
         //---                  (A1  A12)-1
         //--- X*A^-1 = (X1 X2)*(       )
         //---                  (     A2)
         RMatrixRightTrsM(m,s1,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         RMatrixGemm(m,s2,s1,-1.0,x,i2,j2,0,a,i1,j1+s1,0,1.0,x,i2,j2+s1);
         RMatrixRightTrsM(m,s2,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2,j2+s1);
         //--- exit the function
         return;
        }
      //--- check
      if(isUpper && optype!=0)
        {
         //---                  (A1'     )-1
         //--- X*A^-1 = (X1 X2)*(        )
         //---                  (A12' A2')
         RMatrixRightTrsM(m,s2,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2,j2+s1);
         RMatrixGemm(m,s1,s2,-1.0,x,i2,j2+s1,0,a,i1,j1+s1,optype,1.0,x,i2,j2);
         RMatrixRightTrsM(m,s1,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         //--- exit the function
         return;
        }
      //--- check
      if(!isUpper && optype==0)
        {
         //---                  (A1     )-1
         //--- X*A^-1 = (X1 X2)*(       )
         //---                  (A21  A2)
         RMatrixRightTrsM(m,s2,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2,j2+s1);
         RMatrixGemm(m,s1,s2,-1.0,x,i2,j2+s1,0,a,i1+s1,j1,0,1.0,x,i2,j2);
         RMatrixRightTrsM(m,s1,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         //--- exit the function
         return;
        }
      //--- check
      if(!isUpper && optype!=0)
        {
         //---                  (A1' A21')-1
         //--- X*A^-1 = (X1 X2)*(        )
         //---                  (     A2')
         RMatrixRightTrsM(m,s1,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         RMatrixGemm(m,s2,s1,-1.0,x,i2,j2,0,a,i1+s1,j1,optype,1.0,x,i2,j2+s1);
         RMatrixRightTrsM(m,s2,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2,j2+s1);
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Same as CMatrixLeftTRSM, but for real matrices                   |
//| OpType may be only 0 or 1.                                       |
//+------------------------------------------------------------------+
static void CAblas::RMatrixLeftTrsM(const int m,const int n,CMatrixDouble &a,
                                    const int i1,const int j1,const bool isUpper,
                                    const bool isUnit,const int optype,
                                    CMatrixDouble &x,const int i2,const int j2)
  {
//--- create variables
   int s1=0;
   int s2=0;
   int bs=AblasBlockSize();
//--- check
   if(m<=bs && n<=bs)
     {
      //--- basic algorithm
      RMatrixLeftTrsM2(m,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
      //--- exit the function
      return;
     }
//--- check
   if(n>=m)
     {
      //--- Split X: op(A)^-1*X = op(A)^-1*(X1 X2)
      AblasSplitLength(x,n,s1,s2);
      RMatrixLeftTrsM(m,s1,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
      RMatrixLeftTrsM(m,s2,a,i1,j1,isUpper,isUnit,optype,x,i2,j2+s1);
     }
   else
     {
      //--- Split A
      AblasSplitLength(a,m,s1,s2);
      //--- check
      if(isUpper && optype==0)
        {
         //---           (A1  A12)-1  ( X1 )
         //--- A^-1*X* = (       )   *(    )
         //---           (     A2)    ( X2 )
         RMatrixLeftTrsM(s2,n,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2+s1,j2);
         RMatrixGemm(s1,n,s2,-1.0,a,i1,j1+s1,0,x,i2+s1,j2,0,1.0,x,i2,j2);
         RMatrixLeftTrsM(s1,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         //--- exit the function
         return;
        }
      //--- check
      if(isUpper && optype!=0)
        {
         //---          (A1'     )-1 ( X1 )
         //--- A^-1*X = (        )  *(    )
         //---          (A12' A2')   ( X2 )
         RMatrixLeftTrsM(s1,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         RMatrixGemm(s2,n,s1,-1.0,a,i1,j1+s1,optype,x,i2,j2,0,1.0,x,i2+s1,j2);
         RMatrixLeftTrsM(s2,n,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2+s1,j2);
         //--- exit the function
         return;
        }
      //--- check
      if(!isUpper && optype==0)
        {
         //---          (A1     )-1 ( X1 )
         //--- A^-1*X = (       )  *(    )
         //---          (A21  A2)   ( X2 )
         RMatrixLeftTrsM(s1,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
         RMatrixGemm(s2,n,s1,-1.0,a,i1+s1,j1,0,x,i2,j2,0,1.0,x,i2+s1,j2);
         RMatrixLeftTrsM(s2,n,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2+s1,j2);
         //--- exit the function
         return;
        }
      //--- check
      if(!isUpper && optype!=0)
        {
         //---          (A1' A21')-1 ( X1 )
         //--- A^-1*X = (        )  *(    )
         //---          (     A2')   ( X2 )
         RMatrixLeftTrsM(s2,n,a,i1+s1,j1+s1,isUpper,isUnit,optype,x,i2+s1,j2);
         RMatrixGemm(s1,n,s2,-1.0,a,i1+s1,j1,optype,x,i2+s1,j2,0,1.0,x,i2,j2);
         RMatrixLeftTrsM(s1,n,a,i1,j1,isUpper,isUnit,optype,x,i2,j2);
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| This subroutine calculates  C=alpha*A*A^H+beta*C or              |
//| C=alpha*A^H*A+beta*C where:                                      |
//| * C is NxN Hermitian matrix given by its upper/lower triangle    |
//| * A is NxK matrix when A*A^H is calculated, KxN matrix otherwise |
//| Additional info:                                                 |
//| * cache-oblivious algorithm is used.                             |
//| * multiplication result replaces C. If Beta=0, C elements are not|
//|   used in calculations (not multiplied by zero - just not        |
//|   referenced)                                                    |
//| * if Alpha=0, A is not used (not multiplied by zero - just not   |
//|   referenced)                                                    |
//| * if both Beta and Alpha are zero, C is filled by zeros.         |
//| INPUT PARAMETERS                                                 |
//|     N       -   matrix size, N>=0                                |
//|     K       -   matrix size, K>=0                                |
//|     Alpha   -   coefficient                                      |
//|     A       -   matrix                                           |
//|     IA      -   submatrix offset                                 |
//|     JA      -   submatrix offset                                 |
//|     OpTypeA -   multiplication type:                             |
//|                 * 0 - A*A^H is calculated                        |
//|                 * 2 - A^H*A is calculated                        |
//|     Beta    -   coefficient                                      |
//|     C       -   matrix                                           |
//|     IC      -   submatrix offset                                 |
//|     JC      -   submatrix offset                                 |
//|     IsUpper -   whether C is upper triangular or lower triangular|
//+------------------------------------------------------------------+
static void CAblas::CMatrixSyrk(const int n,const int k,const double alpha,
                                CMatrixComplex &a,const int ia,const int ja,
                                const int optypea,const double beta,CMatrixComplex &c,
                                const int ic,const int jc,const bool isUpper)
  {
//--- create variables
   al_complex Alpha(alpha,0);
   al_complex Beta(beta,0);
   int        s1=0;
   int        s2=0;
   int        bs=AblasComplexBlockSize();
//--- check
   if(n<=bs && k<=bs)
     {
      //--- basic algorithm
      CMatrixSyrk2(n,k,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
      //--- exit the function
      return;
     }
//--- check
   if(k>=n)
     {
      //--- Split K
      AblasComplexSplitLength(a,k,s1,s2);
      //--- check
      if(optypea==0)
        {
         CMatrixSyrk(n,s1,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
         CMatrixSyrk(n,s2,alpha,a,ia,ja+s1,optypea,1.0,c,ic,jc,isUpper);
        }
      else
        {
         CMatrixSyrk(n,s1,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
         CMatrixSyrk(n,s2,alpha,a,ia+s1,ja,optypea,1.0,c,ic,jc,isUpper);
        }
     }
   else
     {
      //--- Split N
      AblasComplexSplitLength(a,n,s1,s2);
      //--- check
      if(optypea==0 && isUpper)
        {
         CMatrixSyrk(s1,k,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
         CMatrixGemm(s1,s2,k,Alpha,a,ia,ja,0,a,ia+s1,ja,2,Beta,c,ic,jc+s1);
         CMatrixSyrk(s2,k,alpha,a,ia+s1,ja,optypea,beta,c,ic+s1,jc+s1,isUpper);
         //--- exit the function
         return;
        }
      //--- check
      if(optypea==0 && !isUpper)
        {
         CMatrixSyrk(s1,k,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
         CMatrixGemm(s2,s1,k,Alpha,a,ia+s1,ja,0,a,ia,ja,2,Beta,c,ic+s1,jc);
         CMatrixSyrk(s2,k,alpha,a,ia+s1,ja,optypea,beta,c,ic+s1,jc+s1,isUpper);
         //--- exit the function
         return;
        }
      //--- check
      if(optypea!=0 && isUpper)
        {
         CMatrixSyrk(s1,k,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
         CMatrixGemm(s1,s2,k,Alpha,a,ia,ja,2,a,ia,ja+s1,0,Beta,c,ic,jc+s1);
         CMatrixSyrk(s2,k,alpha,a,ia,ja+s1,optypea,beta,c,ic+s1,jc+s1,isUpper);
         //--- exit the function
         return;
        }
      //--- check
      if(optypea!=0 && !isUpper)
        {
         CMatrixSyrk(s1,k,alpha,a,ia,ja,optypea,beta,c,ic,jc,isUpper);
         CMatrixGemm(s2,s1,k,Alpha,a,ia,ja+s1,2,a,ia,ja,0,Beta,c,ic+s1,jc);
         CMatrixSyrk(s2,k,alpha,a,ia,ja+s1,optypea,beta,c,ic+s1,jc+s1,isUpper);
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| This subroutine calculates C = alpha*op1(A)*op2(B) +beta*C where:|
//| * C is MxN general matrix                                        |
//| * op1(A) is MxK matrix                                           |
//| * op2(B) is KxN matrix                                           |
//| * "op" may be identity transformation, transposition, conjugate  |
//| transposition                                                    |
//| Additional info:                                                 |
//| * cache-oblivious algorithm is used.                             |
//| * multiplication result replaces C. If Beta=0, C elements are not|
//|   used in calculations (not multiplied by zero - just not        |
//|   referenced)                                                    |
//| * if Alpha=0, A is not used (not multiplied by zero - just not   |
//|   referenced)                                                    |
//| * if both Beta and Alpha are zero, C is filled by zeros.         |
//| INPUT PARAMETERS                                                 |
//|     N       -   matrix size, N>0                                 |
//|     M       -   matrix size, N>0                                 |
//|     K       -   matrix size, K>0                                 |
//|     Alpha   -   coefficient                                      |
//|     A       -   matrix                                           |
//|     IA      -   submatrix offset                                 |
//|     JA      -   submatrix offset                                 |
//|     OpTypeA -   transformation type:                             |
//|                 * 0 - no transformation                          |
//|                 * 1 - transposition                              |
//|                 * 2 - conjugate transposition                    |
//|     B       -   matrix                                           |
//|     IB      -   submatrix offset                                 |
//|     JB      -   submatrix offset                                 |
//|     OpTypeB -   transformation type:                             |
//|                 * 0 - no transformation                          |
//|                 * 1 - transposition                              |
//|                 * 2 - conjugate transposition                    |
//|     Beta    -   coefficient                                      |
//|     C       -   matrix                                           |
//|     IC      -   submatrix offset                                 |
//|     JC      -   submatrix offset                                 |
//+------------------------------------------------------------------+
static void CAblas::CMatrixGemm(const int m,const int n,const int k,al_complex &alpha,
                                CMatrixComplex &a,const int ia,const int ja,
                                const int optypea,CMatrixComplex &b,const int ib,
                                const int jb,const int optypeb,al_complex &beta,
                                CMatrixComplex &c,const int ic,const int jc)
  {
//--- create variables
   al_complex Beta(1,0);
   int        s1=0;
   int        s2=0;
   int        bs=AblasComplexBlockSize();
//--- check
   if(m<=bs && n<=bs && k<=bs)
     {
      //--- basic algorithm
      CMatrixGemmk(m,n,k,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
      //--- exit the function
      return;
     }
//--- check
   if(m>=n && m>=k)
     {
      //--- A*B = (A1 A2)^T*B
      AblasComplexSplitLength(a,m,s1,s2);
      CMatrixGemm(s1,n,k,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
      //--- check
      if(optypea==0)
         CMatrixGemm(s2,n,k,alpha,a,ia+s1,ja,optypea,b,ib,jb,optypeb,beta,c,ic+s1,jc);
      else
         CMatrixGemm(s2,n,k,alpha,a,ia,ja+s1,optypea,b,ib,jb,optypeb,beta,c,ic+s1,jc);
      //--- exit the function
      return;
     }
//--- check
   if(n>=m && n>=k)
     {
      //---A*B = A*(B1 B2)
      AblasComplexSplitLength(a,n,s1,s2);
      //--- check
      if(optypeb==0)
        {
         CMatrixGemm(m,s1,k,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         CMatrixGemm(m,s2,k,alpha,a,ia,ja,optypea,b,ib,jb+s1,optypeb,beta,c,ic,jc+s1);
        }
      else
        {
         CMatrixGemm(m,s1,k,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         CMatrixGemm(m,s2,k,alpha,a,ia,ja,optypea,b,ib+s1,jb,optypeb,beta,c,ic,jc+s1);
        }
      //--- exit the function
      return;
     }
//--- check
   if(k>=m && k>=n)
     {
      //--- A*B = (A1 A2)*(B1 B2)^T
      AblasComplexSplitLength(a,k,s1,s2);
      //--- check
      if(optypea==0 && optypeb==0)
        {
         CMatrixGemm(m,n,s1,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         CMatrixGemm(m,n,s2,alpha,a,ia,ja+s1,optypea,b,ib+s1,jb,optypeb,Beta,c,ic,jc);
        }
      //--- check
      if(optypea==0 && optypeb!=0)
        {
         CMatrixGemm(m,n,s1,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         CMatrixGemm(m,n,s2,alpha,a,ia,ja+s1,optypea,b,ib,jb+s1,optypeb,Beta,c,ic,jc);
        }
      //--- check
      if(optypea!=0 && optypeb==0)
        {
         CMatrixGemm(m,n,s1,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         CMatrixGemm(m,n,s2,alpha,a,ia+s1,ja,optypea,b,ib+s1,jb,optypeb,Beta,c,ic,jc);
        }
      //--- check
      if(optypea!=0 && optypeb!=0)
        {
         CMatrixGemm(m,n,s1,alpha,a,ia,ja,optypea,b,ib,jb,optypeb,beta,c,ic,jc);
         CMatrixGemm(m,n,s2,alpha,a,ia+s1,ja,optypea,b,ib,jb+s1,optypeb,Beta,c,ic,jc);
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Level 2 variant of CMatrixRightTRSM                              |
//+------------------------------------------------------------------+
static void CAblas::CMatrixRightTrsM2(const int m,const int n,CMatrixComplex &a,
                                      const int i1,const int j1,const bool isUpper,
                                      const bool isUnit,const int optype,
                                      CMatrixComplex &x,const int i2,const int j2)
  {
//--- check
   if(n*m==0)
      return;
//--- create variables
   int        i=0;
   int        j=0;
   al_complex vc=0;
   al_complex vd=0;
   int        i_=0;
   int        i1_=0;
//--- General case
   if(isUpper)
     {
      //--- Upper triangular matrix
      if(optype==0)
        {
         //--- X*A^(-1)
         for(i=0;i<m;i++)
           {
            for(j=0;j<n;j++)
              {
               //--- check
               if(isUnit)
                  vd=1;
               else
                  vd=a[i1+j][j1+j];
               //--- change x
               x[i2+i].Set(j2+j,x[i2+i][j2+j]/vd);
               //--- check
               if(j<n-1)
                 {
                  vc=x[i2+i][j2+j];
                  i1_=j1-j2;
                  for(i_=j2+j+1;i_<j2+n;i_++)
                     x[i2+i].Set(i_,x[i2+i][i_]-vc*a[i1+j][i_+i1_]);
                 }
              }
           }
         //--- exit the function
         return;
        }
      //--- check
      if(optype==1)
        {
         //--- X*A^(-T)
         for(i=0;i<m;i++)
           {
            for(j=n-1;j>=0;j--)
              {
               vc=0;
               vd=1;
               //--- check
               if(j<n-1)
                 {
                  i1_=j1-j2;
                  vc=0.0;
                  for(i_=j2+j+1;i_<j2+n;i_++)
                     vc+=x[i2+i][i_]*a[i1+j][i_+i1_];
                 }
               //--- check
               if(!isUnit)
                  vd=a[i1+j][j1+j];
               //--- change x
               x[i2+i].Set(j2+j,(x[i2+i][j2+j]-vc)/vd);
              }
           }
         //--- exit the function
         return;
        }
      //--- check
      if(optype==2)
        {
         //--- X*A^(-H)
         for(i=0;i<m;i++)
           {
            for(j=n-1;j>=0;j--)
              {
               vc=0;
               vd=1;
               //--- check
               if(j<n-1)
                 {
                  i1_=j1-j2;
                  vc=0.0;
                  for(i_=j2+j+1;i_<j2+n;i_++)
                     vc+=x[i2+i][i_]*CMath::Conj(a[i1+j][i_+i1_]);
                 }
               //--- check
               if(!isUnit)
                  vd=CMath::Conj(a[i1+j][j1+j]);
               //--- change x
               x[i2+i].Set(j2+j,(x[i2+i][j2+j]-vc)/vd);
              }
           }
         //--- exit the function
         return;
        }
     }
   else
     {
      //--- Lower triangular matrix
      if(optype==0)
        {
         //--- X*A^(-1)
         for(i=0;i<m;i++)
           {
            for(j=n-1;j>=0;j--)
              {
               //--- check
               if(isUnit)
                  vd=1;
               else
                  vd=a[i1+j][j1+j];
               //--- change x
               x[i2+i].Set(j2+j,x[i2+i][j2+j]/vd);
               //--- check
               if(j>0)
                 {
                  vc=x[i2+i][j2+j];
                  i1_=j1-j2;
                  for(i_=j2;i_<j2+j;i_++)
                     x[i2+i].Set(i_,x[i2+i][i_]-vc*a[i1+j][i_+i1_]);
                 }
              }
           }
         //--- exit the function
         return;
        }
      if(optype==1)
        {
         //--- X*A^(-T)
         for(i=0;i<m;i++)
           {
            for(j=0;j<n;j++)
              {
               vc=0;
               vd=1;
               //--- check
               if(j>0)
                 {
                  i1_=j1-j2;
                  vc=0.0;
                  for(i_=j2;i_<j2+j;i_++)
                     vc+=x[i2+i][i_]*a[i1+j][i_+i1_];
                 }
               //--- check
               if(!isUnit)
                  vd=a[i1+j][j1+j];
               //--- change x
               x[i2+i].Set(j2+j,(x[i2+i][j2+j]-vc)/vd);
              }
           }
         //--- exit the function
         return;
        }
      if(optype==2)
        {
         //--- X*A^(-H)
         for(i=0;i<m;i++)
           {
            for(j=0;j<n;j++)
              {
               vc=0;
               vd=1;
               //--- check
               if(j>0)
                 {
                  i1_=j1-j2;
                  vc=0.0;
                  for(i_=j2;i_<j2+j;i_++)
                     vc+=x[i2+i][i_]*CMath::Conj(a[i1+j][i_+i1_]);
                 }
               //--- check
               if(!isUnit)
                  vd=CMath::Conj(a[i1+j][j1+j]);
               //--- change x
               x[i2+i].Set(j2+j,(x[i2+i][j2+j]-vc)/vd);
              }
           }
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Level-2 subroutine                                               |
//+------------------------------------------------------------------+
static void CAblas::CMatrixLeftTrsM2(const int m,const int n,CMatrixComplex &a,
                                     const int i1,const int j1,const bool isUpper,
                                     const bool isUnit,const int optype,
                                     CMatrixComplex &x,const int i2,const int j2)
  {
//--- check
   if(n*m==0)
      return;
//--- create variables
   al_complex Beta(1,0);
   int        i=0;
   int        j=0;
   al_complex vc=0;
   al_complex vd=0;
   int        i_=0;
//--- General case
   if(isUpper)
     {
      //--- Upper triangular matrix
      if(optype==0)
        {
         //--- A^(-1)*X
         for(i=m-1;i>=0;i--)
           {
            for(j=i+1;j<=m-1;j++)
              {
               vc=a[i1+i][j1+j];
               //--- change x
               for(i_=j2;i_<j2+n;i_++)
                  x[i2+i].Set(i_,x[i2+i][i_]-vc*x[i2+j][i_]);
              }
            //--- check
            if(!isUnit)
              {
               vd=Beta/a[i1+i][j1+i];
               //--- change x
               for(i_=j2;i_<j2+n;i_++)
                  x[i2+i].Set(i_,vd*x[i2+i][i_]);
              }
           }
         //--- exit the function
         return;
        }
      //--- check
      if(optype==1)
        {
         //--- A^(-T)*X
         for(i=0;i<m;i++)
           {
            //--- check
            if(isUnit)
               vd=1;
            else
               vd=Beta/a[i1+i][j1+i];
            //--- change x
            for(i_=j2;i_<j2+n;i_++)
               x[i2+i].Set(i_,vd*x[i2+i][i_]);
            for(j=i+1;j<=m-1;j++)
              {
               vc=a[i1+i][j1+j];
               for(i_=j2;i_<j2+n;i_++)
                  x[i2+j].Set(i_,x[i2+j][i_]-vc*x[i2+i][i_]);
              }
           }
         //--- exit the function
         return;
        }
      //--- check
      if(optype==2)
        {
         //--- A^(-H)*X
         for(i=0;i<m;i++)
           {
            //--- check
            if(isUnit)
               vd=1;
            else
               vd=Beta/CMath::Conj(a[i1+i][j1+i]);
            //--- change x
            for(i_=j2;i_<j2+n;i_++)
               x[i2+i].Set(i_,vd*x[i2+i][i_]);
            for(j=i+1;j<=m-1;j++)
              {
               vc=CMath::Conj(a[i1+i][j1+j]);
               for(i_=j2;i_<j2+n;i_++)
                  x[i2+j].Set(i_,x[i2+j][i_]-vc*x[i2+i][i_]);
              }
           }
         //--- exit the function
         return;
        }
     }
   else
     {
      //--- Lower triangular matrix
      if(optype==0)
        {
         //--- A^(-1)*X
         for(i=0;i<m;i++)
           {
            //--- change x
            for(j=0;j<i;j++)
              {
               vc=a[i1+i][j1+j];
               for(i_=j2;i_<j2+n;i_++)
                  x[i2+i].Set(i_,x[i2+i][i_]-vc*x[i2+j][i_]);
              }
            //--- check
            if(isUnit)
               vd=1;
            else
               vd=Beta/a[i1+j][j1+j];
            //--- change x
            for(i_=j2;i_<j2+n;i_++)
               x[i2+i].Set(i_,vd*x[i2+i][i_]);
           }
         //--- exit the function
         return;
        }
      if(optype==1)
        {
         //--- A^(-T)*X
         for(i=m-1;i>=0;i--)
           {
            //--- check
            if(isUnit)
               vd=1;
            else
               vd=Beta/a[i1+i][j1+i];
            //--- change x
            for(i_=j2;i_<j2+n;i_++)
               x[i2+i].Set(i_,vd*x[i2+i][i_]);
            for(j=i-1;j>=0;j--)
              {
               vc=a[i1+i][j1+j];
               for(i_=j2;i_<j2+n;i_++)
                  x[i2+j].Set(i_,x[i2+j][i_]-vc*x[i2+i][i_]);
              }
           }
         //--- exit the function
         return;
        }
      if(optype==2)
        {
         //--- A^(-H)*X
         for(i=m-1;i>=0;i--)
           {
            //--- check
            if(isUnit)
               vd=1;
            else
               vd=Beta/CMath::Conj(a[i1+i][j1+i]);
            //--- change x
            for(i_=j2;i_<j2+n;i_++)
               x[i2+i].Set(i_,vd*x[i2+i][i_]);
            for(j=i-1;j>=0;j--)
              {
               vc=CMath::Conj(a[i1+i][j1+j]);
               for(i_=j2;i_<j2+n;i_++)
                  x[i2+j].Set(i_,x[i2+j][i_]-vc*x[i2+i][i_]);
              }
           }
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Level 2 subroutine                                               |
//+------------------------------------------------------------------+
static void CAblas::RMatrixRightTrsM2(const int m,const int n,CMatrixDouble &a,
                                      const int i1,const int j1,const bool isUpper,
                                      const bool isUnit,const int optype,
                                      CMatrixDouble &x,const int i2,const int j2)
  {
//--- check
   if(n*m==0)
      return;
//--- create variables 
   int    i=0;
   int    j=0;
   double vr=0;
   double vd=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(isUpper)
     {
      //--- Upper triangular matrix
      if(optype==0)
        {
         //--- X*A^(-1)
         for(i=0;i<m;i++)
           {
            for(j=0;j<n;j++)
              {
               //--- check
               if(isUnit)
                  vd=1;
               else
                  vd=a[i1+j][j1+j];
               //--- change x
               x[i2+i].Set(j2+j,x[i2+i][j2+j]/vd);
               //--- check
               if(j<n-1)
                 {
                  vr=x[i2+i][j2+j];
                  i1_=j1-j2;
                  //--- change x
                  for(i_=j2+j+1;i_<j2+n;i_++)
                     x[i2+i].Set(i_,x[i2+i][i_]-vr*a[i1+j][i_+i1_]);
                 }
              }
           }
         //--- exit the function
         return;
        }
      //--- check
      if(optype==1)
        {
         //--- X*A^(-T)
         for(i=0;i<m;i++)
           {
            for(j=n-1;j>=0;j--)
              {
               vr=0;
               vd=1;
               //--- check
               if(j<n-1)
                 {
                  i1_=j1-j2;
                  vr=0.0;
                  for(i_=j2+j+1;i_<j2+n;i_++)
                     vr+=x[i2+i][i_]*a[i1+j][i_+i1_];
                 }
               //--- check
               if(!isUnit)
                  vd=a[i1+j][j1+j];
               //--- change x
               x[i2+i].Set(j2+j,(x[i2+i][j2+j]-vr)/vd);
              }
           }
         //--- exit the function
         return;
        }
     }
   else
     {
      //--- Lower triangular matrix
      if(optype==0)
        {
         //--- X*A^(-1)
         for(i=0;i<m;i++)
           {
            for(j=n-1;j>=0;j--)
              {
               //--- check
               if(isUnit)
                  vd=1;
               else
                  vd=a[i1+j][j1+j];
               //--- change x
               x[i2+i].Set(j2+j,x[i2+i][j2+j]/vd);
               //--- check
               if(j>0)
                 {
                  vr=x[i2+i][j2+j];
                  i1_=j1-j2;
                  //--- change x
                  for(i_=j2;i_<j2+j;i_++)
                     x[i2+i].Set(i_,x[i2+i][i_]-vr*a[i1+j][i_+i1_]);
                 }
              }
           }
         //--- exit the function
         return;
        }
      //--- check
      if(optype==1)
        {
         //--- X*A^(-T)
         for(i=0;i<m;i++)
           {
            for(j=0;j<n;j++)
              {
               vr=0;
               vd=1;
               //--- check
               if(j>0)
                 {
                  i1_=j1-j2;
                  vr=0.0;
                  for(i_=j2;i_<j2+j;i_++)
                     vr+=x[i2+i][i_]*a[i1+j][i_+i1_];
                 }
               //--- check
               if(!isUnit)
                  vd=a[i1+j][j1+j];
               //--- change x
               x[i2+i].Set(j2+j,(x[i2+i][j2+j]-vr)/vd);
              }
           }
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Level 2 subroutine                                               |
//+------------------------------------------------------------------+
static void CAblas::RMatrixLeftTrsM2(const int m,const int n,CMatrixDouble &a,
                                     const int i1,const int j1,const bool isUpper,
                                     const bool isUnit,const int optype,
                                     CMatrixDouble &x,const int i2,const int j2)
  {
//--- check
   if(n*m==0)
      return;
//--- create variables
   int    i=0;
   int    j=0;
   double vr=0;
   double vd=0;
   int    i_=0;
//--- check
   if(isUpper)
     {
      //--- Upper triangular matrix
      if(optype==0)
        {
         //--- A^(-1)*X
         for(i=m-1;i>=0;i--)
           {
            for(j=i+1;j<=m-1;j++)
              {
               vr=a[i1+i][j1+j];
               for(i_=j2;i_<j2+n;i_++)
                  x[i2+i].Set(i_,x[i2+i][i_]-vr*x[i2+j][i_]);
              }
            //--- check
            if(!isUnit)
              {
               vd=1/a[i1+i][j1+i];
               for(i_=j2;i_<j2+n;i_++)
                  x[i2+i].Set(i_,vd*x[i2+i][i_]);
              }
           }
         //--- exit the function
         return;
        }
      //--- check
      if(optype==1)
        {
         //--- A^(-T)*X
         for(i=0;i<m;i++)
           {
            //--- check
            if(isUnit)
               vd=1;
            else
               vd=1/a[i1+i][j1+i];
            //--- change x
            for(i_=j2;i_<j2+n;i_++)
               x[i2+i].Set(i_,vd*x[i2+i][i_]);
            for(j=i+1;j<=m-1;j++)
              {
               vr=a[i1+i][j1+j];
               for(i_=j2;i_<j2+n;i_++)
                  x[i2+j].Set(i_,x[i2+j][i_]-vr*x[i2+i][i_]);
              }
           }
         //--- exit the function
         return;
        }
     }
   else
     {
      //--- Lower triangular matrix
      if(optype==0)
        {
         //--- A^(-1)*X
         for(i=0;i<m;i++)
           {
            for(j=0;j<i;j++)
              {
               vr=a[i1+i][j1+j];
               for(i_=j2;i_<j2+n;i_++)
                  x[i2+i].Set(i_,x[i2+i][i_]-vr*x[i2+j][i_]);
              }
            //--- check
            if(isUnit)
               vd=1;
            else
               vd=1/a[i1+j][j1+j];
            for(i_=j2;i_<j2+n;i_++)
               x[i2+i].Set(i_,vd*x[i2+i][i_]);
           }
         //--- exit the function
         return;
        }
      //--- check
      if(optype==1)
        {
         //--- A^(-T)*X
         for(i=m-1;i>=0;i--)
           {
            //--- check
            if(isUnit)
               vd=1;
            else
               vd=1/a[i1+i][j1+i];
            //--- change x
            for(i_=j2;i_<j2+n;i_++)
               x[i2+i].Set(i_,vd*x[i2+i][i_]);
            for(j=i-1;j>=0;j--)
              {
               vr=a[i1+i][j1+j];
               for(i_=j2;i_<j2+n;i_++)
                  x[i2+j].Set(i_,x[i2+j][i_]-vr*x[i2+i][i_]);
              }
           }
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Level 2 subroutine                                               |
//+------------------------------------------------------------------+
static void CAblas::CMatrixSyrk2(const int n,const int k,const double alpha,
                                 const CMatrixComplex &a,const int ia,const int ja,
                                 const int optypea,const double beta,CMatrixComplex &c,
                                 const int ic,const int jc,const bool isUpper)
  {
//--- check
   if((alpha==0 || k==0) && beta==1)
      return;
//--- create variables
   al_complex Alpha(alpha,0);
   al_complex Beta(beta,0);
   al_complex Zero(0,0);
   int        i=0;
   int        j=0;
   int        j1=0;
   int        j2=0;
   al_complex v=0;
   int        i_=0;
   int        i1_=0;
//--- check
   if(optypea==0)
     {
      //--- C = alpha*A*A^H+beta*C
      for(i=0;i<n;i++)
        {
         //--- check
         if(isUpper)
           {
            j1=i;
            j2=n-1;
           }
         else
           {
            j1=0;
            j2=i;
           }
         //--- cycle
         for(j=j1;j<=j2;j++)
           {
            //--- check
            if(alpha!=0 && k>0)
              {
               v=0.0;
               for(i_=ja;i_<=ja+k-1;i_++)
                  v+=a[ia+i][i_]*CMath::Conj(a[ia+j][i_]);
              }
            else
               v=0;
            //--- check
            if(beta==0)
               c[ic+i].Set(jc+j,Alpha*v);
            else
               c[ic+i].Set(jc+j,Beta*c[ic+i][jc+j]+Alpha*v);
           }
        }
      //--- exit the function
      return;
     }
   else
     {
      //--- C = alpha*A^H*A+beta*C
      for(i=0;i<n;i++)
        {
         //--- check
         if(isUpper)
           {
            j1=i;
            j2=n-1;
           }
         else
           {
            j1=0;
            j2=i;
           }
         //--- check
         if(beta==0)
           {
            for(j=j1;j<=j2;j++)
               c[ic+i].Set(jc+j,Zero);
           }
         else
           {
            for(i_=jc+j1;i_<=jc+j2;i_++)
               c[ic+i].Set(i_,Beta*c[ic+i][i_]);
           }
        }
      //--- cycle
      for(i=0;i<k;i++)
        {
         for(j=0;j<n;j++)
           {
            //--- check
            if(isUpper)
              {
               j1=j;
               j2=n-1;
              }
            else
              {
               j1=0;
               j2=j;
              }
            v=Alpha*CMath::Conj(a[ia+i][ja+j]);
            i1_=(ja+j1)-(jc+j1);
            for(i_=jc+j1;i_<=jc+j2;i_++)
               c[ic+j].Set(i_,c[ic+j][i_]+v*a[ia+i][i_+i1_]);
           }
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| GEMM kernel                                                      |
//+------------------------------------------------------------------+
static void CAblas::CMatrixGemmk(const int m,const int n,const int k,al_complex &alpha,
                                 const CMatrixComplex &a,const int ia,const int ja,
                                 const int optypea,const CMatrixComplex &b,const int ib,
                                 const int jb,const int optypeb,al_complex &beta,
                                 CMatrixComplex &c,const int ic,const int jc)
  {
//--- check
   if(m*n==0)
      return;
//--- create variables
   al_complex Zero(0,0);
   int        i=0;
   int        j=0;
   al_complex v=0;
   int        i_=0;
   int        i1_=0;
//--- special case
   if(k==0)
     {
      //--- check
      if(beta!=Zero)
        {
         //--- get c
         for(i=0;i<m;i++)
            for(j=0;j<n;j++)
               c[ic+i].Set(jc+j,beta*c[ic+i][jc+j]);
        }
      else
        {
         //--- get c
         for(i=0;i<m;i++)
            for(j=0;j<n;j++)
               c[ic+i].Set(jc+j,Zero);
        }
      //--- exit the function
      return;
     }
//--- General case
   if(optypea==0 && optypeb!=0)
     {
      //--- A*B'
      for(i=0;i<m;i++)
        {
         for(j=0;j<n;j++)
           {
            //--- check
            if(k==0 || alpha==Zero)
               v=0;
            else
              {
               //--- check
               if(optypeb==1)
                 {
                  i1_=(jb)-(ja);
                  v=0.0;
                  for(i_=ja;i_<=ja+k-1;i_++)
                     v+=a[ia+i][i_]*b[ib+j][i_+i1_];
                 }
               else
                 {
                  i1_=(jb)-(ja);
                  v=0.0;
                  for(i_=ja;i_<=ja+k-1;i_++)
                     v+=a[ia+i][i_]*CMath::Conj(b[ib+j][i_+i1_]);
                 }
              }
            //--- check
            if(beta==Zero)
               c[ic+i].Set(jc+j,alpha*v);
            else
               c[ic+i].Set(jc+j,beta*c[ic+i][jc+j]+alpha*v);
           }
        }
      //--- exit the function
      return;
     }
   if(optypea==0 && optypeb==0)
     {
      //--- A*B
      for(i=0;i<m;i++)
        {
         //--- check
         if(beta!=Zero)
           {
            for(i_=jc;i_<=jc+n-1;i_++)
               c[ic+i].Set(i_,beta*c[ic+i][i_]);
           }
         else
           {
            for(j=0;j<n;j++)
               c[ic+i].Set(jc+j,Zero);
           }
         //--- check
         if(alpha!=Zero)
           {
            for(j=0;j<=k-1;j++)
              {
               v=alpha*a[ia+i][ja+j];
               i1_=(jb)-(jc);
               for(i_=jc;i_<=jc+n-1;i_++)
                  c[ic+i].Set(i_,c[ic+i][i_]+v*b[ib+j][i_+i1_]);
              }
           }
        }
      //--- exit the function
      return;
     }
//--- check
   if(optypea!=0 && optypeb!=0)
     {
      //--- A'*B'
      for(i=0;i<m;i++)
        {
         for(j=0;j<n;j++)
           {
            //--- check
            if(alpha==Zero)
               v=0;
            else
              {
               //--- check
               if(optypea==1)
                 {
                  //--- check
                  if(optypeb==1)
                    {
                     i1_=(jb)-(ia);
                     v=0.0;
                     for(i_=ia;i_<=ia+k-1;i_++)
                        v+=a[i_][ja+i]*b[ib+j][i_+i1_];
                    }
                  else
                    {
                     i1_=(jb)-(ia);
                     v=0.0;
                     for(i_=ia;i_<=ia+k-1;i_++)
                        v+=a[i_][ja+i]*CMath::Conj(b[ib+j][i_+i1_]);
                    }
                 }
               else
                 {
                  //--- check
                  if(optypeb==1)
                    {
                     i1_=(jb)-(ia);
                     v=0.0;
                     for(i_=ia;i_<=ia+k-1;i_++)
                        v+=CMath::Conj(a[i_][ja+i])*b[ib+j][i_+i1_];
                    }
                  else
                    {
                     i1_=(jb)-(ia);
                     v=0.0;
                     for(i_=ia;i_<=ia+k-1;i_++)
                        v+=CMath::Conj(a[i_][ja+i])*CMath::Conj(b[ib+j][i_+i1_]);
                    }
                 }
              }
            //--- check
            if(beta==Zero)
               c[ic+i].Set(jc+j,alpha*v);
            else
               c[ic+i].Set(jc+j,beta*c[ic+i][jc+j]+alpha*v);
           }
        }
      //--- exit the function
      return;
     }
//--- check
   if(optypea!=0 && optypeb==0)
     {
      //--- A'*B
      if(beta==Zero)
        {
         for(i=0;i<m;i++)
            for(j=0;j<n;j++)
               c[ic+i].Set(jc+j,Zero);
        }
      else
        {
         for(i=0;i<m;i++)
            for(i_=jc;i_<=jc+n-1;i_++)
               c[ic+i].Set(i_,beta*c[ic+i][i_]);
        }
      //--- check
      if(alpha!=Zero)
        {
         for(j=0;j<=k-1;j++)
            for(i=0;i<m;i++)
              {
               //--- check
               if(optypea==1)
                  v=alpha*a[ia+j][ja+i];
               else
                  v=alpha*CMath::Conj(a[ia+j][ja+i]);
               i1_=(jb)-(jc);
               for(i_=jc;i_<=jc+n-1;i_++)
                  c[ic+i].Set(i_,c[ic+i][i_]+v*b[ib+j][i_+i1_]);
              }
        }
     }
//--- exit the function
   return;
  }
//+------------------------------------------------------------------+
//| Orthogonal factorizations                                        |
//+------------------------------------------------------------------+
class COrtFac
  {
private:
   //--- private methods
   static void       RMatrixQRBaseCase(CMatrixDouble &a,const int m,const int n,double &work[],double &t[],double &tau[]);
   static void       RMatrixLQBaseCase(CMatrixDouble &a,const int m,const int n,double &work[],double &t[],double &tau[]);
   static void       CMatrixQRBaseCase(CMatrixComplex &a,const int m,const int n,al_complex &work[],al_complex &t[],al_complex &tau[]);
   static void       CMatrixLQBaseCase(CMatrixComplex &a,const int m,const int n,al_complex &work[],al_complex &t[],al_complex &tau[]);
   static void       RMatrixBlockReflector(CMatrixDouble &a,double &tau[],const bool columnwisea,const int lengtha,const int blocksize,CMatrixDouble &t,double &work[]);
   static void       CMatrixBlockReflector(CMatrixComplex &a,al_complex &tau[],const bool columnwisea,const int lengtha,const int blocksize,CMatrixComplex &t,al_complex &work[]);
public:
   //--- constructor, destructor
                     COrtFac(void);
                    ~COrtFac(void);
   //--- real matrix
   static void       RMatrixQR(CMatrixDouble &a,const int m,const int n,double &tau[]);
   static void       RMatrixLQ(CMatrixDouble &a,const int m,const int n,double &tau[]);
   static void       RMatrixQRUnpackQ(CMatrixDouble &a,const int m,const int n,double &tau[],const int qcolumns,CMatrixDouble &q);
   static void       RMatrixQRUnpackR(CMatrixDouble &a,const int m,const int n,CMatrixDouble &r);
   static void       RMatrixLQUnpackQ(CMatrixDouble &a,const int m,const int n,double &tau[],const int qrows,CMatrixDouble &q);
   static void       RMatrixLQUnpackL(CMatrixDouble &a,const int m,const int n,CMatrixDouble &l);
   static void       RMatrixBD(CMatrixDouble &a,const int m,const int n,double &tauq[],double &taup[]);
   static void       RMatrixBDUnpackQ(CMatrixDouble &qp,const int m,const int n,double &tauq[],const int qcolumns,CMatrixDouble &q);
   static void       RMatrixBDMultiplyByQ(CMatrixDouble &qp,const int m,const int n,double &tauq[],CMatrixDouble &z,const int zrows,const int zcolumns,const bool fromtheright,const bool dotranspose);
   static void       RMatrixBDUnpackPT(CMatrixDouble &qp,const int m,const int n,double &taup[],const int ptrows,CMatrixDouble &pt);
   static void       RMatrixBDMultiplyByP(CMatrixDouble &qp,const int m,const int n,double &taup[],CMatrixDouble &z,const int zrows,const int zcolumns,const bool fromtheright,const bool dotranspose);
   static void       RMatrixBDUnpackDiagonals(CMatrixDouble &b,const int m,const int n,bool &isupper,double &d[],double &e[]);
   static void       RMatrixHessenberg(CMatrixDouble &a,const int n,double &tau[]);
   static void       RMatrixHessenbergUnpackQ(CMatrixDouble &a,const int n,double &tau[],CMatrixDouble &q);
   static void       RMatrixHessenbergUnpackH(CMatrixDouble &a,const int n,CMatrixDouble &h);
   static void       SMatrixTD(CMatrixDouble &a,const int n,const bool isupper,double &tau[],double &d[],double &e[]);
   static void       SMatrixTDUnpackQ(CMatrixDouble &a,const int n,const bool isupper,double &tau[],CMatrixDouble &q);
   //--- complex matrix
   static void       CMatrixQR(CMatrixComplex &a,const int m,const int n,al_complex &tau[]);
   static void       CMatrixLQ(CMatrixComplex &a,const int m,const int n,al_complex &tau[]);
   static void       CMatrixQRUnpackQ(CMatrixComplex &a,const int m,const int n,al_complex &tau[],const int qcolumns,CMatrixComplex &q);
   static void       CMatrixQRUnpackR(CMatrixComplex &a,const int m,const int n,CMatrixComplex &r);
   static void       CMatrixLQUnpackQ(CMatrixComplex &a,const int m,const int n,al_complex &tau[],const int qrows,CMatrixComplex &q);
   static void       CMatrixLQUnpackL(CMatrixComplex &a,const int m,const int n,CMatrixComplex &l);
   static void       HMatrixTD(CMatrixComplex &a,const int n,const bool isupper,al_complex &tau[],double &d[],double &e[]);
   static void       HMatrixTDUnpackQ(CMatrixComplex &a,const int n,const bool isupper,al_complex &tau[],CMatrixComplex &q);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
COrtFac::COrtFac(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
COrtFac::~COrtFac(void)
  {

  }
//+------------------------------------------------------------------+
//| QR decomposition of a rectangular matrix of size MxN             |
//| Input parameters:                                                |
//|     A   -   matrix A whose indexes range within [0..M-1, 0..N-1].|
//|     M   -   number of rows in matrix A.                          |
//|     N   -   number of columns in matrix A.                       |
//| Output parameters:                                               |
//|     A   -   matrices Q and R in compact form (see below).        |
//|     Tau -   array of scalar factors which are used to form       |
//|             matrix Q. Array whose index ranges within            |
//|             [0.. Min(M-1,N-1)].                                  |
//| Matrix A is represented as A = QR, where Q is an orthogonal      |
//| matrix of size MxM, R - upper triangular (or upper trapezoid)    |
//| matrix of size M x N.                                            |
//| The elements of matrix R are located on and above the main       |
//| diagonal of matrix A. The elements which are located in Tau      |
//| array and below the main diagonal of matrix A are used to form   |
//| matrix Q as follows:                                             |
//| Matrix Q is represented as a product of elementary reflections   |
//| Q = H(0)*H(2)*...*H(k-1),                                        |
//| where k = min(m,n), and each H(i) is in the form                 |
//| H(i) = 1 - tau * v * (v^T)                                       |
//| where tau is a scalar stored in Tau[I]; v - real vector,         |
//| so that v(0:i-1) = 0, v(i) = 1, v(i+1:m-1) stored in             |
//| A(i+1:m-1,i).                                                    |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixQR(CMatrixDouble &a,const int m,const int n,double &tau[])
  {
//--- check
   if(m<=0 || n<=0)
      return;
//--- create variables
   int minmn=MathMin(m,n);
   int blockstart=0;
   int blocksize=0;
   int rowscount=0;
   int i=0;
   int i_=0;
   int i1_=0;
//--- create arrays
   double work[];
   double t[];
   double taubuf[];
//--- create matrix
   CMatrixDouble tmpa;
   CMatrixDouble tmpt;
   CMatrixDouble tmpr;
//--- allocation
   ArrayResizeAL(work,MathMax(m,n)+1);
   ArrayResizeAL(t,MathMax(m,n)+1);
   ArrayResizeAL(tau,minmn);
   ArrayResizeAL(taubuf,minmn);
//--- allocation
   tmpa.Resize(m,CAblas::AblasBlockSize());
   tmpt.Resize(CAblas::AblasBlockSize(),2*CAblas::AblasBlockSize());
   tmpr.Resize(2*CAblas::AblasBlockSize(),n);
//--- Blocked code
   while(blockstart!=minmn)
     {
      //--- Determine block size
      blocksize=minmn-blockstart;
      if(blocksize>CAblas::AblasBlockSize())
         blocksize=CAblas::AblasBlockSize();
      //--- change
      rowscount=m-blockstart;
      //--- QR decomposition of submatrix.
      //--- Matrix is copied to temporary storage to solve
      //--- some TLB issues arising from non-contiguous memory
      //--- access pattern.
      CAblas::RMatrixCopy(rowscount,blocksize,a,blockstart,blockstart,tmpa,0,0);
      RMatrixQRBaseCase(tmpa,rowscount,blocksize,work,t,taubuf);
      CAblas::RMatrixCopy(rowscount,blocksize,tmpa,0,0,a,blockstart,blockstart);
      i1_=-blockstart;
      for(i_=blockstart;i_<=blockstart+blocksize-1;i_++)
         tau[i_]=taubuf[i_+i1_];
      //--- check
      if(blockstart+blocksize<=n-1)
        {
         //--- Update the rest, choose between:
         //--- a) Level 2 algorithm (when the rest of the matrix is small enough)
         //--- b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
         //---    representation for products of Householder transformations',
         //---    by R. Schreiber and C. Van Loan.
         if(n-blockstart-blocksize>=2*CAblas::AblasBlockSize() || rowscount>=4*CAblas::AblasBlockSize())
           {
            //--- Prepare block reflector
            RMatrixBlockReflector(tmpa,taubuf,true,rowscount,blocksize,tmpt,work);
            //--- Multiply the rest of A by Q'.
            //--- Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
            //--- Q' = E + Y*T'*Y' = E + TmpA*TmpT'*TmpA'
            CAblas::RMatrixGemm(blocksize,n-blockstart-blocksize,rowscount,1.0,tmpa,0,0,1,a,blockstart,blockstart+blocksize,0,0.0,tmpr,0,0);
            CAblas::RMatrixGemm(blocksize,n-blockstart-blocksize,blocksize,1.0,tmpt,0,0,1,tmpr,0,0,0,0.0,tmpr,blocksize,0);
            CAblas::RMatrixGemm(rowscount,n-blockstart-blocksize,blocksize,1.0,tmpa,0,0,0,tmpr,blocksize,0,0,1.0,a,blockstart,blockstart+blocksize);
           }
         else
           {
            //--- Level 2 algorithm
            for(i=0;i<blocksize;i++)
              {
               i1_=i-1;
               for(i_=1;i_<=rowscount-i;i_++)
                  t[i_]=tmpa[i_+i1_][i];
               t[1]=1;
               //--- function call
               CReflections::ApplyReflectionFromTheLeft(a,taubuf[i],t,blockstart+i,m-1,blockstart+blocksize,n-1,work);
              }
           }
        }
      //--- change value
      blockstart=blockstart+blocksize;
     }
  }
//+------------------------------------------------------------------+
//| LQ decomposition of a rectangular matrix of size MxN             |
//| Input parameters:                                                |
//|     A   -   matrix A whose indexes range within [0..M-1, 0..N-1].|
//|     M   -   number of rows in matrix A.                          |
//|     N   -   number of columns in matrix A.                       |
//| Output parameters:                                               |
//|     A   -   matrices L and Q in compact form (see below)         |
//|     Tau -   array of scalar factors which are used to form       |
//|             matrix Q. Array whose index ranges within            |
//|             [0..Min(M,N)-1].                                     |
//| Matrix A is represented as A = LQ, where Q is an orthogonal      |
//| matrix of size MxM, L - lower triangular (or lower trapezoid)    |
//| matrix of size M x N.                                            |
//| The elements of matrix L are located on and below the main       |
//| diagonal of matrix A. The elements which are located in Tau      |
//| array and above the main diagonal of matrix A are used to form   |
//| matrix Q as follows:                                             |
//| Matrix Q is represented as a product of elementary reflections   |
//| Q = H(k-1)*H(k-2)*...*H(1)*H(0),                                 |
//| where k = min(m,n), and each H(i) is of the form                 |
//| H(i) = 1 - tau * v * (v^T)                                       |
//| where tau is a scalar stored in Tau[I]; v - real vector, so that |
//| v(0:i-1)=0, v(i) = 1, v(i+1:n-1) stored in A(i,i+1:n-1).         |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixLQ(CMatrixDouble &a,const int m,const int n,double &tau[])
  {
//--- check
   if(m<=0 || n<=0)
      return;
//--- create variables
   int minmn=MathMin(m,n);
   int blockstart=0;
   int blocksize=0;
   int columnscount=0;
   int i=0;
   int i_=0;
   int i1_=0;
//--- create arrays
   double work[];
   double t[];
   double taubuf[];
//--- create matrix
   CMatrixDouble tmpa;
   CMatrixDouble tmpt;
   CMatrixDouble tmpr;
//--- allocation
   ArrayResizeAL(work,MathMax(m,n)+1);
   ArrayResizeAL(t,MathMax(m,n)+1);
   ArrayResizeAL(tau,minmn);
   ArrayResizeAL(taubuf,minmn);
//--- allocation
   tmpa.Resize(CAblas::AblasBlockSize(),n);
   tmpt.Resize(CAblas::AblasBlockSize(),2*CAblas::AblasBlockSize());
   tmpr.Resize(m,2*CAblas::AblasBlockSize());
//--- Blocked code
   while(blockstart!=minmn)
     {
      //--- Determine block size
      blocksize=minmn-blockstart;
      if(blocksize>CAblas::AblasBlockSize())
         blocksize=CAblas::AblasBlockSize();
      //--- change
      columnscount=n-blockstart;
      //--- LQ decomposition of submatrix.
      //--- Matrix is copied to temporary storage to solve
      //--- some TLB issues arising from non-contiguous memory
      //--- access pattern.
      CAblas::RMatrixCopy(blocksize,columnscount,a,blockstart,blockstart,tmpa,0,0);
      RMatrixLQBaseCase(tmpa,blocksize,columnscount,work,t,taubuf);
      CAblas::RMatrixCopy(blocksize,columnscount,tmpa,0,0,a,blockstart,blockstart);
      i1_=-blockstart;
      for(i_=blockstart;i_<=blockstart+blocksize-1;i_++)
        {
         tau[i_]=taubuf[i_+i1_];
        }
      //--- Update the rest, choose between:
      //--- a) Level 2 algorithm (when the rest of the matrix is small enough)
      //--- b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
      //---    representation for products of Householder transformations',
      //---    by R. Schreiber and C. Van Loan.
      if(blockstart+blocksize<=m-1)
        {
         //--- check
         if(m-blockstart-blocksize>=2*CAblas::AblasBlockSize())
           {
            //--- prepare
            RMatrixBlockReflector(tmpa,taubuf,false,columnscount,blocksize,tmpt,work);
            //--- Multiply the rest of A by Q.
            //--- Q = E + Y*T*Y' = E + TmpA'*TmpT*TmpA
            CAblas::RMatrixGemm(m-blockstart-blocksize,blocksize,columnscount,1.0,a,blockstart+blocksize,blockstart,0,tmpa,0,0,1,0.0,tmpr,0,0);
            CAblas::RMatrixGemm(m-blockstart-blocksize,blocksize,blocksize,1.0,tmpr,0,0,0,tmpt,0,0,0,0.0,tmpr,0,blocksize);
            CAblas::RMatrixGemm(m-blockstart-blocksize,columnscount,blocksize,1.0,tmpr,0,blocksize,0,tmpa,0,0,0,1.0,a,blockstart+blocksize,blockstart);
           }
         else
           {
            //--- Level 2 algorithm
            for(i=0;i<blocksize;i++)
              {
               i1_=i-1;
               for(i_=1;i_<=columnscount-i;i_++)
                  t[i_]=tmpa[i][i_+i1_];
               t[1]=1;
               //--- function call
               CReflections::ApplyReflectionFromTheRight(a,taubuf[i],t,blockstart+blocksize,m-1,blockstart+i,n-1,work);
              }
           }
        }
      //--- change value
      blockstart=blockstart+blocksize;
     }
  }
//+------------------------------------------------------------------+
//| QR decomposition of a rectangular complex matrix of size MxN     |
//| Input parameters:                                                |
//|     A   -   matrix A whose indexes range within [0..M-1, 0..N-1] |
//|     M   -   number of rows in matrix A.                          |
//|     N   -   number of columns in matrix A.                       |
//| Output parameters:                                               |
//|     A   -   matrices Q and R in compact form                     |
//|     Tau -   array of scalar factors which are used to form       |
//|             matrix Q. Array whose indexes range within           |
//|             [0.. Min(M,N)-1]                                     |
//| Matrix A is represented as A = QR, where Q is an orthogonal      |
//| matrix of size MxM, R - upper triangular (or upper trapezoid)    |
//| matrix of size MxN.                                              |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      September 30, 1994                                          |
//+------------------------------------------------------------------+
static void COrtFac::CMatrixQR(CMatrixComplex &a,const int m,const int n,al_complex &tau[])
  {
//--- check
   if(m<=0 || n<=0)
      return;
//--- create arrays
   al_complex work[];
   al_complex t[];
   al_complex taubuf[];
//--- create matrix
   CMatrixComplex tmpa;
   CMatrixComplex tmpt;
   CMatrixComplex tmpr;
//--- create variables
   int        minmn=MathMin(m,n);;
   int        blockstart=0;
   int        blocksize=0;
   int        rowscount=0;
   int        i=0;
   int        i_=0;
   int        i1_=0;
   al_complex One(1,0);
   al_complex Alpha(1,0);
   al_complex Beta(0,0);
//--- allocation
   ArrayResizeAL(work,MathMax(m,n)+1);
   ArrayResizeAL(t,MathMax(m,n)+1);
   ArrayResizeAL(tau,minmn);
   ArrayResizeAL(taubuf,minmn);
//--- allocation
   tmpa.Resize(m,CAblas::AblasComplexBlockSize());
   tmpt.Resize(CAblas::AblasComplexBlockSize(),CAblas::AblasComplexBlockSize());
   tmpr.Resize(2*CAblas::AblasComplexBlockSize(),n);
//--- Blocked code
   while(blockstart!=minmn)
     {
      //--- Determine block size
      blocksize=minmn-blockstart;
      if(blocksize>CAblas::AblasComplexBlockSize())
         blocksize=CAblas::AblasComplexBlockSize();
      rowscount=m-blockstart;
      //--- QR decomposition of submatrix.
      //--- Matrix is copied to temporary storage to solve
      //--- some TLB issues arising from non-contiguous memory
      //--- access pattern.
      CAblas::CMatrixCopy(rowscount,blocksize,a,blockstart,blockstart,tmpa,0,0);
      CMatrixQRBaseCase(tmpa,rowscount,blocksize,work,t,taubuf);
      CAblas::CMatrixCopy(rowscount,blocksize,tmpa,0,0,a,blockstart,blockstart);
      i1_=-blockstart;
      for(i_=blockstart;i_<=blockstart+blocksize-1;i_++)
         tau[i_]=taubuf[i_+i1_];
      //--- Update the rest, choose between:
      //--- a) Level 2 algorithm (when the rest of the matrix is small enough)
      //--- b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
      //---    representation for products of Householder transformations',
      //---    by R. Schreiber and C. Van Loan.
      if(blockstart+blocksize<=n-1)
        {
         //--- check
         if(n-blockstart-blocksize>=2*CAblas::AblasComplexBlockSize())
           {
            //--- Prepare block reflector
            CMatrixBlockReflector(tmpa,taubuf,true,rowscount,blocksize,tmpt,work);
            //--- Multiply the rest of A by Q'.
            //--- Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
            //--- Q' = E + Y*T'*Y' = E + TmpA*TmpT'*TmpA'
            CAblas::CMatrixGemm(blocksize,n-blockstart-blocksize,rowscount,Alpha,tmpa,0,0,2,a,blockstart,blockstart+blocksize,0,Beta,tmpr,0,0);
            CAblas::CMatrixGemm(blocksize,n-blockstart-blocksize,blocksize,Alpha,tmpt,0,0,2,tmpr,0,0,0,Beta,tmpr,blocksize,0);
            CAblas::CMatrixGemm(rowscount,n-blockstart-blocksize,blocksize,Alpha,tmpa,0,0,0,tmpr,blocksize,0,0,Alpha,a,blockstart,blockstart+blocksize);
           }
         else
           {
            //--- Level 2 algorithm
            for(i=0;i<blocksize;i++)
              {
               i1_=i-1;
               for(i_=1;i_<=rowscount-i;i_++)
                  t[i_]=tmpa[i_+i1_][i];
               t[1]=One;
               //--- function call
               CComplexReflections::ComplexApplyReflectionFromTheLeft(a,CMath::Conj(taubuf[i]),t,blockstart+i,m-1,blockstart+blocksize,n-1,work);
              }
           }
        }
      //--- change value
      blockstart=blockstart+blocksize;
     }
  }
//+------------------------------------------------------------------+
//| LQ decomposition of a rectangular complex matrix of size MxN     |
//| Input parameters:                                                |
//|     A   -   matrix A whose indexes range within [0..M-1, 0..N-1] |
//|     M   -   number of rows in matrix A.                          |
//|     N   -   number of columns in matrix A.                       |
//| Output parameters:                                               |
//|     A   -   matrices Q and L in compact form                     |
//|     Tau -   array of scalar factors which are used to form       |
//|             matrix Q. Array whose indexes range within           |
//|             [0.. Min(M,N)-1]                                     |
//| Matrix A is represented as A = LQ, where Q is an orthogonal      |
//| matrix of size MxM, L - lower triangular (or lower trapezoid)    |
//| matrix of size MxN.                                              |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      September 30, 1994                                          |
//+------------------------------------------------------------------+
static void COrtFac::CMatrixLQ(CMatrixComplex &a,const int m,const int n,al_complex &tau[])
  {
//--- check
   if(m<=0 || n<=0)
      return;
//--- create arrays
   al_complex work[];
   al_complex t[];
   al_complex taubuf[];
//--- create matrix
   CMatrixComplex tmpa;
   CMatrixComplex tmpt;
   CMatrixComplex tmpr;
//--- create variables
   int        minmn=MathMin(m,n);
   int        blockstart=0;
   int        blocksize=0;
   int        columnscount=0;
   int        i=0;
   int        i_=0;
   int        i1_=0;
   al_complex One(1,0);
   al_complex Alpha(1,0);
   al_complex Beta(0,0);
//--- allocation
   ArrayResizeAL(work,MathMax(m,n)+1);
   ArrayResizeAL(t,MathMax(m,n)+1);
   ArrayResizeAL(tau,minmn);
   ArrayResizeAL(taubuf,minmn);
//--- allocation
   tmpa.Resize(CAblas::AblasComplexBlockSize(),n);
   tmpt.Resize(CAblas::AblasComplexBlockSize(),CAblas::AblasComplexBlockSize());
   tmpr.Resize(m,2*CAblas::AblasComplexBlockSize());
//--- Blocked code
   while(blockstart!=minmn)
     {
      //--- Determine block size
      blocksize=minmn-blockstart;
      if(blocksize>CAblas::AblasComplexBlockSize())
         blocksize=CAblas::AblasComplexBlockSize();
      columnscount=n-blockstart;
      //--- LQ decomposition of submatrix.
      //--- Matrix is copied to temporary storage to solve
      //--- some TLB issues arising from non-contiguous memory
      //--- access pattern.
      CAblas::CMatrixCopy(blocksize,columnscount,a,blockstart,blockstart,tmpa,0,0);
      CMatrixLQBaseCase(tmpa,blocksize,columnscount,work,t,taubuf);
      CAblas::CMatrixCopy(blocksize,columnscount,tmpa,0,0,a,blockstart,blockstart);
      i1_=-blockstart;
      for(i_=blockstart;i_<=blockstart+blocksize-1;i_++)
         tau[i_]=taubuf[i_+i1_];
      //--- Update the rest, choose between:
      //--- a) Level 2 algorithm (when the rest of the matrix is small enough)
      //--- b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
      //---    representation for products of Householder transformations',
      //---    by R. Schreiber and C. Van Loan.
      if(blockstart+blocksize<=m-1)
        {
         //--- check
         if(m-blockstart-blocksize>=2*CAblas::AblasComplexBlockSize())
           {
            //--- Prepare block reflector
            CMatrixBlockReflector(tmpa,taubuf,false,columnscount,blocksize,tmpt,work);
            //--- Multiply the rest of A by Q.
            //--- Q  = E + Y*T*Y'  = E + TmpA'*TmpT*TmpA
            CAblas::CMatrixGemm(m-blockstart-blocksize,blocksize,columnscount,Alpha,a,blockstart+blocksize,blockstart,0,tmpa,0,0,2,Beta,tmpr,0,0);
            CAblas::CMatrixGemm(m-blockstart-blocksize,blocksize,blocksize,Alpha,tmpr,0,0,0,tmpt,0,0,0,Beta,tmpr,0,blocksize);
            CAblas::CMatrixGemm(m-blockstart-blocksize,columnscount,blocksize,Alpha,tmpr,0,blocksize,0,tmpa,0,0,0,Alpha,a,blockstart+blocksize,blockstart);
           }
         else
           {
            //--- Level 2 algorithm
            for(i=0;i<blocksize;i++)
              {
               i1_=i-1;
               for(i_=1;i_<=columnscount-i;i_++)
                  t[i_]=CMath::Conj(tmpa[i][i_+i1_]);
               t[1]=One;
               //--- function call
               CComplexReflections::ComplexApplyReflectionFromTheRight(a,taubuf[i],t,blockstart+blocksize,m-1,blockstart+i,n-1,work);
              }
           }
        }
      //--- change value
      blockstart=blockstart+blocksize;
     }
  }
//+------------------------------------------------------------------+
//| Partial unpacking of matrix Q from the QR decomposition of a     |
//| matrix A                                                         |
//| Input parameters:                                                |
//|     A       -   matrices Q and R in compact form.                |
//|                 Output of RMatrixQR subroutine.                  |
//|     M       -   number of rows in given matrix A. M>=0.          |
//|     N       -   number of columns in given matrix A. N>=0.       |
//|     Tau     -   scalar factors which are used to form Q.         |
//|                 Output of the RMatrixQR subroutine.              |
//|     QColumns -  required number of columns of matrix Q.          |
//|                 M>=QColumns>=0.                                  |
//| Output parameters:                                               |
//|     Q       -   first QColumns columns of matrix Q.              |
//|                 Array whose indexes range within                 |
//|                 [0..M-1, 0..QColumns-1].                         |
//|                 If QColumns=0, the array remains unchanged.      |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixQRUnpackQ(CMatrixDouble &a,const int m,const int n,
                                      double &tau[],const int qcolumns,CMatrixDouble &q)
  {
//--- check
   if(!CAp::Assert(qcolumns<=m,__FUNCTION__+": QColumns>M!"))
      return;
//--- check
   if(m<=0 || n<=0 || qcolumns<=0)
      return;
//--- create arrays
   double work[];
   double t[];
   double taubuf[];
//--- create matrix
   CMatrixDouble tmpa;
   CMatrixDouble tmpt;
   CMatrixDouble tmpr;
//--- create variables
   int minmn=MathMin(m,n);
   int refcnt=MathMin(minmn,qcolumns);
   int blockstart=CAblas::AblasBlockSize()*(refcnt/CAblas::AblasBlockSize());
   int blocksize=refcnt-blockstart;
   int rowscount=0;
   int i=0;
   int j=0;
   int i_=0;
   int i1_=0;
//--- allocation
   q.Resize(m,qcolumns);
//--- identity matrix
   for(i=0;i<m;i++)
     {
      for(j=0;j<qcolumns;j++)
        {
         if(i==j)
            q[i].Set(j,1);
         else
            q[i].Set(j,0);
        }
     }
//--- allocation
   ArrayResizeAL(work,MathMax(m,qcolumns)+1);
   ArrayResizeAL(t,MathMax(m,qcolumns)+1);
   ArrayResizeAL(taubuf,minmn);
//--- allocation
   tmpa.Resize(m,CAblas::AblasBlockSize());
   tmpt.Resize(CAblas::AblasBlockSize(),2*CAblas::AblasBlockSize());
   tmpr.Resize(2*CAblas::AblasBlockSize(),qcolumns);
//--- Blocked code
   while(blockstart>=0)
     {
      rowscount=m-blockstart;
      //--- Copy current block
      CAblas::RMatrixCopy(rowscount,blocksize,a,blockstart,blockstart,tmpa,0,0);
      i1_=blockstart;
      for(i_=0;i_<blocksize;i_++)
         taubuf[i_]=tau[i_+i1_];
      //--- Update, choose between:
      //--- a) Level 2 algorithm (when the rest of the matrix is small enough)
      //--- b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
      //---    representation for products of Householder transformations',
      //---    by R. Schreiber and C. Van Loan.
      if(qcolumns>=2*CAblas::AblasBlockSize())
        {
         //--- Prepare block reflector
         RMatrixBlockReflector(tmpa,taubuf,true,rowscount,blocksize,tmpt,work);
         //--- Multiply matrix by Q.
         //--- Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
         CAblas::RMatrixGemm(blocksize,qcolumns,rowscount,1.0,tmpa,0,0,1,q,blockstart,0,0,0.0,tmpr,0,0);
         CAblas::RMatrixGemm(blocksize,qcolumns,blocksize,1.0,tmpt,0,0,0,tmpr,0,0,0,0.0,tmpr,blocksize,0);
         CAblas::RMatrixGemm(rowscount,qcolumns,blocksize,1.0,tmpa,0,0,0,tmpr,blocksize,0,0,1.0,q,blockstart,0);
        }
      else
        {
         //--- Level 2 algorithm
         for(i=blocksize-1;i>=0;i--)
           {
            i1_=i-1;
            for(i_=1;i_<=rowscount-i;i_++)
              {
               t[i_]=tmpa[i_+i1_][i];
              }
            t[1]=1;
            //--- function call
            CReflections::ApplyReflectionFromTheLeft(q,taubuf[i],t,blockstart+i,m-1,0,qcolumns-1,work);
           }
        }
      //--- change value
      blockstart=blockstart-CAblas::AblasBlockSize();
      blocksize=CAblas::AblasBlockSize();
     }
  }
//+------------------------------------------------------------------+
//| Unpacking of matrix R from the QR decomposition of a matrix A    |
//| Input parameters:                                                |
//|     A       -   matrices Q and R in compact form.                |
//|                 Output of RMatrixQR subroutine.                  |
//|     M       -   number of rows in given matrix A. M>=0.          |
//|     N       -   number of columns in given matrix A. N>=0.       |
//| Output parameters:                                               |
//|     R       -   matrix R, array[0..M-1, 0..N-1].                 |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixQRUnpackR(CMatrixDouble &a,const int m,const int n,CMatrixDouble &r)
  {
//--- check
   if(m<=0 || n<=0)
      return;
//--- create variables
   int i=0;
   int k=MathMin(m,n);
   int i_=0;
//--- allocation
   r.Resize(m,n);
//--- Prepare matrix
   for(i=0;i<n;i++)
      r[0].Set(i,0);
   for(i=1;i<m;i++)
     {
      for(i_=0;i_<n;i_++)
         r[i].Set(i_,r[0][i_]);
     }
//--- get result
   for(i=0;i<k;i++)
     {
      for(i_=i;i_<n;i_++)
         r[i].Set(i_,a[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| Partial unpacking of matrix Q from the LQ decomposition of a     |
//| matrix A                                                         |
//| Input parameters:                                                |
//|     A       -   matrices L and Q in compact form.                |
//|                 Output of RMatrixLQ subroutine.                  |
//|     M       -   number of rows in given matrix A. M>=0.          |
//|     N       -   number of columns in given matrix A. N>=0.       |
//|     Tau     -   scalar factors which are used to form Q.         |
//|                 Output of the RMatrixLQ subroutine.              |
//|     QRows   -   required number of rows in matrix Q. N>=QRows>=0.|
//| Output parameters:                                               |
//|     Q       -   first QRows rows of matrix Q. Array whose indexes|
//|                 range within [0..QRows-1, 0..N-1]. If QRows=0,   |
//|                 the array remains unchanged.                     |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixLQUnpackQ(CMatrixDouble &a,const int m,const int n,
                                      double &tau[],const int qrows,CMatrixDouble &q)
  {
//--- check
   if(!CAp::Assert(qrows<=n,__FUNCTION__+": QRows>N!"))
      return;
//--- check
   if(m<=0 || n<=0 || qrows<=0)
      return;
//--- create arrays
   double work[];
   double t[];
   double taubuf[];
//--- create matrix
   CMatrixDouble tmpa;
   CMatrixDouble tmpt;
   CMatrixDouble tmpr;
//--- create variables
   int minmn=MathMin(m,n);
   int refcnt=MathMin(minmn,qrows);
   int blockstart=CAblas::AblasBlockSize()*(refcnt/CAblas::AblasBlockSize());
   int blocksize=refcnt-blockstart;
   int columnscount=0;
   int i=0;
   int j=0;
   int i_=0;
   int i1_=0;
//--- allocation
   ArrayResizeAL(work,MathMax(m,n)+1);
   ArrayResizeAL(t,MathMax(m,n)+1);
   ArrayResizeAL(taubuf,minmn);
//--- allocation
   tmpa.Resize(CAblas::AblasBlockSize(),n);
   tmpt.Resize(CAblas::AblasBlockSize(),2*CAblas::AblasBlockSize());
   tmpr.Resize(qrows,2*CAblas::AblasBlockSize());
   q.Resize(qrows,n);
//--- identity matrix
   for(i=0;i<=qrows-1;i++)
     {
      for(j=0;j<n;j++)
        {
         if(i==j)
            q[i].Set(j,1);
         else
            q[i].Set(j,0);
        }
     }
//--- Blocked code
   while(blockstart>=0)
     {
      columnscount=n-blockstart;
      //--- Copy submatrix
      CAblas::RMatrixCopy(blocksize,columnscount,a,blockstart,blockstart,tmpa,0,0);
      i1_=blockstart;
      for(i_=0;i_<blocksize;i_++)
         taubuf[i_]=tau[i_+i1_];
      //--- Update matrix, choose between:
      //--- a) Level 2 algorithm (when the rest of the matrix is small enough)
      //--- b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
      //---    representation for products of Householder transformations',
      //---    by R. Schreiber and C. Van Loan.
      if(qrows>=2*CAblas::AblasBlockSize())
        {
         //--- Prepare block reflector
         RMatrixBlockReflector(tmpa,taubuf,false,columnscount,blocksize,tmpt,work);
         //--- Multiply the rest of A by Q'.
         //--- Q'  = E + Y*T'*Y'  = E + TmpA'*TmpT'*TmpA
         CAblas::RMatrixGemm(qrows,blocksize,columnscount,1.0,q,0,blockstart,0,tmpa,0,0,1,0.0,tmpr,0,0);
         CAblas::RMatrixGemm(qrows,blocksize,blocksize,1.0,tmpr,0,0,0,tmpt,0,0,1,0.0,tmpr,0,blocksize);
         CAblas::RMatrixGemm(qrows,columnscount,blocksize,1.0,tmpr,0,blocksize,0,tmpa,0,0,0,1.0,q,0,blockstart);
        }
      else
        {
         //--- Level 2 algorithm
         for(i=blocksize-1;i>=0;i--)
           {
            i1_=i-1;
            for(i_=1;i_<=columnscount-i;i_++)
               t[i_]=tmpa[i][i_+i1_];
            t[1]=1;
            //--- function call
            CReflections::ApplyReflectionFromTheRight(q,taubuf[i],t,0,qrows-1,blockstart+i,n-1,work);
           }
        }
      //--- change value
      blockstart=blockstart-CAblas::AblasBlockSize();
      blocksize=CAblas::AblasBlockSize();
     }
  }
//+------------------------------------------------------------------+
//| Unpacking of matrix L from the LQ decomposition of a matrix A    |
//| Input parameters:                                                |
//|     A    -matrices Q and L in compact form.                      |
//|                 Output of RMatrixLQ subroutine.                  |
//|     M    -number of rows in given matrix A. M>=0.                |
//|     N    -number of columns in given matrix A. N>=0.             |
//| Output parameters:                                               |
//|     L    -matrix L, array[0..M-1,0..N-1].                        |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixLQUnpackL(CMatrixDouble &a,const int m,const int n,CMatrixDouble &l)
  {
//--- check
   if(m<=0 || n<=0)
      return;
//--- create variables
   int i=0;
   int k=0;
   int i_=0;
//--- allocation
   l.Resize(m,n);
//--- Prepare matrix
   for(i=0;i<n;i++)
      l[0].Set(i,0);
   for(i=1;i<m;i++)
     {
      for(i_=0;i_<n;i_++)
         l[i].Set(i_,l[0][i_]);
     }
//--- get result
   for(i=0;i<m;i++)
     {
      k=MathMin(i,n-1);
      for(i_=0;i_<=k;i_++)
         l[i].Set(i_,a[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| Partial unpacking of matrix Q from QR decomposition of a complex |
//| matrix A.                                                        |
//| Input parameters:                                                |
//|     A           -   matrices Q and R in compact form.            |
//|                     Output of CMatrixQR subroutine .             |
//|     M           -   number of rows in matrix A. M>=0.            |
//|     N           -   number of columns in matrix A. N>=0.         |
//|     Tau         -   scalar factors which are used to form Q.     |
//|                     Output of CMatrixQR subroutine .             |
//|     QColumns    -   required number of columns in matrix Q.      |
//|                     M>=QColumns>=0.                              |
//| Output parameters:                                               |
//|     Q           -   first QColumns columns of matrix Q.          |
//|                     Array whose index ranges within [0..M-1,     |
//|                     0..QColumns-1].                              |
//|                     If QColumns=0, array isn't changed.          |
//+------------------------------------------------------------------+
static void COrtFac::CMatrixQRUnpackQ(CMatrixComplex &a,const int m,const int n,
                                      al_complex &tau[],const int qcolumns,CMatrixComplex &q)
  {
//--- check
   if(!CAp::Assert(qcolumns<=m,__FUNCTION__+": QColumns>M!"))
      return;
//--- check
   if(m<=0 || n<=0)
      return;
//--- create arrays
   al_complex work[];
   al_complex t[];
   al_complex taubuf[];
//--- create matrix
   CMatrixComplex tmpa;
   CMatrixComplex tmpt;
   CMatrixComplex tmpr;
//--- create variables
   int        minmn=MathMin(m,n);
   int        refcnt=MathMin(minmn,qcolumns);
   int        blockstart=CAblas::AblasComplexBlockSize()*(refcnt/CAblas::AblasComplexBlockSize());
   int        blocksize=refcnt-blockstart;
   int        rowscount=0;
   int        i=0;
   int        j=0;
   int        i_=0;
   int        i1_=0;
   al_complex One(1,0);
   al_complex Zero(0,0);
//--- allocation
   ArrayResizeAL(work,MathMax(m,n)+1);
   ArrayResizeAL(t,MathMax(m,n)+1);
   ArrayResizeAL(taubuf,minmn);
//--- allocation
   tmpa.Resize(m,CAblas::AblasComplexBlockSize());
   tmpt.Resize(CAblas::AblasComplexBlockSize(),CAblas::AblasComplexBlockSize());
   tmpr.Resize(2*CAblas::AblasComplexBlockSize(),qcolumns);
   q.Resize(m,qcolumns);
   for(i=0;i<m;i++)
     {
      for(j=0;j<qcolumns;j++)
        {
         if(i==j)
            q[i].Set(j,One);
         else
            q[i].Set(j,Zero);
        }
     }
//--- Blocked code
   while(blockstart>=0)
     {
      rowscount=m-blockstart;
      //--- QR decomposition of submatrix.
      //--- Matrix is copied to temporary storage to solve
      //--- some TLB issues arising from non-contiguous memory
      //--- access pattern.
      CAblas::CMatrixCopy(rowscount,blocksize,a,blockstart,blockstart,tmpa,0,0);
      i1_=blockstart;
      for(i_=0;i_<blocksize;i_++)
         taubuf[i_]=tau[i_+i1_];
      //--- Update matrix, choose between:
      //--- a) Level 2 algorithm (when the rest of the matrix is small enough)
      //--- b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
      //---    representation for products of Householder transformations',
      //---    by R. Schreiber and C. Van Loan.
      if(qcolumns>=2*CAblas::AblasComplexBlockSize())
        {
         //--- Prepare block reflector
         CMatrixBlockReflector(tmpa,taubuf,true,rowscount,blocksize,tmpt,work);
         //--- Multiply the rest of A by Q.
         //--- Q  = E + Y*T*Y'  = E + TmpA*TmpT*TmpA'
         CAblas::CMatrixGemm(blocksize,qcolumns,rowscount,One,tmpa,0,0,2,q,blockstart,0,0,Zero,tmpr,0,0);
         CAblas::CMatrixGemm(blocksize,qcolumns,blocksize,One,tmpt,0,0,0,tmpr,0,0,0,Zero,tmpr,blocksize,0);
         CAblas::CMatrixGemm(rowscount,qcolumns,blocksize,One,tmpa,0,0,0,tmpr,blocksize,0,0,One,q,blockstart,0);
        }
      else
        {
         //--- Level 2 algorithm
         for(i=blocksize-1;i>=0;i--)
           {
            i1_=i-1;
            for(i_=1;i_<=rowscount-i;i_++)
               t[i_]=tmpa[i_+i1_][i];
            t[1]=1;
            //--- function call
            CComplexReflections::ComplexApplyReflectionFromTheLeft(q,taubuf[i],t,blockstart+i,m-1,0,qcolumns-1,work);
           }
        }
      //--- change value
      blockstart=blockstart-CAblas::AblasComplexBlockSize();
      blocksize=CAblas::AblasComplexBlockSize();
     }
  }
//+------------------------------------------------------------------+
//| Unpacking of matrix R from the QR decomposition of a matrix A    |
//| Input parameters:                                                |
//|     A       -   matrices Q and R in compact form.                |
//|                 Output of CMatrixQR subroutine.                  |
//|     M       -   number of rows in given matrix A. M>=0.          |
//|     N       -   number of columns in given matrix A. N>=0.       |
//| Output parameters:                                               |
//|     R       -   matrix R, array[0..M-1, 0..N-1].                 |
//+------------------------------------------------------------------+
static void COrtFac::CMatrixQRUnpackR(CMatrixComplex &a,const int m,const int n,CMatrixComplex &r)
  {
//--- check
   if(m<=0 || n<=0)
      return;
//--- create variables
   al_complex Zero(0,0);
   int        i=0;
   int        k=MathMin(m,n);
   int        i_=0;
//--- allocation
   r.Resize(m,n);
//--- Prepare matrix
   for(i=0;i<n;i++)
      r[0].Set(i,Zero);
   for(i=1;i<m;i++)
     {
      for(i_=0;i_<n;i_++)
         r[i].Set(i_,r[0][i_]);
     }
//--- get result
   for(i=0;i<k;i++)
     {
      for(i_=i;i_<n;i_++)
         r[i].Set(i_,a[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| Partial unpacking of matrix Q from LQ decomposition of a complex |
//| matrix A.                                                        |
//| Input parameters:                                                |
//|     A           -   matrices Q and R in compact form.            |
//|                     Output of CMatrixLQ subroutine.              |
//|     M           -   number of rows in matrix A. M>=0.            |
//|     N           -   number of columns in matrix A. N>=0.         |
//|     Tau         -   scalar factors which are used to form Q.     |
//|                     Output of CMatrixLQ subroutine .             |
//|     QRows       -   required number of rows in matrix Q.         |
//|                     N>=QColumns>=0.                              |
//| Output parameters:                                               |
//|     Q           -   first QRows rows of matrix Q.                |
//|                     Array whose index ranges within [0..QRows-1, |
//|                     0..N-1].                                     |
//|                     If QRows=0, array isn't changed.             |
//+------------------------------------------------------------------+
static void COrtFac::CMatrixLQUnpackQ(CMatrixComplex &a,const int m,const int n,
                                      al_complex &tau[],const int qrows,CMatrixComplex &q)
  {
//--- check
   if(m<=0 || n<=0)
      return;
//--- create arrays
   al_complex work[];
   al_complex t[];
   al_complex taubuf[];
//--- create matrix
   CMatrixComplex tmpa;
   CMatrixComplex tmpt;
   CMatrixComplex tmpr;
//--- create variables
   int        minmn=MathMin(m,n);
   int        refcnt=MathMin(minmn,qrows);
   int        blockstart=CAblas::AblasComplexBlockSize()*(refcnt/CAblas::AblasComplexBlockSize());
   int        blocksize=refcnt-blockstart;
   int        columnscount=0;
   int        i=0;
   int        j=0;
   int        i_=0;
   int        i1_=0;
   al_complex One(1,0);
   al_complex Zero(0,0);
//--- allocation
   ArrayResizeAL(work,MathMax(m,n)+1);
   ArrayResizeAL(t,MathMax(m,n)+1);
   ArrayResizeAL(taubuf,minmn);
//--- allocation
   tmpa.Resize(CAblas::AblasComplexBlockSize(),n);
   tmpt.Resize(CAblas::AblasComplexBlockSize(),CAblas::AblasComplexBlockSize());
   tmpr.Resize(qrows,2*CAblas::AblasComplexBlockSize());
   q.Resize(qrows,n);
   for(i=0;i<=qrows-1;i++)
     {
      for(j=0;j<n;j++)
        {
         if(i==j)
            q[i].Set(j,One);
         else
            q[i].Set(j,Zero);
        }
     }
//--- Blocked code
   while(blockstart>=0)
     {
      columnscount=n-blockstart;
      //--- LQ decomposition of submatrix.
      //--- Matrix is copied to temporary storage to solve
      //--- some TLB issues arising from non-contiguous memory
      //--- access pattern.
      CAblas::CMatrixCopy(blocksize,columnscount,a,blockstart,blockstart,tmpa,0,0);
      i1_=blockstart;
      for(i_=0;i_<blocksize;i_++)
         taubuf[i_]=tau[i_+i1_];
      //--- Update matrix, choose between:
      //--- a) Level 2 algorithm (when the rest of the matrix is small enough)
      //--- b) blocked algorithm, see algorithm 5 from  'A storage efficient WY
      //---    representation for products of Householder transformations',
      //---    by R. Schreiber and C. Van Loan.
      if(qrows>=2*CAblas::AblasComplexBlockSize())
        {
         //--- Prepare block reflector
         CMatrixBlockReflector(tmpa,taubuf,false,columnscount,blocksize,tmpt,work);
         //--- Multiply the rest of A by Q'.
         //--- Q'  = E + Y*T'*Y'  = E + TmpA'*TmpT'*TmpA
         CAblas::CMatrixGemm(qrows,blocksize,columnscount,One,q,0,blockstart,0,tmpa,0,0,2,Zero,tmpr,0,0);
         CAblas::CMatrixGemm(qrows,blocksize,blocksize,One,tmpr,0,0,0,tmpt,0,0,2,Zero,tmpr,0,blocksize);
         CAblas::CMatrixGemm(qrows,columnscount,blocksize,One,tmpr,0,blocksize,0,tmpa,0,0,0,One,q,0,blockstart);
        }
      else
        {
         //--- Level 2 algorithm
         for(i=blocksize-1;i>=0;i--)
           {
            i1_=i-1;
            for(i_=1;i_<=columnscount-i;i_++)
               t[i_]=CMath::Conj(tmpa[i][i_+i1_]);
            t[1]=1;
            //--- function call
            CComplexReflections::ComplexApplyReflectionFromTheRight(q,CMath::Conj(taubuf[i]),t,0,qrows-1,blockstart+i,n-1,work);
           }
        }
      //--- change value
      blockstart=blockstart-CAblas::AblasComplexBlockSize();
      blocksize=CAblas::AblasComplexBlockSize();
     }
  }
//+------------------------------------------------------------------+
//| Unpacking of matrix L from the LQ decomposition of a matrix A    |
//| Input parameters:                                                |
//|     A       -   matrices Q and L in compact form.                |
//|                 Output of CMatrixLQ subroutine.                  |
//|     M       -   number of rows in given matrix A. M>=0.          |
//|     N       -   number of columns in given matrix A. N>=0.       |
//| Output parameters:                                               |
//|     L       -   matrix L, array[0..M-1, 0..N-1].                 |
//+------------------------------------------------------------------+
static void COrtFac::CMatrixLQUnpackL(CMatrixComplex &a,const int m,const int n,CMatrixComplex &l)
  {
//--- check
   if(m<=0 || n<=0)
      return;
//--- create variables
   al_complex Zero(0,0);
   int        i=0;
   int        k=0;
   int        i_=0;
//--- allocation
   l.Resize(m,n);
//--- Prepare matrix
   for(i=0;i<n;i++)
      l[0].Set(i,Zero);
   for(i=1;i<m;i++)
     {
      for(i_=0;i_<n;i_++)
         l[i].Set(i_,l[0][i_]);
     }
//--- get result
   for(i=0;i<m;i++)
     {
      k=MathMin(i,n-1);
      for(i_=0;i_<=k;i_++)
         l[i].Set(i_,a[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| Reduction of a rectangular matrix to  bidiagonal form            |
//| The algorithm reduces the rectangular matrix A to  bidiagonal    |
//| form by orthogonal transformations P and Q: A = Q*B*P.           |
//| Input parameters:                                                |
//|     A       -   source matrix. array[0..M-1, 0..N-1]             |
//|     M       -   number of rows in matrix A.                      |
//|     N       -   number of columns in matrix A.                   |
//| Output parameters:                                               |
//|     A       -   matrices Q, B, P in compact form (see below).    |
//|     TauQ    -   scalar factors which are used to form matrix Q.  |
//|     TauP    -   scalar factors which are used to form matrix P.  |
//| The main diagonal and one of the secondary diagonals of matrix A |
//| are replaced with bidiagonal matrix B. Other elements contain    |
//| elementary reflections which form MxM matrix Q and NxN matrix P, |
//| respectively.                                                    |
//| If M>=N, B is the upper bidiagonal MxN matrix and is stored in   |
//| the corresponding elements of matrix A. Matrix Q is represented  |
//| as a product of elementary reflections Q = H(0)*H(1)*...*H(n-1), |
//| where H(i) = 1-tau*v*v'. Here tau is a scalar which is stored in |
//| TauQ[i], and vector v has the following structure: v(0:i-1)=0,   |
//| v(i)=1, v(i+1:m-1) is stored in elements A(i+1:m-1,i).Matrix P is|
//| as follows: P = G(0)*G(1)*...*G(n-2), where G(i) = 1 - tau*u*u'. |
//| Tau is stored in TauP[i], u(0:i)=0, u(i+1)=1, u(i+2:n-1) is      |
//| stored in elements A(i,i+2:n-1).                                 |
//| If M<N, B is the lower bidiagonal MxN matrix and is stored in the|
//| corresponding elements of matrix A. Q = H(0)*H(1)*...*H(m-2),    |
//| where H(i) = 1 - tau*v*v', tau is stored in TauQ, v(0:i)=0,      |
//| v(i+1)=1, v(i+2:m-1) is stored in elements A(i+2:m-1,i).         |
//| P = G(0)*G(1)*...*G(m-1), G(i) = 1-tau*u*u', tau is stored in    |
//| TauP, u(0:i-1)=0, u(i)=1, u(i+1:n-1) is stored in A(i,i+1:n-1).  |
//| EXAMPLE:                                                         |
//| m=6, n=5 (m > n):               m=5, n=6 (m < n):                |
//| (  d   e   u1  u1  u1 )         (  d   u1  u1  u1  u1  u1 )      |
//| (  v1  d   e   u2  u2 )         (  e   d   u2  u2  u2  u2 )      |
//| (  v1  v2  d   e   u3 )         (  v1  e   d   u3  u3  u3 )      |
//| (  v1  v2  v3  d   e  )         (  v1  v2  e   d   u4  u4 )      |
//| (  v1  v2  v3  v4  d  )         (  v1  v2  v3  e   d   u5 )      |
//| (  v1  v2  v3  v4  v5 )                                          |
//| Here vi and ui are vectors which form H(i) and G(i), and d and   |
//| e - are the diagonal and off-diagonal elements of matrix B.      |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixBD(CMatrixDouble &a,const int m,const int n,double &tauq[],double &taup[])
  {
//--- check
   if(n<=0 || m<=0)
      return;
//--- create arrays
   double work[];
   double t[];
//--- create variables
   int    minmn=0;
   int    maxmn=MathMax(m,n);
   int    i=0;
   double ltau=0;
   int    i_=0;
   int    i1_=0;
//--- allocation
   ArrayResizeAL(work,maxmn+1);
   ArrayResizeAL(t,maxmn+1);
//--- check
   if(m>=n)
     {
      ArrayResizeAL(tauq,n);
      ArrayResizeAL(taup,n);
     }
   else
     {
      ArrayResizeAL(tauq,m);
      ArrayResizeAL(taup,m);
     }
//--- check
   if(m>=n)
     {
      //--- Reduce to upper bidiagonal form
      for(i=0;i<n;i++)
        {
         //--- Generate elementary reflector H(i) to annihilate A(i+1:m-1,i)
         i1_=i-1;
         for(i_=1;i_<=m-i;i_++)
            t[i_]=a[i_+i1_][i];
         CReflections::GenerateReflection(t,m-i,ltau);
         tauq[i]=ltau;
         i1_=1-i;
         for(i_=i;i_<m;i_++)
            a[i_].Set(i,t[i_+i1_]);
         t[1]=1;
         //--- Apply H(i) to A(i:m-1,i+1:n-1) from the left
         CReflections::ApplyReflectionFromTheLeft(a,ltau,t,i,m-1,i+1,n-1,work);
         //--- check
         if(i<n-1)
           {
            //--- Generate elementary reflector G(i) to annihilate
            //--- A(i,i+2:n-1)
            i1_=i;
            for(i_=1;i_<n-i;i_++)
               t[i_]=a[i][i_+i1_];
            CReflections::GenerateReflection(t,n-1-i,ltau);
            taup[i]=ltau;
            i1_=-i;
            for(i_=i+1;i_<n;i_++)
               a[i].Set(i_,t[i_+i1_]);
            t[1]=1;
            //--- Apply G(i) to A(i+1:m-1,i+1:n-1) from the right
            CReflections::ApplyReflectionFromTheRight(a,ltau,t,i+1,m-1,i+1,n-1,work);
           }
         else
            taup[i]=0;
        }
     }
   else
     {
      //--- Reduce to lower bidiagonal form
      for(i=0;i<m;i++)
        {
         //--- Generate elementary reflector G(i) to annihilate A(i,i+1:n-1)
         i1_=i-1;
         for(i_=1;i_<=n-i;i_++)
            t[i_]=a[i][i_+i1_];
         CReflections::GenerateReflection(t,n-i,ltau);
         taup[i]=ltau;
         i1_=1-i;
         for(i_=i;i_<n;i_++)
            a[i].Set(i_,t[i_+i1_]);
         t[1]=1;
         //--- Apply G(i) to A(i+1:m-1,i:n-1) from the right
         CReflections::ApplyReflectionFromTheRight(a,ltau,t,i+1,m-1,i,n-1,work);
         //--- check
         if(i<m-1)
           {
            //--- Generate elementary reflector H(i) to annihilate
            //--- A(i+2:m-1,i)
            i1_=i;
            for(i_=1;i_<m-i;i_++)
               t[i_]=a[i_+i1_][i];
            CReflections::GenerateReflection(t,m-1-i,ltau);
            tauq[i]=ltau;
            i1_=-i;
            for(i_=i+1;i_<m;i_++)
               a[i_].Set(i,t[i_+i1_]);
            t[1]=1;
            //--- Apply H(i) to A(i+1:m-1,i+1:n-1) from the left
            CReflections::ApplyReflectionFromTheLeft(a,ltau,t,i+1,m-1,i+1,n-1,work);
           }
         else
            tauq[i]=0;
        }
     }
  }
//+------------------------------------------------------------------+
//| Unpacking matrix Q which reduces a matrix to bidiagonal form.    |
//| Input parameters:                                                |
//|     QP          -   matrices Q and P in compact form.            |
//|                     Output of ToBidiagonal subroutine.           |
//|     M           -   number of rows in matrix A.                  |
//|     N           -   number of columns in matrix A.               |
//|     TAUQ        -   scalar factors which are used to form Q.     |
//|                     Output of ToBidiagonal subroutine.           |
//|     QColumns    -   required number of columns in matrix Q.      |
//|                     M>=QColumns>=0.                              |
//| Output parameters:                                               |
//|     Q           -   first QColumns columns of matrix Q.          |
//|                     Array[0..M-1, 0..QColumns-1]                 |
//|                     If QColumns=0, the array is not modified.    |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixBDUnpackQ(CMatrixDouble &qp,const int m,const int n,
                                      double &tauq[],const int qcolumns,CMatrixDouble &q)
  {
//--- check
   if(!CAp::Assert(qcolumns<=m,__FUNCTION__+": QColumns>M!"))
      return;
//--- check
   if(!CAp::Assert(qcolumns>=0,__FUNCTION__+": QColumns<0!"))
      return;
//--- check
   if(m==0 || n==0 || qcolumns==0)
      return;
//--- create variables
   int i=0;
   int j=0;
//--- allocation
   q.Resize(m,qcolumns);
//--- identity matrix
   for(i=0;i<m;i++)
     {
      for(j=0;j<qcolumns;j++)
        {
         if(i==j)
            q[i].Set(j,1);
         else
            q[i].Set(j,0);
        }
     }
//--- get result
   RMatrixBDMultiplyByQ(qp,m,n,tauq,q,m,qcolumns,false,false);
  }
//+------------------------------------------------------------------+
//| Multiplication by matrix Q which reduces matrix A to bidiagonal  |
//| form.                                                            |
//| The algorithm allows pre- or post-multiply by Q or Q'.           |
//| Input parameters:                                                |
//|     QP          -   matrices Q and P in compact form.            |
//|                     Output of ToBidiagonal subroutine.           |
//|     M           -   number of rows in matrix A.                  |
//|     N           -   number of columns in matrix A.               |
//|     TAUQ        -   scalar factors which are used to form Q.     |
//|                     Output of ToBidiagonal subroutine.           |
//|     Z           -   multiplied matrix.                           |
//|                     array[0..ZRows-1,0..ZColumns-1]              |
//|     ZRows       -   number of rows in matrix Z. If FromTheRight= |
//|                     =False, ZRows=M, otherwise ZRows can be      |
//|                     arbitrary.                                   |
//|     ZColumns    -   number of columns in matrix Z. If            |
//|                     FromTheRight=True, ZColumns=M, otherwise     |
//|                     ZColumns can be arbitrary.                   |
//|     FromTheRight -  pre- or post-multiply.                       |
//|     DoTranspose -   multiply by Q or Q'.                         |
//| Output parameters:                                               |
//|     Z           -   product of Z and Q.                          |
//|                     Array[0..ZRows-1,0..ZColumns-1]              |
//|                     If ZRows=0 or ZColumns=0, the array is not   |
//|                     modified.                                    |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixBDMultiplyByQ(CMatrixDouble &qp,const int m,const int n,
                                          double &tauq[],CMatrixDouble &z,const int zrows,
                                          const int zcolumns,const bool fromtheright,
                                          const bool dotranspose)
  {
//--- check
   if(m<=0 || n<=0 || zrows<=0 || zcolumns<=0)
      return;
//--- check
   if(!CAp::Assert((fromtheright && zcolumns==m) || (!fromtheright && zrows==m),__FUNCTION__+": incorrect Z size!"))
      return;
//--- create variables
   int mx=0;
   int i_=0;
   int i1_=0;
   int i=0;
   int i1=0;
   int i2=0;
   int istep=0;
//--- create arrays
   double v[];
   double work[];
//--- initialization
   mx=MathMax(m,n);
   mx=MathMax(mx,zrows);
   mx=MathMax(mx,zcolumns);
//--- allocation
   ArrayResizeAL(v,mx+1);
   ArrayResizeAL(work,mx+1);
//--- check
   if(m>=n)
     {
      //--- setup
      if(fromtheright)
        {
         i1=0;
         i2=n-1;
         istep=1;
        }
      else
        {
         i1=n-1;
         i2=0;
         istep=-1;
        }
      //--- check
      if(dotranspose)
        {
         i=i1;
         i1=i2;
         i2=i;
         istep=-istep;
        }
      //--- Process
      i=i1;
      do
        {
         i1_=i-1;
         for(i_=1;i_<=m-i;i_++)
            v[i_]=qp[i_+i1_][i];
         v[1]=1;
         //--- check
         if(fromtheright)
            CReflections::ApplyReflectionFromTheRight(z,tauq[i],v,0,zrows-1,i,m-1,work);
         else
            CReflections::ApplyReflectionFromTheLeft(z,tauq[i],v,i,m-1,0,zcolumns-1,work);
         i=i+istep;
        }
      while(i!=i2+istep);
     }
   else
     {
      //--- setup
      if(fromtheright)
        {
         i1=0;
         i2=m-2;
         istep=1;
        }
      else
        {
         i1=m-2;
         i2=0;
         istep=-1;
        }
      //--- check
      if(dotranspose)
        {
         i=i1;
         i1=i2;
         i2=i;
         istep=-istep;
        }
      //--- Process
      if(m-1>0)
        {
         i=i1;
         do
           {
            i1_=i;
            for(i_=1;i_<=m-i-1;i_++)
               v[i_]=qp[i_+i1_][i];
            v[1]=1;
            //--- check
            if(fromtheright)
               CReflections::ApplyReflectionFromTheRight(z,tauq[i],v,0,zrows-1,i+1,m-1,work);
            else
               CReflections::ApplyReflectionFromTheLeft(z,tauq[i],v,i+1,m-1,0,zcolumns-1,work);
            i=i+istep;
           }
         while(i!=i2+istep);
        }
     }
  }
//+------------------------------------------------------------------+
//| Unpacking matrix P which reduces matrix A to bidiagonal form.    |
//| The subroutine returns transposed matrix P.                      |
//| Input parameters:                                                |
//|     QP      -   matrices Q and P in compact form.                |
//|                 Output of ToBidiagonal subroutine.               |
//|     M       -   number of rows in matrix A.                      |
//|     N       -   number of columns in matrix A.                   |
//|     TAUP    -   scalar factors which are used to form P.         |
//|                 Output of ToBidiagonal subroutine.               |
//|     PTRows  -   required number of rows of matrix P^T.           |
//|                 N >= PTRows >= 0.                                |
//| Output parameters:                                               |
//|     PT      -   first PTRows columns of matrix P^T               |
//|                 Array[0..PTRows-1, 0..N-1]                       |
//|                 If PTRows=0, the array is not modified.          |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixBDUnpackPT(CMatrixDouble &qp,const int m,const int n,
                                       double &taup[],const int ptrows,CMatrixDouble &pt)
  {
//--- check
   if(!CAp::Assert(ptrows<=n,__FUNCTION__+": PTRows>N!"))
      return;
//--- check
   if(!CAp::Assert(ptrows>=0,__FUNCTION__+": PTRows<0!"))
      return;
//--- check
   if(m==0 || n==0 || ptrows==0)
      return;
//--- create variables
   int i=0;
   int j=0;
//--- allocation
   pt.Resize(ptrows,n);
//--- prepare
   for(i=0;i<=ptrows-1;i++)
     {
      for(j=0;j<n;j++)
        {
         if(i==j)
            pt[i].Set(j,1);
         else
            pt[i].Set(j,0);
        }
     }
//--- get result
   RMatrixBDMultiplyByP(qp,m,n,taup,pt,ptrows,n,true,true);
  }
//+------------------------------------------------------------------+
//| Multiplication by matrix P which reduces matrix A to bidiagonal  |
//| form.                                                            |
//| The algorithm allows pre- or post-multiply by P or P'.           |
//| Input parameters:                                                |
//|     QP          -   matrices Q and P in compact form.            |
//|                     Output of RMatrixBD subroutine.              |
//|     M           -   number of rows in matrix A.                  |
//|     N           -   number of columns in matrix A.               |
//|     TAUP        -   scalar factors which are used to form P.     |
//|                     Output of RMatrixBD subroutine.              |
//|     Z           -   multiplied matrix.                           |
//|                     Array whose indexes range within             |
//|                     [0..ZRows-1,0..ZColumns-1].                  |
//|     ZRows       -   number of rows in matrix Z. If               |
//|                     FromTheRight=False, ZRows=N, otherwise ZRows |
//|                     can be arbitrary.                            |
//|     ZColumns    -   number of columns in matrix Z. If            |
//|                     FromTheRight=True, ZColumns=N, otherwise     |
//|                     ZColumns can be arbitrary.                   |
//|     FromTheRight -  pre- or post-multiply.                       |
//|     DoTranspose -   multiply by P or P'.                         |
//| Output parameters:                                               |
//|     Z - product of Z and P.                                      |
//|                 Array whose indexes range within                 |
//|                 [0..ZRows-1,0..ZColumns-1]. If ZRows=0 or        |
//|                 ZColumns=0, the array is not modified.           |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixBDMultiplyByP(CMatrixDouble &qp,const int m,const int n,
                                          double &taup[],CMatrixDouble &z,const int zrows,
                                          const int zcolumns,const bool fromtheright,
                                          const bool dotranspose)
  {
//--- check
   if(m<=0 || n<=0 || zrows<=0 || zcolumns<=0)
      return;
//--- check
   if(!CAp::Assert((fromtheright && zcolumns==n) || (!fromtheright && zrows==n),__FUNCTION__+": incorrect Z size!"))
      return;
//--- create arrays
   double v[];
   double work[];
//--- create variables
   int i=0;
   int mx=0;
   int i1=0;
   int i2=0;
   int istep=0;
   int i_=0;
   int i1_=0;
//--- initialization
   mx=MathMax(m,n);
   mx=MathMax(mx,zrows);
   mx=MathMax(mx,zcolumns);
//--- allocation
   ArrayResizeAL(v,mx+1);
   ArrayResizeAL(work,mx+1);
//--- check
   if(m>=n)
     {
      //--- setup
      if(fromtheright)
        {
         i1=n-2;
         i2=0;
         istep=-1;
        }
      else
        {
         i1=0;
         i2=n-2;
         istep=1;
        }
      //--- check
      if(!dotranspose)
        {
         i=i1;
         i1=i2;
         i2=i;
         istep=-istep;
        }
      //--- Process
      if(n-1>0)
        {
         i=i1;
         do
           {
            i1_=i;
            for(i_=1;i_<n-i;i_++)
               v[i_]=qp[i][i_+i1_];
            v[1]=1;
            //--- check
            if(fromtheright)
               CReflections::ApplyReflectionFromTheRight(z,taup[i],v,0,zrows-1,i+1,n-1,work);
            else
               CReflections::ApplyReflectionFromTheLeft(z,taup[i],v,i+1,n-1,0,zcolumns-1,work);
            i=i+istep;
           }
         while(i!=i2+istep);
        }
     }
   else
     {
      //--- setup
      if(fromtheright)
        {
         i1=m-1;
         i2=0;
         istep=-1;
        }
      else
        {
         i1=0;
         i2=m-1;
         istep=1;
        }
      //--- check
      if(!dotranspose)
        {
         i=i1;
         i1=i2;
         i2=i;
         istep=-istep;
        }
      //--- Process
      i=i1;
      do
        {
         i1_=i-1;
         for(i_=1;i_<=n-i;i_++)
            v[i_]=qp[i][i_+i1_];
         v[1]=1;
         //--- check
         if(fromtheright)
            CReflections::ApplyReflectionFromTheRight(z,taup[i],v,0,zrows-1,i,n-1,work);
         else
            CReflections::ApplyReflectionFromTheLeft(z,taup[i],v,i,n-1,0,zcolumns-1,work);
         i=i+istep;
        }
      while(i!=i2+istep);
     }
  }
//+------------------------------------------------------------------+
//| Unpacking of the main and secondary diagonals of bidiagonal      |
//| decomposition of matrix A.                                       |
//| Input parameters:                                                |
//|     B   -   output of RMatrixBD subroutine.                      |
//|     M   -   number of rows in matrix B.                          |
//|     N   -   number of columns in matrix B.                       |
//| Output parameters:                                               |
//|     IsUpper -   True, if the matrix is upper bidiagonal.         |
//|                 otherwise IsUpper is False.                      |
//|     D       -   the main diagonal.                               |
//|                 Array whose index ranges within [0..Min(M,N)-1]. |
//|     E       -   the secondary diagonal (upper or lower, depending|
//|                 on the value of IsUpper).                        |
//|                 Array index ranges within [0..Min(M,N)-1], the   |
//|                 last element is not used.                        |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixBDUnpackDiagonals(CMatrixDouble &b,const int m,const int n,
                                              bool &isupper,double &d[],double &e[])
  {
//--- check
   if(m<=0 || n<=0)
      return;
//--- create variables
   int i=0;
//--- check
   if(m>=n)
      isupper=true;
   else
      isupper=false;
//--- check
   if(isupper)
     {
      //--- allocation
      ArrayResizeAL(d,n);
      ArrayResizeAL(e,n);
      //--- get result
      for(i=0;i<n-1;i++)
        {
         d[i]=b[i][i];
         e[i]=b[i][i+1];
        }
      d[n-1]=b[n-1][n-1];
     }
   else
     {
      //--- allocation
      ArrayResizeAL(d,m);
      ArrayResizeAL(e,m);
      //--- get result
      for(i=0;i<m-1;i++)
        {
         d[i]=b[i][i];
         e[i]=b[i+1][i];
        }
      d[m-1]=b[m-1][m-1];
     }
  }
//+------------------------------------------------------------------+
//| Reduction of a square matrix to  upper Hessenberg form:          |
//| Q'*A*Q = H, where Q is an orthogonal matrix, H - Hessenberg      |
//| matrix.                                                          |
//| Input parameters:                                                |
//|     A       -   matrix A with elements [0..N-1, 0..N-1]          |
//|     N       -   size of matrix A.                                |
//| Output parameters:                                               |
//|     A       -   matrices Q and P in  compact form (see below).   |
//|     Tau     -   array of scalar factors which are used to form   |
//|                 matrix Q.                                        |
//|                 Array whose index ranges within [0..N-2]         |
//| Matrix H is located on the main diagonal, on the lower secondary |
//| diagonal and above the main diagonal of matrix A. The elements   |
//| which are used to form matrix Q are situated in array Tau and    |
//| below the lower secondary diagonal of matrix A as follows:       |
//| Matrix Q is represented as a product of elementary reflections   |
//| Q = H(0)*H(2)*...*H(n-2),                                        |
//| where each H(i) is given by                                      |
//| H(i) = 1 - tau * v * (v^T)                                       |
//| where tau is a scalar stored in Tau[I]; v - is a real vector,    |
//| so that v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) stored in             |
//| A(i+2:n-1,i).                                                    |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      October 31, 1992                                            |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixHessenberg(CMatrixDouble &a,const int n,double &tau[])
  {
//--- check
   if(n<=1)
      return;
//--- check
   if(!CAp::Assert(n>=0,__FUNCTION__+": incorrect N!"))
      return;
//--- create arrays
   double t[];
   double work[];
//--- create variables
   int    i=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- allocation
   ArrayResizeAL(tau,n-1);
   ArrayResizeAL(t,n+1);
   ArrayResizeAL(work,n);
   for(i=0;i<n-1;i++)
     {
      //--- Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
      i1_=i;
      for(i_=1;i_<n-i;i_++)
         t[i_]=a[i_+i1_][i];
      CReflections::GenerateReflection(t,n-i-1,v);
      i1_=-i;
      for(i_=i+1;i_<n;i_++)
         a[i_].Set(i,t[i_+i1_]);
      tau[i]=v;
      t[1]=1;
      //--- Apply H(i) to A(1:ihi,i+1:ihi) from the right
      CReflections::ApplyReflectionFromTheRight(a,v,t,0,n-1,i+1,n-1,work);
      //--- Apply H(i) to A(i+1:ihi,i+1:n) from the left
      CReflections::ApplyReflectionFromTheLeft(a,v,t,i+1,n-1,i+1,n-1,work);
     }
  }
//+------------------------------------------------------------------+
//| Unpacking matrix Q which reduces matrix A to upper Hessenberg    |
//| form                                                             |
//| Input parameters:                                                |
//|     A   -   output of RMatrixHessenberg subroutine.              |
//|     N   -   size of matrix A.                                    |
//|     Tau -   scalar factors which are used to form Q.             |
//|             Output of RMatrixHessenberg subroutine.              |
//| Output parameters:                                               |
//|     Q   -   matrix Q.                                            |
//|             Array whose indexes range within [0..N-1, 0..N-1].   |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixHessenbergUnpackQ(CMatrixDouble &a,const int n,
                                              double &tau[],CMatrixDouble &q)
  {
//--- check
   if(n==0)
      return;
//--- create arrays
   double v[];
   double work[];
//--- create variables
   int i=0;
   int j=0;
   int i_=0;
   int i1_=0;
//--- allocation
   q.Resize(n,n);
//--- allocation
   ArrayResizeAL(v,n);
   ArrayResizeAL(work,n);
//--- identity matrix
   for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
        {
         if(i==j)
            q[i].Set(j,1);
         else
            q[i].Set(j,0);
        }
     }
//--- unpack Q
   for(i=0;i<n-1;i++)
     {
      //--- Apply H(i)
      i1_=i;
      for(i_=1;i_<n-i;i_++)
         v[i_]=a[i_+i1_][i];
      v[1]=1;
      CReflections::ApplyReflectionFromTheRight(q,tau[i],v,0,n-1,i+1,n-1,work);
     }
  }
//+------------------------------------------------------------------+
//| Unpacking matrix H (the result of matrix A reduction to upper    |
//| Hessenberg form)                                                 |
//| Input parameters:                                                |
//|     A   -   output of RMatrixHessenberg subroutine.              |
//|     N   -   size of matrix A.                                    |
//| Output parameters:                                               |
//|     H   -   matrix H. Array whose indexes range within           |
//|     [0..N-1, 0..N-1].                                            |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixHessenbergUnpackH(CMatrixDouble &a,const int n,CMatrixDouble &h)
  {
//--- check
   if(n==0)
      return;
//--- create arrays
   double v[];
   double work[];
//--- create variables
   int i=0;
   int j=0;
   int i_=0;
//--- allocation
   h.Resize(n,n);
//--- get result
   for(i=0;i<n;i++)
     {
      for(j=0;j<=i-2;j++)
         h[i].Set(j,0);
      j=(int)MathMax(0,i-1);
      for(i_=j;i_<n;i_++)
         h[i].Set(i_,a[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| Reduction of a symmetric matrix which is given by its higher or  |
//| lower triangular part to a tridiagonal matrix using orthogonal   |
//| similarity transformation: Q'*A*Q=T.                             |
//| Input parameters:                                                |
//|     A       -   matrix to be transformed                         |
//|                 array with elements [0..N-1, 0..N-1].            |
//|     N       -   size of matrix A.                                |
//|     IsUpper -   storage format. If IsUpper = True, then matrix A |
//|                 is given by its upper triangle, and the lower    |
//|                 triangle is not used and not modified by the     |
//|                 algorithm, and vice versa if IsUpper = False.    |
//| Output parameters:                                               |
//|     A       -   matrices T and Q in  compact form (see lower)    |
//|     Tau     -   array of factors which are forming matrices H(i) |
//|                 array with elements [0..N-2].                    |
//|     D       -   main diagonal of symmetric matrix T.             |
//|                 array with elements [0..N-1].                    |
//|     E       -   secondary diagonal of symmetric matrix T.        |
//|                 array with elements [0..N-2].                    |
//|   If IsUpper=True, the matrix Q is represented as a product of   |
//|   elementary reflectors                                          |
//|      Q = H(n-2) . . . H(2) H(0).                                 |
//|   Each H(i) has the form                                         |
//|      H(i) = I - tau * v * v'                                     |
//|   where tau is a real scalar, and v is a real vector with        |
//|   v(i+1:n-1) = 0, v(i) = 1, v(0:i-1) is stored on exit in        |
//|   A(0:i-1,i+1), and tau in TAU(i).                               |
//|   If IsUpper=False, the matrix Q is represented as a product of  |
//|   elementary reflectors                                          |
//|      Q = H(0) H(2) . . . H(n-2).                                 |
//|   Each H(i) has the form                                         |
//|      H(i) = I - tau * v * v'                                     |
//|   where tau is a real scalar, and v is a real vector with        |
//|   v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) is stored on exit in        |
//|   A(i+2:n-1,i), and tau in TAU(i).                               |
//|   The contents of A on exit are illustrated by the following     |
//|   examples with n = 5:                                           |
//|   if UPLO = 'U':                       if UPLO = 'L':            |
//|     (  d   e   v1  v2  v3 )              (  d                  ) |
//|     (      d   e   v2  v3 )              (  e   d              ) |
//|     (          d   e   v3 )              (  v0  e   d          ) |
//|     (              d   e  )              (  v0  v1  e   d      ) |
//|     (                  d  )              (  v0  v1  v2  e   d  ) |
//|   where d and e denote diagonal and off-diagonal elements of T,  |
//|   and vi denotes an element of the vector defining H(i).         |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      October 31, 1992                                            |
//+------------------------------------------------------------------+
static void COrtFac::SMatrixTD(CMatrixDouble &a,const int n,const bool isupper,
                               double &tau[],double &d[],double &e[])
  {
//--- check
   if(n<=0)
      return;
//--- create arrays
   double t[];
   double t2[];
   double t3[];
//--- create variables
   int    i=0;
   double alpha=0;
   double taui=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- allocation
   ArrayResizeAL(t,n+1);
   ArrayResizeAL(t2,n+1);
   ArrayResizeAL(t3,n+1);
//--- check
   if(n>1)
      ArrayResizeAL(tau,n-1);
   ArrayResizeAL(d,n);
//--- check
   if(n>1)
      ArrayResizeAL(e,n-1);
//--- check
   if(isupper)
     {
      //--- Reduce the upper triangle of A
      for(i=n-2;i>=0;i--)
        {
         //--- Generate elementary reflector H() = E - tau * v * v'
         if(i>=1)
           {
            i1_=-2;
            for(i_=2;i_<=i+1;i_++)
               t[i_]=a[i_+i1_][i+1];
           }
         t[1]=a[i][i+1];
         CReflections::GenerateReflection(t,i+1,taui);
         //--- check
         if(i>=1)
           {
            i1_=2;
            for(i_=0;i_<i;i_++)
               a[i_].Set(i+1,t[i_+i1_]);
           }
         a[i].Set(i+1,t[1]);
         e[i]=a[i][i+1];
         //--- check
         if(taui!=0)
           {
            //--- Apply H from both sides to A
            a[i].Set(i+1,1);
            //--- Compute  x := tau * A * v  storing x in TAU
            i1_=-1;
            for(i_=1;i_<=i+1;i_++)
               t[i_]=a[i_+i1_][i+1];
            CSblas::SymmetricMatrixVectorMultiply(a,isupper,0,i,t,taui,t3);
            i1_=1;
            for(i_=0;i_<=i;i_++)
               tau[i_]=t3[i_+i1_];
            //--- Compute  w := x - 1/2 * tau * (x'*v) * v
            v=0.0;
            for(i_=0;i_<=i;i_++)
               v+=tau[i_]*a[i_][i+1];
            alpha=-(0.5*taui*v);
            for(i_=0;i_<=i;i_++)
               tau[i_]=tau[i_]+alpha*a[i_][i+1];
            //--- Apply the transformation as a rank-2 update:
            //---    A := A - v * w' - w * v'
            i1_=-1;
            for(i_=1;i_<=i+1;i_++)
               t[i_]=a[i_+i1_][i+1];
            i1_=-1;
            for(i_=1;i_<=i+1;i_++)
               t3[i_]=tau[i_+i1_];
            CSblas::SymmetricRank2Update(a,isupper,0,i,t,t3,t2,-1);
            a[i].Set(i+1,e[i]);
           }
         d[i+1]=a[i+1][i+1];
         tau[i]=taui;
        }
      d[0]=a[0][0];
     }
   else
     {
      //--- Reduce the lower triangle of A
      for(i=0;i<n-1;i++)
        {
         //--- Generate elementary reflector H = E - tau * v * v'
         i1_=i;
         for(i_=1;i_<n-i;i_++)
            t[i_]=a[i_+i1_][i];
         CReflections::GenerateReflection(t,n-i-1,taui);
         i1_=-i;
         for(i_=i+1;i_<n;i_++)
            a[i_].Set(i,t[i_+i1_]);
         e[i]=a[i+1][i];
         if(taui!=0)
           {
            //--- Apply H from both sides to A
            a[i+1].Set(i,1);
            //--- Compute  x := tau * A * v  storing y in TAU
            i1_=i;
            for(i_=1;i_<n-i;i_++)
               t[i_]=a[i_+i1_][i];
            CSblas::SymmetricMatrixVectorMultiply(a,isupper,i+1,n-1,t,taui,t2);
            i1_=1-i;
            for(i_=i;i_<n-1;i_++)
               tau[i_]=t2[i_+i1_];
            //--- Compute  w := x - 1/2 * tau * (x'*v) * v
            i1_=1;
            v=0.0;
            for(i_=i;i_<=n-2;i_++)
               v+=tau[i_]*a[i_+i1_][i];
            alpha=-(0.5*taui*v);
            i1_=1;
            for(i_=i;i_<n-1;i_++)
               tau[i_]=tau[i_]+alpha*a[i_+i1_][i];
            //--- Apply the transformation as a rank-2 update:
            //---     A := A - v * w' - w * v'
            i1_=i;
            for(i_=1;i_<n-i;i_++)
               t[i_]=a[i_+i1_][i];
            i1_=i-1;
            for(i_=1;i_<n-i;i_++)
               t2[i_]=tau[i_+i1_];
            CSblas::SymmetricRank2Update(a,isupper,i+1,n-1,t,t2,t3,-1);
            a[i+1].Set(i,e[i]);
           }
         d[i]=a[i][i];
         tau[i]=taui;
        }
      d[n-1]=a[n-1][n-1];
     }
  }
//+------------------------------------------------------------------+
//| Unpacking matrix Q which reduces symmetric matrix to a           |
//| tridiagonal form.                                                |
//| Input parameters:                                                |
//|     A       -   the result of a SMatrixTD subroutine             |
//|     N       -   size of matrix A.                                |
//|     IsUpper -   storage format (a parameter of SMatrixTD         |
//|                 subroutine)                                      |
//|     Tau     -   the result of a SMatrixTD subroutine             |
//| Output parameters:                                               |
//|     Q       -   transformation matrix.                           |
//|                 array with elements [0..N-1, 0..N-1].            |
//+------------------------------------------------------------------+
static void COrtFac::SMatrixTDUnpackQ(CMatrixDouble &a,const int n,const bool isupper,
                                      double &tau[],CMatrixDouble &q)
  {
//--- check
   if(n==0)
      return;
//--- create arrays
   double v[];
   double work[];
//--- create variables
   int i=0;
   int j=0;
   int i_=0;
   int i1_=0;
//--- allocation
   q.Resize(n,n);
//--- allocation
   ArrayResizeAL(v,n+1);
   ArrayResizeAL(work,n);
//--- identity matrix
   for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
        {
         if(i==j)
            q[i].Set(j,1);
         else
            q[i].Set(j,0);
        }
     }
//--- unpack Q
   if(isupper)
     {
      for(i=0;i<n-1;i++)
        {
         //--- Apply H(i)
         i1_=-1;
         for(i_=1;i_<=i+1;i_++)
            v[i_]=a[i_+i1_][i+1];
         v[i+1]=1;
         //--- function call
         CReflections::ApplyReflectionFromTheLeft(q,tau[i],v,0,i,0,n-1,work);
        }
     }
   else
     {
      for(i=n-2;i>=0;i--)
        {
         //--- Apply H(i)
         i1_=i;
         for(i_=1;i_<n-i;i_++)
            v[i_]=a[i_+i1_][i];
         v[1]=1;
         //--- function call
         CReflections::ApplyReflectionFromTheLeft(q,tau[i],v,i+1,n-1,0,n-1,work);
        }
     }
  }
//+------------------------------------------------------------------+
//| Reduction of a Hermitian matrix which is given by its higher or  |
//| lower triangular part to a real tridiagonal matrix using unitary |
//| similarity transformation: Q'*A*Q = T.                           |
//| Input parameters:                                                |
//|     A       -   matrix to be transformed                         |
//|                 array with elements [0..N-1, 0..N-1].            |
//|     N       -   size of matrix A.                                |
//|     IsUpper -   storage format. If IsUpper = True, then matrix A |
//|                 is given by its upper triangle, and the lower    |
//|                 triangle is not used and not modified by the     |
//|                 algorithm, and vice versa if IsUpper = False.    |
//| Output parameters:                                               |
//|     A       -   matrices T and Q in  compact form (see lower)    |
//|     Tau     -   array of factors which are forming matrices H(i) |
//|                 array with elements [0..N-2].                    |
//|     D       -   main diagonal of real symmetric matrix T.        |
//|                 array with elements [0..N-1].                    |
//|     E       -   secondary diagonal of real symmetric matrix T.   |
//|                 array with elements [0..N-2].                    |
//|   If IsUpper=True, the matrix Q is represented as a product of   |
//|   elementary reflectors                                          |
//|      Q = H(n-2) . . . H(2) H(0).                                 |
//|   Each H(i) has the form                                         |
//|      H(i) = I - tau * v * v'                                     |
//|   where tau is a complex scalar, and v is a complex vector with  |
//|   v(i+1:n-1) = 0, v(i) = 1, v(0:i-1) is stored on exit in        |
//|   A(0:i-1,i+1), and tau in TAU(i).                               |
//|   If IsUpper=False, the matrix Q is represented as a product of  |
//|   elementary reflectors                                          |
//|      Q = H(0) H(2) . . . H(n-2).                                 |
//|   Each H(i) has the form                                         |
//|      H(i) = I - tau * v * v'                                     |
//|   where tau is a complex scalar, and v is a complex vector with  |
//|   v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) is stored on exit in        |
//|   A(i+2:n-1,i), and tau in TAU(i).                               |
//|   The contents of A on exit are illustrated by the following     |
//|   examples with n = 5:                                           |
//|   if UPLO = 'U':                       if UPLO = 'L':            |
//|     (  d   e   v1  v2  v3 )              (  d                  ) |
//|     (      d   e   v2  v3 )              (  e   d              ) |
//|     (          d   e   v3 )              (  v0  e   d          ) |
//|     (              d   e  )              (  v0  v1  e   d      ) |
//|     (                  d  )              (  v0  v1  v2  e   d  ) |
//| where d and e denote diagonal and off-diagonal elements of T, and|
//| vi denotes an element of the vector defining H(i).               |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      October 31, 1992                                            |
//+------------------------------------------------------------------+
static void COrtFac::HMatrixTD(CMatrixComplex &a,const int n,const bool isupper,
                               al_complex &tau[],double &d[],double &e[])
  {
//--- check
   if(n<=0)
      return;
//--- create arrays
   al_complex t[];
   al_complex t2[];
   al_complex t3[];
//--- create variables
   al_complex Half(0.5,0);
   al_complex Zero(0,0);
   al_complex _One(-1,0);
   al_complex alpha=0;
   al_complex taui=0;
   al_complex v=0;
   int        i=0;
   int        i_=0;
   int        i1_=0;
   for(i=0;i<n;i++)
     {
      //--- check
      if(!CAp::Assert(a[i][i].im==0))
         return;
     }
//--- allocation
   if(n>1)
     {
      ArrayResizeAL(tau,n-1);
      ArrayResizeAL(e,n-1);
     }
   ArrayResizeAL(d,n);
   ArrayResizeAL(t,n);
   ArrayResizeAL(t2,n);
   ArrayResizeAL(t3,n);
//--- check
   if(isupper)
     {
      //--- Reduce the upper triangle of A
      a[n-1].Set(n-1,a[n-1][n-1].re);
      for(i=n-2;i>=0;i--)
        {
         //--- Generate elementary reflector H = I+1 - tau * v * v'
         alpha=a[i][i+1];
         t[1]=alpha;
         //--- check
         if(i>=1)
           {
            i1_=-2;
            for(i_=2;i_<=i+1;i_++)
               t[i_]=a[i_+i1_][i+1];
           }
         //--- function call
         CComplexReflections::ComplexGenerateReflection(t,i+1,taui);
         //--- check
         if(i>=1)
           {
            i1_=2;
            for(i_=0;i_<i;i_++)
               a[i_].Set(i+1,t[i_+i1_]);
           }
         //--- change values
         alpha=t[1];
         e[i]=alpha.re;
         //--- check
         if(taui!=Zero)
           {
            //--- Apply H(I+1) from both sides to A
            a[i].Set(i+1,1);
            //--- Compute  x := tau * A * v  storing x in TAU
            i1_=-1;
            for(i_=1;i_<=i+1;i_++)
               t[i_]=a[i_+i1_][i+1];
            CHblas::HermitianMatrixVectorMultiply(a,isupper,0,i,t,taui,t2);
            i1_=1;
            for(i_=0;i_<=i;i_++)
               tau[i_]=t2[i_+i1_];
            //--- Compute  w := x - 1/2 * tau * (x'*v) * v
            v=0.0;
            for(i_=0;i_<=i;i_++)
               v+=CMath::Conj(tau[i_])*a[i_][i+1];
            //--- calculation
            alpha=Half*taui*v;
            alpha.re=-alpha.re;
            alpha.im=-alpha.im;
            for(i_=0;i_<=i;i_++)
               tau[i_]=tau[i_]+alpha*a[i_][i+1];
            //--- Apply the transformation as a rank-2 update:
            //---    A := A - v * w' - w * v'
            i1_=-1;
            for(i_=1;i_<=i+1;i_++)
               t[i_]=a[i_+i1_][i+1];
            i1_=-1;
            for(i_=1;i_<=i+1;i_++)
               t3[i_]=tau[i_+i1_];
            CHblas::HermitianRank2Update(a,isupper,0,i,t,t3,t2,_One);
           }
         else
            a[i].Set(i,a[i][i].re);
         //--- change values
         a[i].Set(i+1,e[i]);
         d[i+1]=a[i+1][i+1].re;
         tau[i]=taui;
        }
      d[0]=a[0][0].re;
     }
   else
     {
      //--- Reduce the lower triangle of A
      a[0].Set(0,a[0][0].re);
      for(i=0;i<n-1;i++)
        {
         //--- Generate elementary reflector H = I - tau * v * v'
         i1_=i;
         for(i_=1;i_<n-i;i_++)
            t[i_]=a[i_+i1_][i];
         //--- function call
         CComplexReflections::ComplexGenerateReflection(t,n-i-1,taui);
         i1_=-i;
         for(i_=i+1;i_<n;i_++)
            a[i_].Set(i,t[i_+i1_]);
         e[i]=a[i+1][i].re;
         //--- check
         if(taui!=Zero)
           {
            //--- Apply H(i) from both sides to A(i+1:n,i+1:n)
            a[i+1].Set(i,1);
            //--- Compute  x := tau * A * v  storing y in TAU
            i1_=i;
            for(i_=1;i_<n-i;i_++)
               t[i_]=a[i_+i1_][i];
            CHblas::HermitianMatrixVectorMultiply(a,isupper,i+1,n-1,t,taui,t2);
            i1_=1-i;
            for(i_=i;i_<n-1;i_++)
               tau[i_]=t2[i_+i1_];
            //--- Compute  w := x - 1/2 * tau * (x'*v) * v
            i1_=1;
            v=0.0;
            for(i_=i;i_<n-1;i_++)
               v+=CMath::Conj(tau[i_])*a[i_+i1_][i];
            //--- calculation
            alpha=Half*taui*v;
            alpha.re=-alpha.re;
            alpha.im=-alpha.im;
            i1_=1;
            for(i_=i;i_<n-1;i_++)
               tau[i_]=tau[i_]+alpha*a[i_+i1_][i];
            //--- Apply the transformation as a rank-2 update:
            //--- A := A - v * w' - w * v'
            i1_=i;
            for(i_=1;i_<n-i;i_++)
               t[i_]=a[i_+i1_][i];
            i1_=i-1;
            for(i_=1;i_<n-i;i_++)
               t2[i_]=tau[i_+i1_];
            CHblas::HermitianRank2Update(a,isupper,i+1,n-1,t,t2,t3,_One);
           }
         else
            a[i+1].Set(i+1,a[i+1][i+1].re);
         //--- change values
         a[i+1].Set(i,e[i]);
         d[i]=a[i][i].re;
         tau[i]=taui;
        }
      d[n-1]=a[n-1][n-1].re;
     }
  }
//+------------------------------------------------------------------+
//| Unpacking matrix Q which reduces a Hermitian matrix to a real    |
//| tridiagonal form.                                                |
//| Input parameters:                                                |
//|     A       -   the result of a HMatrixTD subroutine             |
//|     N       -   size of matrix A.                                |
//|     IsUpper -   storage format (a parameter of HMatrixTD         |
//|                 subroutine)                                      |
//|     Tau     -   the result of a HMatrixTD subroutine             |
//| Output parameters:                                               |
//|     Q       -   transformation matrix.                           |
//|                 array with elements [0..N-1, 0..N-1].            |
//+------------------------------------------------------------------+
static void COrtFac::HMatrixTDUnpackQ(CMatrixComplex &a,const int n,const bool isupper,
                                      al_complex &tau[],CMatrixComplex &q)
  {
//--- check
   if(n==0)
      return;
//--- create arrays
   al_complex v[];
   al_complex work[];
//--- create variables
   int i=0;
   int j=0;
   int i_=0;
   int i1_=0;
//--- allocation
   q.Resize(n,n);
//--- allocation
   ArrayResizeAL(v,n+1);
   ArrayResizeAL(work,n+1);
//--- identity matrix
   for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
        {
         if(i==j)
            q[i].Set(j,1);
         else
            q[i].Set(j,0);
        }
     }
//--- unpack Q
   if(isupper)
     {
      for(i=0;i<n-1;i++)
        {
         //--- Apply H(i)
         i1_=-1;
         for(i_=1;i_<=i+1;i_++)
            v[i_]=a[i_+i1_][i+1];
         v[i+1]=1;
         //--- function call
         CComplexReflections::ComplexApplyReflectionFromTheLeft(q,tau[i],v,0,i,0,n-1,work);
        }
     }
   else
     {
      for(i=n-2;i>=0;i--)
        {
         //--- Apply H(i)
         i1_=i;
         for(i_=1;i_<n-i;i_++)
            v[i_]=a[i_+i1_][i];
         v[1]=1;
         //--- function call
         CComplexReflections::ComplexApplyReflectionFromTheLeft(q,tau[i],v,i+1,n-1,0,n-1,work);
        }
     }
  }
//+------------------------------------------------------------------+
//| Base case for real QR                                            |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixQRBaseCase(CMatrixDouble &a,const int m,const int n,
                                       double &work[],double &t[],double &tau[])
  {
//--- create variables
   int    i=0;
   int    k=MathMin(m,n);
   int    minmn=MathMin(m,n);
   double tmp=0;
   int    i_=0;
   int    i1_=0;
//--- calculations
   for(i=0;i<k;i++)
     {
      //--- Generate elementary reflector H(i) to annihilate A(i+1:m,i)
      i1_=i-1;
      for(i_=1;i_<=m-i;i_++)
         t[i_]=a[i_+i1_][i];
      CReflections::GenerateReflection(t,m-i,tmp);
      tau[i]=tmp;
      i1_=1-i;
      for(i_=i;i_<m;i_++)
         a[i_].Set(i,t[i_+i1_]);
      t[1]=1;
      //--- check
      if(i<n)
        {
         //--- Apply H(i) to A(i:m-1,i+1:n-1) from the left
         CReflections::ApplyReflectionFromTheLeft(a,tau[i],t,i,m-1,i+1,n-1,work);
        }
     }
  }
//+------------------------------------------------------------------+
//| Base case for real LQ                                            |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixLQBaseCase(CMatrixDouble &a,const int m,const int n,
                                       double &work[],double &t[],double &tau[])
  {
//--- create variables
   int    i=0;
   int    k=MathMin(m,n);
   int    minmn=MathMin(m,n);
   double tmp=0;
   int    i_=0;
   int    i1_=0;
//--- calculation
   for(i=0;i<k;i++)
     {
      //--- Generate elementary reflector H(i) to annihilate A(i,i+1:n-1)
      i1_=i-1;
      for(i_=1;i_<=n-i;i_++)
         t[i_]=a[i][i_+i1_];
      CReflections::GenerateReflection(t,n-i,tmp);
      tau[i]=tmp;
      i1_=1-i;
      for(i_=i;i_<n;i_++)
         a[i].Set(i_,t[i_+i1_]);
      t[1]=1;
      //--- check
      if(i<n)
        {
         //--- Apply H(i) to A(i+1:m,i:n) from the right
         CReflections::ApplyReflectionFromTheRight(a,tau[i],t,i+1,m-1,i,n-1,work);
        }
     }
  }
//+------------------------------------------------------------------+
//| Base case for complex QR                                         |
//+------------------------------------------------------------------+
static void COrtFac::CMatrixQRBaseCase(CMatrixComplex &a,const int m,const int n,
                                       al_complex &work[],al_complex &t[],al_complex &tau[])
  {
//--- create variables
   int        i=0;
   int        k=MathMin(m,n);
   int        mmi=0;
   int        minmn=MathMin(m,n);
   al_complex tmp=0;
   int        i_=0;
   int        i1_=0;
//--- check
   if(minmn<=0)
      return;
//--- calculation
   for(i=0;i<k;i++)
     {
      //--- Generate elementary reflector H(i) to annihilate A(i+1:m,i)
      mmi=m-i;
      i1_=i-1;
      for(i_=1;i_<=mmi;i_++)
         t[i_]=a[i_+i1_][i];
      //--- function call
      CComplexReflections::ComplexGenerateReflection(t,mmi,tmp);
      tau[i]=tmp;
      i1_=1-i;
      for(i_=i;i_<m;i_++)
         a[i_].Set(i,t[i_+i1_]);
      t[1]=1;
      //--- check
      if(i<n-1)
        {
         //--- Apply H'(i) to A(i:m,i+1:n) from the left
         CComplexReflections::ComplexApplyReflectionFromTheLeft(a,CMath::Conj(tau[i]),t,i,m-1,i+1,n-1,work);
        }
     }
  }
//+------------------------------------------------------------------+
//| Base case for complex LQ                                         |
//+------------------------------------------------------------------+
static void COrtFac::CMatrixLQBaseCase(CMatrixComplex &a,const int m,const int n,
                                       al_complex &work[],al_complex &t[],al_complex &tau[])
  {
//--- create variables
   int        i=0;
   int        minmn=MathMin(m,n);
   al_complex tmp=0;
   int        i_=0;
   int        i1_=0;
//--- check
   if(minmn<=0)
      return;
//--- calculation
   for(i=0;i<=minmn-1;i++)
     {
      //--- Generate elementary reflector H(i)
      //--- NOTE: ComplexGenerateReflection() generates left reflector,
      //--- i.e. H which reduces x by applyiong from the left, but we
      //--- need RIGHT reflector. So we replace H=E-tau*v*v' by H^H,
      //--- which changes v to conj(v).
      i1_=i-1;
      for(i_=1;i_<=n-i;i_++)
         t[i_]=CMath::Conj(a[i][i_+i1_]);
      CComplexReflections::ComplexGenerateReflection(t,n-i,tmp);
      tau[i]=tmp;
      i1_=1-i;
      for(i_=i;i_<n;i_++)
         a[i].Set(i_,CMath::Conj(t[i_+i1_]));
      t[1]=1;
      //--- check
      if(i<m-1)
        {
         //--- Apply H'(i)
         CComplexReflections::ComplexApplyReflectionFromTheRight(a,tau[i],t,i+1,m-1,i,n-1,work);
        }
     }
  }
//+------------------------------------------------------------------+
//| Generate block reflector:                                        |
//| * fill unused parts of reflectors matrix by zeros                |
//| * fill diagonal of reflectors matrix by ones                     |
//| * generate triangular factor T                                   |
//| PARAMETERS:                                                      |
//|     A           -   either LengthA*BlockSize (if ColumnwiseA) or |
//|                     BlockSize*LengthA (if not ColumnwiseA) matrix|
//|                     of elementary reflectors.                    |
//|                     Modified on exit.                            |
//|     Tau         -   scalar factors                               |
//|     ColumnwiseA -   reflectors are stored in rows or in columns  |
//|     LengthA     -   length of largest reflector                  |
//|     BlockSize   -   number of reflectors                         |
//|     T           -   array[BlockSize,2*BlockSize]. Left           |
//|                     BlockSize*BlockSize submatrix stores         |
//|                     triangular factor on exit.                   |
//|     WORK        -   array[BlockSize]                             |
//+------------------------------------------------------------------+
static void COrtFac::RMatrixBlockReflector(CMatrixDouble &a,double &tau[],
                                           const bool columnwisea,const int lengtha,
                                           const int blocksize,CMatrixDouble &t,double &work[])
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- fill beginning of new column with zeros,
//--- load 1.0 in the first non-zero element
   for(k=0;k<blocksize;k++)
     {
      //--- check
      if(columnwisea)
        {
         for(i=0;i<k;i++)
            a[i].Set(k,0);
        }
      else
        {
         for(i=0;i<k;i++)
            a[k].Set(i,0);
        }
      a[k].Set(k,1);
     }
//--- Calculate Gram matrix of A
   for(i=0;i<blocksize;i++)
     {
      for(j=0;j<blocksize;j++)
         t[i].Set(blocksize+j,0);
     }
   for(k=0;k<=lengtha-1;k++)
     {
      for(j=1;j<blocksize;j++)
        {
         //--- check
         if(columnwisea)
           {
            v=a[k][j];
            //--- check
            if(v!=0)
              {
               i1_=-blocksize;
               for(i_=blocksize;i_<=blocksize+j-1;i_++)
                  t[j].Set(i_,t[j][i_]+v*a[k][i_+i1_]);
              }
           }
         else
           {
            v=a[j][k];
            //--- check
            if(v!=0)
              {
               i1_=-blocksize;
               for(i_=blocksize;i_<=blocksize+j-1;i_++)
                  t[j].Set(i_,t[j][i_]+v*a[i_+i1_][k]);
              }
           }
        }
     }
//--- Prepare Y (stored in TmpA) and T (stored in TmpT)
   for(k=0;k<blocksize;k++)
     {
      //--- fill non-zero part of T, use pre-calculated Gram matrix
      i1_=blocksize;
      for(i_=0;i_<k;i_++)
         work[i_]=t[k][i_+i1_];
      for(i=0;i<k;i++)
        {
         v=0.0;
         for(i_=i;i_<k;i_++)
            v+=t[i][i_]*work[i_];
         t[i].Set(k,-(tau[k]*v));
        }
      t[k].Set(k,-tau[k]);
      //--- Rest of T is filled by zeros
      for(i=k+1;i<blocksize;i++)
         t[i].Set(k,0);
     }
  }
//+------------------------------------------------------------------+
//| Generate block reflector (complex):                              |
//| * fill unused parts of reflectors matrix by zeros                |
//| * fill diagonal of reflectors matrix by ones                     |
//| * generate triangular factor T                                   |
//+------------------------------------------------------------------+
static void COrtFac::CMatrixBlockReflector(CMatrixComplex &a,al_complex &tau[],
                                           const bool columnwisea,const int lengtha,
                                           const int blocksize,CMatrixComplex &t,al_complex &work[])
  {
//--- create variables
   int        i=0;
   int        k=0;
   al_complex v=0;
   al_complex tauv=0;
   int        i_=0;
//--- Prepare Y (stored in TmpA) and T (stored in TmpT)
   for(k=0;k<blocksize;k++)
     {
      //--- fill beginning of new column with zeros,
      //--- load 1.0 in the first non-zero element
      if(columnwisea)
        {
         for(i=0;i<k;i++)
            a[i].Set(k,0);
        }
      else
        {
         for(i=0;i<k;i++)
            a[k].Set(i,0);
        }
      a[k].Set(k,1);
      //--- fill non-zero part of T
      for(i=0;i<k;i++)
        {
         //--- check
         if(columnwisea)
           {
            v=0.0;
            for(i_=k;i_<=lengtha-1;i_++)
               v+=CMath::Conj(a[i_][i])*a[i_][k];
           }
         else
           {
            v=0.0;
            for(i_=k;i_<=lengtha-1;i_++)
               v+=a[i][i_]*CMath::Conj(a[k][i_]);
           }
         work[i]=v;
        }
      for(i=0;i<k;i++)
        {
         v=0.0;
         for(i_=i;i_<k;i_++)
            v+=t[i][i_]*work[i_];
         //--- change
         tauv=tau[k]*v;
         tauv.re=-tauv.re;
         tauv.im=-tauv.im;
         t[i].Set(k,tauv);
        }
      //--- change
      tauv=tau[k];
      tauv.re=-tauv.re;
      tauv.im=-tauv.im;
      t[k].Set(k,tauv);
      //--- Rest of T is filled by zeros
      for(i=k+1;i<blocksize;i++)
         t[i].Set(k,0);
     }
  }
//+------------------------------------------------------------------+
//| Eigenvalues and eigenvectors                                     |
//+------------------------------------------------------------------+
class CEigenVDetect
  {
private:
   //--- private methods
   static bool       TtidiagonalEVD(double &d[],double &ce[],const int n,const int zneeded,CMatrixDouble &z);
   static void       TdEVDE2(const double a,const double b,const double c,double &rt1,double &rt2);
   static void       TdEVDEv2(const double a,const double b,const double c,double &rt1,double &rt2,double &cs1,double &sn1);
   static double     TdEVDPythag(const double a,const double b);
   static double     TdEVDExtSign(const double a,const double b);
   static bool       InternalBisectionEigenValues(double &cd[],double &ce[],const int n,int irange,const int iorder,const double vl,const double vu,const int il,const int iu,const double abstol,double &w[],int &m,int &nsplit,int &iblock[],int &isplit[],int &errorcode);
   static void       InternalDStein(const int n,double &d[],double &ce[],const int m,double &cw[],int &iblock[],int &isplit[],CMatrixDouble &z,int &ifail[],int &info);
   static void       TdIninternalDLAGTF(const int n,double &a[],const double lambdav,double &b[],double &c[],double tol,double &d[],int &iin[],int &info);
   static void       TdIninternalDLAGTS(const int n,double &a[],double &b[],double &c[],double &d[],int &iin[],double &y[],double &tol,int &info);
   static void       InternalDLAEBZ(const int ijob,const int nitmax,const int n,const int mmax,const int minp,const double abstol,const double reltol,const double pivmin,double &d[],double &e[],double &e2[],int &nval[],CMatrixDouble &ab,double &c[],int &mout,CMatrixInt &nab,double &work[],int &iwork[],int &info);
   static void       InternalTREVC(CMatrixDouble &t,const int n,const int side,const int howmny,bool &cvselect[],CMatrixDouble &vl,CMatrixDouble &vr,int &m,int &info);
   static void       InternalHsEVDLALN2(const bool ltrans,const int na,const int nw,const double smin,const double ca,CMatrixDouble &a,const double d1,const double d2,CMatrixDouble &b,const double wr,const double wi,bool &rswap4[],bool &zswap4[],CMatrixInt &ipivot44,double &civ4[],double &crv4[],CMatrixDouble &x,double &scl,double &xnorm,int &info);
   static void       InternalHsEVDLADIV(const double a,const double b,const double c,const double d,double &p,double &q);
   static bool       NonSymmetricEVD(CMatrixDouble &ca,const int n,const int vneeded,double &wr[],double &wi[],CMatrixDouble &vl,CMatrixDouble &vr);
   static void       ToUpperHessenberg(CMatrixDouble &a,const int n,double &tau[]);
   static void       UnpackQFromUpperHessenberg(CMatrixDouble &a,const int n,double &tau[],CMatrixDouble &q);
public:
                     CEigenVDetect(void);
                    ~CEigenVDetect(void);
   //--- public methods
   static bool       SMatrixEVD(CMatrixDouble &ca,const int n,const int zneeded,const bool isupper,double &d[],CMatrixDouble &z);
   static bool       SMatrixEVDR(CMatrixDouble &ca,const int n,const int zneeded,const bool isupper,const double b1,const double b2,int &m,double &w[],CMatrixDouble &z);
   static bool       SMatrixEVDI(CMatrixDouble &ca,const int n,const int zneeded,const bool isupper,const int i1,const int i2,double &w[],CMatrixDouble &z);
   static bool       HMatrixEVD(CMatrixComplex &ca,const int n,int zneeded,const bool isupper,double &d[],CMatrixComplex &z);
   static bool       HMatrixEVDR(CMatrixComplex &ca,const int n,int zneeded,bool isupper,const double b1,const double b2,int &m,double &w[],CMatrixComplex &z);
   static bool       HMatrixEVDI(CMatrixComplex &ca,const int n,int zneeded,const bool isupper,const int i1,const int i2,double &w[],CMatrixComplex &z);
   static bool       SMatrixTdEVD(double &d[],double &ce[],const int n,const int zneeded,CMatrixDouble &z);
   static bool       SMatrixTdEVDR(double &d[],double &e[],const int n,const int zneeded,const double a,const double b,int &m,CMatrixDouble &z);
   static bool       SMatrixTdEVDI(double &d[],double &e[],const int n,const int zneeded,const int i1,const int i2,CMatrixDouble &z);
   static bool       RMatrixEVD(CMatrixDouble &ca,const int n,const int vneeded,double &wr[],double &wi[],CMatrixDouble &vl,CMatrixDouble &vr);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CEigenVDetect::CEigenVDetect(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CEigenVDetect::~CEigenVDetect(void)
  {

  }
//+------------------------------------------------------------------+
//| Finding the eigenvalues and eigenvectors of a symmetric matrix   |
//| The algorithm finds eigen pairs of a symmetric matrix by reducing|
//| it to tridiagonal form and using the QL/QR algorithm.            |
//| Input parameters:                                                |
//|     A       -   symmetric matrix which is given by its upper or  |
//|                 lower triangular part.                           |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..N-1].                                |
//|     N       -   size of matrix A.                                |
//|     ZNeeded -   flag controlling whether the eigenvectors are    |
//|                 needed or not.                                   |
//|                 If ZNeeded is equal to:                          |
//|                  * 0, the eigenvectors are not returned;         |
//|                  * 1, the eigenvectors are returned.             |
//|     IsUpper -   storage format.                                  |
//| Output parameters:                                               |
//|     D       -   eigenvalues in ascending order.                  |
//|                 Array whose index ranges within [0..N-1].        |
//|     Z       -   if ZNeeded is equal to:                          |
//|                  * 0, Z hasn?t changed;                          |
//|                  * 1, Z contains the eigenvectors.               |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..N-1].                                |
//|                 The eigenvectors are stored in the matrix        |
//|                 columns.                                         |
//| Result:                                                          |
//|     True, if the algorithm has converged.                        |
//|     False, if the algorithm hasn't converged (rare case).        |
//+------------------------------------------------------------------+
static bool CEigenVDetect::SMatrixEVD(CMatrixDouble &ca,const int n,const int zneeded,
                                      const bool isupper,double &d[],CMatrixDouble &z)
  {
//--- create arrays
   double tau[];
   double e[];
//--- create copy
   CMatrixDouble a;
   a=ca;
//--- check
   if(!CAp::Assert(zneeded==0 || zneeded==1,__FUNCTION__+": incorrect ZNeeded"))
      return(false);
//--- function call
   COrtFac::SMatrixTD(a,n,isupper,tau,d,e);
//--- check
   if(zneeded==1)
     {
      //--- function call
      COrtFac::SMatrixTDUnpackQ(a,n,isupper,tau,z);
     }
//--- return result
   return(SMatrixTdEVD(d,e,n,zneeded,z));
  }
//+------------------------------------------------------------------+
//| Subroutine for finding the eigenvalues (and eigenvectors) of a   |
//| symmetric matrix in a given half open interval (A, B] by using a |
//| bisection and inverse iteration                                  |
//| Input parameters:                                                |
//|     A       -   symmetric matrix which is given by its upper or  |
//|                 lower triangular part. Array [0..N-1, 0..N-1].   |
//|     N       -   size of matrix A.                                |
//|     ZNeeded -   flag controlling whether the eigenvectors are    |
//|                 needed or not.                                   |
//|                 If ZNeeded is equal to:                          |
//|                  * 0, the eigenvectors are not returned;         |
//|                  * 1, the eigenvectors are returned.             |
//|     IsUpperA -  storage format of matrix A.                      |
//|     B1, B2 -    half open interval (B1, B2] to search            |
//|                 eigenvalues in.                                  |
//| Output parameters:                                               |
//|     M       -   number of eigenvalues found in a given           |
//|                 half-interval (M>=0).                            |
//|     W       -   array of the eigenvalues found.                  |
//|                 Array whose index ranges within [0..M-1].        |
//|     Z       -   if ZNeeded is equal to:                          |
//|                  * 0, Z hasn?t changed;                          |
//|                  * 1, Z contains eigenvectors.                   |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..M-1].                                |
//|                 The eigenvectors are stored in the matrix        |
//|                 columns.                                         |
//| Result:                                                          |
//|     True, if successful. M contains the number of eigenvalues in |
//|     the given half-interval (could be equal to 0), W contains the|
//|     eigenvalues, Z contains the eigenvectors (if needed).        |
//|     False, if the bisection method subroutine wasn't able to find|
//|     the eigenvalues in the given interval or if the inverse      |
//|     iteration subroutine wasn't able to find all the             |
//|     corresponding eigenvectors. In that case, the eigenvalues    |
//|     and eigenvectors are not returned, M is equal to 0.          |
//+------------------------------------------------------------------+
static bool CEigenVDetect::SMatrixEVDR(CMatrixDouble &ca,const int n,const int zneeded,
                                       const bool isupper,const double b1,const double b2,
                                       int &m,double &w[],CMatrixDouble &z)
  {
//--- create arrays
   double tau[];
   double e[];
//--- create copy
   CMatrixDouble a;
   a=ca;
//--- initialization
   m=0;
//--- check
   if(!CAp::Assert(zneeded==0 || zneeded==1,__FUNCTION__+": incorrect ZNeeded"))
      return(false);
//--- function call
   COrtFac::SMatrixTD(a,n,isupper,tau,w,e);
//--- check
   if(zneeded==1)
     {
      //--- function call
      COrtFac::SMatrixTDUnpackQ(a,n,isupper,tau,z);
     }
//--- return result
   return(SMatrixTdEVDR(w,e,n,zneeded,b1,b2,m,z));
  }
//+------------------------------------------------------------------+
//| Subroutine for finding the eigenvalues and eigenvectors of a     |
//| symmetric matrix with given indexes by using bisection and       |
//| inverse iteration methods.                                       |
//| Input parameters:                                                |
//|     A       -   symmetric matrix which is given by its upper or  |
//|                 lower triangular part. Array whose indexes range |
//|                 within [0..N-1, 0..N-1].                         |
//|     N       -   size of matrix A.                                |
//|     ZNeeded -   flag controlling whether the eigenvectors are    |
//|                 needed or not.                                   |
//|                 If ZNeeded is equal to:                          |
//|                  * 0, the eigenvectors are not returned;         |
//|                  * 1, the eigenvectors are returned.             |
//|     IsUpperA -  storage format of matrix A.                      |
//|     I1, I2 -    index interval for searching (from I1 to I2).    |
//|                 0 <= I1 <= I2 <= N-1.                            |
//| Output parameters:                                               |
//|     W       -   array of the eigenvalues found.                  |
//|                 Array whose index ranges within [0..I2-I1].      |
//|     Z       -   if ZNeeded is equal to:                          |
//|                  * 0, Z hasn?t changed;                          |
//|                  * 1, Z contains eigenvectors.                   |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..I2-I1].                              |
//|                 In that case, the eigenvectors are stored in the |
//|                 matrix columns.                                  |
//| Result:                                                          |
//|     True, if successful. W contains the eigenvalues, Z contains  |
//|     the eigenvectors (if needed).                                |
//|     False, if the bisection method subroutine wasn't able to find|
//|     the eigenvalues in the given interval or if the inverse      |
//|     iteration subroutine wasn't able to find all the             |
//|     corresponding eigenvectors. In that case, the eigenvalues    |
//|     and eigenvectors are not returned.                           |
//+------------------------------------------------------------------+
static bool CEigenVDetect::SMatrixEVDI(CMatrixDouble &ca,const int n,const int zneeded,
                                       const bool isupper,const int i1,const int i2,
                                       double &w[],CMatrixDouble &z)
  {
//--- create arrays
   double tau[];
   double e[];
//--- create copy
   CMatrixDouble a;
   a=ca;
//--- check
   if(!CAp::Assert(zneeded==0 || zneeded==1,__FUNCTION__+": incorrect ZNeeded"))
      return(false);
//--- function call
   COrtFac::SMatrixTD(a,n,isupper,tau,w,e);
//--- check
   if(zneeded==1)
     {
      //--- function call
      COrtFac::SMatrixTDUnpackQ(a,n,isupper,tau,z);
     }
//--- return result
   return(SMatrixTdEVDI(w,e,n,zneeded,i1,i2,z));
  }
//+------------------------------------------------------------------+
//| Finding the eigenvalues and eigenvectors of a Hermitian matrix   |
//| The algorithm finds eigen pairs of a Hermitian matrix by reducing|
//| it to real tridiagonal form and using the QL/QR algorithm.       |
//| Input parameters:                                                |
//|     A       -   Hermitian matrix which is given by its upper or  |
//|                 lower triangular part.                           |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..N-1].                                |
//|     N       -   size of matrix A.                                |
//|     IsUpper -   storage format.                                  |
//|     ZNeeded -   flag controlling whether the eigenvectors are    |
//|                 needed or not. If ZNeeded is equal to:           |
//|                  * 0, the eigenvectors are not returned;         |
//|                  * 1, the eigenvectors are returned.             |
//| Output parameters:                                               |
//|     D       -   eigenvalues in ascending order.                  |
//|                 Array whose index ranges within [0..N-1].        |
//|     Z       -   if ZNeeded is equal to:                          |
//|                  * 0, Z hasn?t changed;                          |
//|                  * 1, Z contains the eigenvectors.               |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..N-1].                                |
//|                 The eigenvectors are stored in the matrix        |
//|                 columns.                                         |
//| Result:                                                          |
//|     True, if the algorithm has converged.                        |
//|     False, if the algorithm hasn't converged (rare case).        |
//| Note:                                                            |
//|     eigenvectors of Hermitian matrix are defined up to           |
//|     multiplication by a complex number L, such that |L|=1.       |
//+------------------------------------------------------------------+
static bool CEigenVDetect::HMatrixEVD(CMatrixComplex &ca,const int n,int zneeded,
                                      const bool isupper,double &d[],CMatrixComplex &z)
  {
//--- create variables
   int    i=0;
   int    k=0;
   double v=0;
   int    i_=0;
   bool   result;
//--- create arrays
   al_complex tau[];
   double     e[];
   double     work[];
//--- create matrix
   CMatrixDouble  t;
   CMatrixComplex q;
//--- create copy
   CMatrixComplex a;
   a=ca;
//--- check
   if(!CAp::Assert(zneeded==0 || zneeded==1,__FUNCTION__+": incorrect ZNeeded"))
      return(false);
//--- Reduce to tridiagonal form
   COrtFac::HMatrixTD(a,n,isupper,tau,d,e);
//--- check
   if(zneeded==1)
     {
      //--- function call
      COrtFac::HMatrixTDUnpackQ(a,n,isupper,tau,q);
      zneeded=2;
     }
//--- get result
   result=SMatrixTdEVD(d,e,n,zneeded,t);
//--- Eigenvectors are needed
//--- Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
   if(result && zneeded!=0)
     {
      ArrayResizeAL(work,n);
      z.Resize(n,n);
      for(i=0;i<n;i++)
        {
         //--- Calculate real part
         for(k=0;k<n;k++)
            work[k]=0;
         for(k=0;k<n;k++)
           {
            v=q[i][k].re;
            for(i_=0;i_<n;i_++)
               work[i_]=work[i_]+v*t[k][i_];
           }
         //--- get real part
         for(k=0;k<n;k++)
            z[i].SetRe(k,work[k]);
         //--- Calculate imaginary part
         for(k=0;k<n;k++)
            work[k]=0;
         for(k=0;k<n;k++)
           {
            v=q[i][k].im;
            for(i_=0;i_<n;i_++)
               work[i_]=work[i_]+v*t[k][i_];
           }
         //--- get imaginary part
         for(k=0;k<n;k++)
            z[i].SetIm(k,work[k]);
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Subroutine for finding the eigenvalues (and eigenvectors) of a   |
//| Hermitian matrix in a given half-interval (A, B] by using a      |
//| bisection and inverse iteration                                  |
//| Input parameters:                                                |
//|     A       -   Hermitian matrix which is given by its upper or  |
//|                 lower triangular part. Array whose indexes range |
//|                 within [0..N-1, 0..N-1].                         |
//|     N       -   size of matrix A.                                |
//|     ZNeeded -   flag controlling whether the eigenvectors are    |
//|                needed or not. If ZNeeded is equal to:            |
//|                  * 0, the eigenvectors are not returned;         |
//|                  * 1, the eigenvectors are returned.             |
//|     IsUpperA -  storage format of matrix A.                      |
//|     B1, B2 -    half-interval (B1, B2] to search eigenvalues in. |
//| Output parameters:                                               |
//|     M       -   number of eigenvalues found in a given           |
//|                 half-interval, M>=0                              |
//|     W       -   array of the eigenvalues found.                  |
//|                 Array whose index ranges within [0..M-1].        |
//|     Z       -   if ZNeeded is equal to:                          |
//|                  * 0, Z hasn?t changed;                          |
//|                  * 1, Z contains eigenvectors.                   |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..M-1].                                |
//|                 The eigenvectors are stored in the matrix        |
//|                 columns.                                         |
//| Result:                                                          |
//|     True, if successful. M contains the number of eigenvalues    |
//|     in the given half-interval (could be equal to 0), W contains |
//|     the eigenvalues, Z contains the eigenvectors (if needed).    |
//|     False, if the bisection method subroutine wasn't able to find|
//|     the eigenvalues in the given interval or if the inverse      |
//|     iteration subroutine wasn't able to find all the             |
//|     corresponding eigenvectors. In that case, the eigenvalues and|
//|     eigenvectors are not returned, M is equal to 0.              |
//| Note:                                                            |
//|     eigen vectors of Hermitian matrix are defined up to          |
//|     multiplication by a complex number L, such as |L|=1.         |
//+------------------------------------------------------------------+
static bool CEigenVDetect::HMatrixEVDR(CMatrixComplex &ca,const int n,int zneeded,
                                       bool isupper,const double b1,const double b2,
                                       int &m,double &w[],CMatrixComplex &z)
  {
//--- create variables
   int    i=0;
   int    k=0;
   double v=0;
   int    i_=0;
   bool   result;
//--- create arrays
   al_complex tau[];
   double     e[];
   double     work[];
//--- create matrix
   CMatrixComplex q;
   CMatrixDouble  t;
//--- create copy
   CMatrixComplex a;
   a=ca;
//--- initialization
   m=0;
//--- check
   if(!CAp::Assert(zneeded==0 || zneeded==1,__FUNCTION__+": incorrect ZNeeded"))
      return(false);
//--- Reduce to tridiagonal form
   COrtFac::HMatrixTD(a,n,isupper,tau,w,e);
//--- check
   if(zneeded==1)
     {
      //--- function call
      COrtFac::HMatrixTDUnpackQ(a,n,isupper,tau,q);
      zneeded=2;
     }
//--- Bisection and inverse iteration
   result=SMatrixTdEVDR(w,e,n,zneeded,b1,b2,m,t);
//--- Eigenvectors are needed
//--- Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
   if((result && zneeded!=0) && m!=0)
     {
      ArrayResizeAL(work,m);
      z.Resize(n,m);
      for(i=0;i<n;i++)
        {
         //--- Calculate real part
         for(k=0;k<=m-1;k++)
            work[k]=0;
         for(k=0;k<n;k++)
           {
            v=q[i][k].re;
            for(i_=0;i_<m;i_++)
               work[i_]=work[i_]+v*t[k][i_];
           }
         //--- get real part
         for(k=0;k<=m-1;k++)
            z[i].SetRe(k,work[k]);
         //--- Calculate imaginary part
         for(k=0;k<=m-1;k++)
            work[k]=0;
         for(k=0;k<n;k++)
           {
            v=q[i][k].im;
            for(i_=0;i_<m;i_++)
               work[i_]=work[i_]+v*t[k][i_];
           }
         //--- get imaginary part
         for(k=0;k<=m-1;k++)
            z[i].SetIm(k,work[k]);
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Subroutine for finding the eigenvalues and eigenvectors of a     |
//| Hermitian matrix with given indexes by using bisection and       |
//| inverse iteration methods                                        |
//| Input parameters:                                                |
//|     A       -   Hermitian matrix which is given by its upper or  |
//|                 lower triangular part.                           |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..N-1].                                |
//|     N       -   size of matrix A.                                |
//|     ZNeeded -   flag controlling whether the eigenvectors are    |
//|                 needed or not. If ZNeeded is equal to:           |
//|                  * 0, the eigenvectors are not returned;         |
//|                  * 1, the eigenvectors are returned.             |
//|     IsUpperA -  storage format of matrix A.                      |
//|     I1, I2 -    index interval for searching (from I1 to I2).    |
//|                 0 <= I1 <= I2 <= N-1.                            |
//| Output parameters:                                               |
//|     W       -   array of the eigenvalues found.                  |
//|                 Array whose index ranges within [0..I2-I1].      |
//|     Z       -   if ZNeeded is equal to:                          |
//|                  * 0, Z hasn?t changed;                          |
//|                  * 1, Z contains eigenvectors.                   |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..I2-I1].                              |
//|                 In  that  case,  the eigenvectors are stored in  |
//|                 the matrix columns.                              |
//| Result:                                                          |
//|     True, if successful. W contains the eigenvalues, Z contains  |
//|     the eigenvectors (if needed).                                |
//|     False, if the bisection method subroutine wasn't able to find|
//|     the eigenvalues in the given interval or if the inverse      |
//|     corresponding eigenvectors. iteration subroutine wasn't able |
//|     to find all the corresponding  eigenvectors. In that case,   |
//|     the eigenvalues and  eigenvectors are not returned.          |
//| Note:                                                            |
//|     eigen vectors of Hermitian matrix are defined up to          |
//|     multiplication  by a complex number L, such as |L|=1.        |
//+------------------------------------------------------------------+
static bool CEigenVDetect::HMatrixEVDI(CMatrixComplex &ca,const int n,int zneeded,
                                       const bool isupper,const int i1,const int i2,
                                       double &w[],CMatrixComplex &z)
  {
//--- create variables
   int    i=0;
   int    k=0;
   double v=0;
   int    m=0;
   int    i_=0;
   bool   result;
//--- create arrays
   al_complex tau[];
   double     e[];
   double     work[];
//--- create matrix
   CMatrixComplex q;
   CMatrixDouble t;
//--- create copy
   CMatrixComplex a;
   a=ca;
//--- check
   if(!CAp::Assert(zneeded==0 || zneeded==1,__FUNCTION__+": incorrect ZNeeded"))
      return(false);
//--- Reduce to tridiagonal form
   COrtFac::HMatrixTD(a,n,isupper,tau,w,e);
//--- check
   if(zneeded==1)
     {
      //--- function call
      COrtFac::HMatrixTDUnpackQ(a,n,isupper,tau,q);
      zneeded=2;
     }
//--- Bisection and inverse iteration
   result=SMatrixTdEVDI(w,e,n,zneeded,i1,i2,t);
//--- Eigenvectors are needed
//--- Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
   m=i2-i1+1;
//--- check
   if(result && zneeded!=0)
     {
      ArrayResizeAL(work,m);
      z.Resize(n,m);
      for(i=0;i<n;i++)
        {
         //--- Calculate real part
         for(k=0;k<=m-1;k++)
            work[k]=0;
         for(k=0;k<n;k++)
           {
            v=q[i][k].re;
            for(i_=0;i_<m;i_++)
               work[i_]=work[i_]+v*t[k][i_];
           }
         //--- get real part
         for(k=0;k<=m-1;k++)
            z[i].SetRe(k,work[k]);
         //--- Calculate imaginary part
         for(k=0;k<=m-1;k++)
            work[k]=0;
         for(k=0;k<n;k++)
           {
            v=q[i][k].im;
            for(i_=0;i_<m;i_++)
               work[i_]=work[i_]+v*t[k][i_];
           }
         //--- get imaginary part
         for(k=0;k<=m-1;k++)
            z[i].SetIm(k,work[k]);
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Finding the eigenvalues and eigenvectors of a tridiagonal        |
//| symmetric matrix                                                 |
//| The algorithm finds the eigen pairs of a tridiagonal symmetric   |
//| matrix by using an QL/QR algorithm with implicit shifts.         |
//| Input parameters:                                                |
//|     D       -   the main diagonal of a tridiagonal matrix.       |
//|                 Array whose index ranges within [0..N-1].        |
//|     E       -   the secondary diagonal of a tridiagonal matrix.  |
//|                 Array whose index ranges within [0..N-2].        |
//|     N       -   size of matrix A.                                |
//|     ZNeeded -   flag controlling whether the eigenvectors are    |
//|                 needed or not.                                   |
//|                 If ZNeeded is equal to:                          |
//|                  * 0, the eigenvectors are not needed;           |
//|                  * 1, the eigenvectors of a tridiagonal matrix   |
//|                    are multiplied by the square matrix Z. It is  |
//|                    used if the tridiagonal matrix is obtained by |
//|                    the similarity transformation of a symmetric  |
//|                    matrix;                                       |
//|                  * 2, the eigenvectors of a tridiagonal matrix   |
//|                    replace the square matrix Z;                  |
//|                  * 3, matrix Z contains the first row of the     |
//|                    eigenvectors matrix.                          |
//|     Z       -   if ZNeeded=1, Z contains the square matrix by    |
//|                 which the eigenvectors are multiplied.           |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..N-1].                                |
//| Output parameters:                                               |
//|     D       -   eigenvalues in ascending order.                  |
//|                 Array whose index ranges within [0..N-1].        |
//|     Z       -   if ZNeeded is equal to:                          |
//|                  * 0, Z hasn?t changed;                          |
//|                  * 1, Z contains the product of a given matrix   |
//|                    (from the left) and the eigenvectors matrix   |
//|                    (from the right);                             |
//|                  * 2, Z contains the eigenvectors.               |
//|                  * 3, Z contains the first row of the            |
//|                       eigenvectors matrix.                       |
//|                 If ZNeeded<3, Z is the array whose indexes range |
//|                 within [0..N-1, 0..N-1].                         |
//|                 In that case, the eigenvectors are stored in the |
//|                 matrix columns.                                  |
//|                 If ZNeeded=3, Z is the array whose indexes range |
//|                 within [0..0, 0..N-1].                           |
//| Result:                                                          |
//|     True, if the algorithm has converged.                        |
//|     False, if the algorithm hasn't converged.                    |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      September 30, 1994                                          |
//+------------------------------------------------------------------+
static bool CEigenVDetect::SMatrixTdEVD(double &d[],double &ce[],const int n,const int zneeded,CMatrixDouble &z)
  {
//--- create variables
   int  i=0;
   int  i_=0;
   int  i1_=0;
   bool result;
//--- create arrays
   double d1[];
   double e1[];
//--- create matrix
   CMatrixDouble z1;
//--- create copy
   double e[];
   ArrayResizeAL(e,ArraySize(ce));
   ArrayCopy(e,ce);
//--- Prepare 1-based task
   ArrayResizeAL(d1,n+1);
   ArrayResizeAL(e1,n+1);
   i1_=-1;
   for(i_=1;i_<=n;i_++)
      d1[i_]=d[i_+i1_];
//--- check
   if(n>1)
     {
      i1_=-1;
      for(i_=1;i_<n;i_++)
         e1[i_]=e[i_+i1_];
     }
//--- check
   if(zneeded==1)
     {
      z1.Resize(n+1,n+1);
      for(i=1;i<=n;i++)
        {
         i1_=-1;
         for(i_=1;i_<=n;i_++)
            z1[i].Set(i_,z[i-1][i_+i1_]);
        }
     }
//--- Solve 1-based task
   result=TtidiagonalEVD(d1,e1,n,zneeded,z1);
//--- check
   if(!result)
      return(result);
//--- Convert back to 0-based result
   i1_=1;
   for(i_=0;i_<n;i_++)
      d[i_]=d1[i_+i1_];
//--- check
   if(zneeded!=0)
     {
      //--- check
      if(zneeded==1)
        {
         for(i=1;i<=n;i++)
           {
            i1_=1;
            for(i_=0;i_<n;i_++)
               z[i-1].Set(i_,z1[i][i_+i1_]);
           }
         //--- return result
         return(result);
        }
      //--- check
      if(zneeded==2)
        {
         z.Resize(n,n);
         for(i=1;i<=n;i++)
           {
            i1_=1;
            for(i_=0;i_<n;i_++)
               z[i-1].Set(i_,z1[i][i_+i1_]);
           }
         //--- return result
         return(result);
        }
      //--- check
      if(zneeded==3)
        {
         z.Resize(1,n);
         i1_=1;
         for(i_=0;i_<n;i_++)
            z[0].Set(i_,z1[1][i_+i1_]);
         //--- return result
         return(result);
        }
      //--- check
      if(!CAp::Assert(false,__FUNCTION__+": Incorrect ZNeeded!"))
         return(false);
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Subroutine for finding the tridiagonal matrix eigenvalues/vectors|
//| in a given half-interval (A, B] by using bisection and inverse   |
//| iteration.                                                       |
//| Input parameters:                                                |
//|     D       -   the main diagonal of a tridiagonal matrix.       |
//|                 Array whose index ranges within [0..N-1].        |
//|     E       -   the secondary diagonal of a tridiagonal matrix.  |
//|                 Array whose index ranges within [0..N-2].        |
//|     N       -   size of matrix, N>=0.                            |
//|     ZNeeded -   flag controlling whether the eigenvectors are    |
//|                 needed or not. If ZNeeded is equal to:           |
//|                  * 0, the eigenvectors are not needed;           |
//|                  * 1, the eigenvectors of a tridiagonal matrix   |
//|                    are multiplied by the square matrix Z. It is  |
//|                    used if the tridiagonal matrix is obtained by |
//|                    the similarity transformation of a symmetric  |
//|                    matrix.                                       |
//|                  * 2, the eigenvectors of a tridiagonal matrix   |
//|                    replace matrix Z.                             |
//|     A, B    -   half-interval (A, B] to search eigenvalues in.   |
//|     Z       -   if ZNeeded is equal to:                          |
//|                  * 0, Z isn't used and remains unchanged;        |
//|                  * 1, Z contains the square matrix (array whose  |
//|                    indexes range within [0..N-1, 0..N-1]) which  |
//|                    reduces the given symmetric matrix to         |
//|                    tridiagonal form;                             |
//|                  * 2, Z isn't used (but changed on the exit).    |
//| Output parameters:                                               |
//|     D       -   array of the eigenvalues found.                  |
//|                 Array whose index ranges within [0..M-1].        |
//|     M       -   number of eigenvalues found in the given         |
//|                 half-interval (M>=0).                            |
//|     Z       -   if ZNeeded is equal to:                          |
//|                  * 0, doesn't contain any information;           |
//|                  * 1, contains the product of a given NxN matrix |
//|                    Z (from the left) and NxM matrix of the       |
//|                    eigenvectors found (from the right). Array    |
//|                    whose indexes range within [0..N-1, 0..M-1].  |
//|                  * 2, contains the matrix of the eigenvectors    |
//|                    found. Array whose indexes range within       |
//|                    [0..N-1, 0..M-1].                             |
//| Result:                                                          |
//|     True, if successful. In that case, M contains the number of  |
//|     eigenvalues in the given half-interval (could be equal to 0),|
//|     D contains the eigenvalues, Z contains the eigenvectors (if  |
//|     needed). It should be noted that the subroutine changes the  |
//|     size of arrays D and Z.                                      |
//|     False, if the bisection method subroutine wasn't able to find|
//|     the eigenvalues in the given interval or if the inverse      |
//|     iteration subroutine wasn't able to find all the             |
//|     corresponding eigenvectors. In that case, the eigenvalues and|
//|     eigenvectors are not returned, M is equal to 0.              |
//+------------------------------------------------------------------+
static bool CEigenVDetect::SMatrixTdEVDR(double &d[],double &e[],const int n,
                                         const int zneeded,const double a,
                                         const double b,int &m,CMatrixDouble &z)
  {
//--- create variables
   bool   result;
   int    errorcode=0;
   int    nsplit=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    cr=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- create arrays
   int    iblock[];
   int    isplit[];
   int    ifail[];
   double d1[];
   double e1[];
   double w[];
//--- create matrix
   CMatrixDouble z2;
   CMatrixDouble z3;
//--- initialization
   m=0;
//--- check
   if(!CAp::Assert(zneeded>=0 && zneeded<=2,__FUNCTION__+": incorrect ZNeeded!"))
      return(false);
//--- check
   if(b<=a)
     {
      m=0;
      //--- return result
      return(true);
     }
//--- check
   if(n<=0)
     {
      m=0;
      //--- return result
      return(true);
     }
//--- Copy D,E to D1, E1
   ArrayResizeAL(d1,n+1);
   i1_=-1;
   for(i_=1;i_<=n;i_++)
      d1[i_]=d[i_+i1_];
//--- check
   if(n>1)
     {
      ArrayResizeAL(e1,n);
      i1_=-1;
      for(i_=1;i_<n;i_++)
         e1[i_]=e[i_+i1_];
     }
//--- No eigen vectors
   if(zneeded==0)
     {
      //--- get result
      result=InternalBisectionEigenValues(d1,e1,n,2,1,a,b,0,0,-1,w,m,nsplit,iblock,isplit,errorcode);
      //--- check
      if(!result || m==0)
        {
         m=0;
         return(result);
        }
      ArrayResizeAL(d,m);
      i1_=1;
      for(i_=0;i_<m;i_++)
        {
         d[i_]=w[i_+i1_];
        }
      //--- return result
      return(result);
     }
//--- Eigen vectors are multiplied by Z
   if(zneeded==1)
     {
      //--- Find eigen pairs
      result=InternalBisectionEigenValues(d1,e1,n,2,2,a,b,0,0,-1,w,m,nsplit,iblock,isplit,errorcode);
      //--- check
      if(!result || m==0)
        {
         m=0;
         return(result);
        }
      //--- function call
      InternalDStein(n,d1,e1,m,w,iblock,isplit,z2,ifail,cr);
      //--- check
      if(cr!=0)
        {
         m=0;
         //--- get result
         result=false;
         return(result);
        }
      //--- Sort eigen values and vectors
      for(i=1;i<=m;i++)
        {
         k=i;
         for(j=i;j<=m;j++)
           {
            //--- check
            if(w[j]<w[k])
               k=j;
           }
         //--- swap
         v=w[i];
         w[i]=w[k];
         w[k]=v;
         for(j=1;j<=n;j++)
           {
            //--- swap
            v=z2[j][i];
            z2[j].Set(i,z2[j][k]);
            z2[j].Set(k,v);
           }
        }
      //--- Transform Z2 and overwrite Z
      z3.Resize(m+1,n+1);
      for(i=1;i<=m;i++)
         for(i_=1;i_<=n;i_++)
            z3[i].Set(i_,z2[i_][i]);
      for(i=1;i<=n;i++)
        {
         for(j=1;j<=m;j++)
           {
            i1_=1;
            v=0.0;
            for(i_=0;i_<n;i_++)
               v+=z[i-1][i_]*z3[j][i_+i1_];
            z2[i].Set(j,v);
           }
        }
      //--- rewrite
      z.Resize(n,m);
      for(i=1;i<=m;i++)
        {
         i1_=1;
         for(i_=0;i_<n;i_++)
            z[i_].Set(i-1,z2[i_+i1_][i]);
        }
      //--- Store W
      ArrayResizeAL(d,m);
      for(i=1;i<=m;i++)
         d[i-1]=w[i];
      //--- return result
      return(result);
     }
//--- Eigen vectors are stored in Z
   if(zneeded==2)
     {
      //--- Find eigen pairs
      result=InternalBisectionEigenValues(d1,e1,n,2,2,a,b,0,0,-1,w,m,nsplit,iblock,isplit,errorcode);
      //--- check
      if(!result || m==0)
        {
         m=0;
         return(result);
        }
      //--- function call
      InternalDStein(n,d1,e1,m,w,iblock,isplit,z2,ifail,cr);
      //--- check
      if(cr!=0)
        {
         m=0;
         result=false;
         return(result);
        }
      //--- Sort eigen values and vectors
      for(i=1;i<=m;i++)
        {
         k=i;
         for(j=i;j<=m;j++)
           {
            //--- check
            if(w[j]<w[k])
               k=j;
           }
         //--- swap
         v=w[i];
         w[i]=w[k];
         w[k]=v;
         for(j=1;j<=n;j++)
           {
            //--- swap
            v=z2[j][i];
            z2[j].Set(i,z2[j][k]);
            z2[j].Set(k,v);
           }
        }
      //--- Store W
      ArrayResizeAL(d,m);
      for(i=1;i<=m;i++)
         d[i-1]=w[i];
      z.Resize(n,m);
      for(i=1;i<=m;i++)
        {
         i1_=1;
         for(i_=0;i_<n;i_++)
            z[i_].Set(i-1,z2[i_+i1_][i]);
        }
      //--- return result
      return(result);
     }
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Subroutine for finding tridiagonal matrix eigenvalues/vectors    |
//| with given indexes (in ascending order) by using the bisection   |
//| and inverse iteraion.                                            |
//| Input parameters:                                                |
//|     D       -   the main diagonal of a tridiagonal matrix.       |
//|                 Array whose index ranges within [0..N-1].        |
//|     E       -   the secondary diagonal of a tridiagonal matrix.  |
//|                 Array whose index ranges within [0..N-2].        |
//|     N       -   size of matrix. N>=0.                            |
//|     ZNeeded -   flag controlling whether the eigenvectors are    |
//|                 needed or not. If ZNeeded is equal to:           |
//|                  * 0, the eigenvectors are not needed;           |
//|                  * 1, the eigenvectors of a tridiagonal matrix   |
//|                    are multiplied by the square matrix Z. It is  |
//|                    used if the tridiagonal matrix is obtained by |
//|                    the similarity transformation of a symmetric  |
//|                    matrix.                                       |
//|                  * 2, the eigenvectors of a tridiagonal matrix   |
//|                    replace matrix Z.                             |
//|     I1, I2  -   index interval for searching (from I1 to I2).    |
//|                 0 <= I1 <= I2 <= N-1.                            |
//|     Z       -   if ZNeeded is equal to:                          |
//|                  * 0, Z isn't used and remains unchanged;        |
//|                  * 1, Z contains the square matrix (array whose  |
//|                    indexes range within [0..N-1, 0..N-1]) which  |
//|                    reduces the given symmetric matrix to         |
//|                    tridiagonal form;                             |
//|                  * 2, Z isn't used (but changed on the exit).    |
//| Output parameters:                                               |
//|     D       -   array of the eigenvalues found.                  |
//|                 Array whose index ranges within [0..I2-I1].      |
//|     Z       -   if ZNeeded is equal to:                          |
//|                  * 0, doesn't contain any information;           |
//|                  * 1, contains the product of a given NxN matrix |
//|                    Z (from the left) and Nx(I2-I1) matrix of the |
//|                    eigenvectors found (from the right). Array    |
//|                    whose indexes range within [0..N-1, 0..I2-I1].|
//|                  * 2, contains the matrix of the eigenvalues     |
//|                    found. Array whose indexes range within       |
//|                    [0..N-1, 0..I2-I1].                           |
//| Result:                                                          |
//|     True, if successful. In that case, D contains the            |
//|     eigenvalues, Z contains the eigenvectors (if needed).        |
//|     It should be noted that the subroutine changes the size of   |
//|     arrays D and Z.                                              |
//|     False, if the bisection method subroutine wasn't able to find|
//|     the eigenvalues in the given interval or if the inverse      |
//|     iteration subroutine wasn't able to find all the             |
//|     corresponding eigenvectors. In that case, the eigenvalues and|
//|     eigenvectors are not returned.                               |
//+------------------------------------------------------------------+
static bool CEigenVDetect::SMatrixTdEVDI(double &d[],double &e[],const int n,
                                         const int zneeded,const int i1,
                                         const int i2,CMatrixDouble &z)
  {
//--- create variables
   bool   result;
   int    errorcode=0;
   int    nsplit=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    m=0;
   int    cr=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- create arrays
   int    iblock[];
   int    isplit[];
   int    ifail[];
   double w[];
   double d1[];
   double e1[];
//--- create matrix
   CMatrixDouble z2;
   CMatrixDouble z3;
//--- check
   if(!CAp::Assert((0<=i1 && i1<=i2) && i2<n,__FUNCTION__+": incorrect I1/I2!"))
      return(false);
//--- Copy D,E to D1, E1
   ArrayResizeAL(d1,n+1);
   i1_=-1;
   for(i_=1;i_<=n;i_++)
      d1[i_]=d[i_+i1_];
//--- check
   if(n>1)
     {
      ArrayResizeAL(e1,n);
      i1_=-1;
      for(i_=1;i_<n;i_++)
         e1[i_]=e[i_+i1_];
     }
//--- No eigen vectors
   if(zneeded==0)
     {
      result=InternalBisectionEigenValues(d1,e1,n,3,1,0,0,i1+1,i2+1,-1,w,m,nsplit,iblock,isplit,errorcode);
      //--- check
      if(!result)
         return(result);
      //--- check
      if(m!=i2-i1+1)
        {
         result=false;
         return(result);
        }
      ArrayResizeAL(d,m);
      for(i=1;i<=m;i++)
         d[i-1]=w[i];
      //--- return result
      return(result);
     }
//--- Eigen vectors are multiplied by Z
   if(zneeded==1)
     {
      //--- Find eigen pairs
      result=InternalBisectionEigenValues(d1,e1,n,3,2,0,0,i1+1,i2+1,-1,w,m,nsplit,iblock,isplit,errorcode);
      //--- check
      if(!result)
         return(result);
      //--- check
      if(m!=i2-i1+1)
        {
         result=false;
         return(result);
        }
      //--- function call
      InternalDStein(n,d1,e1,m,w,iblock,isplit,z2,ifail,cr);
      //--- check
      if(cr!=0)
        {
         result=false;
         return(result);
        }
      //--- Sort eigen values and vectors
      for(i=1;i<=m;i++)
        {
         k=i;
         for(j=i;j<=m;j++)
           {
            //--- check
            if(w[j]<w[k])
               k=j;
           }
         //--- swap
         v=w[i];
         w[i]=w[k];
         w[k]=v;
         for(j=1;j<=n;j++)
           {
            //--- swap
            v=z2[j][i];
            z2[j].Set(i,z2[j][k]);
            z2[j].Set(k,v);
           }
        }
      //--- Transform Z2 and overwrite Z
      z3.Resize(m+1,n+1);
      for(i=1;i<=m;i++)
         for(i_=1;i_<=n;i_++)
            z3[i].Set(i_,z2[i_][i]);
      for(i=1;i<=n;i++)
        {
         for(j=1;j<=m;j++)
           {
            i1_=1;
            v=0.0;
            for(i_=0;i_<n;i_++)
               v+=z[i-1][i_]*z3[j][i_+i1_];
            z2[i].Set(j,v);
           }
        }
      //--- rewrite z
      z.Resize(n,m);
      for(i=1;i<=m;i++)
        {
         i1_=1;
         for(i_=0;i_<n;i_++)
            z[i_].Set(i-1,z2[i_+i1_][i]);
        }
      //--- Store W
      ArrayResizeAL(d,m);
      for(i=1;i<=m;i++)
         d[i-1]=w[i];
      //--- return result
      return(result);
     }
//--- Eigen vectors are stored in Z
   if(zneeded==2)
     {
      //--- Find eigen pairs
      result=InternalBisectionEigenValues(d1,e1,n,3,2,0,0,i1+1,i2+1,-1,w,m,nsplit,iblock,isplit,errorcode);
      //--- check
      if(!result)
         return(result);
      //--- check
      if(m!=i2-i1+1)
        {
         result=false;
         return(result);
        }
      //--- function call
      InternalDStein(n,d1,e1,m,w,iblock,isplit,z2,ifail,cr);
      //--- check
      if(cr!=0)
        {
         result=false;
         return(result);
        }
      //--- Sort eigen values and vectors
      for(i=1;i<=m;i++)
        {
         k=i;
         for(j=i;j<=m;j++)
           {
            //--- check
            if(w[j]<w[k])
               k=j;
           }
         //--- swap
         v=w[i];
         w[i]=w[k];
         w[k]=v;
         for(j=1;j<=n;j++)
           {
            //--- swap
            v=z2[j][i];
            z2[j].Set(i,z2[j][k]);
            z2[j].Set(k,v);
           }
        }
      //--- Store Z
      z.Resize(n,m);
      for(i=1;i<=m;i++)
        {
         i1_=1;
         for(i_=0;i_<n;i_++)
            z[i_].Set(i-1,z2[i_+i1_][i]);
        }
      //--- Store W
      ArrayResizeAL(d,m);
      for(i=1;i<=m;i++)
         d[i-1]=w[i];
      //--- return result
      return(result);
     }
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| Finding eigenvalues and eigenvectors of a general matrix         |
//| The algorithm finds eigenvalues and eigenvectors of a general    |
//| matrix by using the QR algorithm with multiple shifts. The       |
//| algorithm can find eigenvalues and both left and right           |
//| eigenvectors.                                                    |
//| The right eigenvector is a vector x such that A*x = w*x, and the |
//| left eigenvector is a vector y such that y'*A = w*y' (here y'    |
//| implies a complex conjugate transposition of vector y).          |
//| Input parameters:                                                |
//|     A       -   matrix. Array whose indexes range within         |
//|                 [0..N-1, 0..N-1].                                |
//|     N       -   size of matrix A.                                |
//|     VNeeded -   flag controlling whether eigenvectors are needed |
//|                 or not. If VNeeded is equal to:                  |
//|                  * 0, eigenvectors are not returned;             |
//|                  * 1, right eigenvectors are returned;           |
//|                  * 2, left eigenvectors are returned;            |
//|                  * 3, both left and right eigenvectors are       |
//|                       returned.                                  |
//| Output parameters:                                               |
//|     WR      -   real parts of eigenvalues.                       |
//|                 Array whose index ranges within [0..N-1].        |
//|     WR      -   imaginary parts of eigenvalues.                  |
//|                 Array whose index ranges within [0..N-1].        |
//|     VL, VR  -   arrays of left and right eigenvectors (if they   |
//|                 are needed). If WI[i]=0, the respective          |
//|                 eigenvalue is a real number, and it corresponds  |
//|                 to the column number I of matrices VL/VR. If     |
//|                 WI[i]>0, we have a pair of complex conjugate     |
//|                     numbers with positive and negative imaginary |
//|                     parts: the first eigenvalue WR[i] +          |
//|                     + sqrt(-1)*WI[i]; the second eigenvalue      |
//|                     WR[i+1] + sqrt(-1)*WI[i+1];                  |
//|                     WI[i]>0                                      |
//|                     WI[i+1] = -WI[i] < 0                         |
//|                 In that case, the eigenvector  corresponding to  |
//|                 the first eigenvalue is located in i and i+1     |
//|                 columns of matrices VL/VR (the column number i   |
//|                 contains the real part, and the column number    |
//|                 i+1 contains the imaginary part), and the vector |
//|                 corresponding to the second eigenvalue is a      |
//|                 complex conjugate to the first vector.           |
//|                 Arrays whose indexes range within                |
//|                 [0..N-1, 0..N-1].                                |
//| Result:                                                          |
//|     True, if the algorithm has converged.                        |
//|     False, if the algorithm has not converged.                   |
//| Note 1:                                                          |
//|     Some users may ask the following question: what if WI[N-1]>0?|
//|     WI[N] must contain an eigenvalue which is complex conjugate  |
//|     to the N-th eigenvalue, but the array has only size N?       |
//|     The answer is as follows: such a situation cannot occur      |
//|     because the algorithm finds a pairs of eigenvalues,          |
//|     therefore, if WI[i]>0, I is strictly less than N-1.          |
//| Note 2:                                                          |
//|     The algorithm performance depends on the value of the        |
//|     internal parameter NS of the InternalSchurDecomposition      |
//|     subroutine which defines the number of shifts in the QR      |
//|     algorithm (similarly to the block width in block-matrix      |
//|     algorithms of linear algebra). If you require maximum        |
//|     performance on your machine, it is recommended to adjust     |
//|     this parameter manually.                                     |
//| See also the InternalTREVC subroutine.                           |
//| The algorithm is based on the LAPACK 3.0 library.                |
//+------------------------------------------------------------------+
static bool CEigenVDetect::RMatrixEVD(CMatrixDouble &ca,const int n,const int vneeded,
                                      double &wr[],double &wi[],
                                      CMatrixDouble &vl,CMatrixDouble &vr)
  {
//--- create variables
   int  i=0;
   int  i_=0;
   int  i1_=0;
   bool result;
//--- create arrays
   double wr1[];
   double wi1[];
//--- create matrix
   CMatrixDouble a1;
   CMatrixDouble vl1;
   CMatrixDouble vr1;
//--- create copy
   CMatrixDouble a;
   a=ca;
//--- check
   if(!CAp::Assert(vneeded>=0 && vneeded<=3,__FUNCTION__+": incorrect VNeeded!"))
      return(false);
   a1.Resize(n+1,n+1);
   for(i=1;i<=n;i++)
     {
      i1_=-1;
      for(i_=1;i_<=n;i_++)
         a1[i].Set(i_,a[i-1][i_+i1_]);
     }
//--- get result
   result=NonSymmetricEVD(a1,n,vneeded,wr1,wi1,vl1,vr1);
//--- check
   if(result)
     {
      //--- allocation
      ArrayResizeAL(wr,n);
      ArrayResizeAL(wi,n);
      i1_=1;
      for(i_=0;i_<n;i_++)
         wr[i_]=wr1[i_+i1_];
      i1_=1;
      for(i_=0;i_<n;i_++)
         wi[i_]=wi1[i_+i1_];
      //--- check
      if(vneeded==2 || vneeded==3)
        {
         vl.Resize(n,n);
         for(i=0;i<n;i++)
           {
            i1_=1;
            for(i_=0;i_<n;i_++)
               vl[i].Set(i_,vl1[i+1][i_+i1_]);
           }
        }
      //--- check
      if(vneeded==1 || vneeded==3)
        {
         vr.Resize(n,n);
         for(i=0;i<n;i++)
           {
            i1_=1;
            for(i_=0;i_<n;i_++)
               vr[i].Set(i_,vr1[i+1][i_+i1_]);
           }
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Eigenvalues and eigenvectors in tridiagonal matrix               |
//+------------------------------------------------------------------+
static bool CEigenVDetect::TtidiagonalEVD(double &d[],double &ce[],const int n,
                                          const int zneeded,CMatrixDouble &z)
  {
//--- create variables
   bool   result;
   int    maxit=0;
   int    i=0;
   int    ii=0;
   int    iscale=0;
   int    j=0;
   int    jtot=0;
   int    k=0;
   int    t=0;
   int    l=0;
   int    l1=0;
   int    lend=0;
   int    lendm1=0;
   int    lendp1=0;
   int    lendsv=0;
   int    lm1=0;
   int    lsv=0;
   int    m=0;
   int    mm=0;
   int    mm1=0;
   int    nm1=0;
   int    nmaxit=0;
   int    tmpint=0;
   double anorm=0;
   double b=0;
   double c=0;
   double eps=0;
   double eps2=0;
   double f=0;
   double g=0;
   double p=0;
   double r=0;
   double rt1=0;
   double rt2=0;
   double s=0;
   double safmax=0;
   double safmin=0;
   double ssfmax=0;
   double ssfmin=0;
   double tst=0;
   double tmp=0;
   bool   gotoflag;
   int    zrows=0;
   bool   wastranspose;
   int    i_=0;
//--- create arrays
   double work1[];
   double work2[];
   double workc[];
   double works[];
   double wtemp[];
//--- create copy
   double e[];
   ArrayResizeAL(e,ArraySize(ce));
   ArrayCopy(e,ce);
//--- check
   if(!CAp::Assert(zneeded>=0 && zneeded<=3,__FUNCTION__+": Incorrent ZNeeded"))
      return(false);
//--- check
   if(zneeded<0 || zneeded>3)
      return(false);
//--- initialization
   result=true;
//--- check
   if(n==0)
      return(result);
//--- check
   if(n==1)
     {
      //--- check
      if(zneeded==2 || zneeded==3)
        {
         z.Resize(2,2);
         z[1].Set(1,1);
        }
      //--- return result
      return(result);
     }
//--- initialization
   maxit=30;
//--- allocation
   ArrayResizeAL(wtemp,n+1);
   ArrayResizeAL(work1,n);
   ArrayResizeAL(work2,n);
   ArrayResizeAL(workc,n+1);
   ArrayResizeAL(works,n+1);
//--- Determine the unit roundoff and over/underflow thresholds.
   eps=CMath::m_machineepsilon;
   eps2=CMath::Sqr(eps);
   safmin=CMath::m_minrealnumber;
   safmax=CMath::m_maxrealnumber;
   ssfmax=MathSqrt(safmax)/3;
   ssfmin=MathSqrt(safmin)/eps2;
//--- Here we are using transposition to get rid of column operations
   wastranspose=false;
   zrows=0;
//--- check
   if(zneeded==1)
      zrows=n;
//--- check
   if(zneeded==2)
      zrows=n;
//--- check
   if(zneeded==3)
      zrows=1;
//--- check
   if(zneeded==1)
     {
      wastranspose=true;
      //--- function call
      CBlas::InplaceTranspose(z,1,n,1,n,wtemp);
     }
//--- check
   if(zneeded==2)
     {
      wastranspose=true;
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
//--- check
   if(zneeded==3)
     {
      wastranspose=false;
      z.Resize(2,n+1);
      for(j=1;j<=n;j++)
        {
         //--- check
         if(j==1)
            z[1].Set(j,1);
         else
            z[1].Set(j,0);
        }
     }
//--- initialization
   nmaxit=n*maxit;
   jtot=0;
//--- Determine where the matrix splits and choose QL or QR iteration
//--- for each block, according to whether top or bottom diagonal
//--- element is smaller.
   l1=1;
   nm1=n-1;
   while(true)
     {
      //--- check
      if(l1>n)
         break;
      //--- check
      if(l1>1)
         e[l1-1]=0;
      gotoflag=false;
      m=l1;
      //--- check
      if(l1<=nm1)
        {
         for(m=l1;m<=nm1;m++)
           {
            tst=MathAbs(e[m]);
            //--- check
            if(tst==0.0)
              {
               gotoflag=true;
               //--- break the cycle
               break;
              }
            //--- check
            if(tst<=MathSqrt(MathAbs(d[m]))*MathSqrt(MathAbs(d[m+1]))*eps)
              {
               e[m]=0;
               gotoflag=true;
               //--- break the cycle
               break;
              }
           }
        }
      //--- check
      if(!gotoflag)
         m=n;
      //--- change values
      l=l1;
      lsv=l;
      lend=m;
      lendsv=lend;
      l1=m+1;
      //--- check
      if(lend==l)
         continue;
      //--- Scale submatrix in rows and columns L to LEND
      if(l==lend)
         anorm=MathAbs(d[l]);
      else
        {
         anorm=MathMax(MathAbs(d[l])+MathAbs(e[l]),MathAbs(e[lend-1])+MathAbs(d[lend]));
         for(i=l+1;i<=lend-1;i++)
            anorm=MathMax(anorm,MathAbs(d[i])+MathAbs(e[i])+MathAbs(e[i-1]));
        }
      iscale=0;
      //--- check
      if(anorm==0.0)
         continue;
      //--- check
      if(anorm>(double)(ssfmax))
        {
         iscale=1;
         tmp=ssfmax/anorm;
         tmpint=lend-1;
         for(i_=l;i_<=lend;i_++)
            d[i_]=tmp*d[i_];
         for(i_=l;i_<=tmpint;i_++)
            e[i_]=tmp*e[i_];
        }
      //--- check
      if(anorm<ssfmin)
        {
         iscale=2;
         tmp=ssfmin/anorm;
         tmpint=lend-1;
         for(i_=l;i_<=lend;i_++)
            d[i_]=tmp*d[i_];
         for(i_=l;i_<=tmpint;i_++)
            e[i_]=tmp*e[i_];
        }
      //--- Choose between QL and QR iteration
      if(MathAbs(d[lend])<MathAbs(d[l]))
        {
         lend=lsv;
         l=lendsv;
        }
      //--- check
      if(lend>l)
        {
         //--- QL Iteration
         //--- Look for small subdiagonal element.
         while(true)
           {
            gotoflag=false;
            //--- check
            if(l!=lend)
              {
               lendm1=lend-1;
               for(m=l;m<=lendm1;m++)
                 {
                  tst=CMath::Sqr(MathAbs(e[m]));
                  //--- check
                  if(tst<=eps2*MathAbs(d[m])*MathAbs(d[m+1])+safmin)
                    {
                     gotoflag=true;
                     //--- break the cycle
                     break;
                    }
                 }
              }
            //--- check
            if(!gotoflag)
               m=lend;
            //--- check
            if(m<lend)
               e[m]=0;
            p=d[l];
            //--- check
            if(m!=l)
              {
               //--- If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
               //--- to compute its eigensystem.
               if(m==l+1)
                 {
                  //--- check
                  if(zneeded>0)
                    {
                     //--- function call
                     TdEVDEv2(d[l],e[l],d[l+1],rt1,rt2,c,s);
                     //--- change values
                     work1[l]=c;
                     work2[l]=s;
                     workc[1]=work1[l];
                     works[1]=work2[l];
                     //--- check
                     if(!wastranspose)
                        CRotations::ApplyRotationsFromTheRight(false,1,zrows,l,l+1,workc,works,z,wtemp);
                     else
                        CRotations::ApplyRotationsFromTheLeft(false,l,l+1,1,zrows,workc,works,z,wtemp);
                    }
                  else
                  //--- function call
                     TdEVDE2(d[l],e[l],d[l+1],rt1,rt2);
                  //--- change values
                  d[l]=rt1;
                  d[l+1]=rt2;
                  e[l]=0;
                  l=l+2;
                  //--- check
                  if(l<=lend)
                     continue;
                  //--- break the cycle
                  break;
                 }
               //--- check
               if(jtot==nmaxit)
                  break;
               jtot=jtot+1;
               //--- Form shift.
               g=(d[l+1]-p)/(2*e[l]);
               //--- function call
               r=TdEVDPythag(g,1);
               g=d[m]-p+e[l]/(g+TdEVDExtSign(r,g));
               s=1;
               c=1;
               p=0;
               //--- Inner loop
               mm1=m-1;
               for(i=mm1;i>=l;i--)
                 {
                  f=s*e[i];
                  b=c*e[i];
                  //--- function call
                  CRotations::GenerateRotation(g,f,c,s,r);
                  //--- check
                  if(i!=m-1)
                     e[i+1]=r;
                  g=d[i+1]-p;
                  r=(d[i]-g)*s+2*c*b;
                  p=s*r;
                  d[i+1]=g+p;
                  g=c*r-b;
                  //--- If eigenvectors are desired, then save rotations.
                  if(zneeded>0)
                    {
                     work1[i]=c;
                     work2[i]=-s;
                    }
                 }
               //--- If eigenvectors are desired, then apply saved rotations.
               if(zneeded>0)
                 {
                  for(i=l;i<m;i++)
                    {
                     workc[i-l+1]=work1[i];
                     works[i-l+1]=work2[i];
                    }
                  //--- check
                  if(!wastranspose)
                    {
                     //--- function call
                     CRotations::ApplyRotationsFromTheRight(false,1,zrows,l,m,workc,works,z,wtemp);
                    }
                  else
                    {
                     //--- function call
                     CRotations::ApplyRotationsFromTheLeft(false,l,m,1,zrows,workc,works,z,wtemp);
                    }
                 }
               d[l]=d[l]-p;
               e[l]=g;
               continue;
              }
            //--- Eigenvalue found.
            d[l]=p;
            l=l+1;
            //--- check
            if(l<=lend)
               continue;
            //--- break the cycle
            break;
           }
        }
      else
        {
         //--- QR Iteration
         //--- Look for small superdiagonal element.
         while(true)
           {
            gotoflag=false;
            //--- check
            if(l!=lend)
              {
               lendp1=lend+1;
               for(m=l;m>=lendp1;m--)
                 {
                  tst=CMath::Sqr(MathAbs(e[m-1]));
                  //--- check
                  if(tst<=(double)(eps2*MathAbs(d[m])*MathAbs(d[m-1])+safmin))
                    {
                     gotoflag=true;
                     //--- break the cycle
                     break;
                    }
                 }
              }
            //--- check
            if(!gotoflag)
               m=lend;
            //--- check
            if(m>lend)
               e[m-1]=0;
            p=d[l];
            //--- check
            if(m!=l)
              {
               //--- If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
               //--- to compute its eigensystem.
               if(m==l-1)
                 {
                  //--- check
                  if(zneeded>0)
                    {
                     //--- function call
                     TdEVDEv2(d[l-1],e[l-1],d[l],rt1,rt2,c,s);
                     work1[m]=c;
                     work2[m]=s;
                     workc[1]=c;
                     works[1]=s;
                     //--- check
                     if(!wastranspose)
                       {
                        //--- function call
                        CRotations::ApplyRotationsFromTheRight(true,1,zrows,l-1,l,workc,works,z,wtemp);
                       }
                     else
                       {
                        //--- function call
                        CRotations::ApplyRotationsFromTheLeft(true,l-1,l,1,zrows,workc,works,z,wtemp);
                       }
                    }
                  else
                    {
                     //--- function call
                     TdEVDE2(d[l-1],e[l-1],d[l],rt1,rt2);
                    }
                  d[l-1]=rt1;
                  d[l]=rt2;
                  e[l-1]=0;
                  l=l-2;
                  //--- check
                  if(l>=lend)
                     continue;
                  //--- break the cycle
                  break;
                 }
               //--- check
               if(jtot==nmaxit)
                  break;
               jtot=jtot+1;
               //--- Form shift.
               g=(d[l-1]-p)/(2*e[l-1]);
               //--- function call
               r=TdEVDPythag(g,1);
               g=d[m]-p+e[l-1]/(g+TdEVDExtSign(r,g));
               s=1;
               c=1;
               p=0;
               //--- Inner loop
               lm1=l-1;
               for(i=m;i<=lm1;i++)
                 {
                  f=s*e[i];
                  b=c*e[i];
                  //--- function call
                  CRotations::GenerateRotation(g,f,c,s,r);
                  //--- check
                  if(i!=m)
                     e[i-1]=r;
                  //--- change values
                  g=d[i]-p;
                  r=(d[i+1]-g)*s+2*c*b;
                  p=s*r;
                  d[i]=g+p;
                  g=c*r-b;
                  //--- If eigenvectors are desired, then save rotations.
                  if(zneeded>0)
                    {
                     work1[i]=c;
                     work2[i]=s;
                    }
                 }
               //--- If eigenvectors are desired, then apply saved rotations.
               if(zneeded>0)
                 {
                  mm=l-m+1;
                  for(i=m;i<=l-1;i++)
                    {
                     workc[i-m+1]=work1[i];
                     works[i-m+1]=work2[i];
                    }
                  //--- check
                  if(!wastranspose)
                    {
                     //--- function call
                     CRotations::ApplyRotationsFromTheRight(true,1,zrows,m,l,workc,works,z,wtemp);
                    }
                  else
                    {
                     //--- function call
                     CRotations::ApplyRotationsFromTheLeft(true,m,l,1,zrows,workc,works,z,wtemp);
                    }
                 }
               d[l]=d[l]-p;
               e[lm1]=g;
               continue;
              }
            //--- Eigenvalue found.
            d[l]=p;
            l=l-1;
            //--- check
            if(l>=lend)
               continue;
            //--- break the cycle
            break;
           }
        }
      //--- Undo scaling if necessary
      if(iscale==1)
        {
         tmp=anorm/ssfmax;
         tmpint=lendsv-1;
         for(i_=lsv;i_<=lendsv;i_++)
            d[i_]=tmp*d[i_];
         for(i_=lsv;i_<=tmpint;i_++)
            e[i_]=tmp*e[i_];
        }
      //--- check
      if(iscale==2)
        {
         tmp=anorm/ssfmin;
         tmpint=lendsv-1;
         for(i_=lsv;i_<=lendsv;i_++)
            d[i_]=tmp*d[i_];
         for(i_=lsv;i_<=tmpint;i_++)
            e[i_]=tmp*e[i_];
        }
      //--- Check for no convergence to an eigenvalue after a total
      //--- of N*MAXIT iterations.
      if(jtot>=nmaxit)
        {
         result=false;
         //--- check
         if(wastranspose)
           {
            //--- function call
            CBlas::InplaceTranspose(z,1,n,1,n,wtemp);
           }
         //--- return result
         return(result);
        }
     }
//--- Order eigenvalues and eigenvectors.
   if(zneeded==0)
     {
      //--- Sort
      if(n==1)
         return(result);
      //--- check
      if(n==2)
        {
         //--- check
         if(d[1]>d[2])
           {
            tmp=d[1];
            d[1]=d[2];
            d[2]=tmp;
           }
         //--- return result
         return(result);
        }
      i=2;
      do
        {
         t=i;
         while(t!=1)
           {
            k=t/2;
            //--- check
            if(d[k]>=d[t])
               t=1;
            else
              {
               //--- change values
               tmp=d[k];
               d[k]=d[t];
               d[t]=tmp;
               t=k;
              }
           }
         i=i+1;
        }
      //--- cycle
      while(i<=n);
      i=n-1;
      do
        {
         tmp=d[i+1];
         d[i+1]=d[1];
         d[1]=tmp;
         t=1;
         while(t!=0)
           {
            k=2*t;
            //--- check
            if(k>i)
               t=0;
            else
              {
               //--- check
               if(k<i)
                 {
                  //--- check
                  if(d[k+1]>d[k])
                     k=k+1;
                 }
               //--- check
               if(d[t]>=d[k])
                  t=0;
               else
                 {
                  //--- change values
                  tmp=d[k];
                  d[k]=d[t];
                  d[t]=tmp;
                  t=k;
                 }
              }
           }
         i=i-1;
        }
      while(i>=1);
     }
   else
     {
      //--- Use Selection Sort to minimize swaps of eigenvectors
      for(ii=2;ii<=n;ii++)
        {
         i=ii-1;
         k=i;
         p=d[i];
         for(j=ii;j<=n;j++)
           {
            //--- check
            if(d[j]<p)
              {
               k=j;
               p=d[j];
              }
           }
         //--- check
         if(k!=i)
           {
            d[k]=d[i];
            d[i]=p;
            //--- check
            if(wastranspose)
              {
               for(i_=1;i_<=n;i_++)
                  wtemp[i_]=z[i][i_];
               for(i_=1;i_<=n;i_++)
                  z[i].Set(i_,z[k][i_]);
               for(i_=1;i_<=n;i_++)
                  z[k].Set(i_,wtemp[i_]);
              }
            else
              {
               for(i_=1;i_<=zrows;i_++)
                  wtemp[i_]=z[i_][i];
               for(i_=1;i_<=zrows;i_++)
                  z[i_].Set(i,z[i_][k]);
               for(i_=1;i_<=zrows;i_++)
                  z[i_].Set(k,wtemp[i_]);
              }
           }
        }
      //--- check
      if(wastranspose)
        {
         //--- function call
         CBlas::InplaceTranspose(z,1,n,1,n,wtemp);
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix     |
//|    [  A   B  ]                                                   |
//|    [  B   C  ].                                                  |
//| On return, RT1 is the eigenvalue of larger absolute value, and   |
//| RT2 is the eigenvalue of smaller absolute value.                 |
//|   -- LAPACK auxiliary routine (version 3.0) --                   |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      October 31, 1992                                            |
//+------------------------------------------------------------------+
static void CEigenVDetect::TdEVDE2(const double a,const double b,const double c,
                                   double &rt1,double &rt2)
  {
//--- create variables
   double ab=0;
   double acmn=0;
   double acmx=0;
   double adf=0;
   double df=0;
   double rt=0;
   double sm=0;
   double tb=0;
//--- initialization
   rt1=0;
   rt2=0;
   sm=a+c;
   df=a-c;
   adf=MathAbs(df);
   tb=b+b;
   ab=MathAbs(tb);
//--- check
   if(MathAbs(a)>MathAbs(c))
     {
      acmx=a;
      acmn=c;
     }
   else
     {
      acmx=c;
      acmn=a;
     }
//--- check
   if(adf>ab)
     {
      rt=adf*MathSqrt(1+CMath::Sqr(ab/adf));
     }
   else
     {
      //--- check
      if(adf<ab)
         rt=ab*MathSqrt(1+CMath::Sqr(adf/ab));
      else
        {
         //--- Includes case AB=ADF=0
         rt=ab*MathSqrt(2);
        }
     }
//--- check
   if(sm<0.0)
     {
      rt1=0.5*(sm-rt);
      //--- Order of execution important.
      //--- To get fully accurate smaller eigenvalue,
      //--- next line needs to be executed in higher precision.
      rt2=acmx/rt1*acmn-b/rt1*b;
     }
   else
     {
      //--- check
      if(sm>0.0)
        {
         rt1=0.5*(sm+rt);
         //--- Order of execution important.
         //--- To get fully accurate smaller eigenvalue,
         //--- next line needs to be executed in higher precision.
         rt2=acmx/rt1*acmn-b/rt1*b;
        }
      else
        {
         //--- Includes case RT1 = RT2 = 0
         rt1=0.5*rt;
         rt2=-(0.5*rt);
        }
     }
  }
//+------------------------------------------------------------------+
//| DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric     |
//| matrix                                                           |
//|    [  A   B  ]                                                   |
//|    [  B   C  ].                                                  |
//| On return, RT1 is the eigenvalue of larger absolute value, RT2 is|
//| the eigenvalue of smaller absolute value, and (CS1,SN1) is the   |
//| unit right eigenvector for RT1, giving the decomposition         |
//|    [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]         |
//|    [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].        |
//|   -- LAPACK auxiliary routine (version 3.0) --                   |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      October 31, 1992                                            |
//+------------------------------------------------------------------+
static void CEigenVDetect::TdEVDEv2(const double a,const double b,const double c,
                                    double &rt1,double &rt2,double &cs1,double &sn1)
  {
//--- create variables
   int    sgn1=0;
   int    sgn2=0;
   double ab=0;
   double acmn=0;
   double acmx=0;
   double acs=0;
   double adf=0;
   double cs=0;
   double ct=0;
   double df=0;
   double rt=0;
   double sm=0;
   double tb=0;
   double tn=0;
//--- initialization
   rt1=0;
   rt2=0;
   cs1=0;
   sn1=0;
//--- Compute the eigenvalues
   sm=a+c;
   df=a-c;
   adf=MathAbs(df);
   tb=b+b;
   ab=MathAbs(tb);
//--- check
   if(MathAbs(a)>MathAbs(c))
     {
      acmx=a;
      acmn=c;
     }
   else
     {
      acmx=c;
      acmn=a;
     }
//--- check
   if(adf>ab)
      rt=adf*MathSqrt(1+CMath::Sqr(ab/adf));
   else
     {
      //--- check
      if(adf<ab)
         rt=ab*MathSqrt(1+CMath::Sqr(adf/ab));
      else
        {
         //--- Includes case AB=ADF=0
         rt=ab*MathSqrt(2);
        }
     }
//--- check
   if(sm<0.0)
     {
      rt1=0.5*(sm-rt);
      sgn1=-1;
      //--- Order of execution important.
      //--- To get fully accurate smaller eigenvalue,
      //--- next line needs to be executed in higher precision.
      rt2=acmx/rt1*acmn-b/rt1*b;
     }
   else
     {
      //--- check
      if(sm>0.0)
        {
         rt1=0.5*(sm+rt);
         sgn1=1;
         //--- Order of execution important.
         //--- To get fully accurate smaller eigenvalue,
         //--- next line needs to be executed in higher precision.
         rt2=acmx/rt1*acmn-b/rt1*b;
        }
      else
        {
         //--- Includes case RT1 = RT2 = 0
         rt1=0.5*rt;
         rt2=-(0.5*rt);
         sgn1=1;
        }
     }
//--- Compute the eigenvector
   if(df>=0.0)
     {
      cs=df+rt;
      sgn2=1;
     }
   else
     {
      cs=df-rt;
      sgn2=-1;
     }
   acs=MathAbs(cs);
//--- check
   if(acs>ab)
     {
      ct=-(tb/cs);
      sn1=1/MathSqrt(1+ct*ct);
      cs1=ct*sn1;
     }
   else
     {
      //--- check
      if(ab==0.0)
        {
         cs1=1;
         sn1=0;
        }
      else
        {
         tn=-(cs/tb);
         cs1=1/MathSqrt(1+tn*tn);
         sn1=tn*cs1;
        }
     }
//--- check
   if(sgn1==sgn2)
     {
      tn=cs1;
      cs1=-sn1;
      sn1=tn;
     }
  }
//+------------------------------------------------------------------+
//| Internal routine                                                 |
//+------------------------------------------------------------------+
static double CEigenVDetect::TdEVDPythag(const double a,const double b)
  {
//--- create variables
   double result=0;
//--- check
   if(MathAbs(a)<MathAbs(b))
      result=MathAbs(b)*MathSqrt(1+CMath::Sqr(a/b));
   else
      result=MathAbs(a)*MathSqrt(1+CMath::Sqr(b/a));
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Internal routine                                                 |
//+------------------------------------------------------------------+
static double CEigenVDetect::TdEVDExtSign(const double a,const double b)
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
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static bool CEigenVDetect::InternalBisectionEigenValues(double &cd[],double &ce[],
                                                        const int n,int irange,
                                                        const int iorder,const double vl,
                                                        const double vu,const int il,
                                                        const int iu,const double abstol,
                                                        double &w[],int &m,
                                                        int &nsplit,int &iblock[],
                                                        int &isplit[],int &errorcode)
  {
//--- create variables
   bool   result;
   double fudge=0;
   double relfac=0;
   bool   ncnvrg;
   bool   toofew;
   int    ib=0;
   int    ibegin=0;
   int    idiscl=0;
   int    idiscu=0;
   int    ie=0;
   int    iend=0;
   int    iinfo=0;
   int    im=0;
   int    iin=0;
   int    ioff=0;
   int    iout=0;
   int    itmax=0;
   int    iw=0;
   int    iwoff=0;
   int    j=0;
   int    itmp1=0;
   int    jb=0;
   int    jdisc=0;
   int    je=0;
   int    nwl=0;
   int    nwu=0;
   int    tmpi=0;
   double atoli=0;
   double bnorm=0;
   double gl=0;
   double gu=0;
   double pivmin=0;
   double rtoli=0;
   double safemn=0;
   double tmp1=0;
   double tmp2=0;
   double tnorm=0;
   double ulp=0;
   double wkill=0;
   double wl=0;
   double wlu=0;
   double wu=0;
   double wul=0;
   double scalefactor=0;
   double t=0;
//--- create arrays
   int    idumma[];
   double work[];
   int    iwork[];
   int    ia1s2[];
   double ra1s2[];
   double ra1siin[];
   double ra2siin[];
   double ra3siin[];
   double ra4siin[];
   int    iworkspace[];
   double rworkspace[];
//--- create matrix
   CMatrixDouble ra1s2x2;
   CMatrixInt    ia1s2x2;
   CMatrixDouble ra1siinx2;
   CMatrixInt    ia1siinx2;
//--- create copy
   double d[];
   ArrayResizeAL(d,ArraySize(cd));
   ArrayCopy(d,cd);
//--- create copy
   double e[];
   ArrayResizeAL(e,ArraySize(ce));
   ArrayCopy(e,ce);
//--- initialization
   m=0;
   nsplit=0;
   errorcode=0;
//--- check
   if(n==0)
      return(true);
//--- Get machine constants
//--- NB is the minimum vector length for vector bisection, or 0
//--- if only scalar is to be done.
   fudge=2;
   relfac=2;
   safemn=CMath::m_minrealnumber;
   ulp=2*CMath::m_machineepsilon;
   rtoli=ulp*relfac;
//--- allocation
   ArrayResizeAL(idumma,2);
   ArrayResizeAL(work,4*n+1);
   ArrayResizeAL(iwork,3*n+1);
   ArrayResizeAL(w,n+1);
   ArrayResizeAL(iblock,n+1);
   ArrayResizeAL(isplit,n+1);
   ArrayResizeAL(ia1s2,3);
   ArrayResizeAL(ra1s2,3);
   ArrayResizeAL(ra1siin,n+1);
   ArrayResizeAL(ra2siin,n+1);
   ArrayResizeAL(ra3siin,n+1);
   ArrayResizeAL(ra4siin,n+1);
   ArrayResizeAL(iworkspace,n+1);
   ArrayResizeAL(rworkspace,n+1);
   ra1siinx2.Resize(n+1,3);
   ia1siinx2.Resize(n+1,3);
   ra1s2x2.Resize(3,3);
   ia1s2x2.Resize(3,3);
//--- initialization
   wlu=0;
   wul=0;
//--- Check for Errors
   result=false;
   errorcode=0;
//--- check
   if(irange<=0 || irange>=4)
      errorcode=-4;
//--- check
   if(iorder<=0 || iorder>=3)
      errorcode=-5;
//--- check
   if(n<0)
      errorcode=-3;
//--- check
   if(irange==2 && vl>=vu)
      errorcode=-6;
//--- check
   if(irange==3 && (il<1 || il>MathMax(1,n)))
      errorcode=-8;
//--- check
   if(irange==3 && (iu<MathMin(n,il) || iu>n))
      errorcode=-9;
//--- check
   if(errorcode!=0)
      return(result);
//--- Initialize error flags
   ncnvrg=false;
   toofew=false;
//--- Simplifications:
   if(irange==3 && il==1 && iu==n)
      irange=1;
//--- Special Case when N=1
   if(n==1)
     {
      nsplit=1;
      isplit[1]=1;
      //--- check
      if((irange==2 && vl>=d[1]) || vu<d[1])
         m=0;
      else
        {
         w[1]=d[1];
         iblock[1]=1;
         m=1;
        }
      //--- return result
      return(true);
     }
//--- Scaling
   t=MathAbs(d[n]);
   for(j=1;j<n;j++)
     {
      t=MathMax(t,MathAbs(d[j]));
      t=MathMax(t,MathAbs(e[j]));
     }
   scalefactor=1;
//--- check
   if(t!=0.0)
     {
      //--- check
      if(t>MathSqrt(MathSqrt(CMath::m_minrealnumber))*MathSqrt(CMath::m_maxrealnumber))
         scalefactor=t;
      //--- check
      if(t<MathSqrt(MathSqrt(CMath::m_maxrealnumber))*MathSqrt(CMath::m_minrealnumber))
         scalefactor=t;
      for(j=1;j<n;j++)
        {
         d[j]=d[j]/scalefactor;
         e[j]=e[j]/scalefactor;
        }
      d[n]=d[n]/scalefactor;
     }
//--- Compute Splitting Points
   nsplit=1;
   work[n]=0;
   pivmin=1;
   for(j=2;j<=n;j++)
     {
      tmp1=CMath::Sqr(e[j-1]);
      //--- check
      if(MathAbs(d[j]*d[j-1])*CMath::Sqr(ulp)+safemn>tmp1)
        {
         isplit[nsplit]=j-1;
         nsplit=nsplit+1;
         work[j-1]=0;
        }
      else
        {
         work[j-1]=tmp1;
         pivmin=MathMax(pivmin,tmp1);
        }
     }
   isplit[nsplit]=n;
   pivmin=pivmin*safemn;
//--- Compute Interval and ATOLI
   if(irange==3)
     {
      //--- RANGE='I': Compute the interval containing eigenvalues
      //---     IL through IU.
      //--- Compute Gershgorin interval for entire (split) matrix
      //--- and use it as the initial interval
      gu=d[1];
      gl=d[1];
      tmp1=0;
      for(j=1;j<n;j++)
        {
         //--- change values
         tmp2=MathSqrt(work[j]);
         gu=MathMax(gu,d[j]+tmp1+tmp2);
         gl=MathMin(gl,d[j]-tmp1-tmp2);
         tmp1=tmp2;
        }
      //--- change values
      gu=MathMax(gu,d[n]+tmp1);
      gl=MathMin(gl,d[n]-tmp1);
      tnorm=MathMax(MathAbs(gl),MathAbs(gu));
      gl=gl-fudge*tnorm*ulp*n-fudge*2*pivmin;
      gu=gu+fudge*tnorm*ulp*n+fudge*pivmin;
      //--- Compute Iteration parameters
      itmax=(int)MathCeil((MathLog(tnorm+pivmin)-MathLog(pivmin))/MathLog(2))+2;
      //--- check
      if(abstol<=0.0)
         atoli=ulp*tnorm;
      else
         atoli=abstol;
      //--- change values
      work[n+1]=gl;
      work[n+2]=gl;
      work[n+3]=gu;
      work[n+4]=gu;
      work[n+5]=gl;
      work[n+6]=gu;
      iwork[1]=-1;
      iwork[2]=-1;
      iwork[3]=n+1;
      iwork[4]=n+1;
      iwork[5]=il-1;
      iwork[6]=iu;
      //--- Calling DLAEBZ
      //--- DLAEBZ( 3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E,
      //---    WORK, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT,
      //---    IWORK, W, IBLOCK, IINFO )
      ia1s2[1]=iwork[5];
      ia1s2[2]=iwork[6];
      ra1s2[1]=work[n+5];
      ra1s2[2]=work[n+6];
      ra1s2x2[1].Set(1,work[n+1]);
      ra1s2x2[2].Set(1,work[n+2]);
      ra1s2x2[1].Set(2,work[n+3]);
      ra1s2x2[2].Set(2,work[n+4]);
      ia1s2x2[1].Set(1,iwork[1]);
      ia1s2x2[2].Set(1,iwork[2]);
      ia1s2x2[1].Set(2,iwork[3]);
      ia1s2x2[2].Set(2,iwork[4]);
      //--- function call
      InternalDLAEBZ(3,itmax,n,2,2,atoli,rtoli,pivmin,d,e,work,ia1s2,ra1s2x2,ra1s2,iout,ia1s2x2,w,iblock,iinfo);
      iwork[5]=ia1s2[1];
      iwork[6]=ia1s2[2];
      work[n+5]=ra1s2[1];
      work[n+6]=ra1s2[2];
      work[n+1]=ra1s2x2[1][1];
      work[n+2]=ra1s2x2[2][1];
      work[n+3]=ra1s2x2[1][2];
      work[n+4]=ra1s2x2[2][2];
      iwork[1]=ia1s2x2[1][1];
      iwork[2]=ia1s2x2[2][1];
      iwork[3]=ia1s2x2[1][2];
      iwork[4]=ia1s2x2[2][2];
      //--- check
      if(iwork[6]==iu)
        {
         //--- change values
         wl=work[n+1];
         wlu=work[n+3];
         nwl=iwork[1];
         wu=work[n+4];
         wul=work[n+2];
         nwu=iwork[4];
        }
      else
        {
         //--- change values
         wl=work[n+2];
         wlu=work[n+4];
         nwl=iwork[2];
         wu=work[n+3];
         wul=work[n+1];
         nwu=iwork[3];
        }
      //--- check
      if(nwl<0 || nwl>=n || nwu<1 || nwu>n)
        {
         errorcode=4;
         return(false);
        }
     }
   else
     {
      //--- RANGE='A' or 'V' -- Set ATOLI
      tnorm=MathMax(MathAbs(d[1])+MathAbs(e[1]),MathAbs(d[n])+MathAbs(e[n-1]));
      for(j=2;j<n;j++)
         tnorm=MathMax(tnorm,MathAbs(d[j])+MathAbs(e[j-1])+MathAbs(e[j]));
      //--- check
      if(abstol<=0.0)
         atoli=ulp*tnorm;
      else
         atoli=abstol;
      //--- check
      if(irange==2)
        {
         wl=vl;
         wu=vu;
        }
      else
        {
         wl=0;
         wu=0;
        }
     }
//--- Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
//--- NWL accumulates the number of eigenvalues .le. WL,
//--- NWU accumulates the number of eigenvalues .le. WU
   m=0;
   iend=0;
   errorcode=0;
   nwl=0;
   nwu=0;
   for(jb=1;jb<=nsplit;jb++)
     {
      ioff=iend;
      ibegin=ioff+1;
      iend=isplit[jb];
      iin=iend-ioff;
      //--- check
      if(iin==1)
        {
         //--- check
         if(irange==1 || wl>=d[ibegin]-pivmin)
            nwl=nwl+1;
         //--- check
         if(irange==1 || wu>=d[ibegin]-pivmin)
            nwu=nwu+1;
         //--- check
         if((irange==1 || wl<d[ibegin]-pivmin) && wu>=d[ibegin]-pivmin)
           {
            m=m+1;
            w[m]=d[ibegin];
            iblock[m]=jb;
           }
        }
      else
        {
         //--- General Case -- IIN > 1
         //--- Compute Gershgorin Interval
         //--- and use it as the initial interval
         gu=d[ibegin];
         gl=d[ibegin];
         tmp1=0;
         for(j=ibegin;j<=iend-1;j++)
           {
            //--- change values
            tmp2=MathAbs(e[j]);
            gu=MathMax(gu,d[j]+tmp1+tmp2);
            gl=MathMin(gl,d[j]-tmp1-tmp2);
            tmp1=tmp2;
           }
         //--- change values
         gu=MathMax(gu,d[iend]+tmp1);
         gl=MathMin(gl,d[iend]-tmp1);
         bnorm=MathMax(MathAbs(gl),MathAbs(gu));
         gl=gl-fudge*bnorm*ulp*iin-fudge*pivmin;
         gu=gu+fudge*bnorm*ulp*iin+fudge*pivmin;
         //--- Compute ATOLI for the current submatrix
         if(abstol<=0.0)
            atoli=ulp*MathMax(MathAbs(gl),MathAbs(gu));
         else
            atoli=abstol;
         //--- check
         if(irange>1)
           {
            //--- check
            if(gu<wl)
              {
               nwl=nwl+iin;
               nwu=nwu+iin;
               continue;
              }
            gl=MathMax(gl,wl);
            gu=MathMin(gu,wu);
            //--- check
            if(gl>=gu)
               continue;
           }
         //--- Set Up Initial Interval
         work[n+1]=gl;
         work[n+iin+1]=gu;
         //--- Calling DLAEBZ
         //--- CALL DLAEBZ( 1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
         //---    D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),
         //---    IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM,
         //---    IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
         for(tmpi=1;tmpi<=iin;tmpi++)
           {
            ra1siin[tmpi]=d[ibegin-1+tmpi];
            //--- check
            if(ibegin-1+tmpi<n)
               ra2siin[tmpi]=e[ibegin-1+tmpi];
            //--- change values
            ra3siin[tmpi]=work[ibegin-1+tmpi];
            ra1siinx2[tmpi].Set(1,work[n+tmpi]);
            ra1siinx2[tmpi].Set(2,work[n+tmpi+iin]);
            ra4siin[tmpi]=work[n+2*iin+tmpi];
            rworkspace[tmpi]=w[m+tmpi];
            iworkspace[tmpi]=iblock[m+tmpi];
            ia1siinx2[tmpi].Set(1,iwork[tmpi]);
            ia1siinx2[tmpi].Set(2,iwork[tmpi+iin]);
           }
         //--- function call
         InternalDLAEBZ(1,0,iin,iin,1,atoli,rtoli,pivmin,ra1siin,ra2siin,ra3siin,idumma,ra1siinx2,ra4siin,im,ia1siinx2,rworkspace,iworkspace,iinfo);
         for(tmpi=1;tmpi<=iin;tmpi++)
           {
            //--- change values
            work[n+tmpi]=ra1siinx2[tmpi][1];
            work[n+tmpi+iin]=ra1siinx2[tmpi][2];
            work[n+2*iin+tmpi]=ra4siin[tmpi];
            w[m+tmpi]=rworkspace[tmpi];
            iblock[m+tmpi]=iworkspace[tmpi];
            iwork[tmpi]=ia1siinx2[tmpi][1];
            iwork[tmpi+iin]=ia1siinx2[tmpi][2];
           }
         nwl=nwl+iwork[1];
         nwu=nwu+iwork[iin+1];
         iwoff=m-iwork[1];
         //--- Compute Eigenvalues
         itmax=(int)MathCeil((MathLog(gu-gl+pivmin)-MathLog(pivmin))/MathLog(2))+2;
         //--- Calling DLAEBZ
         //--- CALL DLAEBZ( 2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
         //---    D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),
         //---    IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT,
         //---    IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
         for(tmpi=1;tmpi<=iin;tmpi++)
           {
            ra1siin[tmpi]=d[ibegin-1+tmpi];
            //--- check
            if(ibegin-1+tmpi<n)
               ra2siin[tmpi]=e[ibegin-1+tmpi];
            //--- change values
            ra3siin[tmpi]=work[ibegin-1+tmpi];
            ra1siinx2[tmpi].Set(1,work[n+tmpi]);
            ra1siinx2[tmpi].Set(2,work[n+tmpi+iin]);
            ra4siin[tmpi]=work[n+2*iin+tmpi];
            rworkspace[tmpi]=w[m+tmpi];
            iworkspace[tmpi]=iblock[m+tmpi];
            ia1siinx2[tmpi].Set(1,iwork[tmpi]);
            ia1siinx2[tmpi].Set(2,iwork[tmpi+iin]);
           }
         //--- function call
         InternalDLAEBZ(2,itmax,iin,iin,1,atoli,rtoli,pivmin,ra1siin,ra2siin,ra3siin,idumma,ra1siinx2,ra4siin,iout,ia1siinx2,rworkspace,iworkspace,iinfo);
         for(tmpi=1;tmpi<=iin;tmpi++)
           {
            //--- change values
            work[n+tmpi]=ra1siinx2[tmpi][1];
            work[n+tmpi+iin]=ra1siinx2[tmpi][2];
            work[n+2*iin+tmpi]=ra4siin[tmpi];
            w[m+tmpi]=rworkspace[tmpi];
            iblock[m+tmpi]=iworkspace[tmpi];
            iwork[tmpi]=ia1siinx2[tmpi][1];
            iwork[tmpi+iin]=ia1siinx2[tmpi][2];
           }
         //--- Copy Eigenvalues Into W and IBLOCK
         //--- Use -JB for block number for unconverged eigenvalues.
         for(j=1;j<=iout;j++)
           {
            tmp1=0.5*(work[j+n]+work[j+iin+n]);
            //--- Flag non-convergence.
            if(j>iout-iinfo)
              {
               ncnvrg=true;
               ib=-jb;
              }
            else
               ib=jb;
            for(je=iwork[j]+1+iwoff;je<=iwork[j+iin]+iwoff;je++)
              {
               w[je]=tmp1;
               iblock[je]=ib;
              }
           }
         m=m+im;
        }
     }
//--- If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
//--- If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
   if(irange==3)
     {
      im=0;
      idiscl=il-1-nwl;
      idiscu=nwu-iu;
      //--- check
      if(idiscl>0 || idiscu>0)
        {
         for(je=1;je<=m;je++)
           {
            //--- check
            if(w[je]<=wlu && idiscl>0)
              {
               idiscl=idiscl-1;
              }
            else
              {
               //--- check
               if(w[je]>=wul && idiscu>0)
                  idiscu=idiscu-1;
               else
                 {
                  im=im+1;
                  w[im]=w[je];
                  iblock[im]=iblock[je];
                 }
              }
           }
         m=im;
        }
      //--- check
      if(idiscl>0 || idiscu>0)
        {
         //--- Code to deal with effects of bad arithmetic:
         //--- Some low eigenvalues to be discarded are not in (WL,WLU],
         //--- or high eigenvalues to be discarded are not in (WUL,WU]
         //--- so just kill off the smallest IDISCL/largest IDISCU
         //--- eigenvalues, by simply finding the smallest/largest
         //--- eigenvalue(s).
         //--- (If N(w) is monotone non-decreasing, this should never
         //---  happen.)
         if(idiscl>0)
           {
            wkill=wu;
            for(jdisc=1;jdisc<=idiscl;jdisc++)
              {
               iw=0;
               for(je=1;je<=m;je++)
                 {
                  //--- check
                  if(iblock[je]!=0 && (w[je]<(double)(wkill) || iw==0))
                    {
                     iw=je;
                     wkill=w[je];
                    }
                 }
               iblock[iw]=0;
              }
           }
         //--- check
         if(idiscu>0)
           {
            wkill=wl;
            for(jdisc=1;jdisc<=idiscu;jdisc++)
              {
               iw=0;
               for(je=1;je<=m;je++)
                 {
                  //--- check
                  if(iblock[je]!=0 && (w[je]>(double)(wkill) || iw==0))
                    {
                     iw=je;
                     wkill=w[je];
                    }
                 }
               iblock[iw]=0;
              }
           }
         im=0;
         for(je=1;je<=m;je++)
           {
            //--- check
            if(iblock[je]!=0)
              {
               im=im+1;
               w[im]=w[je];
               iblock[im]=iblock[je];
              }
           }
         m=im;
        }
      //--- check
      if(idiscl<0 || idiscu<0)
         toofew=true;
     }
//--- If ORDER='B', do nothing -- the eigenvalues are already sorted
//---    by block.
//--- If ORDER='E', sort the eigenvalues from smallest to largest
   if(iorder==1 && nsplit>1)
     {
      for(je=1;je<=m-1;je++)
        {
         ie=0;
         tmp1=w[je];
         for(j=je+1;j<=m;j++)
           {
            //--- check
            if(w[j]<tmp1)
              {
               ie=j;
               tmp1=w[j];
              }
           }
         //--- check
         if(ie!=0)
           {
            //--- change values
            itmp1=iblock[ie];
            w[ie]=w[je];
            iblock[ie]=iblock[je];
            w[je]=tmp1;
            iblock[je]=itmp1;
           }
        }
     }
   for(j=1;j<=m;j++)
      w[j]=w[j]*scalefactor;
   errorcode=0;
//--- check
   if(ncnvrg)
      errorcode=errorcode+1;
//--- check
   if(toofew)
      errorcode=errorcode+2;
   result=errorcode==0;
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CEigenVDetect::InternalDStein(const int n,double &d[],double &ce[],
                                          const int m,double &cw[],int &iblock[],
                                          int &isplit[],CMatrixDouble &z,
                                          int &ifail[],int &info)
  {
//--- create variables
   int    maxits=0;
   int    extra=0;
   int    b1=0;
   int    blksiz=0;
   int    bn=0;
   int    gpind=0;
   int    i=0;
   int    iinfo=0;
   int    its=0;
   int    j=0;
   int    j1=0;
   int    jblk=0;
   int    jmax=0;
   int    nblk=0;
   int    nrmchk=0;
   double dtpcrt=0;
   double eps=0;
   double eps1=0;
   double nrm=0;
   double onenrm=0;
   double ortol=0;
   double pertol=0;
   double scl=0;
   double sep=0;
   double tol=0;
   double xj=0;
   double xjm=0;
   double ztr=0;
   bool   tmpcriterion;
   int    ti=0;
   int    i1=0;
   int    i2=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- create arrays
   double work1[];
   double work2[];
   double work3[];
   double work4[];
   double work5[];
   int    iwork[];
//--- create copy
   double e[];
   ArrayResizeAL(e,ArraySize(ce));
   ArrayCopy(e,ce);
//--- create copy
   double w[];
   ArrayResizeAL(w,ArraySize(cw));
   ArrayCopy(w,cw);
//--- initialization
   info=0;
   maxits=5;
   extra=2;
//--- allocation
   ArrayResizeAL(work1,(int)MathMax(n,1)+1);
   ArrayResizeAL(work2,(int)MathMax(n-1,1)+1);
   ArrayResizeAL(work3,(int)MathMax(n,1)+1);
   ArrayResizeAL(work4,(int)MathMax(n,1)+1);
   ArrayResizeAL(work5,(int)MathMax(n,1)+1);
   ArrayResizeAL(iwork,(int)MathMax(n,1)+1);
   ArrayResizeAL(ifail,(int)MathMax(m,1)+1);
   z.Resize((int)MathMax(n,1)+1,(int)MathMax(m,1)+1);
//--- initialization
   gpind=0;
   onenrm=0;
   ortol=0;
   dtpcrt=0;
   xjm=0;
//--- check input parameters
   info=0;
   for(i=1;i<=m;i++)
      ifail[i]=0;
//--- check
   if(n<0)
     {
      info=-1;
      return;
     }
//--- check
   if(m<0 || m>n)
     {
      info=-4;
      return;
     }
   for(j=2;j<=m;j++)
     {
      //--- check
      if(iblock[j]<iblock[j-1])
        {
         info=-6;
         //--- break the cycle
         break;
        }
      //--- check
      if(iblock[j]==iblock[j-1] && w[j]<w[j-1])
        {
         info=-5;
         //--- break the cycle
         break;
        }
     }
//--- check
   if(info!=0)
      return;
//--- check
   if(n==0 || m==0)
      return;
//--- check
   if(n==1)
     {
      z[1].Set(1,1);
      return;
     }
//--- Some preparations
   ti=n-1;
   for(i_=1;i_<=ti;i_++)
      work1[i_]=e[i_];
   ArrayResizeAL(e,n+1);
   for(i_=1;i_<=ti;i_++)
      e[i_]=work1[i_];
   for(i_=1;i_<=m;i_++)
      work1[i_]=w[i_];
   ArrayResizeAL(w,n+1);
   for(i_=1;i_<=m;i_++)
      w[i_]=work1[i_];
//--- Get machine constants.
   eps=CMath::m_machineepsilon;
//--- Compute eigenvectors of matrix blocks.
   j1=1;
   for(nblk=1;nblk<=iblock[m];nblk++)
     {
      //--- Find starting and ending indices of block nblk.
      if(nblk==1)
         b1=1;
      else
         b1=isplit[nblk-1]+1;
      bn=isplit[nblk];
      blksiz=bn-b1+1;
      //--- check
      if(blksiz!=1)
        {
         //--- Compute reorthogonalization criterion and stopping criterion.
         gpind=b1;
         onenrm=MathAbs(d[b1])+MathAbs(e[b1]);
         onenrm=MathMax(onenrm,MathAbs(d[bn])+MathAbs(e[bn-1]));
         for(i=b1+1;i<=bn-1;i++)
            onenrm=MathMax(onenrm,MathAbs(d[i])+MathAbs(e[i-1])+MathAbs(e[i]));
         ortol=0.001*onenrm;
         dtpcrt=MathSqrt(0.1/blksiz);
        }
      //--- Loop through eigenvalues of block nblk.
      jblk=0;
      for(j=j1;j<=m;j++)
        {
         //--- check
         if(iblock[j]!=nblk)
           {
            j1=j;
            //--- break the cycle
            break;
           }
         jblk=jblk+1;
         xj=w[j];
         //--- check
         if(blksiz==1)
           {
            //--- Skip all the work if the block size is one.
            work1[1]=1;
           }
         else
           {
            //--- If eigenvalues j and j-1 are too close, add a relatively
            //--- small perturbation.
            if(jblk>1)
              {
               eps1=MathAbs(eps*xj);
               pertol=10*eps1;
               sep=xj-xjm;
               //--- check
               if(sep<pertol)
                  xj=xjm+pertol;
              }
            its=0;
            nrmchk=0;
            //--- Get random starting vector.
            for(ti=1;ti<=blksiz;ti++)
               work1[ti]=2*CMath::RandomReal()-1;
            //--- Copy the matrix T so it won't be destroyed in factorization.
            for(ti=1;ti<=blksiz-1;ti++)
              {
               work2[ti]=e[b1+ti-1];
               work3[ti]=e[b1+ti-1];
               work4[ti]=d[b1+ti-1];
              }
            work4[blksiz]=d[b1+blksiz-1];
            //--- Compute LU factors with partial pivoting  ( PT = LU )
            tol=0;
            TdIninternalDLAGTF(blksiz,work4,xj,work2,work3,tol,work5,iwork,iinfo);
            //--- Update iteration count.
            do
              {
               its=its+1;
               //--- check
               if(its>maxits)
                 {
                  //--- If stopping criterion was not satisfied, update info and
                  //--- store eigenvector number in array ifail.
                  info=info+1;
                  ifail[info]=j;
                  break;
                 }
               //--- Normalize and scale the righthand side vector Pb.
               v=0;
               for(ti=1;ti<=blksiz;ti++)
                  v=v+MathAbs(work1[ti]);
               scl=blksiz*onenrm*MathMax(eps,MathAbs(work4[blksiz]))/v;
               for(i_=1;i_<=blksiz;i_++)
                  work1[i_]=scl*work1[i_];
               //--- Solve the system LU = Pb.
               TdIninternalDLAGTS(blksiz,work4,work2,work3,work5,iwork,work1,tol,iinfo);
               //--- Reorthogonalize by modified Gram-Schmidt if eigenvalues are
               //--- close enough.
               if(jblk!=1)
                 {
                  //--- check
                  if(MathAbs(xj-xjm)>ortol)
                     gpind=j;
                  //--- check
                  if(gpind!=j)
                    {
                     for(i=gpind;i<j;i++)
                       {
                        i1=b1;
                        i2=b1+blksiz-1;
                        i1_=i1-1;
                        ztr=0.0;
                        for(i_=1;i_<=blksiz;i_++)
                           ztr+=work1[i_]*z[i_+i1_][i];
                        i1_=i1-1;
                        for(i_=1;i_<=blksiz;i_++)
                           work1[i_]=work1[i_]-ztr*z[i_+i1_][i];
                       }
                    }
                 }
               //--- Check the infinity norm of the iterate.
               jmax=CBlas::VectorIdxAbsMax(work1,1,blksiz);
               nrm=MathAbs(work1[jmax]);
               //--- Continue for additional iterations after norm reaches
               //--- stopping criterion.
               tmpcriterion=false;
               //--- check
               if(nrm<dtpcrt)
                  tmpcriterion=true;
               else
                 {
                  nrmchk=nrmchk+1;
                  //--- check
                  if(nrmchk<extra+1)
                     tmpcriterion=true;
                 }
              }
            while(tmpcriterion);
            //--- Accept iterate as jth eigenvector.
            scl=1/CBlas::VectorNorm2(work1,1,blksiz);
            jmax=CBlas::VectorIdxAbsMax(work1,1,blksiz);
            //--- check
            if(work1[jmax]<0.0)
               scl=-scl;
            for(i_=1;i_<=blksiz;i_++)
               work1[i_]=scl*work1[i_];
           }
         for(i=1;i<=n;i++)
            z[i].Set(j,0);
         for(i=1;i<=blksiz;i++)
            z[b1+i-1].Set(j,work1[i]);
         //--- Save the shift to check eigenvalue spacing at next iteration.
         xjm=xj;
        }
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CEigenVDetect::TdIninternalDLAGTF(const int n,double &a[],const double lambdav,
                                              double &b[],double &c[],double tol,
                                              double &d[],int &iin[],int &info)
  {
//--- create variables
   int    k=0;
   double eps=0;
   double mult=0;
   double piv1=0;
   double piv2=0;
   double scale1=0;
   double scale2=0;
   double temp=0;
   double tl=0;
//--- initialization
   info=0;
//--- check
   if(n<0)
     {
      info=-1;
      return;
     }
//--- check
   if(n==0)
      return;
   a[1]=a[1]-lambdav;
   iin[n]=0;
//--- check
   if(n==1)
     {
      //--- check
      if(a[1]==0.0)
         iin[1]=1;
      //--- exit the function
      return;
     }
//--- initialization
   eps=CMath::m_machineepsilon;
   tl=MathMax(tol,eps);
   scale1=MathAbs(a[1])+MathAbs(b[1]);
   for(k=1;k<n;k++)
     {
      a[k+1]=a[k+1]-lambdav;
      scale2=MathAbs(c[k])+MathAbs(a[k+1]);
      //--- check
      if(k<n-1)
         scale2=scale2+MathAbs(b[k+1]);
      //--- check
      if(a[k]==0.0)
         piv1=0;
      else
         piv1=MathAbs(a[k])/scale1;
      //--- check
      if(c[k]==0.0)
        {
         iin[k]=0;
         piv2=0;
         scale1=scale2;
         //--- check
         if(k<n-1)
            d[k]=0;
        }
      else
        {
         piv2=MathAbs(c[k])/scale2;
         //--- check
         if(piv2<=piv1)
           {
            //--- change values
            iin[k]=0;
            scale1=scale2;
            c[k]=c[k]/a[k];
            a[k+1]=a[k+1]-c[k]*b[k];
            //--- check
            if(k<n-1)
               d[k]=0;
           }
         else
           {
            //--- change values
            iin[k]=1;
            mult=a[k]/c[k];
            a[k]=c[k];
            temp=a[k+1];
            a[k+1]=b[k]-mult*temp;
            if(k<n-1)
              {
               d[k]=b[k+1];
               b[k+1]=-(mult*d[k]);
              }
            b[k]=temp;
            c[k]=mult;
           }
        }
      //--- check
      if(MathMax(piv1,piv2)<=tl && iin[n]==0)
         iin[n]=k;
     }
//--- check
   if(MathAbs(a[n])<=scale1*tl && iin[n]==0)
      iin[n]=n;
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CEigenVDetect::TdIninternalDLAGTS(const int n,double &a[],double &b[],
                                              double &c[],double &d[],int &iin[],
                                              double &y[],double &tol,int &info)
  {
//--- create variables
   int    k=0;
   double absak=0;
   double ak=0;
   double bignum=0;
   double eps=0;
   double pert=0;
   double sfmin=0;
   double temp=0;
//--- initialization
   info=0;
//--- check
   if(n<0)
     {
      info=-1;
      return;
     }
//--- check
   if(n==0)
      return;
//--- initialization
   eps=CMath::m_machineepsilon;
   sfmin=CMath::m_minrealnumber;
   bignum=1/sfmin;
//--- check
   if(tol<=0.0)
     {
      tol=MathAbs(a[1]);
      //--- check
      if(n>1)
         tol=MathMax(tol,MathMax(MathAbs(a[2]),MathAbs(b[1])));
      for(k=3;k<=n;k++)
         tol=MathMax(tol,MathMax(MathAbs(a[k]),MathMax(MathAbs(b[k-1]),MathAbs(d[k-2]))));
      tol=tol*eps;
      //--- check
      if(tol==0.0)
         tol=eps;
     }
   for(k=2;k<=n;k++)
     {
      //--- check
      if(iin[k-1]==0)
         y[k]=y[k]-c[k-1]*y[k-1];
      else
        {
         temp=y[k-1];
         y[k-1]=y[k];
         y[k]=temp-c[k-1]*y[k];
        }
     }
   for(k=n;k>=1;k--)
     {
      //--- check
      if(k<=n-2)
         temp=y[k]-b[k]*y[k+1]-d[k]*y[k+2];
      else
        {
         //--- check
         if(k==n-1)
            temp=y[k]-b[k]*y[k+1];
         else
            temp=y[k];
        }
      ak=a[k];
      pert=MathAbs(tol);
      //--- check
      if(ak<0.0)
         pert=-pert;
      while(true)
        {
         absak=MathAbs(ak);
         //--- check
         if(absak<1.0)
           {
            //--- check
            if(absak<sfmin)
              {
               //--- check
               if(absak==0.0 || MathAbs(temp)*sfmin>absak)
                 {
                  ak=ak+pert;
                  pert=2*pert;
                  continue;
                 }
               else
                 {
                  temp=temp*bignum;
                  ak=ak*bignum;
                 }
              }
            else
              {
               //--- check
               if(MathAbs(temp)>absak*bignum)
                 {
                  ak=ak+pert;
                  pert=2*pert;
                  continue;
                 }
              }
           }
         //--- break the cycle
         break;
        }
      y[k]=temp/ak;
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CEigenVDetect::InternalDLAEBZ(const int ijob,const int nitmax,
                                          const int n,const int mmax,const int minp,
                                          const double abstol,const double reltol,
                                          const double pivmin,double &d[],
                                          double &e[],double &e2[],int &nval[],
                                          CMatrixDouble &ab,double &c[],int &mout,
                                          CMatrixInt &nab,double &work[],
                                          int &iwork[],int &info)
  {
//--- create variables
   int    itmp1=0;
   int    itmp2=0;
   int    j=0;
   int    ji=0;
   int    jit=0;
   int    jp=0;
   int    kf=0;
   int    kfnew=0;
   int    kl=0;
   int    klnew=0;
   double tmp1=0;
   double tmp2=0;
//--- initialization
   mout=0;
   info=0;
//--- check
   if(ijob<1 || ijob>3)
     {
      info=-1;
      //--- exit the function
      return;
     }
//--- Initialize NAB
   if(ijob==1)
     {
      //--- Compute the number of eigenvalues in the initial intervals.
      mout=0;
      //--- DIR$ NOVECTOR
      for(ji=1;ji<=minp;ji++)
        {
         for(jp=1;jp<=2;jp++)
           {
            tmp1=d[1]-ab[ji][jp];
            //--- check
            if(MathAbs(tmp1)<pivmin)
               tmp1=-pivmin;
            nab[ji].Set(jp,0);
            //--- check
            if(tmp1<=0.0)
               nab[ji].Set(jp,1);
            for(j=2;j<=n;j++)
              {
               tmp1=d[j]-e2[j-1]/tmp1-ab[ji][jp];
               //--- check
               if(MathAbs(tmp1)<pivmin)
                  tmp1=-pivmin;
               //--- check
               if(tmp1<=0.0)
                  nab[ji].Set(jp,nab[ji][jp]+1);
              }
           }
         mout=mout+nab[ji][2]-nab[ji][1];
        }
      //--- exit the function
      return;
     }
//--- Initialize for loop
//--- KF and KL have the following meaning:
//---   Intervals 1,...,KF-1 have converged.
//---   Intervals KF,...,KL  still need to be refined.
   kf=1;
   kl=minp;
//--- If IJOB=2, initialize C.
//--- If IJOB=3, use the user-supplied starting point.
   if(ijob==2)
     {
      for(ji=1;ji<=minp;ji++)
         c[ji]=0.5*(ab[ji][1]+ab[ji][2]);
     }
//--- Iteration loop
   for(jit=1;jit<=nitmax;jit++)
     {
      //--- Loop over intervals
      //--- Serial Version of the loop
      klnew=kl;
      for(ji=kf;ji<=kl;ji++)
        {
         //--- Compute N(w), the number of eigenvalues less than w
         tmp1=c[ji];
         tmp2=d[1]-tmp1;
         itmp1=0;
         //--- check
         if(tmp2<=pivmin)
           {
            itmp1=1;
            tmp2=MathMin(tmp2,-pivmin);
           }
         //--- A series of compiler directives to defeat vectorization
         //--- for the next loop
         //--- *$PL$ CMCHAR=' '
         //--- CDIR$          NEXTSCALAR
         //--- C$DIR          SCALAR
         //--- CDIR$          NEXT SCALAR
         //--- CVD$L          NOVECTOR
         //--- CDEC$          NOVECTOR
         //--- CVD$           NOVECTOR
         //--- *VDIR          NOVECTOR
         //--- *VOCL          LOOP,SCALAR
         //--- CIBM           PREFER SCALAR
         //--- *$PL$ CMCHAR='*'
         for(j=2;j<=n;j++)
           {
            tmp2=d[j]-e2[j-1]/tmp2-tmp1;
            //--- check
            if(tmp2<=pivmin)
              {
               itmp1=itmp1+1;
               tmp2=MathMin(tmp2,-pivmin);
              }
           }
         //--- check
         if(ijob<=2)
           {
            //--- IJOB=2: Choose all intervals containing eigenvalues.
            //--- Insure that N(w) is monotone
            itmp1=MathMin(nab[ji][2],MathMax(nab[ji][1],itmp1));
            //--- Update the Queue -- add intervals if both halves
            //--- contain eigenvalues.
            if(itmp1==nab[ji][2])
              {
               //--- No eigenvalue in the upper interval:
               //--- just use the lower interval.
               ab[ji].Set(2,tmp1);
              }
            else
              {
               //--- check
               if(itmp1==nab[ji][1])
                 {
                  //--- No eigenvalue in the lower interval:
                  //--- just use the upper interval.
                  ab[ji].Set(1,tmp1);
                 }
               else
                 {
                  //--- check
                  if(klnew<mmax)
                    {
                     //--- Eigenvalue in both intervals -- add upper to queue.
                     klnew=klnew+1;
                     ab[klnew].Set(2,ab[ji][2]);
                     nab[klnew].Set(2,nab[ji][2]);
                     ab[klnew].Set(1,tmp1);
                     nab[klnew].Set(1,itmp1);
                     ab[ji].Set(2,tmp1);
                     nab[ji].Set(2,itmp1);
                    }
                  else
                    {
                     info=mmax+1;
                     //--- exit the function
                     return;
                    }
                 }
              }
           }
         else
           {
            //--- IJOB=3: Binary search.  Keep only the interval
            //--- containing  w  s.t. N(w) = NVAL
            if(itmp1<=nval[ji])
              {
               ab[ji].Set(1,tmp1);
               nab[ji].Set(1,itmp1);
              }
            //--- check
            if(itmp1>=nval[ji])
              {
               ab[ji].Set(2,tmp1);
               nab[ji].Set(2,itmp1);
              }
           }
        }
      kl=klnew;
      //--- Check for convergence
      kfnew=kf;
      for(ji=kf;ji<=kl;ji++)
        {
         tmp1=MathAbs(ab[ji][2]-ab[ji][1]);
         tmp2=MathMax(MathAbs(ab[ji][2]),MathAbs(ab[ji][1]));
         //--- check
         if(tmp1<(double)(MathMax(abstol,MathMax(pivmin,reltol*tmp2))) || nab[ji][1]>=nab[ji][2])
           {
            //--- Converged -- Swap with position KFNEW,
            //--- then increment KFNEW
            if(ji>kfnew)
              {
               tmp1=ab[ji][1];
               tmp2=ab[ji][2];
               itmp1=nab[ji][1];
               itmp2=nab[ji][2];
               //--- change values
               ab[ji].Set(1,ab[kfnew][1]);
               ab[ji].Set(2,ab[kfnew][2]);
               nab[ji].Set(1,nab[kfnew][1]);
               nab[ji].Set(2,nab[kfnew][2]);
               ab[kfnew].Set(1,tmp1);
               ab[kfnew].Set(2,tmp2);
               nab[kfnew].Set(1,itmp1);
               nab[kfnew].Set(2,itmp2);
               //--- check
               if(ijob==3)
                 {
                  itmp1=nval[ji];
                  nval[ji]=nval[kfnew];
                  nval[kfnew]=itmp1;
                 }
              }
            kfnew=kfnew+1;
           }
        }
      kf=kfnew;
      //--- Choose Midpoints
      for(ji=kf;ji<=kl;ji++)
         c[ji]=0.5*(ab[ji][1]+ab[ji][2]);
      //--- If no more intervals to refine, quit.
      if(kf>kl)
         break;
     }
//--- Converged
   info=(int)MathMax(kl+1-kf,0);
   mout=kl;
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee,Univ. of California Berkeley,NAG Ltd.,   |
//|      Courant Institute,Argonne National Lab, and Rice University |
//|      June 30,1999                                                |
//+------------------------------------------------------------------+
static void CEigenVDetect::InternalTREVC(CMatrixDouble &t,const int n,const int side,
                                         const int howmny,bool &cvselect[],CMatrixDouble &vl,
                                         CMatrixDouble &vr,int &m,int &info)
  {
//--- create variables
   bool   allv;
   bool   bothv;
   bool   leftv;
   bool   over;
   bool   pair;
   bool   rightv;
   bool   somev;
   int    i=0;
   int    ierr=0;
   int    ii=0;
   int    ip=0;
   int    iis=0;
   int    j=0;
   int    j1=0;
   int    j2=0;
   int    jnxt=0;
   int    k=0;
   int    ki=0;
   int    n2=0;
   double beta=0;
   double bignum=0;
   double emax=0;
   double ovfl=0;
   double rec=0;
   double remax=0;
   double scl=0;
   double smin=0;
   double smlnum=0;
   double ulp=0;
   double unfl=0;
   double vcrit=0;
   double vmax=0;
   double wi=0;
   double wr=0;
   double xnorm=0;
   bool   skipflag;
   int    k1=0;
   int    k2=0;
   int    k3=0;
   int    k4=0;
   double vt=0;
   int    i_=0;
   int    i1_=0;
//--- create arrays
   double work[];
   double temp[];
   bool   rswap4[];
   bool   zswap4[];
   double civ4[];
   double crv4[];
//--- create matrix
   CMatrixDouble x;
   CMatrixDouble temp11;
   CMatrixDouble temp22;
   CMatrixDouble temp11b;
   CMatrixDouble temp21b;
   CMatrixDouble temp12b;
   CMatrixDouble temp22b;
   CMatrixInt    ipivot44;
//--- create copy
   bool vselect[];
   ArrayResizeAL(vselect,ArraySize(cvselect));
   ArrayCopy(vselect,cvselect);
//--- initialization
   m=0;
   info=0;
//--- allocation
   x.Resize(3,3);
   temp11.Resize(2,2);
   temp11b.Resize(2,2);
   temp21b.Resize(3,2);
   temp12b.Resize(2,3);
   temp22b.Resize(3,3);
   temp22.Resize(3,3);
   ArrayResizeAL(work,3*n+1);
   ArrayResizeAL(temp,n+1);
   ArrayResizeAL(rswap4,5);
   ArrayResizeAL(zswap4,5);
   ArrayResizeAL(civ4,5);
   ArrayResizeAL(crv4,5);
   ipivot44.Resize(5,5);
//--- check
   if(howmny!=1)
     {
      //--- check
      if(side==1 || side==3)
         vr.Resize(n+1,n+1);
      //--- check
      if(side==2 || side==3)
         vl.Resize(n+1,n+1);
     }
//--- Decode and test the input parameters
   bothv=side==3;
   rightv=side==1 || bothv;
   leftv=side==2 || bothv;
   allv=howmny==2;
   over=howmny==1;
   somev=howmny==3;
   info=0;
//--- check
   if(n<0)
     {
      info=-2;
      return;
     }
//--- check
   if(!rightv && !leftv)
     {
      info=-3;
      return;
     }
//--- check
   if((!allv && !over) && !somev)
     {
      info=-4;
      return;
     }
//--- Set M to the number of columns required to store the selected
//--- eigenvectors, standardize the array SELECT if necessary, and
//--- test MM.
   if(somev)
     {
      m=0;
      pair=false;
      for(j=1;j<=n;j++)
        {
         //--- check
         if(pair)
           {
            pair=false;
            vselect[j]=false;
           }
         else
           {
            //--- check
            if(j<n)
              {
               //--- check
               if(t[j+1][j]==0.0)
                 {
                  //--- check
                  if(vselect[j])
                     m=m+1;
                 }
               else
                 {
                  pair=true;
                  //--- check
                  if(vselect[j] || vselect[j+1])
                    {
                     vselect[j]=true;
                     m=m+2;
                    }
                 }
              }
            else
              {
               //--- check
               if(vselect[n])
                  m=m+1;
              }
           }
        }
     }
   else
      m=n;
//--- check
   if(n==0)
      return;
//--- Set the constants to control overflow.
   unfl=CMath::m_minrealnumber;
   ovfl=1/unfl;
   ulp=CMath::m_machineepsilon;
   smlnum=unfl*(n/ulp);
   bignum=(1-ulp)/smlnum;
//--- Compute 1-norm of each column of strictly upper triangular
//--- part of T to control overflow in triangular solver.
   work[1]=0;
   for(j=2;j<=n;j++)
     {
      work[j]=0;
      for(i=1;i<j;i++)
         work[j]=work[j]+MathAbs(t[i][j]);
     }
//--- Index IP is used to specify the real or complex eigenvalue:
//--- IP = 0, real eigenvalue,
//---      1, first of conjugate complex pair: (wr,wi)
//---     -1, second of conjugate complex pair: (wr,wi)
   n2=2*n;
//--- check
   if(rightv)
     {
      //--- Compute right eigenvectors.
      ip=0;
      iis=m;
      for(ki=n;ki>=1;ki--)
        {
         skipflag=false;
         //--- check
         if(ip==1)
            skipflag=true;
         else
           {
            //--- check
            if(ki!=1)
              {
               //--- check
               if(t[ki][ki-1]!=0.0)
                  ip=-1;
              }
            //--- check
            if(somev)
              {
               //--- check
               if(ip==0)
                 {
                  //--- check
                  if(!vselect[ki])
                     skipflag=true;
                 }
               else
                 {
                  //--- check
                  if(!vselect[ki-1])
                     skipflag=true;
                 }
              }
           }
         //--- check
         if(!skipflag)
           {
            //--- Compute the KI-th eigenvalue (WR,WI).
            wr=t[ki][ki];
            wi=0;
            //--- check
            if(ip!=0)
               wi=MathSqrt(MathAbs(t[ki][ki-1]))*MathSqrt(MathAbs(t[ki-1][ki]));
            smin=MathMax(ulp*(MathAbs(wr)+MathAbs(wi)),smlnum);
            //--- check
            if(ip==0)
              {
               //--- Real right eigenvector
               work[ki+n]=1;
               //--- Form right-hand side
               for(k=1;k<=ki-1;k++)
                  work[k+n]=-t[k][ki];
               //--- Solve the upper quasi-triangular system:
               //---   (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
               jnxt=ki-1;
               for(j=ki-1;j>=1;j--)
                 {
                  //--- check
                  if(j>jnxt)
                     continue;
                  j1=j;
                  j2=j;
                  jnxt=j-1;
                  //--- check
                  if(j>1)
                    {
                     //--- check
                     if(t[j][j-1]!=0.0)
                       {
                        j1=j-1;
                        jnxt=j-2;
                       }
                    }
                  //--- check
                  if(j1==j2)
                    {
                     //--- 1-by-1 diagonal block
                     temp11[1].Set(1,t[j][j]);
                     temp11b[1].Set(1,work[j+n]);
                     //--- function call
                     InternalHsEVDLALN2(false,1,1,smin,1,temp11,1.0,1.0,temp11b,wr,0.0,rswap4,zswap4,ipivot44,civ4,crv4,x,scl,xnorm,ierr);
                     //--- Scale X(1,1) to avoid overflow when updating
                     //--- the right-hand side.
                     if(xnorm>1.0)
                       {
                        //--- check
                        if(work[j]>bignum/xnorm)
                          {
                           x[1].Set(1,x[1][1]/xnorm);
                           scl=scl/xnorm;
                          }
                       }
                     //--- Scale if necessary
                     if(scl!=1.0)
                       {
                        k1=n+1;
                        k2=n+ki;
                        for(i_=k1;i_<=k2;i_++)
                           work[i_]=scl*work[i_];
                       }
                     work[j+n]=x[1][1];
                     //--- Update right-hand side
                     k1=1+n;
                     k2=j-1+n;
                     k3=j-1;
                     vt=-x[1][1];
                     i1_=1-k1;
                     for(i_=k1;i_<=k2;i_++)
                       {
                        work[i_]=work[i_]+vt*t[i_+i1_][j];
                       }
                    }
                  else
                    {
                     //--- 2-by-2 diagonal block
                     temp22[1].Set(1,t[j-1][j-1]);
                     temp22[1].Set(2,t[j-1][j]);
                     temp22[2].Set(1,t[j][j-1]);
                     temp22[2].Set(2,t[j][j]);
                     temp21b[1].Set(1,work[j-1+n]);
                     temp21b[2].Set(1,work[j+n]);
                     //--- function call
                     InternalHsEVDLALN2(false,2,1,smin,1.0,temp22,1.0,1.0,temp21b,wr,0,rswap4,zswap4,ipivot44,civ4,crv4,x,scl,xnorm,ierr);
                     //--- Scale X(1,1) and X(2,1) to avoid overflow when
                     //--- updating the right-hand side.
                     if(xnorm>1.0)
                       {
                        beta=MathMax(work[j-1],work[j]);
                        //--- check
                        if(beta>bignum/xnorm)
                          {
                           x[1].Set(1,x[1][1]/xnorm);
                           x[2].Set(1,x[2][1]/xnorm);
                           scl=scl/xnorm;
                          }
                       }
                     //--- Scale if necessary
                     if(scl!=1.0)
                       {
                        k1=1+n;
                        k2=ki+n;
                        for(i_=k1;i_<=k2;i_++)
                           work[i_]=scl*work[i_];
                       }
                     work[j-1+n]=x[1][1];
                     work[j+n]=x[2][1];
                     //--- Update right-hand side
                     k1=1+n;
                     k2=j-2+n;
                     k3=j-2;
                     k4=j-1;
                     vt=-x[1][1];
                     i1_=1-k1;
                     for(i_=k1;i_<=k2;i_++)
                        work[i_]=work[i_]+vt*t[i_+i1_][k4];
                     vt=-x[2][1];
                     i1_=1-k1;
                     for(i_=k1;i_<=k2;i_++)
                        work[i_]=work[i_]+vt*t[i_+i1_][j];
                    }
                 }
               //--- Copy the vector x or Q*x to VR and normalize.
               if(!over)
                 {
                  k1=1+n;
                  k2=ki+n;
                  i1_=k1-1;
                  for(i_=1;i_<=ki;i_++)
                     vr[i_].Set(iis,work[i_+i1_]);
                  //--- function call
                  ii=CBlas::ColumnIdxAbsMax(vr,1,ki,iis);
                  remax=1/MathAbs(vr[ii][iis]);
                  for(i_=1;i_<=ki;i_++)
                     vr[i_].Set(iis,remax*vr[i_][iis]);
                  for(k=ki+1;k<=n;k++)
                     vr[k].Set(iis,0);
                 }
               else
                 {
                  //--- check
                  if(ki>1)
                    {
                     for(i_=1;i_<=n;i_++)
                        temp[i_]=vr[i_][ki];
                     //--- function call
                     CBlas::MatrixVectorMultiply(vr,1,n,1,ki-1,false,work,1+n,ki-1+n,1.0,temp,1,n,work[ki+n]);
                     for(i_=1;i_<=n;i_++)
                        vr[i_].Set(ki,temp[i_]);
                    }
                  //--- function call
                  ii=CBlas::ColumnIdxAbsMax(vr,1,n,ki);
                  remax=1/MathAbs(vr[ii][ki]);
                  for(i_=1;i_<=n;i_++)
                     vr[i_].Set(ki,remax*vr[i_][ki]);
                 }
              }
            else
              {
               //--- Complex right eigenvector.
               //--- Initial solve
               //---     [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
               //---     [ (T(KI,KI-1)   T(KI,KI)   )               ]
               if(MathAbs(t[ki-1][ki])>=MathAbs(t[ki][ki-1]))
                 {
                  work[ki-1+n]=1;
                  work[ki+n2]=wi/t[ki-1][ki];
                 }
               else
                 {
                  work[ki-1+n]=-(wi/t[ki][ki-1]);
                  work[ki+n2]=1;
                 }
               work[ki+n]=0;
               work[ki-1+n2]=0;
               //--- Form right-hand side
               for(k=1;k<=ki-2;k++)
                 {
                  work[k+n]=-(work[ki-1+n]*t[k][ki-1]);
                  work[k+n2]=-(work[ki+n2]*t[k][ki]);
                 }
               //--- Solve upper quasi-triangular system:
               //--- (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
               jnxt=ki-2;
               for(j=ki-2;j>=1;j--)
                 {
                  //--- check
                  if(j>jnxt)
                     continue;
                  j1=j;
                  j2=j;
                  jnxt=j-1;
                  //--- check
                  if(j>1)
                    {
                     //--- check
                     if(t[j][j-1]!=0.0)
                       {
                        j1=j-1;
                        jnxt=j-2;
                       }
                    }
                  //--- check
                  if(j1==j2)
                    {
                     //--- 1-by-1 diagonal block
                     temp11[1].Set(1,t[j][j]);
                     temp12b[1].Set(1,work[j+n]);
                     temp12b[1].Set(2,work[j+n+n]);
                     //--- function call
                     InternalHsEVDLALN2(false,1,2,smin,1.0,temp11,1.0,1.0,temp12b,wr,wi,rswap4,zswap4,ipivot44,civ4,crv4,x,scl,xnorm,ierr);
                     //--- Scale X(1,1) and X(1,2) to avoid overflow when
                     //--- updating the right-hand side.
                     if(xnorm>1.0)
                       {
                        //--- check
                        if(work[j]>bignum/xnorm)
                          {
                           x[1].Set(1,x[1][1]/xnorm);
                           x[1].Set(2,x[1][2]/xnorm);
                           scl=scl/xnorm;
                          }
                       }
                     //--- Scale if necessary
                     if(scl!=1.0)
                       {
                        k1=1+n;
                        k2=ki+n;
                        for(i_=k1;i_<=k2;i_++)
                           work[i_]=scl*work[i_];
                        k1=1+n2;
                        k2=ki+n2;
                        for(i_=k1;i_<=k2;i_++)
                           work[i_]=scl*work[i_];
                       }
                     work[j+n]=x[1][1];
                     work[j+n2]=x[1][2];
                     //--- Update the right-hand side
                     k1=1+n;
                     k2=j-1+n;
                     k3=1;
                     k4=j-1;
                     vt=-x[1][1];
                     i1_=k3-k1;
                     for(i_=k1;i_<=k2;i_++)
                        work[i_]=work[i_]+vt*t[i_+i1_][j];
                     //--- change values
                     k1=1+n2;
                     k2=j-1+n2;
                     k3=1;
                     k4=j-1;
                     vt=-x[1][2];
                     i1_=k3-k1;
                     for(i_=k1;i_<=k2;i_++)
                        work[i_]=work[i_]+vt*t[i_+i1_][j];
                    }
                  else
                    {
                     //--- 2-by-2 diagonal block
                     temp22[1].Set(1,t[j-1][j-1]);
                     temp22[1].Set(2,t[j-1][j]);
                     temp22[2].Set(1,t[j][j-1]);
                     temp22[2].Set(2,t[j][j]);
                     temp22b[1].Set(1,work[j-1+n]);
                     temp22b[1].Set(2,work[j-1+n+n]);
                     temp22b[2].Set(1,work[j+n]);
                     temp22b[2].Set(2,work[j+n+n]);
                     //--- function call
                     InternalHsEVDLALN2(false,2,2,smin,1.0,temp22,1.0,1.0,temp22b,wr,wi,rswap4,zswap4,ipivot44,civ4,crv4,x,scl,xnorm,ierr);
                     //--- Scale X to avoid overflow when updating
                     //--- the right-hand side.
                     if(xnorm>1.0)
                       {
                        beta=MathMax(work[j-1],work[j]);
                        //--- check
                        if(beta>bignum/xnorm)
                          {
                           rec=1/xnorm;
                           x[1].Set(1,x[1][1]*rec);
                           x[1].Set(2,x[1][2]*rec);
                           x[2].Set(1,x[2][1]*rec);
                           x[2].Set(2,x[2][2]*rec);
                           scl=scl*rec;
                          }
                       }
                     //--- Scale if necessary
                     if(scl!=1.0)
                       {
                        for(i_=1+n;i_<=ki+n;i_++)
                           work[i_]=scl*work[i_];
                        for(i_=1+n2;i_<=ki+n2;i_++)
                           work[i_]=scl*work[i_];
                       }
                     //--- change values
                     work[j-1+n]=x[1][1];
                     work[j+n]=x[2][1];
                     work[j-1+n2]=x[1][2];
                     work[j+n2]=x[2][2];
                     //--- Update the right-hand side
                     vt=-x[1][1];
                     i1_=-n;
                     for(i_=n+1;i_<=n+j-2;i_++)
                        work[i_]=work[i_]+vt*t[i_+i1_][j-1];
                     vt=-x[2][1];
                     i1_=-n;
                     for(i_=n+1;i_<=n+j-2;i_++)
                        work[i_]=work[i_]+vt*t[i_+i1_][j];
                     vt=-x[1][2];
                     i1_=-n2;
                     for(i_=n2+1;i_<=n2+j-2;i_++)
                        work[i_]=work[i_]+vt*t[i_+i1_][j-1];
                     vt=-x[2][2];
                     i1_=-n2;
                     for(i_=n2+1;i_<=n2+j-2;i_++)
                        work[i_]=work[i_]+vt*t[i_+i1_][j];
                    }
                 }
               //--- Copy the vector x or Q*x to VR and normalize.
               if(!over)
                 {
                  i1_=n;
                  for(i_=1;i_<=ki;i_++)
                     vr[i_].Set(iis-1,work[i_+i1_]);
                  i1_=n2;
                  for(i_=1;i_<=ki;i_++)
                     vr[i_].Set(iis,work[i_+i1_]);
                  emax=0;
                  for(k=1;k<=ki;k++)
                     emax=MathMax(emax,MathAbs(vr[k][iis-1])+MathAbs(vr[k][iis]));
                  remax=1/emax;
                  //--- copy
                  for(i_=1;i_<=ki;i_++)
                     vr[i_].Set(iis-1,remax*vr[i_][iis-1]);
                  for(i_=1;i_<=ki;i_++)
                     vr[i_].Set(iis,remax*vr[i_][iis]);
                  for(k=ki+1;k<=n;k++)
                     vr[k].Set(iis-1,0);
                  vr[k].Set(iis,0);
                 }
               else
                 {
                  //--- check
                  if(ki>2)
                    {
                     for(i_=1;i_<=n;i_++)
                        temp[i_]=vr[i_][ki-1];
                     //--- function call
                     CBlas::MatrixVectorMultiply(vr,1,n,1,ki-2,false,work,1+n,ki-2+n,1.0,temp,1,n,work[ki-1+n]);
                     for(i_=1;i_<=n;i_++)
                        vr[i_].Set(ki-1,temp[i_]);
                     for(i_=1;i_<=n;i_++)
                        temp[i_]=vr[i_][ki];
                     //--- function call
                     CBlas::MatrixVectorMultiply(vr,1,n,1,ki-2,false,work,1+n2,ki-2+n2,1.0,temp,1,n,work[ki+n2]);
                     for(i_=1;i_<=n;i_++)
                        vr[i_].Set(ki,temp[i_]);
                    }
                  else
                    {
                     vt=work[ki-1+n];
                     //--- copy
                     for(i_=1;i_<=n;i_++)
                        vr[i_].Set(ki-1,vt*vr[i_][ki-1]);
                     vt=work[ki+n2];
                     for(i_=1;i_<=n;i_++)
                        vr[i_].Set(ki,vt*vr[i_][ki]);
                    }
                  emax=0;
                  for(k=1;k<=n;k++)
                     emax=MathMax(emax,MathAbs(vr[k][ki-1])+MathAbs(vr[k][ki]));
                  remax=1/emax;
                  //--- copy
                  for(i_=1;i_<=n;i_++)
                     vr[i_].Set(ki-1,remax*vr[i_][ki-1]);
                  for(i_=1;i_<=n;i_++)
                     vr[i_].Set(ki,remax*vr[i_][ki]);
                 }
              }
            iis=iis-1;
            //--- check
            if(ip!=0)
               iis=iis-1;
           }
         //--- check
         if(ip==1)
            ip=0;
         //--- check
         if(ip==-1)
            ip=1;
        }
     }
//--- check
   if(leftv)
     {
      //--- Compute left eigenvectors.
      ip=0;
      iis=1;
      for(ki=1;ki<=n;ki++)
        {
         skipflag=false;
         //--- check
         if(ip==-1)
            skipflag=true;
         else
           {
            //--- check
            if(ki!=n)
              {
               //--- check
               if(t[ki+1][ki]!=0.0)
                  ip=1;
              }
            //--- check
            if(somev)
              {
               //--- check
               if(!vselect[ki])
                  skipflag=true;
              }
           }
         //--- check
         if(!skipflag)
           {
            //--- Compute the KI-th eigenvalue (WR,WI).
            wr=t[ki][ki];
            wi=0;
            //--- check
            if(ip!=0)
               wi=MathSqrt(MathAbs(t[ki][ki+1]))*MathSqrt(MathAbs(t[ki+1][ki]));
            smin=MathMax(ulp*(MathAbs(wr)+MathAbs(wi)),smlnum);
            //--- check
            if(ip==0)
              {
               //--- Real left eigenvector.
               work[ki+n]=1;
               //--- Form right-hand side
               for(k=ki+1;k<=n;k++)
                  work[k+n]=-t[ki][k];
               //--- Solve the quasi-triangular system:
               //--- (T(KI+1:N,KI+1:N) - WR)'*X = SCALE*WORK
               vmax=1;
               vcrit=bignum;
               jnxt=ki+1;
               for(j=ki+1;j<=n;j++)
                 {
                  //--- check
                  if(j<jnxt)
                     continue;
                  j1=j;
                  j2=j;
                  jnxt=j+1;
                  //--- check
                  if(j<n)
                    {
                     //--- check
                     if(t[j+1][j]!=0.0)
                       {
                        j2=j+1;
                        jnxt=j+2;
                       }
                    }
                  //--- check
                  if(j1==j2)
                    {
                     //--- 1-by-1 diagonal block
                     //--- Scale if necessary to avoid overflow when forming
                     //--- the right-hand side.
                     if(work[j]>vcrit)
                       {
                        rec=1/vmax;
                        for(i_=ki+n;i_<=n+n;i_++)
                          {
                           work[i_]=rec*work[i_];
                          }
                        vmax=1;
                        vcrit=bignum;
                       }
                     i1_=n;
                     vt=0.0;
                     for(i_=ki+1;i_<j;i_++)
                        vt+=t[i_][j]*work[i_+i1_];
                     work[j+n]=work[j+n]-vt;
                     //--- Solve (T(J,J)-WR)'*X = WORK
                     temp11[1].Set(1,t[j][j]);
                     temp11b[1].Set(1,work[j+n]);
                     InternalHsEVDLALN2(false,1,1,smin,1.0,temp11,1.0,1.0,temp11b,wr,0,rswap4,zswap4,ipivot44,civ4,crv4,x,scl,xnorm,ierr);
                     //--- Scale if necessary
                     if(scl!=1.0)
                       {
                        for(i_=ki+n;i_<=n+n;i_++)
                           work[i_]=scl*work[i_];
                       }
                     work[j+n]=x[1][1];
                     vmax=MathMax(MathAbs(work[j+n]),vmax);
                     vcrit=bignum/vmax;
                    }
                  else
                    {
                     //--- 2-by-2 diagonal block
                     //--- Scale if necessary to avoid overflow when forming
                     //--- the right-hand side.
                     beta=MathMax(work[j],work[j+1]);
                     //--- check
                     if(beta>vcrit)
                       {
                        rec=1/vmax;
                        for(i_=ki+n;i_<=n+n;i_++)
                           work[i_]=rec*work[i_];
                        vmax=1;
                        vcrit=bignum;
                       }
                     i1_=n;
                     vt=0.0;
                     for(i_=ki+1;i_<j;i_++)
                        vt+=t[i_][j]*work[i_+i1_];
                     //--- change values
                     work[j+n]=work[j+n]-vt;
                     i1_=n;
                     vt=0.0;
                     for(i_=ki+1;i_<j;i_++)
                        vt+=t[i_][j+1]*work[i_+i1_];
                     work[j+1+n]=work[j+1+n]-vt;
                     //--- Solve
                     //---    [T(J,J)-WR   T(J,J+1)     ]'* X = SCALE*( WORK1 )
                     //---    [T(J+1,J)    T(J+1,J+1)-WR]             ( WORK2 )
                     temp22[1].Set(1,t[j][j]);
                     temp22[1].Set(2,t[j][j+1]);
                     temp22[2].Set(1,t[j+1][j]);
                     temp22[2].Set(2,t[j+1][j+1]);
                     temp21b[1].Set(1,work[j+n]);
                     temp21b[2].Set(1,work[j+1+n]);
                     //--- function call
                     InternalHsEVDLALN2(true,2,1,smin,1.0,temp22,1.0,1.0,temp21b,wr,0,rswap4,zswap4,ipivot44,civ4,crv4,x,scl,xnorm,ierr);
                     //--- Scale if necessary
                     if(scl!=1.0)
                       {
                        for(i_=ki+n;i_<=n+n;i_++)
                           work[i_]=scl*work[i_];
                       }
                     //--- change values
                     work[j+n]=x[1][1];
                     work[j+1+n]=x[2][1];
                     vmax=MathMax(MathAbs(work[j+n]),MathMax(MathAbs(work[j+1+n]),vmax));
                     vcrit=bignum/vmax;
                    }
                 }
               //--- Copy the vector x or Q*x to VL and normalize.
               if(!over)
                 {
                  i1_=n;;
                  for(i_=ki;i_<=n;i_++)
                     vl[i_].Set(iis,work[i_+i1_]);
                  //--- function call
                  ii=CBlas::ColumnIdxAbsMax(vl,ki,n,iis);
                  remax=1/MathAbs(vl[ii][iis]);
                  for(i_=ki;i_<=n;i_++)
                     vl[i_].Set(iis,remax*vl[i_][iis]);
                  for(k=1;k<=ki-1;k++)
                     vl[k].Set(iis,0);
                 }
               else
                 {
                  //--- check
                  if(ki<n)
                    {
                     for(i_=1;i_<=n;i_++)
                        temp[i_]=vl[i_][ki];
                     //--- function call
                     CBlas::MatrixVectorMultiply(vl,1,n,ki+1,n,false,work,ki+1+n,n+n,1.0,temp,1,n,work[ki+n]);
                     for(i_=1;i_<=n;i_++)
                        vl[i_].Set(ki,temp[i_]);
                    }
                  //--- function call
                  ii=CBlas::ColumnIdxAbsMax(vl,1,n,ki);
                  remax=1/MathAbs(vl[ii][ki]);
                  for(i_=1;i_<=n;i_++)
                     vl[i_].Set(ki,remax*vl[i_][ki]);
                 }
              }
            else
              {
               //--- Complex left eigenvector.
               //--- Initial solve:
               //---   ((T(KI,KI)    T(KI,KI+1) )' - (WR - I* WI))*X = 0.
               //---   ((T(KI+1,KI) T(KI+1,KI+1))                )
               if(MathAbs(t[ki][ki+1])>=MathAbs(t[ki+1][ki]))
                 {
                  work[ki+n]=wi/t[ki][ki+1];
                  work[ki+1+n2]=1;
                 }
               else
                 {
                  work[ki+n]=1;
                  work[ki+1+n2]=-(wi/t[ki+1][ki]);
                 }
               work[ki+1+n]=0;
               work[ki+n2]=0;
               //--- Form right-hand side
               for(k=ki+2;k<=n;k++)
                 {
                  work[k+n]=-(work[ki+n]*t[ki][k]);
                  work[k+n2]=-(work[ki+1+n2]*t[ki+1][k]);
                 }
               //--- Solve complex quasi-triangular system:
               //--- ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2
               vmax=1;
               vcrit=bignum;
               jnxt=ki+2;
               for(j=ki+2;j<=n;j++)
                 {
                  //--- check
                  if(j<jnxt)
                     continue;
                  j1=j;
                  j2=j;
                  jnxt=j+1;
                  //--- check
                  if(j<n)
                    {
                     //--- check
                     if(t[j+1][j]!=0.0)
                       {
                        j2=j+1;
                        jnxt=j+2;
                       }
                    }
                  //--- check
                  if(j1==j2)
                    {
                     //--- 1-by-1 diagonal block
                     //--- Scale if necessary to avoid overflow when
                     //--- forming the right-hand side elements.
                     if(work[j]>vcrit)
                       {
                        rec=1/vmax;
                        for(i_=ki+n;i_<=n+n;i_++)
                           work[i_]=rec*work[i_];
                        for(i_=ki+n2;i_<=n+n2;i_++)
                           work[i_]=rec*work[i_];
                        vmax=1;
                        vcrit=bignum;
                       }
                     i1_=n;
                     vt=0.0;
                     for(i_=ki+2;i_<j;i_++)
                        vt+=t[i_][j]*work[i_+i1_];
                     //--- calculation
                     work[j+n]=work[j+n]-vt;
                     i1_=n2;
                     vt=0.0;
                     for(i_=ki+2;i_<j;i_++)
                        vt+=t[i_][j]*work[i_+i1_];
                     work[j+n2]=work[j+n2]-vt;
                     //--- Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2
                     temp11[1].Set(1,t[j][j]);
                     temp12b[1].Set(1,work[j+n]);
                     temp12b[1].Set(2,work[j+n+n]);
                     //--- function call
                     InternalHsEVDLALN2(false,1,2,smin,1.0,temp11,1.0,1.0,temp12b,wr,-wi,rswap4,zswap4,ipivot44,civ4,crv4,x,scl,xnorm,ierr);
                     //--- Scale if necessary
                     if(scl!=1.0)
                       {
                        for(i_=ki+n;i_<=n+n;i_++)
                           work[i_]=scl*work[i_];
                        for(i_=ki+n2;i_<=n+n2;i_++)
                           work[i_]=scl*work[i_];
                       }
                     //--- change values
                     work[j+n]=x[1][1];
                     work[j+n2]=x[1][2];
                     vmax=MathMax(MathAbs(work[j+n]),MathMax(MathAbs(work[j+n2]),vmax));
                     vcrit=bignum/vmax;
                    }
                  else
                    {
                     //--- 2-by-2 diagonal block
                     //--- Scale if necessary to avoid overflow when forming
                     //--- the right-hand side elements.
                     beta=MathMax(work[j],work[j+1]);
                     //--- check
                     if(beta>vcrit)
                       {
                        rec=1/vmax;
                        for(i_=ki+n;i_<=n+n;i_++)
                           work[i_]=rec*work[i_];
                        for(i_=ki+n2;i_<=n+n2;i_++)
                           work[i_]=rec*work[i_];
                        vmax=1;
                        vcrit=bignum;
                       }
                     i1_=n;
                     vt=0.0;
                     for(i_=ki+2;i_<j;i_++)
                        vt+=t[i_][j]*work[i_+i1_];
                     //--- calculation
                     work[j+n]=work[j+n]-vt;
                     i1_=n2;
                     vt=0.0;
                     for(i_=ki+2;i_<j;i_++)
                        vt+=t[i_][j]*work[i_+i1_];
                     //--- calculation
                     work[j+n2]=work[j+n2]-vt;
                     i1_=n;
                     vt=0.0;
                     for(i_=ki+2;i_<j;i_++)
                        vt+=t[i_][j+1]*work[i_+i1_];
                     //--- calculation
                     work[j+1+n]=work[j+1+n]-vt;
                     i1_=n2;
                     vt=0.0;
                     for(i_=ki+2;i_<j;i_++)
                        vt+=t[i_][j+1]*work[i_+i1_];
                     work[j+1+n2]=work[j+1+n2]-vt;
                     //--- Solve 2-by-2 complex linear equation
                     //---   ([T(j,j)   T(j,j+1)  ]'-(wr-i*wi)*I)*X = SCALE*B
                     //---   ([T(j+1,j) T(j+1,j+1)]             )
                     temp22[1].Set(1,t[j][j]);
                     temp22[1].Set(2,t[j][j+1]);
                     temp22[2].Set(1,t[j+1][j]);
                     temp22[2].Set(2,t[j+1][j+1]);
                     temp22b[1].Set(1,work[j+n]);
                     temp22b[1].Set(2,work[j+n+n]);
                     temp22b[2].Set(1,work[j+1+n]);
                     temp22b[2].Set(2,work[j+1+n+n]);
                     //--- function call
                     InternalHsEVDLALN2(true,2,2,smin,1.0,temp22,1.0,1.0,temp22b,wr,-wi,rswap4,zswap4,ipivot44,civ4,crv4,x,scl,xnorm,ierr);
                     //--- Scale if necessary
                     if(scl!=1.0)
                       {
                        for(i_=ki+n;i_<=n+n;i_++)
                           work[i_]=scl*work[i_];
                        for(i_=ki+n2;i_<=n+n2;i_++)
                           work[i_]=scl*work[i_];
                       }
                     //--- change values
                     work[j+n]=x[1][1];
                     work[j+n2]=x[1][2];
                     work[j+1+n]=x[2][1];
                     work[j+1+n2]=x[2][2];
                     vmax=MathMax(MathAbs(x[1][1]),vmax);
                     vmax=MathMax(MathAbs(x[1][2]),vmax);
                     vmax=MathMax(MathAbs(x[2][1]),vmax);
                     vmax=MathMax(MathAbs(x[2][2]),vmax);
                     vcrit=bignum/vmax;
                    }
                 }
               //--- Copy the vector x or Q*x to VL and normalize.
               if(!over)
                 {
                  i1_=n;
                  for(i_=ki;i_<=n;i_++)
                     vl[i_].Set(iis,work[i_+i1_]);
                  i1_=n2;
                  for(i_=ki;i_<=n;i_++)
                     vl[i_].Set(iis+1,work[i_+i1_]);
                  emax=0;
                  for(k=ki;k<=n;k++)
                     emax=MathMax(emax,MathAbs(vl[k][iis])+MathAbs(vl[k][iis+1]));
                  remax=1/emax;
                  //--- copy
                  for(i_=ki;i_<=n;i_++)
                     vl[i_].Set(iis,remax*vl[i_][iis]);
                  for(i_=ki;i_<=n;i_++)
                     vl[i_].Set(iis+1,remax*vl[i_][iis+1]);
                  for(k=1;k<=ki-1;k++)
                    {
                     vl[k].Set(iis,0);
                     vl[k].Set(iis+1,0);
                    }
                 }
               else
                 {
                  //--- check
                  if(ki<n-1)
                    {
                     for(i_=1;i_<=n;i_++)
                        temp[i_]=vl[i_][ki];
                     //--- function call
                     CBlas::MatrixVectorMultiply(vl,1,n,ki+2,n,false,work,ki+2+n,n+n,1.0,temp,1,n,work[ki+n]);
                     for(i_=1;i_<=n;i_++)
                        vl[i_].Set(ki,temp[i_]);
                     for(i_=1;i_<=n;i_++)
                        temp[i_]=vl[i_][ki+1];
                     //--- function call
                     CBlas::MatrixVectorMultiply(vl,1,n,ki+2,n,false,work,ki+2+n2,n+n2,1.0,temp,1,n,work[ki+1+n2]);
                     for(i_=1;i_<=n;i_++)
                        vl[i_].Set(ki+1,temp[i_]);
                    }
                  else
                    {
                     //--- copy
                     vt=work[ki+n];
                     for(i_=1;i_<=n;i_++)
                        vl[i_].Set(ki,vt*vl[i_][ki]);
                     vt=work[ki+1+n2];
                     for(i_=1;i_<=n;i_++)
                        vl[i_].Set(ki+1,vt*vl[i_][ki+1]);
                    }
                  emax=0;
                  for(k=1;k<=n;k++)
                     emax=MathMax(emax,MathAbs(vl[k][ki])+MathAbs(vl[k][ki+1]));
                  remax=1/emax;
                  //--- copy
                  for(i_=1;i_<=n;i_++)
                     vl[i_].Set(ki,remax*vl[i_][ki]);
                  for(i_=1;i_<=n;i_++)
                     vl[i_].Set(ki+1,remax*vl[i_][ki+1]);
                 }
              }
            iis=iis+1;
            //--- check
            if(ip!=0)
               iis=iis+1;
           }
         //--- check
         if(ip==-1)
            ip=0;
         //--- check
         if(ip==1)
            ip=-1;
        }
     }
  }
//+------------------------------------------------------------------+
//| DLALN2 solves a system of the form  (ca A - w D ) X = s B        |
//| or (ca A' - w D) X = s B   with possible scaling ("s") and       |
//| perturbation of A. (A' means A-transpose.)                       |
//| A is an NA x NA real matrix, ca is a real scalar, D is an NA x   |
//| NA real diagonal matrix, w is a real or complex value, and X and |
//| B are NA x 1 matrices -- real if w is real, complex if w is      |
//| complex.  NA may be 1 or 2.                                      |
//| If w is complex, X and B are represented as NA x 2 matrices,     |
//| the first column of each being the real part and the second      |
//| being the imaginary part.                                        |
//| "s" is a scaling factor (.LE. 1), computed by DLALN2, which is   |
//| so chosen that X can be computed without overflow.  X is further |
//| scaled if necessary to assure that norm(ca A - w D)*norm(X) is   |
//| less than overflow.                                              |
//| If both singular values of (ca A - w D) are less than SMIN,      |
//| SMIN*identity will be used instead of (ca A - w D).  If only one |
//| singular value is less than SMIN, one element of (ca A - w D)    |
//| will be perturbed enough to make the smallest singular value     |
//| roughly SMIN. If both singular values are at least SMIN,         |
//| (ca A - w D) will not be perturbed.  In any case, the            |
//| perturbation will be at most some small multiple of max( SMIN,   |
//| ulp*norm(ca A - w D) ).  The singular values are computed by     |
//| infinity-norm approximations, and thus will only be correct to a |
//| factor of 2 or so.                                               |
//| Note: all input quantities are assumed to be smaller than        |
//| overflow by a reasonable factor. (See BIGNUM.)                   |
//|   -- LAPACK auxiliary routine (version 3.0) --                   |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      October 31, 1992                                            |
//+------------------------------------------------------------------+
static void CEigenVDetect::InternalHsEVDLALN2(const bool ltrans,const int na,
                                              const int nw,const double smin,
                                              const double ca,CMatrixDouble &a,
                                              const double d1,const double d2,
                                              CMatrixDouble &b,const double wr,
                                              const double wi,bool &rswap4[],
                                              bool &zswap4[],CMatrixInt &ipivot44,
                                              double &civ4[],double &crv4[],
                                              CMatrixDouble &x,double &scl,
                                              double &xnorm,int &info)
  {
//--- create variables
   int    icmax=0;
   int    j=0;
   double bbnd=0;
   double bi1=0;
   double bi2=0;
   double bignum=0;
   double bnorm=0;
   double br1=0;
   double br2=0;
   double ci21=0;
   double ci22=0;
   double cmax=0;
   double cnorm=0;
   double cr21=0;
   double cr22=0;
   double csi=0;
   double csr=0;
   double li21=0;
   double lr21=0;
   double smini=0;
   double smlnum=0;
   double temp=0;
   double u22abs=0;
   double ui11=0;
   double ui11r=0;
   double ui12=0;
   double ui12s=0;
   double ui22=0;
   double ur11=0;
   double ur11r=0;
   double ur12=0;
   double ur12s=0;
   double ur22=0;
   double xi1=0;
   double xi2=0;
   double xr1=0;
   double xr2=0;
   double tmp1=0;
   double tmp2=0;
//--- initialization
   scl=0;
   xnorm=0;
   info=0;
   zswap4[1]=false;
   zswap4[2]=false;
   zswap4[3]=true;
   zswap4[4]=true;
   rswap4[1]=false;
   rswap4[2]=true;
   rswap4[3]=false;
   rswap4[4]=true;
   ipivot44[1].Set(1,1);
   ipivot44[2].Set(1,2);
   ipivot44[3].Set(1,3);
   ipivot44[4].Set(1,4);
   ipivot44[1].Set(2,2);
   ipivot44[2].Set(2,1);
   ipivot44[3].Set(2,4);
   ipivot44[4].Set(2,3);
   ipivot44[1].Set(3,3);
   ipivot44[2].Set(3,4);
   ipivot44[3].Set(3,1);
   ipivot44[4].Set(3,2);
   ipivot44[1].Set(4,4);
   ipivot44[2].Set(4,3);
   ipivot44[3].Set(4,2);
   ipivot44[4].Set(4,1);
   smlnum=2*CMath::m_minrealnumber;
   bignum=1/smlnum;
   smini=MathMax(smin,smlnum);
//--- initialization
   info=0;
   scl=1;
//--- check
   if(na==1)
     {
      //--- 1 x 1  (i.e., scalar) system   C X = B
      if(nw==1)
        {
         //--- Real 1x1 system.
         //--- C = ca A - w D
         csr=ca*a[1][1]-wr*d1;
         cnorm=MathAbs(csr);
         //--- If | C | < SMINI, use C = SMINI
         if(cnorm<smini)
           {
            csr=smini;
            cnorm=smini;
            info=1;
           }
         //--- Check scaling for  X = B / C
         bnorm=MathAbs(b[1][1]);
         //--- check
         if(cnorm<1.0 && bnorm>1.0)
           {
            //--- check
            if(bnorm>bignum*cnorm)
               scl=1/bnorm;
           }
         //--- Compute X
         x[1].Set(1,b[1][1]*scl/csr);
         xnorm=MathAbs(x[1][1]);
        }
      else
        {
         //--- Complex 1x1 system (w is complex)
         //--- C = ca A - w D
         csr=ca*a[1][1]-wr*d1;
         csi=-(wi*d1);
         cnorm=MathAbs(csr)+MathAbs(csi);
         //--- If | C | < SMINI, use C = SMINI
         if(cnorm<smini)
           {
            csr=smini;
            csi=0;
            cnorm=smini;
            info=1;
           }
         //--- Check scaling for  X = B / C
         bnorm=MathAbs(b[1][1])+MathAbs(b[1][2]);
         //--- check
         if(cnorm<1.0 && bnorm>1.0)
           {
            //--- check
            if(bnorm>bignum*cnorm)
               scl=1/bnorm;
           }
         //--- Compute X
         InternalHsEVDLADIV(scl*b[1][1],scl*b[1][2],csr,csi,tmp1,tmp2);
         x[1].Set(1,tmp1);
         x[1].Set(2,tmp2);
         xnorm=MathAbs(x[1][1])+MathAbs(x[1][2]);
        }
     }
   else
     {
      //--- 2x2 System
      //--- Compute the real part of  C = ca A - w D  (or  ca A' - w D )
      crv4[1+0]=ca*a[1][1]-wr*d1;
      crv4[2+2]=ca*a[2][2]-wr*d2;
      //--- check
      if(ltrans)
        {
         crv4[1+2]=ca*a[2][1];
         crv4[2+0]=ca*a[1][2];
        }
      else
        {
         crv4[2+0]=ca*a[2][1];
         crv4[1+2]=ca*a[1][2];
        }
      //--- check
      if(nw==1)
        {
         //--- Real 2x2 system  (w is real)
         //--- Find the largest element in C
         cmax=0;
         icmax=0;
         for(j=1;j<=4;j++)
           {
            //--- check
            if(MathAbs(crv4[j])>cmax)
              {
               cmax=MathAbs(crv4[j]);
               icmax=j;
              }
           }
         //--- If norm(C) < SMINI, use SMINI*identity.
         if(cmax<smini)
           {
            bnorm=MathMax(MathAbs(b[1][1]),MathAbs(b[2][1]));
            //--- check
            if(smini<1.0 && bnorm>1.0)
              {
               //--- check
               if(bnorm>bignum*smini)
                  scl=1/bnorm;
              }
            //--- change values
            temp=scl/smini;
            x[1].Set(1,temp*b[1][1]);
            x[2].Set(1,temp*b[2][1]);
            xnorm=temp*bnorm;
            info=1;
            //--- exit the function
            return;
           }
         //--- Gaussian elimination with complete pivoting.
         ur11=crv4[icmax];
         cr21=crv4[ipivot44[2][icmax]];
         ur12=crv4[ipivot44[3][icmax]];
         cr22=crv4[ipivot44[4][icmax]];
         ur11r=1/ur11;
         lr21=ur11r*cr21;
         ur22=cr22-ur12*lr21;
         //--- If smaller pivot < SMINI, use SMINI
         if(MathAbs(ur22)<smini)
           {
            ur22=smini;
            info=1;
           }
         //--- check
         if(rswap4[icmax])
           {
            br1=b[2][1];
            br2=b[1][1];
           }
         else
           {
            br1=b[1][1];
            br2=b[2][1];
           }
         br2=br2-lr21*br1;
         bbnd=MathMax(MathAbs(br1*(ur22*ur11r)),MathAbs(br2));
         //--- check
         if(bbnd>1.0 && MathAbs(ur22)<1.0)
           {
            //--- check
            if(bbnd>=bignum*MathAbs(ur22))
               scl=1/bbnd;
           }
         xr2=br2*scl/ur22;
         xr1=scl*br1*ur11r-xr2*(ur11r*ur12);
         //--- check
         if(zswap4[icmax])
           {
            x[1].Set(1,xr2);
            x[2].Set(1,xr1);
           }
         else
           {
            x[1].Set(1,xr1);
            x[2].Set(1,xr2);
           }
         xnorm=MathMax(MathAbs(xr1),MathAbs(xr2));
         //--- Further scaling if  norm(A) norm(X) > overflow
         if(xnorm>1.0 && cmax>1.0)
           {
            //--- check
            if(xnorm>bignum/cmax)
              {
               temp=cmax/bignum;
               x[1].Set(1,temp*x[1][1]);
               x[2].Set(1,temp*x[2][1]);
               xnorm=temp*xnorm;
               scl=temp*scl;
              }
           }
        }
      else
        {
         //--- Complex 2x2 system  (w is complex)
         //--- Find the largest element in C
         civ4[1+0]=-(wi*d1);
         civ4[2+0]=0;
         civ4[1+2]=0;
         civ4[2+2]=-(wi*d2);
         cmax=0;
         icmax=0;
         for(j=1;j<=4;j++)
           {
            //--- check
            if(MathAbs(crv4[j])+MathAbs(civ4[j])>cmax)
              {
               cmax=MathAbs(crv4[j])+MathAbs(civ4[j]);
               icmax=j;
              }
           }
         //--- If norm(C) < SMINI, use SMINI*identity.
         if(cmax<smini)
           {
            bnorm=MathMax(MathAbs(b[1][1])+MathAbs(b[1][2]),MathAbs(b[2][1])+MathAbs(b[2][2]));
            //--- check
            if(smini<1.0 && bnorm>1.0)
              {
               //--- check
               if(bnorm>bignum*smini)
                  scl=1/bnorm;
              }
            //--- change values
            temp=scl/smini;
            x[1].Set(1,temp*b[1][1]);
            x[2].Set(1,temp*b[2][1]);
            x[1].Set(2,temp*b[1][2]);
            x[2].Set(2,temp*b[2][2]);
            xnorm=temp*bnorm;
            info=1;
            //--- exit the function
            return;
           }
         //--- Gaussian elimination with complete pivoting.
         ur11=crv4[icmax];
         ui11=civ4[icmax];
         cr21=crv4[ipivot44[2][icmax]];
         ci21=civ4[ipivot44[2][icmax]];
         ur12=crv4[ipivot44[3][icmax]];
         ui12=civ4[ipivot44[3][icmax]];
         cr22=crv4[ipivot44[4][icmax]];
         ci22=civ4[ipivot44[4][icmax]];
         //--- check
         if(icmax==1 || icmax==4)
           {
            //--- Code when off-diagonals of pivoted C are real
            if(MathAbs(ur11)>MathAbs(ui11))
              {
               temp=ui11/ur11;
               ur11r=1/(ur11*(1+CMath::Sqr(temp)));
               ui11r=-(temp*ur11r);
              }
            else
              {
               temp=ur11/ui11;
               ui11r=-(1/(ui11*(1+CMath::Sqr(temp))));
               ur11r=-(temp*ui11r);
              }
            //--- change values
            lr21=cr21*ur11r;
            li21=cr21*ui11r;
            ur12s=ur12*ur11r;
            ui12s=ur12*ui11r;
            ur22=cr22-ur12*lr21;
            ui22=ci22-ur12*li21;
           }
         else
           {
            //--- Code when diagonals of pivoted C are real
            ur11r=1/ur11;
            ui11r=0;
            lr21=cr21*ur11r;
            li21=ci21*ur11r;
            ur12s=ur12*ur11r;
            ui12s=ui12*ur11r;
            ur22=cr22-ur12*lr21+ui12*li21;
            ui22=-(ur12*li21)-ui12*lr21;
           }
         u22abs=MathAbs(ur22)+MathAbs(ui22);
         //--- If smaller pivot < SMINI, use SMINI
         if(u22abs<smini)
           {
            ur22=smini;
            ui22=0;
            info=1;
           }
         //--- check
         if(rswap4[icmax])
           {
            br2=b[1][1];
            br1=b[2][1];
            bi2=b[1][2];
            bi1=b[2][2];
           }
         else
           {
            br1=b[1][1];
            br2=b[2][1];
            bi1=b[1][2];
            bi2=b[2][2];
           }
         br2=br2-lr21*br1+li21*bi1;
         bi2=bi2-li21*br1-lr21*bi1;
         bbnd=MathMax((MathAbs(br1)+MathAbs(bi1))*(u22abs*(MathAbs(ur11r)+MathAbs(ui11r))),MathAbs(br2)+MathAbs(bi2));
         //--- check
         if(bbnd>1.0 && u22abs<1.0)
           {
            //--- check
            if(bbnd>=bignum*u22abs)
              {
               //--- change values
               scl=1/bbnd;
               br1=scl*br1;
               bi1=scl*bi1;
               br2=scl*br2;
               bi2=scl*bi2;
              }
           }
         //--- function call
         InternalHsEVDLADIV(br2,bi2,ur22,ui22,xr2,xi2);
         xr1=ur11r*br1-ui11r*bi1-ur12s*xr2+ui12s*xi2;
         xi1=ui11r*br1+ur11r*bi1-ui12s*xr2-ur12s*xi2;
         //--- check
         if(zswap4[icmax])
           {
            x[1].Set(1,xr2);
            x[2].Set(1,xr1);
            x[1].Set(2,xi2);
            x[2].Set(2,xi1);
           }
         else
           {
            x[1].Set(1,xr1);
            x[2].Set(1,xr2);
            x[1].Set(2,xi1);
            x[2].Set(2,xi2);
           }
         xnorm=MathMax(MathAbs(xr1)+MathAbs(xi1),MathAbs(xr2)+MathAbs(xi2));
         //--- Further scaling if  norm(A) norm(X) > overflow
         if(xnorm>1.0 && cmax>1.0)
           {
            //--- check
            if(xnorm>bignum/cmax)
              {
               //--- change values
               temp=cmax/bignum;
               x[1].Set(1,temp*x[1][1]);
               x[2].Set(1,temp*x[2][1]);
               x[1].Set(2,temp*x[1][2]);
               x[2].Set(2,temp*x[2][2]);
               xnorm=temp*xnorm;
               scl=temp*scl;
              }
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| performs complex division in  real arithmetic                    |
//|                         a + i*b                                  |
//|              p + i*q = ---------                                 |
//|                         c + i*d                                  |
//| The algorithm is due to Robert L. Smith and can be found         |
//| in D. Knuth, The art of Computer Programming, Vol.2, p.195       |
//|   -- LAPACK auxiliary routine (version 3.0) --                   |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      October 31, 1992                                            |
//+------------------------------------------------------------------+
static void CEigenVDetect::InternalHsEVDLADIV(const double a,const double b,
                                              const double c,const double d,
                                              double &p,double &q)
  {
//--- create variables
   double e=0;
   double f=0;
//--- initialization
   p=0;
   q=0;
//--- check
   if(MathAbs(d)<MathAbs(c))
     {
      //--- get result
      e=d/c;
      f=c+d*e;
      p=(a+b*e)/f;
      q=(b-a*e)/f;
     }
   else
     {
      //--- get result
      e=c/d;
      f=d+c*e;
      p=(b+a*e)/f;
      q=(-a+b*e)/f;
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static bool CEigenVDetect::NonSymmetricEVD(CMatrixDouble &ca,const int n,
                                           const int vneeded,double &wr[],
                                           double &wi[],CMatrixDouble &vl,
                                           CMatrixDouble &vr)
  {
//--- create variables
   bool result;
   int  i=0;
   int  info=0;
   int  m=0;
   int  i_=0;
//--- create arrays
   double tau[];
   bool sel[];
//--- create matrix
   CMatrixDouble s;
//--- create copy
   CMatrixDouble a;
   a=ca;
//--- check
   if(!CAp::Assert(vneeded>=0 && vneeded<=3,__FUNCTION__+": incorrect VNeeded!"))
      return(false);
//--- check
   if(vneeded==0)
     {
      //--- Eigen values only
      ToUpperHessenberg(a,n,tau);
      //--- function call
      CHsSchur::InternalSchurDecomposition(a,n,0,0,wr,wi,s,info);
      //--- get result
      result=info==0;
      //--- return result
      return(result);
     }
//--- Eigen values and vectors
   ToUpperHessenberg(a,n,tau);
//--- function call
   UnpackQFromUpperHessenberg(a,n,tau,s);
//--- function call
   CHsSchur::InternalSchurDecomposition(a,n,1,1,wr,wi,s,info);
//--- get result
   result=info==0;
//--- check
   if(!result)
      return(result);
//--- check
   if(vneeded==1 || vneeded==3)
     {
      vr.Resize(n+1,n+1);
      for(i=1;i<=n;i++)
        {
         for(i_=1;i_<=n;i_++)
            vr[i].Set(i_,s[i][i_]);
        }
     }
//--- check
   if(vneeded==2 || vneeded==3)
     {
      vl.Resize(n+1,n+1);
      for(i=1;i<=n;i++)
        {
         for(i_=1;i_<=n;i_++)
            vl[i].Set(i_,s[i][i_]);
        }
     }
//--- function call
   InternalTREVC(a,n,vneeded,1,sel,vl,vr,m,info);
//--- get result
   result=info==0;
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Upper Hessenberg form                                            |
//+------------------------------------------------------------------+
static void CEigenVDetect::ToUpperHessenberg(CMatrixDouble &a,const int n,double &tau[])
  {
//--- create variables
   int    i=0;
   int    ip1=0;
   int    nmi=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- create arrays
   double t[];
   double work[];
//--- check
   if(!CAp::Assert(n>=0,__FUNCTION__+": incorrect N!"))
      return;
//--- check
   if(n<=1)
      return;
//--- allocation
   ArrayResizeAL(tau,n);
   ArrayResizeAL(t,n+1);
   ArrayResizeAL(work,n+1);
//--- calculations
   for(i=1;i<n;i++)
     {
      //--- Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
      ip1=i+1;
      nmi=n-i;
      i1_=ip1-1;
      for(i_=1;i_<=nmi;i_++)
         t[i_]=a[i_+i1_][i];
      //--- function call
      CReflections::GenerateReflection(t,nmi,v);
      i1_=1-ip1;
      for(i_=ip1;i_<=n;i_++)
         a[i_].Set(i,t[i_+i1_]);
      tau[i]=v;
      t[1]=1;
      //--- Apply H(i) to A(1:ihi,i+1:ihi) from the right
      CReflections::ApplyReflectionFromTheRight(a,v,t,1,n,i+1,n,work);
      //--- Apply H(i) to A(i+1:ihi,i+1:n) from the left
      CReflections::ApplyReflectionFromTheLeft(a,v,t,i+1,n,i+1,n,work);
     }
  }
//+------------------------------------------------------------------+
//| Unpack Q from the matrix of the upper Hessenberg form            |
//+------------------------------------------------------------------+
static void CEigenVDetect::UnpackQFromUpperHessenberg(CMatrixDouble &a,const int n,
                                                      double &tau[],CMatrixDouble &q)
  {
//--- create variables
   int i=0;
   int j=0;
   int ip1=0;
   int nmi=0;
   int i_=0;
   int i1_=0;
//--- create arrays
   double v[];
   double work[];
//--- check
   if(n==0)
      return;
//--- allocation
   q.Resize(n+1,n+1);
   ArrayResizeAL(v,n+1);
   ArrayResizeAL(work,n+1);
//--- identity matrix
   for(i=1;i<=n;i++)
     {
      for(j=1;j<=n;j++)
        {
         //--- check
         if(i==j)
            q[i].Set(j,1);
         else
            q[i].Set(j,0);
        }
     }
//--- unpack Q
   for(i=1;i<n;i++)
     {
      //--- Apply H(i)
      ip1=i+1;
      nmi=n-i;
      i1_=ip1-1;
      for(i_=1;i_<=nmi;i_++)
         v[i_]=a[i_+i1_][i];
      v[1]=1;
      //--- function call
      CReflections::ApplyReflectionFromTheRight(q,tau[i],v,1,n,i+1,n,work);
     }
  }
//+------------------------------------------------------------------+
//| Random matrix generation                                         |
//+------------------------------------------------------------------+
class CMatGen
  {
public:
                     CMatGen(void);
                    ~CMatGen(void);
   //--- public methods
   static void       RMatrixRndOrthogonal(const int n,CMatrixDouble &a);
   static void       RMatrixRndCond(const int n,const double c,CMatrixDouble &a);
   static void       CMatrixRndOrthogonal(const int n,CMatrixComplex &a);
   static void       CMatrixRndCond(const int n,const double c,CMatrixComplex &a);
   static void       SMatrixRndCond(const int n,const double c,CMatrixDouble &a);
   static void       SPDMatrixRndCond(const int n,const double c,CMatrixDouble &a);
   static void       HMatrixRndCond(const int n,const double c,CMatrixComplex &a);
   static void       HPDMatrixRndCond(const int n,const double c,CMatrixComplex &a);
   static void       RMatrixRndOrthogonalFromTheRight(CMatrixDouble &a,const int m,const int n);
   static void       RMatrixRndOrthogonalFromTheLeft(CMatrixDouble &a,const int m,const int n);
   static void       CMatrixRndOrthogonalFromTheRight(CMatrixComplex &a,const int m,const int n);
   static void       CMatrixRndOrthogonalFromTheLeft(CMatrixComplex &a,const int m,const int n);
   static void       SMatrixRndMultiply(CMatrixDouble &a,const int n);
   static void       HMatrixRndMultiply(CMatrixComplex &a,const int n);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMatGen::CMatGen(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMatGen::~CMatGen(void)
  {

  }
//+------------------------------------------------------------------+
//| Generation of a random uniformly distributed (Haar) orthogonal   |
//| matrix                                                           |
//| INPUT PARAMETERS:                                                |
//|     N   -   matrix size, N>=1                                    |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   orthogonal NxN matrix, array[0..N-1,0..N-1]          |
//+------------------------------------------------------------------+
static void CMatGen::RMatrixRndOrthogonal(const int n,CMatrixDouble &a)
  {
//--- create variables
   int i=0;
   int j=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- allocation
   a.Resize(n,n);
   for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
        {
         //--- check
         if(i==j)
            a[i].Set(j,1);
         else
            a[i].Set(j,0);
        }
     }
//--- get result
   RMatrixRndOrthogonalFromTheRight(a,n,n);
  }
//+------------------------------------------------------------------+
//| Generation of random NxN matrix with given condition number and  |
//| norm2(A)=1                                                       |
//| INPUT PARAMETERS:                                                |
//|     N   -   matrix size                                          |
//|     C   -   condition number (in 2-norm)                         |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   random matrix with norm2(A)=1 and cond(A)=C          |
//+------------------------------------------------------------------+
static void CMatGen::RMatrixRndCond(const int n,const double c,CMatrixDouble &a)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double l1=0;
   double l2=0;
//--- check
   if(!CAp::Assert(n>=1&&c>=1.0,__FUNCTION__+": N<1 or C<1!"))
      return;
//--- allocation
   a.Resize(n,n);
//--- check
   if(n==1)
     {
      a[0].Set(0,2*CMath::RandomInteger(2)-1);
      //--- exit the function
      return;
     }
//--- initialization
   l1=0;
   l2=MathLog(1/c);
   for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
         a[i].Set(j,0);
     }
//--- change a
   a[0].Set(0,MathExp(l1));
   for(i=1;i<=n-2;i++)
      a[i].Set(i,MathExp(CMath::RandomReal()*(l2-l1)+l1));
   a[n-1].Set(n-1,MathExp(l2));
//--- function call
   RMatrixRndOrthogonalFromTheLeft(a,n,n);
//--- function call
   RMatrixRndOrthogonalFromTheRight(a,n,n);
  }
//+------------------------------------------------------------------+
//| Generation of a random Haar distributed orthogonal complex matrix|
//| INPUT PARAMETERS:                                                |
//|     N   -   matrix size, N>=1                                    |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   orthogonal NxN matrix, array[0..N-1,0..N-1]          |
//+------------------------------------------------------------------+
static void CMatGen::CMatrixRndOrthogonal(const int n,CMatrixComplex &a)
  {
//--- create variables
   int i=0;
   int j=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return;
//--- allocation
   a.Resize(n,n);
   for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
        {
         //--- check
         if(i==j)
            a[i].Set(j,1);
         else
            a[i].Set(j,0);
        }
     }
//--- get result
   CMatrixRndOrthogonalFromTheRight(a,n,n);
  }
//+------------------------------------------------------------------+
//| Generation of random NxN complex matrix with given condition     |
//| number C and norm2(A)=1                                          |
//| INPUT PARAMETERS:                                                |
//|     N   -   matrix size                                          |
//|     C   -   condition number (in 2-norm)                         |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   random matrix with norm2(A)=1 and cond(A)=C          |
//+------------------------------------------------------------------+
static void CMatGen::CMatrixRndCond(const int n,const double c,CMatrixComplex &a)
  {
//--- create variables
   int        i=0;
   int        j=0;
   double     l1=0;
   double     l2=0;
   al_complex v=0;
//--- object of class
   CHighQualityRandState state;
//--- check
   if(!CAp::Assert(n>=1&&c>=1.0,__FUNCTION__+": N<1 or C<1!"))
      return;
//--- allocation
   a.Resize(n,n);
//--- check
   if(n==1)
     {
      //--- function call
      CHighQualityRand::HQRndRandomize(state);
      //--- function call
      CHighQualityRand::HQRndUnit2(state,v.re,v.im);
      a[0].Set(0,v);
      return;
     }
//--- initialization
   l1=0;
   l2=MathLog(1/c);
   for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
         a[i].Set(j,0);
     }
//--- change values
   a[0].Set(0,MathExp(l1));
   for(i=1;i<=n-2;i++)
      a[i].Set(i,MathExp(CMath::RandomReal()*(l2-l1)+l1));
   a[n-1].Set(n-1,MathExp(l2));
//--- function call
   CMatrixRndOrthogonalFromTheLeft(a,n,n);
//--- function call
   CMatrixRndOrthogonalFromTheRight(a,n,n);
  }
//+------------------------------------------------------------------+
//| Generation of random NxN symmetric matrix with given condition   |
//| number and norm2(A)=1                                            |
//| INPUT PARAMETERS:                                                |
//|     N   -   matrix size                                          |
//|     C   -   condition number (in 2-norm)                         |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   random matrix with norm2(A)=1 and cond(A)=C          |
//+------------------------------------------------------------------+
static void CMatGen::SMatrixRndCond(const int n,const double c,CMatrixDouble &a)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double l1=0;
   double l2=0;
//--- check
   if(!CAp::Assert(n>=1&&c>=1.0,__FUNCTION__+": N<1 or C<1!"))
      return;
//--- allocation
   a.Resize(n,n);
//--- check
   if(n==1)
     {
      a[0].Set(0,2*CMath::RandomInteger(2)-1);
      return;
     }
//--- Prepare matrix
   l1=0;
   l2=MathLog(1/c);
   for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
         a[i].Set(j,0);
     }
//--- initialization
   a[0].Set(0,MathExp(l1));
   for(i=1;i<=n-2;i++)
      a[i].Set(i,(2*CMath::RandomInteger(2)-1)*MathExp(CMath::RandomReal()*(l2-l1)+l1));
   a[n-1].Set(n-1,MathExp(l2));
//--- Multiply
   SMatrixRndMultiply(a,n);
  }
//+------------------------------------------------------------------+
//| Generation of random NxN symmetric positive definite matrix with |
//| given condition number and norm2(A)=1                            |
//| INPUT PARAMETERS:                                                |
//|     N   -   matrix size                                          |
//|     C   -   condition number (in 2-norm)                         |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   random SPD matrix with norm2(A)=1 and cond(A)=C      |
//+------------------------------------------------------------------+
static void CMatGen::SPDMatrixRndCond(const int n,const double c,CMatrixDouble &a)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double l1=0;
   double l2=0;
//--- check
   if(n<=0 || c<1.0)
      return;
//--- allocation
   a.Resize(n,n);
//--- check
   if(n==1)
     {
      a[0].Set(0,1);
      return;
     }
//--- Prepare matrix
   l1=0;
   l2=MathLog(1/c);
   for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
         a[i].Set(j,0);
     }
//--- initialization
   a[0].Set(0,MathExp(l1));
   for(i=1;i<=n-2;i++)
      a[i].Set(i,MathExp(CMath::RandomReal()*(l2-l1)+l1));
   a[n-1].Set(n-1,MathExp(l2));
//--- Multiply
   SMatrixRndMultiply(a,n);
  }
//+------------------------------------------------------------------+
//| Generation of random NxN Hermitian matrix with given condition   |
//| number and norm2(A)=1                                            |
//| INPUT PARAMETERS:                                                |
//|     N   -   matrix size                                          |
//|     C   -   condition number (in 2-norm)                         |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   random matrix with norm2(A)=1 and cond(A)=C          |
//+------------------------------------------------------------------+
static void CMatGen::HMatrixRndCond(const int n,const double c,CMatrixComplex &a)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double l1=0;
   double l2=0;
//--- check
   if(!CAp::Assert(n>=1&&c>=1.0,__FUNCTION__+": N<1 or C<1!"))
      return;
//--- allocation
   a.Resize(n,n);
//--- check
   if(n==1)
     {
      a[0].Set(0,2*CMath::RandomInteger(2)-1);
      return;
     }
//--- Prepare matrix
   l1=0;
   l2=MathLog(1/c);
   for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
         a[i].Set(j,0);
     }
//--- initialization
   a[0].Set(0,MathExp(l1));
   for(i=1;i<=n-2;i++)
      a[i].Set(i,(2*CMath::RandomInteger(2)-1)*MathExp(CMath::RandomReal()*(l2-l1)+l1));
   a[n-1].Set(n-1,MathExp(l2));
//--- Multiply
   HMatrixRndMultiply(a,n);
//--- post-process to ensure that matrix diagonal is real
   for(i=0;i<n;i++)
      a[i].SetIm(i,0);
  }
//+------------------------------------------------------------------+
//| Generation of random NxN Hermitian positive definite matrix with |
//| given condition number and norm2(A)=1                            |
//| INPUT PARAMETERS:                                                |
//|     N   -   matrix size                                          |
//|     C   -   condition number (in 2-norm)                         |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   random HPD matrix with norm2(A)=1 and cond(A)=C      |
//+------------------------------------------------------------------+
static void CMatGen::HPDMatrixRndCond(const int n,const double c,CMatrixComplex &a)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double l1=0;
   double l2=0;
//--- check
   if(n<=0 || c<1.0)
      return;
//--- allocation
   a.Resize(n,n);
//--- check
   if(n==1)
     {
      a[0].Set(0,1);
      return;
     }
//--- Prepare matrix
   l1=0;
   l2=MathLog(1/c);
   for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
         a[i].Set(j,0);
     }
//--- initialization
   a[0].Set(0,MathExp(l1));
   for(i=1;i<=n-2;i++)
      a[i].Set(i,MathExp(CMath::RandomReal()*(l2-l1)+l1));
   a[n-1].Set(n-1,MathExp(l2));
//--- Multiply
   HMatrixRndMultiply(a,n);
//--- post-process to ensure that matrix diagonal is real
   for(i=0;i<n;i++)
      a[i].SetIm(i,0);
  }
//+------------------------------------------------------------------+
//| Multiplication of MxN matrix by NxN random Haar distributed      |
//| orthogonal matrix                                                |
//| INPUT PARAMETERS:                                                |
//|     A   -   matrix, array[0..M-1, 0..N-1]                        |
//|     M, N-   matrix size                                          |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   A*Q, where Q is random NxN orthogonal matrix         |
//+------------------------------------------------------------------+
static void CMatGen::RMatrixRndOrthogonalFromTheRight(CMatrixDouble &a,const int m,
                                                      const int n)
  {
//--- create variables
   double tau=0;
   double lambdav=0;
   int    s=0;
   int    i=0;
   double u1=0;
   double u2=0;
   int i_=0;
//--- create arrays
   double w[];
   double v[];
//--- object of class
   CHighQualityRandState state;
//--- check
   if(!CAp::Assert(n>=1 && m>=1,__FUNCTION__+": N<1 or M<1!"))
      return;
//--- check
   if(n==1)
     {
      tau=2*CMath::RandomInteger(2)-1;
      for(i=0;i<m;i++)
         a[i].Set(0,a[i][0]*tau);
      //--- exit the function
      return;
     }
//--- General case.
//--- First pass.
   ArrayResizeAL(w,m);
   ArrayResizeAL(v,n+1);
//--- function call
   CHighQualityRand::HQRndRandomize(state);
   for(s=2;s<=n;s++)
     {
      //--- Prepare random normal v
      do
        {
         i=1;
         while(i<=s)
           {
            //--- function call
            CHighQualityRand::HQRndNormal2(state,u1,u2);
            v[i]=u1;
            //--- check
            if(i+1<=s)
               v[i+1]=u2;
            i=i+2;
           }
         //--- change values
         lambdav=0.0;
         for(i_=1;i_<=s;i_++)
            lambdav+=v[i_]*v[i_];
        }
      while(lambdav==0.0);
      //--- Prepare and apply reflection
      CReflections::GenerateReflection(v,s,tau);
      v[1]=1;
      //--- function call
      CReflections::ApplyReflectionFromTheRight(a,tau,v,0,m-1,n-s,n-1,w);
     }
//--- Second pass.
   for(i=0;i<n;i++)
     {
      tau=2*CMath::RandomInteger(2)-1;
      for(i_=0;i_<m;i_++)
         a[i_].Set(i,tau*a[i_][i]);
     }
  }
//+------------------------------------------------------------------+
//| Multiplication of MxN matrix by MxM random Haar distributed      |
//| orthogonal matrix                                                |
//| INPUT PARAMETERS:                                                |
//|     A   -   matrix, array[0..M-1, 0..N-1]                        |
//|     M, N-   matrix size                                          |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   Q*A, where Q is random MxM orthogonal matrix         |
//+------------------------------------------------------------------+
static void CMatGen::RMatrixRndOrthogonalFromTheLeft(CMatrixDouble &a,const int m,
                                                     const int n)
  {
//--- create variables
   double tau=0;
   double lambdav=0;
   int    s=0;
   int    i=0;
   int    j=0;
   double u1=0;
   double u2=0;
//--- create arrays
   double w[];
   double v[];
//--- object of class
   CHighQualityRandState state;
   int i_=0;
//--- check
   if(!CAp::Assert(n>=1 && m>=1,__FUNCTION__+": N<1 or M<1!"))
      return;
//--- check
   if(m==1)
     {
      tau=2*CMath::RandomInteger(2)-1;
      for(j=0;j<n;j++)
         a[0].Set(j,a[0][j]*tau);
      //--- exit the function
      return;
     }
//--- General case.
//--- First pass.
   ArrayResizeAL(w,n);
   ArrayResizeAL(v,m+1);
//--- function call
   CHighQualityRand::HQRndRandomize(state);
   for(s=2;s<=m;s++)
     {
      //--- Prepare random normal v
      do
        {
         i=1;
         while(i<=s)
           {
            //--- function call
            CHighQualityRand::HQRndNormal2(state,u1,u2);
            v[i]=u1;
            //--- check
            if(i+1<=s)
               v[i+1]=u2;
            i=i+2;
           }
         //--- change values
         lambdav=0.0;
         for(i_=1;i_<=s;i_++)
            lambdav+=v[i_]*v[i_];
        }
      while(lambdav==0.0);
      //--- Prepare random normal v
      CReflections::GenerateReflection(v,s,tau);
      v[1]=1;
      //--- function call
      CReflections::ApplyReflectionFromTheLeft(a,tau,v,m-s,m-1,0,n-1,w);
     }
//--- Second pass.
   for(i=0;i<m;i++)
     {
      tau=2*CMath::RandomInteger(2)-1;
      for(i_=0;i_<n;i_++)
         a[i].Set(i_,tau*a[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| Multiplication of MxN complex matrix by NxN random Haar          |
//| distributed complex orthogonal matrix                            |
//| INPUT PARAMETERS:                                                |
//|     A   -   matrix, array[0..M-1, 0..N-1]                        |
//|     M, N-   matrix size                                          |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   A*Q, where Q is random NxN orthogonal matrix         |
//+------------------------------------------------------------------+
static void CMatGen::CMatrixRndOrthogonalFromTheRight(CMatrixComplex &a,const int m,const int n)
  {
//--- create variables
   al_complex zero=0;
   al_complex lambdav=0;
   al_complex tau=0;
   int        s=0;
   int        i=0;
   int        i_=0;
//--- create arrays
   al_complex w[];
   al_complex v[];
//--- object of class
   CHighQualityRandState state;
//--- check
   if(!CAp::Assert(n>=1 && m>=1,__FUNCTION__+": N<1 or M<1!"))
      return;
//--- check
   if(n==1)
     {
      //--- function call
      CHighQualityRand::HQRndRandomize(state);
      //--- function call
      CHighQualityRand::HQRndUnit2(state,tau.re,tau.im);
      for(i=0;i<m;i++)
         a[i].Set(0,a[i][0]*tau);
      //--- exit the function
      return;
     }
//--- General case.
//--- First pass.
   ArrayResizeAL(w,m);
   ArrayResizeAL(v,n+1);
//--- function call
   CHighQualityRand::HQRndRandomize(state);
   for(s=2;s<=n;s++)
     {
      //--- Prepare random normal v
      do
        {
         for(i=1;i<=s;i++)
           {
            //--- function call
            CHighQualityRand::HQRndNormal2(state,tau.re,tau.im);
            v[i]=tau;
           }
         //--- change values
         lambdav=0.0;
         for(i_=1;i_<=s;i_++)
            lambdav+=v[i_]*CMath::Conj(v[i_]);
        }
      while(lambdav==zero);
      //--- Prepare and apply reflection
      CComplexReflections::ComplexGenerateReflection(v,s,tau);
      v[1]=1;
      //--- function call
      CComplexReflections::ComplexApplyReflectionFromTheRight(a,tau,v,0,m-1,n-s,n-1,w);
     }
//--- Second pass.
   for(i=0;i<n;i++)
     {
      //--- function call
      CHighQualityRand::HQRndUnit2(state,tau.re,tau.im);
      for(i_=0;i_<m;i_++)
         a[i_].Set(i,tau*a[i_][i]);
     }
  }
//+------------------------------------------------------------------+
//| Multiplication of MxN complex matrix by MxM random Haar          |
//| distributed complex orthogonal matrix                            |
//| INPUT PARAMETERS:                                                |
//|     A   -   matrix, array[0..M-1, 0..N-1]                        |
//|     M, N-   matrix size                                          |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   Q*A, where Q is random MxM orthogonal matrix         |
//+------------------------------------------------------------------+
static void CMatGen::CMatrixRndOrthogonalFromTheLeft(CMatrixComplex &a,const int m,
                                                     const int n)
  {
//--- create variables
   al_complex zero=0;
   al_complex tau=0;
   al_complex lambdav=0;
   int        s=0;
   int        i=0;
   int        j=0;
   int        i_=0;
//--- create arrays
   al_complex w[];
   al_complex v[];
//--- object of class
   CHighQualityRandState state;
//--- check
   if(!CAp::Assert(n>=1 && m>=1,__FUNCTION__+": N<1 or M<1!"))
      return;
//--- check
   if(m==1)
     {
      //--- function call
      CHighQualityRand::HQRndRandomize(state);
      //--- function call
      CHighQualityRand::HQRndUnit2(state,tau.re,tau.im);
      for(j=0;j<n;j++)
         a[0].Set(j,a[0][j]*tau);
      //--- exit the function
      return;
     }
//--- General case.
//--- First pass.
   ArrayResizeAL(w,n);
   ArrayResizeAL(v,m+1);
//--- function call
   CHighQualityRand::HQRndRandomize(state);
   for(s=2;s<=m;s++)
     {
      //--- Prepare random normal v
      do
        {
         for(i=1;i<=s;i++)
           {
            //--- function call
            CHighQualityRand::HQRndNormal2(state,tau.re,tau.im);
            v[i]=tau;
           }
         //--- change values
         lambdav=0.0;
         for(i_=1;i_<=s;i_++)
            lambdav+=v[i_]*CMath::Conj(v[i_]);
        }
      while(lambdav==zero);
      //--- Prepare and apply reflection
      CComplexReflections::ComplexGenerateReflection(v,s,tau);
      v[1]=1;
      //--- function call
      CComplexReflections::ComplexApplyReflectionFromTheLeft(a,tau,v,m-s,m-1,0,n-1,w);
     }
//--- Second pass.
   for(i=0;i<m;i++)
     {
      //--- function call
      CHighQualityRand::HQRndUnit2(state,tau.re,tau.im);
      for(i_=0;i_<n;i_++)
         a[i].Set(i_,tau*a[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| Symmetric multiplication of NxN matrix by random Haar            |
//| distributed orthogonal matrix                                    |
//| INPUT PARAMETERS:                                                |
//|     A   -   matrix, array[0..N-1, 0..N-1]                        |
//|     N   -   matrix size                                          |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   Q'*A*Q, where Q is random NxN orthogonal matrix      |
//+------------------------------------------------------------------+
static void CMatGen::SMatrixRndMultiply(CMatrixDouble &a,const int n)
  {
//--- create variables
   double tau=0;
   double lambdav=0;
   int    s=0;
   int    i=0;
   double u1=0;
   double u2=0;
   int    i_=0;
//--- create arrays
   double w[];
   double v[];
//--- object of class
   CHighQualityRandState state;
//--- General case.
   ArrayResizeAL(w,n);
   ArrayResizeAL(v,n+1);
//--- function call
   CHighQualityRand::HQRndRandomize(state);
   for(s=2;s<=n;s++)
     {
      //--- Prepare random normal v
      do
        {
         i=1;
         while(i<=s)
           {
            //--- function call
            CHighQualityRand::HQRndNormal2(state,u1,u2);
            v[i]=u1;
            //--- check
            if(i+1<=s)
               v[i+1]=u2;
            i=i+2;
           }
         //--- change values
         lambdav=0.0;
         for(i_=1;i_<=s;i_++)
            lambdav+=v[i_]*v[i_];
        }
      while(lambdav==0.0);
      //--- Prepare and apply reflection
      CReflections::GenerateReflection(v,s,tau);
      v[1]=1;
      //--- function call
      CReflections::ApplyReflectionFromTheRight(a,tau,v,0,n-1,n-s,n-1,w);
      //--- function call
      CReflections::ApplyReflectionFromTheLeft(a,tau,v,n-s,n-1,0,n-1,w);
     }
//--- Second pass.
   for(i=0;i<n;i++)
     {
      tau=2*CMath::RandomInteger(2)-1;
      for(i_=0;i_<n;i_++)
         a[i_].Set(i,tau*a[i_][i]);
      for(i_=0;i_<n;i_++)
         a[i].Set(i_,tau*a[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| Hermitian multiplication of NxN matrix by random Haar distributed|
//| complex orthogonal matrix                                        |
//| INPUT PARAMETERS:                                                |
//|     A   -   matrix, array[0..N-1, 0..N-1]                        |
//|     N   -   matrix size                                          |
//| OUTPUT PARAMETERS:                                               |
//|     A   -   Q^H*A*Q, where Q is random NxN orthogonal matrix     |
//+------------------------------------------------------------------+
static void CMatGen::HMatrixRndMultiply(CMatrixComplex &a,const int n)
  {
//--- create variables
   al_complex zero=0;
   al_complex tau=0;
   al_complex lambdav=0;
   int        s=0;
   int        i=0;
   int        i_=0;
//--- create arrays
   al_complex w[];
   al_complex v[];
//--- object of class
   CHighQualityRandState state;
//--- General case.
   ArrayResizeAL(w,n);
   ArrayResizeAL(v,n+1);
//--- function call
   CHighQualityRand::HQRndRandomize(state);
   for(s=2;s<=n;s++)
     {
      //--- Prepare random normal v
      do
        {
         for(i=1;i<=s;i++)
           {
            //--- function call
            CHighQualityRand::HQRndNormal2(state,tau.re,tau.im);
            v[i]=tau;
           }
         //--- change values
         lambdav=0.0;
         for(i_=1;i_<=s;i_++)
            lambdav+=v[i_]*CMath::Conj(v[i_]);
        }
      while(lambdav==zero);
      //--- Prepare and apply reflection
      CComplexReflections::ComplexGenerateReflection(v,s,tau);
      v[1]=1;
      //--- function call
      CComplexReflections::ComplexApplyReflectionFromTheRight(a,tau,v,0,n-1,n-s,n-1,w);
      //--- function call
      CComplexReflections::ComplexApplyReflectionFromTheLeft(a,CMath::Conj(tau),v,n-s,n-1,0,n-1,w);
     }
//--- Second pass.
   for(i=0;i<n;i++)
     {
      //--- function call
      CHighQualityRand::HQRndUnit2(state,tau.re,tau.im);
      for(i_=0;i_<n;i_++)
         a[i_].Set(i,tau*a[i_][i]);
      tau=CMath::Conj(tau);
      for(i_=0;i_<n;i_++)
         a[i].Set(i_,tau*a[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| Triangular factorizations                                        |
//+------------------------------------------------------------------+
class CTrFac
  {
private:
   //--- private methods
   static void       CMatrixLUPRec(CMatrixComplex &a,const int offs,const int m,const int n,int &pivots[],al_complex &tmp[]);
   static void       RMatrixLUPRec(CMatrixDouble &a,const int offs,const int m,const int n,int &pivots[],double &tmp[]);
   static void       CMatrixPLURec(CMatrixComplex &a,const int offs,const int m,const int n,int &pivots[],al_complex &tmp[]);
   static void       RMatrixPLURec(CMatrixDouble &a,const int offs,const int m,const int n,int &pivots[],double &tmp[]);
   static void       CMatrixLUP2(CMatrixComplex &a,const int offs,const int m,const int n,int &pivots[],al_complex &tmp[]);
   static void       RMatrixLUP2(CMatrixDouble &a,const int offs,const int m,const int n,int &pivots[],double &tmp[]);
   static void       CMatrixPLU2(CMatrixComplex &a,const int offs,const int m,const int n,int &pivots[],al_complex &tmp[]);
   static void       RMatrixPLU2(CMatrixDouble &a,const int offs,const int m,const int n,int &pivots[],double &tmp[]);
   static bool       HPDMatrixCholeskyRec(CMatrixComplex &a,const int offs,const int n,const bool isupper,al_complex &tmp[]);
   static bool       HPDMatrixCholesky2(CMatrixComplex &aaa,const int offs,const int n,const bool isupper,al_complex &tmp[]);
   static bool       SPDMatrixCholesky2(CMatrixDouble &aaa,const int offs,const int n,const bool isupper,double &tmp[]);
public:
                     CTrFac(void);
                    ~CTrFac(void);
   //--- public methods
   static void       RMatrixLU(CMatrixDouble &a,const int m,const int n,int &pivots[]);
   static void       CMatrixLU(CMatrixComplex &a,const int m,const int n,int &pivots[]);
   static bool       HPDMatrixCholesky(CMatrixComplex &a,const int n,const bool isupper);
   static bool       SPDMatrixCholesky(CMatrixDouble &a,const int n,const bool isupper);
   static void       RMatrixLUP(CMatrixDouble &a,const int m,const int n,int &pivots[]);
   static void       CMatrixLUP(CMatrixComplex &a,const int m,const int n,int &pivots[]);
   static void       RMatrixPLU(CMatrixDouble &a,const int m,const int n,int &pivots[]);
   static void       CMatrixPLU(CMatrixComplex &a,const int m,const int n,int &pivots[]);
   static bool       SPDMatrixCholeskyRec(CMatrixDouble &a,const int offs,const int n,const bool isupper,double &tmp[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CTrFac::CTrFac(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CTrFac::~CTrFac(void)
  {

  }
//+------------------------------------------------------------------+
//| LU decomposition of a general real matrix with row pivoting      |
//| A is represented as A = P*L*U, where:                            |
//| * L is lower unitriangular matrix                                |
//| * U is upper triangular matrix                                   |
//| * P = P0*P1*...*PK, K=min(M,N)-1,                                |
//|   Pi - permutation matrix for I and Pivots[I]                    |
//| This is cache-oblivous implementation of LU decomposition.       |
//| It is optimized for square matrices. As for rectangular matrices:|
//| * best case - M>>N                                               |
//| * worst case - N>>M, small M, large N, matrix does not fit in CPU|
//|   cache                                                          |
//| INPUT PARAMETERS:                                                |
//|     A       -   array[0..M-1, 0..N-1].                           |
//|     M       -   number of rows in matrix A.                      |
//|     N       -   number of columns in matrix A.                   |
//| OUTPUT PARAMETERS:                                               |
//|     A       -   matrices L and U in compact form:                |
//|                 * L is stored under main diagonal                |
//|                 * U is stored on and above main diagonal         |
//|     Pivots  -   permutation matrix in compact form.              |
//|                 array[0..Min(M-1,N-1)].                          |
//+------------------------------------------------------------------+
static void CTrFac::RMatrixLU(CMatrixDouble &a,const int m,const int n,int &pivots[])
  {
//--- check
   if(!CAp::Assert(m>0,__FUNCTION__+": incorrect M!"))
      return;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- function call
   RMatrixPLU(a,m,n,pivots);
  }
//+------------------------------------------------------------------+
//| LU decomposition of a general complex matrix with row pivoting   |
//| A is represented as A = P*L*U, where:                            |
//| * L is lower unitriangular matrix                                |
//| * U is upper triangular matrix                                   |
//| * P = P0*P1*...*PK, K=min(M,N)-1,                                |
//|   Pi - permutation matrix for I and Pivots[I]                    |
//| This is cache-oblivous implementation of LU decomposition. It is |
//| optimized for square matrices. As for rectangular matrices:      |
//| * best case - M>>N                                               |
//| * worst case - N>>M, small M, large N, matrix does not fit in CPU|
//| cache                                                            |
//| INPUT PARAMETERS:                                                |
//|     A       -   array[0..M-1, 0..N-1].                           |
//|     M       -   number of rows in matrix A.                      |
//|     N       -   number of columns in matrix A.                   |
//| OUTPUT PARAMETERS:                                               |
//|     A       -   matrices L and U in compact form:                |
//|                 * L is stored under main diagonal                |
//|                 * U is stored on and above main diagonal         |
//|     Pivots  -   permutation matrix in compact form.              |
//|                 array[0..Min(M-1,N-1)].                          |
//+------------------------------------------------------------------+
static void CTrFac::CMatrixLU(CMatrixComplex &a,const int m,const int n,int &pivots[])
  {
//--- check
   if(!CAp::Assert(m>0,__FUNCTION__+": incorrect M!"))
      return;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- function call
   CMatrixPLU(a,m,n,pivots);
  }
//+------------------------------------------------------------------+
//| Cache-oblivious Cholesky decomposition                           |
//| The algorithm computes Cholesky decomposition of a Hermitian     |
//| positive - definite matrix. The result of an algorithm is a      |
//| representation of A as A=U'*U or A=L*L' (here X' detones         |
//| conj(X^T)).                                                      |
//| INPUT PARAMETERS:                                                |
//|     A       -   upper or lower triangle of a factorized matrix.  |
//|                 array with elements [0..N-1, 0..N-1].            |
//|     N       -   size of matrix A.                                |
//|     IsUpper -   if IsUpper=True, then A contains an upper        |
//|                 triangle of a symmetric matrix, otherwise A      |
//|                 contains a lower one.                            |
//| OUTPUT PARAMETERS:                                               |
//|     A       -   the result of factorization. If IsUpper=True,    |
//|                 then the upper triangle contains matrix U, so    |
//|                 that A = U'*U, and the elements below the main   |
//|                 diagonal are not modified. Similarly, if         |
//|                 IsUpper = False.                                 |
//| RESULT:                                                          |
//|     If the matrix is positive-definite, the function returns     |
//|     True. Otherwise, the function returns False. Contents of A is|
//|     not determined in such case.                                 |
//+------------------------------------------------------------------+
static bool CTrFac::HPDMatrixCholesky(CMatrixComplex &a,const int n,const bool isupper)
  {
//--- create array
   al_complex tmp[];
//--- check
   if(n<1)
      return(false);
//--- return result
   return(HPDMatrixCholeskyRec(a,0,n,isupper,tmp));
  }
//+------------------------------------------------------------------+
//| Cache-oblivious Cholesky decomposition                           |
//| The algorithm computes Cholesky decomposition of a symmetric     |
//| positive - definite matrix. The result of an algorithm is a      |
//| representation of A as A=U^T*U  or A=L*L^T                       |
//| INPUT PARAMETERS:                                                |
//|     A       -   upper or lower triangle of a factorized matrix.  |
//|                 array with elements [0..N-1, 0..N-1].            |
//|     N       -   size of matrix A.                                |
//|     IsUpper -   if IsUpper=True, then A contains an upper        |
//|                 triangle of a symmetric matrix, otherwise A      |
//|                 contains a lower one.                            |
//| OUTPUT PARAMETERS:                                               |
//|     A       -   the result of factorization. If IsUpper=True,    |
//|                 then the upper triangle contains matrix U, so    |
//|                 that A = U^T*U, and the elements below the main  |
//|                 diagonal are not modified. Similarly, if         |
//|                 IsUpper = False.                                 |
//| RESULT:                                                          |
//|     If the matrix is positive-definite, the function returns     |
//|     True. Otherwise, the function returns False. Contents of A is|
//|     not determined in such case.                                 |
//+------------------------------------------------------------------+
static bool CTrFac::SPDMatrixCholesky(CMatrixDouble &a,const int n,const bool isupper)
  {
//--- create array
   double tmp[];
//--- check
   if(n<1)
      return(false);
//--- return result
   return(SPDMatrixCholeskyRec(a,0,n,isupper,tmp));
  }
//+------------------------------------------------------------------+
//| LUP decomposition of general real matrix                         |
//+------------------------------------------------------------------+
static void CTrFac::RMatrixLUP(CMatrixDouble &a,const int m,const int n,int &pivots[])
  {
//--- create variables
   int    i=0;
   int    j=0;
   double mx=0;
   double v=0;
   int    i_=0;
//--- create array
   double tmp[];
//--- Internal LU decomposition subroutine.
//--- Never call it directly.
   if(!CAp::Assert(m>0,__FUNCTION__+": incorrect M!"))
      return;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- Scale matrix to avoid overflows,
//--- decompose it, then scale back.
   mx=0;
   for(i=0;i<m;i++)
     {
      for(j=0;j<n;j++)
         mx=MathMax(mx,MathAbs(a[i][j]));
     }
//--- check
   if(mx!=0.0)
     {
      v=1/mx;
      //--- change matrix
      for(i=0;i<m;i++)
        {
         for(i_=0;i_<n;i_++)
            a[i].Set(i_,v*a[i][i_]);
        }
     }
//--- allocation
   ArrayResizeAL(pivots,MathMin(m,n));
   ArrayResizeAL(tmp,2*MathMax(m,n));
//--- function call
   RMatrixLUPRec(a,0,m,n,pivots,tmp);
//--- check
   if(mx!=0.0)
     {
      v=mx;
      //--- get result
      for(i=0;i<m;i++)
         for(i_=0;i_<=MathMin(i,n-1);i_++)
            a[i].Set(i_,v*a[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| LUP decomposition of general complex matrix                      |
//+------------------------------------------------------------------+
static void CTrFac::CMatrixLUP(CMatrixComplex &a,const int m,const int n,int &pivots[])
  {
//--- create variables
   int        i=0;
   int        j=0;
   double     mx=0;
   double     v=0;
   al_complex cV;
   int        i_=0;
//--- create array
   al_complex tmp[];
//--- Internal LU decomposition subroutine.
//--- Never call it directly.
   if(!CAp::Assert(m>0,__FUNCTION__+": incorrect M!"))
      return;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- Scale matrix to avoid overflows,
//--- decompose it, then scale back.
   mx=0;
   for(i=0;i<m;i++)
     {
      for(j=0;j<n;j++)
         mx=MathMax(mx,CMath::AbsComplex(a[i][j]));
     }
//--- check
   if(mx!=0.0)
     {
      v=1/mx;
      cV=v;
      //--- change matrix
      for(i=0;i<m;i++)
        {
         for(i_=0;i_<n;i_++)
            a[i].Set(i_,cV*a[i][i_]);
        }
     }
//--- allocation
   ArrayResizeAL(pivots,MathMin(m,n));
   ArrayResizeAL(tmp,2*MathMax(m,n));
//--- function call
   CMatrixLUPRec(a,0,m,n,pivots,tmp);
//--- check
   if(mx!=0.0)
     {
      v=mx;
      cV=v;
      //--- get result
      for(i=0;i<m;i++)
        {
         for(i_=0;i_<=MathMin(i,n-1);i_++)
            a[i].Set(i_,cV*a[i][i_]);
        }
     }
  }
//+------------------------------------------------------------------+
//| PLU decomposition of general real matrix                         |
//+------------------------------------------------------------------+
static void CTrFac::RMatrixPLU(CMatrixDouble &a,const int m,const int n,int &pivots[])
  {
//--- create variables
   int    i=0;
   int    j=0;
   double mx=0;
   double v=0;
   int    i_=0;
//--- create array
   double tmp[];
//--- Internal LU decomposition subroutine.
//--- Never call it directly.
   if(!CAp::Assert(m>0,__FUNCTION__+": incorrect M!"))
      return;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- allocation
   ArrayResizeAL(tmp,2*MathMax(m,n));
   ArrayResizeAL(pivots,MathMin(m,n));
//--- Scale matrix to avoid overflows,
//--- decompose it, then scale back.
   mx=0;
   for(i=0;i<m;i++)
     {
      for(j=0;j<n;j++)
         mx=MathMax(mx,MathAbs(a[i][j]));
     }
//--- check
   if(mx!=0.0)
     {
      v=1/mx;
      //--- change matrix
      for(i=0;i<m;i++)
        {
         for(i_=0;i_<n;i_++)
            a[i].Set(i_,v*a[i][i_]);
        }
     }
//--- function call
   RMatrixPLURec(a,0,m,n,pivots,tmp);
//--- check
   if(mx!=0.0)
     {
      v=mx;
      //--- get result
      for(i=0;i<=MathMin(m,n)-1;i++)
        {
         for(i_=i;i_<n;i_++)
            a[i].Set(i_,v*a[i][i_]);
        }
     }
  }
//+------------------------------------------------------------------+
//| PLU decomposition of general complex matrix                      |
//+------------------------------------------------------------------+
static void CTrFac::CMatrixPLU(CMatrixComplex &a,const int m,const int n,int &pivots[])
  {
//--- create variables
   int        i=0;
   int        j=0;
   double     mx=0;
   al_complex v=0;
   int        i_=0;
//--- create array
   al_complex tmp[];
//--- Internal LU decomposition subroutine.
//--- Never call it directly.
   if(!CAp::Assert(m>0,__FUNCTION__+": incorrect M!"))
      return;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- allocation
   ArrayResizeAL(tmp,2*MathMax(m,n));
   ArrayResizeAL(pivots,MathMin(m,n));
//--- Scale matrix to avoid overflows,
//--- decompose it, then scale back.
   mx=0;
   for(i=0;i<m;i++)
     {
      for(j=0;j<n;j++)
         mx=MathMax(mx,CMath::AbsComplex(a[i][j]));
     }
//--- check
   if(mx!=0.0)
     {
      v=1/mx;
      //--- change matrix
      for(i=0;i<m;i++)
        {
         for(i_=0;i_<n;i_++)
            a[i].Set(i_,v*a[i][i_]);
        }
     }
//--- function call
   CMatrixPLURec(a,0,m,n,pivots,tmp);
//--- check
   if(mx!=0.0)
     {
      v=mx;
      //--- get result
      for(i=0;i<=MathMin(m,n)-1;i++)
        {
         for(i_=i;i_<n;i_++)
            a[i].Set(i_,v*a[i][i_]);
        }
     }
  }
//+------------------------------------------------------------------+
//| Recursive computational subroutine for SPDMatrixCholesky.        |
//| INPUT PARAMETERS:                                                |
//|     A       -   matrix given by upper or lower triangle          |
//|     Offs    -   offset of diagonal block to decompose            |
//|     N       -   diagonal block size                              |
//|     IsUpper -   what half is given                               |
//|     Tmp     -   temporary array; allocated by function, if its   |
//|                 size is too small; can be reused on subsequent   | 
//|                 calls.                                           |
//| OUTPUT PARAMETERS:                                               |
//|     A       -   upper (or lower) triangle contains Cholesky      |
//|                 decomposition                                    |
//| RESULT:                                                          |
//|     True, on success                                             |
//|     False, on failure                                            |
//+------------------------------------------------------------------+
static bool CTrFac::SPDMatrixCholeskyRec(CMatrixDouble &a,const int offs,const int n,
                                         const bool isupper,double &tmp[])
  {
//--- create variables
   bool result;
   int  n1=0;
   int  n2=0;
//--- check
   if(n<1)
      return(false);
//--- prepare bufer
   if(CAp::Len(tmp)<2*n)
      ArrayResizeAL(tmp,2*n);
//--- special cases
   if(n==1)
     {
      //--- check
      if(a[offs][offs]>0.0)
        {
         a[offs].Set(offs,MathSqrt(a[offs][offs]));
         //--- return result
         return(true);
        }
      else
        {
         //--- return result
         return(false);
        }
     }
//--- check
   if(n<=CAblas::AblasBlockSize())
     {
      //--- return result
      return(SPDMatrixCholesky2(a,offs,n,isupper,tmp));
     }
//--- general case: split task in cache-oblivious manner
   result=true;
   CAblas::AblasSplitLength(a,n,n1,n2);
   result=SPDMatrixCholeskyRec(a,offs,n1,isupper,tmp);
//--- check
   if(!result)
     {
      //--- return result
      return(result);
     }
//--- check
   if(n2>0)
     {
      if(isupper)
        {
         //--- function call
         CAblas::RMatrixLeftTrsM(n1,n2,a,offs,offs,isupper,false,1,a,offs,offs+n1);
         //--- function call
         CAblas::RMatrixSyrk(n2,n1,-1.0,a,offs,offs+n1,1,1.0,a,offs+n1,offs+n1,isupper);
        }
      else
        {
         //--- function call
         CAblas::RMatrixRightTrsM(n2,n1,a,offs,offs,isupper,false,1,a,offs+n1,offs);
         //--- function call
         CAblas::RMatrixSyrk(n2,n1,-1.0,a,offs+n1,offs,0,1.0,a,offs+n1,offs+n1,isupper);
        }
      result=SPDMatrixCholeskyRec(a,offs+n1,n2,isupper,tmp);
      //--- check
      if(!result)
        {
         //--- return result
         return(result);
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Recurrent complex LU subroutine.                                 |
//| Never call it directly.                                          |
//+------------------------------------------------------------------+
static void CTrFac::CMatrixLUPRec(CMatrixComplex &a,const int offs,const int m,
                                  const int n,int &pivots[],al_complex &tmp[])
  {
//--- create variables
   int        i=0;
   int        m1=0;
   int        m2=0;
   int        i_=0;
   int        i1_=0;
   al_complex One(1.0,0.0);
   al_complex _One(-1.0,0.0);
//--- Kernel case
   if(MathMin(m,n)<=CAblas::AblasComplexBlockSize())
     {
      CMatrixLUP2(a,offs,m,n,pivots,tmp);
      //--- exit the function
      return;
     }
//--- Preliminary step, make N>=M
//---     ( A1 )
//--- A = (    ), where A1 is square
//---     ( A2 )
//--- Factorize A1, update A2
   if(m>n)
     {
      //--- function call
      CMatrixLUPRec(a,offs,n,n,pivots,tmp);
      for(i=0;i<n;i++)
        {
         i1_=offs+n;
         for(i_=0;i_<m-n;i_++)
            tmp[i_]=a[i_+i1_][offs+i];
         //--- change matrix
         for(i_=offs+n;i_<offs+m;i_++)
            a[i_].Set(offs+i,a[i_][pivots[offs+i]]);
         i1_=-(offs+n);
         for(i_=offs+n;i_<offs+m;i_++)
            a[i_].Set(pivots[offs+i],tmp[i_+i1_]);
        }
      //--- function call
      CAblas::CMatrixRightTrsM(m-n,n,a,offs,offs,true,true,0,a,offs+n,offs);
      //--- exit the function
      return;
     }
//--- Non-kernel case
   CAblas::AblasComplexSplitLength(a,m,m1,m2);
//--- function call
   CMatrixLUPRec(a,offs,m1,n,pivots,tmp);
//--- check
   if(m2>0)
     {
      for(i=0;i<m1;i++)
        {
         //--- check
         if(offs+i!=pivots[offs+i])
           {
            i1_=offs+m1;
            for(i_=0;i_<m2;i_++)
               tmp[i_]=a[i_+i1_][offs+i];
            //--- change matrix
            for(i_=offs+m1;i_<offs+m;i_++)
               a[i_].Set(offs+i,a[i_][pivots[offs+i]]);
            i1_=-(offs+m1);
            for(i_=offs+m1;i_<offs+m;i_++)
               a[i_].Set(pivots[offs+i],tmp[i_+i1_]);
           }
        }
      //--- function call
      CAblas::CMatrixRightTrsM(m2,m1,a,offs,offs,true,true,0,a,offs+m1,offs);
      //--- function call
      CAblas::CMatrixGemm(m-m1,n-m1,m1,_One,a,offs+m1,offs,0,a,offs,offs+m1,0,One,a,offs+m1,offs+m1);
      //--- function call
      CMatrixLUPRec(a,offs+m1,m-m1,n-m1,pivots,tmp);
      for(i=0;i<m2;i++)
        {
         //--- check
         if(offs+m1+i!=pivots[offs+m1+i])
           {
            i1_=offs;
            for(i_=0;i_<m1;i_++)
               tmp[i_]=a[i_+i1_][offs+m1+i];
            //--- change matrix
            for(i_=offs;i_<offs+m1;i_++)
               a[i_].Set(offs+m1+i,a[i_][pivots[offs+m1+i]]);
            i1_=-offs;
            for(i_=offs;i_<offs+m1;i_++)
               a[i_].Set(pivots[offs+m1+i],tmp[i_+i1_]);
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Recurrent real LU subroutine.                                    |
//| Never call it directly.                                          |
//+------------------------------------------------------------------+
static void CTrFac::RMatrixLUPRec(CMatrixDouble &a,const int offs,const int m,
                                  const int n,int &pivots[],double &tmp[])
  {
//--- create variables
   int i=0;
   int m1=0;
   int m2=0;
   int i_=0;
   int i1_=0;
//--- Kernel case
   if(MathMin(m,n)<=CAblas::AblasBlockSize())
     {
      RMatrixLUP2(a,offs,m,n,pivots,tmp);
      //--- exit the function
      return;
     }
//--- Preliminary step, make N>=M
//---     ( A1 )
//--- A = (    ), where A1 is square
//---     ( A2 )
//--- Factorize A1, update A2
   if(m>n)
     {
      //--- function call
      RMatrixLUPRec(a,offs,n,n,pivots,tmp);
      for(i=0;i<n;i++)
        {
         //--- check
         if(offs+i!=pivots[offs+i])
           {
            i1_=offs+n;
            for(i_=0;i_<m-n;i_++)
               tmp[i_]=a[i_+i1_][offs+i];
            //--- change matrix
            for(i_=offs+n;i_<offs+m;i_++)
               a[i_].Set(offs+i,a[i_][pivots[offs+i]]);
            i1_=-(offs+n);
            for(i_=offs+n;i_<offs+m;i_++)
               a[i_].Set(pivots[offs+i],tmp[i_+i1_]);
           }
        }
      //--- function call
      CAblas::RMatrixRightTrsM(m-n,n,a,offs,offs,true,true,0,a,offs+n,offs);
      //--- exit the function
      return;
     }
//--- Non-kernel case
   CAblas::AblasSplitLength(a,m,m1,m2);
//--- function call
   RMatrixLUPRec(a,offs,m1,n,pivots,tmp);
//--- check
   if(m2>0)
     {
      for(i=0;i<m1;i++)
        {
         //--- check
         if(offs+i!=pivots[offs+i])
           {
            i1_=offs+m1;
            for(i_=0;i_<m2;i_++)
               tmp[i_]=a[i_+i1_][offs+i];
            //--- change matrix
            for(i_=offs+m1;i_<offs+m;i_++)
               a[i_].Set(offs+i,a[i_][pivots[offs+i]]);
            i1_=-(offs+m1);
            for(i_=offs+m1;i_<offs+m;i_++)
               a[i_].Set(pivots[offs+i],tmp[i_+i1_]);
           }
        }
      //--- function call
      CAblas::RMatrixRightTrsM(m2,m1,a,offs,offs,true,true,0,a,offs+m1,offs);
      //--- function call
      CAblas::RMatrixGemm(m-m1,n-m1,m1,-1.0,a,offs+m1,offs,0,a,offs,offs+m1,0,1.0,a,offs+m1,offs+m1);
      //--- function call
      RMatrixLUPRec(a,offs+m1,m-m1,n-m1,pivots,tmp);
      for(i=0;i<m2;i++)
        {
         //--- check
         if(offs+m1+i!=pivots[offs+m1+i])
           {
            i1_=offs;
            for(i_=0;i_<m1;i_++)
               tmp[i_]=a[i_+i1_][offs+m1+i];
            //--- change matrix
            for(i_=offs;i_<offs+m1;i_++)
               a[i_].Set(offs+m1+i,a[i_][pivots[offs+m1+i]]);
            i1_=-offs;
            for(i_=offs;i_<offs+m1;i_++)
               a[i_].Set(pivots[offs+m1+i],tmp[i_+i1_]);
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Recurrent complex LU subroutine.                                 |
//| Never call it directly.                                          |
//+------------------------------------------------------------------+
static void CTrFac::CMatrixPLURec(CMatrixComplex &a,const int offs,const int m,const int n,int &pivots[],al_complex &tmp[])
  {
//--- create variables
   int        i=0;
   int        n1=0;
   int        n2=0;
   int        i_=0;
   int        i1_=0;
   al_complex One(1.0,0.0);
   al_complex _One(-1.0,0.0);
//--- Kernel case
   if(MathMin(m,n)<=CAblas::AblasComplexBlockSize())
     {
      //--- function call
      CMatrixPLU2(a,offs,m,n,pivots,tmp);
      //--- exit the function
      return;
     }
//--- Preliminary step, make M>=N.
//--- A = (A1 A2), where A1 is square
//--- Factorize A1, update A2
   if(n>m)
     {
      //--- function call
      CMatrixPLURec(a,offs,m,m,pivots,tmp);
      for(i=0;i<m;i++)
        {
         i1_=offs+m;
         for(i_=0;i_<n-m;i_++)
            tmp[i_]=a[offs+i][i_+i1_];
         //--- change matrix
         for(i_=offs+m;i_<offs+n;i_++)
            a[offs+i].Set(i_,a[pivots[offs+i]][i_]);
         i1_=-(offs+m);
         for(i_=offs+m;i_<offs+n;i_++)
            a[pivots[offs+i]].Set(i_,tmp[i_+i1_]);
        }
      //--- function call
      CAblas::CMatrixLeftTrsM(m,n-m,a,offs,offs,false,true,0,a,offs,offs+m);
      //--- exit the function
      return;
     }
//--- Non-kernel case
   CAblas::AblasComplexSplitLength(a,n,n1,n2);
//--- function call
   CMatrixPLURec(a,offs,m,n1,pivots,tmp);
//--- check
   if(n2>0)
     {
      for(i=0;i<n1;i++)
        {
         //--- check
         if(offs+i!=pivots[offs+i])
           {
            i1_=offs+n1;
            for(i_=0;i_<n2;i_++)
               tmp[i_]=a[offs+i][i_+i1_];
            //--- change matrix
            for(i_=offs+n1;i_<offs+n;i_++)
               a[offs+i].Set(i_,a[pivots[offs+i]][i_]);
            i1_=-(offs+n1);
            for(i_=offs+n1;i_<offs+n;i_++)
               a[pivots[offs+i]].Set(i_,tmp[i_+i1_]);
           }
        }
      //--- function call
      CAblas::CMatrixLeftTrsM(n1,n2,a,offs,offs,false,true,0,a,offs,offs+n1);
      //--- function call
      CAblas::CMatrixGemm(m-n1,n-n1,n1,_One,a,offs+n1,offs,0,a,offs,offs+n1,0,One,a,offs+n1,offs+n1);
      //--- function call
      CMatrixPLURec(a,offs+n1,m-n1,n-n1,pivots,tmp);
      for(i=0;i<n2;i++)
        {
         //--- check
         if(offs+n1+i!=pivots[offs+n1+i])
           {
            i1_=offs;
            for(i_=0;i_<n1;i_++)
               tmp[i_]=a[offs+n1+i][i_+i1_];
            //--- change matrix
            for(i_=offs;i_<offs+n1;i_++)
               a[offs+n1+i].Set(i_,a[pivots[offs+n1+i]][i_]);
            i1_=-offs;
            for(i_=offs;i_<offs+n1;i_++)

               a[pivots[offs+n1+i]].Set(i_,tmp[i_+i1_]);
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Recurrent real LU subroutine.                                    |
//| Never call it directly.                                          |
//+------------------------------------------------------------------+
static void CTrFac::RMatrixPLURec(CMatrixDouble &a,const int offs,const int m,
                                  const int n,int &pivots[],double &tmp[])
  {
//--- create variables
   int i=0;
   int n1=0;
   int n2=0;
   int i_=0;
   int i1_=0;
//--- Kernel case
   if(MathMin(m,n)<=CAblas::AblasBlockSize())
     {
      //--- function call
      RMatrixPLU2(a,offs,m,n,pivots,tmp);
      //--- exit the function
      return;
     }
//--- Preliminary step, make M>=N.
//--- A = (A1 A2), where A1 is square
//--- Factorize A1, update A2
   if(n>m)
     {
      //--- function call
      RMatrixPLURec(a,offs,m,m,pivots,tmp);
      for(i=0;i<m;i++)
        {
         i1_=offs+m;
         for(i_=0;i_<n-m;i_++)
            tmp[i_]=a[offs+i][i_+i1_];
         //--- change matrix
         for(i_=offs+m;i_<offs+n;i_++)
            a[offs+i].Set(i_,a[pivots[offs+i]][i_]);
         i1_=-(offs+m);
         for(i_=offs+m;i_<offs+n;i_++)
            a[pivots[offs+i]].Set(i_,tmp[i_+i1_]);
        }
      //--- function call
      CAblas::RMatrixLeftTrsM(m,n-m,a,offs,offs,false,true,0,a,offs,offs+m);
      //--- exit the function
      return;
     }
//--- Non-kernel case
   CAblas::AblasSplitLength(a,n,n1,n2);
//--- function call
   RMatrixPLURec(a,offs,m,n1,pivots,tmp);
//--- check
   if(n2>0)
     {
      for(i=0;i<n1;i++)
        {
         //--- check
         if(offs+i!=pivots[offs+i])
           {
            i1_=offs+n1;
            for(i_=0;i_<n2;i_++)
               tmp[i_]=a[offs+i][i_+i1_];
            //--- change matrix
            for(i_=offs+n1;i_<offs+n;i_++)
               a[offs+i].Set(i_,a[pivots[offs+i]][i_]);
            i1_=-(offs+n1);
            for(i_=offs+n1;i_<offs+n;i_++)
               a[pivots[offs+i]].Set(i_,tmp[i_+i1_]);
           }
        }
      //--- function call
      CAblas::RMatrixLeftTrsM(n1,n2,a,offs,offs,false,true,0,a,offs,offs+n1);
      //--- function call
      CAblas::RMatrixGemm(m-n1,n-n1,n1,-1.0,a,offs+n1,offs,0,a,offs,offs+n1,0,1.0,a,offs+n1,offs+n1);
      //--- function call
      RMatrixPLURec(a,offs+n1,m-n1,n-n1,pivots,tmp);
      for(i=0;i<n2;i++)
        {
         //--- check
         if(offs+n1+i!=pivots[offs+n1+i])
           {
            i1_=offs;
            for(i_=0;i_<n1;i_++)
               tmp[i_]=a[offs+n1+i][i_+i1_];
            //--- change matrix
            for(i_=offs;i_<offs+n1;i_++)
               a[offs+n1+i].Set(i_,a[pivots[offs+n1+i]][i_]);
            i1_=-offs;
            for(i_=offs;i_<offs+n1;i_++)
               a[pivots[offs+n1+i]].Set(i_,tmp[i_+i1_]);
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Complex LUP kernel                                               |
//+------------------------------------------------------------------+
static void CTrFac::CMatrixLUP2(CMatrixComplex &a,const int offs,const int m,const int n,int &pivots[],al_complex &tmp[])
  {
//--- create variables
   int        i=0;
   int        j=0;
   int        jp=0;
   al_complex s=0;
   int        i_=0;
   int        i1_=0;
   al_complex zero=0;
   al_complex One(1.0,0.0);
//--- check
   if(m==0 || n==0)
      return;
//--- main cycle
   for(j=0;j<=MathMin(m-1,n-1);j++)
     {
      //--- Find pivot, swap columns
      jp=j;
      for(i=j+1;i<n;i++)
        {
         //--- check
         if(CMath::AbsComplex(a[offs+j][offs+i])>CMath::AbsComplex(a[offs+j][offs+jp]))
            jp=i;
        }
      pivots[offs+j]=offs+jp;
      //--- check
      if(jp!=j)
        {
         i1_=offs;
         for(i_=0;i_<m;i_++)
            tmp[i_]=a[i_+i1_][offs+j];
         //--- change matrix
         for(i_=offs;i_<offs+m;i_++)
            a[i_].Set(offs+j,a[i_][offs+jp]);
         i1_=-offs;
         for(i_=offs;i_<offs+m;i_++)
            a[i_].Set(offs+jp,tmp[i_+i1_]);
        }
      //--- LU decomposition of 1x(N-J) matrix
      if(a[offs+j][offs+j]!=zero && j+1<=n-1)
        {
         s=One/a[offs+j][offs+j];
         for(i_=offs+j+1;i_<offs+n;i_++)
            a[offs+j].Set(i_,s*a[offs+j][i_]);
        }
      //--- Update trailing (M-J-1)x(N-J-1) matrix
      if(j<MathMin(m-1,n-1))
        {
         i1_=offs+j+1;
         for(i_=0;i_<m-j-1;i_++)
            tmp[i_]=a[i_+i1_][offs+j];
         i1_=(offs+j+1)-(m);
         for(i_=m;i_<m+n-j-1;i_++)
           {
            tmp[i_]=a[offs+j][i_+i1_];
            tmp[i_].re=-tmp[i_].re;
            tmp[i_].im=-tmp[i_].im;
           }
         //--- function call
         CAblas::CMatrixRank1(m-j-1,n-j-1,a,offs+j+1,offs+j+1,tmp,0,tmp,m);
        }
     }
  }
//+------------------------------------------------------------------+
//| Real LUP kernel                                                  |
//+------------------------------------------------------------------+
static void CTrFac::RMatrixLUP2(CMatrixDouble &a,const int offs,const int m,
                                const int n,int &pivots[],double &tmp[])
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    jp=0;
   double s=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(m==0 || n==0)
      return;
//--- main cycle
   for(j=0;j<=MathMin(m-1,n-1);j++)
     {
      //--- Find pivot, swap columns
      jp=j;
      for(i=j+1;i<n;i++)
        {
         if(MathAbs(a[offs+j][offs+i])>MathAbs(a[offs+j][offs+jp]))
            jp=i;
        }
      pivots[offs+j]=offs+jp;
      //--- check
      if(jp!=j)
        {
         i1_=offs;
         for(i_=0;i_<m;i_++)
            tmp[i_]=a[i_+i1_][offs+j];
         //--- change matrix
         for(i_=offs;i_<offs+m;i_++)
            a[i_].Set(offs+j,a[i_][offs+jp]);
         i1_=-offs;
         for(i_=offs;i_<offs+m;i_++)
            a[i_].Set(offs+jp,tmp[i_+i1_]);
        }
      //--- LU decomposition of 1x(N-J) matrix
      if(a[offs+j][offs+j]!=0.0 && j+1<=n-1)
        {
         s=1/a[offs+j][offs+j];
         for(i_=offs+j+1;i_<offs+n;i_++)
            a[offs+j].Set(i_,s*a[offs+j][i_]);
        }
      //--- Update trailing (M-J-1)x(N-J-1) matrix
      if(j<MathMin(m-1,n-1))
        {
         i1_=offs+j+1;
         for(i_=0;i_<m-j-1;i_++)
            tmp[i_]=a[i_+i1_][offs+j];
         i1_=(offs+j+1)-(m);
         for(i_=m;i_<m+n-j-1;i_++)
            tmp[i_]=-a[offs+j][i_+i1_];
         //--- function call
         CAblas::RMatrixRank1(m-j-1,n-j-1,a,offs+j+1,offs+j+1,tmp,0,tmp,m);
        }
     }
  }
//+------------------------------------------------------------------+
//| Complex PLU kernel                                               |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee,Univ. of California Berkeley,NAG Ltd.,   |
//|      Courant Institute,Argonne National Lab,and Rice University  |
//|      June 30,1992                                                |
//+------------------------------------------------------------------+
static void CTrFac::CMatrixPLU2(CMatrixComplex &a,const int offs,const int m,
                                const int n,int &pivots[],al_complex &tmp[])
  {
//--- create variables
   int        i=0;
   int        j=0;
   int        jp=0;
   al_complex s=0;
   al_complex zero=0;
   al_complex One(1.0,0.0);
   int        i_=0;
   int        i1_=0;
//--- check
   if(m==0 || n==0)
      return;
   for(j=0;j<=MathMin(m-1,n-1);j++)
     {
      //---  Find pivot and test for singularity.
      jp=j;
      for(i=j+1;i<m;i++)
        {
         //--- check
         if(CMath::AbsComplex(a[offs+i][offs+j])>CMath::AbsComplex(a[offs+jp][offs+j]))
            jp=i;
        }
      pivots[offs+j]=offs+jp;
      if(a[offs+jp][offs+j]!=zero)
        {
         //---  Apply the interchange to rows
         if(jp!=j)
           {
            for(i=0;i<n;i++)
              {
               s=a[offs+j][offs+i];
               a[offs+j].Set(offs+i,a[offs+jp][offs+i]);
               a[offs+jp].Set(offs+i,s);
              }
           }
         //--- Compute elements J+1:M of J-th column.
         if(j+1<=m-1)
           {
            s=One/a[offs+j][offs+j];
            for(i_=offs+j+1;i_<offs+m;i_++)
               a[i_].Set(offs+j,s*a[i_][offs+j]);
           }
        }
      //--- check
      if(j<MathMin(m,n)-1)
        {
         //--- Update trailing submatrix.
         i1_=offs+j+1;
         for(i_=0;i_<m-j-1;i_++)
            tmp[i_]=a[i_+i1_][offs+j];
         i1_=(offs+j+1)-(m);
         for(i_=m;i_<m+n-j-1;i_++)
           {
            tmp[i_]=a[offs+j][i_+i1_];
            tmp[i_].re=-tmp[i_].re;
            tmp[i_].im=-tmp[i_].im;
           }
         //--- function call
         CAblas::CMatrixRank1(m-j-1,n-j-1,a,offs+j+1,offs+j+1,tmp,0,tmp,m);
        }
     }
  }
//+------------------------------------------------------------------+
//| Real PLU kernel                                                  |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee,Univ. of California Berkeley,NAG Ltd.,   |
//|      Courant Institute,Argonne National Lab,and Rice University  |
//|      June 30,1992                                                |
//+------------------------------------------------------------------+
static void CTrFac::RMatrixPLU2(CMatrixDouble &a,const int offs,const int m,
                                const int n,int &pivots[],double &tmp[])
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    jp=0;
   double s=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(m==0 || n==0)
      return;
   for(j=0;j<=MathMin(m-1,n-1);j++)
     {
      //--- Find pivot and test for singularity.
      jp=j;
      for(i=j+1;i<m;i++)
        {
         //--- check
         if(MathAbs(a[offs+i][offs+j])>MathAbs(a[offs+jp][offs+j]))
            jp=i;
        }
      pivots[offs+j]=offs+jp;
      //--- check
      if(a[offs+jp][offs+j]!=0.0)
        {
         //--- Apply the interchange to rows
         if(jp!=j)
           {
            for(i=0;i<n;i++)
              {
               s=a[offs+j][offs+i];
               a[offs+j].Set(offs+i,a[offs+jp][offs+i]);
               a[offs+jp].Set(offs+i,s);
              }
           }
         //--- Compute elements J+1:M of J-th column.
         if(j+1<=m-1)
           {
            s=1/a[offs+j][offs+j];
            for(i_=offs+j+1;i_<offs+m;i_++)
               a[i_].Set(offs+j,s*a[i_][offs+j]);
           }
        }
      //--- check
      if(j<MathMin(m,n)-1)
        {
         //--- Update trailing submatrix.
         i1_=offs+j+1;
         for(i_=0;i_<m-j-1;i_++)
            tmp[i_]=a[i_+i1_][offs+j];
         i1_=(offs+j+1)-(m);
         for(i_=m;i_<m+n-j-1;i_++)
            tmp[i_]=-a[offs+j][i_+i1_];
         //--- function call
         CAblas::RMatrixRank1(m-j-1,n-j-1,a,offs+j+1,offs+j+1,tmp,0,tmp,m);
        }
     }
  }
//+------------------------------------------------------------------+
//| Recursive computational subroutine for HPDMatrixCholesky         |
//+------------------------------------------------------------------+
static bool CTrFac::HPDMatrixCholeskyRec(CMatrixComplex &a,const int offs,const int n,
                                         const bool isupper,al_complex &tmp[])
  {
//--- create variables
   bool result;
   int  n1=0;
   int  n2=0;
//--- check
   if(n<1)
      return(false);
//--- prepare bufer
   if(CAp::Len(tmp)<2*n)
      ArrayResizeAL(tmp,2*n);
//--- special cases
   if(n==1)
     {
      //--- check
      if(a[offs][offs].re>0.0)
        {
         a[offs].Set(offs,MathSqrt(a[offs][offs].re));
         result=true;
        }
      else
        {
         result=false;
        }
      //--- return result
      return(result);
     }
//--- check
   if(n<=CAblas::AblasComplexBlockSize())
     {
      result=HPDMatrixCholesky2(a,offs,n,isupper,tmp);
      //--- return result
      return(result);
     }
//--- general case: split task in cache-oblivious manner
   result=true;
//--- function call
   CAblas::AblasComplexSplitLength(a,n,n1,n2);
   result=HPDMatrixCholeskyRec(a,offs,n1,isupper,tmp);
//--- check
   if(!result)
      return(result);
//--- check
   if(n2>0)
     {
      //--- check
      if(isupper)
        {
         //--- function call
         CAblas::CMatrixLeftTrsM(n1,n2,a,offs,offs,isupper,false,2,a,offs,offs+n1);
         //--- function call
         CAblas::CMatrixSyrk(n2,n1,-1.0,a,offs,offs+n1,2,1.0,a,offs+n1,offs+n1,isupper);
        }
      else
        {
         //--- function call
         CAblas::CMatrixRightTrsM(n2,n1,a,offs,offs,isupper,false,2,a,offs+n1,offs);
         //--- function call
         CAblas::CMatrixSyrk(n2,n1,-1.0,a,offs+n1,offs,0,1.0,a,offs+n1,offs+n1,isupper);
        }
      result=HPDMatrixCholeskyRec(a,offs+n1,n2,isupper,tmp);
      //--- check
      if(!result)
         return(result);
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Level-2 Hermitian Cholesky subroutine.                           |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee,Univ. of California Berkeley,NAG Ltd.,   |
//|      Courant Institute,Argonne National Lab,and Rice University  |
//|      February 29,1992                                            |
//+------------------------------------------------------------------+
static bool CTrFac::HPDMatrixCholesky2(CMatrixComplex &aaa,const int offs,
                                       const int n,const bool isupper,al_complex &tmp[])
  {
//--- create variables
   bool       result;
   int        i=0;
   int        j=0;
   double     ajj=0;
   al_complex v=0;
   al_complex cR;
   double     r=0;
   int        i_=0;
   int        i1_=0;
//--- initialization
   result=true;
//--- check
   if(n<0)
      return(false);
//--- check
   if(n==0)
      return(result);
//--- check
   if(isupper)
     {
      //--- Compute the Cholesky factorization A = U'*U.
      for(j=0;j<n;j++)
        {
         //--- Compute U(J,J) and test for non-positive-definiteness.
         v=0.0;
         for(i_=offs;i_<offs+j;i_++)
            v+=CMath::Conj(aaa[i_][offs+j])*aaa[i_][offs+j];
         ajj=(aaa[offs+j][offs+j]-v).re;
         //--- check
         if(ajj<=0.0)
           {
            aaa[offs+j].Set(offs+j,ajj);
            //--- return result
            return(false);
           }
         ajj=MathSqrt(ajj);
         aaa[offs+j].Set(offs+j,ajj);
         //--- Compute elements J+1:N-1 of row J.
         if(j<n-1)
           {
            //--- check
            if(j>0)
              {
               i1_=offs;
               for(i_=0;i_<j;i_++)
                 {
                  tmp[i_]=CMath::Conj(aaa[i_+i1_][offs+j]);
                  tmp[i_].re=-tmp[i_].re;
                  tmp[i_].im=-tmp[i_].im;
                 }
               //--- function call
               CAblas::CMatrixMVect(n-j-1,j,aaa,offs,offs+j+1,1,tmp,0,tmp,n);
               i1_=(n)-(offs+j+1);
               for(i_=offs+j+1;i_<offs+n;i_++)
                  aaa[offs+j].Set(i_,aaa[offs+j][i_]+tmp[i_+i1_]);
              }
            r=1/ajj;
            cR=r;
            //--- change matrix
            for(i_=offs+j+1;i_<offs+n;i_++)
               aaa[offs+j].Set(i_,cR*aaa[offs+j][i_]);
           }
        }
     }
   else
     {
      //--- Compute the Cholesky factorization A = L*L'.
      for(j=0;j<n;j++)
        {
         //--- Compute L(J+1,J+1) and test for non-positive-definiteness.
         v=0.0;
         for(i_=offs;i_<offs+j;i_++)
            v+=CMath::Conj(aaa[offs+j][i_])*aaa[offs+j][i_];
         ajj=(aaa[offs+j][offs+j]-v).re;
         //--- check
         if(ajj<=0.0)
           {
            aaa[offs+j].Set(offs+j,ajj);
            //--- return result
            return(false);
           }
         ajj=MathSqrt(ajj);
         aaa[offs+j].Set(offs+j,ajj);
         //--- Compute elements J+1:N of column J.
         if(j<n-1)
           {
            //--- check
            if(j>0)
              {
               i1_=offs;
               for(i_=0;i_<j;i_++)
                  tmp[i_]=CMath::Conj(aaa[offs+j][i_+i1_]);
               //--- function call
               CAblas::CMatrixMVect(n-j-1,j,aaa,offs+j+1,offs,0,tmp,0,tmp,n);
               for(i=0;i<n-j-1;i++)
                  aaa[offs+j+1+i].Set(offs+j,(aaa[offs+j+1+i][offs+j]-tmp[n+i])/ajj);
              }
            else
              {
               for(i=0;i<n-j-1;i++)
                  aaa[offs+j+1+i].Set(offs+j,aaa[offs+j+1+i][offs+j]/ajj);
              }
           }
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Level-2 Cholesky subroutine                                      |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      February 29, 1992                                           |
//+------------------------------------------------------------------+
static bool CTrFac::SPDMatrixCholesky2(CMatrixDouble &aaa,const int offs,const int n,
                                       const bool isupper,double &tmp[])
  {
//--- create variables
   bool   result;
   int    i=0;
   int    j=0;
   double ajj=0;
   double v=0;
   double r=0;
   int    i_=0;
   int    i1_=0;
//--- initialization
   result=true;
//--- check
   if(n<0)
      return(false);
//--- check
   if(n==0)
      return(result);
//--- check
   if(isupper)
     {
      //--- Compute the Cholesky factorization A = U'*U.
      for(j=0;j<n;j++)
        {
         //--- Compute U(J,J) and test for non-positive-definiteness.
         v=0.0;
         for(i_=offs;i_<offs+j;i_++)
            v+=aaa[i_][offs+j]*aaa[i_][offs+j];
         ajj=aaa[offs+j][offs+j]-v;
         //--- check
         if(ajj<=0.0)
           {
            aaa[offs+j].Set(offs+j,ajj);
            result=false;
            //--- return result
            return(result);
           }
         ajj=MathSqrt(ajj);
         aaa[offs+j].Set(offs+j,ajj);
         //--- Compute elements J+1:N-1 of row J.
         if(j<n-1)
           {
            //--- check
            if(j>0)
              {
               i1_=offs;
               for(i_=0;i_<j;i_++)
                  tmp[i_]=-aaa[i_+i1_][offs+j];
               //--- function call
               CAblas::RMatrixMVect(n-j-1,j,aaa,offs,offs+j+1,1,tmp,0,tmp,n);
               i1_=(n)-(offs+j+1);
               for(i_=offs+j+1;i_<offs+n;i_++)
                  aaa[offs+j].Set(i_,aaa[offs+j][i_]+tmp[i_+i1_]);
              }
            r=1/ajj;
            //--- change matrix
            for(i_=offs+j+1;i_<offs+n;i_++)
               aaa[offs+j].Set(i_,r*aaa[offs+j][i_]);
           }
        }
     }
   else
     {
      //--- Compute the Cholesky factorization A = L*L'.
      for(j=0;j<n;j++)
        {
         //--- Compute L(J+1,J+1) and test for non-positive-definiteness.
         v=0.0;
         for(i_=offs;i_<offs+j;i_++)
            v+=aaa[offs+j][i_]*aaa[offs+j][i_];
         ajj=aaa[offs+j][offs+j]-v;
         //--- check
         if(ajj<=0.0)
           {
            aaa[offs+j].Set(offs+j,ajj);
            //--- return result
            return(false);
           }
         ajj=MathSqrt(ajj);
         aaa[offs+j].Set(offs+j,ajj);
         //--- Compute elements J+1:N of column J.
         if(j<n-1)
           {
            //--- check
            if(j>0)
              {
               i1_=offs;
               for(i_=0;i_<j;i_++)
                  tmp[i_]=aaa[offs+j][i_+i1_];
               //--- function call
               CAblas::RMatrixMVect(n-j-1,j,aaa,offs+j+1,offs,0,tmp,0,tmp,n);
               for(i=0;i<n-j-1;i++)
                  aaa[offs+j+1+i].Set(offs+j,(aaa[offs+j+1+i][offs+j]-tmp[n+i])/ajj);
              }
            else
              {
               for(i=0;i<n-j-1;i++)
                  aaa[offs+j+1+i].Set(offs+j,aaa[offs+j+1+i][offs+j]/ajj);
              }
           }
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Estimate of a matrix condition number                            |
//+------------------------------------------------------------------+
class CRCond
  {
private:
   //--- private methods
   static void       RMatrixRCondTrInternal(CMatrixDouble &a,const int n,const bool isupper,const bool isunit,const bool onenorm,double anorm,double &rc);
   static void       CMatrixRCondTrInternal(CMatrixComplex &a,const int n,const bool isupper,const bool isunit,const bool onenorm,double anorm,double &rc);
   static void       SPDMatrixRCondCholeskyInternal(CMatrixDouble &cha,const int n,const bool isupper,const bool isnormprovided,double anorm,double &rc);
   static void       HPDMatrixRCondCholeskyInternal(CMatrixComplex &cha,const int n,const bool isupper,const bool isnormprovided,double anorm,double &rc);
   static void       RMatrixRCondLUInternal(CMatrixDouble &lua,const int n,const bool onenorm,const bool isanormprovided,double anorm,double &rc);
   static void       CMatrixRCondLUInternal(CMatrixComplex &lua,const int n,const bool onenorm,const bool isanormprovided,double anorm,double &rc);
   static void       RMatrixEstimateNorm(const int n,double &v[],double &x[],int &isgn[],double &est,int &kase);
   static void       CMatrixEstimateNorm(const int n,al_complex &v[],al_complex &x[],double &est,int &kase,int &isave[],double &rsave[]);
   static double     InternalComplexRCondScSum1(al_complex &x[],const int n);
   static int        InternalComplexRCondIcMax1(al_complex &x[],const int n);
   static void       InternalComplexRCondSaveAll(int &isave[],double &rsave[],int &i,int &iter,int &j,int &jlast,int &jump,double &absxi,double &altsgn,double &estold,double &temp);
   static void       InternalComplexRCondLoadAll(int &isave[],double &rsave[],int &i,int &iter,int &j,int &jlast,int &jump,double &absxi,double &altsgn,double &estold,double &temp);
public:
   //--- public methods
   static double     RMatrixRCond1(CMatrixDouble &ca,const int n);
   static double     RMatrixRCondInf(CMatrixDouble &ca,const int n);
   static double     SPDMatrixRCond(CMatrixDouble &ca,const int n,const bool isupper);
   static double     RMatrixTrRCond1(CMatrixDouble &a,const int n,const bool isupper,const bool isunit);
   static double     RMatrixTrRCondInf(CMatrixDouble &a,const int n,const bool isupper,const bool isunit);
   static double     RMatrixLURCond1(CMatrixDouble &lua,const int n);
   static double     RMatrixLURCondInf(CMatrixDouble &lua,const int n);
   static double     SPDMatrixCholeskyRCond(CMatrixDouble &a,const int n,const bool isupper);
   static double     HPDMatrixRCond(CMatrixComplex &ca,const int n,const bool isupper);
   static double     CMatrixRCond1(CMatrixComplex &ca,const int n);
   static double     CMatrixRCondInf(CMatrixComplex &ca,const int n);
   static double     HPDMatrixCholeskyRCond(CMatrixComplex &a,const int n,const bool isupper);
   static double     CMatrixLURCond1(CMatrixComplex &lua,const int n);
   static double     CMatrixLURCondInf(CMatrixComplex &lua,const int n);
   static double     CMatrixTrRCond1(CMatrixComplex &a,const int n,const bool isupper,const bool isunit);
   static double     CMatrixTrRCondInf(CMatrixComplex &a,const int n,const bool isupper,const bool isunit);
   static double     RCondThreshold(void);
  };
//+------------------------------------------------------------------+
//| Estimate of a matrix condition number (1-norm)                   |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| Input parameters:                                                |
//|     A   -   matrix. Array whose indexes range within             |
//|             [0..N-1, 0..N-1].                                    |
//|     N   -   size of matrix A.                                    |
//| Result: 1/LowerBound(cond(A))                                    |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::RMatrixRCond1(CMatrixDouble &ca,const int n)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   double nrm=0;
//--- create arrays
   int    pivots[];
   double t[];
//--- create copy
   CMatrixDouble a;
   a=ca;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- allocation
   ArrayResizeAL(t,n);
//--- fiiling array
   for(i=0;i<n;i++)
      t[i]=0;
   for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
         t[j]=t[j]+MathAbs(a[i][j]);
     }
//--- change values
   nrm=0;
   for(i=0;i<n;i++)
      nrm=MathMax(nrm,t[i]);
//--- function call
   CTrFac::RMatrixLU(a,n,n,pivots);
//--- function call
   RMatrixRCondLUInternal(a,n,true,true,nrm,v);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Estimate of a matrix condition number (infinity-norm).           |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| Input parameters:                                                |
//|     A   -   matrix. Array whose indexes range within             |
//|     [0..N-1, 0..N-1].                                            |
//|     N   -   size of matrix A.                                    |
//| Result: 1/LowerBound(cond(A))                                    |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::RMatrixRCondInf(CMatrixDouble &ca,const int n)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   double nrm=0;
//--- create array
   int pivots[];
//--- create copy
   CMatrixDouble a;
   a=ca;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- change values
   nrm=0;
   for(i=0;i<n;i++)
     {
      v=0;
      for(j=0;j<n;j++)
         v=v+MathAbs(a[i][j]);
      nrm=MathMax(nrm,v);
     }
//--- function call
   CTrFac::RMatrixLU(a,n,n,pivots);
//--- function call
   RMatrixRCondLUInternal(a,n,false,true,nrm,v);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Condition number estimate of a symmetric positive definite       |
//| matrix.                                                          |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| It should be noted that 1-norm and inf-norm of condition numbers |
//| of symmetric matrices are equal, so the algorithm doesn't take   |
//| into account the differences between these types of norms.       |
//| Input parameters:                                                |
//|     A       -   symmetric positive definite matrix which is given|
//|                 by its upper or lower triangle depending on the  |
//|                 value of IsUpper. Array with elements            |
//|                 [0..N-1, 0..N-1].                                |
//|     N       -   size of matrix A.                                |
//|     IsUpper -   storage format.                                  |
//| Result:                                                          |
//|     1/LowerBound(cond(A)), if matrix A is positive definite,     |
//|    -1, if matrix A is not positive definite, and its condition   |
//|     number could not be found by this algorithm.                 |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::SPDMatrixRCond(CMatrixDouble &ca,const int n,const bool isupper)
  {
//--- create variables
   double result=0;
   int    i=0;
   int    j=0;
   int    j1=0;
   int    j2=0;
   double v=0;
   double nrm=0;
//--- create array
   double t[];
//--- create copy
   CMatrixDouble a;
   a=ca;
//--- allocation
   ArrayResizeAL(t,n);
//--- fiiling array
   for(i=0;i<n;i++)
      t[i]=0;
   for(i=0;i<n;i++)
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
      //--- change t
      for(j=j1;j<=j2;j++)
        {
         //--- check
         if(i==j)
            t[i]=t[i]+MathAbs(a[i][i]);
         else
           {
            t[i]=t[i]+MathAbs(a[i][j]);
            t[j]=t[j]+MathAbs(a[i][j]);
           }
        }
     }
//--- change values
   nrm=0;
   for(i=0;i<n;i++)
      nrm=MathMax(nrm,t[i]);
//--- check
   if(CTrFac::SPDMatrixCholesky(a,n,isupper))
     {
      //--- function call
      SPDMatrixRCondCholeskyInternal(a,n,isupper,true,nrm,v);
      //--- get result
      result=v;
     }
   else
      result=-1;
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Triangular matrix: estimate of a condition number (1-norm)       |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| Input parameters:                                                |
//|     A       -   matrix. Array[0..N-1, 0..N-1].                   |
//|     N       -   size of A.                                       |
//|     IsUpper -   True, if the matrix is upper triangular.         |
//|     IsUnit  -   True, if the matrix has a unit diagonal.         |
//| Result: 1/LowerBound(cond(A))                                    |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::RMatrixTrRCond1(CMatrixDouble &a,const int n,
                                      const bool isupper,const bool isunit)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   double nrm=0;
   int    j1=0;
   int    j2=0;
//--- create arrays
   int    pivots[];
   double t[];
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- allocation
   ArrayResizeAL(t,n);
//--- fiiling array
   for(i=0;i<n;i++)
      t[i]=0;
   for(i=0;i<n;i++)
     {
      //--- check
      if(isupper)
        {
         j1=i+1;
         j2=n-1;
        }
      else
        {
         j1=0;
         j2=i-1;
        }
      //--- change t
      for(j=j1;j<=j2;j++)
         t[j]=t[j]+MathAbs(a[i][j]);
      //--- check
      if(isunit)
         t[i]=t[i]+1;
      else
         t[i]=t[i]+MathAbs(a[i][i]);
     }
//--- change values
   nrm=0;
   for(i=0;i<n;i++)
      nrm=MathMax(nrm,t[i]);
//--- function call
   RMatrixRCondTrInternal(a,n,isupper,isunit,true,nrm,v);
//--- return result
   return(v);
  }
//+-------------------------------------------------------------------+
//| Triangular matrix: estimate of a matrix condition number          |
//| (infinity-norm).                                                  |
//| The algorithm calculates a lower bound of the condition number. In|
//| this case, the algorithm does not return a lower bound of the     |
//| condition number, but an inverse number (to avoid an overflow in  |
//| case of a singular matrix).                                       |
//| Input parameters:                                                 |
//|     A   -   matrix. Array whose indexes range within              |
//|             [0..N-1, 0..N-1].                                     |
//|     N   -   size of matrix A.                                     |
//|     IsUpper -   True, if the matrix is upper triangular.          |
//|     IsUnit  -   True, if the matrix has a unit diagonal.          |
//| Result: 1/LowerBound(cond(A))                                     |
//| NOTE:                                                             |
//|     if k(A) is very large, then matrix is assumed degenerate,     |
//|     k(A)=INF, 0.0 is returned in such cases.                      |
//+-------------------------------------------------------------------+
static double CRCond::RMatrixTrRCondInf(CMatrixDouble &a,const int n,
                                        const bool isupper,const bool isunit)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   double nrm=0;
   int    j1=0;
   int    j2=0;
//--- create array
   int pivots[];
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- change values
   nrm=0;
   for(i=0;i<n;i++)
     {
      //--- check
      if(isupper)
        {
         j1=i+1;
         j2=n-1;
        }
      else
        {
         j1=0;
         j2=i-1;
        }
      //--- change v
      v=0;
      for(j=j1;j<=j2;j++)
         v=v+MathAbs(a[i][j]);
      //--- check
      if(isunit)
         v=v+1;
      else
         v=v+MathAbs(a[i][i]);
      nrm=MathMax(nrm,v);
     }
//--- function call
   RMatrixRCondTrInternal(a,n,isupper,isunit,false,nrm,v);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Condition number estimate of a Hermitian positive definite       |
//| matrix.                                                          |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| It should be noted that 1-norm and inf-norm of condition numbers |
//| of symmetric matrices are equal, so the algorithm doesn't take   |
//| into account the differences between these types of norms.       |
//| Input parameters:                                                |
//|     A       -   Hermitian positive definite matrix which is given|
//|                 by its upper or lower triangle depending on the  |
//|                 value of IsUpper. Array with elements            |
//|                 [0..N-1, 0..N-1].                                |
//|     N       -   size of matrix A.                                |
//|     IsUpper -   storage format.                                  |
//| Result:                                                          |
//|     1/LowerBound(cond(A)), if matrix A is positive definite,     |
//|    -1, if matrix A is not positive definite, and its condition   |
//|     number could not be found by this algorithm.                 |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::HPDMatrixRCond(CMatrixComplex &ca,const int n,const bool isupper)
  {
//--- create variables
   double result=0;
   int    i=0;
   int    j=0;
   int    j1=0;
   int    j2=0;
   double v=0;
   double nrm=0;
//--- create array
   double t[];
//--- create copy
   CMatrixComplex a;
   a=ca;
//--- allocation
   ArrayResizeAL(t,n);
//--- fiiling array
   for(i=0;i<n;i++)
      t[i]=0;
   for(i=0;i<n;i++)
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
      //--- change t
      for(j=j1;j<=j2;j++)
        {
         //--- check
         if(i==j)
            t[i]=t[i]+CMath::AbsComplex(a[i][i]);
         else
           {
            t[i]=t[i]+CMath::AbsComplex(a[i][j]);
            t[j]=t[j]+CMath::AbsComplex(a[i][j]);
           }
        }
     }
//--- change values
   nrm=0;
   for(i=0;i<n;i++)
      nrm=MathMax(nrm,t[i]);
//--- check
   if(CTrFac::HPDMatrixCholesky(a,n,isupper))
     {
      //--- function call
      HPDMatrixRCondCholeskyInternal(a,n,isupper,true,nrm,v);
      //--- get result
      result=v;
     }
   else
      result=-1;
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Estimate of a matrix condition number (1-norm)                   |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| Input parameters:                                                |
//|     A   -   matrix. Array whose indexes range within             |
//|             [0..N-1, 0..N-1].                                    |
//|     N   -   size of matrix A.                                    |
//| Result: 1/LowerBound(cond(A))                                    |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::CMatrixRCond1(CMatrixComplex &ca,const int n)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   double nrm=0;
//--- create arrays
   double t[];
   int    pivots[];
//--- create copy
   CMatrixComplex a;
   a=ca;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- allocation
   ArrayResizeAL(t,n);
//--- fiiling array
   for(i=0;i<n;i++)
      t[i]=0;
   for(i=0;i<n;i++)
     {
      for(j=0;j<n;j++)
         t[j]=t[j]+CMath::AbsComplex(a[i][j]);
     }
//--- change values
   nrm=0;
   for(i=0;i<n;i++)
      nrm=MathMax(nrm,t[i]);
//--- function call
   CTrFac::CMatrixLU(a,n,n,pivots);
//--- function call
   CMatrixRCondLUInternal(a,n,true,true,nrm,v);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Estimate of a matrix condition number (infinity-norm).           |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| Input parameters:                                                |
//|     A   -   matrix. Array whose indexes range within             |
//|             [0..N-1, 0..N-1].                                    |
//|     N   -   size of matrix A.                                    |
//| Result: 1/LowerBound(cond(A))                                    |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::CMatrixRCondInf(CMatrixComplex &ca,const int n)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   double nrm=0;
//--- create array
   int pivots[];
//--- create copy
   CMatrixComplex a;
   a=ca;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- change values
   nrm=0;
   for(i=0;i<n;i++)
     {
      v=0;
      for(j=0;j<n;j++)
         v=v+CMath::AbsComplex(a[i][j]);
      nrm=MathMax(nrm,v);
     }
//--- function call
   CTrFac::CMatrixLU(a,n,n,pivots);
//--- function call
   CMatrixRCondLUInternal(a,n,false,true,nrm,v);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Estimate of the condition number of a matrix given by its LU     |
//| decomposition (1-norm)                                           |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| Input parameters:                                                |
//|     LUA         -   LU decomposition of a matrix in compact form.|
//|                     Output of the RMatrixLU subroutine.          |
//|     N           -   size of matrix A.                            |
//| Result: 1/LowerBound(cond(A))                                    |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::RMatrixLURCond1(CMatrixDouble &lua,const int n)
  {
//--- create a variable
   double v=0;
//--- function call
   RMatrixRCondLUInternal(lua,n,true,false,0,v);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Estimate of the condition number of a matrix given by its LU     |
//| decomposition (infinity norm).                                   |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| Input parameters:                                                |
//|     LUA     -   LU decomposition of a matrix in compact form.    |
//|                 Output of the RMatrixLU subroutine.              |
//|     N       -   size of matrix A.                                |
//| Result: 1/LowerBound(cond(A))                                    |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is  assumed  degenerate,  |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::RMatrixLURCondInf(CMatrixDouble &lua,const int n)
  {
//--- create a variable
   double v=0;
//--- function call
   RMatrixRCondLUInternal(lua,n,false,false,0,v);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Condition number estimate of a symmetric positive definite matrix|
//| given by Cholesky decomposition.                                 |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| It should be noted that 1-norm and inf-norm condition numbers of |
//| symmetric matrices are equal, so the algorithm doesn't take into |
//| account the differences between these types of norms.            |
//| Input parameters:                                                |
//|     CD  - Cholesky decomposition of matrix A,                    |
//|           output of SMatrixCholesky subroutine.                  |
//|     N   - size of matrix A.                                      |
//| Result: 1/LowerBound(cond(A))                                    |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::SPDMatrixCholeskyRCond(CMatrixDouble &a,const int n,
                                             const bool isupper)
  {
//--- create a variable
   double v=0;
//--- function call
   SPDMatrixRCondCholeskyInternal(a,n,isupper,false,0,v);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Condition number estimate of a Hermitian positive definite matrix|
//| given by Cholesky decomposition.                                 |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| It should be noted that 1-norm and inf-norm condition numbers of |
//| symmetric matrices are equal, so the algorithm doesn't take into |
//| account the differences between these types of norms.            |
//| Input parameters:                                                |
//|     CD  - Cholesky decomposition of matrix A,                    |
//|           output of SMatrixCholesky subroutine.                  |
//|     N   - size of matrix A.                                      |
//| Result: 1/LowerBound(cond(A))                                    |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::HPDMatrixCholeskyRCond(CMatrixComplex &a,const int n,
                                             const bool isupper)
  {
//--- create a variable
   double v=0;
//--- function call
   HPDMatrixRCondCholeskyInternal(a,n,isupper,false,0,v);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Estimate of the condition number of a matrix given by its LU     |
//| decomposition (1-norm)                                           |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| Input parameters:                                                |
//|     LUA         -   LU decomposition of a matrix in compact form.|
//|                     Output of the CMatrixLU subroutine.          |
//|     N           -   size of matrix A.                            |
//| Result: 1/LowerBound(cond(A))                                    |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::CMatrixLURCond1(CMatrixComplex &lua,const int n)
  {
//--- create a variable
   double v=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- function call
   CMatrixRCondLUInternal(lua,n,true,false,0.0,v);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Estimate of the condition number of a matrix given by its LU     |
//| decomposition (infinity norm).                                   |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| Input parameters:                                                |
//|     LUA     -   LU decomposition of a matrix in compact form.    |
//|                 Output of the CMatrixLU subroutine.              |
//|     N       -   size of matrix A.                                |
//| Result: 1/LowerBound(cond(A))                                    |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::CMatrixLURCondInf(CMatrixComplex &lua,const int n)
  {
//--- create a variable
   double v=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- function call
   CMatrixRCondLUInternal(lua,n,false,false,0.0,v);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Triangular matrix: estimate of a condition number (1-norm)       |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| Input parameters:                                                |
//|     A       -   matrix. Array[0..N-1, 0..N-1].                   |
//|     N       -   size of A.                                       |
//|     IsUpper -   True, if the matrix is upper triangular.         |
//|     IsUnit  -   True, if the matrix has a unit diagonal.         |
//| Result: 1/LowerBound(cond(A))                                    |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::CMatrixTrRCond1(CMatrixComplex &a,const int n,
                                      const bool isupper,const bool isunit)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   double nrm=0;
   int    j1=0;
   int    j2=0;
//--- create arrays
   int    pivots[];
   double t[];
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- allocation
   ArrayResizeAL(t,n);
//--- fiiling array
   for(i=0;i<n;i++)
      t[i]=0;
   for(i=0;i<n;i++)
     {
      //--- check
      if(isupper)
        {
         j1=i+1;
         j2=n-1;
        }
      else
        {
         j1=0;
         j2=i-1;
        }
      //--- change t
      for(j=j1;j<=j2;j++)
         t[j]=t[j]+CMath::AbsComplex(a[i][j]);
      //--- check
      if(isunit)
         t[i]=t[i]+1;
      else
         t[i]=t[i]+CMath::AbsComplex(a[i][i]);
     }
//--- change values
   nrm=0;
   for(i=0;i<n;i++)
      nrm=MathMax(nrm,t[i]);
//--- function call
   CMatrixRCondTrInternal(a,n,isupper,isunit,true,nrm,v);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Triangular matrix: estimate of a matrix condition number         |
//| (infinity-norm).                                                 |
//| The algorithm calculates a lower bound of the condition number.  |
//| In this case, the algorithm does not return a lower bound of the |
//| condition number, but an inverse number (to avoid an overflow in |
//| case of a singular matrix).                                      |
//| Input parameters:                                                |
//|     A   -   matrix. Array whose indexes range within             |
//|             [0..N-1, 0..N-1].                                    |
//|     N   -   size of matrix A.                                    |
//|     IsUpper -   True, if the matrix is upper triangular.         |
//|     IsUnit  -   True, if the matrix has a unit diagonal.         |
//| Result: 1/LowerBound(cond(A))                                    |
//| NOTE:                                                            |
//|     if k(A) is very large, then matrix is assumed degenerate,    |
//|     k(A)=INF, 0.0 is returned in such cases.                     |
//+------------------------------------------------------------------+
static double CRCond::CMatrixTrRCondInf(CMatrixComplex &a,const int n,
                                        const bool isupper,const bool isunit)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   double nrm=0;
   int    j1=0;
   int    j2=0;
//--- create array
   int pivots[];
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- change values
   nrm=0;
   for(i=0;i<n;i++)
     {
      //--- check
      if(isupper)
        {
         j1=i+1;
         j2=n-1;
        }
      else
        {
         j1=0;
         j2=i-1;
        }
      //--- change v
      v=0;
      for(j=j1;j<=j2;j++)
         v=v+CMath::AbsComplex(a[i][j]);
      //--- check
      if(isunit)
         v=v+1;
      else
         v=v+CMath::AbsComplex(a[i][i]);
      nrm=MathMax(nrm,v);
     }
//--- function call
   CMatrixRCondTrInternal(a,n,isupper,isunit,false,nrm,v);
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Threshold for rcond: matrices with condition number beyond this  |
//| threshold are considered singular.                               |
//| Threshold must be far enough from underflow, at least            |
//| Sqr(Threshold)  must be greater than underflow.                  |
//+------------------------------------------------------------------+
static double CRCond::RCondThreshold(void)
  {
//--- return result
   return(MathSqrt(MathSqrt(CMath::m_minrealnumber)));
  }
//+------------------------------------------------------------------+
//| Internal subroutine for condition number estimation              |
//|   -- LAPACK routine (version 3.0)                                |
//|      Univ. of Tennessee,Univ. of California Berkeley,NAG Ltd.,   |
//|      Courant Institute,Argonne National Lab,and Rice University  |
//|      February 29,1992                                            |
//+------------------------------------------------------------------+
static void CRCond::RMatrixRCondTrInternal(CMatrixDouble &a,const int n,
                                           const bool isupper,const bool isunit,
                                           const bool onenorm,double anorm,
                                           double &rc)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    kase=0;
   int    kase1=0;
   int    j1=0;
   int    j2=0;
   double ainvnm=0;
   double maxgrowth=0;
   double s=0;
   bool   mupper;
   bool   mtrans;
   bool   munit;
//--- create arrays
   double ex[];
   double ev[];
   int iwork[];
   double tmp[];
//--- check
   if(onenorm)
      kase1=1;
   else
      kase1=2;
//--- initialization
   rc=0;
   mupper=true;
   mtrans=true;
   munit=true;
//--- allocation
   ArrayResizeAL(iwork,n+1);
   ArrayResizeAL(tmp,n);
//--- prepare parameters for triangular solver
   maxgrowth=1/RCondThreshold();
   s=0;
   for(i=0;i<n;i++)
     {
      //--- check
      if(isupper)
        {
         j1=i+1;
         j2=n-1;
        }
      else
        {
         j1=0;
         j2=i-1;
        }
      //--- change s
      for(j=j1;j<=j2;j++)
         s=MathMax(s,MathAbs(a[i][j]));
      //--- check
      if(isunit)
         s=MathMax(s,1);
      else
         s=MathMax(s,MathAbs(a[i][i]));
     }
//--- check
   if(s==0.0)
      s=1;
   s=1/s;
//--- Scale according to S
   anorm=anorm*s;
//--- Quick return if possible
//--- We assume that ANORM<>0 after this block
   if(anorm==0.0)
      return;
//--- check
   if(n==1)
     {
      rc=1;
      return;
     }
//--- Estimate the norm of inv(A).
   ainvnm=0;
   kase=0;
   while(true)
     {
      //--- function call
      RMatrixEstimateNorm(n,ev,ex,iwork,ainvnm,kase);
      //--- check
      if(kase==0)
         break;
      //--- from 1-based array to 0-based
      for(i=0;i<n;i++)
         ex[i]=ex[i+1];
      //--- multiply by inv(A) or inv(A')
      if(kase==kase1)
        {
         //--- multiply by inv(A)
         if(!CSafeSolve::RMatrixScaledTrSafeSolve(a,s,n,ex,isupper,0,isunit,maxgrowth))
            return;
        }
      else
        {
         //--- multiply by inv(A')
         if(!CSafeSolve::RMatrixScaledTrSafeSolve(a,s,n,ex,isupper,1,isunit,maxgrowth))
            return;
        }
      //--- from 0-based array to 1-based
      for(i=n-1;i>=0;i--)
         ex[i+1]=ex[i];
     }
//--- Compute the estimate of the reciprocal condition number.
   if(ainvnm!=0.0)
     {
      rc=1/ainvnm;
      rc=rc/anorm;
      //--- check
      if(rc<RCondThreshold())
         rc=0;
     }
  }
//+------------------------------------------------------------------+
//| Condition number estimation                                      |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      March 31, 1993                                              |
//+------------------------------------------------------------------+
static void CRCond::CMatrixRCondTrInternal(CMatrixComplex &a,const int n,
                                           const bool isupper,const bool isunit,
                                           const bool onenorm,double anorm,
                                           double &rc)
  {
//--- create variables
   int    kase=0;
   int    kase1=0;
   double ainvnm=0;
   int    i=0;
   int    j=0;
   int    j1=0;
   int    j2=0;
   double s=0;
   double maxgrowth=0;
//--- create arrays
   al_complex ex[];
   al_complex cwork2[];
   al_complex cwork3[];
   al_complex cwork4[];
   int        isave[];
   double     rsave[];
//--- initialization
   rc=0;
//--- check
   if(n<=0)
      return;
//--- check
   if(n==0)
     {
      rc=1;
      return;
     }
//--- allocation
   ArrayResizeAL(cwork2,n+1);
//--- prepare parameters for triangular solver
   maxgrowth=1/RCondThreshold();
   s=0;
   for(i=0;i<n;i++)
     {
      //--- check
      if(isupper)
        {
         j1=i+1;
         j2=n-1;
        }
      else
        {
         j1=0;
         j2=i-1;
        }
      //--- change s
      for(j=j1;j<=j2;j++)
         s=MathMax(s,CMath::AbsComplex(a[i][j]));
      //--- check
      if(isunit)
         s=MathMax(s,1);
      else
         s=MathMax(s,CMath::AbsComplex(a[i][i]));
     }
//--- check
   if(s==0.0)
      s=1;
   s=1/s;
//--- Scale according to S
   anorm=anorm*s;
   if(anorm==0.0)
      return;
//--- Estimate the norm of inv(A).
   ainvnm=0;
//--- check
   if(onenorm)
      kase1=1;
   else
      kase1=2;
   kase=0;
   while(true)
     {
      //--- function call
      CMatrixEstimateNorm(n,cwork4,ex,ainvnm,kase,isave,rsave);
      //--- check
      if(kase==0)
         break;
      //--- from 1-based array to 0-based
      for(i=0;i<n;i++)
         ex[i]=ex[i+1];
      //--- multiply by inv(A) or inv(A')
      if(kase==kase1)
        {
         //--- multiply by inv(A)
         if(!CSafeSolve::CMatrixScaledTrSafeSolve(a,s,n,ex,isupper,0,isunit,maxgrowth))
            return;
        }
      else
        {
         //--- multiply by inv(A')
         if(!CSafeSolve::CMatrixScaledTrSafeSolve(a,s,n,ex,isupper,2,isunit,maxgrowth))
            return;
        }
      //--- from 0-based array to 1-based
      for(i=n-1;i>=0;i--)
         ex[i+1]=ex[i];
     }
//--- Compute the estimate of the reciprocal condition number.
   if(ainvnm!=0.0)
     {
      rc=1/ainvnm;
      rc=rc/anorm;
      //--- check
      if(rc<RCondThreshold())
         rc=0;
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine for condition number estimation              |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      February 29, 1992                                           |
//+------------------------------------------------------------------+
static void CRCond::SPDMatrixRCondCholeskyInternal(CMatrixDouble &cha,const int n,
                                                   const bool isupper,
                                                   const bool isnormprovided,
                                                   double anorm,double &rc)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    kase=0;
   double ainvnm=0;
   double sa=0;
   double v=0;
   double maxgrowth=0;
   int    i_=0;
   int    i1_=0;
//--- create arrays
   double ex[];
   double ev[];
   double tmp[];
   int    iwork[];
//--- check
   if(!CAp::Assert(n>=1))
      return;
//--- allocation
   ArrayResizeAL(tmp,n);
//--- initialization
   rc=0;
//--- prepare parameters for triangular solver
   maxgrowth=1/RCondThreshold();
   sa=0;
//--- check
   if(isupper)
     {
      for(i=0;i<n;i++)
         for(j=i;j<n;j++)
            sa=MathMax(sa,CMath::AbsComplex(cha[i][j]));
     }
   else
     {
      for(i=0;i<n;i++)
         for(j=0;j<=i;j++)
            sa=MathMax(sa,CMath::AbsComplex(cha[i][j]));
     }
//--- check
   if(sa==0.0)
      sa=1;
   sa=1/sa;
//--- Estimate the norm of A
   if(!isnormprovided)
     {
      kase=0;
      anorm=0;
      while(true)
        {
         //--- function call
         RMatrixEstimateNorm(n,ev,ex,iwork,anorm,kase);
         //--- check
         if(kase==0)
            break;
         //--- check
         if(isupper)
           {
            //--- Multiply by U
            for(i=1;i<=n;i++)
              {
               i1_=1;
               v=0.0;
               for(i_=i-1;i_<n;i_++)
                  v+=cha[i-1][i_]*ex[i_+i1_];
               ex[i]=v;
              }
            for(i_=1;i_<=n;i_++)
               ex[i_]=sa*ex[i_];
            //--- Multiply by U'
            for(i=0;i<n;i++)
               tmp[i]=0;
            for(i=0;i<n;i++)
              {
               v=ex[i+1];
               for(i_=i;i_<n;i_++)
                  tmp[i_]=tmp[i_]+v*cha[i][i_];
              }
            //--- change values
            i1_=-1;
            for(i_=1;i_<=n;i_++)
               ex[i_]=tmp[i_+i1_];
            for(i_=1;i_<=n;i_++)
               ex[i_]=sa*ex[i_];
           }
         else
           {
            //--- Multiply by L''
            for(i=0;i<n;i++)
               tmp[i]=0;
            for(i=0;i<n;i++)
              {
               v=ex[i+1];
               for(i_=0;i_<=i;i_++)
                  tmp[i_]=tmp[i_]+v*cha[i][i_];
              }
            //--- change values
            i1_=-1;
            for(i_=1;i_<=n;i_++)
               ex[i_]=tmp[i_+i1_];
            for(i_=1;i_<=n;i_++)
               ex[i_]=sa*ex[i_];
            //--- Multiply by L'
            for(i=n;i>=1;i--)
              {
               i1_=1;
               v=0.0;
               for(i_=0;i_<i;i_++)
                  v+=cha[i-1][i_]*ex[i_+i1_];
               ex[i]=v;
              }
            for(i_=1;i_<=n;i_++)
               ex[i_]=sa*ex[i_];
           }
        }
     }
//--- check
   if(anorm==0.0)
      return;
//--- check
   if(n==1)
     {
      rc=1;
      return;
     }
//--- Estimate the 1-norm of inv(A).
   kase=0;
   while(true)
     {
      //--- function call
      RMatrixEstimateNorm(n,ev,ex,iwork,ainvnm,kase);
      //--- check
      if(kase==0)
         break;
      for(i=0;i<n;i++)
         ex[i]=ex[i+1];
      //--- check
      if(isupper)
        {
         //--- Multiply by inv(U')
         if(!CSafeSolve::RMatrixScaledTrSafeSolve(cha,sa,n,ex,isupper,1,false,maxgrowth))
            return;
         //--- Multiply by inv(U)
         if(!CSafeSolve::RMatrixScaledTrSafeSolve(cha,sa,n,ex,isupper,0,false,maxgrowth))
            return;
        }
      else
        {
         //--- Multiply by inv(L)
         if(!CSafeSolve::RMatrixScaledTrSafeSolve(cha,sa,n,ex,isupper,0,false,maxgrowth))
            return;
         //--- Multiply by inv(L')
         if(!CSafeSolve::RMatrixScaledTrSafeSolve(cha,sa,n,ex,isupper,1,false,maxgrowth))
            return;
        }
      for(i=n-1;i>=0;i--)
         ex[i+1]=ex[i];
     }
//--- Compute the estimate of the reciprocal condition number.
   if(ainvnm!=0.0)
     {
      v=1/ainvnm;
      rc=v/anorm;
      //--- check
      if(rc<RCondThreshold())
         rc=0;
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine for condition number estimation              |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      February 29, 1992                                           |
//+------------------------------------------------------------------+
static void CRCond::HPDMatrixRCondCholeskyInternal(CMatrixComplex &cha,const int n,
                                                   const bool isupper,
                                                   const bool isnormprovided,
                                                   double anorm,double &rc)
  {
//--- create variables
   al_complex Csa;
   int        kase=0;
   double     ainvnm=0;
   al_complex v=0;
   int        i=0;
   int        j=0;
   double     sa=0;
   double     maxgrowth=0;
   int        i_=0;
   int        i1_=0;
//--- create arrays
   int        isave[];
   double     rsave[];
   al_complex ex[];
   al_complex ev[];
   al_complex tmp[];
//--- check
   if(!CAp::Assert(n>=1))
      return;
//--- allocation
   ArrayResizeAL(tmp,n);
//--- initialization
   rc=0;
//--- prepare parameters for triangular solver
   maxgrowth=1/RCondThreshold();
   sa=0;
//--- check
   if(isupper)
     {
      for(i=0;i<n;i++)
         for(j=i;j<n;j++)
            sa=MathMax(sa,CMath::AbsComplex(cha[i][j]));
     }
   else
     {
      for(i=0;i<n;i++)
         for(j=0;j<=i;j++)
            sa=MathMax(sa,CMath::AbsComplex(cha[i][j]));
     }
//--- check
   if(sa==0.0)
      sa=1;
   sa=1/sa;
//--- Estimate the norm of A
   if(!isnormprovided)
     {
      anorm=0;
      kase=0;
      while(true)
        {
         //--- function call
         CMatrixEstimateNorm(n,ev,ex,anorm,kase,isave,rsave);
         //--- check
         if(kase==0)
            break;
         //--- check
         if(isupper)
           {
            //--- Multiply by U
            for(i=1;i<=n;i++)
              {
               i1_=1;
               v=0.0;
               for(i_=i-1;i_<n;i_++)
                  v+=cha[i-1][i_]*ex[i_+i1_];
               ex[i]=v;
              }
            for(i_=1;i_<=n;i_++)
              {
               Csa=sa;
               ex[i_]=Csa*ex[i_];
              }
            //--- Multiply by U'
            for(i=0;i<n;i++)
               tmp[i]=0;
            for(i=0;i<n;i++)
              {
               v=ex[i+1];
               for(i_=i;i_<n;i_++)
                  tmp[i_]=tmp[i_]+v*CMath::Conj(cha[i][i_]);
              }
            //--- change values
            i1_=-1;
            for(i_=1;i_<=n;i_++)
               ex[i_]=tmp[i_+i1_];
            for(i_=1;i_<=n;i_++)
              {
               Csa=sa;
               ex[i_]=Csa*ex[i_];
              }
           }
         else
           {
            //--- Multiply by L''
            for(i=0;i<n;i++)
               tmp[i]=0;
            for(i=0;i<n;i++)
              {
               v=ex[i+1];
               for(i_=0;i_<=i;i_++)
                  tmp[i_]=tmp[i_]+v*CMath::Conj(cha[i][i_]);
              }
            //--- change values
            i1_=-1;
            for(i_=1;i_<=n;i_++)
               ex[i_]=tmp[i_+i1_];
            for(i_=1;i_<=n;i_++)
              {
               Csa=sa;
               ex[i_]=Csa*ex[i_];
              }
            //--- Multiply by L'
            for(i=n;i>=1;i--)
              {
               i1_=1;
               v=0.0;
               for(i_=0;i_<i;i_++)
                  v+=cha[i-1][i_]*ex[i_+i1_];
               ex[i]=v;
              }
            for(i_=1;i_<=n;i_++)
              {
               Csa=sa;
               ex[i_]=Csa*ex[i_];
              }
           }
        }
     }
//--- Quick return if possible
//--- After this block we assume that ANORM<>0
   if(anorm==0.0)
      return;
//--- check
   if(n==1)
     {
      rc=1;
      return;
     }
//--- Estimate the norm of inv(A).
   ainvnm=0;
   kase=0;
   while(true)
     {
      //--- function call
      CMatrixEstimateNorm(n,ev,ex,ainvnm,kase,isave,rsave);
      //--- check
      if(kase==0)
         break;
      for(i=0;i<n;i++)
         ex[i]=ex[i+1];
      //--- check
      if(isupper)
        {
         //--- Multiply by inv(U')
         if(!CSafeSolve::CMatrixScaledTrSafeSolve(cha,sa,n,ex,isupper,2,false,maxgrowth))
            return;
         //--- Multiply by inv(U)
         if(!CSafeSolve::CMatrixScaledTrSafeSolve(cha,sa,n,ex,isupper,0,false,maxgrowth))
            return;
        }
      else
        {
         //--- Multiply by inv(L)
         if(!CSafeSolve::CMatrixScaledTrSafeSolve(cha,sa,n,ex,isupper,0,false,maxgrowth))
            return;
         //--- Multiply by inv(L')
         if(!CSafeSolve::CMatrixScaledTrSafeSolve(cha,sa,n,ex,isupper,2,false,maxgrowth))
            return;
        }
      for(i=n-1;i>=0;i--)
         ex[i+1]=ex[i];
     }
//--- Compute the estimate of the reciprocal condition number.
   if(ainvnm!=0.0)
     {
      rc=1/ainvnm;
      rc=rc/anorm;
      //--- check
      if(rc<RCondThreshold())
         rc=0;
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine for condition number estimation              |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      February 29, 1992                                           |
//+------------------------------------------------------------------+
static void CRCond::RMatrixRCondLUInternal(CMatrixDouble &lua,const int n,
                                           const bool onenorm,
                                           const bool isanormprovided,
                                           double anorm,double &rc)
  {
//--- create variables
   double v=0;
   int    i=0;
   int    j=0;
   int    kase=0;
   int    kase1=0;
   double ainvnm=0;
   double maxgrowth=0;
   double su=0;
   double sl=0;
   bool   mupper;
   bool   mtrans;
   bool   munit;
   int    i_=0;
   int    i1_=0;
//--- create arrays
   double ex[];
   double ev[];
   int    iwork[];
   double tmp[];
//--- check
   if(onenorm)
      kase1=1;
   else
      kase1=2;
//--- initialization
   rc=0;
   mupper=true;
   mtrans=true;
   munit=true;
//--- allocation
   ArrayResizeAL(iwork,n+1);
   ArrayResizeAL(tmp,n);
//--- prepare parameters for triangular solver
   maxgrowth=1/RCondThreshold();
   su=0;
   sl=1;
   for(i=0;i<n;i++)
     {
      for(j=0;j<i;j++)
         sl=MathMax(sl,MathAbs(lua[i][j]));
      for(j=i;j<n;j++)
         su=MathMax(su,MathAbs(lua[i][j]));
     }
//--- check
   if(su==0.0)
      su=1;
   su=1/su;
   sl=1/sl;
//--- Estimate the norm of A
   if(!isanormprovided)
     {
      kase=0;
      anorm=0;
      while(true)
        {
         //--- function call
         RMatrixEstimateNorm(n,ev,ex,iwork,anorm,kase);
         //--- check
         if(kase==0)
            break;
         //--- check
         if(kase==kase1)
           {
            //--- Multiply by U
            for(i=1;i<=n;i++)
              {
               i1_=1;
               v=0.0;
               for(i_=i-1;i_<n;i_++)
                  v+=lua[i-1][i_]*ex[i_+i1_];
               ex[i]=v;
              }
            //--- Multiply by L
            for(i=n;i>=1;i--)
              {
               //--- check
               if(i>1)
                 {
                  i1_=1;
                  v=0.0;
                  for(i_=0;i_<=i-2;i_++)
                     v+=lua[i-1][i_]*ex[i_+i1_];
                 }
               else
                  v=0;
               ex[i]=ex[i]+v;
              }
           }
         else
           {
            //--- Multiply by L'
            for(i=0;i<n;i++)
               tmp[i]=0;
            for(i=0;i<n;i++)
              {
               v=ex[i+1];
               //--- check
               if(i>=1)
                 {
                  for(i_=0;i_<i;i_++)
                     tmp[i_]=tmp[i_]+v*lua[i][i_];
                 }
               tmp[i]=tmp[i]+v;
              }
            //--- change values
            i1_=-1;
            for(i_=1;i_<=n;i_++)
               ex[i_]=tmp[i_+i1_];
            //--- Multiply by U'
            for(i=0;i<n;i++)
               tmp[i]=0;
            for(i=0;i<n;i++)
              {
               v=ex[i+1];
               for(i_=i;i_<n;i_++)
                  tmp[i_]=tmp[i_]+v*lua[i][i_];
              }
            //--- change values
            i1_=-1;
            for(i_=1;i_<=n;i_++)
               ex[i_]=tmp[i_+i1_];
           }
        }
     }
//--- Scale according to SU/SL
   anorm=anorm*su*sl;
//--- Quick return if possible
//--- We assume that ANORM<>0 after this block
   if(anorm==0.0)
      return;
//--- check
   if(n==1)
     {
      rc=1;
      return;
     }
//--- Estimate the norm of inv(A).
   ainvnm=0;
   kase=0;
   while(true)
     {
      //--- function call
      RMatrixEstimateNorm(n,ev,ex,iwork,ainvnm,kase);
      //--- check
      if(kase==0)
         break;
      //--- from 1-based array to 0-based
      for(i=0;i<n;i++)
         ex[i]=ex[i+1];
      //--- multiply by inv(A) or inv(A')
      if(kase==kase1)
        {
         //--- Multiply by inv(L)
         if(!CSafeSolve::RMatrixScaledTrSafeSolve(lua,sl,n,ex,!mupper,0,munit,maxgrowth))
            return;
         //--- Multiply by inv(U)
         if(!CSafeSolve::RMatrixScaledTrSafeSolve(lua,su,n,ex,mupper,0,!munit,maxgrowth))
            return;
        }
      else
        {
         //--- Multiply by inv(U')
         if(!CSafeSolve::RMatrixScaledTrSafeSolve(lua,su,n,ex,mupper,1,!munit,maxgrowth))
            return;
         //--- Multiply by inv(L')
         if(!CSafeSolve::RMatrixScaledTrSafeSolve(lua,sl,n,ex,!mupper,1,munit,maxgrowth))
            return;
        }
      //--- from 0-based array to 1-based
      for(i=n-1;i>=0;i--)
         ex[i+1]=ex[i];
     }
//--- Compute the estimate of the reciprocal condition number.
   if(ainvnm!=0.0)
     {
      rc=1/ainvnm;
      rc=rc/anorm;
      //--- check
      if(rc<RCondThreshold())
         rc=0;
     }
  }
//+------------------------------------------------------------------+
//| Condition number estimation                                      |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      March 31, 1993                                              |
//+------------------------------------------------------------------+
static void CRCond::CMatrixRCondLUInternal(CMatrixComplex &lua,const int n,
                                           const bool onenorm,
                                           const bool isanormprovided,
                                           double anorm,double &rc)
  {
//--- create variables
   int        kase=0;
   int        kase1=0;
   double     ainvnm=0;
   al_complex v=0;
   int        i=0;
   int        j=0;
   double     su=0;
   double     sl=0;
   double     maxgrowth=0;
   int        i_=0;
   int        i1_=0;
//--- create arrays
   al_complex ex[];
   al_complex cwork2[];
   al_complex cwork3[];
   al_complex cwork4[];
   int        isave[];
   double     rsave[];
//--- check
   if(n<=0)
      return;
//--- allocation
   ArrayResizeAL(cwork2,n+1);
//--- initialization
   rc=0;
//--- check
   if(n==0)
     {
      rc=1;
      return;
     }
//--- prepare parameters for triangular solver
   maxgrowth=1/RCondThreshold();
   su=0;
   sl=1;
   for(i=0;i<n;i++)
     {
      for(j=0;j<i;j++)
         sl=MathMax(sl,CMath::AbsComplex(lua[i][j]));
      for(j=i;j<n;j++)
         su=MathMax(su,CMath::AbsComplex(lua[i][j]));
     }
//--- check
   if(su==0.0)
      su=1;
   su=1/su;
   sl=1/sl;
//--- Estimate the norm of SU*SL*A
   if(!isanormprovided)
     {
      anorm=0;
      //--- check
      if(onenorm)
         kase1=1;
      else
         kase1=2;
      kase=0;
      do
        {
         //--- function call
         CMatrixEstimateNorm(n,cwork4,ex,anorm,kase,isave,rsave);
         //--- check
         if(kase!=0)
           {
            //--- check
            if(kase==kase1)
              {
               //--- Multiply by U
               for(i=1;i<=n;i++)
                 {
                  i1_=1;
                  v=0.0;
                  for(i_=i-1;i_<n;i_++)
                     v+=lua[i-1][i_]*ex[i_+i1_];
                  ex[i]=v;
                 }
               //--- Multiply by L
               for(i=n;i>=1;i--)
                 {
                  v=0;
                  //--- check
                  if(i>1)
                    {
                     i1_=1;
                     v=0.0;
                     for(i_=0;i_<=i-2;i_++)
                        v+=lua[i-1][i_]*ex[i_+i1_];
                    }
                  ex[i]=v+ex[i];
                 }
              }
            else
              {
               //--- Multiply by L'
               for(i=1;i<=n;i++)
                  cwork2[i]=0;
               for(i=1;i<=n;i++)
                 {
                  v=ex[i];
                  //--- check
                  if(i>1)
                    {
                     i1_=-1;
                     for(i_=1;i_<i;i_++)
                        cwork2[i_]=cwork2[i_]+v*CMath::Conj(lua[i-1][i_+i1_]);
                    }
                  cwork2[i]=cwork2[i]+v;
                 }
               //--- Multiply by U'
               for(i=1;i<=n;i++)
                  ex[i]=0;
               for(i=1;i<=n;i++)
                 {
                  v=cwork2[i];
                  i1_=-1;
                  for(i_=i;i_<=n;i_++)
                     ex[i_]=ex[i_]+v*CMath::Conj(lua[i-1][i_+i1_]);
                 }
              }
           }
        }
      while(kase!=0);
     }
//--- Scale according to SU/SL
   anorm=anorm*su*sl;
//--- check
   if(anorm==0.0)
      return;
//--- Estimate the norm of inv(A).
   ainvnm=0;
//--- check
   if(onenorm)
      kase1=1;
   else
      kase1=2;
   kase=0;
   while(true)
     {
      //--- function call
      CMatrixEstimateNorm(n,cwork4,ex,ainvnm,kase,isave,rsave);
      //--- check
      if(kase==0)
         break;
      //--- from 1-based array to 0-based
      for(i=0;i<n;i++)
         ex[i]=ex[i+1];
      //--- multiply by inv(A) or inv(A')
      if(kase==kase1)
        {
         //--- Multiply by inv(L)
         if(!CSafeSolve::CMatrixScaledTrSafeSolve(lua,sl,n,ex,false,0,true,maxgrowth))
           {
            rc=0;
            //--- exit the function
            return;
           }
         //--- Multiply by inv(U)
         if(!CSafeSolve::CMatrixScaledTrSafeSolve(lua,su,n,ex,true,0,false,maxgrowth))
           {
            rc=0;
            //--- exit the function
            return;
           }
        }
      else
        {
         //--- Multiply by inv(U')
         if(!CSafeSolve::CMatrixScaledTrSafeSolve(lua,su,n,ex,true,2,false,maxgrowth))
           {
            rc=0;
            //--- exit the function
            return;
           }
         //--- Multiply by inv(L')
         if(!CSafeSolve::CMatrixScaledTrSafeSolve(lua,sl,n,ex,false,2,true,maxgrowth))
           {
            rc=0;
            //--- exit the function
            return;
           }
        }
      //--- from 0-based array to 1-based
      for(i=n-1;i>=0;i--)
         ex[i+1]=ex[i];
     }
//--- Compute the estimate of the reciprocal condition number.
   if(ainvnm!=0.0)
     {
      rc=1/ainvnm;
      rc=rc/anorm;
      //--- check
      if(rc<RCondThreshold())
         rc=0;
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine for matrix norm estimation                   |
//|   -- LAPACK auxiliary routine (version 3.0) --                   |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      February 29, 1992                                           |
//+------------------------------------------------------------------+
static void CRCond::RMatrixEstimateNorm(const int n,double &v[],double &x[],
                                        int &isgn[],double &est,int &kase)
  {
//--- create variables
   int    itmax=0;
   int    i=0;
   double t=0;
   bool   flg;
   int    positer=0;
   int    posj=0;
   int    posjlast=0;
   int    posjump=0;
   int    posaltsgn=0;
   int    posestold=0;
   int    postemp=0;
   int    i_=0;
//--- initialization
   itmax=5;
   posaltsgn=n+1;
   posestold=n+2;
   postemp=n+3;
   positer=n+1;
   posj=n+2;
   posjlast=n+3;
   posjump=n+4;
//--- check
   if(kase==0)
     {
      //--- allocation
      ArrayResizeAL(v,n+4);
      ArrayResizeAL(x,n+1);
      ArrayResizeAL(isgn,n+5);
      //--- change values
      t=1.0/(double)n;
      for(i=1;i<=n;i++)
         x[i]=t;
      kase=1;
      isgn[posjump]=1;
      //--- exit the function
      return;
     }
//--- ................ entry   (jump = 1)
//--- first iteration.  x has been overwritten by a*x.
   if(isgn[posjump]==1)
     {
      //--- check
      if(n==1)
        {
         v[1]=x[1];
         est=MathAbs(v[1]);
         kase=0;
         //--- exit the function
         return;
        }
      //--- change value
      est=0;
      for(i=1;i<=n;i++)
         est=est+MathAbs(x[i]);
      for(i=1;i<=n;i++)
        {
         //--- check
         if(x[i]>=0.0)
            x[i]=1;
         else
            x[i]=-1;
         //--- check
         if(x[i]>0)
            isgn[i]=1;
         //--- check
         if(x[i]<0)
            isgn[i]=-1;
         //--- check
         if(x[i]==0)
            isgn[i]=0;
        }
      kase=2;
      isgn[posjump]=2;
      //--- exit the function
      return;
     }
//--- ................ entry   (jump = 2)
//--- first iteration.  x has been overwritten by trandpose(a)*x.
   if(isgn[posjump]==2)
     {
      isgn[posj]=1;
      for(i=2;i<=n;i++)
        {
         //--- check
         if(MathAbs(x[i])>MathAbs(x[isgn[posj]]))
            isgn[posj]=i;
        }
      isgn[positer]=2;
      //--- main loop - iterations 2,3,...,itmax.
      for(i=1;i<=n;i++)
         x[i]=0;
      x[isgn[posj]]=1;
      kase=1;
      isgn[posjump]=3;
      //--- exit the function
      return;
     }
//--- ................ entry   (jump = 3)
//--- x has been overwritten by a*x.
   if(isgn[posjump]==3)
     {
      for(i_=1;i_<=n;i_++)
         v[i_]=x[i_];
      v[posestold]=est;
      //--- change value
      est=0;
      for(i=1;i<=n;i++)
         est=est+MathAbs(v[i]);
      flg=false;
      for(i=1;i<=n;i++)
        {
         if(((x[i]>=0.0)&&(isgn[i]<0))||((x[i]<0.0)&&(isgn[i]>=0)))
            flg=true;
        }
      //--- repeated sign vector detected, hence algorithm has converged.
      //--- or may be cycling.
      if(!flg || est<=v[posestold])
        {
         v[posaltsgn]=1;
         for(i=1;i<=n;i++)
           {
            x[i]=v[posaltsgn]*(1+(double)(i-1)/(double)(n-1));
            v[posaltsgn]=-v[posaltsgn];
           }
         kase=1;
         isgn[posjump]=5;
         //--- exit the function
         return;
        }
      for(i=1;i<=n;i++)
        {
         //--- check
         if(x[i]>=0.0)
           {
            x[i]=1;
            isgn[i]=1;
           }
         else
           {
            x[i]=-1;
            isgn[i]=-1;
           }
        }
      kase=2;
      isgn[posjump]=4;
      //--- exit the function
      return;
     }
//--- ................ entry   (jump = 4)
//--- x has been overwritten by trandpose(a)*x.
   if(isgn[posjump]==4)
     {
      isgn[posjlast]=isgn[posj];
      isgn[posj]=1;
      for(i=2;i<=n;i++)
        {
         //--- check
         if(MathAbs(x[i])>MathAbs(x[isgn[posj]]))
            isgn[posj]=i;
        }
      //--- check
      if(x[isgn[posjlast]]!=MathAbs(x[isgn[posj]]) && isgn[positer]<itmax)
        {
         isgn[positer]=isgn[positer]+1;
         for(i=1;i<=n;i++)
            x[i]=0;
         x[isgn[posj]]=1;
         kase=1;
         isgn[posjump]=3;
         //--- exit the function
         return;
        }
      //--- iteration complete.  final stage.
      v[posaltsgn]=1;
      for(i=1;i<=n;i++)
        {
         x[i]=v[posaltsgn]*(1+(double)(i-1)/(double)(n-1));
         v[posaltsgn]=-v[posaltsgn];
        }
      kase=1;
      isgn[posjump]=5;
      //--- exit the function
      return;
     }
//--- ................ entry   (jump = 5)
//--- x has been overwritten by a*x.
   if(isgn[posjump]==5)
     {
      v[postemp]=0;
      for(i=1;i<=n;i++)
         v[postemp]=v[postemp]+MathAbs(x[i]);
      v[postemp]=2*v[postemp]/(3*n);
      //--- check
      if(v[postemp]>est)
        {
         for(i_=1;i_<=n;i_++)
            v[i_]=x[i_];
         est=v[postemp];
        }
      kase=0;
      //--- exit the function
      return;
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CRCond::CMatrixEstimateNorm(const int n,al_complex &v[],al_complex &x[],
                                        double &est,int &kase,int &isave[],
                                        double &rsave[])
  {
//--- create variables
   int    itmax=0;
   int    i=0;
   int    iter=0;
   int    j=0;
   int    jlast=0;
   int    jump=0;
   double absxi=0;
   double altsgn=0;
   double estold=0;
   double safmin=0;
   double temp=0;
   int    i_=0;
//--- Executable Statements
   itmax=5;
   safmin=CMath::m_minrealnumber;
//--- check
   if(kase==0)
     {
      //--- allocation
      ArrayResizeAL(v,n+1);
      ArrayResizeAL(x,n+1);
      ArrayResizeAL(isave,5);
      ArrayResizeAL(rsave,4);
      for(i=1;i<=n;i++)
         x[i]=1.0/(double)n;
      kase=1;
      jump=1;
      //--- function call
      InternalComplexRCondSaveAll(isave,rsave,i,iter,j,jlast,jump,absxi,altsgn,estold,temp);
      //--- exit the function
      return;
     }
//--- function call
   InternalComplexRCondLoadAll(isave,rsave,i,iter,j,jlast,jump,absxi,altsgn,estold,temp);
//--- entry   (jump = 1)
//--- first iteration.  x has been overwritten by a*x.
   if(jump==1)
     {
      //--- check
      if(n==1)
        {
         v[1]=x[1];
         //--- function call
         est=CMath::AbsComplex(v[1]);
         kase=0;
         //--- function call
         InternalComplexRCondSaveAll(isave,rsave,i,iter,j,jlast,jump,absxi,altsgn,estold,temp);
         //--- exit the function
         return;
        }
      //--- function call
      est=InternalComplexRCondScSum1(x,n);
      for(i=1;i<=n;i++)
        {
         //--- function call
         absxi=CMath::AbsComplex(x[i]);
         //--- check
         if(absxi>safmin)
            x[i]=x[i]/absxi;
         else
            x[i]=1;
        }
      kase=2;
      jump=2;
      //--- function call
      InternalComplexRCondSaveAll(isave,rsave,i,iter,j,jlast,jump,absxi,altsgn,estold,temp);
      //--- exit the function
      return;
     }
//--- entry   (jump = 2)
//--- first iteration.  x has been overwritten by ctrans(a)*x.
   if(jump==2)
     {
      j=InternalComplexRCondIcMax1(x,n);
      iter=2;
      //--- main loop - iterations 2,3,...,itmax.
      for(i=1;i<=n;i++)
         x[i]=0;
      x[j]=1;
      kase=1;
      jump=3;
      //--- function call
      InternalComplexRCondSaveAll(isave,rsave,i,iter,j,jlast,jump,absxi,altsgn,estold,temp);
      //--- exit the function
      return;
     }
//--- entry   (jump = 3)
//--- x has been overwritten by a*x.
   if(jump==3)
     {
      for(i_=1;i_<=n;i_++)
         v[i_]=x[i_];
      estold=est;
      //--- function call
      est=InternalComplexRCondScSum1(v,n);
      //--- test for cycling.
      if(est<=estold)
        {
         //--- iteration complete.  final stage.
         altsgn=1;
         for(i=1;i<=n;i++)
           {
            x[i]=altsgn*(1+(double)(i-1)/(double)(n-1));
            altsgn=-altsgn;
           }
         kase=1;
         jump=5;
         //--- function call
         InternalComplexRCondSaveAll(isave,rsave,i,iter,j,jlast,jump,absxi,altsgn,estold,temp);
         //--- exit the function
         return;
        }
      for(i=1;i<=n;i++)
        {
         absxi=CMath::AbsComplex(x[i]);
         //--- check
         if(absxi>safmin)
            x[i]=x[i]/absxi;
         else
            x[i]=1;
        }
      kase=2;
      jump=4;
      //--- function call
      InternalComplexRCondSaveAll(isave,rsave,i,iter,j,jlast,jump,absxi,altsgn,estold,temp);
      //--- exit the function
      return;
     }
//--- entry   (jump = 4)
//--- x has been overwritten by ctrans(a)*x.
   if(jump==4)
     {
      jlast=j;
      j=InternalComplexRCondIcMax1(x,n);
      //--- check
      if(CMath::AbsComplex(x[jlast])!=CMath::AbsComplex(x[j]) && iter<itmax)
        {
         iter=iter+1;
         //--- main loop - iterations 2,3,...,itmax.
         for(i=1;i<=n;i++)
            x[i]=0;
         x[j]=1;
         kase=1;
         jump=3;
         //--- function call
         InternalComplexRCondSaveAll(isave,rsave,i,iter,j,jlast,jump,absxi,altsgn,estold,temp);
         //--- exit the function
         return;
        }
      //--- iteration complete.  final stage.
      altsgn=1;
      for(i=1;i<=n;i++)
        {
         x[i]=altsgn*(1+(double)(i-1)/(double)(n-1));
         altsgn=-altsgn;
        }
      kase=1;
      jump=5;
      //--- function call
      InternalComplexRCondSaveAll(isave,rsave,i,iter,j,jlast,jump,absxi,altsgn,estold,temp);
      //--- exit the function
      return;
     }
//--- entry   (jump = 5)
//--- x has been overwritten by a*x.
   if(jump==5)
     {
      temp=2*(InternalComplexRCondScSum1(x,n)/(3*n));
      //--- check
      if(temp>est)
        {
         for(i_=1;i_<=n;i_++)
            v[i_]=x[i_];
         est=temp;
        }
      kase=0;
      //--- function call
      InternalComplexRCondSaveAll(isave,rsave,i,iter,j,jlast,jump,absxi,altsgn,estold,temp);
      //--- exit the function
      return;
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static double CRCond::InternalComplexRCondScSum1(al_complex &x[],const int n)
  {
//--- create variables
   double result=0;
   int    i=0;
//--- get result
   for(i=1;i<=n;i++)
      result=result+CMath::AbsComplex(x[i]);
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static int CRCond::InternalComplexRCondIcMax1(al_complex &x[],const int n)
  {
//--- create variables
   int    result=0;
   int    i=0;
   double m=0;
//--- get result
   result=1;
   m=CMath::AbsComplex(x[1]);
   for(i=2;i<=n;i++)
     {
      //--- check
      if(CMath::AbsComplex(x[i])>m)
        {
         result=i;
         m=CMath::AbsComplex(x[i]);
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CRCond::InternalComplexRCondSaveAll(int &isave[],double &rsave[],
                                                int &i,int &iter,int &j,
                                                int &jlast,int &jump,
                                                double &absxi,double &altsgn,
                                                double &estold,double &temp)
  {
//--- copy
   isave[0]=i;
   isave[1]=iter;
   isave[2]=j;
   isave[3]=jlast;
   isave[4]=jump;
//--- copy
   rsave[0]=absxi;
   rsave[1]=altsgn;
   rsave[2]=estold;
   rsave[3]=temp;
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CRCond::InternalComplexRCondLoadAll(int &isave[],double &rsave[],
                                                int &i,int &iter,int &j,
                                                int &jlast,int &jump,
                                                double &absxi,double &altsgn,
                                                double &estold,double &temp)
  {
//--- get
   i=isave[0];
   iter=isave[1];
   j=isave[2];
   jlast=isave[3];
   jump=isave[4];
//--- get
   absxi=rsave[0];
   altsgn=rsave[1];
   estold=rsave[2];
   temp=rsave[3];
  }
//+------------------------------------------------------------------+
//| Matrix inverse report:                                           |
//| * R1    reciprocal of condition number in 1-norm                 |
//| * RInf  reciprocal of condition number in inf-norm               |
//+------------------------------------------------------------------+
class CMatInvReport
  {
public:
   double            m_r1;
   double            m_rinf;
   //--- constructor, destructor
                     CMatInvReport(void);
                    ~CMatInvReport(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMatInvReport::CMatInvReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMatInvReport::~CMatInvReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Matrix inverse report:                                           |
//| * R1    reciprocal of condition number in 1-norm                 |
//| * RInf  reciprocal of condition number in inf-norm               |
//+------------------------------------------------------------------+
class CMatInvReportShell
  {
private:
   CMatInvReport     m_innerobj;
public:
   //--- constructors, destructor
                     CMatInvReportShell(void);
                     CMatInvReportShell(CMatInvReport &obj);
                    ~CMatInvReportShell(void);
   //--- methods
   double            GetR1(void);
   void              SetR1(double r);
   double            GetRInf(void);
   void              SetRInf(double r);
   CMatInvReport    *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMatInvReportShell::CMatInvReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
CMatInvReportShell::CMatInvReportShell(CMatInvReport &obj)
  {
   m_innerobj.m_r1=obj.m_r1;
   m_innerobj.m_rinf=obj.m_rinf;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMatInvReportShell::~CMatInvReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable r1                             |
//+------------------------------------------------------------------+
double CMatInvReportShell::GetR1(void)
  {
   return(m_innerobj.m_r1);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable r1                            |
//+------------------------------------------------------------------+
void CMatInvReportShell::SetR1(double r)
  {
   m_innerobj.m_r1=r;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable rinf                           |
//+------------------------------------------------------------------+
double CMatInvReportShell::GetRInf(void)
  {
   return(m_innerobj.m_rinf);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable rint                          |
//+------------------------------------------------------------------+
void CMatInvReportShell::SetRInf(double r)
  {
   m_innerobj.m_rinf=r;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMatInvReport *CMatInvReportShell::GetInnerObj(void)
  {
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Inverse matrix                                                   |
//+------------------------------------------------------------------+
class CMatInv
  {
private:
   //--- private methods
   static void       RMatrixTrInverseRec(CMatrixDouble &a,const int offs,const int n,const bool isupper,const bool isunit,double &tmp[],int &info,CMatInvReport &rep);
   static void       CMatrixTrInverseRec(CMatrixComplex &a,const int offs,const int n,const bool isupper,const bool isunit,al_complex &tmp[],int &info,CMatInvReport &rep);
   static void       RMatrixLUInverseRec(CMatrixDouble &a,const int offs,const int n,double &work[],int &info,CMatInvReport &rep);
   static void       CMatrixLUInverseRec(CMatrixComplex &a,const int offs,const int n,al_complex &work[],int &info,CMatInvReport &rep);
   static void       SPDMatrixCholeskyInverseRec(CMatrixDouble &a,const int offs,const int n,const bool isupper,double &tmp[]);
   static void       HPDMatrixCholeskyInverseRec(CMatrixComplex &a,const int offs,const int n,const bool isupper,al_complex &tmp[]);
public:
                     CMatInv(void);
                    ~CMatInv(void);
   //--- public methods
   static void       RMatrixLUInverse(CMatrixDouble &a,int &pivots[],const int n,int &info,CMatInvReport &rep);
   static void       RMatrixInverse(CMatrixDouble &a,const int n,int &info,CMatInvReport &rep);
   static void       SPDMatrixCholeskyInverse(CMatrixDouble &a,const int n,const bool isupper,int &info,CMatInvReport &rep);
   static void       SPDMatrixInverse(CMatrixDouble &a,const int n,const bool isupper,int &info,CMatInvReport &rep);
   static void       RMatrixTrInverse(CMatrixDouble &a,const int n,const bool isupper,const bool isunit,int &info,CMatInvReport &rep);
   static void       CMatrixLUInverse(CMatrixComplex &a,int &pivots[],const int n,int &info,CMatInvReport &rep);
   static void       CMatrixInverse(CMatrixComplex &a,const int n,int &info,CMatInvReport &rep);
   static void       HPDMatrixCholeskyInverse(CMatrixComplex &a,const int n,const bool isupper,int &info,CMatInvReport &rep);
   static void       HPDMatrixInverse(CMatrixComplex &a,const int n,const bool isupper,int &info,CMatInvReport &rep);
   static void       CMatrixTrInverse(CMatrixComplex &a,const int n,const bool isupper,const bool isunit,int &info,CMatInvReport &rep);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMatInv::CMatInv(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMatInv::~CMatInv(void)
  {

  }
//+------------------------------------------------------------------+
//| Inversion of a matrix given by its LU decomposition.             |
//| INPUT PARAMETERS:                                                |
//|     A       -   LU decomposition of the matrix                   |
//|                 (output of RMatrixLU subroutine).                |
//|     Pivots  -   table of permutations                            |
//|                 (the output of RMatrixLU subroutine).            |
//|     N       -   size of matrix A (optional) :                    |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, size is automatically determined |
//|                   from matrix size (A must be square matrix)     |
//| OUTPUT PARAMETERS:                                               |
//|     Info    -   return code:                                     |
//|                 * -3    A is singular, or VERY close to singular.|
//|                         it is filled by zeros in such cases.     |
//|                 *  1    task is solved (but matrix A may be      |
//|                         ill-conditioned, check R1/RInf parameters|
//|                         for condition numbers).                  |
//|     Rep     -   solver report, see below for more info           |
//|     A       -   inverse of matrix A.                             |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..N-1].                                |
//| SOLVER REPORT                                                    |
//| Subroutine sets following fields of the Rep structure:           |
//| * R1        reciprocal of condition number: 1/cond(A), 1-norm.   |
//| * RInf      reciprocal of condition number: 1/cond(A), inf-norm. |
//+------------------------------------------------------------------+
static void CMatInv::RMatrixLUInverse(CMatrixDouble &a,int &pivots[],
                                      const int n,int &info,CMatInvReport &rep)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   double v=0;
//--- create array
   double work[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(pivots)>=n,__FUNCTION__+": len(Pivots)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(a,n,n),__FUNCTION__+": A contains infinite or NaN values!"))
      return;
//--- initialization
   info=1;
   for(i=0;i<n;i++)
     {
      //--- check
      if(pivots[i]>n-1 || pivots[i]<i)
         info=-1;
     }
//--- check
   if(!CAp::Assert(info>0,__FUNCTION__+": incorrect Pivots array!"))
      return;
//--- calculate condition numbers
   rep.m_r1=CRCond::RMatrixLURCond1(a,n);
   rep.m_rinf=CRCond::RMatrixLURCondInf(a,n);
//--- check
   if(rep.m_r1<CRCond::RCondThreshold() || rep.m_rinf<CRCond::RCondThreshold())
     {
      for(i=0;i<n;i++)
        {
         for(j=0;j<n;j++)
            a[i].Set(j,0);
        }
      //--- change values
      rep.m_r1=0;
      rep.m_rinf=0;
      info=-3;
      //--- exit the function
      return;
     }
//--- Call cache-oblivious code
   ArrayResizeAL(work,n);
   RMatrixLUInverseRec(a,0,n,work,info,rep);
//--- apply permutations
   for(i=0;i<n;i++)
     {
      for(j=n-2;j>=0;j--)
        {
         k=pivots[j];
         v=a[i][j];
         a[i].Set(j,a[i][k]);
         a[i].Set(k,v);
        }
     }
  }
//+------------------------------------------------------------------+
//| Inversion of a general matrix.                                   |
//| Input parameters:                                                |
//|     A       -   matrix.                                          |
//|     N       -   size of matrix A (optional) :                    |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, size is automatically determined |
//|                   from matrix size (A must be square matrix)     |
//| Output parameters:                                               |
//|     Info    -   return code, same as in RMatrixLUInverse         |
//|     Rep     -   solver report, same as in RMatrixLUInverse       |
//|     A       -   inverse of matrix A, same as in RMatrixLUInverse |
//| Result:                                                          |
//|     True, if the matrix is not singular.                         |
//|     False, if the matrix is singular.                            |
//+------------------------------------------------------------------+
static void CMatInv::RMatrixInverse(CMatrixDouble &a,const int n,int &info,
                                    CMatInvReport &rep)
  {
//--- create array
   int pivots[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(a,n,n),__FUNCTION__+": A contains infinite or NaN values!"))
      return;
//--- function call
   CTrFac::RMatrixLU(a,n,n,pivots);
//--- function call
   RMatrixLUInverse(a,pivots,n,info,rep);
  }
//+------------------------------------------------------------------+
//| Inversion of a matrix given by its LU decomposition.             |
//| INPUT PARAMETERS:                                                |
//|     A       -   LU decomposition of the matrix                   |
//|                 (output of CMatrixLU subroutine).                |
//|     Pivots  -   table of permutations                            |
//|                 (the output of CMatrixLU subroutine).            |
//|     N       -   size of matrix A (optional) :                    |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, size is automatically determined |
//|                   from matrix size (A must be square matrix)     |
//| OUTPUT PARAMETERS:                                               |
//|     Info    -   return code, same as in RMatrixLUInverse         |
//|     Rep     -   solver report, same as in RMatrixLUInverse       |
//|     A       -   inverse of matrix A, same as in RMatrixLUInverse |
//+------------------------------------------------------------------+
static void CMatInv::CMatrixLUInverse(CMatrixComplex &a,int &pivots[],
                                      const int n,int &info,CMatInvReport &rep)
  {
//--- create variables
   int        i=0;
   int        j=0;
   int        k=0;
   al_complex v=0;
//--- create array
   al_complex work[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(pivots)>=n,__FUNCTION__+": len(Pivots)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteComplexMatrix(a,n,n),__FUNCTION__+": A contains infinite or NaN values!"))
      return;
//--- initialization
   info=1;
   for(i=0;i<n;i++)
     {
      //--- check
      if(pivots[i]>n-1 || pivots[i]<i)
         info=-1;
     }
//--- check
   if(!CAp::Assert(info>0,__FUNCTION__+": incorrect Pivots array!"))
      return;
//--- calculate condition numbers
   rep.m_r1=CRCond::CMatrixLURCond1(a,n);
   rep.m_rinf=CRCond::CMatrixLURCondInf(a,n);
//--- check
   if(rep.m_r1<CRCond::RCondThreshold() || rep.m_rinf<CRCond::RCondThreshold())
     {
      for(i=0;i<n;i++)
        {
         for(j=0;j<n;j++)
            a[i].Set(j,0);
        }
      //--- change values
      rep.m_r1=0;
      rep.m_rinf=0;
      info=-3;
      //--- exit the function
      return;
     }
//--- Call cache-oblivious code
   ArrayResizeAL(work,n);
   CMatrixLUInverseRec(a,0,n,work,info,rep);
//--- apply permutations
   for(i=0;i<n;i++)
     {
      for(j=n-2;j>=0;j--)
        {
         k=pivots[j];
         v=a[i][j];
         a[i].Set(j,a[i][k]);
         a[i].Set(k,v);
        }
     }
  }
//+------------------------------------------------------------------+
//| Inversion of a general matrix.                                   |
//| Input parameters:                                                |
//|     A       -   matrix                                           |
//|     N       -   size of matrix A (optional) :                    |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, size is automatically determined |
//|                   from matrix size (A must be square matrix)     |
//| Output parameters:                                               |
//|     Info    -   return code, same as in RMatrixLUInverse         |
//|     Rep     -   solver report, same as in RMatrixLUInverse       |
//|     A       -   inverse of matrix A, same as in RMatrixLUInverse |
//+------------------------------------------------------------------+
static void CMatInv::CMatrixInverse(CMatrixComplex &a,const int n,int &info,
                                    CMatInvReport &rep)
  {
//--- create array
   int pivots[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteComplexMatrix(a,n,n),__FUNCTION__+": A contains infinite or NaN values!"))
      return;
//--- function call
   CTrFac::CMatrixLU(a,n,n,pivots);
//--- function call
   CMatrixLUInverse(a,pivots,n,info,rep);
  }
//+------------------------------------------------------------------+
//| Inversion of a symmetric positive definite matrix which is given |
//| by Cholesky decomposition.                                       |
//| Input parameters:                                                |
//|     A       -   Cholesky decomposition of the matrix to be       |
//|                 inverted: A=U?*U or A = L*L'.                    |
//|                 Output of  SPDMatrixCholesky subroutine.         |
//|     N       -   size of matrix A (optional) :                    |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, size is automatically determined |
//|                   from matrix size (A must be square matrix)     |
//|     IsUpper -   storage type (optional):                         |
//|                 * if True, symmetric matrix A is given by its    |
//|                   upper triangle, and the lower triangle isn?t   |
//|                   used/changed by function                       |
//|                 * if False, symmetric matrix A is given by its   |
//|                   lower triangle, and the upper triangle isn?t   |
//|                   used/changed  by function                      |
//|                 * if not given, lower half is used.              |
//| Output parameters:                                               |
//|     Info    -   return code, same as in RMatrixLUInverse         |
//|     Rep     -   solver report, same as in RMatrixLUInverse       |
//|     A       -   inverse of matrix A, same as in RMatrixLUInverse |
//+------------------------------------------------------------------+
static void CMatInv::SPDMatrixCholeskyInverse(CMatrixDouble &a,const int n,
                                              const bool isupper,int &info,
                                              CMatInvReport &rep)
  {
//--- create variables
   int  i=0;
   int  j=0;
   bool f;
//--- create array
   double tmp[];
//--- object of class
   CMatInvReport rep2;
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return;
//--- initialization
   info=1;
   f=true;
   for(i=0;i<n;i++)
      f=f && CMath::IsFinite(a[i][i]);
//--- check
   if(!CAp::Assert(f,__FUNCTION__+": A contains infinite or NaN values!"))
      return;
//--- calculate condition numbers
   rep.m_r1=CRCond::SPDMatrixCholeskyRCond(a,n,isupper);
   rep.m_rinf=rep.m_r1;
//--- check
   if(rep.m_r1<CRCond::RCondThreshold() || rep.m_rinf<CRCond::RCondThreshold())
     {//--- check
      if(isupper)
        {
         for(i=0;i<n;i++)
            for(j=i;j<n;j++)
               a[i].Set(j,0);
        }
      else
        {
         for(i=0;i<n;i++)
            for(j=0;j<=i;j++)
               a[i].Set(j,0);
        }
      //--- change values
      rep.m_r1=0;
      rep.m_rinf=0;
      info=-3;
      //--- exit the function
      return;
     }
//--- Inverse
   ArrayResizeAL(tmp,n);
   SPDMatrixCholeskyInverseRec(a,0,n,isupper,tmp);
  }
//+------------------------------------------------------------------+
//| Inversion of a symmetric positive definite matrix.               |
//| Given an upper or lower triangle of a symmetric positive definite|
//| matrix, the algorithm generates matrix A^-1 and saves the upper  |
//| or lower triangle depending on the input.                        |
//| Input parameters:                                                |
//|     A       -   matrix to be inverted (upper or lower triangle). |
//|                 Array with elements [0..N-1,0..N-1].             |
//|     N       -   size of matrix A (optional) :                    |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, size is automatically determined |
//|                   from matrix size (A must be square matrix)     |
//|     IsUpper -   storage type (optional):                         |
//|                 * if True, symmetric matrix A is given by its    |
//|                   upper triangle, and the lower triangle isn?t   |
//|                   used/changed by function                       |
//|                 * if False, symmetric matrix A is given by its   |
//|                   lower triangle, and the upper triangle isn?t   |
//|                   used/changed by function                       |
//|                 * if not given, both lower and upper triangles   |
//|                   must be filled.                                |
//| Output parameters:                                               |
//|     Info    -   return code, same as in RMatrixLUInverse         |
//|     Rep     -   solver report, same as in RMatrixLUInverse       |
//|     A       -   inverse of matrix A, same as in RMatrixLUInverse |
//+------------------------------------------------------------------+
static void CMatInv::SPDMatrixInverse(CMatrixDouble &a,const int n,
                                      const bool isupper,int &info,
                                      CMatInvReport &rep)
  {
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteRTrMatrix(a,n,isupper),__FUNCTION__+": A contains infinite or NaN values!"))
      return;
//--- initialization
   info=1;
//--- check
   if(CTrFac::SPDMatrixCholesky(a,n,isupper))
      SPDMatrixCholeskyInverse(a,n,isupper,info,rep);
   else
      info=-3;
  }
//+------------------------------------------------------------------+
//| Inversion of a Hermitian positive definite matrix which is given |
//| by Cholesky decomposition.                                       |
//| Input parameters:                                                |
//|     A       -   Cholesky decomposition of the matrix to be       |
//|                 inverted: A=U?*U or A = L*L'.                    |
//|                 Output of  HPDMatrixCholesky subroutine.         |
//|     N       -   size of matrix A (optional) :                    |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, size is automatically determined |
//|                   from matrix size (A must be square matrix)     |
//|     IsUpper -   storage type (optional):                         |
//|                 * if True, symmetric matrix A is given by its    |
//|                   upper triangle, and the lower triangle isn?t   |
//|                   used/changed by function                       |
//|                 * if False, symmetric matrix A is given by its   |
//|                   lower triangle, and the upper triangle isn?t   |
//|                   used/changed by function                       |
//|                 * if not given, lower half is used.              |
//| Output parameters:                                               |
//|     Info    -   return code, same as in RMatrixLUInverse         |
//|     Rep     -   solver report, same as in RMatrixLUInverse       |
//|     A       -   inverse of matrix A, same as in RMatrixLUInverse |
//+------------------------------------------------------------------+
static void CMatInv::HPDMatrixCholeskyInverse(CMatrixComplex &a,const int n,
                                              const bool isupper,int &info,
                                              CMatInvReport &rep)
  {
//--- create variables
   int  i=0;
   int  j=0;
   bool f;
//--- create array
   al_complex tmp[];
//--- object of class
   CMatInvReport rep2;
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return;
//--- initialization
   f=true;
   for(i=0;i<n;i++)
      f=(f && CMath::IsFinite(a[i][i].re)) && CMath::IsFinite(a[i][i].im);
//--- check
   if(!CAp::Assert(f,__FUNCTION__+": A contains infinite or NaN values!"))
      return;
   info=1;
//--- calculate condition numbers
   rep.m_r1=CRCond::HPDMatrixCholeskyRCond(a,n,isupper);
   rep.m_rinf=rep.m_r1;
//--- check
   if(rep.m_r1<CRCond::RCondThreshold() || rep.m_rinf<CRCond::RCondThreshold())
     {
      //--- check
      if(isupper)
        {
         for(i=0;i<n;i++)
            for(j=i;j<n;j++)
               a[i].Set(j,0);
        }
      else
        {
         for(i=0;i<n;i++)
            for(j=0;j<=i;j++)
               a[i].Set(j,0);
        }
      //--- change values
      rep.m_r1=0;
      rep.m_rinf=0;
      info=-3;
      //--- exit the function
      return;
     }
//--- Inverse
   ArrayResizeAL(tmp,n);
   HPDMatrixCholeskyInverseRec(a,0,n,isupper,tmp);
  }
//+------------------------------------------------------------------+
//| Inversion of a Hermitian positive definite matrix.               |
//| Given an upper or lower triangle of a Hermitian positive definite|
//| matrix, the algorithm generates matrix A^-1 and saves the upper  |
//| or lower triangle depending on the input.                        |
//| Input parameters:                                                |
//|     A       -   matrix to be inverted (upper or lower triangle). |
//|                 Array with elements [0..N-1,0..N-1].             |
//|     N       -   size of matrix A (optional) :                    |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, size is automatically determined |
//|                   from matrix size (A must be square matrix)     |
//|     IsUpper -   storage type (optional):                         |
//|                 * if True, symmetric matrix A is given by its    |
//|                   upper triangle, and the lower triangle isn?t   |
//|                   used/changed by function                       |
//|                 * if False, symmetric matrix A is given by its   |
//|                   lower triangle, and the upper triangle isn?t   |
//|                   used/changed by function                       |
//|                 * if not given, both lower and upper triangles   |
//|                   must be filled.                                |
//| Output parameters:                                               |
//|     Info    -   return code, same as in RMatrixLUInverse         |
//|     Rep     -   solver report, same as in RMatrixLUInverse       |
//|     A       -   inverse of matrix A, same as in RMatrixLUInverse |
//+------------------------------------------------------------------+
static void CMatInv::HPDMatrixInverse(CMatrixComplex &a,const int n,
                                      const bool isupper,int &info,
                                      CMatInvReport &rep)
  {
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteCTrMatrix(a,n,isupper),__FUNCTION__+": A contains infinite or NaN values!"))
      return;
//--- initialization
   info=1;
//--- check
   if(CTrFac::HPDMatrixCholesky(a,n,isupper))
      HPDMatrixCholeskyInverse(a,n,isupper,info,rep);
   else
      info=-3;
  }
//+------------------------------------------------------------------+
//| Triangular matrix inverse (real)                                 |
//| The subroutine inverts the following types of matrices:          |
//|     * upper triangular                                           |
//|     * upper triangular with unit diagonal                        |
//|     * lower triangular                                           |
//|     * lower triangular with unit diagonal                        |
//| In case of an upper (lower) triangular matrix, the inverse matrix|
//| will also be upper (lower) triangular, and after the end of the  |
//| algorithm, the inverse matrix replaces the source matrix. The    |
//| elements below (above) the main diagonal are not changed by the  |
//| algorithm.                                                       |
//| If the matrix has a unit diagonal, the inverse matrix also has a |
//| unit diagonal, and the diagonal elements are not passed to the   |
//| algorithm.                                                       |
//| Input parameters:                                                |
//|     A       -   matrix, array[0..N-1, 0..N-1].                   |
//|     N       -   size of matrix A (optional) :                    |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, size is automatically determined |
//|                   from matrix size (A must be square matrix)     |
//|     IsUpper -   True, if the matrix is upper triangular.         |
//|     IsUnit  -   diagonal type (optional):                        |
//|                 * if True, matrix has unit diagonal (a[i,i] are  |
//|                   NOT used)                                      |
//|                 * if False, matrix diagonal is arbitrary         |
//|                 * if not given, False is assumed                 |
//| Output parameters:                                               |
//|     Info    -   same as for RMatrixLUInverse                     |
//|     Rep     -   same as for RMatrixLUInverse                     |
//|     A       -   same as for RMatrixLUInverse.                    |
//+------------------------------------------------------------------+
static void CMatInv::RMatrixTrInverse(CMatrixDouble &a,const int n,
                                      const bool isupper,const bool isunit,
                                      int &info,CMatInvReport &rep)
  {
//--- create variables
   int i=0;
   int j=0;
//--- create array
   double tmp[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteRTrMatrix(a,n,isupper),__FUNCTION__+": A contains infinite or NaN values!"))
      return;
//--- initialization
   info=1;
//--- calculate condition numbers
   rep.m_r1=CRCond::RMatrixTrRCond1(a,n,isupper,isunit);
   rep.m_rinf=CRCond::RMatrixTrRCondInf(a,n,isupper,isunit);
//--- check
   if(rep.m_r1<CRCond::RCondThreshold() || rep.m_rinf<CRCond::RCondThreshold())
     {
      for(i=0;i<n;i++)
         for(j=0;j<n;j++)
            a[i].Set(j,0);
      //--- change values
      rep.m_r1=0;
      rep.m_rinf=0;
      info=-3;
      //--- exit the function
      return;
     }
//--- Inverse
   ArrayResizeAL(tmp,n);
   RMatrixTrInverseRec(a,0,n,isupper,isunit,tmp,info,rep);
  }
//+------------------------------------------------------------------+
//| Triangular matrix inverse (complex)                              |
//| The subroutine inverts the following types of matrices:          |
//|     * upper triangular                                           |
//|     * upper triangular with unit diagonal                        |
//|     * lower triangular                                           |
//|     * lower triangular with unit diagonal                        |
//| In case of an upper (lower) triangular matrix, the inverse matrix|
//| will also be upper (lower) triangular, and after the end of the  |
//| algorithm, the inverse matrix replaces the source matrix. The    |
//| elements below (above) the main diagonal are not changed by the  |
//| algorithm.                                                       |
//| If the matrix has a unit diagonal, the inverse matrix also has a |
//| unit diagonal, and the diagonal elements are not passed to the   |
//| algorithm.                                                       |
//| Input parameters:                                                |
//|     A       -   matrix, array[0..N-1, 0..N-1].                   |
//|     N       -   size of matrix A (optional) :                    |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, size is automatically determined |
//|                   from matrix size (A must be square matrix)     |
//|     IsUpper -   True, if the matrix is upper triangular.         |
//|     IsUnit  -   diagonal type (optional):                        |
//|                 * if True, matrix has unit diagonal (a[i,i] are  |
//|                   NOT used)                                      |
//|                 * if False, matrix diagonal is arbitrary         |
//|                 * if not given, False is assumed                 |
//| Output parameters:                                               |
//|     Info    -   same as for RMatrixLUInverse                     |
//|     Rep     -   same as for RMatrixLUInverse                     |
//|     A       -   same as for RMatrixLUInverse.                    |
//+------------------------------------------------------------------+
static void CMatInv::CMatrixTrInverse(CMatrixComplex &a,const int n,
                                      const bool isupper,const bool isunit,
                                      int &info,CMatInvReport &rep)
  {
//--- create variables
   int i=0;
   int j=0;
//--- create array
   al_complex tmp[];
//--- initialization
   info=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteCTrMatrix(a,n,isupper),__FUNCTION__+": A contains infinite or NaN values!"))
      return;
//--- initialization
   info=1;
//--- calculate condition numbers
   rep.m_r1=CRCond::CMatrixTrRCond1(a,n,isupper,isunit);
   rep.m_rinf=CRCond::CMatrixTrRCondInf(a,n,isupper,isunit);
//--- check
   if(rep.m_r1<CRCond::RCondThreshold() || rep.m_rinf<CRCond::RCondThreshold())
     {
      for(i=0;i<n;i++)
        {
         for(j=0;j<n;j++)
            a[i].Set(j,0);
        }
      //--- change values
      rep.m_r1=0;
      rep.m_rinf=0;
      info=-3;
      //--- exit the function
      return;
     }
//--- Inverse
   ArrayResizeAL(tmp,n);
   CMatrixTrInverseRec(a,0,n,isupper,isunit,tmp,info,rep);
  }
//+------------------------------------------------------------------+
//| Triangular matrix inversion, recursive subroutine                |
//+------------------------------------------------------------------+
static void CMatInv::RMatrixTrInverseRec(CMatrixDouble &a,const int offs,
                                         const int n,const bool isupper,
                                         const bool isunit,double &tmp[],
                                         int &info,CMatInvReport &rep)
  {
//--- create variables
   int    n1=0;
   int    n2=0;
   int    i=0;
   int    j=0;
   double v=0;
   double ajj=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(n<1)
     {
      info=-1;
      //--- exit the function
      return;
     }
//--- base case
   if(n<=CAblas::AblasBlockSize())
     {
      //--- check
      if(isupper)
        {
         //--- Compute inverse of upper triangular matrix.
         for(j=0;j<n;j++)
           {
            //--- check
            if(!isunit)
              {
               //--- check
               if(a[offs+j][offs+j]==0.0)
                 {
                  info=-3;
                  //--- exit the function
                  return;
                 }
               a[offs+j].Set(offs+j,1/a[offs+j][offs+j]);
               ajj=-a[offs+j][offs+j];
              }
            else
               ajj=-1;
            //--- Compute elements 1:j-1 of j-th column.
            if(j>0)
              {
               i1_=offs;
               for(i_=0;i_<j;i_++)
                  tmp[i_]=a[i_+i1_][offs+j];
               for(i=0;i<j;i++)
                 {
                  //--- check
                  if(i<j-1)
                    {
                     i1_=(i+1)-(offs+i+1);
                     v=0.0;
                     for(i_=offs+i+1;i_<=offs+j-1;i_++)
                        v+=a[offs+i][i_]*tmp[i_+i1_];
                    }
                  else
                     v=0;
                  //--- check
                  if(!isunit)
                     a[offs+i].Set(offs+j,v+a[offs+i][offs+i]*tmp[i]);
                  else
                     a[offs+i].Set(offs+j,v+tmp[i]);
                 }
               for(i_=offs+0;i_<=offs+j-1;i_++)
                  a[i_].Set(offs+j,ajj*a[i_][offs+j]);
              }
           }
        }
      else
        {
         //--- Compute inverse of lower triangular matrix.
         for(j=n-1;j>=0;j--)
           {
            //--- check
            if(!isunit)
              {
               //--- check
               if(a[offs+j][offs+j]==0.0)
                 {
                  info=-3;
                  //--- exit the function
                  return;
                 }
               a[offs+j].Set(offs+j,1/a[offs+j][offs+j]);
               ajj=-a[offs+j][offs+j];
              }
            else
               ajj=-1;
            //--- check
            if(j<n-1)
              {
               //--- Compute elements j+1:n of j-th column.
               i1_=offs;
               for(i_=j+1;i_<n;i_++)
                  tmp[i_]=a[i_+i1_][offs+j];
               for(i=j+1;i<n;i++)
                 {
                  //--- check
                  if(i>j+1)
                    {
                     i1_=-offs;
                     v=0.0;
                     for(i_=offs+j+1;i_<offs+i;i_++)
                        v+=a[offs+i][i_]*tmp[i_+i1_];
                    }
                  else
                     v=0;
                  //--- check
                  if(!isunit)
                     a[offs+i].Set(offs+j,v+a[offs+i][offs+i]*tmp[i]);
                  else
                     a[offs+i].Set(offs+j,v+tmp[i]);
                 }
               for(i_=offs+j+1;i_<=offs+n-1;i_++)
                  a[i_].Set(offs+j,ajj*a[i_][offs+j]);
              }
           }
        }
      //--- exit the function
      return;
     }
//--- Recursive case
   CAblas::AblasSplitLength(a,n,n1,n2);
//--- check
   if(n2>0)
     {
      //--- check
      if(isupper)
        {
         for(i=0;i<n1;i++)
            for(i_=offs+n1;i_<=offs+n-1;i_++)
               a[offs+i].Set(i_,-1*a[offs+i][i_]);
         //--- function call
         CAblas::RMatrixLeftTrsM(n1,n2,a,offs,offs,isupper,isunit,0,a,offs,offs+n1);
         //--- function call
         CAblas::RMatrixRightTrsM(n1,n2,a,offs+n1,offs+n1,isupper,isunit,0,a,offs,offs+n1);
        }
      else
        {
         for(i=0;i<n2;i++)
            for(i_=offs;i_<=offs+n1-1;i_++)
               a[offs+n1+i].Set(i_,-1*a[offs+n1+i][i_]);
         //--- function call
         CAblas::RMatrixRightTrsM(n2,n1,a,offs,offs,isupper,isunit,0,a,offs+n1,offs);
         //--- function call
         CAblas::RMatrixLeftTrsM(n2,n1,a,offs+n1,offs+n1,isupper,isunit,0,a,offs+n1,offs);
        }
      //--- function call
      RMatrixTrInverseRec(a,offs+n1,n2,isupper,isunit,tmp,info,rep);
     }
//--- function call
   RMatrixTrInverseRec(a,offs,n1,isupper,isunit,tmp,info,rep);
  }
//+------------------------------------------------------------------+
//| Triangular matrix inversion, recursive subroutine                |
//+------------------------------------------------------------------+
static void CMatInv::CMatrixTrInverseRec(CMatrixComplex &a,const int offs,
                                         const int n,const bool isupper,
                                         const bool isunit,al_complex &tmp[],
                                         int &info,CMatInvReport &rep)
  {
//--- create variables
   al_complex One=1.0;
   al_complex _One=-1.0;
   int        n1=0;
   int        n2=0;
   int        i=0;
   int        j=0;
   al_complex v=0;
   al_complex ajj=0;
   int        i_=0;
   int        i1_=0;
//--- check
   if(n<1)
     {
      info=-1;
      //--- exit the function
      return;
     }
//--- base case
   if(n<=CAblas::AblasComplexBlockSize())
     {
      //--- check
      if(isupper)
        {
         //--- Compute inverse of upper triangular matrix.
         for(j=0;j<n;j++)
           {
            //--- check
            if(!isunit)
              {
               //--- check
               if(a[offs+j][offs+j]==0)
                 {
                  info=-3;
                  //--- exit the function
                  return;
                 }
               a[offs+j].Set(offs+j,One/a[offs+j][offs+j]);
               ajj=-a[offs+j][offs+j];
              }
            else
               ajj=-1;
            //--- Compute elements 1:j-1 of j-th column.
            if(j>0)
              {
               i1_=offs;
               for(i_=0;i_<j;i_++)
                  tmp[i_]=a[i_+i1_][offs+j];
               for(i=0;i<j;i++)
                 {
                  //--- check
                  if(i<j-1)
                    {
                     i1_=(i+1)-(offs+i+1);
                     v=0.0;
                     for(i_=offs+i+1;i_<=offs+j-1;i_++)
                        v+=a[offs+i][i_]*tmp[i_+i1_];
                    }
                  else
                     v=0;
                  //--- check
                  if(!isunit)
                     a[offs+i].Set(offs+j,v+a[offs+i][offs+i]*tmp[i]);
                  else
                     a[offs+i].Set(offs+j,v+tmp[i]);
                 }
               for(i_=offs+0;i_<=offs+j-1;i_++)
                  a[i_].Set(offs+j,ajj*a[i_][offs+j]);
              }
           }
        }
      else
        {
         //--- Compute inverse of lower triangular matrix.
         for(j=n-1;j>=0;j--)
           {
            //--- check
            if(!isunit)
              {
               //--- check
               if(a[offs+j][offs+j]==0)
                 {
                  info=-3;
                  //--- exit the function
                  return;
                 }
               a[offs+j].Set(offs+j,One/a[offs+j][offs+j]);
               ajj=-a[offs+j][offs+j];
              }
            else
               ajj=-1;
            //--- check
            if(j<n-1)
              {
               //--- Compute elements j+1:n of j-th column.
               i1_=offs;
               for(i_=j+1;i_<n;i_++)
                  tmp[i_]=a[i_+i1_][offs+j];
               for(i=j+1;i<n;i++)
                 {
                  //--- check
                  if(i>j+1)
                    {
                     i1_=-offs;
                     v=0.0;
                     for(i_=offs+j+1;i_<offs+i;i_++)
                        v+=a[offs+i][i_]*tmp[i_+i1_];
                    }
                  else
                     v=0;
                  //--- check
                  if(!isunit)
                     a[offs+i].Set(offs+j,v+a[offs+i][offs+i]*tmp[i]);
                  else
                     a[offs+i].Set(offs+j,v+tmp[i]);
                 }
               for(i_=offs+j+1;i_<=offs+n-1;i_++)
                  a[i_].Set(offs+j,ajj*a[i_][offs+j]);
              }
           }
        }
      //--- exit the function
      return;
     }
//--- Recursive case
   CAblas::AblasComplexSplitLength(a,n,n1,n2);
//--- check
   if(n2>0)
     {
      //--- check
      if(isupper)
        {
         for(i=0;i<n1;i++)
           {
            for(i_=offs+n1;i_<=offs+n-1;i_++)
               a[offs+i].Set(i_,_One*a[offs+i][i_]);
           }
         //--- function call
         CAblas::CMatrixLeftTrsM(n1,n2,a,offs,offs,isupper,isunit,0,a,offs,offs+n1);
         //--- function call
         CAblas::CMatrixRightTrsM(n1,n2,a,offs+n1,offs+n1,isupper,isunit,0,a,offs,offs+n1);
        }
      else
        {
         for(i=0;i<n2;i++)
           {
            for(i_=offs;i_<=offs+n1-1;i_++)
               a[offs+n1+i].Set(i_,_One*a[offs+n1+i][i_]);
           }
         //--- function call
         CAblas::CMatrixRightTrsM(n2,n1,a,offs,offs,isupper,isunit,0,a,offs+n1,offs);
         //--- function call
         CAblas::CMatrixLeftTrsM(n2,n1,a,offs+n1,offs+n1,isupper,isunit,0,a,offs+n1,offs);
        }
      //--- function call
      CMatrixTrInverseRec(a,offs+n1,n2,isupper,isunit,tmp,info,rep);
     }
//--- function call
   CMatrixTrInverseRec(a,offs,n1,isupper,isunit,tmp,info,rep);
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CMatInv::RMatrixLUInverseRec(CMatrixDouble &a,const int offs,
                                         const int n,double &work[],
                                         int &info,CMatInvReport &rep)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   int    n1=0;
   int    n2=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(n<1)
     {
      info=-1;
      //--- exit the function
      return;
     }
//--- base case
   if(n<=CAblas::AblasBlockSize())
     {
      //--- Form inv(U)
      RMatrixTrInverseRec(a,offs,n,true,false,work,info,rep);
      //--- check
      if(info<=0)
         return;
      //--- Solve the equation inv(A)*L = inv(U) for inv(A).
      for(j=n-1;j>=0;j--)
        {
         //--- Copy current column of L to WORK and replace with zeros.
         for(i=j+1;i<n;i++)
           {
            work[i]=a[offs+i][offs+j];
            a[offs+i].Set(offs+j,0);
           }
         //--- Compute current column of inv(A).
         if(j<n-1)
           {
            for(i=0;i<n;i++)
              {
               i1_=-offs;
               v=0.0;
               for(i_=offs+j+1;i_<=offs+n-1;i_++)
                  v+=a[offs+i][i_]*work[i_+i1_];
               a[offs+i].Set(offs+j,a[offs+i][offs+j]-v);
              }
           }
        }
      //--- exit the function
      return;
     }
//--- Recursive code:
//---         ( L1      )   ( U1  U12 )
//--- A    =  (         ) * (         )
//---         ( L12  L2 )   (     U2  )
//---         ( W   X )
//--- A^-1 =  (       )
//---         ( Y   Z )
   CAblas::AblasSplitLength(a,n,n1,n2);
//--- check
   if(!CAp::Assert(n2>0,__FUNCTION__+": internal error!"))
      return;
//--- X :=inv(U1)*U12*inv(U2)
   CAblas::RMatrixLeftTrsM(n1,n2,a,offs,offs,true,false,0,a,offs,offs+n1);
   CAblas::RMatrixRightTrsM(n1,n2,a,offs+n1,offs+n1,true,false,0,a,offs,offs+n1);
//--- Y :=inv(L2)*L12*inv(L1)
   CAblas::RMatrixLeftTrsM(n2,n1,a,offs+n1,offs+n1,false,true,0,a,offs+n1,offs);
   CAblas::RMatrixRightTrsM(n2,n1,a,offs,offs,false,true,0,a,offs+n1,offs);
//--- W :=inv(L1*U1)+X*Y
   RMatrixLUInverseRec(a,offs,n1,work,info,rep);
//--- check
   if(info<=0)
      return;
//--- function call
   CAblas::RMatrixGemm(n1,n1,n2,1.0,a,offs,offs+n1,0,a,offs+n1,offs,0,1.0,a,offs,offs);
//--- X :=-X*inv(L2)
//--- Y :=-inv(U2)*Y
   CAblas::RMatrixRightTrsM(n1,n2,a,offs+n1,offs+n1,false,true,0,a,offs,offs+n1);
   for(i=0;i<n1;i++)
     {
      for(i_=offs+n1;i_<=offs+n-1;i_++)
         a[offs+i].Set(i_,-1*a[offs+i][i_]);
     }
   CAblas::RMatrixLeftTrsM(n2,n1,a,offs+n1,offs+n1,true,false,0,a,offs+n1,offs);
   for(i=0;i<n2;i++)
     {
      for(i_=offs;i_<=offs+n1-1;i_++)
         a[offs+n1+i].Set(i_,-1*a[offs+n1+i][i_]);
     }
//--- Z :=inv(L2*U2)
   RMatrixLUInverseRec(a,offs+n1,n2,work,info,rep);
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CMatInv::CMatrixLUInverseRec(CMatrixComplex &a,const int offs,
                                         const int n,al_complex &work[],
                                         int &info,CMatInvReport &rep)
  {
//--- create variables
   al_complex One=1.0;
   al_complex _One=-1.0;
   int        i=0;
   int        j=0;
   al_complex v=0;
   int        n1=0;
   int        n2=0;
   int        i_=0;
   int        i1_=0;
//--- check
   if(n<1)
     {
      info=-1;
      //--- exit the function
      return;
     }
//--- base case
   if(n<=CAblas::AblasComplexBlockSize())
     {
      //--- Form inv(U)
      CMatrixTrInverseRec(a,offs,n,true,false,work,info,rep);
      //--- check
      if(info<=0)
         return;
      //--- Solve the equation inv(A)*L = inv(U) for inv(A).
      for(j=n-1;j>=0;j--)
        {
         //--- Copy current column of L to WORK and replace with zeros.
         for(i=j+1;i<n;i++)
           {
            work[i]=a[offs+i][offs+j];
            a[offs+i].Set(offs+j,0);
           }
         //--- Compute current column of inv(A).
         if(j<n-1)
           {
            for(i=0;i<n;i++)
              {
               i1_=-offs;
               v=0.0;
               for(i_=offs+j+1;i_<=offs+n-1;i_++)
                  v+=a[offs+i][i_]*work[i_+i1_];
               a[offs+i].Set(offs+j,a[offs+i][offs+j]-v);
              }
           }
        }
      //--- exit the function
      return;
     }
//--- Recursive code:
//---         ( L1      )   ( U1  U12 )
//--- A    =  (         ) * (         )
//---         ( L12  L2 )   (     U2  )
//---         ( W   X )
//--- A^-1 =  (       )
//---         ( Y   Z )
   CAblas::AblasComplexSplitLength(a,n,n1,n2);
//--- check
   if(!CAp::Assert(n2>0,__FUNCTION__+": internal error!"))
      return;
//--- X :=inv(U1)*U12*inv(U2)
   CAblas::CMatrixLeftTrsM(n1,n2,a,offs,offs,true,false,0,a,offs,offs+n1);
   CAblas::CMatrixRightTrsM(n1,n2,a,offs+n1,offs+n1,true,false,0,a,offs,offs+n1);
//--- Y :=inv(L2)*L12*inv(L1)
   CAblas::CMatrixLeftTrsM(n2,n1,a,offs+n1,offs+n1,false,true,0,a,offs+n1,offs);
   CAblas::CMatrixRightTrsM(n2,n1,a,offs,offs,false,true,0,a,offs+n1,offs);
//--- W :=inv(L1*U1)+X*Y
   CMatrixLUInverseRec(a,offs,n1,work,info,rep);
//--- check
   if(info<=0)
      return;
   CAblas::CMatrixGemm(n1,n1,n2,One,a,offs,offs+n1,0,a,offs+n1,offs,0,One,a,offs,offs);
//--- X :=-X*inv(L2)
//--- Y :=-inv(U2)*Y
   CAblas::CMatrixRightTrsM(n1,n2,a,offs+n1,offs+n1,false,true,0,a,offs,offs+n1);
   for(i=0;i<n1;i++)
     {
      for(i_=offs+n1;i_<=offs+n-1;i_++)
         a[offs+i].Set(i_,_One*a[offs+i][i_]);
     }
   CAblas::CMatrixLeftTrsM(n2,n1,a,offs+n1,offs+n1,true,false,0,a,offs+n1,offs);
   for(i=0;i<n2;i++)
     {
      for(i_=offs;i_<=offs+n1-1;i_++)
         a[offs+n1+i].Set(i_,_One*a[offs+n1+i][i_]);
     }
//--- Z :=inv(L2*U2)
   CMatrixLUInverseRec(a,offs+n1,n2,work,info,rep);
  }
//+------------------------------------------------------------------+
//| Recursive subroutine for SPD inversion.                          |													      |
//+------------------------------------------------------------------+
static void CMatInv::SPDMatrixCholeskyInverseRec(CMatrixDouble &a,const int offs,
                                                 const int n,const bool isupper,
                                                 double &tmp[])
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   int    n1=0;
   int    n2=0;
   int    info2=0;
   int    i_=0;
   int    i1_=0;
//--- object of class
   CMatInvReport rep2;
//--- check
   if(n<1)
      return;
//--- base case
   if(n<=CAblas::AblasBlockSize())
     {
      RMatrixTrInverseRec(a,offs,n,isupper,false,tmp,info2,rep2);
      //--- check
      if(isupper)
        {
         //--- Compute the product U * U'.
         //--- NOTE: we never assume that diagonal of U is real
         for(i=0;i<n;i++)
           {
            //--- check
            if(i==0)
              {
               //--- 1x1 matrix
               a[offs+i].Set(offs+i,CMath::Sqr(a[offs+i][offs+i]));
              }
            else
              {
               //--- (I+1)x(I+1) matrix,
               //--- ( A11  A12 )   ( A11^H        )   ( A11*A11^H+A12*A12^H  A12*A22^H )
               //--- (          ) * (              ) = (                                )
               //--- (      A22 )   ( A12^H  A22^H )   ( A22*A12^H            A22*A22^H )
               //--- A11 is IxI, A22 is 1x1.
               i1_=offs;
               for(i_=0;i_<=i-1;i_++)
                  tmp[i_]=a[i_+i1_][offs+i];
               for(j=0;j<i;j++)
                 {
                  v=a[offs+j][offs+i];
                  i1_=-offs;
                  for(i_=offs+j;i_<offs+i;i_++)
                     a[offs+j].Set(i_,a[offs+j][i_]+v*tmp[i_+i1_]);
                 }
               //--- change values
               v=a[offs+i][offs+i];
               for(i_=offs;i_<offs+i;i_++)
                  a[i_].Set(offs+i,v*a[i_][offs+i]);
               a[offs+i].Set(offs+i,CMath::Sqr(a[offs+i][offs+i]));
              }
           }
        }
      else
        {
         //--- Compute the product L' * L
         //--- NOTE: we never assume that diagonal of L is reall
         for(i=0;i<n;i++)
           {
            //--- check
            if(i==0)
              {
               //--- 1x1 matrix
               a[offs+i].Set(offs+i,CMath::Sqr(a[offs+i][offs+i]));
              }
            else
              {
               //--- (I+1)x(I+1) matrix,
               //--- ( A11^H  A21^H )   ( A11      )   ( A11^H*A11+A21^H*A21  A21^H*A22 )
               //--- (              ) * (          ) = (                                )
               //--- (        A22^H )   ( A21  A22 )   ( A22^H*A21            A22^H*A22 )
               //--- A11 is IxI, A22 is 1x1.
               i1_=offs;
               for(i_=0;i_<=i-1;i_++)
                  tmp[i_]=a[offs+i][i_+i1_];
               for(j=0;j<i;j++)
                 {
                  v=a[offs+i][offs+j];
                  i1_=-offs;
                  for(i_=offs;i_<=offs+j;i_++)
                     a[offs+j].Set(i_,a[offs+j][i_]+v*tmp[i_+i1_]);
                 }
               //--- change values
               v=a[offs+i][offs+i];
               for(i_=offs;i_<offs+i;i_++)
                  a[offs+i].Set(i_,v*a[offs+i][i_]);
               a[offs+i].Set(offs+i,CMath::Sqr(a[offs+i][offs+i]));
              }
           }
        }
      //--- exit the function
      return;
     }
//--- Recursive code: triangular factor inversion merged with
//--- UU' or L'L multiplication
   CAblas::AblasSplitLength(a,n,n1,n2);
//--- form off-diagonal block of trangular inverse
   if(isupper)
     {
      for(i=0;i<n1;i++)
        {
         for(i_=offs+n1;i_<=offs+n-1;i_++)
            a[offs+i].Set(i_,-1*a[offs+i][i_]);
        }
      //--- function call
      CAblas::RMatrixLeftTrsM(n1,n2,a,offs,offs,isupper,false,0,a,offs,offs+n1);
      //--- function call
      CAblas::RMatrixRightTrsM(n1,n2,a,offs+n1,offs+n1,isupper,false,0,a,offs,offs+n1);
     }
   else
     {
      for(i=0;i<n2;i++)
        {
         for(i_=offs;i_<=offs+n1-1;i_++)
            a[offs+n1+i].Set(i_,-1*a[offs+n1+i][i_]);
        }
      //--- function call
      CAblas::RMatrixRightTrsM(n2,n1,a,offs,offs,isupper,false,0,a,offs+n1,offs);
      //--- function call
      CAblas::RMatrixLeftTrsM(n2,n1,a,offs+n1,offs+n1,isupper,false,0,a,offs+n1,offs);
     }
//--- invert first diagonal block
   SPDMatrixCholeskyInverseRec(a,offs,n1,isupper,tmp);
//--- update first diagonal block with off-diagonal block,
//--- update off-diagonal block
   if(isupper)
     {
      //--- function call
      CAblas::RMatrixSyrk(n1,n2,1.0,a,offs,offs+n1,0,1.0,a,offs,offs,isupper);
      //--- function call
      CAblas::RMatrixRightTrsM(n1,n2,a,offs+n1,offs+n1,isupper,false,1,a,offs,offs+n1);
     }
   else
     {
      //--- function call
      CAblas::RMatrixSyrk(n1,n2,1.0,a,offs+n1,offs,1,1.0,a,offs,offs,isupper);
      //--- function call
      CAblas::RMatrixLeftTrsM(n2,n1,a,offs+n1,offs+n1,isupper,false,1,a,offs+n1,offs);
     }
//--- invert second diagonal block
   SPDMatrixCholeskyInverseRec(a,offs+n1,n2,isupper,tmp);
  }
//+------------------------------------------------------------------+
//| Recursive subroutine for HPD inversion.                          |
//+------------------------------------------------------------------+
static void CMatInv::HPDMatrixCholeskyInverseRec(CMatrixComplex &a,const int offs,
                                                 const int n,const bool isupper,
                                                 al_complex &tmp[])
  {
//--- create variables
   al_complex _One=-1.0;
   int        i=0;
   int        j=0;
   al_complex v=0;
   int        n1=0;
   int        n2=0;
   int        info2=0;
   int        i_=0;
   int        i1_=0;
//--- object of class
   CMatInvReport rep2;
//--- check
   if(n<1)
      return;
//--- base case
   if(n<=CAblas::AblasComplexBlockSize())
     {
      CMatrixTrInverseRec(a,offs,n,isupper,false,tmp,info2,rep2);
      //--- check
      if(isupper)
        {
         //--- Compute the product U * U'.
         //--- NOTE: we never assume that diagonal of U is real
         for(i=0;i<n;i++)
           {
            //--- check
            if(i==0)
              {
               //--- 1x1 matrix
               a[offs+i].Set(offs+i,CMath::Sqr(a[offs+i][offs+i].re)+CMath::Sqr(a[offs+i][offs+i].im));
              }
            else
              {
               //--- (I+1)x(I+1) matrix,
               //--- ( A11  A12 )   ( A11^H        )   ( A11*A11^H+A12*A12^H  A12*A22^H )
               //--- (          ) * (              ) = (                                )
               //--- (      A22 )   ( A12^H  A22^H )   ( A22*A12^H            A22*A22^H )
               //--- A11 is IxI, A22 is 1x1.
               i1_=offs;
               for(i_=0;i_<=i-1;i_++)
                  tmp[i_]=CMath::Conj(a[i_+i1_][offs+i]);
               for(j=0;j<i;j++)
                 {
                  v=a[offs+j][offs+i];
                  i1_=-offs;
                  for(i_=offs+j;i_<offs+i;i_++)
                     a[offs+j].Set(i_,a[offs+j][i_]+v*tmp[i_+i1_]);
                 }
               //--- change values
               v=CMath::Conj(a[offs+i][offs+i]);
               for(i_=offs;i_<offs+i;i_++)
                  a[i_].Set(offs+i,v*a[i_][offs+i]);
               a[offs+i].Set(offs+i,CMath::Sqr(a[offs+i][offs+i].re)+CMath::Sqr(a[offs+i][offs+i].im));
              }
           }
        }
      else
        {
         //--- Compute the product L' * L
         //--- NOTE: we never assume that diagonal of L is reall
         for(i=0;i<n;i++)
           {
            //--- check
            if(i==0)
              {
               //--- 1x1 matrix
               a[offs+i].Set(offs+i,CMath::Sqr(a[offs+i][offs+i].re)+CMath::Sqr(a[offs+i][offs+i].im));
              }
            else
              {
               //--- (I+1)x(I+1) matrix,
               //--- ( A11^H  A21^H )   ( A11      )   ( A11^H*A11+A21^H*A21  A21^H*A22 )
               //--- (              ) * (          ) = (                                )
               //--- (        A22^H )   ( A21  A22 )   ( A22^H*A21            A22^H*A22 )
               //--- A11 is IxI, A22 is 1x1.
               i1_=offs;
               for(i_=0;i_<=i-1;i_++)
                  tmp[i_]=a[offs+i][i_+i1_];
               for(j=0;j<i;j++)
                 {
                  v=CMath::Conj(a[offs+i][offs+j]);
                  i1_=-offs;
                  for(i_=offs;i_<=offs+j;i_++)
                     a[offs+j].Set(i_,a[offs+j][i_]+v*tmp[i_+i1_]);
                 }
               v=CMath::Conj(a[offs+i][offs+i]);
               for(i_=offs;i_<offs+i;i_++)
                  a[offs+i].Set(i_,v*a[offs+i][i_]);
               a[offs+i].Set(offs+i,CMath::Sqr(a[offs+i][offs+i].re)+CMath::Sqr(a[offs+i][offs+i].im));
              }
           }
        }
      //--- exit the function
      return;
     }
//--- Recursive code: triangular factor inversion merged with
//--- UU' or L'L multiplication
   CAblas::AblasComplexSplitLength(a,n,n1,n2);
//--- form off-diagonal block of trangular inverse
   if(isupper)
     {
      for(i=0;i<n1;i++)
        {
         for(i_=offs+n1;i_<=offs+n-1;i_++)
            a[offs+i].Set(i_,_One*a[offs+i][i_]);
        }
      //--- function call
      CAblas::CMatrixLeftTrsM(n1,n2,a,offs,offs,isupper,false,0,a,offs,offs+n1);
      //--- function call
      CAblas::CMatrixRightTrsM(n1,n2,a,offs+n1,offs+n1,isupper,false,0,a,offs,offs+n1);
     }
   else
     {
      for(i=0;i<n2;i++)
        {
         for(i_=offs;i_<=offs+n1-1;i_++)
            a[offs+n1+i].Set(i_,_One*a[offs+n1+i][i_]);
        }
      //--- function call
      CAblas::CMatrixRightTrsM(n2,n1,a,offs,offs,isupper,false,0,a,offs+n1,offs);
      //--- function call
      CAblas::CMatrixLeftTrsM(n2,n1,a,offs+n1,offs+n1,isupper,false,0,a,offs+n1,offs);
     }
//--- invert first diagonal block
   HPDMatrixCholeskyInverseRec(a,offs,n1,isupper,tmp);
//--- update first diagonal block with off-diagonal block,
//--- update off-diagonal block
   if(isupper)
     {
      //--- function call
      CAblas::CMatrixSyrk(n1,n2,1.0,a,offs,offs+n1,0,1.0,a,offs,offs,isupper);
      //--- function call
      CAblas::CMatrixRightTrsM(n1,n2,a,offs+n1,offs+n1,isupper,false,2,a,offs,offs+n1);
     }
   else
     {
      //--- function call
      CAblas::CMatrixSyrk(n1,n2,1.0,a,offs+n1,offs,2,1.0,a,offs,offs,isupper);
      //--- function call
      CAblas::CMatrixLeftTrsM(n2,n1,a,offs+n1,offs+n1,isupper,false,2,a,offs+n1,offs);
     }
//--- invert second diagonal block
   HPDMatrixCholeskyInverseRec(a,offs+n1,n2,isupper,tmp);
  }
//+------------------------------------------------------------------+
//| Singular value decomposition of a bidiagonal matrix              |
//+------------------------------------------------------------------+
class CBdSingValueDecompose
  {
private:
   //--- private methods
   static bool       BidiagonalSVDDecompositionInternal(double &d[],double &ce[],const int n,const bool isupper,const bool isfractionalaccuracyrequired,CMatrixDouble &u,const int ustart,const int nru,CMatrixDouble &c,const int cstart,const int ncc,CMatrixDouble &vt,const int vstart,const int ncvt);
   static double     ExtSignBdSQR(const double a,const double b);
   static void       SVD2x2(const double f,const double g,const double h,double &ssmin,double &ssmax);
   static void       SVDV2x2(const double f,const double g,const double h,double &ssmin,double &ssmax,double &snr,double &csr,double &snl,double &csl);
public:
                     CBdSingValueDecompose(void);
                    ~CBdSingValueDecompose(void);
   //--- public methods
   static bool       RMatrixBdSVD(double &d[],double &ce[],const int n,const bool isupper,const bool isfractionalaccuracyrequired,CMatrixDouble &u,const int nru,CMatrixDouble &c,const int ncc,CMatrixDouble &vt,const int ncvt);
   static bool       BidiagonalSVDDecomposition(double &d[],double &ce[],const int n,const bool isupper,const bool isfractionalaccuracyrequired,CMatrixDouble &u,const int nru,CMatrixDouble &c,const int ncc,CMatrixDouble &vt,const int ncvt);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CBdSingValueDecompose::CBdSingValueDecompose(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CBdSingValueDecompose::~CBdSingValueDecompose(void)
  {

  }
//+------------------------------------------------------------------+
//| Singular value decomposition of a bidiagonal matrix (extended    |
//| algorithm)                                                       |
//| The algorithm performs the singular value decomposition of a     |
//| bidiagonal matrix B (upper or lower) representing it as          |
//| B = Q*S*P^T, where Q and P - orthogonal matrices, S - diagonal   |
//| matrix with non-negative elements on the main diagonal, in       |
//| descending order.                                                |
//| The algorithm finds singular values. In addition, the algorithm  |
//| can calculate matrices Q and P (more precisely, not the matrices,|
//| but their product with given matrices U and VT - U*Q and         |
//| (P^T)*VT)). Of course, matrices U and VT can be of any type,     |
//| including identity. Furthermore, the algorithm can calculate Q'*C|
//| (this product is calculated more effectively than U*Q, because   |
//| this calculation operates with rows instead  of matrix columns). |
//| The feature of the algorithm is its ability to find all singular |
//| values including those which are arbitrarily close to 0 with     |
//| relative accuracy close to  machine precision. If the parameter  |
//| IsFractionalAccuracyRequired is set to True, all singular values |
//| will have high relative accuracy close to machine precision. If  |
//| the parameter is set to False, only the biggest singular value   |
//| will have relative accuracy close to machine precision. The      |
//| absolute error of other singular values is equal to the absolute |
//| error of the biggest singular value.                             |
//| Input parameters:                                                |
//|     D       -   main diagonal of matrix B.                       |
//|                 Array whose index ranges within [0..N-1].        |
//|     E       -   superdiagonal (or subdiagonal) of matrix B.      |
//|                 Array whose index ranges within [0..N-2].        |
//|     N       -   size of matrix B.                                |
//|     IsUpper -   True, if the matrix is upper bidiagonal.         |
//|     IsFractionalAccuracyRequired -                               |
//|                 accuracy to search singular values with.         |
//|     U       -   matrix to be multiplied by Q.                    |
//|                 Array whose indexes range within                 |
//|                 [0..NRU-1, 0..N-1].                              |
//|                 The matrix can be bigger, in that case only the  |
//|                 submatrix [0..NRU-1, 0..N-1] will be multiplied  |
//|                 by Q.                                            |
//|     NRU     -   number of rows in matrix U.                      |
//|     C       -   matrix to be multiplied by Q'.                   |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..NCC-1].                              |
//|                 The matrix can be bigger, in that case only the  |
//|                 submatrix [0..N-1, 0..NCC-1] will be multiplied  |
//|                 by Q'.                                           |
//|     NCC     -   number of columns in matrix C.                   |
//|     VT      -   matrix to be multiplied by P^T.                  |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..NCVT-1].                             |
//|                 The matrix can be bigger, in that case only the  |
//|                 submatrix [0..N-1, 0..NCVT-1] will be multiplied |
//|                 by P^T.                                          |
//|     NCVT    -   number of columns in matrix VT.                  |
//| Output parameters:                                               |
//|     D       -   singular values of matrix B in descending order. |
//|     U       -   if NRU>0, contains matrix U*Q.                   |
//|     VT      -   if NCVT>0, contains matrix (P^T)*VT.             |
//|     C       -   if NCC>0, contains matrix Q'*C.                  |
//| Result:                                                          |
//|     True, if the algorithm has converged.                        |
//|     False, if the algorithm hasn't converged (rare case).        |
//| Additional information:                                          |
//|     The type of convergence is controlled by the internal        |
//|     parameter TOL. If the parameter is greater than 0, the       |
//|     singular values will have relative accuracy TOL. If TOL<0,   |
//|     the singular values will have absolute accuracy              |
//|     ABS(TOL)*norm(B). By default, |TOL| falls within the range of|
//|     10*Epsilon and 100*Epsilon, where Epsilon is the machine     |
//|     precision. It is not recommended to use TOL less than        |
//|     10*Epsilon since this will considerably slow down the        |
//|     algorithm and may not lead to error decreasing.              |
//| History:                                                         |
//|     * 31 March, 2007.                                            |
//|         changed MAXITR from 6 to 12.                             |
//|   -- LAPACK routine (version 3.0) --                             |
//|      Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., |
//|      Courant Institute, Argonne National Lab, and Rice University|
//|      October 31, 1999.                                           |
//+------------------------------------------------------------------+
static bool CBdSingValueDecompose::RMatrixBdSVD(double &d[],double &ce[],
                                                const int n,const bool isupper,
                                                const bool isfractionalaccuracyrequired,
                                                CMatrixDouble &u,const int nru,
                                                CMatrixDouble &c,const int ncc,
                                                CMatrixDouble &vt,const int ncvt)
  {
//--- create variables
   bool result;
//--- create arrays
   double d1[];
   double e1[];
   int i_=0;
   int i1_=0;
//--- create copy
   double e[];
   ArrayCopy(e,ce);
//--- allocation
   ArrayResizeAL(d1,n+1);
//--- change values
   i1_=-1;
   for(i_=1;i_<=n;i_++)
      d1[i_]=d[i_+i1_];
//--- check
   if(n>1)
     {
      //--- allocation
      ArrayResizeAL(e1,n);
      //--- change values
      i1_=-1;
      for(i_=1;i_<n;i_++)
         e1[i_]=e[i_+i1_];
     }
//--- get result
   result=BidiagonalSVDDecompositionInternal(d1,e1,n,isupper,isfractionalaccuracyrequired,u,0,nru,c,0,ncc,vt,0,ncvt);
//--- change values
   i1_=1;
   for(i_=0;i_<n;i_++)
      d[i_]=d1[i_+i1_];
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Singular value decomposition of a bidiagonal matrix              |
//+------------------------------------------------------------------+
static bool CBdSingValueDecompose::BidiagonalSVDDecomposition(double &d[],double &ce[],
                                                              const int n,const bool isupper,
                                                              const bool isfractionalaccuracyrequired,
                                                              CMatrixDouble &u,const int nru,
                                                              CMatrixDouble &c,const int ncc,
                                                              CMatrixDouble &vt,const int ncvt)
  {
//--- create copy
   double e[];
   ArrayCopy(e,ce);
//--- return result
   return(BidiagonalSVDDecompositionInternal(d,e,n,isupper,isfractionalaccuracyrequired,u,1,nru,c,1,ncc,vt,1,ncvt));
  }
//+------------------------------------------------------------------+
//| Internal working subroutine for bidiagonal decomposition         |
//+------------------------------------------------------------------+
static bool CBdSingValueDecompose::BidiagonalSVDDecompositionInternal(double &d[],double &ce[],
                                                                      const int n,const bool isupper,
                                                                      const bool isfractionalaccuracyrequired,
                                                                      CMatrixDouble &u,const int ustart,
                                                                      const int nru,CMatrixDouble &c,
                                                                      const int cstart,const int ncc,
                                                                      CMatrixDouble &vt,const int vstart,
                                                                      const int ncvt)
  {
//--- create variables
   bool   result;
   int    i=0;
   int    idir=0;
   int    isub=0;
   int    iter=0;
   int    j=0;
   int    ll=0;
   int    lll=0;
   int    m=0;
   int    maxit=0;
   int    oldll=0;
   int    oldm=0;
   double abse=0;
   double abss=0;
   double cosl=0;
   double cosr=0;
   double cs=0;
   double eps=0;
   double f=0;
   double g=0;
   double h=0;
   double mu=0;
   double oldcs=0;
   double oldsn=0;
   double r=0;
   double shift=0;
   double sigmn=0;
   double sigmx=0;
   double sinl=0;
   double sinr=0;
   double sll=0;
   double smax=0;
   double smin=0;
   double sminl=0;
   double sminlo=0;
   double sminoa=0;
   double sn=0;
   double thresh=0;
   double tol=0;
   double tolmul=0;
   double unfl=0;
   int    maxitr=0;
   bool   matrixsplitflag;
   bool   iterflag;
   bool   rightside;
   bool   fwddir;
   double tmp=0;
   int    mm1=0;
   int    mm0=0;
   bool   bchangedir;
   int    uend=0;
   int    cend=0;
   int    vend=0;
   int    i_=0;
//--- create arrays
   double work0[];
   double work1[];
   double work2[];
   double work3[];
   double utemp[];
   double vttemp[];
   double ctemp[];
   double etemp[];
//--- create copy
   double e[];
   ArrayCopy(e,ce);
//--- initialization
   result=true;
//--- check
   if(n==0)
      return(true);
//--- check
   if(n==1)
     {
      //--- check
      if(d[1]<0.0)
        {
         d[1]=-d[1];
         //--- check
         if(ncvt>0)
            for(i_=vstart;i_<=vstart+ncvt-1;i_++)
               vt[vstart].Set(i_,-1*vt[vstart][i_]);
        }
      //--- return result
      return(result);
     }
//--- initialization
   ll=0;
   oldsn=0;
//--- allocation
   ArrayResizeAL(work0,n);
   ArrayResizeAL(work1,n);
   ArrayResizeAL(work2,n);
   ArrayResizeAL(work3,n);
   uend=ustart+(int)MathMax(nru-1,0);
   vend=vstart+(int)MathMax(ncvt-1,0);
   cend=cstart+(int)MathMax(ncc-1,0);
   ArrayResizeAL(utemp,uend+1);
   ArrayResizeAL(vttemp,vend+1);
   ArrayResizeAL(ctemp,cend+1);
//--- initialization
   maxitr=12;
   rightside=true;
   fwddir=true;
//--- resize E from N-1 to N
   ArrayResizeAL(etemp,n+1);
   for(i=1;i<n;i++)
      etemp[i]=e[i];
   ArrayResizeAL(e,n+1);
   for(i=1;i<n;i++)
      e[i]=etemp[i];
   e[n]=0;
   idir=0;
//--- initialization
   eps=CMath::m_machineepsilon;
   unfl=CMath::m_minrealnumber;
//--- If matrix lower bidiagonal, rotate to be upper bidiagonal
//--- by applying Givens rotations on the left
   if(!isupper)
     {
      for(i=1;i<n;i++)
        {
         //--- function call
         CRotations::GenerateRotation(d[i],e[i],cs,sn,r);
         d[i]=r;
         e[i]=sn*d[i+1];
         d[i+1]=cs*d[i+1];
         work0[i]=cs;
         work1[i]=sn;
        }
      //--- Update singular vectors if desired
      if(nru>0)
        {
         //--- function call
         CRotations::ApplyRotationsFromTheRight(fwddir,ustart,uend,1+ustart-1,n+ustart-1,work0,work1,u,utemp);
        }
      //--- check
      if(ncc>0)
        {
         //--- function call
         CRotations::ApplyRotationsFromTheLeft(fwddir,1+cstart-1,n+cstart-1,cstart,cend,work0,work1,c,ctemp);
        }
     }
//--- Compute singular values to relative accuracy TOL
//--- (By setting TOL to be negative, algorithm will compute
//--- singular values to absolute accuracy ABS(TOL)*norm(input matrix))
   tolmul=MathMax(10,MathMin(100,MathPow(eps,-0.125)));
   tol=tolmul*eps;
//--- check
   if(!isfractionalaccuracyrequired)
      tol=-tol;
//--- Compute approximate maximum, minimum singular values
   smax=0;
   for(i=1;i<=n;i++)
      smax=MathMax(smax,MathAbs(d[i]));
   for(i=1;i<n;i++)
      smax=MathMax(smax,MathAbs(e[i]));
   sminl=0;
//--- check
   if(tol>=0.0)
     {
      //--- Relative accuracy desired
      sminoa=MathAbs(d[1]);
      //--- check
      if(sminoa!=0.0)
        {
         mu=sminoa;
         for(i=2;i<=n;i++)
           {
            mu=MathAbs(d[i])*(mu/(mu+MathAbs(e[i-1])));
            sminoa=MathMin(sminoa,mu);
            //--- check
            if(sminoa==0.0)
               break;
           }
        }
      //--- change values
      sminoa=sminoa/MathSqrt(n);
      thresh=MathMax(tol*sminoa,maxitr*n*n*unfl);
     }
   else
     {
      //--- Absolute accuracy desired
      thresh=MathMax(MathAbs(tol)*smax,maxitr*n*n*unfl);
     }
//--- Prepare for main iteration loop for the singular values
//--- (MAXIT is the maximum number of passes through the inner
//--- loop permitted before nonconvergence signalled.)
   maxit=maxitr*n*n;
   iter=0;
   oldll=-1;
   oldm=-1;
//--- M points to last element of unconverged part of matrix
   m=n;
//--- Begin main iteration loop
   while(true)
     {
      //--- Check for convergence or exceeding iteration count 
      if(m<=1)
         break;
      //--- check
      if(iter>maxit)
         return(false);
      //--- Find diagonal block of matrix to work on
      if(tol<0.0 && MathAbs(d[m])<=thresh)
         d[m]=0;
      //--- change values
      smax=MathAbs(d[m]);
      smin=smax;
      matrixsplitflag=false;
      for(lll=1;lll<=m-1;lll++)
        {
         ll=m-lll;
         abss=MathAbs(d[ll]);
         abse=MathAbs(e[ll]);
         //--- check
         if(tol<0.0 && abss<=thresh)
            d[ll]=0;
         //--- check
         if(abse<=thresh)
           {
            matrixsplitflag=true;
            break;
           }
         //--- change values
         smin=MathMin(smin,abss);
         smax=MathMax(smax,MathMax(abss,abse));
        }
      //--- check
      if(!matrixsplitflag)
         ll=0;
      else
        {
         //--- Matrix splits since E(LL) = 0
         e[ll]=0;
         //--- check
         if(ll==m-1)
           {
            //--- Convergence of bottom singular value, return to top of loop
            m=m-1;
            continue;
           }
        }
      ll=ll+1;
      //--- E(LL) through E(M-1) are nonzero, E(LL-1) is zero
      if(ll==m-1)
        {
         //--- 2 by 2 block, handle separately
         SVDV2x2(d[m-1],e[m-1],d[m],sigmn,sigmx,sinr,cosr,sinl,cosl);
         d[m-1]=sigmx;
         e[m-1]=0;
         d[m]=sigmn;
         //--- Compute singular vectors, if desired
         if(ncvt>0)
           {
            mm0=m+(vstart-1);
            mm1=m-1+(vstart-1);
            //--- swap
            for(i_=vstart;i_<=vend;i_++)
               vttemp[i_]=cosr*vt[mm1][i_];
            for(i_=vstart;i_<=vend;i_++)
               vttemp[i_]=vttemp[i_]+sinr*vt[mm0][i_];
            for(i_=vstart;i_<=vend;i_++)
               vt[mm0].Set(i_,cosr*vt[mm0][i_]);
            for(i_=vstart;i_<=vend;i_++)
               vt[mm0].Set(i_,vt[mm0][i_]-sinr*vt[mm1][i_]);
            for(i_=vstart;i_<=vend;i_++)
               vt[mm1].Set(i_,vttemp[i_]);
           }
         //--- check
         if(nru>0)
           {
            mm0=m+ustart-1;
            mm1=m-1+ustart-1;
            //--- swap
            for(i_=ustart;i_<=uend;i_++)
               utemp[i_]=cosl*u[i_][mm1];
            for(i_=ustart;i_<=uend;i_++)
               utemp[i_]=utemp[i_]+sinl*u[i_][mm0];
            for(i_=ustart;i_<=uend;i_++)
               u[i_].Set(mm0,cosl*u[i_][mm0]);
            for(i_=ustart;i_<=uend;i_++)
               u[i_].Set(mm0,u[i_][mm0]-sinl*u[i_][mm1]);
            for(i_=ustart;i_<=uend;i_++)
               u[i_].Set(mm1,utemp[i_]);
           }
         //--- check
         if(ncc>0)
           {
            mm0=m+cstart-1;
            mm1=m-1+cstart-1;
            //--- swap
            for(i_=cstart;i_<=cend;i_++)
               ctemp[i_]=cosl*c[mm1][i_];
            for(i_=cstart;i_<=cend;i_++)
               ctemp[i_]=ctemp[i_]+sinl*c[mm0][i_];
            for(i_=cstart;i_<=cend;i_++)
               c[mm0].Set(i_,cosl*c[mm0][i_]);
            for(i_=cstart;i_<=cend;i_++)
               c[mm0].Set(i_,c[mm0][i_]-sinl*c[mm1][i_]);
            for(i_=cstart;i_<=cend;i_++)
               c[mm1].Set(i_,ctemp[i_]);
           }
         m=m-2;
         continue;
        }
      //--- If working on new submatrix, choose shift direction
      //--- (from larger end diagonal element towards smaller)
      //--- Previously was
      //---     "if (LL>OLDM) or (M<OLDLL) then"
      //--- fixed thanks to Michael Rolle < m@rolle.name >
      //--- Very strange that LAPACK still contains it.
      bchangedir=false;
      //--- check
      if(idir==1 && MathAbs(d[ll])<1.0E-3*MathAbs(d[m]))
         bchangedir=true;
      //--- check
      if(idir==2 && MathAbs(d[m])<1.0E-3*MathAbs(d[ll]))
         bchangedir=true;
      //--- check
      if(ll!=oldll || m!=oldm || bchangedir)
        {
         //--- check
         if(MathAbs(d[ll])>=MathAbs(d[m]))
           {
            //--- Chase bulge from top(big end) to bottom(small end)
            idir=1;
           }
         else
           {
            //--- Chase bulge from bottom (big end) to top (small end)
            idir=2;
           }
        }
      //--- Apply convergence tests
      if(idir==1)
        {
         //--- Run convergence test in forward direction
         //--- First apply standard test to bottom of matrix
         if(MathAbs(e[m-1])<=MathAbs(tol)*MathAbs(d[m]) || (tol<0.0 && MathAbs(e[m-1])<=thresh))
           {
            e[m-1]=0;
            continue;
           }
         //--- check
         if(tol>=0.0)
           {
            //--- If relative accuracy desired,
            //--- apply convergence criterion forward
            mu=MathAbs(d[ll]);
            sminl=mu;
            iterflag=false;
            for(lll=ll;lll<=m-1;lll++)
              {
               //--- check
               if(MathAbs(e[lll])<=tol*mu)
                 {
                  e[lll]=0;
                  iterflag=true;
                  //--- break the cycle
                  break;
                 }
               //--- change values
               sminlo=sminl;
               mu=MathAbs(d[lll+1])*(mu/(mu+MathAbs(e[lll])));
               sminl=MathMin(sminl,mu);
              }
            //--- check
            if(iterflag)
               continue;
           }
        }
      else
        {
         //--- Run convergence test in backward direction
         //--- First apply standard test to top of matrix
         if(MathAbs(e[ll])<=MathAbs(tol)*MathAbs(d[ll]) || (tol<0.0 && MathAbs(e[ll])<=thresh))
           {
            e[ll]=0;
            continue;
           }
         //--- check
         if(tol>=0.0)
           {
            //--- If relative accuracy desired,
            //--- apply convergence criterion backward
            mu=MathAbs(d[m]);
            sminl=mu;
            iterflag=false;
            for(lll=m-1;lll>=ll;lll--)
              {
               //--- check
               if(MathAbs(e[lll])<=(double)(tol*mu))
                 {
                  e[lll]=0;
                  iterflag=true;
                  //--- break the cycle
                  break;
                 }
               sminlo=sminl;
               mu=MathAbs(d[lll])*(mu/(mu+MathAbs(e[lll])));
               sminl=MathMin(sminl,mu);
              }
            //--- check
            if(iterflag)
               continue;
           }
        }
      //--- change values
      oldll=ll;
      oldm=m;
      //--- Compute shift.  First, test if shifting would ruin relative
      //--- accuracy, and if so set the shift to zero.
      if(tol>=0.0 && n*tol*(sminl/smax)<=MathMax(eps,0.01*tol))
        {
         //--- Use a zero shift to avoid loss of relative accuracy
         shift=0;
        }
      else
        {
         //--- Compute the shift from 2-by-2 block at end of matrix
         if(idir==1)
           {
            sll=MathAbs(d[ll]);
            SVD2x2(d[m-1],e[m-1],d[m],shift,r);
           }
         else
           {
            sll=MathAbs(d[m]);
            SVD2x2(d[ll],e[ll],d[ll+1],shift,r);
           }
         //--- Test if shift negligible, and if so set to zero
         if(sll>0.0)
           {
            //--- check
            if(CMath::Sqr(shift/sll)<eps)
               shift=0;
           }
        }
      //--- Increment iteration count
      iter=iter+m-ll;
      //--- If SHIFT = 0, do simplified QR iteration
      if(shift==0.0)
        {
         //--- check
         if(idir==1)
           {
            //--- Chase bulge from top to bottom
            //--- Save cosines and sines for later singular vector updates
            cs=1;
            oldcs=1;
            for(i=ll;i<m;i++)
              {
               //--- function call
               CRotations::GenerateRotation(d[i]*cs,e[i],cs,sn,r);
               //--- check
               if(i>ll)
                  e[i-1]=oldsn*r;
               //--- function call
               CRotations::GenerateRotation(oldcs*r,d[i+1]*sn,oldcs,oldsn,tmp);
               //--- change values
               d[i]=tmp;
               work0[i-ll+1]=cs;
               work1[i-ll+1]=sn;
               work2[i-ll+1]=oldcs;
               work3[i-ll+1]=oldsn;
              }
            //--- change values
            h=d[m]*cs;
            d[m]=h*oldcs;
            e[m-1]=h*oldsn;
            //--- Update singular vectors
            if(ncvt>0)
              {
               //--- function call
               CRotations::ApplyRotationsFromTheLeft(fwddir,ll+vstart-1,m+vstart-1,vstart,vend,work0,work1,vt,vttemp);
              }
            //--- check
            if(nru>0)
              {
               //--- function call
               CRotations::ApplyRotationsFromTheRight(fwddir,ustart,uend,ll+ustart-1,m+ustart-1,work2,work3,u,utemp);
              }
            //--- check
            if(ncc>0)
              {
               //--- function call
               CRotations::ApplyRotationsFromTheLeft(fwddir,ll+cstart-1,m+cstart-1,cstart,cend,work2,work3,c,ctemp);
              }
            //--- Check for convergence
            if(MathAbs(e[m-1])<=thresh)
               e[m-1]=0;
           }
         else
           {
            //--- Chase bulge from bottom to top
            //--- Save cosines and sines for later singular vector updates
            cs=1;
            oldcs=1;
            for(i=m;i>=ll+1;i--)
              {
               //--- function call
               CRotations::GenerateRotation(d[i]*cs,e[i-1],cs,sn,r);
               //--- check
               if(i<m)
                  e[i]=oldsn*r;
               //--- function call
               CRotations::GenerateRotation(oldcs*r,d[i-1]*sn,oldcs,oldsn,tmp);
               //--- change values
               d[i]=tmp;
               work0[i-ll]=cs;
               work1[i-ll]=-sn;
               work2[i-ll]=oldcs;
               work3[i-ll]=-oldsn;
              }
            //--- change values
            h=d[ll]*cs;
            d[ll]=h*oldcs;
            e[ll]=h*oldsn;
            //--- Update singular vectors
            if(ncvt>0)
              {
               //--- function call
               CRotations::ApplyRotationsFromTheLeft(!fwddir,ll+vstart-1,m+vstart-1,vstart,vend,work2,work3,vt,vttemp);
              }
            //--- check
            if(nru>0)
              {
               //--- function call
               CRotations::ApplyRotationsFromTheRight(!fwddir,ustart,uend,ll+ustart-1,m+ustart-1,work0,work1,u,utemp);
              }
            //--- check
            if(ncc>0)
              {
               //--- function call
               CRotations::ApplyRotationsFromTheLeft(!fwddir,ll+cstart-1,m+cstart-1,cstart,cend,work0,work1,c,ctemp);
              }
            //--- Check for convergence
            if(MathAbs(e[ll])<=thresh)
               e[ll]=0;
           }
        }
      else
        {
         //--- Use nonzero shift
         if(idir==1)
           {
            //--- Chase bulge from top to bottom
            //--- Save cosines and sines for later singular vector updates
            f=(MathAbs(d[ll])-shift)*(ExtSignBdSQR(1,d[ll])+shift/d[ll]);
            g=e[ll];
            for(i=ll;i<m;i++)
              {
               //--- function call
               CRotations::GenerateRotation(f,g,cosr,sinr,r);
               //--- check
               if(i>ll)
                  e[i-1]=r;
               //--- change values
               f=cosr*d[i]+sinr*e[i];
               e[i]=cosr*e[i]-sinr*d[i];
               g=sinr*d[i+1];
               d[i+1]=cosr*d[i+1];
               //--- function call
               CRotations::GenerateRotation(f,g,cosl,sinl,r);
               //--- change values
               d[i]=r;
               f=cosl*e[i]+sinl*d[i+1];
               d[i+1]=cosl*d[i+1]-sinl*e[i];
               //--- check
               if(i<m-1)
                 {
                  g=sinl*e[i+1];
                  e[i+1]=cosl*e[i+1];
                 }
               //--- change values
               work0[i-ll+1]=cosr;
               work1[i-ll+1]=sinr;
               work2[i-ll+1]=cosl;
               work3[i-ll+1]=sinl;
              }
            e[m-1]=f;
            //--- Update singular vectors
            if(ncvt>0)
              {
               //--- function call
               CRotations::ApplyRotationsFromTheLeft(fwddir,ll+vstart-1,m+vstart-1,vstart,vend,work0,work1,vt,vttemp);
              }
            //--- check
            if(nru>0)
              {
               //--- function call
               CRotations::ApplyRotationsFromTheRight(fwddir,ustart,uend,ll+ustart-1,m+ustart-1,work2,work3,u,utemp);
              }
            //--- check
            if(ncc>0)
              {
               //--- function call
               CRotations::ApplyRotationsFromTheLeft(fwddir,ll+cstart-1,m+cstart-1,cstart,cend,work2,work3,c,ctemp);
              }
            //--- Check for convergence
            if(MathAbs(e[m-1])<=thresh)
               e[m-1]=0;
           }
         else
           {
            //--- Chase bulge from bottom to top
            //--- Save cosines and sines for later singular vector updates
            f=(MathAbs(d[m])-shift)*(ExtSignBdSQR(1,d[m])+shift/d[m]);
            g=e[m-1];
            for(i=m;i>=ll+1;i--)
              {
               //--- function call
               CRotations::GenerateRotation(f,g,cosr,sinr,r);
               //--- check
               if(i<m)
                  e[i]=r;
               //--- change values
               f=cosr*d[i]+sinr*e[i-1];
               e[i-1]=cosr*e[i-1]-sinr*d[i];
               g=sinr*d[i-1];
               d[i-1]=cosr*d[i-1];
               //--- function call
               CRotations::GenerateRotation(f,g,cosl,sinl,r);
               //--- change values
               d[i]=r;
               f=cosl*e[i-1]+sinl*d[i-1];
               d[i-1]=cosl*d[i-1]-sinl*e[i-1];
               //--- check
               if(i>ll+1)
                 {
                  g=sinl*e[i-2];
                  e[i-2]=cosl*e[i-2];
                 }
               //--- change values
               work0[i-ll]=cosr;
               work1[i-ll]=-sinr;
               work2[i-ll]=cosl;
               work3[i-ll]=-sinl;
              }
            e[ll]=f;
            //--- Check for convergence
            if(MathAbs(e[ll])<=thresh)
               e[ll]=0;
            //--- Update singular vectors if desired
            if(ncvt>0)
              {
               //--- function call
               CRotations::ApplyRotationsFromTheLeft(!fwddir,ll+vstart-1,m+vstart-1,vstart,vend,work2,work3,vt,vttemp);
              }
            //--- check
            if(nru>0)
              {
               //--- function call
               CRotations::ApplyRotationsFromTheRight(!fwddir,ustart,uend,ll+ustart-1,m+ustart-1,work0,work1,u,utemp);
              }
            //--- check
            if(ncc>0)
              {
               //--- function call
               CRotations::ApplyRotationsFromTheLeft(!fwddir,ll+cstart-1,m+cstart-1,cstart,cend,work0,work1,c,ctemp);
              }
           }
        }
      //--- QR iteration finished, go back and check convergence
      continue;
     }
//--- All singular values converged, so make them positive
   for(i=1;i<=n;i++)
     {
      //--- check
      if(d[i]<0.0)
        {
         d[i]=-d[i];
         //--- Change sign of singular vectors, if desired
         if(ncvt>0)
           {
            for(i_=vstart;i_<=vend;i_++)
               vt[i+vstart-1].Set(i_,-1*vt[i+vstart-1][i_]);
           }
        }
     }
//--- Sort the singular values into decreasing order (insertion sort on
//--- singular values, but only one transposition per singular vector)
   for(i=1;i<n;i++)
     {
      //--- Scan for smallest D(I)
      isub=1;
      smin=d[1];
      for(j=2;j<=n+1-i;j++)
        {
         //--- check
         if(d[j]<=smin)
           {
            isub=j;
            smin=d[j];
           }
        }
      //--- check
      if(isub!=n+1-i)
        {
         //--- Swap singular values and vectors
         d[isub]=d[n+1-i];
         d[n+1-i]=smin;
         //--- check
         if(ncvt>0)
           {
            j=n+1-i;
            //--- swap
            for(i_=vstart;i_<=vend;i_++)
               vttemp[i_]=vt[isub+vstart-1][i_];
            for(i_=vstart;i_<=vend;i_++)
               vt[isub+vstart-1].Set(i_,vt[j+vstart-1][i_]);
            for(i_=vstart;i_<=vend;i_++)
               vt[j+vstart-1].Set(i_,vttemp[i_]);
           }
         if(nru>0)
           {
            j=n+1-i;
            //--- swap
            for(i_=ustart;i_<=uend;i_++)
               utemp[i_]=u[i_][isub+ustart-1];
            for(i_=ustart;i_<=uend;i_++)
               u[i_].Set(isub+ustart-1,u[i_][j+ustart-1]);
            for(i_=ustart;i_<=uend;i_++)
               u[i_].Set(j+ustart-1,utemp[i_]);
           }
         //--- check
         if(ncc>0)
           {
            j=n+1-i;
            //--- swap
            for(i_=cstart;i_<=cend;i_++)
               ctemp[i_]=c[isub+cstart-1][i_];
            for(i_=cstart;i_<=cend;i_++)
               c[isub+cstart-1].Set(i_,c[j+cstart-1][i_]);
            for(i_=cstart;i_<=cend;i_++)
               c[j+cstart-1].Set(i_,ctemp[i_]);
           }
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static double CBdSingValueDecompose::ExtSignBdSQR(const double a,
                                                  const double b)
  {
//--- create a variable
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
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CBdSingValueDecompose::SVD2x2(const double f,const double g,
                                          const double h,double &ssmin,
                                          double &ssmax)
  {
//--- create variables
   double aas=0;
   double at=0;
   double au=0;
   double c=0;
   double fa=0;
   double fhmn=0;
   double fhmx=0;
   double ga=0;
   double ha=0;
//--- initialization
   ssmin=0;
   ssmax=0;
   fa=MathAbs(f);
   ga=MathAbs(g);
   ha=MathAbs(h);
   fhmn=MathMin(fa,ha);
   fhmx=MathMax(fa,ha);
//--- check
   if(fhmn==0.0)
     {
      ssmin=0;
      //--- check
      if(fhmx==0.0)
         ssmax=ga;
      else
         ssmax=MathMax(fhmx,ga)*MathSqrt(1+CMath::Sqr(MathMin(fhmx,ga)/MathMax(fhmx,ga)));
     }
   else
     {
      //--- check
      if(ga<fhmx)
        {
         //--- change values
         aas=1+fhmn/fhmx;
         at=(fhmx-fhmn)/fhmx;
         au=CMath::Sqr(ga/fhmx);
         c=2/(MathSqrt(aas*aas+au)+MathSqrt(at*at+au));
         ssmin=fhmn*c;
         ssmax=fhmx/c;
        }
      else
        {
         au=fhmx/ga;
         //--- check
         if(au==0.0)
           {
            //--- Avoid possible harmful underflow if exponent range
            //--- asymmetric (true SSMIN may not underflow even if
            //--- AU underflows)
            ssmin=fhmn*fhmx/ga;
            ssmax=ga;
           }
         else
           {
            //--- change values
            aas=1+fhmn/fhmx;
            at=(fhmx-fhmn)/fhmx;
            c=1/(MathSqrt(1+CMath::Sqr(aas*au))+MathSqrt(1+CMath::Sqr(at*au)));
            ssmin=fhmn*c*au;
            ssmin=ssmin+ssmin;
            ssmax=ga/(c+c);
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CBdSingValueDecompose::SVDV2x2(const double f,const double g,
                                           const double h,double &ssmin,
                                           double &ssmax,double &snr,
                                           double &csr,double &snl,
                                           double &csl)
  {
//--- create variables
   bool   gasmal;
   bool   swp;
   int    pmax=0;
   double a=0;
   double clt=0;
   double crt=0;
   double d=0;
   double fa=0;
   double ft=0;
   double ga=0;
   double gt=0;
   double ha=0;
   double ht=0;
   double l=0;
   double m=0;
   double mm=0;
   double r=0;
   double s=0;
   double slt=0;
   double srt=0;
   double t=0;
   double temp=0;
   double tsign=0;
   double tt=0;
   double v=0;
//--- initialization
   ssmin=0;
   ssmax=0;
   snr=0;
   csr=0;
   snl=0;
   csl=0;
   ft=f;
   fa=MathAbs(ft);
   ht=h;
   ha=MathAbs(h);
   clt=0;
   crt=0;
   slt=0;
   srt=0;
   tsign=0;
//--- PMAX points to the maximum absolute element of matrix
//--- PMAX = 1 if F largest in absolute values
//--- PMAX = 2 if G largest in absolute values
//--- PMAX = 3 if H largest in absolute values
   pmax=1;
   swp=ha>fa;
//--- check
   if(swp)
     {
      //--- Now FA .ge. HA
      pmax=3;
      temp=ft;
      ft=ht;
      ht=temp;
      temp=fa;
      fa=ha;
      ha=temp;
     }
   gt=g;
   ga=MathAbs(gt);
//--- check
   if(ga==0.0)
     {
      //--- Diagonal matrix
      ssmin=ha;
      ssmax=fa;
      clt=1;
      crt=1;
      slt=0;
      srt=0;
     }
   else
     {
      gasmal=true;
      //--- check
      if(ga>fa)
        {
         pmax=2;
         //--- check
         if(fa/ga<CMath::m_machineepsilon)
           {
            //--- Case of very large GA
            gasmal=false;
            ssmax=ga;
            //--- check
            if(ha>1.0)
              {
               v=ga/ha;
               ssmin=fa/v;
              }
            else
              {
               v=fa/ga;
               ssmin=v*ha;
              }
            //--- change values
            clt=1;
            slt=ht/gt;
            srt=1;
            crt=ft/gt;
           }
        }
      //--- check
      if(gasmal)
        {
         //--- Normal case
         d=fa-ha;
         //--- check
         if(d==fa)
            l=1;
         else
            l=d/fa;
         //--- change values
         m=gt/ft;
         t=2-l;
         mm=m*m;
         tt=t*t;
         s=MathSqrt(tt+mm);
         //--- check
         if(l==0.0)
            r=MathAbs(m);
         else
            r=MathSqrt(l*l+mm);
         //--- change values
         a=0.5*(s+r);
         ssmin=ha/a;
         ssmax=fa*a;
         //--- check
         if(mm==0.0)
           {
            //--- Note that M is very tiny
            if(l==0.0)
               t=ExtSignBdSQR(2,ft)*ExtSignBdSQR(1,gt);
            else
               t=gt/ExtSignBdSQR(d,ft)+m/t;
           }
         else
            t=(m/(s+t)+m/(r+l))*(1+a);
         //--- change values
         l=MathSqrt(t*t+4);
         crt=2/l;
         srt=t/l;
         clt=(crt+srt*m)/a;
         v=ht/ft;
         slt=v*srt/a;
        }
     }
//--- check
   if(swp)
     {
      csl=srt;
      snl=crt;
      csr=slt;
      snr=clt;
     }
   else
     {
      csl=clt;
      snl=slt;
      csr=crt;
      snr=srt;
     }
//--- Correct signs of SSMAX and SSMIN
   if(pmax==1)
      tsign=ExtSignBdSQR(1,csr)*ExtSignBdSQR(1,csl)*ExtSignBdSQR(1,f);
//--- check
   if(pmax==2)
      tsign=ExtSignBdSQR(1,snr)*ExtSignBdSQR(1,csl)*ExtSignBdSQR(1,g);
//--- check
   if(pmax==3)
      tsign=ExtSignBdSQR(1,snr)*ExtSignBdSQR(1,snl)*ExtSignBdSQR(1,h);
//--- get result
   ssmax=ExtSignBdSQR(ssmax,tsign);
   ssmin=ExtSignBdSQR(ssmin,tsign*ExtSignBdSQR(1,f)*ExtSignBdSQR(1,h));
  }
//+------------------------------------------------------------------+
//| Singular value decomposition                                     |
//+------------------------------------------------------------------+
class CSingValueDecompose
  {
public:
   //--- constructor, destructor
                     CSingValueDecompose(void);
                    ~CSingValueDecompose(void);
   //--- method
   static bool       RMatrixSVD(CMatrixDouble &ca,const int m,const int n,const int uneeded,const int vtneeded,const int additionalmemory,double &w[],CMatrixDouble &u,CMatrixDouble &vt);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSingValueDecompose::CSingValueDecompose(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSingValueDecompose::~CSingValueDecompose(void)
  {

  }
//+------------------------------------------------------------------+
//| Singular value decomposition of a rectangular matrix.            |
//| The algorithm calculates the singular value decomposition of a   |
//| matrix of size MxN: A = U * S * V^T                              |
//| The algorithm finds the singular values and, optionally, matrices|
//| U and V^T. The algorithm can find both first min(M,N) columns of |
//| matrix U and rows of matrix V^T (singular vectors), and matrices |
//| U and V^T wholly (of sizes MxM and NxN respectively).            |
//| Take into account that the subroutine does not return matrix V   |
//| but V^T.                                                         |
//| Input parameters:                                                |
//|     A           -   matrix to be decomposed.                     |
//|                     Array whose indexes range within             |
//|                     [0..M-1, 0..N-1].                            |
//|     M           -   number of rows in matrix A.                  |
//|     N           -   number of columns in matrix A.               |
//|     UNeeded     -   0, 1 or 2. See the description of the        |
//|                     parameter U.                                 |
//|     VTNeeded    -   0, 1 or 2. See the description of the        |
//|                     parameter VT.                                |
//|     AdditionalMemory -                                           |
//|                     If the parameter:                            |
//|                      * equals 0, the algorithm doesn?t use       |
//|                        additional memory (lower requirements,    |
//|                        lower performance).                       |
//|                      * equals 1, the algorithm uses additional   |
//|                        memory of size min(M,N)*min(M,N) of real  |
//|                        numbers. It often speeds up the algorithm.|
//|                      * equals 2, the algorithm uses additional   |
//|                        memory of size M*min(M,N) of real numbers.|
//|                        It allows to get a maximum performance.   |
//|                     The recommended value of the parameter is 2. |
//| Output parameters:                                               |
//|     W           -   contains singular values in descending order.|
//|     U           -   if UNeeded=0, U isn't changed, the left      |
//|                     singular vectors are not calculated.         |
//|                     if Uneeded=1, U contains left singular       |
//|                     vectors (first min(M,N) columns of matrix U).|
//|                     Array whose indexes range within             |
//|                     [0..M-1, 0..Min(M,N)-1]. if UNeeded=2, U     |
//|                     contains matrix U wholly. Array whose indexes|
//|                     range within [0..M-1, 0..M-1].               |
//|     VT          -   if VTNeeded=0, VT isn?t changed, the right   |
//|                     singular vectors are not calculated.         |
//|                     if VTNeeded=1, VT contains right singular    |
//|                     vectors (first min(M,N) rows of matrix V^T). |
//|                     Array whose indexes range within             |
//|                     [0..min(M,N)-1, 0..N-1]. if VTNeeded=2, VT   |
//|                     contains matrix V^T wholly. Array whose      |
//|                     indexes range within [0..N-1, 0..N-1].       |
//+------------------------------------------------------------------+
static bool CSingValueDecompose::RMatrixSVD(CMatrixDouble &ca,const int m,
                                            const int n,const int uneeded,
                                            const int vtneeded,
                                            const int additionalmemory,
                                            double &w[],CMatrixDouble &u,
                                            CMatrixDouble &vt)
  {
//--- create variables
   bool result;
   bool isupper;
   int  minmn=0;
   int  ncu=0;
   int  nrvt=0;
   int  nru=0;
   int  ncvt=0;
   int  i=0;
   int  j=0;
//--- create arrays
   double tauq[];
   double taup[];
   double tau[];
   double e[];
   double work[];
//--- create matrix
   CMatrixDouble t2;
//--- create copy
   CMatrixDouble a;
   a=ca;
//--- initialization
   result=true;
//--- check
   if(m==0 || n==0)
      return(true);
//--- check
   if(!CAp::Assert(uneeded>=0 && uneeded<=2,__FUNCTION__+": wrong parameters!"))
      return(false);
//--- check
   if(!CAp::Assert(vtneeded>=0 && vtneeded<=2,__FUNCTION__+": wrong parameters!"))
      return(false);
//--- check
   if(!CAp::Assert(additionalmemory>=0 && additionalmemory<=2,__FUNCTION__+": wrong parameters!"))
      return(false);
//--- initialization
   minmn=MathMin(m,n);
   ArrayResizeAL(w,minmn+1);
   ncu=0;
   nru=0;
//--- check
   if(uneeded==1)
     {
      nru=m;
      ncu=minmn;
      u.Resize(nru,ncu);
     }
//--- check
   if(uneeded==2)
     {
      nru=m;
      ncu=m;
      u.Resize(nru,ncu);
     }
   nrvt=0;
   ncvt=0;
//--- check
   if(vtneeded==1)
     {
      nrvt=minmn;
      ncvt=n;
      vt.Resize(nrvt,ncvt);
     }
//--- check
   if(vtneeded==2)
     {
      nrvt=n;
      ncvt=n;
      vt.Resize(nrvt,ncvt);
     }
//--- M much larger than N
//--- Use bidiagonal reduction with QR-decomposition
   if((double)(m)>(double)(1.6*n))
     {
      //--- check
      if(uneeded==0)
        {
         //-- No left singular vectors to be computed
         COrtFac::RMatrixQR(a,m,n,tau);
         for(i=0;i<n;i++)
           {
            for(j=0;j<i;j++)
               a[i].Set(j,0);
           }
         //--- function call
         COrtFac::RMatrixBD(a,n,n,tauq,taup);
         //--- function call
         COrtFac::RMatrixBDUnpackPT(a,n,n,taup,nrvt,vt);
         //--- function call
         COrtFac::RMatrixBDUnpackDiagonals(a,n,n,isupper,w,e);
         //--- return result
         return(CBdSingValueDecompose::RMatrixBdSVD(w,e,n,isupper,false,u,0,a,0,vt,ncvt));
        }
      else
        {
         //--- Left singular vectors (may be full matrix U) to be computed
         COrtFac::RMatrixQR(a,m,n,tau);
         //--- function call
         COrtFac::RMatrixQRUnpackQ(a,m,n,tau,ncu,u);
         for(i=0;i<n;i++)
           {
            for(j=0;j<i;j++)
               a[i].Set(j,0);
           }
         //--- function call
         COrtFac::RMatrixBD(a,n,n,tauq,taup);
         //--- function call
         COrtFac::RMatrixBDUnpackPT(a,n,n,taup,nrvt,vt);
         //--- function call
         COrtFac::RMatrixBDUnpackDiagonals(a,n,n,isupper,w,e);
         //--- check
         if(additionalmemory<1)
           {
            //--- No additional memory can be used
            COrtFac::RMatrixBDMultiplyByQ(a,n,n,tauq,u,m,n,true,false);
            //--- get result
            result=CBdSingValueDecompose::RMatrixBdSVD(w,e,n,isupper,false,u,m,a,0,vt,ncvt);
           }
         else
           {
            //--- Large U. Transforming intermediate matrix T2
            ArrayResizeAL(work,MathMax(m,n)+1);
            //--- function call
            COrtFac::RMatrixBDUnpackQ(a,n,n,tauq,n,t2);
            //--- function call
            CBlas::CopyMatrix(u,0,m-1,0,n-1,a,0,m-1,0,n-1);
            //--- function call
            CBlas::InplaceTranspose(t2,0,n-1,0,n-1,work);
            //--- get result
            result=CBdSingValueDecompose::RMatrixBdSVD(w,e,n,isupper,false,u,0,t2,n,vt,ncvt);
            //--- function call
            CBlas::MatrixMatrixMultiply(a,0,m-1,0,n-1,false,t2,0,n-1,0,n-1,true,1.0,u,0,m-1,0,n-1,0.0,work);
           }
         //--- return result
         return(result);
        }
     }
//--- N much larger than M
//--- Use bidiagonal reduction with LQ-decomposition
   if((double)(n)>(double)(1.6*m))
     {
      //--- check
      if(vtneeded==0)
        {
         //--- No right singular vectors to be computed
         COrtFac::RMatrixLQ(a,m,n,tau);
         for(i=0;i<=m-1;i++)
           {
            for(j=i+1;j<=m-1;j++)
               a[i].Set(j,0);
           }
         //--- function call
         COrtFac::RMatrixBD(a,m,m,tauq,taup);
         //--- function call
         COrtFac::RMatrixBDUnpackQ(a,m,m,tauq,ncu,u);
         //--- function call
         COrtFac::RMatrixBDUnpackDiagonals(a,m,m,isupper,w,e);
         ArrayResizeAL(work,m+1);
         //--- function call
         CBlas::InplaceTranspose(u,0,nru-1,0,ncu-1,work);
         //--- get result
         result=CBdSingValueDecompose::RMatrixBdSVD(w,e,m,isupper,false,a,0,u,nru,vt,0);
         //--- function call
         CBlas::InplaceTranspose(u,0,nru-1,0,ncu-1,work);
         //--- return result
         return(result);
        }
      else
        {
         //--- Right singular vectors (may be full matrix VT) to be computed
         COrtFac::RMatrixLQ(a,m,n,tau);
         //--- function call
         COrtFac::RMatrixLQUnpackQ(a,m,n,tau,nrvt,vt);
         for(i=0;i<=m-1;i++)
           {
            for(j=i+1;j<=m-1;j++)
               a[i].Set(j,0);
           }
         //--- function call
         COrtFac::RMatrixBD(a,m,m,tauq,taup);
         //--- function call
         COrtFac::RMatrixBDUnpackQ(a,m,m,tauq,ncu,u);
         //--- function call
         COrtFac::RMatrixBDUnpackDiagonals(a,m,m,isupper,w,e);
         ArrayResizeAL(work,MathMax(m,n)+1);
         //--- function call
         CBlas::InplaceTranspose(u,0,nru-1,0,ncu-1,work);
         //--- check
         if(additionalmemory<1)
           {
            //--- No additional memory available
            COrtFac::RMatrixBDMultiplyByP(a,m,m,taup,vt,m,n,false,true);
            //--- get result
            result=CBdSingValueDecompose::RMatrixBdSVD(w,e,m,isupper,false,a,0,u,nru,vt,n);
           }
         else
           {
            //--- Large VT. Transforming intermediate matrix T2
            COrtFac::RMatrixBDUnpackPT(a,m,m,taup,m,t2);
            //--- get result
            result=CBdSingValueDecompose::RMatrixBdSVD(w,e,m,isupper,false,a,0,u,nru,t2,m);
            //--- function call
            CBlas::CopyMatrix(vt,0,m-1,0,n-1,a,0,m-1,0,n-1);
            //--- function call
            CBlas::MatrixMatrixMultiply(t2,0,m-1,0,m-1,false,a,0,m-1,0,n-1,false,1.0,vt,0,m-1,0,n-1,0.0,work);
           }
         //--- function call
         CBlas::InplaceTranspose(u,0,nru-1,0,ncu-1,work);
         //--- return result
         return(result);
        }
     }
//--- M<=N
//--- We can use inplace transposition of U to get rid of columnwise operations
   if(m<=n)
     {
      //--- function call
      COrtFac::RMatrixBD(a,m,n,tauq,taup);
      //--- function call
      COrtFac::RMatrixBDUnpackQ(a,m,n,tauq,ncu,u);
      //--- function call
      COrtFac::RMatrixBDUnpackPT(a,m,n,taup,nrvt,vt);
      //--- function call
      COrtFac::RMatrixBDUnpackDiagonals(a,m,n,isupper,w,e);
      ArrayResizeAL(work,m+1);
      //--- function call
      CBlas::InplaceTranspose(u,0,nru-1,0,ncu-1,work);
      //--- get result
      result=CBdSingValueDecompose::RMatrixBdSVD(w,e,minmn,isupper,false,a,0,u,nru,vt,ncvt);
      //--- function call
      CBlas::InplaceTranspose(u,0,nru-1,0,ncu-1,work);
      //--- return result
      return(result);
     }
//--- Simple bidiagonal reduction
   COrtFac::RMatrixBD(a,m,n,tauq,taup);
//--- function call
   COrtFac::RMatrixBDUnpackQ(a,m,n,tauq,ncu,u);
//--- function call
   COrtFac::RMatrixBDUnpackPT(a,m,n,taup,nrvt,vt);
//--- function call
   COrtFac::RMatrixBDUnpackDiagonals(a,m,n,isupper,w,e);
//--- check
   if(additionalmemory<2 || uneeded==0)
     {
      //--- We cant use additional memory or there is no need in such operations
      result=CBdSingValueDecompose::RMatrixBdSVD(w,e,minmn,isupper,false,u,nru,a,0,vt,ncvt);
     }
   else
     {
      //--- We can use additional memory
      t2.Resize(minmn,m);
      //--- function call
      CBlas::CopyAndTranspose(u,0,m-1,0,minmn-1,t2,0,minmn-1,0,m-1);
      //--- get result
      result=CBdSingValueDecompose::RMatrixBdSVD(w,e,minmn,isupper,false,u,0,t2,m,vt,ncvt);
      //--- function call
      CBlas::CopyAndTranspose(t2,0,minmn-1,0,m-1,u,0,m-1,0,minmn-1);
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Structure which stores state of linear CG solver between         |
//| subsequent calls of FBLSCgIteration(). Initialized with          |
//| FBLSCGCreate().                                                  |
//| USAGE:                                                           |
//| 1. call to FBLSCGCreate()                                        |
//| 2. F:=FBLSCgIteration(State)                                     |
//| 3. if F is False, iterations are over                            |
//| 4. otherwise, fill State.AX with A*x, State.XAX with x'*A*x      |
//| 5. goto 2                                                        |
//| If you want to rerminate iterations, pass zero or negative value |
//| to XAX.                                                          |
//| FIELDS:                                                          |
//|     E1      -   2-norm of residual at the start                  |
//|     E2      -   2-norm of residual at the end                    |
//|     X       -   on return from FBLSCgIteration() it contains     |
//|                 vector for matrix-vector product                 |
//|     AX      -   must be filled with A*x if FBLSCgIteration()     |
//|                 returned True                                    |
//|     XAX     -   must be filled with x'*A*x                       |
//|     XK      -   contains result (if FBLSCgIteration() returned   | 
//|                 False)                                           |
//| Other fields are private and should not be used by outsiders.    |
//+------------------------------------------------------------------+
class CFblsLinCgState
  {
public:
   //--- variables
   double            m_e1;
   double            m_e2;
   double            m_x[];
   double            m_ax[];
   double            m_xax;
   double            m_xk[];
   int               m_n;
   double            m_rk[];
   double            m_rk1[];
   double            m_xk1[];
   double            m_pk[];
   double            m_pk1[];
   double            m_b[];
   RCommState        m_rstate;
   double            m_tmp2[];
   //--- constructor, destructor
                     CFblsLinCgState(void);
                    ~CFblsLinCgState(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CFblsLinCgState::CFblsLinCgState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CFblsLinCgState::~CFblsLinCgState(void)
  {

  }
//+------------------------------------------------------------------+
//| Fast basic linear solutions                                      |
//+------------------------------------------------------------------+
class CFbls
  {
private:
   //--- auxiliary functions for FblsCGiteration
   static void       Func_lbl_rcomm(CFblsLinCgState &state,int n,int k,double rk2,double rk12,double pap,double s,double betak,double v1,double v2);
   static bool       Func_lbl_3(CFblsLinCgState &state,int &n,int &k,double &rk2,double &rk12,double &pap,double &s,double &betak,double &v1,double &v2);
   static bool       Func_lbl_5(CFblsLinCgState &state,int &n,int &k,double &rk2,double &rk12,double &pap,double &s,double &betak,double &v1,double &v2);
public:
   //--- constructor, destructor
                     CFbls(void);
                    ~CFbls(void);
   //--- methods
   static void       FblsCholeskySolve(CMatrixDouble &cha,const double sqrtscalea,const int n,const bool isupper,double &xb[],double &tmp[]);
   static void       FblsSolveCGx(CMatrixDouble &a,const int m,const int n,const double alpha,const double &b[],double &x[],double &buf[]);
   static void       FblsCGCreate(double &x[],double &b[],const int n,CFblsLinCgState &state);
   static bool       FblsCGIteration(CFblsLinCgState &state);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CFbls::CFbls(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CFbls::~CFbls(void)
  {

  }
//+------------------------------------------------------------------+
//| Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.             |
//| This subroutine assumes that:                                    |
//| * A*ScaleA is well scaled                                        |
//| * A is well-conditioned, so no zero divisions or overflow may    | 
//|   occur                                                          |
//| INPUT PARAMETERS:                                                |
//|     CHA     -   Cholesky decomposition of A                      |
//|     SqrtScaleA- square root of scale factor ScaleA               |
//|     N       -   matrix size                                      |
//|     IsUpper -   storage type                                     |
//|     XB      -   right part                                       |
//|     Tmp     -   buffer; function automatically allocates it, if  |
//|                 it is too small. It can be reused if function is |
//|                 called several times.                            |
//| OUTPUT PARAMETERS:                                               |
//|     XB      -   solution                                         |
//| NOTES: no assertion or tests are done during algorithm operation |
//+------------------------------------------------------------------+
static void CFbls::FblsCholeskySolve(CMatrixDouble &cha,const double sqrtscalea,
                                     const int n,const bool isupper,double &xb[],
                                     double &tmp[])
  {
//--- create variables
   int    i=0;
   double v=0;
   int    i_=0;
//--- check
   if(CAp::Len(tmp)<n)
      ArrayResizeAL(tmp,n);
//--- A = L*L' or A=U'*U
   if(isupper)
     {
      //--- Solve U'*y=b first.
      for(i=0;i<n;i++)
        {
         xb[i]=xb[i]/(sqrtscalea*cha[i][i]);
         //--- check
         if(i<n-1)
           {
            v=xb[i];
            for(i_=i+1;i_<n;i_++)
               tmp[i_]=sqrtscalea*cha[i][i_];
            for(i_=i+1;i_<n;i_++)
               xb[i_]=xb[i_]-v*tmp[i_];
           }
        }
      //--- Solve U*x=y then.
      for(i=n-1;i>=0;i--)
        {
         //--- check
         if(i<n-1)
           {
            for(i_=i+1;i_<n;i_++)
               tmp[i_]=sqrtscalea*cha[i][i_];
            v=0.0;
            for(i_=i+1;i_<n;i_++)
               v+=tmp[i_]*xb[i_];
            xb[i]=xb[i]-v;
           }
         xb[i]=xb[i]/(sqrtscalea*cha[i][i]);
        }
     }
   else
     {
      //--- Solve L*y=b first
      for(i=0;i<n;i++)
        {
         //--- check
         if(i>0)
           {
            for(i_=0;i_<i;i_++)
               tmp[i_]=sqrtscalea*cha[i][i_];
            v=0.0;
            for(i_=0;i_<i;i_++)
               v+=tmp[i_]*xb[i_];
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
            for(i_=0;i_<i;i_++)
               tmp[i_]=sqrtscalea*cha[i][i_];
            for(i_=0;i_<i;i_++)
               xb[i_]=xb[i_]-v*tmp[i_];
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Fast basic linear solver: linear SPD CG                          |
//| Solves (A^T*A + alpha*I)*x = b where:                            |
//| * A is MxN matrix                                                |
//| * alpha>0 is a scalar                                            |
//| * I is NxN identity matrix                                       |
//| * b is Nx1 vector                                                |
//| * X is Nx1 unknown vector.                                       |
//| N iterations of linear conjugate gradient are used to solve      |
//| problem.                                                         |
//| INPUT PARAMETERS:                                                |
//|     A   -   array[M,N], matrix                                   |
//|     M   -   number of rows                                       |
//|     N   -   number of unknowns                                   |
//|     B   -   array[N], right part                                 |
//|     X   -   initial approxumation, array[N]                      |
//|     Buf -   buffer; function automatically allocates it, if it   |
//|             is too small. It can be reused if function is called |
//|             several times with same M and N.                     |
//| OUTPUT PARAMETERS:                                               |
//|     X   -   improved solution                                    |
//| NOTES:                                                           |
//| *   solver checks quality of improved solution. If (because of   |
//|     problem condition number, numerical noise, etc.) new solution|
//|     is WORSE than original approximation, then original          |
//|     approximation is returned.                                   |
//| *   solver assumes that both A, B, Alpha are well scaled (i.e.   |
//|     they are less than sqrt(overflow) and greater than           |
//|     sqrt(underflow)).                                            |
//+------------------------------------------------------------------+
static void CFbls::FblsSolveCGx(CMatrixDouble &a,const int m,const int n,
                                const double alpha,const double &b[],
                                double &x[],double &buf[])
  {
//--- create variables
   int    k=0;
   int    offsrk=0;
   int    offsrk1=0;
   int    offsxk=0;
   int    offsxk1=0;
   int    offspk=0;
   int    offspk1=0;
   int    offstmp1=0;
   int    offstmp2=0;
   int    bs=0;
   double e1=0;
   double e2=0;
   double rk2=0;
   double rk12=0;
   double pap=0;
   double s=0;
   double betak=0;
   double v1=0;
   double v2=0;
   int    i_=0;
   int    i1_=0;
//--- Test for special case: B=0
   v1=0.0;
   for(i_=0;i_<n;i_++)
      v1+=b[i_]*b[i_];
//--- check
   if(v1==0.0)
     {
      for(k=0;k<n;k++)
         x[k]=0;
      //--- exit the function
      return;
     }
//--- Offsets inside Buf for:
//--- * R[K], R[K+1]
//--- * X[K], X[K+1]
//--- * P[K], P[K+1]
//--- * Tmp1 - array[M], Tmp2 - array[N]
   offsrk=0;
   offsrk1=offsrk+n;
   offsxk=offsrk1+n;
   offsxk1=offsxk+n;
   offspk=offsxk1+n;
   offspk1=offspk+n;
   offstmp1=offspk1+n;
   offstmp2=offstmp1+m;
   bs=offstmp2+n;
//--- check
   if(CAp::Len(buf)<bs)
      ArrayResizeAL(buf,bs);
//--- x(0)=x
   i1_=-offsxk;
   for(i_=offsxk;i_<=offsxk+n-1;i_++)
      buf[i_]=x[i_+i1_];
//--- r(0)=b-A*x(0)
//--- RK2=r(0)'*r(0)
   CAblas::RMatrixMVect(m,n,a,0,0,0,buf,offsxk,buf,offstmp1);
//--- function call
   CAblas::RMatrixMVect(n,m,a,0,0,1,buf,offstmp1,buf,offstmp2);
//--- change array
   i1_=offsxk-offstmp2;
   for(i_=offstmp2;i_<offstmp2+n;i_++)
      buf[i_]=buf[i_]+alpha*buf[i_+i1_];
//--- change array
   i1_=-offsrk;
   for(i_=offsrk;i_<offsrk+n;i_++)
      buf[i_]=b[i_+i1_];
//--- change array
   i1_=offstmp2-offsrk;
   for(i_=offsrk;i_<offsrk+n;i_++)
      buf[i_]=buf[i_]-buf[i_+i1_];
   rk2=0.0;
   for(i_=offsrk;i_<offsrk+n;i_++)
      rk2+=buf[i_]*buf[i_];
//--- change array
   i1_=offsrk-offspk;
   for(i_=offspk;i_<=offspk+n-1;i_++)
      buf[i_]=buf[i_+i1_];
   e1=MathSqrt(rk2);
//--- cycle
   for(k=0;k<n;k++)
     {
      //--- Calculate A*p(k) - store in Buf[OffsTmp2:OffsTmp2+N-1]
      //--- and p(k)'*A*p(k)  - store in PAP
      //--- If PAP=0, break (iteration is over)
      CAblas::RMatrixMVect(m,n,a,0,0,0,buf,offspk,buf,offstmp1);
      v1=0.0;
      for(i_=offstmp1;i_<=offstmp1+m-1;i_++)
         v1+=buf[i_]*buf[i_];
      v2=0.0;
      for(i_=offspk;i_<=offspk+n-1;i_++)
         v2+=buf[i_]*buf[i_];
      pap=v1+alpha*v2;
      //--- function call
      CAblas::RMatrixMVect(n,m,a,0,0,1,buf,offstmp1,buf,offstmp2);
      i1_=offspk-offstmp2;
      for(i_=offstmp2;i_<offstmp2+n;i_++)
         buf[i_]=buf[i_]+alpha*buf[i_+i1_];
      //--- check
      if(pap==0.0)
         break;
      //--- S=(r(k)'*r(k))/(p(k)'*A*p(k))
      s=rk2/pap;
      //--- x(k+1)=x(k) + S*p(k)
      i1_=offsxk-offsxk1;
      for(i_=offsxk1;i_<=offsxk1+n-1;i_++)
         buf[i_]=buf[i_+i1_];
      i1_=offspk-offsxk1;
      for(i_=offsxk1;i_<=offsxk1+n-1;i_++)
         buf[i_]=buf[i_]+s*buf[i_+i1_];
      //--- r(k+1)=r(k) - S*A*p(k)
      //--- RK12=r(k+1)'*r(k+1)
      //--- Break if r(k+1) small enough (when compared to r(k))
      i1_=offsrk-offsrk1;
      for(i_=offsrk1;i_<=offsrk1+n-1;i_++)
         buf[i_]=buf[i_+i1_];
      i1_=offstmp2-offsrk1;
      for(i_=offsrk1;i_<=offsrk1+n-1;i_++)
         buf[i_]=buf[i_]-s*buf[i_+i1_];
      rk12=0.0;
      for(i_=offsrk1;i_<=offsrk1+n-1;i_++)
         rk12+=buf[i_]*buf[i_];
      //--- check
      if(MathSqrt(rk12)<=100*CMath::m_machineepsilon*MathSqrt(rk2))
        {
         //--- X(k) = x(k+1) before exit -
         //--- - because we expect to find solution at x(k)
         i1_=offsxk1-offsxk;
         for(i_=offsxk;i_<=offsxk+n-1;i_++)
            buf[i_]=buf[i_+i1_];
         break;
        }
      //--- BetaK=RK12/RK2
      //--- p(k+1)=r(k+1)+betak*p(k)
      betak=rk12/rk2;
      i1_=offsrk1-offspk1;
      for(i_=offspk1;i_<=offspk1+n-1;i_++)
         buf[i_]=buf[i_+i1_];
      i1_=offspk-offspk1;
      for(i_=offspk1;i_<=offspk1+n-1;i_++)
         buf[i_]=buf[i_]+betak*buf[i_+i1_];
      //--- r(k) :=r(k+1)
      //--- x(k) :=x(k+1)
      //--- p(k) :=p(k+1)
      i1_=offsrk1-offsrk;
      for(i_=offsrk;i_<offsrk+n;i_++)
         buf[i_]=buf[i_+i1_];
      i1_=offsxk1-offsxk;
      for(i_=offsxk;i_<=offsxk+n-1;i_++)
         buf[i_]=buf[i_+i1_];
      i1_=offspk1-offspk;
      for(i_=offspk;i_<=offspk+n-1;i_++)
         buf[i_]=buf[i_+i1_];
      rk2=rk12;
     }
//--- Calculate E2
   CAblas::RMatrixMVect(m,n,a,0,0,0,buf,offsxk,buf,offstmp1);
//--- function call
   CAblas::RMatrixMVect(n,m,a,0,0,1,buf,offstmp1,buf,offstmp2);
   i1_=offsxk-offstmp2;
   for(i_=offstmp2;i_<offstmp2+n;i_++)
      buf[i_]=buf[i_]+alpha*buf[i_+i1_];
   i1_=-offsrk;
   for(i_=offsrk;i_<offsrk+n;i_++)
      buf[i_]=b[i_+i1_];
   i1_=offstmp2-offsrk;
   for(i_=offsrk;i_<offsrk+n;i_++)
      buf[i_]=buf[i_]-buf[i_+i1_];
   v1=0.0;
   for(i_=offsrk;i_<offsrk+n;i_++)
      v1+=buf[i_]*buf[i_];
   e2=MathSqrt(v1);
//--- Output result (if it was improved)
   if(e2<e1)
     {
      i1_=offsxk;
      for(i_=0;i_<n;i_++)
         x[i_]=buf[i_+i1_];
     }
  }
//+------------------------------------------------------------------+
//| Construction of linear conjugate gradient solver.                |
//| State parameter passed using "var" semantics (i.e. previous state|  
//| is NOT erased). When it is already initialized, we can reause    |
//| prevously allocated memory.                                      |
//| INPUT PARAMETERS:                                                |
//|     X       -   initial solution                                 |
//|     B       -   right part                                       |
//|     N       -   system size                                      |
//|     State   -   structure; may be preallocated, if we want to    |
//|                 reuse memory                                     |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure which is used by FBLSCGIteration() to  |
//|                 store algorithm state between subsequent calls.  |
//| NOTE: no error checking is done; caller must check all           |
//|       parameters, prevent overflows, and so on.                  |
//+------------------------------------------------------------------+
static void CFbls::FblsCGCreate(double &x[],double &b[],const int n,CFblsLinCgState &state)
  {
//--- create variables
   int i_=0;
//--- check
   if(CAp::Len(state.m_b)<n)
      ArrayResizeAL(state.m_b,n);
//--- check
   if(CAp::Len(state.m_rk)<n)
      ArrayResizeAL(state.m_rk,n);
//--- check
   if(CAp::Len(state.m_rk1)<n)
      ArrayResizeAL(state.m_rk1,n);
//--- check
   if(CAp::Len(state.m_xk)<n)
      ArrayResizeAL(state.m_xk,n);
//--- check
   if(CAp::Len(state.m_xk1)<n)
      ArrayResizeAL(state.m_xk1,n);
//--- check
   if(CAp::Len(state.m_pk)<n)
      ArrayResizeAL(state.m_pk,n);
//--- check
   if(CAp::Len(state.m_pk1)<n)
      ArrayResizeAL(state.m_pk1,n);
//--- check
   if(CAp::Len(state.m_tmp2)<n)
      ArrayResizeAL(state.m_tmp2,n);
//--- check
   if(CAp::Len(state.m_x)<n)
      ArrayResizeAL(state.m_x,n);
//--- check
   if(CAp::Len(state.m_ax)<n)
      ArrayResizeAL(state.m_ax,n);
//--- check
   state.m_n=n;
   for(i_=0;i_<n;i_++)
      state.m_xk[i_]=x[i_];
   for(i_=0;i_<n;i_++)
      state.m_b[i_]=b[i_];
//--- allocation
   ArrayResizeAL(state.m_rstate.ia,2);
   ArrayResizeAL(state.m_rstate.ra,7);
   state.m_rstate.stage=-1;
  }
//+------------------------------------------------------------------+
//| Linear CG solver, function relying on reverse communication to   |
//| calculate matrix-vector products.                                |
//| See comments for FBLSLinCGState structure for more info.         |
//+------------------------------------------------------------------+
static bool CFbls::FblsCGIteration(CFblsLinCgState &state)
  {
//--- create variables
   int    n=0;
   int    k=0;
   double rk2=0;
   double rk12=0;
   double pap=0;
   double s=0;
   double betak=0;
   double v1=0;
   double v2=0;
   int    i_=0;
//--- This code initializes locals by:
//--- * random values determined during code
//---   generation - on first subroutine call
//--- * values from previous call - on subsequent calls
   if(state.m_rstate.stage>=0)
     {
      //--- initialization
      n=state.m_rstate.ia[0];
      k=state.m_rstate.ia[1];
      rk2=state.m_rstate.ra[0];
      rk12=state.m_rstate.ra[1];
      pap=state.m_rstate.ra[2];
      s=state.m_rstate.ra[3];
      betak=state.m_rstate.ra[4];
      v1=state.m_rstate.ra[5];
      v2=state.m_rstate.ra[6];
     }
   else
     {
      //--- initialization
      n=-983;
      k=-989;
      rk2=-834;
      rk12=900;
      pap=-287;
      s=364;
      betak=214;
      v1=-338;
      v2=-686;
     }
//--- check
   if(state.m_rstate.stage==0)
     {
      //--- copy
      for(i_=0;i_<n;i_++)
         state.m_rk[i_]=state.m_b[i_];
      //--- calculation
      for(i_=0;i_<n;i_++)
         state.m_rk[i_]=state.m_rk[i_]-state.m_ax[i_];
      //--- change value
      rk2=0.0;
      for(i_=0;i_<n;i_++)
         rk2+=state.m_rk[i_]*state.m_rk[i_];
      //--- copy
      for(i_=0;i_<n;i_++)
         state.m_pk[i_]=state.m_rk[i_];
      state.m_e1=MathSqrt(rk2);
      //--- cycle
      k=0;
      //--- function call
      return(Func_lbl_3(state,n,k,rk2,rk12,pap,s,betak,v1,v2));
     }
//--- check
   if(state.m_rstate.stage==1)
     {
      //--- copy
      for(i_=0;i_<n;i_++)
         state.m_tmp2[i_]=state.m_ax[i_];
      pap=state.m_xax;
      //--- check
      if(!CMath::IsFinite(pap))
         return(Func_lbl_5(state,n,k,rk2,rk12,pap,s,betak,v1,v2));
      //--- check
      if(pap<=0.0)
         return(Func_lbl_5(state,n,k,rk2,rk12,pap,s,betak,v1,v2));
      //--- S=(r(k)'*r(k))/(p(k)'*A*p(k))
      s=rk2/pap;
      //--- x(k+1)=x(k) + S*p(k)
      for(i_=0;i_<n;i_++)
         state.m_xk1[i_]=state.m_xk[i_];
      for(i_=0;i_<n;i_++)
         state.m_xk1[i_]=state.m_xk1[i_]+s*state.m_pk[i_];
      //--- r(k+1)=r(k) - S*A*p(k)
      //--- RK12=r(k+1)'*r(k+1)
      //--- Break if r(k+1) small enough (when compared to r(k))
      for(i_=0;i_<n;i_++)
         state.m_rk1[i_]=state.m_rk[i_];
      for(i_=0;i_<n;i_++)
         state.m_rk1[i_]=state.m_rk1[i_]-s*state.m_tmp2[i_];
      rk12=0.0;
      for(i_=0;i_<n;i_++)
         rk12+=state.m_rk1[i_]*state.m_rk1[i_];
      //--- check
      if(MathSqrt(rk12)<=100*CMath::m_machineepsilon*state.m_e1)
        {
         //--- X(k) = x(k+1) before exit -
         //--- - because we expect to find solution at x(k)
         for(i_=0;i_<n;i_++)
            state.m_xk[i_]=state.m_xk1[i_];
         //--- function call
         return(Func_lbl_5(state,n,k,rk2,rk12,pap,s,betak,v1,v2));
        }
      //--- BetaK=RK12/RK2
      //--- p(k+1)=r(k+1)+betak*p(k)
      //--- NOTE: we expect that BetaK won't overflow because of
      //--- "Sqrt(RK12)<=100*MachineEpsilon*E1" test above.
      betak=rk12/rk2;
      for(i_=0;i_<n;i_++)
         state.m_pk1[i_]=state.m_rk1[i_];
      for(i_=0;i_<n;i_++)
         state.m_pk1[i_]=state.m_pk1[i_]+betak*state.m_pk[i_];
      //--- r(k) :=r(k+1)
      //--- x(k) :=x(k+1)
      //--- p(k) :=p(k+1)
      for(i_=0;i_<n;i_++)
         state.m_rk[i_]=state.m_rk1[i_];
      for(i_=0;i_<n;i_++)
         state.m_xk[i_]=state.m_xk1[i_];
      for(i_=0;i_<n;i_++)
         state.m_pk[i_]=state.m_pk1[i_];
      rk2=rk12;
      k=k+1;
      //--- function call
      return(Func_lbl_3(state,n,k,rk2,rk12,pap,s,betak,v1,v2));
     }
//--- check
   if(state.m_rstate.stage==2)
     {
      //--- copy
      for(i_=0;i_<n;i_++)
         state.m_rk[i_]=state.m_b[i_];
      //--- calculation
      for(i_=0;i_<n;i_++)
         state.m_rk[i_]=state.m_rk[i_]-state.m_ax[i_];
      //--- change value
      v1=0.0;
      for(i_=0;i_<n;i_++)
         v1+=state.m_rk[i_]*state.m_rk[i_];
      state.m_e2=MathSqrt(v1);
      //--- return result
      return(false);
     }
//--- Routine body
//--- prepare locals
   n=state.m_n;
//--- Test for special case: B=0
   v1=0.0;
   for(i_=0;i_<n;i_++)
      v1+=state.m_b[i_]*state.m_b[i_];
//--- check
   if(v1==0.0)
     {
      for(k=0;k<n;k++)
         state.m_xk[k]=0;
      //--- return result
      return(false);
     }
//--- r(0)=b-A*x(0)
//--- RK2=r(0)'*r(0)
   for(i_=0;i_<n;i_++)
      state.m_x[i_]=state.m_xk[i_];
   state.m_rstate.stage=0;
//--- Saving state
   Func_lbl_rcomm(state,n,k,rk2,rk12,pap,s,betak,v1,v2);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for FblsCGiteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static void CFbls::Func_lbl_rcomm(CFblsLinCgState &state,int n,int k,
                                  double rk2,double rk12,double pap,
                                  double s,double betak,double v1,double v2)
  {
//--- save
   state.m_rstate.ia[0]=n;
   state.m_rstate.ia[1]=k;
   state.m_rstate.ra[0]=rk2;
   state.m_rstate.ra[1]=rk12;
   state.m_rstate.ra[2]=pap;
   state.m_rstate.ra[3]=s;
   state.m_rstate.ra[4]=betak;
   state.m_rstate.ra[5]=v1;
   state.m_rstate.ra[6]=v2;
  }
//+------------------------------------------------------------------+
//| Auxiliary function for FblsCGiteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CFbls::Func_lbl_3(CFblsLinCgState &state,int &n,int &k,
                              double &rk2,double &rk12,double &pap,
                              double &s,double &betak,double &v1,double &v2)
  {
//--- check
   if(k>n-1)
      return(Func_lbl_5(state,n,k,rk2,rk12,pap,s,betak,v1,v2));
//--- Calculate A*p(k) - store in State.Tmp2
//--- and p(k)'*A*p(k)  - store in PAP
//--- If PAP=0,break (iteration is over)
   for(int i_=0;i_<n;i_++)
      state.m_x[i_]=state.m_pk[i_];
   state.m_rstate.stage=1;
//--- Saving state
   Func_lbl_rcomm(state,n,k,rk2,rk12,pap,s,betak,v1,v2);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Auxiliary function for FblsCGiteration. Is a product to get rid  |
//| of the operator unconditional jump goto.                         |
//+------------------------------------------------------------------+
static bool CFbls::Func_lbl_5(CFblsLinCgState &state,int &n,int &k,
                              double &rk2,double &rk12,double &pap,
                              double &s,double &betak,double &v1,double &v2)
  {
//--- calculation E2
   for(int i_=0;i_<n;i_++)
      state.m_x[i_]=state.m_xk[i_];
   state.m_rstate.stage=2;
//--- Saving state
   Func_lbl_rcomm(state,n,k,rk2,rk12,pap,s,betak,v1,v2);
//--- return result
   return(true);
  }
//+------------------------------------------------------------------+
//| Determinant calculation                                          |
//+------------------------------------------------------------------+
class CMatDet
  {
public:
                     CMatDet(void);
                    ~CMatDet(void);
   //--- public methods
   static double     RMatrixLUDet(CMatrixDouble &a,int &pivots[],const int n);
   static double     RMatrixDet(CMatrixDouble &ca,const int n);
   static double     SPDMatrixCholeskyDet(CMatrixDouble &a,const int n);
   static double     SPDMatrixDet(CMatrixDouble &ca,const int n,const bool isupper);
   static al_complex CMatrixLUDet(CMatrixComplex &a,int &pivots[],const int n);
   static al_complex CMatrixDet(CMatrixComplex &ca,const int n);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMatDet::CMatDet(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMatDet::~CMatDet(void)
  {

  }
//+------------------------------------------------------------------+
//| Determinant calculation of the matrix given by its LU            |
//| decomposition.                                                   |
//| Input parameters:                                                |
//|     A       -   LU decomposition of the matrix (output of        |
//|                 RMatrixLU subroutine).                           |
//|     Pivots  -   table of permutations which were made during     |
//|                 the LU decomposition.                            |
//|                 Output of RMatrixLU subroutine.                  |
//|     N       -   (optional) size of matrix A:                     |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, automatically determined from    |
//|                   matrix size (A must be square matrix)          |
//| Result: matrix determinant.                                      |
//+------------------------------------------------------------------+
static double CMatDet::RMatrixLUDet(CMatrixDouble &a,int &pivots[],const int n)
  {
//--- create variables
   double result=0;
   int    i=0;
   int    s=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Len(pivots)>=n,__FUNCTION__+": Pivots array is too short!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(a,n,n),__FUNCTION__+": A contains infinite or NaN values!"))
      return(EMPTY_VALUE);
//--- initialization
   result=1;
   s=1;
   for(i=0;i<n;i++)
     {
      result=result*a[i][i];
      //--- check
      if(pivots[i]!=i)
         s=-s;
     }
//--- return result
   return(result*s);
  }
//+------------------------------------------------------------------+
//| Calculation of the determinant of a general matrix               |
//| Input parameters:                                                |
//|     A       -   matrix, array[0..N-1, 0..N-1]                    |
//|     N       -   (optional) size of matrix A:                     |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, automatically determined from    |
//|                   matrix size (A must be square matrix)          |
//| Result: determinant of matrix A.                                 |
//+------------------------------------------------------------------+
static double CMatDet::RMatrixDet(CMatrixDouble &ca,const int n)
  {
//--- create variables
   double result=0;
   int    pivots[];
//--- create copy
   CMatrixDouble a;
   a=ca;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(a,n,n),__FUNCTION__+": A contains infinite or NaN values!"))
      return(EMPTY_VALUE);
//--- function call
   CTrFac::RMatrixLU(a,n,n,pivots);
//--- return result
   return(RMatrixLUDet(a,pivots,n));
  }
//+------------------------------------------------------------------+
//| Determinant calculation of the matrix given by its LU            |
//| decomposition.                                                   |
//| Input parameters:                                                |
//|     A       -   LU decomposition of the matrix (output of        |
//|                 RMatrixLU subroutine).                           |
//|     Pivots  -   table of permutations which were made during     |
//|                 the LU decomposition.                            |
//|                 Output of RMatrixLU subroutine.                  |
//|     N       -   (optional) size of matrix A:                     |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, automatically determined from    |
//|                   matrix size (A must be square matrix)          |
//| Result: matrix determinant.                                      |
//+------------------------------------------------------------------+
static al_complex CMatDet::CMatrixLUDet(CMatrixComplex &a,int &pivots[],const int n)
  {
//--- create variables
   al_complex result=0;
   int        i=0;
   int        s=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Len(pivots)>=n,__FUNCTION__+": Pivots array is too short!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CApServ::IsFiniteComplexMatrix(a,n,n),__FUNCTION__+": A contains infinite or NaN values!"))
      return(EMPTY_VALUE);
//--- initialization
   result=1;
   s=1;
   for(i=0;i<n;i++)
     {
      result=result*a[i][i];
      //--- check
      if(pivots[i]!=i)
         s=-s;
     }
//--- return result
   return(result*s);
  }
//+------------------------------------------------------------------+
//| Calculation of the determinant of a general matrix               |
//| Input parameters:                                                |
//|     A       -   matrix, array[0..N-1, 0..N-1]                    |
//|     N       -   (optional) size of matrix A:                     |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, automatically determined from    |
//|                   matrix size (A must be square matrix)          |
//| Result: determinant of matrix A.                                 |
//+------------------------------------------------------------------+
static al_complex CMatDet::CMatrixDet(CMatrixComplex &ca,const int n)
  {
//--- create variables
   al_complex result=0;
   int        pivots[];
//--- create copy
   CMatrixComplex a;
   a=ca;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CApServ::IsFiniteComplexMatrix(a,n,n),__FUNCTION__+": A contains infinite or NaN values!"))
      return(EMPTY_VALUE);
//--- function call
   CTrFac::CMatrixLU(a,n,n,pivots);
//--- return result
   return(CMatrixLUDet(a,pivots,n));
  }
//+------------------------------------------------------------------+
//| Determinant calculation of the matrix given by the Cholesky      |
//| decomposition.                                                   |
//| Input parameters:                                                |
//|     A       -   Cholesky decomposition,                          |
//|                 output of SMatrixCholesky subroutine.            |
//|     N       -   (optional) size of matrix A:                     |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, automatically determined from    |
//|                   matrix size (A must be square matrix)          |
//| As the determinant is equal to the product of squares of diagonal|
//| elements, it?s not necessary to specify which triangle - lower   |
//| or upper - the matrix is stored in.                              |
//| Result:                                                          |
//|     matrix determinant.                                          |
//+------------------------------------------------------------------+
static double CMatDet::SPDMatrixCholeskyDet(CMatrixDouble &a,const int n)
  {
//--- create variables
   double result=0;
   int    i=0;
   bool   f;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return(EMPTY_VALUE);
//--- initialization
   f=true;
   for(i=0;i<n;i++)
      f=f && CMath::IsFinite(a[i][i]);
//--- check
   if(!CAp::Assert(f,__FUNCTION__+": A contains infinite or NaN values!"))
      return(EMPTY_VALUE);
   result=1;
   for(i=0;i<n;i++)
      result=result*CMath::Sqr(a[i][i]);
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Determinant calculation of the symmetric positive definite       |
//| matrix.                                                          |
//| Input parameters:                                                |
//|     A       -   matrix. Array with elements [0..N-1, 0..N-1].    |
//|     N       -   (optional) size of matrix A:                     |
//|                 * if given, only principal NxN submatrix is      |
//|                   processed and overwritten. other elements are  |
//|                   unchanged.                                     |
//|                 * if not given, automatically determined from    |
//|                   matrix size (A must be square matrix)          |
//|     IsUpper -   (optional) storage type:                         |
//|                 * if True, symmetric matrix A is given by its    |
//|                   upper triangle, and the lower triangle isn?t   |
//|                   used/changed  by function                      |
//|                 * if False, symmetric matrix A is given by its   |
//|                   lower triangle, and the upper triangle isn?t   |
//|                   used/changed by function                       |
//|                 * if not given, both lower and upper triangles   |
//|                   must be filled.                                |
//| Result:                                                          |
//|     determinant of matrix A.                                     |
//|     If matrix A is not positive definite, exception is thrown.   |
//+------------------------------------------------------------------+
static double CMatDet::SPDMatrixDet(CMatrixDouble &ca,const int n,const bool isupper)
  {
//--- create variables
   double result=0;
   bool   b;
//--- create copy
   CMatrixDouble a;
   a=ca;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Rows(a)>=n,__FUNCTION__+": rows(A)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CAp::Cols(a)>=n,__FUNCTION__+": cols(A)<N!"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(CApServ::IsFiniteRTrMatrix(a,n,isupper),__FUNCTION__+": A contains infinite or NaN values!"))
      return(EMPTY_VALUE);
//--- function call
   b=CTrFac::SPDMatrixCholesky(a,n,isupper);
//--- check
   if(!CAp::Assert(b,__FUNCTION__+": A is not SPD!"))
      return(EMPTY_VALUE);
//--- return result
   return(SPDMatrixCholeskyDet(a,n));
  }
//+------------------------------------------------------------------+
//| Generalized symmetric eigensolver                                |
//+------------------------------------------------------------------+
class CSpdGEVD
  {
public:
   //--- constructor, destructor
                     CSpdGEVD(void);
                    ~CSpdGEVD(void);
   //--- methods
   static bool       SMatrixGEVD(CMatrixDouble &ca,const int n,const bool isuppera,CMatrixDouble &b,const bool isupperb,const int zneeded,const int problemtype,double &d[],CMatrixDouble &z);
   static bool       SMatrixGEVDReduce(CMatrixDouble &a,const int n,const bool isuppera,CMatrixDouble &b,const bool isupperb,const int problemtype,CMatrixDouble &r,bool &isupperr);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSpdGEVD::CSpdGEVD(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSpdGEVD::~CSpdGEVD(void)
  {

  }
//+------------------------------------------------------------------+
//| Algorithm for solving the following generalized symmetric        |
//| positive-definite eigenproblem:                                  |
//|     A*x = lambda*B*x (1) or                                      |
//|     A*B*x = lambda*x (2) or                                      |
//|     B*A*x = lambda*x (3).                                        |
//| where A is a symmetric matrix, B - symmetric positive-definite   |
//| matrix. The problem is solved by reducing it to an ordinary      |
//| symmetric eigenvalue problem.                                    |
//| Input parameters:                                                |
//|     A           -   symmetric matrix which is given by its upper |
//|                     or lower triangular part.                    |
//|                     Array whose indexes range within             |
//|                     [0..N-1, 0..N-1].                            |
//|     N           -   size of matrices A and B.                    |
//|     IsUpperA    -   storage format of matrix A.                  |
//|     B           -   symmetric positive-definite matrix which is  |
//|                     given by its upper or lower triangular part. |
//|                     Array whose indexes range within             |
//|                     [0..N-1, 0..N-1].                            |
//|     IsUpperB    -   storage format of matrix B.                  |
//|     ZNeeded     -   if ZNeeded is equal to:                      |
//|                      * 0, the eigenvectors are not returned;     |
//|                      * 1, the eigenvectors are returned.         |
//|     ProblemType -   if ProblemType is equal to:                  |
//|                      * 1, the following problem is solved:       |
//|                           A*x = lambda*B*x;                      |
//|                      * 2, the following problem is solved:       |
//|                           A*B*x = lambda*x;                      |
//|                      * 3, the following problem is solved:       |
//|                           B*A*x = lambda*x.                      |
//| Output parameters:                                               |
//|     D           -   eigenvalues in ascending order.              |
//|                     Array whose index ranges within [0..N-1].    |
//|     Z           -   if ZNeeded is equal to:                      |
//|                      * 0, Z hasn?t changed;                      |
//|                      * 1, Z contains eigenvectors.               |
//|                     Array whose indexes range within             |
//|                     [0..N-1, 0..N-1].                            |
//|                     The eigenvectors are stored in matrix        |
//|                     columns. It should be noted that the         |
//|                     eigenvectors in such problems do not form an |
//|                     orthogonal system.                           |
//| Result:                                                          |
//|     True, if the problem was solved successfully.                |
//|     False, if the error occurred during the Cholesky             |
//|     decomposition of matrix B (the matrix isn?t                  |
//|     positive-definite) or during the work of the iterative       |
//|     algorithm for solving the symmetric eigenproblem.            |
//| See also the GeneralizedSymmetricDefiniteEVDReduce subroutine.   |
//+------------------------------------------------------------------+
static bool CSpdGEVD::SMatrixGEVD(CMatrixDouble &ca,const int n,const bool isuppera,
                                  CMatrixDouble &b,const bool isupperb,const int zneeded,
                                  const int problemtype,double &d[],CMatrixDouble &z)
  {
//--- create variables
   bool   result;
   bool   isupperr;
   int    j1=0;
   int    j2=0;
   int    j1inc=0;
   int    j2inc=0;
   int    i=0;
   int    j=0;
   double v=0;
   int    i_=0;
//--- create matrix
   CMatrixDouble r;
   CMatrixDouble t;
//--- create copy
   CMatrixDouble a;
   a=ca;
//--- Reduce and solve
   result=SMatrixGEVDReduce(a,n,isuppera,b,isupperb,problemtype,r,isupperr);
//--- check
   if(!result)
      return(result);
//--- get result
   result=CEigenVDetect::SMatrixEVD(a,n,zneeded,isuppera,d,t);
//--- check
   if(!result)
      return(result);
//--- Transform eigenvectors if needed
   if(zneeded!=0)
     {
      //--- fill Z with zeros
      z.Resize(n,n);
      for(j=0;j<n;j++)
         z[0].Set(j,0.0);
      for(i=1;i<n;i++)
        {
         for(i_=0;i_<n;i_++)
            z[i].Set(i_,z[0][i_]);
        }
      //--- Setup R properties
      if(isupperr)
        {
         j1=0;
         j2=n-1;
         j1inc=1;
         j2inc=0;
        }
      else
        {
         j1=0;
         j2=0;
         j1inc=0;
         j2inc=1;
        }
      //--- Calculate R*Z
      for(i=0;i<n;i++)
        {
         for(j=j1;j<=j2;j++)
           {
            v=r[i][j];
            for(i_=0;i_<n;i_++)
               z[i].Set(i_,z[i][i_]+v*t[j][i_]);
           }
         j1=j1+j1inc;
         j2=j2+j2inc;
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Algorithm for reduction of the following generalized symmetric   |
//| positive- definite eigenvalue problem:                           |
//|     A*x = lambda*B*x (1) or                                      |
//|     A*B*x = lambda*x (2) or                                      |
//|     B*A*x = lambda*x (3)                                         |
//| to the symmetric eigenvalues problem C*y = lambda*y (eigenvalues |
//| of this and the given problems are the same, and the eigenvectors|
//| of the given problem could be obtained by multiplying the        |
//| obtained eigenvectors by the transformation matrix x = R*y).     |
//| Here A is a symmetric matrix, B - symmetric positive-definite    |
//| matrix.                                                          |
//| Input parameters:                                                |
//|     A           -   symmetric matrix which is given by its upper |
//|                     or lower triangular part.                    |
//|                     Array whose indexes range within             |
//|                     [0..N-1, 0..N-1].                            |
//|     N           -   size of matrices A and B.                    |
//|     IsUpperA    -   storage format of matrix A.                  |
//|     B           -   symmetric positive-definite matrix which is  |
//|                     given by its upper or lower triangular part. |
//|                     Array whose indexes range within             |
//|                     [0..N-1, 0..N-1].                            |
//|     IsUpperB    -   storage format of matrix B.                  |
//|     ProblemType -   if ProblemType is equal to:                  |
//|                      * 1, the following problem is solved:       |
//|                           A*x = lambda*B*x;                      |
//|                      * 2, the following problem is solved:       |
//|                           A*B*x = lambda*x;                      |
//|                      * 3, the following problem is solved:       |
//|                           B*A*x = lambda*x.                      |
//| Output parameters:                                               |
//|     A           -   symmetric matrix which is given by its upper |
//|                     or lower triangle depending on IsUpperA.     |
//|                     Contains matrix C. Array whose indexes range |
//|                     within [0..N-1, 0..N-1].                     |
//|     R           -   upper triangular or low triangular           |
//|                     transformation matrix which is used to obtain|
//|                     the eigenvectors of a given problem as the   |
//|                     product of eigenvectors of C (from the right)|
//|                     and matrix R (from the left). If the matrix  |
//|                     is upper triangular, the elements below the  |
//|                     main diagonal are equal to 0 (and vice versa)|
//|                     Thus, we can perform the multiplication      |
//|                     without taking into account the internal     |
//|                     structure (which is an easier though less    |
//|                     effective way). Array whose indexes range    |
//|                     within [0..N-1, 0..N-1].                     |
//|     IsUpperR    -   type of matrix R (upper or lower triangular).|
//| Result:                                                          |
//|     True, if the problem was reduced successfully.               |
//|     False, if the error occurred during the Cholesky             |
//|         decomposition of matrix B (the matrix is not             |
//|         positive-definite).                                      |
//+------------------------------------------------------------------+
static bool CSpdGEVD::SMatrixGEVDReduce(CMatrixDouble &a,const int n,const bool isuppera,
                                        CMatrixDouble &b,const bool isupperb,const int problemtype,
                                        CMatrixDouble &r,bool &isupperr)
  {
//--- create variables
   bool   result;
   int    i=0;
   int    j=0;
   double v=0;
   int    info=0;
   int    i_=0;
   int    i1_=0;
//--- create matrix
   CMatrixDouble t;
//--- create arrays
   double w1[];
   double w2[];
   double w3[];
//--- object of class
   CMatInvReport rep;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": N<=0!"))
      return(false);
//--- check
   if(!CAp::Assert(problemtype==1 || problemtype==2 || problemtype==3,__FUNCTION__+": incorrect ProblemType!"))
      return(false);
//--- initialization
   result=true;
//--- Problem 1:  A*x = lambda*B*x
//--- Reducing to:
//---     C*y = lambda*y
//---     C = L^(-1) * A * L^(-T)
//---     x = L^(-T) * y
   if(problemtype==1)
     {
      //--- Factorize B in T: B = LL'
      t.Resize(n,n);
      //--- check
      if(isupperb)
        {
         for(i=0;i<n;i++)
            for(i_=i;i_<n;i_++)
               t[i_].Set(i,b[i][i_]);
        }
      else
        {
         for(i=0;i<n;i++)
            for(i_=0;i_<=i;i_++)
               t[i].Set(i_,b[i][i_]);
        }
      //--- check
      if(!CTrFac::SPDMatrixCholesky(t,n,false))
         return(false);
      //--- Invert L in T
      CMatInv::RMatrixTrInverse(t,n,false,false,info,rep);
      //--- check
      if(info<=0)
         return(false);
      //--- Build L^(-1) * A * L^(-T) in R
      ArrayResizeAL(w1,n+1);
      ArrayResizeAL(w2,n+1);
      r.Resize(n,n);
      for(j=1;j<=n;j++)
        {
         //--- Form w2 = A * l'(j) (here l'(j) is j-th column of L^(-T))
         i1_=-1;
         for(i_=1;i_<=j;i_++)
            w1[i_]=t[j-1][i_+i1_];
         //--- function call
         CSblas::SymmetricMatrixVectorMultiply(a,isuppera,0,j-1,w1,1.0,w2);
         //--- check
         if(isuppera)
            CBlas::MatrixVectorMultiply(a,0,j-1,j,n-1,true,w1,1,j,1.0,w2,j+1,n,0.0);
         else
            CBlas::MatrixVectorMultiply(a,j,n-1,0,j-1,false,w1,1,j,1.0,w2,j+1,n,0.0);
         //--- Form l(i)*w2 (here l(i) is i-th row of L^(-1))
         for(i=1;i<=n;i++)
           {
            i1_=1;
            v=0.0;
            for(i_=0;i_<i;i_++)
               v+=t[i-1][i_]*w2[i_+i1_];
            r[i-1].Set(j-1,v);
           }
        }
      // Copy R to A
      for(i=0;i<n;i++)
        {
         for(i_=0;i_<n;i_++)
            a[i].Set(i_,r[i][i_]);
        }
      //--- Copy L^(-1) from T to R and transpose
      isupperr=true;
      for(i=0;i<n;i++)
        {
         for(j=0;j<i;j++)
            r[i].Set(j,0);
        }
      for(i=0;i<n;i++)
        {
         for(i_=i;i_<n;i_++)
            r[i].Set(i_,t[i_][i]);
        }
      //--- return result
      return(result);
     }
//--- Problem 2:  A*B*x = lambda*x
//--- or
//--- problem 3:  B*A*x = lambda*x
//--- Reducing to:
//---     C*y = lambda*y
//---     C = U * A * U'
//---     B = U'* U
   if(problemtype==2||problemtype==3)
     {
      //--- Factorize B in T: B = U'*U
      t.Resize(n,n);
      //--- check
      if(isupperb)
        {
         for(i=0;i<n;i++)
            for(i_=i;i_<n;i_++)
               t[i].Set(i_,b[i][i_]);
        }
      else
        {
         for(i=0;i<n;i++)
            for(i_=i;i_<n;i_++)
               t[i].Set(i_,b[i_][i]);
        }
      //--- check
      if(!CTrFac::SPDMatrixCholesky(t,n,true))
         return(false);
      //--- Build U * A * U' in R
      ArrayResizeAL(w1,n+1);
      ArrayResizeAL(w2,n+1);
      ArrayResizeAL(w3,n+1);
      r.Resize(n,n);
      for(j=1;j<=n;j++)
        {
         //--- Form w2 = A * u'(j) (here u'(j) is j-th column of U')
         i1_=(j-1)-(1);
         for(i_=1;i_<=n-j+1;i_++)
            w1[i_]=t[j-1][i_+i1_];
         //--- function call
         CSblas::SymmetricMatrixVectorMultiply(a,isuppera,j-1,n-1,w1,1.0,w3);
         i1_=(1)-(j);
         for(i_=j;i_<=n;i_++)
            w2[i_]=w3[i_+i1_];
         i1_=(j-1)-(j);
         for(i_=j;i_<=n;i_++)
            w1[i_]=t[j-1][i_+i1_];
         //--- check
         if(isuppera)
            CBlas::MatrixVectorMultiply(a,0,j-2,j-1,n-1,false,w1,j,n,1.0,w2,1,j-1,0.0);
         else
            CBlas::MatrixVectorMultiply(a,j-1,n-1,0,j-2,true,w1,j,n,1.0,w2,1,j-1,0.0);
         //--- Form u(i)*w2 (here u(i) is i-th row of U)
         for(i=1;i<=n;i++)
           {
            i1_=(i)-(i-1);
            v=0.0;
            for(i_=i-1;i_<n;i_++)
               v+=t[i-1][i_]*w2[i_+i1_];
            r[i-1].Set(j-1,v);
           }
        }
      //--- Copy R to A
      for(i=0;i<n;i++)
        {
         for(i_=0;i_<n;i_++)
            a[i].Set(i_,r[i][i_]);
        }
      //--- check
      if(problemtype==2)
        {
         //--- Invert U in T
         CMatInv::RMatrixTrInverse(t,n,true,false,info,rep);
         //--- check
         if(info<=0)
            return(false);
         //--- Copy U^-1 from T to R
         isupperr=true;
         for(i=0;i<n;i++)
           {
            for(j=0;j<i;j++)
               r[i].Set(j,0);
           }
         for(i=0;i<n;i++)
           {
            for(i_=i;i_<n;i_++)
               r[i].Set(i_,t[i][i_]);
           }
        }
      else
        {
         //--- Copy U from T to R and transpose
         isupperr=false;
         for(i=0;i<n;i++)
           {
            for(j=i+1;j<n;j++)
               r[i].Set(j,0);
           }
         for(i=0;i<n;i++)
           {
            for(i_=i;i_<n;i_++)
               r[i_].Set(i,t[i][i_]);
           }
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Inverse and update                                               |
//+------------------------------------------------------------------+
class CInverseUpdate
  {
public:
                     CInverseUpdate(void);
                    ~CInverseUpdate(void);
   //--- methods
   static void       RMatrixInvUpdateSimple(CMatrixDouble &inva,const int n,const int updrow,const int updcolumn,const double updval);
   static void       RMatrixInvUpdateRow(CMatrixDouble &inva,const int n,const int updrow,double &v[]);
   static void       RMatrixInvUpdateColumn(CMatrixDouble &inva,const int n,const int updcolumn,double &u[]);
   static void       RMatrixInvUpdateUV(CMatrixDouble &inva,const int n,double &u[],double &v[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CInverseUpdate::CInverseUpdate(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CInverseUpdate::~CInverseUpdate(void)
  {

  }
//+------------------------------------------------------------------+
//| Inverse matrix update by the Sherman-Morrison formula            |
//| The algorithm updates matrix A^-1 when adding a number to an     |
//| element of matrix A.                                             |
//| Input parameters:                                                |
//|     InvA    -   inverse of matrix A.                             |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..N-1].                                |
//|     N       -   size of matrix A.                                |
//|     UpdRow  -   row where the element to be updated is stored.   |
//|     UpdColumn - column where the element to be updated is stored.|
//|     UpdVal  -   a number to be added to the element.             |
//| Output parameters:                                               |
//|     InvA    -   inverse of modified matrix A.                    |
//+------------------------------------------------------------------+
static void CInverseUpdate::RMatrixInvUpdateSimple(CMatrixDouble &inva,const int n,
                                                   const int updrow,const int updcolumn,
                                                   const double updval)
  {
//--- create arrays
   double t1[];
   double t2[];
//--- create variables
   int    i=0;
   double lambdav=0;
   double vt=0;
   int    i_=0;
//--- check
   if(!CAp::Assert(updrow>=0 && updrow<n,"RMatrixInvUpdateSimple: incorrect UpdRow!"))
      return;
//--- check
   if(!CAp::Assert(updcolumn>=0 && updcolumn<n,"RMatrixInvUpdateSimple: incorrect UpdColumn!"))
      return;
//--- allocation
   ArrayResizeAL(t1,n);
   ArrayResizeAL(t2,n);
//--- T1=InvA * U
   for(i_=0;i_<n;i_++)
      t1[i_]=inva[i_][updrow];
//--- T2=v*InvA
   for(i_=0;i_<n;i_++)
      t2[i_]=inva[updcolumn][i_];
//--- Lambda=v * InvA * U
   lambdav=updval*inva[updcolumn][updrow];
//--- InvA=InvA - correction
   for(i=0;i<n;i++)
     {
      vt=updval*t1[i];
      vt=vt/(1+lambdav);
      for(i_=0;i_<n;i_++)
         inva[i].Set(i_,inva[i][i_]-vt*t2[i_]);
     }
  }
//+------------------------------------------------------------------+
//| Inverse matrix update by the Sherman-Morrison formula            |
//| The algorithm updates matrix A^-1 when adding a vector to a row  |
//| of matrix A.                                                     |
//| Input parameters:                                                |
//|     InvA    -   inverse of matrix A.                             |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..N-1].                                |
//|     N       -   size of matrix A.                                |
//|     UpdRow  -   the row of A whose vector V was added.           |
//|                 0 <= Row <= N-1                                  |
//|     V       -   the vector to be added to a row.                 |
//|                 Array whose index ranges within [0..N-1].        |
//| Output parameters:                                               |
//|     InvA    -   inverse of modified matrix A.                    |
//+------------------------------------------------------------------+
static void CInverseUpdate::RMatrixInvUpdateRow(CMatrixDouble &inva,const int n,
                                                const int updrow,double &v[])
  {
//--- create arrays
   double t1[];
   double t2[];
//--- create variables
   int    i=0;
   int    j=0;
   double lambdav=0;
   double vt=0;
   int    i_=0;
//--- allocation
   ArrayResizeAL(t1,n);
   ArrayResizeAL(t2,n);
//--- T1=InvA * U
   for(i_=0;i_<n;i_++)
      t1[i_]=inva[i_][updrow];
//--- T2=v*InvA
//--- Lambda=v * InvA * U
   for(j=0;j<=n-1;j++)
     {
      vt=0.0;
      for(i_=0;i_<n;i_++)
         vt+=v[i_]*inva[i_][j];
      t2[j]=vt;
     }
   lambdav=t2[updrow];
//--- InvA=InvA - correction
   for(i=0;i<n;i++)
     {
      vt=t1[i]/(1+lambdav);
      for(i_=0;i_<n;i_++)
         inva[i].Set(i_,inva[i][i_]-vt*t2[i_]);
     }
  }
//+------------------------------------------------------------------+
//| Inverse matrix update by the Sherman-Morrison formula            |
//| The algorithm updates matrix A^-1 when adding a vector to a      |
//| column of matrix A.                                              |
//| Input parameters:                                                |
//|     InvA        -   inverse of matrix A.                         |
//|                     Array whose indexes range within             |
//|                     [0..N-1, 0..N-1].                            |
//|     N           -   size of matrix A.                            |
//|     UpdColumn   -   the column of A whose vector U was added.    |
//|                     0 <= UpdColumn <= N-1                        |
//|     U           -   the vector to be added to a column.          |
//|                     Array whose index ranges within [0..N-1].    |
//| Output parameters:                                               |
//|     InvA        -   inverse of modified matrix A.                |
//+------------------------------------------------------------------+
static void CInverseUpdate::RMatrixInvUpdateColumn(CMatrixDouble &inva,const int n,
                                                   const int updcolumn,double &u[])
  {
//--- create arrays
   double t1[];
   double t2[];
//--- create variables
   int    i=0;
   double lambdav=0;
   double vt=0;
   int    i_=0;
//--- allocation
   ArrayResizeAL(t1,n);
   ArrayResizeAL(t2,n);
//--- T1=InvA * U
//--- Lambda=v * InvA * U
   for(i=0;i<n;i++)
     {
      vt=0.0;
      for(i_=0;i_<n;i_++)
         vt+=inva[i][i_]*u[i_];
      t1[i]=vt;
     }
   lambdav=t1[updcolumn];
//--- T2=v*InvA
   for(i_=0;i_<n;i_++)
      t2[i_]=inva[updcolumn][i_];
//--- InvA=InvA - correction
   for(i=0;i<n;i++)
     {
      vt=t1[i]/(1+lambdav);
      for(i_=0;i_<n;i_++)
         inva[i].Set(i_,inva[i][i_]-vt*t2[i_]);
     }
  }
//+------------------------------------------------------------------+
//| Inverse matrix update by the Sherman-Morrison formula            |
//| The algorithm computes the inverse of matrix A+u*v? by using the |
//| given matrix A^-1 and the vectors u and v.                       |
//| Input parameters:                                                |
//|     InvA    -   inverse of matrix A.                             |
//|                 Array whose indexes range within                 |
//|                 [0..N-1, 0..N-1].                                |
//|     N       -   size of matrix A.                                |
//|     U       -   the vector modifying the matrix.                 |
//|                 Array whose index ranges within [0..N-1].        |
//|     V       -   the vector modifying the matrix.                 |
//|                 Array whose index ranges within [0..N-1].        |
//| Output parameters:                                               |
//|     InvA - inverse of matrix A + u*v'.                           |
//+------------------------------------------------------------------+
static void CInverseUpdate::RMatrixInvUpdateUV(CMatrixDouble &inva,const int n,
                                               double &u[],double &v[])
  {
//--- create arrays
   double t1[];
   double t2[];
//--- create variables
   int i=0;
   int j=0;
   double lambdav=0;
   double vt=0;
   int i_=0;
//--- allocation
   ArrayResizeAL(t1,n);
   ArrayResizeAL(t2,n);
//--- T1=InvA * U
//--- Lambda=v * T1
   for(i=0;i<n;i++)
     {
      vt=0.0;
      for(i_=0;i_<n;i_++)
         vt+=inva[i][i_]*u[i_];
      t1[i]=vt;
     }
   lambdav=0.0;
   for(i_=0;i_<n;i_++)
      lambdav+=v[i_]*t1[i_];
//--- T2=v*InvA
   for(j=0;j<=n-1;j++)
     {
      vt=0.0;
      for(i_=0;i_<n;i_++)
         vt+=v[i_]*inva[i_][j];
      t2[j]=vt;
     }
//--- InvA=InvA - correction
   for(i=0;i<n;i++)
     {
      vt=t1[i]/(1+lambdav);
      for(i_=0;i_<n;i_++)
         inva[i].Set(i_,inva[i][i_]-vt*t2[i_]);
     }
  }
//+------------------------------------------------------------------+
//| Schur decomposition                                              |
//+------------------------------------------------------------------+
class CSchur
  {
public:
   //--- constructor, destructor
                     CSchur(void);
                    ~CSchur(void);
   //--- method
   static bool       RMatrixSchur(CMatrixDouble &a,const int n,CMatrixDouble &s);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSchur::CSchur(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSchur::~CSchur(void)
  {

  }
//+------------------------------------------------------------------+
//| Subroutine performing the Schur decomposition of a general matrix|
//| by using the QR algorithm with multiple shifts.                  |
//| The source matrix A is represented as S'*A*S = T, where S is an  |
//| orthogonal matrix (Schur vectors), T - upper quasi-triangular    |
//| matrix (with blocks of sizes 1x1 and 2x2 on the main diagonal).  |
//| Input parameters:                                                |
//|     A   -   matrix to be decomposed.                             |
//|             Array whose indexes range within [0..N-1, 0..N-1].   |
//|     N   -   size of A, N>=0.                                     |
//| Output parameters:                                               |
//|     A   -   contains matrix T.                                   |
//|             Array whose indexes range within [0..N-1, 0..N-1].   |
//|     S   -   contains Schur vectors.                              |
//|             Array whose indexes range within [0..N-1, 0..N-1].   |
//| Note 1:                                                          |
//|     The block structure of matrix T can be easily recognized:    |
//|     since all the elements below the blocks are zeros, the       |
//|     elements a[i+1,i] which are equal to 0 show the block border.|
//| Note 2:                                                          |
//|     The algorithm performance depends on the value of the        |
//|     internal parameter NS of the InternalSchurDecomposition      |
//|     subroutine which defines the number of shifts in the QR      |
//|     algorithm (similarly to the block width in block-matrix      |
//|     algorithms in linear algebra). If you require maximum        |
//|     performance on your machine, it is recommended to adjust     |
//|     this parameter manually.                                     |
//| Result:                                                          |
//|     True,                                                        |
//|         if the algorithm has converged and parameters A and S    |
//|         contain the result.                                      |
//|     False,                                                       |
//|         if the algorithm has not converged.                      |
//| Algorithm implemented on the basis of the DHSEQR subroutine      |
//| (LAPACK 3.0 library).                                            |
//+------------------------------------------------------------------+
static bool CSchur::RMatrixSchur(CMatrixDouble &a,const int n,CMatrixDouble &s)
  {
//--- create variables
   bool result;
   int  info=0;
   int  i=0;
   int  j=0;
//--- create arrays
   double tau[];
   double wi[];
   double wr[];
//--- create matrix
   CMatrixDouble a1;
   CMatrixDouble s1;
//--- Upper Hessenberg form of the 0-based matrix
   COrtFac::RMatrixHessenberg(a,n,tau);
   COrtFac::RMatrixHessenbergUnpackQ(a,n,tau,s);
//--- Convert from 0-based arrays to 1-based,
//--- then call InternalSchurDecomposition
//--- Awkward, of course, but Schur decompisiton subroutine
//--- is too complex to fix it.
   a1.Resize(n+1,n+1);
   s1.Resize(n+1,n+1);
   for(i=1;i<=n;i++)
     {
      for(j=1;j<=n;j++)
        {
         a1[i].Set(j,a[i-1][j-1]);
         s1[i].Set(j,s[i-1][j-1]);
        }
     }
//--- function call
   CHsSchur::InternalSchurDecomposition(a1,n,1,1,wr,wi,s1,info);
   result=info==0;
//--- convert from 1-based arrays to -based
   for(i=1;i<=n;i++)
     {
      for(j=1;j<=n;j++)
        {
         a[i-1].Set(j-1,a1[i][j]);
         s[i-1].Set(j-1,s1[i][j]);
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
