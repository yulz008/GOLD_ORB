//+------------------------------------------------------------------+
//|                                               fasttransforms.mqh |
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
#include "complex.mqh"
#include "alglibinternal.mqh"
//+------------------------------------------------------------------+
//| Fast Fourier Transform                                           |
//+------------------------------------------------------------------+
class CFastFourierTransform
  {
public:
                     CFastFourierTransform(void);
                    ~CFastFourierTransform(void);

   static void       FFTC1D(al_complex &a[],const int n);
   static void       FFTC1DInv(al_complex &a[],const int n);
   static void       FFTR1D(double &a[],const int n,al_complex &f[]);
   static void       FFTR1DInv(al_complex &f[],const int n,double &a[]);
   static void       FFTR1DInternalEven(double &a[],const int n,double &buf[],CFtPlan &plan);
   static void       FFTR1DInvInternalEven(double &a[],const int n,double &buf[],CFtPlan &plan);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CFastFourierTransform::CFastFourierTransform(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CFastFourierTransform::~CFastFourierTransform(void)
  {

  }
//+------------------------------------------------------------------+
//| 1-dimensional complex FFT.                                       |
//| Array size N may be arbitrary number (composite or prime).       |
//| Composite N's are handled with cache-oblivious variation of a    |
//| Cooley-Tukey algorithm. Small prime-factors are transformed using|
//| hard coded codelets (similar to FFTW codelets, but without       |
//| low-level optimization), large prime-factors are handled with    |
//| Bluestein's algorithm.                                           |
//| Fastests transforms are for smooth N's (prime factors are 2, 3,  |
//| 5 only), most fast for powers of 2. When N have prime factors    |
//| larger than these, but orders of magnitude smaller than N,       |
//| computations will be about 4 times slower than for nearby highly |
//| composite N's. When N itself is prime, speed will be 6 times     |
//| lower.                                                           |
//| Algorithm has O(N*logN) complexity for any N (composite or       |
//| prime).                                                          |
//| INPUT PARAMETERS                                                 |
//|     A   -   array[0..N-1] - complex function to be transformed   |
//|     N   -   problem size                                         |
//| OUTPUT PARAMETERS                                                |
//|     A   -   DFT of a input array, array[0..N-1]                  |
//|             A_out[j] = SUM(A_in[k]*exp(-2*pi*sqrt(-1)*j*k/N),    |
//|             k = 0..N-1)                                          |
//+------------------------------------------------------------------+
static void CFastFourierTransform::FFTC1D(al_complex &a[],const int n)
  {
//--- create a variable
   int i=0;
//--- create array
   double buf[];
//--- object of class
   CFtPlan plan;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(a)>=n,__FUNCTION__+": Length(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteComplexVector(a,n),__FUNCTION__+": A contains infinite or NAN values!"))
      return;
//--- Special case: N=1,FFT is just identity transform.
//--- After this block we assume that N is strictly greater than 1.
   if(n==1)
      return;
//--- convert input array to the more convinient format
   ArrayResizeAL(buf,2*n);
   for(i=0;i<=n-1;i++)
     {
      buf[2*i+0]=a[i].re;
      buf[2*i+1]=a[i].im;
     }
//--- Generate plan and execute it.
//--- Plan is a combination of a successive factorizations of N and
//--- precomputed data. It is much like a FFTW plan,but is not stored
//--- between subroutine calls and is much simpler.
   CFtBase::FtBaseGenerateComplexFFtPlan(n,plan);
   CFtBase::FtBaseExecutePlan(buf,0,n,plan);
//--- result
   for(i=0;i<=n-1;i++)
     {
      a[i].re=buf[2*i+0];
      a[i].im=buf[2*i+1];
     }
  }
//+------------------------------------------------------------------+
//| 1-dimensional complex inverse FFT.                               |
//| Array size N may be arbitrary number (composite or prime).       |
//| Algorithm has O(N*logN) complexity for any N (composite or prime)|
//| See FFTC1D() description for more information about algorithm    |
//| performance.                                                     |
//| INPUT PARAMETERS                                                 |
//|     A   -   array[0..N-1] - complex array to be transformed      |
//|     N   -   problem size                                         |
//| OUTPUT PARAMETERS                                                |
//|     A   -   inverse DFT of a input array, array[0..N-1]          |
//|             A_out[j] = SUM(A_in[k]/N*exp(+2*pi*sqrt(-1)*j*k/N),  |
//|             k = 0..N-1)                                          |
//+------------------------------------------------------------------+
static void CFastFourierTransform::FFTC1DInv(al_complex &a[],const int n)
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(a)>=n,__FUNCTION__+": Length(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteComplexVector(a,n),__FUNCTION__+": A contains infinite or NAN values!"))
      return;
//--- Inverse DFT can be expressed in terms of the DFT as
//---     invfft(x)=fft(x')'/N
//--- here x' means conj(x).
   for(i=0;i<=n-1;i++)
      a[i].im=-a[i].im;
//--- function call
   FFTC1D(a,n);
//--- change values
   for(i=0;i<=n-1;i++)
     {
      a[i].re=a[i].re/n;
      a[i].im=-(a[i].im/n);
     }
  }
//+------------------------------------------------------------------+
//| 1-dimensional real FFT.                                          |
//| Algorithm has O(N*logN) complexity for any N (composite or       |
//| prime).                                                          |
//| INPUT PARAMETERS                                                 |
//|     A   -   array[0..N-1] - real function to be transformed      |
//|     N   -   problem size                                         |
//| OUTPUT PARAMETERS                                                |
//|     F   -   DFT of a input array, array[0..N-1]                  |
//|             F[j] = SUM(A[k]*exp(-2*pi*sqrt(-1)*j*k/N),           |
//|             k = 0..N-1)                                          |
//| NOTE:                                                            |
//|     F[] satisfies symmetry property F[k] = conj(F[N-k]), so just |
//| one half of array is usually needed. But for convinience         |
//| subroutine returns full complex array (with frequencies above    |
//| N/2), so its result may be used by other FFT-related subroutines.|
//+------------------------------------------------------------------+
static void CFastFourierTransform::FFTR1D(double &a[],const int n,al_complex &f[])
  {
//--- create variables
   int        i=0;
   int        n2=0;
   int        idx=0;
   al_complex hn=0;
   al_complex hmnc=0;
   al_complex v=0;
   int        i_=0;
//--- create array
   double buf[];
//--- object of class
   CFtPlan plan;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(a)>=n,__FUNCTION__+": Length(A)<N!"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteVector(a,n),__FUNCTION__+": A contains infinite or NAN values!"))
      return;
//--- Special cases:
//--- * N=1,FFT is just identity transform.
//--- * N=2,FFT is simple too
//--- After this block we assume that N is strictly greater than 2
   if(n==1)
     {
      //--- allocation
      ArrayResizeAL(f,1);
      f[0]=a[0];
      //--- exit the function
      return;
     }
//--- check
   if(n==2)
     {
      //--- allocation
      ArrayResizeAL(f,2);
      f[0].re=a[0]+a[1];
      f[0].im=0;
      f[1].re=a[0]-a[1];
      f[1].im=0;
      //--- exit the function
      return;
     }
//--- Choose between odd-size and even-size FFTs
   if(n%2==0)
     {
      //--- even-size real FFT,use reduction to the complex task
      n2=n/2;
      //--- allocation
      ArrayResizeAL(buf,n);
      for(i_=0;i_<=n-1;i_++)
         buf[i_]=a[i_];
      //--- function call
      CFtBase::FtBaseGenerateComplexFFtPlan(n2,plan);
      CFtBase::FtBaseExecutePlan(buf,0,n2,plan);
      //--- allocation
      ArrayResizeAL(f,n);
      //--- calculation
      for(i=0;i<=n2;i++)
        {
         idx=2*(i%n2);
         hn.re=buf[idx+0];
         hn.im=buf[idx+1];
         idx=2*((n2-i)%n2);
         hmnc.re=buf[idx+0];
         hmnc.im=-buf[idx+1];
         v.re=-MathSin(-(2*M_PI*i/n));
         v.im=MathCos(-(2*M_PI*i/n));
         f[i]=hn+hmnc-v*(hn-hmnc);
         f[i].re=0.5*f[i].re;
         f[i].im=0.5*f[i].im;
        }
      //--- copy
      for(i=n2+1;i<=n-1;i++)
         f[i]=CMath::Conj(f[n-i]);
     }
   else
     {
      //--- use complex FFT
      ArrayResizeAL(f,n);
      //--- copy
      for(i=0;i<=n-1;i++)
         f[i]=a[i];
      //--- function call
      FFTC1D(f,n);
     }
  }
//+------------------------------------------------------------------+
//| 1-dimensional real inverse FFT.                                  |
//| Algorithm has O(N*logN) complexity for any N (composite or       |
//| prime).                                                          |
//| INPUT PARAMETERS                                                 |
//|     F   -   array[0..floor(N/2)] - frequencies from forward real |
//|             FFT                                                  |
//|     N   -   problem size                                         |
//| OUTPUT PARAMETERS                                                |
//|     A   -   inverse DFT of a input array, array[0..N-1]          |
//| NOTE:                                                            |
//|     F[] should satisfy symmetry property F[k] = conj(F[N-k]),    |
//| so just one half of frequencies array is needed - elements from 0|
//| to floor(N/2). F[0] is ALWAYS real. If N is even F[floor(N/2)] is|
//| real too. If N is odd, then F[floor(N/2)] has no special         |
//| properties.                                                      |
//| Relying on properties noted above, FFTR1DInv subroutine uses only|
//| elements from 0th to floor(N/2)-th. It ignores imaginary part of |
//| F[0], and in case N is even it ignores imaginary part of         |
//| F[floor(N/2)] too.                                               |
//| When you call this function using full arguments list -          |
//| "FFTR1DInv(F,N,A)"                                               |
//| - you can pass either either frequencies array with N elements or|
//| reduced array with roughly N/2 elements - subroutine will        |
//| successfully transform both.                                     |
//| If you call this function using reduced arguments list -         |
//| "FFTR1DInv(F,A)" - you must pass FULL array with N elements      |
//| (although higher N/2 are still not used) because array size is   |
//| used to automatically determine FFT length                       |
//+------------------------------------------------------------------+
static void CFastFourierTransform::FFTR1DInv(al_complex &f[],const int n,
                                             double &a[])
  {
//--- create a variable
   int i=0;
//--- create arrays
   double     h[];
   al_complex fh[];
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(f)>=(int)MathFloor((double)n/2.0)+1,__FUNCTION__+": Length(F)<Floor(N/2)+1!"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(f[0].re),__FUNCTION__+": F contains infinite or NAN values!"))
      return;
   for(i=1;i<=(int)MathFloor((double)n/2.0)-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(f[i].re) && CMath::IsFinite(f[i].im),__FUNCTION__+": F contains infinite or NAN values!"))
         return;
     }
//--- check
   if(!CAp::Assert(CMath::IsFinite(f[(int)MathFloor((double)n/2.0)].re),__FUNCTION__+": F contains infinite or NAN values!"))
      return;
//--- check
   if(n%2!=0)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(f[(int)MathFloor((double)n/2.0)].im),__FUNCTION__+": F contains infinite or NAN values!"))
         return;
     }
//--- Special case: N=1,FFT is just identity transform.
//--- After this block we assume that N is strictly greater than 1.
   if(n==1)
     {
      //--- allocation
      ArrayResizeAL(a,1);
      a[0]=f[0].re;
      //--- exit the function
      return;
     }
//--- inverse real FFT is reduced to the inverse real FHT,
//--- which is reduced to the forward real FHT,
//--- which is reduced to the forward real FFT.
//--- Don't worry,it is really compact and efficient reduction :)
   ArrayResizeAL(h,n);
   ArrayResizeAL(a,n);
//--- change values
   h[0]=f[0].re;
   for(i=1;i<=(int)MathFloor((double)n/2.0)-1;i++)
     {
      h[i]=f[i].re-f[i].im;
      h[n-i]=f[i].re+f[i].im;
     }
//--- check
   if(n%2==0)
      h[(int)MathFloor((double)n/2.0)]=f[(int)MathFloor((double)n/2.0)].re;
   else
     {
      h[(int)MathFloor((double)n/2.0)]=f[(int)MathFloor((double)n/2.0)].re-f[(int)MathFloor((double)n/2.0)].im;
      h[(int)MathFloor((double)n/2.0)+1]=f[(int)MathFloor((double)n/2.0)].re+f[(int)MathFloor((double)n/2.0)].im;
     }
//--- function call
   FFTR1D(h,n,fh);
//--- change values
   for(i=0;i<=n-1;i++)
      a[i]=(fh[i].re-fh[i].im)/n;
  }
//+------------------------------------------------------------------+
//| Internal subroutine. Never call it directly!                     |
//+------------------------------------------------------------------+
static void CFastFourierTransform::FFTR1DInternalEven(double &a[],const int n,
                                                      double &buf[],CFtPlan &plan)
  {
//--- create variables
   double     x=0;
   double     y=0;
   int        i=0;
   int        n2=0;
   int        idx=0;
   al_complex hn=0;
   al_complex hmnc=0;
   al_complex v=0;
   int        i_=0;
//--- check
   if(!CAp::Assert(n>0 && n%2==0,__FUNCTION__+": incorrect N!"))
      return;
//--- Special cases:
//--- * N=2
//--- After this block we assume that N is strictly greater than 2
   if(n==2)
     {
      //--- change values
      x=a[0]+a[1];
      y=a[0]-a[1];
      a[0]=x;
      a[1]=y;
      //--- exit the function
      return;
     }
//--- even-size real FFT,use reduction to the complex task
   n2=n/2;
//--- copy
   for(i_=0;i_<=n-1;i_++)
      buf[i_]=a[i_];
//--- function call
   CFtBase::FtBaseExecutePlan(buf,0,n2,plan);
//--- change value
   a[0]=buf[0]+buf[1];
   for(i=1;i<=n2-1;i++)
     {
      //--- calculation
      idx=2*(i%n2);
      hn.re=buf[idx+0];
      hn.im=buf[idx+1];
      idx=2*(n2-i);
      hmnc.re=buf[idx+0];
      hmnc.im=-buf[idx+1];
      v.re=-MathSin(-(2*M_PI*i/n));
      v.im=MathCos(-(2*M_PI*i/n));
      v=hn+hmnc-v*(hn-hmnc);
      a[2*i+0]=0.5*v.re;
      a[2*i+1]=0.5*v.im;
     }
//--- change value
   a[1]=buf[0]-buf[1];
  }
//+------------------------------------------------------------------+
//| Internal subroutine. Never call it directly!                     |
//+------------------------------------------------------------------+
static void CFastFourierTransform::FFTR1DInvInternalEven(double &a[],const int n,
                                                         double &buf[],CFtPlan &plan)
  {
//--- create variables
   double x=0;
   double y=0;
   double t=0;
   int    i=0;
   int    n2=0;
//--- check
   if(!CAp::Assert(n>0 && n%2==0,__FUNCTION__+": incorrect N!"))
      return;
//--- Special cases:
//--- * N=2
//--- After this block we assume that N is strictly greater than 2
   if(n==2)
     {
      //--- change values
      x=0.5*(a[0]+a[1]);
      y=0.5*(a[0]-a[1]);
      a[0]=x;
      a[1]=y;
      //--- exit the function
      return;
     }
//--- inverse real FFT is reduced to the inverse real FHT,
//--- which is reduced to the forward real FHT,
//--- which is reduced to the forward real FFT.
//--- Don't worry,it is really compact and efficient reduction :)
   n2=n/2;
   buf[0]=a[0];
//--- calculation
   for(i=1;i<=n2-1;i++)
     {
      x=a[2*i+0];
      y=a[2*i+1];
      buf[i]=x-y;
      buf[n-i]=x+y;
     }
   buf[n2]=a[1];
//--- function call
   FFTR1DInternalEven(buf,n,a,plan);
//--- change values
   a[0]=buf[0]/n;
   t=1.0/(double)n;
//--- calculation
   for(i=1;i<=n2-1;i++)
     {
      x=buf[2*i+0];
      y=buf[2*i+1];
      a[i]=t*(x-y);
      a[n-i]=t*(x+y);
     }
//--- change value
   a[n2]=buf[1]/n;
  }
//+------------------------------------------------------------------+
//| Convolution class                                                |
//+------------------------------------------------------------------+
class CConv
  {
public:
   //--- constructor, destructor
                     CConv(void);
                    ~CConv(void);
   //--- methods
   static void       ConvC1D(al_complex &a[],const int m,al_complex &b[],const int n,al_complex &r[]);
   static void       ConvC1DInv(al_complex &a[],const int m,al_complex &b[],const int n,al_complex &r[]);
   static void       ConvC1DCircular(al_complex &s[],const int m,al_complex &r[],const int n,al_complex &c[]);
   static void       ConvC1DCircularInv(al_complex &a[],const int m,al_complex &b[],const int n,al_complex &r[]);
   static void       ConvR1D(double &a[],const int m,double &b[],const int n,double &r[]);
   static void       ConvR1DInv(double &a[],const int m,double &b[],const int n,double &r[]);
   static void       ConvR1DCircular(double &s[],const int m,double &r[],const int n,double &c[]);
   static void       ConvR1DCircularInv(double &a[],const int m,double &b[],const int n,double &r[]);
   static void       ConvC1DX(al_complex &a[],const int m,al_complex &b[],const int n,const bool circular,int alg,int q,al_complex &r[]);
   static void       ConvR1DX(double &a[],const int m,double &b[],const int n,const bool circular,int alg,int q,double &r[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CConv::CConv(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CConv::~CConv(void)
  {

  }
//+------------------------------------------------------------------+
//| 1-dimensional complex convolution.                               |
//| For given A/B returns conv(A,B) (non-circular). Subroutine can   |
//| automatically choose between three implementations:              |
//| straightforward O(M*N) formula for very small N (or M),          |
//| significantly larger than min(M,N), but O(M*N) algorithm is too  |
//| slow, and general FFT-based formula for cases where two previois |
//| algorithms are too slow.                                         |
//| Algorithm has max(M,N)*log(max(M,N)) complexity for any M/N.     |
//| INPUT PARAMETERS                                                 |
//|     A   -   array[0..M-1] - complex function to be transformed   |
//|     M   -   problem size                                         |
//|     B   -   array[0..N-1] - complex function to be transformed   |
//|     N   -   problem size                                         |
//| OUTPUT PARAMETERS                                                |
//|     R   -   convolution: A*B. array[0..N+M-2].                   |
//| NOTE:                                                            |
//|     It is assumed that A is zero at T<0, B is zero too. If one or|
//| both functions have non-zero values at negative T's, you can     |
//| still use this subroutine - just shift its result                |
//| correspondingly.                                                 |
//+------------------------------------------------------------------+
static void CConv::ConvC1D(al_complex &a[],const int m,al_complex &b[],const int n,
                           al_complex &r[])
  {
//--- check
   if(!CAp::Assert(n>0 && m>0,__FUNCTION__+": incorrect N or M!"))
      return;
//--- normalize task: make M>=N,
//--- so A will be longer that B.
   if(m<n)
     {
      //--- function call
      ConvC1D(b,n,a,m,r);
      return;
     }
//--- function call
   ConvC1DX(a,m,b,n,false,-1,0,r);
  }
//+------------------------------------------------------------------+
//| 1-dimensional complex non-circular deconvolution (inverse of     |
//| ConvC1D()).                                                      |
//| Algorithm has M*log(M)) complexity for any M (composite or prime)|
//| INPUT PARAMETERS                                                 |
//|     A   -   array[0..M-1] - convolved signal, A = conv(R, B)     |
//|     M   -   convolved signal length                              |
//|     B   -   array[0..N-1] - response                             |
//|     N   -   response length, N<=M                                |
//| OUTPUT PARAMETERS                                                |
//|     R   -   deconvolved signal. array[0..M-N].                   |
//| NOTE:                                                            |
//|     deconvolution is unstable process and may result in division |
//| by zero (if your response function is degenerate, i.e. has zero  |
//| Fourier coefficient).                                            |
//| NOTE:                                                            |
//|     It is assumed that A is zero at T<0, B is zero too. If one   |
//| or both functions have non-zero values at negative T's, you can  |
//| still use this subroutine - just shift its result correspondingly|
//+------------------------------------------------------------------+
static void CConv::ConvC1DInv(al_complex &a[],const int m,al_complex &b[],
                              const int n,al_complex &r[])
  {
//--- create variables
   int        i=0;
   int        p=0;
   al_complex c1=0;
   al_complex c2=0;
   al_complex c3=0;
   double     t=0;
//--- create arrays
   double buf[];
   double buf2[];
//--- object of class
   CFtPlan plan;
//--- check
   if(!CAp::Assert((n>0 && m>0) && n<=m,__FUNCTION__+": incorrect N or M!"))
      return;
//--- function call
   p=CFtBase::FtBaseFindSmooth(m);
   CFtBase::FtBaseGenerateComplexFFtPlan(p,plan);
//--- allocation
   ArrayResizeAL(buf,2*p);
//--- copy
   for(i=0;i<=m-1;i++)
     {
      buf[2*i+0]=a[i].re;
      buf[2*i+1]=a[i].im;
     }
//--- make zero
   for(i=m;i<=p-1;i++)
     {
      buf[2*i+0]=0;
      buf[2*i+1]=0;
     }
//--- allocation
   ArrayResizeAL(buf2,2*p);
//--- copy
   for(i=0;i<=n-1;i++)
     {
      buf2[2*i+0]=b[i].re;
      buf2[2*i+1]=b[i].im;
     }
//--- make zero
   for(i=n;i<=p-1;i++)
     {
      buf2[2*i+0]=0;
      buf2[2*i+1]=0;
     }
//--- function call
   CFtBase::FtBaseExecutePlan(buf,0,p,plan);
   CFtBase::FtBaseExecutePlan(buf2,0,p,plan);
//--- calculation
   for(i=0;i<=p-1;i++)
     {
      c1.re=buf[2*i+0];
      c1.im=buf[2*i+1];
      c2.re=buf2[2*i+0];
      c2.im=buf2[2*i+1];
      c3=c1/c2;
      buf[2*i+0]=c3.re;
      buf[2*i+1]=-c3.im;
     }
//--- function call
   CFtBase::FtBaseExecutePlan(buf,0,p,plan);
   t=1.0/(double)p;
//--- allocation
   ArrayResizeAL(r,m-n+1);
//--- change values
   for(i=0;i<=m-n;i++)
     {
      r[i].re=t*buf[2*i+0];
      r[i].im=-(t*buf[2*i+1]);
     }
  }
//+------------------------------------------------------------------+
//| 1-dimensional circular complex convolution.                      |
//| For given S/R returns conv(S,R) (circular). Algorithm has        |
//| linearithmic complexity for any M/N.                             |
//| IMPORTANT:  normal convolution is commutative, i.e. it is        |
//| symmetric - conv(A,B)=conv(B,A). Cyclic convolution IS NOT. One  |
//| function - S - is a signal, periodic function, and another - R - |
//| is a response, non-periodic function with limited length.        |
//| INPUT PARAMETERS                                                 |
//|     S   -   array[0..M-1] - complex periodic signal              |
//|     M   -   problem size                                         |
//|     B   -   array[0..N-1] - complex non-periodic response        |
//|     N   -   problem size                                         |
//| OUTPUT PARAMETERS                                                |
//|     R   -   convolution: A*B. array[0..M-1].                     |
//| NOTE:                                                            |
//|     It is assumed that B is zero at T<0. If it has non-zero      |
//| values at negative T's, you can still use this subroutine - just |
//| shift its result correspondingly.                                |
//+------------------------------------------------------------------+
static void CConv::ConvC1DCircular(al_complex &s[],const int m,al_complex &r[],
                                   const int n,al_complex &c[])
  {
//--- create variables
   int i1=0;
   int i2=0;
   int j2=0;
   int i_=0;
   int i1_=0;
//--- create array
   al_complex buf[];
//--- check
   if(!CAp::Assert(n>0 && m>0,__FUNCTION__+": incorrect N or M!"))
      return;
//--- normalize task: make M>=N,
//--- so A will be longer (at least - not shorter) that B.
   if(m<n)
     {
      //--- allocation
      ArrayResizeAL(buf,m);
      //--- initialization
      for(i1=0;i1<=m-1;i1++)
         buf[i1]=0;
      i1=0;
      //--- calculation
      while(i1<n)
        {
         //--- change values
         i2=MathMin(i1+m-1,n-1);
         j2=i2-i1;
         i1_=i1;
         for(i_=0;i_<=j2;i_++)
            buf[i_]=buf[i_]+r[i_+i1_];
         i1=i1+m;
        }
      //--- function call
      ConvC1DCircular(s,m,buf,m,c);
      //--- exit the function
      return;
     }
//--- function call
   ConvC1DX(s,m,r,n,true,-1,0,c);
  }
//+------------------------------------------------------------------+
//| 1-dimensional circular complex deconvolution (inverse of         |
//| ConvC1DCircular()).                                              |
//| Algorithm has M*log(M)) complexity for any M (composite or prime)|
//| INPUT PARAMETERS                                                 |
//|     A   -   array[0..M-1] - convolved periodic signal,           |
//|             A = conv(R, B)                                       |
//|     M   -   convolved signal length                              |
//|     B   -   array[0..N-1] - non-periodic response                |
//|     N   -   response length                                      |
//| OUTPUT PARAMETERS                                                |
//|     R   -   deconvolved signal. array[0..M-1].                   |
//| NOTE:                                                            |
//|     deconvolution is unstable process and may result in division |
//| by zero (if your response function is degenerate, i.e. has zero  |
//| Fourier coefficient).                                            |
//| NOTE:                                                            |
//|     It is assumed that B is zero at T<0. If it has non-zero      |
//| values at negative T's, you can still use this subroutine - just |
//| shift its result correspondingly.                                |
//+------------------------------------------------------------------+
static void CConv::ConvC1DCircularInv(al_complex &a[],const int m,al_complex &b[],
                                      const int n,al_complex &r[])
  {
//--- create variables
   int        i=0;
   int        i1=0;
   int        i2=0;
   int        j2=0;
   al_complex c1=0;
   al_complex c2=0;
   al_complex c3=0;
   double     t=0;
   int        i_=0;
   int        i1_=0;
//--- create arrays
   double     buf[];
   double     buf2[];
   al_complex cbuf[];
//--- object of class
   CFtPlan plan;
//--- check
   if(!CAp::Assert(n>0 && m>0,__FUNCTION__+": incorrect N or M!"))
      return;
//--- normalize task: make M>=N,
//--- so A will be longer (at least - not shorter) that B.
   if(m<n)
     {
      //--- allocation
      ArrayResizeAL(cbuf,m);
      //--- initialization
      for(i=0;i<=m-1;i++)
         cbuf[i]=0;
      i1=0;
      //--- calculation
      while(i1<n)
        {
         //--- change values
         i2=MathMin(i1+m-1,n-1);
         j2=i2-i1;
         i1_=i1;
         for(i_=0;i_<=j2;i_++)
            cbuf[i_]=cbuf[i_]+b[i_+i1_];
         i1=i1+m;
        }
      //--- function call
      ConvC1DCircularInv(a,m,cbuf,m,r);
      //--- exit the function
      return;
     }
//--- Task is normalized
   CFtBase::FtBaseGenerateComplexFFtPlan(m,plan);
//--- allocation
   ArrayResizeAL(buf,2*m);
//--- copy
   for(i=0;i<=m-1;i++)
     {
      buf[2*i+0]=a[i].re;
      buf[2*i+1]=a[i].im;
     }
//--- allocation
   ArrayResizeAL(buf2,2*m);
//--- copy
   for(i=0;i<=n-1;i++)
     {
      buf2[2*i+0]=b[i].re;
      buf2[2*i+1]=b[i].im;
     }
//--- make zero
   for(i=n;i<=m-1;i++)
     {
      buf2[2*i+0]=0;
      buf2[2*i+1]=0;
     }
//--- function call
   CFtBase::FtBaseExecutePlan(buf,0,m,plan);
   CFtBase::FtBaseExecutePlan(buf2,0,m,plan);
//--- calculation
   for(i=0;i<=m-1;i++)
     {
      c1.re=buf[2*i+0];
      c1.im=buf[2*i+1];
      c2.re=buf2[2*i+0];
      c2.im=buf2[2*i+1];
      c3=c1/c2;
      buf[2*i+0]=c3.re;
      buf[2*i+1]=-c3.im;
     }
//--- function call
   CFtBase::FtBaseExecutePlan(buf,0,m,plan);
   t=1.0/(double)m;
//--- allocation
   ArrayResizeAL(r,m);
//--- change values
   for(i=0;i<=m-1;i++)
     {
      r[i].re=t*buf[2*i+0];
      r[i].im=-(t*buf[2*i+1]);
     }
  }
//+------------------------------------------------------------------+
//| 1-dimensional real convolution.                                  |
//| Analogous to ConvC1D(), see ConvC1D() comments for more details. |
//| INPUT PARAMETERS                                                 |
//|     A   -   array[0..M-1] - real function to be transformed      |
//|     M   -   problem size                                         |
//|     B   -   array[0..N-1] - real function to be transformed      |
//|     N   -   problem size                                         |
//| OUTPUT PARAMETERS                                                |
//|     R   -   convolution: A*B. array[0..N+M-2].                   |
//| NOTE:                                                            |
//|     It is assumed that A is zero at T<0, B is zero too. If one   |
//| or both functions have non-zero values at negative T's, you can  |
//| still use this subroutine - just shift its result correspondingly|
//+------------------------------------------------------------------+
static void CConv::ConvR1D(double &a[],const int m,double &b[],
                           const int n,double &r[])
  {
//--- check
   if(!CAp::Assert(n>0 && m>0,__FUNCTION__+": incorrect N or M!"))
      return;
//--- normalize task: make M>=N,
//--- so A will be longer that B.
   if(m<n)
     {
      //--- function call
      ConvR1D(b,n,a,m,r);
      return;
     }
//--- function call
   ConvR1DX(a,m,b,n,false,-1,0,r);
  }
//+------------------------------------------------------------------+
//| 1-dimensional real deconvolution (inverse of ConvC1D()).         |
//| Algorithm has M*log(M)) complexity for any M (composite or prime)|
//| INPUT PARAMETERS                                                 |
//|     A   -   array[0..M-1] - convolved signal, A = conv(R, B)     |
//|     M   -   convolved signal length                              |
//|     B   -   array[0..N-1] - response                             |
//|     N   -   response length, N<=M                                |
//| OUTPUT PARAMETERS                                                |
//|     R   -   deconvolved signal. array[0..M-N].                   |
//| NOTE:                                                            |
//|     deconvolution is unstable process and may result in division |
//| by zero (if your response function is degenerate, i.e. has zero  |
//| Fourier coefficient).                                            |
//| NOTE:                                                            |
//|     It is assumed that A is zero at T<0, B is zero too. If one or|
//| both functions have non-zero values at negative T's, you can     |
//| still use this subroutine - just shift its result correspondingly|
//+------------------------------------------------------------------+
static void CConv::ConvR1DInv(double &a[],const int m,double &b[],
                              const int n,double &r[])
  {
//--- create variables
   int        i=0;
   int        p=0;
   al_complex c1=0;
   al_complex c2=0;
   al_complex c3=0;
   int        i_=0;
//--- create arrays
   double buf[];
   double buf2[];
   double buf3[];
//--- object of class
   CFtPlan plan;
//--- check
   if(!CAp::Assert((n>0 && m>0) && n<=m,__FUNCTION__+": incorrect N or M!"))
      return;
//--- function call
   p=CFtBase::FtBaseFindSmoothEven(m);
//--- allocation
   ArrayResizeAL(buf,p);
//--- copy
   for(i_=0;i_<=m-1;i_++)
      buf[i_]=a[i_];
//--- make zero
   for(i=m;i<=p-1;i++)
      buf[i]=0;
//--- allocation
   ArrayResizeAL(buf2,p);
   for(i_=0;i_<=n-1;i_++)
      buf2[i_]=b[i_];
//--- make zero
   for(i=n;i<=p-1;i++)
      buf2[i]=0;
//--- allocation
   ArrayResizeAL(buf3,p);
//--- function call
   CFtBase::FtBaseGenerateComplexFFtPlan(p/2,plan);
//--- function call
   CFastFourierTransform::FFTR1DInternalEven(buf,p,buf3,plan);
//--- function call
   CFastFourierTransform::FFTR1DInternalEven(buf2,p,buf3,plan);
//--- change values
   buf[0]=buf[0]/buf2[0];
   buf[1]=buf[1]/buf2[1];
//--- calculation
   for(i=1;i<=p/2-1;i++)
     {
      c1.re=buf[2*i+0];
      c1.im=buf[2*i+1];
      c2.re=buf2[2*i+0];
      c2.im=buf2[2*i+1];
      c3=c1/c2;
      buf[2*i+0]=c3.re;
      buf[2*i+1]=c3.im;
     }
//--- function call
   CFastFourierTransform::FFTR1DInvInternalEven(buf,p,buf3,plan);
//--- allocation
   ArrayResizeAL(r,m-n+1);
//--- copy
   for(i_=0;i_<=m-n;i_++)
      r[i_]=buf[i_];
  }
//+------------------------------------------------------------------+
//| 1-dimensional circular real convolution.                         |
//| Analogous to ConvC1DCircular(), see ConvC1DCircular() comments   |
//| for more details.                                                |
//| INPUT PARAMETERS                                                 |
//|     S   -   array[0..M-1] - real signal                          |
//|     M   -   problem size                                         |
//|     B   -   array[0..N-1] - real response                        |
//|     N   -   problem size                                         |
//| OUTPUT PARAMETERS                                                |
//|     R   -   convolution: A*B. array[0..M-1].                     |
//| NOTE:                                                            |
//|     It is assumed that B is zero at T<0. If it has non-zero      |
//| values at negative T's, you can still use this subroutine - just |
//| shift its result correspondingly.                                |
//+------------------------------------------------------------------+
static void CConv::ConvR1DCircular(double &s[],const int m,double &r[],
                                   const int n,double &c[])
  {
//--- create variables
   int i1=0;
   int i2=0;
   int j2=0;
   int i_=0;
   int i1_=0;
//--- create array
   double buf[];
//--- check
   if(!CAp::Assert(n>0 && m>0,__FUNCTION__+": incorrect N or M!"))
      return;
//--- normalize task: make M>=N,
//--- so A will be longer (at least - not shorter) that B.
   if(m<n)
     {
      //--- allocation
      ArrayResizeAL(buf,m);
      //--- initialization
      for(i1=0;i1<=m-1;i1++)
         buf[i1]=0;
      i1=0;
      //--- calculation
      while(i1<n)
        {
         i2=MathMin(i1+m-1,n-1);
         j2=i2-i1;
         i1_=i1;
         for(i_=0;i_<=j2;i_++)
            buf[i_]=buf[i_]+r[i_+i1_];
         i1=i1+m;
        }
      ConvR1DCircular(s,m,buf,m,c);
      //--- exit the function
      return;
     }
//--- reduce to usual convolution
   ConvR1DX(s,m,r,n,true,-1,0,c);
  }
//+------------------------------------------------------------------+
//| 1-dimensional complex deconvolution (inverse of ConvC1D()).      |
//| Algorithm has M*log(M)) complexity for any M (composite or prime)|
//| INPUT PARAMETERS                                                 |
//|     A   -   array[0..M-1] - convolved signal, A = conv(R, B)     |
//|     M   -   convolved signal length                              |
//|     B   -   array[0..N-1] - response                             |
//|     N   -   response length                                      |
//| OUTPUT PARAMETERS                                                |
//|     R   -   deconvolved signal. array[0..M-N].                   |
//| NOTE:                                                            |
//|     deconvolution is unstable process and may result in division |
//| by zero (if your response function is degenerate, i.e. has zero  |
//| Fourier coefficient).                                            |
//| NOTE:                                                            |
//|     It is assumed that B is zero at T<0. If it has non-zero      |
//| values at negative T's, you can still use this subroutine - just |
//| shift its result correspondingly.                                |
//+------------------------------------------------------------------+
static void CConv::ConvR1DCircularInv(double &a[],const int m,double &b[],
                                      const int n,double &r[])
  {
//--- create variables
   int        i=0;
   int        i1=0;
   int        i2=0;
   int        j2=0;
   al_complex c1=0;
   al_complex c2=0;
   al_complex c3=0;
   int        i_=0;
   int        i1_=0;
//--- create arrays
   double     buf[];
   double     buf2[];
   double     buf3[];
   al_complex cbuf[];
   al_complex cbuf2[];
//--- object of class
   CFtPlan plan;
//--- check
   if(!CAp::Assert(n>0 && m>0,__FUNCTION__+": incorrect N or M!"))
      return;
//--- normalize task: make M>=N,
//--- so A will be longer (at least - not shorter) that B.
   if(m<n)
     {
      //--- allocation
      ArrayResizeAL(buf,m);
      //--- initialization
      for(i=0;i<=m-1;i++)
         buf[i]=0;
      i1=0;
      //--- calculation
      while(i1<n)
        {
         i2=MathMin(i1+m-1,n-1);
         j2=i2-i1;
         i1_=i1;
         for(i_=0;i_<=j2;i_++)
            buf[i_]=buf[i_]+b[i_+i1_];
         i1=i1+m;
        }
      //--- function call
      ConvR1DCircularInv(a,m,buf,m,r);
      //--- exit the function
      return;
     }
//--- Task is normalized
   if(m%2==0)
     {
      //--- size is even,use fast even-size FFT
      ArrayResizeAL(buf,m);
      for(i_=0;i_<=m-1;i_++)
         buf[i_]=a[i_];
      //--- allocation
      ArrayResizeAL(buf2,m);
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         buf2[i_]=b[i_];
      //--- make zero
      for(i=n;i<=m-1;i++)
         buf2[i]=0;
      //--- allocation
      ArrayResizeAL(buf3,m);
      //--- function call
      CFtBase::FtBaseGenerateComplexFFtPlan(m/2,plan);
      //--- function call
      CFastFourierTransform::FFTR1DInternalEven(buf,m,buf3,plan);
      //--- function call
      CFastFourierTransform::FFTR1DInternalEven(buf2,m,buf3,plan);
      //--- change values
      buf[0]=buf[0]/buf2[0];
      buf[1]=buf[1]/buf2[1];
      //--- calculation
      for(i=1;i<=m/2-1;i++)
        {
         c1.re=buf[2*i+0];
         c1.im=buf[2*i+1];
         c2.re=buf2[2*i+0];
         c2.im=buf2[2*i+1];
         c3=c1/c2;
         buf[2*i+0]=c3.re;
         buf[2*i+1]=c3.im;
        }
      //--- function call
      CFastFourierTransform::FFTR1DInvInternalEven(buf,m,buf3,plan);
      //--- allocation
      ArrayResizeAL(r,m);
      //--- copy
      for(i_=0;i_<=m-1;i_++)
         r[i_]=buf[i_];
     }
   else
     {
      //--- odd-size,use general real FFT
      CFastFourierTransform::FFTR1D(a,m,cbuf);
      //--- allocation
      ArrayResizeAL(buf2,m);
      //--- copy
      for(i_=0;i_<=n-1;i_++)
         buf2[i_]=b[i_];
      //--- initialization
      for(i=n;i<=m-1;i++)
         buf2[i]=0;
      //--- function call
      CFastFourierTransform::FFTR1D(buf2,m,cbuf2);
      //--- calculation
      for(i=0;i<=(int)MathFloor((double)m/2.0);i++)
         cbuf[i]=cbuf[i]/cbuf2[i];
      //--- function call
      CFastFourierTransform::FFTR1DInv(cbuf,m,r);
     }
  }
//+------------------------------------------------------------------+
//| 1-dimensional complex convolution.                               |
//| Extended subroutine which allows to choose convolution algorithm.|
//| Intended for internal use, ALGLIB users should call              |
//| ConvC1D()/ConvC1DCircular().                                     |
//| INPUT PARAMETERS                                                 |
//|     A   -   array[0..M-1] - complex function to be transformed   |
//|     M   -   problem size                                         |
//|     B   -   array[0..N-1] - complex function to be transformed   |
//|     N   -   problem size, N<=M                                   |
//|     Alg -   algorithm type:                                      |
//|             *-2     auto-select Q for overlap-add                |
//|             *-1     auto-select algorithm and parameters         |
//|             * 0     straightforward formula for small N's        |
//|             * 1     general FFT-based code                       |
//|             * 2     overlap-add with length Q                    |
//|     Q   -   length for overlap-add                               |
//| OUTPUT PARAMETERS                                                |
//|     R   -   convolution: A*B. array[0..N+M-1].                   |
//+------------------------------------------------------------------+
static void CConv::ConvC1DX(al_complex &a[],const int m,al_complex &b[],const int n,
                            const bool circular,int alg,int q,al_complex &r[])
  {
//--- create variables
   int        i=0;
   int        j=0;
   int        p=0;
   int        ptotal=0;
   int        i1=0;
   int        i2=0;
   int        j1=0;
   int        j2=0;
   al_complex v=0;
   double     ax=0;
   double     ay=0;
   double     bx=0;
   double     by=0;
   double     t=0;
   double     tx=0;
   double     ty=0;
   double     flopcand=0;
   double     flopbest=0;
   int        algbest=0;
   int        i_=0;
   int        i1_=0;
//--- create arrays
   al_complex bbuf[];
   double     buf[];
   double     buf2[];
//--- object of class
   CFtPlan plan;
//--- check
   if(!CAp::Assert(n>0 && m>0,__FUNCTION__+": incorrect N or M!"))
      return;
//--- check
   if(!CAp::Assert(n<=m,__FUNCTION__+": N<M assumption is false!"))
      return;
//--- Auto-select
   if(alg==-1 || alg==-2)
     {
      //--- Initial candidate: straightforward implementation.
      //--- If we want to use auto-fitted overlap-add,
      //--- flop count is initialized by large real number - to force
      //--- another algorithm selection
      algbest=0;
      //--- check
      if(alg==-1)
         flopbest=2*m*n;
      else
         flopbest=CMath::m_maxrealnumber;
      //--- Another candidate - generic FFT code
      if(alg==-1)
        {
         //--- check
         if(circular && CFtBase::FtBaseIsSmooth(m))
           {
            //--- special code for circular convolution of a sequence with a smooth length
            flopcand=3*CFtBase::FtBaseGetFlopEstimate(m)+6*m;
            //--- check
            if(flopcand<flopbest)
              {
               algbest=1;
               flopbest=flopcand;
              }
           }
         else
           {
            //--- general cyclic/non-cyclic convolution
            p=CFtBase::FtBaseFindSmooth(m+n-1);
            flopcand=3*CFtBase::FtBaseGetFlopEstimate(p)+6*p;
            //--- check
            if(flopcand<flopbest)
              {
               algbest=1;
               flopbest=flopcand;
              }
           }
        }
      //--- Another candidate - overlap-add
      q=1;
      ptotal=1;
      while(ptotal<n)
         ptotal=ptotal*2;
      //--- calculation
      while(ptotal<=m+n-1)
        {
         //--- change values
         p=ptotal-n+1;
         flopcand=(int)MathCeil((double)m/(double)p)*(2*CFtBase::FtBaseGetFlopEstimate(ptotal)+8*ptotal);
         //--- check
         if(flopcand<flopbest)
           {
            flopbest=flopcand;
            algbest=2;
            q=p;
           }
         ptotal=ptotal*2;
        }
      //--- change value
      alg=algbest;
      //--- function call
      ConvC1DX(a,m,b,n,circular,alg,q,r);
      //--- exit the function
      return;
     }
//--- straightforward formula for
//--- circular and non-circular convolutions.
//--- Very simple code,no further comments needed.
   if(alg==0)
     {
      //--- Special case: N=1
      if(n==1)
        {
         //--- allocation
         ArrayResizeAL(r,m);
         v=b[0];
         for(i_=0;i_<=m-1;i_++)
            r[i_]=v*a[i_];
         //--- exit the function
         return;
        }
      //--- use straightforward formula
      if(circular)
        {
         //--- circular convolution
         ArrayResizeAL(r,m);
         //--- change values
         v=b[0];
         for(i_=0;i_<=m-1;i_++)
            r[i_]=v*a[i_];
         //--- calculation
         for(i=1;i<=n-1;i++)
           {
            //--- change values
            v=b[i];
            i1=0;
            i2=i-1;
            j1=m-i;
            j2=m-1;
            i1_=j1-i1;
            //--- calculation
            for(i_=i1;i_<=i2;i_++)
               r[i_]=r[i_]+v*a[i_+i1_];
            //--- change values
            i1=i;
            i2=m-1;
            j1=0;
            j2=m-i-1;
            i1_=j1-i1;
            //--- calculation
            for(i_=i1;i_<=i2;i_++)
               r[i_]=r[i_]+v*a[i_+i1_];
           }
        }
      else
        {
         //--- non-circular convolution
         ArrayResizeAL(r,m+n-1);
         //--- initialization
         for(i=0;i<=m+n-2;i++)
            r[i]=0;
         //--- calculation
         for(i=0;i<=n-1;i++)
           {
            v=b[i];
            i1_=-i;
            for(i_=i;i_<=i+m-1;i_++)
               r[i_]=r[i_]+v*a[i_+i1_];
           }
        }
      //--- exit the function
      return;
     }
//--- general FFT-based code for
//--- circular and non-circular convolutions.
//--- First,if convolution is circular,we test whether M is smooth or not.
//--- If it is smooth,we just use M-length FFT to calculate convolution.
//--- If it is not,we calculate non-circular convolution and wrap it arount.
//--- IF convolution is non-circular,we use zero-padding + FFT.
   if(alg==1)
     {
      //--- check
      if(circular && CFtBase::FtBaseIsSmooth(m))
        {
         //--- special code for circular convolution with smooth M
         CFtBase::FtBaseGenerateComplexFFtPlan(m,plan);
         //--- allocation
         ArrayResizeAL(buf,2*m);
         //--- copy
         for(i=0;i<=m-1;i++)
           {
            buf[2*i+0]=a[i].re;
            buf[2*i+1]=a[i].im;
           }
         //--- allocation
         ArrayResizeAL(buf2,2*m);
         //--- copy
         for(i=0;i<=n-1;i++)
           {
            buf2[2*i+0]=b[i].re;
            buf2[2*i+1]=b[i].im;
           }
         //--- make zero
         for(i=n;i<=m-1;i++)
           {
            buf2[2*i+0]=0;
            buf2[2*i+1]=0;
           }
         //--- function call
         CFtBase::FtBaseExecutePlan(buf,0,m,plan);
         CFtBase::FtBaseExecutePlan(buf2,0,m,plan);
         //--- calculation
         for(i=0;i<=m-1;i++)
           {
            ax=buf[2*i+0];
            ay=buf[2*i+1];
            bx=buf2[2*i+0];
            by=buf2[2*i+1];
            tx=ax*bx-ay*by;
            ty=ax*by+ay*bx;
            buf[2*i+0]=tx;
            buf[2*i+1]=-ty;
           }
         //--- function call
         CFtBase::FtBaseExecutePlan(buf,0,m,plan);
         t=1.0/(double)m;
         //--- allocation
         ArrayResizeAL(r,m);
         //--- change values
         for(i=0;i<=m-1;i++)
           {
            r[i].re=t*buf[2*i+0];
            r[i].im=-(t*buf[2*i+1]);
           }
        }
      else
        {
         //--- M is non-smooth,general code (circular/non-circular):
         //--- * first part is the same for circular and non-circular
         //---   convolutions. zero padding,FFTs,inverse FFTs
         //--- * second part differs:
         //---   * for non-circular convolution we just copy array
         //---   * for circular convolution we add array tail to its head
         p=CFtBase::FtBaseFindSmooth(m+n-1);
         CFtBase::FtBaseGenerateComplexFFtPlan(p,plan);
         //--- allocation
         ArrayResizeAL(buf,2*p);
         //--- copy
         for(i=0;i<=m-1;i++)
           {
            buf[2*i+0]=a[i].re;
            buf[2*i+1]=a[i].im;
           }
         //--- make zero
         for(i=m;i<=p-1;i++)
           {
            buf[2*i+0]=0;
            buf[2*i+1]=0;
           }
         //--- allocation
         ArrayResizeAL(buf2,2*p);
         //--- copy
         for(i=0;i<=n-1;i++)
           {
            buf2[2*i+0]=b[i].re;
            buf2[2*i+1]=b[i].im;
           }
         //--- make zero
         for(i=n;i<=p-1;i++)
           {
            buf2[2*i+0]=0;
            buf2[2*i+1]=0;
           }
         //--- function call
         CFtBase::FtBaseExecutePlan(buf,0,p,plan);
         CFtBase::FtBaseExecutePlan(buf2,0,p,plan);
         //--- calculation
         for(i=0;i<=p-1;i++)
           {
            ax=buf[2*i+0];
            ay=buf[2*i+1];
            bx=buf2[2*i+0];
            by=buf2[2*i+1];
            tx=ax*bx-ay*by;
            ty=ax*by+ay*bx;
            buf[2*i+0]=tx;
            buf[2*i+1]=-ty;
           }
         //--- function call
         CFtBase::FtBaseExecutePlan(buf,0,p,plan);
         t=1.0/(double)p;
         //--- check
         if(circular)
           {
            //--- circular,add tail to head
            ArrayResizeAL(r,m);
            //--- copy
            for(i=0;i<=m-1;i++)
              {
               r[i].re=t*buf[2*i+0];
               r[i].im=-(t*buf[2*i+1]);
              }
            //--- change values
            for(i=m;i<=m+n-2;i++)
              {
               r[i-m].re=r[i-m].re+t*buf[2*i+0];
               r[i-m].im=r[i-m].im-t*buf[2*i+1];
              }
           }
         else
           {
            //--- non-circular,just copy
            ArrayResizeAL(r,m+n-1);
            for(i=0;i<=m+n-2;i++)
              {
               r[i].re=t*buf[2*i+0];
               r[i].im=-(t*buf[2*i+1]);
              }
           }
        }
      //--- exit the function
      return;
     }
//--- overlap-add method for
//--- circular and non-circular convolutions.
//--- First part of code (separate FFTs of input blocks) is the same
//--- for all types of convolution. Second part (overlapping outputs)
//--- differs for different types of convolution. We just copy output
//--- when convolution is non-circular. We wrap it around,if it is
//--- circular.
   if(alg==2)
     {
      //--- allocation
      ArrayResizeAL(buf,2*(q+n-1));
      //--- prepare R
      if(circular)
        {
         //--- allocation
         ArrayResizeAL(r,m);
         for(i=0;i<=m-1;i++)
            r[i]=0;
        }
      else
        {
         //--- allocation
         ArrayResizeAL(r,m+n-1);
         for(i=0;i<=m+n-2;i++)
            r[i]=0;
        }
      //--- pre-calculated FFT(B)
      ArrayResizeAL(bbuf,q+n-1);
      for(i_=0;i_<=n-1;i_++)
         bbuf[i_]=b[i_];
      for(j=n;j<=q+n-2;j++)
         bbuf[j]=0;
      //--- function call
      CFastFourierTransform::FFTC1D(bbuf,q+n-1);
      //--- prepare FFT plan for chunks of A
      CFtBase::FtBaseGenerateComplexFFtPlan(q+n-1,plan);
      //--- main overlap-add cycle
      i=0;
      //--- calculation
      while(i<=m-1)
        {
         p=MathMin(q,m-i);
         //--- copy
         for(j=0;j<=p-1;j++)
           {
            buf[2*j+0]=a[i+j].re;
            buf[2*j+1]=a[i+j].im;
           }
         //--- make zero
         for(j=p;j<=q+n-2;j++)
           {
            buf[2*j+0]=0;
            buf[2*j+1]=0;
           }
         //--- function call
         CFtBase::FtBaseExecutePlan(buf,0,q+n-1,plan);
         //--- calculation
         for(j=0;j<=q+n-2;j++)
           {
            ax=buf[2*j+0];
            ay=buf[2*j+1];
            bx=bbuf[j].re;
            by=bbuf[j].im;
            tx=ax*bx-ay*by;
            ty=ax*by+ay*bx;
            buf[2*j+0]=tx;
            buf[2*j+1]=-ty;
           }
         //--- function call
         CFtBase::FtBaseExecutePlan(buf,0,q+n-1,plan);
         t=1.0/(double)(q+n-1);
         //--- check
         if(circular)
           {
            j1=MathMin(i+p+n-2,m-1)-i;
            j2=j1+1;
           }
         else
           {
            j1=p+n-2;
            j2=j1+1;
           }
         //--- change values
         for(j=0;j<=j1;j++)
           {
            r[i+j].re=r[i+j].re+buf[2*j+0]*t;
            r[i+j].im=r[i+j].im-buf[2*j+1]*t;
           }
         //--- change values
         for(j=j2;j<=p+n-2;j++)
           {
            r[j-j2].re=r[j-j2].re+buf[2*j+0]*t;
            r[j-j2].im=r[j-j2].im-buf[2*j+1]*t;
           }
         i=i+p;
        }
      //--- exit the function
      return;
     }
  }
//+------------------------------------------------------------------+
//| 1-dimensional real convolution.                                  |
//| Extended subroutine which allows to choose convolution algorithm.|
//| Intended for internal use, ALGLIB users should call ConvR1D().   |
//| INPUT PARAMETERS                                                 |
//|     A   -   array[0..M-1] - complex function to be transformed   |
//|     M   -   problem size                                         |
//|     B   -   array[0..N-1] - complex function to be transformed   |
//|     N   -   problem size, N<=M                                   |
//|     Alg -   algorithm type:                                      |
//|             *-2     auto-select Q for overlap-add                |
//|             *-1     auto-select algorithm and parameters         |
//|             * 0     straightforward formula for small N's        |
//|             * 1     general FFT-based code                       |
//|             * 2     overlap-add with length Q                    |
//|     Q   -   length for overlap-add                               |
//| OUTPUT PARAMETERS                                                |
//|     R   -   convolution: A*B. array[0..N+M-1].                   |
//+------------------------------------------------------------------+
static void CConv::ConvR1DX(double &a[],const int m,double &b[],const int n,
                            const bool circular,int alg,int q,double &r[])
  {
//--- create variables
   double v=0;
   int    i=0;
   int    j=0;
   int    p=0;
   int    ptotal=0;
   int    i1=0;
   int    i2=0;
   int    j1=0;
   int    j2=0;
   double ax=0;
   double ay=0;
   double bx=0;
   double by=0;
   double tx=0;
   double ty=0;
   double flopcand=0;
   double flopbest=0;
   int    algbest=0;
   int    i_=0;
   int    i1_=0;
//--- create arrays
   double buf[];
   double buf2[];
   double buf3[];
//--- object of class
   CFtPlan plan;
//--- check
   if(!CAp::Assert(n>0 && m>0,__FUNCTION__+": incorrect N or M!"))
      return;
//--- check
   if(!CAp::Assert(n<=m,__FUNCTION__+": N<M assumption is false!"))
      return;
//--- handle special cases
   if(MathMin(m,n)<=2)
      alg=0;
//--- Auto-select
   if(alg<0)
     {
      //--- Initial candidate: straightforward implementation.
      //--- If we want to use auto-fitted overlap-add,
      //--- flop count is initialized by large real number - to force
      //--- another algorithm selection
      algbest=0;
      //--- check
      if(alg==-1)
         flopbest=0.15*m*n;
      else
         flopbest=CMath::m_maxrealnumber;
      //--- Another candidate - generic FFT code
      if(alg==-1)
        {
         //--- check
         if((circular && CFtBase::FtBaseIsSmooth(m)) && m%2==0)
           {
            //--- special code for circular convolution of a sequence with a smooth length
            flopcand=3*CFtBase::FtBaseGetFlopEstimate(m/2)+(double)(6*m)/2.0;
            //--- check
            if(flopcand<flopbest)
              {
               algbest=1;
               flopbest=flopcand;
              }
           }
         else
           {
            //--- general cyclic/non-cyclic convolution
            p=CFtBase::FtBaseFindSmoothEven(m+n-1);
            flopcand=3*CFtBase::FtBaseGetFlopEstimate(p/2)+(double)(6*p)/2.0;
            //--- check
            if(flopcand<flopbest)
              {
               algbest=1;
               flopbest=flopcand;
              }
           }
        }
      //--- Another candidate - overlap-add
      q=1;
      ptotal=1;
      while(ptotal<n)
         ptotal=ptotal*2;
      //--- calculation
      while(ptotal<=m+n-1)
        {
         //--- change values
         p=ptotal-n+1;
         flopcand=(int)MathCeil((double)m/(double)p)*(2*CFtBase::FtBaseGetFlopEstimate(ptotal/2)+1*(ptotal/2));
         //--- check
         if(flopcand<flopbest)
           {
            flopbest=flopcand;
            algbest=2;
            q=p;
           }
         ptotal=ptotal*2;
        }
      alg=algbest;
      //--- function call
      ConvR1DX(a,m,b,n,circular,alg,q,r);
      //--- exit the function
      return;
     }
//--- straightforward formula for
//--- circular and non-circular convolutions.
//--- Very simple code,no further comments needed.
   if(alg==0)
     {
      //--- Special case: N=1
      if(n==1)
        {
         //--- allocation
         ArrayResizeAL(r,m);
         v=b[0];
         for(i_=0;i_<=m-1;i_++)
            r[i_]=v*a[i_];
         //--- exit the function
         return;
        }
      //--- use straightforward formula
      if(circular)
        {
         //--- circular convolution
         ArrayResizeAL(r,m);
         v=b[0];
         for(i_=0;i_<=m-1;i_++)
            r[i_]=v*a[i_];
         //--- calculation
         for(i=1;i<=n-1;i++)
           {
            //--- change values
            v=b[i];
            i1=0;
            i2=i-1;
            j1=m-i;
            j2=m-1;
            i1_=j1-i1;
            //--- calculation
            for(i_=i1;i_<=i2;i_++)
               r[i_]=r[i_]+v*a[i_+i1_];
            //--- change values
            i1=i;
            i2=m-1;
            j1=0;
            j2=m-i-1;
            i1_=j1-i1;
            //--- calculation
            for(i_=i1;i_<=i2;i_++)
               r[i_]=r[i_]+v*a[i_+i1_];
           }
        }
      else
        {
         //--- non-circular convolution
         ArrayResizeAL(r,m+n-1);
         for(i=0;i<=m+n-2;i++)
            r[i]=0;
         //--- calculation
         for(i=0;i<=n-1;i++)
           {
            v=b[i];
            i1_=-i;
            for(i_=i;i_<=i+m-1;i_++)
               r[i_]=r[i_]+v*a[i_+i1_];
           }
        }
      //--- exit the function
      return;
     }
//--- general FFT-based code for
//--- circular and non-circular convolutions.
//--- First,if convolution is circular,we test whether M is smooth or not.
//--- If it is smooth,we just use M-length FFT to calculate convolution.
//--- If it is not,we calculate non-circular convolution and wrap it arount.
//--- If convolution is non-circular,we use zero-padding + FFT.
//--- We assume that M+N-1>2 - we should call small case code otherwise
   if(alg==1)
     {
      //--- check
      if(!CAp::Assert(m+n-1>2,__FUNCTION__+": internal error!"))
         return;
      //--- check
      if((circular && CFtBase::FtBaseIsSmooth(m)) && m%2==0)
        {
         //--- special code for circular convolution with smooth even M
         ArrayResizeAL(buf,m);
         for(i_=0;i_<=m-1;i_++)
            buf[i_]=a[i_];
         //--- allocation
         ArrayResizeAL(buf2,m);
         for(i_=0;i_<=n-1;i_++)
            buf2[i_]=b[i_];
         for(i=n;i<=m-1;i++)
            buf2[i]=0;
         //--- allocation
         ArrayResizeAL(buf3,m);
         //--- function call
         CFtBase::FtBaseGenerateComplexFFtPlan(m/2,plan);
         //--- function call
         CFastFourierTransform::FFTR1DInternalEven(buf,m,buf3,plan);
         //--- function call
         CFastFourierTransform::FFTR1DInternalEven(buf2,m,buf3,plan);
         //--- change values
         buf[0]=buf[0]*buf2[0];
         buf[1]=buf[1]*buf2[1];
         //--- calculation
         for(i=1;i<=m/2-1;i++)
           {
            ax=buf[2*i+0];
            ay=buf[2*i+1];
            bx=buf2[2*i+0];
            by=buf2[2*i+1];
            tx=ax*bx-ay*by;
            ty=ax*by+ay*bx;
            buf[2*i+0]=tx;
            buf[2*i+1]=ty;
           }
         //--- function call
         CFastFourierTransform::FFTR1DInvInternalEven(buf,m,buf3,plan);
         //--- allocation
         ArrayResizeAL(r,m);
         //--- copy
         for(i_=0;i_<=m-1;i_++)
            r[i_]=buf[i_];
        }
      else
        {
         //--- M is non-smooth or non-even,general code (circular/non-circular):
         //--- * first part is the same for circular and non-circular
         //---   convolutions. zero padding,FFTs,inverse FFTs
         //--- * second part differs:
         //---   * for non-circular convolution we just copy array
         //---   * for circular convolution we add array tail to its head
         p=CFtBase::FtBaseFindSmoothEven(m+n-1);
         //--- allocation
         ArrayResizeAL(buf,p);
         for(i_=0;i_<=m-1;i_++)
            buf[i_]=a[i_];
         for(i=m;i<=p-1;i++)
            buf[i]=0;
         //--- allocation
         ArrayResizeAL(buf2,p);
         for(i_=0;i_<=n-1;i_++)
            buf2[i_]=b[i_];
         for(i=n;i<=p-1;i++)
            buf2[i]=0;
         //--- allocation
         ArrayResizeAL(buf3,p);
         //--- function call
         CFtBase::FtBaseGenerateComplexFFtPlan(p/2,plan);
         //--- function call
         CFastFourierTransform::FFTR1DInternalEven(buf,p,buf3,plan);
         //--- function call
         CFastFourierTransform::FFTR1DInternalEven(buf2,p,buf3,plan);
         //--- change values
         buf[0]=buf[0]*buf2[0];
         buf[1]=buf[1]*buf2[1];
         //--- calculation
         for(i=1;i<=p/2-1;i++)
           {
            ax=buf[2*i+0];
            ay=buf[2*i+1];
            bx=buf2[2*i+0];
            by=buf2[2*i+1];
            tx=ax*bx-ay*by;
            ty=ax*by+ay*bx;
            buf[2*i+0]=tx;
            buf[2*i+1]=ty;
           }
         //--- function call
         CFastFourierTransform::FFTR1DInvInternalEven(buf,p,buf3,plan);
         //--- check
         if(circular)
           {
            //--- circular,add tail to head
            ArrayResizeAL(r,m);
            for(i_=0;i_<=m-1;i_++)
               r[i_]=buf[i_];
            //--- check
            if(n>=2)
              {
               i1_=m;
               for(i_=0;i_<=n-2;i_++)
                  r[i_]=r[i_]+buf[i_+i1_];
              }
           }
         else
           {
            //--- non-circular,just copy
            ArrayResizeAL(r,m+n-1);
            for(i_=0;i_<=m+n-2;i_++)
               r[i_]=buf[i_];
           }
        }
      //--- exit the function
      return;
     }
//--- overlap-add method
   if(alg==2)
     {
      //--- check
      if(!CAp::Assert((q+n-1)%2==0,__FUNCTION__+": internal error!"))
         return;
      //--- allocation
      ArrayResizeAL(buf,q+n-1);
      ArrayResizeAL(buf2,q+n-1);
      ArrayResizeAL(buf3,q+n-1);
      //--- function call
      CFtBase::FtBaseGenerateComplexFFtPlan((q+n-1)/2,plan);
      //--- prepare R
      if(circular)
        {
         //--- allocation
         ArrayResizeAL(r,m);
         for(i=0;i<=m-1;i++)
            r[i]=0;
        }
      else
        {
         //--- allocation
         ArrayResizeAL(r,m+n-1);
         for(i=0;i<=m+n-2;i++)
            r[i]=0;
        }
      //--- pre-calculated FFT(B)
      for(i_=0;i_<=n-1;i_++)
         buf2[i_]=b[i_];
      for(j=n;j<=q+n-2;j++)
         buf2[j]=0;
      //--- function call
      CFastFourierTransform::FFTR1DInternalEven(buf2,q+n-1,buf3,plan);
      //--- main overlap-add cycle
      i=0;
      //--- calculation
      while(i<=m-1)
        {
         p=MathMin(q,m-i);
         i1_=i;
         //--- copy
         for(i_=0;i_<=p-1;i_++)
            buf[i_]=a[i_+i1_];
         //--- make zero
         for(j=p;j<=q+n-2;j++)
            buf[j]=0;
         //--- function call
         CFastFourierTransform::FFTR1DInternalEven(buf,q+n-1,buf3,plan);
         //--- change values
         buf[0]=buf[0]*buf2[0];
         buf[1]=buf[1]*buf2[1];
         //--- calculation
         for(j=1;j<=(q+n-1)/2-1;j++)
           {
            ax=buf[2*j+0];
            ay=buf[2*j+1];
            bx=buf2[2*j+0];
            by=buf2[2*j+1];
            tx=ax*bx-ay*by;
            ty=ax*by+ay*bx;
            buf[2*j+0]=tx;
            buf[2*j+1]=ty;
           }
         //--- function call
         CFastFourierTransform::FFTR1DInvInternalEven(buf,q+n-1,buf3,plan);
         //--- check
         if(circular)
           {
            j1=MathMin(i+p+n-2,m-1)-i;
            j2=j1+1;
           }
         else
           {
            j1=p+n-2;
            j2=j1+1;
           }
         //--- change values
         i1_=-i;
         for(i_=i;i_<=i+j1;i_++)
            r[i_]=r[i_]+buf[i_+i1_];
         //--- check
         if(p+n-2>=j2)
           {
            i1_=j2;
            for(i_=0;i_<=p+n-2-j2;i_++)
               r[i_]=r[i_]+buf[i_+i1_];
           }
         i=i+p;
        }
      //--- exit the function
      return;
     }
  }
//+------------------------------------------------------------------+
//| Cross-correlation class                                          |
//+------------------------------------------------------------------+
class CCorr
  {
public:
   //--- constructor, destructor
                     CCorr(void);
                    ~CCorr(void);
   //--- methods
   static void       CorrC1D(al_complex &signal[],const int n,al_complex &pattern[],const int m,al_complex &r[]);
   static void       CorrC1DCircular(al_complex &signal[],const int m,al_complex &pattern[],const int n,al_complex &c[]);
   static void       CorrR1D(double &signal[],const int n,double &pattern[],const int m,double &r[]);
   static void       CorrR1DCircular(double &signal[],const int m,double &pattern[],const int n,double &c[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CCorr::CCorr(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CCorr::~CCorr(void)
  {

  }
//+------------------------------------------------------------------+
//| 1-dimensional complex cross-correlation.                         |
//| For given Pattern/Signal returns corr(Pattern,Signal)            |
//| (non-circular).                                                  |
//| Correlation is calculated using reduction to convolution.        |
//| Algorithm with max(N,N)*log(max(N,N)) complexity is used (see    |
//| ConvC1D() for more info about performance).                      |
//| IMPORTANT:                                                       |
//|     for historical reasons subroutine accepts its parameters in  |
//|     reversed order: CorrC1D(Signal, Pattern) = Pattern x Signal  |
//|     (using traditional definition of cross-correlation, denoting |
//|     cross-correlation as "x").                                   |
//| INPUT PARAMETERS                                                 |
//|     Signal  -   array[0..N-1] - complex function to be           |
//|                 transformed, signal containing pattern           |
//|     N       -   problem size                                     |
//|     Pattern -   array[0..M-1] - complex function to be           |
//|                 transformed, pattern to search withing signal    |
//|     M       -   problem size                                     |
//| OUTPUT PARAMETERS                                                |
//|     R       -   cross-correlation, array[0..N+M-2]:              |
//|                 * positive lags are stored in R[0..N-1],         |
//|                   R[i] = sum(conj(pattern[j])*signal[i+j]        |
//|                 * negative lags are stored in R[N..N+M-2],       |
//|                   R[N+M-1-i] = sum(conj(pattern[j])*signal[-i+j] |
//| NOTE:                                                            |
//|     It is assumed that pattern domain is [0..M-1]. If Pattern is |
//| non-zero on [-K..M-1], you can still use this subroutine, just   |
//| shift result by K.                                               |
//+------------------------------------------------------------------+
static void CCorr::CorrC1D(al_complex &signal[],const int n,al_complex &pattern[],
                           const int m,al_complex &r[])
  {
//--- create variables
   int i=0;
   int i_=0;
   int i1_=0;
//--- create arrays
   al_complex p[];
   al_complex b[];
//--- check
   if(!CAp::Assert(n>0 && m>0,__FUNCTION__+": incorrect N or M!"))
      return;
//--- allocation
   ArrayResizeAL(p,m);
   for(i=0;i<=m-1;i++)
      p[m-1-i]=CMath::Conj(pattern[i]);
//--- function call
   CConv::ConvC1D(p,m,signal,n,b);
//--- allocation
   ArrayResizeAL(r,m+n-1);
   i1_=m-1;
   for(i_=0;i_<=n-1;i_++)
      r[i_]=b[i_+i1_];
//--- check
   if(m+n-2>=n)
     {
      i1_=-n;
      for(i_=n;i_<=m+n-2;i_++)
         r[i_]=b[i_+i1_];
     }
  }
//+------------------------------------------------------------------+
//| 1-dimensional circular complex cross-correlation.                |
//| For given Pattern/Signal returns corr(Pattern,Signal) (circular).|
//| Algorithm has linearithmic complexity for any M/N.               |
//| IMPORTANT:                                                       |
//|     for historical reasons subroutine accepts its parameters in  |
//|     reversed order: CorrC1DCircular(Signal, Pattern) = Pattern x |
//|     Signal (using traditional definition of cross-correlation,   |
//|     denoting cross-correlation as "x").                          |
//| INPUT PARAMETERS                                                 |
//|     Signal  -   array[0..N-1] - complex function to be           |
//|                 transformed, periodic signal containing pattern  |
//|     N       -   problem size                                     |
//|     Pattern -   array[0..M-1] - complex function to be           |
//|                 transformed, non-periodic pattern to search      |
//|                 withing signal                                   |
//|     M       -   problem size                                     |
//| OUTPUT PARAMETERS                                                |
//|     R   -   convolution: A*B. array[0..M-1].                     |
//+------------------------------------------------------------------+
static void CCorr::CorrC1DCircular(al_complex &signal[],const int m,al_complex &pattern[],
                                   const int n,al_complex &c[])
  {
//--- create variables
   int i1=0;
   int i2=0;
   int i=0;
   int j2=0;
   int i_=0;
   int i1_=0;
//--- create arrays
   al_complex p[];
   al_complex b[];
//--- check
   if(!CAp::Assert(n>0 && m>0,__FUNCTION__+": incorrect N or M!"))
      return;
//--- normalize task: make M>=N,
//--- so A will be longer (at least - not shorter) that B.
   if(m<n)
     {
      //--- allocation
      ArrayResizeAL(b,m);
      //--- initialization
      for(i1=0;i1<=m-1;i1++)
         b[i1]=0;
      i1=0;
      //--- calculation
      while(i1<n)
        {
         //--- change values
         i2=MathMin(i1+m-1,n-1);
         j2=i2-i1;
         i1_=i1;
         for(i_=0;i_<=j2;i_++)
            b[i_]=b[i_]+pattern[i_+i1_];
         i1=i1+m;
        }
      //--- function call
      CorrC1DCircular(signal,m,b,m,c);
      //--- exit the function
      return;
     }
//--- Task is normalized
   ArrayResizeAL(p,n);
   for(i=0;i<=n-1;i++)
      p[n-1-i]=CMath::Conj(pattern[i]);
//--- function call
   CConv::ConvC1DCircular(signal,m,p,n,b);
//--- allocation
   ArrayResizeAL(c,m);
   i1_=n-1;
//--- copy
   for(i_=0;i_<=m-n;i_++)
      c[i_]=b[i_+i1_];
//--- check
   if(m-n+1<=m-1)
     {
      i1_=-(m-n+1);
      for(i_=m-n+1;i_<=m-1;i_++)
         c[i_]=b[i_+i1_];
     }
  }
//+------------------------------------------------------------------+
//| 1-dimensional real cross-correlation.                            |
//| For given Pattern/Signal returns corr(Pattern,Signal)            |
//| (non-circular).                                                  |
//| Correlation is calculated using reduction to convolution.        |
//| Algorithm with max(N,N)*log(max(N,N)) complexity is used (see    |
//| ConvC1D() for more info about performance).                      |
//| IMPORTANT:                                                       |
//|     for  historical reasons subroutine accepts its parameters in |
//|     reversed order: CorrR1D(Signal, Pattern) = Pattern x Signal  |
//|     (using  traditional definition of cross-correlation, denoting|
//|     cross-correlation as "x").                                   |
//| INPUT PARAMETERS                                                 |
//|     Signal  -   array[0..N-1] - real function to be transformed, |
//|                 signal containing pattern                        |
//|     N       -   problem size                                     |
//|     Pattern -   array[0..M-1] - real function to be transformed, |
//|                 pattern to search withing signal                 |
//|     M       -   problem size                                     |
//| OUTPUT PARAMETERS                                                |
//|     R       -   cross-correlation, array[0..N+M-2]:              |
//|                 * positive lags are stored in R[0..N-1],         |
//|                   R[i] = sum(pattern[j]*signal[i+j]              |
//|                 * negative lags are stored in R[N..N+M-2],       |
//|                   R[N+M-1-i] = sum(pattern[j]*signal[-i+j]       |
//| NOTE:                                                            |
//|     It is assumed that pattern domain is [0..M-1]. If Pattern is |
//| non-zero on [-K..M-1],  you can still use this subroutine, just  |
//| shift result by K.                                               |
//+------------------------------------------------------------------+
static void CCorr::CorrR1D(double &signal[],const int n,double &pattern[],
                           const int m,double &r[])
  {
//--- create variables
   int i=0;
   int i_=0;
   int i1_=0;
//--- create arrays
   double p[];
   double b[];
//--- check
   if(!CAp::Assert(n>0 && m>0,__FUNCTION__+": incorrect N or M!"))
      return;
//--- allocation
   ArrayResizeAL(p,m);
//--- copy
   for(i=0;i<=m-1;i++)
      p[m-1-i]=pattern[i];
//--- function call
   CConv::ConvR1D(p,m,signal,n,b);
//--- allocation
   ArrayResizeAL(r,m+n-1);
   i1_=m-1;
//--- copy
   for(i_=0;i_<=n-1;i_++)
      r[i_]=b[i_+i1_];
//--- check
   if(m+n-2>=n)
     {
      i1_=-n;
      for(i_=n;i_<=m+n-2;i_++)
         r[i_]=b[i_+i1_];
     }
  }
//+------------------------------------------------------------------+
//| 1-dimensional circular real cross-correlation.                   |
//| For given Pattern/Signal returns corr(Pattern,Signal) (circular).|
//| Algorithm has linearithmic complexity for any M/N.               |
//| IMPORTANT:                                                       |
//|     for  historical reasons subroutine accepts its parameters in |
//|     reversed order: CorrR1DCircular(Signal, Pattern) = Pattern x |
//|     Signal (using traditional definition of cross-correlation,   |
//|     denoting cross-correlation as "x").                          |
//| INPUT PARAMETERS                                                 |
//|     Signal  -   array[0..N-1] - real function to be transformed, |
//|                 periodic signal containing pattern               |
//|     N       -   problem size                                     |
//|     Pattern -   array[0..M-1] - real function to be transformed, |
//|                 non-periodic pattern to search withing signal    |
//|     M       -   problem size                                     |
//| OUTPUT PARAMETERS                                                |
//|     R   -   convolution: A*B. array[0..M-1].                     |
//+------------------------------------------------------------------+
static void CCorr::CorrR1DCircular(double &signal[],const int m,double &pattern[],
                                   const int n,double &c[])
  {
//--- create variables
   int i1=0;
   int i2=0;
   int i=0;
   int j2=0;
   int i_=0;
   int i1_=0;
//--- create arrays
   double p[];
   double b[];
//--- check
   if(!CAp::Assert(n>0 && m>0,__FUNCTION__+": incorrect N or M!"))
      return;
//--- normalize task: make M>=N,
//--- so A will be longer (at least - not shorter) that B.
   if(m<n)
     {
      //--- allocation
      ArrayResizeAL(b,m);
      //--- initialization
      for(i1=0;i1<=m-1;i1++)
         b[i1]=0;
      i1=0;
      //--- calculation
      while(i1<n)
        {
         //--- change values
         i2=MathMin(i1+m-1,n-1);
         j2=i2-i1;
         i1_=i1;
         for(i_=0;i_<=j2;i_++)
            b[i_]=b[i_]+pattern[i_+i1_];
         i1=i1+m;
        }
      //--- function call
      CorrR1DCircular(signal,m,b,m,c);
      //--- exit the function
      return;
     }
//--- Task is normalized
   ArrayResizeAL(p,n);
   for(i=0;i<=n-1;i++)
      p[n-1-i]=pattern[i];
//--- function call
   CConv::ConvR1DCircular(signal,m,p,n,b);
//--- allocation
   ArrayResizeAL(c,m);
   i1_=n-1;
//--- copy
   for(i_=0;i_<=m-n;i_++)
      c[i_]=b[i_+i1_];
//--- check
   if(m-n+1<=m-1)
     {
      i1_=-(m-n+1);
      for(i_=m-n+1;i_<=m-1;i_++)
         c[i_]=b[i_+i1_];
     }
  }
//+------------------------------------------------------------------+
//| Fast Hartley Transform                                           |
//+------------------------------------------------------------------+
class CFastHartleyTransform
  {
public:
   //--- constructor, destructor
                     CFastHartleyTransform(void);
                    ~CFastHartleyTransform(void);
   //--- methods
   static void       FHTR1D(double &a[],const int n);
   static void       FHTR1DInv(double &a[],const int n);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CFastHartleyTransform::CFastHartleyTransform(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CFastHartleyTransform::~CFastHartleyTransform(void)
  {

  }
//+------------------------------------------------------------------+
//| 1-dimensional Fast Hartley Transform.                            |
//| Algorithm has O(N*logN) complexity for any N (composite or prime)|
//| INPUT PARAMETERS                                                 |
//|     A  -  array[0..N-1] - real function to be transformed        |
//|     N  -  problem size                                           |
//| OUTPUT PARAMETERS                                                |
//|     A  -  FHT of a input array, array[0..N-1],                   |
//|           A_out[k]=sum(A_in[j]*(cos(2*pi*j*k/N)+sin(2*pi*j*k/N)),|
//|           j=0..N-1)                                              |
//+------------------------------------------------------------------+
static void CFastHartleyTransform::FHTR1D(double &a[],const int n)
  {
//--- create a variable
   int i=0;
//--- create array
   al_complex fa[];
//--- object of class
   CFtPlan plan;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- Special case: N=1,FHT is just identity transform.
//--- After this block we assume that N is strictly greater than 1.
   if(n==1)
      return;
//--- Reduce FHt to real FFT
   CFastFourierTransform::FFTR1D(a,n,fa);
//--- change values
   for(i=0;i<=n-1;i++)
      a[i]=fa[i].re-fa[i].im;
  }
//+------------------------------------------------------------------+
//| 1-dimensional inverse FHT.                                       |
//| Algorithm has O(N*logN) complexity for any N (composite or prime)|
//| INPUT PARAMETERS                                                 |
//|     A   -   array[0..N-1] - complex array to be transformed      |
//|     N   -   problem size                                         |
//| OUTPUT PARAMETERS                                                |
//|     A   -   inverse FHT of a input array, array[0..N-1]          |
//+------------------------------------------------------------------+
static void CFastHartleyTransform::FHTR1DInv(double &a[],const int n)
  {
//--- create a variable
   int i=0;
//--- check
   if(!CAp::Assert(n>0,__FUNCTION__+": incorrect N!"))
      return;
//--- Special case: N=1,iFHT is just identity transform.
//--- After this block we assume that N is strictly greater than 1.
   if(n==1)
      return;
//--- Inverse FHT can be expressed in terms of the FHT as
//---     invfht(x)=fht(x)/N
   FHTR1D(a,n);
//--- change values
   for(i=0;i<=n-1;i++)
      a[i]=a[i]/n;
  }
//+------------------------------------------------------------------+
