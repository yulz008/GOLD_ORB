//+------------------------------------------------------------------+
//|                                                 dataanalysis.mqh |
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
#include "ap.mqh"
#include "optimization.mqh"
#include "statistics.mqh"
#include "solvers.mqh"
//+------------------------------------------------------------------+
//| Auxiliary class for CBdSS                                        |
//+------------------------------------------------------------------+
class CCVReport
  {
public:
   double            m_relclserror;
   double            m_avgce;
   double            m_rmserror;
   double            m_avgerror;
   double            m_avgrelerror;

                     CCVReport(void);
                    ~CCVReport(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CCVReport::CCVReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CCVReport::~CCVReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Data analysis                                                    |
//+------------------------------------------------------------------+
class CBdSS
  {
private:
   //--- private methods
   static double     XLnY(const double x,const double y);
   static double     GetCV(int &cnt[],const int nc);
   static void       TieAddC(int &c[],int &ties[],const int ntie,const int nc,int &cnt[]);
   static void       TieSubC(int &c[],int &ties[],const int ntie,const int nc,int &cnt[]);
public:
   //--- constructor, destructor
                     CBdSS(void);
                    ~CBdSS(void);
   //--- public methods
   static void       DSErrAllocate(const int nclasses,double &buf[]);
   static void       DSErrAccumulate(double &buf[],double &y[],double &desiredy[]);
   static void       DSErrFinish(double &buf[]);
   static void       DSNormalize(CMatrixDouble &xy,const int npoints,const int nvars,int &info,double &means[],double &sigmas[]);
   static void       DSNormalizeC(CMatrixDouble &xy,const int npoints,const int nvars,int &info,double &means[],double &sigmas[]);
   static double     DSGetMeanMindIstance(CMatrixDouble &xy,const int npoints,const int nvars);
   static void       DSTie(double &a[],const int n,int &ties[],int &tiecount,int &p1[],int &p2[]);
   static void       DSTieFastI(double &a[],int &b[],const int n,int &ties[],int &tiecount,double &bufr[],int &bufi[]);
   static void       DSOptimalSplit2(double &ca[],int &cc[],const int n,int &info,double &threshold,double &pal,double &pbl,double &par,double &pbr,double &cve);
   static void       DSOptimalSplit2Fast(double &a[],int &c[],int &tiesbuf[],int &cntbuf[],double &bufr[],int &bufi[],const int n,const int nc,double alpha,int &info,double &threshold,double &rms,double &cvrms);
   static void       DSSplitK(double &ca[],int &cc[],const int n,const int nc,int kmax,int &info,double &thresholds[],int &ni,double &cve);
   static void       DSOptimalSplitK(double &ca[],int &cc[],const int n,const int nc,int kmax,int &info,double &thresholds[],int &ni,double &cve);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CBdSS::CBdSS(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CBdSS::~CBdSS(void)
  {

  }
//+------------------------------------------------------------------+
//| This set of routines (DSErrAllocate, DSErrAccumulate,            |
//| DSErrFinish) calculates different error functions (classification| 
//| error, cross-entropy, rms, avg, avg.rel errors).                 |
//| 1. DSErrAllocate prepares buffer.                                |
//| 2. DSErrAccumulate accumulates individual errors:                |
//|     * Y contains predicted output (posterior probabilities for   |
//|       classification)                                            |
//|     * DesiredY contains desired output (class number for         |
//|       classification)                                            |
//| 3. DSErrFinish outputs results:                                  |
//|    * Buf[0] contains relative classification error (zero for     |
//|      regression tasks)                                           |
//|    * Buf[1] contains avg. cross-entropy (zero for regression     |
//|      tasks)                                                      |
//|    * Buf[2] contains rms error (regression, classification)      |
//|    * Buf[3] contains average error (regression, classification)  |
//|    * Buf[4] contains average relative error (regression,         |
//|      classification)                                             |
//| NOTES(1):                                                        |
//|     "NClasses>0" means that we have classification task.         |
//|     "NClasses<0" means regression task with -NClasses real       |
//|     outputs.                                                     |
//| NOTES(2):                                                        |
//|     rms. avg, avg.rel errors for classification tasks are        |
//|     interpreted as errors in posterior probabilities with        |
//|     respect to probabilities given by training/test set.         |
//+------------------------------------------------------------------+
static void CBdSS::DSErrAllocate(const int nclasses,double &buf[])
  {
//--- allocation
   ArrayResizeAL(buf,8);
//--- initialization
   buf[0]=0;
   buf[1]=0;
   buf[2]=0;
   buf[3]=0;
   buf[4]=0;
   buf[5]=nclasses;
   buf[6]=0;
   buf[7]=0;
  }
//+------------------------------------------------------------------+
//| See DSErrAllocate for comments on this routine.                  |
//+------------------------------------------------------------------+
static void CBdSS::DSErrAccumulate(double &buf[],double &y[],double &desiredy[])
  {
//--- create variables
   int    nclasses=0;
   int    nout=0;
   int    offs=0;
   int    mmax=0;
   int    rmax=0;
   int    j=0;
   double v=0;
   double ev=0;
//--- initialization
   offs=5;
   nclasses=(int)MathRound(buf[offs]);
//--- check
   if(nclasses>0)
     {
      //--- Classification
      rmax=(int)MathRound(desiredy[0]);
      mmax=0;
      //--- initialization
      for(j=1;j<=nclasses-1;j++)
        {
         //--- check
         if(y[j]>y[mmax])
            mmax=j;
        }
      //--- check
      if(mmax!=rmax)
         buf[0]=buf[0]+1;
      //--- check
      if(y[rmax]>0.0)
         buf[1]=buf[1]-MathLog(y[rmax]);
      else
         buf[1]=buf[1]+MathLog(CMath::m_maxrealnumber);
      //--- calculation
      for(j=0;j<=nclasses-1;j++)
        {
         v=y[j];
         //--- check
         if(j==rmax)
            ev=1;
         else
            ev=0;
         //--- change values
         buf[2]=buf[2]+CMath::Sqr(v-ev);
         buf[3]=buf[3]+MathAbs(v-ev);
         //--- check
         if(ev!=0.0)
           {
            buf[4]=buf[4]+MathAbs((v-ev)/ev);
            buf[offs+2]=buf[offs+2]+1;
           }
        }
      //--- change value
      buf[offs+1]=buf[offs+1]+1;
     }
   else
     {
      //--- Regression
      nout=-nclasses;
      rmax=0;
      //--- initialization
      for(j=1;j<=nout-1;j++)
        {
         //--- check
         if(desiredy[j]>desiredy[rmax])
            rmax=j;
        }
      //--- initialization
      mmax=0;
      for(j=1;j<=nout-1;j++)
        {
         //--- check
         if(y[j]>y[mmax])
            mmax=j;
        }
      //--- check
      if(mmax!=rmax)
         buf[0]=buf[0]+1;
      //--- calculation
      for(j=0;j<=nout-1;j++)
        {
         //--- change values
         v=y[j];
         ev=desiredy[j];
         buf[2]=buf[2]+CMath::Sqr(v-ev);
         buf[3]=buf[3]+MathAbs(v-ev);
         //--- check
         if(ev!=0.0)
           {
            buf[4]=buf[4]+MathAbs((v-ev)/ev);
            buf[offs+2]=buf[offs+2]+1;
           }
        }
      //--- change value
      buf[offs+1]=buf[offs+1]+1;
     }
  }
//+------------------------------------------------------------------+
//| See DSErrAllocate for comments on this routine.                  |
//+------------------------------------------------------------------+
static void CBdSS::DSErrFinish(double &buf[])
  {
//--- create variables
   int nout=0;
   int offs=0;
//--- initialization
   offs=5;
   nout=(int)(MathAbs((int)MathRound(buf[offs])));
//--- check
   if(buf[offs+1]!=0.0)
     {
      //--- change values
      buf[0]=buf[0]/buf[offs+1];
      buf[1]=buf[1]/buf[offs+1];
      buf[2]=MathSqrt(buf[2]/(nout*buf[offs+1]));
      buf[3]=buf[3]/(nout*buf[offs+1]);
     }
//--- check
   if(buf[offs+2]!=0.0)
      buf[4]=buf[4]/buf[offs+2];
  }
//+------------------------------------------------------------------+
//| Normalize                                                        |
//+------------------------------------------------------------------+
static void CBdSS::DSNormalize(CMatrixDouble &xy,const int npoints,
                               const int nvars,int &info,double &means[],
                               double &sigmas[])
  {
//--- create variables
   int    i=0;
   int    j=0;
   double mean=0;
   double variance=0;
   double skewness=0;
   double kurtosis=0;
   int    i_=0;
//--- create array
   double tmp[];
//--- initialization
   info=0;
//--- Test parameters
   if(npoints<=0 || nvars<1)
     {
      info=-1;
      return;
     }
//--- change value
   info=1;
//--- Standartization
   ArrayResizeAL(means,nvars);
   ArrayResizeAL(sigmas,nvars);
   ArrayResizeAL(tmp,npoints);
//--- calculation
   for(j=0;j<=nvars-1;j++)
     {
      //--- copy
      for(i_=0;i_<=npoints-1;i_++)
         tmp[i_]=xy[i_][j];
      //--- function call
      CBaseStat::SampleMoments(tmp,npoints,mean,variance,skewness,kurtosis);
      //--- change values
      means[j]=mean;
      sigmas[j]=MathSqrt(variance);
      //--- check
      if(sigmas[j]==0.0)
         sigmas[j]=1;
      //--- change values
      for(i=0;i<=npoints-1;i++)
         xy[i].Set(j,(xy[i][j]-means[j])/sigmas[j]);
     }
  }
//+------------------------------------------------------------------+
//| Normalize                                                        |
//+------------------------------------------------------------------+
static void CBdSS::DSNormalizeC(CMatrixDouble &xy,const int npoints,
                                const int nvars,int &info,double &means[],
                                double &sigmas[])
  {
//--- create variables
   int    j=0;
   double mean=0;
   double variance=0;
   double skewness=0;
   double kurtosis=0;
   int    i_=0;
//--- create array
   double tmp[];
//--- initialization
   info=0;
//--- Test parameters
   if(npoints<=0 || nvars<1)
     {
      info=-1;
      return;
     }
//--- change value
   info=1;
//--- Standartization
   ArrayResizeAL(means,nvars);
   ArrayResizeAL(sigmas,nvars);
   ArrayResizeAL(tmp,npoints);
   for(j=0;j<=nvars-1;j++)
     {
      //--- copy
      for(i_=0;i_<=npoints-1;i_++)
         tmp[i_]=xy[i_][j];
      //--- function call
      CBaseStat::SampleMoments(tmp,npoints,mean,variance,skewness,kurtosis);
      //--- change values
      means[j]=mean;
      sigmas[j]=MathSqrt(variance);
      //--- check
      if(sigmas[j]==0.0)
         sigmas[j]=1;
     }
  }
//+------------------------------------------------------------------+
//| Method                                                           |
//+------------------------------------------------------------------+
static double CBdSS::DSGetMeanMindIstance(CMatrixDouble &xy,const int npoints,
                                          const int nvars)
  {
//--- create variables
   double result=0;
   int    i=0;
   int    j=0;
   double v=0;
   int    i_=0;
//--- creating arrays
   double tmp[];
   double tmp2[];
//--- Test parameters
   if(npoints<=0 || nvars<1)
      return(0);
//--- Process
   ArrayResizeAL(tmp,npoints);
   for(i=0;i<=npoints-1;i++)
      tmp[i]=CMath::m_maxrealnumber;
//--- allocation
   ArrayResizeAL(tmp2,nvars);
   for(i=0;i<=npoints-1;i++)
     {
      for(j=i+1;j<=npoints-1;j++)
        {
         //--- calculation
         for(i_=0;i_<=nvars-1;i_++)
            tmp2[i_]=xy[i][i_];
         for(i_=0;i_<=nvars-1;i_++)
            tmp2[i_]=tmp2[i_]-xy[j][i_];
         v=0.0;
         for(i_=0;i_<=nvars-1;i_++)
            v+=tmp2[i_]*tmp2[i_];
         //--- change values
         v=MathSqrt(v);
         tmp[i]=MathMin(tmp[i],v);
         tmp[j]=MathMin(tmp[j],v);
        }
     }
//--- get result
   result=0;
   for(i=0;i<=npoints-1;i++)
      result=result+tmp[i]/npoints;
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Method                                                           |
//+------------------------------------------------------------------+
static void CBdSS::DSTie(double &a[],const int n,int &ties[],int &tiecount,
                         int &p1[],int &p2[])
  {
//--- create variables
   int i=0;
   int k=0;
//--- create array
   int tmp[];
//--- initialization
   tiecount=0;
//--- Special case
   if(n<=0)
     {
      tiecount=0;
      return;
     }
//--- Sort A
   CTSort::TagSort(a,n,p1,p2);
//--- Process ties
   tiecount=1;
   for(i=1;i<=n-1;i++)
     {
      //--- check
      if(a[i]!=a[i-1])
         tiecount=tiecount+1;
     }
//--- allocation
   ArrayResizeAL(ties,tiecount+1);
//--- change values
   ties[0]=0;
   k=1;
//--- calculation
   for(i=1;i<=n-1;i++)
     {
      //--- check
      if(a[i]!=a[i-1])
        {
         ties[k]=i;
         k=k+1;
        }
     }
//--- change value
   ties[tiecount]=n;
  }
//+------------------------------------------------------------------+
//| Method                                                           |
//+------------------------------------------------------------------+
static void CBdSS::DSTieFastI(double &a[],int &b[],const int n,int &ties[],
                              int &tiecount,double &bufr[],int &bufi[])
  {
//--- create variables
   int i=0;
   int k=0;
//--- create array
   int tmp[];
//--- initialization
   tiecount=0;
//--- Special case
   if(n<=0)
     {
      tiecount=0;
      return;
     }
//--- Sort A
   CTSort::TagSortFastI(a,b,bufr,bufi,n);
//--- Process ties
   ties[0]=0;
   k=1;
//--- calculation
   for(i=1;i<=n-1;i++)
     {
      //--- check
      if(a[i]!=a[i-1])
        {
         ties[k]=i;
         k=k+1;
        }
     }
//--- change values
   ties[k]=n;
   tiecount=k;
  }
//+------------------------------------------------------------------+
//| Optimal binary classification                                    |
//| Algorithms finds optimal (=with minimal cross-entropy) binary    |
//| partition.                                                       |
//| Internal subroutine.                                             |
//| INPUT PARAMETERS:                                                |
//|     A       -   array[0..N-1], variable                          |
//|     C       -   array[0..N-1], class numbers (0 or 1).           |
//|     N       -   array size                                       |
//| OUTPUT PARAMETERS:                                               |
//|     Info    -   completetion code:                               |
//|                 * -3, all values of A[] are same (partition is   |
//|                   impossible)                                    |
//|                 * -2, one of C[] is incorrect (<0, >1)           |
//|                 * -1, incorrect pararemets were passed (N<=0).   |
//|                 *  1, OK                                         |
//|     Threshold-  partiton boundary. Left part contains values     |
//|                 which are strictly less than Threshold. Right    |
//|                 part contains values which are greater than or   |
//|                 equal to Threshold.                              |
//|     PAL, PBL-   probabilities P(0|v<Threshold) and               |
//|                 P(1|v<Threshold)                                 |
//|     PAR, PBR-   probabilities P(0|v>=Threshold) and              |
//|                 P(1|v>=Threshold)                                |
//|     CVE     -   cross-validation estimate of cross-entropy       |
//+------------------------------------------------------------------+
static void CBdSS::DSOptimalSplit2(double &ca[],int &cc[],const int n,
                                   int &info,double &threshold,double &pal,
                                   double &pbl,double &par,double &pbr,
                                   double &cve)
  {
//--- create variables
   int    i=0;
   int    t=0;
   double s=0;
   int    tiecount=0;
   int    k=0;
   int    koptimal=0;
   double pak=0;
   double pbk=0;
   double cvoptimal=0;
   double cv=0;
//--- creating arrays
   int    ties[];
   int    p1[];
   int    p2[];
   double a[];
   int    c[];
//--- copy
   ArrayCopy(a,ca);
   ArrayCopy(c,cc);
//--- initialization
   info=0;
   threshold=0;
   pal=0;
   pbl=0;
   par=0;
   pbr=0;
   cve=0;
//--- Test for errors in inputs
   if(n<=0)
     {
      info=-1;
      return;
     }
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(c[i]!=0 && c[i]!=1)
        {
         info=-2;
         return;
        }
     }
//--- change value
   info=1;
//--- Tie
   DSTie(a,n,ties,tiecount,p1,p2);
//--- swap
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(p2[i]!=i)
        {
         t=c[i];
         c[i]=c[p2[i]];
         c[p2[i]]=t;
        }
     }
//--- Special case: number of ties is 1.
//--- NOTE: we assume that P[i][j] equals to 0 or 1,
//---       intermediate values are not allowed.
   if(tiecount==1)
     {
      info=-3;
      return;
     }
//--- General case,number of ties > 1
//--- NOTE: we assume that P[i][j] equals to 0 or 1,
//---       intermediate values are not allowed.
   pal=0;
   pbl=0;
   par=0;
   pbr=0;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(c[i]==0)
         par=par+1;
      //--- check
      if(c[i]==1)
         pbr=pbr+1;
     }
//--- change values
   koptimal=-1;
   cvoptimal=CMath::m_maxrealnumber;
   for(k=0;k<=tiecount-2;k++)
     {
      //--- first,obtain information about K-th tie which is
      //--- moved from R-part to L-part
      pak=0;
      pbk=0;
      for(i=ties[k];i<=ties[k+1]-1;i++)
        {
         //--- check
         if(c[i]==0)
            pak=pak+1;
         //--- check
         if(c[i]==1)
            pbk=pbk+1;
        }
      //--- Calculate cross-validation CE
      cv=0;
      cv=cv-XLnY(pal+pak,(pal+pak)/(pal+pak+pbl+pbk+1));
      cv=cv-XLnY(pbl+pbk,(pbl+pbk)/(pal+pak+1+pbl+pbk));
      cv=cv-XLnY(par-pak,(par-pak)/(par-pak+pbr-pbk+1));
      cv=cv-XLnY(pbr-pbk,(pbr-pbk)/(par-pak+1+pbr-pbk));
      //--- Compare with best
      if(cv<cvoptimal)
        {
         cvoptimal=cv;
         koptimal=k;
        }
      //--- update
      pal=pal+pak;
      pbl=pbl+pbk;
      par=par-pak;
      pbr=pbr-pbk;
     }
//--- change values
   cve=cvoptimal;
   threshold=0.5*(a[ties[koptimal]]+a[ties[koptimal+1]]);
   pal=0;
   pbl=0;
   par=0;
   pbr=0;
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(a[i]<threshold)
        {
         //--- check
         if(c[i]==0)
            pal=pal+1;
         else
            pbl=pbl+1;
        }
      else
        {
         //--- check
         if(c[i]==0)
            par=par+1;
         else
            pbr=pbr+1;
        }
     }
//--- change values
   s=pal+pbl;
   pal=pal/s;
   pbl=pbl/s;
   s=par+pbr;
   par=par/s;
   pbr=pbr/s;
  }
//+------------------------------------------------------------------+
//| Optimal partition, internal subroutine. Fast version.            |
//| Accepts:                                                         |
//|     A       array[0..N-1]       array of attributes array[0..N-1]|
//|     C       array[0..N-1]       array of class labels            |
//|     TiesBuf array[0..N]         temporaries (ties)               |
//|     CntBuf  array[0..2*NC-1]    temporaries (counts)             |
//|     Alpha                       centering factor (0<=alpha<=1,   |
//|                                 recommended value - 0.05)        |
//|     BufR    array[0..N-1]       temporaries                      |
//|     BufI    array[0..N-1]       temporaries                      |
//| Output:                                                          |
//|     Info    error code (">0"=OK, "<0"=bad)                       |
//|     RMS     training set RMS error                               |
//|     CVRMS   leave-one-out RMS error                              |
//| Note:                                                            |
//|     content of all arrays is changed by subroutine;              |
//|     it doesn't allocate temporaries.                             |
//+------------------------------------------------------------------+
static void CBdSS::DSOptimalSplit2Fast(double &a[],int &c[],int &tiesbuf[],
                                       int &cntbuf[],double &bufr[],int &bufi[],
                                       const int n,const int nc,double alpha,
                                       int &info,double &threshold,
                                       double &rms,double &cvrms)
  {
//--- create variables
   int    i=0;
   int    k=0;
   int    cl=0;
   int    tiecount=0;
   double cbest=0;
   double cc=0;
   int    koptimal=0;
   int    sl=0;
   int    sr=0;
   double v=0;
   double w=0;
   double x=0;
//--- initialization
   info=0;
   threshold=0;
   rms=0;
   cvrms=0;
//--- Test for errors in inputs
   if(n<=0 || nc<2)
     {
      info=-1;
      return;
     }
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(c[i]<0 || c[i]>=nc)
        {
         info=-2;
         return;
        }
     }
//--- change value
   info=1;
//--- Tie
   DSTieFastI(a,c,n,tiesbuf,tiecount,bufr,bufi);
//--- Special case: number of ties is 1.
   if(tiecount==1)
     {
      info=-3;
      return;
     }
//--- General case,number of ties > 1
   for(i=0;i<=2*nc-1;i++)
      cntbuf[i]=0;
   for(i=0;i<=n-1;i++)
      cntbuf[nc+c[i]]=cntbuf[nc+c[i]]+1;
//--- change values
   koptimal=-1;
   threshold=a[n-1];
   cbest=CMath::m_maxrealnumber;
   sl=0;
   sr=n;
//--- calculation
   for(k=0;k<=tiecount-2;k++)
     {
      //--- first,move Kth tie from right to left
      for(i=tiesbuf[k];i<=tiesbuf[k+1]-1;i++)
        {
         cl=c[i];
         cntbuf[cl]=cntbuf[cl]+1;
         cntbuf[nc+cl]=cntbuf[nc+cl]-1;
        }
      sl=sl+(tiesbuf[k+1]-tiesbuf[k]);
      sr=sr-(tiesbuf[k+1]-tiesbuf[k]);
      //--- Calculate RMS error
      v=0;
      for(i=0;i<=nc-1;i++)
        {
         w=cntbuf[i];
         v=v+w*CMath::Sqr(w/sl-1);
         v=v+(sl-w)*CMath::Sqr(w/sl);
         w=cntbuf[nc+i];
         v=v+w*CMath::Sqr(w/sr-1);
         v=v+(sr-w)*CMath::Sqr(w/sr);
        }
      //--- change value
      v=MathSqrt(v/(nc*n));
      //--- Compare with best
      x=(double)(2*sl)/(double)(sl+sr)-1;
      cc=v*(1-alpha+alpha*CMath::Sqr(x));
      //--- check
      if(cc<cbest)
        {
         //--- store split
         rms=v;
         koptimal=k;
         cbest=cc;
         //--- calculate CVRMS error
         cvrms=0;
         for(i=0;i<=nc-1;i++)
           {
            //--- check
            if(sl>1)
              {
               w=cntbuf[i];
               cvrms=cvrms+w*CMath::Sqr((w-1)/(sl-1)-1);
               cvrms=cvrms+(sl-w)*CMath::Sqr(w/(sl-1));
              }
            else
              {
               w=cntbuf[i];
               cvrms=cvrms+w*CMath::Sqr(1.0/(double)nc-1);
               cvrms=cvrms+(sl-w)*CMath::Sqr(1.0/(double)nc);
              }
            //--- check
            if(sr>1)
              {
               w=cntbuf[nc+i];
               cvrms=cvrms+w*CMath::Sqr((w-1)/(sr-1)-1);
               cvrms=cvrms+(sr-w)*CMath::Sqr(w/(sr-1));
              }
            else
              {
               w=cntbuf[nc+i];
               cvrms=cvrms+w*CMath::Sqr(1.0/(double)nc-1);
               cvrms=cvrms+(sr-w)*CMath::Sqr(1.0/(double)nc);
              }
           }
         //--- change value
         cvrms=MathSqrt(cvrms/(nc*n));
        }
     }
//--- Calculate threshold.
//--- Code is a bit complicated because there can be such
//--- numbers that 0.5(A+B) equals to A or B (if A-B=epsilon)
   threshold=0.5*(a[tiesbuf[koptimal]]+a[tiesbuf[koptimal+1]]);
//--- check
   if(threshold<=a[tiesbuf[koptimal]])
      threshold=a[tiesbuf[koptimal+1]];
  }
//+------------------------------------------------------------------+
//| Automatic non-optimal discretization, internal subroutine.       |
//+------------------------------------------------------------------+
static void CBdSS::DSSplitK(double &ca[],int &cc[],const int n,const int nc,
                            int kmax,int &info,double &thresholds[],int &ni,
                            double &cve)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    j1=0;
   int    k=0;
   int    tiecount=0;
   double v2=0;
   int    bestk=0;
   double bestcve=0;
   double curcve=0;
//--- creating arrays
   int    ties[];
   int    p1[];
   int    p2[];
   int    cnt[];
   int    bestsizes[];
   int    cursizes[];
   double a[];
   int    c[];
//--- copy
   ArrayCopy(a,ca);
   ArrayCopy(c,cc);
//--- initialization
   info=0;
   ni=0;
   cve=0;
//--- Test for errors in inputs
   if((n<=0 || nc<2) || kmax<2)
     {
      info=-1;
      return;
     }
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(c[i]<0 || c[i]>=nc)
        {
         info=-2;
         return;
        }
     }
//--- change value
   info=1;
//--- Tie
   DSTie(a,n,ties,tiecount,p1,p2);
//--- swap
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(p2[i]!=i)
        {
         k=c[i];
         c[i]=c[p2[i]];
         c[p2[i]]=k;
        }
     }
//--- Special cases
   if(tiecount==1)
     {
      info=-3;
      return;
     }
//--- General case:
//--- 0. allocate arrays
   kmax=MathMin(kmax,tiecount);
//--- allocation
   ArrayResizeAL(bestsizes,kmax);
   ArrayResizeAL(cursizes,kmax);
   ArrayResizeAL(cnt,nc);
//--- General case:
//--- 1. prepare "weak" solution (two subintervals,divided at median)
   v2=CMath::m_maxrealnumber;
   j=-1;
   for(i=1;i<=tiecount-1;i++)
     {
      //--- check
      if(MathAbs(ties[i]-0.5*(n-1))<v2)
        {
         v2=MathAbs(ties[i]-0.5*n);
         j=i;
        }
     }
//--- check
   if(!CAp::Assert(j>0,__FUNCTION__+": internal error #1!"))
      return;
//--- change values
   bestk=2;
   bestsizes[0]=ties[j];
   bestsizes[1]=n-j;
   bestcve=0;
//--- calculation
   for(i=0;i<=nc-1;i++)
      cnt[i]=0;
   for(i=0;i<=j-1;i++)
      TieAddC(c,ties,i,nc,cnt);
   bestcve=bestcve+GetCV(cnt,nc);
//--- calculation
   for(i=0;i<=nc-1;i++)
      cnt[i]=0;
   for(i=j;i<=tiecount-1;i++)
      TieAddC(c,ties,i,nc,cnt);
   bestcve=bestcve+GetCV(cnt,nc);
//--- General case:
//--- 2. Use greedy algorithm to find sub-optimal split in O(KMax*N) time
   for(k=2;k<=kmax;k++)
     {
      //--- Prepare greedy K-interval split
      for(i=0;i<=k-1;i++)
         cursizes[i]=0;
      //--- change values
      i=0;
      j=0;
      //--- cycle
      while(j<=tiecount-1 && i<=k-1)
        {
         //--- Rule: I-th bin is empty,fill it
         if(cursizes[i]==0)
           {
            cursizes[i]=ties[j+1]-ties[j];
            j=j+1;
            continue;
           }
         //--- Rule: (K-1-I) bins left,(K-1-I) ties left (1 tie per bin);next bin
         if(tiecount-j==k-1-i)
           {
            i=i+1;
            continue;
           }
         //--- Rule: last bin,always place in current
         if(i==k-1)
           {
            cursizes[i]=cursizes[i]+ties[j+1]-ties[j];
            j=j+1;
            continue;
           }
         //--- Place J-th tie in I-th bin,or leave for I+1-th bin.
         if(MathAbs(cursizes[i]+ties[j+1]-ties[j]-(double)n/(double)k)<MathAbs(cursizes[i]-(double)n/(double)k))
           {
            cursizes[i]=cursizes[i]+ties[j+1]-ties[j];
            j=j+1;
           }
         else
            i=i+1;
        }
      //--- check
      if(!CAp::Assert(cursizes[k-1]!=0 && j==tiecount,__FUNCTION__+": internal error #1"))
         return;
      //--- Calculate CVE
      curcve=0;
      j=0;
      for(i=0;i<=k-1;i++)
        {
         //--- calculation
         for(j1=0;j1<=nc-1;j1++)
            cnt[j1]=0;
         for(j1=j;j1<=j+cursizes[i]-1;j1++)
            cnt[c[j1]]=cnt[c[j1]]+1;
         curcve=curcve+GetCV(cnt,nc);
         j=j+cursizes[i];
        }
      //--- Choose best variant
      if(curcve<bestcve)
        {
         for(i=0;i<=k-1;i++)
            bestsizes[i]=cursizes[i];
         bestcve=curcve;
         bestk=k;
        }
     }
//--- Transform from sizes to thresholds
   cve=bestcve;
   ni=bestk;
//--- allocation
   ArrayResizeAL(thresholds,ni-1);
   j=bestsizes[0];
//--- calculation
   for(i=1;i<=bestk-1;i++)
     {
      thresholds[i-1]=0.5*(a[j-1]+a[j]);
      j=j+bestsizes[i];
     }
  }
//+------------------------------------------------------------------+
//| Automatic optimal discretization, internal subroutine.           |
//+------------------------------------------------------------------+
static void CBdSS::DSOptimalSplitK(double &ca[],int &cc[],const int n,
                                   const int nc,int kmax,int &info,
                                   double &thresholds[],int &ni,double &cve)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    s=0;
   int    jl=0;
   int    jr=0;
   double v2=0;
   int    tiecount=0;
   double cvtemp=0;
   int    k=0;
   int    koptimal=0;
   double cvoptimal=0;
//--- creating arrays
   int    ties[];
   int    p1[];
   int    p2[];
   int    cnt[];
   int    cnt2[];
   double a[];
   int    c[];
//--- create matrix
   CMatrixDouble cv;
   CMatrixInt    splits;
//--- copy
   ArrayCopy(a,ca);
   ArrayCopy(c,cc);
//--- initialization
   info=0;
   ni=0;
   cve=0;
//--- Test for errors in inputs
   if((n<=0 || nc<2) || kmax<2)
     {
      info=-1;
      return;
     }
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(c[i]<0 || c[i]>=nc)
        {
         info=-2;
         return;
        }
     }
//--- change value
   info=1;
//--- Tie
   DSTie(a,n,ties,tiecount,p1,p2);
//--- swap
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(p2[i]!=i)
        {
         k=c[i];
         c[i]=c[p2[i]];
         c[p2[i]]=k;
        }
     }
//--- Special cases
   if(tiecount==1)
     {
      info=-3;
      return;
     }
//--- General case
//--- Use dynamic programming to find best split in O(KMax*NC*TieCount^2) time
   kmax=MathMin(kmax,tiecount);
//--- allocation
   cv.Resize(kmax,tiecount);
   splits.Resize(kmax,tiecount);
   ArrayResizeAL(cnt,nc);
   ArrayResizeAL(cnt2,nc);
//--- calculation
   for(j=0;j<=nc-1;j++)
      cnt[j]=0;
   for(j=0;j<=tiecount-1;j++)
     {
      TieAddC(c,ties,j,nc,cnt);
      splits[0].Set(j,0);
      cv[0].Set(j,GetCV(cnt,nc));
     }
   for(k=1;k<=kmax-1;k++)
     {
      for(j=0;j<=nc-1;j++)
         cnt[j]=0;
      //--- Subtask size J in [K..TieCount-1]:
      //--- optimal K-splitting on ties from 0-th to J-th.
      for(j=k;j<=tiecount-1;j++)
        {
         //--- Update Cnt - let it contain classes of ties from K-th to J-th
         TieAddC(c,ties,j,nc,cnt);
         //--- Search for optimal split point S in [K..J]
         for(i=0;i<=nc-1;i++)
            cnt2[i]=cnt[i];
         cv[k].Set(j,cv[k-1][j-1]+GetCV(cnt2,nc));
         splits[k].Set(j,j);
         //--- calculation
         for(s=k+1;s<=j;s++)
           {
            //--- Update Cnt2 - let it contain classes of ties from S-th to J-th
            TieSubC(c,ties,s-1,nc,cnt2);
            //--- Calculate CVE
            cvtemp=cv[k-1][s-1]+GetCV(cnt2,nc);
            //--- check
            if(cvtemp<cv[k][j])
              {
               cv[k].Set(j,cvtemp);
               splits[k].Set(j,s);
              }
           }
        }
     }
//--- Choose best partition,output result
   koptimal=-1;
   cvoptimal=CMath::m_maxrealnumber;
   for(k=0;k<=kmax-1;k++)
     {
      //--- check
      if(cv[k][tiecount-1]<cvoptimal)
        {
         cvoptimal=cv[k][tiecount-1];
         koptimal=k;
        }
     }
//--- check
   if(!CAp::Assert(koptimal>=0,__FUNCTION__+": internal error #1!"))
      return;
//--- check
   if(koptimal==0)
     {
      //--- Special case: best partition is one big interval.
      //--- Even 2-partition is not better.
      //--- This is possible when dealing with "weak" predictor variables.
      //--- Make binary split as close to the median as possible.
      v2=CMath::m_maxrealnumber;
      j=-1;
      for(i=1;i<=tiecount-1;i++)
        {
         //--- check
         if(MathAbs(ties[i]-0.5*(n-1))<v2)
           {
            v2=MathAbs(ties[i]-0.5*(n-1));
            j=i;
           }
        }
      //--- check
      if(!CAp::Assert(j>0,__FUNCTION__+": internal error #2!"))
         return;
      //--- allocation
      ArrayResizeAL(thresholds,1);
      //--- change values
      thresholds[0]=0.5*(a[ties[j-1]]+a[ties[j]]);
      ni=2;
      cve=0;
      //--- calculation
      for(i=0;i<=nc-1;i++)
         cnt[i]=0;
      for(i=0;i<=j-1;i++)
         TieAddC(c,ties,i,nc,cnt);
      cve=cve+GetCV(cnt,nc);
      for(i=0;i<=nc-1;i++)
         cnt[i]=0;
      for(i=j;i<=tiecount-1;i++)
         TieAddC(c,ties,i,nc,cnt);
      cve=cve+GetCV(cnt,nc);
     }
   else
     {
      //--- General case: 2 or more intervals
      ArrayResizeAL(thresholds,koptimal);
      ni=koptimal+1;
      cve=cv[koptimal][tiecount-1];
      jl=splits[koptimal][tiecount-1];
      jr=tiecount-1;
      //--- calculation
      for(k=koptimal;k>=1;k--)
        {
         thresholds[k-1]=0.5*(a[ties[jl-1]]+a[ties[jl]]);
         jr=jl-1;
         jl=splits[k-1][jl-1];
        }
     }
  }
//+------------------------------------------------------------------+
//| Internal function                                                |
//+------------------------------------------------------------------+
static double CBdSS::XLnY(const double x,const double y)
  {
//--- check
   if(x==0.0)
      return(0);
//--- return result
   return(x*MathLog(y));
  }
//+------------------------------------------------------------------+
//| Internal function,                                               |
//| returns number of samples of class I in Cnt[I]                   |
//+------------------------------------------------------------------+
static double CBdSS::GetCV(int &cnt[],const int nc)
  {
//--- create variables
   double result=0;
   int    i=0;
   double s=0;
//--- calculation
   s=0;
   for(i=0;i<=nc-1;i++)
      s=s+cnt[i];
//--- get result
   result=0;
   for(i=0;i<=nc-1;i++)
      result=result-XLnY(cnt[i],cnt[i]/(s+nc-1));
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Internal function,adds number of samples of class I in tie NTie  |
//| to Cnt[I]                                                        |
//+------------------------------------------------------------------+
static void CBdSS::TieAddC(int &c[],int &ties[],const int ntie,const int nc,
                           int &cnt[])
  {
//--- create a variable
   int i=0;
//--- calculation
   for(i=ties[ntie];i<=ties[ntie+1]-1;i++)
      cnt[c[i]]=cnt[c[i]]+1;
  }
//+------------------------------------------------------------------+
//| Internal function,subtracts number of samples of class I in tie  |
//| NTie to Cnt[I]                                                   |
//+------------------------------------------------------------------+
static void CBdSS::TieSubC(int &c[],int &ties[],const int ntie,const int nc,
                           int &cnt[])
  {
//--- create a variable
   int i=0;
//--- calculation
   for(i=ties[ntie];i<=ties[ntie+1]-1;i++)
      cnt[c[i]]=cnt[c[i]]-1;
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CDForest                                     |
//+------------------------------------------------------------------+
class CDecisionForest
  {
public:
   int               m_nvars;
   int               m_nclasses;
   int               m_ntrees;
   int               m_bufsize;
   double            m_trees[];
   //--- constructor, destructor
                     CDecisionForest(void);
                    ~CDecisionForest(void);
   //--- copy
   void              Copy(CDecisionForest &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDecisionForest::CDecisionForest(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CDecisionForest::~CDecisionForest(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CDecisionForest::Copy(CDecisionForest &obj)
  {
//--- copy variables
   m_nvars=obj.m_nvars;
   m_nclasses=obj.m_nclasses;
   m_ntrees=obj.m_ntrees;
   m_bufsize=obj.m_bufsize;
//--- copy array
   ArrayCopy(m_trees,obj.m_trees);
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CDecisionForest                  |
//+------------------------------------------------------------------+
class CDecisionForestShell
  {
private:
   CDecisionForest   m_innerobj;
public:
   //--- constructors, destructor
                     CDecisionForestShell(void);
                     CDecisionForestShell(CDecisionForest &obj);
                    ~CDecisionForestShell(void);
   //--- method
   CDecisionForest *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDecisionForestShell::CDecisionForestShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CDecisionForestShell::CDecisionForestShell(CDecisionForest &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CDecisionForestShell::~CDecisionForestShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CDecisionForest *CDecisionForestShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CDForest                                     |
//+------------------------------------------------------------------+
class CDFReport
  {
public:
   //--- variables
   double            m_relclserror;
   double            m_avgce;
   double            m_rmserror;
   double            m_avgerror;
   double            m_avgrelerror;
   double            m_oobrelclserror;
   double            m_oobavgce;
   double            m_oobrmserror;
   double            m_oobavgerror;
   double            m_oobavgrelerror;
   //--- constructor, destructor
                     CDFReport(void);
                    ~CDFReport(void);
   //--- copy
   void              Copy(CDFReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDFReport::CDFReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CDFReport::~CDFReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CDFReport::Copy(CDFReport &obj)
  {
//--- copy variables
   m_relclserror=obj.m_relclserror;
   m_avgce=obj.m_avgce;
   m_rmserror=obj.m_rmserror;
   m_avgerror=obj.m_avgerror;
   m_avgrelerror=obj.m_avgrelerror;
   m_oobrelclserror=obj.m_oobrelclserror;
   m_oobavgce=obj.m_oobavgce;
   m_oobrmserror=obj.m_oobrmserror;
   m_oobavgerror=obj.m_oobavgerror;
   m_oobavgrelerror=obj.m_oobavgrelerror;
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CDFReport                        |
//+------------------------------------------------------------------+
class CDFReportShell
  {
private:
   CDFReport         m_innerobj;
public:
   //--- constructors, destructor
                     CDFReportShell(void);
                     CDFReportShell(CDFReport &obj);
                    ~CDFReportShell(void);
   //--- methods
   double            GetRelClsError(void);
   void              SetRelClsError(const double d);
   double            GetAvgCE(void);
   void              SetAvgCE(const double d);
   double            GetRMSError(void);
   void              SetRMSError(const double d);
   double            GetAvgError(void);
   void              SetAvgError(const double d);
   double            GetAvgRelError(void);
   void              SetAvgRelError(const double d);
   double            GetOOBRelClsError(void);
   void              SetOOBRelClsError(const double d);
   double            GetOOBAvgCE(void);
   void              SetOOBAvgCE(const double d);
   double            GetOOBRMSError(void);
   void              SetOOBRMSError(const double d);
   double            GetOOBAvgError(void);
   void              SetOOBAvgError(const double d);
   double            GetOOBAvgRelError(void);
   void              SetOOBAvgRelError(const double d);
   CDFReport        *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDFReportShell::CDFReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CDFReportShell::CDFReportShell(CDFReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CDFReportShell::~CDFReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable relclserror                    |
//+------------------------------------------------------------------+
double CDFReportShell::GetRelClsError(void)
  {
//--- return result
   return(m_innerobj.m_relclserror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable relclserror                   |
//+------------------------------------------------------------------+
void CDFReportShell::SetRelClsError(const double d)
  {
//--- change value
   m_innerobj.m_relclserror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgce                          |
//+------------------------------------------------------------------+
double CDFReportShell::GetAvgCE(void)
  {
//--- return result
   return(m_innerobj.m_avgce);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgce                         |
//+------------------------------------------------------------------+
void CDFReportShell::SetAvgCE(const double d)
  {
//--- change value
   m_innerobj.m_avgce=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable rmserror                       |
//+------------------------------------------------------------------+
double CDFReportShell::GetRMSError(void)
  {
//--- return result
   return(m_innerobj.m_rmserror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable rmserror                      |
//+------------------------------------------------------------------+
void CDFReportShell::SetRMSError(const double d)
  {
//--- change value
   m_innerobj.m_rmserror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgerror                       |
//+------------------------------------------------------------------+
double CDFReportShell::GetAvgError(void)
  {
//--- return result
   return(m_innerobj.m_avgerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgerror                      |
//+------------------------------------------------------------------+
void CDFReportShell::SetAvgError(const double d)
  {
//--- change value
   m_innerobj.m_avgerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgrelerror                    |
//+------------------------------------------------------------------+
double CDFReportShell::GetAvgRelError(void)
  {
//--- return result
   return(m_innerobj.m_avgrelerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgrelerror                   |
//+------------------------------------------------------------------+
void CDFReportShell::SetAvgRelError(const double d)
  {
//--- change value
   m_innerobj.m_avgrelerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable oobrelclserror                 |
//+------------------------------------------------------------------+
double CDFReportShell::GetOOBRelClsError(void)
  {
//--- return result
   return(m_innerobj.m_oobrelclserror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable oobrelclserror                |
//+------------------------------------------------------------------+
void CDFReportShell::SetOOBRelClsError(const double d)
  {
//--- change value
   m_innerobj.m_oobrelclserror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable oobavgce                       |
//+------------------------------------------------------------------+
double CDFReportShell::GetOOBAvgCE(void)
  {
//--- return result
   return(m_innerobj.m_oobavgce);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable oobavgce                      |
//+------------------------------------------------------------------+
void CDFReportShell::SetOOBAvgCE(const double d)
  {
//--- change value
   m_innerobj.m_oobavgce=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable oobrmserror                    |
//+------------------------------------------------------------------+
double CDFReportShell::GetOOBRMSError(void)
  {
//--- return result
   return(m_innerobj.m_oobrmserror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable oobrmserror                   |
//+------------------------------------------------------------------+
void CDFReportShell::SetOOBRMSError(const double d)
  {
//--- change value
   m_innerobj.m_oobrmserror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable oobavgerror                    |
//+------------------------------------------------------------------+
double CDFReportShell::GetOOBAvgError(void)
  {
//--- return result
   return(m_innerobj.m_oobavgerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable oobavgerror                   |
//+------------------------------------------------------------------+
void CDFReportShell::SetOOBAvgError(const double d)
  {
//--- change value
   m_innerobj.m_oobavgerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable oobavgrelerror                 |
//+------------------------------------------------------------------+
double CDFReportShell::GetOOBAvgRelError(void)
  {
//--- return result
   return(m_innerobj.m_oobavgrelerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable oobavgrelerror                |
//+------------------------------------------------------------------+
void CDFReportShell::SetOOBAvgRelError(const double d)
  {
//--- change value
   m_innerobj.m_oobavgrelerror=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CDFReport *CDFReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CDForest                                     |
//+------------------------------------------------------------------+
class CDFInternalBuffers
  {
public:
   //--- arrays
   double            m_treebuf[];
   int               m_idxbuf[];
   double            m_tmpbufr[];
   double            m_tmpbufr2[];
   int               m_tmpbufi[];
   int               m_classibuf[];
   double            m_sortrbuf[];
   double            m_sortrbuf2[];
   int               m_sortibuf[];
   int               m_varpool[];
   bool              m_evsbin[];
   double            m_evssplits[];
   //--- constructor, destructor
                     CDFInternalBuffers(void);
                    ~CDFInternalBuffers(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDFInternalBuffers::CDFInternalBuffers(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CDFInternalBuffers::~CDFInternalBuffers(void)
  {

  }
//+------------------------------------------------------------------+
//| Decision forest class                                            |
//+------------------------------------------------------------------+
class CDForest
  {
private:
   //--- private methods
   static int        DFClsError(CDecisionForest &df,CMatrixDouble &xy,const int npoints);
   static void       DFProcessInternal(CDecisionForest &df,const int offs,double &x[],double &y[]);
   static void       DFBuildTree(CMatrixDouble &xy,const int npoints,const int nvars,const int nclasses,const int nfeatures,const int nvarsinpool,const int flags,CDFInternalBuffers &bufs);
   static void       DFBuildTreeRec(CMatrixDouble &xy,const int npoints,const int nvars,const int nclasses,const int nfeatures,int nvarsinpool,const int flags,int &numprocessed,const int idx1,const int idx2,CDFInternalBuffers &bufs);
   static void       DFSplitC(double &x[],int &c[],int &cntbuf[],const int n,const int nc,const int flags,int &info,double &threshold,double &e,double &sortrbuf[],int &sortibuf[]);
   static void       DFSplitR(double &x[],double &y[],const int n,const int flags,int &info,double &threshold,double &e,double &sortrbuf[],double &sortrbuf2[]);
public:
   //--- class constants
   static const int  m_innernodewidth;
   static const int  m_leafnodewidth;
   static const int  m_dfusestrongsplits;
   static const int  m_dfuseevs;
   static const int  m_dffirstversion;
   //--- constructor, destructor
                     CDForest(void);
                    ~CDForest(void);
   //--- public methods
   static void       DFBuildRandomDecisionForest(CMatrixDouble &xy,const int npoints,const int nvars,const int nclasses,const int ntrees,const double r,int &info,CDecisionForest &df,CDFReport &rep);
   static void       DFBuildRandomDecisionForestX1(CMatrixDouble &xy,const int npoints,const int nvars,const int nclasses,const int ntrees,const int nrndvars,const double r,int &info,CDecisionForest &df,CDFReport &rep);
   static void       DFBuildInternal(CMatrixDouble &xy,const int npoints,const int nvars,const int nclasses,const int ntrees,const int samplesize,const int nfeatures,const int flags,int &info,CDecisionForest &df,CDFReport &rep);
   static void       DFProcess(CDecisionForest &df,double &x[],double &y[]);
   static void       DFProcessI(CDecisionForest &df,double &x[],double &y[]);
   static double     DFRelClsError(CDecisionForest &df,CMatrixDouble &xy,const int npoints);
   static double     DFAvgCE(CDecisionForest &df,CMatrixDouble &xy,const int npoints);
   static double     DFRMSError(CDecisionForest &df,CMatrixDouble &xy,const int npoints);
   static double     DFAvgError(CDecisionForest &df,CMatrixDouble &xy,const int npoints);
   static double     DFAvgRelError(CDecisionForest &df,CMatrixDouble &xy,const int npoints);
   static void       DFCopy(CDecisionForest &df1,CDecisionForest &df2);
   static void       DFAlloc(CSerializer &s,CDecisionForest &forest);
   static void       DFSerialize(CSerializer &s,CDecisionForest &forest);
   static void       DFUnserialize(CSerializer &s,CDecisionForest &forest);
  };
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const int CDForest::m_innernodewidth=3;
const int CDForest::m_leafnodewidth=2;
const int CDForest::m_dfusestrongsplits=1;
const int CDForest::m_dfuseevs=2;
const int CDForest::m_dffirstversion=0;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDForest::CDForest(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CDForest::~CDForest(void)
  {

  }
//+------------------------------------------------------------------+
//| This subroutine builds random decision forest.                   |
//| INPUT PARAMETERS:                                                |
//|     XY          -   training set                                 |
//|     NPoints     -   training set size, NPoints>=1                |
//|     NVars       -   number of independent variables, NVars>=1    |
//|     NClasses    -   task type:                                   |
//|                     * NClasses=1 - regression task with one      |
//|                                    dependent variable            |
//|                     * NClasses>1 - classification task with      |
//|                                    NClasses classes.             |
//|     NTrees      -   number of trees in a forest, NTrees>=1.      |
//|                     recommended values: 50-100.                  |
//|     R           -   percent of a training set used to build      |
//|                     individual trees. 0<R<=1.                    |
//|                     recommended values: 0.1 <= R <= 0.66.        |
//| OUTPUT PARAMETERS:                                               |
//|     Info        -   return code:                                 |
//|                     * -2, if there is a point with class number  |
//|                           outside of [0..NClasses-1].            |
//|                     * -1, if incorrect parameters was passed     |
//|                           (NPoints<1, NVars<1, NClasses<1,       |
//|                           NTrees<1, R<=0 or R>1).                |
//|                     *  1, if task has been solved                |
//|     DF          -   model built                                  |
//|     Rep         -   training report, contains error on a training|
//|                     set and out-of-bag estimates of              |
//|                     generalization error.                        |
//+------------------------------------------------------------------+
static void CDForest::DFBuildRandomDecisionForest(CMatrixDouble &xy,
                                                  const int npoints,
                                                  const int nvars,
                                                  const int nclasses,
                                                  const int ntrees,
                                                  const double r,int &info,
                                                  CDecisionForest &df,
                                                  CDFReport &rep)
  {
//--- create a variable
   int samplesize=0;
//--- initialization
   info=0;
//--- check
   if(r<=0.0 || r>1.0)
     {
      info=-1;
      return;
     }
//--- calculation
   samplesize=(int)(MathMax((int)MathRound(r*npoints),1));
//--- function call
   DFBuildInternal(xy,npoints,nvars,nclasses,ntrees,samplesize,(int)(MathMax(nvars/2,1)),m_dfusestrongsplits+m_dfuseevs,info,df,rep);
  }
//+------------------------------------------------------------------+
//| This subroutine builds random decision forest.                   |
//| This function gives ability to tune number of variables used when|
//| choosing best split.                                             |
//| INPUT PARAMETERS:                                                |
//|     XY          -   training set                                 |
//|     NPoints     -   training set size, NPoints>=1                |
//|     NVars       -   number of independent variables, NVars>=1    |
//|     NClasses    -   task type:                                   |
//|                     * NClasses=1 - regression task with one      |
//|                                    dependent variable            |
//|                     * NClasses>1 - classification task with      |
//|                                    NClasses classes.             |
//|     NTrees      -   number of trees in a forest, NTrees>=1.      |
//|                     recommended values: 50-100.                  |
//|     NRndVars    -   number of variables used when choosing best  |
//|                     split                                        |
//|     R           -   percent of a training set used to build      |
//|                     individual trees. 0<R<=1.                    |
//|                     recommended values: 0.1 <= R <= 0.66.        |
//| OUTPUT PARAMETERS:                                               |
//|     Info        -   return code:                                 |
//|                     * -2, if there is a point with class number  |
//|                           outside of [0..NClasses-1].            |
//|                     * -1, if incorrect parameters was passed     |
//|                           (NPoints<1, NVars<1, NClasses<1,       |
//|                           NTrees<1, R<=0 or R>1).                |
//|                     *  1, if task has been solved                |
//|     DF          -   model built                                  |
//|     Rep         -   training report, contains error on a training| 
//|                     set and out-of-bag estimates of              |
//|                     generalization error.                        |
//+------------------------------------------------------------------+
static void CDForest::DFBuildRandomDecisionForestX1(CMatrixDouble &xy,
                                                    const int npoints,
                                                    const int nvars,
                                                    const int nclasses,
                                                    const int ntrees,
                                                    const int nrndvars,
                                                    const double r,int &info,
                                                    CDecisionForest &df,
                                                    CDFReport &rep)
  {
//--- create a variable
   int samplesize=0;
//--- initialization
   info=0;
//--- check
   if(r<=0.0 || r>1.0)
     {
      info=-1;
      return;
     }
//--- check
   if(nrndvars<=0 || nrndvars>nvars)
     {
      info=-1;
      return;
     }
//--- calculation
   samplesize=(int)(MathMax((int)MathRound(r*npoints),1));
//--- function call
   DFBuildInternal(xy,npoints,nvars,nclasses,ntrees,samplesize,nrndvars,m_dfusestrongsplits+m_dfuseevs,info,df,rep);
  }
//+------------------------------------------------------------------+
//| Class method                                                     |
//+------------------------------------------------------------------+
static void CDForest::DFBuildInternal(CMatrixDouble &xy,const int npoints,
                                      const int nvars,const int nclasses,
                                      const int ntrees,const int samplesize,
                                      const int nfeatures,const int flags,
                                      int &info,CDecisionForest &df,
                                      CDFReport &rep)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    tmpi=0;
   int    lasttreeoffs=0;
   int    offs=0;
   int    ooboffs=0;
   int    treesize=0;
   int    nvarsinpool=0;
   bool   useevs;
   int    oobcnt=0;
   int    oobrelcnt=0;
   double v=0;
   double vmin=0;
   double vmax=0;
   bool   bflag;
   int    i_=0;
   int    i1_=0;
//--- creating arrays
   int    permbuf[];
   double oobbuf[];
   int    oobcntbuf[];
   double x[];
   double y[];
//--- create matrix
   CMatrixDouble xys;
//--- create object of class
   CDFInternalBuffers bufs;
//--- initialization
   info=0;
//--- Test for inputs
   if(npoints<1 || samplesize<1 || samplesize>npoints || nvars<1 || nclasses<1 || ntrees<1 || nfeatures<1)
     {
      info=-1;
      return;
     }
//--- check
   if(nclasses>1)
     {
      for(i=0;i<=npoints-1;i++)
        {
         //--- check
         if((int)MathRound(xy[i][nvars])<0 || (int)MathRound(xy[i][nvars])>=nclasses)
           {
            info=-2;
            return;
           }
        }
     }
//--- change value
   info=1;
//--- Flags
   useevs=flags/m_dfuseevs%2!=0;
//--- Allocate data,prepare header
   treesize=1+m_innernodewidth*(samplesize-1)+m_leafnodewidth*samplesize;
//--- allocation
   ArrayResizeAL(permbuf,npoints);
   ArrayResizeAL(bufs.m_treebuf,treesize);
   ArrayResizeAL(bufs.m_idxbuf,npoints);
   ArrayResizeAL(bufs.m_tmpbufr,npoints);
   ArrayResizeAL(bufs.m_tmpbufr2,npoints);
   ArrayResizeAL(bufs.m_tmpbufi,npoints);
   ArrayResizeAL(bufs.m_sortrbuf,npoints);
   ArrayResizeAL(bufs.m_sortrbuf2,npoints);
   ArrayResizeAL(bufs.m_sortibuf,npoints);
   ArrayResizeAL(bufs.m_varpool,nvars);
   ArrayResizeAL(bufs.m_evsbin,nvars);
   ArrayResizeAL(bufs.m_evssplits,nvars);
   ArrayResizeAL(bufs.m_classibuf,2*nclasses);
   ArrayResizeAL(oobbuf,nclasses*npoints);
   ArrayResizeAL(oobcntbuf,npoints);
   ArrayResizeAL(df.m_trees,ntrees*treesize);
   xys.Resize(samplesize,nvars+1);
   ArrayResizeAL(x,nvars);
   ArrayResizeAL(y,nclasses);
//--- initialization
   for(i=0;i<=npoints-1;i++)
      permbuf[i]=i;
   for(i=0;i<=npoints*nclasses-1;i++)
      oobbuf[i]=0;
   for(i=0;i<=npoints-1;i++)
      oobcntbuf[i]=0;
//--- Prepare variable pool and EVS (extended variable selection/splitting) buffers
//--- (whether EVS is turned on or not):
//--- 1. detect binary variables and pre-calculate splits for them
//--- 2. detect variables with non-distinct values and exclude them from pool
   for(i=0;i<=nvars-1;i++)
      bufs.m_varpool[i]=i;
   nvarsinpool=nvars;
//--- check
   if(useevs)
     {
      for(j=0;j<=nvars-1;j++)
        {
         vmin=xy[0][j];
         vmax=vmin;
         //--- calculation
         for(i=0;i<=npoints-1;i++)
           {
            v=xy[i][j];
            vmin=MathMin(vmin,v);
            vmax=MathMax(vmax,v);
           }
         //--- check
         if(vmin==vmax)
           {
            //--- exclude variable from pool
            bufs.m_varpool[j]=bufs.m_varpool[nvarsinpool-1];
            bufs.m_varpool[nvarsinpool-1]=-1;
            nvarsinpool=nvarsinpool-1;
            continue;
           }
         //--- change value
         bflag=false;
         for(i=0;i<=npoints-1;i++)
           {
            v=xy[i][j];
            //--- check
            if(v!=vmin&&v!=vmax)
              {
               bflag=true;
               break;
              }
           }
         //--- check
         if(bflag)
           {
            //--- non-binary variable
            bufs.m_evsbin[j]=false;
           }
         else
           {
            //--- Prepare
            bufs.m_evsbin[j]=true;
            bufs.m_evssplits[j]=0.5*(vmin+vmax);
            //--- check
            if(bufs.m_evssplits[j]<=vmin)
               bufs.m_evssplits[j]=vmax;
           }
        }
     }
//--- RANDOM FOREST FORMAT
//--- W[0]         -   size of array
//--- W[1]         -   version number
//--- W[2]         -   NVars
//--- W[3]         -   NClasses (1 for regression)
//--- W[4]         -   NTrees
//--- W[5]         -   trees offset
//--- TREE FORMAT
//--- W[Offs]      -   size of sub-array
//---     node info:
//--- W[K+0]       -   variable number        (-1 for leaf mode)
//--- W[K+1]       -   threshold              (class/value for leaf node)
//--- W[K+2]       -   ">=" branch index      (absent for leaf node)
   df.m_nvars=nvars;
   df.m_nclasses=nclasses;
   df.m_ntrees=ntrees;
//--- Build forest
   offs=0;
   for(i=0;i<=ntrees-1;i++)
     {
      //--- Prepare sample
      for(k=0;k<=samplesize-1;k++)
        {
         //--- calculation
         j=k+CMath::RandomInteger(npoints-k);
         tmpi=permbuf[k];
         permbuf[k]=permbuf[j];
         permbuf[j]=tmpi;
         j=permbuf[k];
         for(i_=0;i_<=nvars;i_++)
            xys[k].Set(i_,xy[j][i_]);
        }
      //--- build tree,copy
      DFBuildTree(xys,samplesize,nvars,nclasses,nfeatures,nvarsinpool,flags,bufs);
      //--- calculation
      j=(int)MathRound(bufs.m_treebuf[0]);
      i1_=-offs;
      for(i_=offs;i_<=offs+j-1;i_++)
         df.m_trees[i_]=bufs.m_treebuf[i_+i1_];
      lasttreeoffs=offs;
      offs=offs+j;
      //--- OOB estimates
      for(k=samplesize;k<=npoints-1;k++)
        {
         for(j=0;j<=nclasses-1;j++)
            y[j]=0;
         j=permbuf[k];
         for(i_=0;i_<=nvars-1;i_++)
            x[i_]=xy[j][i_];
         //--- function call
         DFProcessInternal(df,lasttreeoffs,x,y);
         //--- calculation
         i1_=-j*nclasses;
         for(i_=j*nclasses;i_<=(j+1)*nclasses-1;i_++)
            oobbuf[i_]=oobbuf[i_]+y[i_+i1_];
         oobcntbuf[j]=oobcntbuf[j]+1;
        }
     }
   df.m_bufsize=offs;
//--- Normalize OOB results
   for(i=0;i<=npoints-1;i++)
     {
      //--- check
      if(oobcntbuf[i]!=0)
        {
         v=1.0/(double)oobcntbuf[i];
         for(i_=i*nclasses;i_<=i*nclasses+nclasses-1;i_++)
            oobbuf[i_]=v*oobbuf[i_];
        }
     }
//--- Calculate training set estimates
   rep.m_relclserror=DFRelClsError(df,xy,npoints);
   rep.m_avgce=DFAvgCE(df,xy,npoints);
   rep.m_rmserror=DFRMSError(df,xy,npoints);
   rep.m_avgerror=DFAvgError(df,xy,npoints);
   rep.m_avgrelerror=DFAvgRelError(df,xy,npoints);
//--- Calculate OOB estimates.
   rep.m_oobrelclserror=0;
   rep.m_oobavgce=0;
   rep.m_oobrmserror=0;
   rep.m_oobavgerror=0;
   rep.m_oobavgrelerror=0;
   oobcnt=0;
   oobrelcnt=0;
   for(i=0;i<=npoints-1;i++)
     {
      //--- check
      if(oobcntbuf[i]!=0)
        {
         ooboffs=i*nclasses;
         //--- check
         if(nclasses>1)
           {
            //--- classification-specific code
            k=(int)MathRound(xy[i][nvars]);
            tmpi=0;
            for(j=1;j<=nclasses-1;j++)
              {
               //--- check
               if(oobbuf[ooboffs+j]>oobbuf[ooboffs+tmpi])
                  tmpi=j;
              }
            //--- check
            if(tmpi!=k)
               rep.m_oobrelclserror=rep.m_oobrelclserror+1;
            //--- check
            if(oobbuf[ooboffs+k]!=0.0)
               rep.m_oobavgce=rep.m_oobavgce-MathLog(oobbuf[ooboffs+k]);
            else
               rep.m_oobavgce=rep.m_oobavgce-MathLog(CMath::m_minrealnumber);
            //--- calculation
            for(j=0;j<=nclasses-1;j++)
              {
               //--- check
               if(j==k)
                 {
                  rep.m_oobrmserror=rep.m_oobrmserror+CMath::Sqr(oobbuf[ooboffs+j]-1);
                  rep.m_oobavgerror=rep.m_oobavgerror+MathAbs(oobbuf[ooboffs+j]-1);
                  rep.m_oobavgrelerror=rep.m_oobavgrelerror+MathAbs(oobbuf[ooboffs+j]-1);
                  oobrelcnt=oobrelcnt+1;
                 }
               else
                 {
                  rep.m_oobrmserror=rep.m_oobrmserror+CMath::Sqr(oobbuf[ooboffs+j]);
                  rep.m_oobavgerror=rep.m_oobavgerror+MathAbs(oobbuf[ooboffs+j]);
                 }
              }
           }
         else
           {
            //--- regression-specific code
            rep.m_oobrmserror=rep.m_oobrmserror+CMath::Sqr(oobbuf[ooboffs]-xy[i][nvars]);
            rep.m_oobavgerror=rep.m_oobavgerror+MathAbs(oobbuf[ooboffs]-xy[i][nvars]);
            //--- check
            if(xy[i][nvars]!=0.0)
              {
               rep.m_oobavgrelerror=rep.m_oobavgrelerror+MathAbs((oobbuf[ooboffs]-xy[i][nvars])/xy[i][nvars]);
               oobrelcnt=oobrelcnt+1;
              }
           }
         //--- update OOB estimates count.
         oobcnt=oobcnt+1;
        }
     }
//--- check
   if(oobcnt>0)
     {
      //--- change values
      rep.m_oobrelclserror=rep.m_oobrelclserror/oobcnt;
      rep.m_oobavgce=rep.m_oobavgce/oobcnt;
      rep.m_oobrmserror=MathSqrt(rep.m_oobrmserror/(oobcnt*nclasses));
      rep.m_oobavgerror=rep.m_oobavgerror/(oobcnt*nclasses);
      //--- check
      if(oobrelcnt>0)
         rep.m_oobavgrelerror=rep.m_oobavgrelerror/oobrelcnt;
     }
  }
//+------------------------------------------------------------------+
//| Procesing                                                        |
//| INPUT PARAMETERS:                                                |
//|     DF      -   decision forest model                            |
//|     X       -   input vector,  array[0..NVars-1].                |
//| OUTPUT PARAMETERS:                                               |
//|     Y       -   result. Regression estimate when solving         |
//|                 regression task, vector of posterior             |
//|                 probabilities for classification task.           |
//| See also DFProcessI.                                             |
//+------------------------------------------------------------------+
static void CDForest::DFProcess(CDecisionForest &df,double &x[],double &y[])
  {
//--- create variables
   int    offs=0;
   int    i=0;
   double v=0;
   int    i_=0;
//--- Proceed
   if(CAp::Len(y)<df.m_nclasses)
      ArrayResizeAL(y,df.m_nclasses);
//--- initialization
   offs=0;
   for(i=0;i<=df.m_nclasses-1;i++)
      y[i]=0;
   for(i=0;i<=df.m_ntrees-1;i++)
     {
      //--- Process basic tree
      DFProcessInternal(df,offs,x,y);
      //--- Next tree
      offs=offs+(int)MathRound(df.m_trees[offs]);
     }
//--- calculation
   v=1.0/(double)df.m_ntrees;
   for(i_=0;i_<=df.m_nclasses-1;i_++)
      y[i_]=v*y[i_];
  }
//+------------------------------------------------------------------+
//| 'interactive' variant of DFProcess for languages like Python     |
//| which support constructs like "Y = DFProcessI(DF,X)" and         |
//| interactive mode of interpreter                                  |
//| This function allocates new array on each call, so it is         |
//| significantly slower than its 'non-interactive' counterpart, but |
//| it is more convenient when you call it from command line.        |
//+------------------------------------------------------------------+
static void CDForest::DFProcessI(CDecisionForest &df,double &x[],double &y[])
  {
//--- function call
   DFProcess(df,x,y);
  }
//+------------------------------------------------------------------+
//| Relative classification error on the test set                    |
//| INPUT PARAMETERS:                                                |
//|     DF      -   decision forest model                            |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     percent of incorrectly classified cases.                     |
//|     Zero if model solves regression task.                        |
//+------------------------------------------------------------------+
static double CDForest::DFRelClsError(CDecisionForest &df,CMatrixDouble &xy,
                                      const int npoints)
  {
//--- return result
   return((double)DFClsError(df,xy,npoints)/(double)npoints);
  }
//+------------------------------------------------------------------+
//| Average cross-entropy (in bits per element) on the test set      |
//| INPUT PARAMETERS:                                                |
//|     DF      -   decision forest model                            |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     CrossEntropy/(NPoints*LN(2)).                                |
//|     Zero if model solves regression task.                        |
//+------------------------------------------------------------------+
static double CDForest::DFAvgCE(CDecisionForest &df,CMatrixDouble &xy,
                                const int npoints)
  {
//--- create variables
   double result=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    tmpi=0;
   int    i_=0;
//--- creating arrays
   double x[];
   double y[];
//--- allocation
   ArrayResizeAL(x,df.m_nvars);
   ArrayResizeAL(y,df.m_nclasses);
//--- initialization
   result=0;
   for(i=0;i<=npoints-1;i++)
     {
      for(i_=0;i_<=df.m_nvars-1;i_++)
         x[i_]=xy[i][i_];
      //--- function call
      DFProcess(df,x,y);
      //--- check
      if(df.m_nclasses>1)
        {
         //--- classification-specific code
         k=(int)MathRound(xy[i][df.m_nvars]);
         tmpi=0;
         for(j=1;j<=df.m_nclasses-1;j++)
           {
            //--- check
            if(y[j]>(double)(y[tmpi]))
               tmpi=j;
           }
         //--- check
         if(y[k]!=0.0)
            result=result-MathLog(y[k]);
         else
            result=result-MathLog(CMath::m_minrealnumber);
        }
     }
//--- return result
   return(result/npoints);
  }
//+------------------------------------------------------------------+
//| RMS error on the test set                                        |
//| INPUT PARAMETERS:                                                |
//|     DF      -   decision forest model                            |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     root mean square error.                                      |
//|     Its meaning for regression task is obvious. As for           |
//|     classification task,RMS error means error when estimating    |
//|     posterior probabilities.                                     |
//+------------------------------------------------------------------+
static double CDForest::DFRMSError(CDecisionForest &df,CMatrixDouble &xy,
                                   const int npoints)
  {
//--- create variables
   double result=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    tmpi=0;
   int    i_=0;
//--- creating arrays
   double x[];
   double y[];
//--- allocation
   ArrayResizeAL(x,df.m_nvars);
   ArrayResizeAL(y,df.m_nclasses);
//--- initialization
   result=0;
   for(i=0;i<=npoints-1;i++)
     {
      for(i_=0;i_<=df.m_nvars-1;i_++)
         x[i_]=xy[i][i_];
      //--- function call
      DFProcess(df,x,y);
      //--- check
      if(df.m_nclasses>1)
        {
         //--- classification-specific code
         k=(int)MathRound(xy[i][df.m_nvars]);
         tmpi=0;
         for(j=1;j<=df.m_nclasses-1;j++)
           {
            //--- check
            if(y[j]>y[tmpi])
               tmpi=j;
           }
         for(j=0;j<=df.m_nclasses-1;j++)
           {
            //--- check
            if(j==k)
               result=result+CMath::Sqr(y[j]-1);
            else
               result=result+CMath::Sqr(y[j]);
           }
        }
      else
        {
         //--- regression-specific code
         result=result+CMath::Sqr(y[0]-xy[i][df.m_nvars]);
        }
     }
//--- return result
   return(MathSqrt(result/(npoints*df.m_nclasses)));
  }
//+------------------------------------------------------------------+
//| Average error on the test set                                    |
//| INPUT PARAMETERS:                                                |
//|     DF      -   decision forest model                            |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     Its meaning for regression task is obvious. As for           |
//|     classification task, it means average error when estimating  |
//|     posterior probabilities.                                     |
//+------------------------------------------------------------------+
static double CDForest::DFAvgError(CDecisionForest &df,CMatrixDouble &xy,
                                   const int npoints)
  {
//--- create variables
   double result=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    i_=0;
//--- creating arrays
   double x[];
   double y[];
//--- allocation
   ArrayResizeAL(x,df.m_nvars);
   ArrayResizeAL(y,df.m_nclasses);
//--- initialization
   result=0;
   for(i=0;i<=npoints-1;i++)
     {
      //--- copy
      for(i_=0;i_<=df.m_nvars-1;i_++)
         x[i_]=xy[i][i_];
      //--- function call
      DFProcess(df,x,y);
      //--- check
      if(df.m_nclasses>1)
        {
         //--- classification-specific code
         k=(int)MathRound(xy[i][df.m_nvars]);
         for(j=0;j<=df.m_nclasses-1;j++)
           {
            //--- check
            if(j==k)
               result=result+MathAbs(y[j]-1);
            else
               result=result+MathAbs(y[j]);
           }
        }
      else
        {
         //--- regression-specific code
         result=result+MathAbs(y[0]-xy[i][df.m_nvars]);
        }
     }
//--- return result
   return(result/(npoints*df.m_nclasses));
  }
//+------------------------------------------------------------------+
//| Average relative error on the test set                           |
//| INPUT PARAMETERS:                                                |
//|     DF      -   decision forest model                            |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     Its meaning for regression task is obvious. As for           |
//|     classification task, it means average relative error when    |
//|     estimating posterior probability of belonging to the correct |
//|     class.                                                       |
//+------------------------------------------------------------------+
static double CDForest::DFAvgRelError(CDecisionForest &df,CMatrixDouble &xy,
                                      const int npoints)
  {
//--- create variables
   double result=0;
   int    relcnt=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    i_=0;
//--- creating arrays
   double x[];
   double y[];
//--- allocation
   ArrayResizeAL(x,df.m_nvars);
   ArrayResizeAL(y,df.m_nclasses);
//--- initialization
   result=0;
   relcnt=0;
   for(i=0;i<=npoints-1;i++)
     {
      //--- copy
      for(i_=0;i_<=df.m_nvars-1;i_++)
         x[i_]=xy[i][i_];
      //--- function call
      DFProcess(df,x,y);
      //--- check
      if(df.m_nclasses>1)
        {
         //--- classification-specific code
         k=(int)MathRound(xy[i][df.m_nvars]);
         for(j=0;j<=df.m_nclasses-1;j++)
           {
            //--- check
            if(j==k)
              {
               result=result+MathAbs(y[j]-1);
               relcnt=relcnt+1;
              }
           }
        }
      else
        {
         //--- regression-specific code
         if(xy[i][df.m_nvars]!=0.0)
           {
            result=result+MathAbs((y[0]-xy[i][df.m_nvars])/xy[i][df.m_nvars]);
            relcnt=relcnt+1;
           }
        }
     }
//--- check
   if(relcnt>0)
      result=result/relcnt;
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Copying of DecisionForest strucure                               |
//| INPUT PARAMETERS:                                                |
//|     DF1 -   original                                             |
//| OUTPUT PARAMETERS:                                               |
//|     DF2 -   copy                                                 |
//+------------------------------------------------------------------+
static void CDForest::DFCopy(CDecisionForest &df1,CDecisionForest &df2)
  {
//--- create a variable
   int i_=0;
//--- change values
   df2.m_nvars=df1.m_nvars;
   df2.m_nclasses=df1.m_nclasses;
   df2.m_ntrees=df1.m_ntrees;
   df2.m_bufsize=df1.m_bufsize;
//--- allocation
   ArrayResizeAL(df2.m_trees,df1.m_bufsize);
//--- copy
   for(i_=0;i_<=df1.m_bufsize-1;i_++)
      df2.m_trees[i_]=df1.m_trees[i_];
  }
//+------------------------------------------------------------------+
//| Serializer: allocation                                           |
//+------------------------------------------------------------------+
static void CDForest::DFAlloc(CSerializer &s,CDecisionForest &forest)
  {
//--- preparation to serialize
   s.Alloc_Entry();
   s.Alloc_Entry();
   s.Alloc_Entry();
   s.Alloc_Entry();
   s.Alloc_Entry();
   s.Alloc_Entry();
//--- function call
   CApServ::AllocRealArray(s,forest.m_trees,forest.m_bufsize);
  }
//+------------------------------------------------------------------+
//| Serializer: serialization                                        |
//+------------------------------------------------------------------+
static void CDForest::DFSerialize(CSerializer &s,CDecisionForest &forest)
  {
//--- serializetion
   s.Serialize_Int(CSCodes::GetRDFSerializationCode());
   s.Serialize_Int(m_dffirstversion);
   s.Serialize_Int(forest.m_nvars);
   s.Serialize_Int(forest.m_nclasses);
   s.Serialize_Int(forest.m_ntrees);
   s.Serialize_Int(forest.m_bufsize);
//--- function call
   CApServ::SerializeRealArray(s,forest.m_trees,forest.m_bufsize);
  }
//+------------------------------------------------------------------+
//| Serializer: unserialization                                      |
//+------------------------------------------------------------------+
static void CDForest::DFUnserialize(CSerializer &s,CDecisionForest &forest)
  {
//--- create variables
   int i0=0;
   int i1=0;
//--- check correctness of header
   i0=s.Unserialize_Int();
//--- check
   if(!CAp::Assert(i0==CSCodes::GetRDFSerializationCode(),__FUNCTION__+": stream header corrupted"))
      return;
//--- unserializetion
   i1=s.Unserialize_Int();
//--- check
   if(!CAp::Assert(i1==m_dffirstversion,__FUNCTION__+": stream header corrupted"))
      return;
//--- Unserialize data
   forest.m_nvars=s.Unserialize_Int();
   forest.m_nclasses=s.Unserialize_Int();
   forest.m_ntrees=s.Unserialize_Int();
   forest.m_bufsize=s.Unserialize_Int();
//--- function call
   CApServ::UnserializeRealArray(s,forest.m_trees);
  }
//+------------------------------------------------------------------+
//| Classification error                                             |
//+------------------------------------------------------------------+
static int CDForest::DFClsError(CDecisionForest &df,CMatrixDouble &xy,
                                const int npoints)
  {
//--- create variables
   int result=0;
   int i=0;
   int j=0;
   int k=0;
   int tmpi=0;
   int i_=0;
//--- creating arrays
   double x[];
   double y[];
//--- check
   if(df.m_nclasses<=1)
      return(0);
//--- allocation
   ArrayResizeAL(x,df.m_nvars);
   ArrayResizeAL(y,df.m_nclasses);
//--- initialization
   result=0;
   for(i=0;i<=npoints-1;i++)
     {
      //--- copy
      for(i_=0;i_<=df.m_nvars-1;i_++)
         x[i_]=xy[i][i_];
      //--- function call
      DFProcess(df,x,y);
      //--- change values
      k=(int)MathRound(xy[i][df.m_nvars]);
      tmpi=0;
      for(j=1;j<=df.m_nclasses-1;j++)
        {
         //--- check
         if(y[j]>(double)(y[tmpi]))
            tmpi=j;
        }
      //--- check
      if(tmpi!=k)
         result=result+1;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Internal subroutine for processing one decision tree starting at |
//| Offs                                                             |
//+------------------------------------------------------------------+
static void CDForest::DFProcessInternal(CDecisionForest &df,const int offs,
                                        double &x[],double &y[])
  {
//--- create variables
   int k=0;
   int idx=0;
//--- Set pointer to the root
   k=offs+1;
//--- Navigate through the tree
   while(true)
     {
      //--- check
      if(df.m_trees[k]==-1.0)
        {
         //--- check
         if(df.m_nclasses==1)
            y[0]=y[0]+df.m_trees[k+1];
         else
           {
            idx=(int)MathRound(df.m_trees[k+1]);
            y[idx]=y[idx]+1;
           }
         //--- break the cycle
         break;
        }
      //--- check
      if(x[(int)MathRound(df.m_trees[k])]<df.m_trees[k+1])
         k=k+m_innernodewidth;
      else
         k=offs+(int)MathRound(df.m_trees[k+2]);
     }
  }
//+------------------------------------------------------------------+
//| Builds one decision tree. Just a wrapper for the DFBuildTreeRec. |
//+------------------------------------------------------------------+
static void CDForest::DFBuildTree(CMatrixDouble &xy,const int npoints,
                                  const int nvars,const int nclasses,
                                  const int nfeatures,const int nvarsinpool,
                                  const int flags,CDFInternalBuffers &bufs)
  {
//--- create variables
   int numprocessed=0;
   int i=0;
//--- check
   if(!CAp::Assert(npoints>0))
      return;
//--- Prepare IdxBuf. It stores indices of the training set elements.
//--- When training set is being split,contents of IdxBuf is
//--- correspondingly reordered so we can know which elements belong
//--- to which branch of decision tree.
   for(i=0;i<=npoints-1;i++)
      bufs.m_idxbuf[i]=i;
//--- Recursive procedure
   numprocessed=1;
//--- function call
   DFBuildTreeRec(xy,npoints,nvars,nclasses,nfeatures,nvarsinpool,flags,numprocessed,0,npoints-1,bufs);
//--- change values
   bufs.m_treebuf[0]=numprocessed;
  }
//+------------------------------------------------------------------+
//| Builds one decision tree (internal recursive subroutine)         |
//| Parameters:                                                      |
//|     TreeBuf     -   large enough array,at least TreeSize         |
//|     IdxBuf      -   at least NPoints elements                    |
//|     TmpBufR     -   at least NPoints                             |
//|     TmpBufR2    -   at least NPoints                             |
//|     TmpBufI     -   at least NPoints                             |
//|     TmpBufI2    -   at least NPoints+1                           |
//+------------------------------------------------------------------+
static void CDForest::DFBuildTreeRec(CMatrixDouble &xy,const int npoints,
                                     const int nvars,const int nclasses,
                                     const int nfeatures,int nvarsinpool,
                                     const int flags,int &numprocessed,
                                     const int idx1,const int idx2,
                                     CDFInternalBuffers &bufs)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   bool   bflag;
   int    i1=0;
   int    i2=0;
   int    info=0;
   double sl=0;
   double sr=0;
   double w=0;
   int    idxbest=0;
   double ebest=0;
   double tbest=0;
   int    varcur=0;
   double s=0;
   double v=0;
   double v1=0;
   double v2=0;
   double threshold=0;
   int    oldnp=0;
   double currms=0;
   bool   useevs;
//--- these initializers are not really necessary,
//--- but without them compiler complains about uninitialized locals
   tbest=0;
//--- Prepare
   if(!CAp::Assert(npoints>0))
      return;
//--- check
   if(!CAp::Assert(idx2>=idx1))
      return;
   useevs=flags/m_dfuseevs%2!=0;
//--- Leaf node
   if(idx2==idx1)
     {
      bufs.m_treebuf[numprocessed]=-1;
      bufs.m_treebuf[numprocessed+1]=xy[bufs.m_idxbuf[idx1]][nvars];
      numprocessed=numprocessed+m_leafnodewidth;
      //--- exit the function
      return;
     }
//--- Non-leaf node.
//--- Select random variable,prepare split:
//--- 1. prepare default solution - no splitting,class at random
//--- 2. investigate possible splits,compare with default/best
   idxbest=-1;
//--- check
   if(nclasses>1)
     {
      //--- default solution for classification
      for(i=0;i<=nclasses-1;i++)
         bufs.m_classibuf[i]=0;
      s=idx2-idx1+1;
      for(i=idx1;i<=idx2;i++)
        {
         j=(int)MathRound(xy[bufs.m_idxbuf[i]][nvars]);
         bufs.m_classibuf[j]=bufs.m_classibuf[j]+1;
        }
      //--- calculation
      ebest=0;
      for(i=0;i<=nclasses-1;i++)
         ebest=ebest+bufs.m_classibuf[i]*CMath::Sqr(1-bufs.m_classibuf[i]/s)+(s-bufs.m_classibuf[i])*CMath::Sqr(bufs.m_classibuf[i]/s);
      ebest=MathSqrt(ebest/(nclasses*(idx2-idx1+1)));
     }
   else
     {
      //--- default solution for regression
      v=0;
      for(i=idx1;i<=idx2;i++)
         v=v+xy[bufs.m_idxbuf[i]][nvars];
      v=v/(idx2-idx1+1);
      //--- calculation
      ebest=0;
      for(i=idx1;i<=idx2;i++)
         ebest=ebest+CMath::Sqr(xy[bufs.m_idxbuf[i]][nvars]-v);
      ebest=MathSqrt(ebest/(idx2-idx1+1));
     }
//--- change value
   i=0;
//--- cycle
   while(i<=MathMin(nfeatures,nvarsinpool)-1)
     {
      //--- select variables from pool
      j=i+CMath::RandomInteger(nvarsinpool-i);
      k=bufs.m_varpool[i];
      bufs.m_varpool[i]=bufs.m_varpool[j];
      bufs.m_varpool[j]=k;
      varcur=bufs.m_varpool[i];
      //--- load variable values to working array
      //--- apply EVS preprocessing: if all variable values are same,
      //--- variable is excluded from pool.
      //--- This is necessary for binary pre-splits (see later) to work.
      for(j=idx1;j<=idx2;j++)
         bufs.m_tmpbufr[j-idx1]=xy[bufs.m_idxbuf[j]][varcur];
      //--- check
      if(useevs)
        {
         bflag=false;
         v=bufs.m_tmpbufr[0];
         for(j=0;j<=idx2-idx1;j++)
           {
            //--- check
            if(bufs.m_tmpbufr[j]!=v)
              {
               bflag=true;
               break;
              }
           }
         //--- check
         if(!bflag)
           {
            //--- exclude variable from pool,
            //--- go to the next iteration.
            //--- I is not increased.
            k=bufs.m_varpool[i];
            bufs.m_varpool[i]=bufs.m_varpool[nvarsinpool-1];
            bufs.m_varpool[nvarsinpool-1]=k;
            nvarsinpool=nvarsinpool-1;
            continue;
           }
        }
      //--- load labels to working array
      if(nclasses>1)
        {
         for(j=idx1;j<=idx2;j++)
            bufs.m_tmpbufi[j-idx1]=(int)MathRound(xy[bufs.m_idxbuf[j]][nvars]);
        }
      else
        {
         for(j=idx1;j<=idx2;j++)
            bufs.m_tmpbufr2[j-idx1]=xy[bufs.m_idxbuf[j]][nvars];
        }
      //--- calculate split
      if(useevs && bufs.m_evsbin[varcur])
        {
         //--- Pre-calculated splits for binary variables.
         //--- Threshold is already known,just calculate RMS error
         threshold=bufs.m_evssplits[varcur];
         //--- check
         if(nclasses>1)
           {
            //--- classification-specific code
            for(j=0;j<=2*nclasses-1;j++)
               bufs.m_classibuf[j]=0;
            //--- change values
            sl=0;
            sr=0;
            //--- calculation
            for(j=0;j<=idx2-idx1;j++)
              {
               k=bufs.m_tmpbufi[j];
               //--- check
               if(bufs.m_tmpbufr[j]<threshold)
                 {
                  bufs.m_classibuf[k]=bufs.m_classibuf[k]+1;
                  sl=sl+1;
                 }
               else
                 {
                  bufs.m_classibuf[k+nclasses]=bufs.m_classibuf[k+nclasses]+1;
                  sr=sr+1;
                 }
              }
            //--- check
            if(!CAp::Assert(sl!=0.0 && sr!=0.0,__FUNCTION__+": something strange!"))
               return;
            //--- change values
            currms=0;
            //--- calculation
            for(j=0;j<=nclasses-1;j++)
              {
               w=bufs.m_classibuf[j];
               currms=currms+w*CMath::Sqr(w/sl-1);
               currms=currms+(sl-w)*CMath::Sqr(w/sl);
               w=bufs.m_classibuf[nclasses+j];
               currms=currms+w*CMath::Sqr(w/sr-1);
               currms=currms+(sr-w)*CMath::Sqr(w/sr);
              }
            currms=MathSqrt(currms/(nclasses*(idx2-idx1+1)));
           }
         else
           {
            //--- regression-specific code
            sl=0;
            sr=0;
            v1=0;
            v2=0;
            //--- calculation
            for(j=0;j<=idx2-idx1;j++)
              {
               //--- check
               if(bufs.m_tmpbufr[j]<threshold)
                 {
                  v1=v1+bufs.m_tmpbufr2[j];
                  sl=sl+1;
                 }
               else
                 {
                  v2=v2+bufs.m_tmpbufr2[j];
                  sr=sr+1;
                 }
              }
            //--- check
            if(!CAp::Assert(sl!=0.0 && sr!=0.0,__FUNCTION__+": something strange!"))
               return;
            //--- change values
            v1=v1/sl;
            v2=v2/sr;
            currms=0;
            for(j=0;j<=idx2-idx1;j++)
              {
               //--- check
               if(bufs.m_tmpbufr[j]<threshold)
                  currms=currms+CMath::Sqr(v1-bufs.m_tmpbufr2[j]);
               else
                  currms=currms+CMath::Sqr(v2-bufs.m_tmpbufr2[j]);
              }
            currms=MathSqrt(currms/(idx2-idx1+1));
           }
         //--- change value
         info=1;
        }
      else
        {
         //--- Generic splits
         if(nclasses>1)
            DFSplitC(bufs.m_tmpbufr,bufs.m_tmpbufi,bufs.m_classibuf,idx2-idx1+1,nclasses,m_dfusestrongsplits,info,threshold,currms,bufs.m_sortrbuf,bufs.m_sortibuf);
         else
            DFSplitR(bufs.m_tmpbufr,bufs.m_tmpbufr2,idx2-idx1+1,m_dfusestrongsplits,info,threshold,currms,bufs.m_sortrbuf,bufs.m_sortrbuf2);
        }
      //--- check
      if(info>0)
        {
         //--- check
         if(currms<=ebest)
           {
            ebest=currms;
            idxbest=varcur;
            tbest=threshold;
           }
        }
      //--- Next iteration
      i=i+1;
     }
//--- to split or not to split
   if(idxbest<0)
     {
      //--- All values are same,cannot split.
      bufs.m_treebuf[numprocessed]=-1;
      //--- check
      if(nclasses>1)
        {
         //--- Select random class label (randomness allows us to
         //--- approximate distribution of the classes)
         bufs.m_treebuf[numprocessed+1]=(int)MathRound(xy[bufs.m_idxbuf[idx1+CMath::RandomInteger(idx2-idx1+1)]][nvars]);
        }
      else
        {
         //--- Select average (for regression task).
         v=0;
         for(i=idx1;i<=idx2;i++)
            v=v+xy[bufs.m_idxbuf[i]][nvars]/(idx2-idx1+1);
         bufs.m_treebuf[numprocessed+1]=v;
        }
      //--- change value
      numprocessed=numprocessed+m_leafnodewidth;
     }
   else
     {
      //--- we can split
      bufs.m_treebuf[numprocessed]=idxbest;
      bufs.m_treebuf[numprocessed+1]=tbest;
      i1=idx1;
      i2=idx2;
      //--- cycle
      while(i1<=i2)
        {
         //--- Reorder indices so that left partition is in [Idx1..I1-1],
         //--- and right partition is in [I2+1..Idx2]
         if(xy[bufs.m_idxbuf[i1]][idxbest]<tbest)
           {
            i1=i1+1;
            continue;
           }
         //--- check
         if(xy[bufs.m_idxbuf[i2]][idxbest]>=tbest)
           {
            i2=i2-1;
            continue;
           }
         //--- change values
         j=bufs.m_idxbuf[i1];
         bufs.m_idxbuf[i1]=bufs.m_idxbuf[i2];
         bufs.m_idxbuf[i2]=j;
         i1=i1+1;
         i2=i2-1;
        }
      //--- change values
      oldnp=numprocessed;
      numprocessed=numprocessed+m_innernodewidth;
      //--- function call
      DFBuildTreeRec(xy,npoints,nvars,nclasses,nfeatures,nvarsinpool,flags,numprocessed,idx1,i1-1,bufs);
      bufs.m_treebuf[oldnp+2]=numprocessed;
      //--- function call
      DFBuildTreeRec(xy,npoints,nvars,nclasses,nfeatures,nvarsinpool,flags,numprocessed,i2+1,idx2,bufs);
     }
  }
//+------------------------------------------------------------------+
//| Makes split on attribute                                         |
//+------------------------------------------------------------------+
static void CDForest::DFSplitC(double &x[],int &c[],int &cntbuf[],const int n,
                               const int nc,const int flags,int &info,
                               double &threshold,double &e,double &sortrbuf[],
                               int &sortibuf[])
  {
//--- create variables
   int    i=0;
   int    neq=0;
   int    nless=0;
   int    ngreater=0;
   int    q=0;
   int    qmin=0;
   int    qmax=0;
   int    qcnt=0;
   double cursplit=0;
   int    nleft=0;
   double v=0;
   double cure=0;
   double w=0;
   double sl=0;
   double sr=0;
//--- initialization
   info=0;
   threshold=0;
   e=0;
//--- function call
   CTSort::TagSortFastI(x,c,sortrbuf,sortibuf,n);
//--- change values
   e=CMath::m_maxrealnumber;
   threshold=0.5*(x[0]+x[n-1]);
   info=-3;
//--- check
   if(flags/m_dfusestrongsplits%2==0)
     {
      //--- weak splits,split at half
      qcnt=2;
      qmin=1;
      qmax=1;
     }
   else
     {
      //--- strong splits: choose best quartile
      qcnt=4;
      qmin=1;
      qmax=3;
     }
   for(q=qmin;q<=qmax;q++)
     {
      //--- change values
      cursplit=x[n*q/qcnt];
      neq=0;
      nless=0;
      ngreater=0;
      //--- calculation
      for(i=0;i<=n-1;i++)
        {
         //--- check
         if(x[i]<cursplit)
            nless=nless+1;
         //--- check
         if(x[i]==cursplit)
            neq=neq+1;
         //--- check
         if(x[i]>cursplit)
            ngreater=ngreater+1;
        }
      //--- check
      if(!CAp::Assert(neq!=0,__FUNCTION__+": NEq=0,something strange!!!"))
         return;
      //--- check
      if(nless!=0 || ngreater!=0)
        {
         //--- set threshold between two partitions, with
         //--- some tweaking to avoid problems with floating point
         //--- arithmetics.
         //--- The problem is that when you calculates C = 0.5*(A+B) there
         //--- can be no C which lies strictly between A and B (for example,
         //--- there is no floating point number which is
         //--- greater than 1 and less than 1+eps). In such situations
         //--- we choose right side as theshold (remember that
         //--- points which lie on threshold falls to the right side).
         if(nless<ngreater)
           {
            cursplit=0.5*(x[nless+neq-1]+x[nless+neq]);
            nleft=nless+neq;
            //--- check
            if(cursplit<=(double)(x[nless+neq-1]))
               cursplit=x[nless+neq];
           }
         else
           {
            cursplit=0.5*(x[nless-1]+x[nless]);
            nleft=nless;
            //--- check
            if(cursplit<=(double)(x[nless-1]))
               cursplit=x[nless];
           }
         //--- change value
         info=1;
         cure=0;
         //--- calculation
         for(i=0;i<=2*nc-1;i++)
            cntbuf[i]=0;
         for(i=0;i<=nleft-1;i++)
            cntbuf[c[i]]=cntbuf[c[i]]+1;
         for(i=nleft;i<=n-1;i++)
            cntbuf[nc+c[i]]=cntbuf[nc+c[i]]+1;
         //--- change values
         sl=nleft;
         sr=n-nleft;
         v=0;
         //--- calculation
         for(i=0;i<=nc-1;i++)
           {
            w=cntbuf[i];
            v=v+w*CMath::Sqr(w/sl-1);
            v=v+(sl-w)*CMath::Sqr(w/sl);
            w=cntbuf[nc+i];
            v=v+w*CMath::Sqr(w/sr-1);
            v=v+(sr-w)*CMath::Sqr(w/sr);
           }
         cure=MathSqrt(v/(nc*n));
         //--- check
         if(cure<e)
           {
            threshold=cursplit;
            e=cure;
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Makes split on attribute                                         |
//+------------------------------------------------------------------+
static void CDForest::DFSplitR(double &x[],double &y[],const int n,const int flags,
                               int &info,double &threshold,double &e,
                               double &sortrbuf[],double &sortrbuf2[])
  {
//--- create variables
   int    i=0;
   int    neq=0;
   int    nless=0;
   int    ngreater=0;
   int    q=0;
   int    qmin=0;
   int    qmax=0;
   int    qcnt=0;
   double cursplit=0;
   int    nleft=0;
   double v=0;
   double cure=0;
//--- initialization
   info=0;
   threshold=0;
   e=0;
//--- function call
   CTSort::TagSortFastR(x,y,sortrbuf,sortrbuf2,n);
//--- change values
   e=CMath::m_maxrealnumber;
   threshold=0.5*(x[0]+x[n-1]);
   info=-3;
//--- check
   if(flags/m_dfusestrongsplits%2==0)
     {
      //--- weak splits,split at half
      qcnt=2;
      qmin=1;
      qmax=1;
     }
   else
     {
      //--- strong splits: choose best quartile
      qcnt=4;
      qmin=1;
      qmax=3;
     }
//--- calculation
   for(q=qmin;q<=qmax;q++)
     {
      //--- change values
      cursplit=x[n*q/qcnt];
      neq=0;
      nless=0;
      ngreater=0;
      for(i=0;i<=n-1;i++)
        {
         //--- check
         if(x[i]<cursplit)
            nless=nless+1;
         //--- check
         if(x[i]==cursplit)
            neq=neq+1;
         //--- check
         if(x[i]>cursplit)
            ngreater=ngreater+1;
        }
      //--- check
      if(!CAp::Assert(neq!=0,__FUNCTION__+": NEq=0,something strange!!!"))
         return;
      //--- check
      if(nless!=0 || ngreater!=0)
        {
         //--- set threshold between two partitions, with
         //--- some tweaking to avoid problems with floating point
         //--- arithmetics.
         //--- The problem is that when you calculates C = 0.5*(A+B) there
         //--- can be no C which lies strictly between A and B (for example,
         //--- there is no floating point number which is
         //--- greater than 1 and less than 1+eps). In such situations
         //--- we choose right side as theshold (remember that
         //--- points which lie on threshold falls to the right side).
         if(nless<ngreater)
           {
            cursplit=0.5*(x[nless+neq-1]+x[nless+neq]);
            nleft=nless+neq;
            //--- check
            if(cursplit<=(double)(x[nless+neq-1]))
               cursplit=x[nless+neq];
           }
         else
           {
            cursplit=0.5*(x[nless-1]+x[nless]);
            nleft=nless;
            //--- check
            if(cursplit<=(double)(x[nless-1]))
               cursplit=x[nless];
           }
         //--- change value
         info=1;
         cure=0;
         v=0;
         //--- calculation
         for(i=0;i<=nleft-1;i++)
            v=v+y[i];
         v=v/nleft;
         for(i=0;i<=nleft-1;i++)
            cure=cure+CMath::Sqr(y[i]-v);
         v=0;
         for(i=nleft;i<=n-1;i++)
            v=v+y[i];
         v=v/(n-nleft);
         for(i=nleft;i<=n-1;i++)
            cure=cure+CMath::Sqr(y[i]-v);
         cure=MathSqrt(cure/n);
         //--- check
         if(cure<e)
           {
            threshold=cursplit;
            e=cure;
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Middle and clusterization                                        |
//+------------------------------------------------------------------+
class CKMeans
  {
private:
   //--- private method
   static bool       SelectCenterPP(CMatrixDouble &xy,const int npoints,const int nvars,CMatrixDouble &centers,bool &cbusycenters[],const int ccnt,double &d2[],double &p[],double &tmp[]);
public:
   //--- constructor, destructor
                     CKMeans(void);
                    ~CKMeans(void);
   //--- public method
   static void       KMeansGenerate(CMatrixDouble &xy,const int npoints,const int nvars,const int k,const int restarts,int &info,CMatrixDouble &c,int &xyc[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CKMeans::CKMeans(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CKMeans::~CKMeans(void)
  {

  }
//+------------------------------------------------------------------+
//| k-means++ clusterization                                         |
//| INPUT PARAMETERS:                                                |
//|     XY          -   dataset, array [0..NPoints-1,0..NVars-1].    |
//|     NPoints     -   dataset size, NPoints>=K                     |
//|     NVars       -   number of variables, NVars>=1                |
//|     K           -   desired number of clusters, K>=1             |
//|     Restarts    -   number of restarts, Restarts>=1              |
//| OUTPUT PARAMETERS:                                               |
//|     Info        -   return code:                                 |
//|                     * -3, if task is degenerate (number of       |
//|                           distinct points is less than K)        |
//|                     * -1, if incorrect                           |
//|                           NPoints/NFeatures/K/Restarts was passed|
//|                     *  1, if subroutine finished successfully    |
//|     C           -   array[0..NVars-1,0..K-1].matrix whose columns|
//|                     store cluster's centers                      |
//|     XYC         -   array[NPoints], which contains cluster       |
//|                     indexes                                      |
//+------------------------------------------------------------------+
static void CKMeans::KMeansGenerate(CMatrixDouble &xy,const int npoints,
                                    const int nvars,const int k,
                                    const int restarts,int &info,
                                    CMatrixDouble &c,int &xyc[])
  {
//--- create variables
   int    i=0;
   int    j=0;
   double e=0;
   double ebest=0;
   double v=0;
   int    cclosest=0;
   bool   waschanges;
   bool   zerosizeclusters;
   int    pass=0;
   int    i_=0;
   double dclosest=0;
//--- creating arrays
   int    xycbest[];
   double x[];
   double tmp[];
   double d2[];
   double p[];
   int    csizes[];
   bool   cbusy[];
   double work[];
//--- create matrix
   CMatrixDouble ct;
   CMatrixDouble ctbest;
//--- initialization
   info=0;
//--- Test parameters
   if(npoints<k || nvars<1 || k<1 || restarts<1)
     {
      info=-1;
      return;
     }
//--- TODO: special case K=1
//--- TODO: special case K=NPoints
   info=1;
//--- Multiple passes of k-means++ algorithm
   ct.Resize(k,nvars);
   ctbest.Resize(k,nvars);
   ArrayResizeAL(xyc,npoints);
   ArrayResizeAL(xycbest,npoints);
   ArrayResizeAL(d2,npoints);
   ArrayResizeAL(p,npoints);
   ArrayResizeAL(tmp,nvars);
   ArrayResizeAL(csizes,k);
   ArrayResizeAL(cbusy,k);
//--- change value
   ebest=CMath::m_maxrealnumber;
//--- calculation
   for(pass=1;pass<=restarts;pass++)
     {
      //--- Select initial centers  using k-means++ algorithm
      //--- 1. Choose first center at random
      //--- 2. Choose next centers using their distance from centers already chosen
      //--- Note that for performance reasons centers are stored in ROWS of CT,not
      //--- in columns. We'll transpose CT in the end and store it in the C.
      i=CMath::RandomInteger(npoints);
      for(i_=0;i_<=nvars-1;i_++)
         ct[0].Set(i_,xy[i][i_]);
      cbusy[0]=true;
      for(i=1;i<=k-1;i++)
         cbusy[i]=false;
      //--- check
      if(!SelectCenterPP(xy,npoints,nvars,ct,cbusy,k,d2,p,tmp))
        {
         info=-3;
         return;
        }
      //--- Update centers:
      //--- 2. update center positions
      for(i=0;i<=npoints-1;i++)
         xyc[i]=-1;
      //--- cycle
      while(true)
        {
         //--- fill XYC with center numbers
         waschanges=false;
         for(i=0;i<=npoints-1;i++)
           {
            //--- change values
            cclosest=-1;
            dclosest=CMath::m_maxrealnumber;
            for(j=0;j<=k-1;j++)
              {
               //--- calculation
               for(i_=0;i_<=nvars-1;i_++)
                  tmp[i_]=xy[i][i_];
               for(i_=0;i_<=nvars-1;i_++)
                  tmp[i_]=tmp[i_]-ct[j][i_];
               v=0.0;
               for(i_=0;i_<=nvars-1;i_++)
                  v+=tmp[i_]*tmp[i_];
               //--- check
               if(v<dclosest)
                 {
                  cclosest=j;
                  dclosest=v;
                 }
              }
            //--- check
            if(xyc[i]!=cclosest)
               waschanges=true;
            //--- change value
            xyc[i]=cclosest;
           }
         //--- Update centers
         for(j=0;j<=k-1;j++)
            csizes[j]=0;
         for(i=0;i<=k-1;i++)
           {
            for(j=0;j<=nvars-1;j++)
               ct[i].Set(j,0);
           }
         //--- change values
         for(i=0;i<=npoints-1;i++)
           {
            csizes[xyc[i]]=csizes[xyc[i]]+1;
            for(i_=0;i_<=nvars-1;i_++)
               ct[xyc[i]].Set(i_,ct[xyc[i]][i_]+xy[i][i_]);
           }
         zerosizeclusters=false;
         for(i=0;i<=k-1;i++)
           {
            cbusy[i]=csizes[i]!=0;
            zerosizeclusters=zerosizeclusters || csizes[i]==0;
           }
         //--- check
         if(zerosizeclusters)
           {
            //--- Some clusters have zero size - rare,but possible.
            //--- We'll choose new centers for such clusters using k-means++ rule
            //--- and restart algorithm
            if(!SelectCenterPP(xy,npoints,nvars,ct,cbusy,k,d2,p,tmp))
              {
               info=-3;
               return;
              }
            continue;
           }
         //--- copy
         for(j=0;j<=k-1;j++)
           {
            v=1.0/(double)csizes[j];
            for(i_=0;i_<=nvars-1;i_++)
               ct[j].Set(i_,v*ct[j][i_]);
           }
         //--- if nothing has changed during iteration
         if(!waschanges)
            break;
        }
      //--- 3. Calculate E,compare with best centers found so far
      e=0;
      for(i=0;i<=npoints-1;i++)
        {
         for(i_=0;i_<=nvars-1;i_++)
            tmp[i_]=xy[i][i_];
         for(i_=0;i_<=nvars-1;i_++)
            tmp[i_]=tmp[i_]-ct[xyc[i]][i_];
         //--- calculation
         v=0.0;
         for(i_=0;i_<=nvars-1;i_++)
            v+=tmp[i_]*tmp[i_];
         e=e+v;
        }
      //--- check
      if(e<ebest)
        {
         //--- store partition.
         ebest=e;
         //--- function call
         CBlas::CopyMatrix(ct,0,k-1,0,nvars-1,ctbest,0,k-1,0,nvars-1);
         //--- copy
         for(i=0;i<=npoints-1;i++)
            xycbest[i]=xyc[i];
        }
     }
//--- Copy and transpose
   c.Resize(nvars,k);
//--- function call
   CBlas::CopyAndTranspose(ctbest,0,k-1,0,nvars-1,c,0,nvars-1,0,k-1);
//--- copy
   for(i=0;i<=npoints-1;i++)
      xyc[i]=xycbest[i];
  }
//+------------------------------------------------------------------+
//| Select center for a new cluster using k-means++ rule             |
//+------------------------------------------------------------------+
static bool CKMeans::SelectCenterPP(CMatrixDouble &xy,const int npoints,
                                    const int nvars,CMatrixDouble &centers,
                                    bool &cbusycenters[],const int ccnt,
                                    double &d2[],double &p[],double &tmp[])
  {
//--- create variables
   bool   result;
   int    i=0;
   int    j=0;
   int    cc=0;
   double v=0;
   double s=0;
   int    i_=0;
//--- create array
   double busycenters[];
//--- copy
   ArrayCopy(busycenters,cbusycenters);
//--- initialization
   result=true;
//--- calculation
   for(cc=0;cc<=ccnt-1;cc++)
     {
      //--- check
      if(!busycenters[cc])
        {
         //--- fill D2
         for(i=0;i<=npoints-1;i++)
           {
            d2[i]=CMath::m_maxrealnumber;
            for(j=0;j<=ccnt-1;j++)
              {
               //--- check
               if(busycenters[j])
                 {
                  for(i_=0;i_<=nvars-1;i_++)
                     tmp[i_]=xy[i][i_];
                  for(i_=0;i_<=nvars-1;i_++)
                     tmp[i_]=tmp[i_]-centers[j][i_];
                  //--- calculation
                  v=0.0;
                  for(i_=0;i_<=nvars-1;i_++)
                     v+=tmp[i_]*tmp[i_];
                  //--- check
                  if(v<d2[i])
                     d2[i]=v;
                 }
              }
           }
         //--- calculate P (non-cumulative)
         s=0;
         for(i=0;i<=npoints-1;i++)
            s=s+d2[i];
         //--- check
         if(s==0.0)
            return(false);
         //--- change value
         s=1/s;
         for(i_=0;i_<=npoints-1;i_++)
            p[i_]=s*d2[i_];
         //--- choose one of points with probability P
         //--- random number within (0,1) is generated and
         //--- inverse empirical CDF is used to randomly choose a point.
         s=0;
         v=CMath::RandomReal();
         //--- calculation
         for(i=0;i<=npoints-1;i++)
           {
            s=s+p[i];
            //--- check
            if(v<=s || i==npoints-1)
              {
               for(i_=0;i_<=nvars-1;i_++)
                  centers[cc].Set(i_,xy[i][i_]);
               busycenters[cc]=true;
               //--- break the cycle
               break;
              }
           }
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Multiclass Fisher LDA                                            |
//+------------------------------------------------------------------+
class CLDA
  {
public:
   //--- constructor, destructor
                     CLDA(void);
                    ~CLDA(void);
   //--- methods
   static void       FisherLDA(CMatrixDouble &xy,const int npoints,const int nvars,const int nclasses,int &info,double &w[]);
   static void       FisherLDAN(CMatrixDouble &xy,const int npoints,const int nvars,const int nclasses,int &info,CMatrixDouble &w);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLDA::CLDA(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLDA::~CLDA(void)
  {

  }
//+------------------------------------------------------------------+
//| Multiclass Fisher LDA                                            |
//| Subroutine finds coefficients of linear combination which        |
//| optimally separates training set on classes.                     |
//| INPUT PARAMETERS:                                                |
//|     XY          -   training set, array[0..NPoints-1,0..NVars].  |
//|                     First NVars columns store values of          |
//|                     independent variables, next column stores    |
//|                     number of class (from 0 to NClasses-1) which |
//|                     dataset element belongs to. Fractional values|
//|                     are rounded to nearest integer.              |
//|     NPoints     -   training set size, NPoints>=0                |
//|     NVars       -   number of independent variables, NVars>=1    |
//|     NClasses    -   number of classes, NClasses>=2               |
//| OUTPUT PARAMETERS:                                               |
//|     Info        -   return code:                                 |
//|                     * -4, if internal EVD subroutine hasn't      |
//|                           converged                              |
//|                     * -2, if there is a point with class number  |
//|                           outside of [0..NClasses-1].            |
//|                     * -1, if incorrect parameters was passed     |
//|                           (NPoints<0, NVars<1, NClasses<2)       |
//|                     *  1, if task has been solved                |
//|                     *  2, if there was a multicollinearity in    |
//|                           training set, but task has been solved.|
//|     W           -   linear combination coefficients,             |
//|                     array[0..NVars-1]                            |
//+------------------------------------------------------------------+
static void CLDA::FisherLDA(CMatrixDouble &xy,const int npoints,
                            const int nvars,const int nclasses,
                            int &info,double &w[])
  {
//--- create a variable
   int i_=0;
//--- create matrix
   CMatrixDouble w2;
//--- initialization
   info=0;
//--- function call
   FisherLDAN(xy,npoints,nvars,nclasses,info,w2);
//--- check
   if(info>0)
     {
      //--- allocation
      ArrayResizeAL(w,nvars);
      //--- copy
      for(i_=0;i_<=nvars-1;i_++)
         w[i_]=w2[i_][0];
     }
  }
//+------------------------------------------------------------------+
//| N-dimensional multiclass Fisher LDA                              |
//| Subroutine finds coefficients of linear combinations which       |
//| optimally separates                                              |
//| training set on classes. It returns N-dimensional basis whose    |
//| vector are sorted                                                |
//| by quality of training set separation (in descending order).     |
//| INPUT PARAMETERS:                                                |
//|     XY          -   training set, array[0..NPoints-1,0..NVars].  |
//|                     First NVars columns store values of          |
//|                     independent variables, next column stores    |
//|                     number of class (from 0 to NClasses-1) which |
//|                     dataset element belongs to. Fractional values|
//|                     are rounded to nearest integer.              |
//|     NPoints     -   training set size, NPoints>=0                |
//|     NVars       -   number of independent variables, NVars>=1    |
//|     NClasses    -   number of classes, NClasses>=2               |
//| OUTPUT PARAMETERS:                                               |
//|     Info        -   return code:                                 |
//|                     * -4, if internal EVD subroutine hasn't      |
//|                           converged                              |
//|                     * -2, if there is a point with class number  |
//|                           outside of [0..NClasses-1].            |
//|                     * -1, if incorrect parameters was passed     |
//|                           (NPoints<0, NVars<1, NClasses<2)       |
//|                     *  1, if task has been solved                |
//|                     *  2, if there was a multicollinearity in    |
//|                           training set, but task has been solved.|
//|     W           -   basis, array[0..NVars-1,0..NVars-1]          |
//|                     columns of matrix stores basis vectors,      |
//|                     sorted by quality of training set separation | 
//|                     (in descending order)                        |
//+------------------------------------------------------------------+
static void CLDA::FisherLDAN(CMatrixDouble &xy,const int npoints,const int nvars,
                             const int nclasses,int &info,CMatrixDouble &w)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    m=0;
   double v=0;
   int    i_=0;
//--- creating arrays
   int    c[];
   double mu[];
   int    nc[];
   double tf[];
   double d[];
   double d2[];
   double work[];
//--- create matrix
   CMatrixDouble muc;
   CMatrixDouble sw;
   CMatrixDouble st;
   CMatrixDouble z;
   CMatrixDouble z2;
   CMatrixDouble tm;
   CMatrixDouble sbroot;
   CMatrixDouble a;
   CMatrixDouble xyproj;
   CMatrixDouble wproj;
//--- initialization
   info=0;
//--- Test data
   if((npoints<0 || nvars<1) || nclasses<2)
     {
      info=-1;
      return;
     }
   for(i=0;i<=npoints-1;i++)
     {
      //--- check
      if((int)MathRound(xy[i][nvars])<0 || (int)MathRound(xy[i][nvars])>=nclasses)
        {
         info=-2;
         return;
        }
     }
//--- change value
   info=1;
//--- Special case: NPoints<=1
//--- Degenerate task.
   if(npoints<=1)
     {
      info=2;
      //--- allocation
      w.Resize(nvars,nvars);
      //--- initialization
      for(i=0;i<=nvars-1;i++)
        {
         for(j=0;j<=nvars-1;j++)
           {
            //--- check
            if(i==j)
               w[i].Set(j,1);
            else
               w[i].Set(j,0);
           }
        }
      //--- exit the function
      return;
     }
//--- Prepare temporaries
   ArrayResizeAL(tf,nvars);
   ArrayResizeAL(work,MathMax(nvars,npoints)+1);
//--- Convert class labels from reals to integers (just for convenience)
   ArrayResizeAL(c,npoints);
   for(i=0;i<=npoints-1;i++)
      c[i]=(int)MathRound(xy[i][nvars]);
//--- Calculate class sizes and means
   ArrayResizeAL(mu,nvars);
   muc.Resize(nclasses,nvars);
   ArrayResizeAL(nc,nclasses);
   for(j=0;j<=nvars-1;j++)
      mu[j]=0;
   for(i=0;i<=nclasses-1;i++)
     {
      nc[i]=0;
      for(j=0;j<=nvars-1;j++)
         muc[i].Set(j,0);
     }
//--- calculation
   for(i=0;i<=npoints-1;i++)
     {
      for(i_=0;i_<=nvars-1;i_++)
         mu[i_]=mu[i_]+xy[i][i_];
      for(i_=0;i_<=nvars-1;i_++)
         muc[c[i]].Set(i_,muc[c[i]][i_]+xy[i][i_]);
      nc[c[i]]=nc[c[i]]+1;
     }
   for(i=0;i<=nclasses-1;i++)
     {
      v=1.0/(double)nc[i];
      for(i_=0;i_<=nvars-1;i_++)
         muc[i].Set(i_,v*muc[i][i_]);
     }
//--- change values
   v=1.0/(double)npoints;
   for(i_=0;i_<=nvars-1;i_++)
      mu[i_]=v*mu[i_];
//--- Create ST matrix
   st.Resize(nvars,nvars);
   for(i=0;i<=nvars-1;i++)
     {
      for(j=0;j<=nvars-1;j++)
         st[i].Set(j,0);
     }
//--- calculation
   for(k=0;k<=npoints-1;k++)
     {
      for(i_=0;i_<=nvars-1;i_++)
         tf[i_]=xy[k][i_];
      for(i_=0;i_<=nvars-1;i_++)
         tf[i_]=tf[i_]-mu[i_];
      for(i=0;i<=nvars-1;i++)
        {
         v=tf[i];
         for(i_=0;i_<=nvars-1;i_++)
            st[i].Set(i_,st[i][i_]+v*tf[i_]);
        }
     }
//--- Create SW matrix
   sw.Resize(nvars,nvars);
   for(i=0;i<=nvars-1;i++)
     {
      for(j=0;j<=nvars-1;j++)
         sw[i].Set(j,0);
     }
//--- calculation
   for(k=0;k<=npoints-1;k++)
     {
      for(i_=0;i_<=nvars-1;i_++)
         tf[i_]=xy[k][i_];
      for(i_=0;i_<=nvars-1;i_++)
         tf[i_]=tf[i_]-muc[c[k]][i_];
      for(i=0;i<=nvars-1;i++)
        {
         v=tf[i];
         for(i_=0;i_<=nvars-1;i_++)
            sw[i].Set(i_,sw[i][i_]+v*tf[i_]);
        }
     }
//--- Maximize ratio J=(w'*ST*w)/(w'*SW*w).
//--- First,make transition from w to v such that w'*ST*w becomes v'*v:
//---    v=root(ST)*w=R*w
//---    R=root(D)*Z'
//---    w=(root(ST)^-1)*v=RI*v
//---    RI=Z*inv(root(D))
//---    J=(v'*v)/(v'*(RI'*SW*RI)*v)
//---    ST=Z*D*Z'
//---    so we have
//---    J=(v'*v) / (v'*(inv(root(D))*Z'*SW*Z*inv(root(D)))*v)=
//=(v'*v) / (v'*A*v)
   if(!CEigenVDetect::SMatrixEVD(st,nvars,1,true,d,z))
     {
      info=-4;
      return;
     }
//--- allocation
   w.Resize(nvars,nvars);
//--- check
   if(d[nvars-1]<=0.0 || d[0]<=1000*CMath::m_machineepsilon*d[nvars-1])
     {
      //--- Special case: D[NVars-1]<=0
      //--- Degenerate task (all variables takes the same value).
      if(d[nvars-1]<=0.0)
        {
         info=2;
         for(i=0;i<=nvars-1;i++)
           {
            for(j=0;j<=nvars-1;j++)
              {
               //--- check
               if(i==j)
                  w[i].Set(j,1);
               else
                  w[i].Set(j,0);
              }
           }
         //--- exit the function
         return;
        }
      //--- Special case: degenerate ST matrix,multicollinearity found.
      //--- Since we know ST eigenvalues/vectors we can translate task to
      //--- non-degenerate form.
      //--- Let WG is orthogonal basis of the non zero variance subspace
      //--- of the ST and let WZ is orthogonal basis of the zero variance
      //--- subspace.
      //--- Projection on WG allows us to use LDA on reduced M-dimensional
      //--- subspace,N-M vectors of WZ allows us to update reduced LDA
      //--- factors to full N-dimensional subspace.
      m=0;
      for(k=0;k<=nvars-1;k++)
        {
         //--- check
         if(d[k]<=1000*CMath::m_machineepsilon*d[nvars-1])
            m=k+1;
        }
      //--- check
      if(!CAp::Assert(m!=0,__FUNCTION__+": internal error #1"))
         return;
      //--- allocation
      xyproj.Resize(npoints,nvars-m+1);
      //--- function call
      CBlas::MatrixMatrixMultiply(xy,0,npoints-1,0,nvars-1,false,z,0,nvars-1,m,nvars-1,false,1.0,xyproj,0,npoints-1,0,nvars-m-1,0.0,work);
      for(i=0;i<=npoints-1;i++)
         xyproj[i].Set(nvars-m,xy[i][nvars]);
      //--- function call
      FisherLDAN(xyproj,npoints,nvars-m,nclasses,info,wproj);
      //--- check
      if(info<0)
         return;
      //--- function call
      CBlas::MatrixMatrixMultiply(z,0,nvars-1,m,nvars-1,false,wproj,0,nvars-m-1,0,nvars-m-1,false,1.0,w,0,nvars-1,0,nvars-m-1,0.0,work);
      //--- change values
      for(k=nvars-m;k<=nvars-1;k++)
        {
         for(i_=0;i_<=nvars-1;i_++)
            w[i_].Set(k,z[i_][k-nvars+m]);
        }
      info=2;
     }
   else
     {
      //--- General case: no multicollinearity
      tm.Resize(nvars,nvars);
      a.Resize(nvars,nvars);
      //--- function call
      CBlas::MatrixMatrixMultiply(sw,0,nvars-1,0,nvars-1,false,z,0,nvars-1,0,nvars-1,false,1.0,tm,0,nvars-1,0,nvars-1,0.0,work);
      CBlas::MatrixMatrixMultiply(z,0,nvars-1,0,nvars-1,true,tm,0,nvars-1,0,nvars-1,false,1.0,a,0,nvars-1,0,nvars-1,0.0,work);
      //--- change values
      for(i=0;i<=nvars-1;i++)
        {
         for(j=0;j<=nvars-1;j++)
            a[i].Set(j,a[i][j]/MathSqrt(d[i]*d[j]));
        }
      //--- check
      if(!CEigenVDetect::SMatrixEVD(a,nvars,1,true,d2,z2))
        {
         info=-4;
         return;
        }
      //--- calculation
      for(k=0;k<=nvars-1;k++)
        {
         for(i=0;i<=nvars-1;i++)
            tf[i]=z2[i][k]/MathSqrt(d[i]);
         for(i=0;i<=nvars-1;i++)
           {
            v=0.0;
            for(i_=0;i_<=nvars-1;i_++)
               v+=z[i][i_]*tf[i_];
            w[i].Set(k,v);
           }
        }
     }
//--- Post-processing:
//--- * normalization
//--- * converting to non-negative form,if possible
   for(k=0;k<=nvars-1;k++)
     {
      //--- calculation
      v=0.0;
      for(i_=0;i_<=nvars-1;i_++)
         v+=w[i_][k]*w[i_][k];
      v=1/MathSqrt(v);
      for(i_=0;i_<=nvars-1;i_++)
         w[i_].Set(k,v*w[i_][k]);
      v=0;
      for(i=0;i<=nvars-1;i++)
         v=v+w[i][k];
      //--- check
      if(v<0.0)
        {
         for(i_=0;i_<=nvars-1;i_++)
            w[i_].Set(k,-1*w[i_][k]);
        }
     }
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CLinReg                                      |
//+------------------------------------------------------------------+
class CLinearModel
  {
public:
   double            m_w[];
   //--- constructor, destructor
                     CLinearModel(void);
                    ~CLinearModel(void);
   //--- copy
   void              Copy(CLinearModel &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLinearModel::CLinearModel(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLinearModel::~CLinearModel(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CLinearModel::Copy(CLinearModel &obj)
  {
//--- copy array
   ArrayCopy(m_w,obj.m_w);
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CLinearModel                     |
//+------------------------------------------------------------------+
class CLinearModelShell
  {
private:
   CLinearModel      m_innerobj;
public:
   //--- constructors, destructor
                     CLinearModelShell(void);
                     CLinearModelShell(CLinearModel &obj);
                    ~CLinearModelShell(void);
   //--- method
   CLinearModel     *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLinearModelShell::CLinearModelShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CLinearModelShell::CLinearModelShell(CLinearModel &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLinearModelShell::~CLinearModelShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CLinearModel *CLinearModelShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| LRReport structure contains additional information about linear  |
//| model:                                                           |
//| * C             -   covariation matrix, array[0..NVars,0..NVars].|
//|                     C[i,j] = Cov(A[i],A[j])                      |
//| * RMSError      -   root mean square error on a training set     |
//| * AvgError      -   average error on a training set              |
//| * AvgRelError   -   average relative error on a training set     |
//|                     (excluding observations with zero function   |
//|                     value).                                      |
//| * CVRMSError    -   leave-one-out cross-validation estimate of   |
//|                     generalization error. Calculated using fast  |
//|                     algorithm with O(NVars*NPoints) complexity.  |
//| * CVAvgError    -   cross-validation estimate of average error   |
//| * CVAvgRelError -   cross-validation estimate of average relative| 
//|                     error                                        |
//| All other fields of the structure are intended for internal use  |
//| and should not be used outside ALGLIB.                           |
//+------------------------------------------------------------------+
class CLRReport
  {
public:
   //--- variables
   double            m_rmserror;
   double            m_avgerror;
   double            m_avgrelerror;
   double            m_cvrmserror;
   double            m_cvavgerror;
   double            m_cvavgrelerror;
   int               m_ncvdefects;
   //--- array
   int               m_cvdefects[];
   //--- matrix
   CMatrixDouble     m_c;
   //--- constructor, destructor
                     CLRReport(void);
                    ~CLRReport(void);
   //--- copy
   void              Copy(CLRReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLRReport::CLRReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLRReport::~CLRReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CLRReport::Copy(CLRReport &obj)
  {
//--- copy variables
   m_rmserror=obj.m_rmserror;
   m_avgerror=obj.m_avgerror;
   m_avgrelerror=obj.m_avgrelerror;
   m_cvrmserror=obj.m_cvrmserror;
   m_cvavgerror=obj.m_cvavgerror;
   m_cvavgrelerror=obj.m_cvavgrelerror;
   m_ncvdefects=obj.m_ncvdefects;
//--- copy array
   ArrayCopy(m_cvdefects,obj.m_cvdefects);
//--- copy matrix
   m_c=obj.m_c;
  }
//+------------------------------------------------------------------+
//| LRReport structure contains additional information about linear  | 
//| model:                                                           |
//| * C             -   covariation matrix, array[0..NVars,0..NVars].|
//|                     C[i,j]=Cov(A[i],A[j])                        |
//| * RMSError      -   root mean square error on a training set     |
//| * AvgError      -   average error on a training set              |
//| * AvgRelError   -   average relative error on a training set     |
//|                     (excluding observations with zero function   |
//|                     value).                                      |
//| * CVRMSError    -   leave-one-out cross-validation estimate of   |
//|                     generalization error. Calculated using fast  |
//|                     algorithm with O(NVars*NPoints) complexity.  |
//| * CVAvgError    -   cross-validation estimate of average error   |
//| * CVAvgRelError -   cross-validation estimate of average relative| 
//|                     error                                        |
//| All other fields of the structure are intended for internal use  | 
//| and should not be used outside ALGLIB.                           |
//+------------------------------------------------------------------+
class CLRReportShell
  {
private:
   CLRReport         m_innerobj;
public:
   //--- constructors, destructor
                     CLRReportShell(void);
                     CLRReportShell(CLRReport &obj);
                    ~CLRReportShell(void);
   //--- methods
   double            GetRMSError(void);
   void              SetRMSError(const double d);
   double            GetAvgError(void);
   void              SetAvgError(const double d);
   double            GetAvgRelError(void);
   void              SetAvgRelError(const double d);
   double            GetCVRMSError(void);
   void              SetCVRMSError(const double d);
   double            GetCVAvgError(void);
   void              SetCVAvgError(const double d);
   double            GetCVAvgRelError(void);
   void              SetCVAvgRelError(const double d);
   int               GetNCVDEfects(void);
   void              SetNCVDEfects(const int i);
   CLRReport        *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLRReportShell::CLRReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CLRReportShell::CLRReportShell(CLRReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLRReportShell::~CLRReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable rmserror                       |
//+------------------------------------------------------------------+
double CLRReportShell::GetRMSError(void)
  {
//--- return result
   return(m_innerobj.m_rmserror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable rmserror                      |
//+------------------------------------------------------------------+
void CLRReportShell::SetRMSError(const double d)
  {
//--- change value
   m_innerobj.m_rmserror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgerror                       |
//+------------------------------------------------------------------+
double CLRReportShell::GetAvgError(void)
  {
//--- return result
   return(m_innerobj.m_avgerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgerror                      |
//+------------------------------------------------------------------+
void CLRReportShell::SetAvgError(const double d)
  {
//--- change value
   m_innerobj.m_avgerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgrelerror                    |
//+------------------------------------------------------------------+
double CLRReportShell::GetAvgRelError(void)
  {
//--- return result
   return(m_innerobj.m_avgrelerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgrelerror                   |
//+------------------------------------------------------------------+
void CLRReportShell::SetAvgRelError(const double d)
  {
//--- change value
   m_innerobj.m_avgrelerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable cvrmserror                     |
//+------------------------------------------------------------------+
double CLRReportShell::GetCVRMSError(void)
  {
//--- return result
   return(m_innerobj.m_cvrmserror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable cvrmserror                    |
//+------------------------------------------------------------------+
void CLRReportShell::SetCVRMSError(const double d)
  {
//--- change value
   m_innerobj.m_cvrmserror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable cvavgerror                     |
//+------------------------------------------------------------------+
double CLRReportShell::GetCVAvgError(void)
  {
//--- return result
   return(m_innerobj.m_cvavgerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable cvavgerror                    |
//+------------------------------------------------------------------+
void CLRReportShell::SetCVAvgError(const double d)
  {
//--- change value
   m_innerobj.m_cvavgerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable cvavgrelerror                  |
//+------------------------------------------------------------------+
double CLRReportShell::GetCVAvgRelError(void)
  {
//--- return result
   return(m_innerobj.m_cvavgrelerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable cvavgrelerror                 |
//+------------------------------------------------------------------+
void CLRReportShell::SetCVAvgRelError(const double d)
  {
//--- change value
   m_innerobj.m_cvavgrelerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable ncvdefects                     |
//+------------------------------------------------------------------+
int CLRReportShell::GetNCVDEfects(void)
  {
//--- return result
   return(m_innerobj.m_ncvdefects);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable ncvdefects                    |
//+------------------------------------------------------------------+
void CLRReportShell::SetNCVDEfects(const int i)
  {
//--- change value
   m_innerobj.m_ncvdefects=i;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CLRReport *CLRReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Linear regression class                                          |
//+------------------------------------------------------------------+
class CLinReg
  {
private:
   //--- private method
   static void       LRInternal(CMatrixDouble &xy,double &s[],const int npoints,const int nvars,int &info,CLinearModel &lm,CLRReport &ar);
public:
   //--- constant
   static const int  m_lrvnum;
   //--- constructor, destructor
                     CLinReg(void);
                    ~CLinReg(void);
   //--- public methods
   static void       LRBuild(CMatrixDouble &xy,const int npoints,const int nvars,int &info,CLinearModel &lm,CLRReport &ar);
   static void       LRBuildS(CMatrixDouble &xy,double &s[],const int npoints,const int nvars,int &info,CLinearModel &lm,CLRReport &ar);
   static void       LRBuildZS(CMatrixDouble &xy,double &s[],const int npoints,const int nvars,int &info,CLinearModel &lm,CLRReport &ar);
   static void       LRBuildZ(CMatrixDouble &xy,const int npoints,const int nvars,int &info,CLinearModel &lm,CLRReport &ar);
   static void       LRUnpack(CLinearModel &lm,double &v[],int &nvars);
   static void       LRPack(double &v[],const int nvars,CLinearModel &lm);
   static double     LRProcess(CLinearModel &lm,double &x[]);
   static double     LRRMSError(CLinearModel &lm,CMatrixDouble &xy,const int npoints);
   static double     LRAvgError(CLinearModel &lm,CMatrixDouble &xy,const int npoints);
   static double     LRAvgRelError(CLinearModel &lm,CMatrixDouble &xy,const int npoints);
   static void       LRCopy(CLinearModel &lm1,CLinearModel &lm2);
   static void       LRLines(CMatrixDouble &xy,double &s[],const int n,int &info,double &a,double &b,double &vara,double &varb,double &covab,double &corrab,double &p);
   static void       LRLine(CMatrixDouble &xy,const int n,int &info,double &a,double &b);
  };
//+------------------------------------------------------------------+
//| Initialize constant                                              |
//+------------------------------------------------------------------+
const int CLinReg::m_lrvnum=5;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLinReg::CLinReg(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLinReg::~CLinReg(void)
  {

  }
//+------------------------------------------------------------------+
//| Linear regression                                                |
//| Subroutine builds model:                                         |
//|     Y = A(0)*X[0] + ... + A(N-1)*X[N-1] + A(N)                   |
//| and model found in ALGLIB format, covariation matrix, training   | 
//| set errors (rms, average, average relative) and leave-one-out    |
//| cross-validation estimate of the generalization error. CV        |
//| estimate calculated using fast algorithm with O(NPoints*NVars)   |
//| complexity.                                                      |
//| When  covariation  matrix  is  calculated  standard deviations of| 
//| function values are assumed to be equal to RMS error on the      |
//| training set.                                                    |
//| INPUT PARAMETERS:                                                |
//|     XY          -   training set, array [0..NPoints-1,0..NVars]: |
//|                     * NVars columns - independent variables      |
//|                     * last column - dependent variable           |
//|     NPoints     -   training set size, NPoints>NVars+1           |
//|     NVars       -   number of independent variables              |
//| OUTPUT PARAMETERS:                                               |
//|     Info        -   return code:                                 |
//|                     * -255, in case of unknown internal error    |
//|                     * -4, if internal SVD subroutine haven't     |
//|                           converged                              |
//|                     * -1, if incorrect parameters was passed     |
//|                           (NPoints<NVars+2, NVars<1).            |
//|                     *  1, if subroutine successfully finished    |
//|     LM          -   linear model in the ALGLIB format. Use       |
//|                     subroutines of this unit to work with the    |
//|                     model.                                       |
//|     AR          -   additional results                           |
//+------------------------------------------------------------------+
static void CLinReg::LRBuild(CMatrixDouble &xy,const int npoints,const int nvars,
                             int &info,CLinearModel &lm,CLRReport &ar)
  {
//--- create variables
   int    i=0;
   double sigma2=0;
   int    i_=0;
//--- create array
   double s[];
//--- initialization
   info=0;
//--- check
   if(npoints<=nvars+1 || nvars<1)
     {
      info=-1;
      return;
     }
//--- allocation
   ArrayResizeAL(s,npoints);
   for(i=0;i<=npoints-1;i++)
      s[i]=1;
//--- function call
   LRBuildS(xy,s,npoints,nvars,info,lm,ar);
//--- check
   if(info<0)
      return;
//--- calculation
   sigma2=CMath::Sqr(ar.m_rmserror)*npoints/(npoints-nvars-1);
   for(i=0;i<=nvars;i++)
     {
      for(i_=0;i_<=nvars;i_++)
         ar.m_c[i].Set(i_,sigma2*ar.m_c[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| Linear regression                                                |
//| Variant of LRBuild which uses vector of standatd deviations      |
//| (errors in function values).                                     |
//| INPUT PARAMETERS:                                                |
//|     XY          -   training set, array [0..NPoints-1,0..NVars]: |
//|                     * NVars columns - independent variables      |
//|                     * last column - dependent variable           |
//|     S           -   standard deviations (errors in function      |
//|                     values) array[0..NPoints-1], S[i]>0.         |
//|     NPoints     -   training set size, NPoints>NVars+1           |
//|     NVars       -   number of independent variables              |
//| OUTPUT PARAMETERS:                                               |
//|     Info        -   return code:                                 |
//|                     * -255, in case of unknown internal error    |
//|                     * -4, if internal SVD subroutine haven't     |
//|                     converged                                    |
//|                     * -1, if incorrect parameters was passed     |
//|                     (NPoints<NVars+2, NVars<1).                  |
//|                     * -2, if S[I]<=0                             |
//|                     *  1, if subroutine successfully finished    |
//|     LM          -   linear model in the ALGLIB format. Use       |
//|                     subroutines of this unit to work with the    |
//|                     model.                                       |
//|     AR          -   additional results                           |
//+------------------------------------------------------------------+
static void CLinReg::LRBuildS(CMatrixDouble &xy,double &s[],const int npoints,
                              const int nvars,int &info,CLinearModel &lm,
                              CLRReport &ar)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   int    offs=0;
   double mean=0;
   double variance=0;
   double skewness=0;
   double kurtosis=0;
   int    i_=0;
//--- creating arrays
   double x[];
   double means[];
   double sigmas[];
//--- create array
   CMatrixDouble xyi;
//--- initialization
   info=0;
//--- Test parameters
   if(npoints<=nvars+1 || nvars<1)
     {
      info=-1;
      return;
     }
//--- Copy data,add one more column (constant term)
   xyi.Resize(npoints,nvars+2);
   for(i=0;i<=npoints-1;i++)
     {
      for(i_=0;i_<=nvars-1;i_++)
         xyi[i].Set(i_,xy[i][i_]);
      xyi[i].Set(nvars,1);
      xyi[i].Set(nvars+1,xy[i][nvars]);
     }
//--- Standartization
   ArrayResizeAL(x,npoints);
   ArrayResizeAL(means,nvars);
   ArrayResizeAL(sigmas,nvars);
   for(j=0;j<=nvars-1;j++)
     {
      //--- copy
      for(i_=0;i_<=npoints-1;i_++)
         x[i_]=xy[i_][j];
      //--- function call
      CBaseStat::SampleMoments(x,npoints,mean,variance,skewness,kurtosis);
      //--- change values
      means[j]=mean;
      sigmas[j]=MathSqrt(variance);
      //--- check
      if(sigmas[j]==0.0)
         sigmas[j]=1;
      //--- calculation
      for(i=0;i<=npoints-1;i++)
         xyi[i].Set(j,(xyi[i][j]-means[j])/sigmas[j]);
     }
//--- Internal processing
   LRInternal(xyi,s,npoints,nvars+1,info,lm,ar);
//--- check
   if(info<0)
      return;
//--- Un-standartization
   offs=(int)MathRound(lm.m_w[3]);
   for(j=0;j<=nvars-1;j++)
     {
      //--- Constant term is updated (and its covariance too,
      //--- since it gets some variance from J-th component)
      lm.m_w[offs+nvars]=lm.m_w[offs+nvars]-lm.m_w[offs+j]*means[j]/sigmas[j];
      v=means[j]/sigmas[j];
      for(i_=0;i_<=nvars;i_++)
         ar.m_c[nvars].Set(i_,ar.m_c[nvars][i_]-v*ar.m_c[j][i_]);
      for(i_=0;i_<=nvars;i_++)
         ar.m_c[i_].Set(nvars,ar.m_c[i_][nvars]-v*ar.m_c[i_][j]);
      //--- J-th term is updated
      lm.m_w[offs+j]=lm.m_w[offs+j]/sigmas[j];
      v=1/sigmas[j];
      for(i_=0;i_<=nvars;i_++)
         ar.m_c[j].Set(i_,v*ar.m_c[j][i_]);
      for(i_=0;i_<=nvars;i_++)
         ar.m_c[i_].Set(j,v*ar.m_c[i_][j]);
     }
  }
//+------------------------------------------------------------------+
//| Like LRBuildS, but builds model                                  |
//|     Y=A(0)*X[0] + ... + A(N-1)*X[N-1]                            |
//| i.m_e. with zero constant term.                                  |
//+------------------------------------------------------------------+
static void CLinReg::LRBuildZS(CMatrixDouble &xy,double &s[],const int npoints,
                               const int nvars,int &info,CLinearModel &lm,
                               CLRReport &ar)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double v=0;
   int    offs=0;
   double mean=0;
   double variance=0;
   double skewness=0;
   double kurtosis=0;
   int    i_=0;
//--- creating arrays
   double x[];
   double c[];
//--- create matrix
   CMatrixDouble xyi;
//--- initialization
   info=0;
//--- Test parameters
   if(npoints<=nvars+1 || nvars<1)
     {
      info=-1;
      return;
     }
//--- Copy data,add one more column (constant term)
   xyi.Resize(npoints,nvars+2);
   for(i=0;i<=npoints-1;i++)
     {
      for(i_=0;i_<=nvars-1;i_++)
         xyi[i].Set(i_,xy[i][i_]);
      xyi[i].Set(nvars,0);
      xyi[i].Set(nvars+1,xy[i][nvars]);
     }
//--- Standartization: unusual scaling
   ArrayResizeAL(x,npoints);
   ArrayResizeAL(c,nvars);
   for(j=0;j<=nvars-1;j++)
     {
      for(i_=0;i_<=npoints-1;i_++)
         x[i_]=xy[i_][j];
      //--- function call
      CBaseStat::SampleMoments(x,npoints,mean,variance,skewness,kurtosis);
      //--- check
      if(MathAbs(mean)>MathSqrt(variance))
        {
         //--- variation is relatively small,it is better to
         //--- bring mean value to 1
         c[j]=mean;
        }
      else
        {
         //--- variation is large,it is better to bring variance to 1
         if(variance==0.0)
            variance=1;
         c[j]=MathSqrt(variance);
        }
      for(i=0;i<=npoints-1;i++)
         xyi[i].Set(j,xyi[i][j]/c[j]);
     }
//--- Internal processing
   LRInternal(xyi,s,npoints,nvars+1,info,lm,ar);
//--- check
   if(info<0)
      return;
//--- Un-standartization
   offs=(int)MathRound(lm.m_w[3]);
   for(j=0;j<=nvars-1;j++)
     {
      //--- J-th term is updated
      lm.m_w[offs+j]=lm.m_w[offs+j]/c[j];
      v=1/c[j];
      for(i_=0;i_<=nvars;i_++)
         ar.m_c[j].Set(i_,v*ar.m_c[j][i_]);
      for(i_=0;i_<=nvars;i_++)
         ar.m_c[i_].Set(j,v*ar.m_c[i_][j]);
     }
  }
//+------------------------------------------------------------------+
//| Like LRBuild but builds model                                    |
//|     Y=A(0)*X[0] + ... + A(N-1)*X[N-1]                            |
//| i.m_e. with zero constant term.                                  |
//+------------------------------------------------------------------+
static void CLinReg::LRBuildZ(CMatrixDouble &xy,const int npoints,const int nvars,
                              int &info,CLinearModel &lm,CLRReport &ar)
  {
//--- create variables
   int i=0;
   double sigma2=0;
   int i_=0;
//--- create array
   double s[];
//--- initialization
   info=0;
//--- check
   if(npoints<=nvars+1 || nvars<1)
     {
      info=-1;
      return;
     }
//--- allocation
   ArrayResizeAL(s,npoints);
   for(i=0;i<=npoints-1;i++)
      s[i]=1;
//--- function call
   LRBuildZS(xy,s,npoints,nvars,info,lm,ar);
//--- check
   if(info<0)
      return;
//--- calculation
   sigma2=CMath::Sqr(ar.m_rmserror)*npoints/(npoints-nvars-1);
   for(i=0;i<=nvars;i++)
     {
      for(i_=0;i_<=nvars;i_++)
         ar.m_c[i].Set(i_,sigma2*ar.m_c[i][i_]);
     }
  }
//+------------------------------------------------------------------+
//| Unpacks coefficients of linear model.                            |
//| INPUT PARAMETERS:                                                |
//|     LM          -   linear model in ALGLIB format                |
//| OUTPUT PARAMETERS:                                               |
//|     V           -   coefficients,array[0..NVars]                 |
//|                     constant term (intercept) is stored in the   |
//|                     V[NVars].                                    |
//|     NVars       -   number of independent variables (one less    |
//|                     than number of coefficients)                 |
//+------------------------------------------------------------------+
static void CLinReg::LRUnpack(CLinearModel &lm,double &v[],int &nvars)
  {
//--- create variables
   int offs=0;
   int i_=0;
   int i1_=0;
//--- initialization
   nvars=0;
//--- check
   if(!CAp::Assert((int)MathRound(lm.m_w[1])==m_lrvnum,__FUNCTION__+": Incorrect LINREG version!"))
      return;
//--- change values
   nvars=(int)MathRound(lm.m_w[2]);
   offs=(int)MathRound(lm.m_w[3]);
//--- allocation
   ArrayResizeAL(v,nvars+1);
//--- calculation
   i1_=offs;
   for(i_=0;i_<=nvars;i_++)
      v[i_]=lm.m_w[i_+i1_];
  }
//+------------------------------------------------------------------+
//| "Packs" coefficients and creates linear model in ALGLIB format   |
//| (LRUnpack reversed).                                             |
//| INPUT PARAMETERS:                                                |
//|     V           -   coefficients, array[0..NVars]                |
//|     NVars       -   number of independent variables              |
//| OUTPUT PAREMETERS:                                               |
//|     LM          -   linear model.                                |
//+------------------------------------------------------------------+
static void CLinReg::LRPack(double &v[],const int nvars,CLinearModel &lm)
  {
//--- create variables
   int offs=0;
   int i_=0;
   int i1_=0;
//--- allocation
   ArrayResizeAL(lm.m_w,5+nvars);
//--- change values
   offs=4;
   lm.m_w[0]=4+nvars+1;
   lm.m_w[1]=m_lrvnum;
   lm.m_w[2]=nvars;
   lm.m_w[3]=offs;
//--- calculation
   i1_=-offs;
   for(i_=offs;i_<=offs+nvars;i_++)
      lm.m_w[i_]=v[i_+i1_];
  }
//+------------------------------------------------------------------+
//| Procesing                                                        |
//| INPUT PARAMETERS:                                                |
//|     LM      -   linear model                                     |
//|     X       -   input vector, array[0..NVars-1].                 |
//| Result:                                                          |
//|     value of linear model regression estimate                    |
//+------------------------------------------------------------------+
static double CLinReg::LRProcess(CLinearModel &lm,double &x[])
  {
//--- create variables
   double v=0;
   int    offs=0;
   int    nvars=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(!CAp::Assert((int)MathRound(lm.m_w[1])==m_lrvnum,__FUNCTION__+": Incorrect LINREG version!"))
      return(EMPTY_VALUE);
//--- change values
   nvars=(int)MathRound(lm.m_w[2]);
   offs=(int)MathRound(lm.m_w[3]);
   i1_=offs;
   v=0.0;
//--- calculation
   for(i_=0;i_<=nvars-1;i_++)
      v+=x[i_]*lm.m_w[i_+i1_];
//--- return result
   return(v+lm.m_w[offs+nvars]);
  }
//+------------------------------------------------------------------+
//| RMS error on the test set                                        |
//| INPUT PARAMETERS:                                                |
//|     LM      -   linear model                                     |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     root mean square error.                                      |
//+------------------------------------------------------------------+
static double CLinReg::LRRMSError(CLinearModel &lm,CMatrixDouble &xy,
                                  const int npoints)
  {
//--- create variables
   double result=0;
   int    i=0;
   double v=0;
   int    offs=0;
   int    nvars=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(!CAp::Assert((int)MathRound(lm.m_w[1])==m_lrvnum,__FUNCTION__+": Incorrect LINREG version!"))
      return(EMPTY_VALUE);
//--- change values
   nvars=(int)MathRound(lm.m_w[2]);
   offs=(int)MathRound(lm.m_w[3]);
   result=0;
//--- calculation
   for(i=0;i<=npoints-1;i++)
     {
      i1_=offs;
      v=0.0;
      for(i_=0;i_<=nvars-1;i_++)
         v+=xy[i][i_]*lm.m_w[i_+i1_];
      v=v+lm.m_w[offs+nvars];
      result=result+CMath::Sqr(v-xy[i][nvars]);
     }
//--- return result
   return(MathSqrt(result/npoints));
  }
//+------------------------------------------------------------------+
//| Average error on the test set                                    |
//| INPUT PARAMETERS:                                                |
//|     LM      -   linear model                                     |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     average error.                                               |
//+------------------------------------------------------------------+
static double CLinReg::LRAvgError(CLinearModel &lm,CMatrixDouble &xy,
                                  const int npoints)
  {
//--- create variables
   double result=0;
   int    i=0;
   double v=0;
   int    offs=0;
   int    nvars=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(!CAp::Assert((int)MathRound(lm.m_w[1])==m_lrvnum,__FUNCTION__+": Incorrect LINREG version!"))
      return(EMPTY_VALUE);
//--- initialization
   nvars=(int)MathRound(lm.m_w[2]);
   offs=(int)MathRound(lm.m_w[3]);
   result=0;
//--- calculation
   for(i=0;i<=npoints-1;i++)
     {
      i1_=offs;
      v=0.0;
      for(i_=0;i_<=nvars-1;i_++)
         v+=xy[i][i_]*lm.m_w[i_+i1_];
      v=v+lm.m_w[offs+nvars];
      result=result+MathAbs(v-xy[i][nvars]);
     }
//--- return result
   return(result/npoints);
  }
//+------------------------------------------------------------------+
//| RMS error on the test set                                        |
//| INPUT PARAMETERS:                                                |
//|     LM      -   linear model                                     |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     average relative error.                                      |
//+------------------------------------------------------------------+
static double CLinReg::LRAvgRelError(CLinearModel &lm,CMatrixDouble &xy,
                                     const int npoints)
  {
//--- create variables
   double result=0;
   int    i=0;
   int    k=0;
   double v=0;
   int    offs=0;
   int    nvars=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(!CAp::Assert((int)MathRound(lm.m_w[1])==m_lrvnum,__FUNCTION__+": Incorrect LINREG version!"))
      return(EMPTY_VALUE);
//--- initialization
   nvars=(int)MathRound(lm.m_w[2]);
   offs=(int)MathRound(lm.m_w[3]);
   result=0;
   k=0;
//--- calculation
   for(i=0;i<=npoints-1;i++)
     {
      //--- check
      if(xy[i][nvars]!=0.0)
        {
         i1_=offs;
         v=0.0;
         for(i_=0;i_<=nvars-1;i_++)
            v+=xy[i][i_]*lm.m_w[i_+i1_];
         v=v+lm.m_w[offs+nvars];
         //--- get result
         result=result+MathAbs((v-xy[i][nvars])/xy[i][nvars]);
         k=k+1;
        }
     }
//--- check
   if(k!=0)
      result=result/k;
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Copying of LinearModel strucure                                  |
//| INPUT PARAMETERS:                                                |
//|     LM1 -   original                                             |
//| OUTPUT PARAMETERS:                                               |
//|     LM2 -   copy                                                 |
//+------------------------------------------------------------------+
static void CLinReg::LRCopy(CLinearModel &lm1,CLinearModel &lm2)
  {
//--- create variables
   int k=0;
   int i_=0;
//--- initialization
   k=(int)MathRound(lm1.m_w[0]);
//--- allocation
   ArrayResizeAL(lm2.m_w,k);
//--- copy
   for(i_=0;i_<=k-1;i_++)
      lm2.m_w[i_]=lm1.m_w[i_];
  }
//+------------------------------------------------------------------+
//| Class method                                                     |
//+------------------------------------------------------------------+
static void CLinReg::LRLines(CMatrixDouble &xy,double &s[],const int n,
                             int &info,double &a,double &b,double &vara,
                             double &varb,double &covab,double &corrab,
                             double &p)
  {
//--- create variables
   int    i=0;
   double ss=0;
   double sx=0;
   double sxx=0;
   double sy=0;
   double stt=0;
   double e1=0;
   double e2=0;
   double t=0;
   double chi2=0;
//--- initialization
   info=0;
   a=0;
   b=0;
   vara=0;
   varb=0;
   covab=0;
   corrab=0;
   p=0;
//--- check
   if(n<2)
     {
      info=-1;
      return;
     }
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if((double)(s[i])<=0.0)
        {
         info=-2;
         return;
        }
     }
//--- change value
   info=1;
//--- Calculate S,SX,SY,SXX
   ss=0;
   sx=0;
   sy=0;
   sxx=0;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      t=CMath::Sqr(s[i]);
      ss=ss+1/t;
      sx=sx+xy[i][0]/t;
      sy=sy+xy[i][1]/t;
      sxx=sxx+CMath::Sqr(xy[i][0])/t;
     }
//--- Test for condition number
   t=MathSqrt(4*CMath::Sqr(sx)+CMath::Sqr(ss-sxx));
   e1=0.5*(ss+sxx+t);
   e2=0.5*(ss+sxx-t);
//--- check
   if(MathMin(e1,e2)<=1000*CMath::m_machineepsilon*MathMax(e1,e2))
     {
      info=-3;
      return;
     }
//--- Calculate A,B
   a=0;
   b=0;
   stt=0;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      t=(xy[i][0]-sx/ss)/s[i];
      b=b+t*xy[i][1]/s[i];
      stt=stt+CMath::Sqr(t);
     }
   b=b/stt;
   a=(sy-sx*b)/ss;
//--- Calculate goodness-of-fit
   if(n>2)
     {
      chi2=0;
      for(i=0;i<=n-1;i++)
         chi2=chi2+CMath::Sqr((xy[i][1]-a-b*xy[i][0])/s[i]);
      //--- function call
      p=CIncGammaF::IncompleteGammaC((double)(n-2)/(double)2,chi2/2);
     }
   else
      p=1;
//--- Calculate other parameters
   vara=(1+CMath::Sqr(sx)/(ss*stt))/ss;
   varb=1/stt;
   covab=-(sx/(ss*stt));
   corrab=covab/MathSqrt(vara*varb);
  }
//+------------------------------------------------------------------+
//| Class method                                                     |
//+------------------------------------------------------------------+
static void CLinReg::LRLine(CMatrixDouble &xy,const int n,int &info,
                            double &a,double &b)
  {
//--- create variables
   int    i=0;
   double vara=0;
   double varb=0;
   double covab=0;
   double corrab=0;
   double p=0;
//--- create array
   double s[];
//--- initialization
   info=0;
   a=0;
   b=0;
//--- check
   if(n<2)
     {
      info=-1;
      return;
     }
//--- allocation
   ArrayResizeAL(s,n);
   for(i=0;i<=n-1;i++)
      s[i]=1;
//--- function call
   LRLines(xy,s,n,info,a,b,vara,varb,covab,corrab,p);
  }
//+------------------------------------------------------------------+
//| Internal linear regression subroutine                            |
//+------------------------------------------------------------------+
static void CLinReg::LRInternal(CMatrixDouble &xy,double &s[],const int npoints,
                                const int nvars,int &info,CLinearModel &lm,
                                CLRReport &ar)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    ncv=0;
   int    na=0;
   int    nacv=0;
   double r=0;
   double p=0;
   double epstol=0;
   int    offs=0;
   int    i_=0;
   int    i1_=0;
//--- creating arrays
   double b[];
   double sv[];
   double t[];
   double svi[];
   double work[];
//--- create matrix
   CMatrixDouble a;
   CMatrixDouble u;
   CMatrixDouble vt;
   CMatrixDouble vm;
   CMatrixDouble xym;
//--- create objects of classes
   CLRReport    ar2;
   CLinearModel tlm;
//--- initialization
   info=0;
   epstol=1000;
//--- Check for errors in data
   if(npoints<nvars || nvars<1)
     {
      info=-1;
      return;
     }
   for(i=0;i<=npoints-1;i++)
     {
      //--- check
      if(s[i]<=0.0)
        {
         info=-2;
         return;
        }
     }
//--- change value
   info=1;
//--- Create design matrix
   a.Resize(npoints,nvars);
   ArrayResizeAL(b,npoints);
//--- calculation
   for(i=0;i<=npoints-1;i++)
     {
      r=1/s[i];
      for(i_=0;i_<=nvars-1;i_++)
         a[i].Set(i_,r*xy[i][i_]);
      b[i]=xy[i][nvars]/s[i];
     }
//--- Allocate W:
//--- W[0]     array size
//--- W[1]     version number,0
//--- W[2]     NVars (minus 1,to be compatible with external representation)
//--- W[3]     coefficients offset
   ArrayResizeAL(lm.m_w,4+nvars);
   offs=4;
   lm.m_w[0]=4+nvars;
   lm.m_w[1]=m_lrvnum;
   lm.m_w[2]=nvars-1;
   lm.m_w[3]=offs;
//--- Solve problem using SVD:
//--- 0. check for degeneracy (different types)
//--- 1. A = U*diag(sv)*V'
//--- 2. T = b'*U
//--- 3. w = SUM((T[i]/sv[i])*V[..,i])
//--- 4. cov(wi,wj) = SUM(Vji*Vjk/sv[i]^2,K=1..M)
//--- see $15.4 of "Numerical Recipes in C" for more information
   ArrayResizeAL(t,nvars);
   ArrayResizeAL(svi,nvars);
   ar.m_c.Resize(nvars,nvars);
   vm.Resize(nvars,nvars);
//--- check
   if(!CSingValueDecompose::RMatrixSVD(a,npoints,nvars,1,1,2,sv,u,vt))
     {
      info=-4;
      return;
     }
//--- check
   if(sv[0]<=0.0)
     {
      //--- Degenerate case: zero design matrix.
      for(i=offs;i<=offs+nvars-1;i++)
         lm.m_w[i]=0;
      //--- change values
      ar.m_rmserror=LRRMSError(lm,xy,npoints);
      ar.m_avgerror=LRAvgError(lm,xy,npoints);
      ar.m_avgrelerror=LRAvgRelError(lm,xy,npoints);
      ar.m_cvrmserror=ar.m_rmserror;
      ar.m_cvavgerror=ar.m_avgerror;
      ar.m_cvavgrelerror=ar.m_avgrelerror;
      ar.m_ncvdefects=0;
      //--- allocation
      ArrayResizeAL(ar.m_cvdefects,nvars);
      ar.m_c.Resize(nvars,nvars);
      for(i=0;i<=nvars-1;i++)
        {
         for(j=0;j<=nvars-1;j++)
            ar.m_c[i].Set(j,0);
        }
      //--- exit the function
      return;
     }
//--- check
   if(sv[nvars-1]<=epstol*CMath::m_machineepsilon*sv[0])
     {
      //--- Degenerate case,non-zero design matrix.
      //--- We can leave it and solve task in SVD least squares fashion.
      //--- Solution and covariance matrix will be obtained correctly,
      //--- but CV error estimates - will not. It is better to reduce
      //--- it to non-degenerate task and to obtain correct CV estimates.
      for(k=nvars;k>=1;k--)
        {
         //--- check
         if(sv[k-1]>epstol*CMath::m_machineepsilon*sv[0])
           {
            //--- Reduce
            xym.Resize(npoints,k+1);
            for(i=0;i<=npoints-1;i++)
              {
               for(j=0;j<=k-1;j++)
                 {
                  //--- calculation
                  r=0.0;
                  for(i_=0;i_<=nvars-1;i_++)
                     r+=xy[i][i_]*vt[j][i_];
                  xym[i].Set(j,r);
                 }
               xym[i].Set(k,xy[i][nvars]);
              }
            //--- Solve
            LRInternal(xym,s,npoints,k,info,tlm,ar2);
            //--- check
            if(info!=1)
               return;
            //--- Convert back to un-reduced format
            for(j=0;j<=nvars-1;j++)
               lm.m_w[offs+j]=0;
            for(j=0;j<=k-1;j++)
              {
               r=tlm.m_w[offs+j];
               i1_=-offs;
               for(i_=offs;i_<=offs+nvars-1;i_++)
                  lm.m_w[i_]=lm.m_w[i_]+r*vt[j][i_+i1_];
              }
            //--- change values
            ar.m_rmserror=ar2.m_rmserror;
            ar.m_avgerror=ar2.m_avgerror;
            ar.m_avgrelerror=ar2.m_avgrelerror;
            ar.m_cvrmserror=ar2.m_cvrmserror;
            ar.m_cvavgerror=ar2.m_cvavgerror;
            ar.m_cvavgrelerror=ar2.m_cvavgrelerror;
            ar.m_ncvdefects=ar2.m_ncvdefects;
            //--- allocation
            ArrayResizeAL(ar.m_cvdefects,nvars);
            for(j=0;j<=ar.m_ncvdefects-1;j++)
               ar.m_cvdefects[j]=ar2.m_cvdefects[j];
            //--- allocation
            ar.m_c.Resize(nvars,nvars);
            ArrayResizeAL(work,nvars+1);
            //--- function calls
            CBlas::MatrixMatrixMultiply(ar2.m_c,0,k-1,0,k-1,false,vt,0,k-1,0,nvars-1,false,1.0,vm,0,k-1,0,nvars-1,0.0,work);
            CBlas::MatrixMatrixMultiply(vt,0,k-1,0,nvars-1,true,vm,0,k-1,0,nvars-1,false,1.0,ar.m_c,0,nvars-1,0,nvars-1,0.0,work);
            //--- exit the function
            return;
           }
        }
      //--- change value
      info=-255;
      //--- exit the function
      return;
     }
//--- change values
   for(i=0;i<=nvars-1;i++)
     {
      //--- check
      if(sv[i]>epstol*CMath::m_machineepsilon*sv[0])
         svi[i]=1/sv[i];
      else
         svi[i]=0;
     }
//--- change values
   for(i=0;i<=nvars-1;i++)
      t[i]=0;
//--- change values
   for(i=0;i<=npoints-1;i++)
     {
      r=b[i];
      for(i_=0;i_<=nvars-1;i_++)
         t[i_]=t[i_]+r*u[i][i_];
     }
   for(i=0;i<=nvars-1;i++)
      lm.m_w[offs+i]=0;
//--- calculation
   for(i=0;i<=nvars-1;i++)
     {
      r=t[i]*svi[i];
      i1_=-offs;
      for(i_=offs;i_<=offs+nvars-1;i_++)
         lm.m_w[i_]=lm.m_w[i_]+r*vt[i][i_+i1_];
     }
//--- calculation
   for(j=0;j<=nvars-1;j++)
     {
      r=svi[j];
      for(i_=0;i_<=nvars-1;i_++)
         vm[i_].Set(j,r*vt[j][i_]);
     }
//--- calculation
   for(i=0;i<=nvars-1;i++)
     {
      for(j=i;j<=nvars-1;j++)
        {
         r=0.0;
         for(i_=0;i_<=nvars-1;i_++)
            r+=vm[i][i_]*vm[j][i_];
         ar.m_c[i].Set(j,r);
         ar.m_c[j].Set(i,r);
        }
     }
//--- Leave-1-out cross-validation error.
//--- NOTATIONS:
//--- A            design matrix
//--- A*x = b      original linear least squares task
//--- U*S*V'       SVD of A
//--- ai           i-th row of the A
//--- bi           i-th element of the b
//--- xf           solution of the original LLS task
//--- Cross-validation error of i-th element from a sample is
//--- calculated using following formula:
//---     ERRi = ai*xf - (ai*xf-bi*(ui*ui'))/(1-ui*ui')     (1)
//--- This formula can be derived from normal equations of the
//--- original task
//---     (A'*A)x = A'*b                                    (2)
//--- by applying modification (zeroing out i-th row of A) to (2):
//---     (A-ai)'*(A-ai) = (A-ai)'*b
//--- and using Sherman-Morrison formula for updating matrix inverse
//--- NOTE 1: b is not zeroed out since it is much simpler and
//--- does not influence final result.
//--- NOTE 2: some design matrices A have such ui that 1-ui*ui'=0.
//--- Formula (1) can't be applied for such cases and they are skipped
//--- from CV calculation (which distorts resulting CV estimate).
//--- But from the properties of U we can conclude that there can
//--- be no more than NVars such vectors. Usually
//--- NVars << NPoints, so in a normal case it only slightly
//--- influences result.
   ncv=0;
   na=0;
   nacv=0;
   ar.m_rmserror=0;
   ar.m_avgerror=0;
   ar.m_avgrelerror=0;
   ar.m_cvrmserror=0;
   ar.m_cvavgerror=0;
   ar.m_cvavgrelerror=0;
   ar.m_ncvdefects=0;
//--- allocation
   ArrayResizeAL(ar.m_cvdefects,nvars);
   for(i=0;i<=npoints-1;i++)
     {
      //--- Error on a training set
      i1_=offs;
      r=0.0;
      for(i_=0;i_<=nvars-1;i_++)
         r+=xy[i][i_]*lm.m_w[i_+i1_];
      //--- change values
      ar.m_rmserror=ar.m_rmserror+CMath::Sqr(r-xy[i][nvars]);
      ar.m_avgerror=ar.m_avgerror+MathAbs(r-xy[i][nvars]);
      //--- check
      if(xy[i][nvars]!=0.0)
        {
         ar.m_avgrelerror=ar.m_avgrelerror+MathAbs((r-xy[i][nvars])/xy[i][nvars]);
         na=na+1;
        }
      //--- Error using fast leave-one-out cross-validation
      p=0.0;
      for(i_=0;i_<=nvars-1;i_++)
         p+=u[i][i_]*u[i][i_];
      //--- check
      if(p>1-epstol*CMath::m_machineepsilon)
        {
         ar.m_cvdefects[ar.m_ncvdefects]=i;
         ar.m_ncvdefects=ar.m_ncvdefects+1;
         continue;
        }
      //--- change values
      r=s[i]*(r/s[i]-b[i]*p)/(1-p);
      ar.m_cvrmserror=ar.m_cvrmserror+CMath::Sqr(r-xy[i][nvars]);
      ar.m_cvavgerror=ar.m_cvavgerror+MathAbs(r-xy[i][nvars]);
      //--- check
      if(xy[i][nvars]!=0.0)
        {
         ar.m_cvavgrelerror=ar.m_cvavgrelerror+MathAbs((r-xy[i][nvars])/xy[i][nvars]);
         nacv=nacv+1;
        }
      ncv=ncv+1;
     }
//--- check
   if(ncv==0)
     {
      //--- Something strange: ALL ui are degenerate.
      //--- Unexpected...
      info=-255;
      //--- exit the function
      return;
     }
//--- change values
   ar.m_rmserror=MathSqrt(ar.m_rmserror/npoints);
   ar.m_avgerror=ar.m_avgerror/npoints;
//--- check
   if(na!=0)
      ar.m_avgrelerror=ar.m_avgrelerror/na;
   ar.m_cvrmserror=MathSqrt(ar.m_cvrmserror/ncv);
   ar.m_cvavgerror=ar.m_cvavgerror/ncv;
//--- check
   if(nacv!=0)
      ar.m_cvavgrelerror=ar.m_cvavgrelerror/nacv;
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CMLPBase                                     |
//+------------------------------------------------------------------+
class CMultilayerPerceptron
  {
public:
   //--- variables
   int               m_hlnetworktype;
   int               m_hlnormtype;
   //--- arrays
   int               m_hllayersizes[];
   int               m_hlconnections[];
   int               m_hlneurons[];
   int               m_structinfo[];
   double            m_weights[];
   double            m_columnmeans[];
   double            m_columnsigmas[];
   double            m_neurons[];
   double            m_dfdnet[];
   double            m_derror[];
   double            m_x[];
   double            m_y[];
   double            m_nwbuf[];
   int               m_integerbuf[];
   //--- matrix
   CMatrixDouble     m_chunks;
   //--- constructor, destructor
                     CMultilayerPerceptron(void);
                    ~CMultilayerPerceptron(void);
   //--- copy
   void              Copy(CMultilayerPerceptron &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMultilayerPerceptron::CMultilayerPerceptron(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMultilayerPerceptron::~CMultilayerPerceptron(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMultilayerPerceptron::Copy(CMultilayerPerceptron &obj)
  {
//--- copy variables
   m_hlnetworktype=obj.m_hlnetworktype;
   m_hlnormtype=obj.m_hlnormtype;
//--- copy arrays
   ArrayCopy(m_hllayersizes,obj.m_hllayersizes);
   ArrayCopy(m_hlconnections,obj.m_hlconnections);
   ArrayCopy(m_hlneurons,obj.m_hlneurons);
   ArrayCopy(m_structinfo,obj.m_structinfo);
   ArrayCopy(m_weights,obj.m_weights);
   ArrayCopy(m_columnmeans,obj.m_columnmeans);
   ArrayCopy(m_columnsigmas,obj.m_columnsigmas);
   ArrayCopy(m_neurons,obj.m_neurons);
   ArrayCopy(m_dfdnet,obj.m_dfdnet);
   ArrayCopy(m_derror,obj.m_derror);
   ArrayCopy(m_x,obj.m_x);
   ArrayCopy(m_y,obj.m_y);
   ArrayCopy(m_nwbuf,obj.m_nwbuf);
   ArrayCopy(m_integerbuf,obj.m_integerbuf);
//--- copy matrix
   m_chunks=obj.m_chunks;
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CMultilayerPerceptron            |
//+------------------------------------------------------------------+
class CMultilayerPerceptronShell
  {
private:
   CMultilayerPerceptron m_innerobj;
public:
   //--- constructors, destructor
                     CMultilayerPerceptronShell(void);
                     CMultilayerPerceptronShell(CMultilayerPerceptron &obj);
                    ~CMultilayerPerceptronShell(void);
   //--- method
   CMultilayerPerceptron *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMultilayerPerceptronShell::CMultilayerPerceptronShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMultilayerPerceptronShell::CMultilayerPerceptronShell(CMultilayerPerceptron &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMultilayerPerceptronShell::~CMultilayerPerceptronShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMultilayerPerceptron *CMultilayerPerceptronShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Multilayer perceptron class                                      |
//+------------------------------------------------------------------+
class CMLPBase
  {
private:
   //--- private methods
   static void       AddInputLayer(const int ncount,int &lsizes[],int &ltypes[],int &lconnfirst[],int &lconnlast[],int &lastproc);
   static void       AddBiasedSummatorLayer(const int ncount,int &lsizes[],int &ltypes[],int &lconnfirst[],int &lconnlast[],int &lastproc);
   static void       AddActivationLayer(const int functype,int &lsizes[],int &ltypes[],int &lconnfirst[],int &lconnlast[],int &lastproc);
   static void       AddZeroLayer(int &lsizes[],int &ltypes[],int &lconnfirst[],int &lconnlast[],int &lastproc);
   static void       HLAddInputLayer(CMultilayerPerceptron &network,int &connidx,int &neuroidx,int &structinfoidx,int nin);
   static void       HLAddOutputLayer(CMultilayerPerceptron &network,int &connidx,int &neuroidx,int &structinfoidx,int &weightsidx,const int k,const int nprev,const int nout,const bool iscls,const bool islinearout);
   static void       HLAddHiddenLayer(CMultilayerPerceptron &network,int &connidx,int &neuroidx,int &structinfoidx,int &weightsidx,const int k,const int nprev,const int ncur);
   static void       FillHighLevelInformation(CMultilayerPerceptron &network,const int nin,const int nhid1,const int nhid2,const int nout,const bool iscls,const bool islinearout);
   static void       MLPCreate(const int nin,const int nout,int &lsizes[],int &ltypes[],int &lconnfirst[],int &lconnlast[],const int layerscount,const bool isclsnet,CMultilayerPerceptron &network);
   static void       MLPHessianBatchInternal(CMultilayerPerceptron &network,CMatrixDouble &xy,const int ssize,const bool naturalerr,double &e,double &grad[],CMatrixDouble &h);
   static void       MLPInternalCalculateGradient(CMultilayerPerceptron &network,double &neurons[],double &weights[],double &derror[],double &grad[],const bool naturalerrorfunc);
   static void       MLPChunkedGradient(CMultilayerPerceptron &network,CMatrixDouble &xy,const int cstart,const int csize,double &e,double &grad[],const bool naturalerrorfunc);
   static double     SafeCrossEntropy(const double t,const double z);
public:
   //--- variables
   static const int  m_mlpvnum;
   static const int  m_mlpfirstversion;
   static const int  m_nfieldwidth;
   static const int  m_hlconm_nfieldwidth;
   static const int  m_hlm_nfieldwidth;
   static const int  m_chunksize;
   //--- constructor, destructor
                     CMLPBase(void);
                    ~CMLPBase(void);
   //--- public methods
   static void       MLPCreate0(const int nin,const int nout,CMultilayerPerceptron &network);
   static void       MLPCreate1(const int nin,const int nhid,const int nout,CMultilayerPerceptron &network);
   static void       MLPCreate2(const int nin,const int nhid1,const int nhid2,const int nout,CMultilayerPerceptron &network);
   static void       MLPCreateB0(const int nin,const int nout,const double b,double d,CMultilayerPerceptron &network);
   static void       MLPCreateB1(const int nin,const int nhid,const int nout,const double b,double d,CMultilayerPerceptron &network);
   static void       MLPCreateB2(const int nin,const int nhid1,const int nhid2,const int nout,const double b,double d,CMultilayerPerceptron &network);
   static void       MLPCreateR0(const int nin,const int nout,const double a,const double b,CMultilayerPerceptron &network);
   static void       MLPCreateR1(const int nin,const int nhid,const int nout,const double a,const double b,CMultilayerPerceptron &network);
   static void       MLPCreateR2(const int nin,const int nhid1,const int nhid2,const int nout,const double a,const double b,CMultilayerPerceptron &network);
   static void       MLPCreateC0(const int nin,const int nout,CMultilayerPerceptron &network);
   static void       MLPCreateC1(const int nin,const int nhid,const int nout,CMultilayerPerceptron &network);
   static void       MLPCreateC2(const int nin,const int nhid1,const int nhid2,const int nout,CMultilayerPerceptron &network);
   static void       MLPCopy(CMultilayerPerceptron &network1,CMultilayerPerceptron &network2);
   static void       MLPSerializeOld(CMultilayerPerceptron &network,double &ra[],int &rlen);
   static void       MLPUnserializeOld(double &ra[],CMultilayerPerceptron &network);
   static void       MLPRandomize(CMultilayerPerceptron &network);
   static void       MLPRandomizeFull(CMultilayerPerceptron &network);
   static void       MLPInitPreprocessor(CMultilayerPerceptron &network,CMatrixDouble &xy,const int ssize);
   static void       MLPProperties(CMultilayerPerceptron &network,int &nin,int &nout,int &wcount);
   static bool       MLPIsSoftMax(CMultilayerPerceptron &network);
   static int        MLPGetLayersCount(CMultilayerPerceptron &network);
   static int        MLPGetLayerSize(CMultilayerPerceptron &network,const int k);
   static void       MLPGetInputScaling(CMultilayerPerceptron &network,const int i,double &mean,double &sigma);
   static void       MLPGetOutputScaling(CMultilayerPerceptron &network,const int i,double &mean,double &sigma);
   static void       MLPGetNeuronInfo(CMultilayerPerceptron &network,const int k,const int i,int &fkind,double &threshold);
   static double     MLPGetWeight(CMultilayerPerceptron &network,const int k0,const int i0,const int k1,const int i1);
   static void       MLPSetInputScaling(CMultilayerPerceptron &network,const int i,const double mean,double sigma);
   static void       MLPSetOutputScaling(CMultilayerPerceptron &network,const int i,const double mean,double sigma);
   static void       MLPSetNeuronInfo(CMultilayerPerceptron &network,const int k,const int i,const int fkind,const double threshold);
   static void       MLPSetWeight(CMultilayerPerceptron &network,const int k0,const int i0,const int k1,const int i1,const double w);
   static void       MLPActivationFunction(double net,const int k,double &f,double &df,double &d2f);
   static void       MLPProcess(CMultilayerPerceptron &network,double &x[],double &y[]);
   static void       MLPProcessI(CMultilayerPerceptron &network,double &x[],double &y[]);
   static double     MLPError(CMultilayerPerceptron &network,CMatrixDouble &xy,const int ssize);
   static double     MLPErrorN(CMultilayerPerceptron &network,CMatrixDouble &xy,const int ssize);
   static int        MLPClsError(CMultilayerPerceptron &network,CMatrixDouble &xy,const int ssize);
   static double     MLPRelClsError(CMultilayerPerceptron &network,CMatrixDouble &xy,const int npoints);
   static double     MLPAvgCE(CMultilayerPerceptron &network,CMatrixDouble &xy,const int npoints);
   static double     MLPRMSError(CMultilayerPerceptron &network,CMatrixDouble &xy,const int npoints);
   static double     MLPAvgError(CMultilayerPerceptron &network,CMatrixDouble &xy,const int npoints);
   static double     MLPAvgRelError(CMultilayerPerceptron &network,CMatrixDouble &xy,const int npoints);
   static void       MLPGrad(CMultilayerPerceptron &network,double &x[],double &desiredy[],double &e,double &grad[]);
   static void       MLPGradN(CMultilayerPerceptron &network,double &x[],double &desiredy[],double &e,double &grad[]);
   static void       MLPGradBatch(CMultilayerPerceptron &network,CMatrixDouble &xy,const int ssize,double &e,double &grad[]);
   static void       MLPGradNBatch(CMultilayerPerceptron &network,CMatrixDouble &xy,const int ssize,double &e,double &grad[]);
   static void       MLPHessianNBatch(CMultilayerPerceptron &network,CMatrixDouble &xy,const int ssize,double &e,double &grad[],CMatrixDouble &h);
   static void       MLPHessianBatch(CMultilayerPerceptron &network,CMatrixDouble &xy,const int ssize,double &e,double &grad[],CMatrixDouble &h);
   static void       MLPInternalProcessVector(int &structinfo[],double &weights[],double &columnmeans[],double &columnsigmas[],double &neurons[],double &dfdnet[],double &x[],double &y[]);
   static void       MLPAlloc(CSerializer &s,CMultilayerPerceptron &network);
   static void       MLPSerialize(CSerializer &s,CMultilayerPerceptron &network);
   static void       MLPUnserialize(CSerializer &s,CMultilayerPerceptron &network);
  };
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const int CMLPBase::m_mlpvnum=7;
const int CMLPBase::m_mlpfirstversion=0;
const int CMLPBase::m_nfieldwidth=4;
const int CMLPBase::m_hlconm_nfieldwidth=5;
const int CMLPBase::m_hlm_nfieldwidth=4;
const int CMLPBase::m_chunksize=32;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMLPBase::CMLPBase(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMLPBase::~CMLPBase(void)
  {

  }
//+------------------------------------------------------------------+
//| Creates  neural  network  with  NIn  inputs,  NOut outputs,      |
//| without hidden layers, with linear output layer. Network weights | 
//| are filled with small random values.                             |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCreate0(const int nin,const int nout,
                                 CMultilayerPerceptron &network)
  {
//--- create variables
   int layerscount=0;
   int lastproc=0;
//--- creating arrays
   int lsizes[];
   int ltypes[];
   int lconnfirst[];
   int lconnlast[];
//--- initialization
   layerscount=4;
//--- Allocate arrays
   ArrayResizeAL(lsizes,layerscount);
   ArrayResizeAL(ltypes,layerscount);
   ArrayResizeAL(lconnfirst,layerscount);
   ArrayResizeAL(lconnlast,layerscount);
//--- Layers
   AddInputLayer(nin,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nout,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(-5,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- Create
   MLPCreate(nin,nout,lsizes,ltypes,lconnfirst,lconnlast,layerscount,false,network);
//--- function call
   FillHighLevelInformation(network,nin,0,0,nout,false,true);
  }
//+------------------------------------------------------------------+
//| Same as MLPCreate0, but with one hidden layer (NHid neurons) with|
//| non-linear activation function. Output layer is linear.          |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCreate1(const int nin,const int nhid,const int nout,
                                 CMultilayerPerceptron &network)
  {
//--- create variables
   int layerscount=0;
   int lastproc=0;
//--- creating arrays
   int lsizes[];
   int ltypes[];
   int lconnfirst[];
   int lconnlast[];
//--- create variables
   layerscount=7;
//--- Allocate arrays
   ArrayResizeAL(lsizes,layerscount);
   ArrayResizeAL(ltypes,layerscount);
   ArrayResizeAL(lconnfirst,layerscount);
   ArrayResizeAL(lconnlast,layerscount);
//--- Layers
   AddInputLayer(nin,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nhid,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nout,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(-5,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- Create
   MLPCreate(nin,nout,lsizes,ltypes,lconnfirst,lconnlast,layerscount,false,network);
//--- function call
   FillHighLevelInformation(network,nin,nhid,0,nout,false,true);
  }
//+------------------------------------------------------------------+
//| Same as MLPCreate0,but with two hidden layers (NHid1 and NHid2   |
//| neurons) with non-linear activation function. Output layer is    |
//| linear.                                                          |
//|  $ALL                                                            |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCreate2(const int nin,const int nhid1,const int nhid2,
                                 const int nout,CMultilayerPerceptron &network)
  {
//--- create variables
   int layerscount=0;
   int lastproc=0;
//--- creating arrays
   int lsizes[];
   int ltypes[];
   int lconnfirst[];
   int lconnlast[];
//--- initialization
   layerscount=10;
//--- Allocate arrays
   ArrayResizeAL(lsizes,layerscount);
   ArrayResizeAL(ltypes,layerscount);
   ArrayResizeAL(lconnfirst,layerscount);
   ArrayResizeAL(lconnlast,layerscount);
//--- Layers
   AddInputLayer(nin,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nhid1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nhid2,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nout,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(-5,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- Create
   MLPCreate(nin,nout,lsizes,ltypes,lconnfirst,lconnlast,layerscount,false,network);
//--- function call
   FillHighLevelInformation(network,nin,nhid1,nhid2,nout,false,true);
  }
//+------------------------------------------------------------------+
//| Creates neural network with NIn inputs, NOut outputs, without    |
//| hidden layers with non-linear output layer. Network weights are  | 
//| filled with small random values.                                 |
//| Activation function of the output layer takes values:            |
//|     (B, +INF), if D>=0                                           |
//| or                                                               |
//|     (-INF, B), if D<0.                                           |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCreateB0(const int nin,const int nout,const double b,
                                  double d,CMultilayerPerceptron &network)
  {
//--- create variables
   int layerscount=0;
   int lastproc=0;
   int i=0;
//--- creating arrays
   int lsizes[];
   int ltypes[];
   int lconnfirst[];
   int lconnlast[];
//--- initialization
   layerscount=4;
//--- check
   if(d>=0.0)
      d=1;
   else
      d=-1;
//--- Allocate arrays
   ArrayResizeAL(lsizes,layerscount);
   ArrayResizeAL(ltypes,layerscount);
   ArrayResizeAL(lconnfirst,layerscount);
   ArrayResizeAL(lconnlast,layerscount);
//--- Layers
   AddInputLayer(nin,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nout,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(3,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- Create
   MLPCreate(nin,nout,lsizes,ltypes,lconnfirst,lconnlast,layerscount,false,network);
//--- function call
   FillHighLevelInformation(network,nin,0,0,nout,false,false);
//--- Turn on ouputs shift/scaling.
   for(i=nin;i<=nin+nout-1;i++)
     {
      network.m_columnmeans[i]=b;
      network.m_columnsigmas[i]=d;
     }
  }
//+------------------------------------------------------------------+
//| Same as MLPCreateB0 but with non-linear hidden layer.            |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCreateB1(const int nin,const int nhid,const int nout,
                                  const double b,double d,CMultilayerPerceptron &network)
  {
//--- create variables
   int layerscount=0;
   int lastproc=0;
   int i=0;
//--- creating arrays
   int lsizes[];
   int ltypes[];
   int lconnfirst[];
   int lconnlast[];
   layerscount=7;
//--- check
   if(d>=0.0)
      d=1;
   else
      d=-1;
//--- Allocate arrays
   ArrayResizeAL(lsizes,layerscount);
   ArrayResizeAL(ltypes,layerscount);
   ArrayResizeAL(lconnfirst,layerscount);
   ArrayResizeAL(lconnlast,layerscount);
//--- Layers
   AddInputLayer(nin,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nhid,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nout,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(3,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- Create
   MLPCreate(nin,nout,lsizes,ltypes,lconnfirst,lconnlast,layerscount,false,network);
//--- function call
   FillHighLevelInformation(network,nin,nhid,0,nout,false,false);
//--- Turn on ouputs shift/scaling.
   for(i=nin;i<=nin+nout-1;i++)
     {
      network.m_columnmeans[i]=b;
      network.m_columnsigmas[i]=d;
     }
  }
//+------------------------------------------------------------------+
//| Same as MLPCreateB0 but with two non-linear hidden layers.       |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCreateB2(const int nin,const int nhid1,const int nhid2,
                                  const int nout,const double b,double d,
                                  CMultilayerPerceptron &network)
  {
//--- create variables
   int layerscount=0;
   int lastproc=0;
   int i=0;
//--- creating arrays
   int lsizes[];
   int ltypes[];
   int lconnfirst[];
   int lconnlast[];
//--- initialization
   layerscount=10;
//--- check
   if(d>=0.0)
      d=1;
   else
      d=-1;
//--- Allocate arrays
   ArrayResizeAL(lsizes,layerscount);
   ArrayResizeAL(ltypes,layerscount);
   ArrayResizeAL(lconnfirst,layerscount);
   ArrayResizeAL(lconnlast,layerscount);
//--- Layers
   AddInputLayer(nin,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nhid1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nhid2,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nout,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(3,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- Create
   MLPCreate(nin,nout,lsizes,ltypes,lconnfirst,lconnlast,layerscount,false,network);
//--- function call
   FillHighLevelInformation(network,nin,nhid1,nhid2,nout,false,false);
//--- Turn on ouputs shift/scaling.
   for(i=nin;i<=nin+nout-1;i++)
     {
      network.m_columnmeans[i]=b;
      network.m_columnsigmas[i]=d;
     }
  }
//+------------------------------------------------------------------+
//| Creates  neural  network  with  NIn  inputs,  NOut outputs,      |
//| without hidden layers with non-linear output layer. Network      |
//| weights are filled with small random values. Activation function | 
//| of the output layer takes values [A,B].                          |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCreateR0(const int nin,const int nout,const double a,
                                  const double b,CMultilayerPerceptron &network)
  {
//--- create variables
   int layerscount=0;
   int lastproc=0;
   int i=0;
//--- creating arrays
   int lsizes[];
   int ltypes[];
   int lconnfirst[];
   int lconnlast[];
//--- initialization
   layerscount=1+3;
//--- Allocate arrays
   ArrayResizeAL(lsizes,layerscount);
   ArrayResizeAL(ltypes,layerscount);
   ArrayResizeAL(lconnfirst,layerscount);
   ArrayResizeAL(lconnlast,layerscount);
//--- Layers
   AddInputLayer(nin,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nout,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- Create
   MLPCreate(nin,nout,lsizes,ltypes,lconnfirst,lconnlast,layerscount,false,network);
//--- function call
   FillHighLevelInformation(network,nin,0,0,nout,false,false);
//--- Turn on outputs shift/scaling.
   for(i=nin;i<=nin+nout-1;i++)
     {
      network.m_columnmeans[i]=0.5*(a+b);
      network.m_columnsigmas[i]=0.5*(a-b);
     }
  }
//+------------------------------------------------------------------+
//| Same as MLPCreateR0,but with non-linear hidden layer.            |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCreateR1(const int nin,const int nhid,const int nout,
                                  const double a,const double b,
                                  CMultilayerPerceptron &network)
  {
//--- create variables
   int layerscount=0;
   int lastproc=0;
   int i=0;
//--- creating arrays
   int lsizes[];
   int ltypes[];
   int lconnfirst[];
   int lconnlast[];
//--- initialization
   layerscount=7;
//--- Allocate arrays
   ArrayResizeAL(lsizes,layerscount);
   ArrayResizeAL(ltypes,layerscount);
   ArrayResizeAL(lconnfirst,layerscount);
   ArrayResizeAL(lconnlast,layerscount);
//--- Layers
   AddInputLayer(nin,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nhid,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nout,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- Create
   MLPCreate(nin,nout,lsizes,ltypes,lconnfirst,lconnlast,layerscount,false,network);
//--- function call
   FillHighLevelInformation(network,nin,nhid,0,nout,false,false);
//--- Turn on outputs shift/scaling.
   for(i=nin;i<=nin+nout-1;i++)
     {
      network.m_columnmeans[i]=0.5*(a+b);
      network.m_columnsigmas[i]=0.5*(a-b);
     }
  }
//+------------------------------------------------------------------+
//| Same as MLPCreateR0,but with two non-linear hidden layers.       |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCreateR2(const int nin,const int nhid1,const int nhid2,
                                  const int nout,const double a,const double b,
                                  CMultilayerPerceptron &network)
  {
//--- create variables
   int layerscount=0;
   int lastproc=0;
   int i=0;
//--- creating arrays
   int lsizes[];
   int ltypes[];
   int lconnfirst[];
   int lconnlast[];
//--- initialization
   layerscount=10;
//--- Allocate arrays
   ArrayResizeAL(lsizes,layerscount);
   ArrayResizeAL(ltypes,layerscount);
   ArrayResizeAL(lconnfirst,layerscount);
   ArrayResizeAL(lconnlast,layerscount);
//--- Layers
   AddInputLayer(nin,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nhid1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nhid2,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nout,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- Create
   MLPCreate(nin,nout,lsizes,ltypes,lconnfirst,lconnlast,layerscount,false,network);
//--- function call
   FillHighLevelInformation(network,nin,nhid1,nhid2,nout,false,false);
//--- Turn on outputs shift/scaling.
   for(i=nin;i<=nin+nout-1;i++)
     {
      network.m_columnmeans[i]=0.5*(a+b);
      network.m_columnsigmas[i]=0.5*(a-b);
     }
  }
//+------------------------------------------------------------------+
//| Creates classifier network with NIn inputs and NOut possible     |
//| classes.                                                         |
//| Network contains no hidden layers and linear output layer with   |
//| SOFTMAX-normalization (so outputs sums up to 1.0 and converge to |
//| posterior probabilities).                                        |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCreateC0(const int nin,const int nout,
                                  CMultilayerPerceptron &network)
  {
//--- create variables
   int layerscount=0;
   int lastproc=0;
//--- creating arrays
   int lsizes[];
   int ltypes[];
   int lconnfirst[];
   int lconnlast[];
//--- check
   if(!CAp::Assert(nout>=2,__FUNCTION__+": NOut<2!"))
      return;
//--- initialization
   layerscount=4;
//--- Allocate arrays
   ArrayResizeAL(lsizes,layerscount);
   ArrayResizeAL(ltypes,layerscount);
   ArrayResizeAL(lconnfirst,layerscount);
   ArrayResizeAL(lconnlast,layerscount);
//--- Layers
   AddInputLayer(nin,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nout-1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddZeroLayer(lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- Create
   MLPCreate(nin,nout,lsizes,ltypes,lconnfirst,lconnlast,layerscount,true,network);
//--- function call
   FillHighLevelInformation(network,nin,0,0,nout,true,true);
  }
//+------------------------------------------------------------------+
//| Same as MLPCreateC0,but with one non-linear hidden layer.        |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCreateC1(const int nin,const int nhid,const int nout,
                                  CMultilayerPerceptron &network)
  {
//--- create variables
   int layerscount=0;
   int lastproc=0;
//--- creating arrays
   int lsizes[];
   int ltypes[];
   int lconnfirst[];
   int lconnlast[];
//--- check
   if(!CAp::Assert(nout>=2,__FUNCTION__+": NOut<2!"))
      return;
//--- initialization
   layerscount=7;
//--- Allocate arrays
   ArrayResizeAL(lsizes,layerscount);
   ArrayResizeAL(ltypes,layerscount);
   ArrayResizeAL(lconnfirst,layerscount);
   ArrayResizeAL(lconnlast,layerscount);
//--- Layers
   AddInputLayer(nin,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nhid,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nout-1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddZeroLayer(lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- Create
   MLPCreate(nin,nout,lsizes,ltypes,lconnfirst,lconnlast,layerscount,true,network);
//--- function call
   FillHighLevelInformation(network,nin,nhid,0,nout,true,true);
  }
//+------------------------------------------------------------------+
//| Same as MLPCreateC0, but with two non-linear hidden layers.      |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCreateC2(const int nin,const int nhid1,const int nhid2,
                                  const int nout,CMultilayerPerceptron &network)
  {
//--- create variables
   int layerscount=0;
   int lastproc=0;
//--- creating arrays
   int lsizes[];
   int ltypes[];
   int lconnfirst[];
   int lconnlast[];
//--- check
   if(!CAp::Assert(nout>=2,__FUNCTION__+": NOut<2!"))
      return;
//--- initialization
   layerscount=10;
//--- Allocate arrays
   ArrayResizeAL(lsizes,layerscount);
   ArrayResizeAL(ltypes,layerscount);
   ArrayResizeAL(lconnfirst,layerscount);
   ArrayResizeAL(lconnlast,layerscount);
//--- Layers
   AddInputLayer(nin,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nhid1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nhid2,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddActivationLayer(1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddBiasedSummatorLayer(nout-1,lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- function call
   AddZeroLayer(lsizes,ltypes,lconnfirst,lconnlast,lastproc);
//--- Create
   MLPCreate(nin,nout,lsizes,ltypes,lconnfirst,lconnlast,layerscount,true,network);
//--- function call
   FillHighLevelInformation(network,nin,nhid1,nhid2,nout,true,true);
  }
//+------------------------------------------------------------------+
//| Copying of neural network                                        |
//| INPUT PARAMETERS:                                                |
//|     Network1 -   original                                        |
//| OUTPUT PARAMETERS:                                               |
//|     Network2 -   copy                                            |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCopy(CMultilayerPerceptron &network1,
                              CMultilayerPerceptron &network2)
  {
//--- copy
   network2.m_hlnetworktype=network1.m_hlnetworktype;
   network2.m_hlnormtype=network1.m_hlnormtype;
//--- function calls
   CApServ::CopyIntegerArray(network1.m_hllayersizes,network2.m_hllayersizes);
   CApServ::CopyIntegerArray(network1.m_hlconnections,network2.m_hlconnections);
   CApServ::CopyIntegerArray(network1.m_hlneurons,network2.m_hlneurons);
   CApServ::CopyIntegerArray(network1.m_structinfo,network2.m_structinfo);
   CApServ::CopyRealArray(network1.m_weights,network2.m_weights);
   CApServ::CopyRealArray(network1.m_columnmeans,network2.m_columnmeans);
   CApServ::CopyRealArray(network1.m_columnsigmas,network2.m_columnsigmas);
   CApServ::CopyRealArray(network1.m_neurons,network2.m_neurons);
   CApServ::CopyRealArray(network1.m_dfdnet,network2.m_dfdnet);
   CApServ::CopyRealArray(network1.m_derror,network2.m_derror);
   CApServ::CopyRealArray(network1.m_x,network2.m_x);
   CApServ::CopyRealArray(network1.m_y,network2.m_y);
   CApServ::CopyRealMatrix(network1.m_chunks,network2.m_chunks);
   CApServ::CopyRealArray(network1.m_nwbuf,network2.m_nwbuf);
   CApServ::CopyIntegerArray(network1.m_integerbuf,network2.m_integerbuf);
  }
//+------------------------------------------------------------------+
//| Serialization of MultiLayerPerceptron strucure                   |
//| INPUT PARAMETERS:                                                |
//|     Network -   original                                         |
//| OUTPUT PARAMETERS:                                               |
//|     RA      -   array of real numbers which stores network,      |
//|                 array[0..RLen-1]                                 |
//|     RLen    -   RA lenght                                        |
//+------------------------------------------------------------------+
static void CMLPBase::MLPSerializeOld(CMultilayerPerceptron &network,
                                      double &ra[],int &rlen)
  {
//--- create variables
   int i=0;
   int ssize=0;
   int ntotal=0;
   int nin=0;
   int nout=0;
   int wcount=0;
   int sigmalen=0;
   int offs=0;
   int i_=0;
   int i1_=0;
//--- initialization
   rlen=0;
//--- Unload info
   ssize=network.m_structinfo[0];
   nin=network.m_structinfo[1];
   nout=network.m_structinfo[2];
   ntotal=network.m_structinfo[3];
   wcount=network.m_structinfo[4];
//--- check
   if(MLPIsSoftMax(network))
      sigmalen=nin;
   else
      sigmalen=nin+nout;
//---  RA format:
//---      LEN         DESRC.
//---      1           RLen
//---      1           version (MLPVNum)
//---      1           StructInfo size
//---      SSize       StructInfo
//---      WCount      Weights
//---      SigmaLen    ColumnMeans
//---      SigmaLen    ColumnSigmas
   rlen=3+ssize+wcount+2*sigmalen;
//--- allocation
   ArrayResizeAL(ra,rlen);
//--- change values
   ra[0]=rlen;
   ra[1]=m_mlpvnum;
   ra[2]=ssize;
//--- calculation
   offs=3;
   for(i=0;i<=ssize-1;i++)
      ra[offs+i]=network.m_structinfo[i];
//--- calculation
   offs=offs+ssize;
   i1_=-offs;
   for(i_=offs;i_<=offs+wcount-1;i_++)
      ra[i_]=network.m_weights[i_+i1_];
//--- calculation
   offs=offs+wcount;
   i1_=-offs;
   for(i_=offs;i_<=offs+sigmalen-1;i_++)
      ra[i_]=network.m_columnmeans[i_+i1_];
//--- calculation
   offs=offs+sigmalen;
   i1_=-offs;
   for(i_=offs;i_<=offs+sigmalen-1;i_++)
      ra[i_]=network.m_columnsigmas[i_+i1_];
   offs=offs+sigmalen;
  }
//+------------------------------------------------------------------+
//| Unserialization of MultiLayerPerceptron strucure                 |
//| INPUT PARAMETERS:                                                |
//|     RA      -   real array which stores network                  |
//| OUTPUT PARAMETERS:                                               |
//|     Network -   restored network                                 |
//+------------------------------------------------------------------+
static void CMLPBase::MLPUnserializeOld(double &ra[],CMultilayerPerceptron &network)
  {
//--- create variables
   int i=0;
   int ssize=0;
   int ntotal=0;
   int nin=0;
   int nout=0;
   int wcount=0;
   int sigmalen=0;
   int offs=0;
   int i_=0;
   int i1_=0;
//--- check
   if(!CAp::Assert((int)MathRound(ra[1])==m_mlpvnum,__FUNCTION__+": incorrect array!"))
      return;
//--- Unload StructInfo from IA
   offs=3;
   ssize=(int)MathRound(ra[2]);
//--- allocation
   ArrayResizeAL(network.m_structinfo,ssize);
   for(i=0;i<=ssize-1;i++)
      network.m_structinfo[i]=(int)MathRound(ra[offs+i]);
   offs=offs+ssize;
//--- Unload info from StructInfo
   ssize=network.m_structinfo[0];
   nin=network.m_structinfo[1];
   nout=network.m_structinfo[2];
   ntotal=network.m_structinfo[3];
   wcount=network.m_structinfo[4];
//--- check
   if(network.m_structinfo[6]==0)
      sigmalen=nin+nout;
   else
      sigmalen=nin;
//--- Allocate space for other fields
   ArrayResizeAL(network.m_weights,wcount);
   ArrayResizeAL(network.m_columnmeans,sigmalen);
   ArrayResizeAL(network.m_columnsigmas,sigmalen);
   ArrayResizeAL(network.m_neurons,ntotal);
   network.m_chunks.Resize(3*ntotal+1,m_chunksize);
   ArrayResizeAL(network.m_nwbuf,MathMax(wcount,2*nout));
   ArrayResizeAL(network.m_dfdnet,ntotal);
   ArrayResizeAL(network.m_x,nin);
   ArrayResizeAL(network.m_y,nout);
   ArrayResizeAL(network.m_derror,ntotal);
//--- Copy parameters from RA
   i1_=offs;
   for(i_=0;i_<=wcount-1;i_++)
      network.m_weights[i_]=ra[i_+i1_];
//--- calculation
   offs=offs+wcount;
   i1_=offs;
   for(i_=0;i_<=sigmalen-1;i_++)
      network.m_columnmeans[i_]=ra[i_+i1_];
//--- calculation
   offs=offs+sigmalen;
   i1_=offs;
   for(i_=0;i_<=sigmalen-1;i_++)
      network.m_columnsigmas[i_]=ra[i_+i1_];
   offs=offs+sigmalen;
  }
//+------------------------------------------------------------------+
//| Randomization of neural network weights                          |
//+------------------------------------------------------------------+
static void CMLPBase::MLPRandomize(CMultilayerPerceptron &network)
  {
//--- create variables
   int i=0;
   int nin=0;
   int nout=0;
   int wcount=0;
//--- function call
   MLPProperties(network,nin,nout,wcount);
//--- change values
   for(i=0;i<=wcount-1;i++)
      network.m_weights[i]=CMath::RandomReal()-0.5;
  }
//+------------------------------------------------------------------+
//| Randomization of neural network weights and standartisator       |
//+------------------------------------------------------------------+
static void CMLPBase::MLPRandomizeFull(CMultilayerPerceptron &network)
  {
//--- create variables
   int i=0;
   int nin=0;
   int nout=0;
   int wcount=0;
   int ntotal=0;
   int istart=0;
   int offs=0;
   int ntype=0;
//--- function call
   MLPProperties(network,nin,nout,wcount);
//--- initialization
   ntotal=network.m_structinfo[3];
   istart=network.m_structinfo[5];
//--- Process network
   for(i=0;i<=wcount-1;i++)
      network.m_weights[i]=CMath::RandomReal()-0.5;
   for(i=0;i<=nin-1;i++)
     {
      network.m_columnmeans[i]=2*CMath::RandomReal()-1;
      network.m_columnsigmas[i]=1.5*CMath::RandomReal()+0.5;
     }
//--- check
   if(!MLPIsSoftMax(network))
     {
      for(i=0;i<=nout-1;i++)
        {
         offs=istart+(ntotal-nout+i)*m_nfieldwidth;
         ntype=network.m_structinfo[offs+0];
         //--- check
         if(ntype==0)
           {
            //--- Shifts are changed only for linear outputs neurons
            network.m_columnmeans[nin+i]=2*CMath::RandomReal()-1;
           }
         //--- check
         if(ntype==0 || ntype==3)
           {
            //--- Scales are changed only for linear or bounded outputs neurons.
            //--- Note that scale randomization preserves sign.
            network.m_columnsigmas[nin+i]=MathSign(network.m_columnsigmas[nin+i])*(1.5*CMath::RandomReal()+0.5);
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine.                                             |
//+------------------------------------------------------------------+
static void CMLPBase::MLPInitPreprocessor(CMultilayerPerceptron &network,
                                          CMatrixDouble &xy,const int ssize)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    jmax=0;
   int    nin=0;
   int    nout=0;
   int    wcount=0;
   int    ntotal=0;
   int    istart=0;
   int    offs=0;
   int    ntype=0;
   double s=0;
//--- creating arrays
   double means[];
   double sigmas[];
//--- function call
   MLPProperties(network,nin,nout,wcount);
//--- initialization
   ntotal=network.m_structinfo[3];
   istart=network.m_structinfo[5];
//--- Means/Sigmas
   if(MLPIsSoftMax(network))
      jmax=nin-1;
   else
      jmax=nin+nout-1;
//--- allocation
   ArrayResizeAL(means,jmax+1);
   ArrayResizeAL(sigmas,jmax+1);
//--- calculation
   for(j=0;j<=jmax;j++)
     {
      //--- means
      means[j]=0;
      for(i=0;i<=ssize-1;i++)
         means[j]=means[j]+xy[i][j];
      means[j]=means[j]/ssize;
      //--- sigmas
      sigmas[j]=0;
      for(i=0;i<=ssize-1;i++)
         sigmas[j]=sigmas[j]+CMath::Sqr(xy[i][j]-means[j]);
      sigmas[j]=MathSqrt(sigmas[j]/ssize);
     }
//--- Inputs
   for(i=0;i<=nin-1;i++)
     {
      network.m_columnmeans[i]=means[i];
      network.m_columnsigmas[i]=sigmas[i];
      //--- check
      if(network.m_columnsigmas[i]==0.0)
         network.m_columnsigmas[i]=1;
     }
//--- Outputs
   if(!MLPIsSoftMax(network))
     {
      for(i=0;i<=nout-1;i++)
        {
         offs=istart+(ntotal-nout+i)*m_nfieldwidth;
         ntype=network.m_structinfo[offs+0];
         //--- Linear outputs
         if(ntype==0)
           {
            network.m_columnmeans[nin+i]=means[nin+i];
            network.m_columnsigmas[nin+i]=sigmas[nin+i];
            //--- check
            if(network.m_columnsigmas[nin+i]==0.0)
               network.m_columnsigmas[nin+i]=1;
           }
         //--- Bounded outputs (half-interval)
         if(ntype==3)
           {
            s=means[nin+i]-network.m_columnmeans[nin+i];
            //--- check
            if(s==0.0)
               s=MathSign(network.m_columnsigmas[nin+i]);
            //--- check
            if(s==0.0)
               s=1.0;
            //--- change value
            network.m_columnsigmas[nin+i]=MathSign(network.m_columnsigmas[nin+i])*MathAbs(s);
            //--- check
            if((double)(network.m_columnsigmas[nin+i])==0.0)
               network.m_columnsigmas[nin+i]=1;
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Returns information about initialized network: number of inputs, | 
//| outputs, weights.                                                |
//+------------------------------------------------------------------+
static void CMLPBase::MLPProperties(CMultilayerPerceptron &network,int &nin,
                                    int &nout,int &wcount)
  {
//--- change values
   nin=network.m_structinfo[1];
   nout=network.m_structinfo[2];
   wcount=network.m_structinfo[4];
  }
//+------------------------------------------------------------------+
//| Tells whether network is SOFTMAX-normalized (i.m_e. classifier)  |
//| or not.                                                          |
//+------------------------------------------------------------------+
static bool CMLPBase::MLPIsSoftMax(CMultilayerPerceptron &network)
  {
//--- check
   if(network.m_structinfo[6]==1)
      return(true);
//--- return result
   return(false);
  }
//+------------------------------------------------------------------+
//| This function returns total number of layers (including input,   |
//| hidden and output layers).                                       |
//+------------------------------------------------------------------+
static int CMLPBase::MLPGetLayersCount(CMultilayerPerceptron &network)
  {
//--- return result
   return(CAp::Len(network.m_hllayersizes));
  }
//+------------------------------------------------------------------+
//| This function returns size of K-th layer.                        |
//| K=0 corresponds to input layer, K=CNT-1 corresponds to output    |
//| layer.                                                           |
//| Size of the output layer is always equal to the number of        |
//| outputs, although when we have softmax-normalized network, last  | 
//| neuron doesn't have any connections - it is just zero.           |
//+------------------------------------------------------------------+
static int CMLPBase::MLPGetLayerSize(CMultilayerPerceptron &network,
                                     const int k)
  {
//--- check
   if(!CAp::Assert(k>=0 && k<CAp::Len(network.m_hllayersizes),__FUNCTION__+": incorrect layer index"))
      return(-1);
//--- return result
   return(network.m_hllayersizes[k]);
  }
//+------------------------------------------------------------------+
//| This function returns offset/scaling coefficients for I-th input | 
//| of the network.                                                  |
//| INPUT PARAMETERS:                                                |
//|     Network     -   network                                      |
//|     I           -   input index                                  |
//| OUTPUT PARAMETERS:                                               |
//|     Mean        -   mean term                                    |
//|     Sigma       -   sigma term,guaranteed to be nonzero.         |
//| I-th input is passed through linear transformation               |
//|     IN[i]=(IN[i]-Mean)/Sigma                                     |
//| before feeding to the network                                    |
//+------------------------------------------------------------------+
static void CMLPBase::MLPGetInputScaling(CMultilayerPerceptron &network,
                                         const int i,double &mean,
                                         double &sigma)
  {
//--- initialization
   mean=0;
   sigma=0;
//--- check
   if(!CAp::Assert(i>=0 && i<network.m_hllayersizes[0],__FUNCTION__+": incorrect (nonexistent) I"))
      return;
//--- change values
   mean=network.m_columnmeans[i];
   sigma=network.m_columnsigmas[i];
//--- check
   if(sigma==0.0)
      sigma=1;
  }
//+------------------------------------------------------------------+
//| This function returns offset/scaling coefficients for I-th output|
//| of the network.                                                  |
//| INPUT PARAMETERS:                                                |
//|     Network     -   network                                      |
//|     I           -   input index                                  |
//| OUTPUT PARAMETERS:                                               |
//|     Mean        -   mean term                                    |
//|     Sigma       -   sigma term, guaranteed to be nonzero.        |
//| I-th output is passed through linear transformation              |
//|     OUT[i] = OUT[i]*Sigma+Mean                                   |
//| before returning it to user. In case we have SOFTMAX-normalized  |
//| network, we return (Mean,Sigma)=(0.0,1.0).                       |
//+------------------------------------------------------------------+
static void CMLPBase::MLPGetOutputScaling(CMultilayerPerceptron &network,
                                          const int i,double &mean,
                                          double &sigma)
  {
//--- initialization
   mean=0;
   sigma=0;
//--- check
   if(!CAp::Assert(i>=0 && i<network.m_hllayersizes[CAp::Len(network.m_hllayersizes)-1],__FUNCTION__+": incorrect (nonexistent) I"))
      return;
//--- check
   if(network.m_structinfo[6]==1)
     {
      //--- change values
      mean=0;
      sigma=1;
     }
   else
     {
      //--- change values
      mean=network.m_columnmeans[network.m_hllayersizes[0]+i];
      sigma=network.m_columnsigmas[network.m_hllayersizes[0]+i];
     }
  }
//+------------------------------------------------------------------+
//| This function returns information about Ith neuron of Kth layer  |
//| INPUT PARAMETERS:                                                |
//|     Network     -   network                                      |
//|     K           -   layer index                                  |
//|     I           -   neuron index (within layer)                  |
//| OUTPUT PARAMETERS:                                               |
//|     FKind       -   activation function type (used by            |
//|                     MLPActivationFunction()) this value is zero  |
//|                     for input or linear neurons                  |
//|     Threshold   -   also called offset, bias                     |
//|                     zero for input neurons                       |
//| NOTE: this function throws exception if layer or neuron with     |
//| given index do not exists.                                       |
//+------------------------------------------------------------------+
static void CMLPBase::MLPGetNeuronInfo(CMultilayerPerceptron &network,
                                       const int k,const int i,int &fkind,
                                       double &threshold)
  {
//--- create variables
   int ncnt=0;
   int istart=0;
   int highlevelidx=0;
   int activationoffset=0;
//--- initialization
   fkind=0;
   threshold=0;
   ncnt=CAp::Len(network.m_hlneurons)/m_hlm_nfieldwidth;
   istart=network.m_structinfo[5];
//--- search
   network.m_integerbuf[0]=k;
   network.m_integerbuf[1]=i;
//--- function call
   highlevelidx=CApServ::RecSearch(network.m_hlneurons,m_hlm_nfieldwidth,2,0,ncnt,network.m_integerbuf);
//--- check
   if(!CAp::Assert(highlevelidx>=0,__FUNCTION__+": incorrect (nonexistent) layer or neuron index"))
      return;
//--- 1. find offset of the activation function record in the
   if(network.m_hlneurons[highlevelidx*m_hlm_nfieldwidth+2]>=0)
     {
      activationoffset=istart+network.m_hlneurons[highlevelidx*m_hlm_nfieldwidth+2]*m_nfieldwidth;
      fkind=network.m_structinfo[activationoffset+0];
     }
   else
      fkind=0;
//--- check
   if(network.m_hlneurons[highlevelidx*m_hlm_nfieldwidth+3]>=0)
      threshold=network.m_weights[network.m_hlneurons[highlevelidx*m_hlm_nfieldwidth+3]];
   else
      threshold=0;
  }
//+------------------------------------------------------------------+
//| This function returns information about connection from I0-th    |
//| neuron of K0-th layer to I1-th neuron of K1-th layer.            |
//| INPUT PARAMETERS:                                                |
//|     Network     -   network                                      |
//|     K0          -   layer index                                  |
//|     I0          -   neuron index (within layer)                  |
//|     K1          -   layer index                                  |
//|     I1          -   neuron index (within layer)                  |
//| RESULT:                                                          |
//|     connection weight (zero for non-existent connections)        |
//| This function:                                                   |
//| 1. throws exception if layer or neuron with given index do not   |
//|    exists.                                                       |
//| 2. returns zero if neurons exist, but there is no connection     |
//|    between them                                                  |
//+------------------------------------------------------------------+
static double CMLPBase::MLPGetWeight(CMultilayerPerceptron &network,
                                     const int k0,const int i0,
                                     const int k1,const int i1)
  {
//--- create variables
   double result=0;
   int    ccnt=0;
   int    highlevelidx=0;
//--- initialization
   ccnt=CAp::Len(network.m_hlconnections)/m_hlconm_nfieldwidth;
//--- check params
   if(!CAp::Assert(k0>=0 && k0<CAp::Len(network.m_hllayersizes),__FUNCTION__+": incorrect (nonexistent) K0"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(i0>=0 && i0<network.m_hllayersizes[k0],__FUNCTION__+": incorrect (nonexistent) I0"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(k1>=0 && k1<CAp::Len(network.m_hllayersizes),__FUNCTION__+": incorrect (nonexistent) K1"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(i1>=0 && i1<network.m_hllayersizes[k1],__FUNCTION__+": incorrect (nonexistent) I1"))
      return(EMPTY_VALUE);
//--- search
   network.m_integerbuf[0]=k0;
   network.m_integerbuf[1]=i0;
   network.m_integerbuf[2]=k1;
   network.m_integerbuf[3]=i1;
//--- function call
   highlevelidx=CApServ::RecSearch(network.m_hlconnections,m_hlconm_nfieldwidth,4,0,ccnt,network.m_integerbuf);
//--- check
   if(highlevelidx>=0)
      result=network.m_weights[network.m_hlconnections[highlevelidx*m_hlconm_nfieldwidth+4]];
   else
      result=0;
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| This function sets offset/scaling coefficients for I-th input of | 
//| the network.                                                     |
//| INPUT PARAMETERS:                                                |
//|     Network     -   network                                      |
//|     I           -   input index                                  |
//|     Mean        -   mean term                                    |
//|     Sigma       -   sigma term (if zero,will be replaced by 1.0) |
//| NTE: I-th input is passed through linear transformation          |
//|     IN[i]=(IN[i]-Mean)/Sigma                                     |
//| before feeding to the network. This function sets Mean and Sigma.|
//+------------------------------------------------------------------+
static void CMLPBase::MLPSetInputScaling(CMultilayerPerceptron &network,
                                         const int i,const double mean,
                                         double sigma)
  {
//--- check
   if(!CAp::Assert(i>=0 && i<network.m_hllayersizes[0],__FUNCTION__+": incorrect (nonexistent) I"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(mean),__FUNCTION__+": infinite or NAN Mean"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(sigma),__FUNCTION__+": infinite or NAN Sigma"))
      return;
//--- check
   if(sigma==0.0)
      sigma=1;
//--- change values
   network.m_columnmeans[i]=mean;
   network.m_columnsigmas[i]=sigma;
  }
//+------------------------------------------------------------------+
//| This function sets offset/scaling coefficients for I-th output of| 
//| the network.                                                     |
//| INPUT PARAMETERS:                                                |
//|     Network     -   network                                      |
//|     I           -   input index                                  |
//|     Mean        -   mean term                                    |
//|     Sigma       -   sigma term (if zero, will be replaced by 1.0)|
//| OUTPUT PARAMETERS:                                               |
//| NOTE: I-th output is passed through linear transformation        |
//|     OUT[i] = OUT[i]*Sigma+Mean                                   |
//| before returning it to user. This function sets Sigma/Mean. In   | 
//| case we have SOFTMAX-normalized network, you can not set (Sigma, |
//| Mean) to anything other than(0.0,1.0) - this function will throw |
//| exception.                                                       |
//+------------------------------------------------------------------+
static void CMLPBase::MLPSetOutputScaling(CMultilayerPerceptron &network,
                                          const int i,const double mean,
                                          double sigma)
  {
//--- check
   if(!CAp::Assert(i>=0 && i<network.m_hllayersizes[CAp::Len(network.m_hllayersizes)-1],__FUNCTION__+": incorrect (nonexistent) I"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(mean),__FUNCTION__+": infinite or NAN Mean"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(sigma),__FUNCTION__+": infinite or NAN Sigma"))
      return;
//--- check
   if(network.m_structinfo[6]==1)
     {
      //--- check
      if(!CAp::Assert(mean==0.0,__FUNCTION__+": you can not set non-zero Mean term for classifier network"))
         return;
      //--- check
      if(!CAp::Assert(sigma==1.0,__FUNCTION__+": you can not set non-unit Sigma term for classifier network"))
         return;
     }
   else
     {
      //--- check
      if(sigma==0.0)
         sigma=1;
      //--- change values
      network.m_columnmeans[network.m_hllayersizes[0]+i]=mean;
      network.m_columnsigmas[network.m_hllayersizes[0]+i]=sigma;
     }
  }
//+------------------------------------------------------------------+
//| This function modifies information about Ith neuron of Kth layer |
//| INPUT PARAMETERS:                                                |
//|     Network     -   network                                      |
//|     K           -   layer index                                  |
//|     I           -   neuron index (within layer)                  |
//|     FKind       -   activation function type (used by            |
//|                     MLPActivationFunction()) this value must be  |
//|                     zero for input neurons (you can not set      |
//|                     activation function for input neurons)       |
//|     Threshold   -   also called offset, bias                     |
//|                     this value must be zero for input neurons    |
//|                     (you can not set threshold for input neurons)|
//| NOTES:                                                           |
//| 1. this function throws exception if layer or neuron with given  |
//|    index do not exists.                                          |
//| 2. this function also throws exception when you try to set       |
//|    non-linear activation function for input neurons (any kind    |
//|    of network) or for output neurons of classifier network.      |
//| 3. this function throws exception when you try to set non-zero   |
//|    threshold for input neurons (any kind of network).            |
//+------------------------------------------------------------------+
static void CMLPBase::MLPSetNeuronInfo(CMultilayerPerceptron &network,
                                       const int k,const int i,
                                       const int fkind,const double threshold)
  {
//--- create variables
   int ncnt=0;
   int istart=0;
   int highlevelidx=0;
   int activationoffset=0;
//--- check
   if(!CAp::Assert(CMath::IsFinite(threshold),__FUNCTION__+": infinite or NAN Threshold"))
      return;
//--- convenience vars
   ncnt=CAp::Len(network.m_hlneurons)/m_hlm_nfieldwidth;
   istart=network.m_structinfo[5];
//--- search
   network.m_integerbuf[0]=k;
   network.m_integerbuf[1]=i;
//--- function call
   highlevelidx=CApServ::RecSearch(network.m_hlneurons,m_hlm_nfieldwidth,2,0,ncnt,network.m_integerbuf);
//--- check
   if(!CAp::Assert(highlevelidx>=0,__FUNCTION__+": incorrect (nonexistent) layer or neuron index"))
      return;
//--- activation function
   if(network.m_hlneurons[highlevelidx*m_hlm_nfieldwidth+2]>=0)
     {
      activationoffset=istart+network.m_hlneurons[highlevelidx*m_hlm_nfieldwidth+2]*m_nfieldwidth;
      network.m_structinfo[activationoffset+0]=fkind;
     }
   else
     {
      //--- check
      if(!CAp::Assert(fkind==0,__FUNCTION__+": you try to set activation function for neuron which can not have one"))
         return;
     }
//--- Threshold
   if(network.m_hlneurons[highlevelidx*m_hlm_nfieldwidth+3]>=0)
      network.m_weights[network.m_hlneurons[highlevelidx*m_hlm_nfieldwidth+3]]=threshold;
   else
     {
      //--- check
      if(!CAp::Assert(threshold==0.0,__FUNCTION__+": you try to set non-zero threshold for neuron which can not have one"))
         return;
     }
  }
//+------------------------------------------------------------------+
//| This function modifies information about connection from I0-th   |
//| neuron of K0-th layer to I1-th neuron of K1-th layer.            |
//| INPUT PARAMETERS:                                                |
//|     Network     -   network                                      |
//|     K0          -   layer index                                  |
//|     I0          -   neuron index (within layer)                  |
//|     K1          -   layer index                                  |
//|     I1          -   neuron index (within layer)                  |
//|     W           -   connection weight (must be zero for          |
//|                     non-existent connections)                    |
//| This function:                                                   |
//| 1. throws exception if layer or neuron with given index do not   |
//|    exists.                                                       |
//| 2. throws exception if you try to set non-zero weight for        |
//|    non-existent connection                                       |
//+------------------------------------------------------------------+
static void CMLPBase::MLPSetWeight(CMultilayerPerceptron &network,const int k0,
                                   const int i0,const int k1,
                                   const int i1,const double w)
  {
//--- create variables
   int ccnt=0;
   int highlevelidx=0;
//--- initialization
   ccnt=CAp::Len(network.m_hlconnections)/m_hlconm_nfieldwidth;
//--- check params
   if(!CAp::Assert(k0>=0 && k0<CAp::Len(network.m_hllayersizes),__FUNCTION__+": incorrect (nonexistent) K0"))
      return;
//--- check
   if(!CAp::Assert(i0>=0 && i0<network.m_hllayersizes[k0],__FUNCTION__+": incorrect (nonexistent) I0"))
      return;
//--- check
   if(!CAp::Assert(k1>=0 && k1<CAp::Len(network.m_hllayersizes),__FUNCTION__+": incorrect (nonexistent) K1"))
      return;
//--- check
   if(!CAp::Assert(i1>=0 && i1<network.m_hllayersizes[k1],__FUNCTION__+": incorrect (nonexistent) I1"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(w),__FUNCTION__+": infinite or NAN weight"))
      return;
//--- search
   network.m_integerbuf[0]=k0;
   network.m_integerbuf[1]=i0;
   network.m_integerbuf[2]=k1;
   network.m_integerbuf[3]=i1;
//--- function call
   highlevelidx=CApServ::RecSearch(network.m_hlconnections,m_hlconm_nfieldwidth,4,0,ccnt,network.m_integerbuf);
//--- check
   if(highlevelidx>=0)
      network.m_weights[network.m_hlconnections[highlevelidx*m_hlconm_nfieldwidth+4]]=w;
   else
     {
      //--- check
      if(!CAp::Assert(w==0.0,__FUNCTION__+": you try to set non-zero weight for non-existent connection"))
         return;
     }
  }
//+------------------------------------------------------------------+
//| Neural network activation function                               |
//| INPUT PARAMETERS:                                                |
//|     NET         -   neuron input                                 |
//|     K           -   function index (zero for linear function)    |
//| OUTPUT PARAMETERS:                                               |
//|     F           -   function                                     |
//|     DF          -   its derivative                               |
//|     D2F         -   its second derivative                        |
//+------------------------------------------------------------------+
static void CMLPBase::MLPActivationFunction(double net,const int k,double &f,
                                            double &df,double &d2f)
  {
//--- create variables
   double net2=0;
   double arg=0;
   double root=0;
   double r=0;
//--- initialization
   f=0;
   df=0;
   d2f=0;
//--- check
   if(k==0 || k==-5)
     {
      f=net;
      df=1;
      d2f=0;
      //--- exit the function
      return;
     }
//--- check
   if(k==1)
     {
      //--- TanH activation function
      if(MathAbs(net)<100.0)
         f=MathTanh(net);
      else
         f=MathSign(net);
      //--- change values
      df=1-CMath::Sqr(f);
      d2f=-(2*f*df);
      //--- exit the function
      return;
     }
//--- check
   if(k==3)
     {
      //--- EX activation function
      if(net>=0.0)
        {
         //--- change values
         net2=net*net;
         arg=net2+1;
         root=MathSqrt(arg);
         f=net+root;
         r=net/root;
         df=1+r;
         d2f=(root-net*r)/arg;
        }
      else
        {
         //--- change values
         f=MathExp(net);
         df=f;
         d2f=f;
        }
      //--- exit the function
      return;
     }
//--- check
   if(k==2)
     {
      //--- calculation
      f=MathExp(-CMath::Sqr(net));
      df=-(2*net*f);
      d2f=-(2*(f+df*net));
      //--- exit the function
      return;
     }
//--- change values
   f=0;
   df=0;
   d2f=0;
  }
//+------------------------------------------------------------------+
//| Procesing                                                        |
//| INPUT PARAMETERS:                                                |
//|     Network -   neural network                                   |
//|     X       -   input vector,  array[0..NIn-1].                  |
//| OUTPUT PARAMETERS:                                               |
//|     Y       -   result. Regression estimate when solving         |
//|                 regression task, vector of posterior             |
//|                 probabilities for classification task.           |
//| See also MLPProcessI                                             |
//+------------------------------------------------------------------+
static void CMLPBase::MLPProcess(CMultilayerPerceptron &network,double &x[],
                                 double &y[])
  {
//--- check
   if(CAp::Len(y)<network.m_structinfo[2])
      ArrayResizeAL(y,network.m_structinfo[2]);
//--- function call
   MLPInternalProcessVector(network.m_structinfo,network.m_weights,network.m_columnmeans,network.m_columnsigmas,network.m_neurons,network.m_dfdnet,x,y);
  }
//+------------------------------------------------------------------+
//| 'interactive' variant of MLPProcess for languages like Python    | 
//| which support constructs like "Y = MLPProcess(NN,X)" and         |
//| interactive mode of the interpreter                              |
//| This function allocates new array on each call, so it is         |
//| significantly slower than its 'non-interactive' counterpart,     |
//| but it is more convenient when you call it from command line.    |
//+------------------------------------------------------------------+
static void CMLPBase::MLPProcessI(CMultilayerPerceptron &network,double &x[],
                                  double &y[])
  {
//--- function call
   MLPProcess(network,x,y);
  }
//+------------------------------------------------------------------+
//| Error function for neural network,internal subroutine.           |
//+------------------------------------------------------------------+
static double CMLPBase::MLPError(CMultilayerPerceptron &network,
                                 CMatrixDouble &xy,const int ssize)
  {
//--- create variables
   double result=0;
   int    i=0;
   int    k=0;
   int    nin=0;
   int    nout=0;
   int    wcount=0;
   double e=0;
   int    i_=0;
   int    i1_=0;
//--- function call
   MLPProperties(network,nin,nout,wcount);
//--- calculation
   for(i=0;i<=ssize-1;i++)
     {
      for(i_=0;i_<=nin-1;i_++)
         network.m_x[i_]=xy[i][i_];
      //--- function call
      MLPProcess(network,network.m_x,network.m_y);
      //--- check
      if(MLPIsSoftMax(network))
        {
         //--- class labels outputs
         k=(int)MathRound(xy[i][nin]);
         //--- check
         if(k>=0 && k<nout)
            network.m_y[k]=network.m_y[k]-1;
        }
      else
        {
         //--- real outputs
         i1_=nin;
         for(i_=0;i_<=nout-1;i_++)
            network.m_y[i_]=network.m_y[i_]-xy[i][i_+i1_];
        }
      //--- calculation
      e=0.0;
      for(i_=0;i_<=nout-1;i_++)
         e+=network.m_y[i_]*network.m_y[i_];
      result=result+e/2;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Natural error function for neural network,internal subroutine.   |
//+------------------------------------------------------------------+
static double CMLPBase::MLPErrorN(CMultilayerPerceptron &network,
                                  CMatrixDouble &xy,const int ssize)
  {
//--- create variables
   double result=0;
   int    i=0;
   int    k=0;
   int    nin=0;
   int    nout=0;
   int    wcount=0;
   double e=0;
   int    i_=0;
   int    i1_=0;
//--- function call
   MLPProperties(network,nin,nout,wcount);
//--- calculation
   for(i=0;i<=ssize-1;i++)
     {
      //--- Process vector
      for(i_=0;i_<=nin-1;i_++)
         network.m_x[i_]=xy[i][i_];
      //--- function call
      MLPProcess(network,network.m_x,network.m_y);
      //--- Update error function
      if(network.m_structinfo[6]==0)
        {
         //--- Least squares error function
         i1_=nin;
         for(i_=0;i_<=nout-1;i_++)
            network.m_y[i_]=network.m_y[i_]-xy[i][i_+i1_];
         //--- calculation
         e=0.0;
         for(i_=0;i_<=nout-1;i_++)
            e+=network.m_y[i_]*network.m_y[i_];
         result=result+e/2;
        }
      else
        {
         //--- Cross-entropy error function
         k=(int)MathRound(xy[i][nin]);
         //--- check
         if(k>=0 && k<nout)
            result=result+SafeCrossEntropy(1,network.m_y[k]);
        }
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Classification error                                             |
//+------------------------------------------------------------------+
static int CMLPBase::MLPClsError(CMultilayerPerceptron &network,
                                 CMatrixDouble &xy,const int ssize)
  {
//--- create variables
   int result=0;
   int i=0;
   int j=0;
   int nin=0;
   int nout=0;
   int wcount=0;
   int nn=0;
   int ns=0;
   int nmax=0;
   int i_=0;
//--- creating arrays
   double workx[];
   double worky[];
//--- function call
   MLPProperties(network,nin,nout,wcount);
//--- allocation
   ArrayResizeAL(workx,nin);
   ArrayResizeAL(worky,nout);
//--- calculation
   for(i=0;i<=ssize-1;i++)
     {
      //--- Process
      for(i_=0;i_<=nin-1;i_++)
         workx[i_]=xy[i][i_];
      //--- function call
      MLPProcess(network,workx,worky);
      //--- Network version of the answer
      nmax=0;
      for(j=0;j<=nout-1;j++)
        {
         //--- check
         if(worky[j]>worky[nmax])
            nmax=j;
        }
      nn=nmax;
      //--- Right answer
      if(MLPIsSoftMax(network))
         ns=(int)MathRound(xy[i][nin]);
      else
        {
         nmax=0;
         for(j=0;j<=nout-1;j++)
           {
            //--- check
            if(xy[i][nin+j]>xy[i][nin+nmax])
               nmax=j;
           }
         ns=nmax;
        }
      //--- compare
      if(nn!=ns)
         result=result+1;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Relative classification error on the test set                    |
//| INPUT PARAMETERS:                                                |
//|     Network -   network                                          |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     percent of incorrectly classified cases. Works both for      |
//|     classifier networks and general purpose networks used as     |
//|     classifiers.                                                 |
//+------------------------------------------------------------------+
static double CMLPBase::MLPRelClsError(CMultilayerPerceptron &network,
                                       CMatrixDouble &xy,const int npoints)
  {
//--- return result
   return((double)MLPClsError(network,xy,npoints)/(double)npoints);
  }
//+------------------------------------------------------------------+
//| Average cross-entropy (in bits per element) on the test set      |
//| INPUT PARAMETERS:                                                |
//|     Network -   neural network                                   |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     CrossEntropy/(NPoints*LN(2)).                                |
//|     Zero if network solves regression task.                      |
//+------------------------------------------------------------------+
static double CMLPBase::MLPAvgCE(CMultilayerPerceptron &network,CMatrixDouble &xy,
                                 const int npoints)
  {
//--- create variables
   double result=0;
   int    nin=0;
   int    nout=0;
   int    wcount=0;
//--- check
   if(MLPIsSoftMax(network))
     {
      //--- function call
      MLPProperties(network,nin,nout,wcount);
      //--- get result
      result=MLPErrorN(network,xy,npoints)/(npoints*MathLog(2));
     }
   else
      result=0;
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| RMS error on the test set                                        |
//| INPUT PARAMETERS:                                                |
//|     Network -   neural network                                   |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     root mean square error.                                      |
//|     Its meaning for regression task is obvious. As for           |
//|     classification task,RMS error means error when estimating    |
//|     posterior probabilities.                                     |
//+------------------------------------------------------------------+
static double CMLPBase::MLPRMSError(CMultilayerPerceptron &network,
                                    CMatrixDouble &xy,const int npoints)
  {
//--- create variables
   int nin=0;
   int nout=0;
   int wcount=0;
//--- function call
   MLPProperties(network,nin,nout,wcount);
//--- return result
   return(MathSqrt(2*MLPError(network,xy,npoints)/(npoints*nout)));
  }
//+------------------------------------------------------------------+
//| Average error on the test set                                    |
//| INPUT PARAMETERS:                                                |
//|     Network -   neural network                                   |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     Its meaning for regression task is obvious. As for           |
//|     classification task,it means average error when estimating   |
//|     posterior probabilities.                                     |
//+------------------------------------------------------------------+
static double CMLPBase::MLPAvgError(CMultilayerPerceptron &network,
                                    CMatrixDouble &xy,const int npoints)
  {
//--- create variables
   double result=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    nin=0;
   int    nout=0;
   int    wcount=0;
   int    i_=0;
//--- function call
   MLPProperties(network,nin,nout,wcount);
//--- calculation
   for(i=0;i<=npoints-1;i++)
     {
      for(i_=0;i_<=nin-1;i_++)
         network.m_x[i_]=xy[i][i_];
      //--- function call
      MLPProcess(network,network.m_x,network.m_y);
      //--- check
      if(MLPIsSoftMax(network))
        {
         //--- class labels
         k=(int)MathRound(xy[i][nin]);
         for(j=0;j<=nout-1;j++)
           {
            //--- check
            if(j==k)
               result=result+MathAbs(1-network.m_y[j]);
            else
               result=result+MathAbs(network.m_y[j]);
           }
        }
      else
        {
         //--- real outputs
         for(j=0;j<=nout-1;j++)
            result=result+MathAbs(xy[i][nin+j]-network.m_y[j]);
        }
     }
//--- return result
   return(result/(npoints*nout));
  }
//+------------------------------------------------------------------+
//| Average relative error on the test set                           |
//| INPUT PARAMETERS:                                                |
//|     Network -   neural network                                   |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     Its meaning for regression task is obvious. As for           |
//|     classification task, it means average relative error when    |
//|     estimating posterior probability of belonging to the correct |
//|     class.                                                       |
//+------------------------------------------------------------------+
static double CMLPBase::MLPAvgRelError(CMultilayerPerceptron &network,
                                       CMatrixDouble &xy,const int npoints)
  {
//--- create variables
   double result=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    lk=0;
   int    nin=0;
   int    nout=0;
   int    wcount=0;
   int    i_=0;
//--- function call
   MLPProperties(network,nin,nout,wcount);
//--- initialization
   result=0;
   k=0;
//--- calculation
   for(i=0;i<=npoints-1;i++)
     {
      for(i_=0;i_<=nin-1;i_++)
         network.m_x[i_]=xy[i][i_];
      //--- function call
      MLPProcess(network,network.m_x,network.m_y);
      //--- check
      if(MLPIsSoftMax(network))
        {
         //--- class labels
         lk=(int)MathRound(xy[i][nin]);
         for(j=0;j<=nout-1;j++)
           {
            //--- check
            if(j==lk)
              {
               result=result+MathAbs(1-network.m_y[j]);
               k=k+1;
              }
           }
        }
      else
        {
         //--- real outputs
         for(j=0;j<=nout-1;j++)
           {
            //--- check
            if(xy[i][nin+j]!=0.0)
              {
               result=result+MathAbs(xy[i][nin+j]-network.m_y[j])/MathAbs(xy[i][nin+j]);
               k=k+1;
              }
           }
        }
     }
//--- check
   if(k!=0)
      result=result/k;
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Gradient calculation                                             |
//| INPUT PARAMETERS:                                                |
//|     Network -   network initialized with one of the network      |
//|                 creation funcs                                   |
//|     X       -   input vector, length of array must be at least   |
//|                 NIn                                              |
//|     DesiredY-   desired outputs, length of array must be at least|
//|                 NOut                                             |
//|     Grad    -   possibly preallocated array. If size of array is |
//|                 smaller than WCount, it will be reallocated. It  |
//|                 is recommended to reuse previously allocated     |
//|                 array to reduce allocation overhead.             |
//| OUTPUT PARAMETERS:                                               |
//|     E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)   |
//|     Grad    -   gradient of E with respect to weights of network,|
//|                 array[WCount]                                    |
//+------------------------------------------------------------------+
static void CMLPBase::MLPGrad(CMultilayerPerceptron &network,double &x[],
                              double &desiredy[],double &e,double &grad[])
  {
//--- create variables
   int i=0;
   int nout=0;
   int ntotal=0;
//--- Alloc
   if(CAp::Len(grad)<network.m_structinfo[4])
      ArrayResizeAL(grad,network.m_structinfo[4]);
//--- Prepare dError/dOut,internal structures
   MLPProcess(network,x,network.m_y);
//--- initialization
   nout=network.m_structinfo[2];
   ntotal=network.m_structinfo[3];
   e=0;
//--- change values
   for(i=0;i<=ntotal-1;i++)
      network.m_derror[i]=0;
   for(i=0;i<=nout-1;i++)
     {
      network.m_derror[ntotal-nout+i]=network.m_y[i]-desiredy[i];
      e=e+CMath::Sqr(network.m_y[i]-desiredy[i])/2;
     }
//--- gradient
   MLPInternalCalculateGradient(network,network.m_neurons,network.m_weights,network.m_derror,grad,false);
  }
//+------------------------------------------------------------------+
//| Gradient calculation (natural error function is used)            |
//| INPUT PARAMETERS:                                                |
//|     Network -   network initialized with one of the network      |
//|                 creation funcs                                   |
//|     X       -   input vector, length of array must be at least   | 
//|                 NIn                                              |
//|     DesiredY-   desired outputs, length of array must be at least| 
//|                 NOut                                             |
//|     Grad    -   possibly preallocated array. If size of array is |
//|                 smaller than WCount, it will be reallocated. It  |
//|                 is recommended to reuse previously allocated     |
//|                 array to reduce allocation overhead.             |
//| OUTPUT PARAMETERS:                                               |
//|     E       -   error function, sum-of-squares for regression    |
//|                 networks, cross-entropy for classification       |
//|                 networks.                                        |
//|     Grad    -   gradient of E with respect to weights of network,| 
//|                 array[WCount]                                    |
//+------------------------------------------------------------------+
static void CMLPBase::MLPGradN(CMultilayerPerceptron &network,double &x[],
                               double &desiredy[],double &e,double &grad[])
  {
//--- create variables
   double s=0;
   int    i=0;
   int    nout=0;
   int    ntotal=0;
//--- initialization
   e=0;
//--- Alloc
   if(CAp::Len(grad)<network.m_structinfo[4])
      ArrayResizeAL(grad,network.m_structinfo[4]);
//--- Prepare dError/dOut,internal structures
   MLPProcess(network,x,network.m_y);
//--- change values
   nout=network.m_structinfo[2];
   ntotal=network.m_structinfo[3];
   for(i=0;i<=ntotal-1;i++)
      network.m_derror[i]=0;
   e=0;
//--- check
   if(network.m_structinfo[6]==0)
     {
      //--- Regression network,least squares
      for(i=0;i<=nout-1;i++)
        {
         network.m_derror[ntotal-nout+i]=network.m_y[i]-desiredy[i];
         e=e+CMath::Sqr(network.m_y[i]-desiredy[i])/2;
        }
     }
   else
     {
      //--- Classification network,cross-entropy
      s=0;
      for(i=0;i<=nout-1;i++)
         s=s+desiredy[i];
      for(i=0;i<=nout-1;i++)
        {
         network.m_derror[ntotal-nout+i]=s*network.m_y[i]-desiredy[i];
         e=e+SafeCrossEntropy(desiredy[i],network.m_y[i]);
        }
     }
//--- gradient
   MLPInternalCalculateGradient(network,network.m_neurons,network.m_weights,network.m_derror,grad,true);
  }
//+------------------------------------------------------------------+
//| Batch gradient calculation for a set of inputs/outputs           |
//| INPUT PARAMETERS:                                                |
//|     Network -   network initialized with one of the network      |
//|                 creation funcs                                   |
//|     XY      -   set of inputs/outputs; one sample = one row;     |
//|                 first NIn columns contain inputs,                |
//|                 next NOut columns - desired outputs.             |
//|     SSize   -   number of elements in XY                         |
//|     Grad    -   possibly preallocated array. If size of array is | 
//|                 smaller than WCount, it will be reallocated. It  | 
//|                 is recommended to reuse previously allocated     |
//|                 array to reduce allocation overhead.             |
//| OUTPUT PARAMETERS:                                               |
//|     E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)   |
//|     Grad    -   gradient of E with respect to weights of network,| 
//|                 array[WCount]                                    |
//+------------------------------------------------------------------+
static void CMLPBase::MLPGradBatch(CMultilayerPerceptron &network,
                                   CMatrixDouble &xy,const int ssize,
                                   double &e,double &grad[])
  {
//--- create variables
   int i=0;
   int nin=0;
   int nout=0;
   int wcount=0;
//--- function call
   MLPProperties(network,nin,nout,wcount);
//--- initialization
   for(i=0;i<=wcount-1;i++)
      grad[i]=0;
   e=0;
   i=0;
//--- calculation
   while(i<=ssize-1)
     {
      MLPChunkedGradient(network,xy,i,MathMin(ssize,i+m_chunksize)-i,e,grad,false);
      i=i+m_chunksize;
     }
  }
//+------------------------------------------------------------------+
//| Batch gradient calculation for a set of inputs/outputs           |
//| (natural error function is used)                                 |
//| INPUT PARAMETERS:                                                |
//|     Network -   network initialized with one of the network      |
//|                 creation funcs                                   |
//|     XY      -   set of inputs/outputs; one sample=one row;       |
//|                 first NIn columns contain inputs,                |
//|                 next NOut columns - desired outputs.             |
//|     SSize   -   number of elements in XY                         |
//|     Grad    -   possibly preallocated array. If size of array is |
//|                 smaller than WCount, it will be reallocated. It  |
//|                 is recommended to reuse previously allocated     |
//|                 array to reduce allocation overhead.             |
//| OUTPUT PARAMETERS:                                               |
//|     E       -   error function, sum-of-squares for regression    |
//|                 networks, cross-entropy for classification       |
//|                 networks.                                        |
//|     Grad    -   gradient of E with respect to weights of network,|
//|                 array[WCount]                                    |
//+------------------------------------------------------------------+
static void CMLPBase::MLPGradNBatch(CMultilayerPerceptron &network,
                                    CMatrixDouble &xy,const int ssize,
                                    double &e,double &grad[])
  {
//--- create variables
   int i=0;
   int nin=0;
   int nout=0;
   int wcount=0;
//--- function call
   MLPProperties(network,nin,nout,wcount);
//--- initialization
   for(i=0;i<=wcount-1;i++)
      grad[i]=0;
   e=0;
   i=0;
//--- calculation
   while(i<=ssize-1)
     {
      MLPChunkedGradient(network,xy,i,MathMin(ssize,i+m_chunksize)-i,e,grad,true);
      i=i+m_chunksize;
     }
  }
//+------------------------------------------------------------------+
//| Batch Hessian calculation (natural error function) using         |
//| R-algorithm. Internal subroutine.                                |
//|      Hessian calculation based on R-algorithm described in       |
//|      "Fast Exact Multiplication by the Hessian",                 |
//|      B. A. Pearlmutter,                                          |
//|      Neural Computation, 1994.                                   |
//+------------------------------------------------------------------+
static void CMLPBase::MLPHessianNBatch(CMultilayerPerceptron &network,
                                       CMatrixDouble &xy,const int ssize,
                                       double &e,double &grad[],
                                       CMatrixDouble &h)
  {
//--- initialization
   e=0;
//--- function call
   MLPHessianBatchInternal(network,xy,ssize,true,e,grad,h);
  }
//+------------------------------------------------------------------+
//| Batch Hessian calculation using R-algorithm.                     |
//| Internal subroutine.                                             |
//|      Hessian calculation based on R-algorithm described in       |
//|      "Fast Exact Multiplication by the Hessian",                 |
//|      B. A. Pearlmutter,                                          |
//|      Neural Computation, 1994.                                   |
//+------------------------------------------------------------------+
static void CMLPBase::MLPHessianBatch(CMultilayerPerceptron &network,
                                      CMatrixDouble &xy,const int ssize,
                                      double &e,double &grad[],
                                      CMatrixDouble &h)
  {
//--- initialization
   e=0;
//--- function call
   MLPHessianBatchInternal(network,xy,ssize,false,e,grad,h);
  }
//+------------------------------------------------------------------+
//| Internal subroutine, shouldn't be called by user.                |
//+------------------------------------------------------------------+
static void CMLPBase::MLPInternalProcessVector(int &structinfo[],double &weights[],
                                               double &columnmeans[],
                                               double &columnsigmas[],
                                               double &neurons[],
                                               double &dfdnet[],
                                               double &x[],double &y[])
  {
//--- create variables
   int    i=0;
   int    n1=0;
   int    n2=0;
   int    w1=0;
   int    w2=0;
   int    ntotal=0;
   int    nin=0;
   int    nout=0;
   int    istart=0;
   int    offs=0;
   double net=0;
   double f=0;
   double df=0;
   double d2f=0;
   double mx=0;
   bool   perr;
   int    i_=0;
   int    i1_=0;
//--- Read network geometry
   nin=structinfo[1];
   nout=structinfo[2];
   ntotal=structinfo[3];
   istart=structinfo[5];
//--- Inputs standartisation and putting in the network
   for(i=0;i<=nin-1;i++)
     {
      //--- check
      if(columnsigmas[i]!=0.0)
         neurons[i]=(x[i]-columnmeans[i])/columnsigmas[i];
      else
         neurons[i]=x[i]-columnmeans[i];
     }
//--- Process network
   for(i=0;i<=ntotal-1;i++)
     {
      offs=istart+i*m_nfieldwidth;
      //--- check
      if(structinfo[offs+0]>0 || structinfo[offs+0]==-5)
        {
         //--- Activation function
         MLPActivationFunction(neurons[structinfo[offs+2]],structinfo[offs+0],f,df,d2f);
         //--- change values
         neurons[i]=f;
         dfdnet[i]=df;
         continue;
        }
      //--- check
      if(structinfo[offs+0]==0)
        {
         //--- Adaptive summator
         n1=structinfo[offs+2];
         n2=n1+structinfo[offs+1]-1;
         w1=structinfo[offs+3];
         w2=w1+structinfo[offs+1]-1;
         i1_=(n1)-(w1);
         net=0.0;
         //--- calculation
         for(i_=w1;i_<=w2;i_++)
            net+=weights[i_]*neurons[i_+i1_];
         neurons[i]=net;
         dfdnet[i]=1.0;
         continue;
        }
      //--- check
      if(structinfo[offs+0]<0)
        {
         perr=true;
         //--- check
         if(structinfo[offs+0]==-2)
           {
            //--- input neuron,left unchanged
            perr=false;
           }
         //--- check
         if(structinfo[offs+0]==-3)
           {
            //--- "-1" neuron
            neurons[i]=-1;
            perr=false;
           }
         //--- check
         if(structinfo[offs+0]==-4)
           {
            //--- "0" neuron
            neurons[i]=0;
            perr=false;
           }
         //--- check
         if(!CAp::Assert(!perr,__FUNCTION__+": internal error - unknown neuron type!"))
            return;
         continue;
        }
     }
//--- Extract result
   i1_=ntotal-nout;
   for(i_=0;i_<=nout-1;i_++)
      y[i_]=neurons[i_+i1_];
//--- Softmax post-processing or standardisation if needed
   if(!CAp::Assert(structinfo[6]==0 || structinfo[6]==1,__FUNCTION__+": unknown normalization type!"))
      return;
//--- check
   if(structinfo[6]==1)
     {
      //--- Softmax
      mx=y[0];
      for(i=1;i<=nout-1;i++)
         mx=MathMax(mx,y[i]);
      //--- calculation
      net=0;
      for(i=0;i<=nout-1;i++)
        {
         y[i]=MathExp(y[i]-mx);
         net=net+y[i];
        }
      for(i=0;i<=nout-1;i++)
         y[i]=y[i]/net;
     }
   else
     {
      //--- Standardisation
      for(i=0;i<=nout-1;i++)
         y[i]=y[i]*columnsigmas[nin+i]+columnmeans[nin+i];
     }
  }
//+------------------------------------------------------------------+
//| Serializer: allocation                                           |
//+------------------------------------------------------------------+
static void CMLPBase::MLPAlloc(CSerializer &s,CMultilayerPerceptron &network)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    fkind=0;
   double threshold=0;
   double v0=0;
   double v1=0;
   int    nin=0;
   int    nout=0;
//--- initialization
   nin=network.m_hllayersizes[0];
   nout=network.m_hllayersizes[CAp::Len(network.m_hllayersizes)-1];
//--- preparation to serialize
   s.Alloc_Entry();
   s.Alloc_Entry();
   s.Alloc_Entry();
//--- function call
   CApServ::AllocIntegerArray(s,network.m_hllayersizes,-1);
   for(i=1;i<=CAp::Len(network.m_hllayersizes)-1;i++)
     {
      for(j=0;j<=network.m_hllayersizes[i]-1;j++)
        {
         //--- function call
         MLPGetNeuronInfo(network,i,j,fkind,threshold);
         //--- preparation to serialize
         s.Alloc_Entry();
         s.Alloc_Entry();
         for(k=0;k<=network.m_hllayersizes[i-1]-1;k++)
            s.Alloc_Entry();
        }
     }
   for(j=0;j<=nin-1;j++)
     {
      //--- function call
      MLPGetInputScaling(network,j,v0,v1);
      //--- preparation to serialize
      s.Alloc_Entry();
      s.Alloc_Entry();
     }
   for(j=0;j<=nout-1;j++)
     {
      //--- function call
      MLPGetOutputScaling(network,j,v0,v1);
      //--- preparation to serialize
      s.Alloc_Entry();
      s.Alloc_Entry();
     }
  }
//+------------------------------------------------------------------+
//| Serializer: serialization                                        |
//+------------------------------------------------------------------+
static void CMLPBase::MLPSerialize(CSerializer &s,CMultilayerPerceptron &network)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    fkind=0;
   double threshold=0;
   double v0=0;
   double v1=0;
   int    nin=0;
   int    nout=0;
//--- change values
   nin=network.m_hllayersizes[0];
   nout=network.m_hllayersizes[CAp::Len(network.m_hllayersizes)-1];
//--- serializetion
   s.Serialize_Int(CSCodes::GetMLPSerializationCode());
   s.Serialize_Int(m_mlpfirstversion);
   s.Serialize_Bool(MLPIsSoftMax(network));
//--- function call
   CApServ::SerializeIntegerArray(s,network.m_hllayersizes,-1);
   for(i=1;i<=CAp::Len(network.m_hllayersizes)-1;i++)
     {
      for(j=0;j<=network.m_hllayersizes[i]-1;j++)
        {
         //--- function call
         MLPGetNeuronInfo(network,i,j,fkind,threshold);
         //--- serializetion
         s.Serialize_Int(fkind);
         s.Serialize_Double(threshold);
         for(k=0;k<=network.m_hllayersizes[i-1]-1;k++)
            s.Serialize_Double(MLPGetWeight(network,i-1,k,i,j));
        }
     }
   for(j=0;j<=nin-1;j++)
     {
      //--- function call
      MLPGetInputScaling(network,j,v0,v1);
      //--- serializetion
      s.Serialize_Double(v0);
      s.Serialize_Double(v1);
     }
   for(j=0;j<=nout-1;j++)
     {
      //--- function call
      MLPGetOutputScaling(network,j,v0,v1);
      //--- serializetion
      s.Serialize_Double(v0);
      s.Serialize_Double(v1);
     }
  }
//+------------------------------------------------------------------+
//| Serializer: unserialization                                      |
//+------------------------------------------------------------------+
static void CMLPBase::MLPUnserialize(CSerializer &s,CMultilayerPerceptron &network)
  {
//--- create variables
   int    i0=0;
   int    i1=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    fkind=0;
   double threshold=0;
   double v0=0;
   double v1=0;
   int    nin=0;
   int    nout=0;
   bool   issoftmax;
//--- create array
   int    layersizes[];
//--- check correctness of header
   i0=s.Unserialize_Int();
//--- check
   if(!CAp::Assert(i0==CSCodes::GetMLPSerializationCode(),__FUNCTION__+": stream header corrupted"))
      return;
//--- unserializetion
   i1=s.Unserialize_Int();
//--- check
   if(!CAp::Assert(i1==m_mlpfirstversion,__FUNCTION__+": stream header corrupted"))
      return;
//--- Create network
   issoftmax=s.Unserialize_Bool();
//--- function call
   CApServ::UnserializeIntegerArray(s,layersizes);
//--- check
   if(!CAp::Assert((CAp::Len(layersizes)==2 || CAp::Len(layersizes)==3) || CAp::Len(layersizes)==4,__FUNCTION__+": too many hidden layers!"))
      return;
//--- change values
   nin=layersizes[0];
   nout=layersizes[CAp::Len(layersizes)-1];
//--- check
   if(CAp::Len(layersizes)==2)
     {
      //--- check
      if(issoftmax)
         MLPCreateC0(layersizes[0],layersizes[1],network);
      else
         MLPCreate0(layersizes[0],layersizes[1],network);
     }
//--- check
   if(CAp::Len(layersizes)==3)
     {
      //--- check
      if(issoftmax)
         MLPCreateC1(layersizes[0],layersizes[1],layersizes[2],network);
      else
         MLPCreate1(layersizes[0],layersizes[1],layersizes[2],network);
     }
//--- check
   if(CAp::Len(layersizes)==4)
     {
      //--- check
      if(issoftmax)
         MLPCreateC2(layersizes[0],layersizes[1],layersizes[2],layersizes[3],network);
      else
         MLPCreate2(layersizes[0],layersizes[1],layersizes[2],layersizes[3],network);
     }
//--- Load neurons and weights
   for(i=1;i<=CAp::Len(layersizes)-1;i++)
     {
      for(j=0;j<=layersizes[i]-1;j++)
        {
         //--- unserializetion
         fkind=s.Unserialize_Int();
         threshold=s.Unserialize_Double();
         //--- function call
         MLPSetNeuronInfo(network,i,j,fkind,threshold);
         //--- unserializetion
         for(k=0;k<=layersizes[i-1]-1;k++)
           {
            v0=s.Unserialize_Double();
            //--- function call
            MLPSetWeight(network,i-1,k,i,j,v0);
           }
        }
     }

//
//--- Load standartizator
//
   for(j=0;j<=nin-1;j++)
     {
      //--- unserializetion
      v0=s.Unserialize_Double();
      v1=s.Unserialize_Double();
      //--- function call
      MLPSetInputScaling(network,j,v0,v1);
     }
   for(j=0;j<=nout-1;j++)
     {
      //--- unserializetion
      v0=s.Unserialize_Double();
      v1=s.Unserialize_Double();
      //--- function call
      MLPSetOutputScaling(network,j,v0,v1);
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine: adding new input layer to network           |
//+------------------------------------------------------------------+
static void CMLPBase::AddInputLayer(const int ncount,int &lsizes[],
                                    int &ltypes[],int &lconnfirst[],
                                    int &lconnlast[],int &lastproc)
  {
//--- change values
   lsizes[0]=ncount;
   ltypes[0]=-2;
   lconnfirst[0]=0;
   lconnlast[0]=0;
   lastproc=0;
  }
//+------------------------------------------------------------------+
//| Internal subroutine: adding new summator layer to network        |
//+------------------------------------------------------------------+
static void CMLPBase::AddBiasedSummatorLayer(const int ncount,int &lsizes[],
                                             int &ltypes[],int &lconnfirst[],
                                             int &lconnlast[],int &lastproc)
  {
//--- change values
   lsizes[lastproc+1]=1;
   ltypes[lastproc+1]=-3;
   lconnfirst[lastproc+1]=0;
   lconnlast[lastproc+1]=0;
   lsizes[lastproc+2]=ncount;
   ltypes[lastproc+2]=0;
   lconnfirst[lastproc+2]=lastproc;
   lconnlast[lastproc+2]=lastproc+1;
   lastproc=lastproc+2;
  }
//+------------------------------------------------------------------+
//| Internal subroutine: adding new summator layer to network        |
//+------------------------------------------------------------------+
static void CMLPBase::AddActivationLayer(const int functype,int &lsizes[],
                                         int &ltypes[],int &lconnfirst[],
                                         int &lconnlast[],int &lastproc)
  {
//--- check
   if(!CAp::Assert(functype>0 || functype==-5,__FUNCTION__+": incorrect function type"))
      return;
//--- change values
   lsizes[lastproc+1]=lsizes[lastproc];
   ltypes[lastproc+1]=functype;
   lconnfirst[lastproc+1]=lastproc;
   lconnlast[lastproc+1]=lastproc;
   lastproc=lastproc+1;
  }
//+------------------------------------------------------------------+
//| Internal subroutine: adding new zero layer to network            |
//+------------------------------------------------------------------+
static void CMLPBase::AddZeroLayer(int &lsizes[],int &ltypes[],
                                   int &lconnfirst[],int &lconnlast[],
                                   int &lastproc)
  {
//--- change values
   lsizes[lastproc+1]=1;
   ltypes[lastproc+1]=-4;
   lconnfirst[lastproc+1]=0;
   lconnlast[lastproc+1]=0;
   lastproc=lastproc+1;
  }
//+------------------------------------------------------------------+
//| This routine adds input layer to the high-level description of   |
//| the network.                                                     |
//| It modifies Network.HLConnections and Network.HLNeurons and      |
//| assumes that these arrays have enough place to store data.       |
//| It accepts following parameters:                                 |
//|     Network     -   network                                      |
//|     ConnIdx     -   index of the first free entry in the         |
//|                     HLConnections                                |
//|     NeuroIdx    -   index of the first free entry in the         |
//|                     HLNeurons                                    |
//|     StructInfoIdx-  index of the first entry in the low level    |
//|                     description of the current layer (in the     |
//|                     StructInfo array)                            |
//|     NIn         -   number of inputs                             |
//| It modified Network and indices.                                 |
//+------------------------------------------------------------------+
static void CMLPBase::HLAddInputLayer(CMultilayerPerceptron &network,
                                      int &connidx,int &neuroidx,
                                      int &structinfoidx,int nin)
  {
//--- create variables
   int i=0;
   int offs=0;
//--- initialization
   offs=m_hlm_nfieldwidth*neuroidx;
//--- change values
   for(i=0;i<=nin-1;i++)
     {
      network.m_hlneurons[offs+0]=0;
      network.m_hlneurons[offs+1]=i;
      network.m_hlneurons[offs+2]=-1;
      network.m_hlneurons[offs+3]=-1;
      offs=offs+m_hlm_nfieldwidth;
     }
//--- change values
   neuroidx=neuroidx+nin;
   structinfoidx=structinfoidx+nin;
  }
//+------------------------------------------------------------------+
//| This routine adds output layer to the high-level description of  |
//| the network.                                                     |
//| It modifies Network.HLConnections and Network. HLNeurons and     |
//| assumes that these arrays have enough place to store data. It    |
//| accepts following parameters:                                    |
//|     Network     -   network                                      |
//|     ConnIdx     -   index of the first free entry in the         |
//|                     HLConnections                                |
//|     NeuroIdx    -   index of the first free entry in the         |
//|                     HLNeurons                                    |
//|     StructInfoIdx-  index of the first entry in the low level    |
//|                     description of the current layer (in the     |
//|                     StructInfo array)                            |
//|     WeightsIdx  -   index of the first entry in the Weights      |
//|                     array which corresponds to the current layer |
//|     K           -   current layer index                          |
//|     NPrev       -   number of neurons in the previous layer      |
//|     NOut        -   number of outputs                            |
//|     IsCls       -   is it classifier network?                    |
//|     IsLinear    -   is it network with linear output?            |
//| It modified Network and ConnIdx/NeuroIdx/StructInfoIdx/WeightsIdx|
//+------------------------------------------------------------------+
static void CMLPBase::HLAddOutputLayer(CMultilayerPerceptron &network,
                                       int &connidx,int &neuroidx,
                                       int &structinfoidx,int &weightsidx,
                                       const int k,const int nprev,
                                       const int nout,const bool iscls,
                                       const bool islinearout)
  {
//--- create variables
   int i=0;
   int j=0;
   int neurooffs=0;
   int connoffs=0;
//--- check
   if(!CAp::Assert((iscls && islinearout) || !iscls,__FUNCTION__+": internal error"))
      return;
//--- initialization
   neurooffs=m_hlm_nfieldwidth*neuroidx;
   connoffs=m_hlconm_nfieldwidth*connidx;
//--- check
   if(!iscls)
     {
      //--- Regression network
      for(i=0;i<=nout-1;i++)
        {
         //--- change values
         network.m_hlneurons[neurooffs+0]=k;
         network.m_hlneurons[neurooffs+1]=i;
         network.m_hlneurons[neurooffs+2]=structinfoidx+1+nout+i;
         network.m_hlneurons[neurooffs+3]=weightsidx+nprev+(nprev+1)*i;
         neurooffs=neurooffs+m_hlm_nfieldwidth;
        }
      for(i=0;i<=nprev-1;i++)
        {
         for(j=0;j<=nout-1;j++)
           {
            //--- change values
            network.m_hlconnections[connoffs+0]=k-1;
            network.m_hlconnections[connoffs+1]=i;
            network.m_hlconnections[connoffs+2]=k;
            network.m_hlconnections[connoffs+3]=j;
            network.m_hlconnections[connoffs+4]=weightsidx+i+j*(nprev+1);
            connoffs=connoffs+m_hlconm_nfieldwidth;
           }
        }
      //--- change values
      connidx=connidx+nprev*nout;
      neuroidx=neuroidx+nout;
      structinfoidx=structinfoidx+2*nout+1;
      weightsidx=weightsidx+nout*(nprev+1);
     }
   else
     {
      //--- Classification network
      for(i=0;i<=nout-2;i++)
        {
         //--- change values
         network.m_hlneurons[neurooffs+0]=k;
         network.m_hlneurons[neurooffs+1]=i;
         network.m_hlneurons[neurooffs+2]=-1;
         network.m_hlneurons[neurooffs+3]=weightsidx+nprev+(nprev+1)*i;
         neurooffs=neurooffs+m_hlm_nfieldwidth;
        }
      //--- change values
      network.m_hlneurons[neurooffs+0]=k;
      network.m_hlneurons[neurooffs+1]=i;
      network.m_hlneurons[neurooffs+2]=-1;
      network.m_hlneurons[neurooffs+3]=-1;
      for(i=0;i<=nprev-1;i++)
        {
         for(j=0;j<=nout-2;j++)
           {
            //--- change values
            network.m_hlconnections[connoffs+0]=k-1;
            network.m_hlconnections[connoffs+1]=i;
            network.m_hlconnections[connoffs+2]=k;
            network.m_hlconnections[connoffs+3]=j;
            network.m_hlconnections[connoffs+4]=weightsidx+i+j*(nprev+1);
            connoffs=connoffs+m_hlconm_nfieldwidth;
           }
        }
      //--- change values
      connidx=connidx+nprev*(nout-1);
      neuroidx=neuroidx+nout;
      structinfoidx=structinfoidx+nout+2;
      weightsidx=weightsidx+(nout-1)*(nprev+1);
     }
  }
//+------------------------------------------------------------------+
//| This routine adds hidden layer to the high-level description of  |
//| the network.                                                     |
//| It modifies Network.HLConnections and Network.HLNeurons and      |
//| assumes that these arrays have enough place to store data. It    | 
//| accepts following parameters:                                    |
//|     Network     -   network                                      |
//|     ConnIdx     -   index of the first free entry in the         |
//|                     HLConnections                                |
//|     NeuroIdx    -   index of the first free entry in the         |
//|                     HLNeurons                                    |
//|     StructInfoIdx-  index of the first entry in the low level    |
//|                     description of the current layer (in the     |
//|                     StructInfo array)                            |
//|     WeightsIdx  -   index of the first entry in the Weights      |
//|                     array which corresponds to the current layer |
//|     K           -   current layer index                          |
//|     NPrev       -   number of neurons in the previous layer      |
//|     NCur        -   number of neurons in the current layer       |
//| It modified Network and ConnIdx/NeuroIdx/StructInfoIdx/WeightsIdx|
//+------------------------------------------------------------------+
static void CMLPBase::HLAddHiddenLayer(CMultilayerPerceptron &network,
                                       int &connidx,int &neuroidx,
                                       int &structinfoidx,int &weightsidx,
                                       const int k,const int nprev,
                                       const int ncur)
  {
//--- create variables
   int i=0;
   int j=0;
   int neurooffs=0;
   int connoffs=0;
//--- change values
   neurooffs=m_hlm_nfieldwidth*neuroidx;
   connoffs=m_hlconm_nfieldwidth*connidx;
   for(i=0;i<=ncur-1;i++)
     {
      //--- change values
      network.m_hlneurons[neurooffs+0]=k;
      network.m_hlneurons[neurooffs+1]=i;
      network.m_hlneurons[neurooffs+2]=structinfoidx+1+ncur+i;
      network.m_hlneurons[neurooffs+3]=weightsidx+nprev+(nprev+1)*i;
      neurooffs=neurooffs+m_hlm_nfieldwidth;
     }
   for(i=0;i<=nprev-1;i++)
     {
      for(j=0;j<=ncur-1;j++)
        {
         //--- change values
         network.m_hlconnections[connoffs+0]=k-1;
         network.m_hlconnections[connoffs+1]=i;
         network.m_hlconnections[connoffs+2]=k;
         network.m_hlconnections[connoffs+3]=j;
         network.m_hlconnections[connoffs+4]=weightsidx+i+j*(nprev+1);
         connoffs=connoffs+m_hlconm_nfieldwidth;
        }
     }
//--- change values
   connidx=connidx+nprev*ncur;
   neuroidx=neuroidx+ncur;
   structinfoidx=structinfoidx+2*ncur+1;
   weightsidx=weightsidx+ncur*(nprev+1);
  }
//+------------------------------------------------------------------+
//| This function fills high level information about network created | 
//| using internal MLPCreate() function.                             |
//| This function does NOT examine StructInfo for low level          | 
//| information, it just expects that network has following          |
//| structure:                                                       |
//|     input neuron            \                                    |
//|     ...                      | input layer                       |
//|     input neuron            /                                    |
//|     "-1" neuron             \                                    |
//|     biased summator          |                                   |
//|     ...                      |                                   |
//|     biased summator          | hidden layer(s), if there are     |
//|     activation function      | exists any                        |
//|     ...                      |                                   |
//|     activation function     /                                    |
//|     "-1" neuron            \                                     |
//|     biased summator         | output layer:                      |
//|     ...                     | * we have NOut summators/activators|
//|     biased summator         |   for regression networks          |
//|     activation function     | * we have only NOut-1 summators and| 
//|     ...                     |   no activators for classifiers    |
//|     activation function     | * we have "0" neuron only when we  |
//|     "0" neuron              /   have classifier                  |
//+------------------------------------------------------------------+
static void CMLPBase::FillHighLevelInformation(CMultilayerPerceptron &network,
                                               const int nin,const int nhid1,
                                               const int nhid2,const int nout,
                                               const bool iscls,const bool islinearout)
  {
//--- create variables
   int idxweights=0;
   int idxstruct=0;
   int idxneuro=0;
   int idxconn=0;
//--- check
   if(!CAp::Assert((iscls && islinearout) || !iscls,__FUNCTION__+": internal error"))
      return;
//--- Preparations common to all types of networks
   idxweights=0;
   idxneuro=0;
   idxstruct=0;
   idxconn=0;
   network.m_hlnetworktype=0;
//--- network without hidden layers
   if(nhid1==0)
     {
      //--- allocation
      ArrayResizeAL(network.m_hllayersizes,2);
      //--- change values
      network.m_hllayersizes[0]=nin;
      network.m_hllayersizes[1]=nout;
      //--- check
      if(!iscls)
        {
         //--- allocation
         ArrayResizeAL(network.m_hlconnections,m_hlconm_nfieldwidth*nin*nout);
         ArrayResizeAL(network.m_hlneurons,m_hlm_nfieldwidth*(nin+nout));
         network.m_hlnormtype=0;
        }
      else
        {
         //--- allocation
         ArrayResizeAL(network.m_hlconnections,m_hlconm_nfieldwidth*nin*(nout-1));
         ArrayResizeAL(network.m_hlneurons,m_hlm_nfieldwidth*(nin+nout));
         network.m_hlnormtype=1;
        }
      //--- function call
      HLAddInputLayer(network,idxconn,idxneuro,idxstruct,nin);
      //--- function call
      HLAddOutputLayer(network,idxconn,idxneuro,idxstruct,idxweights,1,nin,nout,iscls,islinearout);
      //--- exit the function
      return;
     }
//--- network with one hidden layers
   if(nhid2==0)
     {
      //--- allocation
      ArrayResizeAL(network.m_hllayersizes,3);
      //--- change values
      network.m_hllayersizes[0]=nin;
      network.m_hllayersizes[1]=nhid1;
      network.m_hllayersizes[2]=nout;
      //--- check
      if(!iscls)
        {
         //--- allocation
         ArrayResizeAL(network.m_hlconnections,m_hlconm_nfieldwidth*(nin*nhid1+nhid1*nout));
         ArrayResizeAL(network.m_hlneurons,m_hlm_nfieldwidth*(nin+nhid1+nout));
         network.m_hlnormtype=0;
        }
      else
        {
         //--- allocation
         ArrayResizeAL(network.m_hlconnections,m_hlconm_nfieldwidth*(nin*nhid1+nhid1*(nout-1)));
         ArrayResizeAL(network.m_hlneurons,m_hlm_nfieldwidth*(nin+nhid1+nout));
         network.m_hlnormtype=1;
        }
      //--- function call
      HLAddInputLayer(network,idxconn,idxneuro,idxstruct,nin);
      //--- function call
      HLAddHiddenLayer(network,idxconn,idxneuro,idxstruct,idxweights,1,nin,nhid1);
      //--- function call
      HLAddOutputLayer(network,idxconn,idxneuro,idxstruct,idxweights,2,nhid1,nout,iscls,islinearout);
      //--- exit the function
      return;
     }
//--- Two hidden layers
   ArrayResizeAL(network.m_hllayersizes,4);
//--- change values
   network.m_hllayersizes[0]=nin;
   network.m_hllayersizes[1]=nhid1;
   network.m_hllayersizes[2]=nhid2;
   network.m_hllayersizes[3]=nout;
//--- check
   if(!iscls)
     {
      //--- allocation
      ArrayResizeAL(network.m_hlconnections,m_hlconm_nfieldwidth*(nin*nhid1+nhid1*nhid2+nhid2*nout));
      ArrayResizeAL(network.m_hlneurons,m_hlm_nfieldwidth*(nin+nhid1+nhid2+nout));
      network.m_hlnormtype=0;
     }
   else
     {
      //--- allocation
      ArrayResizeAL(network.m_hlconnections,m_hlconm_nfieldwidth*(nin*nhid1+nhid1*nhid2+nhid2*(nout-1)));
      ArrayResizeAL(network.m_hlneurons,m_hlm_nfieldwidth*(nin+nhid1+nhid2+nout));
      network.m_hlnormtype=1;
     }
//--- function call
   HLAddInputLayer(network,idxconn,idxneuro,idxstruct,nin);
//--- function call
   HLAddHiddenLayer(network,idxconn,idxneuro,idxstruct,idxweights,1,nin,nhid1);
//--- function call
   HLAddHiddenLayer(network,idxconn,idxneuro,idxstruct,idxweights,2,nhid1,nhid2);
//--- function call
   HLAddOutputLayer(network,idxconn,idxneuro,idxstruct,idxweights,3,nhid2,nout,iscls,islinearout);
  }
//+------------------------------------------------------------------+
//| Internal subroutine.                                             |
//+------------------------------------------------------------------+
static void CMLPBase::MLPCreate(const int nin,const int nout,int &lsizes[],
                                int &ltypes[],int &lconnfirst[],int &lconnlast[],
                                const int layerscount,const bool isclsnet,
                                CMultilayerPerceptron &network)
  {
//--- create variables
   int i=0;
   int j=0;
   int ssize=0;
   int ntotal=0;
   int wcount=0;
   int offs=0;
   int nprocessed=0;
   int wallocated=0;
//--- creating arrays
   int localtemp[];
   int lnfirst[];
   int lnsyn[];
//--- Check
   if(!CAp::Assert(layerscount>0,__FUNCTION__+": wrong parameters!"))
      return;
//--- check
   if(!CAp::Assert(ltypes[0]==-2,__FUNCTION__+": wrong LTypes[0] (must be -2)!"))
      return;
   for(i=0;i<=layerscount-1;i++)
     {
      //--- check
      if(!CAp::Assert(lsizes[i]>0,__FUNCTION__+": wrong LSizes!"))
         return;
      //--- check
      if(!CAp::Assert(lconnfirst[i]>=0 &&(lconnfirst[i]<i || i==0),__FUNCTION__+": wrong LConnFirst!"))
         return;
      //--- check
      if(!CAp::Assert(lconnlast[i]>=lconnfirst[i]&&(lconnlast[i]<i || i==0),__FUNCTION__+": wrong LConnLast!"))
         return;
     }
//--- Build network geometry
   ArrayResizeAL(lnfirst,layerscount);
   ArrayResizeAL(lnsyn,layerscount);
//--- initialization
   ntotal=0;
   wcount=0;
//--- calculation
   for(i=0;i<=layerscount-1;i++)
     {
      //--- Analyze connections.
      //--- This code must throw an assertion in case of unknown LTypes[I]
      lnsyn[i]=-1;
      //--- check
      if(ltypes[i]>=0 || ltypes[i]==-5)
        {
         lnsyn[i]=0;
         for(j=lconnfirst[i];j<=lconnlast[i];j++)
            lnsyn[i]=lnsyn[i]+lsizes[j];
        }
      else
        {
         //--- check
         if((ltypes[i]==-2 || ltypes[i]==-3) || ltypes[i]==-4)
            lnsyn[i]=0;
        }
      //--- check
      if(!CAp::Assert(lnsyn[i]>=0,__FUNCTION__+": internal error #0!"))
         return;
      //--- Other info
      lnfirst[i]=ntotal;
      ntotal=ntotal+lsizes[i];
      //--- check
      if(ltypes[i]==0)
         wcount=wcount+lnsyn[i]*lsizes[i];
     }
   ssize=7+ntotal*m_nfieldwidth;
//--- Allocate
   ArrayResizeAL(network.m_structinfo,ssize);
   ArrayResizeAL(network.m_weights,wcount);
//--- check
   if(isclsnet)
     {
      //--- allocation
      ArrayResizeAL(network.m_columnmeans,nin);
      ArrayResizeAL(network.m_columnsigmas,nin);
     }
   else
     {
      //--- allocation
      ArrayResizeAL(network.m_columnmeans,nin+nout);
      ArrayResizeAL(network.m_columnsigmas,nin+nout);
     }
//--- allocation
   ArrayResizeAL(network.m_neurons,ntotal);
   network.m_chunks.Resize(3*ntotal+1,m_chunksize);
   ArrayResizeAL(network.m_nwbuf,MathMax(wcount,2*nout));
   ArrayResizeAL(network.m_integerbuf,4);
   ArrayResizeAL(network.m_dfdnet,ntotal);
   ArrayResizeAL(network.m_x,nin);
   ArrayResizeAL(network.m_y,nout);
   ArrayResizeAL(network.m_derror,ntotal);
//--- Fill structure: global info
   network.m_structinfo[0]=ssize;
   network.m_structinfo[1]=nin;
   network.m_structinfo[2]=nout;
   network.m_structinfo[3]=ntotal;
   network.m_structinfo[4]=wcount;
   network.m_structinfo[5]=7;
//--- check
   if(isclsnet)
      network.m_structinfo[6]=1;
   else
      network.m_structinfo[6]=0;
//--- Fill structure: neuron connections
   nprocessed=0;
   wallocated=0;
//--- calculation
   for(i=0;i<=layerscount-1;i++)
     {
      for(j=0;j<=lsizes[i]-1;j++)
        {
         offs=network.m_structinfo[5]+nprocessed*m_nfieldwidth;
         network.m_structinfo[offs+0]=ltypes[i];
         //--- check
         if(ltypes[i]==0)
           {
            //--- Adaptive summator:
            //--- * connections with weights to previous neurons
            network.m_structinfo[offs+1]=lnsyn[i];
            network.m_structinfo[offs+2]=lnfirst[lconnfirst[i]];
            network.m_structinfo[offs+3]=wallocated;
            wallocated=wallocated+lnsyn[i];
            nprocessed=nprocessed+1;
           }
         //--- check
         if(ltypes[i]>0 || ltypes[i]==-5)
           {
            //--- Activation layer:
            //--- * each neuron connected to one (only one) of previous neurons.
            //--- * no weights
            network.m_structinfo[offs+1]=1;
            network.m_structinfo[offs+2]=lnfirst[lconnfirst[i]]+j;
            network.m_structinfo[offs+3]=-1;
            nprocessed=nprocessed+1;
           }
         //--- check
         if((ltypes[i]==-2 || ltypes[i]==-3) || ltypes[i]==-4)
            nprocessed=nprocessed+1;
        }
     }
//--- check
   if(!CAp::Assert(wallocated==wcount,__FUNCTION__+": internal error #1!"))
      return;
//--- check
   if(!CAp::Assert(nprocessed==ntotal,__FUNCTION__+": internal error #2!"))
      return;
//--- Fill weights by small random values
//--- Initialize means and sigmas
   for(i=0;i<=wcount-1;i++)
      network.m_weights[i]=CMath::RandomReal()-0.5;
   for(i=0;i<=nin-1;i++)
     {
      network.m_columnmeans[i]=0;
      network.m_columnsigmas[i]=1;
     }
//--- check
   if(!isclsnet)
     {
      for(i=0;i<=nout-1;i++)
        {
         network.m_columnmeans[nin+i]=0;
         network.m_columnsigmas[nin+i]=1;
        }
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine for Hessian calculation.                     |
//| WARNING!!! Unspeakable math far beyong human capabilities :)     |
//+------------------------------------------------------------------+
static void CMLPBase::MLPHessianBatchInternal(CMultilayerPerceptron &network,
                                              CMatrixDouble &xy,const int ssize,
                                              const bool naturalerr,double &e,
                                              double &grad[],CMatrixDouble &h)
  {
//--- create variables
   int    nin=0;
   int    nout=0;
   int    wcount=0;
   int    ntotal=0;
   int    istart=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    kl=0;
   int    offs=0;
   int    n1=0;
   int    n2=0;
   int    w1=0;
   int    w2=0;
   double s=0;
   double t=0;
   double v=0;
   double et=0;
   bool   bflag;
   double f=0;
   double df=0;
   double d2f=0;
   double deidyj=0;
   double mx=0;
   double q=0;
   double z=0;
   double s2=0;
   double expi=0;
   double expj=0;
   int    i_=0;
   int    i1_=0;
//--- creating arrays
   double x[];
   double desiredy[];
   double gt[];
   double zeros[];
//--- create matrix
   CMatrixDouble rx;
   CMatrixDouble ry;
   CMatrixDouble rdx;
   CMatrixDouble rdy;
//--- initialization
   e=0;
//--- function call
   MLPProperties(network,nin,nout,wcount);
//--- initialization
   ntotal=network.m_structinfo[3];
   istart=network.m_structinfo[5];
//--- Prepare
   ArrayResizeAL(x,nin);
   ArrayResizeAL(desiredy,nout);
   ArrayResizeAL(zeros,wcount);
   ArrayResizeAL(gt,wcount);
   rx.Resize(ntotal+nout,wcount);
   ry.Resize(ntotal+nout,wcount);
   rdx.Resize(ntotal+nout,wcount);
   rdy.Resize(ntotal+nout,wcount);
//--- initialization
   e=0;
   for(i=0;i<=wcount-1;i++)
      zeros[i]=0;
   for(i_=0;i_<=wcount-1;i_++)
      grad[i_]=zeros[i_];
   for(i=0;i<=wcount-1;i++)
     {
      for(i_=0;i_<=wcount-1;i_++)
         h[i].Set(i_,zeros[i_]);
     }
//--- Process
   for(k=0;k<=ssize-1;k++)
     {
      //--- Process vector with MLPGradN.
      //--- Now Neurons,DFDNET and DError contains results of the last run.
      for(i_=0;i_<=nin-1;i_++)
         x[i_]=xy[k][i_];
      //--- check
      if(MLPIsSoftMax(network))
        {
         //--- class labels outputs
         kl=(int)MathRound(xy[k][nin]);
         for(i=0;i<=nout-1;i++)
           {
            //--- check
            if(i==kl)
               desiredy[i]=1;
            else
               desiredy[i]=0;
           }
        }
      else
        {
         //--- real outputs
         i1_=nin;
         for(i_=0;i_<=nout-1;i_++)
            desiredy[i_]=xy[k][i_+i1_];
        }
      //--- check
      if(naturalerr)
         MLPGradN(network,x,desiredy,et,gt);
      else
         MLPGrad(network,x,desiredy,et,gt);
      //--- grad,error
      e=e+et;
      for(i_=0;i_<=wcount-1;i_++)
         grad[i_]=grad[i_]+gt[i_];
      //--- Hessian.
      //--- Forward pass of the R-algorithm
      for(i=0;i<=ntotal-1;i++)
        {
         offs=istart+i*m_nfieldwidth;
         for(i_=0;i_<=wcount-1;i_++)
            rx[i].Set(i_,zeros[i_]);
         for(i_=0;i_<=wcount-1;i_++)
            ry[i].Set(i_,zeros[i_]);
         //--- check
         if(network.m_structinfo[offs+0]>0 || network.m_structinfo[offs+0]==-5)
           {
            //--- Activation function
            n1=network.m_structinfo[offs+2];
            for(i_=0;i_<=wcount-1;i_++)
               rx[i].Set(i_,ry[n1][i_]);
            //--- calculation
            v=network.m_dfdnet[i];
            for(i_=0;i_<=wcount-1;i_++)
               ry[i].Set(i_,v*rx[i][i_]);
            continue;
           }
         //--- check
         if(network.m_structinfo[offs+0]==0)
           {
            //--- Adaptive summator
            n1=network.m_structinfo[offs+2];
            n2=n1+network.m_structinfo[offs+1]-1;
            w1=network.m_structinfo[offs+3];
            w2=w1+network.m_structinfo[offs+1]-1;
            //--- calculation
            for(j=n1;j<=n2;j++)
              {
               v=network.m_weights[w1+j-n1];
               for(i_=0;i_<=wcount-1;i_++)
                  rx[i].Set(i_,rx[i][i_]+v*ry[j][i_]);
               rx[i].Set(w1+j-n1,rx[i][w1+j-n1]+network.m_neurons[j]);
              }
            for(i_=0;i_<=wcount-1;i_++)
               ry[i].Set(i_,rx[i][i_]);
            continue;
           }
         //--- check
         if(network.m_structinfo[offs+0]<0)
           {
            bflag=true;
            //--- check
            if(network.m_structinfo[offs+0]==-2)
              {
               //--- input neuron,left unchanged
               bflag=false;
              }
            //--- check
            if(network.m_structinfo[offs+0]==-3)
              {
               //--- "-1" neuron,left unchanged
               bflag=false;
              }
            //--- check
            if(network.m_structinfo[offs+0]==-4)
              {
               //--- "0" neuron,left unchanged
               bflag=false;
              }
            //--- check
            if(!CAp::Assert(!bflag,__FUNCTION__+": internal error - unknown neuron type!"))
               return;
            continue;
           }
        }
      //--- Hessian. Backward pass of the R-algorithm.
      //--- Stage 1. Initialize RDY
      for(i=0;i<=ntotal+nout-1;i++)
        {
         for(i_=0;i_<=wcount-1;i_++)
            rdy[i].Set(i_,zeros[i_]);
        }
      //--- check
      if(network.m_structinfo[6]==0)
        {
         //--- Standardisation.
         //--- In context of the Hessian calculation standardisation
         //--- is considered as additional layer with weightless
         //--- activation function:
         //--- F(NET) :=Sigma*NET
         //--- So we add one more layer to forward pass,and
         //--- make forward/backward pass through this layer.
         for(i=0;i<=nout-1;i++)
           {
            n1=ntotal-nout+i;
            n2=ntotal+i;
            //--- Forward pass from N1 to N2
            for(i_=0;i_<=wcount-1;i_++)
               rx[n2].Set(i_,ry[n1][i_]);
            v=network.m_columnsigmas[nin+i];
            for(i_=0;i_<=wcount-1;i_++)
               ry[n2].Set(i_,v*rx[n2][i_]);
            //--- Initialization of RDY
            for(i_=0;i_<=wcount-1;i_++)
               rdy[n2].Set(i_,ry[n2][i_]);
            //--- Backward pass from N2 to N1:
            //--- 1. Calculate R(dE/dX).
            //--- 2. No R(dE/dWij) is needed since weight of activation neuron
            //---    is fixed to 1. So we can update R(dE/dY) for
            //---    the connected neuron (note that Vij=0,Wij=1)
            df=network.m_columnsigmas[nin+i];
            for(i_=0;i_<=wcount-1;i_++)
               rdx[n2].Set(i_,df*rdy[n2][i_]);
            for(i_=0;i_<=wcount-1;i_++)
               rdy[n1].Set(i_,rdy[n1][i_]+rdx[n2][i_]);
           }
        }
      else
        {
         //--- Softmax.
         //--- Initialize RDY using generalized expression for ei'(yi)
         //--- (see expression (9) from p. 5 of "Fast Exact Multiplication by the Hessian").
         //--- When we are working with softmax network,generalized
         //--- expression for ei'(yi) is used because softmax
         //--- normalization leads to ei,which depends on all y's
         if(naturalerr)
           {
            //--- softmax + cross-entropy.
            //--- We have:
            //--- S=sum(exp(yk)),
            //--- ei=sum(trn)*exp(yi)/S-trn_i
            //--- j=i:   d(ei)/d(yj)=T*exp(yi)*(S-exp(yi))/S^2
            //--- j<>i:  d(ei)/d(yj)=-T*exp(yi)*exp(yj)/S^2
            t=0;
            for(i=0;i<=nout-1;i++)
               t=t+desiredy[i];
            mx=network.m_neurons[ntotal-nout];
            //--- calculation
            for(i=0;i<=nout-1;i++)
               mx=MathMax(mx,network.m_neurons[ntotal-nout+i]);
            s=0;
            for(i=0;i<=nout-1;i++)
              {
               network.m_nwbuf[i]=MathExp(network.m_neurons[ntotal-nout+i]-mx);
               s=s+network.m_nwbuf[i];
              }
            //--- calculation
            for(i=0;i<=nout-1;i++)
              {
               for(j=0;j<=nout-1;j++)
                 {
                  //--- check
                  if(j==i)
                    {
                     deidyj=t*network.m_nwbuf[i]*(s-network.m_nwbuf[i])/CMath::Sqr(s);
                     for(i_=0;i_<=wcount-1;i_++)
                        rdy[ntotal-nout+i].Set(i_,rdy[ntotal-nout+i][i_]+deidyj*ry[ntotal-nout+i][i_]);
                    }
                  else
                    {
                     deidyj=-(t*network.m_nwbuf[i]*network.m_nwbuf[j]/CMath::Sqr(s));
                     for(i_=0;i_<=wcount-1;i_++)
                        rdy[ntotal-nout+i].Set(i_,rdy[ntotal-nout+i][i_]+deidyj*ry[ntotal-nout+j][i_]);
                    }
                 }
              }
           }
         else
           {
            //--- For a softmax + squared error we have expression
            //--- far beyond human imagination so we dont even try
            //--- to comment on it. Just enjoy the code...
            //--- P.S. That's why "natural error" is called "natural" -
            //--- compact beatiful expressions,fast code....
            mx=network.m_neurons[ntotal-nout];
            for(i=0;i<=nout-1;i++)
               mx=MathMax(mx,network.m_neurons[ntotal-nout+i]);
            //--- calculation
            s=0;
            s2=0;
            for(i=0;i<=nout-1;i++)
              {
               network.m_nwbuf[i]=MathExp(network.m_neurons[ntotal-nout+i]-mx);
               s=s+network.m_nwbuf[i];
               s2=s2+CMath::Sqr(network.m_nwbuf[i]);
              }
            //--- calculation
            q=0;
            for(i=0;i<=nout-1;i++)
               q=q+(network.m_y[i]-desiredy[i])*network.m_nwbuf[i];
            for(i=0;i<=nout-1;i++)
              {
               //--- change values
               z=-q+(network.m_y[i]-desiredy[i])*s;
               expi=network.m_nwbuf[i];
               for(j=0;j<=nout-1;j++)
                 {
                  expj=network.m_nwbuf[j];
                  //--- check
                  if(j==i)
                     deidyj=expi/CMath::Sqr(s)*((z+expi)*(s-2*expi)/s+expi*s2/CMath::Sqr(s));
                  else
                     deidyj=expi*expj/CMath::Sqr(s)*(s2/CMath::Sqr(s)-2*z/s-(expi+expj)/s+(network.m_y[i]-desiredy[i])-(network.m_y[j]-desiredy[j]));
                  for(i_=0;i_<=wcount-1;i_++)
                     rdy[ntotal-nout+i].Set(i_,rdy[ntotal-nout+i][i_]+deidyj*ry[ntotal-nout+j][i_]);
                 }
              }
           }
        }
      //--- Hessian. Backward pass of the R-algorithm
      //--- Stage 2. Process.
      for(i=ntotal-1;i>=0;i--)
        {
         //--- Possible variants:
         //--- 1. Activation function
         //--- 2. Adaptive summator
         //--- 3. Special neuron
         offs=istart+i*m_nfieldwidth;
         //--- check
         if(network.m_structinfo[offs+0]>0 || network.m_structinfo[offs+0]==-5)
           {
            n1=network.m_structinfo[offs+2];
            //--- First,calculate R(dE/dX).
            MLPActivationFunction(network.m_neurons[n1],network.m_structinfo[offs+0],f,df,d2f);
            v=d2f*network.m_derror[i];
            for(i_=0;i_<=wcount-1;i_++)
               rdx[i].Set(i_,df*rdy[i][i_]);
            for(i_=0;i_<=wcount-1;i_++)
               rdx[i].Set(i_,rdx[i][i_]+v*rx[i][i_]);
            //--- No R(dE/dWij) is needed since weight of activation neuron
            //--- is fixed to 1.
            //--- So we can update R(dE/dY) for the connected neuron.
            //--- (note that Vij=0,Wij=1)
            for(i_=0;i_<=wcount-1;i_++)
               rdy[n1].Set(i_,rdy[n1][i_]+rdx[i][i_]);
            continue;
           }
         //--- check
         if(network.m_structinfo[offs+0]==0)
           {
            //--- Adaptive summator
            n1=network.m_structinfo[offs+2];
            n2=n1+network.m_structinfo[offs+1]-1;
            w1=network.m_structinfo[offs+3];
            w2=w1+network.m_structinfo[offs+1]-1;
            //--- First,calculate R(dE/dX).
            for(i_=0;i_<=wcount-1;i_++)
               rdx[i].Set(i_,rdy[i][i_]);
            //--- Then,calculate R(dE/dWij)
            for(j=w1;j<=w2;j++)
              {
               v=network.m_neurons[n1+j-w1];
               for(i_=0;i_<=wcount-1;i_++)
                  h[j].Set(i_,h[j][i_]+v*rdx[i][i_]);
               //--- calculation
               v=network.m_derror[i];
               for(i_=0;i_<=wcount-1;i_++)
                  h[j].Set(i_,h[j][i_]+v*ry[n1+j-w1][i_]);
              }
            //--- And finally,update R(dE/dY) for connected neurons.
            for(j=w1;j<=w2;j++)
              {
               v=network.m_weights[j];
               for(i_=0;i_<=wcount-1;i_++)
                  rdy[n1+j-w1].Set(i_,rdy[n1+j-w1][i_]+v*rdx[i][i_]);
               rdy[n1+j-w1].Set(j,rdy[n1+j-w1][j]+network.m_derror[i]);
              }
            continue;
           }
         //--- check
         if(network.m_structinfo[offs+0]<0)
           {
            bflag=false;
            //--- check
            if((network.m_structinfo[offs+0]==-2 || network.m_structinfo[offs+0]==-3) || network.m_structinfo[offs+0]==-4)
              {
               //--- Special neuron type,no back-propagation required
               bflag=true;
              }
            //--- check
            if(!CAp::Assert(bflag,__FUNCTION__+": unknown neuron type!"))
               return;
            continue;
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//| Network must be processed by MLPProcess on X                     |
//+------------------------------------------------------------------+
static void CMLPBase::MLPInternalCalculateGradient(CMultilayerPerceptron &network,
                                                   double &neurons[],
                                                   double &weights[],
                                                   double &derror[],
                                                   double &grad[],
                                                   const bool naturalerrorfunc)
  {
//--- create variables
   int    i=0;
   int    n1=0;
   int    n2=0;
   int    w1=0;
   int    w2=0;
   int    ntotal=0;
   int    istart=0;
   int    nin=0;
   int    nout=0;
   int    offs=0;
   double dedf=0;
   double dfdnet=0;
   double v=0;
   double fown=0;
   double deown=0;
   double net=0;
   double mx=0;
   bool   bflag;
   int    i_=0;
   int    i1_=0;
//--- Read network geometry
   nin=network.m_structinfo[1];
   nout=network.m_structinfo[2];
   ntotal=network.m_structinfo[3];
   istart=network.m_structinfo[5];
//--- Pre-processing of dError/dOut:
//--- from dError/dOut(normalized) to dError/dOut(non-normalized)
   if(!CAp::Assert(network.m_structinfo[6]==0 || network.m_structinfo[6]==1,__FUNCTION__+": unknown normalization type!"))
      return;
//--- check
   if(network.m_structinfo[6]==1)
     {
      //--- Softmax
      if(!naturalerrorfunc)
        {
         mx=network.m_neurons[ntotal-nout];
         for(i=0;i<=nout-1;i++)
            mx=MathMax(mx,network.m_neurons[ntotal-nout+i]);
         net=0;
         for(i=0;i<=nout-1;i++)
           {
            network.m_nwbuf[i]=MathExp(network.m_neurons[ntotal-nout+i]-mx);
            net=net+network.m_nwbuf[i];
           }
         //--- calculation
         i1_=-(ntotal-nout);
         v=0.0;
         for(i_=ntotal-nout;i_<=ntotal-1;i_++)
            v+=network.m_derror[i_]*network.m_nwbuf[i_+i1_];
         for(i=0;i<=nout-1;i++)
           {
            fown=network.m_nwbuf[i];
            deown=network.m_derror[ntotal-nout+i];
            network.m_nwbuf[nout+i]=(-v+deown*fown+deown*(net-fown))*fown/CMath::Sqr(net);
           }
         for(i=0;i<=nout-1;i++)
            network.m_derror[ntotal-nout+i]=network.m_nwbuf[nout+i];
        }
     }
   else
     {
      //--- Un-standardisation
      for(i=0;i<=nout-1;i++)
         network.m_derror[ntotal-nout+i]=network.m_derror[ntotal-nout+i]*network.m_columnsigmas[nin+i];
     }
//--- Backpropagation
   for(i=ntotal-1;i>=0;i--)
     {
      //--- Extract info
      offs=istart+i*m_nfieldwidth;
      //--- check
      if(network.m_structinfo[offs+0]>0 || network.m_structinfo[offs+0]==-5)
        {
         //--- Activation function
         dedf=network.m_derror[i];
         dfdnet=network.m_dfdnet[i];
         derror[network.m_structinfo[offs+2]]=derror[network.m_structinfo[offs+2]]+dedf*dfdnet;
         continue;
        }
      //--- check
      if(network.m_structinfo[offs+0]==0)
        {
         //--- Adaptive summator
         n1=network.m_structinfo[offs+2];
         n2=n1+network.m_structinfo[offs+1]-1;
         w1=network.m_structinfo[offs+3];
         w2=w1+network.m_structinfo[offs+1]-1;
         dedf=network.m_derror[i];
         dfdnet=1.0;
         v=dedf*dfdnet;
         i1_=n1-w1;
         //--- calculation
         for(i_=w1;i_<=w2;i_++)
            grad[i_]=v*neurons[i_+i1_];
         i1_=w1-n1;
         for(i_=n1;i_<=n2;i_++)
            derror[i_]=derror[i_]+v*weights[i_+i1_];
         continue;
        }
      //--- check
      if(network.m_structinfo[offs+0]<0)
        {
         bflag=false;
         //--- check
         if((network.m_structinfo[offs+0]==-2 || network.m_structinfo[offs+0]==-3) || network.m_structinfo[offs+0]==-4)
           {
            //--- Special neuron type,no back-propagation required
            bflag=true;
           }
         //--- check
         if(!CAp::Assert(bflag,__FUNCTION__+": unknown neuron type!"))
            return;
         continue;
        }
     }
  }
//+------------------------------------------------------------------+
//| Internal subroutine, chunked gradient                            |
//+------------------------------------------------------------------+
static void CMLPBase::MLPChunkedGradient(CMultilayerPerceptron &network,
                                         CMatrixDouble &xy,const int cstart,
                                         const int csize,double &e,
                                         double &grad[],const bool naturalerrorfunc)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    kl=0;
   int    n1=0;
   int    n2=0;
   int    w1=0;
   int    w2=0;
   int    c1=0;
   int    c2=0;
   int    ntotal=0;
   int    nin=0;
   int    nout=0;
   int    offs=0;
   double f=0;
   double df=0;
   double d2f=0;
   double v=0;
   double s=0;
   double fown=0;
   double deown=0;
   double net=0;
   double lnnet=0;
   double mx=0;
   bool   bflag;
   int    istart=0;
   int    ineurons=0;
   int    idfdnet=0;
   int    iderror=0;
   int    izeros=0;
   int    i_=0;
   int    i1_=0;
//--- Read network geometry,prepare data
   nin=network.m_structinfo[1];
   nout=network.m_structinfo[2];
   ntotal=network.m_structinfo[3];
   istart=network.m_structinfo[5];
   c1=cstart;
   c2=cstart+csize-1;
   ineurons=0;
   idfdnet=ntotal;
   iderror=2*ntotal;
   izeros=3*ntotal;
   for(j=0;j<=csize-1;j++)
      network.m_chunks[izeros].Set(j,0);
//--- Forward pass:
//--- 1. Load inputs from XY to Chunks[0:NIn-1,0:CSize-1]
//--- 2. Forward pass
   for(i=0;i<=nin-1;i++)
     {
      for(j=0;j<=csize-1;j++)
        {
         //--- check
         if(network.m_columnsigmas[i]!=0.0)
            network.m_chunks[i].Set(j,(xy[c1+j][i]-network.m_columnmeans[i])/network.m_columnsigmas[i]);
         else
            network.m_chunks[i].Set(j,xy[c1+j][i]-network.m_columnmeans[i]);
        }
     }
   for(i=0;i<=ntotal-1;i++)
     {
      offs=istart+i*m_nfieldwidth;
      //--- check
      if(network.m_structinfo[offs+0]>0 || network.m_structinfo[offs+0]==-5)
        {
         //--- Activation function:
         //--- * calculate F vector,F(i)=F(NET(i))
         n1=network.m_structinfo[offs+2];
         for(i_=0;i_<=csize-1;i_++)
            network.m_chunks[i].Set(i_,network.m_chunks[n1][i_]);
         for(j=0;j<=csize-1;j++)
           {
            //--- function call
            MLPActivationFunction(network.m_chunks[i][j],network.m_structinfo[offs+0],f,df,d2f);
            //--- change values
            network.m_chunks[i].Set(j,f);
            network.m_chunks[idfdnet+i].Set(j,df);
           }
         continue;
        }
      //--- check
      if(network.m_structinfo[offs+0]==0)
        {
         //--- Adaptive summator:
         //--- * calculate NET vector,NET(i)=SUM(W(j,i)*Neurons(j),j=N1..N2)
         n1=network.m_structinfo[offs+2];
         n2=n1+network.m_structinfo[offs+1]-1;
         w1=network.m_structinfo[offs+3];
         w2=w1+network.m_structinfo[offs+1]-1;
         //--- calculation
         for(i_=0;i_<=csize-1;i_++)
            network.m_chunks[i].Set(i_,network.m_chunks[izeros][i_]);
         for(j=n1;j<=n2;j++)
           {
            v=network.m_weights[w1+j-n1];
            for(i_=0;i_<=csize-1;i_++)
               network.m_chunks[i].Set(i_,network.m_chunks[i][i_]+v*network.m_chunks[j][i_]);
           }
         continue;
        }
      //--- check
      if(network.m_structinfo[offs+0]<0)
        {
         bflag=false;
         //--- check
         if(network.m_structinfo[offs+0]==-2)
           {
            //--- input neuron,left unchanged
            bflag=true;
           }
         //--- check
         if(network.m_structinfo[offs+0]==-3)
           {
            //--- "-1" neuron
            for(k=0;k<=csize-1;k++)
               network.m_chunks[i].Set(k,-1);
            bflag=true;
           }
         //--- check
         if(network.m_structinfo[offs+0]==-4)
           {
            //--- "0" neuron
            for(k=0;k<=csize-1;k++)
               network.m_chunks[i].Set(k,0);
            bflag=true;
           }
         //--- check
         if(!CAp::Assert(bflag,__FUNCTION__+": internal error - unknown neuron type!"))
            return;
         continue;
        }
     }
//--- Post-processing,error,dError/dOut
   for(i=0;i<=ntotal-1;i++)
     {
      for(i_=0;i_<=csize-1;i_++)
         network.m_chunks[iderror+i].Set(i_,network.m_chunks[izeros][i_]);
     }
//--- check
   if(!CAp::Assert(network.m_structinfo[6]==0 || network.m_structinfo[6]==1,__FUNCTION__+": unknown normalization type!"))
      return;
//--- check
   if(network.m_structinfo[6]==1)
     {
      //--- Softmax output,classification network.
      //--- For each K=0..CSize-1 do:
      //--- 1. place exp(outputs[k]) to NWBuf[0:NOut-1]
      //--- 2. place sum(exp(..)) to NET
      //--- 3. calculate dError/dOut and place it to the second block of Chunks
      for(k=0;k<=csize-1;k++)
        {
         //--- Normalize
         mx=network.m_chunks[ntotal-nout][k];
         for(i=1;i<=nout-1;i++)
            mx=MathMax(mx,network.m_chunks[ntotal-nout+i][k]);
         net=0;
         for(i=0;i<=nout-1;i++)
           {
            network.m_nwbuf[i]=MathExp(network.m_chunks[ntotal-nout+i][k]-mx);
            net=net+network.m_nwbuf[i];
           }
         //--- Calculate error function and dError/dOut
         if(naturalerrorfunc)
           {
            //--- Natural error func.
            s=1;
            lnnet=MathLog(net);
            kl=(int)MathRound(xy[cstart+k][nin]);
            //--- calculation
            for(i=0;i<=nout-1;i++)
              {
               //--- check
               if(i==kl)
                  v=1;
               else
                  v=0;
               network.m_chunks[iderror+ntotal-nout+i].Set(k,s*network.m_nwbuf[i]/net-v);
               e=e+SafeCrossEntropy(v,network.m_nwbuf[i]/net);
              }
           }
         else
           {
            //--- Least squares error func
            //--- Error,dError/dOut(normalized)
            kl=(int)MathRound(xy[cstart+k][nin]);
            for(i=0;i<=nout-1;i++)
              {
               //--- check
               if(i==kl)
                  v=network.m_nwbuf[i]/net-1;
               else
                  v=network.m_nwbuf[i]/net;
               network.m_nwbuf[nout+i]=v;
               e=e+CMath::Sqr(v)/2;
              }
            //--- From dError/dOut(normalized) to dError/dOut(non-normalized)
            i1_=-nout;
            v=0.0;
            for(i_=nout;i_<=2*nout-1;i_++)
               v+=network.m_nwbuf[i_]*network.m_nwbuf[i_+i1_];
            //--- calculation
            for(i=0;i<=nout-1;i++)
              {
               fown=network.m_nwbuf[i];
               deown=network.m_nwbuf[nout+i];
               network.m_chunks[iderror+ntotal-nout+i].Set(k,(-v+deown*fown+deown*(net-fown))*fown/CMath::Sqr(net));
              }
           }
        }
     }
   else
     {
      //--- Normal output,regression network
      //--- For each K=0..CSize-1 do:
      //--- 1. calculate dError/dOut and place it to the second block of Chunks
      for(i=0;i<=nout-1;i++)
        {
         for(j=0;j<=csize-1;j++)
           {
            v=network.m_chunks[ntotal-nout+i][j]*network.m_columnsigmas[nin+i]+network.m_columnmeans[nin+i]-xy[cstart+j][nin+i];
            network.m_chunks[iderror+ntotal-nout+i].Set(j,v*network.m_columnsigmas[nin+i]);
            e=e+CMath::Sqr(v)/2;
           }
        }
     }
//--- Backpropagation
   for(i=ntotal-1;i>=0;i--)
     {
      //--- Extract info
      offs=istart+i*m_nfieldwidth;
      //--- check
      if(network.m_structinfo[offs+0]>0 || network.m_structinfo[offs+0]==-5)
        {
         //--- Activation function
         n1=network.m_structinfo[offs+2];
         for(k=0;k<=csize-1;k++)
            network.m_chunks[iderror+i].Set(k,network.m_chunks[iderror+i][k]*network.m_chunks[idfdnet+i][k]);
         for(i_=0;i_<=csize-1;i_++)
            network.m_chunks[iderror+n1].Set(i_,network.m_chunks[iderror+n1][i_]+network.m_chunks[iderror+i][i_]);
         continue;
        }
      //--- check
      if(network.m_structinfo[offs+0]==0)
        {
         //--- "Normal" activation function
         n1=network.m_structinfo[offs+2];
         n2=n1+network.m_structinfo[offs+1]-1;
         w1=network.m_structinfo[offs+3];
         w2=w1+network.m_structinfo[offs+1]-1;
         //--- calculation
         for(j=w1;j<=w2;j++)
           {
            v=0.0;
            for(i_=0;i_<=csize-1;i_++)
               v+=network.m_chunks[n1+j-w1][i_]*network.m_chunks[iderror+i][i_];
            grad[j]=grad[j]+v;
           }
         //--- calculation
         for(j=n1;j<=n2;j++)
           {
            v=network.m_weights[w1+j-n1];
            for(i_=0;i_<=csize-1;i_++)
               network.m_chunks[iderror+j].Set(i_,network.m_chunks[iderror+j][i_]+v*network.m_chunks[iderror+i][i_]);
           }
         continue;
        }
      //--- check
      if(network.m_structinfo[offs+0]<0)
        {
         bflag=false;
         //--- check
         if((network.m_structinfo[offs+0]==-2 || network.m_structinfo[offs+0]==-3) || network.m_structinfo[offs+0]==-4)
           {
            //--- Special neuron type,no back-propagation required
            bflag=true;
           }
         //--- check
         if(!CAp::Assert(bflag,__FUNCTION__+": unknown neuron type!"))
            return;
         continue;
        }
     }
  }
//+------------------------------------------------------------------+
//| Returns T*Ln(T/Z), guarded against overflow/underflow.           |
//| Internal subroutine.                                             |
//+------------------------------------------------------------------+
static double CMLPBase::SafeCrossEntropy(const double t,const double z)
  {
//--- create variables
   double result=0;
   double r=0;
//--- check
   if(t==0.0)
      result=0;
   else
     {
      //--- check
      if(MathAbs(z)>1.0)
        {
         //--- Shouldn't be the case with softmax,
         //--- but we just want to be sure.
         if(t/z==0.0)
            r=CMath::m_minrealnumber;
         else
            r=t/z;
        }
      else
        {
         //--- Normal case
         if(z==0.0 || MathAbs(t)>=CMath::m_maxrealnumber*MathAbs(z))
            r=CMath::m_maxrealnumber;
         else
            r=t/z;
        }
      //--- get result
      result=t*MathLog(r);
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CLogit                                       |
//+------------------------------------------------------------------+
class CLogitModel
  {
public:
   double            m_w[];
   //--- constructor, destructor
                     CLogitModel(void);
                    ~CLogitModel(void);
   //--- copy
   void              Copy(CLogitModel &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLogitModel::CLogitModel(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLogitModel::~CLogitModel(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CLogitModel::Copy(CLogitModel &obj)
  {
//--- copy array
   ArrayCopy(m_w,obj.m_w);
  }
//+------------------------------------------------------------------+
//| This class is a shell for class CLogitModel                      |
//+------------------------------------------------------------------+
class CLogitModelShell
  {
private:
   CLogitModel       m_innerobj;
public:
   //--- constructors, destructor
                     CLogitModelShell(void);
                     CLogitModelShell(CLogitModel &obj);
                    ~CLogitModelShell(void);
   //--- method
   CLogitModel      *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLogitModelShell::CLogitModelShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CLogitModelShell::CLogitModelShell(CLogitModel &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLogitModelShell::~CLogitModelShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CLogitModel *CLogitModelShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Auxiliary class for CLogit                                       |
//+------------------------------------------------------------------+
class CLogitMCState
  {
public:
   //--- variables
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
                     CLogitMCState(void);
                    ~CLogitMCState(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLogitMCState::CLogitMCState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLogitMCState::~CLogitMCState(void)
  {

  }
//+------------------------------------------------------------------+
//| MNLReport structure contains information about training process: |
//| * NGrad     -   number of gradient calculations                  |
//| * NHess     -   number of Hessian calculations                   |
//+------------------------------------------------------------------+
class CMNLReport
  {
public:
   //--- variables
   int               m_ngrad;
   int               m_nhess;
   //--- constructor, destructor
                     CMNLReport(void);
                    ~CMNLReport(void);
   //--- copy
   void              Copy(CMNLReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMNLReport::CMNLReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMNLReport::~CMNLReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMNLReport::Copy(CMNLReport &obj)
  {
//--- copy variables
   m_ngrad=obj.m_ngrad;
   m_nhess=obj.m_nhess;
  }
//+------------------------------------------------------------------+
//| MNLReport structure contains information about training process: |
//| * NGrad     -   number of gradient calculations                  |
//| * NHess     -   number of Hessian calculations                   |
//+------------------------------------------------------------------+
class CMNLReportShell
  {
private:
   CMNLReport        m_innerobj;
public:
   //--- constructors, destructor
                     CMNLReportShell(void);
                     CMNLReportShell(CMNLReport &obj);
                    ~CMNLReportShell(void);
   //--- methods
   int               GetNGrad(void);
   void              SetNGrad(const int i);
   int               GetNHess(void);
   void              SetNHess(const int i);
   CMNLReport       *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMNLReportShell::CMNLReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMNLReportShell::CMNLReportShell(CMNLReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMNLReportShell::~CMNLReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable ngrad                          |
//+------------------------------------------------------------------+
int CMNLReportShell::GetNGrad(void)
  {
//--- return result
   return(m_innerobj.m_ngrad);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable ngrad                         |
//+------------------------------------------------------------------+
void CMNLReportShell::SetNGrad(const int i)
  {
//--- change value
   m_innerobj.m_ngrad=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable nhess                          |
//+------------------------------------------------------------------+
int CMNLReportShell::GetNHess(void)
  {
//--- return result
   return(m_innerobj.m_nhess);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable nhess                         |
//+------------------------------------------------------------------+
void CMNLReportShell::SetNHess(const int i)
  {
//--- change value
   m_innerobj.m_nhess=i;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMNLReport *CMNLReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Class logit model                                                |
//+------------------------------------------------------------------+
class CLogit
  {
private:
   //--- private methods
   static void       MNLIExp(double &w[],double &x[]);
   static void       MNLAllErrors(CLogitModel &lm,CMatrixDouble &xy,const int npoints,double &relcls,double &avgce,double &rms,double &avg,double &avgrel);
   static void       MNLMCSrch(const int n,double &x[],double &f,double &g[],double &s[],double &stp,int &info,int &nfev,double &wa[],CLogitMCState &state,int &stage);
   static void       MNLMCStep(double &stx,double &fx,double &dx,double &sty,double &fy,double &dy,double &stp,const double fp,const double dp,bool &brackt,const double stmin,const double stmax,int &info);
public:
   //--- variables
   static const double m_xtol;
   static const double m_ftol;
   static const double m_gtol;
   static const int  m_maxfev;
   static const double m_stpmin;
   static const double m_stpmax;
   static const int  m_logitvnum;
   //--- constructor, destructor
                     CLogit(void);
                    ~CLogit(void);
   //--- public methods
   static void       MNLTrainH(CMatrixDouble &xy,const int npoints,const int nvars,const int nclasses,int &info,CLogitModel &lm,CMNLReport &rep);
   static void       MNLProcess(CLogitModel &lm,double &x[],double &y[]);
   static void       MNLProcessI(CLogitModel &lm,double &x[],double &y[]);
   static void       MNLUnpack(CLogitModel &lm,CMatrixDouble &a,int &nvars,int &nclasses);
   static void       MNLPack(CMatrixDouble &a,const int nvars,const int nclasses,CLogitModel &lm);
   static void       MNLCopy(CLogitModel &lm1,CLogitModel &lm2);
   static double     MNLAvgCE(CLogitModel &lm,CMatrixDouble &xy,const int npoints);
   static double     MNLRelClsError(CLogitModel &lm,CMatrixDouble &xy,const int npoints);
   static double     MNLRMSError(CLogitModel &lm,CMatrixDouble &xy,const int npoints);
   static double     MNLAvgError(CLogitModel &lm,CMatrixDouble &xy,const int npoints);
   static double     MNLAvgRelError(CLogitModel &lm,CMatrixDouble &xy,const int ssize);
   static int        MNLClsError(CLogitModel &lm,CMatrixDouble &xy,const int npoints);
  };
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const double CLogit::m_xtol=100*CMath::m_machineepsilon;
const double CLogit::m_ftol=0.0001;
const double CLogit::m_gtol=0.3;
const int    CLogit::m_maxfev=20;
const double CLogit::m_stpmin=1.0E-2;
const double CLogit::m_stpmax=1.0E5;
const int    CLogit::m_logitvnum=6;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLogit::CLogit(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLogit::~CLogit(void)
  {

  }
//+------------------------------------------------------------------+
//| This subroutine trains logit model.                              |
//| INPUT PARAMETERS:                                                |
//|     XY          -   training set, array[0..NPoints-1,0..NVars]   |
//|                     First NVars columns store values of          |
//|                     independent variables, next column stores    |
//|                     number of class (from 0 to NClasses-1) which |
//|                     dataset element belongs to. Fractional values|
//|                     are rounded to nearest integer.              |
//|     NPoints     -   training set size, NPoints>=1                |
//|     NVars       -   number of independent variables, NVars>=1    |
//|     NClasses    -   number of classes, NClasses>=2               |
//| OUTPUT PARAMETERS:                                               |
//|     Info        -   return code:                                 |
//|                     * -2, if there is a point with class number  |
//|                           outside of [0..NClasses-1].            |
//|                     * -1, if incorrect parameters was passed     |
//|                           (NPoints<NVars+2, NVars<1, NClasses<2).|
//|                     *  1, if task has been solved                |
//|     LM          -   model built                                  |
//|     Rep         -   training report                              |
//+------------------------------------------------------------------+
static void CLogit::MNLTrainH(CMatrixDouble &xy,const int npoints,
                              const int nvars,const int nclasses,
                              int &info,CLogitModel &lm,CMNLReport &rep)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    k=0;
   int    ssize=0;
   bool   allsame;
   int    offs=0;
   double threshold=0;
   double wminstep=0;
   double decay=0;
   int    wdim=0;
   int    expoffs=0;
   double v=0;
   double s=0;
   int    nin=0;
   int    nout=0;
   int    wcount=0;
   double e=0;
   bool   spd;
   double wstep=0;
   int    mcstage=0;
   int    mcinfo=0;
   int    mcnfev=0;
   int    solverinfo=0;
   int    i_=0;
   int    i1_=0;
//--- creating arrays
   double g[];
   double x[];
   double y[];
   double wbase[];
   double wdir[];
   double work[];
//--- create matrix
   CMatrixDouble h;
//--- create objects of classes
   CLogitMCState         mcstate;
   CDenseSolverReport    solverrep;
   CMultilayerPerceptron network;
//--- initialization
   info=0;
   threshold=1000*CMath::m_machineepsilon;
   wminstep=0.001;
   decay=0.001;
//--- Test for inputs
   if((npoints<nvars+2 || nvars<1) || nclasses<2)
     {
      info=-1;
      return;
     }
   for(i=0;i<=npoints-1;i++)
     {
      //--- check
      if((int)MathRound(xy[i][nvars])<0 || (int)MathRound(xy[i][nvars])>=nclasses)
        {
         info=-2;
         return;
        }
     }
//--- change value
   info=1;
//--- Initialize data
   rep.m_ngrad=0;
   rep.m_nhess=0;
//--- Allocate array
   wdim=(nvars+1)*(nclasses-1);
   offs=5;
   expoffs=offs+wdim;
   ssize=5+(nvars+1)*(nclasses-1)+nclasses;
//--- allocation
   ArrayResizeAL(lm.m_w,ssize);
//--- change values
   lm.m_w[0]=ssize;
   lm.m_w[1]=m_logitvnum;
   lm.m_w[2]=nvars;
   lm.m_w[3]=nclasses;
   lm.m_w[4]=offs;
//--- Degenerate case: all outputs are equal
   allsame=true;
   for(i=1;i<=npoints-1;i++)
     {
      //--- check
      if((int)MathRound(xy[i][nvars])!=(int)MathRound(xy[i-1][nvars]))
         allsame=false;
     }
//--- check
   if(allsame)
     {
      for(i=0;i<=(nvars+1)*(nclasses-1)-1;i++)
         lm.m_w[offs+i]=0;
      //--- change values
      v=-(2*MathLog(CMath::m_minrealnumber));
      k=(int)MathRound(xy[0][nvars]);
      //--- check
      if(k==nclasses-1)
        {
         for(i=0;i<=nclasses-2;i++)
            lm.m_w[offs+i*(nvars+1)+nvars]=-v;
        }
      else
        {
         for(i=0;i<=nclasses-2;i++)
           {
            //--- check
            if(i==k)
               lm.m_w[offs+i*(nvars+1)+nvars]=v;
            else
               lm.m_w[offs+i*(nvars+1)+nvars]=0;
           }
        }
      //--- exit the function
      return;
     }
//--- General case.
//--- Prepare task and network. Allocate space.
   CMLPBase::MLPCreateC0(nvars,nclasses,network);
//--- function call
   CMLPBase::MLPInitPreprocessor(network,xy,npoints);
//--- function call
   CMLPBase::MLPProperties(network,nin,nout,wcount);
   for(i=0;i<=wcount-1;i++)
      network.m_weights[i]=(2*CMath::RandomReal()-1)/nvars;
//--- allocation
   ArrayResizeAL(g,wcount);
   h.Resize(wcount,wcount);
   ArrayResizeAL(wbase,wcount);
   ArrayResizeAL(wdir,wcount);
   ArrayResizeAL(work,wcount);
//--- First stage: optimize in gradient direction.
   for(k=0;k<=wcount/3+10;k++)
     {
      //--- Calculate gradient in starting point
      CMLPBase::MLPGradNBatch(network,xy,npoints,e,g);
      v=0.0;
      for(i_=0;i_<=wcount-1;i_++)
         v+=network.m_weights[i_]*network.m_weights[i_];
      //--- change value
      e=e+0.5*decay*v;
      for(i_=0;i_<=wcount-1;i_++)
         g[i_]=g[i_]+decay*network.m_weights[i_];
      rep.m_ngrad=rep.m_ngrad+1;
      //--- Setup optimization scheme
      for(i_=0;i_<=wcount-1;i_++)
         wdir[i_]=-g[i_];
      v=0.0;
      for(i_=0;i_<=wcount-1;i_++)
         v+=wdir[i_]*wdir[i_];
      //--- change values
      wstep=MathSqrt(v);
      v=1/MathSqrt(v);
      for(i_=0;i_<=wcount-1;i_++)
         wdir[i_]=v*wdir[i_];
      mcstage=0;
      //--- function call
      MNLMCSrch(wcount,network.m_weights,e,g,wdir,wstep,mcinfo,mcnfev,work,mcstate,mcstage);
      //--- cycle
      while(mcstage!=0)
        {
         //--- function call
         CMLPBase::MLPGradNBatch(network,xy,npoints,e,g);
         v=0.0;
         for(i_=0;i_<=wcount-1;i_++)
            v+=network.m_weights[i_]*network.m_weights[i_];
         //--- change value
         e=e+0.5*decay*v;
         for(i_=0;i_<=wcount-1;i_++)
            g[i_]=g[i_]+decay*network.m_weights[i_];
         rep.m_ngrad=rep.m_ngrad+1;
         //--- function call
         MNLMCSrch(wcount,network.m_weights,e,g,wdir,wstep,mcinfo,mcnfev,work,mcstate,mcstage);
        }
     }
//--- Second stage: use Hessian when we are close to the minimum
   while(true)
     {
      //--- Calculate and update E/G/H
      CMLPBase::MLPHessianNBatch(network,xy,npoints,e,g,h);
      v=0.0;
      for(i_=0;i_<=wcount-1;i_++)
         v+=network.m_weights[i_]*network.m_weights[i_];
      //--- change value
      e=e+0.5*decay*v;
      for(i_=0;i_<=wcount-1;i_++)
         g[i_]=g[i_]+decay*network.m_weights[i_];
      for(k=0;k<=wcount-1;k++)
         h[k].Set(k,h[k][k]+decay);
      rep.m_nhess=rep.m_nhess+1;
      //--- Select step direction
      //--- NOTE: it is important to use lower-triangle Cholesky
      //--- factorization since it is much faster than higher-triangle version.
      spd=CTrFac::SPDMatrixCholesky(h,wcount,false);
      //--- function call
      CDenseSolver::SPDMatrixCholeskySolve(h,wcount,false,g,solverinfo,solverrep,wdir);
      spd=solverinfo>0;
      //--- check
      if(spd)
        {
         //--- H is positive definite.
         //--- Step in Newton direction.
         for(i_=0;i_<=wcount-1;i_++)
            wdir[i_]=-1*wdir[i_];
         spd=true;
        }
      else
        {
         //--- H is indefinite.
         //--- Step in gradient direction.
         for(i_=0;i_<=wcount-1;i_++)
            wdir[i_]=-g[i_];
         spd=false;
        }
      //--- Optimize in WDir direction
      v=0.0;
      for(i_=0;i_<=wcount-1;i_++)
         v+=wdir[i_]*wdir[i_];
      //--- change values
      wstep=MathSqrt(v);
      v=1/MathSqrt(v);
      for(i_=0;i_<=wcount-1;i_++)
         wdir[i_]=v*wdir[i_];
      mcstage=0;
      //--- function call
      MNLMCSrch(wcount,network.m_weights,e,g,wdir,wstep,mcinfo,mcnfev,work,mcstate,mcstage);
      //--- cycle
      while(mcstage!=0)
        {
         //--- function call
         CMLPBase::MLPGradNBatch(network,xy,npoints,e,g);
         v=0.0;
         for(i_=0;i_<=wcount-1;i_++)
            v+=network.m_weights[i_]*network.m_weights[i_];
         //--- change value
         e=e+0.5*decay*v;
         for(i_=0;i_<=wcount-1;i_++)
            g[i_]=g[i_]+decay*network.m_weights[i_];
         rep.m_ngrad=rep.m_ngrad+1;
         //--- function call
         MNLMCSrch(wcount,network.m_weights,e,g,wdir,wstep,mcinfo,mcnfev,work,mcstate,mcstage);
        }
      //--- check
      if(spd && ((mcinfo==2 || mcinfo==4) || mcinfo==6))
         break;
     }
//--- Convert from NN format to MNL format
   i1_=-offs;
   for(i_=offs;i_<=offs+wcount-1;i_++)
      lm.m_w[i_]=network.m_weights[i_+i1_];
   for(k=0;k<=nvars-1;k++)
     {
      for(i=0;i<=nclasses-2;i++)
        {
         s=network.m_columnsigmas[k];
         //--- check
         if(s==0.0)
            s=1;
         //--- change values
         j=offs+(nvars+1)*i;
         v=lm.m_w[j+k];
         lm.m_w[j+k]=v/s;
         lm.m_w[j+nvars]=lm.m_w[j+nvars]+v*network.m_columnmeans[k]/s;
        }
     }
//--- calculation
   for(k=0;k<=nclasses-2;k++)
      lm.m_w[offs+(nvars+1)*k+nvars]=-lm.m_w[offs+(nvars+1)*k+nvars];
  }
//+------------------------------------------------------------------+
//| Procesing                                                        |
//| INPUT PARAMETERS:                                                |
//|     LM      -   logit model, passed by non-constant reference    |
//|                 (some fields of structure are used as temporaries|
//|                 when calculating model output).                  |
//|     X       -   input vector,  array[0..NVars-1].                |
//|     Y       -   (possibly) preallocated buffer; if size of Y is  |
//|                 less than NClasses, it will be reallocated.If it |
//|                 is large enough, it is NOT reallocated, so we    |
//|                 can save some time on reallocation.              |
//| OUTPUT PARAMETERS:                                               |
//|     Y       -   result, array[0..NClasses-1]                     |
//|                 Vector of posterior probabilities for            |
//|                 classification task.                             |
//+------------------------------------------------------------------+
static void CLogit::MNLProcess(CLogitModel &lm,double &x[],double &y[])
  {
//--- create variables
   int    nvars=0;
   int    nclasses=0;
   int    offs=0;
   int    i=0;
   int    i1=0;
   double s=0;
//--- check
   if(!CAp::Assert(lm.m_w[1]==m_logitvnum,__FUNCTION__+": unexpected model version"))
      return;
//--- initialization
   nvars=(int)MathRound(lm.m_w[2]);
   nclasses=(int)MathRound(lm.m_w[3]);
   offs=(int)MathRound(lm.m_w[4]);
//--- function call
   MNLIExp(lm.m_w,x);
   s=0;
//--- calculation
   i1=offs+(nvars+1)*(nclasses-1);
   for(i=i1;i<=i1+nclasses-1;i++)
      s=s+lm.m_w[i];
//--- check
   if(CAp::Len(y)<nclasses)
      ArrayResizeAL(y,nclasses);
//--- change values
   for(i=0;i<=nclasses-1;i++)
      y[i]=lm.m_w[i1+i]/s;
  }
//+------------------------------------------------------------------+
//| 'interactive' variant of MNLProcess for languages like Python    |
//| which support constructs like "Y=MNLProcess(LM,X)" and           |
//| interactive mode of the interpreter                              |
//| This function allocates new array on each call, so it is         |
//| significantly slower than its 'non-interactive' counterpart,     |
//| but it is  more  convenient when you call it from command line.  |
//+------------------------------------------------------------------+
static void CLogit::MNLProcessI(CLogitModel &lm,double &x[],double &y[])
  {
//--- function call
   MNLProcess(lm,x,y);
  }
//+------------------------------------------------------------------+
//| Unpacks coefficients of logit model. Logit model have form:      |
//|     P(class=i) = S(i) / (S(0) + S(1) + ... +S(M-1))              |
//|     S(i) = Exp(A[i,0]*X[0] + ... + A[i,N-1]*X[N-1] + A[i,N]),    |
//|            when i<M-1                                            |
//|     S(M-1) = 1                                                   |
//| INPUT PARAMETERS:                                                |
//|     LM          -   logit model in ALGLIB format                 |
//| OUTPUT PARAMETERS:                                               |
//|     V           -   coefficients, array[0..NClasses-2,0..NVars]  |
//|     NVars       -   number of independent variables              |
//|     NClasses    -   number of classes                            |
//+------------------------------------------------------------------+
static void CLogit::MNLUnpack(CLogitModel &lm,CMatrixDouble &a,int &nvars,
                              int &nclasses)
  {
//--- create variables
   int offs=0;
   int i=0;
   int i_=0;
   int i1_=0;
//--- initialization
   nvars=0;
   nclasses=0;
//--- check
   if(!CAp::Assert(lm.m_w[1]==m_logitvnum,__FUNCTION__+": unexpected model version"))
      return;
//--- initialization
   nvars=(int)MathRound(lm.m_w[2]);
   nclasses=(int)MathRound(lm.m_w[3]);
   offs=(int)MathRound(lm.m_w[4]);
//--- allocation
   a.Resize(nclasses-1,nvars+1);
//--- calculation
   for(i=0;i<=nclasses-2;i++)
     {
      i1_=offs+i*(nvars+1);
      for(i_=0;i_<=nvars;i_++)
         a[i].Set(i_,lm.m_w[i_+i1_]);
     }
  }
//+------------------------------------------------------------------+
//| "Packs" coefficients and creates logit model in ALGLIB format    |
//| (MNLUnpack reversed).                                            |
//| INPUT PARAMETERS:                                                |
//|     A           -   model (see MNLUnpack)                        |
//|     NVars       -   number of independent variables              |
//|     NClasses    -   number of classes                            |
//| OUTPUT PARAMETERS:                                               |
//|     LM          -   logit model.                                 |
//+------------------------------------------------------------------+
static void CLogit::MNLPack(CMatrixDouble &a,const int nvars,const int nclasses,
                            CLogitModel &lm)
  {
//--- create variables
   int offs=0;
   int i=0;
   int wdim=0;
   int ssize=0;
   int i_=0;
   int i1_=0;
//--- initialization
   wdim=(nvars+1)*(nclasses-1);
   offs=5;
   ssize=5+(nvars+1)*(nclasses-1)+nclasses;
//--- allocation
   ArrayResizeAL(lm.m_w,ssize);
//--- initialization
   lm.m_w[0]=ssize;
   lm.m_w[1]=m_logitvnum;
   lm.m_w[2]=nvars;
   lm.m_w[3]=nclasses;
   lm.m_w[4]=offs;
//--- calculation
   for(i=0;i<=nclasses-2;i++)
     {
      i1_=-(offs+i*(nvars+1));
      for(i_=offs+i*(nvars+1);i_<=offs+i*(nvars+1)+nvars;i_++)
         lm.m_w[i_]=a[i][i_+i1_];
     }
  }
//+------------------------------------------------------------------+
//| Copying of LogitModel strucure                                   |
//| INPUT PARAMETERS:                                                |
//|     LM1 -   original                                             |
//| OUTPUT PARAMETERS:                                               |
//|     LM2 -   copy                                                 |
//+------------------------------------------------------------------+
static void CLogit::MNLCopy(CLogitModel &lm1,CLogitModel &lm2)
  {
//--- create variables
   int k=0;
   int i_=0;
//--- initialization
   k=(int)MathRound(lm1.m_w[0]);
//--- allocation
   ArrayResizeAL(lm2.m_w,k);
//--- copy
   for(i_=0;i_<=k-1;i_++)
      lm2.m_w[i_]=lm1.m_w[i_];
  }
//+------------------------------------------------------------------+
//| Average cross-entropy (in bits per element) on the test set      |
//| INPUT PARAMETERS:                                                |
//|     LM      -   logit model                                      |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     CrossEntropy/(NPoints*ln(2)).                                |
//+------------------------------------------------------------------+
static double CLogit::MNLAvgCE(CLogitModel &lm,CMatrixDouble &xy,const int npoints)
  {
//--- create variables
   double result=0;
   int    nvars=0;
   int    nclasses=0;
   int    i=0;
   int    i_=0;
//--- creating arrays
   double workx[];
   double worky[];
//--- check
   if(!CAp::Assert(lm.m_w[1]==m_logitvnum,__FUNCTION__+": unexpected model version"))
      return(EMPTY_VALUE);
//--- initialization
   nvars=(int)MathRound(lm.m_w[2]);
   nclasses=(int)MathRound(lm.m_w[3]);
//--- allocation
   ArrayResizeAL(workx,nvars);
   ArrayResizeAL(worky,nclasses);
//--- calculation
   for(i=0;i<=npoints-1;i++)
     {
      //--- check
      if(!CAp::Assert((int)MathRound(xy[i][nvars])>=0 &&(int)MathRound(xy[i][nvars])<nclasses,__FUNCTION__+": incorrect class number!"))
         return(EMPTY_VALUE);
      //--- Process
      for(i_=0;i_<=nvars-1;i_++)
         workx[i_]=xy[i][i_];
      //--- function call
      MNLProcess(lm,workx,worky);
      //--- check
      if(worky[(int)MathRound(xy[i][nvars])]>0.0)
         result=result-MathLog(worky[(int)MathRound(xy[i][nvars])]);
      else
         result=result-MathLog(CMath::m_minrealnumber);
     }
//--- return result
   return(result/(npoints*MathLog(2)));
  }
//+------------------------------------------------------------------+
//| Relative classification error on the test set                    |
//| INPUT PARAMETERS:                                                |
//|     LM      -   logit model                                      |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     percent of incorrectly classified cases.                     |
//+------------------------------------------------------------------+
static double CLogit::MNLRelClsError(CLogitModel &lm,CMatrixDouble &xy,
                                     const int npoints)
  {
//--- return result
   return((double)MNLClsError(lm,xy,npoints)/(double)npoints);
  }
//+------------------------------------------------------------------+
//| RMS error on the test set                                        |
//| INPUT PARAMETERS:                                                |
//|     LM      -   logit model                                      |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     root mean square error (error when estimating posterior      |
//|     probabilities).                                              |
//+------------------------------------------------------------------+
static double CLogit::MNLRMSError(CLogitModel &lm,CMatrixDouble &xy,
                                  const int npoints)
  {
//--- create variables
   double relcls=0;
   double avgce=0;
   double rms=0;
   double avg=0;
   double avgrel=0;
//--- check
   if(!CAp::Assert((int)MathRound(lm.m_w[1])==m_logitvnum,__FUNCTION__+": Incorrect MNL version!"))
      return(EMPTY_VALUE);
//--- function call
   MNLAllErrors(lm,xy,npoints,relcls,avgce,rms,avg,avgrel);
//--- return result
   return(rms);
  }
//+------------------------------------------------------------------+
//| Average error on the test set                                    |
//| INPUT PARAMETERS:                                                |
//|     LM      -   logit model                                      |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     average error (error when estimating posterior               |
//|     probabilities).                                              |
//+------------------------------------------------------------------+
static double CLogit::MNLAvgError(CLogitModel &lm,CMatrixDouble &xy,
                                  const int npoints)
  {
//--- create variables
   double relcls=0;
   double avgce=0;
   double rms=0;
   double avg=0;
   double avgrel=0;
//--- check
   if(!CAp::Assert((int)MathRound(lm.m_w[1])==m_logitvnum,__FUNCTION__+": Incorrect MNL version!"))
      return(EMPTY_VALUE);
//--- function call
   MNLAllErrors(lm,xy,npoints,relcls,avgce,rms,avg,avgrel);
//--- return result
   return(avg);
  }
//+------------------------------------------------------------------+
//| Average relative error on the test set                           |
//| INPUT PARAMETERS:                                                |
//|     LM      -   logit model                                      |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     average relative error (error when estimating posterior      |
//|     probabilities).                                              |
//+------------------------------------------------------------------+
static double CLogit::MNLAvgRelError(CLogitModel &lm,CMatrixDouble &xy,
                                     const int ssize)
  {
//--- create variables
   double relcls=0;
   double avgce=0;
   double rms=0;
   double avg=0;
   double avgrel=0;
//--- check
   if(!CAp::Assert((int)MathRound(lm.m_w[1])==m_logitvnum,__FUNCTION__+": Incorrect MNL version!"))
      return(EMPTY_VALUE);
//--- function call
   MNLAllErrors(lm,xy,ssize,relcls,avgce,rms,avg,avgrel);
//--- return result
   return(avgrel);
  }
//+------------------------------------------------------------------+
//| Classification error on test set = MNLRelClsError*NPoints        |
//+------------------------------------------------------------------+
static int CLogit::MNLClsError(CLogitModel &lm,CMatrixDouble &xy,const int npoints)
  {
//--- create variables
   int result=0;
   int nvars=0;
   int nclasses=0;
   int i=0;
   int j=0;
   int nmax=0;
   int i_=0;
//--- creating arrays
   double workx[];
   double worky[];
//--- check
   if(!CAp::Assert(lm.m_w[1]==m_logitvnum,__FUNCTION__+": unexpected model version"))
      return(-1);
//--- initialization
   nvars=(int)MathRound(lm.m_w[2]);
   nclasses=(int)MathRound(lm.m_w[3]);
//--- allocation
   ArrayResizeAL(workx,nvars);
   ArrayResizeAL(worky,nclasses);
//--- calculation
   for(i=0;i<=npoints-1;i++)
     {
      //--- Process
      for(i_=0;i_<=nvars-1;i_++)
         workx[i_]=xy[i][i_];
      //--- function call
      MNLProcess(lm,workx,worky);
      //--- Logit version of the answer
      nmax=0;
      for(j=0;j<=nclasses-1;j++)
        {
         //--- check
         if(worky[j]>worky[nmax])
            nmax=j;
        }
      //--- compare
      if(nmax!=(int)MathRound(xy[i][nvars]))
         result=result+1;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Internal subroutine. Places exponents of the anti-overflow       |
//| shifted internal linear outputs into the service part of the W   |
//| array.                                                           |
//+------------------------------------------------------------------+
static void CLogit::MNLIExp(double &w[],double &x[])
  {
//--- create variables
   int    nvars=0;
   int    nclasses=0;
   int    offs=0;
   int    i=0;
   int    i1=0;
   double v=0;
   double mx=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(!CAp::Assert(w[1]==m_logitvnum,__FUNCTION__+": unexpected model version"))
      return;
//--- initialization
   nvars=(int)MathRound(w[2]);
   nclasses=(int)MathRound(w[3]);
   offs=(int)MathRound(w[4]);
//--- calculation
   i1=offs+(nvars+1)*(nclasses-1);
   for(i=0;i<=nclasses-2;i++)
     {
      //--- change values
      i1_=-(offs+i*(nvars+1));
      v=0.0;
      for(i_=offs+i*(nvars+1);i_<=offs+i*(nvars+1)+nvars-1;i_++)
         v+=w[i_]*x[i_+i1_];
      w[i1+i]=v+w[offs+i*(nvars+1)+nvars];
     }
//--- change values
   w[i1+nclasses-1]=0;
   mx=0;
//--- calculation
   for(i=i1;i<=i1+nclasses-1;i++)
      mx=MathMax(mx,w[i]);
   for(i=i1;i<=i1+nclasses-1;i++)
      w[i]=MathExp(w[i]-mx);
  }
//+------------------------------------------------------------------+
//| Calculation of all types of errors                               |
//+------------------------------------------------------------------+
static void CLogit::MNLAllErrors(CLogitModel &lm,CMatrixDouble &xy,
                                 const int npoints,double &relcls,
                                 double &avgce,double &rms,double &avg,
                                 double &avgrel)
  {
//--- create variables
   int nvars=0;
   int nclasses=0;
   int i=0;
   int i_=0;
//--- creating arrays
   double buf[];
   double workx[];
   double y[];
   double dy[];
//--- initialization
   relcls=0;
   avgce=0;
   rms=0;
   avg=0;
   avgrel=0;
//--- check
   if(!CAp::Assert((int)MathRound(lm.m_w[1])==m_logitvnum,__FUNCTION__+": Incorrect MNL version!"))
      return;
//--- initialization
   nvars=(int)MathRound(lm.m_w[2]);
   nclasses=(int)MathRound(lm.m_w[3]);
//--- allocation
   ArrayResizeAL(workx,nvars);
   ArrayResizeAL(y,nclasses);
   ArrayResizeAL(dy,1);
//--- function call
   CBdSS::DSErrAllocate(nclasses,buf);
   for(i=0;i<=npoints-1;i++)
     {
      for(i_=0;i_<=nvars-1;i_++)
         workx[i_]=xy[i][i_];
      //--- function call
      MNLProcess(lm,workx,y);
      //--- change value
      dy[0]=xy[i][nvars];
      //--- function call
      CBdSS::DSErrAccumulate(buf,y,dy);
     }
//--- function call
   CBdSS::DSErrFinish(buf);
//--- change values
   relcls=buf[0];
   avgce=buf[1];
   rms=buf[2];
   avg=buf[3];
   avgrel=buf[4];
  }
//+------------------------------------------------------------------+
//| The purpose of mcsrch is to find a step which satisfies a        |
//| sufficient decrease condition and a curvature condition.         |
//| At each stage the subroutine updates an interval of uncertainty  |
//| with endpoints stx and sty. The interval of uncertainty is       |
//| initially chosen so that it contains a minimizer of the modified |
//| function                                                         |
//|     f(x+stp*s) - f(x) - ftol*stp*(gradf(x)'s).                   |
//| If a step is obtained for  which the modified function has a     |
//| nonpositive function value and nonnegative derivative, then the  |
//| interval of uncertainty is chosen so that it contains a minimizer|
//| of f(x+stp*s).                                                   |
//| The algorithm is designed to find a step which satisfies the     |
//| sufficient decrease condition                                    |
//|     f(x+stp*s) .le. f(x) + ftol*stp*(gradf(x)'s),                |
//| and the curvature condition                                      |
//|     abs(gradf(x+stp*s)'s)) .le. gtol*abs(gradf(x)'s).            |
//| If ftol is less than gtol and if, for example, the function is   |
//| bounded below, then there is always a step which satisfies both  |
//| conditions. If no step can be found which satisfies both         |
//| conditions, then the algorithm usually stops when rounding       |
//| errors prevent further progress. In this case stp only satisfies |
//| the sufficient decrease condition.                               |
//| Parameters descriprion                                           |
//| N is a positive integer input variable set to the number of      |
//| variables.                                                       |
//| X is an array of length n. on input it must contain the base     |
//| point for the line search. on output it contains x+stp*s.        |
//| F is a variable. on input it must contain the value of f at x. On|
//| output it contains the value of f at x + stp*s.                  |
//| G is an array of length n. on input it must contain the gradient |
//| of f at x. On output it contains the gradient of f at x + stp*s. |
//| s is an input array of length n which specifies the search       |
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
//|               is at most xtol.                                   |
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
static void CLogit::MNLMCSrch(const int n,double &x[],double &f,double &g[],
                              double &s[],double &stp,int &info,int &nfev,
                              double &wa[],CLogitMCState &state,int &stage)
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
         //--- CHECK THE INPUT PARAMETERS FOR ERRORS.
         if(n<=0 || stp<=0.0 || m_ftol<0.0 || m_gtol<zero || m_xtol<zero || m_stpmin<zero || m_stpmax<m_stpmin || m_maxfev<=0)
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
         state.m_width=m_stpmax-m_stpmin;
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
         //--- NEXT
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
         if(stp>m_stpmax)
            stp=m_stpmax;
         //--- check
         if(stp<m_stpmin)
            stp=m_stpmin;
         //--- if an unusual termination is to occur then let
         //--- stp be the lowest point obtained so far.
         if((((state.m_brackt && (stp<=state.m_stmin || stp>=state.m_stmax)) || nfev>=m_maxfev-1) || state.m_infoc==0) || (state.m_brackt && state.m_stmax-state.m_stmin<=m_xtol*state.m_stmax))
            stp=state.m_stx;
         //--- evaluate the function and gradient at stp
         //--- and compute the directional derivative.
         for(i_=0;i_<=n-1;i_++)
            x[i_]=wa[i_];
         for(i_=0;i_<=n-1;i_++)
            x[i_]=x[i_]+stp*s[i_];
         //--- next
         stage=4;
         return;
        }
      //--- check
      if(stage==4)
        {
         info=0;
         nfev=nfev+1;
         v=0.0;
         //--- calculation
         for(i_=0;i_<=n-1;i_++)
            v+=g[i_]*s[i_];
         state.m_dg=v;
         state.m_ftest1=state.m_finit+stp*state.m_dgtest;
         //--- test for convergence.
         if((state.m_brackt && (stp<=state.m_stmin || stp>=state.m_stmax)) || state.m_infoc==0)
            info=6;
         //--- check
         if((stp==m_stpmax &&f<=state.m_ftest1) &&state.m_dg<=state.m_dgtest)
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
         if(f<=state.m_ftest1 && MathAbs(state.m_dg)<=-(m_gtol*state.m_dginit))
            info=1;
         //--- check for termination.
         if(info!=0)
           {
            stage=0;
            return;
           }
         //--- in the first stage we seek a step for which the modified
         //--- function has a nonpositive value and nonnegative derivative.
         if((state.m_stage1 && f<=state.m_ftest1) && state.m_dg>=MathMin(m_ftol,m_gtol)*state.m_dginit)
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
            MNLMCStep(state.m_stx,state.m_fxm,state.m_dgxm,state.m_sty,state.m_fym,state.m_dgym,stp,state.m_fm,state.m_dgm,state.m_brackt,state.m_stmin,state.m_stmax,state.m_infoc);
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
            MNLMCStep(state.m_stx,state.m_fx,state.m_dgx,state.m_sty,state.m_fy,state.m_dgy,stp,f,state.m_dg,state.m_brackt,state.m_stmin,state.m_stmax,state.m_infoc);
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
//| Auxiliary function for MNLMCSrch                                 |
//+------------------------------------------------------------------+
static void CLogit::MNLMCStep(double &stx,double &fx,double &dx,double &sty,
                              double &fy,double &dy,double &stp,const double fp,
                              const double dp,bool &brackt,const double stmin,
                              const double stmax,int &info)
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
   if(((brackt && (stp<=MathMin(stx,sty) || stp>=MathMax(stx,sty))) || dx*(stp-stx)>=0.0) || stmax<stmin)
      return;
//--- determine if the derivatives have opposite sign.
   sgnd=dp*(dx/MathAbs(dx));
//--- first case. a higher function value.
//--- the minimum is bracketed. if the cubic step is closer
//--- to stx than the quadratic step,the cubic step is taken,
//--- else the average of the cubic and quadratic steps is taken.
   if(fp>fx)
     {
      //--- change value
      info=1;
      bound=true;
      theta=3*(fx-fp)/(stp-stx)+dx+dp;
      s=MathMax(MathAbs(theta),MathMax(MathAbs(dx),MathAbs(dp)));
      gamma=s*MathSqrt(CMath::Sqr(theta/s)-dx/s*(dp/s));
      //--- check
      if(stp<stx)
         gamma=-gamma;
      //--- change value
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
      brackt=true;
     }
   else
     {
      //--- check
      if(sgnd<0.0)
        {
         //--- second case. a lower function value and derivatives of
         //--- opposite sign. the minimum is bracketed. if the cubic
         //--- step is closer to stx than the quadratic (secant) step,
         //--- the cubic step is taken,else the quadratic step is taken.
         info=2;
         bound=false;
         theta=3*(fx-fp)/(stp-stx)+dx+dp;
         s=MathMax(MathAbs(theta),MathMax(MathAbs(dx),MathAbs(dp)));
         gamma=s*MathSqrt(CMath::Sqr(theta/s)-dx/s*(dp/s));
         //--- check
         if(stp>stx)
            gamma=-gamma;
         //--- change values
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
         brackt=true;
        }
      else
        {
         //--- check
         if(MathAbs(dp)<MathAbs(dx))
           {
            //--- third case. a lower function value,derivatives of the
            //--- same sign,and the magnitude of the derivative decreases.
            //--- the cubic step is only used if the cubic tends to infinity
            //--- in the direction of the step or if the minimum of the cubic
            //--- is beyond stp. otherwise the cubic step is defined to be
            //--- either stpmin or stpmax. the quadratic (secant) step is also
            //--- computed and if the minimum is bracketed then the the step
            //--- closest to stx is taken,else the step farthest away is taken.
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
            if(r<0.0 && gamma!=0.0)
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
            if(brackt)
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
            //--- same sign,and the magnitude of the derivative does
            //--- not decrease. if the minimum is not bracketed,the step
            //--- is either stpmin or stpmax,else the cubic step is taken.
            info=4;
            bound=false;
            //--- check
            if(brackt)
              {
               //--- change values
               theta=3*(fp-fy)/(sty-stp)+dy+dp;
               s=MathMax(MathAbs(theta),MathMax(MathAbs(dy),MathAbs(dp)));
               gamma=s*MathSqrt(CMath::Sqr(theta/s)-dy/s*(dp/s));
               //--- check
               if(stp>sty)
                  gamma=-gamma;
               //--- change values
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
      sty=stp;
      fy=fp;
      dy=dp;
     }
   else
     {
      //--- check
      if(sgnd<0.0)
        {
         sty=stx;
         fy=fx;
         dy=dx;
        }
      //--- change values
      stx=stp;
      fx=fp;
      dx=dp;
     }
//--- compute the new step and safeguard it.
   stpf=MathMin(stmax,stpf);
   stpf=MathMax(stmin,stpf);
   stp=stpf;
//--- check
   if(brackt && bound)
     {
      //--- check
      if(sty>stx)
         stp=MathMin(stx+0.66*(sty-stx),stp);
      else
         stp=MathMax(stx+0.66*(sty-stx),stp);
     }
  }
//+------------------------------------------------------------------+
//| This structure is a MCPD (Markov Chains for Population Data)     |
//| solver. You should use ALGLIB functions in order to work with    |
//| this object.                                                     |
//+------------------------------------------------------------------+
class CMCPDState
  {
public:
   //--- variables
   int               m_n;
   int               m_npairs;
   int               m_ccnt;
   double            m_regterm;
   int               m_repinneriterationscount;
   int               m_repouteriterationscount;
   int               m_repnfev;
   int               m_repterminationtype;
   CMinBLEICState    m_bs;
   CMinBLEICReport   m_br;
   //--- arrays
   int               m_states[];
   int               m_ct[];
   double            m_pw[];
   double            m_tmpp[];
   double            m_effectivew[];
   double            m_effectivebndl[];
   double            m_effectivebndu[];
   int               m_effectivect[];
   double            m_h[];
   //--- matrices
   CMatrixDouble     m_data;
   CMatrixDouble     m_ec;
   CMatrixDouble     m_bndl;
   CMatrixDouble     m_bndu;
   CMatrixDouble     m_c;
   CMatrixDouble     m_priorp;
   CMatrixDouble     m_effectivec;
   CMatrixDouble     m_p;
   //--- constructor, destructor
                     CMCPDState(void);
                    ~CMCPDState(void);
   //--- copy
   void              Copy(CMCPDState &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMCPDState::CMCPDState(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMCPDState::~CMCPDState(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMCPDState::Copy(CMCPDState &obj)
  {
//--- copy variables
   m_n=obj.m_n;
   m_npairs=obj.m_npairs;
   m_ccnt=obj.m_ccnt;
   m_regterm=obj.m_regterm;
   m_repinneriterationscount=obj.m_repinneriterationscount;
   m_repouteriterationscount=obj.m_repouteriterationscount;
   m_repnfev=obj.m_repnfev;
   m_repterminationtype=obj.m_repterminationtype;
   m_bs.Copy(obj.m_bs);
   m_br.Copy(obj.m_br);
//--- copy arrays
   ArrayCopy(m_states,obj.m_states);
   ArrayCopy(m_ct,obj.m_ct);
   ArrayCopy(m_pw,obj.m_pw);
   ArrayCopy(m_tmpp,obj.m_tmpp);
   ArrayCopy(m_effectivew,obj.m_effectivew);
   ArrayCopy(m_effectivebndl,obj.m_effectivebndl);
   ArrayCopy(m_effectivebndu,obj.m_effectivebndu);
   ArrayCopy(m_effectivect,obj.m_effectivect);
   ArrayCopy(m_h,obj.m_h);
//--- copy matrices
   m_data=obj.m_data;
   m_ec=obj.m_ec;
   m_bndl=obj.m_bndl;
   m_bndu=obj.m_bndu;
   m_c=obj.m_c;
   m_priorp=obj.m_priorp;
   m_effectivec=obj.m_effectivec;
   m_p=obj.m_p;
  }
//+------------------------------------------------------------------+
//| This structure is a MCPD (Markov Chains for Population Data)     |
//| solver.                                                          |
//| You should use ALGLIB functions in order to work with this object|
//+------------------------------------------------------------------+
class CMCPDStateShell
  {
private:
   CMCPDState        m_innerobj;
public:
   //--- constructors, destructor
                     CMCPDStateShell(void);
                     CMCPDStateShell(CMCPDState &obj);
                    ~CMCPDStateShell(void);
   //--- method
   CMCPDState       *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMCPDStateShell::CMCPDStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMCPDStateShell::CMCPDStateShell(CMCPDState &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMCPDStateShell::~CMCPDStateShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMCPDState *CMCPDStateShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| This structure is a MCPD training report:                        |
//|     InnerIterationsCount    -   number of inner iterations of the|
//|                                 underlying optimization algorithm|
//|     OuterIterationsCount    -   number of outer iterations of the|
//|                                 underlying optimization algorithm|
//|     NFEV                    -   number of merit function         |
//|                                 evaluations                      |
//|     TerminationType         -   termination type                 |
//|                                 (same as for MinBLEIC optimizer, |
//|                                 positive values denote success,  | 
//|                                 negative ones - failure)         |
//+------------------------------------------------------------------+
class CMCPDReport
  {
public:
   //--- variables
   int               m_inneriterationscount;
   int               m_outeriterationscount;
   int               m_nfev;
   int               m_terminationtype;
   //--- constructor, destructor
                     CMCPDReport(void);
                    ~CMCPDReport(void);
   //--- copy
   void              Copy(CMCPDReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMCPDReport::CMCPDReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMCPDReport::~CMCPDReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMCPDReport::Copy(CMCPDReport &obj)
  {
//--- copy variables
   m_inneriterationscount=obj.m_inneriterationscount;
   m_outeriterationscount=obj.m_outeriterationscount;
   m_nfev=obj.m_nfev;
   m_terminationtype=obj.m_terminationtype;
  }
//+------------------------------------------------------------------+
//| This structure is a MCPD training report:                        |
//|     InnerIterationsCount    -   number of inner iterations of the|
//|                                 underlying optimization algorithm|
//|     OuterIterationsCount    -   number of outer iterations of the|
//|                                 underlying optimization algorithm|
//|     NFEV                    -   number of merit function         |
//|                                 evaluations                      |
//|     TerminationType         -   termination type                 |
//|                                 (same as for MinBLEIC optimizer, |
//|                                 positive values denote success,  | 
//|                                 negative ones - failure)         |
//+------------------------------------------------------------------+
class CMCPDReportShell
  {
private:
   CMCPDReport       m_innerobj;
public:
   //--- constructors, destructor
                     CMCPDReportShell(void);
                     CMCPDReportShell(CMCPDReport &obj);
                    ~CMCPDReportShell(void);
   //--- methods
   int               GetInnerIterationsCount(void);
   void              SetInnerIterationsCount(const int i);
   int               GetOuterIterationsCount(void);
   void              SetOuterIterationsCount(const int i);
   int               GetNFev(void);
   void              SetNFev(const int i);
   int               GetTerminationType(void);
   void              SetTerminationType(const int i);
   CMCPDReport      *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMCPDReportShell::CMCPDReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMCPDReportShell::CMCPDReportShell(CMCPDReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMCPDReportShell::~CMCPDReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable inneriterationscount           |
//+------------------------------------------------------------------+
int CMCPDReportShell::GetInnerIterationsCount(void)
  {
//--- return result
   return(m_innerobj.m_inneriterationscount);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable inneriterationscount          |
//+------------------------------------------------------------------+
void CMCPDReportShell::SetInnerIterationsCount(const int i)
  {
//--- change value
   m_innerobj.m_inneriterationscount=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable outeriterationscount           |
//+------------------------------------------------------------------+
int CMCPDReportShell::GetOuterIterationsCount(void)
  {
//--- return result
   return(m_innerobj.m_outeriterationscount);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable outeriterationscount          |
//+------------------------------------------------------------------+
void CMCPDReportShell::SetOuterIterationsCount(const int i)
  {
//--- change value
   m_innerobj.m_outeriterationscount=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable nfev                           |
//+------------------------------------------------------------------+
int CMCPDReportShell::GetNFev(void)
  {
//--- return result
   return(m_innerobj.m_nfev);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable nfev                          |
//+------------------------------------------------------------------+
void CMCPDReportShell::SetNFev(const int i)
  {
//--- change value
   m_innerobj.m_nfev=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable terminationtype                |
//+------------------------------------------------------------------+
int CMCPDReportShell::GetTerminationType(void)
  {
//--- return result
   return(m_innerobj.m_terminationtype);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable terminationtype               |
//+------------------------------------------------------------------+
void CMCPDReportShell::SetTerminationType(const int i)
  {
//--- change value
   m_innerobj.m_terminationtype=i;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMCPDReport *CMCPDReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Markov chains class                                              |
//+------------------------------------------------------------------+
class CMarkovCPD
  {
private:
   //--- private method
   static void       MCPDInit(const int n,const int entrystate,const int exitstate,CMCPDState &s);
public:
   //--- constant
   static const double m_xtol;
   //--- constructor, destructor
                     CMarkovCPD(void);
                    ~CMarkovCPD(void);
   //--- public methods
   static void       MCPDCreate(const int n,CMCPDState &s);
   static void       MCPDCreateEntry(const int n,const int entrystate,CMCPDState &s);
   static void       MCPDCreateExit(const int n,const int exitstate,CMCPDState &s);
   static void       MCPDCreateEntryExit(const int n,const int entrystate,const int exitstate,CMCPDState &s);
   static void       MCPDAddTrack(CMCPDState &s,CMatrixDouble &xy,const int k);
   static void       MCPDSetEC(CMCPDState &s,CMatrixDouble &ec);
   static void       MCPDAddEC(CMCPDState &s,const int i,const int j,const double c);
   static void       MCPDSetBC(CMCPDState &s,CMatrixDouble &bndl,CMatrixDouble &bndu);
   static void       MCPDAddBC(CMCPDState &s,const int i,const int j,double bndl,double bndu);
   static void       MCPDSetLC(CMCPDState &s,CMatrixDouble &c,int &ct[],const int k);
   static void       MCPDSetTikhonovRegularizer(CMCPDState &s,const double v);
   static void       MCPDSetPrior(CMCPDState &s,CMatrixDouble &cpp);
   static void       MCPDSetPredictionWeights(CMCPDState &s,double &pw[]);
   static void       MCPDSolve(CMCPDState &s);
   static void       MCPDResults(CMCPDState &s,CMatrixDouble &p,CMCPDReport &rep);
  };
//+------------------------------------------------------------------+
//| Initialize constant                                              |
//+------------------------------------------------------------------+
const double CMarkovCPD::m_xtol=1.0E-8;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMarkovCPD::CMarkovCPD(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMarkovCPD::~CMarkovCPD(void)
  {

  }
//+------------------------------------------------------------------+
//| DESCRIPTION:                                                     |
//| This function creates MCPD (Markov Chains for Population Data)   |
//| solver.                                                          |
//| This solver can be used to find transition matrix P for          |
//| N-dimensional prediction problem where transition from X[i] to   |
//|     X[i+1] is modelled as X[i+1] = P*X[i]                        |
//| where X[i] and X[i+1] are N-dimensional population vectors       |
//| (components of each X are non-negative), and P is a N*N          |
//| transition matrix (elements of   are non-negative, each column   |
//| sums to 1.0).                                                    |
//| Such models arise when when:                                     |
//| * there is some population of individuals                        |
//| * individuals can have different states                          |
//| * individuals can transit from one state to another              |
//| * population size is constant, i.e. there is no new individuals  |
//|   and no one leaves population                                   |
//| * you want to model transitions of individuals from one state    |
//|   into another                                                   |
//| USAGE:                                                           |
//| Here we give very brief outline of the MCPD. We strongly         |
//| recommend you to read examples in the ALGLIB Reference Manual    |
//| and to read ALGLIB User Guide on data analysis which is          |
//| available at http://www.alglib.net/dataanalysis/                 |
//| 1. User initializes algorithm state with MCPDCreate() call       |
//| 2. User adds one or more tracks -  sequences of states which     |
//|    describe evolution of a system being modelled from different  |
//|    starting conditions                                           |
//| 3. User may add optional boundary, equality and/or linear        |
//|    constraints on the coefficients of P by calling one of the    |
//|    following functions:                                          |
//|    * MCPDSetEC() to set equality constraints                     |
//|    * MCPDSetBC() to set bound constraints                        |
//|    * MCPDSetLC() to set linear constraints                       |
//| 4. Optionally, user may set custom weights for prediction errors |
//|    (by default, algorithm assigns non-equal, automatically chosen|
//|    weights for errors in the prediction of different components  |
//|    of X). It can be done with a call of                          |
//|    MCPDSetPredictionWeights() function.                          |
//| 5. User calls MCPDSolve() function which takes algorithm state   |
//|    and pointer (delegate, etc.) to callback function which       |
//|    calculates F/G.                                               |
//| 6. User calls MCPDResults() to get solution                      |
//| INPUT PARAMETERS:                                                |
//|     N       -   problem dimension, N>=1                          |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure stores algorithm state                 |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDCreate(const int n,CMCPDState &s)
  {
//--- check
   if(!CAp::Assert(n>=1,"MCPDCreate: N<1"))
      return;
//--- function call
   MCPDInit(n,-1,-1,s);
  }
//+------------------------------------------------------------------+
//| DESCRIPTION:                                                     |
//| This function is a specialized version of MCPDCreate() function, |
//| and we recommend  you  to read comments for this function for    |
//| general information about MCPD solver.                           |
//| This function creates MCPD (Markov Chains for Population Data)   |
//| solver for "Entry-state" model, i.e. model where transition from |
//| X[i] to X[i+1] is modelled as                                    |
//|     X[i+1] = P*X[i]                                              |
//| where                                                            |
//|     X[i] and X[i+1] are N-dimensional state vectors              |
//|     P is a N*N transition matrix                                 |
//| and  one  selected component of X[] is called "entry" state and  |
//| is treated in a special way:                                     |
//|     system state always transits from "entry" state to some      |
//|     another state                                                |
//|     system state can not transit from any state into "entry"     |
//|     state                                                        |
//| Such conditions basically mean that row of P which corresponds to|
//| "entry" state is zero.                                           |
//| Such models arise when:                                          |
//| * there is some population of individuals                        |
//| * individuals can have different states                          |
//| * individuals can transit from one state to another              |
//| * population size is NOT constant -  at every moment of time     |
//|   there is some (unpredictable) amount of "new" individuals,     |
//|   which can transit into one of the states at the next turn, but |
//|   still no one leaves population                                 |
//| * you want to model transitions of individuals from one state    |
//|   into another                                                   |
//| * but you do NOT want to predict amount of "new" individuals     |
//|   because it does not depends on individuals already present     |
//|   (hence system can not transit INTO entry state - it can only   |
//|   transit FROM it).                                              |
//| This model is discussed in more details in the ALGLIB User Guide |
//| (see http://www.alglib.net/dataanalysis/ for more data).         |
//| INPUT PARAMETERS:                                                |
//|     N       -   problem dimension, N>=2                          |
//|     EntryState- index of entry state, in 0..N-1                  |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure stores algorithm state                 |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDCreateEntry(const int n,const int entrystate,
                                        CMCPDState &s)
  {
//--- check
   if(!CAp::Assert(n>=2,__FUNCTION__+": N<2"))
      return;
//--- check
   if(!CAp::Assert(entrystate>=0,__FUNCTION__+": EntryState<0"))
      return;
//--- check
   if(!CAp::Assert(entrystate<n,__FUNCTION__+": EntryState>=N"))
      return;
//--- function call
   MCPDInit(n,entrystate,-1,s);
  }
//+------------------------------------------------------------------+
//| DESCRIPTION:                                                     |
//| This function is a specialized version of MCPDCreate() function, |
//| and we recommend  you  to read comments for this function for    |
//| general information about MCPD solver.                           |
//| This function creates MCPD (Markov Chains for Population Data)   |
//| solver for "Exit-state" model, i.e. model where transition from  |
//| X[i] to X[i+1] is modelled as                                    |
//|     X[i+1] = P*X[i]                                              |
//| where                                                            |
//|     X[i] and X[i+1] are N-dimensional state vectors              |
//|     P is a N*N transition matrix                                 |
//| and  one  selected component of X[] is called "exit" state and   |
//| is treated in a special way:                                     |
//|     system state can transit from any state into "exit" state    |
//|     system state can not transit from "exit" state into any other|
//|     state transition operator discards "exit" state (makes it    |
//|     zero at each turn)                                           |
//| Such conditions basically mean that column of P which            |
//| corresponds to "exit" state is zero. Multiplication by such P    |
//| may decrease sum of vector components.                           |
//| Such models arise when:                                          |
//| * there is some population of individuals                        |
//| * individuals can have different states                          |
//| * individuals can transit from one state to another              |
//| * population size is NOT constant - individuals can move into    |
//|   "exit" state and leave population at the next turn, but there  |
//|   are no new individuals                                         |
//| * amount of individuals which leave population can be predicted  |
//| * you want to model transitions of individuals from one state    |
//|   into another (including transitions into the "exit" state)     |
//| This model is discussed in more details in the ALGLIB User Guide |
//| (see http://www.alglib.net/dataanalysis/ for more data).         |
//| INPUT PARAMETERS:                                                |
//|     N       -   problem dimension, N>=2                          |
//|     ExitState-  index of exit state, in 0..N-1                   |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure stores algorithm state                 |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDCreateExit(const int n,const int exitstate,
                                       CMCPDState &s)
  {
//--- check
   if(!CAp::Assert(n>=2,__FUNCTION__+": N<2"))
      return;
//--- check
   if(!CAp::Assert(exitstate>=0,__FUNCTION__+": ExitState<0"))
      return;
//--- check
   if(!CAp::Assert(exitstate<n,__FUNCTION__+": ExitState>=N"))
      return;
//--- function call
   MCPDInit(n,-1,exitstate,s);
  }
//+------------------------------------------------------------------+
//| DESCRIPTION:                                                     |
//| This function is a specialized version of MCPDCreate() function, |
//| and we recommend you to read comments for this function for      |
//| general information about MCPD solver.                           |
//| This function creates MCPD (Markov Chains for Population Data)   |
//| solver for "Entry-Exit-states" model, i.e. model where transition|
//| from X[i] to X[i+1] is modelled as                               |
//|     X[i+1] = P*X[i]                                              |
//| where                                                            |
//|     X[i] and X[i+1] are N-dimensional state vectors              |
//|     P is a N*N transition matrix                                 |
//| one selected component of X[] is called "entry" state and is a   |
//| treated in special way:                                          |
//|     system state always transits from "entry" state to some      |
//|     another state                                                |
//|     system state can not transit from any state into "entry"     |
//|     state                                                        |
//| and another one component of X[] is called "exit" state and is   |
//| treated in a special way too:                                    |
//|     system state can transit from any state into "exit" state    |
//|     system state can not transit from "exit" state into any other|
//|     state transition operator discards "exit" state (makes it    |
//|     zero at each turn)                                           |
//| Such conditions basically mean that:                             |
//|     row of P which corresponds to "entry" state is zero          |
//|     column of P which corresponds to "exit" state is zero        |
//| Multiplication by such P may decrease sum of vector components.  |
//| Such models arise when:                                          |
//| * there is some population of individuals                        |
//| * individuals can have different states                          |
//| * individuals can transit from one state to another              |
//| * population size is NOT constant                                |
//| * at every moment of time there is some (unpredictable) amount   |
//|   of "new" individuals, which can transit into one of the states |
//|   at the next turn                                               |
//| * some individuals can move (predictably) into "exit" state      |
//|   and leave population at the next turn                          |
//| * you want to model transitions of individuals from one state    |
//|   into another, including transitions from the "entry" state and |
//|   into the "exit" state.                                         |
//| * but you do NOT want to predict amount of "new" individuals     |
//|   because it does not depends on individuals already present     |
//|   (hence system can not transit INTO entry state - it can only   |
//|   transit FROM it).                                              |
//| This model is discussed  in  more  details  in  the ALGLIB User  |
//| Guide (see http://www.alglib.net/dataanalysis/ for more data).   |
//| INPUT PARAMETERS:                                                |
//|     N       -   problem dimension, N>=2                          |
//|     EntryState- index of entry state, in 0..N-1                  |
//|     ExitState-  index of exit state, in 0..N-1                   |
//| OUTPUT PARAMETERS:                                               |
//|     State   -   structure stores algorithm state                 |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDCreateEntryExit(const int n,const int entrystate,
                                            const int exitstate,CMCPDState &s)
  {
//--- check
   if(!CAp::Assert(n>=2,__FUNCTION__+": N<2"))
      return;
//--- check
   if(!CAp::Assert(entrystate>=0,__FUNCTION__+": EntryState<0"))
      return;
//--- check
   if(!CAp::Assert(entrystate<n,__FUNCTION__+": EntryState>=N"))
      return;
//--- check
   if(!CAp::Assert(exitstate>=0,__FUNCTION__+": ExitState<0"))
      return;
//--- check
   if(!CAp::Assert(exitstate<n,__FUNCTION__+": ExitState>=N"))
      return;
//--- check
   if(!CAp::Assert(entrystate!=exitstate,__FUNCTION__+": EntryState=ExitState"))
      return;
//--- function call
   MCPDInit(n,entrystate,exitstate,s);
  }
//+------------------------------------------------------------------+
//| This function is used to add a track - sequence of system states |
//| at the different moments of its evolution.                       |
//| You may add one or several tracks to the MCPD solver. In case you|
//| have several tracks, they won't overwrite each other. For        |
//| example, if you pass two tracks, A1-A2-A3 (system at t=A+1, t=A+2|
//| and t=A+3) and B1-B2-B3, then solver will try to model           |
//| transitions from t=A+1 to t=A+2, t=A+2 to t=A+3, t=B+1 to t=B+2, |
//| t=B+2 to t=B+3. But it WONT mix these two tracks - i.e. it wont  |
//| try to model transition from t=A+3 to t=B+1.                     |
//| INPUT PARAMETERS:                                                |
//|     S       -   solver                                           |
//|     XY      -   track, array[K,N]:                               |
//|                 * I-th row is a state at t=I                     |
//|                 * elements of XY must be non-negative (exception |
//|                   will be thrown on negative elements)           |
//|     K       -   number of points in a track                      |
//|                 * if given, only leading K rows of XY are used   |
//|                 * if not given, automatically determined from    |
//|                   size of XY                                     |
//| NOTES:                                                           |
//| 1. Track may contain either proportional or population data:     |
//|    * with proportional data all rows of XY must sum to 1.0, i.e. |
//|      we have proportions instead of absolute population values   |
//|    * with population data rows of XY contain population counts   |
//|      and generally do not sum to 1.0 (although they still must be|
//|      non-negative)                                               |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDAddTrack(CMCPDState &s,CMatrixDouble &xy,
                                     const int k)
  {
//--- create variables
   int    i=0;
   int    j=0;
   int    n=0;
   double s0=0;
   double s1=0;
//--- initialization
   n=s.m_n;
//--- check
   if(!CAp::Assert(k>=0,__FUNCTION__+": K<0"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(xy)>=n,__FUNCTION__+": Cols(XY)<N"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(xy)>=k,__FUNCTION__+": Rows(XY)<K"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(xy,k,n),__FUNCTION__+": XY contains infinite or NaN elements"))
      return;
   for(i=0;i<=k-1;i++)
     {
      for(j=0;j<=n-1;j++)
        {
         //--- check
         if(!CAp::Assert((double)(xy[i][j])>=0.0,__FUNCTION__+": XY contains negative elements"))
            return;
        }
     }
//--- check
   if(k<2)
      return;
//--- check
   if(CAp::Rows(s.m_data)<s.m_npairs+k-1)
      CApServ::RMatrixResize(s.m_data,MathMax(2*CAp::Rows(s.m_data),s.m_npairs+k-1),2*n);
//--- calculation
   for(i=0;i<=k-2;i++)
     {
      s0=0;
      s1=0;
      for(j=0;j<=n-1;j++)
        {
         //--- check
         if(s.m_states[j]>=0)
            s0=s0+xy[i][j];
         //--- check
         if(s.m_states[j]<=0)
            s1=s1+xy[i+1][j];
        }
      //--- check
      if(s0>0.0 && s1>0.0)
        {
         for(j=0;j<=n-1;j++)
           {
            //--- check
            if(s.m_states[j]>=0)
               s.m_data[s.m_npairs].Set(j,xy[i][j]/s0);
            else
               s.m_data[s.m_npairs].Set(j,0.0);
            //--- check
            if(s.m_states[j]<=0)
               s.m_data[s.m_npairs].Set(n+j,xy[i+1][j]/s1);
            else
               s.m_data[s.m_npairs].Set(n+j,0.0);
           }
         //--- change value
         s.m_npairs=s.m_npairs+1;
        }
     }
  }
//+------------------------------------------------------------------+
//| This function is used to add equality constraints on the elements|
//| of the transition matrix P.                                      |
//| MCPD solver has four types of constraints which can be placed    |
//| on P:                                                            |
//| * user-specified equality constraints (optional)                 |
//| * user-specified bound constraints (optional)                    |
//| * user-specified general linear constraints (optional)           |
//| * basic constraints (always present):                            |
//|   * non-negativity: P[i,j]>=0                                    |
//|   * consistency: every column of P sums to 1.0                   |
//| Final constraints which are passed to the underlying optimizer   |
//| are calculated as intersection of all present constraints. For   |
//| example, you may specify boundary constraint on P[0,0] and       |
//| equality one:                                                    |
//|     0.1<=P[0,0]<=0.9                                             |
//|     P[0,0]=0.5                                                   |
//| Such combination of constraints will be silently reduced to their|
//| intersection, which is P[0,0]=0.5.                               |
//| This function can be used to place equality constraints on       |
//| arbitrary subset of elements of P. Set of constraints is         |
//| specified by EC, which may contain either NAN's or finite numbers|
//| from [0,1]. NAN denotes absence of constraint, finite number     |
//| denotes equality constraint on specific element of P.            |
//| You can also use MCPDAddEC() function which allows to ADD        |
//| equality constraint for one element of P without changing        |
//| constraints for other elements.                                  |
//| These functions (MCPDSetEC and MCPDAddEC) interact as follows:   |
//| * there is internal matrix of equality constraints which is      |
//|   stored in the MCPD solver                                      |
//| * MCPDSetEC() replaces this matrix by another one (SET)          |
//| * MCPDAddEC() modifies one element of this matrix and leaves     |
//|   other ones unchanged (ADD)                                     |
//| * thus MCPDAddEC() call preserves all modifications done by      |
//|   previous calls, while MCPDSetEC() completely discards all      |
//|   changes done to the equality constraints.                      |
//| INPUT PARAMETERS:                                                |
//|     S       -   solver                                           |
//|     EC      -   equality constraints, array[N,N]. Elements of EC |
//|                 can be either NAN's or finite numbers from [0,1].|
//|                 NAN denotes absence of constraints, while finite |
//|                 value denotes equality constraint on the         |
//|                 corresponding element of P.                      |
//| NOTES:                                                           |
//| 1. infinite values of EC will lead to exception being thrown.    |
//| Values less than 0.0 or greater than 1.0 will lead to error code |
//| being returned after call to MCPDSolve().                        |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDSetEC(CMCPDState &s,CMatrixDouble &ec)
  {
//--- create variables
   int i=0;
   int j=0;
   int n=0;
//--- initialization
   n=s.m_n;
//--- check
   if(!CAp::Assert(CAp::Cols(ec)>=n,__FUNCTION__+": Cols(EC)<N"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(ec)>=n,__FUNCTION__+": Rows(EC)<N"))
      return;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      for(j=0;j<=n-1;j++)
        {
         //--- check
         if(!CAp::Assert(CMath::IsFinite(ec[i][j]) || CInfOrNaN::IsNaN(ec[i][j]),"MCPDSetEC: EC containts infinite elements"))
            return;
         s.m_ec[i].Set(j,ec[i][j]);
        }
     }
  }
//+------------------------------------------------------------------+
//| This function is used to add equality constraints on the elements|
//| of the transition matrix P.                                      |
//| MCPD solver has four types of constraints which can be placed    |
//| on P:                                                            |
//| * user-specified equality constraints (optional)                 |
//| * user-specified bound constraints (optional)                    |
//| * user-specified general linear constraints (optional)           |
//| * basic constraints (always present):                            |
//|   * non-negativity: P[i,j]>=0                                    |
//|   * consistency: every column of P sums to 1.0                   |
//| Final constraints which are passed to the underlying optimizer   |
//| are calculated as intersection of all present constraints. For   |
//| example, you may specify boundary constraint on P[0,0] and       |
//| equality one:                                                    |
//|     0.1<=P[0,0]<=0.9                                             |
//|     P[0,0]=0.5                                                   |
//| Such combination of constraints will be silently reduced to their|
//| intersection, which is P[0,0]=0.5.                               |
//| This function can be used to ADD equality constraint for one     |
//| element of P without changing constraints for other elements.    |
//| You can also use MCPDSetEC() function which allows you to specify|
//| arbitrary set of equality constraints in one call.               |
//| These functions (MCPDSetEC and MCPDAddEC) interact as follows:   |
//| * there is internal matrix of equality constraints which is      |
//|   stored in the MCPD solver                                      |
//| * MCPDSetEC() replaces this matrix by another one (SET)          |
//| * MCPDAddEC() modifies one element of this matrix and leaves     |
//|   other ones unchanged (ADD)                                     |
//| * thus MCPDAddEC() call preserves all modifications done by      |
//|   previous calls, while MCPDSetEC() completely discards all      |
//|   changes done to the equality constraints.                      |
//| INPUT PARAMETERS:                                                |
//|     S       -   solver                                           |
//|     I       -   row index of element being constrained           |
//|     J       -   column index of element being constrained        |
//|     C       -   value (constraint for P[I,J]). Can be either NAN |
//|                 (no constraint) or finite value from [0,1].      |
//| NOTES:                                                           |
//| 1. infinite values of C will lead to exception being thrown.     |
//| Values less than 0.0 or greater than 1.0 will lead to error code |
//| being returned after call to MCPDSolve().                        |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDAddEC(CMCPDState &s,const int i,const int j,
                                  const double c)
  {
//--- check
   if(!CAp::Assert(i>=0,__FUNCTION__+": I<0"))
      return;
//--- check
   if(!CAp::Assert(i<s.m_n,__FUNCTION__+": I>=N"))
      return;
//--- check
   if(!CAp::Assert(j>=0,__FUNCTION__+": J<0"))
      return;
//--- check
   if(!CAp::Assert(j<s.m_n,__FUNCTION__+": J>=N"))
      return;
//--- check
   if(!CAp::Assert(CInfOrNaN::IsNaN(c) || CMath::IsFinite(c),"MCPDAddEC: C is not finite number or NAN"))
      return;
   s.m_ec[i].Set(j,c);
  }
//+------------------------------------------------------------------+
//| This function is used to add bound constraints on the elements   |
//| of the transition matrix P.                                      |
//| MCPD solver has four types of constraints which can be placed    |
//| on P:                                                            |
//| * user-specified equality constraints (optional)                 |
//| * user-specified bound constraints (optional)                    |
//| * user-specified general linear constraints (optional)           |
//| * basic constraints (always present):                            |
//|   * non-negativity: P[i,j]>=0                                    |
//|   * consistency: every column of P sums to 1.0                   |
//| Final constraints which are passed to the underlying optimizer   |
//| are calculated as intersection of all present constraints. For   |
//| example, you may specify boundary constraint on P[0,0] and       |
//| equality one:                                                    |
//|     0.1<=P[0,0]<=0.9                                             |
//|     P[0,0]=0.5                                                   |
//| Such combination of constraints will be silently reduced to their|
//| intersection, which is P[0,0]=0.5.                               |
//| This function can be used to place bound constraints on arbitrary|
//| subset of elements of P. Set of constraints is specified by      |
//| BndL/BndU matrices, which may contain arbitrary combination of   |
//| finite numbers or infinities (like -INF<x<=0.5 or 0.1<=x<+INF).  |
//| You can also use MCPDAddBC() function which allows to ADD bound  |
//| constraint for one element of P without changing constraints for |
//| other elements.                                                  |
//| These functions (MCPDSetBC and MCPDAddBC) interact as follows:   |
//| * there is internal matrix of bound constraints which is stored  |
//|   in the MCPD solver                                             |
//| * MCPDSetBC() replaces this matrix by another one (SET)          |
//| * MCPDAddBC() modifies one element of this matrix and leaves     |
//|   other ones unchanged (ADD)                                     |
//| * thus MCPDAddBC() call preserves all modifications done by      |
//|   previous calls, while MCPDSetBC() completely discards all      |
//|   changes done to the equality constraints.                      |
//| INPUT PARAMETERS:                                                |
//|     S       -   solver                                           |
//|     BndL    -   lower bounds constraints, array[N,N]. Elements of|
//|                 BndL can be finite numbers or -INF.              |
//|     BndU    -   upper bounds constraints, array[N,N]. Elements of|
//|                 BndU can be finite numbers or +INF.              |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDSetBC(CMCPDState &s,CMatrixDouble &bndl,
                                  CMatrixDouble &bndu)
  {
//--- create variables
   int i=0;
   int j=0;
   int n=0;
//--- initialization
   n=s.m_n;
//--- check
   if(!CAp::Assert(CAp::Cols(bndl)>=n,__FUNCTION__+": Cols(BndL)<N"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(bndl)>=n,__FUNCTION__+": Rows(BndL)<N"))
      return;
//--- check
   if(!CAp::Assert(CAp::Cols(bndu)>=n,__FUNCTION__+": Cols(BndU)<N"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(bndu)>=n,__FUNCTION__+": Rows(BndU)<N"))
      return;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      for(j=0;j<=n-1;j++)
        {
         //--- check
         if(!CAp::Assert(CMath::IsFinite(bndl[i][j]) || CInfOrNaN::IsNegativeInfinity(bndl[i][j]),"MCPDSetBC: BndL containts NAN or +INF"))
            return;
         //--- check
         if(!CAp::Assert(CMath::IsFinite(bndu[i][j]) || CInfOrNaN::IsPositiveInfinity(bndu[i][j]),"MCPDSetBC: BndU containts NAN or -INF"))
            return;
         //--- change values
         s.m_bndl[i].Set(j,bndl[i][j]);
         s.m_bndu[i].Set(j,bndu[i][j]);
        }
     }
  }
//+------------------------------------------------------------------+
//| This function is used to add bound constraints on the elements   |
//| of the transition matrix P.                                      |
//| MCPD solver has four types of constraints which can be placed    |
//| on P:                                                            |
//| * user-specified equality constraints (optional)                 |
//| * user-specified bound constraints (optional)                    |
//| * user-specified general linear constraints (optional)           |
//| * basic constraints (always present):                            |
//|   * non-negativity: P[i,j]>=0                                    |
//|   * consistency: every column of P sums to 1.0                   |
//| Final constraints which are passed to the underlying optimizer   |
//| are calculated as intersection of all present constraints. For   |
//| example, you may specify boundary constraint on P[0,0] and       |
//| equality one:                                                    |
//|     0.1<=P[0,0]<=0.9                                             |
//|     P[0,0]=0.5                                                   |
//| Such combination of constraints will be silently reduced to their|
//| intersection, which is P[0,0]=0.5.                               |
//| This function can be used to ADD bound constraint for one element|
//| of P without changing constraints for other elements.            |
//| You can also use MCPDSetBC() function which allows to place bound|
//| constraints on arbitrary subset of elements of P. Set of         |
//| constraints is specified  by  BndL/BndU matrices, which may      |
//| contain arbitrary combination of finite numbers or infinities    |
//| (like -INF<x<=0.5 or 0.1<=x<+INF).                               |
//| These functions (MCPDSetBC and MCPDAddBC) interact as follows:   |
//| * there is internal matrix of bound constraints which is stored  |
//|   in the MCPD solver                                             |
//| * MCPDSetBC() replaces this matrix by another one (SET)          |
//| * MCPDAddBC() modifies one element of this matrix and leaves     |
//|   other ones unchanged (ADD)                                     |
//| * thus MCPDAddBC() call preserves all modifications done by      |
//|   previous calls, while MCPDSetBC() completely discards all      |
//|   changes done to the equality constraints.                      |
//| INPUT PARAMETERS:                                                |
//|     S       -   solver                                           |
//|     I       -   row index of element being constrained           |
//|     J       -   column index of element being constrained        |
//|     BndL    -   lower bound                                      |
//|     BndU    -   upper bound                                      |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDAddBC(CMCPDState &s,const int i,const int j,
                                  double bndl,double bndu)
  {
//--- check
   if(!CAp::Assert(i>=0,__FUNCTION__+": I<0"))
      return;
//--- check
   if(!CAp::Assert(i<s.m_n,__FUNCTION__+": I>=N"))
      return;
//--- check
   if(!CAp::Assert(j>=0,__FUNCTION__+": J<0"))
      return;
//--- check
   if(!CAp::Assert(j<s.m_n,__FUNCTION__+": J>=N"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(bndl) || CInfOrNaN::IsNegativeInfinity(bndl),"MCPDAddBC: BndL is NAN or +INF"))
      return;
//--- check
   if(!CAp::Assert(CMath::IsFinite(bndu) || CInfOrNaN::IsPositiveInfinity(bndu),"MCPDAddBC: BndU is NAN or -INF"))
      return;
//--- change values
   s.m_bndl[i].Set(j,bndl);
   s.m_bndu[i].Set(j,bndu);
  }
//+------------------------------------------------------------------+
//| This function is used to set linear equality/inequality          |
//| constraints on the elements of the transition matrix P.          |
//| This function can be used to set one or several general linear   |
//| constraints on the elements of P. Two types of constraints are   |
//| supported:                                                       |
//| * equality constraints                                           |
//| * inequality constraints (both less-or-equal and                 |
//|   greater-or-equal)                                              |
//| Coefficients of constraints are specified by matrix C (one of the|
//| parameters). One row of C corresponds to one constraint.         |
//| Because transition matrix P has N*N elements, we need N*N columns|
//| to store all coefficients  (they  are  stored row by row), and   |
//| one more column to store right part - hence C has N*N+1 columns. |
//| Constraint kind is stored in the CT array.                       |
//| Thus, I-th linear constraint is                                  |
//|     P[0,0]*C[I,0] + P[0,1]*C[I,1] + .. + P[0,N-1]*C[I,N-1] +     |
//|         + P[1,0]*C[I,N] + P[1,1]*C[I,N+1] + ... +                |
//|         + P[N-1,N-1]*C[I,N*N-1]  ?=?  C[I,N*N]                   |
//| where ?=? can be either "=" (CT[i]=0), "<=" (CT[i]<0) or ">="    |
//| (CT[i]>0).                                                       |
//| Your constraint may involve only some subset of P (less than N*N |
//| elements).                                                       |
//| For example it can be something like                             |
//|     P[0,0] + P[0,1] = 0.5                                        |
//| In this case you still should pass matrix  with N*N+1 columns,   |
//| but all its elements (except for C[0,0], C[0,1] and C[0,N*N-1])  |
//| will be zero.                                                    |
//| INPUT PARAMETERS:                                                |
//|     S       -   solver                                           |
//|     C       -   array[K,N*N+1] - coefficients of constraints     |
//|                 (see above for complete description)             |
//|     CT      -   array[K] - constraint types                      |
//|                 (see above for complete description)             |
//|     K       -   number of equality/inequality constraints, K>=0: |
//|                 * if given, only leading K elements of C/CT are  |
//|                   used                                           |
//|                 * if not given, automatically determined from    |
//|                   sizes of C/CT                                  |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDSetLC(CMCPDState &s,CMatrixDouble &c,int &ct[],
                                  const int k)
  {
//--- create variables
   int i=0;
   int j=0;
   int n=0;
//--- initialization
   n=s.m_n;
//--- check
   if(!CAp::Assert(CAp::Cols(c)>=n*n+1,__FUNCTION__+": Cols(C)<N*N+1"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(c)>=k,__FUNCTION__+": Rows(C)<K"))
      return;
//--- check
   if(!CAp::Assert(CAp::Len(ct)>=k,__FUNCTION__+": Len(CT)<K"))
      return;
//--- check
   if(!CAp::Assert(CApServ::IsFiniteMatrix(c,k,n*n+1),__FUNCTION__+": C contains infinite or NaN values!"))
      return;
//--- function call
   CApServ::RMatrixSetLengthAtLeast(s.m_c,k,n*n+1);
//--- function call
   CApServ::IVectorSetLengthAtLeast(s.m_ct,k);
//--- calculation
   for(i=0;i<=k-1;i++)
     {
      for(j=0;j<=n*n;j++)
         s.m_c[i].Set(j,c[i][j]);
      s.m_ct[i]=ct[i];
     }
   s.m_ccnt=k;
  }
//+------------------------------------------------------------------+
//| This function allows to tune amount of Tikhonov regularization   |
//| being applied to your problem.                                   |
//| By default, regularizing term is equal to r*||P-prior_P||^2,     |
//| where r is a small non-zero value,  P is transition matrix,      |
//| prior_P is identity matrix, ||X||^2 is a sum of squared elements |
//| of X.                                                            |
//| This function allows you to change coefficient r. You can also   |
//| change prior values with MCPDSetPrior() function.                |
//| INPUT PARAMETERS:                                                |
//|     S      -   solver                                            |
//|     V      -   regularization  coefficient, finite non-negative  |
//|                value. It is not recommended to specify zero      |
//|                value unless you are pretty sure that you want it.|
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDSetTikhonovRegularizer(CMCPDState &s,const double v)
  {
//--- check
   if(!CAp::Assert(CMath::IsFinite(v),__FUNCTION__+": V is infinite or NAN"))
      return;
//--- check
   if(!CAp::Assert(v>=0.0,__FUNCTION__+": V is less than zero"))
      return;
//--- change value
   s.m_regterm=v;
  }
//+------------------------------------------------------------------+
//| This function allows to set prior values used for regularization |
//| of your problem.                                                 |
//| By default, regularizing term is equal to r*||P-prior_P||^2,     |
//| where r is a small non-zero value,  P is transition matrix,      |
//| prior_P is identity matrix, ||X||^2 is a sum of squared elements |
//| of X.                                                            |
//| This function allows you to change prior values prior_P. You can |
//| also change r with MCPDSetTikhonovRegularizer() function.        |
//| INPUT PARAMETERS:                                                |
//|     S       -   solver                                           |
//|     PP      -   array[N,N], matrix of prior values:              |
//|                 1. elements must be real numbers from [0,1]      |
//|                 2. columns must sum to 1.0.                      |
//|                 First property is checked (exception is thrown   |
//|                 otherwise), while second one is not              |
//|                 checked/enforced.                                |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDSetPrior(CMCPDState &s,CMatrixDouble &cpp)
  {
//--- create variables
   int i=0;
   int j=0;
   int n=0;
//--- create copy of matrices
   CMatrixDouble pp;
   pp=cpp;
//--- initialization
   n=s.m_n;
//--- check
   if(!CAp::Assert(CAp::Cols(pp)>=n,__FUNCTION__+": Cols(PP)<N"))
      return;
//--- check
   if(!CAp::Assert(CAp::Rows(pp)>=n,__FUNCTION__+": Rows(PP)<K"))
      return;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      for(j=0;j<=n-1;j++)
        {
         //--- check
         if(!CAp::Assert(CMath::IsFinite(pp[i][j]),__FUNCTION__+": PP containts infinite elements"))
            return;
         //--- check
         if(!CAp::Assert(pp[i][j]>=0.0 && pp[i][j]<=1.0,__FUNCTION__+": PP[i][j] is less than 0.0 or greater than 1.0"))
            return;
         //--- change value
         s.m_priorp[i].Set(j,pp[i][j]);
        }
     }
  }
//+------------------------------------------------------------------+
//| This function is used to change prediction weights               |
//| MCPD solver scales prediction errors as follows                  |
//|     Error(P) = ||W*(y-P*x)||^2                                   |
//| where                                                            |
//|     x is a system state at time t                                |
//|     y is a system state at time t+1                              |
//|     P is a transition matrix                                     |
//|     W is a diagonal scaling matrix                               |
//| By default, weights are chosen in order to minimize relative     |
//| prediction error instead of absolute one. For example, if one    |
//| component of state is about 0.5 in magnitude and another one is  |
//| about 0.05, then algorithm will make corresponding weights equal |
//| to 2.0 and 20.0.                                                 |
//| INPUT PARAMETERS:                                                |
//|     S       -   solver                                           |
//|     PW      -   array[N], weights:                               |
//|                 * must be non-negative values (exception will be |
//|                 thrown otherwise)                                |
//|                 * zero values will be replaced by automatically  |
//|                 chosen values                                    |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDSetPredictionWeights(CMCPDState &s,double &pw[])
  {
//--- create variables
   int i=0;
   int n=0;
//--- initialization
   n=s.m_n;
//--- check
   if(!CAp::Assert(CAp::Len(pw)>=n,__FUNCTION__+": Length(PW)<N"))
      return;
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(!CAp::Assert(CMath::IsFinite(pw[i]),__FUNCTION__+": PW containts infinite or NAN elements"))
         return;
      //--- check
      if(!CAp::Assert(pw[i]>=0.0,__FUNCTION__+": PW containts negative elements"))
         return;
      //--- change value
      s.m_pw[i]=pw[i];
     }
  }
//+------------------------------------------------------------------+
//| This function is used to start solution of the MCPD problem.     |
//| After return from this function, you can use MCPDResults() to get|
//| solution and completion code.                                    |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDSolve(CMCPDState &s)
  {
//--- create variables
   int    n=0;
   int    npairs=0;
   int    ccnt=0;
   int    i=0;
   int    j=0;
   int    k=0;
   int    k2=0;
   double v=0;
   double vv=0;
   int    i_=0;
   int    i1_=0;
//--- initialization
   n=s.m_n;
   npairs=s.m_npairs;
//--- init fields of S
   s.m_repterminationtype=0;
   s.m_repinneriterationscount=0;
   s.m_repouteriterationscount=0;
   s.m_repnfev=0;
   for(k=0;k<=n-1;k++)
     {
      for(k2=0;k2<=n-1;k2++)
         s.m_p[k].Set(k2,CInfOrNaN::NaN());
     }
//--- Generate "effective" weights for prediction and calculate preconditioner
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(s.m_pw[i]==0.0)
        {
         //--- change values
         v=0;
         k=0;
         for(j=0;j<=npairs-1;j++)
           {
            //--- check
            if(s.m_data[j][n+i]!=0.0)
              {
               v=v+s.m_data[j][n+i];
               k=k+1;
              }
           }
         //--- check
         if(k!=0)
            s.m_effectivew[i]=k/v;
         else
            s.m_effectivew[i]=1.0;
        }
      else
         s.m_effectivew[i]=s.m_pw[i];
     }
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      for(j=0;j<=n-1;j++)
         s.m_h[i*n+j]=2*s.m_regterm;
     }
//--- calculation
   for(k=0;k<=npairs-1;k++)
     {
      for(i=0;i<=n-1;i++)
        {
         for(j=0;j<=n-1;j++)
            s.m_h[i*n+j]=s.m_h[i*n+j]+2*CMath::Sqr(s.m_effectivew[i])*CMath::Sqr(s.m_data[k][j]);
        }
     }
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      for(j=0;j<=n-1;j++)
        {
         //--- check
         if(s.m_h[i*n+j]==0.0)
            s.m_h[i*n+j]=1;
        }
     }
//--- Generate "effective" BndL/BndU
   for(i=0;i<=n-1;i++)
     {
      for(j=0;j<=n-1;j++)
        {
         //--- Set default boundary constraints.
         //--- Lower bound is always zero,upper bound is calculated
         //--- with respect to entry/exit states.
         s.m_effectivebndl[i*n+j]=0.0;
         //--- check
         if(s.m_states[i]>0 || s.m_states[j]<0)
            s.m_effectivebndu[i*n+j]=0.0;
         else
            s.m_effectivebndu[i*n+j]=1.0;
         //--- Calculate intersection of the default and user-specified bound constraints.
         //--- This code checks consistency of such combination.
         if(CMath::IsFinite(s.m_bndl[i][j]) && s.m_bndl[i][j]>s.m_effectivebndl[i*n+j])
            s.m_effectivebndl[i*n+j]=s.m_bndl[i][j];
         //--- check
         if(CMath::IsFinite(s.m_bndu[i][j]) && s.m_bndu[i][j]<s.m_effectivebndu[i*n+j])
            s.m_effectivebndu[i*n+j]=s.m_bndu[i][j];
         //--- check
         if(s.m_effectivebndl[i*n+j]>s.m_effectivebndu[i*n+j])
           {
            s.m_repterminationtype=-3;
            return;
           }
         //--- Calculate intersection of the effective bound constraints
         //--- and user-specified equality constraints.
         //--- This code checks consistency of such combination.
         if(CMath::IsFinite(s.m_ec[i][j]))
           {
            //--- check
            if(s.m_ec[i][j]<s.m_effectivebndl[i*n+j] || s.m_ec[i][j]>s.m_effectivebndu[i*n+j])
              {
               s.m_repterminationtype=-3;
               return;
              }
            //--- change values
            s.m_effectivebndl[i*n+j]=s.m_ec[i][j];
            s.m_effectivebndu[i*n+j]=s.m_ec[i][j];
           }
        }
     }
//--- Generate linear constraints:
//--- * "default" sums-to-one constraints (not generated for "exit" states)
   CApServ::RMatrixSetLengthAtLeast(s.m_effectivec,s.m_ccnt+n,n*n+1);
//--- function call
   CApServ::IVectorSetLengthAtLeast(s.m_effectivect,s.m_ccnt+n);
   ccnt=s.m_ccnt;
   for(i=0;i<=s.m_ccnt-1;i++)
     {
      for(j=0;j<=n*n;j++)
         s.m_effectivec[i].Set(j,s.m_c[i][j]);
      s.m_effectivect[i]=s.m_ct[i];
     }
//--- calculation
   for(i=0;i<=n-1;i++)
     {
      //--- check
      if(s.m_states[i]>=0)
        {
         for(k=0;k<=n*n-1;k++)
            s.m_effectivec[ccnt].Set(k,0);
         for(k=0;k<=n-1;k++)
            s.m_effectivec[ccnt].Set(k*n+i,1);
         //--- change values
         s.m_effectivec[ccnt].Set(n*n,1.0);
         s.m_effectivect[ccnt]=0;
         ccnt=ccnt+1;
        }
     }
//--- create optimizer
   for(i=0;i<=n-1;i++)
     {
      for(j=0;j<=n-1;j++)
         s.m_tmpp[i*n+j]=1.0/(double)n;
     }
//--- function calls
   CMinBLEIC::MinBLEICRestartFrom(s.m_bs,s.m_tmpp);
   CMinBLEIC::MinBLEICSetBC(s.m_bs,s.m_effectivebndl,s.m_effectivebndu);
   CMinBLEIC::MinBLEICSetLC(s.m_bs,s.m_effectivec,s.m_effectivect,ccnt);
   CMinBLEIC::MinBLEICSetInnerCond(s.m_bs,0,0,m_xtol);
   CMinBLEIC::MinBLEICSetOuterCond(s.m_bs,m_xtol,1.0E-5);
   CMinBLEIC::MinBLEICSetPrecDiag(s.m_bs,s.m_h);
//--- solve problem
   while(CMinBLEIC::MinBLEICIteration(s.m_bs))
     {
      //--- check
      if(!CAp::Assert(s.m_bs.m_needfg,__FUNCTION__+": internal error"))
         return;
      //--- check
      if(s.m_bs.m_needfg)
        {
         //--- Calculate regularization term
         s.m_bs.m_f=0.0;
         vv=s.m_regterm;
         for(i=0;i<=n-1;i++)
           {
            for(j=0;j<=n-1;j++)
              {
               s.m_bs.m_f=s.m_bs.m_f+vv*CMath::Sqr(s.m_bs.m_x[i*n+j]-s.m_priorp[i][j]);
               s.m_bs.m_g[i*n+j]=2*vv*(s.m_bs.m_x[i*n+j]-s.m_priorp[i][j]);
              }
           }
         //--- calculate prediction error/gradient for K-th pair
         for(k=0;k<=npairs-1;k++)
           {
            for(i=0;i<=n-1;i++)
              {
               i1_=(0)-(i*n);
               v=0.0;
               for(i_=i*n;i_<=i*n+n-1;i_++)
                  v+=s.m_bs.m_x[i_]*s.m_data[k][i_+i1_];
               vv=s.m_effectivew[i];
               s.m_bs.m_f=s.m_bs.m_f+CMath::Sqr(vv*(v-s.m_data[k][n+i]));
               for(j=0;j<=n-1;j++)
                 {
                  s.m_bs.m_g[i*n+j]=s.m_bs.m_g[i*n+j]+2*vv*vv*(v-s.m_data[k][n+i])*s.m_data[k][j];
                 }
              }
           }
         //--- continue
         continue;
        }
     }
//--- function call
   CMinBLEIC::MinBLEICResultsBuf(s.m_bs,s.m_tmpp,s.m_br);
   for(i=0;i<=n-1;i++)
     {
      for(j=0;j<=n-1;j++)
         s.m_p[i].Set(j,s.m_tmpp[i*n+j]);
     }
//--- change values
   s.m_repterminationtype=s.m_br.m_terminationtype;
   s.m_repinneriterationscount=s.m_br.m_inneriterationscount;
   s.m_repouteriterationscount=s.m_br.m_outeriterationscount;
   s.m_repnfev=s.m_br.m_nfev;
  }
//+------------------------------------------------------------------+
//| MCPD results                                                     |
//| INPUT PARAMETERS:                                                |
//|     State   -   algorithm state                                  |
//| OUTPUT PARAMETERS:                                               |
//|     P      -   array[N,N], transition matrix                     |
//|     Rep    -   optimization report. You should check Rep.        |
//|                TerminationType in order to distinguish successful|
//|                termination from unsuccessful one. Speaking short,|
//|                positive values denote success, negative ones are |
//|                failures. More information about fields of this   |
//|                structure  can befound in the comments on         |
//|                MCPDReport datatype.                              |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDResults(CMCPDState &s,CMatrixDouble &p,
                                    CMCPDReport &rep)
  {
//--- create variables
   int i=0;
   int j=0;
//--- allocation
   p.Resize(s.m_n,s.m_n);
//--- copy
   for(i=0;i<=s.m_n-1;i++)
     {
      for(j=0;j<=s.m_n-1;j++)
         p[i].Set(j,s.m_p[i][j]);
     }
//--- change values
   rep.m_terminationtype=s.m_repterminationtype;
   rep.m_inneriterationscount=s.m_repinneriterationscount;
   rep.m_outeriterationscount=s.m_repouteriterationscount;
   rep.m_nfev=s.m_repnfev;
  }
//+------------------------------------------------------------------+
//| Internal initialization function                                 |
//+------------------------------------------------------------------+
static void CMarkovCPD::MCPDInit(const int n,const int entrystate,
                                 const int exitstate,CMCPDState &s)
  {
//--- create variables
   int i=0;
   int j=0;
//--- check
   if(!CAp::Assert(n>=1,__FUNCTION__+": N<1"))
      return;
//--- initialization
   s.m_n=n;
//--- allocation
   ArrayResizeAL(s.m_states,n);
   for(i=0;i<=n-1;i++)
      s.m_states[i]=0;
//--- check
   if(entrystate>=0)
      s.m_states[entrystate]=1;
//--- check
   if(exitstate>=0)
      s.m_states[exitstate]=-1;
//--- initialization
   s.m_npairs=0;
   s.m_regterm=1.0E-8;
   s.m_ccnt=0;
//--- allocation
   s.m_p.Resize(n,n);
   s.m_ec.Resize(n,n);
   s.m_bndl.Resize(n,n);
   s.m_bndu.Resize(n,n);
   ArrayResizeAL(s.m_pw,n);
   s.m_priorp.Resize(n,n);
   ArrayResizeAL(s.m_tmpp,n*n);
   ArrayResizeAL(s.m_effectivew,n);
   ArrayResizeAL(s.m_effectivebndl,n*n);
   ArrayResizeAL(s.m_effectivebndu,n*n);
   ArrayResizeAL(s.m_h,n*n);
//--- change values
   for(i=0;i<=n-1;i++)
     {
      for(j=0;j<=n-1;j++)
        {
         s.m_p[i].Set(j,0.0);
         s.m_priorp[i].Set(j,0.0);
         s.m_bndl[i].Set(j,CInfOrNaN::NegativeInfinity());
         s.m_bndu[i].Set(j,CInfOrNaN::PositiveInfinity());
         s.m_ec[i].Set(j,CInfOrNaN::NaN());
        }
      s.m_pw[i]=0.0;
      s.m_priorp[i].Set(i,1.0);
     }
//--- allocation
   s.m_data.Resize(1,2*n);
   for(i=0;i<=2*n-1;i++)
      s.m_data[0].Set(i,0.0);
   for(i=0;i<=n*n-1;i++)
      s.m_tmpp[i]=0.0;
//--- function call
   CMinBLEIC::MinBLEICCreate(n*n,s.m_tmpp,s.m_bs);
  }
//+------------------------------------------------------------------+
//| Training report:                                                 |
//|     * NGrad     - number of gradient calculations                |
//|     * NHess     - number of Hessian calculations                 |
//|     * NCholesky - number of Cholesky decompositions              |
//+------------------------------------------------------------------+
class CMLPReport
  {
public:
   //--- variables
   int               m_ngrad;
   int               m_nhess;
   int               m_ncholesky;
   //--- constructor, destructor
                     CMLPReport(void);
                    ~CMLPReport(void);
   //--- copy
   void              Copy(CMLPReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMLPReport::CMLPReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMLPReport::~CMLPReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMLPReport::Copy(CMLPReport &obj)
  {
//--- copy variables
   m_ngrad=obj.m_ngrad;
   m_nhess=obj.m_nhess;
   m_ncholesky=obj.m_ncholesky;
  }
//+------------------------------------------------------------------+
//| Training report:                                                 |
//|     * NGrad     - number of gradient calculations                |
//|     * NHess     - number of Hessian calculations                 |
//|     * NCholesky - number of Cholesky decompositions              |
//+------------------------------------------------------------------+
class CMLPReportShell
  {
private:
   CMLPReport        m_innerobj;
public:
   //--- constructors, destructor
                     CMLPReportShell(void);
                     CMLPReportShell(CMLPReport &obj);
                    ~CMLPReportShell(void);
   //--- methods
   int               GetNGrad(void);
   void              SetNGrad(const int i);
   int               GetNHess(void);
   void              SetNHess(const int i);
   int               GetNCholesky(void);
   void              SetNCholesky(const int i);
   CMLPReport       *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMLPReportShell::CMLPReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMLPReportShell::CMLPReportShell(CMLPReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMLPReportShell::~CMLPReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable ngrad                          |
//+------------------------------------------------------------------+
int CMLPReportShell::GetNGrad(void)
  {
//--- return result
   return(m_innerobj.m_ngrad);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable ngrad                         |
//+------------------------------------------------------------------+
void CMLPReportShell::SetNGrad(const int i)
  {
//--- change value
   m_innerobj.m_ngrad=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable nhess                          |
//+------------------------------------------------------------------+
int CMLPReportShell::GetNHess(void)
  {
//--- return result
   return(m_innerobj.m_nhess);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable nhess                         |
//+------------------------------------------------------------------+
void CMLPReportShell::SetNHess(const int i)
  {
//--- change value
   m_innerobj.m_nhess=i;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable ncholesky                      |
//+------------------------------------------------------------------+
int CMLPReportShell::GetNCholesky(void)
  {
//--- return result
   return(m_innerobj.m_ncholesky);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable ncholesky                     |
//+------------------------------------------------------------------+
void CMLPReportShell::SetNCholesky(const int i)
  {
//--- change value
   m_innerobj.m_ncholesky=i;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMLPReport *CMLPReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Cross-validation estimates of generalization error               |
//+------------------------------------------------------------------+
class CMLPCVReport
  {
public:
   //--- variables
   double            m_relclserror;
   double            m_avgce;
   double            m_rmserror;
   double            m_avgerror;
   double            m_avgrelerror;
   //--- constructor, destructor
                     CMLPCVReport(void);
                    ~CMLPCVReport(void);
   //--- copy
   void              Copy(CMLPCVReport &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMLPCVReport::CMLPCVReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMLPCVReport::~CMLPCVReport(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMLPCVReport::Copy(CMLPCVReport &obj)
  {
//--- copy variables
   m_relclserror=obj.m_relclserror;
   m_avgce=obj.m_avgce;
   m_rmserror=obj.m_rmserror;
   m_avgerror=obj.m_avgerror;
   m_avgrelerror=obj.m_avgrelerror;
  }
//+------------------------------------------------------------------+
//| Cross-validation estimates of generalization error               |
//+------------------------------------------------------------------+
class CMLPCVReportShell
  {
private:
   CMLPCVReport      m_innerobj;
public:
   //--- constructors, destructor
                     CMLPCVReportShell(void);
                     CMLPCVReportShell(CMLPCVReport &obj);
                    ~CMLPCVReportShell(void);
   //--- methods
   double            GetRelClsError(void);
   void              SetRelClsError(const double d);
   double            GetAvgCE(void);
   void              SetAvgCE(const double d);
   double            GetRMSError(void);
   void              SetRMSError(const double d);
   double            GetAvgError(void);
   void              SetAvgError(const double d);
   double            GetAvgRelError(void);
   void              SetAvgRelError(const double d);
   CMLPCVReport     *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMLPCVReportShell::CMLPCVReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMLPCVReportShell::CMLPCVReportShell(CMLPCVReport &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMLPCVReportShell::~CMLPCVReportShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Returns the value of the variable relclserror                    |
//+------------------------------------------------------------------+
double CMLPCVReportShell::GetRelClsError(void)
  {
//--- return result
   return(m_innerobj.m_relclserror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable relclserror                   |
//+------------------------------------------------------------------+
void CMLPCVReportShell::SetRelClsError(const double d)
  {
//--- change value
   m_innerobj.m_relclserror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgce                          |
//+------------------------------------------------------------------+
double CMLPCVReportShell::GetAvgCE(void)
  {
//--- return result
   return(m_innerobj.m_avgce);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgce                         |
//+------------------------------------------------------------------+
void CMLPCVReportShell::SetAvgCE(const double d)
  {
//--- change value
   m_innerobj.m_avgce=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable rmserror                       |
//+------------------------------------------------------------------+
double CMLPCVReportShell::GetRMSError(void)
  {
//--- return result
   return(m_innerobj.m_rmserror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable rmserror                      |
//+------------------------------------------------------------------+
void CMLPCVReportShell::SetRMSError(const double d)
  {
//--- change value
   m_innerobj.m_rmserror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgerror                       |
//+------------------------------------------------------------------+
double CMLPCVReportShell::GetAvgError(void)
  {
//--- return result
   return(m_innerobj.m_avgerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgerror                      |
//+------------------------------------------------------------------+
void CMLPCVReportShell::SetAvgError(const double d)
  {
//--- change value
   m_innerobj.m_avgerror=d;
  }
//+------------------------------------------------------------------+
//| Returns the value of the variable avgrelerror                    |
//+------------------------------------------------------------------+
double CMLPCVReportShell::GetAvgRelError(void)
  {
//--- return result
   return(m_innerobj.m_avgrelerror);
  }
//+------------------------------------------------------------------+
//| Changing the value of the variable avgrelerror                   |
//+------------------------------------------------------------------+
void CMLPCVReportShell::SetAvgRelError(const double d)
  {
//--- change value
   m_innerobj.m_avgrelerror=d;
  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMLPCVReport *CMLPCVReportShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Training neural networks                                         |
//+------------------------------------------------------------------+
class CMLPTrain
  {
private:
   //--- private methods
   static void       MLPKFoldCVGeneral(CMultilayerPerceptron &n,CMatrixDouble &xy,const int npoints,const double decay,const int restarts,const int foldscount,const bool lmalgorithm,const double wstep,const int maxits,int &info,CMLPReport &rep,CMLPCVReport &cvrep);
   static void       MLPKFoldSplit(CMatrixDouble &xy,const int npoints,const int nclasses,const int foldscount,const bool stratifiedsplits,int &folds[]);
public:
   //--- constant
   static const double m_mindecay;
   //--- constructor, destructor
                     CMLPTrain(void);
                    ~CMLPTrain(void);
   //--- public methods
   static void       MLPTrainLM(CMultilayerPerceptron &network,CMatrixDouble &xy,const int npoints,double decay,const int restarts,int &info,CMLPReport &rep);
   static void       MLPTrainLBFGS(CMultilayerPerceptron &network,CMatrixDouble &xy,const int npoints,double decay,const int restarts,const double wstep,int maxits,int &info,CMLPReport &rep);
   static void       MLPTrainES(CMultilayerPerceptron &network,CMatrixDouble &trnxy,const int trnsize,CMatrixDouble &valxy,const int valsize,const double decay,const int restarts,int &info,CMLPReport &rep);
   static void       MLPKFoldCVLBFGS(CMultilayerPerceptron &network,CMatrixDouble &xy,const int npoints,const double decay,const int restarts,const double wstep,const int maxits,const int foldscount,int &info,CMLPReport &rep,CMLPCVReport &cvrep);
   static void       MLPKFoldCVLM(CMultilayerPerceptron &network,CMatrixDouble &xy,const int npoints,const double decay,const int restarts,int foldscount,int &info,CMLPReport &rep,CMLPCVReport &cvrep);
  };
//+------------------------------------------------------------------+
//| Initialize constant                                              |
//+------------------------------------------------------------------+ 
const double CMLPTrain::m_mindecay=0.001;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMLPTrain::CMLPTrain(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMLPTrain::~CMLPTrain(void)
  {

  }
//+------------------------------------------------------------------+
//| Neural network training using modified Levenberg-Marquardt with  |
//| exact Hessian calculation and regularization. Subroutine trains  |
//| neural network with restarts from random positions. Algorithm is |
//| well suited for small                                            |
//| and medium scale problems (hundreds of weights).                 |
//| INPUT PARAMETERS:                                                |
//|     Network     -   neural network with initialized geometry     |
//|     XY          -   training set                                 |
//|     NPoints     -   training set size                            |
//|     Decay       -   weight decay constant, >=0.001               |
//|                     Decay term 'Decay*||Weights||^2' is added to |
//|                     error function.                              |
//|                     If you don't know what Decay to choose, use  |
//|                     0.001.                                       |
//|     Restarts    -   number of restarts from random position, >0. |
//|                     If you don't know what Restarts to choose,   |
//|                     use 2.                                       |
//| OUTPUT PARAMETERS:                                               |
//|     Network     -   trained neural network.                      |
//|     Info        -   return code:                                 |
//|                     * -9, if internal matrix inverse subroutine  |
//|                           failed                                 |
//|                     * -2, if there is a point with class number  |
//|                           outside of [0..NOut-1].                |
//|                     * -1, if wrong parameters specified          |
//|                           (NPoints<0, Restarts<1).               |
//|                     *  2, if task has been solved.               |
//|     Rep         -   training report                              |
//+------------------------------------------------------------------+
static void CMLPTrain::MLPTrainLM(CMultilayerPerceptron &network,CMatrixDouble &xy,
                                  const int npoints,double decay,const int restarts,
                                  int &info,CMLPReport &rep)
  {
//--- create variables
   int    nin=0;
   int    nout=0;
   int    wcount=0;
   double lmm_ftol=0;
   double lmsteptol=0;
   int    i=0;
   int    k=0;
   double v=0;
   double e=0;
   double enew=0;
   double xnorm2=0;
   double stepnorm=0;
   bool   spd;
   double nu=0;
   double lambdav=0;
   double lambdaup=0;
   double lambdadown=0;
   int    pass=0;
   double ebest=0;
   int    invinfo=0;
   int    solverinfo=0;
   int    i_=0;
//--- creating arrays
   double g[];
   double d[];
   double x[];
   double y[];
   double wbase[];
   double wdir[];
   double wt[];
   double wx[];
   double wbest[];
//--- create matrix
   CMatrixDouble h;
   CMatrixDouble hmod;
   CMatrixDouble z;
//--- objects of classes
   CMinLBFGSReport    internalrep;
   CMinLBFGSState     state;
   CMatInvReport      invrep;
   CDenseSolverReport solverrep;
//--- initialization
   info=0;
//--- function call
   CMLPBase::MLPProperties(network,nin,nout,wcount);
//--- initialization
   lambdaup=10;
   lambdadown=0.3;
   lmm_ftol=0.001;
   lmsteptol=0.001;
//--- Test for inputs
   if(npoints<=0 || restarts<1)
     {
      info=-1;
      return;
     }
//--- check
   if(CMLPBase::MLPIsSoftMax(network))
     {
      for(i=0;i<=npoints-1;i++)
        {
         //--- check
         if((int)MathRound(xy[i][nin])<0 || (int)MathRound(xy[i][nin])>=nout)
           {
            info=-2;
            return;
           }
        }
     }
//--- change values
   decay=MathMax(decay,m_mindecay);
   info=2;
//--- Initialize data
   rep.m_ngrad=0;
   rep.m_nhess=0;
   rep.m_ncholesky=0;
//--- General case.
//--- Prepare task and network. Allocate space.
   CMLPBase::MLPInitPreprocessor(network,xy,npoints);
//--- allocation
   ArrayResizeAL(g,wcount);
   h.Resize(wcount,wcount);
   hmod.Resize(wcount,wcount);
   ArrayResizeAL(wbase,wcount);
   ArrayResizeAL(wdir,wcount);
   ArrayResizeAL(wbest,wcount);
   ArrayResizeAL(wt,wcount);
   ArrayResizeAL(wx,wcount);
//--- initialization
   ebest=CMath::m_maxrealnumber;
//--- Multiple passes
   for(pass=1;pass<=restarts;pass++)
     {
      //--- Initialize weights
      CMLPBase::MLPRandomize(network);
      //--- First stage of the hybrid algorithm: LBFGS
      for(i_=0;i_<=wcount-1;i_++)
         wbase[i_]=network.m_weights[i_];
      //--- function calls
      CMinLBFGS::MinLBFGSCreate(wcount,(int)(MathMin(wcount,5)),wbase,state);
      CMinLBFGS::MinLBFGSSetCond(state,0,0,0,(int)(MathMax(25,wcount)));
      while(CMinLBFGS::MinLBFGSIteration(state))
        {
         //--- gradient
         for(i_=0;i_<=wcount-1;i_++)
            network.m_weights[i_]=state.m_x[i_];
         //--- function call
         CMLPBase::MLPGradBatch(network,xy,npoints,state.m_f,state.m_g);
         //--- weight decay
         v=0.0;
         for(i_=0;i_<=wcount-1;i_++)
            v+=network.m_weights[i_]*network.m_weights[i_];
         state.m_f=state.m_f+0.5*decay*v;
         for(i_=0;i_<=wcount-1;i_++)
            state.m_g[i_]=state.m_g[i_]+decay*network.m_weights[i_];
         //--- next iteration
         rep.m_ngrad=rep.m_ngrad+1;
        }
      //--- function call
      CMinLBFGS::MinLBFGSResults(state,wbase,internalrep);
      for(i_=0;i_<=wcount-1;i_++)
         network.m_weights[i_]=wbase[i_];
      //--- Second stage of the hybrid algorithm: LM
      //--- Initialize H with identity matrix,
      //--- G with gradient,
      //--- E with regularized error.
      CMLPBase::MLPHessianBatch(network,xy,npoints,e,g,h);
      v=0.0;
      for(i_=0;i_<=wcount-1;i_++)
         v+=network.m_weights[i_]*network.m_weights[i_];
      //--- change values
      e=e+0.5*decay*v;
      for(i_=0;i_<=wcount-1;i_++)
         g[i_]=g[i_]+decay*network.m_weights[i_];
      for(k=0;k<=wcount-1;k++)
         h[k].Set(k,h[k][k]+decay);
      //--- change values
      rep.m_nhess=rep.m_nhess+1;
      lambdav=0.001;
      nu=2;
      //--- cycle
      while(true)
        {
         //--- 1. HMod=H+lambda*I
         //--- 2. Try to solve (H+Lambda*I)*dx=-g.
         //---    Increase lambda if left part is not positive definite.
         for(i=0;i<=wcount-1;i++)
           {
            for(i_=0;i_<=wcount-1;i_++)
               hmod[i].Set(i_,h[i][i_]);
            hmod[i].Set(i,hmod[i][i]+lambdav);
           }
         //--- function call
         spd=CTrFac::SPDMatrixCholesky(hmod,wcount,true);
         rep.m_ncholesky=rep.m_ncholesky+1;
         //--- check
         if(!spd)
           {
            lambdav=lambdav*lambdaup*nu;
            nu=nu*2;
            continue;
           }
         //--- function call
         CDenseSolver::SPDMatrixCholeskySolve(hmod,wcount,true,g,solverinfo,solverrep,wdir);
         //--- check
         if(solverinfo<0)
           {
            lambdav=lambdav*lambdaup*nu;
            nu=nu*2;
            continue;
           }
         for(i_=0;i_<=wcount-1;i_++)
            wdir[i_]=-1*wdir[i_];
         //--- Lambda found.
         //--- 1. Save old w in WBase
         //--- 1. Test some stopping criterions
         //--- 2. If error(w+wdir)>error(w),increase lambda
         for(i_=0;i_<=wcount-1;i_++)
            network.m_weights[i_]=network.m_weights[i_]+wdir[i_];
         xnorm2=0.0;
         for(i_=0;i_<=wcount-1;i_++)
            xnorm2+=network.m_weights[i_]*network.m_weights[i_];
         //--- change value
         stepnorm=0.0;
         for(i_=0;i_<=wcount-1;i_++)
            stepnorm+=wdir[i_]*wdir[i_];
         stepnorm=MathSqrt(stepnorm);
         //--- function call
         enew=CMLPBase::MLPError(network,xy,npoints)+0.5*decay*xnorm2;
         //--- check
         if(stepnorm<lmsteptol*(1+MathSqrt(xnorm2)))
            break;
         //--- check
         if(enew>e)
           {
            lambdav=lambdav*lambdaup*nu;
            nu=nu*2;
            continue;
           }
         //--- Optimize using inv(cholesky(H)) as preconditioner
         CMatInv::RMatrixTrInverse(hmod,wcount,true,false,invinfo,invrep);
         //--- check
         if(invinfo<=0)
           {
            //--- if matrix can't be inverted then exit with errors
            //--- TODO: make WCount steps in direction suggested by HMod
            info=-9;
            return;
           }
         //--- calculation
         for(i_=0;i_<=wcount-1;i_++)
            wbase[i_]=network.m_weights[i_];
         for(i=0;i<=wcount-1;i++)
            wt[i]=0;
         //--- function calls
         CMinLBFGS::MinLBFGSCreateX(wcount,wcount,wt,1,0.0,state);
         CMinLBFGS::MinLBFGSSetCond(state,0,0,0,5);
         while(CMinLBFGS::MinLBFGSIteration(state))
           {
            //--- gradient
            for(i=0;i<=wcount-1;i++)
              {
               v=0.0;
               for(i_=i;i_<=wcount-1;i_++)
                  v+=state.m_x[i_]*hmod[i][i_];
               network.m_weights[i]=wbase[i]+v;
              }
            //--- function call
            CMLPBase::MLPGradBatch(network,xy,npoints,state.m_f,g);
            for(i=0;i<=wcount-1;i++)
               state.m_g[i]=0;
            for(i=0;i<=wcount-1;i++)
              {
               v=g[i];
               for(i_=i;i_<=wcount-1;i_++)
                  state.m_g[i_]=state.m_g[i_]+v*hmod[i][i_];
              }
            //--- weight decay
            //--- grad(x'*x)=A'*(x0+A*t)
            v=0.0;
            for(i_=0;i_<=wcount-1;i_++)
               v+=network.m_weights[i_]*network.m_weights[i_];
            state.m_f=state.m_f+0.5*decay*v;
            for(i=0;i<=wcount-1;i++)
              {
               v=decay*network.m_weights[i];
               for(i_=i;i_<=wcount-1;i_++)
                  state.m_g[i_]=state.m_g[i_]+v*hmod[i][i_];
              }
            //--- next iteration
            rep.m_ngrad=rep.m_ngrad+1;
           }
         //--- function call
         CMinLBFGS::MinLBFGSResults(state,wt,internalrep);
         //--- Accept new position.
         //--- Calculate Hessian
         for(i=0;i<=wcount-1;i++)
           {
            v=0.0;
            for(i_=i;i_<=wcount-1;i_++)
               v+=wt[i_]*hmod[i][i_];
            network.m_weights[i]=wbase[i]+v;
           }
         //--- function call
         CMLPBase::MLPHessianBatch(network,xy,npoints,e,g,h);
         v=0.0;
         for(i_=0;i_<=wcount-1;i_++)
            v+=network.m_weights[i_]*network.m_weights[i_];
         //--- change value
         e=e+0.5*decay*v;
         for(i_=0;i_<=wcount-1;i_++)
            g[i_]=g[i_]+decay*network.m_weights[i_];
         for(k=0;k<=wcount-1;k++)
            h[k].Set(k,h[k][k]+decay);
         rep.m_nhess=rep.m_nhess+1;
         //--- Update lambda
         lambdav=lambdav*lambdadown;
         nu=2;
        }
      //--- update WBest
      v=0.0;
      for(i_=0;i_<=wcount-1;i_++)
         v+=network.m_weights[i_]*network.m_weights[i_];
      //--- change value
      e=0.5*decay*v+CMLPBase::MLPError(network,xy,npoints);
      //--- check
      if(e<ebest)
        {
         ebest=e;
         for(i_=0;i_<=wcount-1;i_++)
            wbest[i_]=network.m_weights[i_];
        }
     }
//--- copy WBest to output
   for(i_=0;i_<=wcount-1;i_++)
      network.m_weights[i_]=wbest[i_];
  }
//+------------------------------------------------------------------+
//| Neural network training using L-BFGS algorithm with              |
//| regularization. Subroutine trains neural network with restarts   |
//| from random positions. Algorithm is well suited for problems of  |
//| any dimensionality (memory requirements and step complexity are  |
//| linear by weights number).                                       |
//| INPUT PARAMETERS:                                                |
//|     Network    -   neural network with initialized geometry      |
//|     XY         -   training set                                  |
//|     NPoints    -   training set size                             |
//|     Decay      -   weight decay constant, >=0.001                |
//|                    Decay term 'Decay*||Weights||^2' is added to  |
//|                    error function.                               |
//|                    If you don't know what Decay to choose, use   |
//|                    0.001.                                        |
//|     Restarts   -   number of restarts from random position, >0.  |
//|                    If you don't know what Restarts to choose,    |
//|                    use 2.                                        |
//|     WStep      -   stopping criterion. Algorithm stops if step   |
//|                    size is less than WStep. Recommended          |
//|                    value - 0.01. Zero step size means stopping   |
//|                    after MaxIts iterations.                      |
//|     MaxIts     -   stopping criterion. Algorithm stops after     |
//|                    MaxIts iterations (NOT gradient calculations).|
//|                    Zero MaxIts means stopping when step is       |
//|                    sufficiently small.                           |
//| OUTPUT PARAMETERS:                                               |
//|     Network     -   trained neural network.                      |
//|     Info        -   return code:                                 |
//|                     * -8, if both WStep=0 and MaxIts=0           |
//|                     * -2, if there is a point with class number  |
//|                           outside of [0..NOut-1].                |
//|                     * -1, if wrong parameters specified          |
//|                           (NPoints<0, Restarts<1).               |
//|                     *  2, if task has been solved.               |
//|     Rep         -   training report                              |
//+------------------------------------------------------------------+
static void CMLPTrain::MLPTrainLBFGS(CMultilayerPerceptron &network,
                                     CMatrixDouble &xy,const int npoints,
                                     double decay,const int restarts,
                                     const double wstep,int maxits,
                                     int &info,CMLPReport &rep)
  {
//--- create variables
   int    i=0;
   int    pass=0;
   int    nin=0;
   int    nout=0;
   int    wcount=0;
   double e=0;
   double v=0;
   double ebest=0;
//--- creating arrays
   double w[];
   double wbest[];
//--- create objects of classes
   CMinLBFGSReport internalrep;
   CMinLBFGSState  state;
   int i_=0;
//--- initialization
   info=0;
//--- Test inputs,parse flags,read network geometry
   if(wstep==0.0 && maxits==0)
     {
      info=-8;
      return;
     }
//--- check
   if(((npoints<=0 || restarts<1) || wstep<0.0) || maxits<0)
     {
      info=-1;
      return;
     }
//--- function call
   CMLPBase::MLPProperties(network,nin,nout,wcount);
//--- check
   if(CMLPBase::MLPIsSoftMax(network))
     {
      for(i=0;i<=npoints-1;i++)
        {
         //--- check
         if((int)MathRound(xy[i][nin])<0 || (int)MathRound(xy[i][nin])>=nout)
           {
            info=-2;
            return;
           }
        }
     }
//--- change values
   decay=MathMax(decay,m_mindecay);
   info=2;
//--- Prepare
   CMLPBase::MLPInitPreprocessor(network,xy,npoints);
//--- allocation
   ArrayResizeAL(w,wcount);
   ArrayResizeAL(wbest,wcount);
//--- initialization
   ebest=CMath::m_maxrealnumber;
//--- Multiple starts
   rep.m_ncholesky=0;
   rep.m_nhess=0;
   rep.m_ngrad=0;
   for(pass=1;pass<=restarts;pass++)
     {
      //--- Process
      CMLPBase::MLPRandomize(network);
      for(i_=0;i_<=wcount-1;i_++)
         w[i_]=network.m_weights[i_];
      //--- function calls
      CMinLBFGS::MinLBFGSCreate(wcount,(int)(MathMin(wcount,10)),w,state);
      CMinLBFGS::MinLBFGSSetCond(state,0.0,0.0,wstep,maxits);
      while(CMinLBFGS::MinLBFGSIteration(state))
        {
         for(i_=0;i_<=wcount-1;i_++)
            network.m_weights[i_]=state.m_x[i_];
         //--- function call
         CMLPBase::MLPGradNBatch(network,xy,npoints,state.m_f,state.m_g);
         v=0.0;
         for(i_=0;i_<=wcount-1;i_++)
            v+=network.m_weights[i_]*network.m_weights[i_];
         state.m_f=state.m_f+0.5*decay*v;
         for(i_=0;i_<=wcount-1;i_++)
            state.m_g[i_]=state.m_g[i_]+decay*network.m_weights[i_];
         rep.m_ngrad=rep.m_ngrad+1;
        }
      //--- function call
      CMinLBFGS::MinLBFGSResults(state,w,internalrep);
      for(i_=0;i_<=wcount-1;i_++)
         network.m_weights[i_]=w[i_];
      //--- Compare with best
      v=0.0;
      for(i_=0;i_<=wcount-1;i_++)
         v+=network.m_weights[i_]*network.m_weights[i_];
      //--- change value
      e=CMLPBase::MLPErrorN(network,xy,npoints)+0.5*decay*v;
      //--- check
      if(e<ebest)
        {
         for(i_=0;i_<=wcount-1;i_++)
            wbest[i_]=network.m_weights[i_];
         ebest=e;
        }
     }
//--- The best network
   for(i_=0;i_<=wcount-1;i_++)
      network.m_weights[i_]=wbest[i_];
  }
//+------------------------------------------------------------------+
//| Neural network training using early stopping (base algorithm -   |
//| L-BFGS with regularization).                                     |
//| INPUT PARAMETERS:                                                |
//|     Network     -   neural network with initialized geometry     |
//|     TrnXY       -   training set                                 |
//|     TrnSize     -   training set size                            |
//|     ValXY       -   validation set                               |
//|     ValSize     -   validation set size                          |
//|     Decay       -   weight decay constant, >=0.001               |
//|                     Decay term 'Decay*||Weights||^2' is added to |
//|                     error function.                              |
//|                     If you don't know what Decay to choose, use  |
//|                     0.001.                                       |
//|     Restarts    -   number of restarts from random position, >0. |
//|                     If you don't know what Restarts to choose,   |
//|                     use 2.                                       |
//| OUTPUT PARAMETERS:                                               |
//|     Network     -   trained neural network.                      |
//|     Info        -   return code:                                 |
//|                     * -2, if there is a point with class number  |
//|                           outside of [0..NOut-1].                |
//|                     * -1, if wrong parameters specified          |
//|                           (NPoints<0, Restarts<1, ...).          |
//|                     *  2, task has been solved, stopping         |
//|                           criterion met - sufficiently small     |
//|                           step size. Not expected (we use EARLY  |
//|                           stopping) but possible and not an error|
//|                     *  6, task has been solved, stopping         |
//|                           criterion  met - increasing of         |
//|                           validation set error.                  |
//|     Rep         -   training report                              |
//| NOTE:                                                            |
//| Algorithm stops if validation set error increases for a long     |
//| enough or step size is small enought (there are task where       |
//| validation set may decrease for eternity). In any case solution  |
//| returned corresponds to the minimum of validation set error.     |
//+------------------------------------------------------------------+
static void CMLPTrain::MLPTrainES(CMultilayerPerceptron &network,
                                  CMatrixDouble &trnxy,const int trnsize,
                                  CMatrixDouble &valxy,const int valsize,
                                  const double decay,const int restarts,
                                  int &info,CMLPReport &rep)
  {
//--- create variables
   int    i=0;
   int    pass=0;
   int    nin=0;
   int    nout=0;
   int    wcount=0;
   double e=0;
   double v=0;
   double ebest=0;
   int    itbest=0;
   double wstep=0;
   int    i_=0;
//--- creating arrays
   double w[];
   double wbest[];
   double wfinal[];
   double efinal=0;
//--- objects of classes
   CMinLBFGSReport internalrep;
   CMinLBFGSState state;
//--- initialization
   info=0;
   wstep=0.001;
//--- Test inputs,parse flags,read network geometry
   if(((trnsize<=0 || valsize<=0) || restarts<1) || decay<0.0)
     {
      info=-1;
      return;
     }
//--- function call
   CMLPBase::MLPProperties(network,nin,nout,wcount);
//--- check
   if(CMLPBase::MLPIsSoftMax(network))
     {
      for(i=0;i<=trnsize-1;i++)
        {
         //--- check
         if((int)MathRound(trnxy[i][nin])<0 || (int)MathRound(trnxy[i][nin])>=nout)
           {
            info=-2;
            return;
           }
        }
      for(i=0;i<=valsize-1;i++)
        {
         //--- check
         if((int)MathRound(valxy[i][nin])<0 || (int)MathRound(valxy[i][nin])>=nout)
           {
            info=-2;
            return;
           }
        }
     }
//--- change value
   info=2;
//--- Prepare
   CMLPBase::MLPInitPreprocessor(network,trnxy,trnsize);
//--- allocation
   ArrayResizeAL(w,wcount);
   ArrayResizeAL(wbest,wcount);
   ArrayResizeAL(wfinal,wcount);
//--- initialization
   efinal=CMath::m_maxrealnumber;
   for(i=0;i<=wcount-1;i++)
      wfinal[i]=0;
//--- Multiple starts
   rep.m_ncholesky=0;
   rep.m_nhess=0;
   rep.m_ngrad=0;
//--- calculation
   for(pass=1;pass<=restarts;pass++)
     {
      //--- Process
      CMLPBase::MLPRandomize(network);
      //--- change values
      ebest=CMLPBase::MLPError(network,valxy,valsize);
      for(i_=0;i_<=wcount-1;i_++)
         wbest[i_]=network.m_weights[i_];
      //--- change values
      itbest=0;
      for(i_=0;i_<=wcount-1;i_++)
         w[i_]=network.m_weights[i_];
      //--- function calls
      CMinLBFGS::MinLBFGSCreate(wcount,(int)(MathMin(wcount,10)),w,state);
      CMinLBFGS::MinLBFGSSetCond(state,0.0,0.0,wstep,0);
      CMinLBFGS::MinLBFGSSetXRep(state,true);
      while(CMinLBFGS::MinLBFGSIteration(state))
        {
         //--- Calculate gradient
         for(i_=0;i_<=wcount-1;i_++)
            network.m_weights[i_]=state.m_x[i_];
         //--- function call
         CMLPBase::MLPGradNBatch(network,trnxy,trnsize,state.m_f,state.m_g);
         v=0.0;
         for(i_=0;i_<=wcount-1;i_++)
            v+=network.m_weights[i_]*network.m_weights[i_];
         state.m_f=state.m_f+0.5*decay*v;
         for(i_=0;i_<=wcount-1;i_++)
            state.m_g[i_]=state.m_g[i_]+decay*network.m_weights[i_];
         rep.m_ngrad=rep.m_ngrad+1;
         //--- Validation set
         if(state.m_xupdated)
           {
            for(i_=0;i_<=wcount-1;i_++)
               network.m_weights[i_]=w[i_];
            //--- function call
            e=CMLPBase::MLPError(network,valxy,valsize);
            //--- check
            if(e<ebest)
              {
               ebest=e;
               for(i_=0;i_<=wcount-1;i_++)
                  wbest[i_]=network.m_weights[i_];
               itbest=internalrep.m_iterationscount;
              }
            //--- check
            if(internalrep.m_iterationscount>30 && (double)(internalrep.m_iterationscount)>(double)(1.5*itbest))
              {
               info=6;
               break;
              }
           }
        }
      //--- function call
      CMinLBFGS::MinLBFGSResults(state,w,internalrep);
      //--- Compare with final answer
      if(ebest<efinal)
        {
         for(i_=0;i_<=wcount-1;i_++)
            wfinal[i_]=wbest[i_];
         efinal=ebest;
        }
     }
//--- The best network
   for(i_=0;i_<=wcount-1;i_++)
      network.m_weights[i_]=wfinal[i_];
  }
//+------------------------------------------------------------------+
//| Cross-validation estimate of generalization error.               |
//| Base algorithm - L-BFGS.                                         |
//| INPUT PARAMETERS:                                                |
//|     Network     -   neural network with initialized geometry.    |
//|                     Network is not changed during                |
//|                     cross-validation - it is used only as a      |
//|                     representative of its architecture.          |
//|     XY          -   training set.                                |
//|     SSize       -   training set size                            |
//|     Decay       -   weight  decay, same as in MLPTrainLBFGS      |
//|     Restarts    -   number of restarts, >0.                      |
//|                     restarts are counted for each partition      |
//|                     separately, so total number of restarts will |
//|                     be Restarts*FoldsCount.                      |
//|     WStep       -   stopping criterion, same as in MLPTrainLBFGS |
//|     MaxIts      -   stopping criterion, same as in MLPTrainLBFGS |
//|     FoldsCount  -   number of folds in k-fold cross-validation,  |
//|                     2<=FoldsCount<=SSize.                        |
//|                     recommended value: 10.                       |
//| OUTPUT PARAMETERS:                                               |
//|     Info        -   return code, same as in MLPTrainLBFGS        |
//|     Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS  |
//|     CVRep       -   generalization error estimates               |
//+------------------------------------------------------------------+
static void CMLPTrain::MLPKFoldCVLBFGS(CMultilayerPerceptron &network,
                                       CMatrixDouble &xy,const int npoints,
                                       const double decay,const int restarts,
                                       const double wstep,const int maxits,
                                       const int foldscount,int &info,
                                       CMLPReport &rep,CMLPCVReport &cvrep)
  {
//--- initialization
   info=0;
//--- function call
   MLPKFoldCVGeneral(network,xy,npoints,decay,restarts,foldscount,false,wstep,maxits,info,rep,cvrep);
  }
//+------------------------------------------------------------------+
//| Cross-validation estimate of generalization error.               |
//| Base algorithm - Levenberg-Marquardt.                            |
//| INPUT PARAMETERS:                                                |
//|     Network     -   neural network with initialized geometry.    |
//|                     Network is not changed during                |
//|                     cross-validation - it is used only as a      |
//|                     representative of its architecture.          |
//|     XY          -   training set.                                |
//|     SSize       -   training set size                            |
//|     Decay       -   weight  decay, same as in MLPTrainLBFGS      |
//|     Restarts    -   number of restarts, >0.                      |
//|                     restarts are counted for each partition      |
//|                     separately, so total number of restarts will |
//|                     be Restarts*FoldsCount.                      |
//|     FoldsCount  -   number of folds in k-fold cross-validation,  |
//|                     2<=FoldsCount<=SSize.                        |
//|                     recommended value: 10.                       |
//| OUTPUT PARAMETERS:                                               |
//|     Info        -   return code, same as in MLPTrainLBFGS        |
//|     Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS  |
//|     CVRep       -   generalization error estimates               |
//+------------------------------------------------------------------+
static void CMLPTrain::MLPKFoldCVLM(CMultilayerPerceptron &network,
                                    CMatrixDouble &xy,const int npoints,
                                    const double decay,const int restarts,
                                    int foldscount,int &info,CMLPReport &rep,
                                    CMLPCVReport &cvrep)
  {
//--- initialization
   info=0;
//--- function call
   MLPKFoldCVGeneral(network,xy,npoints,decay,restarts,foldscount,true,0.0,0,info,rep,cvrep);
  }
//+------------------------------------------------------------------+
//| Internal cross-validation subroutine                             |
//+------------------------------------------------------------------+
static void CMLPTrain::MLPKFoldCVGeneral(CMultilayerPerceptron &n,
                                         CMatrixDouble &xy,const int npoints,
                                         const double decay,const int restarts,
                                         const int foldscount,const bool lmalgorithm,
                                         const double wstep,const int maxits,
                                         int &info,CMLPReport &rep,
                                         CMLPCVReport &cvrep)
  {
//--- create variables
   int i=0;
   int fold=0;
   int j=0;
   int k=0;
   int nin=0;
   int nout=0;
   int rowlen=0;
   int wcount=0;
   int nclasses=0;
   int tssize=0;
   int cvssize=0;
   int relcnt=0;
   int i_=0;
//--- creating arrays
   int    folds[];
   double x[];
   double y[];
//--- create matrix
   CMatrixDouble cvset;
   CMatrixDouble testset;
//--- creating arrays
   CMultilayerPerceptron network;
   CMLPReport            internalrep;
//--- initialization
   info=0;
//--- Read network geometry,test parameters
   CMLPBase::MLPProperties(n,nin,nout,wcount);
//--- check
   if(CMLPBase::MLPIsSoftMax(n))
     {
      nclasses=nout;
      rowlen=nin+1;
     }
   else
     {
      nclasses=-nout;
      rowlen=nin+nout;
     }
//--- check
   if((npoints<=0 || foldscount<2) || foldscount>npoints)
     {
      info=-1;
      return;
     }
//--- function call
   CMLPBase::MLPCopy(n,network);
//--- K-fold out cross-validation.
//--- First,estimate generalization error
   testset.Resize(npoints,rowlen);
   cvset.Resize(npoints,rowlen);
   ArrayResizeAL(x,nin);
   ArrayResizeAL(y,nout);
//--- function call
   MLPKFoldSplit(xy,npoints,nclasses,foldscount,false,folds);
//--- change values
   cvrep.m_relclserror=0;
   cvrep.m_avgce=0;
   cvrep.m_rmserror=0;
   cvrep.m_avgerror=0;
   cvrep.m_avgrelerror=0;
   rep.m_ngrad=0;
   rep.m_nhess=0;
   rep.m_ncholesky=0;
   relcnt=0;
//--- calculation
   for(fold=0;fold<=foldscount-1;fold++)
     {
      //--- Separate set
      tssize=0;
      cvssize=0;
      for(i=0;i<=npoints-1;i++)
        {
         //--- check
         if(folds[i]==fold)
           {
            for(i_=0;i_<=rowlen-1;i_++)
               testset[tssize].Set(i_,xy[i][i_]);
            tssize=tssize+1;
           }
         else
           {
            for(i_=0;i_<=rowlen-1;i_++)
               cvset[cvssize].Set(i_,xy[i][i_]);
            cvssize=cvssize+1;
           }
        }
      //--- Train on CV training set
      if(lmalgorithm)
         MLPTrainLM(network,cvset,cvssize,decay,restarts,info,internalrep);
      else
         MLPTrainLBFGS(network,cvset,cvssize,decay,restarts,wstep,maxits,info,internalrep);
      //--- check
      if(info<0)
        {
         //--- change values
         cvrep.m_relclserror=0;
         cvrep.m_avgce=0;
         cvrep.m_rmserror=0;
         cvrep.m_avgerror=0;
         cvrep.m_avgrelerror=0;
         //--- exit the function
         return;
        }
      //--- change values
      rep.m_ngrad=rep.m_ngrad+internalrep.m_ngrad;
      rep.m_nhess=rep.m_nhess+internalrep.m_nhess;
      rep.m_ncholesky=rep.m_ncholesky+internalrep.m_ncholesky;
      //--- Estimate error using CV test set
      if(CMLPBase::MLPIsSoftMax(network))
        {
         //--- classification-only code
         cvrep.m_relclserror=cvrep.m_relclserror+CMLPBase::MLPClsError(network,testset,tssize);
         cvrep.m_avgce=cvrep.m_avgce+CMLPBase::MLPErrorN(network,testset,tssize);
        }
      //--- calculation
      for(i=0;i<=tssize-1;i++)
        {
         for(i_=0;i_<=nin-1;i_++)
            x[i_]=testset[i][i_];
         //--- function call
         CMLPBase::MLPProcess(network,x,y);
         //--- check
         if(CMLPBase::MLPIsSoftMax(network))
           {
            //--- Classification-specific code
            k=(int)MathRound(testset[i][nin]);
            for(j=0;j<=nout-1;j++)
              {
               //--- check
               if(j==k)
                 {
                  //--- change values
                  cvrep.m_rmserror=cvrep.m_rmserror+CMath::Sqr(y[j]-1);
                  cvrep.m_avgerror=cvrep.m_avgerror+MathAbs(y[j]-1);
                  cvrep.m_avgrelerror=cvrep.m_avgrelerror+MathAbs(y[j]-1);
                  relcnt=relcnt+1;
                 }
               else
                 {
                  //--- change values
                  cvrep.m_rmserror=cvrep.m_rmserror+CMath::Sqr(y[j]);
                  cvrep.m_avgerror=cvrep.m_avgerror+MathAbs(y[j]);
                 }
              }
           }
         else
           {
            //--- Regression-specific code
            for(j=0;j<=nout-1;j++)
              {
               cvrep.m_rmserror=cvrep.m_rmserror+CMath::Sqr(y[j]-testset[i][nin+j]);
               cvrep.m_avgerror=cvrep.m_avgerror+MathAbs(y[j]-testset[i][nin+j]);
               //--- check
               if(testset[i][nin+j]!=0.0)
                 {
                  cvrep.m_avgrelerror=cvrep.m_avgrelerror+MathAbs((y[j]-testset[i][nin+j])/testset[i][nin+j]);
                  relcnt=relcnt+1;
                 }
              }
           }
        }
     }
//--- check
   if(CMLPBase::MLPIsSoftMax(network))
     {
      cvrep.m_relclserror=cvrep.m_relclserror/npoints;
      cvrep.m_avgce=cvrep.m_avgce/(MathLog(2)*npoints);
     }
//--- change values
   cvrep.m_rmserror=MathSqrt(cvrep.m_rmserror/(npoints*nout));
   cvrep.m_avgerror=cvrep.m_avgerror/(npoints*nout);
   cvrep.m_avgrelerror=cvrep.m_avgrelerror/relcnt;
   info=1;
  }
//+------------------------------------------------------------------+
//| Subroutine prepares K-fold split of the training set.            |
//| NOTES:                                                           |
//|     "NClasses>0" means that we have classification task.         |
//|     "NClasses<0" means regression task with -NClasses real       |
//|     outputs.                                                     |
//+------------------------------------------------------------------+
static void CMLPTrain::MLPKFoldSplit(CMatrixDouble &xy,const int npoints,
                                     const int nclasses,const int foldscount,
                                     const bool stratifiedsplits,int &folds[])
  {
//--- create variables
   int i=0;
   int j=0;
   int k=0;
//--- test parameters
   if(!CAp::Assert(npoints>0,__FUNCTION__+": wrong NPoints!"))
      return;
//--- check
   if(!CAp::Assert(nclasses>1 || nclasses<0,__FUNCTION__+": wrong NClasses!"))
      return;
//--- check
   if(!CAp::Assert(foldscount>=2 && foldscount<=npoints,__FUNCTION__+" wrong FoldsCount!"))
      return;
//--- check
   if(!CAp::Assert(!stratifiedsplits,__FUNCTION__+": stratified splits are not supported!"))
      return;
//--- Folds
   ArrayResizeAL(folds,npoints);
   for(i=0;i<=npoints-1;i++)
      folds[i]=i*foldscount/npoints;
//--- calculation
   for(i=0;i<=npoints-2;i++)
     {
      j=i+CMath::RandomInteger(npoints-i);
      //--- check
      if(j!=i)
        {
         k=folds[i];
         folds[i]=folds[j];
         folds[j]=k;
        }
     }
  }
//+------------------------------------------------------------------+
//| Neural networks ensemble                                         |
//+------------------------------------------------------------------+
class CMLPEnsemble
  {
public:
   //--- variables
   int               m_ensemblesize;
   int               m_nin;
   int               m_nout;
   int               m_wcount;
   bool              m_issoftmax;
   bool              m_postprocessing;
   int               m_serializedlen;
   //--- arrays
   int               m_structinfo[];
   double            m_weights[];
   double            m_columnmeans[];
   double            m_columnsigmas[];
   double            m_serializedmlp[];
   double            m_tmpweights[];
   double            m_tmpmeans[];
   double            m_tmpsigmas[];
   double            m_neurons[];
   double            m_dfdnet[];
   double            m_y[];
   //--- constructor, destructor
                     CMLPEnsemble(void);
                    ~CMLPEnsemble(void);
   //--- copy
   void              Copy(CMLPEnsemble &obj);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMLPEnsemble::CMLPEnsemble(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMLPEnsemble::~CMLPEnsemble(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy                                                             |
//+------------------------------------------------------------------+
void CMLPEnsemble::Copy(CMLPEnsemble &obj)
  {
//--- copy variables
   m_ensemblesize=obj.m_ensemblesize;
   m_nin=obj.m_nin;
   m_nout=obj.m_nout;
   m_wcount=obj.m_wcount;
   m_issoftmax=obj.m_issoftmax;
   m_postprocessing=obj.m_postprocessing;
   m_serializedlen=obj.m_serializedlen;
//--- copy arrays
   ArrayCopy(m_structinfo,obj.m_structinfo);
   ArrayCopy(m_weights,obj.m_weights);
   ArrayCopy(m_columnmeans,obj.m_columnmeans);
   ArrayCopy(m_columnsigmas,obj.m_columnsigmas);
   ArrayCopy(m_serializedmlp,obj.m_serializedmlp);
   ArrayCopy(m_tmpweights,obj.m_tmpweights);
   ArrayCopy(m_tmpmeans,obj.m_tmpmeans);
   ArrayCopy(m_tmpsigmas,obj.m_tmpsigmas);
   ArrayCopy(m_neurons,obj.m_neurons);
   ArrayCopy(m_dfdnet,obj.m_dfdnet);
   ArrayCopy(m_y,obj.m_y);
  }
//+------------------------------------------------------------------+
//| Neural networks ensemble                                         |
//+------------------------------------------------------------------+
class CMLPEnsembleShell
  {
private:
   CMLPEnsemble      m_innerobj;
public:
   //--- constructors, destructor
                     CMLPEnsembleShell(void);
                     CMLPEnsembleShell(CMLPEnsemble &obj);
                    ~CMLPEnsembleShell(void);
   //--- method
   CMLPEnsemble     *GetInnerObj(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMLPEnsembleShell::CMLPEnsembleShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy constructor                                                 |
//+------------------------------------------------------------------+
CMLPEnsembleShell::CMLPEnsembleShell(CMLPEnsemble &obj)
  {
//--- copy
   m_innerobj.Copy(obj);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMLPEnsembleShell::~CMLPEnsembleShell(void)
  {

  }
//+------------------------------------------------------------------+
//| Return object of class                                           |
//+------------------------------------------------------------------+
CMLPEnsemble *CMLPEnsembleShell::GetInnerObj(void)
  {
//--- return result
   return(GetPointer(m_innerobj));
  }
//+------------------------------------------------------------------+
//| Neural networks ensemble                                         |
//+------------------------------------------------------------------+
class CMLPE
  {
private:
   //--- private methods
   static void       MLPEAllErrors(CMLPEnsemble &ensemble,CMatrixDouble &xy,const int npoints,double &relcls,double &avgce,double &rms,double &avg,double &avgrel);
   static void       MLPEBaggingInternal(CMLPEnsemble &ensemble,CMatrixDouble &xy,const int npoints,const double decay,const int restarts,const double wstep,const int maxits,const bool lmalgorithm,int &info,CMLPReport &rep,CMLPCVReport &ooberrors);
public:
   //--- class constants
   static const int  m_mlpntotaloffset;
   static const int  m_mlpevnum;
   //--- constructor, destructor
                     CMLPE(void);
                    ~CMLPE(void);
   //--- public methods
   static void       MLPECreate0(const int nin,const int nout,const int ensemblesize,CMLPEnsemble &ensemble);
   static void       MLPECreate1(const int nin,const int nhid,const int nout,const int ensemblesize,CMLPEnsemble &ensemble);
   static void       MLPECreate2(const int nin,const int nhid1,const int nhid2,const int nout,const int ensemblesize,CMLPEnsemble &ensemble);
   static void       MLPECreateB0(const int nin,const int nout,const double b,const double d,const int ensemblesize,CMLPEnsemble &ensemble);
   static void       MLPECreateB1(const int nin,const int nhid,const int nout,const double b,const double d,const int ensemblesize,CMLPEnsemble &ensemble);
   static void       MLPECreateB2(const int nin,const int nhid1,const int nhid2,const int nout,const double b,const double d,const int ensemblesize,CMLPEnsemble &ensemble);
   static void       MLPECreateR0(const int nin,const int nout,const double a,const double b,const int ensemblesize,CMLPEnsemble &ensemble);
   static void       MLPECreateR1(const int nin,const int nhid,const int nout,const double a,const double b,const int ensemblesize,CMLPEnsemble &ensemble);
   static void       MLPECreateR2(const int nin,const int nhid1,const int nhid2,const int nout,const double a,const double b,const int ensemblesize,CMLPEnsemble &ensemble);
   static void       MLPECreateC0(const int nin,const int nout,const int ensemblesize,CMLPEnsemble &ensemble);
   static void       MLPECreateC1(const int nin,const int nhid,const int nout,const int ensemblesize,CMLPEnsemble &ensemble);
   static void       MLPECreateC2(const int nin,const int nhid1,const int nhid2,const int nout,const int ensemblesize,CMLPEnsemble &ensemble);
   static void       MLPECreateFromNetwork(CMultilayerPerceptron &network,const int ensemblesize,CMLPEnsemble &ensemble);
   static void       MLPECopy(CMLPEnsemble &ensemble1,CMLPEnsemble &ensemble2);
   static void       MLPESerialize(CMLPEnsemble &ensemble,double &ra[],int &rlen);
   static void       MLPEUnserialize(double &ra[],CMLPEnsemble &ensemble);
   static void       MLPERandomize(CMLPEnsemble &ensemble);
   static void       MLPEProperties(CMLPEnsemble &ensemble,int &nin,int &nout);
   static bool       MLPEIsSoftMax(CMLPEnsemble &ensemble);
   static void       MLPEProcess(CMLPEnsemble &ensemble,double &x[],double &y[]);
   static void       MLPEProcessI(CMLPEnsemble &ensemble,double &x[],double &y[]);
   static double     MLPERelClsError(CMLPEnsemble &ensemble,CMatrixDouble &xy,const int npoints);
   static double     MLPEAvgCE(CMLPEnsemble &ensemble,CMatrixDouble &xy,const int npoints);
   static double     MLPERMSError(CMLPEnsemble &ensemble,CMatrixDouble &xy,const int npoints);
   static double     MLPEAvgError(CMLPEnsemble &ensemble,CMatrixDouble &xy,const int npoints);
   static double     MLPEAvgRelError(CMLPEnsemble &ensemble,CMatrixDouble &xy,const int npoints);
   static void       MLPEBaggingLM(CMLPEnsemble &ensemble,CMatrixDouble &xy,const int npoints,const double decay,const int restarts,int &info,CMLPReport &rep,CMLPCVReport &ooberrors);
   static void       MLPEBaggingLBFGS(CMLPEnsemble &ensemble,CMatrixDouble &xy,const int npoints,const double decay,const int restarts,const double wstep,const int maxits,int &info,CMLPReport &rep,CMLPCVReport &ooberrors);
   static void       MLPETrainES(CMLPEnsemble &ensemble,CMatrixDouble &xy,const int npoints,const double decay,const int restarts,int &info,CMLPReport &rep);
  };
//+------------------------------------------------------------------+
//| Initialize constants                                             |
//+------------------------------------------------------------------+
const int CMLPE::m_mlpntotaloffset=3;
const int CMLPE::m_mlpevnum=9;
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMLPE::CMLPE(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMLPE::~CMLPE(void)
  {

  }
//+------------------------------------------------------------------+
//| Like MLPCreate0, but for ensembles.                              |
//+------------------------------------------------------------------+
static void CMLPE::MLPECreate0(const int nin,const int nout,const int ensemblesize,
                               CMLPEnsemble &ensemble)
  {
//--- object of class
   CMultilayerPerceptron net;
//--- function call
   CMLPBase::MLPCreate0(nin,nout,net);
//--- function call
   MLPECreateFromNetwork(net,ensemblesize,ensemble);
  }
//+------------------------------------------------------------------+
//| Like MLPCreate1, but for ensembles.                              |
//+------------------------------------------------------------------+
static void CMLPE::MLPECreate1(const int nin,const int nhid,const int nout,
                               const int ensemblesize,CMLPEnsemble &ensemble)
  {
//--- object of class
   CMultilayerPerceptron net;
//--- function call
   CMLPBase::MLPCreate1(nin,nhid,nout,net);
//--- function call
   MLPECreateFromNetwork(net,ensemblesize,ensemble);
  }
//+------------------------------------------------------------------+
//| Like MLPCreate2, but for ensembles.                              |
//+------------------------------------------------------------------+
static void CMLPE::MLPECreate2(const int nin,const int nhid1,const int nhid2,
                               const int nout,const int ensemblesize,
                               CMLPEnsemble &ensemble)
  {
//--- object of class
   CMultilayerPerceptron net;
//--- function call
   CMLPBase::MLPCreate2(nin,nhid1,nhid2,nout,net);
//--- function call
   MLPECreateFromNetwork(net,ensemblesize,ensemble);
  }
//+------------------------------------------------------------------+
//| Like MLPCreateB0, but for ensembles.                             |
//+------------------------------------------------------------------+
static void CMLPE::MLPECreateB0(const int nin,const int nout,const double b,
                                const double d,const int ensemblesize,
                                CMLPEnsemble &ensemble)
  {
//--- object of class
   CMultilayerPerceptron net;
//--- function call
   CMLPBase::MLPCreateB0(nin,nout,b,d,net);
//--- function call
   MLPECreateFromNetwork(net,ensemblesize,ensemble);
  }
//+------------------------------------------------------------------+
//| Like MLPCreateB1, but for ensembles.                             |
//+------------------------------------------------------------------+
static void CMLPE::MLPECreateB1(const int nin,const int nhid,const int nout,
                                const double b,const double d,const int ensemblesize,
                                CMLPEnsemble &ensemble)
  {
//--- object of class
   CMultilayerPerceptron net;
//--- function call
   CMLPBase::MLPCreateB1(nin,nhid,nout,b,d,net);
//--- function call
   MLPECreateFromNetwork(net,ensemblesize,ensemble);
  }
//+------------------------------------------------------------------+
//| Like MLPCreateB2, but for ensembles.                             |
//+------------------------------------------------------------------+
static void CMLPE::MLPECreateB2(const int nin,const int nhid1,const int nhid2,
                                const int nout,const double b,const double d,
                                const int ensemblesize,CMLPEnsemble &ensemble)
  {
//--- object of class
   CMultilayerPerceptron net;
//--- function call
   CMLPBase::MLPCreateB2(nin,nhid1,nhid2,nout,b,d,net);
//--- function call
   MLPECreateFromNetwork(net,ensemblesize,ensemble);
  }
//+------------------------------------------------------------------+
//| Like MLPCreateR0, but for ensembles.                             |
//+------------------------------------------------------------------+
static void CMLPE::MLPECreateR0(const int nin,const int nout,const double a,
                                const double b,const int ensemblesize,
                                CMLPEnsemble &ensemble)
  {
//--- object of class
   CMultilayerPerceptron net;
//--- function call
   CMLPBase::MLPCreateR0(nin,nout,a,b,net);
//--- function call
   MLPECreateFromNetwork(net,ensemblesize,ensemble);
  }
//+------------------------------------------------------------------+
//| Like MLPCreateR1, but for ensembles.                             |
//+------------------------------------------------------------------+
static void CMLPE::MLPECreateR1(const int nin,const int nhid,const int nout,
                                const double a,const double b,
                                const int ensemblesize,CMLPEnsemble &ensemble)
  {
//--- object of class
   CMultilayerPerceptron net;
//--- function call
   CMLPBase::MLPCreateR1(nin,nhid,nout,a,b,net);
//--- function call
   MLPECreateFromNetwork(net,ensemblesize,ensemble);
  }
//+------------------------------------------------------------------+
//| Like MLPCreateR2, but for ensembles.                             |
//+------------------------------------------------------------------+
static void CMLPE::MLPECreateR2(const int nin,const int nhid1,const int nhid2,
                                const int nout,const double a,const double b,
                                const int ensemblesize,CMLPEnsemble &ensemble)
  {
//--- object of class
   CMultilayerPerceptron net;
//--- function call
   CMLPBase::MLPCreateR2(nin,nhid1,nhid2,nout,a,b,net);
//--- function call
   MLPECreateFromNetwork(net,ensemblesize,ensemble);
  }
//+------------------------------------------------------------------+
//| Like MLPCreateC0, but for ensembles.                             |
//+------------------------------------------------------------------+
static void CMLPE::MLPECreateC0(const int nin,const int nout,const int ensemblesize,
                                CMLPEnsemble &ensemble)
  {
//--- object of class
   CMultilayerPerceptron net;
//--- function call
   CMLPBase::MLPCreateC0(nin,nout,net);
//--- function call
   MLPECreateFromNetwork(net,ensemblesize,ensemble);
  }
//+------------------------------------------------------------------+
//| Like MLPCreateC1, but for ensembles.                             |
//+------------------------------------------------------------------+
static void CMLPE::MLPECreateC1(const int nin,const int nhid,const int nout,
                                const int ensemblesize,CMLPEnsemble &ensemble)
  {
//--- object of class
   CMultilayerPerceptron net;
//--- function call
   CMLPBase::MLPCreateC1(nin,nhid,nout,net);
//--- function call
   MLPECreateFromNetwork(net,ensemblesize,ensemble);
  }
//+------------------------------------------------------------------+
//| Like MLPCreateC2, but for ensembles.                             |
//+------------------------------------------------------------------+
static void CMLPE::MLPECreateC2(const int nin,const int nhid1,const int nhid2,
                                const int nout,const int ensemblesize,
                                CMLPEnsemble &ensemble)
  {
//--- object of class
   CMultilayerPerceptron net;
//--- function call
   CMLPBase::MLPCreateC2(nin,nhid1,nhid2,nout,net);
//--- function call
   MLPECreateFromNetwork(net,ensemblesize,ensemble);
  }
//+------------------------------------------------------------------+
//| Creates ensemble from network. Only network geometry is copied.  |
//+------------------------------------------------------------------+
static void CMLPE::MLPECreateFromNetwork(CMultilayerPerceptron &network,
                                         const int ensemblesize,
                                         CMLPEnsemble &ensemble)
  {
//--- create variables
   int i=0;
   int ccount=0;
   int i_=0;
   int i1_=0;
//--- check
   if(!CAp::Assert(ensemblesize>0,__FUNCTION__+": incorrect ensemble size!"))
      return;
//--- network properties
   CMLPBase::MLPProperties(network,ensemble.m_nin,ensemble.m_nout,ensemble.m_wcount);
//--- check
   if(CMLPBase::MLPIsSoftMax(network))
      ccount=ensemble.m_nin;
   else
      ccount=ensemble.m_nin+ensemble.m_nout;
//--- change values
   ensemble.m_postprocessing=false;
   ensemble.m_issoftmax=CMLPBase::MLPIsSoftMax(network);
   ensemble.m_ensemblesize=ensemblesize;
//--- structure information
   ArrayResizeAL(ensemble.m_structinfo,network.m_structinfo[0]);
//--- copy
   for(i=0;i<=network.m_structinfo[0]-1;i++)
      ensemble.m_structinfo[i]=network.m_structinfo[i];
//--- weights,means,sigmas
   ArrayResizeAL(ensemble.m_weights,ensemblesize*ensemble.m_wcount);
   ArrayResizeAL(ensemble.m_columnmeans,ensemblesize*ccount);
   ArrayResizeAL(ensemble.m_columnsigmas,ensemblesize*ccount);
//--- calculation
   for(i=0;i<=ensemblesize*ensemble.m_wcount-1;i++)
      ensemble.m_weights[i]=CMath::RandomReal()-0.5;
//--- calculation
   for(i=0;i<=ensemblesize-1;i++)
     {
      i1_=-(i*ccount);
      for(i_=i*ccount;i_<=(i+1)*ccount-1;i_++)
         ensemble.m_columnmeans[i_]=network.m_columnmeans[i_+i1_];
      i1_=-(i*ccount);
      for(i_=i*ccount;i_<=(i+1)*ccount-1;i_++)
         ensemble.m_columnsigmas[i_]=network.m_columnsigmas[i_+i1_];
     }
//--- serialized part
   CMLPBase::MLPSerializeOld(network,ensemble.m_serializedmlp,ensemble.m_serializedlen);
//--- temporaries,internal buffers
   ArrayResizeAL(ensemble.m_tmpweights,ensemble.m_wcount);
   ArrayResizeAL(ensemble.m_tmpmeans,ccount);
   ArrayResizeAL(ensemble.m_tmpsigmas,ccount);
   ArrayResizeAL(ensemble.m_neurons,ensemble.m_structinfo[m_mlpntotaloffset]);
   ArrayResizeAL(ensemble.m_dfdnet,ensemble.m_structinfo[m_mlpntotaloffset]);
   ArrayResizeAL(ensemble.m_y,ensemble.m_nout);
  }
//+------------------------------------------------------------------+
//| Copying of MLPEnsemble strucure                                  |
//| INPUT PARAMETERS:                                                |
//|     Ensemble1 -   original                                       |
//| OUTPUT PARAMETERS:                                               |
//|     Ensemble2 -   copy                                           |
//+------------------------------------------------------------------+
static void CMLPE::MLPECopy(CMLPEnsemble &ensemble1,CMLPEnsemble &ensemble2)
  {
//--- create variables
   int i=0;
   int ssize=0;
   int ccount=0;
   int ntotal=0;
   int i_=0;
//--- Unload info
   ssize=ensemble1.m_structinfo[0];
//--- check
   if(ensemble1.m_issoftmax)
      ccount=ensemble1.m_nin;
   else
      ccount=ensemble1.m_nin+ensemble1.m_nout;
//--- change value
   ntotal=ensemble1.m_structinfo[m_mlpntotaloffset];
//--- Allocate space
   ArrayResizeAL(ensemble2.m_structinfo,ssize);
   ArrayResizeAL(ensemble2.m_weights,ensemble1.m_ensemblesize*ensemble1.m_wcount);
   ArrayResizeAL(ensemble2.m_columnmeans,ensemble1.m_ensemblesize*ccount);
   ArrayResizeAL(ensemble2.m_columnsigmas,ensemble1.m_ensemblesize*ccount);
   ArrayResizeAL(ensemble2.m_tmpweights,ensemble1.m_wcount);
   ArrayResizeAL(ensemble2.m_tmpmeans,ccount);
   ArrayResizeAL(ensemble2.m_tmpsigmas,ccount);
   ArrayResizeAL(ensemble2.m_serializedmlp,ensemble1.m_serializedlen);
   ArrayResizeAL(ensemble2.m_neurons,ntotal);
   ArrayResizeAL(ensemble2.m_dfdnet,ntotal);
   ArrayResizeAL(ensemble2.m_y,ensemble1.m_nout);
//--- Copy
   ensemble2.m_nin=ensemble1.m_nin;
   ensemble2.m_nout=ensemble1.m_nout;
   ensemble2.m_wcount=ensemble1.m_wcount;
   ensemble2.m_ensemblesize=ensemble1.m_ensemblesize;
   ensemble2.m_issoftmax=ensemble1.m_issoftmax;
   ensemble2.m_postprocessing=ensemble1.m_postprocessing;
   ensemble2.m_serializedlen=ensemble1.m_serializedlen;
//--- copy
   for(i=0;i<=ssize-1;i++)
      ensemble2.m_structinfo[i]=ensemble1.m_structinfo[i];
   for(i_=0;i_<=ensemble1.m_ensemblesize*ensemble1.m_wcount-1;i_++)
      ensemble2.m_weights[i_]=ensemble1.m_weights[i_];
   for(i_=0;i_<=ensemble1.m_ensemblesize*ccount-1;i_++)
      ensemble2.m_columnmeans[i_]=ensemble1.m_columnmeans[i_];
   for(i_=0;i_<=ensemble1.m_ensemblesize*ccount-1;i_++)
      ensemble2.m_columnsigmas[i_]=ensemble1.m_columnsigmas[i_];
   for(i_=0;i_<=ensemble1.m_serializedlen-1;i_++)
      ensemble2.m_serializedmlp[i_]=ensemble1.m_serializedmlp[i_];
  }
//+------------------------------------------------------------------+
//| Serialization of MLPEnsemble strucure                            |
//| INPUT PARAMETERS:                                                |
//|     Ensemble-   original                                         |
//| OUTPUT PARAMETERS:                                               |
//|     RA      -   array of real numbers which stores ensemble,     |
//|                 array[0..RLen-1]                                 |
//|     RLen    -   RA lenght                                        |
//+------------------------------------------------------------------+
static void CMLPE::MLPESerialize(CMLPEnsemble &ensemble,double &ra[],int &rlen)
  {
//--- create variables
   int i=0;
   int ssize=0;
   int ntotal=0;
   int ccount=0;
   int hsize=0;
   int offs=0;
   int i_=0;
   int i1_=0;
//--- initialization
   rlen=0;
   hsize=13;
   ssize=ensemble.m_structinfo[0];
//--- check
   if(ensemble.m_issoftmax)
      ccount=ensemble.m_nin;
   else
      ccount=ensemble.m_nin+ensemble.m_nout;
//--- change values
   ntotal=ensemble.m_structinfo[m_mlpntotaloffset];
   rlen=hsize+ssize+ensemble.m_ensemblesize*ensemble.m_wcount+2*ccount*ensemble.m_ensemblesize+ensemble.m_serializedlen;
//---  RA format:
//---  [0]     RLen
//---  [1]     Version (MLPEVNum)
//---  [2]     EnsembleSize
//---  [3]     NIn
//---  [4]     NOut
//---  [5]     WCount
//---  [6]     IsSoftmax 0/1
//---  [7]     PostProcessing 0/1
//---  [8]     sizeof(StructInfo)
//---  [9]     NTotal (sizeof(Neurons),sizeof(DFDNET))
//---  [10]    CCount (sizeof(ColumnMeans),sizeof(ColumnSigmas))
//---  [11]    data offset
//---  [12]    SerializedLen
//---  [..]    StructInfo
//---  [..]    Weights
//---  [..]    ColumnMeans
//---  [..]    ColumnSigmas
   ArrayResizeAL(ra,rlen);
//--- change values
   ra[0]=rlen;
   ra[1]=m_mlpevnum;
   ra[2]=ensemble.m_ensemblesize;
   ra[3]=ensemble.m_nin;
   ra[4]=ensemble.m_nout;
   ra[5]=ensemble.m_wcount;
//--- check
   if(ensemble.m_issoftmax)
      ra[6]=1;
   else
      ra[6]=0;
//--- check
   if(ensemble.m_postprocessing)
      ra[7]=1;
   else
      ra[7]=9;
//--- change values
   ra[8]=ssize;
   ra[9]=ntotal;
   ra[10]=ccount;
   ra[11]=hsize;
   ra[12]=ensemble.m_serializedlen;
//--- copy
   offs=hsize;
   for(i=offs;i<=offs+ssize-1;i++)
      ra[i]=ensemble.m_structinfo[i-offs];
//--- copy
   offs=offs+ssize;
   i1_=-offs;
   for(i_=offs;i_<=offs+ensemble.m_ensemblesize*ensemble.m_wcount-1;i_++)
      ra[i_]=ensemble.m_weights[i_+i1_];
//--- copy
   offs=offs+ensemble.m_ensemblesize*ensemble.m_wcount;
   i1_=-offs;
   for(i_=offs;i_<=offs+ensemble.m_ensemblesize*ccount-1;i_++)
      ra[i_]=ensemble.m_columnmeans[i_+i1_];
//--- copy
   offs=offs+ensemble.m_ensemblesize*ccount;
   i1_=-offs;
   for(i_=offs;i_<=offs+ensemble.m_ensemblesize*ccount-1;i_++)
      ra[i_]=ensemble.m_columnsigmas[i_+i1_];
//--- copy
   offs=offs+ensemble.m_ensemblesize*ccount;
   i1_=-offs;
   for(i_=offs;i_<=offs+ensemble.m_serializedlen-1;i_++)
      ra[i_]=ensemble.m_serializedmlp[i_+i1_];
   offs=offs+ensemble.m_serializedlen;
  }
//+------------------------------------------------------------------+
//| Unserialization of MLPEnsemble strucure                          |
//| INPUT PARAMETERS:                                                |
//|     RA      -   real array which stores ensemble                 |
//| OUTPUT PARAMETERS:                                               |
//|     Ensemble-   restored structure                               |
//+------------------------------------------------------------------+
static void CMLPE::MLPEUnserialize(double &ra[],CMLPEnsemble &ensemble)
  {
//--- create variables
   int i=0;
   int ssize=0;
   int ntotal=0;
   int ccount=0;
   int hsize=0;
   int offs=0;
   int i_=0;
   int i1_=0;
//--- check
   if(!CAp::Assert((int)MathRound(ra[1])==m_mlpevnum,__FUNCTION__+": incorrect array!"))
      return;
//--- load info
   hsize=13;
   ensemble.m_ensemblesize=(int)MathRound(ra[2]);
   ensemble.m_nin=(int)MathRound(ra[3]);
   ensemble.m_nout=(int)MathRound(ra[4]);
   ensemble.m_wcount=(int)MathRound(ra[5]);
   ensemble.m_issoftmax=(int)MathRound(ra[6])==1;
   ensemble.m_postprocessing=(int)MathRound(ra[7])==1;
   ssize=(int)MathRound(ra[8]);
   ntotal=(int)MathRound(ra[9]);
   ccount=(int)MathRound(ra[10]);
   offs=(int)MathRound(ra[11]);
   ensemble.m_serializedlen=(int)MathRound(ra[12]);
//---  Allocate arrays
   ArrayResizeAL(ensemble.m_structinfo,ssize);
   ArrayResizeAL(ensemble.m_weights,ensemble.m_ensemblesize*ensemble.m_wcount);
   ArrayResizeAL(ensemble.m_columnmeans,ensemble.m_ensemblesize*ccount);
   ArrayResizeAL(ensemble.m_columnsigmas,ensemble.m_ensemblesize*ccount);
   ArrayResizeAL(ensemble.m_tmpweights,ensemble.m_wcount);
   ArrayResizeAL(ensemble.m_tmpmeans,ccount);
   ArrayResizeAL(ensemble.m_tmpsigmas,ccount);
   ArrayResizeAL(ensemble.m_neurons,ntotal);
   ArrayResizeAL(ensemble.m_dfdnet,ntotal);
   ArrayResizeAL(ensemble.m_serializedmlp,ensemble.m_serializedlen);
   ArrayResizeAL(ensemble.m_y,ensemble.m_nout);
//--- load data
   for(i=offs;i<=offs+ssize-1;i++)
      ensemble.m_structinfo[i-offs]=(int)MathRound(ra[i]);
//--- copy
   offs=offs+ssize;
   i1_=offs;
   for(i_=0;i_<=ensemble.m_ensemblesize*ensemble.m_wcount-1;i_++)
      ensemble.m_weights[i_]=ra[i_+i1_];
//--- copy
   offs=offs+ensemble.m_ensemblesize*ensemble.m_wcount;
   i1_=offs;
   for(i_=0;i_<=ensemble.m_ensemblesize*ccount-1;i_++)
      ensemble.m_columnmeans[i_]=ra[i_+i1_];
//--- copy
   offs=offs+ensemble.m_ensemblesize*ccount;
   i1_=offs;
   for(i_=0;i_<=ensemble.m_ensemblesize*ccount-1;i_++)
      ensemble.m_columnsigmas[i_]=ra[i_+i1_];
//--- copy
   offs=offs+ensemble.m_ensemblesize*ccount;
   i1_=offs;
   for(i_=0;i_<=ensemble.m_serializedlen-1;i_++)
      ensemble.m_serializedmlp[i_]=ra[i_+i1_];
   offs=offs+ensemble.m_serializedlen;
  }
//+------------------------------------------------------------------+
//| Randomization of MLP ensemble                                    |
//+------------------------------------------------------------------+
static void CMLPE::MLPERandomize(CMLPEnsemble &ensemble)
  {
//--- create a variable
   int i=0;
//--- calculation
   for(i=0;i<=ensemble.m_ensemblesize*ensemble.m_wcount-1;i++)
      ensemble.m_weights[i]=CMath::RandomReal()-0.5;
  }
//+------------------------------------------------------------------+
//| Return ensemble properties (number of inputs and outputs).       |
//+------------------------------------------------------------------+
static void CMLPE::MLPEProperties(CMLPEnsemble &ensemble,int &nin,int &nout)
  {
//--- change values
   nin=ensemble.m_nin;
   nout=ensemble.m_nout;
  }
//+------------------------------------------------------------------+
//| Return normalization type (whether ensemble is SOFTMAX-normalized|
//| or not).                                                         |
//+------------------------------------------------------------------+
static bool CMLPE::MLPEIsSoftMax(CMLPEnsemble &ensemble)
  {
//--- return result
   return(ensemble.m_issoftmax);
  }
//+------------------------------------------------------------------+
//| Procesing                                                        |
//| INPUT PARAMETERS:                                                |
//|     Ensemble-   neural networks ensemble                         |
//|     X       -   input vector,  array[0..NIn-1].                  |
//|     Y       -   (possibly) preallocated buffer; if size of Y is  |
//|                 less than NOut, it will be reallocated. If it is |
//|                 large enough, it is NOT reallocated, so we can   |
//|                 save some time on reallocation.                  |
//| OUTPUT PARAMETERS:                                               |
//|     Y       -   result. Regression estimate when solving         |
//|                 regression task, vector of posterior             |
//|                 probabilities for classification task.           |
//+------------------------------------------------------------------+
static void CMLPE::MLPEProcess(CMLPEnsemble &ensemble,double &x[],double &y[])
  {
//--- create variables
   int    i=0;
   int    es=0;
   int    wc=0;
   int    cc=0;
   double v=0;
   int    i_=0;
   int    i1_=0;
//--- check
   if(CAp::Len(y)<ensemble.m_nout)
      ArrayResizeAL(y,ensemble.m_nout);
//--- initialization
   es=ensemble.m_ensemblesize;
   wc=ensemble.m_wcount;
//--- check
   if(ensemble.m_issoftmax)
      cc=ensemble.m_nin;
   else
      cc=ensemble.m_nin+ensemble.m_nout;
//--- initialization
   v=1.0/(double)es;
   for(i=0;i<=ensemble.m_nout-1;i++)
      y[i]=0;
//--- calculation
   for(i=0;i<=es-1;i++)
     {
      i1_=i*wc;
      for(i_=0;i_<=wc-1;i_++)
         ensemble.m_tmpweights[i_]=ensemble.m_weights[i_+i1_];
      i1_=i*cc;
      for(i_=0;i_<=cc-1;i_++)
         ensemble.m_tmpmeans[i_]=ensemble.m_columnmeans[i_+i1_];
      i1_=i*cc;
      for(i_=0;i_<=cc-1;i_++)
         ensemble.m_tmpsigmas[i_]=ensemble.m_columnsigmas[i_+i1_];
      //--- function call
      CMLPBase::MLPInternalProcessVector(ensemble.m_structinfo,ensemble.m_tmpweights,ensemble.m_tmpmeans,ensemble.m_tmpsigmas,ensemble.m_neurons,ensemble.m_dfdnet,x,ensemble.m_y);
      //--- change values
      for(i_=0;i_<=ensemble.m_nout-1;i_++)
         y[i_]=y[i_]+v*ensemble.m_y[i_];
     }
  }
//+------------------------------------------------------------------+
//| 'interactive' variant of MLPEProcess for languages like Python   |
//| which support constructs like "Y = MLPEProcess(LM,X)" and        |
//| interactive mode of the interpreter                              |
//| This function allocates new array on each call, so it is         |
//| significantly slower than its 'non-interactive' counterpart, but |
//| it is more convenient when you call it from command line.        |
//+------------------------------------------------------------------+
static void CMLPE::MLPEProcessI(CMLPEnsemble &ensemble,double &x[],double &y[])
  {
//--- function call
   MLPEProcess(ensemble,x,y);
  }
//+------------------------------------------------------------------+
//| Relative classification error on the test set                    |
//| INPUT PARAMETERS:                                                |
//|     Ensemble-   ensemble                                         |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     percent of incorrectly classified cases.                     |
//|     Works both for classifier betwork and for regression networks|
//|     which are used as classifiers.                               |
//+------------------------------------------------------------------+
static double CMLPE::MLPERelClsError(CMLPEnsemble &ensemble,CMatrixDouble &xy,
                                     const int npoints)
  {
//--- create variables
   double relcls=0;
   double avgce=0;
   double rms=0;
   double avg=0;
   double avgrel=0;
//--- function call
   MLPEAllErrors(ensemble,xy,npoints,relcls,avgce,rms,avg,avgrel);
//--- return result
   return(relcls);
  }
//+------------------------------------------------------------------+
//| Average cross-entropy (in bits per element) on the test set      |
//| INPUT PARAMETERS:                                                |
//|     Ensemble-   ensemble                                         |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     CrossEntropy/(NPoints*LN(2)).                                |
//|     Zero if ensemble solves regression task.                     |
//+------------------------------------------------------------------+
static double CMLPE::MLPEAvgCE(CMLPEnsemble &ensemble,CMatrixDouble &xy,
                               const int npoints)
  {
//--- create variables
   double relcls=0;
   double avgce=0;
   double rms=0;
   double avg=0;
   double avgrel=0;
//--- function call
   MLPEAllErrors(ensemble,xy,npoints,relcls,avgce,rms,avg,avgrel);
//--- return result
   return(avgce);
  }
//+------------------------------------------------------------------+
//| RMS error on the test set                                        |
//| INPUT PARAMETERS:                                                |
//|     Ensemble-   ensemble                                         |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     root mean square error.                                      |
//|     Its meaning for regression task is obvious. As for           |
//|     classification task RMS error means error when estimating    |
//|     posterior probabilities.                                     |
//+------------------------------------------------------------------+
static double CMLPE::MLPERMSError(CMLPEnsemble &ensemble,CMatrixDouble &xy,
                                  const int npoints)
  {
//--- create variables
   double relcls=0;
   double avgce=0;
   double rms=0;
   double avg=0;
   double avgrel=0;
//--- function call
   MLPEAllErrors(ensemble,xy,npoints,relcls,avgce,rms,avg,avgrel);
//--- return result
   return(rms);
  }
//+------------------------------------------------------------------+
//| Average error on the test set                                    |
//| INPUT PARAMETERS:                                                |
//|     Ensemble-   ensemble                                         |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     Its meaning for regression task is obvious. As for           |
//|     classification task it means average error when estimating   |
//|     posterior probabilities.                                     |
//+------------------------------------------------------------------+
static double CMLPE::MLPEAvgError(CMLPEnsemble &ensemble,CMatrixDouble &xy,
                                  const int npoints)
  {
//--- create variables
   double result=0;
   double relcls=0;
   double avgce=0;
   double rms=0;
   double avg=0;
   double avgrel=0;
//--- function call
   MLPEAllErrors(ensemble,xy,npoints,relcls,avgce,rms,avg,avgrel);
//--- return result
   return(avg);
  }
//+------------------------------------------------------------------+
//| Average relative error on the test set                           |
//| INPUT PARAMETERS:                                                |
//|     Ensemble-   ensemble                                         |
//|     XY      -   test set                                         |
//|     NPoints -   test set size                                    |
//| RESULT:                                                          |
//|     Its meaning for regression task is obvious. As for           |
//|     classification task it means average relative error when     |
//|     estimating posterior probabilities.                          |
//+------------------------------------------------------------------+
static double CMLPE::MLPEAvgRelError(CMLPEnsemble &ensemble,CMatrixDouble &xy,
                                     const int npoints)
  {
//--- create variables
   double result=0;
   double relcls=0;
   double avgce=0;
   double rms=0;
   double avg=0;
   double avgrel=0;
//--- function call
   MLPEAllErrors(ensemble,xy,npoints,relcls,avgce,rms,avg,avgrel);
//--- return result
   return(avgrel);
  }
//+------------------------------------------------------------------+
//| Training neural networks ensemble using  bootstrap  aggregating  |
//| (bagging).                                                       |
//| Modified Levenberg-Marquardt algorithm is used as base training  |
//| method.                                                          |
//| INPUT PARAMETERS:                                                |
//|     Ensemble    -   model with initialized geometry              |
//|     XY          -   training set                                 |
//|     NPoints     -   training set size                            |
//|     Decay       -   weight decay coefficient, >=0.001            |
//|     Restarts    -   restarts, >0.                                |
//| OUTPUT PARAMETERS:                                               |
//|     Ensemble    -   trained model                                |
//|     Info        -   return code:                                 |
//|                     * -2, if there is a point with class number  |
//|                           outside of [0..NClasses-1].            |
//|                     * -1, if incorrect parameters was passed     |
//|                           (NPoints<0, Restarts<1).               |
//|                     *  2, if task has been solved.               |
//|     Rep         -   training report.                             |
//|     OOBErrors   -   out-of-bag generalization error estimate     |
//+------------------------------------------------------------------+
static void CMLPE::MLPEBaggingLM(CMLPEnsemble &ensemble,CMatrixDouble &xy,
                                 const int npoints,const double decay,
                                 const int restarts,int &info,CMLPReport &rep,
                                 CMLPCVReport &ooberrors)
  {
//--- initialization
   info=0;
//--- function call
   MLPEBaggingInternal(ensemble,xy,npoints,decay,restarts,0.0,0,true,info,rep,ooberrors);
  }
//+------------------------------------------------------------------+
//| Training neural networks ensemble using  bootstrap aggregating   |
//| (bagging). L-BFGS algorithm is used as base training method.     |
//| INPUT PARAMETERS:                                                |
//|     Ensemble    -   model with initialized geometry              |
//|     XY          -   training set                                 |
//|     NPoints     -   training set size                            |
//|     Decay       -   weight decay coefficient, >=0.001            |
//|     Restarts    -   restarts, >0.                                |
//|     WStep       -   stopping criterion, same as in MLPTrainLBFGS |
//|     MaxIts      -   stopping criterion, same as in MLPTrainLBFGS |
//| OUTPUT PARAMETERS:                                               |
//|     Ensemble    -   trained model                                |
//|     Info        -   return code:                                 |
//|                     * -8, if both WStep=0 and MaxIts=0           |
//|                     * -2, if there is a point with class number  |
//|                           outside of [0..NClasses-1].            |
//|                     * -1, if incorrect parameters was passed     |
//|                           (NPoints<0, Restarts<1).               |
//|                     *  2, if task has been solved.               |
//|     Rep         -   training report.                             |
//|     OOBErrors   -   out-of-bag generalization error estimate     |
//+------------------------------------------------------------------+
static void CMLPE::MLPEBaggingLBFGS(CMLPEnsemble &ensemble,CMatrixDouble &xy,
                                    const int npoints,const double decay,
                                    const int restarts,const double wstep,
                                    const int maxits,int &info,CMLPReport &rep,
                                    CMLPCVReport &ooberrors)
  {
//--- initialization
   info=0;
//--- function call
   MLPEBaggingInternal(ensemble,xy,npoints,decay,restarts,wstep,maxits,false,info,rep,ooberrors);
  }
//+------------------------------------------------------------------+
//| Training neural networks ensemble using early stopping.          |
//| INPUT PARAMETERS:                                                |
//|     Ensemble    -   model with initialized geometry              |
//|     XY          -   training set                                 |
//|     NPoints     -   training set size                            |
//|     Decay       -   weight decay coefficient, >=0.001            |
//|     Restarts    -   restarts, >0.                                |
//| OUTPUT PARAMETERS:                                               |
//|     Ensemble    -   trained model                                |
//|     Info        -   return code:                                 |
//|                     * -2, if there is a point with class number  |
//|                           outside of [0..NClasses-1].            |
//|                     * -1, if incorrect parameters was passed     |
//|                           (NPoints<0, Restarts<1).               |
//|                     *  6, if task has been solved.               |
//|     Rep         -   training report.                             |
//|     OOBErrors   -   out-of-bag generalization error estimate     |
//+------------------------------------------------------------------+
static void CMLPE::MLPETrainES(CMLPEnsemble &ensemble,CMatrixDouble &xy,
                               const int npoints,const double decay,
                               const int restarts,int &info,CMLPReport &rep)
  {
//--- create variables
   int i=0;
   int k=0;
   int ccount=0;
   int pcount=0;
   int trnsize=0;
   int valsize=0;
   int tmpinfo=0;
   int i_=0;
   int i1_=0;
//--- create matrix
   CMatrixDouble trnxy;
   CMatrixDouble valxy;
//--- objects of classes
   CMultilayerPerceptron network;
   CMLPReport            tmprep;
//--- initialization
   info=0;
//--- check
   if((npoints<2 || restarts<1) || decay<0.0)
     {
      info=-1;
      return;
     }
//--- check
   if(ensemble.m_issoftmax)
     {
      for(i=0;i<=npoints-1;i++)
        {
         //--- check
         if((int)MathRound(xy[i][ensemble.m_nin])<0 || (int)MathRound(xy[i][ensemble.m_nin])>=ensemble.m_nout)
           {
            info=-2;
            return;
           }
        }
     }
//--- change value
   info=6;
//--- allocate
   if(ensemble.m_issoftmax)
     {
      ccount=ensemble.m_nin+1;
      pcount=ensemble.m_nin;
     }
   else
     {
      ccount=ensemble.m_nin+ensemble.m_nout;
      pcount=ensemble.m_nin+ensemble.m_nout;
     }
//--- allocation
   trnxy.Resize(npoints,ccount);
   valxy.Resize(npoints,ccount);
//--- function call
   CMLPBase::MLPUnserializeOld(ensemble.m_serializedmlp,network);
//--- change values
   rep.m_ngrad=0;
   rep.m_nhess=0;
   rep.m_ncholesky=0;
//--- train networks
   for(k=0;k<=ensemble.m_ensemblesize-1;k++)
     {
      //--- Split set
      do
        {
         trnsize=0;
         valsize=0;
         for(i=0;i<=npoints-1;i++)
           {
            //--- check
            if(CMath::RandomReal()<0.66)
              {
               //--- Assign sample to training set
               for(i_=0;i_<=ccount-1;i_++)
                  trnxy[trnsize].Set(i_,xy[i][i_]);
               trnsize=trnsize+1;
              }
            else
              {
               //--- Assign sample to validation set
               for(i_=0;i_<=ccount-1;i_++)
                  valxy[valsize].Set(i_,xy[i][i_]);
               valsize=valsize+1;
              }
           }
        }
      while(!(trnsize!=0 && valsize!=0));
      //--- Train
      CMLPTrain::MLPTrainES(network,trnxy,trnsize,valxy,valsize,decay,restarts,tmpinfo,tmprep);
      //--- check
      if(tmpinfo<0)
        {
         info=tmpinfo;
         return;
        }
      //--- save results
      i1_=-(k*ensemble.m_wcount);
      for(i_=k*ensemble.m_wcount;i_<=(k+1)*ensemble.m_wcount-1;i_++)
         ensemble.m_weights[i_]=network.m_weights[i_+i1_];
      i1_=-(k*pcount);
      for(i_=k*pcount;i_<=(k+1)*pcount-1;i_++)
         ensemble.m_columnmeans[i_]=network.m_columnmeans[i_+i1_];
      i1_=-(k*pcount);
      for(i_=k*pcount;i_<=(k+1)*pcount-1;i_++)
         ensemble.m_columnsigmas[i_]=network.m_columnsigmas[i_+i1_];
      //--- change values
      rep.m_ngrad=rep.m_ngrad+tmprep.m_ngrad;
      rep.m_nhess=rep.m_nhess+tmprep.m_nhess;
      rep.m_ncholesky=rep.m_ncholesky+tmprep.m_ncholesky;
     }
  }
//+------------------------------------------------------------------+
//| Calculation of all types of errors                               |
//+------------------------------------------------------------------+
static void CMLPE::MLPEAllErrors(CMLPEnsemble &ensemble,CMatrixDouble &xy,
                                 const int npoints,double &relcls,
                                 double &avgce,double &rms,
                                 double &avg,double &avgrel)
  {
//--- create variables
   int i=0;
   int i_=0;
   int i1_=0;
//--- creating arrays
   double buf[];
   double workx[];
   double y[];
   double dy[];
//--- initialization
   relcls=0;
   avgce=0;
   rms=0;
   avg=0;
   avgrel=0;
//--- allocation
   ArrayResizeAL(workx,ensemble.m_nin);
   ArrayResizeAL(y,ensemble.m_nout);
//--- check
   if(ensemble.m_issoftmax)
     {
      //--- allocation
      ArrayResizeAL(dy,1);
      //--- function call
      CBdSS::DSErrAllocate(ensemble.m_nout,buf);
     }
   else
     {
      //--- allocation
      ArrayResizeAL(dy,ensemble.m_nout);
      //--- function call
      CBdSS::DSErrAllocate(-ensemble.m_nout,buf);
     }
//--- calculation
   for(i=0;i<=npoints-1;i++)
     {
      for(i_=0;i_<=ensemble.m_nin-1;i_++)
         workx[i_]=xy[i][i_];
      //--- function call
      MLPEProcess(ensemble,workx,y);
      //--- check
      if(ensemble.m_issoftmax)
         dy[0]=xy[i][ensemble.m_nin];
      else
        {
         i1_=ensemble.m_nin;
         for(i_=0;i_<=ensemble.m_nout-1;i_++)
            dy[i_]=xy[i][i_+i1_];
        }
      //--- function call
      CBdSS::DSErrAccumulate(buf,y,dy);
     }
//--- function call
   CBdSS::DSErrFinish(buf);
//--- change values
   relcls=buf[0];
   avgce=buf[1];
   rms=buf[2];
   avg=buf[3];
   avgrel=buf[4];
  }
//+------------------------------------------------------------------+
//| Internal bagging subroutine.                                     |
//+------------------------------------------------------------------+
static void CMLPE::MLPEBaggingInternal(CMLPEnsemble &ensemble,CMatrixDouble &xy,
                                       const int npoints,const double decay,
                                       const int restarts,const double wstep,
                                       const int maxits,const bool lmalgorithm,
                                       int &info,CMLPReport &rep,CMLPCVReport &ooberrors)
  {
//--- create variables
   int nin=0;
   int nout=0;
   int ccnt=0;
   int pcnt=0;
   int i=0;
   int j=0;
   int k=0;
   double v=0;
   int i_=0;
   int i1_=0;
//--- creating arrays
   bool   s[];
   int    oobcntbuf[];
   double x[];
   double y[];
   double dy[];
   double dsbuf[];
//--- create matrix
   CMatrixDouble xys;
   CMatrixDouble oobbuf;
//--- objects of classes
   CMLPReport            tmprep;
   CMultilayerPerceptron network;
//--- initialization
   info=0;
//--- Test for inputs
   if((!lmalgorithm && wstep==0.0) && maxits==0)
     {
      info=-8;
      return;
     }
//--- check
   if(((npoints<=0 || restarts<1) || wstep<0.0) || maxits<0)
     {
      info=-1;
      return;
     }
//--- check
   if(ensemble.m_issoftmax)
     {
      for(i=0;i<=npoints-1;i++)
        {
         //--- check
         if((int)MathRound(xy[i][ensemble.m_nin])<0 || (int)MathRound(xy[i][ensemble.m_nin])>=ensemble.m_nout)
           {
            info=-2;
            return;
           }
        }
     }
//--- allocate temporaries
   info=2;
   rep.m_ngrad=0;
   rep.m_nhess=0;
   rep.m_ncholesky=0;
   ooberrors.m_relclserror=0;
   ooberrors.m_avgce=0;
   ooberrors.m_rmserror=0;
   ooberrors.m_avgerror=0;
   ooberrors.m_avgrelerror=0;
   nin=ensemble.m_nin;
   nout=ensemble.m_nout;
//--- check
   if(ensemble.m_issoftmax)
     {
      ccnt=nin+1;
      pcnt=nin;
     }
   else
     {
      ccnt=nin+nout;
      pcnt=nin+nout;
     }
//--- allocation
   xys.Resize(npoints,ccnt);
   ArrayResizeAL(s,npoints);
   oobbuf.Resize(npoints,nout);
   ArrayResizeAL(oobcntbuf,npoints);
   ArrayResizeAL(x,nin);
   ArrayResizeAL(y,nout);
//--- check
   if(ensemble.m_issoftmax)
      ArrayResizeAL(dy,1);
   else
      ArrayResizeAL(dy,nout);
//--- initialization
   for(i=0;i<=npoints-1;i++)
     {
      for(j=0;j<=nout-1;j++)
         oobbuf[i].Set(j,0);
     }
   for(i=0;i<=npoints-1;i++)
      oobcntbuf[i]=0;
//--- function call
   CMLPBase::MLPUnserializeOld(ensemble.m_serializedmlp,network);
//--- main bagging cycle
   for(k=0;k<=ensemble.m_ensemblesize-1;k++)
     {
      //--- prepare dataset
      for(i=0;i<=npoints-1;i++)
         s[i]=false;
      for(i=0;i<=npoints-1;i++)
        {
         j=CMath::RandomInteger(npoints);
         s[j]=true;
         for(i_=0;i_<=ccnt-1;i_++)
            xys[i].Set(i_,xy[j][i_]);
        }
      //--- train
      if(lmalgorithm)
         CMLPTrain::MLPTrainLM(network,xys,npoints,decay,restarts,info,tmprep);
      else
         CMLPTrain::MLPTrainLBFGS(network,xys,npoints,decay,restarts,wstep,maxits,info,tmprep);
      //--- check
      if(info<0)
         return;
      //--- save results
      rep.m_ngrad=rep.m_ngrad+tmprep.m_ngrad;
      rep.m_nhess=rep.m_nhess+tmprep.m_nhess;
      rep.m_ncholesky=rep.m_ncholesky+tmprep.m_ncholesky;
      //--- copy
      i1_=-(k*ensemble.m_wcount);
      for(i_=k*ensemble.m_wcount;i_<=(k+1)*ensemble.m_wcount-1;i_++)
         ensemble.m_weights[i_]=network.m_weights[i_+i1_];
      //--- copy
      i1_=-(k*pcnt);
      for(i_=k*pcnt;i_<=(k+1)*pcnt-1;i_++)
         ensemble.m_columnmeans[i_]=network.m_columnmeans[i_+i1_];
      //--- copy
      i1_=-(k*pcnt);
      for(i_=k*pcnt;i_<=(k+1)*pcnt-1;i_++)
         ensemble.m_columnsigmas[i_]=network.m_columnsigmas[i_+i1_];
      //--- OOB estimates
      for(i=0;i<=npoints-1;i++)
        {
         //--- check
         if(!s[i])
           {
            for(i_=0;i_<=nin-1;i_++)
               x[i_]=xy[i][i_];
            //--- function call
            CMLPBase::MLPProcess(network,x,y);
            //--- change value
            for(i_=0;i_<=nout-1;i_++)
               oobbuf[i].Set(i_,oobbuf[i][i_]+y[i_]);
            oobcntbuf[i]=oobcntbuf[i]+1;
           }
        }
     }
//--- OOB estimates
   if(ensemble.m_issoftmax)
     {
      //--- function call
      CBdSS::DSErrAllocate(nout,dsbuf);
     }
   else
     {
      //--- function call
      CBdSS::DSErrAllocate(-nout,dsbuf);
     }
   for(i=0;i<=npoints-1;i++)
     {
      //--- check
      if(oobcntbuf[i]!=0)
        {
         v=1.0/(double)oobcntbuf[i];
         for(i_=0;i_<=nout-1;i_++)
            y[i_]=v*oobbuf[i][i_];
         //--- check
         if(ensemble.m_issoftmax)
            dy[0]=xy[i][nin];
         else
           {
            i1_=nin;
            for(i_=0;i_<=nout-1;i_++)
               dy[i_]=v*xy[i][i_+i1_];
           }
         //--- function call
         CBdSS::DSErrAccumulate(dsbuf,y,dy);
        }
     }
//--- function call
   CBdSS::DSErrFinish(dsbuf);
//--- change values
   ooberrors.m_relclserror=dsbuf[0];
   ooberrors.m_avgce=dsbuf[1];
   ooberrors.m_rmserror=dsbuf[2];
   ooberrors.m_avgerror=dsbuf[3];
   ooberrors.m_avgrelerror=dsbuf[4];
  }
//+------------------------------------------------------------------+
//| Principal components analysis                                    |
//+------------------------------------------------------------------+
class CPCAnalysis
  {
public:
   //--- constructor, destructor
                     CPCAnalysis(void);
                    ~CPCAnalysis(void);
   //--- method
   static void       PCABuildBasis(CMatrixDouble &x,const int npoints,const int nvars,int &info,double &s2[],CMatrixDouble &v);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CPCAnalysis::CPCAnalysis(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CPCAnalysis::~CPCAnalysis(void)
  {

  }
//+------------------------------------------------------------------+
//| Principal components analysis                                    |
//| Subroutine builds orthogonal basis where first axis corresponds  | 
//| to direction with maximum variance, second axis maximizes        |
//| variance in subspace orthogonal to first axis and so on.         |
//| It should be noted that, unlike LDA, PCA does not use class      |
//| labels.                                                          |
//| INPUT PARAMETERS:                                                |
//|     X           -   dataset, array[0..NPoints-1,0..NVars-1].     |
//|                     matrix contains ONLY INDEPENDENT VARIABLES.  |
//|     NPoints     -   dataset size, NPoints>=0                     |
//|     NVars       -   number of independent variables, NVars>=1    |
//| OUTPUT PARAMETERS:                                               |
//|     Info        -   return code:                                 |
//|                     * -4, if SVD subroutine haven't converged    |
//|                     * -1, if wrong parameters has been passed    |
//|                           (NPoints<0, NVars<1)                   |
//|                     *  1, if task is solved                      |
//|     S2          -   array[0..NVars-1]. variance values           |
//|                     corresponding to basis vectors.              |
//|     V           -   array[0..NVars-1,0..NVars-1]                 |
//|                     matrix, whose columns store basis vectors.   |
//+------------------------------------------------------------------+
static void CPCAnalysis::PCABuildBasis(CMatrixDouble &x,const int npoints,
                                       const int nvars,int &info,double &s2[],
                                       CMatrixDouble &v)
  {
//--- create variables
   int    i=0;
   int    j=0;
   double mean=0;
   double variance=0;
   double skewness=0;
   double kurtosis=0;
   int    i_=0;
//--- creating arrays
   double m[];
   double t[];
//--- create matrix
   CMatrixDouble a;
   CMatrixDouble u;
   CMatrixDouble vt;
//--- initialization
   info=0;
//--- Check input data
   if(npoints<0 || nvars<1)
     {
      info=-1;
      return;
     }
//--- change value
   info=1;
//--- Special case: NPoints=0
   if(npoints==0)
     {
      //--- allocation
      ArrayResizeAL(s2,nvars);
      v.Resize(nvars,nvars);
      //--- initialization
      for(i=0;i<=nvars-1;i++)
         s2[i]=0;
      for(i=0;i<=nvars-1;i++)
        {
         for(j=0;j<=nvars-1;j++)
           {
            //--- check
            if(i==j)
               v[i].Set(j,1);
            else
               v[i].Set(j,0);
           }
        }
      //--- exit the function
      return;
     }
//--- Calculate means
   ArrayResizeAL(m,nvars);
   ArrayResizeAL(t,npoints);
   for(j=0;j<=nvars-1;j++)
     {
      for(i_=0;i_<=npoints-1;i_++)
         t[i_]=x[i_][j];
      //--- function call
      CBaseStat::SampleMoments(t,npoints,mean,variance,skewness,kurtosis);
      m[j]=mean;
     }
//--- Center,apply SVD,prepare output
   a.Resize(MathMax(npoints,nvars),nvars);
//--- calculation
   for(i=0;i<=npoints-1;i++)
     {
      for(i_=0;i_<=nvars-1;i_++)
         a[i].Set(i_,x[i][i_]);
      for(i_=0;i_<=nvars-1;i_++)
         a[i].Set(i_,a[i][i_]-m[i_]);
     }
   for(i=npoints;i<=nvars-1;i++)
     {
      for(j=0;j<=nvars-1;j++)
         a[i].Set(j,0);
     }
//--- check
   if(!CSingValueDecompose::RMatrixSVD(a,MathMax(npoints,nvars),nvars,0,1,2,s2,u,vt))
     {
      info=-4;
      return;
     }
//--- check
   if(npoints!=1)
     {
      for(i=0;i<=nvars-1;i++)
         s2[i]=CMath::Sqr(s2[i])/(npoints-1);
     }
//--- allocation
   v.Resize(nvars,nvars);
//--- function call
   CBlas::CopyAndTranspose(vt,0,nvars-1,0,nvars-1,v,0,nvars-1,0,nvars-1);
  }
//+------------------------------------------------------------------+
