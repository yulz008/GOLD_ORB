//+------------------------------------------------------------------+
//|                                             specialfunctions.mqh |
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
#include "alglibinternal.mqh"
//+------------------------------------------------------------------+
//| Gamma function                                                   |
//+------------------------------------------------------------------+
class CGammaFunc
  {
private:
   static double     GammaStirlFunc(double x);

public:
                     CGammaFunc(void);
                    ~CGammaFunc(void);
   static double     GammaFunc(double x);
   static double     LnGamma(double x,double &sgngam);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CGammaFunc::CGammaFunc(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CGammaFunc::~CGammaFunc(void)
  {

  }
//+------------------------------------------------------------------+
//| Stirling function                                                |
//+------------------------------------------------------------------+
static double CGammaFunc::GammaStirlFunc(double x)
  {
//--- create variables
   double y;
   double w;
   double v;
   double stir;
//--- initialization
   w=1/x;
   stir=7.87311395793093628397E-4;
   stir=-2.29549961613378126380E-4+w*stir;
   stir=-2.68132617805781232825E-3+w*stir;
   stir=3.47222221605458667310E-3+w*stir;
   stir=8.33333333333482257126E-2+w*stir;
   w=1+w*stir;
   y=MathExp(x);
//--- check
   if(x>143.01608)
     {
      v=MathPow(x,0.5*x-0.25);
      y=v*(v/y);
     }
   else
      y=MathPow(x,x-0.5)/y;
//--- return result
   return(2.50662827463100050242*y*w);
  }
//+------------------------------------------------------------------+
//| Gamma function                                                   |
//| Input parameters:                                                |
//|     X   -   argument                                             |
//| Domain:                                                          |
//|     0 < X < 171.6                                                |
//|     -170 < X < 0, X is not an integer.                           |
//| Relative error:                                                  |
//|  arithmetic   domain     # trials      peak         rms          |
//|     IEEE    -170,-33      20000       2.3e-15     3.3e-16        |
//|     IEEE     -33,  33     20000       9.4e-16     2.2e-16        |
//|     IEEE      33, 171.6   20000       2.3e-15     3.2e-16        |
//+------------------------------------------------------------------+
static double CGammaFunc::GammaFunc(double x)
  {
//--- create variables
   double p;
   double pp;
   double q=MathAbs(x);
   double qq;
   double z;
   int    i;
   double sgngam=1;
//--- check
   if(q>33.0)
     {
      //--- check
      if(x<0.0)
        {
         p=(int)MathFloor(q);
         i=(int)MathRound(p);
         //--- check
         if(i%2==0)
            sgngam=-1;
         z=q-p;
         //--- check
         if(z>0.5)
           {
            p++;
            z=q-p;
           }
         //--- calculation z
         z=q*MathSin(M_PI*z);
         z=MathAbs(z);
         z=M_PI/(z*GammaStirlFunc(q));
        }
      else
         z=GammaStirlFunc(x);
      //--- return result
      return(sgngam*z);
     }
//--- q<=33
   z=1;
   while(x>=3)
     {
      x--;
      z*=x;
     }
//---  x<0
   while(x<0)
     {
      //--- check
      if(x>-0.000000001)
         return(z/((1+0.5772156649015329*x)*x));
      z/=x;
      x++;
     }
//--- 0<x<2
   while(x<2)
     {
      //--- check
      if(x<0.000000001)
         return(z/((1+0.5772156649015329*x)*x));
      z=z/x;
      x=x+1.0;
     }
//--- check
   if(x==2)
      return(z);
//--- calculation
   x-=2.0;
   pp=1.60119522476751861407E-4;
   pp=1.19135147006586384913E-3+x*pp;
   pp=1.04213797561761569935E-2+x*pp;
   pp=4.76367800457137231464E-2+x*pp;
   pp=2.07448227648435975150E-1+x*pp;
   pp=4.94214826801497100753E-1+x*pp;
   pp=9.99999999999999996796E-1+x*pp;
   qq=-2.31581873324120129819E-5;
   qq=5.39605580493303397842E-4+x*qq;
   qq=-4.45641913851797240494E-3+x*qq;
   qq=1.18139785222060435552E-2+x*qq;
   qq=3.58236398605498653373E-2+x*qq;
   qq=-2.34591795718243348568E-1+x*qq;
   qq=7.14304917030273074085E-2+x*qq;
   qq=1.00000000000000000320+x*qq;
//--- return result
   return(z*pp/qq);
  }
//+------------------------------------------------------------------+
//| Natural logarithm of gamma function                              |
//| Input parameters:                                                |
//|     X       -   argument                                         |
//| Result:                                                          |
//|     logarithm of the absolute value of the Gamma(X).             |
//| Output parameters:                                               |
//|     SgnGam  -   sign(Gamma(X))                                   |
//| Domain:                                                          |
//|     0 < X < 2.55e305                                             |
//|     -2.55e305 < X < 0, X is not an integer.                      |
//| ACCURACY:                                                        |
//| arithmetic      domain        # trials     peak         rms      |
//|    IEEE    0, 3                 28000     5.4e-16     1.1e-16    |
//|    IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17    |
//| The error criterion was relative when the function magnitude     |
//| was greater than one but absolute when it was less than one.     |
//| The following test used the relative error criterion, though     |
//| at certain points the relative error could be much higher than   |
//| indicated.                                                       |
//|    IEEE    -200, -4             10000     4.8e-16     1.3e-16    |
//+------------------------------------------------------------------+
static double CGammaFunc::LnGamma(double x,double &sgngam)
  {
//--- create variables
   double a;
   double b;
   double c;
   double p;
   double q;
   double u;
   double w;
   double z;
   int    i;
   double logpi;
   double ls2pi;
   double tmp=0;
//--- initialization
   sgngam=1;
   logpi=1.14472988584940017414;
   ls2pi=0.91893853320467274178;
//--- check
   if(x<-34.0)
     {
      q=-x;
      w=LnGamma(q,tmp);
      p=(int)MathFloor(q);
      i=(int)MathRound(p);
      //--- check
      if(i%2==0)
         sgngam=-1;
      else
         sgngam=1;
      z=q-p;
      //--- check
      if(z>0.5)
        {
         p++;
         z=p-q;
        }
      z=q*MathSin(M_PI*z);
      //--- return result
      return(logpi-MathLog(z)-w);
     }
//--- x>=-34
   if(x<13)
     {
      z=1;
      p=0;
      u=x;
      //--- cycle
      while(u>=3)
        {
         p--;
         u=x+p;
         z*=u;
        }
      //--- cycle
      while(u<2)
        {
         z/=u;
         p++;
         u=x+p;
        }
      //--- check
      if(z<0)
        {
         sgngam=-1;
         z=-z;
        }
      else
         sgngam=1;
      //--- check
      if(u==2)
         return(MathLog(z));
      //--- calculation
      p-=2;
      x+=p;
      b=-1378.25152569120859100;
      b=-38801.6315134637840924+x*b;
      b=-331612.992738871184744+x*b;
      b=-1162370.97492762307383+x*b;
      b=-1721737.00820839662146+x*b;
      b=-853555.664245765465627+x*b;
      c=1;
      c=-351.815701436523470549+x*c;
      c=-17064.2106651881159223+x*c;
      c=-220528.590553854454839+x*c;
      c=-1139334.44367982507207+x*c;
      c=-2532523.07177582951285+x*c;
      c=-2018891.41433532773231+x*c;
      p=x*b/c;
      //--- return result
      return(MathLog(z)+p);
     }
//--- x>=13
   q=(x-0.5)*MathLog(x)-x+ls2pi;
//--- check
   if(x>100000000)
      return(q);
//--- change value
   p=1/(x*x);
//--- check
   if(x>=1000.0)
      q+=((7.9365079365079365079365*0.0001*p-2.7777777777777777777778*0.001)*p+0.0833333333333333333333)/x;
   else
     {
      a=8.11614167470508450300*0.0001;
      a=-(5.95061904284301438324*0.0001)+p*a;
      a=7.93650340457716943945*0.0001+p*a;
      a=-(2.77777777730099687205*0.001)+p*a;
      a=8.33333333333331927722*0.01+p*a;
      q+=a/x;
     }
//--- return result
   return(q);
  }
//+------------------------------------------------------------------+
//| Normal distribution                                              |
//+------------------------------------------------------------------+
class CNormalDistr
  {
public:
   //--- constructor, destructor
                     CNormalDistr(void);
                    ~CNormalDistr(void);
   //--- methods
   static double     ErrorFunction(double x);
   static double     ErrorFunctionC(double x);
   static double     NormalDistribution(const double x);
   static double     InvErF(const double e);
   static double     InvNormalDistribution(double y0);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNormalDistr::CNormalDistr(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNormalDistr::~CNormalDistr(void)
  {

  }
//+------------------------------------------------------------------+
//| Error function                                                   |
//| The integral is                                                  |
//|                           x                                      |
//|                            -                                     |
//|                 2         | |          2                         |
//|   erf(x)  =  --------     |    exp( - t  ) dt.                   |
//|              sqrt(pi)   | |                                      |
//|                          -                                       |
//|                           0                                      |
//| For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise      |
//| erf(x) = 1 - erfc(x).                                            |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0,1         30000       3.7e-16     1.0e-16         |
//+------------------------------------------------------------------+
static double CNormalDistr::ErrorFunction(double x)
  {
//--- create variables
   double xsq=0;
   double s=0;
   double p=0;
   double q=0;
//--- initialization
   s=MathSign(x);
//--- change value
   x=MathAbs(x);
//--- check
   if(x<0.5)
     {
      //--- calculation
      xsq=x*x;
      p=0.007547728033418631287834;
      p=0.288805137207594084924010+xsq*p;
      p=14.3383842191748205576712+xsq*p;
      p=38.0140318123903008244444+xsq*p;
      p=3017.82788536507577809226+xsq*p;
      p=7404.07142710151470082064+xsq*p;
      p=80437.3630960840172832162+xsq*p;
      q=0.0;
      q=1.00000000000000000000000+xsq*q;
      q=38.0190713951939403753468+xsq*q;
      q=658.070155459240506326937+xsq*q;
      q=6379.60017324428279487120+xsq*q;
      q=34216.5257924628539769006+xsq*q;
      q=80437.3630960840172826266+xsq*q;
      //--- return result
      return(s*1.1283791670955125738961589031*x*p/q);
     }
//--- check
   if(x>=10)
      return(s);
//--- return result
   return(s*(1-ErrorFunctionC(x)));
  }
//+------------------------------------------------------------------+
//| Complementary error function                                     |
//|  1 - erf(x) =                                                    |
//|                           inf.                                   |
//|                             -                                    |
//|                  2         | |          2                        |
//|   erfc(x)  =  --------     |    exp( - t  ) dt                   |
//|               sqrt(pi)   | |                                     |
//|                           -                                      |
//|                            x                                     |
//| For small x, erfc(x) = 1 - erf(x); otherwise rational            |
//| approximations are computed.                                     |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0,26.6417   30000       5.7e-14     1.5e-14         |
//+------------------------------------------------------------------+
static double CNormalDistr::ErrorFunctionC(double x)
  {
//--- create variables
   double p=0;
   double q=0;
//--- check
   if(x<0.0)
      return(2-ErrorFunctionC(-x));
//--- check
   if(x<0.5)
      return(1.0-ErrorFunction(x));
//--- check
   if(x>=10)
      return(0);
//--- calculation
   p=0.0;
   p=0.5641877825507397413087057563+x*p;
   p=9.675807882987265400604202961+x*p;
   p=77.08161730368428609781633646+x*p;
   p=368.5196154710010637133875746+x*p;
   p=1143.262070703886173606073338+x*p;
   p=2320.439590251635247384768711+x*p;
   p=2898.0293292167655611275846+x*p;
   p=1826.3348842295112592168999+x*p;
   q=1.0;
   q=17.14980943627607849376131193+x*q;
   q=137.1255960500622202878443578+x*q;
   q=661.7361207107653469211984771+x*q;
   q=2094.384367789539593790281779+x*q;
   q=4429.612803883682726711528526+x*q;
   q=6089.5424232724435504633068+x*q;
   q=4958.82756472114071495438422+x*q;
   q=1826.3348842295112595576438+x*q;
//--- return result
   return(MathExp(-CMath::Sqr(x))*p/q);
  }
//+------------------------------------------------------------------+
//| Normal distribution function                                     |
//| Returns the area under the Gaussian probability density          |
//| function, integrated from minus infinity to x:                   |
//|                            x                                     |
//|                             -                                    |
//|                   1        | |          2                        |
//|    ndtr(x)  = ---------    |    exp( - t /2 ) dt                 |
//|               sqrt(2pi)  | |                                     |
//|                           -                                      |
//|                          -inf.                                   |
//|             =  ( 1 + erf(z) ) / 2                                |
//|             =  erfc(z) / 2                                       |
//| where z = x/sqrt(2). Computation is via the functions            |
//| erf and erfc.                                                    |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE     -13,0        30000       3.4e-14     6.7e-15         |
//+------------------------------------------------------------------+
static double CNormalDistr::NormalDistribution(const double x)
  {
//--- return result
   return(0.5*(ErrorFunction(x/1.41421356237309504880)+1));
  }
//+------------------------------------------------------------------+
//| Inverse of the error function                                    |
//+------------------------------------------------------------------+
static double CNormalDistr::InvErF(const double e)
  {
//--- return result
   return(InvNormalDistribution(0.5*(e+1))/MathSqrt(2));
  }
//+------------------------------------------------------------------+
//| Inverse of Normal distribution function                          |
//| Returns the argument, x, for which the area under the            |
//| Gaussian probability density function (integrated from           |
//| minus infinity to x) is equal to y.                              |
//| For small arguments 0 < y < exp(-2), the program computes        |
//| z = sqrt( -2.0 * log(y) );  then the approximation is            |
//| x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).                       |
//| There are two rational functions P/Q, one for 0 < y < exp(-32)   |
//| and the other for y up to exp(-2).  For larger arguments,        |
//| w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).       |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain        # trials      peak         rms        |
//|    IEEE     0.125, 1        20000       7.2e-16     1.3e-16      |
//|    IEEE     3e-308, 0.135   50000       4.6e-16     9.8e-17      |
//+------------------------------------------------------------------+
static double CNormalDistr::InvNormalDistribution(double y0)
  {
//--- check
   if(y0<=0)
      return(-CMath::m_maxrealnumber);
//--- check
   if(y0>=1)
      return(CMath::m_maxrealnumber);
//--- create variables
   double expm2=0.13533528323661269189;
   double s2pi=2.50662827463100050242;
   double x;
   double y=y0;
   double z;
   double y2;
   double x0;
   double x1;
   int    code=1;
   double p0;
   double q0;
   double p1;
   double q1;
   double p2;
   double q2;
//--- check
   if(y>1.0-expm2)
     {
      y=1.0-y;
      code=0;
     }
//--- check
   if(y>expm2)
     {
      y-=0.5;
      y2=y*y;
      p0=-59.9633501014107895267;
      p0=98.0010754185999661536+y2*p0;
      p0=-56.6762857469070293439+y2*p0;
      p0=13.9312609387279679503+y2*p0;
      p0=-1.23916583867381258016+y2*p0;
      q0=1;
      q0=1.95448858338141759834+y2*q0;
      q0=4.67627912898881538453+y2*q0;
      q0=86.3602421390890590575+y2*q0;
      q0=-225.462687854119370527+y2*q0;
      q0=200.260212380060660359+y2*q0;
      q0=-82.0372256168333339912+y2*q0;
      q0=15.9056225126211695515+y2*q0;
      q0=-1.18331621121330003142+y2*q0;
      x=y+y*y2*p0/q0;
      x*=s2pi;
      //--- return result
      return(x);
     }
//--- calculation
   x=MathSqrt(-(2.0*MathLog(y)));
   x0=x-MathLog(x)/x;
   z=1.0/x;
//--- check
   if(x<8.0)
     {
      p1=4.05544892305962419923;
      p1=31.5251094599893866154+z*p1;
      p1=57.1628192246421288162+z*p1;
      p1=44.0805073893200834700+z*p1;
      p1=14.6849561928858024014+z*p1;
      p1=2.18663306850790267539+z*p1;
      p1=-(1.40256079171354495875*0.1)+z*p1;
      p1=-(3.50424626827848203418*0.01)+z*p1;
      p1=-(8.57456785154685413611*0.0001)+z*p1;
      q1=1;
      q1=15.7799883256466749731+z*q1;
      q1=45.3907635128879210584+z*q1;
      q1=41.3172038254672030440+z*q1;
      q1=15.0425385692907503408+z*q1;
      q1=2.50464946208309415979+z*q1;
      q1=-(1.42182922854787788574*0.1)+z*q1;
      q1=-(3.80806407691578277194*0.01)+z*q1;
      q1=-(9.33259480895457427372*0.0001)+z*q1;
      x1=z*p1/q1;
     }
   else
     {
      p2=3.23774891776946035970;
      p2=6.91522889068984211695+z*p2;
      p2=3.93881025292474443415+z*p2;
      p2=1.33303460815807542389+z*p2;
      p2=2.01485389549179081538*0.1+z*p2;
      p2=1.23716634817820021358*0.01+z*p2;
      p2=3.01581553508235416007*0.0001+z*p2;
      p2=2.65806974686737550832*0.000001+z*p2;
      p2=6.23974539184983293730*0.000000001+z*p2;
      q2=1;
      q2=6.02427039364742014255+z*q2;
      q2=3.67983563856160859403+z*q2;
      q2=1.37702099489081330271+z*q2;
      q2=2.16236993594496635890*0.1+z*q2;
      q2=1.34204006088543189037*0.01+z*q2;
      q2=3.28014464682127739104*0.0001+z*q2;
      q2=2.89247864745380683936*0.000001+z*q2;
      q2=6.79019408009981274425*0.000000001+z*q2;
      x1=z*p2/q2;
     }
   x=x0-x1;
//--- check
   if(code!=0)
      x=-x;
//--- return result
   return(x);
  }
//+------------------------------------------------------------------+
//| Incomplete gamma function                                        |
//+------------------------------------------------------------------+
class CIncGammaF
  {
public:
   //--- constructor, destructor
                     CIncGammaF(void);
                    ~CIncGammaF(void);
   //--- methods
   static double     IncompleteGamma(const double a,const double x);
   static double     IncompleteGammaC(const double a,const double x);
   static double     InvIncompleteGammaC(const double a,const double y0);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CIncGammaF::CIncGammaF(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CIncGammaF::~CIncGammaF(void)
  {

  }
//+------------------------------------------------------------------+
//| Incomplete gamma integral                                        |
//| The function is defined by                                       |
//|                           x                                      |
//|                            -                                     |
//|                   1       | |  -t  a-1                           |
//|  igam(a,x)  =   -----     |   e   t   dt.                        |
//|                  -      | |                                      |
//|                 | (a)    -                                       |
//|                           0                                      |
//| In this implementation both arguments must be positive.          |
//| The integral is evaluated by either a power series or            |
//| continued fraction expansion, depending on the relative          |
//| values of a and x.                                               |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0,30       200000       3.6e-14     2.9e-15         |
//|    IEEE      0,100      300000       9.9e-14     1.5e-14         |
//+------------------------------------------------------------------+
static double CIncGammaF::IncompleteGamma(const double a,const double x)
  {
//--- check
   if(x<=0 || a<=0)
      return(0);
//--- check
   if(x>1 && x>a)
      return(1-IncompleteGammaC(a,x));
//--- create variables
   double igammaepsilon=0.000000000000001;
   double ans;
   double ax;
   double c;
   double r;
   double tmp=0;
//--- change value
   ax=a*MathLog(x)-x-CGammaFunc::LnGamma(a,tmp);
//--- check
   if(ax<-709.78271289338399)
      return(0);
//--- change values
   ax=MathExp(ax);
   r=a;
   c=1;
   ans=1;
//--- cycle
   do
     {
      r++;
      c*=x/r;
      ans+=c;
     }
   while(c/ans>igammaepsilon);
//--- return result
   return(ans*ax/a);
  }
//+------------------------------------------------------------------+
//| Complemented incomplete gamma integral                           |
//| The function is defined by                                       |
//|  igamc(a,x)   =   1 - igam(a,x)                                  |
//|                            inf.                                  |
//|                              -                                   |
//|                     1       | |  -t  a-1                         |
//|               =   -----     |   e   t   dt.                      |
//|                    -      | |                                    |
//|                   | (a)    -                                     |
//|                             x                                    |
//| In this implementation both arguments must be positive.          |
//| The integral is evaluated by either a power series or            |
//| continued fraction expansion, depending on the relative          |
//| values of a and x.                                               |
//| ACCURACY:                                                        |
//| Tested at random a, x.                                           |
//|                a         x                      Relative error:  |
//| arithmetic   domain   domain     # trials      peak         rms  |
//|    IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15|
//|    IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15|
//+------------------------------------------------------------------+
static double CIncGammaF::IncompleteGammaC(const double a,const double x)
  {
//--- check
   if(x<=0 || a<=0)
      return(1);
//--- check
   if(x<1 || x<a)
      return(1-IncompleteGamma(a,x));
//--- create variables
   double igammaepsilon=0.000000000000001;
   double igammabignumber=4503599627370496.0;
   double igammabignumberinv=2.22044604925031308085*0.0000000000000001;
   double ans=0;
   double ax=0;
   double c=0;
   double yc=0;
   double r=0;
   double t=0;
   double y=0;
   double z=0;
   double pk=0;
   double pkm1=0;
   double pkm2=0;
   double qk=0;
   double qkm1=0;
   double qkm2=0;
   double tmp=0;
//--- change value
   ax=a*MathLog(x)-x-CGammaFunc::LnGamma(a,tmp);
//--- check
   if(ax<-709.78271289338399)
      return(0);
//--- change values
   ax=MathExp(ax);
   y=1-a;
   z=x+y+1;
   c=0;
   pkm2=1;
   qkm2=x;
   pkm1=x+1;
   qkm1=z*x;
   ans=pkm1/qkm1;
//--- cycle
   do
     {
      c++;
      y++;
      z+=2;
      yc=y*c;
      pk=pkm1*z-pkm2*yc;
      qk=qkm1*z-qkm2*yc;
      //--- check
      if(qk!=0)
        {
         r=pk/qk;
         t=MathAbs((ans-r)/r);
         ans=r;
        }
      else
         t=1;
      //--- change values
      pkm2=pkm1;
      pkm1=pk;
      qkm2=qkm1;
      qkm1=qk;
      //--- check
      if(MathAbs(pk)>igammabignumber)
        {
         pkm2=pkm2*igammabignumberinv;
         pkm1=pkm1*igammabignumberinv;
         qkm2=qkm2*igammabignumberinv;
         qkm1=qkm1*igammabignumberinv;
        }
     }
   while(t>igammaepsilon);
//--- return result
   return(ans*ax);
  }
//+------------------------------------------------------------------+
//| Inverse of complemented imcomplete gamma integral                |
//| Given p, the function finds x such that                          |
//|  igamc( a, x ) = p.                                              |
//| Starting with the approximate value                              |
//|         3                                                        |
//|  x = a t                                                         |
//|  where                                                           |
//|  t = 1 - d - ndtri(p) sqrt(d)                                    |
//| and                                                              |
//|  d = 1/9a,                                                       |
//| the routine performs up to 10 Newton iterations to find the      |
//| root of igamc(a,x) - p = 0.                                      |
//| ACCURACY:                                                        |
//| Tested at random a, p in the intervals indicated.                |
//|                a        p                      Relative error:   |
//| arithmetic   domain   domain     # trials      peak        rms   |
//|    IEEE     0.5,100   0,0.5       100000       1.0e-14    1.7e-15|
//|    IEEE     0.01,0.5  0,0.5       100000       9.0e-14    3.4e-15|
//|    IEEE    0.5,10000  0,0.5        20000       2.3e-13    3.8e-14|
//+------------------------------------------------------------------+
static double CIncGammaF::InvIncompleteGammaC(const double a,const double y0)
  {
//--- create variables
   double igammaepsilon=0;
   double iinvgammabignumber=0;
   double x0=0;
   double x1=0;
   double x=0;
   double yl=0;
   double yh=0;
   double y=0;
   double d=0;
   double lgm=0;
   double dithresh=0;
   int    i=0;
   int    dir=0;
   double tmp=0;
//--- initialization
   igammaepsilon=0.000000000000001;
   iinvgammabignumber=4503599627370496.0;
   x0=iinvgammabignumber;
   yl=0;
   x1=0;
   yh=1;
   dithresh=5*igammaepsilon;
   d=1/(9*a);
   y=1-d-CNormalDistr::InvNormalDistribution(y0)*MathSqrt(d);
   x=a*y*y*y;
   lgm=CGammaFunc::LnGamma(a,tmp);
   i=0;
//--- cycle
   while(i<10)
     {
      //--- check
      if(x>x0 || x<x1)
        {
         d=0.0625;
         break;
        }
      //--- function call
      y=IncompleteGammaC(a,x);
      //--- check
      if(y<yl || y>yh)
        {
         d=0.0625;
         break;
        }
      //--- check
      if(y<y0)
        {
         x0=x;
         yl=y;
        }
      else
        {
         x1=x;
         yh=y;
        }
      //--- change value
      d=(a-1)*MathLog(x)-x-lgm;
      //--- check
      if(d<(double)(-709.78271289338399))
        {
         d=0.0625;
         break;
        }
      //--- change values
      d=-MathExp(d);
      d=(y-y0)/d;
      //--- check
      if(MathAbs(d/x)<igammaepsilon)
         return(x);
      //--- change values
      x=x-d;
      i=i+1;
     }
//--- check
   if(x0==iinvgammabignumber)
     {
      //--- check
      if(x<=0.0)
         x=1;
      //--- cycle
      while(x0==iinvgammabignumber)
        {
         x=(1+d)*x;
         y=IncompleteGammaC(a,x);
         //--- check
         if(y<y0)
           {
            x0=x;
            yl=y;
            break;
           }
         d=d+d;
        }
     }
//--- change values
   d=0.5;
   dir=0;
   i=0;
//--- cycle
   while(i<400)
     {
      x=x1+d*(x0-x1);
      y=IncompleteGammaC(a,x);
      lgm=(x0-x1)/(x1+x0);
      //--- check
      if(MathAbs(lgm)<dithresh)
         break;
      //--- change value
      lgm=(y-y0)/y0;
      //--- check
      if(MathAbs(lgm)<dithresh)
         break;
      //--- check
      if(x<=0.0)
         break;
      //--- check
      if(y>=y0)
        {
         x1=x;
         yh=y;
         //--- check
         if(dir<0)
           {
            dir=0;
            d=0.5;
           }
         else
           {
            //--- check
            if(dir>1)
               d=0.5*d+0.5;
            else
               d=(y0-yl)/(yh-yl);
           }
         dir=dir+1;
        }
      else
        {
         x0=x;
         yl=y;
         //--- check
         if(dir>0)
           {
            dir=0;
            d=0.5;
           }
         else
           {
            //--- check
            if(dir<-1)
               d=0.5*d;
            else
               d=(y0-yl)/(yh-yl);
           }
         dir=dir-1;
        }
      i=i+1;
     }
//--- return result
   return(x);
  }
//+------------------------------------------------------------------+
//| Airy function                                                    |
//+------------------------------------------------------------------+
class CAiryF
  {
public:
   //--- constructor, destructor
                     CAiryF(void);
                    ~CAiryF(void);
   //--- method
   static void       Airy(const double x,double &ai,double &aip,double &bi,double &bip);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CAiryF::CAiryF(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CAiryF::~CAiryF(void)
  {

  }
//+------------------------------------------------------------------+
//| Airy function                                                    |
//| Solution of the differential equation                            |
//| y"(x) = xy.                                                      |
//| The function returns the two independent solutions Ai, Bi        |
//| and their first derivatives Ai'(x), Bi'(x).                      |
//| Evaluation is by power series summation for small x,             |
//| by rational minimax approximations for large x.                  |
//| ACCURACY:                                                        |
//| Error criterion is absolute when function <= 1, relative         |
//| when function > 1, except * denotes relative error criterion.    |
//| For large negative x, the absolute error increases as x^1.5.     |
//| For large positive x, the relative error increases as x^1.5.     |
//| Arithmetic  domain   function  # trials      peak         rms    |
//| IEEE        -10, 0     Ai        10000       1.6e-15     2.7e-16 |
//| IEEE          0, 10    Ai        10000       2.3e-14*    1.8e-15*|
//| IEEE        -10, 0     Ai'       10000       4.6e-15     7.6e-16 |
//| IEEE          0, 10    Ai'       10000       1.8e-14*    1.5e-15*|
//| IEEE        -10, 10    Bi        30000       4.2e-15     5.3e-16 |
//| IEEE        -10, 10    Bi'       30000       4.9e-15     7.3e-16 |
//+------------------------------------------------------------------+
static void CAiryF::Airy(const double x,double &ai,double &aip,double &bi,
                         double &bip)
  {
//--- create variables
   double z=0;
   double zz=0;
   double t=0;
   double f=0;
   double g=0;
   double uf=0;
   double ug=0;
   double k=0;
   double zeta=0;
   double theta=0;
   int    domflg=0;
   double c1=0;
   double c2=0;
   double sqrt3=0;
   double sqpii=0;
   double afn=0;
   double afd=0;
   double agn=0;
   double agd=0;
   double apfn=0;
   double apfd=0;
   double apgn=0;
   double apgd=0;
   double an=0;
   double ad=0;
   double apn=0;
   double apd=0;
   double bn16=0;
   double bd16=0;
   double bppn=0;
   double bppd=0;
//--- initialization
   ai=0;
   aip=0;
   bi=0;
   bip=0;
   sqpii=5.64189583547756286948E-1;
   c1=0.35502805388781723926;
   c2=0.258819403792806798405;
   sqrt3=1.732050807568877293527;
   domflg=0;
//--- check
   if(x>25.77)
     {
      //--- initialization
      ai=0;
      aip=0;
      bi=CMath::m_maxrealnumber;
      bip=CMath::m_maxrealnumber;
      //--- exit the function
      return;
     }
//--- check
   if(x<-2.09)
     {
      //--- calculation
      domflg=15;
      t=MathSqrt(-x);
      zeta=-(2.0*x*t/3.0);
      t=MathSqrt(t);
      k=sqpii/t;
      z=1.0/zeta;
      zz=z*z;
      //--- calculation
      afn=-1.31696323418331795333E-1;
      afn=afn*zz-6.26456544431912369773E-1;
      afn=afn*zz-6.93158036036933542233E-1;
      afn=afn*zz-2.79779981545119124951E-1;
      afn=afn*zz-4.91900132609500318020E-2;
      afn=afn*zz-4.06265923594885404393E-3;
      afn=afn*zz-1.59276496239262096340E-4;
      afn=afn*zz-2.77649108155232920844E-6;
      afn=afn*zz-1.67787698489114633780E-8;
      //--- calculation
      afd=1.00000000000000000000E0;
      afd=afd*zz+1.33560420706553243746E1;
      afd=afd*zz+3.26825032795224613948E1;
      afd=afd*zz+2.67367040941499554804E1;
      afd=afd*zz+9.18707402907259625840E0;
      afd=afd*zz+1.47529146771666414581E0;
      afd=afd*zz+1.15687173795188044134E-1;
      afd=afd*zz+4.40291641615211203805E-3;
      afd=afd*zz+7.54720348287414296618E-5;
      afd=afd*zz+4.51850092970580378464E-7;
      uf=1.0+zz*afn/afd;
      //--- calculation
      agn=1.97339932091685679179E-2;
      agn=agn*zz+3.91103029615688277255E-1;
      agn=agn*zz+1.06579897599595591108E0;
      agn=agn*zz+9.39169229816650230044E-1;
      agn=agn*zz+3.51465656105547619242E-1;
      agn=agn*zz+6.33888919628925490927E-2;
      agn=agn*zz+5.85804113048388458567E-3;
      agn=agn*zz+2.82851600836737019778E-4;
      agn=agn*zz+6.98793669997260967291E-6;
      agn=agn*zz+8.11789239554389293311E-8;
      agn=agn*zz+3.41551784765923618484E-10;
      //--- calculation
      agd=1.00000000000000000000E0;
      agd=agd*zz+9.30892908077441974853E0;
      agd=agd*zz+1.98352928718312140417E1;
      agd=agd*zz+1.55646628932864612953E1;
      agd=agd*zz+5.47686069422975497931E0;
      agd=agd*zz+9.54293611618961883998E-1;
      agd=agd*zz+8.64580826352392193095E-2;
      agd=agd*zz+4.12656523824222607191E-3;
      agd=agd*zz+1.01259085116509135510E-4;
      agd=agd*zz+1.17166733214413521882E-6;
      agd=agd*zz+4.91834570062930015649E-9;
      //--- calculation
      ug=z*agn/agd;
      theta=zeta+0.25*M_PI;
      f=MathSin(theta);
      g=MathCos(theta);
      ai=k*(f*uf-g*ug);
      bi=k*(g*uf+f*ug);
      apfn=1.85365624022535566142E-1;
      apfn=apfn*zz+8.86712188052584095637E-1;
      apfn=apfn*zz+9.87391981747398547272E-1;
      apfn=apfn*zz+4.01241082318003734092E-1;
      apfn=apfn*zz+7.10304926289631174579E-2;
      apfn=apfn*zz+5.90618657995661810071E-3;
      apfn=apfn*zz+2.33051409401776799569E-4;
      apfn=apfn*zz+4.08718778289035454598E-6;
      apfn=apfn*zz+2.48379932900442457853E-8;
      //--- calculation
      apfd=1.00000000000000000000E0;
      apfd=apfd*zz+1.47345854687502542552E1;
      apfd=apfd*zz+3.75423933435489594466E1;
      apfd=apfd*zz+3.14657751203046424330E1;
      apfd=apfd*zz+1.09969125207298778536E1;
      apfd=apfd*zz+1.78885054766999417817E0;
      apfd=apfd*zz+1.41733275753662636873E-1;
      apfd=apfd*zz+5.44066067017226003627E-3;
      apfd=apfd*zz+9.39421290654511171663E-5;
      apfd=apfd*zz+5.65978713036027009243E-7;
      uf=1.0+zz*apfn/apfd;
      //--- calculation
      apgn=-3.55615429033082288335E-2;
      apgn=apgn*zz-6.37311518129435504426E-1;
      apgn=apgn*zz-1.70856738884312371053E0;
      apgn=apgn*zz-1.50221872117316635393E0;
      apgn=apgn*zz-5.63606665822102676611E-1;
      apgn=apgn*zz-1.02101031120216891789E-1;
      apgn=apgn*zz-9.48396695961445269093E-3;
      apgn=apgn*zz-4.60325307486780994357E-4;
      apgn=apgn*zz-1.14300836484517375919E-5;
      apgn=apgn*zz-1.33415518685547420648E-7;
      apgn=apgn*zz-5.63803833958893494476E-10;
      //--- calculation
      apgd=1.00000000000000000000E0;
      apgd=apgd*zz+9.85865801696130355144E0;
      apgd=apgd*zz+2.16401867356585941885E1;
      apgd=apgd*zz+1.73130776389749389525E1;
      apgd=apgd*zz+6.17872175280828766327E0;
      apgd=apgd*zz+1.08848694396321495475E0;
      apgd=apgd*zz+9.95005543440888479402E-2;
      apgd=apgd*zz+4.78468199683886610842E-3;
      apgd=apgd*zz+1.18159633322838625562E-4;
      apgd=apgd*zz+1.37480673554219441465E-6;
      apgd=apgd*zz+5.79912514929147598821E-9;
      ug=z*apgn/apgd;
      k=sqpii*t;
      aip=-(k*(g*uf+f*ug));
      bip=k*(f*uf-g*ug);
      //--- exit the function
      return;
     }
//--- check
   if(x>=(double)(2.09))
     {
      domflg=5;
      t=MathSqrt(x);
      zeta=2.0*x*t/3.0;
      g=MathExp(zeta);
      t=MathSqrt(t);
      k=2.0*t*g;
      z=1.0/zeta;
      //--- calculation
      an=3.46538101525629032477E-1;
      an=an*z+1.20075952739645805542E1;
      an=an*z+7.62796053615234516538E1;
      an=an*z+1.68089224934630576269E2;
      an=an*z+1.59756391350164413639E2;
      an=an*z+7.05360906840444183113E1;
      an=an*z+1.40264691163389668864E1;
      an=an*z+9.99999999999999995305E-1;
      ad=5.67594532638770212846E-1;
      ad=ad*z+1.47562562584847203173E1;
      ad=ad*z+8.45138970141474626562E1;
      ad=ad*z+1.77318088145400459522E2;
      ad=ad*z+1.64234692871529701831E2;
      ad=ad*z+7.14778400825575695274E1;
      ad=ad*z+1.40959135607834029598E1;
      ad=ad*z+1.00000000000000000470E0;
      //--- calculation
      f=an/ad;
      ai=sqpii*f/k;
      k=-(0.5*sqpii*t/g);
      apn=6.13759184814035759225E-1;
      apn=apn*z+1.47454670787755323881E1;
      apn=apn*z+8.20584123476060982430E1;
      apn=apn*z+1.71184781360976385540E2;
      apn=apn*z+1.59317847137141783523E2;
      apn=apn*z+6.99778599330103016170E1;
      apn=apn*z+1.39470856980481566958E1;
      apn=apn*z+1.00000000000000000550E0;
      apd=3.34203677749736953049E-1;
      apd=apd*z+1.11810297306158156705E1;
      apd=apd*z+7.11727352147859965283E1;
      apd=apd*z+1.58778084372838313640E2;
      apd=apd*z+1.53206427475809220834E2;
      apd=apd*z+6.86752304592780337944E1;
      apd=apd*z+1.38498634758259442477E1;
      apd=apd*z+9.99999999999999994502E-1;
      f=apn/apd;
      aip=f*k;
      //--- check
      if(x>(double)(8.3203353))
        {
         //--- calculation
         bn16=-2.53240795869364152689E-1;
         bn16=bn16*z+5.75285167332467384228E-1;
         bn16=bn16*z-3.29907036873225371650E-1;
         bn16=bn16*z+6.44404068948199951727E-2;
         bn16=bn16*z-3.82519546641336734394E-3;
         bd16=1.00000000000000000000E0;
         bd16=bd16*z-7.15685095054035237902E0;
         bd16=bd16*z+1.06039580715664694291E1;
         bd16=bd16*z-5.23246636471251500874E0;
         bd16=bd16*z+9.57395864378383833152E-1;
         bd16=bd16*z-5.50828147163549611107E-2;
         //--- calculation
         f=z*bn16/bd16;
         k=sqpii*g;
         bi=k*(1.0+f)/t;
         bppn=4.65461162774651610328E-1;
         bppn=bppn*z-1.08992173800493920734E0;
         bppn=bppn*z+6.38800117371827987759E-1;
         bppn=bppn*z-1.26844349553102907034E-1;
         bppn=bppn*z+7.62487844342109852105E-3;
         bppd=1.00000000000000000000E0;
         bppd=bppd*z-8.70622787633159124240E0;
         bppd=bppd*z+1.38993162704553213172E1;
         bppd=bppd*z-7.14116144616431159572E0;
         bppd=bppd*z+1.34008595960680518666E0;
         bppd=bppd*z-7.84273211323341930448E-2;
         f=z*bppn/bppd;
         bip=k*t*(1.0+f);
         //--- exit the function
         return;
        }
     }
//--- change values
   f=1.0;
   g=x;
   t=1.0;
   uf=1.0;
   ug=x;
   k=1.0;
   z=x*x*x;
//--- cycle
   while(t>CMath::m_machineepsilon)
     {
      //--- calculation
      uf=uf*z;
      k=k+1.0;
      uf=uf/k;
      ug=ug*z;
      k=k+1.0;
      ug=ug/k;
      uf=uf/k;
      f=f+uf;
      k=k+1.0;
      ug=ug/k;
      g=g+ug;
      t=MathAbs(uf/f);
     }
//--- change values
   uf=c1*f;
   ug=c2*g;
//--- check
   if(domflg%2==0)
      ai=uf-ug;
//--- check
   if(domflg/2%2==0)
      bi=sqrt3*(uf+ug);
//--- change values
   k=4.0;
   uf=x*x/2.0;
   ug=z/3.0;
   f=uf;
   g=1.0+ug;
   uf=uf/3.0;
   t=1.0;
//--- cycle
   while(t>CMath::m_machineepsilon)
     {
      //--- calculation
      uf=uf*z;
      ug=ug/k;
      k=k+1.0;
      ug=ug*z;
      uf=uf/k;
      f=f+uf;
      k=k+1.0;
      ug=ug/k;
      uf=uf/k;
      g=g+ug;
      k=k+1.0;
      t=MathAbs(ug/g);
     }
//--- change values
   uf=c1*f;
   ug=c2*g;
//--- check
   if(domflg/4%2==0)
      aip=uf-ug;
//--- check
   if(domflg/8%2==0)
      bip=sqrt3*(uf+ug);
  }
//+------------------------------------------------------------------+
//| Bessel function                                                  |
//+------------------------------------------------------------------+
class CBessel
  {
private:
   //--- private methods
   static void       BesselMFirstCheb(const double c,double &b0,double &b1,double &b2);
   static void       BesselMNextCheb(const double x,const double c,double &b0,double &b1,double &b2);
   static void       BesselM1FirstCheb(const double c,double &b0,double &b1,double &b2);
   static void       BesselM1NextCheb(const double x,const double c,double &b0,double &b1,double &b2);
   static void       BesselAsympt0(const double x,double &pzero,double &qzero);
   static void       BesselAsympt1(const double x,double &pzero,double &qzero);
public:
   //--- constructor, destructor
                     CBessel(void);
                    ~CBessel(void);
   //--- public methods
   static double     BesselJ0(double x);
   static double     BesselJ1(double x);
   static double     BesselJN(int n,double x);
   static double     BesselY0(double x);
   static double     BesselY1(double x);
   static double     BesselYN(int n,double x);
   static double     BesselI0(double x);
   static double     BesselI1(double x);
   static double     BesselK0(double x);
   static double     BesselK1(double x);
   static double     BesselKN(int nn,double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CBessel::CBessel(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CBessel::~CBessel(void)
  {

  }
//+------------------------------------------------------------------+
//| Bessel function of order zero                                    |
//| Returns Bessel function of order zero of the argument.           |
//| The domain is divided into the intervals [0, 5] and              |
//| (5, infinity). In the first interval the following rational      |
//| approximation is used:                                           |
//|        2         2                                               |
//| (w - r  ) (w - r  ) P (w) / Q (w)                                |
//|       1         2    3       8                                   |
//|            2                                                     |
//| where w = x  and the two r's are zeros of the function.          |
//| In the second interval, the Hankel asymptotic expansion          |
//| is employed with two rational functions of degree 6/6            |
//| and 7/7.                                                         |
//| ACCURACY:                                                        |
//|                      Absolute error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0, 30       60000       4.2e-16     1.1e-16         |
//+------------------------------------------------------------------+
static double CBessel::BesselJ0(double x)
  {
//--- create variables
   double xsq=0;
   double nn=0;
   double pzero=0;
   double qzero=0;
   double p1=0;
   double q1=0;
//--- check
   if(x<0.0)
      x=-x;
//--- check
   if(x>8.0)
     {
      //--- function call
      BesselAsympt0(x,pzero,qzero);
      nn=x-M_PI/4;
      //--- return result
      return(MathSqrt(2/M_PI/x)*(pzero*MathCos(nn)-qzero*MathSin(nn)));
     }
//--- calculation
   xsq=CMath::Sqr(x);
   p1=26857.86856980014981415848441;
   p1=-40504123.71833132706360663322+xsq*p1;
   p1=25071582855.36881945555156435+xsq*p1;
   p1=-8085222034853.793871199468171+xsq*p1;
   p1=1434354939140344.111664316553+xsq*p1;
   p1=-136762035308817138.6865416609+xsq*p1;
   p1=6382059341072356562.289432465+xsq*p1;
   p1=-117915762910761053603.8440800+xsq*p1;
   p1=493378725179413356181.6813446+xsq*p1;
   q1=1.0;
   q1=1363.063652328970604442810507+xsq*q1;
   q1=1114636.098462985378182402543+xsq*q1;
   q1=669998767.2982239671814028660+xsq*q1;
   q1=312304311494.1213172572469442+xsq*q1;
   q1=112775673967979.8507056031594+xsq*q1;
   q1=30246356167094626.98627330784+xsq*q1;
   q1=5428918384092285160.200195092+xsq*q1;
   q1=493378725179413356211.3278438+xsq*q1;
//--- return result
   return(p1/q1);
  }
//+------------------------------------------------------------------+
//| Bessel function of order one                                     |
//| Returns Bessel function of order one of the argument.            |
//| The domain is divided into the intervals [0, 8] and              |
//| (8, infinity). In the first interval a 24 term Chebyshev         |
//| expansion is used. In the second, the asymptotic                 |
//| trigonometric representation is employed using two               |
//| rational functions of degree 5/5.                                |
//| ACCURACY:                                                        |
//|                      Absolute error:                             |
//| arithmetic   domain      # trials      peak         rms          |
//|    IEEE      0, 30       30000       2.6e-16     1.1e-16         |
//+------------------------------------------------------------------+
static double CBessel::BesselJ1(double x)
  {
//--- create variables
   double result=0;
   double s=0;
   double xsq=0;
   double nn=0;
   double pzero=0;
   double qzero=0;
   double p1=0;
   double q1=0;
//--- calculation
   s=MathSign(x);
//--- check
   if(x<0.0)
      x=-x;
//--- check
   if(x>8.0)
     {
      //--- function call
      BesselAsympt1(x,pzero,qzero);
      nn=x-3*M_PI/4;
      result=MathSqrt(2/M_PI/x)*(pzero*MathCos(nn)-qzero*MathSin(nn));
      //--- check
      if(s<0.0)
         result=-result;
      //--- return result
      return(result);
     }
//--- calculation
   xsq=CMath::Sqr(x);
   p1=2701.122710892323414856790990;
   p1=-4695753.530642995859767162166+xsq*p1;
   p1=3413234182.301700539091292655+xsq*p1;
   p1=-1322983480332.126453125473247+xsq*p1;
   p1=290879526383477.5409737601689+xsq*p1;
   p1=-35888175699101060.50743641413+xsq*p1;
   p1=2316433580634002297.931815435+xsq*p1;
   p1=-66721065689249162980.20941484+xsq*p1;
   p1=581199354001606143928.050809+xsq*p1;
   q1=1.0;
   q1=1606.931573481487801970916749+xsq*q1;
   q1=1501793.594998585505921097578+xsq*q1;
   q1=1013863514.358673989967045588+xsq*q1;
   q1=524371026216.7649715406728642+xsq*q1;
   q1=208166122130760.7351240184229+xsq*q1;
   q1=60920613989175217.46105196863+xsq*q1;
   q1=11857707121903209998.37113348+xsq*q1;
   q1=1162398708003212287858.529400+xsq*q1;
//--- return result
   return(s*x*p1/q1);
  }
//+------------------------------------------------------------------+
//| Bessel function of integer order                                 |
//| Returns Bessel function of order n, where n is a                 |
//| (possibly negative) integer.                                     |
//| The ratio of jn(x) to j0(x) is computed by backward              |
//| recurrence.  First the ratio jn/jn-1 is found by a               |
//| continued fraction expansion.  Then the recurrence               |
//| relating successive orders is applied until j0 or j1 is          |
//| reached.                                                         |
//| If n = 0 or 1 the routine for j0 or j1 is called                 |
//| directly.                                                        |
//| ACCURACY:                                                        |
//|                      Absolute error:                             |
//| arithmetic   range      # trials      peak         rms           |
//|    IEEE      0, 30        5000       4.4e-16     7.9e-17         |
//| Not suitable for large n or x. Use jv() (fractional order)       |
//| instead.                                                         |
//+------------------------------------------------------------------+
static double CBessel::BesselJN(int n,double x)
  {
//--- create variables
   double result=0;
   double pkm2=0;
   double pkm1=0;
   double pk=0;
   double xk=0;
   double r=0;
   double ans=0;
   int    k=0;
   int    sg=0;
//--- check
   if(n<0)
     {
      n=-n;
      //--- check
      if(n%2==0)
         sg=1;
      else
         sg=-1;
     }
   else
      sg=1;
//--- check
   if(x<0.0)
     {
      //--- check
      if(n%2!=0)
         sg=-sg;
      x=-x;
     }
//--- check
   if(n==0)
      return(sg*BesselJ0(x));
//--- check
   if(n==1)
      return(sg*BesselJ1(x));
//--- check
   if(n==2)
     {
      //--- check
      if(x==0.0)
         return(0);
      else
         return(sg*(2.0*BesselJ1(x)/x-BesselJ0(x)));
     }
//--- check
   if(x<CMath::m_machineepsilon)
      return(0);
//--- change values
   k=53;
   pk=2*(n+k);
   ans=pk;
   xk=x*x;
//--- cycle
   do
     {
      pk=pk-2.0;
      ans=pk-xk/ans;
      k=k-1;
     }
   while(k!=0);
//--- change values
   ans=x/ans;
   pk=1.0;
   pkm1=1.0/ans;
   k=n-1;
   r=2*k;
//--- cycle
   do
     {
      pkm2=(pkm1*r-pk*x)/x;
      pk=pkm1;
      pkm1=pkm2;
      r=r-2.0;
      k=k-1;
     }
   while(k!=0);
//--- check
   if(MathAbs(pk)>MathAbs(pkm1))
      ans=BesselJ1(x)/pk;
   else
      ans=BesselJ0(x)/pkm1;
//--- return result
   return(sg*ans);
  }
//+------------------------------------------------------------------+
//| Bessel function of the second kind, order zero                   |
//| Returns Bessel function of the second kind, of order             |
//| zero, of the argument.                                           |
//| The domain is divided into the intervals [0, 5] and              |
//| (5, infinity). In the first interval a rational approximation    |
//| R(x) is employed to compute                                      |
//|   y0(x)  = R(x)  +   2 * log(x) * j0(x) / PI.                    |
//| Thus a call to j0() is required.                                 |
//| In the second interval, the Hankel asymptotic expansion          |
//| is employed with two rational functions of degree 6/6            |
//| and 7/7.                                                         |
//| ACCURACY:                                                        |
//|  Absolute error, when y0(x) < 1; else relative error:            |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0, 30       30000       1.3e-15     1.6e-16         |
//+------------------------------------------------------------------+
static double CBessel::BesselY0(double x)
  {
//--- create variables
   double nn=0;
   double xsq=0;
   double pzero=0;
   double qzero=0;
   double p4=0;
   double q4=0;
//--- check
   if(x>8.0)
     {
      //--- function call
      BesselAsympt0(x,pzero,qzero);
      nn=x-M_PI/4;
      //--- return result
      return(MathSqrt(2/M_PI/x)*(pzero*MathSin(nn)+qzero*MathCos(nn)));
     }
//--- calculation
   xsq=CMath::Sqr(x);
   p4=-41370.35497933148554125235152;
   p4=59152134.65686889654273830069+xsq*p4;
   p4=-34363712229.79040378171030138+xsq*p4;
   p4=10255208596863.94284509167421+xsq*p4;
   p4=-1648605817185729.473122082537+xsq*p4;
   p4=137562431639934407.8571335453+xsq*p4;
   p4=-5247065581112764941.297350814+xsq*p4;
   p4=65874732757195549259.99402049+xsq*p4;
   p4=-27502866786291095837.01933175+xsq*p4;
   q4=1.0;
   q4=1282.452772478993804176329391+xsq*q4;
   q4=1001702.641288906265666651753+xsq*q4;
   q4=579512264.0700729537480087915+xsq*q4;
   q4=261306575504.1081249568482092+xsq*q4;
   q4=91620380340751.85262489147968+xsq*q4;
   q4=23928830434997818.57439356652+xsq*q4;
   q4=4192417043410839973.904769661+xsq*q4;
   q4=372645883898616588198.9980+xsq*q4;
//--- return result
   return(p4/q4+2/M_PI*BesselJ0(x)*MathLog(x));
  }
//+------------------------------------------------------------------+
//| Bessel function of second kind of order one                      |
//| Returns Bessel function of the second kind of order one          |
//| of the argument.                                                 |
//| The domain is divided into the intervals [0, 8] and              |
//| (8, infinity). In the first interval a 25 term Chebyshev         |
//| expansion is used, and a call to j1() is required.               |
//| In the second, the asymptotic trigonometric representation       |
//| is employed using two rational functions of degree 5/5.          |
//| ACCURACY:                                                        |
//|                      Absolute error:                             |
//| arithmetic   domain      # trials      peak         rms          |
//|    IEEE      0, 30       30000       1.0e-15     1.3e-16         |
//+------------------------------------------------------------------+
static double CBessel::BesselY1(double x)
  {
//--- create variables
   double nn=0;
   double xsq=0;
   double pzero=0;
   double qzero=0;
   double p4=0;
   double q4=0;
//--- check
   if(x>8.0)
     {
      //--- function call
      BesselAsympt1(x,pzero,qzero);
      nn=x-3*M_PI/4;
      //--- return result
      return(MathSqrt(2/M_PI/x)*(pzero*MathSin(nn)+qzero*MathCos(nn)));
     }
//--- calculation
   xsq=CMath::Sqr(x);
   p4=-2108847.540133123652824139923;
   p4=3639488548.124002058278999428+xsq*p4;
   p4=-2580681702194.450950541426399+xsq*p4;
   p4=956993023992168.3481121552788+xsq*p4;
   p4=-196588746272214065.8820322248+xsq*p4;
   p4=21931073399177975921.11427556+xsq*p4;
   p4=-1212297555414509577913.561535+xsq*p4;
   p4=26554738314348543268942.48968+xsq*p4;
   p4=-99637534243069222259967.44354+xsq*p4;
   q4=1.0;
   q4=1612.361029677000859332072312+xsq*q4;
   q4=1563282.754899580604737366452+xsq*q4;
   q4=1128686837.169442121732366891+xsq*q4;
   q4=646534088126.5275571961681500+xsq*q4;
   q4=297663212564727.6729292742282+xsq*q4;
   q4=108225825940881955.2553850180+xsq*q4;
   q4=29549879358971486742.90758119+xsq*q4;
   q4=5435310377188854170800.653097+xsq*q4;
   q4=508206736694124324531442.4152+xsq*q4;
//--- return result
   return(x*p4/q4+2/M_PI*(BesselJ1(x)*MathLog(x)-1/x));
  }
//+------------------------------------------------------------------+
//| Bessel function of second kind of integer order                  |
//| Returns Bessel function of order n, where n is a                 |
//| (possibly negative) integer.                                     |
//| The function is evaluated by forward recurrence on               |
//| n, starting with values computed by the routines                 |
//| y0() and y1().                                                   |
//| If n = 0 or 1 the routine for y0 or y1 is called                 |
//| directly.                                                        |
//| ACCURACY:                                                        |
//|                      Absolute error, except relative             |
//|                      when y > 1:                                 |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0, 30       30000       3.4e-15     4.3e-16         |
//+------------------------------------------------------------------+
static double CBessel::BesselYN(int n,double x)
  {
//--- create variables
   int    i=0;
   double a=0;
   double b=0;
   double tmp=0;
   double s=0;
//--- initialization
   s=1;
//--- check
   if(n<0)
     {
      n=-n;
      //--- check
      if(n%2!=0)
         s=-1;
     }
//--- check
   if(n==0)
      return(BesselY0(x));
//--- check
   if(n==1)
      return(s*BesselY1(x));
//--- calculation
   a=BesselY0(x);
   b=BesselY1(x);
   for(i=1;i<=n-1;i++)
     {
      tmp=b;
      b=2*i/x*b-a;
      a=tmp;
     }
//--- return result
   return(s*b);
  }
//+------------------------------------------------------------------+
//| Modified Bessel function of order zero                           |
//| Returns modified Bessel function of order zero of the            |
//| argument.                                                        |
//| The function is defined as i0(x) = j0( ix ).                     |
//| The range is partitioned into the two intervals [0,8] and        |
//| (8, infinity).  Chebyshev polynomial expansions are employed     |
//| in each interval.                                                |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0,30        30000       5.8e-16     1.4e-16         |
//+------------------------------------------------------------------+
static double CBessel::BesselI0(double x)
  {
//--- create variables
   double y=0;
   double v=0;
   double z=0;
   double b0=0;
   double b1=0;
   double b2=0;
//--- check
   if(x<0.0)
      x=-x;
//--- check
   if(x<=8.0)
     {
      y=x/2.0-2.0;
      //--- function calls
      BesselMFirstCheb(-4.41534164647933937950E-18,b0,b1,b2);
      BesselMNextCheb(y,3.33079451882223809783E-17,b0,b1,b2);
      BesselMNextCheb(y,-2.43127984654795469359E-16,b0,b1,b2);
      BesselMNextCheb(y,1.71539128555513303061E-15,b0,b1,b2);
      BesselMNextCheb(y,-1.16853328779934516808E-14,b0,b1,b2);
      BesselMNextCheb(y,7.67618549860493561688E-14,b0,b1,b2);
      BesselMNextCheb(y,-4.85644678311192946090E-13,b0,b1,b2);
      BesselMNextCheb(y,2.95505266312963983461E-12,b0,b1,b2);
      BesselMNextCheb(y,-1.72682629144155570723E-11,b0,b1,b2);
      BesselMNextCheb(y,9.67580903537323691224E-11,b0,b1,b2);
      BesselMNextCheb(y,-5.18979560163526290666E-10,b0,b1,b2);
      BesselMNextCheb(y,2.65982372468238665035E-9,b0,b1,b2);
      BesselMNextCheb(y,-1.30002500998624804212E-8,b0,b1,b2);
      BesselMNextCheb(y,6.04699502254191894932E-8,b0,b1,b2);
      BesselMNextCheb(y,-2.67079385394061173391E-7,b0,b1,b2);
      BesselMNextCheb(y,1.11738753912010371815E-6,b0,b1,b2);
      BesselMNextCheb(y,-4.41673835845875056359E-6,b0,b1,b2);
      BesselMNextCheb(y,1.64484480707288970893E-5,b0,b1,b2);
      BesselMNextCheb(y,-5.75419501008210370398E-5,b0,b1,b2);
      BesselMNextCheb(y,1.88502885095841655729E-4,b0,b1,b2);
      BesselMNextCheb(y,-5.76375574538582365885E-4,b0,b1,b2);
      BesselMNextCheb(y,1.63947561694133579842E-3,b0,b1,b2);
      BesselMNextCheb(y,-4.32430999505057594430E-3,b0,b1,b2);
      BesselMNextCheb(y,1.05464603945949983183E-2,b0,b1,b2);
      BesselMNextCheb(y,-2.37374148058994688156E-2,b0,b1,b2);
      BesselMNextCheb(y,4.93052842396707084878E-2,b0,b1,b2);
      BesselMNextCheb(y,-9.49010970480476444210E-2,b0,b1,b2);
      BesselMNextCheb(y,1.71620901522208775349E-1,b0,b1,b2);
      BesselMNextCheb(y,-3.04682672343198398683E-1,b0,b1,b2);
      BesselMNextCheb(y,6.76795274409476084995E-1,b0,b1,b2);
      //--- calculation
      v=0.5*(b0-b2);
      //--- return result
      return(MathExp(x)*v);
     }
//--- change value
   z=32.0/x-2.0;
//--- function calls
   BesselMFirstCheb(-7.23318048787475395456E-18,b0,b1,b2);
   BesselMNextCheb(z,-4.83050448594418207126E-18,b0,b1,b2);
   BesselMNextCheb(z,4.46562142029675999901E-17,b0,b1,b2);
   BesselMNextCheb(z,3.46122286769746109310E-17,b0,b1,b2);
   BesselMNextCheb(z,-2.82762398051658348494E-16,b0,b1,b2);
   BesselMNextCheb(z,-3.42548561967721913462E-16,b0,b1,b2);
   BesselMNextCheb(z,1.77256013305652638360E-15,b0,b1,b2);
   BesselMNextCheb(z,3.81168066935262242075E-15,b0,b1,b2);
   BesselMNextCheb(z,-9.55484669882830764870E-15,b0,b1,b2);
   BesselMNextCheb(z,-4.15056934728722208663E-14,b0,b1,b2);
   BesselMNextCheb(z,1.54008621752140982691E-14,b0,b1,b2);
   BesselMNextCheb(z,3.85277838274214270114E-13,b0,b1,b2);
   BesselMNextCheb(z,7.18012445138366623367E-13,b0,b1,b2);
   BesselMNextCheb(z,-1.79417853150680611778E-12,b0,b1,b2);
   BesselMNextCheb(z,-1.32158118404477131188E-11,b0,b1,b2);
   BesselMNextCheb(z,-3.14991652796324136454E-11,b0,b1,b2);
   BesselMNextCheb(z,1.18891471078464383424E-11,b0,b1,b2);
   BesselMNextCheb(z,4.94060238822496958910E-10,b0,b1,b2);
   BesselMNextCheb(z,3.39623202570838634515E-9,b0,b1,b2);
   BesselMNextCheb(z,2.26666899049817806459E-8,b0,b1,b2);
   BesselMNextCheb(z,2.04891858946906374183E-7,b0,b1,b2);
   BesselMNextCheb(z,2.89137052083475648297E-6,b0,b1,b2);
   BesselMNextCheb(z,6.88975834691682398426E-5,b0,b1,b2);
   BesselMNextCheb(z,3.36911647825569408990E-3,b0,b1,b2);
   BesselMNextCheb(z,8.04490411014108831608E-1,b0,b1,b2);
//--- calculation
   v=0.5*(b0-b2);
//--- return result
   return(MathExp(x)*v/MathSqrt(x));
  }
//+------------------------------------------------------------------+
//| Modified Bessel function of order one                            |
//| Returns modified Bessel function of order one of the             |
//| argument.                                                        |
//| The function is defined as i1(x) = -i j1( ix ).                  |
//| The range is partitioned into the two intervals [0,8] and        |
//| (8, infinity).  Chebyshev polynomial expansions are employed     |
//| in each interval.                                                |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0, 30       30000       1.9e-15     2.1e-16         |
//+------------------------------------------------------------------+
static double CBessel::BesselI1(double x)
  {
//--- create variables
   double y=0;
   double z=0;
   double v=0;
   double b0=0;
   double b1=0;
   double b2=0;
//--- initialization
   z=MathAbs(x);
//--- check
   if(z<=8.0)
     {
      //--- initialization
      y=z/2.0-2.0;
      //--- function call
      BesselM1FirstCheb(2.77791411276104639959E-18,b0,b1,b2);
      BesselM1NextCheb(y,-2.11142121435816608115E-17,b0,b1,b2);
      BesselM1NextCheb(y,1.55363195773620046921E-16,b0,b1,b2);
      BesselM1NextCheb(y,-1.10559694773538630805E-15,b0,b1,b2);
      BesselM1NextCheb(y,7.60068429473540693410E-15,b0,b1,b2);
      BesselM1NextCheb(y,-5.04218550472791168711E-14,b0,b1,b2);
      BesselM1NextCheb(y,3.22379336594557470981E-13,b0,b1,b2);
      BesselM1NextCheb(y,-1.98397439776494371520E-12,b0,b1,b2);
      BesselM1NextCheb(y,1.17361862988909016308E-11,b0,b1,b2);
      BesselM1NextCheb(y,-6.66348972350202774223E-11,b0,b1,b2);
      BesselM1NextCheb(y,3.62559028155211703701E-10,b0,b1,b2);
      BesselM1NextCheb(y,-1.88724975172282928790E-9,b0,b1,b2);
      BesselM1NextCheb(y,9.38153738649577178388E-9,b0,b1,b2);
      BesselM1NextCheb(y,-4.44505912879632808065E-8,b0,b1,b2);
      BesselM1NextCheb(y,2.00329475355213526229E-7,b0,b1,b2);
      BesselM1NextCheb(y,-8.56872026469545474066E-7,b0,b1,b2);
      BesselM1NextCheb(y,3.47025130813767847674E-6,b0,b1,b2);
      BesselM1NextCheb(y,-1.32731636560394358279E-5,b0,b1,b2);
      BesselM1NextCheb(y,4.78156510755005422638E-5,b0,b1,b2);
      BesselM1NextCheb(y,-1.61760815825896745588E-4,b0,b1,b2);
      BesselM1NextCheb(y,5.12285956168575772895E-4,b0,b1,b2);
      BesselM1NextCheb(y,-1.51357245063125314899E-3,b0,b1,b2);
      BesselM1NextCheb(y,4.15642294431288815669E-3,b0,b1,b2);
      BesselM1NextCheb(y,-1.05640848946261981558E-2,b0,b1,b2);
      BesselM1NextCheb(y,2.47264490306265168283E-2,b0,b1,b2);
      BesselM1NextCheb(y,-5.29459812080949914269E-2,b0,b1,b2);
      BesselM1NextCheb(y,1.02643658689847095384E-1,b0,b1,b2);
      BesselM1NextCheb(y,-1.76416518357834055153E-1,b0,b1,b2);
      BesselM1NextCheb(y,2.52587186443633654823E-1,b0,b1,b2);
      //--- calculation
      v=0.5*(b0-b2);
      z=v*z*MathExp(z);
     }
   else
     {
      //--- initialization
      y=32.0/z-2.0;
      //--- function calls
      BesselM1FirstCheb(7.51729631084210481353E-18,b0,b1,b2);
      BesselM1NextCheb(y,4.41434832307170791151E-18,b0,b1,b2);
      BesselM1NextCheb(y,-4.65030536848935832153E-17,b0,b1,b2);
      BesselM1NextCheb(y,-3.20952592199342395980E-17,b0,b1,b2);
      BesselM1NextCheb(y,2.96262899764595013876E-16,b0,b1,b2);
      BesselM1NextCheb(y,3.30820231092092828324E-16,b0,b1,b2);
      BesselM1NextCheb(y,-1.88035477551078244854E-15,b0,b1,b2);
      BesselM1NextCheb(y,-3.81440307243700780478E-15,b0,b1,b2);
      BesselM1NextCheb(y,1.04202769841288027642E-14,b0,b1,b2);
      BesselM1NextCheb(y,4.27244001671195135429E-14,b0,b1,b2);
      BesselM1NextCheb(y,-2.10154184277266431302E-14,b0,b1,b2);
      BesselM1NextCheb(y,-4.08355111109219731823E-13,b0,b1,b2);
      BesselM1NextCheb(y,-7.19855177624590851209E-13,b0,b1,b2);
      BesselM1NextCheb(y,2.03562854414708950722E-12,b0,b1,b2);
      BesselM1NextCheb(y,1.41258074366137813316E-11,b0,b1,b2);
      BesselM1NextCheb(y,3.25260358301548823856E-11,b0,b1,b2);
      BesselM1NextCheb(y,-1.89749581235054123450E-11,b0,b1,b2);
      BesselM1NextCheb(y,-5.58974346219658380687E-10,b0,b1,b2);
      BesselM1NextCheb(y,-3.83538038596423702205E-9,b0,b1,b2);
      BesselM1NextCheb(y,-2.63146884688951950684E-8,b0,b1,b2);
      BesselM1NextCheb(y,-2.51223623787020892529E-7,b0,b1,b2);
      BesselM1NextCheb(y,-3.88256480887769039346E-6,b0,b1,b2);
      BesselM1NextCheb(y,-1.10588938762623716291E-4,b0,b1,b2);
      BesselM1NextCheb(y,-9.76109749136146840777E-3,b0,b1,b2);
      BesselM1NextCheb(y,7.78576235018280120474E-1,b0,b1,b2);
      //--- calculation
      v=0.5*(b0-b2);
      z=v*MathExp(z)/MathSqrt(z);
     }
//--- check
   if(x<0.0)
      z=-z;
//--- return result
   return(z);
  }
//+------------------------------------------------------------------+
//| Modified Bessel function, second kind, order zero                |
//| Returns modified Bessel function of the second kind              |
//| of order zero of the argument.                                   |
//| The range is partitioned into the two intervals [0,8] and        |
//| (8, infinity).  Chebyshev polynomial expansions are employed     |
//| in each interval.                                                |
//| ACCURACY:                                                        |
//| Tested at 2000 random points between 0 and 8.  Peak absolute     |
//| error (relative when K0 > 1) was 1.46e-14; rms, 4.26e-15.        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0, 30       30000       1.2e-15     1.6e-16         |
//+------------------------------------------------------------------+
static double CBessel::BesselK0(double x)
  {
//--- create variables
   double result=0;
   double y=0;
   double z=0;
   double v=0;
   double b0=0;
   double b1=0;
   double b2=0;
//--- check
   if(!CAp::Assert(x>0.0,__FUNCTION__+": x<=0"))
      return(EMPTY_VALUE);
//--- check
   if(x<=2.0)
     {
      //--- initialization
      y=x*x-2.0;
      //--- function calls
      BesselMFirstCheb(1.37446543561352307156E-16,b0,b1,b2);
      BesselMNextCheb(y,4.25981614279661018399E-14,b0,b1,b2);
      BesselMNextCheb(y,1.03496952576338420167E-11,b0,b1,b2);
      BesselMNextCheb(y,1.90451637722020886025E-9,b0,b1,b2);
      BesselMNextCheb(y,2.53479107902614945675E-7,b0,b1,b2);
      BesselMNextCheb(y,2.28621210311945178607E-5,b0,b1,b2);
      BesselMNextCheb(y,1.26461541144692592338E-3,b0,b1,b2);
      BesselMNextCheb(y,3.59799365153615016266E-2,b0,b1,b2);
      BesselMNextCheb(y,3.44289899924628486886E-1,b0,b1,b2);
      BesselMNextCheb(y,-5.35327393233902768720E-1,b0,b1,b2);
      v=0.5*(b0-b2);
      v=v-MathLog(0.5*x)*BesselI0(x);
     }
   else
     {
      //--- initialization
      z=8.0/x-2.0;
      //--- function call
      BesselMFirstCheb(5.30043377268626276149E-18,b0,b1,b2);
      BesselMNextCheb(z,-1.64758043015242134646E-17,b0,b1,b2);
      BesselMNextCheb(z,5.21039150503902756861E-17,b0,b1,b2);
      BesselMNextCheb(z,-1.67823109680541210385E-16,b0,b1,b2);
      BesselMNextCheb(z,5.51205597852431940784E-16,b0,b1,b2);
      BesselMNextCheb(z,-1.84859337734377901440E-15,b0,b1,b2);
      BesselMNextCheb(z,6.34007647740507060557E-15,b0,b1,b2);
      BesselMNextCheb(z,-2.22751332699166985548E-14,b0,b1,b2);
      BesselMNextCheb(z,8.03289077536357521100E-14,b0,b1,b2);
      BesselMNextCheb(z,-2.98009692317273043925E-13,b0,b1,b2);
      BesselMNextCheb(z,1.14034058820847496303E-12,b0,b1,b2);
      BesselMNextCheb(z,-4.51459788337394416547E-12,b0,b1,b2);
      BesselMNextCheb(z,1.85594911495471785253E-11,b0,b1,b2);
      BesselMNextCheb(z,-7.95748924447710747776E-11,b0,b1,b2);
      BesselMNextCheb(z,3.57739728140030116597E-10,b0,b1,b2);
      BesselMNextCheb(z,-1.69753450938905987466E-9,b0,b1,b2);
      BesselMNextCheb(z,8.57403401741422608519E-9,b0,b1,b2);
      BesselMNextCheb(z,-4.66048989768794782956E-8,b0,b1,b2);
      BesselMNextCheb(z,2.76681363944501510342E-7,b0,b1,b2);
      BesselMNextCheb(z,-1.83175552271911948767E-6,b0,b1,b2);
      BesselMNextCheb(z,1.39498137188764993662E-5,b0,b1,b2);
      BesselMNextCheb(z,-1.28495495816278026384E-4,b0,b1,b2);
      BesselMNextCheb(z,1.56988388573005337491E-3,b0,b1,b2);
      BesselMNextCheb(z,-3.14481013119645005427E-2,b0,b1,b2);
      BesselMNextCheb(z,2.44030308206595545468E0,b0,b1,b2);
      //--- calculation
      v=0.5*(b0-b2);
      v=v*MathExp(-x)/MathSqrt(x);
     }
//--- return result
   return(v);
  }
//+------------------------------------------------------------------+
//| Modified Bessel function, second kind, order one                 |
//| Computes the modified Bessel function of the second kind         |
//| of order one of the argument.                                    |
//| The range is partitioned into the two intervals [0,2] and        |
//| (2, infinity).  Chebyshev polynomial expansions are employed     |
//| in each interval.                                                |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0, 30       30000       1.2e-15     1.6e-16         |
//+------------------------------------------------------------------+
static double CBessel::BesselK1(double x)
  {
//--- create variables
   double result=0;
   double y=0;
   double z=0;
   double v=0;
   double b0=0;
   double b1=0;
   double b2=0;
//--- initialization
   z=0.5*x;
//--- check
   if(!CAp::Assert(z>0.0,__FUNCTION__+": domain error"))
      return(EMPTY_VALUE);
//--- check
   if(x<=2.0)
     {
      //--- initialization
      y=x*x-2.0;
      //--- function calls
      BesselM1FirstCheb(-7.02386347938628759343E-18,b0,b1,b2);
      BesselM1NextCheb(y,-2.42744985051936593393E-15,b0,b1,b2);
      BesselM1NextCheb(y,-6.66690169419932900609E-13,b0,b1,b2);
      BesselM1NextCheb(y,-1.41148839263352776110E-10,b0,b1,b2);
      BesselM1NextCheb(y,-2.21338763073472585583E-8,b0,b1,b2);
      BesselM1NextCheb(y,-2.43340614156596823496E-6,b0,b1,b2);
      BesselM1NextCheb(y,-1.73028895751305206302E-4,b0,b1,b2);
      BesselM1NextCheb(y,-6.97572385963986435018E-3,b0,b1,b2);
      BesselM1NextCheb(y,-1.22611180822657148235E-1,b0,b1,b2);
      BesselM1NextCheb(y,-3.53155960776544875667E-1,b0,b1,b2);
      BesselM1NextCheb(y,1.52530022733894777053E0,b0,b1,b2);
      v=0.5*(b0-b2);
      //--- get result
      result=MathLog(z)*BesselI1(x)+v/x;
     }
   else
     {
      //--- initialization
      y=8.0/x-2.0;
      //--- function calls
      BesselM1FirstCheb(-5.75674448366501715755E-18,b0,b1,b2);
      BesselM1NextCheb(y,1.79405087314755922667E-17,b0,b1,b2);
      BesselM1NextCheb(y,-5.68946255844285935196E-17,b0,b1,b2);
      BesselM1NextCheb(y,1.83809354436663880070E-16,b0,b1,b2);
      BesselM1NextCheb(y,-6.05704724837331885336E-16,b0,b1,b2);
      BesselM1NextCheb(y,2.03870316562433424052E-15,b0,b1,b2);
      BesselM1NextCheb(y,-7.01983709041831346144E-15,b0,b1,b2);
      BesselM1NextCheb(y,2.47715442448130437068E-14,b0,b1,b2);
      BesselM1NextCheb(y,-8.97670518232499435011E-14,b0,b1,b2);
      BesselM1NextCheb(y,3.34841966607842919884E-13,b0,b1,b2);
      BesselM1NextCheb(y,-1.28917396095102890680E-12,b0,b1,b2);
      BesselM1NextCheb(y,5.13963967348173025100E-12,b0,b1,b2);
      BesselM1NextCheb(y,-2.12996783842756842877E-11,b0,b1,b2);
      BesselM1NextCheb(y,9.21831518760500529508E-11,b0,b1,b2);
      BesselM1NextCheb(y,-4.19035475934189648750E-10,b0,b1,b2);
      BesselM1NextCheb(y,2.01504975519703286596E-9,b0,b1,b2);
      BesselM1NextCheb(y,-1.03457624656780970260E-8,b0,b1,b2);
      BesselM1NextCheb(y,5.74108412545004946722E-8,b0,b1,b2);
      BesselM1NextCheb(y,-3.50196060308781257119E-7,b0,b1,b2);
      BesselM1NextCheb(y,2.40648494783721712015E-6,b0,b1,b2);
      BesselM1NextCheb(y,-1.93619797416608296024E-5,b0,b1,b2);
      BesselM1NextCheb(y,1.95215518471351631108E-4,b0,b1,b2);
      BesselM1NextCheb(y,-2.85781685962277938680E-3,b0,b1,b2);
      BesselM1NextCheb(y,1.03923736576817238437E-1,b0,b1,b2);
      BesselM1NextCheb(y,2.72062619048444266945E0,b0,b1,b2);
      v=0.5*(b0-b2);
      //--- get result
      result=MathExp(-x)*v/MathSqrt(x);
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Modified Bessel function, second kind, integer order             |
//| Returns modified Bessel function of the second kind              |
//| of order n of the argument.                                      |
//| The range is partitioned into the two intervals [0,9.55] and     |
//| (9.55, infinity).  An ascending power series is used in the      |
//| low range, and an asymptotic expansion in the high range.        |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0,30        90000       1.8e-8      3.0e-10         |
//| Error is high only near the crossover point x = 9.55             |
//| between the two expansions used.                                 |
//+------------------------------------------------------------------+
static double CBessel::BesselKN(int nn,double x)
  {
//--- create variables
   double k=0;
   double kf=0;
   double nk1f=0;
   double nkf=0;
   double zn=0;
   double t=0;
   double s=0;
   double z0=0;
   double z=0;
   double ans=0;
   double fn=0;
   double pn=0;
   double pk=0;
   double zmn=0;
   double tlg=0;
   double tox=0;
   int    i=0;
   int    n=0;
   double eul=0;
//--- initialization
   eul=5.772156649015328606065e-1;
//--- check
   if(nn<0)
      n=-nn;
   else
      n=nn;
//--- check
   if(!CAp::Assert(n<=31,__FUNCTION__+": overflow"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(x>0.0,__FUNCTION__+":domain error"))
      return(EMPTY_VALUE);
//--- check
   if(x<=9.55)
     {
      //--- initialization
      ans=0.0;
      z0=0.25*x*x;
      fn=1.0;
      pn=0.0;
      zmn=1.0;
      tox=2.0/x;
      //--- check
      if(n>0)
        {
         //--- calculation
         pn=-eul;
         k=1.0;
         for(i=1;i<=n-1;i++)
           {
            pn=pn+1.0/k;
            k=k+1.0;
            fn=fn*k;
           }
         zmn=tox;
         //--- check
         if(n==1)
            ans=1.0/x;
         else
           {
            //--- change values
            nk1f=fn/n;
            kf=1.0;
            s=nk1f;
            z=-z0;
            zn=1.0;
            //--- calculation
            for(i=1;i<=n-1;i++)
              {
               nk1f=nk1f/(n-i);
               kf=kf*i;
               zn=zn*z;
               t=nk1f*zn/kf;
               s=s+t;
               //--- check
               if(!CAp::Assert(CMath::m_maxrealnumber-MathAbs(t)>MathAbs(s),__FUNCTION__+": overflow"))
                  return(EMPTY_VALUE);
               //--- check
               if(!CAp::Assert(!(tox>1.0 && CMath::m_maxrealnumber/tox<zmn),__FUNCTION__+": overflow"))
                  return(EMPTY_VALUE);
               zmn=zmn*tox;
              }
            //--- change values
            s=s*0.5;
            t=MathAbs(s);
            //--- check
            if(!CAp::Assert(!(zmn>1.0 && CMath::m_maxrealnumber/zmn<t),__FUNCTION__+": overflow"))
               return(EMPTY_VALUE);
            //--- check
            if(!CAp::Assert(!(t>1.0 && CMath::m_maxrealnumber/t<zmn),__FUNCTION__+": overflow"))
               return(EMPTY_VALUE);
            ans=s*zmn;
           }
        }
      //--- change values
      tlg=2.0*MathLog(0.5*x);
      pk=-eul;
      //--- check
      if(n==0)
        {
         pn=pk;
         t=1.0;
        }
      else
        {
         pn=pn+1.0/n;
         t=1.0/fn;
        }
      s=(pk+pn-tlg)*t;
      k=1.0;
      //--- cycle
      do
        {
         t=t*(z0/(k*(k+n)));
         pk=pk+1.0/k;
         pn=pn+1.0/(k+n);
         s=s+(pk+pn-tlg)*t;
         k=k+1.0;
        }
      while(MathAbs(t/s)>CMath::m_machineepsilon);
      s=0.5*s/zmn;
      //--- check
      if(n%2!=0)
        {
         s=-s;
        }
      ans=ans+s;
      //--- return result
      return(ans);
     }
//--- check
   if(x>MathLog(CMath::m_maxrealnumber))
      return(0);
//--- change values
   k=n;
   pn=4.0*k*k;
   pk=1.0;
   z0=8.0*x;
   fn=1.0;
   t=1.0;
   s=t;
   nkf=CMath::m_maxrealnumber;
   i=0;
//--- cycle
   do
     {
      z=pn-pk*pk;
      t=t*z/(fn*z0);
      nk1f=MathAbs(t);
      //--- check
      if(i>=n && nk1f>(double)(nkf))
         break;
      //--- change values
      nkf=nk1f;
      s=s+t;
      fn=fn+1.0;
      pk=pk+2.0;
      i=i+1;
     }
   while(MathAbs(t/s)>CMath::m_machineepsilon);
//--- return result
   return(MathExp(-x)*MathSqrt(M_PI/(2.0*x))*s);
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CBessel::BesselMFirstCheb(const double c,double &b0,double &b1,
                                      double &b2)
  {
//--- change values
   b0=c;
   b1=0.0;
   b2=0.0;
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CBessel::BesselMNextCheb(const double x,const double c,double &b0,
                                     double &b1,double &b2)
  {
//--- change values
   b2=b1;
   b1=b0;
   b0=x*b1-b2+c;
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CBessel::BesselM1FirstCheb(const double c,double &b0,double &b1,
                                       double &b2)
  {
//--- change values
   b0=c;
   b1=0.0;
   b2=0.0;
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CBessel::BesselM1NextCheb(const double x,const double c,double &b0,
                                      double &b1,double &b2)
  {
//--- change values
   b2=b1;
   b1=b0;
   b0=x*b1-b2+c;
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CBessel::BesselAsympt0(const double x,double &pzero,double &qzero)
  {
//--- create variables
   double xsq=0;
   double p2=0;
   double q2=0;
   double p3=0;
   double q3=0;
//--- initialization
   pzero=0;
   qzero=0;
//--- calculation
   xsq=64.0/(x*x);
   p2=0.0;
   p2=2485.271928957404011288128951+xsq*p2;
   p2=153982.6532623911470917825993+xsq*p2;
   p2=2016135.283049983642487182349+xsq*p2;
   p2=8413041.456550439208464315611+xsq*p2;
   p2=12332384.76817638145232406055+xsq*p2;
   p2=5393485.083869438325262122897+xsq*p2;
   q2=1.0;
   q2=2615.700736920839685159081813+xsq*q2;
   q2=156001.7276940030940592769933+xsq*q2;
   q2=2025066.801570134013891035236+xsq*q2;
   q2=8426449.050629797331554404810+xsq*q2;
   q2=12338310.22786324960844856182+xsq*q2;
   q2=5393485.083869438325560444960+xsq*q2;
   p3=-0.0;
   p3=-4.887199395841261531199129300+xsq*p3;
   p3=-226.2630641933704113967255053+xsq*p3;
   p3=-2365.956170779108192723612816+xsq*p3;
   p3=-8239.066313485606568803548860+xsq*p3;
   p3=-10381.41698748464093880530341+xsq*p3;
   p3=-3984.617357595222463506790588+xsq*p3;
   q3=1.0;
   q3=408.7714673983499223402830260+xsq*q3;
   q3=15704.89191515395519392882766+xsq*q3;
   q3=156021.3206679291652539287109+xsq*q3;
   q3=533291.3634216897168722255057+xsq*q3;
   q3=666745.4239319826986004038103+xsq*q3;
   q3=255015.5108860942382983170882+xsq*q3;
//--- get result
   pzero=p2/q2;
   qzero=8*p3/q3/x;
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CBessel::BesselAsympt1(const double x,double &pzero,double &qzero)
  {
//--- create variables
   double xsq=0;
   double p2=0;
   double q2=0;
   double p3=0;
   double q3=0;
//--- initialization
   pzero=0;
   qzero=0;
//--- calculation
   xsq=64.0/(x*x);
   p2=-1611.616644324610116477412898;
   p2=-109824.0554345934672737413139+xsq*p2;
   p2=-1523529.351181137383255105722+xsq*p2;
   p2=-6603373.248364939109255245434+xsq*p2;
   p2=-9942246.505077641195658377899+xsq*p2;
   p2=-4435757.816794127857114720794+xsq*p2;
   q2=1.0;
   q2=-1455.009440190496182453565068+xsq*q2;
   q2=-107263.8599110382011903063867+xsq*q2;
   q2=-1511809.506634160881644546358+xsq*q2;
   q2=-6585339.479723087072826915069+xsq*q2;
   q2=-9934124.389934585658967556309+xsq*q2;
   q2=-4435757.816794127856828016962+xsq*q2;
   p3=35.26513384663603218592175580;
   p3=1706.375429020768002061283546+xsq*p3;
   p3=18494.26287322386679652009819+xsq*p3;
   p3=66178.83658127083517939992166+xsq*p3;
   p3=85145.16067533570196555001171+xsq*p3;
   p3=33220.91340985722351859704442+xsq*p3;
   q3=1.0;
   q3=863.8367769604990967475517183+xsq*q3;
   q3=37890.22974577220264142952256+xsq*q3;
   q3=400294.4358226697511708610813+xsq*q3;
   q3=1419460.669603720892855755253+xsq*q3;
   q3=1819458.042243997298924553839+xsq*q3;
   q3=708712.8194102874357377502472+xsq*q3;
//--- get result
   pzero=p2/q2;
   qzero=8*p3/q3/x;
  }
//+------------------------------------------------------------------+
//| Beta function                                                    |
//+------------------------------------------------------------------+
class CBetaF
  {
public:
   //--- constructor, destructor
                     CBetaF(void);
                    ~CBetaF(void);
   //--- method
   static double     Beta(const double a,const double b);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CBetaF::CBetaF(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CBetaF::~CBetaF(void)
  {

  }
//+------------------------------------------------------------------+
//| Beta function                                                    |
//|                   -     -                                        |
//|                  | (a) | (b)                                     |
//| beta( a, b )  =  -----------.                                    |
//|                     -                                            |
//|                    | (a+b)                                       |
//| For large arguments the logarithm of the function is             |
//| evaluated using lgam(), then exponentiated.                      |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE       0,30       30000       8.1e-14     1.1e-14         |
//+------------------------------------------------------------------+
static double CBetaF::Beta(const double a,const double b)
  {
//--- create variables
   double y=0;
   double sg=0;
   double s=0;
//--- initialization
   sg=1;
//--- check
   if(!CAp::Assert(a>0.0 || a!=(double)((int)MathFloor(a)),__FUNCTION__+": overflow"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(b>0.0 || b!=(double)((int)MathFloor(b)),__FUNCTION__+": overflow"))
      return(EMPTY_VALUE);
   y=a+b;
//--- check
   if(MathAbs(y)>171.624376956302725)
     {
      //--- calculation
      y=CGammaFunc::LnGamma(y,s);
      sg=sg*s;
      y=CGammaFunc::LnGamma(b,s)-y;
      sg=sg*s;
      y=CGammaFunc::LnGamma(a,s)+y;
      sg=sg*s;
      //--- check
      if(!CAp::Assert(y<=(double)(MathLog(CMath::m_maxrealnumber)),__FUNCTION__+": overflow"))
         return(EMPTY_VALUE);
      //--- return result
      return(sg*MathExp(y));
     }
//--- function call
   y=CGammaFunc::GammaFunc(y);
//--- check
   if(!CAp::Assert(y!=0.0,"Overflow in Beta"))
      return(EMPTY_VALUE);
//--- check
   if(a>b)
     {
      //--- calculation
      y=CGammaFunc::GammaFunc(a)/y;
      y=y*CGammaFunc::GammaFunc(b);
     }
   else
     {
      //--- calculation
      y=CGammaFunc::GammaFunc(b)/y;
      y=y*CGammaFunc::GammaFunc(a);
     }
//--- return result
   return(y);
  }
//+------------------------------------------------------------------+
//| Incomplete beta function                                         |
//+------------------------------------------------------------------+
class CIncBetaF
  {
private:
   //--- private methods
   static double     IncompleteBetaFracExpans(const double a,const double b,const double x,const double big,const double biginv);
   static double     IncompleteBetaFracExpans2(const double a,const double b,const double x,const double big,const double biginv);
   static double     IncompleteBetaPowSeries(const double a,const double b,const double x,const double maxgam);
public:
   //--- constructor, destructor
                     CIncBetaF(void);
                    ~CIncBetaF(void);
   //--- public methods
   static double     IncompleteBeta(double a,double b,double x);
   static double     InvIncompleteBeta(const double a,double b,double y);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CIncBetaF::CIncBetaF(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CIncBetaF::~CIncBetaF(void)
  {

  }
//+------------------------------------------------------------------+
//| Incomplete beta integral                                         |
//| Returns incomplete beta integral of the arguments, evaluated     |
//| from zero to x.  The function is defined as                      |
//|                  x                                               |
//|     -            -                                               |
//|    | (a+b)      | |  a-1     b-1                                 |
//|  -----------    |   t   (1-t)   dt.                              |
//|   -     -     | |                                                |
//|  | (a) | (b)   -                                                 |
//|                 0                                                |
//| The domain of definition is 0 <= x <= 1.  In this                |
//| implementation a and b are restricted to positive values.        |
//| The integral from x to 1 may be obtained by the symmetry         |
//| relation                                                         |
//|    1 - incbet( a, b, x )  =  incbet( b, a, 1-x ).                |
//| The integral is evaluated by a continued fraction expansion      |
//| or, when b*x is small, by a power series.                        |
//| ACCURACY:                                                        |
//| Tested at uniformly distributed random points (a,b,x) with a and |
//| b in "domain" and x between 0 and 1.                             |
//|                                        Relative error            |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0,5         10000       6.9e-15     4.5e-16         |
//|    IEEE      0,85       250000       2.2e-13     1.7e-14         |
//|    IEEE      0,1000      30000       5.3e-12     6.3e-13         |
//|    IEEE      0,10000    250000       9.3e-11     7.1e-12         |
//|    IEEE      0,100000    10000       8.7e-10     4.8e-11         |
//| Outputs smaller than the IEEE gradual underflow threshold        |
//| were excluded from these statistics.                             |
//+------------------------------------------------------------------+
static double CIncBetaF::IncompleteBeta(double a,double b,double x)
  {
//--- check
   if(!CAp::Assert(a>0 && b>0,__FUNCTION__+": the error variable"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(x>=0 && x<=1,__FUNCTION__+": the error variable"))
      return(EMPTY_VALUE);
//--- check
   if(x==0)
      return(0);
//--- check
   if(x==1)
      return(1);
//--- create variables
   double t;
   double xc;
   double w;
   double y;
   int    flag=0;
   double sg=0;
   double big=4.503599627370496e15;
   double biginv=2.22044604925031308085e-16;
   double maxgam=171.624376956302725;
   double minlog=MathLog(CMath::m_minrealnumber);
   double maxlog=MathLog(CMath::m_maxrealnumber);

//--- check
   if(b*x<=1.0 && x<=0.95)
      return(IncompleteBetaPowSeries(a,b,x,maxgam));
//--- change value
   w=1.0-x;
//--- check
   if(x>a/(a+b))
     {
      flag=1;
      t=a;
      a=b;
      b=t;
      xc=x;
      x=w;
     }
   else
      xc=w;
//--- check
   if(flag==1 && b*x<=1.0 && x<=0.95)
     {
      t=IncompleteBetaPowSeries(a,b,x,maxgam);
      //--- check
      if(t<=CMath::m_machineepsilon)
         return(1.0-CMath::m_machineepsilon);
      else
         return(1.0-t);
     }
//--- change value
   y=x*(a+b-2.0)-(a-1.0);
//--- check
   if(y<0.0)
      w=IncompleteBetaFracExpans(a,b,x,big,biginv);
   else
      w=IncompleteBetaFracExpans2(a,b,x,big,biginv)/xc;
//--- change values
   y=a*MathLog(x);
   t=b*MathLog(xc);
//--- check
   if(a+b<maxgam && MathAbs(y)<maxlog && MathAbs(t)<maxlog)
     {
      t=MathPow(xc,b);
      t*=MathPow(x,a);
      t/=a;
      t*=w;
      t*=CGammaFunc::GammaFunc(a+b)/(CGammaFunc::GammaFunc(a)*CGammaFunc::GammaFunc(b));
      //--- check
      if(flag==1)
        {
         //--- check
         if(t<=CMath::m_machineepsilon)
            return(1.0-CMath::m_machineepsilon);
         else
            return(1.0-t);
        }
      else
         return(t);
     }
//--- change values
   y+=t+CGammaFunc::LnGamma(a+b,sg)-CGammaFunc::LnGamma(a,sg)-CGammaFunc::LnGamma(b,sg);
   y+=MathLog(w/a);
//--- check
   if(y<minlog)
      t=0.0;
   else
      t=MathExp(y);
//--- check
   if(flag==1)
     {
      if(t<=CMath::m_machineepsilon)
         t=1.0-CMath::m_machineepsilon;
      else
         t=1.0-t;
     }
//--- return result
   return(t);
  }
//+------------------------------------------------------------------+
//| Inverse of imcomplete beta integral                              |
//| Given y, the function finds x such that                          |
//|  incbet( a, b, x ) = y .                                         |
//| The routine performs interval halving or Newton iterations to    |
//| find the root of incbet(a,b,x) - y = 0.                          |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//|                x     a,b                                         |
//| arithmetic   domain  domain  # trials    peak       rms          |
//|    IEEE      0,1    .5,10000   50000    5.8e-12   1.3e-13        |
//|    IEEE      0,1   .25,100    100000    1.8e-13   3.9e-15        |
//|    IEEE      0,1     0,5       50000    1.1e-12   5.5e-15        |
//| With a and b constrained to half-integer or integer values:      |
//|    IEEE      0,1    .5,10000   50000    5.8e-12   1.1e-13        |
//|    IEEE      0,1    .5,100    100000    1.7e-14   7.9e-16        |
//| With a = .5, b constrained to half-integer or integer values:    |
//|    IEEE      0,1    .5,10000   10000    8.3e-11   1.0e-11        |
//+------------------------------------------------------------------+
static double CIncBetaF::InvIncompleteBeta(const double a,double b,double y)
  {
//--- check
   if(!CAp::Assert(y>=0 && y<=1,__FUNCTION__+": the error variable"))
      return(EMPTY_VALUE);
//--- check
   if(y==0)
      return(0);
//--- check
   if(y==1.0)
      return(1);
//--- create variables
   double aaa=0;
   double bbb=0;
   double y0=0;
   double d=0;
   double yyy=0;
   double x=0;
   double x0=0.0;
   double x1=1.0;
   double lgm=0;
   double yp=0;
   double di=0;
   double dithresh=0;
   double yl=0.0;
   double yh=1.0;
   double xt=0;
   int    i=0;
   int    rflg=0;
   int    dir=0;
   int    nflg=0;
   double s=0;
   int    mainlooppos=0;
   int    ihalve=1;
   int    ihalvecycle=2;
   int    newt=3;
   int    newtcycle=4;
   int    breaknewtcycle=5;
   int    breakihalvecycle=6;
//--- main
   while(true)
     {
      //--- start
      if(mainlooppos==0)
        {
         if(a<=1.0 || b<=1.0)
           {
            dithresh=1.0e-6;
            rflg=0;
            aaa=a;
            bbb=b;
            y0=y;
            x=aaa/(aaa+bbb);
            yyy=IncompleteBeta(aaa,bbb,x);
            mainlooppos=ihalve;
            continue;
           }
         else
            dithresh=1.0e-4;
         //--- change value
         yp=-CNormalDistr::InvNormalDistribution(y);
         //--- check
         if(y>0.5)
           {
            rflg=1;
            aaa=b;
            bbb=a;
            y0=1.0-y;
            yp=-yp;
           }
         else
           {
            rflg=0;
            aaa=a;
            bbb=b;
            y0=y;
           }
         //--- change values
         lgm=(yp*yp-3.0)/6.0;
         x=2.0/(1.0/(2.0*aaa-1.0)+1.0/(2.0*bbb-1.0));
         d=yp*MathSqrt(x+lgm)/x-(1.0/(2.0*bbb-1.0)-1.0/(2.0*aaa-1.0))*(lgm+5.0/6.0-2.0/(3.0*x));
         d=2.0*d;
         //--- check
         if(d<MathLog(CMath::m_minrealnumber))
           {
            x=0;
            break;
           }
         //--- change values
         x=aaa/(aaa+bbb*MathExp(d));
         yyy=IncompleteBeta(aaa,bbb,x);
         yp=(yyy-y0)/y0;
         //--- check
         if(MathAbs(yp)<0.2)
           {
            mainlooppos=newt;
            continue;
           }
         mainlooppos=ihalve;
         continue;
        }
      //--- check
      if(mainlooppos==ihalve)
        {
         dir=0;
         di=0.5;
         i=0;
         mainlooppos=ihalvecycle;
         continue;
        }
      //--- check
      if(mainlooppos==ihalvecycle)
        {
         //--- check
         if(i<=99)
           {
            //--- check
            if(i!=0)
              {
               x=x0+di*(x1-x0);
               //--- check
               if(x==1.0)
                  x=1.0-CMath::m_machineepsilon;
               //--- check
               if(x==0.0)
                 {
                  di=0.5;
                  x=x0+di*(x1-x0);
                  //--- check
                  if(x==0.0)
                     break;
                 }
               //--- change values
               yyy=IncompleteBeta(aaa,bbb,x);
               yp=(x1-x0)/(x1+x0);
               //--- check
               if(MathAbs(yp)<dithresh)
                 {
                  mainlooppos=newt;
                  continue;
                 }
               yp=(yyy-y0)/y0;
               //--- check
               if(MathAbs(yp)<dithresh)
                 {
                  mainlooppos=newt;
                  continue;
                 }
              }
            //--- check
            if(yyy<y0)
              {
               x0=x;
               yl=yyy;
               //--- check
               if(dir<0)
                 {
                  dir=0;
                  di=0.5;
                 }
               else
                 {
                  //--- check
                  if(dir>3)
                     di=1.0-(1.0-di)*(1.0-di);
                  else
                    {
                     //--- check
                     if(dir>1)
                        di=0.5*di+0.5;
                     else
                        di=(y0-yyy)/(yh-yl);
                    }
                 }
               dir++;
               //--- check
               if(x0>0.75)
                 {
                  //--- check
                  if(rflg==1)
                    {
                     rflg=0;
                     aaa=a;
                     bbb=b;
                     y0=y;
                    }
                  else
                    {
                     rflg=1;
                     aaa=b;
                     bbb=a;
                     y0=1.0-y;
                    }
                  //--- change values
                  x=1.0-x;
                  yyy=IncompleteBeta(aaa,bbb,x);
                  x0=0.0;
                  yl=0.0;
                  x1=1.0;
                  yh=1.0;
                  mainlooppos=ihalve;
                  continue;
                 }
              }
            else
              {
               x1=x;
               //--- check
               if(rflg==1 && x1<CMath::m_machineepsilon)
                 {
                  x=0.0;
                  break;
                 }
               yh=yyy;
               //--- check
               if(dir>0)
                 {
                  dir=0;
                  di=0.5;
                 }
               else
                 {
                  //--- check
                  if(dir<-3)
                     di*=di;
                  else
                    {
                     //--- check
                     if(dir<-1)
                        di*=0.5;
                     else
                        di=(yyy-y0)/(yh-yl);
                    }
                 }
               dir--;
              }
            i++;
            mainlooppos=ihalvecycle;
            continue;
           }
         else
           {
            mainlooppos=breakihalvecycle;
            continue;
           }
        }
      //--- check
      if(mainlooppos==breakihalvecycle)
        {
         //--- check
         if(x0>=1.0)
           {
            x=1.0-CMath::m_machineepsilon;
            break;
           }
         //--- check
         if(x<=0.0)
           {
            x=0.0;
            break;
           }
         mainlooppos=newt;
         continue;
        }
      //--- check
      if(mainlooppos==newt)
        {
         if(nflg!=0)
            break;
         //--- change values
         nflg=1;
         lgm=CGammaFunc::LnGamma(aaa+bbb,s)-CGammaFunc::LnGamma(aaa,s)-CGammaFunc::LnGamma(bbb,s);
         i=0;
         mainlooppos=newtcycle;
         continue;
        }
      //--- check
      if(mainlooppos==newtcycle)
        {
         //--- check
         if(i<=7)
           {
            //--- check
            if(i!=0)
               yyy=IncompleteBeta(aaa,bbb,x);
            //--- check
            if(yyy<yl)
              {
               x=x0;
               yyy=yl;
              }
            else
              {
               //--- check
               if(yyy>yh)
                 {
                  x=x1;
                  yyy=yh;
                 }
               else
                 {
                  //--- check
                  if(yyy<y0)
                    {
                     x0=x;
                     yl=yyy;
                    }
                  else
                    {
                     x1=x;
                     yh=yyy;
                    }
                 }
              }
            //--- check
            if(x==1.0 || x==0.0)
              {
               mainlooppos=breaknewtcycle;
               continue;
              }
            //--- change value
            d=(aaa-1.0)*MathLog(x)+(bbb-1.0)*MathLog(1.0-x)+lgm;
            //--- check
            if(d<MathLog(CMath::m_minrealnumber))
               break;
            //--- check
            if(d>MathLog(CMath::m_maxrealnumber))
              {
               mainlooppos=breaknewtcycle;
               continue;
              }
            //--- change values
            d=MathExp(d);
            d=(yyy-y0)/d;
            xt=x-d;
            //--- check
            if(xt<=x0)
              {
               yyy=(x-x0)/(x1-x0);
               xt=x0+0.5*yyy*(x-x0);
               //--- check
               if(xt<=0.0)
                 {
                  mainlooppos=breaknewtcycle;
                  continue;
                 }
              }
            //--- check
            if(xt>=x1)
              {
               yyy=(x1-x)/(x1-x0);
               xt=x1-0.5*yyy*(x1-x);
               //--- check
               if(xt>=1.0)
                 {
                  mainlooppos=breaknewtcycle;
                  continue;
                 }
              }
            x=xt;
            //--- check
            if(MathAbs(d/x)<128.0*CMath::m_machineepsilon)
               break;
            //--- change values
            i=i+1;
            mainlooppos=newtcycle;
            continue;
           }
         else
           {
            mainlooppos=breaknewtcycle;
            continue;
           }
        }
      //--- check
      if(mainlooppos==breaknewtcycle)
        {
         dithresh=256.0*CMath::m_machineepsilon;
         mainlooppos=ihalve;
         continue;
        }
     }
//--- get result
   if(rflg!=0)
     {
      if(x<=CMath::m_machineepsilon)
         x=1.0-CMath::m_machineepsilon;
      else
         x=1.0-x;
     }
//--- return result
   return(x);
  }
//+------------------------------------------------------------------+
//| Continued fraction expansion #1 for incomplete beta integral     |
//+------------------------------------------------------------------+
static double CIncBetaF::IncompleteBetaFracExpans(const double a,const double b,
                                                  const double x,const double big,
                                                  const double biginv)
  {
//--- create variables
   double xk=0;
   double pk=0;
   double pkm1=1.0;
   double pkm2=0.0;
   double qk=0;
   double qkm1=1.0;
   double qkm2=1.0;
   double k1=a;
   double k2=a+b;
   double k3=a;
   double k4=a+1.0;
   double k5=1.0;
   double k6=b-1.0;
   double k7=k4;
   double k8=a+2.0;
   double r=1.0;
   double t=0;
   double ans=1.0;
   double thresh=3.0*CMath::m_machineepsilon;
   int    n=0;
//--- cycle
   do
     {
      xk=-(x*k1*k2/(k3*k4));
      pk=pkm1+pkm2*xk;
      qk=qkm1+qkm2*xk;
      pkm2=pkm1;
      pkm1=pk;
      qkm2=qkm1;
      qkm1=qk;
      xk=x*k5*k6/(k7*k8);
      pk=pkm1+pkm2*xk;
      qk=qkm1+qkm2*xk;
      pkm2=pkm1;
      pkm1=pk;
      qkm2=qkm1;
      qkm1=qk;
      //--- check
      if(qk!=0)
         r=pk/qk;
      //--- check
      if(r!=0)
        {
         t=MathAbs((ans-r)/r);
         ans=r;
        }
      else
         t=1.0;
      //--- check
      if(t<thresh)
         break;
      //--- change values
      k1+=1.0;
      k2+=1.0;
      k3+=2.0;
      k4+=2.0;
      k5+=1.0;
      k6-=1.0;
      k7+=2.0;
      k8+=2.0;
      //--- check
      if(MathAbs(qk)+MathAbs(pk)>big)
        {
         pkm2*=biginv;
         pkm1*=biginv;
         qkm2*=biginv;
         qkm1*=biginv;
        }
      //--- check
      if(MathAbs(qk)<biginv || MathAbs(pk)<biginv)
        {
         pkm2*=big;
         pkm1*=big;
         qkm2*=big;
         qkm1*=big;
        }
      n++;
     }
   while(n!=300);
//--- return result
   return(ans);
  }
//+------------------------------------------------------------------+
//| Continued fraction expansion #2                                  |
//| for incomplete beta integral                                     |
//+------------------------------------------------------------------+
static double CIncBetaF::IncompleteBetaFracExpans2(const double a,const double b,
                                                   const double x,const double big,
                                                   const double biginv)
  {
//--- create variables
   double xk=0;
   double pk=0;
   double pkm1=1.0;
   double pkm2=0.0;
   double qk=0;
   double qkm1=1.0;
   double qkm2=1.0;
   double k1=a;
   double k2=b-1.0;
   double k3=a;
   double k4=a+1.0;
   double k5=1.0;
   double k6=a+b;
   double k7=a+1.0;
   double k8=a+2.0;
   double r=1.0;
   double t=0;
   double ans=1.0;
   double z=x/(1.0-x);
   double thresh=3.0*CMath::m_machineepsilon;
   int    n=0;
//--- cycle
   do
     {
      xk=-(z*k1*k2/(k3*k4));
      pk=pkm1+pkm2*xk;
      qk=qkm1+qkm2*xk;
      pkm2=pkm1;
      pkm1=pk;
      qkm2=qkm1;
      qkm1=qk;
      xk=z*k5*k6/(k7*k8);
      pk=pkm1+pkm2*xk;
      qk=qkm1+qkm2*xk;
      pkm2=pkm1;
      pkm1=pk;
      qkm2=qkm1;
      qkm1=qk;
      //--- check
      if(qk!=0)
         r=pk/qk;
      //--- check
      if(r!=0)
        {
         t=MathAbs((ans-r)/r);
         ans=r;
        }
      else
         t=1.0;
      //--- check
      if(t<thresh)
         break;
      //--- change values
      k1+=1.0;
      k2-=1.0;
      k3+=2.0;
      k4+=2.0;
      k5+=1.0;
      k6+=1.0;
      k7+=2.0;
      k8+=2.0;
      //--- check
      if(MathAbs(qk)+MathAbs(pk)>big)
        {
         pkm2*=biginv;
         pkm1*=biginv;
         qkm2*=biginv;
         qkm1*=biginv;
        }
      //--- check
      if(MathAbs(qk)<biginv || MathAbs(pk)<biginv)
        {
         pkm2*=big;
         pkm1*=big;
         qkm2*=big;
         qkm1*=big;
        }
      n++;
     }
   while(n!=300);
//--- return result
   return(ans);
  }
//+------------------------------------------------------------------+
//| Power series for incomplete beta integral.                       |
//| Use when b*x is small and x not too close to 1.                  |
//+------------------------------------------------------------------+
static double CIncBetaF::IncompleteBetaPowSeries(const double a,const double b,
                                                 const double x,const double maxgam)
  {
//--- create variables
   double s=0.0;
   double u=(1.0-b)*x;
   double t=u;
   double v=u/(a+1.0);
   double n=2.0;
   double t1=v;
   double ai=1.0/a;
   double z=CMath::m_machineepsilon*ai;
   double sg=0;
//--- cycle
   while(MathAbs(v)>z)
     {
      u=(n-b)*x/n;
      t=t*u;
      v=t/(a+n);
      s=s+v;
      n=n+1.0;
     }
//--- change values
   s+=t1;
   s+=ai;
   u=a*MathLog(x);
//--- check
   if(a+b<maxgam && MathAbs(u)<MathLog(CMath::m_maxrealnumber))
     {
      t=CGammaFunc::GammaFunc(a+b)/(CGammaFunc::GammaFunc(a)*CGammaFunc::GammaFunc(b));
      s*=t*MathPow(x,a);
     }
   else
     {
      t=CGammaFunc::LnGamma(a+b,sg)-CGammaFunc::LnGamma(a,sg)-CGammaFunc::LnGamma(b,sg)+u+MathLog(s);
      //--- check
      if(t<MathLog(CMath::m_minrealnumber))
         s=0.0;
      else
         s=MathExp(t);
     }
//--- return result
   return(s);
  }
//+------------------------------------------------------------------+
//| Binomial distribution                                            |
//+------------------------------------------------------------------+
class CBinomialDistr
  {
public:
   //--- constructor, destructor
                     CBinomialDistr(void);
                    ~CBinomialDistr(void);
   //--- methods
   static double     BinomialDistribution(const int k,const int n,const double p);
   static double     BinomialComplDistribution(const int k,const int n,const double p);
   static double     InvBinomialDistribution(const int k,const int n,const double y);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CBinomialDistr::CBinomialDistr(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CBinomialDistr::~CBinomialDistr(void)
  {

  }
//+------------------------------------------------------------------+
//| Binomial distribution                                            |
//| Returns the sum of the terms 0 through k of the Binomial         |
//| probability density:                                             |
//|   k                                                              |
//|   --  ( n )   j      n-j                                         |
//|   >   (   )  p  (1-p)                                            |
//|   --  ( j )                                                      |
//|  j=0                                                             |
//| The terms are not summed directly; instead the incomplete        |
//| beta integral is employed, according to the formula              |
//| y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p ).                   |
//| The arguments must be positive, with p ranging from 0 to 1.      |
//| ACCURACY:                                                        |
//| Tested at random points (a,b,p), with p between 0 and 1.         |
//|               a,b                     Relative error:            |
//| arithmetic  domain     # trials      peak         rms            |
//|  For p between 0.001 and 1:                                      |
//|    IEEE     0,100       100000      4.3e-15     2.6e-16          |
//+------------------------------------------------------------------+
static double CBinomialDistr::BinomialDistribution(const int k,const int n,
                                                   const double p)
  {
//--- check
   if(!CAp::Assert(p>=0 && p<=1,__FUNCTION__+": the eroor variable"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(k>=-1 && k<=n,__FUNCTION__+": the eroor variable"))
      return(EMPTY_VALUE);
//--- check
   if(k==-1)
      return(0);
//--- check
   if(k==n)
      return(1);
//--- create variables
   double dk;
   double dn=n-k;
//--- check
   if(k==0)
      dk=MathPow(1.0-p,dn);
   else
     {
      dk=k+1;
      dk=CIncBetaF::IncompleteBeta(dn,dk,1.0-p);
     }
//--- return result
   return(dk);
  }
//+------------------------------------------------------------------+
//| Complemented binomial distribution                               |
//| Returns the sum of the terms k+1 through n of the Binomial       |
//| probability density:                                             |
//|   n                                                              |
//|   --  ( n )   j      n-j                                         |
//|   >   (   )  p  (1-p)                                            |
//|   --  ( j )                                                      |
//|  j=k+1                                                           |
//| The terms are not summed directly; instead the incomplete        |
//| beta integral is employed, according to the formula              |
//| y = bdtrc( k, n, p ) = incbet( k+1, n-k, p ).                    |
//| The arguments must be positive, with p ranging from 0 to 1.      |
//| ACCURACY:                                                        |
//| Tested at random points (a,b,p).                                 |
//|               a,b                     Relative error:            |
//| arithmetic  domain     # trials      peak         rms            |
//|  For p between 0.001 and 1:                                      |
//|    IEEE     0,100       100000      6.7e-15     8.2e-16          |
//|  For p between 0 and .001:                                       |
//|    IEEE     0,100       100000      1.5e-13     2.7e-15          |
//+------------------------------------------------------------------+
static double CBinomialDistr::BinomialComplDistribution(const int k,
                                                        const int n,
                                                        const double p)
  {
//--- check
   if(!CAp::Assert(p>=0 && p<=1,__FUNCTION__+": the eroor variable"))
      return(EMPTY_VALUE);
//--- check
   if(!CAp::Assert(k>=-1 && k<=n,__FUNCTION__+": the eroor variable"))
      return(EMPTY_VALUE);
//--- check
   if(k==-1)
      return(0);
//--- check
   if(k==n)
      return(1);
//--- create variables
   double dk;
   double dn=n-k;
//--- check
   if(k==0)
      dk=MathPow(1.0-p,dn);
   else
     {
      dk=k+1;
      dk=CIncBetaF::IncompleteBeta(dn,dk,1.0-p);
     }
//--- return result
   return(dk);
  }
//+------------------------------------------------------------------+
//| Inverse binomial distribution                                    |
//| Finds the event probability p such that the sum of the           |
//| terms 0 through k of the Binomial probability density            |
//| is equal to the given cumulative probability y.                  |
//| This is accomplished using the inverse beta integral             |
//| function and the relation                                        |
//| 1 - p = incbi( n-k, k+1, y ).                                    |
//| ACCURACY:                                                        |
//| Tested at random points (a,b,p).                                 |
//|               a,b                     Relative error:            |
//| arithmetic  domain     # trials      peak         rms            |
//|  For p between 0.001 and 1:                                      |
//|    IEEE     0,100       100000      2.3e-14     6.4e-16          |
//|    IEEE     0,10000     100000      6.6e-12     1.2e-13          |
//|  For p between 10^-6 and 0.001:                                  |
//|    IEEE     0,100       100000      2.0e-12     1.3e-14          |
//|    IEEE     0,10000     100000      1.5e-12     3.2e-14          |
//+------------------------------------------------------------------+
static double CBinomialDistr::InvBinomialDistribution(const int k,const int n,
                                                      const double y)
  {
//--- create variables
   double dk=0;
   double dn=0;
   double p=0;
//--- check
   if(!CAp::Assert(k>=0 && k<n,__FUNCTION__+": domain error"))
      return(EMPTY_VALUE);
//--- initialization
   dn=n-k;
//--- check
   if(k==0)
     {
      //--- check
      if(y>0.8)
         p=-CNearUnitYUnit::NUExp1m(CNearUnitYUnit::NULog1p(y-1.0)/dn);
      else
         p=1.0-MathPow(y,1.0/dn);
     }
   else
     {
      dk=k+1;
      //--- function call
      p=CIncBetaF::IncompleteBeta(dn,dk,0.5);
      //--- check
      if(p>0.5)
         p=CIncBetaF::InvIncompleteBeta(dk,dn,1.0-y);
      else
         p=1.0-CIncBetaF::InvIncompleteBeta(dn,dk,y);
     }
//--- return result
   return(p);
  }
//+------------------------------------------------------------------+
//| Chebyshev polynomials                                            |
//+------------------------------------------------------------------+
class CChebyshev
  {
public:
   //--- constructor, destructor
                     CChebyshev(void);
                    ~CChebyshev(void);
   //--- methods
   static double     ChebyshevCalculate(const int r,const int n,const double x);
   static double     ChebyshevSum(double &c[],const int r,const int n,const double x);
   static void       ChebyshevCoefficients(const int n,double &c[]);
   static void       FromChebyshev(double &a[],const int n,double &b[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CChebyshev::CChebyshev(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CChebyshev::~CChebyshev(void)
  {

  }
//+------------------------------------------------------------------+
//| Calculation of the value of the Chebyshev polynomials of the     |
//| first and second kinds.                                          |
//| Parameters:                                                      |
//|     r   -   polynomial kind, either 1 or 2.                      |
//|     n   -   degree, n>=0                                         |
//|     x   -   argument, -1 <= x <= 1                               |
//| Result:                                                          |
//|     the value of the Chebyshev polynomial at x                   |
//+------------------------------------------------------------------+
static double CChebyshev::ChebyshevCalculate(const int r,const int n,
                                             const double x)
  {
//--- create variables
   double result=0;
   int    i=0;
   double a=0;
   double b=0;
//--- Prepare A and B
   if(r==1)
     {
      a=1;
      b=x;
     }
   else
     {
      a=1;
      b=2*x;
     }
//--- Special cases: N=0 or N=1
   if(n==0)
      return(a);
//--- check
   if(n==1)
      return(b);
//--- General case: N>=2
   for(i=2;i<=n;i++)
     {
      result=2*x*b-a;
      a=b;
      b=result;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Summation of Chebyshev polynomials using Clenshaw?s recurrence   |
//| formula.                                                         |
//| This routine calculates                                          |
//|     c[0]*T0(x) + c[1]*T1(x) + ... + c[N]*TN(x)                   |
//| or                                                               |
//|     c[0]*U0(x) + c[1]*U1(x) + ... + c[N]*UN(x)                   |
//| depending on the R.                                              |
//| Parameters:                                                      |
//|     r   -   polynomial kind, either 1 or 2.                      |
//|     n   -   degree, n>=0                                         |
//|     x   -   argument                                             |
//| Result:                                                          |
//|     the value of the Chebyshev polynomial at x                   |
//+------------------------------------------------------------------+
static double CChebyshev::ChebyshevSum(double &c[],const int r,const int n,
                                       const double x)
  {
//--- create variables
   double result=0;
   double b1=0;
   double b2=0;
   int    i=0;
//--- initialization
   b1=0;
   b2=0;
//--- calculation
   for(i=n;i>=1;i--)
     {
      result=2*x*b1-b2+c[i];
      b2=b1;
      b1=result;
     }
//--- check
   if(r==1)
      result=-b2+x*b1+c[0];
   else
      result=-b2+2*x*b1+c[0];
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Representation of Tn as C[0] + C[1]*X + ... + C[N]*X^N           |
//| Input parameters:                                                |
//|     N   -   polynomial degree, n>=0                              |
//| Output parameters:                                               |
//|     C   -   coefficients                                         |
//+------------------------------------------------------------------+
static void CChebyshev::ChebyshevCoefficients(const int n,double &c[])
  {
//--- create a variable
   int i=0;
//--- allocation
   ArrayResizeAL(c,n+1);
   for(i=0;i<=n;i++)
      c[i]=0;
//--- check
   if(n==0 || n==1)
      c[n]=1;
   else
     {
      c[n]=MathExp((n-1)*MathLog(2));
      //--- calculation
      for(i=0;i<=n/2-1;i++)
         c[n-2*(i+1)]=-(c[n-2*i]*(n-2*i)*(n-2*i-1)/4/(i+1)/(n-i-1));
     }
  }
//+------------------------------------------------------------------+
//| Conversion of a series of Chebyshev polynomials to a power       |
//| series.                                                          |
//| Represents A[0]*T0(x) + A[1]*T1(x) + ... + A[N]*Tn(x) as         |
//| B[0] + B[1]*X + ... + B[N]*X^N.                                  |
//| Input parameters:                                                |
//|     A   -   Chebyshev series coefficients                        |
//|     N   -   degree, N>=0                                         |
//| Output parameters                                                |
//|     B   -   power series coefficients                            |
//+------------------------------------------------------------------+
static void CChebyshev::FromChebyshev(double &a[],const int n,double &b[])
  {
//--- create variables
   int    i=0;
   int    k=0;
   double e=0;
   double d=0;
//--- allocation
   ArrayResizeAL(b,n+1);
//--- initialization
   for(i=0;i<=n;i++)
      b[i]=0;
//--- calculation
   d=0;
   i=0;
   do
     {
      k=i;
      //--- cycle
      do
        {
         e=b[k];
         b[k]=0;
         //--- check
         if(i<=1 && k==i)
            b[k]=1;
         else
           {
            //--- check
            if(i!=0)
               b[k]=2*d;
            //--- check
            if(k>i+1)
               b[k]=b[k]-b[k-2];
           }
         //--- change values
         d=e;
         k=k+1;
        }
      while(k<=n);
      //--- change values
      d=b[i];
      e=0;
      k=i;
      //--- calculation
      while(k<=n)
        {
         e=e+b[k]*a[k];
         k=k+2;
        }
      //--- change values
      b[i]=e;
      i=i+1;
     }
   while(i<=n);
  }
//+------------------------------------------------------------------+
//| Chi-square distribution                                          |
//+------------------------------------------------------------------+
class CChiSquareDistr
  {
public:
   //--- constructor, destructor
                     CChiSquareDistr(void);
                    ~CChiSquareDistr(void);
   //--- methods
   static double     ChiSquareDistribution(const double v,const double x);
   static double     ChiSquareComplDistribution(const double v,const double x);
   static double     InvChiSquareDistribution(const double v,const double y);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CChiSquareDistr::CChiSquareDistr(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CChiSquareDistr::~CChiSquareDistr(void)
  {

  }
//+------------------------------------------------------------------+
//| Chi-square distribution                                          |
//| Returns the area under the left hand tail (from 0 to x)          |
//| of the Chi square probability density function with              |
//| v degrees of freedom.                                            |
//|                                   x                              |
//|                                    -                             |
//|                        1          | |  v/2-1  -t/2               |
//|  P( x | v )   =   -----------     |   t      e     dt            |
//|                    v/2  -       | |                              |
//|                   2    | (v/2)   -                               |
//|                                   0                              |
//| where x is the Chi-square variable.                              |
//| The incomplete gamma integral is used, according to the          |
//| formula                                                          |
//| y = chdtr( v, x ) = igam( v/2.0, x/2.0 ).                        |
//| The arguments must both be positive.                             |
//| ACCURACY:                                                        |
//| See incomplete gamma function                                    |
//+------------------------------------------------------------------+
static double CChiSquareDistr::ChiSquareDistribution(const double v,
                                                     const double x)
  {
//--- check
   if(!CAp::Assert(x>=0 && v>=1,__FUNCTION__+": the error variable"))
      return(EMPTY_VALUE);
//--- return result
   return(CIncGammaF::IncompleteGamma(v/2.0,x/2.0));
  }
//+------------------------------------------------------------------+
//| Complemented Chi-square distribution                             |
//| Returns the area under the right hand tail (from x to            |
//| infinity) of the Chi square probability density function         |
//| with v degrees of freedom:                                       |
//|                                  inf.                            |
//|                                    -                             |
//|                        1          | |  v/2-1  -t/2               |
//|  P( x | v )   =   -----------     |   t      e     dt            |
//|                    v/2  -       | |                              |
//|                   2    | (v/2)   -                               |
//|                                   x                              |
//| where x is the Chi-square variable.                              |
//| The incomplete gamma integral is used, according to the          |
//| formula                                                          |
//| y = chdtr( v, x ) = igamc( v/2.0, x/2.0 ).                       |
//| The arguments must both be positive.                             |
//| ACCURACY:                                                        |
//| See incomplete gamma function                                    |
//+------------------------------------------------------------------+
static double CChiSquareDistr::ChiSquareComplDistribution(const double v,
                                                          const double x)
  {
//--- check
   if(!CAp::Assert(x>=0.0 && v>=1.0,__FUNCTION__+": domain error"))
      return(EMPTY_VALUE);
//--- return result
   return(CIncGammaF::IncompleteGammaC(v/2.0,x/2.0));
  }
//+------------------------------------------------------------------+
//| Inverse of complemented Chi-square distribution                  |
//| Finds the Chi-square argument x such that the integral           |
//| from x to infinity of the Chi-square density is equal            |
//| to the given cumulative probability y.                           |
//| This is accomplished using the inverse gamma integral            |
//| function and the relation                                        |
//|    x/2 = igami( df/2, y );                                       |
//| ACCURACY:                                                        |
//| See inverse incomplete gamma function                            |
//+------------------------------------------------------------------+
static double CChiSquareDistr::InvChiSquareDistribution(const double v,
                                                        const double y)
  {
//--- check
   if(!CAp::Assert((y>=0.0 && y<=1.0) && v>=1.0,__FUNCTION__+": domain error"))
      return(EMPTY_VALUE);
//--- return result
   return(2*CIncGammaF::InvIncompleteGammaC(0.5*v,y));
  }
//+------------------------------------------------------------------+
//| Dawson integral                                                  |
//+------------------------------------------------------------------+
class CDawson
  {
public:
   //--- constructor, destructor
                     CDawson(void);
                    ~CDawson(void);
   //--- method
   static double     DawsonIntegral(double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDawson::CDawson(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CDawson::~CDawson(void)
  {

  }
//+------------------------------------------------------------------+
//| Dawson's Integral                                                |
//| Approximates the integral                                        |
//|                             x                                    |
//|                             -                                    |
//|                      2     | |        2                          |
//|  dawsn(x)  =  exp( -x  )   |    exp( t  ) dt                     |
//|                          | |                                     |
//|                           -                                      |
//|                           0                                      |
//| Three different rational approximations are employed, for        |
//| the intervals 0 to 3.25; 3.25 to 6.25; and 6.25 up.              |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0,10        10000       6.9e-16     1.0e-16         |
//+------------------------------------------------------------------+
static double CDawson::DawsonIntegral(double x)
  {
//--- create variables
   double result=0;
   double x2=0;
   double y=0;
   int    sg=0;
   double an=0;
   double ad=0;
   double bn=0;
   double bd=0;
   double cn=0;
   double cd=0;
//--- initialization
   sg=1;
//--- check
   if(x<0.0)
     {
      sg=-1;
      x=-x;
     }
//--- check
   if(x<3.25)
     {
      //--- calculation
      x2=x*x;
      an=1.13681498971755972054E-11;
      an=an*x2+8.49262267667473811108E-10;
      an=an*x2+1.94434204175553054283E-8;
      an=an*x2+9.53151741254484363489E-7;
      an=an*x2+3.07828309874913200438E-6;
      an=an*x2+3.52513368520288738649E-4;
      an=an*x2+-8.50149846724410912031E-4;
      an=an*x2+4.22618223005546594270E-2;
      an=an*x2+-9.17480371773452345351E-2;
      an=an*x2+9.99999999999999994612E-1;
      ad=2.40372073066762605484E-11;
      ad=ad*x2+1.48864681368493396752E-9;
      ad=ad*x2+5.21265281010541664570E-8;
      ad=ad*x2+1.27258478273186970203E-6;
      ad=ad*x2+2.32490249820789513991E-5;
      ad=ad*x2+3.25524741826057911661E-4;
      ad=ad*x2+3.48805814657162590916E-3;
      ad=ad*x2+2.79448531198828973716E-2;
      ad=ad*x2+1.58874241960120565368E-1;
      ad=ad*x2+5.74918629489320327824E-1;
      ad=ad*x2+1.00000000000000000539E0;
      y=x*an/ad;
      //--- return result
      return(sg*y);
     }
   x2=1.0/(x*x);
//--- check
   if(x<6.25)
     {
      //--- calculation
      bn=5.08955156417900903354E-1;
      bn=bn*x2-2.44754418142697847934E-1;
      bn=bn*x2+9.41512335303534411857E-2;
      bn=bn*x2-2.18711255142039025206E-2;
      bn=bn*x2+3.66207612329569181322E-3;
      bn=bn*x2-4.23209114460388756528E-4;
      bn=bn*x2+3.59641304793896631888E-5;
      bn=bn*x2-2.14640351719968974225E-6;
      bn=bn*x2+9.10010780076391431042E-8;
      bn=bn*x2-2.40274520828250956942E-9;
      bn=bn*x2+3.59233385440928410398E-11;
      bd=1.00000000000000000000E0;
      bd=bd*x2-6.31839869873368190192E-1;
      bd=bd*x2+2.36706788228248691528E-1;
      bd=bd*x2-5.31806367003223277662E-2;
      bd=bd*x2+8.48041718586295374409E-3;
      bd=bd*x2-9.47996768486665330168E-4;
      bd=bd*x2+7.81025592944552338085E-5;
      bd=bd*x2-4.55875153252442634831E-6;
      bd=bd*x2+1.89100358111421846170E-7;
      bd=bd*x2-4.91324691331920606875E-9;
      bd=bd*x2+7.18466403235734541950E-11;
      y=1.0/x+x2*bn/(bd*x);
      //--- return result
      return(sg*0.5*y);
     }
//--- check
   if(x>1.0E9)
      return(sg*0.5/x);
//--- calculation
   cn=-5.90592860534773254987E-1;
   cn=cn*x2+6.29235242724368800674E-1;
   cn=cn*x2-1.72858975380388136411E-1;
   cn=cn*x2+1.64837047825189632310E-2;
   cn=cn*x2-4.86827613020462700845E-4;
   cd=1.00000000000000000000E0;
   cd=cd*x2-2.69820057197544900361E0;
   cd=cd*x2+1.73270799045947845857E0;
   cd=cd*x2-3.93708582281939493482E-1;
   cd=cd*x2+3.44278924041233391079E-2;
   cd=cd*x2-9.73655226040941223894E-4;
   y=1.0/x+x2*cn/(cd*x);
//--- return result
   return(sg*0.5*y);
  }
//+------------------------------------------------------------------+
//| Elliptic integral                                                |
//+------------------------------------------------------------------+
class CElliptic
  {
public:
   //--- constructor, destructor
                     CElliptic(void);
                    ~CElliptic(void);
   //--- methods
   static double     EllipticIntegralK(const double m);
   static double     EllipticIntegralKhighPrecision(const double m1);
   static double     IncompleteEllipticIntegralK(double phi,const double m);
   static double     EllipticIntegralE(double m);
   static double     IncompleteEllipticIntegralE(const double phi,const double m);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CElliptic::CElliptic(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CElliptic::~CElliptic(void)
  {

  }
//+------------------------------------------------------------------+
//| Complete elliptic integral of the first kind                     |
//| Approximates the integral                                        |
//|            pi/2                                                  |
//|             -                                                    |
//|            | |                                                   |
//|            |           dt                                        |
//| K(m)  =    |    ------------------                               |
//|            |                   2                                 |
//|          | |    sqrt( 1 - m sin t )                              |
//|           -                                                      |
//|            0                                                     |
//| using the approximation                                          |
//|     P(x)  -  log x Q(x).                                         |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE       0,1        30000       2.5e-16     6.8e-17         |
//+------------------------------------------------------------------+
static double CElliptic::EllipticIntegralK(const double m)
  {
//--- return result
   return(EllipticIntegralKhighPrecision(1.0-m));
  }
//+------------------------------------------------------------------+
//| Complete elliptic integral of the first kind                     |
//| Approximates the integral                                        |
//|            pi/2                                                  |
//|             -                                                    |
//|            | |                                                   |
//|            |           dt                                        |
//| K(m)  =    |    ------------------                               |
//|            |                   2                                 |
//|          | |    sqrt( 1 - m sin t )                              |
//|           -                                                      |
//|            0                                                     |
//| where m = 1 - m1, using the approximation                        |
//|     P(x)  -  log x Q(x).                                         |
//| The argument m1 is used rather than m so that the logarithmic    |
//| singularity at m = 1 will be shifted to the origin; this         |
//| preserves maximum accuracy.                                      |
//| K(0) = pi/2.                                                     |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE       0,1        30000       2.5e-16     6.8e-17         |
//+------------------------------------------------------------------+
static double CElliptic::EllipticIntegralKhighPrecision(const double m1)
  {
//--- create variables
   double result=0;
   double p=0;
   double q=0;
//--- check
   if(m1<=CMath::m_machineepsilon)
      result=1.3862943611198906188E0-0.5*MathLog(m1);
   else
     {
      //--- calculation
      p=1.37982864606273237150E-4;
      p=p*m1+2.28025724005875567385E-3;
      p=p*m1+7.97404013220415179367E-3;
      p=p*m1+9.85821379021226008714E-3;
      p=p*m1+6.87489687449949877925E-3;
      p=p*m1+6.18901033637687613229E-3;
      p=p*m1+8.79078273952743772254E-3;
      p=p*m1+1.49380448916805252718E-2;
      p=p*m1+3.08851465246711995998E-2;
      p=p*m1+9.65735902811690126535E-2;
      p=p*m1+1.38629436111989062502E0;
      q=2.94078955048598507511E-5;
      q=q*m1+9.14184723865917226571E-4;
      q=q*m1+5.94058303753167793257E-3;
      q=q*m1+1.54850516649762399335E-2;
      q=q*m1+2.39089602715924892727E-2;
      q=q*m1+3.01204715227604046988E-2;
      q=q*m1+3.73774314173823228969E-2;
      q=q*m1+4.88280347570998239232E-2;
      q=q*m1+7.03124996963957469739E-2;
      q=q*m1+1.24999999999870820058E-1;
      q=q*m1+4.99999999999999999821E-1;
      //--- get result
      result=p-q*MathLog(m1);
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Incomplete elliptic integral of the first kind F(phi|m)          |
//| Approximates the integral                                        |
//|                phi                                               |
//|                 -                                                |
//|                | |                                               |
//|                |           dt                                    |
//| F(phi_\m)  =   |    ------------------                           |
//|                |                   2                             |
//|              | |    sqrt( 1 - m sin t )                          |
//|               -                                                  |
//|                0                                                 |
//| of amplitude phi and modulus m, using the arithmetic -           |
//| geometric mean algorithm.                                        |
//| ACCURACY:                                                        |
//| Tested at random points with m in [0, 1] and phi as indicated.   |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE     -10,10       200000      7.4e-16     1.0e-16         |
//+------------------------------------------------------------------+
static double CElliptic::IncompleteEllipticIntegralK(double phi,const double m)
  {
//--- create variables
   double a=0;
   double b=0;
   double c=0;
   double e=0;
   double temp=0;
   double pio2=0;
   double t=0;
   double k=0;
   int    d=0;
   int    md=0;
   int    s=0;
   int    npio2=0;
//--- initialization
   pio2=1.57079632679489661923;
//--- check
   if(m==0.0)
      return(phi);
//--- initialization
   a=1-m;
//--- check
   if(a==0.0)
      return(MathLog(MathTan(0.5*(pio2+phi))));
//--- initialization
   npio2=(int)MathFloor(phi/pio2);
//--- check
   if(npio2%2!=0)
      npio2=npio2+1;
//--- check
   if(npio2!=0)
     {
      k=EllipticIntegralK(1-a);
      phi=phi-npio2*pio2;
     }
   else
      k=0;
//--- check
   if(phi<0.0)
     {
      phi=-phi;
      s=-1;
     }
   else
      s=0;
//--- calculation
   b=MathSqrt(a);
   t=MathTan(phi);
//--- check
   if(MathAbs(t)>10)
     {
      e=1.0/(b*t);
      //--- check
      if(MathAbs(e)<10)
        {
         e=MathArctan(e);
         //--- check
         if(npio2==0)
            k=EllipticIntegralK(1-a);
         //--- function call
         temp=k-IncompleteEllipticIntegralK(e,m);
         //--- check
         if(s<0)
            temp=-temp;
         //--- return result
         return(temp+npio2*k);
        }
     }
//--- change values
   a=1.0;
   c=MathSqrt(m);
   d=1;
   md=0;
//--- cycle
   while(MathAbs(c/a)>CMath::m_machineepsilon)
     {
      //--- calculation
      temp=b/a;
      phi=phi+MathArctan(t*temp)+md*M_PI;
      md=(int)((phi+pio2)/M_PI);
      t=t*(1.0+temp)/(1.0-temp*t*t);
      c=0.5*(a-b);
      temp=MathSqrt(a*b);
      a=0.5*(a+b);
      b=temp;
      d=d+d;
     }
//--- change value
   temp=(MathArctan(t)+md*M_PI)/(d*a);
//--- check
   if(s<0)
      temp=-temp;
//--- return result
   return(temp+npio2*k);
  }
//+------------------------------------------------------------------+
//| Complete elliptic integral of the second kind                    |
//| Approximates the integral                                        |
//|            pi/2                                                  |
//|             -                                                    |
//|            | |                 2                                 |
//| E(m)  =    |    sqrt( 1 - m sin t ) dt                           |
//|          | |                                                     |
//|           -                                                      |
//|            0                                                     |
//| using the approximation                                          |
//|      P(x)  -  x log x Q(x).                                      |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE       0, 1       10000       2.1e-16     7.3e-17         |
//+------------------------------------------------------------------+
static double CElliptic::EllipticIntegralE(double m)
  {
//--- create variables
   double p=0;
   double q=0;
//--- check
   if(!CAp::Assert(m>=0.0 && m<=1.0,__FUNCTION__+": m<0 or m>1"))
      return(EMPTY_VALUE);
//--- change value
   m=1-m;
//--- check
   if(m==0.0)
      return(1);
//--- calculation
   p=1.53552577301013293365E-4;
   p=p*m+2.50888492163602060990E-3;
   p=p*m+8.68786816565889628429E-3;
   p=p*m+1.07350949056076193403E-2;
   p=p*m+7.77395492516787092951E-3;
   p=p*m+7.58395289413514708519E-3;
   p=p*m+1.15688436810574127319E-2;
   p=p*m+2.18317996015557253103E-2;
   p=p*m+5.68051945617860553470E-2;
   p=p*m+4.43147180560990850618E-1;
   p=p*m+1.00000000000000000299E0;
   q=3.27954898576485872656E-5;
   q=q*m+1.00962792679356715133E-3;
   q=q*m+6.50609489976927491433E-3;
   q=q*m+1.68862163993311317300E-2;
   q=q*m+2.61769742454493659583E-2;
   q=q*m+3.34833904888224918614E-2;
   q=q*m+4.27180926518931511717E-2;
   q=q*m+5.85936634471101055642E-2;
   q=q*m+9.37499997197644278445E-2;
   q=q*m+2.49999999999888314361E-1;
//--- return result
   return(p-q*m*MathLog(m));
  }
//+------------------------------------------------------------------+
//| Incomplete elliptic integral of the second kind                  |
//| Approximates the integral                                        |
//|                phi                                               |
//|                 -                                                |
//|                | |                                               |
//|                |                   2                             |
//| E(phi_\m)  =   |    sqrt( 1 - m sin t ) dt                       |
//|                |                                                 |
//|              | |                                                 |
//|               -                                                  |
//|                0                                                 |
//| of amplitude phi and modulus m, using the arithmetic -           |
//| geometric mean algorithm.                                        |
//| ACCURACY:                                                        |
//| Tested at random arguments with phi in [-10, 10] and m in        |
//| [0, 1].                                                          |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE     -10,10      150000       3.3e-15     1.4e-16         |
//+------------------------------------------------------------------+
static double CElliptic::IncompleteEllipticIntegralE(const double phi,
                                                     const double m)
  {
//--- create variables
   double pio2=0;
   double a=0;
   double b=0;
   double c=0;
   double e=0;
   double temp=0;
   double lphi=0;
   double t=0;
   double ebig=0;
   int d=0;
   int md=0;
   int npio2=0;
   int s=0;
//--- initialization
   pio2=1.57079632679489661923;
//--- check
   if(m==0.0)
      return(phi);
//--- change values
   lphi=phi;
   npio2=(int)MathFloor(lphi/pio2);
//--- check
   if(npio2%2!=0)
      npio2=npio2+1;
   lphi=lphi-npio2*pio2;
//--- check
   if(lphi<0.0)
     {
      lphi=-lphi;
      s=-1;
     }
   else
      s=1;
//--- change values
   a=1.0-m;
   ebig=EllipticIntegralE(m);
//--- check
   if(a==0.0)
     {
      temp=MathSin(lphi);
      //--- check
      if(s<0)
         temp=-temp;
      //--- return result
      return(temp+npio2*ebig);
     }
//--- calculation
   t=MathTan(lphi);
   b=MathSqrt(a);
//--- Thanks to Brian Fitzgerald <fitzgb@mml0.meche.rpi.edu>
//--- for pointing out an instability near odd multiples of pi/2
   if(MathAbs(t)>10)
     {
      //--- Transform the amplitude
      e=1.0/(b*t);
      //--- ... but avoid multiple recursions.
      if(MathAbs(e)<10)
        {
         e=MathArctan(e);
         temp=ebig+m*MathSin(lphi)*MathSin(e)-IncompleteEllipticIntegralE(e,m);
         //--- check
         if(s<0)
            temp=-temp;
         //--- return result
         return(temp+npio2*ebig);
        }
     }
//--- change values
   c=MathSqrt(m);
   a=1.0;
   d=1;
   e=0.0;
   md=0;
//--- cycle
   while(MathAbs(c/a)>CMath::m_machineepsilon)
     {
      //--- calculation
      temp=b/a;
      lphi=lphi+MathArctan(t*temp)+md*M_PI;
      md=(int)((lphi+pio2)/M_PI);
      t=t*(1.0+temp)/(1.0-temp*t*t);
      c=0.5*(a-b);
      temp=MathSqrt(a*b);
      a=0.5*(a+b);
      b=temp;
      d=d+d;
      e=e+c*MathSin(lphi);
     }
//--- change value
   temp=ebig/EllipticIntegralK(m);
   temp=temp*((MathArctan(t)+md*M_PI)/(d*a));
   temp=temp+e;
//--- check
   if(s<0)
      temp=-temp;
//--- return result
   return(temp+npio2*ebig);
  }
//+------------------------------------------------------------------+
//| Exponential integral                                             |
//+------------------------------------------------------------------+
class CExpIntegrals
  {
public:
   //--- constructor, destructor
                     CExpIntegrals(void);
                    ~CExpIntegrals(void);
   //--- methods
   static double     ExponentialIntegralEi(const double x);
   static double     ExponentialIntegralEn(const double x,const int n);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CExpIntegrals::CExpIntegrals(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CExpIntegrals::~CExpIntegrals(void)
  {

  }
//+------------------------------------------------------------------+
//| Exponential integral Ei(x)                                       |
//|               x                                                  |
//|                -     t                                           |
//|               | |   e                                            |
//|    Ei(x) =   -|-   ---  dt .                                     |
//|             | |     t                                            |
//|              -                                                   |
//|             -inf                                                 |
//| Not defined for x <= 0.                                          |
//| See also expn.c.                                                 |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE       0,100       50000      8.6e-16     1.3e-16         |
//+------------------------------------------------------------------+
static double CExpIntegrals::ExponentialIntegralEi(const double x)
  {
//--- create variables
   double eul=0;
   double f=0;
   double f1=0;
   double f2=0;
   double w=0;
//--- initialization
   eul=0.5772156649015328606065;
//--- check
   if(x<=0.0)
      return(0);
//--- check
   if(x<2.0)
     {
      //--- calculation
      f1=-5.350447357812542947283;
      f1=f1*x+218.5049168816613393830;
      f1=f1*x-4176.572384826693777058;
      f1=f1*x+55411.76756393557601232;
      f1=f1*x-331338.1331178144034309;
      f1=f1*x+1592627.163384945414220;
      f2=1.000000000000000000000;
      f2=f2*x-52.50547959112862969197;
      f2=f2*x+1259.616186786790571525;
      f2=f2*x-17565.49581973534652631;
      f2=f2*x+149306.2117002725991967;
      f2=f2*x-729494.9239640527645655;
      f2=f2*x+1592627.163384945429726;
      f=f1/f2;
      //--- return result
      return(eul+MathLog(x)+x*f);
     }
//--- check
   if(x<4.0)
     {
      //--- calculation
      w=1/x;
      f1=1.981808503259689673238E-2;
      f1=f1*w-1.271645625984917501326;
      f1=f1*w-2.088160335681228318920;
      f1=f1*w+2.755544509187936721172;
      f1=f1*w-4.409507048701600257171E-1;
      f1=f1*w+4.665623805935891391017E-2;
      f1=f1*w-1.545042679673485262580E-3;
      f1=f1*w+7.059980605299617478514E-5;
      f2=1.000000000000000000000;
      f2=f2*w+1.476498670914921440652;
      f2=f2*w+5.629177174822436244827E-1;
      f2=f2*w+1.699017897879307263248E-1;
      f2=f2*w+2.291647179034212017463E-2;
      f2=f2*w+4.450150439728752875043E-3;
      f2=f2*w+1.727439612206521482874E-4;
      f2=f2*w+3.953167195549672482304E-5;
      f=f1/f2;
      //--- return result
      return(MathExp(x)*w*(1+w*f));
     }
//--- check
   if(x<8.0)
     {
      //--- calculation
      w=1/x;
      f1=-1.373215375871208729803;
      f1=f1*w-7.084559133740838761406E-1;
      f1=f1*w+1.580806855547941010501;
      f1=f1*w-2.601500427425622944234E-1;
      f1=f1*w+2.994674694113713763365E-2;
      f1=f1*w-1.038086040188744005513E-3;
      f1=f1*w+4.371064420753005429514E-5;
      f1=f1*w+2.141783679522602903795E-6;
      f2=1.000000000000000000000;
      f2=f2*w+8.585231423622028380768E-1;
      f2=f2*w+4.483285822873995129957E-1;
      f2=f2*w+7.687932158124475434091E-2;
      f2=f2*w+2.449868241021887685904E-2;
      f2=f2*w+8.832165941927796567926E-4;
      f2=f2*w+4.590952299511353531215E-4;
      f2=f2*w+-4.729848351866523044863E-6;
      f2=f2*w+2.665195537390710170105E-6;
      f=f1/f2;
      //--- return result
      return(MathExp(x)*w*(1+w*f));
     }
//--- check
   if(x<16.0)
     {
      //--- calculation
      w=1/x;
      f1=-2.106934601691916512584;
      f1=f1*w+1.732733869664688041885;
      f1=f1*w-2.423619178935841904839E-1;
      f1=f1*w+2.322724180937565842585E-2;
      f1=f1*w+2.372880440493179832059E-4;
      f1=f1*w-8.343219561192552752335E-5;
      f1=f1*w+1.363408795605250394881E-5;
      f1=f1*w-3.655412321999253963714E-7;
      f1=f1*w+1.464941733975961318456E-8;
      f1=f1*w+6.176407863710360207074E-10;
      f2=1.000000000000000000000;
      f2=f2*w-2.298062239901678075778E-1;
      f2=f2*w+1.105077041474037862347E-1;
      f2=f2*w-1.566542966630792353556E-2;
      f2=f2*w+2.761106850817352773874E-3;
      f2=f2*w-2.089148012284048449115E-4;
      f2=f2*w+1.708528938807675304186E-5;
      f2=f2*w-4.459311796356686423199E-7;
      f2=f2*w+1.394634930353847498145E-8;
      f2=f2*w+6.150865933977338354138E-10;
      f=f1/f2;
      //--- return result
      return(MathExp(x)*w*(1+w*f));
     }
//--- check
   if(x<32.0)
     {
      //--- calculation
      w=1/x;
      f1=-2.458119367674020323359E-1;
      f1=f1*w-1.483382253322077687183E-1;
      f1=f1*w+7.248291795735551591813E-2;
      f1=f1*w-1.348315687380940523823E-2;
      f1=f1*w+1.342775069788636972294E-3;
      f1=f1*w-7.942465637159712264564E-5;
      f1=f1*w+2.644179518984235952241E-6;
      f1=f1*w-4.239473659313765177195E-8;
      f2=1.000000000000000000000;
      f2=f2*w-1.044225908443871106315E-1;
      f2=f2*w-2.676453128101402655055E-1;
      f2=f2*w+9.695000254621984627876E-2;
      f2=f2*w-1.601745692712991078208E-2;
      f2=f2*w+1.496414899205908021882E-3;
      f2=f2*w-8.462452563778485013756E-5;
      f2=f2*w+2.728938403476726394024E-6;
      f2=f2*w-4.239462431819542051337E-8;
      f=f1/f2;
      //--- return result
      return(MathExp(x)*w*(1+w*f));
     }
//--- check
   if(x<64)
     {
      //--- calculation
      w=1/x;
      f1=1.212561118105456670844E-1;
      f1=f1*w-5.823133179043894485122E-1;
      f1=f1*w+2.348887314557016779211E-1;
      f1=f1*w-3.040034318113248237280E-2;
      f1=f1*w+1.510082146865190661777E-3;
      f1=f1*w-2.523137095499571377122E-5;
      f2=1.000000000000000000000;
      f2=f2*w-1.002252150365854016662;
      f2=f2*w+2.928709694872224144953E-1;
      f2=f2*w-3.337004338674007801307E-2;
      f2=f2*w+1.560544881127388842819E-3;
      f2=f2*w-2.523137093603234562648E-5;
      f=f1/f2;
      //--- return result
      return(MathExp(x)*w*(1+w*f));
     }
//--- calculation
   w=1/x;
   f1=-7.657847078286127362028E-1;
   f1=f1*w+6.886192415566705051750E-1;
   f1=f1*w-2.132598113545206124553E-1;
   f1=f1*w+3.346107552384193813594E-2;
   f1=f1*w-3.076541477344756050249E-3;
   f1=f1*w+1.747119316454907477380E-4;
   f1=f1*w-6.103711682274170530369E-6;
   f1=f1*w+1.218032765428652199087E-7;
   f1=f1*w-1.086076102793290233007E-9;
   f2=1.000000000000000000000;
   f2=f2*w-1.888802868662308731041;
   f2=f2*w+1.066691687211408896850;
   f2=f2*w-2.751915982306380647738E-1;
   f2=f2*w+3.930852688233823569726E-2;
   f2=f2*w-3.414684558602365085394E-3;
   f2=f2*w+1.866844370703555398195E-4;
   f2=f2*w-6.345146083130515357861E-6;
   f2=f2*w+1.239754287483206878024E-7;
   f2=f2*w-1.086076102793126632978E-9;
   f=f1/f2;
//--- return result
   return(MathExp(x)*w*(1+w*f));
  }
//+------------------------------------------------------------------+
//| Exponential integral En(x)                                       |
//| Evaluates the exponential integral                               |
//|                 inf.                                             |
//|                   -                                              |
//|                  | |   -xt                                       |
//|                  |    e                                          |
//|      E (x)  =    |    ----  dt.                                  |
//|       n          |      n                                        |
//|                | |     t                                         |
//|                 -                                                |
//|                  1                                               |
//| Both n and x must be nonnegative.                                |
//| The routine employs either a power series, a continued           |
//| fraction, or an asymptotic formula depending on the              |
//| relative values of n and x.                                      |
//| ACCURACY:                                                        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0, 30       10000       1.7e-15     3.6e-16         |
//+------------------------------------------------------------------+
static double CExpIntegrals::ExponentialIntegralEn(const double x,const int n)
  {
//--- create variables
   double result=0;
   double r=0;
   double t=0;
   double yk=0;
   double xk=0;
   double pk=0;
   double pkm1=0;
   double pkm2=0;
   double qk=0;
   double qkm1=0;
   double qkm2=0;
   double psi=0;
   double z=0;
   int    i=0;
   int    k=0;
   double big=0;
   double eul=0;
//--- initialization
   eul=0.57721566490153286060;
   big=1.44115188075855872*MathPow(10,17);
//--- check
   if(((n<0 || x<0.0) || x>170) || (x==0.0 && n<2))
      return(-1);
//--- check
   if(x==0.0)
      return(1.0/(double)(n-1));
//--- check
   if(n==0)
      return(MathExp(-x)/x);
//--- check
   if(n>5000)
     {
      //--- calculation
      xk=x+n;
      yk=1/(xk*xk);
      t=n;
      result=yk*t*(6*x*x-8*t*x+t*t);
      result=yk*(result+t*(t-2.0*x));
      result=yk*(result+t);
      result=(result+1)*MathExp(-x)/xk;
      //--- return result
      return(result);
     }
//--- check
   if(x<=1.0)
     {
      //--- calculation
      psi=-eul-MathLog(x);
      for(i=1;i<=n-1;i++)
         psi=psi+1.0/(double)i;
      z=-x;
      xk=0;
      yk=1;
      pk=1-n;
      //--- check
      if(n==1)
         result=0.0;
      else
         result=1.0/pk;
      do
        {
         //--- change values
         xk=xk+1;
         yk=yk*z/xk;
         pk=pk+1;
         //--- check
         if(pk!=0.0)
            result=result+yk/pk;
         //--- check
         if(result!=0.0)
            t=MathAbs(yk/result);
         else
            t=1;
        }
      while(t>=CMath::m_machineepsilon);
      t=1;
      for(i=1;i<=n-1;i++)
         t=t*z/i;
      //--- return result
      return(psi*t-result);
     }
   else
     {
      //--- calculation
      k=1;
      pkm2=1;
      qkm2=x;
      pkm1=1.0;
      qkm1=x+n;
      result=pkm1/qkm1;
      do
        {
         k=k+1;
         //--- check
         if(k%2==1)
           {
            yk=1;
            xk=n+(double)(k-1)/2.0;
           }
         else
           {
            yk=x;
            xk=(double)k/2.0;
           }
         //--- change values
         pk=pkm1*yk+pkm2*xk;
         qk=qkm1*yk+qkm2*xk;
         //--- check
         if(qk!=0.0)
           {
            r=pk/qk;
            t=MathAbs((result-r)/r);
            result=r;
           }
         else
            t=1;
         //--- change values
         pkm2=pkm1;
         pkm1=pk;
         qkm2=qkm1;
         qkm1=qk;
         //--- check
         if(MathAbs(pk)>big)
           {
            //--- change values
            pkm2=pkm2/big;
            pkm1=pkm1/big;
            qkm2=qkm2/big;
            qkm1=qkm1/big;
           }
        }
      while(t>=CMath::m_machineepsilon);
      result=result*MathExp(-x);
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| F distribution                                                   |
//+------------------------------------------------------------------+
class CFDistr
  {
public:
   //--- constructor, destructor
                     CFDistr(void);
                    ~CFDistr(void);
   //--- methods
   static double     FDistribution(const int a,const int b,const double x);
   static double     FComplDistribution(const int a,const int b,const double x);
   static double     InvFDistribution(const int a,const int b,const double y);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CFDistr::CFDistr(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CFDistr::~CFDistr(void)
  {

  }
//+------------------------------------------------------------------+
//| F distribution                                                   |
//| Returns the area from zero to x under the F density              |
//| function (also known as Snedcor's density or the                 |
//| variance ratio density).  This is the density                    |
//| of x = (u1/df1)/(u2/df2), where u1 and u2 are random             |
//| variables having Chi square distributions with df1               |
//| and df2 degrees of freedom, respectively.                        |
//| The incomplete beta integral is used, according to the           |
//| formula                                                          |
//| P(x) = incbet( df1/2, df2/2, (df1*x/(df2 + df1*x) ).             |
//| The arguments a and b are greater than zero, and x is            |
//| nonnegative.                                                     |
//| ACCURACY:                                                        |
//| Tested at random points (a,b,x).                                 |
//|                x     a,b                     Relative error:     |
//| arithmetic  domain  domain     # trials      peak         rms    |
//|    IEEE      0,1    0,100       100000      9.8e-15     1.7e-15  |
//|    IEEE      1,5    0,100       100000      6.5e-15     3.5e-16  |
//|    IEEE      0,1    1,10000     100000      2.2e-11     3.3e-12  |
//|    IEEE      1,5    1,10000     100000      1.1e-11     1.7e-13  |
//+------------------------------------------------------------------+
static double CFDistr::FDistribution(const int a,const int b,const double x)
  {
//--- check
   if(!CAp::Assert(a>=1 && b>=1 && x>=0,__FUNCTION__+": the error variable"))
      return(EMPTY_VALUE);
//--- create a variable
   double w=0;
//--- calculation
   w=a*x;
   w/=b+w;
//--- return result
   return(CIncBetaF::IncompleteBeta(0.5*a,0.5*b,w));
  }
//+------------------------------------------------------------------+
//| Complemented F distribution                                      |
//| Returns the area from x to infinity under the F density          |
//| function (also known as Snedcor's density or the                 |
//| variance ratio density).                                         |
//|                      inf.                                        |
//|                       -                                          |
//|              1       | |  a-1      b-1                           |
//| 1-P(x)  =  ------    |   t    (1-t)    dt                        |
//|            B(a,b)  | |                                           |
//|                     -                                            |
//|                      x                                           |
//| The incomplete beta integral is used, according to the           |
//| formula                                                          |
//| P(x) = incbet( df2/2, df1/2, (df2/(df2 + df1*x) ).               |
//| ACCURACY:                                                        |
//| Tested at random points (a,b,x) in the indicated intervals.      |
//|                x     a,b                     Relative error:     |
//| arithmetic  domain  domain     # trials      peak         rms    |
//|    IEEE      0,1    1,100       100000      3.7e-14     5.9e-16  |
//|    IEEE      1,5    1,100       100000      8.0e-15     1.6e-15  |
//|    IEEE      0,1    1,10000     100000      1.8e-11     3.5e-13  |
//|    IEEE      1,5    1,10000     100000      2.0e-11     3.0e-12  |
//+------------------------------------------------------------------+
static double CFDistr::FComplDistribution(const int a,const int b,const double x)
  {
//--- create a variable
   double w=0;
//--- check
   if(!CAp::Assert((a>=1 && b>=1) && x>=0.0,__FUNCTION__+": domain error"))
      return(EMPTY_VALUE);
//--- calculation
   w=b/(b+a*x);
//--- return result
   return(CIncBetaF::IncompleteBeta(0.5*b,0.5*a,w));
  }
//+------------------------------------------------------------------+
//| Inverse of complemented F distribution                           |
//| Finds the F density argument x such that the integral            |
//| from x to infinity of the F density is equal to the              |
//| given probability p.                                             |
//| This is accomplished using the inverse beta integral             |
//| function and the relations                                       |
//|      z = incbi( df2/2, df1/2, p )                                |
//|      x = df2 (1-z) / (df1 z).                                    |
//| Note: the following relations hold for the inverse of            |
//| the uncomplemented F distribution:                               |
//|      z = incbi( df1/2, df2/2, p )                                |
//|      x = df2 z / (df1 (1-z)).                                    |
//| ACCURACY:                                                        |
//| Tested at random points (a,b,p).                                 |
//|              a,b                     Relative error:             |
//| arithmetic  domain     # trials      peak         rms            |
//|  For p between .001 and 1:                                       |
//|    IEEE     1,100       100000      8.3e-15     4.7e-16          |
//|    IEEE     1,10000     100000      2.1e-11     1.4e-13          |
//|  For p between 10^-6 and 10^-3:                                  |
//|    IEEE     1,100        50000      1.3e-12     8.4e-15          |
//|    IEEE     1,10000      50000      3.0e-12     4.8e-14          |
//+------------------------------------------------------------------+
static double CFDistr::InvFDistribution(const int a,const int b,const double y)
  {
//--- create variables
   double result=0;
   double w=0;
//--- check
   if(!CAp::Assert(((a>=1 && b>=1) && y>0.0) && y<=1.0,__FUNCTION__+": domain error"))
      return(EMPTY_VALUE);
//--- Compute probability for x=0.5
   w=CIncBetaF::IncompleteBeta(0.5*b,0.5*a,0.5);
//--- If that is greater than y,then the solution w < .5
//--- Otherwise,solve at 1-y to remove cancellation in (b - b*w)
   if(w>y || y<0.001)
     {
      //--- calculation
      w=CIncBetaF::InvIncompleteBeta(0.5*b,0.5*a,y);
      result=(b-b*w)/(a*w);
     }
   else
     {
      //--- calculation
      w=CIncBetaF::InvIncompleteBeta(0.5*a,0.5*b,1.0-y);
      result=b*w/(a*(1.0-w));
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Fresnel integral                                                 |
//+------------------------------------------------------------------+
class CFresnel
  {
public:
   //--- constructor, destructor
                     CFresnel(void);
                    ~CFresnel(void);
   //--- method
   static void       FresnelIntegral(double x,double &c,double &s);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CFresnel::CFresnel(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CFresnel::~CFresnel(void)
  {

  }
//+------------------------------------------------------------------+
//| Fresnel integral                                                 |
//| Evaluates the Fresnel integrals                                  |
//|           x                                                      |
//|           -                                                      |
//|          | |                                                     |
//| C(x) =   |   cos(pi/2 t**2) dt,                                  |
//|        | |                                                       |
//|         -                                                        |
//|          0                                                       |
//|           x                                                      |
//|           -                                                      |
//|          | |                                                     |
//| S(x) =   |   sin(pi/2 t**2) dt.                                  |
//|        | |                                                       |
//|         -                                                        |
//|          0                                                       |
//| The integrals are evaluated by a power series for x < 1.         |
//| For x >= 1 auxiliary functions f(x) and g(x) are employed        |
//| such that                                                        |
//| C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )       |
//| S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )       |
//| ACCURACY:                                                        |
//|  Relative error.                                                 |
//| Arithmetic  function   domain     # trials      peak        rms  |
//|   IEEE       S(x)      0, 10       10000       2.0e-15    3.2e-16|
//|   IEEE       C(x)      0, 10       10000       1.8e-15    3.3e-16|
//+------------------------------------------------------------------+
static void CFresnel::FresnelIntegral(double x,double &c,double &s)
  {
//--- create variables
   double xxa=0;
   double f=0;
   double g=0;
   double cc=0;
   double ss=0;
   double t=0;
   double u=0;
   double x2=0;
   double sn=0;
   double sd=0;
   double cn=0;
   double cd=0;
   double fn=0;
   double fd=0;
   double gn=0;
   double gd=0;
   double mpi=0;
   double mpio2=0;
//--- initialization
   mpi=3.14159265358979323846;
   mpio2=1.57079632679489661923;
   xxa=x;
   x=MathAbs(xxa);
   x2=x*x;
//--- check
   if(x2<2.5625)
     {
      //--- calculation
      t=x2*x2;
      sn=-2.99181919401019853726E3;
      sn=sn*t+7.08840045257738576863E5;
      sn=sn*t-6.29741486205862506537E7;
      sn=sn*t+2.54890880573376359104E9;
      sn=sn*t-4.42979518059697779103E10;
      sn=sn*t+3.18016297876567817986E11;
      sd=1.00000000000000000000E0;
      sd=sd*t+2.81376268889994315696E2;
      sd=sd*t+4.55847810806532581675E4;
      sd=sd*t+5.17343888770096400730E6;
      sd=sd*t+4.19320245898111231129E8;
      sd=sd*t+2.24411795645340920940E10;
      sd=sd*t+6.07366389490084639049E11;
      cn=-4.98843114573573548651E-8;
      cn=cn*t+9.50428062829859605134E-6;
      cn=cn*t-6.45191435683965050962E-4;
      cn=cn*t+1.88843319396703850064E-2;
      cn=cn*t-2.05525900955013891793E-1;
      cn=cn*t+9.99999999999999998822E-1;
      cd=3.99982968972495980367E-12;
      cd=cd*t+9.15439215774657478799E-10;
      cd=cd*t+1.25001862479598821474E-7;
      cd=cd*t+1.22262789024179030997E-5;
      cd=cd*t+8.68029542941784300606E-4;
      cd=cd*t+4.12142090722199792936E-2;
      cd=cd*t+1.00000000000000000118E0;
      s=MathSign(xxa)*x*x2*sn/sd;
      c=MathSign(xxa)*x*cn/cd;
      //--- exit the function
      return;
     }
//--- check
   if(x>36974.0)
     {
      //--- calculation
      c=MathSign(xxa)*0.5;
      s=MathSign(xxa)*0.5;
      //--- exit the function
      return;
     }
//--- calculation
   x2=x*x;
   t=mpi*x2;
   u=1/(t*t);
   t=1/t;
   fn=4.21543555043677546506E-1;
   fn=fn*u+1.43407919780758885261E-1;
   fn=fn*u+1.15220955073585758835E-2;
   fn=fn*u+3.45017939782574027900E-4;
   fn=fn*u+4.63613749287867322088E-6;
   fn=fn*u+3.05568983790257605827E-8;
   fn=fn*u+1.02304514164907233465E-10;
   fn=fn*u+1.72010743268161828879E-13;
   fn=fn*u+1.34283276233062758925E-16;
   fn=fn*u+3.76329711269987889006E-20;
   fd=1.00000000000000000000E0;
   fd=fd*u+7.51586398353378947175E-1;
   fd=fd*u+1.16888925859191382142E-1;
   fd=fd*u+6.44051526508858611005E-3;
   fd=fd*u+1.55934409164153020873E-4;
   fd=fd*u+1.84627567348930545870E-6;
   fd=fd*u+1.12699224763999035261E-8;
   fd=fd*u+3.60140029589371370404E-11;
   fd=fd*u+5.88754533621578410010E-14;
   fd=fd*u+4.52001434074129701496E-17;
   fd=fd*u+1.25443237090011264384E-20;
   gn=5.04442073643383265887E-1;
   gn=gn*u+1.97102833525523411709E-1;
   gn=gn*u+1.87648584092575249293E-2;
   gn=gn*u+6.84079380915393090172E-4;
   gn=gn*u+1.15138826111884280931E-5;
   gn=gn*u+9.82852443688422223854E-8;
   gn=gn*u+4.45344415861750144738E-10;
   gn=gn*u+1.08268041139020870318E-12;
   gn=gn*u+1.37555460633261799868E-15;
   gn=gn*u+8.36354435630677421531E-19;
   gn=gn*u+1.86958710162783235106E-22;
   gd=1.00000000000000000000E0;
   gd=gd*u+1.47495759925128324529E0;
   gd=gd*u+3.37748989120019970451E-1;
   gd=gd*u+2.53603741420338795122E-2;
   gd=gd*u+8.14679107184306179049E-4;
   gd=gd*u+1.27545075667729118702E-5;
   gd=gd*u+1.04314589657571990585E-7;
   gd=gd*u+4.60680728146520428211E-10;
   gd=gd*u+1.10273215066240270757E-12;
   gd=gd*u+1.38796531259578871258E-15;
   gd=gd*u+8.39158816283118707363E-19;
   gd=gd*u+1.86958710162783236342E-22;
//--- get result
   f=1-u*fn/fd;
   g=t*gn/gd;
   t=mpio2*x2;
   cc=MathCos(t);
   ss=MathSin(t);
   t=mpi*x;
   c=0.5+(f*ss-g*cc)/t;
   s=0.5-(f*cc+g*ss)/t;
   c=c*MathSign(xxa);
   s=s*MathSign(xxa);
  }
//+------------------------------------------------------------------+
//| Hermite polynomial                                               |
//+------------------------------------------------------------------+
class CHermite
  {
public:
   //--- constructor, destructor
                     CHermite(void);
                    ~CHermite(void);
   //--- methods
   static double     HermiteCalculate(const int n,const double x);
   static double     HermiteSum(double &c[],const int n,const double x);
   static void       HermiteCoefficients(const int n,double &c[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CHermite::CHermite(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CHermite::~CHermite(void)
  {

  }
//+------------------------------------------------------------------+
//| Calculation of the value of the Hermite polynomial.              |
//| Parameters:                                                      |
//|     n   -   degree, n>=0                                         |
//|     x   -   argument                                             |
//| Result:                                                          |
//|     the value of the Hermite polynomial Hn at x                  |
//+------------------------------------------------------------------+
static double CHermite::HermiteCalculate(const int n,const double x)
  {
//--- create variables
   double result=0;
   int    i=0;
   double a=0;
   double b=0;
//--- Prepare A and B
   a=1;
   b=2*x;
//--- Special cases: N=0 or N=1
   if(n==0)
      return(a);
//--- check
   if(n==1)
      return(b);
//--- General case: N>=2
   for(i=2;i<=n;i++)
     {
      result=2*x*b-2*(i-1)*a;
      a=b;
      b=result;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Summation of Hermite polynomials using Clenshaw?s recurrence     |
//| formula.                                                         |
//| This routine calculates                                          |
//|     c[0]*H0(x) + c[1]*H1(x) + ... + c[N]*HN(x)                   |
//| Parameters:                                                      |
//|     n   -   degree, n>=0                                         |
//|     x   -   argument                                             |
//| Result:                                                          |
//|     the value of the Hermite polynomial at x                     |
//+------------------------------------------------------------------+
static double CHermite::HermiteSum(double &c[],const int n,const double x)
  {
//--- create variables
   double result=0;
   double b1=0;
   double b2=0;
   int    i=0;
//--- initialization
   b1=0;
   b2=0;
//--- calculation
   for(i=n;i>=0;i--)
     {
      result=2*(x*b1-(i+1)*b2)+c[i];
      b2=b1;
      b1=result;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Representation of Hn as C[0] + C[1]*X + ... + C[N]*X^N           |
//| Input parameters:                                                |
//|     N   -   polynomial degree, n>=0                              |
//| Output parameters:                                               |
//|     C   -   coefficients                                         |
//+------------------------------------------------------------------+
static void CHermite::HermiteCoefficients(const int n,double &c[])
  {
//--- create a variable
   int i=0;
//--- allocation
   ArrayResizeAL(c,n+1);
//--- initialization
   for(i=0;i<=n;i++)
      c[i]=0;
//--- calculation
   c[n]=MathExp(n*MathLog(2));
   for(i=0;i<=n/2-1;i++)
      c[n-2*(i+1)]=-(c[n-2*i]*(n-2*i)*(n-2*i-1)/4/(i+1));
  }
//+------------------------------------------------------------------+
//| Jacobian elliptic functions                                      |
//+------------------------------------------------------------------+
class CJacobianElliptic
  {
public:
   //--- constructor, destructor
                     CJacobianElliptic(void);
                    ~CJacobianElliptic(void);
   //--- method
   static void       JacobianEllipticFunctions(const double u,const double m,double &sn,double &cn,double &dn,double &ph);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CJacobianElliptic::CJacobianElliptic(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CJacobianElliptic::~CJacobianElliptic(void)
  {

  }
//+------------------------------------------------------------------+
//| Jacobian Elliptic Functions                                      |
//| Evaluates the Jacobian elliptic functions sn(u|m), cn(u|m),      |
//| and dn(u|m) of parameter m between 0 and 1, and real             |
//| argument u.                                                      |
//| These functions are periodic, with quarter-period on the         |
//| real axis equal to the complete elliptic integral                |
//| ellpk(1.0-m).                                                    |
//| Relation to incomplete elliptic integral:                        |
//| If u = ellik(phi,m), then sn(u|m) = sin(phi),                    |
//| and cn(u|m) = cos(phi).  Phi is called the amplitude of u.       |
//| Computation is by means of the arithmetic-geometric mean         |
//| algorithm, except when m is within 1e-9 of 0 or 1.  In the       |
//| latter case with m close to 1, the approximation applies         |
//| only for phi < pi/2.                                             |
//| ACCURACY:                                                        |
//| Tested at random points with u between 0 and 10, m between       |
//| 0 and 1.                                                         |
//|            Absolute error (* = relative error):                  |
//| arithmetic   function   # trials      peak         rms           |
//|    IEEE      phi         10000       9.2e-16*    1.4e-16*        |
//|    IEEE      sn          50000       4.1e-15     4.6e-16         |
//|    IEEE      cn          40000       3.6e-15     4.4e-16         |
//|    IEEE      dn          10000       1.3e-12     1.8e-14         |
//|  Peak error observed in consistency check using addition         |
//| theorem for sn(u+v) was 4e-16 (absolute).  Also tested by        |
//| the above relation to the incomplete elliptic integral.          |
//| Accuracy deteriorates when u is large.                           |
//+------------------------------------------------------------------+
static void CJacobianElliptic::JacobianEllipticFunctions(const double u,
                                                         const double m,
                                                         double &sn,double &cn,
                                                         double &dn,double &ph)
  {
//--- create variables
   double ai=0;
   double b=0;
   double phi=0;
   double t=0;
   double twon=0;
   int    i=0;
//--- create arrays
   double a[];
   double c[];
//--- initialization
   sn=0;
   cn=0;
   dn=0;
   ph=0;
//--- check
   if(!CAp::Assert(m>=0.0 && m<=1.0,__FUNCTION__+": m<0 or m>1"))
      return;
//--- allocation
   ArrayResizeAL(a,9);
   ArrayResizeAL(c,9);
//--- check
   if(m<1.0e-9)
     {
      //--- calculation
      t=MathSin(u);
      b=MathCos(u);
      ai=0.25*m*(u-t*b);
      sn=t-ai*b;
      cn=b+ai*t;
      ph=u-ai;
      dn=1.0-0.5*m*t*t;
      //--- exit the function
      return;
     }
//--- check
   if(m>=0.9999999999)
     {
      //--- calculation
      ai=0.25*(1.0-m);
      b=MathCosh(u);
      t=MathTanh(u);
      phi=1.0/b;
      twon=b*MathSinh(u);
      sn=t+ai*(twon-u)/(b*b);
      ph=2.0*MathArctan(MathExp(u))-1.57079632679489661923+ai*(twon-u)/b;
      ai=ai*t*phi;
      cn=phi-ai*(twon-u);
      dn=phi+ai*(twon+u);
      //--- exit the function
      return;
     }
//--- change values
   a[0]=1.0;
   b=MathSqrt(1.0-m);
   c[0]=MathSqrt(m);
   twon=1.0;
   i=0;
//--- cycle
   while(MathAbs(c[i]/a[i])>CMath::m_machineepsilon)
     {
      //--- check
      if(i>7)
        {
         //--- check
         if(!CAp::Assert(false,__FUNCTION__+": overflow"))
            return;
         //--- break the cycle
         break;
        }
      //--- calculation
      ai=a[i];
      i=i+1;
      c[i]=0.5*(ai-b);
      t=MathSqrt(ai*b);
      a[i]=0.5*(ai+b);
      b=t;
      twon=twon*2.0;
     }
   phi=twon*a[i]*u;
//--- cycle
   do
     {
      //--- calculation
      t=c[i]*MathSin(phi)/a[i];
      b=phi;
      phi=(MathArcsin(t)+phi)/2.0;
      i=i-1;
     }
   while(i!=0);
//--- get result
   sn=MathSin(phi);
   t=MathCos(phi);
   cn=t;
   dn=t/MathCos(phi-b);
   ph=phi;
  }
//+------------------------------------------------------------------+
//| Laguerre polynomial                                              |
//+------------------------------------------------------------------+
class CLaguerre
  {
public:
   //--- constructor, destructor
                     CLaguerre(void);
                    ~CLaguerre(void);
   //--- methods
   static double     LaguerreCalculate(const int n,const double x);
   static double     LaguerreSum(double &c[],const int n,const double x);
   static void       LaguerreCoefficients(const int n,double &c[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLaguerre::CLaguerre(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLaguerre::~CLaguerre(void)
  {

  }
//+------------------------------------------------------------------+
//| Calculation of the value of the Laguerre polynomial.             |
//| Parameters:                                                      |
//|     n   -   degree, n>=0                                         |
//|     x   -   argument                                             |
//| Result:                                                          |
//|     the value of the Laguerre polynomial Ln at x                 |
//+------------------------------------------------------------------+
static double CLaguerre::LaguerreCalculate(const int n,const double x)
  {
//--- create variables
   double result=0;
   double a=0;
   double b=0;
   double i=0;
//--- initialization
   result=1;
   a=1;
   b=1-x;
//--- check
   if(n==1)
      result=b;
//--- initialization
   i=2;
//--- cycle
   while(i<=n)
     {
      //--- calculation
      result=((2*i-1-x)*b-(i-1)*a)/i;
      a=b;
      b=result;
      i=i+1;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Summation of Laguerre polynomials using Clenshaw?s recurrence    |
//| formula.                                                         |
//| This routine calculates c[0]*L0(x) + c[1]*L1(x) + ... +          |
//| +  c[N]*LN(x)                                                    |
//| Parameters:                                                      |
//|     n   -   degree, n>=0                                         |
//|     x   -   argument                                             |
//| Result:                                                          |
//|     the value of the Laguerre polynomial at x                    |
//+------------------------------------------------------------------+
static double CLaguerre::LaguerreSum(double &c[],const int n,const double x)
  {
//--- create variables
   double result=0;
   double b1=0;
   double b2=0;
   int    i=0;
//--- initialization
   b1=0;
   b2=0;
   result=0;
//--- calculation
   for(i=n;i>=0;i--)
     {
      result=(2*i+1-x)*b1/(i+1)-(i+1)*b2/(i+2)+c[i];
      b2=b1;
      b1=result;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Representation of Ln as C[0] + C[1]*X + ... + C[N]*X^N           |
//| Input parameters:                                                |
//|     N   -   polynomial degree, n>=0                              |
//| Output parameters:                                               |
//|     C   -   coefficients                                         |
//+------------------------------------------------------------------+
static void CLaguerre::LaguerreCoefficients(const int n,double &c[])
  {
//--- create a variable
   int i=0;
//--- allocation
   ArrayResizeAL(c,n+1);
//--- initialization
   c[0]=1;
//--- calculation
   for(i=0;i<=n-1;i++)
      c[i+1]=-(c[i]*(n-i)/(i+1)/(i+1));
  }
//+------------------------------------------------------------------+
//| Legendre polynomial                                              |
//+------------------------------------------------------------------+
class CLegendre
  {
public:
   //--- constructor, destructor
                     CLegendre(void);
                    ~CLegendre(void);
   //--- methods
   static double     LegendreCalculate(const int n,const double x);
   static double     LegendreSum(double &c[],const int n,const double x);
   static void       LegendreCoefficients(const int n,double &c[]);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CLegendre::CLegendre(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CLegendre::~CLegendre(void)
  {

  }
//+------------------------------------------------------------------+
//| Calculation of the value of the Legendre polynomial Pn.          |
//| Parameters:                                                      |
//|     n   -   degree, n>=0                                         |
//|     x   -   argument                                             |
//| Result:                                                          |
//|     the value of the Legendre polynomial Pn at x                 |
//+------------------------------------------------------------------+
static double CLegendre::LegendreCalculate(const int n,const double x)
  {
//--- create variables
   double result=0;
   double a=0;
   double b=0;
   int    i=0;
//--- initialization
   result=1;
   a=1;
   b=x;
//--- check
   if(n==0)
      return(a);
//--- check
   if(n==1)
      return(b);
//--- calculation
   for(i=2;i<=n;i++)
     {
      result=((2*i-1)*x*b-(i-1)*a)/i;
      a=b;
      b=result;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Summation of Legendre polynomials using Clenshaw?s recurrence    |
//| formula.                                                         |
//| This routine calculates                                          |
//|     c[0]*P0(x) + c[1]*P1(x) + ... + c[N]*PN(x)                   |
//| Parameters:                                                      |
//|     n   -   degree, n>=0                                         |
//|     x   -   argument                                             |
//| Result:                                                          |
//|     the value of the Legendre polynomial at x                    |
//+------------------------------------------------------------------+
static double CLegendre::LegendreSum(double &c[],const int n,const double x)
  {
//--- create variables
   double result=0;
   double b1=0;
   double b2=0;
   int    i=0;
//--- initialization
   b1=0;
   b2=0;
   result=0;
//--- calculation
   for(i=n;i>=0;i--)
     {
      result=(2*i+1)*x*b1/(i+1)-(i+1)*b2/(i+2)+c[i];
      b2=b1;
      b1=result;
     }
//--- return result
   return(result);
  }
//+------------------------------------------------------------------+
//| Representation of Pn as C[0] + C[1]*X + ... + C[N]*X^N           |
//| Input parameters:                                                |
//|     N   -   polynomial degree, n>=0                              |
//| Output parameters:                                               |
//|     C   -   coefficients                                         |
//+------------------------------------------------------------------+
static void CLegendre::LegendreCoefficients(const int n,double &c[])
  {
//--- create a variable
   int i=0;
//--- allocation
   ArrayResizeAL(c,n+1);
//--- initialization
   for(i=0;i<=n;i++)
      c[i]=0;
   c[n]=1;
//--- calculation
   for(i=1;i<=n;i++)
      c[n]=c[n]*(n+i)/2/i;
   for(i=0;i<=n/2-1;i++)
      c[n-2*(i+1)]=-(c[n-2*i]*(n-2*i)*(n-2*i-1)/2/(i+1)/(2*(n-i)-1));
  }
//+------------------------------------------------------------------+
//| Poisson distribution                                             |
//+------------------------------------------------------------------+
class CPoissonDistr
  {
public:
   //--- constructor, destructor
                     CPoissonDistr(void);
                    ~CPoissonDistr(void);
   //--- methods
   static double     PoissonDistribution(const int k,const double m);
   static double     PoissonComplDistribution(const int k,const double m);
   static double     InvPoissonDistribution(const int k,const double y);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CPoissonDistr::CPoissonDistr(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CPoissonDistr::~CPoissonDistr(void)
  {

  }
//+------------------------------------------------------------------+
//| Poisson distribution                                             |
//| Returns the sum of the first k+1 terms of the Poisson            |
//| distribution:                                                    |
//|   k         j                                                    |
//|   --   -m  m                                                     |
//|   >   e    --                                                    |
//|   --       j!                                                    |
//|  j=0                                                             |
//| The terms are not summed directly; instead the incomplete        |
//| gamma integral is employed, according to the relation            |
//| y = pdtr( k, m ) = igamc( k+1, m ).                              |
//| The arguments must both be positive.                             |
//| ACCURACY:                                                        |
//| See incomplete gamma function                                    |
//+------------------------------------------------------------------+
static double CPoissonDistr::PoissonDistribution(const int k,const double m)
  {
//--- check
   if(!CAp::Assert(k>=0 && m>0.0,__FUNCTION__+": domain error"))
      return(EMPTY_VALUE);
//--- return result
   return(CIncGammaF::IncompleteGammaC(k+1,m));
  }
//+------------------------------------------------------------------+
//| Complemented Poisson distribution                                |
//| Returns the sum of the terms k+1 to infinity of the Poisson      |
//| distribution:                                                    |
//|  inf.       j                                                    |
//|   --   -m  m                                                     |
//|   >   e    --                                                    |
//|   --       j!                                                    |
//|  j=k+1                                                           |
//| The terms are not summed directly; instead the incomplete        |
//| gamma integral is employed, according to the formula             |
//| y = pdtrc( k, m ) = igam( k+1, m ).                              |
//| The arguments must both be positive.                             |
//| ACCURACY:                                                        |
//| See incomplete gamma function                                    |
//+------------------------------------------------------------------+
static double CPoissonDistr::PoissonComplDistribution(const int k,const double m)
  {
//--- check
   if(!CAp::Assert(k>=0 && m>0.0,__FUNCTION__+": domain error"))
      return(EMPTY_VALUE);
//--- return result
   return(CIncGammaF::IncompleteGamma(k+1,m));
  }
//+------------------------------------------------------------------+
//| Inverse Poisson distribution                                     |
//| Finds the Poisson variable x such that the integral              |
//| from 0 to x of the Poisson density is equal to the               |
//| given probability y.                                             |
//| This is accomplished using the inverse gamma integral            |
//| function and the relation                                        |
//|    m = igami( k+1, y ).                                          |
//| ACCURACY:                                                        |
//| See inverse incomplete gamma function                            |
//+------------------------------------------------------------------+
static double CPoissonDistr::InvPoissonDistribution(const int k,const double y)
  {
//--- check
   if(!CAp::Assert((k>=0 && y>=0.0) && y<1.0,__FUNCTION__+": domain error"))
      return(EMPTY_VALUE);
//--- return result
   return(CIncGammaF::InvIncompleteGammaC(k+1,y));
  }
//+------------------------------------------------------------------+
//| Psi function                                                     |
//+------------------------------------------------------------------+
class CPsiF
  {
public:
   //--- constructor, destructor
                     CPsiF(void);
                    ~CPsiF(void);
   //--- method
   static double     Psi(double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CPsiF::CPsiF(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CPsiF::~CPsiF(void)
  {

  }
//+------------------------------------------------------------------+
//| Psi (digamma) function                                           |
//|              d      -                                            |
//|   psi(x)  =  -- ln | (x)                                         |
//|              dx                                                  |
//| is the logarithmic derivative of the gamma function.             |
//| For integer x,                                                   |
//|                   n-1                                            |
//|                    -                                             |
//| psi(n) = -EUL  +   >  1/k.                                       |
//|                    -                                             |
//|                   k=1                                            |
//| This formula is used for 0 < n <= 10.  If x is negative, it      |
//| is transformed to a positive argument by the reflection          |
//| formula  psi(1-x) = psi(x) + pi cot(pi x).                       |
//| For general positive x, the argument is made greater than 10     |
//| using the recurrence  psi(x+1) = psi(x) + 1/x.                   |
//| Then the following asymptotic expansion is applied:              |
//|                           inf.   B                               |
//|                            -      2k                             |
//| psi(x) = log(x) - 1/2x -   >   -------                           |
//|                            -        2k                           |
//|                           k=1   2k x                             |
//| where the B2k are Bernoulli numbers.                             |
//| ACCURACY:                                                        |
//|    Relative error (except absolute when |psi| < 1):              |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE      0,30        30000       1.3e-15     1.4e-16         |
//|    IEEE      -30,0       40000       1.5e-15     2.2e-16         |
//+------------------------------------------------------------------+
static double CPsiF::Psi(double x)
  {
//--- create variables
   double result=0;
   double p=0;
   double q=0;
   double nz=0;
   double s=0;
   double w=0;
   double y=0;
   double z=0;
   double polv=0;
   int    i=0;
   int    n=0;
   int    negative=0;
//--- initialization
   negative=0;
   nz=0.0;
//--- check
   if(x<=0.0)
     {
      negative=1;
      q=x;
      p=(int)MathFloor(q);
      //--- check
      if(p==q)
        {
         //--- check
         if(!CAp::Assert(false,__FUNCTION__+": singularity in Psi(x)"))
            return(EMPTY_VALUE);
         //--- return result
         return(CMath::m_maxrealnumber);
        }
      nz=q-p;
      //--- check
      if(nz!=0.5)
        {
         //--- check
         if(nz>0.5)
           {
            p=p+1.0;
            nz=q-p;
           }
         //--- change value
         nz=M_PI/MathTan(M_PI*nz);
        }
      else
         nz=0.0;
      x=1.0-x;
     }
//--- check
   if(x<=10.0 && x==(double)((int)MathFloor(x)))
     {
      y=0.0;
      n=(int)MathFloor(x);
      //--- calculation
      for(i=1;i<=n-1;i++)
        {
         w=i;
         y=y+1.0/w;
        }
      y=y-0.57721566490153286061;
     }
   else
     {
      //--- change values
      s=x;
      w=0.0;
      //--- cycle
      while(s<10.0)
        {
         w=w+1.0/s;
         s=s+1.0;
        }
      //--- check
      if(s<1.0E17)
        {
         //--- calculation
         z=1.0/(s*s);
         polv=8.33333333333333333333E-2;
         polv=polv*z-2.10927960927960927961E-2;
         polv=polv*z+7.57575757575757575758E-3;
         polv=polv*z-4.16666666666666666667E-3;
         polv=polv*z+3.96825396825396825397E-3;
         polv=polv*z-8.33333333333333333333E-3;
         polv=polv*z+8.33333333333333333333E-2;
         y=z*polv;
        }
      else
         y=0.0;
      //--- change value
      y=MathLog(s)-0.5/s-y-w;
     }
//--- check
   if(negative!=0)
      y=y-nz;
//--- return result
   return(y);
  }
//+------------------------------------------------------------------+
//| Student's t distribution                                         |
//+------------------------------------------------------------------+
class CStudenttDistr
  {
public:
   //--- constructor, destructor
                     CStudenttDistr(void);
                    ~CStudenttDistr(void);
   //--- methods
   static double     StudenttDistribution(const int k,const double t);
   static double     InvStudenttDistribution(const int k,double p);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CStudenttDistr::CStudenttDistr(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CStudenttDistr::~CStudenttDistr(void)
  {

  }
//+------------------------------------------------------------------+
//| Student's t distribution                                         |
//| Computes the integral from minus infinity to t of the Student    |
//| t distribution with integer k > 0 degrees of freedom:            |
//|                                      t                           |
//|                                      -                           |
//|                                     | |                          |
//|              -                      |         2   -(k+1)/2       |
//|             | ( (k+1)/2 )           |  (     x   )               |
//|       ----------------------        |  ( 1 + --- )        dx     |
//|                     -               |  (      k  )               |
//|       sqrt( k pi ) | ( k/2 )        |                            |
//|                                   | |                            |
//|                                    -                             |
//|                                   -inf.                          |
//| Relation to incomplete beta integral:                            |
//|        1 - stdtr(k,t) = 0.5 * incbet( k/2, 1/2, z )              |
//| where                                                            |
//|        z = k/(k + t**2).                                         |
//| For t < -2, this is the method of computation.  For higher t,    |
//| a direct method is derived from integration by parts.            |
//| Since the function is symmetric about t=0, the area under the    |
//| right tail of the density is found by calling the function       |
//| with -t instead of t.                                            |
//| ACCURACY:                                                        |
//| Tested at random 1 <= k <= 25.  The "domain" refers to t.        |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE     -100,-2      50000       5.9e-15     1.4e-15         |
//|    IEEE     -2,100      500000       2.7e-15     4.9e-17         |
//+------------------------------------------------------------------+
static double CStudenttDistr::StudenttDistribution(const int k,const double t)
  {
//--- check
   if(!CAp::Assert(k>0,__FUNCTION__+": the error variable"))
      return(EMPTY_VALUE);
//--- check
   if(t==0)
      return(0.5);
//--- create variables
   double x=0;
   double rk=0;
   double z=0;
   double f=0;
   double tz=0;
   double p=0;
   double xsqk=0;
   int    j=0;
//--- check
   if(t<-2.0)
     {
      rk=k;
      z=rk/(rk+t*t);
      //--- return result
      return(0.5*CIncBetaF::IncompleteBeta(0.5*rk,0.5,z));
     }
//--- check
   if(t<0)
      x=-t;
   else
      x=t;
//--- change values
   rk=k;
   z=1.0+x*x/rk;
//--- check
   if(k%2!=0)
     {
      xsqk=x/MathSqrt(rk);
      p=MathArctan(xsqk);
      //--- check
      if(k>1)
        {
         f=1.0;
         tz=1.0;
         j=3;
         //--- cycle
         while(j<k-1 && tz/f>CMath::m_machineepsilon)
           {
            //--- calculation
            tz*=(j-1)/(z*j);
            f+=tz;
            j+=2;
           }
         p+=f*xsqk/z;
        }
      p*=2.0/M_PI;
     }
   else
     {
      f=1.0;
      tz=1.0;
      j=2;
      //--- cycle
      while(j<k-1 && tz/f>CMath::m_machineepsilon)
        {
         //--- calculation
         tz*=(j-1)/(z*j);
         f+=tz;
         j+=2;
        }
      p=f*x/MathSqrt(z*rk);
     }
//--- check
   if(t<0)
      p=-p;
//--- return result
   return(0.5+0.5*p);
  }
//+------------------------------------------------------------------+
//| Functional inverse of Student's t distribution                   |
//| Given probability p, finds the argument t such that stdtr(k,t)   |
//| is equal to p.                                                   |
//| ACCURACY:                                                        |
//| Tested at random 1 <= k <= 100.  The "domain" refers to p:       |
//|                      Relative error:                             |
//| arithmetic   domain     # trials      peak         rms           |
//|    IEEE    .001,.999     25000       5.7e-15     8.0e-16         |
//|    IEEE    10^-6,.001    25000       2.0e-12     2.9e-14         |
//+------------------------------------------------------------------+
static double CStudenttDistr::InvStudenttDistribution(const int k,double p)
  {
//--- create variables
   double    t=0;
   double    rk=k;
   double    z=0;
   int       rflg=0;
   CIncBetaF ibetaf;
//--- check
   if(!CAp::Assert(k>0 && p>0 && p<1,__FUNCTION__+": the error variable"))
      return(EMPTY_VALUE);
//--- check
   if(p>0.25 && p<0.75)
     {
      //--- check
      if(p==0.5)
         return(0);
      //--- change values
      z=1.0-2.0*p;
      z=CIncBetaF::InvIncompleteBeta(0.5,0.5*rk,MathAbs(z));
      t=MathSqrt(rk*z/(1.0-z));
      //--- check
      if(p<0.5)
         t=-t;
      //--- return result
      return(t);
     }
//--- 0<p<=0.25 or 0.75<=p<1
   rflg=-1;
//--- check
   if(p>=0.5)
     {
      p=1.0-p;
      rflg=1;
     }
//--- change value
   z=CIncBetaF::InvIncompleteBeta(0.5*rk,0.5,2.0*p);
//--- check
   if(CMath::m_maxrealnumber*z<rk)
      return(rflg*CMath::m_maxrealnumber);
   t=MathSqrt(rk/z-rk);
//--- return result
   return(rflg*t);
  }
//+------------------------------------------------------------------+
//| Trigonometric integrals                                          |
//+------------------------------------------------------------------+
class CTrigIntegrals
  {
private:
   //--- private method
   static void       ChebIterationShiChi(const double x,const double c,double &b0,double &b1,double &b2);
public:
   //--- constructor, destructor
                     CTrigIntegrals(void);
                    ~CTrigIntegrals(void);
   //--- public methods
   static void       SineCosineIntegrals(double x,double &si,double &ci);
   static void       HyperbolicSineCosineIntegrals(double x,double &shi,double &chi);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CTrigIntegrals::CTrigIntegrals(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CTrigIntegrals::~CTrigIntegrals(void)
  {

  }
//+------------------------------------------------------------------+
//| Sine and cosine integrals                                        |
//| Evaluates the integrals                                          |
//|                          x                                       |
//|                          -                                       |
//|                         |  cos t - 1                             |
//|   Ci(x) = eul + ln x +  |  --------- dt,                         |
//|                         |      t                                 |
//|                        -                                         |
//|                         0                                        |
//|             x                                                    |
//|             -                                                    |
//|            |  sin t                                              |
//|   Si(x) =  |  ----- dt                                           |
//|            |    t                                                |
//|           -                                                      |
//|            0                                                     |
//| where eul = 0.57721566490153286061 is Euler's constant.          |
//| The integrals are approximated by rational functions.            |
//| For x > 8 auxiliary functions f(x) and g(x) are employed         |
//| such that                                                        |
//| Ci(x) = f(x) sin(x) - g(x) cos(x)                                |
//| Si(x) = pi/2 - f(x) cos(x) - g(x) sin(x)                         |
//| ACCURACY:                                                        |
//|    Test interval = [0,50].                                       |
//| Absolute error, except relative when > 1:                        |
//| arithmetic   function   # trials      peak         rms           |
//|    IEEE        Si        30000       4.4e-16     7.3e-17         |
//|    IEEE        Ci        30000       6.9e-16     5.1e-17         |
//+------------------------------------------------------------------+
static void CTrigIntegrals::SineCosineIntegrals(double x,double &si,double &ci)
  {
//--- create variables
   double z=0;
   double c=0;
   double s=0;
   double f=0;
   double g=0;
   int    sg=0;
   double sn=0;
   double sd=0;
   double cn=0;
   double cd=0;
   double fn=0;
   double fd=0;
   double gn=0;
   double gd=0;
//--- initialization
   si=0;
   ci=0;
//--- check
   if(x<0.0)
     {
      sg=-1;
      x=-x;
     }
   else
      sg=0;
//--- check
   if(x==0.0)
     {
      si=0;
      ci=-CMath::m_maxrealnumber;
      //--- exit the function
      return;
     }
//--- check
   if(x>1.0E9)
     {
      si=1.570796326794896619-MathCos(x)/x;
      ci=MathSin(x)/x;
      //--- exit the function
      return;
     }
//--- check
   if(x<=4.0)
     {
      //--- calculation
      z=x*x;
      sn=-8.39167827910303881427E-11;
      sn=sn*z+4.62591714427012837309E-8;
      sn=sn*z-9.75759303843632795789E-6;
      sn=sn*z+9.76945438170435310816E-4;
      sn=sn*z-4.13470316229406538752E-2;
      sn=sn*z+1.00000000000000000302E0;
      sd=2.03269266195951942049E-12;
      sd=sd*z+1.27997891179943299903E-9;
      sd=sd*z+4.41827842801218905784E-7;
      sd=sd*z+9.96412122043875552487E-5;
      sd=sd*z+1.42085239326149893930E-2;
      sd=sd*z+9.99999999999999996984E-1;
      s=x*sn/sd;
      cn=2.02524002389102268789E-11;
      cn=cn*z-1.35249504915790756375E-8;
      cn=cn*z+3.59325051419993077021E-6;
      cn=cn*z-4.74007206873407909465E-4;
      cn=cn*z+2.89159652607555242092E-2;
      cn=cn*z-1.00000000000000000080E0;
      cd=4.07746040061880559506E-12;
      cd=cd*z+3.06780997581887812692E-9;
      cd=cd*z+1.23210355685883423679E-6;
      cd=cd*z+3.17442024775032769882E-4;
      cd=cd*z+5.10028056236446052392E-2;
      cd=cd*z+4.00000000000000000080E0;
      c=z*cn/cd;
      //--- check
      if(sg!=0)
         s=-s;
      //--- calculation
      si=s;
      ci=0.57721566490153286061+MathLog(x)+c;
      //--- exit the function
      return;
     }
//--- change values
   s=MathSin(x);
   c=MathCos(x);
   z=1.0/(x*x);
//--- check
   if(x<8.0)
     {
      //--- calculation
      fn=4.23612862892216586994E0;
      fn=fn*z+5.45937717161812843388E0;
      fn=fn*z+1.62083287701538329132E0;
      fn=fn*z+1.67006611831323023771E-1;
      fn=fn*z+6.81020132472518137426E-3;
      fn=fn*z+1.08936580650328664411E-4;
      fn=fn*z+5.48900223421373614008E-7;
      fd=1.00000000000000000000E0;
      fd=fd*z+8.16496634205391016773E0;
      fd=fd*z+7.30828822505564552187E0;
      fd=fd*z+1.86792257950184183883E0;
      fd=fd*z+1.78792052963149907262E-1;
      fd=fd*z+7.01710668322789753610E-3;
      fd=fd*z+1.10034357153915731354E-4;
      fd=fd*z+5.48900252756255700982E-7;
      f=fn/(x*fd);
      gn=8.71001698973114191777E-2;
      gn=gn*z+6.11379109952219284151E-1;
      gn=gn*z+3.97180296392337498885E-1;
      gn=gn*z+7.48527737628469092119E-2;
      gn=gn*z+5.38868681462177273157E-3;
      gn=gn*z+1.61999794598934024525E-4;
      gn=gn*z+1.97963874140963632189E-6;
      gn=gn*z+7.82579040744090311069E-9;
      gd=1.00000000000000000000E0;
      gd=gd*z+1.64402202413355338886E0;
      gd=gd*z+6.66296701268987968381E-1;
      gd=gd*z+9.88771761277688796203E-2;
      gd=gd*z+6.22396345441768420760E-3;
      gd=gd*z+1.73221081474177119497E-4;
      gd=gd*z+2.02659182086343991969E-6;
      gd=gd*z+7.82579218933534490868E-9;
      g=z*gn/gd;
     }
   else
     {
      //--- calculation
      fn=4.55880873470465315206E-1;
      fn=fn*z+7.13715274100146711374E-1;
      fn=fn*z+1.60300158222319456320E-1;
      fn=fn*z+1.16064229408124407915E-2;
      fn=fn*z+3.49556442447859055605E-4;
      fn=fn*z+4.86215430826454749482E-6;
      fn=fn*z+3.20092790091004902806E-8;
      fn=fn*z+9.41779576128512936592E-11;
      fn=fn*z+9.70507110881952024631E-14;
      fd=1.00000000000000000000E0;
      fd=fd*z+9.17463611873684053703E-1;
      fd=fd*z+1.78685545332074536321E-1;
      fd=fd*z+1.22253594771971293032E-2;
      fd=fd*z+3.58696481881851580297E-4;
      fd=fd*z+4.92435064317881464393E-6;
      fd=fd*z+3.21956939101046018377E-8;
      fd=fd*z+9.43720590350276732376E-11;
      fd=fd*z+9.70507110881952025725E-14;
      f=fn/(x*fd);
      gn=6.97359953443276214934E-1;
      gn=gn*z+3.30410979305632063225E-1;
      gn=gn*z+3.84878767649974295920E-2;
      gn=gn*z+1.71718239052347903558E-3;
      gn=gn*z+3.48941165502279436777E-5;
      gn=gn*z+3.47131167084116673800E-7;
      gn=gn*z+1.70404452782044526189E-9;
      gn=gn*z+3.85945925430276600453E-12;
      gn=gn*z+3.14040098946363334640E-15;
      gd=1.00000000000000000000E0;
      gd=gd*z+1.68548898811011640017E0;
      gd=gd*z+4.87852258695304967486E-1;
      gd=gd*z+4.67913194259625806320E-2;
      gd=gd*z+1.90284426674399523638E-3;
      gd=gd*z+3.68475504442561108162E-5;
      gd=gd*z+3.57043223443740838771E-7;
      gd=gd*z+1.72693748966316146736E-9;
      gd=gd*z+3.87830166023954706752E-12;
      gd=gd*z+3.14040098946363335242E-15;
      g=z*gn/gd;
     }
   si=1.570796326794896619-f*c-g*s;
//--- check
   if(sg!=0)
      si=-si;
//--- get result
   ci=f*s-g*c;
  }
//+------------------------------------------------------------------+
//| Hyperbolic sine and cosine integrals                             |
//| Approximates the integrals                                       |
//|                            x                                     |
//|                            -                                     |
//|                           | |   cosh t - 1                       |
//|   Chi(x) = eul + ln x +   |    -----------  dt,                  |
//|                         | |          t                           |
//|                          -                                       |
//|                          0                                       |
//|               x                                                  |
//|               -                                                  |
//|              | |  sinh t                                         |
//|   Shi(x) =   |    ------  dt                                     |
//|            | |       t                                           |
//|             -                                                    |
//|             0                                                    |
//| where eul = 0.57721566490153286061 is Euler's constant.          |
//| The integrals are evaluated by power series for x < 8            |
//| and by Chebyshev expansions for x between 8 and 88.              |
//| For large x, both functions approach exp(x)/2x.                  |
//| Arguments greater than 88 in magnitude return MAXNUM.            |
//| ACCURACY:                                                        |
//| Test interval 0 to 88.                                           |
//|                      Relative error:                             |
//| arithmetic   function  # trials      peak         rms            |
//|    IEEE         Shi      30000       6.9e-16     1.6e-16         |
//|        Absolute error, except relative when |Chi| > 1:           |
//|    IEEE         Chi      30000       8.4e-16     1.4e-16         |
//+------------------------------------------------------------------+
static void CTrigIntegrals::HyperbolicSineCosineIntegrals(double x,
                                                          double &shi,
                                                          double &chi)
  {
//--- create variables
   double k=0;
   double z=0;
   double c=0;
   double s=0;
   double a=0;
   int    sg=0;
   double b0=0;
   double b1=0;
   double b2=0;
//--- initialization
   shi=0;
   chi=0;
//--- check
   if(x<0.0)
     {
      sg=-1;
      x=-x;
     }
   else
      sg=0;
//--- check
   if(x==0.0)
     {
      shi=0;
      chi=-CMath::m_maxrealnumber;
      //--- exit the function
      return;
     }
//--- check
   if(x<8.0)
     {
      //--- initialization
      z=x*x;
      a=1.0;
      s=1.0;
      c=0.0;
      k=2.0;
      do
        {
         //--- calculation
         a=a*z/k;
         c=c+a/k;
         k=k+1.0;
         a=a/k;
         s=s+a/k;
         k=k+1.0;
        }
      while(MathAbs(a/s)>=CMath::m_machineepsilon);
      s=s*x;
     }
   else
     {
      //--- check
      if(x<18.0)
        {
         //--- calculation
         a=(576.0/x-52.0)/10.0;
         k=MathExp(x)/x;
         b0=1.83889230173399459482E-17;
         b1=0.0;
         //--- function calls
         ChebIterationShiChi(a,-9.55485532279655569575E-17,b0,b1,b2);
         ChebIterationShiChi(a,2.04326105980879882648E-16,b0,b1,b2);
         ChebIterationShiChi(a,1.09896949074905343022E-15,b0,b1,b2);
         ChebIterationShiChi(a,-1.31313534344092599234E-14,b0,b1,b2);
         ChebIterationShiChi(a,5.93976226264314278932E-14,b0,b1,b2);
         ChebIterationShiChi(a,-3.47197010497749154755E-14,b0,b1,b2);
         ChebIterationShiChi(a,-1.40059764613117131000E-12,b0,b1,b2);
         ChebIterationShiChi(a,9.49044626224223543299E-12,b0,b1,b2);
         ChebIterationShiChi(a,-1.61596181145435454033E-11,b0,b1,b2);
         ChebIterationShiChi(a,-1.77899784436430310321E-10,b0,b1,b2);
         ChebIterationShiChi(a,1.35455469767246947469E-9,b0,b1,b2);
         ChebIterationShiChi(a,-1.03257121792819495123E-9,b0,b1,b2);
         ChebIterationShiChi(a,-3.56699611114982536845E-8,b0,b1,b2);
         ChebIterationShiChi(a,1.44818877384267342057E-7,b0,b1,b2);
         ChebIterationShiChi(a,7.82018215184051295296E-7,b0,b1,b2);
         ChebIterationShiChi(a,-5.39919118403805073710E-6,b0,b1,b2);
         ChebIterationShiChi(a,-3.12458202168959833422E-5,b0,b1,b2);
         ChebIterationShiChi(a,8.90136741950727517826E-5,b0,b1,b2);
         ChebIterationShiChi(a,2.02558474743846862168E-3,b0,b1,b2);
         ChebIterationShiChi(a,2.96064440855633256972E-2,b0,b1,b2);
         ChebIterationShiChi(a,1.11847751047257036625E0,b0,b1,b2);
         //--- calculation
         s=k*0.5*(b0-b2);
         b0=-8.12435385225864036372E-18;
         b1=0.0;
         //--- function calls
         ChebIterationShiChi(a,2.17586413290339214377E-17,b0,b1,b2);
         ChebIterationShiChi(a,5.22624394924072204667E-17,b0,b1,b2);
         ChebIterationShiChi(a,-9.48812110591690559363E-16,b0,b1,b2);
         ChebIterationShiChi(a,5.35546311647465209166E-15,b0,b1,b2);
         ChebIterationShiChi(a,-1.21009970113732918701E-14,b0,b1,b2);
         ChebIterationShiChi(a,-6.00865178553447437951E-14,b0,b1,b2);
         ChebIterationShiChi(a,7.16339649156028587775E-13,b0,b1,b2);
         ChebIterationShiChi(a,-2.93496072607599856104E-12,b0,b1,b2);
         ChebIterationShiChi(a,-1.40359438136491256904E-12,b0,b1,b2);
         ChebIterationShiChi(a,8.76302288609054966081E-11,b0,b1,b2);
         ChebIterationShiChi(a,-4.40092476213282340617E-10,b0,b1,b2);
         ChebIterationShiChi(a,-1.87992075640569295479E-10,b0,b1,b2);
         ChebIterationShiChi(a,1.31458150989474594064E-8,b0,b1,b2);
         ChebIterationShiChi(a,-4.75513930924765465590E-8,b0,b1,b2);
         ChebIterationShiChi(a,-2.21775018801848880741E-7,b0,b1,b2);
         ChebIterationShiChi(a,1.94635531373272490962E-6,b0,b1,b2);
         ChebIterationShiChi(a,4.33505889257316408893E-6,b0,b1,b2);
         ChebIterationShiChi(a,-6.13387001076494349496E-5,b0,b1,b2);
         ChebIterationShiChi(a,-3.13085477492997465138E-4,b0,b1,b2);
         ChebIterationShiChi(a,4.97164789823116062801E-4,b0,b1,b2);
         ChebIterationShiChi(a,2.64347496031374526641E-2,b0,b1,b2);
         ChebIterationShiChi(a,1.11446150876699213025E0,b0,b1,b2);
         c=k*0.5*(b0-b2);
        }
      else
        {
         //--- check
         if(x<=88.0)
           {
            //--- calculation
            a=(6336.0/x-212.0)/70.0;
            k=MathExp(x)/x;
            b0=-1.05311574154850938805E-17;
            b1=0.0;
            //--- function calls
            ChebIterationShiChi(a,2.62446095596355225821E-17,b0,b1,b2);
            ChebIterationShiChi(a,8.82090135625368160657E-17,b0,b1,b2);
            ChebIterationShiChi(a,-3.38459811878103047136E-16,b0,b1,b2);
            ChebIterationShiChi(a,-8.30608026366935789136E-16,b0,b1,b2);
            ChebIterationShiChi(a,3.93397875437050071776E-15,b0,b1,b2);
            ChebIterationShiChi(a,1.01765565969729044505E-14,b0,b1,b2);
            ChebIterationShiChi(a,-4.21128170307640802703E-14,b0,b1,b2);
            ChebIterationShiChi(a,-1.60818204519802480035E-13,b0,b1,b2);
            ChebIterationShiChi(a,3.34714954175994481761E-13,b0,b1,b2);
            ChebIterationShiChi(a,2.72600352129153073807E-12,b0,b1,b2);
            ChebIterationShiChi(a,1.66894954752839083608E-12,b0,b1,b2);
            ChebIterationShiChi(a,-3.49278141024730899554E-11,b0,b1,b2);
            ChebIterationShiChi(a,-1.58580661666482709598E-10,b0,b1,b2);
            ChebIterationShiChi(a,-1.79289437183355633342E-10,b0,b1,b2);
            ChebIterationShiChi(a,1.76281629144264523277E-9,b0,b1,b2);
            ChebIterationShiChi(a,1.69050228879421288846E-8,b0,b1,b2);
            ChebIterationShiChi(a,1.25391771228487041649E-7,b0,b1,b2);
            ChebIterationShiChi(a,1.16229947068677338732E-6,b0,b1,b2);
            ChebIterationShiChi(a,1.61038260117376323993E-5,b0,b1,b2);
            ChebIterationShiChi(a,3.49810375601053973070E-4,b0,b1,b2);
            ChebIterationShiChi(a,1.28478065259647610779E-2,b0,b1,b2);
            ChebIterationShiChi(a,1.03665722588798326712E0,b0,b1,b2);
            //--- calculation
            s=k*0.5*(b0-b2);
            b0=8.06913408255155572081E-18;
            b1=0.0;
            //--- function calls
            ChebIterationShiChi(a,-2.08074168180148170312E-17,b0,b1,b2);
            ChebIterationShiChi(a,-5.98111329658272336816E-17,b0,b1,b2);
            ChebIterationShiChi(a,2.68533951085945765591E-16,b0,b1,b2);
            ChebIterationShiChi(a,4.52313941698904694774E-16,b0,b1,b2);
            ChebIterationShiChi(a,-3.10734917335299464535E-15,b0,b1,b2);
            ChebIterationShiChi(a,-4.42823207332531972288E-15,b0,b1,b2);
            ChebIterationShiChi(a,3.49639695410806959872E-14,b0,b1,b2);
            ChebIterationShiChi(a,6.63406731718911586609E-14,b0,b1,b2);
            ChebIterationShiChi(a,-3.71902448093119218395E-13,b0,b1,b2);
            ChebIterationShiChi(a,-1.27135418132338309016E-12,b0,b1,b2);
            ChebIterationShiChi(a,2.74851141935315395333E-12,b0,b1,b2);
            ChebIterationShiChi(a,2.33781843985453438400E-11,b0,b1,b2);
            ChebIterationShiChi(a,2.71436006377612442764E-11,b0,b1,b2);
            ChebIterationShiChi(a,-2.56600180000355990529E-10,b0,b1,b2);
            ChebIterationShiChi(a,-1.61021375163803438552E-9,b0,b1,b2);
            ChebIterationShiChi(a,-4.72543064876271773512E-9,b0,b1,b2);
            ChebIterationShiChi(a,-3.00095178028681682282E-9,b0,b1,b2);
            ChebIterationShiChi(a,7.79387474390914922337E-8,b0,b1,b2);
            ChebIterationShiChi(a,1.06942765566401507066E-6,b0,b1,b2);
            ChebIterationShiChi(a,1.59503164802313196374E-5,b0,b1,b2);
            ChebIterationShiChi(a,3.49592575153777996871E-4,b0,b1,b2);
            ChebIterationShiChi(a,1.28475387530065247392E-2,b0,b1,b2);
            ChebIterationShiChi(a,1.03665693917934275131E0,b0,b1,b2);
            c=k*0.5*(b0-b2);
           }
         else
           {
            //--- check
            if(sg!=0)
               shi=-CMath::m_maxrealnumber;
            else
               shi=CMath::m_maxrealnumber;
            chi=CMath::m_maxrealnumber;
            //--- exit the function
            return;
           }
        }
     }
//--- check
   if(sg!=0)
      s=-s;
//--- get result
   shi=s;
   chi=0.57721566490153286061+MathLog(x)+c;
  }
//+------------------------------------------------------------------+
//| Internal subroutine                                              |
//+------------------------------------------------------------------+
static void CTrigIntegrals::ChebIterationShiChi(const double x,const double c,
                                                double &b0,double &b1,double &b2)
  {
//--- change values
   b2=b1;
   b1=b0;
   b0=x*b1-b2+c;
  }
//+------------------------------------------------------------------+
