//+------------------------------------------------------------------+
//|                                                       Normal.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"

const static double normal_cdf_a[5]=
  {
   2.2352520354606839287E00,1.6102823106855587881E02,
   1.0676894854603709582E03,1.8154981253343561249E04,
   6.5682337918207449113E-2
  };
const static double normal_cdf_b[4]=
  {
   4.7202581904688241870E01,9.7609855173777669322E02,
   1.0260932208618978205E04,4.5507789335026729956E04
  };
//--- coefficients for approximation in second interval
const static double normal_cdf_c[9]=
  {
   3.9894151208813466764E-1,8.8831497943883759412E00,
   9.3506656132177855979E01,5.9727027639480026226E02,
   2.4945375852903726711E03,6.8481904505362823326E03,
   1.1602651437647350124E04,9.8427148383839780218E03,
   1.0765576773720192317E-8
  };
const static double normal_cdf_d[8]=
  {
   2.2266688044328115691E01,2.3538790178262499861E02,
   1.5193775994075548050E03,6.4855582982667607550E03,
   1.8615571640885098091E04,3.4900952721145977266E04,
   3.8912003286093271411E04,1.9685429676859990727E04
  };
//--- coefficients for approximation in third interval
const static double normal_cdf_p[6]=
  {
   2.1589853405795699E-1,1.274011611602473639E-1,
   2.2235277870649807E-2,1.421619193227893466E-3,
   2.9112874951168792E-5,2.307344176494017303E-2
  };
const static double normal_cdf_q[5]=
  {
   1.28426009614491121E00,4.68238212480865118E-1,
   6.59881378689285515E-2,3.78239633202758244E-3,
   7.29751555083966205E-5
  };

//--- coefficients for p close to 0.5
const double normal_q_a0 = 3.3871328727963666080;
const double normal_q_a1 = 1.3314166789178437745E+2;
const double normal_q_a2 = 1.9715909503065514427E+3;
const double normal_q_a3 = 1.3731693765509461125E+4;
const double normal_q_a4 = 4.5921953931549871457E+4;
const double normal_q_a5 = 6.7265770927008700853E+4;
const double normal_q_a6 = 3.3430575583588128105E+4;
const double normal_q_a7 = 2.5090809287301226727E+3;
const double normal_q_b1 = 4.2313330701600911252E+1;
const double normal_q_b2 = 6.8718700749205790830E+2;
const double normal_q_b3 = 5.3941960214247511077E+3;
const double normal_q_b4 = 2.1213794301586595867E+4;
const double normal_q_b5 = 3.9307895800092710610E+4;
const double normal_q_b6 = 2.8729085735721942674E+4;
const double normal_q_b7 = 5.2264952788528545610E+3;
//--- coefficients for p not close to 0, 0.5 or 1
const double normal_q_c0 = 1.42343711074968357734;
const double normal_q_c1 = 4.63033784615654529590;
const double normal_q_c2 = 5.76949722146069140550;
const double normal_q_c3 = 3.64784832476320460504;
const double normal_q_c4 = 1.27045825245236838258;
const double normal_q_c5 = 2.41780725177450611770E-1;
const double normal_q_c6 = 2.27238449892691845833E-2;
const double normal_q_c7 = 7.74545014278341407640E-4;
const double normal_q_d1 = 2.05319162663775882187;
const double normal_q_d2 = 1.67638483018380384940;
const double normal_q_d3 = 6.89767334985100004550E-1;
const double normal_q_d4 = 1.48103976427480074590E-1;
const double normal_q_d5 = 1.51986665636164571966E-2;
const double normal_q_d6 = 5.47593808499534494600E-4;
const double normal_q_d7 = 1.05075007164441684324E-9;
//--- coefficients for p near 0 or 1.
const double normal_q_e0 = 6.65790464350110377720E0;
const double normal_q_e1 = 5.46378491116411436990E0;
const double normal_q_e2 = 1.78482653991729133580E0;
const double normal_q_e3 = 2.96560571828504891230E-1;
const double normal_q_e4 = 2.65321895265761230930E-2;
const double normal_q_e5 = 1.24266094738807843860E-3;
const double normal_q_e6 = 2.71155556874348757815E-5;
const double normal_q_e7 = 2.01033439929228813265E-7;
const double normal_q_f1 = 5.99832206555887937690E-1;
const double normal_q_f2 = 1.36929880922735805310E-1;
const double normal_q_f3 = 1.48753612908506148525E-2;
const double normal_q_f4 = 7.86869131145613259100E-4;
const double normal_q_f5 = 1.84631831751005468180E-5;
const double normal_q_f6 = 1.42151175831644588870E-7;
const double normal_q_f7 = 2.04426310338993978564E-15;
//+------------------------------------------------------------------+
//| Normal probability density function (PDF)                        |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Normal distribution with parameters mu and sigma.         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| mu         : Mean                                                |
//| sigma      : Standard deviation (sigma>0)                        |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityNormal(const double x,const double mu,const double sigma,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check sigma
   if(sigma<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;

//--- prepare argument
   double y=(x-mu)/sigma;
//--- check it
   if(!MathIsValidNumber(y))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check overflow
   y=MathAbs(y);
   if(y>=2*MathSqrt(DBL_MAX))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- return density
   return TailLogValue(M_1_SQRT_2PI*MathExp(-0.5*y*y)/sigma,true,log_mode);
  }
//+------------------------------------------------------------------+
//| Normal probability density function (PDF)                        |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Normal distribution with parameters mu and sigma.         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| mu         : Mean                                                |
//| sigma      : Standard deviation (sigma>0)                        |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityNormal(const double x,const double mu,const double sigma,int &error_code)
  {
   return MathProbabilityDensityNormal(x,mu,sigma,false,error_code);
  }
//+------------------------------------------------------------------+
//| Normal probability density function (PDF)                        |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Normal distribution with parameters mu and sigma             |
//| for values in x.                                                 |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Mean                                                |
//| sigma      : Standard deviation (sigma>0)                        |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityNormal(const double &x[],const double mu,const double sigma,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
      return false;
//--- check sigma
   if(sigma<=0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);

   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(!MathIsValidNumber(x_arg))
         return false;

      //--- prepare argument and check it
      double y=(x_arg-mu)/sigma;
      if(!MathIsValidNumber(y))
         return false;

      //--- check overflow
      y=MathAbs(y);
      if(y>=2*MathSqrt(DBL_MAX))
         return false;

      //--- calculate density
      result[i]=TailLogValue(M_1_SQRT_2PI*MathExp(-0.5*y*y)/sigma,true,log_mode);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Normal probability density function (PDF)                        |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Normal distribution with parameters mu and sigma             |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Mean                                                |
//| sigma      : Standard deviation (sigma>0)                        |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityNormal(const double &x[],const double mu,const double sigma,double &result[])
  {
   return MathProbabilityDensityNormal(x,mu,sigma,false,result);
  }
//+------------------------------------------------------------------+
//| Normal cumulative distribution function (CDF)                    |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Normal distribution with parameters mu and sigma        |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| mu         : Mean                                                |
//| sigma      : Standard deviation (must be positive)               |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Normal cumulative distribution function with    |
//| parameters mu and sigma, evaluated at x.                         |
//+------------------------------------------------------------------+
//| Comment from original FORTRAN code                               |
//| http://www.netlib.org/toms-2014-06-10/639                        |
//| http://www.netlib.org/toms-2014-06-10/715                        |
//|                                                                  |
//| This function evaluates the normal distribution function:        |
//|                                                                  |
//|                              / x                                 |
//|                     1       |       -t*t/2                       |
//|          P(x) = ----------- |      e       dt                    |
//|                 sqrt(2 pi)  |                                    |
//|                             /-oo                                 |
//|                                                                  |
//| The main computation evaluates near-minimax approximations       |
//| derived from those in "Rational Chebyshev approximations for     |
//| the error function" by W. J. Cody, Math. Comp., 1969, 631-637.   |
//| This transportable program uses rational functions that          |
//| theoretically approximate the normal distribution function to    |
//| at least 18 significant decimal digits.  The accuracy achieved   |
//| depends on the arithmetic system, the compiler, the intrinsic    |
//| functions, and proper selection of the machine-dependent         |
//| constants.                                                       |
//|                                                                  |
//| Author:                                                          |
//| W. J. Cody, Mathematics and Computer Science Division            |
//| Argonne National Laboratory, Argonne, IL 60439                   |
//+------------------------------------------------------------------+
double MathCumulativeDistributionNormal(const double x,const double mu,const double sigma,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check sigma
   if(sigma<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- prepare argument
   double xx=(x-mu)/sigma;
//--- mathematical constants
//--- sqrpi = 1 / sqrt(2*pi), root32 = sqrt(32), and
//--- thrsh is the argument for which anorm = 0.75.
   const double sqrpi=1.0/MathSqrt(2*M_PI);
   const double thrsh = 0.66291e0;
   const double root32= MathSqrt(32);
//--- machine-dependent constants
//--- data eps/5.96e-8/,xlow/-12.949e0/,xuppr/5.768e0/
   const double eps=1.11e-16;
   const double xlow=-37.519;
   const double xuppr=8.572;
   int k;
//---
   double xsq=0.0;
   double y=MathAbs(xx);
   double xnum=0.0;
   double xden=0.0;
   double cdf=0.0;
   double del=0.0;
//---
   if(y<=thrsh)
     {
      //--- evaluate for  |x| <= 0.66291
      if(y>eps)
         xsq=xx*xx;

      xnum = normal_cdf_a[4] * xsq;
      xden = xsq;
      for(k=0; k<3; k++)
        {
         xnum=(xnum+normal_cdf_a[k])*xsq;
         xden=(xden+normal_cdf_b[k])*xsq;
        }
      cdf = xx*(xnum+normal_cdf_a[3])/(xden+normal_cdf_b[3]);
      cdf = 0.5 + cdf;
     }
   else
   if(y<=root32)
     {
      //--- evaluate for 0.66291 <= |x| <= sqrt(32)
      xnum = normal_cdf_c[8]*y;
      xden = y;
      for(k=0; k<7; k++)
        {
         xnum=(xnum+normal_cdf_c[k])*y;
         xden=(xden+normal_cdf_d[k])*y;
        }
      cdf=(xnum+normal_cdf_c[7])/(xden+normal_cdf_d[7]);
      xsq=int(y*16)/16;
      del=(y-xsq)*(y+xsq);
      cdf=MathExp(-xsq*xsq*0.5)*MathExp(-del*0.5)*cdf;
      if(xx>0.0) cdf=1.0-cdf;
     }
//--- evaluate for |x| > sqrt(32)
   else
     {
      cdf=0.0;
      if((xx>=xlow) && (xx<xuppr))
        {
         xsq=1.0/(xx*xx);
         xnum = normal_cdf_p[5]*xsq;
         xden = xsq;
         for(k=0; k<3; k++)
           {
            xnum=(xnum+normal_cdf_p[k])*xsq;
            xden=(xden+normal_cdf_q[k])*xsq;
           }
         cdf=xsq*(xnum+normal_cdf_p[4])/(xden+normal_cdf_q[4]);
         cdf=(sqrpi-cdf)/y;
         xsq=int(xx*16)/16;
         del=(xx-xsq)*(xx+xsq);
         cdf=MathExp(-xsq*xsq*0.5)*MathExp(-del*0.5)*cdf;
        }
      if(xx>0.0) cdf=1.0-cdf;
     }
//--- take into account round-off errors for probability
   return TailLogValue(MathMin(cdf,1.0),tail,log_mode);
  }
//+------------------------------------------------------------------+
//| Normal cumulative distribution function (CDF)                    |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Normal distribution with parameters mu and sigma        |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| mu         : Mean                                                |
//| sigma      : Standard deviation (must be positive)               |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Normal cumulative distribution function with    |
//| parameters mu and sigma, evaluated at x.                         |
//+------------------------------------------------------------------+
double MathCumulativeDistributionNormal(const double x,const double mu,const double sigma,int &error_code)
  {
   return MathCumulativeDistributionNormal(x,mu,sigma,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Normal cumulative distribution function (CDF)                    |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Normal distribution with parameters mu and sigma             |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Mean                                                |
//| sigma      : Standard deviation (must be positive)               |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionNormal(const double &x[],const double mu,const double sigma,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
      return false;
//--- check sigma
   if(sigma<=0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      //--- prepare argument
      double xx=(x_arg-mu)/sigma;
      //--- mathematical constants
      //--- sqrpi = 1 / sqrt(2*pi), root32 = sqrt(32), and
      //--- thrsh is the argument for which anorm = 0.75.
      const double sqrpi=1.0/MathSqrt(2*M_PI);
      const double thrsh = 0.66291e0;
      const double root32= MathSqrt(32);
      //--- machine-dependent constants
      //--- data eps/5.96e-8/,xlow/-12.949e0/,xuppr/5.768e0/
      const double eps=1.11e-16;
      const double xlow=-37.519;
      const double xuppr=8.572;
      int k;
      //---
      double xsq=0.0;
      double y=MathAbs(xx);
      double xnum=0.0;
      double xden=0.0;
      double cdf=0.0;
      double del=0.0;
      //---
      if(y<=thrsh)
        {
         //--- evaluate for  |x| <= 0.66291
         if(y>eps)
            xsq=xx*xx;

         xnum = normal_cdf_a[4] * xsq;
         xden = xsq;
         for(k=0; k<3; k++)
           {
            xnum=(xnum+normal_cdf_a[k])*xsq;
            xden=(xden+normal_cdf_b[k])*xsq;
           }
         cdf = xx*(xnum+normal_cdf_a[3])/(xden+normal_cdf_b[3]);
         cdf = 0.5 + cdf;
        }
      else
      if(y<=root32)
        {
         //--- evaluate for 0.66291 <= |x| <= sqrt(32)
         xnum = normal_cdf_c[8]*y;
         xden = y;
         for(k=0; k<7; k++)
           {
            xnum=(xnum+normal_cdf_c[k])*y;
            xden=(xden+normal_cdf_d[k])*y;
           }
         cdf=(xnum+normal_cdf_c[7])/(xden+normal_cdf_d[7]);
         xsq=int(y*16)/16;
         del=(y-xsq)*(y+xsq);
         cdf=MathExp(-xsq*xsq*0.5)*MathExp(-del*0.5)*cdf;
         if(xx>0.0) cdf=1.0-cdf;
        }
      //--- evaluate for |x| > sqrt(32)
      else
        {
         cdf=0.0;
         if((xx>=xlow) && (xx<xuppr))
           {
            xsq=1.0/(xx*xx);
            xnum = normal_cdf_p[5]*xsq;
            xden = xsq;
            for(k=0; k<3; k++)
              {
               xnum=(xnum+normal_cdf_p[k])*xsq;
               xden=(xden+normal_cdf_q[k])*xsq;
              }
            cdf=xsq*(xnum+normal_cdf_p[4])/(xden+normal_cdf_q[4]);
            cdf=(sqrpi-cdf)/y;
            xsq=int(xx*16)/16;
            del=(xx-xsq)*(xx+xsq);
            cdf=MathExp(-xsq*xsq*0.5)*MathExp(-del*0.5)*cdf;
           }
         if(xx>0.0) cdf=1.0-cdf;
        }
      //--- take into account round-off errors for probability
      result[i]=TailLogValue(MathMin(cdf,1.0),tail,log_mode);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Normal cumulative distribution function (CDF)                    |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Normal distribution with parameters mu and sigma             |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Mean                                                |
//| sigma      : Standard deviation (must be positive)               |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionNormal(const double &x[],const double mu,const double sigma,double &result[])
  {
   return MathCumulativeDistributionNormal(x,mu,sigma,true,false,result);
  }
//+------------------------------------------------------------------+
//| Normal distribution quantile function (inverse CDF)              |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Normal distribution with parameters mu and sigma |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| mu          : Mean                                               |
//| sigma       : Standard deviation (must be positive)              |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode,if true it calculates for Log values|
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Normal distribution with parameters mu and sigma.         |
//+------------------------------------------------------------------+
//| Comment from original FORTRAN code                               |
//| http://www1.fpl.fs.fed.us/ni241.f                                |
//| Produces the normal deviate Z corresponding to a given lower     |
//| tail area of P; Z is accurate to about 1 part in 10**16.         |
//| Wichura, M.J. (1988). Algorithm AS 241: The Percentage Points of |
//| the Normal Distribution. Applied Statistics, v.37, N3, 477-484.  |
//+------------------------------------------------------------------+
double MathQuantileNormal(const double probability,const double mu,const double sigma,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(probability) || !MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check sigma  
   if(sigma<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

//--- calculate real probability
   double prob=TailLogProbability(probability,tail,log_mode);
//--- check probability range
   if(prob<0.0 || prob>1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- f(0)=-infinity
   if(prob==0.0)
     {
      error_code=ERR_RESULT_INFINITE;
      return QNEGINF;
     }
//--- f(1)=+infinity
   if(prob==1.0)
     {
      error_code=ERR_RESULT_INFINITE;
      return QPOSINF;
     }

   error_code=ERR_OK;

   double q=prob-0.5;
   double r=0;
   double ppnd16=0.0;
//---
   if(MathAbs(q)<=0.425)
     {
      r=0.180625-q*q;
      ppnd16=q*(((((((normal_q_a7*r+normal_q_a6)*r+normal_q_a5)*r+normal_q_a4)*r+normal_q_a3)*r+normal_q_a2)*r+normal_q_a1)*r+normal_q_a0)/
             (((((((normal_q_b7*r+normal_q_b6)*r+normal_q_b5)*r+normal_q_b4)*r+normal_q_b3)*r+normal_q_b2)*r+normal_q_b1)*r+1.0);
      //---
      error_code=ERR_OK;
      return mu+sigma*ppnd16;
     }
   else
     {
      if(q<0.0)
         r=prob;
      else
         r=1.0-prob;
      //---
      r=MathSqrt(-MathLog(r));
      //---
      if(r<=5.0)
        {
         r=r-1.6;
         ppnd16=(((((((normal_q_c7*r+normal_q_c6)*r+normal_q_c5)*r+normal_q_c4)*r+normal_q_c3)*r+normal_q_c2)*r+normal_q_c1)*r+normal_q_c0)/
                (((((((normal_q_d7*r+normal_q_d6)*r+normal_q_d5)*r+normal_q_d4)*r+normal_q_d3)*r+normal_q_d2)*r+normal_q_d1)*r+1.0);
        }
      else
        {
         r=r-5.0;
         ppnd16=(((((((normal_q_e7*r+normal_q_e6)*r+normal_q_e5)*r+normal_q_e4)*r+normal_q_e3)*r+normal_q_e2)*r+normal_q_e1)*r+normal_q_e0)/
                (((((((normal_q_f7*r+normal_q_f6)*r+normal_q_f5)*r+normal_q_f4)*r+normal_q_f3)*r+normal_q_f2)*r+normal_q_f1)*r+1.0);
        }
      //---
      if(q<0.0)
         ppnd16=-ppnd16;
     }
//--- return rescaled/shifted value
   return mu+sigma*ppnd16;
  }
//+------------------------------------------------------------------+
//| Normal distribution quantile function (inverse CDF)              |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of Normal distribution with parameters mu and sigma     |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| mu          : Mean                                               |
//| sigma       : Standard deviation (must be positive)              |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Normal distribution with parameters mu and sigma.         |
//+------------------------------------------------------------------+
double MathQuantileNormal(const double probability,const double mu,const double sigma,int &error_code)
  {
   return MathQuantileNormal(probability,mu,sigma,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Normal distribution quantile function (inverse CDF)              |
//+------------------------------------------------------------------+
//| The function calculates  the inverse cumulative distribution     |
//| function of the Normal distribution with parameters mu and sigma |
//| for the probability values from array.                           |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| mu          : Mean                                               |
//| sigma       : Standard deviation (must be positive)              |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileNormal(const double &probability[],const double mu,const double sigma,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
      return false;
//--- check sigma  
   if(sigma<0)
      return false;

   int data_count=ArraySize(probability);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);

//--- case sigma==0
   if(sigma==0.0)
     {
      for(int i=0; i<data_count; i++)
         result[i]=mu;
      return true;
     }

   for(int i=0; i<data_count; i++)
     {
      //--- calculate real probability
      double prob=TailLogProbability(probability[i],tail,log_mode);
      //--- check probability range
      if(prob<0.0 || prob>1.0)
         return false;

      //--- f(0)=-infinity, f(1)=+infinity
      if(prob==0.0 || prob==1.0)
        {
         if(prob==0.0)
            result[i]=QNEGINF;
         else
            result[i]=QPOSINF;
        }
      else
        {
         double q=prob-0.5;
         double r=0;
         double ppnd16=0.0;
         //---
         if(MathAbs(q)<=0.425)
           {
            r=0.180625-q*q;
            ppnd16=q*(((((((normal_q_a7*r+normal_q_a6)*r+normal_q_a5)*r+normal_q_a4)*r+normal_q_a3)*r+normal_q_a2)*r+normal_q_a1)*r+normal_q_a0)/
                   (((((((normal_q_b7*r+normal_q_b6)*r+normal_q_b5)*r+normal_q_b4)*r+normal_q_b3)*r+normal_q_b2)*r+normal_q_b1)*r+1.0);
            //--- set rescaled/shifted value
            result[i]=mu+sigma*ppnd16;
           }
         else
           {
            if(q<0.0)
               r=prob;
            else
               r=1.0-prob;
            //---
            r=MathSqrt(-MathLog(r));
            //---
            if(r<=5.0)
              {
               r=r-1.6;
               ppnd16=(((((((normal_q_c7*r+normal_q_c6)*r+normal_q_c5)*r+normal_q_c4)*r+normal_q_c3)*r+normal_q_c2)*r+normal_q_c1)*r+normal_q_c0)/
                      (((((((normal_q_d7*r+normal_q_d6)*r+normal_q_d5)*r+normal_q_d4)*r+normal_q_d3)*r+normal_q_d2)*r+normal_q_d1)*r+1.0);
              }
            else
              {
               r=r-5.0;
               ppnd16=(((((((normal_q_e7*r+normal_q_e6)*r+normal_q_e5)*r+normal_q_e4)*r+normal_q_e3)*r+normal_q_e2)*r+normal_q_e1)*r+normal_q_e0)/
                      (((((((normal_q_f7*r+normal_q_f6)*r+normal_q_f5)*r+normal_q_f4)*r+normal_q_f3)*r+normal_q_f2)*r+normal_q_f1)*r+1.0);
              }
            //---
            if(q<0.0)
               ppnd16=-ppnd16;
           }
         //--- set rescaled/shifted value
         result[i]=mu+sigma*ppnd16;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Normal distribution quantile function (inverse CDF)              |
//+------------------------------------------------------------------+
//| The function calculates  the inverse cumulative distribution     |
//| function of the Normal distribution with parameters mu and sigma |
//| for the probability values from array.                           |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| mu          : Mean                                               |
//| sigma       : Standard deviation (must be positive)              |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileNormal(const double &probability[],const double mu,const double sigma,double &result[])
  {
   return MathQuantileNormal(probability,mu,sigma,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Normal distribution                      |
//+------------------------------------------------------------------+
//| Compute the random variable from the Normal distribution         |
//| with given mean mu and standard deviation sigma.                 |
//|                                                                  |
//| Arguments:                                                       |
//| mu          : Mean                                               |
//| sigma       : Standard deviation (must be positive)              |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Normal distribution.                       |
//+------------------------------------------------------------------+
double MathRandomNormal(const double mu,const double sigma,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check sigma  
   if(sigma<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//---
   if(sigma==0.0)
      return mu;
//--- generate random number
   double rnd=MathRandomNonZero();
//--- return normal random using quantile
   return MathQuantileNormal(rnd,mu,sigma,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Random variate from the Normal distribution                      |
//+------------------------------------------------------------------+
//| Generates random variables from the Normal distribution with     |
//| parameters mu and sigma.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| mu         : Mean                                                |
//| sigma      : Standard deviation (must be positive)               |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomNormal(const double mu,const double sigma,const int data_count,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
      return false;
//--- check sigma  
   if(sigma<0)
      return false;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   if(sigma==0.0)
     {
      for(int i=0; i<data_count; i++)
         result[i]=mu;
      return true;
     }
   int err_code=0;
   for(int i=0; i<data_count; i++)
      result[i]=MathRandomNonZero();
//--- return normal random array using quantile
   return MathQuantileNormal(result,mu,sigma,result);
  }
//+------------------------------------------------------------------+
//| Normal distribution moments                                      |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of the Normal            |
//| distribution with parameters mu and sigma.                       |
//|                                                                  |
//| Arguments:                                                       |
//| mu         : Mean                                                |
//| sigma      : Standard deviation (sigma>0)                        |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsNormal(const double mu,const double sigma,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- check sigma
   if(sigma<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- calculate moments
   mean    =mu;
   variance=MathPow(sigma,2);
   skewness=0;
   kurtosis=0;
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
