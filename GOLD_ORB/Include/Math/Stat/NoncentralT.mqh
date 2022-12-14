//+------------------------------------------------------------------+
//|                                                  NoncentralT.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"
#include "T.mqh"
#include "Normal.mqh"

//+------------------------------------------------------------------+
//| Noncentral T probability density function (PDF)                  |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Noncentral T distribution with parameters nu and delta.   |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| nu         : Degrees of freedom                                  |
//| delta      : Noncentrality parameter                             |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityNoncentralT(const double x,const double nu,const double delta,const bool log_mode,int &error_code)
  {
//--- return T
   if(delta==0.0)
      return MathProbabilityDensityT(x,nu,error_code);
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(nu) || !MathIsValidNumber(delta))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check nu
   if(nu!=MathRound(nu) || nu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//---
   double nu_1=nu+1.0;
   double nu_1_half=nu_1*0.5;
   double log_nu=MathLog(nu);
   double factor1=MathExp(-0.5*(MathLog(M_PI)+log_nu)-(delta*delta)*0.5-MathGammaLog(nu*0.5)+nu_1_half*log_nu);

   double nu_xx=nu+x*x;
   double log_nu_xx=MathLog(nu_xx);
   double factor2=MathExp(-nu_1_half*log_nu_xx);
//---
   const int max_terms=500;
   double pwr=1.0;
   double pwr_factor=x*delta*M_SQRT2;
   double pwr_nuxx=1.0;
   double pwr_nuxx_factor=1.0/MathSqrt(nu_xx);
   double pwr_gamma=1.0;
   int    j=0;
   double pdf=0.0;
   while(j<max_terms)
     {
      if(j>0)
        {
         pwr_nuxx*=pwr_nuxx_factor;
         pwr_gamma/=j;
         pwr*=pwr_factor;
        }
      double t=pwr*pwr_gamma*pwr_nuxx*MathGamma((nu_1+j)*0.5);
      pdf+=t;
      //--- check precision
      if((t/(pdf+10E-10))<10E-20)
         break;
      j++;
     }
//--- check convergence
   if(j<max_terms)
      return TailLogValue(factor1*factor2*pdf,true,log_mode);
   else
     {
      error_code=ERR_NON_CONVERGENCE;
      return QNaN;
     }
  }
//+------------------------------------------------------------------+
//| Noncentral T probability density function (PDF)                  |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Noncentral T distribution with parameters nu and delta.   |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| nu         : Degrees of freedom                                  |
//| delta      : Noncentrality parameter                             |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityNoncentralT(const double x,const double nu,const double delta,int &error_code)
  {
   return MathProbabilityDensityNoncentralT(x,nu,delta,false,error_code);
  }
//+------------------------------------------------------------------+
//| Noncentral T probability density function (PDF)                  |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Noncentral T distribution with parameters nu and delta       |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| delta      : Noncentrality parameter                             |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityNoncentralT(const double &x[],const double nu,const double delta,const bool log_mode,double &result[])
  {
//--- return T
   if(delta==0.0)
      return MathProbabilityDensityT(x,nu,log_mode,result);
//--- check NaN
   if(!MathIsValidNumber(nu) || !MathIsValidNumber(delta))
      return false;
//--- check nu
   if(nu!=MathRound(nu) || nu<=0.0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;
//--- common factors
   double nu_1=nu+1.0;
   double nu_1_half=nu_1*0.5;
   double log_nu=MathLog(nu);
   double factor1=MathExp(-0.5*(MathLog(M_PI)+log_nu)-(delta*delta)*0.5-MathGammaLog(nu*0.5)+nu_1_half*log_nu);
   double factor_delta=delta*M_SQRT2;
   const int max_terms=500;

   int error_code=0;
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(!MathIsValidNumber(x_arg))
         return false;

      double nu_xx=nu+x_arg*x_arg;
      double log_nu_xx=MathLog(nu_xx);
      double factor2=MathExp(-nu_1_half*log_nu_xx);
      //---
      double pwr=1.0;
      double pwr_factor=x_arg*factor_delta;
      double pwr_nuxx=1.0;
      double pwr_nuxx_factor=1.0/MathSqrt(nu_xx);
      double pwr_gamma=1.0;
      int    j=0;
      double pdf=0.0;
      while(j<max_terms)
        {
         if(j>0)
           {
            pwr_nuxx*=pwr_nuxx_factor;
            pwr_gamma/=j;
            pwr*=pwr_factor;
           }
         double t=pwr*pwr_gamma*pwr_nuxx*MathGamma((nu_1+j)*0.5);
         pdf+=t;
         //--- check precision
         if((t/(pdf+10E-10))<10E-20)
            break;
         j++;
        }
      //--- check convergence
      if(j<max_terms)
         result[i]=TailLogValue(factor1*factor2*pdf,true,log_mode);
      else
         return false;
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncentral T probability density function (PDF)                  |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Noncentral T distribution with parameters nu and delta       |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| delta      : Noncentrality parameter                             |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityNoncentralT(const double &x[],const double nu,const double delta,double &result[])
  {
   return MathProbabilityDensityNoncentralT(x,nu,delta,false,result);
  }
//+------------------------------------------------------------------+
//| Noncentral T cumulative distribution function (CDF)              |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Noncentral T distribution with parameters nu and delta  |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| nu         : Degrees of freedom                                  |
//| delta      : Noncentrality parameter                             |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Noncentral T cumulative distribution function   |
//| with parameters nu and delta, evaluated at x.                    |
//+------------------------------------------------------------------+
//| The cdf is calculated using Algorithm 7.2 from                   |
//| Denise Benton, K. Krishnamoorthy,                                |
//| "Computing discrete mixtures of continuous distributions:        |
//| noncentral chisquare, noncentral t and the distribution          |
//| of the square of the sample multiple correlation coeffcient      |
//| Computational Statistics & Data Analysis 43 (2003) 249-267.      |
//+------------------------------------------------------------------+
double MathCumulativeDistributionNoncentralT(const double x,const double nu,const double delta,const bool tail,const bool log_mode,int &error_code)
  {
//--- return T
   if(delta==0.0)
      return MathCumulativeDistributionT(x,nu,error_code);
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(nu) || !MathIsValidNumber(delta))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check nu (must be positive integer)
   if(nu!=MathRound(nu) || nu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- special case
   if(nu==1.0)
     {
      error_code=ERR_OK;
      return TailLogValue(0.5+MathArctan(x)*M_1_PI,tail,log_mode);
     }

   error_code=ERR_OK;

//--- error tolerance
   const double errtol=10E-25;
//--- maximum number of iterations
   const int max_iterations=1000;

   double xx,del;
   double ptermf,qtermf,ptermb,qtermb,error;
//--- if t<0, then the transformation in (3.2) must be used; 
   if(x<0.0)
     {
      xx=-x;
      del=-delta;
     }
   else
     {
      xx=x;
      del=delta;
     }

   int err_code=0;
//--- compute the normal cdf at (-del):
   double cdf_normal=MathMin(MathCumulativeDistributionNormal(-del,0,1,err_code),1.0);

   if(xx==0.0)
      return TailLogValue(cdf_normal,tail,log_mode);

   double y=xx*xx/(nu+xx*xx);
   double dels=0.5*del*del;
//--- k = integral part of (dels)
   int k=(int)dels;
   double a=k+0.5;
   double c=k+1.0;
   double b=0.5*nu;
//--- initialization to compute the Pk's:
   double pkf = MathExp(-dels+k*MathLog(dels)-MathGammaLog(k+1.0));
   double pkb = pkf;
//--- initialization to compute the Qk's:
   double qkf = MathExp(-dels+k*MathLog(dels)-MathGammaLog(k+1.5));
   double qkb = qkf;
//--- compute the incomplete beta function associated with the Pk:
//--- pbetaf = beta distribution at (y; a; b)
   double pbetaf=MathBetaIncomplete(y,a,b);
   double pbetab=pbetaf;
//--- compute the incomplete beta function associated with the Qk:
//--- qbetaf = beta distribution at (y; c; b)
   double qbetaf=MathBetaIncomplete(y,c,b);
   double qbetab=qbetaf;
//--- initialization to compute the incomplete beta functions associated with the Pi's recursively:
   double pgamf=MathExp(MathGammaLog(a+b-1.0)-MathGammaLog(a)-MathGammaLog(b)+(a-1.0)*MathLog(y)+b*MathLog(1.0-y));
   double pgamb=pgamf*y*(a+b-1.0)/a;
//--- initialization to compute the incomplete beta functions associated with the Qi's recursively:
   double qgamf=MathExp(MathGammaLog(c+b-1.0)-MathGammaLog(c)-MathGammaLog(b)+(c-1.0)*MathLog(y)+b*MathLog(1.0-y));
   double qgamb=qgamf*y*(c+b-1.0)/c;
//--- compute the remainder of the Poisson probabilities:
   double rempois=1.0-pkf;
   double delosq2=del/M_SQRT2;
   double sum=pkf*pbetaf+delosq2*qkf*qbetaf;
   double cons=0.5*(1.0+0.5*MathAbs(delta));
   int j=0;
   for(;;)
     {
      j++;
      pgamf*=y*(a+b+j-2.0)/(a+j-1.0);
      pbetaf-=pgamf;
      pkf*=dels/(k+j);
      ptermf=pkf*pbetaf;
      qgamf*=y*(c+b+j-2.0)/(c+j-1.0);
      qbetaf-=qgamf;
      qkf*=dels/(k+j+0.5);
      qtermf=qkf*qbetaf;
      double term=ptermf+delosq2*qtermf;
      sum+=term;
      error=rempois*cons*pbetaf;
      rempois-=pkf;
      //--- do forward and backward computations k times or until convergence:
      if(j<=k)
        {
         pgamb*=(a-j+1.0)/(y*(a+b-j));
         pbetab+=pgamb;
         pkb=(k-j+1.0)*pkb/dels;
         ptermb=pkb*pbetab;
         qgamb*=(c-j+1.0)/(y*(c+b-j));
         qbetab+=qgamb;
         qkb=(k-j+1.5)*qkb/dels;
         qtermb=qkb*qbetab;
         term=ptermb+delosq2*qtermb;
         sum+=term;
         rempois-=pkb;
         if(rempois<=errtol || j>=max_iterations)
            break;
        }
      else
        {
         if(error<=errtol || j>=max_iterations)
            break;
        }
     }
//--- check convergence   
   if(j<max_iterations)
     {
      double cdf=0.5*sum+cdf_normal;
      //--- if x is negative
      if(x<0)
         cdf=1.0-cdf;
      //--- take into account round-off errors for probability
      return TailLogValue(MathMin(cdf,1.0),tail,log_mode);
     }
   else
     {
      error_code=ERR_NON_CONVERGENCE;
      return QNaN;
     }
  }
//+------------------------------------------------------------------+
//| Noncentral T cumulative distribution function (CDF)              |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Noncentral T distribution with parameters nu and delta  |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| nu         : Degrees of freedom                                  |
//| delta      : Noncentrality parameter                             |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Noncentral T cumulative distribution function   |
//| with parameters nu and delta, evaluated at x.                    |
//+------------------------------------------------------------------+
double MathCumulativeDistributionNoncentralT(const double x,const double nu,const double delta,int &error_code)
  {
   return MathCumulativeDistributionNoncentralT(x,nu,delta,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Noncentral T cumulative distribution function (CDF)              |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Noncentral T distribution with parameters nu and delta       |
//| for values in x.                                                 |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| delta      : Noncentrality parameter                             |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
//| The cdf is calculated using Algorithm 7.2 from                   |
//| Denise Benton, K. Krishnamoorthy,                                |
//| "Computing discrete mixtures of continuous distributions:        |
//| noncentral chisquare, noncentral t and the distribution          |
//| of the square of the sample multiple correlation coeffcient      |
//| Computational Statistics & Data Analysis 43 (2003) 249-267.      |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionNoncentralT(const double &x[],const double nu,const double delta,const bool tail,const bool log_mode,double &result[])
  {
//--- return T
   if(delta==0.0)
      return MathCumulativeDistributionT(x,nu,tail,log_mode,result);
//--- check NaN
   if(!MathIsValidNumber(nu) || !MathIsValidNumber(delta))
      return false;

//--- check nu (must be positive integer)
   if(nu!=MathRound(nu) || nu<=0.0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   ArrayResize(result,data_count);
//--- special case
   if(nu==1.0)
     {
      for(int i=0; i<data_count; i++)
        {
         double x_arg=x[i];
         if(!MathIsValidNumber(x_arg))
            return false;
         result[i]=TailLogValue(0.5+MathArctan(x_arg)*M_1_PI,tail,log_mode);
        }
      return true;
     }

   int err_code=0;
//--- compute the normal cdf at (-del):
   double cdf_normal=MathMin(MathCumulativeDistributionNormal(-delta,0,1,err_code),1.0);
//--- error tolerance
   const double errtol=10E-25;
//--- maximum number of iterations
   const int max_iterations=1000;

   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(!MathIsValidNumber(x_arg))
         return false;

      double xx,del;
      double ptermf,qtermf,ptermb,qtermb,error;
      //--- if t<0, then the transformation in (3.2) must be used; 
      if(x_arg<0.0)
        {
         xx=-x_arg;
         del=-delta;
        }
      else
        {
         xx=x_arg;
         del=delta;
        }

      if(xx==0.0)
         result[i]=TailLogValue(cdf_normal,tail,log_mode);
      else
        {
         double y=xx*xx/(nu+xx*xx);
         double dels=0.5*del*del;
         //--- k = integral part of (dels)
         int k=(int)dels;
         double a=k+0.5;
         double c=k+1.0;
         double b=0.5*nu;
         //--- initialization to compute the Pk's:
         double pkf = MathExp(-dels+k*MathLog(dels)-MathGammaLog(k+1.0));
         double pkb = pkf;
         //--- initialization to compute the Qk's:
         double qkf = MathExp(-dels+k*MathLog(dels)-MathGammaLog(k+1.5));
         double qkb = qkf;
         //--- compute the incomplete beta function associated with the Pk:
         //--- pbetaf = beta distribution at (y; a; b)
         double pbetaf=MathBetaIncomplete(y,a,b);
         double pbetab=pbetaf;
         //--- compute the incomplete beta function associated with the Qk:
         //--- qbetaf = beta distribution at (y; c; b)
         double qbetaf=MathBetaIncomplete(y,c,b);
         double qbetab=qbetaf;
         //--- initialization to compute the incomplete beta functions associated with the Pi's recursively:
         double pgamf=MathExp(MathGammaLog(a+b-1.0)-MathGammaLog(a)-MathGammaLog(b)+(a-1.0)*MathLog(y)+b*MathLog(1.0-y));
         double pgamb=pgamf*y*(a+b-1.0)/a;
         //--- initialization to compute the incomplete beta functions associated with the Qi's recursively:
         double qgamf=MathExp(MathGammaLog(c+b-1.0)-MathGammaLog(c)-MathGammaLog(b)+(c-1.0)*MathLog(y)+b*MathLog(1.0-y));
         double qgamb=qgamf*y*(c+b-1.0)/c;
         //--- compute the remainder of the Poisson probabilities:
         double rempois=1.0-pkf;
         double delosq2=del/M_SQRT2;
         double sum=pkf*pbetaf+delosq2*qkf*qbetaf;
         double cons=0.5*(1.0+0.5*MathAbs(delta));
         int j=0;
         for(;;)
           {
            j++;
            pgamf*=y*(a+b+j-2.0)/(a+j-1.0);
            pbetaf-=pgamf;
            pkf*=dels/(k+j);
            ptermf=pkf*pbetaf;
            qgamf*=y*(c+b+j-2.0)/(c+j-1.0);
            qbetaf-=qgamf;
            qkf*=dels/(k+j+0.5);
            qtermf=qkf*qbetaf;
            double term=ptermf+delosq2*qtermf;
            sum+=term;
            error=rempois*cons*pbetaf;
            rempois-=pkf;
            //--- do forward and backward computations k times or until convergence:
            if(j<=k)
              {
               pgamb*=(a-j+1.0)/(y*(a+b-j));
               pbetab+=pgamb;
               pkb=(k-j+1.0)*pkb/dels;
               ptermb=pkb*pbetab;
               qgamb*=(c-j+1.0)/(y*(c+b-j));
               qbetab+=qgamb;
               qkb=(k-j+1.5)*qkb/dels;
               qtermb=qkb*qbetab;
               term=ptermb+delosq2*qtermb;
               sum+=term;
               rempois-=pkb;
               if(rempois<=errtol || j>=max_iterations)
                  break;
              }
            else
              {
               if(error<=errtol || j>=max_iterations)
                  break;
              }
           }
         //--- check convergence   
         if(j<max_iterations)
           {
            double cdf=0.5*sum+cdf_normal;
            //--- if x is negative
            if(x_arg<0)
               cdf=1.0-cdf;
            //--- take into account round-off errors for probability
            result[i]=TailLogValue(MathMin(cdf,1.0),tail,log_mode);
           }
         else
            return false;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncentral T cumulative distribution function (CDF)              |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Noncentral T distribution with parameters nu and delta       |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| delta      : Noncentrality parameter                             |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionNoncentralT(const double &x[],const double nu,const double delta,double &result[])
  {
   return MathCumulativeDistributionNoncentralT(x,nu,delta,true,false,result);
  }
//+------------------------------------------------------------------+
//| The Noncentral T distribution quantile function (inverse CDF)    |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of Noncentral T distribution with parameters nu and     |
//| delta for the desired probability.                               |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| nu          : Degrees of freedom                                 |
//| delta       : Noncentrality parameter                            |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of Noncentral T-distribution with parameters nu and delta.       |
//+------------------------------------------------------------------+
double MathQuantileNoncentralT(const double probability,const double nu,const double delta,const bool tail,const bool log_mode,int &error_code)
  {
//--- return T if delta==0
   if(delta==0.0)
      return MathQuantileT(probability,nu,error_code);
//--- check NaN
   if(!MathIsValidNumber(nu) || !MathIsValidNumber(delta))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check nu
   if(nu!=MathRound(nu) || nu<0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check delta
   if(delta<0.0)
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

   if(prob==0.0 || prob==1.0)
     {
      error_code=ERR_RESULT_INFINITE;
      if(prob==0.0)
         return QNEGINF;
      else
         return QPOSINF;
     }
//--- coefficients for pdf and cdf
   double sqr_delta_half=delta*delta*0.5;
   double nu_half=nu*0.5;
   double log_sqr_delta_half=MathLog(sqr_delta_half);
   double exp_sqr_delta_half=MathExp(-sqr_delta_half);
   double sqrt_sqr_delta_half=MathSqrt(sqr_delta_half);
   double nu_1=nu+1.0;
   double nu_1_half=nu_1*0.5;
   double log_nu=MathLog(nu);
   double factor1=MathExp(-0.5*(MathLog(M_PI)+log_nu)-sqr_delta_half-MathGammaLog(nu*0.5)+nu_1_half*log_nu);

//--- probability for -delta
   int err_code=0;
   double p0=MathCumulativeDistributionNormal(-delta,0.0,1.0,err_code);
//---
   error_code=ERR_OK;
   double precision=10E-20;
   const int max_iterations=150;
   int iterations=0;
   double x=0.5;
   double h=1.0;
   double h_min=precision;
   const int max_terms=500;
//--- Newton iterations
   while(iterations<max_iterations)
     {
      //--- check convergence 
      if((MathAbs(h)>h_min && MathAbs(h)>MathAbs(h_min*x))==false)
         break;

      //--- calculate pdf
      double x_arg_sqr=x*x;
      double nu_xx=nu+x_arg_sqr;
      double log_nu_xx=MathLog(nu_xx);
      double factor2=MathExp(-nu_1_half*log_nu_xx);
      double pwr=1.0;
      double pwr_factor=x*delta*M_SQRT2;
      double pwr_nuxx=1.0;
      double pwr_nuxx_factor=1.0/MathSqrt(nu_xx);
      double pwr_gamma=1.0;
      int    j=0;
      double pdf=0.0;
      while(j<max_terms)
        {
         if(j>0)
           {
            pwr_nuxx*=pwr_nuxx_factor;
            pwr_gamma/=j;
            pwr*=pwr_factor;
           }
         double t=pwr*pwr_gamma*pwr_nuxx*MathGamma((nu_1+j)*0.5);
         pdf+=t;
         //--- check precision
         if((t/(pdf+10E-10))<10E-20)
            break;
         j++;
        }
      //--- check convergence
      if(j>max_terms)
         return false;
      pdf=factor1*factor2*pdf;

      //--- calculate cdf
      double t=(x*x)/(nu+x*x);
      double sum1 = 0.0;
      double sum2 = 0.0;
      //---
      double pwr_coef1=1.0;
      double pwr_coef2=1.0;
      double fact=1.0;
      j=0;
      if(x!=0)
        {
         while(j<max_terms)
           {
            if(j>0)
              {
               pwr_coef1*=sqrt_sqr_delta_half;
               pwr_coef2*=sqr_delta_half;
               fact/=j;
              }
            double coef=1.0/MathGamma(j*0.5+1.0);
            //--- term1: t between 0 and x
            double t1=pwr_coef1*coef*MathBetaIncomplete(t,(j+1.0)*0.5,nu_half);
            sum1+=t1;
            //--- term2: t between x and -x
            double t2=pwr_coef2*fact*MathBetaIncomplete(t,(j+0.5),nu_half);
            sum2+=t2;
            //--- check precision
            if((MathAbs(t1/(sum1+10E-10))<precision) && (MathAbs(t2/(sum2+10E-10))<precision))
               break;
            j++;
           }
        }
      //--- check convergence   
      if(j>max_terms)
         return false;

      //--- compute probability for positive x
      double cdf=p0+exp_sqr_delta_half*sum1*0.5;
      //--- compute probability for negative x
      if(x<0)
         cdf=cdf-exp_sqr_delta_half*sum2;
      //--- take into account round-off errors for probability
      cdf=MathMin(cdf,1.0);

      //--- calculate ratio
      h=(cdf-prob)/pdf;

      double x_new=x-h;
      if(x_new<0.0)
         x_new=x*0.1;
      else
      if(x_new>1.0)
         x_new=1.0-(1.0-x)*0.1;

      if(MathAbs(x_new-x)<10E-15)
         break;

      x=x_new;
      iterations++;
     }
//--- check convergence
   if(iterations<max_iterations)
      return x;
   else
     {
      error_code=ERR_NON_CONVERGENCE;
      return QNaN;
     }
  }
//+------------------------------------------------------------------+
//| Noncentral T distribution quantile function (inverse CDF)        |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Noncentral T distribution with parameters nu     |
//| and delta for the desired probability.                           |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| nu          : Degrees of freedom                                 |
//| delta       : Noncentrality parameter                            |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of Noncentral T-distribution with parameters nu and delta.       |
//+------------------------------------------------------------------+
double MathQuantileNoncentralT(const double probability,const double nu,const double delta,int &error_code)
  {
   return MathQuantileNoncentralT(probability,nu,delta,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Noncentral T distribution quantile function (inverse CDF)        |
//+------------------------------------------------------------------+
//| The function calculates  the inverse cumulative distribution     |
//| function of the Noncentral T distribution with parameters nu     |
//| and delta for values from the probability[] array.               |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| nu          : Degrees of freedom                                 |
//| delta       : Noncentrality parameter                            |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode,if true it calculates for Log values|
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileNoncentralT(const double &probability[],const double nu,const double delta,const bool tail,const bool log_mode,double &result[])
  {
//--- return T if delta==0
   if(delta==0.0)
      return MathQuantileT(probability,nu,tail,log_mode,result);
//--- check NaN
   if(!MathIsValidNumber(nu) || !MathIsValidNumber(delta))
      return false;
//--- check nu
   if(nu!=MathRound(nu) || nu<0.0)
      return false;
//--- check delta
   if(delta<0.0)
      return false;

   int data_count=ArraySize(probability);
   if(data_count==0)
      return false;

//--- coefficients for pdf and cdf
   double sqr_delta_half=delta*delta*0.5;
   double nu_half=nu*0.5;
   double log_sqr_delta_half=MathLog(sqr_delta_half);
   double exp_sqr_delta_half=MathExp(-sqr_delta_half);
   double sqrt_sqr_delta_half=MathSqrt(sqr_delta_half);
   double nu_1=nu+1.0;
   double nu_1_half=nu_1*0.5;
   double log_nu=MathLog(nu);
   double factor1=MathExp(-0.5*(MathLog(M_PI)+log_nu)-sqr_delta_half-MathGammaLog(nu*0.5)+nu_1_half*log_nu);

//--- probability for -delta
   int err_code=0;
   double p0=MathCumulativeDistributionNormal(-delta,0.0,1.0,err_code);
//---
   double precision=10E-20;
   double h_min=precision;
   const int max_iterations=50;
   const int max_terms=500;

   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      //--- calculate real probability
      double prob=TailLogProbability(probability[i],tail,log_mode);

      //--- check probability range
      if(prob<0.0 || prob>1.0)
         return false;

      if(prob==0.0 || prob==1.0)
        {
         if(prob==0.0)
            result[i]=QNEGINF;
         else
            result[i]=QPOSINF;
        }
      else
        {
         double x=0.5;
         double h=1.0;
         int iterations=0;
         //--- Newton iterations
         while(iterations<max_iterations)
           {
            //--- check convergence 
            if((MathAbs(h)>h_min && MathAbs(h)>MathAbs(h_min*x))==false)
               break;

            //--- calculate pdf and cdf
            //double cdf=MathCumulativeDistributionNoncentralT(x,nu,delta,err_code);
            //double pdf=MathProbabilityDensityNoncentralT(x,nu,delta,err_code);

            //--- calculate pdf
            double x_arg_sqr=x*x;
            double nu_xx=nu+x_arg_sqr;
            double log_nu_xx=MathLog(nu_xx);
            double factor2=MathExp(-nu_1_half*log_nu_xx);
            double pwr=1.0;
            double pwr_factor=x*delta*M_SQRT2;
            double pwr_nuxx=1.0;
            double pwr_nuxx_factor=1.0/MathSqrt(nu_xx);
            double pwr_gamma=1.0;
            int    j=0;
            double pdf=0.0;
            while(j<max_terms)
              {
               if(j>0)
                 {
                  pwr_nuxx*=pwr_nuxx_factor;
                  pwr_gamma/=j;
                  pwr*=pwr_factor;
                 }
               double t=pwr*pwr_gamma*pwr_nuxx*MathGamma((nu_1+j)*0.5);
               pdf+=t;
               //--- check precision
               if((t/(pdf+10E-10))<10E-20)
                  break;
               j++;
              }
            //--- check convergence
            if(j>max_terms)
               return false;

            pdf=factor1*factor2*pdf;

            //--- calculate cdf
            double t=(x*x)/(nu+x*x);
            double sum1 = 0.0;
            double sum2 = 0.0;
            //---
            double pwr_coef1=1.0;
            double pwr_coef2=1.0;
            double fact=1.0;
            j=0;
            if(x!=0)
              {
               while(j<max_terms)
                 {
                  if(j>0)
                    {
                     pwr_coef1*=sqrt_sqr_delta_half;
                     pwr_coef2*=sqr_delta_half;
                     fact/=j;
                    }
                  double coef=1.0/MathGamma(j*0.5+1.0);
                  //--- term1: t between 0 and x
                  double t1=pwr_coef1*coef*MathBetaIncomplete(t,(j+1.0)*0.5,nu_half);
                  sum1+=t1;
                  //--- term2: t between x and -x
                  double t2=pwr_coef2*fact*MathBetaIncomplete(t,(j+0.5),nu_half);
                  sum2+=t2;
                  //--- check precision
                  if((MathAbs(t1/(sum1+10E-10))<precision) && (MathAbs(t2/(sum2+10E-10))<precision))
                     break;
                  j++;
                 }
              }
            //--- check convergence   
            if(j>max_terms)
               return false;

            //--- compute probability for positive x
            double cdf=p0+exp_sqr_delta_half*sum1*0.5;
            //--- compute probability for negative x
            if(x<0)
               cdf=cdf-exp_sqr_delta_half*sum2;
            //--- take into account round-off errors for probability
            cdf=MathMin(cdf,1.0);

            //--- calculate ratio
            h=(cdf-prob)/pdf;

            double x_new=x-h;
            if(x_new<0.0)
               x_new=x*0.1;
            else
            if(x_new>1.0)
               x_new=1.0-(1.0-x)*0.1;

            if(MathAbs(x_new-x)<10E-15)
               break;

            x=x_new;
            iterations++;
           }
         //--- check convergence
         if(iterations<max_iterations)
            result[i]=x;
         else
            return false;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncentral T distribution quantile function (inverse CDF)        |
//+------------------------------------------------------------------+
//| The function calculates  the inverse cumulative distribution     |
//| function of the Noncentral T distribution with parameters nu     |
//| and delta for values from the probability[] array.               |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| nu          : Degrees of freedom                                 |
//| delta       : Noncentrality parameter                            |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileNoncentralT(const double &probability[],const double nu,const double delta,double &result[])
  {
   return MathQuantileNoncentralT(probability,nu,delta,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Noncentral T distribution                |
//+------------------------------------------------------------------+
//| Computes the random variable from the Noncentral T distribution  |
//| with parameters nu and delta.                                    |
//|                                                                  |
//| Arguments:                                                       |
//| nu          : Degrees of freedom                                 |
//| delta       : Noncentrality parameter                            |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Noncentral T distribution.                 |
//+------------------------------------------------------------------+
double MathRandomNoncentralT(const double nu,const double delta,int &error_code)
  {
//--- return T if delta==0
   if(delta==0.0)
      return MathRandomT(nu,error_code);
//--- check NaN
   if(!MathIsValidNumber(nu) || !MathIsValidNumber(delta))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(nu!=MathRound(nu) || nu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check delta
   if(delta<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   int err_code=0;
//--- generate normal and chisquare random variables 
   double rnd_normal=delta+MathRandomNormal(0,1,err_code);
   double rnd_nchi=MathRandomGamma(nu*0.5,2.0,err_code);
//--- calculate ratio
   double rnd_value=0;
   if(rnd_nchi!=0)
      rnd_value=rnd_normal/MathSqrt(rnd_nchi/nu);
   return rnd_value;
  }
//+------------------------------------------------------------------+
//| Random variate from the Noncentral T distribution                |
//+------------------------------------------------------------------+
//| Generates random variables from the Noncentral T distribution    |
//| with parameters nu and delta.                                    |
//|                                                                  |
//| Arguments:                                                       |
//| nu         : Degrees of freedom                                  |
//| delta      : Noncentrality parameter                             |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomNoncentralT(const double nu,const double delta,const int data_count,double &result[])
  {
//--- return T if delta==0
   if(delta==0.0)
      return MathRandomT(nu,data_count,result);
//--- check NaN
   if(!MathIsValidNumber(nu) || !MathIsValidNumber(delta))
      return false;
//--- check arguments
   if(nu!=MathRound(nu) || nu<=0.0)
      return false;
//--- check delta
   if(delta<0)
      return false;

   int err_code=0;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      //--- generate normal and chisquare random variables 
      double rnd_normal=delta+MathRandomNormal(0,1,err_code);
      double rnd_nchi=MathRandomGamma(nu*0.5,2.0,err_code);
      //--- calculate ratio
      double rnd_value=0;
      if(rnd_nchi!=0)
        {
         rnd_value=rnd_normal/MathSqrt(rnd_nchi/nu);
         result[i]=rnd_value;
        }
      else
         return false;
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncentral T distribution moments                                |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of the Noncentral T      |
//| distribution with parameters nu and delta.                       |
//|                                                                  |
//| Arguments:                                                       |
//| nu         : Degrees of freedom                                  |
//| delta      : Noncentrality parameter                             |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
double MathMomentsNoncentralT(const double nu,const double delta,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- if delta==0, calc moments for T
   if(delta==0)
      return(MathMomentsT(nu,mean,variance,skewness,kurtosis,error_code));
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(nu) || !MathIsValidNumber(delta))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- check nu
   if(nu!=MathRound(nu) || nu<0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }
//--- check delta
   if(delta<0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- calculate moments
//--- mean
   if(nu>1)
      mean=delta*MathSqrt(nu)*MathGamma((nu-1)*0.5)/(MathSqrt(2)*MathGamma(nu*0.5));
//--- delta^2
   double delta_sqr=delta*delta;
//--- 1/((nu-3)*(nu-2))
   double nu32=1/((nu-3)*(nu-2));
   if(nu>2)
      variance=((delta_sqr+1)*nu)/(nu-2)-MathPow(mean,2);
//--- skewness
   if(nu>3)
     {
      skewness=-2*variance;
      skewness+= nu*(delta_sqr+2*nu-3)*nu32;
      skewness*= mean*MathPow(variance,-1.5);
     }
//--- kurtosis
   if(nu>4)
     {
      kurtosis=-3*variance;
      kurtosis+= nu*(delta_sqr*(nu+1)+3*(3*nu-5))*nu32;
      kurtosis*= -MathPow(mean,2);
      kurtosis+= MathPow(nu,2)*(MathPow(delta,4)+6*delta_sqr+3)/((nu-4)*(nu-2));
      kurtosis*= MathPow(variance,-2);
      kurtosis-=3;
     }
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
