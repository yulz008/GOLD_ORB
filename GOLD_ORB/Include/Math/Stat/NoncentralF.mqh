//+------------------------------------------------------------------+
//|                                                  NoncentralF.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"
#include "F.mqh"
#include "Gamma.mqh"
#include "NoncentralBeta.mqh"

//+------------------------------------------------------------------+
//| Noncentral-F probability density function (PDF)                  |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Noncentral-F distribution with parameters nu1,nu2,sigma.  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| sigma      : Noncentrality parameter                             |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityNoncentralF(const double x,const double nu1,const double nu2,const double sigma,const bool log_mode,int &error_code)
  {
//--- return F if sigma==0
   if(sigma==0.0)
      return MathProbabilityDensityF(x,nu1,nu2,error_code);
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(nu1) || !MathIsValidNumber(nu2) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(nu1!=MathRound(nu1) || nu2!=MathRound(nu2) || nu1<=0 || nu2<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(x<=0.0)
      return TailLog0(true,log_mode);
//--- factors
   double nu1_half=nu1*0.5;
   double nu2_half=nu2*0.5;
   double nu12_half=nu1_half+nu2_half;
   double lambda=sigma*0.5;
   double coef_lambda=MathExp(-lambda);
   double nu_coef=nu1/nu2;
   double g=x*nu_coef;
   double pwr_g=MathExp((nu1_half-1)*MathLog(g));
   double g1=g+1.0;
   double pwr_g1=MathExp(-nu12_half*MathLog(g1));
   double pwr_lambda=1.0;
   double fact_mult=1.0;
//--- initial value for recurrent calculation   
   double r_beta=MathBeta(nu1_half,nu2_half);
//--- direct calculation of the sum
   int    max_terms=100;
   int    j=0;
   double pdf=0;
   while(j<max_terms)
     {
      if(j>0)
        {
         pwr_g*=g;
         pwr_lambda*=lambda;
         fact_mult/=j;
         pwr_g1/=g1;
         double jm1=j-1;
         r_beta*=((nu1_half+jm1)/(nu12_half+jm1));
        }
      double dp=pwr_g*pwr_g1*coef_lambda*pwr_lambda*fact_mult/r_beta;
      pdf+=dp;
      if(dp/(pdf+10E-10)<10E-14)
         break;
      j++;
     }
//--- check convergence
   if(j<max_terms)
      return TailLogValue(pdf*nu_coef,true,log_mode);
   else
     {
      error_code=ERR_NON_CONVERGENCE;
      return QNaN;
     }
  }
//+------------------------------------------------------------------+
//| Noncentral-F probability density function (PDF)                  |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Noncentral-F distribution with parameters nu1,nu2,sigma.  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| sigma      : Noncentrality parameter                             |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityNoncentralF(const double x,const double nu1,const double nu2,const double sigma,int &error_code)
  {
   return MathProbabilityDensityNoncentralF(x,nu1,nu2,sigma,false,error_code);
  }
//+------------------------------------------------------------------+
//| Noncentral-F probability density function (PDF)                  |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Noncentral F distribution with parameters nu1, nu2 and sigma |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| sigma      : Noncentrality parameter                             |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityNoncentralF(const double &x[],const double nu1,const double nu2,const double sigma,const bool log_mode,double &result[])
  {
//--- return F if sigma==0
   if(sigma==0.0)
      return MathProbabilityDensityF(x,nu1,nu2,log_mode,result);
//--- check NaN
   if(!MathIsValidNumber(nu1) || !MathIsValidNumber(nu2) || !MathIsValidNumber(sigma))
      return false;
//--- check arguments
   if(nu1!=MathRound(nu1) || nu2!=MathRound(nu2) || nu1<=0 || nu2<=0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   const int max_terms=100;
//--- common factors
   double nu1_half=nu1*0.5;
   double nu2_half=nu2*0.5;
   double nu12_half=nu1_half+nu2_half;
   double lambda=sigma*0.5;
   double coef_lambda=MathExp(-lambda);
   double nu_coef=nu1/nu2;
   double r_beta0=MathBeta(nu1_half,nu2_half);
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(x_arg<=0.0)
         result[i]=TailLog0(true,log_mode);
      else
        {
         double g=x_arg*nu_coef;
         double g1=g+1.0;
         //--- initial values for recurrent calculation   
         double pwr_g=MathExp((nu1_half-1)*MathLog(g));
         double pwr_g1=MathExp(-nu12_half*MathLog(g1));
         double pwr_lambda=1.0;
         double fact_mult=1.0;
         double r_beta=r_beta0;
         //--- direct calculation of the sum
         int    j=0;
         double pdf=0;
         while(j<max_terms)
           {
            if(j>0)
              {
               pwr_g*=g;
               pwr_lambda*=lambda;
               fact_mult/=j;
               pwr_g1/=g1;
               double jm1=j-1;
               r_beta*=((nu1_half+jm1)/(nu12_half+jm1));
              }
            double dp=pwr_g*pwr_g1*coef_lambda*pwr_lambda*fact_mult/r_beta;
            pdf+=dp;
            if(dp/(pdf+10E-10)<10E-14)
               break;
            j++;
           }
         //--- check convergence
         if(j<max_terms)
            result[i]=TailLogValue(pdf*nu_coef,true,log_mode);
         else
            return false;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncentral-F probability density function (PDF)                  |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Noncentral F distribution with parameters nu1, nu2 and sigma |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| sigma      : Noncentrality parameter                             |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityNoncentralF(const double &x[],const double nu1,const double nu2,const double sigma,double &result[])
  {
   return MathProbabilityDensityNoncentralF(x,nu1,nu2,sigma,false,result);
  }
//+------------------------------------------------------------------+
//| Noncentral F cumulative distribution function (CDF)              |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from Noncentral F distribution with parameters nu1,nu2,sigma     |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| sigma      : Noncentrality parameter                             |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Noncentral F cumulative distribution function   |
//| with parameters nu1,nu2,sigma, evaluated at x.                   |
//+------------------------------------------------------------------+
double MathCumulativeDistributionNoncentralF(const double x,const double nu1,const double nu2,const double sigma,const bool tail,const bool log_mode,int &error_code)
  {
//--- return F if sigma==0
   if(sigma==0.0)
      return MathCumulativeDistributionF(x,nu1,nu2,error_code);
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(nu1) || !MathIsValidNumber(nu2) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(nu1!=MathRound(nu1) || nu2!=MathRound(nu2) || nu1<=0 || nu2<=0 || x<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(x<=0)
      return TailLog0(tail,log_mode);
//--- calculate cdf using Noncentral Beta distribution
   double arg=(nu1/nu2)*x;
   return MathCumulativeDistributionNoncentralBeta(arg/(1.0+arg),nu1*0.5,nu2*0.5,sigma,tail,log_mode,error_code);
  }
//+------------------------------------------------------------------+
//| Noncentral F cumulative distribution function (CDF)              |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from Noncentral F distribution with parameters nu1,nu2,sigma     |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| sigma      : Noncentrality parameter                             |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Noncentral F cumulative distribution function   |
//| with parameters nu1,nu2,sigma, evaluated at x.                   |
//+------------------------------------------------------------------+
double MathCumulativeDistributionNoncentralF(const double x,const double nu1,const double nu2,const double sigma,int &error_code)
  {
   return MathCumulativeDistributionNoncentralF(x,nu1,nu2,sigma,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Noncentral F cumulative distribution function (CDF)              |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Noncentral Fl distribution with parameters nu1,nu2 and sigma |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| sigma      : Noncentrality parameter                             |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionNoncentralF(const double &x[],const double nu1,const double nu2,const double sigma,const bool tail,const bool log_mode,double &result[])
  {
//--- return F if sigma==0
   if(sigma==0.0)
      return MathCumulativeDistributionF(x,nu1,nu2,tail,log_mode,result);
//--- check NaN
   if(!MathIsValidNumber(nu1) || !MathIsValidNumber(nu2) || !MathIsValidNumber(sigma))
      return false;
//--- check arguments
   if(nu1!=MathRound(nu1) || nu2!=MathRound(nu2) || nu1<=0 || nu2<=0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

//--- common constants
   int error_code=0;
   double nu1_half=nu1*0.5;
   double nu2_half=nu2*0.5;
   ArrayResize(result,data_count);

   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(x_arg<=0)
         result[i]=TailLog0(tail,log_mode);
      else
        {
         //--- calculate cdf using Noncentral Beta distribution
         double arg=(nu1/nu2)*x_arg;
         result[i]=MathCumulativeDistributionNoncentralBeta(arg/(1.0+arg),nu1_half,nu2_half,sigma,tail,log_mode,error_code);
         //--- check result
         if(error_code!=ERR_OK)
            return false;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncentral F cumulative distribution function (CDF)              |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Noncentral Fl distribution with parameters nu1,nu2 and sigma |
//| for values in x.                                                 |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| sigma      : Noncentrality parameter                             |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionNoncentralF(const double &x[],const double nu1,const double nu2,const double sigma,double &result[])
  {
   return MathCumulativeDistributionNoncentralF(x,nu1,nu2,sigma,true,false,result);
  }
//+------------------------------------------------------------------+
//| Noncentral F distribution quantile function (inverse CDF)        |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of Noncentral F distribution with parameters nu1,nu2    |
//| and sigma for the desired probability.                           |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| nu1         : Numerator degrees of freedom                       |
//| nu2         : Denominator degrees of freedom                     |
//| sigma       : Noncentrality parameter                            |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse Noncentral F cumulative distribution    |
//| function with parameters nu1,nu2,sigma, evaluated at x.          |
//+------------------------------------------------------------------+
double MathQuantileNoncentralF(const double probability,const double nu1,const double nu2,const double sigma,const bool tail,const bool log_mode,int &error_code)
  {
   if(log_mode==true && probability==QNEGINF)
      return 0.0;
//--- return F if sigma==0
   if(sigma==0.0)
      return MathQuantileF(probability,nu1,nu2,tail,log_mode,error_code);
//--- check NaN
   if(!MathIsValidNumber(probability) || !MathIsValidNumber(nu1) || !MathIsValidNumber(nu2) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(nu1!=MathRound(nu1) || nu2!=MathRound(nu2) || nu1<=0 || nu2<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check sigma
   if(sigma<0.0)
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

   if(prob==1.0)
     {
      error_code=ERR_RESULT_INFINITE;
      return QPOSINF;
     }
   error_code=ERR_OK;
   if(prob==0.0)
      return 0.0;
//---
   int    max_iterations=50;
   int    iterations=0;
//--- initial values
   double h=1.0;
   double h_min=10E-10;
   double x=0.5;
   int    err_code=0;
//--- Newton iterations
   while(iterations<max_iterations)
     {
      //--- check convegence     
      if((MathAbs(h)>h_min && MathAbs(h)>MathAbs(h_min*x))==false)
         break;
      //--- calculate pdf and cdf
      double pdf=MathProbabilityDensityNoncentralF(x,nu1,nu2,sigma,err_code);
      double cdf=MathCumulativeDistributionNoncentralF(x,nu1,nu2,sigma,err_code);
      //--- calculate ratio
      h=(cdf-prob)/pdf;
      //---
      double x_new=x-h;
      //--- check x
      if(x_new<0.0)
         x_new=x*0.1;
      else
      if(x_new>1.0)
         x_new=1.0-(1.0-x)*0.1;

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
//| Noncentral F distribution quantile function (inverse CDF)        |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of Noncentral F distribution with parameters nu1,nu2    |
//| and sigma for the desired probability.                           |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| nu1         : Numerator degrees of freedom                       |
//| nu2         : Denominator degrees of freedom                     |
//| sigma       : Noncentrality parameter                            |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse Noncentral F cumulative distribution    |
//| function with parameters nu1,nu2,sigma, evaluated at x.          |
//+------------------------------------------------------------------+
double MathQuantileNoncentralF(const double probability,const double nu1,const double nu2,const double sigma,int &error_code)
  {
   return MathQuantileNoncentralF(probability,nu1,nu2,sigma,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Noncentral F distribution quantile function (inverse CDF)        |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of Noncentral F distribution with parameters nu1,nu2    |
//| and sigma for values from the probability[] array.               |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| nu1         : Numerator degrees of freedom                       |
//| nu2         : Denominator degrees of freedom                     |
//| sigma       : Noncentrality parameter                            |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileNoncentralF(const double &probability[],const double nu1,const double nu2,const double sigma,const bool tail,const bool log_mode,double &result[])
  {
//--- return F if sigma==0
   if(sigma==0.0)
      return MathQuantileF(probability,nu1,nu2,tail,log_mode,result);
//--- check NaN
   if(!MathIsValidNumber(nu1) || !MathIsValidNumber(nu2) || !MathIsValidNumber(sigma))
      return false;
//--- check arguments
   if(nu1!=MathRound(nu1) || nu2!=MathRound(nu2) || nu1<=0 || nu2<=0)
      return false;
//--- check sigma
   if(sigma<0.0)
      return false;

   int data_count=ArraySize(probability);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      //--- calculate real probability
      double prob=TailLogProbability(probability[i],tail,log_mode);
      //--- check probability range
      if(prob<0.0 || prob>1.0)
         return false;

      if(prob==1.0)
         result[i]=QPOSINF;
      else
      if(prob==0.0)
         result[i]=0.0;
      else
        {
         int    max_iterations=50;
         int    iterations=0;
         //--- initial values
         double h=1.0;
         double h_min=10E-10;
         double x=0.5;
         int    err_code=0;
         //--- Newton iterations
         while(iterations<max_iterations)
           {
            //--- check convegence     
            if((MathAbs(h)>h_min && MathAbs(h)>MathAbs(h_min*x))==false)
               break;
            //--- calculate pdf and cdf
            double pdf=MathProbabilityDensityNoncentralF(x,nu1,nu2,sigma,err_code);
            double cdf=MathCumulativeDistributionNoncentralF(x,nu1,nu2,sigma,err_code);
            //--- calculate ratio
            h=(cdf-prob)/pdf;
            //---
            double x_new=x-h;
            //--- check x
            if(x_new<0.0)
               x_new=x*0.1;
            else
            if(x_new>1.0)
               x_new=1.0-(1.0-x)*0.1;

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
//| Noncentral F distribution quantile function (inverse CDF)        |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of Noncentral F distribution with parameters nu1,nu2    |
//| and sigma for values from the probability[] array.               |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| nu1         : Numerator degrees of freedom                       |
//| nu2         : Denominator degrees of freedom                     |
//| sigma       : Noncentrality parameter                            |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileNoncentralF(const double &probability[],const double nu1,const double nu2,const double sigma,double &result[])
  {
   return MathQuantileNoncentralF(probability,nu1,nu2,sigma,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Noncentral F-distribution                |
//+------------------------------------------------------------------+
//| Compute the random variable from the Noncentral F-distribution   |
//| with parameters nu1, nu2 and sigma.                              |
//|                                                                  |
//| Arguments:                                                       |
//| nu1         : Numerator degrees of freedom                       |
//| nu2         : Denominator degrees of freedom                     |
//| sigma       : Noncentrality parameter                            |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Noncentral F-distribution.                 |
//+------------------------------------------------------------------+
double MathRandomNoncentralF(const double nu1,const double nu2,const double sigma,int &error_code)
  {
//--- return F if sigma==0
   if(sigma==0.0)
      return MathRandomF(nu1,nu2,error_code);
//--- check NaN
   if(!MathIsValidNumber(nu1) || !MathIsValidNumber(nu2) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(nu1!=MathRound(nu1) || nu2!=MathRound(nu2) || nu1<=0 || nu2<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check sigma
   if(sigma<0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- calculate using noncentral chisquare and chisquare distributions
   double num=MathRandomNoncentralChiSquare(nu1,sigma,error_code)*nu2;
   double den=MathRandomGamma(nu2*0.5,2.0,error_code)*nu1;
   if(den!=0)
      return num/den;
   else
     {
      error_code=ERR_NON_CONVERGENCE;
      return QNaN;
     }
  }
//+------------------------------------------------------------------+
//| Random variate from the Noncentral F distribution                |
//+------------------------------------------------------------------+
//| Generates random variables from the Noncentral F distribution    |
//| with parameters nu1, nu2 and sigma.                              |
//|                                                                  |
//| Arguments:                                                       |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| sigma      : Noncentrality parameter                             |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomNoncentralF(const double nu1,const double nu2,const double sigma,const int data_count,double &result[])
  {
//--- return F if sigma==0
   if(sigma==0.0)
      return MathRandomF(nu1,nu2,data_count,result);
//--- check NaN
   if(!MathIsValidNumber(nu1) || !MathIsValidNumber(nu2) || !MathIsValidNumber(sigma))
      return false;
//--- check arguments
   if(nu1!=MathRound(nu1) || nu2!=MathRound(nu2) || nu1<=0 || nu2<=0)
      return false;
//--- check sigma
   if(sigma<0.0)
      return false;
   int error_code=0;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      //--- calculate using noncentral chisquare and chisquare distributions
      double num=MathRandomNoncentralChiSquare(nu1,sigma,error_code)*nu2;
      double den=MathRandomGamma(nu2*0.5,2.0,error_code)*nu1;
      if(den!=0)
         result[i]=num/den;
      else
         return false;
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncentral F distribution moments                                |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of the Noncental F       |
//| distribution with parameters nu1,nu2 and sigma.                  |
//|                                                                  |
//| Arguments:                                                       |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| sigma      : Noncentrality parameter                             |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsNoncentralF(const double nu1,const double nu2,const double sigma,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- if sigma==0, calc moments for F
   if(sigma==0)
      return MathMomentsF(nu1,nu2,mean,variance,skewness,kurtosis,error_code);
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(nu1) || !MathIsValidNumber(nu2) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- check arguments
   if(nu1!=MathRound(nu1) || nu2!=MathRound(nu2) || nu1<=0 || nu2<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }
//--- check sigma
   if(sigma<0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- calculate moments
   if(nu2>2)
      mean=nu2*(nu1+sigma)/(nu1*(nu2-2));
//--- variance
   if(nu2>4)
      variance=2*MathPow(nu2/nu1,2)*((nu2-2)*(nu1+2*sigma)+MathPow(nu1+sigma,2))/((nu2-4)*MathPow(nu2-2,2));
//--- factors
   double sigma_sqr=MathPow(sigma,2);
   double sigma_cube=sigma_sqr*sigma;
   double nu12m2=(nu1+nu2-2);
   double nu2p10=(nu2+10);
//--- skewness
   if(nu2>6)
     {
      skewness=2*M_SQRT2*MathSqrt(nu2-4);
      skewness*=(nu12m2*(6*sigma_sqr+(2*nu1+nu2-2)*(3*sigma+nu1))+2*sigma_cube);
      skewness/=(nu2-6);
      skewness/=MathPow(nu12m2*(2*sigma+nu1)+sigma_sqr,1.5);
     }
//--- kurtosis
   if(nu2>8)
     {
      double coef=nu2p10*(MathPow(nu1,2)+nu1*(nu2-2))+4*MathPow(nu2-2,2);
      kurtosis=1;
      kurtosis=3*(nu2-4);
      kurtosis*=(nu12m2*(coef*(4*sigma+nu1)+nu2p10*(4*sigma_cube+2*sigma_sqr*(3*nu1+2*nu2-4)))+nu2p10*MathPow(sigma,4));
      kurtosis/=(nu2-8)*(nu2-6);
      kurtosis/=MathPow((nu12m2*(2*sigma+nu1)+sigma_sqr),2);
      kurtosis-=3;
     }
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
