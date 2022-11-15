//+------------------------------------------------------------------+
//|                                                    Lognormal.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"
#include "Normal.mqh"

//+------------------------------------------------------------------+
//| Lognormal density function (PDF)                                 |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Lognormal distribution with parameters mu and sigma.      |
//|                                                                  |
//| f(x,mu,sigma)=[1/(x*sigma*sqrt(2pi)]*exp(-(ln(x)-mu)/(2*sigma^2))|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| mu         : Log mean                                            |
//| sigma      : Log standard deviation                              |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityLognormal(const double x,const double mu,const double sigma,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
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
//--- check x
   if(x<=0.0)
      return TailLog0(true,log_mode);
//--- check case sigma==0
   if(sigma==0)
     {
      if(MathLog(MathAbs(x))==mu)
        {
         error_code=ERR_RESULT_INFINITE;
         return QPOSINF;
        }
      else
         return TailLog0(true,log_mode);
     }
//--- prepare argument
   double y=(MathLog(x)-mu)/sigma;
//--- check argument
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
//--- return lognormal density
   return TailLogValue(M_1_SQRT_2PI*MathExp(-0.5*y*y)/(x*sigma),true,log_mode);
  }
//+------------------------------------------------------------------+
//| Lognormal density function (PDF)                                 |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Lognormal distribution with parameters mu and sigma          |
//| for values in x.                                                 |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Log mean                                            |
//| sigma      : Log standard deviation                              |
//| log_mode   : Logarithm mode flag,if true it calculates Log values|
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityLognormal(const double &x[],const double mu,const double sigma,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
      return false;
//--- check sigma
   if(sigma<0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);

//--- check case sigma==0
   if(sigma==0)
     {
      for(int i=0; i<data_count; i++)
        {
         if(MathLog(MathAbs(x[i]))==mu)
            result[i]=QPOSINF;
         else
            result[i]=TailLog0(true,log_mode);
         return true;
        }
     }

   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(!MathIsValidNumber(x_arg))
         return false;

      //--- check x 
      if(x_arg<=0.0)
         result[i]=TailLog0(true,log_mode);
      else
        {
         //--- prepare argument
         double y=(MathLog(x_arg)-mu)/sigma;
         //--- check argument
         if(!MathIsValidNumber(y))
            return false;
         //--- check overflow
         y=MathAbs(y);
         if(y>=2*MathSqrt(DBL_MAX))
            return false;
         //--- return lognormal density
         result[i]=TailLogValue(M_1_SQRT_2PI*MathExp(-0.5*y*y)/(x_arg*sigma),true,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Lognormal density function (PDF)                                 |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Lognormal distribution with parameters mu and sigma          |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Log mean                                            |
//| sigma      : Log standard deviation                              |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityLognormal(const double &x[],const double mu,const double sigma,double &result[])
  {
   return MathProbabilityDensityLognormal(x,mu,sigma,false,result);
  }
//+------------------------------------------------------------------+
//| Lognormal density function (PDF)                                 |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Lognormal distribution with parameters mu and sigma.      |
//|                                                                  |
//| f(x,mu,sigma)=[1/(x*sigma*sqrt(2pi)]*exp(-(ln(x)-mu)/(2*sigma^2))|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| mu         : Log mean                                            |
//| sigma      : Log standard deviation                              |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityLognormal(const double x,const double mu,const double sigma,int &error_code)
  {
   return MathProbabilityDensityLognormal(x,mu,sigma,false,error_code);
  }
//+------------------------------------------------------------------+
//| Lognormal cumulative distribution function (CDF)                 |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Lognormal distribution with parameters mu and sigma     |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| mu         : Log mean                                            |
//| sigma      : Log standard deviation                              |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Lognormal cumulative distribution function      |
//| with parameters mu and sigma, evaluated at x.                    |
//+------------------------------------------------------------------+
double MathCumulativeDistributionLognormal(const double x,const double mu,const double sigma,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
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
//--- check x
   if(x<=0.0)
      return TailLog0(tail,log_mode);
//--- return lognormal cdf using Normal cdf
   return MathCumulativeDistributionNormal(MathLog(x),mu,sigma,tail,log_mode,error_code);
  }
//+------------------------------------------------------------------+
//| Lognormal cumulative distribution function (CDF)                 |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Lognormal distribution with parameters mu and sigma     |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| mu         : Log mean                                            |
//| sigma      : Log standard deviation                              |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Lognormal cumulative distribution function      |
//| with parameters mu and sigma, evaluated at x.                    |
//+------------------------------------------------------------------+
double MathCumulativeDistributionLognormal(const double x,const double mu,const double sigma,int &error_code)
  {
   return MathCumulativeDistributionLognormal(x,mu,sigma,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Lognormal cumulative distribution function (CDF)                 |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Lognormal distribution with parameters mu and sigma          |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Log mean                                            |
//| sigma      : Log standard deviation                              |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionLognormal(const double &x[],const double mu,const double sigma,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
      return false;
//--- check sigma
   if(sigma<0)
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

      //--- check x
      if(x_arg<=0.0)
         result[i]=TailLog0(tail,log_mode);
      else
      //--- return lognormal cdf using Normal cdf
         result[i]=MathCumulativeDistributionNormal(MathLog(x_arg),mu,sigma,tail,log_mode,error_code);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Lognormal cumulative distribution function (CDF)                 |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Lognormal distribution with parameters mu and sigma          |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Log mean                                            |
//| sigma      : Log standard deviation                              |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionLognormal(const double &x[],const double mu,const double sigma,double &result[])
  {
   return MathCumulativeDistributionLognormal(x,mu,sigma,true,false,result);
  }
//+------------------------------------------------------------------+
//| Lognormal distribution quantile function (inverse CDF)           |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Lognormal distribution with parameters mu        |
//| and sigma for the desired probability.                           |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| mu          : Log mean                                           |
//| sigma       : Log standard deviation                             |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode,if true it calculates for Log values|
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The quantile value of the Lognormal distribution.                |
//+------------------------------------------------------------------+
double MathQuantileLognormal(const double probability,const double mu,const double sigma,const bool tail,const bool log_mode,int &error_code)
  {
   if(log_mode==true && probability==QNEGINF)
     {
      error_code=ERR_OK;
      return 0.0;
     }
//--- check NaN
   if(!MathIsValidNumber(probability) || !MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
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

//--- calculate real probability
   double prob=TailLogProbability(probability,tail,log_mode);
//--- check probability range
   if(prob<0.0 || prob>1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

//--- special cases exp(a+b-+infinity)
   if(prob==0.0 || prob==1.0)
     {
      if(sigma==0.0)
        {
         error_code=ERR_OK;
         return MathExp(mu);
        }
      else
      if(prob==0.0)
        {
         if(sigma>0)
           {
            error_code=ERR_OK;
            return 0.0;
           }
         else
         if(sigma<0)
           {
            error_code=ERR_RESULT_INFINITE;
            return QPOSINF;
           }
        }
      else
        {
         if(sigma<0)
           {
            error_code=ERR_OK;
            return 0.0;
           }
         else
         if(sigma>0)
           {
            error_code=ERR_RESULT_INFINITE;
            return QPOSINF;
           }
        }
     }
//--- return lognormal quantile using Normal distribution
   return MathExp(MathQuantileNormal(prob,mu,sigma,error_code));
  }
//+------------------------------------------------------------------+
//| Lognormal distribution quantile function (inverse CDF)           |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of Lognormal distribution with parameters mu and sigma  |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| mu          : Log mean                                           |
//| sigma       : Log standard deviation                             |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The quantile value of the Lognormal distribution.                |
//+------------------------------------------------------------------+
double MathQuantileLognormal(const double probability,const double mu,const double sigma,int &error_code)
  {
   return MathQuantileLognormal(probability,mu,sigma,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Lognormal distribution quantile function (inverse CDF)           |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Lognormal distribution with parameters mu and    |
//| sigma for values from the probability[] array.                   |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| mu          : Log mean                                           |
//| sigma       : Log standard deviation                             |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode,if true it calculates for Log values|
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileLognormal(const double &probability[],const double mu,const double sigma,const bool tail,const bool log_mode,double &result[])
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
   for(int i=0; i<data_count; i++)
     {
      //--- calculate real probability
      double prob=TailLogProbability(probability[i],tail,log_mode);
      //--- check probability range
      if(prob<0.0 || prob>1.0)
         return false;

      //--- special cases exp(a+b-+infinity)
      if(prob==0.0 || prob==1.0)
        {
         if(sigma==0.0)
            result[i]=MathExp(mu);
         else
         if(prob==0.0)
           {
            if(sigma>0)
               result[i]=0.0;
            else
            if(sigma<0)
               result[i]=QPOSINF;
           }
         else
           {
            if(sigma<0)
               result[i]=0.0;
            else
            if(sigma>0)
               result[i]=QPOSINF;
           }
        }
      else
      //--- calculate lognormal quantile using Normal distribution
         result[i]=MathExp(MathQuantileNormal(prob,mu,sigma,error_code));
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Lognormal distribution quantile function (inverse CDF)           |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Lognormal distribution with parameters mu and    |
//| sigma for values from the probability[] array.                   |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| mu          : Log mean                                           |
//| sigma       : Log standard deviation                             |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileLognormal(const double &probability[],const double mu,const double sigma,double &result[])
  {
   return MathQuantileLognormal(probability,mu,sigma,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Lognormal distribution                   |
//+------------------------------------------------------------------+
//| Computes the random variable from the Lognormal distribution     |
//| with parameters mu and sigma.                                    |
//|                                                                  |
//| Arguments:                                                       |
//| mu          : Log mean                                           |
//| sigma       : Log standard deviation                             |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Lognormal distribution.                    |
//+------------------------------------------------------------------+
double MathRandomLognormal(const double mu,const double sigma,int &error_code)
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
//--- generate random number
   double rnd=MathRandomNonZero();
//--- 
   rnd=MathQuantileNormal(rnd,mu,sigma,true,false,error_code);
   return MathExp(rnd);
  }
//+------------------------------------------------------------------+
//| Random variate from the Lognormal distribution                   |
//+------------------------------------------------------------------+
//| Generates random variables from the Lognormal distribution       |
//| with parameters mu and sigma.                                    |
//|                                                                  |
//| Arguments:                                                       |
//| mu         : Log mean                                            |
//| sigma      : Log standard deviation                              |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomLognormal(const double mu,const double sigma,const int data_count,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
      return false;
//--- check sigma  
   if(sigma<0)
      return false;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   int err_code=0;
   for(int i=0; i<data_count; i++)
      result[i]=MathRandomNonZero();
//--- return normal random array using quantile
   MathQuantileNormal(result,mu,sigma,result);
   return MathExp(result);
  }
//+------------------------------------------------------------------+
//| Lognormal distribution moments                                   |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of the Lognormal         |
//| distribution with parameters mu and sigma.                       |
//|                                                                  |
//| Arguments:                                                       |
//| mu         : Log mean                                            |
//| sigma      : Log standard deviation                              |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsLognormal(const double mu,const double sigma,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
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
   if(sigma<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- sigma squared
   double sigma_sqr=sigma*sigma;
   double exp_sigma_sqr=MathExp(sigma_sqr);
//--- calculate moments
   mean    =MathExp(mu+sigma_sqr*0.5);
   variance=(exp_sigma_sqr-1.0)*MathExp(2*mu+sigma_sqr);
   skewness=MathSqrt(exp_sigma_sqr-1.0)*(exp_sigma_sqr+2.0);
   kurtosis=3*MathPowInt(exp_sigma_sqr,2)+2*MathPowInt(exp_sigma_sqr,3)+MathPowInt(exp_sigma_sqr,4)-3-3;
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
