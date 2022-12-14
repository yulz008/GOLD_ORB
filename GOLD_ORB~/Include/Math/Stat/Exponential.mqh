//+------------------------------------------------------------------+
//|                                                  Exponential.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"

//+------------------------------------------------------------------+
//| Exponential density function (PDF)                               |
//+------------------------------------------------------------------+
//| The function returns the probability density function of         |
//| the Exponential distribution with parameter mu.                  |
//|                    f(x,mu)=(1/mu)*exp(-x/mu)                     |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| mu         : Mean                                                |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityExponential(const double x,const double mu,const bool log_mode,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(x) || !MathIsValidNumber(mu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- mu must be positive
   if(mu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- check x
   if(x<0.0)
      return TailLog0(true,log_mode);
//--- calculate lambda;
   double lambda=1.0/mu;
   if(log_mode==true)
      return MathLog(lambda*MathExp(-x*lambda));
//--- return density
   return lambda*MathExp(-x*lambda);
  }
//+------------------------------------------------------------------+
//| Exponential density function (PDF)                               |
//+------------------------------------------------------------------+
//| The function returns the probability density function of         |
//| the Exponential distribution with parameter mu.                  |
//|                    f(x,mu)=(1/mu)*exp(-x/mu)                     |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| mu         : Mean                                                |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityExponential(const double x,const double mu,int &error_code)
  {
   return MathProbabilityDensityExponential(x,mu,false,error_code);
  }
//+------------------------------------------------------------------+
//| Exponential density function (PDF)                               |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Exponential distribution with parameter mu for values in x.  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Mean                                                |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityExponential(const double &x[],const double mu,const bool log_mode,double &result[])
  {
//--- check parameters
   if(!MathIsValidNumber(mu))
      return false;
//--- mu must be positive
   if(mu<=0.0)
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

      if(x_arg<0.0)
         result[i]=TailLog0(true,log_mode);
      else
        {
         //--- calculate lambda;
         double lambda=1.0/mu;
         if(log_mode==true)
            result[i]=MathLog(lambda*MathExp(-x_arg*lambda));
         else
            result[i]=lambda*MathExp(-x_arg*lambda);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Exponential density function (PDF)                               |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Exponential distribution with parameter mu for values in x.  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityExponential(const double &x[],const double mu,double &result[])
  {
   return MathProbabilityDensityExponential(x,mu,false,result);
  }
//+------------------------------------------------------------------+
//| Exponential cumulative distribution function (CDF)               |
//+------------------------------------------------------------------+
//| The function returns the cumulative distribution functin of the  |
//| Exponential distribution with parameter mu, evaluated at x.      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| mu         : Mean                                                |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Exponential cumulative distribution function    |
//| with parameter mu, evaluated at x.                               |
//+------------------------------------------------------------------+
double MathCumulativeDistributionExponential(const double x,const double mu,const bool tail,const bool log_mode,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(x) || !MathIsValidNumber(mu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check mu
   if(mu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- check x
   if(x<0.0)
      return TailLog0(tail,log_mode);
//--- calculate cdf and take into account round-off errors for probability
   double result=MathMin(1.0-MathExp(-x/mu),1.0);
   return TailLogValue(result,tail,log_mode);
  }
//+------------------------------------------------------------------+
//| Exponential cumulative distribution function (CDF)               |
//+------------------------------------------------------------------+
//| The function returns the cumulative distribution function of the |
//| Exponential distribution with parameter mu, evaluated at x.      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| mu         : Mean                                                |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Exponential cumulative distribution function    |
//| with parameter mu, evaluated at x.                               |
//+------------------------------------------------------------------+
double MathCumulativeDistributionExponential(const double x,const double mu,int &error_code)
  {
   return MathCumulativeDistributionExponential(x,mu,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Exponential cumulative distribution function (CDF)               |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Exponential distribution with parameter mu for values in x[].|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Mean                                                |
//| b          : Scale                                               |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionExponential(const double &x[],const double mu,const bool tail,const bool log_mode,double &result[])
  {
//--- check parameters
   if(!MathIsValidNumber(mu))
      return false;
//--- check mu
   if(mu<=0.0)
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

      if(x_arg<0.0)
         result[i]=TailLog0(tail,log_mode);
      else
        {
         //--- calculate cdf and take into account round-off errors for probability
         double cdf=MathMin(1.0-MathExp(-x_arg/mu),1.0);
         result[i]=TailLogValue(cdf,tail,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Exponential cumulative distribution function (CDF)               |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distirbution function of  |
//| the Exponential distribution with parameter mu for values in x[].|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Mean                                                |
//| b          : Scale                                               |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionExponential(const double &x[],const double mu,double &result[])
  {
   return MathCumulativeDistributionExponential(x,mu,true,false,result);
  }
//+------------------------------------------------------------------+
//| Exponential distribution quantile function (inverse CDF)         |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Exponential distribution with parameter mu       |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| mu          : Mean                                               |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Exponential distribution with parameter mu.               |
//+------------------------------------------------------------------+
double MathQuantileExponential(const double probability,const double mu,const bool tail,const bool log_mode,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(mu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- mu must be positive
   if(mu<=0.0)
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

   error_code=ERR_OK;
//--- check zero probability case
   if(prob==0.0)
      return 0.0;
   else
   if(prob==1.0)
      return QPOSINF;
//--- return quantile
   return -mu*MathLog(1.0-prob);
  }
//+------------------------------------------------------------------+
//| Exponential distribution quantile function (inverse CDF)         |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of Exponential distribution with parameter mu           |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| mu          : Mean                                               |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Exponential distribution with parameter mu.               |
//+------------------------------------------------------------------+
double MathQuantileExponential(const double probability,const double mu,int &error_code)
  {
   return MathQuantileExponential(probability,mu,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Exponential distribution quantile function (inverse CDF)         |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Exponential distribution with parameter mu       |
//| for values from the probability[] array.                         |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| mu          : Mean                                               |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileExponential(const double &probability[],const double mu,const bool tail,const bool log_mode,double &result[])
  {
//--- check parameters
   if(!MathIsValidNumber(mu))
      return false;
//--- mu must be positive
   if(mu<=0.0)
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

      //--- check zero probability case
      if(prob==0.0)
         result[i]=0.0;
      else
      if(prob==1.0)
         result[i]=QPOSINF;
      else
         result[i]=-mu*MathLog(1.0-prob);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Exponential distribution quantile function (inverse CDF)         |
//+------------------------------------------------------------------+
//| The function calculates  the inverse cumulative distribution     |
//| function of the Exponential distribution with parameter mu       |
//| for values from the probability[] array.                         |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| mu          : Mean                                               |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileExponential(const double &probability[],const double mu,double &result[])
  {
   return MathQuantileExponential(probability,mu,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Exponential distribution                 |
//+------------------------------------------------------------------+
//| Compute the random variable from the Exponential distribution    |
//| with parameter mu using simple inversion method.                 |
//|                                                                  |
//| Arguments:                                                       |
//| mu          : Mean                                               |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Exponential distribution.                  |
//|                                                                  |
//| Reference:                                                       |
//| Devroye L. "Non-uniform random variate generation",Springer,1986.|
//+------------------------------------------------------------------+
double MathRandomExponential(const double mu,int &error_code)
  {
//--- check mu
   if(!MathIsValidNumber(mu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- mu must be positive
   if(mu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- generate random number
   double rnd=MathRandomNonZero();
//--- return variate using quantile
   return -mu*MathLog(1.0-rnd);
  }
//+------------------------------------------------------------------+
//| Random variate from the Exponential distribution                 |
//+------------------------------------------------------------------+
//| Generates random variables from the Exponential distribution     |
//| with parameter mu.                                               |
//|                                                                  |
//| Arguments:                                                       |
//| mu         : Mean                                                |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomExponential(const double mu,const int data_count,double &result[])
  {
//--- check mu
   if(!MathIsValidNumber(mu))
      return false;
//--- mu must be positive
   if(mu<=0.0)
      return false;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      //--- generate random number
      double rnd=MathRandomNonZero();
      result[i]=-mu*MathLog(1.0-rnd);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Exponential distribution moments                                 |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of the Exponential       |
//| distribution with parameter mu.                                  |
//|                                                                  |
//| Arguments:                                                       |
//| mu         : Mean                                                |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsExponential(const double mu,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(mu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- mu must be positive
   if(mu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- calculate moments
   mean    =mu;
   variance=mu*mu;
   skewness=2;
   kurtosis=6;
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
