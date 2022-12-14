//+------------------------------------------------------------------+
//|                                                     Logistic.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"

//+------------------------------------------------------------------+
//| Logistic distribution density function (PDF)                     |
//+------------------------------------------------------------------+
//| The function returns the probability density function of         |
//| the Logistic distribution with parameters mu and sigma.          |
//| f(x,mu,sigma)=exp[-(x-mu)/sigma]/(sigma*(exp[-(x-mu)/sigma])^2)  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| mu         : Mean                                                |
//| sigma      : Scale parameter                                     |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityLogistic(const double x,const double mu,const double sigma,const bool log_mode,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(x) || !MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check sigma
   if(sigma<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- prepare argument
   double y=(x-mu)/sigma;
//--- check result
   if(!MathIsValidNumber(y))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;

//--- calculate exponents
   double e=MathExp(-y);
   double e1=(1+e);
   double pdf=e/(sigma*(e1*e1));
   if(log_mode==true)
      return MathLog(pdf);
//--- return logistic density
   return pdf;
  }
//+------------------------------------------------------------------+
//| Logistic distribution density function (PDF)                     |
//+------------------------------------------------------------------+
//| The function returns the probability density function of         |
//| the Logistic distribution with parameters mu and sigma.          |
//| f(x,mu,sigma)=exp[-(x-mu)/sigma]/(sigma*(exp[-(x-mu)/sigma])^2)  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| mu         : Mean                                                |
//| sigma      : Scale parameter                                     |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityLogistic(const double x,const double mu,const double sigma,int &error_code)
  {
   return MathProbabilityDensityLogistic(x,mu,sigma,false,error_code);
  }
//+------------------------------------------------------------------+
//| Logistic distribution density function (PDF)                     |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Logistic distribution with parameters mu and sigma           |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Mean                                                |
//| sigma      : Scale parameter                                     |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityLogistic(const double &x[],const double mu,const double sigma,const bool log_mode,double &result[])
  {
//--- check parameters
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
      return false;
//--- check sigma
   if(sigma<=0.0)
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

      //--- prepare argument
      double y=(x_arg-mu)/sigma;
      //--- check result
      if(!MathIsValidNumber(y))
         return false;

      //--- calculate exponents
      double e=MathExp(-y);
      double e1=(1+e);
      double pdf=e/(sigma*(e1*e1));
      if(log_mode==true)
         result[i]=MathLog(pdf);
      else
         result[i]=pdf;
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Logistic distribution density function (PDF)                     |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Logistic distribution with parameters mu and sigma for       |
//| values from x[] array.                                           |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Mean                                                |
//| sigma      : Scale parameter                                     |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityLogistic(const double &x[],const double mu,const double sigma,double &result[])
  {
   return MathProbabilityDensityLogistic(x,mu,sigma,false,result);
  }
//+------------------------------------------------------------------+
//| Logistic cumulative distribution function (CDF)                  |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Logistic distribution with parameters mu and sigma      |
//| is less than or equal to x.                                      |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| mu         : Mean                                                |
//| sigma      : Scale parameter                                     |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Logistic cumulative distribution function       |
//|           F(x,mu,sigma)=1/(1+exp[-(x-mu)/sigma])                 |
//| with parameters mu and sigma, evaluated at x.                    |
//+------------------------------------------------------------------+
double MathCumulativeDistributionLogistic(const double x,const double mu,double sigma,const bool tail,const bool log_mode,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(x) || !MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check sigma
   if(sigma<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- prepare argument
   double y=(x-mu)/sigma;
//--- check result
   if(!MathIsValidNumber(y))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- calculate cdf and take into account round-off errors for probability
   double result=1.0/(1.0+MathExp(-y));
   return TailLogValue(MathMin(result,1.0),tail,log_mode);
  }
//+------------------------------------------------------------------+
//| Logistic cumulative distribution function (CDF)                  |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Logistic distribution with parameters mu and sigma      |
//| is less than or equal to x.                                      |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| mu         : Mean                                                |
//| sigma      : Scale parameter                                     |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Logistic cumulative distribution function       |
//|           F(x,mu,sigma)=1/(1+exp[-(x-mu)/sigma])                 |
//| with parameters mu and sigma, evaluated at x.                    |
//+------------------------------------------------------------------+
double MathCumulativeDistributionLogistic(const double x,const double mu,double sigma,int &error_code)
  {
   return MathCumulativeDistributionLogistic(x,mu,sigma,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Logistic cumulative distribution function (CDF)                  |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Logistic distribution with parameters mu and sigma for       |
//| values from the x[] array.                                       |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Mean                                                |
//| sigma      : Scale parameter                                     |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionLogistic(const double &x[],const double mu,const double sigma,const bool tail,const bool log_mode,double &result[])
  {
//--- check parameters
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
      return false;
//--- check sigma
   if(sigma<=0.0)
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

      //--- prepare argument
      double y=(x_arg-mu)/sigma;
      //--- check result
      if(!MathIsValidNumber(y))
         return false;
      //--- calculate cdf and take into account round-off errors for probability
      double cdf=MathMin(1.0/(1.0+MathExp(-y)),1.0);
      result[i]=TailLogValue(cdf,tail,log_mode);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Logistic cumulative distribution function (CDF)                  |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Logistic distribution with parameters mu and sigma for       |
//| values from the x[] array.                                       |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Mean                                                |
//| sigma      : Scale parameter                                     |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionLogistic(const double &x[],const double mu,const double sigma,double &result[])
  {
   return MathCumulativeDistributionLogistic(x,mu,sigma,true,false,result);
  }
//+------------------------------------------------------------------+
//| Logistic distribution quantile function (inverse CDF)            |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Logistic distribution with parameters mu         |
//| and sigma for the desired probability.                           |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| mu          : Mean                                               |
//| sigma       : Scale parameter                                    |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode,if true it calculates for Log values|
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//|               Q(p,mu,sigma)= mu+sigma*log(p/(1-p))               |
//| of the Logistic distribution with parameters mu and sigma.       |
//+------------------------------------------------------------------+
double MathQuantileLogistic(const double probability,const double mu,const double sigma,const bool tail,const bool log_mode,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(probability) || !MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
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

   if(prob==0.0 || prob==1.0)
     {
      if(sigma==0.0)
        {
         error_code=ERR_OK;
         return mu;
        }
      else
        {
         error_code=ERR_RESULT_INFINITE;
         if(prob==0.0)
            return QNEGINF;
         else
            return QPOSINF;
        }
     }

   error_code=ERR_OK;
//--- calculate quantile
   double q=MathLog(prob/(1.0-prob));
//--- return rescaled/shifted quantile
   return mu+sigma*q;
  }
//+------------------------------------------------------------------+
//| Logistic distribution quantile function (inverse CDF)            |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Logistic distribution with parameters mu         |
//| and sigma for the desired probability.                           |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| mu          : Mean                                               |
//| sigma       : Scale parameter                                    |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Logistic distribution with parameters mu and sigma.       |
//+------------------------------------------------------------------+
double MathQuantileLogistic(const double probability,const double mu,const double sigma,int &error_code)
  {
   return MathQuantileLogistic(probability,mu,sigma,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Logistic distribution quantile function (inverse CDF)            |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Logistic distribution with parameters mu and     |
//| sigma for values from the probability[] array.                   |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| mu          : Mean                                               |
//| sigma       : Scale parameter                                    |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileLogistic(const double &probability[],const double mu,const double sigma,const bool tail,const bool log_mode,double &result[])
  {
//--- check parameters
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
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

      if(prob==0.0 || prob==1.0)
        {
         if(sigma==0.0)
            result[i]=mu;
         else
           {
            if(prob==0.0)
               result[i]=QNEGINF;
            else
               result[i]=QPOSINF;
           }
        }
      else
        {
         //--- calculate quantile
         double q=MathLog(prob/(1.0-prob));
         //--- rescaled/shifted quantile
         result[i]=mu+sigma*q;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Logistic distribution quantile function (inverse CDF)            |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Logistic distribution with parameters mu and     |
//| sigma for values from the probability[] array.                   |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| mu          : Mean                                               |
//| sigma       : Scale parameter                                    |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileLogistic(const double &probability[],const double mu,const double sigma,double &result[])
  {
   return MathQuantileLogistic(probability,mu,sigma,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Logistic distribution                    |
//+------------------------------------------------------------------+
//| Compute the random variable from the Logistic distribution       |
//| with parameters mu and sigma.                                    |
//|                                                                  |
//| Arguments:                                                       |
//| mu          : Mean                                               |
//| sigma       : Scale parameter                                    |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Logistic distribution.                     |
//+------------------------------------------------------------------+
double MathRandomLogistic(const double mu,const double sigma,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check sigma
   if(sigma<0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- check sigma
   if(sigma==0.0)
      return mu;
//--- generate random number
   double rnd=MathRandomNonZero();
//--- return value
   return mu+sigma*MathLog(rnd/(1.0-rnd));
  }
//+------------------------------------------------------------------+
//| Random variate from the Logistic distribution                    |
//+------------------------------------------------------------------+
//| Generates random variables from the Logistic distribution        |
//| with parameters mu and sigma.                                    |
//|                                                                  |
//| Arguments:                                                       |
//| mu         : Mean                                                |
//| sigma      : Scale parameter                                     |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomLogistic(const double mu,const double sigma,const int data_count,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(mu) || !MathIsValidNumber(sigma))
      return false;
//--- check sigma
   if(sigma<0.0)
      return false;

//--- prepare output array
   ArrayResize(result,data_count);
//--- check sigma
   if(sigma==0.0)
     {
      for(int i=0; i<data_count; i++)
         result[i]=mu;
      return true;
     }
//--- calculate random variables
   for(int i=0; i<data_count; i++)
     {
      //--- generate random number
      double rnd=MathRandomNonZero();
      //--- calculate logistic random number
      result[i]=mu+sigma*MathLog(rnd/(1.0-rnd));
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Logistic distribution moments                                    |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of the Logistic          |
//| distribution with parameters mu and sigma.                       |
//|                                                                  |
//| Arguments:                                                       |
//| mu         : Mean parameter                                      |
//| sigma      : Scale parameter                                     |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsLogistic(const double mu,const double sigma,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
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
   if(sigma<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }
   error_code=ERR_OK;
//--- calculate moments
   mean    =mu;
   variance=MathPow(M_PI*sigma,2)/3.0;
   skewness=0;
   kurtosis=(21.0/5.0)-3;
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
