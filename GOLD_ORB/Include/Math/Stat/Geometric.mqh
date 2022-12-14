//+------------------------------------------------------------------+
//|                                                    Geometric.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"

//+------------------------------------------------------------------+
//| Geometric mass function (PDF)                                    |
//+------------------------------------------------------------------+
//| The function returns the probability mass function of the        |
//| Geometric distribution with parameter p.                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| p          : Probability parameter                               |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability mass evaluated at x.                             |
//+------------------------------------------------------------------+
double MathProbabilityDensityGeometric(const double x,const double p,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check probability
   if(p<=0.0 || p>1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check x
   if(x!=MathRound(x))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(x<0)
      return TailLog0(true,log_mode);

   if(p==1.0)
     {
      if(x==0.0)
         return TailLog1(true,log_mode);
      else
         return TailLog0(true,log_mode);
     }
//--- return geometric density
   return TailLogValue(p*MathPow(1.0-p,x),true,log_mode);
  }
//+------------------------------------------------------------------+
//| Geometric mass function (PDF)                                    |
//+------------------------------------------------------------------+
//| The function returns the probability mass function of            |
//| the Geometric distribution with parameter p.                     |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| p          : Probability parameter                               |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability mass evaluated at x.                             |
//+------------------------------------------------------------------+
double MathProbabilityDensityGeometric(const double x,const double p,int &error_code)
  {
   return MathProbabilityDensityGeometric(x,p,false,error_code);
  }
//+------------------------------------------------------------------+
//| Geometric mass function (PDF)                                    |
//+------------------------------------------------------------------+
//| The function calculates the probability mass function of         |
//| the Geometric distribution with parameter p for values in x[].   |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| p          : Probability parameter                               |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityGeometric(const double &x[],const double p,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(p))
      return false;
//--- check probability
   if(p<=0.0 || p>1.0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);

//--- special case p==1.0
   if(p==1.0)
     {
      for(int i=0; i<data_count; i++)
        {
         if(x[i]==0.0)
            result[i]=TailLog1(true,log_mode);
         else
            result[i]=TailLog0(true,log_mode);
        }
      return true;
     }

   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(!MathIsValidNumber(x_arg))
         return false;

      if(x_arg!=MathRound(x_arg))
         return false;

      if(x_arg<0)
         result[i]=TailLog0(true,log_mode);
      else
        {
         double pdf=p*MathPow(1.0-p,x_arg);
         result[i]=TailLogValue(pdf,true,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Geometric mass function (PDF)                                    |
//+------------------------------------------------------------------+
//| The function calculates the probability mass function of         |
//| the Geometric distribution with parameter p for values in x[].   |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| p          : Probability parameter                               |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityGeometric(const double &x[],const double p,double &result[])
  {
   return MathProbabilityDensityGeometric(x,p,false,result);
  }
//+------------------------------------------------------------------+
//| Geometric cumulative distribution function (CDF)                 |
//+------------------------------------------------------------------+
//| The function returns the cumulative distribution function of     |
//| the Geometric distribution with parameter p.                     |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| p          : Probability parameter                               |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode,if true it calculates Log values     |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Geometric cumulative distribution function      |
//| with parameter p, evaluated at x.                                |
//+------------------------------------------------------------------+
double MathCumulativeDistributionGeometric(const double x,const double p,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check probability range
   if(p<=0.0 || p>1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
   error_code=ERR_OK;
//--- check x
   if(x<0)
      return TailLog0(true,log_mode);
//--- check p
   if(p==1.0)
     {
      if(x==0.0)
         return TailLog1(true,log_mode);
      else
         return TailLog0(true,log_mode);
     }
//--- calculate cdf and take into account round-off errors for probability
   double cdf=1.0-MathPow(1.0-p,x+1.0);
   return TailLogValue(MathMin(cdf,1.0),tail,log_mode);
  }
//+------------------------------------------------------------------+
//| Geometric cumulative distribution function (CDF)                 |
//+------------------------------------------------------------------+
//| The function returns the cumulative distribution function of     |
//| the Geometric distribution with parameter p.                     |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| p          : Probability parameter                               |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Geometric cumulative distribution function      |
//| with parameter p, evaluated at x.                                |
//+------------------------------------------------------------------+
double MathCumulativeDistributionGeometric(const double x,const double p,int &error_code)
  {
   return MathCumulativeDistributionGeometric(x,p,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Geometric cumulative distribution function (CDF)                 |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Geometric distribution with parameter p for values in x[].   |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| p          : Probability parameter                               |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode,if true it calculates Log values     |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionGeometric(const double &x[],const double p,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(p))
      return false;
//--- check probability range
   if(p<=0.0 || p>1.0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);

//--- special case p==1.0
   if(p==1.0)
     {
      for(int i=0; i<data_count; i++)
        {
         if(x[i]==0.0)
            result[i]=TailLog1(true,log_mode);
         else
            result[i]=TailLog0(true,log_mode);
        }
      return true;
     }

   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(!MathIsValidNumber(x_arg))
         return false;

      if(x_arg<0)
         result[i]=TailLog0(true,log_mode);
      else
         result[i]=TailLogValue(MathMin(1.0-MathPow(1.0-p,x_arg+1.0),1.0),tail,log_mode);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Geometric cumulative distribution function (CDF)                 |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Geometric distribution with parameter p for values in x[].   |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| p          : Probability parameter                               |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionGeometric(const double &x[],const double p,double &result[])
  {
   return MathCumulativeDistributionGeometric(x,p,true,false,result);
  }
//+------------------------------------------------------------------+
//| Geometric distribution quantile function (inverse CDF)           |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of Geometric distribution with parameter p for the      |
//| desired probability.                                             |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| p           : Probability parameter                              |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode,if true it calculates for Log values|
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Geometric inverse cumulative distribution       |
//| function with parameter p, evaluated at probability.             |
//+------------------------------------------------------------------+
double MathQuantileGeometric(const double probability,const double p,const bool tail,const bool log_mode,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check p range
   if(p<=0.0 || p>=1.0)
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
//--- +infinity
   if(prob==1.0)
      return QPOSINF;
   if(prob==0.0)
      return 0.0;

   double res=MathCeil(-1.0+MathLog(1.0-prob)/MathLog(1.0-p)-1e-12);
   if(res<0)
      res=0;
//--- return quantile
   return res;
  }
//+------------------------------------------------------------------+
//| Geometric distribution quantile function (inverse CDF)           |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Geometric distribution with parameter p          |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| p           : Probability parameter                              |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Geometric quantile function for probability.    |
//+------------------------------------------------------------------+
double MathQuantileGeometric(const double probability,const double p,int &error_code)
  {
   return MathQuantileGeometric(probability,p,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Geometric distribution quantile function (inverse CDF)           |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Geometric distribution with parameter p          |
//| for values form the probability[] array.                         |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| p           : Probability parameter                              |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileGeometric(const double &probability[],const double p,const bool tail,const bool log_mode,double &result[])
  {
//--- check parameters
   if(!MathIsValidNumber(p))
      return false;
//--- check p range
   if(p<=0.0 || p>=1.0)
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

      //--- +infinity
      if(prob==1.0)
         result[i]=QPOSINF;
      if(prob==0.0)
         result[i]=0.0;
      else
        {
         double res=MathCeil(-1.0+MathLog(1.0-prob)/MathLog(1.0-p)-1e-12);
         if(res<0)
            res=0;
         result[i]=res;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Geometric distribution quantile function (inverse CDF)           |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Geometric distribution with parameter p          |
//| for values from the probability[] array.                         |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| p           : Probability parameter                              |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileGeometric(const double &probability[],const double p,double &result[])
  {
   return MathQuantileGeometric(probability,p,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Geometric distribution                   |
//+------------------------------------------------------------------+
//| Computes the random variable from the Geometric distribution     |
//| with parameter p.                                                |
//|                                                                  |
//| Arguments:                                                       |
//| p          : Probability parameter                               |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Geometric distribution.                    |
//+------------------------------------------------------------------+
double MathRandomGeometric(const double p,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check probability range
   if(p<0.0 || p>1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- generate random number
   double rnd=MathRandomNonZero();
   double res=MathCeil(-1.0+MathLog(rnd)/MathLog(1.0-p)-1e-12);
   if(res<0)
      res=0;
   return res;
  }
//+------------------------------------------------------------------+
//| Random variate from the Geometric distribution                   |
//+------------------------------------------------------------------+
//| Generates random variables from the Geometric distribution with  |
//| parameter p.                                                     |
//|                                                                  |
//| Arguments:                                                       |
//| p          : Probability parameter                               |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomGeometric(const double p,const int data_count,double &result[])
  {
//--- check parameters
   if(!MathIsValidNumber(p))
      return false;
//--- check probability range
   if(p<0.0 || p>1.0)
      return false;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      //--- generate random number
      double rnd=MathRandomNonZero();
      double res=MathCeil(-1.0+MathLog(rnd)/MathLog(1.0-p)-1e-12);
      if(res<0)
         res=0;
      result[i]=res;
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Geometric distribution moments                                   |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of Geometric             |
//| distribution with parameter p.                                   |
//|                                                                  |
//| Arguments:                                                       |
//| p          : Probability parameter                               |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsGeometric(const double p,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- check probability range
   if(p<=0.0 || p>=1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return(false);
     }

   error_code=ERR_OK;
//--- calculate moments
   mean    =(1.0/p)-1;
   variance=(1.0-p)/(p*p);
   skewness=(2.0-p)/MathSqrt(1.0-p);
   kurtosis=(p*p-6*p+6)/(1-p);
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
