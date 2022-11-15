//+------------------------------------------------------------------+
//|                                                      Uniform.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"

//+------------------------------------------------------------------+
//| Uniform probability density function (PDF)                       |
//+------------------------------------------------------------------+
//| The function returns the probability density function of the     |
//| Uniform distribution with parameters a and b.                    |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| a          : Lower endpoint (minimum)                            |
//| b          : Upper endpoint (maximum)                            |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityUniform(const double x,const double a,const double b,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check range
   if(b<=a)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- check ranges
   if(x>=a && x<=b)
      return TailLogValue(1.0/(b-a),true,log_mode);
//--- otherwise 0
   return TailLog0(true,log_mode);
  }
//+------------------------------------------------------------------+
//| Uniform probability density function (PDF)                       |
//+------------------------------------------------------------------+
//| The function returns the probability density function of the     |
//| Uniform distribution with parameters a and b.                    |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| a          : Lower endpoint (minimum)                            |
//| b          : Upper endpoint (maximum)                            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityUniform(const double x,const double a,const double b,int &error_code)
  {
   return MathProbabilityDensityUniform(x,a,b,false,error_code);
  }
//+------------------------------------------------------------------+
//| Uniform probability density function (PDF)                       |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of the  |
//| Uniform distribution with parameters a and b for values in x[].  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Lower endpoint (minimum)                            |
//| b          : Upper endpoint (maximum)                            |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityUniform(const double &x[],const double a,const double b,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- check range
   if(b<=a)
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

      if(x_arg>=a && x_arg<=b)
         result[i]=TailLogValue(1.0/(b-a),true,log_mode);
      else
         result[i]=TailLog0(true,log_mode);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Uniform probability density function (PDF)                       |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of the  |
//| Uniform distribution with parameters a and b for values in x[].  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| mu         : Mean                                                |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityUniform(const double &x[],const double a,const double b,double &result[])
  {
   return MathProbabilityDensityUniform(x,a,b,false,result);
  }
//+------------------------------------------------------------------+
//| Uniform cumulative distribution function (CDF)                   |
//+------------------------------------------------------------------+
//| The function returns the cumulative distribution function        |
//| of the Uniform distribution with parameters a and b.             |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| a          : Lower endpoint (minimum)                            |
//| b          : Upper endpoint (maximum)                            |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode flag,if true it calculates Log values|
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Uniform cumulative distribution function with   |
//| parameters a and b, evaluated at x.                              |
//+------------------------------------------------------------------+
double MathCumulativeDistributionUniform(const double x,const double a,const double b,const bool tail,const bool log_mode,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(x) || !MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check ranges
   if(b<a)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;

   if(x>=a && x<=b)
      return TailLogValue(MathMin((x-a)/(b-a),1.0),tail,log_mode);

   if(x>b)
      return TailLog1(tail,log_mode);
   return TailLog0(tail,log_mode);
  }
//+------------------------------------------------------------------+
//| Uniform cumulative distribution function (CDF)                   |
//+------------------------------------------------------------------+
//| The function returns the cumulative distribution function of     |
//| the Uniform distribution with parameters a and b.                |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| a          : Lower endpoint (minimum)                            |
//| b          : Upper endpoint (maximum)                            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Uniform cumulative distribution function with   |
//| parameters a and b, evaluated at x.                              |
//+------------------------------------------------------------------+
double MathCumulativeDistributionUniform(const double x,const double a,const double b,int &error_code)
  {
   return MathCumulativeDistributionUniform(x,a,b,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Uniform cumulative distribution function (CDF)                   |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Uniform distribution with parameters a and b for values in x.|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Mean                                                |
//| b          : Scale                                               |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode flag,if true it calculates Log values|
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionUniform(const double &x[],const double a,const double b,const bool tail,const bool log_mode,double &result[])
  {
//--- check parameters
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- check ranges
   if(b<a)
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

      if(x_arg>=a && x_arg<=b)
         result[i]=TailLogValue(MathMin((x_arg-a)/(b-a),1.0),tail,log_mode);
      else
        {
         if(x_arg>b)
            result[i]=TailLog1(tail,log_mode);
         else
            result[i]=TailLog0(tail,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Uniform cumulative distribution function (CDF)                   |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Uniform distribution with parameters a and b for values in x.|
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
bool MathCumulativeDistributionUniform(const double &x[],const double a,const double b,double &result[])
  {
   return MathCumulativeDistributionUniform(x,a,b,true,false,result);
  }
//+------------------------------------------------------------------+
//| Uniform distribution quantile function (inverse CDF)             |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Uniform distribution with parameters a and b     |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| a           : Lower endpoint (minimum)                           |
//| b           : Upper endpoint (maximum)                           |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode,if true it calculates for Log values|
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of Uniform distribution with parameters a and b.                 |
//+------------------------------------------------------------------+
double MathQuantileUniform(const double probability,const double a,const double b,const bool tail,const bool log_mode,int &error_code)
  {
   if(log_mode==true)
     {
      if(probability==QNEGINF)
         return 0.0;
     }
//--- check parameters
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check bounds
   if(b<a)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
   if(b==a)
      return a;

//--- calculate real probability
   double prob=TailLogProbability(probability,tail,log_mode);
//--- check probability range
   if(prob<0.0 || prob>1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;

   if(prob==0.0)
      return a;
   else
   if(prob==1.0)
      return b;

//--- return quantile
   return a+prob*(b-a);
  }
//+------------------------------------------------------------------+
//| Uniform distribution quantile function (inverse CDF)             |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Uniform distribution with parameters a and b     |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| a           : Lower endpoint (minimum)                           |
//| b           : Upper endpoint (maximum)                           |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of Uniform  distribution with parameters a and b.                |
//+------------------------------------------------------------------+
double MathQuantileUniform(const double probability,const double a,const double b,int &error_code)
  {
   return MathQuantileUniform(probability,a,b,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Uniform distribution quantile function (inverse CDF)             |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of Uniform distribution with parameters a and b         |
//| for values from the probability[] array.                         |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| a           : Lower endpoint (minimum)                           |
//| b           : Upper endpoint (maximum)                           |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileUniform(const double &probability[],const double a,const double b,const bool tail,const bool log_mode,double &result[])
  {
//--- check parameters
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- check ranges
   if(b<a)
      return false;

   int data_count=ArraySize(probability);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      if(log_mode==true && probability[i]==QNEGINF)
         result[i]=0;
      else
        {
         //--- calculate real probability
         double prob=TailLogProbability(probability[i],tail,log_mode);
         //--- check probability range
         if(prob<0.0 || prob>1.0)
            return false;

         //--- check bounds
         if(b==a)
            result[i]=a;
         else
         if(prob==0.0)
            result[i]=a;
         else
         if(prob==1.0)
            result[i]=b;
         else
         //--- quantile
            result[i]=(a+prob*(b-a));
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Uniform distribution quantile function (inverse CDF)             |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of Uniform distribution with parameters a and b         |
//| for values from the probability[] array.                         |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| a           : Lower endpoint (minimum)                           |
//| b           : Upper endpoint (maximum)                           |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileUniform(const double &probability[],const double a,const double b,double &result[])
  {
   return MathQuantileUniform(probability,a,b,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Uniform distribution                     |
//+------------------------------------------------------------------+
//| Computes the random variable from the Uniform distribution       |
//| with parameters a and b.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| a           : Lower endpoint (minimum)                           |
//| b           : Upper endpoint (maximum)                           |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with uniform distribution.                      |
//+------------------------------------------------------------------+
double MathRandomUniform(const double a,const double b,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check upper bound
   if(b<a)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- check ranges
   if(a==b)
      return a;
//---
   return a+MathRandomNonZero()*(b-a);
  }
//+------------------------------------------------------------------+
//| Random variate from the Uniform distribution                     |
//+------------------------------------------------------------------+
//| Generates random variables from the Uniform distribution with    |
//| parameters a and b.                                              |
//|                                                                  |
//| Arguments:                                                       |
//| a          : Lower endpoint (minimum)                            |
//| b          : Upper endpoint (maximum)                            |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomUniform(const double a,const double b,const int data_count,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- check upper bound
   if(b<a)
      return false;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      //--- generate random number
      double rnd=MathRandomNonZero();
      result[i]=a+rnd*(b-a);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Uniform distribution moments                                     |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of the Uniform           |
//| distribution with parameters a and b.                            |
//|                                                                  |
//| Arguments:                                                       |
//| a          : Lower endpoint (minimum)                            |
//| b          : Upper endpoint (maximum)                            |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsUniform(const double a,const double b,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- check range
   if(b<=a)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- calculate moments
   mean    =0.5*(a+b);
   variance=MathPow(b-a,2)/12;
   skewness=0;
   kurtosis=-3+9.0/5.0;
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
