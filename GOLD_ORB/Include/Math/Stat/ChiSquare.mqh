//+------------------------------------------------------------------+
//|                                                    ChiSquare.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"
#include "Gamma.mqh"

//+------------------------------------------------------------------+
//| Chi-Square density function (PDF)                                |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Chi-Square distribution with parameter nu.                |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| nu         : Degrees of freedom                                  |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityChiSquare(const double x,const double nu,const bool log_mode,int &error_code)
  {
//--- check arguments
   if(!MathIsValidNumber(x) || !MathIsValidNumber(nu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- nu must be positive integer
   if(nu<=0 || nu!=MathRound(nu))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(x<=0.0)
      return TailLog0(true,log_mode);

//--- calculate using Gamma density
   double pdf=MathProbabilityDensityGamma(x,nu*0.5,2.0,error_code);
   if(log_mode==true)
      return MathLog(pdf);
   return pdf;
  }
//+------------------------------------------------------------------+
//| Chi-Square density function (PDF)                                |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Chi-Square distribution with parameter nu.                |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| nu         : Degrees of freedom                                  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityChiSquare(const double x,const double nu,int &error_code)
  {
   return MathProbabilityDensityChiSquare(x,nu,false,error_code);
  }
//+------------------------------------------------------------------+
//| Chi-Square density function (PDF)                                |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of the  |
//| ChiSquare distribution with parameter nu for values in x[] array.|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+ 
bool MathProbabilityDensityChiSquare(const double &x[],const double nu,const bool log_mode,double &result[])
  {
//--- check arguments
   if(!MathIsValidNumber(nu))
      return false;
//--- nu must be positive integer
   if(nu<=0 || nu!=MathRound(nu))
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

      if(x_arg<=0.0)
         result[i]=TailLog0(true,log_mode);
      else
        {
         //--- calculate using Gamma density
         double pdf=MathProbabilityDensityGamma(x_arg,nu*0.5,2.0,error_code);
         if(log_mode==true)
            result[i]=MathLog(pdf);
         else
            result[i]=pdf;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Chi-Square density function (PDF)                                |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of the  |
//| ChiSquare distribution with parameter nu for values in x[] array.|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+ 
bool MathProbabilityDensityChiSquare(const double &x[],const double nu,double &result[])
  {
   return MathProbabilityDensityChiSquare(x,nu,false,result);
  }
//+------------------------------------------------------------------+
//| Chi-Square cumulative distribution function (CDF)                |
//+------------------------------------------------------------------+
//| The function returns the cumulative distribution function of the |
//| Chi-Square distribution with given nu, evaluated at x.           |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| nu         : Degrees of freedom                                  |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of Chi-Square cumulative distribution function with    |
//| parameter nu, evaluated at x.                                    |
//+------------------------------------------------------------------+
double MathCumulativeDistributionChiSquare(const double x,const double nu,const bool tail,const bool log_mode,int &error_code)
  {
//--- check x
   if(!MathIsValidNumber(x) || !MathIsValidNumber(nu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- nu must be positive integer
   if(nu<=0 || nu!=MathRound(nu))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(x<=0.0)
      return TailLog0(true,log_mode);
//---- calculate using Gamma distribution
   return MathCumulativeDistributionGamma(x,nu*0.5,2.0,tail,log_mode,error_code);
  }
//+------------------------------------------------------------------+
//| Chi-Square cumulative distribution function (CDF)                |
//+------------------------------------------------------------------+
//| The function returns the cumulative distribution function of the |
//| Chi-Square distribution with given nu, evaluated at x.           |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| nu         : Degrees of freedom                                  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of Chi-Square cumulative distribution function with    |
//| parameter nu, evaluated at x.                                    |
//+------------------------------------------------------------------+
double MathCumulativeDistributionChiSquare(const double x,const double nu,int &error_code)
  {
   return MathCumulativeDistributionChiSquare(x,nu,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Chi-Square cumulative distribution function (CDF)                |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Chi-Square distribution with parameter nu for values in x[]. |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionChiSquare(const double &x[],const double nu,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu))
      return false;
//--- nu must be positive integer
   if(nu<=0 || nu!=MathRound(nu))
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

      if(x_arg<=0.0)
         result[i]=TailLog0(true,log_mode);
      else
        {
         double cdf=MathCumulativeDistributionGamma(x_arg,nu*0.5,2.0,true,false,error_code);
         result[i]=TailLogValue(cdf,tail,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Chi-Square cumulative distribution function (CDF)                |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Chi-Square distribution with parameter nu for values in x[]. |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionChiSquare(const double &x[],const double nu,double &result[])
  {
   return MathCumulativeDistributionChiSquare(x,nu,true,false,result);
  }
//+------------------------------------------------------------------+
//| Chi-Square distribution quantile function (inverse CDF)          |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Chi-Square distribution with parameter nu        |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| nu          : Degrees of freedom                                 |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode,if true it calculates for Log values|
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Chi-Square distribution with parameter nu.                |
//+------------------------------------------------------------------+
double MathQuantileChiSquare(const double probability,const double nu,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(nu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- nu must be positive
   if(nu<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- nu must be integer
   if(nu!=MathRound(nu))
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
   if(prob==0.0)
      return 0.0;

   if(prob==1.0)
      return QPOSINF;

//---- calculate quantile using Gamma distribution
   return MathQuantileGamma(prob,nu*0.5,2.0,error_code);
  }
//+------------------------------------------------------------------+
//| Chi-Square distribution quantile function (inverse CDF)          |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Chi-Square distribution with parameter nu        |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| nu          : Degrees of freedom                                 |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Chi-Square distribution with parameter nu.                |
//+------------------------------------------------------------------+
double MathQuantileChiSquare(const double probability,const double nu,int &error_code)
  {
   return MathQuantileChiSquare(probability,nu,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Chi-Square distribution quantile function (inverse CDF)          |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Chi-Square distribution with parameter nu        |
//| for values from the probability[] array.                         |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| nu          : Degrees of freedom                                 |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileChiSquare(const double &probability[],const double nu,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu))
      return false;
//--- nu must be positive
   if(nu<=0)
      return false;
//--- nu must be integer
   if(nu!=MathRound(nu))
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

      if(prob==0.0)
         result[i]=0.0;
      else
      if(prob==1.0)
         result[i]=QPOSINF;
      else
        {
         //--- calculate using Gamma distribution
         result[i]=MathQuantileGamma(prob,nu*0.5,2.0,error_code);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Chi-Square distribution quantile function (inverse CDF)          |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Chi-Square distribution with parameter nu        |
//| for values from the probability[] array.                         |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| nu          : Degrees of freedom                                 |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileChiSquare(const double &probability[],const double nu,double &result[])
  {
   return MathQuantileChiSquare(probability,nu,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Chi-Square distribution                  |
//+------------------------------------------------------------------+
//| Computes the random variable from the Chi-Square distribution    |
//| with parameter nu.                                               |
//|                                                                  |
//| Arguments:                                                       |
//| nu          : Degrees of freedom                                 |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Chi-Square distribution.                   |
//+------------------------------------------------------------------+
double MathRandomChiSquare(const double nu,int &error_code)
  {
//--- NaN
   if(!MathIsValidNumber(nu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- nu must be integer
   if(nu!=MathRound(nu))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- nu must be positive
   if(nu<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- return gamma(nu/2,2)
   return MathRandomGamma(nu*0.5,2.0,error_code);
  }
//+------------------------------------------------------------------+
//| Random variate from Chi-Square distribution                      |
//+------------------------------------------------------------------+
//| Generates random variables from the Chi-Square distribution      |
//| with parameter nu.                                               |
//|                                                                  |
//| Arguments:                                                       |
//| nu         : Degrees of freedom                                  |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomChiSquare(const double nu,const int data_count,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu))
      return false;
//--- nu must be integer
   if(nu!=MathRound(nu))
      return false;
//--- nu must be positive
   if(nu<=0)
      return false;
   int error_code=0;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      //--- generate Gamma random number
      result[i]=MathRandomGamma(nu*0.5,2.0,error_code);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Chi-Square distribution moments                                  |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of Chi-Square            |
//| distribution with parameter nu.                                  |
//|                                                                  |
//| Arguments:                                                       |
//| nu         : Degrees of freedom                                  |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsChiSquare(const double nu,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(nu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- nu must be positive integer
   if(nu<=0 || nu!=MathRound(nu))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- calculate moments
   mean    =nu;
   variance=2*nu;
   skewness=MathSqrt(8/nu);
   kurtosis=12/nu;
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
