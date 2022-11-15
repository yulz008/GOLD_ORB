//+------------------------------------------------------------------+
//|                                                       Cauchy.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"

//+------------------------------------------------------------------+
//| Cauchy density function (PDF)                                    |
//+------------------------------------------------------------------+
//| Computes the value of the Cauchy probability density function    |
//| with parameters a and b at the desired quantile x.               |
//|                                                                  |
//|              f(x,a,b)= 1/(pi*b*(1.0+((x-a)/b)^2)                 |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| a          : Mean                                                |
//| b          : Scale                                               |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityCauchy(const double x,const double a,const double b,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check scale
   if(b<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- prepare argument
   double y=(x-a)/b;
//--- check result
   if(!MathIsValidNumber(y))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(log_mode==true)
      return -MathLog(M_PI*b*(1.0+y*y));
//--- return Cauchy density
   return 1.0/(M_PI*b*(1.0+y*y));
  }
//+------------------------------------------------------------------+
//| Cauchy density function (PDF)                                    |
//+------------------------------------------------------------------+
//| Computes the value of the Cauchy probability density function    |
//| with parameters a and b at the desired quantile x.               |
//|                                                                  |
//|              f(x,a,b)= 1/(pi*b*(1.0+((x-a)/b)^2)                 |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| a          : Mean                                                |
//| b          : Scale                                               |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityCauchy(const double x,const double a,const double b,int &error_code)
  {
   return MathProbabilityDensityCauchy(x,a,b,false,error_code);
  }
//+------------------------------------------------------------------+
//| Cauchy density function (PDF)                                    |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| Cauchy distribution with parameters a and b for values           |
//| from x[] array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Mean                                                |
//| b          : Scale                                               |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityCauchy(const double &x[],const double a,const double b,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- check scale
   if(b<=0.0)
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
      double y=(x_arg-a)/b;
      if(log_mode==true)
         result[i]=-MathLog(M_PI*b*(1.0+y*y));
      else
         result[i]=(1.0/(M_PI*b*(1.0+y*y)));
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Cauchy density function (PDF)                                    |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| Cauchy distribution with parameters a and b for values           |
//| in x[] array.                                                    |
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
bool MathProbabilityDensityCauchy(const double &x[],const double a,const double b,double &result[])
  {
   return MathProbabilityDensityCauchy(x,a,b,false,result);
  }
//+------------------------------------------------------------------+
//| Cauchy cumulative distribution function (CDF)                    |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Cauchy distribution with parameters a and b             |
//| is less than or equal to x.                                      |
//|              F(x,a,b)=(1/2)+(1/pi)*arctan((x-a)/b)               |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| a          : Mean                                                |
//| b          : Scale                                               |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode flag,if true it calculates Log values|
//| error_code : Variable for error code                             |
//|                                                                  |
//| The value of the Cauchy cumulative distribution function with    |
//| parameters a and b, evaluated at x.                              |
//+------------------------------------------------------------------+
double MathCumulativeDistributionCauchy(const double x,const double a,const double b,const bool tail,const bool log_mode,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(x) || !MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check scale
   if(b<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- calculate argument
   double y=(x-a)/b;
//--- check result
   if(!MathIsValidNumber(y))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
   error_code=ERR_OK;
//--- calculate probability and take into account round-off errors 
   double cdf=0;
   if(y>-1.0)
      cdf=MathMin(0.5+M_1_PI*MathArctan(y),1.0);
   else
      cdf=MathMin(M_1_PI*MathArctan(-1/y),1.0);

   return TailLogValue(cdf,tail,log_mode);
  }
//+------------------------------------------------------------------+
//| Cauchy cumulative distribution function (CDF)                    |
//+------------------------------------------------------------------+
//| The function returns the cumulative distribution function of     |
//| the Cauchy distribution with parameters a and b, evaluated at x. |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| a          : Mean                                                |
//| b          : Scale                                               |
//| error_code : Variable for error code                             |
//|                                                                  |
//| The value of the Cauchy cumulative distribution function with    |
//| parameters a and b, evaluated at x.                              |
//+------------------------------------------------------------------+
double MathCumulativeDistributionCauchy(const double x,const double a,const double b,int &error_code)
  {
   return MathCumulativeDistributionCauchy(x,a,b,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Cauchy cumulative distribution function (CDF)                    |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Cauchy distribution with parameters a and b for values from  |
//| x[] array.                                                       |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Mean                                                |
//| b          : Scale                                               |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode flag,if true it calculates Log values|
//| error_code : Variable for error code                             |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionCauchy(const double &x[],const double a,const double b,const bool tail,const bool log_mode,double &result[])
  {
//--- check parameters
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- check scale
   if(b<=0.0)
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

      //--- calculate argument
      double y=(x_arg-a)/b;
      //--- check result
      if(!MathIsValidNumber(y))
         return false;
      //--- calculate probability and take into account round-off errors
      double cdf=0;
      if(y>-1.0)
         cdf=MathMin(0.5+M_1_PI*MathArctan(y),1.0);
      else
         cdf=MathMin(M_1_PI*MathArctan(-1/y),1.0);

      result[i]=TailLogValue(cdf,tail,log_mode);
     }

   return true;
  }
//+------------------------------------------------------------------+
//| Cauchy cumulative distribution function (CDF)                    |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function     |
//| of the Cauchy distribution with parameters a and b for values    |
//| from x[] array.                                                  |
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
bool MathCumulativeDistributionCauchy(const double &x[],const double a,const double b,double &result[])
  {
   return MathCumulativeDistributionCauchy(x,a,b,true,false,result);
  }
//+------------------------------------------------------------------+
//| Cauchy distribution quantile function (inverse CDF)              |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Cauchy distribution with parameters a and b      |
//| for the desired probability.                                     |
//|                   Q(p,a,b)=a+b*tan*(pi*(p-1/2))                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| a           : Mean                                               |
//| b           : Scale                                              |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Cauchy distribution with parameters a and b.              |
//+------------------------------------------------------------------+
double MathQuantileCauchy(const double probability,const double a,const double b,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(probability) || !MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check scale
   if(b<=0.0)
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

//--- f(1)= + infinity
   if(prob==1.0)
     {
      error_code=ERR_RESULT_INFINITE;
      return QPOSINF;
     }
//--- f(0)= - infinity
   if(prob==0.0)
     {
      error_code=ERR_RESULT_INFINITE;
      return QNEGINF;
     }
   error_code=ERR_OK;
//--- return quantile
   return a+b*MathTan(M_PI*(prob-0.5));
  }
//+------------------------------------------------------------------+
//| Cauchy distribution quantile function (inverse CDF)              |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Cauchy distribution with parameters a and b      |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| a           : Mean                                               |
//| b           : Scale                                              |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Cauchy distribution with parameters a and b.              |
//+------------------------------------------------------------------+
double MathQuantileCauchy(const double probability,const double a,const double b,int &error_code)
  {
   return MathQuantileCauchy(probability,a,b,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Cauchy distribution quantile function (inverse CDF)              |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Cauchy distribution with parameters a and b      |
//| for the probability values from array.                           |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| a           : Mean                                               |
//| b           : Scale                                              |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileCauchy(const double &probability[],const double a,const double b,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- check scale
   if(b<=0.0)
      return false;

   int data_count=ArraySize(probability);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      if(!MathIsValidNumber(probability[i]))
         return false;

      //--- calculate real probability
      double prob=TailLogProbability(probability[i],tail,log_mode);
      //--- check probability range
      if(prob<0.0 || prob>1.0)
         return false;

      //--- f(1)= + infinity
      if(prob==1.0)
         result[i]=QPOSINF;
      else
      //--- f(0)= - infinity
      if(prob==0.0)
         result[i]=QNEGINF;
      else
         result[i]=a+b*MathTan(M_PI*(prob-0.5));
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Cauchy distribution quantile function (inverse CDF)              |
//+------------------------------------------------------------------+
bool MathQuantileCauchy(const double &probability[],const double a,const double b,double &result[])
  {
   return MathQuantileCauchy(probability,a,b,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Cauchy distribution                      |
//+------------------------------------------------------------------+
//| Compute the random variable from the Cauchy distribution         |
//| with parameters a and b.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| a          : Mean                                                |
//| b          : Scale                                               |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Cauchy distribution.                       |
//+------------------------------------------------------------------+
double MathRandomCauchy(const double a,const double b,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check scale
   if(b<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- check scale=0
   if(b==0.0)
      return a;
//--- generate random number
   double rnd=MathRandomNonZero();
//--- return result
   return a+b*MathTan(M_PI*(rnd-0.5));
  }
//+------------------------------------------------------------------+
//| Random variate from the Cauchy distribution                      |
//+------------------------------------------------------------------+
//| Generates random variables from the Cauchy distribution with     |
//| parameters a and b.                                              |
//|                                                                  |
//| Arguments:                                                       |
//| a          : Mean                                                |
//| b          : Scale                                               |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomCauchy(const double a,const double b,const int data_count,double &result[])
  {
//--- check parameters
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- check scale
   if(b<0)
      return false;

//--- prepare output array
   ArrayResize(result,data_count);

//--- check scale=0
   if(b==0.0)
     {
      for(int i=0; i<data_count; i++)
         result[i]=a;
     }
   else
//--- calculate random values
   for(int i=0; i<data_count; i++)
     {
      //--- generate random number
      double rnd=MathRandomNonZero();
      result[i]=a+b*MathTan(M_PI*(rnd-0.5));
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Cauchy distribution moments                                      |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of Cauchy distribution   |
//| with parameters a and b.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| a          : Mean                                                |
//| b          : Scale                                               |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsCauchy(const double a,const double b,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
   error_code=ERR_OK;
//--- set theoretical values for moments (undefined)
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
