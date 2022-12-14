//+------------------------------------------------------------------+
//|                                                      Weibull.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"

//+------------------------------------------------------------------+
//| Weibull probability density function (PDF)                       |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Weibull distribution with parameters a and b.             |
//|           f(x,a,b)=[(a/b)*(x/b)^(a-1)]*exp(-(x/b)^a)             |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| a          : Shape parameter of the distribution (a>0)           |
//| b          : Scale parameter of the distribution (b>0)           |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityWeibull(const double x,const double a,const double b,const bool log_mode,int &error_code)
  {
//--- f(-infinity)=f(infinity)=0
   if(x==QPOSINF || x==QNEGINF)
     {
      error_code=ERR_OK;
      return TailLog0(true,log_mode);
     }
//--- check parameters
   if(!MathIsValidNumber(x) || !MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- a and b must be positive
   if(a<=0 || b<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- check x
   if(x<=0)
      return TailLog0(true,log_mode);
//--- calculate factor
   double pwr=MathPow(x/b,a-1);
   double pdf=(a/b)*pwr*MathExp(-(x/b)*pwr);
   if(log_mode==true)
      return MathLog(pdf);
//--- return density
   return pdf;
  }
//+------------------------------------------------------------------+
//| Weibull probability density function (PDF)                       |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Weibull distribution with parameters a and b.             |
//|           f(x,a,b)=[(a/b)*(x/b)^(a-1)]*exp(-(x/b)^a)             |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| a          : Shape parameter of the distribution (a>0)           |
//| b          : Scale parameter of the distribution (b>0)           |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityWeibull(const double x,const double a,const double b,int &error_code)
  {
   return MathProbabilityDensityWeibull(x,a,b,false,error_code);
  }
//+------------------------------------------------------------------+
//| Weibull probability density function (PDF)                       |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of the  |
//| Weibull distribution with parameters a and b for values in x[].  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Shape parameter of the distribution (a>0)           |
//| b          : Scale parameter of the distribution (b>0)           |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityWeibull(const double &x[],const double a,const double b,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- a and b must be positive
   if(a<=0 || b<=0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      //--- f(-infinity)=f(infinity)=0
      if(x_arg==QPOSINF || x_arg==QNEGINF)
         result[i]=TailLog0(true,log_mode);
      else
      if(x_arg<=0)
                result[i]=TailLog0(true,log_mode);
      else
        {
         //--- calculate factor
         double pwr=MathPow(x_arg/b,a-1);
         double pdf=(a/b)*pwr*MathExp(-(x_arg/b)*pwr);
         if(log_mode==true)
            result[i]=MathLog(pdf);
         else
            result[i]=pdf;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Weibull probability density function (PDF)                       |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of the  |
//| Weibull distribution with parameters a and b for values in x[].  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Shape parameter of the distribution (a>0)           |
//| b          : Scale parameter of the distribution (b>0)           |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityWeibull(const double &x[],const double a,const double b,double &result[])
  {
   return MathProbabilityDensityWeibull(x,a,b,false,result);
  }
//+------------------------------------------------------------------+
//| Weibull cumulative distribution function (CDF)                   |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Weibull distribution with parameters a and b            |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| a          : Shape parameter of the distribution (a>0)           |
//| b          : Scale parameter of the distribution (b>0)           |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Weibull cumulative distribution function        |
//|                    F(a,b)=1-exp(-(x/b)^a)                        |
//| with parameters a and b, evaluated at x.                         |
//+------------------------------------------------------------------+
double MathCumulativeDistributionWeibull(const double x,const double a,const double b,const bool tail,const bool log_mode,int &error_code)
  {
//--- f(-infinity)=0
   if(x==QNEGINF)
     {
      error_code=ERR_OK;
      return TailLog0(tail,log_mode);
     }
//--- f(+infinity)=1
   if(x==QPOSINF)
     {
      error_code=ERR_OK;
      return TailLog1(tail,log_mode);
     }
//--- check parameters
   if(!MathIsValidNumber(x) || !MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- a and b must be positive
   if(a<=0 || b<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- check x
   if(x<=0)
      return TailLog0(tail,log_mode);
//--- calculate probability and take into account round-off errors
   double cdf=MathMin(1.0-MathExp(-MathPow(x/b,a)),1.0);
   return TailLogValue(cdf,tail,log_mode);
  }
//+------------------------------------------------------------------+
//| Weibull cumulative distribution function (CDF)                   |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Weibull distribution with parameters a and b            |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| a          : Shape parameter of the distribution (a>0)           |
//| b          : Scale parameter of the distribution (b>0)           |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Weibull cumulative distribution function        |
//|                    F(a,b)=1-exp(-(x/b)^a)                        |
//| with parameters a and b, evaluated at x.                         |
//+------------------------------------------------------------------+
double MathCumulativeDistributionWeibull(const double x,const double a,const double b,int &error_code)
  {
   return MathCumulativeDistributionWeibull(x,a,b,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Weibull cumulative distribution function (CDF)                   |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Weibull distribution with parameters a and b for values in x.|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Shape parameter of the distribution (a>0)           |
//| b          : Scale parameter of the distribution (b>0)           |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionWeibull(const double &x[],const double a,const double b,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- a and b must be positive
   if(a<=0 || b<=0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      //--- f(-infinity)=0, f(+infinity)=1
      if(x_arg==QNEGINF)
         result[i]=TailLog0(tail,log_mode);
      else
      //--- f(+infinity)=1
      if(x_arg==QPOSINF)
         result[i]=TailLog1(tail,log_mode);
      else
      //--- check x
      if(x_arg<=0)
                result[i]=TailLog0(tail,log_mode);
      else
        {
         //--- calculate probability and take into account round-off errors
         double cdf=MathMin(1.0-MathExp(-MathPow(x_arg/b,a)),1.0);
         result[i]=TailLogValue(cdf,tail,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Weibull cumulative distribution function (CDF)                   |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Weibull distribution with parameters a and b for values in x.|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Shape parameter of the distribution (a>0)           |
//| b          : Scale parameter of the distribution (b>0)           |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionWeibull(const double &x[],const double a,const double b,double &result[])
  {
   return MathCumulativeDistributionWeibull(x,a,b,true,false,result);
  }
//+------------------------------------------------------------------+
//| Weibull distribution quantile function (inverse CDF)             |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of Weibull distribution                                 |
//|                Q(p,a,b)=b*((-ln(1-p)))^(1/a)                     |
//| with parameters a and b for the desired probability.             |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| a           : Shape parameter of the distribution (a>0)          |
//| b           : Scale parameter of the distribution (b>0)          |
//| tail        : Flag to calculate for lower tail                   |
//| log_mode    : Logarithm mode,if true it calculates for Log values|
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Weibull distribution with parameters a and b.             |
//+------------------------------------------------------------------+
double MathQuantileWeibull(const double probability,const double a,const double b,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- a and b must be positive
   if(a<=0 || b<=0)
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
//--- f(1)=+infinity
   if(prob==1.0)
     {
      error_code=ERR_RESULT_INFINITE;
      return QPOSINF;
     }

   error_code=ERR_OK;
//--- f(0)=0
   if(prob==0.0)
      return 0.0;
//--- return quantile
   return b*MathPow(-MathLog(1.0-prob),1.0/a);
  }
//+------------------------------------------------------------------+
//| Weibull distribution quantile function (inverse CDF)             |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Weibull distribution                             |
//|                Q(p,a,b)=b*((-ln(1-p)))^(1/a)                     |
//| with parameters a and b for the desired probability.             |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| a           : Shape parameter of the distribution (a>0)          |
//| b           : Scale parameter of the distribution (b>0)          |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Weibull distribution with parameters a and b.             |
//+------------------------------------------------------------------+
double MathQuantileWeibull(const double probability,const double a,const double b,int &error_code)
  {
   return MathQuantileWeibull(probability,a,b,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Weibull distribution quantile function (inverse CDF)             |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Weibull distribution with parameters a and b     |
//| for the probability values from array.                           |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| a           : Shape parameter of the distribution (a>0)          |
//| b           : Scale parameter of the distribution (b>0)          |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileWeibull(const double &probability[],const double a,const double b,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- a and b must be positive
   if(a<=0 || b<=0)
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

      //--- f(1)=+infinity
      if(prob==1.0)
         result[i]=QPOSINF;
      //--- f(0)=0
      if(prob==0.0)
         result[i]=0.0;
      else
      //--- calc quantile
         result[i]=b*MathPow(-MathLog(1.0-prob),1.0/a);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Weibull distribution quantile function (inverse CDF)             |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Weibull distribution with parameters a and b     |
//| for values from the probability[] array.                         |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| a           : Shape parameter of the distribution (a>0)          |
//| b           : Scale parameter of the distribution (b>0)          |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileWeibull(const double &probability[],const double a,const double b,double &result[])
  {
   return MathQuantileWeibull(probability,a,b,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Weibull distribution                     |
//+------------------------------------------------------------------+
//| Computes the random variable from the Weibull distribution       |
//| with shape a and scale b.                                        |
//|                                                                  |
//| Arguments:                                                       |
//| a          : Shape parameter of the distribution (a>0)           |
//| b          : Scale parameter of the distribution (b>0)           |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Weibull distribution.                      |
//+------------------------------------------------------------------+
double MathRandomWeibull(const double a,const double b,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- a and b must be positive
   if(a<=0 || b<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- generate random number
   double rnd=MathRandomNonZero();
   return b*MathPow(-MathLog(rnd),1.0/a);
  }
//+------------------------------------------------------------------+
//| Random variate from the Weibull distribution                     |
//+------------------------------------------------------------------+
//| Generates random variables from the Weibull distribution with    |
//| parameters a and b.                                              |
//|                                                                  |
//| Arguments:                                                       |
//| a          : Shape parameter of the distribution (a>0)           |
//| b          : Scale parameter of the distribution (b>0)           |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomWeibull(const double a,const double b,const int data_count,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- a and b must be positive
   if(a<=0 || b<=0)
      return false;

//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      //--- generate random number
      double rnd=MathRandomNonZero();
      result[i]=b*MathPow(-MathLog(rnd),1.0/a);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Weibull distribution moments                                     |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of the Weibull           |
//| distribution with parameters a and b.                            |
//|                                                                  |
//| Arguments:                                                       |
//| a          : Shape parameter of the distribution (a>0)           |
//| b          : Scale parameter of the distribution (b>0)           |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsWeibull(const double a,const double b,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
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
//--- a and b must be positive
   if(a<=0 || b<=0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- Gamma function values
   double g1 = MathGamma(1+1.0/a);
   double g2 = MathGamma(1+2.0/a);
   double g3 = MathGamma(1+3.0/a);
   double g4 = MathGamma(1+4.0/a);
//--- calculate moments
   mean    =b*g1;
   variance=b*b*g2-MathPow(g1,2);
   skewness=(2*g1*g1*g1-3*g1*g2+g3)*MathPow(g2-g1*g1,-1.5);
   kurtosis=(-6*MathPow(g1,4)+12*MathPow(g1,2)*g2-3*MathPow(g2,2)-4*g1*g3+g4)*MathPow(g2-g1*g1,-2);
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
