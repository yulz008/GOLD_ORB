//+------------------------------------------------------------------+
//|                                                            F.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"
#include "Beta.mqh"
#include "ChiSquare.mqh"

//+------------------------------------------------------------------+
//| F-density function (PDF)                                         |
//+------------------------------------------------------------------+
//| The function returns the probability density function of the     |
//| F-distribution with parameters nu1 and nu2.                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityF(const double x,const double nu1,const double nu2,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(nu1) || !MathIsValidNumber(nu2))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(nu1!=MathRound(nu1) || nu1!=MathRound(nu1) || nu1<1 || nu2<1)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- check x
   if(x<=0)
      return TailLog0(true,log_mode);
//--- calculate F density
   double value=MathPow((nu1/nu2),nu1*0.5)*MathPow(x,(nu1-2)*0.5)/MathBeta(nu1*0.5,nu2*0.5);
   value=value*MathPow(1.0+(nu1/nu2)*x,-(nu1+nu2)*0.5);
   if(log_mode==true)
      return MathLog(value);
//--- return F density
   return value;
  }
//+------------------------------------------------------------------+
//| F-density function (PDF)                                         |
//+------------------------------------------------------------------+
//| The function returns the probability density function of the     |
//| F-distribution with parameters nu1 and nu2.                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityF(const double x,const double nu1,const double nu2,int &error_code)
  {
   return MathProbabilityDensityF(x,nu1,nu2,false,error_code);
  }
//+------------------------------------------------------------------+
//| F-density function (PDF)                                         |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of the  |
//| F distribution with parameters nu1 and nu2 for values in x[].    |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityF(const double &x[],const double nu1,const double nu2,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu1) || !MathIsValidNumber(nu2))
      return false;
//--- check arguments
   if(nu1!=MathRound(nu1) || nu1!=MathRound(nu1) || nu1<1 || nu2<1)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(x_arg<=0)
         result[i]=TailLog0(true,log_mode);
      else
        {
         //--- calculate F density
         double value=MathPow((nu1/nu2),nu1*0.5)*MathPow(x_arg,(nu1-2)*0.5)/MathBeta(nu1*0.5,nu2*0.5);
         value=value*MathPow(1.0+(nu1/nu2)*x_arg,-(nu1+nu2)*0.5);
         if(log_mode==true)
            result[i]=MathLog(value);
         else
            result[i]=value;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| F-density function (PDF)                                         |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of the  |
//| F distribution with parameters nu1 and nu2 for values in x[].    |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityF(const double &x[],const double nu1,const double nu2,double &result[])
  {
   return MathProbabilityDensityF(x,nu1,nu2,false,result);
  }
//+------------------------------------------------------------------+
//| F cumulative distribution function (CDF)                         |
//+------------------------------------------------------------------+
//| The function returns the cumulative distribution function of the |
//| F-distribution with given nu1 and nu2.                           |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the F cumulative distribution function with         |
//| parameters nu1 and nu2, evaluated at x.                          |
//+------------------------------------------------------------------+
double MathCumulativeDistributionF(const double x,const double nu1,const double nu2,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(nu1) || !MathIsValidNumber(nu2))
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
//--- check x
   if(x<=0)
      return TailLog0(tail,log_mode);
//--- calculate cdf using incomplete Beta and take into account round-off errors for probability
   double cdf=MathMin(1.0-MathBetaIncomplete(nu2/(nu2+nu1*x),nu2*0.5,nu1*0.5),1.0);
   return TailLogValue(cdf,tail,log_mode);
  }
//+------------------------------------------------------------------+
//| F cumulative distribution function (CDF)                         |
//+------------------------------------------------------------------+
//| The function returns the cumulative distribution function of the |
//| F-distribution with given nu1 and nu2.                           |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the F cumulative distribution function with         |
//| parameters nu1 and nu2, evaluated at x.                          |
//+------------------------------------------------------------------+
double MathCumulativeDistributionF(const double x,const double nu1,const double nu2,int &error_code)
  {
   return MathCumulativeDistributionF(x,nu1,nu2,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| F cumulative distribution function (CDF)                         |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the F distribution with parameters nu1 and nu2 for values in x[].|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionF(const double &x[],const double nu1,const double nu2,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu1) || !MathIsValidNumber(nu2))
      return false;
//--- check arguments
   if(nu1!=MathRound(nu1) || nu2!=MathRound(nu2) || nu1<=0 || nu2<=0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      //--- check x
      if(x_arg<=0)
         result[i]=TailLog0(tail,log_mode);
      else
        {
         //--- calculate cdf using incomplete Beta and take into account round-off errors for probability
         double cdf=MathMin(1.0-MathBetaIncomplete(nu2/(nu2+nu1*x_arg),nu2*0.5,nu1*0.5),1.0);
         result[i]=TailLogValue(cdf,tail,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| F cumulative distribution function (CDF)                         |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the F distribution with parameters nu1 and nu2 for values in x[].|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionF(const double &x[],const double nu1,const double nu2,double &result[])
  {
   return MathCumulativeDistributionF(x,nu1,nu2,true,false,result);
  }
//+------------------------------------------------------------------+
//| F-distribution quantile function (inverse CDF)                   |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of F-distribution with parameters nu1 and nu2           |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| nu1         : Numerator degrees of freedom                       |
//| nu2         : Denominator degrees of freedom                     |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of F-distribution with parameters nu1 and nu2.                   |
//+------------------------------------------------------------------+
double MathQuantileF(const double probability,const double nu1,const double nu2,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(nu1) || !MathIsValidNumber(nu2))
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

//--- calculate real probability
   double prob=TailLogProbability(probability,tail,log_mode);
//--- check probability range
   if(prob<0.0 || prob>1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check case probability==1
   if(prob==1.0)
     {
      error_code=ERR_RESULT_INFINITE;
      return QPOSINF;
     }

   error_code=ERR_OK;
   if(prob==0.0)
      return 0.0;
//--- calculate quantile using Beta distribution
   double qBeta=MathQuantileBeta(1.0-prob,nu2*0.5,nu1*0.5,error_code);
//--- return quantile;
   return (nu2/qBeta-nu2)/nu1;
  }
//+------------------------------------------------------------------+
//| F-distribution quantile function (inverse CDF)                   |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of F-distribution with parameters nu1 and nu2           |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| nu1         : Numerator degrees of freedom                       |
//| nu2         : Denominator degrees of freedom                     |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of F-distribution with parameters nu1 and nu2.                   |
//+------------------------------------------------------------------+
double MathQuantileF(const double probability,const double nu1,const double nu2,int &error_code)
  {
   return MathQuantileF(probability,nu1,nu2,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| F-distribution quantile function (inverse CDF)                   |
//+------------------------------------------------------------------+
//| The function calculates  the inverse cumulative distribution     |
//| function of the F distribution with parameters nu1 and nu2       |
//| for values from the probability[] array.                         |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| nu1         : Numerator degrees of freedom                       |
//| nu2         : Denominator degrees of freedom                     |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileF(const double &probability[],const double nu1,const double nu2,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu1) || !MathIsValidNumber(nu2))
      return false;
//--- check arguments
   if(nu1!=MathRound(nu1) || nu2!=MathRound(nu2) || nu1<=0 || nu2<=0)
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

      //--- check case probability==1,0
      if(prob==0.0)
         result[i]=0.0;
      else
      if(prob==1.0)
         result[i]=QPOSINF;
      else
        {
         //--- calculate quantile using Beta distribution
         double qBeta=MathQuantileBeta(1.0-prob,nu2*0.5,nu1*0.5,error_code);
         result[i]=(nu2/qBeta-nu2)/nu1;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| F-distribution quantile function (inverse CDF)                   |
//+------------------------------------------------------------------+
//| The function calculates  the inverse cumulative distribution     |
//| function of the F distribution with parameters nu1 and nu2       |
//| for values from probability[] array.                             |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| nu1         : Numerator degrees of freedom                       |
//| nu2         : Denominator degrees of freedom                     |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileF(const double &probability[],const double nu1,const double nu2,double &result[])
  {
   return MathQuantileF(probability,nu1,nu2,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the F-distribution                           |
//+------------------------------------------------------------------+
//| Compute the random variable from F-distribution                  |
//| with parameters nu1 and nu2.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| nu1         : Numerator degrees of freedom                       |
//| nu2         : Denominator degrees of freedom                     |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with F-distribution.                            |
//+------------------------------------------------------------------+
double MathRandomF(const double nu1,const double nu2,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(nu1) || !MathIsValidNumber(nu2))
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
//--- random F=ChiSquare(nu1)*nu2/ChiSquare(nu2)*nu1;
   double xnum = MathRandomGamma(nu1*0.5,1.0,error_code)*nu2;
   double xden = MathRandomGamma(nu2*0.5,1.0,error_code)*nu1;
//---
   double value=0.0;
   if(xden!=0)
      value= xnum/xden;
   else
     {
      error_code=ERR_NON_CONVERGENCE;
      value=QNaN;
     }
//--- return random F
   return value;
  }
//+------------------------------------------------------------------+
//| Random variate from the F distribution                           |
//+------------------------------------------------------------------+
//| Generates random variables from the F distribution with          |
//| parameters nu1 and nu2.                                          |
//|                                                                  |
//| Arguments:                                                       |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomF(const double nu1,const double nu2,const int data_count,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu1) || !MathIsValidNumber(nu2))
      return false;
//--- check arguments
   if(nu1!=MathRound(nu1) || nu2!=MathRound(nu2) || nu1<=0 || nu2<=0)
      return false;

//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      int error_code=0;
      //--- random F=ChiSquare(nu1)*nu2/ChiSquare(nu2)*nu1;
      double xnum = MathRandomGamma(nu1*0.5,1.0,error_code)*nu2;
      double xden = MathRandomGamma(nu2*0.5,1.0,error_code)*nu1;
      //---
      double value=0.0;
      if(xden!=0)
         value= xnum/xden;
      else
        {
         error_code=ERR_NON_CONVERGENCE;
         value=QNaN;
        }
      //--- random F
      result[i]=value;
     }
   return true;
  }
//+------------------------------------------------------------------+
//| F-distribution moments                                           |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of F-distribution        |
//| with parameters nu1 and nu2.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| nu1        : Numerator degrees of freedom                        |
//| nu2        : Denominator degrees of freedom                      |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsF(const double nu1,const double nu2,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(nu1) || !MathIsValidNumber(nu2))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- check arguments
   if(nu1!=MathRound(nu1) || nu1!=MathRound(nu1) || nu1<1 || nu2<1)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- calculate moments
   if(nu2>2)
      mean=nu2/(nu2-2);
   if(nu2>4)
      variance=2*nu2*nu2*(nu1+nu2-2)/(nu1*(nu2-2)*(nu2-2)*(nu2-4));
   if(nu2>6)
      skewness=2*MathSqrt(2)*MathSqrt(nu2-4)*(2*nu1+nu2-2)/(MathSqrt(nu1*(nu1+nu2-2))*(nu2-6));
   if(nu2>8)
      kurtosis=12*(nu1*(5*nu2-22)*(nu1+nu2-2)+(nu2-4)*(nu2-2)*(nu2-2))/(nu1*(nu2-8)*(nu2-6)*(nu1+nu2-2));
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
