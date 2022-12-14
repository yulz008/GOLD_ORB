//+------------------------------------------------------------------+
//|                                               Hypergeometric.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"

//+------------------------------------------------------------------+
//| Hypergeometric probability mass function (PDF)                   |
//+------------------------------------------------------------------+
//| The function returns the probability mass function               |
//| of the Hypergeometric distribution with parameters m,n,k.        |
//|                f(x,m,k,n)=C(k,x)*C(m-k,n-x)/C(m,n)               |
//| where binomial coefficient C(n,k)=n!/(k!*(n-k)!                  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired number of objects                       |
//| m          : Size of the population                              |
//| k          : Number of items with the desired characteristic     |
//|              in the population                                   |
//| n          : Number of samples drawn                             |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability mass function, evaluated at x.                   |
//+------------------------------------------------------------------+
double MathProbabilityDensityHypergeometric(const double x,const double m,const double k,const double n,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(m) || !MathIsValidNumber(k) || !MathIsValidNumber(n))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return(QNaN);
     }
//--- m,k,n must be integer
   if(m!=MathRound(m) || k!=MathRound(k) || n!=MathRound(n))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return(QNaN);
     }
//--- m,k,n must be positive
   if(m<0 || k<0 || n<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return(QNaN);
     }
//--- check ranges
   if(n>m || k>m)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return(QNaN);
     }

   error_code=ERR_OK;
//--- check ranges
   if(x>n)
      return TailLog0(true,log_mode);
   if(x>k || m-k-n+x+1<=0)
      return TailLog0(true,log_mode);
//--- calculate log binomial coefficients
   double log_pdf=MathBinomialCoefficientLog(k,x)+MathBinomialCoefficientLog(m-k,n-x)-MathBinomialCoefficientLog(m,n);
   if(log_mode==true)
      return log_pdf;
//--- return hypergeometric density
   return MathExp(log_pdf);
  }
//+------------------------------------------------------------------+
//| Hypergeometric probability mass function (PDF)                   |
//+------------------------------------------------------------------+
//| The function returns the probability mass function               |
//| of the Hypergeometric distribution with parameters m,n,k.        |
//|                f(x,m,k,n)=C(k,x)*C(m-k,n-x)/C(m,n)               |
//| where binomial coefficient C(n,k)=n!/(k!*(n-k)!                  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired number of objects                       |
//| m          : Size of the population                              |
//| k          : Number of items with the desired characteristic     |
//|              in the population                                   |
//| n          : Number of samples drawn                             |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability mass function, evaluated at x.                   |
//+------------------------------------------------------------------+
double MathProbabilityDensityHypergeometric(const double x,const double m,const double k,const double n,int &error_code)
  {
   return MathProbabilityDensityHypergeometric(x,m,k,n,false,error_code);
  }
//+------------------------------------------------------------------+
//| Hypergeometric probability mass function (PDF)                   |
//+------------------------------------------------------------------+
//| The function calculates the probability mass function of the     |
//| Hypergeometric distribution with parameter m,k,n for values in x.|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| x          : The desired number of objects                       |
//| m          : Size of the population                              |
//| k          : Number of items with the desired characteristic     |
//|              in the population                                   |
//| n          : Number of samples drawn                             |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityHypergeometric(const double &x[],const double m,const double k,const double n,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(m) || !MathIsValidNumber(k) || !MathIsValidNumber(n))
      return false;
//--- m,k,n must be integer
   if(m!=MathRound(m) || k!=MathRound(k) || n!=MathRound(n))
      return false;
//--- m,k,n must be positive
   if(m<0 || k<0 || n<0)
      return false;
//--- check ranges
   if(n>m || k>m)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   double m_k=m-k;
   int error_code=0;
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(!MathIsValidNumber(x_arg))
         return false;

      if(x_arg<0 || x_arg!=MathRound(x_arg))
         return false;

      if(x_arg>n)
         result[i]=TailLog0(true,log_mode);
      else
      //--- check ranges
      if(x_arg>k || m_k-n+x_arg+1<=0)
         result[i]=TailLog0(true,log_mode);
      else
        {
         //--- calculate log binomial coefficients
         double log_pdf=MathBinomialCoefficientLog(k,x_arg)+MathBinomialCoefficientLog(m_k,n-x_arg)-MathBinomialCoefficientLog(m,n);
         if(log_mode==true)
            result[i]=log_pdf;
         else
            result[i]=MathExp(log_pdf);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Hypergeometric probability mass function (PDF)                   |
//+------------------------------------------------------------------+
//| The function calculates the probability mass function of the     |
//| Hypergeometric distribution with parameter m,k,n for values in x.|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| x          : The desired number of objects                       |
//| m          : Size of the population                              |
//| k          : Number of items with the desired characteristic     |
//|              in the population                                   |
//| n          : Number of samples drawn                             |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityHypergeometric(const double &x[],const double m,const double k,const double n,double &result[])
  {
   return MathProbabilityDensityHypergeometric(x,m,k,n,false,result);
  }
//+------------------------------------------------------------------+
//| Hypergeometric cumulative distribution function (CDF)            |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Hypergeometric distribution with parameters m,n,k       |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired number of objects                       |
//| m          : Size of the population                              |
//| k          : Number of items with the desired characteristic     |
//|              in the population                                   |
//| n          : Number of samples drawn                             |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode,if true it calculates Log values     |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Hypergeometric cumulative distribution function |
//| with parameters m,n,k, evaluated at x.                           |
//+------------------------------------------------------------------+
//| Based on algorithm by John Burkardt                              |
//+------------------------------------------------------------------+
double MathCumulativeDistributionHypergeometric(const double x,const double m,const double k,const double n,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(m) || !MathIsValidNumber(k) || !MathIsValidNumber(n))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- m,k,n,x must be integer
   if(m!=MathRound(m) || k!=MathRound(k) || n!=MathRound(n) || x!=MathRound(x))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- m,k,n,x must be positive
   if(m<0 || k<0 || n<0 || x<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check ranges
   if(n>m || k>m)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(x>=n || x>=k)
      return TailLog1(tail,log_mode);
//--- calculate cdf
   double pdf = MathExp(MathBinomialCoefficientLog(m-k,n)-MathBinomialCoefficientLog(m,n));
   double cdf = pdf;
   double coef=m-k-n+1;
   for(int j=0; j<=x-1; j++)
     {
      pdf = pdf*(k-j)*(n-j)/((j+1)*(coef+j));
      cdf = cdf + pdf;
     }
   return TailLogValue(MathMin(cdf,1.0),tail,log_mode);
  }
//+------------------------------------------------------------------+
//| Hypergeometric cumulative distribution function (CDF)            |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Hypergeometric distribution with parameters m,n,k       |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired number of objects                       |
//| m          : Size of the population                              |
//| k          : Number of items with the desired characteristic     |
//|              in the population                                   |
//| n          : Number of samples drawn                             |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Hypergeometric cumulative distribution function |
//| with parameters m,n,k, evaluated at x.                           |
//+------------------------------------------------------------------+
double MathCumulativeDistributionHypergeometric(const double x,const double m,const double k,const double n,int &error_code)
  {
   return MathCumulativeDistributionHypergeometric(x,m,k,n,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Hypergeometric cumulative distribution function (CDF)            |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Hypergeometric distribution with parameters m,k,n for        |
//| the values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| m          : Size of the population                              |
//| k          : Number of items with the desired characteristic     |
//|              in the population                                   |
//| n          : Number of samples drawn                             |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode,if true it calculates Log values     |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
//| Based on algorithm by John Burkardt                              |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionHypergeometric(const double &x[],const double m,const double k,const double n,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(m) || !MathIsValidNumber(k) || !MathIsValidNumber(n))
      return false;
//--- m,k,n,x must be integer
   if(m!=MathRound(m) || k!=MathRound(k) || n!=MathRound(n))
      return false;
//--- m,k,n must be positive
   if(m<0 || k<0 || n<0)
      return false;
//--- check ranges
   if(n>m || k>m)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   double coef=m-k-n+1;

   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(!MathIsValidNumber(x_arg))
         return false;

      //--- x must be positive and integer
      if(x_arg<0 || x_arg!=MathRound(x_arg))
         return false;

      //--- check ranges
      if(x_arg>=n || x_arg>=k)
         result[i]=TailLog1(tail,log_mode);
      else
        {
         //--- calculate cdf
         double pdf = MathExp(MathBinomialCoefficientLog(m-k,n)-MathBinomialCoefficientLog(m,n));
         double cdf = pdf;
         for(int j=0; j<=x_arg-1; j++)
           {
            pdf = pdf*(k-j)*(n-j)/((j+1)*(coef+j));
            cdf = cdf + pdf;
           }
         result[i]=TailLogValue(MathMin(cdf,1.0),tail,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Hypergeometric cumulative distribution function (CDF)            |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Hypergeometric distribution with parameters m,k,n for        |
//| the values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| m          : Size of the population                              |
//| k          : Number of items with the desired characteristic     |
//|              in the population                                   |
//| n          : Number of samples drawn                             |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionHypergeometric(const double &x[],const double m,const double k,const double n,double &result[])
  {
   return MathCumulativeDistributionHypergeometric(x,m,k,n,true,false,result);
  }
//+------------------------------------------------------------------+
//| Hypergeometric distribution quantile function (inverse CDF)      |
//+------------------------------------------------------------------+
//| Computes the inverse cumulative distribution function of the     |
//| Hypergeometric distribution with parameters m,n,k for the        |
//| desired probability.                                             |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The probability                                    |
//| m           : Size of the population                             |
//| k           : Number of items with the desired characteristic    |
//|               in the population                                  |
//| n           : Number of samples drawn                            |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode,if true it calculates for Log values|
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The smallest value x, such that the hypergeometric CDF(x)        |
//| equals or exceeds the desired probability.                       |
//+------------------------------------------------------------------+
//| Based on algorithm by John Burkardt                              |
//+------------------------------------------------------------------+
double MathQuantileHypergeometric(const double probability,const double m,const double k,const double n,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(probability) || !MathIsValidNumber(m) || !MathIsValidNumber(k) || !MathIsValidNumber(n))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- m,k,n,x must be integer
   if(m!=MathRound(m) || k!=MathRound(k) || n!=MathRound(n))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- m,k,n must be positive
   if(m<0 || k<0 || n<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check ranges
   if(n>m || k>m)
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
//--- check probability
   if(prob==0)
      return 0.0;
   if(prob==1.0)
      return QPOSINF;

   int max_terms=1000;
   prob*=1-1000*DBL_EPSILON;
   double m_k=m-k;
   double pdf = MathExp(MathBinomialCoefficientLog(m_k,n)-MathBinomialCoefficientLog(m,n));
   double cdf = pdf;
   double coef=m_k-n+1;
   int j=0;
   while(cdf<prob && j<max_terms)
     {
      pdf = pdf*(k-j)*(n-j)/((j+1)*(coef+j));
      cdf = cdf + pdf;
      j++;
     }
   return j;
  }
//+------------------------------------------------------------------+
//| Hypergeometric distribution quantile function (inverse CDF)      |
//+------------------------------------------------------------------+
//| Computes the inverse cumulative distribution function of the     |
//| Hypergeometric distribution with parameters m,n,k for the        |
//| desired probability.                                             |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The probability                                    |
//| m           : Size of the population                             |
//| k           : Number of items with the desired characteristic    |
//|               in the population                                  |
//| n           : Number of samples drawn                            |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The smallest value x, such that the hypergeometric CDF(x)        |
//| equals or exceeds the desired probability.                       |
//+------------------------------------------------------------------+
double MathQuantileHypergeometric(const double probability,const double m,const double k,const double n,int &error_code)
  {
   return MathQuantileHypergeometric(probability,m,k,n,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Hypergeometric distribution quantile function (inverse CDF)      |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Hypergeometric distribution with parameters      |
//| m,k,n for values from the probability[] array.                   |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| m           : Size of the population                             |
//| k           : Number of items with the desired characteristic    |
//|               in the population                                  |
//| n           : Number of samples drawn                            |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode,if true it calculates for Log values|
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
//| Based on algorithm by John Burkardt                              |
//+------------------------------------------------------------------+
bool MathQuantileHypergeometric(const double &probability[],const double m,const double k,const double n,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(m) || !MathIsValidNumber(k) || !MathIsValidNumber(n))
      return false;
//--- m,k,n,x must be integer
   if(m!=MathRound(m) || k!=MathRound(k) || n!=MathRound(n))
      return false;
//--- m,k,n must be positive
   if(m<0 || k<0 || n<0)
      return false;
//--- check ranges
   if(n>m || k>m)
      return false;

   int data_count=ArraySize(probability);
   if(data_count==0)
      return false;

   int max_terms=1000;
   double m_k=m-k;
   double pdf0= MathExp(MathBinomialCoefficientLog(m_k,n)-MathBinomialCoefficientLog(m,n));
   double coef=m_k-n+1;

   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      //--- calculate real probability
      double prob=TailLogProbability(probability[i],tail,log_mode);
      //--- check probability range
      if(prob<0.0 || prob>1.0)
         return false;

      //--- check probability
      if(prob==0.0)
         result[i]=0.0;
      else
      if(prob==1.0)
         result[i]=QPOSINF;
      else
        {
         prob*=1-1000*DBL_EPSILON;
         double pdf = pdf0;
         double cdf = pdf;
         int j=0;
         while(cdf<prob && j<max_terms)
           {
            pdf = pdf*(k-j)*(n-j)/((j+1)*(coef+j));
            cdf = cdf + pdf;
            j++;
           }
         result[i]=j;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Hypergeometric distribution quantile function (inverse CDF)      |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Hypergeometric distribution with parameters      |
//| m,k,n for values from the probability[] array.                   |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| m           : Size of the population                             |
//| k           : Number of items with the desired characteristic    |
//|               in the population                                  |
//| n           : Number of samples drawn                            |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileHypergeometric(const double &probability[],const double m,const double k,const double n,double &result[])
  {
   return MathQuantileHypergeometric(probability,m,k,n,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Hypergeometric distribution             |
//+------------------------------------------------------------------+
//| Compute the random variable from the Hypergeometric distribution |
//| with parameters m,n,k.                                           |
//|                                                                  |
//| Arguments:                                                       |
//| m           : Size of the population                             |
//| k           : Number of items with the desired characteristic    |
//|               in the population                                  |
//| n           : Number of samples drawn                            |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Hypergeometric distribution.               |
//+------------------------------------------------------------------+
//| Author:  John Burkardt                                           |
//|                                                                  |
//| Reference:                                                       |
//| Jerry Banks, editor, Handbook of Simulation,                     |
//| Engineering and Management Press Books, 1998, page 165.          |
//+------------------------------------------------------------------+
double MathRandomHypergeometric(const double m,const double k,const double n,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(m) || !MathIsValidNumber(k) || !MathIsValidNumber(n))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- m,k,n,x must be integer
   if(m!=MathRound(m) || k!=MathRound(k) || n!=MathRound(n))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- m,k,n must be positive
   if(m<0 || k<0 || n<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check ranges
   if(n>m || k>m)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- generate random number
   double prob=MathRandomNonZero();
   prob*=1-1000*DBL_EPSILON;
   int max_terms=1000;
   double m_k=m-k;
   double coef=m_k-n+1;
   double pdf= MathExp(MathBinomialCoefficientLog(m_k,n)-MathBinomialCoefficientLog(m,n));
   double cdf= pdf;
   int j=0;
   while(cdf<prob && j<max_terms)
     {
      pdf = pdf*(k-j)*(n-j)/((j+1)*(coef+j));
      cdf = cdf + pdf;
      j++;
     }
   return j;
  }
//+------------------------------------------------------------------+
//| Random variate from the Hypergeometric distribution              |
//+------------------------------------------------------------------+
//| Generates random variables from the Hypergeometric distribution  |
//| with parameters m,k,n.                                           |
//|                                                                  |
//| Arguments:                                                       |
//| m          : Size of the population                              |
//| k          : Number of items with the desired characteristic     |
//|              in the population                                   |
//| n          : Number of samples drawn                             |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomHypergeometric(const double m,const double k,const double n,const int data_count,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(m) || !MathIsValidNumber(k) || !MathIsValidNumber(n))
      return false;
//--- m,k,n,x must be integer
   if(m!=MathRound(m) || k!=MathRound(k) || n!=MathRound(n))
      return false;
//--- m,k,n must be positive
   if(m<0 || k<0 || n<0)
      return false;
//--- check ranges
   if(n>m || k>m)
      return false;
//--- prepare coefficients
   int max_terms=1000;
   double m_k=m-k;
   double coef=m_k-n+1;
   double pdf0= MathExp(MathBinomialCoefficientLog(m_k,n)-MathBinomialCoefficientLog(m,n));
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      //--- generate random number
      double prob=MathRandomNonZero();
      prob*=1-1000*DBL_EPSILON;
      //--- calculate using quantile
      double pdf = pdf0;
      double cdf = pdf;
      int j=0;
      while(cdf<prob && j<max_terms)
        {
         pdf = pdf*(k-j)*(n-j)/((j+1)*(coef+j));
         cdf = cdf + pdf;
         j++;
        }
      result[i]=j;
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Hypergeometric distribution moments                              |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of the Hypergeometric    |
//| distribution with parameters m,n,k.                              |
//|                                                                  |
//| Arguments:                                                       |
//| m          : Size of the population                              |
//| k          : Number of items with the desired characteristic     |
//|              in the population                                   |
//| n          : Number of samples drawn                             |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsHypergeometric(const double m,const double k,const double n,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(m) || !MathIsValidNumber(k) || !MathIsValidNumber(n))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- m,k,n must be integer
   if(m!=MathRound(m) || k!=MathRound(k) || n!=MathRound(n))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }
//--- m,k,n must be positive
   if(m<0 || k<0 || n<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }
//--- check ranges
   if(n>m || k>m)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- calculate moments
   mean    =n*k/m;
   variance=k*n*(1-k/m)*(m-n)/(m*(m-1));
   skewness=MathSqrt(m-1)*(m-2*k)*(m-2*n)/((m-2)*MathSqrt(k*n*(m-k)*(m-n)));
   kurtosis=(m-1)*m*m/(k*n*(m-3)*(m-2)*(m-k)*(m-n));
   kurtosis*=3*k*(m-k)*(m*m*(n-2)-m*n*n+6*n*(m-n))/(m*m)-6*n*(m-n)+m*(m+1);
   kurtosis-=3;
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
