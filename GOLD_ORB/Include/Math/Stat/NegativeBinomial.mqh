//+------------------------------------------------------------------+
//|                                             NegativeBinomial.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"
#include "Gamma.mqh"
#include "Poisson.mqh"

//+------------------------------------------------------------------+
//| Negative Binomial probability mass function (PDF)                |
//+------------------------------------------------------------------+
//| The function returns the probability mass function               |
//| of the Negative Binomial distribution with parameters r and p.   |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| r          : Number of successes                                 |
//| p          : Probability of success                              |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability mass evaluated at x.                             |
//+------------------------------------------------------------------+
double MathProbabilityDensityNegativeBinomial(const double x,const double r,const double p,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(r) || !MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(r!=MathRound(r) || r<1.0 || p<0.0 || p>1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(x<0.0)
      return TailLog0(true,log_mode);
//--- calculate gamma factor for the density
   double coef=MathRound(MathExp(MathGammaLog(r+x)-MathGammaLog(x+1.0)-MathGammaLog(r)));
//--- return density
   return TailLogValue(coef*MathPow(p,r)*MathPow(1.0-p,x),true,log_mode);
  }
//+------------------------------------------------------------------+
//| Negative Binomial probability mass function (PDF)                |
//+------------------------------------------------------------------+
//| The function returns the probability mass function               |
//| of the Negative Binomial distribution with parameters r and p.   |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| r          : Number of successes                                 |
//| p          : Probability of success                              |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability mass evaluated at x.                             |
//+------------------------------------------------------------------+
double MathProbabilityDensityNegativeBinomial(const double x,const double r,const double p,int &error_code)
  {
   return MathProbabilityDensityNegativeBinomial(x,r,p,false,error_code);
  }
//+------------------------------------------------------------------+
//| Negative Binomial probability mass function (PDF)                |
//+------------------------------------------------------------------+
//| The function calculates the probability mass function            |
//| of the Negative Binomial distribution with parameters r and p    |
//| for values from x[] array.                                       |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| r          : Number of successes                                 |
//| p          : Probability of success                              |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityNegativeBinomial(const double &x[],const double r,const double p,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(r) || !MathIsValidNumber(p))
      return false;
//--- check arguments
   if(r!=MathRound(r) || r<1.0 || p<0.0 || p>1.0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   double power_p_r=MathPow(p,r);
   double log_gamma_r=MathGammaLog(r);
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
         //--- calculate pdf
         double pdf=power_p_r*MathPow(1.0-p,x_arg)*MathRound(MathExp(MathGammaLog(r+x_arg)-MathGammaLog(x_arg+1.0)-log_gamma_r));
         result[i]=TailLogValue(pdf,true,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Negative Binomial probability mass function (PDF)                |
//+------------------------------------------------------------------+
//| The function calculates the probability mass function            |
//| of the Negative Binomial distribution with parameters r and p    |
//| for values from x[] array.                                       |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| r          : Number of successes                                 |
//| p          : Probability of success                              |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityNegativeBinomial(const double &x[],const double r,const double p,double &result[])
  {
   return MathProbabilityDensityNegativeBinomial(x,r,p,false,result);
  }
//+------------------------------------------------------------------+
//| Negative Binomial cumulative distribution function (CDF)         |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Negative Binomial distribution with parameters r and p  |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x           : The desired quantile                               |
//| r           : Number of successes                                |
//| p           : Probability of success                             |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Negative Binomial cumulative distribution       |
//| function with parameters r and p, evaluated at x.                |
//+------------------------------------------------------------------+
double MathCumulativeDistributionNegativeBinomial(const double x,const double r,double p,const bool tail,const bool log_mode,int error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(r) || !MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(r!=MathRound(r) || r<1.0 || p<0.0 || p>1.0 || x<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(x<0.0)
      return TailLog0(tail,log_mode);
   int err_code=0;
//--- calculate max term of the sum
   int max_j=(int)MathFloor(x);
   double p1=1.0-p;
//--- initial factors
   double factor1=MathFactorial((int)r-1);
   double factor2=1.0;
   double factor_p=1.0;
   double factor_r=1.0/factor1;
   double power_p_r=MathPowInt(p,int(r))*factor_r;
   double cdf=0.0;
   for(int j=0; j<=max_j; j++)
     {
      if(j>0)
        {
         factor1*=(j+1);
         factor2*=j;
         factor_p*=p1;
        }
      double pdf=power_p_r*factor1*factor_p/factor2;
      cdf+=pdf;
     }
//--- take into account round-off errors for probability
   return TailLogValue(MathMin(cdf,1.0),tail,log_mode);
  }
//+------------------------------------------------------------------+
//| Negative Binomial cumulative distribution function (CDF)         |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Negative Binomial distribution with parameters r and p  |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x           : The desired quantile                               |
//| r           : Number of successes                                |
//| p           : Probability of success                             |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Negative Binomial cumulative distribution       |
//| function with parameters r and p, evaluated at x.                |
//+------------------------------------------------------------------+
double MathCumulativeDistributionNegativeBinomial(const double x,const double r,double p,int error_code)
  {
   return MathCumulativeDistributionNegativeBinomial(x,r,p,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Negative Binomial cumulative distribution function (CDF)         |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function     |
//| of the Negative Binomial distribution with parameters r and p    |
//| for values from x[] array.                                       |
//|                                                                  |
//| Arguments:                                                       |
//| x        : Array with random variables                           |
//| r        : Number of successes                                   |
//| p        : Probability of success                                |
//| tail     : Flag to calculate lower tail                          |
//| log_mode : Logarithm mode, if true it calculates Log values      |
//| result   : Array with calculated values                          |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionNegativeBinomial(const double &x[],const double r,double p,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(r) || !MathIsValidNumber(p))
      return false;
//--- check arguments
   if(r!=MathRound(r) || r<1.0 || p<0.0 || p>1.0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;
//--- common factors 
   double fact1=MathFactorial((int)r-1);
   double factor_r=1.0/fact1;
   double power_p_r=MathPowInt(p,int(r))*factor_r;
   double p1=1.0-p;
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
         int err_code=0;
         //--- calculate max term of the sum
         int max_j=(int)MathFloor(x_arg);
         //--- initial factors
         double factor1=fact1;
         double factor2=1.0;
         double factor_p=1.0;
         double cdf=0.0;
         for(int j=0; j<=max_j; j++)
           {
            if(j>0)
              {
               factor1*=(j+1);
               factor2*=j;
               factor_p*=p1;
              }
            double pdf=power_p_r*factor1*factor_p/factor2;
            cdf+=pdf;
           }
         //--- take into account round-off errors for probability
         result[i]=TailLogValue(MathMin(cdf,1.0),tail,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Negative Binomial cumulative distribution function (CDF)         |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function     |
//| of the Negative Binomial distribution with parameters r and p    |
//| for values from x[] array.                                       |
//|                                                                  |
//| Arguments:                                                       |
//| x      : Array with random variables                             |
//| r      : Number of successes                                     |
//| p      : Probability of success                                  |
//| result : Array with calculated values                            |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionNegativeBinomial(const double &x[],const double r,double p,double &result[])
  {
   return MathCumulativeDistributionNegativeBinomial(x,r,p,true,false,result);
  }
//+------------------------------------------------------------------+
//| Negative Binomial distribution quantile function (inverse CDF)   |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Negative Binomial distribution with parameters   |
//| r and p for the desired probability.                             |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| r           : Number of successes                                |
//| p           : Probability of success                             |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Negative Binomial distribution with parameters r and p.   |
//+------------------------------------------------------------------+
double MathQuantileNegativeBinomial(const double probability,const double r,const double p,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(probability) || !MathIsValidNumber(r) || !MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(r!=MathRound(r) || r<1.0 || p<0.0 || p>1.0)
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
//--- check cases p=0 and p=1
   if(prob==1.0)
     {
      error_code=ERR_RESULT_INFINITE;
      return QPOSINF;
     }
   error_code=ERR_OK;
   if(prob==0.0)
      return 0.0;

   int max_terms=1000;
   int err_code=0;
//--- factors 
   double fact1=MathFactorial((int)r-1);
   double factor_r=1.0/fact1;
   double power_p_r=MathPowInt(p,int(r))*factor_r;
   double p1=1.0-p;
//--- initial factors
   double factor1=fact1;
   double factor2=1.0;
   double factor_p=1.0;
   double cdf=0.0;
   int j=0;
   while(cdf<prob && j<max_terms)
     {
      if(j>0)
        {
         factor1*=(j+1);
         factor2*=j;
         factor_p*=p1;
        }
      double pdf=power_p_r*factor1*factor_p/factor2;
      cdf+=pdf;
      j++;
     }
//--- check convergence
   if(j<max_terms)
     {
      if(j==0)
         return 0;
      else
         return j-1;
     }
   else
     {
      error_code=ERR_NON_CONVERGENCE;
      return 0;
     }
  }
//+------------------------------------------------------------------+
//| Negative Binomial distribution quantile function (inverse CDF)   |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Negative Binomial distribution with parameters   |
//| r and p for the desired probability.                             |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| r           : Number of successes                                |
//| p           : Probability of success                             |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Negative Binomial distribution with parameters r and p.   |
//+------------------------------------------------------------------+
double MathQuantileNegativeBinomial(const double probability,const double r,const double p,int &error_code)
  {
   return MathQuantileNegativeBinomial(probability,r,p,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Negative Binomial distribution quantile function (inverse CDF)   |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Negative Binomial distribution with parameters   |
//| r and p for values form the probability[] array.                 |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| r           : Number of successes                                |
//| p           : Probability of success                             |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileNegativeBinomial(const double &probability[],const double r,const double p,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(r) || !MathIsValidNumber(p))
      return false;
//--- check arguments
   if(r!=MathRound(r) || r<1.0 || p<0.0 || p>1.0)
      return false;

   int data_count=ArraySize(probability);
   if(data_count==0)
      return false;
//--- common factors
   double fact1=MathFactorial((int)r-1);
   double factor_r=1.0/fact1;
   double power_p_r=MathPowInt(p,int(r))*factor_r;
   double p1=1.0-p;
   int max_terms=500;
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
         double factor1=fact1;
         double factor2=1.0;
         double factor_p=1.0;
         double cdf=0.0;
         int j=0;
         while(cdf<prob && j<max_terms)
           {
            if(j>0)
              {
               factor1*=(j+1);
               factor2*=j;
               factor_p*=p1;
              }
            double pdf=power_p_r*factor1*factor_p/factor2;
            cdf+=pdf;
            j++;
           }
         if(j<max_terms)
           {
            if(j==0)
               result[i]=0;
            else
               result[i]=j-1;
           }
         else
            return false;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Negative Binomial distribution quantile function (inverse CDF)   |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Negative Binomial distribution with parameters   |
//| r and p for values from the probability[] array.                 |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| r           : Number of successes                                |
//| p           : Probability of success                             |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileNegativeBinomial(const double &probability[],const double r,const double p,double &result[])
  {
   return MathQuantileNegativeBinomial(probability,r,p,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Negative Binomial distribution           |
//+------------------------------------------------------------------+
//| Computes the random variable from the Negative Binomial          |
//| distribution with parameters r and p.                            |
//|                                                                  |
//| Arguments:                                                       |
//| r           : Number of successes                                |
//| p           : Probability of success                             |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Negative Binomial distribution.            |
//+------------------------------------------------------------------+
double MathRandomNegativeBinomial(const double r,const double p,int error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(r) || !MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(r<=0.0 || p<=0.0 || p>=1.0)
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
   double r_gamma=MathRandomGamma(r,(1-p)/p);
   return MathRandomPoisson(r_gamma,error_code);
  }
//+------------------------------------------------------------------+
//| Random variate from the Negative Binomial distribution            |
//+------------------------------------------------------------------+
//| Generates random variables from the Negative Binomial            |
//| distribution with parameters r and p.                            |
//|                                                                  |
//| Arguments:                                                       |
//| r          : Number of successes                                 |
//| p          : Probability of success                              |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomNegativeBinomial(const double r,const double p,const int data_count,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(r) || !MathIsValidNumber(p))
      return false;
//--- check arguments
   if(r<=0.0 || p<=0.0 || p>=1.0)
      return false;

   double p_coef=(1-p)/p;
   int error_code=0;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double r_gamma=MathRandomGamma(r,p_coef);
      result[i]=MathRandomPoisson(r_gamma,error_code);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Negative Binomial distribution moments                           |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of Negative Binomial     |
//| distribution with parameters r and p.                            |
//|                                                                  |
//| Arguments:                                                       |
//| r          : Number of successes                                 |
//| p          : Probability of success                              |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsNegativeBinomial(const double r,double p,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(r) || !MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- check arguments
   if(r!=MathRound(r) || r<1.0 || p<=0.0 || p>=1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- calculate moments
   mean    =r*(1.0-p)/p;
   variance=mean/p;
   skewness=(2.0-p)/MathSqrt((r*(1.0-p)));
   kurtosis=(p*p-6*p+6)/(r*(1.0-p));
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
