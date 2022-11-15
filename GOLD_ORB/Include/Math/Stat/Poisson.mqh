//+------------------------------------------------------------------+
//|                                                      Poisson.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"
#include "Gamma.mqh"

//+------------------------------------------------------------------+
//| Poisson probability mass function (PDF)                          |
//+------------------------------------------------------------------+
//| The function returns the probability mass function               |
//| of the Poisson distribution with parameter lambda.               |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| lambda     : Mean                                                |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability mass evaluated at x.                             |
//+------------------------------------------------------------------+
double MathProbabilityDensityPoisson(const double x,const double lambda,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(lambda))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- lambda must be positive, x must be integer
   if(lambda<=0.0 || x!=MathRound(x))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- check x
   if(x<0.0)
      return TailLog0(true,log_mode);

//--- calculate log pdf using LogGamma
   double log_pdf=-lambda+x*MathLog(lambda)-MathGammaLog(x+1.0);
   if(log_mode)
      return log_pdf;
//--- return density
   return MathExp(log_pdf);
  }
//+------------------------------------------------------------------+
//| Poisson probability mass function (PDF)                          |
//+------------------------------------------------------------------+
//| The function returns the probability mass function               |
//| of the Poisson distribution with parameter lambda.               |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| lambda     : Mean                                                |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability mass evaluated at x.                             |
//+------------------------------------------------------------------+
double MathProbabilityDensityPoisson(const double x,const double lambda,int &error_code)
  {
   return MathProbabilityDensityPoisson(x,lambda,false,error_code);
  }
//+------------------------------------------------------------------+
//| Poisson probability mass function (PDF)                          |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Poisson distribution with parameter lambda for values in x[].|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| lambda     : Mean                                                |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityPoisson(const double &x[],const double lambda,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(lambda))
      return false;
//--- lambda must be positive
   if(lambda<=0.0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(!MathIsValidNumber(x_arg))
         return false;

      if(x_arg!=MathRound(x_arg))
         return false;

      //--- check x
      if(x_arg<0.0)
         result[i]=TailLog0(true,log_mode);
      else
        {
         //--- calculate log pdf using LogGamma
         double log_pdf=-lambda+x_arg*MathLog(lambda)-MathGammaLog(x_arg+1.0);
         if(log_mode)
            result[i]=log_pdf;
         else
            result[i]=MathExp(log_pdf);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Poisson probability mass function (PDF)                          |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Poisson distribution with parameter lambda for values in x[].|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| lambda     : Mean                                                |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityPoisson(const double &x[],const double lambda,double &result[])
  {
   return MathProbabilityDensityPoisson(x,lambda,false,result);
  }
//+------------------------------------------------------------------+
//| Poisson cumulative distribution function (CDF)                   |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from Poisson distribution with parameter lambda                  |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| lambda     : Mean                                                |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Poisson cumulative distribution function with   |
//| parameter lambda, evaluated at x.                                |
//+------------------------------------------------------------------+
double MathCumulativeDistributionPoisson(const double x,const double lambda,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(lambda))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- lambda must be positive, x must be integer
   if(lambda<=0.0 || x!=MathRound(x))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- check x
   if(x<0.0)
      return TailLog0(tail,log_mode);
   int err_code=0;

   int t=(int)MathFloor(x+10e-10);
   double cdf=MathCumulativeDistributionGamma(lambda,t+1,1,false,false,err_code);
   return TailLogValue(cdf,tail,log_mode);
  }
//+------------------------------------------------------------------+
//| Poisson cumulative distribution function (CDF)                   |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from Poisson distribution with parameter lambda                  |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| lambda     : Mean                                                |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Poisson cumulative distribution function with   |
//| parameter lambda, evaluated at x.                                |
//+------------------------------------------------------------------+
double MathCumulativeDistributionPoisson(const double x,const double lambda,int &error_code)
  {
   return MathCumulativeDistributionPoisson(x,lambda,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Poisson cumulative distribution function (CDF)                   |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Poisson distribution with parameter lambda for values in x[].|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| lambda     : Mean                                                |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionPoisson(const double &x[],const double lambda,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(lambda))
      return false;
//--- lambda must be positive
   if(lambda<=0.0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);
   double coef_lambda=MathExp(-lambda);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(x_arg!=MathRound(x_arg))
         return false;

      if(x_arg<0.0)
         result[i]=TailLog0(tail,log_mode);
      else
        {
         int err_code=0;
         int t=(int)MathFloor(x_arg+10e-10);
         double cdf=MathCumulativeDistributionGamma(lambda,t+1,1,false,false,err_code);
         result[i]=TailLogValue(cdf,tail,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Poisson cumulative distribution function (CDF)                   |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Poisson distribution with parameter lambda for values in x[].|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| lambda     : Mean                                                |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionPoisson(const double &x[],const double lambda,double &result[])
  {
   return MathCumulativeDistributionPoisson(x,lambda,true,false,result);
  }
//+------------------------------------------------------------------+
//| Poisson distribution quantile function (inverse CDF)             |
//+------------------------------------------------------------------+
//| Computes the inverse cumulative distribution function of the     |
//| Poisson distribution with parameter lambda for the desired       |
//| probability.                                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| lambda      : Mean                                               |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode,if true it calculates for Log values|
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Poisson distribution with parameter lambda.               |
//+------------------------------------------------------------------+
double MathQuantilePoisson(const double probability,const double lambda,const bool tail,const bool log_mode,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(probability) || !MathIsValidNumber(lambda))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- lambda must be positive
   if(lambda<=0.0)
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
//--- check
   if(prob==1.0)
     {
      error_code=ERR_RESULT_INFINITE;
      return QPOSINF;
     }
   error_code=ERR_OK;
   if(prob==0.0)
      return 0.0;

   prob*=1-1000*DBL_EPSILON;
   int    err_code=0;
   int    j=0;
   const int max_terms=500;
   double coef_lambda=MathExp(-lambda);
   double pwr_lambda=1.0;
   double inverse_fact=1.0;
   double sum=0;
//--- direct calculation of the quantile
   while(sum<prob && j<max_terms)
     {
      if(j>0)
        {
         pwr_lambda*=lambda;
         inverse_fact/=j;
        }
      sum+=coef_lambda*pwr_lambda*inverse_fact;
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
      error_code=ERR_RESULT_INFINITE;
      return QPOSINF;
     }
  }
//+------------------------------------------------------------------+
//| Poisson distribution quantile function (inverse CDF)             |
//+------------------------------------------------------------------+
//| Computes the inverse cumulative distribution function of the     |
//| Poisson distribution with parameter lambda for the desired       |
//| probability.                                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| lambda      : Mean                                               |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of Poisson distribution with parameter lambda.                   |
//+------------------------------------------------------------------+
double MathQuantilePoisson(const double probability,const double lambda,int &error_code)
  {
   return MathQuantilePoisson(probability,lambda,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Poisson distribution quantile function (inverse CDF)             |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Poisson distribution with parameter lambda       |
//| for values from the probability[] array.                         |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| lambda      : Mean                                               |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantilePoisson(const double &probability[],const double lambda,const bool tail,const bool log_mode,double &result[])
  {
//--- NaN
   if(!MathIsValidNumber(lambda))
      return false;
//--- lambda must be positive
   if(lambda<=0.0)
      return false;

   int data_count=ArraySize(probability);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);
   double coef_lambda=MathExp(-lambda);

   for(int i=0; i<data_count; i++)
     {
      //--- calculate real probability
      double prob=TailLogProbability(probability[i],tail,log_mode);

      //--- check probability range
      if(prob<0.0 || prob>1.0)
         return false;
      else
      if(prob==1.0)
         result[i]=QPOSINF;
      if(prob==0.0)
         result[i]=0;
      else
        {
         prob*=1-1000*DBL_EPSILON;
         int    err_code=0;
         int    j=0;
         double sum=0.0;
         const int max_terms=500;
         double pwr_lambda=1.0;
         double inverse_fact=1.0;
         //--- direct calculation
         while(sum<prob && j<max_terms)
           {
            if(j>0)
              {
               pwr_lambda*=lambda;
               inverse_fact/=j;
              }
            sum+=coef_lambda*pwr_lambda*inverse_fact;
            j++;
           }
         //--- check convergence
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
//| Poisson distribution quantile function (inverse CDF)             |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the Poisson distribution with parameter lambda       |
//| for values from the probability[] array.                         |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| lambda      : Mean                                               |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantilePoisson(const double &probability[],const double lambda,double &result[])
  {
   return MathQuantilePoisson(probability,lambda,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Poisson distribution                     |
//+------------------------------------------------------------------+
//| Compute the random variable from the Poisson distribution        |
//| with parameter lambda.                                           |
//|                                                                  |
//| Arguments                                                        |
//| lambda      : Mean                                               |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Poisson distribution.                      |
//+------------------------------------------------------------------+
//| Original FORTRAN77 version by Barry Brown, James Lovato.         |
//| C version by John Burkardt.                                      |
//|                                                                  |
//| Reference:                                                       |
//| Joachim Ahrens, Ulrich Dieter, "Computer Generation of Poisson   |
//| "Deviates From Modified Normal Distributions",                   |
//| ACM Transactions on Mathematical Software,                       |
//| Volume 8, Number 2, June 1982, pages 163-179.                    |
//+------------------------------------------------------------------+
double MathRandomPoisson(const double lambda)
  {
   const double a0 = -0.5;
   const double a1 =  0.3333333;
   const double a2 = -0.2500068;
   const double a3 =  0.2000118;
   const double a4 = -0.1661269;
   const double a5 =  0.1421878;
   const double a6 = -0.1384794;
   const double a7 =  0.1250060;
   int    kflag;
   double fk=0,difmuk=0;
   double e=0,fx,fy,g,p0,px,py,p,q,s,t,u=0,v,x,xx;
   int    value=0;
//--- start new table and calculate P0
   if(lambda<10.0)
     {
      int m=MathMax(1,(int)(lambda));
      p = MathExp(-lambda);
      q = p;
      p0= p;
      //---  uniform sample for inversion method
      for(;;)
        {
         u=MathRandomNonZero();
         value=0;

         if(u<=p0)
            return value;
         //--- creation of new Poisson probabilities
         for(int k=1; k<=35; k++)
           {
            p=p*lambda/double(k);
            q=q+p;
            if(u<=q)
              {
               value=k;
               return value;
              }
           }
        }
     }
   else
     {
      s=MathSqrt(lambda);
      double d=6.0*lambda*lambda;
      int l=(int)(lambda-1.1484);
      //--- generate normal deviate
      double f,x1,x2,r2;
      do
        {
         x1=2.0*MathRandomNonZero()-1.0;
         x2=2.0*MathRandomNonZero()-1.0;
         r2=x1*x1+x2*x2;
        }
      while(r2>=1.0 || r2==0.0);
      //--- Box-Muller transform
      f=MathSqrt(-2.0*MathLog(r2)/r2);
      double snorm=f*x2;
      //--- normal sample
      g=lambda+s*snorm;

      if(0.0<=g)
        {
         value=(int)(g);
         //--- immediate acceptance if large enough
         if(l<=value)
            return value;
         //--- squeeze acceptance
         fk=(double)(value);
         difmuk=lambda-fk;
         u=MathRandomNonZero();
         //---
         if(difmuk*difmuk*difmuk<=d*u)
            return value;
        }
      //--- preparation for steps P and Q
      double omega=0.3989423/s;
      double b1 = 0.04166667/lambda;
      double b2 = 0.3*b1*b1;
      double c3 = 0.1428571*b1*b2;
      double c2 = b2 - 15.0*c3;
      double c1 = b1 - 6.0*b2 + 45.0*c3;
      double c0 = 1.0 - b1 + 3.0*b2 - 15.0*c3;
      double c=0.1069/lambda;
      double del=0;

      if(0.0<=g)
        {
         kflag=0;

         if(value<10)
           {
            px = -lambda;
            py = MathPow(lambda,value)/MathFactorial(value);
           }
         else
           {
            del = 0.8333333E-01/fk;
            del = del - 4.8*del*del*del;
            v=difmuk/fk;

            if(0.25<MathAbs(v))
              {
               px=fk*MathLog(1.0+v)-difmuk-del;
              }
            else
              {
               px=fk*v*v*(((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0)-del;
              }
            py=0.3989423/MathSqrt(fk);
           }
         x=(0.5-difmuk)/s;
         xx = x * x;
         fx = -0.5 * xx;
         fy = omega*(((c3*xx+c2)*xx+c1)*xx+c0);

         if(kflag<=0)
           {
            if(fy-u*fy<=py*MathExp(px-fx))
               return value;
           }
         else
           {
            if(c*MathAbs(u)<=py*MathExp(px+e)-fy*MathExp(fx+e))
               return value;
           }
        }
      //--- exponential sample
      for(;;)
        {
         double rnd=MathRandomNonZero();
         e=-MathLog(1.0-rnd);

         u=2.0*MathRandomNonZero()-1.0;
         if(u<0.0)
            t=1.8-MathAbs(e);
         else
            t=1.8+MathAbs(e);

         if(t<=-0.6744)
            continue;

         value=(int)(lambda+s*t);
         fk=(double)(value);
         difmuk=lambda-fk;

         kflag=1;
         //--- calculation of PX, PY, FX, FY
         if(value<10)
           {
            px = -lambda;
            py = MathPow(lambda,value)/MathFactorial(value);
           }
         else
           {
            del = 0.8333333E-01/fk;
            del = del - 4.8*del*del*del;
            v=difmuk/fk;

            if(0.25<MathAbs(v))
               px=fk*MathLog(1.0+v)-difmuk-del;
            else
               px=fk*v*v*(((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v+a0)-del;

            py=0.3989423/MathSqrt(fk);
           }

         x=(0.5-difmuk)/s;
         xx = x*x;
         fx = -0.5*xx;
         fy = omega*(((c3*xx+c2)*xx+c1)*xx+c0);

         if(kflag<=0)
           {
            if(fy-u*fy<=py*MathExp(px-fx))
               return value;
           }
         else
           {
            if(c*MathAbs(u)<=py*MathExp(px+e)-fy*MathExp(fx+e))
               return value;
           }
        }
     }
   return value;
  }
//+------------------------------------------------------------------+
//| Random variate from the Poisson distribution                     |
//+------------------------------------------------------------------+
//| Compute the random variable from the Poisson distribution        |
//| with parameter lambda.                                           |
//|                                                                  |
//| Arguments                                                        |
//| lambda      : Mean                                               |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Poisson distribution.                      |
//+------------------------------------------------------------------+
double MathRandomPoisson(const double lambda,int &error_code)
  {
//--- check parameters
   if(!MathIsValidNumber(lambda))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- lambda must be positive
   if(lambda<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
   error_code=ERR_OK;
   return MathRandomPoisson(lambda);
  }
//+------------------------------------------------------------------+
//| Random variate from the Poisson distribution                     |
//+------------------------------------------------------------------+
//| Generates random variables from the Poisson distribution         |
//| with parameter lambda.                                           |
//|                                                                  |
//| Arguments:                                                       |
//| lambda     : Mean                                                |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomPoisson(const double lambda,const int data_count,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(lambda))
      return false;
//--- lambda must be positive
   if(lambda<=0.0)
      return false;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      result[i]=MathRandomPoisson(lambda);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Poisson distribution moments                                     |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of the Poisson           |
//| distribution with parameter lambda.                              |
//|                                                                  |
//| Arguments:                                                       |
//| lambda     : Mean                                                |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsPoisson(const double lambda,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(lambda))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- lambda must be positive
   if(lambda<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- calculate moments
   mean    =lambda;
   variance=lambda;
   skewness=MathPow(lambda,-0.5);
   kurtosis=1.0/lambda;
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
