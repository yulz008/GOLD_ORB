//+------------------------------------------------------------------+
//|                                                        Gamma.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Normal.mqh"

const double DoubleEpsilon=1.11022302462515654042E-16;
const double LogMax=7.09782712893383996732E2;
//+------------------------------------------------------------------+
//| Inverse of the incomplete Gamma integral                         |
//+------------------------------------------------------------------+
double MathInverseGammaIncomplete(const double a,const double y)
  {
//--- bound the solution
   double x0 = DBL_MAX;
   double yl = 0;
   double x1 = 0;
   double yh = 1.0;
   double dithresh=5.0*DoubleEpsilon;
//--- approximation to inverse function
   double d=1.0/(9.0*a);
   int err_code=0;
   double q_normal=MathQuantileNormal(y,0,1,true,false,err_code);
   double yy=(1.0-d-q_normal*MathSqrt(d));
   double x=a*yy*yy*yy;
   double lgm=MathGammaLog(a);
   for(int i=0; i<10; i++)
     {
      if(x>x0 || x<x1)
         break;
      yy=1.0-MathGammaIncomplete(x,a);
      if(yy<yl || yy>yh)
         break;
      if(yy<y)
        {
         x0 = x;
         yl = yy;
        }
      else
        {
         x1 = x;
         yh = yy;
        }
      //--- compute the derivative of the function at this point
      d=(a-1.0)*MathLog(x)-x-lgm;
      if(d<-LogMax)
         break;
      d=-MathExp(d);
      //--- compute the step to the next approximation of x
      d=(yy-y)/d;
      if(MathAbs(d/x)<DoubleEpsilon)
         return (x);
      x=x-d;
     }
//--- resort to interval halving if Newton iteration did not converge. 
   d=0.0625;
   if(x0==DBL_MAX)
     {
      if(x<=0.0)
         x=1.0;
      while(x0==DBL_MAX && MathIsValidNumber(x))
        {
         x=(1.0+d)*x;
         yy=1.0-MathGammaIncomplete(x,a);
         if(yy<y)
           {
            x0 = x;
            yl = yy;
            break;
           }
         d=d+d;
        }
     }
   d=0.5;
   double dir=0;
   for(int i=0; i<400; i++)
     {
      double t=x1+d *(x0-x1);
      if(!MathIsValidNumber(t))
         break;
      x=t;
      yy=1.0-MathGammaIncomplete(x,a);
      lgm=(x0-x1)/(x1+x0);
      if(MathAbs(lgm)<dithresh)
         break;
      lgm=(yy-y)/y;
      if(MathAbs(lgm)<dithresh)
         break;
      if(x<=0.0)
         break;
      if(yy>=y)
        {
         x1 = x;
         yh = yy;
         if(dir<0)
           {
            dir=0;
            d=0.5;
           }
         else
         if(dir>1)
            d=0.5*d+0.5;
         else
            d=(y-yl)/(yh-yl);
         dir+=1;
        }
      else
        {
         x0 = x;
         yl = yy;
         if(dir>0)
           {
            dir=0;
            d=0.5;
           }
         else
         if(dir<-1)
            d=0.5*d;
         else
            d=(y-yl)/(yh-yl);
         dir-=1;
        }
     }
   if(x==0.0 || !MathIsValidNumber(x))
     {
      Print("Errors in an arithmetic, casting, or conversion operation.");
      return(QNaN);
     }
//---
   return(x);
  }
//+------------------------------------------------------------------+
//| Gamma probability density function (PDF)                         |
//+------------------------------------------------------------------+
//| The function returns the probability density function of         |
//| of the Gamma distribution with shape parameters a and b.         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| a          : Shape                                               |
//| b          : Scale                                               |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityGamma(const double x,const double a,const double b,const bool log_mode,int &error_code)
  {
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
//--- check negative x
   if(x<=0)
      return TailLog0(true,log_mode);
//--- calculate log Gamma density
   double log_result=(a-1.0)*MathLog(x)-(x/b)-MathGammaLog(a)-a*MathLog(b);
   if(log_mode==true)
      return(log_result);
//--- return Gamma density
   return MathExp(log_result);
  }
//+------------------------------------------------------------------+
//| Gamma probability density function (PDF)                         |
//+------------------------------------------------------------------+
//| The function returns the probability density function of         |
//| of the Gamma distribution with shape parameters a and b.         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| a          : Shape                                               |
//| b          : Scale                                               |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityGamma(const double x,const double a,const double b,int &error_code)
  {
   return MathProbabilityDensityGamma(x,a,b,false,error_code);
  }
//+------------------------------------------------------------------+
//| Gamma probability density function (PDF)                         |
//+------------------------------------------------------------------+
//| The function calculates the Gamma probability density function   |
//| with parameters a and b for values in x[] array.                 |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Shape                                               |
//| b          : Scale                                               |
//| log_mode   : Logarithm mode flag,if true it calculates Log values|
//| result     : Output array for calculated values                  |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathProbabilityDensityGamma(const double &x[],const double a,const double b,const bool log_mode,double &result[])
  {
//--- check parameters
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

      if(!MathIsValidNumber(x_arg))
         return false;

      if(x_arg>0)
        {
         //--- calculate log Gamma density
         double log_result=(a-1.0)*MathLog(x_arg)-(x_arg/b)-MathGammaLog(a)-a*MathLog(b);
         if(log_mode==true)
            result[i]=log_result;
         else
            result[i]=MathExp(log_result);
        }
      else
         result[i]=TailLog0(true,log_mode);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Gamma probability density function (PDF)                         |
//+------------------------------------------------------------------+
//| The function calculates the Gamma probability density function   |
//| with parameters a and b for values from x[] array.               |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Shape                                               |
//| b          : Scale                                               |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathProbabilityDensityGamma(const double &x[],const double a,const double b,double &result[])
  {
   return MathProbabilityDensityGamma(x,a,b,false,result);
  }
//+------------------------------------------------------------------+
//| Gamma cumulative distribution function (CDF)                     |
//+------------------------------------------------------------------+
//| The function returns the cumulative distribution function of the |
//| Gamma distribution with parameters a and b.                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| a          : Shape                                               |
//| b          : Scale                                               |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode flag,if true it calculates Log values|
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Gamma cumulative distribution function          |
//| with parameters a and b, evaluated at x.                         |
//+------------------------------------------------------------------+
double MathCumulativeDistributionGamma(const double x,const double a,const double b,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
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
//--- calculate probability using Incomplete Gamma function and take into account round-off errors
   double cdf=MathMin(MathGammaIncomplete(x/b,a),1.0);
   return TailLogValue(cdf,tail,log_mode);
  }
//+------------------------------------------------------------------+
//| Gamma cumulative distribution function (CDF)                     |
//+------------------------------------------------------------------+
//| The function returns the cumulative distribution function of the |
//| Gamma distribution with parameters a and b.                      |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| a          : Shape                                               |
//| b          : Scale                                               |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Gamma cumulative distribution function          |
//| with parameters a and b, evaluated at x.                         |
//+------------------------------------------------------------------+
double MathCumulativeDistributionGamma(const double x,const double a,const double b,int &error_code)
  {
   return MathCumulativeDistributionGamma(x,a,b,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Gamma cumulative distribution function (CDF)                     |
//+------------------------------------------------------------------+
//| The function calculates the values of the Gamma cumulative       |
//| distribution function with given a and b for values in x[] array.|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Shape                                               |
//| b          : Scale                                               |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode flag,if true it calculates Log values|
//| resut      : Output array for calculated values                  |
//|                                                                  |
//| Return value:                                                    |
//| true if successul, otherwise false.                              |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionGamma(const double &x[],const double a,const double b,const bool tail,const bool log_mode,double &result[])
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

      if(!MathIsValidNumber(x_arg))
         return false;

      if(x_arg<=0)
         result[i]=TailLog0(tail,log_mode);
      else
        {
         //--- calculate probability using Incomplete Gamma function and take into account round-off errors
         double cdf=MathMin(MathGammaIncomplete(x_arg/b,a),1.0);
         result[i]=TailLogValue(cdf,tail,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Gamma cumulative distribution function (CDF)                     |
//+------------------------------------------------------------------+
//| The function calculates the values of the Gamma cumulative       |
//| distribution function with given a and b for values in x[] array.|
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : Shape                                               |
//| b          : Scale                                               |
//| result     : Output array for calculated values                  |
//|                                                                  |
//| Return value:                                                    |
//| true if successul, otherwise false.                              |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionGamma(const double &x[],const double a,const double b,double &result[])
  {
   return MathCumulativeDistributionGamma(x,a,b,true,false,result);
  }
//+------------------------------------------------------------------+
//| Gamma distribution quantile function (inverse CDF)               |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Gamma distribution with parameters a and b       |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| a           : Shape                                              |
//| b           : Scale                                              |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Gamma distribution with parameters a and b.               |
//+------------------------------------------------------------------+
double MathQuantileGamma(const double probability,const double a,const double b,const bool tail,const bool log_mode,int &error_code)
  {
//--- case log probability==-inf
   if(log_mode==true && probability==QNEGINF)
     {
      error_code=ERR_OK;
      return 0.0;
     }
//--- check NaN
   if(!MathIsValidNumber(probability) || !MathIsValidNumber(a) || !MathIsValidNumber(b))
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
   error_code=ERR_OK;
//--- case probability==0
   if(prob==0.0)
      return 0.0;
//--- case probability==1
   if(prob==1.0)
     {
      error_code=ERR_RESULT_INFINITE;
      return QPOSINF;
     }
//--- calculate quantile
   double quantile=MathInverseGammaIncomplete(a,1.0-prob)*b;
   if(!MathIsValidNumber(quantile))
      error_code=ERR_NON_CONVERGENCE;
//---
   return(quantile);
  }
//+------------------------------------------------------------------+
//| Gamma distribution quantile function (inverse CDF)               |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Gamma distribution with parameters a and b       |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| a           : Shape                                              |
//| b           : Scale                                              |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Gamma distribution with parameters a and b.               |
//+------------------------------------------------------------------+
double MathQuantileGamma(const double probability,const double a,const double b,int &error_code)
  {
   return MathQuantileGamma(probability,a,b,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Gamma distribution quantile function (inverse CDF)               |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Gamma distribution with parameters a and b       |
//| for values from the probability[] array.                         |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| a           : Shape                                              |
//| b           : Scale                                              |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode,if true it calculates for Log values|
//| result      : Output array for calculated values                 |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileGamma(const double &probability[],const double a,const double b,const bool tail,const bool log_mode,double &result[])
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

   const double eps=10E-18;
   double max_h=MathSqrt(eps);
   const int max_iterations=1000;

   for(int i=0; i<data_count; i++)
     {
      //--- calculate real probability
      double prob=TailLogProbability(probability[i],tail,log_mode);

      if(!MathIsValidNumber(prob))
         return false;

      //--- check probability range
      if(prob<0.0 || prob>1.0)
         return false;

      //--- case probability==0
      if(prob==0.0)
         result[i]=0.0;
      else
      //--- case probability==1
      if(prob==1.0)
         result[i]=QPOSINF;
      else
        {
         double quantile=MathInverseGammaIncomplete(a,1.0-prob)*b;
         if(MathIsValidNumber(quantile))
            result[i]=quantile;
         else
            return false;
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Gamma distribution quantile function (inverse CDF)               |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Gamma distribution with parameters a and b       |
//| for the desired probability.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| a           : Shape                                              |
//| b           : Scale                                              |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileGamma(const double &probability[],const double a,const double b,double &result[])
  {
   return MathQuantileGamma(probability,a,b,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Gamma distribution                       |
//+------------------------------------------------------------------+
//| Compute the random variable from the Gamma distribution          |
//| with parameters a and b.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| a          : Shape                                               |
//| b          : Scale                                               |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Gamma distribution.                        |
//+------------------------------------------------------------------+
//| Author: Robert Kern                                              |
//+------------------------------------------------------------------+
double MathRandomGamma(const double a,const double b)
  {
   double bb,c,U,V,X=0,Y;
//--- check shape
   if(a==1.0)
     {
      //--- exponential
      return -MathLog(1.0-MathRandomNonZero());
     }
   else
   if(a<1.0)
     {
      for(;;)
        {
         U=MathRandomNonZero();
         //--- exponential
         V=-MathLog(1.0-MathRandomNonZero());

         if(U<=1.0-a)
           {
            X=MathPow(U,1.0/a);
            if(X<=V)
               return b*X;
           }
         else
           {
            Y = -MathLog((1-U)/a);
            X = MathPow(1.0 - a + a*Y, 1.0/a);
            if(X<=(V+Y))
               return(b*X);
           }
        }
     }
   else
     {
      bb= a-1.0/3.0;
      c = 1.0/MathSqrt(9*bb);
      for(;;)
        {
         do
           {
            //--- generate normal random variate
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
            X=f*x2;

            V=1.0+c*X;
           }
         while(V<=0.0);

         V = V*V*V;
         U = MathRandomNonZero();

         if(U<1.0-0.0331*(X*X)*(X*X))
            return(bb*V*b);

         if(MathLog(U)<0.5*X*X+bb*(1.0-V+MathLog(V)))
            return(bb*V*b);
        }
     }
   return(X*b);
  }
//+------------------------------------------------------------------+
//| Random variate from the Gamma distribution                       |
//+------------------------------------------------------------------+
//| Compute the random variable from the Gamma distribution          |
//| with parameters a and b.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| a          : Shape                                               |
//| b          : Scale                                               |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Gamma distribution.                        |
//+------------------------------------------------------------------+
double MathRandomGamma(const double a,const double b,int &error_code)
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
   return MathRandomGamma(a,b);
  }
//+------------------------------------------------------------------+
//| Random variate from the Gamma distribution                       |
//+------------------------------------------------------------------+
//| The function generates random variables from the Gamma           |
//| distribution with parameters a and b.                            |
//|                                                                  |
//| Arguments:                                                       |
//| a           : First shape parameter (a>0)                        |
//| b           : Second shape parameter (b>0)                       |
//| data_count  : Number of values needed                            |
//| result      : Output array for random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomGamma(const double a,const double b,const int data_count,double &result[])
  {
   if(data_count<=0)
      return false;
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- a and b must be positive
   if(a<=0 || b<=0)
      return false;
//--- prepare output array and calculate values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
      result[i]=MathRandomGamma(a,b);
   return true;
  }
//+------------------------------------------------------------------+
//| Gamma distribution moments                                       |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of Gamma distribution    |
//| with parameters a and b.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| a          : Shape                                               |
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
bool MathMomentsGamma(const double a,const double b,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
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
//--- calculate moments
   mean    =a*b;
   variance=a*b*b;
   skewness=2/MathSqrt(a);
   kurtosis=6/a;
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
