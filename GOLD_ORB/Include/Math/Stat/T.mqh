//+------------------------------------------------------------------+
//|                                                            T.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"
#include "Gamma.mqh"

//+------------------------------------------------------------------+
//| T probability density function (PDF)                             |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the T-distribution with parameter nu.                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| nu         : Degrees of freedom                                  |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityT(const double x,const double nu,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(nu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check nu
   if(nu!=MathRound(nu) || nu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- calculate T density
   double pdf=MathExp(MathGammaLog((nu+1.0)*0.5)-MathGammaLog(nu*0.5));
   pdf=pdf/(MathSqrt(nu*M_PI)*MathPow(1+x*x/nu,(nu+1.0)*0.5));
//--- return density
   return TailLogValue(pdf,true,log_mode);
  }
//+------------------------------------------------------------------+
//| T probability density function (PDF)                             |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the T-distribution with parameter nu.                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| nu         : Degrees of freedom                                  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityT(const double x,const double nu,int &error_code)
  {
   return MathProbabilityDensityT(x,nu,false,error_code);
  }
//+------------------------------------------------------------------+
//| T probability density function (PDF)                             |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the T distribution with parameter nu for values in x[] array.    |
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
bool MathProbabilityDensityT(const double &x[],const double nu,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu))
      return false;
//--- check nu
   if(nu!=MathRound(nu) || nu<=0.0)
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

      error_code=ERR_OK;
      //--- calculate T density
      double pdf=MathExp(MathGammaLog((nu+1.0)*0.5)-MathGammaLog(nu*0.5));
      pdf=pdf/(MathSqrt(nu*M_PI)*MathPow(1+x_arg*x_arg/nu,(nu+1.0)*0.5));
      //--- return density
      result[i]=TailLogValue(pdf,true,log_mode);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| T probability density function (PDF)                             |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the T distribution with parameter nu for values in x[] array.    |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityT(const double &x[],const double nu,double &result[])
  {
   return MathProbabilityDensityT(x,nu,false,result);
  }
//+------------------------------------------------------------------+
//| T cumulative distribution function (CDF)                         |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation from    |
//| T-distribution with parameter nu is less than or equal to x.     |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| nu         : Degrees of freedom                                  |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the T cumulative distribution function with         |
//| parameter nu, evaluated at x.                                    |
//+------------------------------------------------------------------+
double MathCumulativeDistributionT(const double x,const double nu,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(nu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check nu (must be positive integer)
   if(nu!=MathRound(nu) || nu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- special case
   if(nu==1.0)
      return TailLogValue(0.5+MathArctan(x)/M_PI,tail,log_mode);
//--- otherwise
   if(x==0)
      return TailLogValue(0.5,tail,log_mode);
//--- calculate pdf using incomplete Beta
   double cdf=1.0-MathBetaIncomplete(nu/(nu+x*x),nu*0.5,0.5);
   cdf=(1.0-cdf)*0.5;
//--- check x
   if(x>0.0)
      cdf=1.0-cdf;
//--- take into account round-off errors for probability
   return TailLogValue(MathMin(cdf,1.0),tail,log_mode);
  }
//+------------------------------------------------------------------+
//| T cumulative distribution function (CDF)                         |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation from    |
//| T-distribution with parameter nu is less than or equal to x.     |
//|                                                                  |
//| Arguments:                                                       |
//| x          : The desired quantile                                |
//| nu         : Degrees of freedom                                  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of T cumulative distribution function with parameter   |
//| nu, evaluated at x.                                              |
//+------------------------------------------------------------------+
double MathCumulativeDistributionT(const double x,const double nu,int &error_code)
  {
   return MathCumulativeDistributionT(x,nu,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| T cumulative distribution function (CDF)                         |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the T distribution with parameter nu for values in x.            |
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
bool MathCumulativeDistributionT(const double &x[],const double nu,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu))
      return false;
//--- check nu (must be positive integer)
   if(nu!=MathRound(nu) || nu<=0.0)
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

      //--- special case
      if(nu==1.0)
         result[i]=TailLogValue(0.5+MathArctan(x_arg)/M_PI,tail,log_mode);
      else
      //--- otherwise
      if(x_arg==0)
         result[i]=TailLogValue(0.5,tail,log_mode);
      else
        {
         //--- calculate pdf using incomplete Beta
         double cdf=1.0-MathBetaIncomplete(nu/(nu+x_arg*x_arg),nu*0.5,0.5);
         cdf=(1.0-cdf)*0.5;
         //--- check x
         if(x_arg>0.0)
            cdf=1.0-cdf;
         //--- take into account round-off errors for probability
         result[i]=TailLogValue(MathMin(cdf,1.0),tail,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| T cumulative distribution function (CDF)                         |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the T distribution with parameter nu for values in x[] array.    |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| nu         : Degrees of freedom                                  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionT(const double &x[],const double nu,double &result[])
  {
   return MathCumulativeDistributionT(x,nu,true,false,result);
  }
//+------------------------------------------------------------------+
//| T distribution quantile function (inverse CDF)                   |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the T distribution with parameter nu for the desired |
//| probability.                                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| nu          : Degrees of freedom                                 |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the T-distribution with parameter nu.                         |
//+------------------------------------------------------------------+
double MathQuantileT(const double probability,const double nu,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(probability) || !MathIsValidNumber(nu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check nu
   if(nu!=MathRound(nu) || nu<0.0)
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

//--- check cases when probability==0 or 1
   if(prob==0.0 || prob==1.0)
     {
      error_code=ERR_RESULT_INFINITE;
      //---
      if(prob==0.0)
         return(QNEGINF);
      else
         return(QPOSINF);
     }

   error_code=ERR_OK;
//--- special case nu=1
   if(nu==1.0)
      return MathTan(M_PI*(prob-0.5));
//--- special case
   if(prob==0.5)
      return 0.0;
//---
   int    max_iterations=50;
   int    iterations=0;
//--- initial values
   double h=1.0;
   double h_min=10E-20;
   double x=0.5;
   int    err_code=0;
//--- Newton iterations
   while(iterations<max_iterations)
     {
      //--- check convegence
      if((MathAbs(h)>h_min && MathAbs(h)>MathAbs(h_min*x))==false)
         break;
      //--- calculate pdf and cdf
      double pdf=MathProbabilityDensityT(x,nu,err_code);
      double cdf=MathCumulativeDistributionT(x,nu,err_code);
      //--- calculate ratio
      h=(cdf-prob)/pdf;
      //---
      double x_new=x-h;
      //--- check x
      if(x_new<0.0)
         x_new=x*0.1;
      else
      if(x_new>1.0)
         x_new=1.0-(1.0-x)*0.1;

      x=x_new;

      iterations++;
     }
//--- check convergence
   if(iterations<max_iterations)
      return x;
   else
     {
      error_code=ERR_NON_CONVERGENCE;
      return QNaN;
     }
  }
//+------------------------------------------------------------------+
//| T distribution quantile function (inverse CDF)                   |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the T distribution with parameter nu for the desired |
//| probability.                                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| nu          : Degrees of freedom                                 |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the T-distribution with parameter nu.                         |
//+------------------------------------------------------------------+
double MathQuantileT(const double probability,const double nu,int &error_code)
  {
   return MathQuantileT(probability,nu,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| T distribution quantile function (inverse CDF)                   |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the T distribution with parameter nu for             |
//| values from the probability[] array.                             |
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
bool MathQuantileT(const double &probability[],const double nu,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu))
      return false;
//--- check nu
   if(nu!=MathRound(nu) || nu<0.0)
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

      //--- special case p=0.5
      if(prob==0.5)
         result[i]=0.0;
      else
      if(prob==0.0)
         result[i]=QNEGINF;
      else
      if(prob==1.0)
         result[i]=QPOSINF;
      else
        {
         //--- special case nu=1
         if(nu==1.0)
            result[i]=MathTan(M_PI*(prob-0.5));
         else
           {
            int    max_iterations=50;
            int    iterations=0;
            //--- initial values
            double h=1.0;
            double h_min=10E-18;
            double x=0.5;
            int    err_code=0;
            //--- Newton iterations
            while(iterations<max_iterations)
              {
               //--- check convegence     
               if((MathAbs(h)>h_min && MathAbs(h)>MathAbs(h_min*x))==false)
                  break;
               //--- calculate pdf and cdf
               double pdf=MathProbabilityDensityT(x,nu,err_code);
               double cdf=MathCumulativeDistributionT(x,nu,err_code);
               //--- calculate ratio
               h=(cdf-prob)/pdf;
               //---
               double x_new=x-h;
               //--- check x
               if(x_new<0.0)
                  x_new=x*0.1;
               else
               if(x_new>1.0)
                  x_new=1.0-(1.0-x)*0.1;
                  
               if (MathAbs(x_new-x)<10E-15)
               break;

               x=x_new;

               iterations++;
              }
            //--- check convergence
            if(iterations<max_iterations)
               result[i]=x;
            else
               return false;
           }
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| T distribution quantile function (inverse CDF)                   |
//+------------------------------------------------------------------+
//| The function calculates the inverse cumulative distribution      |
//| function of the T distribution with parameter nu for             |
//| values from the probability[] array.                             |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| nu          : Degrees of freedom                                 |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileT(const double &probability[],const double nu,double &result[])
  {
   return MathQuantileT(probability,nu,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the T distribution                           |
//+------------------------------------------------------------------+
//| Computes the random variable from the T distribution             |
//| with parameter nu.                                               |
//|                                                                  |
//| Arguments:                                                       |
//| nu          : Degrees of freedom                                 |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with T distribution.                            |
//+------------------------------------------------------------------+
double MathRandomT(const double nu,int error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(nu))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check arguments
   if(nu!=MathRound(nu) || nu<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- calculate normal random variable using Box-Muller transform
   double x1,x2,r2;
   do
     {
      x1=2.0*MathRandomNonZero()-1.0;
      x2=2.0*MathRandomNonZero()-1.0;
      r2=x1*x1+x2*x2;
     }
   while(r2>=1.0 || r2==0.0);
//--- generate normal and gamma random variables 
   double rnd_normal=x2*MathSqrt(-2.0*MathLog(r2)/r2);
   double rnd_gamma=MathRandomGamma(nu*0.5,1,error_code);
//--- calculate ratio
   double result=0;
   if(rnd_gamma!=0)
      result=MathSqrt(nu*0.5)*rnd_normal/MathSqrt(rnd_gamma);
   return(result);
  }
//+------------------------------------------------------------------+
//| Random variate from the T distribution                           |
//+------------------------------------------------------------------+
//| Generates random variables from the T distribution with          |
//| parameter nu.                                                    |
//|                                                                  |
//| Arguments:                                                       |
//| nu         : Degrees of freedom                                  |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomT(const double nu,const int data_count,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(nu))
      return false;
//--- check arguments
   if(nu!=MathRound(nu) || nu<=0.0)
      return false;

   int error_code=0;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      //--- calculate normal random variable using Box-Muller transform
      double x1,x2,r2;
      do
        {
         x1=2.0*MathRandomNonZero()-1.0;
         x2=2.0*MathRandomNonZero()-1.0;
         r2=x1*x1+x2*x2;
        }
      while(r2>=1.0 || r2==0.0);
      //--- generate normal and gamma random variables 
      double rnd_normal=x2*MathSqrt(-2.0*MathLog(r2)/r2);
      double rnd_gamma=MathRandomGamma(nu*0.5,1,error_code);
      //--- calculate ratio
      double rnd=0;
      if(rnd_gamma!=0)
         rnd=MathSqrt(nu*0.5)*rnd_normal/MathSqrt(rnd_gamma);
      result[i]=rnd;
     }
   return true;
  }
//+------------------------------------------------------------------+
//| T distribution moments                                           |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of the T distribution    |
//| with parameter nu.                                               |
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
double MathMomentsT(const double nu,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
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
      return QNaN;
     }
//--- check nu
   if(nu!=MathRound(nu) || nu<0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- calculate moments
   mean=0;
   if(nu>2)
      variance=nu/(nu-2);
   skewness=0;
   if(nu>4)
      kurtosis=6/(nu-4);
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
