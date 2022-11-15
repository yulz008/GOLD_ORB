//+------------------------------------------------------------------+
//|                                                     Binomial.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"
#include "Beta.mqh"

//+------------------------------------------------------------------+
//| Binomial probability mass function (PDF)                         |
//+------------------------------------------------------------------+
//| The function returns the Binomial probability mass function      |
//| with parameters n and p at x.                                    |
//|                                                                  |
//|             f(x,n,p)= C(n,x)*(p^x)*(1-p)^(n-x)                   |
//|                                                                  |
//| where binomial coefficient C(n,k)=n!/(k!*(n-k)!)                 |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Integer random variable                             |
//| n          : Number of trials                                    |
//| p          : Probability of success for each trial               |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability mass function evaluated at x.                    |
//+------------------------------------------------------------------+
double MathProbabilityDensityBinomial(const double x,const double n,const double p,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(n) || !MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check n
   if(n<0 || n!=MathRound(n))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check p range
   if(p<0.0 || p>1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- case p=0
   if(p==0.0 || p==1.0)
      return TailLog0(true,log_mode);
//--- check x range
   if(x<0 || x>n)
      return TailLog0(true,log_mode);

   double log_result=MathGammaLog(n+1.0)-MathGammaLog(x+1.0)-MathGammaLog(n-x+1.0)+x*MathLog(p)+(n-x)*MathLog(1.0-p);
   if(log_mode==true)
      return log_result;
//--- return probability mass
   return MathExp(log_result);
  }
//+------------------------------------------------------------------+
//| Binomial probability mass function (PDF)                         |
//+------------------------------------------------------------------+
//| The function returns the Binomial probability mass function      |
//| with parameters n and p at x.                                    |
//|                                                                  |
//|             f(x,n,p)= C(n,x)*(p^x)*(1-p)^(n-x)                   |
//|                                                                  |
//| where binomial coefficient C(n,k)=n!/(k!*(n-k)!)                 |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Integer random variable                             |
//| n          : Number of trials                                    |
//| p          : Probability of success for each trial               |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability mass function evaluated at x.                    |
//+------------------------------------------------------------------+
double MathProbabilityDensityBinomial(const double x,const double n,const double p,int &error_code)
  {
   return MathProbabilityDensityBinomial(x,n,p,false,error_code);
  }
//+------------------------------------------------------------------+
//| Binomial probability mass function (PDF)                         |
//+------------------------------------------------------------------+
//| The function calculates the Binomial probability mass function   |
//| with parameters n and p for values in x[] array.                 |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with integer random variables                 |
//| n          : Number of trials                                    |
//| p          : Probability of success for each trial               |
//| log_mode   : Logarithm mode flag,if true it calculates Log values|
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathProbabilityDensityBinomial(const double &x[],const double n,const double p,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(n) || !MathIsValidNumber(p))
      return false;
//--- check n
   if(n<0 || n!=MathRound(n))
      return false;
//--- check p range
   if(p<0.0 || p>1.0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);
//--- case p=0 or p=1
   if(p==0.0 || p==1.0)
     {
      for(int i=0; i<data_count; i++)
         result[i]=TailLog0(true,log_mode);
      return true;
     }
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];
      if(MathIsValidNumber(x_arg) && x_arg==MathRound(x_arg))
        {
         //--- check x range
         if(x_arg<0 || x_arg>n)
            result[i]=TailLog0(true,log_mode);
         else
           {
            double log_result=MathGammaLog(n+1.0)-MathGammaLog(x_arg+1.0)-MathGammaLog(n-x_arg+1.0)+x_arg*MathLog(p)+(n-x_arg)*MathLog(1.0-p);
            if(log_mode==true)
               result[i]=log_result;
            else
               result[i]=MathExp(log_result);
           }
        }
      else
         return false;
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Binomial probability mass function (PDF)                         |
//+------------------------------------------------------------------+
//| The function calculates the Binomial probability mass function   |
//| with parameters n and p for values in x[] array.                 |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with integer random variables                 |
//| n          : Number of trials                                    |
//| p          : Probability of success for each trial               |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathProbabilityDensityBinomial(const double &x[],const double n,const double p,double &result[])
  {
   return MathProbabilityDensityBinomial(x,n,p,false,result);
  }
//+------------------------------------------------------------------+
//| Binomial cumulative distribution function (CDF)                  |
//+------------------------------------------------------------------+
//| The function returns the value of the Binomial cumulative        |
//| distribution function with given n and p at the desired x.       |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Integer random variable                             |
//| n          : Number of trials                                    |
//| p          : Probability of success for each trial               |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode flag,if true it calculates Log values|
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The cumulative distribution function evaluated at x.             |
//+------------------------------------------------------------------+
double MathCumulativeDistributionBinomial(const double x,const double n,double p,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(x) || !MathIsValidNumber(n) || !MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check n
   if(n<0 || n!=MathRound(n))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check probability
   if(p<0.0 || p>1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- case p==0
   if(p==0.0)
     {
      if(x>=0)
         return TailLog1(tail,log_mode);
      else
         return TailLog0(tail,log_mode);
     }
//--- case p==1
   if(p==1.0)
     {
      if(x>n)
         return TailLog1(tail,log_mode);
      else
         return TailLog0(tail,log_mode);
     }
//--- x must be>=0
   if(x<0)
      return TailLog0(tail,log_mode);
//--- check x
   if(x>n)
      return TailLog1(tail,log_mode);
   int err_code=0;
//--- calculate using Beta distribution and correct round-off errors
   double result=MathMin(1.0-MathCumulativeDistributionBeta(p,x+1.0,n-x,err_code),1.0);
//--- return result depending on arguments
   return TailLogValue(result,tail,log_mode);
  }
//+------------------------------------------------------------------+
//| Binomial cumulative distribution function (CDF)                  |
//+------------------------------------------------------------------+
//| The function returns the value of the Binomial cumulative        |
//| distribution function with given n and p at the desired x.       |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Integer random variable                             |
//| n          : Number of trials                                    |
//| p          : Probability of success for each trial               |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The cumulative distribution function evaluated at x.             |
//+------------------------------------------------------------------+
double MathCumulativeDistributionBinomial(const double x,const double n,double p,int &error_code)
  {
   return MathCumulativeDistributionBinomial(x,n,p,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Binomial cumulative distribution function (CDF)                  |
//+------------------------------------------------------------------+
//| The function returns the value of the Binomial cumulative        |
//| distribution function with given n and p at the desired x.       |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with integer random variables                 |
//| n          : Number of trials                                    |
//| p          : Probability of success for each trial               |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode flag,if true it calculates Log values|
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionBinomial(const double &x[],const double n,double p,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(n) || !MathIsValidNumber(p))
      return false;
//--- check n
   if(n<0 || n!=MathRound(n))
      return false;
//--- check probability
   if(p<0.0 || p>1.0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);

//--- case p=0 and p==1
   if(p==0.0 || p==1.0)
     {
      if(p==0.0)
        {
         for(int i=0; i<data_count; i++)
           {
            if(x[i]>=0)
               result[i]=TailLog1(tail,log_mode);
            else
               result[i]=TailLog0(tail,log_mode);
           }
        }
      else
      //--- p==1.0
        {
         for(int i=0; i<data_count; i++)
           {
            if(x[i]>n)
               result[i]=TailLog1(tail,log_mode);
            else
               result[i]=TailLog0(tail,log_mode);
           }
        }
      return true;
     }

   int err_code=0;
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];
      if(MathIsValidNumber(x_arg))
        {
         if(x_arg<0)
            result[i]=TailLog0(tail,log_mode);
         else
         if(x_arg>n)
            result[i]=TailLog1(tail,log_mode);
         else
           {
            double value=MathMin(1.0-MathCumulativeDistributionBeta(p,x_arg+1.0,n-x_arg,err_code),1.0);
            //--- calculate result depending on arguments
            result[i]=TailLogValue(value,tail,log_mode);
           }
        }
      else
         return false;
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Binomial cumulative distribution function (CDF)                  |
//+------------------------------------------------------------------+
//| The function returns the value of the Binomial cumulative        |
//| distribution function with given n and p for values              |
//| from x[] array.                                                  |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with integer random variables                 |
//| n          : Number of trials                                    |
//| p          : Probability of success for each trial               |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionBinomial(const double &x[],const double n,double p,double &result[])
  {
   return MathCumulativeDistributionBinomial(x,n,p,true,false,result);
  }
//+------------------------------------------------------------------+
//| Binomial distribution quantile function (inverse CDF)            |
//+------------------------------------------------------------------+
//| The function returns the value of the inverse Binomial cumulative|
//| distribution function with parameters n and p for the desired    |
//| probability.                                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| n           : Number of trials                                   |
//| p           : Probability of success for each trial              |
//| tail        : Lower tail flag (lower tail of probability used)   |
//| log_mode    : Logarithm mode flag (log probability used)         |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Binomial distribution with parameters n and p.            |
//+------------------------------------------------------------------+
double MathQuantileBinomial(const double probability,const double n,const double p,const bool tail,const bool log_mode,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(probability) || !MathIsValidNumber(n) || !MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check n
   if(n<0 || n!=MathRound(n))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check p range
   if(p<0.0 || p>1.0)
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
   int iterations=0;
   const int max_iterations=1000;
//--- direct cdf calculation
   double sum=MathProbabilityDensityBinomial(0,n,p,false,error_code);
   while(sum<prob && iterations<max_iterations)
     {
      sum+=MathProbabilityDensityBinomial(iterations,n,p,false,error_code);
      iterations++;
     }
//--- check convergence
   if(iterations<max_iterations)
     {
      if(iterations==0)
         return 0.0;
      else
         return iterations-1;
     }
   else
     {
      error_code=ERR_RESULT_INFINITE;
      return QPOSINF;
     }
  }
//+------------------------------------------------------------------+
//| Binomial distribution quantile function (inverse CDF)            |
//+------------------------------------------------------------------+
//| The function returns the value of the inverse Binomial cumulative|
//| distribution function with parameters n and p for the desired    |
//| probability.                                                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| n           : Number of trials                                   |
//| p           : Probability of success for each trial              |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function        |
//| of the Binomial distribution with parameters n and p.            |
//+------------------------------------------------------------------+
double MathQuantileBinomial(const double probability,const double n,const double p,int &error_code)
  {
   return MathQuantileBinomial(probability,n,p,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Binomial distribution quantile function (inverse CDF)            |
//+------------------------------------------------------------------+
//| The function calculates the value of the inverse Binomial        |
//| cumulative distribution function with parameters n and p for     |
//| the probability values from probability[] array.                 |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probability values                      |
//| n           : Number of trials                                   |
//| p           : Probability of success for each trial              |
//| tail        : Lower tail flag (lower tail of probability used)   |
//| log_mode    : Logarithm mode flag (log probability used)         |
//| result      : Output array with quantile values                  |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+ 
bool MathQuantileBinomial(const double &probability[],const double n,const double p,const bool tail,const bool log_mode,double &result[])
  {
//--- check NaN
   if(!MathIsValidNumber(n) || !MathIsValidNumber(p))
      return false;
//--- check n
   if(n<0 || n!=MathRound(n))
      return false;
//--- check p range
   if(p<0.0 || p>1.0)
      return false;

   int data_count=ArraySize(probability);
   if(data_count==0)
      return false;

   int error_code=0;
   ArrayResize(result,data_count);

   const int max_iterations=1000;
   for(int i=0; i<data_count; i++)
     {
      //--- calculate real probability
      double prob=TailLogProbability(probability[i],tail,log_mode);

      if(MathIsValidNumber(prob))
        {
         //--- check probability range
         if(prob<0.0 || prob>1.0)
            return false;

         int iterations=0;
         //--- direct cdf calculation
         double sum=MathProbabilityDensityBinomial(0,n,p,false,error_code);
         while(sum<prob && iterations<max_iterations)
           {
            sum+=MathProbabilityDensityBinomial(iterations,n,p,false,error_code);
            iterations++;
           }
         //--- check convergence
         if(iterations<max_iterations)
           {
            if(iterations==0)
               result[i]=0;
            else
               result[i]=iterations-1;
           }
         else
            return false;
        }
      else
         return false;
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Binomial distribution quantile function (inverse CDF)            |
//+------------------------------------------------------------------+
//| The function returns the value of the inverse Binomial cumulative|
//| distribution function with parameters n and p for the desired    |
//| probability values from probability[] array.                     |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probability values                      |
//| n           : Number of trials                                   |
//| p           : Probability of success for each trial              |
//| result      : Output array with quantile values                  |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+ 
bool MathQuantileBinomial(const double &probability[],const double n,const double p,double &result[])
  {
   return MathQuantileBinomial(probability,n,p,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from Binomial distribution                        |
//+------------------------------------------------------------------+
//| This procedure generates a single random deviate from Binomial   |
//| distribution whose number of trials is N and whose probability   |
//| of an event in each trial is p.                                  |
//|                                                                  |
//| Input parameters:                                                |
//| n : Number of binomial trials from which a random deviate        |
//|     will be generated.                                           |
//| p : The probability of an event in each trial of the binomial    |
//|     distribution from which a random deviate is to be generated. |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Binomial distribution.                     |
//|                                                                  |
//| Reference:                                                       |
//| Voratas Kachitvichyanukul, Bruce Schmeiser,                      |
//| "Binomial Random Variate Generation", Communications of the ACM, |
//| Volume 31, Number 2, February 1988, pages 216-222.               |
//|                                                                  |
//| Original FORTRAN77 version by Barry Brown, James Lovato.         |
//| C version by John Burkardt.                                      |
//+------------------------------------------------------------------+
double MathRandomBinomial(const double n,const double p)
  {
   int    ix,ix1,mp;
   double f,g,qn,r,t,u,v,w,w2,x,z;
   int    value=0;
   int    n1=(int)n;

   double pp= MathMin(p,1.0-p);
   double q = 1.0-pp;
   double xnp=(double)(n1)*pp;

   if(xnp<30.0)
     {
      qn= MathPow(q,n1);
      r = pp/q;
      g = r*(double)(n1+1);

      for(;;)
        {
         ix= 0;
         f = qn;
         u = MathRandomNonZero();

         for(;;)
           {
            if(u<f)
              {
               if(0.5<p)
                 {
                  ix=n1-ix;
                 }
               value=ix;
               return value;
              }
            if(110<ix)
               break;
            u=u-f;
            ix=ix+1;
            f=f*(g/(double)(ix)-r);
           }
        }
     }
   double ffm=xnp+pp;
   int m=int(ffm);
   double fm=m;
   double xnpq=xnp*q;
   double p1 = (int)(2.195*MathSqrt(xnpq)-4.6*q)+0.5;
   double xm = fm + 0.5;
   double xl = xm - p1;
   double xr = xm + p1;
   double c=0.134+20.5/(15.3+fm);
   double al=(ffm-xl)/(ffm-xl*pp);
   double xll=al*(1.0+0.5*al);
   al=(xr-ffm)/(xr*q);
   double xlr= al*(1.0 + 0.5*al);
   double p2 = p1*(1.0 + c + c);
   double p3 = p2 + c/xll;
   double p4 = p3 + c/xlr;
//--- generate a variate
   for(;;)
     {
      u = MathRandomNonZero()*p4;
      v = MathRandomNonZero();
      //--- triangle
      if(u<p1)
        {
         ix=int(xm-p1*v+u);
         if(0.5<p)
            ix=n1-ix;

         value=ix;
         return value;
        }
      //--- parallelogram
      if(u<=p2)
        {
         x = xl+(u - p1)/c;
         v = v*c + 1.0 - MathAbs(xm-x)/p1;

         if(v<=0.0 || 1.0<v)
            continue;
         ix=int(x);
        }
      else
      if(u<=p3)
        {
         ix=int(xl+MathLog(v)/xll);
         if(ix<0)
            continue;
         v=v*(u-p2)*xll;
        }
      else
        {
         ix=int(xr-MathLog(v)/xlr);
         if(n1<ix)
            continue;
         v=v*(u-p3)*xlr;
        }
      int k=MathAbs(ix-m);

      if(k<=20 || xnpq/2.0-1.0<=k)
        {
         f = 1.0;
         r = pp/q;
         g = (n1+1)*r;

         if(m<ix)
           {
            mp=m+1;
            for(int i=mp; i<=ix; i++)
               f=f*(g/i-r);
           }
         else
         if(ix<m)
           {
            ix1=ix+1;
            for(int i=ix1; i<=m; i++)
               f=f/(g/i-r);
           }

         if(v<=f)
           {
            if(0.5<p)
               ix=n1-ix;

            value=ix;
            return value;
           }
        }
      else
        {
         double amaxp=(k/xnpq)*((k*(k/3.0+0.625)+0.1666666666666)/xnpq+0.5);
         double ynorm=-double((k*k)/(2.0*xnpq));
         double alv=MathLog(v);

         if(alv<ynorm-amaxp)
           {
            if(0.5<p)
               ix=n1-ix;

            value=ix;
            return value;
           }

         if(ynorm+amaxp<alv)
            continue;

         double x1 = double(ix+1);
         double f1 = fm + 1.0;
         z = (double)(n1+1) - fm;
         w = (double)(n1-ix+1);
         double z2 = z * z;
         double x2 = x1 * x1;
         double f2 = f1 * f1;
         w2=w*w;

         t=xm*MathLog(f1/x1)+(n1-m+0.5)*MathLog(z/w)+(double)(ix-m)*MathLog(w*pp/(x1*q))
           +(13860.0 -(462.0 -(132.0 -(99.0-140.0/f2)/f2)/f2)/f2)/f1/166320.0
           +(13860.0 -(462.0 -(132.0 -(99.0-140.0/z2)/z2)/z2)/z2)/z/166320.0
           +(13860.0 -(462.0 -(132.0 -(99.0-140.0/x2)/x2)/x2)/x2)/x1/166320.0
           +(13860.0 -(462.0 -(132.0 -(99.0-140.0/w2)/w2)/w2)/w2)/w/166320.0;

         if(alv<=t)
           {
            if(0.5<p)
               ix=n1-ix;

            value=ix;
            return value;
           }
        }
     }
   return value;
  }
//+------------------------------------------------------------------+
//| Random variate from Binomial distribution                        |
//+------------------------------------------------------------------+
//| The function returns random deviate from Binomial distribution   |
//| with parameters n and p.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| n           : Number of trials                                   |
//| p           : Probability of success for each trial              |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Binomial distribution.                     |
//+------------------------------------------------------------------+
double MathRandomBinomial(const double n,const double p,int &error_code)
  {
//--- check NaN
   if(!MathIsValidNumber(n) || !MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- check n
   if(n<=0 || n!=MathRound(n))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }
//--- check probability
   if(p<=0 || p>=1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- return binomial random value
   return MathRandomBinomial(n,p);
  }
//+------------------------------------------------------------------+
//| Random variate from Binomial distribution                        |
//+------------------------------------------------------------------+
//| The function generates random variables from Binomial            |
//| distribution with parameters n and p.                            |
//|                                                                  |
//| Arguments:                                                       |
//| n           : Number of trials                                   |
//| p           : Probability of success for each trial              |
//| data_count  : Number of values needed                            |
//| result      : Output array with random values                    |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomBinomial(const double n,const double p,const int data_count,double &result[])
  {
   if(data_count<=0)
      return false;
//--- check NaN
   if(!MathIsValidNumber(n) || !MathIsValidNumber(p))
      return false;
//--- check n
   if(n<=0 || n!=MathRound(n))
      return false;
//--- check probability
   if(p<=0 || p>=1.0)
      return false;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
      result[i]=MathRandomBinomial(n,p);
   return true;
  }
//+------------------------------------------------------------------+
//| Binomial distriburion moments                                    |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of the Binomial          |
//| distribution with parameters n and p.                            |
//|                                                                  |
//| Arguments:                                                       |
//| n          : Number of trials                                    |
//| p          : Probability of success for each trial               |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
bool MathMomentsBinomial(const double n,const double p,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(n) || !MathIsValidNumber(p))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- check n
   if(n<0 || n!=MathRound(n))
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }
//--- check p range
   if(p<=0.0 || p>=1.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- prepare factors
   double np=n*p;
   double one_mp=(1.0-p);
//--- calculate moments
   mean    =np;
   variance=np*one_mp;
   skewness=(1-2*p)/MathSqrt(variance);
   kurtosis=(1-6*p*one_mp)/variance;
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
