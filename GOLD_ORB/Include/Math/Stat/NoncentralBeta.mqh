//+------------------------------------------------------------------+
//|                                               NoncentralBeta.mqh |
//|                   Copyright 2016-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#include "Math.mqh"
#include "Beta.mqh"
#include "NoncentralChiSquare.mqh"

//+------------------------------------------------------------------+
//| Noncental Beta density function (PDF)                            |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Noncental Beta distribution with parameters a,b,lambda    |
//|               Infinity                                           |
//| f(x,a,b,lambda)=Sum [p(k)*x^(a+k-1)*(1-x)^(b-1)]/Beta(a+k,b)     |
//|                 k=0                                              |
//|                                                                  |
//| where p(k)=(1/k!)*exp(-lambda/2)*(lambda/2)^k,                   |
//|       Beta(a,b)=Gamma(a)*Gamma(b)/Gamma(a+b)                     |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| a          : First shape parameter                               |
//| b          : Second shape parameter                              |
//| lambda     : Noncentrality parameter                             |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityNoncentralBeta(const double x,const double a,const double b,const double lambda,const bool log_mode,int &error_code)
  {
//--- if lambda==0, return Beta density
   if(lambda==0.0)
      return MathProbabilityDensityBeta(x,a,b,error_code);
//--- check parameters
   if(!MathIsValidNumber(x) || !MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- a,b,lambda must be positive
   if(a<=0.0 || b<=0.0 || lambda<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(x<=0.0 || x>=1.0)
      return TailLog0(true,log_mode);
//--- factors
   double lambda_half=lambda*0.5;
   double fact_mult=1.0;
   double pwr_lambda_half=1.0;
   double pwr_x=MathExp((a-1.0)*MathLog(x));
   double r_beta=MathBeta(a,b);
   double pdf=0;
//--- direct sum calculation  
   for(int j=0;; j++)
     {
      if(j>0)
        {
         pwr_x*=x;
         pwr_lambda_half*=lambda_half;
         fact_mult/=j;
         double jm1=j-1;
         r_beta*=((a+jm1)/(a+b+jm1));
        }
      double term=pwr_x*fact_mult*pwr_lambda_half/r_beta;
      //---
      if(term<10E-18)
         break;
      pdf+=term;
     }
//--- calculate density coef
   pdf*=MathExp((b-1.0)*MathLog(1.0-x))*MathExp(-lambda_half);
//--- return density
   return TailLogValue(pdf,true,log_mode);
  }
//+------------------------------------------------------------------+
//| Noncental Beta density function (PDF)                            |
//+------------------------------------------------------------------+
//| The function returns the probability density function            |
//| of the Noncental Beta distribution with parameters a,b,lambda.   |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Random variable                                     |
//| a          : First shape parameter                               |
//| b          : Second shape parameter                              |
//| lambda     : Noncentrality parameter                             |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The probability density evaluated at x.                          |
//+------------------------------------------------------------------+
double MathProbabilityDensityNoncentralBeta(const double x,const double a,const double b,const double lambda,int &error_code)
  {
   return MathProbabilityDensityNoncentralBeta(x,a,b,lambda,false,error_code);
  }
//+------------------------------------------------------------------+
//| Noncental Beta density function (PDF)                            |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Noncentral Beta distribution with parameters a,b,lambda      |
//| for values in x.                                                 |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : First shape parameter                               |
//| b          : Second shape parameter                              |
//| lambda     : Noncentrality parameter                             |
//| tail       : Flag to calculate lower tail                        |
//| log_mode   : Logarithm mode flag, if true it returns Log values  |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityNoncentralBeta(const double &x[],const double a,const double b,const double lambda,const bool log_mode,double &result[])
  {
//--- if lambda==0, return Beta density
   if(lambda==0.0)
      return MathProbabilityDensityBeta(x,a,b,log_mode,result);
//--- check parameters
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- a,b,lambda must be positive
   if(a<=0.0 || b<=0.0 || lambda<0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;
//--- common factors
   double lambda_half=lambda*0.5;
   double exp_lambda_half=MathExp(-lambda_half);
   double r_beta0=MathBeta(a,b);
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(x_arg<=0.0 || x_arg>=1.0)
         result[i]=TailLog0(true,log_mode);
      else
        {
         double fact_mult=1.0;
         double pwr_lambda_half=1.0;
         double pwr_x=MathExp((a-1.0)*MathLog(x_arg));
         double r_beta=r_beta0;
         double pdf=0;
         for(int j=0;; j++)
           {
            if(j>0)
              {
               pwr_x*=x_arg;
               pwr_lambda_half*=lambda_half;
               fact_mult/=j;
               double jm1=j-1;
               r_beta*=((a+jm1)/(a+b+jm1));
              }
            double term=pwr_x*fact_mult*pwr_lambda_half/r_beta;
            //---
            if(term<10E-18)
               break;
            pdf+=term;
           }
         //--- calculate density coef
         pdf*=MathExp((b-1.0)*MathLog(1.0-x_arg))*exp_lambda_half;
         result[i]=TailLogValue(pdf,true,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncental Beta density function (PDF)                            |
//+------------------------------------------------------------------+
//| The function calculates the probability density function of      |
//| the Noncentral Beta distribution with parameters a,b,lambda      |
//| for values in x.                                                 |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : First shape parameter                               |
//| b          : Second shape parameter                              |
//| lambda     : Noncentrality parameter                             |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+  
bool MathProbabilityDensityNoncentralBeta(const double &x[],const double a,const double b,const double lambda,double &result[])
  {
   return MathProbabilityDensityNoncentralBeta(x,a,b,lambda,false,result);
  }
//+------------------------------------------------------------------+
//| Noncental Beta cumulative distribution function (CDF)            |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Noncental Beta distribution with parameters a,b,lambda  |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Input parameters:                                                |
//| x          : The desired quantile                                |
//| a          : First shape parameter                               |
//| b          : Second shape parameter                              |
//| lambda     : Noncentrality parameter                             |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Noncental Beta cumulative distribution function |
//| with parameters a,b,lambda, evaluated at x.                      |
//|                                                                  |
//|                          Infinity                                |
//|            F(x,a,b,lambda)=Sum p(k)*Ix(a+k,b)                    |
//|                            k=0                                   |
//|                                                                  |
//| where p(k)=(1/k!)*exp(-lambda/2)*(lambda/2)^k,                   |
//|       Ix(a,b) - incomplete Beta function                         |
//|                                                                  |
//| Author: John Burkardt                                            | 
//|                                                                  |
//| Reference:                                                       |
//| Harry Posten,"An Effective Algorithm for the Noncentral Beta     |
//| Distribution Function", The American Statistician,               |
//| Volume 47, Number 2, May 1993, pages 129-131.                    |
//+------------------------------------------------------------------+
double MathCumulativeDistributionNoncentralBeta(const double x,const double a,const double b,const double lambda,const bool tail,const bool log_mode,int &error_code)
  {
//--- if lambda==0, return Beta CDF
   if(lambda==0.0)
      return MathCumulativeDistributionBeta(x,a,b,error_code);
//--- check parameters
   if(!MathIsValidNumber(x) || !MathIsValidNumber(a) || !MathIsValidNumber(b))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- a,b,lambda must be positive
   if(a<=0.0 || b<=0.0 || lambda<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
   if(x<=0.0)
      return TailLog0(tail,log_mode);
   if(x>=1.0)
      return TailLog1(tail,log_mode);

   const int max_terms=100;
   double c=lambda*0.5;
   double x0 = int(MathMax(c - 5*MathSqrt(c), 0));
   double a0 = a + x0;
   double beta = MathGammaLog(a0) + MathGammaLog(b) - MathGammaLog(a0+b);
   double temp = MathBetaIncomplete(x, a0, b);
   double gx=MathExp(a0*MathLog(x)+b*MathLog(1-x)-beta-MathLog(a0));

   double q=0;
   if(a0>a)
      q=MathExp(-c+x0*MathLog(c)-MathGammaLog(x0+1));
   else
      q=MathExp(-c);

   double sumq=1-q;
   double betanc=q*temp;
   double ab=a+b;
   int j=0;
   for(;;)
     {
      j++;
      temp-=gx;
      gx*=x*(ab+j-1)/(a+j);
      q*=c/j;
      sumq-=q;
      betanc+=temp*q;
      double err=(temp-gx)*sumq;
      if(j>max_terms || err<1E-18)
         break;
     }
   double cdf=MathMin(betanc,1.0);
   return TailLogValue(cdf,tail,log_mode);
  }
//+------------------------------------------------------------------+
//| Noncental Beta cumulative distribution function (CDF)            |
//+------------------------------------------------------------------+
//| The function returns the probability that an observation         |
//| from the Noncental Beta distribution with parameters a,b,lambda  |
//| is less than or equal to x.                                      |
//|                                                                  |
//| Input parameters:                                                |
//| x          : The desired quantile                                |
//| a          : First shape parameter                               |
//| b          : Second shape parameter                              |
//| lambda     : Noncentrality parameter                             |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| The value of the Noncental Beta cumulative distribution function |
//| with parameters a,b,lambda, evaluated at x.                      |
//+------------------------------------------------------------------+
double MathCumulativeDistributionNoncentralBeta(const double x,const double a,const double b,const double lambda,int &error_code)
  {
   return MathCumulativeDistributionNoncentralBeta(x,a,b,lambda,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Noncental Beta cumulative distribution function (CDF)            |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Noncentral Beta distribution with parameters a,b,lambda      |
//| for values in x.                                                 |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : First shape parameter                               |
//| b          : Second shape parameter                              |
//| lambda     : Noncentrality parameter                             |
//| log_mode   : Logarithm mode, if true it calculates Log values    |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionNoncentralBeta(const double &x[],const double a,const double b,const double lambda,const bool tail,const bool log_mode,double &result[])
  {
//--- if lambda==0, return Beta CDF
   if(lambda==0.0)
      return MathCumulativeDistributionBeta(x,a,b,tail,log_mode,result);
//--- check parameters
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b))
      return false;
//--- a,b,lambda must be positive
   if(a<=0.0 || b<=0.0 || lambda<0)
      return false;

   int data_count=ArraySize(x);
   if(data_count==0)
      return false;

   const int max_terms=100;
   ArrayResize(result,data_count);
   for(int i=0; i<data_count; i++)
     {
      double x_arg=x[i];

      if(x_arg<=0.0)
         result[i]=TailLog0(tail,log_mode);
      if(x_arg>=1.0)
         result[i]=TailLog1(tail,log_mode);
      else
        {
         double c=lambda*0.5;
         double x0 = int(MathMax(c - 5*MathSqrt(c), 0));
         double a0 = a + x0;
         double beta = MathGammaLog(a0) + MathGammaLog(b) - MathGammaLog(a0+b);
         double temp = MathBetaIncomplete(x_arg, a0, b);
         double gx=MathExp(a0*MathLog(x_arg)+b*MathLog(1-x_arg)-beta-MathLog(a0));

         double q=0;
         if(a0>a)
            q=MathExp(-c+x0*MathLog(c)-MathGammaLog(x0+1));
         else
            q=MathExp(-c);

         double sumq=1-q;
         double betanc=q*temp;
         int j=0;
         double ab=a+b;
         for(;;)
           {
            j++;
            temp-=gx;
            gx*=x_arg*(ab+j-1)/(a+j);
            q*=c/j;
            sumq-=q;
            betanc+=temp*q;
            double err=(temp-gx)*sumq;
            if(j>max_terms || err<1E-18)
               break;
           }
         double cdf=MathMin(betanc,1.0);
         result[i]=TailLogValue(cdf,tail,log_mode);
        }
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncental Beta cumulative distribution function (CDF)            |
//+------------------------------------------------------------------+
//| The function calculates the cumulative distribution function of  |
//| the Noncentral Beta distribution with parameters a,b,lambda      |
//| for values in x[] array.                                         |
//|                                                                  |
//| Arguments:                                                       |
//| x          : Array with random variables                         |
//| a          : First shape parameter                               |
//| b          : Second shape parameter                              |
//| lambda     : Noncentrality parameter                             |
//| result     : Array with calculated values                        |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathCumulativeDistributionNoncentralBeta(const double &x[],const double a,const double b,const double lambda,double &result[])
  {
   return MathCumulativeDistributionNoncentralBeta(x,a,b,lambda,true,false,result);
  }
//+------------------------------------------------------------------+
//| Noncental Beta distribution quantile function (inverse CDF)      |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Noncental Beta distribution with parameters a,b  |
//| and lambda for the desired probability.                          |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| a           : First shape parameter                              |
//| b           : Second shape parameter                             |
//| lambda      : Noncentrality parameter                            |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function of     |
//| of Noncental Beta distribution with parameters a,b and lambda.   |
//+------------------------------------------------------------------+
double MathQuantileNoncentralBeta(const double probability,const double a,const double b,const double lambda,const bool tail,const bool log_mode,int &error_code)
  {
   if(log_mode==true && probability==QNEGINF)
      return 0.0;
   if(log_mode==false && probability==0)
      return 0.0;
//--- if lambda==0, return beta quantile
   if(lambda==0.0)
      return MathQuantileBeta(probability,a,b,error_code);
//--- check parameters
   if(!MathIsValidNumber(probability) || !MathIsValidNumber(a) || !MathIsValidNumber(b) || !MathIsValidNumber(lambda))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- a,b,lambda must be positive
   if(a<=0.0 || b<=0.0 || lambda<0)
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
//--- check probabilty
   if(prob==0.0)
      return 0.0;
   if(prob==1.0)
      return 1.0;

   double lambda_half=lambda*0.5;
   double lambda_half_log=MathLog(lambda_half);
   double lambda_half_sqrt=MathSqrt(lambda_half);
   double lambda_half_exp=MathExp(-lambda_half);

   double x0=int(MathMax(lambda_half-5*lambda_half_sqrt,0));
   double b_gamma_log=MathGammaLog(b);
   double eps=10E-18;
   double h_min=MathSqrt(eps);

//double lambda_half=lambda*0.5;
   double r_beta0=MathBeta(a,b);

   int    err_code=0;
   double x=0.5;
   double h=1.0;
   const int max_terms=100;
//--- Newton iterations
   const int max_iterations=50;
   int iterations=0;
   while(iterations<max_iterations)
     {
      //--- check convergence
      if((MathAbs(h)>h_min*MathAbs(x) && MathAbs(h)>h_min)==false)
         break;

      //--- calculate PDF
      double pdf=0;
      if(x<=0.0 || x>=1.0)
         pdf=0;
      else
        {
         double fact_mult=1.0;
         double pwr_lambda_half=1.0;
         double pwr_x=MathExp((a-1.0)*MathLog(x));
         double r_beta=r_beta0;
         //--- direct sum calculation  
         for(int j=0;; j++)
           {
            if(j>0)
              {
               pwr_x*=x;
               pwr_lambda_half*=lambda_half;
               fact_mult/=j;
               double jm1=j-1;
               r_beta*=((a+jm1)/(a+b+jm1));
              }
            double term=pwr_x*fact_mult*pwr_lambda_half/r_beta;
            //---
            if(term<10E-18)
               break;
            pdf+=term;
           }
         //--- calculate density coef
         pdf*=MathExp((b-1.0)*MathLog(1.0-x))*lambda_half_exp;
        }

      //--- calculate CDF
      double cdf=0;
      if(x<=0.0)
         cdf=0;
      if(x>=1.0)
         cdf=1;
      else
        {
         double a0=a+x0;
         double beta = MathGammaLog(a0) + b_gamma_log - MathGammaLog(a0+b);
         double temp = MathBetaIncomplete(x, a0, b);
         double gx=MathExp(a0*MathLog(x)+b*MathLog(1-x)-beta-MathLog(a0));

         double q=0;
         if(a0>a)
            q=MathExp(-lambda_half+x0*lambda_half_log-MathGammaLog(x0+1));
         else
            q=lambda_half_exp;

         double sumq=1-q;
         double betanc=q*temp;
         int j=0;
         double ab=a+b;
         for(;;)
           {
            j++;
            temp-=gx;
            gx*=x*(ab+j-1)/(a+j);
            q*=lambda_half/j;
            sumq-=q;
            betanc+=temp*q;
            double err=(temp-gx)*sumq;
            if(j>max_terms || err<1E-18)
               break;
           }
         cdf=MathMin(betanc,1.0);
        }

      //--- calculate ratio
      h=(cdf-prob)/pdf;

      double x_new=x-h;
      if(x_new<0.0)
         x_new=x*0.1;
      else
      if(x_new>1.0)
         x_new=1.0-(1-x)*0.1;

      if(MathAbs(x_new-x)<10E-16)
         break;
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
   return x;
  }
//+------------------------------------------------------------------+
//| Noncental Beta distribution quantile function (inverse CDF)      |
//+------------------------------------------------------------------+
//| The function returns the inverse cumulative distribution         |
//| function of the Noncental Beta distribution with parameters a, b |
//| and lambda for the desired probability.                          |
//|                                                                  |
//| Arguments:                                                       |
//| probability : The desired probability                            |
//| a           : First shape parameter                              |
//| b           : Second shape parameter                             |
//| lambda      : Noncentrality parameter                            |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The value of the inverse cumulative distribution function of     |
//| of Noncental Beta distribution with parameters a,b and lambda.   |
//+------------------------------------------------------------------+
double MathQuantileNoncentralBeta(const double probability,const double a,const double b,const double lambda,int &error_code)
  {
   return MathQuantileNoncentralBeta(probability,a,b,lambda,true,false,error_code);
  }
//+------------------------------------------------------------------+
//| Noncental Beta distribution quantile function (inverse CDF)      |
//+------------------------------------------------------------------+
//| The function calculates  the inverse cumulative distribution     |
//| function of the Noncentral Beta distribution with parameter a,b  |
//| lambda for the probability values from array.                    |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| a           : First shape parameter                              |
//| b           : Second shape parameter                             |
//| lambda      : Noncentrality parameter                            |
//| tail        : Flag to calculate lower tail                       |
//| log_mode    : Logarithm mode, if true it calculates Log values   |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileNoncentralBeta(const double &probability[],const double a,const double b,const double lambda,const bool tail,const bool log_mode,double &result[])
  {
//--- if lambda==0, return beta quantile
   if(lambda==0.0)
      return MathQuantileBeta(probability,a,b,tail,log_mode,result);
//--- check parameters
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b) || !MathIsValidNumber(lambda))
      return false;
//--- a,b,lambda must be positive
   if(a<=0.0 || b<=0.0 || lambda<0)
      return false;

   int data_count=ArraySize(probability);
   if(data_count==0)
      return false;

   int err_code=0;
   ArrayResize(result,data_count);

   double lambda_half=lambda*0.5;
   double lambda_half_log=MathLog(lambda_half);
   double lambda_half_sqrt=MathSqrt(lambda_half);
   double lambda_half_exp=MathExp(-lambda_half);
   double r_beta0=MathBeta(a,b);

   double x0=int(MathMax(lambda_half-5*lambda_half_sqrt,0));
   double b_gamma_log=MathGammaLog(b);
   const double eps=10E-18;
   double h_min=MathSqrt(eps);
   const int max_terms=100;

   for(int i=0; i<data_count; i++)
     {
      //--- calculate real probability
      double prob=TailLogProbability(probability[i],tail,log_mode);

      if(!MathIsValidNumber(prob))
         return false;

      if(prob==0.0)
         result[i]=0.0;
      else
      if(prob==1.0)
         result[i]=1.0;
      else
        {
         double x=0.5;
         double h=1.0;
         //--- Newton iterations
         const int max_iterations=50;
         int iterations=0;
         while(iterations<max_iterations)
           {
            //--- check convergence
            if((MathAbs(h)>h_min*MathAbs(x) && MathAbs(h)>h_min)==false)
               break;

            //--- calculate PDF
            double pdf=0;
            if(x<=0.0 || x>=1.0)
               pdf=0;
            else
              {
               double fact_mult=1.0;
               double pwr_lambda_half=1.0;
               double pwr_x=MathExp((a-1.0)*MathLog(x));
               double r_beta=r_beta0;
               //--- direct sum calculation  
               for(int j=0;; j++)
                 {
                  if(j>0)
                    {
                     pwr_x*=x;
                     pwr_lambda_half*=lambda_half;
                     fact_mult/=j;
                     double jm1=j-1;
                     r_beta*=((a+jm1)/(a+b+jm1));
                    }
                  double term=pwr_x*fact_mult*pwr_lambda_half/r_beta;
                  //---
                  if(term<10E-18)
                     break;
                  pdf+=term;
                 }
               //--- calculate density coef
               pdf*=MathExp((b-1.0)*MathLog(1.0-x))*lambda_half_exp;
              }

            //--- calculate CDF
            double cdf=0;
            if(x<=0.0)
               cdf=0;
            if(x>=1.0)
               cdf=1;
            else
              {
               double a0=a+x0;
               double beta = MathGammaLog(a0) + b_gamma_log - MathGammaLog(a0+b);
               double temp = MathBetaIncomplete(x, a0, b);
               double gx=MathExp(a0*MathLog(x)+b*MathLog(1-x)-beta-MathLog(a0));

               double q=0;
               if(a0>a)
                  q=MathExp(-lambda_half+x0*lambda_half_log-MathGammaLog(x0+1));
               else
                  q=lambda_half_exp;

               double sumq=1-q;
               double betanc=q*temp;
               int j=0;
               double ab=a+b;
               for(;;)
                 {
                  j++;
                  temp-=gx;
                  gx*=x*(ab+j-1)/(a+j);
                  q*=lambda_half/j;
                  sumq-=q;
                  betanc+=temp*q;
                  double err=(temp-gx)*sumq;
                  if(j>max_terms || err<1E-18)
                     break;
                 }
               cdf=MathMin(betanc,1.0);
              }

            //--- calculate ratio
            h=(cdf-prob)/pdf;

            double x_new=x-h;
            if(x_new<0.0)
               x_new=x*0.1;
            else
            if(x_new>1.0)
               x_new=1.0-(1-x)*0.1;

            if(MathAbs(x_new-x)<10E-16)
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
   return true;
  }
//+------------------------------------------------------------------+
//| Noncental Beta distribution quantile function (inverse CDF)      |
//+------------------------------------------------------------------+
//| The function calculates  the inverse cumulative distribution     |
//| function of the Noncentral Beta distribution with parameter a,b  |
//| lambda for the probability values from array.                    |
//|                                                                  |
//| Arguments:                                                       |
//| probability : Array with probabilities                           |
//| a           : First shape parameter                              |
//| b           : Second shape parameter                             |
//| lambda      : Noncentrality parameter                            |
//| result      : Array with calculated values                       |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathQuantileNoncentralBeta(const double &probability[],const double a,const double b,const double lambda,double &result[])
  {
   return MathQuantileNoncentralBeta(probability,a,b,lambda,true,false,result);
  }
//+------------------------------------------------------------------+
//| Random variate from the Noncentral Beta distribution             |
//+------------------------------------------------------------------+
//| Compute the random variable from the Noncentral Beta             |
//| distribution with parameters a,b and lambda.                     |
//|                                                                  |
//| Arguments:                                                       |
//| a           : First shape parameter                              |
//| b           : Second shape parameter                             |
//| lambda      : Noncentrality parameter                            |
//| error_code  : Variable for error code                            |
//|                                                                  |
//| Return value:                                                    |
//| The random value with Noncentral Beta distribution.              |
//+------------------------------------------------------------------+
double MathRandomNoncentralBeta(const double a,const double b,const double lambda,int &error_code)
  {
//--- if lambda==0, return beta random variate
   if(lambda==0.0)
      return MathRandomBeta(a,b,error_code);
//--- check parameters
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b) || !MathIsValidNumber(lambda))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return QNaN;
     }
//--- a,b,lambda must be positive
   if(a<=0.0 || b<=0.0 || lambda<0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return QNaN;
     }

   error_code=ERR_OK;
//--- generate random number using Noncentral ChiSquare
   double chi1=MathRandomNoncentralChiSquare(2*a,2*lambda,error_code);
   double chi2=MathRandomNoncentralChiSquare(2*b,2*lambda,error_code);
   return chi1/(chi1+chi2);
  }
//+------------------------------------------------------------------+
//| Random variate from the Noncentral Beta distribution              |
//+------------------------------------------------------------------+
//| Generates random variables from the Noncentral Beta distribution |
//| with parameters a,b, lambda.                                     |
//|                                                                  |
//| Arguments:                                                       |
//| a          : First shape parameter                               |
//| b          : Second shape parameter                              |
//| lambda     : Noncentrality parameter                             |
//| data_count : Number of values needed                             |
//| result     : Output array with random values                     |
//|                                                                  |
//| Return value:                                                    |
//| true if successful, otherwise false.                             |
//+------------------------------------------------------------------+
bool MathRandomNoncentralBeta(const double a,const double b,const double lambda,const int data_count,double &result[])
  {
//--- if lambda==0, return beta random variate
   if(lambda==0.0)
      return MathRandomBeta(a,b,data_count,result);
//--- check parameters
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b) || !MathIsValidNumber(lambda))
      return false;
//--- a,b,lambda must be positive
   if(a<=0.0 || b<=0.0 || lambda<0)
      return false;

   int error_code=0;
//--- prepare output array and calculate random values
   ArrayResize(result,data_count);
   double lambda2=lambda*2;
   double a2=a*2;
   double b2=b*2;
   for(int i=0; i<data_count; i++)
     {
      //--- generate random number using Noncentral ChiSquare
      double chi1=MathRandomNoncentralChiSquare(a2,lambda2,error_code);
      double chi2=MathRandomNoncentralChiSquare(b2,lambda2,error_code);
      result[i]=chi1/(chi1+chi2);
     }
   return true;
  }
//+------------------------------------------------------------------+
//| Noncental Beta distribution moments                              |
//+------------------------------------------------------------------+
//| The function calculates 4 first moments of the Noncental Beta    |
//| distribution with parameters a,b and lambda.                     |
//|                                                                  |
//| Arguments:                                                       |
//| a          : First shape parameter                               |
//| b          : Second shape parameter                              |
//| lambda     : Noncentrality parameter                             |
//| mean       : Variable for mean value (1st moment)                |
//| variance   : Variable for variance value (2nd moment)            |
//| skewness   : Variable for skewness value (3rd moment)            |
//| kurtosis   : Variable for kurtosis value (4th moment)            |
//| error_code : Variable for error code                             |
//|                                                                  |
//| Return value:                                                    |
//| true if moments calculated successfully, otherwise false.        |
//+------------------------------------------------------------------+
double MathMomentsNoncentralBeta(const double a,const double b,const double lambda,double &mean,double &variance,double &skewness,double &kurtosis,int &error_code)
  {
//--- default values
   mean    =QNaN;
   variance=QNaN;
   skewness=QNaN;
   kurtosis=QNaN;
//--- check NaN
   if(!MathIsValidNumber(a) || !MathIsValidNumber(b) || !MathIsValidNumber(lambda))
     {
      error_code=ERR_ARGUMENTS_NAN;
      return false;
     }
//--- a and b must be positive
   if(a<=0.0 || b<=0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }
//--- check lambda
   if(lambda<0.0)
     {
      error_code=ERR_ARGUMENTS_INVALID;
      return false;
     }

   error_code=ERR_OK;
//--- prepare coefficients
   double lambda_half=lambda*0.5;
//--- hypergeometric function values   
   double f1=MathHypergeometric2F2(a+1,a+b,a,a+b+1,lambda_half);
   double f2=MathHypergeometric2F2(a+2,a+b,a,a+b+2,lambda_half);
   double f3=MathHypergeometric2F2(a+3,a+b,a,a+b+3,lambda_half);
   double f4=MathHypergeometric2F2(a+4,a+b,a,a+b+4,lambda_half);
//--- exponents
   double exp_lambda_half=MathExp(-lambda_half);
   double exp_lambda=MathPow(exp_lambda_half,2);
//--- factors
   double aab=a/(a+b);
   double aab2=MathPow(aab,2);
   double ab1=(a+1)/(a+b+1);
   double ab2=(a+2)/(a+b+2);
   double ab3=(a+3)/(a+b+3);
//--- calculate moments
   mean=aab*exp_lambda_half*f1;
   double mean2=MathPow(mean,2);
   variance=aab*ab1*exp_lambda_half*f2-mean2;
   skewness=(2*MathPow(mean,3)+exp_lambda_half*aab*ab1*(-3*mean*f2+ab2*f3))*MathPow(variance,-1.5);
   kurtosis=-3+(-3*MathPow(mean,4)+exp_lambda*f1*aab2*(6*mean*ab1*f2-4*ab1*ab2*f3)+aab*ab1*ab2*ab3*exp_lambda_half*f4)*MathPow(aab*ab1*exp_lambda_half*f2-mean2,-2);
//--- successful
   return true;
  }
//+------------------------------------------------------------------+
