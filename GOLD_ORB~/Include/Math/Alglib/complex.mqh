//+------------------------------------------------------------------+
//|                                                      complex.mqh |
//|            Copyright 2003-2012 Sergey Bochkanov (ALGLIB project) |
//|                   Copyright 2012-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
//| Implementation of ALGLIB library in MetaQuotes Language 5        |
//|                                                                  |
//| The features of the library include:                             |
//| - Linear algebra (direct algorithms, EVD, SVD)                   |
//| - Solving systems of linear and non-linear equations             |
//| - Interpolation                                                  |
//| - Optimization                                                   |
//| - FFT (Fast Fourier Transform)                                   |
//| - Numerical integration                                          |
//| - Linear and nonlinear least-squares fitting                     |
//| - Ordinary differential equations                                |
//| - Computation of special functions                               |
//| - Descriptive statistics and hypothesis testing                  |
//| - Data analysis - classification, regression                     |
//| - Implementing linear algebra algorithms, interpolation, etc.    |
//|   in high-precision arithmetic (using MPFR)                      |
//|                                                                  |
//| This file is free software; you can redistribute it and/or       |
//| modify it under the terms of the GNU General Public License as   |
//| published by the Free Software Foundation (www.fsf.org); either  |
//| version 2 of the License, or (at your option) any later version. |
//|                                                                  |
//| This program is distributed in the hope that it will be useful,  |
//| but WITHOUT ANY WARRANTY; without even the implied warranty of   |
//| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the     |
//| GNU General Public License for more details.                     |
//+------------------------------------------------------------------+
//+------------------------------------------------------------------+
//| Complex numbers                                                  |
//+------------------------------------------------------------------+
struct al_complex
  {
public:
   double            re; // real part
   double            im; // imaginary part

public:
                     al_complex(void);
                     al_complex(const double x);
                     al_complex(const double x,const double y);
                    ~al_complex(void);
   //--- operations
   void              Copy(const al_complex &rhs);
   bool              Eq(const al_complex &lhs,const al_complex &rhs);
   bool              NotEq(const al_complex &lhs,const al_complex &rhs);
   al_complex        Add(const al_complex &lhs,const al_complex &rhs);
   al_complex        Sub(const al_complex &lhs,const al_complex &rhs);
   al_complex        Mul(const al_complex &lhs,const al_complex &rhs);
   al_complex        Div(const al_complex &lhs,const al_complex &rhs);
   al_complex        Conjugate(void);
   //--- overloading
   void              operator=(const double rhs);
   void              operator=(const al_complex &rhs);
   void              operator+=(const al_complex &rhs);
   void              operator-=(const al_complex &rhs);
   bool              operator==(const al_complex &rhs);
   bool              operator==(const double rhs);
   bool              operator!=(const al_complex &rhs);
   bool              operator!=(const double rhs);
   al_complex        operator+(const al_complex &rhs);
   al_complex        operator+(const double rhs);
   al_complex        operator+(void);
   al_complex        operator-(const al_complex &rhs);
   al_complex        operator-(const double rhs);
   al_complex        operator-(void);
   al_complex        operator*(const al_complex &rhs);
   al_complex        operator*(const double rhs);
   al_complex        operator/(const al_complex &rhs);
   al_complex        operator/(const double rhs);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
al_complex::al_complex(void): re(0),im(0)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with one parameter                                   |
//+------------------------------------------------------------------+
al_complex::al_complex(const double x): re(x),im(0)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with two parameters                                  |
//+------------------------------------------------------------------+
al_complex::al_complex(const double x,const double y): re(x),im(y)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
al_complex::~al_complex(void)
  {

  }
//+------------------------------------------------------------------+
//| Copy complex                                                     |
//+------------------------------------------------------------------+
void al_complex::Copy(const al_complex &rhs)
  {
   re=rhs.re;
   im=rhs.im;
  }
//+------------------------------------------------------------------+
//| Comparison (==)                                                  |
//+------------------------------------------------------------------+
bool al_complex::Eq(const al_complex &lhs,const al_complex &rhs)
  {
//--- comparison
   if(lhs.re==rhs.re && lhs.im==rhs.im)
      return(true);
//--- numbers are not equal
   return(false);
  }
//+------------------------------------------------------------------+
//| Comparison (!=)                                                  |
//+------------------------------------------------------------------+
bool al_complex::NotEq(const al_complex &lhs,const al_complex &rhs)
  {
//--- comparison
   if(lhs.re!=rhs.re || lhs.im!=rhs.im)
      return(true);
//--- numbers are equal
   return(false);
  }
//+------------------------------------------------------------------+
//| Sum                                                              |
//+------------------------------------------------------------------+
al_complex al_complex::Add(const al_complex &lhs,const al_complex &rhs)
  {
   al_complex res;
//--- sum
   res.re=lhs.re+rhs.re;
   res.im=lhs.im+rhs.im;
//--- return result
   return(res);
  }
//+------------------------------------------------------------------+
//| Subtraction                                                      |
//+------------------------------------------------------------------+
al_complex al_complex::Sub(const al_complex &lhs,const al_complex &rhs)
  {
   al_complex res;
//--- subtraction
   res.re=lhs.re-rhs.re;
   res.im=lhs.im-rhs.im;
//--- return result
   return(res);
  }
//+------------------------------------------------------------------+
//| Multiplication                                                   |
//+------------------------------------------------------------------+
al_complex al_complex::Mul(const al_complex &lhs,const al_complex &rhs)
  {
   al_complex res;
//--- multiplication
   res.re=lhs.re*rhs.re-lhs.im*rhs.im;
   res.im=lhs.re*rhs.im+lhs.im*rhs.re;
//--- return result
   return(res);
  }
//+------------------------------------------------------------------+
//| Division                                                          |
//+------------------------------------------------------------------+
al_complex al_complex::Div(const al_complex &lhs,const al_complex &rhs)
  {
//--- empty complex value
   al_complex res(EMPTY_VALUE,EMPTY_VALUE);
//--- check
   if(rhs.re==0 && rhs.im==0)
     {
      Print(__FUNCTION__+": number is zero");
      return(res);
     }
//--- create variables
   double e;
   double f;
//--- division
   if(MathAbs(rhs.im)<MathAbs(rhs.re))
     {
      e=rhs.im/rhs.re;
      f=rhs.re+rhs.im*e;
      res.re=(lhs.re+lhs.im*e)/f;
      res.im=(lhs.im-lhs.re*e)/f;
      //--- return result
      return(res);
     }
   e=rhs.re/rhs.im;
   f=rhs.im+rhs.re*e;
   res.re=(lhs.im+lhs.re*e)/f;
   res.im=(-lhs.re+lhs.im*e)/f;
//--- return result
   return(res);
  }
//+------------------------------------------------------------------+
//| Conjugate                                                        |
//+------------------------------------------------------------------+
al_complex al_complex::Conjugate(void)
  {
//--- conjugate
   al_complex res(re,-im);
//--- return result
   return res;
  }
//+------------------------------------------------------------------+
//| Overloading (=)                                                  |
//+------------------------------------------------------------------+
void al_complex::operator=(const double rhs)
  {
   re=rhs;
   im=0;
  }
//+------------------------------------------------------------------+
//| Overloading (=)                                                  |
//+------------------------------------------------------------------+
void al_complex::operator=(const al_complex &rhs)
  {
   this.Copy(rhs);
  }
//+------------------------------------------------------------------+
//| Overloading (+=)                                                 |
//+------------------------------------------------------------------+
void al_complex::operator+=(const al_complex &rhs)
  {
   re+=rhs.re;
   im+=rhs.im;
  }
//+------------------------------------------------------------------+
//| Overloading (-=)                                                 |
//+------------------------------------------------------------------+
void al_complex::operator-=(const al_complex &rhs)
  {
   re-=rhs.re;
   im-=rhs.im;
  }
//+------------------------------------------------------------------+
//| Overloading (==)                                                 |
//+------------------------------------------------------------------+
bool al_complex::operator==(const al_complex &rhs)
  {
   return(Eq(this,rhs));
  }
//+------------------------------------------------------------------+
//| Overloading (==)                                                 |
//+------------------------------------------------------------------+
bool al_complex::operator==(const double rhs)
  {
   al_complex r(rhs,0);
//--- return result
   return(Eq(this,r));
  }
//+------------------------------------------------------------------+
//| Overloading (!=)                                                 |
//+------------------------------------------------------------------+
bool al_complex::operator!=(const al_complex &rhs)
  {
   return(NotEq(this,rhs));
  }
//+------------------------------------------------------------------+
//| Overloading (!=)                                                 |
//+------------------------------------------------------------------+
bool al_complex::operator!=(const double rhs)
  {
   al_complex r(rhs,0);
//--- return result
   return(NotEq(this,r));
  }
//+------------------------------------------------------------------+
//| Overloading of binary (+)                                        |
//+------------------------------------------------------------------+
al_complex al_complex::operator+(const al_complex &rhs)
  {
   return(Add(this,rhs));
  }
//+------------------------------------------------------------------+
//| Overloading of binary (+)                                        |
//+------------------------------------------------------------------+
al_complex al_complex::operator+(const double rhs)
  {
   al_complex r(rhs,0);
//--- return result
   return(Add(this,r));
  }
//+------------------------------------------------------------------+
//| Overloading of unary (+)                                         |
//+------------------------------------------------------------------+
al_complex al_complex::operator+(void)
  {
//--- return result
   return(this);
  }
//+------------------------------------------------------------------+
//| Overloading of binary (-)                                        |
//+------------------------------------------------------------------+
al_complex al_complex::operator-(const al_complex &rhs)
  {
   return(Sub(this,rhs));
  }
//+------------------------------------------------------------------+
//| Overloading of binary (-)                                        |
//+------------------------------------------------------------------+
al_complex al_complex::operator-(const double rhs)
  {
   al_complex r(rhs,0);
//--- return result
   return(Sub(this,r));
  }
//+------------------------------------------------------------------+
//| Overloading of unary (-)                                         |
//+------------------------------------------------------------------+
al_complex al_complex::operator-(void)
  {
   al_complex c(-this.re,-this.im);
//--- return result
   return(c);
  }
//+------------------------------------------------------------------+
//| Overloading (*)                                                  |
//+------------------------------------------------------------------+
al_complex al_complex::operator*(const al_complex &rhs)
  {
   return(Mul(this,rhs));
  }
//+------------------------------------------------------------------+
//| Overloading (*)                                                  |
//+------------------------------------------------------------------+
al_complex al_complex::operator*(const double rhs)
  {
   al_complex r(rhs,0);
//--- return result
   return(Mul(this,r));
  }
//+------------------------------------------------------------------+
//| Overloading (/)                                                  |
//+------------------------------------------------------------------+
al_complex al_complex::operator/(const al_complex &rhs)
  {
   return(Div(this,rhs));
  }
//+------------------------------------------------------------------+
//| Overloading (/)                                                  |
//+------------------------------------------------------------------+
al_complex al_complex::operator/(const double rhs)
  {
   al_complex r(rhs,0);
//--- return result
   return(Div(this,r));
  }
//+------------------------------------------------------------------+
