//+------------------------------------------------------------------+
//|                                           membershipfunction.mqh |
//|                   Copyright 2015-2017, MetaQuotes Software Corp. |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
//| Implementation of Fuzzy library in MetaQuotes Language 5         |
//|                                                                  |
//| The features of the library include:                             |
//| - Create Mamdani fuzzy model                                     |
//| - Create Sugeno fuzzy model                                      |
//| - Normal membership function                                     |
//| - Triangular membership function                                 |
//| - Trapezoidal membership function                                |
//| - Constant membership function                                   |
//| - Defuzzification method of center of gravity (COG)              |
//| - Defuzzification method of bisector of area (BOA)               |
//| - Defuzzification method of mean of maxima (MeOM)                |
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
#include <Arrays\List.mqh>
//+------------------------------------------------------------------+
//| Purpose: creating membership functions.                          |
//+------------------------------------------------------------------+
//+------------------------------------------------------------------+
//| Types of membership functions composition                        |
//+------------------------------------------------------------------+
enum MfCompositionType
  {
   MinMF,                               // Minumum of functions
   MaxMF,                               // Maximum of functions
   ProdMF,                              // Production of functions
   SumMF                                // Sum of functions
  };
//+------------------------------------------------------------------+
//| The base class of all classes of membership functions            |
//+------------------------------------------------------------------+
class IMembershipFunction : public CObject
  {
public:
   //--- method evaluate value of the membership function
   virtual double    GetValue(const double x)=NULL;
  };
//+------------------------------------------------------------------+
//| Gaussian combination membership function                         |
//+------------------------------------------------------------------+
class CNormalCombinationMembershipFunction : public IMembershipFunction
  {
private:
   double            m_b1;            // Parametr b1: coordinate of the minimum membership function
   double            m_sigma1;        // Parametr sigma1: concentration factor of the left path of function 
   double            m_b2;            // Parametr b2: coordinate of the maximum membership function
   double            m_sigma2;        // Parametr sigma2: concentration factor of the rigth path of function

public:
                     CNormalCombinationMembershipFunction(void);
                     CNormalCombinationMembershipFunction(const double b1,const double sigma1,const double b2,const double sigma2);
                    ~CNormalCombinationMembershipFunction(void);
   //--- methods gets or sets the parametrs
   void              B1(const double b1)           { m_b1=b1;          }
   double            B1(void)                      { return(m_b1);     }
   void              Sigma1(const double sigma1)   { m_sigma1=sigma1;  }
   double            Sigma1(void)                  { return(m_sigma1); }
   void              B2(const double b2)           { m_b2=b2;          }
   double            B2(void)                      { return(m_b2);     }
   void              Sigma2(const double sigma2)   { m_sigma2=sigma2;  }
   double            Sigma2(void)                  { return(m_sigma2); }
   //--- method gets the argument (x axis value)
   double            GetValue(const double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CNormalCombinationMembershipFunction::CNormalCombinationMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+
CNormalCombinationMembershipFunction::CNormalCombinationMembershipFunction(const double b1,const double sigma1,const double b2,const double sigma2)
  {
   m_b1=b1;
   m_sigma1=sigma1;
   m_b2=b2;
   m_sigma2=sigma2;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNormalCombinationMembershipFunction::~CNormalCombinationMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Get argument (x axis value)                                      |
//+------------------------------------------------------------------+
double CNormalCombinationMembershipFunction::GetValue(const double x)
  {
   if(m_b1<=m_b2)
     {
      if(x<m_b1)
        {
         //--- return result
         return (exp((x - m_b1) * (x - m_b1) / ( -2.0 * m_sigma1 * m_sigma1)));
        }
      else if(x>m_b2)
        {
         //--- return result
         return (exp((x - m_b2) * (x - m_b2) / ( -2.0 * m_sigma2 * m_sigma2)));
        }
      else
        {
         //--- m_b1 <= x && x <= m_b2
         //--- return result 
         return (1);
        }
     }
   if(m_b1>m_b2)
     {
      if(x<m_b2)
        {
         //--- return result
         return (exp((x - m_b1) * (x - m_b1) / ( -2.0 * m_sigma1 * m_sigma1)));
        }
      else if(x>m_b1)
        {
         //--- return result
         return (exp((x - m_b2) * (x - m_b2) / ( -2.0 * m_sigma2 * m_sigma2)));
        }
      else
        {
         //--- m_b1 <= x && x <= m_b2
         //--- return result 
         return ( exp((x - m_b1) * (x - m_b1) / ( -2.0 * m_sigma1 * m_sigma1)) * exp((x - m_b2) * (x - m_b2) / ( -2.0 * m_sigma2 * m_sigma2)) );
        }
     }
//--- m_b1 == m_b2
//--- return result  
   return (m_b1);
  }
//+------------------------------------------------------------------+
//| Generalized bell-shaped membership function                      |
//+------------------------------------------------------------------+
class CGeneralizedBellShapedMembershipFunction : public IMembershipFunction
  {
private:
   double            m_a;             // Parametr a: the concentration factor of the membership function
   double            m_b;             // Parametr b: coefficients slope of the membership function
   double            m_c;             // Parametr c:  the maximum coordinate of the membership function

public:
                     CGeneralizedBellShapedMembershipFunction(void);
                     CGeneralizedBellShapedMembershipFunction(const double a,const double b,const double c);
                    ~CGeneralizedBellShapedMembershipFunction(void);
   //--- methods gets or sets the parametrs
   void              A(const double a) { m_a=a;       }
   double            A(void)           { return(m_a); }
   void              B(const double b) { m_b=b;       }
   double            B(void)           { return(m_b); }
   void              C(const double c) { m_c=c;       }
   double            C(void)           { return(m_c); }
   //--- method gets the argument (x axis value)
   double            GetValue(const double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CGeneralizedBellShapedMembershipFunction::CGeneralizedBellShapedMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+                     
CGeneralizedBellShapedMembershipFunction::CGeneralizedBellShapedMembershipFunction(const double a,const double b,const double c)
  {
   m_a=a;
   m_b=b;
   m_c=c;
  }
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
CGeneralizedBellShapedMembershipFunction::~CGeneralizedBellShapedMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Get argument (x axis value)                                      |
//+------------------------------------------------------------------+
double CGeneralizedBellShapedMembershipFunction::GetValue(const double x)
  {
//--- return result
   return (1 / (1 + pow( fabs((x - m_a) / m_c) , 2 * m_b )));
  }
//+------------------------------------------------------------------+
//| S-shaped membership function                                     |
//+------------------------------------------------------------------+
class CS_ShapedMembershipFunction : public IMembershipFunction
  {
private:
   double            m_a;             // Parametr a: beginning of the interval increases        
   double            m_b;             // Parametr b: end of the interval increases  

public:
                     CS_ShapedMembershipFunction(void);
                     CS_ShapedMembershipFunction(const double a,const double b);
                    ~CS_ShapedMembershipFunction(void);
   //--- methods gets or sets the parametrs  
   void              A(const double a) { m_a=a;       }
   double            A(void)           { return(m_a); }
   void              B(const double b) { m_b=b;       }
   double            B(void)           { return(m_b); }
   //--- method gets the argument (x axis value)
   double            GetValue(const double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CS_ShapedMembershipFunction::CS_ShapedMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+
CS_ShapedMembershipFunction::CS_ShapedMembershipFunction(const double a,const double b)
  {
   m_a=a;
   m_b=b;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CS_ShapedMembershipFunction::~CS_ShapedMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Get argument (x axis value)                                      |
//+------------------------------------------------------------------+ 
double CS_ShapedMembershipFunction::GetValue(const double x)
  {
   if(x<=m_a)
     {
      //--- return result
      return (0.0);
     }
   else if(m_a<x && x<=(m_a+m_b)/2.0)
     {
      //--- return result
      return (2.0*((x-m_a)/(m_b-m_a))*((x-m_a)/(m_b-m_a)));
     }
   else if((m_a+m_b)/2.0<x && x<m_b)
     {
      //--- return result
      return (1.0 - 2.0*((x - m_b)/(m_b - m_a))*((x - m_b)/(m_b - m_a)));
     }
   else
     {//--- x >= m_b
      //--- return result
      return (1.0);
     }
  }
//+------------------------------------------------------------------+
//| Z-shaped membership function                                     |
//+------------------------------------------------------------------+
class CZ_ShapedMembershipFunction : public IMembershipFunction
  {
private:
   double            m_a;             // Parametr a: beginning of the interval decreasing       
   double            m_b;             // Parametr b: end of the interval decreasing 

public:
                     CZ_ShapedMembershipFunction(void);
                     CZ_ShapedMembershipFunction(const double a,const double b);
                    ~CZ_ShapedMembershipFunction(void);
   //--- methods gets or sets the parametrs  
   void              A(const double a) { m_a=a;       }
   double            A(void)           { return(m_a); }
   void              B(const double b) { m_b=b;       }
   double            B(void)           { return(m_b); }
   //--- method gets the argument (x axis value)
   double            GetValue(const double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CZ_ShapedMembershipFunction::CZ_ShapedMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+
CZ_ShapedMembershipFunction::CZ_ShapedMembershipFunction(const double a,const double b)
  {
   m_a=a;
   m_b=b;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CZ_ShapedMembershipFunction::~CZ_ShapedMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Get argument (x axis value)                                      |
//+------------------------------------------------------------------+
double CZ_ShapedMembershipFunction::GetValue(const double x)
  {
   if(x<=m_a)
     {
      //--- return result
      return (1.0);
     }
   else if(m_a<x && x<=(m_a+m_b)/2.0)
     {
      //--- return result
      return (1.0 - 2.0*((x - m_a)/(m_b - m_a))*((x - m_a)/(m_b - m_a)));
     }
   else if((m_a+m_b)/2.0<=x && x<=m_b)
     {
      //--- return result
      return (2.0*((x-m_b)/(m_b-m_a))*((x-m_b)/(m_b-m_a)));
     }
   else
     {//--- x >= m_b
      //--- return result
      return (0.0);
     }
  }
//+------------------------------------------------------------------+
//| P-shaped membership function                                     |
//+------------------------------------------------------------------+
class CP_ShapedMembershipFunction : public IMembershipFunction
  {
private:
   double            m_a;             // Parametr a: carrier fuzzy set
   double            m_d;             // Parametr d: carrier fuzzy set
   double            m_b;             // Parametr b: the core of a fuzzy set
   double            m_c;             // Parametr c: the core of a fuzzy set

public:
                     CP_ShapedMembershipFunction(void);
                     CP_ShapedMembershipFunction(const double a,const double b,const double c,const double d);
                    ~CP_ShapedMembershipFunction(void);
   //--- methods gets or sets the parametrs  
   void              A(const double a) { m_a=a;       }
   double            A(void)           { return(m_a); }
   void              D(const double d) { m_d=d;       }
   double            D(void)           { return(m_d); }
   void              B(const double b) { m_b=b;       }
   double            B(void)           { return(m_b); }
   void              C(const double c) { m_c=c;       }
   double            C(void)           { return(m_c); }
   //--- method gets the argument (x axis value)
   double            GetValue(const double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CP_ShapedMembershipFunction::CP_ShapedMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+
CP_ShapedMembershipFunction::CP_ShapedMembershipFunction(const double a,const double b,const double c,const double d)
  {
   m_a=a;
   m_d=d;
   m_b=b;
   m_c=c;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CP_ShapedMembershipFunction::~CP_ShapedMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Get argument (x axis value)                                      |
//+------------------------------------------------------------------+     
double CP_ShapedMembershipFunction::GetValue(const double x)
  {
   if(x<=m_a)
     {
      return(0.0);
     }
   else
   if(m_a<x && x<=(m_a+m_b)/2.0)
     {
      return(2.0 *((x-m_a)/(m_b-m_a))*((x-m_a)/(m_b-m_a)));
     }
   else
   if((m_a+m_b)/2.0<x && x<m_b)
     {
      return(1.0-2.0 *((x-m_b)/(m_b-m_a))*((x-m_b)/(m_b-m_a)));
     }
   else
   if(m_b<=x && x<=m_c)
     {
      return(1.0);
     }
   else
   if(m_c<x && x<=(m_c+m_d)/2.0)
     {
      return(1.0-2.0 *((x-m_c)/(m_d-m_c))*((x-m_c)/(m_d-m_c)));
     }
   else
   if((m_c+m_d)/2.0<x && x<m_d)
     {
      return(2.0 *((x-m_d)/(m_d-m_c))*((x-m_d)/(m_d-m_c)));
     }
   else
     {
      return(0.0);
     }
  }
//+------------------------------------------------------------------+
//| Sigmoidal membership function                                    |
//+------------------------------------------------------------------+
class CSigmoidalMembershipFunction : public IMembershipFunction
  {
private:
   double            m_a;            // Parametr a1: the slope coefficient of membership functions
   double            m_c;            // Parametr c1: coordinate of the inflection of membership function

public:
                     CSigmoidalMembershipFunction(void);
                     CSigmoidalMembershipFunction(const double a,const double c);
                    ~CSigmoidalMembershipFunction(void);
   //--- methods gets or sets the parametrs  
   void              A(const double a) { m_a=a;        }
   double            A(void)           { return (m_a); }
   void              C(const double c) { m_c=c;        }
   double            C(void)           { return (m_c); }
   //--- method gets the argument (x axis value)
   double            GetValue(const double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSigmoidalMembershipFunction::CSigmoidalMembershipFunction(void)
  {
   m_a = 0.0;
   m_c = 0.0;
  }
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+
CSigmoidalMembershipFunction::CSigmoidalMembershipFunction(const double a,const double c)
  {
   m_a = a;
   m_c = c;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CSigmoidalMembershipFunction::~CSigmoidalMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Get argument (x axis value)                                      |
//+------------------------------------------------------------------+
double CSigmoidalMembershipFunction::GetValue(const double x)
  {
//--- return result
   return (1/(1+exp(-m_a*(x-m_c))));
  }
//+------------------------------------------------------------------+
//| Product of two sigmoidal membership functions                    |
//+------------------------------------------------------------------+
class CProductTwoSigmoidalMembershipFunctions : public IMembershipFunction
  {
private:
   double            m_a1;            // Parametr a1: the slope coefficient of the first functions
   double            m_c1;            // Parametr c1: coordinate of the inflection of the first function
   double            m_a2;            // Parametr a1: the slope coefficient of the second functions
   double            m_c2;            // Parametr c2: coordinate of the inflection of the second function

public:
                     CProductTwoSigmoidalMembershipFunctions(void);
                     CProductTwoSigmoidalMembershipFunctions(const double a1,const double c1,const double a2,const double c2);
                    ~CProductTwoSigmoidalMembershipFunctions(void);
   //--- methods gets or sets the parametrs  
   void              A1(const double a1)  { m_a1=a1;       }
   double            A1(void)             { return (m_a1); }
   void              C1(const double c1)  { m_c1=c1;       }
   double            C1(void)             { return (m_c1); }
   void              A2(const double a2)  { m_a2=a2;       }
   double            A2(void)             { return (m_a2); }
   void              C2(const double c2)  { m_c2=c2;       }
   double            C2(void)             { return (m_c2); }
   //--- method gets the argument (x axis value)
   double            GetValue(const double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CProductTwoSigmoidalMembershipFunctions::CProductTwoSigmoidalMembershipFunctions(void)
  {
   m_a1 = 0.0;
   m_a2 = 0.0;
   m_c1 = 0.0;
   m_c2 = 0.0;
  }
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+
CProductTwoSigmoidalMembershipFunctions::CProductTwoSigmoidalMembershipFunctions(const double a1,const double c1,const double a2,const double c2)
  {
   m_a1 = a1;
   m_a2 = a2;
   m_c1 = c1;
   m_c2 = c2;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CProductTwoSigmoidalMembershipFunctions::~CProductTwoSigmoidalMembershipFunctions(void)
  {

  }
//+------------------------------------------------------------------+
//| Get argument (x axis value)                                      |
//+------------------------------------------------------------------+
double CProductTwoSigmoidalMembershipFunctions::GetValue(const double x)
  {
   double first_equation=1/(1+exp(-m_a1*(x-m_c1)));
   double second_equation=1/(1+exp(-m_a2*(x-m_c2)));
//--- return result
   return (first_equation * second_equation);
  }
//+------------------------------------------------------------------+
//| Difference between two sigmoidal functions membership function   |
//+------------------------------------------------------------------+
class CDifferencTwoSigmoidalMembershipFunction : public IMembershipFunction
  {
private:
   double            m_a1;            // Parametr a1: the slope coefficient of the first functions
   double            m_c1;            // Parametr c1: coordinate of the inflection of the first function
   double            m_a2;            // Parametr a1: the slope coefficient of the second functions
   double            m_c2;            // Parametr c2: coordinate of the inflection of the second function

public:
                     CDifferencTwoSigmoidalMembershipFunction(void);
                     CDifferencTwoSigmoidalMembershipFunction(const double a1,const double c1,const double a2,const double c2);
                    ~CDifferencTwoSigmoidalMembershipFunction(void);
   //--- methods gets or sets the parametrs  
   void              A1(const double a1)  { m_a1=a1;      }
   double            A1(void)             { return(m_a1); }
   void              C1(const double c1)  { m_c1=c1;      }
   double            C1(void)             { return(m_c1); }
   void              A2(const double a2)  { m_a2=a2;      }
   double            A2(void)             { return(m_a2); }
   void              C2(const double c2)  { m_c2=c2;      }
   double            C2(void)             { return(m_c2); }
   //--- method gets the argument (x axis value)
   double            GetValue(const double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDifferencTwoSigmoidalMembershipFunction::CDifferencTwoSigmoidalMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+
CDifferencTwoSigmoidalMembershipFunction::CDifferencTwoSigmoidalMembershipFunction(const double a1,const double c1,const double a2,const double c2)
  {
   m_a1 = a1;
   m_a2 = a2;
   m_c1 = c1;
   m_c2 = c2;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CDifferencTwoSigmoidalMembershipFunction::~CDifferencTwoSigmoidalMembershipFunction()
  {

  }
//+------------------------------------------------------------------+
//| Get argument (x axis value)                                      |
//+------------------------------------------------------------------+
double CDifferencTwoSigmoidalMembershipFunction::GetValue(const double x)
  {
   double first_equation=1/(1+exp(-m_a1*(x-m_c1)));
   double second_equation=1/(1+exp(-m_a2*(x-m_c2)));
//--- return result
   return (first_equation - second_equation);
  }
//+------------------------------------------------------------------+
//| Trapezoidal-shaped membership function                           |
//+------------------------------------------------------------------+ 
class CTrapezoidMembershipFunction : public IMembershipFunction
  {
private:
   double            m_x1;            // Parametr x1: the first point on the abscissa
   double            m_x2;            // Parametr x2: the second point on the abscissa
   double            m_x3;            // Parametr x3: the third point on the abscissa
   double            m_x4;            // Parametr x4: the fourth point on the abscissa

public:
                     CTrapezoidMembershipFunction(void);
                     CTrapezoidMembershipFunction(const double x1,const double x2,const double x3,const double x4);
                    ~CTrapezoidMembershipFunction(void);
   //--- methods gets or sets the parametrs  
   void              X1(const double x)   { m_x1=x;       }
   double            X1(void)             { return(m_x1); }
   void              X2(const double x)   { m_x2=x;       }
   double            X2(void)             { return(m_x2); }
   void              X3(const double x)   { m_x3=x;       }
   double            X3(void)             { return(m_x3); }
   void              X4(const double x)   { m_x4=x;       }
   double            X4(void)             { return(m_x4); }
   //--- method gets the argument (x axis value)
   double            GetValue(const double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CTrapezoidMembershipFunction::CTrapezoidMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+
CTrapezoidMembershipFunction::CTrapezoidMembershipFunction(const double x1,const double x2,const double x3,const double x4)
  {
   if(!(x1<=x2 && x2<=x3 && x3<=x4))
     {
      Print("Incorrect parameters! It is necessary to re-initialize them.");
     }
   else
     {
      m_x1 = x1;
      m_x2 = x2;
      m_x3 = x3;
      m_x4 = x4;
     }
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CTrapezoidMembershipFunction::~CTrapezoidMembershipFunction()
  {

  }
//+------------------------------------------------------------------+
//| Get argument (x axis value)                                      |
//+------------------------------------------------------------------+ 
double CTrapezoidMembershipFunction::GetValue(const double x)
  {
   double result=0;
   if(x==m_x1 && x==m_x2)
     {
      result=1.0;
     }
   else if(x==m_x3 && x==m_x4)
     {
      result=1.0;
     }
   else if(x<=m_x1 || x>=m_x4)
     {
      result=0;
     }
   else if((x>=m_x2) && (x<=m_x3))
     {
      result=1;
     }
   else if((x>m_x1) && (x<m_x2))
     {
      result=(x/(m_x2-m_x1)) -(m_x1/(m_x2-m_x1));
     }
   else
     {
      result=(-x/(m_x4-m_x3))+(m_x4/(m_x4-m_x3));
     }
//--- return result
   return (result);
  }
//+------------------------------------------------------------------+
//| Gaussian curve membership function                               |
//+------------------------------------------------------------------+  
class CNormalMembershipFunction : public IMembershipFunction
  {
private:
   double            m_b;             // Parametr b: coordinate of the maximum membership function
   double            m_sigma;         // Parametr sigma: concentration factor of the membership function

public:
                     CNormalMembershipFunction(void);
                     CNormalMembershipFunction(const double b,const double sigma);
                    ~CNormalMembershipFunction(void);
   //--- methods gets or sets the parametrs  
   void              B(const double b)          { m_b=b;           }
   double            B(void)                    { return(m_b);     }
   void              Sigma(const double sigma)  { m_sigma=sigma;   }
   double            Sigma(void)                { return(m_sigma); }
   //--- method gets the argument (x axis value)
   double            GetValue(const double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+ 
CNormalMembershipFunction::CNormalMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+ 
CNormalMembershipFunction::CNormalMembershipFunction(const double b,const double sigma)
  {
   m_b=b;
   m_sigma=sigma;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CNormalMembershipFunction::~CNormalMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Get argument (x axis value)                                      |
//+------------------------------------------------------------------+ 
double CNormalMembershipFunction::GetValue(const double x)
  {
//--- return result
   return (exp(-(x - m_b) * (x - m_b) / (2.0 * m_sigma * m_sigma)));
  }
//+------------------------------------------------------------------+
//| Triangular-shaped membership function                            |
//+------------------------------------------------------------------+
class CTriangularMembershipFunction : public IMembershipFunction
  {
private:
   double            m_x1;            // Parametr x1: the first point on the abscissa
   double            m_x2;            // Parametr x2: the second point on the abscissa
   double            m_x3;            // Parametr x3: the third point on the abscissa

public:
                     CTriangularMembershipFunction(void);
                     CTriangularMembershipFunction(const double x1,const double x2,const double x3);
                    ~CTriangularMembershipFunction(void);
   //--- methods gets or sets the parametrs  
   void              X1(const double x)   { m_x1=x;       }
   double            X1(void)             { return(m_x1); }
   void              X2(const double x)   { m_x2=x;       }
   double            X2(void)             { return(m_x2); }
   void              X3(const double x)   { m_x3=x;       }
   double            X3(void)             { return(m_x3); }
   //--- method convert triangular membership function to normal
   CNormalMembershipFunction *ToNormalMF(void);
   //--- method gets the argument (x axis value)
   double            GetValue(const double x);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CTriangularMembershipFunction::CTriangularMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+
CTriangularMembershipFunction::CTriangularMembershipFunction(const double x1,const double x2,const double x3)
  {
   if(!(x1<=x2 && x2<=x3))
     {
      Print("Incorrect parameters! It is necessary to re-initialize them.");
      return;
     }
   m_x1=x1;
   m_x2=x2;
   m_x3=x3;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CTriangularMembershipFunction::~CTriangularMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Approximately converts triangular membership function to normal  |
//+------------------------------------------------------------------+ 
CNormalMembershipFunction *CTriangularMembershipFunction::ToNormalMF(void)
  {
   double b=m_x2;
   double sigma25=(m_x3-m_x1)/2.0;
   double sigma=sigma25/2.5;
//--- return normal membership function
   return new CNormalMembershipFunction(b, sigma);
  }
//+------------------------------------------------------------------+
//| Get argument (x axis value)                                      |
//+------------------------------------------------------------------+ 
double CTriangularMembershipFunction::GetValue(const double x)
  {
   double result=0;
   if(x==m_x1 && x==m_x2)
     {
      result=1.0;
     }
   else if(x==m_x2 && x==m_x3)
     {
      result=1.0;
     }
   else if(x<=m_x1 || x>=m_x3)
     {
      result=0;
     }
   else if(x==m_x2)
     {
      result=1;
     }
   else if((x>m_x1) && (x<m_x2))
     {
      result=(x/(m_x2-m_x1)) -(m_x1/(m_x2-m_x1));
     }
   else
     {
      result=(-x/(m_x3-m_x2))+(m_x3/(m_x3-m_x2));
     }
//--- return result
   return (result);
  }
//+------------------------------------------------------------------+
//| Constant membership function                                     |
//+------------------------------------------------------------------+  
class CConstantMembershipFunction : public IMembershipFunction
  {
private:
   double            m_constValue;     // Value of the function at all points

public:
                     CConstantMembershipFunction();
                     CConstantMembershipFunction(const double constValue);
                    ~CConstantMembershipFunction();
   //--- method gets the argument (x axis value)
   double GetValue(const double x) { return(m_constValue); }
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CConstantMembershipFunction::CConstantMembershipFunction(void)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+
CConstantMembershipFunction::CConstantMembershipFunction(const double constValue)
  {
   if(constValue<0.0 || constValue>1.0)
     {
      Print("Incorrect parameter! It is necessary to re-initialize them.");
     }
   m_constValue=constValue;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CConstantMembershipFunction::~CConstantMembershipFunction(void)
  {

  }
//+--------------------------------------------------------------------------------------+
//| Composition of several membership functions represened as single membership function |
//+--------------------------------------------------------------------------------------+ 
class CCompositeMembershipFunction : public IMembershipFunction
  {
private:
   CList            *m_mfs;           // List of membership functions
   MfCompositionType m_composType;    // Composite Type

public:
                     CCompositeMembershipFunction(MfCompositionType composType);
                     CCompositeMembershipFunction(MfCompositionType composType,IMembershipFunction *mf1,IMembershipFunction *mf2);
                     CCompositeMembershipFunction(MfCompositionType composType,CList *mfs);
                    ~CCompositeMembershipFunction(void);
   //--- methods gets or sets the parametrs  
   CList            *MembershipFunctions(void)                 { return(m_mfs);        }
   MfCompositionType CompositionType(void)                     { return(m_composType); }
   void              CompositionType(MfCompositionType value)  { m_composType=value;   }
   //--- method gets the argument (x axis value)

   //+------------------------------------------------------------------+
   //| Get argument (x axis value)                                      |
   //+------------------------------------------------------------------+
   double            GetValue(double const x);
private:
   //--- composition of the membership functions
   double            Compose(const double val1,const double val2);
  };
//+------------------------------------------------------------------+
//| First constructor with parameters                                |
//+------------------------------------------------------------------+
CCompositeMembershipFunction::CCompositeMembershipFunction(MfCompositionType composType)
  {
   m_mfs=new CList;
   m_composType=composType;
  }
//+------------------------------------------------------------------+
//| Second constructor with parameters                               |
//+------------------------------------------------------------------+
CCompositeMembershipFunction::CCompositeMembershipFunction(MfCompositionType composType,IMembershipFunction *mf1,IMembershipFunction *mf2)
  {
   m_mfs=new CList;
   m_mfs.Add(mf1);
   m_mfs.Add(mf2);
   m_composType=composType;
  }
//+------------------------------------------------------------------+
//| Third constructor with parameters                                |
//+------------------------------------------------------------------+
CCompositeMembershipFunction::CCompositeMembershipFunction(MfCompositionType composType,CList *mfs)
  {
   m_mfs=mfs;
   m_composType=composType;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+     
CCompositeMembershipFunction::~CCompositeMembershipFunction()
  {
   m_mfs.FreeMode(false);
   delete m_mfs;
  }
//+------------------------------------------------------------------+
//| Get argument (x axis value)                                      |
//+------------------------------------------------------------------+
double CCompositeMembershipFunction::GetValue(double const x)
  {
   if(m_mfs.Total()==0)
     {
      //--- return result
      return 0.0;
     }
   else if(m_mfs.Total()==1)
     {
      IMembershipFunction *fun=m_mfs.GetNodeAtIndex(0);
      //--- return result
      return fun.GetValue(x);
     }
   else
     {
      IMembershipFunction *fun=m_mfs.GetNodeAtIndex(0);
      double result=fun.GetValue(x);
      for(int i=1; i<m_mfs.Total(); i++)
        {
         fun=m_mfs.GetNodeAtIndex(i);
         result=Compose(result,fun.GetValue(x));
        }
      //--- return result
      return (result);
     }
  }
//+------------------------------------------------------------------+
//| The composition of the membership functions                      |
//+------------------------------------------------------------------+
double CCompositeMembershipFunction::Compose(const double val1,const double val2)
  {
   switch(m_composType)
     {
      case MaxMF:
         //--- return result of composition
         return fmax(val1, val2);
      case MinMF:
         //--- return result of composition
         return fmin(val1, val2);
      case ProdMF:
         //--- return result of composition
         return (val1 * val2);
      case SumMF:
         //--- return result of composition
         return (val1 + val2);
      default:
        {
         Print("Incorrect type of composition");
         //--- return 
         return (NULL);
        }
     }
  }
//+------------------------------------------------------------------+
