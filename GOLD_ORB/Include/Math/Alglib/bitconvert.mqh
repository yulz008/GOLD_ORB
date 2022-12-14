//+------------------------------------------------------------------+
//|                                                   bitconvert.mqh |
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
#include "arrayresize.mqh"
//+------------------------------------------------------------------+
//| Converting numbers into an array of bits and vice versa          |
//+------------------------------------------------------------------+
class BitConverter
  {
public:
                     BitConverter(void);
                    ~BitConverter(void);

   static void       GetBytes(const int d,uchar &bytes[]);
   static void       GetBytes(const double d,uchar &bytes[]);
   static int        ToInt32(uchar &bytes[]);
   static double     ToDouble(uchar &bytes[]);
   static bool       IsLittleEndian(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
BitConverter::BitConverter(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
BitConverter::~BitConverter(void)
  {

  }
//+------------------------------------------------------------------+
//| Converting integer to a byte array                               |
//+------------------------------------------------------------------+
static void BitConverter::GetBytes(const int d,uchar &bytes[])
  {
//--- create variables
   int x;
   int q;
   int r;
   int i;
   int div;
//--- allocation
   ArrayResizeAL(bytes,4);
   for(i=0;i<4;i++)
     {
      //--- check
      if(d>=0)
         bytes[i]=0;
      else
         bytes[i]=255;
     }
//--- initialization
   q=-1;
   r=-1;
   i=3;
   div=256*256*256;
//--- check
   if(d<0)
      x=~d;
   else
      x=d;
//--- converting number
   while(i!=-1)
     {
      //--- quotient
      q=x/div;
      //--- remainder of division
      r=x%div;
      //--- get byte
      if(d>=0)
         bytes[i]+=(uchar)q;
      else
         bytes[i]-=(uchar)q;
      //--- the next iteration is reduced divisor
      x=r;
      div=div/256;
      i--;
     }
  }
//+------------------------------------------------------------------+
//| Converting double to a byte array                                |
//+------------------------------------------------------------------+
static void BitConverter::GetBytes(const double d,uchar &bytes[])
  {
//--- module
   double abs_d=MathAbs(d);
//--- number without its fractional
   double floor_d=MathFloor(abs_d);
//--- fractional part
   double fractional_d=abs_d-floor_d;
//--- variable will store the degree
   int    power;
//--- exponent shift
   double exp_shift;
//--- create variables
   int    k;
   int    j;
   uchar  u;
   double step;
   double f;
//--- abs_d as bits
   bool abs_d_to_bitArray[];
//--- d as bits in format IEEE 754
   bool d_to_bitArray[];
//--- allocation
   ArrayResizeAL(d_to_bitArray,64);
   ArrayResizeAL(bytes,8);
//--- initialization
   power=0;
//--- for integer part
   while(1)
     {
      //--- if the number is less than or equal floor_d, we increase the degree
      //--- if 2^power > floor_d, then maximal number < floor_d - it 2^(power-1)
      if(floor_d>=MathPow(2,power))
         power++;
      else
         break;
     }
//--- get power-1
   power--;
//--- if power=-1, then floor_d=0
//--- find a negative power for the fractional part
   if(power==-1)
     {
      power=0;
      while(1)
        {
         //--- the same principle as above
         if(fractional_d<MathPow(2,power))
            power--;
         else
            break;
        }
     }
//--- convert decimal to binary
   j=0;
   step=power;
   while(abs_d!=0.0)
     {
      f=MathPow(2,step);
      //--- check
      if(f>abs_d)
        {
         //--- if abs_d < f - this bit is zero
         ArrayResizeAL(abs_d_to_bitArray,j+1);
         abs_d_to_bitArray[j]=0;
         j++;
         step-=1;
        }
      else
        {
         //--- if abs_d >= f, then this bit is equal one
         //--- reduction abs_d
         abs_d-=f;
         ArrayResizeAL(abs_d_to_bitArray,j+1);
         abs_d_to_bitArray[j]=1;
         j++;
         step-=1;
        }
     }
//--- according to IEEE 754,
//--- zero bit determines sign of the number, 0 -> '+', 1 -> '-'.
   if(d>=0)
      d_to_bitArray[0]=0;
   else
      d_to_bitArray[0]=1;
//--- offset input
   exp_shift=1023+power;
//--- bits from the first and 11 are reserved for the shifted exponential
   j=1;
   for(int i=10;i>=0;i--)
     {
      if(MathPow(2,i)>exp_shift)
         d_to_bitArray[j]=0;
      else
        {
         d_to_bitArray[j]=1;
         //--- reduction
         exp_shift-=MathPow(2,i);
        }
      j++;
     }
//--- Get the length of the array of the binary representation of abs_d
   k=ArraySize(abs_d_to_bitArray);
   j=1;
//--- Bits from 12 to 63 are filled with binary representation of abs_b
//--- the first element abs_d_to_bitArray is always 1
   for(int i=12;i<64;i++)
     {
      if(j<k)
        {
         d_to_bitArray[i]=abs_d_to_bitArray[i-11];
         j++;
        }
      else
         d_to_bitArray[i]=0;
     }
//--- reverse
   ArrayReverse(d_to_bitArray);
//--- converting bit array to byte array
   for(int i=0;i<8;i++)
     {
      u=0;
      //--- get byte
      for(int t=0;t<8;t++)
         u+=(uchar)(d_to_bitArray[i*8+t]*MathPow(2,t));
      //--- save byte
      bytes[i]=u;
     }
  }
//+------------------------------------------------------------------+
//| Converting byte array to a integer                               |
//+------------------------------------------------------------------+
static int BitConverter::ToInt32(uchar &bytes[])
  {
//--- get size
   int size=ArraySize(bytes);
//--- create variables
   int d=0;
   int mul=256*256*256;
//--- get number
   for(int i=size-1;i>=0;i--)
     {
      d+=bytes[i]*mul;
      mul=mul/256;
     }
//--- return result
   return(d);
  }
//+------------------------------------------------------------------+
//| Converting byte array to a double                                |
//+------------------------------------------------------------------+
static double BitConverter::ToDouble(uchar &bytes[])
  {
//--- create variables
   int    s;
//--- exponent shift
   int    e=0;
//--- mantissa
   double m=0;
//--- array of bits in IEEE 754
   bool bits[];
   ArrayResizeAL(bits,64);
//--- get array of bits from array of bytes
   for(int i=0;i<8;i++)
     {
      for(int j=7;j>=0;j--)
        {
         //--- if 2 in power >, bits[i*8+j]=0, else bits[i*8+j]=0
         if(MathPow(2,j)>bytes[i])
            bits[i*8+j]=0;
         else
           {
            bits[i*8+j]=1;
            //--- reduction
            bytes[i]-=(uchar)MathPow(2,j);
           }
        }
     }
//--- search bits with 1
   bool allzero=true;
   for(int i=0;i<64;i++)
      if(bits[i]==1)
         allzero=false;
//--- if all bits are 0, then number is 0
   if(allzero==true)
      return(0.0);
//--- reverse array
   ArrayReverse(bits);
//--- s-the first bit, determines sign of the number
   s=bits[0];
//--- calculation exponent shift
   for(int i=10;i>=0;i--)
      e+=(int)(bits[11-i]*MathPow(2,i));
//--- get mantissa
   for(int i=0;i<52;i++)
      m+=bits[12+i]*MathPow(2,-1-i);
//--- return result
   return(MathPow(-1,s)*MathPow(2,e-1023)*(1+m));
  }
//+------------------------------------------------------------------+
//| Byte ordering (forward, backward)                                |
//+------------------------------------------------------------------+
static bool BitConverter::IsLittleEndian(void)
  {
//--- forward
   return(true);
  }
//+------------------------------------------------------------------+
//| Array reverse (uchar)                                            |
//+------------------------------------------------------------------+
void ArrayReverse(uchar &array[])
  {
//--- size
   int size=ArraySize(array);
//--- half of size
   int half=size/2;
//--- create a variable
   uchar temp;
//--- reverse
   for(int i=0;i<half;i++)
     {
      temp=array[i];
      array[i]=array[size-1-i];
      array[size-1-i]=temp;
     }
  }
//+------------------------------------------------------------------+
//| Array reverse (bool)                                             |
//+------------------------------------------------------------------+
void ArrayReverse(bool &array[])
  {
//--- size
   int size=ArraySize(array);
//--- half of size
   int half=size/2;
//--- create a variable
   bool temp;
//--- reverse
   for(int i=0;i<half;i++)
     {
      temp=array[i];
      array[i]=array[size-1-i];
      array[size-1-i]=temp;
     }
  }
//+------------------------------------------------------------------+
//| Get string from char array                                       |
//+------------------------------------------------------------------+
string GetSelectionString(char &buf[],int startIndex,int lenght)
  {
//--- create array
   ushort res[];
   int    size=ArraySize(buf);
   ArrayResizeAL(res,size);
//--- converting
   for(int i=0;i<size;i++)
      res[i]=(ushort)buf[i];
//--- return result
   return(ShortArrayToString(res,startIndex,lenght));
  }
//+------------------------------------------------------------------+
//| Get sign of number                                               |
//+------------------------------------------------------------------+
double MathSign(const double x)
  {
//--- if x>0
   if(x>0)
      return(1);
//--- if ?==0
   if(x==0)
      return(0);
//--- ?<0
   return(-1);
  }
#ifndef __MQL5__
//+------------------------------------------------------------------+
//| Hyperbolic sine                                                  |
//+------------------------------------------------------------------+
double MathSinh(const double x)
  {
//--- return result
   return((MathPow(M_E,x)-MathPow(M_E,-x))/2);
  }
//+------------------------------------------------------------------+
//| Hyperbolic cosine                                                |
//+------------------------------------------------------------------+
double MathCosh(const double x)
  {
//--- return result
   return((MathPow(M_E,x)+MathPow(M_E,-x))/2);
  }
//+------------------------------------------------------------------+
//| Hyperbolic tangent                                               |
//+------------------------------------------------------------------+
double MathTanh(const double x)
  {
//--- return result
   return(MathSinh(x)/MathCosh(x));
  }
#endif
//+------------------------------------------------------------------+
//| Structure stores a variable of type double                       |
//+------------------------------------------------------------------+
union UDoubleValue
  {
   double            value;
   long              bits;

                     UDoubleValue(double dbl):value(dbl) { }
                     UDoubleValue(long bit_value):bits(bit_value) { }
  };
//+------------------------------------------------------------------+
//| Work with infinity and NaN                                       |
//+------------------------------------------------------------------+
class CInfOrNaN
  {
public:
   //--- constructor, destructor
                     CInfOrNaN(void);
                    ~CInfOrNaN(void);
   //--- checks
   static bool       IsPositiveInfinity(const double x);
   static bool       IsNegativeInfinity(const double x);
   static bool       IsInfinity(const double x);
   static bool       IsNaN(const double x);
   //--- generation values
   static double     PositiveInfinity(void);
   static double     NegativeInfinity(void);
   static double     NaN(void);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CInfOrNaN::CInfOrNaN(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CInfOrNaN::~CInfOrNaN(void)
  {

  }
//+------------------------------------------------------------------+
//| Check for +inf                                                   |
//+------------------------------------------------------------------+
static bool CInfOrNaN::IsPositiveInfinity(const double x)
  {
   UDoubleValue val=x;
//--- check
   return(val.bits==0x7FF0000000000000);
  }
//+------------------------------------------------------------------+
//| Check for -inf                                                   |
//+------------------------------------------------------------------+
static bool CInfOrNaN::IsNegativeInfinity(const double x)
  {
   UDoubleValue val=x;
//--- check
   return(val.bits==0xFFF0000000000000);
  }
//+------------------------------------------------------------------+
//| Check for +-inf                                                  |
//+------------------------------------------------------------------+
static bool CInfOrNaN::IsInfinity(const double x)
  {
   UDoubleValue val=x;
//--- check
   return(val.bits==0x7FF0000000000000 || val.bits==0xFFF0000000000000);
  }
//+------------------------------------------------------------------+
//| Check for NaN                                                    |
//+------------------------------------------------------------------+
static bool CInfOrNaN::IsNaN(const double x)
  {
//--- check
   if(MathIsValidNumber(x))
     {
      //--- is valid number
      return(false);
     }
//--- check
   if(IsInfinity(x))
     {
      //--- +-inf
      return(false);
     }
//--- is NaN
   return(true);
  }
//+------------------------------------------------------------------+
//| Return +inf                                                      |
//+------------------------------------------------------------------+
static double CInfOrNaN::PositiveInfinity(void)
  {
   UDoubleValue val(0x7FF0000000000000);
   return(val.value);
  }
//+------------------------------------------------------------------+
//| Return -inf                                                      |
//+------------------------------------------------------------------+
static double CInfOrNaN::NegativeInfinity(void)
  {
   UDoubleValue val(0xFFF0000000000000);
   return(val.value);
  }
//+------------------------------------------------------------------+
//| Return NaN                                                       |
//+------------------------------------------------------------------+
static double CInfOrNaN::NaN(void)
  {
   UDoubleValue val(0x7FFFFFFFFFFFFFFF);
   return(val.value);
  }
//+------------------------------------------------------------------+
