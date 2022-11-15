//+------------------------------------------------------------------+
//|                                                       matrix.mqh |
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
#include "complex.mqh"
#include "arrayresize.mqh"
//+------------------------------------------------------------------+
//| Rows (double)                                                    |
//+------------------------------------------------------------------+
class CRowDouble
  {
private:
   double            m_array[];

public:
                     CRowDouble(void);
                    ~CRowDouble(void);
   //--- methods
   int               Size(void) const;
   void              Resize(const int n);
   void              Set(const int i,const double d);
   //--- overloading
   double            operator[](const int i) const;
   void              operator=(const double &array[]);
   void              operator=(const CRowDouble &r);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CRowDouble::CRowDouble(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CRowDouble::~CRowDouble(void)
  {

  }
//+------------------------------------------------------------------+
//| Row size                                                         |
//+------------------------------------------------------------------+
int CRowDouble::Size(void) const
  {
   return(ArraySize(m_array));
  }
//+------------------------------------------------------------------+
//| Resize                                                           |
//+------------------------------------------------------------------+
void CRowDouble::Resize(const int n)
  {
   ArrayResizeAL(m_array,n);
  }
//+------------------------------------------------------------------+
//| Set value                                                        |
//+------------------------------------------------------------------+
void CRowDouble::Set(const int i,const double d)
  {
   m_array[i]=d;
  }
//+------------------------------------------------------------------+
//| Indexing operator                                                |
//+------------------------------------------------------------------+
double CRowDouble::operator[](const int i) const
  {
   return(m_array[i]);
  }
//+------------------------------------------------------------------+
//| Overloading (=)                                                  |
//+------------------------------------------------------------------+
void CRowDouble::operator=(const double &array[])
  {
   int size=ArraySize(array);
//--- check
   if(size==0)
      return;
//--- filling array
   ArrayResizeAL(m_array,size);
   for(int i=0;i<size;i++)
      m_array[i]=array[i];
  }
//+------------------------------------------------------------------+
//| Overloading (=)                                                  |
//+------------------------------------------------------------------+
void CRowDouble::operator=(const CRowDouble &r)
  {
   int size=r.Size();
//--- check
   if(size==0)
      return;
//--- filling array
   ArrayResizeAL(m_array,size);
   for(int i=0;i<size;i++)
      m_array[i]=r[i];
  }
//+------------------------------------------------------------------+
//| Rows (int)                                                       |
//+------------------------------------------------------------------+
class CRowInt
  {
private:
   int               m_array[];

public:
                     CRowInt(void);
                    ~CRowInt(void);
   //--- methods
   int               Size(void) const;
   void              Resize(const int n);
   void              Set(const int i,const int d);
   //--- overloading
   int               operator[](const int i) const;
   void              operator=(const int &array[]);
   void              operator=(const CRowInt &r);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CRowInt::CRowInt(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CRowInt::~CRowInt(void)
  {

  }
//+------------------------------------------------------------------+
//| Row size                                                         |
//+------------------------------------------------------------------+
int CRowInt::Size(void) const
  {
   return(ArraySize(m_array));
  }
//+------------------------------------------------------------------+
//| Resize                                                           |
//+------------------------------------------------------------------+
void CRowInt::Resize(const int n)
  {
   ArrayResizeAL(m_array,n);
  }
//+------------------------------------------------------------------+
//| Set value                                                        |
//+------------------------------------------------------------------+
void CRowInt::Set(const int i,const int d)
  {
   m_array[i]=d;
  }
//+------------------------------------------------------------------+
//| Indexing operator                                                |
//+------------------------------------------------------------------+
int CRowInt::operator[](const int i) const
  {
   return(m_array[i]);
  }
//+------------------------------------------------------------------+
//| Overloading (=)                                                  |
//+------------------------------------------------------------------+
void CRowInt::operator=(const int &array[])
  {
   int size=ArraySize(array);
//--- check
   if(size==0)
      return;
//--- filling array
   ArrayResizeAL(m_array,size);
   for(int i=0;i<size;i++)
      m_array[i]=array[i];
  }
//+------------------------------------------------------------------+
//| Overloading (=)                                                  |
//+------------------------------------------------------------------+
void CRowInt::operator=(const CRowInt &r)
  {
   int size=r.Size();
//--- check
   if(size==0)
      return;
//--- filling array
   ArrayResizeAL(m_array,size);
   for(int i=0;i<size;i++)
      m_array[i]=r[i];
  }
//+------------------------------------------------------------------+
//| Rows (complex)                                                   |
//+------------------------------------------------------------------+
class CRowComplex
  {
private:
   //--- array
   al_complex        m_array[];
public:
   //--- constructor, destructor
                     CRowComplex(void);
                    ~CRowComplex(void);
   //--- methods
   int               Size(void) const;
   void              Resize(const int n);
   void              Set(const int i,const al_complex &d);
   void              Set(const int i,const double d);
   void              SetRe(const int i,const double d);
   void              SetIm(const int i,const double d);
   //--- overloading
   al_complex        operator[](const int i) const;
   void              operator=(const al_complex &array[]);
   void              operator=(const CRowComplex &r);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CRowComplex::CRowComplex(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CRowComplex::~CRowComplex(void)
  {

  }
//+------------------------------------------------------------------+
//| Row size                                                         |
//+------------------------------------------------------------------+
int CRowComplex::Size(void) const
  {
   return(ArraySize(m_array));
  }
//+------------------------------------------------------------------+
//| Resize                                                           |
//+------------------------------------------------------------------+
void CRowComplex::Resize(const int n)
  {
   ArrayResizeAL(m_array,n);
  }
//+------------------------------------------------------------------+
//| Set value                                                        |
//+------------------------------------------------------------------+
void CRowComplex::Set(const int i,const al_complex &d)
  {
   m_array[i]=d;
  }
//+------------------------------------------------------------------+
//| Set value                                                        |
//+------------------------------------------------------------------+
void CRowComplex::Set(const int i,const double d)
  {
   al_complex c(d,0);
   m_array[i]=c;
  }
//+------------------------------------------------------------------+
//| Set real part                                                    |
//+------------------------------------------------------------------+
void CRowComplex::SetRe(const int i,const double d)
  {
   m_array[i].re=d;
  }
//+------------------------------------------------------------------+
//| Set imaginary part                                               |
//+------------------------------------------------------------------+
void CRowComplex::SetIm(const int i,const double d)
  {
   m_array[i].im=d;
  }
//+------------------------------------------------------------------+
//| Indexing operator                                                |
//+------------------------------------------------------------------+
al_complex CRowComplex::operator[](const int i) const
  {
   return(m_array[i]);
  }
//+------------------------------------------------------------------+
//| Overloading (=)                                                  |
//+------------------------------------------------------------------+
void CRowComplex::operator=(const al_complex &array[])
  {
   int size=ArraySize(array);
//--- check
   if(size==0)
      return;
//--- filling array
   ArrayResizeAL(m_array,size);
   for(int i=0;i<size;i++)
      m_array[i]=array[i];
  }
//+------------------------------------------------------------------+
//| Overloading (=)                                                  |
//+------------------------------------------------------------------+
void CRowComplex::operator=(const CRowComplex &r)
  {
   int size=r.Size();
//--- check
   if(size==0)
      return;
//--- filling array
   ArrayResizeAL(m_array,size);
   for(int i=0;i<size;i++)
      m_array[i]=r[i];
  }
//+------------------------------------------------------------------+
//| Matrix (double)                                                  |
//+------------------------------------------------------------------+
class CMatrixDouble
  {
private:
   //--- array
   CRowDouble        m_rows[];
public:
   //--- constructors, destructor
                     CMatrixDouble(void);
                     CMatrixDouble(const int rows);
                     CMatrixDouble(const int rows,const int cols);
                    ~CMatrixDouble(void);
   //--- methods
   int               Size(void) const;
   void              Resize(const int n,const int m);
   //--- overloading
   CRowDouble       *operator[](const int i) const;
   void              operator=(const CMatrixDouble &m);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMatrixDouble::CMatrixDouble(void)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with one parameter                                   |
//+------------------------------------------------------------------+
CMatrixDouble::CMatrixDouble(const int rows)
  {
   ArrayResizeAL(m_rows,rows);
  }
//+------------------------------------------------------------------+
//| Constructor with two parameters                                  |
//+------------------------------------------------------------------+
CMatrixDouble::CMatrixDouble(const int rows,const int cols)
  {
   ArrayResizeAL(m_rows,rows);
   for(int i=0;i<rows;i++)
      m_rows[i].Resize(cols);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMatrixDouble::~CMatrixDouble(void)
  {

  }
//+------------------------------------------------------------------+
//| Get size                                                         |
//+------------------------------------------------------------------+
int CMatrixDouble::Size(void) const
  {
   return(ArraySize(m_rows));
  }
//+------------------------------------------------------------------+
//| Resize                                                           |
//+------------------------------------------------------------------+
void CMatrixDouble::Resize(const int n,const int m)
  {
//--- check
   if(n<0 || m<0)
      return;
//--- change sizes
   ArrayResizeAL(m_rows,n);
   for(int i=0;i<n;i++)
      m_rows[i].Resize(m);
  }
//+------------------------------------------------------------------+
//| Indexing operator                                                |
//+------------------------------------------------------------------+
CRowDouble *CMatrixDouble::operator[](const int i) const
  {
   return(GetPointer(m_rows[i]));
  }
//+------------------------------------------------------------------+
//| Overloading (=)                                                  |
//+------------------------------------------------------------------+
void CMatrixDouble::operator=(const CMatrixDouble &m)
  {
   int r=m.Size();
//--- check
   if(r==0)
      return;
   int c=m[0].Size();
//--- check
   if(c==0)
      return;
//--- change size
   ArrayResizeAL(m_rows,r);
   for(int i=0;i<r;i++)
      m_rows[i].Resize(c);
//--- copy
   for(int i=0;i<r;i++)
      m_rows[i]=m[i];
  }
//+------------------------------------------------------------------+
//| Matrix (int)                                                     |
//+------------------------------------------------------------------+
class CMatrixInt
  {
private:
   //--- array
   CRowInt           m_rows[];
public:
   //--- constructors, destructor
                     CMatrixInt(void);
                     CMatrixInt(const int rows);
                     CMatrixInt(const int rows,const int cols);
                    ~CMatrixInt(void);
   //--- methods
   int               Size(void) const;
   void              Resize(const int n,const int m);
   //--- overloading
   CRowInt          *operator[](const int i) const;
   void              operator=(const CMatrixInt &m);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMatrixInt::CMatrixInt(void)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with one parameter                                   |
//+------------------------------------------------------------------+
CMatrixInt::CMatrixInt(const int rows)
  {
   ArrayResizeAL(m_rows,rows);
  }
//+------------------------------------------------------------------+
//| Constructor with two parameters                                  |
//+------------------------------------------------------------------+
CMatrixInt::CMatrixInt(const int rows,const int cols)
  {
   ArrayResizeAL(m_rows,rows);
   for(int i=0;i<rows;i++)
      m_rows[i].Resize(cols);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMatrixInt::~CMatrixInt(void)
  {

  }
//+------------------------------------------------------------------+
//| Get size                                                         |
//+------------------------------------------------------------------+
int CMatrixInt::Size(void) const
  {
   return(ArraySize(m_rows));
  }
//+------------------------------------------------------------------+
//| Resize                                                           |
//+------------------------------------------------------------------+
void CMatrixInt::Resize(const int n,const int m)
  {
//--- check
   if(n<0 || m<0)
      return;
//--- change sizes
   ArrayResizeAL(m_rows,n);
   for(int i=0;i<n;i++)
      m_rows[i].Resize(m);
  }
//+------------------------------------------------------------------+
//| Indexing operator                                                |
//+------------------------------------------------------------------+
CRowInt *CMatrixInt::operator[](const int i) const
  {
   return(GetPointer(m_rows[i]));
  }
//+------------------------------------------------------------------+
//| Overloading (=)                                                  |
//+------------------------------------------------------------------+
void CMatrixInt::operator=(const CMatrixInt &m)
  {
   int r=m.Size();
//--- check
   if(r==0)
      return;
   int c=m[0].Size();
//--- check
   if(c==0)
      return;
//--- change size
   ArrayResizeAL(m_rows,r);
   for(int i=0;i<r;i++)
      m_rows[i].Resize(c);
//--- copy
   for(int i=0;i<r;i++)
      m_rows[i]=m[i];
  }
//+------------------------------------------------------------------+
//| Matrix (complex)                                                 |
//+------------------------------------------------------------------+
class CMatrixComplex
  {
private:
   //--- array
   CRowComplex       m_rows[];
public:
   //--- constructors, destructor
                     CMatrixComplex(void);
                     CMatrixComplex(const int rows);
                     CMatrixComplex(const int rows,const int cols);
                    ~CMatrixComplex(void);
   //--- methods
   int               Size(void) const;
   void              Resize(const int n,const int m);
   //--- overloading
   CRowComplex      *operator[](const int i) const;
   void              operator=(const CMatrixComplex &m);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMatrixComplex::CMatrixComplex(void)
  {

  }
//+------------------------------------------------------------------+
//| Constructor with one parameter                                   |
//+------------------------------------------------------------------+
CMatrixComplex::CMatrixComplex(const int rows)
  {
   ArrayResizeAL(m_rows,rows);
  }
//+------------------------------------------------------------------+
//| Constructor with two parameters                                  |
//+------------------------------------------------------------------+
CMatrixComplex::CMatrixComplex(const int rows,const int cols)
  {
   ArrayResizeAL(m_rows,rows);
   for(int i=0;i<rows;i++)
      m_rows[i].Resize(cols);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CMatrixComplex::~CMatrixComplex(void)
  {

  }
//+------------------------------------------------------------------+
//| Get size                                                         |
//+------------------------------------------------------------------+
int CMatrixComplex::Size(void) const
  {
   return(ArraySize(m_rows));
  }
//+------------------------------------------------------------------+
//| Resize                                                           |
//+------------------------------------------------------------------+
void CMatrixComplex::Resize(const int n,const int m)
  {
//--- check
   if(n<0 || m<0)
      return;
//--- change sizes
   ArrayResizeAL(m_rows,n);
   for(int i=0;i<n;i++)
      m_rows[i].Resize(m);
  }
//+------------------------------------------------------------------+
//| Indexing operator                                                |
//+------------------------------------------------------------------+
CRowComplex *CMatrixComplex::operator[](const int i) const
  {
   return(GetPointer(m_rows[i]));
  }
//+------------------------------------------------------------------+
//| Overloading (=)                                                  |
//+------------------------------------------------------------------+
void CMatrixComplex::operator=(const CMatrixComplex &m)
  {
   int r=m.Size();
//--- check
   if(r==0)
      return;
   int c=m[0].Size();
//--- check
   if(c==0)
      return;
//--- change size
   ArrayResizeAL(m_rows,r);
   for(int i=0;i<r;i++)
      m_rows[i].Resize(c);
//--- copy
   for(int i=0;i<r;i++)
      m_rows[i]=m[i];
  }
//+------------------------------------------------------------------+
