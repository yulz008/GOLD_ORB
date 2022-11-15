//+------------------------------------------------------------------+
//|                                                    fuzzyterm.mqh |
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
#include "MembershipFunction.mqh"
#include "Helper.mqh"
//+------------------------------------------------------------------+
//| Purpose: creating fuzzy term.                                    |
//+------------------------------------------------------------------+
//+------------------------------------------------------------------+
//| Fuzzy or linguistic term.                                        |
//+------------------------------------------------------------------+
class CFuzzyTerm : public CNamedValueImpl
  {
private:
   IMembershipFunction *m_mf;         // The membership function of the term

public:
                     CFuzzyTerm(const string name,IMembershipFunction *mf);
                    ~CFuzzyTerm(void);
   //--- method to check type
   virtual bool      IsTypeOf(EnType type)   { return(type==TYPE_CLASS_FuzzyTerm); }
   //--- method gets the membership function initially associated with the term
   IMembershipFunction *MembershipFunction() { return(m_mf); }
  };
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+
CFuzzyTerm::CFuzzyTerm(const string name,IMembershipFunction *mf)
  {
   CNamedValueImpl::Name(name);
   m_mf=mf;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+     
CFuzzyTerm::~CFuzzyTerm(void)
  {
   if(CheckPointer(m_mf)==POINTER_DYNAMIC)
     {
      delete m_mf;
     }
  }
//+------------------------------------------------------------------+
