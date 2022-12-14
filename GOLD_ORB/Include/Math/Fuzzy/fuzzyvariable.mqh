//+------------------------------------------------------------------+
//|                                                fuzzyvariable.mqh |
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
#include "FuzzyTerm.mqh"
//+------------------------------------------------------------------+
//| Purpose: creating fuzzy variable                                 |
//+------------------------------------------------------------------+
//+------------------------------------------------------------------+
//| Fuzzy or linguistic variable.                                    |
//+------------------------------------------------------------------+
class CFuzzyVariable : public CNamedVariableImpl
  {
private:
   double            m_min;         // Minimum value of the variable
   double            m_max;         // Maximum value of the variable
   CList            *m_terms;       // List of terms in a variable

public :
                     CFuzzyVariable(const string name,const double min,const double max);
                    ~CFuzzyVariable(void);
   //--- method to check type
   virtual bool      IsTypeOf(EnType type) { return(type==TYPE_CLASS_FuzzyVariable); }
   //--- methods gets or sets parameters of varriable
   void              Max(const double max)   { m_max=max;      }
   double            Max(void)               { return (m_max); }
   void              Min(const double min)   { m_min=min;      }
   double            Min(void)               { return (m_min); }
   //--- methods gets or sets the terms
   CList            *Terms() { return(m_terms); }
   void              Terms(CList *terms) {  m_terms=terms;  }
   //--- add fuzzy term  
   void              AddTerm(CFuzzyTerm *term);
   //--- get membership function by name 
   CFuzzyTerm       *GetTermByName(const string name);
   //--- overload
   CList            *Values() { return(m_terms); }
  };
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+   
CFuzzyVariable::CFuzzyVariable(const string name,const double min,const double max)
  {
   CNamedVariableImpl::Name(name);
   m_terms=new CList();
   if(min>max)
     {
      Print("Incorrect parameters! Maximum value must be greater than minimum one.");
     }
   else
     {
      m_min = min;
      m_max = max;
     }
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+     
CFuzzyVariable::~CFuzzyVariable(void)
  {
   delete m_terms;
  }
//+------------------------------------------------------------------+
//| Add fuzzy term to list terms in a variable                       |
//+------------------------------------------------------------------+     
void CFuzzyVariable::AddTerm(CFuzzyTerm *term)
  {
   m_terms.Add(term);
  }
//+------------------------------------------------------------------+
//| Get membership function (term) by name                           |
//+------------------------------------------------------------------+
CFuzzyTerm *CFuzzyVariable::GetTermByName(const string name)
  {
   for(int i=0; i<m_terms.Total(); i++)
     {
      CFuzzyTerm *term=m_terms.GetNodeAtIndex(i);
      if(term.Name()==name)
        {
         //--- return fuzzy term
         return (term);
        }
     }
   Print("Term with the same name can not be found!");
//--- return 
   return (NULL);
  }
//+------------------------------------------------------------------+
