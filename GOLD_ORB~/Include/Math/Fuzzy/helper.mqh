//+------------------------------------------------------------------+
//|                                                       helper.mqh |
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
#include "InferenceMethod.mqh"
//+------------------------------------------------------------------+
//| Purpose: Analysis of the fuzzy rules                             |
//+------------------------------------------------------------------+
//+------------------------------------------------------------------+
//| This class must be implemented by values in parsable rules       |
//+------------------------------------------------------------------+
class INamedValue : public CObject
  {
public:
   //--- method to check type
   virtual bool      IsTypeOf(EnType type) { return(type==TYPE_CLASS_INamedValue); }
   //--- methods gets or sets the name
   virtual string    Name(void)              { return(""); }
   virtual void      Name(const string name) {             }
  };
//+------------------------------------------------------------------+
//| This class must be implemented by values in parsable rules       |
//+------------------------------------------------------------------+
class INamedVariable : public INamedValue
  {
public:
   //--- method to check type
   virtual bool      IsTypeOf(EnType type) { return(type==TYPE_CLASS_INamedValue); }
   //--- get list of values that belongs to the variable
   virtual CList    *Values(void) { return(NULL); }
  };
//+------------------------------------------------------------------+
//| Named variable                                                   |
//+------------------------------------------------------------------+ 
class CNamedVariableImpl : public INamedVariable
  {
private:
   string            m_name;        // Name of the variable  

public:
   //--- method to check type 
   virtual bool      IsTypeOf(EnType type) { return(type==TYPE_CLASS_NamedVariableImpl); }
   //--- methods gets or sets varriable name
   virtual void      Name(const string name);
   virtual string    Name(void) { return(m_name); }
   //--- get list of values that belongs to the variable
   virtual CList    *Values(void) { return(NULL); }
  };
//+------------------------------------------------------------------+
//| Set variable name                                                |
//+------------------------------------------------------------------+  
void CNamedVariableImpl::Name(const string name)
  {
   if(!CNameHelper::IsValidName(name))
     {
      Print("Invalid variable name.");
     }
   m_name=name;
  }
//+------------------------------------------------------------------+
//| Named value of variable                                          |
//+------------------------------------------------------------------+ 
class CNamedValueImpl : public INamedValue
  {
private:
   string            m_name;        // Name of the value 

public:
   //--- method to check type 
   virtual bool      IsTypeOf(EnType type) { return(type==TYPE_CLASS_NamedVariableImpl); }
   //--- methods gets or sets varriable name
   virtual void      Name(const string name);
   virtual string    Name(void) { return(m_name); }
  };
//+------------------------------------------------------------------+
//| Set variable name                                                |
//+------------------------------------------------------------------+
void CNamedValueImpl::Name(const string name)
  {
   if(!CNameHelper::IsValidName(name))
     {
      Print("Invalid term name.");
     }
   m_name=name;
  }
//+------------------------------------------------------------------+
//| Keywords:                                                        |
//+------------------------------------------------------------------+
static string  KEYWORDS[]={ "if","then","is","and","or","not","(",")","slightly","somewhat","very","extremely" }; // Keywords in rules
//+------------------------------------------------------------------+
//| Class NameHelper checks the availability of names                |
//+------------------------------------------------------------------+
class CNameHelper
  {
public:
   //+------------------------------------------------------------------+
   //| Check the name of variable/term                                  |
   //+------------------------------------------------------------------+
   static bool IsValidName(const string name)
     {
      //--- Empty names are not allowed
      if(StringLen(name)==0)
        {
         //--- return false
         return (false);
        }

      for(int i=0; i<StringLen(name); i++)
        {
         //--- Only letters, numbers or '_' are allowed
         char s=(char) StringGetCharacter(name,i);
         if(s!='_' && !(s>=48 && s<=57)    // Not numbers and symbol '_'
            && !( s >= 65 && s <= 90 )     // Not capital letters
            && !( s >= 97 && s <= 122 ))   // Not letters
           {
            //--- return false
            return (false);
           }
        }
      //--- Identifier cannot be a keword
      for(int i=0; i<ArraySize(KEYWORDS); i++)
        {
         if(name==KEYWORDS[i])
           {
            //--- return false
            return (false);
           }
        }
      //--- return true
      return (true);
     }
  };
//+------------------------------------------------------------------+
