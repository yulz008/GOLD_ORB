//+------------------------------------------------------------------+
//|                                                   dictionary.mqh |
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
#include <Object.mqh>
#include <Arrays\List.mqh>
#include "RuleParser.mqh"
//+------------------------------------------------------------------+
//| Gets the value associated with the specified key in the CList    |
//| Where key - string, value - CObject                              |
//+------------------------------------------------------------------+
bool TryGetValue(CList *list,string key,CObject *&value)
  {
   for(int i=0; i<list.Total(); i++)
     {
      CDictionary_String_Obj *pair=list.GetNodeAtIndex(i);
      if(pair.Key()==key)
        {
         value=pair.Value();
         return (true);
        }
     }
   return (false);
  }
//+------------------------------------------------------------------+
//| Removes a range of elements from a list of CList                 |
//+------------------------------------------------------------------+
void RemoveRange(CArrayObj &list,const int index,const int count)
  {
   for(int i=0; i<count; i++)
     {
      list.Delete(index);
     }
  }
//+------------------------------------------------------------------+
//| It creates a shallow copy of a range of elements                 |
//| from the original list of CList                                  |
//+------------------------------------------------------------------+   
CArrayObj *GetRange(CArrayObj *list,const int index,const int count)
  {
   CArrayObj *new_list=new CArrayObj;
   for(int i=0; i<count; i++)
     {
      new_list.Add(list.At(i+index));
     }
   return (new_list);
  }
//+------------------------------------------------------------------+
//| Dictionary: Object - Object                                      |
//+------------------------------------------------------------------+
class CDictionary_Obj_Obj : public CObject
  {
private:
   CObject          *m_key;
   CObject          *m_value;

public:
                     CDictionary_Obj_Obj(void);
                    ~CDictionary_Obj_Obj(void);
   //--- methods gets or sets the value
   CObject          *Key() { return(m_key);   }
   void              Key(CObject *key) { m_key=key;       }
   //--- methods gets or sets the key    
   CObject          *Value()                 { return(m_value); }
   void              Value(CObject *value)   { m_value=value;   }
   //--- method sets the value and key
   void              SetAll(CObject *key,CObject *value);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDictionary_Obj_Obj::CDictionary_Obj_Obj(void)
  {
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+     
CDictionary_Obj_Obj::~CDictionary_Obj_Obj()
  {
  }
//+------------------------------------------------------------------+
//| Sets the value and key                                           |
//+------------------------------------------------------------------+
void CDictionary_Obj_Obj::SetAll(CObject *key,CObject *value)
  {
   m_key=key;
   m_value=value;
  }
//+------------------------------------------------------------------+
//| Dictionary: String - Object                                      |
//+------------------------------------------------------------------+
class CDictionary_String_Obj : public CObject
  {
private:
   string            m_key;
   CObject          *m_value;

public:
                     CDictionary_String_Obj(void);
                    ~CDictionary_String_Obj(void);
   //--- methods gets or sets the value
   string            Key()                   { return(m_key);   }
   void              Key(const string key)   { m_key=key;       }
   //--- methods gets or sets the key    
   CObject          *Value()                 { return(m_value); }
   void              Value(CObject *value)   { m_value=value;   }
   //--- method sets the value and key
   void              SetAll(const string key,CObject *value);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDictionary_String_Obj::CDictionary_String_Obj(void)
  {
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+     
CDictionary_String_Obj::~CDictionary_String_Obj()
  {
   if(CheckPointer(m_value)==POINTER_DYNAMIC)
      delete m_value;
  }
//+------------------------------------------------------------------+
//| Sets the value and key                                           |
//+------------------------------------------------------------------+
void CDictionary_String_Obj::SetAll(const string key,CObject *value)
  {
   m_key=key;
   m_value=value;
  }
//+------------------------------------------------------------------+
//| Dictionary: Object - Double                                      |
//+------------------------------------------------------------------+
class CDictionary_Obj_Double : public CObject
  {
private:
   CObject          *m_key;
   double            m_value;

public:
                     CDictionary_Obj_Double(void);
                    ~CDictionary_Obj_Double(void);
   //--- methods gets or sets the value
   CObject          *Key() { return(m_key);   }
   void              Key(CObject *key) { m_key=key;       }
   //--- methods gets or sets the key    
   double            Value()                    { return(m_value); }
   void              Value(const double value)  { m_value=value;   }
   //--- method sets the value and key
   void              SetAll(CObject *key,const double value);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CDictionary_Obj_Double::CDictionary_Obj_Double(void)
  {
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+     
CDictionary_Obj_Double::~CDictionary_Obj_Double()
  {
  }
//+------------------------------------------------------------------+
//| Sets the value and key                                           |
//+------------------------------------------------------------------+
void CDictionary_Obj_Double::SetAll(CObject *key,const double value)
  {
   m_key=key;
   m_value=value;
  }
//+------------------------------------------------------------------+
