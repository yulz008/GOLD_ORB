//+------------------------------------------------------------------+
//|                                                    fuzzyrule.mqh |
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
#include "FuzzyVariable.mqh"
#include "InferenceMethod.mqh"
//+------------------------------------------------------------------+
//| Purpose: Creating fuzzy rules                                    |
//+------------------------------------------------------------------+
//+------------------------------------------------------------------+
//| And/Or operator type                                             |
//+------------------------------------------------------------------+ 
enum OperatorType
  {
   And,                               // And operator
   Or                                 // Or operator
  };
//+------------------------------------------------------------------+
//| Hedge modifiers                                                  |
//+------------------------------------------------------------------+ 
enum HedgeType
  {
   None,                              // None
   Slightly,                          // Cube root
   Somewhat,                          // Square root
   Very,                              // Square
   Extremely                          // Cube
  };
//+------------------------------------------------------------------+
//| Class of CConditions used in the 'if' expression                 |
//+------------------------------------------------------------------+
class ICondition : public CObject
  {
public:
   //--- method to check type 
   virtual bool      IsTypeOf(EnCondition type) { return(type==TYPE_CLASS_ICondition); }
  };
//+------------------------------------------------------------------+
//| Single condition                                                 |
//+------------------------------------------------------------------+
class CSingleCondition : public ICondition
  {
private:
   INamedVariable   *m_var;        // Type of variable
   INamedValue      *m_term;       // Type of value
   bool              m_not;        // Is MF inverted 

public:
                     CSingleCondition(void);
                     CSingleCondition(INamedVariable *var,INamedValue *term);
                     CSingleCondition(INamedVariable *var,INamedValue *term,bool not);
                    ~CSingleCondition(void);
   //--- methods gets or sets the varriable
   INamedVariable   *Var(void) { return(m_var); }
   void              Var(INamedVariable *value) { m_var=value;   }
   //--- methods gets or sets mark "Is MF inverted"
   bool              Not(void)      { return(m_not); }
   void              Not(bool not)  { m_not=not;     }
   //--- methods gets or sets term in expression
   INamedValue      *Term(void) { return(m_term); }
   void              Term(INamedValue *value) { m_term=value;   }
   //--- method to check type
   virtual bool      IsTypeOf(EnCondition type) { return(type==TYPE_CLASS_SingleCondition); }
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CSingleCondition::CSingleCondition(void)
  {
   m_var = NULL;
   m_not = false;
   m_term=NULL;
  };
//+------------------------------------------------------------------+
//| First constructor with parameters                                |
//+------------------------------------------------------------------+
CSingleCondition::CSingleCondition(INamedVariable *var,INamedValue *term)
  {
   m_var=var;
   m_term=term;
  }
//+------------------------------------------------------------------+
//| Second constructor with parameters                               |
//+------------------------------------------------------------------+
CSingleCondition::CSingleCondition(INamedVariable *var,INamedValue *term,bool not)
  {
   m_var=var;
   m_term=term;
   m_not=not;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+ 
CSingleCondition::~CSingleCondition(void)
  {
   if(CheckPointer(m_var)==POINTER_DYNAMIC)
      delete m_var;
   if(CheckPointer(m_term)==POINTER_DYNAMIC)
      delete m_term;
  }
//+------------------------------------------------------------------+
//| Condition of fuzzy rule for the both Mamdani and Sugeno systems  |
//+------------------------------------------------------------------+
class CFuzzyCondition : public CSingleCondition
  {
private:
   HedgeType         m_hedge;         // hedge type

public:
                     CFuzzyCondition(CFuzzyVariable *var,CFuzzyTerm *term,bool not);
                     CFuzzyCondition(CFuzzyVariable *var,CFuzzyTerm *term,bool not,HedgeType hedge);
                     CFuzzyCondition(CFuzzyVariable *var,CFuzzyTerm *term);
                    ~CFuzzyCondition(void);
   //--- methods gets or sets the hedge type
   HedgeType         Hedge(void)             { return (m_hedge); }
   void              Hedge(HedgeType value)  { m_hedge=value;    }
   //--- method to check type 
   virtual bool      IsTypeOf(EnCondition type) { return(type==TYPE_CLASS_FuzzyCondition); }
  };
//+------------------------------------------------------------------+
//| First constructor with parameters                                |
//+------------------------------------------------------------------+
CFuzzyCondition::CFuzzyCondition(CFuzzyVariable *var,CFuzzyTerm *term,bool not)
  {
   CSingleCondition::Var(var);
   CSingleCondition::Term(term);
   CSingleCondition::Not(not);
   m_hedge=None;
  }
//+------------------------------------------------------------------+
//| Second constructor with parameters                               |
//+------------------------------------------------------------------+     
CFuzzyCondition::CFuzzyCondition(CFuzzyVariable *var,CFuzzyTerm *term,bool not,HedgeType hedge)
  {
   CSingleCondition::Var(var);
   CSingleCondition::Term(term);
   CSingleCondition::Not(not);
   m_hedge=hedge;
  }
//+------------------------------------------------------------------+
//| Thrid constructor with parameters                                |
//+------------------------------------------------------------------+       
CFuzzyCondition::CFuzzyCondition(CFuzzyVariable *var,CFuzzyTerm *term)

  {
   CSingleCondition::Var(var);
   CSingleCondition::Term(term);
   CSingleCondition::Not(false);
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CFuzzyCondition::~CFuzzyCondition(void)
  {

  }
//+------------------------------------------------------------------+
//| Several CConditions linked by or/and operators                   |
//+------------------------------------------------------------------+
class CConditions : public ICondition
  {
private:
   bool              m_not;           // Default : false
   OperatorType      m_op;            // Type of operator. Default : And 
   CList            *m_conditions;    // List of CConditions

public:
                     CConditions(void);
                    ~CConditions(void);
   //--- methods gets or sets the mark "Is MF inverted"
   bool              Not(void)         { return(m_not); }
   void              Not(bool value)   { m_not=value;   }
   //--- methods gets or sets operator that links expressions (and/or)
   OperatorType      Op(void)                { return (m_op); }
   void              Op(OperatorType value)  { m_op=value;    }
   //--- method gets the list of CConditions (single or multiples)
   CList            *ConditionsList(void) { return(m_conditions); }
   //--- method to check type
   virtual bool      IsTypeOf(EnCondition type) { return(type==TYPE_CLASS_Conditions); }
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CConditions::CConditions(void)
  {
   m_not=false;
   m_op = And;
   m_conditions=new CList;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+ 
CConditions::~CConditions(void)
  {
   delete m_conditions;
  }
//+------------------------------------------------------------------+
//| Class used by rule parser                                        |
//+------------------------------------------------------------------+
class IParsableRule : public CObject
  {
public:
   //--- methods gets or sets the condition (IF) part of the rule
   virtual CConditions *Condition(void)                  { return(NULL); }
   virtual void      Condition(CConditions *value)       {  }
   //--- methods gets or sets the conclusion (THEN) part of the rule
   virtual  CSingleCondition *Conclusion(void)           { return(NULL); }
   virtual void      Conclusion(CSingleCondition *value) {  }
   //--- method to check type
   virtual bool      IsTypeOf(EnRule type) { return(type==TYPE_CLASS_IParsableRule); }
  };
//+------------------------------------------------------------------+
//| Implements common functionality of fuzzy rules                   |
//+------------------------------------------------------------------+
class CGenericFuzzyRule : public IParsableRule
  {
private:
   CConditions      *m_generic_condition; // Generic path of condition   

public:
                     CGenericFuzzyRule(void);
                    ~CGenericFuzzyRule(void);
   //--- methods gets or sets the condition (IF) part of the rule
   CConditions      *Condition(void)               { return(m_generic_condition); }
   void              Condition(CConditions *value) { m_generic_condition=value;   }
   //--- methods create a single condition
   CFuzzyCondition *CreateCondition(CFuzzyVariable *var,CFuzzyTerm *term);
   CFuzzyCondition *CreateCondition(CFuzzyVariable *var,CFuzzyTerm *term,bool not);
   CFuzzyCondition *CreateCondition(CFuzzyVariable *var,CFuzzyTerm *term,bool not,HedgeType hedge);
   //--- methods gets or sets the conclusion (THEN) part of the rule
   virtual CSingleCondition *Conclusion(void)            {  return(NULL);  }
   virtual void      Conclusion(CSingleCondition *value) {  }
   //--- method to check type
   virtual bool      IsTypeOf(EnRule type) { return(type==TYPE_CLASS_GenericFuzzyRule); }
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CGenericFuzzyRule::CGenericFuzzyRule(void)
  {
   m_generic_condition=new CConditions();
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+ 
CGenericFuzzyRule::~CGenericFuzzyRule(void)
  {
   delete m_generic_condition;
  }
//+------------------------------------------------------------------+
//| Create a single condition(1)                                     |
//+------------------------------------------------------------------+ 
CFuzzyCondition *CGenericFuzzyRule::CreateCondition(CFuzzyVariable *var,CFuzzyTerm *term)
  {
//--- return fuzzy condition
   return new CFuzzyCondition(var, term);
  }
//+------------------------------------------------------------------+
//| Create a single condition(2)                                     |
//+------------------------------------------------------------------+ 
CFuzzyCondition *CGenericFuzzyRule::CreateCondition(CFuzzyVariable *var,CFuzzyTerm *term,bool not)
  {
//--- return fuzzy condition 
   return new CFuzzyCondition(var, term, not);
  }
//+------------------------------------------------------------------+
//| Create a single condition(3)                                     |
//+------------------------------------------------------------------+ 
CFuzzyCondition *CGenericFuzzyRule::CreateCondition(CFuzzyVariable *var,CFuzzyTerm *term,bool not,HedgeType hedge)
  {
//--- return fuzzy condition
   return new CFuzzyCondition(var, term, not, hedge);
  }
//+------------------------------------------------------------------+
//| Fuzzy rule for Mamdani fuzzy system.                             |
//| NOTE: a rule cannot be created directly, only via                |
//| MamdaniFuzzySystem::EmptyRule or MamdaniFuzzySystem::ParseRule   |
//+------------------------------------------------------------------+
class CMamdaniFuzzyRule : public CGenericFuzzyRule
  {
private:
   CSingleCondition *m_mamdani_conclusion;   // Mamdani conclusion
   double            m_weight;               // Weight of Mamdani rule

public:
                     CMamdaniFuzzyRule(void);
                    ~CMamdaniFuzzyRule(void);
   //--- methods gets or sets the conclusion (THEN) part of the rule
   CSingleCondition *Conclusion(void)                    { return(m_mamdani_conclusion); }
   void              Conclusion(CSingleCondition *value) { m_mamdani_conclusion=value;   }
   //--- methods gets or sets the rule weight 
   double            Weight(void)               { return(m_weight); }
   void              Weight(const double value) { m_weight=value;   }
   //--- method to check type
   virtual bool      IsTypeOf(EnRule type) { return(type==TYPE_CLASS_MamdaniFuzzyRule); }
  };
//+---------------------------------------------------------------+
//| Constructor without parameters                                |
//+---------------------------------------------------------------+
CMamdaniFuzzyRule::CMamdaniFuzzyRule(void)
  {
   m_weight=1.0;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+ 
CMamdaniFuzzyRule::~CMamdaniFuzzyRule(void)
  {
   if(CheckPointer(m_mamdani_conclusion)==POINTER_DYNAMIC)
      delete m_mamdani_conclusion;
  }
//+------------------------------------------------------------------+
//| Fuzzy rule for Sugeno fuzzy system                               |
//| NOTE: a rule cannot be created directly, only via                |
//| SugenoFuzzySystem::EmptyRule or SugenoFuzzySystem::ParseRule     |
//+------------------------------------------------------------------+
class CSugenoFuzzyRule : public CGenericFuzzyRule
  {
private:
   CSingleCondition *m_sugeno_conclusion; // Sugeno conclusion

public:
                     CSugenoFuzzyRule(void);
                    ~CSugenoFuzzyRule(void);
   //--- methods gets or sets the conclusion (THEN) part of the rule
   CSingleCondition *Conclusion(void)                    { return(m_sugeno_conclusion); }
   void              Conclusion(CSingleCondition *value) { m_sugeno_conclusion=value;   }
   //--- method to check type
   virtual bool      IsTypeOf(EnRule type) { return(type==TYPE_CLASS_SugenoFuzzyRule); }
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+-------------------------- ---------------------------------------+
CSugenoFuzzyRule::CSugenoFuzzyRule(void)
  {

  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+ 
CSugenoFuzzyRule::~CSugenoFuzzyRule(void)
  {
   if(CheckPointer(m_sugeno_conclusion)==POINTER_DYNAMIC)
      delete m_sugeno_conclusion;
  }
//+------------------------------------------------------------------+
