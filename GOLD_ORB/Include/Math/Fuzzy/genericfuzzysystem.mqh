//+------------------------------------------------------------------+
//|                                           genericfuzzysystem.mqh |
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
#include <Arrays\ArrayObj.mqh>
#include "FuzzyRule.mqh"
#include "InferenceMethod.mqh"
#include "Dictionary.mqh"
//+------------------------------------------------------------------+
//| Purpose: Creating generic fuzzy system                           |
//+------------------------------------------------------------------+
//+------------------------------------------------------------------+
//| Common functionality of Mamdani and Sugeno fuzzy systems         |
//+------------------------------------------------------------------+
class CGenericFuzzySystem
  {
private:
   CList            *m_input;          // List of input fuzzy variables
   AndMethod         m_and_method;     // And method from InferenceMethod
   OrMethod          m_or_method;      // Or method from InferenceMethod

protected:
                     CGenericFuzzySystem(void);
                    ~CGenericFuzzySystem(void);
public:
   //--- method gets the input linguistic variables
   CList            *Input(void) { return(m_input); }
   //--- method gets or sets the type of "And method"
   void              AndMethod(AndMethod value) { m_and_method=value;    }
   AndMethod         AndMethod(void)            { return (m_and_method); }
   //--- method gets or sets the type of "Or method"
   void              OrMethod(OrMethod value)   { m_or_method=value;     }
   OrMethod          OrMethod(void)             { return (m_or_method);  }
   //--- method gets the varriable by name
   CFuzzyVariable *InputByName(const string name);
   //--- common steps of calculating
   CList *Fuzzify(CList *inputValues);
protected:
   double EvaluateCondition(ICondition *condition,CList *fuzzifiedInput);
   double EvaluateConditionPair(const double cond1,const double cond2,OperatorType op);
private:
   bool ValidateInputValues(CList *inputValues,string &msg);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CGenericFuzzySystem::CGenericFuzzySystem(void)
  {
   m_input=new CList;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+     
CGenericFuzzySystem::~CGenericFuzzySystem(void)
  {
   if(CheckPointer(m_input)==POINTER_DYNAMIC)
     {
      delete m_input;
     }
  }
//+------------------------------------------------------------------+
//| Get input linguistic variable by its name                        |
//+------------------------------------------------------------------+ 
CFuzzyVariable *CGenericFuzzySystem::InputByName(const string name)
  {
   CList *result=CGenericFuzzySystem::Input();
   for(int i=0; i<result.Total(); i++)
     {
      CFuzzyVariable *var=result.GetNodeAtIndex(i);
      if(var.Name()==name)
        {
         //--- return fuzzy variable 
         return (var);
        }
     }
   Print("The variable with that name is not found");
//--- return 
   return (NULL);
  }
//+------------------------------------------------------------------+
//| Fuzzify input                                                    |
//+------------------------------------------------------------------+ 
CList *CGenericFuzzySystem::Fuzzify(CList *inputValues)
  {
//--- Validate input
   string msg;
   if(!ValidateInputValues(inputValues,msg))
     {
      Print(msg);
      //--- return 
      return (NULL);
     }
//--- Fill results list
   CList *result=new CList;
   for(int i=0; i<Input().Total(); i++)
     {
      CFuzzyVariable *var=Input().GetNodeAtIndex(i);
      double value=NULL;
      for(int k=0; k<inputValues.Total(); k++)
        {
         CDictionary_Obj_Double *p_vd=inputValues.GetNodeAtIndex(i);
         CFuzzyVariable *v=p_vd.Key();
         if(p_vd.Key()==var)
           {
            value=p_vd.Value();
            break;
           }
        }
      CList *resultForVar=new CList;
      for(int j=0; j<var.Terms().Total(); j++)
        {
         CDictionary_Obj_Double *p_vd=new CDictionary_Obj_Double;
         CFuzzyTerm *term=var.Terms().GetNodeAtIndex(j);
         p_vd.SetAll(term,term.MembershipFunction().GetValue(value));
         resultForVar.Add(p_vd);
        }
      CDictionary_Obj_Obj *p_vl=new CDictionary_Obj_Obj;
      p_vl.SetAll(var,resultForVar);
      result.Add(p_vl);
     }
//--- return result
   return (result);
  }
//+------------------------------------------------------------------+
//| Evaluate fuzzy condition (or conditions)                         |
//+------------------------------------------------------------------+ 
double CGenericFuzzySystem::EvaluateCondition(ICondition *condition,CList *fuzzifiedInput)
  {
   double result=0.0;
   ICondition *IC;
   if(condition.IsTypeOf(TYPE_CLASS_Conditions))
     {
      CConditions *conds=condition;
      if(conds.ConditionsList().Total()==0)
        {
         Print("Inner exception.");
        }
      else if(conds.ConditionsList().Total()==1)
        {
         IC=conds.ConditionsList().GetNodeAtIndex(0);
         result=EvaluateCondition(IC,fuzzifiedInput);
        }
      else
        {
         IC=conds.ConditionsList().GetNodeAtIndex(0);
         result=EvaluateCondition(IC,fuzzifiedInput);
         for(int i=1; i<conds.ConditionsList().Total(); i++)
           {
            IC=conds.ConditionsList().GetNodeAtIndex(i);
            double cond2=EvaluateCondition(IC,fuzzifiedInput);;
            result=EvaluateConditionPair(result,cond2,conds.Op());
           }
        }
      if(conds.Not())
        {
         result=1.0-result;
        }
      //--- return result 
      return (result);
     }
   else if(condition.IsTypeOf(TYPE_CLASS_FuzzyCondition))
     {
      CFuzzyCondition *cond=condition;
      CDictionary_Obj_Obj *p_vl;
      CDictionary_Obj_Double *p_td;
      for(int i=0; i<fuzzifiedInput.Total(); i++)
        {
         p_vl=fuzzifiedInput.GetNodeAtIndex(i);
         if(p_vl.Key()==cond.Var())
           {
            CList *list=p_vl.Value();
            for(int j=0; j<list.Total(); j++)
              {
               p_td=list.GetNodeAtIndex(j);
               if(p_td.Key()==cond.Term())
                 {
                  break;
                 }
              }
            break;
           }
        }
      result=p_td.Value();
      switch(cond.Hedge())
        {
         case Slightly:
            //--- Cube root
            result=pow(result,1.0/3.0);
            break;
         case Somewhat:
            result=sqrt(result);
            break;
         case Very:
            result=result*result;
            break;
         case Extremely:
            result=result*result*result;
            break;
         default:
            break;
        }
      if(cond.Not())
        {
         result=1.0-result;
        }
      //--- return result 
      return (result);
     }
   else
     {
      Print("Internal error.");
      //--- return 
      return (NULL);
     }
  }
//+------------------------------------------------------------------+
//| Evaluate fuzzy condition (or conditions)                         |
//+------------------------------------------------------------------+ 
double CGenericFuzzySystem::EvaluateConditionPair(const double cond1,const double cond2,OperatorType op)
  {
   if(op==And)
     {
      if(CGenericFuzzySystem::AndMethod()==MinAnd)
        {
         //--- return evaluate condition
         return fmin(cond1, cond2);
        }
      else if(CGenericFuzzySystem::AndMethod()==ProductionAnd)
        {
         //--- return evaluate condition 
         return (cond1 * cond2);
        }
      else
        {
         Print("Internal error.");
         //--- return  
         return(NULL);
        }
     }
   else if(op==Or)
     {
      if(CGenericFuzzySystem::OrMethod()==MaxOr)
        {
         //--- return evaluate condition 
         return fmax(cond1, cond2);
        }
      else if(CGenericFuzzySystem::OrMethod()==ProbabilisticOr)
        {
         //--- return evaluate condition 
         return (cond1 + cond2 - cond1 * cond2);
        }
      else
        {
         Print("Internal error.");
         //--- return 
         return (NULL);
        }
     }
   else
     {
      Print("Internal error.");
      //--- return 
      return (NULL);
     }
  }
//+------------------------------------------------------------------+
//|  Validate input values                                           |
//+------------------------------------------------------------------+
bool CGenericFuzzySystem::ValidateInputValues(CList *inputValues,string &msg)
  {
   msg=NULL;
   if(inputValues.Total()!=Input().Total())
     {
      msg="Input values count is incorrect.";
      //--- return false
      return (false);
     }
   bool contain;
   for(int i=0; i<Input().Total(); i++)
     {
      CFuzzyVariable *var=Input().GetNodeAtIndex(i);
      contain=false;
      for(int j=0; j<inputValues.Total();j++)
        {
         CDictionary_Obj_Double *p_vd=inputValues.GetNodeAtIndex(j);
         if(p_vd.Key()==var)
           {
            contain=true;
            double val=p_vd.Value();
            if(val<var.Min() || val>var.Max())
              {
               msg=StringFormat("Value for the %s variable is out of range.",var.Name());
               //--- return false
               return (false);
              }
           }
        }
      if(contain==false)
        {
         msg=StringFormat("Value for the %s variable does not present.",var.Name());
         //--- return false
         return (false);
        }
     }
//--- return true 
   return (true);
  }
//+------------------------------------------------------------------+
