//+------------------------------------------------------------------+
//|                                           mandanifuzzysystem.mqh |
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
#include <Arrays\ArrayDouble.mqh>
#include "GenericFuzzySystem.mqh"
#include "InferenceMethod.mqh"
#include "RuleParser.mqh"
#include "FuzzyRule.mqh"
//+------------------------------------------------------------------+
//| Purpose: Creating Mamdani fuzzy system                           |
//+------------------------------------------------------------------+
//+------------------------------------------------------------------+
//| Mamdani fuzzy inference system                                   |
//+------------------------------------------------------------------+
class CMamdaniFuzzySystem : public CGenericFuzzySystem
  {
private:
   CList            *m_output;              // List of fuzzy variable
   CList            *m_rules;               // List of Mamdani fuzzy rule
   ImplicationMethod m_impl_method;         // Implication method
   AggregationMethod m_aggr_method;         // Aggregation method
   DefuzzificationMethod m_defuzz_method;   // Defuzzification method

public:
                     CMamdaniFuzzySystem(void);
                    ~CMamdaniFuzzySystem(void);
   //--- method gets the output linguistic variables
   CList             *Output() { return(m_output); }
   //--- method gets the fuzzy rule
   CList             *Rules() { return(m_rules);  }
   //--- methods gets or sets the implication method
   ImplicationMethod ImplicationMethod()                                { return (m_impl_method);   }
   void              ImplicationMethod(ImplicationMethod value)         { m_impl_method=value;      }
   //--- methods gets or sets the aggregation method
   AggregationMethod AggregationMethod()                                { return (m_aggr_method);   }
   void              AggregationMethod(AggregationMethod value)         { m_aggr_method=value;      }
   //--- methods gets or sets the defuzzification method
   DefuzzificationMethod DefuzzificationMethod()                        { return (m_defuzz_method); }
   void              DefuzzificationMethod(DefuzzificationMethod value) { m_defuzz_method=value;    }
   //--- maethod gets the variable by name
   CFuzzyVariable   *OutputByName(const string name);
   //--- create a new rule
   CMamdaniFuzzyRule *EmptyRule();
   //--- parse rule
   CMamdaniFuzzyRule *ParseRule(const string rule);
   //--- method for calculate result
   CList            *Calculate(CList *inputValues);
   CList            *EvaluateConditions(CList *fuzzifiedInput);
   CList            *Implicate(CList *conditions);
   CList            *Aggregate(CList *conclusions);
   CList            *Defuzzify(CList *fuzzyResult);
   double            Defuzzify(IMembershipFunction *mf,const double min,const double max);
  };
//+------------------------------------------------------------------+
//| Constructor without parameters                                   |
//+------------------------------------------------------------------+
CMamdaniFuzzySystem::CMamdaniFuzzySystem(void)
  {
   m_output= new CList;
   m_rules = new CList;
   m_impl_method = MinIpm;                // Implication method default is Min
   m_aggr_method =MaxAgg;                 // Aggregation method default is Max
   m_defuzz_method = CentroidDef;         // Defuzzification method default is Centroid   
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+     
CMamdaniFuzzySystem::~CMamdaniFuzzySystem(void)
  {
   delete m_output;
   delete m_rules;
  }
//+------------------------------------------------------------------+
//| Get output linguistic variable by its name                       |
//+------------------------------------------------------------------+
CFuzzyVariable *CMamdaniFuzzySystem::OutputByName(const string name)
  {
   for(int i=0; i<m_output.Total(); i++)
     {
      CFuzzyVariable *var=m_output.GetNodeAtIndex(i);
      if(var.Name()==name)
        {
         //--- return varriable
         return (var);
        }
     }
   Print("Variable with that name is not found");
//--- return 
   return (NULL);
  }
//+------------------------------------------------------------------+
//| Create new empty rule                                            |
//+------------------------------------------------------------------+
CMamdaniFuzzyRule *CMamdaniFuzzySystem::EmptyRule()
  {
//--- return empty rule
   return new CMamdaniFuzzyRule();
  }
//+------------------------------------------------------------------+
//| Parse rule from the string                                       |
//+------------------------------------------------------------------+
CMamdaniFuzzyRule *CMamdaniFuzzySystem::ParseRule(const string rule)
  {
//--- return Mamdani fuzzy rule
   return CRuleParser::Parse(rule, EmptyRule(), Input(), Output());
  }
//+------------------------------------------------------------------+
//| Calculate output values                                          |
//+------------------------------------------------------------------+
CList *CMamdaniFuzzySystem::Calculate(CList *inputValues)
  {
//--- There should be one rule as minimum
   if(m_rules.Total()==0)
     {
      Print("There should be one rule as minimum.");
      //--- return 
      return (NULL);
     }
//--- Fuzzification step
   CList *fuzzifiedInput=Fuzzify(inputValues);
//--- Evaluate the conditions
   CList *evaluatedConditions=EvaluateConditions(fuzzifiedInput);
//--- Do implication for each rule
   CList *implicatedConclusions=Implicate(evaluatedConditions);
//--- Aggrerate the results
   CList *fuzzyResult=Aggregate(implicatedConclusions);
//--- Defuzzify the result
   CList *result=Defuzzify(fuzzyResult);
//--- 
   delete fuzzyResult;
   for(int i=0; i<implicatedConclusions.Total(); i++)
     {
      CDictionary_Obj_Obj *pair=implicatedConclusions.GetNodeAtIndex(i);
      CCompositeMembershipFunction *composite=pair.Value();
      delete composite.MembershipFunctions().GetNodeAtIndex(0);
      delete composite;
     }
   delete implicatedConclusions;
   delete evaluatedConditions;
   for(int i=0; i<fuzzifiedInput.Total(); i++)
     {
      CDictionary_Obj_Obj *pair=fuzzifiedInput.GetNodeAtIndex(i);
      delete pair.Value();
     }
   delete fuzzifiedInput;
//--- return result
   return (result);
  }
//+------------------------------------------------------------------+
//| Evaluate conditions                                              |
//+------------------------------------------------------------------+
CList *CMamdaniFuzzySystem::EvaluateConditions(CList *fuzzifiedInput)
  {
   CList *result=new CList;
   for(int i=0; i<Rules().Total(); i++)
     {
      CDictionary_Obj_Double *p_rd=new CDictionary_Obj_Double;
      CMamdaniFuzzyRule *rule=Rules().GetNodeAtIndex(i);
      p_rd.SetAll(rule,EvaluateCondition(rule.Condition(),fuzzifiedInput));
      result.Add(p_rd);
     }
//--- return result
   return (result);
  }
//+------------------------------------------------------------------+
//| Implicate rule results                                           |
//+------------------------------------------------------------------+
CList *CMamdaniFuzzySystem::Implicate(CList *conditions)
  {
   CList *conclusions=new CList;
   for(int i=0; i<conditions.Total(); i++)
     {
      CDictionary_Obj_Double *p_rd=conditions.GetNodeAtIndex(i);
      CMamdaniFuzzyRule *rule=p_rd.Key();
      MfCompositionType compType;
      switch(m_impl_method)
        {
         case MinIpm :
           {
            compType=MinMF;
            break;
           }
         case ProductionImp :
           {
            compType=ProdMF;
            break;
           }
         default :
           {
            Print("Internal error.");
            //---
            return (NULL);
           }
        }
      CFuzzyTerm *val=rule.Conclusion().Term();
      IMembershipFunction *first_fun=new CConstantMembershipFunction(p_rd.Value());
      IMembershipFunction *second_fun=val.MembershipFunction();
      CCompositeMembershipFunction *resultMF=new CCompositeMembershipFunction(compType,first_fun,second_fun);
      CDictionary_Obj_Obj *p_rf=new CDictionary_Obj_Obj;
      p_rf.SetAll(rule,resultMF);
      conclusions.Add(p_rf);
     }
//--- return conclusions
   return (conclusions);
  }
//+------------------------------------------------------------------+
//| Aggregate results                                                |
//+------------------------------------------------------------------+
CList *CMamdaniFuzzySystem::Aggregate(CList *conclusions)
  {
   CList *fuzzyResult=new CList;
   for(int i=0; i<Output().Total(); i++)
     {
      CFuzzyVariable *var=Output().GetNodeAtIndex(i);
      CList *mfList=new CList;
      for(int j=0; j<conclusions.Total(); j++)
        {
         CDictionary_Obj_Obj *p_rf=conclusions.GetNodeAtIndex(j);
         CMamdaniFuzzyRule *rule=p_rf.Key();
         if(rule.Conclusion().Var()==var)
           {
            mfList.Add(p_rf.Value());
           }
        }
      MfCompositionType composType;
      switch(m_aggr_method)
        {
         case MaxAgg:
            composType=MaxMF;
            break;
         case SumAgg:
            composType=SumMF;
            break;
         default:
           {
            Print("Internal exception.");
            //--- return 
            return (NULL);
           }
        }
      CDictionary_Obj_Obj *p_vf=new CDictionary_Obj_Obj;
      CCompositeMembershipFunction *func=new CCompositeMembershipFunction(composType,mfList);
      p_vf.SetAll(var,func);
      fuzzyResult.Add(p_vf);
     }
//--- return  result
   return (fuzzyResult);
  }
//+------------------------------------------------------------------+
//| Calculate crisp result for each rule                             |
//+------------------------------------------------------------------+
CList *CMamdaniFuzzySystem::Defuzzify(CList *fuzzyResult)
  {
   CList *crispResult=new CList;
   for(int i=0; i<fuzzyResult.Total(); i++)
     {
      CDictionary_Obj_Double *p_vd=new CDictionary_Obj_Double;
      CDictionary_Obj_Obj *p_vf=fuzzyResult.GetNodeAtIndex(i);
      CFuzzyVariable *var=p_vf.Key();
      p_vd.SetAll(var,Defuzzify(p_vf.Value(),var.Min(),var.Max()));
      crispResult.Add(p_vd);
     }
//--- return  result
   return (crispResult);
  }
//+------------------------------------------------------------------+
//| Helpers                                                          |
//+------------------------------------------------------------------+
double CMamdaniFuzzySystem::Defuzzify(IMembershipFunction *mf,const double min,const double max)
  {
   if(m_defuzz_method==CentroidDef)
     {
      int k=50;                        // The function is divided into "k" steps 
      double step=(max-min)/k;         // Calculate the step function
      //+------------------------------------------------------------------+
      //| Calculate a center of gravity as integral                        |
      //+------------------------------------------------------------------+
      double ptLeft=0.0;
      double ptCenter= 0.0;
      double ptRight = 0.0;
      double valLeft=0.0;
      double valCenter= 0.0;
      double valRight = 0.0;
      double val2Left=0.0;
      double val2Center= 0.0;
      double val2Right = 0.0;
      double numerator=0.0;
      double denominator=0.0;
      for(int i=0; i<k; i++)
        {
         if(i==0)
           {
            ptRight=min;
            valRight=mf.GetValue(ptRight);
            val2Right=ptRight*valRight;
           }
         ptLeft=ptRight;
         ptCenter= min+step *((double)i+0.5);
         ptRight = min+step *(i+1);
         valLeft=valRight;
         valCenter= mf.GetValue(ptCenter);
         valRight = mf.GetValue(ptRight);
         val2Left=val2Right;
         val2Center= ptCenter * valCenter;
         val2Right = ptRight * valRight;
         numerator+=step *(val2Left+4*val2Center+val2Right)/3.0;
         denominator+=step *(valLeft+4*valCenter+valRight)/3.0;
        }
      delete mf;
      if(denominator!=0)
        {
         //--- return result
         return (numerator / denominator);
        }
      else
        {
         //--- return NAN
         return (MathLog(-1));
        }
     }
   else if(m_defuzz_method==BisectorDef)
     {
      //+-------------------------------------------------------------------------------------+
      //| The method Bisector consists in finding the point on the abscissa,                  |
      //| which divides the area under the curve of the membership function in two equal parts|
      //+-------------------------------------------------------------------------------------+   
      double Area=0.0;                 // The area under the function
      int k=50;                        // The function is divided into "k" steps 
      double now=min;                  // The current position
      for(int i=0; i<k; i++)
        {
         Area+=mf.GetValue(now);
         now=now+(max-min)/k;
        }
      now=min;
      double halfArea=fabs(Area/2-mf.GetValue(min));
      Area=0.0;
      while(true)
        {
         Area+=mf.GetValue(now);
         if(Area>=halfArea)
           {
            break;
           }
         now=now+(max-min)/k;
        }
      delete mf;
      //--- return result
      return (now);
     }
   else if(m_defuzz_method==AverageMaximumDef)
     {
      //+------------------------------------------------------------------------------------------+
      //| AverageMaximum method is the arithmetic mean of all the maxima of the membership function|
      //+------------------------------------------------------------------------------------------+    
      double sum_max=0;                // Sum of local maxima
      double count_max=0;              // Count of local maxima
      int k=50;                        // The function is divided into "k" steps 
      double now=min;                  // The current position
      double step=(max-min)/k;         // Calculate the step function
      for(int i=1; i<k; i++)
        {
         double point_1 = mf.GetValue(now);
         double point_0 = mf.GetValue(now - step);
         double point_2 = mf.GetValue(now + step);
         //--- check the first element
         if(i==1)
           {
            if(mf.GetValue(min)>mf.GetValue(min+step))
              {
               sum_max+=mf.GetValue(min);
               count_max++;
              }
           }
         //--- check the second element  
         if(i==k-1)
           {
            if(mf.GetValue(max)>mf.GetValue(max-step))
              {
               sum_max+=mf.GetValue(max);
               count_max++;
              }
           }
         //--- check all the other elements   
         if((point_1>point_0) && (point_1>point_2))
           {
            sum_max+=point_1;
            count_max++;
           }
        }
      if(count_max==0)
        {
         delete mf;
         //--- return result
         return (0);
        }
      else
        {
         delete mf;
         //--- return result
         return (sum_max/count_max);
        }
     }
   else if(m_defuzz_method==LargestMaximumDef)
     {
      CArrayDouble *local_max=new CArrayDouble; // Array of all local maximum
      double result;                   // Result of defuzzification method
      int k=50;                        // The function is divided into "k" steps 
      double now=min;                  // The current position
      double step=(max-min)/k;         // Calculate the step function
      for(int i=1; i<k; i++)
        {
         double point_1 = mf.GetValue(now);
         double point_0 = mf.GetValue(now - step);
         double point_2 = mf.GetValue(now + step);
         //--- check the first element
         if(i==1)
           {
            if(mf.GetValue(min)>mf.GetValue(min+step))
              {
               local_max.Add(mf.GetValue(min));
              }
           }
         //--- check the second element  
         if(i==k-1)
           {
            if(mf.GetValue(max)>mf.GetValue(max-step))
              {
               local_max.Add(mf.GetValue(max));
              }
           }
         //--- check all the other elements   
         if((point_1>point_0) && (point_1>point_2))
           {
            local_max.Add(point_1);
           }
         now+=step;
        }
      result=local_max.At(0);
      for(int i=0; i<local_max.Total(); i++)
        {
         if(result<=local_max.At(i))
           {
            result=local_max.At(i);
           }
        }
      now=min;
      while(true)
        {
         if(mf.GetValue(now)==result)
           {
            break;
           }
         now+=step;
        }
      delete local_max;
      delete mf;
      //--- return result
      return (now);
     }
   else if(m_defuzz_method==SmallestMaximumDef)
     {
      CArrayDouble *local_max=new CArrayDouble; // Array of all local maximum
      double result;                   // Result of defuzzification method
      int k=50;                        // The function is divided into "k" steps 
      double now=min;                  // The current position
      double step=(max-min)/k;         // Calculate the step function
      for(int i=1; i<k; i++)
        {
         double point_1 = mf.GetValue(now);
         double point_0 = mf.GetValue(now - step);
         double point_2 = mf.GetValue(now + step);
         //--- check the first element
         if(i==1)
           {
            if(mf.GetValue(min)>mf.GetValue(min+step))
              {
               local_max.Add(mf.GetValue(min));
              }
           }
         //--- check the second element  
         if(i==k-1)
           {
            if(mf.GetValue(max)>mf.GetValue(max-step))
              {
               local_max.Add(mf.GetValue(max));
              }
           }
         //--- check all the other elements   
         if((point_1>point_0) && (point_1>point_2))
           {
            local_max.Add(point_1);
           }
         now+=step;
        }
      result=local_max.At(0);
      for(int i=0; i<local_max.Total(); i++)
        {
         if(result>=local_max.At(i))
           {
            result=local_max.At(i);
           }
        }
      now=min;
      while(true)
        {
         if(mf.GetValue(now)==result)
           {
            break;
           }
         now+=step;
        }
      delete local_max;
      delete mf;
      //--- return result
      return (now);
     }
   else
     {
      Print("Internal exception.");
      delete mf;
      //--- return 
      return (0);
     }
  }
//+------------------------------------------------------------------+
