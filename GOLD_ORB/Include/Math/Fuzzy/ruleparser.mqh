//+------------------------------------------------------------------+
//|                                                   ruleparser.mqh |
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
#include "Dictionary.mqh"
#include "FuzzyRule.mqh"
#include "Helper.mqh"
#include "InferenceMethod.mqh"
#include "SugenoVariable.mqh"
//+------------------------------------------------------------------+
//| Purpose: Analysis of the fuzzy rules                             |
//+------------------------------------------------------------------+
//+------------------------------------------------------------------+
//| Base class for all expression                                    |
//+------------------------------------------------------------------+
class IExpression : public CObject
  {
public:
   //--- method gets the text of expression
   virtual string    Text(void)=NULL;
   //--- method to check type
   virtual bool      IsTypeOf(EnLexem type) { return(type==TYPE_CLASS_IExpression); }
  };
//+------------------------------------------------------------------+
//| Class for creating lexem                                         |
//+------------------------------------------------------------------+
class CLexem : public IExpression
  {
public:
   //--- method gets the text of expression
   virtual string    Text(void)=NULL;
   //--- method to check type
   virtual bool      IsTypeOf(EnLexem type) { return(type==TYPE_CLASS_Lexem); }
  };
//+------------------------------------------------------------------+
//| Class condition expression                                       |
//+------------------------------------------------------------------+     
class CConditionExpression : public IExpression
  {
private:
   CArrayObj        *m_expressions;    // List of expression
   CFuzzyCondition *m_condition;       // Fuzzy condition

public:
                     CConditionExpression(CArrayObj *expressions,CFuzzyCondition *condition);
                    ~CConditionExpression(void);
   //--- methods gets or sets array of expression                    
   CArrayObj        *Expressions(void)             { return (m_expressions); }
   void              Expressions(CArrayObj *value) { m_expressions=value;    }
   //--- methods gets or sets fuzzy condition
   CFuzzyCondition *Condition(void)                { return m_condition;     }
   void Condition(CFuzzyCondition *value)          { m_condition=value;      }
   //--- method gets the text of expressions
   string            Text(void);
   //--- method to check type  
   virtual bool      IsTypeOf(EnLexem type) { return(type==TYPE_CLASS_ConditionExpression); }
  };
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+
CConditionExpression::CConditionExpression(CArrayObj *expressions,CFuzzyCondition *condition)
  {
   m_expressions=expressions;
   m_condition=condition;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CConditionExpression::~CConditionExpression(void)
  {

  }
//+------------------------------------------------------------------+
//| Convert all expressions to text(string)                          |
//+------------------------------------------------------------------+  
string CConditionExpression::Text(void)
  {
   string sb;
   for(int i=0; i<m_expressions.Total(); i++)
     {
      IExpression *ex=m_expressions.At(i);
      string str=ex.Text();
      StringAdd(sb,str);
     }
//--- return text of all expressions
   return ((string)sb);
  }
//+------------------------------------------------------------------+
//| Class keyword lexem                                              |
//+------------------------------------------------------------------+       
class CKeywordLexem : public CLexem
  {
private:
   string            m_name;          // Name of keyword lexem

public:
                     CKeywordLexem(const string name);
                    ~CKeywordLexem(void);
   //--- method gets the text of lexem     
   string            Text(void) { return(m_name); }
   //--- method to check type   
   virtual bool      IsTypeOf(EnLexem type) { return(type==TYPE_CLASS_KeywordLexem); }
  };
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+   
CKeywordLexem::CKeywordLexem(const string name)
  {
   m_name=name;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CKeywordLexem::~CKeywordLexem(void)
  {

  }
//+------------------------------------------------------------------+
//| Class for handling variable lexem                                |
//+------------------------------------------------------------------+      
class CVarLexem : public CLexem
  {
private:
   INamedVariable   *m_var;           // Variable type
   bool              m_input;         // Confirmation of the variable

public:
                     CVarLexem(INamedVariable *var,bool in);
                    ~CVarLexem(void);
   //--- methods gets or sets the varriable                         
   INamedVariable   *Var(void)                  { return (m_var); }
   void              Var(INamedVariable *var)   { m_var=var;      }
   //--- method gets the text of lexem         
   string            Text(void) { return(m_var.Name()); }
   //--- methods gets or sets mark is input varriable          
   bool              Input(void)       { return (m_input); }
   void              Input(bool value) { m_input=value;    }
   //--- method to check type  
   virtual bool      IsTypeOf(EnLexem type) { return(type==TYPE_CLASS_VarLexem); }
  };
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+   
CVarLexem::CVarLexem(INamedVariable *var,bool in)
  {
   m_var=var;
   m_input=in;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+     
CVarLexem::~CVarLexem(void)
  {

  }
//+------------------------------------------------------------------+
//| Class alternative lexem                                          |
//+------------------------------------------------------------------+       
class IAltLexem : public CLexem
  {
public:
   //--- methods gets or sets alternative lexem
   virtual IAltLexem *Alternative(void)=NULL;
   virtual void      Alternative(IAltLexem *value)=NULL;
   //--- method to check type  
   virtual bool      IsTypeOf(EnLexem type) { return(type==TYPE_CLASS_AltLexem); }
  };
//+------------------------------------------------------------------+
//| Class term lexem                                                 |
//+------------------------------------------------------------------+
class CTermLexem : public IAltLexem
  {
private:
   INamedValue      *m_term;          // Value type
   IAltLexem        *m_alternative;   // Alternative lexem
   bool              m_input;

public:
                     CTermLexem(INamedValue *term,bool in);
                    ~CTermLexem(void);
   //--- methods gets or sets the term        
   INamedValue      *Term(void)                 { return(m_term); }
   void              Term(INamedValue *value)   { m_term=value;   }
   //--- method gets the text of lexem 
   string            Text(void) { return m_term.Name(); }
   //--- methods gets or sets the alternative lexem
   IAltLexem        *Alternative(void)             { return(m_alternative); }
   void              Alternative(IAltLexem *value) { m_alternative=value;   }
   //--- method to check type  
   virtual bool      IsTypeOf(EnLexem type) { return(type==TYPE_CLASS_TermLexem); }
  };
//+------------------------------------------------------------------+
//| Constructor with parameters                                      |
//+------------------------------------------------------------------+      
CTermLexem::CTermLexem(INamedValue *term,bool in)
  {
   m_term=term;
   m_input=in;
  }
//+------------------------------------------------------------------+
//| Destructor                                                       |
//+------------------------------------------------------------------+
CTermLexem::~CTermLexem(void)
  {
   if(CheckPointer(m_alternative)==POINTER_DYNAMIC)
      delete m_alternative;
  }
//+------------------------------------------------------------------+
//| Class responsible for parsing                                    |
//+------------------------------------------------------------------+
class CRuleParser : public INamedVariable
  {
public:
   //--- methods of parsing
   static IParsableRule *Parse(const string rule,IParsableRule *emptyRule,CList *in,CList *out);
private:
   static CList     *BuildLexemsList(CList *in,CList *out);
   static void       BuildLexemsList(INamedVariable *var,bool in,CList *lexems);
   static CArrayObj *ParseLexems(const string rule,CList *lexems);
   static CArrayObj *ExtractSingleCondidtions(CArrayObj *conditionExpression,CList *in,CList *lexems);
   static CConditions *ParseConditions(CArrayObj *conditionExpression,CList *in,CList *lexems);
   static int        FindPairBracket(CArrayObj *expressions,CList *lexems);
   static ICondition *ParseConditionsRecurse(CArrayObj *expressions,CList *lexems);
   static CSingleCondition *ParseConclusion(CArrayObj *conditionExpression,CList *out,CList *lexems);
  };
//+------------------------------------------------------------------+
//| Parse                                                            |
//+------------------------------------------------------------------+
static IParsableRule *CRuleParser::Parse(const string rule,IParsableRule *emptyRule,CList *in,CList *out)
  {
   if(StringLen(rule)==0)
     {
      Print("Rule cannot be empty.");
      //--- return 
      return (NULL);
     }
//--- Surround brakes with spaces, remove double spaces
   string sb=NULL;
   char ch;
   for(int i=0; i<StringLen(rule); i++)
     {
      ch=(char)StringGetCharacter(rule,i);
      if((ch=='(') || (ch==')'))
        {
         if(StringLen(sb)>0 && StringGetCharacter(sb,StringLen(sb)-1)==' ')
           {
            //--- Do not duplicate spaces
           }
         else
           {
            sb+=CharToString(' ');
           }
         sb+=CharToString(ch);
         sb+=CharToString(' ');
        }
      else
        {
         if(ch==' ' && StringLen(sb)>0 && StringGetCharacter(sb,StringLen(sb)-1)==' ')
           {
            // Do not duplicate spaces
           }
         else
           {
            sb+=CharToString(ch);
           }
        }
     }
//--- Remove spaces
//+------------------------------------------------------------------+
//| Use conditional compilation to determine                         |
//| the type of program MQL4 or MQL5 because they have               |
//| different realization of StringTrimRight() and StringTrimLeft()  |
//+------------------------------------------------------------------+      
#ifdef __MQL5__
   StringTrimRight(sb);
   StringTrimLeft(sb);
#else 
#ifdef  __MQL4__
   sb=StringTrimRight(sb);
   sb=StringTrimLeft(sb);
#endif
#endif
   string prepRule=sb;
//--- Build lexems dictionary
   CList *lexemsDict=BuildLexemsList(in,out);
//--- At first we parse lexems
   CArrayObj *expressions=ParseLexems(prepRule,lexemsDict);
   if(expressions.Total()==0)
     {
      Print("No valid identifiers found.");
      //--- return 
      return (NULL);
     }
//--- Find condition & conclusion parts part  
   CDictionary_String_Obj *p_so;
   for(int i=0; i<lexemsDict.Total(); i++)
     {
      p_so=lexemsDict.GetNodeAtIndex(i);
      if(p_so.Key()=="if")
        {
         break;
        }
     }
   if(expressions.At(0)!=p_so.Value())
     {
      Print("'if' should be the first identifier.");
      //--- return 
      return (NULL);
     }
   int thenIndex=-1;
   for(int i=1; i<expressions.Total(); i++)
     {
      for(int j=0; j<lexemsDict.Total(); j++)
        {
         p_so=lexemsDict.GetNodeAtIndex(j);
         if(p_so.Key()=="then"){break;}
        }
      if(expressions.At(i)==p_so.Value())
        {
         thenIndex=i;
         break;
        }
     }
   if(thenIndex==-1)
     {
      Print("'then' identifier not found.");
      //--- return 
      return (NULL);
     }
   int conditionLen=thenIndex-1;
   if(conditionLen<1)
     {
      Print("Condition part of the rule not found.");
      //--- return 
      return (NULL);
     }
   int conclusionLen=expressions.Total()-thenIndex-1;
   if(conclusionLen<1)
     {
      Print("Conclusion part of the rule not found.");
      //--- return 
      return (NULL);
     }
   CArrayObj *conditionExpressions=GetRange(expressions,1,conditionLen);
   CArrayObj *conclusionExpressions=GetRange(expressions,thenIndex+1,conclusionLen);
   CConditions *conditions=ParseConditions(conditionExpressions,in,lexemsDict);
   delete emptyRule.Condition();
   emptyRule.Condition(conditions);
   CSingleCondition *conclusion=ParseConclusion(conclusionExpressions,out,lexemsDict);
   if(emptyRule.IsTypeOf(TYPE_CLASS_MamdaniFuzzyRule))
     {
      CMamdaniFuzzyRule *mamdani_rule=emptyRule;
      mamdani_rule.Conclusion(conclusion);
      emptyRule=mamdani_rule;
     }
   if(emptyRule.IsTypeOf(TYPE_CLASS_SugenoFuzzyRule))
     {
      CSugenoFuzzyRule *sugeno_rule=emptyRule;
      sugeno_rule.Conclusion(conclusion);
      emptyRule=sugeno_rule;
     }
//--- return empty rule
   conditionExpressions.FreeMode(false);
   delete conditionExpressions;
   conclusionExpressions.FreeMode(false);
   delete conclusionExpressions;
   expressions.FreeMode(false);
   delete expressions;
   delete lexemsDict;
   return (emptyRule);
  }
//+------------------------------------------------------------------+
//| Build lexems list                                                |
//+------------------------------------------------------------------+
static CList *CRuleParser::BuildLexemsList(CList *in,CList *out)
  {
   CList *lexems=new CList();
   for(int i=0;i<ArraySize(KEYWORDS); i++)
     {
      string keyword=KEYWORDS[i];
      CKeywordLexem *keywordLexem=new CKeywordLexem(keyword);
      CDictionary_String_Obj *p_so=new CDictionary_String_Obj;
      p_so.SetAll(keywordLexem.Text(),keywordLexem);
      lexems.Add(p_so);
     }
   for(int i=0; i<in.Total(); i++)
     {
      INamedVariable *var=in.GetNodeAtIndex(i);
      BuildLexemsList(var,true,lexems);
     }
   for(int i=0; i<out.Total(); i++)
     {
      INamedVariable *var=out.GetNodeAtIndex(i);
      BuildLexemsList(var,false,lexems);
     }
//--- return lexems
   return (lexems);
  }
//+------------------------------------------------------------------+
//| Build lexems list                                                |
//+------------------------------------------------------------------+
static void CRuleParser::BuildLexemsList(INamedVariable *var,bool in,CList *lexems)
  {
   CVarLexem *varLexem=new CVarLexem(var,in);
   CDictionary_String_Obj *p_so_var=new CDictionary_String_Obj;
   p_so_var.SetAll(varLexem.Text(),varLexem);
   lexems.Add(p_so_var);
   for(int i=0; i<var.Values().Total(); i++)
     {
      INamedValue *term=var.Values().GetNodeAtIndex(i);
      CTermLexem *termLexem=new CTermLexem(term,in);
      CLexem *foundLexem;
      bool contain=false;
      for(int j=0; j<lexems.Total(); j++)
        {
         CDictionary_String_Obj *p_so=lexems.GetNodeAtIndex(j);
         foundLexem=p_so.Value();
         if(p_so.Key()==termLexem.Text())
           {
            contain=true;
            break;
           }
        }
      if(contain==false)
        {
         CDictionary_String_Obj *p_so_val=new CDictionary_String_Obj;
         p_so_val.SetAll(termLexem.Text(),termLexem);
         lexems.Add(p_so_val);
        }
      else
        {
         if(foundLexem.IsTypeOf(TYPE_CLASS_TermLexem))
           {
            //--- There can be more than one terms with the same name.
            //--- TODO: But only if they belong to defferent variables.
            CTermLexem *foundTermLexem=foundLexem;
            while(foundTermLexem.Alternative()!=NULL)
              {
               foundTermLexem=foundTermLexem.Alternative();
              }
            foundTermLexem.Alternative(termLexem);
           }
         else
           {
            //--- Only terms of different vatiables can have the same name
            Print("Found more than one lexems with the same name");
           }
        }
     }
  }
//+------------------------------------------------------------------+
//| Parse lexems                                                     |
//+------------------------------------------------------------------+
static CArrayObj *CRuleParser::ParseLexems(const string rule,CList *lexems)
  {
   CArrayObj *expressions=new CArrayObj;
   string words[];
   StringSplit(rule,' ',words);
   int index=0;
   for(int i=0; i<ArraySize(words); i++)
     {
      string word=words[i];
      CObject *lexem=NULL;
      if(TryGetValue(lexems,word,lexem))
        {
         expressions.Add(lexem);
        }
      else
        {
         Print(StringFormat("Unknown identifier : %s",word));
         //--- return  NULL
         return(NULL);
        }
     }
//--- return expressions
   return(expressions);
  }
//+------------------------------------------------------------------+
//| Extract single condidtions                                       |
//+------------------------------------------------------------------+
static CArrayObj *CRuleParser::ExtractSingleCondidtions(CArrayObj *conditionExpression,CList *in,CList *lexems)
  {
   CArrayObj *copyExpressions=conditionExpression;
   CArrayObj *expressions=new CArrayObj;
   int index=0;
   while(copyExpressions.Total()-index>0)
     {
      IExpression *expr0=copyExpressions.At(index);
      if(expr0.IsTypeOf(TYPE_CLASS_VarLexem))
        {
         //--- Parse variable
         CVarLexem *varLexem=copyExpressions.At(index);
         if(copyExpressions.Total()<3)
           {
            Print(StringFormat("Condition strated with '%s' is incorrect.",varLexem.Text()));
            //--- return 
            return(NULL);
           }
         if(varLexem.Input()==false)
           {
            Print("The variable in condition part must be an input variable.");
            //--- return 
            return(NULL);
           }
         //--- Parse 'is' lexem
         CLexem *exprIs=copyExpressions.At(index+1);
         CDictionary_String_Obj *p_so;
         for(int i=0;i<lexems.Total();i++)
           {
            p_so=lexems.GetNodeAtIndex(i);
            if(p_so.Key()=="is")
              {
               break;
              }
           }
         if(exprIs!=p_so.Value())
           {
            Print(StringFormat("'is' keyword must go after '%s' identifier.",varLexem.Text()));
            //--- return 
            return(NULL);
           }
         //--- Parse 'not' lexem (if exists)
         int cur=2;
         bool not=false;
         for(int i=0;i<lexems.Total();i++)
           {
            p_so=lexems.GetNodeAtIndex(i);
            if(p_so.Key()=="not")
              {
               break;
              }
           }
         if(copyExpressions.At(cur+index)==p_so.Value())
           {
            not=true;
            cur++;
            if(copyExpressions.Total()-index<=cur)
              {
               Print("Error at 'not' in condition part of the rule.");
               //--- return 
               return(NULL);
              }
           }
         //--- "slightly"
         //--- "somewhat"
         //--- "very"
         //--- "extremely"
         //--- Parse hedge modifier (if exists)
         HedgeType hedge=None;
         for(int i=0;i<lexems.Total();i++)
           {
            p_so=lexems.GetNodeAtIndex(i);
            if(p_so.Key()=="slightly")
              {
               if(copyExpressions.At(cur+index)==p_so.Value())
                 {
                  hedge=Slightly;
                  break;
                 }
              }
            else if(p_so.Key()=="somewhat")
              {
               if(copyExpressions.At(cur+index)==p_so.Value())
                 {
                  hedge=Somewhat;
                  break;
                 }
              }
            else if(p_so.Key()=="very")
              {
               if(copyExpressions.At(cur+index)==p_so.Value())
                 {
                  hedge=Very;
                  break;
                 }
              }
            else if(p_so.Key()=="extremely")
              {
               if(copyExpressions.At(cur+index)==p_so.Value())
                 {
                  hedge=Extremely;
                  break;
                 }
              }
           }
         if(hedge!=None)
           {
            cur++;
            if(copyExpressions.Total()-index<=cur)
              {
               Print("Error in condition part of the rule.");
               //--- return 
               return(NULL);
              }
           }
         //--- Parse term
         CLexem *exprTerm=copyExpressions.At(cur+index);
         if(!exprTerm.IsTypeOf(TYPE_CLASS_TermLexem))
           {
            Print(StringFormat("Wrong identifier '%s' in conditional part of the rule.",exprTerm.Text()));
            //--- return 
            return(NULL);
           }
         IAltLexem *altLexem=exprTerm;
         CTermLexem *termLexem=NULL;
         do
           {
            if(!altLexem.IsTypeOf(TYPE_CLASS_TermLexem))
              {
               continue;
              }
            termLexem=altLexem;
            if(varLexem.Var().Values().IndexOf(termLexem.Term())==-1)
              {
               termLexem=NULL;
               continue;
              }
           }
         while((altLexem=altLexem.Alternative())!=NULL && termLexem==NULL);
         if(termLexem==NULL)
           {
            Print(StringFormat("Wrong identifier '%s' in conditional part of the rule.",exprTerm.Text()));
            //--- return 
            return(NULL);
           }
         //--- Add new condition expression
         CFuzzyCondition *condition=new CFuzzyCondition(varLexem.Var(),termLexem.Term(),not,hedge);
         CArrayObj *cutExpression=GetRange(copyExpressions,index,cur+1);
         expressions.Add(new CConditionExpression(cutExpression,condition));
         cutExpression.FreeMode(false);
         delete cutExpression;
         index=index+cur+1;
        }
      else
        {
         CDictionary_String_Obj *p_so;
         IExpression *expr=copyExpressions.At(index);
         bool contain=false;
         for(int i=0; i<lexems.Total(); i++)
           {
            p_so=lexems.GetNodeAtIndex(i);
            if(p_so.Key()=="and" || p_so.Key()=="or" || p_so.Key()=="(" || p_so.Key()==")")
              {
               if(expr==p_so.Value())
                 {
                  contain=true;
                  break;
                 }
              }
           }
         if(contain==true)
           {
            expressions.Add(copyExpressions.At(index));
            index=index+1;
           }
         else
           {
            CLexem *unknownLexem=expr;
            Print(StringFormat("CLexem '%s' found at the wrong place in condition part of the rule.",unknownLexem.Text()));
            //--- return 
            return (NULL);
           }
        }
     }
//--- return expressions
   return (expressions);
  }
//+------------------------------------------------------------------+
//| Parse conditions                                                 |
//+------------------------------------------------------------------+
static CConditions *CRuleParser::ParseConditions(CArrayObj *conditionExpression,CList *in,CList *lexems)
  {
//--- Extract single conditions
   CArrayObj *expressions=ExtractSingleCondidtions(conditionExpression,in,lexems);
   if(expressions.Total()==0)
     {
      Print("No valid conditions found in conditions part of the rule.");
      //--- return 
      return (NULL);
     }
   ICondition *cond=ParseConditionsRecurse(expressions,lexems);
   delete expressions;
   if(cond.IsTypeOf(TYPE_CLASS_Conditions))
     {
      //--- return conditions
      return (cond);
     }
   else
     {
      delete cond;
      delete expressions;
      CConditions *conditions=new CConditions();
      //--- return conditions
      return (conditions);
     }
  }
//+------------------------------------------------------------------+
//| Find pair bracket                                                |
//+------------------------------------------------------------------+
static int CRuleParser::FindPairBracket(CArrayObj *expressions,CList *lexems)
  {
//--- Assume that '(' stands at first place
   int bracketsOpened=1;
   int closeBracket=-1;
   CDictionary_String_Obj *p_so_open;
   CDictionary_String_Obj *p_so_close;
   for(int i=0; i<lexems.Total(); i++)
     {
      CDictionary_String_Obj *p_so=lexems.GetNodeAtIndex(i);
      if(p_so.Key()=="(")
        {
         p_so_open=p_so;
        }
      if(p_so.Key()==")")
        {
         p_so_close=p_so;
        }
     }
   for(int i=1; i<expressions.Total(); i++)
     {
      if(expressions.At(i)==p_so_open.Value())
        {
         bracketsOpened++;
        }
      else if(expressions.At(i)==p_so_close.Value())
        {
         bracketsOpened--;
         if(bracketsOpened==0)
           {
            closeBracket=i;
            break;
           }
        }
     }
//--- return index of bracket
   return (closeBracket);
  }
//+------------------------------------------------------------------+
//| Parse conditions recurse                                         |
//+------------------------------------------------------------------+
static ICondition *CRuleParser::ParseConditionsRecurse(CArrayObj *expressions,CList *lexems)
  {
   if(expressions.Total()<1)
     {
      Print("Empty condition found.");
      //--- return 
      return(NULL);
     }
   CDictionary_String_Obj *p_so;
   for(int i=0; i<lexems.Total();i++)
     {
      p_so=lexems.GetNodeAtIndex(i);
      if(p_so.Key()=="(")
        {
         break;
        }
     }
   IExpression *expr=expressions.At(0);
   if(expressions.At(0)==p_so.Value() && FindPairBracket(expressions,lexems)==expressions.Total())
     {
      //--- Remove extra brackets
      //-- return  conditions
      CArrayObj *cutExpression=GetRange(expressions,1,expressions.Total()-2);
      ICondition *cond=ParseConditionsRecurse(cutExpression,lexems);
      cutExpression.FreeMode(false);
      delete cutExpression;
      return cond;
     }
   else if(expressions.Total()==1 && expr.IsTypeOf(TYPE_CLASS_ConditionExpression))
     {
      //-- return single conditions
      CConditionExpression *condExp=expressions.At(0);
      return condExp.Condition();
     }
   else
     {
      //--- Parse list of one level conditions connected by or/and
      CArrayObj *copyExpressions=expressions;
      CConditionExpression *condExp;
      CConditions *conds=new CConditions();
      int index=0;
      bool setOrAnd=false;
      while(copyExpressions.Total()-index>0)
        {
         ICondition *cond=NULL;
         for(int i=0; i<lexems.Total();i++)
           {
            p_so=lexems.GetNodeAtIndex(i);
            if(p_so.Key()=="(")
              {
               break;
              }
           }
         if(copyExpressions.At(0)==p_so.Value())
           {
            //--- Find pair bracket
            int closeBracket=FindPairBracket(copyExpressions,lexems);
            if(closeBracket==-1)
              {
               Print("Parenthesis error.");
               //--- return 
               return (NULL);
              }
            CArrayObj *cutExpression=GetRange(copyExpressions,index+1,closeBracket-1);
            cond=ParseConditionsRecurse(cutExpression,lexems);
            cutExpression.FreeMode(false);
            delete cutExpression;
            index=index+closeBracket+1;
           }
         else if(expr.IsTypeOf(TYPE_CLASS_ConditionExpression))
           {
            condExp=copyExpressions.At(index);
            cond=condExp.Condition();
            index=index+1;
           }
         else
           {
            condExp=copyExpressions.At(index);
            Print(StringFormat("Wrong expression in condition part at '%s'",condExp.Text()));
            //--- return 
            return (NULL);
           }
         //--- And condition to the list
         conds.ConditionsList().Add(cond);
         p_so=NULL;
         CDictionary_String_Obj *p_so_and;
         CDictionary_String_Obj *p_so_or;
         for(int i=0; i<lexems.Total(); i++)
           {
            p_so=lexems.GetNodeAtIndex(i);
            if(p_so.Key()=="and")
              {
               p_so_and=p_so;
              }
            if(p_so.Key()=="or")
              {
               p_so_or=p_so;
              }
           }
         if(copyExpressions.Total()-index>0)
           {
            if((copyExpressions.At(index)==p_so_and.Value() && p_so_and.Key()=="and") || (copyExpressions.At(index)==p_so_or.Value() && p_so_or.Key()=="or"))
              {
               if(copyExpressions.Total()-index<2)
                 {
                  condExp=copyExpressions.At(index);
                  Print(StringFormat("Error at %s in condition part.",condExp.Text()));
                  //--- return 
                  return (NULL);
                 }
               //--- Set and/or for conditions list
               OperatorType newOp=NULL;
               if(copyExpressions.At(index)==p_so_and.Value() && p_so_and.Key()=="and")
                 {
                  newOp=And;
                 }
               if(copyExpressions.At(index)==p_so_or.Value() && p_so_or.Key()=="or")
                 {
                  newOp=Or;
                 }
               if(setOrAnd)
                 {
                  if(conds.Op()!=newOp)
                    {
                     Print("At the one nesting level cannot be mixed and/or operations.");
                     //--- return 
                     return (NULL);
                    }
                 }
               else
                 {
                  conds.Op(newOp);
                  setOrAnd=true;
                 }
               index=index+1;
              }
            else
              {
               string str;
               condExp=copyExpressions.At(index);
               str=condExp.Text();
               condExp=copyExpressions.At(index+1);
               Print(StringFormat("%s cannot goes after %s",str,condExp.Text()));
               //--- return 
               return (NULL);
              }
           }
        }
      //--- return conditions
      return (conds);
     }
  }
//+------------------------------------------------------------------+
//| Parse conclusion                                                 |
//+------------------------------------------------------------------+
static CSingleCondition *CRuleParser::ParseConclusion(CArrayObj *conditionExpression,CList *out,CList *lexems)
  {
   CArrayObj *copyExpression=conditionExpression;
//--- Remove extra brackets
   CDictionary_String_Obj *p_so;
   CDictionary_String_Obj *p_so_open;
   CDictionary_String_Obj *p_so_close;
   for(int i=0; i<lexems.Total(); i++)
     {
      p_so=lexems.GetNodeAtIndex(i);
      if(p_so.Key()=="(")
        {
         p_so_open=p_so;
        }
      if(p_so.Key()==")")
        {
         p_so_close=p_so;
        }
     }
   int index=0;
   int Total=copyExpression.Total();
   while(Total>=2 && (copyExpression.At(index)==p_so_open.Value() && copyExpression.At(index+conditionExpression.Total()-1)==p_so_close.Value()))
     {
      index=index+1;
      Total=Total-2;
     }
   if(Total!=3)
     {
      Print("Conclusion part of the rule should be in form: 'variable is term'");
      //--- return 
      return (NULL);
     }
//--- Parse variable
   CLexem *exprVariable=copyExpression.At(index);
   if(!exprVariable.IsTypeOf(TYPE_CLASS_VarLexem))
     {
      Print(StringFormat("Wrong identifier '%s' in conclusion part of the rule.",exprVariable.Text()));
      //--- return 
      return (NULL);
     }
   CVarLexem *varLexem=exprVariable;
   if(varLexem.Input()==true)
     {
      Print("The variable in conclusion part must be an output variable.");
      //--- return 
      return (NULL);
     }
//--- Parse 'is' lexem
   CLexem *exprIs=copyExpression.At(index+1);
   for(int i=0; i<lexems.Total(); i++)
     {
      p_so=lexems.GetNodeAtIndex(i);
      if(p_so.Key()=="is")
        {
         break;
        }
     }
   if(exprIs!=p_so.Value())
     {
      Print(StringFormat("'is' keyword must go after %s identifier.",varLexem.Text()));
      //--- return 
      return (NULL);
     }
//--- Parse term
   CLexem *exprTerm=copyExpression.At(index+2);
   if(!exprTerm.IsTypeOf(TYPE_CLASS_TermLexem))
     {
      Print(StringFormat("Wrong identifier '%s' in conclusion part of the rule.",exprTerm.Text()));
     }
   IAltLexem *altLexem=exprTerm;
   CTermLexem *termLexem=NULL;
   do
     {
      if(!altLexem.IsTypeOf(TYPE_CLASS_TermLexem))
        {
         continue;
        }
      termLexem=altLexem;
      if(varLexem.Var().Values().IndexOf(termLexem.Term())==-1)
        {
         termLexem=NULL;
         continue;
        }
     }
   while((altLexem=altLexem.Alternative())!=NULL && termLexem==NULL);
   if(termLexem==NULL)
     {
      Print(StringFormat("Wrong identifier '%s' in conclusion part of the rule.",exprTerm.Text()));
      //--- return 
      return (NULL);
     }
//--- Return fuzzy rule's conclusion
   INamedVariable *var=varLexem.Var();
   INamedValue *term=termLexem.Term();
//--- return single condition
   return new CSingleCondition(var, term, false);
  }
//+------------------------------------------------------------------+
