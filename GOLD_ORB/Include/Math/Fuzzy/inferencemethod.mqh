//+------------------------------------------------------------------+
//|                                              inferencemethod.mqh |
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
//+------------------------------------------------------------------+
//| Purpose: Contains a number of enumerations,                      |
//| for the convenience of working with other files                  |
//+------------------------------------------------------------------+
//+------------------------------------------------------------------+
//| And evaluating method                                            |
//+------------------------------------------------------------------+
enum AndMethod
  {
   MinAnd,                            // Minimum: min(a, b)
   ProductionAnd                      // Production: a * b
  };
//+------------------------------------------------------------------+
//| Or evaluating method                                             |
//+------------------------------------------------------------------+
enum OrMethod
  {
   MaxOr,                             // Maximum: max(a, b)
   ProbabilisticOr                    // Probabilistic OR: a + b - a * b
  };
//+------------------------------------------------------------------+
//| Fuzzy implication method                                         |
//+------------------------------------------------------------------+
enum ImplicationMethod
  {
   MinIpm,                            // Truncation of output fuzzy set
   ProductionImp                      // Scaling of output fuzzy set
  };
//+------------------------------------------------------------------+
//| Aggregation method for membership functions                      |
//+------------------------------------------------------------------+
enum AggregationMethod
  {
   MaxAgg,                            // Maximum of rule outpus
   SumAgg                             // Sum of rule output 
  };
//+------------------------------------------------------------------+
//| Defuzzification method                                           |
//+------------------------------------------------------------------+
enum DefuzzificationMethod
  {
   CentroidDef,                       // Center of area of fuzzy result MF
   BisectorDef,                       // The point divides the area under the MF into two equal
   AverageMaximumDef,                 // Arithmetic mean of all the maxima of the MF
   LargestMaximumDef,                 // The largest of the maxima of the membership function 
   SmallestMaximumDef                 // The smallest of the maxima of the membership function
  };
//+------------------------------------------------------------------+
//| Type of varriable and term                                       |
//+------------------------------------------------------------------+
enum EnType
  {
   TYPE_CLASS_INamedValue,            // Base class
   TYPE_CLASS_INamedVariable,         // INamedVariable : INamedValue
   TYPE_CLASS_NamedVariableImpl,      // NamedVariableImpl : INamedVariable
   TYPE_CLASS_NamedValueImpl,         // NamedValueImpl : INamedValue
   TYPE_CLASS_FuzzyTerm,              // FuzzyTerm : NamedValueImpl
   TYPE_CLASS_FuzzyVariable,          // FuzzyVariable : NamedVariableImpl
   TYPE_CLASS_SugenoVariable,         // SugenoVariable : NamedVariableImpl
   TYPE_CLASS_ISugenoFunction,        // ISugenoFunction : NamedValueImpl
   TYPE_CLASS_LinearSugenoFunction    // LinearSugenoFunction : ISugenoFunction
  };
//+------------------------------------------------------------------+
//| Type of expression                                               |
//+------------------------------------------------------------------+
enum EnLexem
  {
   TYPE_CLASS_IExpression,            // Base class
   TYPE_CLASS_Lexem,                  // Lexem : IExpression
   TYPE_CLASS_ConditionExpression,    // ConditionExpression : IExpression   
   TYPE_CLASS_VarLexem,               // VarLexem : Lexem
   TYPE_CLASS_KeywordLexem,           // KeywordLexem : Lexem
   TYPE_CLASS_AltLexem,               // AltLexem : Lexem
   TYPE_CLASS_TermLexem               // TermLexem : AltLexem 
  };
//+------------------------------------------------------------------+
//| Type of condition                                                |
//+------------------------------------------------------------------+
enum EnCondition
  {
   TYPE_CLASS_ICondition,             // Base class 
   TYPE_CLASS_Conditions,             // Conditions : ICondition
   TYPE_CLASS_SingleCondition,        // SingleCondition : ICondition
   TYPE_CLASS_FuzzyCondition          // FuzzyCondition : SingleCondition      
  };
//+------------------------------------------------------------------+
//| Type of rule                                                     |
//+------------------------------------------------------------------+ 
enum EnRule
  {
   TYPE_CLASS_IParsableRule,          // Base class
   TYPE_CLASS_GenericFuzzyRule,       // GenericFuzzyRule : IParsableRule
   TYPE_CLASS_MamdaniFuzzyRule,       // MamdaniFuzzyRule : GenericFuzzyRule
   TYPE_CLASS_SugenoFuzzyRule         // SugenoFuzzyRule : GenericFuzzyRule
  };
//+------------------------------------------------------------------+
