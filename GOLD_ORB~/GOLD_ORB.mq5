//+------------------------------------------------------------------+
//|                                               MQL5 Practice1.mq5 |
//|                                              Playground Inc 2021 |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "Playground Inc 2021"
#property link      "https://www.mql5.com"
#property version   "1.00"

/* 
Author: Ulysses O. Andulte
Date Created: 10/26/2022
*/


//Include Files
#include "Include\Trade.mqh"
#include "Include\TrailingStops.mqh"
#include "Include\price_action.mqh"
#include "Include\Indicators.mqh"
#include "Include\TradeVirtual.mqh"
#include "Include\TrailingStopsVirtual.mqh"
#include "Include\MoneyManagement.mqh"
#include "Include\Math\Stat\Normal.mqh"
#include "Include\RiskManagement.mqh"



//Input Variables
input group  "SymBolInformation"
input int StartOfTradingHour_ServerTime = 1;

input group  "Trade Management"
input int TakeProfit =1200;
input int StopLoss =400;
input int MaxTradePerDay =2;
input bool LongPosition = true;
input bool ShortPosition = true;

input group "Trail Management"
input bool EnableTrail=true;

input group "Risk Management"
input double MaxEquityDrawdownPercent = 10;
input double MaxRiskPerTradePercent = 1;
input double FixedVolume = 0.1;


input group "Advanced Equity Monitoring Module"
input bool SlopeDetection = false;
input int LossStreakCounter = 0;

input group "Indicators"
input int PriceActionORB_CandleComposition = 3;


//Global Variables
bool execute_trade;
double capital;
int indicator_2;
bool indicator_3;
double TradeVolume = FixedVolume;


//Class objects

//Trade management Module
//++++++++++++++++++++++++
CTrade trade; //a class for executing orders on the server
CTrailing trail; //a class for trail stop

//Indicators Module
//++++++++++++++++++++++++
Price_Action pa;// an indicator class for price action
CiMA MA100; //an indicator class for moving averages


//Virtual Trading Environment Module
//++++++++++++++++++++++++++++++++++++
CTradeVirtual tradevirtual;// a class for executing orders on virtual trade environment
CTrailingVirtual trailvirtual; //a class for trail stop virtual
VirtualTradeInfo VTrade; //a class for storing virtual information: details on  position, deals and closed trades




//Working Code
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
// Event handler: Initialization
int OnInit()
  {
  
   //Initialize Price action object for default or user input
   pa.Init();
   pa.candle_composition = PriceActionORB_CandleComposition;
   pa.trades_per_day = MaxTradePerDay;
   pa.StartOfTradinghour_servertime = StartOfTradingHour_ServerTime;
   

   //Initialize virtual trade environment
   tradevirtual.Init(VTrade);

   //******************************
   //Extra Variables (test cases)
   execute_trade = true;
   capital = AccountInfoDouble(ACCOUNT_BALANCE);
   return("Initialization Success");
  }


//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
// Event handler: Execute each tick that will arrive from the server
void OnTick()
  {
   
   
   RiskManagementModule();
   if(EnableTrail) TrailModule();
   if(pa.new_candle_check2())
     {
      IndicatorModule();
      ExecuteOrders();
     }
  }

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
//End of Program




/////////////////////////////////////////////////////////////////////
//                     Functions
/////////////////////////////////////////////////////////////////////

//+------------------------------------------------------------------+
//| //Risk Management Module     
//|      
//|  **MonitoringVirtualPosition - Monitors and Closes virtual position 
//|                                and update virtualtrade information
//|  
                              
//+------------------------------------------------------------------+
void RiskManagementModule(void)

  {
  
   //Virtual Equity Monitoring
   MonitorVirtualPostion(VTrade); 

   //Setting Up Equity Trail, will stop executing order once hit, off by default if input is "0"
   if(MaxEquityDrawdownPercent!=0)
     {
      if(AccountInfoDouble(ACCOUNT_BALANCE)>capital)
         capital = AccountInfoDouble(ACCOUNT_BALANCE);
      if(100*((AccountInfoDouble(ACCOUNT_BALANCE) - capital)/capital) < -1* MaxEquityDrawdownPercent)
         printf(100*((AccountInfoDouble(ACCOUNT_BALANCE) - capital)/capital));
         execute_trade = false;
     }


   //This module will detect if Lossing streak ended depending on the input integer and if equity is recovering, upward
   if(SlopeDetection || LossStreakCounter!=0)
     {
      bool LossStreak_flag = LossStreakCounter(VTrade,3); //Lossing Streak Detection
      bool Slope_Equity_Flag = CheckSlope(VTrade,12); // Slope Equity Monitoring

      if(LossStreak_flag)
         execute_trade = false;
      if(Slope_Equity_Flag == true && LossStreak_flag == false)
         execute_trade = true;
     }
   
   //Dynamic Position Sizing Relative to port size, or FixedVolume for default if no input in MaxRiskPerTradePercent
   TradeVolume = MoneyManagement(_Symbol,FixedVolume,MaxRiskPerTradePercent,StopLoss);


  }



//+------------------------------------------------------------------+
//| //Trail Stop Module                                              |
//+------------------------------------------------------------------+
void TrailModule(void)

  {
   int i = 0;
   int total_positions = PositionsTotal();


   //RealPort Trail Module, will loop on all open positions to check if trail is hit
   if(total_positions!=0)
      for(i=0 ; i <= total_positions -1; i++)
        {
         trail.TrailingStop(PositionGetTicket(i),700,100,10); //  set 700 trail stop below sa TP then min profit is  100 for secure profits, 10 is the step size
        }


   //Virtual Port Trail Module, will loop on all open positions to check if trail is hit
   int j=0;
   int total_positions_virtual = ArraySize(VTrade.position);
   if(total_positions_virtual!=0)
      for(j=0 ; j <= total_positions_virtual -1; j++)
        {
         trailvirtual.TrailingStop(VTrade,j,700,100,10); //  set 700 trail stop below sa TP then min profit is  100 for secure profits, 10 is the step size
        }
  }




//+------------------------------------------------------------------+
//| //Indicators Module                                              |
//+------------------------------------------------------------------+
void IndicatorModule(void)
  {


   //Price action indicator
   //outputs "11" for Long position signal and "10" for Short position signal
   indicator_2 = pa.Open_Range_Breakout(); 


   // Moving Average indicator
   MA100.Init(_Symbol,PERIOD_CURRENT,100,0,MODE_SMA,PRICE_CLOSE); 
   double ma = MA100.Main(0); // get the value of the latest ma value wrt to latest candle
   iClose(_Symbol,_Period,1) > ma ? indicator_3 = true:indicator_3 = false; //compare the value with the candle

  }



//+------------------------------------------------------------------+
//| //Trade Execution Module                                          |
//+------------------------------------------------------------------+

void ExecuteOrders(void)

  {

//Buy/Sell Order: //Execute buy/sell orders given the indicators and user inputs if its enabled
   if(indicator_2 == 11 && LongPosition)
     {

      tradevirtual.Buy(VTrade,_Symbol,TradeVolume,StopLoss,TakeProfit);
      if(execute_trade)//if this is false then the equity hits its maximum draw down as input by the user
         trade.Buy(_Symbol,TradeVolume,StopLoss,TakeProfit);
     }


   if(indicator_2 == 10 && ShortPosition)
     {
      tradevirtual.Sell(VTrade,_Symbol,TradeVolume,StopLoss,TakeProfit);
      if(execute_trade)//if this is false then the equity hits its maximum draw down as input by the user
         trade.Sell(_Symbol,TradeVolume,StopLoss,TakeProfit);

     }

  }


//+------------------------------------------------------------------+
