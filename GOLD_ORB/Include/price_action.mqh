//+------------------------------------------------------------------+
//|                                                 price_action.mqh |
//|                                              Playground Inc 2021 |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "Playground Inc 2021"
#property link      "https://www.mql5.com"
#property version   "1.00"




#include "TradeVirtual.mqh"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
class Price_Action
  {
protected:


   //
   int                      price_action_return_flag;


   //
   bool                     short_position_flag;
   bool                     long_position_flag;
   bool                     trade_today;

   //
   bool                     candle_reset_resistance;
   bool                     candle_reset_support;
   int                      candle_counter_support;
   int                      candle_counter_resistance;


   //
   int                      m_period_flags;
   MqlDateTime              m_last_tick_time;
   datetime                 initial_time_new_trading_day;
   MqlDateTime              time;
   MqlTick                  tick;
   int                      init_time;


   //
   bool                     new_trading_day_flag;
   bool                     first_candle;

   //
   CandleInfo               previous_candle;
   CandleInfo               range;




   //protected function
   void                     GetCandleInfo(CandleInfo &candle);
   void                     ResetVariables();
   void                     FirstCandleUpdateSnR(CandleInfo &range, CandleInfo &previous_candle);
   void                     CandleUpdateSnR(CandleInfo &range, CandleInfo &previous_candle);
   int                      GetBuySellSignal(CandleInfo &range, CandleInfo &previous_candle);
   void                     UpdateFlags(void);
   void                     UpdateDrawingObjects(void);
   int                      TimeframesFlags(MqlDateTime &time);
   void                     TimeframeAdd(ENUM_TIMEFRAMES period);


public:

   //variables
   int                      candle_composition;
   int                      trades_per_day;
   int                      StartOfTradinghour_servertime;

   //public functions
   bool                     new_candle_check(void);
   bool                     new_candle_check2(void);
   int                      Open_Range_Breakout();
   void                     Init();


  };


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void Price_Action::GetCandleInfo(CandleInfo &candle)
  {


//assign OPEN,CLOSE, HIGH, LOW

//previous candle is bullish
   if(iClose(_Symbol,_Period,1)> iOpen(_Symbol,_Period,1))

     {

      candle.body_high =iClose(_Symbol,_Period,1);
      candle.body_low =iOpen(_Symbol,_Period,1);
      candle.wick_high =iHigh(_Symbol,_Period,1);
      candle.wick_low= iLow(_Symbol,_Period,1);
      candle.direction = true;

     }

//previous candle is bearish
   else

     {

      candle.body_high =iOpen(_Symbol,_Period,1);
      candle.body_low =iClose(_Symbol,_Period,1);
      candle.wick_high =iHigh(_Symbol,_Period,1);
      candle.wick_low= iLow(_Symbol,_Period,1);
      candle.direction = false;
     }

  }

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void Price_Action::ResetVariables()
  {

   range.wick_high = 0;
   range.wick_low = 0;
   candle_counter_support = 0;
   candle_counter_resistance = 0;
   trade_today = false;
   candle_reset_resistance = false;
   candle_reset_support = false;
   previous_candle.body_high =0;
   previous_candle.body_low =0;
   previous_candle.wick_high =0;
   previous_candle.wick_low=0;
   previous_candle.direction = NULL;
   short_position_flag = false;
   long_position_flag = false;
   trade_today = true;
   initial_time_new_trading_day = TimeCurrent()- 2*PeriodSeconds(_Period);
   
  }

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void Price_Action::FirstCandleUpdateSnR(CandleInfo &range, CandleInfo &previous_candle)

  {


   if(fabs((previous_candle.body_high/_Point) - (previous_candle.wick_high/_Point)) > 500)

     {

      range.wick_high = previous_candle.body_high;
      range.body_high = previous_candle.body_high;


     }

//Not long wick
   else

     {
      range.wick_high = previous_candle.wick_high;
      range.body_high = previous_candle.body_high;


     }

//Long wick scenario on support update
   if(fabs((previous_candle.body_low/_Point) - (previous_candle.wick_low/_Point)) > 500)

     {
      range.wick_low =  previous_candle.body_low;
      range.body_low = previous_candle.body_low;


     }

//Not long wick
   else
     {
      range.wick_low =  previous_candle.wick_low;
      range.body_low = previous_candle.body_low;



     }

  }



//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void Price_Action::CandleUpdateSnR(CandleInfo &range, CandleInfo &previous_candle)
  {
//update range.wick_high and range.wick_low on the incoming candles
   if(range.wick_high > 0 && range.wick_low > 0 && (candle_counter_resistance <= candle_composition) && (candle_counter_resistance > 0))

     {

      //this will update new resistance given that new candle establishes new high AND significant change (>1) difference AND bullish candle
      if(previous_candle.wick_high > range.wick_high && (fabs(previous_candle.wick_high-range.wick_high)>0.1)  && previous_candle.body_high > range.body_high && fabs(previous_candle.body_high - range.body_high)> 0.1)


        {

         //Long wick scenario on resistance update
         if(fabs((previous_candle.body_high/_Point) - (previous_candle.wick_high/_Point)) > 500)

           {
            range.wick_high = previous_candle.body_high;
            range.body_high = previous_candle.body_high;

            candle_counter_resistance = 1;
            candle_reset_resistance = true;

           }

         else

           {
            range.wick_high = previous_candle.wick_high;
            range.body_high = previous_candle.body_high;

            candle_counter_resistance = 1;
            candle_reset_resistance = true;

           }


        }

     }


   if(range.wick_high > 0 && range.wick_low > 0 && (candle_counter_support <= candle_composition) && (candle_counter_support > 0))

     {


      //this will update new support given that new candle establishes new low AND significant change (>1) difference AND bearish candle

      if(previous_candle.wick_low < range.wick_low && (fabs(previous_candle.wick_low-range.wick_low)>0.1) &&  fabs(previous_candle.body_low - range.body_low)> 0.1 && previous_candle.body_low < range.body_low)

        {


         //Long wick scenario on support update

         if(fabs((previous_candle.body_low/_Point) - (previous_candle.wick_low/_Point)) > 500)

           {
            range.wick_low =  previous_candle.body_low;
            range.body_low = previous_candle.body_low;


            candle_counter_support = 1;
            candle_reset_support = true;
            //printf("candle_reset_support");


           }


         else
           {
            range.wick_low =  previous_candle.wick_low;
            range.body_low = previous_candle.body_low;


            candle_counter_support = 1;
            candle_reset_support = true;
            //printf("candle_reset_support");


           }

        }

     }


  }


//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
int Price_Action::GetBuySellSignal(CandleInfo &range, CandleInfo &previous_candle)
  {

//execute buy orders

//buy order condition: Resistance_candle_counter greater than candle_composition (box range) and previous candle is bullish
   if((candle_counter_resistance> candle_composition) && previous_candle.direction == true)

     {
      //2nd condition a breakout from the range resistance (wick_high)
      if(previous_candle.body_high > range.wick_high)

        {
         //3rd condition is to ensure that EA will only open position once/twice per day depending on user setting..
         if(trade_today == true && long_position_flag == false)
           {

            //Update flags
            if(trades_per_day==2)
               long_position_flag = true;
            if(trades_per_day==1)
               trade_today = false;
            return(11); //return BuyLong Signal

           }
        }


     }



//execute sell orders
//sell order condition: Support_candle_counter greater than candle_composition (box range) and previous candle is bearish
   if((candle_counter_support > candle_composition)&& previous_candle.direction == false)

     {
      //2nd condition a breakdown from the range low (wick_lowh)
      if(previous_candle.body_low < range.wick_low)

        {
         //3rd condition is to ensure that EA will only open position once/twice per day depending on user setting..
         if(trade_today == true && short_position_flag == false)
           {

            //Update flags
            if(trades_per_day==2)
               short_position_flag = true;
            if(trades_per_day==1)
               trade_today = false;
            return(10); //return SellShort Signal

           }

        }
     }


   return (0); //return "0" if no signal was generated..

  }



//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void Price_Action::UpdateFlags(void)

  {

     {


      if(long_position_flag == true && short_position_flag == true)

        {
         trade_today = false;
        }


      //resets candle counter resistance/support or increments it
      if(candle_reset_resistance == true)
        {
         candle_counter_resistance= 1;
         candle_reset_resistance = false;
        }

      else
        {
         candle_counter_resistance=candle_counter_resistance+1;
        }


      if(candle_reset_support == true)
        {
         candle_counter_support= 1;
         candle_reset_support = false;
        }

      else
        {
         candle_counter_support=candle_counter_support+1;
        }

     }

  }




//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void Price_Action::UpdateDrawingObjects(void)
  {

   if(trade_today)
     {
     
     
     
      // Line Support and Resistance:
      //Box Object Create Resistance
           // ObjectCreate(0,"My Line",OBJ_HLINE,0,0,range.wick_high);
      //ObjectCreate(0,"My Line3",OBJ_HLINE,0,0,range.body_high);
      //ObjectSetInteger(0,"My Line3",OBJPROP_COLOR,clrBlue);

      //Box Object Create Support
          //ObjectCreate(0,"My Line1",OBJ_HLINE,0,0,range.wick_low);
      //ObjectCreate(0,"My Line2",OBJ_HLINE,0,0,range.body_low);
      //ObjectSetInteger(0,"My Line2",OBJPROP_COLOR,clrBlue);
      
      
      // Box Support and Resistance
      initial_time_new_trading_day = initial_time_new_trading_day;  // offset left border - time
      datetime current_time = TimeCurrent() + 2*PeriodSeconds(_Period); // offset right border - time
      
      ObjectCreate(_Symbol,"Rectangle",OBJ_RECTANGLE,0,initial_time_new_trading_day, range.wick_high,current_time,range.wick_low);

     }
  }



//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
int Price_Action::TimeframesFlags(MqlDateTime &time)

  {
//--- set flags for all timeframes
   int result=OBJ_ALL_PERIODS;
//--- if first check, then setting flags all timeframes
   if(m_last_tick_time.min==-1)
      return(result);
//--- check change time
   if(time.min==m_last_tick_time.min &&
      time.hour==m_last_tick_time.hour &&
      time.day==m_last_tick_time.day &&
      time.mon==m_last_tick_time.mon)
      return(OBJ_NO_PERIODS);
//--- new month?
   if(time.mon!=m_last_tick_time.mon)
      return(result);
//--- reset the "new month" flag
   result^=OBJ_PERIOD_MN1;
//--- new day?
   if(time.day!=m_last_tick_time.day)
      return(result);
//--- reset the "new day" and "new week" flags
   result^=OBJ_PERIOD_D1+OBJ_PERIOD_W1;
//--- temporary variables to speed up working with structures
   int curr,delta;
//--- new hour?
   curr=time.hour;
   delta=curr-m_last_tick_time.hour;
   if(delta!=0)
     {
      if(curr%2>=delta)
         result^=OBJ_PERIOD_H2;
      if(curr%3>=delta)
         result^=OBJ_PERIOD_H3;
      if(curr%4>=delta)
         result^=OBJ_PERIOD_H4;
      if(curr%6>=delta)
         result^=OBJ_PERIOD_H6;
      if(curr%8>=delta)
         result^=OBJ_PERIOD_H8;
      if(curr%12>=delta)
         result^=OBJ_PERIOD_H12;
      return(result);
     }
//--- reset all flags for hour timeframes
   result^=OBJ_PERIOD_H1+OBJ_PERIOD_H2+OBJ_PERIOD_H3+OBJ_PERIOD_H4+OBJ_PERIOD_H6+OBJ_PERIOD_H8+OBJ_PERIOD_H12;
//--- new minute?
   curr=time.min;
   delta=curr-m_last_tick_time.min;
   if(delta!=0)
     {
      if(curr%2>=delta)
         result^=OBJ_PERIOD_M2;
      if(curr%3>=delta)
         result^=OBJ_PERIOD_M3;
      if(curr%4>=delta)
         result^=OBJ_PERIOD_M4;
      if(curr%5>=delta)
         result^=OBJ_PERIOD_M5;
      if(curr%6>=delta)
         result^=OBJ_PERIOD_M6;
      if(curr%10>=delta)
         result^=OBJ_PERIOD_M10;
      if(curr%12>=delta)
         result^=OBJ_PERIOD_M12;
      if(curr%15>=delta)
         result^=OBJ_PERIOD_M15;
      if(curr%20>=delta)
         result^=OBJ_PERIOD_M20;
      if(curr%30>=delta)
         result^=OBJ_PERIOD_M30;
     }
//--- result
   return(result);

  }


//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void Price_Action::TimeframeAdd(ENUM_TIMEFRAMES period)

  {
   switch(period)
     {
      case PERIOD_M1:
         m_period_flags|=OBJ_PERIOD_M1;
         break;
      case PERIOD_M2:
         m_period_flags|=OBJ_PERIOD_M2;
         break;
      case PERIOD_M3:
         m_period_flags|=OBJ_PERIOD_M3;
         break;
      case PERIOD_M4:
         m_period_flags|=OBJ_PERIOD_M4;
         break;
      case PERIOD_M5:
         m_period_flags|=OBJ_PERIOD_M5;
         break;
      case PERIOD_M6:
         m_period_flags|=OBJ_PERIOD_M6;
         break;
      case PERIOD_M10:
         m_period_flags|=OBJ_PERIOD_M10;
         break;
      case PERIOD_M12:
         m_period_flags|=OBJ_PERIOD_M12;
         break;
      case PERIOD_M15:
         m_period_flags|=OBJ_PERIOD_M15;
         break;
      case PERIOD_M20:
         m_period_flags|=OBJ_PERIOD_M20;
         break;
      case PERIOD_M30:
         m_period_flags|=OBJ_PERIOD_M30;
         break;
      case PERIOD_H1:
         m_period_flags|=OBJ_PERIOD_H1;
         break;
      case PERIOD_H2:
         m_period_flags|=OBJ_PERIOD_H2;
         break;
      case PERIOD_H3:
         m_period_flags|=OBJ_PERIOD_H3;
         break;
      case PERIOD_H4:
         m_period_flags|=OBJ_PERIOD_H4;
         break;
      case PERIOD_H6:
         m_period_flags|=OBJ_PERIOD_H6;
         break;
      case PERIOD_H8:
         m_period_flags|=OBJ_PERIOD_H8;
         break;
      case PERIOD_H12:
         m_period_flags|=OBJ_PERIOD_H12;
         break;
      case PERIOD_D1:
         m_period_flags|=OBJ_PERIOD_D1;
         break;
      case PERIOD_W1:
         m_period_flags|=OBJ_PERIOD_W1;
         break;
      case PERIOD_MN1:
         m_period_flags|=OBJ_PERIOD_MN1;
         break;
      default:
         m_period_flags=WRONG_VALUE;
         break;
     }


  }

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
bool  Price_Action::new_candle_check(void)
  {


   TimeCurrent(time);


//checks if the new tick came from a new candle

   if(time.hour != init_time)
     {

      init_time = time.hour;
      return(true);//new candle detected

     }

   else
     {

      return(false);
     }

  }


//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+

bool  Price_Action::new_candle_check2(void)

  {
//MqlDateTime time;
   TimeCurrent(time);
   if(m_period_flags!=WRONG_VALUE && m_period_flags!=0)
      if((m_period_flags&TimeframesFlags(time))==0)
        {
         return(false);
        }
   m_last_tick_time=time;
//--- refresh indicators
   return(true);

  }




//+------------------------------------------------------------------+
int  Price_Action::Open_Range_Breakout()
//+------------------------------------------------------------------+
  {

   int y =0;

//get previous candle information then store it to a private member structure "previous_candle"
   GetCandleInfo(previous_candle);


   if(time.hour == StartOfTradinghour_servertime + 1) //flag to be used for timing of 1:00:00 start of trading day
      first_candle = true;


//1st candle of the day, this will only execute at candle open during the  first hour of the trading day at 1:00:00 or 6am at PH time
   if(time.hour == StartOfTradinghour_servertime && first_candle == true)
     {

      ResetVariables(); //reset all flags/variables: counters and range high and range low values at new trading day 6am
      
      first_candle = false;//this is a flag so conditions above will not be executed
      new_trading_day_flag = true;//flag for "2nd candle" and the succedding candles.
     }
//2nd candle of the day
   else
      if(new_trading_day_flag == true)
        {
         FirstCandleUpdateSnR(range,previous_candle); //computes 1st candle range at start of trading day
         new_trading_day_flag = false;
         UpdateFlags();
         UpdateDrawingObjects();
        }

      //Sucedding Candles
      else
        {
         CandleUpdateSnR(range,previous_candle); //Updates SNR ranges for the upcoming candles after the 2nd candle
         UpdateFlags(); //update affected flags
         UpdateDrawingObjects();
        }

   y = GetBuySellSignal(range,previous_candle); //Generate buy/sell signals returns "11" for buy long and  "10" for sell short

   return (y);

  }



//+------------------------------------------------------------------+
void  Price_Action::Init()
  {

   TimeframeAdd(Period());
   first_candle = true;
   candle_composition = 3;
   trades_per_day =2;
  }
//+------------------------------------------------------------------+
