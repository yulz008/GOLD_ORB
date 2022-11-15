//+------------------------------------------------------------------+
//|                                                 TradeVirtual.mqh |
//|                                              Playground Inc 2021 |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "Playground Inc 2021"
#property link      "https://www.mql5.com"
#property version   "1.00"


void MonitorVirtualPostion(VirtualTradeInfo &VTrade)

  {

   int i = 0;
   int total_positions = ArraySize(VTrade.position);
//printf(ArraySize(VTrade.position));


   if(total_positions!=0)
     {

      for(i=0 ; i <= total_positions -1; i++)

        {


         if(VTrade.position[i].type == "long" && SymbolInfoDouble(_Symbol,SYMBOL_BID) >= VTrade.position[i].tp &&VTrade.position[i].tp > 0)

           {

            //update deals
            tradevirtual.TakeProfit(VTrade,i);

            //update close trade
            //Update position
            VTrade.position[i].comment = "closed";
           }

         //Print("Stoploss = ",  VTrade.position[i].sl);
         //Print("Current_Price = ", SymbolInfoDouble(_Symbol,SYMBOL_BID));
         if(VTrade.position[i].type == "long" && SymbolInfoDouble(_Symbol,SYMBOL_BID) <= VTrade.position[i].sl)

           {

            //update deals
            tradevirtual.StopLoss(VTrade,i);

            //update close trade
            //Update position
            VTrade.position[i].comment = "closed";

           }


         if(VTrade.position[i].type == "short" && SymbolInfoDouble(_Symbol,SYMBOL_ASK) <= VTrade.position[i].tp && VTrade.position[i].tp > 0)

           {
            //update deals
            tradevirtual.TakeProfit(VTrade,i);

            //update close trade
            //Update position
            VTrade.position[i].comment = "closed";

           }


         if(VTrade.position[i].type == "short" && SymbolInfoDouble(_Symbol,SYMBOL_ASK) >= VTrade.position[i].sl)

           {

            //update deals
            tradevirtual.StopLoss(VTrade,i);


            //update close trade
            //Update position
            VTrade.position[i].comment = "closed";

           }

        }


      //Update Position container - resizing
      bool position_update = false;
      i=0;
      do

        {
         total_positions = ArraySize(VTrade.position);
         for(i=0 ; i <= total_positions -1; i++)

           {

            if(total_positions ==0)
               position_update = true;

            if(VTrade.position[i].comment == "closed")
              {
               //remove this element on the array
               ArrayRemove(VTrade.position,i,1);


               if(i== total_positions-1)

                 {
                  position_update = true;
                 }

               break;

              }

            if(i== total_positions-1)

              {
               position_update = true;
              }

           }
        }
      while(position_update == false);

     }

  }
//+------------------------------------------------------------------+

//+------------------------------------------------------------------+

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
bool CheckSlope(VirtualTradeInfo &VTrade, int offset_length)
//+------------------------------------------------------------------+
  {


   int deal_size = ArraySize(VTrade.deals);
   int index = deal_size -1;
//int offset_length = 5;

   if(deal_size > offset_length)
     {

      if((VTrade.deals[index].balance - VTrade.deals[index-(offset_length)].balance)>0)

        {


         return(true);
        }
      else
         return(false);
     }

   else
      return(false);

  }
//+------------------------------------------------------------------+
bool LossStreakCounter(VirtualTradeInfo &VTrade, double streak)
  {
   
   
   int index = ArraySize(VTrade.closetrades)-1;
   double container[];
   ArrayResize(container,streak);
   int i =0;

   
   if(ArraySize(VTrade.closetrades) >= streak)
     {

      for(i = index ; i> index - streak; i--)

        {

         if(VTrade.closetrades[i].result == "WIN")
            container[index - i] = 0;

         if(VTrade.closetrades[i].result == "LOSS")
            container[index - i] = 1;

        }

      if(MathSum(container) == streak)
        {
         return(true);
        }
      else
         return(false);
     }

   else
      return (false);
      
     

  }
//+------------------------------------------------------------------+

