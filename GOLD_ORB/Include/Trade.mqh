//+------------------------------------------------------------------+
//|                                                        Trade.mqh |
//|                                              Playground Inc 2021 |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "Playground Inc 2021"
#property link      "https://www.mql5.com"
#property version   "1.00"



#define MAX_RETRIES 5  // Max retries on error
#define RETRY_DELAY 3000 // Retry delay in ms

#include "errordescription.mqh"



//+------------------------------------------------------------------+
//| CTrade Class - Open, Close and Modify Orders                                                           |
//+------------------------------------------------------------------+


class CTrade
  {
protected:

   //create an object "request" inheriting the class MqlTradeRequest
   MqlTradeRequest   request;

   bool              OpenPosition(string pSymbol, ENUM_ORDER_TYPE pType, double pVolume, double pStop = 0, double pProfit = 0, string pComment = NULL);
   bool              OpenPending(string pSymbol, ENUM_ORDER_TYPE pType, double pVolume, double pPrice, double pStop = 0, double pProfit = 0, double pStopLimit = 0, datetime pExpiration = 0, string pComment = NULL);
   void              LogTradeRequest();

   ulong             magicNumber;
   ulong             deviation;
   ENUM_ORDER_TYPE_FILLING fillType;

public:
   MqlTradeResult    result;

   bool              Buy(string pSymbol, double pVolume, double pStop = 0, double pProfit = 0, string pComment = NULL);
   bool              Sell(string pSymbol, double pVolume, double pStop = 0, double pProfit = 0, string pComment = NULL);
   void              MagicNumber(ulong pMagic);
   void              Deviation(ulong pDeviation);
   void              FillType(ENUM_ORDER_TYPE_FILLING pFill);

  };





// Open position
// 110 - 120 Code explanation on "Expert Advisor Programming for MetaTrader5""
bool CTrade::OpenPosition(string pSymbol, ENUM_ORDER_TYPE pType, double pVolume, double pStop = 0, double pProfit = 0, string pComment = NULL)
  {
   ZeroMemory(request);
   ZeroMemory(result);

   request.action = TRADE_ACTION_DEAL;
   request.symbol = pSymbol;
   request.type = pType;
   request.comment = pComment;
   request.deviation = 50;
   request.type_filling = ORDER_FILLING_IOC;
   request.volume = pVolume;

 


// Order loop
   int retryCount = 0;
   int checkCode = 0;

   do
     {
      //assign market price to local variable request.price
      if(pType == ORDER_TYPE_BUY)
        {
         request.price = SymbolInfoDouble(pSymbol,SYMBOL_ASK);
         if(pStop > 0)
            request.sl = request.price - (pStop * _Point);
         if(pProfit > 0)
            request.tp = request.price + (pProfit * _Point);


        }
      else
         if(pType == ORDER_TYPE_SELL)

           {
            request.price = SymbolInfoDouble(pSymbol,SYMBOL_BID);
            if(pStop > 0)
               request.sl = request.price + (pStop * _Point);
            if(pProfit > 0)
               request.tp = request.price - (pProfit * _Point);


           }
      bool sent = OrderSend(request,result);//send order

      checkCode = CheckReturnCode(result.retcode); //stores retcode on local variabe checkcode

      if(checkCode == CHECK_RETCODE_OK)
         break; // exits the loop after order is confirmed on checkcode.
      else
         if(checkCode == CHECK_RETCODE_ERROR) //error handling for failed order
           {
            string errDesc =  (result.retcode);
            Alert("Open market order: Error ",result.retcode," - ",errDesc);
            LogTradeRequest();
            break;
           }
         else
           {
            Print("Server error detected, retrying...");
            Sleep(RETRY_DELAY);
            retryCount++;
           }
     }
   while(retryCount < MAX_RETRIES);

   if(retryCount >= MAX_RETRIES)
     {
      string errDesc = TradeServerReturnCodeDescription(result.retcode);
      Alert("Max retries exceeded: Error ",result.retcode," - ",errDesc);
     }

   string orderType = CheckOrderType(pType);

   string errDesc = TradeServerReturnCodeDescription(result.retcode);
   Print("Open ",orderType," order #",result.order,": ",result.retcode," - ",errDesc,", Volume: ",result.volume,", Price: ",result.price,", Bid: ",result.bid,", Ask: ",result.ask);



//Returns true if order is successfull and return false if its not
   if(checkCode == CHECK_RETCODE_OK)
     {
      Comment(orderType," position opened at ",result.price," on ",pSymbol);
      return(true);
     }
   else
      return(false);
  }



// Open pending order
bool CTrade::OpenPending(string pSymbol,ENUM_ORDER_TYPE pType,double pVolume,double pPrice,double pStop=0.000000,double pProfit=0.000000,double pStopLimit = 0,datetime pExpiration=0,string pComment=NULL)
  {
   ZeroMemory(request);
   ZeroMemory(result);

   request.action = TRADE_ACTION_PENDING;
   request.symbol = pSymbol;
   request.type = pType;
   request.sl = pStop;
   request.tp = pProfit;
   request.comment = pComment;
   request.price = pPrice;
   request.volume = pVolume;
   request.stoplimit = pStopLimit;
   request.deviation = deviation;
   request.type_filling = ORDER_FILLING_IOC;
   request.magic = magicNumber;

   if(pExpiration > 0)
     {
      request.expiration = pExpiration;
      request.type_time = ORDER_TIME_SPECIFIED;
     }
   else
      request.type_time = ORDER_TIME_GTC;

// Order loop
   int retryCount = 0;
   int checkCode = 0;

   do
     {
      bool sent = OrderSend(request,result);

      checkCode = CheckReturnCode(result.retcode);

      if(checkCode == CHECK_RETCODE_OK)
         break;
      else
         if(checkCode == CHECK_RETCODE_ERROR)
           {
            string errDesc = TradeServerReturnCodeDescription(result.retcode);
            Alert("Open pending order: Error ",result.retcode," - ",errDesc);
            LogTradeRequest();
            break;
           }
         else
           {
            Print("Server error detected, retrying...");
            Sleep(RETRY_DELAY);
            retryCount++;
           }
     }
   while(retryCount < MAX_RETRIES);

   if(retryCount >= MAX_RETRIES)
     {
      string errDesc = TradeServerReturnCodeDescription(result.retcode);
      Alert("Max retries exceeded: Error ",result.retcode," - ",errDesc);
     }

   string orderType = CheckOrderType(pType);
   string errDesc = TradeServerReturnCodeDescription(result.retcode);

   Print("Open ",orderType," order #",result.order,": ",result.retcode," - ",errDesc,", Volume: ",result.volume,", Price: ",request.price,
         ", Bid: ",SymbolInfoDouble(pSymbol,SYMBOL_BID),", Ask: ",SymbolInfoDouble(pSymbol,SYMBOL_ASK),", SL: ",request.sl,", TP: ",request.tp,
         ", Stop Limit: ",request.stoplimit,", Expiration: ",request.expiration);

   if(checkCode == CHECK_RETCODE_OK)
     {
      Comment(orderType," order opened at ",request.price," on ",pSymbol);
      return(true);
     }
   else
      return(false);
  }




//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
void CTrade::LogTradeRequest()
  {
   Print("MqlTradeRequest - action:",request.action,", comment:",request.comment,", deviation:",request.deviation,", expiration:",request.expiration,", magic:",request.magic,", order:",request.order,", position:",request.position,", position_by:",request.position_by,", price:",request.price,", ls:",request.sl,", stoplimit:",request.stoplimit,", symbol:",request.symbol,", tp:",request.tp,", type:",request.type,", type_z:",request.type_filling,", type_time:",request.type_time,", volume:",request.volume);
   Print("MqlTradeResult - ask:",result.ask,", bid:",result.bid,", comment:",result.comment,", deal:",result.deal,", order:",result.order,", price:",result.price,", request_id:",result.request_id,", retcode:",result.retcode,", retcode_external:",result.retcode_external,", volume:",result.volume);
  }


// Trade opening shortcuts
bool CTrade::Buy(string pSymbol,double pVolume,double pStop=0.000000,double pProfit=0.000000,string pComment=NULL)
  {
   bool success = OpenPosition(pSymbol,ORDER_TYPE_BUY,pVolume,pStop,pProfit,pComment);
   return(success);
  }

//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
bool CTrade::Sell(string pSymbol,double pVolume,double pStop=0.000000,double pProfit=0.000000,string pComment=NULL)
  {
   bool success = OpenPosition(pSymbol,ORDER_TYPE_SELL,pVolume,pStop,pProfit,pComment);
   return(success);
  }




// Set magic number
void CTrade::MagicNumber(ulong pMagic)
  {
   magicNumber = pMagic;
  }


// Set deviation
void CTrade::Deviation(ulong pDeviation)
  {
   deviation = pDeviation;
  }


// Set fill type
void CTrade::FillType(ENUM_ORDER_TYPE_FILLING pFill)
  {
   fillType = pFill;
  }


// Return code check
int CheckReturnCode(uint pRetCode)
  {
   int status;
   switch(pRetCode)
     {
      case TRADE_RETCODE_REQUOTE:
      case TRADE_RETCODE_CONNECTION:
      case TRADE_RETCODE_PRICE_CHANGED:
      case TRADE_RETCODE_TIMEOUT:
      case TRADE_RETCODE_PRICE_OFF:
      case TRADE_RETCODE_REJECT:
      case TRADE_RETCODE_ERROR:

         status = CHECK_RETCODE_RETRY;
         break;

      case TRADE_RETCODE_DONE:
      case TRADE_RETCODE_DONE_PARTIAL:
      case TRADE_RETCODE_PLACED:
      case TRADE_RETCODE_NO_CHANGES:

         status = CHECK_RETCODE_OK;
         break;

      default:
         status = CHECK_RETCODE_ERROR;
     }

   return(status);
  }




//+------------------------------------------------------------------+
//| Miscellaneous Functions & Enumerations                                                            |
//+------------------------------------------------------------------+


enum ENUM_CHECK_RETCODE
  {
   CHECK_RETCODE_OK,
   CHECK_RETCODE_ERROR,
   CHECK_RETCODE_RETRY
  };


//+------------------------------------------------------------------+
//|                                                                  |
//+------------------------------------------------------------------+
string CheckOrderType(ENUM_ORDER_TYPE pType)
  {
   string orderType;
   if(pType == ORDER_TYPE_BUY)
      orderType = "buy";
   else
      if(pType == ORDER_TYPE_SELL)
         orderType = "sell";
      else
         if(pType == ORDER_TYPE_BUY_STOP)
            orderType = "buy stop";
         else
            if(pType == ORDER_TYPE_BUY_LIMIT)
               orderType = "buy limit";
            else
               if(pType == ORDER_TYPE_SELL_STOP)
                  orderType = "sell stop";
               else
                  if(pType == ORDER_TYPE_SELL_LIMIT)
                     orderType = "sell limit";
                  else
                     if(pType == ORDER_TYPE_BUY_STOP_LIMIT)
                        orderType = "buy stop limit";
                     else
                        if(pType == ORDER_TYPE_SELL_STOP_LIMIT)
                           orderType = "sell stop limit";
                        else
                           orderType = "invalid order type";
   return(orderType);
  }
//+------------------------------------------------------------------+




