//+------------------------------------------------------------------+
//|                                                TrailingStops.mqh |
//|                                              Playground Inc 2021 |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "Playground Inc 2021"
#property link      "https://www.mql5.com"
#property version   "1.00"



#include "errordescription.mqh"
#include "Trade.mqh"


//+------------------------------------------------------------------+
//| Trailing Stop Class                                              |
//+------------------------------------------------------------------+


class CTrailing
{
	protected:
		MqlTradeRequest request;
		
	public:
		MqlTradeResult result;
		
		bool TrailingStop(string pSymbol, int pTrailPoints, int pMinProfit = 0, int pStep = 10);
		bool TrailingStop(string pSymbol, double pTrailPrice, int pMinProfit = 0, int pStep = 10);
		
		bool TrailingStop(ulong pTicket, int pTrailPoints, int pMinProfit = 0, int pStep = 10);
		bool TrailingStop(ulong pTicket, double pTrailPrice, int pMinProfit = 0, int pStep = 10);
		
		bool BreakEven(string pSymbol, int pBreakEven, int pLockProfit=0);
		bool BreakEven(ulong pTicket, int pBreakEven, int pLockProfit=0);
};


// Trailing stop (points)
bool CTrailing::TrailingStop(string pSymbol,int pTrailPoints,int pMinProfit=0,int pStep=10)
{
	if(PositionSelect(pSymbol) == true && pTrailPoints > 0)
	{
		request.action = TRADE_ACTION_SLTP;
		request.symbol = pSymbol;
		
		long posType = PositionGetInteger(POSITION_TYPE);
		double currentStop = PositionGetDouble(POSITION_SL);
		double openPrice = PositionGetDouble(POSITION_PRICE_OPEN);
		
		double point = SymbolInfoDouble(pSymbol,SYMBOL_POINT);
		int digits = (int)SymbolInfoInteger(pSymbol,SYMBOL_DIGITS);
		
		if(pStep < 10) pStep = 10;
		double step = pStep * point;
		
		double minProfit = pMinProfit * point;
		double trailStop = pTrailPoints * point;
		currentStop = NormalizeDouble(currentStop,digits);
		
		double trailStopPrice;
		double currentProfit;
		
		// Order loop
		int retryCount = 0;
		int checkRes = 0;
		
		do 
		{
			if(posType == POSITION_TYPE_BUY)
			{
				trailStopPrice = SymbolInfoDouble(pSymbol,SYMBOL_BID) - trailStop;
				trailStopPrice = NormalizeDouble(trailStopPrice,digits);
				currentProfit = SymbolInfoDouble(pSymbol,SYMBOL_BID) - openPrice;
				
				if(trailStopPrice > currentStop + step && currentProfit >= minProfit)
				{
					request.sl = trailStopPrice;
					bool sent = OrderSend(request,result);
				}
				else return(false);
			}
			else if(posType == POSITION_TYPE_SELL)
			{
				trailStopPrice = SymbolInfoDouble(pSymbol,SYMBOL_ASK) + trailStop;
				trailStopPrice = NormalizeDouble(trailStopPrice,digits);
				currentProfit = openPrice - SymbolInfoDouble(pSymbol,SYMBOL_ASK);
				
				if((trailStopPrice < currentStop - step || currentStop == 0) && currentProfit >= minProfit)
				{	
					request.sl = trailStopPrice;
					bool sent = OrderSend(request,result);
				}
				else return(false);
			}
			
			checkRes = CheckReturnCode(result.retcode);
		
			if(checkRes == CHECK_RETCODE_OK) break;
			else if(checkRes == CHECK_RETCODE_ERROR)
			{
				string errDesc = TradeServerReturnCodeDescription(result.retcode);
				Alert("Trailing stop: Error ",result.retcode," - ",errDesc);
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
		
		string errDesc = TradeServerReturnCodeDescription(result.retcode);
		Print("Trailing stop: ",result.retcode," - ",errDesc,", Old SL: ",currentStop,", New SL: ",request.sl,", Bid: ",SymbolInfoDouble(pSymbol,SYMBOL_BID),", Ask: ",SymbolInfoDouble(pSymbol,SYMBOL_ASK),", Stop Level: ",SymbolInfoInteger(pSymbol,SYMBOL_TRADE_STOPS_LEVEL));
		
		if(checkRes == CHECK_RETCODE_OK) return(true);
		else return(false);
	}
	
	else return(false);
}


// Trailing stop (price)
bool CTrailing::TrailingStop(string pSymbol,double pTrailPrice,int pMinProfit=0,int pStep=10)
{
	if(PositionSelect(pSymbol) == true && pTrailPrice > 0)
	{
		request.action = TRADE_ACTION_SLTP;
		request.symbol = pSymbol;
		
		long posType = PositionGetInteger(POSITION_TYPE);
		double currentStop = PositionGetDouble(POSITION_SL);
		double openPrice = PositionGetDouble(POSITION_PRICE_OPEN);
		
		double point = SymbolInfoDouble(pSymbol,SYMBOL_POINT);
		int digits = (int)SymbolInfoInteger(pSymbol,SYMBOL_DIGITS);
		
		if(pStep < 10) pStep = 10;
		double step = pStep * point;
		double minProfit = pMinProfit * point;
				
		currentStop = NormalizeDouble(currentStop,digits);
		pTrailPrice = NormalizeDouble(pTrailPrice,digits);
		
		double currentProfit;
		
		int retryCount = 0;
		int checkRes = 0;
		
		double bid = 0, ask = 0;
		
		do 
		{
			if(posType == POSITION_TYPE_BUY)
			{
				bid = SymbolInfoDouble(pSymbol,SYMBOL_BID);
				currentProfit = bid - openPrice;
				if(pTrailPrice > currentStop + step && currentProfit >= minProfit) 
				{
					request.sl = pTrailPrice;
					bool sent = OrderSend(request,result);
				}
				else return(false);
			}
			else if(posType == POSITION_TYPE_SELL)
			{
				ask = SymbolInfoDouble(pSymbol,SYMBOL_ASK);
				currentProfit = openPrice - ask;
				if((pTrailPrice < currentStop - step || currentStop == 0) && currentProfit >= minProfit)
				{
					request.sl = pTrailPrice;
					bool sent = OrderSend(request,result);
				}
				else return(false);
			}
			
			checkRes = CheckReturnCode(result.retcode);
		
			if(checkRes == CHECK_RETCODE_OK) break;
			else if(checkRes == CHECK_RETCODE_ERROR)
			{
				string errDesc = TradeServerReturnCodeDescription(result.retcode);
				Alert("Trailing stop: Error ",result.retcode," - ",errDesc);
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
		
		string errDesc = TradeServerReturnCodeDescription(result.retcode);
		Print("Trailing stop: ",result.retcode," - ",errDesc,", Old SL: ",currentStop,", New SL: ",request.sl,", Bid: ",bid,", Ask: ",ask,", Stop Level: ",SymbolInfoInteger(pSymbol,SYMBOL_TRADE_STOPS_LEVEL));
		
		if(checkRes == CHECK_RETCODE_OK) return(true);
		else return(false);
	}
	else return(false);
}


// Trailing stop (points, hedging orders)
bool CTrailing::TrailingStop(ulong pTicket,int pTrailPoints,int pMinProfit=0,int pStep=10)
{
	if(PositionSelectByTicket(pTicket) == true && pTrailPoints > 0)
	{
		request.action = TRADE_ACTION_SLTP;
		request.position = pTicket;
		
		long posType = PositionGetInteger(POSITION_TYPE);
		double currentStop = PositionGetDouble(POSITION_SL);
		double openPrice = PositionGetDouble(POSITION_PRICE_OPEN);
		string symbol = PositionGetString(POSITION_SYMBOL);
		
		double point = SymbolInfoDouble(symbol,SYMBOL_POINT);
		int digits = (int)SymbolInfoInteger(symbol,SYMBOL_DIGITS);
		
		if(pStep < 10) pStep = 10;
		double step = pStep * point;
		
		double minProfit = pMinProfit * point;
		double trailStop = pTrailPoints * point;
		currentStop = NormalizeDouble(currentStop,digits);
		
		double trailStopPrice;
		double currentProfit;
		
		// Order loop
		int retryCount = 0;
		int checkRes = 0;
		
		do 
		{
			if(posType == POSITION_TYPE_BUY)
			{
				trailStopPrice = SymbolInfoDouble(symbol,SYMBOL_BID) - trailStop;
				trailStopPrice = NormalizeDouble(trailStopPrice,digits);
				currentProfit = SymbolInfoDouble(symbol,SYMBOL_BID) - openPrice;
				
				if(trailStopPrice > currentStop + step && currentProfit >= minProfit)
				{
					request.sl = trailStopPrice;
					bool sent = OrderSend(request,result);
				}
				else return(false);
			}
			else if(posType == POSITION_TYPE_SELL)
			{
				trailStopPrice = SymbolInfoDouble(symbol,SYMBOL_ASK) + trailStop;
				trailStopPrice = NormalizeDouble(trailStopPrice,digits);
				currentProfit = openPrice - SymbolInfoDouble(symbol,SYMBOL_ASK);
				
				if((trailStopPrice < currentStop - step || currentStop == 0) && currentProfit >= minProfit)
				{	
					request.sl = trailStopPrice;
					bool sent = OrderSend(request,result);
				}
				else return(false);
			}
			
			checkRes = CheckReturnCode(result.retcode);
		
			if(checkRes == CHECK_RETCODE_OK) break;
			else if(checkRes == CHECK_RETCODE_ERROR)
			{
				string errDesc = TradeServerReturnCodeDescription(result.retcode);
				Alert("Trailing stop: Error ",result.retcode," - ",errDesc);
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
		
		string errDesc = TradeServerReturnCodeDescription(result.retcode);
		Print("Trailing stop: ",result.retcode," - ",errDesc,", #",pTicket,", Old SL: ",currentStop,", New SL: ",request.sl,", Bid: ",SymbolInfoDouble(symbol,SYMBOL_BID),", Ask: ",SymbolInfoDouble(symbol,SYMBOL_ASK),", Stop Level: ",SymbolInfoInteger(symbol,SYMBOL_TRADE_STOPS_LEVEL));
		
		if(checkRes == CHECK_RETCODE_OK) return(true);
		else return(false);
	}
	
	else return(false);
}


// Trailing stop (price, hedging orders)
bool CTrailing::TrailingStop(ulong pTicket,double pTrailPrice,int pMinProfit=0,int pStep=10)
{
	if(PositionSelectByTicket(pTicket) == true && pTrailPrice > 0)
	{
		request.action = TRADE_ACTION_SLTP;
		request.position = pTicket;
		
		long posType = PositionGetInteger(POSITION_TYPE);
		double currentStop = PositionGetDouble(POSITION_SL);
		double openPrice = PositionGetDouble(POSITION_PRICE_OPEN);
		string symbol = PositionGetString(POSITION_SYMBOL);
		
		double point = SymbolInfoDouble(symbol,SYMBOL_POINT);
		int digits = (int)SymbolInfoInteger(symbol,SYMBOL_DIGITS);
		
		if(pStep < 10) pStep = 10;
		double step = pStep * point;
		double minProfit = pMinProfit * point;
				
		currentStop = NormalizeDouble(currentStop,digits);
		pTrailPrice = NormalizeDouble(pTrailPrice,digits);
		
		double currentProfit;
		
		int retryCount = 0;
		int checkRes = 0;
		
		double bid = 0, ask = 0;
		
		do 
		{
			if(posType == POSITION_TYPE_BUY)
			{
				bid = SymbolInfoDouble(symbol,SYMBOL_BID);
				currentProfit = bid - openPrice;
				if(pTrailPrice > currentStop + step && currentProfit >= minProfit) 
				{
					request.sl = pTrailPrice;
					bool sent = OrderSend(request,result);
				}
				else return(false);
			}
			else if(posType == POSITION_TYPE_SELL)
			{
				ask = SymbolInfoDouble(symbol,SYMBOL_ASK);
				currentProfit = openPrice - ask;
				if((pTrailPrice < currentStop - step || currentStop == 0) && currentProfit >= minProfit)
				{
					request.sl = pTrailPrice;
					bool sent = OrderSend(request,result);
				}
				else return(false);
			}
			
			checkRes = CheckReturnCode(result.retcode);
		
			if(checkRes == CHECK_RETCODE_OK) break;
			else if(checkRes == CHECK_RETCODE_ERROR)
			{
				string errDesc = TradeServerReturnCodeDescription(result.retcode);
				Alert("Trailing stop: Error ",result.retcode," - ",errDesc);
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
		
		string errDesc = TradeServerReturnCodeDescription(result.retcode);
		Print("Trailing stop: ",result.retcode," - ",errDesc,", #",pTicket,", Old SL: ",currentStop,", New SL: ",request.sl,", Bid: ",bid,", Ask: ",ask,", Stop Level: ",SymbolInfoInteger(symbol,SYMBOL_TRADE_STOPS_LEVEL));
		
		if(checkRes == CHECK_RETCODE_OK) return(true);
		else return(false);
	}
	else return(false);
}



// Break even stop
bool CTrailing::BreakEven(string pSymbol,int pBreakEven,int pLockProfit=0)
{
	if(PositionSelect(pSymbol) == true && pBreakEven > 0)
	{
		request.action = TRADE_ACTION_SLTP;
		request.symbol = pSymbol;
		
		long posType = PositionGetInteger(POSITION_TYPE);
		double currentSL = PositionGetDouble(POSITION_SL);
		double openPrice = PositionGetDouble(POSITION_PRICE_OPEN);
		
		double point = SymbolInfoDouble(pSymbol,SYMBOL_POINT);
		int digits = (int)SymbolInfoInteger(pSymbol,SYMBOL_DIGITS);
		
		double breakEvenStop;
		double currentProfit;
		
		int retryCount = 0;
		int checkRes = 0;
		
		double bid = 0, ask = 0;
		
		do 
		{
			if(posType == POSITION_TYPE_BUY)
			{
				bid = SymbolInfoDouble(pSymbol,SYMBOL_BID);
				breakEvenStop = openPrice + (pLockProfit * point);
				currentProfit = bid - openPrice;
				
				breakEvenStop = NormalizeDouble(breakEvenStop, digits);
				currentProfit = NormalizeDouble(currentProfit, digits);
				
				if(currentSL < breakEvenStop && currentProfit >= pBreakEven * point) 
				{
					request.sl = breakEvenStop;
					bool sent = OrderSend(request,result);
				}
				else return(false);
			}
			else if(posType == POSITION_TYPE_SELL)
			{
				ask = SymbolInfoDouble(pSymbol,SYMBOL_ASK);
				breakEvenStop = openPrice - (pLockProfit * point);
				currentProfit = openPrice - ask;
				
				breakEvenStop = NormalizeDouble(breakEvenStop, digits);
				currentProfit = NormalizeDouble(currentProfit, digits);
				
				if((currentSL > breakEvenStop || currentSL == 0) && currentProfit >= pBreakEven * point)
				{
					request.sl = breakEvenStop;
					bool sent = OrderSend(request,result);
				}
				else return(false);
			}
			
			checkRes = CheckReturnCode(result.retcode);
		
			if(checkRes == CHECK_RETCODE_OK) break;
			else if(checkRes == CHECK_RETCODE_ERROR)
			{
				string errDesc = TradeServerReturnCodeDescription(result.retcode);
				Alert("Break even stop: Error ",result.retcode," - ",errDesc);
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
		
		string errDesc = TradeServerReturnCodeDescription(result.retcode);
		Print("Break even stop: ",result.retcode," - ",errDesc,", SL: ",request.sl,", Bid: ",bid,", Ask: ",ask,", Stop Level: ",SymbolInfoInteger(pSymbol,SYMBOL_TRADE_STOPS_LEVEL));
		
		if(checkRes == CHECK_RETCODE_OK) return(true);
		else return(false);
	}
	else return(false);
}


// Break even stop (hedging orders)
bool CTrailing::BreakEven(ulong pTicket,int pBreakEven,int pLockProfit=0)
{
	if(PositionSelectByTicket(pTicket) == true && pBreakEven > 0)
	{
		request.action = TRADE_ACTION_SLTP;
		request.position = pTicket;
		
		long posType = PositionGetInteger(POSITION_TYPE);
		double currentSL = PositionGetDouble(POSITION_SL);
		double openPrice = PositionGetDouble(POSITION_PRICE_OPEN);
		string symbol = PositionGetString(POSITION_SYMBOL);
		
		double point = SymbolInfoDouble(symbol,SYMBOL_POINT);
		int digits = (int)SymbolInfoInteger(symbol,SYMBOL_DIGITS);
		
		double breakEvenStop;
		double currentProfit;
		
		int retryCount = 0;
		int checkRes = 0;
		
		double bid = 0, ask = 0;
		
		do 
		{
			if(posType == POSITION_TYPE_BUY)
			{
				bid = SymbolInfoDouble(symbol,SYMBOL_BID);
				breakEvenStop = openPrice + (pLockProfit * point);
				currentProfit = bid - openPrice;
				
				breakEvenStop = NormalizeDouble(breakEvenStop, digits);
				currentProfit = NormalizeDouble(currentProfit, digits);
				
				if(currentSL < breakEvenStop && currentProfit >= pBreakEven * point) 
				{
					request.sl = breakEvenStop;
					bool sent = OrderSend(request,result);
				}
				else return(false);
			}
			else if(posType == POSITION_TYPE_SELL)
			{
				ask = SymbolInfoDouble(symbol,SYMBOL_ASK);
				breakEvenStop = openPrice - (pLockProfit * point);
				currentProfit = openPrice - ask;
				
				breakEvenStop = NormalizeDouble(breakEvenStop, digits);
				currentProfit = NormalizeDouble(currentProfit, digits);
				
				if((currentSL > breakEvenStop || currentSL == 0) && currentProfit >= pBreakEven * point)
				{
					request.sl = breakEvenStop;
					bool sent = OrderSend(request,result);
				}
				else return(false);
			}
			
			checkRes = CheckReturnCode(result.retcode);
		
			if(checkRes == CHECK_RETCODE_OK) break;
			else if(checkRes == CHECK_RETCODE_ERROR)
			{
				string errDesc = TradeServerReturnCodeDescription(result.retcode);
				Alert("Break even stop: Error ",result.retcode," - ",errDesc);
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
		
		string errDesc = TradeServerReturnCodeDescription(result.retcode);
		Print("Break even stop: ",result.retcode," - ",errDesc,", #",pTicket,", SL: ",request.sl,", Bid: ",bid,", Ask: ",ask,", Stop Level: ",SymbolInfoInteger(symbol,SYMBOL_TRADE_STOPS_LEVEL));
		
		if(checkRes == CHECK_RETCODE_OK) return(true);
		else return(false);
	}
	else return(false);
}