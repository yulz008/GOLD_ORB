# Introduction

GOLD_ORB is an Expert Advisor for XAUUSD (GOLD) and uses 1HR TF written in MQL5.

## About the Project

This project aims to develop a working trading bot (EA) which utilizes price action as buy/sell signal. This project also caters the basics and structure of creating a
trading bot making it as a guide for beginners on creating their own trading bot. The project is designed so that it is customizable for developers to add new features
and strategies.


## Strategy

We will first start with our strategy, GOLD_ORB uses Open Range Breakout strategy to generate buy and sell signals (for opening position long/short)

"The opening range is high and low for a given period after the market opens. This period is generally the first 30 or 60 minutes of trading. " The opening hour of the 
market is associated with big trading volumes and volatility. This time of the trading
session provides many trading opportunities. In this way, traders use the opening range to set the entry points and to predict and forecast the price action for the
day.

Once the final range is established we just need to wait for the price to breakout/breakdown on its range. Buy long for breakouts and Short Sell for breakdown, these 
are the signals to be generated. 

For commodities such as GOLD(XAU/USD) market open is the start of trading hour which is 1:02 (server time).

Since GOLD_ORB uses 1HR TF, 1 hour after the market opens, the "initial" range will be calculated. Which is nothing but the high and low of the 1st candle of 
day on the H1 timeframe. This will become the initial low(support) and high(resistance) of the range. This is not yet the "final" range, the high and low will be 
updated on the incoming candles. The EA is configured to wait for at least 3 "candle composition" or upcoming candles must consolidate within the initial range so  the 
range can be considered "final".  Otherwise the EA will just update the new support or resistance on the incoming candles once it establishes a new low or new high. 
Once the minimum "candle compostion" or consolidation is meet the range is now considered final and if  a breakout/breakdown happens on the "final" range it will 
generatea buy/sell signal. The "candle compostion" can be adjusted depending on the user prefrence this will be discussed later in detail.



### Strategy Logic

//execute only at the beginning of each new candle
1. 1HR after market open, get initial range low and high of the candle. (initial support and resistance)
2. at succeding candles, update range high and low else do nothing 
3. generate buy/sell signal at breakout/breakdown on final range.
4. repeat



// 1 HR Timeframe of XAU/USD, the vertical lines signifies the beginning of trading hours (1:02 server time)
![EA Working](https://user-images.githubusercontent.com/117939069/206975769-a6170f95-bf78-4efc-b967-593984513111.gif)



![image](https://user-images.githubusercontent.com/117939069/206984297-14a01ae2-9019-4608-84c9-70bf85ed253a.png)

### Strategy Class
By using OOP, a class is created with the strategy stated above, the class "price_action" is located  at 

GOLD_ORB_v1/GOLD_ORB/Include/price_action.mqh

The class Price_Action contains two significant public members: new_candle_check2 and  Open_Range_Breakout. Take note that these members can only be process under the 
event handler of OnTick() and OnInit() since these functions reads and anaylyze the data tick by tick.  

**new_candle_check2**

The new_candle_check2 function as the name implies checks if the present tick came from the opening of a new candle. It will return a boolean type, if it detects the 
present of new candle it will return "1" and "0" otherwise. This is an important function since we want our EA to know the beginning of each candle. Since our EA will 
only analyze and execute its logic at the beginning of each new candle. We first wait for the candle to conclude before doing any processing to eliminate market noise.


**Open_Range_Breakout**

This is the heart of our EA, as this is function for the Stategy Logic stated above. This will process only at the beggining of each candle. It will return an integer 
type. This will generate the buy and sell signals of our EA on opening position. It will return "11" for "buy-long" signal, "10" for "short-sell" signal and "0" if no 
signal generated.

From here you can create your own class depending on your strategy. If you use price action a lot like me, then this class will be very helpfull to start
with.


## Algorithm Structure
![image](https://user-images.githubusercontent.com/117939069/201953919-a6dd6e05-4918-42ac-b340-c46b8edf67a4.png)


## Features

### Trade Management
The Trade Management include adjustable to TakeProfit, StopLoss, MaxTradePerDay
  - TakeProfit
  - StopLoss
  - MaxTradePerDay
  - LongPosition
  - ShortPosition
 
### Trail Management
The Trail management includes the trail module, this is not customizable on the inputs, since as per my backtest, this setting yields the greater result compared to
other input combination. Once the gain of the position reaches to 700 points this will activate the trail to 100points at minimum.
  
### Risk Management
One of the very important feature of this bot is the Risk Management, this automatically calculates the position size relative to the Account size.
  - MaxEquityDrawdownPercent
  - RiskPerTrade
  - FixedLot
 
### Indicators
For this EA, only 1 indicator is available which is the price action. The "Candle Composition" can be adjusted, depending on the consoldidation period the user prefer

### Advanced Equity Monitoring Module

## Backtest

![image](https://user-images.githubusercontent.com/117939069/201954432-3b38daf8-e183-4cfa-8384-3f88e2d8fc1c.png)

![image](https://user-images.githubusercontent.com/117939069/201954567-2c30a2c0-ec65-4bf6-83ac-59ddee0d188b.png)


## How to use the EA

## Revisions

## Disclaimer



