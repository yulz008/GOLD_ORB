# Introduction

GOLD_ORB is an Expert Advisor for XAUUSD (GOLD) and uses 1HR TF written in MQL5.

## About the Project

This project aims to develop a working trading bot (EA) which utilizes price action as buy/sell signal. This project also caters the basics and structure of creating a
trading bot making it as a guide for beginners on creating their own trading bot. The project is designed so that it is customizable for developers to add new features
and strategies.


## Strategy

We will first start with our strategy, GOLD_ORB uses Open Range Breakout strategy to generate buy and sell signals (for opening position long/short)

>"The opening range is high and low for a given period after the market opens. This period is generally the first 30 or 60 minutes of trading. " The opening hour of 
the market is associated with big trading volumes and volatility. This time of the trading session provides many trading opportunities. In this way, traders use the 
opening range to set the entry points and to predict and forecast the price action for the
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
![image](https://user-images.githubusercontent.com/117939069/207205822-ced98100-df1f-4d9e-ab33-8c27225311fb.png)

Let us first cover the backbone and algorithm structure of GOLD_ORB. This structure is very common to every MQL program. We have the 
preprocessor directives followed by the inputs and global variables, then the object creation and lastly we have the event handlers.

**Preprocessor Directives**
  - this are used to set properties, define constants, include files and import functions/libraries. This are declared on top of each MQL program. For this EA the 
   #include directive is used to call files containing several class and structure that will be used in the program.

**Inputs and Global Varibles**
   - inputs are parameters accessibale to the user once the program is compiled. This value can be customized depending on the user inputs. The global variables are 
   additional variables used to interact input and output of different classes and functions
   
**Classes** 
   - This section creates the object from a specific class. It is important to note that before we can create the object we must first include the file where the class defintion in located. The following are the classes that will be used on the EA. Some of this classes are found on the book.
   
      - CTrade - for excuting orders
      - CTrailing - for trail stops
      - Price_Action - a class for the strategy
      - CiMA - a class for time series indicator: Moving average
      - CTradeVirtual - for executing orders virtually
      - CTrailingVirtual - for virtual trail stops
      - VirtualTradeInfo - Virtual Environment: information on the virtualtrades

**Event handler**
  -  Event handlers are the method by which and MQL5 runs. It executes whenever a certain event occurs.[2] 
  -  For this EA two event handlers where used:
     - OnInit() - runs once the program starts, this is were we can initialize our function or objects
     - OnTick() - runs whenever the new tick is received by the EA. 
 
 **OnTick()**   
 
 This contains the code everytime a price change occur. Most of the functionality of this EA is located under this event handler. For organizing purposes the block of codes are categorize depending on its functionality which includes the Risk Management Group, Trail Management Group, Indicators Group and the Trade Management Group.
 
 - Risk Management Group -  includes 2 subcategory modules the "Essential" and "Advance" modules.
  
    Essential
     - Dynamic Position Sizing Module - If this module is enabled, the position size per trade varies dependently on the acount size.
     - Drawdown Monitoring Module - this monitors the acount balance, once the maximum drawdown is reach by the EA, it will automaticaly stops trading. Until the user modifies it back.

     Advaced(optional)
     
     - Equity Slope detection - detects if the equity balance is trending upward or downward. 
     - Losing Streak Module - detects losing streak and stop the EA on trading. The losing streak number is set by the user.
     - Virtual Equity Monitoring Module - this is enables by default. This copy the trades of the EA and execute it virtually. If the EA stopped from trading due to lossing streak or drawdown hit the EA will continue to trade virtually. The EA can be configured to resume trading in real account together with Equity Slope Detection and Losing Streak Module. 
   
   
 - Trail Mangement Group 
    - The Trail management includes the trail module, this is not customizable on the inputs, since as per my backtest, this setting yields the greater result compared to other input combination. Once the gain of the position reaches to 700 points this will activate the trail to 100points at minimum. Uses the class CTrail.
    

- Indicators Group

    - The indicator group for this EA includes the Price_Action class only as discuss above, other indicators such time series indicators (e.g. Moving Averages, RSI, and MACD) can be used together with the Price_Action class as long as the new time serires class is included and the object is created prior eg. CiMA
     
- Trade Management Group

    - This is where the buying and selling took place to open the position on the EA. It uses the CTrade Class discuss on the Algorithm structure. Take note that when we open a position the TakeProfit and Stoploss order is sent together with our "opening a position order". This module is heavily dependent on the output of the Indicator group. Without the signal generated by the indicators, no buy or sell orders will be sent to the server. 
    
## Features

### Trade Management
The Trade Management include adjustable to TakeProfit, StopLoss, MaxTradePerDay
  - TakeProfit
  - StopLoss
  - MaxTradePerDay
  - LongPosition
  - ShortPosition
 
### Trail Management
To enable Trail
  
### Risk Management
One of the very important feature of this bot is the Risk Management, this automatically calculates the position size relative to the Account size.
  - MaxEquityDrawdownPercent 
  - RiskPerTradePercent
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



