# Introduction

GOLD_ORB is an Expert Advisor for XAUUSD (GOLD) on 1HR TF, written in MQL5. 

MQL5 is a C++ based programming language.




## About the Project

This project aims to develop a working trading bot (EA) which utilizes price action as buy/sell signal. This project also caters the structure of creating a
trading bot, making it as a guide for beginners on creating their own trading bot. The project is designed so that it is customizable for developers to add new features and strategies. 

>If you're a complete beginner, before you proceed it is best for you to read the book about the  basic of MQL5:  
>Expert Advisor Programming for MetaTrader 5 by Andrew Young (see reference [2])

However, if you are already familiar with MQL5 or with other trading bot framework, you can proceed on the next outline. Otherwise, if things begin to sound confusing then you might want to go back and familiarize first with the code structure and basic classes of MQL5. 

## Strategy

We will first start with our strategy, GOLD_ORB uses an Open Range Breakout strategy in 1HR timeframe to generate buy and sell signals (for opening position long/short)

>The opening range is high and low for a given period after the market opens. This period is generally the first 30 or 60 minutes of trading. The opening hour of the market is associated with big trading volumes and volatility. This time of the trading session provides many trading opportunities. In this way, traders use the opening range to set the entry points and to predict and forecast the price action for the day. [1]

Once the final range is established we just need to wait for the price to breakout/breakdown on its range. Buy long for breakouts and Short Sell for breakdown, these 
are the signals to be generated. 

For commodities such as GOLD(XAUUSD) the market open around 1:02 (server time) which is the start of trading hour.

Since GOLD_ORB uses 1HR TF, 1 hour after the market opens, the "initial" range will be calculated. This is nothing but the high and low of the 1st candle of day on the H1 timeframe. This will become the initial low(support) and high(resistance) of the range. This is not yet the "final" range, the high and low will be updated on the incoming candles. The EA is configured to wait for at least 3 "candle composition" or upcoming candles must consolidate within the initial range so the range can be considered "final".  Otherwise, the EA will just update the new support or resistance on the incoming candles if it has a new low or new high.

Once the minimum "candle composition" is meet the range is now considered final. And if a breakout/breakdown happens on the "final" range it will generate buy/sell signal. The "candle composition" can be adjusted depending on the user preference but by default it set as “3”.



### Strategy Logic

//execute only at the beginning of each new candle
1. 1 hour after market open, get initial range low and high of the candle. (Initial support and resistance)
2. at succeeding candles update range high and low, else do nothing.
3. on final range, generate buy/sell signal at breakout/breakdown.
4. repeat



// H1 Timeframe of XAUUSD, the vertical lines signifies the beginning of trading hours (1:02 server time)
![EA Working](https://user-images.githubusercontent.com/117939069/206975769-a6170f95-bf78-4efc-b967-593984513111.gif)


![image](https://user-images.githubusercontent.com/117939069/206984297-14a01ae2-9019-4608-84c9-70bf85ed253a.png)

### Strategy Class
By using OOP, a class is created with the strategy stated above, the class "price_action" is located at 

GOLD_ORB_v1/GOLD_ORB/Include/price_action.mqh

The class Price_Action contains two significant public members: new_candle_check2 and  Open_Range_Breakout. Take note that these members can only be process under the 
event handler of OnTick() and OnInit() since these functions reads and analyze the data tick by tick.  

**new_candle_check2**

The new_candle_check2 function as the name implies checks if the present tick came from the opening of a new candle. It will return a Boolean type. If it detects the 
present of new candle it will return "true" and "false" otherwise. This is an important function since we want our EA to know the beginning of each candle. Our EA will 
only analyze and execute its logic at the beginning of each new candle. We first wait for the candle to conclude before doing any processing to eliminate market noise.


**Open_Range_Breakout**

This is the heart of our EA, this is the function for the Strategy Logic discuss earlier. This will process only at the beginning of each candle. It will return an integer type. This will generate the buy and sell signals of our EA for opening a position. It will return "11" for buy-long signal, "10" for short-sell signal and "0" if no signal is generated.

From here you can create your own class depending on your strategy. If you use price action a lot like me, then this class will be very helpful to start
with. You can also used time series indicators such as the  Moving Average and RSI.  MQL5 has its built-in class for this time series indicator. 



##  
## Algorithm Structure
![image](https://user-images.githubusercontent.com/117939069/207205822-ced98100-df1f-4d9e-ab33-8c27225311fb.png)

Let us first cover the backbone and algorithm structure of GOLD_ORB. This structure is very common to every MQL program. We have the 
preprocessor directives followed by the inputs and global variables, then the object creation and lastly we have the event handlers.

**Preprocessor Directives**
  - this are used to set properties, define constants, include files and import functions/libraries. This are declared on top of each MQL program. For this EA the 
   #include directive is used to call files containing several class and structure that will be used in the program.

**Inputs and Global Varibles**
   - inputs are parameters accessible to the user once the program is compiled. This value can be customized depending on the user inputs. The global variables are 
   additional variables used to interact with the input and output of different classes and functions within the program.
   
**Classes** 

   -  All includes files(mqh) are located here GOLD_ORB_v1/GOLD_ORB/Include
   - This section creates the object from a specific class. It is important to note that before we can create the object, we must first include the file where the class definition in located. The following are the classes that will be used on the EA. Some of these classes and functions are found on the book, and are reused on this EA to save time from developing. (denoted in *)  
   
 
      - *CTrade - for excuting orders
      - *CTrailing - for trail stops
      - *CiMA - a class for time series indicator: Moving average
      - Price_Action - a class for open range breakout strategy
      - CTradeVirtual - for executing orders virtually
      - CTrailingVirtual - for virtual trail stops
      - VirtualTradeInfo - Virtual Environment: information on the virtualtrades


**Event handler**
  
  Event handlers are the method by which and MQL5 runs. It executes whenever a certain event occurs.[2] 
  For this EA two event handlers where used:
     - OnInit() - runs once the program starts, this is were we can initialize our function or objects
     - OnTick() - runs whenever the new tick is received by the EA. 
 
 **OnTick()**   
 
 This contains the code every time a price change occur. Most of the functionality of this EA is located under this event handler. For organizing purposes ,the block of codes are categorize depending on its functionality which includes the Risk Management Group, Trail Management Group, Indicators Group and the Trade Management Group.
 
 - Risk Management Group -  includes 2 subcategory modules: the "Essential" and "Advance" modules.
  
    Essential
     - Dynamic Position Sizing Module - If this module is enabled, in accordance with the input, the position size per trade varies dependently on the account size. 
     - Drawdown Monitoring Module - this monitors the account balance, once the maximum set drawdown is reach by the EA, it will automatically stop trading. It will remain stop until the user reset the EA and resumes it back again.

     Advaced(optional)
     
     - Equity Slope detection - detects if the equity balance is trending upward or downward. 
     - Losing Streak Module - detects losing streak and stop the EA on trading. The losing streak number is set by the user.
     - Virtual Equity Monitoring Module - this is enabled by default. This copies the trades of the EA and execute it virtually. If the EA stopped from trading due to losing streak or drawdown hit the EA will continue to trade virtually. The EA can be configured to resume trading in real account together with Equity Slope Detection and Losing Streak Module configuration.
   
   
 - Trail Mangement Group 
    - The Trail management includes the trail module, this is not customizable on the inputs, since as per my backtest on XAAUSD, this setting yields the greater result compared to other input combination. Once the gain of the position reaches to 700 points this will activate the trail to 100points at minimum. Uses the class CTrailing.
    

- Indicators Group

    - The indicator group for this EA includes the Price_Action class only as discuss earlier, other indicators such time series indicators (e.g. Moving Averages, RSI, and MACD) can be used together with the Price_Action class as long as the new time series class is included and the object is created prior (eg. CiMA)
     
- Trade Management Group

    - This is where the buying and selling took place to open the position on the EA. It uses the CTrade Class discuss on the Algorithm structure. Take note that when we open a position the TakeProfit and Stoploss order is sent together with our "opening a position order". This module is heavily dependent on the output of the Indicator group. Without the signal generated by the indicators, no buy or sell orders will be sent to the server. 


### GOLD_ORB EA Complete Logic

//After Compiling the source code
1. Read inputs
2. Oninit():Initialize objects, functions, and global variables
3. Ontick():
    - If Enabled, Compute Position Size
    - If Enabled, Check Losing Streak 
    - If Enabled, Check if max drawdown is hit, stop trading if true
    - If Enabled, check if trailstop is hit on existing position
    - If new candle, Check Price_Action for signal
    - If Price_Action has signal, execute orders
    - All actions above are reflected in the virtual environment
 
 
##  
## EA Features and Inputs

![image](https://user-images.githubusercontent.com/117939069/207226726-5f667083-e449-4149-ab4f-bb1827c0ac23.png)


### SymbolInformation
Contains the StartOFTradingHour_ServerTime for XAUUSD the default value for this is "1". 

### Trade Management
The Trade Management include customizable TakeProfit, StopLoss, MaxTradePerDay
  - TakeProfit      - integer input: default "1200"
  - StopLoss        - integer input: default "400"
  - MaxTradePerDay  - integer input: default "2" 
    - The maximum trade to open a position per day is 2: Long-1 and Short-1. If only 1 trade per day is desired, input must set to "1"
 
  - LongPosition    - default "True"
  - ShortPosition   - default "True"
 
### Trail Management
To enable Trail it must be set to "True", "False" if otherwise.
  
### Risk Management
One of the very important features of this bot is the Risk Management, the EA can automatically calculate the position size relative to the account size.
  - MaxEquityDrawdownPercent
      - integer input: default 10 
      - A "10" on the input means 10%, this means a 10% drawdown from the account will stop the EA from trading.
  - RiskPerTradePercent      
      - integer input: default 1
      - A "1" on the input means that for every trade, the EA will only risk 1% of the present account balance
  - FixedLot  
      - double input: default 0.1
      - if the user want to set a fixed lot per trade, then this could be adjusted here, RiskPerTradePercent must be disable or set to "0"
 
### Advanced Equity Monitoring Module
  - SlopeDetection - default "False"
  - LossStreakCounter - integer input: default: 0
 
### Indicators
For this EA, only one indicator is available: price action(PA). The "candle composition" can be adjusted, depending on the consolidation period the user prefer. Default value is "3". 


### Broker and Hedging Mode

The broker I used to test and run this EA is ICMarkets MT5 (Real Account Server). Hedging mode is enabled on the broker which means that different open position will be place on each buy/sell order on the same trading symbol (eg XAU/USD GOLD). You will notice once the EA is running 3-7+ position is open with different direction but under the same symbol. This is called hedging mode which is the opposite of netting mode


## Backtest and Results

![image](https://user-images.githubusercontent.com/117939069/207229689-e2f05309-f80a-4fa1-86f3-731bcd2b735b.png)
![image](https://user-images.githubusercontent.com/117939069/207229753-0e2b3829-2c2b-4a06-ba0c-622e8e6b7da9.png)


## How to use the EA

### Backtesting
1. After copying the project folder to MQL5 directory. Perform recalibration test to ensure its functionality. See Test Cases folder
2. EA can be used for backtesting other TimePeriod, Timeframes, Pairs and Etc. this are all customizable. Please note the starting trading time for different pairs. Different pairs open on different server time.

### Forward Testing
1. Attach the MQL program to chart
2. Load default inputs
3. You will notice on the upper right side of the chart the name of the program if EA is already running.
![image](https://user-images.githubusercontent.com/117939069/207310994-b30402b1-078c-495b-a964-523818ab933b.png)


## Reference
[1] “Opening range breakouts - what is orb trading strategy? - elearnmarkets.” [Online]. Available: https://www.elearnmarkets.com/blog/how-to-trade-opening-range-breakouts/. [Accessed: 13-Dec-2022]. 

[2] A. R. Young, Expert Advisor Programming for MetaTrader 5: Creating Automated Trading Systems in the MQL5 language. Nashville, TN: Edgehill Publishing, 2013. 

[3] “MQL5 reference - how to use algorithmic/automated trading language for metatrader 5,” MQL5. [Online]. Available: https://www.mql5.com/en/docs. [Accessed: 13-Dec-2022]. 


## Disclaimer

This code is just meant to be used for learning. In no way do I promise profitable trading outcomes. Do not risk money that you can't afford to lose because the authors and any affiliates assume no responsibility for your trading results. This strategy DO NOT come with ANY warranty, thus there may be flaws in the code. Investments are risky by nature! Future outcomes cannot be predicted based on past performance!

