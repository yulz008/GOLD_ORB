Introduction

GOLD_ORB is an Expert Advisor for XAUUSD (GOLD) and uses 1HR TF written in MQL5.

About the Project

This project aims to develop a working trading bot (EA) which utilizes price action as buy/sell signal. This project also caters the basics and structure of creating a
trading bot making it as a guide for beginners on creating their own trading bot. The project is designed so that it is customizable for developers to add new features
and strategies.


Strategy

GOLD_ORB uses Open Range Breakout strategy to generate buy and sell signals (for opening position long/short)

"The opening range is high and low for a given period after the market opens. This period is generally the first 30 or 60 minutes of trading. It is one most important
chart patterns to make money in the stock market." The opening hour of the market is associated with big trading volumes and volatility. This time of the trading
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
generatea buy/sell signal. The "candle compostion" can be adjusted depending on the user prefrence. This will be discussed later



This is best illustrated below:

1. 1HR after market open, get initial range low and high of the candle. (initial support and resistance)
2. at succeding candles, update range high and low else do nothing 
3. generate buy/sell signal at breakout/breaokdown on final range.
4. Repeat




[INSERT GIF that shows the working EA updates the support and resistnace]


Range Breakout: Long Position

![image](https://user-images.githubusercontent.com/117939069/201953535-3fc70a14-5b7f-4648-a80c-0160accc31aa.png)


Range Breakdown: Short Position

![image](https://user-images.githubusercontent.com/117939069/201955988-b487401c-3458-40eb-b4ae-1de67ef04795.png)


Algorithm Structure
![image](https://user-images.githubusercontent.com/117939069/201953919-a6dd6e05-4918-42ac-b340-c46b8edf67a4.png)


Backtest Result

![image](https://user-images.githubusercontent.com/117939069/201954432-3b38daf8-e183-4cfa-8384-3f88e2d8fc1c.png)

![image](https://user-images.githubusercontent.com/117939069/201954567-2c30a2c0-ec65-4bf6-83ac-59ddee0d188b.png)

