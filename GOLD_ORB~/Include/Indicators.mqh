//+------------------------------------------------------------------+
//|                                                   Indicators.mqh |
//|                                              Playground Inc 2021 |
//|                                             https://www.mql5.com |
//+------------------------------------------------------------------+
#property copyright "Playground Inc 2021"
#property link      "https://www.mql5.com"
#property version   "1.00"


#define MAX_COUNT 100


//+------------------------------------------------------------------+
//| Base Class                                                       |
//+------------------------------------------------------------------+

class CIndicator
{
	protected:
		int handle;
		double main[];
		
	public:
		CIndicator(void);
		double Main(int pShift=0);
		void Release();
		virtual int Init() { return(handle); }
};

CIndicator::CIndicator(void)
{
	ArraySetAsSeries(main,true);
}

double CIndicator::Main(int pShift=0)
{
	CopyBuffer(handle,0,0,MAX_COUNT,main);
	double value = NormalizeDouble(main[pShift],_Digits);
	return(value);
}

void CIndicator::Release(void)
{
	IndicatorRelease(handle);
}


//+------------------------------------------------------------------+
//| Moving Average                                                   |
//+------------------------------------------------------------------+

/*

CiMA MA;

sinput string MA;		// Moving Average
input int MAPeriod = 10;
input ENUM_MA_METHOD MAMethod = 0;
input int MAShift = 0;
input ENUM_APPLIED_PRICE MAPrice = PRICE_CLOSE;

MA.Init(_Symbol,_Period,MAPeriod,MAShift,MAMethod,MAPrice);

MA.Main()

*/

class CiMA : public CIndicator
{
	public:
		int Init(string pSymbol,ENUM_TIMEFRAMES pTimeframe,int pMAPeriod,int pMAShift,ENUM_MA_METHOD pMAMethod,ENUM_APPLIED_PRICE pMAPrice);
};

int CiMA::Init(string pSymbol,ENUM_TIMEFRAMES pTimeframe,int pMAPeriod,int pMAShift,ENUM_MA_METHOD pMAMethod,ENUM_APPLIED_PRICE pMAPrice)
{
	handle = iMA(pSymbol,pTimeframe,pMAPeriod,pMAShift,pMAMethod,pMAPrice);
	return(handle);
}


//+------------------------------------------------------------------+
//| RSI                                                              |
//+------------------------------------------------------------------+

/*

CiRSI RSI;

sinput string RS;	// RSI
input int RSIPeriod = 8;
input ENUM_APPLIED_PRICE RSIPrice = PRICE_CLOSE;

RSI.Init(_Symbol,_Period,RSIPeriod,RSIPrice);

RSI.Main()

*/




class CiRSI : public CIndicator
{
	public:
		int Init(string pSymbol, ENUM_TIMEFRAMES pTimeframe, int pRSIPeriod, ENUM_APPLIED_PRICE pRSIPrice);
};

int CiRSI::Init(string pSymbol, ENUM_TIMEFRAMES pTimeframe, int pRSIPeriod, ENUM_APPLIED_PRICE pRSIPrice)
{
	handle = iRSI(pSymbol,pTimeframe,pRSIPeriod,pRSIPrice);
	return(handle);
}


//+------------------------------------------------------------------+
//| Stochastic                                                       |
//+------------------------------------------------------------------+

/*

CiStochastic Stoch;

sinput string STO;	// Stochastic
input int KPeriod = 10;
input int DPeriod = 3;
input int Slowing = 3;
input ENUM_MA_METHOD StochMethod = MODE_SMA;
input ENUM_STO_PRICE StochPrice = STO_LOWHIGH;

Stoch.Init(_Symbol,_Period,KPeriod,DPeriod,Slowing,StochMethod,StochPrice);

Stoch.Main()
Stoch.Signal()

*/




//+------------------------------------------------------------------+
//| Bollinger Bands                                                  |
//+------------------------------------------------------------------+

/*

CiBollinger Bands;

sinput string BB;		// Bollinger Bands
input int BandsPeriod = 20;
input int BandsShift = 0;
input double BandsDeviation = 2;
input ENUM_APPLIED_PRICE BandsPrice = PRICE_CLOSE; 

Bands.Init(_Symbol,_Period,BandsPeriod,BandsShift,BandsDeviation,BandsPrice);

Bands.Upper()
Bands.Lower()

*/




//+------------------------------------------------------------------+
//| Blank Indicator Class Templates                                  |
//+------------------------------------------------------------------+

/* 

Replace _INDNAME_ with the name of the indicator.
Replace _INDFUNC_ with the name of the correct technical indicator function.
Add appropriate input parameters (...) to Init() function.
Rename Buffer1(), Buffer2(), etc. to something user-friendly.
Add or remove buffer arrays and functions as necessary.



// Single Buffer Indicator

class Ci_INDNAME_ : public CIndicator
{
	public:
		int Init(string pSymbol, ENUM_TIMEFRAMES pTimeframe, ... );
}; 


int Ci_INDNAME_::Init(string pSymbol,ENUM_TIMEFRAMES pTimeframe, ... )
{
	handle = _INDFUNC_(pSymbol,pTimeframe, ... );
	return(handle);
}



// Multi-Buffer Indicator

class Ci_INDNAME_ : public CIndicator
{
	private:
		double buffer1[];
		double buffer2[];
		
	public:
		int Init(string pSymbol,ENUM_TIMEFRAMES pTimeframe, ... );
		double Buffer1(int pShift=0);
		double Buffer2(int pShift=0);
		Ci_INDNAME_();
}; 


Ci_INDNAME_::Ci_INDNAME_()
{
	ArraySetAsSeries(buffer1,true);
	ArraySetAsSeries(buffer2,true);
}


int Ci_INDNAME_::Init(string pSymbol,ENUM_TIMEFRAMES pTimeframe,...)
{
	handle = _INDFUNC_(pSymbol,pTimeframe,...);
	return(handle);
}


double Ci_INDNAME_::Buffer1(int pShift=0)
{
	CopyBuffer(handle,1,0,MAX_COUNT,buffer1);
	double value = NormalizeDouble(buffer1[pShift],_Digits);
	return(value); 
} 


double Ci_INDNAME_::Buffer2(int pShift=0)
{
	CopyBuffer(handle,1,0,MAX_COUNT,buffer2);
	double value = NormalizeDouble(buffer2[pShift],_Digits);
	return(value); 
} 


*/