# YahooVIXCalc
Calculates a real-time "Volatility index" for any stock or ETF with option chains on Yahoo Finance, following the calculation methodology
of the CBOE VIX index. Note that the methodology is not exactly the same, but the differences should be negligible, and these numbers are
not meant to be interpreted in a vacuum or replicate any published numbers. For instance, the CBOE uses a seconds methodology to calculate 
times to maturity, whereas I a just using day counts, along with the same interest rate for both near and next terms contracts, etc. 

The code offers a real time VIX-type calculation for the stock or ETF chosen as "ticker." However, note that you need to pick an 
underlying with enough contracts such that there is a reasonable close cluster around 30 days. For instance, suppose at the time of writing,
(5/31/2018), I would like to see what a small-cap VIX looks like. Maybe the Vanguard ETF with ticker VB comes to mind as I hold this in my portfolio.
However, the closest VB options to 30 days are June 15th and July 20th, which are quite far and actually fail the CBOE's methodology
that the contracts should be >23 days and <37 days. Instead, I can choose IWM, which is a smaller ETF but has more actively traded contracts, including 
for this little example a June 29th marturity and a July 6th maturity, which are perfect for the calculation.

Example:
Say I'd like to get a real-time VIX calc for large caps, small caps, and gold. I could call

[MFIV_LC] = YahooMFIVCalc('SPY');
[MFIV_SC] = YahooMFIVCalc('IWM');
[MFIV_GLD] = YahooMFIVCalc('GLD');

**Note that you could technically use the code for indices, i.e. make the ticker '^...' however, you need to be very careful here
and probably manually check that the right contracts are available on Yahoo. For instance, the option chains on Yahoo for ^SPX are not
any good, i.e. they clearly aren't the SPX traditional product contracts available from the CBOE. There are 2 overwhelmingly likely reasons for
error with running the code - connection errors and an underlying that does not have good chains available on Yahoo.
