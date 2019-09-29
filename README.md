# EquityOptionPricing_OOP
European Option Pricing and Greeks Calculation for Non dividend paying stocks


#This program performs 3 functions. It performs Equity Vanilla Option Pricing, generates greeks and performs volatility smile fitting based on user inputted Vol/Strike grid.

--------------------------------------------------------------------------------------------------------------------------------
Function 1: European Option Pricing for Non dividend paying stocks
--------------------------------------------------------------------------------------------------------------------------------
The module involves implementation of Black Scholes pricing for user inputted Vanilla Equity Option trade economics - Spot Price (S), Strike Price (K), Volatility (V), Risk Free Rate of Return (r), Time to Expiry (T) and Option Type ('P' for Put, 'C' for Call). The option premium is calculated as per below.

Call Option Price = S * N(d1) - K * exp(-r * T) * N(d2)

Put Option Price = K * exp(-r * T) * N(-d2) - S * N(-d1)

where d1 = ( ln(S / K) + (r + (V * V)/2) * T) / (V * sqrt(T))

d2 = ( ln(S / K) + (r - (V * V)/2) * T) / (V * sqrt(T)) = d1 - V * sqrt(T)

N(x) represents cumulative probability function for a standardised normal distribution

The code uses an approximation of the cumulative probability function for the standard normal distribution. Please note that volatility and risk free rate of return should be input as fractions (Ex: 0.05 for 5%). The code calculates the option premium and prints on console.

------------------------------------------------------------------------------------------------------------------------------------
Function 2: Vanilla Option Greeks Calculation
------------------------------------------------------------------------------------------------------------------------------------
The module calculates the following option greeks.
1> Delta: Change in option premium with respect to unit change in Spot Price (S) 
2> Gamma: Change in option delta with respect to unit change in Spot Price (S) 
3> Vega : Change in option premium with respect to unit change in Volatility (V) 
4> Theta: Change in option premium with respect to unit change in time to expiry (T) 
5> Rho : Change in option premium with respect to unit change in risk free rate of return (r)

-----------------------------------------------------------------------------------------------------------------------------------
Function 3: Volatility Smile Fitting
-----------------------------------------------------------------------------------------------------------------------------------
The module aims to arrive at best fit for Volatility Surface based on user inputted Volatility/ Strike grid and returns a pointer to the array of calibrated parameters - Convexity, Skew and Constant. The module expresses smile surface (volatility as a function of Strike) as Volatility = Convexity * Strike^2 + Skew * Strike + Constant and then tries to arrive at best fit by changing the parameters Convexity, Skew, Constant by a fixed amount (tweak amount) in each iteration. Standard error is calculated for each iteration. The iteration for which the standard error is minimal, is considered as the best fit for the volatility smile and can be used to quote implied volatility for range of strikes. Please note that module calibrates well within a close range of ATM Strike. However, for deep out of/ in the money strikes, the fit may not be that optimal meaning the quadratic behaviour breaks at higher strikes. However, the module does fairly well for close to ATM Strikes and for lesser dispersion of Volatility across the strike grid.
The user inputs Implied Volatility against various different strikes in the system. The module arrives at the best fit as per the above and prints the calibrated parameters on the console along with the standard error. The tweak amount in each iteration can be lowered to arrive at better fits however at the cost of increased computation time. Hence the user needs to balance between accuracy of the fit and compute time to arrive at best fit.

---------------------------------------------------End of Documentation-----------------------------------------------------------------
