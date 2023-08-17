##### ARMA Basics by Christopher Thomson

## Basics in constructing ARMA models using zoo and tseries R packages 

library(zoo)
library(tseries)

read.csv("^IXIC.csv")

colnames(Nasdaq) <- c("Date","Open","High","Low","Close","Adjusted Close","Volume")

nsdf <- read.csv("^IXIC.csv") ### Extracts S&P500 Data
S_t <- nsdf$Adj.Close         ### Daily Ending Prices
LN_t <- log(S_t[1:249]/S_t[2:250]) ### Computes Daily log returns days 1-251 of trading

S_ts <- ts(S_t, start = 1, frequency = 1) ### Creates time series of raw prices 
ln_ts <- ts(LN_t, start = 1, frequency = 1) ### Creates time series of log returns


ferrari <- read.csv("RACE.csv") ### Read in Ferrari data
ferr_S <- ferrari$Adj.Close ### Raw closing price data
ln_f <- log(ferr_S[1:249]/ferr_S[2:250]) ### Convert to log-returns
ln_fts <- ts(ln_f, start = 1, frequency = 1) ### Create time series of log returns



plot(S_ts) ### Price plot
plot(ln_ts) ### Log-return time series plot
plot(ln_fts) ### Log-tertun time series plot of Ferrari

acf(ln_ts, lag.max = 20)
### Computes autocorrelation plot
acf(ln_fts, lag.max = 30)

Box.test(ln_ts, lag = 5, type = "Ljung-Box") ### Simultaneous test for multiple lagged values
Box.test(ln_ts, lag = 50, type = "Ljung-Box") ### For some CI of p-value, null hypothesis is not one k-lag value indicates correlation
### A test statistic with a p-value less than or equal to the threshold indicates reject null, and that on more lags exhibits autocorrelation

Box.test(ln_fts, lag = 25, type = "Ljung-Box") ### Computes ac in ferrari series


adf.test(ln_ts) ### Dickey Fuller Unit Root Test, indicates stationarity via characteristic polynomial root test (tseries package)


### Forecasting Outward
install.packages("forecast")
library(forecast)
auto.arima(ln_ts, max.P = 25, max.Q = 20, ic = "aic")

auto.arima(ln_fts, max.p = 25, max.q = 25, ic = "aic") ### Estimates best fitting ARIMA
mu_ferrari <- mean(ln_fts) ### mean of TS


arima(ln_fts, order = c(2,0,0))
fit_ferrari <- arima(ln_fts, order = c(2,0,0))
forecast <- predict(fit_ferrari, 30)
print(forecast$pred)

plot(forecast$pred)


