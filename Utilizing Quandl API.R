#### Time Series AR + MA Practice 
library(Quandl)

# Imports equity price data of shares of Boston Scientific trading on the Hong Kong Stock Exchange
BostonScientific <- Quandl("HKEX/06669", type = "xts")

BostonScientific_DF <- as.data.frame(BostonScientific)
write.csv(BostonScientific_DF,"/Users/christopherthomson/Desktop")

# Computing ACF of daily ending prices
acf(BostonScientific[,1], lag.max = 50, type = "correlation") 
# Daily closing prices are heavily correlated 

# Recall that log(P_t/P_t-1) = ln(P_t) - ln(P_t-1), therefore diff(log(.)) will give 
BS_logrets <- na.omit(diff(log(BostonScientific[,1])))
names(BS_logrets) <- c("Logarithmic Returns")


acf(BS_logrets, lag.max = 35, type = "correlation")
# No autocorrelation beyond the 2nd lag, barely significant with the first.

Box.test(BS_logrets, lag = 20, type = "Ljung-Box")
Box.test(BS_logrets, lag = 20, type = "Box-Pierce")
# Both tests indicate autocorrelation occurs at at least lag value between 1 and 20
mean(BS_logrets) # -.0004 
sd(BS_logrets) # 0.03989 over past 20 years of daily trading on Hong Kong Exchange


ts.plot(BS_logrets)
abline(h = mean(BS_logrets), col = "red")
abline(h = max(BS_logrets), col = "blue")
max(BS_logrets)

partac <- pacf(BS_logrets, lag.max = 25, plot = FALSE) ### Note statistical signifcance at lag 10

arima(BS_logrets, order = c(10,0,0))

((0.941*1.05) - 1)/-0.4

