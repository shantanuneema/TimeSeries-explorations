
library('astsa')
library('CombMSC')
library('forecast')
library('TSA')

sundata <- read.csv("SN_y_tot_V2.0.csv",header = FALSE,sep=";")
trdata <- ts(sundata[1:311,2])
tsdata <- ts(sundata[-c(1:311),2],start = 312, end = 316)

# Test of Stationarity
adf.test(sundata[,2])
# Results are satisfactory

!is.null(tbats(sundata[,2])$seasonal)
# No seasonality

# Let's plot ACF and PACF
acf2(trdata,50)


ts.plot(diff(sundata[,2],differences = 2))
acf2(diff(trdata,differences = 3),50)

f <- 6
round(sarima.for(trdata,n.ahead = f,2,0,1,0,0,0,0)$pred,1)
mat <- NULL
for (I in 1:7) {
  for (J in 1:7) {
    future <- sarima.for(trdata,n.ahead = f,I,0,J,0,1,0,11)
    lines(tsdata,lty="dashed",lwd=2,col='blue')
    RMSE <- sqrt(sum((tsdata - future$pred)^2)/f)
    mat <- rbind(mat,c(RMSE,I,J,mean(future$se)))
    cat(I,J,"\n")
  }
}
plot(mat[,1],type='l')
mat[order(mat[,1])[1:8],]

data1 <- diff(trans,differences = 2) # d = 2
data2 <- diff(data1,11) # D = 1

answ <- NULL
for (K in 1:4){
  li <- 295+(K-1)*5
  trdata <- ts(sundata[1:li,2])
  tsdata <- ts(sundata[-c(1:li),2],start = li+1)
  f <- 316-li
  for (I in 1:3) {
    for (J in 1:2) {
      model <- sarima(trdata,I,2,J,1,1,0,11)
      future <- sarima.for(trdata,n.ahead = f,I,2,J,1,1,0,11)
      lines(tsdata,lty="dashed",lwd=2,col='blue')
      RMSE <- sqrt(sum((tsdata - future$pred)^2)/6)
      answ <- rbind(answ,c(I,J,K,RMSE))
    }
  }
}

answ

lamb <- BoxCox.lambda(sundata[,2])
trans <- BoxCox(trdata,0.5)
adf.test(trans)


hist(trans,breaks=50)
hist(trdata,breaks=25)

acf2(trans,100)
d1 <- diff(trans,11)
acf2(d1,100)
d2 <- diff(d1,100)

model <- arima(trdata,order=c(9,0,0))
future <- sarima.for(trdata,n.ahead = 5,2,0,1,0,0,0,0)
lines(tsdata,lty="dashed",lwd=2,col='blue')
sarima(trdata,8,0,0,0,0,0,0)


