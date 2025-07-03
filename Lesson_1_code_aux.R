# log 
v<-c(1,10,100,200)
w<-log(v)
v
w


indprod <- read.table("Data/ind.dat")
indprod <- ts(indprod,start=1949,frequency=12)
plot(indprod,xlab="Year",ylab="Industrial Production")

### Movie: decompose 
movie <- TRUE
delay <- 1
window <- 20
n <- length(indprod)/12 # number of years
n
if(movie) {
  for(t in 1:(n-window +1))
  {
    Sys.sleep(delay)
    subsamp <- indprod[((t-1)*12+1):((t-1+window)*12)] # extraction directe 
    newyear <- 1948 + t
    plot(ts(subsamp,start=newyear,frequency=12),ylab="")
  } }

# auto reg 

model2 <- lm(pop[-1] ~ pop[-n]) # same length variables Xt et Xt-1
summary(model2)
plot(ts(model2$residuals))


cor(pop[-1],pop[-n])
plot(pop[-1],pop[-n],xlab="X Past",ylab="X Present")
