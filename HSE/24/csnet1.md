# Code from Random 

## Current value of a changing quantity in the computer

Different quantities in the computer are unpredictably changing and could serve as a source of randomness. Most of them are very difficult to access from high-level programming languages such as R. A good candidate is the system time
```
> t <- Sys.time()
> t
[1] "2022-08-14 23:17:59 CEST"
> cat(t,"\n")
1660511880 
> class(t)
[1] "POSIXct" "POSIXt" 
> s <- Sys.time()
> s
[1] "2022-08-14 23:19:42 CEST"
> cat(s,"\n")
1660511982 
> e <- Sys.time()
> e
[1] "2022-08-14 23:20:12 CEST"
> cat(e,"\n")
1660512013 
> e-s
Time difference of 30.37726 secs
> as.numeric(s)
[1] 1660511982
> as.numeric(e)
[1] 1660512013
> as.numeric(e)-as.numeric(s)
[1] 30.37726
> T <- "1970-01-01 00:00:14"
> t <- as.POSIXct(T,tz="GMT")
> as.numeric(t)
[1] 14
```
In R, the Sys.time() counts the seconds from January 1, 1970. For many possible applications, the changes in seconds are too slow.

Recently two R libraries were developed that provide access to time at the level of nanoseconds: [[https://github.com/eddelbuettel/nanotime|nanotime]] and [[https://dirk.eddelbuettel.com/code/rcpp.cctz.html|RcppCCTZ]]. They are based on the 64 bit representation of integers.
```
> install.packages("RcppCCTZ")
> install.packages("nanotime")
```
Let's try!
```
> library(nanotime)
> library(bit64)
> t <- Sys.time()
> as.numeric(t)
[1] 1660442526
> v <- nanotime(t)
> as.integer64(v)
integer64
[1] 1660442525739342000
```
The last three places are 0 - Windows provides time stamps at the level of microseconds.

We could use Sys.time() values for throwing dice in an interactive game. Changing fast are last 3 or 4 places.
```
> k <- 6
> repeat{
+   a <- readline("Throw - press ENTER (yes) or s and ENTER (stop)")
+   if(a=="s") break
+   t <- nanotime(Sys.time()); n <- as.integer64(t)
+   d <- as.integer((n%/%1000)%%k+1); print(n); print(d)
+ }
Throw - press ENTER (yes) or s and ENTER (stop)
integer64
[1] 1660517375129460000
[1] 1
Throw - press ENTER (yes) or s and ENTER (stop)
integer64
[1] 1660517379681713000
[1] 6
Throw - press ENTER (yes) or s and ENTER (stop)
integer64
[1] 1660517381783372000
[1] 3
Throw - press ENTER (yes) or s and ENTER (stop)
integer64
[1] 1660517384964824000
[1] 3
Throw - press ENTER (yes) or s and ENTER (stop)
integer64
[1] 1660517387146669000
[1] 2
Throw - press ENTER (yes) or s and ENTER (stop)
integer64
[1] 1660517391795778000
[1] 5
...
```
A more serious application is a random initialization of some parameter in a program.

## Wichmann and Hill's generator

```
> WichHill <- function(n){
+    output <- numeric(n)
+    seed <- get(".WHseed",env=.WH)
+    x <- seed[1]; y <- seed[2]; z <- seed[3]
+    for(i in 1:n){
+       x <- (171*x) %% 30269
+       y <- (172*y) %% 30307
+       z <- (170*z) %% 30323
+       output[i] <- (x/30269 + y/30307 + z/30323) %% 1.0
+    }
+    assign(".WHseed",c(x,y,z),env=.WH)
+    output
+ }
> 
> .WH <- new.env(); assign(".WHseed",1:3,env=.WH)
> WichHill(10)
 [1] 0.03381877 0.77754189 0.05273525 0.74462407 0.49036219 0.98285437
 [7] 0.80915099 0.71338138 0.80102091 0.98958603
> r <- WichHill(10000)
> hist(r,prob=TRUE)
> r <- WichHill(100000)
> hist(r,prob=TRUE)
```


## Basic generators in R

```
RNGkind()
set.seed(2018)
.Random.seed[1:6]
runif(6)
.Random.seed[1:6]
set.seed(2018)
runif(6)
RNGkind("Super")
RNGkind()
.Random.seed[1:6]
```

```
> RNGkind()
[1] "Mersenne-Twister" "Inversion"       
> set.seed(2018)
> .Random.seed[1:6]
[1]        403        624  -79855522 -803040953  -27922212 -118944723
> runif(6)
[1] 0.33615347 0.46372327 0.06058539 0.19743361 0.47431419 0.30104860
> .Random.seed[1:6]
[1]        403          6 -881521081  953668654 1212651912  902275211
> set.seed(2018)
> runif(6)
[1] 0.33615347 0.46372327 0.06058539 0.19743361 0.47431419 0.30104860
> RNGkind("Super")
> RNGkind()
[1] "Super-Duper" "Inversion"  
> .Random.seed[1:6]
[1]         402  1572731790 -1300846921          NA          NA          NA
```

## Uniform and Bernoulli distribution

```
> random <- function() {return(runif(1,0,1))}
> dice <- function(n=6) {return(1+trunc(n*random()))}
> Bernoulli <- function(p) {if(random()<=p) return(1) else return(0)}
> n <- 10; s <- numeric(n)
> for(i in 1:n) s[i] <- random()
> s
 [1] 0.29062721 0.70364315 0.25889652 0.22473082 0.69339126 0.13130646
 [7] 0.04117978 0.53640553 0.39368591 0.56187369
> n <- 20; s <- numeric(n)
> for(i in 1:n) s[i] <- dice(3)
> s
 [1] 1 3 2 2 2 3 2 2 2 2 3 3 2 1 1 1 1 2 2 3
> for(i in 1:n) s[i] <- Bernoulli(1/3)
> s
 [1] 0 1 0 0 0 1 0 0 0 0 1 0 0 0 0 0 1 1 0 0
> table(s)
s
 0  1 
15  5 
> n <- 10000; s <- numeric(n)
> for(i in 1:n) s[i] <- Bernoulli(1/3)
> table(s)
s
   0    1 
6672 3328 
```

## Geometric distribution


```
> geometric <- function(p){
+   if (p>=1) return(1)
+   if (p<=0) return(Inf)
+   return(trunc(log(1-random())/log(1-p))+1)
+ }
> n <- 10; s <- numeric(n)
> for(i in 1:n) s[i] <- geometric(1/3)
> s
 [1] 3 2 2 1 4 4 2 5 1 3
> for(i in 1:n) s[i] <- geometric(0.1)
> s
 [1] 21 11 12  2 12 10 17  2  2  6
> n <- 100000; s <- numeric(n)
> for(i in 1:n) s[i] <- geometric(0.1)
> t <- table(s)
> x <- as.numeric(names(t))
> plot(x,t,pch=16,cex=0.7)
``` 

## Tabelaric distribution


```
> tabelaR <- function(p){
+   r <- random(); k <- 0;
+   while (r >= 0) {k <- k+1; r <- r - p[k]}
+   return(names(p)[k])
+ }
>
> f <- c(1,2,3,2,1); names(f) <- c("mon","tue","wed","thu","fri")
> p <- f/sum(f)
> s <- vector("character",10)
> for (i in 1:10) s[i] <- tabelaR(p)
> s
 [1] "fri" "tue" "wed" "wed" "tue" "wed" "tue" "wed" "tue" "wed"
> n <- 100000; s <- vector("character",n)
> for (i in 1:n) s[i] <- tabelaR(t)
> table(s)
s
  fri   mon   thu   tue   wed 
11127 11114 22403 22015 33341 
```

## Poisson distribution


```
> PoissonRnd <- function(lambda){
+   k <- 0; p <- exp(-lambda); r <- random()-p
+   while(r > 0) {k <- k+1; p <- p*lambda/k; r <- r-p}
+   return(k)
+ }
> 
> n <- 10000; s <- numeric(n)
> for(i in 1:n) s[i] <- PoissonRnd(0.3)
> table(s)
s
   0    1    2    3    4    5 
7381 2251  333   33    1    1 
> z <- rpois(n,0.3)
> table(z)
z
   0    1    2    3    4 
7430 2204  324   38    4 
```


## Uniform continuous distribution

```
> uniformRnd <- function(a,b) {a+random()*(b-a)}
> n <- 100000; s <- numeric(n)
> for(i in 1:n) s[i] <- uniformRnd(75,100)
> hist(s,prob=TRUE)
```

## (Negative) exponential distribution

```
> exponRnd <- function(lambda) {return(-log(1-random())/lambda)} 
> n <- 100000; s <- numeric(n)
> for(i in 1:n) s[i] <- exponRnd(1.5)
> hist(s,prob=TRUE); lines(density(s),col="red")
```

## Cauchy distribution


```
> cauchyRnd <- function(a=1,b=0) {tan(pi*(random()-0.5))/a + b}
> n <- 20; s <- numeric(n)
> for(i in 1:n) s[i] <- cauchyRnd()
> s
 [1]  1.4517597 -2.3231091 -1.2913811 -1.6924863  8.4182520 -4.6988607
 [7]  1.4703372  0.1900819  0.2629615  2.9184602 -0.1944157  6.6547234
[13] -0.3956465 -0.4921238  0.4071398 -8.2433149 -0.6530219  0.2205136
[19]  0.1589820 -0.5175836
```

## von Neumann rejection method
Triangle { (0,0), (1,1), (2,0) }
```
> vonNeumann <- function(a,b,G,g,...){
+   repeat { s <- a + random()*(b-a)
+     if (G*random()<=g(s,...)) return(s)
+   }
+ }
> triang <- function(x){if (x>2) return(0); if (x<0) return(0);
+   if (x<1) return(x) else return(2-x)
+ }
> 
> n <- 100000; s <- numeric(n)
> for(i in 1:n) s[i] <- vonNeumann(-0.5,2.5,1,triang)
> hist(s,prob=TRUE); lines(density(s),col="red",lw=2)
```
## Normal or Gaussian distribution 

```
> gaussRnd <- function(st=1,m=0,s=1){
+   if (st==0){new <- TRUE; .Gauss <<- list(m=m,s=s,n=new,r=0)}
+   else {m <- .Gauss$m; s <- .Gauss$s; new <- .Gauss$n }
+   if (new) {p <- sqrt(-2*log(random())); q <- 2*pi*random()
+     x <- p*sin(q); .Gauss$r <<- p*cos(q)
+   } else  x <- .Gauss$r
+   .Gauss$n <<- !new
+   return(m+s*x)
+ }
> n <- 100000; s <- numeric(n); t <- gaussRnd(0,175,10)
> for(i in 1:n) s[i] <- gaussRnd()
> hist(s,prob=TRUE); lines(density(s),col="red",lw=2)
```

## Multidimensional normal distribution

``` 
> multinormal <- function(T){return(t(t(T)%*%rnorm(dim(T)[1])))}
> data(longley); R <- cor(longley); T <- chol(R); pairs(longley)
> n <- 2000; s <- NULL;
> for (i in 1:n) s <- rbind(s,multinormal(T))
> pairs(s)
```

## Random 3D direction


```
> dir3D<-function(){return(c(2*pi*random(),acos(1-2*random())))}
> vector3D <- function(dir){
+   return( c( sin(dir[2])*cos(dir[1]), sin(dir[2])*sin(dir[1]), cos(dir[2]) ))
+ }
> points3D <- function(n){
+   cat(file="sphere.net",c("*vertices ",n,"\n"))
+   for (i in 1:n) {
+     cat(file="sphere.net",i,i,vector3D(dir3D()),"\n",append=TRUE)
+   }
+ }
> setwd("C:/Users/batagelj/Documents/papers/2018/moskva/NetR/nets/test")
> points3D(500)
```
 
To view the file ''sphere.net'' use Pajek.


## Levy flights

```
> direction <- function() {phi <- 2*pi*random(); return(c(cos(phi),sin(phi)))}
>
> # Brown motion
> x <- c(0,0); n <- 200; s <- matrix(0,nrow=n,ncol=2); s[1,] <- x
> for (i in 2:n) {x <- x + direction(); s[i,] <- x}
> plot(s,pch=16,col="red"); lines(s)
>
> # Levy flights
> x <- c(0,0); n <- 500; s <- matrix(0,nrow=n,ncol=2); s[1,] <- x
> for (i in 2:n) {x <- x + cauchyRnd(3)*direction(); s[i,] <- x}
> plot(s,pch=16,col="red"); lines(s)
> 
> x <- c(0,0); n <- 500; s <- matrix(0,nrow=n,ncol=2); s[1,] <- x
> for (i in 2:n) {x <- x + cauchyRnd(0.1)*direction(); s[i,] <- x}
> plot(s,pch=16,col="red"); lines(s)
```

## Uniform distribution in multidimensional ellipsoid

```
> elipsoid <- function(T){
+   m <- dim(T)[1]; r <- random()^(1/m)
+   x <- rnorm(m); y <- x*r/sqrt(sum(x^2))
+   return(t(t(T) %*% y))
+ }
> 
> R <-
+   c( 4.0000000, 1.1395235, 1.775876, 0.7753723,
+      1.1395235, 4.0000000, 1.065415, 0.7881692,
+      1.7758761, 1.0654147, 4.000000, 1.5841046,
+      0.7753723, 0.7881692, 1.584105, 4.0000000 )
> m <- 4; dim(R) <- c(m,m); T <- chol(R)
> n <- 500; s <- matrix(0,nrow=n,ncol=m)
> for (i in 1:n) s[i,] <- elipsoid(T)
> pairs(s)
>
> R <- matrix(0,4,4); diag(R) <- c(1,2,3,4); T <- chol(R)
> n <- 500; s <- matrix(0,nrow=n,ncol=m)
> for (i in 1:n) s[i,] <- elipsoid(T)
> pairs(s)
```

## Shuffle

```
> shuffle <- function(x){
+   n <- length(x)
+   for ( i in n:2 ){j <- dice(i)
+     t <- x[i]; x[i] <- x[j]; x[j] <- t }
+   return(x)
+ }
> x <- c("⇐","⇑","⇒","⇓","⇖","⇗","⇘","⇙","●") 
> for(i in 1:5) cat(x <- shuffle(x),"\n")
⇒ ⇘ ⇗ ⇖ ⇙ ● ⇓ ⇑ ⇐ 
⇐ ⇖ ⇒ ⇓ ⇑ ⇗ ● ⇙ ⇘ 
● ⇑ ⇐ ⇓ ⇙ ⇒ ⇘ ⇗ ⇖ 
⇖ ⇐ ⇒ ⇗ ● ⇑ ⇓ ⇘ ⇙ 
⇗ ⇒ ⇓ ● ⇐ ⇖ ⇙ ⇑ ⇘ 
```
 
## Sample

```
> sampleN <- function(n,m){
+   k <- 0; s <- integer(m)
+   for(i in 0:n) if ((n-i)*random() < m-k) {
+     k <- k+1; s[k] <- i
+     if (k==m) return(s)
+   }
+ }
> for(i in 1:5) cat(sampleN(100,15),"\n")
2 14 19 25 26 33 39 45 53 58 62 73 87 92 94 
6 18 20 21 31 33 34 36 40 49 54 69 73 74 98 
9 12 27 41 44 45 48 49 53 54 71 77 78 82 97 
12 24 33 36 37 44 49 52 54 58 67 72 73 79 84 
1 5 15 26 49 53 54 58 61 74 76 84 92 93 97 
```

## Number π 


```
> k <- 0; n <- 2000; u <- integer(n)
> for(i in 1:n){
+   if (random()^2+random()^2 < 1) k <- k+1;
+   u[i] <- k }
> p <- u/seq(n)*4
> plot(p,pch=20,cex=0.8)
> lines(c(-10,n),c(pi,pi),col="red")
```

## Cube

Average distance between two random points in a unit cube.

{{ru:hse:snet:pics:answercube.png}}

[[http://mathworld.wolfram.com/HypercubeLinePicking.html|MathWorld]]
```
> n <- 100000; t <- 0.6617
> x <- matrix(runif(3*n),byrow=FALSE,nrow=n)
> y <- matrix(runif(3*n),byrow=FALSE,nrow=n) 
> a <- cumsum(d<-sqrt(rowSums((x-y)^2)))/seq(n)
> plot(1:n,a,pch=16,cex=0.5)
> lines(c(-n/20,n),c(t,t),col="red",lw=2)
> 
> I <- a[n]; lambda <- 3
> sigma <- sqrt(sum((d-I)^2)/n) 
> delta <- lambda*sigma/sqrt(n)
> cat("I =",I,"within error",delta,"with probability 0.997\n") 
I = 0.6619123 within error 0.001121217 with probability 0.997
> plot(1:n,a,pch=16,cex=0.5,ylim=c(0.66,0.665))
> lines(c(-n/20,n),c(t,t),col="red",lw=2)
```
## Sandokan 

```
> sandokan <- function(n,m=-1,r=100){
+   if (m<0) m <- n
+   s <- integer(r)
+   for ( i in 1:r ){
+     album <- rep(TRUE,n)
+     different <- 0; cards <-0
+     while (different < m) {
+       k <- dice(n); cards <- cards + 1
+       if (album[k]) {
+         album[k] <- FALSE; different <- different + 1
+       }
+     }
+     s[i] <- cards
+   }
+   return(s)
+ }
> 
> sandokanT <- function(n,m=-1){
+   if (m<0) m <- n
+   return(n*sum(1/((n-m+1):n)))
+ }
> 
> picture1 <- function(t,s,...){
+   n <- length(s)
+   mu <- cumsum(s)/1:n
+   sigma <- sqrt(cumsum(s^2)/1:n - mu^2)
+   plot(s,pch=16,cex=0.2,...)
+   lines(1:n,mu,col="red")
+   lines(1:n,mu+sigma,col="blue")
+   lines(1:n,mu-sigma,col="blue")
+   abline(h=t,col="goldenrod")
+ }
> 
> picture2 <- function(s,breaks=50,...){
+   n <- length(s); t <- NULL 
+   for(x in min(s):max(s)) {t <- c(t,5*n*dnorm(x,mean(s),sqrt(var(s))))}
+   hist(s,breaks=breaks,...)
+   lines(min(s):max(s),t,col="red",lw=2)
+ }
>
> (teo <- sandokanT(400,350))
[1] 828.2897
>
> s <- sandokan(400,350,2000)
> head(s)
[1] 866 837 747 772 815 801
> picture1(teo,s)
> picture2(s)
>
> pdf("sandokan.pdf")
> picture1(teo,s,ylim=c(820,835))
> dev.off()
```

## Monte Carlo in statistics

Variates X1 and X2 have standard normal distribution. Estimate the expected value of their absolute difference.

```
> N <- 10000
> x1 <- rnorm(N); x2 <- rnorm(N); y <- abs(x1-x2)
> print(theta.hat <- mean(y))
[1] 1.14084
> print(se.theta <- sd(y)/sqrt(N))
[1] 0.008639172
> print(theta <- 2/sqrt(pi))
[1] 1.128379
> print(se <- sqrt((2-4/pi)/N))
[1] 0.008525025
>
> N <- 1000000
> x1 <- rnorm(N); x2 <- rnorm(N); y <- abs(x1-x2)
> print(theta.hat <- mean(y))
[1] 1.128559
> print(se.theta <- sd(y)/sqrt(N))
[1] 0.0008520579
> print(theta <- 2/sqrt(pi))
[1] 1.128379
> print(se <- sqrt((2-4/pi)/N))
[1] 0.0008525025
```

## Resampling

Approximating a distribution F by the sample
```
> n <- 20
> mu <- 175; sigma <- 10; a <- 130; b <- 220
> x <- sort(rnorm(n,mean=mu,sd=sigma))
> y <- (1:n)/n
> tit <- paste("N(",mu,",",sigma,"), n=",n,sep="")
> curve(pnorm(x,mean=mu,sd=sigma),a,b,lwd=2,main=tit,ylab="p")
> points(x,y,pch=16,col="red")
```

## Bootstrap

```
> library(bootstrap)
> plot(law82$LSAT,law82$GPA)
> points(law$LSAT,law$GPA,pch=16,col='Red')
> print(cor(law82$LSAT,law82$GPA))
[1] 0.7599979
> print(theta.hat <- cor(law$LSAT,law$GPA))
[1] 0.7763745
> N <- 2000
> n <- nrow(law)
> theta.star <- numeric(N)
> for(k in 1:N){
+   i <- sample(1:n,size=n,replace=TRUE)
+   L <- law$LSAT[i];G <- law$GPA[i]
+   theta.star[k] <- cor(L,G)
+ }
> hist(theta.star,prob=TRUE)
> print(theta.star.mean <- mean(theta.star))
[1] 0.7712286
> print(bias <- theta.star.mean - theta.hat)
[1] -0.005145919
> print(se.theta.star <- sd(theta.star))
[1] 0.1302567
```

<hr/>

[HSE](../2024.md); [Docs](docs.md)






