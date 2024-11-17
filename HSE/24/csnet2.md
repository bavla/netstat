# Code for Matrices




## Operations

```
> r <- 2
> '%Pf+%' <- function(a,b){return(min(a,b))}
> '%Pf*%' <- function(a,b){return((a^r+b^r)^(1/r))}
> 3 %Pf+% 4
[1] 3
> 3 %Pf*% 4
[1] 5
```

## Semiring library

[[https://github.com/bavla/semirings/|Github / semirings]]

### Library

```
> Sr.set <- function(sr="combinatorial",r=2){
+   if(sr=="combinatorial"){
+     Sr.zero <<- 0; Sr.one <<- 1; Sr.absorption <<- FALSE
+     '%(+)%' <<- function(a,b) return(a+b)
+     '%(.)%' <<- function(a,b) return(a*b)
+     Sr.star <<- function(a) if(a==1) return(Inf) else
+       if(a==Inf) return(Inf) else return(1/(1-a))          
+   } else if(sr=="shortpaths"){
+     Sr.zero <<- Inf; Sr.one <<- 0; Sr.absorption <<- TRUE
+     '%(+)%' <<- min
+     '%(.)%' <<- function(a,b) return(a+b)
+     Sr.star <<- function(a) return(0)           
+   } else if(sr=="logical"){
+     Sr.zero <<- 0; Sr.one <<- 1; Sr.absorption <<- TRUE
+     '%(+)%' <<- max
+     '%(.)%' <<- min
+     Sr.star <<- function(a) return(1)           
+   } else if(sr=="maxmin"){
+     Sr.zero <<- 0; Sr.one <<- Inf; Sr.absorption <<- TRUE
+     '%(+)%' <<- max
+     '%(.)%' <<- min
+     Sr.star <<- function(a) return(Inf)           
+   } else if(sr=="log"){
+     Sr.zero <<- Inf; Sr.one <<- 0; Sr.absorption <<- FALSE
+     '%(+)%' <<- function(a,b) return(-log(exp(-a)+exp(-b)))
+     '%(.)%' <<- function(a,b) return(a+b)
+   } else if(sr=="pathfinder"){
+     Sr.zero <<- Inf; Sr.one <<- 0; Sr.absorption <<- TRUE
+     '%(+)%' <<- min
+     '%(.)%' <<- function(a,b){return((a^r+b^r)^(1/r))}
+     Sr.r <<- r
+     Sr.star <<- function(a) return(0)           
+   } else {cat("unknown semiring\n"); return(NULL)}
+   Sr.type <<- sr
+ }
>  
> Sr.set()
>
> Sr.get <- function() return(Sr.sr)
>
> Sr.Zero <- function(n) return(matrix(Sr.zero,nrow=n,ncol=n))
>
> Sr.One <- function(n){ I <- Sr.Zero(n); diag(I) <- Sr.one; return(I)}
> 
> Sr.adapt <- function(A) {P <- A; P[P==0] <- Sr.zero; return(P)}
> 
> Sr.binary <- function(A) {B <- A; B[B!=Sr.zero] <- Sr.one; return(B)}
>
> Sr.select <- function(n,s) {S <- rep(Sr.zero,n); S[s] <- Sr.one; return(S)}
>
> Sr.Perm <- function(A,p) return(A[p,p])
>
> Sr.perminv <- function(p) {q <- p; q[p] <- 1:length(p); return(q)}
>
> '%[+]%' <- Sr.Plus <- function(A,B){
+   if(is.null(dim(A))) {na <- length(A); C <- numeric(na)
+     for(i in 1:na) C[i] <- A[i] %(+)% B[i]; return(C)}
+   na <- nrow(A); ma <- ncol(A) 
+   C <- matrix(0,nrow=na,ncol=ma)
+   for(i in 1:na) for(j in 1:ma) C[i,j] <- A[i,j] %(+)% B[i,j]
+   return(C)
+ }
> 
> '%[.]%' <- Sr.Times <- function(A,B){ vector <- FALSE
+   if(is.null(dim(A))) {A <- matrix(A,ncol=length(A),nrow=1); vector <- TRUE}
+   if(is.null(dim(B))) {B <- matrix(B,nrow=length(B),ncol=1); vector <- TRUE}
+   na <- nrow(A); ma <- ncol(A); nb <- nrow(B); mb <- ncol(B)
+   if(ma!=nb) {cat("incompatible matrices\n"); return(NULL)}
+   C <- matrix(Sr.zero,nrow=na,ncol=mb)
+   for(i in 1:na) for(j in 1:mb) { s <- Sr.zero
+     for(k in 1:ma) s <- s %(+)% (A[i,k] %(.)% B[k,j])
+     C[i,j] <- s
+   }
+   if(vector) return(as.vector(C)) else return(C)
+ }
>
> Sr.Power <- function(A,k){
+   n <- nrow(A); Mt <- Sr.One(n)
+   rownames(Mt) <- colnames(Mt) <- rownames(A)
+   if (k > 0) { i <- k; Ms <- A
+     repeat {
+       if ((i %% 2) == 1) { Mt <- Mt %[.]% Ms }
+       i <- i %/% 2; if (i == 0) break
+       Ms <- Ms %[.]% Ms
+     }
+   }
+   return(Mt)
+ }
>
> .spw <- function(i,k,M){
+   n <- nrow(M)
+   if (i == 1) return(list(A=M %[+]% Sr.One(n), B=M %[.]% M))
+   if (i %% 2 == 0) {
+     AB <- .spw(i-1,k,M)
+     list(A=AB$A %[+]% AB$B, B=if(i<k){AB$B %[.]% M}else{NULL})
+   } else {
+     AB <- .spw(i %/% 2,k,M)
+     P <- AB$B %[+]% Sr.One(n)
+     list(A=P %[.]% AB$A, B=if(i<k){AB$B %[.]% AB$B}else{NULL})
+   }
+ }
> 
> Sr.Sumpow <- function(A,k){ n <- nrow(A)
+   if (k == 0) return(Sr.One(n))
+   if (k == 1) return(Sr.One(n) %[+]% A)
+   return(.spw(k,k,A)$A)
+ }
>
> Sr.Closure <- function(A){
+   if(!Sr.absorption) {cat("Semiring is not absorptive\n"); return(NULL)}
+   na <- nrow(A); C <- A
+   for(k in 1:na) {for(i in 1:na) for(j in 1:na)
+     C[i,j] <- C[i,j] %(+)% (C[i,k] %(.)% Sr.star(C[k,k]) %(.)% C[k,j])
+     C[k,k] <- Sr.one %(+)% C[k,k] 
+   }
+   return(C)
+ }
> 
> Sr.save.net <- function(fnet,A){
+   net <- file(fnet,"w")
+   n <- nrow(A); rn <- rownames(A)
+   cat("*vertices",n,"\n",file=net)   
+   for (i in 1:n) cat(i," \"",rn[i],"\"\n",file=net,sep="")
+   cat("*matrix\n",file=net)
+   for (i in 1:n) cat(A[i,],"\n",file=net)
+   close(net)
+ }
```

### Examples


```
> setwd("C:/Users/batagelj/Documents/papers/2018/moskva/NetR/nets")
> source("https://github.com/bavla/semirings/raw/master/semirings.R")
> Sr.set("shortpaths")
> Sr.Zero(5)
> Sr.One(5)
> Sr.select(9,4)
> Sr.set("log")
> 3 %(+)% 4
> 1 %(+)% 1
> 100 %(+)% 100
> 1 %(+)% 100
> 0.5 %(+)% 0.6
> -3 %(+)% 4
> -3 %(+)% 0
> 0 %(+)% 0
> 100 %(+)% Inf
> Sr.set("shortpaths")
> a <- c(1,0,5,2); b <- c(3,2,0,4)
> cbind(a,b,c=a %[+]% b)
> cbind(a,b,c=a %[.]% b)
```

<img src="https://github.com/bavla/netstat/blob/master/HSE/24/pics/semi.png" width=400>

```
> w <- c(
+   0, 2, 3,  5, 0, 0,  0, 0, 0,
+   0, 0, 0,  2, 4, 0,  0, 0, 0,
+   0, 0, 0,  4, 0, 0,  3, 0, 0,
+   0, 0, 0,  0, 3, 0,  2, 0, 0,
+   0, 0, 0,  0, 0, 5,  0, 2, 0,
+   0, 0, 0,  1, 0, 3,  5, 0, 0,
+   5, 0, 0,  0, 0, 0,  0, 0, 0,
+   0, 0, 0,  0, 0, 0,  0, 0, 0,
+   0, 0, 0,  0, 0, 2,  0, 0, 0  )
> W <- matrix(w,byrow=TRUE,nrow=9)
> colnames(W) <- rownames(W) <- paste("v",1:9,sep="")
> Sr.save.net("semi.mat",W)
> Sr.set("logical")
> B <- Sr.binary(W)
> Sr.Power(B,21)
> Sr.set("combinatorial")
> B <- Sr.binary(W)
> Sr.Sumpow(B,5)
> Q <- Sr.One(9)
> (Q <- Q %[.]% B)
> (s <- Sr.select(9,4))
> s
[1] 0 0 0 1 0 0 0 0 0
> (s <- s %[.]% B)
[1] 0 0 0 0 1 0 1 0 0
> (s <- s %[.]% B)
[1] 1 0 0 0 0 1 0 1 0
> (s <- s %[.]% B)
[1] 0 1 1 2 0 1 1 0 0
> (s <- s %[.]% B)
[1] 1 0 0 3 3 1 4 0 0
> Sr.set("shortpaths")
> W <- Sr.adapt(matrix(w,byrow=TRUE,nrow=9))
> colnames(W) <- rownames(W) <- paste("v",1:9,sep="")
> P <- Sr.Closure(W)
> P
    v1  v2  v3  v4  v5  v6  v7 v8  v9
v1   0   2   3   4   6  11   6  8 Inf
v2   9   0  12   2   4   9   4  6 Inf
v3   8  10   0   4   7  12   3  9 Inf
v4   7   9  10   0   3   8   2  5 Inf
v5  13  15  16   6   0   5   8  2 Inf
v6   8  10  11   1   4   0   3  6 Inf
v7   5   7   8   9  11  16   0 13 Inf
v8 Inf Inf Inf Inf Inf Inf Inf  0 Inf
v9  10  12  13   3   6   2   5  8   0
```
<html><!--
==== Canibals ====

Two missionaries and two cannibals reach the left bank of a river simultaneously and all of them want to cross the river with the aid of the single available boat which carries no more than two persons. On no bank of the river may the cannibals outnumber the missionaries otherwise the former jeopardize the latter.

<code>
> # Canibals
> # states (m,s) numbers on the left bank: 1 - initial state, 7 - final state
> Can <- matrix(c(
+   0, 1, 1, 1, 1, 0, 0,
+   0, 0, 1, 1, 0, 1, 0,
+   0, 0, 0, 0, 0, 0, 1,
+   0, 0, 0, 0, 0, 1, 1,
+   0, 0, 0, 0, 0, 1, 1,
+   0, 0, 0, 0, 0, 0, 1,
+   0, 0, 0, 0, 0, 0, 0 ),
+   byrow=TRUE,nrow=7)
> colnames(Can) <- rownames(Can) <- 
+   c("(2,2)","(2,1)","(2,0)","(1,1)","(0,2)","(0,1)","(0,0)")
</code>
--></html>

[[notes:net:allp|All paths]] in Python.

## Betweenness 

[[https://github.com/bavla/semirings/|Github / Betweenness]]

```
> mat.geodesics <- function(m)
+ { n <- nrow(m)
+   md <- m; md[m==0] <- Inf
+   mc <- m; mc[m>0] <- 1
+   for(k in 1:n) { for(u in 1:n){ for(v in 1:n){
+     dst <- md[u,k] + md[k,v]
+     if(md[u,v] >= dst) {
+       cnt <- mc[u,k]*mc[k,v];
+       if (md[u,v] == dst) {mc[u,v] <- mc[u,v] + cnt }
+       else{ md[u,v] <- dst; mc[u,v] <- cnt }
+     }
+   }}}
+   return(list(dis=md,cnt=mc))
+ }
> 
> vec.Closeness <- function(dis)
+ { n <- nrow(dis); return((n-1)/rowSums(dis)) }
> 
> vec.Betweenness <- function(dis,cnt){
+   n <- nrow(dis); bw <- rep(0,n)
+   for (v in 1:n) {
+     b <- 0
+     for(u in 1:n) { for(w in 1:n) {
+       if((cnt[u,w] > 0) && (u != w) && (u != v) && (v != w) &&
+             ((dis[u,v] + dis[v,w]) == dis[u,w]))
+         {b <- b + cnt[u,v]*cnt[v,w] / cnt[u,w]}
+     }}
+     bw[v] <- b/((n-1)*(n-2))
+   }
+   return(bw)
+ }
> 
> vec.betweenness <- function(m)
+ { mt <- mat.geodesics(m); return(vec.Betweenness(mt$dis,m$tcnt)) }
```

```
> library(sna)
> setwd("C:/Users/batagelj/Documents/papers/2018/moskva/NetR/nets")
> source("https://github.com/bavla/semirings/raw/master/betweenness.R")
> g <- c(
+  0, 1, 1, 0, 0, 0, 0, 0,
+  0, 0, 0, 1, 1, 0, 0, 0,
+  0, 0, 0, 1, 1, 0, 0, 0,
+  0, 0, 0, 0, 0, 1, 1, 0,
+  0, 0, 0, 0, 1, 1, 1, 0,
+  0, 0, 0, 0, 0, 0, 0, 1,
+  0, 0, 0, 0, 0, 0, 0, 1,
+  0, 1, 0, 0, 0, 0, 0, 0  )
> n <- 1:8
> g <- matrix(nrow=8,byrow=T,dimnames=list(n,n),data=g)
> g
  1 2 3 4 5 6 7 8
1 0 1 1 0 0 0 0 0
2 0 0 0 1 1 0 0 0
3 0 0 0 1 1 0 0 0
4 0 0 0 0 0 1 1 0
5 0 0 0 0 1 1 1 0
6 0 0 0 0 0 0 0 1
7 0 0 0 0 0 0 0 1
8 0 1 0 0 0 0 0 0
> gplot(g,displaylabels=TRUE)
>
> G <- mat.geodesics(g)
> G
$dis
    1 2   3 4 5 6 7 8
1 Inf 1   1 2 2 3 3 4
2 Inf 4 Inf 1 1 2 2 3
3 Inf 4 Inf 1 1 2 2 3
4 Inf 3 Inf 4 4 1 1 2
5 Inf 3 Inf 4 1 1 1 2
6 Inf 2 Inf 3 3 4 4 1
7 Inf 2 Inf 3 3 4 4 1
8 Inf 1 Inf 2 2 3 3 4

$cnt
  1 2 3 4 5 6 7 8
1 0 1 1 2 2 4 4 8
2 0 4 0 1 1 2 2 4
3 0 4 0 1 1 2 2 4
4 0 2 0 2 2 1 1 2
5 0 2 0 2 1 1 1 2
6 0 1 0 1 1 2 2 1
7 0 1 0 1 1 2 2 1
8 0 1 0 1 1 2 2 4

> b <- vec.betweenness(g)
> b
[1] 0.00000000 0.34523810 0.05952381 0.16666667 0.16666667 0.11904762 0.11904762 0.30952381
```

## Pitts Russian Rivers

Unpack in ''/rivers'' the data from [[pajek:data:conv:pitts|Data set]].
```
> library(sna)
> setwd("C:/Users/batagelj/Documents/papers/2018/moskva/NetR/nets/rivers")
> source("https://github.com/bavla/semirings/raw/master/betweenness.R")
> pitts <- read.paj("PittsGeo.net")
> list.vertex.attributes(pitts)
[1] "na"           "vertex.names" "x"            "y"            "z"           
> R <- as.matrix(pitts); nam <- rownames(R)
> x <- pitts %v% "x"; y <- pitts %v% "y"
```

{{ru:hse:snet:pics:brjansk.jpg?200}} {{ru:hse:snet:pics:dorogobuz.jpg?200}} {{ru:hse:snet:pics:elec.jpg?200}} {{ru:hse:snet:pics:murom.jpg?200}} 

{{ru:hse:snet:pics:smolensk.jpg?200}} {{ru:hse:snet:pics:velikij_novgorod.jpg?200}} {{ru:hse:snet:pics:vjazma.jpg?200}} {{ru:hse:snet:pics:vladimir.jpg?200}}

### Binary

```
> Geo <- mat.geodesics(R)
> Gcl <- vec.Closeness(Geo$dis); names(Gcl) <- nam
> Gbw <- vec.Betweenness(Geo$dis,Geo$cnt); names(Gbw) <- nam
> q <- order(Gbw,decreasing=TRUE)
> cbind(bw=Gbw[q],cl=Gcl[q])
                                 bw        cl
 Kolomna                0.351141367 0.3584906
 Moskva                 0.341112917 0.3486239
 Ksnyatin               0.296721534 0.3166667
 Kozelsk                0.214065569 0.3220339
 Dorogobuzh             0.189087584 0.3089431
 Tver                   0.159225767 0.2923077
 Vladimir               0.157159791 0.2814815
 Vyazma                 0.141979950 0.3304348
 Bryansk                0.116639572 0.2878788
 Dedoslavl              0.112091038 0.2857143
 Mtsensk                0.102492718 0.2657343
 Smolensk               0.101303936 0.2533333
 B                      0.087743006 0.2533333
 Nizhniy_Novgorod       0.080844002 0.2331288
 A                      0.074114340 0.2835821
 Mozhaysk               0.064383933 0.3166667
 Pereslavl              0.056187767 0.2794118
 Uglich                 0.046941679 0.2467532
 Vyshny_Volochyok       0.043859649 0.2375000
 Isady-Ryazan           0.041488857 0.2405063
 Kursk                  0.028473210 0.2345679
 Vitebsk                0.027738265 0.2183908
 Volokolamsk            0.024425930 0.2714286
 Suzdal                 0.024276908 0.2303030
 Chernigov              0.022569938 0.2076503
 Novgorod               0.021811285 0.2134831
 Novgorod_Severskiy     0.021052632 0.2317073
 Kiev                   0.020697013 0.2171429
 Dubok                  0.020104315 0.2389937
 Murom                  0.019440493 0.2196532
 Pronsk                 0.016121385 0.2360248
 C                      0.015291607 0.2159091
 Karachev               0.014485538 0.2389937
 Yelets                 0.012185870 0.2209302
 Tula                   0.009554291 0.2585034
 Rostov                 0.004267425 0.2065217
 Dmitrov                0.004267425 0.2483660
 Yaroslavl              0.001422475 0.2043011
 Bolgar                 0.000000000 0.1900000
> plot(Gbw,Gcl,pch=16,col="red",xlim=c(-0.02,0.35),main="Pitts - binary")
> text(Gbw,Gcl,nam,cex=0.5)
```
### Distance

```
> mat.dist <- function(R,x,y){
+   n <- nrow(R); D <- R
+   for(i in 1:n) { for(j in 1:n) {
+     if(R[i,j]>0) D[i,j] <- sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2) }}
+   return(D)
+ }
> 
> Deo <- mat.geodesics(mat.dist(R,x,y))
> Dcl <- vec.Closeness(Deo$dis); names(Dcl) <- nam
> Dbw <- vec.Betweenness(Deo$dis,Deo$cnt); names(Dbw) <- nam
> p <- order(Dbw,decreasing=TRUE)
> cbind(bw=Dbw[p],cl=Dcl[p])
                                 bw         cl
 Moskva                 0.238975818 0.21250622
 Ksnyatin               0.217638691 0.19895479
 Kolomna                0.207681366 0.20971455
 Dorogobuzh             0.170697013 0.18185725
 Kozelsk                0.159317212 0.19084617
 Mozhaysk               0.145092461 0.20713976
 Dedoslavl              0.133712660 0.19798117
 Tver                   0.132290185 0.19174094
 Vyazma                 0.129445235 0.19372181
 Bryansk                0.109530583 0.17079453
 B                      0.106685633 0.17333355
 Volokolamsk            0.106685633 0.19863966
 Tula                   0.088193457 0.19367514
 Vladimir               0.086770982 0.15823979
 A                      0.071123755 0.17294217
 Isady-Ryazan           0.065433855 0.16560027
 Smolensk               0.064011380 0.15353368
 Nizhniy_Novgorod       0.059743954 0.11880757
 C                      0.054054054 0.15657771
 Novgorod_Severskiy     0.049786629 0.13593519
 Suzdal                 0.048364154 0.15546934
 Pronsk                 0.048364154 0.16105212
 Mtsensk                0.045519203 0.16559230
 Pereslavl              0.045519203 0.18198257
 Murom                  0.045519203 0.13205730
 Uglich                 0.042674253 0.17449043
 Vyshny_Volochyok       0.027027027 0.15057193
 Dubok                  0.024182077 0.16612955
 Dmitrov                0.022759602 0.18560002
 Chernigov              0.021337127 0.11004909
 Karachev               0.018492176 0.15798282
 Vitebsk                0.017069701 0.12114010
 Rostov                 0.008534851 0.15590427
 Kiev                   0.004267425 0.10347625
 Yelets                 0.004267425 0.15427963
 Novgorod               0.001422475 0.10667918
 Yaroslavl              0.001422475 0.14758213
 Kursk                  0.000000000 0.13815674
 Bolgar                 0.000000000 0.07328784
> plot(Dbw,Dcl,pch=16,col="red",xlim=c(-0.02,0.24),main="Pitts - distance")
> text(Dbw,Dcl,nam,cex=0.5)
```

### Labels in Cyrillic

```
> A <- read.table("pitts.nam",header=FALSE,skip=2)
> namCy <- as.vector(A$V2)
> Encoding(namCy) <- "UTF-8"
> plot(Dbw,Dcl,pch=16,col="red",xlim=c(-0.02,0.24),main="Pitts - distance")
> text(Dbw,Dcl,namCy,cex=0.5)
```

### Closeness and betweenness in Sna

Sna does not take weights into account.
```
> Dis <- as.network(D)
> library(sna)
> (Rc <- closeness(pitts))
 [1] 0.2159091 0.2209302 0.2567568 0.2196532 0.2099448 0.2345679 0.2375000 0.2923077 0.2420382 0.3275862 0.3140496 0.3362832
[13] 0.2878788 0.2968750 0.2405063 0.3220339 0.2500000 0.2065217 0.2087912 0.2567568 0.2183908 0.2331288 0.2857143 0.2360248
[25] 0.1919192 0.2435897 0.2389937 0.2420382 0.2235294 0.2695035 0.2620690 0.2900763 0.2835821 0.3653846 0.3551402 0.3220339
[37] 0.2516556 0.2753623 0.2222222
> (Dc <- closeness(Dis))
 [1] 0.2159091 0.2209302 0.2567568 0.2196532 0.2099448 0.2345679 0.2375000 0.2923077 0.2420382 0.3275862 0.3140496 0.3362832
[13] 0.2878788 0.2968750 0.2405063 0.3220339 0.2500000 0.2065217 0.2087912 0.2567568 0.2183908 0.2331288 0.2857143 0.2360248
[25] 0.1919192 0.2435897 0.2389937 0.2420382 0.2235294 0.2695035 0.2620690 0.2900763 0.2835821 0.3653846 0.3551402 0.3220339
[37] 0.2516556 0.2753623 0.2222222
> (Rb <- betweenness(pitts))
 [1]  30.66667  39.00000 142.43333  29.10000  31.73333  29.60000  40.03333 163.99524  20.36667 300.97619 265.85714 199.62381
[13] 104.20476 223.87143  61.66667 417.19048  66.00000   2.00000   6.00000 123.36667  21.50000  34.13333 220.96667 113.66667
[25]   0.00000  58.33333  22.66667  28.26667  17.13333 144.10476  13.43333 157.60000  79.00000 493.70476 479.60476  90.52381
[37]   6.00000  34.34286  27.33333
> (Db <- betweenness(Dis))
 [1]  30.66667  39.00000 142.43333  29.10000  31.73333  29.60000  40.03333 163.99524  20.36667 300.97619 265.85714 199.62381
[13] 104.20476 223.87143  61.66667 417.19048  66.00000   2.00000   6.00000 123.36667  21.50000  34.13333 220.96667 113.66667
[25]   0.00000  58.33333  22.66667  28.26667  17.13333 144.10476  13.43333 157.60000  79.00000 493.70476 479.60476  90.52381
[37]   6.00000  34.34286  27.33333
```

## Transitions

```
> P <- rbind(
+   c(0.0, 0.8, 0.0, 0.0, 0.2),
+   c(0.0, 0.0, 0.5, 0.0, 0.5),
+   c(0.0, 0.0, 0.3, 0.7, 0.0),
+   c(0.0, 0.0, 0.0, 0.0, 1.0),
+   c(1.0, 0.0, 0.0, 0.0, 0.0) )
> rownames(P) <- colnames(P) <- paste("S",1:5,sep="")
> P
   S1  S2  S3  S4  S5
S1  0 0.8 0.0 0.0 0.2
S2  0 0.0 0.5 0.0 0.5
S3  0 0.0 0.3 0.7 0.0
S4  0 0.0 0.0 0.0 1.0
S5  1 0.0 0.0 0.0 0.0
> library(sna)
> plot.sociomatrix(P)
> gplot(P,displaylabels=TRUE,diag=TRUE,edge.curve=0.1,loop.cex=4,edge.lwd=10)
> x <- c(10000,0,0,0,0)
> y <- x; k <- 0
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
0 : 10000 0 0 0 0 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
1 : 0 8000 0 0 2000 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
2 : 2000 0 4000 0 4000 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
3 : 4000 1600 1200 2800 400 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
4 : 400 3200 1160 840 4400 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
5 : 4400 320 1948 812 2520 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
6 : 2520 3520 744.4 1363.6 1852 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
7 : 1852 2016 1983.32 521.08 3627.6 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
8 : 3627.6 1481.6 1602.996 1388.324 1899.48 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
9 : 1899.48 2902.08 1221.699 1122.097 2854.644 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
10 : 2854.644 1519.584 1817.55 855.1892 2953.033 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
11 : 2953.033 2283.715 1305.057 1272.285 2185.91 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
12 : 2185.91 2362.427 1533.375 913.5398 3004.749 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
13 : 3004.749 1748.728 1641.226 1073.362 2531.935 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
14 : 2531.935 2403.799 1366.732 1148.858 2548.676 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
15 : 2548.676 2025.548 1611.919 956.7122 2857.145 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
16 : 2857.145 2038.941 1496.35 1128.343 2479.221 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
17 : 2479.221 2285.716 1468.375 1047.445 2719.243 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
18 : 2719.243 1983.377 1583.37 1027.863 2686.147 
...
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
51 : 2651.384 2121.438 1515.048 1060.649 2651.481 
> cat(k,":",y,"\n"); k <- k+1; y <- y%*%P
52 : 2651.481 2121.107 1515.233 1060.534 2651.645 
>
> P1 <- cbind(P,0)
> P1 <- rbind(P1,0)
> rownames(P1)[6] <- colnames(P1)[6] <- "S6"
> P1[4,5] <- 0.6; P1[4,6] <- 0.4
> P1
   S1  S2  S3  S4  S5  S6
S1  0 0.8 0.0 0.0 0.2 0.0
S2  0 0.0 0.5 0.0 0.5 0.0
S3  0 0.0 0.3 0.7 0.0 0.0
S4  0 0.0 0.0 0.0 0.6 0.4
S5  1 0.0 0.0 0.0 0.0 0.0
S6  0 0.0 0.0 0.0 0.0 0.0
> plot.sociomatrix(P1)
> gplot(P1,displaylabels=TRUE,diag=TRUE,edge.curve=0.1,loop.cex=4,edge.lwd=10)
```

<img src="https://github.com/bavla/netstat/blob/master/HSE/24/pics/p1.png" width=400>

Rows with sum different from 1; zero sum. 

Kelvin Lancaster: Mathematical Economics.
[[https://www.amazon.com/Mathematical-Economics-Dover-Computer-Science/dp/0486653919|Dover]], 2011

```
> P1[5,1] <- 1.18
> P1
> plot.sociomatrix(P1)
> gplot(P1,displaylabels=TRUE,diag=TRUE,loop.cex=4,edge.lwd=10)
> x <- c(10000,0,0,0,0,0)
> y <- x; k <- 0
> cat(k,":",y,"=",sum(y),"\n"); k <- k+1; y <- y%*%P1
```

<img src="https://github.com/bavla/netstat/blob/master/HSE/24/pics/q.png" width=400>

```
> Q <- rbind(
+   c(0.0, 0.0, 0.2, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0), 
+   c(0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 
+   c(0.0, 0.0, 0.0, 0.0, 0.3, 0.7, 0.0, 0.0, 0.0), 
+   c(0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.8, 0.0, 0.0), 
+   c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0), 
+   c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0), 
+   c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0), 
+   c(0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), 
+   c(0.6, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) )
> 
> rownames(Q) <- colnames(Q) <- paste("Q",1:9,sep="")
> write.dl(Q,"./test/Q.net",vertex.lab=rownames(Q),matrix.lab="Periodic transition matrix")
> plot.sociomatrix(Q)
> gplot(Q,displaylabels=TRUE,edge.lwd=10)
> x <- c(10000,0,0,0,0,0,0,0,0)
> y <- x; k <- 0
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
0 10000 0 0 0 0 0 0 0 0 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
1 0 0 2000 8000 0 0 0 0 0 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
2 0 0 0 0 600 3000 6400 0 0 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
3 0 0 0 0 0 0 0 7000 3000 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
4 5300 4700 0 0 0 0 0 0 0 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
5 0 0 5760 4240 0 0 0 0 0 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
6 0 0 0 0 1728 4880 3392 0 0 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
7 0 0 0 0 0 0 0 5120 4880 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
8 5488 4512 0 0 0 0 0 0 0 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
9 0 0 5609.6 4390.4 0 0 0 0 0 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
10 0 0 0 0 1682.88 4804.8 3512.32 0 0 
...
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
19 0 0 0 0 0 0 0 5192.312 4807.688 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
20 5480.769 4519.231 0 0 0 0 0 0 0 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
21 0 0 5615.385 4384.615 0 0 0 0 0 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
22 0 0 0 0 1684.615 4807.692 3507.692 0 0 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
23 0 0 0 0 0 0 0 5192.308 4807.692 
> cat(k,y,"\n"); k <- k+1; y <- y%*%Q
24 5480.769 4519.231 0 0 0 0 0 0 0 
```

## Markov chains

### Data 

Pics: [[http://vladowiki.fmf.uni-lj.si/lib/exe/fetch.php?media=ru:hse:snet:pics:roberts1.jpg|Roberts1]]; [[http://vladowiki.fmf.uni-lj.si/lib/exe/fetch.php?media=ru:hse:snet:pics:roberts2.jpg|Roberts2]]
 
```
> setwd("C:/Users/batagelj/Documents/papers/2018/moskva/NetR/nets")
> library(sna)
> # pasture
> # F.S. Roberts: Discrete mathematical models, p.266
> P <- c(  3/5 , 3/10,  0 , 1/10,
+          1/10, 2/5 , 1/2,  0  ,
+          3/4 ,  0  , 1/5, 1/20,
+           0  ,  0  ,  0 ,  1   )
> n <- c('soil','grass','cattle','outside')
> P <- matrix(nrow=4,byrow=T,dimnames=list(n,n),data=P)
> 
> # summer weather
> # F.S. Roberts: Discrete mathematical models, p.267
> S <- c( 1/3, 1/2, 1/6,
+         1/2, 1/3, 1/6,
+         1/3, 1/3, 1/3 )
> n <- c('hot','moderate','cool')
> S <- matrix(nrow=3,byrow=T,dimnames=list(n,n),data=S)
> 
> # gambler's ruin
> # F.S. Roberts: Discrete mathematical models, p.268
> p <- 1/3
> G <- c( 1 , 0 , 0 , 0, 0,
+        1-p, 0 , p , 0, 0,
+         0 ,1-p, 0 , p, 0,
+         0 , 0 ,1-p, 0, p,
+         0 , 0 , 0 , 0, 1 )
> n <- c('$0','$1','$2','$3','$4')
> G <- matrix(nrow=5,byrow=T,dimnames=list(n,n),data=G)
> 
> P
        soil grass cattle outside
soil    0.60   0.3    0.0    0.10
grass   0.10   0.4    0.5    0.00
cattle  0.75   0.0    0.2    0.05
outside 0.00   0.0    0.0    1.00
> S
               hot  moderate      cool
hot      0.3333333 0.5000000 0.1666667
moderate 0.5000000 0.3333333 0.1666667
cool     0.3333333 0.3333333 0.3333333
> G
          $0        $1        $2        $3        $4
$0 1.0000000 0.0000000 0.0000000 0.0000000 0.0000000
$1 0.6666667 0.0000000 0.3333333 0.0000000 0.0000000
$2 0.0000000 0.6666667 0.0000000 0.3333333 0.0000000
$3 0.0000000 0.0000000 0.6666667 0.0000000 0.3333333
$4 0.0000000 0.0000000 0.0000000 0.0000000 1.0000000
```

[[..:snet:dat:mc|Additional data]]

### Computing



```
> setwd("C:/Users/batagelj/Documents/papers/2018/moskva/NetR/nets")
> library(sna)
> p <- c(0,0,1)
> q <- p
> q
[1] 0 0 1
> (q <- q %*% S) 
           hot  moderate      cool
[1,] 0.3333333 0.3333333 0.3333333
> (q <- q %*% S) 
           hot  moderate      cool
[1,] 0.3888889 0.3888889 0.2222222
> (q <- q %*% S) 
           hot  moderate      cool
[1,] 0.3981481 0.3981481 0.2037037
> (q <- q %*% S) 
           hot  moderate      cool
[1,] 0.3996914 0.3996914 0.2006173
> (q <- q %*% S) 
           hot  moderate      cool
[1,] 0.3999486 0.3999486 0.2001029
> (q <- q %*% S) 
           hot  moderate      cool
[1,] 0.3999914 0.3999914 0.2000171
> (q <- q %*% S) 
           hot  moderate      cool
[1,] 0.3999986 0.3999986 0.2000029
> (q <- q %*% S) 
           hot  moderate      cool
[1,] 0.3999998 0.3999998 0.2000005
> (q <- q %*% S) 
     hot moderate      cool
[1,] 0.4      0.4 0.2000001
> (q <- q %*% S) 
     hot moderate cool
[1,] 0.4      0.4  0.2
> (q <- q %*% S) 
     hot moderate cool
[1,] 0.4      0.4  0.2
```

```
# fixed point of stable chain
n <- nrow(S)
A <- t(S)-diag(n); A[n,] <- rep(1,n)
b <- rep(0,n); b[n] <- 1
w <- solve(A,b)
W <- matrix(nrow=n,ncol=n,byrow=T,w,dimnames=dimnames(S))
Z <- solve(diag(n) - (S - W))
E <- (diag(n) - Z + matrix(nrow=n,ncol=n,byrow=T,diag(Z))) %*% diag(1/w)
dimnames(E) <- dimnames(S)
> A
                hot   moderate      cool
hot      -0.6666667  0.5000000 0.3333333
moderate  0.5000000 -0.6666667 0.3333333
cool      1.0000000  1.0000000 1.0000000
> w
     hot moderate     cool 
     0.4      0.4      0.2 
> W
         hot moderate cool
hot      0.4      0.4  0.2
moderate 0.4      0.4  0.2
cool     0.4      0.4  0.2
> Z
                 hot    moderate  cool
hot       0.94857143  0.09142857 -0.04
moderate  0.09142857  0.94857143 -0.04
cool     -0.08000000 -0.08000000  1.16
> E
              hot moderate cool
hot      2.500000 2.142857    6
moderate 2.142857 2.500000    6
cool     2.571429 2.571429    5


> max(abs(S %*% W - W %*% S))
[1] 8.326673e-17
> max(abs(S %*% W - W))
[1] 5.551115e-17
```

```
source("https://github.com/bavla/semirings/raw/master/semirings.R")
gplot(G,displaylabels=TRUE,diag=TRUE,loop.cex=3)
n <- 5; m <- 2
G <- Sr.Perm(G,c(1,5,2,3,4))
Q <- G[(m+1):n,(m+1):n]
N <- solve(diag(n-m)-Q)
R <- G[(m+1):n,1:m]
B <- N %*% R

> G
          $0        $4        $1        $2        $3
$0 1.0000000 0.0000000 0.0000000 0.0000000 0.0000000
$4 0.0000000 1.0000000 0.0000000 0.0000000 0.0000000
$1 0.6666667 0.0000000 0.0000000 0.3333333 0.0000000
$2 0.0000000 0.0000000 0.6666667 0.0000000 0.3333333
$3 0.0000000 0.3333333 0.0000000 0.6666667 0.0000000
> Q
          $1        $2        $3
$1 0.0000000 0.3333333 0.0000000
$2 0.6666667 0.0000000 0.3333333
$3 0.0000000 0.6666667 0.0000000
> N
    $1  $2  $3
$1 1.4 0.6 0.2
$2 1.2 1.8 0.6
$3 0.8 1.2 1.4
> R
          $0        $4
$1 0.6666667 0.0000000
$2 0.0000000 0.0000000
$3 0.0000000 0.3333333
> B
          $0         $4
$1 0.9333333 0.06666667
$2 0.8000000 0.20000000
$3 0.5333333 0.46666667
```

#### Frog

<img src="https://github.com/bavla/netstat/blob/master/HSE/24/pics/frog.png" width=400>

```
> library(sna)
> MC2 <- c(
+  0 ,  1 ,  0 ,  0 ,  0 ,  0 ,  0 ,
+  0 , 1/2, 1/2,  0 ,  0 ,  0 ,  0 ,
+ 1/2,  0 , 1/2,  0 ,  0 ,  0 ,  0 ,
+  0 ,  0 , 1/4, 1/2, 1/4,  0 ,  0 ,
+  0 ,  0 ,  0 ,  0 ,  0 , 1/2, 1/2,
+  0 ,  0 ,  0 ,  1 ,  0 ,  0 ,  0 ,
+  0 ,  0 ,  0 ,  0 ,  0 ,  0 ,  1   )
> frog <- matrix(MC2,nrow=7,byrow=TRUE)
> rownames(frog)<- colnames(frog) <- paste("L",1:7,sep="")
> frog
    L1  L2   L3  L4   L5  L6  L7
L1 0.0 1.0 0.00 0.0 0.00 0.0 0.0
L2 0.0 0.5 0.50 0.0 0.00 0.0 0.0
L3 0.5 0.0 0.50 0.0 0.00 0.0 0.0
L4 0.0 0.0 0.25 0.5 0.25 0.0 0.0
L5 0.0 0.0 0.00 0.0 0.00 0.5 0.5
L6 0.0 0.0 0.00 1.0 0.00 0.0 0.0
L7 0.0 0.0 0.00 0.0 0.00 0.0 1.0
> gplot(frog,displaylabels=TRUE,diag=TRUE,loop.cex=3)
> mat.shrink <- function(A,C,lc){
+   n <- nrow(A); other <- setdiff(1:n,C); m <- length(other)
+   R <- matrix(0,nrow=m+1,ncol=m+1)
+   rownames(R) <- colnames(R) <- c(lc,rownames(A)[other])
+   R[2:(m+1),2:(m+1)] <- A[other,other]; R[1,1] <- 1
+   R[2:(m+1),1] <- rowSums(A[other,C])
+   return(R)
+ }
> (FR <- mat.shrink(frog,c(1,2,3),"A1"))
     A1  L4   L5  L6  L7
A1 1.00 0.0 0.00 0.0 0.0
L4 0.25 0.5 0.25 0.0 0.0
L5 0.00 0.0 0.00 0.5 0.5
L6 0.00 1.0 0.00 0.0 0.0
L7 0.00 0.0 0.00 0.0 1.0
> gplot(FR,displaylabels=TRUE,diag=TRUE,loop.cex=3)
> source("https://github.com/bavla/semirings/raw/master/semirings.R")
> G <- Sr.Perm(FR,c(1,5,2,3,4))
> G
      A  L7  L4   L5  L6
A  1.00 0.0 0.0 0.00 0.0
L7 0.00 1.0 0.0 0.00 0.0
L4 0.25 0.0 0.5 0.25 0.0
L5 0.00 0.5 0.0 0.00 0.5
L6 0.00 0.0 1.0 0.00 0.0
> n <- 5; m <- 2
> Q <- G[(m+1):n,(m+1):n]
> N <- solve(diag(n-m)-Q)
> R <- G[(m+1):n,1:m]
> B <- N %*% R
> N
         L4        L5        L6
L4 2.666667 0.6666667 0.3333333
L5 1.333333 1.3333333 0.6666667
L6 2.666667 0.6666667 1.3333333
> R
      A  L7
L4 0.25 0.0
L5 0.00 0.5
L6 0.00 0.0
> B
           A        L7
L4 0.6666667 0.3333333
L5 0.3333333 0.6666667
L6 0.6666667 0.3333333 
```

R package [[https://cran.r-project.org/web/packages/markovchain/|markovchain]].


### To do

  * Add row/col names to results 
  * Implement shrinking of absorbing states in MC (done June 2, 2018)
  * Implement matrix/vector operations based on semiring. Reprogram [[mat|Mat]] library. {{pub:pdf:relcalc.pdf|RelCalc}}
  * Transition matrix from a given text determined by consecutive letters pairs. Generate random text. 




<hr/>

[HSE](../2024.md); [Docs](doc.md)






