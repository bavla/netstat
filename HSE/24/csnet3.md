# Basic models

Results of fire experiments:
```
> x <- c(16,24,34,41,46,51,56,59,61,65,75)
> y <- c(0.5,0.6,0.6,1,1.5,3,12.3,22.9,75.3,94.6,99.3)
> plot(x,y,pch=16,type="b")

> density <- c(10,25,35,45,50,55,57,59,60,61,62,65,75,90)
> burned <- c(0.6,0.6,0.7,1.3,2.8,3.5,22.5,32.1,71.9,83.6,89.3,94.6,99.3,100)
> plot(density,burned,type="b",col="red",lwd=2)
```

## Random graphs in sna

  * rgnm - Erdös-Rényi: Draw Density-Conditioned Random Graphs
  * rgraph - Gilbert: Generate Bernoulli Random Graphs
  * rguman - Draw Dyad Census-Conditioned Random Graphs
  * rgnmix - Mixing-conditioned random graphs
  * rgws - The Watts-Strogatz rewiring model

Melissa Clarkson: [[http://www.melissaclarkson.com/resources/R_guides/documents/random_graphs_Ver1.pdf|Random graphs in sna]]

### "Home made"

```
> setwd("C:/Users/batagelj/Documents/papers/2018/moskva/NetR/nets")
> library(sna)
> s <- rbinom(
+   size = 1,  # 1 binomial trial -- a Bernoulli trial
+   n = 50 ^ 2,  # number of 0's or 1's to generate
+   p = .05 )
> g <- network(matrix(s, nrow = 50, ncol = 50))
> plot(g)
>
> n <- 50; p <- 0.05
> s <- numeric(n*n); i <- 0
> repeat{
+   i <- i+1+rgeom(1,p)
+   if(i > n*n) break
+   s[i] <- 1
+ }
> g <- network(matrix(s,nrow=n,ncol=n))
> plot(g)
```
Problems: loops, for large networks use list of links.

### sna <-> igraph 

```
> ig <- igraph::erdos.renyi.game(30,0.1)
> ig
IGRAPH 35441a0 U--- 30 42 -- Erdos renyi (gnp) graph
+ attr: name (g/c), type (g/c), loops (g/l), p (g/n)
+ edges from 35441a0:
 [1]  1-- 6  1-- 8  7-- 8  5-- 9  8--12  1--14 11--14  2--15  5--16 10--16
[11]  2--17  4--20 13--20  4--21 15--21 17--21 11--22 16--22 21--22  3--23
[21]  7--23 10--23 14--23 22--23 14--24 16--24 22--24  4--25 14--25 22--25
[31]  5--26 12--26 23--26 11--27  5--28  5--29 15--29 22--29 27--29  2--30
[41] 21--30 24--30
> ga <- intergraph::asNetwork(ig) 
> gplot(ga)
> gi <- intergraph::asIgraph(g) 
> gi
IGRAPH 4616d77 D--- 50 113 -- 
+ attr: na (v/l), vertex.names (v/n), na (e/l)
+ edges from 4616d77:
 [1] 20-> 1 29-> 1 41-> 1 50-> 1 37-> 2 44-> 3  6-> 4 20-> 4 28-> 4 39-> 4
[11] 47-> 4 29-> 5 38-> 5 16-> 6  5-> 7 25-> 7  3-> 8  2->10 49->11 29->12
[21] 30->12 36->12 20->13 26->13 28->13 30->13  8->14  5->15 17->15 21->15
[31] 30->15 34->15 37->16  4->17 15->17  7->18 41->18 14->19 17->19 47->19
[41]  6->20  7->20 14->20 33->20 33->21 38->21 47->21 10->22 37->22  7->23
[51] 16->23  5->24  8->24 37->24  1->25 31->25 13->26 21->27 24->27  3->28
[61] 14->28 16->28 30->28 46->28 50->28 23->29 38->29 39->29  7->30 12->30
[71] 41->30 47->30  8->31 45->31 26->32 30->32 44->32 47->33 35->34 39->35
+ ... omitted several edges
```

### rgnm

Fixed number of links

```
> setwd("C:/Users/batagelj/Documents/papers/2018/moskva/NetR/nets")
> library(sna)
> help(rgnm)
starting httpd help server ... done
> Erdos.gph <- rgnm(1, 100, 200, mode="digraph")
> class(Erdos.gph)
[1] "matrix"
> Erdos.net <- network(Erdos.gph)
> Erdos.net
 Network attributes:
  vertices = 100 
  directed = TRUE 
  hyper = FALSE 
  loops = FALSE 
  multiple = FALSE 
  bipartite = FALSE 
  total edges= 200 
    missing edges= 0 
    non-missing edges= 200 

 Vertex attribute names: 
    vertex.names 

No edge attributes
> plot(Erdos.net)
```

#### Degree distribution

```
> g <- rgnm(1,1000,4000,mode="graph",return.as.edgelist=TRUE,diag=FALSE)
> d <- degree(g,gmode="graph")
> t <- table(d)
> p <- t/sum(t)
> m <- mean(d); s <- sd(d)
> m
[1] 8
> s
[1] 2.874415
> plot(p)
> x <- 1:18
> y <- dnorm(x,mean=m,sd=s)
> lines(x,y,col="red",lw=2)
> z <- dpois(x,m)
> lines(x,z,col="blue",lw=2)
```

#### Random acyclic network

```
> n <- 50; p <- 0.05
> m <- matrix(0,nrow=n,ncol=n)
> s <- rbinom(size=1,n=n*(n-1)/2,p=p)
> m[upper.tri(m, diag = FALSE)] <- s
> g <- network(m)
> plot.sociomatrix(g,cex=0.4,main="random acyclic network")
```

#### Random two-mode network

[[https://www.rdocumentation.org/packages/sna/versions/2.4/topics/gplot|gplot]]

```
> n1 <- 30; n2 <- 20; p <- 0.2; n <- n1+n2
> m <- matrix(0,nrow=n,ncol=n)
> s <- rbinom(size=1,n=n1*n2,p=p)
> m[1:n1,(n1+1):n] <- s 
> g <- network(m,directed=FALSE)
> plot.sociomatrix(g,cex=0.4,main="random two-mode network")
> b <- c(rep("red",n1),rep("blue",n2))
> gplot(g,gmode="graph",vertex.cex=1,vertex.col=b,mode="kamadakawai",
+ main="random two-mode network")
```


### rgraph

Bernoulli graphs - link selected with probability p
```
> help(rgraph)
> Erdos2.net <- network(rgraph(100,1,tprob=0.02,mode="graph"))
> gplot(Erdos2.net,gmode="graph")
> gden(Erdos2.net)
[1] 0.02232323
> 
> Erdos.gphs <- rgraph(100, m=500, tprob=0.02)  # 500 random graphs
> hist(gden(Erdos.gphs))
```

#### Giant component

```
> r <- 100; p <- 1:r/r/r*2; s <- numeric(r)
> for(i in 1:r) 
+   s[i] <- max(component.dist(rgraph(r,1,tprob=p[i]),connected="weak")$csize)
> plot(p,s,type="p",pch=16,cex=0.5,col="red",main="Giant component")
> lines(p,s)
>
> k <- 500; r <- 100; p <- 1:r/r/r*2; s <- m <- M <- numeric(r)
> for(i in 1:100) { t <- numeric(k)
+   for(j in 1:k) 
+     t[j] <- max(component.dist(rgraph(100,1,tprob=p[i]),connected="weak")$csize)
+   s[i] <- mean(t); m[i] <- min(t); M[i] <- max(t)
+ }
> plot(p,s,type="p",pch=16,cex=0.5,col="red",main="Giant component")
> lines(p,s); lines(p,m,col="blue"); lines(p,M,col="blue")
```

##### Size of the giant component 


```
> f <- 0.5; d <- 3
> for(i in 1:10) {cat(i,f,"\n"); f <- 1 - exp(-d*f)}
1 1 
2 0.9502129 
...
9 0.9404798 
10 0.9404798 
> d <- 1:50/10; s <- numeric(50); n <- integer(50)
> for(r in 1:50){ f <- 0.5; k <- 0
+   repeat{fo <- f; k <- k+1; f <- 1 - exp(-d[r]*f); if(abs(fo-f)<0.00001) break}
+   s[r] <- f; n[r] <- k 
+ }
> plot(d,s,type="b",pch=16,cex=0.7)
> plot(d,n)
> plot(d,n,ylim=c(0,50))
```

##### Eigenvalues

```
> gER <- rgraph(500,1,tprob=0.06,mode="graph")
> # gplot(gER,gmode="graph")
> la <- eigen(gER,only.values=TRUE)
> hist(la$values,breaks=15)
> summary(la$values)
```

### rguman

  * m - mutual dyads
  * a - asymmetric dyads
  * n - null dyads
```
> help(rguman)
> rman <- rguman(1, 20, mut=50, asym = 100, null = 40, method="exact")
> dyad.census(rman)
     Mut Asym Null
[1,]  50  100   40
> rman <- rguman(1, 20, mut=10, asym = 20, null = 160, method="exact")
> plot(network(rman))
> rman <- rguman(1, 20, mut=50, asym = 100, null = 40, method="probability")
> dyad.census(rman)
     Mut Asym Null
[1,]  46   97   47
```

### configuration model

Chung-Lu method
```
> library(statnet)
> data(samplk)
> gplot(samplk3, gmode="digraph")
> id <- degree(samplk3,cmode="indegree")
> od <- degree(samplk3,cmode="outdegree")
> id
 [1] 4 6 3 4 6 2 5 2 4 0 2 6 2 2 2 1 2 3
> od
 [1] 3 3 4 3 3 3 3 3 3 4 3 3 3 3 3 3 3 3
> m <- sum(id)
> m
[1] 56
> n <- network.size(samplk3)
> C <- matrix(0,nrow=n,ncol=n)
> rownames(C) <- colnames(C) <- paste("v",1:n,sep="")
> for(u in 1:n) for(v in 1:n) if(u!=v) if(runif(1) < od[u]*id[v]/m) C[u,v] <- 1
> (ic <- degree(C,cmode="indegree"))
 [1] 2 6 4 7 7 2 4 1 5 0 1 6 1 2 4 0 2 6
> (oc <- degree(C,cmode="outdegree"))
 [1] 6 2 3 6 3 4 2 5 3 4 3 1 6 3 1 5 1 2
> sum(ic)
[1] 60
> gplot(C,displaylabels=TRUE) 
```

#### Molloy-Reed method

```
> b <- c(); for(i in 1:n) b <- c(b,rep(i,od[i]))
> e <- c(); for(i in 1:n) e <- c(e,rep(i,id[i])) 
> m <- length(b); M <- matrix(0,nrow=n,ncol=n)
> rownames(M) <- colnames(M) <- paste("v",1:n,sep="") 
> I <- matrix(nrow=m,ncol=2)
> I[,1] <- sample(b,m); I[,2] <- sample(e,m)
> M[I] <- 1; diag(M) <- 0
> sum(M)
[1] 53
> gplot(M,displaylabels=TRUE)
```

### rgnmix


```
> help(rgnmix)
> type <- c(rep(1,4),rep(2,3),rep(3,5))
> type
 [1] 1 1 1 1 2 2 2 3 3 3 3 3
> mix <- rbind(c(9,3,0),c(2,0,10),c(11,0,14))
> mix
     [,1] [,2] [,3]
[1,]    9    3    0
[2,]    2    0   10
[3,]   11    0   14
> mix.gph <- rgnmix(1,type,mix,method="exact")
> mix.gph
      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
 [1,]    0    0    1    1    0    0    1    0    0     0     0     0
 [2,]    1    0    1    1    0    0    0    0    0     0     0     0
 [3,]    1    1    0    1    1    0    0    0    0     0     0     0
 [4,]    0    1    0    0    0    0    1    0    0     0     0     0
 [5,]    0    0    0    0    0    0    0    0    1     1     0     1
 [6,]    0    1    0    0    0    0    0    0    1     1     1     1
 [7,]    1    0    0    0    0    0    0    1    1     0     1     0
 [8,]    0    1    1    0    0    0    0    0    1     1     1     0
 [9,]    1    1    1    0    0    0    0    1    0     0     0     1
[10,]    0    1    1    0    0    0    0    1    1     0     0     1
[11,]    1    0    0    1    0    0    0    0    1     1     0     0
[12,]    0    0    1    1    0    0    0    1    1     1     1     0
> plot.sociomatrix(mix.gph)
```


### rgws


The Watts-Strogatz rewiring model (Small world)
```
> # cases, #nodes, dim, deg, prob
> ws.net <- network(rgws(1,60,1,5,0),directed=FALSE)
> gplot(ws.net,gmode="graph")
> ws.net <- network(rgws(1,60,1,5,0.1),directed=FALSE)
> gplot(ws.net,gmode="graph")
> ws.net <- network(rgws(1,60,1,5,0.2),directed=FALSE)
> gplot(ws.net,gmode="graph")
```

The libraries sna / network / statnet do not contain functions for determining diameter, average path length and (average) clustering coefficient; neither do support the generation of scale-free networks.

```
> ws.net <- network(rgws(1,12,1,2,0.3))
> plot(ws.net)
> (D <- geodist(ws.net))
$counts
      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
 [1,]    1    1    1    2    1    5    3    2    1     1     1     6
 [2,]    1    1    1    1    2    2    1    1    3     5     5     2
 [3,]    1    1    1    1    1    3    2    1    1     1     1     4
 [4,]    2    1    1    1    1    2    1    1    2     4     4     2
 [5,]    1    2    1    1    1    1    1    3    1     1     1     2
 [6,]    5    2    3    2    1    1    1    1    3     1     1     1
 [7,]    3    1    2    1    1    1    1    1    2     1     1     1
 [8,]    2    1    1    1    3    1    1    1    1     2     2     1
 [9,]    1    3    1    2    1    3    2    1    1     1     1     1
[10,]    1    5    1    4    1    1    1    2    1     1     2     1
[11,]    1    5    1    4    1    1    1    2    1     2     1     1
[12,]    6    2    4    2    2    1    1    1    1     1     1     1

$gdist
      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
 [1,]    0    1    1    2    2    4    3    3    3     4     4     5
 [2,]    1    0    1    1    2    3    2    2    3     4     4     4
 [3,]    1    1    0    1    1    3    2    2    2     3     3     4
 [4,]    2    1    1    0    1    2    1    1    2     3     3     3
 [5,]    2    2    1    1    0    2    1    2    1     2     2     3
 [6,]    4    3    3    2    2    0    1    1    2     1     1     1
 [7,]    3    2    2    1    1    1    0    1    2     2     2     2
 [8,]    3    2    2    1    2    1    1    0    1     2     2     2
 [9,]    3    3    2    2    1    2    2    1    0     1     1     2
[10,]    4    4    3    3    2    1    2    2    1     0     2     1
[11,]    4    4    3    3    2    1    2    2    1     2     0     2
[12,]    5    4    4    3    3    1    2    2    2     1     2     0
> (L <- sum(D$gdist)/nrow(D$gdist)/(nrow(D$gdist)-1))
[1] 2.121212
> (diam <- max(D$gdist))
[1] 5
```
Function ''gtrans'' computes the global clustering coefficient.

sna has two functions to apply WS rewiring  to a given graph: ''rewire.ws'' and ''rewire.ud''.

```
> g<-matrix(0,20,20)
> g[1,] <- 1
> gplot(g,mode="circle")
> gr <- rewire.ws(g,0.5)
> gplot(gr,mode="circle")
> gr <- rewire.ws(g,1)
> gplot(gr,mode="circle")
```



## Random graphs in igraph 

[[http://igraph.org/c/doc/igraph-Generators.html|Graph generators]]; [[http://igraph.org/c/doc/igraph-Generators.html#idm470936067488|Random]]

  * erdos.renyi.game — Generates a random (Erdos-Renyi) graph.
  * degree.sequence.game — Generates a random graph with a given degree 
  * watts.strogatz.game — The Watts-Strogatz small-world model
  * rewire — Rewire the edges of a graph with constant probability
  * barabasi.game — Generates a graph based on the Barabási-Albert model
  * grg.game — Generating geometric random graphs
  * growing.random.game	— Growing random graph generation


and some others.


### erdos.renyi.game

Erdös-Rényi random graphs
https://www.rdocumentation.org/packages/igraph/versions/0.1.1/topics/erdos.renyi.game
```
> g <- erdos.renyi.game(1000, 1/1000)
> degree.distribution(g)
[1] 0.350 0.364 0.206 0.060 0.017 0.002 0.001
```

### degree.sequence.game

https://www.rdocumentation.org/packages/igraph/versions/0.1.1/topics/degree.sequence.game

```
> g2 <- degree.sequence.game(1:10, 10:1)
> degree(g2, mode="out")
 [1]  1  2  3  4  5  6  7  8  9 10
> degree(g2, mode="in")
 [1] 10  9  8  7  6  5  4  3  2  1
> plot(g2)
> G2 <- simplify(g2)
> plot(G2)
```





### watts.strogatz.game

Small world
```
# R script: small_world.R
library(igraph)
g <- watts.strogatz.game(1, 100, 5, 0.05)
plot(g,layout=layout.circle)
average.path.length(g)
transitivity(g, type="average")
```

[[https://stackoverflow.com/questions/48853610/average-clustering-coefficient-of-a-network-igraph|average clustering coefficient]]

```
> # Global clustering coefficient
> transitivity(g)
[1] 0.04639175
> # Average clustering coefficient
> transitivity(g, type = "average")
[1] 0.04577047
> # The same as above
> mean(transitivity(g, type = "local"), na.rm = TRUE)
[1] 0.04577047
```

```
> p <- 0.0001
> for(i in 1:100) {q <- p[i]*1.35; if(q>1) break; p <- c(p,q)}
> k <- length(p); nr <- 50
> t <- a <- numeric(k)
> for(i in 1:k){ tt <- aa <- numeric(nr)
+   for(r in 1:nr) {g <- watts.strogatz.game(1,300,7,p[i])
+     aa[r] <- average.path.length(g)
+     tt[r] <- transitivity(g,type="average")}
+   a[i] <- mean(aa); t[i] <- mean(tt) }
> Ma <- max(a); Mt <- max(t); ma <- min(a); mt <- min(t)
> x <- (a-ma)/(Ma-ma); y <- (t-mt)/(Mt-mt)
> plot(p,x,pch=16,log="x",col="blue",ylim=c(0,1))
> points(p,y,pch=16,cex=0.7,col="red")
```

### rewire

[[http://igraph.org/r/doc/rewire.html|rewire]]

```
> g <- graph.lattice( length=100, dim=1, nei=5 )
> average.path.length(g)
[1] 7.141414
> gr <- rewire(g,each_edge(p = .1, loops = FALSE))
> plot(gr,layout=layout_in_circle)
> average.path.length(gr)
[1] 2.493131
> g3 <- rewire(g,each_edge(p = .3, loops = FALSE))
> plot(g3,layout=layout_in_circle)
> average.path.length(g3)
[1] 2.288081
```

### barabasi.game

[[https://www.rdocumentation.org/packages/igraph/versions/0.1.1/topics/barabasi.game|barabasi.game]]

```
> n <- 10000
> g <- barabasi.game(n,m=2)
> plot(degree.distribution(g),log="xy")
> d <- degree(g)
> mean(d)
[1] 3.9994
> sd(d)
[1] 13.3639
> p <- mean(d)/(n-1)
> er <- erdos.renyi.game(n,p)
> der <- degree(er)
> mean(der)
[1] 3.9498
> sd(der)
[1] 1.978957
> lines(degree.distribution(er),col="red",lw=2)
```

```
> n <- 100000
> g <- igraph::barabasi.game(n,m=2)
> d <- igraph::degree(g)
> t <- table(d)
> x <- as.numeric(names(t))
> plot(x,t,pch=16,cex=0.6,xlab='deg',ylab='num',main='BA')
> plot(x,t,log="xy",pch=16,cex=0.6,xlab='deg',ylab='num',main="BA logs")
> len <- length(t)
> z <- t[len:1]; y <- cumsum(z); z <- y[len:1]
> plot(x,z,log="xy",pch=16,cex=0.6,xlab='deg',ylab='num',main="BA cum logs")
```

#### SN5 citation network input degrees

[[http://vlado.fmf.uni-lj.si/pub/networks/data/WoS/SN5.zip|SN5]] ("social network*" AND SO=(Social networks)) plus most frequently cited works plus around 100 SNA researchers (collected January 2008).
 
The function plfit(x) fits the distribution p(x) = (alpha-1)/x<sub>min</sub> (x/x<sub>min</sub>)<sup>-alpha</sup>. It returns estimates of alpha and x_min and D - the Kolmogorov-Smirnov  goodness-of-fit statistic.

Install libraries VGAM (zeta function) and R.matlab .
```
> source("https://sites.santafe.edu/~aaronc/powerlaws/plfit.r")
> infile <- "https://raw.githubusercontent.com/bavla/biblio/master/dat/indegCite.vec"
> d <- read.table(infile,skip=1,header=FALSE)
> t <- table(d)
> x <- as.numeric(names(t))
> plot(x,t,log="xy",pch=16,main="SN5 cite input degrees / logs")
> len <- length(t)
> z <- t[len:1]; y <- cumsum(z); z <- y[len:1]
> plot(x,z,log="xy",pch=16,cex=0.6,xlab='deg',ylab='num',main="BA cum logs")
> D <- as.vector(d$V1); D <- D[D>0]
> r <- plfit(D)
> names(r)
[1] "xmin"  "alpha" "D"    
> r$alpha
[1] 2.45
> r$xmin
[1] 3
> r$D
[1] 0.006133737
> plot(x,t,log="xy",pch=16,main="SN5 cite input degrees / logs")
> b <- (r$alpha-1)*r$xmin**(r$alpha-1)
> abline(b-r$xmin+1,-r$alpha,col="red",lw=2)
```
 

<hr/>

[HSE](../2024.md); [Docs](doc.md)





