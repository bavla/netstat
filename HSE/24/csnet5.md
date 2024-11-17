# ====== Statistics ======





===== Combinatorics =====

<code>
> combn(letters[1:4], 2)
     [,1] [,2] [,3] [,4] [,5] [,6]
[1,] "a"  "a"  "a"  "b"  "b"  "c" 
[2,] "b"  "c"  "d"  "c"  "d"  "d" 
> combn(letters[1:4], 3)
     [,1] [,2] [,3] [,4]
[1,] "a"  "a"  "a"  "b" 
[2,] "b"  "b"  "c"  "c" 
[3,] "c"  "d"  "d"  "d" 
> utils:::menuInstallPkgs()
package ‘combinat’ successfully unpacked and MD5 sums checked
> library(combinat)
> permn(letters[1:4])
[[1]]
[1] "a" "b" "c" "d"

[[2]]
[1] "a" "b" "d" "c"

[[3]]
[1] "a" "d" "b" "c"
</code>
https://www.topcoder.com/blog/generating-combinations/

https://www.geeksforgeeks.org/heaps-algorithm-for-generating-permutations/

https://www.topcoder.com/blog/generating-permutations/


<code>
> startSubsets <- function(n){.n <<- n; .s <<- -1}
> 
> nextSubset <- function(){
+   .s <<- .s+1
+   S <- as.logical(intToBits(.s)[1:(.n+1)])
+   if(S[.n+1]) {.s <<- -1; return(NULL)}
+   return(S[1:.n])
+ }

#.Machine$integer.max

> startSubsets(3)
> repeat {S <- nextSubset(); if(is.null(S)) break; cat("{",letters[1:3][S],'}\n')}
{  }
{ a }
{ b }
{ a b }
{ c }
{ a c }
{ b c }
{ a b c }
</code>

[[https://www4.uwsp.edu/math/nwodarz/Math209Files/209-0809F-L10-Section06_03-AlgorithmsForGeneratingPermutationsAndCombinations-Notes.pdf|Permutations]] in increasing lexicographic order.

<code>
> startPerm <- function(n){.n <<- n; .s <<- 1:n; .first <<- TRUE }
> startPerm <- function(n){.n <<- n; .s <<- 1:n; .first <<- TRUE }
> nextPerm <- function(){
+   if(.first){ .first <<- FALSE; return(.s) }
+   m <- .n-1
+   while(.s[m]>.s[m+1]) m <- m-1
+   k <- .n
+   while(.s[m]>.s[k]) k <- k-1 
+   .t <<- .s[m]; .s[m] <<- .s[k]; .s[k] <<- .t
+   p <- m+1; q <- .n
+   while(p<q) {.t <<- .s[p]; .s[p] <<- .s[q]; .s[q] <<- .t; p <- p+1; q <- q-1 }
+   return(.s)
+ }
> 
> startPerm(4)
> tryCatch({for(i in 1:30) cat(i,nextPerm(),"\n")}, error=function(err){})
1 1 2 3 4 
2 1 2 4 3 
3 1 3 2 4 
4 1 3 4 2 
5 1 4 2 3 
6 1 4 3 2 
7 2 1 3 4 
8 2 1 4 3 
9 2 3 1 4 
10 2 3 4 1 
11 2 4 1 3 
12 2 4 3 1 
13 3 1 2 4 
14 3 1 4 2 
15 3 2 1 4 
16 3 2 4 1 
17 3 4 1 2 
18 3 4 2 1 
19 4 1 2 3 
20 4 1 3 2 
21 4 2 1 3 
22 4 2 3 1 
23 4 3 1 2 
24 4 3 2 1 
NULL
> 

</code>

https://blogs.msdn.microsoft.com/gpalem/2013/03/28/make-vectorize-your-friend-in-r/

===== Autocorrelation =====

==== Simple example ====

<code>
> library(statnet)
> g <- c(
+   0, 1, 0,  1, 0, 0,  0, 0, 0,
+   0, 0, 1,  0, 1, 0,  0, 0, 0,
+   0, 0, 0,  0, 0, 1,  0, 0, 0,
+   0, 0, 0,  0, 1, 0,  1, 0, 0,
+   0, 0, 0,  0, 0, 1,  0, 1, 0,
+   0, 0, 0,  0, 0, 0,  0, 0, 1,
+   0, 0, 0,  0, 0, 0,  0, 1, 0,
+   0, 0, 0,  0, 0, 0,  0, 0, 1,
+   0, 0, 0,  0, 0, 0,  0, 0, 0  )
> G <- symmetrize(matrix(g,byrow=TRUE,nrow=9))
> rownames(G) <- colnames(G) <- c("A","B","C","D","E","F","G","H","I")
> plot.sociomatrix(G)
> gplot(G,gmode="graph",displaylabels=TRUE)
> a1 <- c( 1, 2, 3, 2, 3, 4, 3, 4, 5 )
> gplot(G,gmode="graph",displaylabels=TRUE,vertex.cex=a1)
> nacf(G,a1,type="geary",mode="graph")[2]
        1 
0.3333333 
> nacf(G,a1,type="moran",mode="graph")[2]
  1 
0.5 
> a2 <- c(3, 4, 3, 4, 3, 2, 1, 2, 5 )
> nacf(G,a2,type="geary",mode="graph")[2]
1 
1 
> nacf(G,a2,type="moran",mode="graph")[2]
    1 
-0.25 
> a3 <- c(4, 1, 4, 2, 5, 2, 3, 3, 3 )
> nacf(G,a3,type="geary",mode="graph")[2]
       1 
1.833333 
> nacf(G,a3,type="moran",mode="graph")[2]
     1 
-0.875 
</code>

==== Florentine families ====

<code>
> library(statnet)
> data(florentine)
> fw <- flomarriage %v% "wealth"
> gplot(flomarriage,displaylabels=TRUE,vertex.cex=0.025*fw,
+ gmode="graph",main='Florentine families')
> I <- nacf(flomarriage,fw,type="moran",mode="graph")
> I
          0           1           2           3           4           5 
 1.00000000 -0.31073529  0.06531299 -0.06045322 -0.06267282  0.01039729 
          6           7           8           9          10          11 
 0.00000000  0.00000000  0.00000000  0.00000000  0.00000000  0.00000000 
         12          13          14          15 
 0.00000000  0.00000000  0.00000000  0.00000000 
> C <- nacf(flomarriage,fw,type="geary",mode="graph")
> C
          0           1           2           3           4           5 
0.000000000 1.683607336 0.811642725 0.899673648 0.826831324 0.006496392 
          6           7           8           9          10          11 
0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 
         12          13          14          15 
0.000000000 0.000000000 0.000000000 0.000000000 
> I[2]
         1 
-0.3107353 
> C[2]
       1 
1.683607 
</code>





===== CUG =====

<code>
> library(statnet)
> data(florentine)
> rSize <- cug.test(flobusiness,centralization,
+   FUN.arg=list(FUN=betweenness),mode="graph",cmode="size")
> rEdges <- cug.test(flobusiness,centralization,
+   FUN.arg=list(FUN=betweenness),mode="graph",cmode="edges")
> rDyad <- cug.test(flobusiness,centralization,
+   FUN.arg=list(FUN=betweenness),mode="graph",cmode="dyad.census")
> names(rSize)
[1] "obs.stat" "rep.stat" "mode"     "diag"     "cmode"    "plteobs"  "pgteobs" 
[8] "reps"    
> rSize

Univariate Conditional Uniform Graph Test

Conditioning Method: size 
Graph Type: graph 
Diagonal Used: FALSE 
Replications: 1000 

Observed Value: 0.2057143 
Pr(X>=Obs): 0.001 
Pr(X<=Obs): 0.999 
> 
> # Aggregate results
> Betweenness <- c(rSize$obs.stat,rEdges$obs.stat,rDyad$obs.stat)
> PctGreater <- c(rSize$pgteobs,rEdges$pgteobs,rDyad$pgteobs)
> PctLess <- c(rSize$plteobs,rEdges$plteobs,rDyad$plteobs)
> report <- cbind(Betweenness, PctGreater, PctLess)
> rownames(report) <- c("Size","Edges","Dyads")
> report
      Betweenness PctGreater PctLess
Size    0.2057143      0.001   0.999
Edges   0.2057143      0.713   0.289
Dyads   0.2057143      0.738   0.263
> gplot(flobusiness,gmode="graph")
> gplot(flobusiness,gmode="graph",displaylabels=TRUE)
> 
> par(mfrow=c(1,3))
> plot(rSize, main="Betweenness \nConditioned on Size" )
> plot(rEdges, main="Betweenness \nConditioned on Edges" )
> plot(rDyad, main="Betweenness \nConditioned on Dyads" )
> par(mfrow=c(1,1))
</code>

===== QAP =====

<code>
> help(qaptest)
> gcor(flobusiness,flomarriage)
[1] 0.3718679
> (rCor <- qaptest(list(flobusiness,flomarriage),gcor, 
+   g1=1,g2=2,reps=1000))

QAP Test Results

Estimated p-values:
        p(f(perm) >= f(d)): 0.001 
        p(f(perm) <= f(d)): 1 

> 
> plot(rCor, xlim=c(-0.25, 0.4))
</code>
===== Preparing data =====

<html><!--
[[http://thomasgrund.weebly.com/teaching.html|Grund]]
--></html>
C:\Users\batagelj\Documents\papers\2018\moskva\NetR\doc\R\grund\


https://raw.githubusercontent.com/bavla/netstat/master/data/LGang/

<html><!--
http://thomasgrund.weebly.com/teaching.html
--></html>
The LGang data includes attributes and the co-offending network of a youth gang. Data were collected by James Densley (Metropolitan University) and Thomas Grund (UCD).

<code>
> setwd("C:/Users/batagelj/Documents/papers/2018/moskva/NetR/nets")
> library(sna)
> list.files()
[1] "campnet" "grund"   "padgett" "rivers"  "test"   
> # LG.mat <- read.csv("./grund/LGang-matrix.csv",header=TRUE)
> LG.mat <- read.csv("https://raw.githubusercontent.com/bavla/netstat/master/data/LGang/LGang-matrix.csv",header=TRUE)
> head(LG.mat)
  ID X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16 X17 X18 X19 X20 X21 X22 X23 X24 X25 X26
1  1  0  1  1  2  1  1  2  3  2   2   3   1   0   0   0   0   1   1   0   2   2   3   1   0   1   0
...
> LG.mat <- LG.mat[,-1]  # remove ID
> LG.mat[LG.mat==1]<-0   # recode
> LG.mat[LG.mat>=2]<-1
> LG.net <- network(as.matrix(LG.mat),directed=FALSE) # convert to network
> # LG.attrib <- read.csv("./grund/LGang-attributes.csv")
> LG.attrib <- read.csv("https://raw.githubusercontent.com/bavla/netstat/master/data/LGang/LGang-attributes.csv")
> head(LG.attrib)
  ID year age  birthplace residence arrests convictions prison music rank rankrev core
1  1 1989  20 West Africa         0      16           4      1     1    1       5    1
2  2 1989  20   Caribbean         0      16           7      1     0    2       4    1
3  3 1990  19   Caribbean         0      12           4      1     0    2       4    1
4  4 1988  21   Caribbean         0       8           1      0     0    2       4    1
5  5 1985  24   Caribbean         0      11           3      0     0    2       4    1
6  6 1984  25          UK         1      17          10      0     0    2       4    1
> LG.net %v% "ethnicity" <- as.character(LG.attrib$birthplace)
> LG.net %v% "prison" <- LG.attrib$prison
> plot(LG.net)
> plot(LG.net,vertex.col="ethnicity")
> plot(LG.net,vertex.col="prison")
> prisonsize<- 1 + LG.net %v% "prison"  # display as vertex size
> plot(LG.net, vertex.cex=prisonsize, vertex.col="ethnicity")
</code>

===== Changes of transitivity with density =====

<code>
> # changes for different densities
> (d<-seq(from=0.05,to=0.85,by=0.1))
[1] 0.05 0.15 0.25 0.35 0.45 0.55 0.65 0.75 0.85
> den <-lapply(d, function(x) gden(rgraph(100,m=50, tprob= x)))
> class(den)
[1] "list"
> length(den)
[1] 9
> class(den[[1]])
[1] "numeric"
> den[[1]]
 [1] 0.05161616 0.05282828 0.04979798 0.04727273 0.05161616 0.05111111 0.05252525 0.04969697 0.04777778
[10] 0.04969697 0.04919192 0.05181818 0.05424242 0.04939394 0.04909091 0.05000000 0.04727273 0.04626263
[19] 0.04959596 0.05474747 0.04858586 0.04888889 0.05181818 0.05171717 0.04868687 0.05060606 0.05343434
[28] 0.04585859 0.04777778 0.04969697 0.05242424 0.05161616 0.04989899 0.04929293 0.05090909 0.04898990
[37] 0.04808081 0.05161616 0.04616162 0.04949495 0.05020202 0.04737374 0.04737374 0.04656566 0.05171717
[46] 0.05080808 0.05040404 0.04939394 0.05313131 0.05060606
> library(beanplot)
> beanplot(den, names = d)
> trans <-lapply(d, function(x) gtrans(rgraph(100,m=50, tprob= x), measure="weakcensus"))
> beanplot(trans, names = d)
> # changes of betweenness with density
> d <- seq(0.05,0.5,by=0.05)
> n <- 100
> r <- 50
> 
> bw <- lapply(d, 
+     function(x){ temp <- vector('list',r)
+       for (i in 1:r) temp[[i]] <- betweenness(rgraph(n,tprob=x),rescale=TRUE)
+       return(temp)
+   })
> beanplot(dat, names=d, what=c(1,1,1,0))
</code>

===== Testing =====
==== Transitivity ====

Are there fewer unclosed triangles in the LGang network than one would expect by chance, given the size and number of edges that are there?

<code>
> (LG.trans <- gtrans(LG.net))
[1] 0.3635371
> (LG.density <- gden(LG.net))
[1] 0.092942
> (LG.size <- network.size(LG.net))
[1] 54
> (LG.nedges <- network.edgecount(LG.net)) 
[1] 133
> ERgphs <- rgnm(n=200, nv=LG.size, m=LG.nedges, mode="graph")
> ERgph1 <- network(ERgphs[1,,],directed=F)
> plot(ERgph1)
> plot(LG.net)
> gden(ERgph1)
[1] 0.092942
> gden(LG.net)
[1] 0.092942
> ERtrans <- gtrans(ERgphs)
> plot(density(ERtrans))
> plot(density(ERtrans),xlim=c(0,0.38))
> abline(v=LG.trans,col="red",lw=2)
> # packed test
> cg <- cug.test(LG.net,"gtrans",mode="graph",cmode=c("edges"),reps=10000)
> cg

Univariate Conditional Uniform Graph Test

Conditioning Method: edges 
Graph Type: graph 
Diagonal Used: FALSE 
Replications: 10000 

Observed Value: 0.3635371 
Pr(X>=Obs): 0 
Pr(X<=Obs): 1 

> plot(cg)
</code>

The probability for a random network with the same number of edges as the LGang network to have more triangles than the LGang network is practically 0. So, relative to a random network with the same amount of edges there would appear to be high clustering in the LGang network.

<code>
> cg2 <- cug.test(LG.net,gtrans,cmode="dyad.census")
> cg2

Univariate Conditional Uniform Graph Test

Conditioning Method: dyad.census 
Graph Type: digraph 
Diagonal Used: FALSE 
Replications: 1000 

Observed Value: 0.3635371 
Pr(X>=Obs): 0 
Pr(X<=Obs): 1 

> plot(cg2)
</code>

===== Permutation test QAP =====

<code>
> ety <- get.vertex.attribute(LG.net,"ethnicity")
> ety
 [1] "West Africa" "Caribbean"   "Caribbean"   "Caribbean"   "Caribbean"   "UK"          "East Africa"
 [8] "West Africa" "West Africa" "West Africa" "West Africa" "UK"          "UK"          "UK"         
...
> Mety <- 1 * outer(ety, ety, "==") # 1 * for conversion to integer
> Mety[1:10,1:10]
      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
 [1,]    1    0    0    0    0    0    0    1    1     1
 [2,]    0    1    1    1    1    0    0    0    0     0
 [3,]    0    1    1    1    1    0    0    0    0     0
 [4,]    0    1    1    1    1    0    0    0    0     0
 [5,]    0    1    1    1    1    0    0    0    0     0
 [6,]    0    0    0    0    0    1    0    0    0     0
 [7,]    0    0    0    0    0    0    1    0    0     0
 [8,]    1    0    0    0    0    0    0    1    1     1
 [9,]    1    0    0    0    0    0    0    1    1     1
[10,]    1    0    0    0    0    0    0    1    1     1
> gcor(LG.net,Mety)
[1] 0.1249278
> LG.perm <- network(rmperm(LG.net),directed=F)

> gden(LG.perm)
[1] 0.092942
> gden(LG.net)
[1] 0.092942
> list.vertex.attributes(LG.net)
[1] "ethnicity"    "na"           "prison"       "vertex.names"
> list.vertex.attributes(LG.perm)
[1] "na"           "vertex.names"
> LG.net %v% "vertex.names"
 [1] "X1"  "X2"  "X3"  "X4"  "X5"  "X6"  "X7"  "X8"  "X9"  "X10" "X11" "X12" "X13" "X14" "X15" "X16" "X17"
[18] "X18" "X19" "X20" "X21" "X22" "X23" "X24" "X25" "X26" "X27" "X28" "X29" "X30" "X31" "X32" "X33" "X34"
[35] "X35" "X36" "X37" "X38" "X39" "X40" "X41" "X42" "X43" "X44" "X45" "X46" "X47" "X48" "X49" "X50" "X51"
[52] "X52" "X53" "X54"
> LG.perm %v% "vertex.names"
 [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
[35] 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54
> LG.perm %v% "vertex.names" <- LG.net %v% "vertex.names"
> plot(LG.net,displaylabels=TRUE)
> plot(LG.perm,displaylabels=TRUE)
> gcor(LG.net, LG.perm)
[1] 0.03845129
> gcor(LG.perm, Mety)
[1] -0.03857968
> qap <- qaptest(list(LG.net,Mety),g1=1,g2=2,FUN=gcor,reps=5000)
> summary(qap)

QAP Test Results

Estimated p-values:
        p(f(perm) >= f(d)): 2e-04 
        p(f(perm) <= f(d)): 0.9998 

Test Diagnostics:
        Test Value (f(d)): 0.1249278 
        Replications: 5000 
        Distribution Summary:
                Min:     -0.09659847 
                1stQ:    -0.02275637 
                Med:     -0.001658628 
                Mean:    0.0001832054 
                3rdQ:    0.01943912 
                Max:     0.1354767 

> plot(qap)
</code>
Our actual test-statistic (here, correlation) is much higher than what we would see if we were to scramble the network.

<code>
> Netr1<-rgraph(100, tprob=0.1)
> Netr2<-rgraph(100, tprob=0.5)
> plot(network(Netr1))
> 
> qap<-qaptest(list(Netr1,Netr2),g1=1,g2=2,FUN=gcor,reps=5000)
> summary(qap)

QAP Test Results

Estimated p-values:
        p(f(perm) >= f(d)): 0.3958 
        p(f(perm) <= f(d)): 0.6132 

Test Diagnostics:
        Test Value (f(d)): 0.002507697 
        Replications: 5000 
        Distribution Summary:
                Min:     -0.03736669 
                1stQ:    -0.006954022 
                Med:     -0.0001956514 
                Mean:    -0.0001938942 
                3rdQ:    0.006562719 
                Max:     0.03562371 

> plot(qap)
</code>
Our observed correlation falls right into the distribution of correlations we would expect by chance, controlling for the underlying structure. Hence, the correlation between the two random networks is not significant.

===== Regression =====

<code>
> ties <- as.matrix(LG.net)
> diag(ties) <- NA
> ties <- ties[upper.tri(ties)]
> Maty <- Mety
> Maty <- Maty[upper.tri(Maty)]
> dyads <- data.frame(ties,Maty)
> dyads[1:10,]
   ties Maty
1     0    0
2     0    0
3     1    1
4     1    0
5     0    1
6     1    1
7     0    0
8     0    1
9     1    1
10    1    1
> results <- glm(ties~Maty,family="binomial",data=dyads)
> summary(results)

Call:
glm(formula = ties ~ Maty, family = "binomial", data = dyads)

Deviance Residuals: 
    Min       1Q   Median       3Q      Max  
-0.5679  -0.5679  -0.3794  -0.3794   2.3096  

Coefficients:
            Estimate Std. Error z value Pr(>|z|)    
(Intercept)  -2.5953     0.1239 -20.946  < 2e-16 ***
Maty          0.8523     0.1844   4.622  3.8e-06 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 885.19  on 1430  degrees of freedom
Residual deviance: 864.48  on 1429  degrees of freedom
AIC: 868.48

Number of Fisher Scoring iterations: 5

>
</code>
 The null hypothesis is that networks are random with the right number of edges (so the dyads are independent). The x1 estimate (=beta coefficient) is positive  at 0.85 and this effect is significant with p = 4.14...e-06 = very small p. Basically, the regression says that a tie between two gang members is Exp(b) = 2.34 times more likely when the two gang members have the same ethnicity. There is clear evidence for ethnic homophily in co-offending.

<code>
help(netlogit)
> dyQap <- netlogit(LG.net,Mety,mode="graph",nullhyp="qapy",reps=100)
> dyQap

Network Logit Model

Coefficients:
            Estimate   Exp(b)     Pr(<=b) Pr(>=b) Pr(>=|b|)
(intercept) -2.5952547 0.07462687 0.91    0.09    0.91     
x1           0.8522854 2.34500000 1.00    0.00    0.00     

Goodness of Fit Statistics:

Null deviance: 1983.787 on 1431 degrees of freedom
Residual deviance: 864.4812 on 1429 degrees of freedom
Chi-Squared test of fit improvement:
         1119.306 on 2 degrees of freedom, p-value 0 
AIC: 868.4812   BIC: 879.0135 
Pseudo-R^2 Measures:
        (Dn-Dr)/(Dn-Dr+dfn): 0.4388909 
        (Dn-Dr)/Dn: 0.5642268 

> 


</code>
==== Interactive coordinates ====

<code>
> xy <- gplot(LG.net,displaylabels=TRUE,label.cex=0.7,usearrows=FALSE,interactive=TRUE)
> gplot(LG.net,coord=xy,displaylabels=TRUE,label.cex=0.7,usearrows=FALSE)
</code>


======  ======
\\



