# ====== Markov chains data ======

===== Black Friday =====

[[https://www.countbayesie.com/blog/2015/11/21/the-black-friday-puzzle-understanding-markov-chains|The Black Friday Puzzle]]

{{ru:hse:snet:pics:blackfridaygph.png?350}}
<code>
> bf <- c(
+ 0.5 , 0.1, 0.1, 0.05, 0.2,
+ 0.2 , 0.3, 0.2, 0.15, 0.1,
+ 0.15, 0.2, 0.3, 0.3 , 0.1,
+ 0.1 , 0.3, 0.2, 0.4 , 0.1,
+ 0.05, 0.1, 0.2, 0.1 , 0.5 )
> BF <- matrix(bf,nrow=5)
> rownames(BF) <- colnames(BF) <- c("Books", "Children",
+   "Puzzles", "Toys", "Music")
> 
> BF
         Books Children Puzzles Toys Music
Books     0.50     0.20    0.15  0.1  0.05
Children  0.10     0.30    0.20  0.3  0.10
Puzzles   0.10     0.20    0.30  0.2  0.20
Toys      0.05     0.15    0.30  0.4  0.10
Music     0.20     0.10    0.10  0.1  0.50
</code>
