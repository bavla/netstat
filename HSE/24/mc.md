# Markov chains data

## Black Friday

[The Black Friday Puzzle](https://www.countbayesie.com/blog/2015/11/21/the-black-friday-puzzle-understanding-markov-chains)

<img src="https://github.com/bavla/netstat/blob/master/HSE/24/pics/blackFridayGph.png" width=350>

```
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
```


<hr/>

[HSE](../2024.md); [Docs](doc.md); [Matrices](csnet2.md)

