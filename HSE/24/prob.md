# Problems 

Solve the problem with your number with the Monte Carlo method. If you are able to solve it also theoretically compare both results.

  1. We have three containers. The first contains n white balls, the second contains  n  black balls, and in the third contains  n  red balls. Randomly take one ball from each container and then place the ball from the first container into the second container, the ball from the second to the third, and the ball from the third into the first container. What is the expected number of balls of each color in containers after  k  steps? Solve for n = 5, 10, 25, 50, 100 and k = 1, ..., 10, 20, 50, 100.
  1. The players  A  and  B  have at the beginning 12 coins each. They play with three dices. After each throw of dices (who throws does not matter), if the sum of the points is equal to 11, the player A gives a coin to B; but if the sum of the points is 14, the player B  gives the coin to  A. A player who first collects all the coins wins. How many throws take the average game? What is the probability that the player A wins?
  1. Two decks of n  cards are shuffled and we take away the top  m cards. What is the expected number of pairs in the remainder? Solve for n = 52 and all possible values of m. [Hint 2](#hint-2)
  1. n  people are present at the meeting. What is the probability of having at least two birthdays on the same day of the year? What is the probability of three people? Solve for n = 10, 50, 100.
  1. n  people are present at the meeting. What is the probability of having at least two birthdays on the same day of the year? What is the probability for two pairs, three pairs (each pair for its own date)? Solve for n = 10, 50, 100.
  1. n hunters came to the inn and put off their hats. When they were leaving, the lights went out. That's why everyone took some hat. What's the probability that nobody took his own hat? Solve for n = 5, 10, 15, 20, 50, 100.
  1. The player travels along the tape, which consists of  n  consecutive cells. He starts in cell 1. At each step, a coin (number/face) is thrown. If the face falls, he is advancing for one cell, otherwise by two. What is the probability that he will enter the cell k, k = 1, .., n ? What is the average number of steps up to the end of the game, if we require the player to end the game in the field  n (if the move is impossible he stays in the current cell)? Solve it for n = 100.
  1. Given is a convex shape L. What is the probability that the randomly selected four points in L will determine a convex quadrilateral? Solve for L = circle. [Hint 1](#hint-1)
  1. Given is a convex shape L. What is the probability that the randomly selected four points in L will determine a convex quadrilateral? Solve for L = equilateral triangle. [Hint 1](#hint-1)
  1. Given is a convex shape L. What is the probability that the randomly selected four points in L will determine a convex quadrilateral? Solve for L = square. [Hint 1](#hint-1)
  1. Given is a convex shape L. What is the probability that the randomly selected four points in L will determine a convex quadrilateral? Solve for L = semicircle. [Hint 1](#hint-1)
  1. Eight boys and seven girls bought each his/her own ticket in the same row at the cinema. What is the average number of pairs in a row? For example: there are 9 pairs in the BBGGBBGBGBGBBGG sequence.
  1. Which event is more likely
    *  to throw at least one 6 in six rounds; or
    *  to throw at least three 6 in eighteen rounds?
  1. The drunk starts to melt at random with the same probability in directions S, N, E, W. What is the probability that it will be in n, n = 1, ..., 50 steps returned to the starting point?
  1. A  1m  long stick is randomly cut into two parts. What is the average length of shorter side? What is the average ratio between the length of the shorter and longer side?
  1. A  1m  long stick is randomly cut into three parts. What is the probability that the parts form the sides of the triangle?
  1. Shuffle the deck of 52 cards. Then we drop the cards from the top of the deck until we throw off the first ace. How many cards are discarded on average? [Hint 2](#hint-2)
  1. In the confectionery in the dough, n nuts are mixed. From the dough, m cookies are made. How many nuts must at least be used on average so that there will be at most p percent of cookies without nuts in them? (m = 100, p = 5)
  1. In the confectionery in the dough, n nuts are mixed. From the dough, m cookies are made. How many nuts must at least be used on average so that at least q percent of cookies contain at least 2 nuts? (m = 100, q = 90)
  1. Points A and B are uniformly distributed within the unit circle. Determine the probability that the length of the line AB is shorter than the radius of the circle r = 1.
  1. Repeat in R the [Buffon's needle](https://en.wikipedia.org/wiki/Buffon%27s_needle_problem) computation of the number π.
  1. In the unit square [0,1]x[0,1] a pair of random line segments are selected. What is the probability that they intersect? [1](https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/), [2](https://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect)
  1. Shuffle the deck of 52 cards. How many cards separate the first and the second ace in the deck on average? [Hint 2](#hint-2)


## Hint 1 

The points A, B, C, D form a convex quadrilateral iff no point is inside the triangle determined by the remaining three points. [1](https://stackoverflow.com/questions/2049582/how-to-determine-if-a-point-is-in-a-2d-triangle), [2](http://totologic.blogspot.com/2014/01/accurate-point-in-triangle-test.html), [3](https://stackoverflow.com/questions/9513107/find-if-4-points-form-a-quadrilateral)

Another criterion is that one of the pairs of segments  (AB,CD), (AD,BC), (AC,BD) intersect. See also problem 23.	

## Hint 2 

A [standard deck of cards](https://en.wikipedia.org/wiki/Standard_52-card_deck) can be created as follows
<code>
> suits <- c("♣", "♦", "♥", "♠")
> ranks <- c(2:9,"X", "J", "Q", "K", "A")
> deck <- as.vector(outer(ranks,suits,"paste",sep=""))
> deck
 [1] "2♣" "3♣" "4♣" "5♣" "6♣" "7♣" "8♣" "9♣" "X♣" "J♣" "Q♣" "K♣" "A♣"
[14] "2♦" "3♦" "4♦" "5♦" "6♦" "7♦" "8♦" "9♦" "X♦" "J♦" "Q♦" "K♦" "A♦"
[27] "2♥" "3♥" "4♥" "5♥" "6♥" "7♥" "8♥" "9♥" "X♥" "J♥" "Q♥" "K♥" "A♥"
[40] "2♠" "3♠" "4♠" "5♠" "6♠" "7♠" "8♠" "9♠" "X♠" "J♠" "Q♠" "K♠" "A♠"
</code>
To get two characters card codes I replaced the number 10 with the roman numeral X.

To shuffle the deck use
<code>
> deck <- sample(deck)
> deck
 [1] "8♦" "7♣" "K♣" "5♠" "2♣" "9♦" "6♦" "X♥" "J♣" "5♥" "K♥" "3♥" "Q♣"
[14] "J♥" "4♣" "9♥" "5♣" "2♥" "A♦" "J♠" "A♥" "Q♠" "4♦" "9♠" "A♣" "6♥"
[27] "3♣" "3♠" "8♠" "8♥" "4♥" "X♦" "A♠" "6♠" "5♦" "Q♥" "8♣" "7♥" "K♦"
[40] "7♦" "7♠" "9♣" "4♠" "K♠" "6♣" "3♦" "Q♦" "X♠" "2♦" "2♠" "X♣" "J♦"
</code>
and to find the position of selected card(s) use
<code>
> grep("Q♥",deck)
[1] 36
> grep("A",deck)
[1] 19 21 25 33
</code>

<HTML><!--
https://demonstrations.wolfram.com/ConcaveRandomQuadrilateralsFromFourPointsInADisk/

17. 53/5 = 10.6 
   https://www.quora.com/What-is-the-expected-number-of-cards-that-need-to-be-turned-over-in-a-regular-52-card-deck-in-order-to-see-the-first-ace
--></HTML>

\\

[[ru:hse:snet22:stu|Snet projects]] 

