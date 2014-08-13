Control Flow and loops in R ############################

Control Flow
============

The standard if else
--------------------

.. code:: r

    p.test <- function(p) {
        if (p <= 0.05) 
            print("yeah!!!!") else if (p >= 0.9) 
            print("high!!!!") else print("somewhere in the middle")
    }

Now pick a number and put it in ``p.test``

.. code:: r

    p.test(0.5)

::

    ## [1] "somewhere in the middle"

ifelse()
========

A better and vectorized way of doing this is ``ifelse(test, yes, no)``
function. ifelse() is far more useful as it is vectorized.

.. code:: r

    p.test.2 <- function(p) {
        ifelse(p <= 0.05, print("yippee"), print("bummer, man"))
    }

Test this with the following sequence. See what happens if you use
``if`` vs. ``ifelse()``.

.. code:: r

    x <- runif(10, 0, 1)
    x

::

    ##  [1] 0.27332 0.14155 0.89000 0.07041 0.79419 0.25013 0.02324 0.86766
    ##  [9] 0.41114 0.56165

Now try it with ``p.test()`` (uses ``if``).

.. code:: r

    p.test(x)

::

    ## Warning: the condition has length > 1 and only the first element will be used
    ## Warning: the condition has length > 1 and only the first element will be used

::

    ## [1] "somewhere in the middle"

Now try it with ``p.test.2()``

.. code:: r

    p.test.2(x)

::

    ## [1] "yippee"
    ## [1] "bummer, man"

::

    ##  [1] "bummer, man" "bummer, man" "bummer, man" "bummer, man" "bummer, man"
    ##  [6] "bummer, man" "yippee"      "bummer, man" "bummer, man" "bummer, man"

Other vectorized ways of control flow.
======================================

There are many times that you may think you need to use an if with
(iterating with a for loop... see below), or ifelse, but there may be
far better ways.

For instance, say you are doing some simulations for a power analysis,
and you want to know how often your simulation gives you a p-value less
than 0.05.

.. code:: r

    p.1000 <- runif(n = 1000, min = 0, max = 1)

The line above generates 1000 random values between 0-1, which we will
pretend are our p-values for differential expression from our
simulation.

You may try and count how often it less than 0.05

.. code:: r

    p.ifelse <- ifelse(p.1000 < 0.05, 1, 0)  # If it is less than 0.05, then you get a 1, otherwise 0. 

Our approximate false positives. Should be close to 0.05

.. code:: r

    sum(p.ifelse)/length(p.1000)

::

    ## [1] 0.059

In R, think index!
~~~~~~~~~~~~~~~~~~

However the best and fastest way to accomplish this is to use the index,
by setting up the Boolean (TRUE/FALSE) in the index of the vector.

.. code:: r

    length(p.1000[p.1000 < 0.05])/length(p.1000)

::

    ## [1] 0.059

Same number, faster and simpler computation.

Simple loops
============

while() function..
------------------

I tend to avoid these, so you will not see them much here

.. code:: r

    i <- 1
    while (i <= 10) {
        print(i)
        i <- i + 0.5
    }

::

    ## [1] 1
    ## [1] 1.5
    ## [1] 2
    ## [1] 2.5
    ## [1] 3
    ## [1] 3.5
    ## [1] 4
    ## [1] 4.5
    ## [1] 5
    ## [1] 5.5
    ## [1] 6
    ## [1] 6.5
    ## [1] 7
    ## [1] 7.5
    ## [1] 8
    ## [1] 8.5
    ## [1] 9
    ## [1] 9.5
    ## [1] 10

for loop
========

If I run a loop I most often use ``for(){}`` automatically iterates
across a list (in this case the sequence from 1:10).

.. code:: r

    for (i in 1:10) {
        print(i)
    }

::

    ## [1] 1
    ## [1] 2
    ## [1] 3
    ## [1] 4
    ## [1] 5
    ## [1] 6
    ## [1] 7
    ## [1] 8
    ## [1] 9
    ## [1] 10

If you do not want to use integers, how might you do it using the for()?

.. code:: r

    for (i in seq(from = 1, to = 5, by = 0.5)) {
        print(i)
    }

::

    ## [1] 1
    ## [1] 1.5
    ## [1] 2
    ## [1] 2.5
    ## [1] 3
    ## [1] 3.5
    ## [1] 4
    ## [1] 4.5
    ## [1] 5

behavior of strings.
^^^^^^^^^^^^^^^^^^^^

Using strings is a bit more involved in R, compared to other languages.
For instance the following does not do what you want::

.. code:: r

    for (letter in "word") {
        print(letter)
    }

::

    ## [1] "word"

(try letters for a hoot.)

Instead in R, we have to split the word "word" into single characters
using strsplit(), i.e::

.. code:: r

    strsplit("word", split = "")

::

    ## [[1]]
    ## [1] "w" "o" "r" "d"

So for the for loop we would do the following:
==============================================

.. code:: r

    for (letter in strsplit("word", split = "")) {
        print(letter)
    }

::

    ## [1] "w" "o" "r" "d"

More avoiding loops
===================

Many would generate random numbers like so.

.. code:: r

    for (i in 1:100) {
        print(rnorm(n = 1, mean = 0, sd = 1))
    }

::

    ## [1] -0.1837
    ## [1] -0.9313
    ## [1] 1.648
    ## [1] -0.6964
    ## [1] 0.2112
    ## [1] 0.3441
    ## [1] 1.036
    ## [1] 0.7439
    ## [1] 0.5859
    ## [1] -0.6087
    ## [1] -0.4014
    ## [1] 1.44
    ## [1] -0.3906
    ## [1] -1.861
    ## [1] -0.739
    ## [1] -1.204
    ## [1] 0.07794
    ## [1] -1.65
    ## [1] 1.261
    ## [1] 0.6753
    ## [1] 0.6736
    ## [1] 0.3238
    ## [1] -1.316
    ## [1] 0.2965
    ## [1] 1.499
    ## [1] 0.4326
    ## [1] 0.4488
    ## [1] 0.8873
    ## [1] -1.304
    ## [1] -0.347
    ## [1] 0.3491
    ## [1] 0.24
    ## [1] 0.1425
    ## [1] -0.2785
    ## [1] -0.5072
    ## [1] -1.775
    ## [1] -0.04051
    ## [1] 0.9452
    ## [1] 0.3322
    ## [1] -0.01994
    ## [1] -0.2308
    ## [1] -0.4053
    ## [1] -0.5685
    ## [1] -1.631
    ## [1] -0.1484
    ## [1] 0.434
    ## [1] 1.653
    ## [1] 1.57
    ## [1] 0.1308
    ## [1] -1.059
    ## [1] -0.7157
    ## [1] -0.8316
    ## [1] 0.06561
    ## [1] 0.8243
    ## [1] 0.1841
    ## [1] 1.048
    ## [1] 0.1612
    ## [1] -0.9553
    ## [1] -0.7569
    ## [1] -0.288
    ## [1] -1.837
    ## [1] 0.7301
    ## [1] -2.103
    ## [1] -1.869
    ## [1] -1.298
    ## [1] -1.077
    ## [1] -0.2139
    ## [1] -0.9419
    ## [1] 0.4694
    ## [1] -1.344
    ## [1] -0.08514
    ## [1] -2.055
    ## [1] -0.803
    ## [1] -0.7281
    ## [1] 1.778
    ## [1] -1.116
    ## [1] 1.33
    ## [1] 0.1535
    ## [1] -2.897
    ## [1] 0.7305
    ## [1] 1.228
    ## [1] 1.697
    ## [1] -0.8183
    ## [1] -1.013
    ## [1] -0.634
    ## [1] -0.942
    ## [1] -0.3395
    ## [1] 0.1396
    ## [1] 1.022
    ## [1] 0.9868
    ## [1] -0.7778
    ## [1] 1.075
    ## [1] -0.1029
    ## [1] 0.2644
    ## [1] 0.01165
    ## [1] 0.8025
    ## [1] -1.24
    ## [1] -0.8865
    ## [1] 0.981
    ## [1] 0.5333

We are cycling through and generating one random number at each
iteration. Look at the indices, and you can see we keep generating
vectors of length 1.

better/cleaner/faster to generate them all at one time

.. code:: r

    rnorm(n = 100, mean = 0, sd = 1)

::

    ##   [1] -0.08683 -1.55262 -1.16909  0.30451 -1.14555  0.76682  0.12643
    ##   [8] -0.61174 -0.29103 -0.10707 -0.03397 -0.05926  0.27294  1.32693
    ##  [15] -0.53284  1.83234  0.43959 -0.88991  0.25383  0.96709 -0.23210
    ##  [22] -1.00190 -1.32289  1.80030  1.15272 -1.82907  0.75989  1.35966
    ##  [29]  0.53943  0.01429 -0.58707 -0.11886 -0.70367 -2.38988  0.08033
    ##  [36] -0.22795 -0.62166 -0.19832 -1.95990 -0.85127  0.94236  0.37771
    ##  [43]  0.32617 -0.08393 -0.54506 -2.58781 -0.58433  0.20985 -0.41613
    ##  [50]  0.60527  0.51713  1.57950 -0.61079 -0.28564 -0.16444  0.55007
    ##  [57]  0.57258  0.58513 -0.86728 -0.81185 -0.29333 -1.23935  0.46169
    ##  [64] -1.53586 -0.32583  0.17629 -0.85579  1.04989  1.22120  1.53359
    ##  [71] -2.37276  1.44393  1.47506  0.40110 -0.10157  0.35485 -0.72068
    ##  [78] -1.27910  0.63152 -0.65216  1.60160  0.27109  0.50904 -1.00531
    ##  [85]  0.76743 -0.78954 -0.01159  1.06944  1.15661 -0.91031  1.54919
    ##  [92] -0.84334  2.19994  0.26716  0.02081  0.53577  0.07840 -0.79387
    ##  [99] -1.18941  1.24745

What if we wanted to put all of these numbers in a vector?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**The not advisable approach**

First we initialize a vector to store all of the numbers. Why do we
initialize this vector first?

.. code:: r

    n <- 1e+05
    x <- rep(NA, n)

The step above creates a vector of n NA's. They will be replaced sequentially with the random numbers as we generate them (using a function like the above one).
================================================================================================================================================================

.. code:: r

    head(x)

::

    ## [1] NA NA NA NA NA NA

Now we run the for loop.

.. code:: r

    for (i in 1:n) {
        x[i] <- rnorm(n = 1, mean = 0, sd = 1)
    }

for each ``i`` in the index, one number is generated, and placed in x

.. code:: r

    head(x)

::

    ## [1]  0.2848 -0.5432  1.1391 -1.0901  0.8515  0.5490

However this is computationally inefficient in R. Which has vectorized
operations.

.. code:: r

    system.time(

    for (i in 1:n){
        x[i] <- rnorm(n=1, mean=0, sd=1)})

::

    ##    user  system elapsed 
    ##   0.562   0.023   0.584

We can also use the replicate function to do the same thing. Easier
syntax to write.

.. code:: r

    system.time(z <- replicate(n, rnorm(n = 1, mean = 0, sd = 1)))

::

    ##    user  system elapsed 
    ##   0.561   0.035   0.841

This is ~20% faster.

The way to do it
^^^^^^^^^^^^^^^^

However, since R is vectorized, both of the will be far slower than:

.. code:: r

    system.time(y <- rnorm(n, 0, 1))

::

    ##    user  system elapsed 
    ##   0.010   0.000   0.011

About 65 times faster than the for loop

The general rule in R is that loops are slower than the apply family of
functions (for small to medium data sets, not true for very large data)
which are slower than vectorized computations.
