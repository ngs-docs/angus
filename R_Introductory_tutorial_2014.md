R Tutorial for NGS2014
######################

In this tutorial I will introduce `R`, a programming language that is currently the *lingua franca* for data analysis (although `python` has many powerful data analysis tools through `numpy`, `scipy` and other libraries.

## What is **R**?
R is a bit of funny programming language, with a [funny history](http://en.wikipedia.org/wiki/R_(programming_language)). While many researchers with a statistical bent really appreciate aspects of the syntax, and how it works, there are a great number of idiosyncracies, and "features" of the language that [drive many programmers crazy](http://www.burns-stat.com/pages/Tutor/R_inferno.pdf). It is also known to be a bit of a memory hog, so big data sets **can** be an issue to work with if you are not aware of the appropriate approaches to deal with. Despite this, there are some features of `R` that make it great for doing data analyses of all kinds, and there are a huge number of libraries of tools covering all aspects of statistics, genomics, bioinformatics, etc... In many cases new methods are associated with at least a basic implementation in `R` if not a full blown library.

## Installing R

Many of you may already have R installed on your local machine. If you would like to install R on your computer just go to http://www.r-project.org/ click download and follow for your OS. For linux machines, it is best to instead use your package manager (`yum` or `apt-get` for instance).

If you would like to run this on an Amazon EC2 instance, [set up and log into your instance](http://angus.readthedocs.org/en/2014/day1.html) as you did in the earlier tutorial.

R generally does not come pre-installed on most Linux distributions including Ubuntu (Debian) Linux, which we are using on our Amazon EC2 instance but it is very easy to install:

```
apt-get install r-base
```

You may wish to [link and mount your instance to dropbox](http://angus.readthedocs.org/en/2014/amazon/installing-dropbox.html) as we did in the earlier tutorials.

It is also important to point out that `R` has been packaged in a number of different easy to use ways. For instance, many scientists really like [R-studio](http://www.rstudio.com/) which is also a free easy to use IDE that makes it easier for you to implement *version control* and do reproducible research. I personally do not use `R` via R-studio, but feel free to try it out on your local machine if you want to.

For the remainder of the tutorial I am assuming that we will all be running `R` on your Amazon EC2 instance though.

#R

## What is R, really....
While many of you are used to seeing `R` through some GUI (whether on windows, OSX,etc).. It is fundamentally just a command line program like what you have been using at the command line or when you call `python`, etc.. The GUI versions just add some additional functionality.

At your shell prompt just type:
```
R
```

You now have a new command prompt (likely `>`) which means you are in `R`. Indeed you probably see something that looks like:
```
##  R version 3.0.1 (2013-05-16) -- "Good Sport"
##  Copyright (C) 2013 The R Foundation for Statistical Computing
## Platform: x86_64-apple-darwin10.8.0 (64-bit)

##  R is free software and comes with ABSOLUTELY NO WARRANTY.
##  You are welcome to redistribute it under certain conditions.
##  Type 'license()' or 'licence()' for distribution details.
##
##   Natural language support but running in an English locale
##
##  R is a collaborative project with many contributors.
##  Type 'contributors()' for more information and
##  'citation()' on how to cite R or R packages in publications.

##  Type 'demo()' for some demos, 'help()' for on-line help, or
## 'help.start()' for an HTML browser interface to help.
##  Type 'q()' to quit R.
##  
##  >
```

### Where to find stuff in GUI R (navigating GUI-R)
I will quickly (on the Mac OSX version) show you some of navigating the GUI.

## How to close `R`
To quit R:
```r
q() 
```
**NOTE to MAC (OSX) users**: This may no longer work (R V2.11.+ ) from the Mac R GUI... See below in OSX specific notes..

**NOTE**: For the moment when it asks you to save workspace image, say no.

##R Basics
Ok, return to `R` (type `R` at unix command prompt) so we can begin the tutorial.

I stole this table from [Josh Herr's `R` tutorial](https://github.com/jrherr/quick_basic_R_tutorial/blob/master/R_tutorial.md):

Table 1 - Important R Symbols

| Symbol        | Meaning / Use  |
|:-------------:|:--------------:|
| >      | Prompt for a new command line, you do NOT type this, it appears in the console  |
| +      | Continuation of a previously started command, you do NOT type this      |
| # | Comment in code; R does not "see" this text -- important to explain your computational train of thought          |
| <- or = (= not prefered) |    set a value to a name     |
 

**Note**: Comments in `R` are performed by using `#`. Anything followed by the number sign is ignored by `R`.

For the purposes of this tutorial we will be typing most things at the prompt. However this is annoying and difficult to save. So it is generally best practice to write things in a script editor. The `R` GUIs have some (OSX version with syntax highlighting as does gedit). So does R Studio. If we have time I will show you how to set up syntax highlighting on `nano` as well.
 
##R as a calculator

Let's start by adding `2+2`.
```r
2 + 2
```

This will produce the output:
```
## [1] 4
```

### creating variables in R

Since `R` is a programming language we can easily generate variables to store all sorts of things.

```r
y = 2
```
When you create a variable like this, it does not provide any immediate output.

However when you type `y` and press return:
```r
y
```

You will see the output:
```
## [1] 2
```

Note that the `[1]` is just an index for keeping track where the answer was put.  It actually means that it is the first element in a vector (more on this in a few minutes).

R is case SENSITIVE. That is `y` & `Y` are not the same variable.

```r
y
```

```
## [1] 2
```

```r
Y
```

```
## Error: object 'Y' not found
```

We can generate variables and do computations with them.
```r
x = 3
x + y
```

```
## [1] 5
```

```r
z <- x + y
```


You will notice that sometimes I am using '=' and sometimes '<-'.  These are called *assignment operators*. In most instances they are equivalent.
 The '<-' is preferred in R, and can be used anywhere. You can look at the help file (more on this in a second) to try to parse the difference
```r
`?`("=")
```

We may want to ask whether a variable that we have computed equals something in particular for this we need to use '==' not '=' (one equals is an assignment, two means 'equal to')

```r
x == 3
```

```
## [1] TRUE
```

```r
x == 4
```

```
## [1] FALSE
```

```r
x == y
```

```
## [1] FALSE
```

What happens if we write?

```r
x = y
```

We have now assigned the current value of `y` (2) to `x`. This also shows you that you can overwrite a variable assignment. This is powerful, but also means you need to be very careful that you are actually using the value you think you are.

### standard mathematical operators

Operator `*` for multiplication.

```r
2 * 3
```

```
## [1] 6
```

For division `/`.
```r
6/3
```

```
## [1] 2
```

Operator for exponents `^`. Can also use `**`
```r
3^2
```

```
## [1] 9
```

```r
3**2  # same as above
```

```
## [1] 9
```
You can use `^0.5` or `sqrt()` function for square roots.

```r
9^0.5
```

```
## [1] 3
```

```r
sqrt(9)
```

```
## [1] 3
```

# to raise something to e^some exponent

```r
exp(2)  # this is the performing e^2
```

```
## [1] 7.389
```

For natural log (base e), use `log()`.

```r
log(2.7)
```

```
## [1] 0.9933
```

To raise to an arbitrary base use `log( , base)` like the following:

```r
log(2.7, 10)  # base 10
```

```
## [1] 0.4314
```

You can also use `log10()` or `log2()` for base 10 or base 2.

While these are all standard operators (except `<-`) common to most programming languages, it is a lot to remember. 


# A bit on data structures in `R`

`R` **is vectorized**.  It can do its operations on vectors. Indeed there is no data structure in `R` for scalar values at all. This means when you assign `y <- 2` you are in fact generating a vector of length 1. 

```r
a <- c(2, 6, 4, 5)
b <- c(2, 2, 2, 1)
```

The `c` is short for concatenate you can add the elements of the vectors together. Keep in mind *how* it is doing its addition this way (element-by-element).

```r
a + b
```

```
## [1] 4 8 6 6
```

Or multiply the elements of the vector together (this is **NOT** vector multiplication, but element-by-element multiplication. For vector multiplication (inner & outer products) there are specific operators, i.e. `%*%`.

```r
a * b
```

```
## [1]  4 12  8  5
```

If you want to take vectors `a` and `b` (imagine these are each columns of read counts for RNAseq from different samples) and put them together in a single vector you use the `c()` (concatenate) function. I just realized that I am using `c` as a variable name (for the new vector), and the function is `c()`. This is entirely by accident and they have no relation to one another.

```r
c <- c(a, b)
```

How might you make a vector that repeats vector `a` 3 times?

MANY MANY operations can be vectorized, and R is really good at it!!! Vectors are very useful as way of storing data relevant for all sorts of data analysis.

## GETTING HELP in R

There are a number of places where you can get help with R directly from the console. For instance, what if we want to get help on the function `lm()` which is what is used to fit all sorts of *linear models* like regression, ANOVA, ANCOVA etc.. 

```r
?lm
```
This brings up a description of the function 'lm'


For operators we need to use quotes
```r
?'*' # for help for operators use quotes
```

sometimes you will need to use `help.search('lm')` This brings up all references to the `lm()` function in packages and commands in `R`.  We will talk about packages later.


You can also use `RSiteSearch('lm')`. This is quite a comprehensive search that covers R functions, contributed packages and R-help postings.  It is very useful but uses the web (so you need to be online).

You can also use the html version of help using `help.start()`.

Or if using the GUI, just go to the help menu!
 
## Simple functions in base R

`R` has many functions, and indeed everything (including the operators) are functions (even though the details of the function call are sometimes hidden). As we will see it is very easy to also write new functions.

Here are some examples of the useful functions in base `R`.

You can find out the length of the new vector using `length()`:

```r
length(c) 
```

```
## [1] 8
```


`length() ` is an example of a pre-built function in R. Most things in R revolve around using functions to do something, or extract something. We will write our own simple functions soon.

Here are some more common ones that you may use

```r
mean(c)
```

```
## [1] 3
```

```r
sum(c)
```

```
## [1] 24
```

standard deviation
```r
sd(c)  
```

```
## [1] 1.773
```

variance
```r
var(c)
```

```
## [1] 3.143
```

Correlation coefficient.

```r
cor(a, b) 
```

This gets the Pearson correlation (there are options to change this to other types of correlations, among the arguments for this function.....).

```
## [1] -0.2928
```

Say we want to keep the mean of `c` for later computation we can assign it to a variable

```r
mean_c <- mean(c)
```

We can look at the underlying code of the function (although some times it is buried, in these cases).

We can use a function like `mean()` to add up all of the elements of the vector.

We can also join the two vectors together to make a matrix of numbers.

```r
d <- cbind(a, b)
d
```

```
##      a b
## [1,] 2 2
## [2,] 6 2
## [3,] 4 2
## [4,] 5 1
```

We can double check that we really made a matrix:
```r
is.matrix(d) 
```

This sets up a 'Boolean'. In other words when we ask 'is d a matrix'? it answers TRUE or FALSE.


```
## [1] TRUE
```

```r
mode(d)
```

```
## [1] "numeric"
```

```r
class(d)
```

```
## [1] "matrix"
```
While the mode of d is still numeric (the most basic "atomic" type of the data), the class of the object we have created is a matrix.

Exercise: Make a new vector q that goes a,b,a*b

## Objects in R, classes of objects, mode of objects.

R is an object-oriented language. Everything in R is considered an object. Each object has one or more attributes (which we do not generally need to worry about, but useful for programming.)  Most objects in R have an attribute which is the 'class' of the object, which is what we will usually care about. R has a bunch of useful classes for statistical programming. `bioconductor` also has expanded the classes for objects of use for bioinformatics and genomics. We will see these more this afternoon and in the tutorials next week.


#### Mode of the object. 
The most basic (atomic) feature of an object. NOTE this does not mean the 'mode' of a distribution.

```r
mode(c)
```

```
## [1] "numeric"
```

```r
class(c)  # class of the object
```

```
## [1] "numeric"
```

In this case for the vector `c` the mode and class of `c` is the same. This is not always going to be the case as we see below. 

```r
# typeof(c) # internal representation of type, rarely of interest
```

Let's look at some other objects we have generated.

```
mode(mean_c) 
```

```
## [1] "numeric"
```

```r
class(mean_c)  #
```

```
## [1] "numeric"
```

Despite what we have seen up to now, mode and class are not the same thing. mode tells the basic 'structures' for the objects. integer, numeric (vector), character, logical (TRUE,FALSE) are examples of the atomic structures. 


 There are many different classes of objects each with their own attributes.  The basic ones that we will see are numeric, character, matrix, data.frame, formula, array, list & factor. It is relatively straightforward (but we will not discuss it here) to extend classes as you need.


We can also make vectors that store values that are not numeric.

```r
cities <- c("Okemos", "E.Lansing", "Toronto", "Montreal")
class(cities)
```

```
## [1] "character"
```

```r
mode(cities)
```

```
## [1] "character"
```

Let's use one of the built-in functions that we used before to look at the "length" of `cities`.

```r
length(cities)
```

```
## [1] 4
```

This tells us how many strings we have in the object 'cities' not the length of the string! To get the number of characters in each string we use `nchar()`.


```r
nchar(cities)  # This tells us how many characters we have for each string.
```

```
## [1] 6 9 7 8
```

So if we just do this
```r
q = "okemos"
length(q)
```

```
## [1] 1
```

```r
nchar(q)
```

```
## [1] 6
```

Exercise: How would you compute the total number of characters for all of cities in the object `cities`?


Let's create a second vector storing the names of rivers in each city.


```r
rivers <- c("Red Cedar", "Red Cedar", "Don Valley", "Sainte-Laurent")
cities_rivers <- cbind(cities, rivers)
cities_rivers
```

```
##      cities      rivers          
## [1,] "Okemos"    "Red Cedar"     
## [2,] "E.Lansing" "Red Cedar"     
## [3,] "Toronto"   "Don Valley"    
## [4,] "Montreal"  "Sainte-Laurent"
```

```r
class(cities_rivers)
```

```
## [1] "matrix"
```

```r
mode(cities_rivers)
```

```
## [1] "character"
```

In this above example we have made a matrix, but filled with characters, not numerical values.

####Formula

Another type of object we will need for this workshop (well a variant of it) is called formula. Not surprisingly this is used generally to generate a formula for a statistical model we want to fit.

```r
model_1 <- y ~ x1 + x2 + x1:x2  
model_1
```
This is just the model formula, and we HAVE NOT FIT ANY MODEL YET!!!!!! It just tells us the model we want to fit. That is the object `model_1` has not yet been 'evaluated'.

```
## y ~ x1 + x2 + x1:x2
## <environment: 0x100b504b0>
```

```r
# typeof(model_1)
class(model_1)
```

```
## [1] "formula"
```

```r
terms(model_1)  # also see all.names() and all.vars
```

```
## y ~ x1 + x2 + x1:x2
## attr(,"variables")
## list(y, x1, x2)
## attr(,"factors")
##    x1 x2 x1:x2
## y   0  0     0
## x1  1  0     1
## x2  0  1     1
## attr(,"term.labels")
## [1] "x1"    "x2"    "x1:x2"
## attr(,"order")
## [1] 1 1 2
## attr(,"intercept")
## [1] 1
## attr(,"response")
## [1] 1
## attr(,".Environment")
## <environment: 0x100b504b0>
```

### Factors

Let's make a new vector that will store some read counts of transcript expression values (from a single transcript). The first four are from samples of a "wild type" genotype and the second four samples from a "mutant" genotype. 


```r
counts_transcript_a <- c(250, 157, 155, 300, 125, 100, 153, 175)
```

We may be interested in comparing and contrasting these. So we need to make a variable that stores the information on the two different genotypes. We do this using one of the underlying classes in `R`, called `factors`.

We will make a new variable for the mutant and wild type using the `gl()` (generate levels) function. There are other ways of generating levels for a categorical treatment variable (and R usually can figure this out), but this will work for here.

```r
genotype <- gl(n=2, k=4, labels = c("wild_type", "mutant"))
```

This is our first time using arguments within a function. For the `gl()` the `n=2` means we have two treatment levels. `k=4` means we 4 samples within each. 

One obvious thing we might want to do is make a single object with both variables (genotype, and counts_transcript_a). 

Your first thought (given that we just did this a few minutes ago) might be to make a matrix from this.

```r
expression_data <- matrix(counts_transcript_a, genotype)
```

But as you see you get an error message

```
## Error in matrix(counts_transcript_a, genotype) : 
##  non-numeric matrix extent
```

Why? Because for an object of class `matrix` all of the atomic types in it need to be the same. Basically we can not combine numeric and factors. 

The way we do this is to instead create a `data.frame`. A `data.frame` is a particular representation of another type of object (`list`) that allows for the storage of heterogeneous data types. Effectively you can think about it like a spreadsheet for purposes of analysis (for now anyways.)

```r
expression_data <- data.frame(counts_transcript_a, genotype)
expression_data
```

When you import most of your data, and do most analyses, you will more often than not be using a data frame.


## Workspaces, and objects in them

R stores variables, datafiles, functions, vectors, etc in what is called the Workspace. This contains all of the items that you can access directly within your R session.  You can list all of the objects in your
# workspace using:

```r
ls()
```

```
##  [1] "a"             "b"             "c"             "cities"       
##  [5] "cities_rivers" "d"             "mean_c"        "model_1"      
##  [9] "q"             "rivers"        "x"             "y"            
## [13] "z"
```

If you want to remove a particular variable (say x) use the rm() function

```r
rm(x)
```
you could remove multiple objects

```r
rm(x, y, z)
```

```
## Warning: object 'x' not found
```
If you want to remove all of the objects in your workspace
`rm(list = ls())`. We will learn what this means later, but basically we are making a list that contains all of the objects found by performing `ls()`, and then removing everything in that list.

###Saving the workspace.

Some people like to save their workspaces, not only because it contains all of the commands they have written, but also all of the objects they have created during that session.  I personally do not do this unless I have created objects that have taken a long time to compute. Instead I just save the scripts I write.

However if you write your commands directly at the console (like we have been doing, but really you should not do) without a script editor, you should save your workspaces.

```r
save.image('file_name')
```

Which will save it to your current working directory (you can find that using `getwd()`).

If you want to load it again
```r
load('file_name.RData')
```

You will need to have the correct working directory set, which I will show you how to do shortly.

## SCRIPT!

Writing everything at the console can be a bit annoying, so we will use a script editor.

In Mac OS X I personally find the built-in script editor useful You can highlight the text in the script editor and press command (apple) + return to send it to the R console. Or place the cursor at the end of the line
that you want to submit to R with command+ return. It also provides syntax highlighting, and shows the syntax & options for functions.

However, for those of you are under the spell of Bill Gates....... While the basic script editor for windows does not have much functionality, many people have written excellent script editors. The base script editor in windows will submit a line with ctrl-R(???). There are many windows script editors with syntax highlighting (such as Tinn-R).

For a list of some http://www.sciviews.org/_rgui/ In general we will save R scripts with the extension `.R`

Also there is R-studio.


let's type something into our new script

```r
x <- c(3, 6, 6, 7)
```

now highlight that line and press ctrl+r (windows), or apple key + return (mac). This should send the highlighted portion to R.


```r
x <- c(2, 2, 2, 2)
y <- c(3, 3, 3, 3)
z <- cbind(x, y)
z
```

```
##      x y
## [1,] 2 3
## [2,] 2 3
## [3,] 2 3
## [4,] 2 3
```

We have also just used the function, `cbind`, which stands for column bind. This will create a new object stitching them together as columns.


## Writing our own functions in R

We have now used a few built-in functions in R (there are many).  Anything where you use `()` is a function in `R`. Like I mentioned, pretty much everything you do in `R` is actually a call to a function.

However, we will often want to compute something for which there is no pre-built function. Thankfully it is very easy to write our own functions in R. You should
definitely get in the habit of doing so.

Functions have the following format:

```
aFunction <- function(input variable 1, input variable 2, argument, etc...) {expressions to calculate}
```
This is abstract so let me give you a real example. For our read counts, we want to compute the standard error of the mean (a measure of sampling uncertainty), which is ~equal to the sd/sqrt(sample size). How might we do it?

We want to compute it for the numeric vector of read counts.

```r
counts_transcript_a
```

We could do it by hand

```r
sd_a <- sd(counts_transcript_a)
sample_a <- length(counts_transcript_a)
sd_a/sqrt(sample_a)
```

```
## [1] 23.35857
```
(notice that the last value was printed, as we did not store it in a variable).

Or we could do it in one line
```r
sd(a)/sqrt(length(a))  # notice the function within a function
```

```
## [1] 23.35857
```

But we can also do it so that we can use any vector input we wanted by writing a function.

```r
StdErr <- function(vector) {
    sd(vector)/sqrt(length(vector))
}
```

Now type `StdErr`
```r
StdErr
```

```
## function(vector) {
##     sd(vector)/sqrt(length(vector))
## }
## <environment: 0x100b504b0>
```

This just repeats the function. If you want to edit the function just type `edit(StdErr)`.

Let's use our new function
```r
StdErr(counts_transcript_a)
```

```
## [1] 23.35857
```
But now we can use this for any set of samples that we need to calculate the SEM. In this case transcript 'b'.

```r
counts_transcript_b <- c(75, 85, 82, 79, 77, 83, 96, 62)
StdErr(counts_transcript_b)
```

```
## [1] 3.414452
```

Exercise: Write your own function to do something simple, like calculate the co-efficient of variation (CV) which is the sd/mean.

It takes some practice but learning to write small discrete functions can be
extremely helpful for R.

One thing to keep in mind, is that it is very easy to call one function from within another.  It is generally considered good practice to write functions that do one thing, and one thing only. It is way easier to find problems (debug).

One of the great things that we can (and need to do) often is to compute the mean (or SEM, etc..) for each transcript by sample. We can do this for the data set we made

```r
expression_data
```

I will not explain it at the moment, but we can use one of the apply like functions to compute the mean and SEM for each genotype (works the same whether it is 2 or 2000).

```r
with(expression_data, tapply(X=counts_transcript_a, INDEX=genotype, FUN=mean))
```

```
## wild_type    mutant 
##  215.50    138.25 
```

And then for the SEM

```r
with(expression_data, tapply(X=counts_transcript_a, INDEX=genotype, FUN=StdErr))
```

```
# wild_type    mutant 
#  35.83876  16.34715
``` 
The `with()` just makes it a bit easier to utilize the variables within the `expression_data` object. There are a number of other ways of doing it, primarily using the `$` (extract) operator (for lists including data frames). You will also see people use the `attach()`. Avoid using attach at all costs.


## Using source() to load your functions

One of the great things about writing simple functions, is that once you have them working, you can keep using them over and over.  However, it is generally a pain(and bad practice) to have to include the text of the function in every script you write (what if there is a bug...).  Instead, R has a function source() which allows you to 'load' a script that contains functions you have written (and other options you may want), so that you can use them.

We will likely see `source()` in later tutorials, so watch for it.

## Regular Sequences

Sometimes we want regular sequences or to create objects of repeated numbers or characters. R makes this easy.

 If you want to create regular sequences of integers by units of 1

```r
one_to_20 <- 1:20
one_to_20
```

```
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
```

```r
twenty_to_1 <- 20:1
twenty_to_1
```

```
##  [1] 20 19 18 17 16 15 14 13 12 11 10  9  8  7  6  5  4  3  2  1
```

for other more complicated sequences, use the `seq()` function

```r
seq1 <- seq(from = 1, to = 20, by = 0.5)
seq1
```

```
##  [1]  1.0  1.5  2.0  2.5  3.0  3.5  4.0  4.5  5.0  5.5  6.0  6.5  7.0  7.5
## [15]  8.0  8.5  9.0  9.5 10.0 10.5 11.0 11.5 12.0 12.5 13.0 13.5 14.0 14.5
## [29] 15.0 15.5 16.0 16.5 17.0 17.5 18.0 18.5 19.0 19.5 20.0
```

or

```r
seq1 <- seq(1, 20, 0.5)
seq1
```

```
##  [1]  1.0  1.5  2.0  2.5  3.0  3.5  4.0  4.5  5.0  5.5  6.0  6.5  7.0  7.5
## [15]  8.0  8.5  9.0  9.5 10.0 10.5 11.0 11.5 12.0 12.5 13.0 13.5 14.0 14.5
## [29] 15.0 15.5 16.0 16.5 17.0 17.5 18.0 18.5 19.0 19.5 20.0
```

This shows that for default options (in the correct order) you do not need to specify things like 'from' or 'by'

Exercise: Make a sequence from -10 to 10 by units of 2

What if you want to repeat a number or character a set number of times?

```r
many_2 <- rep(2, times = 20)
```

Works for characters as well

```r
many_a <- rep("a", times = 10)
```

We can even use this to combine vectors

```r
seq_rep <- rep(20:1, times = 2)
seq_rep
```

```
##  [1] 20 19 18 17 16 15 14 13 12 11 10  9  8  7  6  5  4  3  2  1 20 19 18
## [24] 17 16 15 14 13 12 11 10  9  8  7  6  5  4  3  2  1
```

What if you wanted to repeat a sequence of numbers (1,2,3) 3 times?
```r
rep_3_times <- rep(c(1, 2, 3), times = 3)
# or
rep(1:3, times = 3)
```

```
## [1] 1 2 3 1 2 3 1 2 3
```

What if we wanted to perform this to create a matrix

```r
matrix(rep(20:1, 4), 20, 4)
```

```
##       [,1] [,2] [,3] [,4]
##  [1,]   20   20   20   20
##  [2,]   19   19   19   19
##  [3,]   18   18   18   18
##  [4,]   17   17   17   17
##  [5,]   16   16   16   16
##  [6,]   15   15   15   15
##  [7,]   14   14   14   14
##  [8,]   13   13   13   13
##  [9,]   12   12   12   12
## [10,]   11   11   11   11
## [11,]   10   10   10   10
## [12,]    9    9    9    9
## [13,]    8    8    8    8
## [14,]    7    7    7    7
## [15,]    6    6    6    6
## [16,]    5    5    5    5
## [17,]    4    4    4    4
## [18,]    3    3    3    3
## [19,]    2    2    2    2
## [20,]    1    1    1    1
```


## Indexing, extracting values and subsetting from the objects we have created

Often we will want to extract certain elements from a vector, list or matrix. Sometimes this will be a single number, sometimes a whole row or column.

We index in R using `[ ]` (square brackets).

**NOTE** `R` indexes starting with 1, not 0!

```r
a <- 1:20
b <- 5 * a
a
```

```
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
```

```r
b
```

```
##  [1]   5  10  15  20  25  30  35  40  45  50  55  60  65  70  75  80  85
## [18]  90  95 100
```

```r
length(a)
```

```
## [1] 20
```

```r
length(b)
```

```
## [1] 20
```

If we want to extract the 5th element from `a`.
```r
a[5]
```

```
## [1] 5
```

If we want to extract the 5th and 7th element from 'b'

```r
b[c(5, 7)]
```

```
## [1] 25 35
```


IF we want to extract the fifth through 10th element from 'b'
```r
b[5:10]
```

```
## [1] 25 30 35 40 45 50
```

How about if we want all but the 20th element of 'a'?
```r
a[-20]
```

```
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
```

Indexing can also be used when we want all elements greater than (less than etc...) a certain value. Under the hood this is generating a logical/boolean (T vs. F).

```r
b[b > 20]
```

```
##  [1]  25  30  35  40  45  50  55  60  65  70  75  80  85  90  95 100
```

Or between certain numbers

```r
b[b > 20 & b < 80]
```

```
##  [1] 25 30 35 40 45 50 55 60 65 70 75
```

Exercise: generate a vector with 20 elements create a 'sub' vector that has elements 1:5, 16:20 create a 'sub' vector with odd elements 1,3,5,...,19.


### Indexing for matrices

```r
c <- a + b
q_matrix <- cbind(a, b, c)
q_matrix
```

`cbind()` 'binds' column vectors together into a matrix (also see `rbind()`).

```
##        a   b   c
##  [1,]  1   5   6
##  [2,]  2  10  12
##  [3,]  3  15  18
##  [4,]  4  20  24
##  [5,]  5  25  30
##  [6,]  6  30  36
##  [7,]  7  35  42
##  [8,]  8  40  48
##  [9,]  9  45  54
## [10,] 10  50  60
## [11,] 11  55  66
## [12,] 12  60  72
## [13,] 13  65  78
## [14,] 14  70  84
## [15,] 15  75  90
## [16,] 16  80  96
## [17,] 17  85 102
## [18,] 18  90 108
## [19,] 19  95 114
## [20,] 20 100 120
```

(What happens if we ask for the length of q.matrix?...)

```r
length(q.matrix)
````

We can instead ask for number of rows or columns

```r
nrow(q_matrix)
```

```
## [1] 20
```

```r
ncol(q_matrix)
```

```
## [1] 3
```

Or just use `dim(q_matrix)` (dimensions) to get both # rows and # columns. R always specifies in row by column format.

```
## [1] 20  3
```

Say we want to extract the element from the 3rd row of the second column (b)?


```r
q_matrix[3, 2]
```

```
##  b 
## 15
```

How about if we want to extract the entire third row?

```r
q_matrix[3, ]
```

```
##  a  b  c 
##  3 15 18
```


We can also pull things out by name

```r
q_matrix[, "c"]  
```

```
##  [1]   6  12  18  24  30  36  42  48  54  60  66  72  78  84  90  96 102
## [18] 108 114 120
```
This is an example of indexing via 'key' instead of numerical order.

### Accessing values in objects


The at `@` is used to extract the contents of a slot in an object.. We will not use it much for this class, but it is essential for object oriented programming in R (S4) objects. objectName@slotName.

More often we will use the dollar sign `$`, which is used to extract elements of an object of class list (including data frames).. We will use this a lot to extract information from objects (such as information from our models, like co-efficients) object.name$element.name.

For more information `?'$'`.

## Where to go from here?
There are a huge number of resources for R. Everyone has favorite tutorials and books. Here are but a few.

[I have a few screencasts](http://beaconcourse.pbworks.com/w/page/62520939/CSE845_2013_Screencasts) that you can access. I also have a number of [tutorials](http://beaconcourse.pbworks.com/w/browse/#view=ViewFolder&param=2014_August). I have way more R resources for a graduate class I teach in computational statistical approaches which I am happy to share as well.

The R site also has access to [numerous tutorials and books](http://cran.r-project.org/doc/manuals/) or [other documents](http://www.r-project.org/other-docs.html).

For more advanced programming check out [Hadley Wickham's online book](http://adv-r.had.co.nz/).

Here is a reasonably decent [R wikibook](http://en.wikibooks.org/wiki/R_Programming/Advanced_programming)

I really like the book [art of R programming](http://www.nostarch.com/artofr.htm).


##A few advanced topics... For your own amusement (not nescessary for now, but helps for more advanced R programming).

### Setting attributes of objects.

Objects have attributes. The one we have thought about most is the class  of the object, which tells us (and R) how to think about the object, and how it can be used or manipulated (methods). We have also looked at dim() which is another attribute Here is a list of common ones: class, comment, dim, dimnames, names, row.names and tsp.

We can set attributes of objects in easy ways like

```r
x <- 4:6
names(x) <- c("observation_1", "observation_2", "observation_3")
x
```

```
## observation_1 observation_2 observation_3 
##             4             5             6
```

You can see the attributes in a bunch of ways
```r
str(x)
```

```
##  Named int [1:3] 4 5 6
##  - attr(*, "names")= chr [1:3] "observation_1" "observation_2" "observation_3"
```

```r
attributes(x)
```

```
## $names
## [1] "observation_1" "observation_2" "observation_3"
```

Same as above, but we will be able to use this to set attributes of the object x as well
```r
attr(x, "names") 
```

```
## [1] "observation_1" "observation_2" "observation_3"
```

```r
y <- cbind(1:5, 11:15)
attributes(y)
```

```
## $dim
## [1] 5 2
```

```r
colnames(y) <- c("vec1", "vec2")
comment(y) <- c("the first column is pretend data", "the second column is yet more pretend data ")
str(y)
```

```
##  int [1:5, 1:2] 1 2 3 4 5 11 12 13 14 15
##  - attr(*, "dimnames")=List of 2
##   ..$ : NULL
##   ..$ : chr [1:2] "vec1" "vec2"
##  - attr(*, "comment")= chr [1:2] "the first column is pretend data" "the second column is yet more pretend data "
```

```r
attributes(y)
```

```
## $dim
## [1] 5 2
## 
## $dimnames
## $dimnames[[1]]
## NULL
## 
## $dimnames[[2]]
## [1] "vec1" "vec2"
## 
## 
## $comment
## [1] "the first column is pretend data"           
## [2] "the second column is yet more pretend data "
```

```r


### Generic functions and methods 
Calling a function like summary() will do very different things for different object classes.  We will use this call a lot for data frames and output from statistical models, etc..

```r
summary(x)  # numeric vector
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     4.0     4.5     5.0     5.0     5.5     6.0
```

```r
summary(string_1)  # character string
```

```
##    Length     Class      Mode 
##         1 character character
```

The call to `summary()` is generic, which first looks at the class of the object, and then uses a class specific method to generate a summary of x.

```r
summary(x)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     4.0     4.5     5.0     5.0     5.5     6.0
```

```r
summary.default(x)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     4.0     4.5     5.0     5.0     5.5     6.0
```

but..

```r
summary.lm(x)  # Since this was looking for an object of class lm
```

```
## Error: $ operator is invalid for atomic vectors
```

To see all of the methods used when you call the generic `summary()` for S3 classes.

```r
methods(summary) 
```

```
##  [1] summary.aov             summary.aovlist        
##  [3] summary.aspell*         summary.connection     
##  [5] summary.data.frame      summary.Date           
##  [7] summary.default         summary.ecdf*          
##  [9] summary.factor          summary.glm            
## [11] summary.infl            summary.lm             
## [13] summary.loess*          summary.manova         
## [15] summary.matrix          summary.mlm            
## [17] summary.nls*            summary.packageStatus* 
## [19] summary.PDF_Dictionary* summary.PDF_Stream*    
## [21] summary.POSIXct         summary.POSIXlt        
## [23] summary.ppr*            summary.prcomp*        
## [25] summary.princomp*       summary.proc_time      
## [27] summary.srcfile         summary.srcref         
## [29] summary.stepfun         summary.stl*           
## [31] summary.table           summary.tukeysmooth*   
## 
##    Non-visible functions are asterisked
```

## Syntax style guide
Generally it is advisable to use a consistent way of scripting. For any given programming language there are syntax style guide. The [Style guide for my class](https://www.msu.edu/~idworkin/ZOL851_style_guide.html). You can also check out the [R style guide from Google](http://google-styleguide.googlecode.com/svn/trunk/google-r-style.html).


## Random bits
Note about using q() on the Mac R GUI in v2.11.+ The programming team decided the default behaviour was potentially 'dangerous', and people may lose their files, so they have changed it to command + q to quit instead. If you are an old-fogey like me and like to use q(), you have a couple of options. base::q() # This will work, but it is annoying.

you can set your `.Rprofile` to have the following line.
```r
options(RGUI.base.quit=T) 
```
and the next time you run R the old `q()` will work.

If you do not know how to create or edit .Rprofile, come speak with me...

## session info
The R session information (including the OS info, R version and all
packages used):

```r
sessionInfo()
```

```
## R version 3.0.1 (2013-05-16)
## Platform: x86_64-apple-darwin10.8.0 (64-bit)
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] knitr_1.5
## 
## loaded via a namespace (and not attached):
## [1] evaluate_0.5.1 formatR_0.10   stringr_0.6.2  tools_3.0.1
```

```r
Sys.time()
```

```
## [1] "2014-08-07 10:39:49 EDT"
```


### A few thoughts about indexing for those used to Python (If you don't, ignore this)


# R indexing begins at 1 (not 0 like Python) Negative values of indexes in R
# mean something very different. for instance
a[-1]  # this removes the first element of a, and prints out all of the remaining elements.
```

```
##  [1]  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
```

```r

# As far as I know all classes of objects are mutable, which means you can
# write over the name of the objects, values within the objects, and
# slots....

# Indexing on a character string does not work in R
string_1 <- "hello world"
string_1[1]
```

```
## [1] "hello world"
```

```r

# instead you need to use the substr() function
substr(x = string_1, start = 1, stop = 1)
```

```
## [1] "h"
```

```r

# similarly
length(string_1)  # this gives an output of 1
```

```
## [1] 1
```

```r
nchar(string_1)  # this gives the 11 characters
```

```
## [1] 11
```

# TOC
# Section 1: What is R;  R at the console; quiting R
# Section 2: R basics; R as a calculator; assigning variables; vectorized computation in R
# Section 3: pre-built functions in R
# Section 4: Objects, classes, modes - Note: should I add attributes?
# Section 5: The R workspace; listing objects, removing objects (should I add attach and detach?)
# Section 6: Getting Help in R
# Section 7: Using A script editor for R
# Section 8: Writing simple functions in R
# Section 8b: Using source() to call a set of functions 
# Section 9: Regular sequences in R
# Section 10: Extracting (and replacing), indexing & subsetting (using the index). Can also be used for sorting.

#### Advanced stuff to learn on your own...
# .....  setting attributes of objects.... (names, class, dim )
# .....  environments   (see ?environments)

