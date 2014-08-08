
### reading in data, combining, parsing, etc...

# we read in data files using a family of 'read.' functions...
# in this case we have our data stored as comma separated values (.csv) so we'll use 'read.csv'

genes_fpkm  <- read.csv( "~/Desktop/genes_fpkm_table.csv" )

#...will spit the contents of the file to your console. Not super-effective :-(
# a better idea would eb to make an R object out of it

genes_fpkm <- read.csv( "~/Desktop/genes_fpkm_table.csv" )

# Alternatively, we can point R at the directory that we want use like this:

setwd( "~/Desktop" )

# which makes reading in data (and writing things out; see below) easier to write

genes_attr <- read.csv( "genes.attr_table.csv" )

#Looking at your data
# R is usually smart enough to make something sensible out of your data, in this case you ought to have a data.frame
class( genes_fpkm )

# which we can learn more about using 'str'
# note that we can see how R has guessed (hopefully) that the 1st column was labels and the following columns are numeric


str( genes_fpkm )

# And let's look at a summary of the data

summary(genes_fpkm)

head(genes_fpkm)

# Any missing data
sum(is.na(genes_fpkm))


# And the gene attributes.
# note that this data is all non-numeric: we have meta-data to apply to the sequencing numbers from the first dataset

class( genes_attr )
str( genes_attr )

# so, we have 2 related data.frames. We could keep them separate and R would be happy to work with them, but it is ugly
# e.g. we could get combined data & meta-data for a row (row 14 here) like this:

cbind( genes_attr[ genes_attr$tracking_id == genes_fpkm$tracking_id[ 14 ], ], genes_fpkm[ 14, c(1, 4:5) ] )

# ...but that's no fun. It's often easier to combine the two datasets in R's memory. There are a few ways to do this.

# Situation 1: "all my data are beautifully arranged"
# there are the same no. rows

dim( genes_attr )
dim( genes_fpkm )

# ...and the tracking_ids are in the same order

length( intersect( genes_attr$tracking_id, genes_fpkm$tracking_id ))  # same IDs

#This is effectively doing this...
(genes_attr$tracking_id == genes_fpkm$tracking_id)    # same order

# we can simply glom the two datasets together

genes_data <- cbind( genes_attr, genes_fpkm )
str( genes_data )

genes_data <- cbind( genes_attr, genes_fpkm[, 2:ncol(genes_fpkm) ] ) # if you don't want to double up the IDs

# Situation 2: "my data is <100% perfect" (i.e. real-world)

# lets reorder 1 dataset
genes_attr_2 <- genes_attr[ order( genes_attr$locus ), ]

# now we will use the binary operator %in$ for the match() to match columns.
length( intersect( genes_attr_2$tracking_id, genes_fpkm$tracking_id ))    # same IDs
(genes_attr_2$tracking_id == genes_fpkm$tracking_id)        # different order!?!?!?


# in this case we want to combine rows *only* where the tracking_ids are confirmed to be the same
# as always with R, there are a few ways to achieve in the given goal

genes_data_2 <- data.frame( genes_attr_2,
                            genes_fpkm[ genes_fpkm$tracking_id %in% genes_attr_2$tracking_id , 2:ncol(genes_fpkm) ] )

# or

genes_data_3 <- merge( genes_attr_2, genes_fpkm, by.x="tracking_id", by.y="tracking_id" )


str( genes_data_2 )
str( genes_data_3 )


########

### Doing interesting things with our data!

# look at the str (or names) of the new data.frame

names( genes_data_2 )

# you'll notice that we have biological replicates, e.g. Hyb_vg_0 & Hyb_vg_1
# we might be interested in how repeateble our measurements are

plot( genes_data_2$Hyb_vg_0, genes_data_2$Hyb_vg_1 )
plot( log2(genes_data_2$Hyb_vg_0), log2(genes_data_2$Hyb_vg_1) )

cor( genes_data_2$Hyb_vg_0, genes_data_2$Hyb_vg_1 )

# we seem to be in good shape for Hyb_vg, how about we get fancy

plot( genes_data_2$Hyb_vg_0, genes_data_2$Hyb_vg_1, xlim=c( 0, 60000 ), ylim=c( 0, 60000 ), 
    xlab="rep_0", ylab="rep_1", main= paste( "correlation=", 
    eval(round( cor( genes_data_2$Hyb_vg_0, genes_data_2$Hyb_vg_1 ), 2 )) ) )

abline( a=0, b=1, col="red" )

# how about the other genotypes?

# let's make some easy-to-use data.frames

genes_data_2_rep_0 <- genes_data_2[, seq( 5, 15, 2) ]
genes_data_2_rep_1 <- genes_data_2[, seq( 6, 16, 2) ]

names( genes_data_2_rep_1 )
names( genes_data_2_rep_0 )

# as an excercise; why not use the help functions to work out what I'm doing here?

par( mfrow=c( 2, 3 ) )

for ( raptor in 1:6 ) {
  plot( log10( genes_data_2_rep_0[, raptor] ), log10( genes_data_2_rep_1[, raptor] ),
        main= paste( substr( names( genes_data_2_rep_0 )[ raptor ], 1, 6 ) ),
        xlab= paste( names( genes_data_2_rep_0 )[ raptor ] ),
        ylab= paste( names( genes_data_2_rep_1 )[ raptor ] ) )
  abline(  a=0, b=1, col="red" )
  text( -2, max( log10( genes_data_2_rep_1[, raptor]) - 1 ),
        paste( "cor=", eval(round( cor( genes_data_2_rep_0[, raptor], genes_data_2_rep_1[, raptor] ), 2 )) ),
        cex= 1.5, col= "red")
}

par( mfrow= c(1, 1) )

# maybe we are interested in differences among chromosomes...
# this information is currently part of the 'locus' factor

str( genes_data_2 )
genes_data_2$locus[ 1:5 ]

# lets make a new classifying factor to identify data by chromosome

genes_data_2$chromo <- factor( sub( ':[[:digit:]]+-[[:digit:]]+', "", genes_data_2$locus ))

# ...with which we can interrogate our data

genes_data_2$all_fpkm <- rowMeans( genes_data_2[, 5:16 ] )

plot( genes_data_2$chromo, genes_data_2$all_fpkm  )

# are these chromosomes really different?

chromo_mod <- lm( all_fpkm ~ chromo, data= genes_data_2 )
summary( chromo_mod )

# there appears to be some evidence that the mitochondria and some of the 'Het' reads are a bit weird compared to the others...
# let's say that we're only interested in chromosomes X, 2 & 3 (which I'll call the shorthand genome)

levels( genes_data_2$chromo )

shorthand_chromos <- factor( c( "X", "2L", "2R", "3L", "3R" ))

genes_data_shorthand <- genes_data_2[ genes_data_2$chromo %in% shorthand_chromos , ]

dim( genes_data_2 )
dim( genes_data_shorthand )

# what percentage of rows did we lose?

nrow( genes_data_shorthand ) / nrow( genes_data_2 ) * 100


