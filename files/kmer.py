#this will be analyze.py
# to run, type:
# python analyze.py

# this script relies on the sys package, must import
import sys 

#data for this script:
#http://athyra.idyll.org/~t/transfer/aa.txt
#http://athyra.idyll.org/~t/transfer/bb.txt
#http://athyra.idyll.org/~t/transfer/cc.txt
#http://athyra.idyll.org/~t/transfer/dd.txt


#print "hello world!"

# set open list
aalist = []

# opens the file, supplied as the first argument after script filename on command line when this script is run,
# reads each line
# adds each line of data to list after stripping white space at the end of each line
for line in open (sys.argv[1]):
    line = line.strip()
    aalist.append(line)

# these are commented out as tests, may want to look at later?
#start = aalist[0]
#print start

# set open dictionary
aadict = {}

#kmers of 4 or 10
# this value of k is variable, change this later to see what happens with large or small kmer
k=10

# goes through each data entry in aalist

for read in aalist:
# splits each data entry into words of length k
    for j in range(0,len(read) - k + 1):
# adds each kmer into dictionary with count
# checks to see if it already exists in dictionary 
# if so, adds count value
        kmer = read[j:j+k]
        aadict[kmer] = aadict.get(kmer,0)+1
        #print j, kmer, aadict[kmer]
        
# prints all kmers and abundance counts
for kmer in aadict:
    print aadict[kmer], kmer
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

