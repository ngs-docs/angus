
echo "Enter file list"
read filelist
echo "Enter search term"
read sequence

 for filenames in `${filelist}`   #Matt thought calling the thing bananas was funny
        do grep -i -B 1 ${sequence} ${filenames} | #search for  ThgI
        grep ">" | #search for header lines
        cut -f 1 -d " " | #Only keep sequence name
        cut -c 2- # get rid of ">"
 done