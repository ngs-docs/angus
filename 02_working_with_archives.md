
###Working with Compressed files

As previously mentioned, genomics data files tend to be *large*. Since larger files are slower and more costly to move around, you will often encounter files that have been *compressed* to save time/space/money. The two most commonly encountered types of compressed files are Zip archives (e.g. `filename.zip`), Gzip archives (e.g. `filename.gz`) and Tarballs (e.g. `filename.tar` or `filename.tar.gz`).

Once you've convinced yourself that the file you have is the file that you *ought* to have, the next thing that you'll want to do is unzip it (a.k.a. uncompress or decompress or extract). You can unzip your .zip archive using the `unzip` program:

```bash
unzip <filename.zip>
```

If you don't want to extract everything, but rather check the contents, you can view what a zip contains using the `-l` flag ('list'):

```bash
unzip -l <filename.zip>
```

When you want to go in the other direction and make your own archive the command is simply `zip`. It works like this:

```bash
zip <mynewarchive.zip> <myfirstfile.txt> <mysecondfile.sam>
```

Note that you can also use the `-r` flag (recursive) to zip up a folder and all its contents, including subfolders like so:

```bash
zip -r <myproject.zip> myproject/
```

If you have been sent a big bundle of data as a gzip archive, then happily the same procedure applies for viewing and extracting as with zip archives, but with the `gunzip` program:

```bash
gunzip -l <bundle.gz>
gunzip <bundle.gz>
```

Things are slightly different (read 'complex') if you encounter a tarball: `thisfile.tar` or `thatfile.tar.gz` or `tacofile.tgz`.

![]( http://imgs.xkcd.com/comics/tar.png "I'm so sorry..." )

You can view the contents of tarballs using the `tar` program:

```bash
tar -tf <thisfile.tar>
tar -ztvf <thatfile.tar.gz>
tar -ztvf <tacofile.tgz>
```

...and extract them like this:

```bash
tar -xf <thisfile.tar>
tar -zxvf <thatfile.tar.gz>
tar -zxvf <tacofile.tgz>
```

Other types of compressed files and archives do exist, but these are the most common. 
