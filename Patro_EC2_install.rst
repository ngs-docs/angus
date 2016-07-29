First, install the "build tools" (compilers etc. that may be needed)

'> sudo apt-get update`

`> sudo apt-get install build-essential`

To make use of linuxbrew, we'll need Ruby:

`> sudo apt-get install git ruby`

Now, we can get linuxbrew with the following command:

`> ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"`

when prompted, hit `RETURN`.  To make the software we'll install via linuxbrew accessable we have 
to place the directory where linuxbrew installs programs into our "path".  The following sequence of 
commands will do this:

`> echo 'export PATH="/home/ec2-user/.linuxbrew/bin:$PATH"' >>~/.bash_profile`

`> echo 'export MANPATH="/home/ec2-user/.linuxbrew/share/man:$MANPATH"' >>~/.bash_profile`

`> echo 'export INFOPATH="/home/ec2-user/.linuxbrew/share/info:$INFOPATH"' >>~/.bash_profile`

The programs we're interested in installing are part of hombrew-science.  We can "tap" the science keg (;P) as follows:

`> brew tap homebrew/science`

Now, we can install the STAR aligner like so:

`> brew install homebrew/science/rna-star`

Since the latest (pre-release) salmon is not yet a binary available in linuxbrew, we'll grab a pre-compiled binary directly.
We can download it using `wget` like so:

`> wget --no-check-certificate 'https://drive.google.com/uc?export=download&id=0B3iS9-xjPftjaFQwYUlvQnN0UFU' -O Salmon-v0.7.0.tgz`

and we can untar and unzip the resulting file with the following command:

`> tar xzf Salmon-v0.7.0.tgz`

Finally, so that we can simply type `salmon` to execute salmon, we'll add the appropriate directory to our path variable again.

`> echo 'export PATH="SalmonBeta-0.7.0-pre-july27_CentOS5/bin:$PATH"' >>~/.bash_profile`

The path will be automatically set when you login, but we want the changes to take effect now, so we must "source" the 
file containing all of the commands that we created:

`> source ~/.bash_profile`
