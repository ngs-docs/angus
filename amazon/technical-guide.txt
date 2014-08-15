Technical guide to the ANGUS course
===================================

Packages we install
-------------------

Install::

   apt-get update
   apt-get -y install screen git curl gcc make g++ python-dev unzip \
           default-jre pkg-config libncurses5-dev r-base-core \
           r-cran-gplots python-matplotlib sysstat samtools python-pip \
           ipython-notebook

Creating your own AMI
---------------------

Boot up the machine you want (Ubuntu 14.04 LTS PV is what we used in
2014), then run the apt-get commands above.  Shut it down, and then
follow these instructions:

  http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/creating-an-ami-ebs.html

I usually name them something like '2014-08-03 ubuntu 14.04 angus',
with 'v1', 'v2', 'v3' on the end of 03 (e.g. 2014-08-03.v1) just so I
can sort them more easily.

At the end of this, you should be able to boot it up and run the
apt-get commands, above, and have them complete instantly (beacuse
it's up to date)

Be sure to make it public! (Actions... Modify Image permissions under
AMI menu).

