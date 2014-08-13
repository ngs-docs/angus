===========================================================
Logging into your new instance "in the cloud" (Mac version)
===========================================================

OK, so you've created a running computer.  How do you get to it?

The main thing you'll need is the network name of your new computer.
To retrieve this, go to the instance view and click on the instance,
and find the "Public DNS".  This is the public name of your computer
on the Internet.

Copy this name, and connect to that computer with ssh under the username
'root', as follows.

First, find your private key file; it's the .pem file you downloaded
when starting up your EC2 instance.  It should be in your Downloads
folder.  Move it onto your desktop and rename it to 'amazon.pem'.

Next, start Terminal (in Applications... Utilities...) and type::

  chmod og-rwx ~/Desktop/amazon.pem

to set the permissions on the private key file to "closed to all evildoers".

Then type::

  ssh -i ~/Desktop/amazon.pem ubuntu@ec2-???-???-???-???.compute-1.amazonaws.com

(but you have to replace the stuff after the '@' sign with the name of the host).

Here, you're logging in as user 'ubuntu' to the machine
'ec2-174-129-122-189.compute-1.amazonaws.com' using the authentication
key located in 'amazon.pem' on your Desktop.

You should now see a text line that starts with something like
``ubuntu@ip-10-235-34-223:~$``.  You're in!  Now type::

   sudo bash
   cd /root

to switch into superuser mode (see: http://xkcd.com/149/) and go to your
home directory.

This is where the rest of the tutorials will start!

If you have Dropbox, you should now visit :doc:`installing-dropbox`.

You might also want to read about :doc:`terminating-your-instance`.

To log out, type::

   exit
   logout

or just close the window.
