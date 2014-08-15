Installing Dropbox on your EC2 machine
======================================

**IMPORTANT: Dropbox will sync everything you have to your EC2 machine, so
if you are already using Dropbox for a lot of stuff, you might want to 
create a separate Dropbox account just for the course.**

Start at the login prompt on your EC2 machine::

  sudo bash
  cd /root

Then, grab the latest dropbox installation package for Linux::

   wget -O dropbox.tar.gz "http://www.dropbox.com/download/?plat=lnx.x86_64"

Unpack it::

   tar -xvzf dropbox.tar.gz

Make the Dropbox directory on /mnt and link it in::

   mkdir /mnt/Dropbox
   ln -fs /mnt/Dropbox /root

and then run it, configuring it to put stuff in /mnt::

   HOME=/mnt /root/.dropbox-dist/dropboxd &

When you get a message saying "this client is not linked to any account",
copy/paste the URL into browser and go log in.  Voila!

Your directory '/root/Dropbox', or, equivalently, '/mnt/Dropbox', is now
linked to your Dropbox account!
