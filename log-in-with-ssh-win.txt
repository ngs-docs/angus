===============================================================
Logging into your new instance "in the cloud" (Windows version)
===============================================================

Download Putty and Puttygen from here: http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html

Generate a ppk file from your pem file
======================================

(You only need to do this once!)

Open puttygen; select "Load".

.. image:: images/win-puttygen.png
   :width: 90%

Find and load your '.pem' file; it's probably in your Downloads
folder.  Note, you have to select 'All files' on the bottom.

.. image:: images/win-puttygen-2.png
   :width: 90%

Load it.

.. image:: images/win-puttygen-3.png
   :width: 90%

Now, "save private key".  Put it somewhere easy to find.

.. image:: images/win-puttygen-4.png
   :width: 90%

Logging into your EC2 instance with Putty
=========================================

Open up putty, and enter your hostname into the Host Name box.

.. image:: images/win-putty-1.png
   :width: 90%

Now, go find the 'SSH' section and enter your ppk file (generated above
by puttygen).  Then select 'Open'.

.. image:: images/win-putty-2.png
   :width: 90%

Log in as "root".

.. image:: images/win-putty-3.png
   :width: 90%

Declare victory!

.. image:: images/win-putty-4.png
   :width: 90%

