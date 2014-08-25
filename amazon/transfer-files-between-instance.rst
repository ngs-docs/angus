=====================================
Transferring Files between your laptop and Amazon instance 
=====================================

For linux/Unix/Mac system, we can use a command-line tool "scp" to transfer 
files between your laptop and Amazon instance. Also we can use a GUI tool 
"FileZilla" to do the transfer, which is more user-friendly.

Using scp to transfer data
-------------

"scp" means "secure copy", which can copy files between computers on a network.
You can use this tool in a Terminal on a Unix/Linux/Mac system.

To upload a file from your laptop to Amazon instance::

  $scp -i ~/Desktop/amazon.pem ~/Desktop/MS115.fa  ubuntu@ec2-54-166-128-20.compute-1.amazonaws.com:~/data/
  
This command will upload a file - MS115.fa in your ~/Desktop/ folder of 
your laptop to folder ~/data/ on an Amazon instance. Note you still need to 
use the private key you used to connect to the Amazon instance with ssh. (In
this example, it is the amazon.pem file in ~/Desktop/.

Note: You need to make sure that the user "ubuntu" has the permission to 
write in the target directory. In this example, if ~/data/ was created by user 
"ubuntu", it should be fine. 

Similarly, to download a file from Amazon instance to your laptop::

  $scp -i ~/Desktop/amazon.pem ubuntu@ec2-54-166-128-20.compute-1.amazonaws.com:/data/ecoli_ref-5m-trim.fastq.gz ~/Download/
  
This command will download a file /data/ecoli_ref-5m-trim.fastq.gz from 
Amazon instance to your ~/Download folder in your laptop.

Note: You can use asterisk(*) to download multiple files, like *.fasta.gz.


Using FileZilla to transfer data
---------------

If you want a more user-friendly tool to transfer data, FileZilla is a 
good choice. It is free, it supports Windows/Linux/Mac systems, and it has a
good user interface. It supports FTP, SFTP and other file transfer
protocols. 

Firstly, go to 'https://filezilla-project.org/' and click "Download FileZilla 
Client" button to download it.

The interface of FileZilla is like this:

.. image:: images/filezilla_interface.png
   :width: 50%

If you want to use FileZila to upload to or download data from a normal 
FTP server if you have the user and password, just put the information in the 
"Host", "Username", "Password" box and connect. 
However for Amazon instance, we use key-pair to log in 
instead of password for better safety. So it is a little bit more complicated
to configure.

Open "Settings" and click "SFTP":

.. image:: images/filezilla_sftp.png
   :width: 50%

Click "Add keyfile...":

.. image:: images/add_keyfile.png
   :width: 50%

Then select the ".pem" file you used to connect to Amazon instance with ssh.

.. image:: images/select_pem.png
   :width: 50%

There is a dialog box to ask you if you want to convert the ".pem" file 
into a supported format. Click "Yes".

.. image:: images/convert_pem.png
   :width: 50%


Name it with extension as ".ppk" and save it.

.. image:: images/save_ppk.png
   :width: 50%

You will see the a private key has been added.

.. image:: images/show_key.png
   :width: 50%

Close "Settings" and go back to the main interface.

Click button to open the site manager.

.. image:: images/site_manager.png
   :width: 50%

Click "New Site".

.. image:: images/new_site.png
   :width: 50%

Put the Amazon instance URL like ec2-54-166-128-20.compute-1.amazonaws.com 
in the "Host" box. Set "Protocol" as "SFTP", "Logon Type" as "Normal", 
"User" as "ubuntu" and leave "Password" as blank. Then click "Connect".

.. image:: images/site_config.png
   :width: 50%

There will be a dialogue box to ask you about "Unknown host key", just click 
"Ok". 

.. image:: images/host_key_diag.png
   :width: 50%

All right. Now you have logged in the Amazon instance. You can drag and drop
to transfer the files between the remote machine and your local laptop.

.. image:: images/final_filezilla.png
   :width: 50%

