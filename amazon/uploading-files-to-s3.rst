=====================================
Uploading files to Amazon S3 to share
=====================================

"Amazon Simple Storage Service(S3) is storage for the Internet." It is 
designed to make web-scale computing easier for developers with a
highly scalable, reliable, fast, inexpensive data storage infrastructure. 
Companies like Dropbox, Pinterest, Tumblr store their data on the S3 servers.

(Learn more from http://docs.aws.amazon.com/AmazonS3/latest/dev/Welcome.html)

For personal users like us, S3 is also a reliable, inexpensive service to 
store data or share data around the world. 

This tutorial will instruct you to upload files to Amazon S3 and share them.

Uploading files to Amazon S3
---------------------------------------------------------------

Go to Amazon S3 Console: https://console.aws.amazon.com/s3/

Click "Create Bucket".

In the dialog box, in the Bucket Name box, enter a bucket name.
Here a "bucket" is like a "folder" in concept.
The bucket name must be unique across all existing bucket names in Amazon S3.
You can not change the name after you create a bucket.
Also, the bucket name you chose here will be visible in the URL
that points to the objects(files) stored in the bucket. So make sure the 
bucket name you choose is appropriate.

Leave the "Region" as default and click "Create".

.. image:: images/create_bucket.png
   :width: 50%

Next you can click the bucket we just created, and it is an empty one.

Now we can add files into this bucket(folder) by clicking "Upload".

.. image:: images/upload_object.png
   :width: 50%
   
After you select the files and upload them successfully, you can see them
in the current bucket. Right click the file, you can manipulate it in many ways 
you like, like "Rename", "Delete", "Download", and others.

Out of them, "Make Public" is what you want to click to share the file.

.. image:: images/s3_public.png
   :width: 50%
   
When it is done, highlight the file you just shared and click "Properties" 
on the upper right. You will see the link of this file, like 
https://s3.amazonaws.com/msu_angus/MS115.bam in this example. 

.. image:: images/s3_link.png
   :width: 50%


So that's the public URL of that file. You can share this link to the person
you want to share this file with.


Downloading files from Amazon S3
---------------------------------------------------------------

You can use any internet browser to download the file from Amazon S3 directly
with the public URL.

In command-line environment, you can use curl to download the file::

  $curl -O https://s3.amazonaws.com/msu_angus/MS115.bam
  
This command wil download the file into current directory.





