=====================================
Using Amazon EBS Snapshots for sharing and backing up data
=====================================

Now you have the Amazon EBS Volume to store the data for your instance to use. 
But it can only be attached to the EC2 instance you created. 
If you want to share the EBS Volume with the data with other colleagues so 
they can use the data on the EC2 instance they created, you need to create
an Amazon EBS Snapshot and share it with other Amazon EC2 user. You can also
create Amazon EBS Snapshots periodically to backup your Amazon EBS Volume.

In this tutorial, you will be instructed to create an Amazon EBS Snapshot
from an Volume, to share the EBS Snapshot with others and to restore an 
Volume from an Snapshot.


Prerequisites
-------------

This tutorial assumes you've completed the EC2 tutorial to set up an
Amazon EBS Volume.

Creating an Amazon EBS Snapshot from a Volume
---------------------------------------------------------------


Firstly open the Amazon EC2 console at 'https://console.aws.amazon.com/ec2'
 and make sure it says North Virginia in the upper right.

At the AWS Management Console, on the left menu bar, click "Volumes".



Here you can see the 20GiB volume you created in the tutorial 
"Storing data persistently with Amazon EBS Volumes". 

.. image:: images/volume_view.png
   :width: 50%
   
This tutorial will guide you to create a Amazon EBS Snapshot for this volume.
In this example,
the "Volume ID" of that volume is "vol-180eef57". Record this ID, we will 
use it later.


Next, on the left menu bar, click "Snapshots" and click "Create Snapshot".

.. image:: images/create_snapshot.png
   :width: 50%

Choose the Volume we want to create Snapshot for. (The ID is vol-180eef57 for
this example, as we recored in last step.) You can enter the information 
for "Name" and "Description" as you like or just leave them blank.

Then click "Create".

.. image:: images/snapshot_option.png
   :width: 50%

Ok. Now you have created a Snapshot from that 20G Volume. 

Sharing an Amazon EBS Snapshot
---------------------------------------------------------------

For now, the snapshot we just created is private, which can only be viewed
and used by yourself. 

To make it public or share it with other Amazon EC2 user:

Right click the snapshot and click "Modify Snapshot Permissions".

.. image:: images/snapshot_share.png
   :width: 50%
   

If you want to make this snapshot public, which means any Amazon EC2 user can 
have access to this snapshot and all the data in it, click "Public".

If you just want to share this snapshot with specific person, like your 
colleague, keep the option as "Private" but put the AWS Account Number (which
can be acquired by checking "Account Settings") of the
person in the box and click "Add Permission".

.. image:: images/modify_snapshot_permission.png
   :width: 50%

Now you can share the "Snapshot ID" (snap-a8a00b07 in this example)
to your colleague so they can restore
an EBS Volume from this snapshot and use the data in it.


Restoring an Amazon EBS Volume from a Snapshot
---------------------------------------------------------------

Your colleague shares an Amazon EBS Snapshot with you and you want to use
the data on it. You did something terribly wrong to the data on your Volume
and you want to restore the data from the backup Snapshot of that Volume. 
Under these circumstances, you want to restore an EBS Volume from a Snapshot.

It is similar to how you create a new EBS Volume. The only difference is that 
in the dialog box after you click "Create Volume" button, input the 
"Snapshot ID" of the snapshot you want to restore. Similarly, also select the 
zone in which your instance is running. 

.. image:: images/restore_snapshot.png
   :width: 50%

Ok, now you have the volume available to attach to your running instance.
