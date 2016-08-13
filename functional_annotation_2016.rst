================================================
Attach a Volume
================================================
Lets add a volume - as a reminder, this is sort of like plugging in an external harddrive to your laptop. Its adding a set of data, in this case, the databases that we are going to use for comparison for functional annotation.

This requires a few new steps in our tried and true Amazon EC2 instance protocol.

1. Choose AMI - Still going to choose Ubuntu server 14.04. Click Select.
2. Choose c4.2xlarge - Click Next:Configure Instance.
3. **New Step** For Subnet, select the one that ends in "us-east-1c". Click Next: Add Storage.
4. Change Size 8Gb to 100Gb.
5. Click review and launch.

Now go to the EC2 dashboard and click on Volumes.  
1. Create a new volume
2. Make the size 20gb 
3. For the availability zone, enter us-east-1c
4. For snapshot id, enter ``snap-6df6088a``. A little drop down should say "Angus2016DammitDBs", select the dropdown option.
4. Click create

Now we need to attach the volume to the instance. While still in the Volume part of the dashboard:
1. Actions-> Attach Volume
2. Click on instance, it should give you a single option, the instance you are already running
3. Click attach

Now you can SSH into your instance as normal.


================================================
Make your Volume available
================================================

Start by seeing where you Volume is hanging out
::
lsblk

