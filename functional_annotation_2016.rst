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

You will get results that look like:

.. image:: figures/lsblk.jpg
   :width: 50%

This tells you information you need to mount the volume. In this case, the mount point is /dev/xvdf. Change this if your mount point is different.
::
sudo mkdir /mnt/dammit_dbs
sudo mount /dev/xvdf /mnt/dammit_dbs/

Check out all the files you now have
::
ls /mnt/dammit_dbs/

Now start the install process.  First install linux brew, as described `here <http://angus.readthedocs.io/en/2016/linuxbrew_install.html>`__.

Next, we have a handy script stored on the volume with the rest of the install commands.
::
cp /mnt/dammit_dbs/install.sh .

Take a look at it
::
cat install.sh

And run it (this will take ~10 minutes)
::
./install.sh

