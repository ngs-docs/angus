================================================
Attach a Volume
================================================
Lets add a volume - as a reminder, this is sort of like plugging in an external harddrive to your laptop. Its adding a set of data, in this case, the databases that we are going to use for comparison for functional annotation.

This requires a few new steps in our tried and true Amazon EC2 instance protocol.

1. Choose AMI - Still going to choose Ubuntu server 14.04. Click Select.
2. Choose c4.2xlarge - Click Next:Configure Instance.
3. **New Step** For Subnet, select the one that ends in "us-east-1c". Click Next: Add Storage.
4. Change Size 8Gb to 100Gb.
5. **New Step** Click add Volume. Enter snapshot 'snap-6df6088a'. Make it 100Gb as well. And change to /dev/sdf
5. Click review and launch.

Now you can SSH into your instance as normal.

Make your Volume available. We made the mount point is /dev/xvdf. 
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

