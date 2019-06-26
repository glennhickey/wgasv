#!/bin/bash

# This assumes an i3.4xlarge instance with one EBS
# The volume names are hardcoded to the defaults (now).
# they can be checked with: sudo fdisk -l
# The EBS volume is only formatted if unititialized
#
# Todo: would like to figure out how to mount all local disks
# to one volume. 
#
# The idea is to have all the installed tools and stuff on the ebs volume
# and big working data on the ephemeral volumes
#
sudo mkdir /data1
sudo mkfs.ext4 /dev/nvme0n1
sudo mount -t ext4 /dev/nvme0n1 /data1
echo "/dev/nvme0n1 /data1 auto noatime 0 0" | sudo tee -a /etc/fstab

sudo mkdir /data2
sudo mkfs.ext4 /dev/nvme1n1
sudo mount -t ext4 /dev/nvme1n1 /data2
echo "/dev/nvme1n1 /data1 auto noatime 0 0" | sudo tee -a /etc/fstab

sudo mkdir /ebs1
if [ "$(sudo file -s /dev/xvdb)" == "/dev/xvdb: data" ]
then
	 sudo mkfs.ext4 /dev/xvdb
fi
sudo mount -t ext4 /dev/xvdb /ebs1
echo "/dev/xvdb /ebs1 auto noatime 0 0" | sudo tee -a /etc/fstab

df -h
