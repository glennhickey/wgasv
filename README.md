### SV graphs from genome alignments

Worflow is as follows:
* Create an i3.4xlarge instance using the blue Launch Instance button from the EC2 console
* In Step 1, select this AMI **Ubuntu Server 18.04 LTS (HVM), SSD Volume Type - ami-005bdb005fb00e791 (64-bit x86)**
* Select an **i3.4xlarge** insteance type in Step 2.  This isn't so important, but other instance types may require hacking the device names in the mounting script below a little bit.
* Click "Configure Instance Details" and tick off the "Request Spot Instances" box in the Purchasing Options section near the top.
* Click "Next: Add Storage".  And add an EBS volume if desired (note there's lots of ephemeral space on this node type).  If you already have a volume, don't do anything here and attach it to your instance later via the EBS->Volumes menu. 
* Launch the instance. Once it's running, add the Owner and Name tags from the instances menu. 
* Copy the instance's public IP from the Description part of the instance information (when clicking on it).

Once the instance is successfully create:
```
ssh ubuntu@<public ip>
git clone https://github.com/glennhickey/wgasv.git
cd wgasv

# mount the disks
cd cloud
./setup-volumes-i3.4xlarge.sh
# make sure the df output looks sensible

# install the software
./setup-cloud.sh
```
From here on, you need to make sure you're using somewhere on `/data1/` or `/data2/` for scratch space and general storage.  You may need to run something like `sudo chmod -R 777 /data/1` to get permissions.  

#### Notes
* You need to log out then back in again after having run `setup-cloud.sh` before using Docker.
* Getting Kyotocabinet to compile required the following manual patch for me

```
diff --git a/kcdbext.h b/kcdbext.h
index 001c09a..93c612d 100644
--- a/kcdbext.h
+++ b/kcdbext.h
@@ -1278,7 +1278,7 @@ class IndexDB {
     if (omode_ == 0) {
       set_error(_KCCODELINE_, BasicDB::Error::INVALID, "not opened");
       *sp = 0;
-      return false;
+      return NULL;
     }
     if (!cache_) return db_.get(kbuf, ksiz, sp);
     size_t dvsiz = 0;
```

* I have not gotten toil-vg to run docker images due to some weird permission thing (they run okay on command line).  But the scripts should build everything needed to run most commands without Docker
