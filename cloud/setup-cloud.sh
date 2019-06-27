#!/bin/bash
#
# This will install a bunch of graph-genome tools and dependencies
# via github and apt.  By default it will install everything to /ebs1
# (see setup-volumes-i3.4xlarge.sh for how /ebs1 is mounted)
# but that can be changed via search-replace below
#
# cactus and vg will need activation of their respective virtualenvs to run
# the other tools need to be in the path (see exports at bottom)
#
# This process takes about 10 minutes or so.  You need to click OK at the
# beginning for the apt upgrade, and re-log back in after if you want to use Docker
#
# Tested on the "Ubuntu Server 18.04 LTS (HVM), SSD Volume Type - ami-005bdb005fb00e791 (64-bit x86)" image
#
sudo apt update
sudo apt upgrade -qyy
sudo apt install -y make pkg-config emacs virtualenv
sudo apt install -y samtools
sudo apt install -y tabix
sudo apt install -y libssl-dev
sudo apt install -y docker.io
sudo usermod -aG docker $USER
printf "\n\nNeed to log out then in for Docker to work!!\n\n"

mkdir /ebs1

cd /ebs1
git clone https://github.com/vgteam/vg.git --recursive
cd vg
make get-deps
make -j $(getconf _NPROCESSORS_ONLN)

cd /ebs1
git clone https://github.com/vcflib/vcflib.git --recursive
cd vcflib
make -j $(getconf _NPROCESSORS_ONLN)

cd /ebs1
git clone https://github.com/vgteam/toil-vg.git
cd toil-vg
virtualenv toilvenv
source toilvenv/bin/activate
make prepare
make develop

cd /ebs1
git clone https://github.com/ekg/seqwish.git --recursive
cd seqwish
cmake -H. -Bbuild && cmake --build build -- -j $(getconf _NPROCESSORS_ONLN)

cd /ebs1
git clone git://github.com/samtools/htslib.git
git clone git://github.com/samtools/bcftools.git
cd bcftools
make

cd /ebs1
git clone https://github.com/lh3/minimap2
cd minimap2
make

# these packages aren't in apt in ubuntu 18.04.  but if we move to ubuntu 16.04, we can't
# build vg due to cmake version issues.  so we install by hand
cd /ebs1
mkdir local
git clone https://github.com/cloudflarearchive/kyotocabinet.git --recursive
cd kyotocabinet
./configure --prefix=/ebs1/local
make -j $(getconf _NPROCESSORS_ONLN) && make install

cd /ebs1
git clone https://github.com/cloudflarearchive/kyototycoon.git --recursive
cd kyototycoon
./configure --prefix=/ebs1/local --with-kc=/ebs1/local/
make -j $(getconf _NPROCESSORS_ONLN) && make install

cd /ebs1
git clone https://github.com/ComparativeGenomicsToolkit/cactus.git --recursive
cd cactus
sudo apt install -y python-pip
pip install virtualenv
virtualenv cactus_env
source cactus_env/bin/activate
pip install --upgrade toil[aws,mesos]
pip install --upgrade .
ttPrefix=/ebs1/local
export kyotoTycoonIncl="-I${ttPrefix}/include -DHAVE_KYOTO_TYCOON=1"
export kyotoTycoonLib="-L${ttPrefix}/lib -Wl,-rpath,${ttPrefix}/lib -lkyototycoon -lkyotocabinet -lz -lbz2 -lpthread -lm -lstdc++"
sudo apt-get install -y git gcc g++ build-essential python-dev zlib1g-dev libtokyocabinet-dev libkyotocabinet-dev wget valgrind libbz2-dev libhiredis-dev pkg-config
make

cd /ebs1
git clone https://github.com/joelarmstrong/repeatMaskerPipeline.git

cd /ebs1
git clone https://github.com/ComparativeGenomicsToolkit/hal2vg.git --recursive
cd hal2vg
export PATH="/ebs1/cactus/submodules/hdf5/bin:$PATH"
export h5prefix="-prefix=/ebs1/cactus/submodules/hdf5"
make

export PATH="/ebs1/vcflib/bin:$PATH"
export PATH="/ebs1/seqwish/bin:$PATH"
export PATH="/ebs1/minimap2:$PATH"
export PATH="/ebs1/bcftools:$PATH"
export PATH="/ebs1/hal2vg:$PATH"



