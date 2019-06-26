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
sudo apt update
sudo apt upgrade -qyy
#cmake https://apt.kitware.com/
sudo apt-get update 
sudo apt-get install -y apt-transport-https ca-certificates gnupg software-properties-common wget
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | sudo apt-key add -
sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
sudo apt-get update
sudo apt-get install -y kitware-archive-keyring
sudo apt-key --keyring /etc/apt/trusted.gpg del C1F34CDD40CD72DA
sudo apt-get install -y cmake
#
sudo apt install -y make pkg-config emacs
sudo apt install -y samtools
sudo apt install -y tabix
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

cd /ebs1
git clone https://github.com/ComparativeGenomicsToolkit/cactus.git --recursive
cd cactus
sudo apt install -y python-pip
pip install virtualenv
virtualenv cactus_env
source cactus_env/bin/activate
pip install --upgrade toil[aws,mesos]
pip install --upgrade .
sudo apt-get install -y git gcc g++ build-essential python-dev zlib1g-dev libkyototycoon-dev libtokyocabinet-dev libkyotocabinet-dev wget valgrind libbz2-dev libhiredis-dev pkg-config
make

export PATH="/ebs1/vcflib/bin:$PATH"
export PATH="/ebs1/seqwish/bin:$PATH"
export PATH="/ebs1/minimap2:$PATH"
export PATH="/ebs1/bcftools:$PATH"



