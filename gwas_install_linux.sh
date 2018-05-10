#!/bin/bash

mkdir gtdata
cd gtdata
wget https://www.cog-genomics.org/static/bin/plink180410/plink_linux_x86_64.zip
unzip plink_linux_x86_64.zip -d .
wget http://zzz.bwh.harvard.edu/plink/dist/hapmap_JPT_CHB_r23a.zip 
unzip hapmap_JPT_CHB_r23a.zip -d . 

