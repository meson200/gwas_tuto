#!/bin/bash

mkdir gtdata
cd gtdata
curl https://www.cog-genomics.org/static/bin/plink180410/plink_mac.zip -o plink_mac.zip
tar --strip-components=1 -xvzf plink_mac.zip 
curl http://zzz.bwh.harvard.edu/plink/dist/hapmap_JPT_CHB_r23a.zip -o hapmap_asian.zip
tar -xvzf hapmap_asian.zip 


