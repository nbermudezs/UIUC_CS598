#!/usr/bin/env bash

mkdir data/raw
mkdir data/raw/iRefIndex
mkdir data/raw/MultiNet
mkdir data/raw/HINT
mkdir data/raw/TCGA

echo "Downloading iRefIndex data for Homo sapiens"
curl -L -k "https://uofi.box.com/shared/static/i6c1xwe03nd6brkqg0lwdz4glho9mfr6.zip" --output data/raw/iRefIndex/9606.mitab.zip
download_md5=$(md5 data/raw/iRefIndex/9606.mitab.zip | cut -d\  -f4)
ref_md5=$(cat data/iRefIndex.md5)

if download_md5!=ref_md5; then
    echo "Downloaded iRefIndex data has wrong MD5 checksum"
    exit 1
fi


echo "Downloading Multinet data"
curl -L -k "https://uofi.box.com/shared/static/bseep4fprvxrgw8wgmzykfs9oumlch1n.txt" --output data/raw/MultiNet/Multinet.interactions.txt
curl -L -k "https://uofi.box.com/shared/static/9wg67jp2ikta4r8l9cs7ip56bpr4rjai.txt" --output data/raw/MultiNet/Multinet.interactions.network_presence.txt


echo "Downloading HINT data"
curl -L -k "https://uofi.box.com/shared/static/mb26vwyt45a4o09lsonsicvwzg8anypb.txt" --output data/raw/HINT/HomoSapiens_binary_hq.txt
curl -L -k "https://uofi.box.com/shared/static/vz198ectq2htmgslbnkct4zz4xwapnhz.txt" --output data/raw/HINT/HomoSapiens_cocomp_hq.txt
curl -L -k "https://uofi.box.com/shared/static/318ywc0fc9dc6ip35y3udyavsk0345ro.txt" --output data/raw/HINT/HomoSapiens_lcc_hq.txt
curl -L -k "https://uofi.box.com/shared/static/emtf55y3hyar5e6s3b4uckn38zdr91l9.txt" --output data/raw/HINT/HomoSapiens_lcb_hq.txt


echo "Downloading TCGA Pancancer data"
curl -L -k "https://uofi.box.com/shared/static/meb2mims3sag5ebd0iui4j1efr5mdg27" --output data/raw/TCGA/pancancer12.txt


echo "Downloading TCGA MAFs"
curl -L -k "https://uofi.box.com/shared/static/2hzczmzv64udzoabvwge31rla944znho.zip" --output data/raw/TCGA/mafs.zip
unzip data/raw/TCGA/mafs.zip -d data/raw/TCGA/
rm data/raw/TCGA/mafs.zip
