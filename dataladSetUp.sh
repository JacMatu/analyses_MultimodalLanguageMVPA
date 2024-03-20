#!/usr/bin/env bash

root_dir=${PWD}
raw_dir=${root_dir}/inputs/raw
derivatives_dir=${root_dir}/outputs/derivatives
#preproc_dir=${derivatives_dir}/bidspm-preproc
#stats_dir=${derivatives_dir}/bidspm-stats
#roi_dir=${derivatives_dir}/bidspm-roi
#cosmo_dir=${derivatives_dir}/cosmo-mvpa

# get url of the gin repos from config
source dataladConfig.sh

# install raw dataset
datalad install -d . -s "${URL_RAW}" "${raw_dir}"

# install derivatives dataset with submodules (hopefully)

datalad install -r -d . -s "${URL_DER}" "${derivatives_dir}"

cd "${derivatives_dir}"

datalad push --to origin -r

cd "${root_dir}"

datalad save . -m "add code and folders to set subdatasets"

datalad push --to origin

echo "############################"
echo "# DATALAD IS READY TO WORK #"
echo "############################"
