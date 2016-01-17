#!/bin/bash

# Get the data
#export datasrc=/hdfs/store/user/cepeda/
#export jobid=MiniAODSIMv2-Spring15-25ns_LFV_October13
export datasrc=/hdfs/store/user/taroni
export jobid=MiniAODv2_2fb_v3
#export jobid=test
export afile=`find $datasrc/$jobid | grep root | head -n 1`

## Build the cython wrappers
#rake "make_wrapper[$afile, mt/final/Ntuple, MuTauTreeOLD]"
#rake "make_wrapper[$afile, mm/final/Ntuple, MuMuTree]"
#rake "make_wrapper[$afile, mmt/final/Ntuple, MuMuTauTree]"

#ls *pyx | sed "s|pyx|so|" | xargs rake 
#echo "finishing compilation" 
#bash compileTree.txt

#rake "meta:getinputs[$jobid, $datasrc,mt/metaInfo, mt/summedWeights]"
rake "meta:getmeta[inputs/$jobid, mt/metaInfo, 13, mt/summedWeights]"


