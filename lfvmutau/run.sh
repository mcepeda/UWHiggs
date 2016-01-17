#!/bin/bash

set -o nounset
set -o errexit

export jobidpu=MiniAODv2_2fb_v3
export jobid=MiniAODv2_2fb_v3
export isZTauTau=0
export isRealData=true
rake mutauData
export isRealData=false
rake mutauMC
export isZTauTau=1
rake ztautauMC

#rake mutauMC
#rake mutauMCNOPU
#export jobid=MiniAODSIMv2-Spring15-25ns_LFV_October13
#rake signal

#rake analyzeSpring2015
#rake muons
export isRealData=false
#rake mumutau
export isRealData=true
#rake mumutauData
