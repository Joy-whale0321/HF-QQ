#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n new  # setup sPHENIX environment in the singularity container shell. Note the shell is bash by default

# Additional commands for my local environment
export SPHENIX=/sphenix/u/jzhang1
export MYINSTALL=$SPHENIX/install

# Setup MYINSTALL to local directory and run sPHENIX setup local script
# to adjust PATH, LD LIBRARY PATH, ROOT INCLUDE PATH, etc
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

echo "sPHENIX environment setup finished"

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`
echo running: $this_script $*

nEvents=$1
InTrackDst=$2
InTrackPath=$3
InCaloDst=$4
InCaloPath=$5
OutPrefix=$6
OutPath=$7
Index=${8}
StepSize=${9}

# echo "OFFLINE_MAIN=$OFFLINE_MAIN"
# echo "$LD_LIBRARY_PATH" | grep -o 'release_[^/]*/' | head -1

root.exe -q -l Fun4All_Quarkonium_pp.C\($nEvents,\"${InTrackDst}\",\"${InTrackPath}\",\"${InCaloDst}\",\"${InCaloPath}\",\"${OutPrefix}\",\"${OutPath}\",$Index,$StepSize\)

echo Script done
