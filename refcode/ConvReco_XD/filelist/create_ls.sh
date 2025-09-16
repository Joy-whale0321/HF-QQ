#!/bin/bash

localpath=$PWD
#runs=(53046 53079 53081 53195 53494 53513 53517 53530 53531 53532 53571 53578 53580 53581 53590 53630 53631 53632 53652 53687 53716 53742 53744 53756) #green runs, both good Tracking and calo runs
runs=(53018 53194 53196 53197 53586 53587 53741 53743 53783 53871 53876 53877) #yellow runs

for ((k=0; k<${#runs[@]}; k++))
do
  tagversion=ana505_2024p023_v001
  dsttype=DST_TRKR_CLUSTER
  runspecies=run2pp
  runtype=physics
  run=${runs[$k]}
  dirStart=${run:0:3}
  dirEnd=$(($dirStart + 1))
  directory=/sphenix/lustre01/sphnxpro/production/${runspecies}/${runtype}/${tagversion}/${dsttype}/run_000${dirStart}00_000${dirEnd}00/dst/
  fileHeader=${dsttype}_${runspecies}_${tagversion}-000${run}-
  outlist=${localpath}/dst_trkr_cluster_run2pp-000${run}.list
  echo "printing ${directory}/${fileHeader}*.root -> ${outlist}"
  cd ${directory}
  ls ${fileHeader}*.root > ${outlist}
done

for ((k=0; k<${#runs[@]}; k++))
do
  tagversion=ana506_2024p023_v001
  dsttype=DST_TRKR_TRACKS
  runspecies=run2pp
  runtype=physics
  run=${runs[$k]}
  dirStart=${run:0:3}
  dirEnd=$(($dirStart + 1))
  directory=/sphenix/lustre01/sphnxpro/production/${runspecies}/${runtype}/${tagversion}/${dsttype}/run_000${dirStart}00_000${dirEnd}00/dst/
  fileHeader=${dsttype}_${runspecies}_${tagversion}-000${run}-
  outlist=${localpath}/dst_trkr_tracks_run2pp-000${run}.list
  echo "printing ${directory}/${fileHeader}*.root -> ${outlist}"
  cd ${directory}
  ls ${fileHeader}*.root > ${outlist}
done
