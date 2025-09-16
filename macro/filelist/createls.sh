#!/bin/bash

localpath=$PWD

# 逐行读取 run.list
while read run || [[ -n "$run" ]]; do
    # 如果是空行就跳过
    [[ -z "$run" ]] && continue

    # 补齐 8 位
    run=$(printf "%08d" "$run")

    # 根据 run number 算目录范围，例如 run=00053877 对应 run_00053800_00053900
    lower_bound=$(( (10#$run / 100) * 100 ))
    upper_bound=$(( lower_bound + 100 ))
    formatted_lower=$(printf "%08d" "$lower_bound")
    formatted_upper=$(printf "%08d" "$upper_bound")

    # CLUSTER dst files output list creation ------------------
    outfile=${localpath}/dst_trkr_cluster_run2pp-${run}.list
    > "$outfile"

    tagversion=ana505_2024p023_v001
    dsttype=DST_TRKR_CLUSTER
    runspecies=run2pp
    runtype=physics
    directory="/sphenix/lustre01/sphnxpro/production/${runspecies}/${runtype}/${tagversion}/${dsttype}/run_${formatted_lower}_${formatted_upper}/dst/"

    fileHeader=${dsttype}_${runspecies}_${tagversion}-${run}-
    echo "printing ${directory}/${fileHeader}*.root -> ${outfile}"
    cd ${directory}
    ls ${fileHeader}*.root > ${outfile}
done < run.list

# 逐行读取 run.list
while read run || [[ -n "$run" ]]; do
    # 如果是空行就跳过
    [[ -z "$run" ]] && continue
    # 补齐 8 位
    run=$(printf "%08d" "$run")

    # 根据 run number 算目录范围，例如 run=00053877 对应 run_00053800_00053900
    lower_bound=$(( (10#$run / 100) * 100 ))
    upper_bound=$(( lower_bound + 100 ))
    formatted_lower=$(printf "%08d" "$lower_bound")
    formatted_upper=$(printf "%08d" "$upper_bound")

    # TRACKS dst files output list creation ------------------
    outfile=${localpath}/dst_trkr_tracks_run2pp-${run}.list
    > "$outfile"

    tagversion=ana506_2024p023_v001
    dsttype=DST_TRKR_TRACKS
    runspecies=run2pp
    runtype=physics
    directory="/sphenix/lustre01/sphnxpro/production/${runspecies}/${runtype}/${tagversion}/${dsttype}/run_${formatted_lower}_${formatted_upper}/dst/"

    fileHeader=${dsttype}_${runspecies}_${tagversion}-${run}-
    echo "printing ${directory}/${fileHeader}*.root -> ${outfile}" 
    cd ${directory}
    ls ${fileHeader}*.root > ${outfile}

    echo "生成 $outfile, 包含 $(wc -l < "$outfile") 个文件"
done < run.list
