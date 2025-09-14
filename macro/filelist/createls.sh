#!/bin/bash

tag="DST_TRKR_TRACKS_run2pp_ana506_2024p023_v001"

# 逐行读取 run.list
while read run; do
    # 补齐 8 位
    run=$(printf "%08d" "$run")
    outfile="dst_trkr_tracks_run2pp-${run}.list"

    # 清空或创建文件
    > "$outfile"

    # 根据 run number 算目录范围，例如 run=00053877 对应 run_00053800_00053900
    lower_bound=$(( (10#$run / 100) * 100 ))
    upper_bound=$(( lower_bound + 100 ))
    formatted_lower=$(printf "%08d" "$lower_bound")
    formatted_upper=$(printf "%08d" "$upper_bound")

    dir="/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana506_2024p023_v001/DST_TRKR_TRACKS/run_${formatted_lower}_${formatted_upper}/dst"

    # 查找匹配的所有文件
    for f in ${dir}/${tag}-${run}-*.root; do
        if [[ -f "$f" ]]; then
            echo "$f" >> "$outfile"
        fi
    done

    echo "生成 $outfile, 包含 $(wc -l < "$outfile") 个文件"
done < run.list
