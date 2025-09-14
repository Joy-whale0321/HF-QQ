#!/bin/bash

# 目标目录和 tag
dir="/sphenix/lustre01/sphnxpro/production/run2pp/physics/ana506_2024p023_v001/DST_TRKR_TRACKS/run_00053800_00053900/dst/"
tag="DST_TRKR_TRACKS_run2pp_ana506_2024p023_v001"

# 逐行读取 run.list
while read run; do
    # 补齐 8 位
    run=$(printf "%08d" "$run")
    outfile="dst_trkr_tracks_run2pp-${run}.list"

    # 清空或创建文件
    > "$outfile"

    # 查找匹配的所有文件
    for f in ${dir}/${tag}-${run}-*.root; do
        if [[ -f "$f" ]]; then
            echo "$f" >> "$outfile"
        fi
    done

    echo "生成 $outfile, 包含 $(wc -l < $outfile) 个文件"
done < run.list
