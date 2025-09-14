#!/bin/bash

runs=(53046) #test run
# runs=(53046 53079 53081 53195 53494 53513 53517 53530 53531 53532 53571 53578 53580 53581 53590 53630 53631 53632 53652 53687 53716 53742 53744 53756) #green runs, both good Tracking and calo runs
#runs=(53018 53194 53196 53197 53586 53587 53741 53743 53783 53871 53876 53877) #yellow runs

# set the ratio of number of calo files over number of trkr files
ratio_nCalo_over_nTrkr=10

# loop over the selected runs
for ((k=0; k<${#runs[@]}; k++))
do
    # input list file settings
    input_trkr_tracks_file="dst_trkr_tracks_run2pp-000${runs[$k]}.list"
    # input_trkr_cluster_file="dst_trkr_cluster_run2pp-000${runs[$k]}.list"
    # input_calo_file="dst_calo_run2pp-000${runs[$k]}.list"
    input_calo_file="dst_calo-000${runs[$k]}.list"
    # input_calo_file="dst_calofitting-000${runs[$k]}.list"

    # check if input files exist, skip missing runs and segments files
    if [ ! -f "$input_trkr_tracks_file" ] || [ ! -f "$input_calo_file" ];  then
        continue
    fi

    # synchronize(sync) list file settings
    output_trkr_tracks_cluster_calo_file="dst_sync_trkr_calo_run2pp-000${runs[$k]}.list"

    # clear output file, if not exist, create a new one
    > "$output_trkr_tracks_cluster_calo_file"

    # read trkr tracks file and build synchronization map
    # remove old trkr_tracks_map and create a new one
    unset trkr_tracks_map 
    declare -A trkr_tracks_map
    # initialize variables
    trkr_tracks_prefix=0
    trkr_tracks_run=0
    # while ... done < file, read file line by line from input_trkr_tracks_file
    while IFS= read -r trkr_tracks_line; do        
        # cut -d'-' -f1、-f2、-f3 abstract the prefix、run number、segment number from the line(filename)
        trkr_tracks_prefix=$(echo "$trkr_tracks_line" | cut -d'-' -f1)
        trkr_tracks_run=$(echo "$trkr_tracks_line" | cut -d'-' -f2)
        trkr_tracks_segment=$(echo "$trkr_tracks_line" | cut -d'-' -f3 | cut -d'.' -f1)
        # eg.
        # read DST_TRKR_SEED_run2pp_ana468_2024p012_v002-00053741-00001.root
        # prefix: DST_TRKR_SEED_run2pp_ana468_2024p012_v002
        # run number: 00053741
        # segment number: 00001

        echo "[DEBUG] trkr_tracks_line=$trkr_tracks_line -> segment=$trkr_tracks_segment"

        # store in synchronization map, key is segment
        trkr_tracks_map["$trkr_tracks_segment"]="$trkr_tracks_line"
    done < "$input_trkr_tracks_file"

    # read calo file and build synchronization map (workflow similar as above)
    unset calo_map
    declare -A calo_map
    calo_prefix=0
    calo_run=0
    while IFS= read -r calo_line; do
        calo_prefix=$(echo "$calo_line" | cut -d'-' -f1)
        calo_run=$(echo "$calo_line" | cut -d'-' -f2)
        calo_segment=$(echo "$calo_line" | cut -d'-' -f3 | cut -d'.' -f1)

        echo "[DEBUG] calo_line=$calo_line -> segment=$calo_segment"

        # store in synchronization map, key is segment
        calo_map["$calo_segment"]="$calo_line"
    done < "$input_calo_file"

    # loop over all tracker segments, write out only when calo can also match with increased order
    for trkr_segment in $(printf "%s\n" "${!trkr_tracks_map[@]}" | sort -n); do
        matched_calo_segment=$(printf "%05d" $((10#$trkr_segment / ratio_nCalo_over_nTrkr)))
        # 调试输出可保留观察次序
        echo "[DEBUG-ORDERED] trkr_segment=${trkr_segment} -> calo_segment=${matched_calo_segment}"
        if [[ -n "${calo_map[$matched_calo_segment]}" ]]; then
            echo "$trkr_segment" "$matched_calo_segment" >> "$output_trkr_tracks_cluster_calo_file"
        fi
    done

    njobs=`cat $output_trkr_tracks_cluster_calo_file | wc -l`
    echo "Run ${runs[$k]} match done! Writing to $output_trkr_tracks_cluster_calo_file. ${njobs} Jobs in total"

done
