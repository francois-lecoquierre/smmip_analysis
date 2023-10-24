#!/bin/bash

# this script generates a pileup file for a single position, for all bam files in a list

# Parse the command-line arguments
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    --ref)
    ref="$2"
    shift
    shift
    ;;
    --chr)
    chr="$2"
    shift
    shift
    ;;
    --pos)
    pos="$2"
    shift
    shift
    ;;
    --bam_list)
    bam_list="$2"
    shift
    shift
    ;;
    --output_folder)
    output_folder="$2"
    shift
    shift
    ;;
    *)    # unknown option
    echo "Unknown option: $key"
    exit 1
    ;;
esac
done

# Check if all required arguments are provided
if [ -z "$ref" ] || [ -z "$chr" ] || [ -z "$pos" ] || [ -z "$bam_list" ] || [ -z "$output_folder" ]; then
    echo "Missing required arguments"
    echo "Usage: pileup_on_1_pos.sh --ref <reference> --chr <chromosome> --pos <position> --bam_list <bam_list> --output_folder <output_folder>"
    exit 1
fi

variant=$chr-$pos
echo $variant
pileup=$output_folder$variant.pileup

# Generate header
lines=$(cat $bam_list)
output="chr\tpos\tref\t"
for line in $lines; do
    IFS="/" read -ra parts <<< "$line"
    file="${parts[-1]}"
    output="${output}${file}_count\t${file}_value\t${file}_quality\t"
    output=$(echo $output | sed 's/\t$//')
    echo -e $output > $pileup
done

# Add pileup data
samtools mpileup --verbosity 2 -f "$ref" -r "$chr:$pos-$pos" -d 100000 -b $bam_list >> "$pileup"

echo "*** Saved pileup data to : $pileup"