#!/bin/bash

#This script will generate all pairwise combinations between populations - the result is 
#a bunch of bed files under /users/akuntzle/scratch/IBD_archaic_shared/pairwise_comparisons_hmmix/{pop1_pop2}

#There's also a file at that directory called pop1_pop2.txt which has all the individuals and how much they share
ml bedops 
# Get the populations from the command line arguments
pop1="$1"
pop2="$2"
input_dir="$3"
output_dir="$4"

echo "$pop1 $pop2 $input_dir $output_dir"

# Paths to population directories
dir1="$input_dir/$pop1"
dir2="$input_dir/$pop2"

# Create a result directory for the pair
pair_output_dir="$output_dir/${pop1}_${pop2}"
mkdir -p "$pair_output_dir"

# Initialize total intersection and total combinations variables
total_intersection=0
total_combinations=0

# Get all bed files in each population folder
bed_files1=$(find "$dir1" -name "*.bed")
bed_files2=$(find "$dir2" -name "*.bed")

individuals_file="$output_dir/${pop1}_{$pop2}.txt"
> "$individuals_file"

# Compute intersection for each combination of BED files
for bed1 in $bed_files1; do
    for bed2 in $bed_files2; do
        # Output file for the intersection
        id1=$(basename "$bed1" .bed)
        id2=$(basename "$bed2" .bed)

        output_file="$pair_output_dir/intersection_${id1}_${id2}.bed.gz"

        # Run bedops intersection command
        bedops --intersect "$bed1" "$bed2" | gzip > "$output_file"

        # Get the number of base pairs in the intersection
        intersect_base_pairs=0
        intersect_base_pairs=$(zcat "$output_file" | awk '{sum += ($3 - $2)} END {print sum}')

        echo -e "${id1}\t${id2}\t${intersect_base_pairs}" >> "$individuals_file"

        total_intersection=$((total_intersection + intersect_base_pairs))
        total_combinations=$((total_combinations + 1))

        rm "$output_file"
    done
done

# Calculate the average for this pair of populations
if ((total_combinations > 0)); then
    avg_pairwise_sharing=$(echo "scale=6; $total_intersection / $total_combinations" | bc)

    # Save the result to a file with the name pop1_pop2_pairwise_average_sharing.txt
    output_file="$output_dir/${pop1}_${pop2}_pairwise_average_sharing.txt"
    echo -e "$pop1\t$pop2\t$avg_pairwise_sharing" > "$output_file"

    echo "Pairwise average sharing result saved to $output_file"
else
    echo "No intersections found between $pop1 and $pop2"
fi
