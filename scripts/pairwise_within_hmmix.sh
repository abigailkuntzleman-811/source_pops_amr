#!/bin/bash

ml bedops 

# Get the populations from the command line arguments
pop="$1"
input_dir="$2"
rootoutdir="$3"

# Paths to population directories
dir="$input_dir/$pop"

# Create a result directory for the pair
output_dir="$rootoutdir/${pop}"
mkdir -p "$output_dir"

# Initialize total intersection and total combinations variables
total_intersection=0
total_combinations=0

# Get all bed files in the population folder
bed_files=$(find "$dir" -name "*.bed")

individuals_file="$rootoutdir/all_$pop.txt"
> "$individuals_file"

# Compute intersection for each combination of BED files
for bed1 in $bed_files; do
    for bed2 in $bed_files; do
        if [ "$bed1" != "$bed2" ] && [[ "$bed1" > "$bed2" ]]; then
            # Output file for the intersection
            output_file="$output_dir/intersection_$(basename "$bed1")_$(basename "$bed2").bed.gz"

            # Run bedops intersection command
            bedops --intersect "$bed1" "$bed2" | gzip > "$output_file"

            # Get the number of base pairs in the intersection
            intersect_base_pairs=$(zcat "$output_file" | awk '{sum += $3 - $2} END {print sum}')

            # rm "output_file"
            total_intersection=$((total_intersection + intersect_base_pairs))
            total_combinations=$((total_combinations + 1))

            id1=$(basename "$bed1" .bed)
            id2=$(basename "$bed2" .bed)

            echo -e "${id1}\t${id2}\t${intersect_base_pairs}" >> "$individuals_file"
            rm "$output_file"
        fi
    done
done

# # Calculate the average for this pair of populations
if ((total_combinations > 0)); then
    avg_pairwise_sharing=$(echo "scale=6; $total_intersection / $total_combinations" | bc)

#     # Save the result to a file with the name pop1_pop2_pairwise_average_sharing.txt
    output_file="$rootoutdir/${pop}_pairwise_average_sharing.txt"
    echo -e "$pop\t$avg_pairwise_sharing" > "$output_file"

    echo "Pairwise average sharing result saved to $output_file"
else
    echo "No intersections found between $pop1 and $pop2"
fi