# Motif Grange object make
# The positions of the motifs were taken out and made into Grange objects
#!/bin/bash

# The value of specified to match the fourth column
target_motifs=("GB48655_M06252_2.00")  # Here to add your value to match

# Specifies the path to two BED files
bed_file1="Or_LOC102656221_promoter_motif.bed"
bed_file2="Or_LOC102656904_promoter_motif.bed"

# Specifies the path to the output file
output_file="LOC552079_motifs.bed"

# Clear the output file (if present)
> "$output_file"

# iterate over each motif we want to match
for motif in "${target_motifs[@]}"; do
    # We look for matching lines in the first BED file and append them to the output file
    awk -v motif="$motif" -F '\t' -v OFS='\t' '{ if ($4 == motif) { split($1, region, ":|-"); $1 = region[1]; $2 = $2 + region[2]; $3 = $3 + region[2]; print } }' "$bed_file1" >> "$output_file"

    # We look for matching lines in the second BED file and append them to the output file
    awk -v motif="$motif" -F '\t' -v OFS='\t' '{ if ($4 == motif) { split($1, region, ":|-"); $1 = region[1]; $2 = $2 + region[2]; $3 = $3 + region[2]; print } }' "$bed_file2" >> "$output_file"
done


# The positions of the motifs were taken out and made into Grange objects
#!/bin/bash

# The value of specified to match the fourth column
target_motifs=("GB45715_M03994_2.00")  # Here to add your value to match

# Specifies the path to two BED files
bed_file1="Or_LOC102656221_promoter_motif.bed"
bed_file2="Or_LOC102656904_promoter_motif.bed"

# Specifies the path to the output file
output_file="LOC411133_motifs.bed"

# Clear the output file (if present)
> "$output_file"

# iterate over each motif we want to match
for motif in "${target_motifs[@]}"; do
    # We look for matching lines in the first BED file and append them to the output file
    awk -v motif="$motif" -F '\t' -v OFS='\t' '{ if ($4 == motif) { split($1, region, ":|-"); $1 = region[1]; $2 = $2 + region[2]; $3 = $3 + region[2]; print } }' "$bed_file1" >> "$output_file"

    # We look for matching lines in the second BED file and append them to the output file
    awk -v motif="$motif" -F '\t' -v OFS='\t' '{ if ($4 == motif) { split($1, region, ":|-"); $1 = region[1]; $2 = $2 + region[2]; $3 = $3 + region[2]; print } }' "$bed_file2" >> "$output_file"
done

# The positions of the motifs were taken out and made into Grange objects
#!/bin/bash

# The value of specified to match the fourth column
target_motifs=("GB50892_M03837_2.00")  # Here to add your value to match

# Specifies the path to two BED files
bed_file1="Or_LOC102656221_promoter_motif.bed"
bed_file2="Or_LOC102656904_promoter_motif.bed"

# Specifies the path to the output file
output_file="LOC100577980_motifs.bed"

# Clear the output file (if present)
> "$output_file"

# iterate over each motif we want to match
for motif in "${target_motifs[@]}"; do
    # We look for matching lines in the first BED file and append them to the output file
    awk -v motif="$motif" -F '\t' -v OFS='\t' '{ if ($4 == motif) { split($1, region, ":|-"); $1 = region[1]; $2 = $2 + region[2]; $3 = $3 + region[2]; print } }' "$bed_file1" >> "$output_file"

    # We look for matching lines in the second BED file and append them to the output file
    awk -v motif="$motif" -F '\t' -v OFS='\t' '{ if ($4 == motif) { split($1, region, ":|-"); $1 = region[1]; $2 = $2 + region[2]; $3 = $3 + region[2]; print } }' "$bed_file2" >> "$output_file"
done


