#!/bin/bash

# Output file header
echo -e "Sample\tHDM_Read_Count" > hdm_read_counts.tsv

for file in *.blastn; do
    sample=$(basename "$file" .blastn)

    # Sort by qseqid and descending bitscore (column 12 is bitscore)
    sort -k1,1 -k12,12nr "$file" > tmp_sorted.txt

    # Parse and count only reads whose top bitscore hit(s) are all Dermatophagoides
    awk -F'\t' '
    {
        qseqid = $1
        bitscore = $12
        ssciname = $15

        # If this is a new qseqid, evaluate the previous one
        if (qseqid != prev_qseqid) {
            # Count previous qseqid if all top hits were Dermatophagoides
            if (NR > 1 && top_is_derma == 1) {
                count++
            }

            # Start tracking new qseqid
            prev_qseqid = qseqid
            max_bitscore = bitscore
            top_is_derma = (ssciname ~ /Dermatophagoides farinae/) ? 1 : 0
        }
        else {
            # If bitscore still equals top bitscore, keep checking ssciname
            if (bitscore == max_bitscore && top_is_derma == 1 && ssciname !~ /Dermatophagoides farinae/) {
                top_is_derma = 0
            }
        }
    }
    END {
        # Check last qseqid
        if (top_is_derma == 1) {
            count++
        }
        print "'$sample'\t" count
    }
    ' tmp_sorted.txt >> hdm_read_counts.tsv

    rm tmp_sorted.txt
done

echo "Done. Output saved to hdm_read_counts.tsv"

