#!/bin/bash

# ========================
# Argument processing
# ========================

threads=16  # Default value for threads
evalue_filter=1e-4  # Default value for evalue filter

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --output_dir) output_dir="$2"; shift ;;
        --fasta_input) fasta_input="$2"; shift ;;
        --tree) tree="$2"; shift ;;
        --correspondance_file) correspondance_file="$2"; shift ;;
        --database) database="$2"; shift ;;
        --threads) threads="$2"; shift ;;
        --levels-up) levels_up="$2"; shift ;;
        --evalue_filter) evalue_filter="$2"; shift ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
    shift
done

# Check mandatory parameters
if [[ -z "$output_dir" || -z "$fasta_input" || -z "$tree" || -z "$correspondance_file" || -z "$database" ]]; then
    echo "Usage: sbatch Date.sh --output_dir <path> --fasta_input <file> --tree <file> --correspondance_file <file> --database <file> [--threads <number>] [--evalue_filter <value>] [--levels-up <number>]" >&2
    exit 1
fi

# ========================
# Internal parameters
# ========================

diamond_bin="/store/EQUIPES/BIM/MEMBERS/antoine.grislain/diamond_2.0.13"
NR_db="/datas/NR/nr_2.0.13.dmnd"  # NR database used in step 4

mkdir -p "$output_dir"
short_name=$(basename "$output_dir")
levels_up=${levels_up:-1}  # Default value: 1

# ========================
# PIPELINE
# ========================

# Step 0: Run BLAST with the chosen database
job0=$(sbatch --parsable --time=182:00:00 --mem=90G --cpus-per-task=$threads --partition=bim --nodelist=node04 \
--wrap="/usr/bin/time -v $diamond_bin blastp \
    -d $database \
    -q $fasta_input \
    -o $output_dir/${short_name}_blast_results.txt \
    -k0 \
    --sensitive \
    --threads $threads \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
    --unal 1 \
    2> $output_dir/${short_name}_blast_timelog.txt")
echo "$job0" > "$output_dir/step0_blastp.jobid"

# Step 1
job1=$(sbatch --parsable --dependency=afterok:$job0 --time=50:00:00 --mem=90G --cpus-per-task=$threads \
--wrap="python3 scripts/extract_taxid_optimized_parallel.py \
-i $output_dir/${short_name}_blast_results.txt \
-o $output_dir/${short_name}_taxids \
-t $output_dir/${short_name}_focale \
--workers $threads \
--evalue_filter $evalue_filter")
echo "$job1" > "$output_dir/step1_extract_taxid_optimized_parallel.jobid"

# Step 2
job2=$(sbatch --parsable --dependency=afterok:$job1 --time=50:00:00 --mem=90G --cpus-per-task=$threads \
--wrap="python3 scripts/find_mrca.py \
-t $tree \
-x $output_dir/${short_name}_taxids \
-o $output_dir/${short_name}_mrca \
-T $tree \
--focale_file $output_dir/${short_name}_focale \
-m species_percentage")
echo "$job2" > "$output_dir/step2_find_mrca.jobid"

# Step 3
job3=$(sbatch --parsable --dependency=afterok:$job2 --time=50:00:00 --mem=90G --cpus-per-task=$threads \
--wrap="python3 scripts/get_nr_subtree.py \
-i $output_dir/${short_name}_mrca \
-o $output_dir/${short_name}_old_mrca \
-n $output_dir/${short_name}_new_mrca \
-t $tree \
-m $correspondance_file \
--levels-up $levels_up")
echo "$job3" > "$output_dir/step3_get_nr_subtree.jobid"

# Step 4
job4=$(sbatch --parsable --dependency=afterok:$job3 --time=50:00:00 --mem=90G --partition=bim --nodelist=node04 --cpus-per-task=$threads \
--wrap="python3 scripts/regroup_blast.py \
-i $output_dir/${short_name}_new_mrca \
-f $fasta_input \
-b $diamond_bin \
-d $NR_db \
-o $output_dir/${short_name}_blast_nr/")
echo "$job4" > "$output_dir/step4_regroup_blast.jobid"

# Step 5
job5=$(sbatch --parsable --dependency=afterok:$job4 --time=50:00:00 --mem=90G --cpus-per-task=$threads \
--wrap="python3 scripts/extract_taxid_optimized_parallel.py \
-i $output_dir/${short_name}_blast_nr/final_blast_results.txt \
-o $output_dir/${short_name}_hybride_taxids \
-t $output_dir/${short_name}_hybride_focale \
--workers $threads \
--evalue_filter $evalue_filter")
echo "$job5" > "$output_dir/step5_extract_taxid_optimized_parallel.jobid"

# Step 6
job6=$(sbatch --parsable --dependency=afterok:$job5 --time=50:00:00 --mem=90G --cpus-per-task=$threads \
--wrap="python3 scripts/get_nr_mrca_parallel.py \
-i $output_dir/${short_name}_hybride_taxids \
-o $output_dir/${short_name}_hybride_mrca \
--focale_file $output_dir/${short_name}_hybride_focale \
-t $threads")
echo "$job6" > "$output_dir/step6_get_nr_mrca_parallel.jobid"

# Step 7
job7=$(sbatch --parsable --dependency=afterok:$job6 --time=50:00:00 --mem=90G --cpus-per-task=$threads \
--wrap="python3 scripts/verification.py \
    --hybride $output_dir/${short_name}_hybride_mrca \
    --new $output_dir/${short_name}_new_mrca \
    -o $output_dir/${short_name}_mrca_verification \
    -t $threads")
echo "$job7" > "$output_dir/step7_verification.jobid"

# Step 8
job8=$(sbatch --dependency=afterok:$job7 --time=01:00:00 --mem=4G --cpus-per-task=1 \
--wrap="cat $output_dir/${short_name}_old_mrca \
$output_dir/${short_name}_hybride_mrca \
> $output_dir/${short_name}_final_mrca")
echo "$job8" > "$output_dir/step8_cat_final_mrca.jobid"
