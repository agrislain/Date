---

# README - Date Pipeline

This pipeline automates a multi-step process to analyze protein sequences, extract taxonomic identifiers, determine the Most Recent Common Ancestor (MRCA), and refine the analysis using the NR database.

## Requirements

Ensure the following dependencies are available:

* 'diamond' (tested with version 2.0.13)
* 'python3'
* 'sbatch' (SLURM scheduler)
* A DIAMOND database ('.dmnd')
* An NR DIAMOND database (e.g., 'nr_2.0.13.dmnd')
* Python scripts located in the 'scripts/' directory

## Required Input Parameters

You must provide the following arguments when launching the pipeline:

* '--output_dir' : Output directory
* '--fasta_input' : Input protein sequences in FASTA format
* '--tree' : Newick-format phylogenetic tree file
* '--correspondance_file' : Taxonomic name mapping file
* '--database' : DIAMOND '.dmnd' database for initial BLAST search

Optional arguments:

* '--threads' : Number of threads to use (default: 16)
* '--evalue_filter' : E-value threshold for filtering BLAST hits (default: 1e-4)
* '--levels-up' : Number of tree levels to go up for MRCA determination (default: 1)

## Launching the Pipeline

'''
sbatch pipeline.sh \
  --output_dir /path/to/results \
  --fasta_input /path/to/sequences.faa \
  --tree /path/to/tree.newick \
  --correspondance_file /path/to/correspondences.txt \
  --database /path/to/database.dmnd \
  [--threads 16] \
  [--evalue_filter 1e-4] \
  [--levels-up 1]
'''

## Step-by-step Execution Details

### Step 0: DIAMOND BLAST against custom database

Command:

'''
diamond blastp \
  -d <database> \
  -q <fasta_input> \
  -o <output_dir>/<prefix>_blast_results.txt \
  -k0 \
  --sensitive \
  --threads <threads> \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
  --unal 1
'''

### Step 1: Extract taxonomic IDs from BLAST results

Command:

'''
python3 scripts/extract_taxid_optimized_parallel.py \
  -i <output_dir>/<prefix>_blast_results.txt \
  -o <output_dir>/<prefix>_taxids \
  -t <output_dir>/<prefix>_focale \
  --workers <threads> \
  --evalue_filter <evalue_filter>
'''

### Step 2: Find MRCA from taxonomic IDs

Command:

'''
python3 scripts/find_mrca.py \
  -t <tree> \
  -x <output_dir>/<prefix>_taxids \
  -o <output_dir>/<prefix>_mrca \
  -T <tree> \
  --focale_file <output_dir>/<prefix>_focale \
  -m species_percentage
'''

### Step 3: Generate NR subtree from MRCA

Command:

'''
python3 scripts/get_nr_subtree.py \
  -i <output_dir>/<prefix>_mrca \
  -o <output_dir>/<prefix>_old_mrca \
  -n <output_dir>/<prefix>_new_mrca \
  -t <tree> \
  -m <correspondance_file> \
  --levels-up <levels_up>
'''

### Step 4: DIAMOND BLAST against NR database

Command:

'''
python3 scripts/regroup_blast.py \
  -i <output_dir>/<prefix>_new_mrca \
  -f <fasta_input> \
  -b <diamond_bin> \
  -d <nr_database> \
  -o <output_dir>/<prefix>_blast_nr/
'''

### Step 5: Extract taxonomic IDs from NR BLAST

Command:

'''
python3 scripts/extract_taxid_optimized_parallel.py \
  -i <output_dir>/<prefix>_blast_nr/final_blast_results.txt \
  -o <output_dir>/<prefix>_hybride_taxids \
  -t <output_dir>/<prefix>_hybride_focale \
  --workers <threads> \
  --evalue_filter <evalue_filter>
'''

### Step 6: Determine MRCA from hybrid taxids

Command:

'''
python3 scripts/get_nr_mrca_parallel.py \
  -i <output_dir>/<prefix>_hybride_taxids \
  -o <output_dir>/<prefix>_hybride_mrca \
  --focale_file <output_dir>/<prefix>_hybride_focale \
  -t <threads>
'''

### Step 7: Merge old and hybrid MRCAs

Command:

'''
cat <output_dir>/<prefix>_old_mrca \
    <output_dir>/<prefix>_hybride_mrca \
    > <output_dir>/<prefix>_final_mrca
'''

## Output Summary

Key output files stored in '<output_dir>':

* '<prefix>_blast_results.txt' : Initial BLAST output
* '<prefix>_taxids', '<prefix>_focale' : Extracted taxonomic data
* '<prefix>_mrca', '<prefix>_old_mrca', '<prefix>_new_mrca' : MRCA intermediate results
* '<prefix>_blast_nr/' : Results from NR BLAST
* '<prefix>_hybride_taxids', '<prefix>_hybride_mrca' : Hybrid analysis outputs
* '<prefix>_final_mrca' : Final merged MRCA result

---
