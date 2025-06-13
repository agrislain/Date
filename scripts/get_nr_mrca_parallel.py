import argparse
from posixpath import basename
import subprocess
import os
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor
import sys

def clean_lineage_term(term):
    """Cleans and extracts the last term of a taxonomic lineage."""
    return term.split(";")[-1].strip()

@lru_cache(maxsize=None)
def cached_taxonkit_lca(taxids, threads=1):
    """Calls TaxonKit LCA with caching and thread support."""
    command = f"taxonkit lca --threads {threads}"
    result = subprocess.run(command, input=taxids, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"TaxonKit LCA error: {result.stderr.strip()}")
    return result.stdout.strip()

@lru_cache(maxsize=None)
def cached_taxonkit_lineage(taxid, threads=4):
    """Calls TaxonKit lineage with caching and thread support."""
    command = f"taxonkit lineage --threads {threads}"
    result = subprocess.run(command, input=taxid, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"TaxonKit lineage error: {result.stderr.strip()}")
    return result.stdout.strip()

def process_line(line, threads, focal_taxid):
    """Processes a line by calling TaxonKit and returns the result."""
    query, taxids = line.strip().split(";")
    # Add focal taxid to the taxids except for "no_hit" or "no_valid_taxid"
    if taxids not in ["no_hit", "no_valid_taxid"]:
        taxids = f"{taxids} {focal_taxid}"

    # Case where all taxids are "no_hit"
    if taxids == "no_hit":
        return f"{query}\tno MRCA\n"

    # Case where all taxids are "no_valid_taxid"
    if taxids == "no_valid_taxid":
        return f"{query}\tsynthetic construct or malformed\n"

    try:
        # Call TaxonKit for valid taxids
        mrca_result = cached_taxonkit_lca(taxids, threads)
        mrca_taxid = mrca_result.split("\t")[-1]
        lineage_result = cached_taxonkit_lineage(mrca_taxid, threads)
        lineage = lineage_result.split("\t")[-1]
        last_term = clean_lineage_term(lineage)
        return f"{query}\t{last_term}\t{mrca_taxid}\n"
    except Exception as e:
        print(f"Error for {query}: {e}")
        return f"{query}\tError\n"

def run_taxonkit_lca_with_last_lineage(input_file, output_file, threads=4, focal_taxid=None):
    """Applies TaxonKit LCA + lineage on each line with caching in parallel."""
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        lines = infile.readlines()

        # Use ThreadPoolExecutor to parallelize line processing
        with ThreadPoolExecutor(max_workers=threads) as executor:
            results = executor.map(lambda line: process_line(line, threads, focal_taxid), lines)

        # Write results to the output file
        outfile.writelines(results)

def main(input_file, focal_taxid, output_file, threads):
    run_taxonkit_lca_with_last_lineage(input_file, output_file, threads, focal_taxid)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Finds the MRCA for each query with a focal taxid and returns the last term of the taxonomic lineage.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Path to the input file (query -> taxids)")
    parser.add_argument("-f", "--focal_taxid", type=str, help="Taxid of the focal species to include")
    parser.add_argument("--focale_file", type=str, help="Path to a file containing the taxid of the focal species")
    parser.add_argument("-o", "--output", type=str, required=True, help="Path to the output file")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use for TaxonKit")
    args = parser.parse_args()

    # Handle focal taxid
    if args.focale_file:
        try:
            with open(args.focale_file, "r", encoding="utf-8") as f:
                focal_taxid = f.read().strip()
                print(f"Focal taxid loaded from file: {focal_taxid}")
        except FileNotFoundError:
            print(f"Error: File '{args.focale_file}' not found.")
            sys.exit(1)
    elif args.focal_taxid:
        focal_taxid = args.focal_taxid
    else:
        print("Error: You must provide either --focal_taxid or --focale_file.")
        sys.exit(1)

    main(args.input, focal_taxid, args.output, args.threads)
