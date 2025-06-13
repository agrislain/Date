#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import argparse
import time
from multiprocessing import Pool
from collections import Counter, defaultdict


def detect_columns(input_file):
    """Detects the columns of the input file and returns the appropriate names."""
    standard_columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                        "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    nr_columns = standard_columns + ["staxids"]

    sample = pd.read_csv(input_file, sep="\t", nrows=1, header=None)
    if sample.shape[1] == len(standard_columns):
        return standard_columns, False
    elif sample.shape[1] == len(nr_columns):
        return nr_columns, True
    else:
        raise ValueError("Invalid input file. Expected 12 or 13 columns.")


def read_chunks_with_query_grouping(file_path, chunksize, columns, usecols, sep="\t"):
    """Reads a file in chunks while grouping rows associated with the same query."""
    leftover = pd.DataFrame()
    for chunk in pd.read_csv(file_path, sep=sep, chunksize=chunksize, header=None, names=columns):
        # Add incomplete rows from the previous chunk
        chunk = pd.concat([leftover, chunk], ignore_index=True)

        # Identify incomplete rows (split in the middle of a query)
        last_query = chunk["qseqid"].iloc[-1]
        incomplete_rows = chunk[chunk["qseqid"] == last_query]

        # Save these rows for the next chunk
        leftover = incomplete_rows

        # Remove incomplete rows from the current chunk
        chunk = chunk[chunk["qseqid"] != last_query]

        # Filter necessary columns
        chunk = chunk[usecols]

        yield chunk

    # Process remaining rows
    if not leftover.empty:
        yield leftover


def process_chunk(chunk, use_taxid, evalue_filter):
    """Processes a chunk and returns a dictionary of taxids by query."""
    query_to_taxids = {}

    # Convert the evalue column to float
    chunk["evalue"] = pd.to_numeric(chunk["evalue"], errors="coerce")

    # Identify rows without hits
    if use_taxid:
        no_hit_queries = chunk[
            (
                (chunk["staxids"].isna() | (chunk["staxids"] == "") | (chunk["staxids"] == "0")) &
                (chunk["evalue"] == -1)
            )
        ]["qseqid"].tolist()
    else:
        no_hit_queries = chunk[
            ((chunk["sseqid"] == "*") & (chunk["evalue"] == -1))
        ]["qseqid"].tolist()

    for query_id in no_hit_queries:
        query_to_taxids[query_id] = ["no_hit"]

    # Remove rows without hits
    if use_taxid:
        chunk = chunk[~(
            (chunk["staxids"].isna() | (chunk["staxids"] == "") | (chunk["staxids"] == "0")) &
            (chunk["evalue"] == -1)
        )]
    else:
        chunk = chunk[~((chunk["sseqid"] == "*") & (chunk["evalue"] == -1))]

    # Filter rows with valid evalues
    all_queries = set(chunk["qseqid"].unique())
    chunk = chunk[chunk["evalue"] < evalue_filter]
    remaining_queries = set(chunk["qseqid"].unique())
    no_hit_queries_after_filter = all_queries - remaining_queries

    for query_id in no_hit_queries_after_filter:
        query_to_taxids[query_id] = ["no_hit"]

    if chunk.empty:
        return query_to_taxids

    # Clean and standardize taxids/sseqid
    if use_taxid:
        chunk["staxids"] = chunk["staxids"].fillna("").astype(str).str.split(";")
        chunk["staxids"] = chunk["staxids"].apply(lambda x: [taxid for taxid in x if taxid != "32630"])
        chunk["staxids"] = chunk["staxids"].apply(lambda x: ";".join(x) if x else "no_valid_taxid")
    else:
        chunk["sseqid"] = chunk["sseqid"].fillna("").astype(str)
        chunk["sseqid"] = chunk["sseqid"].str.strip().str.replace(r"_[^_]*@", "_", regex=True)

    # Group by query_id and collect cleaned taxids
    grouped = chunk.groupby("qseqid")
    for query_id, group in grouped:
        if use_taxid:
            values = group["staxids"].tolist()
            flat_values = [item for sublist in values for item in sublist.split(";") if item.strip()]
        else:
            values = group["sseqid"].tolist()
            flat_values = [item for item in values if item.strip()]
        if flat_values:
            query_to_taxids[query_id] = list(set(flat_values))
        else:
            query_to_taxids[query_id] = ["no_valid_taxid"]

    return query_to_taxids


def extract_taxids_parallel(input_file, output_file, chunksize=100000, num_workers=4, evalue_filter=1e-4):
    """Extracts taxids in parallel."""
    columns, use_taxid = detect_columns(input_file)

    query_to_taxids = defaultdict(list)

    # Use multiprocessing to process chunks in parallel
    with Pool(processes=num_workers) as pool:
        tasks = []
        for chunk in read_chunks_with_query_grouping(
            file_path=input_file,
            chunksize=chunksize,
            columns=columns,
            usecols=["qseqid", "evalue", "staxids"] if use_taxid else ["qseqid", "sseqid", "evalue"],
            sep="\t",
        ):
            tasks.append(pool.apply_async(process_chunk, (chunk, use_taxid, evalue_filter)))

        # Collect results from processes
        for task in tasks:
            result = task.get()
            for query_id, taxids in result.items():
                query_to_taxids[query_id].extend(taxids)

    # Remove duplicates in taxids
    for query_id in query_to_taxids:
        query_to_taxids[query_id] = sorted(set(query_to_taxids[query_id]))

    # Write results to the output file
    with open(output_file, "w") as f:
        for k, v in query_to_taxids.items():
            f.write(f"{k};{' '.join(v)}\n")


def process_top_species(chunk, use_taxid):
    """Processes a chunk to calculate the most frequent species and collect species by query_id."""
    species_counter = Counter()
    query_to_species = {}

    # Standardize the species column
    if use_taxid:
        chunk["species"] = chunk["staxids"].astype(str).str.replace(r"(^|;)32630(;|$)", r"\1\2", regex=True)
    else:
        chunk["species"] = chunk["sseqid"].astype(str).str.replace(r"_[^_]*@", "_", regex=True)

    chunk["species"] = chunk["species"].fillna("unknown")

    # Keep only the first row per query_id
    chunk_first = chunk.drop_duplicates(subset="qseqid", keep="first")

    # Update the global counter and collect species by query_id
    for _, row in chunk_first.iterrows():
        query_id = row["qseqid"]
        species_list = [sp.strip() for sp in row["species"].split(";") if sp.strip()]
        species_counter.update(species_list)
        query_to_species[query_id] = species_list

    return species_counter, query_to_species


def get_top_species_parallel(input_file, chunksize=100000, num_workers=4):
    """Calculates the most frequent species in parallel."""
    columns, use_taxid = detect_columns(input_file)

    global_counter = Counter()
    query_to_species = {}

    # Use multiprocessing to process chunks in parallel
    with Pool(processes=num_workers) as pool:
        tasks = []
        for chunk in pd.read_csv(input_file, sep="\t", chunksize=chunksize, header=None, names=columns,
                                 usecols=["qseqid", "staxids"] if use_taxid else ["qseqid", "sseqid"]):
            tasks.append(pool.apply_async(process_top_species, (chunk, use_taxid)))

        # Collect results from processes
        for task in tasks:
            species_counter, chunk_query_to_species = task.get()
            global_counter.update(species_counter)
            query_to_species.update(chunk_query_to_species)

    # Identify the most frequent species
    most_common_species = global_counter.most_common(1)
    if not most_common_species:
        return "no_species_found"
    global_best = str(most_common_species[0][0])

    # Final calculation
    final_count = Counter()
    for query_id, species_list in query_to_species.items():
        if not species_list:
            continue  # Skip if species_list is empty
        # Ensure comparison is with string
        species_list_str = [str(s) for s in species_list]
        if global_best in species_list_str:
            final_count[global_best] += 1
        else:
            final_count[species_list_str[0]] += 1

    return final_count.most_common(1)[0][0]


def main(input_file, output_file, top_specie_path, chunksize=100000, num_workers=4, evalue_filter=1e-4):
    start = time.perf_counter()

    # Extract taxids
    extract_taxids_parallel(input_file, output_file, chunksize, num_workers, evalue_filter)
    print("Taxid extraction done")

    # Calculate the most frequent species
    top_species = get_top_species_parallel(input_file, chunksize, num_workers)
    with open(top_specie_path, "w") as f:
        f.write(top_species + "\n")
    print(f"Top species is: {top_species}")

    end = time.perf_counter()
    print(f"Elapsed time: {end - start:.2f}s")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parallelized taxid extraction and top species detection")
    parser.add_argument("-i", "--input", required=True, help="Input BLAST file")
    parser.add_argument("-o", "--output", required=True, help="Output file for taxids")
    parser.add_argument("-t", "--top_specie", required=True, help="Output file for top species")
    parser.add_argument("--chunksize", type=int, default=100000, help="Chunksize for reading")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers")
    parser.add_argument("--evalue_filter", type=float, default=1e-4, help="E-value filter threshold")

    args = parser.parse_args()
    main(args.input, args.output, args.top_specie, args.chunksize, args.workers, args.evalue_filter)
