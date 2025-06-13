import argparse
import os
import subprocess
import tempfile

def parse_get_nr_subtree_output(input_file, output_dir):
    """
    Reads the output file from get_nr_subtree line by line and stores results
    in temporary files per MRCA to avoid memory overload.
    """
    os.makedirs(output_dir, exist_ok=True)
    mrca_files = {}

    with open(input_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue  # Skip malformed lines
            qseqid, mrca_name, mrca_taxid = parts[:3]

            # Replace spaces in the MRCA name with underscores
            mrca_safe = mrca_name.replace(" ", "_")

            # Create a file per MRCA if it doesn't already exist
            if mrca_safe not in mrca_files:
                mrca_files[mrca_safe] = open(os.path.join(output_dir, f"mrca_{mrca_safe}.txt"), "w")

            # Write data to the corresponding file
            mrca_files[mrca_safe].write(f"{qseqid}\t{mrca_taxid}\n")

    # Close all open files
    for f in mrca_files.values():
        f.close()

    return list(mrca_files.keys())  # Return the list of processed MRCAs

def get_fasta_sequence(fasta_file, query_id):
    """
    Searches for a specific FASTA sequence in a file without loading the entire file into memory.
    """
    with open(fasta_file, "r") as f:
        record = None
        for line in f:
            if line.startswith(">"):  # Start of a new sequence
                record_id = line.strip()[1:]
                if record_id == query_id:
                    record = ""
                else:
                    if record:
                        return record  # Return the found sequence
                    record = None
            elif record is not None:
                record += line.strip()
        return record  # Return the last found sequence

def run_diamond_blast(mrca, mrca_file, fasta_file, diamond_bin, diamond_db, output_dir):
    """
    Runs a DIAMOND BLAST for a specific MRCA by reading FASTA sequences on demand.
    """
    # Replace spaces in the MRCA name with underscores
    mrca_safe = mrca.replace(" ", "_")

    # Create a temporary FASTA file for the queries
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as fasta_temp:
        with open(mrca_file, "r") as f:
            for line in f:
                query_id, taxon_id = line.strip().split("\t")
                sequence = get_fasta_sequence(fasta_file, query_id)
                if sequence:
                    fasta_temp.write(f">{query_id}\n{sequence}\n")
                else:
                    print(f"Warning: Query ID {query_id} not found in FASTA file.")
        fasta_temp.flush()

    # Run DIAMOND
    output_file = os.path.join(output_dir, f"blast_results_mrca_{mrca_safe}.txt")
    command = (
        f"{diamond_bin} blastp -q {fasta_temp.name} -d {diamond_db} "
        f"--taxonlist {taxon_id} -o \"{output_file}\" -k0 --sensitive --outfmt 6 "
        f"qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids "
        f"--unal 1 "
    )

    print(f"Running DIAMOND for MRCA {mrca_safe}: {command}")
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"DIAMOND error for MRCA {mrca_safe}: {e}")

    # Remove the temporary FASTA file
    os.remove(fasta_temp.name)

def concatenate_blast_outputs(output_dir, final_output_file):
    """
    Combines all DIAMOND output files into a single file.
    """
    with open(final_output_file, "w") as outfile:
        for file_name in sorted(os.listdir(output_dir)):
            if file_name.startswith("blast_results_mrca_") and file_name.endswith(".txt"):
                with open(os.path.join(output_dir, file_name), "r") as infile:
                    outfile.write(infile.read())

    print(f"All results have been combined into {final_output_file}")

def main(input_file, fasta_file, diamond_bin, diamond_db, output_dir):
    """
    Main function: parses the input file, runs DIAMOND sequentially,
    and combines the results.
    """
    os.makedirs(output_dir, exist_ok=True)
    print("Parsing the get_nr_subtree output file...")
    mrca_list = parse_get_nr_subtree_output(input_file, output_dir)

    print("Running DIAMOND BLAST sequentially...")

    for mrca in mrca_list:
        mrca_file = os.path.join(output_dir, f"mrca_{mrca}.txt")
        run_diamond_blast(mrca, mrca_file, fasta_file, diamond_bin, diamond_db, output_dir)

    print("Combining output files...")
    final_output_file = os.path.join(output_dir, "final_blast_results.txt")
    concatenate_blast_outputs(output_dir, final_output_file)

    print("Removing individual DIAMOND output files...")
    for file_name in os.listdir(output_dir):
        if file_name.startswith("blast_results_mrca_") and file_name.endswith(".txt"):
            os.remove(os.path.join(output_dir, file_name))

    print("Removing temporary MRCA files...")
    for file_name in os.listdir(output_dir):
        if file_name.startswith("mrca_") and file_name.endswith(".txt"):
            os.remove(os.path.join(output_dir, file_name))

    print("Processing complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Groups queries by MRCA and runs DIAMOND BLAST.")
    parser.add_argument("-i", "--input", required=True, help="Output file from get_nr_subtree")
    parser.add_argument("-f", "--fasta", required=True, help="FASTA file containing the queries")
    parser.add_argument("-b", "--diamond-bin", required=True, help="Path to the DIAMOND executable")
    parser.add_argument("-d", "--diamond-db", required=True, help="DIAMOND database")
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory for results")

    args = parser.parse_args()

    try:
        main(args.input, args.fasta, args.diamond_bin, args.diamond_db, args.output_dir)
    except Exception as e:
        print(f"Error: {e}")
