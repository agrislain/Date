import argparse
from ete3 import Tree, NCBITaxa
from functools import lru_cache
import subprocess
import re
import csv


def clean_taxon_name(name):
    """Cleans and standardizes taxon names."""
    name = name.split("_", 1)[-1].replace("_", " ")
    name = re.sub(r'\b(sp|Colp-|Cur-|TD-)\b.*', '', name, flags=re.IGNORECASE).strip()
    return name

def load_mapping(mapping_file):
    """Loads the MRCA-to-NR mapping from a file."""
    mapping = {}
    with open(mapping_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 3:
                continue
            mrca_key, taxon_name, taxid = parts
            mapping[mrca_key] = (taxon_name, taxid)
    return mapping

def process_old_mrca(qseqid, mrca_name, mrca_node, old_outfile, mapping_dict):
    """Processes old MRCA entries."""
    try:
        if mrca_name not in mapping_dict:
            return f"{qseqid}\t{mrca_name}\tNot found in mapping\n"
        taxon_name, _ = mapping_dict[mrca_name]
        old_outfile.write(f"{qseqid}\t{taxon_name}\n")
    except Exception as e:
        return f"Error processing old MRCA {mrca_name}: {e}\n"

def process_new_mrca(qseqid, mrca_node, include_topology, new_outfile, mapping_dict):
    """Processes new MRCA entries."""
    try:
        mrca_key = mrca_node.name
        if mrca_key not in mapping_dict:
            return f"{qseqid}\t{mrca_key}\tNot found in mapping\n"

        taxon_name, taxid = mapping_dict[mrca_key]
        if include_topology:
            new_outfile.write(f"{qseqid}\t{taxon_name}\t{taxid}\tNA\n")
        else:
            new_outfile.write(f"{qseqid}\t{taxon_name}\t{taxid}\n")
    except Exception as e:
        return f"Error processing new MRCA: {e}\n"

def process_line(line, original_nodes, ncbi, include_topology, old_outfile, new_outfile, mapping_dict, levels_up):
    """Processes a single line from the input file."""
    try:
        parts = line.strip().split("\t")
        if len(parts) != 3:
            return f"Malformed line: {line.strip()}\n"

        qseqid, mrca_name, _ = parts
        mrca_node = original_nodes.get(mrca_name)

        if not mrca_node:
            return f"MRCA {mrca_name} not found in the tree.\n"

        # Traverse up the specified number of levels
        current_node = mrca_node
        for _ in range(levels_up): 
            if not current_node.up or not current_node.up.name:
                # If unable to traverse further, treat as an old MRCA
                return process_old_mrca(qseqid, mrca_name, mrca_node, old_outfile, mapping_dict)
            current_node = current_node.up

        mrca_plusN_key = current_node.name

        # Check if MRCA+N exists in the mapping
        if mrca_plusN_key in mapping_dict:
            mrca_plusN_taxon, _ = mapping_dict[mrca_plusN_key]
            # If MRCA+N corresponds to "Eukaryota", treat as an old MRCA
            if mrca_plusN_taxon == "Eukaryota":
                return process_old_mrca(qseqid, mrca_name, mrca_node, old_outfile, mapping_dict)

        # Otherwise, treat as a new MRCA
        return process_new_mrca(qseqid, current_node, include_topology, new_outfile, mapping_dict)
    except Exception as e:
        return f"Error processing line {line.strip()}: {e}\n"

def get_nr_subtree_and_species(input_file, original_tree_file, old_output_file, new_output_file, include_topology, mapping_file, levels_up):
    """Processes the input file to associate MRCA+N with NR using a mapping file."""
    try:
        tree = Tree(original_tree_file, format=1)
        original_nodes = {node.name: node for node in tree.traverse()}
        ncbi = NCBITaxa()
        mapping_dict = load_mapping(mapping_file)

        with open(input_file, "r") as f, \
             open(old_output_file, "w") as old_out_f, \
             open(new_output_file, "w") as new_out_f:
            next(f)  # Skip header
            for line in f:
                result = process_line(line, original_nodes, ncbi, include_topology, old_out_f, new_out_f, mapping_dict, levels_up)
                if result:
                    print(result.strip())
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Associates MRCA+N with NR using a mapping file.")
    parser.add_argument("-i", "--input", required=True, help="Input file with qseqid and MRCA")
    parser.add_argument("-t", "--tree", required=True, help="Newick file containing the original tree")
    parser.add_argument("-o", "--old-output", required=True, help="Output file for old taxa")
    parser.add_argument("-n", "--new-output", required=True, help="Output file for new taxa")
    parser.add_argument("-m", "--mapping-file", required=True, help="Mapping file for MRCA->NR")
    parser.add_argument("--include-topology", action="store_true", help="Include topology in the results (placeholder NA)")
    parser.add_argument("--levels-up", type=int, default=1, help="Number of levels to traverse up the tree (default: 1)")

    args = parser.parse_args()

    try:
        get_nr_subtree_and_species(
            args.input, args.tree, args.old_output, args.new_output,
            args.include_topology, args.mapping_file, args.levels_up
        )
    except Exception as e:
        print(f"Error: {e}")
