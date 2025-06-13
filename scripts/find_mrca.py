import argparse
from ete3 import Tree
import csv
import sys
import time
from tqdm import tqdm

csv.field_size_limit(sys.maxsize)

def load_taxid_mapping(taxid_file):
    """Loads the taxid mapping from a file."""
    taxid_mapping = {}
    with open(taxid_file, newline='', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile, delimiter=';')
        for row in reader:
            if len(row) >= 2:
                qseqid = row[0]
                taxids = row[1].split()
                taxid_mapping[qseqid] = taxids
    return taxid_mapping

def assign_internal_node_names(tree):
    """Assigns names to internal nodes if they are not already named."""
    all_annotated = all(node.name for node in tree.traverse() if not node.is_leaf())
    if all_annotated:
        print("All internal nodes are already annotated. No need to rewrite the tree.")
        return False  # Indicates no modifications were made

    counter = 1
    for node in tree.traverse():
        if not node.name:
            node.name = f"Internal_{counter}"
            counter += 1
    return True  # Indicates modifications were made

def find_mrca(tree_file, taxid_file, output_file, output_tree, anchor_specie, method):
    """Finds the Most Recent Common Ancestor (MRCA) for query sequences."""
    tree = Tree(tree_file, format=1)
    modified = assign_internal_node_names(tree)
    taxid_mapping = load_taxid_mapping(taxid_file)
    
    taxon_nodes = {}
    for node in tree.traverse():
        if node.name:
            taxon_id = node.name 
            if taxon_id not in taxon_nodes:
                taxon_nodes[taxon_id] = []
            taxon_nodes[taxon_id].append(node)
    
    results = []
    mrca_count = {}
    anchor_node = tree & anchor_specie  # Get the anchor node for the focal species
    
    start_time = time.time()
    
    for qseqid, taxids in tqdm(taxid_mapping.items(), desc="Processing qseqid", unit="qseqid"):
        # Filter valid taxids
        valid_taxids = [taxid for taxid in taxids if taxid not in {"no_hit", "no_valid_taxid"}]
        
        if not valid_taxids:
            # Handle cases with no valid taxids
            if all(taxid == "no_hit" for taxid in taxids):
                mrca = anchor_node
                metric = "Assigned to focal species (no_hit)"
                results.append([qseqid, mrca.name, metric])
                continue
            
            if all(taxid == "no_valid_taxid" for taxid in taxids):
                mrca = anchor_node
                metric = "Assigned to focal species (no_valid_taxid)"
                results.append([qseqid, mrca.name, metric])
                print(f"Warning: Query {qseqid} hit a synthetic construct or was malformed.")
                continue
        
        taxa = [taxon_nodes[tid] for tid in valid_taxids if tid in taxon_nodes]
        taxa.append([anchor_node])
        
        if taxa:
            flat_taxa = [t for sublist in taxa for t in sublist]
            mrca = tree.get_common_ancestor(flat_taxa)
            
            # Populate observed leaves
            observed_leaves = set()
            for tid in valid_taxids:
                if tid in taxon_nodes:
                    for node in taxon_nodes[tid]:
                        observed_leaves.update(node.get_leaf_names())
                else:
                    print(f"Warning: taxid {tid} not found in taxon_nodes")
            
            # Calculate metrics based on the chosen method
            if method == "species_percentage":
                leaves_under_mrca = set(mrca.get_leaf_names())
                observed_count = len(leaves_under_mrca.intersection(observed_leaves))
                total_count = len(leaves_under_mrca)
                coverage_percentage = (observed_count / total_count) if total_count > 0 else 0
                metric = f"{coverage_percentage:.2f}"
            
            elif method == "node_activation":
                leaves_seen = set()
                path_to_mrca = []
                current_node = anchor_node
                while current_node != mrca:
                    path_to_mrca.append(current_node)
                    current_node = current_node.up
                path_to_mrca.append(mrca)
                activation_vector = []
                
                for node in path_to_mrca:
                    leaves_under_node = set(node.get_leaf_names())
                    new_species = leaves_under_node - leaves_seen
                    observed_in_new_species = len(new_species.intersection(observed_leaves))
                    total_new_species = len(new_species)
                    proportion_new_species = (observed_in_new_species / total_new_species) if total_new_species > 0 else 0
                    activation_vector.append(proportion_new_species)
                    leaves_seen.update(leaves_under_node)
                
                metric = activation_vector
        
        else:
            mrca = None
            metric = "Not found" if method == "species_percentage" else []
        
        results.append([qseqid, mrca.name if mrca else "Not found", metric])
        
        if mrca and mrca.name not in mrca_count:
            mrca_count[mrca.name] = 0
        if mrca:
            mrca_count[mrca.name] += 1
    
    total = len(results)
    most_found = max(mrca_count, key=mrca_count.get) if mrca_count else "None"
    least_found = min(mrca_count, key=mrca_count.get) if mrca_count else "None"
    not_found_count = len([r for r in results if r[1] == "Not found"])
    
    mrca_percentages = {mrca: (count / total) * 100 for mrca, count in mrca_count.items()}
    
    print(f"Found MRCA for {total} query proteins")
    print(f"MRCA most found: {most_found}")
    print(f"MRCA least found: {least_found}")
    print(f"MRCA not found: {not_found_count}")
    print("MRCA percentages:")
    for mrca, percentage in mrca_percentages.items():
        print(f"  {mrca}: {percentage:.2f}%")
    
    if modified:
        tree.write(format=1, format_root_node=True, outfile=output_tree)
    else:
        print("The tree was not modified and will not be rewritten.")
    
    with open(output_file, "w", newline='') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        if method == "species_percentage":
            writer.writerow(["qseqid", "mrca", "coverage_percentage"])
        elif method == "node_activation":
            writer.writerow(["qseqid", "mrca", "activation_vector"])
        writer.writerows(results)
    
    end_time = time.time()
    print(f"Execution time: {end_time - start_time:.2f} seconds")
    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find the Most Recent Common Ancestor (MRCA) for query sequences.")
    parser.add_argument("-t", "--tree_file", required=True, help="Path to the tree file.")
    parser.add_argument("-x", "--taxid_file", required=True, help="Path to the taxid mapping file.")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output file.")
    parser.add_argument("-T", "--output_tree", required=True, help="Path to the output tree file.")
    parser.add_argument("-f", "--focale", help="Anchor species name.")
    parser.add_argument("--focale_file", help="Path to a file containing the anchor species name.")
    parser.add_argument("-m", "--method", required=True, choices=["species_percentage", "node_activation"], help="Method to determine node activation.")

    args = parser.parse_args()

    if args.focale_file:
        try:
            with open(args.focale_file, "r", encoding="utf-8") as f:
                focale = f.read().strip()
                print(f"Focal species loaded from file: {focale}")
        except FileNotFoundError:
            print(f"Error: Focal species file '{args.focale_file}' not found.")
            sys.exit(1)
    elif args.focale:
        focale = args.focale
    else:
        print("Error: You must provide either --focale or --focale_file.")
        sys.exit(1)

    find_mrca(
        tree_file=args.tree_file,
        taxid_file=args.taxid_file,
        output_file=args.output_file,
        output_tree=args.output_tree,
        anchor_specie=focale,
        method=args.method
    )
