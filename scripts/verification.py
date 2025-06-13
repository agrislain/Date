import pandas as pd
import subprocess
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed

def get_taxid(name):
    """Convert a name to taxid using TaxonKit."""
    result = subprocess.run(
        f'echo "{name}" | taxonkit name2taxid',
        capture_output=True,
        text=True,
        shell=True
    )
    taxid = result.stdout.strip().split("\t")[-1]
    return taxid if taxid.isdigit() else None

def get_lineage(taxid):
    """Get the lineage string for a taxid using TaxonKit."""
    result = subprocess.run(
        f'echo "{taxid}" | taxonkit lineage',
        capture_output=True,
        text=True,
        shell=True
    )
    if result.returncode != 0:
        return None
    parts = result.stdout.strip().split("\t")
    return parts[-1] if len(parts) >= 2 else None

def batch_call(func, items, n_threads=8):
    """Run a function in parallel using ThreadPoolExecutor."""
    results = {}
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        future_to_item = {executor.submit(func, item): item for item in items}
        for future in as_completed(future_to_item):
            item = future_to_item[future]
            try:
                results[item] = future.result()
            except Exception:
                results[item] = None
    return results

def main(hybride_file, new_file, output_file, n_threads=8):
    # Adapter la lecture du fichier hybride
    df1 = pd.read_csv(hybride_file, names=["query", "mrca_hybride", "taxid_hybride"], sep="\t")
    df2 = pd.read_csv(new_file, names=["query", "mrca_new", "taxid_new"], sep="\t")
    df = pd.merge(df1, df2, on="query")

    # Nettoyage des taxids : conversion en int puis str pour éviter les .0
    df["taxid_hybride"] = pd.to_numeric(df["taxid_hybride"], errors="coerce").dropna().astype(int).astype(str)
    df["taxid_new"] = pd.to_numeric(df["taxid_new"], errors="coerce").dropna().astype(int).astype(str)

    # On ne garde que les taxids numériques valides
    all_taxids = pd.concat([df["taxid_hybride"], df["taxid_new"]]).dropna().unique()
    all_taxids = [t for t in all_taxids if t.isdigit()]

    taxid2lineage = batch_call(get_lineage, all_taxids, n_threads=n_threads)
    df["lineage_hybride"] = df["taxid_hybride"].map(taxid2lineage)
    df["lineage_new"] = df["taxid_new"].map(taxid2lineage)

    results = []
    for _, row in df.iterrows():
        taxid_hybrid = row["taxid_hybride"]
        taxid_new = row["taxid_new"]
        lineage_hybrid = row["lineage_hybride"]
        lineage_new = row["lineage_new"]
        # Vérifie que les lineages sont bien des chaînes non nulles
        if not taxid_hybrid or not taxid_new:
            continue
        if not isinstance(lineage_hybrid, str) or not isinstance(lineage_new, str):
            continue
        lineage_hybrid_list = lineage_hybrid.split(";")
        lineage_new_list = lineage_new.split(";")
        if len(lineage_hybrid_list) < len(lineage_new_list) or str(lineage_new_list[-1]) not in lineage_hybrid_list:
            results.append({
                "query": row["query"],
                "mrca_hybride": row["mrca_hybride"],
                "mrca_new": row["mrca_new"],
                "taxid_hybride": taxid_hybrid,
                "taxid_new": taxid_new,
                "lineage_hybride": lineage_hybrid,
                "lineage_new": lineage_new
            })
    pd.DataFrame(results).to_csv(output_file, sep="\t", index=False)
    print(f"{len(results)} queries where hybrid MRCA is older than new MRCA written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check that hybrid MRCA is younger (more specific) than new MRCA using TaxonKit lineages (parallelized).")
    parser.add_argument("--hybride", required=True, help="Hybride MRCA file (query<TAB>mrca)")
    parser.add_argument("--new", required=True, help="New NR MRCA file (query<TAB>mrca)")
    parser.add_argument("-o", "--output", required=True, help="Output file for queries where hybrid MRCA is older")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads for parallelization")
    args = parser.parse_args()
    main(args.hybride, args.new, args.output, n_threads=args.threads)
