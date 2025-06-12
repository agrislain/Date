sbatch jobs/Date.sh --output_dir <chemin_du_répertoire_de_sortie> --fasta_input <chemin_du_fichier_fasta> --tree <chemin_du_fichier_arbre> --correspondance_file <chemin_du_fichier_correspondance> --database <chemin_de_la_base_de_données> [--threads <nombre_de_threads>]

------------------------------------------------------------------------------------------------------------------------------------------

sbatch --time=50:00:00 --mem=90G --cpus-per-task=8 --wrap="python3 scripts/extract_taxid_optimized_parallel.py -i diamond_blast_22genomes_vs_264_sensitive/Dmel_CDS_vs_264_sensitive.txt -o drosophila_melanogaster_264/drosophila264_taxids -t drosophila_melanogaster_264/drosophila264_focale --workers 8 --evalue_filter 1e-3"

sbatch --time=50:00:00 --mem=90G --cpus-per-task=8 --wrap="python3 scripts/find_mrca.py -t Dmel_264_sensitive/264_taxa_excavata-rooted_Dmel_internal_nodes.nw -x drosophila_melanogaster_264/drosophila264_taxids -o drosophila_melanogaster_264/drosophila264_mrca -T drosophila_melanogaster_264/264_taxa_excavata-rooted_Dmel_internal_nodes.nw -f EP00099_Drosophila_melanogaster -m species_percentage "

sbatch --time=50:00:00 --mem=90G --cpus-per-task=8 --wrap="python3 scripts/get_nr_subtree.py -i drosophila_melanogaster_264/drosophila264_mrca -o drosophila_melanogaster_264/drosophila264_old_mrca -n drosophila_melanogaster_264/drosophila264_new_mrca -t Dmel_264_sensitive/264_taxa_excavata-rooted_Dmel_internal_nodes.nw -m correspondances/correspondance_264_nr_nodes.tsv "

sbatch --time=50:00:00 --mem=90G --partition=bim --nodelist=node04 --cpus-per-task=8 --wrap="python3 scripts/regroup_blast.py -i drosophila_melanogaster_264/drosophila264_new_mrca -f 22_proteomes_modeles/Dmel_CDS.faa -b /store/USERS/antoine.grislain/diamond_2.0.13 -d /datas/NR/nr_2.0.13.dmnd -o drosophila_melanogaster_264/drosophila264_blast_nr/ "

sbatch --time=50:00:00 --mem=90G --cpus-per-task=8 --wrap="python3 scripts/extract_taxid_optimized_parallel.py -i drosophila_melanogaster_264/drosophila264_blast_nr/final_blast_results.txt -o drosophila_melanogaster_264/drosophila264_hybride_taxids -t drosophila_melanogaster_264/drosophila264_hybride_focale --workers 8 "

sbatch --time=50:00:00 --mem=90G --cpus-per-task=8 --wrap="python3 scripts/get_nr_mrca_parallel.py -i drosophila_melanogaster_264/drosophila264_hybride_taxids -o drosophila_melanogaster_264/drosophila264_hybride_mrca -f 7227 -t 8"

cat drosophila_melanogaster_264/drosophila264_old_mrca drosophila_melanogaster_264/drosophila264_hybride_mrca > drosophila_melanogaster_264/drosophila264_final_mrca