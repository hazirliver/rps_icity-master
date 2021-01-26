import argparse
import pandas as pd

ap = argparse.ArgumentParser(description = "Rename original CLUSTER_ID to pfam names")
ap.add_argument("-f", help = "Filenames suffix", required = True)

opts = ap.parse_args()
suffix = opts.f


clusters = pd.read_csv('icity_pfam' + suffix + '.csv', header = 0)

clusters_dict = dict(zip(clusters.icity_cluster_id, clusters.pfam))

relevance = pd.read_csv('Relevance' + suffix + '.tsv', sep = "\t", header = None, names = ["Cluster_ID","e_size_baits", "e_size_database", "median_dist_to_baits", "icity"])

relevance.Cluster_ID = relevance.Cluster_ID.map(clusters_dict)
relevance.to_csv('relevance_pfam_named' + suffix + '.csv', sep = "\t", index=False)
