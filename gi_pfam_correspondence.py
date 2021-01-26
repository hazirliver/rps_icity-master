import pandas as pd
import argparse

ap = argparse.ArgumentParser(description = "Create gi pfam correspondence table")
ap.add_argument("-f", help = "Filenames suffix", required = True)

opts = ap.parse_args()
suffix = opts.f

icity = pd.read_csv('gi_cluster_id' + suffix + '.csv', header = 0, names = ["icity_cluster_id","gi", "icity_score"], dtype=str)
rps = pd.read_csv('gi_pfam_named' + suffix + '.csv', header = 0, names = ["gi","pfam"], dtype=str)
vicinity = pd.read_csv('VicinityPermissiveClustsLinear' + suffix + '.tsv', header = None, names = ["gi","composition"], sep = "\t", dtype=str)

pfam = pd.merge(rps, vicinity['gi'], on = "gi", how = "inner")
pfam['idx'] = pfam.index
pfam = pfam.drop('gi', axis = 1)

#tmp = icity.drop_duplicates('icity_cluster_id').sort_values('icity_cluster_id')
tmp = icity
tmp[['c','idx']] = tmp.icity_cluster_id.str.split('_', expand=True)
tmp = tmp.drop(['c'], axis = 1)
tmp.idx = tmp.idx.astype('int')
tmp = tmp.sort_values('idx')
tmp.idx = tmp.idx - 1

cluster = pd.merge(tmp, pfam, on='idx').drop(['idx'], axis=1)
cluster.to_csv('icity_pfam' + suffix + '.csv', index = False)