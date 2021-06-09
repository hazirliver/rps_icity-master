import os
import config
import argparse
import pandas as pd
from collections import defaultdict

ap = argparse.ArgumentParser(description = "Run RPS-BLAST")
ap.add_argument("-f", help = "Filenames suffix", required = True)

opts = ap.parse_args()
suffix = opts.f


def run_hmmsearch(eval):
    command_hmmsearch = 'hmmsearch --tblout hmm_out' + suffix + '.txt -E 1e-5 --cpu 24 ./pfam/Pfam-A.hmm ./Vicinity' + suffix + '.faa'
    os.system(command_hmmsearch)

def sort_hmmsearch():

    os.system("Rscript --vanilla parse_hmmer3_tbloutput.R " + suffix)

    df_new = pd.read_csv('./gi_pfam_named' + suffix + '.tsv', sep="\t", header=0)

    pfam_cog_clusters_groups = df_new.groupby("pfam_cog_id")
    clusters = [pfam_cog_clusters_groups.get_group(x) for x in pfam_cog_clusters_groups.groups]

    gi_clustered = defaultdict(list)
    i = 0
    for x in clusters:
        for gi in x["gi"]:
            gi_clustered[i].append(gi)
        i += 1
    to_output = sorted(gi_clustered.values())

    with open('./VicinityPermissiveClustsLinear' + suffix + '.tsv', 'w') as out:
        for cluster in to_output:
            if len(cluster) > 1:
                j = 0
                for gi in cluster:
                    if j == 0:
                        out.write(str(gi) + "\t" + str(gi) + " ")
                        j += 1
                    else:
                        out.write(str(gi) + " ")
                out.write("\n")

def clustering(eval):
    run_hmmsearch(eval)
    sort_hmmsearch()


clustering(config.ICITY_CONFIG_INPUT["RpsBlastEval"])

