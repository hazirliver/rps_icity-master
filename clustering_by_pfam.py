import os
import config
import argparse
import pandas as pd
from collections import defaultdict

ap = argparse.ArgumentParser(description = "Run RPS-BLAST")
ap.add_argument("-f", help = "Filenames suffix", required = True)

opts = ap.parse_args()
suffix = opts.f


def run_rps_blast(eval):
    command_rps_blast = 'rpsblast+ -query ./Vicinity' + suffix + '.faa -db ./Database_pfam/Cog_Pfam -out ./rps_blast_results' + suffix + '.tsv -evalue ' + str(eval) + ' -outfmt \"6 qacc stitle evalue\"'
    os.system(command_rps_blast)

def sort_rps_blast():
    df = pd.read_csv('./rps_blast_results' + suffix + '.tsv', sep="\t", header=None)

    first_column_tmp = (df[0].str.rstrip("|")).str.split("|")
    data = first_column_tmp.to_list()
    names = ["tmp1", "gi", "tmp2", "gb"]
    new_first_col = (pd.DataFrame(data, columns=names)).drop(['tmp1', 'tmp2'], axis=1)

    second_column_tmp = df[1].str.split(",", n=2)
    data = second_column_tmp.to_list()
    names = ["pfam_cog_id", "pfam_cog_smth_add", "pfam_cog_name"]
    new_sec_col = pd.DataFrame(data, columns=names)

    df_tmp = pd.concat([new_first_col, new_sec_col, df[2]], axis=1)
    df_tmp.columns = ["gi", "gb", "pfam_cog_id", "pfam_cog_smth_add", "pfam_cog_name", "eval"]
    df_new = df_tmp.loc[df_tmp.groupby('gi').eval.idxmin()]

    gi_pfam = df_new[["gi", "pfam_cog_id"]]
    gi_pfam.to_csv('gi_pfam_named' + suffix + '.csv', index=False, header=True)

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
    run_rps_blast(eval)
    sort_rps_blast()


clustering(config.ICITY_CONFIG_INPUT["RpsBlastEval"])

