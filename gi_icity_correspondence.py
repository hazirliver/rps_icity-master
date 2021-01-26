import argparse
import pandas as pd

ap = argparse.ArgumentParser(description = "Creating gi-icity correspondence table")
ap.add_argument("-f", help = "Filenames suffix", required = True)

opts = ap.parse_args()
suffix = opts.f

relevance_dict = {}
with open('Relevance' + suffix + '.tsv') as relevance:
    for line in relevance:
        line_splited = line.split("\t")
        relevance_dict[line_splited[0]] = line_splited[4].rstrip()

    relevance = pd.DataFrame(relevance_dict.items(), columns = ['icity_cluster_id', 'icity_score'])

cluster_gis = pd.DataFrame(columns=['icity_cluster_id', 'gi'])
i = 1
while True:
    try:
        with open('./CLUSTERS' + suffix + '/Sorted/CLUSTER_' + str(i) + ".hits_sorted", "r") as hits_sorted:
            for line in hits_sorted:
                cluster_gis = cluster_gis.append({'icity_cluster_id': "CLUSTER_" + str(i), 'gi': line.split("\t")[0]}, ignore_index=True)
            i += 1
    except FileNotFoundError:
            break

to_out = cluster_gis.merge(relevance, on = 'icity_cluster_id', how = 'left')


to_out.to_csv('gi_cluster_id' + suffix + '.csv', index = False)