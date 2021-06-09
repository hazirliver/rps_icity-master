import glob
import os
import argparse

ap = argparse.ArgumentParser(description = "Create icity CLUSTER_ID -- pfam names correspondence table")
ap.add_argument("-f", help = "Filenames suffix", required = True)

opts = ap.parse_args()
suffix = opts.f

GI_PFAM = {}
CLUSTER = {}

# Загружаем все именованные gi
with open('gi_pfam_named' + suffix + '.tsv', "r") as gi_named:
    for _ in range(1):
        next(gi_named)
    for line in gi_named:
        tmp = line.rstrip().split(",")
        GI_PFAM[tmp[0]] = tmp[1]


# Проходимся по всем .ali файлам, содержащим gi, относящиеся к одному кластеру
path_cl = './CLUSTERS' + suffix
for filename in glob.glob(os.path.join(path_cl, '*.ali')):
   with open(os.path.join(os.getcwd(), filename), 'r') as f:
       first_line = f.readline()
       gi = first_line.split("|")[1]

       Cluster_id = filename.rstrip(".ali").split("/")[-1]

       CLUSTER[Cluster_id] = GI_PFAM[gi]

# Запись полученных именованных кластеров в файлы

csv_file = 'clusters' + suffix + '.csv'
try:
    with open(csv_file, 'w') as csvfile:
        for key in CLUSTER.keys():
            csvfile.write("%s,%s\n" % (key, CLUSTER[key]))
except IOError:
    print("I/O error")

