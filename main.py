import os
import argparse

ap = argparse.ArgumentParser(description = "Run main rps-blast clustering icity pipeline")
ap.add_argument("-s", help = "Seeds filename", required = True)
ap.add_argument("-f", help = "Output filenames suffix", required = True)

opts = ap.parse_args()
seeds = opts.s
suffix = opts.f
#
# 3. Selecting 10k neighborhoods
os.system('python3 SelectNeighborhood.py -p ./Database/CDS.pty -s ./' + seeds + ' -o Vicinity' + suffix + '.tsv -d 10000')

# 4. Getting protein IDs from bait neighborhoods
os.system('grep -v "===" Vicinity' + suffix + '.tsv | cut -f1 | sort -u > VicinityIDs' + suffix + '.lst')

# 4.1 Delete cr in file
os.system('sed -i "s/\r$//" VicinityIDs' + suffix + '.lst')

# 5. Getting protein sequences from bait neighborhoods
os.system('blastdbcmd -db ./Database/ProteinsDB -long_seqids > Vicinity' + suffix + '.faa -entry_batch VicinityIDs' + suffix + '.lst')

# 6. Run rps-blast clustering, sort and save results to VicinityPermissiveClustsLinear.tsv
os.system('python3 clustering_by_pfam.py -f ' + suffix)

# 7. Make profiles for the clusters
os.system('python3 ./MakeProfiles.py -f VicinityPermissiveClustsLinear' + suffix + '.tsv -c CLUSTERS' + suffix + '/ -d ./Database/ProteinsDB')

# 8. Running PSI-BLAST to search for generated protein profiles
os.system('python3 ./RunPSIBLAST.py -c CLUSTERS' + suffix + '/ -d ./Database/ProteinsDB')

# 9. Sorting of blast hits between clusters
os.system('python3 ./SortBLASTHitsInMemory.py -c CLUSTERS' + suffix + '/ -o CLUSTERS' + suffix + '/Sorted/ -p ./Database/CDS.pty -i VicinityIDs' + suffix + '.lst -s ' + seeds + ' -v Vicinity' + suffix + '.tsv -z 0.4 -x 0.35 -f ' + suffix)

# 10. Calculating Icity metric
os.system('bash ./CalculateICITY.sh ./CLUSTERS' + suffix + '/Sorted ./Database/ProteinsDB VicinityPermissiveClustsLinear' + suffix + '.tsv Relevance' + suffix + '.tsv')

# 11. Corresponding gi and icity ids clusters
os.system('python3 gi_icity_correspondence.py -f ' + suffix)

# 12. Corresponding gi and pfam_cog names
os.system('python3 gi_pfam_correspondence.py -f ' + suffix)

# 13. Replace icity clusters by pfam_cog names in relevance.tsv
os.system('python3 correcting_relevance.py -f ' + suffix)
