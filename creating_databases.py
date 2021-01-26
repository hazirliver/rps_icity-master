import os

# 1. Creating proteins database for rps-icity
os.system("makeblastdb -in ./Database/ProteinSequencesFasta.faa -out ./Database/ProteinsDB -dbtype prot -parse_seqids")

# 2. Creating rps-blast database for rps-icity
os.system('cd ./Database_pfam/SMPs; makeprofiledb -in Pfam_Cog.pn -title Cog_Pfam -out ../Cog_Pfam -dbtype "rps"')
