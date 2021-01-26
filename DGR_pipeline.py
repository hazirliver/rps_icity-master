import os

# 14. Run main 3 tiems for RT, TR, VR
os.system('python3 main.py -s seeds_RT.tsv -f _RT')
os.system('python3 main.py -s seeds_TR.tsv -f _TR')
os.system('python3 main.py -s seeds_VR.tsv -f _VR')

# 15. Move files to corresponding folders
os.system('move_dgr.sh')