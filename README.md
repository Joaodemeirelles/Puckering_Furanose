# Puckering_Furanose.py
Script used to analyze puckering furanose conformation in Zx and Zy, P and Amp values in GROMACS

Output is various .csv values, being:

Zx.csv, Zy.csv, Zx_Zy.csv, P.csv, Amp.csv, P_Amp.csv

Usage: python Puckering_Furanose.py res_name_on_pdb res_id_number  your_sim.tpr your_sim.xtc

Example: python Puckering_Furanose.py 2DR 20 deoxyribose_100ns.tpr deoxyribose_100ns.xtc

# Dependencies

You need to install (pip or pip3 is a good way) the Biopython package (https://biopython.org/) and math.

