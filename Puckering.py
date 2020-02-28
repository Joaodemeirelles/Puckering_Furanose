# Script used to analyze puckering furanose conformation in Zx and Zy, P and Amp values
# Output is various .csv values, being:
# Zx.csv, Zy.csv, Zx_Zy.csv, P.csv, Amp.csv, P_Amp.csv 

# You need to install Biopython and math python packages

# Usage: python Puckering_Furanose.py res_name_on_pdb res_id_number_on_pdb  your_sim.tpr your_sim.xtc
# Example: python Puckering_Furanose.py 2DR 20 deoxyribose_100ns.tpr deoxyribose_100ns.xtc 


#########################################################   PACKAGES NEEDED     #########################################################

import os
import sys
from collections import defaultdict
from Bio.PDB import *
import math


######################################################### 	FUNCTIONS  	#########################################################

def dumping_frames(res_name, tpr, xtc):

	if not os.path.exists("frames/"):
		os.system("mkdir frames")
	else:
		os.system("rm frames/*")

	os.system("echo {} System | gmx_51X trjconv -s " + tpr + " -f " + xtc + " -sep -skip " + str(skip) + " -o frames/frame_.pdb -pbc mol -center".format(res_name))


##########################################################################################################################################

res_name = sys.argv[1]		# Residue name is the first input of the script
res_id = int(sys.argv[2])	# Second input is residue identifier
tpr = sys.argv[3]		# Third, dynamics tpr
xtc = sys.argv[4]		# Fourth, dynamics xtc
skip = 2
dt = 2*skip			# Definition of timestep, depends on the simulation you're running


#dumping_frames(res_name,tpr, xtc)


current_dir = os.getcwd() # Takes current directory
pdbs = os.listdir(current_dir + "/frames")	# Takes frames directory, saving all pdb files
os.chdir(current_dir + "/frames")

Zx_all = []
Zy_all = []
Ptheta_all = []
Ar_all = []	# These are the list that will store values of Zx, Zy, Ptheta and Ar of all frames

arq = open("Zx.csv","w")
arq.close()
arq = open("Zy.csv","w")
arq.close()
arq = open("Zx_Zy.csv","w")
arq.close()
arq = open("P.csv","w")
arq.close()
arq = open("Amp.csv","w")
arq.close()
pi = math.pi
z =  4*pi/5
zx_1 = 2*math.cos(z) * 57.29
zy_1 = 2*math.sin(z) * 57.29


for i in range(len(pdbs)):
	
    item = pdbs[i]
	
    if ".csv" in item:		# Just in the case there is a left .csv in frames
        continue

    struct = PDBParser() # Creates a pdb object
    structure = struct.get_structure(item,item)
    residue_list = structure.get_residues()
	
    for residue in residue_list:
        	 # In this first part, we test if the atom is from the same residue we want
                

        if residue.get_resname() != res_name:   # If this atom residue name is not the same you want: jump this atom in the loop
            continue
                
        if residue.get_id()[1] != res_id:       # If this atom res_id is different from the same you want: jump this atom in the loop
            continue
	
        atom_list = residue.get_atoms()

        for atom in atom_list:
            if(atom.get_name() == "C1"):

                vectorC1 = atom.get_vector()	# Get vectors for calculating dihedral

            if(atom.get_name() == "C2"):

                vectorC2 = atom.get_vector()	# Get vectors for calculating dihedral

            if(atom.get_name() == "C3"):

                vectorC3 = atom.get_vector()

            if(atom.get_name() == "C4"):

                vectorC4 = atom.get_vector()
		
            if(atom.get_name() == "O4"):

                vectorO4 = atom.get_vector()



    # Calculate torsion dihedrals v1 and v2 to calculate Zx and Zy
    v1 = calc_dihedral(vectorO4,vectorC1,vectorC2,vectorC3) * 57.29 
    v3 = calc_dihedral(vectorO4,vectorC4,vectorC3,vectorC2) * 57.29
    # Calculate Zx and Zy based on v1 and v3

    Zx = ((v1 + v3) / zx_1) * 57.29
    Zy = ((v1 - v3) / zy_1) * 57.29  
    # Now, calculate Ptheta and Ar based on Zx and Zy
    Ptheta = 57.296 * math.atan(Zy/Zx)

    Ar = math.sqrt((Zy**2) + (Zx**2)) * 57.269
	
    # We add this values to it's respective list
	
    Zx_all.append(Zx)
    Zy_all.append(Zy)
    Ptheta_all.append(Ptheta)
    Ar_all.append(Ar)



with open("Zx.csv","a") as file:
    file.write("# Frame, Zx (Degrees)\n")
    for i in range(len(Zx_all)):
        file.write("{}\t{}\n".format(i,Zx_all[i]))

with open("Zy.csv","a") as file:
    file.write("# Frame, Zy (Degrees)\n")
    for i in range(len(Zy_all)):
        file.write("{}\t{}\n".format(i,Zy_all[i]))

with open("Zx_Zy.csv","a") as file:  
    file.write("# Zx (Degrees), Zy (Degrees)\n")
    for i in range(len(Zx_all)):
        file.write("{}\t{}\n".format(Zx_all[i],Zy_all[i]))

with open("P.csv","a") as file:
    file.write("# Frame, P (Degrees)\n")
    for i in range(len(Ptheta_all)):
        file.write("{}\t{}\n".format(i,Ptheta_all[i]))

with open("Amp.csv","a") as file:
    file.write("# Frame, Amp (Degrees)\n")
    for i in range(len(Ar_all)):
        file.write("{}\t{}\n".format(i,Ar_all[i]))

with open("P_Amp.csv","a") as file:
    file.write("# P (Degrees), Amp (Degrees)\n")
    for i in range(len(Ar_all)):
        file.write("{}\t{}\n".format(Ptheta_all[i],Ar_all[i]))

# Mv files to working directory
mvfiles = "mv Zx.csv Zy.csv Zx_Zy.csv P.csv Amp.csv P_Amp.csv {}".format(current_dir)
os.system(mvfiles)
