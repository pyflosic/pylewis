#   Copyright 2019 PyLEWIS developers 
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
# Authors: 
#           Sebastian Schwalbe (SS) 
#           Kai Trepte (KT) 
#
#
# pylewis  
# 
# Lewis structures as input for FLO-SIC 
# - especially for automatic guess creators (fodmc,density-centroids,pyeff etc.) 
#
# Author:       Sebastian Schwalbe 
#		Kai Trepte (Geometrical aspects)
# Date:         20.09.2018 
#		- Inital script creation 
#	        21.09.2018 
#		- Update to work for 3D structures (e.g. C60) 
#		Introduction of two new input parameters (aromatic,dim,addH)
#		- Added the functionality to start from pdb files
#		  xyz files can be converted with Avogadro to pdb files which have the right bond connect information
#		  ase and babel only convert nuclei information and no bonding information 
#		25.09.2018 
#		- Added support for atoms which are one electron radicals (e.g CH3, allyl radical) 
#		- Added support for core electrons from 1-4 electrons (1s, 2s) 
#               25.04.2019 
#               - changed name to pylewis 
# Test:         Institute Laptop 
# Links:        RDkit homepage https://www.rdkit.org/docs/GettingStartedInPython.html#writing-molecules
# Installation: RDkit, under linux sudo apt-get install python-rdkit librdkit1 rdkit-data
# Input:        Similes strings are available in every chemistry database 
#               e.g. @pubchem https://pubchem.ncbi.nlm.nih.gov/compound/anthracene
# TODO:		- SDF files include also information about radicals "RAD" flag 
#		  http://c4.cabrillo.edu/404/ctfile.pdf
#		- for the core electrons of each atom in the system one need to count the bondings 
#		  and remove the electrons to get the number of core electrons 
#		  as well one needs to count for the charge state of this atom and 
#		  add or remove electrons
#		- analyse the charge state may work with formal charges
#		- lone pairs detection 
 
from rdkit import Chem
from rdkit.Chem import AllChem,Draw
from ase.io import read,write
import numpy as np
from ase.atoms import Atom, Atoms

def lewis2xyz(lewis,f='SMILES',aromatic=False,dim=2,addH=True):
        # Input 
        #       lewis 		... 	smiles format string 
	#	aromatic 	... 	is the system aromatic 
	#	dim		...	2d == 2 and 3d == 3 
        # Output 
        #       Bonding information 
        #       Lewis pictures
        #       sdf and xyz files 

	if f == 'pdb':
		m = Chem.MolFromPDBFile(lewis)
	if f == 'xyz':
		struct = read(lewis)
		write('lewis.pdb',struct,'proteindatabank')
		m = Chem.MolFromPDBFile(lewis)
       	if f == 'SMILES':
		m = Chem.MolFromSmiles(lewis)
	
	# Change armatic bonds to alternate single and double bonds 
	if aromatic == True:	
		if m.GetBondWithIdx(0).GetIsAromatic() == True:
       			Chem.Kekulize(m)

	if addH == True:	
		m=Chem.AddHs(m) 
	if dim == 1 or dim == 2: 
                AllChem.Compute2DCoords(m)
	if dim == 3:
                AllChem.EmbedMolecule(m)
		AllChem.UFFOptimizeMolecule(m)

	fods = Atoms()
       	# core electrons  
	for atom in m.GetAtoms():
		idx = atom.GetIdx()
		# atomic number 
		N = m.GetAtomWithIdx(idx).GetAtomicNum()
		# number of valence electrons 
		V = m.GetAtomWithIdx(idx).GetTotalValence()
		# number of radicals 
		R = atom.GetNumRadicalElectrons()
		# formal charge of the atom 
		F = atom.GetFormalCharge()
		# number of core electrons 
		C = N-V-R-F
		# lone pairs  
		#L = 0
		#bonds = atom.GetBonds()
		#for b in range(len(bonds)):
		#	bond_type = bonds[b].GetBondType()
		#	if str(bond_type) == 'SINGLE':
		#		L = L + 1 		
		#	if str(bond_type) == 'DOUBLE':
		#		L = L + 2 
		#	if str(bond_type) == 'TRIPLE':
                #                L = L + 3
		#print(V-L)
		if C > 0:
			pos = m.GetConformer().GetAtomPosition(idx)
			offset = 0.002
			if C == 1:
				fods.append(Atom('X',[pos[0],pos[1],pos[2]]))			
			if C == 2:
				fods.append(Atom('X',[pos[0],pos[1],pos[2]]))
				fods.append(Atom('He',[pos[0]+offset,pos[1]+offset,pos[2]+offset]))
			if C == 3:
				fods.append(Atom('X',[pos[0],pos[1],pos[2]]))
                                fods.append(Atom('He',[pos[0]+offset,pos[1]+offset,pos[2]+offset])) 
				fods.append(Atom('X',[pos[0]+offset,pos[1]+offset,pos[2]+4*offset]))
			if C == 4:
                                fods.append(Atom('X',[pos[0],pos[1],pos[2]]))
                                fods.append(Atom('He',[pos[0]+offset,pos[1]+offset,pos[2]+offset]))
                                fods.append(Atom('X',[pos[0]+offset,pos[1]+offset,pos[2]+40*offset]))
				fods.append(Atom('He',[pos[0]+offset,pos[1]+offset,pos[2]-40*offset]))
			if C == 5:
                                fods.append(Atom('X',[pos[0],pos[1],pos[2]]))
                                fods.append(Atom('He',[pos[0]+offset,pos[1]+offset,pos[2]+offset]))
                                fods.append(Atom('X',[pos[0]+offset,pos[1]+offset,pos[2]+40*offset]))
                                fods.append(Atom('He',[pos[0]+offset,pos[1]+offset,pos[2]-40*offset]))
				fods.append(Atom('X',[pos[0]+40*offset,pos[1]+40*offset,pos[2]+offset]))		
			if C == 6: 
                                fods.append(Atom('X',[pos[0],pos[1],pos[2]]))
                                fods.append(Atom('He',[pos[0]+offset,pos[1]+offset,pos[2]+offset]))
                                
                                fods.append(Atom('X',[pos[0]+40*offset,pos[1]+40*offset,pos[2]+40*offset]))
                                fods.append(Atom('He',[pos[0]+offset,pos[1]-40*offset,pos[2]-40*offset]))
                                fods.append(Atom('X',[pos[0]-40*offset,pos[1]+40*offset,pos[2]-40*offset]))
				fods.append(Atom('He',[pos[0]-40*offset,pos[1]-40*offset,pos[2]+40*offset]))

			if C == 7:
                                fods.append(Atom('X',[pos[0],pos[1],pos[2]]))
                                fods.append(Atom('He',[pos[0]+offset,pos[1]+offset,pos[2]+offset]))

                                fods.append(Atom('X',[pos[0]+40*offset,pos[1]+40*offset,pos[2]+40*offset]))
                                fods.append(Atom('He',[pos[0]+offset,pos[1]-40*offset,pos[2]-40*offset]))
                                fods.append(Atom('X',[pos[0]-40*offset,pos[1]+40*offset,pos[2]-40*offset]))
                                fods.append(Atom('He',[pos[0]-40*offset,pos[1]-40*offset,pos[2]+40*offset]))

				fods.append(Atom('He',[pos[0]+42*offset,pos[1]+42*offset,pos[2]+42*offset]))
			
			if C == 8:
                                fods.append(Atom('X',[pos[0],pos[1],pos[2]]))
                                fods.append(Atom('He',[pos[0]+offset,pos[1]+offset,pos[2]+offset]))

                                fods.append(Atom('X',[pos[0]+40*offset,pos[1]+40*offset,pos[2]+40*offset]))
                                fods.append(Atom('He',[pos[0]+40*offset,pos[1]-40*offset,pos[2]-40*offset]))
                                fods.append(Atom('X',[pos[0]-40*offset,pos[1]+40*offset,pos[2]-40*offset]))
                                fods.append(Atom('He',[pos[0]-40*offset,pos[1]-40*offset,pos[2]+40*offset]))

                                fods.append(Atom('He',[pos[0]+42*offset,pos[1]+42*offset,pos[2]+42*offset]))
				fods.append(Atom('X',[pos[0]+42*offset,pos[1]-42*offset,pos[2]-42*offset]))


	# find radicals 
	for atom in m.GetAtoms():
		rad = atom.GetNumRadicalElectrons()
		#print(rad)
		if rad == 1:
			idx = atom.GetIdx()
			pos = m.GetConformer().GetAtomPosition(idx)
			bonds = atom.GetBonds()
			vec_BA_x = 0
			vec_BA_y = 0 
			vec_BA_z = 0
			BA = []
			for b in range(len(bonds)):
                		ai = bonds[b].GetBeginAtomIdx()
                		aj = bonds[b].GetEndAtomIdx()
                		bond_type = bonds[b].GetBondType()

                		# Positions 
                		pos_ai = m.GetConformer().GetAtomPosition(ai)
                		pos_aj = m.GetConformer().GetAtomPosition(aj)
		
				# adding bonding vectors for 3d 	
				vec_BA_x = vec_BA_x + (pos_aj[0]-pos_ai[0])
				vec_BA_y = vec_BA_y + (pos_aj[1]-pos_ai[1])
				vec_BA_z = vec_BA_z + (pos_aj[2]-pos_ai[2]) 
				# save bonding vectors for 2d 
				BA.append(np.array([(pos_aj[0]-pos_ai[0]),(pos_aj[1]-pos_ai[1]),(pos_aj[2]-pos_ai[2])]))
			
			if dim == 1:
				fods.append(Atom('X',[-1*BA[0][0],-1*BA[0][1],-1*BA[0][2]]))	
			if dim == 2: 
				vec = np.cross(BA[0],BA[1])
				norm = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
				pos = vec/norm 
				fods.append(Atom('X',[pos[0],pos[1],pos[2]]))

			if dim == 3:
				norm = np.sqrt(vec_BA_x**2 + vec_BA_y**2 + vec_BA_z**2)	
				pos = np.array([-1*vec_BA_x,-1*vec_BA_y,-1*vec_BA_z])/norm
				fods.append(Atom('X',[pos[0],pos[1],pos[2]])) 
			


        bonds = m.GetBonds()
        for b in range(len(bonds)):
                ai = bonds[b].GetBeginAtomIdx()
                aj = bonds[b].GetEndAtomIdx()
                bond_type = bonds[b].GetBondType()

                # Positions 
                pos_ai = m.GetConformer().GetAtomPosition(ai)
		pos_aj = m.GetConformer().GetAtomPosition(aj)
                # Bondlength 
                #print('%0.5f' % np.linalg.norm([pos_ai[0]-pos_aj[0],pos_ai[1]-pos_aj[1],pos_ai[2]-pos_aj[2]]))
                # Middle of the vector 
                MP = [(pos_ai[0]+pos_aj[0])/2.,(pos_ai[1]+pos_aj[1])/2.,(pos_ai[2]+pos_aj[2])/2.]
                #print('%0.5f %0.5f %0.5f' %(MP[0],MP[1],MP[2]))
                # Print positions 
                #print('%0.5f %0.5f %0.5f' %(pos_ai[0],pos_ai[1],pos_ai[2]))
                #print('%0.5f %0.5f %0.5f' %(pos_aj[0],pos_aj[1],pos_aj[2]))

                for spin in ['He','X']:
                        if spin == 'He':
                                offset = 0.000
                        if spin == 'X':
                                #offset = 0.002
				offset = 0.02 
                        if str(bond_type) == 'SINGLE':
                                #print('%s %0.5f %0.5f %0.5f' %(spin,MP[0]+offset,MP[1]+offset,MP[2]+offset))
                                fods.append(Atom(spin,[MP[0]+offset,MP[1]+offset,MP[2]+offset]))

                        if str(bond_type) == 'DOUBLE':
                                #print('%s %0.5f %0.5f %0.5f' %(spin,MP[0]+offset,MP[1]+offset,MP[2]+0.5+offset))
                                #print('%s %0.5f %0.5f %0.5f' %(spin,MP[0]+offset,MP[1]+offset,MP[2]-0.5+offset))
                                fods.append(Atom(spin,[MP[0]+offset,MP[1]+offset,MP[2]+0.5+offset]))
                                fods.append(Atom(spin,[MP[0]+offset,MP[1]+offset,MP[2]-0.5+offset]))

                        if str(bond_type) == 'TRIPLE':
                                #print('%s %0.5f %0.5f %0.5f' %(spin,MP[0]+offset,MP[1]+offset,MP[2]+0.5+offset))
                                #print('%s %0.5f %0.5f %0.5f' %(spin,MP[0]+offset,MP[1]+offset,MP[2]-0.5+offset))
                                #print('%s %0.5f %0.5f %0.5f' %(spin,MP[0]-0.5+offset,MP[1]+offset,MP[2]+offset))
                                fods.append(Atom(spin,[MP[0]+offset,MP[1]+offset,MP[2]+0.5+offset]))
                                fods.append(Atom(spin,[MP[0]+offset,MP[1]+offset,MP[2]-0.5+offset]))
                                fods.append(Atom(spin,[MP[0]-0.5+offset,MP[1]+offset,MP[2]+offset]))

                w = Chem.SDWriter('%s.sdf' % 'pylewis')
                w.write(m)
                w.flush()

                # sdf2xyz       
                struct = read('%s.sdf' % 'pylewis')
                struct_fods = Atoms()
                struct_fods.extend(struct)
                struct_fods.extend(fods)

                #print(struct.get_chemical_symbols())
                write('%s.xyz' % 'pylewis',struct,'xyz')
                write('%s_fods.xyz' % 'pylewis',struct_fods,'xyz')

# Test systems 
# N2, some strange test, antracene, benzene, coronene 
#C60 = 'C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4=C3C3=C5C(=C81)C%10=C23'
#lewis = ['N#N','C1=C-C1','C1=CC=C2C=CC=CC2=C1','C1=C-C=C-C=C1','C1=CC2=C3C4=C1C=CC5=C4C6=C(C=C5)C=CC7=C6C3=C(C=C2)C=C7',C60,'C'][0]

# coronene 
#lewis2xyz('C1=CC2=C3C4=C1C=CC5=C4C6=C(C=C5)C=CC7=C6C3=C(C=C2)C=C7',aromatic=True,dim=2)
# N2
#lewis2xyz('N#N',aromatic=False,dim=2)
# C60
#C60 = 'C12=C3C4=C5C6=C1C7=C8C9=C1C%10=C%11C(=C29)C3=C2C3=C4C4=C5C5=C9C6=C7C6=C7C8=C1C1=C8C%10=C%10C%11=C2C2=C3C3=C4C4=C5C5=C%11C%12=C(C6=C95)C7=C1C1=C%12C5=C%11C4=C3C3=C5C(=C81)C%10=C23'
#lewis2xyz(C60,aromatic=True,dim=3)
#lewis2xyz('C60.pdb',f='pdb',aromatic=True,dim=3)
# CH4 
#lewis2xyz('C',aromatic=False,dim=3)
#lewis2xyz('CH4.pdb',f='pdb',aromatic=False,dim=3)
# CH3 radical 
#lewis2xyz('[CH3]',aromatic=False,dim=2)
# allyl radical 
#lewis2xyz('C=C[CH2]',aromatic=False,dim=2)
# C6H6
#lewis2xyz('C1=C-C=C-C=C1',aromatic=True,dim=2)
# CN radical 
#lewis2xyz('[C]#N',aromatic=False,dim=1)
# CN anion 
#lewis2xyz('[C-]#N',aromatic=False,dim=1)
# HCN
#lewis2xyz('C#N',aromatic=False,dim=1)
# NH3 
#lewis2xyz('N',aromatic=False,dim=3)
# BH3 
#lewis2xyz('BH3.pdb',f='pdb',aromatic=False,dim=3,addH=False)
#lewis2xyz('B',aromatic=False,dim=3)
# N3H 
#lewis2xyz('N=[N+]=[N-]',aromatic=False,dim=3)
# H2O 
#lewis2xyz('O',aromatic=False,dim=3)
# HF 
#lewis2xyz('F',aromatic=False,dim=3)
# F2 
#lewis2xyz('FF',aromatic=False,dim=3)
# F3N 
#lewis2xyz('FN(F)F',aromatic=False,dim=3)
# H2O2 
#lewis2xyz('OO',aromatic=False,dim=3)
# LiF
# ATTENTION: Doesnt work for now 
#lewis2xyz('LiF.pdb',f='pdb',aromatic=False,dim=3,addH=False)
# C2H2 
#lewis2xyz('C#C',aromatic=False,dim=3)
# C4H4 
#lewis2xyz('C=C',aromatic=False,dim=2)
