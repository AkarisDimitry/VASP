# *** warning supresion
import warnings, os, time
import urllib.request
warnings.filterwarnings("ignore")

# *** warning supresion
import warnings; warnings.filterwarnings("ignore")

# *** numpy libraries
import numpy as np 	
import matplotlib.pyplot as plt
try:	
	from ase.visualize import view
	from ase import Atoms
	from ase.visualize.plot import plot_atoms
except: print('can not load ASE module. (pipX install ase-atomistics)')

# *** 
from pypdb import * # pip install pypdb
#urllib.request.urlretrieve('http://files.rcsb.org/download/101M.pdb', '101m.pdb')

class PDB(object): # generador de datos
	def __init__(self, path=None, name=None):
		self.path = path 
		self.name = name 

		# PDB metadata
		self.HEADER = None
		self.TITLE  = None
		self.COMPND = None
		self.SOURCE = None
		self.MUTATION = None
		self.ENGINEERED = None

		self.HETNAM = None
		self.FORMUL = None

		# ==== PROTEIN ATOMS ==== # # ==== PROTEIN ATOMS ==== # # ==== PROTEIN ATOMS ==== # # ==== PROTEIN ATOMS ==== #
		self.ATOM   = None
									# 1 -  6        Record name   "ATOM  "
		self.ATOM_num = None 		# 7 - 11        Integer       serial       Atom  serial number.
		self.ATOM_idx = None 		# 13 - 16        Atom          name         Atom name.
									# 17             Character     altLoc       Alternate location indicator.
		self.ATOM_alt = None		# 18 - 20        Residue name  resName      Residue name.
		self.ATOM_ani = None		# 22             Character     chainID      Chain identifier.

		self.ATOM_res = None		# 23 - 26        Integer       resSeq       Residue sequence number.
		self.ATOM_iCode = None 		# 27             AChar         iCode        Code for insertion of residues.
		self.ATOM_xyz = None		# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
									# 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
									# 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		self.ATOM_occupancy = None	# 55 - 60        Real(6.2)     occupancy    Occupancy.

		self.ATOM_tempFactor = None	# 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		self.ATOM_element = None	# 77 - 78        LString(2)    element      Element symbol, right-justified.
		self.ATOM_charge = None		# 79 - 80        LString(2)    charge       Charge  on the atom.

		# ==== HETERO ATOMS ==== # # ==== HETERO ATOMS ==== # # ==== HETERO ATOMS ==== # # ==== HETERO ATOMS ==== #
		self.HETATM   = None
									# 1 -  6        Record name   "ATOM  "
		self.HETATM_num = None 		# 7 - 11        Integer       serial       Atom  serial number.
		self.HETATM_idx = None 		# 13 - 16        Atom          name         Atom name.
									# 17             Character     altLoc       Alternate location indicator.
		self.HETATM_alt = None		# 18 - 20        Residue name  resName      Residue name.
		self.HETATM_ani = None		# 22             Character     chainID      Chain identifier.

		self.HETATM_res = None		# 23 - 26        Integer       resSeq       Residue sequence number.
		self.HETATM_iCode = None 		# 27             AChar         iCode        Code for insertion of residues.
		self.HETATM_xyz = None		# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
									# 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
									# 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		self.HETATM_occupancy = None	# 55 - 60        Real(6.2)     occupancy    Occupancy.

		self.HETATM_tempFactor = None	# 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		self.HETATM_element = None	# 77 - 78        LString(2)    element      Element symbol, right-justified.
		self.HETATM_charge = None		# 79 - 80        LString(2)    charge       Charge  on the atom.


		# ==== HETERO ATOMS ==== # # ==== HETERO ATOMS ==== # # ==== HETERO ATOMS ==== # # ==== HETERO ATOMS ==== #
		self.ATMall   = None
										# 1 -  6        Record name   "ATOM  "
		self.ATMall_num = None 			# 7 - 11        Integer       serial       Atom  serial number.
		self.ATMall_idx = None 			# 13 - 16        Atom          name         Atom name.
										# 17             Character     altLoc       Alternate location indicator.
		self.ATMall_alt = None			# 18 - 20        Residue name  resName      Residue name.
		self.ATMall_ani = None			# 22             Character     chainID      Chain identifier.

		self.ATMall_res 	= None		# 23 - 26        Integer       resSeq       Residue sequence number.
		self.ATMall_iCode 	= None 		# 27             AChar         iCode        Code for insertion of residues.
		self.ATMall_xyz 	= None		# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
										# 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
										# 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		self.ATMall_occupancy = None	# 55 - 60        Real(6.2)     occupancy    Occupancy.

		self.ATMall_tempFactor 	= None	# 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		self.ATMallM_element 	= None	# 77 - 78        LString(2)    element      Element symbol, right-justified.
		self.ATMall_charge   	= None	# 79 - 80        LString(2)    charge       Charge  on the atom.


	def timer(func):
		def wrapper(*args, **kwargs):
			before = time.time()
			r = func(*args, **kwargs)
			name = str(func).split(' ')[1]
			print( f'Resuldo of {name} : {r} (execution time {time.time()-before}s)' ) 
			return r
		return wrapper

	@timer
	def READ(self, path=None, name=None, save=True, v=True):
		path = path if not path is None else self.path 
		name = name if not name is None else self.name

		file = open(f'{path}/{name}')
		HEADER = []
		TITLE  = []
		COMPND = []
		SOURCE = []
		ENGINEERED = False
		MUTATION = False
		HETNAM = [] 
		FORMUL = []

		# ==== PROTEIN ATOMS ==== # # ==== PROTEIN ATOMS ==== # # ==== PROTEIN ATOMS ==== # # ==== PROTEIN ATOMS ==== #
		ATOM   = []
								# 1 -  6        Record name   "ATOM  "
		ATOM_num = [] 			# 7 - 11        Integer       serial       Atom  serial number.
		ATOM_idx = [] 	 		# 13 - 16        Atom          name         Atom name.
								# 17             Character     altLoc       Alternate location indicator.
		ATOM_alt = [] 			# 18 - 20        Residue name  resName      Residue name.
		ATOM_ani = [] 			# 22             Character     chainID      Chain identifier.

		ATOM_res   = [] 		# 23 - 26        Integer       resSeq       Residue sequence number.
		ATOM_iCode = [] 		# 27             AChar         iCode        Code for insertion of residues.
		ATOM_xyz   = [] 		# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
								# 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
								# 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		ATOM_occupancy = [] 	# 55 - 60        Real(6.2)     occupancy    Occupancy.

		ATOM_tempFactor = []	# 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		ATOM_element    = [] 	# 77 - 78        LString(2)    element      Element symbol, right-justified.
		ATOM_charge     = []	# 79 - 80        LString(2)    charge       Charge  on the atom.
		
		# ==== HETERO ATOMS ==== # # ==== HETERO ATOMS ==== # # ==== HETERO ATOMS ==== # # ==== HETERO ATOMS ==== #
		HETATM   = []
								# 1 -  6        Record name   "ATOM  "
		HETATM_num = [] 		# 7 - 11        Integer       serial       Atom  serial number.
		HETATM_idx = [] 		# 13 - 16        Atom          name         Atom name.
								# 17             Character     altLoc       Alternate location indicator.
		HETATM_alt = []			# 18 - 20        Residue name  resName      Residue name.
		HETATM_ani = []			# 22             Character     chainID      Chain identifier.

		HETATM_res = []			# 23 - 26        Integer       resSeq       Residue sequence number.
		HETATM_iCode = [] 		# 27             AChar         iCode        Code for insertion of residues.
		HETATM_xyz = []			# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
								# 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
								# 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		HETATM_occupancy = []	# 55 - 60        Real(6.2)     occupancy    Occupancy.

		HETATM_tempFactor = []	# 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		HETATM_element    = []	# 77 - 78        LString(2)    element      Element symbol, right-justified.
		HETATM_charge     = []	# 79 - 80        LString(2)    charge       Charge  on the atom.

		for line in file:
			vec = [m for m in line.replace('\t',' ').split(' ') if m != '' and m != '\n']

			if vec[0] == 'HEADER': 	HEADER.append(''.join(vec[1:]))
			if vec[0] == 'TITLE': 	TITLE.append(''.join(vec[1:]))
			if vec[0] == 'COMPND': 	COMPND.append(''.join(vec[1:]))
			if vec[0] == 'SOURCE': 	SOURCE.append(''.join(vec[1:]))
			if vec[0] == 'HETNAM': 	HETNAM.append( [line[11:14].strip(), line[15:].strip()] )
			if vec[0] == 'FORMUL': 	FORMUL.append( [line[7:11].strip(), line[12:17].strip(), line[18:].strip() ] )

			if vec[0] == 'COMPND': 
				if vec[2] == 'ENGINEERED': 	ENGINEERED = True
				if vec[2] == 'MUTATION': 	MUTATION = True

			if vec[0] == 'ATOM': 	
				ATOM.append(vec[1:])
															# 1 -  6        Record name   "ATOM  "
				ATOM_num.append( int(line[7:11])) 			# 7 - 11        Integer       serial       Atom  serial number.
				ATOM_idx.append( line[13:16]) 				# 13 - 16        Atom          name         Atom name.
															# 17             Character     altLoc       Alternate location indicator.
				ATOM_alt.append( line[18:20])				# 18 - 20        Residue name  resName      Residue name.
				ATOM_ani.append( line[22])					# 22             Character     chainID      Chain identifier.

				ATOM_res.append( int(line[23:26]) )			# 23 - 26        Integer       resSeq       Residue sequence number.
				ATOM_iCode.append( line[27]) 				# 27             AChar         iCode        Code for insertion of residues.
				ATOM_xyz.append([float(line[31:38]), 		# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
								 float(line[39:46]),		# 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
								 float(line[47:54])])		# 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
				ATOM_occupancy.append( float(line[55:60]))	# 55 - 60        Real(6.2)     occupancy    Occupancy.
				ATOM_tempFactor.append( float(line[61:66]))	# 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
				ATOM_element.append( line[77:78])			# 77 - 78        LString(2)    element      Element symbol, right-justified.
				ATOM_charge.append( line[79:80])			# 79 - 80        LString(2)    charge       Charge  on the atom.

			if vec[0] == 'HETATM': 	
				HETATM.append(vec[1:])
																# 1 -  6        Record name   "ATOM  "
				HETATM_num.append( int(line[7:11])) 			# 7 - 11        Integer       serial       Atom  serial number.
				HETATM_idx.append( line[13:16]) 				# 13 - 16        Atom          name         Atom name.
																# 17             Character     altLoc       Alternate location indicator.
				HETATM_alt.append( line[18:20])					# 18 - 20        Residue name  resName      Residue name.
				HETATM_ani.append( line[22])					# 22             Character     chainID      Chain identifier.

				HETATM_res.append( int(line[23:26]) )			# 23 - 26        Integer       resSeq       Residue sequence number.
				HETATM_iCode.append( line[27]) 					# 27             AChar         iCode        Code for insertion of residues.
				HETATM_xyz.append([float(line[31:38]), 			# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
								 float(line[39:46]),			# 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
								 float(line[47:54])])			# 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
				HETATM_occupancy.append( float(line[55:60]))	# 55 - 60        Real(6.2)     occupancy    Occupancy.
				HETATM_tempFactor.append( float(line[61:66]))	# 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
				HETATM_element.append( 							# 77 - 78        LString(2)    element      Element symbol, right-justified.
										line[75:78].strip()[0]+line[75:78].strip()[1:].lower() )	
				HETATM_charge.append( line[79:80])				# 79 - 80        LString(2)    charge       Charge  on the atom.

		if save:
			self.HEADER = HEADER
			self.TITLE = TITLE
			self.COMPND = COMPND
			self.SOURCE = SOURCE

			self.ENGINEERED = ENGINEERED
			self.MUTATION = MUTATION
			
			self.HETNAM = HETNAM
			self.FORMUL = FORMUL

			# ==== PROTEIN ATOMS ==== # # ==== PROTEIN ATOMS ==== # # ==== PROTEIN ATOMS ==== # # ==== PROTEIN ATOMS ==== #
			self.ATOM = ATOM
													# 1 -  6        Record name   "ATOM  "
			self.ATOM_num = ATOM_num				# 7 - 11        Integer       serial       Atom  serial number.
			self.ATOM_idx = ATOM_idx 				# 13 - 16        Atom          name         Atom name.
													# 17             Character     altLoc       Alternate location indicator.
			self.ATOM_alt = ATOM_alt				# 18 - 20        Residue name  resName      Residue name.
			self.ATOM_ani = ATOM_ani				# 22             Character     chainID      Chain identifier.

			self.ATOM_res = ATOM_res				# 23 - 26        Integer       resSeq       Residue sequence number.
			self.ATOM_iCode = ATOM_iCode 			# 27             AChar         iCode        Code for insertion of residues.
			self.ATOM_xyz = ATOM_xyz				# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
													# 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
													# 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
			self.ATOM_occupancy = ATOM_occupancy	# 55 - 60        Real(6.2)     occupancy    Occupancy.

			self.ATOM_tempFactor = ATOM_tempFactor	# 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
			self.ATOM_element    = ATOM_element		# 77 - 78        LString(2)    element      Element symbol, right-justified.
			self.ATOM_charge     = ATOM_charge		# 79 - 80        LString(2)    charge       Charge  on the atom.

			# ==== HETERO ATOMS ==== # # ==== HETERO ATOMS ==== # # ==== HETERO ATOMS ==== # # ==== HETERO ATOMS ==== #
			self.HETATM = HETATM
														# 1 -  6        Record name   "ATOM  "
			self.HETATM_num = HETATM_num				# 7 - 11        Integer       serial       Atom  serial number.
			self.HETATM_idx = HETATM_idx 				# 13 - 16        Atom          name         Atom name.
														# 17             Character     altLoc       Alternate location indicator.
			self.HETATM_alt = HETATM_alt				# 18 - 20        Residue name  resName      Residue name.
			self.HETATM_ani = HETATM_ani				# 22             Character     chainID      Chain identifier.

			self.HETATM_res   = HETATM_res				# 23 - 26        Integer       resSeq       Residue sequence number.
			self.HETATM_iCode = HETATM_iCode 			# 27             AChar         iCode        Code for insertion of residues.
			self.HETATM_xyz   = HETATM_xyz				# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
														# 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
														# 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
			self.HETATM_occupancy = HETATM_occupancy	# 55 - 60        Real(6.2)     occupancy    Occupancy.

			self.HETATM_tempFactor = HETATM_tempFactor	# 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
			self.HETATM_element    = HETATM_element		# 77 - 78        LString(2)    element      Element symbol, right-justified.
			self.HETATM_charge     = HETATM_charge		# 79 - 80        LString(2)    charge       Charge  on the atom.

	def ATOMS_join(self, ):
		self.ATMall = self.ATOM + self.HETATM
																				# 1 -  6        Record name   "ATOM  "
		self.ATMall_num = self.HETATM_num + self.ATOM_num						# 7 - 11        Integer       serial       Atom  serial number.
		self.ATMall_idx = self.HETATM_idx + self.ATOM_idx						# 13 - 16        Atom          name         Atom name.
																				# 17             Character     altLoc       Alternate location indicator.
		self.ATMall_alt = self.HETATM_alt + self.ATOM_alt						# 18 - 20        Residue name  resName      Residue name.
		self.ATMall_ani = self.HETATM_ani + self.ATOM_ani						# 22             Character     chainID      Chain identifier.

		self.ATMall_res   = self.HETATM_res + self.ATOM_res						# 23 - 26        Integer       resSeq       Residue sequence number.
		self.ATMall_iCode = self.HETATM_iCode + self.ATOM_iCo 					# 27             AChar         iCode        Code for insertion of residues.
		self.ATMall_xyz   = self.HETATM_xyz + self.ATOM_xyz						# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
																				# 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
																				# 47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		self.ATMall_occupancy = self.HETATM_occupancy + self.ATOM_occupancy		# 55 - 60        Real(6.2)     occupancy    Occupancy.

		self.ATMall_tempFactor = self.HETATM_tempFactor + self.ATOM_tempFactor	# 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		self.ATMall_element    = self.HETATM_element + self.ATOM_element		# 77 - 78        LString(2)    element      Element symbol, right-justified.
		self.ATMall_charge     = self.HETATM_charge + self.ATOM_charge			# 79 - 80        LString(2)    charge       Charge  on the atom.

	def plot_img(self, save=False, path='./no_name.png'):
		if type(self.ATOM_xyz) != type(None):
			fig, ax = plt.subplots(2,2)
			plt.subplots_adjust(wspace=0, hspace=0)
			atoms = self.get_ase()
			rot=['0x,0y,0z', '-90x,0y,0z', '0x,90y,0z', '0x,0y,90z']
			for itr, a in enumerate(ax.flatten()):
			    plot_atoms(atoms, a, rotation=rot[itr], show_unit_cell=True, scale=0.5)
			    a.axis('off')
			    #a.set_facecolor('black')
			    a.autoscale()
			fig.set_facecolor('white')
			if save:
				fig.savefig(f'{path}' , dpi=100, pad_inches=0.1, bbox_inches='tight', horizontalalignment='right') 

	def plot(self, save=False, path='./no_name.png'):
		if type(self.ATOM_xyz) != type(None):
			atoms = self.get_ase()
			view(atoms)#, viewer='VMD')

			if save:
				fig.savefig(f'{path}' , dpi=100, pad_inches=0.1, bbox_inches='tight', horizontalalignment='right') 

	def get_ase(self, postions=None, element=None, pbc=[0,0,0]):
		element  = element  if not element  is None else self.ATOM_element
		postions = postions if not postions is None else self.ATOM_xyz

		return Atoms(	element,
             positions=	postions,
             #cell     =	self.cell,
             pbc      =	[0, 0, 0])

	def is_mutated(self, ): return self.MUTATION
	def is_engenierated(self, ): return self.ENGINEERED

	def has_heteroatom(self, name=None, match='aprox'):
		if match == 'exact':
			return any([ any([name in atom for atom in ligand]) for ligand in PDB1.FORMUL])
		elif match == 'aprox':
			return any([ any([name.lower() in atom.lower() for atom in ligand]) for ligand in PDB1.FORMUL])

PDB1 = PDB(path='.', name='6zys.pdb')
PDB1.READ()

print( PDB1.has_heteroatom('HOH') )



