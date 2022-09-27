# *** warning supresion
import warnings, os, time, argparse
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

try:	from src import Logs
except:	
	try: import Logs as Logs
	except: print('WARNING :: Set.import_libraries() :: can not import ORR ')

# *** 
#from pypdb import * # pip install pypdb
#urllib.request.urlretrieve('http://files.rcsb.org/download/101M.pdb', '101m.pdb')

class PDB(object): # generador de datos
	def __init__(self, path=None, name=None):
		self.path = path 
		self.name = name 

		# PDB metadata
		self.HEADER = []			# HEADER contains the following fields: classification="PROTEIN" (or imported value), date, idCode="NONE" (or imported value).
		self.TITLE  = []			# TITLE, SOURCE, KEYWDS, EXPDTA: The imported value is exported. Default: "NULL".
		self.COMPND = []	 		# COMPND: The imported value is exported. Default: "MOLECULE:name", where "name" is the molecule name.
		self.SOURCE = []			# COMPND: The imported value is exported. Default: "MOLECULE:name", where "name" is the molecule name.
		self.AUTHOR = []			# AUTHOR: The imported value is exported. Default: "Marvin".
		self.REMARK = [] 			# REMARK In entries where REMARK 0 is included as described above, REMARK 900 will also reflect the reuse of existing experimental data. 
									# https://www.wwpdb.org/documentation/file-format-content/format33/remarks1.html
									# http://www.bmsc.washington.edu/CrystaLinks/man/pdb/part_31.html
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

		self.connectivity = None

	def timer(func):
		def wrapper(*args, **kwargs):
			before = time.time()
			r = func(*args, **kwargs)
			name = str(func).split(' ')[1]
			print( f'Resuldo of {name} : {r} (execution time {time.time()-before}s)' ) 
			return r
		return wrapper

	@timer
	@Logs.LogDecorator()
	def READ(self, path=None, name=None, save=True, v=True):
		# https://www.wwpdb.org/documentation/file-format
		# https://files.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
		path = path if not path is None else self.path 
		name = name if not name is None else self.name

		file = open(f'{path}/{name}')
		HEADER = []			# HEADER contains the following fields: classification="PROTEIN" (or imported value), date, idCode="NONE" (or imported value).
		TITLE  = []			# TITLE, SOURCE, KEYWDS, EXPDTA: The imported value is exported. Default: "NULL".
		COMPND = [] 		# COMPND: The imported value is exported. Default: "MOLECULE:name", where "name" is the molecule name.
		SOURCE = []			# COMPND: The imported value is exported. Default: "MOLECULE:name", where "name" is the molecule name.
		AUTHOR = []			# AUTHOR: The imported value is exported. Default: "Marvin".
		REMARK = [] 		# REMARK In entries where REMARK 0 is included as described above, REMARK 900 will also reflect the reuse of existing experimental data. 
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

			if vec[0] == 'HEADER': 	HEADER.append(''.join(vec[1:]))													# HEADER contains the following fields: classification="PROTEIN" (or imported value), date, idCode="NONE" (or imported value).
			if vec[0] == 'TITLE': 	TITLE.append(''.join(vec[1:]))													# TITLE, SOURCE, KEYWDS, EXPDTA: The imported value is exported. Default: "NULL".
			if vec[0] == 'COMPND': 	COMPND.append(''.join(vec[1:])) 												# COMPND: The imported value is exported. Default: "MOLECULE:name", where "name" is the molecule name.
			if vec[0] == 'SOURCE': 	SOURCE.append(''.join(vec[1:]))													# COMPND: The imported value is exported. Default: "MOLECULE:name", where "name" is the molecule name.
			if vec[0] == 'AUTHOR': 	AUTHOR.append(''.join(vec[1:]))													# AUTHOR: The imported value is exported. Default: "Marvin".
			if vec[0] == 'REMARK':  REMARK.append(' '.join(vec[1:])) 												# # REMARK In entries where REMARK 0 is included as described above, REMARK 900 will also reflect the reuse of existing experimental data. 
			if vec[0] == 'HETNAM': 	HETNAM.append( [line[11:14].strip(), line[15:].strip()] )
			if vec[0] == 'FORMUL': 	FORMUL.append( [line[7:11].strip(), line[12:17].strip(), line[18:].strip() ] )

			if vec[0] == 'COMPND': 
				try:
					if vec[2] == 'ENGINEERED': 	ENGINEERED = True
					if vec[2] == 'MUTATION': 	MUTATION = True
				except: pass

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

		if v: print(f' >> Read :: END :: {path}/{name}')

		if save:
			self.HEADER = HEADER
			self.TITLE = TITLE
			self.COMPND = COMPND
			self.SOURCE = SOURCE
			self.AUTHOR = AUTHOR
			self.REMARK = REMARK

			self.ENGINEERED = ENGINEERED
			self.MUTATION = MUTATION
			
			self.HETNAM = HETNAM
			self.FORMUL = FORMUL

			# ==== PROTEIN ATOMS ==== # # ==== PROTEIN ATOMS ==== # # ==== PROTEIN ATOMS ==== # # ==== PROTEIN ATOMS ==== #
			self.ATOM = ATOM
													# 1 -  6        Record name   "ATOM  "
			self.ATOM_num = ATOM_num				# 7 - 11        Integer       serial       Atom  serial number.
			self.ATOM_idx = np.array(ATOM_idx) 		# 13 - 16        Atom          name         Atom name.
													# 17             Character     altLoc       Alternate location indicator.
			self.ATOM_alt = ATOM_alt				# 18 - 20        Residue name  resName      Residue name.
			self.ATOM_ani = ATOM_ani				# 22             Character     chainID      Chain identifier.
 
			self.ATOM_res   = ATOM_res				# 23 - 26        Integer       resSeq       Residue sequence number.
			self.ATOM_iCode = ATOM_iCode 			# 27             AChar         iCode        Code for insertion of residues.
			self.ATOM_xyz   = np.array(ATOM_xyz)				# 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
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

	@Logs.LogDecorator()
	def export_PDB(self, name, atoms=None, atoms_names_list=None, connectivity=None, sing=False, v=True):
		if v: print(f' >> Export as PDB >> {name}')

		atoms_names_list 	= atoms_names_list 	if not atoms_names_list is None else self.ATOM_idx
		atoms 				= atoms 			if not atoms 			is None else self.ATOM_xyz
		name  				= name if name[-4:].lower() == '.pdb' else name+'.pdb' if not name is None else self.name
		connectivity		= connectivity  	if not connectivity		is None else self.connectivity
		
		remark_dic =	 {	'HEARDER'	:	self.HEADER	, # HEADER contains the following fields: classification="PROTEIN" (or imported value), date, idCode="NONE" (or imported value).
							'TITLE'		:	self.TITLE 	, # TITLE, SOURCE, KEYWDS, EXPDTA: The imported value is exported. Default: "NULL".
							'COMPND'	:	self.COMPND	, # COMPND: The imported value is exported. Default: "MOLECULE:name", where "name" is the molecule name.
							'SOURCE'	:	self.SOURCE	, # COMPND: The imported value is exported. Default: "MOLECULE:name", where "name" is the molecule name.
														  # AUTHOR: The imported value is exported. Default: "Marvin".
							'AUTHOR'	:	 ' GENERATED BY cyclopentanoperhydrophenanthre@gmail.com  ' if sing else self.AUTHOR	, 
							'REMARK'	:	self.REMARK	, # # REMARK In entries where REMARK 0 is included as described above, REMARK 900 will also reflect the reuse of existing experimental data. 
						}

		dataout = open(f'{name}', 'w')

		for key, annotation in remark_dic.items():
			for a in annotation:
				annotation_str = a.replace('\n', '')
				print(annotation_str)
				if len(annotation_str) > 2: dataout.write("%-8s  %-70s\n" % (key, annotation_str)) # REMARK

		for i, pos in enumerate(atoms):		#loop over different atoms
			S = "%-4s  %5d %-4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6s%6s      %-4s%4s%2s\n" % ('ATOM', int(i+1), atoms_names_list[i], ' ', 
																							  'MOL', '', '1', '', pos[0], pos[1], pos[2],
																							  '1.00', '0.00', '', '', '')
			dataout.write(S) #ATOM

		if not self.connectivity is None and type(np.array(self.connectivity)) == np.ndarray and self.connectivity.shape[0] > 2:
			for c1, c in enumerate(connectivity):
				dataout.write(f'CONECT{int(c[0]):>5}' + f'{int(c[1]):>5}'*c[2] + '\n')
		else:
			for i1, pos1 in enumerate(atoms):		#loop over different atoms
				for i2, pos2 in enumerate(atoms):
					if  i1>i2 and np.linalg.norm(pos1-pos2) < 1.6:
						dataout.write(f'CONECT{int(i1+1):>5}{int(i2+1):>5}\n')

		dataout.write(f'END\n')
		dataout.close()
		return True

	def add_PDB(self, PDB2:object, save:bool=True) -> bool:
		# add new PDB data from PDB2 to self.object
		self.add_atoms( name 		= PDB2.ATOM_idx, 
						position  	= PDB2.ATOM_xyz )

		return True

	def add_atoms(self, name:str=None, position:list=None, 
						atoms_names_list:np.ndarray=None, atoms:np.ndarray=None,
						save:bool=True) -> dict:
		
		atoms_names_list 	= atoms_names_list 	 if not atoms_names_list 	is None else self.ATOM_idx
		atoms 				= np.array(atoms) 	 if not atoms 				is None else self.ATOM_xyz
		position 			= np.array(position) if type(position) 			is list else position
		name 				= np.array([name]) 	 if type(name) 				is str  else name
		
		print(atoms.shape)
		if len(atoms.shape) == 2:	atoms = np.concatenate( (atoms, position ), axis=0 ) 
		else:						atoms = np.concatenate( (atoms, position[np.newaxis,:] ), axis=0 ) 
		
		atoms_names_list = np.concatenate( (atoms_names_list, name), axis=0 ) 

		if save: 
			self.ATOM_xyz = atoms
			self.ATOM_idx = atoms_names_list
	
		return {
				'names'		:	atoms_names_list,
				'position'	:	atoms 
				}
	
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

	def is_mutated(self, ) -> bool: return self.MUTATION

	def is_engenierated(self, ): return self.ENGINEERED

	def has_heteroatom(self, name=None, match='aprox'):
		if match == 'exact':
			return any([ any([name in atom for atom in ligand]) for ligand in PDB1.FORMUL])
		elif match == 'aprox':
			return any([ any([name.lower() in atom.lower() for atom in ligand]) for ligand in PDB1.FORMUL])

	def cicle(self, distance_criteria=2, ):
		def linalg_norm(data):
		    a, b = data[0]
		    return np.linalg.norm(a - b, axis=1)

		def linalg_norm_T(data):
		    a, b = data[1]
		    return np.linalg.norm(a - b, axis=0)

		def sqrt_sum(data):
			a, b = data[0]
			return np.sqrt(np.sum((a - b) ** 2, axis=1))

		def sqrt_sum_T(data):
		    a, b = data[1]
		    return np.sqrt(np.sum((a - b) ** 2, axis=0))

		def scipy_distance(data):
		    a, b = data[0]
		    return list(map(distance.euclidean, a, b))

		def sqrt_einsum(data):
			a, b = data[0]
			a_min_b = a - b
			return np.sqrt(np.einsum("ij,ij->i", a_min_b, a_min_b))

		def exp(include, actual, last, rings=[]):
			xyz = self.ATOM_xyz[actual, :]
			for i1, n1 in enumerate(self.ATOM_xyz):
				if i1 != actual and i1 != last and np.linalg.norm(n1-xyz) < distance_criteria:
					if i1 in include:
						Sring = sorted(include[include.index(i1):])
						if not Sring in rings:
							rings.append( Sring )
					else:
						rings = exp( include=include+[i1], actual=i1, last=actual, rings=rings)
			return rings

		rings = []
		for i, n in enumerate(self.ATOM_xyz):
			rings = exp(include=[], actual=i, last=i, rings=rings)

		return rings 

	def add_center_cicle(self, path, name, allow_cicle_len:list=[5,6], 
							append:bool=True, v:bool=True):
		rings = self.cicle()
		if v: print( rings )
		if append:
			for ri, r in enumerate(rings):
				if len(r) in allow_cicle_len:
					mn = np.mean(self.ATOM_xyz[[r], :], axis=1 )[0] 
					PDB1.add_atoms( 'C', mn )
			return True

		else:
			file 	 = open(f'{path}/{name}')
			file_out = open(f'{path}/out_{name}', 'w')
			for line in file:
				vec = [m for m in line.replace('\t',' ').split(' ') if m != '' and m != '\n']

				if vec[0] == 'ENDROOT\n':# and int(vec[1]) == self.ATOM_num[-1]:
					allowed_rings_number = 0
					for ri, r in enumerate(rings):
						if len(r) in allow_cicle_len:
							allowed_rings_number += 1
							mn = np.mean(self.ATOM_xyz[[r], :], axis=1 )[0]
							S = "ATOM  %5d %2s   D10     1  %8.3f%8.3f%8.3f  1.00  0.00     0.000 AC\n" % (self.ATOM_num[-1]+allowed_rings_number, 'CM', mn[0],mn[1],mn[2])
							file_out.write(S)
				file_out.write(line)

			file_out.close()
			file.close()

			return True

	# *******************************************************************************************************************************************************************
	# * === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === *
	# *******************************************************************************************************************************************************************
def main(argv):
	# === organize arg === #
	inputfile  = argv['input']
	outputfile = argv['output']
	outputfile = 'outfiles.pdb' if outputfile is None else outputfile
	task 	   = argv['task']
	v 	  	   = True
	
	path  = ['/'.join(inputs.split('/')[:-1]) for inputs in inputfile]
	names = [inputs.split('/')[-1] for inputs in inputfile]

	if task == 'join':
		# console INPUT example 
		# python3 PDB.py -t join -i /home/akaris/Documents/temporal/denise/3spu_A_sinOH.pdb, /home/akaris/Documents/temporal/denise/3spu_A_sinOH.pdb -o denise.pdb

		# === Make data holder === #
		PDB1 = PDB(path=path[0], name=names[0])
		PDB1.READ()

		# === Make data holder === #
		PDB2 = PDB(path=path[1], name=names[1])
		PDB2.READ()

		PDB1.add_PDB(PDB2)
		PDB1.export_PDB( f'{outputfile}' )	

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-t','--task', help="task to accomplish \n   read  :  Read files from dataset  \n    summarise  :  resume data from dataset ",
	                    type=str, default='read', required=True)

	parser.add_argument('-o','--output', help="name of data output file",
	                    type=str, default=None, required=False)

	parser.add_argument('-i','--input', help="File list",
	                    type=str, default='all', nargs='+', required=False)

	args = vars(parser.parse_args())
	main(args)

'''
mypath = '/home/akaris/Documents/code/VASP/v4.6/files/POSCAR/ligands'
from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

for n in onlyfiles:
	if not 'dum' in n:
		print(n)
		PDB1 = PDB(path=mypath, name=n)
		PDB1.READ()
		
		PDB1.add_center_cicle(path=mypath, name=n)
		PDB1.export_PDB(name=mypath + '/' + n +'_exp')

print( PDB1.has_heteroatom('HOH') )
'''


