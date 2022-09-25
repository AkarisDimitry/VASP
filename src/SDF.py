# *** warning supresion
import warnings, os
warnings.filterwarnings("ignore")

# *** warning supresion
import warnings; warnings.filterwarnings("ignore")
import argparse

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

class SDF(object):
	def __init__(self, name=None, atoms=None, atoms_list=None):
		self.name = name
		self.path = None
		self.file_name = None

		self.N = None 		# N total number of atoms :: len(self.atoms)
		self.scale = None 	# scale factor

		self.atoms = atoms 			# np.array(3,N)
		self.atoms_number = None 	# [n(Fe), n(N), n(C), n(H) ]

		self.contrains = None
		
		self.atoms_names_list = None	# [Fe, N, N, N, N, C, C, C, C, H]
		self.atoms_names_ID = None  	# [Fe, N, C, H] 
		self.atoms_names_full = None 	# FeFeFeNNNNNNNCCCCCCCCCCCCCCCHHHHHHHHHHHHHHHH

		self.cell = None 
		self.distances = None # distance matrix

		self.plot_color = [
			'#DC143C', # 	crimson 			#DC143C 	(220,20,60)
			'#ADFF2F', #	green yellow 		#ADFF2F 	(173,255,47)
			'#40E0D0', #	turquoise 			#40E0D0 	(64,224,208)
			'#FF8C00', #  	dark orange 		#FF8C00 	(255,140,0)
			'#BA55D3', #	medium orchid 		#BA55D3 	(186,85,211)
			'#1E90FF', #	dodger blue 		#1E90FF 	(30,144,255)
			'#FF1493', #	deep pink 			#FF1493 	(255,20,147)
			'#8B4513', #	saddle brown 		#8B4513 	(139,69,19)
			'#FFD700', #	gold 				#FFD700 	(255,215,0)
			'#808000', #	Olive 				#808000 	(128,128,0)
			'#808080', #	Gray 				#808080 	(128,128,128)
			'#FF00FF', #	Magenta / Fuchsia 	#FF00FF 	(255,0,255)
			'#00FFFF', #	Cyan / Aqua 		#00FFFF 	(0,255,255)
			'#000000', #	Black 				#000000 	(0,0,0)
							] 

	@Logs.LogDecorator()
	def isnum(self, n:int):
		'''
		 ------------------ Define if n is or not a number ------------------ # 
		 n     :   VAR     :   VAR to check if it is a numerical VAR
		 return :  BOOL    : True/False
		1.      The float() function takes a single parameter, n, which is the number to be converted.
		2.      If n is not a number, the function returns False.
		3.      If n is a number, the function returns the float() converted value.
		'''
		try: float(n); return True
		except: return False

	@Logs.LogDecorator()
	def load(self, file_name:str=None, save:bool=True):
		# https://chem.libretexts.org/Courses/Intercollegiate_Courses/Cheminformatics_OLCC_(2019)/02%3A_Representing_Small_Molecules_on_Computers/2.05%3A_Structural_Data_Files
		if type(file_name) == str and file_name != None: 
			self.file_name = file_name

		# ------------------------ LOAD data from POSCAR ------------------------ #
		# name              :   STR     :   path+name of file to load 
		# ------------------------------- #     # ------------------------------- #
		# atom_name         :   LIST    :   List with all atoms names with-out repeat
		# atom_numbers      :   LIST    :   List with the number of atoms of each specie
		# atom_position     :   LIST    :   contains XYZ data from all atoms divided by species
		# atoms  :   N-MAT   :   Numpy array with all XYZ data
		# cell              :   N-MAT   :   Numpy array with cell parameters
		# N                 :   INT     :   Integer with total number of atoms
		# ------------------------------- #    # ------------------------------- #
		try: 	f = open(self.file_name,'r')
		except: 	print('ERROR :: POSCAR.load() :: missing file {}'.format(file_name) ); 	return

		# variables preset #
		atoms    			= [] 	# np.array(3,N)
		atoms_names_list 	= []  	# [Fe, N, N, N, N, C, C, C, C, H]
		atoms_names_ID 		= []	# [Fe, N, C, H] 
		atoms_number 		= []	# [n(Fe), n(N), n(C), n(H) ]
		contrains 			= []	# [Xbool, Ybool, Zbool]
		connectivity  		= []	# [atom0, atom1, multiplicity]
		N, C = 0, 0

		for i, n in enumerate(f):
			vec = [m for m in n.replace('\t',' ').split(' ') if m != '' and m != '\n']

			if i == 0: name = n[:-2]
			if i == 1: made_file = n[:-2]
			if i == 3: N, C = int(vec[0]), int(vec[1]) 

			if i >= 4 and i < N+4:
				atoms.append( [ float(vec[j]) for j in range(3)] )
				atoms_names_list.append( vec[3] )
				contrains.append([True, True, True] )

			if i >= N+4 and i < C+N+4:
				connectivity.append( [int(vec[j]) for j in range(3)] )

		atoms_names_ID 	 = list(set(atoms_names_list))								# [Fe, N, C, H] 
		atoms_number 	 = [atoms_names_list.count(n) for n in atoms_names_ID ] 	# [n(Fe), n(N), n(C), n(H) ]
		atoms_names_full = ''.join(atoms_names_list)								# FeFeFeNNNNNNNCCCCCCCCCCCCCCCHHHHHHHHHHHHHHHH

		cell = (np.max( atoms, axis=0 ) - np.min( atoms, axis=0 ))*10
		print( type(np.array(connectivity)) == np.ndarray )

		if save:
			self.N 					= np.array(N)
			self.C 					= np.array(C)

			self.atoms 				= np.array(atoms)
			self.atoms_number 		= np.array(atoms_number)

			self.atoms_names_list 	= np.array(atoms_names_list)
			self.atoms_names_full 	= np.array(atoms_names_full)
			self.atoms_names_ID 	= np.array(atoms_names_ID)
			self.contrains 			= np.array(contrains)

			self.connectivity		= np.array(connectivity)
			self.cell 				= np.array([[cell[0], 0, 0], [0, cell[1], 0], [0, 0, cell[2]]])

		try: f.close()
		except: 	print('ERROR :: POSCAR.load() :: can NOT close file {}'.format(file_name) ); 	return

		return True

	@Logs.LogDecorator()
	def export_PDB(self, name, atoms=None, atoms_names_list=None, connectivity=None, v=True):
		if v: print(f' >> Export as PDB >> {name}')

		atoms_names_list 	= atoms_names_list 	if not atoms_names_list is None else self.atoms_names_list
		atoms 				= atoms 			if not atoms 			is None else self.atoms
		name  				= name  			if not name  			is None else self.name
		connectivity		= connectivity  	if not connectivity		is None else self.connectivity

		dataout = open(f'{name}.pdb', 'w')
		for i, pos in enumerate(atoms):		#loop over different atoms
			S = "ATOM  %5d %2s   MOL     1  %8.3f%8.3f%8.3f  1.00  0.00\n" % (int(i+1), atoms_names_list[i], pos[0], pos[1], pos[2])
			dataout.write(S) #ATOM

		if not self.connectivity is None and type(np.array(self.connectivity)) == np.ndarray and self.connectivity.shape[0] > 2:
			for c1, c in enumerate(connectivity):
				dataout.write(f'CONECT{int(c[0]):>5}' + f'{int(c[1]):>5}'*c[2] + '\n')
		else:
			for i1, pos1 in enumerate(atoms):		#loop over different atoms
				for i2, pos2 in enumerate(atoms):
					if  i1>i2 and np.linalg.norm(pos1-pos2) < 1.6:
						dataout.write(f'CONECT{int(i1+1):>5}{int(i2+1):>5}\n')

		dataout.close()
		return True
			
	@Logs.LogDecorator()
	def plot(self, save=False, path='./no_name.png'):
		if type(self.atoms) != type(None):
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

	@Logs.LogDecorator()
	def summary(self, ):
		POSCAR_resumen = self.resume()

		if type(self.file_name) == str and 'CONTCAR' in self.file_name:		toprint = '            *--> [CONTCAR] ::\n'
		if type(self.file_name) == str and 'POSCAR' in self.file_name: 		toprint = '            *--> [POSCAR]  ::\n'
		else:																toprint = '            *--> [P/C CAR] ::\n'

		for osr in POSCAR_resumen['list']:	toprint += '    |         \t{}:{}\n'.format(osr[0], osr[1])
		print(toprint)
		
		return None

	@Logs.LogDecorator()
	def resume(self, v=0):
		try:
			info_POSCAR = { 	'CELL': ';  '.join([ ','.join([ '{:.2f}'.format(y) for y in x]) for  x in self.cell ]), }
			self.relevant_distances()
			for i, (key_atom, value_atom) in enumerate( self.relevant_distances_dict_extend.items()):
				info_POSCAR[key_atom] = value_atom
			#'Distances': '; '.join([ '; '.join([ '; '.join([ '{}{}-{}{}:{:.3}'.format(key_atom, key_atom2, value_atom2['names'][k], value_atom2['index'][k], dist) for k, dist in enumerate(value_atom2['distances'])]) for j, (key_atom2, value_atom2) in enumerate( value_atom.items())]) for i, (key_atom, value_atom) in enumerate( self.relevant_distances(v=0).items()) ]) }
		except: info_POSCAR = {}
		
		return {'list':[[key, value] for i, (key, value) in enumerate( info_POSCAR.items())], 'dict':info_POSCAR}


	def direct_to_cartesian(self, *arg, **karg):
		return True

	@Logs.LogDecorator()
	def operations(self, operation={'name':'shake', 'intencity':0.1}, 
				atoms=None, atoms_number=None, atoms_names_ID=None, atoms_names_list=None, atoms_names_full=None,
				cell=None, contrains=None, save=False, v=True):
		atoms_names_ID = np.array(atoms_names_ID) if type(atoms_names_ID) != type(None) else np.array(self.atoms_names_ID)
		atoms_names_list = np.array(atoms_names_list) if type(atoms_names_list) != type(None) else np.array(self.atoms_names_list)
		
		atoms = np.array(atoms) if type(atoms) != type(None) else np.array(self.atoms)
		atoms_number = np.array(atoms_number) if type(atoms_number) != type(None) else np.array(self.atoms_number)
		atoms_names_full = atoms_names_full if not atoms_names_full is None else self.atoms_names_full
		
		contrains = np.array(contrains) if type(contrains) != type(None) else np.array(self.contrains)
		
		cell = np.array(cell) if type(cell) != type(None) else np.array(self.cell)

		atom_N = atoms.shape[0]

		def random_shake(intencity):	
			return atoms + np.random.rand(atom_N, 3)*intencity

		def move(atoms2move, distance, ):
			displace = np.zeros_like(atoms)
			displace[atoms2move,:] += distance
			return atoms + displace

		def replicate(R):
			if v: print(f' >> Replicate {atoms.shape[0]}:ATOMS >> ({R[0]}, {R[1]}, {R[2]}) >> {atoms.shape[0]*R[0]*R[1]*R[2]}:ATOMS')
			#self.atoms      		= atoms # np.array(3,N)
			#self.atoms_number 		= None 	# [n(Fe), n(N), n(C), n(H) ]	
			#self.atoms_names_list 	= None	# [Fe, N, N, N, N, C, C, C, C, H]
			#self.atoms_names_ID 	= None  # [Fe, N, C, H] 
			#self.atoms_names_full 	= None 	# FeFeFeNNNNNNNCCCCCCCCCCCCCCCHHHHHHHHHHHHHHHH
			rep_atoms_names_list 	= []
			rep_atoms 				= []
			rep_atoms_number 		= [ an*R[0]*R[1]*R[2] for an in atoms_number ]	 
			rep_atoms_names_full = ''.join([aid*int(an) for an, aid in zip(rep_atoms_number, atoms_names_ID) ])
			rep_contrains = []

			for i, atom in enumerate(atoms):
				for Rx in range(R[0]):
					for Ry in range(R[1]):
						for Rz in range(R[2]):
							rep_atoms.append( np.array(atom + Rx*cell[0] + Ry*cell[1] + Rz*cell[2] ) )
							rep_atoms_names_list.append( atoms_names_list[i] )
							rep_contrains.append(contrains[i])
			rep_cell = np.array([R[0]*cell[0], R[1]*cell[1], R[2]*cell[2]])

			return np.array(rep_atoms), rep_atoms_names_list, np.array(rep_contrains), rep_atoms_names_full, rep_atoms_number, rep_cell

		if operation['name'] == 'shake':
			atoms = random_shake(operation['intencity'])

		if operation['name'] == 'move':
			atoms = move(operation['atoms2move'], operation['distance'])

		if operation['name'] == 'replicate':
			atoms, atoms_names_list, contrains, atoms_names_full, atoms_number, cell= replicate(operation['replicate'])

		if save: 
			self.atoms_names_ID = atoms_names_ID
			self.atoms_names_list = atoms_names_list
			self.atoms = atoms
			self.atoms_number = atoms_number
			self.atoms_names_full = atoms_names_full

			self.contrains = contrains

			self.cell = cell

		return atoms

	@Logs.LogDecorator()
	def get_ase(self, ):
		return Atoms(	self.atoms_names_list,
			 positions=	self.atoms,
			 cell     =	self.cell,
			 pbc      =	[1, 1, 0])

	@Logs.LogDecorator()
	def get_distance_matrix(self, atoms_position=None, save:bool=True ):
		N = self.N
		cell = self.cell
		pos = self.atoms[:N, :] if type(atoms_position) == type(None) else atoms_position 

		# === arange position data === #
		v1a = np.outer(np.ones(N), pos[:,0])
		v2a = np.outer(np.ones(N), pos[:,1])
		v3a = np.outer(np.ones(N), pos[:,2])

		# === arange position data === #
		v1b = np.outer( pos[:,0], np.ones(N) )
		v2b = np.outer( pos[:,1], np.ones(N) )
		v3b = np.outer( pos[:,2], np.ones(N) )

		# === Evaluate distance per component === #
		dif1 = v1a - v1b 
		dif2 = v2a - v2b 
		dif3 = v3a - v3b

		# === Chech first neighbour  === #
		condition1 = np.abs(dif1) > cell[0][0]/2
		condition2 = np.abs(dif2) > cell[1][1]/2
		condition3 = np.abs(dif3) > cell[2][2]/2

		dif1_edited = dif1  + condition1*cell[0][0]
		dif2_edited = dif2  + condition2*cell[1][1]
		dif3_edited = dif3  + condition3*cell[2][2]

		# === First neighbour distance VECTOR  === #
		vector_dist = np.array([dif1_edited, dif2_edited, dif3_edited])
		dist = np.linalg.norm(vector_dist, axis=0)
		# === First neighbour distance VERSOR  === #
		versor_dist = vector_dist / dist

		'''
		# === Filter atoms by distance distances  === #
		dist_filter = np.tril(dist) # amtrix triangular inferior
		filt = np.logical_and(dist_filter<cut_dist[1], dist_filter>cut_dist[0])
	
		# === Calcular angulos === #
		unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
		unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
		dot_product = np.dot(unit_vector_1, unit_vector_2)
		angle = np.arccos(dot_product)
		'''

		if save:
			self.distance_matrix = dist
			self.versor_dist = versor_dist

		return dist

	# *******************************************************************************************************************************************************************
	# * === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === ARGPARSE === *
	# *******************************************************************************************************************************************************************
def main(argv):
	# === organize arg === #
	inputfile  		= argv['input']
	outputfile 		= argv['output']
	outputfile 		= inputfile+'.o' if type(outputfile) == type(None) else outputfile
	task 	   		= argv['task']
	replicateX 	    = argv['replicateX']
	replicateY 	    = argv['replicateY']
	replicateZ 	    = argv['replicateZ']
	incell 	   		= True if argv['incell'] in ['T', 'True', 'TRUE', 1, '1'] else False
	output_format 	= argv['format']

	# === Make data holder === #
	sdf = SDF()
	# load file #
	print(f'Input  >> {inputfile} ')
	path = '/'.join(inputfile.split('/')[:-1])
	sdf.load(file_name = inputfile)

	if replicateX+replicateY+replicateZ > 3:
		path = '/'.join(inputfile.split('/')[:-1])
		sdf.load(file_name = inputfile)
		sdf.operations( operation={'name':'replicate', 'replicate':[replicateX,replicateY,replicateZ]}, save=True )

	if incell:
		sdf.atoms = sdf.get_molecule(dist_filter=2.0,)

	if task == 'export':

		if   output_format == 'POSCAR':		sdf.export(outputfile)
		elif output_format == 'PDB':		sdf.export_PDB(outputfile)



	print(f'Input  >> {inputfile} ')
	print(f'OUTPUT >> {outputfile}')

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-t','--task', help="task to accomplish :: read-Read files from dataset :: summarise-resume data from dataset :: export-export data ",
	                    type=str, default='export', required=False)
	
	parser.add_argument('-i','--input', help="[STR] path to inputfile || eg. -i ./",
	                    type=str, default='', required=True)

	parser.add_argument('-o','--output', help="name of data output file",
	                    type=str, default=None, required=False)

	parser.add_argument('-rx','--replicateX', help="[INT] number of cell replica in y axes || eg. -rx 1",
	                    type=int, default=1, required=False)

	parser.add_argument('-ry','--replicateY', help="[INT] number of cell replica in y axes || eg. -ry 1",
	                    type=int, default=1, required=False)

	parser.add_argument('-rz','--replicateZ', help="[INT] number of cell replica in y axes || eg. -rz 1",
	                    type=int, default=1, required=False)

	parser.add_argument('-incell','--incell', help="[BOOL] join atoms into one molecule || eg. -ry TRUE",
	                    type=str, default='False', required=False)

	parser.add_argument('-f','--format', help="[STR] number of cell replica in y axes || eg. -f PDB -f POSCAR",
	                    type=str, default='POSCAR', required=False)

	parser.add_argument('-v','--verbosity', help="verbosity",
	                    type=str, default=1, required=False)
	
	# ---- cookbook ---- #
	# python3 
	args = vars(parser.parse_args())
	main(args)

