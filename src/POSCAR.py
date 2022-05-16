# *** warning supresion
import warnings, os
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

class POSCAR(object):
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

	def isnum(self, n):
		# ------------------ Define if n is or not a number ------------------ # 
		# n     :   VAR     :   VAR to check if it is a numerical VAR
		# return :  BOOL    : True/False
		try: float(n); return True
		except: return False

	def load(self, file_name=None, save=True):
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
		cell 				= np.zeros((3,3)) 
		N 					= 0 
		
		atoms 				= []
		atoms_number 		= [] 

		atoms_names_list 		= [] 
		atoms_names_ID 		= [] 
		atoms_names_full 	= ''

		direct 				= 0
		scale 				= 0
		contrains 			= []

		for i, n in enumerate(f):
			vec = [m for m in n.replace('\t',' ').split(' ') if m != '' and m != '\n']
			try: 
				if vec[-1][-1:] == '\n': vec[-1] = vec[-1][:-1]
			except: pass

			# **** LOAD SCALE FACTOR**** #
			if i == 1: scale = float(vec[0])

			# **** LOAD CELL **** #
			for m in range(3):
				if i == 2+m:  cell[m,0]=float(vec[0]);  cell[m,1]=float(vec[1]); cell[m,2]=float(vec[2])

			# **** LOAD ATOMS names (without repetition) **** #
			if i == 5:
				for m in vec: atoms_names_ID.append(m)

			# **** LOAD ATOMS names (with repetition) **** #
			if i == 6:
				for j, m in enumerate(vec):  
					N += int(m); 									# total number of atoms
					atoms_number.append( int(m) );  					# number of each specie
					atoms_names_full += str(atoms_names_ID[j])*int(m)
					atoms_names_list += [atoms_names_ID[j]]*int(m)
				atoms = np.zeros((N, 3))

			# **** LOAD ATOMS positions **** #	
			if i == 8 and vec[0]=='Direct': 		direct = True
			elif i == 8 and not vec[0]=='Direct':	direct = False

			# **** LOAD ATOMS contrains **** #	
			if i >= 9:
				if i < N+9: 
					atoms[i-9,:] = np.array( [vec[0], vec[1], vec[2]] )
					contrains.append([True if 'T' in n else False for n in [vec[3], vec[4], vec[5]]] )

		# ---- store loaded data in OBJ ---- #
		if save:
			self.cell 				= np.array(cell)
			self.N 					= N
			self.scale 				= scale

			self.atoms 				= np.array(atoms)
			self.atoms_number 		= atoms_number

			self.contrains 			= contrains

			self.atoms_names_list 		= atoms_names_list
			self.atoms_names_full 	= atoms_names_full
			self.atoms_names_ID 	= atoms_names_ID

			if direct: 	self.coordenate = 'Direct'
			else:		self.coordenate = 'Cartesian'
		
		return 	np.array(cell), N, np.array(atoms), atoms_number,  contrains, atoms_names_list, atoms_names_full, atoms_names_ID

	def direct_to_cartesian(self, atoms=None, cell=None, save=True, force=False, criteria=False, v=True):
		if not type(atoms) == np.ndarray: 	atoms = self.atoms
		elif type(atoms) == list:			atoms = np.array(atoms)

		if not type(cell) == np.ndarray: 	cell = self.cell
		elif type(cell) == list:			cell = np.array(cell)

		criteria = True if (atoms<= 1.0).all()  else False

		if self.coordenate == 'Direct' or force or criteria:
			try:
				atoms = np.dot(atoms, cell)	
			except: 
				if v: print('ERROR :: POSCAR.direct_to_cartesian() :: can not convert direct to cartesian' )
				return

			if save: self.atoms = atoms
			self.coordenate = 'Cartesian'
		else:
			pass
			#if v: print('WARNNING :: POSCAR.direct_to_cartesian() :: can not convert direct to cartesian' )

		return atoms

	def operations(self, operation={'name':'shake', 'intencity':0.1}, 
				atoms=None, atoms_number=None, atoms_names_ID=None, atoms_names_list=None, cell=None, save=False):
		atoms_names_ID = np.array(atoms_names_ID) if type(atoms_names_ID) != type(None) else np.array(self.atoms_names_ID)
		atoms_names_list = np.array(atoms_names_list) if type(atoms_names_list) != type(None) else np.array(self.atoms_names_list)
		
		atoms = np.array(atoms) if type(atoms) != type(None) else np.array(self.atoms)
		atoms_number = np.array(atoms_number) if type(atoms_number) != type(None) else np.array(self.atoms_number)
		
		cell = np.array(cell) if type(cell) != type(None) else np.array(self.cell)

		atom_N = atoms.shape[0]

		def random_shake(intencity):	
			return atoms + np.random.rand(atom_N, 3)*intencity

		def move(atoms2move, distance, ):
			displace = np.zeros_like(atoms)
			displace[atoms2move,:] += distance
			return atoms + displace

		if operation['name'] == 'shake':
			atoms = random_shake(operation['intencity'])

		if operation['name'] == 'move':
			atoms = move(operation['atoms2move'], operation['distance'])

		if save: 
			self.atoms_names_ID = atoms_names_ID
			self.atoms_names_list = atoms_names_list
			self.atoms = atoms
			self.atoms_number = atoms_number
			self.cell = cell

		return atoms

	def data_augmentation(self, path=None, operation={'name':'shake', 'intencity':0.1}, size=10,
						atoms=None, atoms_number=None, atoms_names_ID=None, atoms_names_list=None, cell=None, ):

		atoms_names_ID = np.array(atoms_names_ID) if type(atoms_names_ID) != type(None) else np.array(self.atoms_names_ID)
		atoms_names_list = np.array(atoms_names_list) if type(atoms_names_list) != type(None) else np.array(self.atoms_names_list)
		
		atoms = np.array(atoms) if type(atoms) != type(None) else np.array(self.atoms)
		atoms_number = np.array(atoms_number) if type(atoms_number) != type(None) else np.array(self.atoms_number)
		
		cell = np.array(cell) if type(cell) != type(None) else np.array(self.cell)

		atom_N = atoms.shape[0]

		self.make_subdirectories( path )

		for n in range(size):
			self.make_subdirectories( path+'/sample_{}'.format(n) )
			atoms = self.operations( operation=operation, atoms=atoms, atoms_number=atoms_number, 
									atoms_names_ID=atoms_names_ID, atoms_names_list=atoms_names_list, cell=cell, save=False)
			self.export(name_export=path+'/sample_{}/POSCAR'.format(n), atoms=atoms)

	def make_subdirectories(self, path):
		# define the name of the directory to be created
		try:
		    os.makedirs(path)
		except OSError:
		    print ("Creation of the directory %s failed" % path)
		else:
		    print ("Successfully created the directory %s" % path)

	def get_postion(self, ion, ):
		# gives all the position of ion in self.atoms_names_list  
		# ------------------------------------------------------
		# ion  ::	STR 	:: ion name str 
		if type(self.atoms_names_list) == list:
			return list(filter(lambda x: (x >= 0), [ i if n == ion else -1 for i, n in enumerate( self.atoms_names_list ) ] ))
		else:
			return []
			
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

	def plot_distance(self, ax=None, save=False, path='./no_name.png'):
		if type(self.atoms) != type(None):
			self.relevant_distances()
			if ax==None:	fig, ax = plt.subplots()
			else:			fig = ax.get_figure()

			ref_dic = {}
			# == grab reference == #
			for i, (key_atoms, value_atoms) in enumerate( self.relevant_distances_dict.items()):
				for j, (key_atom, value_atom) in enumerate(value_atoms.items()):
					for k, name in enumerate(set(value_atom['names'])): 
						if not f'{key_atoms}{name}' in ref_dic:
							ref_dic[f'{key_atoms}{name}'] = len(ref_dic)

			# == plot scatter == #
			for i, (key_atoms, value_atoms) in enumerate( self.relevant_distances_dict.items()):
				for j, (key_atom, value_atom) in enumerate(value_atoms.items()):
					for k, name in enumerate(value_atom['names']): 
						
						ax.plot( ref_dic[f'{key_atoms}{name}'], value_atom['distances'][k], 'o', color=self.plot_color[ref_dic[f'{key_atoms}{name}']] )
						ax.annotate(f'{key_atoms}-{name}', (ref_dic[f'{key_atoms}{name}'], value_atom['distances'][k]))

					#legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')
			ax.set_title('Distances')
			ax.set_ylabel('distance (A)')
			ax.set_xlabel('specie')

			if save:
				fig.savefig(f'{path}' , dpi=100, pad_inches=0.1, bbox_inches='tight', horizontalalignment='right') 

	def plot_embedding(self, atomic_embedding=None, atoms_names_ID=None, atoms_names_list=None):
		atoms_names_ID = np.array(atoms_names_ID) if type(atoms_names_ID) != type(None) else np.array(self.atoms_names_ID)
		atoms_names_list = np.array(atoms_names_list) if type(atoms_names_list) != type(None) else np.array(self.atoms_names_list)
		atomic_embedding = np.array(atomic_embedding) if type(atomic_embedding) != type(None) else np.array(self.atomic_embedding)

		atomic_embedding_sum = np.sum(atomic_embedding, axis=2)
		plt.hist(atomic_embedding_sum, )

	def export(self, name_export=None, 
				atoms=None, atoms_number=None, atoms_names_ID=None, atoms_names_list=None, cell=None, ):

		atoms_names_ID = np.array(atoms_names_ID) if type(atoms_names_ID) != type(None) else np.array(self.atoms_names_ID)
		atoms_names_list = np.array(atoms_names_list) if type(atoms_names_list) != type(None) else np.array(self.atoms_names_list)
		
		atoms = np.array(atoms) if type(atoms) != type(None) else np.array(self.atoms)
		atoms_number = np.array(atoms_number) if type(atoms_number) != type(None) else np.array(self.atoms_number)
		
		cell = np.array(cell) if type(cell) != type(None) else np.array(self.cell)

		if not type(name_export) == str:
			print('ERROR :: POSCAR.export() :: need a valid name ')

		file = open(name_export, 'w')

		file.write('Comment_line \n') # coment line
		file.write('{}\n'.format(self.scale) ) # factor scale

		file.write('{:>18.15f}\t{:>18.15f}\t{:>18.15f}\n'.format(cell[0][0], cell[0][1], cell[0][2]) )	# cell vector 00
		file.write('{:>18.15f}\t{:>18.15f}\t{:>18.15f}\n'.format(cell[1][0], cell[1][1], cell[1][2]))	# cell vector 01
		file.write('{:>18.15f}\t{:>18.15f}\t{:>18.15f}\n'.format(cell[2][0], cell[2][1], cell[2][2]))	# cell vector 02
		file.write('    '.join(atoms_names_ID) +' \n' )	# atom type list
		file.write('    '.join( [str(n) for n in atoms_number] ) +' \n' )	# number of each atom type

		file.write('Selective dynamics \n')	# selective dynamics

		if self.coordenate in ['direct', 'Direct'] :		file.write('Direct\n')	# direct
		elif self.coordenate in ['cartesian', 'Cartesian']:	file.write('Cartesian\n')	# direct
		else: 												file.write('Direct\n')	# direct

		for i, atom in enumerate(atoms):
			file.write(	'\t'+'\t'.join( [ '{:>18.15f}'.format(n) for n in atom] ) +
						'\t'+'  '.join( [ 'T' if n else 'F' for n in self.contrains[i] ] ) +
						'\n' )

		file.write('Comment_line \n')	# atom N

	def inyect_data(self, obj):
		try: self.__dict__ = obj.__dict__.copy()
		except: print(' ERROR  :: code X :: DATA.inyect_data() :: can not inject data into {}'.format(str(obj)))

	def periodic_distance_matrix(self, ):
		N = self.N
		cell = self.cell
		pos = self.atoms

		v1a = np.dot(np.ones(N)[:,np.newaxis], pos[:,0,np.newaxis].T)
		v2a = np.dot(np.ones(N)[:,np.newaxis], pos[:,1,np.newaxis].T)
		v3a = np.dot(np.ones(N)[:,np.newaxis], pos[:,2,np.newaxis].T)

		v1b = np.dot( pos[:,0,np.newaxis], np.ones(N)[np.newaxis, :] )
		v2b = np.dot( pos[:,1,np.newaxis], np.ones(N)[np.newaxis, :] )
		v3b = np.dot( pos[:,2,np.newaxis], np.ones(N)[np.newaxis, :] )

		dif1 = np.abs(v1a - v1b) 
		dif2 = np.abs(v2a - v2b) 
		dif3 = np.abs(v3a - v3b) 

		condition1 = dif1 > cell[0][0]/2
		condition2 = dif2 > cell[1][1]/2
		condition3 = dif3 > cell[2][2]/2

		dif1_edited = dif1  - condition1*cell[0][0]
		dif2_edited = dif2  - condition2*cell[1][1]
		dif3_edited = dif3  - condition3*cell[2][2]

		dist = (dif1_edited**2 + dif2_edited**2 + dif3_edited**2)**0.5

		return dist

	def internal_distance(self, metric='norm2', atoms=None):
		return np.array([ np.linalg.norm((self.atoms - n[np.newaxis, :]), axis=1) for n in self.atoms[self.has(atoms)] ])

	def periodic_distance(self, metric='norm2', atoms=None):
		distance_list = []

		if not type(self.atoms) == type(None):	# check if atoms position have been loaded

			if type(atoms) == list and len(atoms)>0 and type(atoms[0]) == str:
				# atoms = ['Fe', 'O']
				for n in self.atoms[self.has(atoms)]:
					a = self.atoms - n[np.newaxis, :] 
					a[ a[:,0] > self.cell[0,0]/2] -= [self.cell[0,0], 0, 0]
					a[ a[:,1] > self.cell[1,1]/2] -= [0, self.cell[1,1], 0]
					a[ a[:,2] > self.cell[2,2]/2] -= [0, 0, self.cell[2,2]]

					a[ a[:,0] < -self.cell[0,0]/2] += [self.cell[0,0], 0, 0]
					a[ a[:,1] < -self.cell[1,1]/2] += [0, self.cell[1,1], 0]
					a[ a[:,2] < -self.cell[2,2]/2] += [0, 0, self.cell[2,2]]

					distance_list.append( np.linalg.norm(a, axis=1) )

			if type(atoms) == int:
				# atoms = 1
				a = self.atoms - self.atoms[atoms][np.newaxis, :] 

				a[ a[:,0] > self.cell[0,0]/2] -= [self.cell[0,0], 0, 0]
				a[ a[:,1] > self.cell[1,1]/2] -= [0, self.cell[1,1], 0]
				a[ a[:,2] > self.cell[2,2]/2] -= [0, 0, self.cell[2,2]]

				a[ a[:,0] < -self.cell[0,0]/2] += [self.cell[0,0], 0, 0]
				a[ a[:,1] < -self.cell[1,1]/2] += [0, self.cell[1,1], 0]
				a[ a[:,2] < -self.cell[2,2]/2] += [0, 0, self.cell[2,2]]

				distance_list = np.linalg.norm(a, axis=1)

			return np.array(distance_list), np.array(self.atoms_names_list)[self.has(atoms)]

		else:	return None, None

	def closest(self, atom, distances_filter='min'):
		# return distance, atoms index
		distances, names = self.periodic_distance(metric='norm2', atoms=atom) 
		
		if distances_filter == 'min':
			if type(distances) != type(None):
				return np.min(distances[distances>0]), np.argmin(distances[distances>0])
			else:
				return None, None

		elif self.isnum(distances_filter) :
			if type(distances) != type(None):
				return distances[distances<distances_filter], np.arange( len(self.atoms_names_list) )[distances<distances_filter]
			else:
				return None, None

		else:
			print('ERROR :: can NOT identify "distances_filter" parameter. try use: distances_filter="min",  distances_filter=2.1')
			return None

	def relevant_atoms_position(self, relevance=2, save=True, v=0):
		relevant_mask = self.relevant_atoms_mask()
		return self.atoms[relevant_mask]

	def relevant_atoms_mask(self, relevance=3, save=True, v=0):
		relevant_mask = [ self.atoms_names_list.count(specie) <= relevance for specie in self.atoms_names_list ]
		if save: self.relevant_mask = relevant_mask
		return relevant_mask

	def relevant_distances(self, relevance=2, relevance_distance=2.5, save=True, v=0):
		# 'relevance' is the maximum number of atoms of the n-th species allowed in order to be considered as a relevant specie
		relevant_distances_dict = {} 
		if type(self.atoms_names_ID) != type(None) and type(self.atoms_names_list) != type(None) and  type(self.atoms) != type(None):
			self.direct_to_cartesian( criteria=True )
			for atom_name in self.atoms_names_ID:
				if self.atoms_names_list.count(atom_name) <= relevance:
					relevant_distances_dict[atom_name] = {}
					index_list = list(filter(lambda x: (x >= 0 ), [index if specie == atom_name else -1 for index, specie in enumerate(self.atoms_names_list)]))
					for atom_index in index_list:
						distances, distances_index = self.closest(atom=atom_index, distances_filter=relevance_distance) 
						distances_index = distances_index[distances>0]
						distances = distances[distances>0]
						relevant_distances_dict[atom_name][atom_index] = {'distances': distances, 'index':distances_index, 'names':[self.atoms_names_list[n] for n in distances_index]}
		else:
			print('WARNNING :: POSCAR.relevant_distances() :: type(self.atoms_names_ID) != type(None) and type(self.atoms_names_list) != type(None) and  type(self.atoms) != type(None)' )

		if v >= 1:
			print(' Relevant species: ' )
			for i, n in enumerate(relevant_distances_dict.keys()):
				print( ' ({0}) {1} ( in system {2} )'.format(i, n, len(relevant_distances_dict[n].keys()) ) )
				for j, m in enumerate(relevant_distances_dict[n]):
					print( ' **** {0} idex {1} '.format(n, m) )
					for k, o in enumerate(relevant_distances_dict[n][m]['index']):
						print( ' **** **** d({0}_{1}, {2}_{3} \t : \t {4:0.3})  '.format(n, m, relevant_distances_dict[n][m]['names'][k], relevant_distances_dict[n][m]['index'][k], relevant_distances_dict[n][m]['distances'][k]) )

		if save: 
			self.relevant_distances_dict = relevant_distances_dict
			self.relevant_distances_dict_extend = {'{}{}-{}{}:'.format(key_atom, key_atom2, value_atom2['names'][k], value_atom2['index'][k]) : '{:.3}'.format(dist)   for i, (key_atom, value_atom) in enumerate( self.relevant_distances_dict.items()) for j, (key_atom2, value_atom2) in enumerate( value_atom.items()) for k, dist in enumerate(value_atom2['distances'])}
			self.relevant_distances_text = '; '.join([ '; '.join([ '; '.join([ '{}{}-{}{}:{:.3}'.format(key_atom, key_atom2, value_atom2['names'][k], value_atom2['index'][k], dist) for k, dist in enumerate(value_atom2['distances'])]) for j, (key_atom2, value_atom2) in enumerate( value_atom.items())]) for i, (key_atom, value_atom) in enumerate( self.relevant_distances_dict.items() ) ])

		return relevant_distances_dict

	def isnum(self, n):
		# ------------------ Define if n is or not a number ------------------ # 
		# n     :   VAR     :   VAR to check if it is a numerical VAR
		# return :  BOOL    : True/False
		try: float(n); return True
		except: return False

	def has(self, atoms, v=False):	
		try:		return np.array([ True if n in atoms else False for n in self.atoms_names_list])
		except: 	
			if v: print('ERROR :: POSCAR.has :: can NOT get atoms from self.atoms_names_list')	

	def summary(self, ):
		POSCAR_resumen = self.resume()

		if type(self.file_name) == str and 'CONTCAR' in self.file_name:		toprint = '            *--> [CONTCAR] ::\n'
		if type(self.file_name) == str and 'POSCAR' in self.file_name: 		toprint = '            *--> [POSCAR]  ::\n'
		else:																toprint = '            *--> [P/C CAR] ::\n'

		for osr in POSCAR_resumen['list']:	toprint += '    |         \t{}:{}\n'.format(osr[0], osr[1])
		print(toprint)
		
		return None

	def resume(self, v=0):
		try:
			info_POSCAR = { 	'CELL': ';  '.join([ ','.join([ '{:.2f}'.format(y) for y in x]) for  x in self.cell ]), }
			self.relevant_distances()
			for i, (key_atom, value_atom) in enumerate( self.relevant_distances_dict_extend.items()):
				info_POSCAR[key_atom] = value_atom
			#'Distances': '; '.join([ '; '.join([ '; '.join([ '{}{}-{}{}:{:.3}'.format(key_atom, key_atom2, value_atom2['names'][k], value_atom2['index'][k], dist) for k, dist in enumerate(value_atom2['distances'])]) for j, (key_atom2, value_atom2) in enumerate( value_atom.items())]) for i, (key_atom, value_atom) in enumerate( self.relevant_distances(v=0).items()) ]) }
		except: info_POSCAR = {}
		
		return {'list':[[key, value] for i, (key, value) in enumerate( info_POSCAR.items())], 'dict':info_POSCAR}

	def get_first_neighbour_distance_matrix(self, atoms_position=None, save=True ):

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

		dif1_edited = dif1  - condition1*cell[0][0]
		dif2_edited = dif2  - condition2*cell[1][1]
		dif3_edited = dif3  - condition3*cell[2][2]

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

	def get_embedding(self, distance_matrix=None, atoms_names_ID=None, atoms_names_list=None, embedding={'partition':200, 'min distance':0.5, 'max distance':3.0, 'range distance':0.05}, save=True):

		def gaussian(mu, sigma, n, ): return np.e**(-1.0/2 * ((mu-np.arange(n))/sigma)**2) * 0.3989422804/sigma

		# == define some variables == #
		Dmin, Dmax ,Dsig,N = embedding['min distance'], embedding['max distance'], embedding['range distance'], embedding['partition']
		Nsig = Dsig * N / (Dmax-Dmin)

		atoms_names_ID = np.array(atoms_names_ID) if type(atoms_names_ID) != type(None) else np.array(self.atoms_names_ID)
		atoms_names_list = np.array(atoms_names_list) if type(atoms_names_list) != type(None) else np.array(self.atoms_names_list)
		distance_matrix = np.array(distance_matrix) if type(distance_matrix) != type(None) else np.array(self.distance_matrix)
		# == asing memory space == #
		# atomic_embedding = { i1:{a2:np.zeros(N) for a2 in atoms_names_in} for i1, a1 in enumerate(atoms_names) } #
		atomic_embedding = np.zeros( (atoms_names_list.shape[0], atoms_names_ID.shape[0], N) )

		# == Filter by atomic type == #
		atomic_filter = np.zeros( (atoms_names_list.shape[0], atoms_names_list.shape[0], atoms_names_ID.shape[0])  )
		for i, a in enumerate(atoms_names_ID): atomic_filter[atoms_names_list == a, atoms_names_list == a, i] = 1

		# == Filter by distance == #
		dist_emb = ((distance_matrix - Dmin) / (Dmax-Dmin) * N).astype(int)
		dist_emb[ dist_emb<0 ] = 0
		dist_emb[ dist_emb>N+Nsig*3 ] = 0

		# == Apply filters == #
		v1a = np.tensordot(dist_emb, atomic_filter, axes=1)

		# == Filter by distance == #
		for i1, a1 in enumerate(atoms_names_list):
			for i2, a2 in enumerate(atoms_names_ID):
				unique, counts = np.unique(v1a[i1,:,i2], return_counts=True) 
				for u, c in zip(unique[1:], counts[1:]):
					atomic_embedding[i1][i2][:] += gaussian(u, Nsig, N)*c

		if save: self.atomic_embedding = atomic_embedding
		return atomic_embedding # [ atoms , ID, N ]

	def get_ase(self, ):
		return Atoms(	self.atoms_names_list,
             positions=	self.atoms,
             cell     =	self.cell,
             pbc      =	[1, 1, 0])


'''
# *** EG - DATA augmentation *** #
pos = POSCAR()
pos.load(file_name = '/home/akaris/Documents/code/VASP/v4.1/files/POSCAR/FePC/FePC_O2_Au/POSCAR')
pos.summary() 
pos.plot_distance()
plt.show()
'''


'''
pos.plot(save=True)

pos.direct_to_cartesian()
#pos.export( 'POSCAR' )
#eroro
pos.data_augmentation(	path='/home/akaris/Documents/temporal/Free_standing/FePC_O2', 
						operation={	'name':'move', 
									'atoms2move':[True if n in [147, 148] else False for n in range( len(pos.atoms) )], 
									'distance':[0, 0, +0.10]},
						size=10)
'''

'''
pos.data_augmentation(path='/home/akaris/Documents/temporal/Free_standing/FePC', operation={'name':'shake', 'parameters':{'intencity':0.02}})

# *** EG *** #
# eg. relevant_distances
pos = POSCAR()
pos.load(file_name = '/home/akaris/Documents/code/VASP/v3.5/files/POSCAR/FeBeBeAu/FePC_5Bz_OH')
#pos.periodic_distance_matrix()
pos.direct_to_cartesian()

dist = pos.get_first_neighbour_distance_matrix()
embedding = pos.get_embedding()
pos.plot_embedding()
plt.show()

pos.plot( name_file= 'files/POSCAR' )
pos.export(name_export='afsldkmlasd')
'''