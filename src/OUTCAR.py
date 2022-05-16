# *** warning supresion
import warnings
warnings.filterwarnings("ignore")

# *** warning supresion
import warnings; warnings.filterwarnings("ignore")

# *** numpy libraries
try:
	import numpy as np
except: print('WARNING :: DATA.import_libraries() :: can not import numpy ')

# *** matplotlib libraries
try:
	import matplotlib.pyplot as plt
	import matplotlib.axes as ax
	import matplotlib.patches as patches
except:	print('WARNING :: Set.import_libraries() :: can not import matplotlib ')


class OUTCAR(object):
	def __init__(self, name=None, ):
		self.name = name
		self.file_name = None

		self.POTCARs = None
		self.POTCAR = None
		self.POSCAR = None
		self.E_fermi = None
		self.NIONS = None

		self.atoms_names_list = None	# [Fe, N, N, N, N, C, C, C, C, H]
		self.atoms_names_ID = None  	# [Fe, N, C, H] 

		self.E = None#

		self.orbitals = {'s':0, 'px':1, 'py':2, 'pz':3,
					'dxy':4, 'dyz':5, 'dz2':6, 'dxz':7, 'dx2':8, 'tot':9 } #    s     py     pz     px    dxy    dyz    dz2    dxz    dx2

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
							# ------- REPEAT -------- #
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
 


	def load(self, file_name=None):
		if file_name != None: self.file_name = file_name

		try: 	f = open(self.file_name,'r')
		except: 	print('ERROR :: OUTCAR.load() :: missing file {}'.format(file_name) ); 	return

		self.POTCAR = []# POTCAR information 
		self.POTCAR_full = []# POTCAR information 
		self.POSCAR = ''# POSCAR name
		self.atoms_names_list = []	# [Fe, N, N, N, N, C, C, C, C, H]

		self.E = []		# total Energy
		self.Edisp = [] # E VDW for D3 correction
		self.positions = []
		self.cell = []
		self.vibrations = []

		self.total_charge = []
		self.magnetization = []
		self.IRdisplacement = []

		self.total_force = []

		CL = 0
		TC, MG, TF, IR = 0, 0, 0, 0 # mark for 
		#	TC : total charge
		#	MG : magnetization
		# 	TF : Total force
		# 	IR : Atoms displace in IR calculation

		for i, n in enumerate(f):
			vec = [m for m in n.split(' ') if m != '']

			# === POTCAR === #
			if 'POTCAR' in n: 
				try: 	
					if not n.split(':')[1] in self.POTCAR_full:
						self.POTCAR_full.append(n.split(':')[1])

					if not n.split(':')[1].strip().split(' ')[1] in self.POTCAR:
						self.POTCAR.append( n.split(':')[1].strip().split(' ')[1] )

				except: pass

			if 'POSCAR' in n: 
				try: self.POSCAR = n.split(':')[1]
				except: pass
			if 'E-fermi' in n : self.E_fermi = vec[2]
			if 'NIONS' in n: self.NIONS = vec[11]

			# === ION type === #
			if 'ions per type' in n:
				for ion_n, ion in enumerate(vec[4:]):
					self.atoms_names_list += [self.POTCAR[ion_n]]*int(ion) 	# [Fe, N, N, N, N, C, C, C, C, H]
				self.atoms_names_ID = self.POTCAR   	# [Fe, N, C, H] 


			# --- store CELL --- # dim = steps, 3, 3
			if 'direct lattice vectors' in n: CL = 4; cell = []
			if int(CL) > 0 :
				CL -= 1
				if int(CL) < 4:
					try: cell.append( [float(vec[0]),float(vec[1]),float(vec[2])] )
					except: pass
				if CL == 0: self.cell.append( np.array(cell))

			# --- store CHAGE --- # dim = steps, ions, 4
			if 'total charge ' in n: TC = int(self.NIONS)+4; total_charge = []
			if int(TC) > 0 : 
				TC -= 1
				if int(TC) < int(self.NIONS):
					try: total_charge.append( [float(vec[1]),float(vec[2]),float(vec[3]),float(vec[4])] )
					except: pass
				if TC == 0: self.total_charge.append( np.array(total_charge))

			# --- store CHAGE --- # dim = steps, ions, 4
			if 'magnetization (x)' in n: MG = int(self.NIONS)+4; magnetization = []
			if int(MG) > 0 :
				MG -= 1
				if int(MG) < int(self.NIONS):
					try: magnetization.append( [float(vec[1]),float(vec[2]),float(vec[3]),float(vec[4])] )
					except: pass
				if MG == 0: self.magnetization.append( np.array(magnetization))
			
			# --- store IR displacement --- # dim = steps, ions, 4
			if '2PiTHz' in n: IR = int(self.NIONS)+2; IRdisplacement = []; self.vibrations.append([float(vec[3]),float(vec[5]),float(vec[7]),float(vec[9])])
			if int(IR) > 0 :
				IR -= 1
				if int(IR) < int(self.NIONS):
					try: IRdisplacement.append( [	float(vec[0]),float(vec[1]),float(vec[2]),
													float(vec[3]),float(vec[4]),float(vec[5])] )
					except: pass
				if IR == 0: self.IRdisplacement.append( np.array(IRdisplacement))


			# --- store CHAGE --- # dim = steps, 1
			if 'Edisp' in n: self.Edisp.append(vec[-1][:-1])

			# --- store FORCE --- # dim = steps-1, ions, 6
			if 'TOTAL-FORCE' in n: TF = int(self.NIONS)+2; total_force = []

			if 'TOTAL-FORCE' in n: TF = int(self.NIONS)+2; total_force = []

			if int(TF) > 0 :
				TF -= 1
				if int(TF) < int(self.NIONS):
					try: total_force.append( [float(vec[0]),float(vec[1]),float(vec[2]),float(vec[3]),float(vec[4]),float(vec[5])] )
					except: pass
				if TF == 0: self.total_force.append( np.array(total_force))

			if 'energy  without entropy=' in n:
				self.E.append( float(vec[-1]) )

		self.total_charge 		= np.array(self.total_charge)
		self.magnetization 		= np.array(self.magnetization)
		self.total_force 		= np.array(self.total_force)
		self.IRdisplacement 	= np.array(self.IRdisplacement)
		self.vibrations		 	= np.array(self.vibrations)
		self.cell 				= np.array(self.cell)
		self.atoms_names_list 	= self.atoms_names_list
		self.Edisp 				= np.array(self.Edisp)

		#self.summary()

	def save_trajectory_step(self, atoms_position=None, atoms_names_list=None, atoms_names_ID=None, 
									file_name=None, append_data=True, traj_step=None, traj_time=None,save=True):

		atoms_names_ID = np.array(atoms_names_ID) if type(atoms_names_ID) != type(None) else np.array(self.atoms_names_ID)
		atoms_names_list = np.array(atoms_names_list) if type(atoms_names_list) != type(None) else np.array(self.atoms_names_list)
		traj_step = np.array(traj_step) if type(traj_step) != type(None) else np.array(0)
		traj_time = np.array(traj_time) if type(traj_time) != type(None) else np.array(0)

		path = file_name.split('/')[:-1]
		path = file_name.split('/')[:-1]
		name = file_name.split('/')[-1]

		if append_data: 	file = open('{}'.format(file_name), 'a' )
		else:  				file = open('{}'.format(file_name), 'rw')

		atom_number = atoms_position.shape[0]
		file.write('{}\n'.format(atom_number) )
		file.write(' TIme: {0} fs -- step count: {1}\n'.format(traj_time, traj_step) )
		for atom_n in range(atom_number):
			file.write('{} \t {:.7} \t {:.7} \t {} \t 0 \t 0\n'.format(atoms_names_list[atom_n], atoms_position[atom_n,0], atoms_position[atom_n,1], atoms_position[atom_n,2]  ) )
		
		file.close()
		return None

	def save_IR(self, modes='all', filtering=True):
		if modes == 'all': 
			modes = self.vibrations.shape[0]

		for i, n in enumerate(self.vibrations):
			file = open(f'mode_{i}.xsf', 'w' )
			file.write('CRYSTAL\n' )
			file.write('PRIMVEC\n' )
			file.write(f'  {self.cell[0][0][0]}   {self.cell[0][0][1]}   {self.cell[0][0][2]}\n' )
			file.write(f'  {self.cell[0][1][0]}   {self.cell[0][1][1]}   {self.cell[0][1][2]}\n' )
			file.write(f'  {self.cell[0][2][0]}   {self.cell[0][2][1]}   {self.cell[0][2][2]}\n' )
			file.write('PRIMCOORD\n' )
			file.write(f' {self.IRdisplacement.shape[1]}\n' )
			if filtering:
				maxs = np.max( np.linalg.norm(self.IRdisplacement[i, :, 3:], axis=1) )
				filt = maxs*0.3 > np.linalg.norm(self.IRdisplacement[i, :, 3:], axis=1)
				self.IRdisplacement[i, :, 3][filt] = 0
				self.IRdisplacement[i, :, 4][filt] = 0
				self.IRdisplacement[i, :, 5][filt] = 0
			
			for j in range(self.IRdisplacement.shape[1]):
				str_list = '   '.join([str(o) for o in self.IRdisplacement[i, j, :]])
				file.write(f'{self.atoms_names_list[j]}  {str_list}\n' )
			file.close()

	def integrity(self):
		try: 	self.POTCAR
		except: self.POTCAR = None

		try: 	self.Edisp
		except: self.Edisp = None

		try: 	self.cell
		except: self.cell = None

		try: 	self.total_force
		except: self.total_force = None

		try: 	self.total_charge
		except: self.total_charge = None

		try: 	self.magnetization
		except: self.magnetization = None

		return True

	def summary(self, ):
		OUTCAR_resumen = self.resume()

		toprint  =  '    |       *-->[OUTCAR]\n'
		toprint += f'    |           |\n'
		toprint += f'    |           *-> (POTCAR):{self.POTCAR}\n'
		toprint += f'    |           *-> (Edisp):{self.Edisp[-1]} ({np.mean(np.array(self.Edisp, dtype=np.float64))}+-{np.std(np.array(self.Edisp, dtype=np.float64))}) \n'
		#toprint +=  '    |           *-> (Cell):{0}'.format(';'.join(self.cell[-1])) + '}\n'
		#toprint += f'    |           *-> (Magnetization):{self.magnetization[-1]}\n'
		unique_atoms_names, count_atoms_names = np.unique(np.array(self.atoms_names_list), return_counts=True) 
		unique_filter = [ n in unique_atoms_names[count_atoms_names<=3 ] for n in self.atoms_names_list ]
		for i, n in enumerate(self.atoms_names_list): 
			if unique_filter[i]:
				toprint += f'    |           |     {n}:  Mag {self.magnetization[-1][i,:]}\n'
				toprint += f'    |           |          Chg {self.total_charge[-1][i,:]} \n'
		toprint += f'    |           *-> (Fmax):{np.max(self.total_force[-1][:,3:])}\n'
		print(toprint)
		
		return None

	def resume(self, ):
		self.integrity()

		info_OUTCAR = { 'POTCAR':self.POTCAR, 					'Edisp':self.Edisp[-1] if len(self.Edisp) > 0 else 0, 
						'cell':self.cell[-1], 					'total_force':self.total_force[-1], 	
						'total_charge':self.total_charge[-1], 	'magnetization':self.magnetization[-1],
						'atoms_names_list':self.atoms_names_list}

		return {'list':[[key, value] for i, (key, value) in enumerate( info_OUTCAR.items())], 'dict':info_OUTCAR}

	def plot_forces(self, where=None, figure=None):
		try:
			if figure == None:
				plt.hist(self.total_force[:,3], alpha=0.9, color=self.plot_color[0])
				plt.plot(self.total_force[:,3],-1+np.zeros(self.total_force[:,3].shape[0]), 'o' , color=self.plot_color[0])
				plt.hist(self.total_force[:,4], alpha=0.9, color=self.plot_color[1])
				plt.plot(self.total_force[:,4],-2+np.zeros(self.total_force[:,4].shape[0]), 'o' , color=self.plot_color[1])
				plt.hist(self.total_force[:,5], alpha=0.9, color=self.plot_color[2])	
				plt.plot(self.total_force[:,5],-3+np.zeros(self.total_force[:,5].shape[0]), 'o' , color=self.plot_color[2])
				plt.xlabel('Force mod')
		except Exception as e:
			print('*'*10 + str('error :: OUTCAR NOT loaded :: ') + '*'*10)
		else:
			pass
		finally:
			pass
	
	def get_magnetization(self, atoms_list=None, step=None, ):
		# return magnetization of atoms in 'atoms_list' from 'step' 
		step = -1 if step == None else step
		atoms_list = np.array(atoms_list) if type(atoms_list) == list else atoms_list

		try:
			if type(atoms_list) is np.ndarray:	return self.magnetization[step][atoms_list]
			else: 								return self.magnetization[step]
		except: return None

	def get_charge(self, atoms_list=None, step=None, ):
		# return total_charge of atoms in 'atoms_list' from 'step' 
		step = -1 if step == None else step
		atoms_list = np.array(atoms_list) if type(atoms_list) == list else atoms_list

		try:
			if type(atoms_list) is np.ndarray:	return self.total_charge[step][atoms_list]
			else: 								return self.total_charge[step]
		except: return None

	def get_dipolarmoment(self, atoms_list=None, step=None):
		pass

	def get_first_neighbour_distance_matrix(self, atoms_position=None, cell=None, save=True ):

		N = int(self.NIONS)
		cell = self.cell[0] if type(cell) == type(None) else cell
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
			self.vector_dist = vector_dist
			self.versor_dist = versor_dist

		return dist, vector_dist

	def get_embedding(self, distance_matrix=None, atoms_names_ID=None, atoms_names_list=None, embedding_parameters={'partition':200, 'min distance':0.5, 'max distance':3.0, 'range distance':0.05}, save=True):

		def gaussian(mu, sigma, n, ): return np.e**(-1.0/2 * ((mu-np.arange(n))/sigma)**2) * 0.3989422804/sigma

		# == define some variables == #
		Dmin, Dmax ,Dsig, N = embedding_parameters['min distance'], embedding_parameters['max distance'], embedding_parameters['range distance'], embedding_parameters['partition']
		Nsig = 200 * Dsig * N / (Dmax-Dmin)

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


	def get_topologic_embedding(self, distance_matrix=None, vector_dist=None, atoms_names_ID=None, atoms_names_list=None, embedding_parameters={'partition':200, 'min distance':0.5, 'max distance':3.0, 'range distance':0.05}, save=True):

		# == define some variables == #
		Dmin, Dmax ,Dsig, embedding_size = embedding_parameters['min distance'], embedding_parameters['max distance'], embedding_parameters['range distance'], embedding_parameters['partition']

		atoms_names_ID = np.array(atoms_names_ID) if type(atoms_names_ID) != type(None) else np.array(self.atoms_names_ID)
		atoms_names_list = np.array(atoms_names_list) if type(atoms_names_list) != type(None) else np.array(self.atoms_names_list)
		distance_matrix = np.array(distance_matrix) if type(distance_matrix) != type(None) else np.array(self.distance_matrix)
		vector_dist = np.array(vector_dist) if type(vector_dist) != type(None) else np.array(self.vector_dist) 

		# == asing memory space == #
		distance_matrix[distance_matrix==0] = np.inf
		atomic_embedding = np.zeros( (atoms_names_list.shape[0], atoms_names_ID.shape[0], 3, embedding_size) )

		# == Filter by distance == #
		for atom_n, atom in enumerate(atoms_names_list):
			for atom_id_n, atom_id in enumerate(atoms_names_ID):
				atomic_embedding[atom_n, atom_id_n, :] = vector_dist[:, atom_n, atoms_names_list==atom_id][:,distance_matrix[atom_n, atoms_names_list==atom_id].argsort()][:, :embedding_size]
		
		if save: self.atomic_embedding = atomic_embedding
		return atomic_embedding # [ atoms , ID, 3, embeding_size ]

	def get_all_embedding(self, embedding_parameters, ):
		# === calls the embedding getter === #
		return self.get_embeddings()

	def get_embeddings(self, atoms_names_ID=None, atoms_names_list=None, dimention='topologic',
			embedding_parameters={'partition':40, 'min distance':0.5, 'max distance':3.0, 'range distance':0.05}, save=True):

		# === Set input parameters === #
		atoms_names_ID = np.array(atoms_names_ID) if type(atoms_names_ID) != type(None) else np.array(self.atoms_names_ID)
		atoms_names_list = np.array(atoms_names_list) if type(atoms_names_list) != type(None) else np.array(self.atoms_names_list)

		# === alocate mem === #
		embedding_size = len(self.POTCAR)*embedding_parameters['partition'] if dimention == '2D' else len(self.POTCAR)*embedding_parameters['partition']*3
		embedding_vector = np.zeros((self.total_force.shape[0], self.total_force.shape[1], embedding_size))
		force_vector = np.zeros((self.total_force.shape[0], self.total_force.shape[1], 3))
		partition = embedding_parameters['partition']

		# === Get the embedding for each STEP === #
		for step, postion_forces in enumerate(self.total_force):
			if self.total_force.shape[0] > 0:
				distance_matrix, vector_dist = self.get_first_neighbour_distance_matrix(atoms_position=postion_forces[:, :3], save=True )
				embedding = self.get_embedding(		distance_matrix=distance_matrix, atoms_names_ID=atoms_names_ID, atoms_names_list=atoms_names_list, 
													embedding_parameters=embedding_parameters) if dimention == '2D' else self.get_topologic_embedding(
													distance_matrix=distance_matrix, vector_dist=vector_dist,
													atoms_names_ID=atoms_names_ID, atoms_names_list=atoms_names_list, 
													embedding_parameters=embedding_parameters)
				force_vector[step, :] = postion_forces[:, 3:]  
				embedding_vector[step, :] = np.reshape(embedding, [self.total_force.shape[1], embedding_size] , order='C') 
		
		# === Store X and Y for neural regresor === #
		X = { atomID:embedding_vector[:, atoms_names_list==atomID, :].reshape(-1,embedding_size) for atomID in atoms_names_ID }
		Y = { atomID:force_vector[:, atoms_names_list==atomID].reshape(-1,3) for atomID in atoms_names_ID }

		# === Store data === #
		if save: 
			self.embedding_vector = embedding_vector
			self.force_vector = force_vector
			self.ff_X = X
			self.ff_Y = Y

		return embedding_vector, force_vector, X, Y

	def get_NNforces(self, position=None, model=None, atoms_names_ID=None, atoms_names_list=None, dimention='topologic',
		embedding_parameters={'partition':200, 'min distance':0.5, 'max distance':3.0, 'range distance':0.05}, save=True):
		
		atoms_names_ID = np.array(atoms_names_ID) if type(atoms_names_ID) != type(None) else np.array(self.atoms_names_ID)
		atoms_names_list = np.array(atoms_names_list) if type(atoms_names_list) != type(None) else np.array(self.atoms_names_list)
		partition = embedding_parameters['partition']
		embedding_size = len(self.POTCAR)*embedding_parameters['partition'] if dimention == '2D' else len(self.POTCAR)*embedding_parameters['partition']*3
		
		distance_matrix, vector_dist = self.get_first_neighbour_distance_matrix(atoms_position=position, save=True )

		embedding = self.get_embedding(		distance_matrix=distance_matrix, atoms_names_ID=atoms_names_ID, atoms_names_list=atoms_names_list, 
													embedding_parameters=embedding_parameters) if dimention == '2D' else self.get_topologic_embedding(
													distance_matrix=distance_matrix, vector_dist=vector_dist,
													atoms_names_ID=atoms_names_ID, atoms_names_list=atoms_names_list, 
													embedding_parameters=embedding_parameters)
		forces = np.zeros_like(position)
		for atom_n, atom in enumerate(atoms_names_list):
			#if atom =='Fe':
			#	# borrar !
			#	plt.plot( np.reshape(embedding[atom_n,:,:], [embedding.shape[1]*partition] , order='C') )
			#	plt.show()
			forces[atom_n, :] = model[atom].predict( [embedding[atom_n,:,:].flatten(order='C')] )

		return forces 

'''
# === OUTCAR IR === #
OC = OUTCAR()
OC.load('/home/akaris/Documents/code/VASP/v4.1/files/OUTCAR/FeTPP/IR/rmFeNC/OUTCAR')
OC.summary()
OC.save_IR()

print(	OC.total_charge.shape , OC.IRdisplacement.shape,
	OC.vibrations)
print( OC.total_charge[0,:,-1] - OC.total_charge[80,:,-1] )
def gaussian(n, mu, s):
	return np.e**( -(mu-n)**2/s )

n = np.linspace(1,2000,2000)
s = 50
data = np.zeros(2000)
for m in OC.vibrations:
	data += gaussian(n, m[2], s)
plt.plot( data )
plt.show()
# === OUTCAR IR === #

OC.get_embeddings(embedding_parameters={'partition':5, 'min distance':0.5, 'max distance':7.0, 'range distance':0.001})
'''

'''
import scipy.io
from sklearn.neural_network import MLPRegressor
import pickle

make_embedding = False
if make_embedding:
	# ========== GET embedding ========== #
	OC.get_embeddings(embedding_parameters={'partition':400, 'min distance':0.5, 'max distance':7.0, 'range distance':0.001})

	# ========== LOAD data ========== #
	#** Ru **# #** Ru **# #** Ru **# #** Ru **# #** Ru **# #** Ru **# #** Ru **# #** Ru **# #** Ru **# #** Ru **# 
	X_Ru_train, y_Ru_train = OC.X['Ru'], OC.Y['Ru']
	np.savetxt('/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/Ru/X_train', X_Ru_train)
	np.savetxt('/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/Ru/y_train', y_Ru_train)
	# ========== Train model NN ========== #
	regrRu = MLPRegressor(random_state=3, activation='relu', #early_stopping=True, validation_fraction=0.10,
						hidden_layer_sizes=(200, 100, 100, 100, 100, 100, 200, ),max_iter=1000, #learning_rate='adaptive', warm_start=True,
						#hidden_layer_sizes=(  400, 540, 540, 400,  ),max_iter=5000,
						momentum=0.9, verbose=True, n_iter_no_change=100, ).fit(X_Ru_train, y_Ru_train)

	# ========== SAVE model NN ========== #
	pkl_filename = "/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/Ru/NNmodel_Ru.pkl"
	with open(pkl_filename, 'wb') as file:
		pickle.dump(regrRu, file)

	#** N **# #** N **# #** N **# #** N **# #** N **# #** N **# #** N **# #** N **# #** N **# #** N **# #** N **# 
	X_N_train, y_N_train = OC.X['N'], OC.Y['N']
	np.savetxt('/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/N/X_train', X_N_train)
	np.savetxt('/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/N/y_train', y_N_train)
	# ========== Train model NN ========== #
	regrN = MLPRegressor(random_state=3, activation='relu', #early_stopping=True, validation_fraction=0.10,
						hidden_layer_sizes=(200, 100, 100, 100, 100, 100, 200, ),max_iter=1000, #learning_rate='adaptive', warm_start=True,
						#hidden_layer_sizes=(  400, 540, 540, 400,  ),max_iter=5000,
						momentum=0.9, verbose=True, n_iter_no_change=100, ).fit(X_N_train, y_N_train)

	# ========== SAVE model NN ========== #
	pkl_filename = "/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/N/NNmodel_N.pkl"
	with open(pkl_filename, 'wb') as file:
	    pickle.dump(regrN, file)

	#** C **# #** C **# #** C **# #** C **# #** C **# #** C **# #** C **# #** C **# #** C **# #** C **# #** C **# 
	X_C_train, y_C_train = OC.X['C'], OC.Y['C']
	np.savetxt('/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/C/X_train', X_C_train)
	np.savetxt('/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/C/y_train', y_C_train)
	# ========== Train model NN ========== #
	regrC = MLPRegressor(random_state=3, activation='relu', #early_stopping=True, validation_fraction=0.10,
						hidden_layer_sizes=(200, 100, 100, 100, 100, 100, 200, ),max_iter=1000, #learning_rate='adaptive', warm_start=True,
						#hidden_layer_sizes=(  400, 540, 540, 400,  ),max_iter=5000,
						momentum=0.9, verbose=True, n_iter_no_change=100, ).fit(X_C_train, y_C_train)

	# ========== SAVE model NN ========== #
	pkl_filename = "/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/C/NNmodel_C.pkl"
	with open(pkl_filename, 'wb') as file:
	    pickle.dump(regrC, file)

	#** H **# #** H **# #** H **# #** H **# #** H **# #** H **# #** H **# #** H **# #** H **# #** H **# #** H **# 
	X_H_train, y_H_train = OC.X['H'], OC.Y['H']
	np.savetxt('/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/H/X_train', X_H_train)
	np.savetxt('/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/H/y_train', y_H_train)
	# ========== Train model NN ========== #
	regrH = MLPRegressor(random_state=3, activation='relu', #early_stopping=True, validation_fraction=0.10,
						hidden_layer_sizes=(200, 100, 100, 100, 100, 100, 200, ),max_iter=1000, #learning_rate='adaptive', warm_start=True,
						#hidden_layer_sizes=(  400, 540, 540, 400,  ),max_iter=5000,
						momentum=0.9, verbose=True, n_iter_no_change=100, ).fit(X_H_train, y_H_train)

	# ========== SAVE model NN ========== #
	pkl_filename = "/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/H/NNmodel_H.pkl"
	with open(pkl_filename, 'wb') as file:
	    pickle.dump(regrH, file)

	#** Au **# #** Au **# #** Au **# #** Au **# #** Au **# #** Au **# #** Au **# #** Au **# #** Au **# #** Au **# 
	X_Au_train, y_Au_train = OC.X['Au'], OC.Y['Au']
	np.savetxt('/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/Au/X_train', X_Au_train)
	np.savetxt('/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/Au/y_train', y_Au_train)
	# ========== Train model NN ========== #
	regrAu = MLPRegressor(random_state=3, activation='relu', #early_stopping=True, validation_fraction=0.10,
						hidden_layer_sizes=(200, 100, 100, 100, 100, 100, 200, ),max_iter=1000, #learning_rate='adaptive', warm_start=True,
						#hidden_layer_sizes=(  400, 540, 540, 400,  ),max_iter=5000,
						momentum=0.9, verbose=True, n_iter_no_change=100, ).fit(X_Au_train, y_Au_train)

	# ========== SAVE model NN ========== #
	pkl_filename = "/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/Au/NNmodel_Au.pkl"
	with open(pkl_filename, 'wb') as file:
	    pickle.dump(regrAu, file)

load_embedding = True
if load_embedding:
	model = {}
	pkl_filename = "/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/Ru/NNmodel_Ru.pkl"
	with open(pkl_filename, 'rb') as file:
		model['Ru'] = pickle.load(file)

	pkl_filename = "/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/N/NNmodel_N.pkl"
	with open(pkl_filename, 'rb') as file:
		model['N'] = pickle.load(file)

	pkl_filename = "/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/C/NNmodel_C.pkl"
	with open(pkl_filename, 'rb') as file:
		model['C'] = pickle.load(file)

	pkl_filename = "/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/H/NNmodel_H.pkl"
	with open(pkl_filename, 'rb') as file:
		model['H'] = pickle.load(file)

	pkl_filename = "/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/Au/NNmodel_Au.pkl"
	with open(pkl_filename, 'rb') as file:
		model['Au'] = pickle.load(file)
'''

'''
X_C_train = np.loadtxt('/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/C/X_train')
y_C_train = np.loadtxt('/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/C/y_train')
forces = model['C'].predict( X_C_train )
plt.plot(forces, y_C_train,'o')
plt.show()
error
'''

'''
# run dinamic #
path = '/home/akaris/Documents/code/VASP/v3.5/files/dataset/force_trainning/path4'
#OC.seve_trajectory_step(atoms_position=position, file_name=path, append_data=True, traj_step=None, traj_time=None,save=True)
position = OC.total_force[5,:,:3]
velocidad = np.zeros_like(position)
force = np.zeros_like(position)
dt = 0.1
for n in range(200):
	forces = OC.get_NNforces(position=position, model=model, embedding_parameters={'partition':400, 'min distance':0.5, 'max distance':7.0, 'range distance':0.001}, save=True)
	velocidad[:4,:] += forces[:4,:]*dt 
	position  += velocidad*dt 
	OC.save_trajectory_step(atoms_position=position, file_name=path, append_data=True, traj_step=None, traj_time=None,save=True)
	print(n)
	print(forces[0,:])
'''


'''	
print(OC.NIONS ) 
print(OC.atoms_names_list)
print(OC.POTCAR )
'''

