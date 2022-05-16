print('Python 3.X - https://www.python.org/download/releases/3.0/')
###################################
# Step 1 || Load python libraries #
###################################
# *** warning supresion
import warnings, time
warnings.filterwarnings("ignore")

# *** numeric libraries *** #
try:
	import numpy as np
except: print('ERROR :: DATA.import_libraries() :: can not import numpy ')

try:
	import matplotlib.pyplot as plt
	import matplotlib.axes as ax
	import matplotlib.patches as patches
except:	print('ERROR :: System.import_libraries() :: can not import matplotlib ')

from scipy.signal import savgol_filter

try:
	from sklearn.decomposition import PCA
	from sklearn.preprocessing import StandardScaler
	from sklearn import linear_model
	from sklearn.model_selection import cross_val_predict
	from sklearn.metrics import mean_squared_error, r2_score
except:	print('ERROR :: DATA.import_libraries() :: can not import sklearn')

try:
	from ase import Atoms
	from ase.visualize import view
except:	print('ERROR :: DATA.import_libraries() :: can not import ase ')


try:
	from os import path
	import itertools, operator, logging, time, copy, pickle, os.path
except:  print('ERROR :: DATA.import_libraries() :: can not import itertools, operator, logging, time, copy, pickle or os')

# *** load own libraries *** #
try:	from src import POSCAR
except:	
	try: import POSCAR as POSCAR
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

try:	from src import DOSCAR
except:	
	try: import DOSCAR as DOSCAR
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')
try:	from src import OUTCAR
except:	
	try: import OUTCAR as OUTCAR
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

try:	from src import OSZICAR
except:	
	try: import OSZICAR as OSZICAR
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

try:	from src import ORR
except:	
	try: import ORR as ORR
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

class System(object):
	def __init__(self, 
					name=None	, description=None	,
					POSCAR=None	, INCAR=None		, CONTCAR=None	, EIGENVAL=None,
					DOSCAR=None	, PROCAR=None		, OSZICAR=None 	, PARCHGCAR=None,
					CHGCAR=None	, OUTCAR=None		, POTCAR=None,
					ZPE=None	, S=None):
		
		self.name = name
		self.description = description

		self.POSCAR = POSCAR
		self.INCAR = INCAR
		self.POTCAR = POTCAR
		self.CONTCAR = CONTCAR
		self.CHGCAR = CHGCAR
		self.DOSCAR = DOSCAR
		self.PROCAR = PROCAR
		self.EIGENVAL = EIGENVAL
		self.PARCHGCAR = PARCHGCAR
		self.OSZICAR =  OSZICAR
		self.OUTCAR = OUTCAR

		self.ZPE    = ZPE   # Zero point energy #
		self.S      = S     # Entropy #

		self.white_list = [
							'POSCAR', 'INCAR', 'CONTCAR', 'CHGCAR', 'DOSCAR', 'POTCAR',
							'PROCAR', 'EIGENVAL', 'PARCHGCAR', 'OSZICAR', 'OUTCAR']

		self.loaded = []
		self.loaded_dict = {
			'POSCAR':self.POSCAR,	'CONTCAR':self.CONTCAR,
			'INCAR':self.INCAR,		'OSZICAR':self.OSZICAR,
			'OUTCAR':self.OUTCAR,	'EIGENVAL':self.EIGENVAL,
			'DOSCAR':self.DOSCAR,	'PROCAR':self.PROCAR,
		}

	def load(self, path=None, file='POSCAR', v=True):
		# LOAD {file} from {path}
		if path == '' or path == None:  
			print('ERROR :: System.load() :: can not LOAD path: {} '.format(path) )

		if not type(file) == str or not file in self.white_list:    
			print('WARNNING :: System.load() :: can not LOAD file: {} '.format(file) )
		
		def load_POSCAR(path, v=True):
			try:
				poscar = POSCAR.POSCAR()
				poscar.load(file_name = path, v=v)
				self.POSCAR = poscar
				self.loaded.append('POSCAR')
				if v>0: self.POSCAR.summary()

			except: 
				if v==2: print('ERROR :: System.load() :: can not LOAD POSCAR')
		
		def load_POTCAR(path, v=True):
			try:
				potcar = POTCAR.POTCAR()
				potcar.load(file_name = path)
				self.POTCAR = potcar
				self.loaded.append('POSCAR')
			except: 
				if v==2: print('ERROR :: System.load() :: can not LOAD POTCAR')
		
		def load_CONTCAR(path, v=True):
			try:
				contcar = POSCAR.POSCAR()
				contcar.load(file_name = path)
				self.CONTCAR = contcar
				self.loaded.append('CONTCAR')
				if v>0: self.CONTCAR.summary()
			except: 
				if v==2: print('ERROR :: System.load() :: can not LOAD CONTCAR')

		def load_OUTCAR(path, v=True):
			try:
				outcar = OUTCAR.OUTCAR()
				outcar.load(file_name = path)
				self.OUTCAR = outcar 
				self.loaded.append('OUTCAR')
			except: 
				if v==2: print('ERROR :: System.load() :: can not LOAD OUTCAR')

		def load_OSZICAR(path, v=True):
			try:
				oszicar = OSZICAR.OSZICAR()
				oszicar.load(file_name = path)
				self.OSZICAR = oszicar
				self.loaded.append('OSZICAR')
				if v>0: self.OSZICAR.summary()
			except:
				if v==2: print('ERROR :: System.load() :: can not LOAD OSZICAR')

		def load_DOSCAR(path, v=True):
			try:
				doscar = DOSCAR.DOSCAR()
				#doscar.load(file_name = path) !!!!!!!!!!!!!!!!!!!!!!!!
				self.DOSCAR = doscar
				self.loaded.append('DOSCAR')
			except:
				if v==2: print('ERROR :: System.load() :: can not LOAD DOSCAR')

		if      file == 'POSCAR':   load_POSCAR(path,  v=v)
		if      file == 'CONTCAR':  load_CONTCAR(path, v=v)
		elif    file == 'OUTCAR':   load_OUTCAR(path,  v=v)
		elif    file == 'OSZICAR':  load_OSZICAR(path, v=v)
		elif    file == 'DOSCAR':   load_DOSCAR(path,  v=v)
		elif    file == 'POTCAR':   load_POTCAR(path,  v=v)

	def load_all(self, path, v=True):
		for file_load in self.white_list:		
			self.load(path='{}/{}'.format(path, file_load), file=file_load, v=v)

	def get_embeddings(self, path=None, embedding_parameters=None, save=False, v=True):
		# === Get data from OUTCAR === #
		embedding_parameters = {'partition':5, 'min distance':0.5, 'max distance':7.0, 'range distance':0.001} if type(embedding_parameters) == type(None) else embedding_parameters

		if v:
			t1 = time.time() 
			print('Getting embedding :: {}'.format(embedding_parameters))
			print(f'Processing :: { np.array(self.OUTCAR.E).shape[0] } steps' )

		embedding_vector, force_vector, ff_X, ff_Y = self.OUTCAR.get_embeddings(embedding_parameters=embedding_parameters)

		if v:
			resume_str = f'Time {time.time()-t1} \n'
			for i, (key, value) in enumerate( ff_X.items()):
				resume_str += f'\t\t({i}) Ion:{key} Samples:{value.shape[0]} Dimention:{value.shape[1]} \n'
			print(resume_str)

		# === Store data === #
		if save: 
			self.embedding_vector = embedding_vector
			self.force_vector = force_vector
			self.ff_X = ff_X
			self.ff_Y = ff_Y
			self.embedding_dict = {'embedding_vector':embedding_vector, 'force_vector':force_vector, 
				'ff_X':ff_X, 'ff_Y':ff_Y}
				
		return {'embedding_vector':embedding_vector, 'force_vector':force_vector, 
				'ff_X':ff_X, 'ff_Y':ff_Y}

	def summary(self, ):
		#toprint = 'Loaded files : '
		#for ld in self.loaded:	toprint += '\t{};'.format(ld)
		#print(toprint)
		
		if 'CONTCAR' in self.loaded and self.CONTCAR == object:	self.CONTCAR.summary() 
		if 'POTCAR' in self.loaded 	and self.POTCAR == object:	self.POTCAR.summary() 
		if 'INCAR' in self.loaded 	and self.INCAR == object:	self.INCAR.summary() 
		if 'OSZICAR' in self.loaded and self.OSZICAR == object:	self.OSZICAR.summary() 
		if 'OUTCAR' in self.loaded 	and self.OUTCAR == object:	self.OUTCAR.summary() 

		return	None





