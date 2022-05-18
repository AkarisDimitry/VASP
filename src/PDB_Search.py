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

class PDB_Search(object): # generador de datos
	def __init__(self, path=None, name=None):
		self.path = path 
		self.name = name 

		self.result = None

	def timer(func):
		def wrapper(*args, **kwargs):
			before = time.time()
			r = func(*args, **kwargs)
			name = str(func).split(' ')[1]

			print( f'{name} : {len(r) if type(r) == list else r} (execution time {time.time()-before}s)' ) 
			return r
		return wrapper

	@timer
	def search(self, query, save=True, v=True):
		result = Query(query).search()
		
		if v:	print(f'About {len(result)} results for {query}')
		if save: 
			self.query = query
			self.result = result
		
		return result

	@timer
	def download(self, search=None, path=None, v=True):
		search = search if not search is None else self.result
		path = path if not path is None else self.path

		if not os.path.isdir(f'download'):  os.makedirs(f'download')
		if not os.path.isdir(f'download/{self.query}'):  os.makedirs(f'download/{self.query}')

		for name in search:
			if not os.path.isfile(f'download/{self.query}/{name}.pdb'):
				if v: print(f'Downloading {self.query} {name}...')
				urllib.request.urlretrieve(f'http://files.rcsb.org/download/{name}.pdb', f'download/{self.query}/{name}.pdb')
			else: 
				if v: print(f'Allready dowload {self.query}-{name}...')

		return True

PDBs = PDB_Search()
S1 = PDBs.search(query="NDM")
S2 = PDBs.search(query='New Delhi "metallo-beta-lactamase"')
S3 = PDBs.search(query='INHIBITOR')

for i, n in enumerate( set(S1).intersection(set(S2)).intersection(set(S3)) ):
	print(i, n)

#PDBs.download()

