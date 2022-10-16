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

class CHGCAR(object):
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

		self.CHARGE = None
		self.CHARGE_size = None

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
	def load_poscar(self, ):
		pass

	@Logs.LogDecorator()
	def load_charge(self, file_name=None, save=True):
		file_name = file_name if not file_name is None else self.file_name

		check = [0, 0, 0, 0]
		read = [0, 0]
		augmentation_occupancies = []

		charge_size 		= None
		charge 				= None
		magnetization_size 	= None
		magnetization		= None

		with open(file_name) as file:
			for i, n in enumerate(file):
				vec = [m for m in n.replace('\t',' ').split(' ') if m != '' and m != '\n']

				# prepare for read structure

				# prepare for read CHG  total charge density (spin up + spin down) 
				if check[0] == 0 and i > 5 and len(vec) == 0: check = [1, i, 0]
					
				elif check[0] == 1 and int(check[2]/5) + 1 + check[1] == i:
					# Reading total charge density (spin up + spin down) from CHGCAR
					print( f'>> Reading total charge density [{int(vec[0])} {int(vec[1])} {int(vec[2])}] (spin up + spin down) from CHGCAR')
					charge_size = [ int(vec[0]), int(vec[1]), int(vec[2]) ]
					charge = np.zeros((int(vec[0]), int(vec[1]), int(vec[2])))
					# read augmentation after len(augmentation_occupancies)%5 lines 
					check, read = [2, i, int(vec[0]) * int(vec[1]) * int(vec[2])], [1, i+1]

				elif check[0] == -1 and int(check[2]/5) + 1 + check[1] == i:
					# Reading magnetization density (spin up - spin down) from CHGCAR
					print( f'>> Reading magnetization density [{int(vec[0])} {int(vec[1])} {int(vec[2])}] (spin up - spin down) from CHGCAR')
					magnetization_size = [ int(vec[0]), int(vec[1]), int(vec[2]) ]
					magnetization = np.zeros((int(vec[0]), int(vec[1]), int(vec[2])))
					# read augmentation after len(augmentation_occupancies)%5 lines 
					check, read = [2, i, int(vec[0]) * int(vec[1]) * int(vec[2])], [-1, i+1]
					print(f'>> Reading augmentation occupancies... ')

				elif check[0] == 2 and int(check[2]/5) + 1 + check[1] == i:

					if vec[0] == 'augmentation':
						# Reading CHG  total charge density (spin up + spin down) 
						augmentation_occupancies.append([])
						# read augmentation after len(augmentation_occupancies)%5 lines 
						check, read = [2, i if int(vec[3])%5 == 0 else i+1, int(vec[3])], [2, i]

					else:	
						# Reading CHG  total charge density (spin up + spin down) 
						print(f'>> Reading augmentation occupancies... ')
						# go to Reading magnetization density after len(augmentation_occupancies)%5 lines 
						check, read = [-1, i-1 if len(augmentation_occupancies)%5 == 0 else i, len(augmentation_occupancies)], [2, i]

				if read[0] == 1:
					linear_index = i-read[1]
					if linear_index >= 0:
						try:
							index = [[	(linear_index*5+n)%charge_size[0],
										int((linear_index*5+n)/charge_size[0])%charge_size[1],
										int((linear_index*5+n)/(charge_size[0]*charge_size[1]))%charge_size[2]
									] for n in range(5)]
							for j, v in enumerate(vec):	charge[index[j][0], index[j][1], index[j][2]] = float(v)
						except:	print(f'ERROR :: CHGCAR.load_charge() :: error while reading line {i} \n >> {vec}')

				elif read[0] == -1:
					linear_index = i-read[1]
					if linear_index >= 0:
						try:
							index = [[	(linear_index*5+n)%magnetization_size[0],
										int((linear_index*5+n)/magnetization_size[0])%magnetization_size[1],
										int((linear_index*5+n)/(magnetization_size[0]*magnetization_size[1]))%magnetization_size[2]
									] for n in range(5)]
							for j, v in enumerate(vec):	magnetization[index[j][0], index[j][1], index[j][2]] = float(v)
						except:	print(f'ERROR :: CHGCAR.load_charge() :: error while reading line {i} \n >> {vec}')

		print(f'Reading compleat {file_name}')
		if save:
			self.magnetization_size 		= magnetization_size
			self.magnetization 				= magnetization
			self.charge_size 				= charge_size
			self.charge 					= charge

			self.augmentation_occupancies 	= augmentation_occupancies

	@Logs.LogDecorator()
	def plot(self, plotter='mlab'):

		if plotter == 'mlab':
			from mayavi import mlab
			values = self.charge

			mlab.contour3d(values) #, contours=[+67.924])
			#mlab.savefig('surface.obj')

		elif plotter == 'plotly':
			import plotly.graph_objects as go
			X, Y, Z = np.mgrid[-5:5:self.charge_size[0]*1j, -5:5:self.charge_size[1]*1j, -5:5:self.charge_size[2]*1j]
			values = self.charge

			fig = go.Figure(data=go.Isosurface(
			    x=X.flatten(),
			    y=Y.flatten(),
			    z=Z.flatten(),
			    value=values.flatten(),
			    isomin=10,
			    isomax=40,
			    caps=dict(x_show=False, y_show=False)
			    ))
			fig.show()
if __name__ == "__main__":
	# === Charge diference plot === #
	chg = CHGCAR()
	chg.load_charge(		file_name='/home/akaris/Documents/code/VASP/v4.6/files/CHGCAR/CoFeTPyPCo-CoFeTPyP[Co]/CHGCAR')

	chg_missCo = CHGCAR()
	chg_missCo.load_charge(	file_name='/home/akaris/Documents/code/VASP/v4.6/files/CHGCAR/CoFeTPyPCo-CoFeTPyP[Co]/CHGCAR_missFe')

	chg_Co = CHGCAR()
	chg_Co.load_charge(		file_name='/home/akaris/Documents/code/VASP/v4.6/files/CHGCAR/CoFeTPyPCo-CoFeTPyP[Co]/CHGCAR_Fe')

	chg2 = CHGCAR()
	chg2.charge = chg.charge - chg_missCo.charge - chg_Co.charge
	
	print('max-min', np.max(chg2.charge), np.min(chg2.charge))
	chg2.plot()
	plt.show()

	asdfas
	# === Charge diference plot === #
	chg = CHGCAR()
	chg.load_charge(		file_name='/home/akaris/Documents/code/VASP/v4.6/files/CHGCAR/CHGCAR')

	chg_missCo = CHGCAR()
	chg_missCo.load_charge(	file_name='/home/akaris/Documents/code/VASP/v4.6/files/CHGCAR/CHGCAR_missCo')

	chg_Co = CHGCAR()
	chg_Co.load_charge(		file_name='/home/akaris/Documents/code/VASP/v4.6/files/CHGCAR/CHGCAR_Co')

	chg2 = CHGCAR()
	chg2.charge = chg.charge - chg_missCo.charge - chg_Co.charge
	
	print('max-min', np.max(chg2.charge), np.min(chg2.charge))
	chg2.plot()
	plt.show()



