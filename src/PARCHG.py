import numpy as np 	
import matplotlib.pyplot as plt
try:	from ase.visualize import view
except: print('can not load ASE module. (pipX install ase-atomistics)')
import plotly.graph_objects as go

class PARCHG(object):
	def __init__(self, name=None):
		self.name = name
		self.path = None
		self.file_name = None

		self.atoms = None # np.array(3,N)
		self.atoms_list = None # [np.array(3,n(Fe)), np.array(3,n(N)), np.array(3,n(C)), np.array(3,n(H))]
		self.n = None # N
		self.atoms_names_list = None # kpoint
		self.atoms_names = None # [Fe, N, C, H]
		self.atoms_number = None # [n(Fe), n(N), n(C), n(H) ]

		self.cell = None 

		self.CHG		= None
		self.CHG_grid	= None
		self.MAG 		= None	
		self.MAG_grid	= None	

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

	def load(self, file_name=None):
		if file_name != None: self.file_name = file_name

		# ------------------------ LOAD data from POSCAR ------------------------ #
		# name              :   STR     :   path+name of file to load 
		# ------------------------------- #     # ------------------------------- #
		# atom_name         :   LIST    :   List with all atoms names with-out repeat
		# atom_numbers      :   LIST    :   List with the number of atoms of each specie
		# atom_position     :   LIST    :   contains XYZ data from all atoms divided by species
		# atom_allposition  :   N-MAT   :   Numpy array with all XYZ data
		# cell              :   N-MAT   :   Numpy array with cell parameters
		# N                 :   INT     :   Integer with total number of atoms
	    # ------------------------------- #    # ------------------------------- #
		file = open(self.file_name, 'r')
		cell = np.zeros((3,3)) ; atom_name = [] ; atom_numbers = [] ; N = 0 ; atom_position = [] ; atom_allposition = []; direct = 0

		UP = 0
		DOWN = 0
		UP_data = []
		for i, n in enumerate(file):
			vec = [m for m in n.split(' ') if m != '' and m != '\n']
			
			if len(vec) > 0 and vec[-1][-1:] == '\n': vec[-1] = vec[-1][:-1]

			print(i, vec)

			if i == 0: 			self.name = n

			for m in range(3):
				if i == 2+m:  cell[m,0]=float(vec[0]);  cell[m,1]=float(vec[1]); cell[m,2]=float(vec[2])
			if i == 5:
				for m in vec: atom_name.append(m)
			if i == 6:
				for m in vec:  N += int(m); atom_numbers.append(int(m)); atom_position.append( np.zeros((int(m), 3)) )
				atom_allposition = np.zeros((N, 3))

			if i == 8 and vec[0]=='Direct': direct = 1

			if UP == 1: 		self.grid = np.array([int(vec[0]), int(vec[1]), int(vec[2])]); self.end_poscar = int(i); break

			if i > 5 and vec == []: UP = 1

			if UP == 0 and DOWN == 0 and i > 8 :
				var = 9
				for j, m in enumerate(atom_numbers):
					if i >= var and i < var+m:  
						atom_position[j][i-var,0]=float(vec[0]) ;  atom_position[j][i-var,1]=float(vec [1]) ;   atom_position[j][i-var,2]=float(vec[2])
						atom_allposition[i-9,0]=float(vec[0])   ;  atom_allposition[i-9,1]=float(vec [1])   ;   atom_allposition[i-9,2]=float(vec[2])
					var += m

			if i > 100: break

		if direct==1: atom_allposition = np.dot(atom_allposition, cell)

		self.cell = np.array(cell)
		self.n = N
		self.atoms = np.array(atom_allposition)
		self.atoms_list = atom_position
		self.atoms_names = atom_name
		self.atoms_names_list = [n for n, m in zip(atom_name, atom_numbers) for m in range(m)] 

		grid_len = int(self.grid[0]*self.grid[1]*self.grid[2]/10)
		print( grid_len)

		try:
			self.CHG = np.loadtxt(fname=self.file_name, delimiter=None, skiprows=self.end_poscar+1, max_rows=grid_len).T
			print('CHARGE grid readed, grid size : {} {},'.format(self.CHG.shape[0], self.CHG.shape[1]))
		except:
			print('ERROR :: PARCHG.load() :: can NOT read CHARGE grid from {}'.format(file_name) )

		self.CHG_grid = np.reshape(self.CHG, self.grid, order='F')

		try:
			self.MAG = np.loadtxt(fname=self.file_name, delimiter=None, skiprows=grid_len+self.end_poscar+3).T
			print('MAGNETIZATION grid readed, grid size : {} {},'.format(self.MAG.shape[0], self.MAG.shape[1] ))
		except:
			print('ERROR :: PARCHG.load() :: can NOT read MAGNETIZATION grid from {}'.format(file_name) )

		self.MAG_grid = np.reshape(self.MAG, self.grid, order='F')

		print(self.MAG_grid.shape)

		return atom_name, atom_numbers, atom_position, np.array(atom_allposition), np.array(cell), N


	def direct_copy(self, file_name=None):
		if file_name != None: self.file_name = file_name
		file_up 	= open('{}_CHARGE.vasp'.format(self.file_name), 'w')
		file_down 	= open('{}_MAGENTIZATION.vasp'.format(self.file_name), 'w')

		file		= open(self.file_name, 'r')

		UP = 1
		DOWN = 1
		UP_data = []
		for i, line in enumerate(file):
			vec = [m for m in line.split(' ') if m != '' and m != '\n']
	
			if UP:			file_up.write(line)
			if DOWN:		file_down.write(line)
			if vec == []:
				if UP and DOWN: DOWN = 0
				elif UP: 		DOWN = 1; UP = 0

		file_up.close()
		file_down.close()


	def plot(self, file_name=None):

		X, Y, Z = np.mgrid[
							-0:10:complex(0, self.grid[0]),
							-0:10:complex(0, self.grid[1]), 
							-0:10:complex(0, self.grid[2])]


		values = self.CHG_grid + self.MAG_grid

		spacing = 4
		values 	= values[0::spacing,0::spacing,0::spacing]
		X 		= X[0::spacing,0::spacing,0::spacing]
		Y 		= Y[0::spacing,0::spacing,0::spacing]
		Z 		= Z[0::spacing,0::spacing,0::spacing]

		
		# Create figure
		fig = go.Figure()

		# Add traces, one for each slider step
		for step in [100, 50, 30, 10, 7, 5, 3, 1, 0.5, 0.1]:
		    fig.add_trace(
		        go.Isosurface(
					    x=X.flatten(),
					    y=Y.flatten(),
					    z=Z.flatten(),
					    value=values.flatten(),
					    colorscale='BlueRed',
					    isomin=step,
					    isomax=step,
					    surface_count=1,
					    caps=dict(x_show=False, y_show=False)
					    ))

		# Make 10th trace visible
		fig.data[0].visible = True

		# Create and add slider
		steps = []
		for i in range(len(fig.data)):
		    step = dict(
		        method="restyle",
		        args=["visible", [False] * len(fig.data)],
		    )
		    step["args"][1][i] = True  # Toggle i'th trace to "visible"
		    steps.append(step)

		sliders = [dict(
		    active=0,
		    currentvalue={"prefix": "ISOVALUE: "},
		    pad={"t": 6},
		    steps=steps
		)]

		fig.update_layout(
		    sliders=sliders
		)

		fig.show()


	def plot2(self, file_name=None):

		X, Y, Z = np.mgrid[
							-0:10:complex(0, self.grid[0]),
							-0:10:complex(0, self.grid[1]), 
							-0:10:complex(0, self.grid[2])]


		values = self.CHG_grid

		spacing = 5
		values 	= values[0::spacing,0::spacing,0::spacing]
		X 		= X[0::spacing,0::spacing,0::spacing]
		Y 		= Y[0::spacing,0::spacing,0::spacing]
		Z 		= Z[0::spacing,0::spacing,0::spacing]

		fig = go.Figure(data=go.Isosurface(
		    x=X.flatten(),
		    y=Y.flatten(),
		    z=Z.flatten(),
		    value=values.flatten(),
		    colorscale='BlueRed',
		    isomin=0.0002,
		    isomax=0.002,
		    surface_count=3,
		    caps=dict(x_show=False, y_show=False)
		    ))
		print(X.shape, values.shape)
		fig.show()

#parchg = PARCHG()
#parchg.direct_copy(file_name='files/orbitales/save/PARCHG.0088.ALLK')
#parchg.load('files/orbitales/save/PARCHG.0089.ALLK')
#parchg.plot()








