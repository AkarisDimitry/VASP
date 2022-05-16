# ------------------------------------- #
# ---  XML XML XML XML XML XML XML  --- #
# ------------------------------------- #
print('Python 3.X - https://www.python.org/download/releases/3.0/')
#####################################
# 		 Load python libraries 		#
#####################################

# === === === import libraries === === === #
# *** warning supresion
import warnings
warnings.filterwarnings("ignore")

# *** warning supresion
import warnings; warnings.filterwarnings("ignore")

# *** python libraries
try:
	import itertools, operator, logging, time, copy, pickle, datetime, os, sys, argparse
	#os.chmod('.', 777)
except:  print('WARNING :: XML.import_libraries() :: can not import itertools, operator, logging, time, copy, pickle or os')

# *** numpy libraries
try:
	import numpy as np
except: print('WARNING :: XML.import_libraries() :: can not import numpy ')


# XML #
import openpyxl
from openpyxl.utils import get_column_letter
from openpyxl.styles import Font, Border, Side, Color, PatternFill, NamedStyle
from openpyxl.chart import BarChart, Series, Reference

# PANDAs #
import pandas as pd

# PLOT #
# *** matplotlib libraries 
try:
	# seaborn #
	import seaborn as sns

	# matplotlib #
	import matplotlib.pyplot as plt
	import matplotlib.axes as ax
	import matplotlib.patches as patches

	# ASE #
	from ase import Atoms
	from ase.visualize import view
	from ase.visualize.plot import plot_atoms
	from ase.lattice.cubic import FaceCenteredCubic
except:	print('WARNING :: XML.import_libraries() :: can not import matplotlib ')


try:	from src import POSCAR
except:	
	try: import POSCAR as POSCAR
	except: print('WARNING :: XML.import_libraries() :: can not import POSCAR ')

try:	from src import DOSCAR
except:	
	try: import DOSCAR as DOSCAR
	except: print('WARNING :: XML.import_libraries() :: can not import DOSCAR ')

try:	from src import OUTCAR
except:	
	try: import OUTCAR as OUTCAR
	except: print('WARNING :: XML.import_libraries() :: can not import OUTCAR ')

try:	from src import OSZICAR
except:	
	try: import OSZICAR as OSZICAR
	except: print('WARNING :: XML.import_libraries() :: can not import OSZICAR ')

try:	from src import Data
except:	
	try: import Data as Data
	except: print('WARNING :: XML.import_libraries() :: can not import Data ')

try:	from src import System
except:	
	try: import System as System
	except: print('WARNING :: XML.import_libraries() :: can not import System ')

try:	from src import ORR
except:	
	try: import ORR as ORR
	except: print('WARNING :: XML.import_libraries() :: can not import ORR ')

try:
	from scipy.signal import savgol_filter

	# NON linear models #
	from sklearn.neural_network import MLPRegressor
	from sklearn.datasets import make_regression
	from sklearn.model_selection import train_test_split
	from sklearn.model_selection import ShuffleSplit # or StratifiedShuffleSplit
	# linear models #
	from sklearn.decomposition import PCA
	from sklearn.preprocessing import StandardScaler
	from sklearn import linear_model
	from sklearn.model_selection import cross_val_predict
	from sklearn.metrics import mean_squared_error, r2_score

except:	print('ERROR :: DATA.import_libraries() :: can not import sklearn')

class XML(object):
	def __init__(self, 	name=None, sheet=None, 
						DATA=None, path=None, ):
		self.name  = name
		self.sheet = []
		self.path  = path

		# == DATA == #
		self.set = DATA
		
		self.ORR_array = None
		self.ORR_df    = None
		self.ORR_relevant = [	'G1_ORR','G2_ORR','G3_ORR','G4_ORR','limiting_step_ORR',
								'Eabs_O','Gabs_O','Eabs_OH','Gabs_OH','Eabs_OOH','Gabs_OOH',
								'overpotencial_ORR_4e']
		
		self.Magnetization_array = None
		self.Magnetization_df    = None

		self.Distance_array = None
		self.Distance_df    = None
		
		self.WB = openpyxl.Workbook()

	def make_path(self, path):
		if not os.path.isdir(path):  	 			
			os.makedirs(path)
			return True

		return False

	def creat_styles(self, ):
		highlight = NamedStyle(name="highlight")
		highlight.font = Font(bold=True, size=20)
		bd = Side(style='thick', color="000000")
		highlight.border = Border(left=bd, top=bd, right=bd, bottom=bd)
		self.WB.add_named_style(highlight)

	def plot_bars(self, data, x, y, hue, 
					   Title='', path=None, name='plot.png',
					    ):

		if type(data) != type(None) and len(data) > 0:
			ax = sns.barplot(
			    x=x, 
			    y=y, 
			    hue=hue, 
			    data=data, 
			    ci="sd", 
			    edgecolor="black",
			    errcolor="black",
			    errwidth=1.5,
			    capsize = 0.1,
			    alpha=0.5
			)
			sns.stripplot(
			    x=x, 
			    y=y, 
			    hue=hue, 
			    data=data, dodge=True, alpha=0.6, ax=ax
			)

			handles, labels = ax.get_legend_handles_labels()
			ax.legend(handles, labels, title=hue, bbox_to_anchor=(1, 1.02), loc='upper left')
			ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
			ax.set_title(Title)
			fig = ax.get_figure()
			if not os.path.exists(f"{path}/{name}"):
				fig.savefig(f"{path}/{name}" , dpi=100, pad_inches=0.1, bbox_inches='tight', horizontalalignment='right') 
			plt.clf()

		return ax

	def plot_add(self, path, name, anchor, sheet, height=400, width=400):
		img = openpyxl.drawing.image.Image(f"{path}/{name}")
		img.anchor = f'{get_column_letter(i*5+1)}{line*20+1}'
		img.height = 400 # insert image height in pixels as float or int (e.g. 305.5)
		img.width  = 400 # insert image width in pixels as float or int (e.g. 405.8)
		sheet.add_image(img)

		return True

	def summary(self, ):
		# fill the workbook with all possible data from DATA.py #

		# ====== Create folthers ====== #
		self.make_path(f'{self.path}/plots/geometrias/')
		self.make_path(f'{self.path}/plots/ORR/')
		self.make_path(f'{self.path}/plots/metanalisis/name1')
		self.make_path(f'{self.path}/plots/metanalisis/name2')
		self.make_path(f'{self.path}/plots/metanalisis/name3')
		self.make_path(f'{self.path}/plots/metanalisis/magnetization')
		self.make_path(f'{self.path}/plots/metanalisis/magnetization/System')
		self.make_path(f'{self.path}/plots/metanalisis/magnetization/Functiona/')
		self.make_path(f'{self.path}/plots/metanalisis/distance')
		self.make_path(f'{self.path}/plots/metanalisis/distance/System')
		self.make_path(f'{self.path}/plots/metanalisis/distance/Functiona/')

		self.get_ORR_array()
		self.get_ORR_df()

		self.get_Magnetization_array()
		self.get_Magnetization_DF()

		self.get_Distance_array()
		self.get_Distance_DF()



		self.plot_bars(data, x, y, hue, Title='', path=None, name='plot.png')
		self.plot_add(path, name, anchor, sheet, height=400, width=400)

	def summary_system(self, ):
		for set_n, (key_data, data) in enumerate(self.set.items()):
			if not os.path.isdir(f'{path}/plots/systemas/{key_data}'):  	  os.makedirs(f'{path}/plots/systemas/{key_data}')	

			# == make new sheet == (with an acceptable name)
			sheet_name = key_data.replace('*','x') 
			sheet_name = sheet_name.replace('[','L') 
			sheet_name = sheet_name.replace(']','L') 

			ws = wb.create_sheet(sheet_name) # insert at the end (default)
			sheet = wb.get_sheet_by_name(sheet_name)

			# == set some sheer variables ==
			sheet.column_dimensions['A'].width = 30
			head_space = 37 # space from first line
			line = 1+head_space 
			for system_n, (key_system, system) in enumerate(data.system.items()):
				# ======================= PLOT ATOMS ======================= #
				if system.CONTCAR != None:	
					try:
						#if not os.path.exists(f'{path}/plots/systemas/{key_data}/{key_system}.png'):
							#system.CONTCAR.plot(save=True, path=f'{path}/plots/systemas/{key_data}/{key_system}.png')
						img = openpyxl.drawing.image.Image(f'{path}/plots/systemas/{key_data}/{key_system}.png')
						img.anchor = f'A{line}'
						sheet.add_image(img)
					except: pass

	def get_ORR_array(self, save=True):
		# set conteiners #
		ORR_relevant = self.ORR_relevant
		ORR_array   = []
		for set_n, (key_data, data) in enumerate(self.set.items()):
			# -- generate OBJ and alalize the ORR reaction -- #
			data.ORR()
			ORR = data.reaction('ORR')

			# -- generate plot -- #
			if not os.path.exists(f"{self.path}/plots/ORR/{key_data}.png"):
				ORR.plot(folder=f'{self.path}/plots/ORR/.', name=key_data)

			# -- store DATA in vector-- #
			if type(ORR.ORR) == dict:	 # store ORR info of all systems #
				vector_labels  = [key_data, key_data.split('_')[0], key_data.split('_')[1], '_'.join(key_data.replace('__','_').split('_')[2:]) ] 
				vector_numeric = [ORR.ORR[label] for label in ORR_relevant]

			# == ORR_array data append == #
			ORR_array.append(vector_numeric+vector_labels)
		ORR_array = np.array(ORR_array) # store ORR info of all systems #

		if save:	
			self.ORR_array = ORR_array
			
		return ORR_array

	def get_ORR_DF(self, ORR_array=None, ORR_relevant=None, save=True):
		ORR_array = ORR_array if type(ORR_array) == type(None) else self.ORR_array
		ORR_relevant = ORR_relevant if type(ORR_relevant) == type(None) else self.ORR_relevant

		# -- make ORR dataframe -- #
		ORR_df  	= pd.DataFrame(data=ORR_matrix,  columns=ORR_relevant+['FullID', 'System', 'Ligand', 'Cell', 'Functional']) # store ORR info of all systems #
		ORR_df[ORR_relevant] = ORR_df[ORR_relevant].apply(pd.to_numeric)
		if save:	self.ORR_df = ORR_df
		return ORR_df

	def get_Magnetization_array(self, ):
		# set conteiners #
		Magnetization_array = []
		
		# iterate each data conteiner #
		for set_n, (key_data, data) in enumerate(self.set.items()):
			vector_labels  = [key_data, key_data.split('_')[0], key_data.split('_')[1], '_'.join(key_data.replace('__','_').split('_')[2:]) ] 

			# iterate each system  #
			for system_n, (key_system, system) in enumerate(data.system.items()):
				key_subtype = key_system.split('_')[-1] 
				try:
					if key_subtype in ['*', '*O', '*OH', '*OOH']:
						# -- store DATA in vector-- #
						atoms_names =  np.array(system.OUTCAR.atoms_names_list)[system.CONTCAR.relevant_atoms_mask()] 
						atoms_magnetization = np.abs(system.OUTCAR.resume()['dict']['magnetization'][system.CONTCAR.relevant_atoms_mask()][:,-1])
						Magnetization_array += [ vector_labels+[key_subtype, atoms_names[atom_i], atom_i, atoms_magnetization[atom_i]] for atom_i, atom_name in enumerate(atoms_names)  ]	
				except: 
					print('ERROR :: XML.get_Magnetization_array() :: Can not Extrar magnetization data from dataset')
		
		if save: self.Magnetization_array = Magnetization_array
			
		return Magnetization_array

	def get_Magnetization_DF(self, Magnetization_array=None, save=True):
		Magnetization_array = Magnetization_array if type(Magnetization_array) == type(None) else self.Magnetization_array
		Magnetization_array = Magnetization_array if type(Magnetization_array) == type(None) else self.Magnetization_array

		# -- make MARNETIZATION dataframe -- #
		Magnetization_array = np.array(Magnetization_array)
		Magnetization_df  	 = pd.DataFrame(data=Magnetization_array,  columns= ['FullID', 'System', 'cell', 'Functional']+['System', 'Atom', 'number', 'Magnetization'] ) # store ORR info of all systems #
		Magnetization_df['Magnetization'] = Magnetization_df['Magnetization'].apply(pd.to_numeric)
		
		if save: self.Magnetization_df = Magnetization_df
		
		return Magnetization_df

	def get_Distance_array(self, save=True ):
		# set conteiners #
		Distance_array = []

		# iterate each data conteiner #
		for set_n, (key_data, data) in enumerate(self.set.items()):
			vector_labels  = [key_data, key_data.split('_')[0], key_data.split('_')[1], '_'.join(key_data.replace('__','_').split('_')[2:]) ] 
			
			# iterate each system  #
			for system_n, (key_system, system) in enumerate(data.system.items()):
				key_subtype = key_system.split('_')[-1] 
				try:
					if key_subtype in ['*', '*O', '*OH', '*OOH']:
						# correct CONTCAR ID with OUTCAR data # 
						system.CONTCAR.atoms_names_list = system.OUTCAR.atoms_names_list
	
						# -- store DATA in vector-- #
						Distance_vec = [[ vector_labels+[key_subtype, atom1, atom2, prop['distances'][pos]] for pos, atom2 in enumerate(prop['names']) ] for atom1, atom1N in system.CONTCAR.relevant_distances().items()  for num, prop in atom1N.items()  ] 
						for n in Distance_vec:
							Distance_array += n
				except: pass

		if save: self.Distance_array = Distance_array
			
		return Distance_array

	def get_Distance_DF(self, Distance_array=None, save=True):
		# set conteiners #
		Distance_array = Distance_array if type(Distance_array) == type(None) else self.Distance_array
		
		# -- make DISTANCE dataframe -- #
		Distance_array = np.array(Distance_array)
		Distance_df  	 = pd.DataFrame(data=Distance_array,  columns= ['FullID', 'Name', 'cell', 'Functional']+['System', 'Atom1', 'Atom2', 'Distance'] ) # store ORR info of all systems #
		Distance_df['Distance'] = Distance_df['Distance'].apply(pd.to_numeric)

		if save: self.Distance_df = Distance_df

		return Distance_df
		

	def add_magnetization_sheet(self, Magnetization_df=None, path=None):
		path = path if type(path) != type(None) else self.path
		# == Magnetization Per ATOM/System/FUNCTIONAL == #  # == Magnetization Per ATOM/System/FUNCTIONAL == #  # == Magnetization Per ATOM/System/FUNCTIONAL == # 
		# -- Magnetization vs Functional -- #
		line = -2
		ws = wb.create_sheet(f'Mag_functional') 
		MAGsheet = wb.get_sheet_by_name(f'Mag_functional')
	
		for name_i in pd.unique(Magnetization_df['name1']):
			line += 1

			for atom_i in pd.unique(Magnetization_df['Atom']):
				line += 1

				for i, sys_i in enumerate(pd.unique(Magnetization_df['System'])):
					# == filter data == # 
					data  = Magnetization_df[Magnetization_df['Magnetization'] < 5][Magnetization_df['Atom'] == atom_i][Magnetization_df['System'] == sys_i][Magnetization_df['Name'] == name_i]

					if len(data) > 0:			
						self.plot_bars(data, x="Magnetization", y="Magnetization", hue="Atom", 
										Title=f'Systems {sys_i} | Atom Type {atom_i} |  {name_i} ', 
										path =f"{path}/plots/metanalisis/magnetization/Functional", 
										name =f"bars_plot_{sys_i}_{atom_i}_{name_i}.png")

						self.plot_add(	path  =f"{path}/plots/metanalisis/magnetization/Functional", 
										name  =f"bars_plot_{sys_i}_{atom_i}_{name_i}.png", 
										anchor=f'{get_column_letter(i*5+1)}{line*20+1}', 
										sheet =MAGsheet, height=400, width=400)						
				
		# == Magnetization Per ATOM/System == #	# == Magnetization Per ATOM/System == # # == Magnetization Per ATOM/System == #		
		# Magnetization vs System #	
		line = -1
		ws = wb.create_sheet(f'Mag_sys') # insert at the end (default)
		MAGsheet = wb.get_sheet_by_name(f'Mag_sys')

		for atom_i in pd.unique(Magnetization_df['Atom']):
			line += 1

			for i, sys_i in enumerate(pd.unique(Magnetization_df['System'])):
				data = Magnetization_df[Magnetization_df['Magnetization'] < 5][Magnetization_df['Atom'] == atom_i][Magnetization_df['System'] == sys_i]


				if len(data) > 0:			
					self.plot_bars(data, x="Magnetization", y="Magnetization", hue="Atom", 
									Title=f'Atom Type {atom_i} |  {name_i} ', 
									path =f"{path}/plots/metanalisis/magnetization/System", 
									name =f"bars_plot_{sys_i}_{atom_i}.png")

					self.plot_add(	path  =f"{path}/plots/metanalisis/magnetization/System", 
									name  =f"bars_plot_{sys_i}_{atom_i}.png", 
									anchor=f'{get_column_letter(i*5+1)}{line*20+1}', 
									sheet =MAGsheet, height=400, width=400)						
			
	def add_distance_sheet(self, Magnetization_df=None):
		# == Magnetization Per ATOM/System/FUNCTIONAL == #  # == Magnetization Per ATOM/System/FUNCTIONAL == #  # == Magnetization Per ATOM/System/FUNCTIONAL == # 

		# -- Magnetization vs Functional -- #
		line = -2
		ws = wb.create_sheet(f'Mag_functional') 
		MAGsheet = wb.get_sheet_by_name(f'Mag_functional')
	
		for name_i in pd.unique(Magnetization_df['name1']):
			line += 1

			for atom_i in pd.unique(Magnetization_df['Atom']):
				line += 1

				for i, sys_i in enumerate(pd.unique(Magnetization_df['System'])):
					# == filter data == # 
					data  = Magnetization_df[Magnetization_df['Magnetization'] < 5][Magnetization_df['Atom'] == atom_i][Magnetization_df['System'] == sys_i][Magnetization_df['Name'] == name_i]

					if len(data) > 0:			
						self.plot_bars(data, x="Magnetization", y="Magnetization", hue="Atom", 
										Title=f'Systems {sys_i} | Atom Type {atom_i} |  {name_i} ', 
										path =f"{path}/plots/metanalisis/magnetization/Functional", 
										name =f"bars_plot_{sys_i}_{atom_i}_{name_i}.png")

						self.plot_add(	path  =f"{path}/plots/metanalisis/magnetization/Functional", 
										name  =f"bars_plot_{sys_i}_{atom_i}_{name_i}.png", 
										anchor=f'{get_column_letter(i*5+1)}{line*20+1}', 
										sheet =MAGsheet, height=400, width=400)						
				
		# == Magnetization Per ATOM/System == #	# == Magnetization Per ATOM/System == # # == Magnetization Per ATOM/System == #		
		# Magnetization vs System #	
		line = -1
		ws = wb.create_sheet(f'Mag_sys') # insert at the end (default)
		MAGsheet = wb.get_sheet_by_name(f'Mag_sys')

		for atom_i in pd.unique(Magnetization_df['Atom']):
			line += 1

			for i, sys_i in enumerate(pd.unique(Magnetization_df['System'])):
				data = Magnetization_df[Magnetization_df['Magnetization'] < 5][Magnetization_df['Atom'] == atom_i][Magnetization_df['System'] == sys_i]


				if len(data) > 0:			
					self.plot_bars(data, x="Magnetization", y="Magnetization", hue="Atom", 
									Title=f'Atom Type {atom_i} |  {name_i} ', 
									path =f"{path}/plots/metanalisis/magnetization/System", 
									name =f"bars_plot_{sys_i}_{atom_i}.png")

					self.plot_add(	path  =f"{path}/plots/metanalisis/magnetization/System", 
									name  =f"bars_plot_{sys_i}_{atom_i}.png", 
									anchor=f'{get_column_letter(i*5+1)}{line*20+1}', 
									sheet =MAGsheet, height=400, width=400)						
			
						
