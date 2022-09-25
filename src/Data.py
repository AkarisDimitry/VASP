print('Python 3.X - https://www.python.org/download/releases/3.0/')
###################################
# Step 1 || Load python libraries #
###################################
# *** warning supresion
import multiprocessing as mp
import warnings
warnings.filterwarnings("ignore")

# *** numeric libraries *** #
try:
	import numpy as np
except: print('WARNING :: DATA.import_libraries() :: can not import numpy ')

import scipy.io

try:
	from sklearn.decomposition import PCA
	from sklearn.preprocessing import StandardScaler
	from sklearn import linear_model
	from sklearn.model_selection import cross_val_predict
	from sklearn.metrics import mean_squared_error, r2_score
except:	print('ERROR :: DATA.import_libraries() :: can not import sklearn')

try:
	from os import path
	import itertools, operator, logging, time, copy, pickle
except:  print('WARNING :: DATA.import_libraries() :: can not import itertools, operator, logging, time, copy, pickle or os')


try:
	import matplotlib.pyplot as plt
	import matplotlib.axes as ax
	import matplotlib.patches as patches
except:	print('WARNING :: DATA.import_libraries() :: can not import matplotlib ')

try:	from src import POSCAR
except:	
	try: import POSCAR as POSCAR
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

try:	from src import DOSCAR
except:	
	try: import DOSCAR as DOSCAR
	except: print('WARNING :: Set.import_libraries() :: can not import DOSCAR ')

try:	from src import OUTCAR
except:	
	try: import OUTCAR as OUTCAR
	except: print('WARNING :: Set.import_libraries() :: can not import OUTCAR ')

try:	from src import OSZICAR
except:	
	try: import OSZICAR as OSZICAR
	except: print('WARNING :: Set.import_libraries() :: can not import OSZICAR ')

try:	from src import System
except:	
	try: import System as System
	except: print('WARNING :: Set.import_libraries() :: can not import System ')

try:	from src import ORR as ORR
except:	
	try: import ORR as ORR
	except: print('WARNING :: Set.import_libraries() :: can not import ORR ')

from scipy.signal import savgol_filter

try:
	from ase import Atoms
	from ase.visualize import view
except:	print('ERROR :: DATA.import_libraries() :: can not import ase ')

class Data(object):
	def __init__(self, name=None, path=None):
		self.name = name
		self.system = {}
		self.path = path

		self.analysis_results = {}

	def load(self, data:dict=None, files:str='all', v:bool:True ) -> list:
		return [self.load_system(name=name, path=path, files=files, v=v) for name, path in data.items() ]

	def load_system_list(self, name_list=None, path_list=None, files='all', v=True) -> list:
		return [ self.load_system(name=name, path=path, files=files, v=v) for name, path in zip(name_list, path_list) ]

	def load_system(self, name:str=None, path:str=None, files='all', v:bool=True) -> object:
		if name == None or not type(name) == str: 	name = 'no_name'
		if not type(path) == str: path = self.path

		system = System.System(name=name)

		if type(files) == str and files == 'all':
			system.load_all(path, v=v)
		elif type(files) == list:
			for file_load in files:		
				system.load(path='{}/{}'.format(path, file_load), file=file_load, v=v)

		self.system[str(name)] = system 

		return self.system

	def summary(self, feature={'ORR':True} ):
		for key, system in self.system.items():
			print( '    *=== >> * {} '.format(key) )
			print( '    |       | ' )
			system.summary()

		# ------- ORR ------- #
		if 'ORR' in feature and feature['ORR']:
			try: self.analysis_results
			except: self.reaction(reaction='ORR')
			print( '            * == ORR == ' )
			if 'ORR' in self.analysis_results:
				try:
					self.analysis_results['ORR'].summarise_steps()
					self.analysis_results['ORR'].summarise_absortion()
				except:
					print('WARNING :: Data.summary() :: can NOT print summarise ORR')

	def reaction(self, reaction='ORR',):
		if reaction == 'ORR':	reaction = self.ORR()
		
		return reaction

	def ORR(self, ZPE=None, S=None, Gb=0, thermodimanic_corrections=True):
		try: self.analysis_results
		except: self.analysis_results = {}

		ZPE = ZPE if type(ZPE) == dict else {'sys':0, 'O':0.07, 'OH':0.35, 'OOH':0.43, 'O2':0.00, 'H2':0.27, 'H2O':0.56} 
		S = S if type(S) == dict else 		{'sys':0, 'O':0.0, 'OH':0.0, 'OOH':0.0, 'O2':0.64, 'H2':0.41, 'H2O':0.67}

		sys = self.system[list(filter(lambda x: (x[-1:] == '*'), self.system))[0]]
		if sys.ZPE == None:	sys.ZPE	= ZPE['sys'] # if None use default value # 
		if sys.S == None:	sys.S 	= S['sys'] # if None use default value # 

		sys_O = self.system[list(filter(lambda x: (x[-2:] == '*O'), self.system))[0]]
		if sys_O.ZPE == None:	sys_O.ZPE	= ZPE['O'] # if None use default value # 
		if sys_O.S == None:	sys_O.S 		= S['O'] # if None use default value # 

		sys_OH = self.system[list(filter(lambda x: (x[-3:] == '*OH'), self.system))[0]]
		if sys_OH.ZPE == None:	sys_OH.ZPE	= ZPE['OH'] # if None use default value # 
		if sys_OH.S == None:	sys_OH.S 	= S['OH'] # if None use default value # 

		sys_OOH = self.system[list(filter(lambda x: (x[-4:] == '*OOH'), self.system))[0]]
		if sys_OOH.ZPE == None:	sys_OOH.ZPE	= ZPE['OOH'] # if None use default value # 
		if sys_OOH.S == None:	sys_OOH.S 	= S['OOH'] # if None use default value # 

		try:
			sys_O2 = self.system[list(filter(lambda x: (x[-2:] == 'O2'), self.system))[0]]
			if sys_O2.ZPE == None:	sys_O2.ZPE	= ZPE['O2'] # if None use default value # 
			if sys_O2.S == None:	sys_O2.S 	= S['O2'] # if None use default value # 
		except: sys_O2 = None

		sys_H2 = self.system[list(filter(lambda x: (x[-2:] == 'H2'), self.system))[0]]
		if sys_H2.ZPE == None:	sys_H2.ZPE	= ZPE['H2'] # if None use default value # 
		if sys_H2.S == None:	sys_H2.S 	= S['H2'] # if None use default value # 

		sys_H2O = self.system[list(filter(lambda x: (x[-3:] == 'H2O'), self.system))[0]]
		if sys_H2O.ZPE == None:	sys_H2O.ZPE	= ZPE['H2O'] # if None use default value # 
		if sys_H2O.S == None:	sys_H2O.S 	= S['H2O'] # if None use default value # 

		try:
			self.ORR_analysis = ORR.OxigenReaction( name='Oxigen_Reaction', 	sys=sys, sys_O=sys_O, sys_OH=sys_OH, sys_OOH=sys_OOH, sys_O2=sys_O2,
																H2O=sys_H2O, H2=sys_H2, T=298)
			self.ORR_analysis.calculate( Gb=Gb, thermodimanic_corrections=thermodimanic_corrections)
			
			self.analysis_results['ORR'] = self.ORR_analysis

		except OSError as err:
			print('ERROR :: ORR_analysis :: Data.ORR() ')
			print("OS error: {0}".format(err))
		except ValueError:
			print('ERROR :: ORR_analysis :: Data.ORR() ')
			print("Could not convert data to an integer.")
		except:
			print('ERROR :: ORR_analysis :: Data.ORR() ')
			print("Unexpected error:", self.exc_info()[0])
			raise

		return self.ORR_analysis

	def plot(self,	data, 
					color=				{'frame':'#444444', 'zero_line':'#222222'}, 
					figure=				{'fig_n':0}, 
					label=				{}, 
					xlabel=				{'font_size':10}, 
					ylabel=				{'font_size':10}, 
					plot_limits=		{}, 
					ticks=				{'font_size':10}, 
					step_dimentions= 	{'step_width':1, 'step_sparce':0},  # ancho, separacion
					steps_names=		{'show':True, 'font_size':30, 'color':'#333333'}, 
					delta=				{'show': True, 'font_size':10},
					OP_plot=			{'show':True, 'font_size':10}, 
					plot=				{'lw': 4},
					save=				{},
					v_lines=			{'show':True,}):
		'''
		color=				{'h_lines':color, 'v_lines':color, 'zero_line':color, 
							'x_text':color, 'y_text':color, 'label_text':color,
							'frame':color}, 
		figure=				{'fig':fig, 'ax':ax}, 
		label=				{'title':str, 'font_size':40, 'color':color}, 
		xlabel=				{'label':'', 'font_size':40, 'color':'#111111', 'backgroundcolor':'#ffffff', 'weight':'roman'}, 
		ylabel=				{'label':str, 'font_size':40, 'color':color}, 
		system_name=		None, 
		plot_limits=		{'Y':[min, max], 'x':[min, max]}, 
		ticks=				{'font_size':18}, 
		step_dimentions= 	{'step_width':1, 'step_sparce':0},  # ancho, separacion
		steps_names=		{'show':True, 'label':[],'font_size':30, 'color':color}, 
		delta=				{'show': False, 'font_size':40, 'color':'#111111', 'backgroundcolor':'#ffffff', 'weight':'roman'},
		OP_plot=			{'show':False, 'font_size':30, 'backgroundcolor':30, 'color':30, 
							'weight':30, 'size':30, 'verticalalignment':30, 'alpha':30, 
							'facecolor':30, 'font_size':30, 'font_size':30, }, 
		plot=				{'lw': 7},
		save=				{'name':'', 'folder':'', 'ext':'', 'dpi':'' },
		v_lines=			{'show':True, 'linestyle':'--', 'linewidth':1.0}
		'''

		# ---- data read ---- # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
		if not type(data) is np.ndarray:
			try: data = np.array(data)
			except: print('ERROR :: code 0?? :: data plot need DATA as argument ')

		# ---- Figure configuration ---- # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
		fig, ax = plt.subplots() if type(figure) is dict else [figure['fig'], figure['ax']]
		fig, ax = [figure['fig'], figure['ax']]

		# ---- some global parameters ---- # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
		try:	step_width, step_sparce = step_dimentions['step_width'], step_dimentions['step_sparce']
		except: print('ERROR :: key error [step_width] or [step_sparce]')
		step_com = step_width + step_sparce
		step_number = data.shape[0]

		MIN, MAX = np.min(data), np.max(data) 
		DELTA = MAX - MIN

		color = { key:'#EE5555' if not key in color else color[key] for key in ['h_lines', 'v_lines', 'zero_line', 'x_text', 'y_text', 'label_text', 'frame'] }  

		label['title'] = label['title'] if 'title' in label else '' 

		v_lines = { key:value if not key in v_lines else v_lines[key] for (key,value) in {'show':True, 'linestyle':'-', 'linewidth':4.0, 'color':'#EE5555'}.items()}
		plot['lw'] = 4.0 if not 'lw' in plot else plot['lw']

		ticks['x_font_size'] = 20.0 if not 'x_font_size' in ticks else ticks['x_font_size']
		ticks['y_font_size'] = 15.0 if not 'y_font_size' in ticks else ticks['y_font_size']

		# ---- FRAME ---- # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
		for spine in ax.spines: ax.spines[spine].set_linewidth(7); ax.spines[spine].set_color(color['frame'])

		# ---- PLOT limits ---- # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
		if 'X' in plot_limits:	ax.set_xlim(plot_limits['X'][0], plot_limits['X'][1])
		else: 					ax.set_xlim((0, step_number*(step_com) ))

		if 'Y' in plot_limits:	ax.set_ylim( [plot_limits['Y'][0], plot_limits['Y'][1]] )
		else:					ax.set_ylim( [MIN-0.5*DELTA, MAX+0.5*DELTA] )

		# ---- reference lvl ---- # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
		ax.plot([-1, step_number*(step_com)+step_width], [0,0], '--', lw = 4, color=color['zero_line'], alpha=0.5)

		# ---- Y-TICKS ---- # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
		y_ticks = ticks['yticks'] if 'yticks' in ticks else [ float('{:.2}'.format(n)) for n in np.linspace(start=MIN-0.5*DELTA, stop=MAX+0.5*DELTA, num=10)] 
		ax.set_yticks(y_ticks, y_ticks)
		ax.set_yticklabels(labels=y_ticks, fontdict = { 'family': 'serif', 'color':  'darkred', 'weight': 'normal', 'size': ticks['y_font_size'],})

		# ---- X-TICKS ---- # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
		x_ticks = ticks['xticks'] if 'xticks' in ticks else [float('{:.2}'.format(i*step_com+step_width/2)) for i, n in enumerate(data)]  
		x_ticks_labels = ticks['x_ticks_labels'] if 'x_ticks_labels' in ticks else [i+1 for i, n in enumerate(data)] 
		ax.set_xticks(x_ticks, x_ticks_labels)
		ax.set_xticklabels(labels=x_ticks, fontdict = { 'family': 'serif', 'color':  'darkred', 'weight': 'normal', 'size': ticks['x_font_size'],})

		# ---- PLOT ---- # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
		# ---- PLOT horizontal lines ---- #
		patches.Rectangle((1,1),2,2, color='#AAAAAA')
		
		steps_names = steps_names if type(steps_names) == dict else {} 
		
		steps_names['label'] = steps_names['label'] if 'show' in steps_names and steps_names['show'] and 'label' in steps_names and len(steps_names['label']) >= step_number else np.arange(1, step_number+1)
		steps_names['color'] = steps_names['color'] if 'color' in steps_names else '#333333' 
		steps_names = steps_names if type(steps_names) == dict else {} 

		for i, value in enumerate(data):
			ax.plot([i*step_com, i*step_com+step_width], [value, value], color=color['h_lines'], lw=plot['lw'] )

			if 'show' in steps_names and steps_names['show']:
				text_font_size = steps_names['font_size'] if 'font_size' in steps_names else 30

				ax.text(x=i*step_com+step_width/2 , y=value+DELTA*0.1, s=str(steps_names['label'][i]),
					            backgroundcolor='#FFFFFF', color=steps_names['color'], weight='roman', horizontalalignment='center',
					            size=text_font_size, alpha=0.9,
					            bbox=dict(facecolor='red', alpha=0.0)) 

		# ---- PLOT vertical lines ---- #
		if 'show' in v_lines and v_lines['show']:
			for n in range(data.shape[0]-1):
				ax.plot([(n+1)*(step_width)+n*step_sparce, (n+1)*(step_width+step_sparce)], [data[n], data[n+1]], 
						color=v_lines['color'], 
						linewidth = v_lines['linewidth'], 
						linestyle= v_lines['linestyle'])

		# ---- TITLE ---- #---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
		label  = { key:value if not key in label else label[key] for (key,value) in {'label':'', 'font_size':40, 'color':'#111111', 'backgroundcolor':'#ffffff', 'weight':'roman'}.items()}
		xlabel = { key:value if not key in xlabel else xlabel[key] for (key,value) in {'label':'', 'font_size':40, 'color':'#111111', 'backgroundcolor':'#ffffff', 'weight':'roman'}.items()}
		ylabel = { key:value if not key in ylabel else ylabel[key] for (key,value) in {'label':'', 'font_size':40, 'color':'#111111', 'backgroundcolor':'#ffffff', 'weight':'roman'}.items()}

		ax.set_title('{}'.format(label['label']), backgroundcolor='#ffffff', color=label['color'], weight='roman', size=label['font_size'], pad=30)

		# --- Axis label --- #---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
		ax.set_xlabel(str(xlabel['label']), backgroundcolor='#ffffff', color=xlabel['color'], weight='roman', size=xlabel['font_size'], labelpad=10) 	
		ax.set_ylabel(str(ylabel['label']), backgroundcolor='#ffffff', color=ylabel['color'], weight='roman', size=ylabel['font_size'], labelpad=20)

		# ---- REFERENCES ---- #
		#for i, n in enumerate(range(data.shape[1])):
		#	plt.figtext(0.84 , 0.76+0.05*i, '                                                                                         ',
		#	            backgroundcolor=color[i], color='black', weight='roman',
		#	            size=5)
		#	plt.figtext(0.90, 0.75+0.05*i, system_name[i],
		#	            backgroundcolor='#ffffff', color='black', weight='roman',
		#	            size=26, )

		# ---- DELTA value ---- # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
		delta = { key:value if not key in delta else delta[key] for (key,value) in {'show':True, 'label':'', 'font_size':5, 'color':'#111111', 'backgroundcolor':'#ffffff', 'weight':'roman', 'alpha':0.0}.items()}

		if 'show' in delta and delta['show']:
			for i in range(step_number-1):
				t = ax.text(i*(step_com)+1.2*step_width, (data[i]+data[i+1])/2 , '{0:> 2.2f}'.format(data[i+1]-data[i]),
								backgroundcolor=delta['backgroundcolor'], color=delta['color'], weight=delta['weight'],
								size=delta['font_size']) # transform=ax.transAxes
				t.set_bbox(dict(facecolor=color['v_lines'], alpha=0.0, edgecolor=color['v_lines'] ))

		# ---- overpotencial value PLOT ---- # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
		OP_plot = { key:value if not key in OP_plot else OP_plot[key] for (key,value) in {'show':True, 'label':'', 'font_size':10, 'color':'#111111', 'backgroundcolor':'#ffffff', 'weight':'roman', 'alpha':0.0, 'verticalalignment':'center'}.items()}

		if 'show' in OP_plot and OP_plot['show']:
			OP_j = np.max(-(data[:-1] - data[1:]))
			OP_j_arg = np.argmax(-(data[:-1] - data[1:]))

			X0 = OP_j_arg*(step_com)+step_width*1.1
			Y0, Y1 = data[OP_j_arg+1], data[OP_j_arg+1] - OP_j

			t = ax.text(X0 + 0.1*step_com, (data[OP_j_arg]+data[OP_j_arg+1])/2 , '{0:> 2.3f}'.format(data[OP_j_arg+1]-data[OP_j_arg]),
							backgroundcolor=OP_plot['backgroundcolor'], color=OP_plot['color'], weight=OP_plot['weight'],
							size=OP_plot['font_size'], verticalalignment=OP_plot['verticalalignment'], alpha=OP_plot['alpha'],
							bbox=dict(facecolor='red', alpha=0.0)) # transform=ax.transAxes

			ax.annotate(s='', xy=(X0, Y0), xytext=(X0, Y1), size=20, arrowprops=dict(arrowstyle='<->', color=color['v_lines'], lw=2, ), color=color['v_lines'], alpha=0.3 )

		# ---- save ---- # ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----  ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 		
		try:
			save['name'] = save['name'] if 'name' in save else 'saved_fig'+list(filter(lambda x: (x[-1:] == '*'), self.system))[0]
			save['folder'] = save['folder'] if 'folder' in save else '.'
			save['ext'] = save['ext'] if 'ext' in save else 'png'
			save['dpi'] = save['dpi'] if 'dpi' in save else 300

			ax.figure.savefig('{}/{}.{}'.format(save['folder'], save['name'], save['ext']), bbox_inches='tight', dpi=save['dpi'])
		except:
			print('ERROR :: code 0?? :: can NOT save figure ')

	def get_embeddings(self, path=None, embedding_parameters=None, processing_pool={'multiprocessing':True, 'cores':'check'}, save=True):
		'''
		| Multi-args   Concurrence    Blocking     Ordered-results
		---------------------------------------------------------------------
		Pool.map          | no           yes            yes          yes
		Pool.map_async    | no           yes            no           yes
		Pool.apply        | yes          no             yes          no
		Pool.apply_async  | yes          yes            no           no
		Pool.starmap      | yes          yes            yes          yes
		Pool.starmap_async| yes          yes            no           no
		'''

		def log_result(result):    result_list.append(result)

		# === Get data from OUTCAR === #
		embedding_parameters = {'partition':400, 'min distance':0.5, 'max distance':7.0, 'range distance':0.001} if embedding_parameters == None else embedding_parameters

		# === Processing pool segmentation initialization === # 
		if processing_pool['multiprocessing']: 
			if type(processing_pool['cores']) == int:
				pool = mp.Pool(processing_pool['cores'])
			if type(processing_pool['cores']) == str:
				pool = mp.Pool(mp.cpu_count())
				print(f' Number of cores detected : {mp.cpu_count()}')
		
		result_list = []

		for i, (key, system) in enumerate(self.system.items()):
			# ===  embedding segmentation by core === # 
			if processing_pool['multiprocessing']:
				if len(np.array(system.OUTCAR.E).shape) > 0 and np.array(system.OUTCAR.E).shape[0] < 400:
					pool.apply_async(system.get_embeddings, args = [embedding_parameters], callback = log_result )	
				else: pass

			else:
				# === one core processing === # 
				if len(np.array(system.OUTCAR.E).shape) > 0 and np.array(system.OUTCAR.E).shape[0] < 400:
					embedding_dict = system.get_embeddings(embedding_parameters=embedding_parameters)
					
					embedding_vector, force_vector, ff_X, ff_Y = embedding_dict['embedding_vector'], embedding_dict['force_vector'], embedding_dict['ff_X'], embedding_dict['ff_Y']
					result_list.append(embedding_dict)
				else: pass
		'''
		=== Eg FeTPyP 79 steps ===
		Getting embedding :: {'partition': 400, 'min distance': 0.5, 'max distance': 7.0, 'range distance': 0.001}
		Processing :: 79 steps
		Time 9.389496564865112 
			(0) Ion:Fe Samples:78 Dimention:1600 
			(1) Ion:N Samples:624 Dimention:1600 
			(2) Ion:C Samples:2496 Dimention:1600 
			(3) Ion:H Samples:1248 Dimention:1600 

		embedding_vector.shape :	(78, 57, 4, 400) 	[steps, atoms, atomsID, embedding size] 
		force_vector :				(78, 57, 3) 		[steps, atoms, 3D]
		ff_X : 						(2496, 1600) 		[steps*atoms_i, atomsID*embedding size_i]
		ff_Y : 						(2496, 3)			[steps*atoms_i, 3D]
		
		'''

		# === synchronizing core date === # 
		if processing_pool['multiprocessing']:
			pool.close()
			pool.join()

		# === Store data === #
		if save: 	
			self.embedding_result = result_list
			if path != None:
				pass

		return result_list

	def get_absortion(self, analysis='ORR', **karg):
		if analysis=='ORR':
			return self.ORR_analysis.get_absortion(**karg)
		return None

	def exc_info(self, ): 
		return ['Can not properlly calculate ORR parameters.']


