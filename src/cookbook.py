###################################
# Step 1 || Load python libraries #
###################################
# *** warning supresion
import warnings
warnings.filterwarnings("ignore")

# *** numeric libraries *** #
import numpy as np
import scipy.io
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn.model_selection import ShuffleSplit # or StratifiedShuffleSplit
from sklearn.decomposition import PCA

from sklearn import manifold

# *** graph libraries *** #
try:
	import matplotlib.pyplot as plt
	from mpl_toolkits import mplot3d
	import matplotlib as mpl
except: 
	print('WARNNING :: main_simulation.py :: can NOT correctly load "matplotlib" libraries')
	print('Install by: ( pip3 install matplotlib )')
	

# *** python common libraries *** #
import logging, operator, pickle, os, random

# *** load own libraries *** #
try:	from src import POSCAR
except:	
	try: from POSCAR import POSCAR
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

try:	from src import DOSCAR
except:	
	try: from DOSCAR import DOSCAR
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

try:	from src import PROCAR
except:	
	try: from PROCAR import PROCAR
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

try:	from src import OUTCAR
except:	
	try: from OUTCAR import OUTCAR
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

try:	from src import OSZICAR
except:	
	try: from OSZICAR import OSZICAR
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

try:	from src import System
except:	
	try: from System import System
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

try:	from src import Set
except:	
	try:	from Set import Set
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

try:	from src import Data
except:	
	try: from Data import Data
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

try:	from src import ORR
except:	
	try: from ORR import OxigenReaction
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')

#########################################################
# ORR eg || FePC/FePCAu/FePCBzAu  CoPC/CoPCAu/CoPCBzAu  #
#########################################################
def ORR_PC_plot():
	# ============== CoPC  ============== # 
	orr = OxigenReaction() # ========= CoPC FREE (VOID) ========= # 
	orr.calculate(sys={'E':-422.133,'ZPE':0.0,'S':0.0}, sys_O={'E':-426.500,'ZPE':0.07,'S':0.0},
				 sys_OH={'E':-432.0393,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-436.660,'ZPE':0.43,'S':0.0}, sys_O2=None, 
							H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0)
	orr.plot_integrated('CoPC FREE :: Ea={:.3}'.format(-422.133+421.979))
	orr.summarise_steps()
	orr.summarise_absortion()

	# ============== CoPC + 5Bz + Au ============== # 
	orr = OxigenReaction() # ========= CoPC + 5Bz + Au (Au site) ========= # 
	orr.calculate(sys={'E':-1120.709,'ZPE':0.0,'S':0.0}, sys_O={'E':-1125.299,'ZPE':0.07,'S':0.0},
				 sys_OH={'E':-1130.573,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-1135.187,'ZPE':0.43,'S':0.0}, sys_O2=None, 
							H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0)
	orr.plot_integrated('CoPC + 5Bz + Au (Au site) :: Ea={:.3}'.format(-1120.709+421.979+696.631) )
	orr.summarise_steps()
	orr.summarise_absortion()

	orr = OxigenReaction() # ========= CoPC + 5Bz + Au (C site) ========= # 
	orr.calculate(sys={'E':-1120.579,'ZPE':0.0,'S':0.0}, sys_O={'E':-1124.962,'ZPE':0.07,'S':0.0},
				 sys_OH={'E':-1130.464,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-1135.090,'ZPE':0.43,'S':0.0}, sys_O2=None, 
							H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0)
	orr.plot_integrated('CoPC + 5Bz + Au (C site) :: Ea={:.3}'.format(-1120.579+421.979+696.631) )
	orr.summarise_steps()
	orr.summarise_absortion()

	orr = OxigenReaction() # ========= CoPC + 5Bz + Au (Bz site) ========= # 
	orr.calculate(sys={'E':-1120.389,'ZPE':0.0,'S':0.0}, sys_O={'E':-1125.100,'ZPE':0.07,'S':0.0},
				 sys_OH={'E':-1130.270,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-1134.892,'ZPE':0.43,'S':0.0}, sys_O2=None, 
							H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0)
	orr.plot_integrated('CoPC + 5Bz + Au (Bz site) :: Ea={:.3}'.format(-1120.389+421.979+696.631) )
	orr.summarise_steps()
	orr.summarise_absortion()

	# ============== CoPC + Au ============== # 
	orr = OxigenReaction() # ========= CoPC + Au (HOLLOW site) ========= #
	orr.calculate(sys={'E':-739.3821,'ZPE':0.0,'S':0.0}, sys_O={'E':-743.649,'ZPE':0.07,'S':0.0}, 
				sys_OH={'E':-749.105,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-753.758,'ZPE':0.43,'S':0.0}, sys_O2=None, 
							H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0*4)
	orr.plot_integrated('CoPC + Au (HOLLOW site) :: Ea={:.3}'.format(-739.3821+421.979+312.16) )
	orr.summarise_steps()

	orr = OxigenReaction() # ========= CoPC + Au (TOP site) ========= #
	orr.calculate(sys={'E':-739.1643,'ZPE':0.0,'S':0.0}, sys_O={'E':-743.473,'ZPE':0.07,'S':0.0}, 
				sys_OH={'E':-748.899,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-753.566,'ZPE':0.43,'S':0.0}, sys_O2=None, 
							H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0*4)
	orr.plot_integrated('CoPC + Au ((TOP site) :: Ea={:.3}'.format(-739.1643+421.979+312.16) )
	orr.summarise_steps()
	orr.summarise_absortion()

	orr = OxigenReaction() # ========= CoPC + Au (BRIDGE site) ========= #
	orr.calculate(sys={'E':-739.336,'ZPE':0.0,'S':0.0}, sys_O={'E':-743.624,'ZPE':0.07,'S':0.0}, 
				sys_OH={'E':-749.081,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-753.727,'ZPE':0.43,'S':0.0}, sys_O2=None, 
							H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0*4)
	orr.plot_integrated('CoPC + Au (BRIDGE site) :: Ea={:.3}'.format(-739.336+421.979+312.16) )
	orr.summarise_steps()
	orr.summarise_absortion()

	orr.summarise_absortion()

	plt.show()


	# ============== FePC + Bz + Au ============== # 
	orr = OxigenReaction() # ========= FePC + 5Bz + Au (Bz site) ========= # 
	orr.calculate(sys={'E':-1121.632,'ZPE':0.0,'S':0.0}, sys_O={'E':-1127.478,'ZPE':0.07,'S':0.0},
				 sys_OH={'E':-1131.878,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-1136.452,'ZPE':0.43,'S':0.0}, sys_O2=None, 
				 			H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0)
	orr.plot_integrated('FePC + 5Bz + Au (Bz site) :: Ea={:.3}'.format(-1121.632+423.126+696.631) )
	orr.summarise_steps()
	orr.summarise_absortion()

	orr = OxigenReaction() # ========= FePC + 5Bz + Au (C site) ========= # 
	orr.calculate(sys={'E':-1121.803,'ZPE':0.0,'S':0.0}, sys_O={'E':-1127.666,'ZPE':0.07,'S':0.0},
				 sys_OH={'E':-1132.0515,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-1136.632,'ZPE':0.43,'S':0.0}, sys_O2=None, 
				 			H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0)
	orr.plot_integrated('FePC + 5Bz + Au (C site) :: Ea={:.3}'.format(-1121.803+423.126+696.631) )
	orr.summarise_steps()
	orr.summarise_absortion()
	orr = OxigenReaction() # ========= FePC + 5Bz + Au (Au site) ========= # 
	orr.calculate(sys={'E':-1121.9378,'ZPE':0.0,'S':0.0}, sys_O={'E':-1127.7830,'ZPE':0.07,'S':0.0},
				 sys_OH={'E':-1132.2013,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-1136.7592,'ZPE':0.43,'S':0.0}, sys_O2=None, 
				 			H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0)
	orr.plot_integrated('FePC + 5Bz + Au (Au site) :: Ea={:.3}'.format(-1121.9378+423.126+696.631) )
	orr.summarise_steps()
	orr.summarise_absortion()

	# ============== FePC + Au ============== # 
	orr = OxigenReaction() # ========= FePC + Au (HOLLOW site) ========= #
	orr.calculate(sys={'E':-740.5996,'ZPE':0.0,'S':0.0}, sys_O={'E':-746.2118,'ZPE':0.07,'S':0.0}, 
				  sys_OH={'E':-750.7263,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-755.3344,'ZPE':0.43,'S':0.0}, sys_O2=None, 
							H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0*4)
	orr.plot_integrated('FePC + Au (HOLLOW site) :: Ea={:.3}'.format(-740.5996+423.126+312.16) )
	orr.summarise_steps()
	orr.summarise_absortion()

	orr = OxigenReaction() # ========= FePC + Au (BRIDGE site) ========= #
	orr.calculate(sys={'E':-740.5352,'ZPE':0.0,'S':0.0}, sys_O={'E':-746.1820,'ZPE':0.07,'S':0.0}, 
				  sys_OH={'E':-750.698,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-755.3063,'ZPE':0.43,'S':0.0}, sys_O2=None, 
							H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0*4)
	orr.plot_integrated('FePC + Au (BRIDGE site) :: Ea={:.3}'.format(-740.5352+423.126+312.16) )
	orr.summarise_steps()
	orr.summarise_absortion()

	orr = OxigenReaction() # ========= FePC + Au (TOP site) ========= #
	orr.calculate(sys={'E':-740.3873,'ZPE':0.0,'S':0.0}, sys_O={'E':-746.022,'ZPE':0.07,'S':0.0}, 
				  sys_OH={'E':-750.495,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-755.306,'ZPE':0.43,'S':0.0}, sys_O2=None, 
							H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0*4)
	orr.plot_integrated( 'FePC + Au (TOP site) :: Ea={:.3}'.format(-740.387+423.126+312.16) )
	orr.summarise_steps()
	orr.summarise_absortion()

	# ============== FePC ============== # 
	orr = OxigenReaction() # ========= FePC + Free (VOID) ========= #
	orr.calculate(sys={'E':-423.2843,'ZPE':0.0,'S':0.0}, sys_O={'E':-429.238,'ZPE':0.07,'S':0.0}, 
				  sys_OH={'E':-433.647,'ZPE':0.35,'S':0.0}, sys_OOH={'E':-438.148,'ZPE':0.43,'S':0.0}, sys_O2=None, 
				 			 H2O={'E':-14.213,'ZPE':0.56,'S':0.67}, H2={'E':-6.76,'ZPE':0.27,'S':0.41}, Gb=-0.0*4)
	orr.plot_integrated('FePC + Free (VOID) :: Ec=-0.158')
	orr.summarise_steps()
	orr.summarise_absortion()

	plt.show()

######################################
# DOSCAR eg || READ and plot DOSCAR  #
######################################
def DOSCAR_plot():
	# buscamos comparar la E de interaccion M-O en terminos de las PDOS para los
	# sistemas FePC+Au y CoPC+Au. 
	doscar_Co = DOSCAR()
	doscar_Fe = DOSCAR()
	# --- CoPC --- # LOAD files
	doscar_Co.load( '/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/CoPC+Au/DOSCAR' )
	doscar_Fe.load( '/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/FePC+Au/DOSCAR' )

	# ** SIGMA ** # PLOT
	fig, ax = plt.subplots()
	doscar_Fe.plot(ion=[0], orbital=['d_z2-r2',], color=(0.8,0.3,0.3), ax=ax,
						 sum_plot=True, add_fermi=True, legend_plot=True) # sum_plot=True,
	doscar_Co.plot(ion=[0], orbital=['d_z2-r2',], color=(0.3,0.8,0.8), ax=ax,
						 sum_plot=True, add_fermi=True, legend_plot=True) # sum_plot=True,

	fig, ax = plt.subplots()
	doscar_Fe.plot(ion=[0], orbital=['d_x2-y2',], color=(0.8,0.3,0.3), ax=ax,
						 sum_plot=True, add_fermi=True, legend_plot=True) # sum_plot=True,
	doscar_Co.plot(ion=[0], orbital=['d_x2-y2',], color=(0.3,0.8,0.8), ax=ax,
						 sum_plot=True, add_fermi=True, legend_plot=True) # sum_plot=True,

	fig, ax = plt.subplots()
	doscar_Fe.plot(ion=[0], orbital=['d_xz',], color=(0.8,0.3,0.3), ax=ax,
						 sum_plot=True, add_fermi=True, legend_plot=True) # sum_plot=True,
	doscar_Co.plot(ion=[0], orbital=['d_xz',], color=(0.3,0.8,0.8), ax=ax,
						 sum_plot=True, add_fermi=True, legend_plot=True) # sum_plot=True,

	fig, ax = plt.subplots()
	doscar_Fe.plot(ion=[0], orbital=['d_xy',], color=(0.8,0.3,0.3), ax=ax,
						 sum_plot=True, add_fermi=True, legend_plot=True) # sum_plot=True,
	doscar_Co.plot(ion=[0], orbital=['d_xy',], color=(0.3,0.8,0.8), ax=ax,
						 sum_plot=True, add_fermi=True, legend_plot=True) # sum_plot=True,

	fig, ax = plt.subplots()
	doscar_Fe.plot(ion=[0], orbital=['d_yz',], color=(0.8,0.3,0.3), ax=ax,
						 sum_plot=True, add_fermi=True, legend_plot=True) # sum_plot=True,
	doscar_Co.plot(ion=[0], orbital=['d_yz',], color=(0.3,0.8,0.8), ax=ax,
						 sum_plot=True, add_fermi=True, legend_plot=True) # sum_plot=True,

	plt.show()

	# ** PI ** # PLOT
	# - Co - #
	fig = plt.figure(1)
	doscar.plot(ion=[0], orbital=['d_yz', 'd_xz',], figure=fig, color=1, 
						 sum_plot=True, add_fermi=True, legend_plot=True) # sum_plot=True,
	fig.axes[0].set_title('CoPC+O (mirando Co) d_z2-r2 (RED) / d_yz, d_xz (GREEN)')
	# ** ** SHOW
	plt.show()

DOSCAR_plot()
asd
#########################################
# BAND plot eg || READ and plot PROCAR  #
#########################################
def BAND_plot(file):
	'''
	BAND_plot(file='files/BAND/Si_diamond/rmNONSC/PROCAR')
	plot band analisys from PROCAR
	file 	:: 	STR 	:: 		Path to PROCAR
	'''
	PC = PROCAR() 
	PC.load(file)
	PC.plot()
	plt.show()
	PC.summary()

#########################################
# Read - EDIT - export POSCAR/CONTCAR   #
#########################################
def POSCAR2PDB(file):
	# Read POSCAR  -->  Set molecule to cell  -->  replicate cell  -->  SAVE as PDB and export as POSCAR
	pos = POSCAR()
	pos.load(file_name = '/home/akaris/Documents/code/VASP/v4.2/files/dataset/CoFeTPyP/Cellsize/FeTPyP/POSCAR/FePTyP_11_1295_PROPELER')

	pos.atoms = pos.get_molecule()

	pos.operations( operation={'name':'replicate', 'replicate':[1,1,1]}, save=True )
	pos.export_PDB('/home/akaris/Documents/code/Blender/VASP_BLENDER/vasp-files/FeTPyP/FePTyP_11_1295_PROPELER')
	pos.export('/home/akaris/Documents/code/Blender/VASP_BLENDER/vasp-files/FeTPyP/POSCAR')

	return pos

#########################################
# Scale relation || data analysis	    #
#########################################
def ORR_TPyP_analysis(load, save ):

	def load_features(file):
		dataset = Set()
		dataset.load_data( filename=file)
		features = dataset.extract_features(feature= { 
							'magnetization'	: {'atoms'  : {'name':['Fe', 'Co', 'Cu', 'Mg']},},
							'charge'		: {'atoms'  : {'name':['Fe', 'Co', 'Cu', 'Mg']},},
							'ORR' : [	'overpotencial_ORR_4e', 'Gabs_OOH', 'Eabs_OOH', 'Gabs_OH', 
										'Eabs_OH', 'Gabs_O', 'Eabs_O', 'G1_ORR', 'G2_ORR', 'G3_ORR', 'G4_ORR', 'limiting_step_ORR'],
							}, v=True )
		ORR_data, names = [], []
		for n in features:
			if features[n]['ORR'] != [None]:
				ORR_data.append(  features[n]['ORR'] )
				names.append( n )

		return np.array(ORR_data), names

	def ScaleRelation_plot(x, y, names, category, category_name, 
							interpolation='Linear', plot_name=True):
		N = x.shape[0]
		color = [	[0,1,1], [1,0,0], [0,0,1], [0,0,0], [0,1,0], [1,1,0], [1,0,1], [0.7,0.4,0.4], [0.3,0.8,0.3], [0.3,0.3,0.8],
					[0,1,1], [1,0,0], [0,0,1], [0,0,0], [0,1,0], [1,1,0], [1,0,1], [0.7,0.4,0.4], [0.3,0.8,0.3], [0.3,0.3,0.8]]
		xmin, xmax = np.min(x), np.max(x)
		ymin, ymax = np.min(y), np.max(y)

		# === polyfit === #
		polyfit = np.polyfit( x, y, 1, full=True)
		linear = np.array([xmin*0.95,xmax*1.05])*polyfit[0][0]+polyfit[0][1]
		print('polyfit :y = {:.4} x + {:.4} ; RMSE = {:.4}'.format(polyfit[0][0], polyfit[0][1], polyfit[1][0]))
		
		# === polyfit === #
		for n in range(0, 20):
			try:
				polyfit = np.polyfit( x[category==n], y[category==n], 1, full=True)
				linear = np.array([xmin*0.95,xmax*1.05])*polyfit[0][0]+polyfit[0][1]
				print('{} :y = {:.4} x + {:.4} ; RMSE = {:.4}'.format(category_name[n-1], polyfit[0][0], polyfit[0][1], polyfit[1][0]))
				plt.plot(np.array([xmin*0.95,xmax*1.05]), linear, color=color[n] )
				plt.plot(x[category==n], y[category==n], 'o', color=color[n], label='{} :y = {:.4} x + {:.4} ; RMSE = {:.4}'.format(category_name[n-1], polyfit[0][0], polyfit[0][1], polyfit[1][0]))
			except: pass

		plt.legend(); plt.xlabel('Gabs_OH'); plt.ylabel('Gabs_OOH')

		if plot_name:
			for i in range(N):
				for c in range(20):
					if   category[i] == c: 	plt.text(x[i], y[i], str(names[i]), bbox=dict(facecolor=color[c], alpha=0.1))			
	
	def ScaleRelation_compare_plot(x, y, z, names, category, category_name, 
							interpolation='Linear', plot_name=True):
		# category_name = ['DF286R', 'DF2', 'DF', 'opt86b', 'optPBE', 'D3-BJ', 'D3']
		fig, ax = plt.subplots(2, 3, sharex=True, sharey=True, figsize=(18, 13), dpi=80, facecolor='w', edgecolor='k') # ax.shape = [2,3]
		print(x[category==1])
		color1 = (0.95,0.4,0.4)
		color2 = (0.4,0.4,0.4)
		for i,j,n in [	[0,0,5], [0,1,0], [0,2,1], 
						[1,0,2], [1,1,3], [1,2,4], ]:

			ax[i][j].set_xlim([0.4, 1.8])
			ax[i][j].set_ylim([1, 5])
			ax[i][j].set_title( category_name[n], fontsize=20 )
			ax[i][j].grid(True, color="black", alpha=.1, linewidth=1, linestyle=":")
			ax[i][j].tick_params(axis="y",direction="in", which="major")
			ax[i][j].tick_params(axis="x",direction="in", which="major")
			ax[i][j].tick_params(axis="y",direction="in", which="minor")
			ax[i][j].tick_params(axis="x",direction="in", which="minor")
			ax[i][j].xaxis.set_major_locator(plt.MultipleLocator(0.5))
			ax[i][j].xaxis.set_minor_locator(plt.MultipleLocator(0.1))
			ax[i][j].yaxis.set_major_locator(plt.MultipleLocator(1.0))
			ax[i][j].yaxis.set_minor_locator(plt.MultipleLocator(0.2))

			xmin, xmax = 0.4, 1.8
			polyfit = np.polyfit( y[category==n+1], x[category==n+1], 1, full=True)
			linear = np.array([xmin*0.95,xmax*1.05])*polyfit[0][0]+polyfit[0][1]
			print('{} :y = {:.4} x + {:.4} ; RMSE = {:.4}'.format(category_name[n-1], polyfit[0][0], polyfit[0][1], polyfit[1][0]))
			ax[i][j].plot(np.array([xmin*0.95,xmax*1.05]), linear, color=color1, alpha=0.5)
			ax[i][j].plot(x[category==n], y[category==n], 'o', color=color1, label=r'$\Delta G_{OOH^*} =$' +'{:.3}'.format(polyfit[0][0]) + r'$\Delta G_{OH^*}$' + '+ {:.3}; NRMSE = {:.3}'.format(polyfit[0][1], polyfit[1][0]) )
			ax[i][j].plot(np.array([xmin*0.95,xmax*1.05]), np.array([xmin*0.95,xmax*1.05])*1+3.2, ':', color=color1, alpha=0.5)

			xmin, xmax = 0.4, 1.8
			polyfit = np.polyfit( y[category==n+1], z[category==n+1], 1, full=True)
			linear = np.array([xmin*0.95,xmax*1.05])*polyfit[0][0]+polyfit[0][1]
			print('{} :y = {:.4} x + {:.4} ; RMSE = {:.4}'.format(category_name[n-1], polyfit[0][0], polyfit[0][1], polyfit[1][0]))
			ax[i][j].plot(np.array([xmin*0.95,xmax*1.05]), linear, color=color2, alpha=0.5)
			ax[i][j].plot(x[category==n], y[category==n], 'o', color=color2, label=r'$\Delta G_{O^*} =$' +'{:.3}'.format(polyfit[0][0]) + r'$\Delta G_{OH^*}$' + '+ {:.3}; NRMSE = {:.3}'.format(polyfit[0][1], polyfit[1][0]) )
			ax[i][j].legend()
			ax[i][j].plot(np.array([xmin*0.95,xmax*1.05]), np.array([xmin*0.95,xmax*1.05])*2, ':',color=color2, alpha=0.5)

			m = ["o", "^", "s", "*", "D", "P"]
			for s, system in enumerate(['FeTPyP', 'CoTPyP', 'Fe*TPyPCo', 'FeTPyPCo*', 'Co*FeTPyPCo', 'CoFeTPyPCo*']):
				ax[i][j].scatter( y[category==n+1][s], x[category==n+1][s], color=color1, marker=m[s], )
				ax[i][j].scatter( y[category==n+1][s], z[category==n+1][s], color=color2, marker=m[s], )
				ax[i][j].text(y[category==n+1][s], x[category==n+1][s], names[category==n+1][s].split('_')[0], bbox=dict(facecolor=color1, alpha=0.0), )	
				ax[i][j].text(y[category==n+1][s], z[category==n+1][s], names[category==n+1][s].split('_')[0], bbox=dict(facecolor=color1, alpha=0.0), )	



		ax[1][1].set_xlabel(r'$\Delta G_{b(OH)}$', fontsize=20); ax[1][0].set_ylabel(r'$\Delta G_{b}$', fontsize=20)
		plt.savefig('test.png', dpi=300, figsize=(30, 14) )
		plt.show()

	def PCA_analysis():
		# PCA ANALYSIS #
		X = np.array([x, y]).T
		mean, std = np.mean(X,axis=0), np.std(X,axis=0)
		X = (X - mean)/std
		pca = PCA(n_components=2, svd_solver='full')
		pca.fit( X )
		print( pca.inverse_transform( [1,0] ), mean, pca.inverse_transform( [1,0] )+mean)
		v1, v2 = pca.inverse_transform( [1,0] ), pca.inverse_transform( [0,1] )
		plt.plot( [mean[0], v1[0]+mean[0]], [mean[1], v1[1]+mean[1]], lw=3, color=[0,0,0] )
		plt.plot( [mean[0], v2[0]+mean[0]], [mean[1], v2[1]+mean[1]], lw=3, color=[0,0,0] )
		PCA(n_components=2)
		print('explained_variance_ratio_', pca.explained_variance_ratio_)
		print('singular_values_', pca.singular_values_)
		plt.figure(2)
		Xt = pca.transform(X)
		for n in range(0, 20):
			try:	plt.plot(Xt[category==n, 0], Xt[category==n, 1], 'o', color=color[n], label='{} :y = {:.4} x + {:.4} ; RMSE = {:.4}'.format(category_name[n-1], polyfit[0][0], polyfit[0][1], polyfit[1][0]))
			except: pass

	if load['file'] == 'pkl':
		ORR_data, names = load_features( '{}/{}'.format(load['path'], load['name']) )

	elif load['file'] == 'npy':
		X = 	np.load('{}/{}.npy'.format(load['path'], load['name']), )
		f = 			   open('{}/names'.format(load['path'],), 'r')
		names = [ n[:-1] for n in f ]
		f.close()

	if save['state']:
		np.save('{}/{}'.format(save['path']), ORR_data)
		f = open('{}/names'.format(save['path']), 'w')
		for n in names: 	f.write(n+'\n')
		f.close()

	# ======== FILTER ======== #
	def cond_principal(name):
		if not 'Ru' in name and not 'MnPc' in name and not 'CrTPyP' in name and not 'MgPc' in name and not 'CuPc' in name :	return True 	# FeTPyP			
		else:					return False 	# TPyPCo* corr

	def cond_principal(name):
		if '13.70' in name and ('FeTPyP_' in name or 'Fe*TPyPCo_' in name or 'FeTPyPCo*_' in name or 'Co*FeTPyPCo_' in name or 
			'CoFeTPyPCo*_' in name) and not '551' in name and not '331' in name:	return True 	# FeTPyP			
		else:				return False 	# TPyPCo* corr

	# filter by value #
	x, y, z = X[:,1], X[:,3], X[:,5] 
	names = np.array(names)
	'''
	y = y[x<1.9]
	names = np.array(names)[x<1.9]
	x = x[x<1.9]

	y = y[x>0.5]
	names = np.array(names)[x>0.5]
	x = x[x>0.5]
	'''

	# filter by name #
	index_filter =  np.array([ cond_principal(n) for i, n in enumerate(names) ]) 
	x, y, z, names = x[index_filter], y[index_filter], z[index_filter], names[index_filter]
	
	# ======== INDEX ======== #
	def cond(name):
		if not 'Co' in name:	return 1	# FeTPyP			
		elif 'Fe*' in name:		return 1	# FeTPyP
		elif 'Co*Fe' in name:	return 2	# TPyPCo* sust
		elif 'Co*' in name:		return 3	# TPyPCo* corr
		else:					return 2 	# TPyPCo* sust		

	def cond_functional(name):
		if   'B86' in name:		return 1	# DF2_B86			
		elif 'DF2' in name:		return 2	# DF2
		elif 'DF' in name:		return 3	# DF
		elif 'OPT86' in name:	return 4	# OPT86
		elif 'OPT' in name:		return 5	# OPT_PBE
		elif 'D3BJ' in name:	return 6	# D3BJ	
		elif 'D3' in name:		return 7	# D3
		else: return 5

	def cond_metal(name):
		if   'Fe' in name:		return 1	# DF2_B86			
		elif 'Co' in name:		return 2	# DF2
		elif 'Cr' in name:		return 3	# DF
		elif 'Cu' in name:		return 4	# OPT86
		elif 'Mg' in name:		return 5	# OPT_PBE
		elif 'Mn' in name:		return 6	# D3BJ	
		elif 'Ni' in name:		return 7	# D3
		elif 'Rh' in name:		return 8	# D3
		elif 'Ru' in name:		return 9	# D3
		else: return 10

	def cond_metal(name):
		if   'FeTPyP' in name:		return 1	# DF2_B86			
		elif 'CoTPyP' in name:		return 2	# DF2
		elif 'CrTPyP' in name:		return 3	# DF
		elif 'CuTPyP' in name:		return 4	# OPT86
		elif 'MgTPyP' in name:		return 5	# OPT_PBE
		elif 'MnTPyP' in name:		return 6	# D3BJ	
		elif 'NiTPyP' in name:		return 7	# D3
		elif 'RhTPyP' in name:		return 8	# D3
		elif 'RuTPyP' in name:		return 9	# D3
		elif 'FePc' in name:	    return 10	# DF2_B86			
		elif 'CoPc' in name:		return 11	# DF2
		elif 'CrPc' in name:		return 12	# DF
		elif 'CuPc' in name:		return 13	# OPT86
		elif 'MgPc' in name:		return 14	# OPT_PBE
		elif 'MnPc' in name:		return 15	# D3BJ	
		elif 'NiPc' in name:		return 16	# D3
		elif 'RhPc' in name:		return 17	# D3
		elif 'RuPc' in name:		return 18	# D3
		else: return 10

	index_category = np.array([ cond_functional(n) for i, n in enumerate(names) ]) # index condicion de categoria
	category_name=[ 'FeTPyP', 'CoTPyP', 'CrTPyP', 'CuTPyP', 'MgTPyP', 'MnTPyP', 'NiTPyP', 'RhTPyP', 'RuTPyP', 
					'FePc', 'CoPc', 'CrPc', 'CuPc', 'MgPc', 'MnPc', 'NiPc', 'RhPc', 'RuPc',]
	category_name=['DF286R', 'DF2', 'DF', 'opt86b', 'optPBE', 'D3-BJ', 'D3', '', '', '', '']

	ScaleRelation_compare_plot(x, y, z, names, index_category, interpolation='Linear', category_name=category_name)
	plt.show()
	error
	ScaleRelation_plot(x, y, names, index_category, interpolation='Linear', category_name=category_name)
	
	#plt.ylim([0,7])
	#plt.xlim([-1,2.5])

