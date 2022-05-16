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
	try: from ORR import ORR
	except: print('WARNING :: Set.import_libraries() :: can not import POSCAR ')



######################################
# DOSCAR eg || READ and plot DOSCAR  #
######################################
def DOSCAR_plot():
	# buscamos comparar la E de interaccion M-O en terminos de las PDOS para los
	# sistemas FePC+Au y CoPC+Au. 
	doscar = DOSCAR()
	# --- CoPC --- # LOAD files
	doscar.load( './files/PDOS/CoPC+O+Au/DOSCAR' )

	# ** SIGMA ** # PLOT
	# - Co - #
	fig = plt.figure(1)
	doscar.plot(ion=[0], orbital=['d_z2-r2',], figure=fig, color=0, 
						 sum_plot=True, add_fermi=True, legend_plot=True) # sum_plot=True,
	# ** PI ** # PLOT
	# - Co - #
	fig = plt.figure(1)
	doscar.plot(ion=[0], orbital=['d_yz', 'd_xz',], figure=fig, color=1, 
						 sum_plot=True, add_fermi=True, legend_plot=True) # sum_plot=True,
	fig.axes[0].set_title('CoPC+O (mirando Co) d_z2-r2 (RED) / d_yz, d_xz (GREEN)')
	# ** ** SHOW
	plt.show()

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

load = {'file':'pkl', 'path':'../files/dataset/MeTPyP/DastaSet_TPyP.pkl',}
save = {'state':True, 'path':'../files/dataset/MeTPyP/ORR_only'}

# === load from pkl and save to npy === #
#ORR_TPyP_analysis(	load = {'file':'pkl', 'path':'../files/dataset/MeTPyP', 		'name':'DastaSet_TPyP.pkl'},
#					save = {'state':True, 'path':'../files/dataset/MeTPyP/ORR_only','name':'ORR_data'} )
# === load from npy === #
ORR_TPyP_analysis(	load = {'file':'npy', 'path':'../files/dataset/MeTPyP/ORR_only', 		'name':'ORR_data'},
					save = {'state':False,} )
error
# === load from npy === #
ORR_TPyP_analysis(	load = {'file':'npy', 'path':'../files/dataset/Metales', 		'name':'NNy_metales'},
					save = {'state':False,} )
error


#########################################
# SET plk load || load 				    #
#########################################
def load_features(file):
	dataset = Set()
	dataset.load_data( filename=file)
	save_path = '/home/akaris/Documents/code/VASP/v3.0/files/ORR/summary070921'
	#dataset.summary(xml=True, path=save_path)
	#error
	#dataset.save_plot_orr( path=save_path )
	features = dataset.extract_features(feature= { 
						#'PDOS': {	'config' : {'start':-10.0, 'end':5.0, 'point':1000},
						#			'atoms'  : {'name':['Fe', 'Co', 'Cu', 'Mg']},  
						#			'orbital': [9,10,11,12,13,14,15,16,17,18,]  }, # [9,10,11,12,13,14,15,16,17,18]
						'magnetization'	: {'atoms'  : {'name':['Fe', 'Co', 'Cu', 'Mg']},},
						'charge'		: {'atoms'  : {'name':['Fe', 'Co', 'Cu', 'Mg']},},
						'ORR' : [	'overpotencial_ORR_4e', 'Gabs_OOH', 'Eabs_OOH', 'Gabs_OH', 
									'Eabs_OH', 'Gabs_O', 'Eabs_O', 'G1_ORR', 'G2_ORR', 'G3_ORR', 'G4_ORR', 'limiting_step_ORR'],
						}, v=True )

	NNx = []
	NNy = []
	NNz = []
	names = []
	print(features)
	for n in features:
		for m in features[n]:
			NNy.append(  features[n]['ORR'] )
			names.append( n )
			if 'PDOS12312' in m:
				try:
					if type(features[n]['ORR'][0]) != type(None) and type(np.reshape(features[n][m]['PDOS'], (1000*10))) != type(None):
						NNx.append(  np.concatenate((features[n][m]['magnetization'][0], features[n][m]['charge'][0][0,:], np.reshape(features[n][m]['PDOS'], (1000*10)))) )
						NNy.append( features[n]['ORR'] )
						names.append( n )
				except: pass


	#random.shuffle(NNx)
	#random.shuffle(NNy)
	return NNx, NNy, NNz, names


# ========== set path ========== #
path_TPyP = 	'../files/dataset/MeTPyP/DastaSet_TPyP.pkl'
path_metales = 	'../files/dataset/Metales/metales.pkl'
path_PC = 		'../files/dataset/MePC/M1PC.pkl'

'''
# ========== Load pkl selected features and save npy data ========== #
NNx, NNy, NNz, names = load_features(path_metales)
np.save('../files/dataset/Metales/NNx_metales', NNx)
np.save('../files/dataset/Metales/NNy_metales', NNy)
np.save('../files/dataset/Metales/NNz_metales', NNz)
f = open('../files/dataset/Metales/names_metales', 'w')
for n in names: 	f.write(n+'\n')
f.close()

NNx, NNy, NNz, names = load_features(path_TPyP)
np.save('../files/dataset/MeTPyP/NNx_TPyP', NNx)
np.save('../files/dataset/MeTPyP/NNy_TPyP', NNy)
np.save('../files/dataset/MeTPyP/NNz_TPyP', NNz)
f = open('../files/dataset/MeTPyP/names_TPyP', 'w')
for n in names: 	f.write(n+'\n')
f.close()

NNx, NNy, NNz, names = load_features(path_PC)
np.save('../files/dataset/MePC/NNx_PC', NNx)
np.save('../files/dataset/MePC/NNy_PC', NNy)
np.save('../files/dataset/MePC/NNz_PC', NNz)
f = open('../files/dataset/MePC/names_PC', 'w')
for n in names: 	f.write(n+'\n')
f.close()
'''


# ========== Load npy data ========== #
NNx_metales = 	np.load('../files/dataset/Metales/NNx_metales.npy', )
NNy_metales = 	np.load('../files/dataset/Metales/NNy_metales.npy', )
NNz_metales = 	np.load('../files/dataset/Metales/NNz_metales.npy', )
f = 			   open('../files/dataset/Metales/names_metales', 'r')
names_metales = [ n[:-1] for n in f ]
f.close()

'''
NNx_TPyP = 		np.load('../files/dataset/MeTPyP/NNx_TPyP.npy', )
NNy_TPyP = 		np.load('../files/dataset/MeTPyP/NNy_TPyP.npy', )
NNz_TPyP = 		np.load('../files/dataset/MeTPyP/NNz_TPyP.npy', )
f = 			   open('../files/dataset/MeTPyP/names_TPyP', 'r')
names_TPyP = [ n[:-1] for n in f ]
f.close()

NNx_PC = 		np.load('../files/dataset/MePC/NNx_PC.npy', )
NNy_PC =		np.load('../files/dataset/MePC/NNy_PC.npy', )
NNz_PC =		np.load('../files/dataset/MePC/NNz_PC.npy', )
f = 			   open('../files/dataset/MePC/names_PC', 'r')
names_PC = [ n[:-1] for n in f ]
f.close()
'''

# ========== Reshape data ========== #
#NNx = np.concatenate( (NNx_TPyP, NNx_PC, NNx_metales ) )
#NNy = np.concatenate( (NNy_TPyP, NNy_PC, NNy_metales ) )
NNx = NNx_metales
NNy = NNy_metales

NN_names = np.array(names_metales)

print(NNy.shape, NNx.shape, NN_names.shape)
NNx = NNx[:,4000:6000]
#NNx = NNx[:,4000:5000] - NNx[:,5000:6000]
NNy = NNy[:,:2]

NN_names = NN_names[ NNy[:,0]<2 ]
NNx = NNx[ NNy[:,0]<2,: ]
NNy = NNy[ NNy[:,0]<2,: ]

'''
# ========== Manifold ========== #
n_neighbors = 3
n_components = 2
LLE = manifold.LocallyLinearEmbedding(n_neighbors, n_components, eigen_solver='auto',method='modified')
LLE = manifold.TSNE(n_components=n_components, init='pca', random_state=0)
Y = LLE.fit_transform(NNx)
colores = [ [1,0,0] if 'Fe' in n else [0,0,1] if 'Co' in n  else  [0,0,1] if 'Ni' in n else [0,1,0] for n in NN_names ]
plt.scatter(Y[:, 0], Y[:, 1], color=colores, cmap=plt.cm.Spectral)
plt.show()
'''

# ========== shuffle split ========== #
#X_train, X_test, y_train, y_test = train_test_split(NNx, NNy, test_size=0.05)
sss = ShuffleSplit(n_splits=1, test_size=0.2)
sss.get_n_splits(NNx, NNy)
train_index, test_index = next(sss.split(NNx, NNy)) 

X_train, X_test = NNx[train_index], NNx[test_index] 
y_train, y_test = NNy[train_index], NNy[test_index]

# ========== SELECT DATA ========== #
#def cond(n): return 'Fe' in n or 'Ni' in n or 'Co' in n   or 'Cu' in n   or 'Mn' in n  or 'Cr' in n 
def cond(n):  return 'Co' in n or 'Ni' in n  or 'Cu' in n  or 'Mn' in n or 'Fe' in n or 'Rh' in n or 'Fe' in n 
def cond(n):  return ('TPyP' in n and 'Ni' in n) or ('TPyP' in n and 'Mn' in n) or ('TPyP' in n and 'Cr' in n) or ('TPyP' in n and 'Cu' in n)  or ('Ru' in n and 'Pc' in n)
def cond(n):  return ('TPyP' in n and 'Ni' in n) or ('TPyP' in n and 'Co' in n) or ('TPyP' in n and 'Cu' in n) 

def cond_train(n):  return  not ('Pc' in n  or 'Fe' in n or 'Ni' in n or 'Mg' in n) 
def cond_test(n):  return  'Fe' in n or 'Ni' in n


train_index = list(filter(lambda x: x>=0, [ i if cond_train(n) else -1 for i, n in enumerate(NN_names)]     ))
test_index = list(filter(lambda x: x>=0, [ i if cond_test(n) else -1 for i, n in enumerate(NN_names)]  ))

X_train, X_test = NNx[train_index], NNx[test_index] 
y_train, y_test = NNy[train_index], NNy[test_index]
print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)

# data mining add reverse upside dpwm PDOS
# ========== Train NN ========== #
regr = MLPRegressor(random_state=3, activation='relu', #early_stopping=True, validation_fraction=0.10,
					hidden_layer_sizes=(1000, 1000),max_iter=1000, #learning_rate='adaptive', warm_start=True,
					#hidden_layer_sizes=(  400, 540, 540, 400,  ),max_iter=5000,
					momentum=0.9, verbose=True, n_iter_no_change=20, ).fit(X_train, y_train)


# ========== NN regresion ========== #
NNp_test =  regr.predict(X_test)
NNp =  regr.predict(X_train)
NNp_all =  regr.predict(NNx)

text = ['overpotencial_ORR_4e', 'Gabs_OOH', 'Eabs_OOH', 'Gabs_OH', 'Eabs_OH', 
'Gabs_O', 'Eabs_O', 'G1_ORR', 'G2_ORR', 'G3_ORR', 'G4_ORR', 'limiting_step_ORR']

# ========== Regresion plot ========== #
for n in range(1):
	plt.figure(n)
	for i, m in enumerate(NNp_all[:,n]):		
		if abs(NNy[i,n] - NNp_all[i,n]) > 0.0:
			plt.text(NNy[i,n], NNp_all[i,n], NN_names[i], bbox=dict(facecolor='red', alpha=0.04))

	plt.title( text[n] )
	plt.plot(y_train[:,n], NNp[:,n], 'o')
	plt.plot(y_test[:,n], NNp_test[:,n], 'o')


	plt.plot( [0,10], [0,10] )
	plt.plot( [0,10], [0.5,10.5] )
	plt.plot( [0,10], [-0.5,9.5] )
	if n == 1 or n ==2:
		plt.xlim(3, 7)
		plt.ylim(3, 7)
	else:
		plt.xlim(0, 4)
		plt.ylim(0, 4)

'''
plt.figure(10)
plt.plot(y_train[:,n], np.max(NNp[:,7:12], axis=1) , 'o')
plt.plot(y_test[:,n],  np.max(NNp_test[:,7:12], axis=1) , 'o')

plt.figure(11)
plt.hist(NNy[:,0] - NNp_all[:,0], bins=32)


limiting_step_ORR_ytrain =  np.argmax( NNy[:,7:12], axis=1 ) 
limiting_step_ORR_NNp =  np.argmax( NNp_all[:,7:12], axis=1 ) 

print(' ***** TRAIN *****')
for i, m in enumerate(train_index): print('{} : {} : {:.3} : {} : {} '.format(limiting_step_ORR_ytrain[m], limiting_step_ORR_NNp[m],  NNy[m,0], NNp_all[m,0], NN_names[m],) )
print(' ***** TEST *****')
for i, m in enumerate(test_index): print('{} : {} : {:.3} : {:.3} : {} '.format(limiting_step_ORR_ytrain[m], limiting_step_ORR_NNp[m],  NNy[m,0], NNp_all[m,0], NN_names[m],) )
'''
print( np.sum(np.abs(y_test[:,0]-NNp_test[:,0])) )
print( np.sum(np.abs(y_test[:,0]-NNp_test[:,0])) )
plt.show()

error
'''
print('Total number of samples:{}'.format(limiting_step_ORR_NNp.shape[0] ))
c = np.sum(limiting_step_ORR_ytrain == limiting_step_ORR_NNp) 
print('Confusion vector : ')

print('correctly\t\t\t{:}\t\t\tIncorrectly\t\t\t{:}'.format(c, limiting_step_ORR_NNp.shape[0]-c) )
print('correctly\t\t\t{:.3}\t\t\tIncorrectly\t\t\t{:.3}'.format(c/limiting_step_ORR_NNp.shape[0], 
											(limiting_step_ORR_NNp.shape[0]-c)/limiting_step_ORR_NNp.shape[0]) )

# ------------------------------------- # # ------------------------------------- #
limiting_step_ORR_ytrain =  np.argmax( y_test[:,7:12], axis=1 ) 
limiting_step_ORR_NNp =  np.argmax( NNp_test[:,7:12], axis=1 ) 

print('Total number of samples:{}'.format(limiting_step_ORR_NNp.shape[0] ))
c = np.sum(limiting_step_ORR_ytrain == limiting_step_ORR_NNp) 
print('Confusion vector : ')

print('correctly\t\t\t{:}\t\t\tIncorrectly\t\t\t{:}'.format(c, limiting_step_ORR_NNp.shape[0]-c) )
print('correctly\t\t\t{:.3}\t\t\tIncorrectly\t\t\t{:.3}'.format(c/limiting_step_ORR_NNp.shape[0], 
											(limiting_step_ORR_NNp.shape[0]-c)/limiting_step_ORR_NNp.shape[0]) )
'''

plt.show()

error
CoTPyP 		= [ 0.57, 	 0.53, 	 0.16, 	 0.22, 	 0.21, 	 0.17, 	 0.20, 	]
FeTPyP      = [  0.50, 	 0.51, 	 0.67, 	 0.69, 	 0.52, 	 0.53, 	 0.59, 	]
CoFeTPyP    = [  0.83 ,	 0.84, 	 0.35,   0.25, 	 0.53,	 0.54, 	 0.35,	]
Fe0TPyPCo   = [ 0.59, 	 0.60, 	 0.77, 	 0.75, 	 0.64, 	 0.65, 	 0.68,  ]

FeTPyPCo0   = [  0.30 ,	 0.30,	 0.60, 	 0.70 ,	 0.47, 	 0.55, 	 0.56, 	]
Co0FeTPyPCo = [  0.86 ,	 0.86, 	 0.40,   0.28, 	 0.58, 	 0.52, 	 0.41, 	]
CoFeTPyPCo0 = [  0.36, 	 0.35, 	 0.63, 	 0.72 ,	 0.46, 	 0.46, 	 0.59, 	]

plotlist = [CoTPyP, FeTPyP, CoFeTPyP, Fe0TPyPCo, FeTPyPCo0, Co0FeTPyPCo, CoFeTPyPCo0 ]
plotname = ['CoTPyP', 'FeTPyP', 'CoFeTPyP', 'Fe0TPyPCo', 'FeTPyPCo0', 'Co0FeTPyPCo', 'CoFeTPyPCo0']
fig, ax = plt.subplots()
for i, n in enumerate(plotlist):
	order = [1,0,4,5,6,2,3]
	ax.plot( np.array(n)[order], '-o', label=plotname[i] )
legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')

ax.set_xticklabels(['','D3bj', 'D3', 'DF2B86', 'OPTB86', 'OPTPBE', 'DF', 'DF2'])

plt.show()
error
