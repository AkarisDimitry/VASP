#import yahoo_fin.stock_info as yf
#import requests, ftplib, io, re, json, datetime, time, Logspp
import pandas as pd
import random, Logs
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression

# -------------------------------------------------------------------------------- #

# -------------------------------------------------------------------------------- #

class QMlearning(object):
	def __init__(self, name:str=None,
						X:np.ndarray=None, Y:np.ndarray=None, Z:np.ndarray=None, 
						N:np.ndarray=None, ):
		self._name = name

		self._X = X
		self._Y = Y
		self._Z = Z

		self._N = N
		# --- NN --- #
		self.act_func = 'tanh'
		self.neural_regresor_model = None

		self.cnames = {
'aliceblue':            '#F0F8FF','antiquewhite':         '#FAEBD7','aqua':                 '#00FFFF','aquamarine':           '#7FFFD4',
'azure':                '#F0FFFF','beige':                '#F5F5DC','bisque':               '#FFE4C4','black':                '#000000',
'blanchedalmond':       '#FFEBCD','blue':                 '#0000FF','blueviolet':           '#8A2BE2','brown':                '#A52A2A',
'burlywood':            '#DEB887','cadetblue':            '#5F9EA0','chartreuse':           '#7FFF00','chocolate':            '#D2691E',
'coral':                '#FF7F50','cornflowerblue':       '#6495ED','cornsilk':             '#FFF8DC','crimson':              '#DC143C',
'cyan':                 '#00FFFF','darkblue':             '#00008B','darkcyan':             '#008B8B','darkgoldenrod':        '#B8860B',
'darkgray':             '#A9A9A9','darkgreen':            '#006400','darkkhaki':            '#BDB76B','darkmagenta':          '#8B008B',
'darkolivegreen':       '#556B2F','darkorange':           '#FF8C00','darkorchid':           '#9932CC','darkred':              '#8B0000',
'darksalmon':           '#E9967A','darkseagreen':         '#8FBC8F','darkslateblue':        '#483D8B','darkslategray':        '#2F4F4F',
'darkturquoise':        '#00CED1','darkviolet':           '#9400D3','deeppink':             '#FF1493','deepskyblue':          '#00BFFF',
'dimgray':              '#696969','dodgerblue':           '#1E90FF','firebrick':            '#B22222','floralwhite':          '#FFFAF0',
'forestgreen':          '#228B22','fuchsia':              '#FF00FF','gainsboro':            '#DCDCDC','ghostwhite':           '#F8F8FF',
'gold':                 '#FFD700','goldenrod':            '#DAA520','gray':                 '#808080','green':                '#008000',
'greenyellow':          '#ADFF2F','honeydew':             '#F0FFF0','hotpink':              '#FF69B4','indianred':            '#CD5C5C',
'indigo':               '#4B0082','ivory':                '#FFFFF0','khaki':                '#F0E68C','lavender':             '#E6E6FA',
'lavenderblush':        '#FFF0F5','lawngreen':            '#7CFC00','lemonchiffon':         '#FFFACD','lightblue':            '#ADD8E6',
'lightcoral':           '#F08080','lightcyan':            '#E0FFFF','lightgoldenrodyellow': '#FAFAD2','lightgreen':           '#90EE90',
'lightgray':            '#D3D3D3','lightpink':            '#FFB6C1','lightsalmon':          '#FFA07A','lightseagreen':        '#20B2AA',
'lightskyblue':         '#87CEFA','lightslategray':       '#778899','lightsteelblue':       '#B0C4DE','lightyellow':          '#FFFFE0',
'lime':                 '#00FF00','limegreen':            '#32CD32','linen':                '#FAF0E6','magenta':              '#FF00FF',
'maroon':               '#800000','mediumaquamarine':     '#66CDAA','mediumblue':           '#0000CD','mediumorchid':         '#BA55D3',
'mediumpurple':         '#9370DB','mediumseagreen':       '#3CB371','mediumslateblue':      '#7B68EE','mediumspringgreen':    '#00FA9A',
'mediumturquoise':      '#48D1CC','mediumvioletred':      '#C71585','midnightblue':         '#191970','mintcream':            '#F5FFFA',
'mistyrose':            '#FFE4E1','moccasin':             '#FFE4B5','navajowhite':          '#FFDEAD','navy':                 '#000080',
'oldlace':              '#FDF5E6','olive':                '#808000','olivedrab':            '#6B8E23','orange':               '#FFA500',
'orangered':            '#FF4500','orchid':               '#DA70D6','palegoldenrod':        '#EEE8AA','palegreen':            '#98FB98',
'paleturquoise':        '#AFEEEE','palevioletred':        '#DB7093','papayawhip':           '#FFEFD5','peachpuff':            '#FFDAB9',
'peru':                 '#CD853F','pink':                 '#FFC0CB','plum':                 '#DDA0DD','powderblue':           '#B0E0E6',
'purple':               '#800080','red':                  '#FF0000','rosybrown':            '#BC8F8F','royalblue':            '#4169E1',
'saddlebrown':          '#8B4513','salmon':               '#FA8072','sandybrown':           '#FAA460','seagreen':             '#2E8B57',
'seashell':             '#FFF5EE','sienna':               '#A0522D','silver':               '#C0C0C0','skyblue':              '#87CEEB',
'slateblue':            '#6A5ACD','slategray':            '#708090','snow':                 '#FFFAFA','springgreen':          '#00FF7F',
'steelblue':            '#4682B4','tan':                  '#D2B48C','teal':                 '#008080','thistle':              '#D8BFD8',
'tomato':               '#FF6347','turquoise':            '#40E0D0','violet':               '#EE82EE','wheat':                '#F5DEB3',
'white':                '#FFFFFF','whitesmoke':           '#F5F5F5','yellow':               '#FFFF00','yellowgreen':          '#9ACD32'}

		self.functional_metadata = {
			'D3_' 		: {'marker'	: 'o'},
			'D3BJ' 		: {'marker'	: 'v'},
			'DF_' 		: {'marker'	: '1'},
			'DF2_' 		: {'marker'	: 's'},
			'DF2B86' 	: {'marker'	: '+'},
			'OPT86' 	: {'marker'	: 'P'},
			'OPTPBE' 	: {'marker'	: '*'},
		}

		self.atomic_metadata = {
		    # Name [period, column,                  valence]
		    'Al' : [    3,      13,      [3]                , 0     ,	'aliceblue'			],
		    'Bi' : [    6,      15,      [3,5]              , 1     ,	'azure'				],
		    'Ca' : [    4,       2,      [2]                , 2     ,	'blanchedalmond'	],
		    'Cd' : [    5,      12,      [2]                , 3     ,	'coral'				],
		    'Co' : [    4,       9,      [2,3]              , 4     ,	'cyan'				],
		    'Cr' : [    4,       6,      [2,3,6]            , 5     ,	'darkgray'			],
		    'Cu' : [    4,      11,      [1,2]              , 6     ,	'darkolivegreen'	],
		    'Ir' : [    6,       9,      [1,2]              , 7     ,	'darkturquoise'		],
		    'Fe' : [    4,       8,      [2,3]              , 8     ,	'dimgray'			],
		    'Mg' : [    3,       2,      [2]                , 9     ,	'indigo'			],
		    'Mn' : [    4,       7,      [2,3,4,6,7]        , 10    ,	'lavenderblush'		],
		    'Mo' : [    5,       6,      [2,3,4,5,6]        , 11    ,	'lightskyblue'		],
		    'Ni' : [    4,      10,      [2,3]              , 12    ,	'maroon'			],
		    'Pb' : [    6,      14,      [2,4]              , 13    ,	'mediumturquoise'	],
		    'Pd' : [    5,      10,      [2,4]              , 14    ,	'skyblue'			],
		    'Pt' : [    6,      10,      [2,4]              , 15    ,	'orangered'			],
		    'Ru' : [    5,       8,      [2,3,4,5,6,7,8]    , 16    ,	'paleturquoise'		],
		    'Rh' : [    5,       9,      [2,3,4,5,6]        , 17    ,	'purple'			],
		    'Sc' : [    4,       3,      [3]                , 18    ,	'saddlebrown'		],
		    'Ti' : [    4,       4,      [2,3,4]            , 19    ,	'plum'				],
		    'Tc' : [    5,       7,      [2,3,4]            , 20    ,	'tomato'			],
		     'V' : [    4,       5,      [2,3,4,5]          , 21    ,	'teal'				],
		    'Zn' : [    4,      12,      [2,3,4,5]          , 22    ,	'slategray'			],
		}
		for key, item in self.atomic_metadata.items():	item.append( random.choice(list(self.cnames)) )

	# ------------------- name ------------------- #
	@property
	def name(self, ) -> str:
		if self._name is None:	return 'Unknow'	
		else:					return self._name

	@name.setter
	def name(self, name:str):	self._name = name

	@name.deleter
	def name(self, ):	del self._name

	# ------------------- X ------------------- #
	@property
	def X(self, ) -> np.ndarray:
		if self._X is None:	return 'Unknow'	
		else:					return self._X

	@X.setter
	def name(self, X:np.ndarray):	self._X = X

	@X.deleter
	def X(self, ):	del self._X

	# ------------------- X ------------------- #
	@property
	def Y(self, ) -> np.ndarray:
		if self._Y is None:	return 'Unknow'	
		else:					return self._Y

	@Y.setter
	def name(self, Y:np.ndarray):	self._Y = Y

	@Y.deleter
	def Y(self, ):	del self._Y


	# ------------------- X ------------------- #
	@property
	def Z(self, ) -> np.ndarray:
		if self._Z is None:	return 'Unknow'	
		else:					return self._Z

	@Z.setter
	def name(self, Z:np.ndarray):	self._Z = Z

	@Z.deleter
	def Z(self, ):	del self._Z


	# ------------------- N ------------------- #
	@property
	def N(self, ) -> np.ndarray:
		if self._N is None:	return 'Unknow'	
		else:					return self._N

	@N.setter
	def name(self, N:np.ndarray):	self._N = N

	@N.deleter
	def N(self, ):	del self._N

	def load(self, file:str=None) -> np.ndarray:
		if file[-3:] == '.npy':	return np.load(file)
		else:					
			try:	return np.loadtxt(file)
			except:	return np.loadtxt(file, dtype=str)

	@Logs.LogDecorator()
	def load_folder(self, path:str=None, 
						save:bool=True, v:bool=True) -> list:

		try:
			from os import listdir
			from os.path import isfile, join
		except:
			print('>> (!) ERROR :: can not import os :: try pip3 install os')

		loaded_dict = {}
		for n in [f for f in listdir(path) if isfile(join(path, f))]:

			if n[-5:] == 'X.dat':
				if v: print( f' >> Loading X.dat :: from {path}/{n}')
				loaded_dict['X'] = self.load(file=f'{path}/{n}')
				if v: print( f' (*) Loaded X.dat :: {loaded_dict["X"].shape}')
				if save: 
					if self._X is None: self._X = loaded_dict['X']
					else:				self._X = np.concatenate( (self._X, loaded_dict['X']), axis=0 )

			if n[-5:] == 'Y.dat':
				if v: print( f' >> Loading Y.dat :: from {path}/{n}')
				loaded_dict['Y'] = self.load(file=f'{path}/{n}')
				if v: print( f' (*) Loaded Y.dat :: {loaded_dict["Y"].shape}')
				if save: 
					if self._Y is None: self._Y = loaded_dict['Y']
					else:				self._Y = np.concatenate( (self._Y, loaded_dict['Y']), axis=0 )
			if n[-5:] == 'Z.dat':
				if v: print( f' >> Loading Z.dat :: from {path}/{n}')
				loaded_dict['Z'] = self.load(file=f'{path}/{n}')
				if v: print( f' (*) Loaded Z.dat :: {loaded_dict["Z"].shape}')
				if save: 
					if self._Z is None: self._Z = loaded_dict['Z']
					else:				self._Z = np.concatenate( (self._Z, loaded_dict['Z']), axis=0 )

			if n[-5:] == 'N.dat':
				if v: print( f' >> Loading N.dat :: from {path}/{n}')
				loaded_dict['N'] = self.load(file=f'{path}/{n}')
				if v: print( f' (*) Loaded N.dat :: {loaded_dict["N"].shape}')
				if save: 
					if self._N is None: self._N = loaded_dict['N']
					else:				self._N = np.concatenate( (self._N, loaded_dict['N']) )
		return loaded_dict

	def filter(self, filters:list=None, X:np.ndarray=None, Y:np.ndarray=None, Z:np.ndarray=None, N:np.ndarray=None,
					v:bool=True, save:bool=True):
		X = X if not X is None else self.X
		Y = Y if not Y is None else self.Y
		Z = Z if not Z is None else self.Z
		N = N if not N is None else self.N

		for f in filters:
			if f['type'] == 'N':
				if f['condition'] == 'whitelist':
					filter_array = np.array([ np.any( np.array([ m in n for m in f['list']],dtype=bool) ) == True for n in N])

		if save:
			self._X = X[filter_array]
			try:	self._X = X[filter_array]
			except: pass

			try:	self._Y = Y[filter_array]
			except: pass

			try:	self._Z = Z[filter_array]
			except: pass

			try:	self._N = N[filter_array]
			except: pass

		return True

	def ORR_analysis(self, X:np.ndarray=None, Y:np.ndarray=None, Z:np.ndarray=None, N:np.ndarray=None,
						curve_fit:bool=True,
						degree:int=1, label:list=None, save:bool=True, v:bool=True, text:bool=False):
		X = X if not X is None else self.X
		Y = Y if not Y is None else self.Y
		Z = Z if not Z is None else self.Z
		N = N if not N is None else self.N
		label = label if not label is None else [ 	'overpotencial_ORR_4e', 
													'Gabs_OOH', 'Eabs_OOH', 
													'Gabs_OH',  'Eabs_OH', 
													'Gabs_O', 'Eabs_O', 
													'G1_ORR', 'G2_ORR', 'G3_ORR', 'G4_ORR', 
													'limiting_step_ORR']

		def band_center(band, Emin=0, Emax=1 ):  return np.sum(band*np.linspace(Emin, Emax, num=band.shape[0]))/np.sum(band)
		samples, Xdim = X.shape
		samples, Ydim = Y.shape

		fig, ([ax11, ax12, ax13], [ax21, ax22, ax23], [ax31, ax32, ax33]) = plt.subplots(3, 3, figsize=(10, 10) )
		fig.tight_layout()
		ax11.set(xlabel='Gabs_OOH (eV)'	, ylabel='overpotencial_ORR_4e  (eV)'	)
		ax12.set(xlabel='Gabs_OH (eV)'	, ylabel='overpotencial_ORR_4e  (eV)'	)
		ax13.set(xlabel='Gabs_O (eV)'	, ylabel='overpotencial_ORR_4e  (eV)'	)
		ax21.set(xlabel='G1_ORR (eV)'	, ylabel='overpotencial_ORR_4e  (eV)'	)
		ax22.set(xlabel='G2_ORR (eV)'	, ylabel='overpotencial_ORR_4e  (eV)'	)
		ax23.set(xlabel='G3_ORR (eV)'	, ylabel='overpotencial_ORR_4e  (eV)'	)
		ax31.set(xlabel='Gabs_OH (eV)'	, ylabel='Gabs_O  (eV)'					)
		ax32.set(xlabel='Gabs_OOH (eV)'	, ylabel='Gabs_OH  (eV)'				)
		ax33.set(xlabel='Gabs_OOH (eV)'	, ylabel='Gabs_O  (eV)'					)

		for n in range(samples):
			for key, metadata in self.atomic_metadata.items():
				if key in N[n]: break

			for key, marker in self.functional_metadata.items():
				if key in N[n]: marker=marker['marker']; break

			ax11.plot( X[n,label.index('Gabs_OOH')], X[n,label.index('overpotencial_ORR_4e')], marker=marker , color=metadata[4] )
			if text: ax11.text( X[n,label.index('Gabs_OOH')], X[n,label.index('overpotencial_ORR_4e')], N[n], color=metadata[4] )
			
			ax12.plot( X[n,label.index('Gabs_OH')], X[n,label.index('overpotencial_ORR_4e')], marker=marker , color=metadata[4] )
			if text: ax12.text( X[n,label.index('Gabs_OH')], X[n,label.index('overpotencial_ORR_4e')], N[n], color=metadata[4] )

			ax13.plot( X[n,label.index('Gabs_O')], X[n,label.index('overpotencial_ORR_4e')], marker=marker , color=metadata[4] )
			if text: ax13.text( X[n,label.index('Gabs_O')], X[n,label.index('overpotencial_ORR_4e')], N[n], color=metadata[4] )


			ax21.plot( X[n,label.index('G1_ORR')], X[n,label.index('overpotencial_ORR_4e')], marker=marker , color=metadata[4] )
			if text: ax21.text( X[n,label.index('G1_ORR')], X[n,label.index('overpotencial_ORR_4e')], N[n], color=metadata[4] )

			ax22.plot( X[n,label.index('G2_ORR')], X[n,label.index('overpotencial_ORR_4e')], marker=marker , color=metadata[4] )
			if text: ax22.text( X[n,label.index('G2_ORR')], X[n,label.index('overpotencial_ORR_4e')], N[n], color=metadata[4] )

			ax23.plot( X[n,label.index('G3_ORR')], X[n,label.index('overpotencial_ORR_4e')], marker=marker , color=metadata[4] )
			if text: ax23.text( X[n,label.index('G3_ORR')], X[n,label.index('overpotencial_ORR_4e')], N[n], color=metadata[4] )


			ax31.plot( X[n,label.index('Gabs_OH')], X[n,label.index('Gabs_O')], marker=marker , color=metadata[4] )
			if text: ax31.text( X[n,label.index('Gabs_OH')], X[n,label.index('Gabs_O')], N[n], color=metadata[4] )

			ax32.plot( X[n,label.index('Gabs_OOH')], X[n,label.index('Gabs_OH')], marker=marker , color=metadata[4] )
			if text: ax32.text( X[n,label.index('Gabs_OOH')], X[n,label.index('Gabs_OH')], N[n], color=metadata[4] )

			ax33.plot( X[n,label.index('Gabs_OOH')], X[n,label.index('Gabs_O')], marker=marker , color=metadata[4] )
			if text: ax33.text( X[n,label.index('Gabs_OOH')], X[n,label.index('Gabs_O')], N[n], color=metadata[4] )

		if curve_fit:
			from scipy.optimize import curve_fit

			def modelo(x, a, b, c, d):
				Y = np.zeros_like(x)
				Y[x<=c] = a*x[x<=c] + b 
				Y[x> c] = d*x[x> c] + (a-d)*c+b
				return Y

			Xmin, Xmax = np.min(X[:,label.index('Gabs_OOH')])*0.9, np.max(X[:,label.index('Gabs_OOH')])*1.1
			Xnew = np.linspace(Xmin, Xmax, 100)
			popt, pcov = curve_fit(modelo, X[:,label.index('Gabs_OOH')], X[:,label.index('overpotencial_ORR_4e')], p0=[1, 5, X[np.argmin(X[:,label.index('overpotencial_ORR_4e')]),label.index('Gabs_OOH')], -1])
			ss_res = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - modelo(X[:,label.index('Gabs_OOH')], *popt))**2  )	# Suma de los cuadrados de los residuos
			ss_tot = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - np.mean(X[:,label.index('overpotencial_ORR_4e')]) )**2  ) 		# Suma total de cuadrados
			ax11.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}Gabs_OOH+{popt[1]:.2} \n {popt[3]:.2}Gabs_OOH+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
			ax11.legend(loc='best')

			Xmin, Xmax = np.min(X[:,label.index('Gabs_OH')])*0.9, np.max(X[:,label.index('Gabs_OH')])*1.1
			Xnew = np.linspace(Xmin, Xmax, 100)
			popt, pcov = curve_fit(modelo, X[:,label.index('Gabs_OH')], X[:,label.index('overpotencial_ORR_4e')], p0=[-1,8, X[np.argmin(X[:,label.index('overpotencial_ORR_4e')]),label.index('Gabs_OH')] ,1])
			ss_res = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - modelo(X[:,label.index('Gabs_OH')], *popt))**2  )	# Suma de los cuadrados de los residuos
			ss_tot = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - np.mean(X[:,label.index('overpotencial_ORR_4e')]) )**2  ) 		# Suma total de cuadrados
			ax12.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}Gabs_OH+{popt[1]:.2} \n {popt[3]:.2}Gabs_OH+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
			ax12.legend(loc='best')

			Xmin, Xmax = np.min(X[:,label.index('Gabs_O')])*0.9, np.max(X[:,label.index('Gabs_O')])*1.1
			Xnew = np.linspace(Xmin, Xmax, 100)
			popt, pcov = curve_fit(modelo, X[:,label.index('Gabs_O')], X[:,label.index('overpotencial_ORR_4e')], p0=[-1,8, X[np.argmin(X[:,label.index('overpotencial_ORR_4e')]),label.index('Gabs_O')] ,1])
			ss_res = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - modelo(X[:,label.index('Gabs_O')], *popt))**2  )	# Suma de los cuadrados de los residuos
			ss_tot = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - np.mean(X[:,label.index('overpotencial_ORR_4e')]) )**2  ) 		# Suma total de cuadrados
			ax13.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}Gabs_O+{popt[1]:.2} \n {popt[3]:.2}Gabs_O+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
			ax13.legend(loc='best')


			Xmin, Xmax = np.min(X[:,label.index('G1_ORR')])*0.9, np.max(X[:,label.index('G1_ORR')])*1.1
			Xnew = np.linspace(Xmin, Xmax, 100)
			popt, pcov = curve_fit(modelo, X[:,label.index('G1_ORR')], X[:,label.index('overpotencial_ORR_4e')], p0=[-1,8, X[np.argmin(X[:,label.index('overpotencial_ORR_4e')]),label.index('G1_ORR')] ,1])
			ss_res = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - modelo(X[:,label.index('G1_ORR')], *popt))**2  )	# Suma de los cuadrados de los residuos
			ss_tot = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - np.mean(X[:,label.index('overpotencial_ORR_4e')]) )**2  ) 		# Suma total de cuadrados
			ax21.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}G1_ORR+{popt[1]:.2} \n {popt[3]:.2}G1_ORR+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
			ax21.legend(loc='best')

			Xmin, Xmax = np.min(X[:,label.index('G2_ORR')])*0.9, np.max(X[:,label.index('G2_ORR')])*1.1
			Xnew = np.linspace(Xmin, Xmax, 100)
			popt, pcov = curve_fit(modelo, X[:,label.index('G2_ORR')], X[:,label.index('overpotencial_ORR_4e')], p0=[-1,0, X[np.argmin(X[:,label.index('overpotencial_ORR_4e')]),label.index('G2_ORR')] ,1])
			ss_res = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - modelo(X[:,label.index('G2_ORR')], *popt))**2  )	# Suma de los cuadrados de los residuos
			ss_tot = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - np.mean(X[:,label.index('overpotencial_ORR_4e')]) )**2  ) 		# Suma total de cuadrados
			ax22.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}G2_ORR+{popt[1]:.2} \n {popt[3]:.2}G2_ORR+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
			ax22.legend(loc='best')

			Xmin, Xmax = np.min(X[:,label.index('G3_ORR')])*0.9, np.max(X[:,label.index('G3_ORR')])*1.1
			Xnew = np.linspace(Xmin, Xmax, 100)
			popt, pcov = curve_fit(modelo, X[:,label.index('G3_ORR')], X[:,label.index('overpotencial_ORR_4e')], p0=[-1,8, X[np.argmin(X[:,label.index('overpotencial_ORR_4e')]),label.index('G3_ORR')] ,1])
			ss_res = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - modelo(X[:,label.index('G3_ORR')], *popt))**2  )	# Suma de los cuadrados de los residuos
			ss_tot = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - np.mean(X[:,label.index('overpotencial_ORR_4e')]) )**2  ) 		# Suma total de cuadrados
			ax23.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}G3_ORR+{popt[1]:.2} \n {popt[3]:.2}G3_ORR+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
			ax23.legend(loc='best')


			Xmin, Xmax = np.min(X[:,label.index('Gabs_OH')])*0.9, np.max(X[:,label.index('Gabs_OH')])*1.1
			Xnew = np.linspace(Xmin, Xmax, 100)
			z = np.polyfit(X[:,label.index('Gabs_OH')], X[:,label.index('Gabs_O')] , degree)
			poli11 = np.poly1d(z)
			ss_res = np.sum( (X[:,label.index('Gabs_O')] - poli11(X[:,label.index('Gabs_OH')]))**2  )	# Suma de los cuadrados de los residuos
			ss_tot = np.sum( (X[:,label.index('Gabs_O')] - np.mean(X[:,label.index('Gabs_OH')]) )**2  ) 		# Suma total de cuadrados
			ax31.plot( Xnew, 	poli11(Xnew), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{z[0]:.2}Gabs_OH+{z[1]:.2} \n R2={1-(ss_res/ss_tot):.2}' )
			ax31.legend(loc='best')

			Xmin, Xmax = np.min(X[:,label.index('Gabs_OOH')])*0.9, np.max(X[:,label.index('Gabs_OOH')])*1.1
			Xnew = np.linspace(Xmin, Xmax, 100)
			z = np.polyfit(X[:,label.index('Gabs_OOH')], X[:,label.index('Gabs_OH')] , degree)
			poli11 = np.poly1d(z)
			ss_res = np.sum( (X[:,label.index('Gabs_OH')] - poli11(X[:,label.index('Gabs_OOH')]))**2  )	# Suma de los cuadrados de los residuos
			ss_tot = np.sum( (X[:,label.index('Gabs_OH')] - np.mean(X[:,label.index('Gabs_OOH')]) )**2  ) 		# Suma total de cuadrados
			ax32.plot( Xnew, 	poli11(Xnew), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{z[0]:.2}Gabs_OOH+{z[1]:.2} \n R2={1-(ss_res/ss_tot):.2}' )
			ax32.legend(loc='best')

			Xmin, Xmax = np.min(X[:,label.index('Gabs_OOH')])*0.9, np.max(X[:,label.index('Gabs_OOH')])*1.1
			Xnew = np.linspace(Xmin, Xmax, 100)
			z = np.polyfit(X[:,label.index('Gabs_OOH')], X[:,label.index('Gabs_O')] , degree)
			poli11 = np.poly1d(z)
			ss_res = np.sum( (X[:,label.index('Gabs_O')] - poli11(X[:,label.index('Gabs_OOH')]))**2  )	# Suma de los cuadrados de los residuos
			ss_tot = np.sum( (X[:,label.index('Gabs_O')] - np.mean(X[:,label.index('Gabs_OOH')]) )**2  ) 		# Suma total de cuadrados
			ax33.plot( Xnew, 	poli11(Xnew), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{z[0]:.2}Gabs_OOH+{z[1]:.2} \n R2={1-(ss_res/ss_tot):.2}' )
			ax33.legend(loc='best')

		#plt.show()
		period_N, column_N, valence_N = [], [], []
		for n in N:
			for key, item in self.atomic_metadata.items():
				if key in n: 
					period_N.append(item[0]), column_N.append(item[1]), valence_N.append(np.mean(item[3]))
		period_N, column_N, valence_N = np.array(period_N), np.array(column_N), np.array(valence_N)

		orbitals = [XYd, XYu, YZd, YZu, Z2d, Z2u, XZd, XZu, Y2d, Y2u] = [ Y[:, y*500+8:(y+1)*500+8] for y in range(10) ]
		LCobitals = [	XYd+XYu, 				YZd+YZu, 						Z2d+Z2u, 
						XZd+XZu, 				Y2d+Y2u, 						(XYd+XYu)+(Y2d+Y2u), 
						(YZd+YZu)+(XZd+XZu), 	(YZd+YZu)+(XZd+XZu)+(Z2d+Z2u), 	XYd+XYu+YZd+YZu+Z2d+Z2u+XZd+XZu+Y2d+Y2u,
						]

		LCnames = [	'XYd+XYu', 				'YZd+YZu', 							'Z2d+Z2u', 
					'XZd+XZu', 				'Y2d+Y2u', 							'(XYd+XYu)+(Y2d+Y2u)', 
					'(YZd+YZu)+(XZd+XZu)', 	'(YZd+YZu)+(XZd+XZu)+(Z2d+Z2u)', 	'XYd+XYu+YZd+YZu+Z2d+Z2u+XZd+XZu+Y2d+Y2u',
					'period_N', 			'column_N', 						'valence_N'		]
		
		# BC.shape = (samples, factors)    period_N, 				column_N, 						valence_N	
		for pp in range(4):
			BC = np.array([ [band_center(orbital[n, :], -5, 5)  for n in range(samples)] for orbital in LCobitals]).T 
			BC = np.concatenate( (BC, np.array([period_N, column_N, valence_N]).T), axis=1)

			fig, ([ax11, ax12, ax13, ax14], [ax21, ax22, ax23, ax24], [ax31, ax32, ax33, ax34]) = plt.subplots(3, 4, figsize=(13, 10) )
			ax11.set(xlabel=f'BC_{LCnames[pp*3+0]}', ylabel='overpotencial_ORR_4e  (eV)')
			ax12.set(xlabel=f'BC_{LCnames[pp*3+0]}', ylabel='Gabs_O    (eV)')
			ax13.set(xlabel=f'BC_{LCnames[pp*3+0]}', ylabel='Gabs_OH   (eV)')
			ax14.set(xlabel=f'BC_{LCnames[pp*3+0]}', ylabel='Gabs_OOH  (eV)')

			ax21.set(xlabel=f'BC_{LCnames[pp*3+1]}', ylabel='overpotencial_ORR_4e  (eV)')
			ax22.set(xlabel=f'BC_{LCnames[pp*3+1]}', ylabel='Gabs_O    (eV)')
			ax23.set(xlabel=f'BC_{LCnames[pp*3+1]}', ylabel='Gabs_OH   (eV)')
			ax24.set(xlabel=f'BC_{LCnames[pp*3+1]}', ylabel='Gabs_OOH  (eV)')

			ax31.set(xlabel=f'BC_{LCnames[pp*3+2]}', ylabel='overpotencial_ORR_4e  (eV)')
			ax32.set(xlabel=f'BC_{LCnames[pp*3+2]}', ylabel='Gabs_O    (eV)')
			ax33.set(xlabel=f'BC_{LCnames[pp*3+2]}', ylabel='Gabs_OH   (eV)')
			ax34.set(xlabel=f'BC_{LCnames[pp*3+2]}', ylabel='Gabs_OOH  (eV)')

			fig.tight_layout()

			for n in range(samples):
				for key, metadata in self.atomic_metadata.items():
					if key in N[n]: break

				for key, marker in self.functional_metadata.items():
					if key in N[n]: marker=marker['marker']; break

				ax11.plot( BC[n,pp*3+0], X[n,label.index('overpotencial_ORR_4e')], marker=marker , color=metadata[4] )
				if text: ax11.text( BC[n,pp*3+0], X[n,label.index('overpotencial_ORR_4e')], N[n], color=metadata[4] )

				ax12.plot( BC[n,pp*3+0], X[n,label.index('Gabs_O')], marker=marker , color=metadata[4] )
				if text: ax12.text( BC[n,pp*3+0], X[n,label.index('Gabs_O')], N[n], color=metadata[4] )

				ax13.plot( BC[n,pp*3+0], X[n,label.index('Gabs_OH')], marker=marker , color=metadata[4] )
				if text: ax13.text( BC[n,pp*3+0], X[n,label.index('Gabs_OH')], N[n], color=metadata[4] )

				ax14.plot( BC[n,pp*3+0], X[n,label.index('Gabs_OOH')], marker=marker , color=metadata[4] )
				if text: ax14.text( BC[n,pp*3+0], X[n,label.index('Gabs_OOH')], N[n], color=metadata[4] )


				ax21.plot( BC[n,pp*3+1], X[n,label.index('overpotencial_ORR_4e')], marker=marker , color=metadata[4] )
				if text: ax21.text( BC[n,pp*3+1], X[n,label.index('overpotencial_ORR_4e')], N[n], color=metadata[4] )

				ax22.plot( BC[n,pp*3+1], X[n,label.index('Gabs_O')], marker=marker , color=metadata[4] )
				if text: ax22.text( BC[n,pp*3+1], X[n,label.index('Gabs_O')], N[n], color=metadata[4] )

				ax23.plot( BC[n,pp*3+1], X[n,label.index('Gabs_OH')], marker=marker , color=metadata[4] )
				if text: ax23.text( BC[n,pp*3+1], X[n,label.index('Gabs_OH')], N[n], color=metadata[4] )

				ax24.plot( BC[n,pp*3+1], X[n,label.index('Gabs_OOH')], marker=marker , color=metadata[4] )
				if text: ax24.text( BC[n,pp*3+1], X[n,label.index('Gabs_OOH')], N[n], color=metadata[4] )


				ax31.plot( BC[n,pp*3+2], X[n,label.index('overpotencial_ORR_4e')], marker=marker , color=metadata[4] )
				if text: ax31.text( BC[n,pp*3+2], X[n,label.index('overpotencial_ORR_4e')], N[n], color=metadata[4] )

				ax32.plot( BC[n,pp*3+2], X[n,label.index('Gabs_O')], marker=marker , color=metadata[4] )
				if text: ax32.text( BC[n,pp*3+2], X[n,label.index('Gabs_O')], N[n], color=metadata[4] )

				ax33.plot( BC[n,pp*3+2], X[n,label.index('Gabs_OH')], marker=marker , color=metadata[4] )
				if text: ax33.text( BC[n,pp*3+2], X[n,label.index('Gabs_OH')], N[n], color=metadata[4] )

				ax34.plot( BC[n,pp*3+2], X[n,label.index('Gabs_OOH')], marker=marker , color=metadata[4] )
				if text: ax34.text( BC[n,pp*3+2], X[n,label.index('Gabs_OOH')], N[n], color=metadata[4] )

			if curve_fit:
				from scipy.optimize import curve_fit

				def modelo(x, a, b, c, d):
					Y = np.zeros_like(x)
					Y[x<=c] = a*x[x<=c] + b 
					Y[x> c] = d*x[x> c] + (a-d)*c+b
					return Y


				try: 
					Xmin, Xmax = np.min(BC[:,pp*3+0][~np.isnan(BC[:,pp*3+0])])*0.9, np.max(BC[:,pp*3+0][~np.isnan(BC[:,pp*3+0])])*1.1
					Xnew = np.linspace(Xmin, Xmax, 100)
					popt, pcov = curve_fit(modelo, BC[:,pp*3+0][~np.isnan(BC[:,pp*3+0])], X[:,label.index('overpotencial_ORR_4e')][~np.isnan(BC[:,pp*3+0])], p0=[-1, 5, BC[np.argmin(X[:,label.index('overpotencial_ORR_4e')]),pp*3+0], 1])
					ss_res = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - modelo(BC[:,pp*3+0], *popt))**2  )	# Suma de los cuadrados de los residuos
					ss_tot = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - np.mean(X[:,label.index('overpotencial_ORR_4e')]) )**2  ) 		# Suma total de cuadrados
					ax11.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}{LCnames[pp*3+0]}+{popt[1]:.2} \n {popt[3]:.2}{LCnames[pp*3+0]}+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
					ax11.legend(loc='best')
				except: pass 

				try:
					Xmin, Xmax = np.min(BC[:,pp*3+0][~np.isnan(BC[:,pp*3+0])])*0.9, np.max(BC[:,pp*3+0][~np.isnan(BC[:,pp*3+0])])*1.1
					Xnew = np.linspace(Xmin, Xmax, 100)
					popt, pcov = curve_fit(modelo, BC[:,pp*3+0][~np.isnan(BC[:,pp*3+0])], X[:,label.index('Gabs_O')][~np.isnan(BC[:,pp*3+0])], p0=[-1, 5, BC[np.argmin(X[:,label.index('Gabs_O')]),pp*3+0], -1])
					ss_res = np.sum( (X[:,label.index('Gabs_O')] - modelo(BC[:,pp*3+0], *popt))**2  )	# Suma de los cuadrados de los residuos
					ss_tot = np.sum( (X[:,label.index('Gabs_O')] - np.mean(X[:,label.index('Gabs_O')]) )**2  ) 		# Suma total de cuadrados
					ax12.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}{LCnames[pp*3+0]}+{popt[1]:.2} \n {popt[3]:.2}{LCnames[pp*3+0]}+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
					ax12.legend(loc='best')
				except: 
					pass
					
				try:
					Xmin, Xmax = np.min(BC[:,pp*3+0][~np.isnan(BC[:,pp*3+0])])*0.9, np.max(BC[:,pp*3+0][~np.isnan(BC[:,pp*3+0])])*1.1
					Xnew = np.linspace(Xmin, Xmax, 100)
					popt, pcov = curve_fit(modelo, BC[:,pp*3+0][~np.isnan(BC[:,pp*3+0])], X[:,label.index('Gabs_OH')][~np.isnan(BC[:,pp*3+0])], p0=[1, 5, BC[np.argmin(X[:,label.index('Gabs_OH')]),pp*3+0], -1])
					ss_res = np.sum( (X[:,label.index('Gabs_OH')] - modelo(BC[:,pp*3+0], *popt))**2  )	# Suma de los cuadrados de los residuos
					ss_tot = np.sum( (X[:,label.index('Gabs_OH')] - np.mean(X[:,label.index('Gabs_OH')]) )**2  ) 		# Suma total de cuadrados
					ax13.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}{LCnames[pp*3+0]}+{popt[1]:.2} \n {popt[3]:.2}{LCnames[pp*3+0]}+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
					ax13.legend(loc='best')
				except: 
					pass
					
				try:
					Xmin, Xmax = np.min(BC[:,pp*3+0][~np.isnan(BC[:,pp*3+0])])*0.9, np.max(BC[:,pp*3+0][~np.isnan(BC[:,pp*3+0])])*1.1
					Xnew = np.linspace(Xmin, Xmax, 100)
					popt, pcov = curve_fit(modelo, BC[:,pp*3+0][~np.isnan(BC[:,pp*3+0])], X[:,label.index('Gabs_OOH')][~np.isnan(BC[:,pp*3+0])], p0=[1, 5, BC[np.argmin(X[:,label.index('Gabs_OOH')]),pp*3+0], -1])
					ss_res = np.sum( (X[:,label.index('Gabs_OOH')] - modelo(BC[:,pp*3+0], *popt))**2  )	# Suma de los cuadrados de los residuos
					ss_tot = np.sum( (X[:,label.index('Gabs_OOH')] - np.mean(X[:,label.index('Gabs_OOH')]) )**2  ) 		# Suma total de cuadrados
					ax14.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}{LCnames[pp*3+0]}+{popt[1]:.2} \n {popt[3]:.2}{LCnames[pp*3+0]}+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
					ax14.legend(loc='best')
				except: 
					pass
					

				try:
					Xmin, Xmax = np.min(BC[:,pp*3+1][~np.isnan(BC[:,pp*3+1])])*0.9, np.max(BC[:,pp*3+1][~np.isnan(BC[:,pp*3+1])])*1.1
					Xnew = np.linspace(Xmin, Xmax, 100)
					popt, pcov = curve_fit(modelo, BC[:,pp*3+1][~np.isnan(BC[:,pp*3+1])], X[:,label.index('overpotencial_ORR_4e')][~np.isnan(BC[:,pp*3+1])], p0=[1, 5, BC[np.argmin(X[:,label.index('overpotencial_ORR_4e')]),pp*3+1], -1])
					ss_res = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - modelo(BC[:,pp*3+1], *popt))**2  )	# Suma de los cuadrados de los residuos
					ss_tot = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - np.mean(X[:,label.index('overpotencial_ORR_4e')]) )**2  ) 		# Suma total de cuadrados
					ax21.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}{LCnames[pp*3+1]}+{popt[1]:.2} \n {popt[3]:.2}{LCnames[pp*3+0]}+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
					ax21.legend(loc='best')
				except: 
					pass
					
				try:
					Xmin, Xmax = np.min(BC[:,pp*3+1][~np.isnan(BC[:,pp*3+1])])*0.9, np.max(BC[:,pp*3+1][~np.isnan(BC[:,pp*3+1])])*1.1
					Xnew = np.linspace(Xmin, Xmax, 100)
					popt, pcov = curve_fit(modelo, BC[:,pp*3+1][~np.isnan(BC[:,pp*3+1])], X[:,label.index('Gabs_O')][~np.isnan(BC[:,pp*3+1])], p0=[1, 5, BC[np.argmin(X[:,label.index('Gabs_O')]),pp*3+1], -1])
					ss_res = np.sum( (X[:,label.index('Gabs_O')] - modelo(BC[:,pp*3+1], *popt))**2  )	# Suma de los cuadrados de los residuos
					ss_tot = np.sum( (X[:,label.index('Gabs_O')] - np.mean(X[:,label.index('Gabs_O')]) )**2  ) 		# Suma total de cuadrados
					ax22.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}{LCnames[pp*3+1]}+{popt[1]:.2} \n {popt[3]:.2}{LCnames[pp*3+0]}+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
					ax22.legend(loc='best')
				except: 
					pass
					
				try:
					Xmin, Xmax = np.min(BC[:,pp*3+1][~np.isnan(BC[:,pp*3+1])])*0.9, np.max(BC[:,pp*3+1][~np.isnan(BC[:,pp*3+1])])*1.1
					Xnew = np.linspace(Xmin, Xmax, 100)
					popt, pcov = curve_fit(modelo, BC[:,pp*3+1][~np.isnan(BC[:,pp*3+1])], X[:,label.index('Gabs_OH')][~np.isnan(BC[:,pp*3+1])], p0=[1, 5, BC[np.argmin(X[:,label.index('Gabs_OH')]),pp*3+1], -1])
					ss_res = np.sum( (X[:,label.index('Gabs_OH')] - modelo(BC[:,pp*3+1], *popt))**2  )	# Suma de los cuadrados de los residuos
					ss_tot = np.sum( (X[:,label.index('Gabs_OH')] - np.mean(X[:,label.index('Gabs_OH')]) )**2  ) 		# Suma total de cuadrados
					ax23.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}{LCnames[pp*3+1]}+{popt[1]:.2} \n {popt[3]:.2}{LCnames[pp*3+0]}+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
					ax23.legend(loc='best')
				except: 
					pass
					
				try:
					Xmin, Xmax = np.min(BC[:,pp*3+1][~np.isnan(BC[:,pp*3+1])])*0.9, np.max(BC[:,pp*3+1][~np.isnan(BC[:,pp*3+1])])*1.1
					Xnew = np.linspace(Xmin, Xmax, 100)
					popt, pcov = curve_fit(modelo, BC[:,pp*3+1][~np.isnan(BC[:,pp*3+1])], X[:,label.index('Gabs_OOH')][~np.isnan(BC[:,pp*3+1])], p0=[1, 5, BC[np.argmin(X[:,label.index('Gabs_OOH')]),pp*3+1], -1])
					ss_res = np.sum( (X[:,label.index('Gabs_OOH')] - modelo(BC[:,pp*3+1], *popt))**2  )	# Suma de los cuadrados de los residuos
					ss_tot = np.sum( (X[:,label.index('Gabs_OOH')] - np.mean(X[:,label.index('Gabs_OOH')]) )**2  ) 		# Suma total de cuadrados
					ax24.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}{LCnames[pp*3+1]}+{popt[1]:.2} \n {popt[3]:.2}{LCnames[pp*3+0]}+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
					ax24.legend(loc='best')
				except: 
					pass
					

				try:
					Xmin, Xmax = np.min(BC[:,pp*3+2][~np.isnan(BC[:,pp*3+2])])*0.9, np.max(BC[:,pp*3+2][~np.isnan(BC[:,pp*3+2])])*1.1
					Xnew = np.linspace(Xmin, Xmax, 100)
					popt, pcov = curve_fit(modelo, BC[:,pp*3+2][~np.isnan(BC[:,pp*3+2])], X[:,label.index('overpotencial_ORR_4e')][~np.isnan(BC[:,pp*3+2])], p0=[1, 5, BC[np.argmin(X[:,label.index('overpotencial_ORR_4e')]),pp*3+2], -1])
					ss_res = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - modelo(BC[:,pp*3+2], *popt))**2  )	# Suma de los cuadrados de los residuos
					ss_tot = np.sum( (X[:,label.index('overpotencial_ORR_4e')] - np.mean(X[:,label.index('overpotencial_ORR_4e')]) )**2  ) 		# Suma total de cuadrados
					ax31.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}{LCnames[pp*3+2]}+{popt[1]:.2} \n {popt[3]:.2}{LCnames[pp*3+0]}+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
					ax31.legend(loc='best')
				except: 
					pass
					
				try:
					Xmin, Xmax = np.min(BC[:,pp*3+2][~np.isnan(BC[:,pp*3+2])])*0.9, np.max(BC[:,pp*3+2][~np.isnan(BC[:,pp*3+2])])*1.1
					Xnew = np.linspace(Xmin, Xmax, 100)
					popt, pcov = curve_fit(modelo, BC[:,pp*3+2][~np.isnan(BC[:,pp*3+2])], X[:,label.index('Gabs_O')][~np.isnan(BC[:,pp*3+2])], p0=[1, 5, BC[np.argmin(X[:,label.index('Gabs_O')]),pp*3+2], -1])
					ss_res = np.sum( (X[:,label.index('Gabs_O')] - modelo(BC[:,pp*3+2], *popt))**2  )	# Suma de los cuadrados de los residuos
					ss_tot = np.sum( (X[:,label.index('Gabs_O')] - np.mean(X[:,label.index('Gabs_O')]) )**2  ) 		# Suma total de cuadrados
					ax32.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}{LCnames[pp*3+2]}+{popt[1]:.2} \n {popt[3]:.2}{LCnames[pp*3+0]}+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
					ax32.legend(loc='best')
				except: 
					pass
					
				try:
					Xmin, Xmax = np.min(BC[:,pp*3+2][~np.isnan(BC[:,pp*3+2])])*0.9, np.max(BC[:,pp*3+2][~np.isnan(BC[:,pp*3+2])])*1.1
					Xnew = np.linspace(Xmin, Xmax, 100)
					popt, pcov = curve_fit(modelo, BC[:,pp*3+2][~np.isnan(BC[:,pp*3+2])], X[:,label.index('Gabs_OH')][~np.isnan(BC[:,pp*3+2])], p0=[1, 5, BC[np.argmin(X[:,label.index('Gabs_OH')]),pp*3+2], -1])
					ss_res = np.sum( (X[:,label.index('Gabs_OH')] - modelo(BC[:,pp*3+2], *popt))**2  )	# Suma de los cuadrados de los residuos
					ss_tot = np.sum( (X[:,label.index('Gabs_OH')] - np.mean(X[:,label.index('Gabs_OH')]) )**2  ) 		# Suma total de cuadrados
					ax33.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}{LCnames[pp*3+2]}+{popt[1]:.2} \n {popt[3]:.2}{LCnames[pp*3+0]}+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
					ax33.legend(loc='best')
				except: 
					pass

				try:
					Xmin, Xmax = np.min(BC[:,pp*3+2][~np.isnan(BC[:,pp*3+2])])*0.9, np.max(BC[:,pp*3+2][~np.isnan(BC[:,pp*3+2])])*1.1
					Xnew = np.linspace(Xmin, Xmax, 100)
					popt, pcov = curve_fit(modelo, BC[:,pp*3+2][~np.isnan(BC[:,pp*3+2])], X[:,label.index('Gabs_OOH')][~np.isnan(BC[:,pp*3+2])], p0=[1, 5, BC[np.argmin(X[:,label.index('Gabs_OOH')]),pp*3+2], -1])
					ss_res = np.sum( (X[:,label.index('Gabs_OOH')] - modelo(BC[:,pp*3+2], *popt))**2  )	# Suma de los cuadrados de los residuos
					ss_tot = np.sum( (X[:,label.index('Gabs_OOH')] - np.mean(X[:,label.index('Gabs_OOH')]) )**2  ) 		# Suma total de cuadrados
					ax34.plot(Xnew, modelo(Xnew, *popt), marker=None, linestyle='dashed', color=(0.9,0.3,0.3), label=f'{popt[0]:.2}{LCnames[pp*3+2]}+{popt[1]:.2} \n {popt[3]:.2}{LCnames[pp*3+0]}+{(popt[0]-popt[3])*popt[2]:.2} \n cut at {popt[2]:.2} | R2={1-(ss_res/ss_tot):.2}')
					ax34.legend(loc='best')
				except: 
					pass

		plt.show()

	def neural_regresor_fit(self, X:np.ndarray=None, Y:np.ndarray=None, Z:np.ndarray=None, N:np.ndarray=None,
								label:str=None, act_func:str='tanh', 
								v:bool=True, save:bool=True):
		
		X 			= X 		if not X is None 		else self.X
		Y 			= Y 		if not Y is None 		else self.Y
		Z 			= Z 		if not Z is None 		else self.Z
		N 			= N 		if not N is None 		else self.N
		act_func 	= act_func 	if not act_func is None else self.act_func

		label = label if not label is None else [ 	'overpotencial_ORR_4e', 
													'Gabs_OOH', 'Eabs_OOH', 
													'Gabs_OH',  'Eabs_OH', 
													'Gabs_O', 	'Eabs_O', 
													'G1_ORR', 	'G2_ORR', 'G3_ORR', 'G4_ORR', 
													'limiting_step_ORR']

		def get_model_l0(width:int):
			my_l2 = l2(0.03)
			model = Sequential()
			model.add(Reshape((width, 1), input_shape=(width,)))
			model.add(Conv1D(filters=40, kernel_size=100, activation='tanh', kernel_regularizer=my_l2))
			model.add(MaxPooling1D(40))
			model.add(Flatten())
			model.add(Dense(30, activation='relu', kernel_regularizer=my_l2))
			model.add(Dense(30, activation='relu', kernel_regularizer=my_l2))
			model.add(Dense(30, activation='relu', kernel_regularizer=my_l2))
			model.add(Dense(1, activation='linear', kernel_regularizer=my_l2))
			model.compile(loss='mean_squared_error', optimizer='adam', metrics=[coeff_determination, ])

			return model

		def get_model_l2(width:int):
			my_l2 = l2(0.04)
			model = Sequential()
			model.add(Reshape((width, 1), input_shape=(width,)))
			model.add(Conv1D(filters=10, kernel_size=100, activation=act_func, kernel_regularizer=my_l2))
			model.add(MaxPooling1D(100))
			model.add(Flatten())
			model.add(Dense(60, activation='relu', kernel_regularizer=my_l2))
			model.add(Dense(1, activation='linear', kernel_regularizer=my_l2))
			model.compile(loss='mean_squared_error', optimizer='adam', metrics=[coeff_determination, ])

			return model

		def get_model_l1(width:int):
			my_l2 = l2(0.04)
			model = Sequential()
			model.add(Reshape((width, 1), input_shape=(width,)))
			model.add(Conv1D(filters=40, kernel_size=50, activation=act_func, kernel_regularizer=my_l2))
			model.add(MaxPooling1D(60))
			model.add(Flatten())
			model.add(Dense(40, activation='relu', kernel_regularizer=my_l2))
			model.add(Dense(1, activation='linear', kernel_regularizer=my_l2))
			model.compile(loss='mean_squared_error', optimizer='adam', metrics=[coeff_determination, ])

			return model

		def coeff_determination(y_true, y_pred):
			SS_res = K.sum(K.square(y_true-y_pred))
			SS_tot = K.sum(K.square(y_true - K.mean(y_true)))

			return 1 - SS_res/SS_tot

		try:	
			import keras.backend as K
			import keras.callbacks
			from keras.models import Sequential
			from keras.layers import Conv1D, Flatten, BatchNormalization, MaxPooling1D, Dense, Reshape
			from keras.regularizers import l2, l1
			from sklearn.model_selection import train_test_split
			import seaborn as sns
			sns.set_style('ticks')
		except: print('(!) ERROR :: can NOT import lib')

		# === MODEL data === #
		print( '>> == Data shape ==' )
		print( f' (OK) X : {X.shape}' )
		print( f' (OK) Y : {Y.shape}' )
		print( f' (OK) N : {N.shape}' )

		Y = Y[:,3]
		orbitals = [XYd, XYu, YZd, YZu, Z2d, Z2u, XZd, XZu, Y2d, Y2u] = [ X[:, x*500+8:(x+1)*500+8] for x in range(10) ]
		LCobitals = (	Z2d+Z2u,	
						(YZd+YZu)+(XZd+XZu), 	(YZd+YZu)+(XZd+XZu)+(Z2d+Z2u), 	XYd+XYu+YZd+YZu+Z2d+Z2u+XZd+XZu+Y2d+Y2u,
						)
		#LCobitals = tuple( [XYd+XYu, Y2d+Y2u] )
		#X = np.concatenate( ((YZd+YZu), (XZd+XZu), (Z2d+Z2u)), axis=1 )
		X = np.concatenate( LCobitals, axis=1 )
		
		# === set MODEL === #
		model = get_model_l0(X.shape[1])
		model.summary()
		X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1)
		history =  model.fit(X_train, y_train, epochs=50, validation_data=(X_test, y_test), verbose=2)

		# === PLOT data === #
		# regresion ECM vs epoch #
		fig, ([ax11, ax12, ax13]) = plt.subplots(1, 3, figsize=(15, 5) )
		fig.tight_layout()

		ax11.plot( history.history['loss'] )
		ax11.set(xlabel='epoch'	, ylabel='ECM', title='ECM vs. epochs'	)

		# data vs PREDICT #
		y_regr = model.predict(X_train)
		ax12.plot(y_train, y_regr, 'o', color=self.cnames['tomato'], label='trainnig data')
		y_regr = model.predict(X_test)
		ax12.plot(y_test, y_regr, 'o', color=self.cnames['royalblue'], label='test data')
		ax12.set(xlabel='y data'	, ylabel='y prediction', title='Regresion Curve'	)
		ss_res = np.sum( (np.concatenate( (y_train, y_test) ) - np.concatenate( (model.predict(X_train)[:,0], model.predict(X_test)[:,0]) ) )**2  )	# Suma de los cuadrados de los residuos
		ss_tot = np.sum( (np.concatenate( (y_train, y_test) ) - np.mean(np.concatenate( (y_train, y_test) )) )**2  ) 		# Suma total de cuadrados
		ax12.plot( [ 0.0, 5], [   0, 5.0] , ':' , alpha=0.5, c=(0,0,0), label=f'R2 = {1-(ss_res/ss_tot):.2}')
		ax12.plot( [ 0.0, 5], [ 0.5, 5.5] , '-' , alpha=0.5, c=(1,0,0), label='500meV error')
		ax12.plot( [ 0.0, 5], [-0.5, 4.5] , '-' , alpha=0.5, c=(1,0,0))
		ax12.legend(loc='best')

		# RMSE hist #
		ax13.hist( np.concatenate( (y_train-model.predict(X_train)[:,0], y_test-model.predict(X_test)[:,0]) ) )
		ax13.set(xlabel='nominal - prediction'	, ylabel='ECM', title='Residues'	)

		plt.show()

		# save model #
		if save: self.neural_regresor_model = model

		return True

	def neural_regresor_predict(self, X:np.ndarray=None, Y:np.ndarray=None, Z:np.ndarray=None, N:np.ndarray=None,
								act_func:str='tanh',
								v:bool=True):
		
		X 			= X 		if not X is None 		else self.X
		Y 			= Y 		if not Y is None 		else self.Y
		Z 			= Z 		if not Z is None 		else self.Z
		N 			= N 		if not N is None 		else self.N
		act_func 	= act_func 	if not act_func is None else self.act_func

		def coeff_determination(y_true, y_pred):
			SS_res = K.sum(K.square(y_true-y_pred))
			SS_tot = K.sum(K.square(y_true - K.mean(y_true)))

			return 1 - SS_res/SS_tot

		try:	
			import keras.backend as K
			import keras.callbacks
			from keras.models import Sequential
			from keras.layers import Conv1D, Flatten, BatchNormalization, MaxPooling1D, Dense, Reshape
			from keras.regularizers import l2, l1
			from sklearn.model_selection import train_test_split
			import seaborn as sns
			sns.set_style('ticks')
		except: print('(!) ERROR :: can NOT import lib')

		# === MODEL data === #
		print( '>> == Data shape ==' )
		print( f' (OK) X : {X.shape}' )
		print( f' (OK) Y : {Y.shape}' )
		print( f' (OK) N : {N.shape}' )

		Y = Y[:,3]
		orbitals = [XYd, XYu, YZd, YZu, Z2d, Z2u, XZd, XZu, Y2d, Y2u] = [ X[:, x*500+8:(x+1)*500+8] for x in range(10) ]
		LCobitals = (	
						Z2d+Z2u,	(YZd+YZu)+(XZd+XZu), 	(YZd+YZu)+(XZd+XZu)+(Z2d+Z2u), 	
						XYd+XYu+YZd+YZu+Z2d+Z2u+XZd+XZu+Y2d+Y2u,
						)
		#LCobitals = tuple( [XYd+XYu, Y2d+Y2u] )
		#X = np.concatenate( ((YZd+YZu), (XZd+XZu), (Z2d+Z2u)), axis=1 )
		X = np.concatenate( LCobitals, axis=1 )

		# === load MODEL === #
		model = self.neural_regresor_model

		# === PLOT data === #
		fig, ([ax11, ax12]) = plt.subplots(1, 2, figsize=(10, 5) )
		ax11.set(xlabel='y data'	, ylabel='y prediction', title='Regresion Curve'	)
		y_regr = model.predict(X)
		Xmin, Xmax = np.min(Y), np.max(Y)
		Ymin, Ymax = np.min(y_regr), np.max(y_regr)
		ss_res = np.sum( (Y - y_regr[:,0])**2  )			# Suma de los cuadrados de los residuos
		ss_tot = np.sum( (Y - np.mean(Y))**2  ) 		# Suma total de cuadrados
		print(Y.shape, y_regr.shape, Y - np.mean(Y))
		ax11.plot(Y, y_regr, 'o', color='r')
		ax11.plot( [ Xmin, Xmax], [ Xmin, Xmax] , ':' , alpha=0.5, c=(0,0,0), label=f'R2 = {1-(ss_res/ss_tot):.2}')
		ax11.plot( [ Xmin, Xmax], [ Xmin+0.5, Xmax+0.5] , '-' , alpha=0.5, c=(1,0,0)) 
		ax11.plot( [ Xmin, Xmax], [ Xmin-0.5, Xmax-0.5] , '-' , alpha=0.5, c=(1,0,0)) 
		ax11.legend(loc='best')	
		for i, n in enumerate(Y-y_regr[:,0]):
			if abs(n) > 0.5:	ax11.text( Y[i], y_regr[i,0], N[i] , alpha=0.6, c=(0,0,0))

		# RMSE hist #
		ax12.hist( Y-model.predict(X)[:,0] )
		ax12.set(xlabel='nominal - prediction'	, ylabel='ECM', title='Residues'	)

		plt.show()

# ==== QM00 ==== # # ==== QM00 ==== # # ==== QM00 ==== #
qm = QMlearning()

# --> LOAD <-- # 
#qm.load_folder('/home/akaris/Documents/code/VASP/v4.6/files/dataset/Metales/pkl/dataset/MPC1')
#qm.load_folder('/home/akaris/Documents/code/VASP/v4.6/files/dataset/Metales/pkl/dataset/MPC2')
qm.load_folder('/home/akaris/Documents/code/VASP/v4.6/files/dataset/Metales/pkl/dataset/MPCwoAu')

#qm.load_folder('/home/akaris/Documents/code/VASP/v4.6/files/dataset/Metales/pkl/dataset/MTPyP1')
#qm.load_folder('/home/akaris/Documents/code/VASP/v4.6/files/dataset/Metales/pkl/dataset/MTPyP2')
#qm.load_folder('/home/akaris/Documents/code/VASP/v4.6/files/dataset/Metales/pkl/dataset/MTPyPwoAu')


# --> FILTER <-- # 
qm.filter( filters=[{ 'type':'N', 'condition':'whitelist', 'list':['Sc','Ti','V','Cr','Mn',  'Fe','Co','Ni','Cu','Zn',] }] )
#qm.filter( filters=[{ 'type':'N', 'condition':'whitelist', 'list':['Sc','Ti','V','Cr','Mn',  'Co','Cu','Zn',] }] )
#qm.filter( filters=[{ 'type':'N', 'condition':'whitelist', 'list':['DF2',] }] )

# --> NNmodel <-- # 
#qm.neural_regresor_fit( X=qm.Y, Y=qm.X )

# --> Data analysis <-- # 
qm.ORR_analysis(text=False)

# ==== QM01 ==== # # ==== QM01 ==== # # ==== QM01 ==== #
qm1 = QMlearning()

# --> LOAD <-- # 
#qm1.load_folder('/home/akaris/Documents/code/VASP/v4.6/files/dataset/Metales/pkl/dataset/MPC1')
#qm1.load_folder('/home/akaris/Documents/code/VASP/v4.6/files/dataset/Metales/pkl/dataset/MPC2')
qm1.load_folder('/home/akaris/Documents/code/VASP/v4.6/files/dataset/Metales/pkl/dataset/MTPyP1')
qm1.load_folder('/home/akaris/Documents/code/VASP/v4.6/files/dataset/Metales/pkl/dataset/MTPyP2')

# --> FILTER <-- # 
qm1.filter( filters=[{ 'type':'N', 'condition':'whitelist', 'list':['Fe','Ni',] }] )
qm1.filter( filters=[{ 'type':'N', 'condition':'whitelist', 'list':['D3'] }] )

# --> NNmodel <-- # 
qm1.neural_regresor_model = qm.neural_regresor_model
qm1.neural_regresor_predict( X=qm1.Y, Y=qm1.X )



plt.show()

# plot vs centro de banda


