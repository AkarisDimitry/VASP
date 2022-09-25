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

# *** import python libraries
import re

class OSZICAR(object):
	def __init__(self, name=None, ):
		self.name = name
		self.file_name = None

		self.ionic_step = None
		self.electronic_step = None

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

	def E(self, ): 
		try:		return self.ionic_step[2][-1]
		except:		return None

	def isnum(self, n):
		# ------------------ Define if n is or not a number ------------------ # 
		# n     :   VAR     :   VAR to check if it is a numerical VAR
		# return :  BOOL    : True/False
		try: float(n); return True
		except: return False

	def load(self, file_name=None, save=True, low_mem=False):
		if file_name != None: self.file_name = file_name

		try: 	f = open(self.file_name,'r')
		except: 	print('ERROR :: OSZICAR.load() :: missing file {}'.format(file_name) ); 	return

		ionic_step, ionic_step_dict, electronic_step = [], {}, [[]]
		E, F, mag, N= [], [], [], 0
		variables = {	'F':'F', 	'E0':'E', 	'E=':'dE','mag':'mag', 	'EK':'EK', 	'T':'T', 	'SP':'SP', 	'SK':'SK'  }

		for i, line in enumerate(f):
			vec = [float(m) if self.isnum(m) else m for m in line.split(' ') if m != '']

			if self.isnum(vec[0]):
				dict_ionic = dict(filter(lambda elem: type(elem[0]) != type(None), { i[0]:i[1] for i in [ [variables[list(filter(lambda x: type(x)!=type(None), [key if type(v) == str and key in v else None for key, value in variables.items() ]))[0]] , vec[i+1]] if sum([1 if type(v) == str and key in v else 0 for key, value in variables.items() ])>0 else [None, None] for i, v in enumerate(vec) ] }.items() ))
				dict_ionic['N'] = int(vec[0])

				ionic_step_dict[i] = dict_ionic
				N += 1
				F += [dict_ionic['F']]
				E += [dict_ionic['E']]
				mag += [dict_ionic['mag']]

				try:
					ionic_step.append( list(filter(lambda x: self.isnum(x), vec )) )
					if not low_mem: 
						electronic_step.append([])
				except: print('ERROR :: OSZICAR.load() :: error while READING OSZICAR' )
				
			elif 'DAV' in vec[0] and not low_mem:
				try:	electronic_step[-1].append( list(filter(lambda x: self.isnum(x), vec )) )
				except: print('ERROR :: OSZICAR.load() :: error while READING OSZICAR' )

		# close the file 
		try: 		f.close()
		except: 	print('ERROR :: OSZICAR.load() :: can NOT close file {}'.format(file_name) ); 	return

		if save:
			self.N = N
			self.E = E
			self.F = F
			self.mag = mag

			self.ionic_step_dict 	= []

			self.ionic_step 		= ionic_step
			self.electronic_step	= electronic_step[:-1]
		
		return ionic_step

	def get_E(self, ): return self.ionic_step[-1][1]

	def summary(self, ): 
		OSZICAR_resumen = self.resume()

		toprint = '    |       *-->[OSZICAR] ::\n'
		for osr in OSZICAR_resumen['list']:	toprint += '    |         \t{}:{}\n'.format(osr[0], osr[1] )
		print(toprint)
		
		return None

	def resume(self, ):
		try:
			if type(self.ionic_step) == list and len(self.ionic_step)>0:	
				if len(self.ionic_step[-1]) > 3:
					info_OSZICAR = { 'N': len(self.ionic_step), 'E':self.E[-1], 'MAG':self.ionic_step[-1][-1], }
				else:
					info_OSZICAR = { 'N': len(self.ionic_step), 'E':self.E[-1] }
			else: 								
				info_OSZICAR = { 'N': 0, 'ionic_step':0, }
		except OSError as err:
			print("OS error: {0}".format(err))
			info_OSZICAR = { 'N': 0, 'ionic_step':0, }
		except ValueError:
			print("Could not convert data to an integer.")
			info_OSZICAR = { 'N': 0, 'ionic_step':0, }
		except:
			print("Unexpected error:", sys.exc_info()[0])
			info_OSZICAR = { 'N': 0, 'ionic_step':0, }

		return {'list': [[key, value] for i, (key, value) in enumerate( info_OSZICAR.items())], 'dict':info_OSZICAR}

'''
# === Eg. load OSZICAR === #
oszicar = OSZICAR()
oszicar.load('/home/akaris/Documents/code/VASP/v4.1/files/OSZICAR/OSZICAR_FePC_IBRION2')
oszicar.summary()
'''



