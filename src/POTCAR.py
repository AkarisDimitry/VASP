import numpy as np 	
import matplotlib.pyplot as plt
try:	
	from ase.visualize import view
	from ase import Atoms
	from ase.visualize.plot import plot_atoms
except: print('can not load ASE module. (pipX install ase-atomistics)')

class POTCAR(object):
	def __init__(self, name=None, atoms=None, atoms_list=None):
		self.name = name
		self.path = None
		self.file_name = None

		self.tag_list = [	'TITEL', 'LULTRA', 'IUNSCR', 'RPACOR', 'POMASS', 'ZVAL', 
							'RCORE', 'RWIGS', 'ENMAX', 'ENMIN', 'RCLOC', 'LCOR', 
							'LPAW', 'EAUG', 'RMAX', 'RAUG', 'RDEP', 'RDEPT']
		self.tag_dict = {key:[] for key in self.tag_list}
		'''
				'TITEL':	[], 
				'LULTRA':	[], 
				'IUNSCR':	[], 
				'RPACOR':	[],
				'POMASS':	[], 
				'ZVAL':		[], 
				'RWIGS':	[], 
				'ENMAX':	[], 
				'ENMIN':	[], 
				'RCLOC':	[], 
				'LCOR':		[], 
				'LPAW':		[], 
				'EAUG':		[], 
				'RMAX':		[], 
				'RAUG':		[], 
				'RDEP':		[], 
				'RDEPT':	[], 
							}
		'''

		self.TITEL = self.tag_dict['TITEL'] 	# eg: ('H') 	| TITEL
		self.LULTRA = self.tag_dict['LULTRA']   # eg: (F) 		| use ultrasoft PP ?
		self.IUNSCR = self.tag_dict['IUNSCR']  	# eg: (0) 		| unscreen: 0-lin 1-nonlin 2-no
		self.RPACOR = self.tag_dict['RPACOR']   # eg: (0.000) 	| partial core radius
		self.POMASS = self.tag_dict['POMASS'] 	# eg: (1.000) 	| partial core radius
		self.ZVAL   = self.tag_dict['ZVAL']  	# eg: (1.000) 	| mass and valenz
		self.RCORE  = self.tag_dict['RCORE'] 	# eg: (1.100) 	| outmost cutoff radius
		self.RWIGS  = self.tag_dict['RWIGS'] 	# eg: (0.700) 	| RWIGS  =    0.370    wigner-seitz radius (au A)
		self.ENMAX  = self.tag_dict['ENMAX']  	# eg: (250.000) | eV
		self.ENMIN  = self.tag_dict['ENMIN']  	# eg: (200.000) | eV
		self.RCLOC  = self.tag_dict['RCLOC']   	# eg: (0.701) 	| cutoff for local pot
		self.LCOR   = self.tag_dict['LCOR']    	# eg: (T) 		| correct aug charges
		self.LPAW   = self.tag_dict['LPAW']    	# eg: (T) 		| paw PP
		self.EAUG   = self.tag_dict['EAUG']
		self.RMAX   = self.tag_dict['RMAX']    	# eg: (1.123) 	| core radius for proj-oper
		self.RAUG   = self.tag_dict['RAUG']    	# eg: (1.200) 	| factor for augmentation sphere
		self.RDEP   = self.tag_dict['RDEP']    	# eg: (1.112) 	| radius for radial grids
		self.RDEPT  = self.tag_dict['RDEPT']   	# eg: (0.926) 	| core radius for aug-charge


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

	def load(self, file_name=None, save=True, v=0):
		if type(file_name) == str and file_name != None: 
			self.file_name = file_name

		try: 		f = open(self.file_name,'r')
		except: 	print('ERROR :: POTCAR.load() :: missing file {}'.format(file_name) ); 	pass

		for n in self.tag_list:
			if type(self.tag_dict[n]) != list: 
				try:self.tag_dict[n] = []
				except: print('ERROR :: POTCAR.load() :: error KEY ({}) in self.tag_dict'.format(n) ); pass
			else: pass

		
		for i, n in enumerate(f): # read the file {file_name} line by line
			vec = [m for m in n.split(' ') if m != '' and m != '\n'] # clean the n line
			if len(vec) > 0 and vec[0] in self.tag_dict.keys():
				try:
					if v>1: print('reading {} '.format(vec[0])) # verbose
					if vec[0] == 'TITEL': 	# TITEL case
						if v>0: print('new ATOM detected {} '.format(vec[3])) 				# verbose
						self.tag_dict[vec[0]].append( vec[3] ) 								# TITLE read case
					elif self.isnum(vec[2]):  self.tag_dict[vec[0]].append( float(vec[2]) ) # numeric read case
					else:  self.tag_dict[vec[0]].append( vec[2] )							# NON numeric read case
				except: print('ERROR :: POTCAR.load() :: error while reading ({}) in file {}'.format(vec[0], file_name) ); pass

		if v>0: print('SUCCESSFUL READ POTCAR :: from file {}'.format(file_name) ) # verbose
		if v>1: print('{}'.format(self.tag_dict) ) # verbosereturn self.tag_dict

	def summary(self, ): 
		potcar_resumen = self.resume()

		print('[POTCAR]\t\t' + ','.join(potcar_resumen))

	def resume(self, ):
		return 	self.TITEL

# *** eg. usage:
potcar = POTCAR()
potcar.load(file_name='/home/akaris/Documents/code/VASP/v4.0/files/POTCARs', v=False)
potcar.summary()
''' 
'''

