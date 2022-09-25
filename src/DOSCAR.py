import numpy as np 	
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import argparse

try:	from src import Logs
except:	
	try: import Logs as Logs
	except: print('WARNING :: Set.import_libraries() :: can not import ORR ')

try:		from scipy.interpolate import interp1d
except: 	print('WARNNING :: can NOT load scipy.interpolate :: function cut will be disabled')

class DOSCAR(object):
	def __init__(self, name=None, plot_ions=None, plot_orbitals=None):
		self.name = name
		self.file_name = None

		self.n_E = None
		self.n_oins = None
		self.n_orb = None

		self.fermi = None
		self.E = None # kpoint
		self.E_total = None # kpoint

		self.orbital_name = {
							0: 'E',

							1: 's(u)',
							2: 's(d)',

							3: 'py(u)',
							4: 'py(d)',
							5: 'pz(u)',
							6: 'pz(d)',
							7: 'p1(u)',
							8: 'p1(d)',

							9: 'd_xy(u)',
							10:'d_xy(d)',
							11:'d_yz(u)',
							12:'d_yz(d)',
							13:'d_z2-r2(u)',
							14:'d_z2-r2(d)',
							15:'d_xz(u)',
							16:'d_xz(d)',
							17:'d_x2-y2(u)',
							18:'d_x2-y2(d)',
							}

		self.orbital_index = {
							'all'	 :[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
							's(u)'	:	[1],
							's(d)'	:	[2],
							's'		:	[1,2],
							
							'py(u)'	:	[3],
							'py(d)'	:	[4],
							'py'	:	[3, 4],
							'pz(u)'	:	[5],
							'pz(d)'	:	[6],
							'pz'	:	[5, 6],
							'px(u)'	:	[7],
							'px(d)'	:	[8],
							'px'	:	[7, 8],
							'p'		:	[3, 4, 5, 6, 7, 8],

							'd_xy(u)'	:	[9],
							'd_xy(d)'	:	[10],
							'd_xy'		:	[9, 10],
							'd_yz(u)'	:	[11],
							'd_yz(d)'	:	[12],
							'd_yz'	 	:	[11, 12],
							'd_z2-r2(u)':	[13],
							'd_z2-r2(d)':	[14],
							'd_z2-r2'	:	[13, 14],
							'd_xz(u)'	:	[15],
							'd_xz(d)'	:	[16],
							'd_xz'		:	[15, 16],
							'd_x2-y2(u)':	[17],
							'd_x2-y2(d)':	[18],
							'd_x2-y2'	:	[17, 18],
							'd'			:	[9, 10, 11, 12, 13, 14, 15, 16, 17, 18],
							'z'			:	[5, 6, 11, 12, 13, 14, 15, 16]	
							}

		self.plot_ions = plot_ions
		self.plot_orbitals = plot_orbitals

		self.color = [
			'#DC143C', # 	crimson 			#DC143C 	(220,20,60)
			'#ADFF2F', #	green yellow 		#ADFF2F 	(173,255,47)
			'#40E0D0', #	turquoise 			#40E0D0 	(64,224,208)
			'#FF8C00', #  	dark orange 		#FF8C00 	(255,140,0)
			'#BA55D3', #	medium orchid 		#BA55D3 	(186,85,211)
			'#1E90FF', #	DOSCAR_FeTPyPFe_11_saddleger blue 		#1E90FF 	(30,144,255)
			'#FF1493', #	deep pink 			#FF1493 	(255,20,147)
			'#8B4513', #	saddle brown 		#8B4513 	(139,69,19)
			'#FFD700', #	gold 				#FFD700 	(255,215,0)
			'#808000', #	Olive 				#808000 	(128,128,0)
			'#808080', #	Gray 				#808080 	(128,128,128)
			'#FF00FF', #	Magenta / Fuchsia 	#FF00FF 	(255,0,255)
			'#00FFFF', #	Cyan / Aqua 		#00FFFF 	(0,255,255)
			'#000000', #	Black 				#000000 	(0,0,0)

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
		f = open(self.file_name,'r')

		self.n_E = 6
		self.E = [] 
		self.E_total = [] 
		for i, n in enumerate(f):
			vec = [ m for m in n.split(' ') if m != '']
			if i == 5: 	self.n_E = float(vec[2]) ; self.fermi = float(vec[3])
			
			if i > 5 and i < self.n_E+1+5: self.E_total.append([float(m) for m in vec])

			if i-5%(self.n_E+1) > 0 and int((i-5)/(self.n_E+1)) == 0: self.E_total.append([float(m) for m in vec])

			if int(i-5)%(self.n_E+1) == 0 and int((i-5)/(self.n_E+1)) > 0: self.E.append([])
			if int(i-5)%int(self.n_E+1) > 0 and int((i-5)/(self.n_E+1)) > 0: self.E[-1].append([float(m) for m in vec])

		self.E = np.array(self.E)
		self.n_oins, var, self.n_orb = self.E.shape
		print( f' >> Load complete :: data shape {self.E.shape} :: {file_name}')
		
		f.close()

		return self.E

	def export(self, file_name=None, fermi=True):
		if file_name == None: file_name = 'default'
		export_data = []

		if add_fermi: 	fermi = self.fermi 	# define fermi reference
		else: 			fermi = 0

		for i in ion:
			if i == 0: export_data.append( up_down_coef*self.E[0,:,o] )
			sum_acumulated = None
			for o in orbital:
				up_down_coef = o%2*2-1
				if sum_plot: # cumulated plot 
					if type(sum_acumulated) is np.ndarray:		sum_acumulated += up_down_coef*self.E[i,:,o]
					else: 										sum_acumulated =  up_down_coef*self.E[i,:,o]
				else:		export_data.append( self.E[i,:,0]-ferm )
			if sum_plot:	export_data.append( sum_acumulated )

		export_data = np.array(export_data)
		print(export_data.shape)
		return export_data

	def cut(self, ion=0, orbital=0, 
		start=-1.0, end=1.0, point=100, interpolation='cubic', fermi=True, endpoint=True):
		# ----------------------------------------------------------------------------------------------------- #
		# PLOT pDOS from DOSCAR result (VASP calculation)
		# ----------------------------------------------------------------------------------------------------- #
		# ion 			: 	LIST 	: [0< int <= ion number]  	eg. [1,2,3,4]	: 	selected ion to plot
		# orbital 		: 	LIST 	: [0< int <= 18]  			eg. [1,2,3,4]	: 	selected orbital to plot
		# start			: 	FLOAT 	: MIN value of the Energy interval  {start, end}
		# end 			:	FLOAT 	: END value of the Energy interval 	{start, end}
		# point			:	INT 	: number of interpolation points  	
		# interpolation	:	str		: type of interpolation || eg. 'linear', 'cubic' 
		# fermi 		:	BOOL	: add fermi level to Energy || eg. TRUE FALSE 
		# endpoint		:	BOOL	: include end point from interpolation || eg. TRUE FALSE
		# ----------------------------------------------------------------------------------------------------- #

		# set domine variable
		try: self.orbital_index
		except: 
			self.orbital_index = {
							'all'	 :[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18], 'z':	[5, 6, 11, 12, 13, 14, 15, 16],
							's(u)'	:	[1],	's(d)'	:	[2],	's'		:	[1,2],
							'py(u)'	:	[3],'py(d)'	:	[4],'py'	:	[3, 4],	'pz(u)'	:	[5],	'pz(d)'	:	[6],	'pz'	:	[5, 6],	'px(u)'	:	[7],	'px(d)'	:	[8],	'px'	:	[7, 8],	'p'		:	[3, 4, 5, 6, 7, 8],
							'd_xy(u)'	:	[9], 'd_xy(d)'	:	[10],	'd_xy'		:	[9, 10],	'd_yz(u)'	:	[11],	'd_yz(d)'	:	[12],	'd_yz'	 	:	[11, 12],	'd_z2-r2(u)':	[13],	'd_z2-r2(d)':	[14],	'd_z2-r2'	:	[13, 14],'d_xz(u)'	:	[15],	'd_xz(d)'	:	[16],	'd_xz'		:	[15, 16], 'd_x2-y2(u)':	[17],	'd_x2-y2(d)':	[18], 'd_x2-y2'	:	[17, 18],	'd':	[9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
							}

		if fermi: 		x = self.E[ion,:,0] - self.fermi
		else:	 		x = self.E[ion,:,0]

		orbital = orbital if type(orbital) == int else self.orbital_index[orbital][0] if type(orbital) == str else 1

		# set predictor variable
		y = self.E[ion,:,orbital]
		# train the predictor 
		f = interp1d(x, y, kind=interpolation)
		# make the prediction
		x0 = np.linspace(start, end, num=point, endpoint=True)
		# return result
		return f(x0)

	def plot(self, 	ion:list=None, orbital:list=None, ax=None, color=None, 
					positive_plot=False, sum_plot=False, add_fermi=True, legend_plot=True, spin_invertion=False, path=None, save=False):
		# ----------------------------------------------------------------------------------------------------- #
		# PLOT pDOS from DOSCAR result (VASP calculation)
		# ----------------------------------------------------------------------------------------------------- #
		# ion 		: 	LIST 	: [0< int <= ion number]  	eg. [1,2,3,4]	: 	selected ion to plot
		# orbital 	: 	LIST : [0< int <= 18]  			eg. [1,2,3,4]	: 	selected orbital to plot5
		# figure	: 	OBJ 	: matplot lib figure obj
		# color 	:	TUPLE 	: selected color 
		# sum_plot	:	BOOL	: Plot the total sum of all the orbitals 	
		# add_fermi	:	BOOL	: add fermi E to the energy 
		# legend_plot:	BOOL	: plot figure legend
		# ----------------------------------------------------------------------------------------------------- #

		try: 	self.orbital_index
		except: self.orbital_index = {
							'all'	 :[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18], 'z':	[5, 6, 11, 12, 13, 14, 15, 16],
							's(u)'	 :[1],  's(d)'	 :[2],  's'	  :[1,2],
							'py(u)'	: [3],'py(d)'	:	[4],'py'	:	[3, 4],	'pz(u)'	:	[5],	'pz(d)'	:	[6],	'pz'	:	[5, 6],	'px(u)'	:	[7],	'px(d)'	:	[8],	'px'	:	[7, 8],	'p'		:	[3, 4, 5, 6, 7, 8],
							'd_xy(u)':[9],  'd_xy(d)':[10],     'd_xy'      :[9, 10], 'd_yz(u)'   :[11],  'd_yz(d)':[12],      'd_yz':[11, 12],  'd_z2-r2(u)':[13],  'd_z2-r2(d)':[14],  'd_z2-r2':[13, 14],  'd_xz(u)':[15],	
							'd_xz(d)':[16],	'd_xz'	 :[15, 16], 'd_x2-y2(u)':[17],    'd_x2-y2(d)':[18],  'd_x2-y2':[17, 18],  'd'   :[9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
							}

		if ion != None and ion != []: pass
		elif self.plot_ions != None and self.plot_ions != []: ion = self.plot_ions
		else: ion = [0]

		if orbital != None and orbital != []: pass
		elif self.plot_orbitals != None and self.plot_orbitals != []: orbital = self.plot_orbitals
		else: orbital = [1,2]	

		orbital = orbital if type(orbital[0]) == int else [ m for n in orbital if n in self.orbital_index for m in self.orbital_index[n] ] if type(orbital[0]) == str else [1, 2] if type(orbital) == list and len(orbital) > 0 else [1, 2]
		orbital = np.array(orbital)


		if add_fermi: 	fermi = self.fermi 	# define fermi reference
		else: 			fermi = 0

		X = self.E[0,:,0]

		orbital_up = [ not True^spin_invertion	if o%2==1 else not False^spin_invertion for o in orbital]
		orbital_down = [ not False^spin_invertion	if o%2==1 else not True^spin_invertion for o in orbital]
		
		if ax is None: fig, ax = plt.subplots()

		for i in ion:
			for j, (o1, o2) in enumerate(zip(orbital[orbital_up], orbital[orbital_down])):
				if  sum_plot:
					E_up   =  np.sum(self.E[i,:,o1], axis=0)
					E_down = -np.sum(self.E[i,:,o2], axis=0)
					ax.plot(X-fermi, E_up.T,   label='UP',   color=color, alpha=0.8)
					ax.plot(X-fermi, E_down.T, label='DOWN', color=color, alpha=0.8 )

				elif positive_plot:
					E_up   =  np.sum(self.E[i,:,orbital[orbital_up]], axis=0)
					E_down =  np.sum(self.E[i,:,orbital[orbital_down]], axis=0)
					ax.plot(X-fermi, E_up.T + E_down.T,   label='UP + DOWN',   color=color, alpha=0.8)

				else:
					print( [self.orbital_name[o] for o in orbital[orbital_up] ] )
					E_up   =  self.E[i,:,o1]
					E_down = -self.E[i,:,o2]
					ax.plot(X-fermi, E_up.T,   color=self.color[j], label=f'{self.orbital_name[o1]}' )
					ax.plot(X-fermi, E_down.T, color=self.color[j], label=f'{self.orbital_name[o2]}' )
					
		Emax, Emin = np.max(E_up), np.min(E_down)   	# evaluate min/max in X #
		ax.plot( (0,0), (Emin, Emax), ls='--', color=[0.3,0.3,0.3], lw=2 , alpha=0.1 )
		if legend_plot:	legend = ax.legend(shadow=True, fontsize='x-large')
		ax.set_xlim((-4, 4))

		if save:
			fig = ax.get_figure()
			fig.savefig(f'{path}' , dpi=400, pad_inches=0.1, bbox_inches='tight', horizontalalignment='right') 

	@Logs.LogDecorator()
	def plot_orbital_diagram(self, 	ion:list=None, orbital:list=None, add_fermi:bool=True, 
									ax:object=None, save:bool=True, v:bool=True, path:str=None):
		try: 	self.orbital_index
		except: self.orbital_index = {
							'all'	 :[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],  'z':	[5, 6, 11, 12, 13, 14, 15, 16],
							's(u)'	 :[1],  's(d)'	 :[2],  's'	  :[1,2],
							'py(u)'	:	[3],'py(d)'	:	[4],'py'	:	[3, 4],	'pz(u)'	:	[5],	'pz(d)'	:	[6],	'pz'	:	[5, 6],	'px(u)'	:	[7],	'px(d)'	:	[8],	'px'	:	[7, 8],	'p'		:	[3, 4, 5, 6, 7, 8],
							'd_xy(u)':[9],  'd_xy(d)':[10],     'd_xy'      :[9, 10], 'd_yz(u)'   :[11],  'd_yz(d)':[12],      'd_yz':[11, 12],  'd_z2-r2(u)':[13],  'd_z2-r2(d)':[14],  'd_z2-r2':[13, 14],  'd_xz(u)':[15],	
							'd_xz(d)':[16],	'd_xz'	 :[15, 16], 'd_x2-y2(u)':[17],    'd_x2-y2(d)':[18],  'd_x2-y2':[17, 18],  'd'   :[9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
							} # dxy, dyz, dz2, dxz, dx2y2

		orbitals = orbital if type(orbital[0]) == int else [ m for n in orbital if n in self.orbital_index for m in self.orbital_index[n] ] if type(orbital[0]) == str else [1, 2] if type(orbital) == list and len(orbital) > 0 else [1, 2]
		orbitals = np.array(orbitals)
		if v: print( '-'*100 )

		if ax is None: fig, ax = plt.subplots()

		ax.set_xlim(-2, 5)
		ax.plot( (0, 3), (0, 0), ':', color=(0.4,0.4,0.4), alpha=0.5 )
		ax.set_ylabel('Energy (eV)')
		ax.set_title('')
		occupancy = self.get_occupancy(ion=ion, orbital=orbital, add_fermi=add_fermi, v=False)
		for i in ion:
			E, E_up, E_down = self.get_data(ion=[i], orbital=orbital, sum_plot=False, add_fermi=add_fermi )

			for j, (o, o_occ) in enumerate(occupancy[i].items()):

				print('[{:<10}(UP)  ] E:{:>8.3}eV, ocupated:{:>8.4} || Unocupated:{:>8.4} '.format(o,  E[np.argmax(E_up[j,:])],   o_occ['ocupated']['up'], o_occ['Unocupated']['up']) )
				print('[{:<10}(DOWN)] E:{:>8.3}eV, ocupated:{:>8.4} || Unocupated:{:>8.4} '.format(o,  E[np.argmin(E_down[j,:])], o_occ['ocupated']['down'], o_occ['Unocupated']['down']) )
				
				ax.plot( (0, 1), (E[np.argmin(E_down[j,:])], E[np.argmin(E_down[j,:])]), color=self.color[j] )
				ax.plot( (2, 3), (E[np.argmax(E_up[j,:])],   E[np.argmax(E_up[j,:])]), color=self.color[j] )
				ax.text(-1, E[np.argmin(E_down[j,:])], o,)
				ax.text( 4, E[np.argmax(E_up[j,:])], o,)	
					
			if save:
				fig = ax.get_figure()
				fig.savefig(f'{path}_{i}' , dpi=400, pad_inches=0.1, bbox_inches='tight', horizontalalignment='right') 

	def get_data(self, ion=None, orbital=None, sum_plot=False, add_fermi=True, spin_invertion=False, save=False):
		try: 	self.orbital_index
		except: self.orbital_index = {
							'all'	 :[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],  'z':	[5, 6, 11, 12, 13, 14, 15, 16],
							's(u)'	 :[1],  's(d)'	 :[2],  's'	  :[1,2],
							'py(u)'	:	[3],'py(d)'	:	[4],'py'	:	[3, 4],	'pz(u)'	:	[5],	'pz(d)'	:	[6],	'pz'	:	[5, 6],	'px(u)'	:	[7],	'px(d)'	:	[8],	'px'	:	[7, 8],	'p'		:	[3, 4, 5, 6, 7, 8],
							'd_xy(u)':[9],  'd_xy(d)':[10],     'd_xy'      :[9, 10], 'd_yz(u)'   :[11],  'd_yz(d)':[12],      'd_yz':[11, 12],  'd_z2-r2(u)':[13],  'd_z2-r2(d)':[14],  'd_z2-r2':[13, 14],  'd_xz(u)':[15],	
							'd_xz(d)':[16],	'd_xz'	 :[15, 16], 'd_x2-y2(u)':[17],    'd_x2-y2(d)':[18],  'd_x2-y2':[17, 18],  'd'   :[9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
							}

		if ion != None and ion != []: pass
		else: ion = [0]

		if orbital != None and orbital != []: pass
		else: orbital = [1,2]	

		orbital = orbital if type(orbital[0]) == int else [ m for n in orbital if n in self.orbital_index for m in self.orbital_index[n] ] if type(orbital[0]) == str else [1, 2] if type(orbital) == list and len(orbital) > 0 else [1, 2]
		orbital = np.array(orbital)

		E = self.E[0,:,0]-self.fermi if add_fermi else self.E[0,:,0]

		orbital_up = [ not True^spin_invertion	if o%2==1 else not False^spin_invertion for o in orbital]
		orbital_down = [ not False^spin_invertion	if o%2==1 else not True^spin_invertion for o in orbital]

		for i in ion:
			if  sum_plot:
				E_up   =  np.sum(self.E[i,:,orbital[orbital_up]], axis=0)
				E_down = -np.sum(self.E[i,:,orbital[orbital_down]], axis=0)
			else:
				E_up   =  self.E[i,:,orbital[orbital_up]]
				E_down = -self.E[i,:,orbital[orbital_down]]
				
		return E, E_up, E_down

	def get_center_of_band(self, ion=None, orbital=None, sum_plot=False, add_fermi=True, save=True):
		# === get data === # 
		E, E_up, E_down = self.get_data(ion=ion, orbital=orbital, sum_plot=sum_plot, add_fermi=add_fermi )

		# === get center_band === # 
		center_band = np.sum( (E_up-E_down) *E * (E[1]-E[0])) / np.sum((E_up-E_down) * (E[1]-E[0]) )
		
		# === save data into obj === # 
		if save: self.center_band = center_band
		
		return center_band 


	def get_band_dispercion(self, ion=None, orbital=None, sum_plot=False, add_fermi=True, save=True):
		# === get data === # 
		E, E_up, E_down = self.get_data(ion=ion, orbital=orbital, sum_plot=sum_plot, add_fermi=add_fermi )

		# === get center_band === # 
		center_band = np.sum( (E_up-E_down) *E * (E[1]-E[0])) / np.sum((E_up-E_down) * (E[1]-E[0]) )
		
		# === save data into obj === # 
		if save: self.center_band = center_band
		
		return center_band 

	@Logs.LogDecorator()
	def get_occupancy(self, ion:list=None, orbital:list=None,  add_fermi:bool=True, save:bool=True, v:bool=True, path:str=None):
		try: 	self.orbital_index
		except: self.orbital_index = {
							'all'	 :[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],  'z':	[5, 6, 11, 12, 13, 14, 15, 16],
							's(u)'	 :[1],  's(d)'	 :[2],  's'	  :[1,2],
							'py(u)'	:	[3],'py(d)'	:	[4],'py'	:	[3, 4],	'pz(u)'	:	[5],	'pz(d)'	:	[6],	'pz'	:	[5, 6],	'px(u)'	:	[7],	'px(d)'	:	[8],	'px'	:	[7, 8],	'p'		:	[3, 4, 5, 6, 7, 8],
							'd_xy(u)':[9],  'd_xy(d)':[10],     'd_xy'      :[9, 10], 'd_yz(u)'   :[11],  'd_yz(d)':[12],      'd_yz':[11, 12],  'd_z2-r2(u)':[13],  'd_z2-r2(d)':[14],  'd_z2-r2':[13, 14],  'd_xz(u)':[15],	
							'd_xz(d)':[16],	'd_xz'	 :[15, 16], 'd_x2-y2(u)':[17],    'd_x2-y2(d)':[18],  'd_x2-y2':[17, 18],  'd'   :[9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
							} # dxy, dyz, dz2, dxz, dx2y2

		orbitals = orbital if type(orbital[0]) == int else [ m for n in orbital if n in self.orbital_index for m in self.orbital_index[n] ] if type(orbital[0]) == str else [1, 2] if type(orbital) == list and len(orbital) > 0 else [1, 2]
		orbitals = np.array(orbitals)
		if v: print( '-'*100 )

		# === responce dict === # 
		occupancy = {}

		file = open(f'{path}_occupancy.dat', 'w') 
		for i in ion:
			# === responce dict === # 
			occupancy[i] = {}

			# === get data === # # === get data === #  # === get data === # 
			E, E_up, E_down = self.get_data(ion=[i], orbital=orbital, sum_plot=False, add_fermi=add_fermi )

			for j, (Eiu, Eid) in enumerate( zip(E_up, E_down) ):
				# === get occupancy === # # === get occupancy === #  # === get occupancy === # 

				# ***** OCUPATE states ***** # 
				E_ocupate				= E[E<0]
				# -- UP -- #
				state_occupancy_up   	= Eiu[E<0]
				state_occupancy_up_mean = np.sum(state_occupancy_up)*(E[1]-E[0])
				state_center_up			= np.sum( state_occupancy_up * E_ocupate * (E[1]-E[0])) / np.sum(state_occupancy_up * (E[1]-E[0]) )
				# -- DOWN -- #
				state_occupancy_down 		= Eid[E<0]
				state_occupancy_down_mean 	= np.sum(state_occupancy_down)*(E[1]-E[0])
				state_center_down			= np.sum( state_occupancy_down * E_ocupate * (E[1]-E[0])) / np.sum(state_occupancy_down * (E[1]-E[0]) )

				# ***** OCUPATE states ***** # 
				E_unocupate 				= E[E>0]
				# -- UP -- #
				state_unoccupancy_up    	= Eiu[E>0]
				state_unoccupancy_up_mean 	= np.sum(state_unoccupancy_up)*(E[1]-E[0])
				state_center_up				= np.sum( state_unoccupancy_up * E_unocupate * (E[1]-E[0])) / np.sum(state_unoccupancy_up * (E[1]-E[0]) )
				# -- DOWN -- #
				state_unoccupancy_down  	= Eid[E>0]
				state_unoccupancy_down_mean = np.sum(state_unoccupancy_down)*(E[1]-E[0])
				state_center_down			= np.sum( state_unoccupancy_down * E_unocupate * (E[1]-E[0])) / np.sum(state_unoccupancy_down * (E[1]-E[0]) )

				# === responce dict === # 
				occupancy[i][self.orbital_name[orbitals[j*2]][:-3]] = { 'ocupated':   {'up':state_occupancy_up_mean,   'down':state_occupancy_down_mean, 	'total':state_occupancy_up_mean-state_occupancy_down_mean},
																	 'Unocupated': {'up':state_unoccupancy_up_mean, 'down':state_unoccupancy_down_mean, 'total':state_unoccupancy_up_mean-state_unoccupancy_down_mean} }

				if v: print( '[{:<10}] Ocupated   (Up/down/total) : {:<5.3f}/{:<5.3f}/{:<5.3f} \n             Unocupated (Up/down/total) : {:<5.3f}/{:<5.3f}/{:<5.3f}\n'.format( self.orbital_name[orbitals[j*2]][:-3],
										state_occupancy_up_mean, state_occupancy_down_mean, state_occupancy_up_mean-state_occupancy_down_mean,
										state_unoccupancy_up_mean, state_unoccupancy_down_mean, state_unoccupancy_up_mean-state_unoccupancy_down_mean))
				
				# === save data into obj === # 
				if save: 
					file.write( '[{:<10}] Ocupated   (Up/down/total) : {:<5.3f}/{:<5.3f}/{:<5.3f} \n             Unocupated (Up/down/total) : {:<5.3f}/{:<5.3f}/{:<5.3f}\n'.format( self.orbital_name[orbitals[j*2]][:-3],
										state_occupancy_up_mean, state_occupancy_down_mean, state_occupancy_up_mean-state_occupancy_down_mean,
										state_unoccupancy_up_mean, state_unoccupancy_down_mean, state_unoccupancy_up_mean-state_unoccupancy_down_mean) )
			file.write( '-'*100+'\n' )

		if v: print( '-'*100 )
		file.close()

		return occupancy 
		
'''
doscar = DOSCAR()
doscar.load( file_name='/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/FePC+Au/DOSCAR' )

doscar.plot_orbital_diagram(ion=[0], orbital=['d'], add_fermi=True)
#doscar.plot( orbital=['d'] )
doscar.get_occupancy(ion=[0], orbital=['d'], add_fermi=True)


doscar.load( file_name='/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/CoPC+Au/DOSCAR' )
#doscar.get_occupancy(ion=[0], orbital=['d'], add_fermi=True)
doscar.plot_orbital_diagram(ion=[0], orbital=['d'], add_fermi=True)
plt.show()

'''



'''
doscar = DOSCAR()
doscar.load( file_name='/home/akaris/Documents/code/VASP/v4.6/files/PDOS/TPyP/catalisis/FeTPyPAu/DOSCAR' )
doscar.get_occupancy(ion=[64], orbital=['d'], add_fermi=True)


doscar.load( file_name='/home/akaris/Documents/code/VASP/v4.6/files/PDOS/TPyP/catalisis/FeTPyP[Co]Au/DOSCAR' )
doscar.get_occupancy(ion=[0], orbital=['d'], add_fermi=True)

doscar.load( file_name='/home/akaris/Documents/code/VASP/v4.6/files/PDOS/TPyP/catalisis/FeTPyPCoAu/DOSCAR' )
doscar.get_occupancy(ion=[1], orbital=['d'], add_fermi=True)


'''
'''
# FeTPyP[Co]
a = np.array([
# Up Ocupado, Down ocupado, Up libre, Down libre, 
		[0.921, -0.936, 0.0178, -0.0162], # dxy
		[0.129, -0.938, 0.832, -0.0386], # dyz
		[0.713, -0.923, 0.278, -0.0896], # dz2 
		[0.115, -0.94 , 0.847, -0.0381], # dxz
		[0.329, -0.437, 0.715, -0.612], ]) # dx2y2

# FeTPyPCo
b = np.array([
# Up Ocupado, Down ocupado, Up libre, Down libre, 
		[ -0.869, 0.925, -0.0556, 0.0139,], # dxy
		[ -0.458, 0.919,  -0.494, 0.045,], # dyz
		[ -0.208, 0.909, -0.763, 0.0917,], # dz2
		[ -0.403, 0.919, -0.549, 0.0456,], # dxz 
		[ -0.316, 0.407, -0.717, 0.632,], ])*-1 # dx2y2

# FeTPyP
c = np.array([
# Up Ocupado, Down ocupado, Up libre, Down libre, 
		[0.716, -0.764, 0.256, -0.219], # dxy
		[0.253, -0.934, 0.705, -0.039], # dyz 
		[0.671, -0.921, 0.314, -0.0863], # dz2 
		[0.23, -0.933, 0.731, -0.0437], # dxz
		[0.54, -0.615, 0.466, -0.399], ]) # dx2y2
 '''
# dxy, dyz, dz2, dxz, dx2y2
# Sin Co - Con Co

# dyz gana 0.33 
# dz2 pierde 0.5 
# dxz gana 0.30



def main(argv):
	# === organize arg === #
	inputfile  = argv['input']
	outputfile = argv['output']
	outputfile = 'PDOS.png' if type(outputfile) == type(None) else outputfile
	task 	   = argv['task']
	ion	   	   = argv['ion']
	orbital    = argv['orbital']

	# === Make data holder === #
	doscar = DOSCAR()

	if task == 'plot':
		figure, ax = plt.subplots(1)
		doscar.load( file_name=inputfile )
		ax.set_xlim(-4, 4)
		doscar.plot(ion=ion, orbital=orbital, ax=ax, color=(.8, .3, .3), 
							positive_plot=True, sum_plot=False, add_fermi=True, 
							legend_plot=True, save=True, path=f'{inputfile}_{orbital}_{ion}_PDOS.png')
		plt.show()

	if task == 'dat':
		doscar.load( file_name=inputfile )
		E, E_up, E_down = doscar.get_data(ion=ion, orbital=orbital, sum_plot=True, add_fermi=True )
		np.savetxt(f'{inputfile}_{orbital}_{ion}.dat', np.array([E, E_up, E_down, (E_up+np.abs(E_down))]).T, fmt='%10f')
	
	if task == 'band_center' or task == 'bc' :
		doscar.load( file_name=inputfile )
		center_band = doscar.get_center_of_band(ion=ion, orbital=orbital, sum_plot=True, add_fermi=True )
		print(f'>> Centro de banda: {center_band} \n Ion : {ion} \n Orbital : {orbital}')

	if task == 'analisys':
		figure, ax = plt.subplots(1)
		doscar.load( file_name=inputfile )
		center_band = doscar.get_center_of_band(ion=ion, orbital=orbital, sum_plot=True, add_fermi=True )

		doscar.plot(ion=ion, orbital=orbital, ax=ax, color=(.8, .3, .3), 
							positive_plot=False, sum_plot=False, add_fermi=True, 
							legend_plot=True, save=True, path=f'{inputfile}_{orbital}_{ion}_{center_band:.3}eV.png')

		doscar.plot_orbital_diagram(ion=ion, orbital=orbital, add_fermi=True, save=True, path=f'{inputfile}_orbital_diagram')

		E, E_up, E_down = doscar.get_data(ion=ion, orbital=orbital, sum_plot=True, add_fermi=True )
		
		doscar.get_occupancy(ion=ion, orbital=orbital, add_fermi=True, path=f'{inputfile}')

		np.savetxt(f'{inputfile}_{orbital}_{ion}_{center_band:.3}eV.dat', np.array([E, E_up, E_down, (E_up+np.abs(E_down))]).T, fmt='%10f')


		plt.show()

	print(f'Input  >> {inputfile} ')
	print(f'OUTPUT >> {outputfile}')


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-t','--task', help="task to accomplish \n   read  :  Read files from dataset  \n    summarise  :  resume data from dataset ",
	                    type=str, default='read', required=True)
	
	parser.add_argument('-i','--input', help="path to inputfile",
	                    type=str, default='', required=True)

	parser.add_argument('-o','--output', help="name of data output file",
	                    type=str, default=None, required=False)

	parser.add_argument('-ion','--ion', action='append', help="Number of the ion (starting from 1)",
	                    type=int, default=None, required=False)

	parser.add_argument('-orbital','--orbital', action='append', help="Orbital to plot",
	                    type=str, default=None, required=False)

	args = vars(parser.parse_args())
	main(args)
