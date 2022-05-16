import numpy as np 	
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
							'd':	[9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
							}

		self.plot_ions = plot_ions
		self.plot_orbitals = plot_orbitals

		self.plot_color = [
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
		print(self.E.shape)
		self.n_oins, var, self.n_orb = self.E.shape
		
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
							's(u)'	:	[1],	's(d)'	:	[2],	's'		:	[1,2],
							'px(u)'	:	[3],'px(d)'	:	[4],'px'	:	[3, 4],	'py(u)'	:	[5],	'py(d)'	:	[6],	'py'	:	[5, 6],	'pz(u)'	:	[7],	'pz(d)'	:	[8],	'pz'	:	[7, 8],	'p'		:	[3, 4, 5, 6, 7, 8],
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

	def plot(self, 	ion=None, orbital=None, figure=None, color=None, 
					sum_plot=False, add_fermi=True, legend_plot=True, spin_invertion=False, path=None, save=False):
		# ----------------------------------------------------------------------------------------------------- #
		# PLOT pDOS from DOSCAR result (VASP calculation)
		# ----------------------------------------------------------------------------------------------------- #
		# ion 		: 	LIST 	: [0< int <= ion number]  	eg. [1,2,3,4]	: 	selected ion to plot
		# orbital 	: 	LIST : [0< int <= 18]  			eg. [1,2,3,4]	: 	selected orbital to plot
		# figure	: 	OBJ 	: matplot lib figure obj
		# color 	:	TUPLE 	: selected color 
		# sum_plot	:	BOOL	: Plot the total sum of all the orbitals 	
		# add_fermi	:	BOOL	: add fermi E to the energy 
		# legend_plot:	BOOL	: plot figure legend
		# ----------------------------------------------------------------------------------------------------- #

		try: 	self.orbital_index
		except: self.orbital_index = {
							's(u)'	 :[1],  's(d)'	 :[2],  's'	  :[1,2],
							'px(u)'	 :[3],  'px(d)'  :[4],  'px'  :[3, 4],  'py(u)'  :[5],   'py(d)'  :[6],   'py'  :[5, 6],    'pz(u)'	:	[7],	'pz(d)'	:	[8],	'pz'	:	[7, 8],	'p'		:	[3, 4, 5, 6, 7, 8],
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
		#fig, ax = plt.subplots()
		for i in ion:
			if  sum_plot:
				E_up   =  np.sum(self.E[i,:,orbital[orbital_up]], axis=0)
				E_down = -np.sum(self.E[i,:,orbital[orbital_down]], axis=0)
				figure.plot(X-fermi, E_up.T,   label='UP',   color=color, alpha=0.8)
				figure.plot(X-fermi, E_down.T, label='DOWN', color=color, alpha=0.8 )
			else:
				E_up   =  self.E[i,:,orbital[orbital_up]]
				E_down = -self.E[i,:,orbital[orbital_down]]
				figure.plot(X-fermi, E_up.T, label='UP' )
				figure.plot(X-fermi, E_down.T, label='DOWN')
				
		Emax, Emin = np.max(E_up), np.min(E_down)   	# evaluate min/max in X #
		figure.plot( (0,0), (Emin, Emax), ls='--', color=[0.3,0.3,0.3], lw=2  )
		legend = ax.legend(shadow=True, fontsize='x-large')

		if save:
			fig = axs.get_figure()
			fig.savefig(f'{path}' , dpi=100, pad_inches=0.1, bbox_inches='tight', horizontalalignment='right') 

def main(argv):
	# === organize arg === #
	inputfile  = argv['input']
	outputfile = argv['output']
	outputfile = 'PDOS.png' if type(outputfile) == type(None) else outputfile
	task 	   = argv['task']
	ion	   = argv['ion']

	# === Make data holder === #
	dataset = Set()

	if task == 'd':
		figure, ax = plt.subplots(1)
		self.load( file_name=inputfile )
		self.plot(ion=[ion], orbital=['d' ], figure=ax, color=(.8, .3, .3), 
							sum_plot=True, add_fermi=True, legend_plot=True, save=True, path=outputfile)
		
	print(f'Input  >> {inputfile} ')
	print(f'OUTPUT >> {outputfile}')

'''
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	
	parser.add_argument('-t','--task', help="task to accomplish \n   read  :  Read files from dataset  \n    summarise  :  resume data from dataset ",
	                    type=str, default='read', required=True)
	
	parser.add_argument('-i','--input', help="path to inputfile",
	                    type=str, default='', required=True)

	parser.add_argument('-o','--output', help="name of data output file",
	                    type=str, default=None, required=False)

	parser.add_argument('-ion','--ion', help="Number of the ion (starting from 1)",
	                    type=str, default=None, required=False)

	args = vars(parser.parse_args())
	main(args)
'''

# --- eg. (load DOSCAR file)--- #
dos = DOSCAR()
figure, ax = plt.subplots(1)



dos.load( file_name='/home/akaris/Documents/code/VASP/v4.1/files/PDOS/FeTPyPFe/DOSCAR771' )
dos.plot(ion=[0], orbital=['d' ], figure=ax, color=(.9, .2, .9), 
					sum_plot=True, add_fermi=True, legend_plot=True)

dos.load( file_name='/home/akaris/Documents/code/VASP/v4.1/files/PDOS/FeTPyPFe/DOSCAR771' )
dos.plot(ion=[1], orbital=['d' ], figure=ax, color=(.9, .2, .2), 
					sum_plot=True, add_fermi=True, legend_plot=True)


plt.show()




dos.plot(ion=[0], orbital=['d' ], figure=ax, color=(.3, .3, .8), 
					sum_plot=True, add_fermi=True, legend_plot=True)

'''
plt.show()

figure, ax = plt.subplots(1)
figure, ax = plt.subplots(1)
dos.load( file_name='/home/akaris/Documents/code/VASP/v4.1/files/PDOS/FeTPyPFe/DOSCAR771' )
dos.plot(ion=[2,3,4,5], orbital=['d'], figure=ax, color=(.8, .3, .3), 
					sum_plot=True, add_fermi=True, legend_plot=True)


dos.plot(ion=[0], orbital=['d'], figure=ax, color=(.3, .3, .8), 
					sum_plot=True, add_fermi=True, legend_plot=True)
plt.show()


dos.load( file_name='/home/akaris/Documents/code/VASP/v4.1/files/PDOS/FeTPyPFe/DOSCAR111' )
figure, ax = plt.subplots(1)
dos.plot(ion=[0], orbital=['d'], figure=ax, color=(.3, .3, .8), 
					sum_plot=True, add_fermi=True, legend_plot=True)

dos.load( file_name='/home/akaris/Documents/code/VASP/v4.1/files/PDOS/FeTPyPFe/DOSCAR331' )
dos.plot(ion=[0], orbital=['d'], figure=ax, color=(.8, .3, .8), 
					sum_plot=True, add_fermi=True, legend_plot=True)

dos.load( file_name='/home/akaris/Documents/code/VASP/v4.1/files/PDOS/FeTPyPFe/DOSCAR551' )
dos.plot(ion=[0], orbital=['d'], figure=ax, color=(.8, .3, .3), 
					sum_plot=True, add_fermi=True, legend_plot=True)

dos.load( file_name='/home/akaris/Documents/code/VASP/v4.1/files/PDOS/FeTPyPFe/DOSCAR771' )
dos.plot(ion=[0], orbital=['d'], figure=ax, color=(.3, .8, .3), 
					sum_plot=True, add_fermi=True, legend_plot=True)

plt.show()


'''



'''

dos.load( file_name='files/PDOS/DOSCAR_FeTPyPFe_11_plannar' )
dos.plot(ion=[0], orbital=['d'], figure=figure, color=1, 
					sum_plot=True, add_fermi=True, legend_plot=True)

dos.load( file_name='files/PDOS/DOSCAR_FeTPyPFe_22_saddle' )
dos.plot(ion=[257], orbital=['d'], figure=figure, color=2, 
					sum_plot=True, add_fermi=True, legend_plot=True, spin_invertion=True)

plt.ylim([-3.3, 3.3])
plt.savefig( fname='files/PDOS/PDOS111122.png', dpi=1200 )

figure = plt.figure(2)

dos.load( file_name='files/PDOS/DOSCAR_FeTPyPFe_22_saddle' )
dos.plot(ion=[257], orbital=['d'], figure=figure, color=0, 
					sum_plot=True, add_fermi=True, legend_plot=True, spin_invertion=True)

dos.plot(ion=[256], orbital=['d'], figure=figure, color=1, 
					sum_plot=True, add_fermi=True, legend_plot=True)

dos.plot(ion=[258], orbital=['d'], figure=figure, color=2, 
					sum_plot=True, add_fermi=True, legend_plot=True, spin_invertion=True)

dos.plot(ion=[259], orbital=['d'], figure=figure, color=3, 
					sum_plot=True, add_fermi=True, legend_plot=True)

plt.ylim([-3.3, 3.3])
plt.savefig( fname='files/PDOS/PDOS22saddle.png', dpi=1200 )

figure = plt.figure(3)

dos.load( file_name='files/PDOS/DOSCAR_FeTPyPFe_22_propeller' )
dos.plot(ion=[0], orbital=['d'], figure=figure, color=0, 
					sum_plot=True, add_fermi=True, legend_plot=True, spin_invertion=True)

dos.plot(ion=[1], orbital=['d'], figure=figure, color=1, 
					sum_plot=True, add_fermi=True, legend_plot=True, spin_invertion=True)

dos.plot(ion=[2], orbital=['d'], figure=figure, color=2, 
					sum_plot=True, add_fermi=True, legend_plot=True)

dos.plot(ion=[3], orbital=['d'], figure=figure, color=3, 
					sum_plot=True, add_fermi=True, legend_plot=True)

plt.ylim([-3.3, 3.3])
plt.savefig( fname='files/PDOS/PDOS22propeller.png', dpi=1200 )

figure = plt.figure(4)

dos.load( file_name='files/PDOS/DOSCAR_FeTPyPFe_22_propeller' )
dos.plot(ion=[0], orbital=['d'], figure=figure, color=0, 
					sum_plot=True, add_fermi=True, legend_plot=True, spin_invertion=True)

dos.plot(ion=[1], orbital=['d'], figure=figure, color=1, 
					sum_plot=True, add_fermi=True, legend_plot=True, spin_invertion=True)

dos.load( file_name='files/PDOS/DOSCAR_FeTPyPFe_11_saddle' )
dos.plot(ion=[64], orbital=['d'], figure=figure, color=0, 
					sum_plot=True, add_fermi=True, legend_plot=True)

dos.load( file_name='files/PDOS/DOSCAR_FeTPyPFe_11_plannar' )
dos.plot(ion=[0], orbital=['d'], figure=figure, color=1, 
					sum_plot=True, add_fermi=True, legend_plot=True)

plt.ylim([-3.3, 3.3])
plt.savefig( fname='files/PDOS/PDOS111122propeller.png', dpi=1200 )

plt.show()
'''

# --- // --- #

'''
dos = DOSCAR()
dos.load( file_name='files/PDOS/DOSCAR' )

dos.plot(ion=[64], orbital=[17,18],legend_plot=True)
plt.show()
'''

'''
path_be = '/home/akaris/Documents/code/VASP/v3.5/files/PDOS/FePC+BE+Au'
path_au = '/home/akaris/Documents/code/VASP/v3.5/files/PDOS/FePC+Au'

fig, ax = plt.subplots()
dos = DOSCAR()
dos.load( file_name='{}/DOSCAR'.format(path_be) )
dos.plot(ion=[0], orbital=['d'], sum_plot=True, legend_plot=True,spin_invertion=True,figure=ax)

dos1 = DOSCAR()
dos1.load( file_name='{}/DOSCAR'.format(path_au) )
dos1.plot(ion=[0], orbital=['d'], sum_plot=True, legend_plot=True,spin_invertion=False,figure=ax)
plt.show()
'''