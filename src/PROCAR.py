import numpy as np 	
import matplotlib.pyplot as plt

class PROCAR(object):
	def __init__(self, name=None, plot_ions=None, plot_spins=None, plot_orbitals=None, ):
		self.name = name
		self.file_name = None
		self.n_electrons = None
		self.KPOINTS = None # kpoint
		self.weight = None # kpoint

		self.E = None # [Kpoints, band] --> Energy
		self.occupation = None # [Kpoints, band] --> Ocupation
		self.P = None # [Kpoint, bands, ions, orbital] --> <Y|Phi>

		self.n_kpoints = None
		self.n_band = None
		self.n_oins = None

		self.orbitals = {'$s$':0, '$px$':1, '$py$':2, '$pz$':3,
					r'$dxy$':4, r'$dyz':5, r'$dz^{2}$':6, r'$dxz$':7, r'$dx^{2}$':8, r'tot':9 } #    s     py     pz     px    dxy    dyz    dz2    dxz    dx2

		self.plot_ions = plot_ions
		self.plot_spins = plot_spins
		self.plot_orbitals = plot_orbitals

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

	def isnum(self, n):
		# ------------------ Define if n is or not a number ------------------ # 
		# n     :   VAR     :   VAR to check if it is a numerical VAR
		# return :  BOOL    : True/False
		try: float(n); return True
		except: return False

	def load(self, file_name=None):
		if file_name != None: self.file_name = file_name
		f = open(self.file_name,'r')

		self.KPOINTS = [ ] 
		self.E = [ ]
		self.occupation = []
		self.P = [ ]

		for i, n in enumerate(f):
			vec = [ m for m in n.split(' ') if m != '']

			if vec[0] == 'k-point': self.KPOINTS.append([vec[3],vec[4],vec[5]]); self.E.append([]); self.occupation.append([]); self.P.append([])
			if vec[0] == 'band': self.E[-1].append(float(vec[4])); self.P[-1].append([]) 
			if self.isnum(vec[0]): self.P[-1][-1].append( [float(m) for i, m in enumerate(vec) if i>0] ) 

		self.KPOINTS = np.array(self.KPOINTS)  ;  self.E = np.array(self.E) 
		self.occupation = np.array(self.occupation) ; self.P = np.array(self.P) 

		self.n_kpoints = self.E.shape[0]/2 ; self.n_band = self.E.shape[1] ; self.n_oins = self.P.shape[2]

	def plot(self, plot='band', parameters={} ):

		if  plot == 'band':
			parameters = { key:value if not key in parameters else parameters[key] for key, value in {'spin':None, 'ion':None, 'orbital':None, 'figure':None}.items() }
			self.plot_band( spin=parameters['spin'], ion=parameters['ion'], orbital=parameters['orbital'], figure=parameters['figure'] )
		else:
			print(' WARNNING :: PROCAR.plot() :: Can not identify plot type *{}*'.format(plot) )
		return None

	def plot_band(self, spin=None, ion=None, orbital=None, figure=None):
		#if not self.isnum(orbital): orbital=self.orbitals[orbital]
		if ion != None and ion != []: pass
		elif self.plot_ions != None and self.plot_ions != []: ion = self.plot_ions
		else: ion = [0]

		if spin != None and spin != []: pass
		elif self.plot_spins != None and self.plot_spins != []: spin = self.plot_spins
		else: spin = [0]

		if orbital != None and orbital != []: pass
		elif self.plot_orbitals != None and self.plot_orbitals != []: orbital = self.plot_orbitals
		else: orbital = [0]	

		if not figure == None: figure = figure.add_subplot(111)
		 
		for n in range(self.n_band):
			C = np.array([1.0, 0.0, 0.0])

			for s_i, s in enumerate(spin):
				for i_i, i in enumerate(ion):
					for o_i, o in enumerate(orbital):
						if [s_i,i_i,o_i] != [0,0,0]: P += np.array(self.P[int(s*self.n_kpoints):int((s+1)*self.n_kpoints), n, i, o])
						else: P = np.array(self.P[int(s*self.n_kpoints):int((s+1)*self.n_kpoints), n, i, o])
			try:P = np.where(P < 1, P, 1.0)
			except: pass
			
			if figure == None:
				print(self.P.shape)
				print(self.n_kpoints*2)

				plt.figure(1), plt.scatter(np.arange(0,self.n_kpoints*2,1), self.E[:,n] ) # s=P*50		, color=np.dot(np.array([P]).T, np.array([C]))
				plt.figure(1), plt.plot( self.E[:,n], ls='-', markersize=2, color=np.mean(P)*C )
			else: 
				figure.scatter(np.arange(0,30,1), self.E[:30,n], s=P*50, color=np.dot(np.array([P]).T, np.array([C])) )
				figure.plot( self.E[:30,n], ls='-', markersize=2, color=np.mean(P)*C )

	def plot_E(self, spin=None, ion=None, orbital=None, figure=None, place=None):
		#if not self.isnum(orbital): orbital=self.orbitals[orbital]
		if ion != None and ion != []: pass
		elif self.plot_ions != None and self.plot_ions != []: ion = self.plot_ions
		else: ion = [0]

		if spin != None and spin != []: pass
		elif self.plot_spins != None and self.plot_spins != []: spin = self.plot_spins
		else: spin = [0]

		if orbital != None and orbital != []: pass
		elif self.plot_orbitals != None and self.plot_orbitals != []: orbital = self.plot_orbitals
		else: orbital = [0]	

		if place == None: place = 0 
		else: pass

		if not figure == None: figure = figure.add_subplot(111)
		 
		C1 = np.array([0.8,0.3,0.3])
		C2 = np.array([0.2,0.2,0.2])

		dic = {0:'s', 1:'px', 2:'py', 3:'pz',
					4:'dxy', 5:'dyz', 6:'dz2', 7:'dxz', 8:'dx2', 9:'tot' }
		for n in range(self.P.shape[1]):
			if figure == None: 
				color = (C1-C2)*self.P[0,n,0,9]/np.max(self.P[0,:,0,9]) + C2
				#plt.figure(1), plt.plot((0,1), (self.E[0,n],self.E[0,n]), lw=2, color='#888888' ) 
				plt.figure(1), plt.plot((place+0,place+1), (self.E[0,n],self.E[0,n]), lw=2, color=color ) 
				text = ''
				for o in range(9):
					if self.P[0, n, 0, o] > 0.1:
						text += " {} = {} ||".format(dic[o], self.P[0, n, 0, o])
					if len(text)>0: plt.text(place, self.E[0,n], text,
			    			        	color=color, weight='heavy', 
		    	        				size=16, )
			else: 
				color = (C1-C2)*self.P[0,n,0,9]/np.max(self.P[0,:,0,9]) + C2
				#plt.figure(1), plt.plot((0,1), (self.E[0,n],self.E[0,n]), lw=2, color='#888888' ) 
				figure.plot((place+0,place+1), (self.E[0,n],self.E[0,n]), lw=2, color=color ) 
				text = ''
				for o in range(9):
					if self.P[0, n, 0, o] > 0.1:
						text += " {} = {} ||".format(dic[o], self.P[0, n, 0, o])
					if len(text)>0: figure.text(place, self.E[0,n], text,
			    			        	color=color, weight='heavy', 
		    	        				size=16, )

	def summary(self, vebosity=1):
		print( '*'*10+' PROCAR_summary '+'*'*10)
		print('KPOINTs number: {}'.format(self.n_kpoints) )
		print('BANDs number: {}'.format(self.n_band) )
		print('IONs number: {}'.format(self.n_oins) )
		#print('ORBITALs number: {}'.format(self.n_kpoints) )

		dic = {	0:r'$s$', 1:r'$p_{x}$', 2:r'$p_{z}$', 3:r'$p_{y}$',
				4:r'$d_{xy}$', 5:r'$d_{yz}$', 6:r'$d_{z^{2}}$', 7:r'$d_{xz}$', 8:r'$d_{x^{2}}$', 9:r'tot' }

		if vebosity > 0:
			print('[self.E] = [{} (KPOINTS_UP+KPOINTS_down), {} (bands)]'.format(self.n_kpoints, self.n_band) )
			print('[self.P] = [{} (KPOINTS_UP+KPOINTS_down), {} (bands), {} (ions), {} (orbitals+1)]'.format(self.n_kpoints, self.n_band, self.n_oins, 10))

		for spin in range(self.P.shape[0]):
			print('****** SPIN {} ******'.format(spin) )
			for molecular_orbital in range(self.P.shape[1]):
				tot =  np.sum(self.P[spin, molecular_orbital, :, :], axis=0)
				print(r'\begin{tabular}{ |p{2cm}||p{2cm}|p{2cm}|p{2cm}|p{2cm}|p{2cm}|  }')
				print(r'\hline')
				print( r'\multicolumn{4}{|c|}{'+r' MOLECULAR ORBITAL {} S:{:.1f}\% P:{:.1f}\% D:{:.1f}\% :: {:.3f}  '.format(molecular_orbital+1, 100*tot[0]/tot[9], 100*(tot[1]+tot[2]+tot[3])/tot[9], 100*(tot[4]+tot[5]+tot[6]+tot[7]+tot[8])/tot[9], tot[9]) + r'} \\' )
				print(r'\hline')
				for atom in range(self.P.shape[2]):
					if self.P[spin, molecular_orbital, atom, 9] > 0.01:
						str_info = 'ATOM {}: '.format(atom)

						for atomic_orbital in range(self.P.shape[3]):
							if  100*self.P[spin, molecular_orbital, atom, atomic_orbital]/tot[9] > 0.2 and atomic_orbital != 9:
								str_info += r'& {}:{:.1f}\% '.format( dic[atomic_orbital],  100*self.P[spin, molecular_orbital, atom, atomic_orbital]/tot[9] )
						print(str_info+r' \\')	
				print(r'\hline')
				print(r'\end{tabular}')

# How to .. 
# of electrons, # of k-points, #of bands 

'''
PC = PROCAR() 
PC.load('/home/akaris/Documents/code/VASP/v4.1/files/PDOS/FeTPyPFe/PROCAR_111')
PC.plot()
plt.show()
PC.summary()
'''




