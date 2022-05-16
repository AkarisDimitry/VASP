import numpy as np 	
import matplotlib.pyplot as plt

class EIGENVAL(object):
	def __init__(self, name=None):
		self.name = name
		self.n_electrons = None
		self.n_kpoints = None
		self.n_bands = None
		self.bands = None
		self.kpoints = None
		self.k_distance = None

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
							
	def load(self, name=None):
		if name != None: self.name = name
		f_eigenval = open(self.name,'r')

		self.bands, self.kpoints=[], []
		var = -1
		for i, n in enumerate(f_eigenval):

			vec = [float(m) for m in n.split(' ') if self.isnum(m) ] 
			if i == 5: self.bands = np.zeros((int(vec[1]), int(vec[2]), 4))
			if len(vec) == 4 and i>5: self.kpoints.append(vec); var+=1
			if len(vec) == 5 and i>5: self.bands[var,int(vec[0]-1), :] = vec[1:]

		self.k_distance = np.zeros((len(self.kpoints)))
		var = 0
		for n in range(len(self.k_distance)-1): 
			var += ((self.kpoints[n][0]-self.kpoints[n+1][0])**2+(self.kpoints[n][1]-self.kpoints[n+1][1])**2+(self.kpoints[n][2]-self.kpoints[n+1][2])**2)**0.5
			self.k_distance[n+1] = var

	def isnum(self, n):
		# ------------------ Define if n is or not a number ------------------ # 
		# n     :   VAR     :   VAR to check if it is a numerical VAR
		# return :  BOOL    : True/False
		try: float(n); return True
		except: return False

	def summary(self, ):
		pass

	def plot(self, ):
		plt.plot(self.k_distance, self.bands[:, :, 0])
		plt.show()

# How to .. 
#EV = EIGENVAL('EIGENVAL')
#EV.load()
#EV.plot()
 # of electrons, # of k-points, #of bands 