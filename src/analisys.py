import numpy as np
import matplotlib.pyplot as plt

names = ['CoPC', 'CoPC_Au (HOLLOW)', 'CoPC_Au (TOP)', 'CoPC_Au (BRIDGE)', 'CoPC_Bz_Au (Au)', 'CoPC_Bz_Au (Bz)', 'CoPC_Bz_Au (C)',
		'FePC', 'FePC_Au (HOLLOW)', 'FePC_Au (TOP)', 'FePC_Au (BRIDGE)', 'FePC_Bz_Au (Au)', 'FePC_Bz_Au (Bz)', 'FePC_Bz_Au (C)' ]




import functional_base as fb

#doscar = DOSCAR()
#doscar.load( file_name='/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/CoPC_Au/DOSCAR_CoPCAu_HOLLOW' )
#doscar.load( file_name='/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/CoPC_Bz_Au/DOSCAR_CoPCBzAu_C' )
#E, E_up, E_down = doscar.get_data(ion=[0], orbital=['d'], sum_plot=True, add_fermi=True )

path_list = [
	'/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/CoPC/DOSCAR_CoPC_[\'z\']_[0]_-0.934eV.dat',
	
	'/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/CoPC_Au/DOSCAR_CoPCAu_BRIDGE_[\'z\']_[0]_-0.382eV.dat',
	'/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/CoPC_Au/DOSCAR_CoPCAu_HOLLOW_[\'z\']_[0]_-0.508eV.dat',
	'/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/CoPC_Au/DOSCAR_CoPCAu_TOP_[\'z\']_[0]_-0.404eV.dat',
	
	'/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/CoPC_Bz_Au/DOSCAR_CoPCBzAu_Au_[\'z\']_[0]_-0.704eV.dat',
	'/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/CoPC_Bz_Au/DOSCAR_CoPCBzAu_Be_[\'z\']_[0]_-0.749eV.dat',
	'/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/CoPC_Bz_Au/DOSCAR_CoPCBzAu_C_[\'z\']_[0]_-0.837eV.dat',

	'/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/FePC/DOSCAR_FePC_[\'z\']_[0]_-0.946eV.dat',
	
	'/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/FePC_Au/DOSCAR_FePCAu_BRIDGE_[\'z\']_[0]_-0.667eV.dat',
	'/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/FePC_Au/DOSCAR_FePCAu_HOLLOW_[\'z\']_[0]_-0.645eV.dat',

	
	'/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/FePC_Bz_Au/DOSCAR_FePCBzAu_Au_[\'z\']_[0]_-0.635eV.dat',
	'/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/FePC_Bz_Au/DOSCAR_FePCBzAu_Be_[\'z\']_[0]_-0.525eV.dat',
	'/home/akaris/Documents/code/VASP/v4.6/files/PDOS/PC/Catalisis/FePC_Bz_Au/DOSCAR_FePCBzAu_C_[\'z\']_[0]_-0.637eV.dat',
	]

def DOSCAR2PEM():
	for i, n in enumerate(path_list):
		path_particular = '/'.join( n.split('/')[:-1])
		data = np.loadtxt(n)

		base = fb.BASE() 
		base.PEM(data=data[:,1], verbosity=True, iterations=30, sigma=np.arange(10, 50, 0.5))

		#base.plot_PEM( )
		np.savetxt( f'{path_particular}/PEM_coef_{i}.dat', base.get_PEM_coef() )

#DOSCAR2PEM()

OP = [  #-0.934,  -0.382, -0.508, -0.404, -0.837, -0.749, -0.704, 
		#-0.946, -0.667, -0.645,  -0.645, -0.525, -0.525, -0.525]
		-1.18,  -0.756, -0.736, -0.740, -

BandX = [  0.509, 0.660, 0.634, 0.645, 0.558, 0.533,  0.525,
 		   0.551, 0.445,  0.487, 0.444, 0.542, 0.560, 0.574]

plt.plot(OP, BandX, 'o')
plt.show()
asdf

data_mean = []
for i, n in enumerate(path_list):
	path_particular = '/'.join( n.split('/')[:-1]) + f'/PEM_coef_{i}.dat'
	data = np.loadtxt(path_particular)
	print( np.mean(data, axis=0) )
	data_mean.append( np.mean(data, axis=0) )

print(data_mean)
data_mean = np.array(data_mean)
plt.plot(OP[:7], data_mean[:7, 2], 'o')
plt.plot(OP[7:], data_mean[7:, 2], 'o')
plt.show()


OP = [  -0.934, -0.508, -0.404, -0.382, -0.837, -0.749, -0.704, 
 		-0.946, -0.645,         -0.667, -0.525, -0.525, -0.525] # falta FePCAuTOP

BandX = [  0.509, 0.660, 0.634, 0.645, 0.558, 0.533,  0.525,
 		   0.551, 0.445,         0.444, 0.542, 0.560, 0.574]

plt.plot(OP[:7], BandX[:7], 'o')
plt.plot(OP[7:], BandX[7:], 'o')
plt.show()

asfd

