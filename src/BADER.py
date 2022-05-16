# *** warning supresion
import warnings
warnings.filterwarnings("ignore")

# *** numeric libraries *** #
import numpy as np
import scipy.io

# *** graph libraries *** #
try:
	import matplotlib.pyplot as plt
	from mpl_toolkits import mplot3d
	import matplotlib as mpl
except: 
	print('WARNNING :: main_simulation.py :: can NOT correctly load "matplotlib" libraries')
	print('Install by: ( pip3 install matplotlib )')
	

# *** python common libraries
import logging, operator, pickle, os

# load own libraries #
try:
	from src import DOSCAR, PROCAR, POSCAR, ORR
except: print('Can NOT load VASP libraries')
from src import DOSCAR, PROCAR, POSCAR, ORR

'''
http://theory.cm.utexas.edu/henkelman/code/bader/

Code: Bader Charge Analysis
News
02/08/20 - Version 1.04 Released
Bug fix for printing Bader volumes.
Introduction
Richard Bader, from McMaster University, developed an intuitive way of dividing molecules into atoms. 
His definition of an atom is based purely on the electronic charge density. Bader uses what are called 
zero flux surfaces to divide atoms. A zero flux surface is a 2-D surface on which the charge density 
is a minimum perpendicular to the surface. Typically in molecular systems, the charge density reaches a 
minimum between atoms and this is a natural place to separate atoms from each other.

Besides being an intuitive scheme for visualizing atoms in molecules, Bader's definition is often useful 
for charge analysis. For example, the charge enclosed within the Bader volume is a good approximation to 
the total electronic charge of an atom. The charge distribution can be used to determine multipole moments 
of interacting atoms or molecules. Bader's analysis has also been used to define the hardness of atoms, 
which can be used to quantify the cost of removing charge from an atom.

Program Overview
We have developed a fast algorithm for doing Bader's analysis on a charge density grid. The program 
(see below) can read in charge densities in the VASP CHGCAR format, or the Gaussian CUBE format. The 
program outputs the total charge associated with each atom, and the zero flux surfaces defining the Bader 
volumes.

Note for VASP users
One major issue with the charge density (CHGCAR) files from the VASP code is that they only contain the 
valance charge density. The Bader analysis assumes that charge density maxima are located at atomic centers 
(or at pseudoatoms). Aggressive pseudopotentials remove charge from atomic centers where it is both expensive 
to calculate and irrelevant for the important bonding properties of atoms.
VASP contains a module (aedens) which allows for the core charge to be written out from PAW calculations. 
This module is included in vasp version 4.6.31 08Feb07 and later. By adding the LAECHG=.TRUE. to the INCAR file, 
the core charge is written to AECCAR0 and the valance charge to AECCAR2. These two charge density files can be 
summed using the chgsum.pl script:

  chgsum.pl AECCAR0 AECCAR2
The total charge will be written to CHGCAR_sum.
The bader analysis can then be done on this total charge density file:

  bader CHGCAR -ref CHGCAR_sum
One finally note is that you need a fine fft grid to accurately reproduce the correct total core charge. It is 
essential to do a few calculations, increasing NG(X,Y,Z)F until the total charge is correct.

'''

def cookbook():
	m1 = np.loadtxt( 'files/BADER/FePC/ACF.dat' )
	m2 = np.loadtxt( 'files/BADER/FePC+Au/ACF.dat' )
	m1m2 = m1[:57,4] - m2[:57,4]
	print( 'FePC (total): {0} , FePC+Au (total sin Au): {1} '.format(np.sum(m1[:57,4]), np.sum(m2[:57,4])),)
	print( 'FePC (total) - FePC+Au (total sin Au): {0} '.format(np.sum(m1[:57,4]) - np.sum(m2[:57,4])) )

	m3 = np.loadtxt( 'files/BADER/CoPC/ACF.dat' )
	m4 = np.loadtxt( 'files/BADER/CoPC+Au/ACF.dat' )
	m3m4 = m3[:57,4] - m4[:57,4]
	print( 'CoPC (total): {0} , CoPC+Au (total sin Au): {1} '.format(np.sum(m3[:57,4]), np.sum(m4[:57,4])),)
	print( 'CoPC (total) - CoPC+Au (total sin Au): {0} '.format(np.abs(np.sum(m3[:57,4])) - np.abs(np.sum(m4[:57,4]))),)

	print(np.abs(np.sum(m1m2[:])) - np.abs(np.sum(m3m4[:])) )
	print(np.abs(np.sum(np.abs(m1m2[:]))) - np.abs(np.sum(np.abs(m3m4[:]))) )
	print( (np.sum(m1[:57,4])-np.sum(m2[:57,4])) + (np.sum(m3[:57,4])-np.sum(m4[:57,4])) ) 
	error
	for i, n in enumerate(['Co']*1+['N']*8+['C']*32+['H']*16):
		print('Atom {}_{} ((FePC   )-(FePC+Au)) - ((CoPC   )-(CoPC+Au)) = {:.4f}'.format(n, i, m1m2[i]-m3m4[i]) )
		print('Atom {}_{} ( {:.4f}  - {:.4f}  ) - ({:.4f}   - {:.4f}  ) = {:.4f}'.format(n, i, m1[i,4], m2[i,4], m3[i,4], m4[i,4], m1m2[i]-m3m4[i]) )
		print('Atom {}_{} ( {:.4f}             -   {:.4f}  )            = {:.4f}'.format(n, i, m1m2[i], m3m4[i], m1m2[i]-m3m4[i]) )

	plt.figure(5)
	plt.plot( m1[:57,4] , '-o')
	plt.plot( m2[:57,4] , '-o')
	plt.title('FePC (blue), FePC+Au(orange)')
	plt.xlabel('Atom number')
	plt.ylabel('Valence charge')

	plt.figure(6)
	plt.plot( m3[:57,4] , '-o')
	plt.plot( m4[:57,4] , '-o')
	plt.title('CoPC (blue), CoPC+Au(orange)')
	plt.xlabel('Atom number')
	plt.ylabel('Valence charge')


	plt.figure(1)
	plt.plot( m1m2[:] , 'r-o')
	plt.plot( m3m4[:] , '-o', c=(0,0,0))
	plt.title(' FePC-FePC+Au(red), CoPC-CoPC+Au(black)')
	plt.xlabel('Atom number')
	plt.ylabel('Valence charge')

	plt.figure(2)
	plt.plot( m1[:57,4], m2[:57,4], 'o', )
	plt.plot( m3[:57,4], m4[:57,4], 'o')
	plt.title(' FePC diff system (red), CoPC diff system (black)')
	plt.xlabel('[FePC or CoPC] valence charge')
	plt.ylabel('[FePC+Au or CoPC+Au] valence charge')


	plt.figure(3)
	plt.plot( m1m2[:]-m3m4[:], '-o')
	plt.title(' Systems diff system ((FePC)-(FePC+Au)) - ((CoPC)-(CoPC+Au))')
	plt.xlabel('Atom number')
	plt.ylabel('((FePC)-(FePC+Au)) - ((CoPC)-(CoPC+Au))')

	plt.show()