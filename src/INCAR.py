import numpy as np 	
import matplotlib.pyplot as plt

class INTCAR(object):
	def __init__(self, name=None, plot_ions=None, plot_orbitals=None):
		self.name = name
		self.file_name = None

# "Here is a short overview of all parameters currently supported. Parameters which are used frequently are emphasized.
		self.help = {
'NGX':	'FFT mesh for orbitals (Sec. 6.3,6.11)',
'NGY': 	'FFT mesh for orbitals (Sec. 6.3,6.11)',
'NGZ': 	'FFT mesh for orbitals (Sec. 6.3,6.11)',
'NGXF':	'FFT mesh for charges (Sec. 6.3,6.11)',
'NGYF':	'FFT mesh for charges (Sec. 6.3,6.11)',
'NGZF': 	'FFT mesh for charges (Sec. 6.3,6.11)',
'NBANDS': 	'number of bands included in the calculation (Sec. 6.5)',
'NBLK': 	'blocking for some BLAS calls (Sec. 6.6)',
'SYSTEM': 	'name of System',
'NWRITE': 	'verbosity write-flag (how much is written)',
'ISTART': 	'startjob: 0-new 1-cont 2-samecut',
'ICHARG': 	'charge: 1-file 2-atom 10-const',
'ISPIN' :	'spin polarized calculation (2-yes 1-no)',
'MAGMOM': 	'initial mag moment / atom',
'INIWAV': 	'initial electr wf. : 0-lowe 1-rand',
'ENCUT' :	'energy cutoff in eV',
'PREC' :	'precession: medium, high or low : VASP.4.5 also: normal, accurate',
'NELM': 		'nr. of electronic steps',
'NELMIN': 	'nr. of electronic steps',
'NELMDL': 	'nr. of electronic steps',
'EDIFF': 	'stopping-criterion for electronic upd.',
'EDIFFG': 	'stopping-criterion for ionic upd.',
'NSW': 	'number of steps for ionic upd.',
'NBLOCK': 'number of steps for ionic upd.',
'KBLOCK':	'inner block; outer block',
'IBRION': 	'ionic relaxation: 0-MD 1-quasi-New 2-CG',
'ISIF': 	'calculate stress and what to relax',
'IWAVPR': 	'prediction of wf.: 0-non 1-charg 2-wave 3-comb',
'ISYM' :	'symmetry: 0-nonsym 1-usesym',
'SYMPREC': 	'precession in symmetry routines',
'LCORR': 	'Harris-correction to forces',
'POTIM': 	'time-step for ion-motion (fs)',
'TEBEG':  'temperature during run',
'TEEND': 	'temperature during run',
'SMASS': 	'Nose mass-parameter (am)',
'NPACO': 	'distance and nr. of slots for P.C.',
'APACO': 	'distance and nr. of slots for P.C.',
'POMASS': 	'mass of ions in am',
'ZVAL': 	'ionic valence',
'RWIGS': 	'Wigner-Seitz radii',
'NELECT': 	'total number of electrons',
'NUPDOWN': 	'fix spin moment to specified value',
'EMIN':	'energy-range for DOSCAR file',
'EMAX': 	'energy-range for DOSCAR file',
'ISMEAR': 'part. occupancies: -5 Blöchl -4-tet -1-fermi 0-gaus >0 MP',
'SIGMA':  'broadening in eV -4-tet -1-fermi 0-gaus',
'ALGO': 	'algorithm: Normal (Davidson) | Fast | Very_Fast (RMM-DIIS)',
'IALGO': 	'algorithm: use only 8 (CG) or 48 (RMM-DIIS)',
'LREAL': 	'non-local projectors in real space',
'ROPT': 	'number of grid points for non-local proj in real space',
'GGA': 	'xc-type: e.g. PE AM or 91',
'VOSKOWN': 	'use Vosko, Wilk, Nusair interpolation',
'DIPOL': 	'center of cell for dipol',
'AMIX':	'tags for mixing',
'BMIX': 	'tags for mixing',
'WEIMIN':	'special control tags',
'EBREAK': 'special control tags',
'DEPER': 	'special control tags',
'TIME': 	'special control tag',
'LWAVE':	'create WAVECAR/CHGCAR/LOCPOT',
'LCHARG':	'create WAVECAR/CHGCAR/LOCPOT',
'LVTOT':	'create WAVECAR/CHGCAR/LOCPOT',
'LVHAR': 	'create WAVECAR/CHGCAR/LOCPOT',
'LELF': 	'create ELFCAR',
'LORBIT': 	'create PROOUT',
'NPAR': 	'parallelization over bands',
'LSCALAPACK' :	'switch off scaLAPACK',
'LSCALU': 	'switch of LU decomposition',
'LASYNC': 'overlap communcation with calculations'}

		self.attr_dic = {}

	def isnum(self, n):
		# ------------------ Define if n is or not a number ------------------ # 
		# n     :   VAR     :   VAR to check if it is a numerical VAR
		# return :  BOOL    : True/False
		try: float(n); return True
		except: return False

	def var_assing(self, var_name, var_value):
		if not var_name in self.attr_dic.keys():
			self.attr_dic[var_name] = var_value
			setattr(self, var_name, var_value)

	def load(self, file_name=None):
		if file_name != None: self.file_name = file_name
		f = open(self.file_name,'r')

		for i, n in enumerate(f):
			try:
				vec = [m for m in n.split(' ') if m != '' and m != '\n']
				if len(vec) > 2 and not '#' in vec[0]:	
					if self.isnum(vec[2]): var_name, var_value = vec[0], float(vec[2])
					else: var_name, var_value = vec[0], vec[2].split('\n')[0] 
					self.var_assing(var_name, var_value)
				if len(vec) > 1 and not '#' in vec[0] and '=' in vec[0]: 
					if self.isnum(vec[1]): var_name, var_value = vec[0], float(vec[1])
					else: var_name, var_value = vec[0].split('=')[0], vec[1].split('\n')[0] 
					self.var_assing(var_name, var_value)
			except: print('can not read line {} :: '.format(i+1), n )

	def view(self, ):
		for n in self.attr_dic.keys():
			if n in self.help.keys():	print( '{} : {} : {}'.format(n, self.attr_dic[n], self.help[n]  ) )
			else:	print( '{} : {} : unknow '.format(n, self.attr_dic[n]) )
			
	def summary(self, ):
		self.view()

# How to...		
#IN = INTCAR('INCAR')
#IN.load('INCAR')
#print(IN.view())






