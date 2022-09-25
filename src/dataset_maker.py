import datetime
import numpy as np

try:    from src import Logs
except: 
	try: import Logs as Logs
	except: print('WARNING :: Set.import_libraries() :: can not import ORR ')

class DATASET(object):
	def __init__(self,  name:str=None, path:str=None, file:str=None, 
						data_name_list:list=None, system_name_list:list=None, system_path_list:list=None, 
						set_dict:dict=None):
		self.name = name
		self.path = path
		self.file = file

		self.data_name_list   = data_name_list      # [ str, ... ]
		self.system_name_list = system_name_list      # [ [str], [str], ... ]
		self.system_path_list = system_path_list    # [ [str], [str], ... ]

		self.set_dict = set_dict

	@Logs.LogDecorator()
	def save_dataset(self, set_dict:dict=None, file:str=None ) -> bool:
		set_dict = set_dict if not set_dict is None else self.set_dict

		with open(f'{file}', 'w') as ds_file:
			ds_file.write(f'#  === >> {datetime.datetime.now()} Data set generated with dataset_maker.py V0.2 === \n')
			ds_file.write(f'#  === >> {datetime.datetime.now()} {file} === \n')
			ds_file.write('#  \n')
			ds_file.write('Set = {  \n')

			for d_key, d_item in set_dict.items():
				ds_file.write('     %-20s : {  \n' % f'"{d_key}"' )
				for s_key, s_item in d_item.items():
					ds_file.write('          %-20s : %s, \n' % (f'"{s_key}"', f'"{s_item}"') )
				ds_file.write('     },  \n'  )
			ds_file.write('}  \n'  )

		return True

	@Logs.LogDecorator()
	def get_set_dict(self, data_name_list:list=None, system_name_list:list=None, system_path_list:list=None) -> dict:
		set_dict = {}
		for d in data_name_list:
			if not d in set_dict: set_dict[d] = {}

			for sn, sp in zip(system_name_list, system_path_list):
				if not sn in set_dict[d]:
					set_dict[d][sn] = sp

		return set_dict

	@Logs.LogDecorator()
	def get_set_list(self, set_dict:dict=None, save:bool=True) -> list:
		set_dict = set_dict if not set_dict is None else self.set_dict

		data_name_list   = [  d_key  for d_key, d_item in set_dict.items() ]
		system_name_list = [[ s_key  for s_key, s_item in   d_item.items() ] for d_key, d_item in set_dict.items() ]
		system_path_list = [[ s_item for s_key, s_item in   d_item.items() ] for d_key, d_item in set_dict.items() ]
		
		if save:
			self.data_name_list   = data_name_list
			self.system_name_list = system_name_list
			self.system_path_list = system_path_list

		return data_name_list, system_name_list, system_path_list

if __name__ == "__main__":
	# ==== Bimetalic site M ==== #
	metals           = ['Co', 'Cr', 'Cu', 'Fe', 'Mn', 'Ni', 'Sc', 'Ti', 'V', 'Zn',]
	H2H2O = {   'H2'    : '/home/busnengo/Juan.Manuel/alone/H2',
				'H2O'   : '/home/busnengo/Juan.Manuel/alone/H2O', }

	data_dict = { f'{M1}TPyP{M2}' : {**{ss1 : f'/home/busnengo/PAULA/FePC/Metal/M1TPyPM2/M2{M2}/{M1}{ss2}' for ss1, ss2 in zip(['*', '*O', '*OH', '*OOH', 'PDOS'], ['', '/O', '/OH', '/OOH', '/PDOS']) }, **H2H2O} for M1 in metals for M2 in metals }
	dataset = DATASET( set_dict=data_dict )
	dataset.get_set_list()

	dataset.save_dataset(file='/home/akaris/Documents/code/VASP/v4.6/files/dataset/M1TPyPM2/dataset_BiMetals_M1TPyPM2Au_catalysis_site1.py')

	# ==== Bimetalic site C ==== #
	metals           = ['Co', 'Cr', 'Cu', 'Fe', 'Mn', 'Ni', 'Sc', 'Ti', 'V', 'Zn',]
	H2H2O = {   'H2'    : '/home/busnengo/Juan.Manuel/alone/H2',
				'H2O'   : '/home/busnengo/Juan.Manuel/alone/H2O', }

	data_dict = { f'{M1}TPyP{M2}' : {**{ss1 : f'/home/busnengo/PAULA/FePC/Metal/M1TPyPM2/M2{M2}/{M1}{ss2}' for ss1, ss2 in zip(['*', '*O', '*OH', '*OOH', 'PDOS'], ['', '/Os', '/OHs', '/OOHs', '/PDOS']) }, **H2H2O} for M1 in metals for M2 in metals }
	dataset = DATASET( set_dict=data_dict )
	dataset.get_set_list()

	dataset.save_dataset(file='/home/akaris/Documents/code/VASP/v4.6/files/dataset/M1TPyPM2/dataset_BiMetals_M1TPyPM2Au_catalysis_site2.py')

	# ==== Bimetalic site M FREESTANDING ==== #
	metals           = ['Co', 'Cr', 'Cu', 'Fe', 'Mn', 'Ni', 'Sc', 'Ti', 'V', 'Zn',]
	H2H2O = {   'H2'    : '/home/busnengo/Juan.Manuel/alone/H2',
				'H2O'   : '/home/busnengo/Juan.Manuel/alone/H2O', }

	data_dict = { f'{M1}TPyP{M2}' : {**{ss1 : f'/home/busnengo/PAULA/FePC/Metal/M1TPyPM2/M2{M2}/{M1}{ss2}' for ss1, ss2 in zip(['*', '*O', '*OH', '*OOH', 'PDOS'], ['/FREE', '/O/FREE', '/OH/FREE', '/OOH/FREE', ]) }, **H2H2O} for M1 in metals for M2 in metals }
	dataset = DATASET( set_dict=data_dict )
	dataset.get_set_list()

	dataset.save_dataset(file='/home/akaris/Documents/code/VASP/v4.6/files/dataset/M1TPyPM2/dataset_BiMetals_M1TPyPM2FREE_catalysis_site1.py')

	# ==== Bimetalic site C FREESTANDING ==== #
	metals           = ['Co', 'Cr', 'Cu', 'Fe', 'Mn', 'Ni', 'Sc', 'Ti', 'V', 'Zn',]
	H2H2O = {   'H2'    : '/home/busnengo/Juan.Manuel/alone/H2',
				'H2O'   : '/home/busnengo/Juan.Manuel/alone/H2O', }

	data_dict = { f'{M1}TPyP{M2}' : {**{ss1 : f'/home/busnengo/PAULA/FePC/Metal/M1TPyPM2/M2{M2}/{M1}{ss2}' for ss1, ss2 in zip(['*', '*O', '*OH', '*OOH', 'PDOS'], ['/FREE', '/Os/FREE', '/OHs/FREE', '/OOHs/FREE', ]) }, **H2H2O} for M1 in metals for M2 in metals }
	dataset = DATASET( set_dict=data_dict )
	dataset.get_set_list()

	dataset.save_dataset(file='/home/akaris/Documents/code/VASP/v4.6/files/dataset/M1TPyPM2/dataset_BiMetals_M1TPyPM2FREE_catalysis_site2.py')