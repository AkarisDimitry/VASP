import datetime
import numpy as np

def dataset_maker(name, name_list, subname_list, path_list):

    with open(f'{name}.py', 'w') as ds_file:
        ds_file.write(f'#  === >> {datetime.datetime.now()} Data set generated with dataset_maker.py V0.0 === \n')
        ds_file.write(f'#  === >> {datetime.datetime.now()} {name} === \n')
        ds_file.write('#  \n')
        ds_file.write('Set = {  \n')
        for nl, snl, pl in zip(name_list, subname_list, path_list):
            ds_file.write(f'\t"{nl}" :{{ \n')
            for snl_, pl_ in zip(snl, pl):

                ds_file.write(f'\t\t"{snl_}"\t\t : \t\t "{pl_}", \n')
            ds_file.write(f'\t}}, \n')
        ds_file.write(f'}} \n')

def cookbook():
    vdw_list = ['rmBJvdw', 'rmDF', 'rmDF2', 'rmDF2_B86', 'rmOPT86', 'rmOPT_PBE']
    ref1_path = '/home/busnengo/PAULA/FePC/Metal/Fe/O2/sample'
    N = 10

    job_list        = [ f'height_{i:02d}'       for i in range(N) ]
    subjob_list     = [ vdw_list                for i in range(N) ]
    path_list       = [[f'{ref1_path}_{i}/{n}' for n in vdw_list]   for i in range(N) ]

    dataset_maker(  'O2_absortion', 
                    job_list,
                    subjob_list,
                    path_list
                    )

'''
'''

