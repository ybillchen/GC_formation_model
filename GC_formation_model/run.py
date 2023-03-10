import numpy as np

from .form import form
from .offset import offset
from .assign import assign
# from .get_tid import get_tid
# from .disrupt import disrupt

def run(params):

    print('########## GC formation and evolution model ##########')

    allcat_name = params['allcat_base'] + '_%.1f_%.1f.txt'%(
        params['p2'], params['p3'])

    run_params = params
    run_params['allcat_name'] = allcat_name

    run_params['cosmo'] = astro_utils.cosmo(h=run_params['h100'], 
        omega_baryon=run_params['Ob'], omega_matter=run_params['Om'])

    # form(run_params)
    # offset(run_params)
    # assign(run_params)
    # get_tid(run_params)
    # disrupt(run_params)

if __name__ == '__main__':
    run(params)