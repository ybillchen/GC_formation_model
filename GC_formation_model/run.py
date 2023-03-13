# Licensed under BSD-3-Clause License - see LICENSE

import numpy as np

from . import logo
from . import astro_utils
from .form import form
from .offset import offset
from .assign import assign
from .get_tid import get_tid
from .evolve import evolve

__all__ = ['run']

def run(params):

    if params['verbose']:
        logo.print_logo()
        logo.print_version()
        print('\nRuning model on %d halo(s).'%len(params['subs']))

    allcat_name = params['allcat_base'] + '_p2%g_p3%g.txt'%(
        params['p2'], params['p3'])

    run_params = params
    run_params['allcat_name'] = allcat_name

    run_params['cosmo'] = astro_utils.cosmo(h=run_params['h100'], 
        omega_baryon=run_params['Ob'], omega_matter=run_params['Om'])

    form(run_params)
    offset(run_params)
    assign(run_params)
    get_tid(run_params)
    evolve(run_params)

    if params['verbose']:
        print('\nModel was run on %d halo(s).\n'%len(params['subs']))