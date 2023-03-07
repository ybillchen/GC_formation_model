import numpy as np

from params import params
import form
import offset
# import assign
# import get_tid
# import disrupt

def run(params):

    allcat_name = params['allcat_base'] + '_%.1f_%.1f.txt'%(
        params['p2'], params['p3'])

    run_params = params
    run_params['allcat_name'] = allcat_name

    form.form(run_params)
    offset.offset(run_params)
    # assign.assign(run_params)
    # get_tid.get_tid(run_params)
    # disrupt.disrupt(run_params)

if __name__ == '__main__':
    run(params)