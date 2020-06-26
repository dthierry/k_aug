#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

from pyomo.environ import *
from pyomo.opt import SolverFactory, ProblemFormat
from shutil import copyfile
"""Example taken from the sipopt manual, also runs with k_aug
   make sure k_aug is in the path
   make sure dot_sens is in the path
   make sure ipopt is in the path
   make sure ipopt_sens is in the path
please check
This illustrates how to use k_aug with the dsdp_mode for sensitivity
"""

__author__ = 'David Thierry'  #: @2018, rev1 @2020


def main():
    #: Declare Model
    m = ConcreteModel()

    m.i = Set(initialize=[1, 2, 3])

    init_vals = {1:25E+07, 2:0.0, 3:0.0}
    #: Variables
    m.x = Var(m.i, initialize=init_vals)
    #: Objective
    m.oF = Objective(expr=m.x[1]**2 +m.x[2]**2 + m.x[3]**2, 
            sense=minimize)

    m.p1 = Var(initialize=5.0)
    m.p2 = Var(initialize=1.0)
    m.p1_0 = Param(initialize=5.0)
    m.p2_0 = Param(initialize=1.0)

    #: Constraints
    m.c1 = Constraint(expr=6.0 * m.x[1] + 3.0 * m.x[2] + 2.0 * m.x[3] - m.p1 == 0.0)
    m.c2 = Constraint(expr=m.p2 * m.x[1] + m.x[2] - m.x[3] - 1.0 == 0.0)
    m.c1p = Constraint(expr=m.p1 - m.p1_0 == 0.0)
    m.c2p = Constraint(expr=m.p2 - m.p2_0 == 0.0)

    #: ipopt suffixes
    m.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
    m.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
    m.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
    m.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
    m.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)

    #: sipopt suffix
    #m.red_hessian = Suffix(direction=Suffix.EXPORT)
    m.sens_state_0 = Suffix(direction=Suffix.EXPORT)
    m.sens_state_1 = Suffix(direction=Suffix.EXPORT)
    m.sens_state_value_1 = Suffix(direction=Suffix.EXPORT)
    m.sens_init_constr = Suffix(direction=Suffix.EXPORT)
    m.sens_sol_state_1 = Suffix(direction=Suffix.IMPORT)

    m.p1.set_suffix_value(m.sens_state_0, 1)
    m.p1.set_suffix_value(m.sens_state_1, 1)
    m.p1.set_suffix_value(m.sens_state_value_1, 4.5)

    m.p2.set_suffix_value(m.sens_state_0, 2)
    m.p2.set_suffix_value(m.sens_state_1, 2)
    m.p2.set_suffix_value(m.sens_state_value_1, 1.0)

    m.c1p.set_suffix_value(m.sens_init_constr, 1)
    m.c2p.set_suffix_value(m.sens_init_constr, 1)

    ipopt = SolverFactory('ipopt')
    sipopt = SolverFactory('ipopt_sens')

    #: write some options for ipopt sens
    with open('ipopt.opt', 'w') as f:
        f.write('output_file my_ouput.txt\n')
        f.write('run_sens yes\n')
        f.write('n_sens_steps 1\n')
        f.close()

    # see the current suffixes
    for i in m.component_objects(Suffix):
        i.pprint()

    #: First try sipopt to compare
    sipopt.solve(m, tee=True)


    #: A file that will contain some interesting results
    with open("myfile.txt", "w") as f:
        f.write("sipopt [nominal]\n\n")
        for i in m.component_objects(Var):
            i.display(ostream=f)
        f.write("\n\nsipopt [update]\n\n")
        m.sens_sol_state_1.display(ostream=f)

    #: Now, let us use k_aug
    #: k_aug -> compute sensitivities
    kaug = SolverFactory('k_aug')
    #: dotsens -> perform sensitivity step
    dotsens = SolverFactory('dot_sens')

    #: K_AUG SUFFIXES
    
    m.dcdp = Suffix(direction=Suffix.EXPORT)  #: the dummy constraints tag (integer >0) (do not duplicate value)
    m.DeltaP = Suffix(direction=Suffix.EXPORT)  #: the parameter change (float)
    m.var_order = Suffix(direction=Suffix.EXPORT)  #: Important variables (optional) (integer >0) (do not duplicate value)
    #: it is important to use numbers >0 bc ampl uses 0 as NULL
    m.c1p.set_suffix_value(m.dcdp, 1)
    m.c2p.set_suffix_value(m.dcdp, 2)  #: this has to be different, in general an arbitrary value although the dxdp_.dat file will have this order
    #: make sure the order is consistent i.e. 1, 2 and 3. E.g. not 1, 1 and 2 (wrong!) (also this is optional)
    m.x[1].set_suffix_value(m.var_order, 1)
    m.x[2].set_suffix_value(m.var_order, 2)
    # m.x[3].set_suffix_value(m.var_order, 3)  #: we could have all, a subset or none at all
    m.c1p.set_suffix_value(m.DeltaP, 0.5)
    m.c2p.set_suffix_value(m.DeltaP, 0.0)
    #: please check the dsdp_in_.in file generated !
    #: please check the dxdp_.dat file generated if var_order was set! (if you want to see the sens_matrix!!)
    

    #: we need to clear out this file first
    with open('ipopt.opt', 'w') as f:
        f.close()

    ipopt.solve(m, tee=True)


    m.ipopt_zL_in.update(m.ipopt_zL_out)  #: important!
    m.ipopt_zU_in.update(m.ipopt_zU_out)  #: important!    
    #: k_aug
    print('k_aug \n\n\n')
    kaug.options['dsdp_mode'] = ""  #: sensitivity mode!
    kaug.solve(m, tee=True)
    dotsens = SolverFactory("dot_sens")
    dotsens.options["dsdp_mode"] = ""
    dotsens.solve(m, tee=True)

    with open("myfile.txt", "a") as f:
        f.write("\n\n\nk_aug [dsdp_mode!!]\n\n")
        for i in m.component_objects(Var):
            i.display(ostream=f)
    #: note that this mode updates all varables at once

if __name__ == "__main__":
    main()

