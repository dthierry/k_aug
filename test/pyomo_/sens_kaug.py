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
https://github.com/coin-or/Ipopt/blob/master/Ipopt/contrib/sIPOPT/examples/redhess_ampl/red_hess.run"""

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
    m.dof_v = Suffix(direction=Suffix.EXPORT)  #: variable tag takes values of 0 or 1
    m.npdp = Suffix(direction=Suffix.EXPORT)  #: constraint tag takes values of float

    #: make sure ||dof_v|| <= nvars - mcons
    #: we can update several variables as long as the previous line is true


    #: Tag the constraint that contains the parameter
    #: The value here indicates the DeltaP of the parameter in each constraint (1 per equation)
    m.c1p.set_suffix_value(m.npdp, 0.5)
    m.c2p.set_suffix_value(m.npdp, 0.0)

    #: we need to clear out this file first
    with open('ipopt.opt', 'w') as f:
        f.close()
    #: k_aug
    print('k_aug \n\n\n')

    # m.write('problem.nl', format=ProblemFormat.nl)
    # f = open('./kaug_sens.res', 'w')
    # f.close()

    #: one by one sensitivity update
    for i in range(1, 4):
        m.dof_v.clear()  #: clear the suffix
        m.x[i].set_suffix_value(m.dof_v, 1)
        ipopt.solve(m)  #: re solve just in case
        m.ipopt_zL_in.update(m.ipopt_zL_out)  #: update bound multipliers
        m.ipopt_zU_in.update(m.ipopt_zU_out)
        #: k_aug call
        kaug.solve(m, tee=True)
        dotsens.solve(m, tee=True)  #: sens step
        # with open('./dot_out.out', 'r') as f:
        #    lines = f.readlines()
        #    f.close()
        #with open('./kaug_sens.res', 'a') as f:
        #    f.write(lines[0])
        #    f.close()
        #copyfile('./dot_out.out', 'res_' + str(i))
        #: write some interesting results
        with open("myfile.txt", "a") as f:
            f.write("\n\nk_aug at position: {}\n\n".format(str(i)))
            for i in m.component_objects(Var):
                i.display(ostream=f)

    #: one by one sensitivity update (version 2: no re-call ipopt)
    ipopt.solve(m)
    for i in range(1, 4):
        m.dof_v.clear()  #: clear the suffix
        m.x[i].set_suffix_value(m.dof_v, 1)
        # ipopt.solve(m)  #: re solve just in case
        m.ipopt_zL_in.update(m.ipopt_zL_out)  #: update bound multipliers
        m.ipopt_zU_in.update(m.ipopt_zU_out)
        #: k_aug call
        kaug.solve(m, tee=True)
        dotsens.solve(m, tee=True)  #: sens step
        # with open('./dot_out.out', 'r') as f:
        #    lines = f.readlines()
        #    f.close()
        #with open('./kaug_sens.res', 'a') as f:
        #    f.write(lines[0])
        #    f.close()
        #copyfile('./dot_out.out', 'res_' + str(i))
        #: write some interesting results
        with open("myfile.txt", "a") as f:
            f.write("\n\nk_aug at position version 2: {}\n\n".format(str(i)))
            for i in m.component_objects(Var):
                i.display(ostream=f)

if __name__ == "__main__":
    main()

