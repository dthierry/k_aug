#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

from pyomo.environ import *
from pyomo.opt import SolverFactory, ProblemFormat

"""Example taken from the sipopt manual
please check
https://github.com/coin-or/Ipopt/blob/master/Ipopt/contrib/sIPOPT/examples/redhess_ampl/red_hess.run"""

__author__ = 'David Thierry'  #: @2018

#: Declare Model
m = ConcreteModel()

m.i = Set(initialize=[1, 2, 3])

init_vals = {1:25E+07, 2:0.0, 3:0.0}
#: Variables
m.x = Var(m.i, initialize=init_vals)
#: Objective
m.oF = Objective(expr=(m.x[1] - 1.0)**2 + exp(m.x[2] - 2.0)**2 + (m.x[3] - 3.0)**2, sense=minimize)
m.oF.pprint()
#: Constraints
m.c1 = Constraint(expr=m.x[1] + 2 * m.x[2] + 3 * m.x[3] == 0.0)

#: ipopt suffixes  REQUIRED FOR K_AUG!
m.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
m.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
m.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
m.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
m.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)

#: sipopt suffix
m.red_hessian = Suffix(direction=Suffix.EXPORT)


m.x[2].set_suffix_value(m.red_hessian, 1)
m.x[3].set_suffix_value(m.red_hessian, 2)
ipopt = SolverFactory('ipopt')
sipopt = SolverFactory('ipopt_sens')  #: we can skip sipopt this is just for comparision

kaug = SolverFactory('k_aug')
#: K_AUG SUFFIXES
m.dof_v = Suffix(direction=Suffix.EXPORT)  #: SUFFIX FOR K_AUG
m.rh_name = Suffix(direction=Suffix.IMPORT)  #: SUFFIX FOR K_AUG AS WELL
#: rh_name will tell us which position the corresponding variable has on the reduced hessian text file.
#: be sure to declare the suffix value (order)
m.x[2].set_suffix_value(m.dof_v, 1)
m.x[3].set_suffix_value(m.dof_v, 1)
#: number of dof_v has to be less or equal than m_equality constraints

m.npdp = Suffix(direction=Suffix.EXPORT)
m.c1.set_suffix_value(m.npdp, 1)



kaug.options["compute_inv"] = ""  #: if the reduced hessian is desired.
#: please check the inv_.in file if the compute_inv option was used

#: write some options for ipopt sens
with open('ipopt.opt', 'w') as f:
    f.write('compute_red_hessian yes\n')  #: computes the reduced hessian (sens_ipopt)
    f.write('output_file my_ouput.txt\n')
    f.write('rh_eigendecomp yes\n')
    f.close()
#: Solve
sipopt.solve(m, tee=True)    #: we can skip sipopt this is just for comparision
with open('ipopt.opt', 'w') as f:
    f.close()

ipopt.solve(m, tee=True)  #: Do not skip this solve!

m.ipopt_zL_in.update(m.ipopt_zL_out)
m.ipopt_zU_in.update(m.ipopt_zU_out)

#: k_aug
print('k_aug \n\n\n')
m.write('problem.nl', format=ProblemFormat.nl)
kaug.solve(m, tee=True)  #: always call k_aug AFTER ipopt.
#: please check inv_.in for the inverse of the RH
#: check result_red_hess.txt for the actual RH
