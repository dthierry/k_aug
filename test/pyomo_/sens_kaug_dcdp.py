#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

from pyomo.environ import *
from pyomo.opt import SolverFactory, ProblemFormat
from shutil import copyfile
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

#: ipopt suffixes  REQUIRED FOR K_AUG!
m.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
m.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
m.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
m.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
m.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)

#: sipopt suffix
#:m.red_hessian = Suffix(direction=Suffix.EXPORT)


#.x[2].set_suffix_value(m.red_hessian, 1)
#.x[3].set_suffix_value(m.red_hessian, 2)
ipopt = SolverFactory('ipopt')
sipopt = SolverFactory('ipopt_sens')
kaug = SolverFactory('k_aug', executable="../../bin/k_aug")
#: K_AUG SUFFIXES
#: rh_name will tell us which position the corresponding variable has on the reduced hessian text file.
#: be sure to declare the suffix value (order)

m.dcdp = Suffix(direction=Suffix.EXPORT)
m.DeltaP = Suffix(direction=Suffix.EXPORT)
m.var_order = Suffix(direction=Suffix.EXPORT)

m.c1p.set_suffix_value(m.dcdp, 1)
m.c2p.set_suffix_value(m.dcdp, 2)
m.x[1].set_suffix_value(m.var_order, 1)
m.x[2].set_suffix_value(m.var_order, 2)

#: please check the inv_.in file if the compute_inv option was used

#: write some options for ipopt sens
with open('ipopt.opt', 'w') as f:
    f.close()

#: Solve
# sipopt.solve(m, tee=True)
#with open('ipopt.opt', 'w') as f:
#    f.close()

ipopt.solve(m, tee=True)
m.ipopt_zL_in.update(m.ipopt_zL_out)
m.ipopt_zU_in.update(m.ipopt_zU_out)
#: k_aug
kaug.options['dsdp_mode'] = ""
kaug.solve(m, tee=True)

