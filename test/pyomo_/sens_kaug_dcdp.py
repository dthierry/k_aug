#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

from pyomo.environ import *
from pyomo.opt import SolverFactory, ProblemFormat
from shutil import copyfile
"""Example taken from the sipopt manual
please check
https://github.com/coin-or/Ipopt/blob/master/Ipopt/contrib/sIPOPT/examples/

This illustrates how to use k_aug with the dsdp_mode for sensitivity"""

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
#: Dummy variables
m.p1 = Var(initialize=5.0)
m.p2 = Var(initialize=1.0)
#: Parameters variables
m.p1_0 = Param(initialize=5.0)
m.p2_0 = Param(initialize=1.0)

#: Constraints
m.c1 = Constraint(expr=6.0 * m.x[1] + 3.0 * m.x[2] + 2.0 * m.x[3] - m.p1 == 0.0)
m.c2 = Constraint(expr=m.p2 * m.x[1] + m.x[2] - m.x[3] - 1.0 == 0.0)
#: Dummy Constraint REQUIRED!
m.c1p = Constraint(expr=m.p1 - m.p1_0 == 0.0)
m.c2p = Constraint(expr=m.p2 - m.p2_0 == 0.0)

#: Ipopt suffixes  REQUIRED FOR K_AUG!
m.dual = Suffix(direction=Suffix.IMPORT_EXPORT)
m.ipopt_zL_out = Suffix(direction=Suffix.IMPORT)
m.ipopt_zU_out = Suffix(direction=Suffix.IMPORT)
m.ipopt_zL_in = Suffix(direction=Suffix.EXPORT)
m.ipopt_zU_in = Suffix(direction=Suffix.EXPORT)

ipopt = SolverFactory('ipopt')
sipopt = SolverFactory('ipopt_sens')
kaug = SolverFactory('k_aug', executable="../../bin/k_aug")
#: K_AUG SUFFIXES
m.dcdp = Suffix(direction=Suffix.EXPORT)  #: the dummy constraints
# m.DeltaP = Suffix(direction=Suffix.EXPORT)  #:
m.var_order = Suffix(direction=Suffix.EXPORT)  #: Important variables (primal)

m.c1p.set_suffix_value(m.dcdp, 1)
m.c2p.set_suffix_value(m.dcdp, 2)
#: make sure the order is consistent i.e. 1, 2 and 3. E.g. not 1, 1 and 2 (wrong!)
m.x[1].set_suffix_value(m.var_order, 1)
m.x[2].set_suffix_value(m.var_order, 2)
# m.x[3].set_suffix_value(m.var_order, 3)  #: we could have all, a subset or none at all

#: please check the dsdp_in_.in file generated !
#: please check the dxdp_.dat file generated if var_order was set!

#: Clear this file
with open('ipopt.opt', 'w') as f:
    f.close()

ipopt.solve(m, tee=True)
m.ipopt_zL_in.update(m.ipopt_zL_out)  #: important!
m.ipopt_zU_in.update(m.ipopt_zU_out)  #: important!
#: k_aug
kaug.options['dsdp_mode'] = ""  #: sensitivity mode!
kaug.solve(m, tee=True)

