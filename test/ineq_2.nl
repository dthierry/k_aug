g3 1 1 0	# problem unknown
 2 3 1 2 0 	# vars, constraints, objectives, ranges, eqns
 0 0 0 0 0 0	# nonlinear constrs, objs; ccons: lin, nonlin, nd, nzlb
 0 0	# network constraints: nonlinear, linear
 0 0 0 	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0 	# discrete variables: binary, integer, nonlinear (b,c,o)
 4 2 	# nonzeros in Jacobian, obj. gradient
 12 8	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
C0	#Time
n0
C1	#Limit[bands]
n0
C2	#Limit[coils]
n0
O0 1	#Total_Profit
n0
x2	# initial guess
0 6000.000059999374
1 1400.0000140003547
r	#3 ranges (rhs's)
1 40
0 0.0 6000.0
0 0.0 4000.0
b	#2 bounds (on variables)
3
3
k1	#intermediate Jacobian column lengths
2
J0 2
0 0.005
1 0.007142857142857143
J1 1
0 1
J2 1
1 1
G0 2
0 25
1 30
