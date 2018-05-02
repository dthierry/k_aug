g3 0 1 0	# problem silly
 5 4 1 0 4	# vars, constraints, objectives, ranges, eqns
 1 1	# nonlinear constraints, objectives
 0 0	# network constraints: nonlinear, linear
 2 4 1	# nonlinear vars in constraints, objectives, both
 0 0 0 1	# linear network variables; functions; arith, flags
 0 0 0 0 0	# discrete variables: binary, integer, nonlinear (b,c,o)
 10 3	# nonzeros in Jacobian, gradients
 0 0	# max name lengths: constraints, variables
 0 0 0 0 0	# common exprs: b,c,o,c1,o1
S0 2 sens_state_0
1 2
4 1
S0 2 sens_state_1
1 2
4 1
S4 2 sens_state_value_1
1 1
4 4.5
S1 2 sens_init_constr
2 1
3 1
S1 2 dcdp
2 1
3 1
b
2 0
3
2 0
2 0
3
x5
0 0.632653048979739
1 1
2 0.3877551216324176
3 0.02040817061215653
4 5
r
4 1
4 0
4 5
4 1
d4
0 0.28571437923136433
1 0.16326528580108643
2 0.16326528580108643
3 -0.18075807371461186
C0
o2
v1
v0
C1
n0
C2
n0
C3
n0
O0 0
o54
3
o5
v0
n2
o5
v2
n2
o5
v3
n2
k4
2
4
6
8
J0 4
0 0
1 0
2 1
3 -1
J1 4
0 6
2 3
3 2
4 -1
J2 1
4 1
J3 1
1 1
G0 3
0 0
2 0
3 0
