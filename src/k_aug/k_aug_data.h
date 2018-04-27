/*
// Created by dav0 on 4/24/18.
*/

#ifndef K_AUG_DATA_H
#define K_AUG_DATA_H

typedef struct inertia_strategy_params{
    double km;
    double kp;
    double kbp;
    double kc;
    double dcb;
    double d_w0;
    double dmin;
    double dmax;
} inertia_params;

typedef struct inertia_strategy_options{
    int no_inertia;
    int always_perturb_jacobian;
    int jacobian_perturbed;
} inertia_options;

typedef struct inertia_perturbations{
    double d_w;
    double d_c;
    double d_w_last;
    double d_c_last;
    int jacobian_perturbed;
} inertia_perts;

typedef struct linear_solver_options{
    int want_accurate;
    int max_refinement_steps;
    int max_inertia_steps;
    int ls_internal_scaling;
    double residual_ratio_max;
    int max_memory_al_steps;
    double pivot_tol0;
    double pivtol_max;
} linsol_opts;


#endif /*K_AUG_DATA_H*/
