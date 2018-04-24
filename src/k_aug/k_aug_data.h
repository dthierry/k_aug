//
// Created by dav0 on 4/24/18.
//

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


#endif //K_AUG_DATA_H
