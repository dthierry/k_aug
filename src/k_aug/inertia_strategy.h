//
// Created by dav0 on 4/24/18.
//

#ifndef INERTIA_STRATEGY_H
#define INERTIA_STRATEGY_H

#include "k_aug_data.h"

int
inertia_strategy(int *row_strt, double *a, int nvar, int ncon, int n_eig, double *d_w, double *d_c, double *d_w_last,
                 double *d_c_last, inertia_params i_parm, inertia_options i_opts, int *try_n, double log10mu, int *pert_pivot, int *jac_pert);


#endif //INERTIA_STRATEGY_H
