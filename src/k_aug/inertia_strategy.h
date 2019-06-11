/*
    *Created by dav0 on 4/24/18.
*/

#ifndef INERTIA_STRATEGY_H
#define INERTIA_STRATEGY_H

#include "k_aug_data.h"

int
inertia_strategy(const int *row_strt, double *a, int nvar, int ncon, int n_eig, inertia_perts *i_pert,
                 inertia_params i_parm,
                 inertia_options *i_opts, int *try_n, double log10mu, int *pert_pivot);


#endif /*INERTIA_STRATEGY_H*/
