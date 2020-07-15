/*
    *Created by dav0 on 4/24/18.
*/

#ifndef INERTIA_STRATEGY_H
#define INERTIA_STRATEGY_H

#include "k_aug_data.h"
#include "asl.h"

fint
inertia_strategy(const fint *row_strt, double *a, fint nvar, fint ncon, fint n_eig, inertia_perts *i_pert,
                 inertia_params i_parm,
                 inertia_options *i_opts, fint *try_n, double log10mu, fint *pert_pivot);


#endif /*INERTIA_STRATEGY_H*/
