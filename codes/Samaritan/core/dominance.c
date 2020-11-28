/*
 * dominance.c:
 *  This file contains the functions to perform the dominance check, constraint violation is considered.
 *
 * Authors:
 *  Ke Li <k.li@exeter.ac.uk>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2017 Ke Li, Renzhi Chen
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "../header/global.h"

/* Check dominace relation bewteen two solutions without considering the constraint violation. */
int check_dominance (individual_real *a, individual_real *b)
{
    int i;
    int flag1, flag2;

    flag1 = flag2 = 0;
    for (i = 0; i < number_objective; i++)
    {
        if (a->obj[i] < b->obj[i])
            flag1 = 1;
        else
        {
            if (a->obj[i] > b->obj[i])
                flag2 = 1;
        }
    }
    if (flag1 == 1 && flag2 == 0)
        return (1);
    else
    {
        if (flag1 == 0 && flag2 == 1)
            return (-1);
        else
            return (0);
    }
}

/* Check dominance relation between two solutions with the consideration of the constraint violation. */
int check_dominance_cv (individual_real *a, individual_real *b)
{
    int i;
    int flag1, flag2;

    flag1 = flag2 = 0;

    if (a->cv < 0 && b->cv < 0)
    {
        if (a->cv > b->cv)
            return (1);
        else if (a->cv < b->cv)
            return (-1);
        else
            return (0);
    }
    else
    {
        if (a->cv < 0 && b->cv > -EPS)
            return (-1);
        else if (a->cv > -EPS && b->cv < 0)
            return (1);
        else
            return check_dominance (a, b);
    }
}
