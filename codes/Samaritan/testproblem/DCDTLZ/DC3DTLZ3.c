/*
 * DC3-DTLZ3.c
 *
 * Authors:
 *  Ke Li <k.li@exeter.ac.uk>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *
 * Copyright (c) 2018 Renzhi Chen, Ke Li
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

#include "../../header/problems.h"

void dc3dtlz3 (individual_real *ind)
{
    int i, j, k, aux;
    double a, b;    // parameters of constraints
    double gx, fsum, re;
    double *xreal, *obj;

    obj   = ind->obj;
    xreal = ind->xreal;

    gx = 0.0;
    k  = number_variable - number_objective + 1;
    for (i = number_variable - k; i < number_variable; i++)
        gx += pow ((xreal[i] - 0.5), 2.0) - cos (20.0 * PI * (xreal[i] - 0.5));
    gx = 100.0 * (k + gx);

    for (i = 0; i < number_objective; i++)
        obj[i] = 1.0 + gx;

    for (i = 0; i < number_objective; i++)
    {
        for (j = 0; j < number_objective - (i + 1); j++)
            obj[i] *= cos (PI * 0.5 * xreal[j]);
        if (i != 0)
        {
            aux     = number_objective - (i + 1);
            obj[i] *= sin (PI * 0.5 * xreal[aux]);
        }
    }

    a  = 5;
    b  = 0.5;

    re = cos (a * PI * gx) - b;
    for (i = 0; i < number_variable - k; i++)
    {
        if (cos (a * PI * xreal[i]) - b < re)
            re = cos (a * PI * xreal[i]) - b;
    }

    if (re > 0) re = 0;

    ind->cv = re;
}
