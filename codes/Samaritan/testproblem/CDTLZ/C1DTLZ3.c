/*
 * C1-DTLZ3.c
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Copyright (c) 2017 Renzhi Chen, Ke Li
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

void c1dtlz3 (individual_real *ind)
{
    int i, j, k, aux;
    double gx, fsum, re, r; // r determines the outermost feasible boundary's radius
    double *xreal, *obj;

    if (number_objective == 3)
        r = 9.0;
    else if (number_objective == 5)
        r = 12.5;
    else if (number_objective == 8)
        r = 12.5;
    else if (number_objective == 10)
        r = 15.0;
    else if (number_objective == 15)
        r = 15.0;
    else
        r = 9.0;

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

    // evaluate the constraint violation
    fsum = 0;
    for (i = 0; i < number_objective; i++)
        fsum = obj[i] * obj[i] + fsum;
    re = (fsum - 16.0) * (fsum - r * r);
    if (re >= 0.0)
        re = 0.0;

    ind->cv = re;

    return;
}
