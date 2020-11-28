/*
 * utility.c:
 *  This file contains the functions to facilitate some common usages.
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

# include "../header/utility.h"

/* Calculate the L2-norm of a vector */
double norm_vector (double *a)
{
    int i;
    double sum;

    sum = 0;
    for (i = 0; i < number_objective; i++)
        sum += a[i] * a[i];

    return sqrt (sum);
}

/* Calculate the Euclidean distance between two points */
double euclidian_distance (double *a, double *b, int dimension)
{
    int i;
    double distance;

    distance = 0.0;
    for(i = 0; i < dimension; i++)
        distance += (a[i] - b[i]) * (a[i] - b[i]);

    return sqrt(distance);
}


/* Build up multi-level directories */
void _mkdir (const char *dir)
{
    char tmp[256];
    char *p = NULL;
    size_t len;

    snprintf (tmp, sizeof(tmp), "%s", dir);
    len = strlen (tmp);
    if (tmp[len - 1] == '/')
        tmp[len - 1] = 0;
    for (p = tmp + 1; *p; p++)
    {
        if (*p == '/')
        {
            *p = 0;
            mkdir (tmp, S_IRWXU);
            *p = '/';
        }
    }
    mkdir (tmp, S_IRWXU);
}

/* Calculate the combinatorial number for n choose k */
int combination (int n, int k)
{
    int i;

    if (n < k)
        return -1;
    double ans = 1;
    for (i = k + 1; i <= n; i++)
    {
        ans = ans * i;
        ans = ans / (double) (i - k);
    }

    return (int) ans;
}

/* Shuffle the 'perm' array */
void random_permutation (int *perm, int size)
{
    int i, num, start;
    int *index, *flag;

    index = malloc (size * sizeof(int));
    flag  = malloc (size * sizeof(int));
    for (i = 0; i < size; i++)
    {
        index[i] = i;
        flag[i]  = 1;
    }

    num = 0;
    while (num < size)
    {
        start = rnd (0, size - 1);
        while (1)
        {
            if (flag[start])
            {
                perm[num] = index[start];
                flag[start] = 0;
                num++;
                break;
            }
            if (start == (size - 1))
                start = 0;
            else
                start++;
        }
    }

    free (index);
    free (flag);
    return;
}

/* Update the current ideal point */
void update_ideal_point (individual_real *individual)
{
    int i;

    for (i = 0; i < number_objective; i++)
        if (individual->obj[i] < ideal_point[i])
            ideal_point[i] = individual->obj[i];

    return;
}


/* Update the current nadir point */
void update_nadir_point (individual_real *individual)
{
    int i;

    for (i = 0; i < number_objective; i++)
        if (individual->obj[i] > nadir_point[i])
            nadir_point[i] = individual->obj[i];

    return;
}

/* Test for duplicate solutions: 1->duplicate; -1->not duplicate */
int is_duplicate (individual_real *ind1, individual_real *ind2)
{
    int i;
    int counter;

    counter = 0;
    for (i = 0; i < number_objective; i++)
        if (ind1->obj[i] == ind2->obj[i])
            counter++;
    if (counter == number_objective)
        return 1;
    else
        return -1;
}
