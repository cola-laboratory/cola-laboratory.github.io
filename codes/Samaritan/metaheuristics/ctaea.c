/*
 * ctaea.c:
 *  This file implements the main procedures of C-TAEA. It is based on the following reference:
 *
 * K. Li, R. Chen, G. Fu and X. Yao "Two-Archive Evolutionary Algorithm for Constrained Multi-Objective Optimization",
 * IEEE Trans. Evol. Comput. accepted for publication, July 2018.
 *
 * Note that this is the version 1.0, more comprehensive functions will be augmented later.
 *
 * Authors:
 *  Ke Li <k.li@exeter.ac.uk>
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
 *
 * Copyright (c) 2018 Ke Li, Renzhi Chen
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

# include "../header/metaheuristics.h"

int c_pool_count;
int c_max_density;
double rho_CA, rho_DA;
int *c_CA_density, *c_CA_location, *c_DA_density, *c_DA_location;

/* Routine for binary tournament */
individual_real *ctaea_tournament (individual_real *ind1, individual_real *ind2)
{
    int flag;

    if (ind1->cv > -EPS && ind2->cv > -EPS)
    {
        flag = check_dominance (ind1, ind2);
        if (flag == 1)
            return (ind1);
        else if (flag == -1)
            return (ind2);
        else
        {
            if ((randomperc()) <= 0.5)
                return (ind1);
            else
                return (ind2);
        }
    }
    else if (ind1->cv > -EPS && ind2->cv < 0)
        return (ind1);
    else if (ind1->cv < 0 && ind2->cv > -EPS)
        return (ind2);
    else
    {
        if ((randomperc()) <= 0.5)
            return (ind1);
        else
            return (ind2);
    }
}

/* Calculate perpendicular distance (inspired by the PBI) from a point to a given vector */
double perpendicular_distance (double *p, double *vector)
{
    int i;
    double d1, d2, nl;

    d1 = d2 = nl = 0.0;
    for (i = 0; i < number_objective; i++)
    {
        d1 += (p[i] - ideal_point[i]) * vector[i];
        nl += vector[i] * vector[i];
    }
    nl = sqrt (nl);
    d1 = fabs (d1) / nl;

    for (i = 0; i < number_objective; i++)
    {
        d2 += ((p[i] - ideal_point[i]) - d1 * (vector[i] / nl)) *
                ((p[i] - ideal_point[i]) - d1 * (vector[i] / nl));
    }
    d2 = sqrt (d2);

    return d2;
}

/* Association procedure in C-TAEA */
void taea_association (population_real *pop, list *pool, int *density, int *location)
{
    int i, j;
    int min_idx;
    double distance, min_distance;
    double diff, maxFun, feval;

    list *temp;

    c_max_density = 0;
    for (i = 0; i < number_weight; i++)
        density[i] = 0;

    temp = pool->child;
    do {
        i = temp->index;
        min_distance = perpendicular_distance (pop->ind[i].obj, lambda[0]);

        min_idx = 0;
        for (j = 1; j < number_weight; j++)
        {
            distance = perpendicular_distance (pop->ind[i].obj, lambda[j]);
            if (distance < min_distance)
            {
                min_distance = distance;
                min_idx      = j;
            }
        }
        density[min_idx]++;
        location[i] = min_idx;
        if (density[min_idx] > c_max_density)
            c_max_density = density[min_idx];

        maxFun = -1.0;
        for (j = 0; j < number_objective; j++)
        {
            diff = fabs (pop->ind[i].obj[j] - ideal_point[j]);
            if (lambda[min_idx][j] == 0)
                feval = diff / 0.0001;
            else
                feval = diff / lambda[min_idx][j];
            if (feval > maxFun)
                maxFun = feval;
        }
        pop->ind[i].fitness = maxFun;

        temp = temp->child;
    } while (temp != NULL);

    return;
}

/* Check dominance relation according to a new 2-obj problem, i.e. CV and ASF */
int CA_dominance (individual_real *a, individual_real *b)
{
    int flag1, flag2;

    flag1 = flag2 = 0;

    if (-a->cv < -b->cv)
        flag1 = 1;
    else
    {
        if (-a->cv > -b->cv)
            flag2 = 1;

    }
    if (a->fitness < b->fitness)
        flag1 = 1;
    else
    {
        if (a->fitness > b->fitness)
            flag2 = 1;
    }

    if (flag1 == 1 && flag2 == 0)
        return (1);
    else if (flag1 == 0 && flag2 == 1)
        return (-1);
    else
        return (0);
}

/* Fill the CA with good infeasible solutions. Here we formulate a new two-objective problem to do so. */
void fill_CA_nd (list *pool, population_real *mixed_pop, population_real *pop, int idx)
{
    int i;
    int flag;
    int end;
    int front_size, archive_size;

    list *elite;
    list *temp1, *temp2;

    struct double_with_index *array;

    elite = (list *) malloc (sizeof(list));
    elite->index  = -1;
    elite->parent = NULL;
    elite->child  = NULL;

    archive_size = idx;
    do {
        temp1 = pool->child;
        insert (elite, temp1->index);
        front_size = 1;
        temp2 = elite->child;
        temp1 = del (temp1);
        temp1 = temp1->child;
        do {
            temp2 = elite->child;
            if (temp1 == NULL)
                break;
            do {
                end  = 0;
                flag = CA_dominance (&(mixed_pop->ind[temp1->index]), &(mixed_pop->ind[temp2->index]));
                if (flag == 1)
                {
                    insert (pool, temp2->index);
                    temp2 = del (temp2);
                    front_size--;
                    temp2 = temp2->child;
                }
                if (flag == 0)
                    temp2 = temp2->child;
                if (flag == -1)
                    end = 1;
            } while (end != 1 && temp2 != NULL);
            if (flag == 0 || flag == 1)
            {
                insert (elite, temp1->index);
                front_size++;
                temp1 = del (temp1);
            }
            temp1 = temp1->child;
        } while (temp1 != NULL);
        if ((archive_size + front_size) <= popsize)
        {
            temp2 = elite->child;
            do {
                copy_ind (&(mixed_pop->ind[temp2->index]), &(pop->ind[idx]));
                archive_size++;
                temp2 = temp2->child;
                idx++;
                if (idx == popsize)
                    break;
            } while (temp2 != NULL);
        }
        else
        {
            i = 0;
            array = malloc (sizeof(struct double_with_index) * front_size);
            temp2 = elite->child;
            do {
                array[i].idx = temp2->index;
                array[i].x   = mixed_pop->ind[temp2->index].cv;
                temp2        = temp2->child;
                i++;
            } while (temp2 != NULL);
            qsort (array, front_size, sizeof (struct double_with_index), double_with_index_smaller_cmp);

            i = 0;
            while (idx < popsize)
            {
                copy_ind (&(mixed_pop->ind[array[i].idx]), &(pop->ind[idx]));
                archive_size++;
                idx++;
                i++;
            }
            free (array);
        }

        temp2 = elite->child;
        do {
            temp2 = del (temp2);
            temp2 = temp2->child;
        } while (elite->child !=NULL);
    } while (archive_size < popsize);

    // garbage collection
    while (elite != NULL)
    {
        temp1 = elite;
        elite = elite->child;
        free (temp1);
    }

    return;
}

/* 'selected_pool' corresponds {F1+F2+...Fl}, where Fl is the last acceptable front */
void taea_assign_rank (population_real *pop, int size, list *selected_pool)
{
    int i, j;
    int flag, end;
    int rank     = 0;
    int cur_size = 0;
    int request  = popsize;
    c_pool_count = 0;

    list *pool;
    list *elite;
    list *temp1, *temp2;

    pool  = (list *) malloc (sizeof(list));
    elite = (list *) malloc (sizeof(list));
    pool->index   = -1;
    pool->parent  = NULL;
    pool->child   = NULL;
    temp1 = selected_pool->child;
    temp2 = pool;
    while (temp1 != NULL)
    {
        insert (temp2, temp1->index);
        temp1 = temp1->child;
        temp2 = temp2->child;
    }
    temp1 = selected_pool->child;    // empty the 'selected' list
    while (temp1 != NULL)
    {
        temp2 = temp1->child;
        del (temp1);
        temp1 = temp2;
    }
    elite->index  = -1;
    elite->parent = NULL;
    elite->child  = NULL;

    i = 0;
    do
    {
        temp1 = pool->child;
        insert (elite, temp1->index);
        temp1 = del (temp1);
        temp1 = temp1->child;
        do
        {
            temp2 = elite->child;
            if (temp1 == NULL)
                break;

            do
            {
                end  = 0;
                flag = check_dominance (&(pop->ind[temp1->index]), &(pop->ind[temp2->index]));
                if (flag == 1)
                {
                    insert (pool, temp2->index);
                    temp2 = del (temp2);
                    temp2 = temp2->child;
                }
                if (flag == 0)
                {
                    temp2 = temp2->child;
                }
                if (flag == -1)
                {
                    end = 1;
                }
            }
            while (end != 1 && temp2 != NULL);
            if (flag == 0 || flag == 1)
            {
                insert (elite, temp1->index);
                temp1 = del (temp1);
            }
            temp1 = temp1->child;
        }
        while (temp1 != NULL);

        // copy each level into candidate
        if (cur_size < request)
        {
            temp2 = elite->child;
            do
            {
                insert (selected_pool, temp2->index);
                c_pool_count++;
                pop->ind[temp2->index].rank = rank;
                cur_size++;
                temp2 = temp2->child;
            }
            while (temp2 != NULL);
            rank++;
        }
        else if (cur_size < size)
        {
            temp2 = elite->child;
            do
            {
                pop->ind[temp2->index].rank = rank;
                cur_size++;
                temp2 = temp2->child;
            }
            while (temp2 != NULL);
            rank++;
        }

        temp2 = elite->child;
        do {
            temp2 = del (temp2);
            temp2 = temp2->child;
        }
        while (elite->child !=NULL);
    }
    while (cur_size < size);

    // garbage collection
    while (pool != NULL)
    {
        temp1 = pool;
        pool  = pool->child;
        free (temp1);
    }
    while (elite != NULL)
    {
        temp1 = elite;
        elite = elite->child;
        free (temp1);
    }
}

/* Get the proportion of non-dominated solutions in CA and DA */
void calculate_nd_proportion (population_real *CA, population_real *DA)
{
    int i, j;
    int flag_CA, flag_DA;
    int nd_CA, nd_DA;

    population_real *mixed_pop;
    mixed_pop = (population_real *) malloc (sizeof(population_real));
    allocate_memory_pop (mixed_pop, 2 * popsize);

    merge (CA, DA, mixed_pop); // combine CA and DA to form a hybrid population
    nd_CA = nd_DA = 0;

    // count the number of non-dominated points in the CA
    for (i = 0; i < popsize; i++)
    {
        for (j = 0; j < 2 * popsize; j++)
        {
            flag_CA = check_dominance (&(CA->ind[i]), &(mixed_pop->ind[j]));
            if (flag_CA == -1)
                break;
        }
        if (flag_CA != -1)
            nd_CA++;
    }

    // count the number of non-dominated points in the DA
    for (i = 0; i < popsize; i++)
    {
        for (j = 0; j < 2 * popsize; j++)
        {
            flag_DA = check_dominance (&(DA->ind[i]), &(mixed_pop->ind[j]));
            if (flag_DA == -1)
                break;
        }
        if (flag_DA != -1)
            nd_DA++;
    }

    rho_CA = (double) nd_CA / (nd_CA + nd_DA);
    rho_DA = (double) nd_DA / (nd_CA + nd_DA);

    deallocate_memory_pop (mixed_pop, popsize);

    return;
}

/* Identify the best solution in among the solutions stored in the 'pool' */
int find_best (population_real *pop, list *pool)
{
    int min_idx;
    int flag, end;

    double min_fitness;

    list *elite;
    list *temp1, *temp2;

    elite = (list *) malloc (sizeof(list));

    elite->index  = -1;
    elite->parent = NULL;
    elite->child  = NULL;

    // find the current non-dominated solutions
    temp1 = pool->child;
    insert (elite, temp1->index);
    temp1 = del (temp1);
    temp1 = temp1->child;
    do {
        temp2 = elite->child;
        if (temp1 == NULL)
            break;

        do {
            end  = 0;
            flag = check_dominance (&(pop->ind[temp1->index]), &(pop->ind[temp2->index]));
            if (flag == 1)
            {
                insert (pool, temp2->index);
                temp2 = del (temp2);
                temp2 = temp2->child;
            }
            if (flag == 0)
            {
                temp2 = temp2->child;
            }
            if (flag == -1)
            {
                end = 1;
            }
        } while (end != 1 && temp2 != NULL);
        if (flag == 0 || flag == 1)
        {
            insert (elite, temp1->index);
            temp1 = del (temp1);
        }
        temp1 = temp1->child;
    } while (temp1 != NULL);

    // find the best solution from the non-dominated set
    temp2       = elite->child;
    min_fitness = pop->ind[temp2->index].fitness;
    min_idx     = temp2->index;
    while (temp2->child != NULL)
    {
        temp2 = temp2->child;
        if (pop->ind[temp2->index].fitness < min_fitness)
        {
            min_fitness = pop->ind[temp2->index].fitness;
            min_idx     = temp2->index;
        }
    }

    // garbage collection
    while (elite != NULL)
    {
        temp1 = elite;
        elite = elite->child;
        free (temp1);
    }

    return min_idx;
}

/* Update mechanism of the CA */
void CA_selection (population_real *mixed_pop, population_real *CA)
{
    int i, j, k;
    int cur_idx, delete_idx, temp_idx;
    int num_feasible;

    double temp_distance, close_distance, worst_fitness;

    list *loop1, *loop2;
    list *temp, *temp_feasible, *temp_infeasible, *delete_ind;
    list *feasible_pool, *infeasible_pool;

    feasible_pool = (list *) malloc (sizeof(list));
    feasible_pool->index  = -1;
    feasible_pool->parent = NULL;
    feasible_pool->child  = NULL;

    infeasible_pool = (list *) malloc (sizeof(list));
    infeasible_pool->index  = -1;
    infeasible_pool->parent = NULL;
    infeasible_pool->child  = NULL;

    // find the feasible solutions in the 'mixed_pop'
    temp_feasible   = feasible_pool;
    temp_infeasible = infeasible_pool;
    num_feasible    = 0;
    for (i = 0; i < 2 * popsize; i++)
    {
        if (mixed_pop->ind[i].cv > -EPS)
        {
            num_feasible++;
            insert (temp_feasible, i);
            temp_feasible = temp_feasible->child;
        }
        else
        {
            insert (temp_infeasible, i);
            temp_infeasible = temp_infeasible->child;
        }
    }

    if (num_feasible > popsize)
    {
        taea_assign_rank (mixed_pop, num_feasible, feasible_pool);    // 'feasible_pool' corresponds {F1+F2+...Fl}, where Fl is the last acceptable front
        taea_association (mixed_pop, feasible_pool, c_CA_density, c_CA_location);    // associate solutions in 'feasible_pool' to the corresponding sub-regions
        while (c_pool_count > popsize)
        {
            worst_fitness = -1.0;
            for (i = 0; i < number_weight; i++)
            {
                if (c_CA_density[i] == c_max_density)
                {
                    temp = feasible_pool->child;
                    do {
                        cur_idx = temp->index;
                        if (c_CA_location[cur_idx] == i && mixed_pop->ind[cur_idx].fitness > worst_fitness)
                        {
                            delete_idx = cur_idx;
                            delete_ind = temp;
                            worst_fitness = mixed_pop->ind[cur_idx].fitness;

                            // find the nearest neighbhor distance to 'cur_idx'
                            loop1 = feasible_pool->child;
                            close_distance = INFINITY;
                            do {
                                temp_idx = loop1->index;
                                if (temp_idx != cur_idx)
                                {
                                    temp_distance = euclidian_distance (mixed_pop->ind[cur_idx].obj, mixed_pop->ind[temp_idx].obj, number_objective);
                                    if (temp_distance < close_distance)
                                        close_distance = temp_distance;
                                }
                                loop1 = loop1->child;
                            } while (loop1 != NULL);

                            // find another pair of solutions which have smaller distance than 'close_distance'
                            loop1 = feasible_pool->child;
                            do {
                                j = loop1->index;
                                loop2 = feasible_pool->child;
                                do {
                                    k = loop2->index;
                                    if (j != k)
                                    {
                                        temp_distance = euclidian_distance (mixed_pop->ind[j].obj, mixed_pop->ind[k].obj, number_objective);
                                        if (temp_distance < close_distance && c_CA_location[j] == i && c_CA_location[k] == i)
                                        {
                                            close_distance = temp_distance;
                                            if (mixed_pop->ind[j].fitness > mixed_pop->ind[k].fitness)
                                            {
                                                delete_idx = j;
                                                delete_ind = loop1;
                                            }
                                            else if (mixed_pop->ind[j].fitness < mixed_pop->ind[k].fitness)
                                            {
                                                delete_idx = k;
                                                delete_ind = loop2;
                                            }
                                            else
                                            {
                                                if ((randomperc()) <= 0.5)
                                                {
                                                    delete_idx = j;
                                                    delete_ind = loop1;
                                                }
                                                else
                                                {
                                                    delete_idx = k;
                                                    delete_ind = loop2;
                                                }
                                            }
                                        }
                                        else if (temp_distance == close_distance && c_CA_location[j] == i && c_CA_location[k] == i)
                                        {
                                            if (mixed_pop->ind[j].fitness > mixed_pop->ind[k].fitness)
                                            {
                                                if (mixed_pop->ind[j].fitness > worst_fitness)
                                                {
                                                    delete_idx = j;
                                                    delete_ind = loop1;
                                                }
                                            }
                                            else if (mixed_pop->ind[j].fitness < mixed_pop->ind[k].fitness)
                                            {
                                                if (mixed_pop->ind[k].fitness > worst_fitness)
                                                {
                                                    delete_idx = k;
                                                    delete_ind = loop2;
                                                }
                                            }
                                            else
                                            {
                                                if (mixed_pop->ind[j].fitness > worst_fitness)
                                                {
                                                    if ((randomperc()) <= 0.5)
                                                    {
                                                        delete_idx = j;
                                                        delete_ind = loop1;
                                                    }
                                                    else
                                                    {
                                                        delete_idx = k;
                                                        delete_ind = loop2;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    loop2 = loop2->child;
                                } while (loop2 != NULL);
                                loop1 = loop1->child;
                            } while (loop1 != NULL);
                            worst_fitness = mixed_pop->ind[delete_idx].fitness;
                        }
                        temp = temp->child;
                    } while (temp != NULL);
                }
            }
            c_pool_count--;
            del (delete_ind);
            c_CA_density[c_CA_location[delete_idx]]--;
            c_CA_location[delete_idx] = -1;

            // update 'c_max_density'
            c_max_density = c_CA_density[0];
            for (j = 1; j < number_weight; j++)
            {
                if (c_max_density < c_CA_density[j])
                    c_max_density = c_CA_density[j];
            }
        }

        // copy the remaining solutions in 'feasible_pool' to the new CA
        j    = 0;
        temp = feasible_pool->child;
        do {
            i = temp->index;
            copy_ind (&(mixed_pop->ind[i]), &(CA->ind[j]));
            j++;
            temp = temp->child;
        } while (temp != NULL);
    }
    else    // num_feasible <= popsize
    {
        // copy all feasible solutions to the new CA
        i    = 0;
        temp = feasible_pool->child;
        while (temp != NULL)
        {
            copy_ind (&(mixed_pop->ind[temp->index]), &(CA->ind[i]));
            i++;
            temp = temp->child;
        }

        if (i < popsize)
        {
            // fill the gap with promising infeasible solutions
            taea_association (mixed_pop, infeasible_pool, c_CA_density, c_CA_location);
            fill_CA_nd (infeasible_pool, mixed_pop, CA, i);
        }
    }

    /* Garbage collection */
    while (feasible_pool != NULL)
    {
        temp          = feasible_pool;
        feasible_pool = feasible_pool->child;
        free (temp);
    }
    while (infeasible_pool != NULL)
    {
        temp            = infeasible_pool;
        infeasible_pool = infeasible_pool->child;
        free (temp);
    }

    return;
}

/* Update mechanism of the DA */
void DA_selection (population_real *mixed_pop, population_real *CA, population_real *DA)
{
    int i, j, k;
    int itr;
    int best_idx;
    int updated_size;

    list *pool, *temp_archive, *temp;

    pool         = (list *) malloc (sizeof(list));
    pool->index  = -1;
    pool->parent = NULL;
    pool->child  = NULL;

    temp_archive         = (list *) malloc (sizeof(list));
    temp_archive->index  = -1;
    temp_archive->parent = NULL;
    temp_archive->child  = NULL;

    // association procedure for the CA
    temp = pool;
    for (i = 0; i < popsize; i++)
    {
        insert (temp, i);
        temp = temp->child;
    }
    taea_association (CA, pool, c_CA_density, c_CA_location);

    // clear the 'pool' list
    while (pool->child != NULL)
    {
        temp = pool;
        pool = pool->child;
        free (temp);
    }

    // association procedure for the mixed population, i.e. DA + offspring
    temp = pool;
    for (i = 0; i < 2 * popsize; i++)
    {
        insert (temp, i);
        temp = temp->child;
    }
    taea_association (mixed_pop, pool, c_DA_density, c_DA_location);

    itr = 0;
    updated_size = 0;
    while (updated_size < popsize)
    {
        itr++;
        for (i = 0; i < number_weight; i++)
        {
            if (c_CA_density[i] < itr)
            {
                for (j = 0; j < (itr - c_CA_density[i]); j++)
                {
                    if (c_DA_density[i] == 0)
                        break;
                    else    // c_DA_density[i] >= 1
                    {
                        // find solutions of 'mixed_pop' (DA + offspring) located in the i-th sub-region
                        temp = pool->child;
                        do {
                            k = temp->index;
                            if (c_DA_location[k] == i)
                                insert (temp_archive, k);
                            temp = temp->child;
                        } while (temp != NULL);

                        // find the best solution in the i-th sub-region
                        best_idx = find_best (mixed_pop, temp_archive);

                        // clear the 'temp_archive' list
                        while (temp_archive->child != NULL)
                        {
                            temp_archive = temp_archive->child;
                            free (temp_archive->parent);
                        }
                        temp_archive->index  = -1;
                        temp_archive->parent = NULL;
                        temp_archive->child  = NULL;

                        // remove the selected 'best_idx' solution
                        temp = pool->child;
                        do {
                            k = temp->index;
                            if (best_idx == k)
                            {
                                del (temp);
                                break;
                            }
                            temp = temp->child;
                        } while (temp != NULL);
                        c_DA_density[c_DA_location[best_idx]]--;
                        c_DA_location[best_idx] = -1;

                        copy_ind (&(mixed_pop->ind[best_idx]), &(DA->ind[updated_size]));
                        updated_size++;

                        if (updated_size == popsize)
                        {
                            while (pool != NULL)
                            {
                                temp = pool;
                                pool = pool->child;
                                free (temp);
                            }
                            while (temp_archive != NULL)
                            {
                                temp         = temp_archive;
                                temp_archive = temp_archive->child;
                                free (temp);
                            }

                            return;
                        }
                    }
                }
            }
        }
    }
}

/* Crossover operation in C-TAEA */
void crossover_taea (population_real *CA, population_real *DA, population_real *offspring_pop)
{
    int i, temp, rand;
    int *a1, *a2;
    individual_real *parent1, *parent2;

    a1 = (int *) malloc (popsize * sizeof(int));
    a2 = (int *) malloc (popsize * sizeof(int));
    for (i = 0; i < popsize; i++)
        a1[i] = a2[i] = i;

    for (i = 0; i < popsize; i++)
    {
        rand     = rnd (i, popsize - 1);
        temp     = a1[rand];
        a1[rand] = a1[i];
        a1[i]    = temp;
        rand     = rnd (i, popsize - 1);
        temp     = a2[rand];
        a2[rand] = a2[i];
        a2[i]    = temp;
    }

    for (i = 0; i < popsize; i += 2)
    {
        if (i == popsize - 1)
            i--;
        if (rho_CA > rho_DA)
            parent1 = ctaea_tournament (&CA->ind[a1[i]], &CA->ind[a1[i + 1]]);
        else
            parent1 = ctaea_tournament (&DA->ind[a1[i]], &DA->ind[a1[i + 1]]);
        if (rndreal (0, 1) < rho_CA)
            parent2 = ctaea_tournament (&CA->ind[a2[i]], &CA->ind[a2[i + 1]]);
        else
            parent2 = ctaea_tournament (&DA->ind[a2[i]], &DA->ind[a2[i + 1]]);

        sbx_crossover (parent1, parent2, &offspring_pop->ind[i], &offspring_pop->ind[i + 1]);
    }

    free (a1);
    free (a2);

    return;
}

/* main procedure of C-TAEA */
void CTAEA (population_real *pop, population_real *offspring_pop, population_real *mixed_pop)
{
    int i;
    int generation;

    population_real *CA, *DA;

    CA = pop;
    DA = (population_real *) malloc (sizeof(population_real));
    allocate_memory_pop (DA, popsize);

    generation       = 1;
    evaluation_count = 0;
    printf ("|\tThe %d run\t|\t1%%\t|", run_index);

    // initialization process
    initialize_uniform_weight (NULL);
    c_CA_density  = malloc (sizeof(int) * number_weight);
    c_CA_location = malloc (sizeof(int) * 2 * popsize);
    c_DA_density  = malloc (sizeof(int) * number_weight);
    c_DA_location = malloc (sizeof(int) * 2 * popsize);

    print_error (number_weight != popsize, 1, "Number of weight vectors must be equal to the population size!");

    initialize_population_real (CA);
    evaluate_population (CA);
    initialize_idealpoint (CA);
    for (i = 0; i < popsize; i++)
        copy_ind (&(CA->ind[i]), &(DA->ind[i]));
    track_evolution (pop, generation, 0);

    while (evaluation_count < max_evaluation)
    {
        print_progress ();
        calculate_nd_proportion (CA, DA);
        crossover_taea (CA, DA, offspring_pop);
        mutation_real (offspring_pop);
        evaluate_population (offspring_pop);

        // update the ideal point
        for (i = 0; i < popsize; i++)
            update_ideal_point (&(offspring_pop->ind[i]));

        // update the CA
        merge (CA, offspring_pop, mixed_pop);
        CA_selection (mixed_pop, CA);

        // update the DA
        merge (DA, offspring_pop, mixed_pop);
        DA_selection (mixed_pop, CA, DA);

        generation++;
        track_evolution (CA, generation, evaluation_count >= max_evaluation);
    }
    merge (CA, DA, mixed_pop);
    CA_selection (mixed_pop, CA);

    track_evolution (CA, generation, evaluation_count >= max_evaluation);

    // garbage collection
    deallocate_memory_pop (DA, popsize);
    free (c_CA_density);
    free (c_DA_density);
    free (c_CA_location);
    free (c_DA_location);

    return;
}
