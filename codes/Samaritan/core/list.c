/*
 * list.c:
 *  This file contains the functions to perform linked list operations.
 *
 * Authors:
 *  Renzhi Chen <rxc332@cs.bham.ac.uk>
 *  Ke Li <k.li@exeter.ac.uk>
 *
 * Institution:
 *  Computational Optimization and Data Analytics (CODA) Group @ University of Exeter
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

#include "../header/global.h"

/* Insert an element X into the list at location specified by NODE */
void insert (list *node, int x)
{
    list *temp;

    if (node == NULL)
    {
        printf("\n EE: asked to enter after a NULL pointer, hence exiting \n");
        exit(1);
    }

    temp = (list *) malloc (sizeof(list));
    temp->index  = x;
    temp->child  = node->child;
    temp->parent = node;
    if (node->child != NULL)
    {
        node->child->parent = temp;
    }
    node->child = temp;

    return;
}

/* Delete the node NODE from the list */
list* del (list *node)
{
    list *temp;

    if (node == NULL)
    {
        printf("\n Error!!! asked to delete a NULL pointer, hence exiting \n");
        exit(1);
    }

    temp = node->parent;
    temp->child = node->child;
    if (temp->child != NULL)
    {
        temp->child->parent = temp;
    }
    free (node);

    return (temp);
}

