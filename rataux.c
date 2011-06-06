/* rataux.c
 * 12 July 2000
 * auxiliary routines for rational arithmetic
 */

#include <stdio.h>
        /* printf, sprintf      */
#include <math.h>
        /* floor()              */
#include <limits.h>
        /* INT_MAX,  INT_MIN    */
                                  
#include "rat.h"
#include "rataux.h"

Rat contfract(double x, int accuracy)
{
    int n0, n1, d0, d1 ;
    double xfl, nnext, dnext;
    Rat result;

    xfl = floor(x);
    n0 = 1;
    d0 = 0;
    n1 = (int) xfl ;
    d1 = 1;

    while(1)
        {
        if (x < xfl + 0.5 / INT_MAX)    /* next inverse too large */
            break;
        x = 1 / (x - xfl);
        xfl = floor(x);
        dnext = d1 * xfl + d0 ;
        nnext = n1 * xfl + n0 ;
        if ( dnext > accuracy || nnext > INT_MAX || nnext < INT_MIN )
            break ;
        d0 = d1 ;
        d1 = (int) dnext ;
        n0 = n1 ;
        n1 = (int) nnext ;
        }
    result.num = n1 ;
    result.den = d1 ;
    return result;
}

void ratscalarprod(int dim, Rat *vec1, Rat *vec2, Rat *result)
{
    int j;
    double sum = (double) result->num / (double) result->den ;

    for (j=0; j < dim; j++)
        sum += ((double) vec1[j].num * (double) vec2[j].num) /
	       ((double) vec1[j].den * (double) vec2[j].den ) ;
    *result = contfract (sum, INT_MAX);
}
