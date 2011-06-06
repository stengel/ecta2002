/* rataux.h
 * 12 July 2000
 * auxiliary routines for rational arithmetic
 */

/* #include before: "rat.h" */

/* return the rational number  a/b  that approximates  x  best
 * among all such rationals with  1 <= b <= accuracy
 */ 
Rat contfract(double x, int accuracy);

/* add to  *result  
 * the scalar product of the vectors  vec1  and  vec2  of dimension  dim,  
 * approximated as good as possible
 * note:  *result  must be initialized to zero to get proper scalar product
 */
void ratscalarprod(int dim, Rat *vec1, Rat *vec2, Rat *result);

