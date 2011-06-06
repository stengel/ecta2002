/* rat.h
 * typedef Bool
 * computing with rationals
 * 22 Apr 2000
 */

typedef int Bool;        /* Boolean value 0/1                   */

typedef struct 
    {
    int        num;    /* numerator    */
    int        den;    /* denominator  */
    }
    Rat;

/* returns sum  a+b, normalized                         */
Rat ratadd (Rat a, Rat b);

/* returns quotient  a/b, normalized                    */
Rat ratdiv (Rat a, Rat b);

/* converts integer i to rational                       */
Rat ratfromi(int i);

/* computes gcd of integers  a  and  b,  0 if both 0    */
int ratgcd(int a, int b);

/* returns Boolean condition that a > b                 */
Bool ratgreat (Rat a, Rat b);

/* returns quotient  1/a, normalized only if  a  is     */
Rat ratinv (Rat a);

/* returns Boolean condition that a==b
 * a, b are assumed to be normalized
 */
Bool ratiseq (Rat a, Rat b);

/* returns product  a*b, normalized                     */
Rat ratmult (Rat a, Rat b);

/* returns -a, normalized only if a normalized          */
Rat ratneg (Rat a);

/* normalizes (make den>0, =1 if num==0)
 * and reduces by  gcd(num,den)
 */
Rat ratreduce (Rat a);

/* converts rational  r  to string  s, omit den 1
 * s  must be sufficiently long to contain result
 * returns length of string
 */
int rattoa (Rat r, char *s);

/* converts rational  a  to  double                     */
double rattodouble (Rat a);
