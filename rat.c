/* rat.c 
 * computing with rationals
 * 27 Apr 2000
 */

#include <stdio.h>
	/* sprintf	*/
#include "rat.h"

#include "mp.h"
	/* for improved precision in  ratadd(a, b)	*/

Rat ratadd (Rat a, Rat b)
{
    /*
    a.num = a.num * b.den + a.den * b.num;
    a.den *= b.den;
    return ratreduce(a);
    */

    mp num, den, x, y;

    itomp (a.num, num) ;
    itomp (a.den, den) ;
    itomp (b.num, x) ;
    itomp (b.den, y) ;
    mulint (y, num, num);
    mulint (den, x, x);
    linint (num, 1, x, 1);
    mulint (y, den, den);
    reduce(num, den) ;
    mptoi( num, &a.num, 1 );
    mptoi( den, &a.den, 1 );
    return a ; 
}

Rat ratdiv (Rat a, Rat b)
{
    return ratmult(a, ratinv(b) );
}

Rat ratfromi(int i)
{
    Rat tmp;
    tmp.num = i;
    tmp.den = 1;
    return tmp;
}

int ratgcd(int a, int b)
{
    int c;
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    if (a < b) { c=a; a=b; b=c; }
    while (b != 0)
        {
        c = a % b;
        a = b ;
        b = c ;
        }
    return a;
}

Rat ratinv (Rat a)
{
    int x;

    x = a.num ;
    a.num = a.den ;
    a.den = x ;
    return a;
}
Bool ratiseq (Rat a, Rat b)
{
    return (a.num == b.num && a.den == b.den);
}

Bool ratgreat (Rat a, Rat b)
{
    Rat c = ratadd(a, ratneg(b));
    return (c.num > 0);
}

Rat ratmult (Rat a, Rat b)
{
    int x;

    /* avoid overflow in intermediate product by cross-cancelling first
     */
    x = a.num ;
    a.num = b.num ;
    b.num = x ;
    a = ratreduce(a);
    b = ratreduce(b);
    a.num *= b.num;
    a.den *= b.den;
    return ratreduce(a);        /* a  or  b  might be non-normalized    s*/
}

Rat ratneg (Rat a)
        /* returns -a                                           */
{
    a.num = - a.num;
    return  a;
}

Rat ratreduce (Rat a)
{
    if (a.num == 0)
        a.den = 1;
    else
        {
        int div;
        if (a.den < 0)
            {
            a.den = -a.den;
            a.num = -a.num;
            }
        div = ratgcd(a.den, a.num);
        a.num = a.num/div;
        a.den = a.den/div;
        }
    return a;
}

int rattoa (Rat r, char *s)
{
    int l, a;
    l = sprintf(s, "%d", r.num);
    if (r.den != 1)
        {
        a = sprintf(s+l, "/%d", r.den);
        l += a + 1;
        }
    return l;
}

double rattodouble(Rat a)
{
    return (double) a.num / (double) a.den ;
}
