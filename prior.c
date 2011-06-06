/* prior.c
 * 12 July 2000
 * prior generation
 */

#include <stdio.h>
        /* printf, sprintf      */
#include <stdlib.h>
        /* rand(), RAND_MAX     */
#include <math.h>
        /* floor()              */
#include <limits.h>
        /* INT_MAX,  INT_MIN    */

#include "prior.h"
                                  
#include "rat.h"
#include "rataux.h"	/* contfract()	*/
#include "treedef.h"
#include "sfnf.h"
	/* realplan[PLAYERS], behavtorealplan()	*/
#include "seqform.h"
	/* outbehavstrat	*/

/* generate the prior, stored in  moves[]->behavprob,
 * where each move has equal probability
 */ 
static void gencentroid(void)
{
    Move c;
    int pl;
    for (pl=1; pl < PLAYERS; pl++)
        for (c = firstmove[pl]+1; c < firstmove [pl+1]; c++)
            {
            c->behavprob.num = 1;
            c->behavprob.den = c->atiset->nmoves;
            }
}

void genprior(Flagsprior flags)
{
    int pl;
    Iset h;

    if (0 == flags.seed)
	{
    	gencentroid();
	return ;
	}
    /* generate random priors for all information sets	*/
    srand(FIRSTPRIORSEED + flags.seed); 
    for (pl=1; pl < PLAYERS; pl++)
        for (h = firstiset[pl]; h < firstiset [pl+1]; h++)
	    if ( h->nmoves > 2)
	    	{
		fprintf(stderr, "Sorry, only binary info sets so far.\n") ; 
		exit(1) ;
		}
            else 
                {
	        Rat a;
	        double x;

	        x = rand() / (double) RAND_MAX;
	        a = contfract( x, flags.accuracy) ;
	        /* make sure to get a properly mixed prior,
	         * unless  flags.accuracy == 1,
	  	 * in which case we have a random pure strategy
		 * because this statement flips 0 to 1 and vice versa
		 */
	        if (a.num == 0)
	            {
		    a.num = 1 ;
		    a.den = flags.accuracy;
		    }
                else if (a.den == 1)  	/* "else" for pure strategy	*/
	            {
		    a.num = flags.accuracy - 1 ;
		    a.den = flags.accuracy;
		    }
                h->move0->behavprob = a ;
                ((h->move0)+1)->behavprob = ratadd(ratfromi(1), ratneg(a)) ;
                }
}

void outprior(void)
{
    int pl;

    printf("------Prior behavior strategies player 1, 2:\n");
    for ( pl = 1; pl < PLAYERS; pl++ )
        {
	behavtorealprob(pl) ;
	realplanfromprob(pl, realplan[pl]);
	outbehavstrat(pl, realplan[pl], 1);
	}
}
