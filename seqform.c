/* seqform.c
 * 27 Apr 2000
 */

#include <stdio.h>
#include <stdlib.h>
	/* free()       */
#include "alloc.h"
	/* CALLOC(n,s), TALLOC(n,type), T2ALLOC(ptr,nrows,ncols,type)   */
	/* FREE2(ptr,nrows)                                             */
#include "col.h"
#include "rat.h"
	/* typedef Rat                  */
	/* ratadd, ratfromi, ratmult    */
#include "lemke.h"
#include "treedef.h"
#include "sfnf.h"

#include "seqform.h"

/* global variables for sequence form   */

Payvec **sfpay;
int **sfconstr[PLAYERS];

void allocsf(void)
{
    static int oldnseqs1 = 0;
    static int oldconstrows[PLAYERS] = {0, 0, 0};
    int pl, i, j;
    int nrows;
    
    /* payoff matrices, two players only here, init to pay 0        */
    FREE2(sfpay, oldnseqs1);
    oldnseqs1 = nseqs[1];
    T2ALLOC (sfpay, nseqs[1], nseqs[2], Payvec);
    for (i=0; i<nseqs[1]; i++)
	for (j=0; j<nseqs[2]; j++)
	    for (pl=1; pl < PLAYERS; pl++)
		sfpay[i][j][pl-1] = ratfromi(0);
    /* constraint matrices, any number of players           */
    /* sfconstr[0] stays unused                             */
    for (pl=1; pl < PLAYERS; pl++)
	{
	FREE2(sfconstr[pl], oldconstrows[pl]);
	oldconstrows[pl] = nrows = nisets[pl]+1;   /* extra row for seq 0  */
	T2ALLOC (sfconstr[pl], nrows, nseqs[pl], int);
	}
}       /* end of allocsf()     */

void gensf(void)
{
    int pl, i, j;
    Outcome z;
    allocsf();
    
    behavtorealprob(0);     /* get realization probabilities of leaves      */
    
    /* sf payoff matrices                   */
    for (z=outcomes; z < lastoutcome; z++)
	{
	Node u = z->whichnode;
	i = u->defseq[1] - firstmove[1];
	j = u->defseq[2] - firstmove[2];
	for (pl=1; pl < PLAYERS; pl++)
	    sfpay[i][j][pl-1] = ratadd(sfpay[i][j][pl-1],
		    ratmult(u->defseq[0]->realprob, z->pay[pl-1]) );
	}
    /* sf constraint matrices, sparse fill  */
    for (pl=1; pl < PLAYERS; pl++)
	{
	sfconstr[pl][0][0] = 1;     /* empty sequence                       */
	for (i=0; i < nisets[pl]; i++)
	    sfconstr[pl][i+1][(firstiset[pl]+i)->seqin - firstmove[pl]] = -1;
	for (j=1; j < nseqs[pl]; j++)
	    sfconstr[pl][(firstmove[pl]+j)->atiset - firstiset[pl]+1][j] = 1;
	}
}       /* end of  gensf()              */



void sflcp(void)
{
    int i;

    gensf();
    setlcp( nseqs[1] + nisets[2]+1 + nseqs[2] + nisets[1]+1 );
    /* fill  M  */
    /* -A       */
    payratmatcpy(sfpay, 0, 1, 0, nseqs[1], nseqs[2], 
	lcpM, 0, nseqs[1] + nisets[2]+1);
    /* -E\T     */
    intratmatcpy(sfconstr[1], 1, 1, nisets[1]+1, nseqs[1], 
	lcpM, 0, nseqs[1] + nisets[2]+1 + nseqs[2]);
    /* F        */
    intratmatcpy(sfconstr[2], 0, 0, nisets[2]+1, nseqs[2], 
	lcpM, nseqs[1], nseqs[1] + nisets[2]+1 );

    /* -B\T     */
    payratmatcpy(sfpay, 1, 1, 1, nseqs[1], nseqs[2], 
	lcpM, nseqs[1] + nisets[2]+1, 0);
    /* -F\T     */
    intratmatcpy(sfconstr[2], 1, 1, nisets[2]+1, nseqs[2], 
	lcpM, nseqs[1] + nisets[2]+1, nseqs[1] );
    /* E        */
    intratmatcpy(sfconstr[1], 0, 0, nisets[1]+1, nseqs[1], 
	lcpM, nseqs[1] + nisets[2]+1 + nseqs[2], 0 );
    /* define RHS q,  using special shape of SF constraints RHS e,f     */
    for (i = 0; i < lcpdim; i++)
	rhsq[i] = ratfromi(0);
    rhsq[ nseqs[1] ] = ratneg(ratfromi(1));
    rhsq[ nseqs[1]  + nisets[2]+1 + nseqs[2] ] = ratneg(ratfromi(1));
} 

void realplanfromprob(int pl, Rat *rplan)
{
    int i;

    for (i = 0; i < nseqs[pl]; i++)
         rplan[i] = (firstmove[pl] + i)->realprob;
}

Bool iseqrealplantoprob(int pl, Rat *rplan, Bool bcomplain)
{
    int i;
    char s[MAXSTRL];
    Bool isok=1;

    for (i = 0; i < nseqs[pl]; i++)
        if (! ratiseq(rplan[i], (firstmove[pl] + i)->realprob))
	    {
	    isok = 0;
            if (bcomplain)
                {
                seqtoa(firstmove[pl] + i, pl, s);
                printf ("Player %d, move %s has realprob ", pl, s);
                rattoa( (firstmove[pl] + i)->realprob, s);
                printf ("%s", s);
                rattoa( rplan[i], s);
                printf (", but should be %s\n", s);
                }
            else
                break;
            }
    return isok ;
}

int  propermixisets(int pl, Rat *rplan)
{
    int mix = 0;
    int i;
    Move c;
    Iset h;

    for (h = firstiset[pl]; h < firstiset[pl+1]; h++)
        for (c = h->move0, i=0; i < h->nmoves; c++, i++)
            if ( rplan[ c - firstmove[pl] ].num != 0  &&
                !ratiseq( rplan[ c - firstmove[pl] ],
                          rplan[ h->seqin - firstmove[pl]]) )
                {
                mix++ ;
                break ;
                }
    return mix;
}

void outrealplan(int pl, Rat *rplan)
{
    int i;
    char s[MAXSTRL];

    colset(nseqs[pl]) ;
    for (i = 0; i < nseqs[pl]; i++)
        {
        seqtoa(firstmove[pl] + i, pl, s);
        colpr(s);
        }
    for (i = 0; i < nseqs[pl]; i++)
        {
        rattoa(rplan[i], s);
        colpr(s);
        }
   colout();
}

void outbehavstrat(int pl, Rat *rplan, Bool bnewline)
{
    char s[MAXSTRL];
    int i;
    Move c;
    Iset h;
    Rat rprob, bprob;

    for (h = firstiset[pl]; h < firstiset[pl+1]; h++)
        for (c = h->move0, i=0; i < h->nmoves; c++, i++)
            {
            rprob = rplan[ c - firstmove[pl] ];
            if ( rprob.num != 0)
                {
                movetoa(c, pl, s);
                printf(" %s", s);
                bprob = ratdiv( rprob, rplan[ h->seqin - firstmove[pl]]);
                if (!ratiseq(bprob, ratfromi(1) ) )
                    {
                    rattoa(bprob, s);
                    printf(":%s", s);
                    }
                }
            }
    if (bnewline)
        printf("\n");
}

void sfprint(void)
{
    char s[MAXSTRL];
    int i, j, k;
    printf("SF payoffs and constraints:\n");
    
    colset(nseqs[2]+2+nisets[1]);
    colleft(0);
    colpr("pay1");
    for (j=0; j<nseqs[2]; j++)
	{
        seqtoa(firstmove[2] + j, 2, s);    colpr(s);
	}
    colnl();
    colpr("pay2");
    for (j=0; j<nseqs[2]; j++)
	colpr(" ");
    colpr("cons1");
    for (k=1; k<=nisets[1]; k++)
	/* headers for constraint matrix pl1, printed on right of payoffs   */
	colpr((firstiset[1]+k-1)->name);
    for (i=0; i<nseqs[1]; i++)
	{
	/* payoffs player 1 */
	seqtoa(firstmove[1] + i, 1, s);    colpr(s);
	for (j=0; j<nseqs[2]; j++)
	    {
	    rattoa(sfpay[i][j][0], s);     colpr(s);
	    }
	colnl();
	/* payoffs player 2 */
	colpr("");
	for (j=0; j<nseqs[2]; j++)
	    {
	    rattoa(sfpay[i][j][1], s);     colpr(s);
	    }
	/* constraints player 1 */
	for (k=0; k<=nisets[1]; k++)
	    colipr(sfconstr[1][k][i]);    
	colnl();
	}
    /* constraints player 2 */
    for (k=0; k<=nisets[2]; k++)
	{
	colnl();
	if (k==0)
	    colpr("cons2");
	else
	    colpr((firstiset[2]+k-1)->name);
	for (j=0; j<nseqs[2]; j++)
	    colipr(sfconstr[2][k][j]);
	colnl();
	}
    colout();
}       /* end of  sfprint()            */

