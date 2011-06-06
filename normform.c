/* normform.c
 * 12 Jul 2000
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
	/*  strcpy      */

#include "alloc.h"
	/* CALLOC(n,s), TALLOC(n,type), T2ALLOC(ptr,nrows,ncols,type)   */
	/* FREE2(ptr,nrows)                                             */
#include "col.h"
#include "rat.h"
#include "lemke.h"

#include "mp.h"         /* for  mixedtorealplan ()      */

#include "treedef.h"
#include "sfnf.h"

#include "normform.h"

/* global variables for normal form             */
/* nf payoffs [row][col]                        */
Payvec **nfpay;
/* number of pure strategies for each player
 * set by  gennf() via numstratsnfpre(pl)
 */
int  nstrats[PLAYERS]; 
Move *movetuple[PLAYERS];       /* movetuple encoding pure strategy     */   
Rat  *mixedstrat[PLAYERS];      /* mixed strategy (vector of probs)     */ 

void allocnf(void)
{
    int pl, i, j;
    /* payoff matrices, two players only here, init to pay 0        */
    T2ALLOC (nfpay, nstrats[1], nstrats[2], Payvec);
    for (i=0; i<nstrats[1]; i++)
	for (j=0; j<nstrats[2]; j++)
	    for (pl=1; pl < PLAYERS; pl++)
		nfpay[i][j][pl-1] = ratfromi(0);
    for (pl=0; pl < PLAYERS; pl++)    
	{
	movetuple[pl] = TALLOC (nisets[pl], Move);
	mixedstrat[pl] = TALLOC (nstrats[pl], Rat);
	}
}

int numstratsnfpre(int pl)
{
    Iset h;
    Move c;
    for (h = firstiset[pl+1]-1; h >= firstiset[pl]; h--)
	h->ncontin = 0;
    for (c = firstmove[pl+1]-1; ; c--)
	{
	c->ncompat = 1;
	for (h = firstiset[pl+1]-1; h >= firstiset[pl]; h--)
	    if (h->seqin == c)
		/* h one of possibly many parallel info sets    */
		{
		h->prefact = c->ncompat;
		c->ncompat *= h->ncontin ;
		}
	if (c == firstmove[pl]) 
	    break;
	c->offset = c->atiset->ncontin;
	c->atiset->ncontin += c->ncompat;
	}
    return c->ncompat;
}       /* end of  numstratsnfpre(pl)        */

int seqtostratlist (Move seq, int pl, int *list)
{
    int nl;         /* number of list elements      */
    int i, j, off, right, left;
    /* initialize, usually  nl = 1  if  seq leads to leaf   */
    nl = seq->ncompat;
    for (i=0; i < nl; i++)
	list[i] = i;
    while (seq != firstmove[pl])    /* seq not empty sequence       */
	{
	off = seq->offset;
	for (i=0; i<nl; i++)
	    list[i] += off;
	right = seq->atiset->prefact;
	if (right>1)        
	    /* seq->atiset  has parallel isets to its right */
	    {
	    for (i=nl-1; i >= 0; i--)
		for (j = right-1; j >=0 ; j--)
		    list[i * right + j] = list[i] * right + j;
	    nl *= right;
	    }
	right *= seq->atiset->ncontin;      /* multiply with current iset   */
	left = seq->atiset->seqin->ncompat / right;
	
	/* debug only (integer divisiblity) */
	if (left * right != seq->atiset->seqin->ncompat)
	    printf("not an integer division!!!\n");
    
	/* left>1  means  seq->atiset  has parallel isets to its left       */
	for (j=1; j < left; j++)
	    for (i=0; i < nl; i++)
		list[j * nl + i] = j * right + list[i];
	nl *= left;
	seq = seq->atiset->seqin;
	}
    return nl;
}       /* end of  int seqtostratlist (seq, pl, list)   */

void gennf(void)
{
    int pl, i, j;
    int  *slist[PLAYERS];   /* list of strategies for nf generation         */
    int nl[PLAYERS];        /* no of strategies compatible with leave seq   */
    Outcome z;
    Payvec  v;
    Payvec *nfrow;
    int     nfcolpos;
    /* determine nf size and allocate       */
    for (pl=0; pl < PLAYERS; pl++)
	nstrats[pl] = numstratsnfpre(pl);
    allocnf();
	
    for (pl=0; pl < PLAYERS; pl++)
	slist[pl] = TALLOC (nstrats[pl], int);
    
    behavtorealprob(0);     /* get realization probabilities of leaves      */
    
    /* nf payoff matrices                   */
    for (z=outcomes; z < lastoutcome; z++)
	{
	Node u = z->whichnode;
	for (pl=1; pl < PLAYERS; pl++)
	    {
	    v[pl-1] = ratmult(u->defseq[0]->realprob, z->pay[pl-1] );
	    nl[pl] = seqtostratlist (u->defseq[pl], pl, slist[pl]);
	    }
	for (i=0; i < nl[1]; i++)
	    {
	    nfrow = nfpay[slist[1][i]];
	    for (j=0; j < nl[2]; j++)
		{
		nfcolpos = slist[2][j];
		for (pl=1; pl < PLAYERS; pl++)
		    nfrow[nfcolpos][pl-1] =
		    ratadd(nfrow[nfcolpos][pl-1], v[pl-1]);
		}
	    }
	}   /* end of  for all outcomes  z  */
    for (pl=0; pl < PLAYERS; pl++)
	free(slist[pl]);
}       /* end of  gennf()              */

void nflcp(void)
{
    int i;

    gennf();
    setlcp( nstrats[1] + 1 + nstrats[2] + 1 );
    /* fill  M  */
    /* -A       */
    payratmatcpy(nfpay, 0, 1, 0, nstrats[1], nstrats[2],
	lcpM, 0, nstrats[1] + 1);
    /* -E\T     */
    for (i = 0; i < nstrats[1]; i++)
	lcpM[i][nstrats[1] + 1 + nstrats[2]] = ratneg(ratfromi(1));
    /* F        */
    for (i = 0; i < nstrats[2]; i++)
	lcpM[nstrats[1]][nstrats[1] + 1 + i] = ratfromi(1);
    /* -B\T     */
    payratmatcpy(nfpay, 1, 1, 1, nstrats[1], nstrats[2],
	lcpM, nstrats[1] + 1, 0);
    /* -F\T     */
    for (i = 0; i < nstrats[2]; i++)
	lcpM[ nstrats[1] + 1 + i] [nstrats[1]] = ratneg(ratfromi(1));
    /* E        */
    for (i = 0; i < nstrats[1]; i++)
	lcpM[ nstrats[1] + 1 + nstrats[2]] [i] = ratfromi(1);

    /* define RHS q     */
    for (i = 0; i < lcpdim; i++)
	rhsq[i] = ratfromi(0);
    rhsq[ nstrats[1] ] = ratneg(ratfromi(1));
    rhsq[ nstrats[1]  + 1 + nstrats[2] ] = ratneg(ratfromi(1));
}

void strattomovetuple (int strat, int pl)
{
    int i, sinh, up;
    Move c;
    Iset h;
    int  *saux;     /* aux no strategies at each iset               */   
    saux = TALLOC (nisets[pl], int);
    
    for (i=0; i<nisets[pl]; i++)
	{
	h = firstiset[pl] + i;
	if (h->seqin == firstmove[pl])
	    {
	    sinh = (strat / h->prefact) % h->ncontin ;
	    c = h->move0;   
	    while (sinh < c->offset)        /* offsets decrease with move   */
		c++;
	    movetuple[pl][i] = c;
	    saux[i] = sinh - c->offset;
	    }
	else
	    {
	    up = h->seqin->atiset - firstiset[pl];  /* parent info set no.  */
	    if (movetuple[pl][up] != h->seqin)
		movetuple[pl][i] = NULL;            /* i  is irrelevant     */
	    else
		{
		sinh = (saux[up] / h->prefact) % h->ncontin ;
		c = h->move0;   
		while (sinh < c->offset)    /* offsets decrease with move   */
		    c++;
		movetuple[pl][i] = c;
		saux[i] = sinh - c->offset;
		}
	    }
	}
    free (saux);
}       /* end of  strattomovetuple (strat, pl)                 */

void behavtomixed(int pl)
{
    Rat r;
    int i, strat;
    for (strat = 0; strat < nstrats[pl]; strat++)
	{
	strattomovetuple(strat, pl);
	/* multiply behavior probabilities      */
	r = ratfromi(1);
	for (i=0; i<nisets[pl]; i++)
	    if(movetuple[pl][i]) 
		r = ratmult(r, movetuple[pl][i]->behavprob) ;
		/* check if mixed strat probabilities themselves
		 * troublesome.  Does not seem so.
		 * ...........comment out code..........
		{
		char s[MAXSTRL];
		r = ratmult(r, movetuple[pl][i]->behavprob) ;
                rattoa(r, s);
		printf("Debug: mixed prob %2d is %s\n", i, s) ;
		}
		 */
	mixedstrat[pl][strat] = r;
	}
}       /* end of   behavtomixed(pl)            */

void mixedtorealplan(int pl, Rat *mixed, Rat *rplan)
{
    int r ;
    int i, len;
    int *list;
    mp num, den;        /* numerator and denominator of sum     */
    mp x, y;

    list = TALLOC(nstrats[pl], int);
    
    /* for test purposes, we even compute the probability
     * of the empty sequence. Otherwise, we would start like this:
     *       rplan[0] = ratfromi(1);  
     *       for (r = 1 ....
     */
    for (r = 0; r < nseqs[pl]; r++)
	{
        itomp(0, num);
        itomp(1, den);
        len = seqtostratlist(firstmove[pl] + r, pl, list);
	for (i=0; i<len; i++)
            {
            itomp(mixed[list[i]].den, y);
            mulint (y, num, num);
            itomp(mixed[list[i]].num, x);
            mulint (den, x, x);
            linint (num, 1, x, 1);
            mulint (y, den, den);
            reduce(num, den) ;
            }
        /* mptoi( num, &(rplan[r].num), 1 ); */
        if (mptoi( num, &(rplan[r].num), 1 ))
	    printf(" numerator in mixed strategy probability %3d overflown\n",
	   	   i);
        /* mptoi( den, &(rplan[r].den), 1 ); */
        if (mptoi( den, &(rplan[r].den), 1 ))
	    printf(" denominator in mixed strategy probability %3d overflown\n",
	   	   i);
        }
    free(list);
}

int movetupletoa (int pl, char* s)
{
    int i, pos=0;
    if (nisets[pl] == 0)
        {
        strcpy(s, "()");
        return 2;
        }
    for (i=0; i<nisets[pl]; i++)
	pos += movetoa(movetuple[pl][i], pl, s+pos);
    return pos;
}       /* end of  int movetupletoa (pl, *s)            */

int  supportsize(int pl, Rat *mixed)
{
    int i;
    int suppsize = 0;

    for (i = 0; i < nstrats[pl]; i++)
        if (mixed[i].num != 0)
            suppsize++ ;
    return suppsize; 
}

void outmixed(int pl, Rat *mixed, Bool bnewline)
{
    int i;
    char s[MAXSTRL];
    Rat prob;

    for (i = 0; i < nstrats[pl]; i++)
        {
        prob = mixed[i] ;
        if ( prob.num != 0)
            {
            strattomovetuple (i, pl);
            movetupletoa(pl, s);
            printf(" %s", s);
            if (!ratiseq(prob, ratfromi(1) ) )
                {
                rattoa(prob, s);
                printf(":%s", s);
                }
            }
        }
    if (bnewline)
        printf("\n");
}

void twolineoutmixed(int pl, Rat *mixed)
{
    int i;
    char s[MAXSTRL];
    int *support;
    int suppsize = 0;


    support = TALLOC( nstrats[pl], int);

    for (i = 0; i < nstrats[pl]; i++)
        if (mixed[i].num != 0)
	    support[suppsize++] = i ;
    colset(suppsize) ;
    for (i = 0; i < suppsize; i++)
        {
        strattomovetuple (support[i], pl);
        movetupletoa(pl, s);
        colpr(s);
        }
    for (i = 0; i < suppsize; i++)
        {
        rattoa(mixed[support[i]], s);
        colpr(s);
        }
   colout();
   free(support);
}

void nfprint(void)
{
    char s[MAXSTRL];
    int  i, j;
    printf("NF payoffs:\n");
    
    colset(nstrats[2]+1);
    colleft(0);
    colpr("pay1");
    for (j=0; j<nstrats[2]; j++)
	{
	strattomovetuple(j, 2);
	movetupletoa(2, s);    colpr(s);
	}
    colpr("pay2");
    colnl();
    for (i=0; i<nstrats[1]; i++)
	{
	strattomovetuple(i, 1);
	movetupletoa(1, s);    colpr(s);
	/* payoffs player 1 */
	for (j=0; j<nstrats[2]; j++)
	    {
	    rattoa(nfpay[i][j][0], s);        colpr(s);
	    }
	/* payoffs player 2 */
	colpr("");
	for (j=0; j<nstrats[2]; j++)
	    {
	    rattoa(nfpay[i][j][1], s);        colpr(s);
	    }
	colnl();
	}
    colout();
}       /* end of  nfprint()    */

void compatstrats(int pl)
{
    char s[MAXSTRL];
    Move seq;
    int i,j,nl;
    int **compattable;
    int  *slist;    /* list of strategies   */
    
    T2ALLOC(compattable, nstrats[pl], nseqs[pl], int);
    slist = TALLOC (nstrats[pl], int);
    
    for (seq = firstmove[pl]; seq < firstmove[pl+1]; seq++)
	{
	j = seq - firstmove[pl];
	nl = seqtostratlist(seq, pl, slist);
	for (i=0; i<nl; i++)
	    compattable [slist[i]] [j] = 1;
	}
    
    printf("\nStrategies compatible with sequences for player %d:\n", pl);
    
    colset(nseqs[pl]+1);
    colleft(0);
    colpr("");
    for (seq = firstmove[pl]; seq < firstmove[pl+1]; seq++)
	{
	seqtoa(seq, pl, s);    colpr(s);
	}
    for (i=0; i<nstrats[pl]; i++)
	{
	strattomovetuple(i, pl);
	movetupletoa(pl, s);   colpr(s);
	for (j=0; j<nseqs[pl]; j++)
	    if (compattable[i][j])
		colpr("x");
	    else
		colpr("");
	}
    colout();
    FREE2(compattable, nstrats[pl]);
    free(slist);
}       /* end of  compatstrats(pl)    */

