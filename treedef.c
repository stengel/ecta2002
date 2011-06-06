/* treedef.c
 * 19 Apr 2000
 * defining a game tree
 */
#include <stdio.h>
#include <stdlib.h>
	/* RAND_MAX     */
#include "alloc.h"
	/* CALLOC(n,s), TALLOC(n,type), T2ALLOC(ptr,nrows,ncols,type)   */
	/* FREE2(ptr,nrows)                                             */
#include "col.h"
#include "rat.h"
	/* typedef Rat                  */
	/* ratadd, ratfromi, ratmult    */

#include "treedef.h"

/* ------------- global variables ------------------------------------- */
char an1[PLAYERS] = { '!', 'A', 'a' };
char an2[PLAYERS] = { '/', 'Z', 'z' };

/* game tree                                    */
Node    nodes;          /* nodes of game tree   */
Node    root;           /* &nodes[ROOT]         */
Iset    isets;          /* information sets     */
Move    moves;          /* moves & sequences    */
Outcome outcomes;       /* outcomes             */

Node lastnode;
Iset firstiset[PLAYERS+1];
Move firstmove[PLAYERS+1];
Outcome lastoutcome;

int nseqs[PLAYERS];    
int nisets[PLAYERS];

static Payvec maxpay;


int nextrandpay (void)
{
    return rand() / (RAND_MAX / MAXRANDPAY) ;
}

void alloctree(int nn, int ni, int nm, int no)
{
    free(nodes);
    free(isets);
    free(moves);
    free(outcomes);

    nodes = TALLOC(nn, struct node); 
    lastnode = nodes + nn;
    
    isets = TALLOC(ni, struct iset);
    firstiset[PLAYERS] = isets + ni;
    
    moves = TALLOC(nm, struct move);
    firstmove[PLAYERS] = moves + nm;
    
    outcomes = TALLOC(no, struct outcome);
    lastoutcome = outcomes + no;
}       /* end of alloctree(nn, ni, nm, no)             */

Bool genseqin(void)
{
    Bool isnotok = 0;
    int pl;
    Iset  h, lasth;
    Node u;
    Move seq;

    /* set  nseqs[], nisets[]               */
    for (pl=0; pl < PLAYERS; pl++)
	{
	nseqs[pl]   = firstmove[pl+1] - firstmove[pl];
	nisets[pl] = firstiset[pl+1] - firstiset[pl];
	}
    /* set seqin for all isets to NULL      */
    lasth = firstiset[PLAYERS];
    for (h = isets; h < lasth; h++)
	h->seqin = NULL;
    
    for (pl=0; pl < PLAYERS; pl++)
	root->defseq[pl] = firstmove[pl];
    root->iset->seqin = firstmove[root->iset->player];
		     
    for (u = root+1; u < lastnode; u++)
	{
	if (u->father >= u)
	    /* tree is not topologically sorted     */
	    {
	    isnotok = 1;
	    printf("tree not topologically sorted: father %d ", 
		    u->father - nodes);
	    printf("is larger than node %d itself.\n", u - nodes);
	    }
    
	/* update sequence triple, new only for move leading to  u  */
	for (pl=0; pl < PLAYERS; pl++)
	    u->defseq[pl] = u->father->defseq[pl];
	u->defseq[u->reachedby->atiset->player] = u->reachedby;
	
	/* update sequence for iset, check perfect recall           */
	if (! (u->terminal))
	    {
	    h = u->iset;
	    seq = u->defseq[h->player];
	    if (h->seqin == NULL)
		h->seqin = seq;
	    else if (seq != h->seqin)
		/* not the same as last sequence leading to info set    */
		{
		isnotok = 1;
		/* need output routines for isets, moves, later         */
		printf("imperfect recall in info set no. %d ", h-isets);
		printf("named %s\n", h->name);
		printf("different sequences no. %d,", seq-moves);
		printf(" %d\n", h->seqin-moves);
		}
	    }       /* end of "u decision node"     */
	}           /* end of "for all nodes u"     */
    return isnotok;
}       /* end of  Bool genseqin()               */

void maxpayminusone(Bool bprint)
{
    char s[MAXSTRL];
    int pm;
    Outcome z;
    Payvec addtopay;

    for (pm=0; pm < PLAYERS-1; pm++)
	{
	maxpay[pm] = ratfromi(MINUSINFTY) ;
	for (z = outcomes; z < lastoutcome; z++)
	    if (ratgreat(z->pay[pm], maxpay[pm]))
		maxpay[pm] = z->pay[pm];
	if (bprint)     /* comment to stdout    */
	    {
	    rattoa(maxpay[pm], s);
	    printf("Player %d's maximum payoff is %s , ", pm+1, s);
	    printf("normalize to -1.\n");
	    }
	addtopay[pm] = ratneg( ratadd(maxpay[pm], ratfromi(1)));
	for (z = outcomes; z < lastoutcome; z++)
	    z->pay[pm] = ratadd( z->pay[pm], addtopay[pm]) ;
	}
}

void autoname(void)
{
    int pl, anbase, max, digits, i, i1, j;
    Iset h;

    for (pl=0; pl<PLAYERS; pl++)    /* name isets of player pl      */
	{
	max = anbase = an2[pl]-an1[pl]+1;
	for (digits = 1; max < nisets[pl]; max *= anbase, digits++)
	    ;
	if (digits >= NAMECHARS)
	    {
	    printf("Too many isets (%d) of player %d.  ", nisets[pl], pl);
	    printf("change NAMECHARS to %d or larger\n", digits+1);
	    exit(1);
	    }
	for (i=0; i < nisets[pl]; i++)
	    {
	    i1 = i;
	    h = (firstiset[pl]+i);
	    h->name[digits]='\0';
	    for (j = digits - 1; j>=0; j--)
		{
		h->name[j] = an1[pl] + (i1 % anbase);
		i1 /= anbase;
		}
	    }
	}
}       /* end of  autoname()   */

int movetoa (Move c, int pl, char *s)
{
    if (c == NULL)
	return sprintf(s, "*");
    if (c == firstmove[pl])
	return sprintf(s, "()");
    return sprintf(s, "%s%d", c->atiset->name, c - c->atiset->move0);
}       /* end of  int movetoa (c, pl, *s)      */

int seqtoa (Move seq, int pl, char *s)
{
    int len;
    if (seq == NULL)
	return sprintf(s, "*");
    if (seq == firstmove[pl])
	return sprintf(s, ".");
    len = seqtoa (seq->atiset->seqin, pl, s);       /* recursive call       */
    return len + movetoa(seq, pl, s+len);
}       /* end of  int seqtoa (seq, pl, *s)     */


void rawtreeprint(void)
{
    char s[MAXSTRL];
    int pl;
    Node u;
    Iset h;
    Move c;
    
    /* printing nodes       */
    colset(6 + PLAYERS-1 + PLAYERS);
    colleft(0);
    colpr("node");
    colpr("leaf");
    colpr("iset");
    colpr("father");
    colpr("reachedby");
    colpr("outcome");
    for (pl=1; pl < PLAYERS; pl++)
	{
	sprintf(s,"pay%d", pl);    colpr(s);
	}
    for (pl=0; pl < PLAYERS; pl++)
	{
	sprintf(s,"isp%d", pl);    colpr(s);
	}
    colnl();
    for (u=root; u < lastnode; u++)
	{
	colipr(u - nodes);
	colipr(u->terminal);
	colipr(u->iset - isets);
	colipr(u->father - nodes);
	colipr(u->reachedby - moves);
	colipr(u->outcome - outcomes);
	for (pl=1; pl < PLAYERS; pl++)
	    if (u->terminal)
		{
		rattoa( u->outcome->pay[pl-1], s);
		colpr(s);
		}
	    else
		colpr("");
	for (pl=0; pl < PLAYERS; pl++)
	    colipr(u->defseq[pl] - moves);
	}
    colout();
    printf("\n");
    /* printing isets       */
    colset(8);
    colleft(0);
    colpr("iset");
    colpr("player");
    colpr("nmoves");
    colpr("move0");
    colpr("name");
    colpr("seqin");
    colpr("ncontin");
    colpr("prefact");
    colnl();
    pl=0;
    for (h=isets; h < firstiset[PLAYERS]; h++)
	{
	while (h==firstiset[pl])
	    {
	    sprintf(s,"pl%d:", pl);         colpr(s);
	    colnl();
	    pl++ ;
	    }
	colipr(h-isets);
	colipr(h->player);
	colipr(h->nmoves);
	colipr(h->move0 - moves);
	colpr(h->name);
	colipr(h->seqin - moves);
	colipr(h->ncontin);
	colipr(h->prefact);
	}
    colout();
    printf("\n");
    /* printing moves       */
    colset(9);
    colleft(0);
    colleft(1);
    colleft(2);
    colpr("move");
    colpr("name");
    colpr("seqname");
    colpr("atiset");
    colpr("behavprob");
    colpr("realprob");
    colpr("redsfcol");
    colpr("ncompat");
    colpr("offset");
    colnl();
    pl=0;
    for (c=moves; c < firstmove[PLAYERS]; c++)
	{
	while (c==firstmove[pl])
	    {
	    sprintf(s,"pl%d:", pl);         colpr(s);
	    colnl();
	    pl++ ;
	    }
	/* pl is now the NEXT possible player       */
	colipr(c - moves);
	movetoa(c, pl-1, s);    colpr(s);
	seqtoa(c, pl-1, s);     colpr(s);
	colipr(c->atiset - isets);
	rattoa(c->behavprob, s);    colpr(s);
	rattoa(c->realprob, s);     colpr(s);
	colipr(c->redsfcol);
	colipr(c->ncompat);
	colipr(c->offset);
	}
    colout();
}       /* end of  rawtreeprint()       */
