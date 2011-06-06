/* rsf.c
 * reduced sequence form routines
 * 11 Apr 2000
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

/* global variables for reduced sequence form   */

int irreddim[PLAYERS];
int redsfdim[PLAYERS];
int ** redsfconstr[PLAYERS];
int *  redsfrhs[PLAYERS];
int ** realplfromredsf[PLAYERS];
int *  realplconst[PLAYERS];


void genredsf(int pl)
{
    int i, j;       /* row, column counters */
    int cdim;       /* dimension of c-matrix, = no. redundant vars          */
    int drow;       /* last filled row of d-matrix for irredundant vars     */
    Iset h, lasthinc;
    Move seq;
    /* inverse of matrix of dependent variables, local      */
    int ** depvarinv;
    /* rowinfo [0..cdim) = 
     *      row with 1 of c-matrix of redundant vars 
     * rowinfo [drow..nisets) =  
     *      col with -1 of d-matrix of irred vars
     * rowinfo (nisets<i<nseqs) = col of  depvarinv[][] 
     *      (i = seq->redsfcol)
     */
    int * rowinfo;
    
    T2ALLOC(depvarinv, nisets[pl]+1, nisets[pl]+1, int);
    cdim = 0;    
    drow = nisets[pl]+1;
    rowinfo = TALLOC(nseqs[pl], int);

    /* gather data for c- and d- matrix, c filled from above,
     * d from below
     */ 
    /* initialize the  redsfcol  field of moves of  pl  to -1       */
    for (seq = firstmove[pl]; seq < firstmove[pl+1]; seq++)
	seq->redsfcol = -1;
    
    /* for all isets  h  of player  pl  */
    for (h = firstiset[pl]; h < firstiset[pl+1]; h++)
	if (h->seqin->redsfcol == -1)
	    /* this insequence not yet looked at            */
	    {
	    lasthinc = h;   /* last h in c-matrix, for last dependent var   */
	    h->seqin->redsfcol = cdim ;
	    if (cdim) /* only row 0 of c-matrix should be empty sequence    */
		      /* this segfaults if the empty sequence is elsewhere  */
		rowinfo [cdim] = h->seqin->atiset->seqin->redsfcol ;
	    cdim++ ;
	    }
	else
	    /* iset  h  has insequence already looked at    */
	    {
	    --drow ;
	    h->move0->redsfcol = drow ;     /* identity matrix next to d    */
	    rowinfo [drow] = h->seqin->redsfcol ;   /*  -1 entry of d       */
	    }
    /* all isets processed, assert drow==cdim */
    if (drow != cdim)
	{
	fprintf(stderr, "Error: for player %d, cdim %d not equal to drow %d\n",
			 pl, cdim, drow);
	exit(1);
	}

    /* generate  c^-1 = depvarinv[][]  from c- and d-matrix         */
    /* c-matrix inverse and next to it: 0 matrix                    */
    for (i=0; i<cdim; i++)
	/* c-matrix rows (redundant vars)           */
	{
	for (j=0; j<i; j++)
	    depvarinv[i][j] = 0 ;
	depvarinv[i][i] = -1 ;      /* c is upper triang with -1 on diag    */
	/* the following uses rowinfo[j] giving the only 1 of c in col j    */
	for (j=i+1; j<cdim; j++)
	    depvarinv[i][j] = depvarinv[i] [ rowinfo[j] ] ;
	for (j=cdim; j < nisets[pl]; j++)   /* except for very last col     */
	    depvarinv[i][j] = 0 ;           /* 0 matrix on the right of  c  */
	}
    for (i = cdim; i < nisets[pl] ; i++)
	/* d-matrix rows (irredundant vars)                                 */
	/* inverse is  -d c^-1,  note  d  has neg unit row vectors          */
	{
	for (j=0; j<nisets[pl]; j++)   /* copy row  rowinfo[i]  of  c^-1    */
	    depvarinv[i][j] = depvarinv[ rowinfo[i] ] [ j ] ;
	depvarinv[i][i] = 1 ;       /* identity on right of c-cols          */
	}
    /* now last row and col of  depvarinv[][]
     * the extra col corresponds to the first move at  lasthinc,
     * assuming  cdim > 0, that is,  nisets[pl] > 0,
     * otherwise that column will be the only (= empty) seq of  pl
     */
    seq = (nisets[pl]) ? lasthinc->move0 : firstmove[pl] ;
    seq->redsfcol = nisets[pl] ;
    /* last row of  depvarinv[][] = row of depvarinv corresp to the equatn  */
    /* "empty seq has prob 1", usually empty seq is col 0 of  c             */
    /* assert the latter                                                    */
    if (firstmove[pl]->redsfcol != 0)
	{
	fprintf(stderr, "Error: for player %d, empty seq in col %d,",
			 pl, firstmove[pl]->redsfcol);
	fprintf(stderr, " and not in col 0 of c\n");
	exit(1);
	}
    for (j=0; j < nisets[pl]; j++)   
	depvarinv [ nisets[pl] ] [j] = -depvarinv [0] [j] ;
    /* last col of  depvarinv[][] = col of depvarinv for  lasthinc->move0   */
    /* the latter is in row  cdim-1  of  c = col of c^-1 (=depvarinv)       */
    for (i=0; i < nisets[pl]; i++)
	depvarinv [i] [ nisets[pl] ] = -depvarinv [i] [cdim-1] ;
    /* update main part via rank-1-product of this last column and row  */
    for (i=0; i < nisets[pl]; i++)   
	for (j=0; j < nisets[pl]; j++)   
	    depvarinv [i][j] += 
		depvarinv [i] [ nisets[pl] ] * depvarinv [ nisets[pl] ] [j] ;
    /* lower right element  */
    depvarinv [ nisets[pl] ] [ nisets[pl] ] = 1 ;
    /* depvarinv[][] computed now */

    /* find  redsfcol  for unused sequences = lcp  variables            */
    drow = nisets[pl] + 1 ;
    for (seq = firstmove[pl]; seq < firstmove[pl+1]; seq++)
	if (seq->redsfcol == -1)
	    {
	    seq->redsfcol = drow ;
	    rowinfo[drow] = seq->atiset->seqin->redsfcol ;
	    drow++ ;
	    }
    /* assert drow==nseqs[pl]                                           */
    if (drow != nseqs[pl])
	{
	fprintf(stderr, "Error: for player %d wrong no. %d of redsfcols,",
			 pl, drow);
	fprintf(stderr, "instead of %d.\n", nseqs[pl] );
	exit(1);
	}
    
    /* free old data            */
    free (redsfrhs[pl]);
    FREE2(redsfconstr[pl], irreddim[pl]);
    
    /* assign new dimensions    */
    irreddim[pl] = nisets[pl] + 1 - cdim ;
    redsfdim[pl] = nseqs[pl] - nisets[pl] - 1 ;
    /* allocate and fill reduced seq form constraint right hand side    */
    redsfrhs[pl] = TALLOC(irreddim[pl], int);      
    for (i=0; i < irreddim[pl]; i++)   
	redsfrhs[pl][i] = depvarinv [cdim+i] [nisets[pl]];
    /* allocate and fill reduced seq form constraints                   */
    T2ALLOC(redsfconstr[pl], irreddim[pl], redsfdim[pl], int);      
    for (i=0; i < irreddim[pl]; i++)   
	for (j=0; j < redsfdim[pl]; j++)   
	    redsfconstr[pl][i][j] = 
		depvarinv [cdim+i] [ rowinfo [nisets[pl]+1+j] ];
		
    /* allocate and fill  realplfromredsf/const                         */
    realplconst[pl] = TALLOC(nseqs[pl], int);      
    T2ALLOC(realplfromredsf[pl], nseqs[pl], redsfdim[pl], int);      
    

    for (seq = firstmove[pl]; seq < firstmove[pl+1]; seq++)
    for (j=0; j < redsfdim[pl]; j++)   
	if ( (i = seq->redsfcol ) > nisets[pl] ) 
	    /* seq is a free variable, identity matrix  */
	    {
	    realplconst[pl][seq - firstmove[pl]] = 0;
	    realplfromredsf[pl][seq - firstmove[pl]] [j] 
		= (i == nisets[pl]+1+j ) ? 1 : 0 ;
	    }
	else
	    /* seq is a dependent variable, use neg column of inverse   */
	    {
	    realplconst[pl][seq - firstmove[pl]] 
		= depvarinv [i] [ nisets[pl] ];
	    realplfromredsf[pl][seq - firstmove[pl]] [j] 
		= - depvarinv [i] [ rowinfo [ nisets[pl]+1+j ] ];
	    }
    /* done with constructing matrices, free local matrices             */
    free( rowinfo );
    FREE2(depvarinv, nisets[pl]+1); 
}       /* end of  genredsf(pl)                                         */

void setrsfcomplconstr (int paytopl, Rat ** intoM,
	int rowoffset, int coloffset, Bool negpaymatrix,
	Rat * rhsvec, Bool negrhs)
{
    int othpl = 3 - paytopl;   /* other player         */
    int i,j,k, a;
    Rat s;
    Rat * tmprow;   /* temporary row for triple matrix product comptn   */
    tmprow = TALLOC( nseqs[othpl] , Rat);
    
    /* first blockrow  of dim   redsfdim[paytopl]       */
    for (i = 0; i < redsfdim[paytopl]; i++)
	{
	/* compute  tmprow (with index  j ) as  
	 * if  paytopl==1:  row  i  of  K\T A, 
	 * if  paytopl==2:  row  i  of  L\T B\T
	 */
	for (j = 0; j < nseqs[othpl]; j++)
	    {
	    s = ratfromi(0);
	    for (k = 0; k < nseqs[paytopl]; k++)
		/* a = entry of  K\T  resp. L\T         */
		if ((a = realplfromredsf[paytopl][k][i] ) != 0)
		    s = ratadd( s, ratmult( ratfromi(a),
		      (paytopl==1) ? sfpay[k][j][0] : sfpay[j][k][1] )) ;
	    tmprow[j] = s;
	    }
	/* compute i,k entry of 
	 * if  paytopl==1:  K\T A L   = Ahat  
	 * if  paytopl==2:  L\T B\T K = Bhat \T
	 * comments from now only for  paytopl==1
	 */
	for (k = 0; k < redsfdim[othpl]; k++)
	    {
	    s = ratfromi(0);
	    for (j = 0; j < nseqs[othpl]; j++)
		if ( (a = realplfromredsf [othpl][j][k] ) != 0)
		    s = ratadd( s, ratmult( ratfromi(a), tmprow[j] )) ;
	    /* fill matrix entry with   -Ahat   */
	    intoM [ rowoffset + i] [ coloffset + k ] = 
		negpaymatrix ? ratneg(s) : s ;
	    }
	/* set matrix entry for  Ehat\T         */
	for (j = 0; j < irreddim[paytopl]; j++)
	    {    
	    s = ratfromi( redsfconstr [paytopl][j][i] );
	    intoM [ rowoffset + i] [ coloffset + redsfdim[othpl] + j] = 
		negpaymatrix ? s : ratneg(s) ;
	    }
	/* set LCP rhs entry  i  of  K\T A l = ahat     */
	s = ratfromi(0);
	for (j = 0; j < nseqs[othpl]; j++)
	    s = ratadd( s, ratmult( tmprow[j],
		ratfromi ( realplconst [othpl] [j] ) ) ) ;
	rhsvec[i] = negrhs ? ratneg(s) : s ;
	}

    /* second blockrow  of dim   irreddim [othpl]       */
    for (i = 0; i < irreddim [othpl]; i++)
	{
	/* set LCP matrix entry for  -Fhat              */
	for (j = 0; j < redsfdim[othpl]; j++)
	    {
	    s = ratfromi( redsfconstr [othpl][i][j] );
	    intoM [ rowoffset + redsfdim[paytopl] + i ] [ coloffset + j] =
		negpaymatrix ? ratneg(s) : s ;
	    }
	/* set LCP rhs entry  i  of  -fhat              */
	s = ratfromi( redsfrhs [othpl][i] );
	rhsvec[ redsfdim[paytopl] + i] = negrhs ? s : ratneg(s) ;
	}
    free(tmprow);
}


void rsflcp(void)
{
    genredsf(1);
    genredsf(2);
    setlcp( redsfdim[1] + irreddim[2] + redsfdim[2] + irreddim[1] );
    setrsfcomplconstr(1, lcpM, 0, redsfdim[1] + irreddim[2], 1,
		      rhsq, 1);
    setrsfcomplconstr(2, lcpM, redsfdim[1] + irreddim[2], 0, 1,
		      rhsq + redsfdim[1] + irreddim[2], 1);
}

