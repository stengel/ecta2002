/* interface.c
 * 4 June 2000
 * interface with Audet/Hansen enumeration algorithm
 */

#include <stdio.h>

#include "rat.h"
#include "treedef.h"    /* PLAYERS              */
#include "normform.h"   /* */
#include "seqform.h"    /* sfpay, sfconstr      */
#include "interface.h"

static FILE *fa, *fb, *fe, *ff, *fmatrices;

static void fileerr(char *s)
{
    fprintf(stderr, "Couldn't open file  %s  for writing, abort\n", s);
    exit (1);
}

/* bnf: Open files "matrices"
 * bsf: Open files "A", "B", "E", "F" 
 * for interfacing, abort if error
 */
static void openfiles(Bool bnf, Bool bsf)
{
    if (bnf)
        {
        fmatrices = fopen ("matrices", "w") ;
        if (fmatrices == NULL) 
            fileerr("matrices");
	}
    if (bsf)
	{
        fa = fopen ("A", "w") ;
        if (fa == NULL) 
            fileerr("A");
        fb = fopen ("B", "w") ;
        if (fb == NULL) 
            fileerr("B");
        fe = fopen ("E", "w") ;
        if (fe == NULL) 
            fileerr("E");
        ff = fopen ("F", "w") ;
        if (ff == NULL) 
            fileerr("F");
	}
}

/* bnf: Close files "matrices"
 * bsf: Close files "A", "B", "E", "F" 
 */
static void closefiles(Bool bnf, Bool bsf)
{
    if (bnf)
        fclose(fmatrices);
    if (bsf)
        {
	fclose(fa);
	fclose(fb);
	fclose(fe);
	fclose(ff);
	}
}

void interface(Bool bnf, Bool bsf)
{
    int i, j;

    openfiles(bnf, bsf) ;
    if (bnf)
        {
        fprintf(fmatrices, "%d  %d\n\n", nstrats[1], nstrats[2]) ;
    	/* payoffs player 1 */
        for (i=0; i<nstrats[1]; i++)
	    {
            for (j=0; j<nstrats[2]; j++)
                fprintf(fmatrices, " %8.4f",
                            rattodouble(nfpay[i][j][0]));
            fprintf(fmatrices, "\n");
	    }
	fprintf(fmatrices, "\n");
    	/* payoffs player 2 */
        for (i=0; i<nstrats[1]; i++)
	    {
            for (j=0; j<nstrats[2]; j++)
                fprintf(fmatrices, " %8.4f",
                            rattodouble(nfpay[i][j][1]));
            fprintf(fmatrices, "\n");
	    }
	}
    if (bsf)
        {
        /* sparse representation in files:  row, column, entry,
         * for nonzero entries
         */
        for (i=0; i<nseqs[1]; i++)
            {
            for (j=0; j<nseqs[2]; j++)
                {
                /* payoffs player 1 */
                if (sfpay[i][j][0].num)
                    fprintf(fa, "%3d %3d %8.4f\n",
                            i, j, rattodouble(sfpay[i][j][0]));
                /* payoffs player 2 */
                if (sfpay[i][j][1].num)
                    fprintf(fb, "%3d %3d %8.4f\n",
                        i, j, rattodouble(sfpay[i][j][1]));
                }
            }
        /* constraints player 1 */
        for (i = 0; i <= nisets[1]; i++)
            for (j=0; j<nseqs[1]; j++)
                if (sfconstr[1][i][j])
                    fprintf(fe, "%3d %3d %2d\n", i, j, sfconstr[1][i][j]);
        /* constraints player 2 */
        for (i = 0; i <= nisets[2]; i++)
            for (j=0; j<nseqs[2]; j++)
                if (sfconstr[2][i][j])
                    fprintf(ff, "%3d %3d %2d\n", i, j, sfconstr[2][i][j]);
        }
    closefiles(bnf, bsf) ;
}
