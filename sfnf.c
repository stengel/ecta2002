/* sfnf.c
 * 12 July 2000
 * common routines for SF, RSF, NF
 */

#include <stdio.h>
#include <stdlib.h>
        /* free()       */ 

#include "alloc.h"
#include "rat.h"
#include "rataux.h"        /* ratscalarprod        */

#include "lemke.h"
#include "treedef.h"

#include "seqform.h"
#include "normform.h"
#include "rsf.h"

#include "sfnf.h"

Rat *realplan[PLAYERS];

void allocrealplan (Rat * realpl[PLAYERS])
{
    int pl;

    for (pl=0; pl < PLAYERS; pl++)
        realpl[pl] = TALLOC (nseqs[pl], Rat);
}

void freerealplan(Rat * realpl[PLAYERS])
{
    int pl;

    for (pl=0; pl < PLAYERS; pl++)
        free(realpl[pl]);
}

void behavtorealprob (int pl)
{
    Move c;
    Move lastmove = firstmove[pl+1];
    firstmove[pl]->realprob = ratfromi(1);  /* empty seq has probability 1  */
    for (c = firstmove[pl]+1; c < lastmove; c++)
        c->realprob = ratmult(c->behavprob, c->atiset->seqin->realprob);
}

Bool eqrealplans(int pl, Rat *rplan1, Rat *rplan2, Bool bcomplain)
{
    int i;
    char s[MAXSTRL];
    Bool isok = 1;

    for (i = 0; i < nseqs[pl]; i++)
        if (! ratiseq (rplan1[i], rplan2[i]) )
            {
            isok = 0;
            if (bcomplain)
                {
                seqtoa(firstmove[pl] + i, pl, s);
                printf ("Player %d, sequence %s has realprob ", pl, s);
                rattoa( rplan1[i], s);
                printf ("%s", s);
                rattoa( rplan2[i], s);
                printf (", but should be %s\n", s);
                }
            else
                break;
            }
    return isok ;
}

void payratmatcpy(Payvec ** frommatr, int plminusone, Bool bnegate,
        Bool btranspfrommatr, int nfromrows, int nfromcols, 
        Rat ** targetmatr, int targrowoffset, int targcoloffset)
{
    int i,j;
    for (i=0; i < nfromrows; i++)
        for (j=0; j < nfromcols; j++)
            if (btranspfrommatr)
                targetmatr[j + targrowoffset][i + targcoloffset]
                = bnegate ? ratneg(frommatr[i][j][plminusone]) : frommatr[i][j][plminusone] ;
            else 
                targetmatr[i + targrowoffset][j + targcoloffset]
                = bnegate ? ratneg(frommatr[i][j][plminusone]) : frommatr[i][j][plminusone] ;
}

void intratmatcpy(int ** frommatr, Bool bnegate,
        Bool btranspfrommatr, int nfromrows, int nfromcols, 
        Rat ** targetmatr, int targrowoffset, int targcoloffset)
{
    int i,j;
    for (i=0; i < nfromrows; i++)
        for (j=0; j < nfromcols; j++)
            if (btranspfrommatr)
                targetmatr[j + targrowoffset][i + targcoloffset]
                = ratfromi(bnegate ? -frommatr[i][j] : frommatr[i][j]);
            else
                targetmatr[i + targrowoffset][j + targcoloffset]
                = ratfromi(bnegate ? -frommatr[i][j] : frommatr[i][j]);
}

void covvector(int whichform)
{
    int i, j, dim1, dim2, offset;

    behavtorealprob(1);
    behavtorealprob(2);
    switch(whichform)
        {
        case SFORM:
            dim1 = nseqs[1];
            dim2 = nseqs[2];
            offset = dim1 + 1 + nisets[2] ;
            break;
        case NFORM:
            dim1 = nstrats[1];
            dim2 = nstrats[2];
            offset = dim1 + 1 ;
            behavtomixed(1);
            behavtomixed(2);
            break;
        case RSFORM:
            dim1 = redsfdim[1];
            dim2 = redsfdim[2];
            offset = dim1 + irreddim[2];
            fprintf(stderr, "RSF covvector not yet implemented.\n");
            exit(1);
            break;
        default:
            fprintf(stderr,
                "whichform %d for covvector illegal.\n", whichform);
            exit(1);
        }
    /* covering vector  = -rhsq */
    for (i = 0; i < lcpdim; i++)
        vecd[i] = ratneg(rhsq[i]);
    /* first blockrow += -Aq    */
    for (i = 0; i < dim1; i++)
        if (whichform == SFORM)
            for (j=0; j < dim2; j++)
                vecd[i] = ratadd( vecd[i], ratmult( lcpM [i] [offset + j],
                          (firstmove[2] + j)->realprob)) ;
        else if (whichform == NFORM)
            ratscalarprod(dim2, lcpM [i] + offset, mixedstrat[2], &(vecd[i]));
        /* RSF yet to be done*/  
    /* third blockrow += -B\T p */
    for (i = offset; i < offset + dim2; i++)
        if (whichform == SFORM)
            for (j=0; j < dim1; j++)
                vecd[i] = ratadd( vecd[i], ratmult( lcpM [i] [j],
                          (firstmove[1] + j)->realprob));
        else if (whichform == NFORM)
            ratscalarprod(dim1, lcpM [i], mixedstrat[1], &(vecd[i]));
        /* RSF yet to be done*/  
}


void showeq(int whichform, Bool bshortequil, int docuseed)
{
    int offset;

    switch(whichform)
        {
        case SFORM:
            offset = nseqs[1] + 1 + nisets[2] ;
            /*  non-sparse printing 
            printf("Equilibrium realization plan player 1:\n");
            outrealplan(1, solz);
            printf("Equilibrium realization plan player 2:\n");
            outrealplan(2, solz + offset);
            */
            if (bshortequil)
                printf("BEQ>%4d<1>", docuseed);
            else
                printf("......Equilibrium behavior strategies player 1, 2:\n");
            outbehavstrat(1, solz, !bshortequil);
            if (bshortequil)
                printf(" <2>");
            outbehavstrat(2, solz + offset, 1);
            break;
        case NFORM:
            offset = nstrats[1] + 1;
            if (bshortequil)
                printf("MEQ>%4d<1>", docuseed);
            else
                printf("______Equilibrium mixed strategies player 1, 2:\n");
            outmixed(1, solz, !bshortequil);
            if (bshortequil)
                printf(" <2>");
            outmixed(2, solz + offset, 1);
            break;
        case RSFORM:
            fprintf(stderr, "RSF not yet implemented.\n");
            exit(1);
            break;
        default:
            abort();
        }
}

