/* gambit.c
 * interfaces to the GAMBIT program
 * 12 June 2000
 */

#include <stdio.h>

#include "rat.h"
#include "treedef.h"
#include "lemke.h"	/* solz	*/
#include "seqform.h"

#include "gambit.h" 

static FILE *fgintf;		/* GAMBIT interface file		*/
static Bool bprefixcomma ; 	/* print a comma before the next equil	*/

void gambopenfile(char *filename) 
{
    fgintf = fopen(filename, "w") ;
    if (fgintf == NULL)
       {
       fprintf(stderr, "Couldn't open file  %s  for writing, abort\n",
      	       filename);
       exit (1);
       }
    /*	comment for Include files only
    fprintf(fgintf, "// GAMBIT interface file \"%s\" produced by bintree\n",
               filename);
    */
    fprintf(fgintf, "{ ");
    bprefixcomma = 0 ; 
}

void gambclosefile(void) 
{
    fprintf(fgintf, " }\n");
    fclose(fgintf) ;
}

void gambshoweq(void) 
{
    int offset = nseqs[1] + 1 + nisets[2] ;
    
    if (bprefixcomma)
        fprintf(fgintf, ",\n\n");
    bprefixcomma = 1 ;
    fprintf(fgintf, "{ ");
    glistbehavstrat(fgintf, 1, solz) ;
    fprintf(fgintf, ", ");
    glistbehavstrat(fgintf, 2, solz + offset) ;
    fprintf(fgintf, " }");
}

void glistbehavstrat(FILE *fp, int pl, Rat *rplan)
{
    char s[MAXSTRL];
    int i;
    Move c;
    Iset h;
    Rat rinprob, bprob;

    fprintf (fp, "{ ");    /*	opening brace for behavstrat	*/
    for (h = firstiset[pl]; h < firstiset[pl+1]; )
        /* to work even for player without isets */
        {
	fprintf (fp, "{");   
	rinprob = rplan[ h->seqin - firstmove[pl] ];
        for (c = h->move0, i=0; ; c++)	/* at least one move necessary	*/
            {
            if ( rinprob.num == 0)
                {
		bprob.num = 1;
		bprob.den = h->nmoves;
		}
	    else
                bprob = ratdiv( rplan[ c - firstmove[pl]], rinprob);
	    rattoa(bprob, s);
	    fprintf(fp, "%s", s);
	    i++ ;
	    if (i == h->nmoves) 
		break ;
	    fprintf(fp, ",");
            }
	fprintf (fp, "}");    
        h++ ;
        if (h == firstiset[pl+1])
            break ;
        fprintf(fp, ", ");
	}
    fprintf (fp, " }");    /*	closing brace for behavstrat	*/
}

