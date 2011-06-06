/* leaves.c
 * 2 June 2000
 */

#include <stdio.h>	/* fgets		*/
#include <string.h>	/* strstr, strlen	*/
#include <stdlib.h>	/* atoi        		*/

#include "rat.h"
#include "treedef.h"
#include "seqform.h"	/* outbehavstrat	*/
#include "sfnf.h"	/* realplan, allocrealplan()	*/

#include "leaves.h"

/* test if input properly processed by generating what should be an echo */
static void testreadequil (int docuseed)
{
    printf("BEQ>%4d<1>", docuseed);
    outbehavstrat(1, realplan[1], 0);
    printf(" <2>");
    outbehavstrat(2, realplan[2], 1);
}

void leavesfrominput (void)
{
    char line[MAXLINELENGTH];
    char *linepos ;
    int docuseed ;

    allocrealplan(realplan);
    printf("Processing input lines:\n");
    while ( fgets(line, MAXLINELENGTH, stdin) != NULL )
        if ( (linepos = strstr(line, "BEQ>") ) != NULL)
            {
	    docuseed = atoi( linepos + 4) ;
	    if (oklinetobehavprob (linepos) )
	        {
		behavtorealprob(0);
		behavtorealprob(1);
		behavtorealprob(2);
		realplanfromprob(1, realplan[1]);
		realplanfromprob(2, realplan[2]);
		/*
		testreadequil (docuseed) ;
		*/
		outleavesreached (docuseed);
		}
	    }
    freerealplan(realplan);
}

Bool oklinetobehavprob (char *line)
{
    int pl;
    char s[MAXSTRL] ;
    char *currpos ;
    Move c;

    for (pl = 1; pl < PLAYERS; pl++)
        {
	/* identify the substring "<1>" for player 1,  "<2>" for player 2 */
	sprintf(s, "<%d>", pl);
	currpos = strstr(line, s) ;
	if (currpos == NULL)   /* not a line in proper format	*/
	    return 0; 
        line = currpos + 4 ;

        for (c = firstmove[pl] + 1 ; c < firstmove[pl+1]; c++)
	    {
            c->behavprob = ratfromi(0) ;
	    movetoa(c, pl, s) ; 
	    currpos = strstr(line, s) ;
	    if (currpos != NULL)
	        {
		int l = strlen(s) ;
		if ( currpos[ l ] == ':' )
		    {
		    c->behavprob.num = atoi(currpos + l + 1) ;
	            currpos = strstr(currpos + l + 2, "/") ;
		    c->behavprob.den = atoi(currpos + 1) ;
		    }
                else
		    c->behavprob = ratfromi(1) ;
		line = currpos ;
		}
	    }
        }
    return 1 ;
}

void outleavesreached (int docuseed)
{
    char s[MAXSTRL];
    Outcome o;
    Rat reachprob ;

    printf("SUP>%4d<>", docuseed);
    for (o = outcomes; o < lastoutcome ; o++ )
    	{
	reachprob = o->whichnode->defseq[0]->realprob ;
	reachprob = ratmult( reachprob, o->whichnode->defseq[1]->realprob) ;
	reachprob = ratmult( reachprob, o->whichnode->defseq[2]->realprob) ;
	if ( reachprob.num )  /* nonzero probability of reaching leaf	*/
	    {
	    rattoa(reachprob, s); 
	    printf(" %d~%s", o - outcomes, s) ;
	    }
	}
    printf("\n") ;
}

