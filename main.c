/*  main.c
 *  tracing procedure algorithm
 *  1 May 2001: extend with example
 *
 *  options:
 *      -A  #   accuracy for prior generation, maximally  MAXACCURACY
 *      -b      both normal and sequence form
 *      -c      complementary pivoting steps shown
 *      -d      degeneracy statistics in lexmin ratio test
 *      -e      give equilibrium  (when -m #)
 *      -E      equilibrium leaves from input only
 *      -i      interface with Audet/Hansen enumeration, NF/SF as chosen
 *      -g      print raw game data
 *      -G file interface with GAMBIT via file (re-writes for each game,
 *              so don't use with -m )
 *      -l  #   bintree with # levels (e.g. #=3, any number in
 *              MINLEVEL .. MAXLEVEL is allowed)
 *              default (no -l  option): solve tracingexample.
 *		negative number: -1      solve forward induction example
 *      -m  #   process multiple games, requires  -l option
 *              -m 1 implies quiet mode
 *      -M  #   multiple priors per game
 *      -n      compute with normal form (default: sequence form)
 *      -o      output LCP
 *      -O      output prior
 *      -p -- # # # # ...     replace the payoffs at leaves by # # # # 
 *              (leaf 0/pl 1, leaf 0/pl 2, leaf 1/pl 1, leaf 1/pl 2 ...)
 *              which should be NEGATIVE since invoked AFTER re-normalization;
 *              this must be the LAST option on the command line
 *      -r      compute with RSF         (not yet implemented)
 *      -s  #   payoff seed 
 *      -S  #   prior seed 
 *      -t      tableaus at every pivoting step
 */

#include <stdio.h>
#include <string.h>	/* strcpy		*/
#include <stdlib.h>	/* atoi(), free()       */
#include <ctype.h>	/* isprint()            */
/* #include <unistd.h> */
#include "getopt.h"
	/* getopt(), optarg, optopt, optind             */
#include <time.h>
	/* clock_t, clock(), CLOCKS_PER_SEC     	*/
#include <limits.h>
        /* INT_MAX,  INT_MIN    */

#include "alloc.h"
#include "rat.h"
#include "lemke.h"
#include "mp.h"         /* record_digits, DIG2DEC()     */
#include "treedef.h"
#include "treegen.h"
#include "sfnf.h"
#include "seqform.h"
#include "normform.h"
#include "prior.h"
#include "leaves.h"
#include "gambit.h"
#include "interface.h"

#define MINLEVEL 1      
#define MAXLEVEL 10
#define FILENAMELENGTH 50
#define CLOCKUNITSPERSECOND 1000.0
#define SCLOCKUNITS "millisecs"
#define WHICHFORMS 2      /* currently only SF (0), NF (1)      */
#define REPEATHEADER 20   /* repeat header if more games than this */

/* global variables for generating and documenting computation  */
static  Flagsprior    fprior;
static  Bool boutlcp = 0;       /* output LCP       (-o option) */
static  Bool boutprior = 0;     /* output prior     (-O option) */
static  Bool bcomment = 0;      /* complementary pivoting steps */
static  Bool bequil = 1;        /* output equilibrium           */
static  Bool bshortequil = 0;   /* output equilibrium shortly   */
static  Bool bleavesonly = 0;  	/* equilibrium leaves only   	*/
static  Bool binterface = 0;  	/* interface with enumeration	*/
static  Bool bgambit = 0;  	/* interface with gambit	*/
/* GAMBIT interface file, option parameter	*/
static  char gintfname[FILENAMELENGTH] = "dummyname" ;
static  Flagsrunlemke flemke;

static  int timeused [WHICHFORMS], sumtimeused [WHICHFORMS];
static  int pivots   [WHICHFORMS], sumpivots   [WHICHFORMS];
static  int lcpsize  [WHICHFORMS];
static  int mpdigits [WHICHFORMS], summpdigits [WHICHFORMS];
static  int eqsize [PLAYERS] [WHICHFORMS], sumeqsize [PLAYERS] [WHICHFORMS];
static  Bool agreenfsf [PLAYERS] ;

/* initialize   sumarray[WHICHFORMS]  to zero   */
void settozero(int *sumarray)
{
    int i;
    for (i = 0; i < WHICHFORMS; i++)
	sumarray[i] = 0;
}

/* returns processor SCLOCKUNITS since the last call to
 * stopwatch() and prints them to stdout if  bprint==1
 */
int stopwatch(Bool bprint)
{
    static clock_t time;
    double x;

    x = (double) (clock()) - (double) time;
    if (x < 0)
    	x += 2 * (double) INT_MAX;
    x /= ((double) CLOCKS_PER_SEC / CLOCKUNITSPERSECOND) ;
    if (bprint)
	printf("time elapsed [%s] %4.0f\n", SCLOCKUNITS, x);
    time = clock();
    return (int) x;
}

/* informs about tree size              */
void infotree()
{
    int pl;
    printf("\nGame tree has %d nodes, ", lastnode - root);
    printf("of which %d are terminal nodes.\n", lastoutcome - outcomes);
    for (pl = 0; pl < PLAYERS; pl++)
	{
	printf("    Player %d has ", pl);
	printf("%3d information sets, ", firstiset[pl+1] - firstiset[pl]);
	printf("%3d moves in total\n", firstmove[pl+1] - firstmove[pl] - 1); 
	}
}

/* informs about normal form, compute and set  lcpsize[NFORM]      */
void infonf()
{
    int dim[PLAYERS];
    int pl;
    
    for (pl = 1; pl < PLAYERS; pl++)
	dim[pl] = numstratsnfpre(pl);
    lcpsize [NFORM] = dim[1] + dim[2] + 2 ;
    printf("Normal form LCP dimension is %d\n", lcpsize [NFORM] );
    for (pl = 1; pl < PLAYERS; pl++)
	{
	printf("    Player %d has ", pl);
	printf("%5d RNF pure strategies\n", dim[pl]); 
	}
}

/* informs about sequence form, set  lcpsize[SFORM]    */
void infosf()
{
    int pl;
    
    lcpsize [SFORM] = nseqs[1] + nisets[2]+1 + nseqs[2] + nisets[1]+1 ;
    printf("Sequence form LCP dimension is %d\n", lcpsize [SFORM] );
    for (pl = 1; pl < PLAYERS; pl++)
	{
	printf("    Player %d has ", pl);
	printf("%3d sequences, ", nseqs[pl]);
	printf("subject to %3d constraints\n", nisets[pl]+1); 
	}
} 

/* give header columns for result information via  inforesult(...)      */
void inforesultheader (Bool bsf, Bool bnf)
{
    printf("PRIOR/PAY| ");
    if (bnf)
	{
        printf("NORMAL FORM          support");
        if (bsf)
            printf(" |agrees| ");
	}
    if (bsf)
        printf("SEQUENCE FORM        mixiset");
    printf("\n");
    printf("Seed/seed| ");
    if (bsf)
	{
	printf("pivot %%n [secs] digs pl1 pl2");
	if (bnf)
            printf(" |pl1,2 | ");
	}
    if (bnf)
	printf("pivot %%n [secs] digs pl1 pl2");
    printf("\n");
}

/* info about results for game with  priorseed  and  (payoff) seed */
void inforesult (Bool bsf, Bool bnf, int priorseed, int seed)
{
    char formatstring[] = "%4d %3.0f %6.2f  %3d %3d %3d" ;
    printf("%4d/%4d| ", priorseed, seed);
    if (bnf)
	{
	printf(formatstring, pivots [NFORM], 
            (double) pivots [NFORM]*100.0 / (double) lcpsize [NFORM],
	    (double) timeused [NFORM] / CLOCKUNITSPERSECOND,
	    mpdigits [NFORM], eqsize [1] [NFORM], eqsize [2] [NFORM]);
        if (bsf)
            printf(" | %s  %s | ",  agreenfsf[1] ? "Y" : "N",
		    agreenfsf[2] ? "Y" : "N");
	}
    if (bsf)
	printf(formatstring, pivots [SFORM], 
            (double) pivots [SFORM]*100.0 / (double) lcpsize [SFORM],
	    (double) timeused [SFORM] / CLOCKUNITSPERSECOND,
	    mpdigits [SFORM], eqsize [1] [SFORM], eqsize [2] [SFORM]);
    printf("\n");
}

/* summary info about results for  m  games     */
void infosumresult (Bool bsf, Bool bnf, int m)
{
    double mm = (double) m;
    char formatstring[] = "%6.1f %3.0f %6.2f %4.1f %3.1f %3.1f" ;
    
    printf("---------| AVERAGES over  %d  games:\n", m);
    if (m > REPEATHEADER)
        inforesultheader (bsf, bnf);
    printf("         ");
    if (bnf)
	{
	printf(formatstring, (double) sumpivots [NFORM]/ mm, 
            (double) sumpivots [NFORM]*100.0 /
                (double) (lcpsize [NFORM] * mm),
	    (double) sumtimeused [NFORM] / (CLOCKUNITSPERSECOND * mm),
	    (double) summpdigits [NFORM] / mm, 
            (double) sumeqsize [1] [NFORM] / mm, 
            (double) sumeqsize [2] [NFORM] / mm);
        if (bsf)
            printf("        ");
	}
    if (bsf)
	printf(formatstring, (double) sumpivots [SFORM]/ mm, 
            (double) sumpivots [SFORM]*100.0 /
                (double) (lcpsize [SFORM] * mm),
	    (double) sumtimeused [SFORM] / (CLOCKUNITSPERSECOND * mm),
	    (double) summpdigits [SFORM] / mm, 
            (double) sumeqsize [1] [SFORM] / mm, 
            (double) sumeqsize [2] [SFORM] / mm);
    printf("\n");
}

/* process game for evaluation
 * for comparison:  call first for  NF  then  SF
 * bnf:  NF is processed, compare result with SF result
 * docuseed:  what seed to output for short equilibrium output
 * realplan[][]  must be allocated
 */
void processgame (int whichform, Bool bnf, int docuseed)
{
    int equilsize;
    int offset;
    int pl;

    if (whichform == NFORM)
        {
        if (bcomment)
            printf("Generating and solving normal form.\n");
        nflcp();
        }
    else if (whichform == SFORM)
        {
        if (bcomment)
            printf("Generating and solving sequence form.\n");
	sflcp();
        }
    else
	abort();        /* no RNF yet   */

    covvector(whichform);
    if (boutlcp)
	outlcp();
    stopwatch(0);
    record_digits = 0;
    runlemke(flemke);
    sumtimeused [whichform] += timeused [whichform] =
	stopwatch(0);
    sumpivots [whichform] += pivots [whichform] =
	pivotcount;
    summpdigits [whichform] += mpdigits [whichform] =
	DIG2DEC(record_digits);
    /* equilibrium size     */
    offset = 0;
    for (pl = 1; pl < PLAYERS; pl++)
	{
	if (whichform == NFORM)
	    {
	    equilsize = supportsize(pl, solz + offset);
	    mixedtorealplan(pl, solz + offset, realplan[pl]);
	    /* the next is  offset  for player 2 */
	    offset = nstrats[1] + 1;
	    }
	else if (whichform == SFORM)
	    {
	    equilsize = propermixisets(pl, solz + offset);
            if (bnf)
                agreenfsf[pl] = eqrealplans(pl, solz + offset,
                                   realplan[pl], bcomment);
	    /* the next is  offset  for player 2 */
	    offset = nseqs[1] + 1 + nisets[2] ;
	    }
	else
	    abort();        /* no RNF yet   */
	sumeqsize [pl] [whichform] +=
	    eqsize [pl] [whichform] = equilsize ;
	}
    if (bequil)
	showeq (whichform, bshortequil, docuseed); 
    if (bgambit && whichform == SFORM)
	gambshoweq(); 
}


int main(int argc, char *argv[])
{
    int levels = 0;     /* which game to process,   (-l option)
			 * 0:  tracing example BvS/Elzen/Talman,
			 * -1: forward induction example
			 * MINLEVEL..MAXLEVEL:  bintree
			 */
  
    int multiplegames  = 0;      /* parameter for    -m option  */
    int multipriors = 0;         /* parameter for    -M option  */
    int seed  = 0;      /* payoff seed for bintree  (-s option) */
    int newpayoffs = 0; /* number of payoffs to be replaced (-p)*/
    int *newp1, *newp2; /* arrays for entering new payoffs      */

    /* whichform currently not used (later for RSF)             */
    int whichform = 0;  /* 0: SF, 1: NF, 2: RSF                 */

    Bool bheadfirst = 0;/* headers first (multiple games)       */
    Bool bsf = 1;       /* Y/N process SF, default              */
    Bool bnf = 0;       /* Y/N process NF           (-n option) */
    Bool bgame = 0;     /* output the raw game tree (-g option) */
    int c ;

    flemke.maxcount   = 0;
    flemke.bdocupivot = 0;
    flemke.binitabl   = 0;
    flemke.bouttabl   = 0;
    flemke.boutsol    = 0;
    flemke.binteract  = 0;
    flemke.blexstats  = 0;

    fprior.seed       = 0 ;
    fprior.accuracy   = DEFAULTACCURACY ;

    /* parse options    */
    while ( (c = getopt (argc, argv, "A:bcdeEgG:il:m:M:noOprs:S:t")) != -1)
	switch (c)
	    {
	    int x; 
	    case 'A':
		x = atoi(optarg);
		if (x <= MAXACCURACY &&  x > 0 )
		    fprior.accuracy = x ;
                else
		    {
		    fprintf(stderr, "Entered accuracy %d for prior ", x);
		    fprintf(stderr, "not in 1..%d, not done.\n", MAXACCURACY);
	            }	
		break;
            case 'b':
		bnf = 1;
		bsf = 1;
		break;
            case 'c':
                bcomment = 1;
		break;
	    case 'd':
		flemke.blexstats   = 1;
		break;
            case 'e':
                bshortequil = 1;
		break;
            case 'E':
                bleavesonly = 1;
		break;
	    case 'g':
		bgame = 1;
		break;
	    case 'G':
		bgambit = 1;
		strcpy(gintfname, optarg);
		break;
	    case 'i':
		binterface = 1;
		break;
	    case 'l':
		levels = atoi(optarg);
		if (levels >= MINLEVEL && levels <= MAXLEVEL)
		    break;
		if (levels == -1)

		    break; 
		fprintf(stderr, "Binary tree level %d ", levels);
		fprintf(stderr, "not in range %d .. %d, not done.\n", 
			MINLEVEL, MAXLEVEL);
		return 1;
	    case 'm':
		multiplegames = atoi(optarg);
		break;
	    case 'M':
		multipriors = atoi(optarg);
		break;
	    case 'n':   
		bnf = 1;
		bsf = 0;
		break;
	    case 'o':
                boutlcp = 1;
		break;
	    case 'O':
                boutprior = 1;
		break;
	    case 'p':
		newpayoffs = 1;
		break;
	    case 'r':
		whichform = RSFORM ;
		break;
	    case 's':
		seed = atoi(optarg);
		break;
	    case 'S':
		x = atoi(optarg);
		if ( x > 0 )
		    fprior.seed = x ;
                else
		    {
		    fprintf(stderr, "Entered Seed %d for prior ", x);
		    fprintf(stderr, "must be positive, ignored.\n");
	            }	
		break;
            case 't':
		flemke.bouttabl   = 1;
                flemke.binitabl   = 1;
		break;
	    case '?':
		if (isprint (optopt))
		    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
		else
		    fprintf (stderr,
			     "Unknown option character `\\x%x'.\n",
			      optopt);
		return 1;
	    default:
		abort ();
	    }
    /* options have been input, amend extras	*/
    if (newpayoffs)
	{
	newpayoffs = (argc - optind)/2;
	printf("Entering %d new payoff pair(s).\n", newpayoffs);
	}
    if (newpayoffs)
	{
	int i;
	newp1 = TALLOC(newpayoffs, int);
	newp2 = TALLOC(newpayoffs, int);
	for (i = 0, c = optind; i < newpayoffs; i++) 
	    {
	    newp1[i] = atoi(argv[c++]);
	    newp2[i] = atoi(argv[c++]);
	    }
	} 
    if (multiplegames)
	{
	if (levels == 0)
	    {
	    printf("Multiple games only for binary tree, ");
	    printf("but will change equilibrium output.\n");
	    multiplegames = 0;
            bequil = 0 ;
	    }
	if (newpayoffs)
	    {
	    printf("Multiple games not with changed payoffs, ");
	    printf("but will change equilibrium output.\n");
	    multiplegames = 1;
	    }
	/* else */
            {
            bequil = bshortequil;
            bheadfirst = !bcomment && !boutlcp && !flemke.bouttabl &&
                !flemke.blexstats && !bgame ;
            }
	}
    if (multipriors > 0)
	{
	/* this would exclude the centroid for multiple priors
        if ( fprior.seed == 0)
            fprior.seed = 1 ; 
	*/
        }
    else
        multipriors = 1 ;
    if (bcomment)
	    {
            flemke.bdocupivot = 1;
            flemke.boutsol    = 1;
	    }

    /* options are parsed and flags set */
    /* document the set options         */
    printf("Options chosen,              [ ] = default:\n");
    printf("    normal form            %s [N],  option -n\n",
	    bnf ? "Y" : "N");
    printf("    sequence form          %s [Y],  option -b : both NF and SF\n",
	    bsf ? "Y" : "N");
    printf("    interface to enumerate %s [N],  option -i , NF/SF as above\n",
            binterface ? "Y" : "N");
    printf("    GAMBIT interface       %s [N],  ", bgambit ? "Y" : "N");
    printf("option -G file, SF only\n");
    printf("    levels binary tree    %2d [0],  option -l # ", levels);
    printf("(default 0: simple example)\n");
    printf("    seed payoffs         %3d [0],  option -s #\n", seed);
    printf("    multiple games       %3d [0],  option -m # ", multiplegames);
    printf("(no equilibria unless option -e)\n");
    printf("    equilibrium one line   %s [N],  option -e\n",
            bshortequil ? "Y" : "N");
    printf("    equil leaves < stdin   %s [N],  option -E\n",
            bleavesonly ? "Y" : "N");
    printf("    payoffs new          %3d [0],  option -p -- # # ... #\n",
	    newpayoffs);
    printf("    Multiple priors     %4d [1],  option -M #\n", multipriors);
    printf("    Accuracy prior      %4d [%d], option -A #\n",
	    fprior.accuracy, DEFAULTACCURACY);
    printf("    Seed prior           %3d [0],  ",
            fprior.seed);
    printf("option -S # (default 0: centroid)\n");
    printf("    Output prior           %s [N],  option -O\n",
            boutprior ? "Y" : "N");
    printf("    game output            %s [N],  option -g\n",
	    bgame ? "Y" : "N");
    printf("    comment LCP pivs & sol %s [N],  option -c\n",
            bcomment ? "Y" : "N");
    printf("    output LCP             %s [N],  option -o\n",
            boutlcp ? "Y" : "N");
    printf("    degeneracy statistics  %s [N],  option -d\n",
	    flemke.blexstats ? "Y" : "N");
    printf("    tableaus               %s [N],  option -t\n",
	    flemke.bouttabl ? "Y" : "N");
                    
    if (levels == 0)
	{
	printf("Solving example from BvS/Elzen/Talman\n");
	tracingexample();
	}
    else
    if (levels == -1)
	{
	printf("Forward induction example\n");
	forwardexample();
	}
    else
	createbintree (levels, seed);

    genseqin();  
    autoname();
    maxpayminusone(bcomment);
    if (newpayoffs)
	{
	int i; 
	if (newpayoffs > lastoutcome-outcomes)
	    {
	    newpayoffs = lastoutcome-outcomes;
	    printf("Only the %d outcomes get new payoffs.\n", newpayoffs);
	    }
	for (i = 0; i < newpayoffs; i++) 
	    {
	    outcomes[i].pay[0] = ratfromi(newp1[i]);
	    printf ("Outcome %2d pay1: %3d, ", i, newp1[i]);
	    outcomes[i].pay[1] = ratfromi(newp2[i]);
	    printf ("pay2: %3d\n", newp2[i]);
	    }
	printf("re-normalize payoffs again.\n");
        maxpayminusone(bcomment);
	free(newp1);
	free(newp2);
	}

    /* game tree is defined, give headline information  */
    infotree();

    if (bleavesonly)
        leavesfrominput();
    else
        {
        if (bnf)
	    infonf();
        if (bsf)
	    infosf();
    
        /* process games                    */
        if (multiplegames == 0)
            multiplegames = 1;	/* simplify counting	*/
        {   
	    int gamecount = 0;
	    int startprior = fprior.seed ;
    
	    allocrealplan(realplan);
            if (bheadfirst) /* otherwise the header is garbled by LCP output */
                inforesultheader (bsf, bnf);
	    while(1)        /* process multiple games       */
	        {
	        int priorcount ;
                /* multiple priors 	*/
		if (bgambit)
		    gambopenfile(gintfname);
	        for (priorcount = 0; priorcount < multipriors; priorcount++)
	            {
	            genprior(fprior);
                    if (bgame)
        	        rawtreeprint();
                    if (boutprior)
        	        outprior();
	            if (bnf)
                        processgame(NFORM, bnf, seed + gamecount); 
	            if (bsf)
                        processgame(SFORM, bnf, seed + gamecount); 
                    if ( ! bheadfirst )
                        inforesultheader (bsf, bnf);
	            inforesult (bsf, bnf, fprior.seed, seed + gamecount);
                    fprior.seed++ ;
	            }
		if (binterface)
		    {
		    interface(bnf, bsf);
		    binterface = 0;
		    }
		if (bgambit)
		    gambclosefile() ;
	        /* next game */
	        gamecount++ ;
	        if (gamecount >= multiplegames)
		    break;
                if (multipriors > 1)
	    	    printf("\n") ;
	        createbintree (levels, seed + gamecount);
	        genseqin();  
	        autoname();
                maxpayminusone(bcomment);
	        fprior.seed = startprior ;
	        }           /* end of processing multiple games     	 */
	    if (multipriors * multiplegames > 1)	/* give averages */
	        infosumresult (bsf, bnf, multipriors * multiplegames);
	    freerealplan(realplan);
        }
        }
    return 0;
}
