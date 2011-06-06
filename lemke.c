/* lemke.c
 * 13 July 2000
 * LCP solver
 */

#include <stdio.h>
#include <stdlib.h>
	/* free()       */ 
#include <string.h>
	/*  strcpy      */

#include "alloc.h"
#include "col.h"
#include "rat.h"

#include "lemke.h"

#include "mp.h"

/* used for tableau:    */
#define Z(i) (i)
#define W(i) (i+n)
	/* VARS   = 0..2n = Z(0) .. Z(n) W(1) .. W(n)           */
	/* ROWCOL = 0..2n,  0 .. n-1: tabl rows (basic vars)    */
	/*                  n .. 2n:  tabl cols  0..n (cobasic) */
#define RHS  (n+1)                   /*  q-column of tableau    */
#define TABCOL(v)  (bascobas[v]-n)   
	/*  v in VARS, v cobasic:  TABCOL(v) is v's tableau col */
	/*  v  basic:  TABCOL(v) < 0,  TABCOL(v)+n   is v's row */

/* LCP input    */
Rat **lcpM;
Rat *rhsq; 
Rat *vecd; 
int lcpdim = 0; /* set in setlcp                */
static int n;   /* LCP dimension as used here   */

/* LCP result   */
Rat  *solz; 
int  pivotcount;

/* tableau:    */
static  mp **A;                 /* tableau                              */
static  int *bascobas;          /* VARS  -> ROWCOL                      */
static  int *whichvar;          /* ROWCOL -> VARS, inverse of bascobas  */
                                                                        
/* scale factors for variables z
 * scfa[Z(0)]   for  d,  scfa[RHS] for  q
 * scfa[Z(1..n)] for cols of  M
 * result variables to be multiplied with these
 */
static  mp *scfa;

static  mp det;                         /* determinant                  */

static  int *lextested, *lexcomparisons;/* statistics for lexminvar     */

static  int *leavecand;
    /* should be local to lexminvar but defined globally for economy    */


/*------------------ error message ----------------*/
void errexit (char *info)
{
    fflush(stdout);
    fprintf(stderr, "Error: %s\n", info);
    fprintf(stderr, "Lemke terminated unsuccessfully.\n");
    exit(1);
}

/* declares */
void assertbasic (int v, const char *info);

/*------------------ memory allocation -------------------------*/
void setlcp (int newn)
{
    if (newn < 1 || newn > MAXLCPDIM)
	{
	fprintf(stderr, "Problem dimension  n= %d not allowed.  ", newn);
	fprintf(stderr, "Minimum  n  is 1, maximum %d.\n", MAXLCPDIM);
	exit(1);
	}
    if (lcpdim > 0)             /* free previously used space   */
	{
	FREE2(lcpM, lcpdim); 
	free(rhsq);
	free(vecd);
	free(solz);
	FREE2(A, lcpdim); 
	free(scfa);
	free(bascobas);
	free(whichvar);
	free(leavecand);
	}
    n = lcpdim = newn;
    /* LCP input/output data    */
    T2ALLOC (lcpM, n, n, Rat);
    rhsq = TALLOC(n, Rat);
    vecd = TALLOC(n, Rat);
    solz = TALLOC(n, Rat);
    /* tableau          */
    T2ALLOC (A, n, n+2, mp);
    scfa = TALLOC (n+2, mp);
    bascobas = TALLOC(2*n+1, int);
    whichvar = TALLOC(2*n+1, int);
    lextested = TALLOC(n+1, int);
    lexcomparisons = TALLOC(n+1, int);
    leavecand = TALLOC(n, int);
	/* should be local to lexminvar but allocated here for economy  */
    /* initialize all LCP entries to zero       */
    {
	int i,j;
	Rat zero = ratfromi(0);
	for (i=0; i<n; i++)
	    {
	    for (j=0; j<n; j++)
		lcpM [i] [j] = zero;
	    vecd [i] = rhsq [i] = zero;
	    }
    }
} /* end of  setlcp(newn)       */

/* asserts that  d >= 0  and not  q >= 0  (o/w trivial sol) 
 * and that q[i] < 0  implies  d[i] > 0
 */
void isqdok (void)
{
    int i;
    int isqpos = 1;
    for (i=0; i<n; i++)
	{
	if (vecd[i].num < 0)
	    {
	    fprintf(stderr, "Covering vector  d[%d] = %d/%d negative\n",
		    i+1, vecd[i].num, vecd[i].den);
	    errexit("Cannot start Lemke.");
	    }
	else if (rhsq[i].num < 0)
	    {
	    isqpos = 0;
	    if (vecd[i].num == 0)
		{
		fprintf(stderr, "Covering vector  d[%d] = 0  ", i+1);
		fprintf(stderr, "where  q[%d] = %d/%d  is negative.\n",
			i+1, rhsq[i].num, rhsq[i].den);
		errexit("Cannot start Lemke.");
		}
	    }
	}   /* end of  for(i=...)   */
    if (isqpos)
	{
	printf("No need to start Lemke since  q>=0. ");
	printf("Trivial solution  z=0.\n");
	exit(0);
	}
}       /* end of  isqdok()     */

/* ------------------- tableau setup ------------------ */
void inittablvars (void)
	/* init tableau variables:                      */
	/* Z(0)...Z(n)  nonbasic,  W(1)...W(n) basic    */
{
    int i;
    for (i=0; i<=n; i++)
	{
	bascobas[Z(i)] = n+i;
	whichvar[n+i]  = Z(i);
	}
    for (i=1; i<=n; i++)
	{
	bascobas[W(i)] = i-1;
	whichvar[i-1]  = W(i);
	}
}       /* end of inittablvars()        */

void filltableau (void)
	/* fill tableau from  M, q, d   */
{
    int i,j;
    int den, num;
    mp tmp, tmp2, tmp3;
    for (j=0; j<=n+1; j++)
	{
	/* compute lcm  scfa[j]  of denominators for  col  j  of  A         */
	itomp(ONE, scfa[j]);
	for (i=0; i<n; i++)
	    {
	    den = (j==0) ? vecd[i].den :
		  (j==RHS) ? rhsq[i].den : lcpM[i][j-1].den ;
	    itomp(den, tmp);
	    lcm(scfa[j], tmp);
	    }
	/* fill in col  j  of  A    */
	for (i=0; i<n; i++)
	    {
	    den = (j==0) ? vecd[i].den :
		  (j==RHS) ? rhsq[i].den : lcpM[i][j-1].den ;
	    num = (j==0) ? vecd[i].num :
		  (j==RHS) ? rhsq[i].num : lcpM[i][j-1].num ;
		/* cols 0..n of  A  contain LHS cobasic cols of  Ax = b     */
		/* where the system is here         -Iw + dz_0 + Mz = -q    */
		/* cols of  q  will be negated after first min ratio test   */
	    /* A[i][j] = num * (scfa[j] / den),  fraction is integral       */
	    itomp (den, tmp);
	    copy (tmp3, scfa[j]);
	    divint(tmp3, tmp, tmp2);        /* divint modifies 1st argument */
	    itomp (num, tmp);
	    mulint(tmp2, tmp, A[i][j]);
	    }
	}   /* end of  for(j=...)   */
    inittablvars();
    itomp (ONE, det);
    changesign(det);
}       /* end of filltableau()         */

/* ---------------- output routines ------------------- */
void outlcp (void)
	/* output the LCP as given      */
{
    int i,j ;
    Rat a;
    char s[LCPSTRL];

    printf("LCP dimension: %d\n", n);
    colset(n + 2);
    for (j=0; j<n; j++)
	colpr("");
    colpr("d");
    colpr("q");
    colnl();
    for (i=0; i<n; i++)
	{
	for (j=0; j<n; j++)
	    {
	    a = lcpM [i] [j];
	    if (a.num == 0)
		colpr(".");
	    else
		{
                rattoa(a, s);
		colpr(s);
		}
	    }
	rattoa( vecd [i], s);
	colpr(s);
	rattoa( rhsq [i], s);
	colpr(s);
	}
    colout();
}

int vartoa(int v, char s[])
	/* create string  s  representing  v  in  VARS,  e.g. "w2"    */
	/* return value is length of that string                      */
{
    if (v > n)
	return sprintf(s, "w%d", v-n);
    else
	return sprintf(s, "z%d", v);
}


void outtabl (void)
	/* output the current tableau, column-adjusted                  */
{
    int i, j;
    char s[INFOSTRINGLENGTH];
    char smp [DIG2DEC(MAX_DIGITS)+2];       /* string to print  mp  into    */
    mptoa (det, smp);
    printf("Determinant: %s\n", smp);                
    colset(n+3);
    colleft(0);
    colpr("var");                   /* headers describing variables */
    for (j=0; j<=n+1; j++)
	{
	if (j==RHS)
	    colpr("RHS");
	else
	    {
	    vartoa(whichvar[j+n], s);
	    colpr(s);
	    } 
	}
    colpr("scfa");                  /* scale factors                */
    for (j=0; j<=n+1; j++)
	{
	if (j==RHS)
	    mptoa(scfa[RHS], smp);
	else if (whichvar[j+n] > n) /* col  j  is some  W           */
	    sprintf(smp, "1");
	else                        /* col  j  is some  Z:  scfa    */
	    mptoa( scfa[whichvar[j+n]], smp);
	colpr(smp);
	}
    colnl();
    for (i=0; i<n; i++)             /* print row  i                 */
	{
	vartoa(whichvar[i], s);
	colpr(s);
	for (j=0; j<=n+1; j++)
	    {
	    mptoa( A[i][j], smp);
	    if (strcmp(smp, "0")==0)
		colpr(".");
	    else
		colpr(smp);
	    }
	}
    colout();
    printf("-----------------end of tableau-----------------\n");
}       /* end of  outtabl()                                    */

/* output the current basic solution            */
void outsol (void)
{
    char s[INFOSTRINGLENGTH];
    char smp [2*DIG2DEC(MAX_DIGITS)+4];  
	    /* string to print 2 mp's  into                 */
    int i, row, pos;
    mp num, den;
    
    colset(n+2);    /* column printing to see complementarity of  w  and  z */
    
    colpr("basis=");
    for (i=0; i<=n; i++) 
	{
	if (bascobas[Z(i)]<n)
	    /*  Z(i) is a basic variable        */
	    vartoa(Z(i), s);
	else if (i>0 && bascobas[W(i)]<n)
	    /*  Z(i) is a basic variable        */
	    vartoa(W(i), s);
	else
	    strcpy (s, "  ");
	colpr(s);
	}
	 
    colpr("z=");
    for (i=0; i<=2*n; i++) 
	{
	if ( (row = bascobas[i]) < n)  /*  i  is a basic variable           */
	    {
	    if (i<=n)       /* printing Z(i)        */
		/* value of  Z(i):  scfa[Z(i)]*rhs[row] / (scfa[RHS]*det)   */
		mulint(scfa[Z(i)], A[row][RHS], num);
	    else            /* printing W(i-n)      */
		/* value of  W(i-n)  is  rhs[row] / (scfa[RHS]*det)         */
		copy(num, A[row][RHS]);
	    mulint(det, scfa[RHS], den);
	    reduce(num, den);
	    pos = mptoa(num, smp);
	    if (!one(den))  /* add the denominator  */
		{
		sprintf(&smp[pos], "/");
		mptoa(den, &smp[pos+1]);
		}
	    colpr(smp);
	    }
	else            /* i is nonbasic    */
	    colpr("0");
	if (i==n)       /* new line since printing slack vars  w  next      */
	    {
	    colpr("w=");
	    colpr("");  /* for complementarity in place of W(0)             */
	    }
	}   /* end of  for (i=...)          */
    colout();
}       /* end of outsol                */

/* current basic solution turned into  solz [0..n-1]
 * note that Z(1)..Z(n)  become indices  0..n-1
 * gives a warning if conversion to ordinary rational fails
 * and returns 1, otherwise 0
 */
Bool notokcopysol (void)
{
    Bool notok = 0;
    int i, row;
    mp num, den;
    
    for (i=1; i<=n; i++) 
	if ( (row = bascobas[i]) < n)  /*  i  is a basic variable */
	    {
	    /* value of  Z(i):  scfa[Z(i)]*rhs[row] / (scfa[RHS]*det)   */
	    mulint(scfa[Z(i)], A[row][RHS], num);
	    mulint(det, scfa[RHS], den);
	    reduce(num, den);
            if ( mptoi(num, &(solz[i-1].num), 1) )
                {
                printf("(Numerator of z%d overflown)\n", i);
		notok = 1;
		}
            if ( mptoi(den, &(solz[i-1].den), 1) )
                {
                printf("(Denominator of z%d overflown)\n", i);
		notok = 1;
		}
	    }
	else            /* i is nonbasic    */
	    solz[i-1] = ratfromi(0);
    return notok;
} /* end of copysol                     */

/* --------------- test output and exception routines ---------------- */       
void assertbasic (int v, const char *info)
	/* assert that  v  in VARS is a basic variable         */
	/* otherwise error printing  info  where               */
{
    char s[INFOSTRINGLENGTH];
    if (bascobas[v] >= n)   
	{
	vartoa(v, s);
	fprintf(stderr, "%s: Cobasic variable %s should be basic.\n", info, s);
	errexit("");
	}
}

void assertcobasic (int v, char *info)
	/* assert that  v  in VARS is a cobasic variable       */
	/* otherwise error printing  info  where               */
{
    char s[INFOSTRINGLENGTH];
    if (TABCOL(v) < 0)   
	{
	vartoa(v, s);
	fprintf(stderr, "%s: Basic variable %s should be cobasic.\n", info, s);
	errexit("");
	}
}

void docupivot (leave, enter)
	/* leave, enter in  VARS.  Documents the current pivot. */
	/* Asserts  leave  is basic and  enter is cobasic.      */
{
    char s[INFOSTRINGLENGTH];
    
    assertbasic(leave, "docupivot");
    assertcobasic(enter, "docupivot");
    
    vartoa(leave, s);
    printf("leaving: %-4s ", s);
    vartoa(enter, s);
    printf("entering: %s\n", s);
}       /* end of  docupivot    */

void raytermination (int enter)
{
    char s[INFOSTRINGLENGTH];
    vartoa(enter, s);
    fprintf(stderr, "Ray termination when trying to enter %s\n", s);
    outtabl();
    printf("Current basis, not an LCP solution:\n");
    outsol();
    errexit("");
}

void testtablvars(void)
	/* test tableau variables: error => msg only, continue  */
{
    int i, j;
    for (i=0; i<=2*n; i++)  /* check if somewhere tableauvars wrong */
	if (bascobas[whichvar[i]]!=i || whichvar[bascobas[i]]!=i)
	    /* found an inconsistency, print everything             */
	    {
	    printf("Inconsistent tableau variables:\n");
	    for (j=0; j<=2*n; j++)
		{
		printf("var j:%3d bascobas:%3d whichvar:%3d ",
			j, bascobas[j], whichvar[j]);
		printf(" b[w[j]]==j: %1d  w[b[j]]==j: %1d\n",
			bascobas[whichvar[j]]==j, whichvar[bascobas[j]]==j);
		}
	    break;          
	    }
}
/* --------------- pivoting and related routines -------------- */

/* complement of  v  in VARS, error if  v==Z(0).
 * this is  W(i) for Z(i)  and vice versa, i=1...n
 */
int complement (int v)
{
    if (v==Z(0))
	errexit("Attempt to find complement of z0.");
    return (v > n) ? Z(v-n) : W(v) ;
}       /* end of  complement (v)     */

/* initialize statistics for minimum ratio test
 */
void initstatistics(void)
{
    int i;
    for (i=0; i<=n; i++)
        lextested[i] = lexcomparisons[i] = 0;
}

/* output statistics of minimum ratio test
 */
void outstatistics(void)
{
    int i;
    char s[LCPSTRL];

    colset(n+2);
    colleft(0);
    colpr("lex-column");
    for (i=0; i<=n; i++)
        colipr(i);
    colnl();
    colpr("times tested");
    for (i=0; i<=n; i++)
        colipr(lextested[i]);
    colpr("% times tested");
    if (lextested[0] > 0)
        {
        colpr("100");
        for (i=1; i<=n; i++)
            {
            sprintf(s, "%2.0f",
                    (double) lextested[i] * 100.0 / (double) lextested[0]);
            colpr(s);
            }
        }
    else
        colnl();
    colpr("avg comparisons");
    for (i=0; i<=n; i++)
        if (lextested[i] > 0)
            {
            sprintf(s, "%1.1f",
                (double) lexcomparisons[i] / (double) lextested[i]);
            colpr(s);
            }
        else
            colpr("-");
    colout();
}

/* returns the leaving variable in  VARS, given by lexmin row, 
 * when  enter  in VARS is entering variable
 * only positive entries of entering column tested
 * boolean  *z0leave  indicates back that  z0  can leave the
 * basis, but the lex-minratio test is performed fully,
 * so the returned value might not be the index of  z0
 */
int lexminvar (int enter, int *z0leave)
{                                                       
    int col, i, j, testcol;
    int numcand;
    
    assertcobasic(enter, "Lexminvar");
    col = TABCOL(enter);
    numcand = 0;
	    /* leavecand [0..numcand-1] = candidates (rows) for leaving var */
    /* start with  leavecand = { i | A[i][col] > 0 }                        */
    for (i=0; i<n; i++)
	if (positive (A[i][col]))
	    leavecand[numcand++] = i;
    if (numcand==0) 
	raytermination(enter);
    if (numcand==1)
        {
        lextested[0]      += 1 ;
        lexcomparisons[0] += 1 ;
        *z0leave = (leavecand[0] == bascobas[Z(0)]);
        }
    for (j = 0; numcand > 1; j++)
        /* as long as there is more than one leaving candidate perform
         * a minimum ratio test for the columns of  j  in RHS, W(1),... W(n)
         * in the tableau.  That test has an easy known result if
         * the test column is basic or equal to the entering variable.
         */
	{
        if (j>n)    /* impossible, perturbed RHS should have full rank  */
	    errexit("lex-minratio test failed");
        lextested[j]      += 1 ;
        lexcomparisons[j] += numcand ;

        testcol = (j==0) ? RHS : TABCOL(W(j)) ;
        if (testcol != col)       /* otherwise nothing will change      */
	{
	if (testcol >= 0)
	    /* not a basic testcolumn: perform minimum ratio tests          */
	    {
	    int sgn;
	    int newnum = 0; 
		    /* leavecand[0..newnum]  contains the new candidates    */
	    for (i=1; i < numcand; i++)
		/* investigate remaining candidates                         */
		{
		sgn = comprod(A[leavecand[0]][testcol], A[leavecand[i]][col], 
			      A[leavecand[i]][testcol], A[leavecand[0]][col]);
		/* sign of  A[l_0,t] / A[l_0,col] - A[l_i,t] / A[l_i,col]   */
		/* note only positive entries of entering column considered */
		if (sgn==0)         /* new ratio is the same as before      */
		    leavecand[++newnum] = leavecand[i];
		else if (sgn==1)    /* new smaller ratio detected           */
		    leavecand[newnum=0] = leavecand[i];
		}
	    numcand = newnum+1;
	    }
	else
            /* testcol < 0: W(j) basic, Eliminate its row from leavecand    */
	    /* since testcol is the  jth  unit column                       */
	    for (i=0; i < numcand; i++)
		if (leavecand[i] == bascobas[W(j)])
		    {
		    leavecand[i] = leavecand[--numcand];
			    /* shuffling of leavecand allowed       */
		    break;
		    }
	}   /* end of  if(testcol != col)                           */
   
	if (j==0)
	    /* seek  z0  among the first-col leaving candidates     */
	    for (i=0; i<numcand; i++)
		if ( (*z0leave = (leavecand[i] == bascobas[Z(0)])) )
		    break;
		    /* alternative, to force z0 leaving the basis:
		     * return whichvar[leavecand[i]];
		     */
        }       /* end of  for ( ... numcand > 1 ... )   */
    return whichvar[leavecand[0]];
}       /* end of lexminvar (col, *z0leave);                        */


/* returns the leaving variable in  VARS  as entered by user,
 * when  enter  in VARS is entering variable
 * only nonzero entries of entering column admitted
 * boolean  *z0leave  indicates back that  z0  has been
 * entered as leaving variable, and then
 * the returned value is the index of  z0
 */
int interactivevar (int enter, int *z0leave)
{                                                       
    char s[INFOSTRINGLENGTH], instring[2];

    int inp, col, var;
    int breject = 1;
    assertcobasic(enter, "interactivevar");
    col = TABCOL(enter);

    vartoa(enter, s);
    printf("   Entering variable (column): %s\n", s);
    while (breject)
	{
	printf("   Leaving row (basic variable z.. or w..), ");
	printf("or 't' for tableau:\n");
	strcpy(instring, "?");
	if (scanf("%1s", instring)==EOF)
	    {
	    printf ("Input terminated too early with EOF\n");
	    exit(1);
	    }
	if ( instring[0] == 't')
	    {
	    printf("\n");
	    outtabl();
	    vartoa(enter, s);
	    printf("   Entering variable (column): %s\n", s);
	    continue;
	    }
	scanf("%d", &inp);
	printf("   You typed %s%d\n", instring, inp);
	if ( (inp < 0) || (inp > n))
	    {
	    printf("Variable index %d outside 0..n=%d\n",
		    inp, n);
	    continue;
	    }
	if ( instring[0] == 'w')
	    {
	    if (inp == 0)
		{
		printf("Variable w0 not allowed\n");
		continue;
		}
	    var = inp + n;
	    }
	else if ( instring[0] == 'z')
	    var = inp;
	else 
	    {
	    printf("Variable not starting with  z  or  w\n");
	    continue;
	    }
	/* var == variable in VARS giving what has been input   */
	if ( bascobas[var] >= n)
	    {
	    vartoa (var, s);
	    printf("Variable %s not basic\n", s);
	    continue;
	    }
	if ( zero( A [bascobas[var]] [col] ) )
	    {
	    vartoa (var, s);
	    printf("Row %s has zero pivot element, not allowed\n", s);
	    continue;
	    }
	breject = 0;    /* now everything ok            */
	}       /* end of  while (breject) for input    */
*z0leave = (var == Z(0));
return var;
}   /* end of  interactivevar (col, *z0leave);          */

void negcol(int col)
	/* negate tableau column  col   */
{
    int i;
    for (i=0; i<n; i++)
	changesign(A[i][col]);
}

void negrow(int row)
	/* negate tableau row.  Used in  pivot()        */
{
    int j;
    for (j=0; j<=n+1; j++)
	if (!zero(A[row][j]))
	    changesign(A[row][j]);
}

/* leave, enter in  VARS  defining  row, col  of  A
 * pivot tableau on the element  A[row][col] which must be nonzero
 * afterwards tableau normalized with positive determinant
 * and updated tableau variables
 */
void pivot (int leave, int enter)
{
    int row, col, i, j;
    int nonzero, negpiv;
    mp pivelt, tmp1, tmp2;
    
    row = bascobas[leave];
    col = TABCOL(enter);
    
    copy (pivelt, A[row][col]);     /* pivelt anyhow later new determinant  */
    negpiv = negative (pivelt);
    if (negpiv)
	changesign(pivelt);
    for (i=0; i<n; i++)
	if (i != row)               /*  A[row][..]  remains unchanged       */
	    {
	    nonzero = !zero(A[i][col]);
	    for (j=0; j<=n+1; j++)      /*  assume here RHS==n+1        */
		if (j != col)
		    /*  A[i,j] =
		       (A[i,j] A[row,col] - A[i,col] A[row,j])/ det     */
		    {
		    mulint (A[i][j], pivelt, tmp1);
		    if (nonzero)
			{
			mulint(A[i][col], A[row][j], tmp2);
			linint(tmp1, 1, tmp2, negpiv ? 1 : -1);
			}
		    divint (tmp1, det, A[i][j]);
		    }
	    /* row  i  has been dealt with, update  A[i][col]  safely   */
	    if (nonzero && !negpiv)
		changesign (A[i][col]);
	    }       /* end of  for (i=...)                              */
    copy(A[row][col], det);
    if (negpiv)
	negrow(row);
    copy(det, pivelt);      /* by construction always positive      */
    
    /* update tableau variables                                     */
    bascobas[leave] = col+n;        whichvar[col+n] = leave;
    bascobas[enter] = row;          whichvar[row]   = enter;
}       /* end of  pivot (leave, enter)                         */

/* ------------------------------------------------------------ */ 
void runlemke(Flagsrunlemke flags)
{
    int leave, enter, z0leave;

    pivotcount = 1;
    initstatistics();

    isqdok();
    /*  printf("LCP seems OK.\n");      */

    filltableau();
    /*  printf("Tableau filled.\n");    */

    if (flags.binitabl)
	{
	printf("After filltableau:\n");
	outtabl();
	}
    
    /* z0 enters the basis to obtain lex-feasible solution      */
    enter = Z(0);                                                   
    leave = flags.binteract ? interactivevar(enter, &z0leave) :
	    lexminvar(enter, &z0leave) ;
    
    /* now give the entering q-col its correct sign             */
    negcol (RHS);   
    
    if (flags.bouttabl) 
	{
	printf("After negcol:\n");
	outtabl();
	}
    while (1)       /* main loop of complementary pivoting                  */
	{
	testtablvars();
	if (flags.bdocupivot)
	    docupivot (leave, enter);
	pivot (leave, enter);
	if (z0leave)
	    break;  /* z0 will have value 0 but may still be basic. Amend?  */
	if (flags.bouttabl) 
	    outtabl();
	enter = complement(leave);
	leave = flags.binteract ? interactivevar(enter, &z0leave) :
		lexminvar(enter, &z0leave) ;
        if (pivotcount++ == flags.maxcount)
	    {
            printf("------- stop after %d pivoting steps --------\n", 
		   flags.maxcount);
	    break;
	    }
	}
    
    if (flags.binitabl)
	{
	printf("Final tableau:\n");
	outtabl();
	}
    if (flags.boutsol)
	outsol();
    if (flags.blexstats)
        outstatistics();
    
    notokcopysol();
}




