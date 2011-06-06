/* lemke.h
 * 16 Apr 2000
 * declarations for lcp solver
 */

/* #include before:  rat.h       */

#define MAXLCPDIM 2000       /* max LCP dimension                       */
#define INFOSTRINGLENGTH 8   /* string naming vars, e.g. "z0", "w187"   */
#define LCPSTRL  60          /* length of string containing LCP entry   */

/* LCP input data                                               */
extern  Rat **lcpM;             /* LCP Matrix                   */
extern  Rat *rhsq;              /* right hand side  q           */
extern  Rat *vecd;              /* LCP covering vector  d       */
extern  int lcpdim;             /* LCP dimension                */

/* LCP result data                                              */
extern  Rat  *solz;             /* LCP solution  z  vector      */
/* no. of Lemke pivot iterations, including the first to pivot z0 in    */
extern  int  pivotcount;

/* allocate and initialize with zero entries an LCP of dimension  n
 * this is the only method changing  lcpdim  
 * exit with error if fails, e.g. if  n  not sensible
 */
void setlcp(int n);

/* output the LCP as given      */
void outlcp (void);

/* flags for  runlemke  */
typedef struct
    {
    int   maxcount  ;   /* max no. of iterations, infinity if 0         */
    int   bdocupivot;   /* Y/N  document pivot step                     */
    int   binitabl ;    /* Y/N  output entire tableau at beginning/end  */
    int   bouttabl  ;   /* Y/N  output entire tableau at each step      */
    int   boutsol   ;   /* Y/N  output solution                         */
    int   binteract ;   /* Y/N  interactive pivoting                    */
    int   blexstats ;   /* Y/N  statistics on lexminratio tests         */
    }
    Flagsrunlemke;

/* solve LCP via Lemke's algorithm,
 * solution in  solz [0..lcpdim-1]
 * exit with error if ray termination
 */
void runlemke(Flagsrunlemke flags);
