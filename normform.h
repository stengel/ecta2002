/* normform.h
 * 27 Apr 2000
 */

/* #include before: "rat.h" "treedef.h"         */

/* global variables for normal form             */
/* nf payoffs [row][col]                        */
extern  Payvec **nfpay;
/* number of pure strategies for each player
 * set by  gennf() via numstratsnfpre(pl)
 */
extern  int  nstrats[PLAYERS]; 
extern  Move *movetuple[PLAYERS];   /* movetuple encoding pure strategy */   
extern  Rat  *mixedstrat[PLAYERS];  /* mixed strategy (vector of probs) */ 

/* allocates normal form payoff matrix
 * sets  all  nf  payoffs to 0
 * allocates movetuple[][]
 * allocates mixedstrat[][]
 * assumes  nstrats[], nisets[]  set
 */
void allocnf(void);

/* computes number of RNF pure strategies for player  pl
 * sets  h->ncontin, h->prefact, c->ncompat, c->offset
 */
int numstratsnfpre(int pl);

/* fills  list  of pure strategies compatible with  seq
 * of player  pl,  returns length of list
 * assumes  list  is allocated size  nstrats[pl]
 */
int seqtostratlist (Move seq, int pl, int *list);

/* generate normal form
 * sets  nstrats[1..PLAYERS-1]
 */
void gennf(void);

/* LCP  (M,q) for the normal form
 * not the covering vector  d
 * generates  nf,  allocates  LCP,  fills LCP
 */
void nflcp(void);

/* converts pure strategy number  strat  of player  pl
 * to its tuple of moves (NULL = irrelevant)
 * in  movetuple[nisets[pl]], assumed  allocated
 * assumes  numstratsnfpre(pl)  has been called
 */
void strattomovetuple (int strat, int pl);

/* converts the behavior strategy of player  pl  given by
 * moves[]->behavprob  to a mixed strategy stored in  mixedstrat[pl]
 */
void behavtomixed(int pl);

/* converts the mixed strategy  mixed  of player  pl,  e.g.
 * mixed = mixedstrat[pl]  to a realization plan   rplan,
 * typically  rplan = realplan[pl]
 */
void mixedtorealplan(int pl, Rat *mixed, Rat *rplan);

/* converts  movetuple[pl]  to string  s
 * s  must be long enough to contain result
 * returns length of string
 */
int  movetupletoa (int pl, char* s);

/* support size of the mixed strategy  mixed,  typically  mixedstrat[pl],
 * of player  pl
 */
int  supportsize(int pl, Rat *mixed);

/* gives the mixed strategy  mixed of player  pl to stdout,
 * in one line, listing only the support
 * bnewline: terminate with \n  
 */
void outmixed(int pl, Rat *mixed, Bool bnewline);

/* gives the mixed strategy  mixed,  typically  mixedstrat[pl],
 * of player  pl,  to stdout, in two lines, listing only the support
 */
void twolineoutmixed(int pl, Rat *mixed);

/*  output normal form          */
void nfprint(void);

/* columns: sequences, row: pure strategies of player  pl       */
void compatstrats(int pl);
