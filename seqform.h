/* seqform.h
 * 27 Apr 2000
 */

/* #include before: "treedef.h"         */

/* global variables for sequence form   */

/* sf payoffs [row][col]                */
extern  Payvec **sfpay;

/* constraint matrix [player][row][col]
 * sacrifice one unused pointer for player 0
 */
extern  int **sfconstr[PLAYERS];

/* allocates sequence form payoff and constraint matrices
 * sets  all  sf  payoffs to 0
 * allocate  realplan[0..PLAYERS-1]
 * assumes  nseqs[], nisets[]  are set via  genseqin()
 * frees old data and keeps track of old dimensions
 */
void allocsf(void);

/* allocate & generate sequence form:  payoff and constraint matrices
 * h->seqin must be defined via  genseqin()
 */
void gensf(void);

/* LCP  (M,q) for the sequence form
 * not the covering vector  d
 * generates  sf,  allocates  LCP,  fills LCP
 */
void sflcp(void);

/* copy for player  pl  moves[]->realprob  to  rplan,
 * usually  rplan = realplan[pl],    assumed allocated
 */
void realplanfromprob(int pl, Rat *rplan);

/* asserts that for player  pl  moves[]->realprob == probvector,
 * if (bcomplain):  prints to stdout all differing positions
 */
Bool iseqrealplantoprob(int pl, Rat *rplan, Bool bcomplain);

/* how many isets of player  pl  have nondeterministic moves
 * in the realization plan  rplan  (e.g. if 0: pure strategy)
 */
int  propermixisets(int pl, Rat *rplan);

/* gives the realization plan  rplan,  typically  realplan[pl],
 * of player  pl,  to stdout, in two lines
 */
void outrealplan(int pl, Rat *rplan);

/* gives the realization plan  rplan  of player  pl
 * as sparse behaviour strategy to stdout, in one line
 * bnewline: terminate with \n  
 */
void outbehavstrat(int pl, Rat *rplan, Bool bnewline);

/*  output sequence form        */
void sfprint(void);
