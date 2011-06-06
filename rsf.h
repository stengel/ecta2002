/* rsf.h
 * reduced sequence form
 * 7 Apr 2000
 */

/* #include before: "rat.h"                                     */

/* global vars for reduced sequence form, see  genredsf(pl)     */
extern  int irreddim[PLAYERS];
extern  int redsfdim[PLAYERS];
extern  int ** redsfconstr[PLAYERS];
extern  int *  redsfrhs[PLAYERS];
extern  int ** realplfromredsf[PLAYERS];
extern  int *  realplconst[PLAYERS];


/* reduced sequence form generation for player  pl :
 * - matrix of reduced sequence form constraints
 *             redsfconstr[pl] [ irreddim[pl] ] [ redsfdim[pl] ]
 *   with RHS     redsfrhs[pl] [ irreddim[pl] ]
 * - matrix to recover realization plans
 *                  realplfromredsf[pl] [ nseqs[pl] ] [ redsfdim[pl] ]
 *   constant added:    realplconst[pl] [ nseqs[pl] ] 
 * - both matrices allocated here (and old data freed before)
 * - dimensions also generated here:
 *   irreddim[pl] : number of irredundant variables
 *   redsfdim[pl] : no. of free vars, used in LCP
 */
void genredsf(int pl);

/* set reduced sequence form complementarity constraints
 * with payoff matrix of player  paytopl  (1 or 2)
 * into the matrix               intoM  
 * starting at offset            intoM [rowoffset] [coloffset]
 * payoff matrix negated if      negpaymatrix
 * and right hand side           rhsvec  (no offset used)
 * rhs is  ahat, -fhat  unless   negrhs
 * settings for LCP:
 * negpaymatrix = negrhs = 1:    intoM = lcpM,  rhsvec = rhsq (+ offsets)
 */
void setrsfcomplconstr (int paytopl, Rat ** intoM,
	int rowoffset, int coloffset, Bool negpaymatrix,
	Rat * rhsvec, Bool negrhs);

/* LCP  (M,q) for the reduced sequence form
 * not the covering vector  d
 * assumes  sf  has been generated
 * generates  rsf,  allocates  LCP,  fills LCP
 */
void rsflcp(void);

