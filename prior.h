/* prior.h
 * 27 Apr 2000
 * prior generation
 */

/* max. allowed denominator for random fractions */
#define MAXACCURACY 1000
#define DEFAULTACCURACY 23	/* max. denominator for random fractions */
#define FIRSTPRIORSEED  500 	/* first seed for prior generation	 */

/* flags for  genprior 	*/
typedef struct
    {
    int   seed 	; 	/* 0: centroid,
                         * >0:  random seed for prior, will be added to
			 * FIRSTPRIORSEED
			 */
    int   accuracy; 	/* largest denominator for random prior,
			 * possibly smaller when only two probabilities
			 * via continued fractions generation.
    			 * default DEFAULTACCURACY
			 */
    }
    Flagsprior; 

/* generate prior, stored in  moves[]->behavprob
 */
void genprior(Flagsprior flags);

/* output prior as pair of behavior strategies
 * realplan[][]  must be allocated 
 */
void outprior(void);
