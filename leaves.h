/* leaves.h
 * 1 June 2000
 */

/* #include before: "rat.h" */

#define MAXLINELENGTH 1000	/* max chars in input line	*/

/* read in lines from stdin, process them one by one by
 * oklinetobehavprob()
 * output leaves in equilibrium path with probabilities
 */ 
void leavesfrominput (void);

/* convert a line previously generated with  showeq()  
 * where bshortequil == 1  to behavior strategies in ->behavprob
 * for the two players.
 * return 1 if the line was syntactically OK, otherwise 0
 */ 
Bool oklinetobehavprob (char *line);

/* output the leaves that are reached by  ->realprob in the form
 * leafno:jointprob 
 * prefixed by  "SUP>(docuseed)<"  to identify the game
 */ 
void outleavesreached (int docuseed);


