/* gambit.h
 * interfaces to the GAMBIT program
 * 14 June 2000
 */

/* include before: <stdio.h>, "rat.h"	*/

/* open file for GAMBIT interfacing with that name, abort if fails	
 * print opening brace into that file
 * set internal Boolean variable to remember comma separation
 */
void gambopenfile(char *filename) ;

/* print closing brace into file for GAMBIT interfacing 
 * close file
 */
void gambclosefile(void) ;

/* print to gambit interface file
 * equilibrium computed in  solz[]  by runlemke()
 * knows by itself if to print comma before that or not
 */
void gambshoweq(void) ;

/* gives the realization plan  rplan  of player  pl
 * as behaviour strategy as list (= behavior strategy) of lists 
 * (= probability distributions) of rationals  with braces and commas
 * to the file pointed to by  fp
 */
void glistbehavstrat(FILE *fp, int pl, Rat *rplan);

