/* col.h                                                        */
/* automatic pretty printing in columns                         */
/* 17 Apr 2000                                                  */
/* Bernhard von Stengel  stengel@maths.lse.ac.uk      */

/* no of bytes of buffer to print into.  Buffer will be 
 * printed and flushed if full (not assumed to happen),
 * suffices for one page of text  
 */
#define COLBUFSIZE 30000  
#define ISTR 20 	/* no of bytes to print an integer      */

/* resetting buffer with  c  columns
 * is assumed to be called before all other routines
 */
void colset(int c);

/* print integer i  into the current column                     */
void colipr(int i);

/* making column  c  in  0..ncols-1  left-adjusted              */
void colleft(int c);

/* terminate current line early.  blank line if in col  0       */
void colnl(void);

/* print out the current buffer, without flushing               */
void colout(void);

/* store string  s  into the current column,
 * updating column width
 */
void colpr(const char *s);
