/* alloc.h
 * 17 Apr 2000
 * memory allocation with error output if fails
 */

#define CALLOC(n,s) xcalloc((size_t) n, (size_t) s,__LINE__,__FILE__)

/* Example: nodes = TALLOC(nn, struct node);
 * necessary type:   struct node *nodes;
 */
#define TALLOC(n,type) (type *) xcalloc((size_t) (n), sizeof(type),__LINE__,__FILE__)

/* T2ALLOC :
 * allocate a 2-dim (nrows,ncols) array of  type  to ptr
 * Example:  T2ALLOC (sfpay, nseqs[1], nseqs[2], Payvec);
 * necessary type:    Payvec **sfpay;
 */

#define T2ALLOC(ptr,nrows,ncols,type) {int i; \
    ptr = TALLOC ((nrows), type *);             \
    for (i=0; i < (nrows); i++) ptr[i] = TALLOC ((ncols), type);}

/* FREE2: free a 2-dim (nrows) array, usually allocated with T2ALLOC */

#define FREE2(ptr,nrows) {int i; \
    for (i=0; i < nrows; i++) free(ptr[i]); free(ptr);}

/* allocate  n  objects of size  s,
 * if fails: error noting line  l  and filename  f
 */
void * xcalloc(size_t n, size_t s, int l, char* f);

