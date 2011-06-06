/* treegen.c
 * 1 May 2001 
 * generating examples of game trees
 */

#include <stdio.h>
	/* NULL	*/
#include <stdlib.h>
        /* srand()      */
#include "rat.h"
#include "treedef.h"

#include "treegen.h"

/* ---------- auxiliary routine for  createbintree ------------ */
static int pow2(int n)   
        /* given  n>=0, returns 2^n     */
{
    int a;
    for (a=1; n>0; n--)
        a *= 2;
    return a;
}       /* end of  pow2(n)              */

void createbintree (int levels, int seed)
{
    int k, l;
    Iset    curriset;
    Move    currmove;
    Outcome currout;
    int nnodes = pow2(levels+1);
    
    alloctree (nnodes, pow2(levels-1), pow2(levels)+3, pow2(levels));
    srand(PAYOFFSEEDPLUS + seed);

    /* creating tree fathers, root has no father                    */
    root = nodes + ROOT;
    root->father = NULL;
    for (k=ROOT+1; k < nnodes; k++)
        nodes[k].father = nodes + k/2;
    
    /* defining isets, moves of player 0,1    */
    firstiset[0] = firstiset[1] = curriset = isets;
    firstmove[0] = moves;
    firstmove[1] = currmove = moves + 1;
    currmove++;             /* leave firstmove[1] for empty sequence */
    for (l=0; l < levels; l+= 2)
        /* the even levels belong to player 1       */
        {
        int lastnodeonlevel = pow2(l+1);
        for ( k=pow2(l); k < lastnodeonlevel; k++)
            {
            nodes[k].iset   = curriset;
            if (k % 2)      /* k is odd     */
                {
                curriset->player = 1;
                curriset->nmoves = 2;
                curriset->move0  = currmove;
                currmove->atiset = curriset;
                currmove ++;
                currmove->atiset = curriset;
                currmove ++;
                curriset ++;
                }
            }
        }
    /* defining isets, moves of player 2      */
    firstiset[2] = curriset;
    firstmove[2] = currmove;
    currmove++;             /* leave firstmove[2] for empty sequence */
    for (l=1; l < levels; l+= 2)
        /* the odd levels belong to player 2        */
        {
        int lastnodeonlevel = pow2(l+1);
        for ( k=pow2(l); k < lastnodeonlevel; k++)
            {
            nodes[k].iset   = curriset;
            if (k % 2)      /* k is odd     */
                {
                curriset->player = 2;
                curriset->nmoves = 2;
                curriset->move0  = currmove;
                currmove->atiset = curriset;
                currmove ++;
                currmove->atiset = curriset;
                currmove ++;
                curriset ++;
                }
            }
        }
    /* consistency check    */
    if (curriset != firstiset[3])
        printf("wrong number %d of info sets!!!\n", (int) (curriset-isets));
    if (currmove != firstmove[3])
        printf("wrong number %d of moves!!!\n", (int) (currmove-moves));
    
    /* creating  reachedby  move except for root    */
    /* even-numbered nodes reached by  move0, odd-numbered by  move0+1      */
    for (k=ROOT+1; k < nnodes; k++)
        nodes[k].reachedby = nodes[k].father->iset->move0 + k % 2 ;
    
    /* outcomes    */
    l = pow2(levels);
    for (k=ROOT; k < l; k++)
        nodes[k].terminal = 0;
    currout = outcomes;
    for (k=pow2(levels); k < nnodes; k++)
        {
        nodes[k].terminal = 1;
        nodes[k].outcome  = currout;
        /* generating  contents of  currout , random payoffs    */
        currout->whichnode = nodes + k;
        currout->pay[0] = ratfromi(nextrandpay());
        currout->pay[1] = ratfromi(nextrandpay());
        currout++;
        }
    /* consistency check    */
    if (currout != lastoutcome)
        printf("wrong number %d of outcomes!!!\n", currout-outcomes );
}       /* end of createbintree(levels, seed)   */


void tracingexample(void)
{
    int pay[2][8] = { {11, 3, 0,  0, 0, 24, 6, 0},
                      { 3, 0, 0, 10, 4,  0, 0, 1} };
    int i;
    Outcome z;
    alloctree(16, 5, 13, 8);
    firstiset[0] = isets;
    firstiset[1] = isets + 1;
    firstiset[2] = isets + 3;
    firstmove[0] = moves;
    firstmove[1] = moves + 3;
    firstmove[2] = moves + 8;
    
    root = nodes + ROOT;
    root->father = NULL;
    nodes[2].father = root;
    nodes[3].father = nodes + 2;
    nodes[4].father = root;
    nodes[5].father = nodes + 3;
    nodes[6].father = nodes + 3;
    nodes[7].father = nodes + 2;
    z = outcomes;
    for (i=8; i<=15; i++)
        {
        nodes[i].father = nodes + i/2;
        nodes[i].terminal = 1;
        nodes[i].outcome  = z;
        z->whichnode = nodes + i;
        z->pay[0] = ratfromi(pay[0][i-8]);
        z->pay[1] = ratfromi(pay[1][i-8]);
        z++;
        }
    nodes[1].iset = firstiset[1];
    nodes[2].iset = firstiset[1]+1;
    nodes[3].iset = firstiset[0];
    nodes[4].iset = firstiset[2];
    nodes[5].iset = firstiset[2];
    nodes[6].iset = firstiset[2]+1;
    nodes[7].iset = firstiset[2]+1;
    
    nodes[2].reachedby = firstmove[1]+2;    /* note empty sequence  */
    nodes[3].reachedby = firstmove[1]+3;
    nodes[4].reachedby = firstmove[1]+1;
    nodes[5].reachedby = firstmove[0]+1;
    nodes[6].reachedby = firstmove[0]+2;
    nodes[7].reachedby = firstmove[1]+4;
    nodes[8].reachedby  = firstmove[2]+1;
    nodes[9].reachedby  = firstmove[2]+2;
    nodes[10].reachedby = firstmove[2]+1;
    nodes[11].reachedby = firstmove[2]+2;
    nodes[12].reachedby = firstmove[2]+3;
    nodes[13].reachedby = firstmove[2]+4;
    nodes[14].reachedby = firstmove[2]+3;
    nodes[15].reachedby = firstmove[2]+4;
    
    isets[0].player = 0;
    isets[1].player = 1;
    isets[2].player = 1;
    isets[3].player = 2;
    isets[4].player = 2;
    
    isets[0].move0 = firstmove[0]+1;
    isets[1].move0 = firstmove[1]+1;
    isets[2].move0 = firstmove[1]+3;
    isets[3].move0 = firstmove[2]+1;
    isets[4].move0 = firstmove[2]+3;
    
    isets[0].nmoves = 2;
    isets[1].nmoves = 2;
    isets[2].nmoves = 2;
    isets[3].nmoves = 2;
    isets[4].nmoves = 2;
    
    moves[1].atiset = firstiset[0];
    moves[2].atiset = firstiset[0];
    moves[4].atiset = firstiset[1];
    moves[5].atiset = firstiset[1];
    moves[6].atiset = firstiset[1]+1;
    moves[7].atiset = firstiset[1]+1;
    moves[9].atiset  = firstiset[2];
    moves[10].atiset = firstiset[2];
    moves[11].atiset = firstiset[2]+1;
    moves[12].atiset = firstiset[2]+1;
    
    moves[1].behavprob.num = 1;
    moves[1].behavprob.den = 2;
    moves[2].behavprob.num = 1;
    moves[2].behavprob.den = 2;
}       /* end of  tracingexample()     */


void forwardexample(void)
{
    int pay[2][8] = { {11, 3, 0,  0, 0, 24, 6, 0},
                      { 3, 0, 0, 10, 4,  0, 0, 1} };
    int i;
    Outcome z;
    alloctree(16, 5, 13, 8);
    firstiset[0] = isets;
    firstiset[1] = isets + 1;
    firstiset[2] = isets + 3;
    firstmove[0] = moves;
    firstmove[1] = moves + 3;
    firstmove[2] = moves + 8;
    
    root = nodes + ROOT;
    root->father = NULL;
    nodes[2].father = root;
    nodes[3].father = nodes + 2;
    nodes[4].father = root;
    nodes[5].father = nodes + 3;
    nodes[6].father = nodes + 3;
    nodes[7].father = nodes + 2;
    z = outcomes;
    for (i=8; i<=15; i++)
        {
        nodes[i].father = nodes + i/2;
        nodes[i].terminal = 1;
        nodes[i].outcome  = z;
        z->whichnode = nodes + i;
        z->pay[0] = ratfromi(pay[0][i-8]);
        z->pay[1] = ratfromi(pay[1][i-8]);
        z++;
        }
    nodes[1].iset = firstiset[1];
    nodes[2].iset = firstiset[1]+1;
    nodes[3].iset = firstiset[0];
    nodes[4].iset = firstiset[2];
    nodes[5].iset = firstiset[2];
    nodes[6].iset = firstiset[2]+1;
    nodes[7].iset = firstiset[2]+1;
    
    nodes[2].reachedby = firstmove[1]+2;    /* note empty sequence  */
    nodes[3].reachedby = firstmove[1]+3;
    nodes[4].reachedby = firstmove[1]+1;
    nodes[5].reachedby = firstmove[0]+1;
    nodes[6].reachedby = firstmove[0]+2;
    nodes[7].reachedby = firstmove[1]+4;
    nodes[8].reachedby  = firstmove[2]+1;
    nodes[9].reachedby  = firstmove[2]+2;
    nodes[10].reachedby = firstmove[2]+1;
    nodes[11].reachedby = firstmove[2]+2;
    nodes[12].reachedby = firstmove[2]+3;
    nodes[13].reachedby = firstmove[2]+4;
    nodes[14].reachedby = firstmove[2]+3;
    nodes[15].reachedby = firstmove[2]+4;
    
    isets[0].player = 0;
    isets[1].player = 1;
    isets[2].player = 1;
    isets[3].player = 2;
    isets[4].player = 2;
    
    isets[0].move0 = firstmove[0]+1;
    isets[1].move0 = firstmove[1]+1;
    isets[2].move0 = firstmove[1]+3;
    isets[3].move0 = firstmove[2]+1;
    isets[4].move0 = firstmove[2]+3;
    
    isets[0].nmoves = 2;
    isets[1].nmoves = 2;
    isets[2].nmoves = 2;
    isets[3].nmoves = 2;
    isets[4].nmoves = 2;
    
    moves[1].atiset = firstiset[0];
    moves[2].atiset = firstiset[0];
    moves[4].atiset = firstiset[1];
    moves[5].atiset = firstiset[1];
    moves[6].atiset = firstiset[1]+1;
    moves[7].atiset = firstiset[1]+1;
    moves[9].atiset  = firstiset[2];
    moves[10].atiset = firstiset[2];
    moves[11].atiset = firstiset[2]+1;
    moves[12].atiset = firstiset[2]+1;
    
    moves[1].behavprob.num = 1;
    moves[1].behavprob.den = 2;
    moves[2].behavprob.num = 1;
    moves[2].behavprob.den = 2;
}       /* end of  forwardexample()     */

