/* SCCS @(#)bsplit.c	1.6 06/06/01 */
/*
** The routine which will find the best split for a node
**
** Input :      node
**              node number
**
** Output:      Fills in the node's
**                      primary splits
**                      competitor splits
**
** ALG: adjusted to compute penalty, best split and improve in
** three separate steps.  Also added rp_parent_objective call
** at the beginning to get scaling factor to trade the improve and
** penalty numbers.
*/
#include "rpart.h"
#include "node.h"
#include <stdio.h>
#include "rpartproto.h"

void bsplit(struct node *me, int nodenum, int variables_used[])
{
    int i, j, k;
    int nc;
    double improve;
    FLOAT split = 0.0; /* in case choose does not set it */
    struct split *tsplit;
    Sint *index;
    int  *which;
    FLOAT *xtemp;  /*these 3 because I got tired of typeing "rp.xtemp", etc*/
    double **ytemp;
    double *wtemp;  
    double penalty;  //added ALG 4/1//2012
    double parent_objective; //ALG 4/30/2012. this is what we scale improve by when looking at penalties

    which = rp.which;
    xtemp = rp.xtemp;
    ytemp = rp.ytemp;
    wtemp = rp.wtemp;

    /*
    ** test out the variables 1 at at time
    */
    me->primary =0;
    for (i=0; i<rp.nvar; i++) /* loop through each predictor variable */
    {
    	//if(i==0){
    	//	Rprintf("---------------CONSIDERING NEW NODE------------------\n");
    	//	Rprintf("%d\n",nodenum);
    	//}

    	index = rp.sorts[i];
    	nc = rp.numcat[i];

    	/* extract x and y data */
    	k=0;
    	for (j=0; j<rp.n; j++){
    		if (index[j] >=0 && which[index[j]]== nodenum) {
    			xtemp[k] = rp.xdata[i][j];
    			ytemp[k] = rp.ydata[index[j]];
    			wtemp[k] = rp.wt[index[j]];
    			k++;
			}
    	} /* end of extracting x & y data loop */



    	if (k==0 || (nc==0 &&  xtemp[0]==xtemp[k-1])) continue;  /*stop if no place to split */

    	/*
    	 * ALG 5/1/2012.
    	 * Compute parent_objective.
    	 */
    	if(i==0){
    		(*rp_parent_objective)(k, ytemp, rp.dummy, &parent_objective, wtemp);
    		//Rprintf("set parent val \n");
    	}
    	//Rprintf("%lf\n",parent_objective);

    	//ALG 4/1/2012: compute penalty for splitting this node on variable i
    	//given variables_used
    	penalty = (*rp_penalty)(i, variables_used);

    	/*ALG 1/21/2012: added last two arguments to fit expanded choose */
    	/*ALG 1/21/2012: added third to last argument */
    	/*ALG 4/11/2012: removed 1/21 args, replaced with pre-computed penalty*/
    	(*rp_choose)(k, ytemp, xtemp, nc, rp.min_node, &improve,
			     &split, rp.csplit, me->risk, wtemp);
    	//Rprintf("%lf\n",improve);

    	//ALG:
    	//combine 'improve', the best objective function value which is some form of
    	//improve = parent_objective - f(left & right objectives)
	//with the 'parent_objective', which is what we normalize improve by.
    	//if the program is set to 'scale by root', we end up ignoring the parent objective value
	// in the function call below.
    	improve = (*rp_improve)(parent_objective,improve,penalty);
    	//Rprintf("%lf\n",improve);

		/*
		** Originally, this just said "if (improve >0)", but rounding
		**  error will sometimes create a zero that's not 0.  Yet we
		**  want to retain invariance to the scale of "improve".
		*/
		if (improve > rp.iscale) {  //adjusted for penalized_improve
			rp.iscale = improve;  /*largest seen so far*/
		}

		if (improve > (rp.iscale * 1e-10)) {  /* improvement is big enough */

			improve /= rp.vcost[i];   /* scale the improvement */
			tsplit = insert_split(&(me->primary), nc, improve, rp.maxpri); //ALG 2/2/2012: decides if this is best so far, if so returns split obj

			if (tsplit !=0) {  /* best split seen so far, so set variables*/
				tsplit->improve = improve;
				tsplit->var_num = i;  /* the variable we split on */
				tsplit->spoint  = split;
				tsplit->count   = k;

				if (nc ==0) {
					tsplit->spoint = split;
					tsplit->csplit[0]= rp.csplit[0];
				}
				else for (k=0; k<nc; k++) tsplit->csplit[k] = rp.csplit[k];
			}
		}/* end of 'if' for big enough improvement */
	}/* end of loop through each predictor variable */
}
