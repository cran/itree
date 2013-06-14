/* SCCS  @(#)partition.c	1.5 01/06/00 */
/*
** The main workhorse of the recursive partitioning module.  When called
**   with a node, it partitions it and then calls itself to partition the
**   children it has created.
** If the node is not splittable (too few people, or complexity is too small)
**   it simply returns.  The routine may not be able to discover that the
**   complexity is too small until after the children have been partitioned,
**   so it needs to check this at the end.
**
** 1/16/2012  ALG: keep track of which variables you've split on.
**
** 4/28/2012 added rp.split_check_offset to ensure that for some methods,
** we do not pay attention to the complexity calculations.
*/
#include <stdio.h>
#include "rpart.h"
#include "node.h"
#include "rpartproto.h"
#include "rpartS.h"
#include <math.h>  //for infinity

int partition(int nodenum, struct node *splitnode, double *sumrisk, int variables_used[]){

    struct node *me;
    double tempcp;
    int i,j;
    double tempcp2;
    double left_risk, right_risk;
    int left_split, right_split;
    double twt;
    int primary_split_var;  /*ALG 1/26/2012: primary split's variable number for this node */
    int vars_used_update[rp.max_depth]; //1/31/2012

    me = splitnode;
    if (nodenum >1) {  //NOT root node, so we need to call the eval fcn
    	j=0; /* number of observations currently in the node */
    	twt =0;  /* weight on the observations currently in the node */
    	for (i=0; i<rp.n; i++){
    		if (rp.which[i] == nodenum) { /*is observation i in this node */
    			rp.wtemp[j]   = rp.wt[i];
    			rp.ytemp[j++] = rp.ydata[i];
    			twt += rp.wt[i];
    		}
    	}

    	(*rp_eval)(j, rp.ytemp, me->response_est, &(me->risk), rp.wtemp);
    	me -> num_obs = j;
    	me -> sum_wt  = twt;
    	tempcp = me->risk;
    	if (tempcp > me->complexity)  tempcp = me->complexity;
	}
    else { //root node, and we already did all this stuff in rpart.c
    	tempcp = me->risk;
	}

    /*
    ** Can I quit now ?
    */
    if (me->num_obs < rp.min_split  ||  tempcp <= rp.alpha	|| nodenum > rp.maxnode) {
    	me->complexity =  rp.alpha;
    	me->leftson = (struct node *)0;
    	me->rightson= (struct node *)0;
    	*sumrisk = me->risk;
    	return(0);
	}

    /*
    ** Guess I have to do the split
    */
    bsplit(me, nodenum, variables_used);  /* ALG 1/24/2012: this figures out where to split, does the splitting! */

    if (me->primary ==0) {
		/*
		** This is rather rare -- but I couldn't find a split worth doing
		*/
		me->complexity = rp.alpha;
		me->leftson = (struct node *)0;
		me->rightson= (struct node *)0;
		*sumrisk = me->risk;
		return(0);
	}

    if (rp.maxsur>0) (void)surrogate(me, nodenum);
    else  me->surrogate =0;
    nodesplit(me, nodenum); /*ALG 1/26/2012: set the parameters of the 'split' struct */

    //ALG 1/31/2012: new array of variables used
    primary_split_var = (me->primary)->var_num;

    //ALG 3/18/2012: new way to store variables used
    int set_new;
    set_new = 0;
    for(i=0; i< rp.max_depth; i++){
    	int varnum;
    	varnum = variables_used[i];
    	if(varnum >= 0 || set_new == 1){  //we already added the new split var OR we're recording previous
    		vars_used_update[i] = varnum;
    	}
    	else{
    		vars_used_update[i] = primary_split_var;
    		set_new = 1;
    	}
    }

    /*
    ** split the leftson
    */
    me->leftson = (struct node *)CALLOC(1, nodesize);
    (me->leftson)->complexity = tempcp - rp.alpha;
    left_split = partition(2*nodenum, me->leftson, &left_risk, vars_used_update);

    /*
    ** Update my estimate of cp, and split the right son.
    */
    tempcp = (me->risk - left_risk) / (left_split +1);
    tempcp2 =(me->risk - (me->leftson)->risk);
    if (tempcp < tempcp2) tempcp = tempcp2;
    if (tempcp > me->complexity)  tempcp = me->complexity;

    me->rightson = (struct node *) CALLOC(1, nodesize);
    (me->rightson)->complexity = tempcp - rp.alpha;
    right_split = partition(1+2*nodenum, me->rightson, &right_risk, vars_used_update);


    /*
    ** Now calculate my actual C.P., which depends on children nodes, and
    **  on grandchildren who do not collapse before the children.
    ** The calculation is done assuming that I am the top node of the
    **  whole tree, an assumption to be fixed up later.
    */
    tempcp = (me->risk - (left_risk + right_risk))/
		    (left_split + right_split +1);


    /* Who goes first -- minimum of tempcp, leftson, and rightson */
	if ( (me->rightson)->complexity  > (me->leftson)->complexity ) {
	if (tempcp > (me->leftson)->complexity) {
		/* leftson collapses first */
		left_risk = (me->leftson)->risk;
		left_split =0;

		tempcp = (me->risk - (left_risk + right_risk)) /
				(left_split + right_split +1);
		if (tempcp > (me->rightson)->complexity) {
		/* right one goes too */
		right_risk = (me->rightson)->risk;
		right_split=0;
		}
		}
	}
	else if (tempcp > (me->rightson)->complexity) {
	/*right hand child goes first */
	right_split =0;
	right_risk = (me->rightson)->risk;

	tempcp = (me->risk - (left_risk + right_risk)) /
			(left_split + right_split +1);
	if (tempcp > (me->leftson)->complexity) {
		/* left one goes too */
		left_risk = (me->leftson)->risk;
		left_split=0;
		}
	}


    me->complexity= (me->risk - (left_risk + right_risk))/
			(left_split + right_split +1);

    if(!rp.collapse_is_possible){  //one sided criteria
    	double mm;
    	mm =fmin(left_risk, right_risk);
    	me->complexity = (me->risk - mm)/(left_split + right_split +1);
    }


	//3/18/2013: added second condition
    if ((me->complexity <= rp.alpha) && !rp.collapse_is_possible) {
	/*
	** All was in vain!  This node doesn't split after all.
	*/
    	//Rprintf("---collapsing!!!----\n");
		free_tree(me->leftson, 1);
		free_tree(me->rightson,1);
		me->leftson = (struct node *)0;
		me->rightson= (struct node *)0;
		*sumrisk = me->risk;
		return(0);             /*return # of splits */
	}
    else { /*we actually do the split */
        *sumrisk = left_risk + right_risk;
        if(!rp.collapse_is_possible){
        	*sumrisk = fmin(left_risk, right_risk);
        }
        return(left_split +right_split +1);
	}

}
