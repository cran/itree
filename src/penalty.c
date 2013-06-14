/* ALG 4/4/2012
/*
** The routines for computing penalties.
*/
#include <stdio.h>
#include "rpart.h"
#include "rpartS.h"
#include "node.h"
#include "rpartproto.h"

/* ALG 4/4/2012: no penalty*/
double no_penalty(int thisvar, int vars_used[]){
	return(0.0);
}

/* ALG 4/4/2012: equal penalty on any new variable
 * This just penalizes you by alpha if this variable
 * hasn't been used before in this branch.
 */
double newvar_penalty(int thisvar, int vars_used[]){
    /* ALG 3/18/2012  updated method
     * of keeping track of variables_used
    */
	double penalty_this_variable;
    int last_var_used, i;
    i = 0;
    last_var_used = vars_used[0];

    if(last_var_used == -1){  //we're at the root
    	penalty_this_variable = 0;  //so no penalty
    }
    else{  //continue checking variables used til we reach this node
		penalty_this_variable = rp.splitparams[0]; //alpha
        while(last_var_used != -1 && i<rp.max_depth){
        	if(last_var_used == thisvar){ //reset to zero
        		penalty_this_variable = 0;
        	}
            i++;
        	last_var_used = vars_used[i];
        }
    }
    return(penalty_this_variable);
}

/* ALG 4/4/2012: exponential moving average*/
double ema_penalty(int thisvar, int vars_used[]){

    /* ALG 3/18/2012  updated method
     * of keeping track of variables_used
    */
	double penalty_this_variable, temp, temp2, coef;
	int last_var_used, i, j;
	i=0;
	double indic_vec[rp.max_depth]; //vector of indicators. indic_vec[i]=1 if ith split is on this var

    last_var_used = vars_used[0];
    if(last_var_used == -1){  //we're at the root
    	penalty_this_variable = 0;  //so no penalty
    }

    else{
    	//first figure out the indicator vector
		while(last_var_used != -1 && i<rp.max_depth){
			if(last_var_used == thisvar){
				indic_vec[i]=1;
			}
			else{
				indic_vec[i]=0;
			}
			i++;
			last_var_used = vars_used[i];
		}
		i=i-1;

	    //now i = depth of this node = # of previous splits.
	    //Get exponential penalty by going backwards from depth
		temp = 0;
		temp2 = 0;  //sum of EMA coefficients
	    for(j=i; j>=0; j--){
	    	coef = pow(1-rp.splitparams[0],i-j);
	    	temp += indic_vec[j]* coef;
	    	temp2 += coef; //
	    }
	    penalty_this_variable = (temp2-temp)*rp.splitparams[0];
    }
    return(penalty_this_variable);
}
