/* ALG 4/4/2012
/*
** The routines for computing improvement. By definition
** for ALL methods/penalties, we split at variable+splitpt combination
** that yields greatest improve.
**
** These functions are simple but de-coupling them from the penalty and objective
** functions allows us far more flexibility and saves us from lengthy
** switch statements.
*/
#include <stdio.h>
#include "rpart.h"
#include "rpartS.h"
#include "node.h"
#include "rpartproto.h"


double scale_by_root(double parent_objective,double best_objective, double penalty){
	//best should be of the form
	// objective(parent) - objective(left & right children)
	//here we divide by root's objective
	double impr;
	impr = (best_objective/rp.root_objective_scaling)-penalty;
	return(impr);
}

double scale_by_parent(double parent_objective,double best_objective, double penalty){
	//best should be of the form
	// objective(parent) - objective(left & right children)
	//here we divide by parent's objective
	double impr;
	if(best_objective > 0){
		impr = (best_objective/parent_objective)-penalty;
	}
	else{
		impr = -1*penalty;
	}
	return(impr);
}


//double anova_pvalue(double root_risk ,double parent_risk,
//									double best, double penalty){
//	return(best);
//}
