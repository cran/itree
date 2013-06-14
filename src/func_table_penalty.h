/*
 * ALG 4/2/2012:
 * Functions for calculating penalties.
 *
 * Split criteria is a function
 * of penalty fcn AND purity/SS/etc (objective) AND improve fcn. This allows us to generate triples
 * automatically rather than hand-doing each potential group. Combinations that don't
 * make sense as well as defaults are taken care of at the R level before calling C.
*/
#ifndef FLOAT
#define FLOAT float   /*see comments in rpart.h */
#endif


/* ALG 4/4/2012: no penalty*/
extern double no_penalty(int thisvar, int vars_used[]);

/* ALG 4/4/2012: equal penalty on any new variable*/
extern double newvar_penalty(int thisvar, int vars_used[]);

/* ALG 4/4/2012: exponential moving average*/
extern double ema_penalty(int thisvar, int vars_used[]);

static struct {
	double (*penalty_fcn)();  //computes the penalty for splitting this variable at this node
}
	func_table_penalty [] = {
	{no_penalty},  //original
	{newvar_penalty},  //alpha*1(is this variable used in this branch)
	{ema_penalty}   // exponential moving avg, prefers recently used variables
};

#define NUM_PENALTY 3
