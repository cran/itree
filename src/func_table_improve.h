/*
 * ALG 4/2/2012:
 *
 * Functions for calculating improvement. Typically we split on the variable/split pt
 * combination that maximizes 'improve'. Improve is a function of 'best'- the best
 * value of the objective function over all split pts for this variable,
 * the penalty for this variable, and which objective function we're
 * using.
 *
 * This allows us to generate the improve
 * automatically rather than hand-doing each potential group. Combinations that don't
 * make sense as well as defaults are taken care of at the R level before calling C.
 *
 *
*/

/* ALG 4/19/2012 */
extern double scale_by_root(double parent_risk,
									double best, double penalty);

/* ALG 4/19/2012 */
extern double scale_by_parent(double parent_risk,
									double best, double penalty);

/*
extern double anova_scale_by_parent(double root_risk ,double parent_risk,
									double best, double penalty);

extern double anova_scale_by_root(double root_risk ,double parent_risk,
									double best, double penalty);

extern double anova_pvalue(double root_risk ,double parent_risk,
									double best, double penalty);

extern double anova_pvalue(double root_risk ,double parent_risk,
									double best, double penalty);

extern double gini_scale_by_parent(double root_risk ,double parent_risk,
									double best, double penalty);

extern double gini_scale_by_root(double root_risk ,double parent_risk,
									double best, double penalty);

*/

static struct {
	double (*improve_fcn)();  //gets the improvement for this split variable
	}
	func_table_improve [] = {
	{scale_by_parent},
	{scale_by_root}
};

#define NUM_IMPROVE 2
