/* SCCS @(#)func_table.h	1.5 06/06/01  */
/*
 * ALG 4/2/2012.
 * Renamed to "objective".
*
*  The objective function is the statistic we are trying to minimize
*  or maximizing before introducing a penalty.  For instance
*  anova minimizes SS.
**
**  init_split   - Will be called before a tree is started.  May do very
**                  little, but empirical Bayes like methods will need to set
**                  some global variables.
**  choose_split - function to find the best split
**  eval         - function to calculate the response estimate and risk
**  error        - Function that returns the prediction error.
**  num_y        - Number of columns needed to represent y (usually 1)
**
**  ALG 5/1/2012: added functions that return the scale factor for each method
**  when normalizing the 'improve'. Except for purity and extremes, this is the
**  same function eval.  They ensure the penalty term always can be expressed as a
**  number between (0,1) in all contexts.
*/

#ifndef FLOAT
#define FLOAT float   /*see comments in rpart.h */
#endif

extern int anovainit( int n,	    double*y[],  int maxcat, char **error, 
		      double *param, int *size,   int who,    double *wt);
extern int anovainit_interp( int n,	    double*y[],  int maxcat, char **error,
		      double *param, int *size,   int who,    double *wt);
extern int anovainit_interp2( int n,	    double*y[],  int maxcat, char **error,
		      double *param, int *size,   int who,    double *wt);
extern int anovainit_interp3( int n,	    double*y[],  int maxcat, char **error,
		      double *param, int *size,   int who,    double *wt);

//ALG 3/3/2012
extern int anovainit_extremes( int n,	    double*y[],  int maxcat, char **error,
		      double *param, int *size,   int who,    double *wt);

//ALG 3/21/2012
extern int purity_regression_init( int n,	    double*y[],  int maxcat, char **error,
		      double *param, int *size,   int who,    double *wt);

//ALG 9/19/2012
extern int purity_classification_init( int n,	    double*y[],  int maxcat, char **error,
		      double *param, int *size,   int who,    double *wt);

extern int classification_extremes_init( int n,	    double*y[],  int maxcat, char **error,
		      double *param, int *size,   int who,    double *wt);


//ALG 3/25/2012
extern int anovainit_intExp( int n,	    double*y[],  int maxcat, char **error,
		      double *param, int *size,   int who,    double *wt);

extern int poissoninit(int n,	     double*y[],  int maxcat, char **error, 
		       double *parm, int *size,   int who,    double *wt);
extern int    giniinit(int n,	     double*y[],  int maxcat, char **error, 
		       double *parm, int *size,   int who,    double *wt);
extern int usersplit_init(int n,     double*y[],  int maxcat, char **error, 
		       double *parm, int *size,   int who,    double *wt);

extern void anovass   (int n,     double *y[], double *value, double *risk,
                                  double *wt);
extern void anovass_interp   (int n,     double *y[], double *value, double *risk,
                                  double *wt);

extern void anovass_interp2   (int n,     double *y[], double *value, double *risk,
                                  double *wt);

extern void anovass_interp3   (int n,     double *y[], double *value, double *risk,
                                  double *wt);

extern void anova_extremes_eval(int n,     double *y[], double *value, double *risk,
                                  double *wt);

extern void classification_extremes_eval(int n,     double *y[], double *value, double *risk,
                                  double *wt);


extern void purity_regression_eval(int n,     double *y[], double *value, double *risk,
                                  double *wt);

extern void purity_classification_eval(int n,     double *y[], double *value, double *risk,
                                  double *wt);

extern void anovass_intExp(int n,     double *y[], double *value, double *risk,
                                  double *wt);


extern void poissondev(int n,     double *y[], double *value, double *risk,
                                  double *wt);
extern void ginidev   (int n,     double *y[], double *value, double *risk,
                                  double *wt);
extern void usersplit_eval(int n,     double *y[], double *value, double *risk,
                           double *wt);

/*ALG 1/21/2012: added 2 arguments to the end of these (the third line) to
 * match the call for anova_interp  */
/*ALG 1/24/2012: added 3rd to last arg to pass which variable we're looking at */
/* ALG 1/26/2012: changed vars_used to int[] */
/* ALG 4/9/2012:  removed arguments added 1/21 ...that functionality is part of bsplit */
/* ALG 4/11/2012:  added penalty argument.. pre-computed by bsplit */
/* ALG 4/11/2012:  took away penalty argument... done later in bsplit */
extern void anova(    int n,    double *y[],     FLOAT *x,     int nclass, 
                      int edge, double *improve, FLOAT *split, int *csplit, 
                      double myrisk,             double *wt);
extern void poisson(  int n,    double *y[],     FLOAT *x,     int nclass, 
                      int edge, double *improve, FLOAT *split, int *csplit, 
                      double myrisk,             double *wt);
extern void gini(     int n,    double *y[],     FLOAT *x,     int nclass, 
                      int edge, double *improve, FLOAT *split, int *csplit, 
                      double myrisk,             double *wt);
extern void usersplit(int n,    double *y[],     FLOAT *x,     int nclass, 
                      int edge, double *improve, FLOAT *split, int *csplit, 
                      double myrisk,             double *wt);


/* ALG 2/15/2012: interp3. p-value approach */
extern void anova_interp3(int n,    double *y[],     FLOAT *x,     int nclass,
        int edge, double *improve, FLOAT *split, int *csplit,
        double myrisk, double *wt);

/* ALG 3/3/2012: high/low mean approach from Buja & YS Lee */
extern void anova_extremes(int n,    double *y[],     FLOAT *x,     int nclass,
        int edge, double *improve, FLOAT *split, int *csplit,
        double myrisk, double *wt);

/* ALG 9/19/2012: */
extern void classification_extremes(int n,    double *y[],     FLOAT *x,     int nclass,
        int edge, double *improve, FLOAT *split, int *csplit,
        double myrisk, double *wt);


/* ALG 3/21/2012: one sided purity for regression from Buja & YS Lee */
extern void purity_regression(int n,    double *y[],     FLOAT *x,     int nclass,
        int edge, double *improve, FLOAT *split, int *csplit,
        double myrisk, double *wt);

/* ALG 3/21/2012: one sided purity for classification from Buja & YS Lee */
extern void purity_classification(int n,    double *y[],     FLOAT *x,     int nclass,
        int edge, double *improve, FLOAT *split, int *csplit,
        double myrisk, double *wt);

extern double anovapred  (double *y, double *yhat);
extern double  ginipred  (double *y, double *yhat);
extern double poissonpred(double *y, double *yhat);
extern double usersplit_pred(double *y, double *yhat);
extern double purity_classification_pred(double *y, double *yhat);
extern double classification_extremes_pred(double *y, double *yhat);


/*ALG 5/1/2012: these scaling factor for improve, same args as eval functions,
 * but they're put in for consistency*/
extern void anova_parent_objective(int n, double *y[], double *value, double *risk, double *wt);
extern void poisson_parent_objective(int n, double *y[], double *value, double *risk, double *wt);
extern void gini_parent_objective(int n, double *y[], double *value, double *risk, double *wt);
extern void extremes_parent_objective(int n, double *y[], double *value, double *risk, double *wt);
extern void purity_parent_objective(int n, double *y[], double *value, double *risk, double *wt);
extern void usersplit_parent_objective(int n, double *y[], double *value, double *risk, double *wt);
extern void purity_classification_parent_objective(int n, double *y[], double *value, double *risk, double *wt);
extern void extremes_classification_parent_objective(int n, double *y[], double *value, double *risk, double *wt);

/*ALG 5/1/2012 added _parent_objective functions*/
/*ALG 9/19/2012 added classification purity*/
static struct {
	int  (*init_split)();
	void (*choose_split)();
	void (*eval)();
	double (*error)();
	void (*parent_objective)();
	  } func_table_objective [] =
		 {{anovainit, anova, anovass, anovapred,anova_parent_objective},
		  {poissoninit, poisson, poissondev, poissonpred,poisson_parent_objective},
		  {giniinit, gini, ginidev, ginipred,gini_parent_objective},
		  {usersplit_init, usersplit, usersplit_eval, usersplit_pred,usersplit_parent_objective},
		  {anovainit_extremes, anova_extremes, anovass, anovapred,extremes_parent_objective},
		  {purity_regression_init, purity_regression,purity_regression_eval,anovapred,purity_parent_objective},
		  {classification_extremes_init, classification_extremes,
				 classification_extremes_eval, classification_extremes_pred,
				 extremes_classification_parent_objective},
		  {purity_classification_init, purity_classification, purity_classification_eval,
				  purity_classification_pred,purity_classification_parent_objective}
		 };

#define NUM_METHODS 8
