/* SCCS @(#)rpart.h	1.10 06/06/01   */
/*
** commom variables for the rpart routine
*/
#include <R.h>
#include <math.h> //ALG 4/28/2012 because rp.alpha sometimes is set to -Infinity
#undef error

#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("rpart", String)
#else
#define _(String) (String)
#endif

#ifndef FLOAT
#define FLOAT double    /*so I can easily change 'x' to double at some later
			 date, with all the consequences thereof.  Also see
		         node.h  */
#endif
#define LEFT  (-1)     /*used for the variable "extra" in nodes */
#define RIGHT  1
#define MISSING 0

#ifdef MAINRP
#define EXTERN
#else
#define EXTERN extern
#endif
/* As a sop to S, I need to keep the total number of external symbols
**  somewhat smaller.  So, pack most of them all into a structure.
*/
EXTERN struct {
    double complexity; //cp parameter
    double alpha;  //cp parameter * risk at the root. used to see if we can stop splitting.
    double iscale;         /* used to check improvement==0, with error */
    double **ydata;
    FLOAT  **xdata;
    FLOAT  *xtemp;
    double *wt;
    double **ytemp;
    double *wtemp;          /* temp vector of weights */
    double *lwt;
    double *rwt;            /*scratch double vectors, of length ncat */
    double *vcost;          /* variable costs */
    Sint    *numcat;        /* variable type: 0=cont, 1+  =#categories */
    Sint   **sorts;             /* allocated on the fly */
    int    n;              /* total number of subjects  */
    int    num_y;          /* number of y variables */
    int    nvar;           /* number of predictors */
    int    maxpri;
    int    maxsur;   /* max # of primary or surrogate splits to use */
    int    usesurrogate;
    int    num_unique_cp;   /*ALG 1/26/2012: this defines which method to use */
    int    min_node;       /* minimum size for any terminal node */
    int    min_split;      /*minimum size before we attempt a split */
    int    num_resp;       /*length of the response vector */
    int    sur_agree;      /*0 =  my style, 1=CART style */
    int    maxnode;        /*controls the maximum depth of the tree */
    int    *tempvec;       /*to be allocated by the mainline, of length n */
    int    *which;
    int    *csplit;
    int    *left;
    int    *right;
    int    *featvec;  	/* ALG 1/16/2012: vector of length nvar to pass to nodes to keep track of split criteria */
    double *splitparams; /* ALG 1/21/2012: vector with parameters needed for path-dependent split criteria */
    double root_objective_scaling;  //ALG 5/1/2012. Value we scale improve by for the root node. Used for scaling if impscale=root.
    int method_number; /*ALG 3/5/2012: which split method*/
    int impscale_number;/*ALG 4/27/2012*/
    int penalty_number;/*ALG 4/27/2012*/
    int max_depth;  //most number of splits
    double *dummy;   //ALG 7/19/2012. for consistency we need to pass a dummy double pointer sometimes
    //double split_check_offset; //4/28/2012: added so that we can insure cp doesn't matter for some methods
    int collapse_is_possible; //ALG 9/21/2012: usually 1 for TRUE. ensures we don't collapse purity/extremes
    						 //due to 'complexity' value computed as per the usual CART criteria.
    }  rp;

EXTERN struct cptable *cptable_tail;
EXTERN int  (*rp_init)();    /*called to initialize a splitting function */
EXTERN void (*rp_choose)();  /*set to the splitting function */
EXTERN void (*rp_eval)();   /*set to the evaluation routine */
EXTERN double (*rp_error)();     /*set to the prediction error routine */

//ALG 5/1/2012: function that returns objective to scale this node's improve by
EXTERN void (*rp_parent_objective)();

//ALG 4/4/2012: added penalty function
EXTERN double (*rp_penalty)(); /* set to the penalty routine */

//ALG 4/11/2012: improve function function
EXTERN double (*rp_improve)(); /* set to the improve routine */

/*
** The user inputs his complexity parameter as a percentage. and the
**   printout is also scaled in this way.  The book and the computations all
**   have an easier time with absolute cp.  So complex = what the user
**   typed and alpha = complex * (risk of top node) = what is used
**   internally.
** The variable 'complex' in node.h is also on the internal scale.
**
** Categorical variables must be coded as 1,2,3, ..., and there may be
**  missing categories.  The upper limit is determined on the fly.
*/
