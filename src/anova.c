/* SCCS @(#)anova.c	1.8 08/13/01  */
/*
** The four routines for anova splitting
*/
#include <stdio.h>
#include "rpart.h"
#include "rpartS.h"
#include "node.h"
#include "rpartproto.h"

static double *mean, *sums;
static double *wts;
static int *countn;
static int *tsplit;

int anovainit(int n,        double *y[],  int maxcat, char **error, 
	      double *parm, int *size,    int who,    double *wt)
    {
    if (who==1 && maxcat >0) {
	graycode_init0(maxcat);
	countn  = (int *)ALLOC(2*maxcat, sizeof(int));
	tsplit  = countn + maxcat;
	mean   = (double *)ALLOC(3*maxcat, sizeof(double));
	wts    = mean + maxcat;
	sums   = wts + maxcat;
	}
    *size =1;
    return(0);
    }
/*
** The anova evaluation function.  Return the mean and the ss.
*/
void anovass(int n, double *y[], double *value, double *risk,
	     double *wt) {
    int i;
    double temp, twt;
    double nodemean, ss;

    temp =0;
    twt  =0;   		/* sum of the weights */

    for (i=0; i<n; i++) {
    	temp += *y[i] * wt[i];
    	twt  += wt[i];
	}
    nodemean = temp/twt;

    ss =0;
    for (i=0; i<n; i++) {
    	temp = *y[i] - nodemean;
    	ss += temp*temp* wt[i];
	}

    *value = nodemean;
    *risk = ss;
}

/*
 * ALG 5/30/2012. The function that calculates
 * what we scale 'improve' by.
 * For anova it's just the SS, so we repeat the eval function.
 *
*/
void anova_parent_objective(int n, double *y[], double *value, double *parent_objective,
	     double *wt) {
	anovass(n,y,value,parent_objective, wt);
}

/*
** The anova splitting function.  Find that split point in x such that
**  the sum of squares of y within the two groups is decreased as much
**  as possible.  It is not necessary to actually calculate the SS, the
**  improvement involves only means in the two groups.
**
**  ALG 1/21/2012: added the last two arguments to match the anova_interp
**  call
**
**  ALG 4/2/2012: start making flexible to passing all types of penalties.
**  ALG 4/9/2012: add computation of penalty.
**
**  ALG 4/9/2012: removed args added 1/21/2012
**  ALG 4/11/2012: added penalty - (pre-computed)
**  ALG 4/19/2012: took away penalty--- done in bsplit
*/
void anova(int n,    double *y[],     FLOAT *x,     int nclass, 
	   int edge, double *improve, FLOAT *split, int *csplit, 
	   double myrisk, double *wt)
    {
    int i,j;
    double temp;
    double left_sum, right_sum;
    double left_wt, right_wt;
    int    left_n,  right_n;
    double grandmean, best;
    int direction = LEFT;
    int where = 0;

    /*
    ** The improvement of a node is SS - (SS_L + SS_R), where
    **   SS = sum of squares in a node = \sum w_i (x_i - \bar x)^2, where
    ** of course \bar x is a weighted mean \sum w_i x_i / \sum w_i
    ** Using the identity 
    **    \sum w_i(x_ - \bar x)^2 = \sum w_i (x_i-c)^2 - (\sum w_i)(c-\bar x)^2
    ** the improvement = w_l*(left mean - grand mean)^2 
    **                  +w_r*(right mean- grand mean)^2
    ** where w_l is the sum of weights in the left node, w_r similarly. 
    **
    */
    right_wt =0;
    right_n  = n;
    right_sum =0;
    for (i=0; i<n; i++) {
    	right_sum += *y[i] * wt[i];
    	right_wt  += wt[i];
	}
    grandmean = right_sum/right_wt;

    if (nclass==0) {   /* continuous predictor */
        left_sum=0;   /* No data in left branch, to start */
	left_wt =0;   left_n =0;
	right_sum=0;  /*after subracting grand mean, it's zero */
	best    = 0;
	for (i=0; right_n>edge; i++) {
	    left_wt += wt[i];  
	    right_wt -= wt[i];
	    left_n++;
	    right_n--;
	    temp = (*y[i] - grandmean) * wt[i];
	    left_sum  +=temp;
	    right_sum -=temp;

	    //alg 3/3/2012: check if we're allowed to do this split
	    if (x[i+1] !=x[i] &&  left_n>=edge) {
			temp = left_sum*left_sum/left_wt  +
					right_sum*right_sum/right_wt;
			if (temp > best) {
				best=temp;
				where =i;
				if (left_sum < right_sum) direction = LEFT;
						  else    direction = RIGHT;
				}
	    }
	 }//end loop through each split

	*improve =  best;  //already in the form of parent - best
	if (best>0) {   /* found something */
	    csplit[0] = direction;
	    *split = (x[where] + x[where+1]) /2;
	    }
	} /* end of continuous predictor */

    /* 
    ** Categorical predictor 
    */
    else {
	for (i=0; i<nclass; i++) {
	    sums[i] =0;
	    countn[i]=0;
	    wts[i] =0;
	    }

	/* rank the classes by their mean y value */
	for (i=0; i<n; i++) {
	    j = (int)x[i] -1;
	    countn[j]++;
	    wts[j] += wt[i];
	    sums[j] += (*y[i] - grandmean) * wt[i];
	    }
	for (i=0; i<nclass; i++)  {
	    if (countn[i] >0) {
		tsplit[i] = RIGHT;
		mean[i] = sums[i]/ wts[i];
		}
	    else tsplit[i] = 0;
	    }
	graycode_init2(nclass, countn, mean);

	/*
	** Now find the split that we want
	*/
	left_wt =0; 
	left_sum=0; right_sum=0;  
	left_n = 0; 
	best =0;
	where =0;
	while((j=graycode()) < nclass) {
	    tsplit[j] = LEFT;
	    left_n += countn[j];
	    right_n-= countn[j];
	    left_wt += wts[j];
	    right_wt-= wts[j];
	    left_sum += sums[j];
	    right_sum-= sums[j];
	    if (left_n>=edge  &&  right_n>=edge) {
		temp = left_sum *left_sum /left_wt  + 
		       right_sum*right_sum/right_wt;
		if (temp > best) {
		    best=temp;
		    if ((left_sum/left_wt) > (right_sum/right_wt)) {
			for (i=0; i<nclass; i++) csplit[i] = -tsplit[i];
			}
		    else {
			for (i=0; i<nclass; i++) csplit[i] = tsplit[i];
			}
		    }
		}
	    }
	// *improve = (*rp_improve)(root_risk,parent_risk,best,penalty);
	*improve = best;
	} /* end of categorical predictor */
}

