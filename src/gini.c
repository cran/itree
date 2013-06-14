/*   @(#)gini.c	1.13 11/08/01    
** The routines for gini-classification
**
** ALG 4/12/2012
** Calculate impurity reduction for classification.
**
** ALG 9/18/2012: despite the name of the class, this implements
** Gini AND cross-entropy. Default is Gini.
*/
#include <math.h>
#include <stdio.h>
#include "rpart.h"
#include "rpartS.h"
#include "rpartproto.h"

static int    numclass;
static double *left,     /*left branch n (weighted)*/
	      *right,
	      **ccnt;
static double *prior,
	      *aprior,   /*altered priors */
	      *freq,     //alg.  freq is supposed to hold the weighted number of observ. for each class
              *loss;      /* loss matrix */
static int    *tsplit,
	      *countn;
static double *awt,
	      *rate;

static double (*impurity)();

//ALG 4/12/2012: Gini: squared error loss on +/-1 - fitted prob
//9/19/2012: don't need p(1-p) b/c adding up all p's over all classes
//just gives 1.
static double gini_impure1(p) double p; {  return(1 - p*p); }

//alg 4/12/2012: Cross-entropy
static double gini_impure2(p)
double p; { if (p==0) return(0.0); else return(-p*log(p)); }

int giniinit(int n,        double **y, int maxcat, char **error, 
	     double *parm, int *size,  int who,    double *wt)
    {
    int i, j, k;
    double temp;

    /* allocate memory  and setup losses */
    if (who==1) {
        numclass =0;   /*number of classes */
        for (i=0; i<n; i++){
        	if (*y[i] > numclass)  numclass = *y[i];
        }

        /* ALG 4/12/2012: Choose which impurity function to use
         */
        if (parm[numclass + numclass*numclass] ==2){
        	impurity = gini_impure2;
        }
        else{
        	impurity = gini_impure1;
        }
 
        left = (double *) ALLOC(numclass*2, sizeof(double));
        right = left+numclass;

        tsplit= (int *) ALLOC(maxcat*2, sizeof(int));
        countn= tsplit + maxcat;

		awt = (double *) ALLOC(maxcat*2, sizeof(double));
		rate= awt +maxcat;

		if (maxcat>0) {
			graycode_init0(maxcat);
			ccnt    = (double **) ALLOC(numclass, sizeof(double *));
			if (ccnt==0) {*error=_("Out of memory"); return(1);}
			ccnt[0] = (double *) ALLOC(numclass*maxcat, sizeof(double));
			if (ccnt[0]==0) {*error=_("Out of memory"); return(1);}
			for (i=1; i<numclass; i++)
			ccnt[i] = ccnt[i-1] + maxcat;
		}

		i = 3*numclass + numclass*numclass;
		prior = (double *) ALLOC(i, sizeof (double));
		if (prior==0) {*error=_("Out of memory"); return(1);}
		aprior = prior + numclass;
		freq   = aprior+ numclass;
		loss   = freq  + numclass;

		for (i=0; i<numclass; i++)  freq[i] =0;
		temp =0;
		for (i=0; i<n; i++) {
			j = *y[i] -1;
			freq[j] += wt[i];
			temp += wt[i];   /*sum total of weights */
	    }
		for (i=0; i<numclass; i++)  freq[i] /=temp;   /*relative frequency */

		temp =0;
		for (i=0; i<numclass; i++) {
			prior[i] = parm[i];
			aprior[i] =0;
			for (j=0; j<numclass; j++) {
				k = numclass*i + j;
				loss[k] = parm[numclass+k];
				temp += loss[k] * prior[i];
				aprior[i] += loss[k] * prior[i];
			}
	    }

		for (i=0; i<numclass; i++) {
			if (freq[i]>0) {  /* watch out for a missing class */
				prior[i] /= freq[i];
				aprior[i] /= (temp * freq[i]);  /* pi_i / n_i */
			}
		}
	}//end if(who==1)
    
	*size = 1 + numclass;
	return(0);
	}

/*
** Compute the predicted response and the classification error
**   This is R(T) (this node's contribution) in the paper.
*/
void ginidev(int n, double **y, double *value, double *risk, double *wt)
    {
    int i, j, max = 0;
    double  temp, dev;

    dev =0;
    /*
    ** count up number in each class
    */
    for (i=0; i<numclass; i++)  freq[i]=0;
    for (i=0; i<n; i++) {
    	j = y[i][0] -1;
    	freq[j] += wt[i];
	}

    /*
    ** Now compute best class and its error
    */
    for (i=0; i<numclass; i++) {  /* assume class i is the predicted class */
		temp =0;

		for (j=0; j<numclass; j++) {
			temp += freq[j] * loss[j*numclass +i] *prior[j];
		}

		//figure out if this is the best class so far.
		if (i==0 || temp < dev) {
			max =i;
			dev = temp;
		}
	}
    //the predicted class is the one with minimum loss.
    value[0] = max +1;    /* remember: external groups start at 1 */

    //here we store for each class, the weighted sum of that class in the node
    //...nothing to do with priors or loss.
    for (i=0; i<numclass; i++) value[i+1] = freq[i];
    *risk  = dev;
}


/*
 * ALG 5/30/2012. The function that calculates
 * how what we scale 'improve' by.
 * For gini it's whatever impurity function we're using, NOT the deviance=miscl error.
 */
void gini_parent_objective(int n, double *y[], double *value, double *parent_objective,
	     double *wt) {

    int i,j;
    double rwt;
    int  rtot;
    double total_ss, temp;

    rwt=0;
    rtot=0;

    //reset right
    for(i=0; i<numclass; i++){
    	right[i] = 0;
    }

    //put everything to the right to start
    for (i=0; i<n; i++) {
    	j = *y[i] -1;   //actual value
    	rwt += aprior[j] * wt[i];    /*altered weight = prior * case_weight */
    	right[j] += wt[i];
		rtot++;
	}

    total_ss =0;
    for (i=0; i<numclass; i++){
    	temp = aprior[i] * right[i]/ rwt;      /* p(class=i, given node A) */
    	total_ss += rwt * (*impurity)(temp);    /* p(A) * I(A). this is total node weight * impur fcn */
    }
    *parent_objective = total_ss;
}


/*
** return the error for a particular point
*/
double ginipred(double *y, double *pred)
    {
    int i, j;
    double temp;
    i = y[0] -1;
    j = *pred -1;
    temp = prior[i]*loss[i*numclass +j];
    return(temp);
}


/*
** The gini splitting function.  Find that split point in x such that
**  the rss within the two groups is decreased as much
**  as possible.
**
**  ALG 4/9/2012: removed args added 1/21/2012
**  ALG 4/30/2012: added parent_objective
*/
void gini(int n,    double *y[],     FLOAT *x,     int numcat,
	  int edge, double *improve, FLOAT *split, int *csplit, double my_risk,
	  double *wt, double *parent_objective)
    { 
    int i,j,k;
    double lwt, rwt;
    int  rtot, ltot;
    int    direction = LEFT, where = 0;
    double total_ss, best, temp, p;

    double lmean, rmean;    /* used to decide direction */

    //reset right & left counts of y classes
    for (i=0; i<numclass; i++) {
    	left[i] =0;
    	right[i]=0;
	}
    lwt =0;  rwt=0;
    rtot=0;  ltot=0;

    //put everything to the right to start
    for (i=0; i<n; i++) {
    	j = *y[i] -1;   //actual value
    	rwt += aprior[j] * wt[i];    /*altered weight = prior * case_weight */
    	right[j] += wt[i];
		rtot++;
	}

    total_ss =0;
    for (i=0; i<numclass; i++){
    	temp = aprior[i] * right[i]/ rwt;      /* p(class=i, given node A) */
    	total_ss += rwt * (*impurity)(temp);    /* p(A) * I(A). this is total node weight * impur fcn */
	}
    best =total_ss; /* total weight of right * impurity of right + 0 *0 */

    /*
    ** at this point we split into 2 disjoint paths
    */
    if (numcat >0) goto categorical;

    //cts predictor
    for (i=0;  rtot >edge; i++) {
	j = *y[i] -1;
	rwt -= aprior[j] * wt[i];
	lwt += aprior[j] * wt[i];
	rtot--;
	ltot++;
	right[j] -= wt[i];
	left[j]  += wt[i];

	if (x[i+1] != x[i] &&  (ltot>=edge)) {
	    temp =0;
            lmean =0; rmean =0;
	    for (j=0; j<numclass; j++) {
		p = aprior[j]*left[j]/lwt;    /* p(j | left) */
		temp += lwt * (*impurity)(p);      /* p(left) * I(left) */
		lmean += p*j;
                p =  aprior[j]*right[j]/rwt;   /* p(j | right) */
		temp += rwt * (*impurity)(p);      /*p(right) * I(right) */
		rmean += p*j;
	        }
	    if (temp < best) {
		best=temp;
		where =i;
		if (lmean < rmean) direction = LEFT;
		    else           direction = RIGHT;
		}
	    }
	}

    //this is the impurity of parent - min(weight(right)*Impurity(right)+weight(left)*Impurity(left)
    //The weights are btwn (0,n) not (0,1).
    *improve =  (total_ss - best);
    if (*improve > 0 ) {   /* found something */
    	csplit[0] = direction;
    	*split = (x[where] + x[where+1]) /2;
	}
    return;

categorical:;
    /*
    ** First collapse the data into a numclass x numcat array
    **  ccnt[i][j] = number of class i obs, category j of the predictor
    */
    for (j=0; j<numcat; j++) {
	awt[j] =0;
	countn[j]=0;
	for (i=0; i<numclass; i++)
	    ccnt[i][j] =0;
	}
    for (i=0; i<n; i++) {
	j = *y[i] -1;
	k = x[i] -1;
	awt[k] += aprior[j] * wt[i];
	countn[k]++;
	ccnt[j][k] += wt[i];
	}

    for (i=0; i<numcat; i++){ 
	if (awt[i]==0) tsplit[i] =0;
	else {
	    rate[i] = ccnt[0][i] / awt[i];   /* a scratch array */
	    tsplit[i]=RIGHT;
	    }
        }

    if (numclass==2) graycode_init2(numcat, countn, rate);
                else graycode_init1(numcat, countn);

    while((i=graycode()) < numcat) {
	/* item i changes groups */
	if (tsplit[i]==LEFT) {
	    tsplit[i]=RIGHT;
	    rwt  += awt[i];
	    lwt -= awt[i];
	    rtot += countn[i];
	    ltot -= countn[i];
	    for (j=0; j<numclass; j++) {
		right[j] += ccnt[j][i];
		left[j]  -= ccnt[j][i];
    	        }
	    }
	else {
	    tsplit[i]=LEFT;
	    rwt -= awt[i];
	    lwt += awt[i];
	    rtot -= countn[i];
	    ltot += countn[i];
	    for (j=0; j<numclass; j++) {
		right[j] -= ccnt[j][i];
		left[j]  += ccnt[j][i];
	    }
	}

	if (ltot>=edge  &&  rtot>=edge) {
	    temp =0;
	    lmean=0; rmean =0;
	    for (j=0; j<numclass; j++) {
		p = aprior[j]*left[j] /lwt;
		temp +=  lwt * (*impurity)(p);
		lmean += p*j;
                p =  aprior[j]*right[j]/rwt;       /* p(j | right) */
		temp += rwt * (*impurity)(p);      /*p(right) * I(right) */
		rmean += p*j;
	        }
	    if (temp < best) {
		best=temp;
		if (lmean < rmean)
			for (j=0; j<numcat; j++) csplit[j] = tsplit[j];
		else
			for (j=0; j<numcat; j++) csplit[j] = -tsplit[j];
	        }
	    }
        }
    *improve = (total_ss - best);
    }

