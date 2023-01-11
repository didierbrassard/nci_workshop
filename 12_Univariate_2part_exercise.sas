 /*************************************************************************/
 /*                                                                       */
 /*                   National Cancer Institute Method                    */
 /*                                                                       */
 /*               Code 12. Univariate, same as 12, TWO PART               */
 /*                                                                       */
 /*                        Author: Didier Brassard                        */
 /*                                                                       */
 /*                           Version 2-Exercise                          */
 /*                               01NOV2022                               */
 /*                                                                       */
 /* EXERCISE 01: Proportion of zero                                       */
 /* EXERCISE 02: Fill in MIXTRAN macro parameters                         */
 /* EXERCISE 03: Calculate variance ratio                                 */
 /* EXERCISE 04: Data manipulation                                        */
 /* EXERCISE 05: Fill in DISTRIB macro parameters                         */
 /*                                                                       */
 /*************************************************************************/

/* Indicate file location to work with relative location*/
	
	/* ## TO DO: Update according to your username or file location */
	
	%let path = /home/didier.brassard.10/nci_workshop/ ;

/* Include NCI macros */
	%include "&path./macros/mixtran_macro_v2.21.sas"; /* Univariate method, boxcox + error model + data manipulation */
	%include "&path./macros/distrib_macro_v2.2.sas"; /* Univariate method, distribution estimation*/

/* Define the `data` and `results` library */
	options dlcreatedir;
	libname data "&path./data/";
	libname results "&path./results/";
	
 /*************************************************************************/
 /*                                                                       */
 /*               Prepare an input data set for NCI macros                */
 /*                                                                       */
 /*************************************************************************/

	/* note: data create in program <01_Basics.sas>. Show 5 observations. */

	proc print data=data.prenci (obs=5) ;run;
	
/* create a copy for analysis */
	data preNCI;
		set data.preNCI;
	run;

 /************************************************/
 /*      Look at distribution on a given day     */
 /************************************************/

	/* ### EXERCISE 01: Assess the proportion of
		individuals who reported 0 gram of milk the first 24-h dietary recall */
	
	/* ### TIP: You can copy the code under <Preliminary analysis: Descriptive statistics>
		in script 01_Basics.sas to perfom this task */
	

 /*************************************************************************/
 /*                                                                       */
 /*                    NCI. MIXTRAN, POOLED, TWO PART                     */
 /*                                                                       */
 /*************************************************************************/

/* Objective: estimate usual intake distribution for an
	episodically-consumed food (milk) among adults 19-30 y and by sex (males, females). */

/*	### EXERCISE 02: in the MIXTRAN call below,
		A) add the dummy covariable for sequence and weekend,
		B) add the dummy covariable for age and sex, 
		C) indicate the error model is a correlated model (two-part, Pr*Amount)
		*/

%MIXTRAN(data          = preNCI, 
		 response      = milk,
		 foodtype      = milk,
		 subject       = adm_rno,
		 repeat        = recallid,
		 covars_prob   = ,
		 covars_amt    = ,
		 outlib        = work,
		 modeltype     = ,
		 lambda        = , 
		 replicate_var = WTS_P, 
		 seq           = ,
		 weekend       = ,
		 titles=0, printlevel=3);
		 
	/* note: model may run up to 1-2 minutes before convergence */

 /************************************************/
 /*              Variance and ratios             */
 /************************************************/

	/* ### EXERCISE 03: in the data step below, calculate the within:between,
		as well as the within:total variance ratios and show the results */

	/* coding tip: parameter estimates from the correlated model (probability, amount, corr(probability, amount))
		are stored in the data named: _param_<foodtype>, i.e., _param_milk in this example */

	data variance_ratio ;
	set work._param_milk; 

	/* ### ADD CALCULATION HERE */
	
	run;
	
 /*************************************************************************/
 /*                                                                       */
 /*                  NCI. MIXTRAN, STRATIFIED, TWO PART                   */
 /*                                                                       */
 /*************************************************************************/

/* ### EXERCISE 04: in a data step, stratify the <preNCI> data according to sex.
	Create one data named <male> when sex=1 and one data named <female> when sex=2 */

/* Stratification of data by sex */

	/* ### ADD DATA STEP HERE */

/* Call the MIXTRAN macro in each strata */
%MIXTRAN(data          = male, 
		 response      = milk,
		 foodtype      = male,
		 subject       = adm_rno,
		 repeat        = recallid,
		 covars_prob   = young ,
		 covars_amt    = young ,
		 outlib        = work,
		 modeltype     = CORR, /* CORR to fit correlated model */
		 lambda        = , 
		 replicate_var = WTS_P, 
		 seq           = seq2 ,
		 weekend       = weekend,
		 titles=0, printlevel=3);
	
%MIXTRAN(data          = female, 
		 response      = milk,
		 foodtype      = female,
		 subject       = adm_rno,
		 repeat        = recallid,
		 covars_prob   = young ,
		 covars_amt    = young ,
		 outlib        = work,
		 modeltype     = CORR, /* CORR to fit correlated model */
		 lambda        = , 
		 replicate_var = WTS_P, 
		 seq           = seq2 ,
		 weekend       = weekend,
		 titles=0, printlevel=3);

 /************************************************/
 /*              Variance and ratios             */
 /************************************************/

	data variance_ratio ;
	set work._param_male (keep=A_VAR: in=a)
		work._param_female (keep=A_VAR: in=b);
	if a then sex=1;
	else if b then sex=2;
	within_between = A_VAR_E / A_VAR_U2;
	within_total = A_VAR_E / (A_VAR_E+A_VAR_U2);
	run;
	
	/* check results */
	proc format;
	value sexfmt 
	1 = "Male"
	2 = "Female";
	run;
	
	proc print data=variance_ratio label;
	title1 "Variance and ratios";
	format sex sexfmt. A_VAR_E A_VAR_U2 within_between within_total 4.2 ;
	id sex;
	label A_VAR_E = "Within-ind. variance"
		  A_VAR_U2 = "Between-ind. variance"
		  within_between = "Within / Between ratio"
		  within_total = "Within / Total ratio"
		  ;
	run;
	title1;

/* Question: what can we conclude regarding error structure in males vs. female? */

 /*************************************************************************/
 /*                                                                       */
 /*          NCI. MIXTRAN data preparation, STRATIFIED, TWO PART          */
 /*                                                                       */
 /*************************************************************************/

/* create a `milk` library for results specific to this analysis */
	options dlcreatedir;
	libname milk "&path./results/nci_milk/";

/* append parameter estimates */
	data milk.param_all ;
	set work._param_male (in=a)
		work._param_female (in=b);
	if a then sex=1;
	else if b then sex=2;
	run;
	
/* append predicted values */
	data milk.pred_all ;
	set work._pred_male (in=a)
		work._pred_female (in=b);
	if a then sex=1;
	else if b then sex=2;
	run;
	
	/* note: there should be as many rows (observations or n)
		as the number of participants */

	proc freq data=preNCI;
	title1 "Confirm sample size";
	table recallid ;
	run;
	title1;

 /*************************************************************************/
 /*                                                                       */
 /*                  NCI. DISTRIB, STRATIFIED, TWO PART                   */
 /*                                                                       */
 /*************************************************************************/


/*	### EXERCISE 05: in the DISTRIB call below,
		A) indicate the error model is a correlated model,
		B) indicate the data for predicted values and parameter estimates,
		C) add the stratification variable in <byvar>
		*/

title1 "Usual intake distribution estimates for milk (grams)";
title2 "Adults 19-30 y from CCHS 2015";

%DISTRIB(call_type  = FULL, /* FULL to simulate usual intake AND estimate distribution */
		 seed       = 1234, /* starting point of random values */
		 nsim_mc    = 250, /* the number of pseudo-individuals */
		 mcsimda    = mcsim, 
		 modeltype  = , 
		 pred       = , 
		 param      = , 
		 outlib     = milk, /* other lib can be specified to save output */
		 cutpoints  = 125 250 500, /* cutpoints for which we want to know Pr(X<x) */
		 ncutpnt    = 3, /* n cutpoints used  */
		 byvar      = , 
		 add_da     = preNCI, /* same as before stratification*/
		 subgroup   = sex,
		 subject    = adm_rno,
		 titles     = 2,
		 food       = milk /* used for naming datasets */
		 );

 /*************************************************************************/
 /*                                                                       */
 /*                      NCI. DISTRIB, VISUALIZATION                      */
 /*                                                                       */
 /*************************************************************************/

/* transpose the descript_<food>_wts_p data */
	proc transpose data=milk.descript_milk_wts_p out=distrib_t prefix=sex_;
	var tpercentile1-tpercentile99 ;
	id sex;
	run;

/* add percentile index, rename variables */
	data distrib_t;
		set distrib_t;
	* remove character from <_NAME_>, change to numeric;
		percentile = input( compress(_NAME_,,'a'), 10.);
	*rename (old_name = new_name);
		rename 
		sex__255 = all
		sex_1 = male
		sex_2 = female;
	run;

/* generate cumulative distribution plot */

title1 "Cumulative distribution of usual milk intake (g/d), adults 19-30 y";
title2 italic "Data from the CCHS 2015 - Nutrition";
	proc sgplot data=distrib_t;
	series x=all y=percentile / legendlabel="Both" lineattrs=(color=black thickness=3);
	series x=male y=percentile / legendlabel="Male" lineattrs=(color=blue thickness=3);
	series x=female y=percentile / legendlabel="Female" lineattrs=(color=salmon thickness=3);
	run;
title1;
title2;

 /*************************************************************************/
 /*                                                                       */
 /*                       NCI. DISTRIB, PROPORTION                        */
 /*                                                                       */
 /*************************************************************************/

/* show proportion of all participants with milk intake >=125, 250 or 500 grams */

/* re-scale proportion from Pr(X<x) to Pr(X>=x) */
	data cutprob_inverse;
	set milk.descript_milk_wts_p ;
	array prop(*) cutprob: ;
	do i=1 to dim(prop);
		prop(i) = 1-prop(i);
	end;
	keep sex numsubjects cutprob: ;
	run;

/* Generate a format for clarity */
	proc format;
	value sexfmt
	-255 = "Males and females, 19-30 y"
	1 = "Males, 19-30 y"
	2 = "Females, 19-30 y";
	run;

/* print results with formatting */
	proc print data=cutprob_inverse label;
	title1 "The proportion of adults 19-30 y from Canada in 2015 with daily usual milk intakes (g) greater or equal to...";
	format sex sexfmt. cutprob: percent6.4;
	id sex numsubjects;
	var cutprob: ;
	label cutprob1 = "Pr(X>=125g)"
		cutprob2 = "Pr(X>=250g)"
		cutprob3 = "Pr(X>=500g)"
		;
	run ;
	title1;
	
/* end of code 12_exercise */

