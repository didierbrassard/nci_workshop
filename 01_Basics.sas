 /*************************************************************************/
 /*                                                                       */
 /*                   National Cancer Institute Method                    */
 /*                                                                       */
 /*    Code 01. Univariate: data preparation, measurement error model     */
 /*                                                                       */
 /*                        Author: Didier Brassard                        */
 /*                                                                       */
 /*                               Version 1                               */
 /*                               31OCT2022                               */
 /*                                                                       */
 /*************************************************************************/

/* Indicate file location to work with relative location*/
	
	/* ## TO DO: Update according to your username or file location */
	
	%let path = /home/didier.brassard.10/nci_workshop/ ;

/* Include NCI macros */
	%include "&path./macros/nlmixed_univariate_macro_v1.2.sas"; /* Univariate method, error model */

/* Define the `data` library */
	options dlcreatedir;
	libname data "&path./data/";

 /*************************************************************************/
 /*                                                                       */
 /*               Prepare an input data set for NCI macros                */
 /*                                                                       */
 /*************************************************************************/

	/* note: complete data pertaining to the 2015 Canadian Community Health Survey -
		Nutrition (Public Use Microdata Files) are available upon request to
		Statistics Canada online: https://www150.statcan.gc.ca/n1/en/catalogue/82M0024X2018001 */
		
/* 1) import preprocessed HS_NCI data using <proc import> */
	PROC IMPORT DATAFILE="&path./data/hsnci_19to30y.xlsx"
		DBMS=XLSX
		OUT=WORK.hsnci  replace;
		GETNAMES=YES;
	RUN;
	
	proc sort data=hsnci;
		by adm_rno;
	run;

/* 2) Additional data manipulation needed for analysis, save */

	/* note: NCI macros always require that recall-specific dietary intakes 
		appear on unique rows, i.e., the data has m records per participant
		where m = max number of 24hr */

	data data.preNCI;
		set hsnci;
	* create indicator variable for sequence of 24-h recall;
	seq2 = 0;
	if recallid=2 then seq2 = 1;
	
	* create indicator for age group (<25 y or >=25 y);
	if not missing(age) then do;
		young = 0;
		if age <25 then young =1 ;
	end;
	
	* create indicator for sex ;
	male=0;
	if sex=1 then male=1;
	
	run;
	
	/* 2.1 confirm data recoding */
	title1 "Prelim. Confirm recoding of covariates";
		proc freq data=data.preNCI;
		table recallid * seq2
			  sex * male
			  age * young
			  ;
		run;
	title1;

  	/* warning: <proc freq> usually not
  		appropriate for raw continuous var. */ 
 
	/* note: for simplicity, the CCHS design is completely ignored (sampling weights, variance).
		Thus, estimates and variance below would not be adequate for a formal analysis. */

 /*************************************************************************/
 /*                                                                       */
 /*             Preliminary analysis: Descriptive statistics              */
 /*                                                                       */
 /*************************************************************************/

/* 1) Distribution of intake by 24-h recall */
	title1 "Prelim. Distribution of X";
	proc sgpanel data=data.preNCI;
	  panelby recallid / rows=2 columns=1 layout=panel headerattrs=(weight=bold);
	  histogram energy ;
	  density energy / type=normal;
	  density energy / type=kernel;
	run;

/* 2) Summary of distribution */
	proc means data=data.preNCI n nmiss mean std min p25 p50 p75 max maxdec=0;
	class recallid;
	var energy;
	run;

/* 3a) How many `zeros` ? */
	
	/* note: useful to define a format to check without creating new variables */
	
	proc format;
		value zerofmt
		. = "Missing"
		low-<0.1 = "Zero"
		0.1-high = "Non-zero"
		;
	run;

	title1 "Prelim. Proportion of zeros for X";
	proc freq data=data.preNCI;
	format energy zerofmt. ;
	table energy *recallid;
	run;

	/* question: what should we do with participants with zero energy intake? */

/* 3b) How many `zeros` ? Example 2, episodic food */
	proc freq data=data.preNCI;
	format egg milk soft_drink zerofmt. ;
	table (egg milk soft_drink) * recallid;
	run;
	title1;
	
	/* question: what should we do when >30% of participants do not consume the food on a given day? */


 /*************************************************************************/
 /*                                                                       */
 /* Preliminary analysis: (Classical) measurement error model assumption  */
 /*                                                                       */
 /*************************************************************************/

	/* Technical note.
		NCI methods uses regression calibration (RC). RC assumes that the
		measured dietary intakes on a 24-h recall follow a classical
		measurement error model. According to this model,
		
		W = X + E
		
		where w = intakes measured via 24-h recall,
			  x = true intakes,
			  E = random errors following N(0,sigma2), iid
		*/

	/* Assumption 1. Measurement errors are normally distributed */
	/* Assumption 2. Measurement error variance does not depend on values of X (`true` intakes)*/
	/* Assumption 3. Measurement error variance does not depend on values of Z (any covariates)*/
	
	/* Reference:
		Carroll RJ, Ruppert D, Stefanski LA, Crainiceanu CM. Measurement Error
		in Nonlinear Models: A Modern Perspective, Second Edition: CRC Press; 2006. */

	/* note: NCI methods mostly take care of these assumptions
		(e.g., using best boxcox transformation). Also, most study uses 24-h recall
		on non-consecutive days which increases the probability that errors are independant. */
	
/* 1) perform some manipulation to create a `wide` data */

	/* 1.1) Remove 24-h recall with zero energy intake from the <preNCI> data */
	data data.preNCI;
		set data.preNCI;
	if energy>0 then output;
	run;
	
	/* note: it is common practice to remove 24hr with zero energy.
		Could also use other cut-off, e.g., 100, 250 kcal or some value that could
		plausibly be reported on a given day. */
	
	/* 1.2 Transpose preNCI */
	proc transpose data=data.preNCI out=energy_wide prefix=energy_r;
	by adm_rno ;
	var energy;
	id recallid;
	copy age sex; *to test assumption 2;
	run;


/* 2) calculate within-subject mean and std of intakes */
	data energy_wide1;
		set energy_wide(where=(not missing(energy_r1))); * delete repeated rows due to copying age;
	
	* remove participant without repeated measures (ie, any missing values);
	if nmiss(of energy_r1-energy_r2)>0 then delete;
	
	* mean of all 24-h dietary recalls;
	energy_r_mean = mean(of energy_r1-energy_r2);
	
	* standard deviation of all 24-h dietary recalls;
	energy_r_std = std(of energy_r1-energy_r2);
	
	* difference among repeated measures ;
	if not missing(energy_r2) then energy_delta2_1 = energy_r2-energy_r1;
		else energy_delta2_1=.;
	run;
	
/* 3) Check assumption 1 - differences among replicates = within-individual (random) errors */
	proc sgplot  data=energy_wide1;
	title1 "Are differences among replicates approximately normally distributed?";
	histogram energy_delta2_1 ;
	density energy_delta2_1 / type=normal ;
	density energy_delta2_1 / type=kernel ;
	run;
	ods select plots ;
	proc univariate data=energy_wide1 normal plot ;
	var energy_delta2_1 ;
	run;
	title1;

/* 4) Check Assumption 2 - Measurement error does not depend on X */

	/* note: StD(W) reflects variance of both `true` intakes and (random) errors */

	proc sgplot data=energy_wide1;
	title1 "Is there any trend between Mean(W) and StD(W)?";
	loess x=energy_r_mean y=energy_r_std / lineattrs=(color=red) markerattrs=(color=gray) ;
	run;
	
	ods select SpearmanCorr;
	proc corr data=energy_wide1 spearman;
	var energy_r_mean;
	with energy_r_std;
	run;
	title1;
	
/* 5) Check Assumption 3 - StD(W) reflects variance of both `true` intakes and (random) errors */
	proc sgplot data=energy_wide1;
	title1 "Is there any trend between StD(W) and Z?";
	loess x=energy_r_std y=age / lineattrs=(color=red) markerattrs=(color=gray);
	run;
	
	ods select SpearmanCorr;
	proc corr data=energy_wide1 spearman;
	var energy_r_std;
	with age;
	run;
	
	/* note: this assumption is required for each covariates, but only age is shown here. */

 /*************************************************************************/
 /*                                                                       */
 /*            NCI. BoxCox + Measurement error model, ONE PART            */
 /*                                                                       */
 /*************************************************************************/

	/* note: for this example, we use the <NLMIxed_Univariate> macro which is
		much more limited than the `full` univariate method (i.e., MIXTRAN+DISTRIB). */

	/* **************************************************************************** */
	/*   Macro Parameters:                                                          */
	/*                                                                              */
	/*      Required Parameters:                                                    */
	/*         data          = name of SAS data set containing the data to be       */
	/*                         analyzed. The data set has multiple observations     */
	/*                         for each subject, one for each reptition of the      */
	/*                         24-hour recall (or other dietary instrument).        */
	/*         subject       = name of the variable that uniquely identifies each   */
	/*                         subject (i.e., ID variable).                         */
	/*         repeat        = name of the variable that indexes repeated           */
	/*                         observations for each subject.                       */
	/*         response      = name of the food/nutrient variable to be modeled     */
	/*                         (24-hour recall variable for the food/nutrient).     */
	/*         modeltype     = model for food/nutrient:                             */
	/*                         to fit the two-part (epsisodic) model, specify       */
	/*                            modeltype = TWOPART                               */
	/*                         to fit the one-part (every day) model, specify       */
	/*                            modeltype = ONEPART                               */
	/*                                                                              */
	/* 	     Optional Parameters:                                                   */
	/* 	        covars_prob   = list of variables that are covariates in the        */
	/* 	                        probability part of the two-part model.             */
	/* 	                        if modeltype=ONEPART, then covars_prob is ignored.  */
	/* 	        covars_amt    = list of variables that are covariates in the        */
	/* 	                        one-part model or the amount part of the            */
	/* 	                        two-part model.                                     */
	/*                                                                              */
	/* **************************************************************************** */
	
/* Call the NCI macro <NLMIxed_Univariate>, no transformation */
title1 "Non-linear mixed measurement error model for one variable";
title2 "W = Energy intake, Z = sequence, weekend, age group, sex - UNTRANSFORMED";
	
	/* note: ignore <he maximum log likelihood was found for either the first or last LAMBDA= value> */
	
%NLMIxed_Univariate (data           = data.preNCI,
                     subject        = adm_rno, /* Unique subject id */
                     repeat         = recallid, /* 24-h recall identifier */
                     response       = energy, /* the dietary variable of interest */
                     modeltype      = ONEPART, /* Ok, since reported mostly by everyone*/
                     covars_prob    = , /* not used for one-part model */
                     covars_amt     = seq2 weekend young male, /* covariates, dummy-coded or continuous variables */
                     lambda         = 1, /* 1 = no transformation */
                     replicate_var  = , /* sampling weight (CCHS) or nothing if not a survey sample */
                     print          = Y, /* Y = show results */
                     ntitle         = 2
                     );
	
	/* question: what do parameters mean? */
	
	/* note: <work.parms_u> data include model paramater estimates */

	/* ************************************************************************ */
	/* Output Data Sets:                                                        */
	/*                                                                          */
	/*   parms_u = data set containing parameter estimates for the model.       */
	/*             parms_u contains the following variables:                    */
	/*                                                                          */
	/*                 A_Intercept = intercept in the amount part of the model. */
	/*                 A_varname   = regression slope for covariate "varname"   */
	/*                               in the amount part of the model.           */
	/*                 A_LogSDe    = Log(Sqrt(Var_e))                           */
	/*                 LogSDu2     = Log(Sqrt(Var_u2))                          */
	/*                 Var_e       = variance of the within-person error in the */
	/*                               amount part of the model.                  */
	/*                 Var_u2      = variance of the random effect in the       */
	/*                               amount part of the model.                  */
	/*                                                                          */
	/* ************************************************************************ */

	proc print data=parms_u ;
	run;

	/* question: what is the within:between ratio? Why is this relevant */
	
	data parms_u_ratio;
		set parms_u;
	within_between = VAR_E/VAR_U2;
	within_total = VAR_E / (VAR_E + VAR_U2) ;
	run;
	
	title1 "Variance ratios";
	proc print data=parms_u_ratio;
	format within_between num4.2 within_total percent6.2;
	var var: within: ;
	run;
	title1;	

/* Same model, but lets now use a log transformation to better satisfy assumptions */
title1 "Non-linear mixed measurement error model for one variable";
title2 "W = Energy intake, Z = sequence, weekend, age group, sex - LOG-TRANSFORMED";
%NLMIxed_Univariate (data           = data.preNCI,
                     subject        = adm_rno,
                     repeat         = recallid,
                     response       = energy,
                     modeltype      = ONEPART,
                     covars_prob    = ,
                     covars_amt     = seq2 weekend young male, /* dummy-coded or continuous variables */
                     lambda         = 0, /* 0=log, 1 =no transformation, or nothing = macro finds best boxcox */
                     replicate_var  = , /* sampling weight (CCHS) or nothing if not a survey sample */
                     print          = Y,
                     ntitle         = 2
                     );

	/* note: parameters are now expressed on the log-transformed scale */
                     
	data parms_u_ratio;
		set parms_u;
	within_between = VAR_E/VAR_U2;
	within_total = VAR_E / (VAR_E + VAR_U2) ;
	run;
	
	title1 "Variance and ratios, log-transformed X";
	proc print data=parms_u_ratio;
	format within_between num4.2 within_total percent6.2;
	var var: within: ;
	run;
	title1;
	
	/* question: the within:between variance ratio is now 2.0 (vs. 1.8 based on untransformed data),
		what does this tell us? */


 /*************************************************************************/
 /*                                                                       */
 /*   Bonus: estimation of within-individual coefficient of variation     */
 /*                                                                       */
 /*************************************************************************/
	
	/* Why? Needed when estimating the plausbility of reported energy intakes */

	/* Technical note: 
		the squared root of the variance is the Standard Deviation (SD) which, when
		experessed on the log scale, is the coefficient of variation  */
		
	/* Reference:
		Cole TJ, Altman DG. Statistics Notes: Percentage differences, symmetry,
	    and natural logarithms. BMJ. 2017;358:j3683. */

	%NLMIxed_Univariate (data           = data.preNCI,
	                     subject        = adm_rno, /* Unique subject id */
	                     repeat         = recallid,
	                     response       = energy,
	                     modeltype      = ONEPART,
	                     covars_prob    = ,
	                     covars_amt     = seq2 , /* why only sequence as covariates ? */
	                     lambda         = 0, /* 0=log */
                     	 replicate_var  = WTS_P, /* sampling weight (CCHS) or nothing if not a survey sample */
	                     print          = N,
	                     ntitle         = 2
	                     );

	data CV_rEI;
		set parms_u;
	sqrt_VAR_E = sqrt(VAR_E);
	CV_rEI = sqrt_VAR_E*100;
	run;
	
	title1 "Within-individual coefficient of variation (CV) in reported energy intake (EI)" ;
	proc print data=CV_rEI;
	format VAR_E sqrt_VAR_E CV_rEI num4.2;
	var VAR_E sqrt_VAR_E CV_rEI;
	run;
	title1;
	
	/* interpretation: we see that in this sample, reported energy intake varies
		approx. 39% across 24-h recalls on average */

 /*************************************************************************/
 /*                                                                       */
 /*     Bonus 2: number of 24-h recalls to average for usual intakes      */
 /*                                                                       */
 /*************************************************************************/

	/* Technical note: within- and between-individual variance are the 
		major contributors of day-to-day differences in intakes (Beaton et al. AJCN 1979).
		If we know the relative contribution of within to between variance, we can
		estimate the number of 24-h recalls that we would need to collect to estimate
		usual intakes within a set margin of error. */

	/* According to the formula by Beaton et al.,
	n = (Z * CV / L) squared
	where  n = number of days, Z = normal deviate, CV = within-individual CV, L = arbitrary limit
	*/
	
	data number_of_24h;
	set CV_rEI(keep=CV_rEI) ;
	z= 1.96;
	cv=CV_rEI ;
	L = 20 ; * lets say, within 20 percent of usual energy intake (eg, 2500 +/- 500 kcal) ;
	n = (z * cv/L)**2 ;
	run;
	
	proc print data=number_of_24h;
	title1 'Number of (raw) 24-h recall needed to estimate usual energy intake within 20%, 95% of the time';
	run;
	title1;

/* end of code 01 */

