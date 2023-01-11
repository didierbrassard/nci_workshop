 /*************************************************************************/
 /*                                                                       */
 /*                   National Cancer Institute Method                    */
 /*                                                                       */
 /*  Code 11. Univariate: measurement error model, MC Sim., distribution  */
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

 /*************************************************************************/
 /*                                                                       */
 /*                NCI. MIXTRAN, POOLED VARIANCE, ONE PART                */
 /*                                                                       */
 /*************************************************************************/

	/* technical note: when the measurement error model is derived based on all
		participants, this assumes that the error structure is common to all
		individuals (Kirkpatrick et al. J Am Dietetic Assoc 2022). This is sometimes
		referred to as `pooled` analysis (vs. stratified) */ 
	
	/* note: not all macro parameters are shown below for simplicity. See
		<&path./macros/mixtran_macro_v2.21.sas> for complete list */

	/*******************************************************************************/
	/*                                                                             */
	/* The syntax for calling the macro is:                                        */
	/*                                                                             */
	/* %MIXTRAN                                                                    */
	/* (data=, response=, foodtype=, subject=, repeat=,                            */
	/*  covars_prob=, covars_amt=, outlib=, modeltype=,                            */
	/*  lambda=, replicate_var=, seq=,                                             */
	/*  weekend=, vargroup=, numvargroups=,                                        */
	/*  start_val1=, start_val2=, start_val3=, vcontrol=,                          */
	/*  nloptions=, titles=, printlevel=,subgroup=)                                */ 
	/*                                                                             */
	/*  where                                                                      */
	/*                                                                             */
	/*   "data"         * Specifies the dataset to be used.                        */
	/*                                                                             */
	/*   "response"     * Specifies the 24 hour recall variable.                   */
	/*                                                                             */
	/*   "foodtype"     * Specifies a name for the analysis, used to               */
	/*                    identify the output data sets.  This value can           */
	/*                    be the same as the response variable.                    */
	/*                                                                             */
	/*   "subject"      * Specifies the variable that uniquely                     */
	/*                    identifies each subject.                                 */
	/*                                                                             */
	/*   "repeat"       * Specifies the variable that indexes repeated             */
	/*                    observations for each subject.                           */
	/*                                                                             */
	/*   "covars_prob"    Specifies a list of covariates for the first             */
	/*                    part of the model that models the probability            */
	/*                    of consumption.  Covariates must be separated            */
	/*                    by spaces.  Interactions must be in the order            */
	/*                    specified by PROC GENMOD.  If the model type             */
	/*                    is "amount" then covars_prob should be left as           */
	/*                    a null string.                                           */
	/*                                                                             */
	/*   "covars_amt"   * Specifies a list of covariates for the second            */
	/*                    part of the model that models the consumption-           */
	/*                    day amount.  Covariates must be separated by             */
	/*                    spaces.  Interactions must be in the order               */
	/*                    specified by PROC GENMOD.                                */
	/*                                                                             */
	/*   "outlib"       * Specifies a directory where output data sets             */
	/*                    are stored.  Outlib can not be null.                     */
	/*                                                                             */
	/*   "modeltype"    * Specifies the model.  The possible values are:           */
	/*                    "null string" = fit correlated model,                    */
	/*                    "corr"        = fit correlated model,                    */
	/*                    "nocorr"      = fit uncorrelated model,                  */
	/*                    "amount"      = fit amount-only model.                   */
	/*                                                                             */
	/*   "lambda"         Specifies a user supplied value for the                  */
	/*                    Box-Cox transformation parameter, lambda.  If            */
	/*                    a value is not supplied, the macro will                  */
	/*                    calculate a value for lambda.                            */
	/*                                                                             */
	/*   "replicate_var"  Specifies the variable to be used in the                 */
	/*                    replicate statement of PROC NLMIXED or the               */
	/*                    freq statement of PROC UNIVARIATE.  The                  */
	/*                    specified variable must be integer valued.               */
	/*                                                                             */
	/*   "seq"            Specifies one or more sequence indicator                 */
	/*                    variables to account for effects due to the              */
	/*                    sequence number of a subject's records.  This            */
	/*                    variable can NOT also appear in covars_prob              */
	/*                    or covars_amt.                                           */
	/*                                                                             */
	/*   "weekend"        Specifies the weekend (Fri.-Sun.) indicator              */
	/*                    variable to account for a weekend effect.  A             */
	/*                    value of 1 represents a Fri.-Sun. record, and            */
	/*                    a value of 0 represents a Mon.-Thurs. record.            */
	/*                    This variable can NOT also appear in                     */
	/*                    covars_prob or covars_amt.                               */
	/*                                                                             */
	/*******************************************************************************/

title1 "Non-linear mixed measurement error model for one variable (POOLED)";
title2 "W = Energy intake, Z = sequence, weekend, age group, sex - BOXCOX TRANSFORMED";

%MIXTRAN(data          = preNCI, 
		 response      = energy, /* the dietary variable of interest */
		 foodtype      = energy, /* used for naming datasets */
		 subject       = adm_rno, /* unique participant identifier */
		 repeat        = recallid, /* 24-h recall index variable */
		 covars_prob   = , /* null for daily consumed food or nutrients */
		 covars_amt    = young male, /* with MIXTRAN: do NOT include sequence or weekend */
		 outlib        = work, /* other lib (e.g., results) can be specified to save output */
		 modeltype     = AMOUNT, /* suitable for daily consumed food or nutrients */
		 lambda        = , /* Manual lambda for transformation or null: macro finds best boxcox */
		 replicate_var = WTS_P, /* sampling weights or none (note: must be integer) */
		 seq           = seq2 , /* MIXTRAN: dedicated parameter to remove sequence effect */
		 weekend       = weekend, /* MIXTRAN: dedicated parameter to weight data according to weekend */
		 titles=2, printlevel=3);
		 
 /************************************************/
 /*            Summary of output data            */
 /************************************************/

	proc print data=work.etas_energy ;
	title1 "Measurement error model equation";
	run;
	title1;
	
	proc print data=work._parmsf2_energy;
	title1 "Measurement error model parameter estimates";
	run;
	title1;

	proc print data=work._param_unc_energy;
	title1 "Measurement error model parameter estimates + additional information";
	title2 "Additional information include min. amount, weighting varialbe (if any), weekend, variance estimates";
	run;
	title1;
	
	proc print data=work._pred_unc_energy(obs=5);
	title1 "Predicted value for each individual (5 observations)";
	run;
	title1;
	
	/* question: what do these values represent? Why _0 and _1 ? */

 /************************************************/
 /*              Variance and ratios             */
 /************************************************/

	/* calculate ratio of within:between and within:total variance */
	data variance_ratio ;
	set work._param_unc_energy(keep=A_VAR:);
	within_between = A_VAR_E / A_VAR_U2;
	within_total = A_VAR_E / (A_VAR_E+A_VAR_U2);
	run;
	
	/* check results */
	proc print data=variance_ratio label;
	title1 "Variance and ratios";
	format A_VAR_E A_VAR_U2 within_between within_total 4.2;
	label A_VAR_E = "Within-ind. variance"
		  A_VAR_U2 = "Between-ind. variance"
		  within_between = "Within / Between ratio"
		  within_total = "Within / Total ratio"
		  ;
	run;
	title1;
	
	/* question:
		what are we looking for?
		What are we expecting? 
		What if we suspect that the error structure may vary according to other characteristics (e.g., sex)
		*/
		
 /*************************************************************************/
 /*                                                                       */
 /*          NCI. MIXTRAN+DISTRIB, STRATIFIED VARIANCE, ONE PART          */
 /*                                                                       */
 /*************************************************************************/

	/* technical note: subgroup-specific measurement error model can be obtained
		by stratifying the data and repeated call for MIXTRAN. */

/* 1) Stratify the data according to sex */
	data male female;
	set preNCI;
	if sex=1 then output male;
	else if sex=2 then output female;
	run;

/* 2.1) Call MIXTRAN with first strata (sex=1, males)*/

title1 "Non-linear mixed measurement error model for one variable (STRATIFIED)";
title2 "W = Energy intake, Z = sequence, weekend, age group - BOXCOX TRANSFORMED";

%MIXTRAN(data          = male,  /* Data for current strata only */
		 response      = energy,
		 foodtype      = male, /* name consistent with data to avoid errors */
		 subject       = adm_rno,
		 repeat        = recallid, 
		 covars_prob   = , 
		 covars_amt    = young /*male*/,  /* question: why isnt the <male> variable used here ?*/
		 outlib        = work, 
		 modeltype     = AMOUNT,
		 lambda        = , 
		 replicate_var = WTS_P, 
		 seq           = seq2 ,
		 weekend       = weekend,
		 titles=2, printlevel=1); /* printlevel=1 to reduce the information printed */
		 
/* 2.2) Call MIXTRAN with first strata (sex=2, females)*/
%MIXTRAN(data          = female,  /* Data for current strata only */
		 response      = energy,
		 foodtype      = female, /* name consistent with data to avoid errors */
		 subject       = adm_rno,
		 repeat        = recallid, 
		 covars_prob   = , 
		 covars_amt    = young /*male*/, 
		 outlib        = work, 
		 modeltype     = AMOUNT,
		 lambda        = , 
		 replicate_var = WTS_P, 
		 seq           = seq2 ,
		 weekend       = weekend,
		 titles=2, printlevel=1); /* printlevel=1 to reduce the information printed */
		
 /************************************************/
 /*              Variance and ratios             */
 /************************************************/

	/* append measurement error model from both strata and
		calculate ratio of within:between and within:total variance */
		
	data variance_ratio ;
	set work._param_unc_male (keep=A_VAR: in=a)
		work._param_unc_female (keep=A_VAR: in=b);
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
	
	/* question: what can we conclude? Should use opt for a POOLED or STRATIFIED approach? */
	
	/* more information on stratification vs pooling:
		Early Experience Analyzing Dietary Intake Data from the Canadian Community Health Survey-Nutrition
		Using the National Cancer Institute (NCI) Method
		https://pubmed.ncbi.nlm.nih.gov/31443191/ 
		*/
		
	/* coding note: had we opted for a stratified approach,
		both predicted values data (_param_unc_male, _param_unc_female) and
		parameters estimates (_param_unc_male, _param_unc_male) should have been
		appended together (`set`) including an index variable (e.g., sex = [1, 2])
		before the next step. */

 /*************************************************************************/
 /*                                                                       */
 /*            NCI. MIXTRAN+DISTRIB, POOLED VARIANCE, ONE PART            */
 /*                                                                       */
 /*************************************************************************/

title1 "Non-linear mixed measurement error model for one variable (POOLED)";
title2 "W = Energy intake, Z = sequence, weekend, age group - BOXCOX TRANSFORMED";

	/* note: MIXTRAN call is the same as above */

%MIXTRAN(data          = preNCI,
		 response      = energy,
		 foodtype      = energy,
		 subject       = adm_rno,
		 repeat        = recallid, 
		 covars_prob   = , 
		 covars_amt    = young male, /* note: male covariable included in this call */
		 outlib        = work, 
		 modeltype     = AMOUNT,
		 lambda        = , 
		 replicate_var = WTS_P, 
		 seq           = seq2 ,
		 weekend       = weekend,
		 titles=2, printlevel=1);

	/*************************************************************************/
	/*                                                                  
	/* The syntax for calling the DISTRIB macro is:                     
	/*                                                                  
	/*       %DISTRIB (call_type=, seed=, nsim_mc=, modeltype=,         
	/*                 pred=, param=,outlib=, cutpoints=, ncutpnt=,     
	/*                 byvar=, subgroup=, add_da=, subject=,          
	/*                 titles=, food=, mcsimda=, 
	/*                 recamt=,recamt_co=,recamt_hi=,
	/*                 wkend_prop=,wkend_mc=);                                    
	/*                                                                  
	/* where:                                                           
	/*                                                                  
	/*  "call_type" * Specifies which parts of the DISTRIB macro should 
	/*                be invoked. (FULL, MC, PC). FULL is the default.
	/*                A null string implies FULL.
	/*                FULL = invoke both the calculation 
	/*                  of the estimated intake amount (using monte carlo 
	/*                  simulation) and of the percentiles (and optionally 
	/*                  the percent not meeting the recommended amount of
	/*                  intake).  
	/*                MC = restrict the macro to calculating 
	/*                  the intake amount (using monte carlo simulation).
	/*                PC = use the intake estimates - calculated 
	/*                  in the MC macro in DISTRIB - 
	/*                  to calculate percentiles and, optionally, the 
	/*                  percent meeting the recommended amount of 
	/*                  intake, and/or cutpoints. 
	/*
	/*  "seed"      * Specifies the seed for the random number          
	/*                generator used for the Monte Carlo simulation of  
	/*                the random effects u1 and u2.  
	/*                Required if the call_type is FULL or MC.                   
	/*                Not used if call_type is PC.                      
	/*                                                                  
	/*  "nsim_mc"   * Specifies the number of repetitions to be used in 
	/*                the Monte Carlo simulation.  For each subject,    
	/*                one record will be output for each repetition. 
	/*                Required if the call_type is FULL or MC.      
	/*                Not used if call_type is PC.                      
	/*                                                                  
	/*  "modeltype" * Specifies the model that was used by the MIXTRAN  
	/*                macro to prepare the data for the DISTRIB macro.  
	/*                The value must be the same as the model declared  
	/*                for the MIXTRAN macro. The default is correlated. 
	/*                The possible values are:  
	/*                  null string  = fit correlated model,             
	/*                  corr         = fit correlated model,             
	/*                  nocorr       = fit uncorrelated model,           
	/*                  amount       = fit amount-only model.            
	/*                                                                  
	/*  "pred"      * Specifies the name of the data set containing     
	/*                predicted values calculated in the MIXTRAN macro.                      
	/*                Required.             
	/*                                                                  
	/*  "param"     * Specifies the name of the data set containing the 
	/*                parameter estimates calculated in the MIXTRAN     
	/*                macro.  
	/*                Required if the call_type is FULL or MC.
	/*                Not used if call_type is PC.                                            
	/*                                                                  
	/*  "outlib"    * Specifies the library reference to which the      
	/*                output data set of distributions will be written. 
	/*                                                                  
	/*  "cutpoints"   Specifies a list of cutoff points separated by a  
	/*                space.                                            
	/*                Not used if call_type is MC.                      
	/*                                                                  
	/*  "ncutpnt"     Specifies the number of cutoff points.  If cutoff 
	/*                points are given, ncutpnt must also be given.     
	/*                Not used if call_type is MC.                      
	/*                                                                  
	/*  "byvar"       Specifies a list of by-variables that are in the  
	/*                data sets "pred" and "param" to indicate that the 
	/*                MIXTRAN model was fit separately for each         
	/*                by-group. The estimates used in the calculation will 
	/*                differ based on the by group, however The DISTRIB 
	/*                macro will ultimately produce estimates of the entire         
	/*                population (not distributions within each by group).          
	/*                To obtain distributions for subpopulations, use the                          
	/*                "subgroup" parameter.                             
	/*                                                                  
	/*  "subgroup"    Specifies one categorical variable used for the   
	/*                calculation of a separate usual intake            
	/*                distribution for each subgroup.  The distribution 
	/*                of usual intake will also be calculated for the   
	/*                overall data set. Requires that the paramater
	/*                add_da be supplied. 
	/*                Not used if call_type is MC.                      
	/*
	/*  "add_da"      The name of the data set containing the subgroup
	/*                variable by which the percentiles are to be
	/*                calculated, and/or the recommended amount variable.
	/*                This data set must include:  
	/*                the ID variable declared in the parameter SUBJECT,
	/*                and one or both of the following variables: 
	/*                the variable named in the parameter SUBGROUP 
	/*                the variable(s) named in the parameter(s) recamt and/or recamt_hi. 
	/*                This parameter is required if either of the parameters
	/*                subgroup or recamt are called.
	/*                Not used if call_type is MC.                      
	/*                                                                  
	/*  "subject" *   Specifies the variable that uniquely identifies   
	/*                each subject. (The ID.)                                    
	/*                                                                    
	/*                                                                  
	/*  "titles"      Specifies the number of title lines to be   
	/*                reserved for the user's titles.      
	/*                The default value is 0.                                                        
	/*                                                                  
	/*  "food"     *  Specifies a name for the analysis, used to name  
	/*                the output data set.                              
	/*                                                                  
	/*  "mcsimda"  *  Specifies the name of the data set containing the 
	/*                intake amount derived from the monte carlo        
	/*                simulations. To read or write the data file from disk include
	/*                a libname. Note: due to simulations (see parameter nsim_mc)  
	/*                the data set can grow quite large.               
	/*                Default value is work.mcsim in which case the data set is
	/*                not saved for later use.                                 
	/*                Required if call_type is PC.                                         
	/*                                                                       */
	/*************************************************************************/	

title1 "Usual intake distribution estimates for energy";
title2 "Adults 19-30 y from CCHS 2015";

%DISTRIB(call_type  = FULL, /* full to simulate usual intake AND estimate distribution */
		 seed       = 1234, /* starting point of random values */
		 nsim_mc    = 250, /* the number of pseudo-individuals */
		 mcsimda    = mcsim, /* name of the data with Monte Carlo Simulations */ 
		 modeltype  = AMOUNT, /* must be consistent with MIXTRAN */
		 pred       = _pred_unc_energy,  /* data with predicted values from MIXTRAN */
		 param      = _param_unc_energy, /* data with parameters from MIXTRAN */
		 outlib     = results, /* other lib can be specified to save output */
		 cutpoints  = 1500 2000 2500 3000 3500, /* cutpoints for which we want to know Pr(X<x) */
		 ncutpnt    = 5, /* 5 cutpoints used for this example */
		 byvar      = , /* had we used the STRATIFIED model, <sex> would go here */
		 add_da     = preNCI, /* same as <data> in MIXTRAN */
		 subgroup   = sex, /* distribution for subgroups used in the error model */
		 subject    = adm_rno,
		 titles     = 2,
		 food       = energy /* used for naming datasets */
		 );
		 
	/* warning: mcsim data can be quite large as they are size of 
		(n * ncsim_mc) = (1901 * 250 ) = 475,250 rows */
		
	/* technical note: with very large sample, using only 100 Monte Carlo simulations 
		per individual may be appropriate, while a larger number of simulations
		(pseudo-individuals) may by used with smaller sample.
		For example, 250 or 500 simulations. */

 /************************************************/
 /*             <MCSIM> data overview            */
 /************************************************/

/* Monte Carlo simulation data */
	proc print data=mcsim (obs=10) label noobs;
	label ADM_RNO = "ADM_RNO: subject id."
		  mcsim_wt = "mcsim_wt: Sampling weight, divided by numsims"
		  numsims = "numsims: Number of simulation per subject"
		  mc_t = "mc_t: Estimated usual intake"
		  mc_a = "mc_a: Estimated usual amount"
		;
	run;
	
 /************************************************/
 /*   <DESCRIPT> data overview (usual intakes)   */
 /************************************************/	

	/* note: data named as: descript_<food>_<replicate_var> */

/* Generate a format for clarity */
	proc format;
	value sexfmt
	-255 = "Males and females, 19-30 y"
	1 = "Males, 19-30 y"
	2 = "Females, 19-30 y";
	run;

/* Print mean usual intakes */
	title1 "Mean energy intake (NCI univariate method)";
	proc print data=results.descript_energy_wts_p ;
	format sex sexfmt. mean_mc_t 4.2;
	id sex;
	var mean_mc_t  ;
	run;
	title1;
	
	/* question: sometimes, the `NCI` mean is slightly
		different than the mean based on raw data. Why? */
	 
/* Print selected percentile estimate */
	title1 "Usual energy intake distribution estimates";
	proc print data = results.descript_energy_wts_p ;
	format sex sexfmt. tpercentile: 4.2 ;
	id sex ;
	var tpercentile1 tpercentile5 tpercentile25 tpercentile50 tpercentile75 tpercentile95 tpercentile99 ;
	run;
	title1;

/* Print cutpoints, i.e., Pr(X<x) */
	title1 "Proportion of individuals 19-30 y below thresholds";
	proc print data = results.descript_energy_wts_p label;
	format sex sexfmt. cut: percent6.4;
	id sex numsubjects;
	var cut:;
	label cutprob1 = "Pr(X<1500)"
		  cutprob2 = "Pr(X<2000)"
		  cutprob3 = "Pr(X<2500)"
		  cutprob4 = "Pr(X<3000)"
		  cutprob5 = "Pr(X<3500)";
	run;
	title1;

 /*************************************************************************/
 /*                                                                       */
 /*             Bonus: Examples of distribution visualization             */
 /*                                                                       */
 /*************************************************************************/

/* note: SAS does not have state-of-the-art graphical capabilities. Nonetheless,
	basic graph can be generated for visualization purpose */
	
/* 1) transpose data to have a long-formatted data */
	proc transpose data=results.descript_energy_wts_p out=distrib_t prefix=GR_;
	var tpercentile1-tpercentile99 ;
	id sex;
	run;

	/* 1.1) Apply some formatting */
	data distrib_t;
	set distrib_t;
	* add numerical value for each percentile;
	percentile = input(compress(_NAME_,,'a'),10.);
	rename GR__255 = GR_0 ;
	run;

/* 2) Histograms, density, cumulative function */
	proc sgplot data=distrib_t;
	title1 "Usual intake distribution (19-30 y)";
	histogram gr_0 ;
	density gr_0 / type=kernel;
	density gr_0 / type=normal ;
	xaxis label="Usual energy intake, kcal";
	run;
	title1; 

	proc sgplot data=distrib_t;
	title1 "Usual intake distribution (19-30 y)";
	density gr_1 /  legendlabel="Males" type=kernel lineattrs=(color=lightblue thickness=3);
	density gr_2 /  legendlabel="Females" type=kernel lineattrs=(color=salmon thickness=3);
	xaxis label="Usual energy intake, kcal";
	run;
	title1; 

	proc sgplot data=distrib_t;
	title1 "Cumulative distribution (19-30 y)";
	series y=percentile x=gr_0 / legendlabel="Both" lineattrs=(color=black thickness=3);
	series y=percentile x=gr_1 / legendlabel="Males" lineattrs=(color=lightblue thickness=3);
	series y=percentile x=gr_2 / legendlabel="Females" lineattrs=(color=salmon thickness=3);
	refline 1500 2000 2500 3000 3500 / label=("1500" "2000" "2500" "3000" "3500") 
		axis=X lineattrs=(pattern=LONGDASH);
	xaxis label="Usual energy intake, kcal";
	run;
	title1; 
	
/* 3) Compare distribution on a given day vs. usual */
	
	/* 3.1) Output distribution on a given day */
	proc univariate data=preNCI noprint;
	where recallid=1;
	var energy ;
	weight WTS_P ;
	output out=_day1 pctlpre=p pctlpts=(1 to 99 by 1) ;
	run;
	
	/* 3.2) transpose from wide to long */
	proc transpose data=_day1 out=day1_w ;
	var p1-p99 ;
	run;
	
	/* 3.3) some formatting */
	data day1_w(keep=energy_day1 percentile);
	set day1_w;
	* add numerical value for each percentile;
	percentile = input(compress(_NAME_,,'a'),10.);
	rename COL1=energy_day1;
	run;
	
	/* 3.4) add usual intake values (all) */
	data usual_vs_day1;
	retain percentile energy_day1 energy_usual ;
	merge distrib_t(keep=percentile GR_0 in=a rename=(GR_0=energy_usual))
		day1_w (in=b ) ;
	by percentile;
	* calculate difference at all percentiles values ;
	diff = energy_usual - energy_day1 ;
	diff_percent = (energy_day1 - energy_usual) / mean(energy_usual, energy_day1);
	run;
	
	/* note: according to theory, lower percentile are UNDER-estimated, while
		upper percentile are OVER-estimated.
		
		Intuitive explanation:
		on a given day, there will always be some individuals that
		ate a *little less* than usual, while there will always be some individuals that are 
		a *little more* than usual  */

	/* 3.5) plot distribution */
	proc sgplot data=usual_vs_day1 ;
	title1 "Energy intake distribution (19-30 y, CCHS-2015)";
	title2 italic "On a given day vs. usual" ; 
	footnote1 italic justify=right "Usual = NCI univariate method";
	density energy_day1 / legendlabel="Given day" type=kernel lineattrs=(color=red thickness=3);
	density energy_usual / legendlabel="NCI/Usual" type=kernel lineattrs=(color=darkgreen thickness=3);
	run;
	title1;
	title2;
	footnote1;

	proc sgplot data=usual_vs_day1 ;
	title1 "Cumulative energy intake distribution (19-30 y, CCHS-2015)";
	title2 italic "On a given day vs. usual" ;
	footnote1 italic justify=right "Usual = NCI univariate method";
	series x=energy_day1 y=percentile / legendlabel="Given day" lineattrs=(color=red thickness=3);
	series x=energy_usual y=percentile / legendlabel="NCI/Usual" lineattrs=(color=darkgreen thickness=3);
	refline 1500 2000 2500 3000 3500 / label=("1500" "2000" "2500" "3000" "3500") 
		axis=X lineattrs=(pattern=LONGDASH);	
	run;
	title1;
	title2;
	footnote1;
	
	proc sgplot data=usual_vs_day1 noautolegend;
	title2 "Compared with USUAL, energy intake estimates on A GIVEN DAY are misestimated by...";
	format diff_percent percent6.4;
	loess y=diff_percent x=percentile /lineattrs=(color=red thickness=3) markerattrs=(color=black);
	xaxis label="Percentile";
	yaxis label="Percentage bias"; 
	refline 0 / axis=y label ="No bias" lineattrs=(pattern=longdash);
	run;
	title1;

 /*************************************************************************/
 /*                                                                       */
 /*           Bonus: estimate distribution for another subgroup           */
 /*                                                                       */
 /*************************************************************************/

/* note: covariates used in the modelling (e.g., age group, sex) or in 
	stratification, had we stratified MIXTRAN, can be used to output distributions.
	Given that the mcsim data is still available, we can call distrib a second time
	with <call_type= PC> to estimate a distribution for a new subgroup. */

title1 "Usual intake distribution estimates for energy";
title2 "Adults 19-30 y from CCHS 2015";

%DISTRIB(call_type  = PC, /* PC to estimate distribution only */
		 mcsimda    = mcsim, /* name of the data with Monte Carlo Simulations */ 
		 modeltype  = AMOUNT, /* must be consistent with MIXTRAN */
		 pred       = _pred_unc_energy,  /* data with predicted values from MIXTRAN */
		 param      = _param_unc_energy, /* data with parameters from MIXTRAN */
		 outlib     = results, /* other lib can be specified to save output */
		 cutpoints  = 1500 2000 2500 3000 3500, /* cutpoints for which we want to know Pr(X<x) */
		 ncutpnt    = 5, /* 5 cutpoints used for this example */
		 add_da     = preNCI, /* same as <data> in MIXTRAN */
		 subgroup   = young, /* distribution for a different subgroup used in the error model */
		 subject    = adm_rno,
		 titles     = 2,
		 food       = energy_age /* used for naming datasets */
		 );

/* warning: it would be incorrect to specify a subgroup NOT considered in MIXTRAN,
	e.g., smoking status */
	
 /*************************************************************************/
 /*                                                                       */
 /*             Bonus: percentile difference across subgroups             */
 /*                                                                       */
 /*************************************************************************/

/* note: we may be interested in difference or ratio among subgroups. For example,
	we might want to look at energy intake difference across the distribution.
	This can be done with some data manipulation as shown below. */
	
/* 1) transpose data to have a long-formatted data */
	proc transpose data=results.descript_energy_wts_p out=distrib_t prefix=GR_;
	var mean_mc_t tpercentile1-tpercentile99 cutprob1-cutprob5;
	id sex;
	run;

/* 2) Some formatting, calculate difference */
	data distrib_t_diff ;
		set distrib_t (rename= (GR__255 = GR_0));
	
	* For percentile values ... ;
	if (index(_NAME_,"percentile")>0) then do;
		* add numerical value for each percentile;
		percentile = input(compress(_NAME_,,'a'),10.);
		
		* difference;
		diff_2_1 = GR_2 - GR_1 ;
	end;
	
	*For cutprobs;
	if (index(_NAME_,"cut")>0) then do;
		* Ratio;
		diff_2_1 = GR_2 / GR_1 ;
	end;
	
	* For mean;
	if _NAME_ ="mean_mc_t" then
		diff_2_1 = GR_2 - GR_1 ;
	run;
	
/* 3) Show brief summary */
	proc print data=distrib_t_diff label;
	title1 "Energy intake difference across distribution percentile between females and males";
	where percentile in (5 25 50 75 95);
	format  GR_1 GR_2 diff_2_1 4.2;
	id percentile ;
	var GR_1 GR_2 diff_2_1 ;
	label GR_1 = "Energy in males, kcal"
		  GR_2 = "Energy in females, kcal"
		  diff_2_1 = "Energy difference in females vs. males, kcal";
	run;
	title1;

	proc print data=distrib_t_diff label;
	title1 "Prevalence ratio of energy < cutpoints between females and males";
	footnote1 italic "Intepretation example: females are 2.4 times more likely to have usual energy intake values below 2000 kcal (cutprob2)";
	where (index(_NAME_,"cut")>0) ;
	format  GR_1 GR_2 percent6.4 diff_2_1 4.2;
	id _NAME_ ;
	var GR_1 GR_2 diff_2_1 ;
	label _NAME_ = 'Cutpoints used in %DISTRIB;'
		  GR_1 = "Pr(X<x) in males"
		  GR_2 = "Pr(X<x) in females"
		  diff_2_1 = "Prevalence ratio";
	run;
	title1;

	
/* 		  cutprob1 = "Pr(X<1500)" */
/* 		  cutprob2 = "Pr(X<2000)" */
/* 		  cutprob3 = "Pr(X<2500)" */
/* 		  cutprob4 = "Pr(X<3000)" */
/* 		  cutprob5 = "Pr(X<3500)" */
	
/* end of code 11 */
		 
