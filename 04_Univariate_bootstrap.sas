 /*************************************************************************/
 /*                                                                       */
 /*                   National Cancer Institute Method                    */
 /*                                                                       */
 /*               Code 14. Univariate, same as 11, ONE PART               */
 /*                                                                       */
 /*                        Author: Didier Brassard                        */
 /*                                                                       */
 /*                               Version 1                               */
 /*                               02NOV2022                               */
 /*                                                                       */
 /*************************************************************************/

/* Indicate file location to work with relative location*/
	
	/* ## TO DO: Update according to your username or file location */
	
	%let path = /home/didier.brassard.10/nci_workshop/ ;

/* Include NCI macros */
	%include "&path./macros/mixtran_macro_v2.21.sas"; /* Univariate method, boxcox + error model + data manipulation */
	%include "&path./macros/distrib_macro_v2.2.sas"; /* Univariate method, distribution estimation*/

/* Define the `data`, `results` and `bootstrap` library */
	options dlcreatedir;
	libname data "&path./data/";
	libname results "&path./results/";
	libname boot "&path./results/bootstraps/";
	
/* warning: variance for survey data should be estimated by looping
	the NCI analysis through the bootstrap replicate weights (not shown).
	However, for demonstration purpose only, we assume the preNCI data is that of
	a simple random sample of the population. Standard error and confidence interval
	obtained in the procedure below are NOT valid for CCHS and only shown
	for demonstration purpose. */
	
 /*************************************************************************/
 /*                                                                       */
 /*           PREPARE AN INPUT DATASET FOR BOOTSTRAP RESAMPLING           */
 /*                                                                       */
 /*************************************************************************/

	/* note: original data create in program <01_Basics.sas>. Show 5 observations. */

	proc print data=data.prenci (obs=5) ;run;
	
/* create a copy for analysis */
	proc sort data=data.preNCI out=preNCI;
	by adm_rno;
	run;

/* 1) we transpose 24-h recall data to `wide` format */

	/* coding note: since we are resampling PARTICIPANT (adm_rno), we have 
	to ensure that 24-hour recalls data appear on the SAME row (i.e., unique row by participant) */


	/* 1.1 keep weekend indicator */
	proc transpose data=preNCI prefix=weekend_r out=_wknd(keep=adm_rno weekend:);
	by adm_rno;
	var weekend ;
	id recallid;
	run;
	
	/* 1.2 dietary variables of interest */
	proc transpose data=preNCI prefix=energy_r out=_var1(keep=adm_rno energy: );
	by adm_rno;
	var energy ;
	id recallid;
	run;

	proc transpose data=preNCI prefix=milk_r out=_var2(keep=adm_rno milk: );
	by adm_rno;
	var milk ;
	id recallid;
	run;
	
	/* 1.3 covariables only */
			
	/* note: here, we keep only covariables that are NOT recall-specific*/
	
	data _covar;
		set preNCI(where=(recallid=1));
	keep adm_rno wts_p age sex young male ;
	run;
	
	/* 1.4 Merge all data in the wide format */
	data preNCI_w;
	merge _covar
		  _wknd 
		  _var1
		  _var2;
	by adm_rno;
	run;
		   
	/* coding tip: at this point, the number of rows in the data
		should equal to the sample size (unique participant) */

 /*************************************************************************/
 /*                                                                       */
 /*                     PERFORM BOOTSTRAP RESAMPLING                      */
 /*                                                                       */
 /*************************************************************************/

	/* technical note: the bootstrap resampling must be consistent with 
		the original sample sampling strategy */

	proc surveyselect data=preNCI_w out=_bootrep seed=1563783 
			method=URS sampsize=1901 outhits reps=30;
	run;

	/* technical note: it is common to use 200+ bootstrap samples.
		For demonstration purpose, only 30 samples are used */

/* append the original sample (no bootstrap) with bootstrap samples */
	data boot.original_and_boot;
	set preNCI_w(in=a) _bootrep ;
	if a then replicate=0; * original sample = replicate 0; 
	run;
	 
	proc sort data=boot.original_and_boot;
	by replicate;
	run;
 
/* add a new participant identifier for each pseudo-subject in bootstrap samples */
	data boot.original_and_boot;
	set boot.original_and_boot;
	by replicate;
	retain replicaterowid;
		if first.replicate then replicaterowid=1;
		else replicaterowid=replicaterowid + 1;
	run;

/* confirm consistent sample size across bootstrap samples */
	proc freq data=boot.original_and_boot;
	title1 "Sample size in original and bootstrap samples (should be the same)";
	table replicate;
	run;
	title1;

 /*************************************************************************/
 /*                                                                       */
 /*                     MACRO FOR BOOTSTRAP ANALYSIS                      */
 /*                                                                       */
 /*************************************************************************/

/* Objective:
	estimate usual energy intake distribution and calculate difference
	between males and females across the distribution */
	
/* coding note: since we must repeat the ENTIRE modelling in each bootstrap sample,
	a macro is used to loop the analysis.*/

%macro boot_univariate_amount(bootlib=,bootdata=,replicfirst=0,repliclast=0,setseed=78112);
%global success;

/* print information to log */
%put &=bootlib;
%put &=bootdata;
%put &=replicfirst;
%put &=repliclast;
%put &=setseed;
%put &sysdate9, &systime ;

 /************************************************/
 /*             Begin bootstrap loop             */
 /************************************************/

%do replicnum = &replicfirst %to &repliclast;
	
	/* show rep number */
	%put # Replicate = &replicnum / &repliclast ;
	
	/* output current replicate data */
	data subj1rec;
		set &bootlib..&bootdata;
	if (replicate=&replicnum) then output;
	run;
	
	proc sort;
		by replicaterowid;
	run;
	
	/* confirm observations in current data >0 */
	data _NULL_;
		if 0 then
			set subj1rec nobs=n;
		call symputx('totobs', n);
		stop;
	run;
	
		/* hardstop if no observations in data */
		%if (&totobs=0) %then %do;
		%put ERROR: &sysmacroname.: there are no observations in <&bootlib..&bootdata> for replicate=&replicnum ;
		%put # Tip: Verify if there are valid observations bootstrap samples with replicate = &replicnum ;
		%return;
		%end;
	
	/************************************************/
	/*          Prepare data for NCI macros         */
	/************************************************/

		/* transpose 24-h recall specific data from WIDE to LONG */
		
			/* note: the NCI method requires a long (or tall) formatted data where
				24-h recall data appears on ROWS */
			
		data datamrec (drop = energy_r1 milk_r1 energy_r2 milk_r2 weekend_r1 weekend_r2 ith) ; * drop recall-specific names ;
			retain replicate replicaterowid recallid;
			set subj1rec;
		* generic var. names ;
			array foodlist0 (*) energy milk ;
		* assign generic var. names to recall specific var. names ;
		  * Recall #1: output data if there are no missing values ;
		if nmiss(of energy_r1 milk_r1 ) = 0 then do;
			* change recall-specific to generic variable names;
			array foodlistr1 (*) energy_r1 milk_r1 ;
				do ith=1 to dim(foodlist0);
					foodlist0(ith) = foodlistr1(ith);
					* reads as generic name = recall-specific ... ;
				end;
			* sequential recall id ;
				recallid = 1;
			* change recall-specific weekend indicator to generic variable name, if need be ;
				weekend = weekend_r1;
			* writes observation for current recall ;
			output;
		end;
		* Recall #2: output data if there are no missing values ;
		if nmiss(of energy_r2 milk_r2 ) = 0 then do;
			* change recall-specific to generic variable names;
			array foodlistr2 (*) energy_r2 milk_r2 ;
				do ith=1 to dim(foodlist0);
					foodlist0(ith) = foodlistr2(ith);
					* reads as generic name = recall-specific ... ;
				end;
			* sequential recall id ;
				recallid = 2;
			* change recall-specific weekend indicator to generic variable name, if need be ;
			weekend = weekend_r2;
			* writes observation for current recall ;
			output;
		end;
		run;
		
		proc sort data=datamrec;
		by replicaterowid recallid;
		run;
		
		data datamrec;
		set datamrec;
		* add sequence of 24-h recall indicators;
		if recallid = 2 then seq2=1;
		else seq2=0;
		run;
	
	/************************************************/
	/*            Call the MIXTRAN macro            */
	/************************************************/

	%MIXTRAN(data      = datamrec,
		 response      = energy,
		 foodtype      = energy,
		 subject       = replicaterowid, /* bootstrap subject identifier (vs. adm_rno) */
		 repeat        = recallid, 
		 covars_prob   = , 
		 covars_amt    = young male,
		 outlib        = work, 
		 modeltype     = AMOUNT,
		 lambda        = , 
		 replicate_var = WTS_P, 
		 seq           = seq2 ,
		 weekend       = weekend,
		 titles=2, printlevel=1);

	/* successful execution of MIXTRAN check */
	%if (&Success=0) %then %goto next_rep;
	
	
	/************************************************/
	/*         Save output data from MIXTRAN        */
	/************************************************/
	/* predicted value data */
	data &bootlib.._pred_energy&replicfirst._to_&repliclast;
	retain replicate;
	set %if (&replicnum=&replicfirst) %then %do;
		_pred_unc_energy(in=a);
		%end;
		%else %if (&replicnum>&replicfirst) %then %do;
		&bootlib.._pred_energy&replicfirst._to_&repliclast
		_pred_unc_energy(in=a)
		%end;
		;
	if a then replicate=&replicnum;
	run;

	/* param data */
	data &bootlib.._param_energy&replicfirst._to_&repliclast;
	retain replicate;
	set %if (&replicnum=&replicfirst) %then %do;
		_param_unc_energy(in=a);
		%end;
		%else %if (&replicnum>&replicfirst) %then %do;
		&bootlib.._param_energy&replicfirst._to_&repliclast
		_param_unc_energy(in=a)
		%end;
		;
	if a then replicate=&replicnum;
	run;

	/************************************************/
	/*            Call the DISTRIB macro            */
	/************************************************/
	
	%let loopseed = %eval(&setseed + &replicnum) ;
	
	%DISTRIB(call_type  = FULL,
			 seed       = &loopseed , /* seed updated for each replicate */
			 nsim_mc    = 250, /* the number of pseudo-individuals */
			 mcsimda    = mcsim,
			 modeltype  = AMOUNT,
			 pred       = _pred_unc_energy,
			 param      = _param_unc_energy,
			 outlib     = work,
			 cutpoints  = 1500 2000 2500 3000 3500,
			 ncutpnt    = 5,
			 byvar      = ,
			 add_da     = datamrec, /* same as <data> in MIXTRAN */
			 subgroup   = sex, /* distribution for subgroups used in the error model */
			 subject    = replicaterowid, /* bootstrap subject identifier (vs. adm_rno) */
			 titles     = 0,
			 food       = energy /* used for naming datasets */
			 );
			 
	/* calculate difference for percentile values */
	proc transpose data=descript_energy_wts_p out=distrib_t prefix=GR_;
	var mean_mc_t tpercentile1-tpercentile99 cutprob1-cutprob5;
	id sex;
	run;

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
	
	/* transpose back */
	proc transpose data=distrib_t_diff out=distrib_w_diff;
	var diff_2_1;
	run;
	
	/* add to percentile estimate data */
	data descript_energy_wts_p;
	set descript_energy_wts_p distrib_w_diff(in=b drop=_NAME_) ;
	if b then sex=9992;
	label cutprob1 = "Pr(X<1500)"
		  cutprob2 = "Pr(X<2000)"
		  cutprob3 = "Pr(X<2500)"
		  cutprob4 = "Pr(X<3000)"
		  cutprob5 = "Pr(X<3500)";
	run;
	
	/* clean temporary */
	proc delete lib=work data= distrib_t distrib_t_diff distrib_w_diff ; run;

	/************************************************/
	/*         Save output data from DISTRIB        */
	/************************************************/
	
	/* distribution data */
	data &bootlib.._distrib_energy&replicfirst._to_&repliclast;
	retain replicate;
	set %if (&replicnum=&replicfirst) %then %do;
		descript_energy_wts_p(in=a);
		%end;
		%else %if (&replicnum>&replicfirst) %then %do;
		&bootlib.._distrib_energy&replicfirst._to_&repliclast
		descript_energy_wts_p(in=a)
		%end;
		;
	if a then replicate=&replicnum;
	run;


%next_rep:

	/* clean temporary data to avoid carry-over */
	proc delete lib=work data=
		datamrec subj1rec
		_param_unc_energy _pred_unc_energy etas_energy _parmsf2_energy
		mcsim descript_energy_wts_p ; 
	run;

 /************************************************/
 /*            End of bootstrap loop             */
 /************************************************/

%end; /* end of the bootstrap loop */

/* end of macro */
%mend boot_univariate_amount;

 /*************************************************************************/
 /*                                                                       */
 /*                      PERFORM BOOTSTRAP ANALYSIS                       */
 /*                                                                       */
 /*************************************************************************/

/* Confirm that macro works as intended on original sample (repicate=0) */

%boot_univariate_amount(bootlib=boot,
				      bootdata=original_and_boot,
				      replicfirst=0,
				      repliclast=0,
				      setseed=78112);

/* Run the analysis for 30 bootstrap samples */
	
	/* tip: it is useful to save the log in a .log file rather than printing */

	FILENAME UILOG "&path./results/bootstraps/log_univariate0_to_30_&sysdate9..log";
	PROC PRINTTO LOG = UILOG; run;
		
	%boot_univariate_amount(bootlib=boot,
					      bootdata=original_and_boot,
					      replicfirst=0,
					      repliclast=30,
					      setseed=78112);
	proc printto;run;
	
	/* note: analysis require approx. 4 minutes */

 /*************************************************************************/
 /*                                                                       */
 /*                          VARIANCE ESTIMATION                          */
 /*                                                                       */
 /*************************************************************************/

/* note: there are many different methods to estimate variance with the bootstrap
	method. For demonstration purpose, the normal approximation is shown */
	
	
/* 1) select some statistics, split original and bootstrap estimates */
	data original boot;
	set boot._distrib_energy0_to_30 ;
	if replicate=0 then output original;
	else output boot;
	keep replicate sex mean_mc_t tpercentile1 tpercentile5 tpercentile10 tpercentile25 tpercentile50 
			tpercentile75 tpercentile90 tpercentile95 tpercentile99 cutprob1-cutprob5 ;
	run;

	proc sort data=boot;
	by sex;
	run;
	 
	proc sort data=original;
	by sex;
	run;
 
/* 2) Count number of non-missing values (i.e., successful bootstrap estimates) */
	proc means data=boot nmiss noprint;
	by sex ;
	var mean_mc_t tpercentile1 tpercentile5 tpercentile10 tpercentile25 tpercentile50 tpercentile75 
	 tpercentile90 tpercentile95 tpercentile99 cutprob1-cutprob5 ;
	output out=_nboot(drop=_TYPE_) nmiss=;
	run;
	 
	data _nboot;
	set _nboot;
	array varlist (*) mean_mc_t tpercentile1 tpercentile5 tpercentile10 tpercentile25 tpercentile50 
	 tpercentile75 tpercentile90 tpercentile95 tpercentile99 cutprob1-cutprob5 ;
	* Number of non-missing bootstrap = total obs. (_FREQ_) - nmiss ;
	do i=1 to dim(varlist);
	varlist(i) = _FREQ_ - (varlist(i));
	end;
	run;
	
	/* 2.1 Transpose for merging purpose */
	proc transpose data=_nboot out=_nboott(rename=(COL1=nboot)) prefix=COL ;
	by sex ;
	var mean_mc_t tpercentile1 tpercentile5 tpercentile10 tpercentile25 tpercentile50 tpercentile75 
	 tpercentile90 tpercentile95 tpercentile99 cutprob1-cutprob5 ;
	run;
 
/* 3) calculate Standard Errors (i.e., Standad Deviation of the sampling distribution) */
	proc means data=boot stddev noprint;
	by sex ;
	var mean_mc_t tpercentile1 tpercentile5 tpercentile10 tpercentile25 tpercentile50 tpercentile75 
	 tpercentile90 tpercentile95 tpercentile99 cutprob1-cutprob5 ;
	output out=_std(drop=_TYPE_) std=;
	run;

	/* 3.1 Transpose for merging purpose */
	proc transpose data=_std out=_stdt(rename=(COL1=se)) prefix=COL ;
	by sex ;
	var mean_mc_t tpercentile1 tpercentile5 tpercentile10 tpercentile25 tpercentile50 tpercentile75 
	 tpercentile90 tpercentile95 tpercentile99 cutprob1-cutprob5 ;
	run;
 
/* 4) transpose estimates for merging purpose */
	proc transpose data=original out=_estimate (rename=(estimate1=estimate)) prefix=estimate;
	by sex ;
	var mean_mc_t tpercentile1 tpercentile5 tpercentile10 tpercentile25 tpercentile50 tpercentile75 
	 tpercentile90 tpercentile95 tpercentile99 cutprob1-cutprob5 ;
	run;

/* 5) merge estimate, standard errors and n */
	data results.distrib_energyf (rename=(_NAME_=name _LABEL_ = label));
	merge _estimate _stdt _nboott;
	by sex ;
	* calculate 95ci, pvalue, coefficient of variation;
	alpha = 1 - (95/100);
	one_minus_half_alpha = 1 - alpha/2;
	if nboot-2 >0 then do;
		t_quant = quantile('t', one_minus_half_alpha, nboot-2);
		lcl95 = estimate - t_quant * se;
		ucl95 = estimate + t_quant * se;
		if (se ne 0) then do;
			tvalue =abs( ( estimate - 0 ) / se );
			format pvalue PVALUE6.4;
			pvalue = 2 * (1 - probt(tvalue, nboot-2 ) );
		end;
		else do;
			tvalue = .;
			pvalue = .;
		end;
	end;
	else do;
		t_quant=.;
		lcl95=.;
		ucl95=.;
		tvalue = .;
		pvalue = .;
	end;
	
	label pvalue = "2-sided pvalue (null H value=0)";
	 
	* coefficient of variation ;
	format cv PERCENT6.1;
	if not missing(estimate) & (estimate ne 0) then cv = abs(se / estimate);
		else cv = .;
	 
	* clear ugly labels;
	label _NAME_= " " _LABEL_ = " " nboot = " ";
	run;

 /*************************************************************************/
 /*                                                                       */
 /*                       PRINT SUMMARY OF RESULTS                        */
 /*                                                                       */
 /*************************************************************************/

/* warning: remember that the bootstrap variance estimation assumed
	that the sample was a simple random sample of the population. However,
	we know this assumption is incorrect, because respondents are from 
	the Canadian Community Health Survey 2015. Therefore, variance estimation
	(standard errors, 95%CI, pvalue) is incorrect. A proper variance estimation
	would have involved weighting each sample with the survey bootstrap replicate weights */

/* Apply formatting for clarity */
	proc format;
	value sexfmt
	-255 = "Males and females, 19-30 y"
	1 = "Males, 19-30 y"
	2 = "Females, 19-30 y"
	9992 = "Diff. females vs. males";
	run;

/* Show estimates */
	title1 'Usual energy intake distribution in 19-30 y (SE and 95%CI assumes SRS)';
	proc print data=results.distrib_energyf;
	by sex;
	format sex sexfmt. estimate se lcl95 ucl95 5.2 cv percent6.4;
	id sex name label nboot;
	var estimate se lcl95 ucl95 cv;
	run;
	title1;

/* Make distribution graph */
	data for_sgplot;
	set results.distrib_energyf;
	*make numerical percentile index;
	if index(name,"percentile")>0 then percentile = input(compress(name,,'a'),10.);
	* keep percentile only;
	if index(name,"percentile")>0 then output;
	run;
	
	
	proc sgplot data=for_sgplot ;
		title1 "Cumulative energy intake distribution (19-30 y, CCHS-2015)";
		footnote1 italic justify=right "NCI univariate method";
		series x=estimate y=percentile / legendlabel="Estimate" lineattrs=(color=black thickness=3);
		series x=lcl95 y=percentile / legendlabel='LCL' lineattrs=(color=black pattern=LONGDASH thickness=2.5);
		series x=ucl95 y=percentile / legendlabel='UCL' lineattrs=(color=black pattern=LONGDASH thickness=2.5);
	where sex=-255; *all;
	run;

/* end of code 14 */
