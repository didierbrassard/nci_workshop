 /*************************************************************************/
 /*                                                                       */
 /*            Code 19. Plausibility of reported energy intakes           */
 /*                                                                       */
 /*                      Huang and McCrory approach                       */
 /*                                                                       */
 /*                        Author: Didier Brassard                        */
 /*                                                                       */
 /*                               Version 1                               */
 /*                               24OCT2024                               */
 /*                                                                       */
 /*************************************************************************/

/* Indicate file location to work with relative location*/
	
	/* ## TO DO: Update according to your username or file location */
	
	%let path = /home/didier.brassard.10/nci_workshop/ ;

/* Include NCI macros */
	%include "&path./macros/nlmixed_univariate_macro_v1.2.sas"; /* Univariate method, error model */
	%include "&path./macros/iom_eer.sas"; /* Univariate method, error model */
	
/* Define the `data`, `results` and `bootstrap` library */
	options dlcreatedir;
	libname data "&path./data/";
	libname results "&path./results/";
	
/* warning: variance for survey data should be estimated via 
	bootstrap replicate weights (not shown). However, for demonstration purpose only,
	we assume the preNCI data is that of a simple random sample of the population. 
	Standard error and confidence interval obtained in the procedure below are
	NOT valid for CCHS and only shown for demonstration purpose. */

 /*************************************************************************/
 /*                                                                       */
 /*                         Prepare required data                         */
 /*                                                                       */
 /*************************************************************************/

/* 1) import preprocessed HS_NCI data using <proc import> */
	PROC IMPORT DATAFILE="&path./data/hsnci_19to30y_hwbw.xlsx"
		DBMS=XLSX
		OUT=WORK.hsnci_hwbw  replace;
		GETNAMES=YES;
	RUN;
	
	proc sort data=hsnci_hwbw;
		by adm_rno;
	run;

/* 2) Additional data manipulation needed for analysis, save */

	/* note: NCI macros always require that recall-specific dietary intakes 
		appear on unique rows, i.e., the data has m records per participant
		where m = max number of 24hr */

	data data.plausible;
		set hsnci_hwbw (keep= adm_rno age sex bodyweight_kg height_m recallid energy);
	* create indicator variable for sequence of 24-h recall;
	seq2 = 0;
	if recallid=2 then seq2 = 1;
	
	* create indicator for age group (<25 y or >=25 y);
	if not missing(age) then do;
		young = 0;
		if age <25 then young =1 ;
	end;
	
	* create indicator for sex ;
	if not missing(sex) then do;
		female=0;
		if sex=2 then female=1;
	end;
	
	* delete participant with missing data;
	if missing(height_m) or missing(bodyweight_kg) then delete;
	
	run;
	
	/* 2.1 confirm data recoding */
	title1 "Prelim. Confirm recoding of covariates";
		proc freq data=data.plausible;
		table recallid * seq2
			  sex * female
			  age * young
			  ;
		run;
	title1;
	
/* 3) Calculate mean energy (i.e. based on all 24-h dietary recall) */

	/* 3.1) create a sorted copy of data */
	proc sort data=data.plausible out=plausible_long;
	by adm_rno;
	run;
	
	title1 "Plausible data - long format";
	proc print data=plausible_long(obs=5); run;
	title1; 
	
	/* 3.2) transpose for long/tall to wide format */
	proc transpose data=plausible_long(keep=adm_rno recallid energy) out=_energy_wide(drop=_:) prefix=energy_r;
	by adm_rno ;
	variable energy ;
	run;
	
	title1 "Plausible data (energy only) - wide format";
	proc print data=_energy_wide(obs=5); run;
	title1;
	
	/* 3.3) generate derived variables needed for calculation of CV */
	data energy_wide;
	set _energy_wide;
	* count number of recall as max. of 2 minus missing values - need to update if more than 2 recalls;
	n_recall = 2- nmiss(energy_r1, energy_r2);
	* calculate mean energy based on recall 1 and 2;
	energy_mean = mean(energy_r1, energy_r2);
	run;
	
	/* delete temporary data*/
	proc datasets lib=work nolist nodetails;
	delete hsnci_hwbw _energy_wide plausible_long; 
	run;
	
	/* note: not mandatory, but can be relevant when dealing with very large data */
	

 /*************************************************************************/
 /*                                                                       */
 /*            Calculate CV - Reported Energy Intake (CV rEI)             */
 /*                                                                       */
 /*************************************************************************/


	/* Technical note: 
		the squared root of the variance is the Standard Deviation (SD) which, when
		experessed on the log scale, is the coefficient of variation  */
		
	/* Reference:
		Cole TJ, Altman DG. Statistics Notes: Percentage differences, symmetry,
	    and natural logarithms. BMJ. 2017;358:j3683. */

	%NLMIxed_Univariate (data           = data.plausible,
	                     subject        = adm_rno, /* Unique subject id */
	                     repeat         = recallid,
	                     response       = energy,
	                     modeltype      = ONEPART,
	                     covars_prob    = ,
	                     covars_amt     = seq2 ,
	                     lambda         = 0, /* 0=log transformation - neeeded here */
                     	 replicate_var  = , /* balancing or sampling weight would go here */
	                     print          = Y,
	                     ntitle         = 2
	                     );
	                     
	/* note: we can ignore the warning because we decided to use 0 as lambda */

	data CV_rEI_all;
		set parms_u;
	sqrt_VAR_E = sqrt(VAR_E);
	CV_rEI = sqrt_VAR_E*100;
	run;
	
	title1 "Within-individual coefficient of variation (CV) in reported energy intake (EI)" ;
	proc print data=CV_rEI_all;
	format VAR_E sqrt_VAR_E CV_rEI num4.2;
	var VAR_E sqrt_VAR_E CV_rEI;
	run;
	title1;
	
	/* interpretation: we see that in this sample, reported energy intake varies
		approx. 39% across 24-h recalls on average */
		
 /************************************************/
 /*        Repeat the analysis for each sex      */
 /************************************************/

/* objective: have one estimate of CV for each sex */

/* check sample size for each */
proc freq data=data.plausible; 
	table sex;
run;

/* Create one data for each */
	data sex1 sex2;
		set data.plausible ;
	if sex=1 then output sex1;
		else output sex2;
	run; 

/* Run univariate measurement error model for one sex */
%NLMIxed_Univariate (data           = sex1,
                     subject        = adm_rno, /* Unique subject id */
                     repeat         = recallid,
                     response       = energy,
                     modeltype      = ONEPART,
                     covars_prob    = ,
                     covars_amt     = seq2 ,
                     lambda         = 0, /* 0=log transformation - neeeded here */
                 	 replicate_var  = , /* balancing or sampling weight would go here */
                     print          = N,
                     ntitle         = 2
                     );
  
/* save output data */
	data CV_rEI_sex1; set parms_u; run; 
                     
/* Run univariate measurement error model for other sex */
%NLMIxed_Univariate (data           = sex2,
                     subject        = adm_rno, /* Unique subject id */
                     repeat         = recallid,
                     response       = energy,
                     modeltype      = ONEPART,
                     covars_prob    = ,
                     covars_amt     = seq2 ,
                     lambda         = 0, /* 0=log transformation - neeeded here */
                 	 replicate_var  = , /* balancing or sampling weight would go here */
                     print          = N,
                     ntitle         = 2
                     );

/* save output data */
	data CV_rEI_sex2; set parms_u; run;
	

 /************************************************/
 /* Append all energy CV, keep relevant columns  */
 /************************************************/

	data cv_rEI;
		set cv_rEI_all(in=all)
		CV_rEI_sex1 (in=sex1)
		CV_rEI_sex2 (in=sex2);
	* calculate CV altogether; 
	CV_rEI = sqrt(VAR_E); 
	*add sex identifier; 
	if all then sex = 0;
	else if sex1 then sex=1;
	else if sex2 then sex=2;
	keep sex CV_rEI ;
	run;
	
	proc print data=cv_rEI;run;

 /*************************************************************************/
 /*                                                                       */
 /*    Calculate Predicted Energy Requirements based on IOM equations     */
 /*                                                                       */
 /*************************************************************************/

/* objective: need a CV with 1 row per participant, i.e., no repeated data */

proc sort data=data.plausible nodupkey out=plausible_1row(drop=recallid seq2 energy);
	by adm_rno;
run;

	/* note: <nodupkey> keeps only the first record of each <by> variable */

/* use the iom_eer macro to calculate pER (see ./macros/iom_eer.sas for details) */
	%iom_eer(
	indat        = plausible_1row,
	 outdat      = plausible_1row_eer,
	 female_flag = female,
	 age         = age,
	 weight      = bodyweight_kg,
	 height      = height_m,
	 pal         = 1 /* note: assuming everyone sedentary, implications?  */
	 ); 


/* calculate CV */
	proc means data=plausible_1row_eer mean std ; 
		class sex; 
		var EER;
		output out=cv_pER mean=pER_mean std=pER_sd;
	run;
	
/* additional formatting for the output data */
	data cv_pER;
		set cv_pER;
	* indicate value for sex= all ; 
	if sex = . then sex =0 ;
	pER_cv = pER_sd / pER_mean;
	rename  _FREQ_ = sample_size; 
	drop _TYPE_ ;
	run;
	

 /*************************************************************************/
 /*                                                                       */
 /*       Calculate average number of 24-h dietary recall completed       */
 /*                                                                       */
 /*************************************************************************/

/* combine average energy data (including n_recall) with plausibility + EER */
	data plausible_1row_eer_energy ; 
		merge energy_wide plausible_1row_eer;
		by adm_rno;
	run;

/* calculate mean number of 24-h dietary recall completed  */
	proc means data=plausible_1row_eer_energy mean ; 
		class sex; 
		var n_recall;
		output out=n_recall_mean mean=n_recall_mean ;
	run;

	/* additional formatting for the output data */
		data n_recall_mean ;
		set n_recall_mean;
			if sex=. then sex=0;
			rename _FREQ_ = sample_size ;
			drop _TYPE_; 
		run;

 /*************************************************************************/
 /*                                                                       */
 /*    Combine all data to calculate uncertainty (confidence limits)      */
 /*                                                                       */
 /*************************************************************************/

data confidence_limit ;
	merge n_recall_mean cv_rEI cv_pER(rename=pER_CV = CV_pER) ;
	by sex;
* estimated normal biological variation in energy requirements (Black and Cole. EJCN 2000); 
	CV_mTEE = 8.2/100 ;
* confidence limits ;
	SD1    = 1*sqrt( (CV_rEI**2) / n_recall_mean + CV_pER**2 + CV_mTEE**2);
	SD1_5  = 1.5*sqrt( (CV_rEI**2) / n_recall_mean + CV_pER**2 + CV_mTEE**2);
	SD2    = 2*sqrt( (CV_rEI**2) / n_recall_mean + CV_pER**2 + CV_mTEE**2);
run;

title1 "Confidence limits for estimating plausibility of reported energy intakes, adults 19-30 y, CCHS-2015";
footnote1 italic "No sampling weights applied, does not reflect population value";

/* define a format for sex - better looking output */
	proc format;
	value sexfmt
		0 = "All"
		1 = "Males"
		2 = "Females";
	run;

	proc print data=confidence_limit ;
	id sex; 
	format sex sexfmt. CV: SD1 percent6.2  ; 
		var sample_size n_recall_mean CV: SD1;
	run;
	title1;
	footnote1;

/* use SD1 to calculate the proportion of under-, plausible, over reporter */
	data _NULL_;
	set confidence_limit; 
		if _N_=1 then call symputx("SD1", SD1);
	run;

%put The value of &=SD1;

	/* note: here, we use the overall value (first row) for all participants.
	Alternatively, we could have merged the <confidence_limit> data by sex with <plausible_1row_eer_energy>
	to use sex-specific values */

 /*************************************************************************/
 /*                                                                       */
 /* Categorize participants based on their ratio and the confidence limit */
 /*                                                                       */
 /*************************************************************************/

/* categorize */ 
	data reporting_status;
		set plausible_1row_eer_energy(rename=EER = pER);
	* input limits ;
	SD1 = &SD1 ;
	LCL = 1 - SD1 ;
	UCL = 1 + SD1 ; 
	
	* mean of 24-h dietary recall vs. predicted energy requirements ratio ;
	if pER not in (0 . ) then rEI_pER_ratio = energy_mean / pER ;
	
	* classify participants;
	if not missing(rEI_pER_ratio) then do;
		if rEI_pER_Ratio < LCL then reporting_status =1 ; 
		else if rEI_pER_Ratio > UCL then reporting_status=3 ;
		else reporting_status=2 ;
	end;
	run;
	
	/*note: numerical values are often easier to handle in further analysis,
	but the value could have been text as well (e.g., <under>, <plausible>, <over>) */

 /*************************************************************************/
 /*                                                                       */
 /*               Show results of reporting status analysis               */
 /*                                                                       */
 /*************************************************************************/

/* make format */
proc format;
value reporting_status
	1 = "Under-reporter"
	2 = "Plausible reporter"
	3 = "Over-reporter";
run;

title1 "Plausibility of reported energy intakes, adults 19-30 y, CCHS-2015";
footnote1 italic "No sampling weights applied, does not reflect population value";

proc freq data=reporting_status ; 
format reporting_status reporting_status. ;
table reporting_status ;
run;

title1;
footnote1;

/* end of code 19 */