 /*************************************************************************/
 /*                                                                       */
 /*                             iom_eer macro                             */
 /*                                                                       */
 /* Calculate Estimated Energy Requirements (EER) based on IOM equations  */
 /*                                                                       */
 /* Author: Didier Brassard                                               */
 /* Contact: didierbrassard.dtp@gmail.com                                 */
 /* Version: 1                                                            */
 /* Date: 25OCT2024                                                       */
 /*                                                                       */
 /*************************************************************************/
 /*                                                                       */
 /* indat       = Input data with all needed variables                    */
 /* outdat      = Output data including EER                               */
 /* female_flag = Sex flag (1 for female, 0 for males).                   */
 /* age         = Age of the individual in years                          */
 /* weight      = Weight of the individual in kilograms                   */
 /* height      = Height of the individual in meters                      */
 /* pal         = Physical act. Level 1 to 4 (sedentary to very active)   */
 /*                                                                       */
 /* Note: Applicable only to individuals aged 9 years +                   */
 /*                                                                       */
 /*************************************************************************/
 

%macro iom_eer(indat=, outdat=, female_flag=, age=, weight=, height=, pal=);
	data &outdat(drop=BMI_macro PAL_macro PAL_factor);
		set &indat;
		PAL_macro=&pal;
		if &AGE < 9 then EER=.;
		BMI_macro=&WEIGHT/(&HEIGHT**2);
		if not missing(BMI_macro) then do;
			if 9<=&AGE<=18 then do;
				if BMI_macro < 25 then do;
					if &female_flag=0 then do;
						if PAL_macro=1 then PAL_factor=1.00;
						if PAL_macro=2 then PAL_factor=1.13;
						if PAL_macro=3 then PAL_factor=1.26;
						if PAL_macro=4 then PAL_factor=1.42;
						EER=113.5 - 61.9*&AGE + PAL_factor*(26.7*&WEIGHT + 903*&HEIGHT);
					end;
					if &female_flag=1 then do;
						if PAL_macro=1 then PAL_factor=1.00;
						if PAL_macro=2 then PAL_factor=1.16;
						if PAL_macro=3 then PAL_factor=1.31;
						if PAL_macro=4 then PAL_factor=1.56;
						EER=160.3 - 30.8*&AGE + PAL_factor*(10*&WEIGHT + 934*&HEIGHT);
					end;
				end; * end of bmi<25 conditional; 
				if BMI_macro >= 25 then do;
					if &female_flag=0 then do;
						if PAL_macro=1 then PAL_factor=1.00;
						if PAL_macro=2 then PAL_factor=1.12;
						if PAL_macro=3 then PAL_factor=1.24;
						if PAL_macro=4 then PAL_factor=1.45;
						EER=-114.1 - 50.9*&AGE + PAL_factor*(19.5*&WEIGHT + 1161.4*&HEIGHT);
					end;
					if &female_flag=1 then do;
						if PAL_macro=1 then PAL_factor=1.00;
						if PAL_macro=2 then PAL_factor=1.18;
						if PAL_macro=3 then PAL_factor=1.35;
						if PAL_macro=4 then PAL_factor=1.60;
						EER=389.2 - 41.2*&AGE + PAL_factor*(15*&WEIGHT + 701.6*&HEIGHT);
					end;
				end; * end of bmi 25+ conditional; 
			end; * end of age <=18 conditional ;
			if &AGE > 18 then do;
				if BMI_macro < 25 then do;
					if &female_flag=0 then do;
						if PAL_macro=1 then PAL_factor=1.00;
						if PAL_macro=2 then PAL_factor=1.11;
						if PAL_macro=3 then PAL_factor=1.25;
						if PAL_macro=4 then PAL_factor=1.48;
						EER=661.8 - 9.53*&AGE + PAL_factor*(15.91*&WEIGHT + 539.6*&HEIGHT);
					end;
					if &female_flag=1 then do;
						if PAL_macro=1 then PAL_factor=1.00;
						if PAL_macro=2 then PAL_factor=1.12;
						if PAL_macro=3 then PAL_factor=1.27;
						if PAL_macro=4 then PAL_factor=1.45;
						EER=354.1 - 6.91*&AGE + PAL_factor*(9.36*&WEIGHT + 726*&HEIGHT);
					end;
				end; *  end of bmi<25 conditional; 
				if BMI_macro >= 25 then do;
					if &female_flag=0 then do;
						if PAL_macro=1 then PAL_factor=1.00;
						if PAL_macro=2 then PAL_factor=1.12;
						if PAL_macro=3 then PAL_factor=1.29;
						if PAL_macro=4 then PAL_factor=1.59;
						EER=1085.6 - 10.08*&AGE + PAL_factor*(13.7*&WEIGHT + 416*&HEIGHT);
					end;
					if &female_flag=1 then do;
						if PAL_macro=1 then PAL_factor=1.00;
						if PAL_macro=2 then PAL_factor=1.16;
						if PAL_macro=3 then PAL_factor=1.27;
						if PAL_macro=4 then PAL_factor=1.44;
						EER=447.6 - 7.95*&AGE + PAL_factor*(11.4*&WEIGHT + 619*&HEIGHT);
					end;
				end; * end of bmi 25+ conditional; 
			end; * end of age 18 +condtional  ;
		end; * end of not missing BMI conditional ;
		else EER=.;
		Label EER="Estimated Energy Requirements (IOM eq)" ;
	run;
%mend iom_eer;
