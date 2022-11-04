
/*************************************************************************/
/*************************************************************************/
/*                                                                       */
/*                            INDIVINT MACRO                             */
/*                                                                       */
/*************************************************************************/
/*                       VERSION 2.3               2/1/2018              */
/*                                                                       */
/*                                                                       */
/* The INDIVINT macro calculates predicted values for regression         */
/* calibration using methods from Kipnis et al. (Biometrics, 2009) and   */
/* using results from an amount-only model or a two-part model fit using */
/* the MIXTRAN macro.  The INDIVINT macro performs adaptive Gaussian     */
/* quadrature to predict usual intake for each individual, and the macro */
/* allows the user to provide a Box-Cox transformation parameter in      */
/* order to calculate the predicted values on a transformed scale.  The  */
/* results from this macro are intended for use in a subsequent          */
/* regression model as discussed by Kipnis et al. (Biometrics, 2009).    */
/*                                                                       */
/* The syntax for calling the INDIVINT macro is:                         */
/*                                                                       */
/* %indivint(model12=, subj1recdata=, recid=, r24vars=, min_amt=,        */
/*           var_u1=, var_u2=, cov_u1u2=, var_e=, lambda=, xbeta1=,      */
/*           xbeta2=, boxcox_t_lamt=, lamt=, dencalc=, denopt=,          */
/*           u1nlmix=, u2nlmix=, titles=, notesprt=);                    */
/*                                                                       */
/* where the parameters are described as follows.                        */
/*                                                                       */
/*  "model12"            Specifies the type of model that was fit prior  */
/*                       to calling this macro.  A value of 1 indicates  */
/*                       that an amount-only model was fit, and a value  */
/*                       of 2 indicates that a two-part model was fit    */
/*                       where part 1 is the probability part of the     */
/*                       model and part 2 is the amount part of the      */
/*                       model.                                          */
/*                                                                       */
/*  "subj1recdata"       Specifies a data set with 1 record for each     */
/*                       individual.                                     */
/*                                                                       */
/*  "recid"              Specifies an identification (ID) variable that  */
/*                       uniquely identifies each individual's record.   */
/*                                                                       */
/*  "r24vars"            Specifies the 24-hour recall variables with     */
/*                       values that are either non-negative or a SAS    */
/*                       missing value if the 24-hour recall is not      */
/*                       available.  Variables must be space delimited   */
/*                       as illustrated in the following example:        */
/*                       "r24vars=r24hr1 r24hr2".                        */
/*                       Note for Advanced Users:  If all 24-hour recall */
/*                       values are missing for each subject, then the   */
/*                       denominator integration should not be           */
/*                       performed, so the "dencalc" macro parameter     */
/*                       should be specified as "dencalc=n".             */
/*                                                                       */
/*  "min_amt"            Specifies a variable that provides the minimum  */
/*                       intake amount.  This value may be selected as   */
/*                       the smallest value among the observed           */
/*                       consumption-day amounts.  Note that the         */
/*                       specified variable provides the same value for  */
/*                       each individual.  This value will be divided in */
/*                       half and used in the calculations for the       */
/*                       numerator integration.                          */
/*                                                                       */
/*  "var_u1"             Specifies a variable that provides the variance */
/*                       estimate for u1, the random effect from the     */
/*                       probability part of the model.  If a variable   */
/*                       is specified, then the macro will use its value */
/*                       as a diagonal entry of the covariance matrix    */
/*                       which is either a 1x1 matrix or a 2x2 matrix    */
/*                       depending on the number of random effects that  */
/*                       are in the model.                               */
/*                                                                       */
/*  "var_u2"             Specifies a variable that provides the variance */
/*                       estimate for u2, the random effect from the     */
/*                       amount part of the model or from an amount-only */
/*                       model.  If a variable is specified, then the    */
/*                       macro will use its value as a diagonal entry of */
/*                       the covariance matrix which is either a 1x1     */
/*                       matrix or a 2x2 matrix depending on the number  */
/*                       of random effects that are in the model.        */
/*                                                                       */
/*  "cov_u1u2"           Specifies a variable that provides the estimate */
/*                       of the covariance(u1, u2) from the two-part     */
/*                       model.  If the two-part model was an            */
/*                       uncorrelated model, then the specified variable */
/*                       should have a value of zero for every           */
/*                       individual's record.                            */
/*                                                                       */
/*  "var_e"              Specifies a variable that provides the variance */
/*                       estimate for e, the within-person error term    */
/*                       from the amount part of the model or from an    */
/*                       amount-only model.                              */
/*                                                                       */
/*  "lambda"             Specifies a variable that provides the estimate */
/*                       of the Box-Cox parameter, lambda, from the      */
/*                       amount part of the model or from an amount-only */
/*                       model.                                          */
/*                                                                       */
/*  "xbeta1"             Specifies a variable that provides the linear   */
/*                       predictor values calculated using the           */
/*                       covariates and estimates of the fixed effects   */
/*                       parameters from the probability part of the     */
/*                       model.                                          */
/*                                                                       */
/*  "xbeta2"             Specifies a variable that provides the linear   */
/*                       predictor values calculated using the           */
/*                       covariates and estimates of the fixed effects   */
/*                       parameters from the amount part of the model or */
/*                       from an amount-only model.                      */
/*                                                                       */
/*  "boxcox_t_lamt"      If "boxcox_t_lamt=y" or "boxcox_t_lamt=Y" then  */
/*                       individual usual intake will be predicted on a  */
/*                       transformed scale where the Box-Cox             */
/*                       transformation is used with the Box-Cox         */
/*                       parameter value provided by the "lamt" macro    */
/*                       parameter.  The default value for               */
/*                       "boxcox_t_lamt" is "n".                         */
/*                                                                       */
/*  "lamt"               Specifies a variable that provides the Box-Cox  */
/*                       parameter value when "boxcox_t_lamt=y" or       */
/*                       "boxcox_t_lamt=Y".  The macro does not allow    */
/*                       the Box-Cox parameter to be negative.           */
/*                                                                       */
/*  "dencalc"            By default, "dencalc=y" so the denominator      */
/*                       integration is performed.                       */
/*                       Note for Advanced Users:  If all 24-hour recall */
/*                       variables are missing for each subject, then    */
/*                       the denominator integration should not be       */
/*                       performed, so the "dencalc" option should be    */
/*                       specified as "dencalc=n".                       */
/*                                                                       */
/*  "denopt"             By default, "denopt=y" so the denominator       */
/*                       optimization is performed as part of the        */
/*                       denominator integration calculations.           */
/*                       Note for Advanced Users:  In some situations    */
/*                       the denominator optimization is redundant       */
/*                       because the empirical Bayes estimates of u1 and */
/*                       u2 are available from the model fitting         */
/*                       software; therefore, in these situations,       */
/*                       setting "denopt=n" or "denopt=N" allows the     */
/*                       macro to skip this optimization step and use    */
/*                       the variables provided by the "u1nlmix" and     */
/*                       "u2nlmix" macro parameters.                     */
/*                                                                       */
/*  "u1nlmix"            Specifies a variable for an Advanced Users      */
/*                       option.  For details, see the description for   */
/*                       the "denopt" macro parameter.                   */
/*                                                                       */
/*  "u2nlmix"            Specifies a variable for an Advanced Users      */
/*                       option.  For details, see the description for   */
/*                       the "denopt" macro parameter.                   */
/*                                                                       */
/*  "titles"             Specifies the number of title lines to be       */
/*                       reserved for the user's titles.  One additional */
/*                       title line is used by the macro.  The default   */
/*                       value is "0".                                   */
/*                                                                       */
/*  "notesprt"           If "notesprt=n" or "notesprt=N" then notes are  */
/*                       not printed to the SAS log.  The default value  */
/*                       for "notesprt" is "y".                          */
/*                                                                       */
/*************************************************************************/


********************************************************************************;
********************************************************************************;
********************************************************************************;


%macro indivint(model12=, subj1recdata=, recid=, r24vars=, min_amt=, 
                var_u1=, var_u2=, cov_u1u2=, var_e=, lambda=, xbeta1=, xbeta2=, 
                boxcox_t_lamt=n, lamt=, dencalc=y, denopt=y, u1nlmix=, u2nlmix=, 
                titles=0, notesprt=y);


  ********************************************************************************;
  ***  Check user input and initialize macro options and variables  **************;
  ********************************************************************************;


  %if (&model12^=1) and (&model12^=2) %then %put Error:  Unknown value for model12.;
  
  %if (&subj1recdata = %str()) %then %put Error:  Specify value for subj1recdata.;
  
  %if (&recid = %str()) %then %put Error:  Specify value for recid.;

  %if (&min_amt = %str()) %then %put Error:  Specify value for min_amt.;

  %if (&var_e = %str()) %then %put Error:  Specify value for var_e.;
  
  %if (&lambda = %str()) %then %put Error:  Specify value for lambda.;
  
  %if (&model12=1) %then %do;
    %if (&xbeta1 ^= %str()) %then %put Error: Amount-only model specified, so xbeta1 is not applicable.;
    %if (&var_u1 ^= %str()) %then %put Error: Amount-only model specified, so var_u1 is not applicable.;
    %if (&cov_u1u2 ^= %str()) %then %put Error: Amount-only model specified, so cov_u1u2 is not applicable.;
    %if (&u1nlmix ^= %str()) %then %put Error: Amount-only model specified, so u1nlmix is not applicable.;
  %end;

  %if (&model12=2) and (&xbeta1 = %str()) %then %put Error:  Specify value for xbeta1.;
  
  %if (&xbeta2 = %str()) %then %put Error:  Specify value for xbeta2.;
  

  %if %upcase(&boxcox_t_lamt)=Y %then %do;
    %if (&lamt = %str()) %then %put Error:  Box-Cox(t,lamt) transformation selected, so lamt is required.;
  %end;
  %else %do;
    %if (&lamt ^= %str()) %then %put Note:  Box-Cox(t,lamt) transformation not selected, so lamt is not applicable.;
  %end;


  %if %upcase(&dencalc)^=Y and %upcase(&dencalc)^=N %then %put Error:  Unknown value for dencalc.;

  %if %upcase(&denopt)^=Y and %upcase(&denopt)^=N %then %put Error:  Unknown value for denopt.;


  %if %upcase(&notesprt)=Y %then %do;
    options notes;
  %end;
  %else %if %upcase(&notesprt)=N %then %do;
    options nonotes;
  %end;
  %else %put Error:  Unknown value for notesprt.;




  %if (&var_u1 = %str()) %then %do;
    
    %let uprob = 0;
    
    %if (&var_u2 = %str()) %then %do;
      %let uamt = 0;
      %put Error:  Variance of the random effect(s) was not provided.;
    %end;
    %else %do;
      %let uamt = u[1];  *** in this situation, the u vector reduces to a scalar, so u = u[1] ***;
      %let randeff = amtonly;

      %if (&cov_u1u2 ^= %str()) %then %put Error:  Only one random effect, so cov_u1u2 is not applicable.;
      
      %if %upcase(&dencalc)=N and %upcase(&denopt=Y) %then %put Error:  Denominator integration not selected, so denominator optimization can not be selected.;
      %else %if %upcase(&dencalc)=Y %then %do;

        %if (%upcase(&denopt)=N) and (&u2nlmix = %str()) %then %put Error:  Denominator optimization not selected, so u2nlmix variable is required.;
        %else %if (%upcase(&denopt)=Y) and (&u2nlmix ^= %str()) %then %put Error:  Denominator optimization selected, so u2nlmix is not applicable.;

      %end;

    %end;
    
  %end;
  %else %do;
    
    %let uprob = u[1];
    
    %if (&var_u2 = %str()) %then %do;
      %let uamt = 0;  *** in this situation, the u vector reduces to a scalar, so u = u[1] ***;
      %let randeff = probonly;

      %if (&cov_u1u2 ^= %str()) %then %put Error:  Only one random effect, so cov_u1u2 is not applicable.;
      
      %if %upcase(&dencalc)=N and %upcase(&denopt=Y) %then %put Error:  Denominator integration not selected, so denominator optimization can not be selected.;
      %else %if %upcase(&dencalc)=Y %then %do;

        %if (%upcase(&denopt)=N) and (&u1nlmix = %str()) %then %put Error:  Denominator optimization not selected, so u1nlmix variable is required.;
        %else %if (%upcase(&denopt)=Y) and (&u1nlmix ^= %str()) %then %put Error:  Denominator optimization selected, so u1nlmix is not applicable.;
        
      %end;

    %end;
    %else %do;
      %let uamt = u[2];   *** in this situation, u = (u[1], u[2]) ***;
      %let randeff = probamt;
      
      %if %upcase(&dencalc)=N and %upcase(&denopt=Y) %then %put Error:  Denominator integration not selected, so denominator optimization can not be selected.;
      %else %if %upcase(&dencalc)=Y %then %do;

        %if (%upcase(&denopt)=N) and (&u1nlmix = %str()) %then %put Error:  Denominator optimization not selected, so u1nlmix variable is required.;
        %else %if (%upcase(&denopt)=Y) and (&u1nlmix ^= %str()) %then %put Error:  Denominator optimization selected, so u1nlmix is not applicable.;
        
        %if (%upcase(&denopt)=N) and (&u2nlmix = %str()) %then %put Error:  Denominator optimization not selected, so u2nlmix variable is required.;
        %else %if (%upcase(&denopt)=Y) and (&u2nlmix ^= %str()) %then %put Error:  Denominator optimization selected, so u2nlmix is not applicable.;

      %end;

    %end;
    
  %end;







  data _null_;
    set &subj1recdata(obs=1);

    %if (&lamt ^= %str()) %then %do;
      if &lamt < 0 then put "Error:  Macro does not allow the lamt variable to have a negative value.";
    %end;

    %if &r24vars=%str() %then %do;
      numr24=0;
    %end;
    %else %do;
      array rvec(*)  &r24vars;
      numr24 = dim(rvec);
    %end;
  
  
    call symput("numr24vars",trim(left(put(numr24, 3.))));
  
  run;

  %put Note:  Number of 24-hour recall variables is &numr24vars;




  *******************************************************************************************;
  *******************************************************************************************;
  
  
  proc format;
    value r24errzf
      0         = "Error: 0 Invalid"
      0< - high = "GT Zero"
      ;
    value r24gtzf
      0< - high = "GT Zero"
      ;
  run;
  
  
  data _subj1rec(drop = 
                    %if &r24vars^=%str() %then %str( kr24vars ); 
                    zero24message);
    set &subj1recdata end=lastrec;
    by &recid;
  
    retain zero24message 0;
  
  
    %if %upcase(&boxcox_t_lamt)^=Y %then %do;
      %let lamt=lamt;
      &lamt=1;  *** No transformation.  Therefore, lamt must be 1 because T##lamt/lamt is used. ***;
    %end;
  
  
  
    %if &r24vars^=%str() %then %do;
  
      array rvec(&numr24vars)    &r24vars;
      array yvec(&numr24vars)    y1-y&numr24vars;   
     
      do kr24vars = 1 to &numr24vars;

        yvec(kr24vars) = rvec(kr24vars); 

        %if &model12=1 %then %do;
          if rvec(kr24vars) = 0 then zero24message=1;  *** set flag so zero value error message is printed for amount-only models;
        %end;

      end;
  
    %end;
  

      
  
    if (lastrec) then do;

      %if &model12=1 %then %do;
        if (zero24message = 1) then put "Error:  Amount-only model specified, so zero is not a valid 24-hour recall value.";
      %end;

    end;
  
    output _subj1rec;
  run;


 
  %if &r24vars^=%str() %then %do;

    proc freq data=_subj1rec;
      tables &r24vars / missing;

      
      %if &model12=1 %then %do;
        format &r24vars r24errzf.;
        title%eval(&titles+1) 'Frequency of 24-Hour Recall Variables.  Amount-Only Model Therefore Zero Not Allowed.';
      %end;
      %else %if &model12=2 %then %do;
        format &r24vars r24gtzf.;
        title%eval(&titles+1) 'Frequency of 24-Hour Recall Variables';
      %end;
      
    run;
    title%eval(&titles+1);

  %end;



  ****************************************************************************;
  ****************************************************************************;
  
  
   
    
  proc iml;
  
    use _subj1rec  var{ &recid, 

                        %if &r24vars^=%str() %then %do;
                          %do hrec = 1 %to &numr24vars;
                            y&hrec,
                          %end;
                        %end;
    
                        %if &model12=2 %then %str(&xbeta1,);
                
                        &xbeta2,

                        %if (&randeff = probamt) %then %do;
                          %if (%upcase(&dencalc)=Y) and (%upcase(&denopt)=N) %then %str(&u1nlmix, &u2nlmix,);
                          &var_u1, &var_u2, &cov_u1u2,
                        %end;
                        %else %if (&randeff = probonly) %then %do;
                          %if (%upcase(&dencalc)=Y) and (%upcase(&denopt)=N) %then %str(&u1nlmix,);
                          &var_u1,
                        %end;
                        %else %if (&randeff = amtonly) %then %do;
                          %if (%upcase(&dencalc)=Y) and (%upcase(&denopt)=N) %then %str(&u2nlmix,);
                          &var_u2,
                        %end;
                        &var_e, &lambda, &lamt, &min_amt
                        };
                





    start gbc_inverse(z_vector, bc_lambda_scalar, minamt_scalar);

      /*** Computes the inverse of the Box-Cox transformation and ensures result >= minimum amount ***/

      if bc_lambda_scalar = 0 then do;
        notneg_gbc_inv = exp(z_vector);
      end;
      else do;
        temp_notneg = (1 + bc_lambda_scalar # z_vector) <> 0;
        notneg_gbc_inv = temp_notneg ## (1 / bc_lambda_scalar);
      end;
      
      gbc_inv = notneg_gbc_inv <> minamt_scalar;

      return (gbc_inv);
        
    finish gbc_inverse;






    start backtransform(n, tran_lambda_j, tran_center_j, tran_scale_j, xbeta_u_j, sigmae_jj, minamt_j);

      /*** When lambda is zero, compute the exact backtransformation for the dietary component ***/
      /*** When lambda is not zero, compute the 9 point backtransformation for the dietary component ***/
      
      /*** Abscissas and weights ***/ 

      x9pt = {-2.1, -1.3, -0.8, -0.5, 0.0, 0.5, 0.8, 1.3, 2.1};
      w9pt_given = {0.063345, 0.080255, 0.070458, 0.159698, 0.252489, 0.159698, 0.070458, 0.080255, 0.063345};
      w9pt = w9pt_given / sum(w9pt_given);


      gstar = j(n, 1, 0);

      if tran_lambda_j = 0 then do;
        notneg_gstar = exp( tran_center_j + tran_scale_j # xbeta_u_j + ((tran_scale_j##2) # sigmae_jj / 2) );
        gstar = notneg_gstar <> minamt_j;
      end;
      else do;
        do qq = 1 to 9;
        
          bc_qq9pt = tran_center_j + tran_scale_j # (xbeta_u_j + x9pt[qq, 1] # (sigmae_jj)##0.5);
          g_inv_qq9pt = gbc_inverse(bc_qq9pt, tran_lambda_j, minamt_j);
          
          gstar = gstar + w9pt[qq, 1] # g_inv_qq9pt;
       
        end;
      end;  

      return (gstar);
        
    finish backtransform;





     
    start mlogintn(u) global(pi, inflate_constant, covmatr, covnrow, &var_e, &lambda, 
                             &lamt, nirep, totcons, 
                             
                             %if &r24vars^=%str() %then %do;
                               %do hrec = 1 %to &numr24vars;
                                 y&hrec, consume&hrec, z&hrec,
                               %end;
                             %end;

                             %if &model12=2 %then %str(&xbeta1,); 
                             &xbeta2, &min_amt);

  
      
      %if &model12=1 %then %do;
      
        logp = 0;  *** amount-only model, so log(p) = log(1) = 0 ***;
      
      %end;
      %else %if &model12=2 %then %do;
  
        xbeta1u1    = &xbeta1+&uprob;
                               
        if xbeta1u1 > 0 then do;
          expneg    = exp(-xbeta1u1);
          if expneg = 0 then expneg=0.000000000001;
          logp      = -log(1+expneg);
          log1_p    = log( expneg/(1+expneg) );
        end;
        else do;
          exppos    = exp(xbeta1u1);
          if exppos = 0 then exppos=0.000000000001;
          logp      = log( exppos/(1+exppos) );
          log1_p    = -log(1+exppos);
        end;
        
      %end;
      

      
      
      nsampsize_in_backtr   = 1;
      tran_center_in_backtr = 0;
      tran_scale_in_backtr  = 1;
      xbeta_u2              = &xbeta2 + &uamt;
      half_min_amt          = 0.5 # &min_amt;
      
      g_star = backtransform(nsampsize_in_backtr, &lambda, tran_center_in_backtr, 
                             tran_scale_in_backtr, xbeta_u2, &var_e, half_min_amt); 


  
      if &lamt=0 then do;
        inflated_prob_amt = inflate_constant # (exp(logp) # g_star);
        gt_temp = log(inflated_prob_amt);
        gt = max(gt_temp, 0.000000000001);
      end;
      else if &lamt>0 then do;
        gt_temp = ((exp(logp) # g_star)##&lamt)/&lamt;
        gt = max(gt_temp, 0.000000000001);
      end;

        
      /* Allow the number of 24-hour recalls to vary between subjects.  Recall observed if yj>=0. */

  
      %if &model12=1 %then %do;

        if nirep ^= totcons then do;
          nocon1_p = .;  *** amount-only model, so only consumption days allowed ***;
        end;
        else nocon1_p = 0;  *** amount-only model ***;

      %end;
      %else %if &model12=2 %then %do;

        nocon1_p = (nirep - totcons) # log1_p;  *** two-part model ***;

      %end;

      
      
      conp = 0;
      %if &r24vars^=%str() %then %do;
        %do krec = 1 %to &numr24vars;
          if y&krec > 0 then conp = conp + ( logp - 0.5#( ( z&krec - &xbeta2 - &uamt )##2 )/&var_e );
        %end;
      %end;


  
      fn      =  - log( gt )
                 - nocon1_p
                 - conp
                 - log( 1 / sqrt( ((2#pi)##covnrow) # det(covmatr)) ) + 0.5#( u * inv(covmatr) * (u`) ) ;
    
      return(fn);
    finish mlogintn;






    start mlogintd(u) global(pi, covmatr, covnrow, &var_e, 
                             nirep, totcons, 
                                                          
                             %if &r24vars^=%str() %then %do;
                               %do hrec = 1 %to &numr24vars;
                                 y&hrec, consume&hrec, z&hrec,
                               %end;
                             %end;
                             
                             %if &model12=2 %then %str(&xbeta1,); 
                             &xbeta2);
  


      %if &model12=1 %then %do;
      
        logp = 0;  *** amount-only model, so log(p) = log(1) = 0 ***;

      %end;
      %else %if &model12=2 %then %do;
  
        xbeta1u1    = &xbeta1+&uprob;
                               
        if xbeta1u1 > 0 then do;
          expneg    = exp(-xbeta1u1);
          if expneg = 0 then expneg=0.000000000001;
          logp      = -log(1+expneg);
          log1_p    = log( expneg/(1+expneg) );
        end;
        else do;
          exppos    = exp(xbeta1u1);
          if exppos = 0 then exppos=0.000000000001;
          logp      = log( exppos/(1+exppos) );
          log1_p    = -log(1+exppos);
        end;
        
      %end;
  
        

      /* Allow the number of 24-hour recalls to vary between subjects.  Recall observed if yj>=0. */

  
      %if &model12=1 %then %do;

        if nirep ^= totcons then do;
          nocon1_p = .;  *** amount-only model, so only consumption days allowed ***;
        end;
        else nocon1_p = 0;  *** amount-only model ***;

      %end;
      %else %if &model12=2 %then %do;

        nocon1_p = (nirep - totcons) # log1_p;  *** two-part model ***;

      %end;


  
      conp = 0;
      %if &r24vars^=%str() %then %do;
        %do krec = 1 %to &numr24vars;
          if y&krec > 0 then conp = conp + ( logp - 0.5#( ( z&krec - &xbeta2 - &uamt )##2 )/&var_e );
        %end;
      %end;


  
      fd      =  - nocon1_p
                 - conp
                 - log( 1 / sqrt( ((2#pi)##covnrow) # det(covmatr)) ) + 0.5#( u * inv(covmatr) * (u`) ) ;
    
      return(fd);
    finish mlogintd;


    
   

  
    options nonotes;
  
  
    *************************************************************************************;
    *************************************************************************************;
  
  
    do data;
     
      read next;
  
      pi = arcos(-1);
      inflate_constant = 10##100;   *** constant to inflate the amount*probability value when lamt=0 ***;

      /* Allow the number of 24-hour recalls to vary between subjects.  Recall observed if yj>=0. */
  
      totcons = 0;
      nirep = 0;

      %if &r24vars^=%str() %then %do;
        %do krec = 1 %to &numr24vars;
          z&krec = .;
          consume&krec = .;
          if y&krec >= 0 then do;

            consume&krec = (y&krec > 0);
            if consume&krec = 1 then do;
              if &lambda = 0 then z&krec = log(y&krec);
              else z&krec = ( y&krec##&lambda - 1 ) / &lambda;
            end;
           
            totcons = totcons + consume&krec;
            nirep = nirep + 1;

          end;
        %end;
      %end;


 
      ************* covariance matrix for random effect(s) ****************;
  
      %if (&randeff = probamt) %then %do;
        row1    = &var_u1 || &cov_u1u2;
        row2    = &cov_u1u2 || &var_u2;
        covmatr = row1 // row2;
      %end;
      %else %if (&randeff = probonly) %then %do;
        covmatr = &var_u1;
      %end;
      %else %if (&randeff = amtonly) %then %do;
        covmatr = &var_u2;
      %end;
  
     
      covnrow =  nrow(covmatr);
  
     
      ***************************************************************************************;
      ************* optimization ************************************************************;
     
      optn = {0 0};
    
      **** optimization for denominator integral ***;
      
      %if %upcase(&dencalc)=Y %then %do;
        %if %upcase(&denopt)=N %then %do;
  
          %if (&randeff = probamt) %then %do;
            ufinald = &u1nlmix || &u2nlmix;                 *** using nlmixed u1hat and u2hat values ***;
          %end;
          %else %if (&randeff = probonly) %then %do;
            ufinald = &u1nlmix;                             *** using nlmixed u1hat values ***;
          %end;
          %else %if (&randeff = amtonly) %then %do;
            ufinald = &u2nlmix;                             *** using nlmixed u2hat values ***;
          %end;
  
          rcd = {99};  *** indicates that nlmixed u1hat and u2hat values are used;
          call nlpfdd(fufinald, gufinald, hufinald, "mlogintd", ufinald);
          ustartn  = ufinald;
          
        %end;
        %else %do;
    
          ustart  = j(1, covnrow, 0.0);
          call nlpqn(rcd, ufinald, "mlogintd", ustart, optn);
          call nlpfdd(fufinald, gufinald, hufinald, "mlogintd", ufinald);
          ustartn  = ufinald;
  
        %end;    
      %end;
      %else %do;
        ustartn  = j(1, covnrow, 0.0);
      %end;
  
      **** optimization for numerator integral ***;
      
      call nlpqn(rcn, ufinaln, "mlogintn", ustartn, optn);
      call nlpfdd(fufinaln, gufinaln, hufinaln, "mlogintn", ufinaln);
   
  
   
        
      ***********************************************************************************************;
      ****  evaluate numerator and denominator integrals and calculate individual usual intake ******;
     
      sumresn = 0;
      %if %upcase(&dencalc)=Y %then %do;   
        sumresd = 0;
      %end;
     
      ****** See Table 25.10 of Abramowitz and Stegun (1972) **************;
  
      
      zn6gh = { - 2.350604973674492   
                - 1.335849074013697  
                - 0.436077411927617  
                  0.436077411927617  
                  1.335849074013697
                  2.350604973674492 };
     
      ezn6sq = exp( zn6gh##2 );
     
      wn6gh = { 0.004530009905509     
                0.1570673203229     
                0.7246295952244    
                0.7246295952244 
                0.1570673203229
                0.004530009905509 };
      
    
  
    
     
      call eigen(eivalsn, eivecsn, hufinaln);
      
      eimhalfn = eivalsn ## (-0.5);
      scalmatn = sqrt(2) # ( eivecsn * diag(eimhalfn) * eivecsn` );
     
  
     
      %if %upcase(&dencalc)=Y %then %do;   
       
        call eigen(eivalsd, eivecsd, hufinald);
      
        eimhalfd = eivalsd ## (-0.5);
        scalmatd = sqrt(2) # ( eivecsd * diag(eimhalfd) * eivecsd` );
     
      %end;
     
     
      do j1 = 1 to 6;
  
        %if &randeff=probamt %then %do;   *** use second loop if two random effects ***;
          do j2 = 1 to 6;
            zvect = zn6gh[j1] // zn6gh[j2];
            wezsq = wn6gh[j1] # ezn6sq[j1] # wn6gh[j2] # ezn6sq[j2];
        %end;
        %else %if (&randeff=probonly) or (&randeff=amtonly) %then %do;   *** only use one loop if one random effect ***;
            zvect = zn6gh[j1];
            wezsq = wn6gh[j1] # ezn6sq[j1];
        %end;
     
            avectn = ufinaln + (scalmatn * zvect)`; 
            sumresn = sumresn + exp(-mlogintn(avectn)) # wezsq;
  
        %if %upcase(&dencalc)=Y %then %do;   
            avectd = ufinald + (scalmatd * zvect)`; 
            sumresd = sumresd + exp(-mlogintd(avectd)) # wezsq;
        %end;
  
        %if &randeff=probamt %then %str( end; );   *** end second loop if two random effects ***;
      end;    
     
  
  
    
      intresn = (2##(covnrow/2)) # (det(hufinaln))##(-0.5) # sumresn;
      
      
      %if (&randeff = probamt) %then %do;
        u1finaln = ufinaln[1];
        u2finaln = ufinaln[2];
      %end;
      %else %if (&randeff = probonly) %then %do;
        u1finaln = ufinaln[1];  *** in this situation, the u vector reduces to a scalar, so u = u[1] ***;
        u2finaln = .;
      %end;
      %else %if (&randeff = amtonly) %then %do;
        u1finaln = .;
        u2finaln = ufinaln[1];  *** in this situation, the u vector reduces to a scalar, so u = u[1] ***;
      %end;
  
         
   
      %if %upcase(&dencalc)=Y %then %do;   
   
        intresd = (2##(covnrow/2)) # (det(hufinald))##(-0.5) # sumresd;
   
        %if (&randeff = probamt) %then %do;
          u1finald = ufinald[1];
          u2finald = ufinald[2];
        %end;
        %else %if (&randeff = probonly) %then %do;
          u1finald = ufinald[1];  *** in this situation, the u vector reduces to a scalar, so u = u[1] ***;
          u2finald = .;
        %end;
        %else %if (&randeff = amtonly) %then %do;
          u1finald = .;
          u2finald = ufinald[1];  *** in this situation, the u vector reduces to a scalar, so u = u[1] ***;
        %end;
   
   

        ******************************************************************;
        ****** E[(T##lamt-1)/lamt|Cond]=E[T##lamt/lamt|Cond]-1/lamt ******;
        ******************************************************************;
        %if %upcase(&boxcox_t_lamt)^=Y %then %do;
          indusint = intresn/intresd;  *** No transformation.  Macro sets lamt to 1 because T##lamt/lamt is used.;
        %end;
        %else %if %upcase(&boxcox_t_lamt)=Y %then %do;
          if &lamt=0 then do;
            indusint = intresn/intresd - log(inflate_constant);
          end;
          else if &lamt>0 then do;
            indusint = intresn/intresd - 1/&lamt;
          end;
        %end;
    
      %end;
      %else %do;
     
        ******************************************************************;
        ****** E[(T##lamt-1)/lamt|Cond]=E[T##lamt/lamt|Cond]-1/lamt ******;
        ******************************************************************;
        %if %upcase(&boxcox_t_lamt)^=Y %then %do;
          indusint = intresn;  *** No transformation.  Macro sets lamt to 1 because T##lamt/lamt is used.;
        %end;
        %else %if %upcase(&boxcox_t_lamt)=Y %then %do;
          if &lamt=0 then do;
            indusint = intresn - log(inflate_constant);
          end;
          else if &lamt>0 then do;
            indusint = intresn - 1/&lamt;
          end;
        %end;
        
      %end;
      
    
  
      resvals   =  &recid || indusint || intresn
                   %if %upcase(&dencalc)=Y %then %str(|| intresd);
                   || nirep || totcons
                   || rcn  || u1finaln  || u2finaln 
                   %if %upcase(&dencalc)=Y %then %str(|| rcd  || u1finald  || u2finald);
                   ;
  
      resmat    =  resmat // resvals ;
  
    end;
  
    *************************************************************************************;
    *************************************************************************************;
    
    %if %upcase(&notesprt)=Y %then %do;
      options notes;
    %end;
  


  
    **print resmat;
  
    varnames = {"&recid" "indusint" "intresn"
                %if %upcase(&dencalc)=Y %then %str("intresd");
                "nirep" "totcons"
                "rcn" "u1finaln" "u2finaln" 
                %if %upcase(&dencalc)=Y %then %str("rcd" "u1finald" "u2finald");
                };
                
    create _resdata from resmat [colname=varnames];
    append from resmat;
    
  quit;



  ****************************************************************************;
  ****************************************************************************;
  
  proc freq data=_resdata;
    tables nirep*totcons / norow nocol nopercent missing;
    title%eval(&titles+1) 'Number of Observed 24-Hour Recalls and Number of 24-Hour Recalls Greater Than 0';
  run;
  title%eval(&titles+1);
  
  
  proc freq data=_resdata;
    tables rcn  %if (%upcase(&dencalc)=Y) and (%upcase(&denopt)=Y) %then rcd;;
    title%eval(&titles+1) 'Return Codes from Optimization Step of Adaptive Gaussian Quadrature.  Positive Value Indicates Successful Termination.';
  run;
  title%eval(&titles+1);
  
  

%mend indivint;  ******** end indivint macro ;


