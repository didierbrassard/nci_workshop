 /****************************************************************************
 *                                                                           *
 * SAS macro NLMixed_Univariate fits a univariate model for a food/nutrient. *
 * The food/nutrient can be episodically consumed or consumed every day.     *
 *                                                                           *
 * Model for episodically consumed foods/nutrients (two-part model):         *
 * For episodically consumed foods/nutrients, the macro fits a two-part      *
 * nonlinear mixed model, where the first part is the probability to         *
 * consume and the second part is the amount consumed on a consumption day.  *
 * The model allows for covariates in each part, includes a random effect    *
 * for each part, and allows the random effects to be correlated.            *
 *                                                                           *
 * Model for foods/nutrients consumed every day (one-part model):            *
 * For foods/nutrients consumed every day, the macro fits a one-part         *
 * nonlinear mixed model of the amount consumed (the probability to consume  *
 * is assumed to be 1). The model allows for covariates and includes a       *
 * random effect.                                                            *
 *                                                                           *
 * For a food/nutrient that is consumed nearly every day by nearly everyone, *
 * so that the number of zero values is small, it may be preferable to use   *
 * the one-part (consumed every day) model, since the two-part model may     *
 * have trouble modeling the probability to consume in such a situation.     *
 *                                                                           *
 * Note, however, that the one-part model requires all responses to be       *
 * greater than zero (zero values are treated as missing values).            *
 * Before fitting the one-part model to a food/nutrient that has some zero   *
 * values, replace the zero values with a small positive value, such as      *
 * half the smallest observed nonzero value.                                 *
 *                                                                           *
 * The macro calls the NLMixed procedure to fit the model.                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * Macro Parameters:                                                         *
 *                                                                           *
 *    Required Parameters:                                                   *
 *       data          = name of SAS data set containing the data to be      *
 *                       analyzed. The data set has multiple observations    *
 *                       for each subject, one for each reptition of the     *
 *                       24-hour recall (or other dietary instrument).       *
 *       subject       = name of the variable that uniquely identifies each  *
 *                       subject (i.e., ID variable).                        *
 *       repeat        = name of the variable that indexes repeated          *
 *                       observations for each subject.                      *
 *       response      = name of the food/nutrient variable to be modeled    *
 *                       (24-hour recall variable for the food/nutrient).    *
 *       modeltype     = model for food/nutrient:                            *
 *                       to fit the two-part (epsisodic) model, specify      *
 *                          modeltype = TWOPART                              *
 *                       to fit the one-part (every day) model, specify      *
 *                          modeltype = ONEPART                              *
 *                                                                           *
 *    Optional Parameters:                                                   *
 *       covars_prob   = list of variables that are covariates in the        *
 *                       probability part of the two-part model.             *
 *                       if modeltype=ONEPART, then covars_prob is ignored.  *
 *       covars_amt    = list of variables that are covariates in the        *
 *                       one-part model or the amount part of the            *
 *                       two-part model.                                     *
 *       link          = link function for the probability part of the two-  *
 *                       part model. to fit a logistic model, specify        *
 *                          link = logit                                     *
 *                       to fit a probit model, specify                      *
 *                          link = probit                                    *
 *                       by default, link = probit.                          *
 *                       if modeltype = ONEPART, then link is ignored.       *
 *       lambda        = Box-Cox transformation parameter for the amount     *
 *                       part of the model. If lambda is not specified,      *
 *                       then it is estimated as part of the model.          *
 *       var_u1        = variance of the random effect in the probability    *
 *                       part of the two-part model.                         *
 *                       If var_u1 is not specified, then it is estimated    *
 *                       as part of the model.                               *
 *                       if modeltype = ONEPART, then var_u1 is ignored.     *
 *       var_u2        = variance of the random effect in the one-part model *
 *                       or the amount part of the two-part model.           *
 *                       If var_u2 is not specified, then it is estimated    *
 *                       as part of the model.                               *
 *       indep_u       = Y if random effects u1 and u2 are independent.      *
 *                     = N if random effects u1 and u2 are dependent.        *
 *                       by default, indep_u = N.                            *
 *                       if modeltype = ONEPART, then indep_u is ignored.    *
 *       replicate_var = name of the sampling weight variable if the data    *
 *                       is from a complex survey with weights.              *
 *                       by default, the macro performs an unweighted        *
 *                       analysis (assumes a simple random sample).          *
 *       nloptions     = options for the NLMixed procedure that are          *
 *                       appended to the PROC NLMIXED statement, e.g.,       *
 *                          nloptions = technique=newrap maxiter=200,        *
 *       init_parms    = name of SAS data set that contains initial          *
 *                       parameter estimates. See the description of output  *
 *                       data set parms_u (below) for further information.   *
 *                       if init_parms is not specified, then the macro      *
 *                       calculates initial parameter estimates.             *
 *       print         = Y to print the output from the model.               *
 *                     = N to supress printing the output from the model.    *
 *                     = V (verbose) to print extra output.                  *
 *                       by default, print = Y.                              *
 *       ntitle        = number of titles defined by the user.               *
 *                       by default, ntitle = 2.                             *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * Output Data Sets:                                                         *
 *                                                                           *
 *   parms_u = data set containing parameter estimates for the model.        *
 *             parms_u contains the following variables:                     *
 *                                                                           *
 *                 A_Intercept = intercept in the amount part of the model.  *
 *                 A_varname   = regression slope for covariate "varname"    *
 *                               in the amount part of the model.            *
 *                 A_LogSDe    = Log(Sqrt(Var_e))                            *
 *                 LogSDu2     = Log(Sqrt(Var_u2))                           *
 *                 Var_e       = variance of the within-person error in the  *
 *                               amount part of the model.                   *
 *                 Var_u2      = variance of the random effect in the        *
 *                               amount part of the model.                   *
 *                                                                           *
 *             if fitting the two-part model, then parms_u also contains     *
 *             the following variables:                                      *
 *                                                                           *
 *                 P_Intercept = intercept in the prob. part of the model.   *
 *                 P_varname   = regression slope for covariate "varname"    *
 *                               in the prob. part of the model.             *
 *                 LogSDu1     = Log(Sqrt(Var_u2))                           *
 *                 z_u1u2      = Fisher transformation of Corr_u1u2:         *
 *                                  z = ln[(1+corr)/(1-corr)] / 2            *
 *                 Var_u1      = variance of the random effect in the        *
 *                               prob. part of the model.                    *
 *                 Cov_u1u2    = covariance of random effects u1 and u2.     *
 *                 Corr_u1u2   = correlation of random effects u1 and u2.    *
 *                                                                           *
 *             note: if specifying initial parameter estimates using the     *
 *                   init_parms option, the init_parms data set should have  *
 *                   the same variables as parms_u, except it should not     *
 *                   include var_e, var_u2, var_u1, cov_u1u2 or corr_u1u2    *
 *                   (these are derived parameters, i.e., functions of the   *
 *                    other parameters).                                     *
 *                                                                           *
 *   pred_x_u = data set containing predicted values for the model.          *
 *              pred_x_u contains all the variables in the input data set,   *
 *              plus the following variable:                                 *
 *                                                                           *
 *                 pred_x_a = predicted mean amount on consumption day.      *
 *                                                                           *
 *              if fitting the two-part model, then pred_x_u also contains   *
 *              the following variable:                                      *
 *                                                                           *
 *                  pred_x_p = predicted probability of consumption.         *
 *                                                                           *
 ****************************************************************************/

%macro NLMixed_Univariate (data           = ,
                        subject        = ,
                        repeat         = ,
                        response       = ,
                        modeltype      = TWOPART,
                        covars_prob    = ,
                        covars_amt     = ,
                        link           = PROBIT,
                        lambda         = ,
                        var_u1         = ,
                        var_u2         = ,
                        indep_u        = N,
                        replicate_var  = ,
                        nloptions      = ,
                        init_parms     = ,
                        print          = Y,
                        ntitle         = 2);

%let modeltype  = %upcase(&modeltype);
%let link       = %upcase(&link);
%let indep_u    = %upcase(%substr(&indep_u,1,1));
%let print      = %upcase(%substr(&print,1,1));

%if (&link      ^= LOGIT)  %then %let link = PROBIT;
%if (&modeltype ^= ONEPART) %then %let modeltype = TWOPART;

%if (&modeltype = ONEPART) %then %do;
  %let link = NONE;
  %let var_u1 = 0;
  %end;
%if (&var_u1 = 0 | &var_u2 = 0) %then %let indep_u = Y;

%if (&print = V) %then %let nloptions = &nloptions itdetails;

%if (&replicate_var ^= %str()) %then %do;
  footnote "Note: Standard Errors are not valid if data are not from a simple random sample";
  %end;

proc sort data=&data;  by &subject &repeat;  run;

data parms_u;
  run;

data pred_x_u;
  run;

data conv_u;
  status = 1;
  run;

 /* determine number of covariates in each part of the model. */
 /* get initial parameter estimate for box-cox transformation parameter lambda. */

data _null_;
  set &data;

  %if (&covars_prob = %str()) %then nx_p = 0%str(;);
  %else %do;
     array _p (*) &covars_prob;
     nx_p = dim(_p);
     %end;

  %if (&covars_amt = %str()) %then nx_a = 0%str(;);
  %else %do;
     array _a (*) &covars_amt;
     nx_a = dim(_a);
     %end;

  call symput("nx_p",trim(left(put(nx_p, 3.))));
  call symput("nx_a",trim(left(put(nx_a, 3.))));
  run;

 /* create indicator variable _consume that equals 0 if response = 0, 1 otherwise.          */

data _data;
  set &data;

  if (&response = 0) then _consume = 0;
  if (&response > 0) then _consume = 1;
  run;

 /*** calculate initial parameter estimates for the probability part of the model. ***/

%if (&init_parms = %str() & &modeltype ^= ONEPART) %then %do;

%if (&link = LOGIT) %then title%eval(&ntitle+1) "Calculate Initial Estimates for Probability Model (Link = Logit)" %str(;);
%else                     title%eval(&ntitle+1) "Calculate Initial Estimates for Probability Model (Link = Probit)" %str(;);

 /* fit probability model without random effects (proc logistic). */

data _conv_p;
  status = 1;
  run;

%if (&print ^= V) %then ods exclude all %str(;);
%if (&print = V) %then ods exclude oddsratios %str(;);

proc logistic data=_data namelen=50;
  model _consume(event="1") = &covars_prob / link=&link;
  ods output ParameterEstimates=_init_p;
  ods output ConvergenceStatus=_conv_p;
  run;

ods exclude none;

 /* check for convergenge. */

data _null_;
  set _conv_p;
  call symput("status",trim(left(put(status, 6.))));
  run;

%if (&status ^= 0) %then %goto exit;

 /* save parameters. */

proc transpose data=_init_p out=_init_p(drop=_name_);
  id variable;
  var estimate;
  run;

data _init_p;
  set _init_p;

  rename intercept = P_Intercept
    %do i = 1 %to &nx_p;
      %let xi = %scan(&covars_prob,&i,%str( ));
      &xi = P_&xi
      %end;
    ;
  run;

 /* fit probability model with random effects (proc nlmixed). */

%if (&var_u1 = %str()) %then %do;

data _init_p1;
  set _init_p;
  LogSDu1 = 0;
  run;

data _conv_p;
  status = 1;
  run;

%if (&print ^= V) %then ods exclude all %str(;);
    
proc nlmixed data=_data &nloptions;
  parms /data=_init_p1;

  %if (&replicate_var ^= %str()) %then replicate &replicate_var%str(;);    /* weight variable */

 /* reparameterize the variance/covariance of u. */

  VAR_U1 = EXP(2*LogSDu1);

  x_p = P_Intercept;
  %do i = 1 %to &nx_p;
    %let xi = %scan(&covars_prob,&i,%str( ));
    x_p = x_p + P_&xi * &xi;
    %end;
  xu_p = x_p + u1;

 /* probit or logit link. */

  %if (&link = LOGIT) %then _p = 1 /(1 + exp(-xu_p)) %str(;);
  %else _p = probnorm(xu_p) %str(;);

  model _consume ~ binary(_p);
           
  random u1 ~ normal(0,var_u1) subject=&subject;

  estimate "VAR_U1" VAR_U1;

  ods output ParameterEstimates=_init_p;
  ods output ConvergenceStatus=_conv_p;
  run;

%if (&print ^= V) %then ods exclude none %str(;);

 /* check for convergenge. */

data _null_;
  set _conv_p;
  call symput("status",trim(left(put(status, 6.))));
  run;

%if (&status ^= 0) %then %goto exit;

 /* save parameters. */

proc transpose data=_init_p out=_init_p(drop=_name_);
  id parameter;
  var estimate;
  run;

%end;  /* %if (&var_u1 = %str()) %then %do; */

title%eval(&ntitle+1);

%end;  /* %if (&init_parms = %str() & &modeltype ^= ONEPART) %then %do */

 /*** calculate initial parameter estimates for the amount part of the model (proc transreg). ***/

%if (&init_parms = %str()) %then %do;

%if (&print ^= V) %then ods exclude all %str(;);
ods output boxcox=_boxcox;
ods output coef=_coef;
ods output fitstatistics=_fitstatistics;

proc transreg data=_data(where=(&response > 0)) plots=none;
  %if (&lambda = %str()) %then %do;
    model boxcox(&response / lambda=0 to 1 by 0.1) = identity(&covars_amt) / ss2;
    %end;
  %else %do;
    model boxcox(&response / lambda=&lambda) = identity(&covars_amt) / ss2;
    %end;
  title%eval(&ntitle+1) "Calculate Initial Parameter Estimates for Amount Model";
  run;

title%eval(&ntitle+1);

ods exclude none;

 /* save parameters. */

data _boxcox (keep=_lambda);
  set _boxcox end=eof;
  retain _lambda _maxlike;
  if (_n_ = 1) then _maxlike = loglike;
  if (loglike >= _maxlike) then do;
    _lambda = lambda;
    _maxlike = loglike;
    end;
  if (eof) then output;
  run;

data _fitstatistics (keep=rmse);
  set _fitstatistics (rename=(value1=rmse));
  if (upcase(label1) = "ROOT MSE") then output; 
  run;

data _coef (keep=_varname coefficient);
  set _coef;
  if (upcase(variable) = "INTERCEPT") then _varname = variable; 
  else _varname = scan(variable,2,"()");
  output;
  run;

proc transpose data=_coef out=_coef(drop=_name_);
  id _varname;
  var coefficient;
  run;

data _init_a (drop=rmse);
  merge _coef _boxcox _fitstatistics;

  %if (&var_u2 = %str(0)) %then %do;
    var_e = rmse**2;
    %end;
  %else %do;
    var_u2 = rmse**2 / 2;
    var_e  = rmse**2 / 2;
    %end;
  
  rename _lambda    = A_Lambda
         intercept = A_Intercept
        %do i = 1 %to &nx_a;
          %let xi = %scan(&covars_amt,&i,%str( ));
          &xi = A_&xi
        %end;
    ;
  run;

 /*** combine initial parameter estimates from probability and amount parts of the model. ***/

data _init_parms (drop=var_e %if (&var_u2 ^= %str(0)) %then var_u2; %if (&lambda ^= %str()) %then a_lambda;);
  merge %if (&modeltype ^= ONEPART) %then _init_p; _init_a;

  A_LogSDe = log(var_e) / 2;
  %if (&var_u2 = %str()) %then LogSDu2 = log(var_u2) / 2 %str(;);
  %if (&indep_u ^= Y) %then Z_u1u2 = 0 %str(;);
  run;

%let init_parms = _init_parms;

%end;  /* %if (&init_parms = %str()) %then %do */

%if (&print = V) %then %do;
  proc print data=&init_parms noobs;
    title%eval(&ntitle+1) "Initial Parameter Estimates for Food/Nutrient Model";
    run;
  %end;

 /*** fit two-part (episodically consumed) or one-part (consumed every day) model (proc nlmixed). ***/

      %if (&link = NONE)  %then title%eval(&ntitle+1) "Food/Nutrient Model"%str(;);
%else %if (&link = LOGIT) %then title%eval(&ntitle+1) "Food/Nutrient Model (Link = Logit)"%str(;);
%else                           title%eval(&ntitle+1) "Food/Nutrient Model (Link = Probit)"%str(;);

data conv_u;
  status = 1;
  run;

%if (&print = N) %then ods exclude all %str(;);

proc nlmixed  data=&data &nloptions;
  parms / data=&init_parms;

  %if (&replicate_var ^= %str()) %then replicate &replicate_var%str(;);    /* survey weight variable */

  %if (&lambda ^= %str()) %then A_Lambda = &lambda%str(;);

  %if (&var_u1 = 0) %then u1 = 0%str(;);
  %if (&var_u2 = 0) %then u2 = 0%str(;);

 /* reparameterize the variance of e. */

  VAR_E = EXP(2 * A_LogSDe);

 /* reparameterize the variance/covariance of u. */

  %if (&var_u1 ^= 0) %then %do;
    %if (&var_u1 = %str()) %then VAR_U1 = exp(2 * LogSDu1)%str(;);
    %else VAR_U1 = &var_u1%str(;);
    %end;
  %if (&var_u2 ^= 0) %then %do;
    %if (&var_u2 = %str()) %then VAR_U2 = exp(2 * LogSDu2)%str(;);
    %else VAR_U2 = &var_u2%str(;);
    %end;
  %if (&var_u1 ^= 0 & &var_u2 ^= 0) %then %do;
    %if (&indep_u = Y) %then %do;
      CORR_U1U2 = 0;
      COV_U1U2  = 0;
      %end;
    %else %do;
      CORR_U1U2 =(EXP(2 * Z_u1u2) - 1) / (EXP(2 * Z_u1u2) + 1);
      COV_U1U2  = CORR_U1U2 * SQRT(VAR_U1 * VAR_U2);
      %end;
    %end;

 /* likelihood for probability of consumption. */

  %if (&modeltype = ONEPART) %then %do;
    ll_b = 0;
    %end;

  %else %do;
    x_p = P_Intercept;
    %do i = 1 %to &nx_p;
      %let xi = %scan(&covars_prob,&i,%str( ));
      x_p = x_p + P_&xi * &xi;
      %end;
    xu_p = x_p + u1;

 /* probit or logit link. */

    %if (&link = PROBIT) %then _p = probnorm(xu_p) %str(;);
    %else _p = 1 /(1 + exp(-xu_p)) %str(;);

    if (&response > 0) then _y = 1;
    else if (&response = 0) then _y = 0;
    else _y = .;
    if (_y ^= .) then ll_b = log((1-_p)**(1-_y)) + log(_p**(_y));
    else ll_b = .;
    %end;       

 /* likelihood for amount consumed on consumption day. */

  x_a = A_Intercept;
  %do i = 1 %to &nx_a;
    %let xi = %scan(&covars_amt,&i,%str( ));
    x_a = x_a + A_&xi * &xi;
    %end;
  xu_a = x_a + u2;

  pi = arcos(-1);

  if (&response > 0) then do; 
    if (A_Lambda = 0) then boxcoxy = log(&response); 
    else boxcoxy = (&response**A_Lambda - 1) / A_Lambda;
    ll_n = -log(sqrt(2 * pi * VAR_E)) - (boxcoxy - xu_a)**2 / (2 * VAR_E) +
           (A_Lambda - 1) * log(&response);
    end;

  %if (&modeltype ^= ONEPART) %then %do;
    else if (&response = 0) then do;
      ll_n = 0;
      end;
    %end;

  else ll_n = .;

  ll = ll_b + ll_n;

  model &response ~ general(ll);

 /* specify random effects. */

  %if (&var_u1 ^= 0 | &var_u2 ^= 0) %then %do;
          %if (&var_u2 = 0) %then random u1 ~ normal(0, VAR_U1);
    %else %if (&var_u1 = 0) %then random u2 ~ normal(0, VAR_U2);
    %else                         random u1 u2 ~ normal([0,0],[VAR_U1, COV_U1U2, VAR_U2]);
                                       subject=&subject /* out=pred_u */;
    %end;

 /* estimate additional parameters. */

  %if (&var_u1 ^= 0) %then estimate "VAR_U1" VAR_U1 %str(;);
  %if (&var_u2 ^= 0) %then estimate "VAR_U2" VAR_U2 %str(;);
  %if (&indep_u ^= Y) %then %do;
    estimate "COV_U1U2" COV_U1U2;
    estimate "CORR_U1U2" CORR_U1U2;
    %end;

  estimate "VAR_E" VAR_E;

  %if (&modeltype ^= ONEPART) %then predict x_p out=pred_x_p%str(;);
  predict x_a out=pred_x_a;
  predict xu_a out=pred_xu_a;

  ods output ParameterEstimates=parms_u;
  ods output AdditionalEstimates=adparms_u;
  ods output FitStatistics=fit_u;
  ods output ConvergenceStatus=conv_u;
  run;

title%eval(&ntitle+1);

%if (&print = N) %then ods exclude none %str(;);

 /* check for convergenge. */

data _null_;
  set conv_u;
  call symput("status",trim(left(put(status, 6.))));
  run;

%if (&status ^= 0) %then %goto exit;

 /*** save parameter estimates and predicted values. ***/

proc transpose data=parms_u out=parms_u(drop=_name_);
  id parameter;
  var estimate;
  run;

proc transpose data=adparms_u out=adparms_u(drop=_name_);
   id label;
   var estimate;
   run;

data parms_u;
  merge parms_u adparms_u;

 /* add fixed parameters. */

  %if (&lambda ^= %str()) %then A_Lambda = &lambda%str(;);

  %if (&modeltype ^= ONEPART & &var_u1 ^= %str()) %then VAR_U1 = &var_u1%str(;);
  %if (&var_u2 ^= %str()) %then VAR_U2 = &var_u2%str(;);
  %if (&modeltype ^= ONEPART & &indep_u = Y) %then %do;
    COV_U1U2  = 0;
    CORR_U1U2 = 0;
    %end;
  run;

data pred_x_u;
  merge &data
        %if (&modeltype ^= ONEPART) %then pred_x_p(keep=&subject &repeat pred rename=(pred=pred_x_p));
                                          pred_x_a(keep=&subject &repeat pred rename=(pred=pred_x_a));
    by &subject &repeat;

  label %if (&modeltype ^= ONEPART) %then pred_x_p = " "; pred_x_a = " ";
  run;

 /*** delete unneeded data sets. */

proc datasets lib=work nolist;
  delete _data _init_a _init_parms adparms_u pred_x_a
        %if (&modeltype ^= ONEPART) %then _init_p _conv_p pred_x_p; %str(;)
  run;
  quit;

%exit:

title%eval(&ntitle+1);

%if (&replicate_var ^= %str()) %then footnote %str(;);

%mend NLMixed_Univariate;
