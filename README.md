National Cancer Institute (NCI) methods workshop
================
Didier Brassard

# Introduction

This workshop is an introduction to random measurement errors in
nutrition research. While it is obvious that a systematic error
(difference) between the “true” value and its measurement can be a
problem, the impact of random errors is often more subtle. However, in
many cases, random errors can be as problematic as systematic errors if
they are ignored.

Before diving into measurement error correction methods, I recommend two
introduction blog posts on the topic: [*‘Statistical concept you should
know’: random and systematic measurement
errors*](https://didierbrassard.github.io/posts/2022/11/blog-post-6/)
and [*Impact of random errors: two nutrition
examples*](https://didierbrassard.github.io/posts/2022/11/blog-post-7/).
Useful reference includes *Using Short-Term Dietary Intake Data to
Address Research Questions Related to Usual Dietary Intake among
Populations and Subpopulations: Assumptions, Statistical Techniques, and
Considerations* by [Kirkpatrick et
al. (2022)](https://pubmed.ncbi.nlm.nih.gov/35283362/) as well as the
*STRATOS guidance document on measurement error and misclassification of
variables in observational epidemiology: Part 1-Basic theory and simple
methods of adjustment* by [Keogh et
al. (2020)](https://pubmed.ncbi.nlm.nih.gov/32246539/)

## Do I need measurement error correction?

Food-frequency questionnaires or screeners reflect **long-term** intakes
and measurement error correction methods - at least for
within-individual random errors - are not needed. On the contrary,
measurement error correction methods may be required when using
**short-term** dietary intake instruments, e.g., food records or 24-hour
dietary recalls.

The flow chart below is a general guide to decide whether measurement
error correction (e.g., with the NCI methods) may be required when using
dietary intakes data from short-term instruments. Please note that in
the context of this repository, measurement error correction refers to
random errors.

``` mermaid

%%{init: {'theme': 'neutral' } }%%
flowchart TD
  Q1("What is the statistic of interest?")
  S1{"Mean<br>only"}
  S2{"Distribution or<br>prevalence"}
  S3{"Regression<br>coefficient"}
  S1-->Q4("Is there a <b>scoring</b> <br>(e.g., diet quality score)?")
  Q4--Yes-->F1
  Q4--No-->F2
  F1{{"Use the<br>population ratio<br>method"}}
  Q1-->S1
  Q2("Are <b>repeated data</b><br>available for n>50?")
  Q1-->S2
  Q1-->S3
  S3-->Q3("Is dietary intake<br>the <b>independent</b> variable?")
  Q3--"Yes"-->Q2
  Q3--"No"-->F2
  F2{{"Measurement error<br>correction not needed"}}
  S2-->Q2
  Q2--Yes-->F3
  F3{{"Apply measurement error<br>correction methods"}}
  Q2--No-->F4
  F4{{"Use raw data,<br>assess bias with simulation,<br> acknowledge limitation"}}

```

For more details about the population ration method, see [*Analyzing the
Canadian Community Health Survey (CCHS) 2015 data with R: mean diet
quality
score*](https://didierbrassard.github.io/posts/2022/12/blog-post-8/)

# NCI univariate method, distribution

SAS codes (01, 11, 12, …) are provided as an introduction to basic
aspects of random measurement errors and application of the NCI
univariate method to generate a distribution of usual dietary intakes.

The *12_Univariate_2part_exercise.sas* code is intended as an exercise
to practice using the univariate method. The answers are provided in
*12_Univariate_2part.sas*.

A companion presentation is [available
here](https://drive.google.com/file/d/1VKxCEDkiGSCbAYVioob-s4krjDJut4uT/view?usp=sharing).

The original NCI macros and further examples are available on the
[Biometry Research Group’s *Software for Measurement Error in Nutrition
Research*
website](https://prevention.cancer.gov/research-groups/biometry/measurement-error-impact/software-measurement-error)

## Schematic overview of the NCI methods

``` mermaid
%%{init: {'theme': 'neutral' } }%%
flowchart TB
subgraph data[<b>Repeated dietary assessment</b>]
    R1("24-hour recall #1")
    R2("24-hour recall #2")
    R3("24-hour recall #J")
    RT("Dietary constituents")
    R1 & R2 & R3-->RT
end
subgraph C[<b>Covariates</b>]
    N(Sequence,<br>weekend,<br>season)
    Z(Subgroups)
end
subgraph A[<b>Assumptions</b>]
    S1(Statistical<br>assumptions)
    S2(Classical error<br>model assumptions)
end
    data --"Box-Cox<br>transformation"-->M[<b>Step 1<br>Measurement error model</b>]
    C-->M
    A-->M
    M-->P1["Model parameters"]
    P1-->P2["Predicted intakes"]
    M-->P3[Within- & between-<br>individual variance]
    P3-->P2
    P1 & P2-->MC["<b>Step 2<br>Monte Carlo simulations</b>"]
    MC--"<I>M</I> simulations<br>per individual"-->PI["<b>Step 3<br>'Usual intakes' distribution<br>calculation among<br>(<I>M*n</I>) pseudo-individuals"]
    
```
