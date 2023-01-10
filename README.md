National Cancer Institute (NCI) methods workshop
================
Didier Brassard

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

## NCI univariate method, distribution

SAS codes (01, 02, … 04) are provided as an introduction to basic
aspects of random measurement errors and application of the NCI
univariate method to generate a distribution of usual dietary intakes.

The *03_Univariate_2part_exercise.sas* code is intended as an exercise
to practice using the univariate method. The answers are provided in
*03_Univariate_2part.sas*.

A companion presentation is [available
here](https://drive.google.com/file/d/1VKxCEDkiGSCbAYVioob-s4krjDJut4uT/view?usp=sharing).

The original NCI macros and further examples are available on the
[Biometry Research Group’s *Software for Measurement Error in Nutrition
Research*
website](https://prevention.cancer.gov/research-groups/biometry/measurement-error-impact/software-measurement-error)

### Schematic overview

``` mermaid
%%{init: {'theme': 'neutral' } }%%
flowchart TB
    subgraph data[<b>Repeated dietary assessment</b>]
        R1("24-hour recall #1")
        R2("24-hour recall #2")
        R3("24-hour recall #J")
    end
    subgraph C[<b>Covariates</b>]
        N(Sequence,<br>weekend, season)
        Z(Subgroups)
    end
    subgraph A[<b>Assumptions</b>]
        S1(Statistical<br>assumptions)
        S2(Classical error<br>model assumptions)
    end
        data -->M[<b>Measurement error model</b>]
        C-->M
        A-->M
        M-->P1["Model parameters"]
        P1-->P2["Predicted intakes"]
        M-->P3[Within- & between-<br>individual variance]
        P3-->P2
        P1 & P2-->MC["<b>Monte Carlo simulations</b>"]
        MC--"<I>M</I> simulations<br>per individual"-->PI["<b>'Usual intakes' among<br>(<I>M*n</I>) pseudo-individuals"]
    
```
