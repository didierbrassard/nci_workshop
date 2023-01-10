National Cancer Institute (NCI) methods workshop
================

This workshop is an introduction to random measurement errors in
nutrition research. While it is obvious that a systematic error
(difference) between the “true” value and its measurement can be a
problem, the impact of random errors is often more subtle. However, in
many cases, random errors can be as problematic as systematic errors if
they are ignored.

Before diving into measurement error correction methods, I recommend two
introduction blog posts on the topic: *[‘Statistical concept you should
know’: random and systematic measurement
errors](https://didierbrassard.github.io/posts/2022/11/blog-post-6/)*
and *[Impact of random errors: two nutrition
examples](https://didierbrassard.github.io/posts/2022/11/blog-post-7/)*.
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

## Schematic overview

``` mermaid
%%{init: {&#39;theme&#39;: &#39;neutral&#39; } }%%
flowchart TB
subgraph data[&lt;b&gt;Repeated dietary assessment&lt;/b&gt;]
    R1(&quot;24-hour recall #1&quot;)
    R2(&quot;24-hour recall #2&quot;)
    R3(&quot;24-hour recall #J&quot;)
end
subgraph C[&lt;b&gt;Covariates&lt;/b&gt;]
    N(Sequence,&lt;br&gt;weekend, season)
    Z(Subgroups)
end
subgraph A[&lt;b&gt;Assumptions&lt;/b&gt;]
    S1(Statistical&lt;br&gt;assumptions)
    S2(Classical error&lt;br&gt;model assumptions)
end
    data --&gt;M[&lt;b&gt;Measurement error model&lt;/b&gt;]
    C--&gt;M
    A--&gt;M
    M--&gt;P1[&quot;Model parameters&quot;]
    P1--&gt;P2[&quot;Predicted intakes&quot;]
    M--&gt;P3[Within- &amp; between-&lt;br&gt;individual variance]
    P3--&gt;P2
    P1 &amp; P2--&gt;MC[&quot;&lt;b&gt;Monte Carlo simulations&lt;/b&gt;&quot;]
    MC--&quot;&lt;I&gt;M&lt;/I&gt; simulations&lt;br&gt;per individual&quot;--&gt;PI[&quot;&lt;b&gt;&#39;Usual intakes&#39; among&lt;br&gt;(&lt;I&gt;M*n&lt;/I&gt;) pseudo-individuals&quot;]

```
