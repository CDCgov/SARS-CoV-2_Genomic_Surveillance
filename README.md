# Genomic Surveillance of SARS-CoV-2 Circulating in the United States

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 

## Project Overview 
The emergence and rapid expansion of multiple [SARS-CoV-2 variants](https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-classifications.html) of interest (VOIs), variants of concern (VOCs), and variants being monitored (VBMs) highlight the need for robust genomic surveillance to monitor circulating viruses and help guide the public health response to the COVID-19 pandemic [(Galloway SE et al. 2021)](https://www.cdc.gov/mmwr/volumes/70/wr/mm7003e2.htm). November 2020, CDC established national genomic surveillance of SARS-CoV-2; activities were scaled up rapidly in early 2021. This repository contains the code and an example dataset to demonstrate the statistical weighting and modeling methods used to generate representative estimates of proportions of variants in the United States. These methods are used to produce the SARS-CoV-2 [variant proportion estimates on CDC's COVID Data tracker](https://covid.cdc.gov/covid-data-tracker/#variant-proportions). The genomic surveillance of SARS-CoV-2 circulating in the United States was described between December 2020-May 2021 by [Paul P et al. (2021)](https://www.cdc.gov/mmwr/volumes/70/wr/mm7023a3.htm) and between June 2021-December 2021 by Lambrou AS, Shirk P, et al. (2022).

## Variant Estimation Methods
 SARS-CoV-2 consensus sequences submitted or tagged for national genomic surveillance are combined, assessed for quality, deduplicated, and analyzed weekly to estimate the proportions of variants circulating at the national, HHS regional, and jurisdiction levels. SARS-CoV-2 variant proportions are estimated weekly for variants of concern (VOCs), variants of interest (VOIs), variants being monitored (VBM), and any other lineages accounting for >1% of sequences nationally during the preceding 12 weeks. All analyses use Pango SARS-CoV-2 lineage nomenclature and sublineages. There are two variant estimation metrics calculated in the corresponding R code: 
 1. <b>Estimated weighted proportions:</b> Weighted analysis using complex survey design methods to produce weekly proportions for the past 3-12 weeks. Because of the lag between sample collection and when sequence data are available, there are not sufficient data to estimate proportions for most recent 1-2 weeks using complex survey design method. 
2. <b>Nowcast model:</b> multinomial regression analysis to produce estimates of proportions for the most recent 1-2 weeks to account for time lag between sample collection and sequence reporting. 

## Weighting and Weighted Proportion Estimates 

Weights are estimated by treating each sequencing laboratory source as a cluster nested with strata (geographic level such as jurisdiction or region), calculated for each week (i.e., to ensure representation weekly). Two sets of weights are estimated: the first, w<sub>p</sub>, for representation among (PCR) test positive individuals; and second, w<sub>i</sub>, for representation among all prevalent infections. Assumptions (modifiable as data become available) in estimating weights are: 
* Each positive (PCR) test is in the sampling frame of one of the source streams, so that for each state, week and source 

<img src="https://latex.codecogs.com/svg.image?w_p&space;=&space;\frac{\mbox{number&space;of&space;positive&space;PCR&space;test&space;results}}{\mbox{number&space;of&space;sequences&space;submitted}}" title="w_p = \frac{\mbox{number of positive PCR test results}}{\mbox{number of sequences submitted}}" />

* Oversampling of S-gene target failure (SGTF) samples by one source results in a reduction in weights of SGTF sequences from that source by a factor that is estimated using a logistic regression model relating the odds of finding an “SGTF variant” by source, geographic level, and week. 
* Estimation of w<sub>i</sub> involves estimating the (unobserved) number of infections from test results. There is no reliable and precise method for this yet, so these weights are subject to considerable uncertainties. Here, a [strategy](https://www.medrxiv.org/content/10.1101/2020.10.07.20208504v2.full) based on test positivity is used:

<img src="https://latex.codecogs.com/svg.image?\frac{\mbox{number&space;of&space;prevalent&space;infections}}{\mbox{number&space;of&space;test&space;positives}}&space;=&space;\sqrt{\frac{\mbox{population&space;of&space;jurisdiction}}{\mbox{number&space;of&space;tests}}}" title="\frac{\mbox{number of prevalent infections}}{\mbox{number of test positives}} = \sqrt{\frac{\mbox{population of jurisdiction}}{\mbox{number of tests}}}" />

If each source (lab) stream is assumed to sample from a base population with the same prevalence of infection as the jurisdiction, it can be shown that the weight specific to each source, based on the test positivity of each source, is: 

<img src="https://latex.codecogs.com/svg.image?w_i&space;=&space;\frac{\mbox{number&space;of&space;prevalent&space;infections&space;in&space;source}}{\mbox{number&space;of&space;source&space;positives}}&space;=&space;\frac{\mbox{number&space;of&space;source&space;positives}}{\mbox{number&space;of&space;source&space;tests}}&space;&space;&space;{\frac{\sqrt{\mbox{population&space;of&space;jurisdiction}\times\mbox{number&space;of&space;tests}}}{\mbox{number&space;of&space;positives}}}" title="w_i = \frac{\mbox{number of prevalent infections in source}}{\mbox{number of source positives}} = \frac{\mbox{number of source positives}}{\mbox{number of source tests}} {\frac{\sqrt{\mbox{population of jurisdiction}\times\mbox{number of tests}}}{\mbox{number of positives}}}" />

The “infection” weight for each sequence is w<sub>p</sub> x w<sub>i</sub> and depends on stratum (geographic level), time (week of collection), and cluster (source). 

Three separate survey designs are used in this analysis: 
* Unweighted to estimate variant prevalence among sequenced samples 
* Weighted for estimation among test positives 
* Weighted for estimation among infections, but using strategies of unproven reliability in this context 

## Nowcast Projection Model 

A simple multivariant model with complete cross-immunity and where the variants (labeled i) differ in transmissibility (&beta;<sub>i</sub>) alone is: 

<img src="https://latex.codecogs.com/svg.image?\frac{dv_{i}}{dt}=\left&space;(&space;\beta_{i}s&space;-&space;\gamma&space;&space;\right&space;)v_{i}" title="\frac{dv_{i}}{dt}=\left ( \beta_{i}s - \gamma \right )v_{i}" />

Here s is the proportion susceptible in the population, and &gamma; is the recovery rate. The formal solution

<img src="https://latex.codecogs.com/svg.image?v_{i}\left&space;(&space;t&space;\right&space;)&space;=&space;v_{i}\left&space;(&space;0&space;\right&space;)e^{\int&space;\left&space;(&space;\beta_{i}s&space;-&space;\gamma&space;\right&space;)dt}" title="v_{i}\left ( t \right ) = v_{i}\left ( 0 \right )e^{\int \left ( \beta_{i}s - \gamma \right )dt}" />

so that 

<img src="https://latex.codecogs.com/svg.image?\frac{v_i(t)}{v_j(t)}&space;=&space;\frac{v_i(0)}{v_j(0)}e^{\int&space;(\beta_i&space;-\beta_j)s\,dt}" title="\frac{v_i(t)}{v_j(t)} = \frac{v_i(0)}{v_j(0)}e^{\int (\beta_i -\beta_j)s\,dt}" />

If, as is often the case, the integrand is a slowly varying function of time, a multinomial log-linear model may be a good fit for the variant proportions 

<img src="https://latex.codecogs.com/svg.image?p_i(t)&space;=&space;\frac{v_i(t)}{\sum_j&space;v_j(t)}" title="p_i(t) = \frac{v_i(t)}{\sum_j v_j(t)}" />

The model coefficients b<sub>0i</sub> and b<sub>1i</sub> are such that 

<img src="https://latex.codecogs.com/svg.image?p_i(t)&space;=&space;\frac{e^{b_{0i}&space;&plus;&space;b_{1i}t}}{\sum_j&space;e^{b_{0j}&space;&plus;&space;b_{1j}t}}" title="p_i(t) = \frac{e^{b_{0i} + b_{1i}t}}{\sum_j e^{b_{0j} + b_{1j}t}}" />

A useful result that follows is 

<img src="https://latex.codecogs.com/svg.image?\frac{d\log&space;p_i}{dt}&space;=&space;b_{1i}&space;-&space;\sum_j&space;p_jb_{1j}" title="\frac{d\log p_i}{dt} = b_{1i} - \sum_j p_jb_{1j}" />

With two variants, the model simplifies to binomial logistic regression. The statistical model used for multinomial regression is 

<img src="https://latex.codecogs.com/svg.image?\mbox{variant}&space;\sim&space;\mbox{collection-week}" title="\mbox{variant} \sim \mbox{collection-week}" />

for the national projections, and 

<img src="https://latex.codecogs.com/svg.image?\mbox{variant}&space;\sim&space;\mbox{collection-week}&space;&plus;&space;\mbox{HHS-region}" title="\mbox{variant} \sim \mbox{collection-week} + \mbox{HHS-region}" />

for the regional projections. 




  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
