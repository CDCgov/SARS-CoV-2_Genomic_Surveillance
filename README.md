# Genomic Surveillance of SARS-CoV-2 Circulating in the United States

**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm). GitHub is not hosted by the CDC but is a third-party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 

## Project Overview 
The emergence and rapid expansion of multiple [SARS-CoV-2 variants](https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-classifications.html) highlights the need for robust genomic surveillance to monitor circulating viruses and help guide the public health response to the COVID-19 pandemic [(Galloway SE et al. 2021)](https://www.cdc.gov/mmwr/volumes/70/wr/mm7003e2.htm). In November 2020, CDC established national genomic surveillance of SARS-CoV-2; activities were scaled up rapidly in early 2021. This repository contains the code and an example dataset to demonstrate the statistical weighting and modeling methods used to generate representative estimates of proportions of variants in the United States. These methods are used to produce the SARS-CoV-2 [variant proportion estimates on CDC's COVID Data tracker](https://covid.cdc.gov/covid-data-tracker/#variant-proportions). These methods have been updated over time to maintain timely and representative estimates with changes in data availability and sequencing and testing practices. Genomic surveillance of SARS-CoV-2 variants circulating in the United States was described between December 2020-May 2021 by [Paul P et al. (2021)](https://www.cdc.gov/mmwr/volumes/70/wr/mm7023a3.htm), between June 2021-January 2022 by [Lambrou AS, Shirk P, et al. (2022)](https://www.cdc.gov/mmwr/volumes/71/wr/mm7106a4.htm?s_cid=mm7106a4_w), and between January 2022-May 2023 by [Ma KC, Shirk P, et al. (2023)](https://www.cdc.gov/mmwr/volumes/72/wr/mm7224a2.htm?s_cid=mm7224a2_w). 

## Variant Estimation Methods
SARS-CoV-2 consensus sequences submitted or tagged for national genomic surveillance are combined, assessed for quality, deduplicated, and analyzed biweekly to estimate the proportions of major variants circulating at the national and HHS regional levels. SARS-CoV-2 variant proportions are estimated every other week for lineages accounting for ≥1% of sequences nationally during any of the 6 two-week periods prior to the most recent 4 weeks. The nowcast model also includes any lineages ≥0.5% during the two-week period prior to the current two-week period. All analyses use Pango SARS-CoV-2 lineage nomenclature and sublineages. There are two variant estimation metrics calculated in the corresponding R code:
1. <b>Estimated weighted proportions:</b> Weighted analysis using complex survey design methods to produce proportions for the 6 most recent two-week periods prior to the most recent 4 weeks. Because of the lag between sample collection and when sequence data are available, there are not sufficient data to estimate proportions for most recent 2 two-week periods using the complex survey design method. 
2. <b>Nowcast model:</b> Survey design-based multinomial regression analysis to produce estimates of proportions for the two most recent two-week periods (i.e., 4 weeks in total) to account for time lag between sample collection and sequence reporting. 

All variant proportion and nowcast analyses use weights to ensure estimates are representative of the U.S. national and HHS regional levels. Our survey design uses laboratory source as clusters nested within state-weeks as strata and weights based on the weekly estimated number of infections represented by each sequence. The variance of all variant proportion and nowcast analyses are adjusted to account for the survey design.

## Weighting and Weighted Proportion Estimates 

Ideally weights would be calculated for each cluster in each stratum of the survey design (i.e., each combination of laboratory-state-week). In the absence of laboratory sampling frame to calculate those weights, we estimate weights for each stratum (i.e., each state-week combination). 

Beginning May 13, 2023, methodological changes were made to genomic surveillance weighting methods following the expiration of the public health emergency declaration in response to declining numbers of cases and sequenced specimens. Regional-level test positivity data from the [National Respiratory and Enteric Virus Surveillance System (NREVSS)](https://www.cdc.gov/surveillance/nrevss/index.html) were used in calculating survey weights, replacing state-level test positivity from the [COVID-19 electronic laboratory reporting (CELR)](https://healthdata.gov/dataset/COVID-19-Diagnostic-Laboratory-Testing-PCR-Testing/j8mb-icvb) data. The total numbers of positive tests still come from CELR, but due to inconsistent availability of state-level data, the number of infections is now estimated at the HHS Region level and distributed among constituent states based on population. NREVSS data are available on [data.cdc.gov](https://data.cdc.gov/Laboratory-Surveillance/Percent-Positivity-of-COVID-19-Nucleic-Acid-Amplif/gvsb-yw6g). Details follow:

We estimate the number of infections represented by each sequence using the methods of [Chiu MA and Ndeffo-Mbah ML (2021)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009374). We use NREVSS testing data for test positivity and CELR data for positive tests to calculate weights at the HHS regional level:

$$
weight_{region} = \frac{\sqrt{ \frac{(number\ of\ positive\ NREVSS\ tests)}{(number\ of\ NREVSS\ tests)} * (population\ reporting) * (number\ of\ positive\ CELR\ tests) }}{(number\ of\ sequences)}
$$

where the "population reporting" includes the populations of states reporting testing data to CELR in a given week. To calculate these regional weights, we aggregate numbers of tests, population, and sequences to the HHS region each week. We then apportion the aggregate weekly regional weights to individual states by their respective total populations: 

$$
weight_{state} = \frac{population_{state}}{population_{region}} * weight_{region}
$$

## Nowcast Projection Model 

A simple multivariant model which assumes cross-immunity in the population to each variant and where the variants (labeled $i$) differ in transmissibility ($\beta_i$) alone is:

$$
\frac{dv_i}{dt} = (\beta_i*s - \gamma) * v_i
$$

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; $v_i$ = number of individuals infected with variant $i$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; $\beta_i$ = transmissibility of variant $i$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; $s$ = proportion of susceptible individuals in the population  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; $\gamma$ = the recovery rate  

The formal solution for the number of individuals infected with variant i at time t is: 

$$
v_i(t) = v_i(0) * e^{ \int{(\beta_i * s - \gamma)dt} }
$$

and the ratio of proportions, $p$, of individuals infected with variants $i$ and $j$ at time $t$ is:

$$
\frac{p_i(t)}{p_j(t)} = \frac{v_i(t)}{v_j(t)} = \frac{v_i(0)}{v_j(0)} * e^{ \int{(\beta_i - \beta_j) * s dt} }
$$

The proportion of variant $i$ at time $t$ is: 

$$
p_i(t) = \frac{v_i(t)}{\sum_j{v_j(t)}}
$$

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; where $\sum_j{v_j(t)}$ = the total number of individuals infected at time $t$,

and the variant proportion can be estimated using a multinomial logistic regression model. If the integral in the function of the ratio of proportions above is linear in time, the multinomial model can use a linear covariate of time (which is what we use). If the integral is non-linear in time, as would be the case if transmissibility or the proportion of susceptible individuals changed rapidly over the time period of data used to fit the multinomial model, then a non-linear covariate of time would be warranted. The multinomial model's estimate of the proportion of variant $i$ at time $t$ is:

$$
p_i(t) = \frac{e^{b_{0i} + b_{1i} * t}}{\sum_j{e^{b_{0j} + b_{1j} * t}}}
$$

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; $b_{0i}$ = intercept for variant $i$  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; $b_{1i}$ = coefficient of the time covariate for variant $i$  

A useful extension of the multinomial model is to estimate lineage growth rates as the rate of change of the log proportion of variant $i$ at time $t$:

$$
\frac{dlogp_i}{dt} = b_{1i} - \sum_j{ p_j(t) * b_{1j} }
$$

With two variants, the model simplifies to binomial logistic regression.

The statistical model used for multinomial regression is:

<p align="center"> variant ~ collection-week </p>

for the national projections, and

<p align="center"> variant ~ collection-week + HHS-region </p>

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