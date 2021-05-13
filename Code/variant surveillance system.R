####################################################################################
# Code accompanying MMWR on Genomic Surveillance 
# Date: 2021-05-11
# Author: Prabasaj Paul (VIG5@CDC.GOV)
####################################################################################

####################################################################################
# Part 1: From data to survey design
####################################################################################

library(survey)
options(survey.adjust.domain.lonely=T,survey.lonely.psu="average") 

# Sequence data, test tallies by state and collection date, state populations: 
load("variant_survey_dat_example.RData")


## Some general parameters:
week0day1 = as.Date("2020-01-05") # First Sunday in 2020
# HHS regions
HHS_reg = list(HHS1=c("CT", "ME", "MA", "NH", "RI", "VT"),
               HHS2=c("NJ", "NY", "PR", "VI"),
               HHS3=c("DE", "DC", "MD", "PA", "VA", "WV"),
               HHS4=c("AL", "FL", "GA", "KY", "MS", "NC", "SC", "TN"),
               HHS5=c("IL", "IN", "MI", "MN", "OH", "WI"),
               HHS6=c("AR", "LA", "NM", "OK", "TX"),
               HHS7=c("IA", "KS", "MO", "NE"),
               HHS8=c("CO", "MT", "ND", "SD", "UT", "WY"),
               HHS9=c("AZ", "CA", "HI", "NV", "AS", "MP", "GU", "MH"),
               HHS10=c("AK", "ID", "OR", "WA")
)
# Territory 2020-07-01 populations from https://en.wikipedia.org/wiki/List_of_states_and_territories_of_the_United_States_by_population (2021-04-19)
# Marshall Islands 2018 estimate https://en.wikipedia.org/wiki/Marshall_Islands (2021-04-19)
pops = rbind(pops, data.frame(STUSAB=c("AS", "GU", "MP", "VI", "MH"), `Total.`=c(49437, 168485, 51433, 106235, 58413)))
hhs = data.frame(STUSAB=toupper(pops$STUSAB))
hhs$HHS = sapply(hhs$STUSAB, grep, HHS_reg)

## Subset by time and place: USA, states with two letter abbreviations, human host, reasonable collection date
# U.S., legitimately coded states, human host; dat may already be subset to US, but retained for backward compatibility
us.dat = subset(dat, primary_country %in% c("United States", "USA") & primary_host=="Human") 
us.dat = subset(us.dat, nchar(primary_state_abv)==2) 
## Get unique records in Hadoop table (2021-04-26)
us.dat = unique(us.dat)
pangolin = unique(pangolin)
baseline = unique(baseline)
# Disambiguate and remove unreasonable dates
us.dat$collection_date = as.Date(us.dat$primary_collection_date)
us.dat$collection_date = with(us.dat, as.Date(ifelse(is.na(collection_date), as.Date(primary_collection_date_dt), collection_date), origin="1970-01-01"))
if ("covv_subm_date" %in% names(us.dat)) us.dat$covv_subm_date = as.Date(us.dat$covv_subm_date)
if ("contractor_receive_date_to_cdc" %in% names(us.dat)) us.dat$contractor_receive_date_to_cdc = as.Date(us.dat$contractor_receive_date_to_cdc)
us.dat$yr_wk = as.character(us.dat$collection_date - as.numeric(strftime(us.dat$collection_date, format="%w"))) 
us.dat$DAY = as.numeric(us.dat$collection_date - week0day1) 
us.dat = subset(us.dat, collection_date >= as.Date("2019-10-01")) 
us.dat = merge(us.dat, pangolin[, c("nt_id", "lineage")], by.x="primary_nt_id", by.y="nt_id", all.x=TRUE)
if ("covv_lineage" %in% names(us.dat)) us.dat$lineage = with(us.dat, ifelse(is.na(lineage), covv_lineage, lineage))
# NS3 identifier in baseline$source:
us.dat = merge(us.dat, baseline[, c("nt_id", "source", "primary_virus_name", "s_mut")],
               by.x=c("primary_nt_id", "primary_virus_name"), by.y=c("nt_id", "primary_virus_name"), all.x=TRUE)
## Add in HHS regions (2021-03-09)
us.dat = merge(us.dat, hhs, by.x="primary_state_abv", by.y="STUSAB", all.x=TRUE)



## Testing data by state and date (HHS Protect) for weighting
# Done: tests = read.csv("surveillance/tests_by_state_collection_date.csv")
names(tests) = gsub("^[^a-zA-Z]+", "", names(tests)) # Odd first name
names(tests) = gsub("^[[:print:]]+state$", "STUSAB", names(tests)) # Clean an odd name
tests$collection_date = as.Date(tests$collection_date)
tests = subset(tests, collection_date >= as.Date("2019-10-01") & collection_date <= Sys.Date())
tests$yr_wk = as.character(tests$collection_date - as.numeric(strftime(tests$collection_date, format="%w"))) 
tests$TOTAL = rowSums(tests[, c("INDETERMINATE", "INVALID", "NEGATIVE", "POSITIVE")], na.rm=TRUE)
tests$POSITIVE = ifelse(is.na(tests$POSITIVE), 0, tests$POSITIVE)
# Aggregate by week for weighting
tests_wk = expand.grid(STUSAB=unique(tests$STUSAB), yr_wk=unique(tests$yr_wk), stringsAsFactors=FALSE)
for (cc in c("POSITIVE", "TOTAL")) {
  tests_wk = merge(tests_wk, data.frame(xtabs(tests[, cc] ~ STUSAB + yr_wk, tests)), all.x=TRUE)
  names(tests_wk) = gsub("[Ff]req", cc, names(tests_wk)) 
}



## Test tally denominator streams
test_tallies_wk = merge(merge(tests_wk, hhs, all.x=TRUE), data.frame(STUSAB=toupper(pops$STUSAB), state_population=pops$Total.), all.x=TRUE)
test_tallies_wk = within(test_tallies_wk, INFECTIONS <- ifelse(TOTAL>0, POSITIVE * sqrt(state_population/TOTAL), 0))
test_tallies_wk$week = as.numeric(as.Date(test_tallies_wk$yr_wk) - week0day1)%/%7
incidence_by_region = merge(
  aggregate(INFECTIONS ~ yr_wk + HHS, data=test_tallies_wk, sum),
  aggregate(state_population ~ yr_wk + HHS, data=test_tallies_wk, sum)
)
incidence_by_region$HHS_INCIDENCE = incidence_by_region$INFECTIONS/incidence_by_region$state_population

## SGTF over-sampling weights 
sgtf.1 = table(us.dat$lineage, paste(us.dat$contractor_vendor_name, us.dat$contractor_targeted_sequencing))
sgtf.vars = rownames(sgtf.1)[ sgtf.1[, grep("dropout$", colnames(sgtf.1))] > sgtf.1[, grep("Illumina $", colnames(sgtf.1))] ]
# Smoothed weights: use set of states and weeks where targeted samples were sequenced
sgtf.sub = unique(subset(us.dat, contractor_targeted_sequencing %in% "Screened for S dropout")[, c("primary_state_abv", "yr_wk")])
sgtf.glm = glm(
  I(lineage %in% sgtf.vars)~ I(contractor_vendor_name %in% "Helix/Illumina") + primary_state_abv + yr_wk,
  family="binomial",
  subset(us.dat, yr_wk %in% sgtf.sub$yr_wk & primary_state_abv %in% sgtf.sub$primary_state_abv))
sgtf.sub$sgtf_weights = predict(sgtf.glm, cbind(contractor_vendor_name=c("Helix/Illumina"), sgtf.sub), type="response")/
  predict(sgtf.glm, cbind(contractor_vendor_name=c("Other"), sgtf.sub), type="response")
test_tallies_wk = merge(test_tallies_wk, sgtf.sub, by.x=c("STUSAB", "yr_wk"), by.y=c("primary_state_abv", "yr_wk"), all.x=TRUE)
test_tallies_wk$sgtf_weights = ifelse(is.na(test_tallies_wk$sgtf_weights), 1, test_tallies_wk$sgtf_weights) # Inconsequential, but avoids disruptions due to NA


## Creating a trimmed down survey dataset
# Removing submission/receive dates for now, for consistency with frozen dataset
svy.dat = data.frame(STUSAB=us.dat$primary_state_abv, HHS=us.dat$HHS, yr_wk=us.dat$yr_wk, DAY=us.dat$DAY,
                     # SUBM_DT=us.dat$covv_subm_date, CDC_DT = us.dat$contractor_receive_date_to_cdc,
                     LAB=toupper(us.dat$source),
                     SGTF_UPSAMPLING=(us.dat$contractor_targeted_sequencing %in% "Screened for S dropout"), SOURCE=us.dat$source,
                     VARIANT=us.dat$lineage, S_MUT=us.dat$s_mut)
svy.dat$LAB[is.na(svy.dat$LAB)] = "OTHER"
svy.dat$week = as.numeric(as.Date(svy.dat$yr_wk) - week0day1)/7
svy.dat = merge(svy.dat, data.frame(STUSAB=toupper(pops$STUSAB), state_population=pops$Total.), all.x=TRUE)



## Weights
# Infections to test-positive ratio
# Geometric mean strategy: https://www.medrxiv.org/content/10.1101/2020.10.07.20208504v2.full-text
#   infections/test positives = sqrt(population/number tested)
svy.dat = merge(svy.dat, cbind(test_tallies_wk[, c("STUSAB", "yr_wk")],
                               I_over_POSITIVE_unadj=sqrt(test_tallies_wk$state_population/test_tallies_wk$TOTAL)), all.x=TRUE)
# State-week totals of total and positive test counts, and population, added in to enable alternate weighting (2021-03-18)
# sgtf_weights merged in to enable alternate weighting (2021-03-19)
svy.dat = merge(svy.dat, test_tallies_wk[, c("STUSAB", "yr_wk", "POSITIVE", "TOTAL", "INFECTIONS", "sgtf_weights")], all.x=TRUE)
svy.dat$sgtf_weights[is.na(svy.dat$sgtf_weights) | !svy.dat$SGTF_UPSAMPLING] = 1
svy.dat = merge(svy.dat, incidence_by_region[, c("HHS", "yr_wk", "HHS_INCIDENCE")], all.x=TRUE)



# Convert any factor to string (useful if svy.dat saved and reloaded)
week0day1 = as.Date("2020-01-05") # First Sunday in 2020
fac2str = sapply(svy.dat, class)
fac2str = names(fac2str[fac2str=="factor"])
for (vv in fac2str) svy.dat[, vv] = as.character(svy.dat[, vv])
# Adding new lab streams (that are not handled by the lab-dependent weighting yet) [updated from just NS3 2021-04-02]
svy.dat$SOURCE = svy.dat$LAB # svy.dat variable redefinition resolves this
src.dat = subset(svy.dat, SOURCE!="OTHER" & !is.na(VARIANT) & VARIANT!="None") # VARIANT exclusions added 2021-04-21
src.dat$sgtf_weights[is.na(src.dat$sgtf_weights) | !src.dat$SGTF_UPSAMPLING] = 1 # svy.dat variable redefinition resolves this
# (Weighted) count of sequences
seq.tbl = with(src.dat, xtabs((1/sgtf_weights) ~ STUSAB + yr_wk))
for (rr in 1:nrow(src.dat)) {
  src.dat$SIMPLE_ADJ_WT[rr] = sqrt(src.dat$state_population[rr]/src.dat$TOTAL[rr]) * 
    src.dat$POSITIVE[rr]/seq.tbl[src.dat$STUSAB[rr], src.dat$yr_wk[rr]] / src.dat$sgtf_weights[rr]
  # Impute by HHS region for states with missing testing data
  if (is.na(src.dat$SIMPLE_ADJ_WT[rr])) {
    src.dat$SIMPLE_ADJ_WT[rr] = 
      src.dat$state_population[rr] * src.dat$HHS_INCIDENCE[rr] / seq.tbl[src.dat$STUSAB[rr], src.dat$yr_wk[rr]] / src.dat$sgtf_weights[rr]
  }
}
# Just to be sure:
src.dat = subset(src.dat, !is.na(SIMPLE_ADJ_WT) & SIMPLE_ADJ_WT < Inf)
src.dat$NORM_WTS = with(src.dat, SIMPLE_ADJ_WT/sum(SIMPLE_ADJ_WT, na.rm=TRUE))
src.dat$NORM_WTS = src.dat$NORM_WTS * sum(!is.na(src.dat$NORM_WTS))

src.dat$FORTNIGHT_END = as.character(week0day1 + src.dat$week%/%2 * 14 + 13)
src.dat$VARIANT = as.character(src.dat$VARIANT) # if saved as factor

svyREG = svydesign(ids=~STUSAB+SOURCE, strata=~HHS, weights= ~ SIMPLE_ADJ_WT, nest=TRUE, data=src.dat)

# Helper wrapper for svyciprop with CI: as vector (p, l, h), formatted string (pp.p (ll.l-hh.h)), or range (ll.l-hh.h)
# Example: myciprop("B.1.1.7", "USA", subset(svyREG, yr_wk=="2021-04-18"))
#          myciprop("L452R", 4, subset(svyREG, yr_wk=="2021-04-18"), mut=TRUE)  
myciprop = function(voc, geoid, svy, str=TRUE, range=FALSE, mut=FALSE, ...) {
  if (mut) {VOC = (unlist(gregexpr(voc, svy$variables$S_MUT)) != -1)} else {VOC = (svy$variables$VARIANT %in% voc)}
  res = svyciprop(~VOC, design = subset(update(svy, VOC=VOC), (STUSAB %in% geoid) | (HHS %in% geoid) | (geoid=="USA")), ...)
  res = c(res, confint(res))
  if (str) { 
    res = round( 100 * res, 1)
    if (range) {res = paste0(res[2], "-", res[3])} else {res = paste0(res[1], " (", res[2], "-", res[3], ")")}
  }
  res 
} 


####################################################################################
# Part 2: Nowcasting model
# Currently uses (normalized) survey weights, but not survey design
# Prediction intervals are conservative
####################################################################################

library(nnet)


n_top = 10 # Top by variant share that must be included in output
n_recent_weeks = 4 # Window for estimates
model_weeks = 16 # Lookback for modeling
week0day1 = get0("week0day1", ifnotfound=as.Date("2020-01-05"))
current_week = as.numeric(Sys.Date() - week0day1)%/%7
# List of variants to track (not just VOC or VOI):
voc = c("B.1.2", "B.1.1.7", "B.1.429", "B.1.596", "B.1.427", "B.1", "B.1.1.519", 
        "B.1.526", "B.1.526.1", "B.1.526.2", "P.2", "B.1.351", "B.1.525", "P.1", "P.2",
        "B.1.617", "B.1.617.1", "B.1.617.2", "B.1.617.3")
# Example of mutation sets to track
moi = c("L452R", "E484K") 
mean_generation_time = 6/7 # weeks; CDC proposed modeling scenarios 2021-03-19
# Estimated degrees of freedom = number of clusters - number of strata
dfs = aggregate(svyREG$cluster$SOURCE, list(USA_or_HHSRegion=svyREG$strata$HHS), function(xx) length(unique(xx)))
dfs = rbind(data.frame(USA_or_HHSRegion="USA", x=sum(dfs$x) - nrow(dfs)), dfs)

us_var = sort(
  prop.table(xtabs(SIMPLE_ADJ_WT ~ VARIANT, subset(src.dat, week >= current_week - n_recent_weeks & VARIANT != "None"))),
  decreasing=TRUE)
us_seq = table(subset(src.dat, week >= current_week - n_recent_weeks)$VARIANT)
us_rank = names(us_var)
all_tops = us_rank[us_rank %in% c(us_rank[1:n_top], voc)] # Ordered by national rank; for display
model_vars = us_var[all_tops] # For multinomial model

src.dat$K_US = sapply(1:nrow(src.dat), function(nn) which(c(names(model_vars), src.dat$VARIANT[nn]) == src.dat$VARIANT[nn])[1])

# National, unadjusted
get_cis = TRUE # Hessian slow!
temp_us = multinom(K_US ~ I(week - current_week), data=subset(src.dat, week >= current_week - model_weeks),
                   weights=NORM_WTS, Hess=get_cis, maxit=1000, trace=FALSE)

# Adjusted by HHS region
get_cis = FALSE # Hessian slow!
temp_hhs = multinom(K_US ~ I(week - current_week) + as.factor(HHS), data=subset(src.dat, week >= current_week - model_weeks), 
                    weights=NORM_WTS, Hess=get_cis, maxit=1000, trace=FALSE)

# Helper function to summarize multinomial model
# Currently needs multinomial regression output from model in one of two forms generated here
# Example: ests = se.multinom(temp_us, current_week)
se.multinom = function(mlm, week, geoid="USA") { 
  # mlm = model output, with Hessian; geoid can be HHS Region (1:10)
  # If mlm is without Hessian, all se's are set to zero
  # mlm not geographically stratified for geoid="USA": ~ (week - current_week)
  # mlm, for HHS Regions: ~ (week - current_week) + HHS (1 as reference level)
  cf = coefficients(mlm)
  if ("Hessian" %in% names(mlm)) {vc = solve(mlm$Hessian)} else {vc=rep(0, length(cf)^2)}
  dim(vc) = rep(rev(dim(cf)), 2) # Rearrange terms to ease pulling out relevant submatrix
  if (geoid==1 | geoid=="USA") {hhs1tail = NULL} else {hhs1tail = geoid + 1}
  indices = c(1, 2, hhs1tail)
  if (geoid==1 | geoid=="USA") {hhs1tail = NULL} else {hhs1tail = 1}
  coeffs = c(1, week - current_week, hhs1tail)
  # b is vector of coefficients of time (week)
  sub.vc = vc[indices,,indices,] 
  y_i = c(0, cf[, indices] %*% coeffs)
  vc.y_i = outer(1:dim(sub.vc)[2], 1:dim(sub.vc)[4], Vectorize(function(i, j) c(coeffs %*% sub.vc[, i, , j] %*% coeffs)))
  vc.y_i = rbind(0, cbind(0, vc.y_i))
  p_i = exp(y_i)/sum(exp(y_i))
  # Taylor series based variance:
  # dp_i/dy_j = delta_ij p_i - p_i * p_j
  dp_dy = diag(p_i) - outer(p_i, p_i, `*`)
  se.p_i = as.vector(sqrt(diag(dp_dy %*% vc.y_i %*% dp_dy)))
  list(p_i=p_i, se.p_i=se.p_i, b_i=c(0, cf[, 2]), se.b_i= c(0, sqrt(diag(vc[2,,2,]))))
}

#########################################################
## 2021-04-28 Prabasaj Paul
## Binomial mixture
## Beta distribution beta(al, be) for p assumed
## Inputs: mu = mean p, sigma = standard deviation of p, N = size
## mu = alpha/(al + be)
## sigma^2 = al * be /((al + be)^2 * (al + be + 1))
#########################################################

dmbinom = function(mu, sigma, size, n) {
  al_plus_be = mu * (1 - mu) / sigma^2 - 1
  al = mu * al_plus_be
  be = al_plus_be - al
  choose(size, n) * exp(lbeta(al + n, be + size - n)  - lbeta(al, be))
}

qmbinom = function(mu, sigma, size, p) {
  if (is.na(sigma) | sigma==0) {qbinom(p, size, mu)} else {findInterval(p, cumsum(dmbinom(mu, sigma, size, 0:size)))}
}

#######################################################

# Example of use
ests = se.multinom(temp_us, current_week, "USA") # Alt call: se.multinom(temp_hhs, current_week, 4) 
df = dfs[dfs[,1]=="USA", 2] # Alt: df = dfs[dfs[,1]==4, 2]
ests = data.frame(Variant=c(names(model_vars), "Other"), Proportion=ests$p_i, se.Proportion=ests$se.p_i)
ci = apply(ests, 1, function(rr) qmbinom(as.numeric(rr["Proportion"]), as.numeric(rr["se.Proportion"]), df, c(0.025, 0.975))/df)
ests$Proportion_lo = ci[1,] 
ests$Proportion_hi = pmax(1/df, ci[2,])
