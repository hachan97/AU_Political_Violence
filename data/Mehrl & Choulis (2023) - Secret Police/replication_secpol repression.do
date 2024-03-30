*****************************************************************************************
** Replication Data
** Marius Mehrl & ioannis Choulis - "Secret Police Organizations and State Repression"
** Journal of Conflict Resolution
*****************************************************************************************

* load data (Change working directory!)

use "C:\...\secret police repression\secret police repression_replication data.dta", clear

keep milex pop gdp_pc effective nonvio_protests year electionyear ccode secretpol_ peergroup_revised peer_countrygroups gwf_party gwf_military gwf_monarchy gwf_personal gwf_nonautocracy theta_mean regtype intrastate interstate attempt cbcount xpers lexclpop PTS_S physint electionyear solschdum

* last variable constructions

gen ln_milex = ln(milex+1)
gen ln_pop = ln(pop+1)
gen ln_gdppc = ln(gdp_pc+1)
gen ln_counterb = ln(effective+1)

label variable peergroup_revised "Secret Polices in peer group countrys"
label variable peer_countrygroups "Peer groupings"

xtset ccode year

*************************************
* Main Analyses (reported in paper) *
*************************************

* Figure 1 - create source data in Stata
preserve

collapse (min) minsp=secretpol_revised (max) maxsp=secretpol_revised, by(ccode)
gen spexist=0
replace spexist=1 if minsp==0 & maxsp==0
replace spexist=2 if minsp==0 & maxsp==1
replace spexist=3 if minsp==1 & maxsp==1
replace spexist=2 if ccode==255
keep ccode spexist
replace ccode = 316 if ccode==315
replace ccode = 260 if ccode==255
replace ccode = 340 if ccode==345
replace ccode = 678 if ccode==679

export delimited using "C:\...\secret police repression\sp_mapsource.csv", replace

restore

* Figure 1 - Load source data and create map in R

* !! SWITCH TO R !!

library(cshapes)
library(tidyverse)
library(classInt)
library(RColorBrewer)

setwd("C:\\...\\secret police repression")

secpol <- read.csv("sp_mapsource.csv")
world <- cshp(date = as.Date("2016-1-1"))

world2 <- left_join(world, secpol, by = c("gwcode"="ccode"))

# generate a color palette 
pal <- brewer.pal(4, "Accent")
# find the class intervals and colors> 
breaks <- classIntervals(world2$spexist, n=4, style="fixed", fixedBreaks=c(-0, 1, 2, 3))
colors <- findColours(breaks, pal)
# create plot and add legend> 

plot(world2["spexist"], main = "", bty="n", col = colors)
par(xpd=TRUE)
legend(x=-0.1, y=0.15, legend=c("Unobserved", "Never", "Variation", "Always"), 
       fill = attr(colors, "palette"), 
       bty = "n", title="Secret Police", ncol=2, cex=0.6)


tiff("spmap.tiff", units="cm", width=16, height=8, res=800)
plot(world2["spexist"], main = "", bty="n", col = colors)
par(xpd=TRUE)
legend(x=-0.1, y=0.15, legend=c("Unobserved", "Never", "Variation", "Always"), 
                                          fill = attr(colors, "palette"), 
       bty = "n", title="Secret Police", ncol=2, cex=0.6)
dev.off()

* !! SWITCH BACK TO STATA !!

* Table 1
gen gwf_categ=.
replace gwf_categ=0 if gwf_party==1
replace gwf_categ=1 if gwf_military==1
replace gwf_categ=2 if gwf_monarchy==1
replace gwf_categ=3 if gwf_personal==1
replace gwf_categ=4 if gwf_categ==. & secretpol_revised!=.
tab secretpol_revised gwf_categ if year<=2010, column

* Figure 2
graph box theta_mean, over(secretpol_, relabel(1 "No Secret Police" 2 "Secret Police")) scheme(tufte) aspect(1) ytitle("Human Rights Score") 


* Table 2
reg theta_mean secretpol_, vce(cluster ccode)
reghdfe theta_mean secretpol_, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc i.regtype, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc  i.regtype intrastate interstate attempt, absorb(ccode) vce(cluster ccode)

* Figure 3
gen secpol_turnedaround=secretpol_
replace secpol_turnedaround=2 if secretpol_==0
replace secpol_turnedaround=secpol_turnedaround-1
tab secpol_turnedaround secretpol_
btscs secpol_turnedaround year ccode, gen(secpolduration)
tab secpolduration

tab cbcount

reghdfe theta_mean i.secretpol_##c.secpolduration##c.secpolduration##c.secpolduration ln_pop ln_gdppc  i.regtype intrastate interstate attempt, absorb(ccode) vce(cluster ccode)
margins, dydx(secretpol_) at(secpolduration=(0(1)45))
marginsplot, scheme(plotplain) yline(0) xlabel(#10) xtitle("Secret Police Duration") ytitle("Change in Human Rights Score") title(" ") name(g1, replace)

reghdfe theta_mean i.secretpol_##c.cbcount##c.cbcount ln_pop ln_gdppc  i.regtype intrastate interstate attempt, absorb(ccode) vce(cluster ccode)
margins, dydx(secretpol_) at(cbcount=(0(1)7))
marginsplot, scheme(plotplain) yline(0) xtitle("Number of Counterweights") ytitle(" ") title(" ") name(g2, replace)

graph combine g1 g2, scheme(plotplain) ycommon

************************************************************
* Additional Analyses (reported in supplementary material) *
************************************************************

* Table A1
quietly reghdfe theta_mean secretpol_ ln_pop ln_gdppc  i.regtype intrastate interstate attempt, absorb(ccode) vce(cluster ccode)

sum theta_mean secretpol_ ln_pop ln_gdppc  i.regtype intrastate interstate attempt secpolduration cbcount peergroup_revised  gwf_military gwf_party gwf_personal xpers nonvio lexclpop ln_milex PTS_S physint electionyear if secretpol_!=. & e(sample)==1

* Tables A3 & A4 - IV models
* Construct average level of human rights violations in peergroup
gen one=1
bysort peer_countrygroups year: egen theta_peersum = sum(theta_mean)
bysort peer_countrygroups year: egen peer_number = sum(one)
gen dropown= theta_peersum-theta_mean
gen theta_peeraverage=(dropown/(peer_number-1))
sum theta_peeraverage
drop one theta_peersum peer_number dropown
xtset ccode year

* Construct average security apparatus fragmentation in peer group
gen one=1
bysort peer_countrygroups year: egen cbcount_peersum = sum(cbcount)
bysort peer_countrygroups year: egen peer_number = sum(one)
gen dropown= cbcount_peersum-cbcount
gen cbcount_peeraverage=(dropown/(peer_number-1))
sum cbcount_peeraverage
drop one cbcount_peersum peer_number dropown
xtset ccode year

sum theta_peeraverage cbcount_peeraverage

* Models 1 & 6
xtlogit secretpol_ peergroup_revised, fe
predict prob0
xtivreg theta_mean  (secretpol_ = prob0), fe vce(cluster ccode)
tab ccode if e(sample)==1, gen(countrydummy)
ivreg2 theta_mean countrydummy2-countrydummy115 (secretpol_ = prob0), cluster(ccode) first
drop prob0 countrydummy*

* Models 2 & 7
xtlogit secretpol_ peergroup_revised ln_pop ln_gdppc i.regtype, fe
predict prob1
xtivreg theta_mean  ln_pop ln_gdppc i.regtype  (secretpol_ = prob1), fe vce(cluster ccode)
tab ccode if e(sample)==1, gen(countrydummy)
ivreg2 theta_mean ln_pop ln_gdppc regtype countrydummy2-countrydummy115 (secretpol_ = prob1), cluster(ccode) first
drop prob1 countrydummy*

* Models 3 & 8
xtlogit secretpol_ peergroup_revised ln_pop ln_gdppc i.regtype intrastate interstate attempt, fe
predict prob2
xtivreg theta_mean  ln_pop ln_gdppc i.regtype intrastate interstate attempt (secretpol_ = prob2), fe vce(cluster ccode)
tab ccode if e(sample)==1, gen(countrydummy)
ivreg2 theta_mean ln_pop ln_gdppc regtype intrastate interstate attempt countrydummy2-countrydummy115 (secretpol_ = prob2), cluster(ccode) first
drop prob2 countrydummy*

* Models 4 & 9
xtlogit secretpol_ peergroup_revised ln_pop ln_gdppc i.regtype intrastate interstate attempt l.theta_peeraverage, fe
predict prob3
xtivreg theta_mean  ln_pop ln_gdppc i.regtype intrastate interstate attempt l.theta_peeraverage (secretpol_ = prob3), fe vce(cluster ccode)
tab ccode if e(sample)==1, gen(countrydummy)
ivreg2 theta_mean ln_pop ln_gdppc regtype intrastate interstate attempt l.theta_peeraverage countrydummy2-countrydummy115 (secretpol_ = prob3), cluster(ccode) first
drop prob3 countrydummy*

* Models 5 & 10
xtlogit secretpol_ peergroup_revised ln_pop ln_gdppc i.regtype intrastate interstate attempt l.cbcount_peeraverage, fe
predict prob4
xtivreg theta_mean  ln_pop ln_gdppc i.regtype intrastate interstate attempt l.cbcount_peeraverage (secretpol_ = prob4), fe vce(cluster ccode)
tab ccode if e(sample)==1, gen(countrydummy)
ivreg2 theta_mean ln_pop ln_gdppc regtype intrastate interstate attempt cbcount l.cbcount_peeraverage countrydummy2-countrydummy80 (secretpol_ = prob4), cluster(ccode) first
drop prob4 countrydummy*


* Table A5 - Year-FE
reghdfe theta_mean secretpol_, absorb(ccode year) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc i.regtype, absorb(ccode year) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc  i.regtype intrastate interstate attempt, absorb(ccode year) vce(cluster ccode)


* Table A6 - lagged DV
reghdfe theta_mean secretpol_ l.theta_mean, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean secretpol_ l.theta_mean ln_pop ln_gdppc i.regtype, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean secretpol_ l.theta_mean ln_pop ln_gdppc i.regtype intrastate interstate attempt, absorb(ccode) vce(cluster ccode)


* Table A7 - Cubic country-trends
gen year_sq=year*year
gen year_cb=year*year_sq

reghdfe theta_mean secretpol_, absorb(i.ccode##c.year i.ccode#c.year_sq i.ccode#c.year_cb) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc i.regtype, absorb(i.ccode##c.year i.ccode#c.year_sq i.ccode#c.year_cb) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc i.regtype intrastate interstate attempt, absorb(i.ccode##c.year i.ccode#c.year_sq i.ccode#c.year_cb) vce(cluster ccode)


* Table A8 - GWF authoritarianism variables
reghdfe theta_mean secretpol_ gwf_military gwf_party gwf_personal, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc gwf_military gwf_party gwf_personal, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc gwf_military gwf_party gwf_personal intrastate interstate attempt, absorb(ccode) vce(cluster ccode)


* Table A9 - GWF Personalization Score
reghdfe theta_mean secretpol_ xpers, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc xpers, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc intrastate interstate attempt xpers, absorb(ccode) vce(cluster ccode)

* Table A10 - Additional controls: NAVCO non-violent protest, Military Expenditures, Counterbalancing, election years 
reghdfe theta_mean secretpol_ ln_pop ln_gdppc  i.regtype intrastate interstate attempt nonvio, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc  i.regtype intrastate interstate attempt lexclpop, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc  i.regtype intrastate interstate attempt ln_milex, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc  i.regtype intrastate interstate attempt c.cbcount##c.cbcount, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean secretpol_ ln_pop ln_gdppc  i.regtype intrastate interstate attempt electionyear, absorb(ccode) vce(cluster ccode)


* Figure A1 - Sensitivitiy to unobserved confounders (Plot color needs to be made white in Graph editor)
sensemakr theta_mean secretpol_ ln_pop ln_gdppc i.regtype intrastate interstate attempt i.ccode, treat(secretpol_) benchmark(intrastate) contourplot kd(1 1.5 2 2.5 3 3.5 4)

* Table A11 - Random Effects Models
xtreg theta_mean secretpol_, re vce(cluster ccode)
xtreg theta_mean secretpol_ ln_pop ln_gdppc i.regtype, re vce(cluster ccode)
xtreg theta_mean secretpol_ ln_pop ln_gdppc  i.regtype intrastate interstate attempt, re vce(cluster ccode)

* Table A12 - Political Terror Scale instead of Farriss Measure (State Department Measure as it has most obs.)
ologit PTS_S secretpol_, vce(cluster ccode)
xtologit PTS_S secretpol_, vce(cluster ccode)
xtologit PTS_S secretpol_ ln_pop ln_gdppc i.regtype, vce(cluster ccode)
xtologit PTS_S secretpol_ ln_pop ln_gdppc  i.regtype intrastate interstate attempt, vce(cluster ccode)

* Table A13 - CIRI Physical Integrity Indicator instead of Farriss Measure 
tab physint
ologit physint secretpol_, vce(cluster ccode)
xtologit physint secretpol_, vce(cluster ccode)
xtologit physint secretpol_ ln_pop ln_gdppc i.regtype, vce(cluster ccode)
xtologit physint secretpol_ ln_pop ln_gdppc  i.regtype intrastate interstate attempt, vce(cluster ccode)

* Table A14 & A15 - Interaction coefficient tables
btscs solschdum year ccode, gen(regduration)
tab regduration

reghdfe theta_mean i.secretpol_##c.secpolduration##c.secpolduration##c.secpolduration, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean i.secretpol_##c.secpolduration##c.secpolduration##c.secpolduration ln_pop ln_gdppc i.regtype, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean i.secretpol_##c.secpolduration##c.secpolduration##c.secpolduration ln_pop ln_gdppc  i.regtype intrastate interstate attempt, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean i.secretpol_##c.secpolduration##c.secpolduration##c.secpolduration ln_pop ln_gdppc  i.regtype intrastate interstate attempt regduration, absorb(ccode) vce(cluster ccode)

reghdfe theta_mean i.secretpol_##c.cbcount##c.cbcount, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean i.secretpol_##c.cbcount##c.cbcount ln_pop ln_gdppc  i.regtype, absorb(ccode) vce(cluster ccode)
reghdfe theta_mean i.secretpol_##c.cbcount##c.cbcount ln_pop ln_gdppc  i.regtype intrastate interstate attempt, absorb(ccode) vce(cluster ccode)
margins, dydx(secretpol_) at(cbcount=(0(1)7))
marginsplot, scheme(plotplain) yline(0) xtitle("Number of Counterweights") ytitle(" ") title(" ") name(g2, replace)


* Table A16 & Figure A2 - Interaction models controlling for protests and election years
reghdfe theta_mean i.secretpol_##c.secpolduration##c.secpolduration##c.secpolduration ln_pop ln_gdppc  i.regtype intrastate interstate attempt nonvio, absorb(ccode) vce(cluster ccode)
margins, dydx(secretpol_) at(secpolduration=(0(1)45))
marginsplot, scheme(plotplain) yline(0) xlabel(#10) xtitle("Control: Protest") ytitle("Change in Human Rights Score") title(" ") name(gr1, replace)

reghdfe theta_mean i.secretpol_##c.secpolduration##c.secpolduration##c.secpolduration ln_pop ln_gdppc  i.regtype intrastate interstate attempt electionyear, absorb(ccode) vce(cluster ccode)
margins, dydx(secretpol_) at(secpolduration=(0(1)45))
marginsplot, scheme(plotplain) yline(0) xlabel(#10) xtitle("Control: Election Year") ytitle(" ") title(" ") name(gr2, replace)

reghdfe theta_mean i.secretpol_##c.cbcount##c.cbcount ln_pop ln_gdppc  i.regtype intrastate interstate attempt nonvio, absorb(ccode) vce(cluster ccode)
margins, dydx(secretpol_) at(cbcount=(0(1)7))
marginsplot, scheme(plotplain) yline(0) xtitle("Control: Protest") title(" ") ytitle("Change in Human Rights Score") name(gr3, replace)

reghdfe theta_mean i.secretpol_##c.cbcount##c.cbcount ln_pop ln_gdppc  i.regtype intrastate interstate attempt electionyear, absorb(ccode) vce(cluster ccode)
margins, dydx(secretpol_) at(cbcount=(0(1)7))
marginsplot, scheme(plotplain) yline(0) xtitle("Control: Election Year") ytitle(" ") title(" ") name(gr4, replace)

graph combine gr1 gr2 gr3 gr4, scheme(plotplain) ycommon


* Figure A3 & A4 - Interactions: Interaction effect linearity
tab regtype, gen(regdummy)
interflex theta_mean secretpol_ secpolduration ln_pop ln_gdppc  regdummy2 intrastate interstate attempt, fe(ccode) cluster(ccode) type(binning) cutoffs(5 10 25 40) xdistr(density) xlabel("Secret Police Duration") ylabel("Human Rights")

interflex theta_mean secretpol_ cbcount ln_pop ln_gdppc  regdummy2 intrastate interstate attempt regduration, fe(ccode) cluster(ccode) type(binning) cutoffs(1 2 3 5) xdistr(density) xlabel("Number of Counterweights") ylabel("Human Rights")


* Table A17 - Cross-Validation
tab ccode, gen(country_fe)
gen anocracy=regtype if regtype!=2
gen lag_theta_mean = l.theta_mean

 * Model including secret police
set seed 1860
forvalues h=1/10 {
preserve

quietly reg theta_mean secretpol_ ln_pop ln_gdppc anocracy intrastate interstate attempt country_fe*, vce(cluster ccode)
xtile group=uniform() if e(sample), nq(4)
drop if group==.
gen yhat=.

forvalues i=1/4 {

quietly reg theta_mean secretpol_ ln_pop ln_gdppc anocracy intrastate interstate attempt country_fe* if group~=`i', vce(cluster ccode)
quietly predict yhat_i

quietly replace yhat=yhat_i if group==`i'
quietly drop yhat_i
}
gen suemod=(theta_mean-yhat)^2
quietly su suemod
local MSPE=r(mean)
display `MSPE'
capture drop group yhat

restore
}

 * Model excluding secret police
 set seed 1860
forvalues h=1/10 {
preserve

quietly reg theta_mean secretpol_ ln_pop ln_gdppc anocracy intrastate interstate attempt country_fe*, vce(cluster ccode)
xtile group=uniform() if e(sample), nq(4)
drop if group==.
gen yhat=.

forvalues i=1/4 {

quietly reg theta_mean  ln_pop ln_gdppc anocracy intrastate interstate attempt country_fe* if group~=`i', vce(cluster ccode)
quietly predict yhat_i

quietly replace yhat=yhat_i if group==`i'
quietly drop yhat_i
}
gen suemod=(theta_mean-yhat)^2
quietly su suemod
local MSPE=r(mean)
display `MSPE'
capture drop group yhat

restore
}

drop country_fe* anocracy lag_theta_mean
