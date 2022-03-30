*-------------------------------------------------------------------------------
* Stata code for replicating Cook (2014)
*-------------------------------------------------------------------------------

*-------------------------------------------------------------------------------
*** Time series regressions
*-------------------------------------------------------------------------------

*-------------------------------------------------------------------------------
*** plot application-year and grant-year on the same graph

use pats_time_series, clear
merge m:1 year using pats_fig2_fig3
drop if year<1870
drop _merge
tsset race year

set scheme plotplain
lab var patgrntpc "Patents by grant year"
lab var pat_appyear_pm "Patents by application year"
tw line pat_appyear_pm patgrntpc year if race==1, title("Black patents per million") xtitle("") legend(pos(6))
* note: pat_appyear_pm is patents for blacks; white patent-by-application-year data is not included

su pat_appyear_pm if race==1
su patgrntpc if race==1
* 1.22 for app, 0.16 for grant

* other variables are same across datasets:
tw line lynchpc lynch_pm year if race==1
tw line riot majorriot year if race==1


*-------------------------------------------------------------------------------
*** Footnote 2
su patgrntpc if race==0
su patgrntpc if race==1
* average value is 425 for whites, 0.16 for Blacks
  * conclusion: patgrntpc is actually a 'per million' variable, not 'per capita'

*-------------------------------------------------------------------------------
*** Table 6 - time series regressions
* robustness check using application-year patents
* (continue with same merged dataset)

gen lpatapppc = log(pat_appyear_pm)

lab var Dllynchpc "Lynchings"
lab var seglaw "Segregation laws"

* patents by grant year, replicating table 6
qui xtscc D.lpatgrntpc riot Dllynchpc DLMRindex  seglaw race  _iyear_1921 INTrace1921  t  DG1899 _iyear_1910 _iyear_1913 _iyear_1928, pool
est sto rep1
qui reg D.lpatgrntpc Dllynchpc riot seglaw DLMRindex _iyear_1921  t DG1899 _iyear_1910 _iyear_1913 _iyear_1928 if race==0, robust
est sto rep2
qui reg D.lpatgrntpc Dllynchpc riot seglaw DLMRindex _iyear_1921  t DG1899 _iyear_1910 _iyear_1913 _iyear_1928 if race==1, robust
est sto rep3

esttab rep1 rep2 rep3, mtitle("Full" "Whites" "Blacks") label replace order(Dllynchpc riot seglaw) keep(riot Dllynchpc seglaw race _iyear_1921 INTrace1921) se r2 star(* 0.1 ** 0.05 *** 0.01) b(%9.3f) title("Patents: grant year")

* patents by grant year, excluding 1940 to match application-year
  * footnote 4
qui xtscc D.lpatgrntpc riot Dllynchpc DLMRindex  seglaw race  _iyear_1921 INTrace1921  t  DG1899 _iyear_1910 _iyear_1913 _iyear_1928  if year<1940, pool
est sto rep_ex1
qui reg D.lpatgrntpc Dllynchpc riot seglaw DLMRindex _iyear_1921  t DG1899 _iyear_1910 _iyear_1913 _iyear_1928 if race==0 & year<1940, robust
est sto rep_ex2
qui reg D.lpatgrntpc Dllynchpc riot seglaw DLMRindex _iyear_1921  t DG1899 _iyear_1910 _iyear_1913 _iyear_1928 if race==1 & year<1940, robust
est sto rep_ex3

esttab rep_ex1 rep_ex2 rep_ex3, mtitle("Full" "Whites" "Blacks") label replace order(Dllynchpc riot seglaw) keep(riot Dllynchpc seglaw race _iyear_1921 INTrace1921) se r2 star(* 0.1 ** 0.05 *** 0.01) b(%9.3f) title("Patents: grant year, excluding 1940")

* patents by application year
qui xtscc D.lpatapppc riot Dllynchpc DLMRindex  seglaw race  _iyear_1921 INTrace1921  t  DG1899 _iyear_1910 _iyear_1913 _iyear_1928 if year<1940, pool
est sto rob1
qui reg D.lpatapppc Dllynchpc riot seglaw DLMRindex _iyear_1921  t DG1899 _iyear_1910 _iyear_1913 _iyear_1928 if race==0 & year<1940, robust
est sto rob2
qui reg D.lpatapppc Dllynchpc riot seglaw DLMRindex _iyear_1921  t DG1899 _iyear_1910 _iyear_1913 _iyear_1928 if race==1 & year<1940, robust
est sto rob3

esttab rob1 rob2 rob3, mtitle("Full" "Whites" "Blacks") label replace order(Dllynchpc riot seglaw) keep(riot Dllynchpc seglaw race _iyear_1921 INTrace1921) se r2 star(* 0.1 ** 0.05 *** 0.01) b(%9.3f) title("Patents: application year")

*-------------------------------------------------------------------------------
*** Panel data regressions
*-------------------------------------------------------------------------------

*-------------------------------------------------------------------------------
* calculating completely balanced panel
* https://en.wikipedia.org/wiki/List_of_U.S._states_by_date_of_admission_to_the_Union
* https://www.statista.com/statistics/1043617/number-us-states-by-year/

use pats_state_regs_wcontrol, clear
labelbook

* don't have state names, so need to manually copy label values
  * order them by entry year (year admitted into US)

/*
pre-1870: 38 including DC
8 DE
39 PA
31 NJ
11 GA
7 CT
22 MA
21 MD
41 SC
30 NH
47 VA
33 NY
34 NC
40 RI
46 VT
18 KY
43 TN
36 OH
19 LA
15 IN
25 MS
14 IL
1 AL
20 ME
26 MO
4 AR
23 MI
10 FL
44 TX
16 IA
50 WI
5 CA
24 MN
38 OR
17 KS
49 WV
29 NV
28 NE
9 DC
* Cook has DC from 1871

post-1870
1876:
6 CO

1889:
35 ND
42 SD
27 MT
48 WA

1890:
13 ID
51 WY

1896:
45 UT

1907:
37 OK

1912:
32 NM
3 AZ

1959:
2 AK
12 HI
*/

preserve
collapse stateno, by(year)

gen complete_panel = 38 if year<1876
replace complete_panel = 39 if year>=1876 & year<1889
replace complete_panel = 43 if year>=1889 & year<1890
replace complete_panel = 45 if year>=1890 & year<1896
replace complete_panel = 46 if year>=1896 & year<1907
replace complete_panel = 47 if year>=1907 & year<1912
replace complete_panel = 49 if year>=1912

collapse (sum) complete_panel
su
* 3210
restore

* complete panel sample size
* (1876-1870)*38 + (1889-1876)*39 + (1890-1889)*43 + (1896-1890)*45 + (1907-1896)*46 + (1912-1907)*47 + (1941-1912)*49


*-------------------------------------------------------------------------------
*** plot number of observations by state and year

use pats_state_regs_AAonly, clear

preserve
collapse (count) count=year, by(stateno)
scatter count stateno, xtitle("State ID") ytitle("") title("Observations by state")
restore


preserve
collapse (count) count=stateno, by(year)
scatter count year, xtitle("") ytitle("") title("Observations by year")
restore

*-------------------------------------------------------------------------------
*** plot number of observations by region, for actual data and complete and balanced panel

* calculating complete and balanced panel by region
use pats_state_regs_AAonly, clear

* fix errors
* state 9 has regmatl=0.33 and regs=1 in 1888; regs=1 for all other years
replace regmatl=0 if year==1888 & stateno==9
* state 14 is regmw=1, except for 1886 when it's regmw=0.5 and regs=0.5
replace regs=0 if year==1886 & stateno==14
replace regmw=1 if year==1886 & stateno==14


gen region = .
replace region = 1 if (regs==1)
replace region = 2 if (regmw==1)
replace region = 3 if (regne==1)
replace region = 4 if (regw==1)
replace region = 5 if (regmatl==1)
label define reg_label 1 "South" 2 "Midwest" 3 "Northeast" 4 "West" 5 "Mid-Atlantic"
label values region reg_label

* from wikipedia, linked above
gen entry_year = .
replace entry_year = 1870 if stateno==8 | stateno==39 | stateno==31 | stateno==11 | stateno==7 | stateno==22 | stateno==21 | stateno==41 | stateno==30 | stateno==47 | stateno==33 | stateno==34 | stateno==40 | stateno==46 | stateno==18 | stateno==43 | stateno==36 | stateno==19 | stateno==15 | stateno==25 | stateno==14 | stateno==1 | stateno==20 | stateno==26 | stateno==4 | stateno==23 | stateno==10 | stateno==44 | stateno==16 | stateno==50 | stateno==5 | stateno==24 | stateno==38 | stateno==17 | stateno==49 | stateno==29 | stateno==28 | stateno==9
replace entry_year = 1876 if stateno==6
replace entry_year = 1889 if stateno==35 | stateno==42 | stateno==27 | stateno==48
replace entry_year = 1890 if stateno==13 | stateno==51
replace entry_year = 1896 if stateno==45
replace entry_year = 1907 if stateno==37
replace entry_year = 1912 if stateno==32 | stateno==3

gen duration = 1941-entry_year

preserve
collapse year, by(stateno region)

collapse (sum) duration, by(region)
bro
* numbers for complete_balanced
restore

preserve
collapse (count) count=stateno, by(region)

* manually grab numbers from browse above
gen complete_balanced = 1028 if region==1
replace complete_balanced = 833 if region==2
replace complete_balanced = 426 if region==3
replace complete_balanced = 639 if region==4
replace complete_balanced = 426 if region==5

set scheme plotplain
lab var count "Actual data"
lab var complete_balanced "Balanced panel"

tw (scatter count region, xlabel(1 "South" 2 "Midwest" 3 "Northeast" 4 "West" 5 "Mid-Atlantic") xtitle("") title("Observations by region") ytitle("Actual data", axis(1)) yaxis(1)) (scatter complete_balanced region, yaxis(2) ytitle("Balanced panel", axis(2))), legend(pos(6))

restore

*-------------------------------------------------------------------------------
*** Compare panel data to time series

use pats_time_series, clear
collapse (sum) riot seglaw if race==1
su
* 35 riots, 290 seglaws

use pats_state_regs_AAonly, clear

collapse (sum) riot seglaw patent
su
* 5 riots, 19.33 seglaws, 702 patents

*-------------------------------------------------------------------------------
*** Heterogeneous effects by region, for Column 1, Table 7

use pats_state_regs_AAonly, clear

lab var lynchrevpc "Lynchings, per 100,000"
lab var riot "Major riots"
lab var seglaw "Segregation laws"
lab var illit "Illiteracy rate"
lab var ind "Industry participation rate"

gen region = .
replace region = 1 if (regs==1)
replace region = 2 if (regmw==1)
replace region = 3 if (regne==1)
replace region = 4 if (regw==1)
replace region = 5 if (regmatl==1)
label define reg_label 1 "South" 2 "Midwest" 3 "Northeast" 4 "West" 5 "Mid-Atlantic"
label values region reg_label

* heterogeneous effects: column 1
est clear
qui xtreg patent lynchrevpc riot seglaw illit  blksh  year1910 year1913 year1928 if region==1, re  vce(cl stateno)
est sto h1
qui xtreg patent lynchrevpc riot seglaw illit  blksh  year1910 year1913 year1928 if region==2, re  vce(cl stateno)
est sto h2
qui xtreg patent lynchrevpc riot seglaw illit  blksh  year1910 year1913 year1928 if region==3, re  vce(cl stateno)
est sto h3
qui xtreg patent lynchrevpc riot seglaw illit  blksh  year1910 year1913 year1928 if region==4, re  vce(cl stateno)
est sto h4
qui xtreg patent lynchrevpc riot seglaw illit  blksh  year1910 year1913 year1928 if region==5, re  vce(cl stateno)
est sto h5
esttab h1 h2 h3 h4 h5, mtitle("South" "Midwest" "Northeast" "West" "Mid-Atlantic") label replace order(lynchrevpc riot seglaw illit ind) keep(lynchrevpc riot seglaw illit) se star(* 0.1 ** 0.05 *** 0.01) b(%9.3f) title("Patents: grant year")

* midatlantic lynching rate is tiny
table region, c(mean lynchrevpc)
su lynchrevpc

* breakdown of riots, seglaws, and patents by region
table region, c(sum riot sum seglaw sum patent)

*-------------------------------------------------------------------------------
*** Table 8

qui xtreg assn lynchrevpc riot seglaw illit blksh ind regs regmw regne regw  grinvent year1910 year1913 year1928, re vce(cl stateno)
est sto t81
qui xtreg mech lynchrevpc riot seglaw illit blksh ind regs regmw regne regw year1910 year1913 year1928, re vce(cl stateno)
est sto t82
qui xtreg elec lynchrevpc riot seglaw illit blksh ind regs regmw regne regw year1910 year1913 year1928, re vce(cl stateno)
est sto t83
qui xtreg patsth lynchrevpc riot seglaw illit blksh ind regs regmw regne regw  year1910 year1913 year1928, re vce(cl stateno)
est sto t84
esttab t81 t82 t83 t84, mtitle("Assigned" "Mechanical" "Electrical" "Southern") label replace order(lynchrevpc riot seglaw illit ind) keep(lynchrevpc riot seglaw illit grinvent) se star(* 0.1 ** 0.05 *** 0.01) b(%9.3f) title("Patents: grant year")
* N is off by one, get slightly different results

*-------------------------------------------------------------------------------
*** Table 9
* Cook's code outputs the negative binomial regression results, but Table 9 reports the average marginal effects
* the code below reports the marginal effects
* the results are very different from Table 9 in the paper

use pats_state_regs_wcontrol.dta, clear
collapse (sum) patent assn mech elec patsth (mean)  lynchpc2  riot seglaw2 illit regmatl regs regmw regne regw estbnumpc, by(stateno race)

***Blacks
preserve

keep if race==1

nbreg patent lynchpc2  riot seglaw2 illit regmw regne regs regw estbnumpc, robust nolog
eststo margin: margins, dydx(lynchpc2 riot seglaw2 illit estbnumpc) post
estimates store r1cw

nbreg assn lynchpc2  riot seglaw2 illit regmw regne regs regw estbnumpc, robust nolog
eststo margin: margins, dydx(lynchpc2 riot seglaw2 illit estbnumpc) post
estimates store r1cx

nbreg mech lynchpc2  riot seglaw2 illit regmw regne regs regw estbnumpc, robust nolog
eststo margin: margins, dydx(lynchpc2 riot seglaw2 illit estbnumpc) post
estimates store r1cy

nbreg elec lynchpc2  riot seglaw2 illit regmw regne regs regw estbnumpc, robust nolog
eststo margin: margins, dydx(lynchpc2 riot seglaw2 illit estbnumpc) post
estimates store r1cz

nbreg patsth lynchpc2  riot seglaw2 illit estbnumpc, robust nolog
eststo margin: margins, dydx(lynchpc2 riot seglaw2 illit estbnumpc) post
estimates store r1caa

esttab r1cw r1cx r1cy r1cz r1caa, nomtitle label replace order(lynchpc2 riot seglaw2 illit estbnumpc) keep(lynchpc2 riot seglaw2 illit estbnumpc) se star(* 0.1 ** 0.05 *** 0.01) b(%9.3f) title("Black patents")

restore

*** Whites
preserve

keep if race==0

nbreg patent lynchpc2  riot seglaw2 illit regmw regne regs regw estbnumpc, robust nolog
eststo margin: margins, dydx(lynchpc2 riot seglaw2 illit estbnumpc) post
estimates store r0cw

nbreg assn lynchpc2  riot seglaw2 illit regmw regne regs regw estbnumpc, robust nolog
eststo margin: margins, dydx(lynchpc2 riot seglaw2 illit estbnumpc) post
estimates store r0cx

nbreg mech lynchpc2  riot seglaw2 illit regmw regne regs regw estbnumpc, robust nolog
eststo margin: margins, dydx(lynchpc2 riot seglaw2 illit estbnumpc) post
estimates store r0cy

nbreg elec lynchpc2  riot seglaw2 illit regmw regne regs regw estbnumpc, robust nolog
eststo margin: margins, dydx(lynchpc2 riot seglaw2 illit estbnumpc) post
estimates store r0cz

nbreg patsth lynchpc2  riot seglaw2 illit estbnumpc, robust nolog
eststo margin: margins, dydx(lynchpc2 riot seglaw2 illit estbnumpc) post
estimates store r0caa

esttab r0cw r0cx r0cy r0cz r0caa, nomtitle label replace order(lynchpc2 riot seglaw2 illit estbnumpc) keep(lynchpc2 riot seglaw2 illit estbnumpc) se star(* 0.1 ** 0.05 *** 0.01) b(%9.3f) title("Black patents")

restore
