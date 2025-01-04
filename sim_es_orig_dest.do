* simulation using size = city id + city-time FE
* staggered timing
* log size
* show diff in diff and event study
* show FWL regression for sd_p0 in event study
* calculate number of city-year observations that are used as both origin and destination
* vary number N of cities

* when a city is observed only as an origin or only as a destination, then Post is collinear with the city FE
* when some cities are origin- or destination-only, and some are both, get partial collinearity
* increasing N reduces the number of cities that are both, which increases the collinearity
	* N=2000: get perfect collinearity

clear
set seed 0

foreach N in 20 100 200 500 1000 2000 {
preserve

* Create 1000 movers
qui set obs 1000
gen id = _n

* Randomly assign origin and destination cities
gen n_city=`N'
gen city1 = ceil(n_city*runiform())
gen city2 = ceil(n_city*runiform())
* if origin and destinations are completely disjoint, then Post is collinear with city FE
// gen city1 = ceil(n_city*runiform(0,1))
// gen city2 = ceil(n_city*runiform(1,2))
qui drop if city1==city2
* drop stayers

* use staggered timing: inventors move in different years
bysort id: gen start_year = ceil(runiform()*40) + 1970

* Create 11 time periods for each mover. Move in period 6.
expand 11
bys id: gen time = _n
gen post = time >= 6
replace time = time - 6

gen year = time + start_year

* Current city
gen city = city1
qui replace city = city2 if post == 1

gegen citypost_fe_tag = tag(post city)
qui gen citypost_fe_temp = ceil(20*runiform()) if citypost_fe_tag==1
gegen citypost_fe = mean(citypost_fe_temp), by(post city)

* city size: city id + time FE
gen size = city + citypost_fe

qui gen logsize = log(size)
gegen logsize_pre_temp = mean(logsize) if post==0, by(id)
gegen logsize_pre = mean(logsize_pre_temp), by(id)
gegen logsize_post_temp = mean(logsize) if post==1, by(id)
gegen logsize_post = mean(logsize_post_temp), by(id)
gen size_diff = logsize_post - logsize_pre

local beta = 0.7
gen y_log= `beta'*logsize + rnormal(0,1)

* individual, time, and city FE
reghdfe y_log c.size_diff##i.post, a(id year city) cluster(city1)

*** event study
forvalues i=0/5 {
    gen m_`i' = (time==-`i')
    gen p_`i' = (time==`i')
    gen sd_m`i' = size_diff*m_`i'
    gen sd_p`i' = size_diff*p_`i'
    lab var sd_m`i' "-`i'"
    lab var sd_p`i' "`i'"
}

cap unab sizediff_all: sd_m* sd_p*
cap unab omit: sd_m1 sd_m0
local inds : list sizediff_all - omit

reghdfe y_log `inds' sd_m1, absorb(id year city) vce(cluster city1)
coefplot, drop(_cons) vert keep(sd_m5 sd_m4 sd_m3 sd_m2 sd_m1 sd_p0 sd_p1 sd_p2 sd_p3 sd_p4 sd_p5) order(sd_m5 sd_m4 sd_m3 sd_m2 sd_m1 sd_p0 sd_p1 sd_p2 sd_p3 sd_p4 sd_p5) omitted xtitle("Years since move") 

* FWL first stage
reghdfe sd_p0 sd_m5 sd_m4 sd_m3 sd_m2  sd_p1 sd_p2 sd_p3 sd_p4 sd_p5, absorb(year id city) vce(cluster city)

* are city-year pairs used as both origin and destination?
collapse (count) count=id, by(city year post)
reshape wide count, i(city year) j(post)
collapse (sum) orig=count0 dest=count1, by(city year)
gen both = orig>0 & dest>0
tab both
di "`N' cities"

restore

}

