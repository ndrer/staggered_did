// Install packages
ssc install blindschemes, replace
ssc install estout, replace	
ssc install ftools, replace
ssc install reghdfe, replace
net install did2s, from("https://raw.githubusercontent.com/kylebutts/did2s_stata/main/ado/") replace
ssc install did_imputation, replace
ssc install event_plot, replace
ssc install did_multiplegt, replace
net install csdid, from ("https://raw.githubusercontent.com/friosavila/csdid_drdid/main/code/") replace
net install github, from("https://haghish.github.io/github/") replace
github install lsun20/eventstudyinteract, replace
github install joshbleiberg/stackedev, replace


// Set font, graph
graph set window fontface "Roboto"
graph set window fontface default
set scheme plotplainblind
set more off
matrix drop _all

clear all
timer clear
set seed 10
global T = 15
global I = 300

set obs `=$I*$T'
gen i = int((_n-1)/$T )+1
gen t = mod((_n-1),$T )+1
tsset i t

// Treatment rollout periods Ei=10..16

gen Ei = ceil(runiform()*7)+$T -6 if t==1
bys i (t): replace Ei = Ei[1]
gen K = t-Ei
gen D = K>=0 & Ei!=.

// Generate outcome with parallel trends and heterogeneous treatment effects

gen tau = cond(D==1, (t-12.5), 0)
gen eps = rnormal()
gen Y = i + 3*t + tau*D + eps

// Generate a panel of 300 units in 15 periods

clear all
timer clear
set seed 10
global T = 15
global I = 300

set obs '=$I*$T'
gen i = int((_n-1)/$T )+1
gen t = mod((_n-1),$T )+1
tsset i t

// Randomly generate treatment rollout periods uniformly across Ei=10..16

gen Ei = ceil(runiform()*7)+$T -6 if t==1
bys i (t): replace Ei = Ei[1]
gen K = t-Ei
gen D = K>=0 & Ei!=.

// Generate the outcome with parallel trends and heterogeneous treatment effects

gen tau = cond(D==1, (t-12.5), 0)
gen eps = rnormal()
gen Y = i + 3*t + tau*D + eps

		
// Borusyak et al. (2021)

did_imputation Y i t Ei, allhorizons pretrends(5)
	
event_plot, default_look graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") title("Borusyak et al. (2021) imputation estimator") xlabel(-5(1)5) name(BJS, replace)) together

estimates store bjs
	
	
// de Chaisemartin and D`haultfoeuille (2020)

did_multiplegt Y i t D, robust_dynamic dynamic(5) placebo(5) longdiff_placebo breps(100) cluster(i)
	
event_plot e(estimates)#e(variances), default_look graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") title("de Chaisemartin and D`haultfoeuille (2020)") xlabel(-5(1)5) name(dCdH, replace)) stub_lag(Effect_#) stub_lead(Placebo_#) together

matrix dcdh_b = e(estimates)
matrix dcdh_v = e(variances)


// Callaway and Sant'Anna (2020)

gen gvar = cond(Ei>15, 0, Ei)

csdid Y, ivar(i) time(t) gvar(gvar) agg(event)

event_plot e(b)#e(V), default_look graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") xlabel(-14(1)5) title("Callaway and Sant'Anna (2020)") name(CS, replace)) stub_lag(T+#) stub_lead(T-#) together

matrix cs_b = e(b)
matrix cs_v = e(V)


// Sun and Abraham (2020)

sum Ei
gen lastcohort = Ei==r(max)
forvalues l = 0/5 {
	gen L`l'event = K==`l'
}
forvalues l = 1/14 {
	gen F`l'event = K==-`l'
}
drop F1event

eventstudyinteract Y L*event F*event, vce(cluster i) absorb(i t) cohort(Ei) control_cohort(lastcohort)

event_plot e(b_iw)#e(V_iw), default_look graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") xlabel(-14(1)5) title("Sun and Abraham (2020)") name(SA, replace)) stub_lag(L#event) stub_lead(F#event) together

matrix sa_b = e(b_iw)
matrix sa_v = e(V_iw)


// Gardner (2021)

did2s Y, first_stage(i.i i.t) second_stage(F*event L*event) treatment(D) cluster(i)

event_plot, default_look stub_lag(L#event) stub_lead(F#event) together graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") xlabel(-14(1)5) title("Gardner (2021)") name(DID2S, replace))
		
matrix did2s_b = e(b)
matrix did2s_v = e(V)

	
// Cengiz et al. (2019) 

gen treat_year=.
replace treat_year=Ei if Ei!=16

gen no_treat= (Ei==16)

cap drop F*event L*event
sum Ei
forvalues l = 0/5 {
	gen L`l'event = K==`l'
	replace L`l'event = 0 if no_treat==1
}
forvalues l = 1/14 {
	gen F`l'event = K==-`l'
	replace F`l'event = 0 if no_treat==1
}
drop F1event

preserve
stackedev Y F*event L*event, cohort(treat_year) time(t) never_treat(no_treat) unit_fe(i) clust_unit(i) 
restore

event_plot e(b)#e(V), default_look graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") xlabel(-14(1)5) title("Cengiz et al. (2019)") name(CDLZ, replace)) stub_lag(L#event) stub_lead(F#event) together	
		
matrix stackedev_b = e(b)
matrix stackedev_v = e(V)
	
// TWFE OLS estimation

reghdfe Y F*event L*event, absorb(i t) vce(cluster i)

event_plot, default_look stub_lag(L#event) stub_lead(F#event) together graph_opt(xtitle("Days since the event") ytitle("OLS coefficients") xlabel(-14(1)5) title("OLS") name(OLS, replace))

estimates store ols

	
	
// Construct vector of true average treatment effects by number of periods since treatment

matrix btrue = J(1,6,.)
matrix colnames btrue = tau0 tau1 tau2 tau3 tau4 tau5
qui forvalues h = 0/5 {
	sum tau if K==`h'
	matrix btrue[1,`h'+1]=r(mean)
}


// Combine all plots using the stored estimates (5 leads and lags around event)

event_plot btrue# bjs dcdh_b#dcdh_v cs_b#cs_v sa_b#sa_v did2s_b#did2s_v stackedev_b#stackedev_v ols, stub_lag(tau# tau# Effect_# T+# L#event L#event L#event L#event) stub_lead(pre# pre# Placebo_# T-# F#event F#event F#event F#event) plottype(scatter) ciplottype(rcap) together perturb(-0.325(0.1)0.325) trimlead(5) noautolegend graph_opt(title("Event study estimators in a simulated panel (300 units, 15 periods)", size(med)) xtitle("Periods since the event", size(small)) ytitle("Average causal effect", size(small)) xlabel(-5(1)5) legend(order(1 "True value" 2 "Borusyak et al." 4 "de Chaisemartin-D`haultfoeuille" 6 "Callaway-Sant'Anna" 8 "Sun-Abraham" 10 "Gardner" 12 "Cengiz et al." 14 "TWFE OLS") rows(2) position(6) region(style(none))) xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) graphregion(color(white)) bgcolor(white) ylabel(, angle(horizontal))) lag_opt1(msymbol(+) color(black)) lag_ci_opt1(color(black)) lag_opt2(msymbol(O) color(cranberry)) lag_ci_opt2(color(cranberry)) lag_opt3(msymbol(Dh) color(navy)) lag_ci_opt3(color(navy)) lag_opt4(msymbol(Th) color(forest_green)) lag_ci_opt4(color(forest_green)) lag_opt5(msymbol(Sh) color(dkorange)) lag_ci_opt5(color(dkorange)) lag_opt6(msymbol(Th) color(blue)) lag_ci_opt6(color(blue)) lag_opt7(msymbol(Dh) color(red)) lag_ci_opt7(color(red)) lag_opt8(msymbol(Oh) color(purple)) lag_ci_opt8(color(purple))
graph export "D:/Installer program/seven_estimators_example_5t.png", replace

// Combine all plots using the stored estimates (5 leads and lags around event)

event_plot btrue# bjs dcdh_b#dcdh_v cs_b#cs_v sa_b#sa_v did2s_b#did2s_v stackedev_b#stackedev_v ols, stub_lag(tau# tau# Effect_# T+# L#event L#event L#event L#event) stub_lead(pre# pre# Placebo_# T-# F#event F#event F#event F#event) plottype(scatter) ciplottype(rcap) together perturb(-0.325(0.1)0.325) trimlead(14) noautolegend graph_opt(title("Event study estimators in a simulated panel (300 units, 15 periods)", size(med)) xtitle("Periods since the event", size(small)) ytitle("Average causal effect", size(small)) xlabel(-14(1)5) legend(order(1 "True value" 2 "Borusyak et al." 4 "de Chaisemartin-D`haultfoeuille" 6 "Callaway-Sant'Anna" 8 "Sun-Abraham" 10 "Gardner" 12 "Cengiz et al." 14 "TWFE OLS") rows(2) position(6) region(style(none))) xline(-0.5, lcolor(gs8) lpattern(dash)) yline(0, lcolor(gs8)) graphregion(color(white)) bgcolor(white) ylabel(, angle(horizontal))) lag_opt1(msymbol(+) color(black)) lag_ci_opt1(color(black)) lag_opt2(msymbol(O) color(cranberry)) lag_ci_opt2(color(cranberry)) lag_opt3(msymbol(Dh) color(navy)) lag_ci_opt3(color(navy)) lag_opt4(msymbol(Th) color(forest_green)) lag_ci_opt4(color(forest_green)) lag_opt5(msymbol(Sh) color(dkorange)) lag_ci_opt5(color(dkorange)) lag_opt6(msymbol(Th) color(blue)) lag_ci_opt6(color(blue)) lag_opt7(msymbol(Dh) color(red)) lag_ci_opt7(color(red)) lag_opt8(msymbol(Oh) color(purple)) lag_ci_opt8(color(purple))
graph export "D:/Installer program//seven_estimators_example_allt.png", replace
