clear all
set seed 1234567

local firms 2000
local periods 100
local groups 50
local obs=`firms'*`periods'
set obs `obs'

gen n=_n
gen firm=ceil(_n/`periods') 
bys firm: gen period=_n
gen group=ceil(_n/(_N/`groups')) 

// Create program to generate sample from data-generating process:
capture program drop panelframe
program define panelframe
	qui sort firm period
	qui gen fe=rnormal() if period==1
	qui replace fe=fe[_n-1] if fe==.
	qui gen x=rnormal()+sqrt(.5)*fe
	qui gen x2=rnormal()
    qui gen u=exp(rnormal())+sqrt(.5)*fe // long-tailed error distibution
	qui gen y=x+u
	qui sort firm period
	drop u fe
end

panelframe

timer clear
timer on 1
robreg m y i.period x x2, ivar(firm) cluster(group) eff(95)
timer off 1

timer on 2
robtwfe m y x x2, ivar(firm) tvar(period) cluster(group) eff(95) 
timer off 2

timer list
