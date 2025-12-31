*! version 1.0.1 20251230 David Veenman

/*
20251230: 1.0.1     Made areg default given faster execution, reghdfe is used for Stata versions below 19
20251224: 1.0.0     First version

Dependencies:
   moremata
   reghdfe
   hdfe
   robreg
*/

program define robtwfe, eclass sortpreserve
	version 15
	syntax [anything] [in] [if], cluster(varlist) ivar(str) tvar(str) eff(real) [tol(real 0) weightvar(str) omitr2]

	capt findfile mf_mm_aqreg.hlp
	if _rc {
		di as error "Program requires the {bf:moremata} package: type {stata ssc install moremata, replace}"
		error 499
	}

	capt findfile reghdfe.ado 
	if (_rc & _caller()<19) {
		di as error "Program requires the {bf:reghdfe} package: type {stata ssc install reghdfe, replace}"
		error 499
	}
	
	local stataversion=_caller()
	
	capt findfile hdfe.ado 
	if _rc {
		di as error "Program requires the {bf:hdfe} package: type {stata ssc install hdfe, replace}"
		error 499
	}

	capt findfile robreg.ado 
	if (_rc & "`omitr2'"=="") {
		di as error "Program requires the {bf:robreg} package: type {stata ssc install robreg, replace}"
		error 499
	}

	marksample touse
		
	tokenize `anything'
	local subcmd `"`1'"'

	local cmdlist "m mm"
	if !`:list subcmd in cmdlist' {
		di as err `"Invalid subcommand: `subcmd'"'
		exit 
	}
		
	macro shift 1
	local depv `"`1'"'
	local varlist `"`*'"'

	// Ensure dv is not a factor variable:
	_fv_check_depvar `depv'
	macro shift 1
	local indepv "`*'"
	
	// Ensure iv list does not contain a factor variable:
    fvexpand `indepv'
    if "`r(fvops)'" == "true" {
        di as err "Independent variable list may not contain factor variables"
        exit 
    }
	else {
		local indepv `r(varlist)'
	}
	
	// Mark out missing observations:
	markout `touse' `depv' `indepv'

	// Check for collinearity:
	_rmcoll `indepv'
	local k_omitted=r(k_omitted)
	
	// Check number of independent variables:
	local varcount=0
	foreach v of local indepv {
		local `varcount++'
	}
	scalar k0 = `varcount'
	
	// Check absorb variables:
	local n1: word count `ivar'
	if (`n1'!=1){
	    di as err "ERROR: Option ivar() may contain only one variable"
		exit		
	}
	local n2: word count `tvar'
	if (`n2'!=1){
	    di as err "ERROR: Option tvar() may contain only one variable"
		exit		
	}
	
	// Check nesting of FE in clusters:
	capture assertnested `cluster' `ivar'
	local nest1=(_rc!=0)
	capture assertnested `cluster' `tvar'
	local nest2=(_rc!=0)
	if (`nest2'==1) {
		local nest1dof = 0
	}
	else {
		local nest1dof = `nest1'
	}
	local nest2dof = `nest2'
	
	// For LAD estimation, expand indepv list to include time indicators:
	local indepv2 "i.`tvar' `indepv'"
	fvexpand `indepv2'
	local indepv2 `r(varlist)'
	local fvcheck `r(fvops)'
	
	// Set tolerance:
	if (`tol'!=0){
		local tolerance=`tol'
	}
	else {
		local tolerance=1e-10
	}	

	// Check efficiency:
	
	if (`eff'<63.7 | `eff'>99.9) {
		di as err "ERROR: Normal efficiency must be between 63.7 and 99.9"
		exit
	}
	
	local nc: word count `cluster'
	if (`nc'>=2){
	    di as err "ERROR: Maximum number of dimensions to cluster on is one"
		exit
	}
	local clusterdim1: word 1 of `cluster'
	local clusterdim2: word 2 of `cluster'
	
	// Create temporary variables: 
	tempvar clus1  
	qui egen double `clus1'=group(`clusterdim1') if `touse'
		
    /////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////
	di ""
	di as text "STEP 1: Estimating LAD and obtaining scale estimate"

	qui sum `depv' if `touse'
    local N=r(N)
	tempvar ivarid tvarid _resid_temp w phi
	qui gen double `phi'=.
	qui egen double `ivarid'=group(`ivar') if `touse'
	qui sum `ivarid'
	local ni=r(max)
	qui egen double `tvarid'=group(`tvar') if `touse'
	qui sum `tvarid'
	local nt=r(max)
	local Kinit: word count `indepv2' 
	local Kinit = `Kinit' - `k_omitted'
    scalar df_initial=`N'-`ni'-(`Kinit'-1)

	local ni_red = (1-`nest1')*`ni' + `nest1dof'
	local nt_red = (1-`nest2')*`nt' + `nest2dof'
	local ni_est = `ni' - (1-`nest1')*`ni' - `nest1dof'
	local nt_est = `nt' - (1-`nest2')*`nt' - `nest2dof'
	local K: word count `indepv' 
	local K = `K' + 1 + `ni_est' + `nt_est'
	
	scalar eff=`eff'
	mata: _lad_initial()

	di as text "STEP 2: Iterating IRWLS"

    local diff=100
	local maxiter=c(maxiter)
    forvalues i=1(1)`maxiter'{
        if `diff'>`tolerance' {
			qui capture drop `_resid_temp'
			if (`stataversion'<19) {
				qui capture reghdfe `depv' `indepv' [aw = `w'] if `touse', absorb(`ivar' `tvar') dof(none) notable nofootnote noheader resid keepsin
				qui ren _reghdfe_resid `_resid_temp' 
			}
			else {
				qui capture areg `depv' `indepv' [aw = `w'] if `touse', absorb(`ivar' `tvar') noabs 
				qui predict `_resid_temp', res 
			}
			matrix b=e(b)            
            if `i'>1 {
                local diff=mreldif(b0,b)
            }
            matrix b0=b
			if (`i'==`maxiter' & `diff'>`tolerance'){
				di as err "ERROR: Convergence not achieved"
				exit
			}
			if (`diff'>`tolerance'){
				qui gen double _z_temp=`_resid_temp'/scale
				mata: _update_weights()
				drop _z_temp
			}
        }
    }

	if ("`weightvar'"!="") {
		capture drop `weightvar'
		if _rc==0 {
			local replaceweightvar "yes"
		}
		gen double `weightvar' = `w' 
	}
	
	qui replace `phi'=1e-20 if `phi'==0 // Ensure that residualized values are also created for phi=0 cases
	qui hdfe `indepv' if `touse' [aw = `phi'], absorb(`ivar' `tvar') gen(_tilde_) keepsin
	local indepvr ""
	foreach v of local indepv {
		local indepvr "`indepvr' _tilde_`v'"
	}

	// For calculation of Pseudo R2:
	if "`omitr2'"=="" {
		mata: _madn()
		matrix mu0=mu0
		local madn=scale0
		qui robreg m `depv' if `touse', eff(95) nor2 nose tol(`tolerance') init(mu0) scale(`madn')
		scalar mu=e(b)[1,1]				
	}
	
	mata: ""
	di as text "STEP 3: Computing standard errors"

	sort `clus1' 
	local cvar "`clus1'"	
	mata: _vce_cluster()    
	local nclusterdim1=mata_nclusters
	local e_df_r=mata_nclusters-1
	matrix beta=b0[.,1..k0]
	matrix Vc=Vclust
	    	
	local factor=(`nclusterdim1'/(`nclusterdim1'-1))*((`N'-1)/(`N'-`K'))
	matrix Vc=`factor'*Vc
    
	ereturn clear
	tempname b V

	matrix colnames Vc=`indepv'
	matrix rownames Vc=`indepv'
    matrix colnames beta=`indepv'	
	
	matrix `b'=beta
	matrix `V'=Vc
	
	ereturn post `b' `V'
	ereturn scalar df_r=`e_df_r'
	ereturn scalar N=`N'
	if "`omitr2'"=="" {
		ereturn scalar r2_p=r2_p
	}
	ereturn scalar N_clust=`nclusterdim1'
	
	di ""
	di in green "Huber M-estimation with `eff'% normal efficiency and twoway fixed effects by " ///
		in yellow "`ivar'" in green " and " in yellow "`tvar'"
	di ""
	di _column(51) in green "Number of obs = " %12.0fc in yellow e(N)
	if "`omitr2'"=="" {
		di _column(51) in green "Pseudo R2" _column(65) "= " %12.4f in yellow e(r2_p)
	}
	
    ereturn display
    
	di "SE clustered by " `nclusterdim1' " clusters in " in yellow "`clusterdim1'" 

	if "`weightvar'"!="" {
		di in green "Robust weights stored in " in yellow "`weightvar'" 	
	}
	
	if "`replaceweightvar'"!="" {
		di in green "Careful: " in yellow "`weightvar'" in green " already existed and now replaced with new data"
	}

	local offset1 = 28 - strlen("`ni'")
	local offset2 = 28 - strlen("`nt'")
	local offset3 = 40 - strlen("`ni_red'")
	local offset4 = 40 - strlen("`nt_red'")
	local offset5 = 52 - strlen("`ni_est'")
	local offset6 = 52 - strlen("`nt_est'")
	if (`nest1'==0) {
		local star1 = "*"
	}
	else {
		local star1 = " "		
	}
	if (`nest2'==0) {
		local star2 = "*"
	}
	else {
		local star2 = " "		
	}
	
	di ""
	di in green "Degrees of freedom used by FE:"
	di "{hline 17}{c TT}{hline 36}{c TRC}"
	di "FE dimension: {col 18}{c |}  Categories - Redundant: {col 55}{c |}"
	di "{hline 17}{c +}{hline 36}{c RT}"
	di in green "`ivar' {col 17} {c |}" _column(`offset1') " `ni'" "   - " _column(`offset3') (1-`nest1')*`ni' + `nest1dof' "   = " _column(`offset5') in yellow `ni' - (1-`nest1')*`ni' - `nest1dof' " `star1' {col 53}{c |}"
	di in green "`tvar' {col 17} {c |}" _column(`offset2') " `nt'" "   - " _column(`offset4') (1-`nest2')*`nt' + `nest2dof' "   = " _column(`offset6') in yellow `nt' - (1-`nest2')*`nt' - `nest2dof' " `star2' {col 53}{c |}"
	di "{hline 17}{c BT}{hline 36}{c BRC}"
	if (`nest1'==0 | `nest2'==0) {
		di in green "* FE nested within cluster; treated as redundant for DoF calculation"
	}
	
	matrix drop beta Vc Vclust b b0
	scalar drop df_initial eff mata_nclusters scale krob
	foreach v of local indepv {
		drop _tilde_`v'
	}
	
end

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
// Mata programs
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

mata:
	mata clear
	void _vce_cluster() {
		st_view(y=., ., st_local("depv"), st_local("touse"))
		st_view(Xr=., ., tokens(st_local("indepvr")), st_local("touse"))
		st_view(r=., ., tokens(st_local("_resid_temp")), st_local("touse"))
		st_view(cvar=., ., tokens(st_local("cvar")), st_local("touse"))
		scale=st_numscalar("scale")		
		krob=st_numscalar("krob")
		omitr2=st_local("omitr2")
		if (omitr2=="") {
			mu=st_numscalar("mu")
		}
		
		// Process input:
		k=cols(Xr)
		z=r:/scale
		psi=mm_huber_psi(z,krob)
		phi=mm_huber_phi(z,krob)			
		
		// Compute VCE:
		XphiXinv=invsym(quadcross(Xr,phi,Xr))
		info=panelsetup(cvar, 1)
        nc=rows(info)
        M=J(k,k,0)
		for(i=1; i<=nc; i++) {
			xi=panelsubmatrix(Xr,i,info)
			psii=panelsubmatrix(psi,i,info)
			M=M+(xi'*psii)*(psii'*xi) 
		}
		
		// Combine:
		Vclust=makesymmetric(scale^2*XphiXinv*M*XphiXinv)
		
		// Export to Stata:
		st_matrix("Vclust",Vclust)
		st_numscalar("mata_nclusters",nc)

		// Compute pseudo-R2:
		if (omitr2=="") {
			z0=(y:-mu):/scale
			rho=mm_huber_rho(z,krob)			
			rho0=mm_huber_rho(z0, krob)
			r2_p=1-(colsum(rho)/colsum(rho0))
			st_numscalar("r2_p", r2_p)
		}
	}
	    
	void _lad_initial() {
		st_view(y=., ., tokens(st_local("depv")), st_local("touse"))
		st_view(X=., ., tokens(st_local("indepv2")), st_local("touse"))
        st_view(ivar=., ., tokens(st_local("ivar")), st_local("touse"))
        df=st_numscalar("df_initial")
		eff=st_numscalar("eff")
        n=rows(X)		
		S = mm_aqreg(y, ivar, X)
        e = mm_aqreg_e(S)
        p = (2*n - df) / (2*n) 
        scale=mm_quantile(abs(e), 1, p) / invnormal(0.75) // For consistency with robreg
        z = e / scale
		krob=mm_huber_k(eff)
		w=mm_huber_w(z, krob)
		st_store(., st_addvar("double", st_local("w")), st_local("touse"), w)
        st_numscalar("scale", scale)
		st_numscalar("krob", krob)
	}
	
	void _update_weights() {
        z = st_data(., "_z_temp")
		eff=st_numscalar("eff")
		krob=mm_huber_k(eff)
		phi=mm_huber_phi(z,krob)			
		w=mm_huber_w(z, krob)
        st_store(., st_local("w"), w)
        st_store(., st_local("phi"), phi)
        printf(".")
    }

	void _madn() {
		st_view(y=., ., tokens(st_local("depv")), st_local("touse"))
		mu0 = mm_median(y)
        scale0 = mm_median(abs(y :- mu0)) / invnormal(0.75)
        st_numscalar("mu0", mu0)
		st_numscalar("scale0", scale0)
	}
    
end
	
	
	