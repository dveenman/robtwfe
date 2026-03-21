# robtwfe

A Stata program to combine robust regression estimation (Huber M) with two-way fixed effects. The program accompanies the Gassen & Veenman (2026) study on ``Estimation Precision and Robust Inference in Archival Research'' (https://ssrn.com/abstract=4975569).

**Comments/feedback welcome**

`robtwfe` is a program that can be used to combine Huber M-estimation with two-way fixed effects (FE). For a firm-time panel dataset, the program mimics `robreg` for Huber M-estimation with one FE dimension (firm) absorbed and the other FE dimension (time) included as indicator variables. Instead of including the time indicator variables, the program leverages (a) the functionality of `reghdfe` and (b) the fact that the iterative reweighting in the robust estimation relies on a sequence of weighted least squares estimations, which can be combined with two-way FE using `reghdfe`. For the first-step quantile regression that is used to obtain the scale estimate, the program implements the MM-QR estimation from [Machado and Santos Silva (2009)](https://www.sciencedirect.com/science/article/pii/S0304407619300648) to combine quantile regression estimation with fixed effects, which [Rios-Avila, Siles, and Canavire-Bacarreza (2024)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4944894) show is easily extended to the multidimensional fixed effects setting. When the second (time) dimension becomes sufficiently large (e.g., >50), the program is substantially faster than `robreg` and provides the same estimates. 

---

Installation:
```
net install robtwfe, replace from(https://raw.githubusercontent.com/dveenman/robtwfe/main/)
```

The program requires `moremata`, `reghdfe`, and `hdfe` to be installed in Stata:
```
ssc inst moremata, replace
ssc inst hdfe, replace
ssc inst reghdfe, replace 
```

---

Syntax:

**robtwfe m depvar indepvars, ivar(str) tvar(str) eff(real) [cluster(varlist) tol(real 0) weightvar(str)]**

 - **depvar** is the dependent variable;
 - **indepvars** is/are the independent variable(s);
 - **ivar()** refers to the variable indicating the unit fixed effect to be absorbed, similar to the option in `robreg`;
 - **tvar()** refers to the variable indicating the time fixed effect to be absorbed;
 - **eff()** refers to the normal efficiency of the robust estimation;
 - **cluster()** refers to the variable specifying the dimension at which standard errors should be clustered (nesting of FE not required); when cluster() is not provided, heteroskedasticity-robust standard errors are computed;
 - **tol(numeric)**: specifies the tolerance for convergence of the iterative reweighted least squares (default 1e-10);
 - **weightvar(string)**: specifies the name of a new variable to be generated with robust regression weights;
 
---

The following examples show the standard output and speed of the program (the program can be tested using the **test_robtwfe.do** file): 

<img width="643" height="665" alt="image" src="https://github.com/user-attachments/assets/c82af52c-429f-4b1c-b8a2-c3d1aef1d9e4" />

<img width="639" height="666" alt="image" src="https://github.com/user-attachments/assets/aa23ea2a-83b9-4c55-8537-90a9b6b51bab" />

---


