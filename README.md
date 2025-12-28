# robtwfe

A Stata program to combine robust regression estimation (Huber M) with two-way fixed effects 

by David Veenman (University of Amsterdam)

**Comments/feedback welcome**

`robtwfe` is a program that can be used to combine Huber M-estimation with two-way fixed effects (FE). For a firm-time panel dataset, the program mimics `robreg` for Huber M-estimation with one FE dimension (firm) absorbed and the other FE dimension (time) included as indicator variables. Instead of including the time indicator variables, the program leverages (a) the functionality of `reghdfe` and (b) the fact that the iterative reweighting in the robust estimation relies on a sequence of weighted least squares estimations, which can be combined with two-way FE using `reghdfe`. When the second (time) dimension becomes sufficiently large (e.g., >50), the program is substantially faster than `robreg` and provides the same estimates. 

---

The program can be installed by simply saving the \*.ado file into your local Stata directory that contains additional ado programs. To identify this folder, type and execute "sysdir" in your Stata command window and go to the folder listed at "PLUS". Make sure to place the program in the folder with the correct starting letter (i.e., the folder named "r") and the file extension is correctly specified as \*.ado.

The program requires `moremata`, `reghdfe`, `hdfe`, and `robreg` to be installed in Stata:
```
ssc inst moremata, replace
ssc inst reghdfe, replace
ssc inst hdfe, replace
ssc inst robreg, replace
```

---

Syntax:

**robtwfe m depvar indepvars, ivar(string) tvar(string) cluster(string) eff(numeric) [options]**

 - **depvar** is the dependent variable;
 - **indepvars** is/are the independent variable(s);
 - **ivar()** refers to the variable indicating the unit fixed effect to be absorbed, similar to the option in `robreg`;
 - **tvar()** refers to the variable indicating the time fixed effect to be absorbed;
 - **cluster()** refers to the variable specifying the dimension at which standard errors should be clustered (nesting of FE not required);
 - **eff()** refers to the normal efficiency of the robust estimation;

---

Optional program options in **[options]**:

- **tol(numeric)**: specifies the tolerance for convergence of the iterative reweighted least squares (default 1e-10);
- **weightvar(string)**: specifies the name of a new variable to be generated with robust regression weights;
- **omitr2**: do not compute and report the pseudo R2 (small speed benefit);

---

The following example shows the standard output from the program (the program can be tested using the **test_robtwfe.do** file): 

<img width="635" height="465" alt="image" src="https://github.com/user-attachments/assets/788fd38c-420e-43e3-99c1-10e8327e29c3" />

---


