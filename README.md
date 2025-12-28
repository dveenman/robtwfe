# robtwfe

by David Veenman (University of Amsterdam)

**Comments/feedback welcome**

`robtwfe` is a program that can be used to combine a Huber M robust regression estimator with two-way fixed effects (FE). For a firm-time panel dataset, it mimics the program `robreg` for Huber M-estimation with one FE dimension (firm) absorbed and the other FE dimension (time) included as indicator variables. Instead of including the time indicator variables, the program leverages (a) the functionality of `reghdfe` and (b) the fact that the iterative reweighting procedure within the robust regression estimation relies on a sequence of weighted least squares regressions, which can be combined with two-way FE. When the second (time) dimension becomes large (>50), the program is substantially faster than `robreg` but provides the same estimates. For a smaller second (time) dimension, `robreg` is typically faster.

---

The program can be installed by simply saving the \*.ado file into your local Stata directory that contains additional ado programs. To identify this folder, type and execute "sysdir" in your Stata command window and go to the folder listed at "PLUS:". Make sure you place the program in the folder with the correct starting letter of the program (i.e., the folder named "r") and the file extension is correctly specified as \*.ado.

The program requires `moremata`, `reghdfe`, `hdfe`, and `robreg` to be installed in Stata:
```
ssc inst moremata, replace
```

---

Syntax:

**robtwfe m depvar indepvars, ivar(string) tvar(string) cluster(string) eff(numeric) [options]**

 - **depvar** is the dependent variable;
 - **indepvars** is/are the independent variable(s);
 - **ivar()** refers to the variable indicating the unit fixed effect to be absorbed, similar to the option in `robreg`;
 - **tvar()** refers to the variable indicating the time fixed effect to be absorbed;
 - **cluster()** refers to the variable specifying the dimension at which standard errors should be clustered;
 - **eff()** refers to the normal efficiency of the robust estimation;

---

Optional program options in **[options]**:

- **tol(numeric)**: specifies the tolerance for the iterative procedures (default 1e-10);
- **weightvar(string)**: specifies name of new variable to be generated with robust regression weights;
- **omitr2**: do not compute and report the pseudo R2 (small speed benefit);

---

The following example shows how the program can be run and how it mimics the `robreg` estimation (the program can be tested using the **test_robtwfe.do** file): 



---


