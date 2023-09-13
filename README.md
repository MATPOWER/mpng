![GitHub Logo](/MPNG_User's_Manual/Figures/MPNG_logo_md.PNG)

# An open-source simulation package for solving Optimal Power and Natural Gas Flow problems

-----------------------------------------------------------------------------------

MPNG is a [MATPOWER][1]-based package for solving optimal power and natural gas flow 
problems. MPNG uses the general user nonlinear constraints capability of MATPOWER to 
model the gas network taking into account gas-fired power generators, storage units, 
wells, power-and-gas-driven compressors, and nodes with stratified demand 
(different market segments get different priorities).


System Requirements
-------------------
* [MATLAB][2] version 7.3 (R2016b) or later.
* [MATPOWER][3] version 7.0 or later.


Installation
------------
Installation and use of MPNG requires familiarity with basic operations of MATLAB and
MATPOWER. It is assumed that the user has installed a functional and operating MATPOWER 
version. In short, installing MPNG is as simple as adding all the distribution files 
to the MATLAB's path. The user could either proceed manually with such an addition, 
or run the quick installer released with the package by opening MATLAB at the 
`<MPNG>` directory and typing:

		install_mpng


Citing MPNG
-----------
We do request that publications derived from the use of MPNG explicitly acknowledge 
that fact by including the following cites:

>   S. García-Marín, W. González-Vanegas and C. E. Murillo-Sánchez, "MPNG: A MATPOWER-Based Tool for Optimal Power and
    Natural Gas Flow Analyses," in IEEE Transactions on Power Systems, 2022, doi: [10.1109/TPWRS.2022.3195684][5]

>   S.García-Marín, W.González-Vanegas, and C.E. Murillo-Sánchez, "MPNG: MATPOWER-Natural Gas,"
    2019. [Online]. Available: [https://github.com/MATPOWER/mpng][4]



---- 
 [1]: https://matpower.org
 [2]: https://www.mathworks.com/
 [3]: https://github.com/MATPOWER/matpower
 [4]: https://github.com/MATPOWER/mpng
 [5]: https://doi.org/10.1109/TPWRS.2022.3195684
