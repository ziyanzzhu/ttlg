# ttlg - electronic structure calculation of twited trilayer graphene

Calculation of the electronic band structure and density of twisted trilayer graphene with two independent twist angles using a momentum space continuum model. All scripts are in MATLAB.

## Citation

For reference of the `ttlg` model, please see and cite the following manuscript: 

"Twisted Trilayer Graphene: a precisely tunable platform for correlated electrons" 

Ziyan Zhu, Stephen Carr, Daniel Massatt, Mitchell Luskin, and Efthimios Kaxiras

Phys. Rev. Lett. 125, 116404



## Contact

Ziyan (Zoe) Zhu: zzhu1 [at] g.harvard.edu

Please contact me with any issues and/or request. 



## Code Descriptions

1. `getRecip.m`: calculate the reciprocal space lattice

2. `Layer.m`: create an object that contains the geometry of 3 monolayers

3. `kDOF_tri.m`: create the k degrees of freedom for a given cutoff

4. `gen_interlayer_terms.m`: construct interlayer Hamiltonian (Koshino et al. 2017 style with w_aa \neq w_ab)

5. `gen_intralayer_terms_dirac.m`: constract intralayer rotated Dirac Hamiltonian 

6. `dos_gauss_smear.m`: calculates the DOS using Gaussian smaering

7. `dos_calc_tri.m`: calculate the ttlg DOS and (optional) save data to folder `/data/`

8. `/geom/` folder contain scripts to calculate the moir\'e of moir\'e lengths and to make Figure 1 of the paper.  


Description of input arguments can be found at the beginning of the individual file. 


## Examples

The following two examples are included. Examples were tested with MATLAB version `MATLAB_R2019b` and `MATLAB_R2020a`. 
To get a more accurate result, increase the value of `k_cutoff` to increase the cutoff radius.  

1. `triG_bands_calc.m`: calculates the band structure at θ12 = 1.3 deg., θ23 = 3.2 deg, output will be saved to folder `/data`

2. `call_dos.m`: calculates the DOS by calling `dos_calc_tri.m` for θ12 = 1.3 deg., θ23 = 2.3, 3.2 deg. 

Outputs will be saved to folder `/data`; default number of parallel workers 4. Parallelize the k-space sampling. Have the option to run on a cluster. 

The DOS is obtained by integrating over the bilayer moir\'e Brillouin zone of L1 and L2 only.  Need to also integrate over the L2 and L3 DOS and overage over the two moir\'e Brillouin zones. `k_cutoff` is set to be 3, resulting in ~1,800 degrees of freedom, and the grid size is 42 x 42. For a more accurate result, need to increase the cutoff radius and adjust the grid sampling. The variable `param` is the Gaussian full-width-half-maximum in eV and needs to be adjusted accordingly. In this example, we set `param = 8e-3`.


## Data Availability

Full density of states dataset used for the paper is available upon request.
