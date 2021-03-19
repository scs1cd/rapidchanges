# rapidchanges

# Workflow
To obtain rates of change at Earth's surface requires running the following codes: 

ylm2glm ->
* brdot -> christensen_criteria (optional)
* paleo -> spikeI -> python notebooks

These notes discuss these codes in turn. 

# Ylm2glm

This codes takes the poloidal spherical harmonic coefficients from the statefiles and produces Gauss coefficients at the surface. 
This code can be found in: 

/nfs/a88/earcd/CODES/PALEO_LIBS

To compile do: make

To run the code, go to a directory with statefiles and run as

ylm2glm 10 100.0 839.0 925.0 > OUT_ylm2glm &

Here the arguments are: number of file to start from, number of magnetic diffusion times to process, Rm for model, Rm for Earth. 

The values of Rm are used to convert the timescale in the simulations, which uses magnetic diffusion time, into dimensional 
units. If you do not mind about this then you can put any numbers here. 

Outputs: 
* gauss_coeffs_surface - columns: l, m, glm. 
* gauss_coeffs_time - columns: dimensionless time, dimensional time (diffusion times). 

# Paleo

This code takes gauss_coeffs_surface and produces the magnetic elements X, Y,Z, I, D, F at Earth's surface on a regular grid. 
This code can be found in: 

/nfs/a88/earcd/CODES/PALEO_LIBS

To compile do: make

To run, go to the directory containing gauss_coeffs_surface and do

paleo < paleo_input.in

The input file paleo_input.in contains comments that describe the inputs. 

There are a number of output files, which will be produced in a new directory:
* surface_ts_dir_00001.dat - time series of directional data: : time, lat, lon, I_dipole, I, Ianomaly=I-Idipole, D
* surface_ts_vec_00001.dat - time series of vector data: time, lat, lon, X, Y, Z, F, F_dipole
* surface_ts_vgp_00001.dat - time series of VGP data (x, and y are the Mollweide projection of the VGP coordinates)
* A number of _ave_ files containing averages. 

# SpikeI

This code takes a file of the form surface_ts_vec_00001.dat and produces rates of change on a grid of points.  
The code can be found in:

/nfs/a88/earcd/GAUSS_SPIKE/DYNAMO_CODES/

To compile do

ifort -o spikeI spikeI.f90

This code produces the following files. All filenames use the same stem, which is based on the input file name. 
If the input file is called LSMOD.2_tllxyz3 then we get
* LSMOD.2_tllxyz3_spike_the_alldt_lr=90_type=3_ns=0001 - rate of change at every point on the grid and every point in time. Columns are lat, lon, dtdt, dvgpdt, dvdmdt, dfdt, time of fastest change at this location
* LSMOD.2_tllxyz3_spike_int_maxts_lr=90_type=3_ns=0001 - data at the location corresponding to location of maximum change in intensity (int), i.e. dfdt. 
* LSMOD.2_tllxyz3_spike_the_maxts_lr=90_type=3_ns=0001 - as above but for theta (\hat{B})
* LSMOD.2_tllxyz3_spike_vdm_maxts_lr=90_type=3_ns=0001 - as above but for VDM (P_V)
* LSMOD.2_tllxyz3_spike_vgp_maxts_lr=90_type=3_ns=0001 - as above but for VGP (\hat{P_V})
* LSMOD.2_tllxyz3_spike_the_lr=90_type=3_ns=0001 - ax change in theta, VGP, VDM and F at each point

Using the notation in Davies & Constable (2020) we have the following mapping: 
* dtdt = d \hat{B} / dt
* dvgpdt = d \hat{P_V} / dt
* dfdf = d B / dt
* dvdmdt = d P_V / dt

# Windowing

In Mound,Davies & Silva 2015 we took the output from ylm2glm and windowed the coefficients in 400-yr and 3000yr windows in order to compare with gufm1 and CALS3k. 

NOTEL: gauss_coeffs_surface contains no header. However output of the windowing, gauss_coeffs_surface_winave, does contain 1 header line. brdot (below) assumed 1header line and writes out a warning. 

# BRDOT

This code produces the inputs to what was lebintdev.m - now Calculate-Hsv. It takes an ascii input file that looks like this: 

** lmax, # timepoints, radius to evaluate SV

32 5470 3485.0

** Name of input file in format l,m, glm/hlm

gauss_coeffs_surface

** axisym (1=yes); EA only (yes=1) g10 only (yes=1); trunc below lmax?

0 0 0 14

and produces an output file that is the input file name appended with "lebedev_lmax=[trunc]" where trunc is the truncation of the original SH file. The output files contains columns: 

phi, theta, x, y, z, abs(z). 

Note that you can use sv.f90 to generate SV coefficients in the same format as Gauss coefficients and feed these into brdot. The result will then be phi, theta and the time derivatives of x, y, z, and abs(z). 

# Calculate-Hsv
Python implementation of Lebedev Integration to determine the "hemisphericity" (or "pacificity") of the SV

Conversion from matlab code, based on pre-calculated quadrature points and weights.

lebintdev.m is the matlab code that is up for conversion

lebedev_XXX.txt are precalculated quadrature points and weightings (rows are east longitude, colatitude, weight)
these have been obtained from https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html

Python Notebook is for work in progress, once a flow is set this should be converted to a standard python script

output is columns of timepoints, Hsv values -- (no header)

# Compliance-SV
Python script to determine the Mound et al (2015) compliance of Hemispheric SV ("H_rho") relative to gufm values, compliance(Hsv) = abs((H_rho_model - gufm_avg)/gufm_stddev)

Reads in the output from Calculate-Hsv

output is columns of timepoints, compliance(Hsv) -- (no header)

# Hsv-summary-stats

Once the Hsv and its compliance have been calculated, this python script will generate some summary stats. 

Input filename needed is the one from Calculate-Hsv

Output file has a headerline and gives: total number of timpoints, total time elapsed, avg value, standard deviation of value
for both Hsv and its compliance (one file for each)
