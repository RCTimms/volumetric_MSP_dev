# Development code for a volumetric implementation of Multiple Sparse Priors
Because who likes a mesh?

**Note that this is development only code and has not yet been tried and tested in any anger.**

Based off of the following papers:
>Friston, K., Harrison, L., Daunizeau, J., Kiebel, S., Phillips, C., Trujillo-Barreto, N., ... & Mattout, J. (2008). Multiple sparse priors for the M/EEG inverse problem. NeuroImage, 39(3), 1104-1120.

>Strobbe, G., van Mierlo, P., De Vos, M., Mijović, B., Hallez, H., Van Huffel, S., ... & Vandenberghe, S. (2014). Multiple sparse volumetric priors for distributed EEG source reconstruction. NeuroImage, 100, 715-724.

## Requirements
- SPM12
- MNI aligned "standard" sourcemodels shipped with FieldTrip. These can be downloaded from here: https://github.com/fieldtrip/fieldtrip/tree/master/template/sourcemodel 

(and optionally the OPM repo if you want to change the number of sensors or array setup, etc.).

If you would like a copy of the D object pertaining to the real data (for examples 7 and 8) then please email r.timms@ucl.ac.uk.

## User guide
The function works by providing the volumetric MSP code with both an SPM D object and a source model struct from FieldTrip. The general work flow of the pipeline is as follows:

1. Creation of SPM D object, **coregistration and forward model calculation are essential**.

2. Creation of FieldTrip source model struct. Currently an MNI-aligned volumetric grid is loaded in from disk. 

3. Call to spm_eeg_invert_classic_volumetric. This takes the arguments:
   - D (SPM data object from step 1)
   - the inversion index (an integer, e.g. 1)
   - the unravelled lead field matrix [N_channels x Nsources*3]. This can be obtained by using the unravel_leadfield function.

## Points of note
- The default in this code is to have one prior per voxel and per orientation, not an amalgamation of voxels into "patches". For instance, if there are 1,000 grid points this would lead to 3 x 1000 = 3,000 priors.

- The defualt is to not apply any temporal filtering to the data in this version of the code.


## Example Scripts
There are some scripts provided in the repository designed to demonstrate the function working, alongside highlighting some of the bonus functionality provided. These are called example_N_volumetric_MSPs.m. A brief description of these scripts are provided here. Common to all datasets is the simulation of a single sine wave in the source space which is epoched into trials.

1. A very basic volumetric MSP example. The source is superficial, the SNR is high and we only reconstruct to one lead field orientation. This is an "easy" source localisation problem.

2. Same as above, but now the (single) source is located deep in frontal cortex, so the biological SNR (i.e. lead field coverage) is a bit weaker.

3. Same as example 1, but now we solve for each of the lead field orientations: x, y and z.

4. Same as example 3, but now we iterate over SNRs (-30 to +10dB) and calculate reconstruction statistics (dipole location error and correlation with ground truth signal).

5. Same as example 3, but now we reconstruct to a 6mm grid. This makes the inference more computationally demanding, as there are more priors created in the MSP algorithm. We also demonstrate the collapsing of the reconstructed timeseries into one component from three orthogonal timeseries.

6. Coming soon...

7. A more realistic example. We use a linear (affine) transform to map the MNI grid into the subject's OPM array. A volumetric grid is used. We simulate for one of the lead field axes. The SNR is very very high. **We normalise the lead field vectors - see what happens if we don't do this step**.

8. This is a very comprehensive script (and hence is a little complicated). Starting with the high level details it: does a non-linear coreg of the MNI grid into the OPM helmet; shows how to deal with an MRI file which is in a non-standard space (e.g. MRI scanner space) and put this into a well-defined one (e.g. CTF MEG space); loads in an atlas; simulates a voxel timeseries in the left hippocampus and runs volumetric MSPs. We do this for just one of the lead field orientations (Y - chosen arbitrarily). The SNR is 0dB. **We normalise the lead fields - feel free to try changing this and see what happens...  **

## Example Results

## Development ideas
1. Clump voxels together to form patches to reduce the computational burden/number of patches. This could be based off of an atlas
2. Add support to export images of activity to niftis
3. ~Collapse timeseries over x, y and z into a single dimension (vector to scalar source solution).~
4. ~Check with a real OPM array. Maybe we need to do a coregistration...~
5. ~Add atlas support to reconstruct to ROIs (basically pure FieldTrip stuff).~
6. Add example result figures
7. Add support for the specification of multiple noise covariance matrices
