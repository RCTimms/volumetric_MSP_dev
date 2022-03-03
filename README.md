# Development code for a volumetric implementation of Multiple Sparse Priors
Because who likes a mesh?

**Note that this is development only code and has not yet been tried and tested in any anger.**

## Requirements
- SPM12
- MNI aligned "standard" sourcemodels shipped with FieldTrip. These can be downloaded from here: https://github.com/fieldtrip/fieldtrip/tree/master/template/sourcemodel 

(and optionally the OPM repo if you want to change the number of sensors or array setup, etc.).

## User guide
The function works by providing the volumetric MSP code with both an SPM D object and a source model struct from FieldTrip. The general work flow of the pipeline is as follows:

1. Creation of SPM D object, **coregistration and forward model calculation (essential)**.
2. Creation of FieldTrip source model struct. Currently an MNI-aligned volumetric grid is loaded in from disk. An additional coregistration step may need to be thrown in here in the future.
3. Call to spm_eeg_invert_classic_volumetric. This takes the arguments:
   - D (SPM data object from step 1)
   - the inversion index (an integer, e.g. 1)
   - the unravelled lead field matrix [N_channels x Nsources*3]. This can be obtained by using the unravel_leadfield function.

## Example scripts

## Points of note

## Development ideas
1. Clump voxels together to form patches to reduce the computational burden/number of patches
2. Add support to export images of activity to niftis
3. Collapse timeseries over x, y and z into a single dimension (vector to scalar source solution).
4. Delete voxels from outside of the volume conductor model
5. Check with a real OPM array. Maybe we need to do a coregistration...
