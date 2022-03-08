%% Volumetric MSPs - Example 8 Script
% A volumetric implementation of MSPs
% This example entails solving a problem which isn't in MNI space. We warp
% the grid as required (and is the most likely situation when dealing with
% real data in SPM).
%
% Highlights the importance of lead field normalisation for when using MSPs
% with OPMs... I think
%
% Does a **non-linear** transform of the MNI grid into the OPM array. This
% assumes that the MRI for the subject is in a sensible, well-defined
% coordinate system, otherwise FieldTrip will freak out and just use a
% linear transformation. We demonstrate how to deal with this issue for an
% MRI which is in a poorly defined space.
%
% Demonstrates the use of an atlas to reconstruct/extract to a region of
% interest (ROI)
%
% Overview:
% 1) The source is in the left hippocampus
% 2) The SNR is high (0dB)
% 3) There are a large number of trials
% 4) We solve for one component of the lead field
% 5) We use a 10mm grid.
%
% Ryan Timms, 2022. If you use this script and ever meet me, please buy me
% a beer. Alternatively pick up some litter.

clc; clear all; close all;
cd D:\volumetric_MSP_dev
addpath('D:\spm')
spm('defaults', 'eeg')
addpath('D:\OPM')
addpath('D:\optitrack');
addpath D:\play

% Run-time options
options=[];

% This time we load in a D object which is in a native OPM/MRI space
try
D=spm_eeg_load('edffsub-001_ses-001_task-phantom_20Hz_run-003_meg.mat')
catch
fname ='D:\PhantomData\sub-001\ses-001\meg\sub-001_ses-001_task-phantom_20Hz_run-003_meg';
S =[];
S.data = [fname,'.bin'];
S.meg=[fname,'.json'];
S.channels=[fname,'.tsv'];
S.positions='D:\PhantomData\sub-001\ses-001\meg\sub-001_ses-001_task-ERN_run-001_positions.tsv';
S.lead=1;
S.sMRI='D:\PhantomData\mmsMQ0484_orig.img';
D = spm_opm_create(S);
end
data=D.ftraw;

% Optional: if your MRI is in a funky or ill-defined space (like native MRI
% scanner space), you will have to manually re-define the space to something
% more standard like CTF or Neuromag space. This can be achieved by
% uncommenting the following code and clicking some buttons:

% S.sMRI='D:\PhantomData\mmsMQ0484_orig.img'
% 
% addpath(genpath('D:\fieldtrip-master\fieldtrip-master\external\freesurfer'))
% mri=ft_read_mri(S.sMRI);
% mri.coordsys='neuromag'; 
% 
% cfg=[];
% cfg.method='interactive';
% ctf_mri=ft_volumerealign(cfg,mri)
% save ctf_mri ctf_mri

% Load in the MRI which is in CTF coords and show this w.r.t the sensors
load ctf_mri
ft_determine_coordsys(ctf_mri,'interactive','no');
ft_plot_sens(D.sensors('MEG'),'facecolor','r')
view([180 0])

% This is obviously wrong. Critically, we need to transform our sensors to
% the same space
ctf_grad=ft_transform_geometry(ctf_mri.transform/ctf_mri.transformorig,data.grad);
ft_determine_coordsys(ctf_mri,'interactive','no');
ft_plot_sens(ctf_grad,'facecolor','r')
view([180 0])

% Template source space
% Load in a grid which is in MNI space
load('D:\fieldtrip-master\fieldtrip-master\template\sourcemodel\standard_sourcemodel3d10mm.mat')
template_grid_MNI=sourcemodel;
clear sourcemodel

% We now non-linearly warp the MNI grid to the subject's brain in the
% well-defined CTF MEG coordinate system. The output of this is the grid
% that we use to reconstruct our data to.
cfg           = [];
cfg.warpmni   = 'yes';
cfg.method = 'basedonmni';
cfg.template  = template_grid_MNI;
cfg.sourcemodel.nonlinear = 'yes';
cfg.mri       = ctf_mri;
template_grid_SS         = ft_prepare_sourcemodel(cfg); % Template grid in subject space

% Visualise to see if this worked
ft_determine_coordsys(ctf_mri,'interactive','no');
ft_plot_sens(ctf_grad,'facecolor','r')
ft_plot_mesh(template_grid_SS.pos(template_grid_SS.inside,:),'vertexcolor','c')
view([180 0])

% Get the headmodel, calculated via SPM. But don't forget the transform to
% CTF space. Note that none of this extra realignment would be necessary if
% we had got the data in a sensible space to start with.
headmodel=D.inv{1}.forward.vol;
headmodel=ft_convert_units(headmodel,'mm'); % Important - Gran's gotchya
ctf_headmodel=ft_transform_geometry(ctf_mri.transform/ctf_mri.transformorig,headmodel);

% Again, visualise to sanity check everything
figure;
ft_plot_mesh(template_grid_SS.pos(template_grid_SS.inside,:),'vertexcolor','r');
hold all
ft_plot_sens(ctf_grad)
ft_plot_headmodel(ft_convert_units(ctf_headmodel,'mm'))

% Compute the leadfields. Note again that everything is in the same CTF space
cfg             = [];
cfg.vol         = ctf_headmodel;
cfg.grid        = template_grid_SS;
cfg.grad        = ctf_grad;
cfg.channel     = data.label;
cfg.normalize   = 'yes';
cfg.reducerank  = 2;
gridLF = ft_prepare_leadfield(cfg);

% Now load in an atlas. We interpolate this with the MNI grid (in MNI
% space). After this interpolation, each voxel index has a one-to-one
% correspondence between the atlas grid, the MNI grid and the grid in
% subject/CTF MEG space. e.g. if voxel 1 in the atlas corresponded to left
% occipital cortex, then template_grid_SS' voxel 1 would also be left
% occipital cortex, etc.
atlas = ft_read_atlas('D:\fieldtrip-master\fieldtrip-master\template\atlas\aal\ROI_MNI_V4.nii');
cfg = [];
cfg.interpmethod = 'nearest';
cfg.parameter = 'tissue';
sourcemodel2 = ft_sourceinterpolate(cfg, atlas,template_grid_MNI);

% CHECK: Did our atlas operation work?
figure; hold all
for lab = 1:length(sourcemodel2.tissuelabel)
    % Find atlas points
    clear atlas_points
    atlas_points    = find(sourcemodel2.tissue==...
        find(contains(sourcemodel2.tissuelabel,atlas.tissuelabel{lab})));
    scatter3(template_grid_MNI.pos(atlas_points,1),template_grid_MNI.pos(atlas_points,2),template_grid_MNI.pos(atlas_points,3),100,'filled')%,'MarkerFaceColor',colz(lab,:))
    
end
axis image

% Extract the leadfield (optional). We now just get the Y component of the
% lead field for each voxel to simplify our lives. We don't need to do this
% - see earlier example scripts where we do a reconstruction to all 
L = gridLF.leadfield(gridLF.inside);
clear L_out
for i=1:size(L,1);
    L_out(:,i)=L{i}(:,2); % Only extract the Y component
end
L=L_out;

% Find all the points in a particular ROI. We first find the tissue
% membership vector from the previously defined atlas. We have to delete
% the voxels which are outside of the brain. We then find all voxels which
% are in ROI, the left hippocampus
atlas_tissue=sourcemodel2.tissue(:); % Get the atlas tissue vector
atlas_tissue(find(~template_grid_MNI.inside))=[]; % Delete bits outside of brain
indx = find(atlas_tissue==37); % Left hippocampus

% Simulate data
% Use the atlas to define the source position
grid_idx=indx(1); % coordinates in x/y/z
source_idx=grid_idx; % Lead field index

% Visualise the source and scalp topo
show_grid_and_ground_truth(gridLF,grid_idx);

% Simulate the data
options.SNRdB=0;
[Y,X_GT,fs,duration,fsig]=simulate_data(L,source_idx);

% Put data into D object
Dnew = clone(D, 'ryan', [D.nchannels 300 50]);
Y0=zeros(size(Dnew));
Y0(D.indchantype('MEG'),:,:)=Y;
Dnew(:,:,:)=Y0;
Dnew.save;

% Off to volumetric MSPs
D=Dnew;
D.inv{1}.inverse=[];D.inv{1}.inverse.type='GS';D.inv{1}.inverse.Han=0;
Dout=spm_eeg_invert_classic_volumetric(D,1,L);

% Show the results
is_vector=0; % Is it a vector lead field? no
log_it=0; % Should we log the power map?
[X]=get_timeseries(Dout,is_vector);

% Show the power and the ground truth location in x, y and z
show_power(X,gridLF,grid_idx,log_it);title('Power');view([-144 45])

% Show the reconstructed timeseries
figure;
t=linspace(0,duration,duration*fs);
hold all
plot(t,X(source_idx,:),'LineWidth',2);
plot(t,X_GT(1:fs),'Linewidth',2)
xlabel('Time (s)');
set(gca,'FontSize',16);
legend('Reconstructed','Ground Truth');
