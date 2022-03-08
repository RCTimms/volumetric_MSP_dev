%% Volumetric MSPs - Example 7 Script
% A volumetric implementation of MSPs
% This example entails solving a problem which isn't in MNI space. We warp
% the grid as required (and is the most likely situation when dealing with
% real data in SPM).
%
% This script makes use of a linear (affine) transform of the grid which is
% MNI space into the OPM helmet. We **DO** normalise the lead fields in
% this example and get 0mm dipole location error (DLE). Feel free to change
% cfg.normalize   = 'yes'; to cfg.normalize   = 'no' and see what happens 
%
% 1) The source is superficial
% 2) The SNR is mind-boggling high (+20dB)
% 3) There are a large number of trials
% 4) We solve for 1 component of the lead field
% 5) We use a 10mm grid.
%
% Ryan Timms, 2022

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
% S.voltype='Magnetic Dipole'; % Important!
S.sMRI='D:\PhantomData\mmsMQ0484_orig.img'; % We don't actually need the MRI, but this helps us with the coreg
D = spm_opm_create(S);
end

% We need to borrow some FieldTrip functionality
% Get the data object
data=ftraw(D);

% Get the headmodel, calculated via SPM
headmodel=D.inv{1}.forward.vol;

% Get the transform which maps us from "some other space" to MNI.
M=D.inv{1}.mesh.Affine;

% The transform we want fromm MNI to some other space is hence the inverse
% of M
load('D:\fieldtrip-master\fieldtrip-master\template\sourcemodel\standard_sourcemodel3d10mm.mat')
sourcemodel=ft_convert_units(sourcemodel,'mm'); % Gotchya!
sourcemodel=ft_transform_geometry(inv(M),sourcemodel);

% Validate that this worked
figure;
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:),'vertexcolor','r');
hold all
ft_plot_sens(D.sensors('MEG'))
ft_plot_headmodel(headmodel)

% Compute the leadfields
cfg             = [];
cfg.vol         = headmodel;
cfg.grid        = sourcemodel;
cfg.grad        = data.grad;
cfg.channel     = data.label;
cfg.normalize   = 'yes';
cfg.reducerank  = 2;
gridLF = ft_prepare_leadfield(cfg);

% Extract the leadfield
% Now we only take the z component of the lead field for simplicity
L = gridLF.leadfield(gridLF.inside);
for i=1:size(L,1);
    L_out(:,i)=L{i}(:,3);
end
L=L_out;

figure;
tmp=L;
tmp=sum(abs(tmp),1);
xyz=gridLF.pos(gridLF.inside,:);
scatter3(xyz(:,1),xyz(:,2),xyz(:,3),100,tmp','filled')

% Choose a source index, any source index.
% Each point has three lead field vectors. E.g. for grid_idx 1 we want lead
% field vectors 1, 2 and 3. For grid_idx 2 we want 4, 5 and 6 etc.
% grid_idx=2400; % coordinates in x/y/z
% source_idx=(grid_idx*3)-2; % Lead field index

grid_idx=1331; % coordinates in x/y/z
source_idx=grid_idx; % Lead field index

% Visualise the source and scalp topo
show_grid_and_ground_truth(gridLF,grid_idx);

% Simulate the data
options.SNRdB=20;
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
is_vector=0; % Is it a vector lead field? YES
log_it=0; % Should we log the power map?
[X]=get_timeseries(Dout,is_vector);

% Show the power and the ground truth location in x
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
