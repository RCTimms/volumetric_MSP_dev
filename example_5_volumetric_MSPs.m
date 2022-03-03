%% Volumetric MSPs - Example 5 Script
% A volumetric implementation of MSPs
% This example solves a slightly more complex source reconstruction
% problem:
% 1) The source is superficial
% 2) The SNR is high
% 3) There are a large number of trials
% 4) We solve for all 3 components of the lead field
% 5) We go HAM and try to do this on a 5mm grid
% 6) Finally, we collapse the timeseries into one PC for the x,y,z
%    timeseries
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
options.sim=1; % Make a D object, otherwise use a pre-existing one

% Make a dummy D object with the sensor array
if options.sim==1;

try
    D=spm_eeg_load('dev_MSPs.mat');
catch
    S =[];
    S.space = 50; % Keep it sparse to start
    S.sMRI=1;
    S.lead=0;
    S.fname='dev_MSPs';
    [D] = spm_opm_sim(S);
end
else
   D=spm_eeg_load('dev_MSPs.mat') 
end
clear S 

% We need to borrow some FieldTrip functionality
% Get the data object
data=ftraw(D);

% Get the headmodel, calculated via SPM
headmodel=D.inv{1}.forward.vol;

% Get the source grid - MNI aligned
load('D:\fieldtrip-master\fieldtrip-master\template\sourcemodel\standard_sourcemodel3d6mm.mat')
sourcemodel=ft_convert_units(sourcemodel,headmodel.unit);

% Compute the leadfields
cfg             = [];
cfg.vol         = headmodel;
cfg.grid        = sourcemodel;
cfg.grad        = data.grad;
cfg.channel     = data.label;
cfg.normalize   = 'no';
cfg.reducerank  = 2;
gridLF = ft_prepare_leadfield(cfg);

% Extract the leadfield
L = gridLF.leadfield(gridLF.inside);

% Unravel the leadfield 
L=unravel_leadfield(L);

% Now we reconstruct activity in all 3 directions, so don't extract any of
% the lead field components
L=L;

% Choose a source index, any source index. We want the same one as in
% example 1. The indexing/logic gets a bit confusing here as the lead
% field index is not the same as the grid index
source_idx=length(L)-2;
grid_idx=length(L)/3;

% Visualise the source and scalp topo
show_grid_and_ground_truth(gridLF,grid_idx);
spm_eeg_plotScalpData(L(:,source_idx),D.coor2D,D.chanlabels);

% Simulate the data
[Y,X_GT,fs,duration,fsig]=simulate_data(L,source_idx);

% Put data into D object
Dnew = clone(D, 'ryan', [D.nchannels 300 50]);
Dnew(:,:,:)=Y;
Dnew.save;

%% Off to volumetric MSPs
D=Dnew;
D.inv{1}.inverse=[];D.inv{1}.inverse.type='GS';D.inv{1}.inverse.Han=0;
Dout=spm_eeg_invert_classic_volumetric(D,1,L);

% Show the results
is_vector=1; % Is it a vector lead field? YES
log_it=0; % Should we log the power map?
[X,X_x,X_y,X_z]=get_timeseries(Dout,is_vector);

% Show the power and the ground truth location in x, y and z
show_power(X_x,gridLF,grid_idx,log_it);title('Power - X');view([-144 45])
show_power(X_y,gridLF,grid_idx,log_it);title('Power - Y');view([-144 45])
show_power(X_z,gridLF,grid_idx,log_it);title('Power - Z');view([-144 45])

% Show the reconstructed timeseries
figure;
t=linspace(0,duration,duration*fs);
hold all
plot(t,X(source_idx,:),'LineWidth',2);
plot(t,X_GT(1:fs),'Linewidth',2)
xlabel('Time (s)');
set(gca,'FontSize',16);
legend('Reconstructed','Ground Truth');

% Collapse the source timeseries
X_collapsed=zeros(size(X_x));
for i=1:size(X_x,1);
    
    % Concatenate x,y,z timeseries
    X_1=[X_x(i,:);X_y(i,:);X_z(i,:)];
    
    % The timeseries must be observations x channels
    X_1=X_1';
    
    PCs=do_PCA(X_1);
    X_collapsed(i,:)=PCs(:,1)';
end



