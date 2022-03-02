clc; clear all; close all;
cd D:\volumetric_MSP_dev
addpath('D:\spm')
spm('defaults', 'eeg')
addpath('D:\OPM')
addpath('D:\optitrack');
addpath D:\play

%%% Run-time options
options=[];
options.sim=1; % Make a D object, otherwise use a pre-existing one
% These pertain to the simulation settings
options.N_trials=50;
options.fs=300;
options.duration=1;
options.fsig=20;

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
%%
% F off to FieldTrip

% Get the data object
data=ftraw(D);

% Get the headmodel, calculated via SPM
headmodel=D.inv{1}.forward.vol;

% Get the source grid - MNI aligned
load('D:\fieldtrip-master\fieldtrip-master\template\sourcemodel\standard_sourcemodel3d10mm.mat')
sourcemodel=ft_convert_units(sourcemodel,headmodel.unit);

ft_plot_headmodel(headmodel)
hold all
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:))

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
close all
%% Unravel the leadfield 
L=unravel_leadfield(L);
L0=L;
%% Start simple - just take the x component of the lead field
L=L0(:,1:3:end);

% Visualise the source and scalp topo
source_idx=1000;
show_grid_and_ground_truth(gridLF,source_idx);
spm_eeg_plotScalpData(L(:,source_idx),D.coor2D,D.chanlabels);

% Simulate and project the source data: this makes a sine wave and adds
% noise
[Dnew,t]=simulate_source_data(options,L,source_idx,D);


% Off to volumetric MSPs
D=Dnew;
D.inv{1}.inverse=[];
D.inv{1}.inverse.type='GS';
D.inv{1}.inverse.Han=0;
Dout=spm_eeg_invert_classic_volumetric(D,1,L);

% Reconstruct the timeseries
X=Dout.inv{1}.inverse.M*Dout.inv{1}.inverse.Y*Dout.inv{1}.inverse.T';

%
power=log(sum(X.^2,2));
figure;
scatter3(gridLF.pos(gridLF.inside,1),gridLF.pos(gridLF.inside,2),gridLF.pos(gridLF.inside,3),100,power,'filled')
hold all;

gridLF.pos(gridLF.inside,:);

hold all;ft_plot_mesh(ans(source_idx,:),'vertexcolor','r')


figure;
t=linspace(0,duration,duration*fs);
hold all
plot(t,X(source_idx,:));




% Just populate a single vertex
x=sin(2*pi*fsig*t);
plot(t,x)

function [D,t]=simulate_source_data(options,L,source_idx,D);
N_trials=options.N_trials;
fs=options.fs;
duration=options.duration;
fsig=options.fsig;

t=linspace(0,N_trials*duration,duration*fs*N_trials);
SNRdB=0;

% Just populate a single vertex
x=sin(2*pi*fsig*t);

% Project
Y=L(:,source_idx)*x;

% Add noise
allchanstd=std(Y');
meanrmssignal=mean(allchanstd);
whitenoise = meanrmssignal.*(10^(-SNRdB/20));
Y=Y+randn(size(Y)).*whitenoise;

% Make it trial wise, plot to sanity check
Y=reshape(Y,[size(Y,1),duration*fs,N_trials]);


% Put the simulated data into the D object


Dnew = clone(D, 'ryan', [D.nchannels 300 50]);
Dnew(:,:,:)=Y;
Dnew.save;
end