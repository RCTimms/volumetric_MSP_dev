%% Volumetric MSPs - Example 4 Script
% A volumetric implementation of MSPs
% This example iterates over SNRs and assesses performance
% 1) The source is superficial
% 2) The SNR varies
% 3) There are a large number of trials
% 4) The orientation of the source is fixed in the x direction; we solve
%    for all 3 components of the lead field
% 5) We calculate the dipole location error and correlation at the ground
%    truth/maximum dipole
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
load('D:\fieldtrip-master\fieldtrip-master\template\sourcemodel\standard_sourcemodel3d10mm.mat')
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

% Get the original D object, so we can come back to this at the start of
% each for loop
D0=D;

% Now loop over SNRs
counter=1;
for snr=[-30,-20,-10,0,10];
    for iterations=1:5;
        
        D=D0;
        % Simulate the data
        options.SNRdB=snr;
        [Y,X_GT,fs,duration,fsig]=simulate_data(L,source_idx,options);
        
        % Put data into D object
        Dnew = clone(D, 'ryan', [D.nchannels 300 50]);
        Dnew(:,:,:)=Y;
        Dnew.save;
        
        % Off to volumetric MSPs
        D=Dnew;
        D.inv{1}.inverse=[];D.inv{1}.inverse.type='GS';D.inv{1}.inverse.Han=0;
        Dout=spm_eeg_invert_classic_volumetric(D,1,L);
        
        % Get the timeseries
        is_vector=1; % Is it a vector lead field? YES
        [X,X_x,X_y,X_z]=get_timeseries(Dout,is_vector);
        
        % Get the dipole location error
        DLE_x(counter,iterations)=get_DLE(X_x,gridLF,grid_idx);
        DLE_y(counter,iterations)=get_DLE(X_y,gridLF,grid_idx);
        DLE_z(counter,iterations)=get_DLE(X_z,gridLF,grid_idx);
        
        % Also get the peak power index so we can do the temporal correlation
        [~,idx_x]=get_DLE(X_x,gridLF,grid_idx);
        [~,idx_y]=get_DLE(X_y,gridLF,grid_idx);
        [~,idx_z]=get_DLE(X_z,gridLF,grid_idx);
        
        corr_x(counter,iterations)=corr(X_x(idx_x,:)',X_GT(1:fs)');
        corr_y(counter,iterations)=corr(X_y(idx_y,:)',X_GT(1:fs)');
        corr_z(counter,iterations)=corr(X_z(idx_z,:)',X_GT(1:fs)');
        
        
        % Show the reconstructed timeseries
        figure;
        t=linspace(0,duration,duration*fs);
        hold all
        plot(t,X(source_idx,:),'LineWidth',2);
        plot(t,X_GT(1:fs),'Linewidth',2)
        xlabel('Time (s)');
        set(gca,'FontSize',16);
        legend('Reconstructed','Ground Truth');
    end
    counter=counter+1;
end

%%% DLE
figure('Position',[0 0 2e3 2e3]);hold all;
plot([-30,-20,-10,0,10],(mean(DLE_x,2)+(std(DLE_x'))')*1e2,'b--')
l1=plot([-30,-20,-10,0,10],mean(DLE_x,2)*1e2,'Color','b','DisplayName','Orientation x (GT)');
plot([-30,-20,-10,0,10],(mean(DLE_x,2)-(std(DLE_x'))')*1e2,'b--')

plot([-30,-20,-10,0,10],(mean(DLE_y,2)+(std(DLE_y'))')*1e2,'r--')
l2=plot([-30,-20,-10,0,10],mean(DLE_y,2)*1e2,'Color','r','DisplayName','Orientation y');
plot([-30,-20,-10,0,10],(mean(DLE_y,2)-(std(DLE_y'))')*1e2,'r--')

plot([-30,-20,-10,0,10],(mean(DLE_z,2)+(std(DLE_z'))')*1e2,'g--')
l3=plot([-30,-20,-10,0,10],mean(DLE_z,2)*1e2,'Color','g','DisplayName','Orientation z');
plot([-30,-20,-10,0,10],(mean(DLE_z,2)-(std(DLE_z'))')*1e2,'g--')

ylim([0 max(ylim)])
legend([l1,l2,l3])
xlabel('SNR (dB)');
ylabel('Dipole location error (mm)')
title('Dipole location over SNRs');
set(gca,'FontSize',16);
%%% /DLE

%%% Correlation at peak voxel
corr_x=abs(corr_x);corr_y=abs(corr_y);corr_z=abs(corr_z);
figure('Position',[0 0 2e3 2e3]);hold all;
plot([-30,-20,-10,0,10],(mean(corr_x,2)+(std(corr_x'))')*1e2,'b--')
l1=plot([-30,-20,-10,0,10],mean(corr_x,2)*1e2,'Color','b','DisplayName','Orientation x (GT)')
plot([-30,-20,-10,0,10],(mean(corr_x,2)-(std(corr_x'))')*1e2,'b--')

plot([-30,-20,-10,0,10],(mean(corr_y,2)+(std(corr_y'))')*1e2,'r--')
l2=plot([-30,-20,-10,0,10],mean(corr_y,2)*1e2,'Color','r','DisplayName','Orientation y')
plot([-30,-20,-10,0,10],(mean(corr_y,2)-(std(corr_y'))')*1e2,'r--')

plot([-30,-20,-10,0,10],(mean(corr_z,2)+(std(corr_z'))')*1e2,'g--')
l3=plot([-30,-20,-10,0,10],mean(corr_z,2)*1e2,'Color','g','DisplayName','Orientation z')
plot([-30,-20,-10,0,10],(mean(corr_z,2)-(std(corr_z'))')*1e2,'g--')

ylim([0 100])
legend([l1,l2,l3])
xlabel('SNR (dB)');
ylabel('Correlation (mm)')
title('Correlation at peak voxel w/ GT');
set(gca,'FontSize',16);
%%% /Correlation at peak voxel