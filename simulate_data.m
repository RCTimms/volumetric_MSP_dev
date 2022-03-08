function [Y,x,fs,duration,fsig]=simulate_data(L,source_idx,options);
try N_trials=options.N_trials; catch, N_trials=50; end
try fs=options.fs; catch, fs=300; end
try duration=options.duration; catch,duration=1; end
try fsig=options.fsig; catch, fsig=20; end
try SNRdB=options.SNRdB; catch, SNRdB=0; end
t=linspace(0,N_trials*duration,duration*fs*N_trials);


% Single source time series
x=sin(2*pi*fsig*t);

% Project
if numel(source_idx)==1; % Single voxel
    Y=L(:,source_idx)*x;
else % Multiple voxels
    Y=L(:,source_idx)*repmat(x,numel(source_idx),1);
end

% Add noise
allchanstd=std(Y');
meanrmssignal=mean(allchanstd);
whitenoise = meanrmssignal.*(10^(-SNRdB/20));
Y=Y+randn(size(Y)).*whitenoise;

% Make it trial wise, plot to sanity check
Y=reshape(Y,[size(Y,1),duration*fs,N_trials]);
mean(Y,3);
end