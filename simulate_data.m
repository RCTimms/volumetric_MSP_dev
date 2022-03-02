function [Y,x]=simulate_data(L,source_idx);
N_trials=50;
fs=300;
duration=1;
fsig=20;
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
mean(Y,3);
end