function [D] = spm_eeg_invert_classic_volumetric(D,val)
%% Volumetric multiple sparse priors
% This version only handles single subject single modality data the removal
% of many scaling factors makes it easier to compare between forward models
%
% Note also that this funtion performs a bit differently from the
% spm_eeg_invert and spm_eeg_invert_classic functions. Namely:
%
% 1) No temporal filtering is carried out to the data by default
% 2) By default, each lead field element has one associated prior (i.e. no
%    "patches" or graph Laplacians are calculated).
% 3) Loreta-like priors/inversions are not (currently) supported.
%
% Ryan Timms and Gareth Barnes, 2023.
%
% Requires:
%           An SPM object, D
%           An inversion value, val
%
% The usual SPM invert shenanigans applies:
% D{i}.inv{val}.inverse:
%     inverse.modality - modality to use in case of multimodal datasets
%     inverse.trials - D.events.types to invert
%     inverse.type   - 'GS' Greedy search on MSPs
%                      'ARD' ARD search on MSPs
%                      'MSP' GS and ARD multiple sparse priors
%                      'IID' minimum norm
%                       'EBB' for empirical bayes beamformer
%     inverse.woi    - time window of interest ([start stop] in ms)
%     inverse.lpf    - band-pass filter - low frequency cut-off (Hz)
%     inverse.hpf    - band-pass filter - high frequency cut-off (Hz)
%     inverse.Han    - switch for Hanning window
%     inverse.Nm     - maximum number of channel modes
%     inverse.Nmax     - maximum number of temporal modes
%     inverse.Nt     - fixed/requested number of temporal modes
%     inverse.Np     - number of sparse priors per hemisphere
%     inverse.sdv    - standard deviations of Gaussian temporal correlation
%     inverse.Qe     - any sensor error components (e.g. empty-room data)
%     inverse.Qe0     - minimum amount of sensor noise power relative to
%                        signal eg 0.1 would correspond to power SNR of 10.0
%     inverse.A       - predefined spatial modes (Nchans*Nmodes) to project
%                       sensor data through
%
% Evaluates:
%
%     inverse.M      - MAP projector (reduced)
%     inverse.J{i}   - Conditional expectation (i conditions) J = M*U*Y
%     inverse.L      - Lead field (reduced UL := U*L)
%     inverse.qC     - spatial covariance
%     inverse.qV     - temporal correlations
%     inverse.T      - temporal projector
%     inverse.U(j)   - spatial projector (j modalities) - derived from data
%     inverse.A      - pre-specified spatial projector
%     inverse.Y{i}   - reduced data (i conditions) UY = UL*J + UE
%     inverse.Is     - Indices of active dipoles
%     inverse.It     - Indices of time bins
%     inverse.Ic{j}  - Indices of good channels (j modalities)
%     inverse.Nd     - number of dipoles
%     inverse.pst    - peristimulus time
%     inverse.dct    - frequency range
%     inverse.F      - log-evidence
%     inverse.VE     - variance explained in spatial/temporal subspaces (%)
%     inverse.R2     - variance in subspaces accounted for by model (%)
%__________________________________________________________________________
%
%
% This version is for single subject single modality analysis and therefore
% contains none of the associated scaling factors. No symmetric priors are
% used in this implementation (just single patches) There is an option for
% a Beamforming prior : inversion type 'EBB' also added new beamforming
% method- using GS rather than ARD- from Juan David Martinez Vargas 'EBBgs'
%
%==========================================================================
% Check to see how many subjects are being fed to the function. This only
% works for a single subject.
Nl = length(D);
if Nl>1
    error('function only defined for a single subject');
end
%==========================================================================

% D - SPM data structure
%==========================================================================
if nargin > 1
    D.val = val;
elseif ~isfield(D, 'val')
    D.val = 1;
end


val=D.val;

inverse   = D.inv{val}.inverse;


%==========================================================================
% Function inputs: assign defaults if arguments are missing
%==========================================================================
try type = inverse.type;   catch, type = 'GS';     end
try hpf  = inverse.hpf;    catch, hpf  = 256;      end 
try lpf  = inverse.lpf;    catch, lpf  = 0;        end
try sdv  = inverse.sdv;    catch, sdv  = 4;        end
try Han  = inverse.Han;    catch, Han  = 1;        end
try woi  = inverse.woi;    catch, woi  = [];       end
try Nmax  = inverse.Nmax;  catch, Nmax  = 512;   end % max number of temporal modes
try Nm   = inverse.Nm;     catch, Nm   = [];       end
try Nt   = inverse.Nt;     catch, Nt   = [];       end %% fixed/requested number of temporal modes
try Ip   = inverse.Ip;     catch, Ip   = [];       end
try QE    = inverse.QE;     catch,  QE=1;          end         %  empty room noise measurement
try Qe0   = inverse.Qe0;     catch, Qe0   = exp(-5);       end  %% set noise floor at 1/100th signal power i.e. assume amplitude SNR of 10
try inverse.A;     catch, inverse.A   = [];       end %% orthogonal channel modes
try no_temporal_filter=inverse.no_temporal_filter; catch, no_temporal_filter=1; end
try complexind=inverse.complexind; catch, complexind=[]; end
%==========================================================================
% Function inputs
%==========================================================================

% Check modalities - this function only works for single modality
%--------------------------------------------------------------------------
modalities = D.inv{val}.forward.modality;
if size(modalities,1)>1
    error('not defined for multiple modalities');
end


% check lead fields and get number of dipoles (Nd) and channels (Nc)
%==========================================================================
fprintf('\nChecking leadfields')
[L,~] = spm_eeg_lgainmat(D);    % Generate/load lead field
Nd=size(L,2); % Number of sources in lead field matrix
Np = Nd; % Number of priors

% Calculate the number of priors
%==========================================================================
if ~isempty(Ip)
    Np   = length(Ip);
else
    Ip=ceil([1:Np]*Nd/Np);
end


Ic  = setdiff(D.indchantype(modalities), badchannels(D));
Nd    = size(L,2);      % Number of dipoles
fprintf(' - done\n')


%==========================================================================
% Spatial projectors: eliminate low SNR spatial modes
%==========================================================================
fprintf('Optimising and aligning spatial modes ...\n')

if isempty(inverse.A) % no spatial modes prespecified
    if isempty(Nm) %% number of modes not specifiedd
        [U,~,~]    = spm_svd((L*L'),exp(-16));
        A     = U';                 % spatial projector A
        UL    = A*L;
        
    else % number of modes pre-specified
        [U,ss,~]    = spm_svd((L*L'),0);
        if length(ss)<Nm
            disp('number available');
            length(ss)
            error('Not this many spatial modes in lead fields');
            
        end
        disp('using preselected number spatial modes !');
        A     = U(:,1:Nm)';                 % spatial projector A
        UL    = A*L;
    end
else %% U was specified in input
    disp('Using pre-specified spatial modes');
    if isempty(Nm)
        error('Need to specify number of spatial modes if U is prespecified');
    end
    
    A=inverse.A;
    UL=A*L;
end

Nm    = size(UL,1);         % Number of spatial projectors

clear ss


Is    = 1:Nd;               % Indices of active dipoles - all of them.
Ns    = length(Is);         % Number of sources, Ns
fprintf('Using %d spatial modes',Nm)
%==========================================================================
% Spatial projectors
%==========================================================================


%==========================================================================
% Time-window of interest (in milliseconds)
%==========================================================================
if isempty(woi)
    w      = 1000*[min(D.time) max(D.time)];
else
    w=woi;
end
It     = (w/1000 - D.timeonset)*D.fsample + 1;
It     = max(1,It(1)):min(It(end), length(D.time));
It     = fix(It);
disp(sprintf('Number of samples %d',length(It)))
%==========================================================================
% Time-window of interest
%==========================================================================

if ~no_temporal_filter
    
    % Peristimulus time
    %----------------------------------------------------------------------
    pst    = 1000*D.time;                   % peristimulus time (ms)
    pst    = pst(It);                       % windowed time (ms)
    dur    = (pst(end) - pst(1))/1000;      % duration (s)
    dct    = (It - It(1))/2/dur;            % DCT frequencies (Hz)
    Nb     = length(It);                    % number of time bins
    
    % Serial correlations
    %----------------------------------------------------------------------
    K      = exp(-(pst - pst(1)).^2/(2*sdv^2)); %% sdv set to 4 by default
    K      = toeplitz(K);
    qV     = sparse(K*K'); %% Samples* samples covariance matrix- assumes smooth iid
    
    % Confounds and temporal subspace
    %----------------------------------------------------------------------
    
    T      = spm_dctmtx(Nb,Nb);
    j      = find( (dct >= lpf) & (dct <= hpf) ); % This is the wrong way round but leave for now for compatibility with spm_eeg_invert
    T      = T(:,j);                    % Apply the filter to discrete cosines
    dct    = dct(j);                    % Frequencies accepted
else
    T=eye(length(It));
    qV=T;
    pst=0;dct=0;
end % if no temp filter

% Hanning window
%----------------------------------------------------------------------

if Han
    W  = sparse(1:Nb,1:Nb,spm_hanning(Nb)); % Use hanning unless specified
else
    W=1;
end


%==========================================================================
% Get trials or conditions
%==========================================================================
try
    trial = D.inv{D.val}.inverse.trials;
catch
    trial = D.condlist;
end
Ntrialtypes=length(trial);
%==========================================================================
% Get trials or conditions
%==========================================================================

%==========================================================================
% Get temporal covariance (Y'*Y) to find temporal modes
%==========================================================================
YY=0;% instantiate value of temporal covariance
N=0; % number of trials used in covariance calculation

badtrialind=D.badtrials;
Ik=[]; %% keep a record of trials used
for j = 1:Ntrialtypes                          % pool over conditions
    c     = D.indtrial(trial{j});     % and trials
    [~,ib]=intersect(c,badtrialind); % remove bad trials ib if there are any
    c=c(setxor(1:length(c),ib));
    Ik=[Ik c];
    Nk    = length(c);
    i=sqrt(-1);
    for k = 1:Nk
        if isempty(complexind)
            data=D(Ic,It,c(k));
        else
            data=squeeze(D(Ic,complexind(1,:),c(k))+i.*D(Ic,complexind(2,:),c(k)));
        end
        Y     = A*data;
        
        YY    = YY + Y'*Y;
        N     = N + 1;
    end
end
YY=YY./N;
%==========================================================================
% Get temporal covariance (Y'*Y) to find temporal modes
%==========================================================================

% Apply any Hanning and filtering
%------------------------------------------------------------------
YY         = W'*YY*W;     % Hanning
YTY         = T'*YY*T;     % Filter


%======================================================================

if isempty(Nt) %% automatically assign appropriate number of temporal modes
    [U, E]  = spm_svd(YTY,exp(-8));          % get temporal modes
    if isempty(U) %% fallback
        warning('nothing found using spm svd, using svd');
        [U, E]  = svd(YTY);          % get temporal modes
    end
    E      = diag(E)/trace(YTY);            % normalise variance
    Nr     = min(length(E),Nmax);           % number of temporal modes
    Nr=max(Nr,1); %% use at least one mode
else %% use predefined number of modes
    [U, E]  = svd(YTY);          % get temporal modes
    E      = diag(E)/trace(YTY);            % normalise variance
    disp('Fixed number of temporal modes');
    Nr=Nt;
end

V      = U(:,1:Nr);                     % temporal modes
VE     = sum(E(1:Nr));                  % variance explained

fprintf('Using %i temporal modes, ',Nr)
fprintf('accounting for %0.2f percent average variance\n',full(100*VE))

% projection and whitening
%----------------------------------------------------------------------
S      = T*V;                           % temporal projector
Vq     = S*pinv(S'*qV*S)*S';            % temporal precision



% get spatial covariance (Y*Y') for Gaussian process model
%======================================================================
% loop over Ntrialtypes trial types
%----------------------------------------------------------------------
UYYU = 0;
AYYA=0;
Nn    =0;                             % number of samples
AY={};
Ntrials=0;

for j = 1:Ntrialtypes
    UY{j} = sparse(0);
    c       = D.indtrial(trial{j});
    [~,ib]=intersect(c,badtrialind); %% remove bad trials ib if there are any
    c=c(setxor(1:length(c),ib));
    Nk    = length(c);
    % loop over epochs
    %------------------------------------------------------------------
    for k = 1:Nk
        % stack (scaled aligned data) over modalities
        %--------------------------------------------------------------
        if isempty(complexind)
            data=D(Ic,It,c(k));
        else
            data=D(Ic,complexind(1,:),c(k))+i.*D(Ic,complexind(2,:),c(k));
        end
        Y       = data*S; %% in temporal subspace
        Y=A*Y; %%  in spatial subspace
        % accumulate first & second-order responses
        %--------------------------------------------------------------
        Nn       = Nn + Nr;         % number of samples
        YY          = Y*Y';                  % and covariance
        Ntrials=Ntrials+1;
        % accumulate statistics (subject-specific)
        %--------------------------------------------------------------
        UY{j}     = UY{j} + Y;           % condition-specific ERP
        UYYU     = UYYU + YY;          % subject-specific covariance
        % and pool for optimisation of spatial priors over subjects
        %--------------------------------------------------------------
        AY{end + 1} = Y;                     % pooled response for MVB
        AYYA        = AYYA    + YY;          % pooled response for ReML
    end
end

AY=spm_cat(AY); %% goes to MVB/GS algorithm

ID    = spm_data_id(AY); %% get a unique ID for these filtered data

% assuming equal noise over subjects (Qe) and modalities AQ
%--------------------------------------------------------------------------
AQeA   = A*QE*A';           % Note that here it is A*A'
Qe{1}  = AQeA/(trace(AQeA)); % it means IID noise in virtual sensor space
Q0          = Qe0*trace(AYYA)*Qe{1}./sum(Nn); %% fixed (min) level of sensor space variance


%==========================================================================
% Step 1: Optimise spatial priors
%==========================================================================

% Create source priors (Qp)

if contains(type,'EBBr')
    reglevel=str2num(type(5:end));
    fprintf('\n Using regularizing beamformer prior to keep %d percent variance\n',reglevel)
    type='EBBr';
end
switch(type)
    case {'MSP','GS','ARD'}
        % create MSP spatial basis set in source space
        %------------------------------------------------------------------
        Qp    = {};
        LQpL  = {};
        for i = 1:Np
            q=sparse(Nd,1);
            q(i)=1;
            Qp{end + 1}.q   = q;
            LQpL{end + 1}.q = UL*q;
        end

    case {'EBBr'}
        % create SMOOTH beamforming prior.
        disp('NB  regularizing EBB algorithm, no smoothing !');
        %------------------------------------------------------------------
        QG=speye(Ns,Ns); %% was surface smoothing but not used in volumetric approach
        [u1,s1,~]=svd(AYYA);
        eigsum=cumsum(diag(s1))./sum(diag(s1));
        
        usecomp=max(find(eigsum<=reglevel/100));
        plot(1:length(s1),eigsum,usecomp,eigsum(usecomp),'*')
        fprintf('\nEBB, Keeping %d components\n',usecomp);
        ayya=u1(:,1:usecomp)'*AYYA*u1(:,1:usecomp);
        
        InvCov = spm_inv(ayya);
        allsource = sparse(Ns,1);
        Sourcepower = sparse(Ns,1);
        for bk = 1:Ns
            q               = QG(:,bk);
            
            smthlead = UL*q;     %% THIS IS WHERE THE SMOOTHNESS GETS ADDED
            smthlead=u1(:,1:usecomp)'*smthlead;
            normpower = 1/(smthlead'*smthlead);
            Sourcepower(bk) = 1/(smthlead'*InvCov*smthlead);
            allsource(bk) = Sourcepower(bk)./normpower;
        end
        allsource = allsource/max(allsource);   % Normalise
        subplot(2,1,1);
        plot(Sourcepower);
        subplot(2,1,2)
        plot(allsource);
        Qp{1} = diag(allsource);
        LQpL{1} = UL*diag(allsource)*UL';
        
    case {'EBB'}
        % create SMOOTH beamforming prior.
        disp('NB  EBB algorithm, no smoothing !');
        %------------------------------------------------------------------
        QG=speye(Ns,Ns); %% was surface smoothing but not used in volumetric approach
        
        
        InvCov = spm_inv(AYYA);
        allsource = sparse(Ns,1);
        Sourcepower = sparse(Ns,1);
        for bk = 1:Ns
            q               = QG(:,bk);
            smthlead = UL*q;     %% THIS IS WHERE THE SMOOTHNESS GETS ADDED
            normpower = 1/(smthlead'*smthlead);
            Sourcepower(bk) = 1/(smthlead'*InvCov*smthlead);
            allsource(bk) = Sourcepower(bk)./normpower;
        end
        allsource = allsource/max(allsource);   % Normalise
        Qp{1} = diag(allsource);
        LQpL{1} = UL*diag(allsource)*UL';
    case {'EBBgs'}
        % create beamforming prior- Juan David- Martinez Vargas
        %------------------------------------------------------------------
        allsource = zeros(Ntrials,Ns);
        for ii = 1:Ntrials
            InvCov = spm_inv(YYep{ii});
            Sourcepower = zeros(Ns,1);
            for bk = 1:Ns
                normpower = 1/(UL(:,bk)'*UL(:,bk));
                Sourcepower(bk) = 1/(UL(:,bk)'*InvCov*UL(:,bk));
                allsource(ii,bk) = Sourcepower(bk)./normpower;
            end
            Qp{ii}.q = allsource(ii,:);
        end
    case {'IID','MMN'}
        % create minimum norm prior
        %------------------------------------------------------------------
        Qp{1}   = speye(Ns,Ns);
        LQpL{1} = UL*UL';
end

fprintf('Using %d spatial source priors provided\n',length(Qp));


% Inverse solution
%==========================================================================
QP     = {};
LQP    = {};
LQPL   = {};

% Get source-level priors
%--------------------------------------------------------------------------
switch(type)
    case {'MSP','GS','EBBgs'}
        % Greedy search over MSPs
        %------------------------------------------------------------------
        Np    = length(Qp);
        Q     = zeros(Ns,Np); %% NB SETTING UP A NEW Q HERE
        for i = 1:Np
            Q(:,i) = Qp{i}.q;
        end
        Q = sparse(Q);
        % Multivariate Bayes (Here is performed the inversion)
        %------------------------------------------------------------------
        MVB   = spm_mvb(AY,UL,[],Q,Qe,16); %% Qe is identity with unit trace
        % Accumulate empirical priors (New set of patches for the second inversion)
        %------------------------------------------------------------------
        % MVB.cp provides the final weights of the hyperparameters
        Qcp           = Q*MVB.cp;
        QP{end + 1}   = sum(Qcp.*Q,2);
        LQP{end + 1}  = (UL*Qcp)*Q';
        LQPL{end + 1} = LQP{end}*UL';
end

switch(type)
    case {'MSP','ARD'}
        % ReML / ARD inversion
        %------------------------------------------------------------------
        [~,h,~,~] = spm_sp_reml(AYYA,[],[Qe LQpL],Nn);
        % Spatial priors (QP)
        %------------------------------------------------------------------
        % h provides the final weights of the hyperparameters
        Ne    = length(Qe);
        Np    = length(Qp);
        hp    = h(Ne + (1:Np));
        qp    = sparse(0);
        for i = 1:Np
            if hp(i) > max(hp)/128
                qp  = qp + hp(i)*Qp{i}.q*Qp{i}.q';
            end
        end
        % Accumulate empirical priors (New set of patches for the second inversion)
        %------------------------------------------------------------------
        QP{end + 1}   = diag(qp);
        LQP{end + 1}  = UL*qp;
        LQPL{end + 1} = LQP{end}*UL';
end

switch(type)
    case {'IID','MMN','EBB','EBBr'}
        % or ReML - ARD (This is where the inference on h is carried out)
        %------------------------------------------------------------------
        [~,h,~,~] = spm_reml_sc(AYYA,[],[Qe LQpL],Nn,-4,16,Q0);
        % Spatial priors (QP)
        %------------------------------------------------------------------
        % h provides the final weights of the hyperparameters
        Ne    = length(Qe);
        Np    = length(Qp);
        
        hp    = h(Ne + (1:Np));
        qp    = sparse(0);
        for i = 1:Np
            qp = qp + hp(i)*Qp{i};
        end
        % Accumulate empirical priors (New set of patches for the second inversion)
        %------------------------------------------------------------------
        QP{end + 1}   = diag(qp);
        LQP{end + 1}  = UL*qp;
        LQPL{end + 1} = LQP{end}*UL';
end


%==========================================================================
% Step 2: Re-estimate for each subject separately (fusing all modalities)
%==========================================================================

fprintf('Inverting subject 1\n')


% re-do ReML (with informative hyperpriors)
%----------------------------------------------------------------------
Np    = length(LQPL);       % Final number of priors
Ne    = length(Qe);         % Sensor noise prior


Q     = [{Q0} LQPL]; %% sensor corvariance prior:  Qe is identity with unit trace, LQPL is in the units of data

if rank(AYYA)~=size(A,1)
    warning('AYYA IS RANK DEFICIENT');
end


[Cy,h,~,F]= spm_reml_sc(AYYA,[],Q,Nn,-4,16,Q0);

% Recalculate F here
Cp    = sparse(0);
LCp   = sparse(0);
hp    = h(Ne + (1:Np));
for j = 1:Np
    Cp  =  Cp + hp(j)*QP{j};
    LCp = LCp + hp(j)*LQP{j};
end

% MAP estimates of instantaneous sources
%======================================================================
% This is equivalent to M = Cp*UL'*inv(Qe + UL*Cp*UL'))
% with Cp the posterior source covariance (with optimal h values)
M     = LCp'/Cy;

% conditional variance (leading diagonal)
% Cq    = Cp - Cp*L'*iC*L*Cp;
%----------------------------------------------------------------------
Cq    = Cp - sum(LCp.*M')';

% evaluate conditional expectation
%----------------------------------------------------------------------
% evaluate conditional expectation (of the sum over trials)
%----------------------------------------------------------------------
SSR   = 0;
SST   = 0;
J     = {};

for j = 1:Ntrialtypes
    % trial-type specific source reconstruction
    %------------------------------------------------------------------
    J{j} = M*UY{j};
    % sum of squares
    %------------------------------------------------------------------
    SSR  = SSR + sum(var((UY{j} - UL*J{j}))); %% changed variance calculation
    SST  = SST + sum(var( UY{j}));
end


% accuracy; signal to noise (over sources)
%======================================================================
R2   = 100*(SST - SSR)/SST;
fprintf('Percent variance explained %.2f (%.2f)\n',full(R2),full(R2*VE));

% Save results
%======================================================================
inverse.type   = type;                 % inverse model
inverse.M      = M;                    % MAP projector (reduced)
inverse.J      = J;                    % Conditional expectation
inverse.Y      = Y;                    % ERP data (reduced)
inverse.L      = UL;                   % Lead-field (reduced)
inverse.qC     = Cq;                   % spatial covariance
inverse.tempU  = U;                    % temporal SVD
inverse.E      = V;                    % temporal modes
inverse.qV     = Vq;                   % temporal correlations
inverse.T      = S;                    % temporal projector
inverse.U      = {A};                    % spatial projector
inverse.Is     = Is;                   % Indices of active dipoles
inverse.It     = It;                   % Indices of time bins
inverse.Ik     =Ik;                    % indices of trials used
try
    inverse.Ic{1}     = Ic;            % Indices of good channels
catch
    inverse.Ic    = Ic;                   % Indices of good channels
end
inverse.Nd     = Nd;                   % number of dipoles
inverse.pst    = pst;                  % peristimulus time
inverse.dct    = dct;                  % frequency range
inverse.F      = F;                    % log-evidence
inverse.ID     = ID;                   % data ID
inverse.R2     = R2;                   % variance explained (reduced)
inverse.VE     = R2*VE;                % variance explained
inverse.woi    = w;                    % time-window inverted
inverse.Ip     =Ip;                    % patch locations
inverse.modality = modalities;         % modalities inverted


% save in struct
%----------------------------------------------------------------------
D.inv{val}.inverse = inverse;
D.inv{val}.method  = 'Imaging';


return