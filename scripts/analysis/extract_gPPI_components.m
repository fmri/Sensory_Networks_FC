function [signal, condition_on, condition_conv, PPI] = extract_gPPI_components(signals, condition_on, HRF)
%EXTRACT_GPPI_COMPONENTS
%
% Inputs:
%       signals: rxs cell array of doubles - preprocessed fMRI BOLD signal where r is ROIs, and s is subjects
%       condition_on: cxs cell array of logicals - indicates times when a condition is
%                     active where s is subjects, and c is conditions (must have at least 2 conditions)
%       HRF: ix1 double - hemodynamic response function to use for
%            deconvolution/convolution where i is timepoints (does not have
%            to be the same as n in signal). Default is canonical HRF
%            calculated below.
%
% Outputs:
%
%
%

%% Deal with Inputs
if nargin==0 % dummy data for testing
    signals = rand(120,2,2);
    signals = squeeze(num2cell(signals,1));
    condition_on = false(120,2,2);
    condition_on(1:60,1,:) = true;
    condition_on(60:120,2,:) = true;
    condition_on = squeeze(num2cell(condition_on,1));
    t = 0:2:20; % TR
    delta = 2.05;                    % shift between the onset of each stimulus and BOLD activity
    tshift = max(t-delta,0);  % onset of the BOLD
    tau = 1.08;          % time constant
    n_HIRF = 3;          % phase delay
    BOLD_baseline = 0; % arbituary, can change according to the subject data
    HRF = ( ((tshift/tau).^(n_HIRF-1)) .* exp(-(tshift/tau)) / (tau*(factorial(n_HIRF-1))) ) + BOLD_baseline;
else
    assert(nargin>=2, 'At least 2 arguments required');
    assert(size(signals,3)==size(condition_on,3), 'Size of dimension 3 in "signal" must equal size of dimension 3 in "condition_on". These should represent number of subjects.')
    assert(size(condition_on,2)>=2, 'Input "condition_on" must have a 2nd dimension of size 2 or more (representing 2 or more conditions');
    assert(size(signal,2)>=2, 'Input "signal" must have a 2nd dimension of size 2 or more (representing 2 or more ROIs');

    if nargin < 3 || isempty(HRF) 
        % Create default canonical HRF
        t = 0:2:20; % TR
        delta = 2.05;                    % shift between the onset of each stimulus and BOLD activity
        tshift = max(t-delta,0);  % onset of the BOLD
        tau = 1.08;          % time constant
        n_HIRF = 3;          % phase delay
        BOLD_baseline = 0; % arbituary, can change according to the subject data
        HRF = ( ((tshift/tau).^(n_HIRF-1)) .* exp(-(tshift/tau)) / (tau*(factorial(n_HIRF-1))) ) + BOLD_baseline;
    end
end

Nsubjs = size(signals,2);
Nconds = size(condition_on,2);
HRF = HRF(:); % Make sure HRF is a column vector

for ss = 1:Nsubjs

    %% Convolve condition_on vectors with HRF
    cond_vecs = condition_on(:,ss); % get condition vectors for this subj
    condition_conv_subj = cellfun(@conv, cond_vecs, num2cell(repmat(HRF,1,Nconds),1), 'UniformOutput', false);
    condition_conv_subj = cellfun(@(x) condition_conv_subj(1:length(x)), cond_vecs); % cut back down to original size

    %% Create PPI vector by 1) deconvolving signal with HRF to get estimated underlying neural signal
    %%                      2) multiply condition vector with neural signal and convolve result with HRF
    signal_subj = signal(:,ss); % get signals for this subj
    % deconvolve
    deconv_signal = cellfun(@(x) deconv(x, HRF, "full", Method="least-squares"), signal_subj, 'UniformOutput', false);
    deconv_signal = cellfun(@(x) deconv_signal(1:length(x)), signal_subj); % cut back down to original size
    % multiply by each condition vector
    deconv_signal_cond = []

end


% Deal with uneven numbers of timepoints in each ROI by padding with 0s and
% padding corresponding condition vectors with 0s?


end

