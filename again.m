%% SECTION 1: Data Loading & Cleanup
clear; clc; close all;

% 1. Load Force Plate Data (.txt files)
% Ensure these files are in your current MATLAB directory
fp_quietEO = load('quietEO.txt');
fp_quietEC = load('quietEC.txt');
fp_ropeEC  = load('quietropeEC.txt');

% Fix the corrupted data: Remove the last 20 rows of zeros from the rope trial
fp_ropeEC(end-19:end, :) = []; 

% 2. Load EMG Data (.data files)
% Assuming the .data files load cleanly as numeric matrices. 
% If they contain text headers, you may need to use importdata('quieteo.data').data instead.
emg_quietEO = load('quieteo.data'); 
emg_quietEC = load('quietec.data');
emg_ropeEC  = load('quietropeec.data');

% Extract SOL (Column 1) and TA (Column 2) based on BioRadio setup
% Quiet Eyes Open
SOL_EO_raw = emg_quietEO(:, 1);
TA_EO_raw  = emg_quietEO(:, 2);

% Quiet Eyes Closed
SOL_EC_raw = emg_quietEC(:, 1);
TA_EC_raw  = emg_quietEC(:, 2);

% Perturbed Eyes Closed
SOL_rope_raw = emg_ropeEC(:, 1);
TA_rope_raw  = emg_ropeEC(:, 2);


%% SECTION 2: COP Calculation & Interpolation
dz = 41.3; % Force plate thickness in mm [cite: 133]

% Anonymous function to calculate COP X (Medial-Lateral) and Y (Anterior-Posterior)
% fp columns: 1=Fx, 2=Fy, 3=Fz, 4=Mx, 5=My, 6=Mz [cite: 132]
calc_ML = @(fp) -(fp(:,5) + fp(:,1)*dz) ./ fp(:,3); 
calc_AP = @(fp)  (fp(:,4) - fp(:,2)*dz) ./ fp(:,3);

% Calculate raw COP for all trials
ML_EO_raw = calc_ML(fp_quietEO);
AP_EO_raw = calc_AP(fp_quietEO);

ML_EC_raw = calc_ML(fp_quietEC);
AP_EC_raw = calc_AP(fp_quietEC);

ML_rope_raw = calc_ML(fp_ropeEC);
AP_rope_raw = calc_AP(fp_ropeEC);

% Interpolation Setup [cite: 143, 144]
fs_fp = 200;
fs_emg = 960;

% Function to upsample force plate data
upsample_fp = @(cop_signal) interp1(0:1/fs_fp:(length(cop_signal)-1)/fs_fp, ...
                                    cop_signal, ...
                                    0:1/fs_emg:(length(cop_signal)-1)/fs_fp, 'spline')';

% Apply interpolation
AP_EO_up = upsample_fp(AP_EO_raw);
ML_EO_up = upsample_fp(ML_EO_raw);

AP_EC_up = upsample_fp(AP_EC_raw);
ML_EC_up = upsample_fp(ML_EC_raw);

AP_rope_up = upsample_fp(AP_rope_raw);
ML_rope_up = upsample_fp(ML_rope_raw);

% %% SECTION 3: Manual Synchronization 
% % Plot 1: Eyes Open
% figure(1);
% plot(SOL_EO_raw); 
% title('Locating  Pulse: SOL Eyes Open');
% xlabel('Index'); ylabel('Amplitude');
% 
% % Plot 2: Eyes Closed
% figure(2);
% plot(SOL_EC_raw); 
% title('Locating Sync Pulse: SOL Eyes Closed');
% xlabel('Index'); ylabel('Amplitude');
% 
% % Plot 3: Perturbed (Rope)
% figure(3);
% plot(SOL_rope_raw); 
% title('Locating Sync Pulse: SOL Perturbed (Rope)');
% xlabel('Index'); ylabel('Amplitude');
% % --> ZOOM IN AND FIND THE INDEX OF THE SPIKE. 
% --> It should be around index 4800 (which is 5 seconds * 960 Hz) [cite: 154]

% STEP 2: Enter your visual findings here:
sync_idx_EO   = 9522; % REPLACE with the actual index from your EO plot
sync_idx_EC   = 8733; % REPLACE with the actual index from your EC plot
sync_idx_rope = 6628; % REPLACE with the actual index from your perturbed plot

% STEP 3: Truncate the EMG data to align with the Force Plate 
SOL_EO_sync = SOL_EO_raw(sync_idx_EO:end);
TA_EO_sync  = TA_EO_raw(sync_idx_EO:end);

SOL_EC_sync = SOL_EC_raw(sync_idx_EC:end);
TA_EC_sync  = TA_EC_raw(sync_idx_EC:end);

SOL_rope_sync = SOL_rope_raw(sync_idx_rope:end);
TA_rope_sync  = TA_rope_raw(sync_idx_rope:end);

% STEP 4: Make sure arrays are the same length before proceeding
min_len_EO = min(length(AP_EO_up), length(SOL_EO_sync));
AP_EO = AP_EO_up(1:min_len_EO); ML_EO = ML_EO_up(1:min_len_EO);
SOL_EO = SOL_EO_sync(1:min_len_EO); TA_EO = TA_EO_sync(1:min_len_EO);

% Repeat length-matching for EC and Rope trials...
min_len_EC = min(length(AP_EC_up), length(SOL_EC_sync));
AP_EC = AP_EC_up(1:min_len_EC); ML_EC = ML_EC_up(1:min_len_EC);
SOL_EC = SOL_EC_sync(1:min_len_EC); TA_EC = TA_EC_sync(1:min_len_EC);

min_len_rope = min(length(AP_rope_up), length(SOL_rope_sync));
AP_rope = AP_rope_up(1:min_len_rope); ML_rope = ML_rope_up(1:min_len_rope);
SOL_rope = SOL_rope_sync(1:min_len_rope); TA_rope = TA_rope_sync(1:min_len_rope);

%% SECTION 4: Filtering
% Filter parameters
[b, a] = butter(4, 5/(fs_emg/2), 'low');

% Function to filter COP (lowpass + zero-mean) [cite: 216, 218, 219]
filter_cop = @(signal) filtfilt(b, a, signal) - mean(filtfilt(b, a, signal));

% Function to filter EMG (rectify + lowpass, NO zero-mean) 
filter_emg = @(signal) filtfilt(b, a, abs(signal));

% Apply to EO data
AP_EO_filt = filter_cop(AP_EO);
ML_EO_filt = filter_cop(ML_EO);
SOL_EO_filt = filter_emg(SOL_EO);
TA_EO_filt  = filter_emg(TA_EO);

% Apply to EC data
AP_EC_filt = filter_cop(AP_EC);
ML_EC_filt = filter_cop(ML_EC);
SOL_EC_filt = filter_emg(SOL_EC);
TA_EC_filt  = filter_emg(TA_EC);

% Apply to Rope data
AP_rope_filt = filter_cop(AP_rope);
SOL_rope_filt = filter_emg(SOL_rope);
TA_rope_filt  = filter_emg(TA_rope);

%% SECTION 5: Quiet Standing Analysis

% Helper function for COP metrics
calc_metrics = @(sig) [mean(abs(sig)), sqrt(mean(sig.^2)), max(sig)-min(sig), sum(abs(diff(sig)))/(length(sig)/fs_emg)];
% Outputs: [MDIST, RDIST, RANGE, MVELO] [cite: 221, 222, 224, 225, 227, 229, 230]

metrics_AP_EO = calc_metrics(AP_EO_filt)
metrics_ML_EO = calc_metrics(ML_EO_filt)

% --- CROSS CORRELATION (CCF) ---
% Discard first 5 seconds (4800 pts) and segment the data [cite: 242]
seg_length = 8192; 
num_segments = 13; 
% We will define a function to do the chunking and averaging
function [avg_ccf, time_shift, CVAL, TS] = compute_ccf(emg, cop, seg_length, num_segments, fs)
    ccf_seg = zeros(num_segments, 2*seg_length - 1);
    for k = 1:num_segments 
        idx_start = (k-1)*seg_length + 4801; 
        idx_end   = idx_start + seg_length - 1;
        
        % Ensure zero mean for CCF calculation [cite: 261]
        e_seg = emg(idx_start:idx_end) - mean(emg(idx_start:idx_end));
        c_seg = cop(idx_start:idx_end) - mean(cop(idx_start:idx_end));
        
        ccf_seg(k, :) = xcorr(e_seg, c_seg, 'coeff'); 
    end
    avg_ccf = mean(ccf_seg); 
    time_shift = (-(seg_length-1):1:(seg_length-1)) / fs; 
    
    [CVAL, max_idx] = max(avg_ccf); 
    TS = time_shift(max_idx); 
end

% Compute CCF for Eyes Open
[ccf_EO_SOL, t_ccf, CVAL_EO_SOL, TS_EO_SOL] = compute_ccf(SOL_EO_filt, AP_EO_filt, seg_length, num_segments, fs_emg)
[ccf_EO_TA,  ~,     CVAL_EO_TA,  TS_EO_TA]  = compute_ccf(TA_EO_filt,  AP_EO_filt, seg_length, num_segments, fs_emg)

% Plot example
figure;
plot(t_ccf, ccf_EO_SOL);
title('Average CCF: AP vs SOL (Eyes Open)');
xlabel('Time Shift (s)'); ylabel('Correlation Coefficient (r)');

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
