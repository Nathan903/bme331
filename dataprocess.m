% Initialization parameters
fs_fp = 200; 
fs_emg = 960; 
dz = 41.3; 

% Filter design [cite: 215]
[b, a] = butter(4, 5/(fs_emg/2), 'low'); 

% File lists
fp_files = {'quietEO.txt', 'quietEC.txt', 'quietropeEC.txt'};
emg_files = {'quieteo.data', 'quietec.data', 'quietropeec.data'};
conditions = {'EO', 'EC', 'Perturbed'};

for i = 1:length(conditions) %[output:group:10c2b333]
    cond = conditions{i};

    % Load Data
    fp_data = load(fp_files{i});
    fp_data(fp_data(:,3) == 0, :) = []; % Remove padding
    
    Fx = fp_data(:, 1); Fy = fp_data(:, 2); Fz = fp_data(:, 3);
    Mx = fp_data(:, 4); My = fp_data(:, 5);
    t_fp = (0:size(fp_data, 1)-1) / fs_fp;

    emg_data = load(emg_files{i});
    SOL_raw = emg_data(:, 1); 
    TA_raw = emg_data(:, 2);
    t_emg = (0:size(emg_data, 1)-1) / fs_emg;

    % COP Calculation [cite: 133, 134, 135]
    ML_raw = -(My + Fx .* dz) ./ Fz; 
    AP_raw = (Mx - Fy .* dz) ./ Fz; 

    % Interpolation [cite: 143]
    AP_up = interp1(t_fp, AP_raw, t_emg, 'linear', 'extrap');
    ML_up = interp1(t_fp, ML_raw, t_emg, 'linear', 'extrap');

    % Synchronization truncation [cite: 146, 153]
    pulse_idx = find(diff(SOL_raw) > max(SOL_raw)*0.5, 1);
    fprintf('Sync pulse index for %s is: %d\n', cond, pulse_idx); %[output:7e73d91d] %[output:38e9e337] %[output:8d996af4]
    % Synchronization truncation (Forced Manual Indices)
    if strcmp(cond, 'EO')
        pulse_idx = 9522;
    elseif strcmp(cond, 'EC')
        pulse_idx = 8733;
    elseif strcmp(cond, 'Perturbed')
        pulse_idx = 6628;
    end

    SOL_sync = SOL_raw(pulse_idx:end);
    TA_sync = TA_raw(pulse_idx:end);
    AP_sync = AP_up(pulse_idx:end)';
    ML_sync = ML_up(pulse_idx:end)';

    min_len = min([length(SOL_sync), length(TA_sync), length(AP_sync)]);
    SOL_sync = SOL_sync(1:min_len);
    TA_sync = TA_sync(1:min_len);
    AP_sync = AP_sync(1:min_len);
    ML_sync = ML_sync(1:min_len);

    % Filtering [cite: 215, 236]
    AP_filt = filtfilt(b, a, AP_sync); 
    ML_filt = filtfilt(b, a, ML_sync);
    AP = AP_filt - mean(AP_filt); % DC removed [cite: 216]
    ML = ML_filt - mean(ML_filt); 

    SOL_rect = abs(SOL_sync); 
    TA_rect = abs(TA_sync); 
    SOL = filtfilt(b, a, SOL_rect); % DC retained [cite: 237]
    TA = filtfilt(b, a, TA_rect); 

    N = length(AP);
    T = N / fs_emg;
    time_axis = (0:N-1) / fs_emg;

    if strcmp(cond, 'EO') || strcmp(cond, 'EC')
        fprintf('\n--- %s Condition ---\n', cond); %[output:5cd31025] %[output:4609f46c]
        
        % Quiet Standing Q1: COP Measures [cite: 221]
        MDIST_AP = mean(abs(AP)); 
        RDIST_AP = sqrt(mean(AP.^2)); 
        RANGE_AP = max(AP) - min(AP); 
        MVELO_AP = sum(abs(diff(AP))) / T; 

        MDIST_ML = mean(abs(ML));
        RDIST_ML = sqrt(mean(ML.^2));
        RANGE_ML = max(ML) - min(ML);
        MVELO_ML = sum(abs(diff(ML))) / T;

        fprintf('AP Measures: MDIST=%.2f, RDIST=%.2f, RANGE=%.2f, MVELO=%.2f\n', MDIST_AP, RDIST_AP, RANGE_AP, MVELO_AP); %[output:6e05e305] %[output:1d3301f7]
        fprintf('ML Measures: MDIST=%.2f, RDIST=%.2f, RANGE=%.2f, MVELO=%.2f\n', MDIST_ML, RDIST_ML, RANGE_ML, MVELO_ML); %[output:87b6f4d3] %[output:1a93c8a8]

        % Plots: COP y vs x [cite: 403]
        fig1 = figure('Name', ['COP XY - ', cond]); %[output:123be958] %[output:7cd1cfe3]
        plot(ML, AP); %[output:123be958] %[output:7cd1cfe3] %[output:1a24ad0d]
        xlabel('ML Displacement'); ylabel('AP Displacement');
        title(['COP Displacement (y vs x) - ', cond]);
        saveas(fig1, sprintf('COP_XY_%s.png', cond));

        % Plots: TA and SOL EMG
        fig2 = figure('Name', ['EMG - ', cond]); %[output:086e977c] %[output:14cc2612]
        subplot(2,1,1); plot(time_axis, TA); title(['TA EMG - ', cond]); ylabel('Amplitude'); %[output:086e977c] %[output:14cc2612]
        subplot(2,1,2); plot(time_axis, SOL); title(['SOL EMG - ', cond]); xlabel('Time (s)'); ylabel('Amplitude');
        saveas(fig2, sprintf('EMG_%s.png', cond));

        % Quiet Standing Q2: Cross-Correlation [cite: 238, 242]
        seg_length = 8192; 
        num_segs = 13; 
        start_idx = 4801; 

        ccf_seg_SOL = zeros(num_segs, 2*seg_length - 1);
        ccf_seg_TA = zeros(num_segs, 2*seg_length - 1);

        for k = 1:num_segs
            idx_a = (k-1)*seg_length + start_idx; 
            idx_b = idx_a + seg_length - 1; 

            ccf_seg_SOL(k,:) = xcorr(SOL(idx_a:idx_b) - mean(SOL(idx_a:idx_b)), AP(idx_a:idx_b) - mean(AP(idx_a:idx_b)), 'coeff'); 
            ccf_seg_TA(k,:) = xcorr(TA(idx_a:idx_b) - mean(TA(idx_a:idx_b)), AP(idx_a:idx_b) - mean(AP(idx_a:idx_b)), 'coeff');
        end

        avg_ccf_SOL = mean(ccf_seg_SOL); %[cite: 244]
        avg_ccf_TA = mean(ccf_seg_TA);
        time_lag = (-(seg_length-1):1:(seg_length-1)) / fs_emg; 

        [CVAL_SOL, max_idx_SOL] = max(avg_ccf_SOL); %[cite: 245]
        TS_SOL = time_lag(max_idx_SOL); 

        [CVAL_TA, max_idx_TA] = max(avg_ccf_TA);
        TS_TA = time_lag(max_idx_TA);

        fprintf('CCF AP x SOL: CVAL=%.4f, TS=%.4fs\n', CVAL_SOL, TS_SOL); %[output:36599289] %[output:63953f97]
        fprintf('CCF AP x TA: CVAL=%.4f, TS=%.4fs\n', CVAL_TA, TS_TA); %[output:37f86e7f] %[output:6b494e87]

        % Plots: Average CCF [cite: 282]
        fig3 = figure('Name', ['CCF - ', cond]); %[output:8cfe3824] %[output:2777df80]
        subplot(2,1,1); plot(time_lag, avg_ccf_SOL); title(['Average CCF: AP vs SOL - ', cond]); xlabel('Time Shift (s)'); ylabel('r'); %[output:8cfe3824] %[output:2777df80]
        subplot(2,1,2); plot(time_lag, avg_ccf_TA); title(['Average CCF: AP vs TA - ', cond]); xlabel('Time Shift (s)'); ylabel('r');
        saveas(fig3, sprintf('CCF_%s.png', cond));

    elseif strcmp(cond, 'Perturbed')
        fprintf('\n--- Perturbed Condition ---\n'); %[output:8d0bc825]
        
        % Perturbed Standing Q1: Aligned Plot [cite: 407]
        fig4 = figure('Name', 'Perturbed Standing Aligned'); %[output:68f742f4]
        ax1 = subplot(3,1,1); plot(time_axis, AP); title('AP Displacement'); ylabel('Displacement'); %[output:68f742f4] %[output:35327de6]
        ax2 = subplot(3,1,2); plot(time_axis, SOL); title('SOL EMG'); ylabel('Amplitude');
        ax3 = subplot(3,1,3); plot(time_axis, TA); title('TA EMG'); xlabel('Time (s)'); ylabel('Amplitude');
        linkaxes([ax1, ax2, ax3], 'x');
        saveas(fig4, 'Perturbed_Aligned.png'); %[output:3aa44002]

        % Perturbed Standing Q2: Visual Inspection & Differentials 
        disp('Select 3 perturbation instances on the SOL plot.');
        disp('Click precisely on the ONSET of the SOL response for each of the 3 pulls.');
        [sol_t, ~] = ginput(3);
        
        disp('Select 3 perturbation instances on the TA plot.');
        disp('Click precisely on the ONSET of the TA response for each of the 3 pulls.');
        [ta_t, ~] = ginput(3);

        ap_max_vals = zeros(1,3);
        t_ap_max = zeros(1,3);
        sol_diff = zeros(1,3);
        ta_diff = zeros(1,3);

        for p = 1:3
            % Define a search window for APmax: 0 to 2 seconds after SOL onset
            window_start = round(sol_t(p) * fs_emg);
            window_end = window_start + (2 * fs_emg); 
            
            % Bound check
            if window_end > N; window_end = N; end
            
            [max_val, max_idx_rel] = max(AP(window_start:window_end));
            ap_max_vals(p) = max_val;
            t_ap_max(p) = time_axis(window_start + max_idx_rel - 1);
            
            sol_diff(p) = sol_t(p) - t_ap_max(p); 
            ta_diff(p) = ta_t(p) - t_ap_max(p);
        end

        mean_APmax = mean(ap_max_vals);
        mean_SOL_diff = mean(sol_diff);
        mean_TA_diff = mean(ta_diff);

        fprintf('Mean APmax: %.4f\n', mean_APmax);
        fprintf('Mean Time Differential (SOL onset - APmax): %.4fs\n', mean_SOL_diff);
        fprintf('Mean Time Differential (TA onset - APmax): %.4fs\n', mean_TA_diff);
    end
end %[output:group:10c2b333]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":40}
%---
%[output:7e73d91d]
%   data: {"dataType":"text","outputData":{"text":"Sync pulse index for EO is: 106\n","truncated":false}}
%---
%[output:5cd31025]
%   data: {"dataType":"text","outputData":{"text":"\n--- EO Condition ---\n","truncated":false}}
%---
%[output:6e05e305]
%   data: {"dataType":"text","outputData":{"text":"AP Measures: MDIST=11.28, RDIST=21.28, RANGE=116.70, MVELO=1.13\n","truncated":false}}
%---
%[output:87b6f4d3]
%   data: {"dataType":"text","outputData":{"text":"ML Measures: MDIST=2.15, RDIST=4.06, RANGE=22.35, MVELO=0.28\n","truncated":false}}
%---
%[output:123be958]
%   data: {"dataType":"image","outputData":{"dataUri":"data:,","height":0,"width":0}}
%---
%[output:086e977c]
%   data: {"dataType":"image","outputData":{"dataUri":"data:,","height":0,"width":0}}
%---
%[output:36599289]
%   data: {"dataType":"text","outputData":{"text":"CCF AP x SOL: CVAL=0.0615, TS=-1.0667s\n","truncated":false}}
%---
%[output:37f86e7f]
%   data: {"dataType":"text","outputData":{"text":"CCF AP x TA: CVAL=0.0674, TS=1.2521s\n","truncated":false}}
%---
%[output:8cfe3824]
%   data: {"dataType":"image","outputData":{"dataUri":"data:,","height":0,"width":0}}
%---
%[output:38e9e337]
%   data: {"dataType":"text","outputData":{"text":"Sync pulse index for EC is: 716\n","truncated":false}}
%---
%[output:8d996af4]
%   data: {"dataType":"text","outputData":{"text":"Sync pulse index for Perturbed is: 47\n","truncated":false}}
%---
%[output:4609f46c]
%   data: {"dataType":"text","outputData":{"text":"\n--- EC Condition ---\n","truncated":false}}
%---
%[output:1d3301f7]
%   data: {"dataType":"text","outputData":{"text":"AP Measures: MDIST=10.84, RDIST=20.66, RANGE=114.76, MVELO=1.24\n","truncated":false}}
%---
%[output:1a93c8a8]
%   data: {"dataType":"text","outputData":{"text":"ML Measures: MDIST=0.23, RDIST=0.44, RANGE=2.54, MVELO=0.17\n","truncated":false}}
%---
%[output:7cd1cfe3]
%   data: {"dataType":"image","outputData":{"dataUri":"data:,","height":0,"width":0}}
%---
%[output:1a24ad0d]
%   data: {"dataType":"warning","outputData":{"text":"Warning: An error occurred while drawing the scene: Could not find node in peer tree during replaceCamera"}}
%---
%[output:14cc2612]
%   data: {"dataType":"image","outputData":{"dataUri":"data:,","height":0,"width":0}}
%---
%[output:2777df80]
%   data: {"dataType":"image","outputData":{"dataUri":"data:,","height":0,"width":0}}
%---
%[output:6b494e87]
%   data: {"dataType":"text","outputData":{"text":"CCF AP x TA: CVAL=0.0888, TS=-0.3833s\n","truncated":false}}
%---
%[output:8d0bc825]
%   data: {"dataType":"text","outputData":{"text":"\n--- Perturbed Condition ---\n","truncated":false}}
%---
%[output:68f742f4]
%   data: {"dataType":"image","outputData":{"dataUri":"data:,","height":0,"width":0}}
%---
%[output:3aa44002]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Error using <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('print', 'C:\\Program Files\\MATLAB\\R2025b\\toolbox\\matlab\\graphics\\graphics\\printing\\print.m', 86)\" style=\"font-weight:bold\">print<\/a> (<a href=\"matlab: opentoline('C:\\Program Files\\MATLAB\\R2025b\\toolbox\\matlab\\graphics\\graphics\\printing\\print.m',86,0)\">line 86<\/a>)\nThere was a problem while generating the output: Export failed.\n\nError in <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('saveas', 'C:\\Program Files\\MATLAB\\R2025b\\toolbox\\matlab\\general\\saveas.m', 181)\" style=\"font-weight:bold\">saveas<\/a> (<a href=\"matlab: opentoline('C:\\Program Files\\MATLAB\\R2025b\\toolbox\\matlab\\general\\saveas.m',181,0)\">line 181<\/a>)\n        print( h, name, ['-d' dev{i}] )\n        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"}}
%---
%[output:63953f97]
%   data: {"dataType":"text","outputData":{"text":"CCF AP x SOL: CVAL=0.0699, TS=-2.2104s\n","truncated":false}}
%---
%[output:35327de6]
%   data: {"dataType":"warning","outputData":{"text":"Warning: An error occurred while drawing the scene: Error while executing frame: TypeError: Cannot read properties of null (reading 'getMaxTextureSize')"}}
%---
