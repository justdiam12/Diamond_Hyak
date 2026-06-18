%%% Directories %%%
main_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Hyak';
data_dir = fullfile(main_dir, "data");
bellhop_save_dir = fullfile(main_dir, 'save', 'BELLHOP');
save_dir = fullfile(main_dir, 'save', 'plots');
depths_dir = fullfile(main_dir, 'data', 'array_depths');
depth_files = dir(fullfile(depths_dir, "*.mat"));


%%% Load Array Depths %%%
S = load(fullfile(depths_dir, 'depths_20200121.mat'));
depths_21 = sort(S.unique_depths);
S = load(fullfile(depths_dir, 'depths_20200124.mat'));
depths_24 = sort(S.unique_depths);
S = load(fullfile(depths_dir, 'depths_20200128.mat'));
depths_28 = sort(S.unique_depths);
S = load(fullfile(depths_dir, 'depths_20200129.mat'));
depths_29 = sort(S.unique_depths);
rd = sort([depths_24 depths_28 depths_29]);


%%% Load the Tide Cycle Array %%%
load(fullfile(data_dir, "tide_height.mat"))


%%% Adjust Array Depths based on Time/Depth Reference %%%
% January 21 and 24, 2020
start_d_21_24 = 42.261;  % Top Hydrophone Depth (m) - Grant's logs
depth_adj_21_24 = start_d_21_24 - depths_21(2);
start_t_21_24 = datetime(2020, 1, 21, 13, 34, 0); % Timestamp from Grant
tide_ref_21_24 = tide_height(find(arms_time == start_t_21_24));

% January 28 and 29, 2020
start_t_28 = datetime(2020, 1, 27, 15, 55, 0); % Time in DJ's Logbook
tide_ref_28 = tide_height(find(arms_time == start_t_28));
start_t_29 = datetime(2020, 1, 28, 15, 24, 0); % Time in DJ's Logbook
tide_ref_29 = tide_height(find(arms_time == start_t_29));


%%% Model Parameters %%%
freqs = 3450:3550;
time = 0:0.0001:2;
source_level = 200; 


%%% Source Signal %%%
% Source Signal Parameters
center_freq = 3500;
fs = 80000;
f0 = center_freq - 50;
f1 = center_freq + 50;

% Source Signal and Time Array
source_time = 0:1/fs:1-1/fs;
pulse = chirp(source_time, ...
              f0, ...
              1, ...
              f1, ...
              'linear');

window = tukeywin(length(pulse), 0.1)';
source_signal = pulse .* window;

%%% Sift Through RAMPE Save Files %%%
bellhop_save_files = dir(fullfile(bellhop_save_dir, "*.mat"));
bellhop_21_files = {};
bellhop_24_files = {};
bellhop_28_files = {};
bellhop_29_files = {};
for f = 1:length(bellhop_save_files)
    if contains(bellhop_save_files(f).name, '20200121')
        bellhop_21_files{end+1} = bellhop_save_files(f).name;
    elseif contains(bellhop_save_files(f).name, '20200124')
        bellhop_24_files{end+1} = bellhop_save_files(f).name;
    elseif contains(bellhop_save_files(f).name, '20200128')
        bellhop_28_files{end+1} = bellhop_save_files(f).name;
    elseif contains(bellhop_save_files(f).name, '20200129')
        bellhop_29_files{end+1} = bellhop_save_files(f).name;
    end
end
bellhop_21_files = sort(bellhop_21_files);
bellhop_24_files = sort(bellhop_24_files);
bellhop_28_files = sort(bellhop_28_files);
bellhop_29_files = sort(bellhop_29_files);


%%% TL and RL at Receiver Depths (Standard Depths) %%%
bellhop_tl_depth_21 = zeros(length(rd), length(bellhop_21_files));
bellhop_tl_depth_24 = zeros(length(rd), length(bellhop_24_files));
bellhop_tl_depth_28 = zeros(length(rd), length(bellhop_28_files));
bellhop_tl_depth_29 = zeros(length(rd), length(bellhop_29_files));
bellhop_rl_depth_21 = zeros(length(rd), length(bellhop_21_files));
bellhop_rl_depth_24 = zeros(length(rd), length(bellhop_24_files));
bellhop_rl_depth_28 = zeros(length(rd), length(bellhop_28_files));
bellhop_rl_depth_29 = zeros(length(rd), length(bellhop_29_files));
rl_depth_avg_21 = zeros(length(bellhop_21_files), 1);
rl_depth_avg_24 = zeros(length(bellhop_24_files), 1);
rl_depth_avg_28 = zeros(length(bellhop_28_files), 1);
rl_depth_avg_29 = zeros(length(bellhop_29_files), 1);

%%% TL and RL at Receiver Depths (Tidal Cycle Adjusted Depths) %%%
bellhop_tl_depth_21_tidal = zeros(length(rd), length(bellhop_21_files));
bellhop_tl_depth_24_tidal = zeros(length(rd), length(bellhop_24_files));
bellhop_tl_depth_28_tidal = zeros(length(rd), length(bellhop_28_files));
bellhop_tl_depth_29_tidal = zeros(length(rd), length(bellhop_29_files));
bellhop_rl_depth_21_tidal = zeros(length(rd), length(bellhop_21_files));
bellhop_rl_depth_24_tidal = zeros(length(rd), length(bellhop_24_files));
bellhop_rl_depth_28_tidal = zeros(length(rd), length(bellhop_28_files));
bellhop_rl_depth_29_tidal = zeros(length(rd), length(bellhop_29_files));
rl_depth_avg_21_tidal = zeros(length(bellhop_21_files), 1);
rl_depth_avg_24_tidal = zeros(length(bellhop_24_files), 1);
rl_depth_avg_28_tidal = zeros(length(bellhop_28_files), 1);
rl_depth_avg_29_tidal = zeros(length(bellhop_29_files), 1);


%%% Time Array %%%
time_21 = {};
time_24 = {};
time_28 = {};
time_29 = {};


%%% Loop through RAMPE Save Files %%%
% January 28, 2020
for s = 1:length(bellhop_28_files)

    type = 1;
    disp("Working on file " + bellhop_28_files(s))

    % Current Datetime
    timestamp_str = erase(bellhop_28_files(s), {'bellhop_time_', '.mat'});
    t_curr = datetime(timestamp_str, 'InputFormat', 'yyyyMMdd_HHmmss');
    time_28{s} = t_curr;

    % Tide Height on Specific Date
    t_idx = find(arms_time == t_curr);
    tide_curr = tide_height(t_idx);

    % Adjust Receiver Depths
    depth_diff =  tide_ref_28 - tide_curr;
    depths_28_tidal = rd - depth_diff;

    % Extract Data
    data = load(fullfile(bellhop_save_dir, bellhop_28_files{s}));
    arrivals = data.arrivals;

    % TL at standard depths (No Tidal Cycle Adjustment)
    tl = tl_vs_depth_bh(arrivals, freqs, time, ...
                     rd, source_signal, source_time);
    bellhop_tl_depth_28(:,s) = tl;
    rl_depth_avg_28(s) = mean(source_level - tl);

    % TL with Tidal Cycle Adjustment
    tl_tidal = tl_vs_depth_bh(arrivals, freqs, time, ...
                     depths_28_tidal, source_signal, source_time);
    bellhop_tl_depth_28_tidal(:,s) = tl_tidal;
    bellhop_rl_depth_28_tidal(:,s) = source_level - tl_tidal;
    rl_depth_avg_28_tidal(s) = mean(source_level - tl_tidal);

end