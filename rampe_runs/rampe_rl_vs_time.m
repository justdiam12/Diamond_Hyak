%%% Directories %%%
main_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Hyak';
data_dir = fullfile(main_dir, "data");
rampe_save_dir = fullfile(main_dir, 'save', 'RAMPE');
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
rampe_save_files = dir(fullfile(rampe_save_dir, "*.mat"));
rampe_21_files = {};
rampe_24_files = {};
rampe_28_files = {};
rampe_29_files = {};
for f = 1:length(rampe_save_files)
    if contains(rampe_save_files(f).name, '20200121')
        rampe_21_files{end+1} = rampe_save_files(f).name;
    elseif contains(rampe_save_files(f).name, '20200124')
        rampe_24_files{end+1} = rampe_save_files(f).name;
    elseif contains(rampe_save_files(f).name, '20200128')
        rampe_28_files{end+1} = rampe_save_files(f).name;
    elseif contains(rampe_save_files(f).name, '20200129')
        rampe_29_files{end+1} = rampe_save_files(f).name;
    end
end
rampe_21_files = sort(rampe_21_files);
rampe_24_files = sort(rampe_24_files);
rampe_28_files = sort(rampe_28_files);
rampe_29_files = sort(rampe_29_files);


%%% TL and RL at Receiver Depths (Standard Depths) %%%
rampe_tl_depth_21 = zeros(length(depths_21), length(rampe_21_files));
rampe_tl_depth_24 = zeros(length(depths_24), length(rampe_24_files));
rampe_tl_depth_28 = zeros(length(depths_28), length(rampe_28_files));
rampe_tl_depth_29 = zeros(length(depths_29), length(rampe_29_files));
rampe_rl_depth_21 = zeros(length(depths_21), length(rampe_21_files));
rampe_rl_depth_24 = zeros(length(depths_24), length(rampe_24_files));
rampe_rl_depth_28 = zeros(length(depths_28), length(rampe_28_files));
rampe_rl_depth_29 = zeros(length(depths_29), length(rampe_29_files));
rl_depth_avg_21 = zeros(length(rampe_21_files), 1);
rl_depth_avg_24 = zeros(length(rampe_24_files), 1);
rl_depth_avg_28 = zeros(length(rampe_28_files), 1);
rl_depth_avg_29 = zeros(length(rampe_29_files), 1);

%%% TL and RL at Receiver Depths (Tidal Cycle Adjusted Depths) %%%
rampe_tl_depth_21_tidal = zeros(length(depths_21), length(rampe_21_files));
rampe_tl_depth_24_tidal = zeros(length(depths_24), length(rampe_24_files));
rampe_tl_depth_28_tidal = zeros(length(depths_28), length(rampe_28_files));
rampe_tl_depth_29_tidal = zeros(length(depths_29), length(rampe_29_files));
rampe_rl_depth_21_tidal = zeros(length(depths_21), length(rampe_21_files));
rampe_rl_depth_24_tidal = zeros(length(depths_24), length(rampe_24_files));
rampe_rl_depth_28_tidal = zeros(length(depths_28), length(rampe_28_files));
rampe_rl_depth_29_tidal = zeros(length(depths_29), length(rampe_29_files));
rl_depth_avg_21_tidal = zeros(length(rampe_21_files), 1);
rl_depth_avg_24_tidal = zeros(length(rampe_24_files), 1);
rl_depth_avg_28_tidal = zeros(length(rampe_28_files), 1);
rl_depth_avg_29_tidal = zeros(length(rampe_29_files), 1);


%%% Time Array %%%
time_21 = {};
time_24 = {};
time_28 = {};
time_29 = {};


%%% Loop through RAMPE Save Files %%%
% January 21, 2020
for s = 1:length(rampe_21_files)

    type = 1;
    disp("Working on file " + rampe_21_files(s))

    % Current Datetime
    timestamp_str = erase(rampe_21_files(s), {'rampe_time_', '.mat'});
    t_curr = datetime(timestamp_str, 'InputFormat', 'yyyyMMdd_HHmmss');
    time_21{s} = t_curr;

    % Tide Height on Specific Date
    t_idx = find(arms_time == t_curr);
    tide_curr = tide_height(t_idx);

    % Adjust Receiver Depths
    depth_diff =  tide_ref_21_24 - tide_curr;
    depths_21_tidal = depths_21 - depth_diff + depth_adj_21_24;

    % Extract Data
    data = load(fullfile(rampe_save_dir, rampe_21_files{s}));
    Depth = data.Depth - 150;
    Depth_mask = (Depth >= -20 & Depth <= 200);
    Depth = Depth(Depth_mask);
    Psv = data.Psv(Depth_mask,:);

    % TL at standard depths (No Tidal Cycle Adjustment)
    tl = tl_vs_depth_ram(Psv, Depth, ...
                     freqs, time, depths_21, ...
                     source_signal, source_time, type);
    rampe_tl_depth_21(:,s) = tl;
    rampe_rl_depth_21(:,s) = source_level - tl;
    rl_depth_avg_21(s) = mean(source_level - tl);

    % TL with Tidal Cycle Adjustment
    tl_tidal = tl_vs_depth_ram(Psv, Depth, ...
                           freqs, time, depths_21_tidal, ...
                           source_signal, source_time, type);
    rampe_tl_depth_21_tidal(:,s) = tl_tidal;
    rampe_rl_depth_21_tidal(:,s) = source_level - tl_tidal;
    rl_depth_avg_21_tidal(s) = mean(source_level - tl_tidal);


end

% January 24, 2020
for s = 1:length(rampe_24_files)

    type = 1;
    disp("Working on file " + rampe_24_files(s))

    % Current Datetime
    timestamp_str = erase(rampe_24_files(s), {'rampe_time_', '.mat'});
    t_curr = datetime(timestamp_str, 'InputFormat', 'yyyyMMdd_HHmmss');
    time_24{s} = t_curr;

    % Tide Height on Specific Date
    t_idx = find(arms_time == t_curr);
    tide_curr = tide_height(t_idx);

    % Adjust Receiver Depths
    depth_diff =  tide_ref_21_24 - tide_curr;
    depths_24_tidal = depths_24 - depth_diff + depth_adj_21_24;

    % Extract Data
    data = load(fullfile(rampe_save_dir, rampe_24_files{s}));
    Depth = data.Depth - 150;
    Depth_mask = (Depth >= -20 & Depth <= 200);
    Depth = Depth(Depth_mask);
    Psv = data.Psv(Depth_mask,:);

    % TL at standard depths (No Tidal Cycle Adjustment)
    tl = tl_vs_depth_ram(Psv, Depth, ...
                     freqs, time, depths_24, ...
                     source_signal, source_time, type);
    rampe_tl_depth_24(:,s) = tl;
    rampe_rl_depth_24(:,s) = source_level - tl;
    rl_depth_avg_24(s) = mean(source_level - tl);

    % TL with Tidal Cycle Adjustment
    tl_tidal = tl_vs_depth_ram(Psv, Depth, ...
                           freqs, time, depths_24_tidal, ...
                           source_signal, source_time, type);
    rampe_tl_depth_24_tidal(:,s) = tl_tidal;
    rampe_rl_depth_24_tidal(:,s) = source_level - tl_tidal;
    rl_depth_avg_24_tidal(s) = mean(source_level - tl_tidal);

end

% January 28, 2020
for s = 1:length(rampe_28_files)

    type = 1;
    disp("Working on file " + rampe_28_files(s))

    % Current Datetime
    timestamp_str = erase(rampe_28_files(s), {'rampe_time_', '.mat'});
    t_curr = datetime(timestamp_str, 'InputFormat', 'yyyyMMdd_HHmmss');
    time_28{s} = t_curr;

    % Tide Height on Specific Date
    t_idx = find(arms_time == t_curr);
    tide_curr = tide_height(t_idx);

    % Adjust Receiver Depths
    depth_diff =  tide_ref_28 - tide_curr;
    depths_28_tidal = depths_28 - depth_diff;

    % Extract Data
    data = load(fullfile(rampe_save_dir, rampe_28_files{s}));
    Depth = data.Depth - 150;
    Depth_mask = (Depth >= -20 & Depth <= 200);
    Depth = Depth(Depth_mask);
    Psv = data.Psv(Depth_mask,:);

    % TL at standard depths (No Tidal Cycle Adjustment)
    tl = tl_vs_depth_ram(Psv, Depth, ...
                     freqs, time, depths_28, ...
                     source_signal, source_time, type);
    rampe_tl_depth_28(:,s) = tl;
    rampe_rl_depth_28(:,s) = source_level - tl;
    rl_depth_avg_28(s) = mean(source_level - tl);

    % TL with Tidal Cycle Adjustment
    tl_tidal = tl_vs_depth_ram(Psv, Depth, ...
                           freqs, time, depths_28_tidal, ...
                           source_signal, source_time, type);
    rampe_tl_depth_28_tidal(:,s) = tl_tidal;
    rampe_rl_depth_28_tidal(:,s) = source_level - tl_tidal;
    rl_depth_avg_28_tidal(s) = mean(source_level - tl_tidal);

end

% January 29, 2020
for s = 1:length(rampe_29_files)

    type = 1;
    disp("Working on file " + rampe_29_files(s))

    % Current Datetime
    timestamp_str = erase(rampe_29_files(s), {'rampe_time_', '.mat'});
    t_curr = datetime(timestamp_str, 'InputFormat', 'yyyyMMdd_HHmmss');
    time_29{s} = t_curr;

    % Tide Height on Specific Date
    t_idx = find(arms_time == t_curr);
    tide_curr = tide_height(t_idx);

    % Adjust Receiver Depths
    depth_diff =  tide_ref_29 - tide_curr;
    depths_29_tidal = depths_29 - depth_diff;

    % Extract Data
    data = load(fullfile(rampe_save_dir, rampe_29_files{s}));
    Depth = data.Depth - 150;
    Depth_mask = (Depth >= -20 & Depth <= 200);
    Depth = Depth(Depth_mask);
    Psv = data.Psv(Depth_mask,:);

    % TL at standard depths (No Tidal Cycle Adjustment)
    tl = tl_vs_depth_ram(Psv, Depth, ...
                     freqs, time, depths_29, ...
                     source_signal, source_time, type);
    rampe_tl_depth_29(:,s) = tl;
    rampe_rl_depth_29(:,s) = source_level - tl;
    rl_depth_avg_29(s) = mean(source_level - tl);

    % TL with Tidal Cycle Adjustment
    tl_tidal = tl_vs_depth_ram(Psv, Depth, ...
                           freqs, time, depths_29_tidal, ...
                           source_signal, source_time, type);
    rampe_tl_depth_29_tidal(:,s) = tl_tidal;
    rampe_rl_depth_29_tidal(:,s) = source_level - tl_tidal;
    rl_depth_avg_29(s) = mean(source_level - tl_tidal);

end


%%% Mean and STDDEV of TL (Standard Depths) %%%
rampe_avg_21 = zeros(length(depths_21), 1);
rampe_avg_24 = zeros(length(depths_24), 1);
rampe_avg_28 = zeros(length(depths_28), 1);
rampe_avg_29 = zeros(length(depths_29), 1);
rampe_std_21 = zeros(length(depths_21), 1);
rampe_std_24 = zeros(length(depths_24), 1);
rampe_std_28 = zeros(length(depths_28), 1);
rampe_std_29 = zeros(length(depths_29), 1);

%%% Mean and STDDEV of TL (Tidal Cycle adjusted Depths) %%%
rampe_avg_21_tidal = zeros(length(depths_21), 1);
rampe_avg_24_tidal = zeros(length(depths_24), 1);
rampe_avg_28_tidal = zeros(length(depths_28), 1);
rampe_avg_29_tidal = zeros(length(depths_29), 1);
rampe_std_21_tidal = zeros(length(depths_21), 1);
rampe_std_24_tidal = zeros(length(depths_24), 1);
rampe_std_28_tidal = zeros(length(depths_28), 1);
rampe_std_29_tidal = zeros(length(depths_29), 1);

for d = 1:length(depths_21)
    % Standard Depths
    rampe_avg_21(d) = mean(rampe_tl_depth_21(d,:));
    rampe_avg_24(d) = mean(rampe_tl_depth_24(d,:));
    rampe_avg_28(d) = mean(rampe_tl_depth_28(d,:));
    rampe_avg_29(d) = mean(rampe_tl_depth_29(d,:));
    rampe_std_21(d) = std(rampe_tl_depth_21(d,:));
    rampe_std_24(d) = std(rampe_tl_depth_24(d,:));
    rampe_std_28(d) = std(rampe_tl_depth_28(d,:));
    rampe_std_29(d) = std(rampe_tl_depth_29(d,:));

    % Tidal Cycle Adjusted Depths %
    rampe_avg_21_tidal(d) = mean(rampe_tl_depth_21(d,:));
    rampe_avg_24_tidal(d) = mean(rampe_tl_depth_24(d,:));
    rampe_avg_28_tidal(d) = mean(rampe_tl_depth_28(d,:));
    rampe_avg_29_tidal(d) = mean(rampe_tl_depth_29(d,:));
    rampe_std_21_tidal(d) = std(rampe_tl_depth_21(d,:));
    rampe_std_24_tidal(d) = std(rampe_tl_depth_24(d,:));
    rampe_std_28_tidal(d) = std(rampe_tl_depth_28(d,:));
    rampe_std_29_tidal(d) = std(rampe_tl_depth_29(d,:));
end


%%% Plot TL vs Depth Figure %%%
fig1 = figure('Position',[100 100 800 1400]);

subplot(1,2,1)
hold on
errorbar(rampe_avg_21, depths_21, rampe_std_21, 'horizontal', ...
    'bo', 'DisplayName', 'January 21', 'LineWidth', 2)
errorbar(rampe_avg_24, depths_24, rampe_std_24, 'horizontal', ...
    'ro', 'DisplayName', 'January 24', 'LineWidth', 2)
errorbar(rampe_avg_28, depths_28, rampe_std_28, 'horizontal', ...
    'go', 'DisplayName', 'January 28', 'LineWidth', 2)
errorbar(rampe_avg_29, depths_29, rampe_std_29, 'horizontal', ...
    'mo', 'DisplayName', 'January 29', 'LineWidth', 2)

set(gca, 'YDir', 'reverse')
set(gca, 'FontSize', 12)
title('No Tidal-Cycle Depth Adjustment', fontsize=12)
xlabel('Average TL (dB re 1 µPa)', fontsize=12)
ylabel('Depth (m)', fontsize=12)
xlim([60 80])
grid on
hold off

subplot(1,2,2)
hold on
errorbar(rampe_avg_21_tidal, depths_21, rampe_std_21_tidal, 'horizontal', ...
    'bo', 'DisplayName', 'January 21', 'LineWidth', 2)
errorbar(rampe_avg_24_tidal, depths_24, rampe_std_24_tidal, 'horizontal', ...
    'ro', 'DisplayName', 'January 24', 'LineWidth', 2)
errorbar(rampe_avg_28_tidal, depths_28, rampe_std_28_tidal, 'horizontal', ...
    'go', 'DisplayName', 'January 28', 'LineWidth', 2)
errorbar(rampe_avg_29_tidal, depths_29, rampe_std_29_tidal, 'horizontal', ...
    'mo', 'DisplayName', 'January 29', 'LineWidth', 2)

set(gca, 'YDir', 'reverse')
set(gca, 'FontSize', 12)
title('Tidal-Cycle Depth Adjusted', fontsize=12)
xlabel('Average TL (dB re 1 µPa)', fontsize=12)
ylabel('Depth (m)', fontsize=12)
xlim([60 80])
legend('Location','northeast','FontSize', 10)
grid on
hold off

exportgraphics(fig1, fullfile(save_dir, "tl_depth_st_and_tidal.png"), "Resolution", 300)
disp("RAMPE TL Plot Complete.")



% %%% Plot RL vs Depth Figure %%%
% Depth Indices Closest to 55m
[~, d_idx_21] = min(abs(depths_21 - 55));
[~, d_idx_24] = min(abs(depths_24 - 55));
[~, d_idx_28] = min(abs(depths_28 - 55));
[~, d_idx_29] = min(abs(depths_29 - 55));


fig2 = figure('Position',[100 100 1600 600]);

subplot(2,4,1)
hold on
plot([time_21{:}], rampe_rl_depth_21(d_idx_21,:), 'b.')
plot([time_21{:}], rl_depth_avg_21, 'k.')
xlim([datetime(2020, 1, 21, 0, 0, 0) datetime(2020, 1, 22, 0, 0, 0)])
ylim([110 150])

hold off

subplot(2,4,2)
hold on
plot([time_24{:}], rampe_rl_depth_24(d_idx_24,:), 'b.')
plot([time_24{:}], rl_depth_avg_24, 'k.')
xlim([datetime(2020, 1, 24, 0, 0, 0) datetime(2020, 1, 25, 0, 0, 0)])
ylim([110 150])
hold off

subplot(2,4,3)
hold on
plot([time_28{:}], rampe_rl_depth_28(d_idx_28,:), 'b.')
plot([time_28{:}], rl_depth_avg_28, 'k.')
xlim([datetime(2020, 1, 28, 0, 0, 0) datetime(2020, 1, 29, 0, 0, 0)])
ylim([110 150])
hold off

subplot(2,4,4)
hold on
plot([time_29{:}], rampe_rl_depth_29(d_idx_29,:), 'b.', 'DisplayName', '55m Channel')
plot([time_29{:}], rl_depth_avg_29, 'k.', 'DisplayName', 'Array Average')
xlim([datetime(2020, 1, 29, 0, 0, 0) datetime(2020, 1, 30, 0, 0, 0)])
ylim([110 150])
legend('Location','northwest','FontSize', 10)
hold off

subplot(2,4,5)
hold on
plot([time_21{:}], rampe_rl_depth_21_tidal(d_idx_21,:), 'b.')
plot([time_21{:}], rl_depth_avg_21_tidal, 'k.')
xlim([datetime(2020, 1, 21, 0, 0, 0) datetime(2020, 1, 22, 0, 0, 0)])
ylim([110 150])
hold off

subplot(2,4,6)
hold on
plot([time_24{:}], rampe_rl_depth_24_tidal(d_idx_24,:), 'b.')
plot([time_24{:}], rl_depth_avg_24_tidal, 'k.')
xlim([datetime(2020, 1, 24, 0, 0, 0) datetime(2020, 1, 25, 0, 0, 0)])
ylim([110 150])
hold off

subplot(2,4,7)
hold on
plot([time_28{:}], rampe_rl_depth_28_tidal(d_idx_28,:), 'b.')
plot([time_28{:}], rl_depth_avg_28_tidal, 'k.')
xlim([datetime(2020, 1, 28, 0, 0, 0) datetime(2020, 1, 29, 0, 0, 0)])
ylim([110 150])
hold off

subplot(2,4,8)
hold on
plot([time_29{:}], rampe_rl_depth_29_tidal(d_idx_29,:), 'b.')
plot([time_29{:}], rl_depth_avg_29_tidal, 'k.')
xlim([datetime(2020, 1, 29, 0, 0, 0) datetime(2020, 1, 30, 0, 0, 0)])
ylim([110 150])
hold off

han = axes(fig2,'visible','off');
han.Title.Visible = 'on';
han.XLabel.Visible = 'on';
han.YLabel.Visible = 'on';

xlabel(han,'Time');
ylabel(han,'Receiver Level (dB re 1 µPa)');
title(han,'RAMPE Receiver Level');
exportgraphics(fig2, fullfile(save_dir, "rl_depth_st_and_tidal.png"), "Resolution", 300)