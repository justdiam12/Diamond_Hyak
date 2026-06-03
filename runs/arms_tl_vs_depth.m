%%% Directories %%%
hyak_dir = '/mmfs1/gscratch/stf/jdiam12/Diamond_Hyak';
rampe_save_dir = fullfile(hyak_dir, 'save', 'RAMPE');
bellhop_save_dir = fullfile(hyak_dir, 'save', 'BELLHOP');
save_dir = fullfile(hyak_dir, 'save', 'plots');
depths_dir = fullfile(hyak_dir, 'data', 'array_depths');
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


%%% Model Parameters %%%
freqs = 3450:3550;
r_range = 5500;
source_level = 195;
time = 0:0.0001:2;


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


%%% Sift through BELLHOP Save Files %%%
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


%%% TL at Receiver Depth over time %%%
rampe_tl_depth_21 = zeros(length(depths_21), length(rampe_21_files)-1);
rampe_tl_depth_24 = zeros(length(depths_24), length(rampe_24_files)-1);
rampe_tl_depth_28 = zeros(length(depths_28), length(rampe_28_files)-1);
rampe_tl_depth_29 = zeros(length(depths_29), length(rampe_29_files)-1);
bellhop_tl_depth_21 = zeros(length(depths_21), length(bellhop_21_files)-1);
bellhop_tl_depth_24 = zeros(length(depths_24), length(bellhop_24_files)-1);
bellhop_tl_depth_28 = zeros(length(depths_28), length(bellhop_28_files)-1);
bellhop_tl_depth_29 = zeros(length(depths_29), length(bellhop_29_files)-1);

c = parcluster('local');

jobid = getenv('SLURM_JOB_ID');
if isempty(jobid)
    jobid = char(java.util.UUID.randomUUID);
end

c.JobStorageLocation = fullfile('/mmfs1/gscratch/stf/jdiam12', ...
    'matlab_job_storage', jobid);

if ~exist(c.JobStorageLocation,'dir')
    mkdir(c.JobStorageLocation)
end

saveProfile(c);

nWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));
if isnan(nWorkers)
    nWorkers = feature('numcores');
end

pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);
end

parpool(c, nWorkers);


%%% Loop through RAMPE Save Files %%%
% January 21, 2020
parfor s = 1:length(rampe_21_files)

    type = 1;
    disp("Working on file " + rampe_21_files(s))

    data = load(fullfile(rampe_save_dir, rampe_21_files{s}));
    Depth = data.Depth - 150;
    Depth_mask = (Depth >= -20 & Depth <= 200);
    Depth = Depth(Depth_mask);
    Psv = data.Psv(Depth_mask,:,:);
    Range = data.Range;

    tl = tl_vs_depth( ...
        source_level, Psv, Depth, Range, ...
        freqs, time, depths_21, r_range, ...
        source_signal, source_time, type);
    rampe_tl_depth_21(:,s) = tl;

end

% January 24, 2020
parfor s = 1:length(rampe_24_files)

    type = 1;
    disp("Working on file " + rampe_24_files(s))

    data = load(fullfile(rampe_save_dir, rampe_24_files{s}));
    Depth = data.Depth - 150;
    Depth_mask = (Depth >= -20 & Depth <= 200);
    Depth = Depth(Depth_mask);
    Psv = data.Psv(Depth_mask,:,:);
    Range = data.Range;

    tl = tl_vs_depth( ...
        source_level, Psv, Depth, Range, ...
        freqs, time, depths_24, r_range, ...
        source_signal, source_time, type);
    rampe_tl_depth_24(:,s) = tl;

end

% January 28, 2020
parfor s = 1:length(rampe_28_files)

    type = 1;
    disp("Working on file " + rampe_28_files(s))

    data = load(fullfile(rampe_save_dir, rampe_28_files{s}));
    Depth = data.Depth - 150;
    Depth_mask = (Depth >= -20 & Depth <= 200);
    Depth = Depth(Depth_mask);
    Psv = data.Psv(Depth_mask,:,:);
    Range = data.Range;

    tl = tl_vs_depth( ...
        source_level, Psv, Depth, Range, ...
        freqs, time, depths_28, r_range, ...
        source_signal, source_time, type);
    rampe_tl_depth_28(:,s) = tl;

end

% January 29, 2020
parfor s = 1:length(rampe_29_files)

    type = 1;
    disp("Working on file " + rampe_29_files(s))

    data = load(fullfile(rampe_save_dir, rampe_29_files{s}));
    Depth = data.Depth - 150;
    Depth_mask = (Depth >= -20 & Depth <= 200);
    Depth = Depth(Depth_mask);
    Psv = data.Psv(Depth_mask,:,:);
    Range = data.Range;

    tl = tl_vs_depth( ...
        source_level, Psv, Depth, Range, ...
        freqs, time, depths_29, r_range, ...
        source_signal, source_time, type);
    rampe_tl_depth_29(:,s) = tl;

end

%%% Mean and STDDEV of TL
rampe_avg_21 = zeros(length(depths_21), 1);
rampe_avg_24 = zeros(length(depths_24), 1);
rampe_avg_28 = zeros(length(depths_28), 1);
rampe_avg_29 = zeros(length(depths_29), 1);
rampe_std_21 = zeros(length(depths_21), 1);
rampe_std_24 = zeros(length(depths_24), 1);
rampe_std_28 = zeros(length(depths_28), 1);
rampe_std_29 = zeros(length(depths_29), 1);

for d = 1:length(depths_21)
    rampe_avg_21(d) = mean(rampe_tl_depth_21(d,:));
    rampe_avg_24(d) = mean(rampe_tl_depth_24(d,:));
    rampe_avg_28(d) = mean(rampe_tl_depth_28(d,:));
    rampe_avg_29(d) = mean(rampe_tl_depth_29(d,:));
    rampe_std_21(d) = std(rampe_tl_depth_21(d,:));
    rampe_std_24(d) = std(rampe_tl_depth_24(d,:));
    rampe_std_28(d) = std(rampe_tl_depth_28(d,:));
    rampe_std_29(d) = std(rampe_tl_depth_29(d,:));
end

%%% Plot RAMPE Figure %%%
figure('Position',[100 100 800 1400])
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
set(gca, 'FontSize', 14)
title('Daily-Averaged TL vs Depth (RAMPE)', fontsize=20)
xlabel('Average TL (dB re 1 µPa)', fontsize=16)
ylabel('Depth (m)', fontsize=16)
xlim([60 80])
legend('Location','northeast','FontSize', 12)
grid on
hold off
exportgraphics(gcf, fullfile(save_dir, "rampe_tl_depth.png"), "Resolution", 300)

disp("RAMPE Plot Complete.")

%%% Loop through BELLHOP Save Files
% January 21, 2020
parfor s = 1:length(bellhop_21_files)

    type = 1;
    disp("Working on file " + bellhop_21_files(s))

    data = load(fullfile(bellhop_save_dir, bellhop_21_files{s}));
    Depth = data.Depth
    Psv = data.Psv;
    Range = data.Range;

    tl = tl_vs_depth( ...
        source_level, Psv, Depth, Range, ...
        freqs, time, depths_21, r_range, ...
        source_signal, source_time, type);
    bellhop_tl_depth_21(:,s) = tl;

end

% January 24, 2020
parfor s = 1:length(bellhop_24_files)

    type = 1;
    disp("Working on file " + bellhop_24_files(s))

    data = load(fullfile(bellhop_save_dir, bellhop_24_files{s}));
    Depth = data.Depth
    Psv = data.Psv;
    Range = data.Range;

    tl = tl_vs_depth( ...
        source_level, Psv, Depth, Range, ...
        freqs, time, depths_24, r_range, ...
        source_signal, source_time, type);
    bellhop_tl_depth_24(:,s) = tl;

end

% January 28, 2020
parfor s = 1:length(bellhop_28_files)

    type = 1;
    disp("Working on file " + bellhop_28_files(s))

    data = load(fullfile(bellhop_save_dir, bellhop_28_files{s}));
    Depth = data.Depth
    Psv = data.Psv;
    Range = data.Range;

    tl = tl_vs_depth( ...
        source_level, Psv, Depth, Range, ...
        freqs, time, depths_28, r_range, ...
        source_signal, source_time, type);
    bellhop_tl_depth_28(:,s) = tl;

end

% January 29, 2020
parfor s = 1:length(bellhop_29_files)

    type = 1;
    disp("Working on file " + bellhop_29_files(s))

    data = load(fullfile(bellhop_save_dir, bellhop_29_files{s}));
    Depth = data.Depth
    Psv = data.Psv;
    Range = data.Range;

    tl = tl_vs_depth( ...
        source_level, Psv, Depth, Range, ...
        freqs, time, depths_29, r_range, ...
        source_signal, source_time, type);
    bellhop_tl_depth_29(:,s) = tl;

end


%%% Mean and STDDEV of TL
bellhop_avg_21 = zeros(length(depths_21), 1);
bellhop_avg_24 = zeros(length(depths_24), 1);
bellhop_avg_28 = zeros(length(depths_28), 1);
bellhop_avg_29 = zeros(length(depths_29), 1);
bellhop_std_21 = zeros(length(depths_21), 1);
bellhop_std_24 = zeros(length(depths_24), 1);
bellhop_std_28 = zeros(length(depths_28), 1);
bellhop_std_29 = zeros(length(depths_29), 1);


for d = 1:length(depths_21)
    bellhop_avg_21(d) = mean(bellhop_tl_depth_21(d,:));
    bellhop_avg_24(d) = mean(bellhop_tl_depth_24(d,:));
    bellhop_avg_28(d) = mean(bellhop_tl_depth_28(d,:));
    bellhop_avg_29(d) = mean(bellhop_tl_depth_29(d,:));
    bellhop_std_21(d) = std(bellhop_tl_depth_21(d,:));
    bellhop_std_24(d) = std(bellhop_tl_depth_24(d,:));
    bellhop_std_28(d) = std(bellhop_tl_depth_28(d,:));
    bellhop_std_29(d) = std(bellhop_tl_depth_29(d,:));
end


%%% Plot BELLHOP %%%
figure('Position',[100 100 800 1400])
hold on

errorbar(bellhop_avg_21, depths_21, bellhop_std_21, 'horizontal', ...
    'bo', 'DisplayName', 'January 21', 'LineWidth', 2)
errorbar(bellhop_avg_24, depths_24, bellhop_std_24, 'horizontal', ...
    'ro', 'DisplayName', 'January 24', 'LineWidth', 2)
errorbar(bellhop_avg_28, depths_28, bellhop_std_28, 'horizontal', ...
    'go', 'DisplayName', 'January 28', 'LineWidth', 2)
errorbar(bellhop_avg_29, depths_29, bellhop_std_29, 'horizontal', ...
    'mo', 'DisplayName', 'January 29', 'LineWidth', 2)

set(gca, 'YDir', 'reverse')
set(gca, 'FontSize', 14)
title('Daily-Averaged TL vs Depth (BELLHOP)', fontsize=20)
xlabel('Average TL (dB re 1 µPa)', fontsize=16)
ylabel('Depth (m)', fontsize=16)
xlim([60 80])
legend('Location','northeast','FontSize', 12)
grid on
hold off
exportgraphics(gcf, fullfile(save_dir, "bellhop_tl_depth.png"), "Resolution", 300)

disp("BELLHOP Plot Complete.")