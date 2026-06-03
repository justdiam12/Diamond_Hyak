%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Animate Transmission Loss vs Depth Over Time
%
% Assumes saved MAT file contains:
%   tl_depth   -> [Ndepth x Ntime]
%   depths     -> receiver depths
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
clc

%%%%%%%%%%%%%%%%%%%%
%%% Load Results %%%
%%%%%%%%%%%%%%%%%%%%
load("rampe_tl.mat")   % <-- change filename

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Optional Time Labels / Datetime %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example:
% time_labels = datetime(2020,1,28,0,0,0) + hours(0:size(tl_depth,2)-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create Figure Window %
%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Color','w', Position=[100 100 400 1200])

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Animation Loop %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
for t = 1:size(tl_depth,2)

    clf

    % Plot TL profile
    plot(tl_depth(:,t), depths, ...
         'LineWidth', 2)

    % Flip depth axis
    set(gca,'YDir','reverse')

    % Labels
    xlabel('Transmission Loss (dB)')
    ylabel('Depth (m)')

    % Axis limits
    xlim([min(tl_depth(:)) max(tl_depth(:))])
    ylim([min(depths) max(depths)])

    grid on

    % Title
    title(sprintf('Transmission Loss vs Depth\nTime Index %d', t))

    % If using datetime labels:
    % title(sprintf('Transmission Loss vs Depth\n%s', ...
    %     datestr(time_labels(t))))

    drawnow

    pause(0.1)

end

v = VideoWriter('tl_depth_animation.mp4','MPEG-4');
v.FrameRate = 2;
open(v);

figure('Color','w', Position=[100  100 400 1200])

for t = 1:size(tl_depth,2)

    clf

    plot(tl_depth(:,t), depths, 'LineWidth',2)

    set(gca,'YDir','reverse')

    xlabel('Transmission Loss (dB)')
    ylabel('Depth (m)')

    xlim([min(tl_depth(:)) max(tl_depth(:))])
    ylim([min(depths) max(depths)])

    grid on

    title(sprintf('Transmission Loss vs Depth\nTime Index %d', t))

    frame = getframe(gcf);
    writeVideo(v, frame);

end

close(v)