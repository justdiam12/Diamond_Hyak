% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program takes the following variables as input:
%           arrivals: Array of Arrival Amplitudes and Times
%           freqs: Frequency Array (Hz)
%           rd: Receiver Depths (m)
%           s_sig: Source Signal (unitless)
%           s_t: Source Time (s)
%
% The outputs of this program are as follows:
%           r_sig: Receiver Signal Array (µPa), 
%           r_t: Receiver Time Array (s)
% 
% Written by: Justin Diamond
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tl] = tl_vs_depth_bh(arrivals,freqs,time, rd,s_sig,s_t)
    
    hold on
    % Loop through depth to extract impusle responses at receiver depths
    for d = 1:length(rd)
        
        arrivals = arrivals;
        
        % Extract Amplitude at Arrival Depth
        H = zeros(length(freqs), 1);
        for f = 1:length(freqs)

            H(f) = 0;

            % Arrival index
            arr = arrivals{f,d};
            for n = 1:length(arr.amp)
                arr_time = arr.time(n);
                phase = arr.phase(n);
                amp = arr.amp(n) * exp(1j * deg2rad(phase));

                H(f) = H(f) + amp ...
                    .* exp(-1i*2*pi*freqs(f)*arr_time);
            end
        
        end

        % Loop through time to extract impulse response
        for t = 1:length(time)
            sig(t) = trapz(freqs, transpose(H) ...
                     .*exp(1i * 2 * pi * freqs * time(t)));
        end
        
        [~, max_idx] = max(abs(real(sig)));

        lower_t = time(max_idx) - 0.05;
        upper_t = time(max_idx) + 0.15;
         
        % plot(time, sig)
        % xline(time(max_idx))

        % Window Source Signal
        t_idx = (time >= lower_t & time <= upper_t);
        time_w = time(t_idx) - lower_t;
        sig_w = sig(t_idx);
        dt_sig = time(2) - time(1);

        % Source Signal
        dt_source = s_t(2) - s_t(1);
        dt = min(dt_sig, dt_source);
        t_start = min([time_w(1), s_t(1)]);
        t_end   = max([time_w(end), s_t(end)]);
        
        % Common Time Grid
        t_common = t_start:dt:t_end;
        sig_i = interp1(time_w, sig_w, t_common, 'linear', 0);
        src_i = interp1(s_t, s_sig, t_common, 'linear', 0);

        % Save Received Signal at Specific Depth
        rec_sig = conv(sig_i, src_i, 'full') * dt;
        r_t = 0:dt:(length(rec_sig)-1)*dt;
        plot(r_t, abs(real(rec_sig)))
        r_idx = (r_t >= 0.1 & r_t <= 1.1);
        r_t = r_t(r_idx);
        rec_sig = rec_sig(r_idx);
        tl(d) = abs(20*log10(rms(rec_sig))); % Units (µPa)
    end
    hold off
end