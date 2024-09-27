function [ecg_analysis] = analyse_ecg(ecg_data,leadn,show_figure, name)
%ANALYSE_ECG Built on 'get_R_addrlist_saa.m', this function performs
%all beat detection and annotation of holter signals

nointerp_index = [];
%% Define ecg_analysis structure

ecg_analysis.Beatstart = [];    % Locations of beats in recording, provided by annotation file
ecg_analysis.Beattimes = [];    % Timestamps of beats in recording, based on holter start time
ecg_analysis.Rloc = [];         % Locations of R peaks
ecg_analysis.Rpeakval = [];     % The yvalue of R peak.
ecg_analysis.Ramp = [];         % Amplitude of R peaks, corresponding to Rloc
ecg_analysis.RR = [];           % Size of RR intervals, corresponding to Rloc
ecg_analysis.RRint = [];        % The RR intervals in milliseconds
ecg_analysis.Qloc = [];         % Locations of Q wave
ecg_analysis.Tloc = [];         % Locations of T peaks,
ecg_analysis.Tamp = [];         % Amplitude of T peaks, corresponding to Tloc
ecg_analysis.Tend = [];         % End of t-wave, for determining QT interval
ecg_analysis.QTpeak = [];           % QT interval, calculated difference between Q location and T end
ecg_analysis.Baseline = [];     % Baseline beat, calculated through an average of samples between Tloc and the next Rloc
ecg_analysis.Potentialevents = [];% Potential events, flagged by QC of beat intervals, an array corresponding to the length of Beatstart (1 = EAD, 2 = DAD)
ecg_analysis.Rampevents = [];   % (1 = Lower than average, 2 = Higher than average, 3 = Missing R wave)
ecg_analysis.Tampevents = [];   % (1 = Lower than average, 2 = Higher than average)
ecg_analysis.Potentialectopics = []; % Potential ectopic beats, corresponding to the beats where the detected T wave is greater than or close to the amplitude of the R wave
ecg_analysis.TPars = []; % To be filled with all data regarding T30, T50 and T70 locations (x and y).
ecg_analysis.QC = []; % TO be filled with 1 (accepted beats) or 0 (rejected beats)


sample_rate = ecg_data.arash_Header.Sampling_Rate;

%% Check for calibration pulses


%% Populate ecg_analysis structure and remove data points from beats other 
% than normal ('N')
% In some cases, data is annotated beyond actual recording length. To
% remove from analysis, mark on annotation array:

raw_ecg = ecg_data.raw_ecgSig(:,leadn);
% filtered_ecg = bandpass(raw_ecg,[0.5 30],sample_rate,'Steepness',0.85,'StopbandAttenuation',60);

%Butterworth:
order = 2;
Fc = [0.5 90];
Fs = 200;
[n,h] = butter(order,Fc/(Fs/2));
filtered_ecg = filtfilt(n,h,raw_ecg);
ecg_analysis.filtered_ecg = filtered_ecg;

% histograms
% % figure();
% % hold on
% % plot(raw_ecg);
% % plot(filtered_ecg);
% % % input('Please press ENTER to continue...');
% % close all


% smooth_ecg = smooth(raw_ecg,ceil(sample_rate*0.1));

if ~isempty(ecg_data.Ann)
    annotations = char(ecg_data.Ann);
    annotations(ecg_data.Rloc>length(raw_ecg)) = 'O'; % 'O' denoting locations [O]utside recording
    mean_rr_int = mean(ecg_data.RR);
    ecg_analysis.RR = ecg_data.RR(annotations=='N');
    ecg_analysis.Beatstart = ecg_data.Rloc(annotations=='N'); 
    ecg_analysis.Potentialevents = zeros(1,numel(ecg_analysis.Beatstart));
else
    ecg_analysis.Beatstart = find_r_loc(filtered_ecg);
    ecg_analysis.RR = diff(ecg_analysis.Beatstart);
    ecg_analysis.Potentialevents = zeros(1,numel(ecg_analysis.Beatstart));
end

%% Beatstart/RR quality control
% Flagging beat intervals that wildly differ from the average of
% the last 5 beat intervals, using coefficient of variance at a threshold
% of 0.2 (stdev/mean)

beatdiffs = [];
CoV_threshold = 0.2;
previouscounts = 10;
event_threshold = .5;

% for index = 1:numel(ecg_analysis.Beatstart)-1
%     if numel(beatdiffs) < previouscounts           % For first count indexes and whilst the cov of these beats are > than threshold, find
%         beatdiffs = [beatdiffs, ecg_analysis.Beatstart(index+1)-ecg_analysis.Beatstart(index)]; 
%     elseif numel(beatdiffs) == previouscounts
%         tempdiffs = [beatdiffs, ecg_analysis.Beatstart(index+1)-ecg_analysis.Beatstart(index)];
%         CoV = std(tempdiffs)/mean(tempdiffs);
%         if CoV > CoV_threshold
%             beatdiffs = [beatdiffs(2:end), ecg_analysis.Beatstart(index+1)-ecg_analysis.Beatstart(index)];
%         else
%             beatdiffs = tempdiffs;
%         end
%     else
%         tempdiffs = [beatdiffs(2:end),ecg_analysis.Beatstart(index+1)-ecg_analysis.Beatstart(index)];
%         diffratio = tempdiffs(end)/mean(beatdiffs);
%         if diffratio <= 1-event_threshold
%             ecg_analysis.Potentialevents(index) = 1; % Short event, e.g. EAD
%         elseif  diffratio >= 1+event_threshold
%             ecg_analysis.Potentialevents(index) = 2; % Long event, e.g. DAD
%         else
%             beatdiffs = tempdiffs;
%         end
%     end     
% end

QRS_estimate = floor(mean(ecg_analysis.RR)/5); % General estimate of QRS complex length


%% Find and populate R peak locations and amplitudes

% r_peak_list = zeros(size(ecg_analysis.Beatstart));
% 
% for i = 1 : length(ecg_analysis.Beatstart)
%     [maxval index] = max(filtered_ecg(ecg_analysis.Beatstart(i)-QRS_estimate:min(ecg_analysis.Beatstart(i)+QRS_estimate)));
%     r_peak_list(i) = index+ecg_analysis.Beatstart(i)-QRS_estimate-1;
% end
% 
% end_r_estimate = floor( 3.0 * mean(r_peak_list-ecg_analysis.Beatstart));
% if  (end_r_estimate <= 6) %???
%     end_r_estimate = QRS_estimate;
% end

r_peak_list = zeros(size(ecg_analysis.Beatstart));

for i = 1 : length(ecg_analysis.Beatstart)
    [maxval, index] = max(filtered_ecg(max(ecg_analysis.Beatstart(i)-QRS_estimate,1):min(ecg_analysis.Beatstart(i)+QRS_estimate,numel(filtered_ecg)))); % Checking if search end does not exceed length of ecg
    r_peak_list(i) = index+ecg_analysis.Beatstart(i)-QRS_estimate-1;
end

end_r_estimate = floor( 3.0 * mean(r_peak_list-ecg_analysis.Beatstart));
if  (end_r_estimate <= 6) %???
    end_r_estimate = QRS_estimate;
end


%% Iterate through r locations and claculate peaks/periods
r_peak_list = zeros(1,size(ecg_analysis.Beatstart,2)-1);
r_min_list = r_peak_list;
q_list = r_peak_list;
r_amplitude = r_peak_list;
r_amp_events = r_peak_list;
t_peak_list = r_peak_list;
t_end_list = r_peak_list;
t_amplitude = r_peak_list;
t_amp_events = r_peak_list;
ectopic_events = r_peak_list;
baseline = r_peak_list;
Rpeakval = r_peak_list;

QC = ones(1,size(ecg_analysis.Beatstart,2)-1);

for i = 1 : length(ecg_analysis.Beatstart)-1

    % snip = filtered_ecg(ecg_analysis.Beatstart(i):ecg_analysis.Beatstart(i+1));
    % baselinedsnip = msbackadj((1:numel(snip))',snip,'WindowSize',numel(snip),'StepSize',numel(snip));
    % plot(baselinedsnip)
    % drawnow
    % pause(0.001);

    %% Finding R peak
    search_samples = filtered_ecg(max(ecg_analysis.Beatstart(i)-QRS_estimate,1):min(ecg_analysis.Beatstart(i)+QRS_estimate,numel(filtered_ecg))); % Adding check if search starts from beginning of the ecg
    [Rmaxval, index] = findpeaks(search_samples, 'MinPeakProminence', 10);
    % Findpeaks can pick up p-wave, but t-peak should be greater, using index trick maxval == max(maxval)
    if isempty(Rmaxval) % Incase r wave is not found/not prominent
        r_amp_events(i) = 3;
        r_peak_list(i) = ecg_analysis.Beatstart(i);
    else
        r_peak_list(i) = max(index(Rmaxval == max(Rmaxval))+ecg_analysis.Beatstart(i)-QRS_estimate-1,index(Rmaxval == max(Rmaxval))); % Adding check if search starts from beginning of the ecg
    end
    Rpeakval(i) = filtered_ecg(r_peak_list(i)); % Hold onto r-peak value for later

    %% Define a search period using upstroke and downstroke of QRS complex, found using signal derivative max and min
    [~, maxIndex] = max(diff(search_samples));
    [~, minIndex] = min(diff(search_samples(maxIndex:end)));
    search_window = minIndex;

    
    %% Finding R minimum
    [minval, index] = min(filtered_ecg(r_peak_list(i):min(r_peak_list(i)+2*search_window)));
    r_min_list(i) = index+r_peak_list(i)-1;


    
    %% Q time
    q_search_start = max(2, r_peak_list(i) - 2*search_window);

    % for j = q_search_start:r_peak_list(i)
    %     if (j < length(raw_ecg)-3 )
    %         delta_1 = raw_ecg(j) - raw_ecg(j+1);
    %         delta_2 = raw_ecg(j+1) - raw_ecg(j+2);
    % 
    %         if ((delta_1 * delta_2) > 0 ) && (abs(delta_1) > (1000/sample_rate)*0.025 * R_amplitude(i)) %200
    %             q_list(i)  = j;
    %             break;
    %         end
    %     else
    %         q_list(i)  = j;
    %         break;
    %     end
    % end

    % Trying different second order differential method of q detection
    % Find first minimum of QRS complex, i.e. when second differential is maximum
    % second_diff = diff(raw_ecg(q_search_start:r_peak_list(i)),2);
    % [maxval, index] = max(second_diff);
    % if ~isempty(index)
    %     q_list(i) = index + q_search_start -1;
    % else
    %     q_list(i) = q_search_start;
    % end

    % Trying different minimum method.
    [~, index] = min(filtered_ecg(q_search_start:r_peak_list(i)));
    q_list(i) = index + q_search_start - 1;



    %% T peak
    % To avoid detecting p wave, the end of the search should be 3/4
    % between r_min and next beat start
    end_t_search = r_min_list(i) + floor((ecg_analysis.Beatstart(i+1)-r_min_list(i))*.7); % TODO: Find more robust detection, will biphasic t_waves be misrepresented? what about negative t-waves???
    t_snip = abs(filtered_ecg(r_min_list(i):end_t_search) - filtered_ecg(r_min_list(i))); % To detect peaks of inverse/biphasic t-waves, set first point of window as baseline and absolute
    [maxval, index] = max(t_snip); % Might not work with lqt2 with camel hump twave
    if ~isempty(index)
        if index == length(t_snip) % Check if found maximum is at the end of the window
            t_peak_list(i) = NaN;

        else
            t_peak_list(i) = index + r_min_list(i) - 1; % minus 1 for index correction.
        end
    else
        t_peak_list(i) = NaN; % Case where beat interval is too small to detect t_peak
    end

    % p = polyfit(r_min_list(i):ecg_analysis.Beatstart(i+1),raw_ecg(r_min_list(i):ecg_analysis.Beatstart(i+1)),5);
    % fitted_curve = polyval(p,r_min_list(i):ecg_analysis.Beatstart(i+1));


    %% Finding baseline of beat
    if ~isnan(t_peak_list(i))
        midpoint_period = ecg_analysis.Beatstart(i+1) - t_peak_list(i);
        midpoint_start = t_peak_list(i) + floor(midpoint_period/4);
        midpoint_end = t_peak_list(i) + floor(midpoint_period*3/4);
        baseline(i) = floor(mean(filtered_ecg(midpoint_start:midpoint_end)));
    else
        baseline(i) = 0;
    end
    
    % % % checks for baseline deviating too much. Problem with it getting
    % % % 'locked' onto a certain value e.g. if the last 5 values are all x,
    % % % then the next value will definitely be more than 2sd away hence will
    % % % become x, and then it gets locked onto x where every iteration will
    % % % be replaced by x.
    % % if i>5
    % %     up_threshold = mean(baseline(i-5:i-1)) + (2*std(baseline(i-5:i-1)));
    % %     down_threshold = mean(baseline(i-5:i-1)) - (2*std(baseline(i-5:i-1)));
    % %     if  (baseline(i) > up_threshold | baseline(i) < down_threshold)
    % %         baseline(i) = round(mean(baseline(i-5:i-1)));
    % %     end
    % % end
    % % % if baseline(i) %<> xSD avg of last 5  baselines,replace with the rounded baseline avg.

    %% T end - calculating using intersection of dv/dt tangent to horizontal from baseline: Temporarily commented out - unreliable in noisy signal
    
    % snip = filtered_ecg(r_peak_list(i):ecg_analysis.Beatstart(i+1)); % Use r_peak as reference
    % search_stop = floor((ecg_analysis.Beatstart(i+1) - t_peak_list(i))*.5); %stop searching just before the p-wave/qrs complex
    % [gradient, index] = min(diff(filtered_ecg(t_peak_list(i):t_peak_list(i)+search_stop)));
    % gradient_point = t_peak_list(i) - r_peak_list(i) + index;
    % if ~isscalar(gradient_point)
    %     input('ERROR HERE!!!!')
    % end
    % tangent_x = linspace(gradient_point - 24, gradient_point + 25,50); % Search for 50 samples around gradient point
    % tangent_y = gradient*(tangent_x-gradient_point) + snip(gradient_point); % Re-arranged point gradient formula
    % baseline_x = 1:numel(snip);
    % baseline_y = zeros(size(baseline_x));
    % baseline_y(:) = baseline(i);
    % 
    % % V Was used for horizontal baseline of signal before baseline
    % % approximation
    % % q_line_x = 1:numel(snip);
    % % q_line_y = zeros(size(q_line_x));
    % % q_line_y(:) = filtered_ecg(q_list(i));
    % 
    % [~, intersections, ~] = find(tangent_y <= baseline(i)); % Find tangent intersection with q-line
    % 
    % if ~isempty(intersections) % Check that an intersection has actually been found!!!!
    %     t_end_list(i) = ecg_analysis.Beatstart(i) + tangent_x(intersections(1)) ;
    % else
    %     [minval, index] = min(diff(filtered_ecg(t_peak_list(i):ecg_analysis.Beatstart(i+1)),3));
    %     t_end_list(i) = index + t_peak_list(i);
    % end

    % plot(snip)
    % hold on
    % plot(diff(snip))
    % plot(tangent_x, tangent_y)
    % plot(baseline_x,baseline_y)
    % plot(gradient_point,snip(gradient_point),'r*')
    % if tangent_x(intersections(1)) < numel(snip)
    %     plot(tangent_x(intersections(1)),snip(tangent_x(intersections(1))),'b*')
    % end
    % drawnow
    % pause(0.001);
    % hold off

    % [minval, index] = min(diff(filtered_ecg(t_peak_list(i):ecg_analysis.Beatstart(i+1)),3));
    % t_end_list(i) = index + t_peak_list(i);


    %% Find R amplitude
    r_amplitude(i) = Rpeakval(i) - baseline(i);

    amplitude_threshold = .5;

    % Check if R_amplitude differ from previous beats
    if i > previouscounts
        rampratio = r_amplitude(i)/mean(r_amplitude(i-previouscounts+1:i-1));
        if rampratio < 1 - amplitude_threshold
            r_amp_events(i) = 1;
        elseif rampratio > 1 + amplitude_threshold
            r_amp_events(i) = 2;
        end
    end
    
    %if (R_amplitude(i) >= 400)             %Commented out
    %    R_amplitude(i) =0;
    %end
  

    %% T amplitude
    if ~isnan(t_peak_list(i))
        t_amplitude(i) =  filtered_ecg(t_peak_list(i)) - baseline(i); % - filtered_ecg(q_list(i));  % Amplitude currently measured from peak to 't-end' (down slope)
        if i > previouscounts
            tampratio = t_amplitude(i)/mean(t_amplitude(i-previouscounts+1:i-1));
            if tampratio < 1 - amplitude_threshold
                t_amp_events(i) = 1;
            elseif tampratio > 1 + amplitude_threshold
                t_amp_events(i) = 2;
            end
        end
  
    else
        t_amplitude(i) = NaN;
        t_amp_events(i) = 3;
    end

    % if t_amplitude(i) <= 0
    %     t_amplitude(i) = NaN;
    % end

    %% Pre-Tpeak and post-Tpeak split. 
    % t30 = Location at 30% from Tpeak height. t50 = Location at 50% of Tpeak height. t70 = Location at 70% from Tpeak height.

    if ~isnan(t_peak_list(i)) | ~isnan(t_amplitude(i))
        end_t_snip = r_min_list(i) + floor((ecg_analysis.Beatstart(i+1)-r_min_list(i))*.7);
        search_snip = r_min_list(i):end_t_snip;

        % IF we have a negative T wave, then prepare the values:
        if filtered_ecg(t_peak_list(i)) < baseline(i) % if ecg_analysis(t_peak_list(i)) < baseline(i), then its a negative t deflection

            t_split.t30_height(i) = baseline(i) - (abs(t_amplitude(i)) .* (1-0.3)); % negative t waves need Baseline to be SUBTRACTION i.e. baseline - (0.7*amp)
            first_above_t30 = find(filtered_ecg(search_snip) <=  t_split.t30_height(i), 1, 'first'); % Negative t waves need to detect first and last values that are less than t?_height
            last_above_t30 = find(filtered_ecg(search_snip) <= t_split.t30_height(i), 1, 'last');

            t_split.t50_height(i) = baseline(i) - (abs(t_amplitude(i)) .* 0.5);
            first_above_t50 = find(filtered_ecg(search_snip) <= t_split.t50_height(i), 1, 'first');
            last_above_t50 = find(filtered_ecg(search_snip) <= t_split.t50_height(i), 1, 'last');

            t_split.t70_height(i) = baseline(i) - (abs(t_amplitude(i)) .* (1-0.7));
            first_above_t70 = find(filtered_ecg(search_snip) <= t_split.t70_height(i), 1, 'first');
            last_above_t70 = find(filtered_ecg(search_snip) <= t_split.t70_height(i), 1, 'last');

        % ELSEIF we have a positive T Wave, then prepare the values.
        elseif filtered_ecg(t_peak_list(i)) > baseline(i)

            t_split.t30_height(i) = baseline(i) + (t_amplitude(i) .* (1-0.3)); % Positive T waves need baseline to be ADDITION!i.e. baseline + (0.7*amp)
            first_above_t30 = find(filtered_ecg(search_snip) >= t_split.t30_height(i), 1, 'first'); % Positive t waves need to detect first and last values that are MORE than t?_height
            last_above_t30 = find(filtered_ecg(search_snip) >= t_split.t30_height(i), 1, 'last');

            t_split.t50_height(i) = baseline(i) + (t_amplitude(i) .* 0.5);
            first_above_t50 = find(filtered_ecg(search_snip) >= t_split.t50_height(i), 1, 'first');
            last_above_t50 = find(filtered_ecg(search_snip) >= t_split.t50_height(i), 1, 'last');

            t_split.t70_height(i) = baseline(i) + (t_amplitude(i) .* (1-0.7));
            first_above_t70 = find(filtered_ecg(search_snip) >= t_split.t70_height(i), 1, 'first');
            last_above_t70 = find(filtered_ecg(search_snip) >= t_split.t70_height(i), 1, 'last');
        end

        if any([first_above_t30, first_above_t50, first_above_t70] <= 1) | any([last_above_t30, last_above_t50, last_above_t70] == length(search_snip))
            nointerp_index(end+1) = i; % If any values are the first/last inside search_snip start a counter.  % disp('No interpolation performed. first/last_above_t?0 is the first/last value in search_snip.')

            t_split.prepeak_t30_loc(i) = first_above_t30 + r_min_list(i) - 1; % Need to do minus 1 for index correction
            t_split.postpeak_t30_loc(i) = last_above_t30 + r_min_list(i) - 1;
            t_split.prepeak_t50_loc(i) = first_above_t50 + r_min_list(i) - 1;
            t_split.postpeak_t50_loc(i) = last_above_t50 + r_min_list(i) - 1;
            t_split.prepeak_t70_loc(i) = first_above_t70 + r_min_list(i) - 1;
            t_split.postpeak_t70_loc(i) = last_above_t70 + r_min_list(i) - 1;
        else
            % elseif first_above t_30 is > 1 and < length(search_snip), do the interpolation
            % Now the values are prepared, interpolate the x values at the t30, t50 and t70 heights.
            t_split.prepeak_t30_loc(i) = interpolate_for_x(t_split.t30_height(i), first_above_t30-1, filtered_ecg(search_snip(first_above_t30-1)), first_above_t30, filtered_ecg(search_snip(first_above_t30))) + r_min_list(i) - 1; % Interpolate function Then adding r_min_list for correct indexing outside of search_snip
            t_split.postpeak_t30_loc(i) = interpolate_for_x(t_split.t30_height(i), last_above_t30, filtered_ecg(search_snip(last_above_t30)), last_above_t30+1, filtered_ecg(search_snip(last_above_t30+1))) + r_min_list(i) - 1; % Minus 1 for index correction
            t_split.prepeak_t50_loc(i) = interpolate_for_x(t_split.t50_height(i), first_above_t50-1, filtered_ecg(search_snip(first_above_t50-1)), first_above_t50, filtered_ecg(search_snip(first_above_t50))) + r_min_list(i) - 1;
            t_split.postpeak_t50_loc(i) = interpolate_for_x(t_split.t50_height(i), last_above_t50, filtered_ecg(search_snip(last_above_t50)), last_above_t50+1, filtered_ecg(search_snip(last_above_t50+1))) + r_min_list(i) - 1;
            t_split.prepeak_t70_loc(i) = interpolate_for_x(t_split.t70_height(i), first_above_t70-1, filtered_ecg(search_snip(first_above_t70-1)), first_above_t70, filtered_ecg(search_snip(first_above_t70))) + r_min_list(i) - 1;
            t_split.postpeak_t70_loc(i) = interpolate_for_x(t_split.t70_height(i), last_above_t70, filtered_ecg(search_snip(last_above_t70)), last_above_t70+1, filtered_ecg(search_snip(last_above_t70+1))) + r_min_list(i) - 1;

        end

    else
        %If t_peak_list(i) or t_amplitude(i) was NaN:
        t_split.prepeak_t30_loc(i) = NaN;
        t_split.postpeak_t30_loc(i) = NaN;
        t_split.prepeak_t50_loc(i) = NaN;
        t_split.postpeak_t50_loc(i) = NaN;
        t_split.prepeak_t70_loc(i) = NaN;
        t_split.postpeak_t70_loc(i) = NaN;
    end
   
%% Area under curves.

    %% Potential ectopic beats:- determined by comparing R wave and respective T wave, if both are similar or T wave is greater, then it is potentially an ectopic beat
    ectopic_threshold = 0.8; % Threshold for how large T wave is before considering it a potential ectopic
    if t_amplitude(i) >= r_amplitude(i)*ectopic_threshold
        ectopic_events(i) = 1;
    end

end

  disp("Tsplit beats not interpolated count = " +  size(nointerp_index, 2))
%% Plot select annotations on FILTERED ecg signal 

events = ecg_analysis.Potentialevents(1:end-1);

if (show_figure == 1)
    t = get_time(1:length(raw_ecg),sample_rate);
    plot(t,filtered_ecg);
    
    % hold on;plot(get_time(ecg_analysis.Beatstart(events == 1),sample_rate),filtered_ecg(ecg_analysis.Beatstart(events == 1)),'r*');
    % hold on;plot(get_time(ecg_analysis.Beatstart(events == 0),sample_rate),filtered_ecg(ecg_analysis.Beatstart(events == 0)),'b*');
    % hold on;plot(get_time(ecg_analysis.Beatstart(ectopic_events == 0),sample_rate),filtered_ecg(ecg_analysis.Beatstart(ectopic_events == 0)),'r*');
    % hold on;plot(get_time(r_peak_list(events == 0),sample_rate),filtered_ecg(r_peak_list(events == 0)),'g*');
    % hold on;plot(get_time(r_min_list(events == 0),sample_rate),filtered_ecg(r_min_list(events == 0)),'g*');
        % hold on;plot(get_time(t_split.prepeak_t50_loc(~isnan(t_split.prepeak_t50_loc)),sample_rate),t_split.t50_height(~isnan(t_split.prepeak_t50_loc)),'k*')
        % hold on;plot(get_time(t_split.postpeak_t50_loc(~isnan(t_split.postpeak_t50_loc)),sample_rate),t_split.t50_height(~isnan(t_split.postpeak_t50_loc)),'m*')
        % hold on;plot(get_time(t_split.postpeak_t30_loc(~isnan(t_split.postpeak_t30_loc)),sample_rate),t_split.t30_height(~isnan(t_split.postpeak_t30_loc)),'r*')
        % hold on;plot(get_time(t_split.prepeak_t30_loc(~isnan(t_split.prepeak_t30_loc)),sample_rate),t_split.t30_height(~isnan(t_split.prepeak_t30_loc)),'b*')
        % hold on;plot(get_time(t_split.prepeak_t70_loc(~isnan(t_split.prepeak_t70_loc)),sample_rate),t_split.t70_height(~isnan(t_split.prepeak_t70_loc)),'g*')
        % hold on;plot(get_time(t_split.postpeak_t70_loc(~isnan(t_split.postpeak_t70_loc)),sample_rate),t_split.t70_height(~isnan(t_split.postpeak_t70_loc)),'c*')
    % hold on;plot(get_time(t_peak_list(events == 0),sample_rate),filtered_ecg(t_peak_list(events == 0)),'k*');
    % hold on;plot(get_time(ecg_analysis.Beatstart(events == 0),sample_rate),filtered_ecg(ecg_analysis.Beatstart(events == 0)),'r*');
    % hold on;plot(get_time(q_list(events == 0),sample_rate),filtered_ecg(q_list(events == 0)),'b*');
    % input('Press ENTER to continue...');
     close all
end




%% Beat Variability Plot
% plot(t,filtered_ecg)
% hold on
% plot(get_time(r_peak_list(1:end-1),sample_rate),diff(r_peak_list(1:end)*(1000/sample_rate)))   %<--- Convert intervals to time in ms
% title('Rpeak intervals next to Filtered ECG')
% % input('Press ENTER to continue...')
% close all

%% Calculate timestamps of beats

holter_starttime = ecg_data.arash_Header.inf.Start_Time;
hrs = holter_starttime(1);
mins = holter_starttime(2);
secs = holter_starttime(3);
start_milliseconds = hrs*60*60*1000 + mins*60*1000 + secs*1000; %% Convert start time to milliseconds from midnight
midnight = 24*60*60*1000; %% midnight in milliseconds

milliseconds_per_frame = 1000/sample_rate;
beattimes = start_milliseconds +  ecg_analysis.Beatstart*milliseconds_per_frame;
beattimes(beattimes>=midnight) = beattimes(beattimes>=midnight)- midnight;

%% Perform QC on beats before passing into structure

% Baseline QC
QC(zscore(abs(baseline))>3) = 0; % Remove beats with extreme baseline calculations (peaks expected to be wildly miscalculated)

% RR Interval QC
rr_freq = sample_rate ./ diff(r_peak_list);
rr_int =  diff(r_peak_list) .* (1000 / sample_rate);
QC(zscore(abs(rr_int))>3) = 0;

% R Amplitude QC

% T Peak QC
QC(isnan(t_peak_list)) = 0;

% T split QC
QC(isnan(t_split.prepeak_t30_loc)) = 0;
QC(isnan(t_split.postpeak_t30_loc)) = 0;
QC(isnan(t_split.prepeak_t50_loc)) = 0;
QC(isnan(t_split.postpeak_t50_loc)) = 0;
QC(isnan(t_split.prepeak_t70_loc)) = 0;
QC(isnan(t_split.postpeak_t70_loc)) = 0;



% %% Heart rate change plot
% heart_rate_change = heart_rate_change_plot(rr_int(QC(1:end-1)==1));

%% Load analysis into structure:

ecg_analysis.Rloc = r_peak_list(QC==1);    
ecg_analysis.Ramp = r_amplitude(QC==1);
ecg_analysis.RR = rr_freq(QC(1:end-1)==1); % Measured in Hz (1Hz = 60bpm, 2Hz = 120bpm). diff(r_peak_list)/sample_rate; 
ecg_analysis.RRint = rr_int(QC(1:end-1)==1); % RR int in ms. (1 sample every 5ms in 200Hz.)
ecg_analysis.Qloc = q_list(QC==1);    
ecg_analysis.Tloc = t_peak_list(QC==1); 
ecg_analysis.Tend = t_end_list(QC==1);
ecg_analysis.Tamp = t_amplitude(QC==1);
ecg_analysis.QTpeak = (ecg_analysis.Tloc - ecg_analysis.Qloc) .* (1000 / sample_rate); % If Tend could be more robust, would be preferable to Q-Tpeak
ecg_analysis.Beattimes = beattimes(QC==1);
ecg_analysis.Baseline = baseline(QC==1);
ecg_analysis.Rampevents = r_amp_events(QC==1);
ecg_analysis.Tampevents = t_amp_events(QC==1);
ecg_analysis.Potentialectopics = ectopic_events(QC==1);
ecg_analysis.TPars.prepeak_t30_loc = t_split.prepeak_t30_loc(QC==1);
ecg_analysis.TPars.postpeak_t30_loc = t_split.postpeak_t30_loc(QC==1);
ecg_analysis.TPars.prepeak_t50_loc = t_split.prepeak_t50_loc(QC==1);
ecg_analysis.TPars.postpeak_t50_loc = t_split.postpeak_t50_loc(QC==1);
ecg_analysis.TPars.prepeak_t70_loc = t_split.prepeak_t70_loc(QC==1);
ecg_analysis.TPars.postpeak_t70_loc = t_split.postpeak_t70_loc(QC==1);
ecg_analysis.QC = QC;






end








function rloc = find_r_loc(ecg)

    % Wavelet transform implemented from:
    % \Holter_fingerprint_code\VCCRI_Sarah_Version\adr_ohara\src\generate_rr.m
    % Original code does not use findpeaks function, uses thresholdng with
    % specific constants depending on the dataset of the holter recording.
    
    cof=cwt(ecg,3,'coif1');

    % R -Peak Identification
    % cofsq=cof.^2;
    % % threshold=mean(cofsq);
    % threshold = max(cofsq)/2;
    % [~, rloc] = findpeaks(cofsq, 'MinPeakHeight', threshold); %

    height_threshold = max(cof)*0.4;
    [~,rloc] = findpeaks(cof, 'MinPeakHeight', height_threshold);
    distance_threshold = mean(diff(rloc))/10; % Add additional thresholding to prevent 'immediate' beats from being picked up
    [~,rloc] = findpeaks(cof, 'MinPeakHeight', height_threshold, 'MinPeakDistance', distance_threshold);

end



function heart_rate_change = heart_rate_change_plot(rr_int)
    window = 10; % Number of previous beats to average interval
    heart_rate_change = zeros(1,length(rr_int)-1);
    heart_rate_change(1:window-1) = diff(rr_int(1:window));
    for i = window:length(rr_int)
        heart_rate_change(i) = rr_int(i) - mean(heart_rate_change(i-(window-1):i-1)); 
    end
   plot(heart_rate_change); 
end



function interp_x_loc = interpolate_for_x(yloc, xinterp1, yinterp1, xinterp2, yinterp2) 
    % interpolation ranges (assuming constant sampling interval).
    % Looking for the interpolated x location at the desired y-value (yloc), given the datapoint before (xinterp1, yinterp1) and after (xinterp2, yinterp2).
    interp_x_loc = xinterp1 + (((xinterp2 - xinterp1) .* (yloc - yinterp1)) ./ (yinterp2-yinterp1)); % Rearranging the formula for Linear interpretation so x = x1 + ( (x2-x1)(y-y1) / (y2-y1) ).
end