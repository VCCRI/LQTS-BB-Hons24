function [ecg_analysis] = LQTS_Program_Function(file,name,named_out_folder, out_folder)

% Following logic of original LQTS_Program_Function, created by saa, but 
% defining ecg_data/ecg_analysis structures and condensing annotation/beat
% detection functions

qtpeak_lower_bound = 200; % PARAMETER: lower bound for qtpeak plot data in ms
qtpeak_upper_bound = 750; % PARAMETER: upper bound for qtpeak plot data in ms
rr_lower_bound = 400;   % PARAMETER: lower bound (higher heart rate) for the plot data in ms.
rr_upper_bound = 2000;   % PARAMETER: Upper bound (lower heart rate) for the plot data in ms.
zoom_mode = 1;          % PARAMETER: zoom_mode = 1 (default) zooms into the RR data range, zoom_mode = 2 does not zoom (400-2000ms)

if zoom_mode ~= 1 & zoom_mode ~= 2
    warning("zoom_mode parameter must be 1 or 2");
    zoom_mode = 1;
else
end

if rr_lower_bound < 400
    rr_lower_bound = 400; % 400ms = 150bpm.
else
end

if rr_upper_bound > 2000
    rr_upper_bound = 2000; % 2000ms = 30bpm.
else
end

%% LOAD IN ARASH FORMATTED ECG
ecg_data = load(file);
sampling_rate = ecg_data.arash_Header.Sampling_Rate;


%% Draw out points of analysis from ecg_data (RR interval, R amplitude, etc)

leadn = 1; % Index of holter lead to analyse
show_figure = 1; % Indicate whether to plot 
ecg_analysis = analyse_ecg(ecg_data,leadn,show_figure, name);

% Save the filtered ecg
figure
plot(ecg_analysis.filtered_ecg)
title(['Filtered ECG - Recording ' name]);
saveas(gcf,strcat(named_out_folder,filesep,'Filtered ECG'));
close all

%% Prepare for analysis:
rr_int = ecg_analysis.RRint;
sec_rr_int = ecg_analysis.RRint ./ 1000;
qtp_int = ecg_analysis.QTpeak;
beat_time_stamp = ecg_analysis.Beattimes;

% Prepare to being able to QC the day (8am-8pm) or night (1am-6am)
day_start = 8*60*60*1000; %% 8 am in milliseconds from midnight
day_end = 20*60*60*1000; %% 8 pm in milliseconds from midnight
night_start = 1*60*60*1000; %% 1 am in milliseconds from midnight
night_end = 6*60*60*1000; %% 6 am in milliseconds from midnight

time = beat_time_stamp((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound));






%% Cycle through 'potential events'

% for index = 1:numel(ecg_analysis.Rampevents)
%     if ecg_analysis.Rampevents(index) == 2
%         plot(ecg_data.raw_ecgSig(ecg_analysis.Beatstart(index):ecg_analysis.Beatstart(index+1),leadn))
%         input('PRESS ENTER TO CONTINUE...')
%     end
% end

%% Find ectopic beats (comapre Ramp to detected Tamp and check if it's a potential ectopic)

% for index = 1:numel(ecg_analysis.Rampevents)
%     if ecg_analysis.Rampevents(index) == 2
%         plot(ecg_data.raw_ecgSig(ecg_analysis.Beatstart(index):ecg_analysis.Beatstart(index+1),leadn))
%         input('PRESS ENTER TO CONTINUE...')
%     end
% end

%% TODO: Find and plot longest QT interval for each frequency band


% f_bins = 0.5:0.1:1.5;
% RR_frequency = ecg_analysis.RR/1000;
% for f = 1:(numel(f_bins)-1)
%     beat_frequency_index =  find(f_bins(f) <= RR_frequency & RR_frequency <= f_bins(f+1));
%     [maxQT, i] = max(ecg_analysis.QTpeak(beat_frequency_index),[],"all");
%     max_QT_i = beat_frequency_index(i);
%     raw_ecg = ecg_data.raw_ecgSig(:,leadn);
%     snip = raw_ecg(ecg_analysis.Beatstart(max_QT_i):ecg_analysis.Beatstart(max_QT_i+1));
%     time_x = 0:(1000/sampling_rate):(numel(snip)-1)*1000/sampling_rate;
%     plot(time_x,snip);
%     hold on
%     title('Max QT of Freq bin: '+string(f_bins(f))+'-'+string(f_bins(f+1))+'Hz')
%     xlabel('Time (ms)');
%     ylabel('Amplitude');
%     r = abs(ecg_analysis.Rloc(max_QT_i)- ecg_analysis.Beatstart(max_QT_i));
%     if r == 0
%         r = 1;
%     end
%     plot(r*1000/sampling_rate,snip(r),'r*');
%     text(r*1000/sampling_rate,snip(r),'R peak')
%     t = abs(ecg_analysis.Tloc(max_QT_i)- ecg_analysis.Beatstart(max_QT_i));
%     plot(t*1000/sampling_rate,snip(t),'b*');
%     text(t*1000/sampling_rate,snip(t),'T peak')
%     annotation('textbox',[.5 .5 .3 .3], 'String',...
%         ['Largest QT interval is '+string(ecg_analysis.QTpeak(max_QT_i))+' ms ',...
%         'at time: '+string(ecg_analysis.Beatstart(max_QT_i))+' ms'],...
%         'FitBoxToText','on');
%     drawnow;
%     % input('PRESS ENTER TO CONTINUE...')
%     hold off
% 
%    close all
% end






% %% Plotting RR/NormRamp Heatmap
% clear total_plot_data day_plot_data night_plot_data;
% 
% if size(ecg_analysis.Ramp,2) == size(ecg_analysis.RRint,2) + 1
%     R_amplitude = ecg_analysis.Ramp(1:end-1);
% elseif size(ecg_analys.Ramp,2) == size(ecg_analysis.RRint,2)
%     R_amplitude = ecg_analysis.Ramp;
% end
% norm_R_amplitude = normalize(R_amplitude, 'center', 'median', 'scale', 'iqr');  % Normalising to the median and IQR
% 
% 
% total_plot_data(:,1) = rr_int(rr_int <= rr_upper_bound & rr_int >= rr_lower_bound);
% total_plot_data(:,2) = norm_R_amplitude(rr_int <= rr_upper_bound & rr_int >= rr_lower_bound);
% 
% total_plot_data(total_plot_data(:,2) < -10 | total_plot_data(:,2) > 10, :) = []; % Removes values outside of defined smoothhist2d range 
% total_plot_data(any(isnan(total_plot_data),2),:) = []; % Removes all rows with at least 1 NaN
% figure("WindowState", "maximized")
% if zoom_mode == 2
%     [inlier_data, ecg_analysis.Feiler.RR_NRamp, RRxNRamp_perimeter] = smoothhist2D(total_plot_data,[400,-10],[2000,10],3,[(2000-400)/5, 200],.05, 'image', 0);  % Save F from smoothhist2d into ecg_analysis struc.
% 
% else
%     [inlier_data, ecg_analysis.Feiler.RR_NRamp, RRxNRamp_perimeter] = smoothhist2D(total_plot_data,[rr_lower_bound,-10],[rr_upper_bound,10],3,[(rr_upper_bound - rr_lower_bound)/5, 200],.05, 'image', 0);  % Save F from smoothhist2d into ecg_analysis struc.
% end
% set(gca, "XDir", "reverse"); % Only for RR intervals. Reverses x axis to be in descending RR (ascending HR).
% colormap('turbo');
% cb = colorbar;
% cb.Title.String = "Beat Density (8-bit color)";
% xlabel('RR interval (ms)');
% ylabel('median-normalised R amplitude');
% title(['RR interval vs R amplitude full plot - Recording ' name]);
% subtitle(strcat(num2str(rr_lower_bound), "<RR<", num2str(rr_upper_bound), " (ms)          Accepted Beats = ", num2str(size(inlier_data,1)), " of ", num2str(size(rr_int,2)), " beats"))
% % % print('-clipboard', '-dmeta')
% % % saveas(gcf,strcat(named_out_folder,filesep,'RR interval vs Normalised R amplitude Full'));
% % % saveas(gcf, strcat(named_out_folder,filesep,'RR interval vs Normalised R amplitude Full.png'));
% inliers.RR_NRamp = inlier_data;
% clear inlier_data
% clear total_plot_data day_plot_data night_plot_data;
% 
% 
% 
% 
% %% Plotting RR/Ramp Heatmap
% clear total_plot_data day_plot_data night_plot_data;
% 
% if size(ecg_analysis.Ramp,2) == size(ecg_analysis.RRint,2) + 1
%     R_amplitude = ecg_analysis.Ramp(1:end-1);
% elseif size(ecg_analys.Ramp,2) == size(ecg_analysis.RRint,2)
%     R_amplitude = ecg_analysis.Ramp;
% end
% 
% 
% % Total
% total_plot_data(:,1) = rr_int(rr_int <= rr_upper_bound & rr_int >= rr_lower_bound);
% total_plot_data(:,2) = R_amplitude(rr_int <= rr_upper_bound & rr_int >= rr_lower_bound);
% 
% total_plot_data(total_plot_data(:,2) < 0 | total_plot_data(:,2) > 350, :) = []; 
% total_plot_data(any(isnan(total_plot_data),2),:) = []; 
% figure("WindowState", "maximized")
% if zoom_mode == 2
%     [inlier_data, ecg_analysis.Feiler.RR_Ramp, RRxRamp_perimeter] = smoothhist2D(total_plot_data,[400,0],[2000,350],3,[(2000-400)/5, 350],.05, 'image', 0);
% else
%     [inlier_data, ecg_analysis.Feiler.RR_Ramp, RRxRamp_perimeter] = smoothhist2D(total_plot_data,[rr_lower_bound,0],[rr_upper_bound,350],3,[(rr_upper_bound - rr_lower_bound)/5, 350],.05, 'image', 0);
% end
% set(gca, "XDir", "reverse"); % Only for RR intervals.
% colormap('turbo');
% %Linear Regression
% hold on
% lm = fitlm(inlier_data(:,1), inlier_data(:,2));
% h = plot(lm); % Plots linear regression line
% h(1).Visible = 0; % Removes the plotted points (since the heatmap has already plotted them)
% legend("", "Linear Fit", "Confidence Bounds", "", "Location", "northeast");
% cb = colorbar;
% cb.Title.String = "Beat Density (8-bit color)";
% xlabel('RR interval (ms)');
% ylabel('R amplitude (units??)');
% title(['RR interval vs R amplitude full plot - Recording ' name]);
% subtitle(strcat(num2str(rr_lower_bound), "<RR<", num2str(rr_upper_bound), " (ms)          Accepted Beats = ", num2str(size(inlier_data,1)), " of ", num2str(size(rr_int,2)), " beats"));

% inliers.RR_Ramp = inlier_data;
% clear lm h eqn inlier_data
% hold off
% %% saveas(gcf, strcat(named_out_folder,filesep,'Regression RR vs Ramp Full.png'));
% % input('Press ENTER to continue...');
% 
% 

% 
% %% Plotting RR/RampTamp Heatmap (See TotalRec of RR/Ramp for explanations)
% clear total_plot_data day_plot_data night_plot_data;
% 
% T_amplitude = ecg_analysis.Tamp;
%     T_amplitude(T_amplitude == 0) = NaN;
% R_amplitude = ecg_analysis.Ramp;
%     R_amplitude(R_amplitude == 0) = NaN;
% RampTamp = R_amplitude ./ T_amplitude;
% 
% 
% 
% 
% 
% % Total
% total_plot_data(:,1) = rr_int((rr_int >=rr_lower_bound) & (rr_int <=rr_upper_bound));
% total_plot_data(:,2) = RampTamp((rr_int >=rr_lower_bound) & (rr_int <=rr_upper_bound));
% 
% total_plot_data(total_plot_data(:,2) < -25 | total_plot_data(:,2) > 25, :) = []; 
% total_plot_data(any(isnan(total_plot_data),2),:) = [];
% figure("WindowState", "maximized")
% if zoom_mode == 2
%     [inlier_data,ecg_analysis.Feiler.RR_RampTamp, RRxRampTamp_perimeter] = smoothhist2D(total_plot_data,[400,-25],[2000,25],3,[(2000-400)/5, 500],.05, 'image', 0);
% else
%     [inlier_data,ecg_analysis.Feiler.RR_RampTamp, RRxRampTamp_perimeter] = smoothhist2D(total_plot_data,[rr_lower_bound,-25],[rr_upper_bound,25],3,[(rr_upper_bound - rr_lower_bound)/5, 500],.05, 'image', 0);
% end
% set(gca, "XDir", "reverse"); % Only for RR intervals. 
% colormap('turbo');
% %Linear Regression
% hold on
% lm = fitlm(inlier_data(:,1), inlier_data(:,2));
% h = plot(lm);
% h(1).Visible = 0;
% legend("", "Linear Fit", "Confidence Bounds", "", "Location", "northeast");
% cb = colorbar;
% cb.Title.String = "Beat Density (8-bit color)"; % Scale is 0 to 255 same as 8-bit colour.
% xlabel('RR interval (ms)');
% ylabel('Ramp/Tamp');
% title(['RR interval vs Ramp/Tamp full plot, Recording ' name]);
% subtitle(strcat(num2str(rr_lower_bound), "<RR<",num2str(rr_upper_bound), " (ms)          Accepted Beats = ", num2str(size(inlier_data,1)), " of ", num2str(size(rr_int,2)), " beats"));
% inliers.RR_RampTamp = inlier_data;
% clear lm h eqn inlier_data
% hold off
% %% saveas(gcf, strcat(named_out_folder,filesep,'Regression RR vs RampTamp Full.png'));
% % input('Press ENTER to continue...');
% 
% 
% 
% 
% 
% 
% %% Plotting Tamp vs Ramp Heatmap (See TotalRec of RR/Ramp for explanations)
% clear total_plot_data day_plot_data night_plot_data;
% T_amplitude = ecg_analysis.Tamp;
%     T_amplitude(T_amplitude == 0) = NaN;
% R_amplitude = ecg_analysis.Ramp;
%     R_amplitude(R_amplitude == 0) = NaN;
% 
% 
% 
% 
% % Total
% total_plot_data(:,1) = T_amplitude((rr_int <= rr_upper_bound) & (rr_int >= rr_lower_bound));
% total_plot_data(:,2) = R_amplitude((rr_int <= rr_upper_bound) & (rr_int >= rr_lower_bound));
% total_plot_data(total_plot_data(:,1) < -100 | total_plot_data(:,1) > 200, :) = []; 
% total_plot_data(total_plot_data(:,2) < -50 | total_plot_data(:,2) > 350, :) = [];
% total_plot_data(any(isnan(total_plot_data),2),:) = []; % Removes all rows with at least 1 NaN
% figure("WindowState", "maximized")
%     [inlier_data, ecg_analysis.Feiler.Tamp_Ramp, TampxRamp_perimeter] = smoothhist2D(total_plot_data,[-100 -50],[200,350],3,[300, 400],.05, 'image', 0);
% colormap('turbo');
% %Linear Regression
% hold on
% lm = fitlm(inlier_data(:,1), inlier_data(:,2));
% h = plot(lm);
% h(1).Visible = 0;
% legend("", "Linear Fit", "Confidence Bounds", "", "Location", "northeast");
% cb = colorbar;
% cb.Title.String = "Beat Density (8-bit color)";
% xlabel('T amplitude (units?)');
% ylabel('R amplitude (units?)');
% title(['T amplitude vs R amplitude full plot, Recording '  name]);
% subtitle(strcat(num2str(rr_lower_bound), "<RR<", num2str(rr_upper_bound), " (ms)          Accepted Beats = ", num2str(size(inlier_data,1)), " of ", num2str(size(rr_int,2)), " beats"))
% % % print('-clipboard', '-dmeta')
% % % saveas(gcf,strcat(named_out_folder,filesep,'T amplitude vs R amplitude Full'));
% % % saveas(gcf, strcat(named_out_folder,filesep,'T amplitude vs R amplitude Full.png'));
% 
% %         if isnan(lm.Rsquared.Adjusted) | isnan(lm.ModelFitVsNullModel.Pvalue)
% %             eqn = sprintf('y = %.2fx + %.2f\nR^2 or p error', lm.Coefficients{2,1}, lm.Coefficients{1,1}); % formatting the equation text)
% %             input("Press ENTER to continue...")
% %         elseif lm.ModelFitVsNullModel.Pvalue >= 0.0001
% %              eqn = sprintf('y = %.2fx + %.3f\nR^2 = %.3f, p = %.3f', lm.Coefficients{2,1}, lm.Coefficients{1,1}, lm.Rsquared.Adjusted, lm.ModelFitVsNullModel.Pvalue);
% %         else
% %             eqn = sprintf('y = %.2fx + %.3f\nR^2 = %.3f, p <0.0001', lm.Coefficients{2,1}, lm.Coefficients{1,1}, lm.Rsquared.Adjusted);
% %         end
% 
% %         text(-90, -10, eqn, 'Color', 'r', 'FontSize', 10); % displaying the equation
% 
% inliers.Tamp_Ramp = inlier_data;
% clear lm h eqn inlier_data
% clear total_plot_data day_plot_data night_plot_data
% hold off
% %% saveas(gcf, strcat(named_out_folder,filesep,'Regression Tamp vs Ramp Full.png'));
% 
% 

%% Plotting RR/TampRamp Heatmap (See TotalRec of RR/Ramp for explanations)
clear total_plot_data day_plot_data night_plot_data;


T_amplitude = ecg_analysis.Tamp;
    T_amplitude(T_amplitude == 0) = NaN;
R_amplitude = ecg_analysis.Ramp;
    R_amplitude(R_amplitude == 0) = NaN;
TampRamp = T_amplitude ./ R_amplitude;


% Total
total_plot_data(:,1) = rr_int((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound));
total_plot_data(:,2) = TampRamp((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound));

total_plot_data(total_plot_data(:,2) < -0.25 | total_plot_data(:,2) > 1.25, :) = []; 
total_plot_data(any(isnan(total_plot_data),2),:) = [];
ecg_analysis.PD.RRxTampRamp = total_plot_data;
figure %("WindowState", "maximized")
if zoom_mode == 2
    [inlier_data, ecg_analysis.Feiler.RR_TampRamp, RRxTampRamp_perimeter] = smoothhist2D(total_plot_data,[400,-0.25],[2000,1.25],3,[(2000-400)/5, 350],.05, 'image', 0);
else
    [inlier_data, ecg_analysis.Feiler.RR_TampRamp, RRxTampRamp_perimeter] = smoothhist2D(total_plot_data,[rr_lower_bound,-0.2],[rr_upper_bound,0.6],3,[(rr_upper_bound - rr_lower_bound)/5, 80],.05, 'image', 0);
end
set(gca, "XDir", "reverse"); % Only for RR intervals.
colormap('turbo');


%Linear Regression
hold on
lm = fitlm(inlier_data(:,1), inlier_data(:,2));
h = plot(lm);
h(1).Visible = 0;
legend("", "Linear Fit", "Confidence Bounds", "", "Location", "northeast");
cb = colorbar;
cb.Title.String = "Beat Density (8-bit color)";
xlabel('RR interval (ms)');
ylabel('Tamp/Ramp');
title(['RR interval vs Tamp/Ramp full plot, Recording ' name]);
subtitle(strcat(num2str(rr_lower_bound), "<RR<", num2str(rr_upper_bound), " (ms)          Accepted Beats = ", num2str(size(inlier_data,1)), " of ", num2str(size(rr_int,2)), " beats"))

inliers.RR_TampRamp = inlier_data;

RRxTampRamp_m = lm.Coefficients.Estimate(2);
RRxTampRamp_b = lm.Coefficients.Estimate(1);
RRxTampRamp_rsquare = lm.Rsquared.Ordinary;

clear lm h eqn inlier_data
hold off
saveas(gcf, strcat(named_out_folder,filesep,'Regression RR vs TampRamp Full.png'));
close all
% input('Press ENTER to continue...');



%% Plotting RR/QTpeak interval Heatmap (See TotalRec of RR/Ramp for explanations)

clear total_plot_data day_plot_data night_plot_data;




% Total
total_plot_data(:,1) = sec_rr_int((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound));
total_plot_data(:,2) = qtp_int((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound));

total_plot_data(total_plot_data(:,2) < qtpeak_lower_bound | total_plot_data(:,2) > qtpeak_upper_bound, :) = []; 
total_plot_data(any(isnan(total_plot_data),2),:) = []; % Removes all rows with at least 1 NaN

ecg_analysis.PD.RRxQTp = total_plot_data;
figure %("WindowState", "maximized")
if zoom_mode == 2
    [inlier_data, ecg_analysis.Feiler.RR_QTp, RRxQTp_perimeter] = smoothhist2D(total_plot_data,[0.4,200],[2,750],3,[(1000 * (2-0.4))/5, (qtpeak_upper_bound-qtpeak_lower_bound)/5], 0.05, 'image', 0);
else
    [inlier_data, ecg_analysis.Feiler.RR_QTp, RRxQTp_perimeter] = smoothhist2D(total_plot_data,[rr_lower_bound/1000,qtpeak_lower_bound],[rr_upper_bound/1000,qtpeak_upper_bound],3,[(rr_upper_bound - rr_lower_bound)/5, (qtpeak_upper_bound-qtpeak_lower_bound)/5],.05, 'image', 0);
end
set(gca, "XDir", "reverse"); % Only for RR intervals.
colormap('turbo');
hold on



%% Fitting a Power function 
[RRxQT.fiteq, RRxQT.gof] = fit(inlier_data(:,1), inlier_data(:,2), 'power1');
eq = ['QTint = ', num2str(RRxQT.fiteq.a), ' * RR^{', num2str(RRxQT.fiteq.b), '}'];
% Getting the p values for the fit
RRxQT.coeffs = coeffvalues(RRxQT.fiteq);
RRxQT.confint = confint(RRxQT.fiteq); % CI
RRxQT.se = (RRxQT.confint(2,:) - RRxQT.confint(2,:)) / (2 * 1.96);% Standard error calculated from the CI
RRxQT.tstats = RRxQT.coeffs ./ RRxQT.se;
RRxQT.df = length(inlier_data(:,1)) - length(RRxQT.coeffs);
RRxQT.pvals = 2 * (1 - tcdf(abs(RRxQT.tstats), RRxQT.df));
%Preparing the CI shadow
xconf = linspace((rr_lower_bound/1000),(rr_upper_bound/1000), ((rr_upper_bound-rr_lower_bound) /5))';
low_ci = RRxQT.confint(1,1) .* (xconf .^ RRxQT.confint(1,2));
high_ci = RRxQT.confint(2,1) .* (xconf .^ RRxQT.confint(2,2));
yconf = [low_ci; flip(high_ci)];
xconf = [xconf; flip(xconf)];

hold on
fill(xconf, yconf, [1 0.8 0.8], "FaceAlpha", '0.75', "EdgeColor", "none"); % Plot 95CI behind the fitted eq
plot(RRxQT.fiteq); % Plot fitted eq

legend("95% Confidence Interval", eq, "Location", "northeast");
cb = colorbar;
cb.Title.String = "Beat Density (8-bit color)";
xlabel('RR interval (s)');
ylabel('Q-Tpeak interval (ms)');
title(['RR interval vs Q-Tpeak interval full plot, Recording ' name]);
subtitle(strcat(num2str(rr_lower_bound), "<RR<", num2str(rr_upper_bound), " (ms)          Accepted Beats = ", num2str(size(inlier_data,1)), " of ", num2str(size(rr_int,2)), " beats"))
% print('-clipboard', '-dmeta')
% saveas(gcf,strcat(named_out_folder,filesep,'RR interval vs Q-Tpeak int Full'));
% saveas(gcf, strcat(named_out_folder,filesep,'RR interval vs Q-Tpeak int Full.png'));

inliers.secRR_QTpeak = inlier_data;
hold off
clear lm h eqn inlier_data eq xconf low_ci high_ci yconf
hold off
saveas(gcf, strcat(named_out_folder,filesep,'Regression RR vs QTpeak int Full.png'));
close all
% input('Press ENTER to continue...');





%% Plotting QTpeak vs Tamp/Ramp
clear total_plot_data day_plot_data night_plot_data;

T_amplitude = ecg_analysis.Tamp;
    T_amplitude(T_amplitude == 0) = NaN;
R_amplitude = ecg_analysis.Ramp;
    R_amplitude(R_amplitude == 0) = NaN;
TampRamp = T_amplitude ./ R_amplitude;



% Total
total_plot_data(:,1) = qtp_int((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound));
total_plot_data(:,2) = TampRamp((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound));
total_plot_data(total_plot_data(:,1) < 200 | total_plot_data(:,1) > 750, :) = []; 
total_plot_data(total_plot_data(:,2) < -0.25 | total_plot_data(:,2) > 1.25, :) = []; 
total_plot_data(any(isnan(total_plot_data),2),:) = [];
ecg_analysis.PD.QTpxTampRamp = total_plot_data;
figure %("WindowState", "maximized")
[inlier_data, ecg_analysis.Feiler.QTp_TampRamp, QTpxTampRamp_perimeter] = smoothhist2D(total_plot_data,[qtpeak_lower_bound,-0.2],[qtpeak_upper_bound,0.6],3,[(qtpeak_upper_bound-qtpeak_lower_bound)/5, 80],.05, 'image', 0);
colormap('turbo');
%Linear Regression
hold on
lm = fitlm(inlier_data(:,1), inlier_data(:,2));
h = plot(lm);
h(1).Visible = 0;
legend("", "Linear Fit", "Confidence Bounds", "", "Location", "northeast")
cb = colorbar;
cb.Title.String = "Beat Density (8-bit color)";
xlabel('Q-Tpeak Interval (ms)');
ylabel('T amplitude/R amplitude');
title(['Q-Tpeak vs Tamp/Ramp full plot, Recording ' name]);
subtitle(strcat(num2str(rr_lower_bound), "<RR<", num2str(rr_upper_bound), " (ms)          Accepted Beats = ", num2str(size(inlier_data,1)), " of ", num2str(size(rr_int,2)), " beats"))
% % print('-clipboard', '-dmeta')
% % saveas(gcf,strcat(named_out_folder,filesep,'QTpeak vs TampRamp Full'));
saveas(gcf, strcat(named_out_folder,filesep,'QTpeak vs TampRamp Full.png'));
close all
inliers.QTpeak_TampRamp = inlier_data;
QTpxTampRamp_m = lm.Coefficients.Estimate(2);
QTpxTampRamp_b = lm.Coefficients.Estimate(1);
QTpxTampRamp_rsquare = lm.Rsquared.Ordinary;
clear lm h eqn inlier_data
hold off








%% Plotting QTpeak/T_amp
clear total_plot_data day_plot_data night_plot_data;

T_amplitude = ecg_analysis.Tamp;


% Total
total_plot_data(:,1) = qtp_int((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound));
total_plot_data(:,2) = T_amplitude((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound));
total_plot_data(total_plot_data(:,1) < 200 | total_plot_data(:,1) > 750, :) = []; 
total_plot_data(total_plot_data(:,2) < -50 | total_plot_data(:,2) > 200, :) = []; 
total_plot_data(any(isnan(total_plot_data),2),:) = []; 
figure("WindowState", "maximized")
[inlier_data, ecg_analysis.Feiler.QTp_Tamp, QTpxTamp_perimeter] = smoothhist2D(total_plot_data,[qtpeak_lower_bound,-50],[qtpeak_upper_bound,200],3,[(qtpeak_upper_bound-qtpeak_lower_bound)/5, 250],.05, 'image', 0);
colormap('turbo');
%Linear Regression
hold on
lm = fitlm(inlier_data(:,1), inlier_data(:,2));
h = plot(lm);
h(1).Visible = 0;
legend("", "Linear Fit", "Confidence Bounds", "", "Location", "northeast")
cb = colorbar;
cb.Title.String = "Beat Density (8-bit color)";
xlabel('Q-Tpeak Interval (ms)');
ylabel('T amplitude (units?)');
title(['Q-Tpeak vs T amplitude full plot, Recording ' name]);
subtitle(strcat(num2str(rr_lower_bound), "<RR<", num2str(rr_upper_bound), " (ms)          Accepted Beats = ", num2str(size(inlier_data,1)), " of ", num2str(size(rr_int,2)), " beats"))

inliers.QTp_Tamp = inlier_data;
clear lm h eqn inlier_data
hold off
% saveas(gcf, strcat(named_out_folder,filesep,'Regression QTpeak int vs Tamp Full.png'));





%% Plotting RR/Tpeak_Td70 interval Heatmap
clear total_plot_data day_plot_data night_plot_data;

Tp_Td70 = (ecg_analysis.TPars.postpeak_t70_loc - ecg_analysis.Tloc) .* (1000 / sampling_rate); % Tpeak to downsloping T70 distance.



% Total
total_plot_data(:,1) = sec_rr_int((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound));
total_plot_data(:,2) = Tp_Td70((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound));

total_plot_data(total_plot_data(:,2) < 0 | total_plot_data(:,2) > 150, :) = []; 
total_plot_data(any(isnan(total_plot_data),2),:) = []; 
ecg_analysis.PD.RRxTpTd70 = total_plot_data;
figure %("WindowState", "maximized")
if zoom_mode == 2
    [inlier_data, ecg_analysis.Feiler.RR_TpTd70, RRxTpTd70_perimeter] = smoothhist2D(total_plot_data,[0.4,0],[2,150],3,[(1000 * (2-0.4))/5, 150], 0.05, 'image', 0);
else
    [inlier_data, ecg_analysis.Feiler.RR_TpTd70, RRxTpTd70_perimeter] = smoothhist2D(total_plot_data,[rr_lower_bound/1000,0],[rr_upper_bound/1000,150],3,[(rr_upper_bound - rr_lower_bound)/5, 150],.05, 'image', 0);
end
set(gca, "XDir", "reverse"); % Only for RR intervals.
colormap('turbo');
hold on


%% Fitting a Power function 
[RRxTpTd70.fiteq, RRxTpTd70.gof] = fit(inlier_data(:,1), inlier_data(:,2), 'power1', "Lower", [0,0]);
eq = ['Tpeak-Td70 int = ', num2str(RRxTpTd70.fiteq.a), ' * RR^{', num2str(RRxTpTd70.fiteq.b), '}'];
% Getting the p values for the fit
RRxTpTd70.coeffs = coeffvalues(RRxTpTd70.fiteq);
RRxTpTd70.confint = confint(RRxTpTd70.fiteq); % CI
RRxTpTd70.se = (RRxTpTd70.confint(2,:) - RRxTpTd70.confint(2,:)) / (2 * 1.96);% Standard error calculated from the CI
RRxTpTd70.tstats = RRxTpTd70.coeffs ./ RRxTpTd70.se;
RRxTpTd70.df = length(inlier_data(:,1)) - length(RRxTpTd70.coeffs);
RRxTpTd70.pvals = 2 * (1 - tcdf(abs(RRxTpTd70.tstats), RRxTpTd70.df));
%Preparing the CI shadow
xconf = linspace((rr_lower_bound/1000),(rr_upper_bound/1000), ((rr_upper_bound-rr_lower_bound) /5))';
low_ci = RRxTpTd70.confint(1,1) .* (xconf .^ RRxTpTd70.confint(1,2));
high_ci = RRxTpTd70.confint(2,1) .* (xconf .^ RRxTpTd70.confint(2,2));
yconf = [low_ci; flip(high_ci)];
xconf = [xconf; flip(xconf)];

hold on
fill(xconf, yconf, [1 0.8 0.8], "FaceAlpha", '0.75', "EdgeColor", "none"); % Plot 95CI behind the fitted eq
plot(RRxTpTd70.fiteq); % Plot fitted eq

legend("95% Confidence Interval", eq, "Location", "northeast");
cb = colorbar;
cb.Title.String = "Beat Density (8-bit color)";
xlabel('RR interval (s)');
ylabel('Tpeak-Td70 interval (ms)');
title(['RR interval vs Tpeak-Td70 interval full plot, Recording ' name]);
subtitle(strcat(num2str(rr_lower_bound), "<RR<", num2str(rr_upper_bound), " (ms)          Accepted Beats = ", num2str(size(inlier_data,1)), " of ", num2str(size(rr_int,2)), " beats"))
inliers.secRR_TpTd70 = inlier_data;
clear lm h eqn inlier_data eq xconf low_ci high_ci yconf
hold off



saveas(gcf, strcat(named_out_folder,filesep,'Regression RR vs Tpeak-Td70 int Full.png'));
close all
% input('Press ENTER to continue...');


clear total_plot_data day_plot_data night_plot_data;







%% HISTOGRAMS

% Histogram for QTpeak
figure('WindowState', "maximized")
histogram(ecg_analysis.QTpeak((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)), [qtpeak_lower_bound:5:qtpeak_upper_bound], 'FaceColor', 'black', 'FaceAlpha', 0.3, "EdgeAlpha", 0.3); % All Beats
hold on
histogram(ecg_analysis.QTpeak((rr_int <= 1000) & (rr_int >=909)), [qtpeak_lower_bound:5:qtpeak_upper_bound], "FaceColor", "blue", "FaceAlpha", 0.4, "EdgeAlpha", 0.4);% Slow Beats
histogram(ecg_analysis.QTpeak((rr_int <= 714) & (rr_int >=666)), [qtpeak_lower_bound:5:qtpeak_upper_bound], "FaceColor", "red", "FaceAlpha", 0.6, "EdgeAlpha", 0.4); % Fast Beats
legend("All beats (30-150bpm)", "Slow Heart Rate (60-66bpm)", "Fast Heart Rate (84-90bpm)", "Location", "northeast")
xlabel("Q-Tpeak interval (ms)");
title(["Q-Tpeak interval histogram, Recording " name]);
% saveas(gcf,strcat(named_out_folder,filesep,'QTpeak interval histogram.png'));
close all
hold off


% Histogram for RR interval (ms)

histogram(ecg_analysis.RRint, [400:40:2000]);
hold on
histogram(ecg_analysis.RRint((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)), [400:40:2000]); % All Beats
xlabel("RR interval (ms)");
title(["RR interval histogram, Recording " name]);
% % saveas(gcf,strcat(named_out_folder,filesep,'RR interval histogram.png'));
hold off

% % Histogram for Ramp
% histogram(ecg_analysis.Ramp, [-50:5:350]);
% hold on
% histogram(ecg_analysis.Ramp((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)), [-50:5:350]); % All Beats
% xlabel("R amplitude (units?)");
% title(["R amplitude histogram, Recording " name]);
% % % saveas(gcf,strcat(named_out_folder,filesep,'R amplitude histogram.png'));
% hold off
% 
% % Histogram for Tamp
% histogram(ecg_analysis.Tamp, [-100:3:200]);
% hold on
% histogram(ecg_analysis.Tamp((rr_int <=rr_upper_bound) & (rr_int>=rr_lower_bound)), [-100:3:200]); % All Beats
% xlabel("T amplitude (units)");
% title(["T amplitude histogram, Recording " name]);
% % % saveas(gcf,strcat(named_out_folder,filesep,'T amplitude histogram.png'));
% hold off


% % Histogram for TampRamp
figure('WindowState', "maximized")
histogram(TampRamp((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)), [-0.25:0.025:1.25], 'FaceColor', 'black', 'FaceAlpha', 0.3, "EdgeAlpha", 0.3); % All Beats
hold on
histogram(TampRamp((rr_int <= 1000) & (rr_int >=909)), [-0.25:0.025:1.25], "FaceColor", "blue", "FaceAlpha", 0.4, "EdgeAlpha", 0.4);% Slow Beats
histogram(TampRamp((rr_int <= 714) & (rr_int >=666)), [-0.25:0.025:1.25], "FaceColor", "red", "FaceAlpha", 0.6, "EdgeAlpha", 0.4);% Fast Beats
legend("All beats (30-150bpm)", "Slow Heart Rate (60-66bpm)", "Fast Heart Rate (84-90bpm)", "Location", "northeast")
xlabel("T amplitude/R amplitude");
title(["T amplitude/R amplitude histogram, Recording " name]);
% saveas(gcf,strcat(named_out_folder,filesep,'TampRamp histogram.png'));
close all
hold off

% % Histogram for TpTd70
figure('WindowState', "maximized")
histogram(Tp_Td70((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)), [0:1.25:150], 'FaceColor', 'black', 'FaceAlpha', 0.3, "EdgeAlpha", 0.3);% All Beats
hold on
histogram(Tp_Td70((rr_int <= 1000) & (rr_int >=909)), [0:1.25:150], "FaceColor", "blue", "FaceAlpha", 0.4, "EdgeAlpha", 0.4);% Slow Beats
histogram(Tp_Td70((rr_int <= 714) & (rr_int >=666)), [0:1.25:150], "FaceColor", "red", "FaceAlpha", 0.6, "EdgeAlpha", 0.4);% Fast Beats
legend("All beats (30-150bpm)", "Slow Heart Rate (60-66bpm)", "Fast Heart Rate (84-90bpm)", "Location", "northeast")
xlabel("T_p_e_a_k to T_d_,_7_0 interval (ms)");
title(["T_p_e_a_k to T_d_,_7_0 interval histogram, Recording " name]);
% saveas(gcf,strcat(named_out_folder,filesep,'Tp to Td70 histogram.png'));
close all
hold off


% In order: RR, QT, TampRamp, TpTd70. Inlier data only.
in_RR = inliers.secRR_QTpeak(:,1) .* 1000;
in_QTp = inliers.secRR_QTpeak(:,2);
in_TampRamp = inliers.RR_TampRamp(:,2);
in_TpTd70 = inliers.secRR_TpTd70(:,2);



%% START REPOL SAVING % If uncommenting, ensure relevant code also uncommented in BATCH_ARASH_FORMAT_ANALYSIS
%% TODO: FIgure out another way to save these so that SSIM_recs later on will not rely on magic index numbers and also such that the pre-assigned zeros for repol_vars is always the right size.
%% If changing the size on this before fixing please also change the zeros size for repol_vars in BATCH_ARASH_FORMAT_ANALYSIS and the locations of the SSIM_rec calls.
% %% Load Items we want to save
% 
% % Col 1: Name of recording
% recording = str2double(name);
% % Col 2: Avg QTpeak
% meanQTpeak = mean(in_QTp);
% % Col 3: Median QTpeak
% medianQTpeak = median(in_QTp);
% % Col 4: a from RRxQTpeak heatmap (ax^b)
% a_from_RRxQT = RRxQT.fiteq.a;
% % Col 5: b from RRxQTpeak heatmap (ax^b)
% b_from_RRxQT = RRxQT.fiteq.b;
% 
% % Col 6: avg Tu30Tp distance
% mean_Tu30Tp = mean(ecg_analysis.Tloc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)) - ecg_analysis.TPars.prepeak_t30_loc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% % Col 7: avg TpTd30 distance
% mean_TpTd30 = mean(ecg_analysis.TPars.postpeak_t30_loc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)) - ecg_analysis.Tloc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% % Col 8: avg Tu50Tp distance
% mean_Tu50Tp = mean(ecg_analysis.Tloc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)) - ecg_analysis.TPars.prepeak_t50_loc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% % Col 9: avg TpTd50 distance
% mean_TpTd50 = mean(ecg_analysis.TPars.postpeak_t50_loc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)) - ecg_analysis.Tloc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% % Col 10: avg Tu70Tp distance
% mean_Tu70Tp = mean(ecg_analysis.Tloc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)) - ecg_analysis.TPars.prepeak_t70_loc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% % Col 11: avg TpTd70 distance
% mean_TpTd70 = mean(ecg_analysis.TPars.postpeak_t70_loc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)) - ecg_analysis.Tloc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% 
% % Col 12: median Tu30Tp distance
% median_Tu30Tp = median(ecg_analysis.Tloc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)) - ecg_analysis.TPars.prepeak_t30_loc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% % Col 13: median TpTd30
% median_TpTd30 = median(ecg_analysis.TPars.postpeak_t30_loc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)) - ecg_analysis.Tloc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% % Col 14: median Tu50Tp distance
% median_Tu50Tp = median(ecg_analysis.Tloc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)) - ecg_analysis.TPars.prepeak_t50_loc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% % Col 15: median TpTd50
% median_TpTd50 = median(ecg_analysis.TPars.postpeak_t50_loc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)) - ecg_analysis.Tloc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% % Col 16: median Tu70Tp distance
% median_Tu70Tp = median(ecg_analysis.Tloc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)) - ecg_analysis.TPars.prepeak_t70_loc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% % Col 17: median TpTd70
% median_TpTd70 = median(ecg_analysis.TPars.postpeak_t70_loc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)) - ecg_analysis.Tloc((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% 
% % Col 18: avg tamp
% meanTamp = mean(ecg_analysis.Tamp((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% % Col 19: median tamp
% medianTamp = median(ecg_analysis.Tamp((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% % Col 20: avg tamp(n+1) - Tamp(n)
% mean_tamp_difference = mean(diff(ecg_analysis.Tamp((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound))));
% % Col 21: median tamp(n+1) - Tamp(n)
% median_tamp_difference = median(diff(ecg_analysis.Tamp((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound))));
% 
% % Col 22: avg rr
% meanRR = mean(in_RR);
% % Col 23: median RR
% medianRR = median(in_RR);
% % Col 24: avg Ramp 
% meanRamp = mean(ecg_analysis.Ramp((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% % Col 25: median Ramp
% medianRamp = median(ecg_analysis.Ramp((rr_int <=rr_upper_bound) & (rr_int >=rr_lower_bound)));
% 
% % Col 26: Number of QC included beats
% QC_included = numel(in_RR);
% % Col 27: Number of QC excluded beats
% QC_excluded = numel(ecg_analysis.QC) - numel(in_RR);
% 
% % Col 28: Average QTint at 60-66bpm
% meanQTpeak_slow = mean(in_QTp((in_RR <=1000) & (in_RR >=909)));
% % Col 29: Median QTint at 60-66bpm
% medianQTpeak_slow = median(in_QTp((in_RR <=1000) & (in_RR >=909)));
% % Col 30: Average QTint at 84-90bpm
% meanQTpeak_fast = mean(in_QTp((in_RR <=714) & (in_RR >=666)));
% % Col 31: Median QTint at 84-90bpm
% medianQTpeak_fast = median(in_QTp((in_RR <=714) & (in_RR >=666)));
% 
% % Col 32: avg Tu30Tp distance at 60-66bpm
% mean_Tu30Tp_slow = mean(ecg_analysis.Tloc((rr_int <= 1000) & (rr_int >= 909)) - ecg_analysis.TPars.prepeak_t30_loc((rr_int <= 1000) & (rr_int >= 909)));
% % Col 33: avg TpTd30 at 60-66bpm
% mean_TpTd30_slow = mean(ecg_analysis.TPars.postpeak_t30_loc((rr_int <= 1000) & (rr_int >= 909)) - ecg_analysis.Tloc((rr_int <= 1000) & (rr_int >= 909)));
% % Col 34: avg Tu50Tp distance at 60-66bpm
% mean_Tu50Tp_slow = mean(ecg_analysis.Tloc((rr_int <= 1000) & (rr_int >= 909)) - ecg_analysis.TPars.prepeak_t50_loc((rr_int <= 1000) & (rr_int >= 909)));
% % Col 35: avg TpTd50 at 60-66bpm
% mean_TpTd50_slow = mean(ecg_analysis.TPars.postpeak_t50_loc((rr_int <= 1000) & (rr_int >= 909)) - ecg_analysis.Tloc((rr_int <= 1000) & (rr_int >= 909)));
% % Col 36: avg Tu70Tp distance at 60-66bpm
% mean_Tu70Tp_slow = mean(ecg_analysis.Tloc((rr_int <= 1000) & (rr_int >= 909)) - ecg_analysis.TPars.prepeak_t70_loc((rr_int <= 1000) & (rr_int >= 909)));
% % Col 37: avg TpTd70 at 60-66bpm
% mean_TpTd70_slow = mean(ecg_analysis.TPars.postpeak_t70_loc((rr_int <= 1000) & (rr_int >= 909)) - ecg_analysis.Tloc((rr_int <= 1000) & (rr_int >= 909)));
% 
% % Col 38: median Tu30Tp distance at 60-66bpm
% median_Tu30Tp_slow = median(ecg_analysis.Tloc((rr_int <= 1000) & (rr_int >= 909)) - ecg_analysis.TPars.prepeak_t30_loc((rr_int <= 1000) & (rr_int >= 909)));
% % Col 39: median TpTd30 at 60-66bpm
% median_TpTd30_slow = median(ecg_analysis.TPars.postpeak_t30_loc((rr_int <= 1000) & (rr_int >= 909)) - ecg_analysis.Tloc((rr_int <= 1000) & (rr_int >= 909)));
% % Col 40: median Tu50Tp distance at 60-66bpm
% median_Tu50Tp_slow = median(ecg_analysis.Tloc((rr_int <= 1000) & (rr_int >= 909)) - ecg_analysis.TPars.prepeak_t50_loc((rr_int <= 1000) & (rr_int >= 909)));
% % Col 41: median TpTd50 at 60-66bpm
% median_TpTd50_slow = median(ecg_analysis.TPars.postpeak_t50_loc((rr_int <= 1000) & (rr_int >= 909)) - ecg_analysis.Tloc((rr_int <= 1000) & (rr_int >= 909)));
% % Col 42: median Tu70Tp distance at 60-66bpm
% median_Tu70Tp_slow = median(ecg_analysis.Tloc((rr_int <= 1000) & (rr_int >= 909)) - ecg_analysis.TPars.prepeak_t70_loc((rr_int <= 1000) & (rr_int >= 909)));
% % Col 43: median TpTd70 at 60-66bpm
% median_TpTd70_slow = median(ecg_analysis.TPars.postpeak_t70_loc((rr_int <= 1000) & (rr_int >= 909)) - ecg_analysis.Tloc((rr_int <= 1000) & (rr_int >= 909)));
% 
% 
% % Col 44: avg Tu30Tp distance at 84-90bpm
% mean_Tu30Tp_fast = mean(ecg_analysis.Tloc((rr_int <= 714) & (rr_int >= 666)) - ecg_analysis.TPars.prepeak_t30_loc((rr_int <= 714) & (rr_int >= 666)));
% % Col 45: avg TpTd30
% mean_TpTd30_fast = mean(ecg_analysis.TPars.postpeak_t30_loc((rr_int <= 714) & (rr_int >= 666)) - ecg_analysis.Tloc((rr_int <= 714) & (rr_int >= 666)));
% % Col 46: avg Tu50Tp distance
% mean_Tu50Tp_fast = mean(ecg_analysis.Tloc((rr_int <= 714) & (rr_int >= 666)) - ecg_analysis.TPars.prepeak_t50_loc((rr_int <= 714) & (rr_int >= 666)));
% % Col 47: avg TpTd50
% mean_TpTd50_fast = mean(ecg_analysis.TPars.postpeak_t50_loc((rr_int <= 714) & (rr_int >= 666)) - ecg_analysis.Tloc((rr_int <= 714) & (rr_int >= 666)));
% % Col 48: avg Tu70Tp distance
% mean_Tu70Tp_fast = mean(ecg_analysis.Tloc((rr_int <= 714) & (rr_int >= 666)) - ecg_analysis.TPars.prepeak_t70_loc((rr_int <= 714) & (rr_int >= 666)));
% % Col 49: avg TpTd70
% mean_TpTd70_fast = mean(ecg_analysis.TPars.postpeak_t70_loc((rr_int <= 714) & (rr_int >= 666)) - ecg_analysis.Tloc((rr_int <= 714) & (rr_int >= 666)));
% 
% % Col 50: median Tu30Tp distance
% median_Tu30Tp_fast = median(ecg_analysis.Tloc((rr_int <= 714) & (rr_int >= 666)) - ecg_analysis.TPars.prepeak_t30_loc((rr_int <= 714) & (rr_int >= 666)));
% % Col 51: median TpTd30
% median_TpTd30_fast = median(ecg_analysis.TPars.postpeak_t30_loc((rr_int <= 714) & (rr_int >= 666)) - ecg_analysis.Tloc((rr_int <= 714) & (rr_int >= 666)));
% % Col 52: median Tu50Tp distance
% median_Tu50Tp_fast = median(ecg_analysis.Tloc((rr_int <= 714) & (rr_int >= 666)) - ecg_analysis.TPars.prepeak_t50_loc((rr_int <= 714) & (rr_int >= 666)));
% % Col 53: median TpTd50
% median_TpTd50_fast = median(ecg_analysis.TPars.postpeak_t50_loc((rr_int <= 714) & (rr_int >= 666)) - ecg_analysis.Tloc((rr_int <= 714) & (rr_int >= 666)));
% % Col 54: median Tu70Tp distance
% median_Tu70Tp_fast = median(ecg_analysis.Tloc((rr_int <= 714) & (rr_int >= 666)) - ecg_analysis.TPars.prepeak_t70_loc((rr_int <= 714) & (rr_int >= 666)));
% % Col 55: median TpTd70
% median_TpTd70_fast = median(ecg_analysis.TPars.postpeak_t70_loc((rr_int <= 714) & (rr_int >= 666)) - ecg_analysis.Tloc((rr_int <= 714) & (rr_int >= 666)));
% 
% % Cols 56-61 TampRamps
% mean_TampRamp = mean(in_TampRamp);
% median_TampRamp = median(in_TampRamp);
% mean_TampRamp_slow = mean(in_TampRamp((inliers.RR_TampRamp(:,1) <=1000) & (inliers.RR_TampRamp(:,1) >=909)));
% median_TampRamp_slow = median(in_TampRamp((inliers.RR_TampRamp(:,1) <=1000) & (inliers.RR_TampRamp(:,1) >=909)));
% mean_TampRamp_fast = mean(in_TampRamp((inliers.RR_TampRamp(:,1) <=714) & (inliers.RR_TampRamp(:,1) >=666)));
% median_TampRamp_fast = median(in_TampRamp((inliers.RR_TampRamp(:,1) <=714) & (inliers.RR_TampRamp(:,1) >=666)));
% 
% % Cols 62-67: TampRamp(n+1) - TampRamp(n)
% mean_TampRamp_difference = mean(diff(in_TampRamp));
% median_TampRamp_difference = median(diff(in_TampRamp));
% mean_TampRamp_difference_slow = mean(diff(in_TampRamp((inliers.RR_TampRamp(:,1) <=1000) & (inliers.RR_TampRamp(:,1) >=909))));
% median_TampRamp_difference_slow = median(diff(in_TampRamp((inliers.RR_TampRamp(:,1) <=1000) & (inliers.RR_TampRamp(:,1) >=909))));
% mean_TampRamp_difference_fast = mean(diff(in_TampRamp((inliers.RR_TampRamp(:,1) <=714) & (inliers.RR_TampRamp(:,1) >=666))));
% median_TampRamp_difference_fast = median(diff(in_TampRamp((inliers.RR_TampRamp(:,1) <=714) & (inliers.RR_TampRamp(:,1) >=666))));
% 
% % Cols 68-73: Inlier Only TpTd70s
% mean_inlier_TpTd70 = mean(in_TpTd70);
% median_inlier_TpTd70 = median(in_TpTd70);
% mean_inlier_TpTd70_slow = mean(in_TpTd70(((inliers.secRR_TpTd70(:,1) .* 1000) <=1000) & ((inliers.secRR_TpTd70(:,1) .* 1000) >=909)));
% median_inlier_TpTd70_slow = median(in_TpTd70(((inliers.secRR_TpTd70(:,1) .* 1000) <=1000) & ((inliers.secRR_TpTd70(:,1) .* 1000) >=909)));
% mean_inlier_TpTd70_fast = mean(in_TpTd70(((inliers.secRR_TpTd70(:,1) .* 1000) <=714) & ((inliers.secRR_TpTd70(:,1) .* 1000) >=666)));
% median_inlier_TpTd70_fast = median(in_TpTd70(((inliers.secRR_TpTd70(:,1) .* 1000) <=714) & ((inliers.secRR_TpTd70(:,1) .* 1000) >=666)));
% 
% % Cols 74-79: RR vs TpTd70 analyses 
% a_from_RRxTpTd70 = RRxTpTd70.fiteq.a;
% b_from_RRxTpTd70 = RRxTpTd70.fiteq.b;
% pval_for_a_RRxTpTd70 = RRxTpTd70.pvals(1);
% pval_for_b_RRxTpTd70 = RRxTpTd70.pvals(2);
% rmse_from_RRxTpTd70 = RRxTpTd70.gof.rmse;
% rsquare_from_RRxTpTd70 = RRxTpTd70.gof.rsquare;
% 
% % Cols 80-83: Extra RR vs QT missed earlier.
% pval_for_a_RRxQT = RRxQT.pvals(1);
% pval_for_b_RRxQT = RRxQT.pvals(2);
% rmse_from_RRxQT = RRxQT.gof.rmse;
% rsquare_from_RRxQT = RRxQT.gof.rsquare;
% 
% % Col 84: 
% RRxTampRamp_perimeter;
% % Col 85:
% RRxQTp_perimeter;
% % Col 86:
% QTpxTampRamp_perimeter;
% % Col 87:
% RRxTpTd70_perimeter;
% 
% % Col 88-90: RRxTampRamp linear regression
% RRxTampRamp_m;
% RRxTampRamp_b;
% RRxTampRamp_rsquare;
% % Col 91-93: QTxTampRamp linear regression
% QTpxTampRamp_m;
% QTpxTampRamp_b;
% QTpxTampRamp_rsquare;
% 
% 
% % append all the above onto a new row
% current_row = [recording, meanQTpeak, medianQTpeak, a_from_RRxQT, b_from_RRxQT, mean_Tu30Tp, mean_TpTd30, mean_Tu50Tp, mean_TpTd50, mean_Tu70Tp, mean_TpTd70, median_Tu30Tp, median_TpTd30, median_Tu50Tp, median_TpTd50, median_Tu70Tp, median_TpTd70, meanTamp, medianTamp, mean_tamp_difference, median_tamp_difference, meanRR, medianRR, meanRamp, medianRamp, QC_included, QC_excluded, meanQTpeak_slow, medianQTpeak_slow, meanQTpeak_fast, medianQTpeak_fast, mean_Tu30Tp_slow, mean_TpTd30_slow, mean_Tu50Tp_slow, mean_TpTd50_slow, mean_Tu70Tp_slow, mean_TpTd70_slow, median_Tu30Tp_slow, median_TpTd30_slow, median_Tu50Tp_slow, median_TpTd50_slow, median_Tu70Tp_slow, median_TpTd70_slow, mean_Tu30Tp_fast, mean_TpTd30_fast, mean_Tu50Tp_fast, mean_TpTd50_fast, mean_Tu70Tp_fast, mean_TpTd70_fast, median_Tu30Tp_fast, median_TpTd30_fast, median_Tu50Tp_fast, median_TpTd50_fast, median_Tu70Tp_fast, median_TpTd70_fast, mean_TampRamp, median_TampRamp, mean_TampRamp_slow, median_TampRamp_slow, mean_TampRamp_fast, median_TampRamp_fast, mean_TampRamp_difference, median_TampRamp_difference, mean_TampRamp_difference_slow, median_TampRamp_difference_slow, mean_TampRamp_difference_fast, median_TampRamp_difference_fast, mean_inlier_TpTd70, median_inlier_TpTd70, mean_inlier_TpTd70_slow, median_inlier_TpTd70_slow, mean_inlier_TpTd70_fast, median_inlier_TpTd70_fast, a_from_RRxTpTd70, b_from_RRxTpTd70, pval_for_a_RRxTpTd70, pval_for_b_RRxTpTd70, rmse_from_RRxTpTd70, rsquare_from_RRxTpTd70, pval_for_a_RRxQT, pval_for_b_RRxQT, rmse_from_RRxQT, rsquare_from_RRxQT, RRxTampRamp_perimeter, RRxQTp_perimeter, QTpxTampRamp_perimeter, RRxTpTd70_perimeter, RRxTampRamp_m, RRxTampRamp_b, RRxTampRamp_rsquare, QTpxTampRamp_m, QTpxTampRamp_b, QTpxTampRamp_rsquare];
% ecg_analysis.summary = current_row;
% 
% % open the repol_vars variable from file
% load(strcat(out_folder,'\repol_savers.mat'), "repol_vars");
% % Append the row onto the file.
% repol_vars = [repol_vars; current_row];
% % Save the .mat
% save(strcat(out_folder,'\repol_savers.mat'), "repol_vars");
% % clear to prevent confusion with next file.
% clear repol_vars current_row


