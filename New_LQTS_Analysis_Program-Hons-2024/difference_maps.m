function [outputArg1,outputArg2] = difference_maps(minuend, subtrahend, minuend_name, subtrahend_name, all_feilers, par, m)
%DIFFERENCE_MAPS Summary of this function goes here

% USER MUST CHOOSE WHICH HEATMAPS ARE PRE AND POST/WHICH heatmaps to
% compare.
% Variable 'par' must be one of the following: RR_NRamp' % 'RR_Ramp' %
% 'RR_RampTamp' % 'Tamp_Ramp' % RR_TampRamp % 'RR_QTp' % 'QTp_Tamp' % 'QTp_TampRamp' % 'RR_TpTd70'
% all_feilers should be a struc with the different par options inside it.
% Each par of all_feilers struc should be a (:,:,n) matrix with n being
% the number of recordings looped through.


%% French flag
modif = fix(m/2);

r = [ones(modif,1); (sqrt(1-((1:modif)/modif).^2))'];
g = [flipud( (sqrt(1-((1:modif)/modif).^2))'); (sqrt(1-((1:modif)/modif).^2))'];
b = [flipud((sqrt(1-((1:modif)/modif).^2))'); ones(modif,1)];
french_map = [r g b];



%% Choose the parameter to be plotted.
if par == "RR_NRamp"
    subtracted_plot_data = all_feilers.RR_NRamp(:,:,minuend) - all_feilers.RR_NRamp(:,:,subtrahend);
    xlab = "RR interval (ms)";
    ylab = "Normalised R amplitude (to median)";
    xticknumber = 9;
    yticknumber = 11;
    xtlab = linspace(400, 2000, xticknumber);
    ytlab = linspace(-10,10, yticknumber);
    xlinsp_tick_max = 320; % x dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.
    ylinsp_tick_max = 200; % y dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.

elseif par == "RR_Ramp"
    subtracted_plot_data = all_feilers.RR_Ramp(:,:,minuend) - all_feilers.RR_Ramp(:,:,subtrahend);
    xlab = "RR interval (ms)";
    ylab = "Filtered R amplitude";
    xticknumber = 9;
    yticknumber = 9;
    xtlab = linspace(400, 2000, xticknumber);
    ytlab = linspace(-50,350, yticknumber);
    xlinsp_tick_max = 320; % x dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.
    ylinsp_tick_max = 350; % y dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.

elseif par == "RR_RampTamp"
    subtracted_plot_data = all_feilers.RR_RampTamp(:,:,minuend) - all_feilers.RR_RampTamp(:,:,subtrahend);
    xlab = "RR interval (ms)";
    ylab = "Ramp/Tamp";
    xticknumber = 9;
    yticknumber = 11;
    xtlab = linspace(400, 2000, xticknumber);
    ytlab = linspace(-25,25, yticknumber);
    xlinsp_tick_max = 320; % x dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.
    ylinsp_tick_max = 500; % y dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.

elseif par == "Tamp_Ramp"
    subtracted_plot_data = all_feilers.Tamp_Ramp(:,:,minuend) - all_feilers.Tamp_Ramp(:,:,subtrahend);
    xlab = "T amplitude (units)";
    ylab = "R amplitude (units)";
    xticknumber = 7;
    yticknumber = 11;
    xtlab = linspace(-100, 200, xticknumber);
    ytlab = linspace(-50,350, yticknumber);
    xlinsp_tick_max = 300; % x dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.
    ylinsp_tick_max = 400; % y dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.

elseif par == "RR_TampRamp"
    subtracted_plot_data = all_feilers.RR_TampRamp(:,:,minuend) - all_feilers.RR_TampRamp(:,:,subtrahend);
    xlab = "RR interval (ms)";
    ylab = "Tamp/Ramp";
    xticknumber = 9;
    yticknumber = 11;
    xtlab = linspace(400, 2000, xticknumber);
    ytlab = linspace(-1,1.5, yticknumber);
    xlinsp_tick_max = 320; % x dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.
    ylinsp_tick_max = 150; % y dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.

elseif par == "RR_QTp"
    subtracted_plot_data = all_feilers.RR_QTp(:,:,minuend) - all_feilers.RR_QTp(:,:,subtrahend);
    xlab = "RR interval (ms)";
    ylab = "QTpeak interval";
    xticknumber = 9;
    yticknumber = 12;
    xtlab = linspace(400, 2000, xticknumber);
    ytlab = linspace(200, 750, yticknumber);
    xlinsp_tick_max = 320; % x dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.
    ylinsp_tick_max = 110; % y dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.

elseif par == "QTp_Tamp"
    subtracted_plot_data = all_feilers.QTp_Tamp(:,:,minuend) - all_feilers.QTp_Tamp(:,:,subtrahend);
    xlab = "QTpeak interval";
    ylab = "T amplitude (units)";
    xticknumber = 12;
    yticknumber = 11;
    xtlab = linspace(200, 750, xticknumber);
    ytlab = linspace(-50,200, yticknumber);
    xlinsp_tick_max = 110; % x dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.
    ylinsp_tick_max = 250; % y dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.

elseif par == "QTp_TampRamp"
    subtracted_plot_data = all_feilers.QTp_TampRamp(:,:,minuend) - all_feilers.QTp_TampRamp(:,:,subtrahend);
    xlab = "QTpeak interval";
    ylab = "Tamp/Ramp";
    xticknumber = 12;
    yticknumber = 11;
    xtlab = linspace(200, 750, xticknumber);
    ytlab = linspace(-10,10, yticknumber);
    xlinsp_tick_max = 110; % x dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.
    ylinsp_tick_max = 150; % y dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.

elseif par == "RR_TpTd70"
    subtracted_plot_data = all_feilers.RR_TpTd70(:,:,minuend) - all_feilers.RR_TpTd70(:,:,subtrahend);
    xlab = "RR interval (ms)";
    ylab = "Tpeak to Td70 interval";
    xticknumber = 9;
    yticknumber = 12;
    xtlab = linspace(400, 2000, xticknumber);
    ytlab = linspace(0, 150, yticknumber);
    xlinsp_tick_max = 320; % x dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.
    ylinsp_tick_max = (150 - 0); % y dimensions of F from the nbins variable in smoothhist2d function. This has to change if you change the nbins for future calls.

else
    input("Please ensure your par input is one of the available options")
end




figure
imagesc(subtracted_plot_data);
drawnow
xlabel(xlab);
ylabel(ylab);
% if xlabel starts with RR, then reverse X axis to be in descending order
if startsWith(xlab, "RR") == 1
    set(gca, 'XDir', 'reverse');
end
set(gca, "YDir", "normal");

xticks(linspace(1, xlinsp_tick_max, xticknumber));
xticklabels(xtlab);
yticks(linspace(1, ylinsp_tick_max, yticknumber)); 
yticklabels(ytlab);
heading = minuend_name + " - " + subtrahend_name;
title(heading);

% Change the colorbar labels to have 0 in the centre (white) and abs(max) as both -ve and +ve lims. Relabel so min = -1 and max = 1.
% cmap = turbo(256);
% white = [1 1 1];
% white_start = -(0.1 .* max(abs(subtracted_plot_data),[],'all'));
% white_end = (0.1 .* max(abs(subtracted_plot_data),[],'all'));
% cmap(white_start:white_end, :) = repmat(white, white_end - white_start + 1, 1);
% colormap(cmap)


cbh = colorbar;
clim([-max(abs(subtracted_plot_data),[],'all'), max(abs(subtracted_plot_data), [], 'all')]);
cbh.Ticks = linspace(-max(abs(subtracted_plot_data),[],'all'), max(abs(subtracted_plot_data), [], 'all'), 9);
cbh.TickLabels = num2cell(-1:0.25:1);
cbh.Title.String = "Portion of maximum intensity";
colormap(redblue)
% saveas(%your save location)

end