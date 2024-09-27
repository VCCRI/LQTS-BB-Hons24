function [sim_out] = SSIM_recs(post, pre, summary_data, plot_data)
% This function finds the SSIM between two recordings specified, as well as
% the displacement vectors.

% Summary_data is a 87 long array in the following order: recording, meanQTpeak, medianQTpeak, a_from_RRxQT, b_from_RRxQT, mean_Tu30Tp, mean_TpTd30, mean_Tu50Tp, mean_TpTd50, mean_Tu70Tp, mean_TpTd70, median_Tu30Tp, median_TpTd30, median_Tu50Tp, median_TpTd50, median_Tu70Tp, median_TpTd70, meanTamp, medianTamp, mean_tamp_difference, median_tamp_difference, meanRR, medianRR, meanRamp, medianRamp, QC_included, QC_excluded, meanQTpeak_slow, medianQTpeak_slow, meanQTpeak_fast, medianQTpeak_fast, mean_Tu30Tp_slow, mean_TpTd30_slow, mean_Tu50Tp_slow, mean_TpTd50_slow, mean_Tu70Tp_slow, mean_TpTd70_slow, median_Tu30Tp_slow, median_TpTd30_slow, median_Tu50Tp_slow, median_TpTd50_slow, median_Tu70Tp_slow, median_TpTd70_slow, mean_Tu30Tp_fast, mean_TpTd30_fast, mean_Tu50Tp_fast, mean_TpTd50_fast, mean_Tu70Tp_fast, mean_TpTd70_fast, median_Tu30Tp_fast, median_TpTd30_fast, median_Tu50Tp_fast, median_TpTd50_fast, median_Tu70Tp_fast, median_TpTd70_fast, mean_TampRamp, median_TampRamp, mean_TampRamp_slow, median_TampRamp_slow, mean_TampRamp_fast, median_TampRamp_fast, mean_TampRamp_difference, median_TampRamp_difference, mean_TampRamp_difference_slow, median_TampRamp_difference_slow, mean_TampRamp_difference_fast, median_TampRamp_difference_fast, mean_inlier_TpTd70, median_inlier_TpTd70, mean_inlier_TpTd70_slow, median_inlier_TpTd70_slow, mean_inlier_TpTd70_fast, median_inlier_TpTd70_fast, a_from_RRxTpTd70, b_from_RRxTpTd70, pval_for_a_RRxTpTd70, pval_for_b_RRxTpTd70, rmse_from_RRxTpTd70, rsquare_from_RRxTpTd70, pval_for_a_RRxQT, pval_for_b_RRxQT, rmse_from_RRxQT, rsquare_from_RRxQT, RRxTampRamp_perimeter, RRxQTp_perimeter, QTpxTampRamp_perimeter, RRxTpTd70_perimeter];
% This will not work, or will break if the order of the summary_data
% changes. You need to modify the below parameters:

RR_index = 23; % summary_data index for median RR
Tamp_index = 57;% summary_data index for  median TampRamp
QTp_index = 3; % summary index for median QTp
TpTd60_index = 69; % summary index for median TpTd70


%% RRxTampRamp

figure
smoothhist2D(plot_data.RRxTampRamp{post}, [400,-0.25],[2000,1.25],3,[(2000-400)/5, 350],.05, 'image', 0);
set(gca, "XDir", "reverse"); % Only for RR intervals. Reverses x axis to be in descending RR (ascending HR).
colormap('turbo');
min = gcf;

figure
smoothhist2D(plot_data.RRxTampRamp{pre}, [400,-0.25],[2000,1.25],3,[(2000-400)/5, 350],.05, 'image', 0);
set(gca, "XDir", "reverse"); % Only for RR intervals. Reverses x axis to be in descending RR (ascending HR).
colormap('turbo');
sub = gcf;

frame1 = getframe(min).cdata;
frame2 = getframe(sub).cdata;

gray1 = rgb2gray(frame1);
gray2 = rgb2gray(frame2);

[ssimval, ssimmap] = ssim(gray1, gray2);
imshow(ssimmap);

xlen = summary_data(:,RR_index,post) - summary_data(:,RR_index,pre);
ylen = summary_data(:,Tamp_index,post) - summary_data(:,Tamp_index,pre);

magnitude = sqrt((xlen)^2 + (ylen)^2);
direction = atan2(ylen, xlen);
direction = rad2deg(direction);

sim_out.RRxTampRamp.ssimval = ssimval;
sim_out.RRxTampRamp.ssimmap = ssimmap;
sim_out.RRxTampRamp.vector_magnitude = magnitude;
sim_out.RRxTampRamp.vector_theta = direction;

close all
clear min sub frame1 frame2 gray1 gray 2 ssimval ssimmap xlen ylen magnitude direction

%% RRxQTp

figure
smoothhist2D(plot_data.RRxQTp{post}, [0.4,200],[2,750],3,[(1000 * (2-0.4))/5, (750-200)/5], 0.05, 'image', 0);
set(gca, "XDir", "reverse"); % Only for RR intervals. Reverses x axis to be in descending RR (ascending HR).
colormap('turbo');
min = gcf;

figure
smoothhist2D(plot_data.RRxQTp{pre}, [0.4,200],[2,750],3,[(1000 * (2-0.4))/5, (750-200)/5], 0.05, 'image', 0);
set(gca, "XDir", "reverse"); % Only for RR intervals. Reverses x axis to be in descending RR (ascending HR).
colormap('turbo');
sub = gcf;

frame1 = getframe(min).cdata;
frame2 = getframe(sub).cdata;

gray1 = rgb2gray(frame1);
gray2 = rgb2gray(frame2);

[ssimval, ssimmap] = ssim(gray1, gray2);
imshow(ssimmap);

xlen = summary_data(:,RR_index,post) - summary_data(:,RR_index,pre);
ylen = summary_data(:,QTp_index,post) - summary_data(:,QTp_index,pre);

magnitude = sqrt((xlen)^2 + (ylen)^2);
direction = atan2(ylen, xlen);
direction = rad2deg(direction);

sim_out.RRxQTp.ssimval = ssimval;
sim_out.RRxQTp.ssimmap = ssimmap;
sim_out.RRxQTp.vector_magnitude = magnitude;
sim_out.RRxQTp.vector_theta = direction;

close all
clear min sub frame1 frame2 gray1 gray 2 ssimval ssimmap xlen ylen magnitude direction




%% QTpxTampRamp
figure
smoothhist2D(plot_data.QTpxTampRamp{post}, [200,-0.25],[750,1.25],3,[(750-200)/5, 150],.05, 'image', 0);
colormap('turbo');
min = gcf;

figure
smoothhist2D(plot_data.QTpxTampRamp{pre}, [200,-0.25],[750,1.25],3,[(750-200)/5, 150],.05, 'image', 0);
colormap('turbo');
sub = gcf;

frame1 = getframe(min).cdata;
frame2 = getframe(sub).cdata;

gray1 = rgb2gray(frame1);
gray2 = rgb2gray(frame2);

[ssimval, ssimmap] = ssim(gray1, gray2);
imshow(ssimmap);

xlen = summary_data(:,QTp_index,post) - summary_data(:,QTp_index,pre);
ylen = summary_data(:,Tamp_index,post) - summary_data(:,Tamp_index,pre);

magnitude = sqrt((xlen)^2 + (ylen)^2);
direction = atan2(ylen, xlen);
direction = rad2deg(direction);

sim_out.QTpxTampRamp.ssimval = ssimval;
sim_out.QTpxTampRamp.ssimmap = ssimmap;
sim_out.QTpxTampRamp.vector_magnitude = magnitude;
sim_out.QTpxTampRamp.vector_theta = direction;

close all
clear min sub frame1 frame2 gray1 gray 2 ssimval ssimmap xlen ylen magnitude direction

%% RRxTpTd70
figure
smoothhist2D(plot_data.RRxTpTd70{post}, [0.4,0],[2,150],3,[(1000 * (2-0.4))/5, 150], 0.05, 'image', 0);
set(gca, "XDir", "reverse"); % Only for RR intervals. Reverses x axis to be in descending RR (ascending HR).
colormap('turbo');
min = gcf;

figure
smoothhist2D(plot_data.RRxTpTd70{pre}, [0.4,0],[2,150],3,[(1000 * (2-0.4))/5, 150], 0.05, 'image', 0);
set(gca, "XDir", "reverse"); % Only for RR intervals. Reverses x axis to be in descending RR (ascending HR).
colormap('turbo');
sub = gcf;

frame1 = getframe(min).cdata;
frame2 = getframe(sub).cdata;

gray1 = rgb2gray(frame1);
gray2 = rgb2gray(frame2);

[ssimval, ssimmap] = ssim(gray1, gray2);
imshow(ssimmap);

xlen = summary_data(:,RR_index,post) - summary_data(:,RR_index,pre);
ylen = summary_data(:,TpTd60_index,post) - summary_data(:,TpTd60_index,pre);

magnitude = sqrt((xlen)^2 + (ylen)^2);
direction = atan2(ylen, xlen);
direction = rad2deg(direction);

sim_out.RRxTpTd70.ssimval = ssimval;
sim_out.RRxTpTd70.ssimmap = ssimmap;
sim_out.RRxTpTd70.vector_magnitude = magnitude;
sim_out.RRxTpTd70.vector_theta = direction;

close all
clear min sub frame1 frame2 gray1 gray 2 ssimval ssimmap xlen ylen magnitude direction

end

