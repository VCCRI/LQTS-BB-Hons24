% Prompts an input folder of input

%% PARAMETERS
diff_map.enter_mapping = 0; % Enter the difference mapping mode. 0 = no, 1 = yes
diff_map.par = 'RR_QTp'; % 'RR_NRamp' % 'RR_Ramp' % 'RR_RampTamp' % 'Tamp_Ramp' % RR_TampRamp % 'RR_QTp' % 'QTp_Tamp' % 'QTp_TampRamp' % 'RR_TpTd70'
diff_map.minuend = 2;% Choose the iteration number that you will subtract from. 1 is first recording, 2 is second recording, 3 is third recording, etc. Default = 2 
diff_map.subtrahend = 1; % Choose the iteration number that will be subtracted from the minuend. 1 means first recording, i.e. minuend - 1st. 1 is 1st recording, 2 is 2nd recording, etc.

% END Parameters

inp_folder = uigetdir('.','Please select input folder...');

out_folder = strcat(inp_folder,'\ProcessedFiles');
mkdir(out_folder);

dirList = dir(strcat(inp_folder,'\*.mat'));
names_list = cell(length(dirList));
errorfiles = [];


% Initiate table for holter parameters saving parametersl
parameters = table();

% preparing difference heatmaps struc - will populate with F from each iteration of smoothhist2d (saved in ecg_analysis.Feiler)to become a multimensional array.
all_feilers = struct();
all_feilers.RR_NRamp = [];
all_feilers.RR_Ramp = [];
all_feilers.RR_RampTamp = [];
all_feilers.Tamp_Ramp = [];
all_feilers.RR_TampRamp = [];
all_feilers.RR_QTp = [];
all_feilers.QTp_Tamp = [];
all_feilers.QTp_TampRamp = [];
all_feilers.RR_TpTd70 = [];



% Prepare for the variable saving (for my figures) % If commented, Need to
% also comment out saving repol section in LQTS_Program_Function
% repol_vars = zeros(1, 93);
% save(strcat(out_folder,'\repol_savers.mat'),"repol_vars");

for i = 1:length(dirList)
    tic
    [~,name,] = fileparts(dirList(i).name);
    named_out_folder = strcat(inp_folder,'\ProcessedFiles\',name);
    mkdir(named_out_folder);
    file = strcat(dirList(i).folder, filesep, dirList(i).name);
%     try
    [ecg_analysis] = LQTS_Program_Function(file, name, named_out_folder, out_folder); %,out_folder);
    
    % all_feilers.RR_NRamp(:,:,i) = ecg_analysis.Feiler.RR_NRamp;
    % all_feilers.RR_Ramp(:,:,i) = ecg_analysis.Feiler.RR_Ramp;
    % all_feilers.RR_RampTamp(:,:,i) = ecg_analysis.Feiler.RR_RampTamp;
    % all_feilers.Tamp_Ramp(:,:,i) = ecg_analysis.Feiler.Tamp_Ramp;
    all_feilers.RR_TampRamp(:,:,i) = ecg_analysis.Feiler.RR_TampRamp;
    all_feilers.RR_QTp(:,:,i) = ecg_analysis.Feiler.RR_QTp;
    % all_feilers.QTp_Tamp(:,:,i) = ecg_analysis.Feiler.QTp_Tamp;
    all_feilers.QTp_TampRamp(:,:,i) = ecg_analysis.Feiler.QTp_TampRamp;
    all_feilers.RR_TpTd70(:,:,i) = ecg_analysis.Feiler.RR_TpTd70;

    summary_data(:,:,i) = ecg_analysis.summary;

    plot_data.RRxTampRamp{i} = ecg_analysis.PD.RRxTampRamp;
    plot_data.RRxQTp{i} = ecg_analysis.PD.RRxQTp;
    plot_data.QTpxTampRamp{i} = ecg_analysis.PD.QTpxTampRamp;
    plot_data.RRxTpTd70{i} = ecg_analysis.PD.RRxTpTd70;

    save(strcat(named_out_folder,filesep,name,'.mat'),'ecg_analysis');
    toc
%     catch
%         disp(['Unable to analyse file number ', name])
%         errorfiles = [errorfiles; {name}];
%     end
    close all
end

%% Feiler analyses for difference maps and TODO: statistical analysis.
if diff_map.enter_mapping == 1
    difference_maps(diff_map.minuend, diff_map.subtrahend, dirList(diff_map.minuend).name, dirList(diff_map.subtrahend).name, all_feilers, diff_map.par, 256);
    savetitle = "difference map " + dirList(diff_map.minuend).name(1:end-4) + " - " + dirList(diff_map.subtrahend).name(1:end-4) + " " + diff_map.par;
    % % saveas(gcf,strcat(named_out_folder,filesep, savetitle));
    % % saveas(gcf,strcat(named_out_folder,filesep, savetitle + ".png"));
  close all
else
end

% extracting SSIM and displacement vectors for selected paired data
%% LQT1
% [similarity.pt_test] = SSIM_recs(2, 1, summary_data, plot_data);
% 
% save(strcat(out_folder, filesep, 'similarity.mat'), "similarity");