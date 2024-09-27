function [handles,day,night,total] = select_n_analysis(ecg_data, leadn, day_starttime,day_endtime,night_starttime,night_endtime)

% leadn=1;
% day_starttime=8;
% day_endtime=20;
% night_starttime=1;
% night_endtime=6;

    holter_starttime = ecg_data.arash_Header.inf.Start_Time;holter_starttime1 = ecg_data.arash_Header.inf.Start_Time;
    holter_starttime = datenum(0,0,0,holter_starttime(1),holter_starttime(2),holter_starttime(3));

    day_start=time_calc(holter_starttime,day_starttime);
    day_end  =time_calc(holter_starttime,day_endtime);
    night_start=time_calc(holter_starttime,night_starttime);
    night_end  =time_calc(holter_starttime,night_endtime);
    
    handles.day_starttime  = day_starttime;
    handles.day_endtime    = day_endtime;
    handles.night_starttime= night_starttime;
    handles.night_endtime  = night_endtime;
    handles.holter_starttime = holter_starttime;
    
    length_rloc    = length(ecg_data.Rloc);
    beat_ind_times = convert_beatind_time(ecg_data,ecg_data.Rloc(1:length_rloc-100));    
    
    %%%%% DO NOT NEED 'process_signal_selection', just replace with all
    %%%%% analysis in 'get_R_addrlist_saa.m' to be rebranded as 'annotate
    %%%%% signal
    total = process_signal_selection_saa(ecg_data, leadn, 1,length_rloc-100); % Analysis of whole signal
    
%     if day_start>holter_starttime
%         daystart2=day_start;
%         day_start=holter_starttime;
% %         day   = process_signal_selection(ecg_data, leadn, find(beat_ind_times >= day_start,1),find(beat_ind_times >= day_end,1));
%         
% %     else
% %          day   = process_signal_selection(ecg_data, leadn, find(beat_ind_times >= day_start,1),find(beat_ind_times >= day_end,1));
%     end

 if day_start>holter_starttime
        day_start2=day_start;
%         day_end2=datenum(0,0,1,holter_starttime1(1),0,0);
        day_end2=beat_ind_times(end);
        day_start=holter_starttime;
    end
    
    
     day   = process_signal_selection_saa(ecg_data, leadn, find(beat_ind_times >= day_start,1),find(beat_ind_times >= day_end,1)); % Analysis of only day
    night = process_signal_selection_saa(ecg_data, leadn, find(beat_ind_times >= night_start,1),find(beat_ind_times >= night_end,1)); % Analysis of only night
% %     nextday   = process_signal_selection(ecg_data, leadn, find(beat_ind_times >= day_start2,1),find(beat_ind_times >= day_end2,1));
%     
    
    h(1)=subplot(2,2,1);
    if (night.signal_len ~= 0)
       axis_range(1,:) =plot_hist(night.rr_freq,night.R_amplitude,'night');
    end
    
    h(2)=subplot(2,2,2);
    if (day.signal_len ~= 0)
       axis_range(2,:) =plot_hist(day.rr_freq,day.R_amplitude,'day');
    end   

    h(3)=subplot(2,2,3);
    if (total.signal_len ~= 0)
       axis_range(3,:) =plot_hist(total.rr_freq,total.R_amplitude,[ecg_data.outfnamestr ' lead' num2str(leadn)]);
       %axis_range(3,:) =plot_hist(day.rr_freq,day.R_amplitude,['day']);
    end
    
    x_min = min(axis_range(:,1));
    x_max = max(axis_range(:,2));
    y_min = min(axis_range(:,3));
    y_max = max(axis_range(:,4));
    
    handles.axis_range = [x_min x_max y_min y_max];
    handles.subplot_handels = h;
    handles.ecg_data = ecg_data;
    handles.leadn = leadn;   
    handles.beat_ind_times = beat_ind_times;
    handles.day_start =day_start;
    handles.day_end =day_end;
    handles.night_start = night_start;
    handles.night_end =night_end;
    handles.total = total;
    
    %GP: TEMPORARY CHANGE FOR HEATMAP CLARITY%%%%%%%%%%%%
    y_max = 1000;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    axis(h,[x_min x_max y_min y_max]);  
    ecg_plot_h = subplot(2,2,4);
    
    handles.ecg_plot_h = ecg_plot_h;
    handles.selected_freq = overlap_freq(day,night);
    
    
    selected_data.night_freq = handles.selected_freq;
    selected_data.day_freq   = handles.selected_freq;    
    
    [day3,night3] = plot_ecg(selected_data,handles);
    sample_rate = handles.ecg_data.arash_Header.Sampling_Rate;

    %% TODO: Is there a reason for the commenting here? Why not plot day sigmoid data?
% %     
% %     if (~isempty(day3.ecg))
% %        day3.sigmoid = fit_t_sigmoid(day3.time,day3.ecg,sample_rate,selected_data.day_freq);
% %        day3.t_wave_data =plot_sigmoid_data(day3,handles,'r');       
% %     end 
% %     
    if (~isempty(night3.ecg))
       night3.sigmoid = fit_t_sigmoid(night3.time,night3.ecg,sample_rate,selected_data.night_freq);
       night3.t_wave_data =plot_sigmoid_data(night3,handles,'k');       
    end
    
%     res = [night3.t_wave_data day3.t_wave_data];
%     res = [day3.t_wave_data];
    
        
    
    
    
    
%     guidata(gcf(), handles);
%     [handles.t_wave_data,day1,night1] = fit_day_night_sigmoids([],[]);
%     
%     guidata(gcf(), handles);    
end



function [timeres] = time_calc(holter_starttime, hour)
   timeres = datenum(0,0,0,hour,0,0);
   if timeres < holter_starttime
       timeres = datenum(0,0,1,hour,0,0);
   end
end



function [ output ] = convert_beatind_time( ecg_data,beat_ind )
%CONVERT_BEATIND_TIME Summary of this function goes here
%   Detailed explanation goes here

    [days,hours,mins,secs] = convert_secs_day_frmt( beat_ind./ecg_data.arash_Header.Sampling_Rate );
    output= datenum(0,0,days,hours,mins,secs);
    
    if isfield(ecg_data.arash_Header,'inf')
        st_time = ecg_data.arash_Header.inf.Start_Time;
        start_time_offset = datenum(0,0,0,st_time(1),st_time(2),st_time(3));
        output= start_time_offset+output;
    end
    
end

function [ day hour min sec] = convert_secs_day_frmt( secs )
%CONVERT_SECS_DAY_FRMT Summary of this function goes here
%   Detailed explanation goes here
   day = floor(secs/(24*3600));
   secs = mod(secs, 24*3600);
   hour = floor(secs/3600);
   secs = mod(secs, 3600);
   min = floor(secs/60);
   sec = mod(secs, 60);
end



function [range] = plot_hist(rr_freq,R_amplitude,title_str)
    R_amplitude=R_amplitude((rr_freq >=0.25) & (rr_freq <=3.5));
    rr_freq    =rr_freq((rr_freq >=0.25) & (rr_freq <=3.5));
    
    range(1) = min(rr_freq);
    range(2) = max(rr_freq);
    range(3) = min(R_amplitude);
    range(4) = max(R_amplitude);
    
    plot_data(:,1) = rr_freq;
    plot_data(:,2) = R_amplitude;
    smoothhist2D(plot_data,5,[1000, 1000],.025);
    title(title_str,'interpreter','none');
end
