function [ ecg ecg_inbins ecg_time] = plot_binaveg_process(name,ecg_data,leadn,out_folder )
%PLOT_BINAVEG_PROCESS showing the beat averaging process 
%   Detailed explanation goes here
    remove_patch = 0;
    freq_list=[0.9,1.0,1.1,1.2,1.3,1.4];
    freq_tol =0.025; 
    freq     =1;
    sp_n     =5;
    sp_m     =3;
    num_bins =100;
  
  
figure;%subplot(sp_n,sp_m,[1,2,3]);
    plot(ecg_data.raw_ecgSig(ecg_data.Rloc(6000):ecg_data.Rloc(6010)));
    grid minor;
    
    if isempty(ecg_data.Rloc)
        disp('error occured in get_ramp_freq: rr_freq is empty');
        return
    end 
    
    selected_R = R_amplitude( (rr_freq >= freq-freq_tol) & (rr_freq <=freq+freq_tol) );
    std_r      = std(selected_R);
    amp_tol=std_r;
    amp        = mean(selected_R);
    
    
    Y = selected_R;
    n=hist(Y,num_bins);
    figure;%subplot(sp_n,sp_m,6);
    gbox_start  = mean(Y)-std(Y)/2;
    gbox_end    = mean(Y)+std(Y)/2;
    if (remove_patch ==0)
        ph=patch([gbox_start gbox_end gbox_end gbox_start gbox_start],...
                 [.1 .1 1.1*max(n) 1.1*max(n) .1],...
                 zeros(1,5));
        set(ph,'facecolor',.65*[1 1 1]);
        set(ph,'edgecolor',.65*[1 1 1]);
        set(ph,'FaceAlpha',0.5);
    end    
    hold on;hist(Y,num_bins);
    ylim([0,1.1*max(n)]);
    xlabel('R_{amp}');
    ylabel('Number of items');
    title(['R_{amp} freq(RR_n)=' num2str(freq) '+/-'  num2str(freq_tol) 'Hz']);
    
      
    selected_RR_ind = find((rr_freq >= freq-freq_tol) & (rr_freq <=freq+freq_tol));
    RR_n1 = selected_RR_ind+1;
    RR_n1(end) = [];
    Y = ecg_data.RR(RR_n1)*(1000/ecg_data.arash_Header.Sampling_Rate);
    n=hist(Y,num_bins);

    figure;%subplot(sp_n,sp_m,12);
    gbox_start  = mean(Y)-std(Y);
    gbox_end    = mean(Y)+std(Y);
    if (remove_patch ==0)
        ph=patch([gbox_start gbox_end gbox_end gbox_start gbox_start],...
                 [.1 .1 max(n) max(n) .1],...
                 zeros(1,5));
        set(ph,'facecolor',.65*[1 1 1]);
        set(ph,'edgecolor',.65*[1 1 1]);
        set(ph,'FaceAlpha',0.5);
    end
    hold on;hist(Y,num_bins);
    title('RR_{n+1}');
    xlabel('Freq');
    ylabel('Number of items');    
    
%     [ slow.ecg slow.ecg_list slow.r_loc slow.time slow.RR_ind] = get_ecgavg_saa( ecg_data, rr_freq, R_amplitude, freq, freq_tol, amp, std_r, leadn,0 );

%     figure;%subplot(sp_n,sp_m,9);  
%     plot_3d_avgecg(slow.ecg_list);    
%     
%      [m,n] = size(slow.ecg_list);
%     avg_sig = mean(slow.ecg_list(:,floor(0.1*n):floor(0.7*n)));
%     
%     figure;    
%     for i = 1:m
%         cur_sig = slow.ecg_list(i,floor(0.1*n):floor(0.7*n));        
%         if (sum((avg_sig-cur_sig).^2)/sum(avg_sig.^2) < 0.05)
%             hold on;plot(slow.time(1:length(cur_sig)), cur_sig,'k');
%         end
%     end
%     hold on;plot( slow.time(1:length(cur_sig)),avg_sig,'r');
%     title(ecg_data.outfnamestr);
    
%     figure;%subplot(sp_n,sp_m,15);  
%     plot( slow.time,slow.ecg,'b');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

col =varycolor(length(freq_list));    icnt =1; Nbeats_inbins=[];%ecg_collection=[];

    % VVV Iterate through frequency bins and plot seperate average beat for
    % each

    for i = 1:length(freq_list)
       
       current_freq = freq_list(i);
       selected_beats = (rr_freq >= current_freq-freq_tol) & (rr_freq <=current_freq+freq_tol);%& (R_amplitude >= amp-std_r) & (R_amplitude <= amp+std_r);

       selected_R = R_amplitude( selected_beats ); %figure,plot(selected_R);end
       %selected_R_sizeratio = (length(selected_R)/num_beats)
       
       Nbeats_inbins{i}= length((selected_R)) ;
       
       
       if (length(selected_R) > 100) %&& (std(selected_R) < 60)
           % try 
               mean_r_amp = mean(selected_R);
               ecg.freq = current_freq;
               [ ecg.rawdata list r_loc ecg.time] = get_ecgavg_saa( ecg_data, rr_freq, R_amplitude, current_freq, freq_tol, mean_r_amp, std(selected_R)/2, leadn,0 );
               ecg.ecg = ecg.rawdata;
               [~, maxindex] = max(ecg.rawdata);
               if (~isempty(ecg.rawdata( ecg.time > ecg.time(maxindex)+150 )))
                   [ecg.sigmoid estimates t_wave sw ac deac] = fit_t_sigmoid_saa(ecg.time,ecg.rawdata,ecg_data.arash_Header.Sampling_Rate,ecg.freq);
%              figure,plot(ac,'Linewidth',2);hold on;plot(deac,'r','Linewidth',2);hold on;plot(sw,'g','Linewidth',2);hold on;plot(t_wave.y,'k','Linewidth',2);hold off
                      
                   %%Polynomial curve fitting %%%%%%%%%%%%%%%%%%%%%%%%%
%         options.samplerate      = 200;
%         options.df              = 10;
%         options.downward_offset = 2;
%         [ t_poly_data ] = fit_t_bifid_polynomial_saa(xdata, ydata, options );
   ecg.sigmoid.t_wave=t_wave;ecg.sigmoid.sw=sw;ecg.sigmoid.ac=ac;ecg.sigmoid.deac=deac;
                   sigmoid_inbins{i}= ecg.sigmoid ;
                   ecg_inbins{i}=ecg.rawdata;
                   ecg_time{i}=ecg.time;
                   Allbeats_inbins{i}= list ;
                   estimates_inbins{i}=estimates;
                   %fit_t_sigmoid(ecg.time,ecg.rawdata,ecg_data.arash_Header.Sampling_Rate,ecg.freq);
                    if (ecg.sigmoid.fit_error < 0.2)
                       ecg_collection(icnt) = ecg;
                     icnt = icnt+1;
%                        if (remove_patch ==0)
%                            gbox_start  = freq_list(i)-freq_tol;
%                            gbox_end    = freq_list(i)+freq_tol;    
%                            ph=patch([gbox_start gbox_end gbox_end gbox_start gbox_start],...
%                                      [.1 .1 max(R_amplitude) max(R_amplitude) .1],...
%                                      zeros(1,5));
%                            set(ph,'facecolor',col(i,:));
%                            set(ph,'edgecolor',col(i,:));
%                            set(ph,'FaceAlpha',0.3);                              
%                        end
                    end    
               end
           % catch
           %    disp(['error occured while analysing sample at ' num2str(current_freq)]);
           else
                   sigmoid_inbins{i}= [] ;
                   ecg_inbins{i}=[];
                   ecg_time{i}=[];
                   Allbeats_inbins{i}=[];
                   estimates_inbins{i}=[];
    end
end
    
    ecg.sigmoids_inbins=sigmoid_inbins;
    ecg.Nbeats_inbins=Nbeats_inbins;
    ecg.Allbeats_inbins=Allbeats_inbins;
    ecg.estimates_inbins =estimates_inbins ;
    
    ecg.ecg_inbins=ecg_inbins;
    ecg.ecg_time=ecg_time;
     
%      close all
%    if size(ecg_collection)~=0  
    figure;
%     
%     set(gcf,'PaperType','A4')
%     set(gcf,'PaperOrientation','landscape')
%     set(gcf,'PaperPositionMode','manual')
%     set(gcf,'PaperUnits','centimeters')
%     set(gcf,'PaperPosition',[0,0,28,20])
    
       
 
    const_max_time  = 1000;
    
%     ecg_inbins=[];
    for i=1: icnt-1
        max_time = max(ecg_collection(i).time);
        if (max_time > const_max_time)
            t_end = const_max_time;
        else
            t_end = max_time;
        end

        selected_time = ecg_collection(i).time <  t_end;
        selected_time1 = ecg_collection(i).time <  t_end+1000; %+200
        hold on;
%         ecg_inbins{i}=ecg_collection(i).rawdata(selected_time1);
%         ecg_time{i}=ecg_collection(i).time(selected_time1);
        plt_h(i) = plot( ecg_collection(i).time(selected_time), ecg_collection(i).rawdata(selected_time),'color',col(i,:));
        title(strcat('Average curves-control-night- ',name));
        plt_legend{i} = ['freq:' num2str(ecg_collection(i).freq)];            
 
    end   

    if exist('plt_h') == 1 % GP: Adding existence check for plt_h
        legend(plt_h,plt_legend); 
    end
%   tachoname=strcat('lqt2_24_',name);  print ('-dpsc2', '-append', tachoname)  %subplot(2,1,2);
%    end

end

