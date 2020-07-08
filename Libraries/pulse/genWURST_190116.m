function [ time, amp, phase_rad ] = genWURST_190116( BWpulse, tpulse, Npts, Norder, fnappend, fndir, flagplot )     
    % generate WURST phase_rad waveforms for Halbach scanner
    % October, 2014
    % 
    % modified: patmcd@mit.edu 18/04/08
    % 
    %  Inputs:
    %   BWpulse     :           Total bandwidth of pulse [Hz]
    %   tpulse      :           Duration of pulse [sec]
    %   Npts        :           Number of time points for pulse
    %   Norder      :           Smoothing parameter for pulse (ie:
    %                           amplitude envelope = (1 - ||cos(a*t)||)^Norder )
    %   fnappend    :           base string for pulse filenames
    %   fndir       :           Directory at which to save files (include
    %                           '/' at end)
    %   flagplot    :           if 0 -> don't plot; else -> generate plots
    %
    %  Outputs:
    %   time        :           vector of time points [sec]
    %   amp         :           vector of pulse amplitude [norm. to unity]
    %   phase       :           vector of pulse phase [rad]
    %
    %
    
    if nargin < 7
        flagplot = 0;
    end
    if nargin < 6
        fndir = './';
    end
    if nargin < 5
        fnappend = '';
    end
    if nargin < 4 || isempty(Norder)
        Norder = 40;
    end
    
    sweep_range = BWpulse;     % sweep range
    tau_pulse = tpulse; % length of pulse in seconds
    N_order = Norder;     % how much to round off edges of boxcar amplitude
    freq_offset = 0;  % not using frequency offset

    np = Npts;   % number of points


    time = linspace(0,tau_pulse,np);


    amp = (1 - abs(cos(pi*time./tau_pulse)).^N_order);

    freq = -(freq_offset + sweep_range/2 - (sweep_range/tau_pulse).*time);

    phase = -(2*pi*((freq_offset + sweep_range/2).*time - (sweep_range./(2.*tau_pulse)).*time.^2));

    phase = phase.'*180/pi;
    phase_rad = phase*pi/180;

%     figure(10),subplot(3,1,1),plot(1000.*time,amp),xlabel('ms','FontSize',18),title('amplitude','FontSize',20)
%     subplot(3,1,2),plot(1000.*time,freq),xlabel('ms','FontSize',18),title('frequency','FontSize',20)
%     subplot(3,1,3),plot(1000.*time,phase),xlabel('ms','FontSize',18),title('phase','FontSize',20)


    out = [amp.' freq.' phase ones(numel(phase),1)];
fid = fopen([fndir 'WURST',num2str_nd(N_order), ...
             '_f',num2str_nd(sweep_range/1000, '%03u'), ...
             'k_t',num2str_nd(1000000*tau_pulse, '%04u'), ...
             'u_N', num2str_nd( np, '%04u'), ...
             fnappend '_amp001.txt'],'w');
         
%             fprintf(fid,'%8.5f  %8.5f  %8.5f   %8.5f\n',[out(:,1) out(:,2) out(:,3) out(:,4)]);
             fprintf(fid,'%8.6f\n',[1*out(:,1)]);
%             fwrite(fid,out(:,1),'double')
            fclose(fid);
            
fid = fopen([fndir 'WURST',num2str_nd(N_order), ...
             '_f',num2str_nd(sweep_range/1000, '%03u'), ...
             'k_t',num2str_nd(1000000*tau_pulse, '%04u'), ...
             'u_N', num2str_nd( np, '%04u'), ...
             fnappend '_amp098.txt'],'w');
         
%             fprintf(fid,'%8.5f  %8.5f  %8.5f   %8.5f\n',[out(:,1) out(:,2) out(:,3) out(:,4)]);
             fprintf(fid,'%8.6f\n',[98*out(:,1)]);
%             fwrite(fid,out(:,1),'double')
            fclose(fid);
            
fid = fopen([fndir 'WURST',num2str_nd(N_order), ...
             '_f',num2str_nd(sweep_range/1000, '%03u'), ...
             'k_t',num2str_nd(1000000*tau_pulse, '%04u'), ...
             'u_N', num2str_nd( np, '%04u'), ...
             fnappend '_amp100.txt'],'w');
         
%             fprintf(fid,'%8.5f  %8.5f  %8.5f   %8.5f\n',[out(:,1) out(:,2) out(:,3) out(:,4)]);
             fprintf(fid,'%8.6f\n',[100*out(:,1)]);
%             fwrite(fid,out(:,1),'double')
            fclose(fid);
            
fid = fopen([fndir 'WURST',num2str_nd(N_order), ...
             '_f',num2str_nd(sweep_range/1000, '%03u'), ...
             'k_t',num2str_nd(1000000*tau_pulse, '%04u'), ...
             'u_N', num2str_nd( np, '%04u'), ...
             fnappend '_phs.txt'],'w');
         
%             fprintf(fid,'%8.5f  %8.5f  %8.5f   %8.5f\n',[out(:,1) out(:,2) out(:,3) out(:,4)]);
             fprintf(fid,'%8.6f\n',[out(:,3)]);
%             fwrite(fid,out(:,1),'double')
            fclose(fid);
            
fid = fopen([fndir 'WURST',num2str_nd(N_order), ...
             '_f',num2str_nd(sweep_range/1000, '%03u'), ...
             'k_t',num2str_nd(1000000*tau_pulse, '%04u'), ...
             'u_N', num2str_nd( np, '%04u'), ...
             fnappend '_freq.txt'],'w');
         
%             fprintf(fid,'%8.5f  %8.5f  %8.5f   %8.5f\n',[out(:,1) out(:,2) out(:,3) out(:,4)]);
             fprintf(fid,'%8.6f\n',[out(:,2)]);
%             fwrite(fid,out(:,1),'double')
            fclose(fid);
fid = fopen([fndir 'WURST',num2str_nd(N_order), ...
             '_f',num2str_nd(sweep_range/1000, '%03u'), ...
             'k_t',num2str_nd(1000000*tau_pulse, '%04u'), ...
             'u_N', num2str_nd( np, '%04u'), ...
             fnappend '_phsn.txt'],'w');
         
%             fprintf(fid,'%8.5f  %8.5f  %8.5f   %8.5f\n',[out(:,1) out(:,2) out(:,3) out(:,4)]);
             fprintf(fid,'%8.6f\n',[-1*out(:,3)]);
%             fwrite(fid,out(:,1),'double')
            fclose(fid);
    


%     save(['WURST',num2str_nd(N_order),'_',num2str_nd(sweep_range/1000),'_',num2str_nd((1000*tau_pulse)) fnappend],'out','-ascii')

%     tnmr_amp_1st180 = 98*amp';
%     tnmr_amp_2nd180 = 100*amp';
    if flagplot==0
        return
    end
    figure; 
    subplot( 1, 3, 1);
    plotsexy(time, amp,'color',0.8*[0 0 1]);
    subplot( 1,3,2 ); 
    plotsexy(time, phase, 'color', 0.8*[1 0 0]);
    subplot( 1,3,3 ); 
    plotsexy(time, freq, 'color', 0.8*[0 1 0]);
    
end