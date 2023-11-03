function [Xd,Xnoise,NoiseLoc] = SDROM(Xn,Nw,T)
% This is a function that shows the audio processing in terms of denoising

% Implementation of SD-ROM method
% AN EFFICIENT METHOD FOR THE REMOVAL OF IMPULSE NOISE FROM SPEECH AND AUDIO SIGNALS
% Charu Chandra, Michael S. Moore and Sanjit IS. Mitra

% Input:
% audio data with impulsive noise: Xn
% window size: Nw (e.g. 5) should be odd number

% T[T1 T2 ...] thresholds according to the window length (size: Dshift)
% T1 threshold (e.g., 4)
% T2 threshold (e.g., 12)
% Every detected impulse is replaced by the ROM, m[k]=(r2[k]+r3[k])/2

% Output
% denoised Xn: Xd
% noise: Xnoise
% locations of impulsive noise: NoiseLoc (to be compared with the true
% locations "noiseIndexes_sorted" from addImplNoise.m

%%%%%%%%%%%
Xd=Xn;
N=length(Xn); % audio data length
% Check for Nw being odd number
if mod(Nw,2)==0   % Even number
    error('Error. \nWindw legnth Nw must be odd number.')
end
Dshift=(Nw-1)/2; % the shift left and right of the current sample
% thresholds should be in mumber according to the window legnth (size: Dshift)
if length(T)~=Dshift
     error('Error. \nThresholds are missing-Please provide (Nw-1)/2 trheshold values in T vector .')
end
NoiseLoc=[];
tt=0;
for i=Dshift+1:N-Dshift
    Xw=[]; wn=[]; rn=[];mu=[];dn=[];
    
    Xw=Xn(i-Dshift:i+Dshift); % windowed signal

    wn=[Xw(1:Dshift) Xw(Dshift+2:end)];

    rn=sort(wn); % ranked ascending order
    % ranked-ordered difference
    mu=(rn(length(rn)/2)+rn((length(rn)/2)+1))/2; %Ranked order mean (ROM)
    for k=1:Dshift
        if Xn(i)<=mu
            dn(k)=rn(k)-Xn(i);
        else
            dn(k)=Xn(i)-rn(Nw-1-k);
        end
    end
    ss=0;
    for k=1:Dshift
        if dn(k)>T(k)
            ss=ss+1;
        end
    end
    if ss>0
        Xd(i)=mu;
        tt=tt+1;
        NoiseLoc(tt)=i;
    end
end
Xnoise=Xn-Xd;

% plots
% % % figure(1)
% % % subplot(311)
% % % plot(Xn)
% % % xlabel('Samples (n)')
% % % ylabel('Amplitude')
% % % title('Original noisy signal')
% % % xlim([1 N])
% % % ylim([min(Xn) max(Xn)])
% % % 
% % % 
% % % subplot(312)
% % % plot(Xd)
% % % xlabel('Samples (n)')
% % % ylabel('Amplitude')
% % % title('SD-ROM denoised signal')
% % % xlim([1 N])
% % % ylim([min(Xn) max(Xn)])
% % % 
% % % 
% % % subplot(313)
% % % plot(Xnoise)
% % % xlabel('Samples (n)')
% % % ylabel('Amplitude')
% % % title('Estimated noise signal')
% % % xlim([1 N])
% % % ylim([min(Xn) max(Xn)])
