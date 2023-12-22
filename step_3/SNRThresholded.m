clear;
close all;
clc;

Range_max = 200; % in meters
Range_resolution = 1; % in meters
Velocity_max = 100; % in m/s
c = 3e8; %speed of light

 
Target_range = 60; % in meters
Target_velocity = 20; % in m/s



%B = c / (2 * Range_resolution);
B = 200e6;

Tsweep = 5.5;
%Tchirp = Tsweep * (2 * Range_max / c) ;
Tchirp = 0.4e-6;
%To find slope of FMCW
slope = B / Tchirp;

fc= 24e9;             %carrier freq

Nd=256;                   % #of doppler cells OR #of sent periods % number of chirps


Nr=512;                  %for length of time OR # of range cells


t=linspace(0,Nd*Tchirp,Nr*Nd); %total time for samples


%Creating the vectors for Tx, Rx and Mix based on the total samples input.
Tx=zeros(1,length(t)); %transmitted signal
Rx=zeros(1,length(t)); %received signal
Mix = zeros(1,length(t)); %beat signal

%Similar vectors for range_covered and time delay.
r_t=zeros(1,length(t));
td=zeros(1,length(t));

for i=1:length(t)         
     
    range_at_t = 2 * (Target_range + (Target_velocity * t(i)))/c;
     
    Tx(i) = cos(2 * pi * (fc * t(i) + slope * t(i)^2 / 2));
    Rx (i)  = cos(2 * pi * (fc * (t(i) - range_at_t) + slope * (t(i) - range_at_t)^2 / 2));
    
   
    Mix(i) = Tx(i).*Rx(i); % .* is used for elementwise multiplication
    
end


% Nr = 1024, Nd = 128
Mix = reshape(Mix, [1024, 128]); 


FFT1 = abs(fft(Mix, Nr));

%To Normalize
FFT_Norm = FFT1./max(FFT1); % dividing every element by max value in the 1st FFT output. 


FFT_single_side = FFT_Norm(1:Nr/2-1);

%plotting the range
figure ('Name','Range from First FFT')

 % *%TODO* :
 % plot FFT output 
plot(FFT_single_side);
axis ([0 200 0 1]);
title('Range from First FFT');
ylabel('Normalized Amplitude');
xlabel('Range (m)');


Mix=reshape(Mix,[Nr,Nd]);

% 2D FFT using the FFT size for both dimensions.
sig_fft2 = fft2(Mix,Nr,Nd);

% Taking just one side of signal from Range dimension.
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10*log10(RDM) ;

%use the surf function to plot the output of 2DFFT and to show axis in both
%dimensions
doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure("Name", 'RDM sans thres') ,surf(doppler_axis,range_axis,RDM);
title('RDM sans thres');
xlabel('V axis');
ylabel('Range axis');
zlabel('Amplitude');

%The below line is uncommented to view the output in different angles
%view(180,0); % chnage the valeus to 45, 90, 180, ..  360 degrees 

%% CFAR implementation


Tr = 10; % Range dimension
Td = 8; % Doppler dimension


Gr = 4; % Range dimension
Gd = 4; % Doppler dimension


SNR_OFFSET = 14;

%To claculate the number of Gaurd and Training cells
% Remove CUT
Gaurd_N = (2 * Gr + 1) * (2 * Gd + 1) - 1;  
Train_N = (2 * Tr + 2 * Gr + 1) * (2 * Td + 2 * Gd + 1) - (Gaurd_N + 1);

% *%TODO* :
% This has been used in the loop
%Create a vector to store noise_level for each iteration on training cells
noise_level = zeros(1,1);


   
for range = Tr + Gr + 1 : (Nr/2) -(Tr+Gr)
    for doppler = Td + Gd + 1 : Nd - (Td+Gd)
        %Create a vector to store noise_level for each iteration on training cells
        noise_level = zeros(1,1);
        
        for i = range - (Tr + Gr) : range + Tr + Gr 
            for j = doppler - (Td + Gd) : doppler + Td + Gd 
                if (abs(range-i) > Gr || abs(range-j) > Gd)
                    noise_level = noise_level + db2pow(RDM(i,j));
                end
            end
        end
        
        threshold = SNR_OFFSET + pow2db(noise_level /(2 * (Td + Gd + 1) * 2 * (Tr + Gr + 1) - (Gr * Gd) - 1));
        
        if (RDM(range,doppler) < threshold)
            RDM(range,doppler) = 0;
        else
            RDM(range,doppler) = 1;
        end
    end
end

RDM(RDM~=0 & RDM~=1) = 0;

figure('Name', '2D CFAR on the Output of RDM'),surf(doppler_axis,range_axis,RDM);
title('RDM filtre');
xlabel('V axis');
ylabel('R axis');
zlabel('Amplitude');
colorbar;


 
