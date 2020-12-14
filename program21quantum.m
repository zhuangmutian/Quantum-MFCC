clc
clear all
% close all

[audio, fs] = audioread('TOther (7).wav');
audiocha1 = audio(:,1);
audiocha2 = audio(:,2);
% figure
% plot(audiocha1);
% sound(audiocha1,fs);
%%
[startpointlocation] =  dlmread('TOther (7) startpointlocation.txt',' ');

% T12 raw
NAME =    'TOther(7)R CH2 Features'  % change name
TOther7RCH2Features = zeros(190,40); % generate zeros matrix
 number = 1;
    startpoint = startpointlocation(number);
    
    input = audiocha2(startpoint:startpoint+floor(0.1*fs)); % which channel
    
    sound(input,fs);
    % % subplot(2,1,2);
    % % plot(audiocha2);
    input = filter([1,-0.98],[1],input);
    % sound(input,fs);

    %%
    % Duration of signal
    signallength = length(input);

    % number of frames
    % numFrames = floor(signallength/framestep);    
 
    %%

    input = input';

    numCoeff = 13;

    % Number of mel banks
    numMel = 40;

    % Number of FFT points
    numFFT = 1024;

    % Create mel matrix for demonstration
    signalScaled = freqToMel(3000, numFFT, fs) * abs(spectrogram(input,numFFT));

    % Compute MFCC matrix 
%     mfccMat = myMFCC(numCoeff, numMel, numFFT, input, fs);
% Function to compute MFCC matrix
% function mfccMatrix = myQuMFCC(numCoeff, numMel, numFFT, x, fs)  

% Arguments:
% numCoeff: desired number of MFCC coefficients
% numMel: number of filters in the mel filter bank
% numFFT: number of FFT points to be used per frame 
% x: input signal
% fs: sample rate

% input = flip(input);
x = input;
hbar = 1.054572e-34;
a = 340*(length(x))/fs;
m = 1e-65;
% m = 1;
% E = x^2;
% n = sqrt(E*2*m*a^2/(pi^2*hhat^2));

% Frame duration in seconds
windowLen = 0.025;

% Number of samples per frame
frameLen = floor(fs*windowLen);

% Duration of signal
L = length(x);

% Frame step in samples (seconds * fs)
frameStep = 0.01 * fs; 

% Maximum number of frames in signal
numFrames = floor(L/frameStep);    

% Matrix to hold cepstral coefficients
mfccMatrix = zeros(numFrames-2, numCoeff);

% Compute half the FFT length
halfFFT = 1+floor(numFFT/2); 

% Filter coefficient for high frequency accentuation
coeff = .95;

% Create mel scale matrix
mels = freqToMel(numMel, numFFT, fs);

    % Window signal, overlap frames, and compute MFCC coefficients
% for i = 1:numFrames-2
        i = 5
        

        % Frame signal
        frame = x((i-1)*frameStep+1:(i-1)*frameStep+frameLen);
        
        Phi = zeros(length(frame),1);

        % Energy
        
%         E = frame.^2;
        E = (flip(frame)).^2;
        n = sqrt(E*2*m*a^2/(pi^2*hbar^2));

        % Time
        t = ((i-1)*frameStep+1:(i-1)*frameStep+frameLen)/fs;
        
        % Position
        P = ((i-1)*frameStep:(i-1)*frameStep+frameLen-1)*340/fs;
%         P = t*340;
        P = a-P;
        % Wave Function
        
        for j = 1:length(frame)
            
            Phi(j) = sqrt(2/a)*sin((n(j)*pi*P(j)/a))*exp(-1*(1i)*(n(j)^2*pi*hbar)/(2*m*a^2)*t(j));
        end
        
        Probability = sum(abs(Phi).^2);
        
      
        % Log of energy in frame
        energy = log10(sum(1.0/frameLen*frame.^2));

        % High pass pre-emphasis filter
        for k=2:length(frame)
            frame(k) = frame(k) - coeff * frame(k-1);
        end

        % Apply window function to frame
        frame = frame.*triang(frameLen); 

        % Compute FFT of the frame
        FFT = fft(frame,frameLen);

        % Take log of non redundant mel scaled frequencies   
        feat = log10(mels.*abs(FFT(1:halfFFT)));

        % Find peaks in spectrum 
        feat = max(feat(:),1e-22);

        % Perform DCT
        c = dct(feat);

        % Replace first coefficient
        c(1) = energy;

        % Retain desired number of coefficients
        coeffs = c(1:numCoeff); 

        % Save current MFCC's in output matrix
%         mfccMatrix(i, :) = Probability*coeffs;
        mfccMatrix(i, :) = coeffs;

% end


    QumfccMat = mfccMatrix';
    mfcc = reshape(QumfccMat,[],1);
    
    TOther7RCH2QuMFCC(:,number) = mfcc;




% 
% startpoint = startpointlocation(number);
%     
% input = audiocha2(startpoint:startpoint+floor(0.1*fs)); % which channel
%     
% a = 340*(length(input))/fs;
% E = input^2;
% n = sqrt(E*2*m*a^2/(pi^2*hhat^2));
% 
% %%
% sound(input,fs);
%     % % subplot(2,1,2);
%     % % plot(audiocha2);
% input = filter([1,-0.98],[1],input);
%     % sound(input,fs);
% 
%     %%
%     % Duration of signal
% signallength = length(input);
% 
%     % number of frames
%     % numFrames = floor(signallength/framestep);    
%     
%     
%     %% Quantum
% E = input^2;
% n = sqrt(E*2*m*a^2/(pi^2*hhat^2));
% 


