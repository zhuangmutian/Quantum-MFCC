clc
clear all
% close all

[audioT11, fs] = audioread('T37.wav');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
audiocha1 = audioT11(:,1);
audiocha2 = audioT11(:,2);

[startpointlocation] =  dlmread('T37 startpointlocation.txt',' ');%%%%%%%%%%%%%%%%%%%%%%%%%
% startpointlocation = flipud(startpointlocation);

NAME =    'T37R CH2 QuMFCC'  % %%%%%%%%%%%%%%%%%%%%%change name

T37RCH2QuMFCC = zeros(104,40); % %%%%%%%%%%%%%%%%generate zeros matrix

for number = 1:40
    startpoint = startpointlocation(number);
    
    
    input = audiocha2(startpoint:startpoint+floor(0.1*fs));%   %%%%%%%%% which channel


    input = filter([1,-0.98],[1],input);
    % Desired number of Mel coefficients
    input = input';

    numCoeff = 13;

    % Number of mel banks
    numMel = 40;

    % Number of FFT points
    numFFT = 1024;

    % Create mel matrix for demonstration
    signalScaled = freqToMel(3000, numFFT, fs) * abs(spectrogram(input,numFFT));

    % Compute MFCC matrix 
    mfccMat = myQuMFCC(numCoeff, numMel, numFFT, input, fs);
    mfccMat = mfccMat';
    mfcc = reshape(mfccMat,[],1);
    
    T37RCH2QuMFCC(:,number) = mfcc; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
save(NAME,'T37RCH2QuMFCC') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%