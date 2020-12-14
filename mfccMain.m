% % Read input signal for demonstration
% [x, fs] = audioread('asthma1.WAV');
% 
% % Initial parameters

% Main program to get MFCC

clc 
clear all
close all

%change name for the ducument here
InPutVoice = 'pneumonia3.wav';
OutPutName = 'pneumonia_MFCC_1';

%read the audio, using 2nd line.
[Breath_orig,Fs] = audioread(InPutVoice);

Breath=Breath_orig(:,2)';


% [startpoint, endpoint, sound]= voicedetect1(Breath);
[sound1]= voicedetect1(Breath);
figure
plot(sound1);
% sound1 = sound1(1:1500);
% Desired number of Mel coefficients
numCoeff = 13;

% Number of mel banks
numMel = 40;

% Number of FFT points
numFFT = 1024;

% Create mel matrix for demonstration
signalScaled = freqToMel(3000, numFFT, Fs) * abs(spectrogram(sound1,numFFT));

% Compute MFCC matrix 
EVmfccMat = myEVMFCC(numCoeff, numMel, numFFT, sound1, Fs);
mfccMat = myMFCC(numCoeff, numMel, numFFT, sound1, Fs);
TVmfccMat = myTVMFCC(numCoeff, numMel, numFFT, sound1, Fs);
% % Display Frequencies on mel scale
% subplot(2,1,1);
% imagesc(signalScaled);
% 
% % Display MFCC Spectrogram
% subplot(2,1,2);
% imagesc(mfccMat);
pneumonia3mfccMat = mfccMat';
pneumonia3EVmfccMat = EVmfccMat';
pneumonia3TVmfccMat = TVmfccMat';
