clc
clear all
% close all
% 
[audio, fs] = audioread('TOther (3).wav');
audiocha1 = audio(:,1);
audiocha2 = audio(:,2);
% figure
% plot(audiocha1);
% % % sound(audiocha1,fs);
% % 
startpoint =  4.752e5 ;
input = audiocha1(startpoint:startpoint+floor(0.*fs));
sound(input,fs);
