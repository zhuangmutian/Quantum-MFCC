clc
clear all

%%K-NN
%% USe loop to calculate
%% Use first all as training data
%% Find miu
load('T12R CH1 Features.mat');
load('T12R CH2 Features.mat');
load('T17R CH1 Features.mat');
load('T17R CH2 Features.mat');
load('T11R CH1 Features.mat');
load('T11R CH2 Features.mat');
load('T18R CH1 Features.mat');
load('T18R CH2 Features.mat');
load('T19R CH1 Features.mat');
load('T19R CH2 Features.mat');
load('TLab(2) CH1 Features.mat');
load('TLab(2) CH2 Features.mat');


% Load training data
% X = [T12RCH1Features T12RCH2Features T17RCH1Features T17RCH2Features];
t1 = [T12RCH1Features T17RCH1Features ];
[rt1,ct1] = size(t1);
miu1 = zeros(rt1,1);
for i = 1:rt1
    miu1(i) = (sum(t1(i,:)))/ct1; 
end

t2 = [T11RCH2Features TLab2RCH2Features];
[rt2,ct2] = size(t1);
miu2 = zeros(rt2,1);
for i = 1:rt2
    miu2(i) = (sum(t2(i,:)))/ct2; 
end

%% Load test data and calculate pdf
%%
test1 = [T12RCH2Features T17RCH2Features T18RCH1Features T19RCH1Features T18RCH2Features T19RCH2Features];
% test1 = [T12RCH2Features T17RCH2Features T18RCH1Features T19RCH1Features]; %T18RCH2Features T19RCH2Features];
[rtest,ctest] = size(test1);
label = ones(1,ctest);
test = test1;

% ctest = 160;

%%
% test2 = [T11RCH2Features T11RCH2Features TLab2RCH1Features TLab2RCH2Features];
% [rtest,ctest] = size(test2);
% label = 2*ones(1,ctest);
% test = test2;

% a = 0;
g = 0;
b = 0;
for number = 1:ctest
    
%     pdf = (1/((sqrt(2*pi))^13))*(1/sqrt(det(sigma)))*exp((-1/2)*((test(:,number) - miu)'*(inv(sigma))*(test(:,number) - miu)));
% pdf = (1/((sqrt(2*pi))^13))*(1/sqrt(det(sigma)))*exp((-1/2)*((test1(:,number) - miu)'*(inv(sigma))*(test1(:,number) - miu)));
    dis1 = sqrt(sum((test(:,number)-miu1).^2));
    dis2 = sqrt(sum((test(:,number)-miu2).^2));
    if dis1 <= dis2*0.3
        class = 1;
    else 
        class = 2;
    end
%     if class == label(number);
%         a = a+1;
%     end
    if class == 1;
        g = g+1;
    else
        b = b+1;
    end
end
%  p = a/ctest
pd =g/ctest

pfa = g/ctest
        


