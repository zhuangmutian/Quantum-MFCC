clc
clear all

%%Classification

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
X = [T12RCH1Features T17RCH1Features ];
[rtraining,ctraining] = size(X);

miu = zeros(rtraining,1);



for i = 1:rtraining
    miu(i) = (sum(X(i,:)))/ctraining;
   
end

%% Find Sigma
sigma = (X-miu)*(X-miu)'
sigma = diag(diag(sigma))
%% Load test data and calculate pdf
%%
% test1 = [T12RCH2Features T17RCH2Features T18RCH1Features T19RCH1Features T18RCH2Features T19RCH2Features];
% [rtest,ctest] = size(test1);
% label = ones(1,ctest);
% test = test1;
% ctest = 160;
%%
test2 = [T11RCH2Features T11RCH2Features TLab2RCH1Features TLab2RCH2Features];
[rtest,ctest] = size(test2);
label = 2*ones(1,ctest);
test = test2;

%%
g = 0;
b = 0;
for number = 1:ctest
    pdf = (1/((sqrt(2*pi))^13))*(1/sqrt(det(sigma)))*exp((-1/2)*((test(:,number) - miu)'*(inv(sigma))*(test(:,number) - miu)));
% pdf = (1/((sqrt(2*pi))^13))*(1/sqrt(det(sigma)))*exp((-1/2)*((test1(:,number) - miu)'*(inv(sigma))*(test1(:,number) - miu)));
    grade = pdf*1e125
    if grade > 2
        class = 1;
    else 
        class = 2;
    end
    if class == 1;
        g = g+1;
    else
        b = b+1;
    end
end

pd =g/ctest

pfa = (ctest-b)/ctest
        


