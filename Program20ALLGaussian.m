clc
clear all

%%Classification
load('Rdata.mat')
%% USe loop to calculate
%% Use first all as training data
%% Find miu
load('T12R CH1 Features.mat');
load('T12R CH2 Features.mat');
load('T17R CH1 Features.mat');
load('T17R CH2 Features.mat');
load('T18R CH1 Features.mat');
load('T18R CH2 Features.mat');
load('T19R CH1 Features.mat');
load('T19R CH2 Features.mat');

load('TLab(2) CH1 Features.mat');
load('TLab(2) CH2 Features.mat');
load('TLab(3) CH1 Features.mat');
load('TLab(3) CH1 Features.mat');
load('T11R CH1 Features.mat');
load('T11R CH2 Features.mat');
load('TOther(3)R CH1 Features.mat')
load('TOther(3)R CH2 Features.mat')


load('T12R CH1 MFCC.mat');
load('T12R CH2 MFCC.mat');
load('T17R CH1 MFCC.mat');
load('T17R CH2 MFCC.mat');
load('T18R CH1 MFCC.mat');
load('T18R CH2 MFCC.mat');
load('T19R CH1 MFCC.mat');
load('T19R CH2 MFCC.mat');

load('T11R CH1 MFCC.mat');
load('T11R CH2 MFCC.mat');
load('TLab(2)R CH1 MFCC.mat');
load('TLab(2)R CH2 MFCC.mat');
load('TLab(3)R CH1 MFCC.mat');
load('TLab(3)R CH2 MFCC.mat');
load('TOther(3)R CH1 MFCC.mat')
load('TOther(3)R CH2 MFCC.mat')
% Load training data
% % X = [T12RCH1Features T12RCH2Features T17RCH1Features T17RCH2Features];
XF = [T12RCH1Features T17RCH1Features T18RCH1Features T19RCH1Features];
[frtraining,fctraining] = size(XF);

fmiu = zeros(frtraining,1);



for i = 1:frtraining
    fmiu(i) = (sum(XF(i,:)))/fctraining;
   
end

% % Find Sigma
fsigma = (XF-fmiu)*(XF-fmiu)';
fsigma = diag(diag(fsigma));
%% Features  
% %% Load test data and calculate pdf
% %% Generate testing matrix and add label
% ftest1 = [T12RCH2Features T17RCH2Features T18RCH1Features T19RCH1Features T18RCH2Features T19RCH2Features];
% [frtest1,fctest1] = size(ftest1);
% flabel1 = ones(1,fctest1);
% 
% ftest2 = [T11RCH2Features T11RCH2Features TLab2RCH1Features TLab2RCH2Features TOther3RCH1Features TOther3RCH2Features];
% [frtest2,fctest2] = size(ftest2);
% flabel2 = 2*ones(1,fctest2);
% 
% ftest = [ftest1 ftest2];
% flabel = [flabel1 flabel2];
% [frtest,fctest] = size(ftest);

ftest = FeatureAllMatrix;
[frtest,fctest] = size(ftest);
fpdfmat = zeros(1,fctest);
flabel = Label;
% %%
for number = 1:fctest
    fpdf = exp((-1/2)*((ftest(:,number) - fmiu)'*(inv(fsigma))*(ftest(:,number) - fmiu)));
    fpdfmat(number) = fpdf;
end
threshold = 0:1e-6:1;
fcorrectmat = zeros(1,length(threshold));    
fPD = zeros(1,length(threshold));
fPFA = zeros(1,length(threshold));
i = 1;
for threshold =0:1e-6:1
    
    fcorrect = 0;    
    fD = 0;
    fFA = 0;

    for number = 1:fctest
%         fpdf = exp((-1/2)*((ftest(:,number) - fmiu)'*(inv(fsigma))*(ftest(:,number) - fmiu)));
%       pdf = (1/((sqrt(2*pi))^13))*(1/sqrt(det(sigma)))*exp((-1/2)*((test(:,number) - miu)'*(inv(sigma))*(test(:,number) - miu)));
%       pdf = (1/((sqrt(2*pi))^13))*(1/sqrt(det(sigma)))*exp((-1/2)*((test1(:,number) - miu)'*(inv(sigma))*(test1(:,number) - miu)));
        fpdf = fpdfmat(number);
        if fpdf > threshold;
            fclass = 1;
        end
        if fpdf <= threshold;
            fclass = 0;
        end
        if fclass == flabel(number)
            fcorrect = fcorrect + 1;
        end
        if flabel(number) == 1 && fclass == 1;
            fD = fD+1;
        end
        if flabel(number) == 0 && fclass == 1;
            fFA = fFA+1;
        end
    end

    fcorrectmat(i) = fcorrect/fctest;
    fPD(i) = fD/1520;
    fPFA(i) = fFA/960;
    
    i = i+1;
end
thresholdx = 0:1e-6:1;
figure
plot(thresholdx,fcorrectmat);
title('Feature correct rate VS threshold')

% figure 
% plot(PD,'b')
% hold on
% plot(PFA,'r')

figure
plot(fPFA,fPD,'b');
hold on
plot(0:0.01:1,0:0.01:1,'r')
title('Feature ROC')
xlabel('PFA')
ylabel('PD')
% 
% % %% 
fMax = max(fcorrectmat)          
[r,fp]=find(fcorrectmat==fMax);
fp = fp(1)

fT = thresholdx(fp)


%% % % % MFCC

XM = [T12RCH1MFCC T17RCH1MFCC T18RCH1MFCC T19RCH1MFCC];
[mrtraining,mctraining] = size(XM);

mmiu = zeros(mrtraining,1);



for i = 1:mrtraining
    mmiu(i) = (sum(XM(i,:)))/mctraining;
   
end

% % Find Sigma
msigma = (XM-mmiu)*(XM-mmiu)';
msigma = diag(diag(msigma));
%% Features  
% %% Load test data and calculate pdf
% %% Generate testing matrix and add label
mtest1 = [T12RCH2MFCC T17RCH2MFCC T18RCH1MFCC T19RCH1MFCC T18RCH2MFCC T19RCH2MFCC];
% [mrtest1,mctest1] = size(mtest1);
% mlabel1 = ones(1,mctest1);
% 
% mtest2 = [T11RCH2MFCC T11RCH2MFCC TLab2RCH1MFCC TLab2RCH2MFCC TOther3RCH1MFCC TOther3RCH2MFCC];
% [mrtest2,mctest2] = size(mtest2);
% mlabel2 = 2*ones(1,mctest2);
% mlabel = [mlabel1 mlabel2];

mtest = MFCCAllMatrix;

[mrtest,mctest] = size(mtest);
mlabel = Label;
mpdfmat = zeros(1,mctest);
% %%
for number = 1:mctest
    mpdf = exp((-1/2)*((mtest(:,number) - mmiu)'*(inv(msigma))*(mtest(:,number) - mmiu)));
    mpdfmat(number) = mpdf;
end
threshold = 0:1e-6:1;
mcorrectmat = zeros(1,length(threshold));    
mPD = zeros(1,length(threshold));
mPFA = zeros(1,length(threshold));
i = 1;
for threshold =0:1e-6:1
    
    mcorrect = 0;    
    mD = 0;
    mFA = 0;

    for number = 1:mctest
%         fpdf = exp((-1/2)*((ftest(:,number) - fmiu)'*(inv(fsigma))*(ftest(:,number) - fmiu)));
%       pdf = (1/((sqrt(2*pi))^13))*(1/sqrt(det(sigma)))*exp((-1/2)*((test(:,number) - miu)'*(inv(sigma))*(test(:,number) - miu)));
%       pdf = (1/((sqrt(2*pi))^13))*(1/sqrt(det(sigma)))*exp((-1/2)*((test1(:,number) - miu)'*(inv(sigma))*(test1(:,number) - miu)));
        mpdf = mpdfmat(number);
        if mpdf > threshold;
            mclass = 1;
        end
        if mpdf <= threshold;
            mclass = 0;
        end
        if mclass == mlabel(number)
            mcorrect = mcorrect + 1;
        end
        if mlabel(number) == 1 && mclass == 1;
            mD = mD+1;
        end
        if mlabel(number) == 0 && mclass == 1;
            mFA = mFA+1;
        end
    end

    mcorrectmat(i) = mcorrect/mctest;
    mPD(i) = mD/1520;
    mPFA(i) = mFA/960;
    
    i = i+1;
end
thresholdx = 0:1e-6:1;
mMax = max(mcorrectmat)          
[r,mp]=find(mcorrectmat==mMax);
mp = mp(1)

mT = thresholdx(mp)


thresholdx = 0:1e-6:1;
figure
plot(thresholdx,mcorrectmat);
title(' MFCC correct rate VS threshold')

% figure 
% plot(mPD,'b')
% hold on
% plot(mPFA,'r')
% title('MFCC ROC')
figure
plot(mPFA,mPD,'b');
% hold on
% plot(0:0.01:1,0:0.01:1,'r')
title('MFCC ROC')
xlabel('PFA')
ylabel('PD')
% 
% 

%% contract
figure
plot(mPFA,mPD,'b');
hold on
plot(fPFA,fPD,'r');
legend('MFCC','Proposed Feature')
% title('MFCC ROC VS Feature')
xlabel('Probability of False Alarm')
ylabel('Probability of Detection')

figure

plot(thresholdx,mcorrectmat,'b');
hold on
plot(thresholdx,fcorrectmat,'r');
legend('MFCC','Proposed Feature')
xlabel('Threshold')
ylabel('Correct Rate')
% title('Correct rate MFCC VS Feature')