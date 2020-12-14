clc
clear all

load('Rdata.mat')
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
%%K-NN
%% USe loop to calculate
%% Use first all as training data
%% Find miu

% % Load training data
%% Select T12  T17 T18 T19 as Native 
%%
% X = [T12RCH1Features T12RCH2Features T17RCH1Features T17RCH2Features];
ft1 = [FeatureAllMatrix(:,1:80)];%, FeatureAllMatrix(:,81:160) ];
% ft1 = [T12RCH1Features T17RCH1Features T18RCH1Features T19RCH1Features];
[frt1,fct1] = size(ft1);
fmiu1 = zeros(frt1,1);
for fi = 1:frt1
    fmiu1(fi) = (sum(ft1(fi,:)))/fct1; 
end
%% Select T11  TLab(2) TLab(3) TOther(3)as Nonative
% ft2 = [FeatureAllMatrix(:,1521:1600), FeatureAllMatrix(:,1681:1760) ];
% ft2 = [T11RCH1Features, TLab2RCH1Features, TLab3RCH1Features, TOther3RCH1Features];
ft2 = [TLab2RCH1Features TLab2RCH2Features];
[frt2,fct2] = size(ft2);
fmiu2 = zeros(frt2,1);
for fi = 1:frt2
    fmiu2(fi) = (sum(ft2(fi,:)))/fct2; 
end
% 
fctest = length(FeatureAllMatrix);
fdis1 = zeros(1,fctest);
fdis2 = zeros(1,fctest);
for number = 1:fctest
    
%     pdf = (1/((sqrt(2*pi))^13))*(1/sqrt(det(sigma)))*exp((-1/2)*((test(:,number) - miu)'*(inv(sigma))*(test(:,number) - miu)));
% pdf = (1/((sqrt(2*pi))^13))*(1/sqrt(det(sigma)))*exp((-1/2)*((test1(:,number) - miu)'*(inv(sigma))*(test1(:,number) - miu)));
    fdis1(number) = sqrt(sum((FeatureAllMatrix(:,number)-fmiu1).^2));
    fdis2(number) = sqrt(sum((FeatureAllMatrix(:,number)-fmiu2).^2));

end
fthreshold = 0:1e-5:4 ;
fN = length(fthreshold);
fPdmat = zeros(1,fN);
fPfamat = zeros(1,fN);
fCorrectratemat = zeros(1,fN);
fi=1;
% 
for fthreshold = 0:1e-5:4
    fcorrect = 0;
    ffa = 0;
    fd = 0;

    
    for number = 1:fctest

        if fdis1(number) <= fdis2(number)*fthreshold
            fclass = 1;
        else 
            fclass = 0;
        end
        if fclass == Label(number);
            fcorrect = fcorrect+1;
        end
        if Label(number) == 0 && fclass == 1
            ffa = ffa+1;
        end
        if Label(number) == 1 && fclass == 1
            fd = fd+1;
        end
    end
    fCorrectratemat(fi) = fcorrect/fctest;
    fPfamat(fi) = ffa/960;
    fPdmat(fi) = fd/1520;
    fi=fi+1;
end
fthreshold = 0:1e-5:4 ;
figure
plot(fthreshold,fCorrectratemat)
xlabel('Threshold')
ylabel('Correct Rate')

figure
plot(fPfamat,fPdmat)

fMax = max(fCorrectratemat)          
[fr,fc]=find(fCorrectratemat==fMax);
fc = fc(1)


fthreshold=fthreshold(fc)

%% Threshold = 1
%  fcorrect = 0;
%     ffa = 0;
%     fd = 0;
% 
%     
%     for number = 1:fctest
% 
%         if fdis1(number) <= fdis2(number)
%             fclass = 1;
%         else 
%             fclass = 0;
%         end
%         if fclass == Label(number);
%             fcorrect = fcorrect+1;
%         end
%         if Label(number) == 0 && fclass == 1
%             ffa = ffa+1;
%         end
%         if Label(number) == 1 && fclass == 1
%             fd = fd+1;
%         end
%     end
%     fCorrectratemat = fcorrect/fctest
%     fPfamat = ffa/960
%     fPdmat = fd/1520

%% MFCC
%% Find miu
%% Select T12 & T17 as Native 
%% Select T11 & TLab(2) as Nonative


mt1 = [MFCCAllMatrix(:,1:80)];% MFCCAllMatrix(:,81:160) ]; 
% mt1 = [ T12RCH1MFCC, T17RCH1MFCC, T18RCH1MFCC, T19RCH1MFCC];

[mrt1,mct1] = size(mt1);
mmiu1 = zeros(mrt1,1);
for mi = 1:mrt1
    mmiu1(mi) = (sum(mt1(mi,:)))/mct1; 
end
% 
% mt2 = [MFCCAllMatrix(:,1521:1600), MFCCAllMatrix(:,1681:1760) ];
% mt2 = [T11RCH1MFCC, TLab2RCH1MFCC, TLab3RCH1MFCC, TOther3RCH1MFCC];
mt2 = [TLab2RCH1MFCC TLab2RCH2MFCC];
[mrt2,mct2] = size(mt2);
mmiu2 = zeros(mrt2,1);
for mi = 1:mrt2
    mmiu2(mi) = (sum(mt2(mi,:)))/mct2; 
end
% 
mctest = length(MFCCAllMatrix);
mdis1 = zeros(1,mctest);
mdis2 = zeros(1,mctest);
for number = 1:mctest
    
%     pdf = (1/((sqrt(2*pi))^13))*(1/sqrt(det(sigma)))*exp((-1/2)*((test(:,number) - miu)'*(inv(sigma))*(test(:,number) - miu)));
% pdf = (1/((sqrt(2*pi))^13))*(1/sqrt(det(sigma)))*exp((-1/2)*((test1(:,number) - miu)'*(inv(sigma))*(test1(:,number) - miu)));
    mdis1(number) = sqrt(sum((MFCCAllMatrix(:,number)-mmiu1).^2));
    mdis2(number) = sqrt(sum((MFCCAllMatrix(:,number)-mmiu2).^2));

end


%
mthreshold = 0:1e-5:4 ;
mN = length(mthreshold);
mPdmat = zeros(1,mN);
mPfamat = zeros(1,mN);
mCorrectratemat = zeros(1,mN);
mi=1;
% 
for mthreshold = 0:1e-5:4
    mcorrect = 0;
    mfa = 0;
    md = 0;

    
    for number = 1:mctest

        if mdis1(number) <= mdis2(number)*mthreshold
            mclass = 1;
        else 
            mclass = 0;
        end
        if mclass == Label(number);
            mcorrect = mcorrect+1;
        end
        if Label(number) == 0 && mclass == 1
            mfa = mfa+1;
        end
        if Label(number) == 1 && mclass == 1
            md = md+1;
        end
    end
    mCorrectratemat(mi) = mcorrect/mctest;
    mPfamat(mi) = mfa/960;
    mPdmat(mi) = md/1520;
    mi=mi+1;
end
mthreshold = 0:1e-5:4 ;
figure
plot(mthreshold,mCorrectratemat,'b')
hold on
plot(mthreshold,fCorrectratemat,'r')
legend('MFCC','Proposed Feature')

xlabel('Threshold')
ylabel('Correct Rate')
figure
plot(mPfamat,mPdmat)
mMax = max(mCorrectratemat)          
[r,mc]=find(mCorrectratemat==mMax);
mc = mc(1)

mthreshold = mthreshold(mc)

figure
plot(mPfamat,mPdmat,'b')
hold on
plot(fPfamat,fPdmat,'r')
legend('MFCC','Proposed Feature')
xlabel('Probability of False Alarm')
ylabel('Probability of Detection')
%% Threshold = 1
% mcorrect = 0;
%     mfa = 0;
%     md = 0;
% 
%     
%     for number = 1:mctest
% 
%         if mdis1(number) <= mdis2(number)
%             mclass = 1;
%         else 
%             mclass = 0;
%         end
%         if mclass == Label(number)
%             mcorrect = mcorrect+1;
%         end
%         if Label(number) == 0 && mclass == 1
%             mfa = mfa+1;
%         end
%         if Label(number) == 1 && mclass == 1
%             md = md+1;
%         end
%     end
%     mCorrectratemat = mcorrect/mctest
%     mPfamat = mfa/960
%     mPdmat = md/1520
% 

