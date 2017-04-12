function CSFI_AUC
% 
%
% Calcurating ROC AUC 
% 
% SO@ACH 2017.4 
%% load data 
% G = readtable('Glc.csv');

T = readtable('Latest20170208.xlsx');

% remove subjects HFA reliability is low
rows =  T.FP< .15 & T.FN<.33 & T.FixLoss_pcnt<.2;
% rows =  T.FP< .15 & T.FixLoss_pcnt<.2;

T2 = T(rows,:);

N = readtable('Normal.csv');

%%
tbl = [10,4; 372, 200];
[h,p,stats] = fishertest(tbl)

%% OAG vs Normal

label = [T2.Type;N.Type];

% ROC curve 
figure; hold on;
[X,Y,t,AUC] = perfcurve(label, [T2.CSFI;N.CSFI_rate/100],'Normal');
plot(Y, X)
1-AUC; % 0.95

[X,Y,t,AUC] = perfcurve(label, [T2.MD30_2;N.MD_30_2],'Normal');
plot(X,Y)
AUC; % 0.90

[X,Y,t,AUC] = perfcurve(label, [T2.cpRNFL;N.cpRNFL],'Normal');
plot(X,Y)
AUC; % 0.96

[X,Y,t,AUC] = perfcurve(label, [T2.wRGC;N.wRGC],'Normal');
plot(X,Y)
AUC; % 0.95

[X,Y,t,AUC] = perfcurve(label, [T2.RGC_HFA;N.RGC_HFA],'Normal');
plot(X,Y,'--')
AUC; % 0.90

[X,Y,t,AUC] = perfcurve(label, [T2.RGC_OCT;N.RGC_OCT],'Normal');
plot(X,Y,'--')
AUC; % 0.95

xlabel('FP rate')
ylabel('TP rate')
title('ROC for diagnosing OAG')
legend({'CSFI','MD','cpRNFL','wRGC','RGC HFA','RGC OCT'})

%% NTG vs POAG
Type = char(T2.Type);
rows = zeros(length(T2.Type),1);

% Pick up POAG
for n = 1: length(T2.Type);
 rows(n) = strcmp(T2.Type(n),'POAG');
 rows    = logical(rows);
end

POAG = T2(rows,:);

% pick up NTG
for n = 1: length(T2.Type);
 rows(n) = strcmp(T2.Type(n),'NTG');
 rows    = logical(rows);
end
NTG  =  T2(rows,:);
