function RGCfrom10_2
%% read data
cd '/Users/shumpei/Google Drive/CSFI/CSFI';
T = readtable('10-2forMatlab2.xlsx');

%% remove subjects HFA reliability is not enough
rows =  T.FP< .15 & T.FN<.33 & T.FixLoss_pcnt<.2;

T2 = T(rows,:);

%% plot 1
figure; hold on;
plot(T2.MD10_2,T2.RGC_OCT,'ob')
plot(T2.MD10_2,T2.RGC_HFA_Adjusted,'or')
plot(T2.MD10_2,T2.wRGC_adjusted,'og')
xlabel('MD value')
ylabel('estimated num of RGC')
legend({'OCT','HFA','weighted'})


%% Bland Altman From matlab exchange

addpath(genpath('BlandAltman'))
label ={'RGC HFA','RGC OCT'};
[rpc, fig, stats] = BlandAltman(T2.RGC_OCT,T2.RGC_HFA_Adjusted,label)

%% lowess Y= MD
figure; hold on;
span = 0.5;
[xx, inds] = sort(T2.RGC_HFA_Adjusted);
yy1 = smooth(T2.RGC_HFA_Adjusted,T2.MD10_2,span,'rloess');
plot(xx,yy1(inds),'r.');
clear xx inds

[xx, inds] = sort(T2.RGC_OCT);
yy2 = smooth(T2.RGC_OCT,T2.MD10_2,span,'rloess');
plot(xx,yy2(inds),'b.')
clear xx inds

[xx, inds] = sort(T2.wRGC);
yy3 = smooth(T2.wRGC_adjusted,T2.MD10_2,span,'rloess');
plot(xx,yy3(inds),'g.')
clear xx inds

xlabel 'RGC count'
ylabel MD
legend({'RGC HFA','RGC OCT','wRGC'})

%% lowess X= MD
figure; hold on;
span = 0.5;
[xx, inds] = sort(T2.MD10_2);
yy1 = smooth(T2.MD10_2, T2.RGC_HFA_Adjusted, span,'rloess');
plot(xx,yy1(inds),'r.');
% clear xx inds

% [xx, inds] = sort(T2.RGC_OCT);
yy2 = smooth(T2.MD10_2,T2.RGC_OCT,span,'rloess');
plot(xx,yy2(inds),'b.')
% clear xx inds

% [xx, inds] = sort(T2.wRGC_adjusted);
yy3 = smooth(T2.MD10_2,T2.wRGC_adjusted,span,'rloess');
plot(xx,yy3(inds),'g.')
clear xx inds

xlabel 'RGC count'
ylabel MD
legend({'RGC HFA','RGC OCT','wRGC'})

