function RGCfrom10_2
%% read data

T = readtable('10-2forMatlab.xlsx');

% %%
% figure; hold on;
% plot(T.RGC_OCT,T.MD10_2,'ob')
% plot(T.RGC_HFA_Adjusted,T.MD10_2,'or')
% plot(T.wRGC_adjusted, T.MD10_2,'og')
% plot(T.wRGC, T.MD10_2,'ob')

%% plot 1
figure; hold on;
plot(T.MD10_2,T.RGC_OCT,'ob')
plot(T.MD10_2,T.RGC_HFA_Adjusted,'or')
plot( T.MD10_2,T.wRGC_adjusted,'og')
xlabel('MD value')
ylabel('estimated num of RGC')
legend({'OCT','HFA','weighted'})

%%
figure
plot(T.RGC_OCT,T.RGC_HFA_Adjusted,'o')
