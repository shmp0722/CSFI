% load data 
% G = readtable('Glc.csv');

T = readtable('Latest20170208.xlsx');


% TF = ismissing(T, {'','.','NA' NaN -99});
% T(any(TF,2),:)
% 
% T = standardizeMissing(T,-99);


% %% take complete data
% TF = ismissing(T);
% T2 = T(~any(TF,2),:);

%%
Var = fieldnames(T);

for ii= 8:length(Var)
    [R(ii),p(ii) ] = corr(T.FixLoss_pcnt,T.(Var{ii}))
    figure;
    plot(T.FixLoss_pcnt,T.(Var{ii}),'o')
    lsline
    title(Var{ii})
end



plot(T.FixLoss_pcnt*100,T.CSFI*100,'o')


%% MD vs VFI

% remove subjects reliability low
rows =  T.FP< .2 & T.FN<.33 & T.FixLoss_pcnt<.2;

T2 = T(rows,:);
summary(T2)

% 
figure; hold on;
c = jet(4);

% plot(T2.CSFI,T2.VFI,'o')

xlabel 'CSFI %'
ylabel 'VFI %'

% fit model
for jj = 1:4 
p = polyfit(T2.CSFI,T2.VFI,jj);
% p = polyfit(T2.CSFI,T2.VFI,2);
% p = polyfit(T2.CSFI,T2.VFI,3);
% p = polyfit(T2.CSFI,T2.VFI,4);

t2 = -0.4:0.01:1;
y2 = polyval(p,t2);

% figure;hold on;
plot(T2.CSFI,T2.VFI,'o',t2,y2)%,'color',c(jj,:))
xlabel 'CSFI %'
ylabel VFI
end
% legend({'1','2','3','4'})

%%
[h,p] = corr(T2.CSFI,T2.VFI)
[h,p] = corr(T2.CSFI,T2.MD30_2)

%%

figure; hold on;
c = jet(4);

% fit 2d 
% p = polyfit(T2.CSFI,T2.MD30_2,2);


% fit model
for jj = 1:4 
p = polyfit(T2.CSFI,T2.MD30_2,jj);
% p = polyfit(T2.CSFI,T2.VFI,2);
% p = polyfit(T2.CSFI,T2.VFI,3);
% p = polyfit(T2.CSFI,T2.VFI,4);

t2 = -0.4:0.01:1;
y2 = polyval(p,t2);

% figure;hold on;
plot(T2.CSFI,T2.MD30_2,'o',t2,y2)%,'color',c(jj,:))
xlabel 'CSFI %'
ylabel MD
end







t2 = 0:1:100;
y2 = polyval(p,t2);

% figure;hold on;
plot(G.CSFI_rate,G.MD_24_2,'o',t2,y2)
xlabel 'CSFI %'
ylabel 'MD'


[h,p] = corr(G.CSFI_rate,G.MD_24_2)
