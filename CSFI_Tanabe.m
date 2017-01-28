function CSFI_Tanabe

%% Load whole data

T = readtable('CSFI_data.csv');
C = readtable('Normal.csv');
N = readtable('NTG.csv');
P = readtable('POAG.csv');
G = readtable('Glc.csv');
%% Age
figure; hold on;
c = jet(20);

box = nan(size(P.age));

box1 = box;
box2 = box;

box1(1:length(C.age))  = C.age;
box2(1:length(N.age)) = N.age;

boxplot([box1,box2, P.age])
hold off

%% Age POAGvs NTG
figure; hold on;
boxplot(G.age,G.Type)
hold off

[p, h] = ranksum(N.age, P.age)

%% CSFI vs VFI
figure; hold on;
c =lines(4);

plot(G.CSFI_rate,G.VFI_rate,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))
xlabel 'CSFI %'
ylabel 'VFI %'

legend('Total n=661')

%% CSFI vs MD
figure; hold on;
c =lines(4);

plot(G.CSFI_rate,G.VFI_rate,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))

plot(G.CSFI_rate,G.MD_30_2,'o','color',[0 0 0],'MarkerFaceColor',c(2,:))

xlabel 'CSFI %'
ylabel 'VFI-MD'

legend({'VFI','MD'})

%% CSFI vs MD subplot
figure; subplot(1,2,1);hold on;
c =lines(4);

plot(G.CSFI_rate,G.VFI_rate,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))

xlabel 'CSFI %'
ylabel 'VFI %'
legend({'VFI'})


subplot(1,2,2)
plot(G.CSFI_rate,G.MD_30_2,'o','color',[0 0 0],'MarkerFaceColor',c(2,:))

xlabel 'CSFI %'
ylabel 'MD'

legend({'MD'})

%% CSFIvs MD
figure; hold on;
plot(G.CSFI_rate,G.MD_30_2,'o','color',[0 0 0],'MarkerFaceColor',c(2,:))

xlabel 'CSFI %'
ylabel 'MD'

legend({'MD'})

%%
[h, p] =corr(G.CSFI_rate,G.VFI_rate)

[h, p] =corr(G.CSFI_rate,G.MD_30_2)


[h, p] =corr(G.VFI_rate,G.MD_30_2)

%%
Y = fieldnames(G);
G.Axial_length =cell2mat(G.Axial_length);
% for kk= 1:length(Y)
for kk = [8,11,16,17,20,25]
    try
    figure;hold on;
    G.(Y{kk})(G.(Y{kk})==0)=nan;
    plot(G.CSFI_rate,G.(Y{kk}),'o','color',[0 0 0],'MarkerFaceColor',c(kk,:))
    
    xlabel 'CSFI %'
    ylabel(Y{kk})
    catch
    end
end

%%
% legend({'MD'})
figure;hold on;
for kk = [8,11,17]
    try
    G.(Y{kk})(G.(Y{kk})==0)=nan;
%     plot(G.CSFI_rate,G.(Y{kk}),'o','color',[0 0 0],'MarkerFaceColor',c(kk,:))
    
    plot(G.CSFI_rate,G.(Y{kk}),'o','color','none','MarkerFaceColor',c(kk,:))

    
    xlabel 'CSFI %'
    ylabel(Y{kk})
    legend({Y{8},Y{11},Y{17}})
    catch
    end
end
%% VFI_rate vs MD
figure; hold on;
c =lines(4);

plot(G.MD_30_2,G.VFI_rate,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))

% plot(G.CSFI_rate,G.MD_30_2,'o','color',[0 0 0],'MarkerFaceColor',c(2,:))

xlabel 'MD'
ylabel 'VFI'

%
lsline

mdl = fitlm(G.MD_30_2,G.VFI_rate)

% ANOVA ??????????????????????????????????????????
anova(mdl,'summary')

% ?????????????????????
coefCI(mdl)

% ??????????????????????????????

[p,F,d] = coefTest(mdl)

% legend({'VFI','MD'})

%% Stage
inds = G.MD_30_2>=-6;
x = G.CSFI_rate(inds);
y = G.VFI_rate(inds);

figure;hold on;
plot(x,y,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))
xlabel CSFI
ylabel VFI

title Early
lsline

mdl = fitlm(x,y)

[h,p] = corr(x, y)

% MD
x = G.CSFI_rate(inds);
y = G.MD_30_2(inds);

figure;hold on;
plot(x,y,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))
xlabel CSFI
ylabel MD

title Early
% lsline

mdl = fitlm(x,y)

[h,p] = corr(x, y)

%% fit 2d VFI
x = G.CSFI_rate(inds);
y = G.VFI_rate(inds);

[p,ErrorEst] = polyfit(x, y,2);
y_fit = polyval(p,x,ErrorEst);

figure; hold on;
% each subject
plot(x,y,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))

% fit
plot(x,y_fit,'.');


% figure;hold on;
% plot(G.CSFI_rate(inds), G.VFI_rate(inds),'o',t2,y2)
xlabel 'CSFI %'
ylabel VFI

% % figure; hold on;
% res = y - y_fit;
% plot(x,res,'o')

%
[y_fit,delta] = polyval(p,x,ErrorEst);

% figure; hold on;
%
% plot(x,y,'o')
% plot(x,y_fit,'g.')

% 95% tile
plot(x,y_fit+2*delta,'g.',...
    x,y_fit-2*delta,'g.');

corrcoef(x,y)

[h,p] = corr(G.CSFI_rate(inds), G.VFI_rate(inds))


%% fit 2d MD
x = G.CSFI_rate(inds);
y = G.MD_30_2(inds);

[p,ErrorEst] = polyfit(x, y,2);
y_fit = polyval(p,x,ErrorEst);

figure; hold on;
% each subject
plot(x,y,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))

% fit
plot(x,y_fit,'.');


% figure;hold on;
% plot(G.CSFI_rate(inds), G.VFI_rate(inds),'o',t2,y2)
xlabel 'CSFI %'
ylabel MD

% % figure; hold on;
% res = y - y_fit;
% plot(x,res,'o')

%
[y_fit,delta] = polyval(p,x,ErrorEst);

% figure; hold on;
%
% plot(x,y,'o')
% plot(x,y_fit,'g.')

% 95% tile
plot(x,y_fit+2*delta,'g.',...
    x,y_fit-2*delta,'g.');

corrcoef(x,y)

[h,p] = corr(G.CSFI_rate(inds), G.VFI_rate(inds))

return
% SO 2017/1/25

%% CSFI vs VFI
% plot(C.CSFI_rate,C.VFI_rate,'o')
% plot(N.CSFI_rate,N.VFI_rate,'o')
% plot(P.CSFI_rate,P.VFI_rate,'o')
figure;
hold on;

plot(C.CSFI_rate,C.VFI_rate,'o','color',c(1,:),'MarkerFaceColor',c(1,:))

plot(G.CSFI_rate,G.VFI_rate,'o','color',[0 0 0],'MarkerFaceColor',c(2,:))


plot(N.CSFI_rate,N.VFI_rate,'o','color',[0 0 0],'MarkerFaceColor',c(2,:))
plot(P.CSFI_rate,P.VFI_rate,'o','color',[0 0 0],'MarkerFaceColor',c(3,:))

xlabel 'CSFI %'
ylabel 'VFI %'

legend({'Healthy','NTG','POAG'})
legend({'NTG','POAG'})


%% NTG vs POAG
figure; hold on;
c =lines(4);

plot(N.CSFI_rate,N.VFI_rate,'o','color',[0 0 0],'MarkerFaceColor',c(2,:))
plot(P.CSFI_rate,P.VFI_rate,'o','color',[0 0 0],'MarkerFaceColor',c(3,:))

xlabel 'CSFI %'
ylabel 'VFI %'

legend({'NTG','POAG'})


%% CSFI ???3????????????????????????
[h,p] = ttest2(C.CSFI_rate,N.CSFI_rate) % S
[h,p] = ttest2(C.CSFI_rate,P.CSFI_rate) % S
[h,p] = ttest2(N.CSFI_rate,P.CSFI_rate) % N.S.

% NTG???POAG?????????????????????

%%
[h,p] = corr(C.age,C.CSFI_rate) % n.s
[h,p] = corr(N.age,N.CSFI_rate) % n.s
[h,p] = corr(P.age,P.CSFI_rate) % S

%% CSFI vs MD
figure; hold on;
c = jet(4);

plot(C.CSFI_rate,C.MD_24_2,'o','color',c(1,:),'MarkerFaceColor',c(1,:))
plot(N.CSFI_rate,N.MD_24_2,'o','color',c(2,:),'MarkerFaceColor',c(2,:))
plot(P.CSFI_rate,P.MD_24_2,'o','color',c(3,:),'MarkerFaceColor',c(3,:))

xlabel 'CSFI %'
ylabel 'MD'

legend({'Healthy','NTG','POAG'})

%% MD_24_2 ???3????????????????????????
[h,p] = ttest2(C.MD_24_2,N.MD_24_2) % S
[h,p] = ttest2(C.MD_24_2,P.MD_24_2) % S
[h,p] = ttest2(N.MD_24_2,P.MD_24_2) % S P= 0.008

% 3????????????????????????

%%


% load data
G = readtable('Glc.csv');


%% MD vs VFI

figure; hold on;
c = jet(4);

% plot(G.CSFI_rate,G.VFI_rate,'o')

xlabel 'CSFI %'
ylabel 'VFI %'

% fit 2d
p = polyfit(G.CSFI_rate,G.VFI_rate,2);

t2 = -0:1:100;
y2 = polyval(p,t2);

% figure;hold on;
plot(G.CSFI_rate,G.VFI_rate,'o',t2,y2)
xlabel 'CSFI %'
ylabel VFI

[h,p] = corr(G.CSFI_rate,G.VFI_rate)

%%

figure; hold on;
c = jet(4);

% fit 2d
p = polyfit(G.CSFI_rate,G.MD_24_2,2);

t2 = 0:1:100;
y2 = polyval(p,t2);

% figure;hold on;
plot(G.CSFI_rate,G.MD_24_2,'o',t2,y2)
xlabel 'CSFI %'
ylabel 'MD'


[h,p] = corr(G.CSFI_rate,G.MD_24_2)
