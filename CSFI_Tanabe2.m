%% load data 
% G = readtable('Glc.csv');

T = readtable('Latest20170208.xlsx');

%% remove subjects HFA reliability is low
rows =  T.FP< .2 & T.FN<.33 & T.FixLoss_pcnt<.2;

T2 = T(rows,:);

figure; hold on;
subplot(1,2,1)
plot(T2.CSFI,T2.MD30_2,'*r') % n = 575
title 'Good HFA n =575 '  
xlabel CSFI
ylabel MD30-2

subplot(1,2,2)
plot(T.CSFI,T.MD30_2,'*') % n = 661
title 'all subject n=661' 
xlabel CSFI
ylabel MD30-2
% NTG = T{rows,'Type'}
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

%%

mean(NTG.age)
nanstd(NTG.age)

mean(POAG.age)
nanstd(POAG.age)

mean(NTG.MD30_2)
mean(POAG.MD30_2)
p = ranksum(POAG.MD30_2,NTG.MD30_2)


mean(NTG.VFI)
mean(POAG.VFI)
std(NTG.VFI)

p = ranksum(POAG.VFI,NTG.VFI)



% same as Man-U test
p = ranksum(POAG.age,NTG.age)
clear p

%% NTG Gender

[tbl,chi2,p,labels] = crosstab(T2.Type,T2.Gender)
[h,p,stats] = fishertest(tbl)


%% NTG vs POAG 
figure; hold on;
plot(NTG.CSFI, NTG.MD30_2,'*r')
plot(POAG.CSFI,POAG.MD30_2,'*b')


title 'NTG vs POAG' 
xlabel CSFI
ylabel MD30-2

legend({'NTG','POAG'})

%% ANCOV
% there is no difference due to glc type
[h,atab,ctab,stats] = aoctool(T2.CSFI, T2.MD30_2,T2.Type);

[h,atab,ctab,stats] = aoctool(T2.CSFI, T2.VFI,T2.Type);

%%
% AIC
figure; hold on;
% scatter(T2.CSFI,T2.MD30_2,10,'*')

plot(NTG.CSFI, NTG.MD30_2,'*r')
plot(POAG.CSFI,POAG.MD30_2,'*b')


AIC = zeros(1,4);
obj = cell(1,4);
for k = 1:4
    figure; hold on;
    plot(NTG.CSFI, NTG.MD30_2,'*r')
    plot(POAG.CSFI,POAG.MD30_2,'*b')
    obj{k} = gmdistribution.fit([T2.CSFI,T2.MD30_2],k);
    AIC(k)= obj{k}.AIC;

h = ezcontour(@(x,y)pdf(obj{k},[x y]));
axis auto
end


% %%
% varfun(@mean,T,'InputVariables','CSFI',...
%     'GroupingVariables','Type')


%% Fix_loss effect something on CSFI
Var = fieldnames(T2);

for ii= 8:length(Var)-1
    [R(ii),p(ii) ] = corr(T2.FixLoss_pcnt,T2.(Var{ii}));

    % showing only significantly correlatied combinations
    if p(ii)<0.05
        figure;
        plot(T2.FixLoss_pcnt,T2.(Var{ii}),'o')
        lsline
        title(Var{ii})
        xlabel 'Fix loss %'
        ylabel(Var{ii})
    end
end

% plot(T.FixLoss_pcnt*100,T.(12),'o')


%% CSFI vs MD ?? Why??
figure; hold on;
c = jet(4);

plot(T2.CSFI,T2.MD30_2,'o')
xlabel 'CSFI %'
ylabel 'MD 30-2'

% fit model
for jj = 1:4 
[p,s, mu] = polyfit(T2.CSFI,T2.MD30_2,jj);
t2 = linspace(-0.4,1);
y2 = polyval(p,t2);

% figure;hold on;
% plot(T2.CSFI,T2.MD30_2,'o',t2,y2)%,'color',c(jj,:))
plot(t2,y2)%,'color',c(jj,:))

% lsline
xlabel 'CSFI %'
ylabel VFI
end
% legend({'1','2','3','4'})
clear p s
%% AIC MD!?
figure; hold on;
scatter(T2.CSFI,T2.MD30_2,10,'*')

AIC = zeros(1,4);
obj = cell(1,4);
for k = 1:4
    obj{k} = gmdistribution.fit([T2.CSFI,T2.MD30_2],k);
    AIC(k)= obj{k}.AIC;
end

% h = ezcontour(@(x,y)pdf(obj{1},[x y]));



%% Fitting exponential did not fit 

[estimates, model] = fitcurvedemo(T2.CSFI,T2.VFI);

figure; hold on;
plot(T2.CSFI,T2.VFI, '*')
% hold on
[sse, FittedCurve] = model(estimates);
plot(T2.CSFI, FittedCurve, 'r')
 
xlabel('T2.CSFI')
ylabel('f(estimates,T2.CSFI)')
title(['Fitting to function ', func2str(model)]);
legend('data', ['fit using ', func2str(model)])
hold off


%%
[h,p] = corr(T2.CSFI,T2.VFI)
[h,p] = corr(T2.CSFI,T2.MD30_2)

%% CSFI vs MD30-2

figure; hold on;
c = lines(4);
plot(T2.CSFI,T2.MD30_2,'*b')


% fit model
for jj = 1:4
p = polyfit(T2.CSFI,T2.MD30_2,jj);

t2 = linspace(-0.4,1);
y2 = polyval(p,t2);

plot(t2,y2,'color',c(jj,:))

xlabel 'CSFI %'
ylabel MD
end
legend({'','2d','3d','4d'})

%% AIC VFI
figure; hold on;
scatter(T2.CSFI,T2.MD30_2,10,'*')

AIC = zeros(1,4);
obj = cell(1,4);
for k = 1:4
    obj{k} = gmdistribution.fit([T2.CSFI,T2.VFI],k);
    AIC(k)= obj{k}.AIC;
end

min(AIC)

h = ezcontour(@(x,y)pdf(obj{4},[x y]));
axis auto

%% staging
inds = T2.MD30_2>-6; 

x = T2.CSFI(inds);
y = T2.MD30_2(inds);

figure;hold on;
plot(x,y,'o','color',[0 0 0])%,'MarkerFaceColor',c(1,:))
xlabel CSFI
ylabel MD

lsline

mdl = fitlm(x,y)

[h,p] = corr(x, y)
title Early


% VFI
y = T2.VFI(inds);

figure;hold on;
plot(x,y,'o','color',[0 0 0])%,'MarkerFaceColor',c(1,:))
xlabel CSFI
ylabel VFI

title Early
lsline

mdl = fitlm(x,y)

[h,p] = corr(x, y)


%% stage 2
inds = T2.MD30_2<-6 & T2.MD30_2>=-12 ; 

x = T2.CSFI(inds);
y = T2.MD30_2(inds);

figure;hold on;
plot(x,y,'o','color',[0 0 0])%,'MarkerFaceColor',c(1,:))
xlabel CSFI
ylabel MD

lsline

mdl = fitlm(x,y)

[h,p] = corr(x, y)
title Middle


% VFI
y = T2.VFI(inds);

figure;hold on;
plot(x,y,'o','color',[0 0 0])%,'MarkerFaceColor',c(1,:))
xlabel CSFI
ylabel VFI

title Middle
lsline

mdl = fitlm(x,y)

[h,p] = corr(x, y)

%% stage 3
inds = T2.MD30_2<-12 ; 

x = T2.CSFI(inds);
y = T2.MD30_2(inds);

figure;hold on;
plot(x,y,'o','color',[0 0 0])%,'MarkerFaceColor',c(1,:))
xlabel CSFI
ylabel MD

lsline

mdl = fitlm(x,y)

[h,p] = corr(x, y)
title Advance


% VFI
y = T2.VFI(inds);

figure;hold on;
plot(x,y,'o','color',[0 0 0])%,'MarkerFaceColor',c(1,:))
xlabel CSFI
ylabel VFI

title Advance
lsline

mdl = fitlm(x,y)

[h,p] = corr(x, y)




%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load whole data

% T = readtable('CSFI_data.csv');
C = readtable('Normal.csv');
N = readtable('NTG.csv');
P = readtable('POAG.csv');
G = readtable('Glc.csv');
%% Age
figure; hold on;
c = jet(100);

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

[p, h] = ranksum(N.MD_30_2, P.MD_30_2)

[p, h] = ranksum(N.VFI_rate, P.VFI_rate)


nanmean(N.MD_30_2)
nanmean(P.MD_30_2)

nanmean(N.VFI_rate)
nanmean(P.VFI_rate)

%% Gender
% NTG f 117: m 209
% POAG 109:226


%% CSFI vs VFI
figure; hold on;
c =lines(100);

plot(G.CSFI_rate,G.VFI_rate,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))
xlabel 'CSFI %'
ylabel 'VFI %'

legend('Total n=661')

%% CSFI vs MD
figure; hold on;
c =lines(100);

plot(G.CSFI_rate,G.VFI_rate,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))
plot(G.CSFI_rate,G.MD_30_2,'o','color',[0 0 0],'MarkerFaceColor',c(2,:))

xlabel 'CSFI %'
ylabel 'VFI-MD'

legend({'VFI','MD'})

%% CSFI vs MD subplot
figure; subplot(1,2,1);hold on;
c =lines(100);

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

%% Axial_length
for jj = 1:length(G.Axial_length)
    G.AL(jj) = str2double(G.Axial_length(jj));
end

Y = fieldnames(G);

% for kk= 1:length(Y)
for kk = [8,11,16,17,20,30]
    try
    figure;hold on;
    plot(G.CSFI_rate,G.(Y{kk}),'o','color',[0 0 0],'MarkerFaceColor',c(kk,:))
    
    xlabel 'CSFI %'
    ylabel(Y{kk})
    catch
    end
    % save the fig
    saveas(gca,Y{kk},'png')
    saveas(gca,Y{kk},'fig')
    saveas(gca,Y{kk},'eps2')
end

!mv *fig Figure/
!mv *png Figure/
!mv *eps Figure/

%% AL
inds = G.AL==0;

figure; hold on;
plot(G.CSFI_rate(~inds),G.AL(~inds),'o','color',[0 0 0],'MarkerFaceColor',c(kk,:));
    lsline
    xlabel 'CSFI %'
    ylabel 'Axial length'

set(gca,'YLim',[21 30])


[h, p] =corr(G.CSFI_rate(~inds),G.AL(~inds))

mdl = fitlm(G.CSFI_rate(~inds),G.AL(~inds))

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

    saveas(gca,'3elements','png')
    saveas(gca,'3elements','fig')
    saveas(gca,'3elements','eps2')
%% cpRNFL vs CSFI
figure; hold on;
x = G.CSFI_rate; y = G.cpRNFL;
plot(x, y,'o')
xlabel 'CSFI %'
ylabel 'cpRNFL'

%% cpRNFL vs MD
figure; hold on;
x = G.MD_30_2; y = G.cpRNFL;
plot(x, y,'o')
xlabel 'MD'
ylabel 'cpRNFL'

%% cpRNFL vs MD
figure; hold on;
x = G.VFI_rate; y = G.cpRNFL;
plot(x, y,'o')
xlabel 'VFI'
ylabel 'cpRNFL'
    
    
%% VFI_rate vs MD
figure; hold on;
c =lines(4);

plot(G.MD_30_2,G.VFI_rate,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))

% plot(G.CSFI_rate,G.MD_30_2,'o','color',[0 0 0],'MarkerFaceColor',c(2,:))

xlabel 'MD'
ylabel 'VFI'

%
lsline

[h,p] = corr(G.MD_30_2,G.VFI_rate)

mdl = fitlm(G.MD_30_2,G.VFI_rate)

% ANOVA ??????????????????????????????????????????
anova(mdl,'summary')

% ?????????????????????
coefCI(mdl)

% ??????????????????????????????

[p,F,d] = coefTest(mdl)

% legend({'VFI','MD'})


%% Stage 1
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


%% stage 2
inds = G.MD_30_2<-6 & G.MD_30_2>=-12 ; 

x = G.CSFI_rate(inds);
y = G.VFI_rate(inds);

figure;hold on;
plot(x,y,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))
xlabel CSFI
ylabel VFI

title Middle
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

title Middle
lsline

mdl = fitlm(x,y)

[h,p] = corr(x, y)

%% stage 3
inds = G.MD_30_2<-12 ; 

x = G.CSFI_rate(inds);
y = G.VFI_rate(inds);

figure;hold on;
plot(x,y,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))
xlabel CSFI
ylabel VFI

title Advance
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

title Advance
lsline

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


%% additional 3
c = jet(15);

for Min = -12:3:30;
    figure;
    inds =G.MD_30_2<=Min+3  & G.MD_30_2>Min ;
    
    x = G.CSFI_rate(inds);
    y = G.VFI_rate(inds);
    
    figure;hold on;
    plot(x,y,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))
    xlabel CSFI
    ylabel VFI
    
    title(sprintf('%d-%d',Min,Min+3))
    lsline
    
%     mdl = fitlm(x,y)
    
    [h,p] = corr(x, y)
    
    legend(sprintf('r = %d',h))
end

%% additional 3
c = jet(15);

for Min = -30:3:3;
    figure;
    inds =G.MD_30_2<=Min+3  & G.MD_30_2>Min ;
    
    x = G.CSFI_rate(inds);
    y = G.MD_30_2(inds);
    
    figure;hold on;
    plot(x,y,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))
    xlabel CSFI
    ylabel MD
    
    title(sprintf('%d-%d',Min,Min+3))
    lsline
    
%     mdl = fitlm(x,y)
    
    [h,p] = corr(x, y)
    
    legend(sprintf('r = %d',h))
end

%%
% MD
x = G.CSFI_rate(inds);
y = G.MD_30_2(inds);

figure;hold on;
plot(x,y,'o','color',[0 0 0],'MarkerFaceColor',c(1,:))
xlabel CSFI
ylabel MD

title Advance
lsline

mdl = fitlm(x,y)

[h,p] = corr(x, y)


