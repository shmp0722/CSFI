%% load data 
% G = readtable('Glc.csv');

T = readtable('Latest20170208.xlsx');

%% remove subjects HFA reliability is low
rows =  T.FP< .15 & T.FN<.33 & T.FixLoss_pcnt<.2;
% rows =  T.FP< .15 & T.FixLoss_pcnt<.2;


T2 = T(rows,:);

figure; hold on;
subplot(1,2,1)
plot(T2.CSFI,T2.MD30_2,'*r') % n = 575
title 'Good HFA n =572 '  
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

%% summary and difference

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

[tbl,chi2,p,labels] = crosstab(T2.Type,T2.Gender);
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

[h,atab,ctab,stats] = aoctool(T2.cpRNFL, T2.cpRNFL,T2.Type);
% [h,atab,ctab,stats] = aoctool(T2.CSFI, T2.VFI,T2.Type);

%% CSFI vs MD
figure; hold on;
% scatter(T2.CSFI,T2.MD30_2,10,'*')

plot(NTG.CSFI, NTG.MD30_2,'*r')
plot(POAG.CSFI,POAG.MD30_2,'*b')
xlabel CSFI
ylabel MD
legend({'NTG','POAG'})
hold off

n = 10;
AIC = zeros(1,n);
obj = cell(1,n);
for k = 1:n
%     figure; hold on;
%     plot(NTG.CSFI, NTG.MD30_2,'*r')
%     plot(POAG.CSFI,POAG.MD30_2,'*b')
    obj{k} = gmdistribution.fit([T2.CSFI,T2.MD30_2],k);
    AIC(k)= obj{k}.AIC;

% h = ezcontour(@(x,y)pdf(obj{k},[x y]));
% axis auto
end
figure; hold on;
plot(AIC,'r.')
line([0 10],[min(AIC),min(AIC)])
title AIC
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


%% AIC MD!?
figure; hold on;
scatter(T2.CSFI,T2.MD30_2,10,'*')
xlabel CSFI
ylabel MD
legend({'OAG'})
hold off;

n = 10;
AIC = zeros(1,n);
obj = cell(1,n);
for k = 1:n
    obj{k} = gmdistribution.fit([T2.CSFI,T2.MD30_2],k);
    AIC(k)= obj{k}.AIC;
end

% AIC
figure; hold on;
plot(AIC,'r.')
line([0 10],[min(AIC),min(AIC)])
title CSFI vs MD
ylabel AIC
xlabel n o function
hold off
% h = ezcontour(@(x,y)pdf(obj{1},[x y]));



%% Fitting exponential did not fit 

[estimates, model] = fitcurvedemo(T2.CSFI,T2.VFI);

figure; hold on;
plot(T2.CSFI,T2.VFI, '*')
% hold on
[sse, FittedCurve] = model(estimates);
plot(T2.CSFI, FittedCurve, 'r')
 
xlabel('CSFI')
ylabel('f(estimates,T2.CSFI)')
title(['Fitting to function ', func2str(model)]);
legend('data', ['fit using ', func2str(model)])
hold off


%% Just correlation 
[h,p] = corrcoef(T2.CSFI,T2.VFI)
[h,p] = corrcoef(T2.CSFI,T2.MD30_2)


%% CSFI vs MD30-2
% fit model
n =5;
for jj = 1:n
figure; hold on;
c = lines(n);
plot(T2.CSFI,T2.MD30_2,'*b')
[p,ErrorEst] = polyfit(T2.CSFI,T2.MD30_2,jj);

t2 = linspace(-0.4,1);
[y2,delta] = polyval(p,t2,ErrorEst);

% 95% c.i.
plot(t2,y2,'color',c(jj,:))
plot(t2,y2+2*delta,'r:')
plot(t2,y2-2*delta,'r:')

xlabel 'CSFI %'
ylabel MD
title(sprintf('FIT %d function with CI',jj))
hold off;
end

%% lowess
figure; hold on;
yy2 = smooth(T2.CSFI,T2.MD30_2,0.3,'rloess');
% plot(yy2)
[xx,ind] = sort(T2.CSFI);
plot(xx,T2.MD30_2(ind),'b.',xx,yy2(ind),'r-')
xlabel 'CSFI %'
ylabel MD
% legend({'POAG','2d','3d','4d'}

figure;
yy3 = smooth(T2.cpRNFL,T2.MD30_2,0.4,'rloess');
[xx,ind] = sort(T2.cpRNFL);
plot(T2.cpRNFL,T2.MD30_2,'b.',xx,yy3(ind),'r-')
xlabel 'cpRNFL'
ylabel MD


%% CSFI vs VFI
% AIC MD
figure; hold on;
scatter(T2.CSFI,T2.VFI,10,'*')
xlabel CSFI
ylabel VFI
legend({'POAG'})
hold off;

n = 10;
AIC = zeros(1,n);
obj = cell(1,n);
for k = 1:n
    obj{k} = gmdistribution.fit([T2.CSFI,T2.VFI],k);
    AIC(k)= obj{k}.AIC;
end

figure; hold on;
plot(AIC,'r.')
line([0 10],[min(AIC),min(AIC)])
% title AIC
xlabel 'N of function'
ylabel AIC
hold off

% 
% AIC
% min(AIC)

% h = ezcontour(@(x,y)pdf(obj{4},[x y]));
% axis auto


% fit model
n = 7;
for jj = 1:n
figure; hold on;
c = lines(n);
plot(T2.CSFI,T2.VFI,'*b')
[p,ErrorEst] = polyfit(T2.CSFI,T2.VFI,jj);

t2 = linspace(-0.4,1);
[y2,delta] = polyval(p,t2,ErrorEst);

% 95% c.i.
plot(t2,y2,'color',c(jj,:))
plot(t2,y2+2*delta,'r:')
plot(t2,y2-2*delta,'r:')

xlabel CSFI
ylabel VFI
title(sprintf('FIT %d function with CI',jj))
hold off;
end

%% lowess
figure; hold on;
yy2 = smooth(T2.CSFI,T2.VFI,0.3,'rloess');
% plot(yy2)
[xx,ind] = sort(T2.CSFI);
plot(xx,T2.VFI(ind),'b.',xx,yy2(ind),'r-')
xlabel 'CSFI %'
ylabel VFI
% legend({'POAG','2d','3d','4d'}

figure;
yy3 = smooth(T2.cpRNFL,T2.VFI,0.4,'rloess');
[xx,ind] = sort(T2.cpRNFL);
plot(T2.cpRNFL,T2.VFI,'b.',xx,yy3(ind),'r-')
xlabel 'cpRNFL'
ylabel VFI

 
%% staging
inds = T2.MD30_2>-6; 
ninds = NTG.MD30_2>-6;
pinds = POAG.MD30_2>-6;

% MD vs CSFI
figure;hold on;
plot(NTG.CSFI(ninds), NTG.MD30_2(ninds),'*r')
plot(POAG.CSFI(pinds),POAG.MD30_2(pinds),'*b')

xlabel CSFI
ylabel MD
title Early
legend({'NTG','POAG'})

% lsline
x = T2.CSFI(inds);
y = T2.MD30_2(inds);
mdl = fitlm(x,y)
[h,p] = corrcoef(x, y)


% VFI vs CSFI
figure;hold on;
% plot(x,y,'o','color',[0 0 0])%,'MarkerFaceColor',c(1,:))
plot(NTG.CSFI(ninds), NTG.VFI(ninds),'*r')
plot(POAG.CSFI(pinds),POAG.VFI(pinds),'*b')
xlabel CSFI
ylabel VFI
title Early
legend({'NTG','POAG'})

% title Early
% lsline
x = T2.CSFI(inds);
y = T2.VFI(inds);
mdl = fitlm(x,y)

[h,p] = corrcoef(x, y)

% RGC vs CSFI
figure;hold on;
% plot(x,y,'o','color',[0 0 0])%,'MarkerFaceColor',c(1,:))
plot(NTG.CSFI(ninds), NTG.RGC_OCT(ninds),'*r')
plot(POAG.CSFI(pinds),POAG.RGC_OCT(pinds),'*b')
xlabel CSFI
ylabel RGC_OCT
title Early
legend({'NTG','POAG'})

% title Early
% lsline
x = T2.CSFI(inds);
y = T2.RGC_OCT(inds);
mdl = fitlm(x,y)

[h,p] = corrcoef(x, y)

%% staging 1 subplot
inds = T2.MD30_2>-6; 
ninds = NTG.MD30_2>-6;
pinds = POAG.MD30_2>-6;

% lets see correlation between MD and others
names = fieldnames(T2);

% corr MD vs all
for kk = 1:length(names)
    try
        [h,p] = corr(T2.MD30_2(inds),T2.(kk)(inds));
        H(kk) = h;P(kk) = p;
    catch
        H(kk) = nan;P(kk) = nan; 
    end
end

names(abs(P)<0.05)

MD = table( names(abs(P)<0.05),H(abs(P)<0.05)',P(abs(P)<0.05)','VariableNames',{'variable','H','P'})

% lets see correlation between VFI and others
names = fieldnames(T2);

% corr VFI vs all
for kk = 1:length(names)
    try
        [h,p] = corr(T2.VFI(inds),T2.(kk)(inds));
        H(kk) = h;P(kk) = p;
    catch
        H(kk) = nan;P(kk) = nan; 
    end
end

names(abs(P)<0.05)

VFI = table( names(abs(P)<0.05),H(abs(P)<0.05)',P(abs(P)<0.05)','VariableNames',{'variable','H','P'})



% MD vs CSFI
figure;hold on;
subplot(1,3,1)
plot(NTG.CSFI(ninds), NTG.MD30_2(ninds),'*r')
plot(POAG.CSFI(pinds),POAG.MD30_2(pinds),'*b')

xlabel CSFI
ylabel MD
title Early
legend({'NTG','POAG'})

% lsline
x = T2.CSFI(inds);
y = T2.MD30_2(inds);
mdl = fitlm(x,y)
[h,p] = corrcoef(x, y)
text

% VFI vs CSFI
% figure;hold on;
subplot(1,3,2)
% plot(x,y,'o','color',[0 0 0])%,'MarkerFaceColor',c(1,:))
plot(NTG.CSFI(ninds), NTG.VFI(ninds),'*r')
plot(POAG.CSFI(pinds),POAG.VFI(pinds),'*b')
xlabel CSFI
ylabel VFI
title Early
legend({'NTG','POAG'})

% title Early
% lsline
x = T2.CSFI(inds);
y = T2.VFI(inds);
mdl = fitlm(x,y)

[h,p] = corrcoef(x, y)

% RGC vs CSFI
subplot(1,3,3)
% plot(x,y,'o','color',[0 0 0])%,'MarkerFaceColor',c(1,:))
plot(NTG.CSFI(ninds), NTG.RGC_OCT(ninds),'*r')
plot(POAG.CSFI(pinds),POAG.RGC_OCT(pinds),'*b')
xlabel CSFI
ylabel 'RGC OCT'
title Early
legend({'NTG','POAG'})

% title Early
% lsline
x = T2.CSFI(inds);
y = T2.RGC_OCT(inds);
mdl = fitlm(x,y)

[h,p] = corrcoef(x, y)


%% stage 2
inds = T2.MD30_2<-6 & T2.MD30_2>=-12 ; 
ninds = NTG.MD30_2<-6 & NTG.MD30_2>=-12 ; 
pinds = POAG.MD30_2<-6 & POAG.MD30_2>=-12 ; 


% lets see correlation between MD and others
names = fieldnames(T2);

% corr MD vs all
for kk = 1:length(names)
    try
        [h,p] = corr(T2.MD30_2(inds),T2.(kk)(inds));
        H(kk) = h;P(kk) = p;
    catch
        H(kk) = nan;P(kk) = nan; 
    end
end

names(abs(P)<0.05)

MD = table( names(abs(P)<0.05),H(abs(P)<0.05)',P(abs(P)<0.05)','VariableNames',{'variable','H','P'})

% lets see correlation between VFI and others
names = fieldnames(T2);

% corr VFI vs all
for kk = 1:length(names)
    try
        [h,p] = corr(T2.VFI(inds),T2.(kk)(inds));
        H(kk) = h;P(kk) = p;
    catch
        H(kk) = nan;P(kk) = nan; 
    end
end

names(abs(P)<0.05)

VFI = table( names(abs(P)<0.05),H(abs(P)<0.05)',P(abs(P)<0.05)','VariableNames',{'variable','H','P'})

%% continue

figure;hold on;

plot(NTG.CSFI(ninds), NTG.MD30_2(ninds),'*r')
plot(POAG.CSFI(pinds),POAG.MD30_2(pinds),'*b')

xlabel CSFI
ylabel MD
title Middle
legend({'NTG','POAG'})


x = T2.CSFI(inds);
y = T2.MD30_2(inds);

mdl = fitlm(x,y)

[h,p] = corrcoef(x, y)
title Middle


% VFI
y = T2.VFI(inds);

figure; hold on;
plot(NTG.CSFI(ninds), NTG.VFI(ninds),'*r')
plot(POAG.CSFI(pinds),POAG.VFI(pinds),'*b')

xlabel CSFI
ylabel VFI
title Middle
legend({'NTG','POAG'})

mdl = fitlm(x,y)

[h,p] = corr(x, y)

%% stage 3
inds = T2.MD30_2<-12 ; 
ninds = NTG.MD30_2<-12 ; 
pinds = POAG.MD30_2<-12 ; 

% lets see correlation between MD and others
names = fieldnames(T2);

% corr MD vs all
for kk = 1:length(names)
    try
        [h,p] = corr(T2.MD30_2(inds),T2.(kk)(inds));
        H(kk) = h;P(kk) = p;
    catch
        H(kk) = nan;P(kk) = nan; 
    end
end

names(abs(P)<0.05)

MD = table( names(abs(P)<0.05),H(abs(P)<0.05)',P(abs(P)<0.05)','VariableNames',{'variable','H','P'})

% lets see correlation between VFI and others
names = fieldnames(T2);

% corr VFI vs all
for kk = 1:length(names)
    try
        [h,p] = corr(T2.VFI(inds),T2.(kk)(inds));
        H(kk) = h;P(kk) = p;
    catch
        H(kk) = nan;P(kk) = nan; 
    end
end

names(abs(P)<0.05)

VFI = table( names(abs(P)<0.05),H(abs(P)<0.05)',P(abs(P)<0.05)','VariableNames',{'variable','H','P'})



%% MD
figure;hold on;

plot(NTG.CSFI(ninds), NTG.MD30_2(ninds),'*r')
plot(POAG.CSFI(pinds),POAG.MD30_2(pinds),'*b')

xlabel CSFI
ylabel MD
title Advanced
legend({'NTG','POAG'})


% VFI
% VFI
y = T2.VFI(inds);

figure; hold on;
plot(NTG.CSFI(ninds), NTG.VFI(ninds),'*r')
plot(POAG.CSFI(pinds),POAG.VFI(pinds),'*b')

xlabel CSFI
ylabel VFI
title Advanced
legend({'NTG','POAG'})

mdl = fitlm(x,y)

[h,p] = corrcoef(x, y)



%% linear regression fit

names = fieldnames(T2);

% corr MD vs all
for kk = 1:length(names)
    try
        [h,p] = corr(T2.MD30_2,T2.(kk));
        H(kk) = h;P(kk) = p;
    catch
        H(kk) = nan;P(kk) = nan; 
    end
end

names(abs(P)<0.05)

MD = table( names(abs(P)<0.05),H(abs(P)<0.05)',P(abs(P)<0.05)','VariableNames',{'variable','H','P'});

%% VFI
% corr VFI vs all
for kk = 1:length(names)
    try
        [h,p] = corr(T2.VFI,T2.(kk));
        H(kk) = h;P(kk) = p;
    catch
        H(kk) = nan;P(kk) = nan; 
    end
end

names(abs(P)<0.05)

VFI = table( names(abs(P)<0.05),H(abs(P)<0.05)',P(abs(P)<0.05)','VariableNames',{'variable','H','P'});


%% corr CSFI vs all
for kk = 1:length(names)
    try
        [h,p] = corr(T2.CSFI,T2.(kk));
        H(kk) = h;P(kk) = p;
    catch
        H(kk) = nan;P(kk) = nan; 
    end
end

names(abs(P)<0.05)

CSFI = table( names(abs(P)<0.05),H(abs(P)<0.05)',P(abs(P)<0.05)','VariableNames',{'variable','H','P'});



%%
tbl = table(T2.CSFI,T2.age,T2.AL,T2.cpRNFL,T2.FixLoss_pcnt,T2.Gender,T2.VFI,T2.MD30_2,...
    T2.RGC_HFA,T2.RGC_OCT,T2.SE,T2.wRGC,T2.AreaOD,...
    'VariableNames',{'CSFI','age','AL','cpRNFL','Fix_Loss','Gender','VFI','MD'...
    'RGC_HFA','RGC_OCT','SE','wRGC','AreaOD'});

lm = fitlm(tbl)

lm = fitlm(tbl,'MD ~ 1 + CSFI + age + cpRN')


%%
tbl = table(T2.CSFI,T2.age,T2.AL,T2.cpRNFL,T2.FixLoss_pcnt,T2.Gender,T2.VFI,T2.MD30_2,T2.AreaOD,...
    'VariableNames',{'CSFI','age','AL','cpRNFL','Fix_Loss','Gender','VFI','MD','AreaOD'});

lm = fitlm(tbl)

lm = fitlm(tbl,'MD ~ 1 + CSFI + age + cpRNFL + Fix_Loss + Gender + VFI')

lm = fitlm(tbl,'MD ~ 1 + CSFI + age + cpRNFL')

%

lm = fitlm(tbl,'CSFI ~ 1 + VFI + age + cpRNFL')

lm = fitlm(tbl,'CSFI ~ 1 + MD + age + cpRNFL')

%% staging
inds = T2.MD30_2>-6; 
Early = T2(inds,:);

inds = T2.MD30_2<-6 & T2.MD30_2>=-12 ; 
Middle = T2(inds,:);

inds = T2.MD30_2<-12 ; 
End  = T2(inds,:);

%% Early
tbl = table(Early.CSFI,Early.age,Early.AL,Early.cpRNFL,Early.FixLoss_pcnt,...
    Early.Gender,Early.VFI,Early.MD30_2,...
    'VariableNames',{'CSFI','age','AL','cpRNFL','Fix_Loss','Gender','VFI','MD'});

lm = fitlm(tbl,'CSFI ~ 1 + VFI + age + cpRNFL')
% anova(lm,'summary')

lm = fitlm(tbl,'CSFI ~ 1 + MD + age + cpRNFL')

%% Middle
tbl = table(Middle.CSFI,Middle.age,Middle.AL,Middle.cpRNFL,Middle.FixLoss_pcnt,...
    Middle.Gender,Middle.VFI,Middle.MD30_2,...
    'VariableNames',{'CSFI','age','AL','cpRNFL','Fix_Loss','Gender','VFI','MD'});

lm = fitlm(tbl,'CSFI ~ 1 + VFI + age + cpRNFL')
lm = fitlm(tbl,'CSFI ~ 1 + MD + age + cpRNFL')

%% Advanced
tbl = table(End.CSFI,End.age,End.AL,End.cpRNFL,End.FixLoss_pcnt,...
    End.Gender,End.VFI,End.MD30_2,...
    'VariableNames',{'CSFI','age','AL','cpRNFL','Fix_Loss','Gender','VFI','MD'});

lm = fitlm(tbl,'CSFI ~ 1 + VFI + age + cpRNFL')
lm = fitlm(tbl,'CSFI ~ 1 + MD + age + cpRNFL')

%% Bland altman
addpath(genpath('BlandAltman'))
label ={'RGC HFA','RGC OCT'};
[rpc, fig, stats] = BlandAltman(T2.RGC_OCT,T2.RGC_HFA,label);

% set(gca , 'FontSize',18)
%% by myself
figure; hold on;
plot(T2.RGC_OCT,T2.RGC_HFA, 'o')
lsline
axis equal
xlabel 'RGC OCT'
ylabel 'RGC HFA'
set(gca, 'FontSize',18, 'XLim',get(gca,'YLim'))
title 'Correlation plot'

[x1] = get(gca,'YLim');
x = linspace( x1(1), x1(2) );

plot(x, x, '--k')
title 'Correlation plot'

%%
[h,p,ci,stats]= ttest(T2.RGC_OCT,T2.RGC_HFA)


% [h,p,ci,stats]= ttest(T2.RGC_OCT-T2.RGC_HFA ,mean([T2.RGC_OCT,T2.RGC_HFA],2))


t = stats.r*sqrt( (stats.N -2)/(1-stats.r2));

%% student t distribution

x = [-10:.1:10];
y1 = tpdf(x,stats.N -2);   % For nu = 5
% y2 = tpdf(x,25);  % For nu = 25
% y3 = tpdf(x,50);  % For nu = 50

figure;
plot(x,y1,'Color','black','LineStyle','-')
hold on
% plot(x,y2,'Color','red','LineStyle','-.')
% plot(x,y3,'Color','blue','LineStyle','--')
% legend({'nu = 5','nu = 25','nu = 50'})
hold off

