%% load data 

T = readtable('Latest20170208.xlsx');

%% remove subjects HFA reliability is low
rows =  T.FP< .15 & T.FN<.33 & T.FixLoss_pcnt<.2;
% rows =  T.FP< .15 & T.FixLoss_pcnt<.2;

T2 = T(rows,:);

%% exponetial to linear

T2.MD30_linear = 10.^T2.MD30_2;

%%
figure; hold on;
% subplot(1,2,1)
scatter(T2.wRGC,T2.MD30_2,40,'filled')
xlabel RGC
ylabel MD-dB

figure; hold on;
% subplot(1,2,2)
scatter(T2.wRGC,T2.MD30_linear,40,'filled')
xlabel RGC
ylabel MD-linear

%% rgcOCT- rgcOCT

figure; hold on;
plot(T.RGC_OCT,T.RGC_HFA,'.')
xlabel RGC_OCT
ylabel RGC_OCT
axis equal
get(gca,'ylim')
set(gca,'xlim',[0 2000],'ylim',[0, 2000])

x = 0:2000;
line(x,x)


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


%% corr MD linear vs all
for kk = 1:length(names)
    try
        [h,p] = corr(T2.MD30_linear(inds),T2.(kk)(inds));
        H(kk) = h;P(kk) = p;
    catch
        H(kk) = nan;P(kk) = nan; 
    end
end

names(abs(P)<0.05)

MD_L = table( names(abs(P)<0.05),H(abs(P)<0.05)',P(abs(P)<0.05)','VariableNames',{'variable','H','P'})