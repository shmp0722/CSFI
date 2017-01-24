function CSFI_Tanabe

%% Load whole data

T = readtable('CSFI_data.csv');
C = readtable('Normal.csv');
N = readtable('NTG.csv');
P = readtable('POAG.csv');

%%
figure; hold on;
c = jet(4);


plot(C.age,C.MD_24_2,'o','color',c(1,:))
plot(N.age,N.MD_24_2,'o','color',c(2,:))
plot(P.age,P.MD_24_2,'o','color',c(3,:))


xlabel age
ylabel MD

%% CSFI vs VFI
figure; hold on;
c = jet(4);

plot(C.CSFI_rate,C.VFI_rate,'o','color',c(1,:),'MarkerFaceColor',c(1,:))
plot(N.CSFI_rate,N.VFI_rate,'o','color',c(2,:),'MarkerFaceColor',c(2,:))
plot(P.CSFI_rate,P.VFI_rate,'o','color',c(3,:),'MarkerFaceColor',c(3,:))

[h,p] = ttest2(C.CSFI_rate,N.CSFI_rate)
[h,p] = ttest2(C.CSFI_rate,P.CSFI_rate)
[h,p] = ttest2(N.CSFI_rate,P.CSFI_rate)


xlabel 'CSFI %'
ylabel 'VFI %'

legend({'Healthy','NTG','POAG'})

%%
ttest(t(:,1),t(:,2));
corr(t(:,1),t(:,2));

figure; hold on;
plot(t(:,1),t(:,2)*100,'o')
xlabel MD
ylabel CSFI

figure; hold on;
plot(t(:,2)*100,t(:,1),'o')
xlabel CSFI
ylabel MD
%% fit 
p = polyfit(t(:,2)*100,t(:,1),2);

t2 = -40:1:100;
y2 = polyval(p,t2);

figure;hold on;
plot(t(:,2)*100,t(:,1),'o',t2,y2)
xlabel CSFI
ylabel MD


%%
cd /media/USB_HDD1/CSFI
ds = csvread('CSFI_data.csv');


%% 
figure; hold on;
plot(age ,VFI_rate,'o')
xlabel Age
ylabel VFI

%%
figure; hold on;
plot(CSFI_rate ,VFI_rate,'o')
xlabel 'CSFI %'
ylabel VFI
[p,s,mu] = polyfit(CSFI_rate ,VFI_rate,2)

%% Soukan VFI vs 
[r,p] = corr(CSFI_rate ,VFI_rate)

[r,p] = corr(RGC_HFA ,VFI_rate)

[r,p] = corr(RGC_OCT ,VFI_rate)

[r,p] = corr(wRGC ,VFI_rate)

[r,p] = corr(wRGC ,VFI_rate)


%%










