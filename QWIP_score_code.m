%%% load Rrs variables
load('rrs_merged_dataset.mat');

Rrs_vis = rrsT(:,52:352);
wave = (400:700);
Rrs_492 = rrsT(:,144);
Rrs_665 = rrsT(:,317);
RRS = table2array(Rrs_vis);

%%% Size up the array and calcualte AVW
[m,n] = size(RRS);
AVW = zeros(m,1);
for i = 1:m
    AVW(i) = sum(RRS(i,:))./sum(RRS(i,:)./wave);
end


%%% Now we calculate the Normalized Difference Index (NDI)
%%% (Rrs_665-Rrs_492)/((Rrs_665+Rrs_492)

[match_492,index_492] = min(abs(wave-492));
[match_665,index_665] = min(abs(wave-665));
avw_poly = (400:640);
NDI = (RRS(:,index_665) - RRS(:,index_492))./(RRS(:,index_665) + RRS(:,index_492));
p=[   -8.399884740300151e-09     1.715532100780679e-05    -1.301670056641901e-02 ...
     4.357837742180596e+00    -5.449532021524279e+02 ];
fit1 = polyval(p,avw_poly);

%%% Just generating an array here of the "predicted" NDI based on the AVW
%%% image values. We then subtract this from the actual NDI that we
%%% calculate above, and get an pixel by pixel QWIP score.
NDI_pred = (p(1)*AVW.^4 + p(2)*AVW.^3 + p(3)*AVW.^2 + p(4)*AVW.^1 + p(5));
QWIP_score = NDI - NDI_pred;
abs_QWIP_score = abs(QWIP_score);
QWIP_flag = abs_QWIP_score>=0.2;


%Water Type

Rrs_665b = table2array(Rrs_vis(:,266));
Rrs_560b = table2array(Rrs_vis(:,161));
Rrs_492b = table2array(Rrs_vis(:,93));

Step1 = (Rrs_665b>Rrs_560b);
Step2 = (Rrs_665b>0.025);
Step3 = (Rrs_560b<Rrs_492b);
ind_600A = (Step1>0 | Step2>0);
ind_500A = (Step1<1 & Step2<1)& Step3<1;
ind_400A = (Step1<1 & Step2<1)& Step3>0;


%%% Generate figure to show QCI index relative to AVW
fit1a = fit1 + 0.1;
fit1b = fit1 - 0.1;
fit2a = fit1 + 0.2;
fit2b = fit1 - 0.2;
fit3a = fit1 + 0.3;
fit3b = fit1 - 0.3;
fit4a = fit1 + 0.4;
fit4b = fit1 - 0.4;
figure(2)
plot1 = plot(AVW(ind_500A),NDI(ind_500A),'og','MarkerSize',1);hold on;
plot(AVW(ind_400A),NDI(ind_400A),'ob','MarkerSize',1);
plot(AVW(ind_600A),NDI(ind_600A),'or','MarkerSize',1);
plot(avw_poly,fit1,'-k','LineWidth',2);
plot(avw_poly,fit1a,'--g','LineWidth',2);
plot(avw_poly,fit1b,'--g','LineWidth',2);
plot(avw_poly,fit2a,'--','LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
plot(avw_poly,fit2b,'--','LineWidth',2,'Color',[0.9290 0.6940 0.1250]);
plot(avw_poly,fit3a,'--','LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
plot(avw_poly,fit3b,'--','LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
plot(avw_poly,fit4a,'-r','LineWidth',2);
plot(avw_poly,fit4b,'-r','LineWidth',2);
xlabel('AVW (nm)','FontSize',16);
ylabel(['NDI (' num2str(wave(index_492)) ',' num2str(wave(index_665)) ')'],'FontSize',16);
ylim([-2.5 2]);
xlim([440 630]);
%%% Forgive my sloppiness here, I often do this little cheat
%%% just to get my MATLAB legends to cooperate...
LH(1) = plot(nan, nan, '--g');
L{1} = 'QWIP ± 0.1';
LH(2) = plot(nan, nan, '--','Color',[0.9290 0.6940 0.125]);
L{2} = 'QWIP ± 0.2';
LH(3) = plot(nan, nan, '--','Color',[0.8500 0.3250 0.0980]);
L{3} = 'QWIP ± 0.3';
LH(4) = plot(nan, nan, '-r');
L{4} = 'QWIP ± 0.4';
leg1=legend(LH, L,'Location','northwest','FontSize',12);hold on;
hold off;


