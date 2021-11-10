mag1=[2 4 10]; 
mag8=[2 5 8];
mag5=[1 4 4];
mags=[mag1;mag8;mag5];
temp1_low=[16,15,28];
temp1_high=[5,26,0];
temp8_low=[15,16,22];
temp8_high=[5,26,1];
temp5_low=[7,24,47];
temp5_high=[3,28,0];

temps=[temp1_low;temp1_high;temp8_low;temp8_high;temp5_low;temp5_high];

bgs_detected_z=[3 3.5 3.4 3.5 2.4 2.3 2.2 2.3 2.2 2.5 2.4 2.9 2.4 2.3 2.5 2.5];
bgs_detected_m=[0.4 0.5 0.4 0.1 0.3 -0.1 -0.6 -0.1 -0.3 0.8 0.0 1.1 -0.4 0.0 0.5 0.7];
bgs_nodetected_z=[3.2 3.7 2.9 2.2 2.9 2.2 2.2 2.4 2.0 2.3 2.2 2.2];
bgs_nodetected_m=[-0.2 -0.3 -0.1 0.2 0.8 -0.2 -0.4 -0.2 -0.5 -0.4 -0.5 -0.5];

bgs_detected=[bgs_detected_z;bgs_detected_m];
bgs_nodetected=[bgs_nodetected_z;bgs_nodetected_m];

numbervec=[1:16];
mymags=[0.5 -0.04 0.2 -0.9 0 0 0.1 -0.1 -0.4 0.6 0 1 -0.7 -0.2 0.4 0.4];

    




%%
figure

plot(numbervec,mymags,'b--x','MarkerSize',16);
hold on
plot(numbervec,bgs_detected_m,'k--x','MarkerSize',16);
title('Magnitude Calculation Comparison');
xlabel('Earthquake No.');
ylabel('Local Magnitude');
labels1{1}='This Study'; labels1{2}='BGS';
legend(labels1)
set(gca,'FontSize',20)
%%
figure
plot(bgs_detected(1,:),bgs_detected(2,:),'k*','MarkerSize',18);
hold on
plot(bgs_nodetected(1,:),bgs_nodetected(2,:),'rx','MarkerSize',18);
labels{1}='Detected'; labels{2}='Not Detected';
legend(labels)
set(gca,'FontSize',20)
title('Algorithm Effectivity Using BGS Detections')
xlabel('Depth (km)');
ylabel('Local Magnitude');
%%

figure
y=bar(mags,'grouped');
title('Cross-Correlation Detected Earthquakes Sorted By Magnitude')
y(1).FaceColor = 'r';
y(2).FaceColor = 'y';
y(3).FaceColor = 'g';
l{1}='Greater than 0.5M'; l{2}='0M to 0.5M'; l{3}='Less than 0M'; 
ylabel('Amount')
set(gca,'FontSize',20)

legend(y,l);
name = {'1.1M Template';'0.8M Template';'0.5M Template'};
set(gca,'xticklabel',name)
%%
figure
d=bar(temps,'grouped');
title('Cross-Correlation Detections')
d(1).FaceColor = 'g';
d(2).FaceColor = 'b';
d(3).FaceColor = 'r';
g{1}='True +'; g{2}='False +'; g{3}='True -'; 
ylabel('Amount')
set(gca,'FontSize',20)
legend(d,g);
name = {'1.1M (0.35)';'1.1M (0.45)';'0.8M (0.35)';'0.8M (0.45)';'0.5M (0.35)';'0.5M (0.45)'};
xlabel('Template and Trigger Coefficient')
set(gca,'xticklabel',name)