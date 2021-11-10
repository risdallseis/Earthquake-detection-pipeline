

%%
%STA/LTA L004 24/10/18 100SR 3.0 thresh
clear
stalta=[24	4	23
25	1	89
26	3	82
27	2	85
29	5	168
];

d=stalta(:,1);% = DateOctober2018;
c=stalta(:,2) ;%= Correct Detections;
f=stalta(:,3) ;%= Incorrect Detections;
figure
plot(d,c)
hold on
plot(d,f)



%CC L004 24/10/18 template: 1.1M  0.44 thresh
cc=[24	2	4
    25	0	0
    26	2	2
    27	3	3
    29	4	6];

cd=cc(:,1);% = DateOctober2018;
cf=cc(:,2) ;%= Correct Detections;
ff=cc(:,3) ;%= Incorrect Detections;

plot(cd,cf)

plot(cd,ff)

%plot either bar graphs with 2bars for each day or confidence graph for
%percentage chance of correct detection