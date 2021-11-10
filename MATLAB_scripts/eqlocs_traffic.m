clear 
h1=worldmap([53.58 53.9],[-3.3 -2.65]);   
%land = shaperead('landareas.shp','UseGeoCoords', true);
%geoshow(land)
xlabel('Latitude');
ylabel('Longitude');

latlim = [54  52];
lonlim = [-3.0 -2.6];
s1=shaperead('GSHHS_h_L1','UseGeoCoords', true);
s2=shaperead('GSHHS_h_L2','UseGeoCoords', true);
s3=shaperead('GSHHS_h_L3','UseGeoCoords', true);
s4=shaperead('GSHHS_h_L4','UseGeoCoords', true);
s5=shaperead('GSHHS_h_L5','UseGeoCoords', true);
s6=shaperead('GSHHS_h_L6','UseGeoCoords', true);
geoshow(s1)
geoshow(s2)
geoshow(s3)
geoshow(s4)
geoshow(s5)
geoshow(s6)
%setm(h1,'FFaceColor','w') % set the frame fill to white
%places = {'Blackpool';'Preston';'Liverpool';'Manchester'}; 
%latp = [53.8169 53.7629 53.4057 53.483959];   %latitude locations of the 3 stations
%lonp = [-3.0371 -2.7035 -2.9931 -2.244644];   %longitude locations of the 3 stations
%l1=textm(latp,lonp,places,'fontsize',14,'color','k'); %the textm function from the mapping toolbox allows text to be inserted on the map using the 3 variables above

placem = ('*');
latm = (53.7871);
lonm = (-2.951620);
l2=textm(latm,lonm,placem,'fontsize',28,'color','r');

placewrite=('fracking');
lonm2= (-2.94);
l3=textm(latm,lonm2,placewrite,'fontsize',16,'color','r');

stationsactive={'L2';'L4';'L6';'L9'};
latsa =[53.84 53.76 53.69 53.78];
lonsa =[-2.84 -2.94 -2.88 -3.04];
l4=textm(latsa,lonsa,stationsactive,'fontsize',12,'color','k');

latr2=(53.8);
lonr2=(-3.28);
txt2= {'University Seismometers';'locations accurate to ± 1km'};
l5=textm(latr2,lonr2,txt2,'fontsize',16,'color','k');

txt3={'Earthquakes over 0.5 magnitude: x'};
l6=textm(53.75, lonr2,txt3,'fontsize',16,'color','r');
txt4={'Earthquakes of magnitudes 0 to 0.5: x'};
l7=textm(53.74,lonr2,txt4,'fontsize',16,'color','b');

txt5={'Earthquakes less than 0 magnitude: x'};
l8=textm(53.73,lonr2,txt5,'fontsize',16,'color','g');

bgs_locs_lat=[53.787 53.785 53.784 53.785 53.788 53.787 53.786 53.786 53.787 53.789 53.788 53.789 53.786 53.788 53.788 53.787];
bgs_locs_lon=[-2.977 -2.971 -2.970 -2.967 -2.965 -2.968 -2.963 -2.966 -2.962 -2.963 -2.963 -2.962 -2.963 -2.964 -2.963 -2.963];
my_locs_lat=[53.725 53.775 53.767 53.766 53.694 53.767 53.768 53.768 53.767 53.765 53.767 53.765 53.768 53.780 53.735 53.765];
my_locs_lon=[-2.657 -2.965 -2.923 -2.918 -2.461 -2.930 -2.919 -2.919 -2.924 -2.910 -2.922 -2.910 -2.929 -3.000 -2.725 -2.912];
onlymy_locs_lat=[53.768 53.767 53.767];
onlymy_locs_lon=[-2.925 -2.924 -2.920];

redquakes_lat=[53.765 53.765];
redquakes_lon=[-2.9103 -2.911];
orangequakes_lat=[53.725 53.767 53.694 53.767 53.768 53.735 53.765];
orangequakes_lon=[-2.657 -2.923 -2.461 -2.93 -2.919 -2.725 -2.912];
greenquakes_lat=[53.775 53.766 53.768 53.767 53.767 53.768 53.78 53.768 53.767 53.767];
greenquakes_lon=[-2.965 -2.918 -2.919 -2.924 -2.922 -2.929 -3 -2.925 -2.924 -2.92];

labeltxt2='x';
for i=1:2
    
    label(i)=textm(redquakes_lat(i), redquakes_lon(i),labeltxt2,'fontsize',19,'color','r');
end

for i=1:7
    labelzz(i)=textm(orangequakes_lat(i),orangequakes_lon(i),labeltxt2,'fontsize',14,'color','b');
end

for i=1:10
    labelzzs(i)=textm(greenquakes_lat(i),greenquakes_lon(i),labeltxt2,'fontsize',12,'color','g');
end
scaleruler on
northarrow('latitude', 53.84, 'longitude', -3.25,'scaleratio', 0.11);
t1=title(h1,'Earthquake Locations and Magnitudes','fontsize',18); %adding a title to figure





