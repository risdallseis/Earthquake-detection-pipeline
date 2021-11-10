clear 
h1=worldmap([53.58 53.9],[-3.3 -2.65]);   
%land = shaperead('landareas.shp','UseGeoCoords', true);


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
places = {'Blackpool';'Preston';'Liverpool';'Manchester'};   %creating a cell array with the 3 stations in order to label them on the map
placem = ('x');
placewrite=('drilling');
stationsactive={'L2';'L4';'L6';'L7';'L8';'L9'};
stationsoff={'L1';'L3';'L5'};
latso=[53.82 53.85 53.77];
lonso=[-2.93 -2.80 -2.87];
txt= {'Offline stations';'locations accurate to ± 1km'};
txt2= {'Online stations';'locations accurate to ± 1km'};
latr=(53.75);
lonr=(-3.25);
latr2=(53.8);
lonr2=(-3.25);
latp = [53.8169 53.7629 53.4057 53.483959];   %latitude locations of the 3 stations
lonp = [-3.0371 -2.7035 -2.9931 -2.244644];   %longitude locations of the 3 stations
latm = (53.7871);
lonm = (-2.951620);
lonm2= (-2.93);
latsa =[53.84 53.76 53.69 53.74 53.60 53.78];
lonsa =[-2.84 -2.94 -2.88 -2.73 -3.04 -3.04];
l1=textm(latp,lonp,places,'fontsize',14,'color','k'); %the textm function from the mapping toolbox allows text to be inserted on the map using the 3 variables above
l2=textm(latsa,lonsa,stationsactive,'fontsize',16,'color','b');
l3=textm(latm,lonm,placem,'fontsize',18,'color','k');
l4=textm(latr,lonr,txt,'fontsize',16,'color','r');
l5=textm(latso,lonso,stationsoff,'fontsize',16,'color','r');
l6=textm(latr2,lonr2,txt2,'fontsize',16,'color','b');
l7=textm(latm,lonm2,placewrite,'fontsize',16,'color','k');
scaleruler on
northarrow('latitude', 53.84, 'longitude', -3.25,'scaleratio', 0.11);
t1=title(h1,'Locations of drilling and seismometer stations','fontsize',18); %adding a title to figure




