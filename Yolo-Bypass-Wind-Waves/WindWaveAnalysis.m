%Code to run the wind analysis, calculate wave heights and run-up at the
%YBEL.
clear
close all
%% Calculate percentage of wind directions for Upper Yolo Bypass
% dir = csvread('e:\CBEC\Projects\18_1011_YB_Wind_Wave\MBK Review Memo\WindDirs.csv');
% bins = 0:22.5:360;
% counts = histcounts(dir,bins);
% f1 = figure(1);
% subplot(121)
% hist = histogram('Categories',{'N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW',...
%     'WSW','W','WNW','NW','NNW'},'bincounts',counts);
% subplot(122)
% theta = dir*pi/180;
% th = theta;
% pr = rose(th,16);
% set(gca,'View',[-90 90],'YDir','reverse');
% 
% %percentage NNW:
% pct = (counts(15)./length(dir))*100;
% fprintf('Percentage of year when wind blows NNW: %0.2f\n',pct)

% %% Import water level data from the DWSE report
% dir = 'e:\CBEC\Projects\19-1034_YBEL\Wind_Wave_Analysis\';
% WSGRR = csvread([dir 'Yolo_DWSE_Plate2_Data_WSGRR.csv']);
% DWSE = csvread([dir 'Yolo_DWSE_Plate2_Data_DWSE.csv']);
% WS57 = csvread([dir 'Yolo_DWSE_Plate2_Data_1957WSProfile.csv']);
% f2 = figure(2);hold on
% p(1) = plot(WSGRR(:,2),WSGRR(:,1),'g','linewidth',1.5);
% p(2) = plot(DWSE(:,2),DWSE(:,1),'b','linewidth',1.5);
% p(3) = plot(WS57(:,2),WS57(:,1),'r','linewidth',1.5);
% grid on
% box on
% xlabel('Levee Station')
% ylabel('Stage, ft [NAVD88]')
% title('Yolo Bypass East Levee Design Profile')
% 
% %% Load station data, find corresponding water level elevations in DWSE
% fid = fopen([dir 'WindWaveStations_v2.csv']);
% wwstations = textscan(fid,'%f%s%f%f%f','Delimiter',',');
% stations = wwstations{5}';
% stationx = wwstations{3}';stationy = wwstations{4}';
% [~,idx] = min(abs(stations-WSGRR(:,2)));
% wse_wsgrr = WSGRR(idx,1);
% [~,idx] = min(abs(stations-DWSE(:,2)));
% wse_dwse = DWSE(idx,1);
% [~,idx] = min(abs(stations-WS57(:,2)));
% wse_ws57 = WS57(idx,1);
% plot(ones(10,1).*stations,ones(1,6).*linspace(27,33,10)','--k')
% 
% legend(p',{'WS GRR (1/200 ACE)';'DWSE (1/200 ACE)';'1957 W.S. Profile'},...
%     'location','eastoutside')
%% Load topobathy data from the cbec YBEL survey (Nov 2019)
% dir = 'e:\CBEC\Projects\19-1034_YBEL\Data\GIS\';
% fid = fopen([dir 'YBEL_cbec_survey.xyz']);disp('Loading topobathy survey data...')
% topobath = textscan(fid,'%f%f%f','Delimiter',' ');
% tbx = topobath{1};tby = topobath{2}; tbz = topobath{3};
% xv = linspace(min(tbx), max(tbx), 20);
% yv = linspace(min(tby), max(tby), 20);disp('Gridding topobathy data...')
% [X,Y] = meshgrid(xv, yv);
% Z = griddata(tbx,tby,tbz,X,Y);
% f3 = figure(3);
% surf(X,Y,Z);
% % Extract depths from the topobathy data for the same station locations
% [~,idx] = min(abs(stationy-tby));
% stationz = tbz(idx);
% hold on
% plot3(stationx,stationy,stationz+2,'marker','o',...
%     'linestyle','none',...
%     'markersize',10,...
%     'markerfacecolor','y',...
%     'markeredgecolor','k')
% text(stationx,stationy,stationz+12,wwstations{2})
% % Plot WSEs from the former plot as lines on the surface
% xs = ones(size(WSGRR))*mean(stationx);
% ys = linspace(min(tby),max(tby),length(WSGRR));
% p1 = plot3(xs,ys,WSGRR(:,1),'g','linewidth',1.5);
% xs = ones(size(DWSE))*mean(stationx);
% ys = linspace(min(tby),max(tby),length(DWSE));
% p2 = plot3(xs,ys,DWSE(:,1),'b','linewidth',1.5);
% xs = ones(size(WS57))*mean(stationx);
% ys = linspace(min(tby),max(tby),length(WS57));
% p3 = plot3(xs,ys,WS57(:,1),'r','linewidth',1.5);
% grid on
% view(-64,9)
% set(gca, 'ZLim',[0 40])
% shading interp
% xlabel('Easting (ft)'),ylabel('Northing (ft)'),zlabel('Elevation (ft) [NAVD88]')
% legend([p1(1)' p2(1)' p3(1)'],{'WS GRR (1/200 ACE)';'DWSE (1/200 ACE)';'1957 W.S. Profile'},...
%     'location','southeast')

%% Calculate wind-wave properties and wave-runup percentage
% Wind data (from wind_wave_eqns.xlsx)
% page 249 of SPM (document page 3-55)
% "3-39"	gH/Ua^2 = 0.283 * tanh [ 0.530 (gd/Ua^2)^(3/4) ] * tanh { (0.00565*(gF/Ua^2)^(1/2))/( tanh [ 0.530* (gd/Ua^2)^(3/4) ]) }
% "3-40"	gT/Ua = 7.54 * tanh [ 0.833 (gd/Ua^2)^(3/8) ] * tanh { (0.0379*(gF/Ua^2)^(1/3))/( tanh [ 0.833* (gd/Ua^2)^(3/8) ]) }
% 3-28a (m/s)	Ua=0.71*U^1.23
%
%   W          WNW      NW        NNW
U = [17.479264 26.15184 24.855424 27.71648]; %50-year extreme wind speed in m/s
h = [20.2 24.3 26.6 29.0]; %WSE from USACE HEC-RAS model 
g = 9.81; %grav constant
Ua = 0.71*U.^1.23; %wind speed at 10m
%Load fetch data from GIS analysis
dir = 'e:\CBEC\Projects\19-1034_YBEL\Wind_Wave_Analysis\';
fid = fopen([dir 'YBEL_FetchLengths_v4.csv']);
fetch = textscan(fid,'%f%s%s%f%f%f%f%f','Delimiter',',','headerlines',1);
f = fetch{4}; %just the fetch lengths
depth = [fetch{5:8}]; %avg water depths along transects
station = fetch{2};stationid = unique(station);
wsname = {'twoyr';'fiveyr';'tenyr';'twentyfiveyr'}; %data structure field names
fn = {'id';'StationID';'WindDirection';'FetchLength';'WindSpeed';...
    'WaterDepth';'WaveHeight';'WavePeriod';'WindSetUp';'WaveRunUp';...
    'TWL';'Taub';'Dn50';'Tauc'}; %field names for the data table
Ws = zeros(length(f),1); %wind setup elevation
Hs = zeros(length(f),1); %wave height
Tp = zeros(length(f),1); %wave period
R2pct = zeros(length(f),1); %2 percent wave runup elevation
Dn50 = zeros(length(f),1); %required median nominal riprap size given wave conditions
Taub = zeros(length(f),1); %wave shear stress due to breaking against levee
Tauc = zeros(length(f),1); %critical shear stress of levee material based on Dn50
for i = 1:length(wsname)
    fprintf('Running analysis for %s water level scenario\n',wsname{i})
    for j = 1:length(stationid)
        idx = strcmp(station,stationid(j));
        F = (f(idx)./3.281)'; %convert freedom units to meters
        d =  (depth(idx,i)./3.281)'; %convert freedom units to meters
        %Wind set-up
        S = ((U.^2).*(F/1000))./(62000*d); %units in m
        
        %Wave height and period: eqns. from USACE SPM 3-39
        gHUasqd = 0.283*tanh(0.53*((g*d)./(Ua.^2)).^(3/4))*tanh((0.00565*((g*F)./(Ua.^2)).^(1/2))/tanh(0.530*((g*d)./(Ua.^2)).^(3/4)));
        gTUa = 7.54*tanh(0.833*((g*d)./(Ua.^2)).^(3/8))*tanh((0.0379*((g*F)./(Ua.^2)).^(1/3))/tanh(0.833*((g*d)./(Ua.^2)).^(3/8)));
        H = (gHUasqd.*(Ua.^2))/g; %wave height in m
        T = (gTUa.*Ua)/g; %wave period in s
        
        %Wave properties
        w = 2*pi./T; %omega
        kh = qkhf(w,d); %k from dispersion relation--code from USGS
        k = kh./d;
        %the following eqns use imperial instead of metric
        dfu = d.*3.281; %ft
        Hfu = H.*3.281; %ft
        kfu = k./3.281; %ft^-1
        lambdafu = 2*pi./kfu; %ft

        %calculate Iribarren number (surf similarity parameter)
        a = 18.4; %slope of levee (arctan(levee height/levee length from toe to top))
        Xi = tand(a)./(sqrt(Hfu./lambdafu));
        A = 1.6;C = 0; %CEM Table V1-5-2 for wind/wave conditions during this study
        gamr = 1; %CEM Table V1-5-3 assuming two rock layers with H/D = 1.5 - 6
        gamb = 1; %Assumed no bermed fetch
        gamh = 1; %Assumed no shallow wave influence
        gambeta = 1; %Assumed waves only approach along fetch cross-section line
        r2pct = Hfu.*(A.*Xi+C)*gamr*gamb*gamh*gambeta;
        
        %estimate the riprap size needed to withstand waves (EM-1110-2-1100, eq. V1-5-85) 
        gamr = 160; %lbf/ft^3 (stone)
        gamw = 62.4; %lbf/ft^3 (fresh water)
        Kd = 16; %V1-5-81
        delta = (gamr/gamw)-1;
        W50 = (gamr*(Hfu.^3))./(Kd*(delta^3)*cotd(a)); %weight of median stone
        dn50 = (W50/gamr).^(1/3); %median nominal size of riprap required to withstand waves (ft)
        
        %compute hydraulic bed shear stress due to wave breaking (URS,
        %2011)
        B = 1; %empirical coefficient = 1
        gfu = 32.2; %ft/s^2
        Hmax = (0.88./kfu).*tanh(0.78*((kfu.*dfu)./0.88)); %wave height at breaking point
        eps = 0.25*(gamw*(1./T)).*(((B.*Hmax).^3)./dfu); %wave energy dissipation rate
        Cg = 0.5*(sqrt((gfu./kfu).*tanh(2*pi*(dfu./lambdafu)))).*(1+((2*(kfu.*dfu))./sinh(2*(kfu.*dfu)))); %wave group speed
        e = 0.075; %wave breaking shear streess efficiency (URS, 2007)
        taub = e*eps./Cg; %shear stress due to wave breaking (lbf/ft^2)
        
        %compute material critical shear stress
        kl = sqrt(1-((sind(a)^2)/(sind(40)^2))); %eq V1-5-129
        SP = 0.047; %shields parameter for riprap
        tauc = kl*SP*((gamr)-(gamw))*dn50; %eq V1-5-128
        
        %Save out data
        Ws(idx) = S.*3.281;
        Hs(idx) = Hfu;
        Tp(idx) = T;
        R2pct(idx) = r2pct;
        Dn50(idx) = dn50.*12; %convert ft to in
        Taub(idx) = taub;
        Tauc(idx) = tauc;
    end
    Umph = repmat(U'.*2.2373,4,1); %convert wind speed to mph
    wse = repmat(h(i),length(f),1); %station depth (ft)
    int = {fetch{1} fetch{2} fetch{3} fetch{4} Umph wse Hs Tp Ws R2pct wse+Ws+R2pct Taub Dn50 Tauc};
    data.(wsname{i}) = cell2struct(int,fn,2);
end
%display data as tables for easy copying to excel
disp('Two-year RI data')
struct2table(data.twoyr)
disp('Five-year RI data')
struct2table(data.fiveyr)
disp('Ten-year RI data')
struct2table(data.tenyr)
disp('Twenty Five-year RI data')
struct2table(data.twentyfiveyr)





