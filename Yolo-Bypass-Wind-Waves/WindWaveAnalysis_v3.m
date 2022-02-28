%Code to run the wind analysis, calculate wave heights and run-up at the
%YBEL. This is version 2 of this script, using only NW, W, and SW winds
clear

%% Calculate wind-wave properties and wave-runup percentage
% Wind data (from WindWaveTables.xlsx) COMMENT OUT IF YOU WANT TO RUN 1HR
% or 10HR data!
%
%   fetch-limited wind duration
%   NW    W    SW
% U = [46 57.6 70.5]; %50-year fetch-limited wind speed in mph
%
%
%   1-hr wind duration
%   NW    W    SW
U = [46.0 57.5 70.5]; %50-year 1-hr extreme wind speed in mph
%
%   10-hr wind duration
%   NW    W    SW
% U = [39.1 48.9 59.9]; %50-year 10-hrs extreme wind speed in mph
%
%
gammar = [0.6 0.6 0.55]; %gamma_r values from Table VI-5-3 for calculation of runup
h = [20.5 20.2 20.1 20; 24.6 24.3 24.2 24.1; 27 26.6 26.5 26.4]; %water surface elevation against levee from HEC-RAS models

%Load fetch data from GIS analysis
dir = 'e:\CBEC\Projects\19-1034_YBEL\Wind_Wave_Analysis\';
fid = fopen([dir 'AMR_2_5_10yr_MaxDepths_v3.csv']);
fetch = textscan(fid,'%f%f%s%s%f%f%f%f%f%f','Delimiter',',','headerlines',1);
f = [fetch{5:7}]; %just the fetch lengths
depth = [fetch{8:10}]; %avg water depths along transects
station = fetch{3};stationid = unique(station);
wsname = {'twoyr';'fiveyr';'tenyr'}; %data structure field names
fn = {'id';'StationID';'WindDirection';'FetchLength';'WindSpeed';...
    'WaterDepth';'WaveHeight';'WavePeriod';'WindSetUp';'WaveRunUp';...
    'TWL';'Bench';'Vellinga';'Tfetch';'Tkmax';'Dn50'}; %field names for the data table
Tfetch = zeros(length(f),1); %Wind duration
Ws = zeros(length(f),1); %wind setup elevation
Hs = zeros(length(f),1); %wave height
Tp = zeros(length(f),1); %wave period
R2pct = zeros(length(f),1); %2 percent wave runup elevation
Dn50 = zeros(length(f),1); %required median nominal riprap size given wave conditions
bench = zeros(length(f),1); %bench length based on FEMA report
Vellinga = zeros(length(f),1); %Vellinga intrusion length
Tkmax = zeros(length(f),1); %Duration to failure for erosion
Umph = zeros(length(f),1); %Wind velocity in mph
for i = 1:length(wsname)
    fprintf('Running analysis for %s water level scenario\n',wsname{i})
    for j = 1:length(stationid)
        idx = strcmp(station,stationid(j));
        F = (f(idx,i))'; %fetch in ft
        d =  (depth(idx,i))'; %depth in ft

        %compute fetch-limited duration (eq. 2.1 from NHC 2011)
        Fm = F./3.281; %fetch in m
        gm = 9.81; %grav constant in SI units
        Um = U./2.237; %wind speed in m/s
        
        tfetch = (77.23*((Fm.^0.67)./((gm^0.34)*(Um.^0.33))))./3600; %in hours
        
        Ua = (0.589*(U.^1.23))*1.46667; %extra term converts mph to ft/s for H and T calcs

        %Wind set-up
        S = ((U.^2).*(F./5280))./(1400.*d); %units in imperial, Eq. C-2 EM 1110-2-1420
        
        %Wave height and period: eqns. from Shore Protection Manual (1984)
        g = 32.2; %grav constant in ft/s^2        
        
%         gHUasqd = 0.283*tanh(0.53*((g*d)./(Ua.^2)).^(3/4)).*tanh((0.00565*((g*F)./(Ua.^2)).^(1/2))./tanh(0.530*((g*d)./(Ua.^2)).^(3/4)));
%         gTUa = 7.54*tanh(0.833*((g*d)./(Ua.^2)).^(3/8)).*tanh((0.0379*((g*F)./(Ua.^2)).^(1/3))./tanh(0.833*((g*d)./(Ua.^2)).^(3/8)));
%         H = gHUasqd.*((Ua.^2)/g); %wave height in ft
%         T = gTUa.*(Ua/g); %wave period in s

        
        %eqns. from Hurdle & Stive (1989)--not using these        
        firstterm = 0.25*((Ua.^2)/g).*tanh(0.6*((g*d)./(Ua.^2)).^(3/4));
        secondterm = (tanh((4.3E-5*((g*F)./(Ua.^2)))./(tanh(0.6*((g*d)./(Ua.^2)).^(3/4))).^2)).^(1/2);  
        H = firstterm.*secondterm;
        
        firstterm = 8.3*((Ua)/g).*tanh(0.76*((g*d)./(Ua.^2)).^(3/8));
        secondterm = (tanh((4.1E-5*((g*F)./(Ua.^2)))./(tanh(0.76*((g*d)./(Ua.^2)).^(3/8))).^3)).^(1/3);
        T = firstterm.*secondterm;
        
        %Wave properties
        lambda = (g.*(T.^2))./(2*pi); %USACE approximation for wavelength
        %estimate deep vs. shallow water behavior
        if any(H./lambda > 0.5)
            disp('Waves are DEEP!')
        elseif any(H./lambda < 0.04)
            disp('Waves are SHALLOW!')
        else
            disp('Waves are TRANSITIONAL!')
        end
        
        %calculate Iribarren number (surf similarity parameter)
        a = 18.4; %slope of levee in degrees        
        S0 = (H./lambda);
        Xi = tand(a)./(sqrt(S0)); %surf similarity parameter V1-5-1
        
        if any(Xi) <= 0.5
            error('Xi < 0.5!')
        end
        A = 1.77; %CEM Table V1-5-2 for wind/wave conditions during this study
        gamr = gammar(i); %CEM Table V1-5-3
        gamb = 1; %Assumed no bermed fetch
        gamh = 1; %Assumed no shallow wave influence
        gambeta = zeros(1,3);
        ifW = strcmp(fetch{4}(idx),'W');
        gambeta(ifW) = 1; %Waves coming from along-transect
        gambeta(~ifW) = 1-0.0022*45; %waves approaching from a 45 deg angle (e.g., SW and NW winds)
        r2pct = 1.77.*H*gamr*gamb*gamh.*gambeta.*Xi;
        
        %estimate the riprap size needed to withstand waves (EM-1110-2-1100, eq. V1-5-85)
        gams = 165; %lbf/ft^3 (stone)
        gamw = 62.4; %lbf/ft^3 (fresh water)
        Kd = 16; %V1-5-81
        delta = (gams/gamw)-1;
        W50 = (gams*(H.^3))./(Kd*(delta^3)*cotd(a)); %weight of median stone (lbs)
        dn50 = (W50/gams).^(1/3); %median nominal size of riprap required to withstand waves (ft)

        %Compute erosion depth using the two methods outlined in the NHC
        %(2011) GRR report
        mZ = cotd(a)*r2pct; %eroded bench length (MK&A method -- FEMA, 2004)
        w = 0.05; %fall velocity of fine sand in m/s that NHC used
        Le = ((0.39*w^0.44)^-1.28)*(H./3.281).^1.28; %in meters
        
        %compute time-to-erosion (eq. 4.3 in NHC study)
        Tpeak = 1.1*T;
        HR = 0.673*((Tpeak.^0.5)./(H.^0.25)).*H;
        ed = 2/12; %2 inches converted to ft
        sf = 3; %factor of safety that NHC used
        Ce = 9.14E-7; %erosion coef of grasses that NHC used
        tkmax = ed./(3600*sf*Ce*(4*HR*tand(a)).^2);
                
        %Save out data
        Tfetch(idx) = tfetch;
        Ws(idx) = S; 
        Hs(idx) = H;
        Tp(idx) = T;
        R2pct(idx) = r2pct;
        Dn50(idx) = dn50.*12; %convert ft to in
        bench(idx) = mZ;
        Vellinga(idx) = Le.*3.281; %convert m to ft
        Tkmax(idx) = tkmax;
        Umph(idx) = U;
    end
    wse = [repmat(h(i,1),3,1); repmat(h(i,2),3,1); repmat(h(i,3),3,1); repmat(h(i,4),3,1)]; %station depth (ft)
    int = {fetch{2} fetch{3} fetch{4} f(:,i) Umph depth(:,i) Hs Tp Ws R2pct wse+Ws+R2pct bench Vellinga Tfetch Tkmax Dn50};
    data.(wsname{i}) = cell2struct(int,fn,2);
end
% display data as tables for easy copying to excel
disp('Two-year RI data')
writetable(struct2table(data.twoyr),'WindWave_twoyear.csv','Delimiter',' ')
type 'WindWave_twoyear.csv'
disp('Five-year RI data')
writetable(struct2table(data.fiveyr),'WindWave_fiveyear.csv','Delimiter',' ')
type 'WindWave_fiveyear.csv'
disp('Ten-year RI data')
writetable(struct2table(data.tenyr),'WindWave_tenyear.csv','Delimiter',' ')
type 'WindWave_tenyear.csv'
