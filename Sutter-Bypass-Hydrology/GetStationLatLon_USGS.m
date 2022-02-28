clear
path = 'e:\CBEC\Projects\20-1006_Sutter_Bypass\Data\';
fid = fopen([path 'USGSGauges.csv']);
data = textscan(fid,'%s%s','delimiter',',');
StationNumber = data{1};
StationName = data{2};
apis = cell(length(StationNumber),1);
lats = zeros(length(StationNumber),1);
lons = zeros(length(StationNumber),1);
for i = 1:length(StationNumber)
    options = weboptions('Timeout',15);
    api1 = ['https://waterdata.usgs.gov/nwis/inventory/?site_no=' StationNumber{i} '&agency_cd=USGS'];
    table1 = webread(api1,options);
    stationid = strfind(table1,['USGS ' StationNumber{i} ' ' StationName{i}]);
    if isempty(stationid)
        continue
    else
    fprintf('Station Number: %s\n',StationNumber{i})
    fprintf('Station Name: %s\n',StationName{i}) 
    %% Extract just the data you need
    
    %Latitude
    latid = strfind(table1,'<dd>Latitude');
    latdms = table1(latid+14:latid+27);
    latdms = regexp(latdms,'[^\dNSEW.]+','split'); %second cell contains the 'degree' symbol (176)!
    latdms = latdms(~cellfun('isempty',latdms));
    latdd = str2double(latdms{1})+(str2double(latdms{3})/60)+(str2double(latdms{4}))/3600;
%     fprintf('Latitude: %0.4f\n',latdd)
    
    %Longitude
    lonid = strfind(table1,'&nbsp; Longitude');
    londms = table1(lonid+17:lonid+30);
    londms = regexp(londms,'[^\dNSEW.]+','split'); %second cell contains the 'degree' symbol (176)!
    londms = londms(~cellfun('isempty',londms));
    londd = -1*(str2double(londms{1})+(str2double(londms{3})/60)+(str2double(londms{4}))/3600);
%     fprintf('Longitude: %0.4f\n',londd)
    
    %Compile data
    apis{i} = api1;
    lats(i) = latdd;
    lons(i) = londd;
    fprintf('\n')
    end
end
datatable = cell2table({apis,lats,lons},'variablenames',{'api','lat','lon'});

