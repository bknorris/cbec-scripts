clear
path = 'e:\CBEC\Projects\20-1006_Sutter_Bypass\Data\';
fid = fopen([path 'WQGauges.csv']);
data = textscan(fid,'%s%s%s','delimiter',',');
StationNumber = data{1};
MapStationName = data{3};
apis = cell(length(StationNumber),1);
lats = zeros(length(StationNumber),1);
lons = zeros(length(StationNumber),1);
for i = 1:length(StationNumber)
    options = weboptions('Timeout',15);
    api1 = ['http://wdl.water.ca.gov/WaterDataLibrary/WaterQualityDataLib.aspx?StationNumber=' StationNumber{i} '&MapStationName=' MapStationName{i}];
    table1 = webread(api1,options);
    stationid = strfind(table1,'WaterQualityStationId=');
    if isempty(stationid)
        continue
    else
    id = table1(stationid+20:stationid+28); %should be a 4-digit number
    id = regexp(id,'[^0-9.]','split');
    id = id(~cellfun('isempty',id));
    fprintf('Station ID: %s\n',id{1})
    fprintf('Station Number: %s\n',StationNumber{i})
    fprintf('Station Name: %s\n',MapStationName{i})
    api2 = ['http://wdl.water.ca.gov/WaterDataLibrary/StationDetailsPopup.aspx?WaterQualityStationId=' id{1}];
    table2 = webread(api2,options);
    
    %% Extract just the data you need
    
    %Latitude
    latid = strfind(table2,'Latitude:</div ></td>');
    latdms = table2(latid+40:latid+60);
    latdms = regexp(latdms,'[^0-9.]','split');
    latdms = latdms(~cellfun('isempty',latdms));
    latdd = str2double(latdms{1})+(str2double(latdms{2})/60)+(str2double(latdms{3}))/3600;
%     fprintf('Latitude: %0.4f\n',latdd)
    
    %Longitude
    lonid = strfind(table2,'Longitude:</div ></td>');
    londms = table2(lonid+40:lonid+70);
    londms = regexp(londms,'[^0-9.]','split');
    londms = londms(~cellfun('isempty',londms));
    londd = -1*(str2double(londms{1})+(str2double(londms{2})/60)+(str2double(londms{3}))/3600);
%     fprintf('Longitude: %0.4f\n',londd)
    
    %Compile data
    apis{i} = api1;
    lats(i) = latdd;
    lons(i) = londd;
    fprintf('\n')
    end
end
datatable = cell2table({apis,lats,lons},'variablenames',{'api','lat','lon'});

