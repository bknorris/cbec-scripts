clear
path = 'e:\CBEC\Projects\20-1006_Sutter_Bypass\Data\';
fid = fopen([path 'FlowGauges.csv']);
data = textscan(fid,'%s%s','delimiter',',');
StationNumber = data{1};
apis = cell(length(StationNumber),1);
lats = zeros(length(StationNumber),1);
lons = zeros(length(StationNumber),1);
for i = 1:length(StationNumber)
    options = weboptions('Timeout',15);
    api1 = ['https://cdec.water.ca.gov/dynamicapp/staMeta?station_id=' StationNumber{i}];
    table1 = webread(api1,options);
    stationid = strfind(table1,'Station ID');
    if isempty(stationid)
        continue
    else
    id = table1(stationid+23:stationid+25); %should be 3 characters
    id = regexp(id,'[^a-zA-Z]','split');
    id = id(~cellfun('isempty',id));
    fprintf('Station ID: %s\n',id{1})
    fprintf('Station Number: %s\n',StationNumber{i})    
    %% Extract just the data you need
    
    %Latitude
    latid = strfind(table1,'Latitude</b></td><td>');
    latdms = table1(latid+21:latid+29);
    latdms = regexp(latdms,'[^0-9.]','split');
    latdms = latdms(~cellfun('isempty',latdms));
    latdd = str2double(latdms{1});
%     fprintf('Latitude: %0.4f\n',latdd)
    
    %Longitude
    lonid = strfind(table1,'Longitude</b></td><td>');
    londms = table1(lonid+22:lonid+32);
    londms = regexp(londms,'[^0-9.]','split');
    londms = londms(~cellfun('isempty',londms));
    londd = -1*(str2double(londms{1}));
%     fprintf('Longitude: %0.4f\n',londd)
    
    %Compile data
    apis{i} = api1;
    lats(i) = latdd;
    lons(i) = londd;
    fprintf('\n')
    end
end
datatable = cell2table({apis,lats,lons},'variablenames',{'api','lat','lon'});

