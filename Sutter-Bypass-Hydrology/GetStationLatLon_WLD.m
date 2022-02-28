clear
path = 'e:\CBEC\Projects\18_1012_Burlingame\HH\Tides\';
tofind = 'A05926';
for i = 1
    options = weboptions('Timeout',15);
    api1 = 'http://wdl.water.ca.gov/WaterDataLibrary/WaterQualityDataLib.aspx?StationNumber=A0510300&MapStationName=FEATHER';
    table1 = webread(api1,options);
    stationid = strfind(table1,'WaterQualityStationId=');
    id = table1(stationid+22:stationid+25); %should be a 4-digit number
    fprintf('Station ID: %s\n',id)
    api2 = ['http://wdl.water.ca.gov/WaterDataLibrary/StationDetailsPopup.aspx?WaterQualityStationId=' id];
    table2 = webread(api2,options);
    
    %% Extract just the data you need
    
    %Latitude
    latid = strfind(table2,'Latitude:</div ></td>');
    latdms = table2(latid+46:latid+57);
    latdms = regexp(latdms,'-','split');
    latdd = str2double(latdms{1})+(str2double(latdms{2})/60)+(str2double(latdms{3}))/3600;
    
    %Longitude
    lonid = strfind(table2,'Longitude:</div ></td>');
    londms = table2(lonid+47:lonid+60);
    londms = regexp(londms,'-','split');
    londd = -1*(str2double(londms{2})+(str2double(londms{3})/60)+(str2double(londms{4}))/3600);
    
    %Compile table of info to copy into excel
    data = table(api1,latdd,londd);
    clipboard('copy',data)
    
    
end