%Analyze Water Level Gauge data for the Pacheco Project to determine the
%time lag between WSE time series.
%
% 03/12/2020 - B.K. Norris - cbec eco-engineering
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%% Load in the instrument data
%These data from 2019_07_31:
filesdir = 'e:\CBEC\Projects\19_1011_Pacheco\Data\Instruments\2019_07_31\';
fname{1} = 'InstrumentData_2019_07_31_C1.csv';
fname{2} = 'InstrumentData_2019_07_31_C4.csv';
fname{3} = 'InstrumentData_2019_07_31_USGS.csv';
instnames = {'C1';'C4';'USGS'}; %these names need to be in the same order as the files in fname{i}
%Read the data into MATLAB
FRMT = {'%s%f';'%s%f';'%s%s%s%f'}; %these format strings need to be in the same order as the files in fname{i}
for i = 1:3
    fid = fopen([filesdir fname{i}]);
    header = fgetl(fid);header = regexp(header,',','split');
    datfile = textscan(fid,FRMT{i},'delimiter',',');
    if i < 3 %The USGS gauge has datetime as the third column, C1 and C4 has datetime as the first column
        whichcol = 1;
    else
        whichcol = 3;
    end
    datenumb = datenum(datfile{whichcol}); %convert date string to matlab date numbers for easy plotting
    data.(instnames{i}) = cell2struct({datenumb datfile{whichcol+1}},{header{whichcol:whichcol+1}},2);
end
%% Plot the data
f1 = figure(1);
p = zeros(3,1);
cc = hsv(3);
hold on
for i = 1:3
    %normalize data so everything plots easily
    data_norm = (data.(instnames{i}).Stage-min(data.(instnames{i}).Stage))./...
        (max(data.(instnames{i}).Stage)-min(data.(instnames{i}).Stage));
    p(i) = plot(data.(instnames{i}).Datetime,data_norm,...
        'linewidth',1.5,...
        'color',cc(i,:));
end
hold off;
leg = legend(p,instnames,'location','northeast');
datetickzoom('x','mm-dd-yy')
ylabel('Normalized WSE [ft]')
xlabel('Date Time')

%% Time series correlation
start = datenum('30-Apr-2019 16:30:00');
stop = datenum(0,0,15,0,0,0);
id1 = data.C1.Datetime>=start & data.C1.Datetime <= start+stop;
x1 = data.C1.Stage(id1);
id2 = data.USGS.Datetime>=start & data.USGS.Datetime <= start+stop;
x2 = data.USGS.Stage(id2);
x1norm = smooth((x1-min(x1))/(max(x1)-min(x1)),5); %normalize data to remove differences in magnitude (stage)
x2norm = smooth((x2-min(x2))/(max(x2)-min(x2)),5);
dt = 15; %min, timestep of WSE data

%plot the data and cross-correlation
f2 = figure(2);
p = zeros(3,1);
sp(1) = subplot(211);
p(1) = plot(data.C1.Datetime(id1),x1norm,'r');hold on
p(2) = plot(data.USGS.Datetime(id2),x2norm,'b');
datetickzoom('x','mm-dd-yy')
ylabel('Normalized WSE [-]')
xlabel('Date Time')
sp(2) = subplot(212);
[cor,lags]=xcorr(x1norm-mean(x1norm),x2norm-mean(x1norm));
p(3) = plot((lags*dt)./60./24,cor,'-k');
xlabel('Time Lag [Days]')
ylabel('Correlation')
[~,I] = max(abs(cor));
timeDiff = (lags(I)*dt)/60/24; %time difference between C1 and USGS in days
fprintf('Time difference between C1 and USGS Gauge: %0.1f days\n',timeDiff)
prettyfigures('box',1)

