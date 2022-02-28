%Read the SC sonde data, do some basic QC and plot results.
clear
close all
dir = 'e:\CBEC\Projects\15-1017-3_Smith_Canal\Data\';
figpath = 'e:\CBEC\Projects\15-1017-3_Smith_Canal\Figures\';
filename = 'Sonde1.csv';
frmt = '%s%s%s%f%f%f%f%f%f%f';
fid = fopen([dir filename]);
header = regexp(fgetl(fid), ',', 'split'); %column names
data = textscan(fid,frmt,'delimiter',',');
sonde = cell2struct(data,header,2); %matlab structure containing fields with csv column headers and data
%% Plot some results
%First, convert datetime from csv into matlab numerical datetime
x = cellfun(@datenum, sonde.DateTime); %if this line doesn't work, check the csv for gaps!
fields = {'Temperature';'Specific Conductivity (mS/cm)';'Salinity (ppt)';'Depth (ft)';...
    'ODO Saturation (%)';'ODO (mg/L)';'Battery (volts)'};
c = 1; %plot counter
for i = 4:10 %skip the first three columns because these contain dates etc
    y = sonde.(header{i});
    y = fillmissing(y,'pchip');
    if strcmp(header{i},'Depth')
        y = y-34; %adjust for atmospheric pressure
        ff = figure(c);
        p(1) = plot(x,y,'b');
    else
        %run the despiking routine
        yfilt = movmedian(y,40); %calculate a running median of ts
        yqc = zeros(length(y),1);
        for j = 1:length(y)
            if y(j)/yfilt(j) < 1.08 && y(j)/yfilt(j) >= 1
                yqc(j) = y(j);
            elseif y(j)/yfilt(j) > 0.92 && y(j)/yfilt(j) <= 1
                yqc(j) = y(j);
            else
                yqc(j) = NaN;
            end
        end
        %fill missing gaps minus leading nans if any
        fnanid = find(~isnan(yqc),1,'first');
        yqc(fnanid:end) = fillmissing(yqc(fnanid:end),'pchip');
        ff = figure(c);
        p(1) = plot(x,y,'b');hold on
        p(2) = plot(x,yfilt,'k');
        p(3) = plot(x,yqc,'r');
        legend(p,{'Raw Data';'45-pt median filter';'Filtered Data'},...
            'location','southeast')
    end
    set(gca,'xlim',[min(x) max(x)],'ylim',[min(y) max(y)])
    datetick('x','mm-dd-yy HH:MM')
    ylabel(fields{c})
    set(ff,'units','inches')
    pos = get(ff,'Position');
    set(ff,'PaperPositionMode','Auto','PaperUnits','Inches',...
        'PaperSize',[pos(3) pos(4)])
    fname = regexprep(filename,'.csv','');
    print(ff,[figpath fname '_' header{i}],'-dpng','-r0')
    c=c+1;
end
