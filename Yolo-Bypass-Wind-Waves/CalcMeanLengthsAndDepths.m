clear
dir = 'e:\CBEC\Projects\19-1034_YBEL\Wind_Wave_Analysis\';
fid = fopen([dir 'AMR_2yrMaxDepth_MATLAB.csv']);
data = textscan(fid,'%f%f%f%s%s%f%f%f','Delimiter',',','headerlines',1);
stations = unique(data{4});
winddir = unique(data{5});
ids = unique(data{3});
lengths = zeros(size(ids));
depths = zeros(size(ids));
stations2 = cell(size(ids));
winddir2 = cell(size(ids));
for i = 1:length(ids) %num of individual fetch lines
    idx = find(data{3}==ids(i));
    ll = data{2}(idx);
    dd = data{8}(idx);
    lengths(i) = ll(find((~isnan(dd) & dd~=0),1,'last'));
    depths(i) = mean(dd((~isnan(dd) & dd~=0)));
    stations2(i) = unique(data{4}(idx));
    winddir2(i) = unique(data{5}(idx));
end
fetchlength = zeros(4,3);
avgdepth = zeros(4,3);
stations3 = cell(4,3);
winddir3 = cell(4,3);
angles = [-9 -6 -3 0 3 6 9];
for i = 1:length(stations)
    idx1 = find(strcmp(stations2,stations{i}));
    for j = 1:length(winddir)
         idx2 = find(strcmp(winddir2,winddir{j}));
         idx3 = intersect(idx1,idx2);
         stations3(i,j) = unique(stations2(idx3));
         winddir3(i,j) = unique(winddir2(idx3));
         fetchlength(i,j) = mean((sum(lengths(idx3).*cosd(angles))./sum(cosd(angles))));
         avgdepth(i,j) = mean(depths(idx3));
    end
end
reshape(stations3,12,1) 
reshape(winddir3,12,1) 
reshape(fetchlength,12,1) 
reshape(avgdepth,12,1)

