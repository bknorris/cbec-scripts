%Code to automate the delination of berms of rasters using the topo toolbox
%https://topotoolbox.wordpress.com/topotoolbox/
%
%
clear

%% Automated Processes: Load DEM, Create Flow Dir Layer
dirc = 'e:\CBEC\Projects\20-1006_Sutter_Bypass\Data\GIS\Rasters\';
rasname = 'DWR_2019_FirstArea2.tif';fname = regexprep(rasname,'.tif','');
t1 = tic;
if ~any(exist([dirc fname '_flowD.mat'],'file'))
DEM = GRIDobj([dirc rasname]);
DEM.Z = DEM.Z.*-1; %first invert the DEM to find ridges
fprintf('DEM loaded in %0.2f minutes\n',toc(t1)/60)

t2 = tic;
DEMf = fillsinks(DEM);

    FD = FLOWobj(DEMf,'preprocess','carve'); %THIS STEP TAKES AWHILE TO RUN!
    save([dirc fname '_flowD.mat'],'FD')
    fprintf('Flow direction analysis completed in in %0.2f minutes\n',toc(t2)/60)
else
    load([dirc fname '_flowD.mat'])
end
A  = flowacc(FD);

%% Edit these values to generate different auto berm IDs
figdir = 'e:\CBEC\Projects\20-1006_Sutter_Bypass\Data\Figures\';
shapedir = 'e:\CBEC\Projects\20-1006_Sutter_Bypass\Data\GIS\Shapes\FlowDirectionAnalysis\';
ww = [25000 50000 100000 150000 200000 250000 300000]; %try a couple of options
f1 = zeros(size(ww));
p2 = zeros(size(ww));
legtext = {size(ww)};
cc = brewermap(length(ww),'Spectral');
close all
for i = 1:length(ww)
    W = ww(i); %set cell size for flow directions
    S = STREAMobj(FD,'minarea',W);
    S2  = klargestconncomps(S);
    S2 = removeshortstreams(S2,W/25); %reduce W by some amount (values of 15 to 30 work well)
%     y = smooth(S2,S2.y,'k',15000,'nstribs',true);
%     x = smooth(S2,S2.x,'k',15000,'nstribs',true);
%     [~,~,xs,ys] = STREAMobj2XY(S2,x,y);
%     S2.y = ys;S2.x = xs;
    
    f1(i) = figure(i);
    p1 = plot(S,'k');hold on
    p2 = plot(S2,'r','linewidth',1.5);
    leg = legend([p1 p2],{'Original';'Reduced'},'location','southeast');
    ylabel('Northing (ft)')
    xlabel('Easting (ft)')
    title(['\bfCell size W = ' sprintf('%0.0f',W)])
    hold off
    
    f2 = figure(length(ww)+1);
    legtext{i} = ['W = ' sprintf('%0.0f',ww(i))];
    p2(i) = plot(S2,'color',cc(i,:),'linewidth',1.5);hold on
    if i == 1
        ylabel('Northing (ft)')
        xlabel('Easting (ft)')
    end
    %     f2 = figure(4);
    %     sp(1) = subplot(121);
    %     plot(S2);
    %     sp(2) = subplot(122);
    %     [Z,x,y] = GRIDobj2mat(DEM);
    %     imagesc(x,flipud(y),Z)
    %     cb = colorbar;
    %     set(cb,'location','eastoutside')
    %     ylabel(cb,'Elevation [ft]')
    %     set(sp(1),'position',[0.15 0.15 0.35 0.75],...
    %         'tickdir','out')
    %     set(sp(2),'position',[0.55 0.15 0.28 0.75],...
    %         'yticklabel',[],...
    %         'tickdir','out')
    %     ylabel(sp(1),'Northing (ft)')
    %     xlabel(sp(1),'Easting (ft)'),xlabel(sp(2),'Easting (ft)')
    %     suplabel(['\bfCell size W = ' sprintf('%0.0f',W)],'t');
    %     export_fig(f2,[figdir 'AutoBermDel_W' sprintf('%0.0f',W) '.png'],'-png','-nocrop')
    % shapedir = 'e:\CBEC\Projects\20-1006_Sutter_Bypass\Data\GIS\Shapes\FlowDirectionAnalysis\';
    MS = STREAMobj2mapstruct(S2);
    save([shapedir fname 'W' sprintf('%0.0f',ww(i)) '_FD.mat'],'MS')
end
leg = legend(legtext,'location','southeast');
prettyfigures('text',11,'labels',11,'box',1)
for i = 1:length(ww)
    export_fig(f1(i),[figdir 'AutoBermReduc_' legtext{i} '.png'],'-png','-nocrop')
end
export_fig(f2,[figdir 'AllBermEsts.png'],'-png','-nocrop')

%% Export Flow Direction data to .xyz for GIS


