% Produces analyses for the Wrinkles project
% input requires a cell array of folders (even if only one folder)


% The macro below  suppresses warnings in the editor about arrays growing 
% in a loop (which is poor coding, but helps readability)
%#ok<*AGROW>

function wrAnalyse (folderlist)

nfolders = length(folderlist);

% Array of traces and means for each participant, plus N for each condition
nDr = 0; nWe = 0; nWr = 0; 
DRYLF = []; WETLF = []; WRILF = []; 
DRYGF = []; WETGF = []; WRIGF = []; 
DRYlag= []; WETlag= []; WRIlag= [];
DRYerr= []; WETerr= []; WRIerr= [];
DRYatt= []; WETatt= []; WRIatt= [];
DRYdec= []; WETdec= []; WRIdec= [];
% Demographics
AgeDr = []; SexDr = []; HandDr = [];
AgeWe = []; SexWe = []; HandWe = [];
AgeWr = []; SexWr = []; HandWr = [];

%% COLLATE DATA
%
% Run through folders, adding data from files to arrays

for f = 1:nfolders
    
    folder = folderlist{f};
    disp(folder)
    dirin = cd;
    cd (folder);
    d = dir('*.mat');
    n = length(d);
    
    for i=1:n
        disp(d(i).name)
        load (d(i).name);
        switch condName
            case 'Dry'
                [G,L,R,lag,n,Att,Dec] = int_getdata(DATA);
                if n>5
                    nDr = nDr+1;
                    DRYLF = [DRYLF; L]; %#ok<*AGROW>
                    DRYGF = [DRYGF; G]; 
                    DRYlag= [DRYlag; lag];
                    DRYerr= [DRYerr; R];
                    DRYatt = [DRYatt; Att];
                    DRYdec = [DRYdec; Dec];
                    AgeDr = [AgeDr; Age];
                    if strcmp(Sex,'F')
                        SexDr = [SexDr; 1];
                    else
                        SexDr = [SexDr; 0];
                    end
                    if strcmp(Hand,'L')
                        HandDr = [HandDr; 1];
                    elseif strcmp(Hand,'R')
                        HandDr = [HandDr; 0];
                    else
                        fprintf('Bad hand %s\n',d(i).name)
                    end
                else
                    % Too few files - skip but print name
                    fprintf('Exc: %s\n', d(i).name)
                end
            case 'Wet'
                [G,L,R,lag,n,Att,Dec] = int_getdata(DATA);
                if n>5
                    nWe = nWe+1;
                    WETLF = [WETLF; L];
                    WETGF = [WETGF; G];
                    WETlag= [WETlag; lag];
                    WETerr= [WETerr; R];
                    WETatt = [WETatt; Att];
                    WETdec = [WETdec; Dec];
                    AgeWe = [AgeWe; Age];
                    if strcmp(Sex,'F')
                        SexWe = [SexWe; 1];
                    else
                        SexWe = [SexWe; 0];
                    end
                    if strcmp(Hand,'L')
                        HandWe = [HandWe; 1];
                    elseif strcmp(Hand,'R')
                        HandWe = [HandWe; 0];
                    else
                        fprintf('Bad hand %s\n',d(i).name)
                    end
                else
                    % Too few files - skip but print name
                    fprintf('Exc: %s\n', d(i).name)
                end
            case 'Wrinkly'
                [G,L,R,lag,n,Att,Dec] = int_getdata(DATA);
                if n>5
                    nWr = nWr+1;
                    WRILF = [WRILF; L];
                    WRIGF = [WRIGF; G];
                    WRIlag= [WRIlag; lag];
                    WRIerr= [WRIerr; R];
                    WRIatt = [WRIatt; Att];
                    WRIdec = [WRIdec; Dec];
                   AgeWr = [AgeWr; Age];
                    if strcmp(Sex,'F')
                        SexWr = [SexWr; 1];
                    else
                        SexWr = [SexWr; 0];
                    end
                    if strcmp(Hand,'L')
                        HandWr = [HandWr; 1];
                    elseif strcmp(Hand,'R')
                        HandWr = [HandWr; 0];
                    else
                        fprintf('Bad hand %s\n',d(i).name)
                    end
                else
                    % Too few files - skip but print name
                    fprintf('Exc: %s\n', d(i).name)
                end
        end
    end
    
    % End of processing, return to original directory
    cd (dirin)
    
end


%% DISPLAY DATA
%
% Show figures and descriptive statistics

% Number of participants
fprintf('\n\n')
fprintf('Numbers:\n  Wet: %d\n  Dry: %d\n  Wri: %d\n\n',nWe,nDr,nWr)
% Sex of participants
fprintf('Sex ratio:\n  Wet: %d F, %d M\n  Dry: %d F, %d M\n  Wri: %d F, %d M\n\n',...
    length(find(SexWe==1)),length(find(SexWe==0)),...
    length(find(SexDr==1)),length(find(SexDr==0)),...
    length(find(SexWr==1)),length(find(SexWr==0)))
% Handedness
fprintf('Handedness:\n  Wet: %d L, %d R\n  Dry: %d L, %d R\n  Wri: %d L, %d R\n\n',...
    length(find(HandWe==1)),length(find(HandWe==0)),...
    length(find(HandDr==1)),length(find(HandDr==0)),...
    length(find(HandWr==1)),length(find(HandWr==0)))
% Ages (NB age could be NaN, so skip these trials
fprintf('Mean Age:\n  Wet: %2.2f\n  Dry: %2.2f\n  Wri: %2.2f\n(Min %2d, Max %2d, Mean %2.2f, SD %2.3f)\n\n',...
    nanmean(AgeWe),nanmean(AgeDr),nanmean(AgeWr),...
    min([AgeWe;AgeDr;AgeWr]),max([AgeWe;AgeDr;AgeWr]),nanmean([AgeWe;AgeDr;AgeWr]),nanstd([AgeWe;AgeDr;AgeWr]))

% Get mean LF and GF for each participant, then snip out a segment of
% interest.
DRYLFm = mean(DRYLF,1);
WETLFm = mean(WETLF,1);
WRILFm = mean(WRILF,1);
DRYGFm = mean(DRYGF,1);
WETGFm = mean(WETGF,1);
WRIGFm = mean(WRIGF,1);
seg = 7001:10000;
DRYLFseg = DRYLF(:,seg);
WETLFseg = WETLF(:,seg);
WRILFseg = WRILF(:,seg);
DRYGFseg = DRYGF(:,seg);
WETGFseg = WETGF(:,seg);
WRIGFseg = WRIGF(:,seg);
% This bit is ugly, but is used to build the error margin (standard error)
% around the traces in the figure
xc = [seg seg(end:-1:1)];
wem = mean(WETGFseg,1);
wes = std(WETGFseg,[],1)./sqrt(size(WETGFseg,1));
wec= [wem+wes wem(end:-1:1)-wes];
wrm = mean(WRIGFseg,1);
wrs = std(WRIGFseg,[],1)./sqrt(size(WRIGFseg,1));
wrc= [wrm+wrs wrm(end:-1:1)-wrs];
drm = mean(DRYGFseg,1);
drs = std(DRYGFseg,[],1)./sqrt(size(DRYGFseg,1));
drc= [drm+drs drm(end:-1:1)-drs];


% plot GF and LF traces, with error margins
figure
plot(int_getTargetTrace(),'k','LineWidth',2)
hold on;
plot ([DRYLFm' WETLFm' WRILFm'],'LineWidth',2);
legend('Target','Dry','Wet','Wrinkled','Location','NorthWest')
fill(xc,drc,'b','FaceAlpha',0.4)
fill(xc,wec,'g','FaceAlpha',0.4)
fill(xc,wrc,'r','FaceAlpha',0.4)
plot ([DRYLFm' WETLFm' WRILFm'],'LineWidth',2);
plot ([DRYGFm' WETGFm' WRIGFm']);
set(gca,'XTick',1000:2000:15000)
set(gca,'YTickLabel',{'0.0','0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0'});
set(gca,'YLim',[0 4.2])
xlabel('Time (ms)')
ylabel('Force (N)')
set(gcf,'Position',[186 49 1122 636])

figure; hold on
line([.5 2.0],[.5 2.0],'color','k',...
    'LineStyle','--','LineWidth',2)
plot (mean(DRYLF(:,1000:end),1),mean(DRYGF(:,1000:end),1),'b','LineWidth',2);
plot (mean(WETLF(:,1000:end),1),mean(WETGF(:,1000:end),1),'g','LineWidth',2);
plot (mean(WRILF(:,1000:end),1),mean(WRIGF(:,1000:end),1),'r','LineWidth',2);
legend('Target','Dry','Wet','Wrinkled','Location','SouthEast')
xlabel('Load force (N)')
ylabel('Grip force (N)')
set (gca, 'YLim', [0 4.25])
set (gca, 'XLim', [0.2 2.25])
set (gca, 'YTick',0:.5:4)
set (gca, 'XTick',0.5:.25:2)


%% STATISTICS
%
% Inferential statistics 
% - assumes Stats Toolbox is available, will generate errors if not

% Grip:Load ratio - one-way analysis of variance
DryPC = mean(DRYGFseg,2)./mean(DRYLFseg,2);
WetPC = mean(WETGFseg,2)./mean(WETLFseg,2);
WriPC = mean(WRIGFseg,2)./mean(WRILFseg,2);
Y = [DryPC; WetPC; WriPC];
G = [ones(size(DryPC)); 2*ones(size(WetPC)); 3*ones(size(WriPC))];
[p,tbl,stats] = anova1(Y,G,'off');
if p>=0.05
    fprintf('No significant effect of condition on GF:LF ratio: F(%d,%d)=%2.3f, p=%1.4f\n\n',...
        tbl{2,3},tbl{3,3},tbl{2,5},p);
else
    fprintf('Significant effect of condition on GF:LF ratio: F(%d,%d)=%2.3f, p=%1.4f\n\n',...
        tbl{2,3},tbl{3,3},tbl{2,5},p);
    stats.gnames = {'Dry';'Wet';'Wrinkly'};
    figure
    c = multcompare(stats)
end


% Correlation of age and GF-LF ratio
A = [AgeDr;AgeWe;AgeWr];
I = find(~isnan(A));
[r,p] = corr(Y(I),A(I));
if p>=0.05
    fprintf('GF:LF ratio does not correlate with age: r(%d)=%1.3f, p=%1.4f\n\n',...
        length(Y)-2,r,p);
else
    fprintf('GF:LF ratio correlates significantly with age: r(%d)=%1.3f, p=%1.4f\n',...
        length(Y)-2,r,p);
    [B,~,~,~,stats] = regress (Y(I),[A(I) ones(size(A(I)))]);
    fprintf('Ratio = %3.3f * Age + %3.3f. R2 = %2.3f\n\n',B(1),B(2),stats(1))
end



% Grip-load lag - one-way analysis of variance, plus one-sample t-test
Y = [DRYlag; WETlag; WRIlag];
G = [ones(size(DRYlag)); 2*ones(size(WETlag)); 3*ones(size(WRIlag))];
[p,tbl,stats] = anova1(Y,G,'off');
if p>=0.05
    fprintf('No significant effect of condition on lag: F(%d,%d)=%2.3f, p=%1.4f\n\n',...
        tbl{2,3},tbl{3,3},tbl{2,5},p);
else
    fprintf('Significant effect of condition on lag: F(%d,%d)=%2.3f, p=%1.4f\n\n',...
        tbl{2,3},tbl{3,3},tbl{2,5},p);
    stats.gnames = {'Dry';'Wet';'Wrinkly'};
    figure
    c = multcompare(stats)
end
[~,p,~,stats] = ttest(Y);
if p>=0.05
    fprintf('Lag not different to zero: t(%d)=%3.3f, p=%1.4f\n\n',...
        stats.df,stats.tstat,p);
else
    fprintf('Lag significantly different to zero: t(%d)=%3.3f, p=%1.4f, mean lag %4.3f ms\n\n',...
        stats.df,stats.tstat,p,mean(Y));
end

% Correlation of age and lag
[r,p] = corr(Y(I),A(I));
if p>=0.05
    fprintf('Lag does not correlate with age: r(%d)=%1.3f, p=%1.4f\n\n',...
        length(Y)-2,r,p);
else
    fprintf('Lag correlates significantly with age: r(%d)=%1.3f, p=%1.4f\n',...
        length(Y)-2,r,p);
    [B,~,~,~,stats] = regress (Y(I),[A(I) ones(size(A(I)))]);
    fprintf('Lag = %3.3f * Age + %3.3f. R2 = %2.3f\n\n',B(1),B(2),stats(1))
end


% Load force error - one-way analysis of variance
Y = [DRYerr; WETerr; WRIerr];
G = [ones(size(DRYerr)); 2*ones(size(WETerr)); 3*ones(size(WRIerr))];
[p,tbl,stats] = anova1(Y,G,'off');
if p>=0.05
    fprintf('No significant effect of condition on error: F(%d,%d)=%2.3f, p=%1.4f, mean r %1.4f\n\n',...
        tbl{2,3},tbl{3,3},tbl{2,5},p,mean(Y));
else
    fprintf('Significant effect of condition on error: F(%d,%d)=%2.3f, p=%1.4f, mean r %1.4f\n\n',...
        tbl{2,3},tbl{3,3},tbl{2,5},p,mean(Y));
    stats.gnames = {'Dry';'Wet';'Wrinkly'};
    figure
    c = multcompare(stats)
end


% Slope of GF rate of change - attack
Ya = [DRYatt; WETatt; WRIatt];
Ga = [ones(size(DRYatt)); 2*ones(size(WETatt)); 3*ones(size(WRIatt))];
disp('--- GF slope')
[p,tbl,stats] = anova1(Ya,Ga,'off')
% figure
[c,m] = multcompare(stats,'display','off')
% int_plotSlopeBars(Ya,Ga);


% Slope of GF rate of change - decay
Yd = [DRYdec; WETdec; WRIdec];
Gd = [ones(size(DRYdec)); 2*ones(size(WETdec)); 3*ones(size(WRIdec))];
disp('--- LF slope')
% figure
[p,tbl,stats] = anova1(Yd,Gd,'off')
% figure
[c,m] = multcompare(stats,'display','off')
% int_plotSlopeBars(Yd,Gd);

int_plotSlopeBars2(Ya,Ga,Yd,Gd);


% t-test for slope???
% [~,p,~,stats] = ttest(Y);
% if p>=0.05
%     fprintf('Lag not different to zero: t(%d)=%3.3f, p=%1.4f\n\n',...
%         stats.df,stats.tstat,p);
% else
%     fprintf('Lag significantly different to zero: t(%d)=%3.3f, p=%1.4f, mean lag %4.3f ms\n\n',...
%         stats.df,stats.tstat,p,mean(Y));
% end
%     
    

%% --- INTERNAL FUNCTIONS

% This function does most of the hard work.
function [GFm,LFm,RMSm,mlag,n,mAtt,mDec] = int_getdata (DATA)
LF = [];
GF = [];
RMS= [];
Att= [];
Dec= [];
[b,a] = butter (2, 20./500, 'low');
n = 0;
lags= [];
APHASE = 4001:6000;
DPHASE = 11001:13000;

for i=1:8
    D = DATA{i};
    lf = D.LF';
    if length(lf)<15000
        lf = [ones(1,(15000-length(lf)))*lf(1) lf];
    end
    if length(lf)>15000
        l = length(lf)-15000;
        lf = lf(1:end-l);
    end
    lff = filtfilt(b,a,lf);
    
    gf = D.GF';
    if length(gf)<15000
        gf = [ones(1,(15000-length(gf)))*gf(1) gf];
    end
    if length(gf)>15000
        l = length(gf)-15000;
        gf = gf(1:end-l);
    end
    gff = filtfilt(b,a,gf);
    
    [y,l]=xcorr(lf,gf,150); % NB max lags 150ms
    
    % Correlation of target and LF. Skips first 1.5 sec.
    T = int_getTargetTrace();
    R = corr(T(1501:end)',lff(1501:end)');
    
    % exclude trial if LF not greater than zero
    XX = 8501:10000;
    if (ttest(lff(XX),0,'Tail','right','Alpha',0.05))
        if (mean(gff)/mean(lff))<10
            % LF is larger, so combine with others
            LF = [LF; lff];
            GF = [GF; gff];
            RMS= [RMS; R];
            [~,mxi] = max(y);
            lags = [lags;l(mxi)];
            n=n+1;
        else
            %
        end
    else
        %
    end
    
    % slope of attack and decay phase - use simple linear regression
    B = regress (gff(APHASE)',[APHASE' ones(size(APHASE'))]);
    Att = [Att; B(1)*1000];
    B = regress (gff(DPHASE)',[DPHASE' ones(size(DPHASE'))]);
    Dec = [Dec; B(1)*1000];

end

LFm = mean(LF,1);
GFm = mean(GF,1);
RMSm= mean(RMS,1);
mlag = mean(lags);
mAtt = mean(Att);
mDec = mean(Dec);
% Att
% Dec


% Build the target load force trace
function T = int_getTargetTrace()
targettrace = [...
    5 * ones(1,3500), ...  % 3.5s ready
    5:(15/2999):20, ...    % 3s ramp up
    20 * ones(1,3999), ... % 4s hold
    20:-(15/3000):5, ...   % 3s ramp down
    5 * ones(1,1500)...    % 1.5s end
    ];
T = targettrace./10;

% Plot bar graph for slopes
function int_plotSlopeBars(allData,group)
    % A few common properties
    xCenter = [1 2 3];% 4 5 7 8 10 11 13 14];
    spread = 0.25; % 0=no spread; 0.5=random spread within box bounds (can be any value)
    txtCentres = [1 2 3];% 4.5 7.5 10.5 13.5];

    % Mean first
    figure; hold on;
    for i = 1:length(allData)
        %plot(rand(size(allData(i)))*spread -(spread/2) + xCenter(i), allData(i), 'k.','linewidth', 2)
        plot(rand(1)*spread -(spread/2) + xCenter(group(i)), allData(i), 'k.','linewidth', 2)
    end
    h = boxplot(allData,group,'Notch','on','positions',xCenter,...
        'Symbol','ko','OutlierSize',4);
    set(h, 'linewidth' ,2,'Color','k')
    set(gca,'XLim',[-.5 4.5])
%      set(gca,'YLim',[13 30])
%     set(gca,'YTick',14:2:28)
    set (gca,'XTick',txtCentres);
    set(gca,'XTickLabel', {'Dry';'Wet';'Wrinkly'})
%     TXT = {'**','**','ns','**','ns'};
%     y = 29;
%     for i=1:5
%         text (txtCentres(i),y,TXT{i},'FontSize',18,'HorizontalAlignment','center',...
%             'VerticalAlignment','middle');
%     end
    ylabel('Mean slope (N/sec)')

    
% Plot bar graph for slopes
function int_plotSlopeBars2(allData1,group1,allData2,group2)
xCenter = [1 2 3];% 4 5 7 8 10 11 13 14];
spread = 0.25; % 0=no spread; 0.5=random spread within box bounds (can be any value)
txtCentres = [1 2 3];% 4.5 7.5 10.5 13.5];
figure
%---
subplot(1,2,1)
hold on;
for i = 1:length(allData1)
    %plot(rand(size(allData(i)))*spread -(spread/2) + xCenter(i), allData(i), 'k.','linewidth', 2)
    plot(rand(1)*spread -(spread/2) + xCenter(group1(i)), allData1(i), 'k.','linewidth', 2)
end
h = boxplot(allData1,group1,'Notch','on','positions',xCenter,...
    'Symbol','ko','OutlierSize',3);
line ([.5 3.5],[.5 .5],'color','k','LineStyle','--','LineWidth',2)
set(h, 'linewidth' ,2,'Color','k')
set(gca,'XLim',[-.5 4.5])
      set(gca,'YLim',[-6 6])
%     set(gca,'YTick',14:2:28)
set (gca,'XTick',txtCentres);
set(gca,'XTickLabel', {'Dry';'Wet';'Wrinkly'})
%     TXT = {'**','**','ns','**','ns'};
%     y = 29;
%     for i=1:5
%         text (txtCentres(i),y,TXT{i},'FontSize',18,'HorizontalAlignment','center',...
%             'VerticalAlignment','middle');
%     end
ylabel('Mean slope (N/sec)')
%---
subplot(1,2,2)
hold on;
for i = 1:length(allData2)
    %plot(rand(size(allData(i)))*spread -(spread/2) + xCenter(i), allData(i), 'k.','linewidth', 2)
    plot(rand(1)*spread -(spread/2) + xCenter(group2(i)), allData2(i), 'k.','linewidth', 2)
end
h = boxplot(allData2,group2,'Notch','on','positions',xCenter,...
    'Symbol','ko','OutlierSize',3);
line ([.5 3.5],[-.5 -.5],'color','k','LineStyle','--','LineWidth',2)
set(h, 'linewidth' ,2,'Color','k')
set(gca,'XLim',[-.5 4.5])
%      set(gca,'YLim',[13 30])
%     set(gca,'YTick',14:2:28)
      set(gca,'YLim',[-6 6])
set (gca,'XTick',txtCentres);
set(gca,'XTickLabel', {'Dry';'Wet';'Wrinkly'})
%     TXT = {'**','**','ns','**','ns'};
%     y = 29;
%     for i=1:5
%         text (txtCentres(i),y,TXT{i},'FontSize',18,'HorizontalAlignment','center',...
%             'VerticalAlignment','middle');
%     end
%ylabel('Mean slope (N/sec)')
