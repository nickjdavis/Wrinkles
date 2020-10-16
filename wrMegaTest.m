% This is some VERY BAD CODE
function wrMegaTest (folderlist)


nfolders = length(folderlist);

% PI = load('persistentInfo');
     PI.LF_LC_offset_V = 0;
           PI.VperN_LF= 0.0275;
     PI.GF_LC_offset_V= 0.0272;
           PI.VperN_GF= 0.0222;




Dr = 0;
We = 0;
Wr = 0;

DRYLF = [];
WETLF = [];
WRILF = [];
DRYGF = [];
WETGF = [];
WRIGF = [];



for f = 1:nfolders
    
    folder = folderlist{f}
    dirin = cd;
    cd (folder);
    d = dir('*.mat');
    n = length(d);
    
    for i=1:n
        
        load (d(i).name);
        switch condName
            case 'Dry'
                Dr = Dr+1;
                [G,L] = int_getdata(DATA,PI);
                DRYLF = [DRYLF; L];
                DRYGF = [DRYGF; G];
            case 'Wet'
                We = We+1;
                [G,L] = int_getdata(DATA,PI);
                WETLF = [WETLF; L];
                WETGF = [WETGF; G];
            case 'Wrinkly'
                Wr = Wr+1;
                [G,L] = int_getdata(DATA,PI);
                WRILF = [WRILF; L];
                WRIGF = [WRIGF; G];
        end
    end
    
    
    cd (dirin)
    
end

% disp([Dr We Wr])

sprintf('Wet: %d\nDry: %d\nWri: %d\n',We,Dr,Wr)

DRYLFm = mean(DRYLF,1);
WETLFm = mean(WETLF,1);
WRILFm = mean(WRILF,1);
DRYGFm = mean(DRYGF,1);
WETGFm = mean(WETGF,1);
WRIGFm = mean(WRIGF,1);

figure
plot ([DRYLFm' WETLFm' WRILFm'],'LineWidth',2);
hold on;
plot ([DRYGFm' WETGFm' WRIGFm']);
legend('Dry','Wet','Wri')

% size(DRYLF)

DRYLFseg = DRYLF(:,8001:10000);
WETLFseg = WETLF(:,8001:10000);
WRILFseg = WRILF(:,8001:10000);
DRYGFseg = DRYGF(:,8001:10000);
WETGFseg = WETGF(:,8001:10000);
WRIGFseg = WRIGF(:,8001:10000);

% size(DRYLFseg)
% m = mean(DRYGFseg,2);
% size(m)
DryPC = mean(DRYGFseg,2)./mean(DRYLFseg,2);
WetPC = mean(WETGFseg,2)./mean(WETLFseg,2);
WriPC = mean(WRIGFseg,2)./mean(WRILFseg,2);
% size(DRYpc)

figure;
barweb ([mean(DryPC) mean(WetPC) mean(WriPC)],...
    [std(DryPC) std(WetPC) std(WriPC)]);
legend('Dry','Wet','Wri')

Y = [DryPC; WetPC; WriPC];
G = [ones(size(DryPC)); 2*ones(size(WetPC)); 3*ones(size(WriPC))];
[p,tbl,stats] = anova1(Y,G)

c = multcompare(stats)

function [GF,LF] = int_getdata (DATA,P)
LF = [];
GF = [];
[b,a] = butter (2, 10./500, 'low');

for i=1:8
    D = DATA{i};
    x = D.LF';
    if length(x)<15000
        x = [ones(1,(15000-length(x)))*x(1) x];
    end
    if length(x)>15000
        l = length(x)-15000;
        x = x(1:end-l);
    end
    xf = filtfilt(b,a,x);
    xfn= xf; %int_VtoN(xf,P.LF_LC_offset_V,P.VperN_LF);
    LF = [LF; xfn];
    
    y = D.GF';
    if length(y)<15000
        y = [ones(1,(15000-length(y)))*y(1) y];
    end
    if length(y)>15000
        l = length(y)-15000;
        y = y(1:end-l);
    end
    yf = filtfilt(b,a,y);
    yfn= yf; %int_VtoN(yf,P.GF_LC_offset_V,P.VperN_GF);
    GF = [GF; yfn];
end


function N = int_VtoN (V,offset,VperN)
N = (V-offset)/VperN;