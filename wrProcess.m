function wrProcess (folder)

dirin = cd;
cd(folder)

d = dir('*.mat');
n = length(d);

WetGF = []; WetLF = []; WetPC = [];
DryGF = []; DryLF = []; DryPC = [];
WriGF = []; WriLF = []; WriPC = [];



[b,a] = butter (2, 10./500, 'low');

for i=1:n
    load(d(i).name);
    
    LL = []; GG = []; PP = [];
    
    for j=1:8
        
        L = DATA{j}.LF';
        if length(L)<15000
            %size(L)
            L = filtfilt(b,a,[ones(1,(15000-length(L)))*L(1) L]);
        end
        LL = [LL; L(1:15000)];
        
        G = DATA{j}.GF';
        if length(G)<15000
            G = filtfilt(b,a,[ones(1,(15000-length(G)))*G(1) G]);
        end
        GG = [GG; G(1:15000)];
    
        PP = [PP; mean(G(8000:12000))/mean(L(8000:12000))];
    end
    
    %disp(condName)
    switch condName
        case 'Dry'
            DryGF = [DryGF; mean(GG,1)];
            DryLF = [DryLF; mean(LL,1)];
            DryPC = [DryPC; mean(PP)];
        case 'Wet'
            WetGF = [WetGF; mean(GG,1)];
            WetLF = [WetLF; mean(LL,1)];
            WetPC = [WetPC; mean(PP)];
        case 'Wrinkly'
            WriGF = [WriGF; mean(GG,1)];
            WriLF = [WriLF; mean(LL,1)];
            WriPC = [WriPC; mean(PP)];
    end
end


% Mean trace
mWtraceGF = mean(WetGF,1); mWtraceLF = mean(WetLF,1);
mDtraceGF = mean(DryGF,1); mDtraceLF = mean(DryLF,1);
mRtraceGF = mean(WriGF,1); mRtraceLF = mean(WriLF,1);


% Percentage overgrip
% bar ([mean(DryPC(6001:12000)) mean(WetPC(6001:12000)) mean(WriPC(6001:12000))]);
% legend('Dry','Wet','Wri')
figure;
barweb ([mean(DryPC) mean(WetPC) mean(WriPC)],...
    [std(DryPC) std(WetPC) std(WriPC)]);
legend('Dry','Wet','Wri')

figure
plot ([mWtraceLF' mWtraceGF' mDtraceLF' mDtraceGF' mRtraceLF' mRtraceGF']);
legend('Wet LF','Wet GF','Dry LF','Dry GF','Wri LF','Wri GF')
xlabel('Time (milliseconds)')
ylabel('Force (not scaled)')



cd (dirin)