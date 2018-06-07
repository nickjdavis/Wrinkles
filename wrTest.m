

W = load ('wrinkles_435');
D = load ('wrinkles_434');


WtraceLF = [];
WtraceGF = [];
DtraceLF = [];
DtraceGF = [];

for i=1:8
    x = W.DATA{i}.LF';
    if length(x)<15000
        size(x)
        x = [ones(1,(15000-length(x)))*x(1) x];
    end
    WtraceLF = [WtraceLF; x(1:15000)];
    x = W.DATA{i}.GF';
    if length(x)<15000
        x = [ones(1,(15000-length(x)))*x(1) x];
    end
    WtraceGF = [WtraceGF; x(1:15000)];
    x = D.DATA{i}.LF';
    if length(x)<15000
        x = [ones(1,(15000-length(x)))*x(1) x];
    end
    DtraceLF = [DtraceLF; x(1:15000)];
    x = D.DATA{i}.GF';
    if length(x)<15000
        x = [ones(1,(15000-length(x)))*x(1) x];
    end
    DtraceGF = [DtraceGF; x(1:15000)];
end

mWtraceLF = mean(WtraceLF,1);
mWtraceGF = mean(WtraceGF,1);
mDtraceLF = mean(DtraceLF,1);
mDtraceGF = mean(DtraceGF,1);

figure
plot ([mWtraceLF' mWtraceGF' mDtraceLF' mDtraceGF']);
legend('Wet LF','Wet GF','Dry LF','Dry GF')

