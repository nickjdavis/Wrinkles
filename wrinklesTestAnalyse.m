% Run through folder to compute GF-LF ratio

function wrinklesTestAnalyse(folder)


dirin = cd;
cd (folder);

d = dir('*.mat');
n = length(d);

wrGF = [];
wrLF = [];
weGF = [];
weLF = [];
drGF = [];
drLF = [];

for i=1:n
    load (d(i).name);
    for j=1:nTrials
        LF = DATA{j}.LF;
        GF = DATA{j}.GF;
        switch(condName)
            case 'Wrinkly'
                wrGF = [wrGF; GF(1:14950)'];
                wrLF = [wrLF; LF(1:14950)'];
            case 'Wet'
                weGF = [weGF; GF(1:14950)'];
                weLF = [weLF; LF(1:14950)'];
            case 'Dry'
                drGF = [drGF; GF(1:14950)'];
                drLF = [drLF; LF(1:14950)'];
        end
    end
end

% size(wGF)
% size(wLF)

% mWrGF = mean(wrGF,1);
% mWrLF = mean(wrLF,1);
% mWeGF = mean(weGF,1);
% mWeLF = mean(weLF,1);
% mDrGF = mean(drGF,1);
% mDrLF = mean(drLF,1);
% plot ([mWrGF' mWrLF' mWeGF' mWeLF' mDrGF' mDrLF']);
% legend('Wrinkly GF','Wrinkly LF','Wet GF','Wet LF','Dry GF','Dry LF')

% Slight cheat for LF
mWrGF = mean(wrGF,1);
mWeGF = mean(weGF,1);
mDrGF = mean(drGF,1);
mLF = mean([wrLF;weLF;drLF],1);
plot ([mWrGF' mWeGF' mDrGF' mLF']);
legend('Wrinkly GF','Wet GF','Dry GF','LF','Location','NorthWest')
xlabel('Time (ms)')
ylabel('Force (??)')

cd (dirin)