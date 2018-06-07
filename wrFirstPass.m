function wrFirstPass (folder)

dirin = cd;
cd(folder)

d = dir('*.mat');
n = length(d);

[b,a] = butter (2, 10./500, 'low');

for i=1:n
    load(d(i).name);
    AcceptedTrials = [];
    for j=1:8
        % plot trial
        LFf = filtfilt(b,a,DATA{j}.LF);
        GFf = filtfilt(b,a,DATA{j}.GF);
        plot ([GFf LFf]); hold on;
        legend ('GF','LF')
        % button for acc/rej
        cA = rectangle ('Position',  [100 1 1000 1],'Curvature',[1 1]);
        cR = rectangle ('Position',[14000 1 1000 1],'Curvature',[1 1]);
        set(cA,'FaceColor','g')
        set(cR,'FaceColor','r')
        % save acc/rej list
        waitforbuttonpress
        hold off
    end
end


cd (dirin)