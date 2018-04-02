% TODO
% Resample LF target trace
% Scale data with cal file
% Finish save


function wrinklesCheck(folder)

[b,a] = butter(2,15/500,'low');

dirin = cd;
cd (folder);
d = dir('*.mat');
n = length(d);

for i=1:n
    load (d(i).name);
    includedTrials = zeros(nTrials,1);
    for j=1:nTrials
        LF = DATA{j}.LF;
        GF = DATA{j}.GF;
        LFf = filtfilt(b,a,LF);
        GFf = filtfilt(b,a,GF);
        % plot data
        plot ([GFf,LFf]);
        hold on
        % get button press
        waitforbuttonpress
        S = get(gcf,'SelectionType');
        switch S
            case 'normal'
                % keep
                includedTrials(j) = 1;
            case 'alt'
                % reject
                includedTrials(j) = 0;
        end
        hold off
    end
    % mean of included trials
    includedTrials
    % save back to file
    
end


cd (dirin)