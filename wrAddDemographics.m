% function to add demographoc data to files

function wrAddDemographics (folder)


% open info file
% first col comes out as 'x___Number' (???)
infofile = 'WrinklesInfo.csv';
T = readtable(infofile);

% run through data files
dirin = cd;
cd (folder);
d = dir('*.mat');
n = length(d);

for i=1:n
    datafile = d(i).name;
    % this is ugly - matches digits (returns string)
    subjectnumberstr = datafile(regexp(datafile,'\d'));
    subjectnumber = str2num(subjectnumberstr);
    x = find(T.x___Number==subjectnumber);
%     disp(T.Age(x))
%     disp(T.Sex(x))
%     disp(T.Hand(x))
    Age = T.Age(x); % NB could be NaN
    Sex = T.Sex(x);
    Hand= T.Hand(x);
    save(datafile, 'Age','Sex','Hand','-append')
end

cd(dirin)