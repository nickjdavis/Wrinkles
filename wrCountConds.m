% count no of people in each cond

function [Dr,We,Wr] = wrCountConds (folder)

dirin = cd;
cd (folder);
d = dir('*.mat');
n = length(d);

Dr = 0;
We = 0;
Wr = 0;

for i=1:n
    load (d(i).name);
    switch condName
        case 'Dry'
            Dr = Dr+1;
        case 'Wet'
            We = We+1;
        case 'Wrinkly'
            Wr = Wr+1;
    end
end


cd (dirin)