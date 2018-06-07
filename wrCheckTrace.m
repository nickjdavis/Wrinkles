% check conflict?

targettrace = [...
    10 * ones(1,3000), ...  % 4s ready
    10:(10/3000):25, ...    % 3s ramp up
    25 * ones(1,5000), ...  % 4s hold
    25:-(10/3000):10, ...   % 3s ramp down
    10 * ones(1,1000)...    % 1s end
    ];


% but in each loop...

Screen('DrawLine', expWin, [0 0 255], mx-500, my+200, floor(mx-(500-500*(4/7.5))), my+200, 2);
Screen('DrawLine', expWin, [0 0 255], floor(mx-(500-500*(4/7.5))), my+200, floor(mx-(500-500*(7/7.5))), my, 2);
Screen('DrawLine', expWin, [0 0 255], floor(mx-(500-500*(7/7.5))), my, floor(mx+(500*(3/7.5))), my, 2);
Screen('DrawLine', expWin, [0 0 255], floor(mx+(500*(3/7.5))), my, floor(mx+(500*(6/7.5))), my+200, 2);
Screen('DrawLine', expWin, [0 0 255], floor(mx+(500*(6/7.5))), my+200, floor(mx+(500*(7/7.5))), my+200, 2);

% ... do these give the save trace?