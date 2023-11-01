function cmd = stepCmd(t, stepTimes, stepValues)
cmd = 0;
if nargin <= 1
    return
end

idx = sum(t >= stepTimes);
if idx == 0
    return
end
cmd = stepValues(idx);
end