function result = CalcMeasures(Y, predY)
% result = [ACC,error_cnt];
if size(Y,2) ~= 1
    Y = Y';
end

if size(predY,2) ~= 1
    predY = predY';
end

% bestMap
predY = bestMap(Y, predY);
if size(Y)~=size(predY)
    predY=predY';
end

error_cnt = sum(Y ~= predY);
AC = length(find(Y == predY))/length(Y);

result = [AC, error_cnt];
