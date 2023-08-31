function [R] = update_R(view,dis,lambda)
R = repmat(1/view,[1,view]);
temp = [];

for i = 1:view
    F_i = sum(sum(dis{i},1));
    temp = [temp exp(-F_i/lambda)];
end

for i = 1:view
    R(i) = temp(i)/sum(temp);
end

end

