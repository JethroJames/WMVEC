function [jaccard] = update_jaccard(view,lambda,R,dis)
part_1 = 0;
for i= 1:view
    part_1 = part_1 + R(i)*sum(sum(dis{i}));
end
part_2 = 0;
for i=1:length(R)
    part_2 = part_2 + lambda .* sum((R(i)+1e-4).*log(R(i)+1e-3));
end

jaccard = part_1 + part_2;
end