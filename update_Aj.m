function [Aj,new_center] = update_Aj(view,cluster,features,Aj,F,center,nbFoc)
% Aj Initialization %
for i=1:view
    new_center = zeros([nbFoc,features(i)]);
    for j=1:length(F)
        if sum(F(j,:)) ~= 0
            temp1 = 0;
            for k=1:cluster
               temp1 = temp1 + F(j,k) .* center{i}(k,:);
            end
            new_center(j,:) = temp1./sum(F(j,:));
        end   
    end
    Aj{i} = new_center;
end
end