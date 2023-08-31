function [centroids] = Centroids_Initialization(X,K)
centroids = [];
for i = 1:K
    if i == 1
        d = [];
        for j = 1:size(X,1)
            center = X(j,:);dis = 0;
            for z = 1:size(X,1)
                % Calculate the distance except for the centroid
                if z~=j
                    dis = dis + norm(center-X(z,:),2);
                end
            end
            d = [d;dis];
        end
        [~,idx] = min(d);
        centroids = [centroids;X(idx,:)];
        X(idx,:) = [];
    else
        if size(centroids,1) == 1
           for j = 1:size(centroids,1)
               center = centroids(j,:);
               dis = sum((X-center).^2,2);
               [~,idx] = max(dis);
               centroids = [centroids;X(idx,:)];
               X(idx,:) = [];
           end
        else 
            Dis = [];Idx = [];
            for j = 1:size(centroids,1)
                center = centroids(j,:);
               dis = sum((X-center).^2,2);
               Dis = [Dis;sum(dis)];
               [~,idx] = max(dis);
               Idx = [Idx;idx];
            end
            [~,idx] = min(Dis);
           centroids = [centroids;X(Idx(idx),:)];
           X(Idx(idx),:) = [];
        end
    end

end
end

