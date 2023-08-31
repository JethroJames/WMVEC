function [V] = update_V(view,cluster,alpha,beta,data,M,R,F_update,features)
V = cell(1,view);
for i = 1:view
    v_i = zeros([cluster,features(i)]);
    for p = 1:features(i)
        B = zeros([cluster,1]); % shape:Cx1
        %% Calculate B
        for j = 1:cluster
            pos = []; % Record the indices
            for k = 1:length(F_update{i})
                if F_update{i}(k,j) == 1
                    pos(end+1) = k;
                end
            end
            B_X = 0;
            for n=1:length(pos)
                card = sum(F_update{i}(pos(n),:));
                r = R(i);
                aj_dis = card^(alpha-1)*r .* data{i}(:,p) .* (M(:,pos(n)).^beta);
                B_X = B_X + sum(aj_dis);
            end
            B(j) = B_X;
        end
        %%  Calculate H
        H = zeros([cluster,cluster]);
        for c=1:cluster
            for k=1:cluster
                loc = [];
                for n = 1:length(F_update{i})
                    if F_update{i}(n,c) == 1 && F_update{i}(n,k)== 1
                        loc(end+1) = n;
                    end
                end
                H_ck = 0;
                for n=1:length(loc)
                    card = sum(F_update{i}(loc(n),:));
                    r = R(i);
                    aj_all_dis = card^(alpha-2) * r .* M(:,loc(n)).^beta;
                    H_ck = H_ck + sum(aj_all_dis);
                end
                H(c,k) = H_ck;
            end
        end
        v = pinv(H) * B;
        v_i(:,p) = v;
    end
    V{i} = v_i;
end
end