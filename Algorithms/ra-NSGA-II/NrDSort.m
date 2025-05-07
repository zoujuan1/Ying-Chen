function [FrontNo,MaxFNo] = NrDSort(PopObj,nSort,Points,W,delta)
% Do non-r-dominated sorting by efficient non-dominated sort (ENS)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    FrontNo = inf(1,size(PopObj,1));
    for i = 1 : size(Points,1)
        FrontNo = min(FrontNo,nrdsort(PopObj,Points(i,:),W(i,:),delta));
    end
    MaxFNo = find(cumsum(hist(FrontNo,1:max(FrontNo)))>=min(nSort,length(FrontNo)),1);
    FrontNo(FrontNo>MaxFNo) = inf;
end

function FrontNo = nrdsort(PopObj,g,w,delta)
% Sort the population according to one preferred point

    [PopObj,~,Loc] = unique(PopObj,'rows');
    % Calculate the weighted Euclidean distance of each solution
    %Dist       = sqrt(((PopObj-repmat(g,size(PopObj,1),1))./repmat(max(PopObj,[],1)-min(PopObj,[],1),size(PopObj,1),1)).^2*(w/sum(w))');
    %DistExtent = max(Dist) - min(Dist);
    
    % ra-dominace
    dist = pdist2(g,PopObj);
    [~,DistTag] = sort(dist);
    Mindist = min(dist);
    MinPop = PopObj(DistTag(1));
    r = Mindist * tan(delta);
    Dist1 = sqrt(sum((PopObj - repmat(g,size(PopObj,1),1)).^2,2));
    Dist2 = sqrt(sum((PopObj - repmat(MinPop,size(PopObj,1),1)).^2,2));
    p = (Dist1 + Dist2 + Mindist)./2;
    Dist = 2 * sqrt(p.*(p-Dist1).*(p-Dist2).*(p-Mindist))./Mindist;
    
    % Sort the population based on their Dist values, so that a solution
    % cannot r-dominate the solutions having smaller Dist values than it
    [Dist,rank] = sort(Dist);
    PopObj      = PopObj(rank,:);
    % Non-r-dominated sorting by ENS
    [N,M]   = size(PopObj);
    FrontNo = inf(1,N);
    MaxFNo  = 0;
    while any(FrontNo==inf)
        MaxFNo = MaxFNo + 1;
        for i = 1 : N
            if FrontNo(i) == inf
                Dominated = false;
                for j = i-1 : -1 : 1
                    if FrontNo(j) == MaxFNo
                        m = 1;
                        while m <= M && PopObj(i,m) >= PopObj(j,m)
                            m = m + 1;
                        end
                        Dominated = m > M;
                        if ~Dominated
                            %Dominated = (Dist(j)-Dist(i))./DistExtent < -delta;
                            Dominated = (Dist(i)-Dist(j)) > r;
                        end
                        if Dominated
                            break;
                        end
                    end
                end
                if ~Dominated
                    FrontNo(i) = MaxFNo;
                end
            end
        end
    end
    FrontNo(rank) = FrontNo;
    FrontNo = FrontNo(Loc);
end