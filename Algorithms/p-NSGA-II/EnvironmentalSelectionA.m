function [Population] = EnvironmentalSelectionA(Population,M,N,Point,delta,gen,dpj_old,dpj)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
     %%
     if gen > 1
       Pop_num = N + N;
     else
        Pop_num = N;
     end
     dist = zeros(1,Pop_num);
     for i = 1 : Pop_num         
         if M == 3
            dist(i) = norm(cross(Point,Population(i).obj))/norm(Point);
        elseif M == 2
            dist(i) = abs(det([Point;Population(i).obj]))/norm(Point);
        else
            dist(i) = pdist2(Population(i).obj,Point);
        end
     end
     dpj = mean(dist);
     if gen == 1
         dpj =  Inf;
     elseif dpj <= delta
         dpj = delta;
     elseif (dpj <= dpj_old) && (dpj > delta)
         dpj = dpj;
     else
         dpj = dpj_old;
     end
     dpj_old = dpj;
     P = zeros(1,Pop_num);
     L = zeros(1,Pop_num);
     count_p = 1;
     count_l = 1;
     for i = 1 : Pop_num
         if dist(i) <= dpj
             P(count_p) = i;
             count_p = count_p + 1;
         else
             L(count_l) = i;
            dist(count_l) = dist(i);
             count_l = count_l + 1;
         end
     end
     count_p = count_p - 1;
     if count_p <= N
         [~,dis] = sort(dist);
         aaa = 1;
         while count_p < N       
             count_p = count_p + 1;
             P(count_p) = dis(aaa);
             aaa = aaa + 1;
         end
         Population = Population(P(1:count_p));
     else  
      Population = Population(P(1:count_p));    
    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population().objs,Population().cons,N);
    Next = FrontNo < MaxFNo;
    
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(Population().objs,FrontNo);
    
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
     end
end