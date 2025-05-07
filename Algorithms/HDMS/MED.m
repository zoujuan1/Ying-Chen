function med = MED(PopObj)
% Calculate the crowding distance of each solution front by front

%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [N,M]    = size(PopObj);
    med = zeros(1,N);
    Med=1000000;
        for i = 1 : length( PopObj)
            TotalDis=0;   
            NearDis=inf;
            Dis=0;
            for j = 1 : length(PopObj)
                if i==j
                    continue
                end
                A=PopObj(i,:)-PopObj(j,:);
                Dis=norm(A,2);
                TotalDis=TotalDis+Dis;
                if Dis<NearDis
                    NearDis=Dis;
                end
            end
            TotalDis;
            NearDis;
            Med=NearDis*TotalDis;
            med(i)=Med;
        end
       
    
end