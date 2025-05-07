function raNSGAII(Global)
% <algorithm> <A>
% r-dominance based NSGA-II
% Points ---     --- Set of preferred points
% W      ---     --- Set of weight vector for each preferred point
% alpha  --- 0.1 --- Non-ra-dominance angle

%------------------------------- Reference --------------------------------
% L. B. Said, S. Bechikh, and K. Ghedira, The r-dominance: A new dominance
% relation for interactive evolutionary multicriteria decision making, IEEE
% Transactions on Evolutionary Computation, 2010, 14(5): 801-818.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [Points,W,alpha] = Global.ParameterSet(zeros(1,Global.M)+0.5,ones(1,Global.M),0.1);

    %% Generate random population
    Population = Global.Initialization();
    FrontNo    = NrDSort(Population.objs,inf,Points,W,alpha);
    CrowdDis   = CrowdingDistance(Population.objs,FrontNo);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
         Offspring  = Global.Variation(Population(MatingPool));
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N,Points,W,alpha);
%          if Global.gen == Global.maxgen
%             filename = strcat(func2str(Global.algorithm),'_',class(Global.problem),'_');
%             out_file=[filename,'.txt'];
%             fid = fopen(out_file,'wt');
%             data = Population.objs;
%             [row,col] = size(data);
%             % ‰≥ˆ÷÷»∫
%             for i = 1 : row
%                 for j = 1 : col
%                     if j == col
%                         fprintf(fid,'%g\n',data(i,j));
%                     else
%                         fprintf(fid,'%g\t',data(i,j));
%                     end
%                 end
%             end
%             fclose(fid);
%         end
    end
end