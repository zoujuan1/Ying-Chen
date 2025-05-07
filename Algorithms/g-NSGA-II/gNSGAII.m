function gNSGAII(Global)
% <algorithm> <A-G>
% g-Dominance: Reference Point Based Dominance for Multiobjective
% Metaheuristics
% Point --- --- Preferred point

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    Point = Global.ParameterSet(zeros(1,Global.M)+0.5);

    %% Generate random population
    Population = Global.Initialization();
    FrontNo    = NDSort(Evaluate(Population.objs,Point),inf);
    CrowdDis   = CrowdingDistance(Population.objs,FrontNo);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N,Point);
        
%         if Global.gen == Global.maxgen
%         filename = strcat(func2str(Global.algorithm),'_',func2str(Global.problem),'_');
%         out_file=[filename,'.txt'];
%         fid = fopen(out_file,'wt');
%         data = Population.objs;
%         [row,col] = size(data);
%         % ‰≥ˆ÷÷»∫
%         for i = 1 : row
%             for j = 1 : col
%                 if j == col
%                     fprintf(fid,'%g\n',data(i,j));
%                 else
%                     fprintf(fid,'%g\t',data(i,j));
%                 end
%             end
%         end
%         fclose(fid);
%         end
    end
end