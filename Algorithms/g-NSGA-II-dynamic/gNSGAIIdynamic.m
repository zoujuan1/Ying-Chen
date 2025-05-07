function gNSGAIIdynamic(Global)
% <algorithm> <A-G>
% g-Dominance: Reference Point Based Dominance for Multiobjective
% Metaheuristics

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    Point = [0.2,0.8];

    %% Generate random population
    Population = Global.Initialization();
    FrontNo    = NDSort(Evaluate(Population.objs,Point),inf);
    CrowdDis   = CrowdingDistance(Population.objs,FrontNo);
    
    filename = strcat(func2str(Global.algorithm),func2str(Global.problem),'_',num2str(Global.gen));
        out_file=['SelfData\',filename,'.txt'];
        fid = fopen(out_file,'wt');
        data = Population.objs;
        [row,col] = size(data);
        %输出种群
        for i = 1 : row
            for j = 1 : col
                if j == col
                    fprintf(fid,'%g\n',data(i,j));
                else
                    fprintf(fid,'%g\t',data(i,j));
                end
            end
        end
        fclose(fid);

    %% Optimization
    while Global.NotTermination(Population)
      %% Update the reference point
       switch rem(ceil(Global.gen/50),7)
            case 0
                Point = [0.2,0.8];
            case 1
                Point = [0.6,0.4];
            case 2
                Point = [0.3,0.7];
            case 3
                Point = [0.7,0.3];
            case 4
                Point = [0.4,0.6];
            case 5
                Point = [0.8,0.2];
            case 6
                Point = [0.5,0.5];
        end
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N,Point);
        
         %% 输出每次变化的种群
        filename = strcat(func2str(Global.algorithm),func2str(Global.problem),'_',num2str(Global.gen));
        out_file=['SelfData\',filename,'.txt'];
        fid = fopen(out_file,'wt');
        data = Population.objs;
        [row,col] = size(data);
        %输出种群
        for i = 1 : row
            for j = 1 : col
                if j == col
                    fprintf(fid,'%g\n',data(i,j));
                else
                    fprintf(fid,'%g\t',data(i,j));
                end
            end
        end
        fclose(fid);
    end
end