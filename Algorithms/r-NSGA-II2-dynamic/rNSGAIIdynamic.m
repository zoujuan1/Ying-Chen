function rNSGAIIdynamic(Global)
% <algorithm> <O-Z>
% The r-Dominance: A New Dominance Relation for Interactive Evolutionary
% Multicriteria Decision Making
% W      ---     --- Set of weight vector for each preferred point
% delta  --- 0.1 --- Non-r-dominance threshold

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    [W,delta] = Global.ParameterSet(ones(1,Global.M),0.1);
    Points = zeros(1, Global.M);
    t =0.1* floor(Global.gen/Global.tao_n);
    Points(1) = (cos(pi*t)).^2 + 0.1;
    for i = 2 : Global.M
        Points(i) = (sin(pi*t)).^2 + 0.1;
    end

    %% Generate random population
    Population = Global.Initialization();
    FrontNo    = NrDSort(Population.objs,inf,Points,W,1-(1-delta)*Global.gen/Global.maxgen);
    CrowdDis   = CrowdingDistance(Population.objs,FrontNo);
    
    %% Optimization
    while Global.NotTermination(Population)    
         %% Update the reference point
       t =0.1* floor(Global.Current()/Global.tao_n);
       Points(1) = (cos(pi*t)).^2 + 0.1;
       for i = 2 : Global.M
           Points(i) = (sin(pi*t)).^2 + 0.1;
       end
       
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        Offspring  = Global.Variation(Population(MatingPool));
        [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N,Points,W,1-(1-delta)*Global.gen/Global.maxgen);
        if rem(Global.Current,Global.tao_n) == 0
          for j = 1 : Global.N
            Population(j) = INDIVIDUAL(Population(j).dec);
          end
        end
        
                %% 输出每次变化的种群
          if rem(Global.Current+1,Global.tao_n) == 0
                filename = strcat(func2str(Global.algorithm),'_',func2str(Global.problem),'_t',num2str(floor(Global.Current()/Global.tao_n)),'_run',num2str(Global.run));
                out_file=[filename,'.txt'];
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
end