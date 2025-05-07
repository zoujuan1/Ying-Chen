function preNSGAII(Global)
% <algorithm> <A>
% Nondominated sorting genetic algorithm II
% point --- --- Preferred point
% r --- 0.1 ---- preferred region

%------------------------------- Reference --------------------------------
% a unified approach for MOEA to solve preference-based MOPs
%--------------------------------------------------------------------------

    %% Parameter setting
    [point,r] = Global.ParameterSet(zeros(1,Global.M)+0.5,0.1);
    %point = [0.2,0.7];
    
    %% Generate random population
    Population = Global.Initialization();
     alpha = 1;  
    beta = 5;
     % 惩罚系数 proportion
        if Global.gen <= 0.7*Global.maxgen     % proportion  = a + b/N*(1+k/T)^cp, cp = log(N(b-a)/b)/log1.6
            proportion = alpha + beta/Global.N * (1 + Global.gen/Global.maxgen)^(log((beta-alpha)*Global.N/beta)/log(1.6));
        else                                % proportion = b
            proportion = beta;
        end        
     %% 距离惩罚
        best_point = point - r;
        worst_point = point + r;
        Disx = zeros(Global.N,1);
        if Global.M == 3
            MatrixA = cross(repmat(best_point - worst_point,Global.N,1),Population.objs - repmat(worst_point,Global.N,1));
            for j = 1 : Global.N
                Disx(j) = norm(MatrixA(j,1:end))/norm(best_point - worst_point);
            end
        elseif Global.M == 2
            MatrixB = Population.objs;
            for j = 1 : Global.N
                Disx(j) = abs(det([best_point-worst_point;MatrixB(j,1:end)-worst_point]))/norm(best_point-worst_point);
            end
        else
            Disx = pdist2(Population.objs,best_point - worst_point);
        end
        MaxDisx = max(Disx);       
        Gx = Population.objs - repmat(point,Global.N,1);
        GX = abs(Gx);
        GX_max = max(GX,[],1);
        for j = 1 : Global.M
            if GX_max(j) == 0
                GX_max(j) = 1;
            end
        end       
        gp = Population.objs - repmat(best_point,Global.N,1);
        ROI = false;    %偏好区域位置
        for k = 1 : Global.N
            if isempty(find(gp(k,1:end) > 0))  %ROI在可行域内
                ROI = true;
                break;
            end
        end          
        for k = 1 : Global.N  % 在ROI内的惩罚为0
            if isempty(find(abs(Gx(k,1:end)) > r))
                GX(k,1:end) = 0;
            end
        end       
        if ROI == true   
            for k = 1 : Global.N
                if isempty(find(Gx(k,1:end) >= 0))  % 对于ROI在可行域，支配参考点的区域惩罚为0
                    GX(k,1:end) = 0;
                end
            end
        end
        % 目标单位化
        ideal_p = min(Population.objs,[],1);
        worst_p = max(Population.objs,[],1);
        FixObj = (Population.objs-repmat(ideal_p,Global.N,1))./(repmat(worst_p - ideal_p,Global.N,1)) + proportion*(GX./repmat(GX_max,Global.N,1)) + 1/proportion*Disx./MaxDisx;
        %FixObj = (Population.objs-repmat(ideal_p,Global.N,1))./(repmat(worst_p - ideal_p,Global.N,1)) + proportion*(GX./repmat(GX_max,Global.N,1));
        %FixObj = (Population.objs-repmat(ideal_p,Global.N,1))./(repmat(worst_p - ideal_p,Global.N,1)) + 5*(GX./repmat(GX_max,Global.N,1)) + 1*Disx./MaxDisx;
        
        [~,FrontNo,CrowdDis] = pEnvironmentalSelection(Population,Global.N,FixObj);
        
       filename = strcat(func2str(Global.algorithm),func2str(Global.problem),'_',num2str(Global.gen));
        out_file=['0.7\',filename,'.txt'];
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
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        
        % 惩罚系数 proportion
        if Global.gen <= 0.7*Global.maxgen     % proportion  = a + b/N*(1+k/T)^cp, cp = log(N(b-a)/b)/log1.6
            proportion = alpha + beta/Global.N * (1 + Global.gen/Global.maxgen)^(log((beta-alpha)*Global.N/beta)/log(1.6));
        else                                % proportion = b
            proportion = beta;
        end
        
     %% 距离惩罚
        if Global.M == 3
            MatrixA = cross(repmat(best_point - worst_point,Global.N,1),Population.objs - repmat(worst_point,Global.N,1));
            for j = 1 : Global.N
                Disx(j) = norm(MatrixA(j,1:end))/norm(best_point - worst_point);
            end
        elseif Global.M == 2
            MatrixB = Population.objs;
            for j = 1 : Global.N
                Disx(j) = abs(det([best_point-worst_point;MatrixB(j,1:end)-worst_point]))/norm(best_point-worst_point);
            end
        else
            Disx = pdist2(Population.objs,best_point - worst_point);
        end
        MaxDisx = max(Disx);
        
        Gx = Population.objs - repmat(point,Global.N,1);
        GX = abs(Gx);
        
        GX_max = max(GX,[],1);
        for j = 1 : Global.M
            if GX_max(j) == 0
                GX_max(j) = 1;
            end
        end
        
        gp = Population.objs - repmat(best_point,Global.N,1);
        ROI = false;    %偏好区域位置
        for k = 1 : Global.N
            if isempty(find(gp(k,1:end) > 0))  %ROI在可行域内
                ROI = true;
                break;
            end
        end
           
        for k = 1 : Global.N  % 在ROI内的惩罚为0
            if isempty(find(abs(Gx(k,1:end)) > r))
                GX(k,1:end) = 0;
            end
        end
        
        if ROI == true   
            for k = 1 : Global.N
                if isempty(find(Gx(k,1:end) >= 0))  % 对于ROI在可行域，支配参考点的区域惩罚为0
                    GX(k,1:end) = 0;
                end
            end
        end
        
        % 目标单位化
        ideal_p = min(Population.objs,[],1);
        worst_p = max(Population.objs,[],1);
        FixObj = (Population.objs-repmat(ideal_p,Global.N,1))./(repmat(worst_p - ideal_p,Global.N,1)) + proportion*(GX./repmat(GX_max,Global.N,1)) + 1/proportion*Disx./MaxDisx;
        %FixObj = (Population.objs-repmat(ideal_p,Global.N,1))./(repmat(worst_p - ideal_p,Global.N,1)) + proportion*(GX./repmat(GX_max,Global.N,1)) ;
        %FixObj = (Population.objs-repmat(ideal_p,Global.N,1))./(repmat(worst_p - ideal_p,Global.N,1)) + 5*(GX./repmat(GX_max,Global.N,1)) + 1*Disx./MaxDisx;
        
        
        Offspring  = Global.Variation(Population(MatingPool));
        
         offspring_size = size(Offspring,2);
        
        offspring_disx = zeros(offspring_size,1);
        if Global.M == 3
            MatrixA = cross(repmat(best_point - worst_point,offspring_size,1),Offspring.objs - repmat(worst_point,offspring_size,1));
            for j = 1 : offspring_size
                offspring_disx(j) = norm(MatrixA(j,1:end))/norm(best_point - worst_point);
            end
        elseif Global.M == 2
            MatrixB = Offspring.objs;
            for j = 1 : offspring_size
                offspring_disx(j) = abs(det([best_point-worst_point;MatrixB(j,1:end)-worst_point]))/norm(best_point-worst_point);
            end
        else
            offspring_disx = pdist2(Offspring.objs,best_point - worst_point);
        end
        
        offspring_Gx = Offspring.objs - repmat(point,offspring_size,1);
        offspring_GX = abs(offspring_Gx);
        for k = 1 : offspring_size  % 在ROI内的惩罚为0
            if isempty(find(abs(offspring_Gx(k,1:end)) > r))
                offspring_GX(k,1:end) = 0;
            end
        end
        if ROI == true
            for k = 1 : offspring_size
                if isempty(find(offspring_Gx(k,1:end) >= 0))  % 对于ROI在可行域，支配参考点的区域惩罚为0
                    offspring_GX(k,1:end) = 0;
                end
            end
        end
        FixObj_offspring = (Offspring.objs-repmat(ideal_p,offspring_size,1))./(repmat(worst_p - ideal_p,offspring_size,1)) + proportion*(offspring_GX./repmat(GX_max,offspring_size,1)) + 1/proportion*offspring_disx./MaxDisx;
        %FixObj_offspring = (Offspring.objs-repmat(ideal_p,offspring_size,1))./(repmat(worst_p - ideal_p,offspring_size,1)) + proportion*(offspring_GX./repmat(GX_max,offspring_size,1));
        %FixObj_offspring = (Offspring.objs-repmat(ideal_p,offspring_size,1))./(repmat(worst_p - ideal_p,offspring_size,1)) + 5*(offspring_GX./repmat(GX_max,offspring_size,1)) + 1*offspring_disx./MaxDisx;
        
        [Population,FrontNo,CrowdDis] = pEnvironmentalSelection([Population,Offspring],Global.N,[FixObj;FixObj_offspring]);
        
      filename = strcat(func2str(Global.algorithm),func2str(Global.problem),'_',num2str(Global.gen));
        out_file=['0.7\',filename,'.txt'];
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