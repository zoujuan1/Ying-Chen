function preMOEAD(Global)
% <algorithm> <A>
% Multiobjective evolutionary algorithm based on decomposition
% type --- 2 --- The type of aggregation function
% point --- --- Preferred point
% r --- 0.2 ---- preferred region
% nr    ---   2 --- Maximum number of solutions replaced by each offspring

%------------------------------- Reference --------------------------------
% a unified approach for MOEA to solve preference-based MOPs
%--------------------------------------------------------------------------

%% Parameter setting
[type,point,r,nr] = Global.ParameterSet(2,zeros(1,Global.M)+0.5,0.2,2);

%% Generate the weight vectors
[W,Global.N] = UniformPoint(Global.N,Global.M);
T = ceil(Global.N/10);

%% Detect the neighbours of each solution
B = pdist2(W,W);
[~,B] = sort(B,2);
B = B(:,1:T);

%% Generate random population
Population = Global.Initialization();
%Z = min(Population.objs,[],1);

alpha = 1;  
beta = 5;
%% Optimization
while Global.NotTermination(Population)

    % 惩罚系数 proportion
    if Global.gen <= 0.6*Global.maxgen     % proportion  = a + b/N*(1+k/T)^cp, cp = log(N(b-a)/b)/log1.6
        proportion = alpha + beta/Global.N * (1 + Global.gen/Global.maxgen)^(log((beta-alpha)*Global.N/beta)/log(1.6));
    else                                % proportion = b
        proportion = beta;
    end
    % proportion = (1) 1/N*(1+k/T)^cp, cp = logN/log1.6;  (2) 1
        
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
        
    
    % For each solution
    for i = 1 : Global.N
        % Choose the parents
        P = B(i,randperm(size(B,2)));
    
        %Offspring = GAhalf(Population(P(1:2)));
       Offspring = Global.Variation(Population([i,P(1:2)]),1,@DE);   %子代个体采用DE生成
        
        offspring_gx = abs(Offspring.obj - point);
        
        if isempty(find(offspring_gx > r)) %子代在ROI内惩罚为0
            offspring_gx(1:end) = 0;
        end
        
        if isempty(find((Offspring.obj - point) > 0))  %子代支配参考点
            offspring_gx(1:end) = 0;
        end
        
        if Global.M == 3
           disx_offspring = norm(cross(best_point - worst_point, Offspring.obj - worst_point))/norm(best_point - worst_point);
        elseif Global.M == 2
            disx_offspring = abs(det([best_point-worst_point;Offspring.objs-worst_point]))/norm(best_point-worst_point);
        else
            disx_offspring = pdist2(Offspring.obj,best_point - worst_point);
        end

        FixObj_offspring = (Offspring.obj-ideal_p)./(worst_p - ideal_p) + proportion*(offspring_gx./GX_max) + 1/proportion*disx_offspring/MaxDisx;
        %FixObj_offspring = (1 + proportion)*((Offspring.obj-ideal_p)./(worst_p - ideal_p) + (offspring_gx./GX_max));
        %FixObj = (1 + proportion)*((Population.objs-repmat(ideal_p,Global.N,1))./(repmat(worst_p - ideal_p,Global.N,1)) + (GX./repmat(GX_max,Global.N,1)));
        
        % Update the ideal point
        Z = min([FixObj;FixObj_offspring],[],1);
        
        %% 聚合函数利用修正化后的obj
        switch type
            case 1
                % PBI approach
                normW   = sqrt(sum(W(P,:).^2,2));
                normP   = sqrt(sum((FixObj(P,1:end)-repmat(Z,T,1)).^2,2));
                normO   = sqrt(sum((FixObj_offspring-Z).^2,2));
                CosineP = sum((FixObj(P,1:end)-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
                CosineO = sum(repmat(FixObj_offspring-Z,T,1).*W(P,:),2)./normW./normO;
                g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
            case 2
                % Tchebycheff approach
                g_old = max(abs(FixObj(P,1:end)-repmat(Z,T,1)).*W(P,:),[],2);
                g_new = max(repmat(abs(FixObj_offspring-Z),T,1).*W(P,:),[],2);
            case 3
                % Tchebycheff approach with normalization
                Zmax  = max(FixObj(P,1:end),[],1);
                g_old = max(abs(FixObj(P,1:end)-repmat(Z,T,1))./repmat(Zmax-Z,T,1).*W(P,:),[],2);
                g_new = max(repmat(abs(FixObj_offspring-Z)./(Zmax-Z),T,1).*W(P,:),[],2);
            case 4
                % Modified Tchebycheff approach
                g_old = max(abs(FixObj(P,1:end)-repmat(Z,T,1))./W(P,:),[],2);
                g_new = max(repmat(abs(FixObj_offspring-Z),T,1)./W(P,:),[],2);
        end
        Population(P(g_old>=g_new)) = Offspring;
        %Population(P(find(g_old>=g_new,nr))) = Offspring;
    end
    
%      if Global.gen == Global.maxgen
%         filename = strcat(func2str(Global.algorithm),'_',func2str(Global.problem),'_');
%         out_file=[filename,'.txt'];
%         fid = fopen(out_file,'wt');
%         data = Population.objs;
%         [row,col] = size(data);
%         %输出种群
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
%      end
end
end