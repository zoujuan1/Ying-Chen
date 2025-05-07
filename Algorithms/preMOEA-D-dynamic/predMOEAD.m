function predMOEAD(Global)
% <algorithm> <A>
% Multiobjective evolutionary algorithm based on decomposition
% type --- 2 --- The type of aggregation function
% point --- --- Preferred point
% r --- 0.1 ---- preferred region

%% Parameter setting
[type,point,r] = Global.ParameterSet(2,[0.2,0.8],0.1);

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

%% some detection segments
e = [];
e_k = 1; %current number of evaluation

alpha = 1;  
beta = 5;

%环境变化
changed = false;


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
    %偏好信息变
        switch rem(ceil(Global.gen/50),7)
            case 0
                point = [0.2,0.8];
            case 1
                point = [0.6,0.4];
            case 2
                point = [0.3,0.7];
            case 3
                point = [0.7,0.3];
            case 4
                point = [0.4,0.6];
            case 5
                point = [0.8,0.2];
            case 6
                point = [0.5,0.5];
        end
    
    % 惩罚系数 proportion
    if Global.gen <= 0.6*Global.maxgen     % proportion  = a + b/N*(1+k/T)^cp, cp = log(N(b-a)/b)/log1.6
        proportion = alpha + beta/Global.N * (1 + Global.gen/Global.maxgen)^(log((beta-alpha)*Global.N/beta)/log(1.6));
    else                                % proportion = b
        proportion = beta;
    end
    % proportion = (1) 1/N*(1+k/T)^cp, cp = logN/log1.6;  (2) 1
        
    % For each solution
    for i = 1 : Global.N
        % Choose the parents
        P = B(i,randperm(size(B,2)));
        
        %% 距离惩罚
        best_point = point - r;
        worst_point = point + r;
        Disx = pdist2(Population.objs,best_point - worst_point);
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
            if isempty(find(Gx(k,1:end) > r))
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
        
        % 修正种群适应度
        FixObj = (Population.objs-repmat(ideal_p,Global.N,1))./(repmat(worst_p - ideal_p,Global.N,1)) + 5*(GX./repmat(GX_max,Global.N,1)) + 0.5*Disx./MaxDisx;
        
        %环境变化判断
        if Global.gen > 10  %前10代不做判断
            if isequal(FixObj(DetePop), DeteObj) ~= 1
                changed = true;
            end
        end
        
         % decide whether an inverse model is used to generate an offspring solution
        fj = FixObj(P,1:end);
        wj = zeros(T,Global.M);
        rj = zeros(1,Global.M);
        k = W(B(i,:));
        for j = 1 : Global.M
            wj(:,j) = k(:,j);
        end
        kx = fj-repmat(mean(fj,1),T,1);
        ky = wj-repmat(mean(wj,1),T,1);
        kup = sum(kx.*ky,1);
        kxx = kx.^2;
        kyy = ky.^2;
        kdown = sum(kxx.*kyy,1);
        rj(1,1:end) = kup(1,1:end)./(kdown(1,1:end).^0.5);
        
        % Generate an offspring
        use_inverse = ~isempty(find(rj>0.7)) || ~isempty(find(rj<-0.7));
        if use_inverse == true  % use inverse model
            X = (Population(P).decs)';
            F = (FixObj(P,1:end))';
            BB = (X*(F'))*inv(F*F');
            diff = diag(FixObj(i,1:end) - Z + 0.01.*ones(1,Global.M));
            Fdesired = FixObj(i,1:end)*diff*diag(normrnd(0,1,1,Global.M));
            while ~IsDominated(Fdesired, FixObj(i,1:end), Global.M)
                Fdesired = FixObj(i,1:end)*diff*diag(normrnd(0,1,1,Global.M));
            end
            x_offspring = (BB*(Fdesired'))';
            Offspring =  INDIVIDUAL.Individual(x_offspring);
            
            e(e_k,:) = Fdesired - Offspring.obj;
            e_k = e_k + 1;
        else
            Offspring = Global.Variation(Population([i,P(1:2)]),1,@DE);
        end
        
        %Offspring = GAhalf(Population(P(1:2)));
        %Offspring = DE(Population(i),Population(P(1)),Population(P(2)));   %子代个体采用DE生成

        offspring_gx = abs(Offspring.obj - point);
        
        if isempty(find(offspring_gx > r)) %子代在ROI内惩罚为0
            offspring_gx(1:end) = 0;
        end
        
        if isempty(find((Offspring.obj - point) > 0))  %子代支配参考点
            offspring_gx(1:end) = 0;
        end

        disx_offspring = pdist2(Offspring.obj,best_point - worst_point);
    
        FixObj_offspring = (Offspring.obj-ideal_p)./(worst_p - ideal_p) + 5*(offspring_gx./GX_max ) + 0.5*disx_offspring/MaxDisx; % 修正子代适应度
        
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
    end
    
    %检测机制   
    Disx = pdist2(Population.objs,best_point - worst_point);
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
        if isempty(find(Gx(k,1:end) > r))
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
    ideal_p = min(Population.objs,[],1);
    worst_p = max(Population.objs,[],1);
    FixObj = (Population.objs-repmat(ideal_p,Global.N,1))./(repmat(worst_p - ideal_p,Global.N,1)) + 5*(GX./repmat(GX_max,Global.N,1)) + 0.5*Disx./MaxDisx;
    
    DetePop = randperm(Global.N,10);
    DeteObj = FixObj(DetePop);
    
%     if changed == true  %确认环境变化
%         for j = 1 : Global.N
%             Population(j) = INDIVIDUAL.Individual(Population(j).dec);
%         end
%     end

    if (rem(Global.gen + 1, 50) == 0) || (rem(Global.gen,50)==0)  %确认环境变化
        for j = 1 : Global.N
            Population(j) = INDIVIDUAL.Individual(Population(j).dec);
        end
    end
    
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