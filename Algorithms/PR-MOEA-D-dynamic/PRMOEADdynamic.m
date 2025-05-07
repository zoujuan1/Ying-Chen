function PRMOEADdynamic(Global)
% <algorithm> <H-N>
% NUMS-MOEA/D: integration of preference in MOEA/D
% propotion --- 0.2 --- the range of ROI
% kind --- 1 --- The type of aggregation function

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

%% Parameter setting
[propotion,kind] = Global.ParameterSet(0.2,1);
point = zeros(1, Global.M);
t =0.1* floor(Global.gen/Global.tao_n);
point(1) = (sin(pi*t)).^2 + 0.1;
for i = 2 : Global.M
    point(i) = (cos(pi*t)).^2 + 0.1;
end

%% Generate the weight vectors
XX = Global.N;
[W,Global.N] = PRPoint(Global.N,Global.M,point,propotion);
% [W,Global.N] = UniformPoint(Global.N,Global.M);  % IM-MOEAD
T = ceil(Global.N/10);

%% Detect the neighbours of each solution
B = pdist2(W,W);
[~,B] = sort(B,2);
B = B(:,1:T);

%% 保留初始种群
Initpop = [Global.Initialization() Global.Initialization() Global.Initialization()];

detect_num= randperm(Global.N,floor(Global.N/10));
Detect_pop = Initpop(:,detect_num);   %检测环境变化的部分个体

%% 搜索离每个向量最近的个体
C = pdist2(Initpop.objs,W);
[~,C] = sort(C,1);
A = C(1,1:end);
for i = 1 : Global.N
    for j = i+1 : Global.N
        if A(j) == A(i)
            A(j) = C(i+1,j);
        end
    end
end

%% Generate random population
Population = Initpop(:,A);
%   Population = Global.Initialization();
Z = min(Population.objs,[],1);

%% 随机选择一个个体作检测个体
DetePop =  Population(unidrnd(Global.N,1));


%% 部分参数或计数设置
PS_old = [];
k1 = [];
Tchebycheff = [];
Tchebycheff_k = 1;


%% Optimization
while Global.NotTermination(Population)
    DF("FDA2",3.5)
%     %% 输出每次变化的种群
%     if rem(Global.Current+1,Global.tao_n) == 0
%         filename = strcat(func2str(Global.algorithm),'_',func2str(Global.problem),'_t',num2str(floor(Global.Current()/Global.tao_n)),'_run',num2str(Global.run));
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
%     end
    
    
    
    %% Update the reference point
    t =0.1* floor(Global.Current()/Global.tao_n);
    point(1) = (sin(pi*t)).^2 + 0.1;
    for i = 2 : Global.M
        point(i) = (cos(pi*t)).^2 + 0.1;
    end
    
    
    %% 检测参考点变化情况
    point_change = false;
    avg_dis = mean(pdist2(point,W));
    if Global.gen >=20
        if avg_dis == dis_old
            point_change = false;
        else
            point_change = true;
        end
    end
    dis_old = avg_dis;
    
    %% 参考点发生变化
    if point_change == true
        %% Update the weight vectors
        [W,Global.N] = PRPoint(XX,Global.M,point,propotion);
        T = ceil(Global.N/10);
        %% Detect the neighbours of each solution
        B = pdist2(W,W);
        [~,B] = sort(B,2);
        B = B(:,1:T);
        
        %             %% 预测模型
        %             % 采用线性回归拟合预测PS
        %             if isempty(PS_old) %第一次预测只用上一次的
        %                 PS = Population.decs;
        %                 x = [ones(Global.N,1) PS(:,2:end) PS(:,2:end).^2];
        %                 y = PS(:,1);
        %                 k = regress(y,x);
        %                 PS_old = Population.decs;
        %                 error = std(PS_old);
        %
        %             else   %之后的预测用前两次的
        %                 PS = PS_old;
        %                 PS1 = Population.decs;
        %                 x = [ones(Global.N,1) PS(:,2:end) PS(:,2:end).^2];
        %                 x1 = [ones(Global.N,1) PS1(:,1:end-1) PS1(:,1:end-1).^2];
        %                 y = PS(:,1);
        %                 k = regress(y,x);
        %                 y1 = PS1(:,end);
        %                 k1 = regress(y1,x1);
        %                 error = std(PS_old);
        %
        %                 PS_old = Population.decs;
        %             end
        %                %% 在变化前更新初始种群
        %                for ii = 1 : Global.N
        %                    Initpop(A(ii)).dec = Population(ii).dec;
        %                    Initpop(A(ii)) = INDIVIDUAL(Initpop(A(ii)).dec);
        %                end
        
        %% 搜索离每个向量最近的个体
        C = pdist2(Initpop.objs,W);
        [~,C] = sort(C,1);
        A = C(1,1:end);
        for i = 1 : Global.N
            for j = i+1 : Global.N
                if A(j) == A(i)
                    A(j) = C(i+1,j);
                end
            end
        end
        Population = Initpop(:,A);
        
        %                 %% 预测更新PS
        %             for i = 1 : Global.N
        %                     ps = Population(i).dec;
        %                     ps(1) = [ones(1,1) ps(:,2:end) ps(:,2:end).^2]*k;
        %                     ps = ps + (diag(normrnd(0,0.01,Global.D,1)*error))';
        %                     ind = INDIVIDUAL(ps);
        %                     if ~isempty(k1)
        %                         ps(end) = [ones(1,1) ps(:,1:end-1) ps(:,1:end-1).^2]*k1;
        %                         ps = ps + (diag(normrnd(0,0.01,Global.D,1)*error))';
        %                         ind1 = INDIVIDUAL(ps);
        %                         if IsBetter(ind1.obj, ind.obj, W, Z)
        %                            ind = ind1;
        %                         end
        %                     end
        %                     if IsBetter(Population(i).obj,ind.obj, W, Z)  %如果预测前个体更优
        %                         Population(i) = Population(i);
        %                     else      % 否则替换个体为预测个体
        %                         Population(i) = ind;
        %                     end
        %             end
    end
    
    % For each solution
    for i = 1 : Global.N
        % Choose the parents
        P = B(i,randperm(size(B,2)));
        
        %% decide whether an inverse model is used to generate an offspring solution
        fj = Population(P).objs;
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
            F = (Population(P).objs)';
            BB = (X*(F'))*inv(F*F');
            diff = diag(Population(i).obj - Z + 0.01.*ones(1,Global.M));
            Fdesired = Population(i).obj*diff*diag(normrnd(0,1,1,Global.M));
            while ~IsDominated(Fdesired, Population(i).obj, Global.M)
                Fdesired = Population(i).obj*diff*diag(normrnd(0,1,1,Global.M));
            end
            x_offspring = (BB*(Fdesired'))';
            Offspring =  INDIVIDUAL(x_offspring);
        else
            Offspring = Global.Variation(Population(P(1:2)),1);
        end
        
        %              % Generate an offspring
        %             Offspring = Global.Variation(Population(P(1:2)),1);
        
        % Update the ideal point
        Z = min(Z,Offspring.obj);
        
        % Update the neighbours
        switch kind
            case 1
                % PBI approach
                normW   = sqrt(sum(W(P,:).^2,2));
                normP   = sqrt(sum((Population(P).objs-repmat(Z,T,1)).^2,2));
                normO   = sqrt(sum((Offspring.obj-Z).^2,2));
                CosineP = sum((Population(P).objs-repmat(Z,T,1)).*W(P,:),2)./normW./normP;
                CosineO = sum(repmat(Offspring.obj-Z,T,1).*W(P,:),2)./normW./normO;
                g_old   = normP.*CosineP + 5*normP.*sqrt(1-CosineP.^2);
                g_new   = normO.*CosineO + 5*normO.*sqrt(1-CosineO.^2);
            case 2
                % Tchebycheff approach
                g_old = max(abs(Population(P).objs-repmat(Z,T,1)).*W(P,:),[],2);
                g_new = max(repmat(abs(Offspring.obj-Z),T,1).*W(P,:),[],2);
            case 3
                % Tchebycheff approach with normalization
                Zmax  = max(Population.objs,[],1);
                g_old = max(abs(Population(P).objs-repmat(Z,T,1))./repmat(Zmax-Z,T,1).*W(P,:),[],2);
                g_new = max(repmat(abs(Offspring.obj-Z)./(Zmax-Z),T,1).*W(P,:),[],2);
            case 4
                % Modified Tchebycheff approach
                g_old = max(abs(Population(P).objs-repmat(Z,T,1))./W(P,:),[],2);
                g_new = max(repmat(abs(Offspring.obj-Z),T,1)./W(P,:),[],2);
        end
        
        %             Dis = pdist2(Population(i).obj,W);
        %             [~,Dis] = sort(Dis,2);
        %             w = W(Dis(1),:);
        %             Tchebycheff(Tchebycheff_k,1) = abs(max((Population(i).obj-Z).*w) - max((Offspring.obj-Z).*w));
        %             Tchebycheff_k = Tchebycheff_k + 1;
        Population(P(g_old>=g_new)) = Offspring;
    end
    
    %         detect = false;
    %         detect_second = false;
    %         %Z检测
    %         if Tchebycheff_k> Global.N
    %           u_global = mean(Tchebycheff(1:Tchebycheff_k-Global.N,1),1);
    %           delta_global = std(Tchebycheff(1:Tchebycheff_k-Global.N,1),0,1);
    %           u_windows = mean(Tchebycheff(Tchebycheff_k-Global.N+1:end,1),1);
    %           delta_windows = std(Tchebycheff(Tchebycheff_k-Global.N+1:end,1),0,1);
    %           SE = (delta_global - delta_windows)./(Global.N^0.5);
    %           zz = (u_windows - u_global)./SE;
    %           if zz > 1
    %               detect = true;
    %           end
    %         end
    %
    %         %二级检测
    %         if detect == true
    %             Newpop = INDIVIDUAL(DetePop.dec);
    %             if Newpop.obj == DetePop.obj
    %                 detect_second = false;
    %             else
    %                 detect_second = true;
    %                 DetePop = Newpop;
    %                 Tchebycheff_k = 1;
    %             end
    %         end
    %         if detect_second == true
    %            for j = 1 : Global.N
    %                 Population(j) = INDIVIDUAL(Population(j).dec);
    %            end
    %            Z = min(Population.objs,[],1);  % reset Z
    %         end
    
    
            if rem(Global.Current,Global.tao_n) == 0
               for j = 1 : Global.N
                    Population(j) = INDIVIDUAL(Population(j).dec);
               end
               Z = min(Population.objs,[],1);  % reset Z
            end
    
end
end