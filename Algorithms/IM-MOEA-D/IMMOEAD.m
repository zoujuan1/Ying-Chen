function IMMOEAD(Global)
% <algorithm> <H-N>
% IM-MOEA/D: A Multiobjective Evolutionary Algorithm Based on Decomposition
% for dynamic MOPs
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
    kind = Global.ParameterSet(1);

    %% Generate the weight vectors
    [W,Global.N] = UniformPoint(Global.N,Global.M);
    
    T = ceil(Global.N/10);

    %% Detect the neighbours of each solution
    B = pdist2(W,W);
    [~,B] = sort(B,2);
    B = B(:,1:T);
    
    %% Generate random population
    Population = Global.Initialization();
    Z = min(Population.objs,[],1);
    
    %% 随机选择一个个体作检测个体
     DetePop =  Population(unidrnd(Global.N,1));

    %% some detection segments
    eval_trigger = 0;
    e = [];
    e_k = 1; %current number of evaluation
    Tchebycheff = [];
    Tchebycheff_k = 1
    PP=cell(1,1000);
    ind=1;
    %% Optimization
    while Global.NotTermination(Population)
        P = PPF(10000,Global);
        if mod(Global.gen,50)
            PP(ind)=P;
            ind=ind+1;
        end
        Draw( PP(ind),'k+','Markeredgecolor',[.1 .4 .6],'Markerfacecolor',[.9 .8 .9]);
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
                
                e(e_k,:) = Fdesired - Offspring.obj;
                e_k = e_k + 1;
            else
               Offspring = Global.Variation(Population(P(1:2)),1);
            end

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
            
             Dis = pdist2(Population(i).obj,W);
            [~,Dis] = sort(Dis,2);
            w = W(Dis(1),:);
            Tchebycheff(Tchebycheff_k,1) = max((Population(i).obj-Z).*w) - max((Offspring.obj-Z).*w);
            Tchebycheff_k = Tchebycheff_k + 1;
            
            Population(P(g_old>=g_new)) = Offspring;        
            
            
%             %% change detection
%             % first-stage
%             detect_first = false;
%             u_global = mean(e,1);
%             delta_global = std(e,0,1);
%             if e_k >= 20
%                 u_windows = mean(e(e_k-19:end,1:end),1);
%                 delta_windows = std(e(e_k-19:end,1:end),0,1);
%                 SE = (delta_global-delta_windows)./10;
%                 zz = (u_windows - u_global)./SE;
%                 if ~isempty(find(zz > 1))
%                     detect_first = true;
%                 end
%             end
%             
%             %second-stage
%             detect_second = false;
%             if e_k - eval_trigger >= 5 && detect_first == true
%                 detect_pop = randperm(Global.N,1);
%                 ps = Population(detect_pop).dec;
%                 detect_ind = INDIVIDUAL(ps);
%                 if detect_ind.obj ~= Population(detect_pop).obj
%                     detect_second = true;
%                 end
%                 eval_trigger = e_k;
%             end
        
%              if rem(Global.Current,Global.tao_n) == 0
%                 detect_second = true;
%              end
%             
%             %% 检测到环境变化
%             if detect_second == true
%                 for j = 1 : Global.N
%                     Population(j) = INDIVIDUAL(Population(j).dec);
%                 end
%                 Z = min(Population.objs,[],1);  % reset Z
%                 break;
%             end  
        end   
        
        % Tchebycheff变化幅度判断环境变化
        detect = false;  
        detect_second = false;
        %Z检测
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
        
        %威尔逊检测
        if Tchebycheff_k> Global.N
            [p,h] = ranksum(Tchebycheff(1:Tchebycheff_k-Global.N,1),Tchebycheff(Tchebycheff_k-Global.N+1:end,1),0.05);
            if h == 1
                detect = true;
            end
        end
        
        %二级检测
        if detect == true
            Newpop = INDIVIDUAL(DetePop.dec);
            if Newpop.obj == DetePop.obj
                detect_second = false;
            else
                detect_second = true;
                DetePop = Newpop;
                Tchebycheff_k = 1;
            end
        end
        
        if detect_second == true
           for j = 1 : Global.N
                Population(j) = INDIVIDUAL(Population(j).dec);
           end
           Z = min(Population.objs,[],1);  % reset Z
        end    
        
        
        
        %% 输出每次环境变化种群
%           if rem(Global.Current+1,Global.tao_n) == 0
%                 filename = strcat(func2str(Global.algorithm),'_',func2str(Global.problem),'_gen',num2str(Global.Current),'_run',num2str(Global.run));
%                 out_file=[filename,'.txt'];
%                 fid = fopen(out_file,'wt');
%                 data = Population.objs;
%                 [row,col] = size(data);
%                 for i = 1 : row
%                     for j = 1 : col
%                         if j == col
%                             fprintf(fid,'%g\n',data(i,j));
%                         else
%                             fprintf(fid,'%g\t',data(i,j));
%                         end
%                     end
%                 end
%                 fprintf(fid,'\n\n\n');
%                 fclose(fid);  
%           end
    end
end

function P = PPF(input,Global)
            f(:,1)    = (0:1/(input-1):1)';
            f(:,2)    = 1-f(:,1).^(0.2 * floor(Global.Current()/Global.tao_n));
            P = f;
end