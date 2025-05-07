function PPSD(Global)
% <algorithm> <A>
% Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    global igd;
    igd=0.1;
    metric1=[];
    T=[];
    global hvd;
    hvd = 0.1;
    metric2 =[];
    %% Generate random population
%     Global.NT=10;
    Population = Global.Initialization();
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);
    %% Optimization
    
    %% Set time step t
    t=0;             %时间序列
    p=3;
    M=23;
    Pt1=Population;  % 上一环境的种群
    Pt2=Population;  % 当前环境的种群
    CenterPoint=[];  % 历史中心点

    PF=Global.PF;
    while Global.NotTermination(Population)
%         Cc=INDIVIDUAL(sum(Population.decs,1)/Global.N);
%         DF(func2str(Global.problem),Global.NT);    %画出对应PF
%         Draw(Global.PF,'ok','Markeredgecolor',[.9 .0 .0],'Markerfacecolor',[.9 .0 .0],'MarkerSize',2);
%         Draw(Cc.objs,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',5);
%         Draw(Pt1.objs,'k+','Markeredgecolor',[.6 .0 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',3);

       
        if Global.evaluated>0.9*Global.evaluation     
             y=0;
        end
        if Change(Global)   % 检测环境变化
            t=t+1           
            Pt2=Population;

%             newPop=[];                  % 预测的种群 
            Ct=sum(Pt2.decs,1)/Global.N;  % 计算中心点
            CenterPoint=[CenterPoint;Ct];       % 将Ct加进历史中心点
            St2=Pt2.decs-sum(Pt2.decs,1)/Global.N;  % 计算当前种群形状
            St1=Pt1.decs-sum(Pt1.decs,1)/Global.N;  % 计算上一环境种群形状
            if t<=6
                 Pp = Global.Initialization();
                 Population=[Pp(1:Global.N/2),Pt2(Global.N/2:Global.N)];
                 Decs=Population.decs;
                 Population=INDIVIDUAL(Decs);
            else

                 
                 S=Manifold_P(St1,St2,Global);    %根据St1和St2两个形状 估计下一代的形状
                 NC=Center_P(CenterPoint,Ct,t,Global,M);  %根据历史中心点 估计下一中心点
                 Decs = S + NC;                    % 得到预测的PS
                 
                 for i = 1 : Global.N              %修正decs
                     for j =1 : Global.D
                               if  Decs(i,j)<Global.lower(j)
                                   Decs(i,j)=0.5*(CenterPoint(t,j)+Global.lower(j));
                               end
                               if  Decs(i,j)>Global.upper(j)
                                   Decs(i,j)=0.5*(CenterPoint(t,j)+Global.upper(j));
                               end  
                     end
                 end
                
                 Population=INDIVIDUAL(Decs);      %根据预测的PS 生成个体，并替换到Population               
            end
            PF=Global.PF;
            Pt1=Population;
            if t>=2      %除去第一次环境变化
                 score1=IGD(Population.objs,Global.PF);      %计算IGD
                 metric1=[metric1,score1];
                 migd=sum(metric1)/length(metric1); %计算MIGD
                 igd=migd;                         %将MIGD通过全局变量igd传给指标函数MIGD()
                 T=[T,t];
                 score2=HV(Population.objs,Global.PF);      %计算HVD
                 metric2=[metric2,score2];
                 mhvd=sum(metric2)/length(metric2); %计算MHVD
                 hvd=mhvd;                         %将MHVD通过全局变量hvd传给指标函数MHVD()
            end
        else
            MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
            Offspring  = Global.Variation(Population(MatingPool));
           [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
        end
    end
end





function S=Manifold_P(St1,St2,Global)
      Distance=pdist2(St1,St2,'euclidean');
      sigma =((sum(min(Distance))/length(St1))^2)/Global.D;     % 正态分布的标准差
      epsilon = normrnd(0, sqrt(sigma),Global.N,Global.D);         
      S=St2+epsilon;
end
function NC=Center_P(CenterPoint,Ct,t,Global,M)
    %% 根据历史中心点建立AR模型，并预测中心点
   ind=t;
   fai=[];
   if t>=23
       ind=23;
   end
   NC=[];
   for i = 1 : Global.D
       fai=[];
       Psai=[];
       for j = 0 : ind-4
            r=[CenterPoint(t-1-j,i),CenterPoint(t-2-j,i),CenterPoint(t-3-j,i)];
            fai=[fai;r];
            Psai=[Psai,CenterPoint(t-j,i)];
       end
       Psai=Psai.';
%        fai0=fai.';
%        fai1=inv(fai0*fai);
%        lambda2=fai1*fai0*Psai
       lambda = fai\Psai;
      
%        lambda = fai\Psai;
       sigma=0;             %初始化误差标准差
       for j = 0 : ind-4
           sigma=sigma+(CenterPoint(t-j,i)-(lambda(1)*CenterPoint(t-j-1,i)+lambda(2)*CenterPoint(t-2-j,i)+lambda(3)*CenterPoint(t-3-j,i))).^2;
       end
       sigma=sigma/(ind-3); %计算误差方差
       epsilon = normrnd(0, sqrt(sigma),1,1);
       xi=(lambda(1)*CenterPoint(t,i))+(lambda(2)*CenterPoint(t-1,i))+(lambda(3)*CenterPoint(t-2,i))+epsilon;
       if  xi<Global.lower(i)
           xi=0.5*(CenterPoint(t-1,i)+Global.lower(i));
       end
       if  xi>Global.upper(i)
           xi=0.5*(CenterPoint(t-1,i)+Global.upper(i));
       end   
       NC=[NC,xi];
   end
   
   
end
function isture = Change(Global)
    isture= rem(Global.Current,Global.tao_n) == 0;
end
