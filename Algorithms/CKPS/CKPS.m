function CKPS(Global)
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
    metric=[];
    
    N1=9     %Knee 点个数
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
    KneesPoints=zeros(N1,Global.D,200);  %初始化历史knee点集合

    PF=Global.PF;
    while Global.NotTermination(Population)
%         Cc=INDIVIDUAL(sum(Population.decs,1)/Global.N);
%         DF(func2str(Global.problem),Global.NT);    %画出对应PF
%         Draw(Global.PF,'ok','Markeredgecolor',[.9 .0 .0],'Markerfacecolor',[.9 .0 .0],'MarkerSize',2);
%         Draw(Cc.objs,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',5);
%         Draw(Pt1.objs,'k+','Markeredgecolor',[.6 .0 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',3);

       
       
        if Change(Global)   % 检测环境变化
            t=t+1      
            Decs=Population.decs;
            Population=INDIVIDUAL(Decs);
            Pt2=Population;

            [FrontNo,MaxFNo] = NDSort(Pt2.objs,Global.N);
            First = find(FrontNo==1);   
            Ct=sum(Pt2(First).decs,1)/length(First);  % 计算非支配层中心点
            CenterPoint=[CenterPoint;Ct];       % 将Ct加进历史中心点

            [P1,KneesPoints]=KneePoint(Pt2(First),Global,t,N1,KneesPoints);               %预测Knee点
            P2=prediction(Pt2(First),Global,t,CenterPoint);              %预测中心点
            Population=RandGenrate(P1,P2,Global);                       %随机产生个体
      
            PF=Global.PF;
            Pt1=Population;
            if t>=2      %除去第一次环境变化
                 score=IGD(Population.objs,PF);      %计算IGD
                 metric=[metric,score];
                 migd=sum(metric)/length(metric); %计算MIGD
                 igd=migd;                         %将MIGD通过全局变量igd传给指标函数MIGD()
            end
        else
%             MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
%             Offspring  = Global.Variation(Population(MatingPool));
              Offspring  = Global.Variation(Population,Global.N,@RMMEDA_operator);
           [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
        end
    end
end

function [P1,KneesPoints]=KneePoint(Pop,Global,t,N1,KneesPoints)    %预测Knee点
    L=length(Pop);                         %得到非支配种群大小
    Pop_Obj=Pop.objs;
    [~,index]=sortrows(Pop_Obj,1);        %对非支配种群 根据第一个目标排序 
    N2=floor(L/N1);                               %将种群分成N1份，每份N2个 个体
    P1=[];
    if N2~=0
        KneesD=zeros(N1,Global.D);
        Knees=[];                             %初始化当前Knee点集合
        for i = 1:N1                        %对每份个体求Knee点
            if i==N1
                Pop2_Obj=Pop_Obj(index((i-1)*N2+1:end),:);   
            else
                Pop2_Obj=Pop_Obj(index((i-1)*N2+1:i*N2),:);  
            end
            if Global.M==2                    %如果是2维目标，确定2个边界点，计算个体点到直线的距离
                [~,index1]=max(Pop_Obj(:,1)); %第一个目标的最大索引
                [~,index2]=max(Pop_Obj(:,2)); %第二个目标的最大索引
                if  Pop_Obj(index1,1)~=Pop_Obj(index2,1) &&  Pop_Obj(index1,2)~=Pop_Obj(index2,2)         %边界点每维都不能相等
                    Point1=Pop_Obj(index1,:); %得到边界点 Point1
                    Point2=Pop_Obj(index2,:); %得到边界点 Point2

                    syms K B
                    eq1=K*Point1(1)-Point1(2)+B==0;   %定义直线方程
                    eq2=K*Point2(1)-Point2(2)+B==0;   %定义直线方程
                    [K,B]=solve(eq1,eq2,K,B);         %解方程

                    double(K);
                    double(B);
                    dis=abs(double(K)* Pop2_Obj(:,1) -  Pop2_Obj(:,2) +double(B))/(sqrt(double(K)^2+1));
                    [~,index3]=min(dis);
                    knee=Pop(index((i-1)*N2+index3)); 
                    Knees=[Knees,knee];                  %把当前区域的knee点 加进knee点集
                end
            elseif Global.M==3               %如果是3维目标，确定3个点，计算个体点到超平面的距离

            end
        end
        if ~isempty(Knees)&&length(Knees)==N1
            KneesD=Knees.decs;
            KneesPoints(:,:,t)=KneesD;         %将当前knee点集合加进历史knee点集合
        end
    
        if t<=6||length(Knees)~=N1
            P1=[];
        else                                   %t大于6时，用AR模型预测下一环境的 knee点集合
            for i = 1:N1 
                 Knees_i=[];                     %第i个区域 的历史knee点集合
                 for j=1:t
                     temp=KneesPoints(i,:,j);
                     Knees_i=[Knees_i;temp];
                 end
          
                 NextKnee_i=Center_P(Knees_i,t,Global);  %预测第i区域 下一环境的knee点
                 P=INDIVIDUAL(NextKnee_i);
                 P1=[P1,P];
            end
        end
    end
end

function P2=prediction(Pt2,Global,t,CenterPoint)    %预测中心点
     Ob=Pt2;
     OP=length(Ob);
     if t>=2
        Decs=Ob.decs+(CenterPoint(t,:)-CenterPoint(t-1,:));               %预测
        for i = 1 : OP              %修正decs
            for j =1 : Global.D
                  if  Decs(i,j)<Global.lower(j)
                          Decs(i,j)=0.5*(CenterPoint(t,j)+Global.lower(j));
                  end
                  if   Decs(i,j)>Global.upper(j)
                          Decs(i,j)=0.5*(CenterPoint(t,j)+Global.upper(j));
                  end  
            end
        end
        P2=INDIVIDUAL(Decs);
     else
        P2=INDIVIDUAL(Ob.decs);
     end
end

function P3=RandGenrate(P1,P2,Global)              %自适应随机产生
    if length(P1)+length(P2)>=Global.N
        [P3,~,~] = EnvironmentalSelection([P1,P2],Global.N);
    else
        P3=[P1,P2];
        for i =1:(Global.N-length(P1)-length(P2))
            Decs=(Global.upper-Global.lower).*rand(1,Global.D);
            P=INDIVIDUAL(Decs);
            P3=[P3,P];
        end
    end
end

function NC=Center_P(CenterPoint,t,Global)
    %% 根据历史点建立AR模型，并预测中心点,这里的CenterPoint是knee点
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
