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
    t=0;             %ʱ������
    p=3;
    M=23;
    Pt1=Population;  % ��һ��������Ⱥ
    Pt2=Population;  % ��ǰ��������Ⱥ
    CenterPoint=[];  % ��ʷ���ĵ�

    PF=Global.PF;
    while Global.NotTermination(Population)
%         Cc=INDIVIDUAL(sum(Population.decs,1)/Global.N);
%         DF(func2str(Global.problem),Global.NT);    %������ӦPF
%         Draw(Global.PF,'ok','Markeredgecolor',[.9 .0 .0],'Markerfacecolor',[.9 .0 .0],'MarkerSize',2);
%         Draw(Cc.objs,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',5);
%         Draw(Pt1.objs,'k+','Markeredgecolor',[.6 .0 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',3);

       
        if Global.evaluated>0.9*Global.evaluation     
             y=0;
        end
        if Change(Global)   % ��⻷���仯
            t=t+1           
            Pt2=Population;

%             newPop=[];                  % Ԥ�����Ⱥ 
            Ct=sum(Pt2.decs,1)/Global.N;  % �������ĵ�
            CenterPoint=[CenterPoint;Ct];       % ��Ct�ӽ���ʷ���ĵ�
            St2=Pt2.decs-sum(Pt2.decs,1)/Global.N;  % ���㵱ǰ��Ⱥ��״
            St1=Pt1.decs-sum(Pt1.decs,1)/Global.N;  % ������һ������Ⱥ��״
            if t<=6
                 Pp = Global.Initialization();
                 Population=[Pp(1:Global.N/2),Pt2(Global.N/2:Global.N)];
                 Decs=Population.decs;
                 Population=INDIVIDUAL(Decs);
            else

                 
                 S=Manifold_P(St1,St2,Global);    %����St1��St2������״ ������һ������״
                 NC=Center_P(CenterPoint,Ct,t,Global,M);  %������ʷ���ĵ� ������һ���ĵ�
                 Decs = S + NC;                    % �õ�Ԥ���PS
                 
                 for i = 1 : Global.N              %����decs
                     for j =1 : Global.D
                               if  Decs(i,j)<Global.lower(j)
                                   Decs(i,j)=0.5*(CenterPoint(t,j)+Global.lower(j));
                               end
                               if  Decs(i,j)>Global.upper(j)
                                   Decs(i,j)=0.5*(CenterPoint(t,j)+Global.upper(j));
                               end  
                     end
                 end
                
                 Population=INDIVIDUAL(Decs);      %����Ԥ���PS ���ɸ��壬���滻��Population               
            end
            PF=Global.PF;
            Pt1=Population;
            if t>=2      %��ȥ��һ�λ����仯
                 score1=IGD(Population.objs,Global.PF);      %����IGD
                 metric1=[metric1,score1];
                 migd=sum(metric1)/length(metric1); %����MIGD
                 igd=migd;                         %��MIGDͨ��ȫ�ֱ���igd����ָ�꺯��MIGD()
                 T=[T,t];
                 score2=HV(Population.objs,Global.PF);      %����HVD
                 metric2=[metric2,score2];
                 mhvd=sum(metric2)/length(metric2); %����MHVD
                 hvd=mhvd;                         %��MHVDͨ��ȫ�ֱ���hvd����ָ�꺯��MHVD()
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
      sigma =((sum(min(Distance))/length(St1))^2)/Global.D;     % ��̬�ֲ��ı�׼��
      epsilon = normrnd(0, sqrt(sigma),Global.N,Global.D);         
      S=St2+epsilon;
end
function NC=Center_P(CenterPoint,Ct,t,Global,M)
    %% ������ʷ���ĵ㽨��ARģ�ͣ���Ԥ�����ĵ�
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
       sigma=0;             %��ʼ������׼��
       for j = 0 : ind-4
           sigma=sigma+(CenterPoint(t-j,i)-(lambda(1)*CenterPoint(t-j-1,i)+lambda(2)*CenterPoint(t-2-j,i)+lambda(3)*CenterPoint(t-3-j,i))).^2;
       end
       sigma=sigma/(ind-3); %��������
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
