function PBDMO(Global)
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
    Score1=[];
    eva=[];
    PopObj=[];
    PopPF=[];
    
    global hvd;
    hvd = 0.1;
    metric2 =[];
    T=[];
    %% Generate random population
%     Global.NT=10;
    Population = Global.Initialization();
    [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);
    %% Optimization
    
    %% Set time step t
    t=0;             %ʱ������
    N2=Global.D;            %����Ҫά�Ĳ�������
    N1=Global.N/(2*N2);     %��Ҫά�Ĳ�����
    CenterPoint=[];  % ��ʷ���ĵ�
   
    Pt1=Population;  % ��һ��������Ⱥ
    Pt2=Population;  % ��ǰ��������Ⱥ
    while Global.NotTermination(Population)
%         Cc=INDIVIDUAL(sum(Population.decs,1)/Global.N);
%         DF(func2str(Global.problem),Global.NT);    %������ӦPF
%         Draw(Global.PF,'ok','Markeredgecolor',[.9 .0 .0],'Markerfacecolor',[.9 .0 .0],'MarkerSize',2);
%         Draw(Cc.objs,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',5);
%         Draw(Pt1.objs,'k+','Markeredgecolor',[.6 .0 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',3);
        
        
        if Global.evaluated>0.9*Global.evaluation
             y=0;
        end
        if Change(Global)  % ��⻷���仯
            k=0;
            t=t+1 
            Pt2=Population;
            
            Decs=Population.decs;
            Population=INDIVIDUAL(Decs);
            
            Ct= sum(Population.decs,1)/Global.N;  % ������Ⱥ���ĵ�
            CenterPoint=[CenterPoint;Ct];         % ��Ct�ӽ���ʷ���ĵ�
            
            P1=predict(Pt2,Pt1,Global);                       %Ԥ������
            P2=Sampling(Pt2,Global,N1,N2);                           %��������
            P3=Shrink(Pt2,Pt1,Global,N1,N2);                        %��������
            Pachive=[P1,P2,P3];     
            [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Pachive,Population],Global.N);
            Pt1=Population;
            if t>=2      %��ȥ��һ�λ����仯
                 score1=IGD(Population.objs,Global.PF);      %����IGD
                 metric1=[metric1,score1];
                 migd=sum(metric1)/length(metric1); %����MIGD
                 igd=migd;                         %��MIGDͨ��ȫ�ֱ���igd����ָ�꺯��MIGD()
                 T=[T,t];
                 Score1=[Score1,log10(score1)];
                 score2=HV(Population.objs,Global.PF);      %����HV
                 metric2=[metric2,score2];
                 mhvd=sum(metric2)/length(metric2); %����MGD
                 hvd=mhvd;                          %��Mhvdͨ��ȫ�ֱ���hvd����ָ�꺯��HV() 
                  temp=Population.objs+0.2*(t-1);
                  PopObj=[PopObj;temp];
                  PopPF=[PopPF;Global.PF+0.2*(t-1)];
            end
        else
           % MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
           % Offspring  = Global.Variation(Population(MatingPool));
            Offspring  = Global.Variation(Population,Global.N,@RMMEDA_operator);
            [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
          
           
        end
    end
end
function  P1 = predict(Pt1,Pt2,Global)       %��̽����
   Ct1= sum(Pt1.decs,1)/Global.N;
   Ct2= sum(Pt2.decs,1)/Global.N;        %���ĵ�
   Dir=Ct2-Ct1;                          %���ĵ�λ������
   [FrontNo,MaxFNo] = NDSort(Pt2.objs,Global.N);
   First = find(FrontNo==1);             %�õ�Pt2�ķ�֧���
   
   Pos1=Pt2(First).decs+0.5*Dir;              %Ԥ������Pos
   Pos2=Pt2(First).decs+Dir;
   Pos3=Pt2(First).decs+1.5*Dir;
   
   Pop1=INDIVIDUAL(Pos1);                %����3������Ⱥ
   Pop2=INDIVIDUAL(Pos2);
   Pop3=INDIVIDUAL(Pos3);
   
   P1=[Pop1,Pop2,Pop3];
end
function P2 = Sampling(Pt2,Global,N1,N2)                                   %��������
   [FrontNo,MaxFNo] = NDSort(Pt2.objs,Global.N);
   First = find(FrontNo==1);             %�õ�Pt2�ķ�֧���
   Decs=Pt2(First).decs;                        %�õ���֧����POS
   [m,index]=max(std(Decs,0,1));                       %��ÿһά���� �ķ�����ط������ֵm����������index
   mainD=index;                               %�������Ϊ ��ά��
   P2=[];
   
   for i = 1 : N1
        for j = 1:N2
            Decs=[];
            for k =1:Global.D
                if k==mainD                              %�����ά������ά��
                    temp=Global.lower(k)+((Global.upper(k)-Global.lower(k))*i)/N1;
                else
                    temp=Global.lower(k)+((Global.upper(k)-Global.lower(k))*j)/N2;
                end
            end
            Decs=[Decs,temp];
            p=INDIVIDUAL(Decs);
            P2=[P2,p];
        end
   end 
  [FrontNo,MaxFNo] = NDSort(P2.objs,Global.N);
  First = find(FrontNo==1); 
  P2=P2(First);
end
function P3 = Shrink(Pt2,Pt1,Global,N1,N2)                        %��������
   [FrontNo,MaxFNo] = NDSort(Pt2.objs,Global.N);
   First = find(FrontNo==1);             %�õ�Pt2�ķ�֧���
   Decs=Pt2(First).decs;                        %�õ�Pt2��֧����POS
   [m,index]=max(std(Decs,0,1));            %��ÿһά���� �ķ�����ط������ֵm����������index
   mainD=index;                               %�������Ϊ ��ά��
   DV2=[];                                   %Pt2ÿһά�Ĺ���ֵ
  
   for i =1:Global.D
        a=Decs(:,i);
        [f,x]=ksdensity(a);             %���ܶȹ���
        [m,ind]=max(f);                  
        xi=x(ind);                         %�����iά����, �������Ĺ���ֵ 
        DV2=[DV2,xi];    
   end
   
   [FrontNo,MaxFNo] = NDSort(Pt1.objs,Global.N);
   First = find(FrontNo==1);             %�õ�Pt1�ķ�֧���
   Decs=Pt1(First).decs;                        %�õ�Pt1��֧����POS
   [m,index]=max(std(Decs,0,1));            %��ÿһά���� �ķ�����ط������ֵm����������index
   mainD=index;                              %�������Ϊ ��ά��
   DV1=[];                                   %Pt1ÿһά�Ĺ���ֵ
   
   for i =1:Global.D
        a=Decs(:,i);
        [f,x]=ksdensity(a);             %���ܶȹ���
        [m,ind]=max(f);                  
        xi=x(ind);                         %�����iά����, �������Ĺ���ֵ 
        DV1=[DV1,xi];    
   end
   DV3=DV2+(DV2-DV1);
   DV3upper=DV3+(DV3-DV2);                 %��һ������ ���������Ͻ�
   DV3lower=DV2;                           %��һ������ ���������½�
   
   P3=[];
    for i = 1 : N1
        for j = 1:N2
            Decs=[];
            for k =1:Global.D
                if k==mainD                              %�����ά������ά��
                    temp=Global.lower(k)+((Global.upper(k)-Global.lower(k))*i)/N1;
                else
                    temp=DV3lower(k)+((DV3upper(k)-DV3lower(k))*j)/N2;
                end
            end
            Decs=[Decs,temp];
            p=INDIVIDUAL(Decs);
            P3=[P3,p];
        end
    end 
  [FrontNo,MaxFNo] = NDSort(P3.objs,Global.N);
  First = find(FrontNo==1); 
  P3=P3(First);
end

function isture = Change(Global)
    isture= rem(Global.Current,Global.tao_n) == 0;
end
