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
    
    N1=9     %Knee �����
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
    KneesPoints=zeros(N1,Global.D,200);  %��ʼ����ʷknee�㼯��

    PF=Global.PF;
    while Global.NotTermination(Population)
%         Cc=INDIVIDUAL(sum(Population.decs,1)/Global.N);
%         DF(func2str(Global.problem),Global.NT);    %������ӦPF
%         Draw(Global.PF,'ok','Markeredgecolor',[.9 .0 .0],'Markerfacecolor',[.9 .0 .0],'MarkerSize',2);
%         Draw(Cc.objs,'ok','Markeredgecolor',[.0 .9 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',5);
%         Draw(Pt1.objs,'k+','Markeredgecolor',[.6 .0 .0],'Markerfacecolor',[.0 .9 .0],'MarkerSize',3);

       
       
        if Change(Global)   % ��⻷���仯
            t=t+1      
            Decs=Population.decs;
            Population=INDIVIDUAL(Decs);
            Pt2=Population;

            [FrontNo,MaxFNo] = NDSort(Pt2.objs,Global.N);
            First = find(FrontNo==1);   
            Ct=sum(Pt2(First).decs,1)/length(First);  % �����֧������ĵ�
            CenterPoint=[CenterPoint;Ct];       % ��Ct�ӽ���ʷ���ĵ�

            [P1,KneesPoints]=KneePoint(Pt2(First),Global,t,N1,KneesPoints);               %Ԥ��Knee��
            P2=prediction(Pt2(First),Global,t,CenterPoint);              %Ԥ�����ĵ�
            Population=RandGenrate(P1,P2,Global);                       %�����������
      
            PF=Global.PF;
            Pt1=Population;
            if t>=2      %��ȥ��һ�λ����仯
                 score=IGD(Population.objs,PF);      %����IGD
                 metric=[metric,score];
                 migd=sum(metric)/length(metric); %����MIGD
                 igd=migd;                         %��MIGDͨ��ȫ�ֱ���igd����ָ�꺯��MIGD()
            end
        else
%             MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
%             Offspring  = Global.Variation(Population(MatingPool));
              Offspring  = Global.Variation(Population,Global.N,@RMMEDA_operator);
           [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
        end
    end
end

function [P1,KneesPoints]=KneePoint(Pop,Global,t,N1,KneesPoints)    %Ԥ��Knee��
    L=length(Pop);                         %�õ���֧����Ⱥ��С
    Pop_Obj=Pop.objs;
    [~,index]=sortrows(Pop_Obj,1);        %�Է�֧����Ⱥ ���ݵ�һ��Ŀ������ 
    N2=floor(L/N1);                               %����Ⱥ�ֳ�N1�ݣ�ÿ��N2�� ����
    P1=[];
    if N2~=0
        KneesD=zeros(N1,Global.D);
        Knees=[];                             %��ʼ����ǰKnee�㼯��
        for i = 1:N1                        %��ÿ�ݸ�����Knee��
            if i==N1
                Pop2_Obj=Pop_Obj(index((i-1)*N2+1:end),:);   
            else
                Pop2_Obj=Pop_Obj(index((i-1)*N2+1:i*N2),:);  
            end
            if Global.M==2                    %�����2άĿ�꣬ȷ��2���߽�㣬�������㵽ֱ�ߵľ���
                [~,index1]=max(Pop_Obj(:,1)); %��һ��Ŀ����������
                [~,index2]=max(Pop_Obj(:,2)); %�ڶ���Ŀ����������
                if  Pop_Obj(index1,1)~=Pop_Obj(index2,1) &&  Pop_Obj(index1,2)~=Pop_Obj(index2,2)         %�߽��ÿά���������
                    Point1=Pop_Obj(index1,:); %�õ��߽�� Point1
                    Point2=Pop_Obj(index2,:); %�õ��߽�� Point2

                    syms K B
                    eq1=K*Point1(1)-Point1(2)+B==0;   %����ֱ�߷���
                    eq2=K*Point2(1)-Point2(2)+B==0;   %����ֱ�߷���
                    [K,B]=solve(eq1,eq2,K,B);         %�ⷽ��

                    double(K);
                    double(B);
                    dis=abs(double(K)* Pop2_Obj(:,1) -  Pop2_Obj(:,2) +double(B))/(sqrt(double(K)^2+1));
                    [~,index3]=min(dis);
                    knee=Pop(index((i-1)*N2+index3)); 
                    Knees=[Knees,knee];                  %�ѵ�ǰ�����knee�� �ӽ�knee�㼯
                end
            elseif Global.M==3               %�����3άĿ�꣬ȷ��3���㣬�������㵽��ƽ��ľ���

            end
        end
        if ~isempty(Knees)&&length(Knees)==N1
            KneesD=Knees.decs;
            KneesPoints(:,:,t)=KneesD;         %����ǰknee�㼯�ϼӽ���ʷknee�㼯��
        end
    
        if t<=6||length(Knees)~=N1
            P1=[];
        else                                   %t����6ʱ����ARģ��Ԥ����һ������ knee�㼯��
            for i = 1:N1 
                 Knees_i=[];                     %��i������ ����ʷknee�㼯��
                 for j=1:t
                     temp=KneesPoints(i,:,j);
                     Knees_i=[Knees_i;temp];
                 end
          
                 NextKnee_i=Center_P(Knees_i,t,Global);  %Ԥ���i���� ��һ������knee��
                 P=INDIVIDUAL(NextKnee_i);
                 P1=[P1,P];
            end
        end
    end
end

function P2=prediction(Pt2,Global,t,CenterPoint)    %Ԥ�����ĵ�
     Ob=Pt2;
     OP=length(Ob);
     if t>=2
        Decs=Ob.decs+(CenterPoint(t,:)-CenterPoint(t-1,:));               %Ԥ��
        for i = 1 : OP              %����decs
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

function P3=RandGenrate(P1,P2,Global)              %����Ӧ�������
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
    %% ������ʷ�㽨��ARģ�ͣ���Ԥ�����ĵ�,�����CenterPoint��knee��
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
