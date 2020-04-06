%********************DYNAMIC OPTIMIZATION*********************************%
%*****   Paper:Collaborative Beamforming for Wireless Sensor Networks*****%
%*****   with Fog Computing: An Dynamic Learning Approach            *****%
%*****    Performance Analysis of CCDF                               *****%
%*****   The Proposed DAL algorithm \eta=10 \eta=20                  *****%
%*****   分析CCDF性能                                                ******%
% 参数: the directions of unintended BS: -150,-70,100
%       M=30 to 80,N>=M/2, \eta=10
% 日期：2019.08.16
%*************************************************************************%
%**********************Network Initial*********************%
clc; 
clear all;
close all; 
MaxV=30;
MinV=-30;
Maxiteration=800;
Best_Sta=zeros(Maxiteration,5);  %运行800次进行统计
Num_N=1;
M=60;                            %总的节点数量

for jj=1:Maxiteration

% XY=zeros(n,2);
%************random node array************
r=fix(MinV+(MaxV-MinV).*rand([1 M]));  %随机产生60个节点极坐标
N_Fai=-pi/2+pi.*rand([1 M]);
% r=cell2mat(struct2cell(load('r.mat')));
% N_Fai=cell2mat(struct2cell(load('N_Fai.mat')));
w=ones(M,1);
%************************
%************location of intended sink node****************
Theta0=pi/2;
%BFai0=pi*5/8;
BFai0=0;
%************************************************************
r=abs(r);
lemda=15;
%**********Unintended direction******************************
% BFai=[pi/4-pi/8,pi/4+pi/8,5*pi/8,-pi/4];
pFai=linspace(-pi,pi,361);
BFai=[pFai(31),pFai(111),pFai(281)];%-150,-70,100
D=3;

%*********************Initialization***********************************
p=0.5.*ones(1,M);   %初始每个节点选择的概率，0.5
%**************************************************************************
N=M/2;              %要求选择的节点数量


MaxEvalNum_T=50;
%***************Proposed DLA -10 ************************
        Tem=0.01*ones(1,M);
        Eval_t=1;
        ut1=zeros(1,M);
        ut0=zeros(1,M);


        SumValue1=zeros(1,M);
        SumValue0=zeros(1,M);

        Act=zeros(1,M);
        rou=0.3;   % important parameter , 0.3 is great.
        E_num=10;  %the number of equlibrium********************* IMPORTANCE************************
        ACTA=zeros(E_num,M);
        Sibe_Enum=zeros(1,E_num);
        for s_e=1:E_num

            for n=1:M
               Act(n)=binornd(1,p(n)); 
            end  
            %Act
            ACTA(s_e,:)=Act;
            %Sidelobe=zeros(1,MaxEvalNum_T);

            Eval_t=0;
        %C_Uti=FitSibelobeBXC(Act,N,r,N_Fai,Theta0,BFai0,BFai,M); %current utility
            while  (Eval_t<MaxEvalNum_T)
                 TemSaveAct=0;
                 for i=1:M
                     TemSaveAct=Act(i);
                     Act(i)=1;
                     ut1(i)=FitSibelobeBXC(Act,M,r,N_Fai,Theta0,BFai0,BFai,D);
                     Act(i)=0;
                     ut0(i)=FitSibelobeBXC(Act,M,r,N_Fai,Theta0,BFai0,BFai,D);
                     Act(i)=TemSaveAct;

                    % SumValue0(i)=(Eval_t-1)/Eval_t*SumValue0(i)+1/Eval_t*ut0(i);
                    % SumValue1(i)=(Eval_t-1)/Eval_t*SumValue1(i)+1/Eval_t*(ut1(i)-ut1(i)/(1.5*N));
                     SumValue0(i)=(1-rou)*SumValue0(i)+rou*ut0(i);
                     SumValue1(i)=(1-rou)*SumValue1(i)+rou*(ut1(i)-ut1(i)/(1.5*M));
                     Temm=randsrc(1,1,[0,1;0.6,0.4]);

                     if Temm==0
                         [Value,vv]=min([SumValue0(i),SumValue1(i)]);
                         ACTA(s_e,i)=vv-1;
                     else
                         ACTA(s_e,i)=Act(i);
                     end
                 end

                 if Act==ACTA(s_e,:)
                     break;
                 end
                 Act=ACTA(s_e,:);
                 Eval_t=Eval_t+1;  
            end
            Sibe_Enum(s_e)=FitSibelobeBXC(ACTA(s_e,:),M,r,N_Fai,Theta0,BFai0,BFai,D);

        end

        [Va,Id]=min(Sibe_Enum);
        Act=ACTA(Id,:);

%****************Proposed DLA -20*****************************************

        Tem=0.01*ones(1,M);
        Eval_t=1;
        ut1=zeros(1,M);
        ut0=zeros(1,M);


        SumValue1=zeros(1,M);
        SumValue0=zeros(1,M);
        Act2=zeros(1,M);
        rou=0.3;   % important parameter , 0.3 is great.
        E_num=20;  %the number of equlibrium********************* IMPORTANCE************************
        ACTA=zeros(E_num,M);
        Sibe_Enum=zeros(1,E_num);
        for s_e=1:E_num

            for n=1:M
               Act2(n)=binornd(1,p(n)); 
            end  
            %Act
            ACTA(s_e,:)=Act2;
            %Sidelobe=zeros(1,MaxEvalNum_T);

            Eval_t=0;
        %C_Uti=FitSibelobeBXC(Act,N,r,N_Fai,Theta0,BFai0,BFai,M); %current utility
            while  (Eval_t<MaxEvalNum_T)
                 TemSaveAct=0;
                 for i=1:M
                     TemSaveAct=Act2(i);
                     Act2(i)=1;
                     ut1(i)=FitSibelobeBXC(Act2,M,r,N_Fai,Theta0,BFai0,BFai,D);
                     Act2(i)=0;
                     ut0(i)=FitSibelobeBXC(Act2,M,r,N_Fai,Theta0,BFai0,BFai,D);
                     Act2(i)=TemSaveAct;

                    % SumValue0(i)=(Eval_t-1)/Eval_t*SumValue0(i)+1/Eval_t*ut0(i);
                    % SumValue1(i)=(Eval_t-1)/Eval_t*SumValue1(i)+1/Eval_t*(ut1(i)-ut1(i)/(1.5*N));
                     SumValue0(i)=(1-rou)*SumValue0(i)+rou*ut0(i);
                     SumValue1(i)=(1-rou)*SumValue1(i)+rou*(ut1(i)-ut1(i)/(1.5*M));
                     Temm=randsrc(1,1,[0,1;0.6,0.4]);

                     if Temm==0
                         [Value,vv]=min([SumValue0(i),SumValue1(i)]);
                         ACTA(s_e,i)=vv-1;
                     else
                         ACTA(s_e,i)=Act2(i);
                     end
                 end

                 if Act2==ACTA(s_e,:)
                     break;
                 end
                 Act2=ACTA(s_e,:);
                 Eval_t=Eval_t+1;  
            end
            Sibe_Enum(s_e)=FitSibelobeBXC(ACTA(s_e,:),M,r,N_Fai,Theta0,BFai0,BFai,D);

        end

        [Va,Id]=min(Sibe_Enum);
        Act2=ACTA(Id,:);


%***************************CEO**********************************

%**************Inital Probility of selecting Action(Transit Power) for each node**************
p=0.5.*ones(1,M);
%**************************************************************************
pFai=linspace(-pi,pi,361);
Un_ID=[31,111,281];
% BFai=[pFai(31),pFai(111),pFai(281)];%-150,-70,110
D=3;
Maxite=100;
L=60;
g=zeros(L,M);
g_cb=zeros(1,L);
rau=0.1;
a=0.8;
for i=1:Maxite
    g_cb=zeros(1,L);
    for l=1:L
        while(1)
            for n=1:M
                g(l,n)=binornd(1,p(n));
            end
          if sum(g(l,:)~=0)>=N
              break;
          end
        end
        g_cb(l)=FitSibelobeBXC(g(l,:),M,r,N_Fai,Theta0,BFai0,BFai,D);
        g_cb(l)=1/D*g_cb(l); 
    end
    Ng_cb=sort(g_cb,'ascend');
    for nn=1:L            %find minium value corrsponding the position :id 
        if Ng_cb(1)==g_cb(nn)
            %N_id=nn;
            N_g=g(nn,:);
            break;
        end
    end    
    rr=Ng_cb(ceil(rau*L));
    sp=find(g_cb<rr);
    if isempty(sp)
        break;
    end
    for m=1:M
        sum_gm=0;
        p_old=p(m);
        for ll=1:length(sp)
            sum_gm=sum_gm+g(sp(ll),m);
        end
        p(m)=1/length(sp)*sum_gm;
        p(m)=a*p(m)+(1-a)*p_old;
    end
end


%*********************radmon node selection******************************%
rg=zeros(L,M);
rg_cb=zeros(1,L);
for i=1:L
    while(1)
          for n=1:M
             rg(i,n)=binornd(1,0.5);
          end
          if sum(rg(i,:)~=0)>=N
             break;
          end
    end
    rg_cb(i)=FitSibelobeBXC(rg(i,:),M,r,N_Fai,Theta0,BFai0,BFai,D);
end
[mv,idd]=min(rg_cb);
zrg=rg(idd,:);
 


%***************************without optimization*************************%
g=zeros(1,M);
while(1)
  for n=1:M
     g(n)=binornd(1,0.5);
  end
  if sum(g(1,:)~=0)>=N
         break;
   end
end

%**************计算非目标基站点的平均INR***********************************%
 AFR=FitnessFunO(N_g.*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0);  %CEO
 ceo=(abs(AFR(31))^2+abs(AFR(111))^2+abs(AFR(281))^2)/D;
 Best_Sta(jj,1)=ceo;
 AFR2=FitnessFunO(Act.*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0); %MY-10
 bxc=(abs(AFR2(31))^2+abs(AFR2(111))^2+abs(AFR2(281))^2)/D
 Best_Sta(jj,2)=bxc;
 
 AFR3=FitnessFunO(Act2.*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0);%MY-20
 bxc2=(abs(AFR3(31))^2+abs(AFR3(111))^2+abs(AFR3(281))^2)/D
 Best_Sta(jj,3)=bxc2;
  
 AFR4=FitnessFunO(zrg.*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0); %RANDOM
 RSA=(abs(AFR4(31))^2+abs(AFR4(111))^2+abs(AFR4(281))^2)/D;
 Best_Sta(jj,4)=RSA;
 AFR5=FitnessFunO(g.*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0);   %WITHOUT
 OWO=(abs(AFR5(31))^2+abs(AFR5(111))^2+abs(AFR5(281))^2)/D;
 Best_Sta(jj,5)=OWO;
 
 
end


    
    