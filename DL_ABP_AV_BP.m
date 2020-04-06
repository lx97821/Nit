%********************DYNAMIC OPTIMIZATION*********************************%
%********Paper:Collaborative Beamforming for Wireless Sensor Networks*****%
%********with Fog Computing: An Dynamic Learning Approach*****************%
%********Performance Analysis of Average baempattern**********************%
%********The Proposed DAL algorithm \eta=10*******************************%
%********验证30个节点的提出方法的平均波瓣性能*******************************%
% 参数: the directions of unintended BS: -150,-70,100
%       M=30,N>=15, \eta=10
% 日期：2019.08.16
%*************************************************************************%
    clc; 
    clear all;
    close all; 
    ceo=0;           %存储CEO 算法的INR值
    noregret=0;

    Mainceo=0;
    Mainnoregret=0;


    MaxV=30;        %节点x,y分布的区别最大最小值
    MinV=-30;
    M=30;           %节点总数
%**************************Load random node array************************%
    r=cell2mat(struct2cell(load('r30.mat')));
    N_Fai=cell2mat(struct2cell(load('N_Fai30.mat')));
    w=ones(M,1);
%************location of intended sink node******************************%

    Theta0=pi/2;   
    BFai0=0;
%************************************************************************%
    r=abs(r);         %极坐标半径
    lemda=15;         %波长
%**********the directions of Unintended direction************************%
    pFai=linspace(-pi,pi,361);
    BFai=[pFai(31),pFai(111),pFai(281)];%-150,-70,110
    D=3;
%***********************
%*********************Initialization***********************************
    p=0.5.*ones(1,M); % 协同波束节点的初始选择，定义基于0，1各1/2概率
    N=M/2;            % 主瓣要求最小N值
    Eval_t=1;         % 迭代变量初始值
    MaxEvalNum_T=50;  % 最大迭代次数，根据节点规模修改
    ceoAVBP=zeros(1,361);
    noregretAVBP=zeros(1,361);
    randomAVBP=zeros(1,361);
    withoutAVBP=zeros(1,361);
 

    real_nums=0;       %用于统计生成的选择节点个数大于N的次数
    for nums=1:60       %运行60次，算法平均波瓣

        ut1=zeros(1,M);%用于存储节点被选择及没被选择的非目标基站的INR
        ut0=zeros(1,M);
        SumValue1=zeros(1,M);%用于存储节点被选择及没被选择的累计值，具体看JSFP
        SumValue0=zeros(1,M);
        Act=zeros(1,M);      %存储M的节点的行为，即1，0；
        rou=0.3;   % important parameter , 0.3 is great.
        E_num=10;  % the number of equlibrium********** IMPORTANCE********
        ACTA=zeros(E_num,M);
        Sibe_Enum=zeros(1,E_num); %用于存储多个均衡的非目标节点INR值
        for s_e=1:E_num
            for n=1:M
               Act(n)=binornd(1,p(n)); %按照1/2的概率初始选择节点
            end  
            ACTA(s_e,:)=Act;           %临时保存当前节点选择情况
            Eval_t=0;
            while  (Eval_t<MaxEvalNum_T)
                 TemSaveAct=0;
                 for i=1:M
                     TemSaveAct=Act(i);  %计算每个节点选择和非选择的效用
                     Act(i)=1;
                     ut1(i)=FitSibelobeBXC(Act,M,r,N_Fai,Theta0,BFai0,BFai,D);
                     Act(i)=0;
                     ut0(i)=FitSibelobeBXC(Act,M,r,N_Fai,Theta0,BFai0,BFai,D);
                     Act(i)=TemSaveAct;
                     SumValue0(i)=(1-rou)*SumValue0(i)+rou*ut0(i);
                     SumValue1(i)=(1-rou)*SumValue1(i)+rou*(ut1(i)-ut1(i)/(3*N/M*M)); %belta=rou*ut1(i)/(3*N/M*M) ;
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
            Sibe_Enum(s_e)=FitSibelobeBXC(ACTA(s_e,:),M,r,N_Fai,Theta0,BFai0,BFai,D); %存储每个均衡的非目标INR值
        end

        [Va,Id]=min(Sibe_Enum);
        Act=ACTA(Id,:);



%+++++++++++++++++++以下是CEO算法++++++++++++++++++++++++++++++++++++++++++%

%+++++++++++++++++++++++++++CEO ALGORITHM+++++++++++++++++++++++++++++++++%
%++++算法描述具体见论文：An Efficient Sensor-Node Selection Algorithm for++%
%+++++++++Sidelobe Control in Collaborative Beamforming+++++++++++++++++++%
        p=0.5.*ones(1,M);
        pFai=linspace(-pi,pi,361);
        D=3;
        Maxite=50;
        L=60;
        g=zeros(L,M);
        g_cb=zeros(1,L);
        rau=0.1;
        a=0.8;
        ceoIterNum=0;
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
                ceoIterNum=i*L;
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
        
        
%***********************Randmon Node Selection****************************%
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
 


 %***********************without optimization*****************************%
        g=zeros(1,M);
        while(1)
          for n=1:M
             g(n)=binornd(1,0.5);
          end
          if sum(g(1,:)~=0)>=N
                 break;
           end
        end      
        
 %**********************计算各算法的非目标节点INR值之和********************%
         
        AFR2=FitnessFunO(N_g.*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0);  %ceo
        Mainceo=abs(AFR2(181)); 

        AFR=FitnessFunO(Act.*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0);   %my
        Mainnoregret=abs(AFR(181));
        
        AFR3=FitnessFunO(zrg.*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0);  %random selection
        
        AFR4=FitnessFunO(g.*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0);    %without optimization
        
%*************************************************************************%
        if (Mainceo>=N&&Mainnoregret>=N) %满足条件下的计算平均波束
            ceoAVBP=ceoAVBP+AFR2;
            noregretAVBP=noregretAVBP+AFR;
            randomAVBP= randomAVBP+AFR3;
            withoutAVBP=withoutAVBP+AFR4;
            real_nums= real_nums+1;
%             ceo=abs(AFR2(31))+abs(AFR2(111))+abs(AFR2(281));
%             noregret=abs(AFR(31))+abs(AFR(111))+abs(AFR(281));
%             ceosum=ceosum+ceo;
%             noregretsum=noregretsum+noregret;
%             if ceo<noregret
%                 ceoNum=ceoNum+1;
%             else
%                 noregretNum=noregretNum+1;
%             end
        end
        
    end
    ceoAVBP=ceoAVBP./real_nums;
    noregretAVBP=noregretAVBP./real_nums;
    randomAVBP=randomAVBP./real_nums;
    withoutAVBP=withoutAVBP./real_nums;

    plot(pFai*(180/pi),10*log10(abs(ceoAVBP).^2),'-.');
    hold on
    plot(pFai*(180/pi),10*log10(abs(noregretAVBP).^2));
    hold on
    plot(pFai*(180/pi),10*log10(abs(randomAVBP).^2));
    hold on
    plot(pFai*(180/pi),10*log10(abs(withoutAVBP).^2));



