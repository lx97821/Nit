%********************DYNAMIC OPTIMIZATION*********************************%
%*****   Paper:Collaborative Beamforming for Wireless Sensor Networks*****%
%*****   with Fog Computing: An Dynamic Learning Approach            *****%
%*****    Performance Analysis of Energy Consumption                 *****%
%*****   The Proposed DAL algorithm \eta=10 \eta=20                  *****%
%*****   计算60个节点的累计消耗能量                                  ******%
%***** 参数: the directions of unintended BS: -150,-70,100          ******%
%***** M=60 ,N>=M/2, \eta=10                                        ******%
%***** \belta=31.623  becasue of l0*log(31.623)=15dB                ******%
%***** 日期：2019.08.20                                             ******%
%*************************************************************************%

    clc; 
    clear all;
    close all; 
    ceo=0;
    noregret=0;                 %统计Average INR值用
    Mainceo=0;
    Mainnoregret=0;
    
    MaxV=30;
    MinV=-30;
    Sta_Num=200;                %每个节点规模运行200次
    
    M=60;                %转化成节点规模
    N=M/2;        
    r=fix(MinV+(MaxV-MinV).*rand([1 M]));  %随机产生M个节点
    N_Fai=-pi/2+pi.*rand([1 M]);
    %***********************Eenergy Parameter****************************%
    Range_Distance=50;      %最大传输50m
    Tr_E=(50+1.3*10^(-6)*(Range_Distance)*10^4)*10^(-9);
    Re_E=50*10^(-9);
  
    %**********Variable Initialization************************************%
    
    %Fog-based -my10-my20 一样
     TotalEnergy_Fog=zeros(1,Sta_Num);
     Handshake_Fog=0;
     Nodestate_Update=0;
     CSI_Acq=0;
     Result_Feedback=0;
     Data_Share_Fog=0;
     Data_Tra_Fog=0;
     
       
    
    
    
    % Traditional Custer-based my-10
    TotalEnergy_Custer=zeros(1,Sta_Num);
    HakeEnergy=0;
    Dynamic_w0=0;
    Dynamic_wM=0;
    Datashare_Energy=0;
    DataTrans_Energy=0;
    
    Total_my_Num=0;
    
    % Traditional Custer-based my-20
    TotalEnergy_Custer2=zeros(1,Sta_Num);
    HakeEnergy2=0;
    Dynamic_w02=0;
    Dynamic_wM2=0;
    Datashare_Energy2=0;
    DataTrans_Energy2=0;
    
    
    
  
    
        
    for jj=1:Sta_Num
        
        Handshake_Fog=Tr_E*(50*8+39*8)+Re_E*(50*8+39*8); 
        TotalEnergy_Fog(jj)=TotalEnergy_Fog(jj)+Handshake_Fog;
        
        
        %*******************Traditional Handshke stage***********************************%
        % Each node receives a CTC and CT message and sends a RTC and ACK message 
        HakeEnergy=Tr_E*(50*8+39*8)+Re_E*(50*8+39*8); 
        TotalEnergy_Custer(jj)=TotalEnergy_Custer(jj)+HakeEnergy;
        %************location of intended sink node****************
        Theta0=pi/2;
        BFai0=0;
        %************************************************************
        r=abs(r);
        lemda=15;      %波长
        %**********the directions of D Unintended direction***********%
        pFai=linspace(-pi,pi,361);
        BFai=[pFai(31),pFai(111),pFai(281)];%-150,-70,110
        D=3;               % 非目标节点的数量
        %***********************
        %*********************DLA-My10-Initialization*****************%
        p=0.5.*ones(1,M);
        Eval_t=1;
        MaxEvalNum_T=50;  %最大迭代数量
        ceoNum=0;         %测试用，用于统计最优解数量
        noregretNum=0;
        noregretIterNum=0;
        
        ut1=zeros(1,M);   %存储节点选择0，1的效用值
        ut0=zeros(1,M);
        
        SumValue1=zeros(1,M);  %存储节点选择0，1的累计效用
        SumValue0=zeros(1,M);
        
        Act=zeros(1,M);
        rou=0.3;              %very important ,具体见JFSP
        E_num=10;             %the number of equlibrium
        ACTA=zeros(E_num,M);
        Sibe_Enum=zeros(1,E_num);
        
        %*********Nodestate Updating*********************%
        %* 2*Node stae ,2*ACK, 
        Nodestate_Update=Tr_E*((50+30)*8*2)+Re_E*(39*8*2);
        TotalEnergy_Fog(jj)=TotalEnergy_Fog(jj)+Nodestate_Update;

        
        %*************CSI Acqu*********************************%
        CSI_Acq=Tr_E*(50+50)*8+ Tr_E*((50+50)*M*8+39*8);
        
         TotalEnergy_Fog(jj)=TotalEnergy_Fog(jj)+CSI_Acq;
        
        for s_e=1:E_num
            
            for n=1:M
                Act(n)=binornd(1,p(n));
            end
            %Act
            ACTA(s_e,:)=Act;
            Sidelobe=zeros(1,MaxEvalNum_T);
            Eval_t=0;
            %*****************Dynamic Energy**************************%
            % w^0 energy ,each node sends Coll-TestData***************%
            % Testdata=100bytes
            Dynamic_w0=Tr_E*((50+50)*8);
            TotalEnergy_Custer(jj)=TotalEnergy_Custer(jj)+Dynamic_w0;
            
            while  (Eval_t<MaxEvalNum_T)
                TemSaveAct=0;
                for i=1:M
                    TemSaveAct=Act(i); %计算每个节点选择和非选择的效用
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
                
                %*****************Dynamic Energy**************************%
                %each node sends M*Coll-TestData+ACK_Utili+NT/ET*************%
                % Testdata=50bytes
                Dynamic_wM=Tr_E*((50+50)*8+(39+M+39)*8);
                TotalEnergy_Custer(jj)=TotalEnergy_Custer(jj)+Dynamic_wM;
                                
                if Act==ACTA(s_e,:)
                    bb=0;
                    Total_my_Num=Total_my_Num+Eval_t;
                    break;
                end
                
                Act=ACTA(s_e,:);
                Eval_t=Eval_t+1;
            end
            Sibe_Enum(s_e)=FitSibelobeBXC(ACTA(s_e,:),M,r,N_Fai,Theta0,BFai0,BFai,D); %存储每个均衡的非目标INR值
            
        end
        [Va,Id]=min(Sibe_Enum);
        Act=ACTA(Id,:);
        
        
        AFR1=FitnessFunO(Act.*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0);
        AFRMY10=(abs(AFR1(31))^2+abs(AFR1(111))^2+abs(AFR1(281))^2)/D;
         %*****************Datashare+DataTrans**************************%
         %each node receives Coll-Data+ Sends ACK, and sends Coll-Data**%
        Datashare_Energy=Re_E*((50+1000)*8)+Tr_E*(39*8);
        DataTrans_Energy=Tr_E*((50+1000)*8);
        TotalEnergy_Fog(jj)=TotalEnergy_Fog(jj)+Datashare_Energy+DataTrans_Energy;
        TotalEnergy_Custer(jj)=TotalEnergy_Custer(jj)+Datashare_Energy+DataTrans_Energy;
        %Total_my_Num
        
        
 %***************DLA-MY-20 Initialization**********************************
 
        %*******************Handshke stage***********************************%
        % Each node receives a CTC and CT message and sends a RTC and ACK message 
        HakeEnergy2=Tr_E*(50*8+39*8)+Re_E*(50*8+39*8); 
        TotalEnergy_Custer2(jj)=TotalEnergy_Custer2(jj)+HakeEnergy2;
        %************location of intended sink node****************
        Theta0=pi/2;
        BFai0=0;
        %************************************************************
        r=abs(r);
        lemda=15;      %波长
        %**********the directions of D Unintended direction***********%
        pFai=linspace(-pi,pi,361);
        BFai=[pFai(31),pFai(111),pFai(281)];%-150,-70,110
        D=3;               % 非目标节点的数量
        %***********************
        %*********************DLA-My10-Initialization*****************%
        p=0.5.*ones(1,M);
        Eval_t=1;
        MaxEvalNum_T=50;  %最大迭代数量
        ceoNum=0;         %测试用，用于统计最优解数量
        noregretNum=0;
        noregretIterNum=0;
        
        ut1=zeros(1,M);   %存储节点选择0，1的效用值
        ut0=zeros(1,M);
        
        SumValue1=zeros(1,M);  %存储节点选择0，1的累计效用
        SumValue0=zeros(1,M);
        
        Act=zeros(1,M);
        rou=0.3;              %very important ,具体见JFSP
        E_num=20;             %the number of equlibrium
        ACTA=zeros(E_num,M);
        Sibe_Enum=zeros(1,E_num);
        for s_e=1:E_num
            
            for n=1:M
                Act(n)=binornd(1,p(n));
            end
            %Act
            ACTA(s_e,:)=Act;
            Sidelobe=zeros(1,MaxEvalNum_T);
            Eval_t=0;
            %*****************Dynamic Energy**************************%
            % w^0 energy ,each node sends Coll-TestData***************%
            % Testdata=100bytes
            Dynamic_w02=Tr_E*((50+50)*8);
            TotalEnergy_Custer2(jj)=TotalEnergy_Custer2(jj)+Dynamic_w0;
            while  (Eval_t<MaxEvalNum_T)
                TemSaveAct=0;
                for i=1:M
                    TemSaveAct=Act(i); %计算每个节点选择和非选择的效用
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
                
%*********************************Dynamic Energy**************************%
%****************each node sends M*Coll-TestData+ACK_Utili+NT/ET**********%
%****************      Testdata=100bytes             *********************%
                Dynamic_wM2=Tr_E*((50+50)*8+(39+M+39)*8);
                TotalEnergy_Custer2(jj)=TotalEnergy_Custer2(jj)+Dynamic_wM2;
                                
                if Act==ACTA(s_e,:)
                    bb=0;
                    break;
                end
                
                Act=ACTA(s_e,:);
                Eval_t=Eval_t+1;
            end
            Sibe_Enum(s_e)=FitSibelobeBXC(ACTA(s_e,:),M,r,N_Fai,Theta0,BFai0,BFai,D); %存储每个均衡的非目标INR值
            
        end
        [Va,Id]=min(Sibe_Enum);
        Act=ACTA(Id,:);
        AFR1=FitnessFunO(Act.*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0);
        AFRMY10=(abs(AFR1(31))^2+abs(AFR1(111))^2+abs(AFR1(281))^2)/D;
%****************************Datashare+DataTrans**************************%
%**********each node receives Coll-Data+ Sends ACK, and sends Coll-Data***%
        Datashare_Energy2=Re_E*((50+1000)*8)+Tr_E*(39*8);
        DataTrans_Energy2=Tr_E*((50+1000)*8);
        TotalEnergy_Custer2(jj)=TotalEnergy_Custer2(jj)+Datashare_Energy2+DataTrans_Energy2;
  
 
    end
    
    AllEnergy=[TotalEnergy_Fog;TotalEnergy_Custer;TotalEnergy_Custer2];
    for tt=1:3 
      for ii=1:Sta_Num
             if ii>=2
         
                 AllEnergy(tt,ii)=AllEnergy(tt,ii)+AllEnergy(tt,ii-1);
             end
        
      end
    end
    
    plot(1:Sta_Num,10*log10(AllEnergy.*10^3));
    hold on 

    NewAllEnergy=zeros(3,20);
    
    for i=1:20
       if i==1
           NewAllEnergy(:,i)=AllEnergy(:,i);
       else
           NewAllEnergy(:,i)=AllEnergy(:,(i-1)*10);
       end
    end
    
    h1=plot(1:10:200,10*log10(NewAllEnergy(1,:).*10^3),'+');
    hold on
    h2=plot(1:10:200,10*log10(NewAllEnergy(2,:).*10^3),'.');
    hold on 
    h3=plot(1:10:200,10*log10(NewAllEnergy(3,:).*10^3),'*');
    
    legend([h1,h2,h3],'The proposed DAL \eta=10 and \eta=20 based on fog computing','The proposed DAL \eta=10 based on cluster achitecture','The proposed DAL \eta=20 based on cluster achitecture');
    