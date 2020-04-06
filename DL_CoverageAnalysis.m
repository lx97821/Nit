%********************DYNAMIC OPTIMIZATION*********************************%
%*****   Paper:Collaborative Beamforming for Wireless Sensor Networks*****%
%*****   with Fog Computing: An Dynamic Learning Approach            *****%
%*****    Performance Analysis of coverage                           *****%
%*****   The Proposed DAL algorithm \eta=10 \eta=20                  *****%
%*****   计算不同N值下的INR值                                        ******%
%***** 参数: the directions of unintended BS: -150,-70,100          ******%
%***** M=30 to 80,N>=M/2, \eta=10                                   ******%
%***** \belta=31.623  becasue of l0*log(31.623)=15dB                ******%
%***** 日期：2019.08.16                                             ******%
%*************************************************************************%
    clc; 
    clear all;
    close all; 
    ceo=0;
    noregret=0;

    Mainceo=0;
    Mainnoregret=0;

    %for kkk=1:100
    MaxV=30;
    MinV=-30;
    M=60;
    N=M/2;
    r=fix(MinV+(MaxV-MinV).*rand([1 M]));
    N_Fai=-pi/2+pi.*rand([1 M]);
    w=ones(M,1);
    %************************
    %************location of intended sink node****************
    Theta0=pi/2;
    BFai0=0;
    %************************************************************
    r=abs(r);
    lemda=15;
    %**********Unintended direction******************************
    pFai=linspace(-pi,pi,361);
    BFai=[pFai(31),pFai(111),pFai(281)];%-150,-70,110
    D=3;
    %***********************
    %*********************Initialization***********************************
    %****************Energy Initialization****************

    p=0.5.*ones(1,M);
    Eval_t=1;
    MaxEvalNum_T=100;
    ceoNum=0;
    noregretNum=0;
    ut1=zeros(1,M);
    ut0=zeros(1,M);


    SumValue1=zeros(1,M);
    SumValue0=zeros(1,M);

    Act=zeros(1,M);
    rou=0.3;
    E_num=10;  %the number of equlibrium********************* IMPORTANCE************************
    Coverage_A=zeros(E_num,MaxEvalNum_T);  %收敛值存储
    ACTA=zeros(E_num,M);
    Sibe_Enum=zeros(1,E_num);
     
    for s_e=1:E_num
        for n=1:M
           Act(n)=binornd(1,p(n)); 
        end  
        ACTA(s_e,:)=Act;
        Sidelobe=zeros(1,MaxEvalNum_T);
        
        Eval_t=0;
        % C_Uti=FitSibelobeBXC(Act,N,r,N_Fai,Theta0,BFai0,BFai,M) %current utility
        while  (Eval_t<MaxEvalNum_T)
            TemSaveAct=0;
            for i=1:M
                TemSaveAct=Act(i);
                Act(i)=1;
                ut1(i)=FitSibelobeBXC(Act,M,r,N_Fai,Theta0,BFai0,BFai,D);
                Act(i)=0;
                ut0(i)=FitSibelobeBXC(Act,M,r,N_Fai,Theta0,BFai0,BFai,D);
                Act(i)=TemSaveAct;
                SumValue0(i)=(1-rou)*SumValue0(i)+rou*ut0(i);
                SumValue1(i)=(1-rou)*SumValue1(i)+rou*(ut1(i)-ut1(i)/(3*M/N*M));
                Temm=randsrc(1,1,[0,1;0.5,0.5]);
                if Temm==0
                    [Value,vv]=min([SumValue0(i),SumValue1(i)]);
                    ACTA(s_e,i)=vv-1;
                else
                    ACTA(s_e,i)=Act(i);
                end
               
            end
            
            %                  if Act==ACTA(s_e,:)
            %                      break;
            %                  end
            Act=ACTA(s_e,:);
            Eval_t=Eval_t+1;
          
            my10AFR=FitnessFunO(Act.*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0);
            Coverage_A(s_e,Eval_t)=(abs(my10AFR(31))^2+abs(my10AFR(111))^2+abs(my10AFR(281))^2)/D;
            
          
        end
        Sibe_Enum(s_e)=FitSibelobeBXC(ACTA(s_e,:).*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0,BFai,D);
        
    end
    
    [Va,Id]=min(Sibe_Enum);
    Act=ACTA(Id,:);


    
    AFR=FitnessFunO(Act.*sqrt(31.623/N),M,r,N_Fai,Theta0,BFai0);    
    figure(1)
    plot(pFai*(180/pi),10*log10(abs(AFR).^2))       
    figure(2)
    plot(1:MaxEvalNum_T,10*log10(Coverage_A(:,1:MaxEvalNum_T)),'LineWidth',1);

    noregretNum=(abs(AFR(31))^2+abs(AFR(111))^2+abs(AFR(281))^2)/D;
% plot(pFai*(180/pi),10*log10(abs(AFR2)));
% 
% 
% hold on
% 
% 
% plot(pFai*(180/pi),10*log10(abs(AFR)));
% 

%end


