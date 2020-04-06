%********************DYNAMIC OPTIMIZATION*********************************%
%*****   Paper:Collaborative Beamforming for Wireless Sensor Networks*****%
%*****   with Fog Computing: An Dynamic Learning Approach            *****%
%*****    Performance Analysis of CCDF                               *****%
%*****   The Proposed DAL algorithm \eta=10 \eta=20                  *****%
%*****   画图CCDF                                                   ******%
%***** 参数: the directions of unintended BS: -150,-70,100          ******%
%***** M=30 to 80,N>=M/2, \eta=10                                   ******%
%***** 日期：2019.08.16                                             ******%
%*************************************************************************%
clc;
clear all;
D=3;
%CEO
M=60;
N=30;
Best_StaN=cell2mat(struct2cell(load('Best_Sta.mat')));


%********************统计CEO*********************************************%
[x1,Num]=size(Best_StaN(:,1));
Amp_size=[0.01:0.001:1];              %INR 的横坐标值
db_Value=10*log10(Amp_size);          %对应的dB功率值
[x,y]=size(db_Value);                 %计算对应的矢量大小
dbPro=zeros(1,y);                     %用于存储大于db_Value中每个INR dB值对应的平均统计值
t=1;
tt=db_Value(1:25:end);
pp=1:25:y;
for j=0.01:0.001:1
    
   dbPro(t)=sum(Best_StaN(:,1)>j)/x1; %计算大于db_Value中每个INR dB值对应的平均统计值
   t=t+1;

end
Chan=find(dbPro<0.001);               %统计值太小的话，赋0       
dbPro(Chan)=0;
semilogy(db_Value,dbPro);             %画图
hold on
h1=semilogy(tt,dbPro(pp),'s');        %再画出相应的标识
hold on
axis([-20,40,0.01,1]);
% c
[zzz,ddd]=min(dbPro(find(dbPro>0.01))); %测试用
db_Value(ddd)



%
%********************MY-10************************************************%
[x1,Num]=size(Best_StaN(:,2));
Amp_size=[0.01:0.001:1];
db_Value=10*log10(Amp_size);
[x,y]=size(db_Value);
dbPro=zeros(1,y);
t=1;
tt=db_Value(1:20:end);
pp=1:20:y;
for j=0.01:0.001:1
    
   dbPro(t)=sum(Best_StaN(:,2)>j)/x1;
   t=t+1;

end

Chan=find(dbPro<0.001);
dbPro(Chan)=0;
hold on;
semilogy(db_Value,dbPro);
hold on
h2=semilogy(tt,dbPro(pp),'+');
[zzz,ddd]=min(dbPro(find(dbPro>0.01)));
db_Value(ddd)


%********************MY-20************************************************%
[x1,Num]=size(Best_StaN(:,1));
Amp_size=[0.01:0.001:1];
db_Value=10*log10(Amp_size);
[x,y]=size(db_Value);
dbPro=zeros(1,y);
t=1;
tt=db_Value(1:20:end);
pp=1:20:y;
for j=0.01:0.001:1
    
   dbPro(t)=sum(Best_StaN(:,3)>j)/x1;
   t=t+1;

end

Chan=find(dbPro<0.001);
dbPro(Chan)=0;
hold on;
semilogy(db_Value,dbPro);
hold on
h3=semilogy(tt,dbPro(pp),'x');
[zzz,ddd]=min(dbPro(find(dbPro>0.01)));
db_Value(ddd)



%********************RANDOM***********************************************%
[x1,Num]=size(Best_StaN(:,4));
Amp_size=[0.01:0.001:19];
db_Value=10*log10(Amp_size);
[x,y]=size(db_Value);
dbPro=zeros(1,y);
t=1;
tt=db_Value(1:500:end);
pp=1:500:y;
for j=0.01:0.001:19
    
   dbPro(t)=sum(Best_StaN(:,4)>j)/x1;
   t=t+1;

end

Chan=find(dbPro<0.001);
dbPro(Chan)=0;
hold on;
semilogy(db_Value,dbPro);
hold on
h4=semilogy(tt,dbPro(pp),'*');
[zzz,ddd]=min(dbPro(find(dbPro>0.01)));
db_Value(ddd)

%********************Without OP*******************************************%
[x1,Num]=size(Best_StaN(:,5));
Amp_size=[0.01:0.001:200];
db_Value=10*log10(Amp_size);
[x,y]=size(db_Value);
dbPro=zeros(1,y);
t=1;
tt=db_Value(1:3500:end);
pp=1:3500:y;
for j=0.01:0.001:200
    
   dbPro(t)=sum(Best_StaN(:,5)>j)/x1;
   t=t+1;

end

Chan=find(dbPro<0.01);
dbPro(Chan)=0;
hold on;
semilogy(db_Value,dbPro);
hold on
h5=semilogy(tt,dbPro(pp),'s');
[zzz,ddd]=min(dbPro(find(dbPro>0.01)));
db_Value(ddd)


%加上legend
legend([h1,h2,h3,h4,h5],'CEO algorithm','The Proposed DAL \eta=10','The Proposed DAL \eta=20','Random node selection','Without optimization');

