%********************DYNAMIC OPTIMIZATION*********************************%
%*****   Paper:Collaborative Beamforming for Wireless Sensor Networks*****%
%*****   with Fog Computing: An Dynamic Learning Approach            *****%
%*****    Performance Analysis of CCDF                               *****%
%*****   The Proposed DAL algorithm \eta=10 \eta=20                  *****%
%*****   ��ͼCCDF                                                   ******%
%***** ����: the directions of unintended BS: -150,-70,100          ******%
%***** M=30 to 80,N>=M/2, \eta=10                                   ******%
%***** ���ڣ�2019.08.16                                             ******%
%*************************************************************************%
clc;
clear all;
D=3;
%CEO
M=60;
N=30;
Best_StaN=cell2mat(struct2cell(load('Best_Sta.mat')));


%********************ͳ��CEO*********************************************%
[x1,Num]=size(Best_StaN(:,1));
Amp_size=[0.01:0.001:1];              %INR �ĺ�����ֵ
db_Value=10*log10(Amp_size);          %��Ӧ��dB����ֵ
[x,y]=size(db_Value);                 %�����Ӧ��ʸ����С
dbPro=zeros(1,y);                     %���ڴ洢����db_Value��ÿ��INR dBֵ��Ӧ��ƽ��ͳ��ֵ
t=1;
tt=db_Value(1:25:end);
pp=1:25:y;
for j=0.01:0.001:1
    
   dbPro(t)=sum(Best_StaN(:,1)>j)/x1; %�������db_Value��ÿ��INR dBֵ��Ӧ��ƽ��ͳ��ֵ
   t=t+1;

end
Chan=find(dbPro<0.001);               %ͳ��ֵ̫С�Ļ�����0       
dbPro(Chan)=0;
semilogy(db_Value,dbPro);             %��ͼ
hold on
h1=semilogy(tt,dbPro(pp),'s');        %�ٻ�����Ӧ�ı�ʶ
hold on
axis([-20,40,0.01,1]);
% c
[zzz,ddd]=min(dbPro(find(dbPro>0.01))); %������
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


%����legend
legend([h1,h2,h3,h4,h5],'CEO algorithm','The Proposed DAL \eta=10','The Proposed DAL \eta=20','Random node selection','Without optimization');

