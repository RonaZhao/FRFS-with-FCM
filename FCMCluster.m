%����ģ��C��ֵ�����ݼ�data��Ϊn��
%�÷�
% 1 [center,U,obj_fcn]=FCMCluster(data,n,options);
% 2 [center,U,obj_fcn]=FCMCluster(data,n);

%���� 
% data    n*m����n����������ÿ��������ά��Ϊm
% n       �����
% options 4*1 ����
%   options(1):�����Ⱦ���U�ļ�Ȩָ��
%   options(2):����������
%   options(3):��������С�仯����������ֹ����
%   options(4):ÿ�ε����Ƿ������Ϣ��־

%���
% center    ��������
% U         �����Ⱦ���
% obj_fun   Ŀ�꺯��ֵ


function [center,U,obj_fun] = FCMCluster(data,n,options)
    if nargin~=2 && nargin~=3
        error('Too many or too few input arguments');
    end 
    data_n=size(data,1);%���ؾ��������
    in_n=size(data,2);
    %disp(data_n)
    %Ĭ�ϲ���
    default_options=[2;100;1e-5;0];

    %��������
      %���ֻ����ǰ����������ѡ��Ĭ�ϵĲ���;�����������С��4������ѡ��Ĭ�ϲ���
      if nargin==2
          options=default_options;
      else
           if length(options)<4
               tmp=default_options;
               tmp(1:length(options))=options;
               options=tmp;
           end 
           nan_index=find(isnan(options)==1);
           options(nan_index)=default_options(nan_index);

           if options(1)<=1
               error('The exponent should be greater than 1!');
           end 
      end 

      %��options �еķ����ֱ�ֵ���ĸ�����
      expo=options(1);
      max_iter=options(2);
      min_impro=options(3);
      display=options(4);

      obj_fun=zeros(max_iter,1);

      %��ʼ��ģ���������
      U=initfcm(n,data_n);
      
      %������
       for i=1:max_iter
           [U,center,obj_fun(i)]=stepfcm(data,U,n,expo);
           if display
               fprintf('FCM:Iteration count=%d,obj_fun=%f\n',i,obj_fun(i));
           end
           %��ֹ�����б�
           if i>1
               if abs(obj_fun(i)-obj_fun(i-1))<min_impro
                   break;
               end 
           end 
       end
       iter_n=i;
       obj_fun(iter_n+1:max_iter)=[];
 end
%%�Ӻ��� ģ�������ʼ��
function U= initfcm(n,data_n)
    U=rand(n,data_n);
    col_sum=sum(U);
    U=U./col_sum(ones(n,1),:);
end


%%�Ӻ��� �𲽾���
function [U_new,center,obj_fun]=stepfcm(data,U,n,expo)
	mf=U.^expo;
	center=mf*data./((ones(size(data,2),1)*sum(mf'))');
	dist=distfcm(center,data);
	obj_fun=sum(sum((dist.^2).*mf));
    tmp=dist.^(-2/(expo-1));
    U_new=tmp./(ones(n,1)*sum(tmp));
end


%%�Ӻ��� �������
function out=distfcm(center,data)
	out=zeros(size(center,1),size(data,1));
	for k=1:size(center,1)
        out(k,:)=sqrt(sum(((data-ones(size(data,1),1)*center(k,:)).^2)',1));
    end
end
            
