%programed by nan in 2018/10/16
%fuzzy rough feature selection
%input:filename of the dataset
%output:the selected subset of features
% data format
% object index |feature|label
%in=input('please input the path of dataset:');
clear;
global data;
tic;
delimiterIn =  ',';%interval of data
headerlinesIn = 1;
in='Data1/wine.csv';
a=1;%%%%%����ʵ���������
b=0;
R1=0;
data = csvread(in,headerlinesIn);
% data = csvread(in,1,R1,[1,R1,5000,3]);%���ļ���������
[r,c]=size(data);%row and column of dataset
data=sortrows(data,c);%ascending order of the label

%====================================
labellist=zeros(1,1000);%partation of label the index
i1=single(2);
labellist(1)=0;
for i=2:r-1%label �������ĳһ��label�����һ��Ԫ�ص�����U/Q�Ļ���
    if (isempty(find(data(i+1,c)==data(i,c),1)))
        labellist(i1)=i;
        i1=i1+1;%sequence of labellist
    end
end

labellist(i1)=r;
labellist=labellist(:,1:i1);
%[~, labelnum] = size(labellist);
% M=clus(labelnum-1);%�����ÿһ������
%���Ա�׼��
v = std(data,0);
M=zeros(r,c-1);
for i=1:c-1
    M(:,i)=data(:,i)/v(i);
end
clearvars data v i1
disp("==========================");
U=cell(1,c-1);%������ϵ����
L=zeros(1,r);%��¼class�ľ��� ���㵽������Ϊ1
C=linspace(1,c-1,c-1);%��������
R=zeros(1,c-1);%
DD=[C;zeros(1,c-1)];%������

 for i=1:c-1
     sum=0;%
        for x=1:r
            for y=1:r
                temp=(a-a*abs(M(x,i)-M(y,i))+b*min(1-M(x,i),1-M(y,i)))./(a-(a-1)*abs(M(x,i)-M(y,i))+b*min(1-M(x,i),1-M(y,i)));
                U{1,i}(x,y)=max(min(temp,1),0);
            %tw-relation
            end
            T=U{1,i}(x,:);%
            maxv=0;
            for y=1:length(labellist)-1%U/Q
                L=zeros(1,r);
                L(:,(labellist(1,y)+1):labellist(1,y+1))=ones(1,labellist(1,y+1)-labellist(1,y));
%                 disp(min(1-T+L,1));%implicator
                if (maxv<min(min(1-T+L,1)))%inf and posȡsup
                    maxv=min(min(1-T+L,1));
               end
            end
             sum=sum+maxv;
        end
        DD(2,i)=sum/r;
 end%��һ����������
 clear M
    CM=ones(r);%t norm matric
    [~,column]=size(C);
    for i=C
        CM=max(CM+U{1,i}-1,0);
    end
    sum=0;
    for x=1:r
         T=CM(x,:);
          maxv=0;
          for y=1:length(labellist)-1%U/Q
                L=zeros(1,r);
                L(:,(labellist(1,y)+1):labellist(1,y+1))=ones(1,labellist(1,y+1)-labellist(1,y));
                %disp(min(1-T+L,1));
                if (maxv<min(min(1-T+L,1)))
                    maxv=min(min(1-T+L,1));
                end
           end
           sum=sum+maxv;
    end
 totalDD=sum./r;    
 [index]=find(DD(2,:)==max(DD(2,:)));
 maxgamma=max(DD(2,:));
 R(index)=index;
 str=['feature=' num2str(index)];
 str2=['gamma=' num2str(maxgamma)];
  disp(str);
  disp(str2);
sum=0;
while(true)
    CM=ones(r);%t norm matric
    C=R(R~=0);
    [~,column]=size(C);
     DD=zeros(2,c-1,'single');
    for i=C
        CM=max(CM+U{1,i}-1,0);
    end
    for i=1:c-1%ÿһ������
        if(R(i)==0)
           
            CM2=max(CM+U{1,i}-1,0);%t-norm
            
            for x=1:r
               T=CM2(x,:);
               maxv=0;
               for y=1:length(labellist)-1%U/Q
                    L=zeros(1,r);
                    L(:,(labellist(1,y)+1):labellist(1,y+1))=ones(1,labellist(1,y+1)-labellist(1,y));
                    %disp(min(1-T+L,1));
                    if (maxv<min(min(1-T+L,1)))
                        maxv=min(min(1-T+L,1));
                    end
               end
               sum=sum+maxv;
            end
         DD(1,i)=i;
        DD(2,i)=sum/r;
        if (sum/r >= 1)
            break
        end
        sum=0;
        end

    end
     if (max(DD(2,:))==maxgamma)||(maxgamma==1)||(nnz(R)==c-1)
            break
     elseif (max(DD(2,:))>=maxgamma)
            [index]=find(DD(2,:)==max(DD(2,:)));
            R(index)=index;
            maxgamma=max(DD(2,:));
            str=['feature=' num2str(index)];
            str2=['3gamma=' num2str(maxgamma)];
            disp(str);
            disp(str2);
     end
end
toc;
R(R==0)=[];
disp(a);
disp(length(R));
disp(R);



    


        
        

