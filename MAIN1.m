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
in='Data1/hillvalley.csv';
a=0.01;%%%%%根据实际情况更改
b=0;
R1=0;
data = csvread(in,headerlinesIn);
% data = csvread(in,1,R1,[1,R1,4500,3]);%从文件导入数据
[r,c]=size(data);%row and column of dataset
data=sortrows(data,c);%ascending order of the label

%====================================
labellist=zeros(1,1000);%partation of label the index
i1=2;
labellist(1)=0;
for i=2:r-1%label 存的数是某一类label的最后一个元素的索引U/Q的划分
    if (isempty(find(data(i+1,c)==data(i,c),1)))
        labellist(i1)=i;
        i1=i1+1;%sequence of labellist
    end
end

labellist(i1)=r;
labellist=labellist(:,1:i1);
[~, labelnum] = size(labellist);
M=clus(labelnum-1);%聚类对每一个特征
disp("==========================");
U=cell(1,c-1);%特征关系矩阵
L=zeros(1,r);%记录class的矩阵 计算到该类则为1
C=linspace(1,c-1,c-1);%特征集合
R=zeros(1,c-1);%
DD=[C;zeros(1,c-1)];%依赖度

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
            for j=1:length(labellist)-1%U/Q
                L=zeros(1,r);
                L(:,(labellist(1,j)+1):labellist(1,j+1))=ones(1,labellist(1,j+1)-labellist(1,j));
%                 disp(min(1-T+L,1));%implicator
                if (maxv<min(min(1-T+L,1)))%inf and pos取sup
                    maxv=min(min(1-T+L,1));
               end
            end
             sum=sum+maxv;
        end
        DD(2,i)=sum/r;
 end%第一遍特征运算
 [index]=find(DD(2,:)==max(DD(2,:)));
 maxgamma=max(DD(2,:));
 R(index)=index;
%  str=['feature=' num2str(index)];
%  str2=['gamma=' num2str(maxgamma)];
%   disp(str);
%   disp(str2);
    
while(true)
     if (maxgamma==1)||(nnz(R)==0)
                break
     end
    CM=ones(r);%t norm matric
    A=R(R~=0);
    [~,column]=size(A);
     DD=zeros(2,c-1);
    for sef=A
        CM=max(CM+U{1,sef}-1,0);
    end
    for i=1:c-1%每一个特征
        if(R(i)==0)
           
            CM2=max(CM+U{1,i}-1,0);%t-norm
            sum=0;
            for x=1:r
               T=CM2(x,:);
               maxv=0;
               for j=1:length(labellist)-1%U/Q
                    L=zeros(1,r);
                    L(:,(labellist(1,j)+1):labellist(1,j+1))=ones(1,labellist(1,j+1)-labellist(1,j));
                    %disp(min(1-T+L,1));
                    if (maxv<min(min(1-T+L,1)))
                        maxv=min(min(1-T+L,1));
                    end
               end
               sum=sum+maxv;
            end
         DD(1,i)=i;
        DD(2,i)=sum/r;
        end
        if (sum/r >= 1)
            break
        end
    end
    if (max(DD(2,:))-maxgamma<(1/r))||(maxgamma==1)||(nnz(R)==c-1)
            break
    elseif (max(DD(2,:))>=maxgamma)
            [index]=find(DD(2,:)==max(DD(2,:)));
            R(index)=index;
            maxgamma=max(DD(2,:));
%             str=['feature=' num2str(index)];
%             str2=['3gamma=' num2str(maxgamma)];
%             disp(str);
%             disp(str2);
    end
end
toc;
R(R==0)=[];
disp(a);
disp(length(R));
disp(R);



    


        
        

