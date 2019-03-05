function out = clus(labelnum)
global data;
[r, c] = size(data);
% cdata=zscore(data);
% data = sortrows(data, c);
% disp(data);
% cdata = data(:,1) ;
% [center,U,obj_fcn] = FCMCluster(cdata,2); 
% M = U(1,:);
% M = M.';%.‘是转置的意思
% [r, c] = size(M);
% fprintf("r=%d,c=%d",r,c);
M = zeros(r,c-1);
X=linspace(1,r,r);
options=[2;100;1e-5;0];
for i=1:c-1
    cdata=data(:,i);
%     cdata=[X.',cdata];
    eva=evalclusters(cdata,'kmeans','DaviesBouldin', 'KList', [1:labelnum]);
%     disp(eva.OptimalK);
    [center,U,~] = fcm(cdata,eva.OptimalK,options); 
%     [seeds] = kpp(cdata,eva.OptimalK);
%     disp(center);
%     disp(seeds);
%     m=center==min(center(:,1));
    MT = max(U,[],1);
    MT = MT.';
    [r, c] = size(MT);
%     fprintf("r=%d,c=%d\n",r,c);
    M(:,i)= MT;
%     C=center;
%     figure(1);
%     subplot(3,5,i)
%     plot(MT,'*r');
%     box off
%     figure(1);
%     subplot(3,5,i);
%     plot(X.',cdata(:,1),'.y');
%     hold on
%     plot(center,'*b');
%     hold off
%     box off
end
out=M;