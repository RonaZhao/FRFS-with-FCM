% Programmed by Javad Rahimipour Anaraki on 23/05/18
% Ph.D. Candidate
% Department of Computer Science
% Memorial University of Newfoundland
% jra066 [AT] mun [DOT] ca | www.cs.mun.ca/~jra066

% This is an implementation of Threshold Quick Reduct algorithm
% Input: dataset name in the same path
% Output: is a subset of selected features
% More info: https://ieeexplore.ieee.org/document/5955425/
clear
global data;
%%================================Data=====================================
strArray = cell(4,1);
strArray{1} = 'Data2/yeast.csv';


for i = 1:1
    data = importdata(strArray{i});
    %========================Preparing the Data============================
    [r, c] = size(data);
    data = sortrows(data, c);
    %===================Calculating overall dependency=====================
    allF = ones(1, c-2);
    totalDD = dependency(allF);
    %==============================Variables===============================
    s = 1;
    true = 1;
    fDD = 0;
    supMat = zeros(r, 1);
    maxDD = zeros(2, c-2);
    md = zeros(r, 3);
    md(:, 1) = data(:, 1);
    maxF = [];
    v = std(data, 0);
    cData = data;
    %====================Calculating Indiscernibilities====================
    pData = [cData(:, 1), cData(:, end)];
    [indI, cls] = IND(pData);
    tmp = zeros(1, length(cls));
    %=========Finding other features referring to dependency dergree=======
    tic
    while(true)
        for f = 2:c-1
            if (~isempty(find(maxF == f, 1)))
                continue;
            end
            moRmF = zeros(1, length(maxF));
            maxF = [maxF, f];
            for x = 1:r
                for nCls = 1:length(cls)
                    part = indI(s:cls(nCls));
                    for y = 1:r
                        
                        if (~isempty(find(data(y,1) == part(1,:), 1)))
                            md(y, 2) = 1;
                        end
                        
                        for h = 1:length(maxF)
                            Mf = maxF(h);
                            fTerm1 = (data(y,Mf)-(data(x,Mf)-v(Mf))) / (data(x,Mf)-(data(x,Mf)-v(Mf)));
                            fTerm2 = ((data(x,Mf)+v(Mf))-data(y,Mf)) / ((data(x,Mf)+v(Mf))-data(x,Mf));
                            moRmF(h) = max(min(fTerm1, fTerm2), 0);
                        end
                        if (length(moRmF) > 1)
                            out = max(moRmF(1)+moRmF(2)-1, 0);
                            for nMoRmF = 3:length(moRmF)
                                out = max(out+moRmF(nMoRmF)-1, 0);
                            end
                        else
                            out = min(moRmF);
                        end
                        
                        md(y, 3) = min(1-out+md(y,2), 1);
                    end
                    tmp(nCls) = min(md(:, 3));
                    md = zeros(r, 3);
                    md(:, 1) = data(:, 1);
                    s = cls(nCls) + 1;
                end
                supMat(x) = max(tmp);
                s = 1;
            end
            
            DD = sum(supMat);
            maxDD(1, f-1) = f;
            maxDD(2, f-1) = DD/r;
            s = 1;
            tmp = zeros(1, nCls);
            md = zeros(r, 3);
            md(:, 1) = data(:, 1);          
            maxF = maxF(:, 1:end-1);
        end
        
        if ((~isempty(maxDD)) && (max(fDD) < max(maxDD(2,:))) && abs((totalDD-max(maxDD(2, :)))) * r >= 1)
            [~, mxF] = max(maxDD(2, :));
            maxF = [maxF, maxDD(1,mxF)];
            fDD = [fDD, maxDD(2,mxF)];
            disp(fDD) %Uncomment this to see the dependency values
        else
            true=0;
        end
        
        if (length(fDD) == (c-2))||(max(fDD) == 1)
            true = 0;
        end
        maxDD = [];
    end
    
    %================================Output================================
    disp(strArray{i});
    disp(['Selected features: ', num2str(maxF-1)]);
    disp('==============================================================');
    toc
end