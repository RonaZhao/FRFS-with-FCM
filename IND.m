% Programmed by Javad Rahimipour Anaraki on 23/05/18
% Ph.D. Candidate
% Department of Computer Science
% Memorial University of Newfoundland
% jra066 [AT] mun [DOT] ca | www.cs.mun.ca/~jra066

% This is an implementation of to find Indiscernibile objects
% Input: dataset
% Output: decerbinle index of samples
% More info: https://ieeexplore.ieee.org/document/5955425/

function [ind, cls] = IND(data)
    ind = [];
    cls = [];
    l = 1;
    [r, c]= size(data);
    while (l<r+1)
        n=1;
        if (length(ind) == length(data))
            break;
        end
        while ((n <= length(ind)) && (l < length(data)))
            if (ind(1, n) == data(l, 1))
                l = l + 1;
                n = 1;
            else
                n = n + 1;
            end
        end
        if (l < r)
            ind = [ind, data(l, 1)];
        end
        for m = (l+1):r
            if (data(l, c) == data(m, c))
                ind = [ind, data(m, 1)];
            end
        end
        cls = [cls; length(ind)];
        l = l + 1;
    end

ind = ind(:, 1:end);