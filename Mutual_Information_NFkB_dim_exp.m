% Mutual Information Calculation (binless method)
 function Mutual_information = Mutual_Information_NFkB_dim_exp(ts,c1)
%%
        knn = 10; % KNN Value
%% entropy Calculation (binless method)
        cond_entropy = 0;
        total_entropy = 0;
    for ns = 1: ts
        X = c1{ns}'*10^3; % convert the NFkB concentration from 1uM to 10^3 pM unit
% conditional entropy Calculation (binless method)
        nobs = size(X,1);
        dimen = size(X,2); % vector dimension
        Vd = (pi^(dimen/2))/gamma((dimen/2)+1);
        [~,d] = knnsearch(X,X,'k',knn);
        %
        factor1 = log2(knn./(nobs*Vd.*(d(:,knn).^dimen)));
        cond_entropy = cond_entropy + (-(sum(factor1)/(nobs*ts)));
% Total entropy Calculation (binless method)
            cc = cell(ts,1);
        for nn = 1: ts
            %
            Y = c1{nn}'*10^3; % convert the NFkB concentration from 1uM to 10^3 pM unit
            [~,dt] = knnsearch(Y,X,'k',knn);
            cc{nn} = dt;
        end
        %
        total_entropy1 = 0;
        for nr = 1: nobs
            sumw = 0;
            tw = 1/ts;
        %
        for nw = 1: ts
            df = cc{nw};
            sumw = sumw + ((knn*tw)/(nobs*Vd*(df(nr,knn)^dimen)));
        end
        factor3 = log2(sumw);
        total_entropy1 = total_entropy1 + factor3;
        end
         factor4 = -(total_entropy1/(nobs*ts));
         total_entropy = total_entropy + factor4;
    end
%% Information transfer (I) between an input(S) and output (R) using this formula
%  I(R;S) = H(R) - H(R/S)
 Mutual_information = total_entropy - cond_entropy;
 end
