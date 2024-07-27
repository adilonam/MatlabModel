    function [f,JacLnS]=rc_newton(del)
    %function [f]=rc_newton(del)
        [S,~,JacS] = RC_shares(del,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Theta2);
        %[S] = RC_shares(del,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Theta2);
        f=-log(SchoolsM.Share)+log(S);
        JacLnS= bsxfun(@rdivide,JacS,S);
    end