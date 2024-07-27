function deltaNew=SolveNewton(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params)

g= @(x)(rc_newton(x));
ops=optimset('Display','off','Jacobian','on','DerivativeCheck','off','FinDiffType','central','maxIter',1000,'TolFun',1e-18);
[deltaNew,~,~,~] = fsolve(g,DeltaM,ops);


if (deltaNew(SchoolsM.NormId==1)~=0)
    %fprintf('Not normed correctly');
    deltaNew=deltaNew-deltaNew(SchoolsM.NormId==1);
end

    function [f,JacLnS]=rc_newton(del)
    %function [f]=rc_newton(del)
        [S,~,JacS] = RC_shares(del,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
        %[S] = RC_shares(del,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Theta2);
        f=-log(SchoolsM.Share)+log(S);
        JacLnS= bsxfun(@rdivide,JacS,S);
    end

end