function [ Delta,  Jac ] = SolveShares(Delta0,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2)

%{

Project:        JMP
Author:         Claudia Allende
Date created:   May 23rd, 2019
Date modified:  May 24th, 2019
Description:    This code Solves shares for delta market by market. Can also return the Jacobian (d delta/ d theta)

%}
    
    Params = GetParams(Estimation,Set,Theta2);
    deltahat{length(Set.marketslist)} = [];
    Jacpart{length(Set.marketslist)} = [];
    
    k=0;
    
    for j = 1:length(Set.marketslist)
        for y= 1:length(Set.years)
            k=k+1;
            i=Set.marketslist(j);
            
            DeltaM          = Delta0(Schools.MarketId==i & Schools.Year==Set.years(y));
            SchoolsM        = Schools(Schools.MarketId==i & Schools.Year==Set.years(y),:);
            DistanceM       = Distance(Schools.MarketId==i & Schools.Year==Set.years(y),Consumers.MarketId==i);
            CweightsMAll    = squeeze(Cweights.All(Consumers.MarketId==i,2,:));
            CweightsMTypes  = squeeze(Cweights.Types(Consumers.MarketId==i,2,:));
            
            if strcmp(Set.method,'Newton')
                deltahat{k} = SolveNewton(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
            elseif strcmp(Set.method,'Squarem')
                deltahat{k} = SolveSquarem2(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
            end
            % this returns the Jacobian once after convergence (for the
            % gradient computation)
            if nargout >1
            [~,~,~,~,~,~,~,Jacpart{k}] =RC_shares(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
            end
        end
    end


Delta   =   real(cat(1,deltahat{:}));

if nargout >1
    Jac     =   real(cat(1,Jacpart{:}));
end
    
end

