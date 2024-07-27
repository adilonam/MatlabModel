    function [gMM,dMM_dTheta,dDeltadThetanorm] = MicroMoments(Moments,Schools,Consumers,Distance,Cweights,Estimation,Set,Params)
        
        MM = [];
        dMM_dTheta = [];
        dDeltadThetanorm = [];
        
        for m = Set.marketslist'
            
            for y = Set.years

                DeltaM = Schools.Delta(Schools.MarketId==m & Schools.Year==y);
                SchoolsM = Schools(Schools.MarketId==m & Schools.Year==y,:);
                DistanceM = Distance(Schools.MarketId==m & Schools.Year==y,Consumers.MarketId==m);
                CweightsMAll = squeeze(Cweights.All(Consumers.MarketId==m,2,:));
                CweightsMTypes = squeeze(Cweights.Types(Consumers.MarketId==m,2,:));
                Schools.Share(Schools.MarketId==m & Schools.Year==y) = RC_shares(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
                SchoolsM = Schools(Schools.MarketId==m & Schools.Year==y,:);
                
                if nargout>=2
                    [~,~,~,~,~,MMm,dMM_dThetam,dDeltadThetanormm] =RC_shares(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
                else
                    [~,~,~,~,~,MMm] =RC_shares(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
                end
                
                MM = [MM;MMm];
                if nargout>=2
                    dMM_dTheta = [dMM_dTheta;dMM_dThetam];
                    dDeltadThetanorm = [dDeltadThetanorm;dDeltadThetanormm];
                end
                
            end
        end
        
        gMM=MM-table2array(Moments.MM(:,5));
        gMM(Moments.MM.SanityCheck==0) = [];
        if nargout>=2
            dMM_dTheta(Moments.MM.SanityCheck==0,:) = [];
        end
        
    end