function [Q,dQ,Q_IV,dQ_IV,Q_MM,dQ_MM,Q_RDM,dQ_RDM]=GMM_ObjEval_NFP(Schools,Consumers,Distance,Cweights,Moments,Estimation,Set,Theta2)

%%%% 1. Market Shares %%%%

Params = GetParams(Estimation,Set,Theta2);

if nargout == 1
    [Schools.Delta] = SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2);
elseif nargout > 1
    [Schools.Delta,dDelta_dTheta] = SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2);
end

%%%% 2. IV Momements %%%%

Y = Schools.Delta;
X = [Schools.Mu Schools.Price Schools.XX_exogenous Estimation.ChainFE Estimation.MarketYearFE];
Z = [Schools.Instruments Schools.XX_exogenous Estimation.ChainFE Estimation.MarketYearFE];

[~,resid,W]=ivregression(Y,X,Z);

Q_IV    = (resid'*Z-Estimation.AdjustmentObjIV)*W*(resid'*Z-Estimation.AdjustmentObjIV)';

if nargout > 1
    dQ_IV   = -2*(dDelta_dTheta'*Z)*W*(resid'*Z)';
end

%%%% 3. MicroMoments %%%%

if nargout == 1
    gMM = MicroMoments(Params);
elseif nargout > 1
    [gMM,dMM_dTheta,dDeltadThetanorm] = MicroMoments(Moments,Schools,Consumers,Distance,Cweights,Estimation,Set,Params);
end

Q_MM    = gMM'*Moments.WMM*gMM-Estimation.AdjustmentObjMM;

if nargout > 1
    dQ_MM   = 2*dMM_dTheta'*Moments.WMM*gMM;
end

%%%% 4. RD Moments %%%%
if strcmp(Set.moments,'moments 2')==1 || strcmp(Set.moments,'moments 3')==1
    if nargout == 1
        gRDM = RDM_Simulation(Schools,Consumers,Distance,Moments,Estimation,Set,Params);
    elseif nargout > 1
        [gRDM,dRDM_dTheta] = RDM_Simulation(Schools,Consumers,Distance,Moments,Estimation,Set,Params,dDeltadThetanorm);
    end
    
    Q_RDM = gRDM'*Moments.WRDM*gRDM;
    
    if nargout > 1
        dQ_RDM = 2*dRDM_dTheta'*Moments.WRDM*gRDM;
    end
end

%%%% 5. MLE Moments %%%%

if strcmp(Set.moments,'moments 3')==1
    if nargout == 1
        gRDM = RDM_Simulation(Schools,Consumers,Distance,Moments,Estimation,Set,Params);
    elseif nargout > 1
        [gRDM,dRDM_dTheta] = RDM_Simulation(Schools,Consumers,Distance,Moments,Estimation,Set,Params,dDeltadThetanorm);
    end
    
    Q_RDM = gRDM'*Moments.WRDM*gRDM;
    
    if nargout > 1
        dQ_RDM = 2*dRDM_dTheta'*Moments.WRDM*gRDM;
    end
end

%%%% 6. Put everything together %%%%

Q=Estimation.W_IV*Q_IV+Estimation.W_MM*Q_MM+Estimation.W_RDM*Q_RDM;

dQ=Estimation.W_IV*dQ_IV+Estimation.W_MM*dQ_MM+Estimation.W_RDM*dQ_RDM;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end