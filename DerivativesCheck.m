function diff = DerivativesCheck(Schools,Consumers,Distance,Cweights,Moments,Estimation,Set,Theta2,m,y,Supply)

Params = GetParams(Estimation,Set,Theta2);

DeltaM = Schools.Delta(Schools.MarketId==m & Schools.Year==y);
SchoolsM = Schools(Schools.MarketId==m & Schools.Year==y,:);
DistanceM = Distance(Schools.MarketId==m & Schools.Year==y,Consumers.MarketId==m);
CweightsMAll = squeeze(Cweights.All(Consumers.MarketId==m,2,:));
CweightsMTypes = squeeze(Cweights.Types(Consumers.MarketId==m,2,:));
Schools.Share(Schools.MarketId==m & Schools.Year==y) = RC_shares(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
SchoolsM = Schools(Schools.MarketId==m & Schools.Year==y,:);

[SM, ~,share_ijvM,dSdDelta,dSdTheta,MM,DMM_dTheta,dDeltadThetanorm] =RC_shares(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);

dS_dDeltaNum=zeros(height(SchoolsM));

for i=1:height(SchoolsM)
    
    DeltaA=DeltaM;
    DeltaB=DeltaM;
    
    DeltaA(i)=DeltaA(i)-Set.step/2;
    DeltaB(i)=DeltaB(i)+Set.step/2;
    
    SharesA = RC_shares(DeltaA,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
    SharesB = RC_shares(DeltaB,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
    
    dS_dDeltaNum(:,i)=((SharesB)-(SharesA))/Set.step;

end

diff(1,1)=max(max(abs(dS_dDeltaNum-dSdDelta)));
diff(2,1)=max(max(abs((dS_dDeltaNum-dSdDelta)./dSdDelta)));

fprintf('\n 1. dS/dDelta: \n Max abs diff numerical gradient and analytical version = %10.3e  \n',diff(1,1))
fprintf(' Max perc diff numerical gradient and analytical version =%10.3e  \n',diff(2,1))

%%%%%  2. Check drivatives wrt Theta %%%%%%%%%%%%%%%%%%%

dS_dThetaNum=zeros(height(SchoolsM),size(Theta2,1));

for i=1:size(Theta2,1)
    
    Theta2A=Theta2;
    Theta2B=Theta2;
    
    Theta2A(i)=Theta2A(i)-Set.step/2;
    Theta2B(i)=Theta2B(i)+Set.step/2;
    
    ParamsA = GetParams(Estimation,Set,Theta2A);
    ParamsB = GetParams(Estimation,Set,Theta2B);
    
    SharesA = RC_shares(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsA);
    SharesB = RC_shares(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsB);
    
    dS_dThetaNum(:,i)=((SharesB)-(SharesA))/Set.step;

end

diff(1,2)=max(max(abs(dS_dThetaNum-dSdTheta)));
diff(2,2)=max(max(abs((dS_dThetaNum-dSdTheta)./dSdTheta)));

fprintf('\n 2. dS/dTheta: \n Max abs diff numerical gradient and analytical version = %10.3e  \n',diff(1,2))
fprintf(' Max perc diff numerical gradient and analytical version =%10.3e \n',diff(2,2))

%%%%  3. Check drivatives wrt quality %%%%%%%%%%%%%%%%%%%

[~,~,dSdq,dSdp] = SupplyMarks(SM,share_ijvM,SchoolsM,CweightsMAll,CweightsMTypes,Estimation,Set,Theta2,'Model1');

dS_dQNUM=zeros(height(SchoolsM),1);

for i=1:height(SchoolsM)
         
    SchoolsA=SchoolsM;
    SchoolsB=SchoolsM;
    
    SchoolsA.Mu(i)=SchoolsM.Mu(i)-Set.step/2;
    SchoolsB.Mu(i)=SchoolsM.Mu(i)+Set.step/2;
    
    S_A =RC_shares(DeltaM,SchoolsA,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
    S_B =RC_shares(DeltaM,SchoolsB,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
    
    dS_dQNUM(i,1)=(S_B(i)-S_A(i))/Set.step;
    
end

diff(1,3)=max(max(abs(dS_dQNUM-diag(dSdq))));
diff(2,3)=max(max(abs((dS_dQNUM-diag(dSdq))./diag(dSdq))));

if Supply == 0
    fprintf('\n 3. dS/dQuality: \n Max abs diff numerical gradient and analytical version = %10.3e  \n',diff(1,3))
else
    fprintf('\n 3.a. dS/dQuality: \n Max abs diff numerical gradient and analytical version = %10.3e  \n',diff(1,3))
end
fprintf(' Max perc diff numerical gradient and analytical version =%10.3e \n',diff(2,3))

%%%%  3.b. Check drivatives wrt quality %%%%%%%%%%%%%%%%%%%

if Supply == 1
    dS_dQNUM=zeros(height(SchoolsM));
    
    for i=1:height(SchoolsM)
        
        SchoolsA=SchoolsM;
        SchoolsB=SchoolsM;
        
        SchoolsA.Mu(i)=SchoolsM.Mu(i)-Set.step/2;
        SchoolsB.Mu(i)=SchoolsM.Mu(i)+Set.step/2;
        
        S_A =RC_shares(DeltaM,SchoolsA,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
        S_B =RC_shares(DeltaM,SchoolsB,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
        
        dS_dQNUM(i,:)=(S_B-S_A)'/Set.step;
        
    end
    
    diff(1,3)=max(max(abs(dS_dQNUM-dSdq)));
    diff(2,3)=max(max(abs((dS_dQNUM-(dSdq))./(dSdq))));
    
    fprintf('\n 3.b. dS/dQuality: \n Max abs diff numerical gradient and analytical version = %10.3e  \n',diff(1,3))
    fprintf(' Max perc diff numerical gradient and analytical version =%10.3e \n',diff(2,3))
end

%%%%  4. Check drivatives wrt price %%%%%%%%%%%%%%%%%%%

Set.step=10^-6;dS_dPNUM=zeros(height(SchoolsM),1);

for i=1:height(SchoolsM)
     
    SchoolsA=SchoolsM;
    SchoolsB=SchoolsM;
    
    SchoolsA.Price(i)=SchoolsM.Price(i)-Set.step/2;
    SchoolsB.Price(i)=SchoolsM.Price(i)+Set.step/2;
    
    S_A =RC_shares(DeltaM,SchoolsA,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
    S_B =RC_shares(DeltaM,SchoolsB,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
    
    dS_dPNUM(i,1)=(S_B(i)-S_A(i))/Set.step;

end

diff(1,4)=max(max(abs(dS_dPNUM-diag(dSdp))));
diff(2,4)=max(max(abs((dS_dPNUM-diag(dSdp))./diag(dSdp))));

if Supply == 0
    fprintf('\n 4. dS/dPrice: \n Max abs diff numerical gradient and analytical version = %10.3e  \n',diff(1,4))
else
    fprintf('\n 4.a. dS/dPrice: \n Max abs diff numerical gradient and analytical version = %10.3e  \n',diff(1,4))
end
fprintf(' Max perc diff numerical gradient and analytical version =%10.3e \n ',diff(2,4))


%%%%  4.b Check drivatives wrt price %%%%%%%%%%%%%%%%%%%

if Supply == 1
    Set.step=10^-6;dS_dPNUM=zeros(height(SchoolsM));
    
    for i=1:height(SchoolsM)
        
        SchoolsA=SchoolsM;
        SchoolsB=SchoolsM;
        
        SchoolsA.Price(i)=SchoolsM.Price(i)-Set.step/2;
        SchoolsB.Price(i)=SchoolsM.Price(i)+Set.step/2;
        
        S_A =RC_shares(DeltaM,SchoolsA,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
        S_B =RC_shares(DeltaM,SchoolsB,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
        
        dS_dPNUM(i,:)=(S_B-S_A)/Set.step;
        
    end
    
    diff(1,4)=max(max(abs(dS_dPNUM-(dSdp))));
    diff(2,4)=max(max(abs((dS_dPNUM-(dSdp))./(dSdp))));
    
    fprintf('\n 4.b. dS/dPrice: \n Max abs diff numerical gradient and analytical version = %10.3e  \n',diff(1,4))
    fprintf(' Max perc diff numerical gradient and analytical version =%10.3e \n ',diff(2,4))
end

%%%%%  5. Check derivatives of Delta wrt theta %%%%%%%%%%%%%%%%%%%


dDelta_dThetaNum=zeros(height(SchoolsM),size(Theta2,1));

for i=1:size(Theta2,1)
    
    Theta2A=Theta2;
    Theta2B=Theta2;
    
    Theta2A(i)=Theta2A(i)-Set.step/2;
    Theta2B(i)=Theta2B(i)+Set.step/2;
    
    ParamsA = GetParams(Estimation,Set,Theta2A);
    ParamsB = GetParams(Estimation,Set,Theta2B);
    
    DeltaM_A=SolveSquarem2(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsA);
    DeltaM_B=SolveSquarem2(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsB);   
    
    dDelta_dThetaNum(:,i)=((DeltaM_B)-(DeltaM_A))/Set.step;

end

diff(1,5)=max(max(abs(dDelta_dThetaNum-dDeltadThetanorm)));
diff(2,5)=max(max(abs((dDelta_dThetaNum-dDeltadThetanorm)./dDeltadThetanorm)));

fprintf('\n 5. dDelta/dTheta: \n Max abs diff numerical and analytical gradient for the Moments = %10.3e  \n',diff(1,5))
fprintf(' Max perc diff numerical and analytical gradient for the Moments =%10.3e \n',diff(2,5))

%%%%%  6. Check derivatives of the moments wrt theta %%%%%%%%%%%%%%%%%%%

dMo_dThetaNum=zeros(length(MM),size(Theta2,1));

for i=1:size(Theta2,1)
    
    Theta2A=Theta2;
    Theta2B=Theta2;
    
    Theta2A(i)=Theta2A(i)-Set.step/2;
    Theta2B(i)=Theta2B(i)+Set.step/2;
    
    ParamsA = GetParams(Estimation,Set,Theta2A);
    ParamsB = GetParams(Estimation,Set,Theta2B);
    
    DeltaM_A=SolveSquarem2(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsA);
    [~,~,~,~,~,MM_A] = RC_shares(DeltaM_A,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsA);
    
    DeltaM_B=SolveSquarem2(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsB);   
    [~,~,~,~,~,MM_B] = RC_shares(DeltaM_B,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsB);
    
    dMo_dThetaNum(:,i)=((MM_B)-(MM_A))/Set.step;

end

diff(1,6)=max(max(abs(dMo_dThetaNum-DMM_dTheta)));
diff(2,6)=max(max(abs((dMo_dThetaNum-DMM_dTheta)./DMM_dTheta)));

fprintf('\n 6. dM/dTheta: \n Max abs diff numerical and analytical gradient for the Moments = %10.3e  \n',diff(1,6))
fprintf(' Max perc diff numerical and analytical gradient for the Moments =%10.3e \n',diff(2,6))

%%%%%  7. Check derivatives of IV component of the Obj function wrt theta %%%%%%%%%%%%%%%%%%%
tic;
[~,dQ,~,dQ_IV,~,dQ_MM,~,dQ_RDM]=GMM_ObjEval_NFP(Schools,Consumers,Distance,Cweights,Moments,Estimation,Set,Theta2);
Time=toc;

fprintf('\n Time to evaluate the obj fun + gradient = %10.2 seconds  \n',Time(1,1))


Y = Schools.Delta;
X = [Schools.Mu Schools.Price Schools.XX_exogenous Estimation.ChainFE Estimation.MarketYearFE];
Z = [Schools.Instruments Schools.XX_exogenous Estimation.ChainFE Estimation.MarketYearFE];

[~,resid,W]=ivregression(Y,X,Z);

JIVnum = zeros(size(resid,1),size(Theta2,1));

for i=1:size(Theta2,1)
    
    Theta2A=Theta2;
    Theta2B=Theta2;
    
    Theta2A(i)=Theta2A(i)-Set.step/2;
    Theta2B(i)=Theta2B(i)+Set.step/2;
               
    DeltaM_A=SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2A);
    [~,residA]=ivregression(DeltaM_A,X,Z);
    
    DeltaM_B=SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2B);
    [~,residB]=ivregression(DeltaM_B,X,Z);
    
    JIVnum(:,i)=(residB-residA)/Set.step;
    
end

dQ_IV_Num=-2*(JIVnum'*Z)*W*(resid'*Z)';

diff(1,7)=max(max(abs(dQ_IV_Num-dQ_IV)));
diff(2,7)=max(max(abs((dQ_IV_Num-dQ_IV)./dQ_IV)));

fprintf('\n 7. dQ_IV/dTheta: \n Max abs diff for numerical and analytical version for IV gradient= %10.3e  \n',diff(1,7))
fprintf(' Max perc diff for numerical and analytical version for IV gradient= %10.3e  \n',diff(2,7))


%%%%%  8. Check derivatives of MM component of the Obj function wrt theta %%%%%%%%%%%%%%%%%%%


gMM = MicroMoments(Moments,Schools,Consumers,Distance,Cweights,Estimation,Set,Params);
J_MM_num = zeros(size(gMM,1),size(Theta2,1));

for i=1:size(Theta2,1)
    
    Theta2A=Theta2;
    Theta2B=Theta2;
    
    Theta2A(i)=Theta2A(i)-Set.step/2;
    Theta2B(i)=Theta2B(i)+Set.step/2;
    
    ParamsA = GetParams(Estimation,Set,Theta2A);
    ParamsB = GetParams(Estimation,Set,Theta2B);
    
    SchoolsA=Schools;
               
    SchoolsA.Delta=SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2A);
    gMMA = MicroMoments(Moments,SchoolsA,Consumers,Distance,Cweights,Estimation,Set,ParamsA);

    SchoolsB=Schools;
                  
    SchoolsB.Delta=SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2B);
    gMMB = MicroMoments(Moments,SchoolsB,Consumers,Distance,Cweights,Estimation,Set,ParamsB);

    
    J_MM_num(:,i)=(gMMB-gMMA)/Set.step;
    
end

dQ_MM_num = 2*J_MM_num'*Moments.WMM*gMM;

diff(1,8)=max(max(abs(dQ_MM-dQ_MM_num)));
diff(2,8)=max(max(abs((dQ_MM-dQ_MM_num)./dQ_MM_num)));

fprintf('\n 8. dQ_MM/dTheta: \n Max abs diff for numerical and analytical version for MM gradient= %10.3e  \n',diff(1,8))
fprintf(' Max perc diff for numerical and analytical version for MM gradient= %10.3e \n',diff(2,8))

%%%%%  9. Check derivatives of RD component of the Obj function wrt theta %%%%%%%%%%%%%%%%%%%
[~,dQ,~,dQ_IV,~,dQ_MM,~,dQ_RDM]=GMM_ObjEval_NFP(Schools,Consumers,Distance,Cweights,Moments,Estimation,Set,Theta2);

gRDM = RDM_Simulation(Schools,Consumers,Distance,Moments,Estimation,Set,Params);
J_RDM_num = zeros(size(gRDM,1),size(Theta2,1));

for i=1:size(Theta2,1)
    
    Theta2A=Theta2;
    Theta2B=Theta2;
    
    Theta2A(i)=Theta2A(i)-Set.step/2;
    Theta2B(i)=Theta2B(i)+Set.step/2;
    
    ParamsA = GetParams(Estimation,Set,Theta2A);
    ParamsB = GetParams(Estimation,Set,Theta2B);
    
    SchoolsA=Schools;
               
    SchoolsA.Delta=SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2A);
%    gRDMA = RDM_Simulation(SchoolsA,Consumers,Distance,Moments,Estimation,Set,Params);
    gRDMA = RDM_Simulation(SchoolsA,Consumers,Distance,Moments,Estimation,Set,ParamsA);

    SchoolsB=Schools;
                  
    SchoolsB.Delta=SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2B);
%    gRDMB = RDM_Simulation(SchoolsB,Consumers,Distance,Moments,Estimation,Set,Params);
    gRDMB = RDM_Simulation(SchoolsB,Consumers,Distance,Moments,Estimation,Set,ParamsB);

    
    J_RDM_num(:,i)=(gRDMB-gRDMA)/Set.step;
    
end

dQ_RDM_num = 2*J_RDM_num'*Moments.WRDM*gRDM;

diff(1,9)=max(max(abs(dQ_RDM-dQ_RDM_num)));
diff(2,9)=max(max(abs((dQ_RDM-dQ_RDM_num)./dQ_RDM_num)));

fprintf('\n 9. dQ_RDM/dTheta: \n Max abs diff for numerical and analytical version for RDM gradient= %10.3e  \n',diff(1,9))
fprintf(' Max perc diff for numerical and analytical version for RDM gradient= %10.3e \n',diff(2,9))
 
%%%%%  9.b. Check derivatives of RD component of the Obj function wrt theta %%%%%%%%%%%%%%%%%%%

[~,~,~,dSdTheta] = RDM_Simulation(Schools,Consumers,Distance,Moments,Estimation,Set,Params,dDeltadThetanorm);
J_RDM_num = zeros(size(gRDM,1),size(Theta2,1));

dSimY{17,1}=[];

%for i=1:size(Theta2,1)
for i=1
    
    Theta2A=Theta2;
    Theta2B=Theta2;
    
    Theta2A(i)=Theta2A(i)-Set.step/2;
    Theta2B(i)=Theta2B(i)+Set.step/2;
    
    ParamsA = GetParams(Estimation,Set,Theta2A);
    ParamsB = GetParams(Estimation,Set,Theta2B);
    
    SchoolsA=Schools;
               
    
    %SchoolsA.Delta=SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2A);
    [~,~,SimYA] = RDM_Simulation(SchoolsA,Consumers,Distance,Moments,Estimation,Set,ParamsA,dDeltadThetanorm);

    SchoolsB=Schools;
                  
    %SchoolsB.Delta=SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2B);
    [~,~,SimYB] = RDM_Simulation(SchoolsB,Consumers,Distance,Moments,Estimation,Set,ParamsB,dDeltadThetanorm);

    
    dSimY{i,1}=(SimYB-SimYA)/Set.step;
    
end

dQ_RDM_num = 2*J_RDM_num'*Moments.WRDM*gRDM;

diff(1,9)=max(max(abs(dQ_RDM-dQ_RDM_num)));
diff(2,9)=max(max(abs((dQ_RDM-dQ_RDM_num)./dQ_RDM_num)));

fprintf('\n 9. dQ_RDM/dTheta: \n Max abs diff for numerical and analytical version for RDM gradient= %10.3e  \n',diff(1,9))
fprintf(' Max perc diff for numerical and analytical version for RDM gradient= %10.3e \n',diff(2,9))
 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [~,dQ,~,dQ_IV,~,dQ_MM,~,dQ_RDM]=GMM_ObjEval_NFP(Schools,Consumers,Distance,Cweights,Moments,Estimation,Set,Theta2);
% 
% [gRDM] = RDM_Simulation(Schools,Consumers,Distance,Moments,Estimation,Set,Params);
% J_RDM_num = zeros(size(gRDM,1),size(Theta2,1));
% 
% for i=1:size(Theta2,1)
%     
%     Theta2A=Theta2;
%     Theta2B=Theta2;
%     
%     Theta2A(i)=Theta2A(i)-Set.step/2;
%     Theta2B(i)=Theta2B(i)+Set.step/2;
%     
%     ParamsA = GetParams(Estimation,Set,Theta2A);
%     ParamsB = GetParams(Estimation,Set,Theta2B);
%     
%     SchoolsA=Schools;
%                
%     SchoolsA.Delta=SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2A);
% %    gRDMA = RDM_Simulation(Schools,Consumers,Distance,Moments,Estimation,Set,ParamsA);
%     [gRDM,~,SSA,~] = RDM_Simulation(SchoolsA,Consumers,Distance,Moments,Estimation,Set,Params,dDeltadThetanorm);
% 
%     SchoolsB=Schools;
%                   
%     SchoolsB.Delta=SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2B);
% %    gRDMB = RDM_Simulation(Schools,Consumers,Distance,Moments,Estimation,Set,ParamsB);
%     [gRDM,~,SSB,~] = RDM_Simulation(SchoolsB,Consumers,Distance,Moments,Estimation,Set,Params,dDeltadThetanorm);
% 
%     
%     J_RDM_num(:,i)=(SSB-SSA)/Set.step;
%     
% end
% 
% dQ_RDM_num = 2*J_RDM_num'*Moments.WRDM*gRDM;
% 
% diff(1,9)=max(max(abs(dQ_RDM-dQ_RDM_num)));
% diff(2,9)=max(max(abs((dQ_RDM-dQ_RDM_num)./dQ_RDM_num)));
% 
% fprintf('\n 9. dQ_RDM/dTheta: \n Max abs diff for numerical and analytical version for RDM gradient= %10.3e  \n',diff(1,9))
% fprintf(' Max perc diff for numerical and analytical version for RDM gradient= %10.3e \n',diff(2,9))
%  
% % %%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%  10. Check derivatives of the Obj function wrt theta %%%%%%%%%%%%%%%%%%%
% 
% gRDM = RDM_Simulation(Schools,Consumers,Distance,Moments,Estimation,Set,Params);
% J_RDM_num = zeros(size(gRDM,1),size(Theta2,1));
% 
% for i=1:size(Theta2,1)
%     
%     Theta2A=Theta2;
%     Theta2B=Theta2;
%     
%     Theta2A(i)=Theta2A(i)-Set.step/2;
%     Theta2B(i)=Theta2B(i)+Set.step/2;
%     
%     ParamsA = GetParams(Estimation,Set,Theta2A);
%     ParamsB = GetParams(Estimation,Set,Theta2B);
%     
%     SchoolsA=Schools;
%                
%     SchoolsA.Delta=SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2A);
%     [~,dQA]=GMM_ObjEval_NFP(SchoolsA,Consumers,Distance,Cweights,Moments,Estimation,Set,Theta2A);
% 
%     SchoolsB=Schools;
%                   
%     SchoolsB.DeltaM_B=SolveShares(Schools.Delta,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2B);
%     [~,dQB]=GMM_ObjEval_NFP(SchoolsB,Consumers,Distance,Cweights,Moments,Estimation,Set,Theta2B);
%     
%     dQ_RDM_num(:,i)=(dQB-dQA)/Set.step;
%     i
% end
% 
% diff(1,10)=max(max(abs(dQ_RDM-dQ_RDM_num)));
% diff(2,10)=max(max(abs((dQ_RDM-dQ_RDM_num)./dQ_RDM_num)));
% 
% fprintf('\n 10. dQ_MM/dTheta: \n Max abs diff for numerical and analytical version for MM gradient= %10.3e  \n',diff(1,10))
% fprintf(' Max perc diff for numerical and analytical version for MM gradient= %10.3e \n',diff(2,10))
%   
        
%%%%%  7. Check derivatives of the moments wrt theta %%%%%%%%%%%%%%%%%%%
% 
% dMo_dThetaNum=zeros(length(MM),size(Theta2,1));
% 
% Params = GetParams(Estimation,Set,Theta2);
% 
% for i=1:size(Theta2,1)
%     
%     Theta2A=Theta2;
%     Theta2B=Theta2;
%     
%     Theta2A(i)=Theta2A(i)-Set.step/2;
%     Theta2B(i)=Theta2B(i)+Set.step/2;
%     
%     ParamsA = GetParams(Estimation,Set,Theta2A);
%     ParamsB = GetParams(Estimation,Set,Theta2B);
%     
%     DeltaM=SolveSquarem2(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
%     
%     DeltaM_A=SolveSquarem2(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsA);
%     [~,~,~,~,~,MM_A,~] = RC_shares(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsA);
%     [~,~,~,~,~,MM_AA,~] = RC_shares(DeltaM_B,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
%     [~,~,~,~,~,MM_AAA,~] = RC_shares(DeltaM_B,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsA);
%     
%     DeltaM_B=SolveSquarem2(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsB);   
%     [~,~,~,~,~,MM_B,~] = RC_shares(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsB);
%     [~,~,~,~,~,MM_BB,~] = RC_shares(DeltaM_A,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
%     [~,~,~,~,~,MM_BBB,~] = RC_shares(DeltaM_A,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,ParamsB);
%     
%     dMo_dsThetaNum(:,i)=((MM_B)-(MM_A))/Set.step;
%     dMo_dsThetaNum2(:,i)=((MM_BB)-(MM_AA))/Set.step;
%     dMo_dsThetaNum3(:,i)=((MM_BBB)-(MM_AAA))/Set.step;
% end
% 
% dif2(1)=max(max(abs(dMo_dThetaNum-DMM_dTheta)));
% dif2(2)=max(max(abs((dMo_dThetaNum-DMM_dTheta)./DMM_dTheta)));
% 
% fprintf('\n Theta: \n Max abs diff numerical and analytical gradient for the Moments = %10.3e  \n',dif2(1))
% fprintf(' Max perc diff numerical and analytical gradient for the Moments =%10.3e \n',dif2(2))
%  
end