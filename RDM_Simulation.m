function [gRDM,dRDM_dTheta,SimY,dSdTheta] =RDM_Simulation(Schools,Consumers,Distance,Moments,Estimation,Set,Params,dDeltadThetanorm)

gRDM = zeros(6,1);
dRDM_dTheta = zeros(6,size(Estimation.BetaQEdu,2));

SimY    = NaN(height(Moments.RDdata),6);

if nargout > 1
    
    dSdTheta{size(Estimation.BetaQEdu,2),1} = [];
    
    for i = 1:size(Estimation.BetaQEdu,2)
        dSdTheta{i,1} = zeros(size(SimY));
    end
end

for m=unique(Moments.RDdata.MarketId)'
    
    for y=unique(Moments.RDdata.Year)'
        
        % By Market, by year
        
        DeltaM = Schools.Delta(Schools.MarketId==m & Schools.Year==y);
        SchoolsM = Schools(Schools.MarketId==m & Schools.Year==y,:);
        DistanceM = Distance(Schools.MarketId==m & Schools.Year==y,Consumers.MarketId==m);
        
        % % % % % % % % % % % % % % % % % % % % % % % % Only for Model 2
        for type=1:4
    
            rdIndex = (Moments.RDdata.MarketId == m & Moments.RDdata.Year == y & Moments.RDdata.Type == type);
            
            nodes = Moments.RDdata.Node(rdIndex);
            
            if isempty(nodes)==0
                
                % Set up interaction between firm Xs and unobservables
                
                UoType  = SchoolsM.Mu*Params.betaK(type,1) + (SchoolsM.Price-0.75*(Moments.RDdata.Treat(rdIndex) == 1)'.*SchoolsM.Price)*Params.betaK(type,2) + DistanceM(:,nodes)*Params.betaK(type,3) + SchoolsM.Ze(:,type)*Params.betaK(type,4) + SchoolsM.Zy(:,type)*Params.betaK(type,5);
                Uv      = SchoolsM.Mu*Params.betai(:,1)' + SchoolsM.Ze(:,type)*Params.betai(:,2)' + SchoolsM.Zy(:,type)*Params.betai(:,3)'; % Random Coefficient Part of Utility
                Uv      = reshape(Uv,size(Uv,1),1,size(Uv,2));
                
                % get maxutility to norm vector of utility to avoid overflow
                maxuij              = max(DeltaM)+max(max(UoType))+max(max(Uv));
                num                 = exp(DeltaM+UoType+Uv-maxuij);
                share_ijv           = bsxfun(@rdivide,num,sum(num,1));
                share_ij            = sum(share_ijv.*Params.w,3);
                
                SimY(rdIndex == 1,:)  = [share_ij'*[SchoolsM.Mu SchoolsM.Price] sum(share_ij.*(SchoolsM.Price-0.75*(Moments.RDdata.Treat(rdIndex) == true)'),1)' sum(share_ij.*DistanceM(:,nodes),1)' share_ij'*[SchoolsM.Ze(:,type) SchoolsM.Zy(:,type)]];
                
                if nargout > 1
                    
                    share_ij = sum(share_ijv.*Params.w,3);
                    
                    Diag = diag(share_ij(:,1));
                    
                    for n=2:size(share_ij,2)
                        
                        Diag = cat(3,Diag,diag(share_ij(:,n)));
                        
                    end
                    
                    dSdDdDdT = squeeze(mmx('mult',(mmx('mult',mmx('mult',permute(-share_ijv,[1 3 2]),diag(squeeze(Params.w(1,1,:)))),permute(share_ijv,[3 1 2])))+Diag,dDeltadThetanorm(Schools.MarketId==m & Schools.Year==y,:)));
                    %dSdDdDdT = zeros(size(squeeze(mmx('mult',(mmx('mult',mmx('mult',permute(-share_ijv,[1 3 2]),diag(squeeze(Params.w(1,1,:)))),permute(share_ijv,[3 1 2])))+Diag,dDeltadThetanorm(Schools.MarketId==m & Schools.Year==y,:)))));
                    
                    tempB =1:size(Estimation.BetaQEdu,2);
                    
                    for i=tempB(not(ismember(tempB,Estimation.TypesTheta{type,1})))
                        
                        dSdTheta{i,1}(rdIndex,:) = [squeeze(dSdDdDdT(:,i,:))'*[SchoolsM.Mu SchoolsM.Price] sum(squeeze(dSdDdDdT(:,i,:)).*(SchoolsM.Price-0.75*(Moments.RDdata.Treat(rdIndex) == true)'),1)' sum(squeeze(dSdDdDdT(:,i,:)).*DistanceM(:,nodes),1)' squeeze(dSdDdDdT(:,i,:))'*[SchoolsM.Ze(:,type) SchoolsM.Zy(:,type)]];
                        
                    end
                    
                    for i=Estimation.TypesTheta{type,1}
                        
                        if Estimation.ThetaMask{i,2}==1
                            
                            if Estimation.ThetaMask{i,1}==0
                                
                                                              
                                dSdThetai       = squeeze(dSdDdDdT(:,i,:)) + sum(share_ijv.*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv,'tn')).*Params.w(1,1,:),3);
                                dSdTheta{i,1}(rdIndex,:) = [dSdThetai'*[SchoolsM.Mu SchoolsM.Price] sum(dSdThetai.*(SchoolsM.Price-0.75*(Moments.RDdata.Treat(rdIndex) == true)'),1)' sum(dSdThetai.*DistanceM(:,nodes),1)' dSdThetai'*[SchoolsM.Ze(:,type) SchoolsM.Zy(:,type)]];
                                
                                
                            elseif Estimation.ThetaMask{i,1}==1
                               
                                dSdThetai = squeeze(dSdDdDdT(:,i,:)) + sum(bsxfun(@times,share_ijv,reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv,'tn')).*Params.w(1,1,:),3);
                                dSdTheta{i,1}(rdIndex,:) = [dSdThetai'*[SchoolsM.Mu SchoolsM.Price] sum(dSdThetai.*(SchoolsM.Price-0.75*(Moments.RDdata.Treat(rdIndex) == true)'),1)' sum(dSdThetai.*DistanceM(:,nodes),1)' dSdThetai'*[SchoolsM.Ze(:,type) SchoolsM.Zy(:,type)]];
                                
                            end
                            
                        elseif Estimation.ThetaMask{i,2}==2
                            
                            if Estimation.ThetaMask{i,1}==0
                                
                                dSdThetai       = squeeze(dSdDdDdT(:,i,:)) +sum(share_ijv.*bsxfun(@minus,(SchoolsM.Price-0.75*(Moments.RDdata.Treat(rdIndex) == 1)'.*SchoolsM.Price),sum(bsxfun(@times,(SchoolsM.Price-0.75*(Moments.RDdata.Treat(rdIndex) == 1)'.*SchoolsM.Price),share_ijv),1)).*Params.w(1,1,:),3);
                                dSdTheta{i,1}(rdIndex,:) = [dSdThetai'*[SchoolsM.Mu SchoolsM.Price] sum(dSdThetai.*(SchoolsM.Price-0.75*(Moments.RDdata.Treat(rdIndex) == true)'),1)' sum(dSdThetai.*DistanceM(:,nodes),1)' dSdThetai'*[SchoolsM.Ze(:,type) SchoolsM.Zy(:,type)]];
                                
                            end
                            
                        elseif Estimation.ThetaMask{i,2}==3
                            
                            if Estimation.ThetaMask{i,1}==0
                                
                                dSdThetai       = squeeze(dSdDdDdT(:,i,:)) +sum(share_ijv.*bsxfun(@minus,DistanceM(:,nodes),sum(bsxfun(@times,DistanceM(:,nodes),share_ijv),1)).*Params.w(1,1,:),3);
                                dSdTheta{i,1}(rdIndex,:) = [dSdThetai'*[SchoolsM.Mu SchoolsM.Price] sum(dSdThetai.*(SchoolsM.Price-0.75*(Moments.RDdata.Treat(rdIndex) == true)'),1)' sum(dSdThetai.*DistanceM(:,nodes),1)' dSdThetai'*[SchoolsM.Ze(:,type) SchoolsM.Zy(:,type)]];
                                
                            end
                            
                        elseif Estimation.ThetaMask{i,2}==4
                            
                            if Estimation.ThetaMask{i,1}==0
                              
                                dSdThetai       = squeeze(dSdDdDdT(:,i,:)) +sum(share_ijv.*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv,'tn')).*Params.w(1,1,:),3);
                                dSdTheta{i,1}(rdIndex,:) = [dSdThetai'*[SchoolsM.Mu SchoolsM.Price] sum(dSdThetai.*(SchoolsM.Price-0.75*(Moments.RDdata.Treat(rdIndex) == true)'),1)' sum(dSdThetai.*DistanceM(:,nodes),1)' dSdThetai'*[SchoolsM.Ze(:,type) SchoolsM.Zy(:,type)]];
                                
                            elseif Estimation.ThetaMask{i,1}==1
                                
                                dSdThetai = squeeze(dSdDdDdT(:,i,:)) + sum(bsxfun(@times,share_ijv,reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv,'tn')).*Params.w(1,1,:),3);
                                dSdTheta{i,1}(rdIndex,:) = [dSdThetai'*[SchoolsM.Mu SchoolsM.Price] sum(dSdThetai.*(SchoolsM.Price-0.75*(Moments.RDdata.Treat(rdIndex) == true)'),1)' sum(dSdThetai.*DistanceM(:,nodes),1)' dSdThetai'*[SchoolsM.Ze(:,type) SchoolsM.Zy(:,type)]];
                                
                            end
                            
                        elseif Estimation.ThetaMask{i,2}==5
                            
                            if Estimation.ThetaMask{i,1}==0
                                
                                dSdThetai       = squeeze(dSdDdDdT(:,i,:)) +sum(share_ijv.*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv,'tn')).*Params.w(1,1,:),3);
                                dSdTheta{i,1}(rdIndex,:) = [dSdThetai'*[SchoolsM.Mu SchoolsM.Price] sum(dSdThetai.*(SchoolsM.Price-0.75*(Moments.RDdata.Treat(rdIndex) == true)'),1)' sum(dSdThetai.*DistanceM(:,nodes),1)' dSdThetai'*[SchoolsM.Ze(:,type) SchoolsM.Zy(:,type)]];
                                
                            elseif Estimation.ThetaMask{i,1}==1
                               
                                dSdThetai = squeeze(dSdDdDdT(:,i,:)) + sum(bsxfun(@times,share_ijv,reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv,'tn')).*Params.w(1,1,:),3);
                                dSdTheta{i,1}(rdIndex,:) = [dSdThetai'*[SchoolsM.Mu SchoolsM.Price] sum(dSdThetai.*(SchoolsM.Price-0.75*(Moments.RDdata.Treat(rdIndex) == true)'),1)' sum(dSdThetai.*DistanceM(:,nodes),1)' dSdThetai'*[SchoolsM.Ze(:,type) SchoolsM.Zy(:,type)]];
                                
                            end
                        end
                    end
                    
                end
            end
        end
        
    end
end


% Try two versions

if strcmp(Set.momentsrd,'min distance')==1
    
    for i = 1:6
        temp = ivregression(SimY(:,i),Moments.RDdata.X,Moments.RDdata.Z);
        gRDM(i,1) = Moments.RDBeta(i,1) - temp(1);
        
        if nargout > 1
            for j = 1:size(Estimation.BetaQEdu,2)
                temp = ivregression(dSdTheta{j,1}(:,i),Moments.RDdata.X,Moments.RDdata.Z);
                dRDM_dTheta(i,j) = - temp(1);
            end
        end
    end
    
elseif strcmp(Set.momentsrd,'moment condition')==1
    
    for i = 1:6
        W = (Moments.RDdata.Z'*Moments.RDdata.Z) \ eye(size(Moments.RDdata.Z,2));
        gRDM(i,1) = (Moments.RDdata.Z'*(SimY(:,i)-Moments.RDdata.X*Moments.RDBeta(1,:)'))'*W*Moments.RDdata.Z'*(SimY(:,i)-Moments.RDdata.X*Moments.RDBeta(1,:)');
        
        if nargout > 1
            for j = 1:size(Estimation.BetaQEdu,2)
                
                %dRDM_dTheta(i,j) = (Moments.RDdata.Z'*(dSdTheta{j,1}(:,i)-Moments.RDdata.X*Moments.RDBeta(1,:)'))'*W*Moments.RDdata.Z'*(dSdTheta{j,1}(:,i)-Moments.RDdata.X*Moments.RDBeta(1,:)');
                dRDM_dTheta(i,j) = 2*(Moments.RDdata.Z'*(dSdTheta{j,1}(:,i)))'*W*Moments.RDdata.Z'*(SimY(:,i)-Moments.RDdata.X*Moments.RDBeta(1,:)');
                
            end
        end
    end
end





end




% % % % % % % % % % % % % % % % % % % % % % % %





