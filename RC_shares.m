function [S, Shares,share_ijv,dSdDelta,dSdTheta,MM,DMM_dTheta,dDeltadThetanorm] =RC_shares(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params)

% Apply Normalization by Market x Year
% if (DeltaM(SchoolsM.NormId==1)~=0)
%     %fprintf('Not normed correctly');
%     DeltaM=DeltaM-DeltaM(SchoolsM.NormId==1);
% end

% By Market, by year

S       = zeros(height(SchoolsM),1);
Shares  = NaN(height(SchoolsM),4);
share_ijv{1,4}=[];



if nargout > 3
    dSdDelta = zeros(height(SchoolsM));
    
    if nargout > 4
        dSdTheta = zeros(height(SchoolsM),size(Estimation.ThetaMask,1));
        
        if nargout > 5
            moM=[];
            %moV=[];
            
            if nargout > 6
                
                dMo_dTheta      = zeros(5*4,size(Estimation.ThetaMask,1));
                dMo_dSdDdDdT    = zeros(5*4,size(Estimation.ThetaMask,1));
            end
        end
    end
end
% % % % % % % % % % % % % % % % % % % % % % % %
for type=1:4
    
    % Set up interaction between firm Xs and unobservables
    if Set.model == 'model 1'
        UoType  = SchoolsM.Mu*Params.betaK(type,1) + SchoolsM.Price*Params.betaK(type,2) + DistanceM*Params.betaK(type,3);
        Uv      = SchoolsM.Mu*Params.betai(:,1)'; % Random Coefficient Part of Utility
        Uv      = reshape(Uv,size(Uv,1),1,size(Uv,2));
    elseif Set.model == 'model 2'
        UoType  = SchoolsM.Mu*Params.betaK(type,1) + SchoolsM.Price*Params.betaK(type,2) + DistanceM*Params.betaK(type,3) + SchoolsM.Ze(:,type)*Params.betaK(type,4) + SchoolsM.Zy(:,type)*Params.betaK(type,5);
        Uv      = SchoolsM.Mu*Params.betai(:,1)' + SchoolsM.Ze(:,type)*Params.betai(:,2)' + SchoolsM.Zy(:,type)*Params.betai(:,3)'; % Random Coefficient Part of Utility
        Uv      = reshape(Uv,size(Uv,1),1,size(Uv,2));
    end
    % get maxutility to norm vector of utility to avoid overflow
    maxuij              = max(DeltaM)+max(max(UoType))+max(max(Uv));
    num                 = exp(DeltaM+UoType+Uv-maxuij);
    share_ijv{1,type}   = bsxfun(@rdivide,num,sum(num,1));
    share_ij            = sum(share_ijv{1,type}.*Params.w,3);
    Shares(:,type)      = sum(share_ij.*CweightsMTypes(:,type)',2);
    S                   = S + sum(share_ij.*CweightsMAll(:,type)',2);
    
    if nargout > 5
        
        if Set.model == 'model 1'
            moMy     = [share_ij'*[SchoolsM.Mu SchoolsM.Price] sum(share_ij.*DistanceM,1)']'*CweightsMTypes(:,type);
            
            %             moVy     = [ std(repmat(SchoolsM.Mu,size(share_ij,2),1),reshape(share_ij.*CweightsMTypes(:,type)',[],1));
            %                 std(repmat(SchoolsM.Price,size(share_ij,2),1),reshape(share_ij.*CweightsMTypes(:,type)',[],1));
            %                 std(reshape(DistanceM,[],1),reshape(share_ij.*CweightsMTypes(:,type)',[],1))];
            moM = [moM ; moMy];
            %             moV = [moV ; moVy];
            
        elseif Set.model == 'model 2'
            moMy     = [share_ij'*[SchoolsM.Mu SchoolsM.Price] sum(share_ij.*DistanceM,1)' share_ij'*[SchoolsM.Ze(:,type) SchoolsM.Zy(:,type)]]'*CweightsMTypes(:,type);
            
            %             moVy     = [ std(repmat(SchoolsM.Mu,size(share_ij,2),1),reshape(share_ij.*CweightsMTypes(:,type)',[],1));
            %                 std(repmat(SchoolsM.Price,size(share_ij,2),1),reshape(share_ij.*CweightsMTypes(:,type)',[],1));
            %                 std(reshape(DistanceM,[],1),reshape(share_ij.*CweightsMTypes(:,type)',[],1));
            %                 std(repmat(SchoolsM.Ze(:,type),size(share_ij,2),1),reshape(share_ij.*CweightsMTypes(:,type)',[],1));
            %                 std(repmat(SchoolsM.Zy(:,type),size(share_ij,2),1),reshape(share_ij.*CweightsMTypes(:,type)',[],1))];
            moM = [moM ; moMy];
            %             moV = [moV ; moVy];
        end
        
        
        
    end
    
    if nargout > 3
        
        if height(SchoolsM)<1400
            dSdDelta= dSdDelta+sum((mmx('mult',mmx('mult',permute(-share_ijv{1,type},[1 3 2]),diag(squeeze(Params.w(1,1,:)))),permute(share_ijv{1,type},[3 1 2]))).*reshape(CweightsMAll(:,type),1,1,[]),3);
        else
            for n=1:length(CweightsMTypes(:,type))
                
                dSdDelta = dSdDelta + -squeeze(share_ijv{1,type}(:,n,:))*diag(squeeze(Params.w(1,1,:)))*squeeze(share_ijv{1,type}(:,n,:))'*CweightsMAll(n,type);
                
            end
            
        end
        
        if nargout > 4
            
            if height(SchoolsM)<1400
                for i=Estimation.TypesTheta{type,1}
                    
                    if Estimation.ThetaMask{i,2}==1
                        
                        if Estimation.ThetaMask{i,1}==0
                            
                            dSdTheta(:,i) = dSdTheta(:,i) + squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv{1,type},'tn'))),CweightsMAll(:,type)))*squeeze(Params.w(1,1,:));
                            
                            if nargout > 6
                                
                                dMo_dTheta(1+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv{1,type},'tn'))),SchoolsM.Mu,'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(2+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv{1,type},'tn'))),SchoolsM.Price,'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(3+(type-1)*5,i) = mmx('mult',squeeze(sum((share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv{1,type},'tn'))).*DistanceM,1)),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(4+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv{1,type},'tn'))),SchoolsM.Ze(:,type),'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(5+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv{1,type},'tn'))),SchoolsM.Zy(:,type),'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                
                            end
                            
                        elseif Estimation.ThetaMask{i,1}==1
                            
                            dSdTheta(:,i) = dSdTheta(:,i) + squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv{1,type},'tn')),CweightsMAll(:,type)))*squeeze(Params.w(1,1,:));
                            
                            if nargout > 6
                                
                                dMo_dTheta(1+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv{1,type},'tn')),SchoolsM.Mu,'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(2+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv{1,type},'tn')),SchoolsM.Price,'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(3+(type-1)*5,i) = squeeze(mmx('mult',sum((bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv{1,type},'tn'))).*DistanceM,1),CweightsMTypes(:,type)))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(4+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv{1,type},'tn')),SchoolsM.Ze(:,type),'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(5+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Mu,mmx('mult',SchoolsM.Mu,share_ijv{1,type},'tn')),SchoolsM.Zy(:,type),'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                
                            end
                            
                        end
                        
                    elseif Estimation.ThetaMask{i,2}==2
                        
                        if Estimation.ThetaMask{i,1}==0
                            dSdTheta(:,i) = dSdTheta(:,i) + squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Price,mmx('mult',SchoolsM.Price,share_ijv{1,type},'tn'))),CweightsMAll(:,type)))*squeeze(Params.w(1,1,:));
                            
                            if nargout > 6
                                
                                dMo_dTheta(1+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Price,mmx('mult',SchoolsM.Price,share_ijv{1,type},'tn'))),SchoolsM.Mu,'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(2+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Price,mmx('mult',SchoolsM.Price,share_ijv{1,type},'tn'))),SchoolsM.Price,'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(3+(type-1)*5,i) = mmx('mult',squeeze(sum((share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Price,mmx('mult',SchoolsM.Price,share_ijv{1,type},'tn'))).*DistanceM,1)),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(4+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Price,mmx('mult',SchoolsM.Price,share_ijv{1,type},'tn'))),SchoolsM.Ze(:,type),'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(5+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Price,mmx('mult',SchoolsM.Price,share_ijv{1,type},'tn'))),SchoolsM.Zy(:,type),'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                
                            end
                            
                        end
                        
                    elseif Estimation.ThetaMask{i,2}==3
                        
                        if Estimation.ThetaMask{i,1}==0
                            
                            dSdTheta(:,i) = dSdTheta(:,i) + squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,DistanceM,sum(bsxfun(@times,DistanceM,share_ijv{1,type}),1))),CweightsMAll(:,type)))*squeeze(Params.w(1,1,:));
                            
                            if nargout > 6
                                
                                dMo_dTheta(1+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,DistanceM,sum(bsxfun(@times,DistanceM,share_ijv{1,type}),1))),SchoolsM.Mu,'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(2+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,DistanceM,sum(bsxfun(@times,DistanceM,share_ijv{1,type}),1))),SchoolsM.Price,'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(3+(type-1)*5,i) = squeeze(mmx('mult',squeeze(sum((share_ijv{1,type}.*bsxfun(@minus,DistanceM,sum(bsxfun(@times,DistanceM,share_ijv{1,type}),1))).*DistanceM,1)),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(4+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,DistanceM,sum(bsxfun(@times,DistanceM,share_ijv{1,type}),1))),SchoolsM.Ze(:,type),'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(5+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,DistanceM,sum(bsxfun(@times,DistanceM,share_ijv{1,type}),1))),SchoolsM.Zy(:,type),'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                
                            end
                            
                        end
                        
                    elseif Estimation.ThetaMask{i,2}==4
                        
                        if Estimation.ThetaMask{i,1}==0
                            
                            dSdTheta(:,i) = dSdTheta(:,i) + squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv{1,type},'tn'))),CweightsMAll(:,type)))*squeeze(Params.w(1,1,:));
                            
                            if nargout > 6
                                
                                dMo_dTheta(1+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv{1,type},'tn'))),SchoolsM.Mu,'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(2+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv{1,type},'tn'))),SchoolsM.Price,'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(3+(type-1)*5,i) = mmx('mult',squeeze(sum((share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv{1,type},'tn'))).*DistanceM,1)),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(4+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv{1,type},'tn'))),SchoolsM.Ze(:,type),'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(5+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv{1,type},'tn'))),SchoolsM.Zy(:,type),'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                
                            end
                            
                        elseif Estimation.ThetaMask{i,1}==1
                            
                            dSdTheta(:,i) = dSdTheta(:,i) + squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv{1,type},'tn')),CweightsMAll(:,type)))*squeeze(Params.w(1,1,:));
                            
                            if nargout > 6
                                
                                dMo_dTheta(1+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv{1,type},'tn')),SchoolsM.Mu,'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(2+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv{1,type},'tn')),SchoolsM.Price,'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(3+(type-1)*5,i) = squeeze(mmx('mult',sum((bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv{1,type},'tn'))).*DistanceM,1),CweightsMTypes(:,type)))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(4+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv{1,type},'tn')),SchoolsM.Ze(:,type),'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(5+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Ze(:,type),mmx('mult',SchoolsM.Ze(:,type),share_ijv{1,type},'tn')),SchoolsM.Zy(:,type),'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                
                            end
                        end
                        
                    elseif Estimation.ThetaMask{i,2}==5
                        
                        if Estimation.ThetaMask{i,1}==0
                            
                            dSdTheta(:,i) = dSdTheta(:,i) + squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv{1,type},'tn'))),CweightsMAll(:,type)))*squeeze(Params.w(1,1,:));
                            
                            if nargout > 6
                                
                                dMo_dTheta(1+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv{1,type},'tn'))),SchoolsM.Mu,'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(2+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv{1,type},'tn'))),SchoolsM.Price,'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(3+(type-1)*5,i) = mmx('mult',squeeze(sum((share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv{1,type},'tn'))).*DistanceM,1)),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(4+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv{1,type},'tn'))),SchoolsM.Ze(:,type),'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(5+(type-1)*5,i) = mmx('mult',squeeze(mmx('mult',(share_ijv{1,type}.*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv{1,type},'tn'))),SchoolsM.Zy(:,type),'tn')),CweightsMTypes(:,type),'tn')'*squeeze(Params.w(1,1,:));
                                
                            end
                            
                        elseif Estimation.ThetaMask{i,1}==1
                            
                            dSdTheta(:,i) = dSdTheta(:,i) + squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv{1,type},'tn')),CweightsMAll(:,type)))*squeeze(Params.w(1,1,:));
                            
                            if nargout > 6
                                
                                dMo_dTheta(1+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv{1,type},'tn')),SchoolsM.Mu,'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(2+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv{1,type},'tn')),SchoolsM.Price,'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(3+(type-1)*5,i) = squeeze(mmx('mult',sum((bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv{1,type},'tn'))).*DistanceM,1),CweightsMTypes(:,type)))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(4+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv{1,type},'tn')),SchoolsM.Ze(:,type),'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                dMo_dTheta(5+(type-1)*5,i) = squeeze(mmx('mult',squeeze(mmx('mult',bsxfun(@times,share_ijv{1,type},reshape(Params.dbetai(:,Estimation.ThetaMask{i,3}),1,1,[])).*bsxfun(@minus,SchoolsM.Zy(:,type),mmx('mult',SchoolsM.Zy(:,type),share_ijv{1,type},'tn')),SchoolsM.Zy(:,type),'tn')),CweightsMTypes(:,type),'tn'))'*squeeze(Params.w(1,1,:));
                                
                            end
                        end
                    end
                    
                    
                end
            else
                
            end
        end
    end
    
end

if nargout > 5
    
    MM = moM;
    % MM = [moM,moV];
    if nargout > 6
        
        dDeltadThetanorm =  zeros(size(S,1),size(Estimation.ThetaMask,1))  ;
        
        dSdDelta= dSdDelta + diag(S);
        %   dDeltadThetanorm = -dSdDelta\dSdTheta;
        %   TempDdDt = -inv(dSdDelta)*dSdTheta;
        %   dDeltadThetanorm = TempDdDt  ;
        %   TempDdDt = (diag(1./S(SchoolsM.NormId==0,:))*dSdDelta(SchoolsM.NormId==0,SchoolsM.NormId==0))\(diag(1./S(SchoolsM.NormId==0,:))*dSdTheta(SchoolsM.NormId==0,:));
        TempDdDt = -inv(dSdDelta(SchoolsM.NormId==0,SchoolsM.NormId==0))*dSdTheta(SchoolsM.NormId==0,:);
        
        dDeltadThetanorm(SchoolsM.NormId==0,:) = TempDdDt  ;
        
        for type=1:4
            
            share_ij = sum(share_ijv{1,type}.*Params.w,3);
            
            Diag = diag(share_ij(:,1));
            
            for n=2:size(share_ij,2)
                
                Diag = cat(3,Diag,diag(share_ij(:,n)));
                
            end
            
            dMo_dSdDdDdT(1+(type-1)*5,:)= sum((mmx('mult',(mmx('mult',(mmx('mult',mmx('mult',permute(-share_ijv{1,type},[1 3 2]),diag(squeeze(Params.w(1,1,:)))),permute(share_ijv{1,type},[3 1 2])))+Diag,dDeltadThetanorm)),SchoolsM.Mu,'tn')).*reshape(CweightsMTypes(:,type),1,1,[]),3)';
            dMo_dSdDdDdT(2+(type-1)*5,:)= sum((mmx('mult',(mmx('mult',(mmx('mult',mmx('mult',permute(-share_ijv{1,type},[1 3 2]),diag(squeeze(Params.w(1,1,:)))),permute(share_ijv{1,type},[3 1 2])))+Diag,dDeltadThetanorm)),SchoolsM.Price,'tn')).*reshape(CweightsMTypes(:,type),1,1,[]),3)';
            dMo_dSdDdDdT(3+(type-1)*5,:)= sum((sum((mmx('mult',(mmx('mult',mmx('mult',permute(-share_ijv{1,type},[1 3 2]),diag(squeeze(Params.w(1,1,:)))),permute(share_ijv{1,type},[3 1 2])))+Diag,dDeltadThetanorm)).*reshape(DistanceM,size(DistanceM,1),1,size(DistanceM,2)),1)).*reshape(CweightsMTypes(:,type),1,1,[]),3);
            dMo_dSdDdDdT(4+(type-1)*5,:)= sum((mmx('mult',(mmx('mult',(mmx('mult',mmx('mult',permute(-share_ijv{1,type},[1 3 2]),diag(squeeze(Params.w(1,1,:)))),permute(share_ijv{1,type},[3 1 2])))+Diag,dDeltadThetanorm)),SchoolsM.Ze(:,type),'tn')).*reshape(CweightsMTypes(:,type),1,1,[]),3)';
            dMo_dSdDdDdT(5+(type-1)*5,:)= sum((mmx('mult',(mmx('mult',(mmx('mult',mmx('mult',permute(-share_ijv{1,type},[1 3 2]),diag(squeeze(Params.w(1,1,:)))),permute(share_ijv{1,type},[3 1 2])))+Diag,dDeltadThetanorm)),SchoolsM.Zy(:,type),'tn')).*reshape(CweightsMTypes(:,type),1,1,[]),3)';
            
        end
        
        DMM_dTheta = (dMo_dSdDdDdT + dMo_dTheta);
        
    end
end


end




% % % % % % % % % % % % % % % % % % % % % % % %





