function [MarkdownQ,MarkupP,dSdq,dSdp] = SupplyMarks(S,share_ijv,SchoolsM,CweightsMAll,CweightsMTypes,Estimation,Set,Theta2,Model)

Params = GetParams(Estimation,Set,Theta2);

dSdq         = zeros(height(SchoolsM));
dSdp         = zeros(height(SchoolsM));

if strcmp(Model,'Model1')
    
    for type=1:4
        
        share_ijv_a=bsxfun(@times,share_ijv{1,type},reshape(Params.betaK(type,1)+Params.betai(:,1)',1,1,[]));
        share_ijv_b=Params.betaK(type,2)*share_ijv{1,type};
        
        if height(SchoolsM)<1400
            
            dSdq= dSdq+sum((mmx('mult',mmx('mult',permute(-share_ijv_a,[1 3 2]),diag(squeeze(Params.w(1,1,:)))),permute(share_ijv{1,type},[3 1 2]))).*reshape(CweightsMAll(:,type),1,1,[]),3) + diag(sum(sum(share_ijv_a.*Params.w,3).*CweightsMAll(:,type)',2));
            dSdp= dSdp+sum((mmx('mult',mmx('mult',permute(-share_ijv_b,[1 3 2]),diag(squeeze(Params.w(1,1,:)))),permute(share_ijv{1,type},[3 1 2]))).*reshape(CweightsMAll(:,type),1,1,[]),3) + diag(sum(sum(share_ijv_b.*Params.w,3).*CweightsMAll(:,type)',2));
        
        else
            
            for n=1:length(CweightsMTypes(:,type))
                
                dSdq = dSdq + -squeeze(share_ijv_a(:,n,:))*diag(squeeze(Params.w(1,1,:)))*squeeze(share_ijv{1,type}(:,n,:))'*CweightsMAll(n,type) + diag(sum(share_ijv_a(:,n,:).*Params.w,3).*CweightsMAll(n,type)');
                dSdp = dSdp + -squeeze(share_ijv_b(:,n,:))*diag(squeeze(Params.w(1,1,:)))*squeeze(share_ijv{1,type}(:,n,:))'*CweightsMAll(n,type) + diag(sum(share_ijv_b(:,n,:).*Params.w,3).*CweightsMAll(n,type)');
                
            end
            
        end
        
    end
    
elseif strcmp(Model,'Model2')
    
end

MarkdownQ   = S./dSdq;
MarkupP     = S./dSdp;

end