function [Params] = GetParams(Estimation,Set,Theta2)

if Set.model == 'model 1'
    
    % Model 1 
    
    BetaQEdu    = Theta2(Estimation.BetaEdu);
    BetaQPoor   = Theta2(Estimation.BetaPoor);     
    AlphaEdu    = Theta2(Estimation.AlphaEdu);      
    AlphaPoor   = Theta2(Estimation.AlphaPoor);     
    LambdaEdu   = Theta2(Estimation.LambdaEdu);    
    LambdaPoor  = Theta2(Estimation.LambdaPoor);
    
    Params.betaK    =   [BetaQPoor          AlphaPoor               LambdaEdu(1)+LambdaPoor ;   
                        0                   0                       LambdaEdu(1)            ;
                        BetaQEdu+BetaQPoor  AlphaEdu+AlphaPoor      LambdaEdu(2)+LambdaPoor ;
                        BetaQEdu            AlphaEdu                LambdaEdu(2)            ];
  
    BetaRC   = Theta2(Estimation.BetaRC);
    
    Params.betai        = [BetaRC*Estimation.drawsN(:,1) zeros(size(Estimation.drawsN(:,1))) zeros(size(Estimation.drawsN(:,1)))];
    Params.dbetai       = Estimation.drawsN(:,1);
    Params.w(1,1,:)     = Estimation.drawsW;
    Params.dbetamask    = [1 1];
    Params.lb           = 0;
    Params.ub           = Inf;
    
elseif Set.model == 'model 2'
    
    % Model 2 

    BetaQEdu     = Theta2(Estimation.BetaQEdu);      
    BetaQPoor    = Theta2(Estimation.BetaQPoor);    
    AlphaEdu     = Theta2(Estimation.AlphaEdu);      
    AlphaPoor    = Theta2(Estimation.AlphaPoor);     
    BetaZeEdu    = Theta2(Estimation.BetaZeEdu);     
    BetaZePoor   = Theta2(Estimation.BetaZePoor);    
    BetaZpEdu    = Theta2(Estimation.BetaZpEdu);    
    BetaZpPoor   = Theta2(Estimation.BetaZpPoor);    
    LambdaEdu    = Theta2(Estimation.LambdaEdu);     
    LambdaPoor   = Theta2(Estimation.LambdaPoor);   
       
    Params.betaK    =   [BetaQPoor          AlphaPoor               LambdaEdu(1)+LambdaPoor     BetaZeEdu(1)+BetaZePoor     BetaZpEdu(1)+BetaZpPoor  ;   
                        0                   0                       LambdaEdu(1)                BetaZeEdu(1)                BetaZpEdu(1)             ;
                        BetaQEdu+BetaQPoor  AlphaEdu+AlphaPoor      LambdaEdu(2)+LambdaPoor     BetaZeEdu(2)+BetaZePoor     BetaZpEdu(2)+BetaZpPoor  ;
                        BetaQEdu            AlphaEdu                LambdaEdu(2)                BetaZeEdu(2)                BetaZpEdu(2)             ];
                 
    BetaQRC     = Theta2(Estimation.BetaQRC);       
    BetaZeRC    = Theta2(Estimation.BetaZeRC);      
    BetaZpRC    = Theta2(Estimation.BetaZpRC);      
    BetaZcorRC  = Theta2(Estimation.BetaZcorRC); 
    
    Params.betai        = [BetaQRC*Estimation.drawsN(:,1), BetaZeRC*Estimation.drawsN(:,2), BetaZcorRC*Estimation.drawsN(:,2) + BetaZpRC*Estimation.drawsN(:,3)];
    Params.dbetai       = Estimation.drawsN(:,1:3);
    Params.w(1,1,:)     = Estimation.drawsW;
    Params.dbetamask    = [1 1; 4 2; 5 2; 4 3];
    Params.lb           = [0      0       -Inf    0];
    Params.ub           = [Inf    Inf     Inf     Inf];
end

end 
