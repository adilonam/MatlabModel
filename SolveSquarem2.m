function [delta2, flag, norm_maxShares, norm_max,  iter]=SolveSquarem2(delta0,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params)


iter=1;

delta0=delta0-delta0(SchoolsM.NormId==1);
norm_maxSharesOld=1;
norm_max=1; norm_maxShares=1;showit=1;normaxold=10; count=0;countM=0;countS=0;flag=0;

while norm_max > Set.tol_inner && norm_maxShares>Set.tol_inner  && iter < 1000

    shares1 =RC_shares(delta0,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
    delta1 =  delta0 + log(SchoolsM.Share) - log(shares1);    
     
    delta2=delta1-delta1(SchoolsM.NormId==1);% Apply Normalization
    
    shares2 =RC_shares(delta2,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
    delta3 =  delta2 + log(SchoolsM.Share) - log(shares2);    
     
    delta4=delta3-delta3(SchoolsM.NormId==1);% Apply Normalization
          
    v=delta4+delta0-(delta2+delta2);  
    r=delta2-delta0;    
    alph=(v'*r)./(v'*v);
    
    delta2=(delta0)-2*alph*r+alph^2*v;
    
    [norm_max, ~]= max(abs(delta2 -delta0 ));
    [norm_maxShares, ~]= max(abs(SchoolsM.Share-shares1));
    
    if iter/100>=showit && iter/1001<showit
        fprintf('Dif of shares of %17.12f and diference of expdelta of %17.12f  at iteration %10.0f   \n ', [norm_maxShares norm_max  iter] )
        showit=showit+1;
    end
    
    % START ROBUSTNESS SECTION
    % If iteration gets stuck, get out after some time.
    if norm_max==normaxold || norm_maxShares==norm_maxSharesOld
        count=count+1;
        %  fprintf('Not updating it seems %4.0f \n',count)
        if count>=11
            flag=1;
            disp('jumping ship')
            break
        end
    end
    % If iteration goes wrong way, get out after some time.
    if norm_max>normaxold
        countS=countS+1;
        if countM>=101
            flag=2;
             disp('jumping ship')
            break
        end
    end
    
    if max(isnan(delta2))>=1
        delta2 = delta0 ;
        flag=3;
        break
    end

    % ----------------------------------------------------------------------------
    % END ROBUSTNESS SECTION
    
    delta0 = (delta2);
    normaxold=norm_max;norm_maxSharesOld=norm_maxShares;
    iter = iter + 1;
end
end






