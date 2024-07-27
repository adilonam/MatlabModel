function Dev = DeviationFoc(x,DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Theta2)

P = x(1:height(SchoolsM));
Q = x(height(SchoolsM)+1:end);

SchoolsM.Price  = P;
SchoolsM.Mu     = Q;

[S,~,~,~,dSdq,dSdp] = RC_shares(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Theta2);

MarkDownQ   =   S./dS_dq; 
MarkUpP     =   S./dS_dp; 

Pnew = Schools.C0+C1*Q+MarkUpP;
Qnew = ;

Dev = [Pnew-P;lQnew-Q];
end