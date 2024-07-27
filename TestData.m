

%     - schools (t): 
%
%         * sIndex [MarketId,Year,SchoolId,SchoolDistr,ChainId,Lat,Lon]
%         * share [Share,ShareType,NumType,]
%         * charact [Price,Mu,Ze,Zy,Delta,Private,Charter,ForProfit,Religious,Emblematic]
%         * Instrumets [Instruments]
%         * NormId [NormId]
% 
%     - consumers (t):
% 
%         * cIndex [marketId,Node,Lat,Lon]
% 
%     - draws (t):
% 
%         * nodes
%         * weights
% 
%     - distance (a,2):
% 
%         * rows: schools, columns: consumer nodes
% 
%     - Cweights (t):
% 
%         * rows: consumers, columns: years, third dim: types
% 
%     - Moments (s), fields:
% 
%         * Moments.MM: Micromoments (MarketId,Year,Type,Charact,Moment,Variance,Nobs,Sanitycheck)
%         * Moments.WMM: MM weighting matrix
%         * Moments.RDM: MarketId,Year,Type,Node,ScoreAc,ScoreLd,Treat

F = importdata([stemWorked filesep '0. Final Model/1. FakeData' filesep 'Firms.csv']);
D = importdata([stemWorked filesep '0. Final Model/1. FakeData' filesep 'Delta.csv']);
C = importdata([stemWorked filesep '0. Final Model/1. FakeData' filesep 'Consumers.csv']);
M = importdata([stemWorked filesep '0. Final Model/1. FakeData' filesep 'MicroMoments.csv']);

% 1. Schools

MarketId    = F(:,1);
Year        = F(:,2)+7;
SchoolId    = F(:,3);
SchoolDistr = F(:,1)+randi([1 5],size(MarketId,1),1);
ChainId     = (ones(size(F(:,9)))-F(:,9)).*(rand(size(F(:,9))) >= 0.7).*randi([1 10],size(F(:,9)));
Lat         = F(:,8);
Lon         = F(:,7);
Share       = F(:,4);
Share(isnan(Share)) = rand(sum(isnan(Share)),1)/100;
ShareType   = repmat(Share(:,1),1,4)+rand(size(repmat(Share(:,1),1,4)))/100;

Price       = rand(size(F(:,6)));
Mu          = randn(size(F(:,5)));

Delta       = D/std(D);
Private     = F(:,9);
Private(1)  = 0;

Price(Private==0) = 0;

Charter     = ones(size(F(:,8),1),1).*(rand(size(F(:,8))) >= 0.9).*(Private==1);


ForProfit   = F(:,11);
ForProfit(Private==0) = 0;

Religious   = F(:,12).*(Private==1);

Emblematic  = ones(size(F(:,8),1),1).*(rand(size(F(:,8))) >= 0.9).*(Private==0);
Instruments = rand(size(F(:,5),1),17);
NormId      = zeros(size(F(:,1)));

markets = unique(MarketId)';
years = unique(Year)';

NumType = zeros(size(ShareType));
ZeB = zeros(size(ShareType));
ZyB = zeros(size(ShareType));

ZeLOO = zeros(size(ShareType));
ZyLOO = zeros(size(ShareType));

for m=markets
    for y=years
        
        idx=find(MarketId == m & Year == y);
        NormId(idx(1)) = 1;
        Share(idx,1) = Share(idx,1)./sum(Share(idx,1),1);
        
        temp = rand(1,4)*10^min(max(floor(rand()*10),3),5);
        
        for type = 1:4
            ShareType(idx,type) = round(ShareType(idx,type)./sum(ShareType(idx,type),1),7);
            NumType(idx,type)   = temp(1,type)*ShareType(idx,type);
        end
        
        % Basic version
        
        ZeB(idx,:)  = repmat((NumType(idx,3) + NumType(idx,4))./(sum(NumType(idx,:),2)),1,4);
        
        ZyB(idx,:)  = repmat((NumType(idx,2) + NumType(idx,4))./(sum(NumType(idx,:),2)),1,4);
 
        % Leave one out version
        
        ZeLOO(idx,1)  = (NumType(idx,3) + NumType(idx,4))./(sum(NumType(idx,:),2)-1);
        ZeLOO(idx,2)  = (NumType(idx,3) + NumType(idx,4))./(sum(NumType(idx,:),2)-1);
        ZeLOO(idx,3)  = (NumType(idx,3) + NumType(idx,4)-1)./(sum(NumType(idx,:),2)-1);
        ZeLOO(idx,4)  = (NumType(idx,3) + NumType(idx,4)-1)./(sum(NumType(idx,:),2)-1);
        
        ZyLOO(idx,1)  = (NumType(idx,2) + NumType(idx,4))./(sum(NumType(idx,:),2)-1);
        ZyLOO(idx,2)  = (NumType(idx,2) + NumType(idx,4)-1)./(sum(NumType(idx,:),2)-1);
        ZyLOO(idx,3)  = (NumType(idx,2) + NumType(idx,4))./(sum(NumType(idx,:),2)-1);
        ZyLOO(idx,4)  = (NumType(idx,2) + NumType(idx,4)-1)./(sum(NumType(idx,:),2)-1);
        
    end
end

Schools = table(MarketId,Year,SchoolId,SchoolDistr,ChainId,Lat,Lon,Share,ShareType,NumType,Price,Mu,ZeB,ZyB,ZeLOO,ZyLOO,Delta,Private,Charter,ForProfit,Religious,Emblematic,Instruments,NormId);

clear MarketId Year SchoolId SchoolDistr ChainId Lat Lon Share ShareType NumType Price Mu Ze Zy Delta Private Charter ForProfit Religious Emblematic Instruments NormId

% 2. Consumers

MarketId    = C(C(:,2)==2011 & C(:,3)==1,1);
Node        = C(C(:,2)==2011 & C(:,3)==1,4);
Lat         = C(C(:,2)==2011 & C(:,3)==1,6);
Lon         = C(C(:,2)==2011 & C(:,3)==1,5);

Consumers = table(MarketId,Node,Lat,Lon);

% 3. Distance

Distance = zeros(size(Schools,1),size(Consumers,1));

for j=1:size(Schools,1)
    cIdx=(Consumers.MarketId==Schools.MarketId(j));
    Node(cIdx)=(1:sum(cIdx,1))';
    Distance(j,cIdx) = deg2km(distance(repmat(Schools.Lat(j),sum(cIdx),1),repmat(Schools.Lon(j),sum(cIdx),1),Consumers.Lat(cIdx),Consumers.Lon(cIdx)))';
end

Consumers.Node = Node;

clear MarketId Node Lat Lon


% 4.  Cweights (t):


% Cweights(:,1,1)= C(C(:,2)==2007 & C(:,3)==1,7);
% Cweights(:,1,2)= C(C(:,2)==2007 & C(:,3)==2,7) + C(C(:,2)==2007 & C(:,3)==3,7);
% Cweights(:,1,3)= C(C(:,2)==2007 & C(:,3)==4,7) + C(C(:,2)==2007 & C(:,3)==5,7);
% Cweights(:,1,4)= C(C(:,2)==2007 & C(:,3)==6,7);
% 
% Cweights(:,2,1)= C(C(:,2)==2011 & C(:,3)==1,7);
% Cweights(:,2,2)= C(C(:,2)==2011 & C(:,3)==2,7) + C(C(:,2)==2011 & C(:,3)==3,7);
% Cweights(:,2,3)= C(C(:,2)==2011 & C(:,3)==4,7) + C(C(:,2)==2011 & C(:,3)==5,7);
% Cweights(:,2,4)= C(C(:,2)==2011 & C(:,3)==6,7);

CC = rand(height(Consumers),length(years),4);
Cweights.All = zeros(size(CC));
Cweights.Types = zeros(size(CC));


for m=markets
    for y=1:size(years,2)
       idx=find(Consumers.MarketId == m);
       Cweights.Types(idx,y,:) = CC(idx,y,:)./sum(CC(idx,y,:),1);
       Cweights.All(idx,y,:) = CC(idx,y,:)./sum(sum(CC(idx,y,:),1));
   end
end

% 5.  Moments (s)

% mm = reshape(repmat(markets,40,1),1,[])';
% yy = reshape(repmat([2014 2018],20,1),1,[])';
% yy = repmat(yy,size(markets,2),1);
% tt = reshape(repmat([1:4],5,1),1,[])';
% tt = repmat(tt,size(markets,2)*2,1);
% mo = repmat([1:5]',size(markets,2)*2*4,1);

MarketId = reshape(repmat(markets,40,1),1,[])';
Year     = reshape(repmat([2014 2018],20,1),1,[])';
Year     = repmat(Year,size(markets,2),1);
Type     = reshape(repmat([1:4],5,1),1,[])';
Type     = repmat(Type,size(markets,2)*2,1);
Charact  = repmat([1:5]',size(markets,2)*2*4,1);
SanityCheck = zeros(size(MarketId,1),1);

Moment      = rand(size(MarketId));
MVar   = rand(size(MarketId));
MobsN       = rand(size(MarketId));

for m=markets
   for y=years
       for t=1:4 
          Moment(MarketId==m & Year==y & Type == t,1)       = randn(5,1);©
          MVar(MarketId==m & Year==y & Type == t,1)         = randi(100,5,1)./100;
          MobsN(MarketId==m & Year==y & Type == t,1)        = rand(1,5)*10^min(max(floor(rand()*10),3),5);
          if y == 2014
              SanityCheck(MarketId==m & Year==y & Type == t,1)  = [1;1;0;1;1];
          else
              SanityCheck(MarketId==m & Year==y & Type == t,1)  = [1;1;1;1;1];
          end
       end  
   end
end

SanityCheck = logical(SanityCheck);
Moments.MM = table(MarketId,Year,Type,Charact,Moment,MVar,MobsN,SanityCheck);
Moments.WMM = diag(ones(size(MarketId,1),1))./size(MarketId,1);

clear MarketId Year Type Charact Moment MVar MobsN

MarketId    = markets(randi(size(markets,2),1000,1))'; 
Year        = years(randi(size(years,2),1000,1))'; 
Type        = randi(4,1000,1);

for i=1:size(MarketId,1)
    
    nodes   = unique(Consumers.Node(Consumers.MarketId==MarketId(i)));
    Node(i,1) = nodes(randi(size(nodes,1),1,1));
    schools = unique(Schools.SchoolId(Schools.MarketId==MarketId(i)));
    School(i,1) = schools(randi(size(schools,1),1,1));
end

ScoreAc     = randn(1000,1);
ScoreLd     = randn(1000,1).*(rand(size(ScoreAc)) >= 0.2);
MiScoreLd   = (ScoreLd==0);
ThresCross  = (ScoreAc>0 & ScoreLd>0);
Treat       = (ScoreAc>0 & ScoreLd>0).*(rand(size(ThresCross)) >= 0.3);

Moments.RDdata = sortrows(table(MarketId,Year,Type,Node,ScoreAc,ScoreLd,MiScoreLd,ThresCross,Treat,School),{'MarketId','Year'});

save([stemWorked filesep '0. Final Model/1. FakeData' filesep 'SchoolData_' datadate '.mat'],'Schools')   
save([stemWorked filesep '0. Final Model/1. FakeData' filesep 'ConsumerData_' datadate '.mat'],'Consumers','Cweights','Distance')
save([stemWorked filesep '0. Final Model/1. FakeData' filesep 'Estimation_' datadate '.mat'],'Moments')

