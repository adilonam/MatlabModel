function [Schools,Consumers,Cweights,Distance,Moments,Estimation,markets]=setupData(Set,stemWorked)





%{

Project:        JMP
Author:         Claudia Allende
Date created:   May 23rd, 2019
Date modified:  May 24th, 2019
Description:    This code sets up the data for estimation 

Input:

   - markets:

        - [xx xx]:  Any selection (markets go from 1 to 102)
        - [991]:    Small markets (under 100)  
        - [992]:    Medium markets (under 300)  
        - [993]:    All markets except Lima (101) 
        - []:       All markets including Lima (102)

    - years:
        - 2014-2018

    - model:

        - 'model 1': {p,q,d} 
        - 'model 2': {p,q,d,z} 

    - moments:

        - 'moments 1': Micro,IV
        - 'moments 2': Micro,IV,RD
        - 'moments 3': Micro,IV,RD,MLE

    - Types:

        Type 1: Poor + Low Educ
        Type 2: Low Educ
        Type 3: Poor + High Educ
        Type 4: High Educ
    
Output: 'Data': Structure with the following fields (t: table; a,n: array,dim; s: struct):    
    
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
%     - Cweights (s):
% 
%         * [rows: consumers, columns: years, third dim: types]
%         * Cweights.All
%         * Cweights.Types
%
%     - Moments (s), fields:
% 
%         * Moments.MM: Micromoments (MarketId,Year,Type,Charact,Moment,Variance,Nobs)
%         * Moments.WMM: MM weighting matrix
%         * Moments.RDM: RD Moments(MarketId,Year,Type,Node,ScoreAc,ScoreLd,Treat)

%}

% ,'shares','sIndex','normID','normXX','MM','MIndex','VMM','NMM','MMnames','Mmicro'


% Load data

load(['Data' filesep 'SchoolData_' Set.datadate '.mat'])   
load(['Data' filesep 'ConsumerData_' Set.datadate '.mat'])
load(['Data' filesep 'Estimation_' Set.datadate '.mat'])

% Choose Markets
if strcmp(Set.markets,'small')==1
    D=tabulate(Schools.MarketId(:,1));
    pick=(D(:,2)/size(Set.years,2)<100 & D(:,2)/size(Set.years,2)>20);
    markets=D(pick,1)';       % small markets
elseif strcmp(Set.markets,'medium')==1
    D=tabulate(Schools.MarketId(:,1));
    pick=(D(:,2)/size(Set.years,2)<300 & D(:,2)/size(Set.years,2)>=20);
    markets=D(pick,1)';          % all markets
elseif strcmp(Set.markets,'large')==1
    D=tabulate(Schools.MarketId(:,1));
    %pick=(D(:,1)~=56 & D(:,2)/size(Set.years,2)>20); % Exclude Lima
    pick=(D(:,1)~=312 & D(:,2)/size(Set.years,2)>20); % Exclude Santiago
    markets=D(pick,1)';
elseif strcmp(Set.markets,'all')==1
    D=tabulate(Schools.MarketId(:,1));
    pick=(D(:,2)/size(Set.years,2)>20);
    markets=D(pick,1)';          % all markets
end

%% Get feasible markets 
% 1) have firm data, have shares and have moments

% Get Row Index for moments
temp_mIndex=[];
for m=markets
    for y=Set.years
        temp_mIndex=[temp_mIndex; find(Moments.MM.MarketId==m & Moments.MM.Year==y & isnan(Moments.MM.Moment)==0)];
    end
end

markets=unique(Moments.MM.MarketId(temp_mIndex)); % update markets to include only markets with moments

% Get Row Index for schools
temp_sIndex=[];
for m=markets'
    for y=Set.years
        temp_sIndex=[temp_sIndex; find(Schools.MarketId==m & Schools.Year==y & (isnan(Schools.Share)==0))];
    end
end

markets=unique(Schools.MarketId(temp_sIndex)); % update markets to include schools restrictions

% Get Row Index for Consumers

temp_cIndex=[];
for m=markets'
    temp_cIndex=[temp_cIndex; find(Consumers.MarketId==m)];
end

% Get Row Index for RD sample

temp_RDIndex=[];
for m=markets'
    temp_RDIndex=[temp_RDIndex; find(Moments.RDdata.MarketId==m)];
end

% Adjust Peers model

if strcmp(Set.peersmodel,'basic')
    Schools = removevars(Schools,{'ZeLOO','ZyLOO'});
    Schools.Properties.VariableNames{'ZeB'} = 'Ze';
    Schools.Properties.VariableNames{'ZyB'} = 'Zy';
elseif strcmp(Set.peersmodel,'leave one out')
    Schools = removevars(Schools,{'ZeB','ZyB'});
    Schools.Properties.VariableNames{'ZeLOO'} = 'Ze';
    Schools.Properties.VariableNames{'ZyLOO'} = 'Zy';
end

% Put the exogenous regressors together

XX_exogenous = [Schools.Private Schools.Charter Schools.ForProfit Schools.Religious Schools.Emblematic];
Schools = addvars(Schools,XX_exogenous,'Before','Private');
Schools = removevars(Schools,{'Private','Charter','ForProfit','Religious','Emblematic'});

%% Now get final group 

Schools         = Schools(temp_sIndex,:);
Consumers       = Consumers(temp_cIndex,:);
Cweights.Types  = Cweights.Types(temp_cIndex,:,:);
Cweights.All    = Cweights.All(temp_cIndex,:,:);
Distance        = Distance(temp_sIndex,temp_cIndex);
Moments.MM      = Moments.MM(temp_mIndex,:);
Moments.WMM     = Moments.WMM(temp_mIndex,temp_mIndex);
Moments.RDdata  = Moments.RDdata(temp_RDIndex,:);

% Schools = renamevars(Schools, 'ZeB', 'Ze');
% Schools = renamevars(Schools, 'ZyB', 'Zy');

%% Standardize Mu, Price, Distance, Ze, and Zy by Market

for m = markets'
    for y = Set.years
        Schools.Mu(Schools.MarketId==m & Schools.Year==y,:) = (Schools.Mu(Schools.MarketId==m & Schools.Year==y,:) - (Schools.Mu(Schools.MarketId==m & Schools.Year==y & Schools.NormId==1,:)))./std(Schools.Mu(Schools.MarketId==m & Schools.Year==y,:));
        Schools.Price(Schools.MarketId==m & Schools.Year==y,:) = (Schools.Price(Schools.MarketId==m & Schools.Year==y,:) - (Schools.Price(Schools.MarketId==m & Schools.Year==y & Schools.NormId==1,:)))./std(Schools.Price(Schools.MarketId==m & Schools.Year==y,:));
        for k=1:4
            Schools.Ze(Schools.MarketId==m & Schools.Year==y,k) = (Schools.Ze(Schools.MarketId==m & Schools.Year==y,k) - (Schools.Ze(Schools.MarketId==m & Schools.Year==y & Schools.NormId==1,k)))./std(Schools.Ze(Schools.MarketId==m & Schools.Year==y,k));
            Schools.Zy(Schools.MarketId==m & Schools.Year==y,k) = (Schools.Zy(Schools.MarketId==m & Schools.Year==y,k) - (Schools.Zy(Schools.MarketId==m & Schools.Year==y & Schools.NormId==1,k)))./std(Schools.Zy(Schools.MarketId==m & Schools.Year==y,k));
        end
        Distance(Schools.MarketId==m & Schools.Year==y,Consumers.MarketId==m)=(Distance(Schools.MarketId==m & Schools.Year==y,Consumers.MarketId==m)-min(abs(Distance(Schools.MarketId==m & Schools.Year==y & Schools.NormId==1,Consumers.MarketId==m))))./std(Distance(Schools.MarketId==m & Schools.Year==y,Consumers.MarketId==m));
    end
end

% Prep IV Regression
 Schools.ChainId(Schools.ChainId==0) = NaN;
 Estimation.ChainFE = dummyvar(Schools.ChainId);
 Estimation.ChainFE(isnan(Estimation.ChainFE)) = 0;
 
 Rfe=[];Cfe=[];
 MarketID=unique(Schools.MarketId);
 TimeID=unique(Schools.Year);
 
 counter=1;
 
 for j=1:length(MarketID)
     for t=1:length(TimeID)
         idm=MarketID(j);
         idt=TimeID(t);
         
         rj=find(idm==Schools.MarketId & idt==Schools.Year);
         Rfe=[Rfe;rj];
         Cfe=[Cfe; counter*ones(length(rj),1)];
         counter=counter+1;
     end
 end
 
Estimation.MarketYearFE=sparse(Rfe,Cfe,ones(length(Rfe),1));

% Fix weights for micro moments

Moments.MM.SanityCheck((abs(Moments.MM.Moment)==0 | Moments.MM.MobsN<30  | (Moments.MM.Charact==3 & Moments.MM.Year==2014)==1))=0;

temp=Moments.MM.MVar(Moments.MM.SanityCheck==1).^(-1);
if sum(temp==Inf)>0
error('Weights for micromoments include Inf')    
end
if sum(temp==0)>0
error('Weights for micromoments include zeros')    
end

for ii=1:5
tempweight=(Moments.MM.Charact(Moments.MM.SanityCheck==1)==ii);
temp(tempweight)=temp(tempweight)/sum(temp(tempweight)); 
end

weightsM=temp;
Moments.WMM=sparse(1:length(Moments.MM.MVar(Moments.MM.SanityCheck==1)),1:length(Moments.MM.MVar(Moments.MM.SanityCheck==1)),weightsM);


%% Prepare RD Data

DataY   = NaN(height(Moments.RDdata),6);

[~,locb1] = ismember(Moments.RDdata.School,Schools.SchoolId);
DataY(:,1) = Schools.Mu(locb1);
DataY(:,2) = Schools.Price(locb1);

DataY(:,4) = Schools.Price(locb1);
DataY(:,5) = Schools.Ze(locb1);
DataY(:,6) = Schools.Zy(locb1);

for m=unique(Moments.RDdata.MarketId)'
    for y=unique(Moments.RDdata.Year)'
        
        % By Market, by year
        
        SchoolsM = Schools(Schools.MarketId==m & Schools.Year==y,:);
        DistanceM = Distance(Schools.MarketId==m & Schools.Year==y,Consumers.MarketId==m);
        
        % % % % % % % % % % % % % % % % % % % % % % % % Only for Model 2
        
        rdIndex = (Moments.RDdata.MarketId == m & Moments.RDdata.Year == y);
        
        nodes = Moments.RDdata.Node(rdIndex);
        
        [~,locb1] = ismember(Moments.RDdata.School,SchoolsM.SchoolId);
        
        locb1 = locb1(locb1>0);
        
        temp = [];
        
        for i=1:size(nodes,1)
            
            temp = [temp;DistanceM(locb1(i),nodes(i)); ];
            
        end
        
        DataY(rdIndex,3) = temp;
    end
end

Moments.RDdata.DataY = DataY;

% RD regression

Moments.RDdata.Z = [Moments.RDdata.ThresCross ones(size(Moments.RDdata.DataY,1),1) Moments.RDdata.ScoreAc Moments.RDdata.ScoreLd Moments.RDdata.MiScoreLd Moments.RDdata.ScoreAc.*Moments.RDdata.ScoreLd Moments.RDdata.ScoreAc.*Moments.RDdata.MiScoreLd ];
Moments.RDdata.X = [Moments.RDdata.Treat ones(size(Moments.RDdata.DataY,1),1) Moments.RDdata.ScoreAc Moments.RDdata.ScoreLd Moments.RDdata.MiScoreLd Moments.RDdata.ScoreAc.*Moments.RDdata.ScoreLd Moments.RDdata.ScoreAc.*Moments.RDdata.MiScoreLd ];

for i = 1:6
     Moments.RDBeta(i,:) = ivregression(DataY(:,i),Moments.RDdata.X,Moments.RDdata.Z)';
end

Moments.WRDM = 1/6*eye(6);

%% Estimation

if Set.model == 'model 1'
    
    [Estimation.drawsN,Estimation.drawsW] = nwspgr('GQN',1,5);
    
    % Model 1 (8 parameters)
    
    Estimation.BetaEdu       = logical([ 1 0 0 0 0 0 0 0 ]);   % 1 parameter
    Estimation.BetaPoor      = logical([ 0 1 0 0 0 0 0 0 ]);   % 1 parameter
    Estimation.AlphaEdu      = logical([ 0 0 1 0 0 0 0 0 ]);   % 1 parameter
    Estimation.AlphaPoor     = logical([ 0 0 0 1 0 0 0 0 ]);   % 1 parameter
    Estimation.LambdaEdu     = logical([ 0 0 0 0 1 1 0 0 ]);   % 2 parameter
    Estimation.LambdaPoor    = logical([ 0 0 0 0 0 0 1 0 ]);   % 1 parameters
    Estimation.BetaRC        = logical([ 0 0 0 0 0 0 0 1 ]);   % 1 parameter
    
%     Estimation.ThetaMask{1,1} = [3,4];
%     Estimation.ThetaMask{2,1} = [1,3];
%     Estimation.ThetaMask{3,1} = [3,4];
%     Estimation.ThetaMask{4,1} = [1,3];
%     Estimation.ThetaMask{5,1} = [1,2];
%     Estimation.ThetaMask{6,1} = [3,4];
%     Estimation.ThetaMask{7,1} = [1,3];
%     Estimation.ThetaMask{8,1} = [1,2,3,4];

    Estimation.ThetaMask{1,1} = 0;
    Estimation.ThetaMask{2,1} = 0;
    Estimation.ThetaMask{3,1} = 0;
    Estimation.ThetaMask{4,1} = 0;
    Estimation.ThetaMask{5,1} = 0;
    Estimation.ThetaMask{6,1} = 0;
    Estimation.ThetaMask{7,1} = 0;
    Estimation.ThetaMask{8,1} = 1;
    
    Estimation.ThetaMask{1,2} = 1;
    Estimation.ThetaMask{2,2} = 1;
    Estimation.ThetaMask{3,2} = 2;
    Estimation.ThetaMask{4,2} = 2;
    Estimation.ThetaMask{5,2} = 3;
    Estimation.ThetaMask{6,2} = 3;
    Estimation.ThetaMask{7,2} = 3;
    Estimation.ThetaMask{8,2} = 1;

    Estimation.TypesTheta{1,1} = [2,4,5,7,8];
    Estimation.TypesTheta{2,1} = [5,8];
    Estimation.TypesTheta{3,1} = [1,2,3,4,6,8];
    Estimation.TypesTheta{4,1} = [1,3,6,8];

elseif Set.model == 'model 2'
    
    [Estimation.drawsN,Estimation.drawsW] = nwspgr('GQN',3,5);
    
    % Model 2 (17 parameters)

    Estimation.BetaQEdu      = logical([ 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameter
    Estimation.BetaQPoor     = logical([ 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameter
    Estimation.AlphaEdu      = logical([ 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameter
    Estimation.AlphaPoor     = logical([ 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]);   % 1 parameter
    Estimation.BetaZeEdu     = logical([ 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0]);   % 2 parameter
    Estimation.BetaZePoor    = logical([ 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]);   % 1 parameters
    Estimation.BetaZpEdu     = logical([ 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0]);   % 2 parameter
    Estimation.BetaZpPoor    = logical([ 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]);   % 1 parameters
    Estimation.LambdaEdu     = logical([ 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0]);   % 2 parameter
    Estimation.LambdaPoor    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]);   % 1 parameters
    Estimation.BetaQRC       = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]);   % 1 parameter
    Estimation.BetaZeRC      = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]);   % 1 parameter
    Estimation.BetaZpRC      = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]);   % 1 parameter
    Estimation.BetaZcorRC    = logical([ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]);   % 1 parameter
    
    Estimation.Theta1Names = {'BetaQ0';'Alpha0';'Private';'Charter';'ForProfit';'Religious';'Emblematic';['ChainFE (' num2str(size(Estimation.ChainFE,2)) ')'];['Market Year FE (' num2str(size(MarketID,1)) 'x' num2str(size(TimeID,1)) ')']};    
    Estimation.Theta1Names = {'BetaQEdu';'BetaQPoor';'AlphaEdu';'AlphaPoor';'BetaZeEdu';'BetaZeEdu';'BetaZePoor';'BetaZpEdu';'BetaZpEdu';'BetaZpPoor';'LambdaEdu';'LambdaEdu';'LambdaPoor';'BetaQRC';'BetaZeRC';'BetaZpRC';'BetaZcorRC'};
   
%     Estimation.ThetaMask{1,1}   = [3,4]; %1
%     Estimation.ThetaMask{2,1}   = [1,3];
%     Estimation.ThetaMask{3,1}   = [3,4];
%     Estimation.ThetaMask{4,1}   = [1,3];
%     Estimation.ThetaMask{5,1}   = [1,2]; %5
%     Estimation.ThetaMask{6,1}   = [3,4];
%     Estimation.ThetaMask{7,1}   = [1,3];
%     Estimation.ThetaMask{8,1}   = [1,2];
%     Estimation.ThetaMask{9,1}   = [3,4];
%     Estimation.ThetaMask{10,1}  = [1,3]; %10
%     Estimation.ThetaMask{11,1}  = [1,2];
%     Estimation.ThetaMask{12,1}  = [3,4];
%     Estimation.ThetaMask{13,1}  = [1,3];
%     Estimation.ThetaMask{14,1}  = [1,2,3,4];
%     Estimation.ThetaMask{15,1}  = [1,2,3,4]; % 15
%     Estimation.ThetaMask{16,1}  = [1,2,3,4];
%     Estimation.ThetaMask{17,1}  = [1,2,3,4];

    Estimation.ThetaMask{1,1}   = 0; %1
    Estimation.ThetaMask{2,1}   = 0;
    Estimation.ThetaMask{3,1}   = 0;
    Estimation.ThetaMask{4,1}   = 0;
    Estimation.ThetaMask{5,1}   = 0; %5
    Estimation.ThetaMask{6,1}   = 0;
    Estimation.ThetaMask{7,1}   = 0;
    Estimation.ThetaMask{8,1}   = 0;
    Estimation.ThetaMask{9,1}   = 0;
    Estimation.ThetaMask{10,1}  = 0; %10
    Estimation.ThetaMask{11,1}  = 0;
    Estimation.ThetaMask{12,1}  = 0;
    Estimation.ThetaMask{13,1}  = 0;
    Estimation.ThetaMask{14,1}  = 1;
    Estimation.ThetaMask{15,1}  = 1; % 15
    Estimation.ThetaMask{16,1}  = 1;
    Estimation.ThetaMask{17,1}  = 1;
    
    Estimation.ThetaMask{1,2}   = 1;
    Estimation.ThetaMask{2,2}   = 1;
    Estimation.ThetaMask{3,2}   = 2;
    Estimation.ThetaMask{4,2}   = 2;
    Estimation.ThetaMask{5,2}   = 4;
    Estimation.ThetaMask{6,2}   = 4;
    Estimation.ThetaMask{7,2}   = 4;
    Estimation.ThetaMask{8,2}   = 5;
    Estimation.ThetaMask{9,2}   = 5;
    Estimation.ThetaMask{10,2}  = 5;
    Estimation.ThetaMask{11,2}  = 3;
    Estimation.ThetaMask{12,2}  = 3;
    Estimation.ThetaMask{13,2}  = 3;
    Estimation.ThetaMask{14,2}  = 1;
    Estimation.ThetaMask{15,2}  = 4;
    Estimation.ThetaMask{16,2}  = 5;
    Estimation.ThetaMask{17,2}  = 5;
    
    Estimation.ThetaMask{14,3}  = 1;
    Estimation.ThetaMask{15,3}  = 2;
    Estimation.ThetaMask{16,3}  = 3;
    Estimation.ThetaMask{17,3}  = 2;

    Estimation.TypesTheta{1,1} = [2,4,5,7,8,10,11,13:17];
    Estimation.TypesTheta{2,1} = [5,8,11,14:17];
    Estimation.TypesTheta{3,1} = [1,2,3,4,6,7,9,10,12:17];
    Estimation.TypesTheta{4,1} = [1,3,6,9,12,14:17];
end


    Estimation.AdjustmentObjMM = 0;
    Estimation.AdjustmentObjIV = 0;

end
  
             
