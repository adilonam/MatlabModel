%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DEMAND ESTIMATION LAB %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

Project:        --
Author:         --
Date created:   May 23rd, 2018
Date modified:  May 24th, 2018
Description:    This code creates fake data and estimates demand model     

%}

%==============================================================================
% Step 0: Paths and globals
%==============================================================================

clear all; clc;

temp=extractBetween(pwd,1,'Desktop');

if strcmp(temp,'your path') == 1 % Paths Claudia Personal
    
    stemGit = [ 'your path' '/CodeUpwork/MatlabModel']; 
   
end


clear temp

global Set
%datadate='2019_05_23';
datadate='Fake_Test';
savedate='2019_05_23';
cd(stemGit);

%==============================================================================
% Step 1: Setup the Basic data
%==============================================================================

% Choose Markets, Years, Model and Moments

%{

Markets:
   
    - [xx xx]:  Any selection (markets go from 1 to 102)
    - [991]:    Small markets (under 100)  
    - [992]:    Medium markets (under 300)  
    - [993]:    All markets except Lima (101) 
    - []:       All markets including Lima (102)

Years:
    - 2014-2018

Models:

    - 'model 1': {p,q,d}, RC only on Q
    - 'model 2': {p,q,d,z}, RC on Q and Ze, Zy

ModelPeers:

    - 'basic'
    - 'leave one out'

ModelRD_Moments:

    - 'min distance'
    - 'moment condition'

Moments:

    - 'moments 1': Micro,IV
    - 'moments 2': Micro,IV,RD
    - 'moments 3': Micro,IV,RD,MLE
%}

% temp

%run TestData.m

Set.markets     =   'small';

%Set.markets     =   'all';
%Set.years       =   [2014:2016 2018];
Set.years       =   [2014 2018];
Set.model       =   'model 2';
Set.peersmodel  =   'leave one out';
Set.momentsrd   =   'moment condition';
Set.moments     =   'moments 2';
Set.method      =   'Squarem';
Set.datadate    =   datadate;

[Schools,Consumers,Cweights,Distance,Moments,Estimation,marketslist]   =   setupData(Set,stemWorked);

Set.marketslist=marketslist;


Set.mW = 0.5;
Set.tol_inner = 10^-12;
Set.tol_outer = 10^-12; 
Estimation.W_IV = 1/3;
Estimation.W_MM = 1/3;
Estimation.W_RDM = 1/3;

if Set.model == 'model 1'

    Set.Th2_0 = [0.5 -0.2 -0.3 0.2 -0.2 0.3 0.5 0.3]';
    
elseif Set.model == 'model 2'
    
    Set.Th2_0 = [0.5 -0.2, -0.3 0.2, 0.8 -0.3 -0.6 0.3, -0.2 0.3, 0.5 0.4 0.6 0.3, 0.1 0.4 0.2]';

end

Theta2 = Set.Th2_0;

Delta0 = ones(size(Schools,1),1);
Schools.Share =zeros (size(Schools.Share));

Params = GetParams(Estimation,Set,Theta2);


for m = marketslist'
    for y = Set.years
        DeltaM = Schools.Delta(Schools.MarketId==m & Schools.Year==y);
        SchoolsM = Schools(Schools.MarketId==m & Schools.Year==y,:);
        DistanceM = Distance(Schools.MarketId==m & Schools.Year==y,Consumers.MarketId==m);
        CweightsMAll = squeeze(Cweights.All(Consumers.MarketId==m,2,:));
        CweightsMTypes = squeeze(Cweights.Types(Consumers.MarketId==m,2,:));
        Schools.Share(Schools.MarketId==m & Schools.Year==y) = RC_shares(DeltaM,SchoolsM,DistanceM,CweightsMAll,CweightsMTypes,Estimation,Set,Params);
      
    end
end


%==============================================================================
% Step 2: Check Derivatives and Time Main Functions
%==============================================================================

%%%%%  1. Time Delta and Check Derivatives of (s) wrt Delta %%%%%%%%%%%%%%%%%%%

% Set.method      =   'Squarem';
% 
% tic;
% delta = SolveShares(Delta0,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2);
% Time(1)=toc;
% 
% tic;
% [delta,Jac] = SolveShares(Delta0,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2);
% Time(2)=toc;

% Set.method      =   'Newton';
% 
% tic;
% delta = SolveShares(Delta0,Schools,Consumers,Distance,Cweights,Estimation,Set,Theta2);
% Time(2)=toc;
% 
% fprintf('\n Time Squarem (seconds) = %5.3g \n Time Newton (seconds) = %5.3g  \n',Time)


%%%%%% 2. Check Derivatives of (s) wrt Delta (use 1 market) %%%%%

Supply=0;
Set.step=10^-6; 

diff = DerivativesCheck(Schools,Consumers,Distance,Cweights,Moments,Estimation,Set,Theta2,6,Set.years(1),Supply);

