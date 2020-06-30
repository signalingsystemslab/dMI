%(c) 2020 Signaling Systems Lab UCLA
%All rights reserved. 
%This MATLAB code package implements the quantification on the time-dependent
%channel capacity for data of time series. It used either
%time-inhomogeneous Markov model or hidden Markov model to learn the
%dynamical patterns of the set of time series. Then, it can reproduce the path
%ensemble by sampling a same amount of time series, and quantify the
%similarity between data and sampling. It further calculate the trajectory
%probability for each time series, and time-dependent channel capacity.

%clear;
cd;
close all;
%% Specify the test dataset or user dataset
clearvars -except Dataset LabelName Data condition ID MainPara OtherPara; 
if exist('Dataset')==1
    MainPara.Dataset=Dataset;
else
    Dataset=0;% Comment later, 0 is test dataset; 1 is NFkB data set; 2 is p53 after Mdmx; 3 is Erk dataset; and 4 is synthetic data. 
    MainPara.Dataset=Dataset;
end

%% Load test dataset or user dataset
if MainPara.Dataset==0    
    if exist('LabelName')~=1
    LabelName='TestData'; %Please give a name to the simulation
    end
    
    if exist('Data')~=1
    Data='data1.mat'; %the file name of the dataset, 'data.mat'. The dataset should contain a set of numeric matrices of single-cell signaling responses, where rows are time points, and columns are single cells. It should be put into the same folder of the package main folder, otherwise please specify the folder path.
    end
    file=[pwd,'\',Data];%
    
    if exist(file)%~=0 && exist('condition')
        load(file);
        if exist('condition')~=1 || length(condition)>length(X)           
            temp1=1:1:length(X);
            condition=cellstr(string(temp1)); %the name label of each stimulus condition.
        end
        if exist('ID')==1
        OtherPara.UserID=ID;
        end
        OtherPara.UserCondition=condition;OtherPara.UserDataName=Data;
    else
        disp('Do not find the file name of data.')
        disp('Please provide the file name of data, or check the folder path of data, etc.');
        return;
    end
    MainPara.LabelName=LabelName;
else
    clearvars -except Dataset LabelName Data condition ID MainPara OtherPara;
    MainPara.Dataset=Dataset;OtherPara.Dataset=Dataset;
end


%% Specify simulation parameters
if ~isfield(MainPara,'Model')
MainPara.Model=1;
end
% Default is 1. 
%1: use Hidden Markov
%2: use Markov
%3: Vector Method
%4: Time-point Method

if ~isfield(OtherPara,'Scancase')
OtherPara.Scancase=0;
end
%0: single case; 
%1: HMM scan fix ratio between hidden and emission states or MM; 
%2: HMM vary hidden states; 
%3: HMM vary emission states;
%4: equal partition within conditions

if ~isfield(OtherPara,'binsize')
OtherPara.binsize=4;%The number of emission states for hidden Markov model
end
if ~isfield(OtherPara,'state')
OtherPara.state=4;%128;%64;%The number of hidden states for hidden Markov model
end

%% Parameters for perturbations
if isfield(MainPara,'CCSelect')~=1
MainPara.CCSelect=0;% First run the model with 0, and then use 1 to select conditions
end

if ~isfield(MainPara,'CCSelectCondi')
MainPara.CCSelectCondi=[1 2];
end


if ~isfield(MainPara,'Mutant')
MainPara.Mutant=0;% Default is 0. For Dataset=1, 1 is WT or 0 is IkB-mutant dataset
end

if ~isfield(MainPara,'TimePointsToUse')
MainPara.TimePointsToUse=143;%72; %Default is 150. The number of time points used 
end

if ~isfield(MainPara,'TotalTimeLength')
MainPara.TotalTimeLength=11.9;%in the unit of hours. Default is 11.9. The length of real time in the unit of hour
end

if ~isfield(MainPara,'permuteData')
MainPara.permuteData=0;% Default is 0. If 1: permute ordering of time points; If 0: no permutation
end
if ~isfield(MainPara,'permuteMode')
MainPara.permuteMode=3;% If 1: descend, 2: ascend, 3: random permutation.
end
if ~isfield(MainPara,'PartialTraining')
MainPara.PartialTraining=0;%Default is 0. If 1: Use part of trajectory in each condition to train the model
end
if ~isfield(OtherPara,'PartialRatio')
OtherPara.PartialRatio=1/2; %If OtherPara.PartialTraining=1:, use this ratio of the total number of trajectories to train HMM
end
OtherPara.PartialTrainTest=1;% For likelihood of test dataset
OtherPara.PartialNum=30;%Number of trajectory as test dataset

%% Setting for plotting
if ~isfield(MainPara,'PlotGoodness')
MainPara.PlotGoodness=0;% Plotting goodness of model training. Default is 1. If 1: plot goodness figures, and save them in a folder.
end
if ~isfield(MainPara,'Plotting')
    MainPara.Plotting=1;% Plotting maximum mutual information. Default is 1. If 1: plot maximum MI figure, and save output as mat and xlsx.
end

%% Setting for the vector method
if MainPara.Model==3
if ~isfield(MainPara,'WindowMode')
MainPara.WindowMode=1;% 1 is Naive scheme; 2 is Intermidate scheme: with using 5 time points
end
if ~isfield(MainPara,'WindowLength')
MainPara.WindowLength=5;
end
end

if MainPara.Model==3 || MainPara.Model==4
    OtherPara.Scancase=0;
end

%% Main script to train model
if OtherPara.Scancase==0
filename=example_HMM(MainPara,OtherPara);
end

%% Parameters for scanning the number of states
if ~isfield(OtherPara,'RatioScan')
OtherPara.RatioScan=2;
end
if ~isfield(OtherPara,'stateScan')
OtherPara.stateScan=2.^[1:6];
end
if ~isfield(OtherPara,'binsizeScan')
OtherPara.binsizeScan=2.^[1:6];
end
if ~isfield(OtherPara,'stateScanFix')
OtherPara.stateScanFix=64;
end
if ~isfield(OtherPara,'binsizeScanFix')
OtherPara.binsizeScanFix=32;
end



%% Conduct scan of states for training goodness or overfitting
if OtherPara.Scancase==4
   % for i=1:length(OtherPara.stateScan)
    if MainPara.Model==1 %HMM
        %OtherPara.state=OtherPara.stateScan(i);
       % OtherPara.binsize=OtherPara.binsizeScanFix;
    filename=example_HMM_Partition(MainPara,OtherPara);
    else
        %OtherPara.state=OtherPara.stateScan(i);
        %OtherPara.binsize=OtherPara.stateScan(i);
    filename=Markov_partition(MainPara,OtherPara);
    end
   % end
end


if OtherPara.Scancase==1 % HMM scan fix ratio between hidden and emission states or MM

    for i=1:length(OtherPara.stateScan)
        if MainPara.Model==1 %HMM
        OtherPara.state=OtherPara.stateScan(i);
        OtherPara.binsize=max(round(OtherPara.stateScan(i)/OtherPara.RatioScan),2);
        filename=example_HMM(MainPara,OtherPara);
       
        elseif MainPara.Model==2 %MM???
        OtherPara.state=OtherPara.stateScan(i);
        OtherPara.binsize=OtherPara.stateScan(i);
        filename=example_HMM(MainPara,OtherPara);
        end
    end
    
elseif OtherPara.Scancase==2 % HMM vary hidden states
    if MainPara.Model==1 %HMM
    for i=1:length(OtherPara.stateScan)
         OtherPara.state=OtherPara.stateScan(i);
        OtherPara.binsize=OtherPara.binsizeScanFix;
        filename=example_HMM(MainPara,OtherPara);
    end
    end
elseif OtherPara.Scancase==3 % HMM vary emission states
    if MainPara.Model==1 %HMM
    for i=1:length(OtherPara.binsizeScan)
        OtherPara.state=OtherPara.stateScanFix;
        OtherPara.binsize=OtherPara.binsizeScan(i);
        filename=example_HMM(MainPara,OtherPara);
    end
    end
    
end

%% Plotting goodness
if OtherPara.Scancase==1 ||  OtherPara.Scancase==2 ||  OtherPara.Scancase==3% MainPara.PlotGoodness==1
Plot_Goodness(MainPara,OtherPara);
end
