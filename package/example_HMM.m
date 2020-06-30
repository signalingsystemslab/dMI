%(c) 2020 Signaling Systems Lab UCLA
%All rights reserved. 
%This MATLAB code package implements the quantification on the time-dependent
%channel capacity for data of time series. It used either
%time-inhomogeneous Markov model or hidden Markov model to learn the
%dynamical patterns of the set of time series. Then, it can reproduce the path
%ensemble by sampling a same amount of time series, and quantify the
%similarity between data and sampling. It further calculate the trajectory
%probability for each time series, and time-dependent channel capacity.

%A detailed description on the methods is given in the main text. A
%guideline for the package is in the README.txt file.

%This script is the main script to train a hidden Markov model, use it to
%sample, evaluate training goodness, and calculate channel-capacity.



%% Specify the dataset and name the folder for output
function filename=example_HMM(MainPara,OtherPara)
filefolder=pwd;
global ID

Dataset=MainPara.Dataset;%1 is our data set and TNFR_BK; 2 is p53 after Mdmx; 3 is Erk dataset; and 4 is synthetic data. 
%MainPara.Mutant=MainPara.Mutant;% Default is 0. For Dataset=1, 1 is WT or 0 is IkB-MainPara.Mutant dataset
if Dataset==0
    LabelName=[MainPara.LabelName];%'User';
    
    
elseif Dataset==1
    LabelName='NFkB';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
    if MainPara.Mutant==1
     LabelName='NFkB_IkBMu';%'New3_HMM_TNFBk_13';%'New3_HMM_IkBMu_6';%: Ade's new full data with IkBMu and TNFBk
    end
elseif Dataset==2
    LabelName='p53';%'p53_New_dataset_2';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
elseif Dataset==3
    LabelName='Erk';%'p53_New_dataset_2';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
elseif Dataset==4
    TrajModes=1;LabelName=['Synthetic_HMM_2_',num2str(TrajModes)];
end
OtherPara.LabelName=LabelName;

%% Parameter input
OtherPara.HMM=MainPara.Model;% Default is 1. If 1: use Hidden Markov
TimePointsToUse=MainPara.TimePointsToUse;%72; %Default is 150. The number of time points used 
permuteData=MainPara.permuteData;% Default is 0. If 1: permute ordering of time points; If 0: no permutation
permuteMode=MainPara.permuteMode;% If 1: descend, 2: ascend, 3: random permutation.
PartialTraining=MainPara.PartialTraining;%Default is 0. If 1: Use part of trajectory in each condition to train the model
% OtherPara.PartialRatio=1/2; %If OtherPara.PartialTraining=1:, use this ratio of the total number of trajectories to train HMM
% 
% OtherPara.binsize=4;%The number of emission states for hidden Markov model
% OtherPara.state=4;%128;%64;%The number of hidden states for hidden Markov model

TempCCSelect=MainPara.CCSelect;%Selective; %Default is 0. If 1: select conditions for channel capacity, and 0 is not 
TempCCSelectCondi=MainPara.CCSelectCondi;
Jackknife=1;% Default is 0. If 1: do Jackknife analysis. 
CalculateTrajProb=1;%Default is 1. If 1: Calculate trajectory probability etc for each trajectory
% Selective=0;%Default is 0. If 1: select conditions for channel capacity, and 0 is not 
%SelectedConditions=[12 11 6 17];%[6 17];%[6 ,12 11 17];%If OtherPara.CCSelect=0, the index of conditions used for channel capacity.
OtherPara.permuteData=permuteData;
OtherPara.permuteMode=permuteMode;
OtherPara.TimePointsToUse=TimePointsToUse;
OtherPara.TotalTimeLength=MainPara.TotalTimeLength;
OtherPara.Dataset=Dataset; %Default is 1. 1 is our data set; 2 is p53 after Mdmx; 3 is ERK; 
OtherPara.HMMFitToAll=1; %Default is 1. If 1: Fit to all traj to get a HMM and use it to sample
OtherPara.HMMAllCondition=0; %Default is 0. If 1: Fit to all condition and use it to get mutual information later
OtherPara.PartialTraining=PartialTraining; 
OtherPara.ModelBasedStatistic=1; %Default is 1. If 1: Calculate model-level trajectory entropy etc
OtherPara.TrajBasedStatistic=CalculateTrajProb;  %Default is 1. If 1: Calculate trajectory entropy etc for each trajectory
OtherPara.CalculateStatistics=1; %Default is 1. If 1: Calculate trajectory probability etc for each trajectory
%OtherPara.TrajBasedNormalize=0; %Default is 0. If 1: Normalize trajectory probability
OtherPara.CalculateDistance=1; %Default is 1. If 1: Calculate the distance between data and sampled trajectory
OtherPara.EnhancedSample=2; %Default is 2. The fold-increase on the number of sampled trajectoreis than data: sample more trajectories to plot heatmap for better visualization
OtherPara.Jackknife=Jackknife; %If 1: Do Jackknife analysis. 
OtherPara.JackknifeRatio=[1:-0.05:0.5]; %If OtherPara.Jackknife=1, the ratio of sub-sampling in Jackknife method
OtherPara.CCSelect=MainPara.CCSelect;%Selective; %Default is 0. If 1: select conditions for channel capacity, and 0 is not 
% if MainPara.Mutant==1
%     SelectedConditions=[3 4];%
% end
OtherPara.CCSelectCondi=MainPara.CCSelectCondi;%SelectedConditions; %If OtherPara.CCSelect=0, the index of conditions used for channel capacity. Select the same as IkB-MainPara.Mutant
OtherPara.PointWiseOptim=1;%Default is 1. If 1: Calculate channel capacity by optimizing at each time point



OtherPara.sampleSize=10000;% If fitting each trajectory to a HMM, i.e., OtherPara.HMMFitToAll=0: the number of sampling to get best fit each traj from HMM
OtherPara.estTR_Num=50;%If fitting each trajectory to a HMM, i.e., OtherPara.HMMFitToAll=0: Number of HMM model to save for each condition

CorrelationToFeatures=0;%Default is 0. If 1: compare the features of data and sampled trajectories
OtherPara.Delayembedding=0;%Default is 0. use Delay embedding or not: 0 no, 1 and 2 differnt use

%% Giving name to the simulation accoridng to the input parameters

if OtherPara.HMM==1  
    LabelName=[LabelName,'_HMM'];
    LabelName=[LabelName,'_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize)];
    if OtherPara.PartialTraining==1
        LabelName=[LabelName,'Partial_',num2str(1/OtherPara.PartialRatio)];
    end
    
    if OtherPara.HMMFitToAll==1
    LabelName=[LabelName];
    end


if permuteData==1
    LabelName=[LabelName,'_Permute_',num2str(permuteMode)];
elseif permuteData==2
    LabelName=[LabelName,'_Permute2_',num2str(permuteMode)];
end

if TimePointsToUse<140
    LabelName=[LabelName,'_Points_',num2str(TimePointsToUse)];
end
end


foldername=[filefolder,'\',LabelName];
filename=[filefolder,'\',LabelName,'.mat'];
OtherPara.figurenamehmmSample=[foldername,'\HMMTrajs_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize)];

if OtherPara.HMM==1
if ~exist(foldername)
mkdir(foldername);
end
end

%% Load the file if already generated and saved. Otherwise, run the model training and sampling, etc.
if exist(filename) 
%     if OtherPara.TrajBasedNormalize==1
%         filename2=[filefolder,'\',LabelName,'_normalized.mat'];
%         if exist(filename2) 
%             load(filename2); filename=filename2;
%         else
%             load(filename); OtherPara.TrajBasedNormalize=1;filename=filename2;
%         end
%     else
%         OtherPara.TrajBasedNormalize=0;
%     end
    disp('Have conducted the model training and mutual information calculation.');
   disp(['File name: ',filename]);
    load(filename);
    if exist('HMM')
    if isfield(HMM,'TrajProb') %If having calculated trajectory probability, jump the procedure.
        OtherPara.TrajBasedStatistic=0;
    end
    end
   RunSimulation=0;
   
else
    RunSimulation=1;
end
%  disp(filename);
%  ddd
% ddddddddd

% OtherPara.CCSelect=Selective; %Default is 0. If 1: select conditions for channel capacity, and 0 is not 
% OtherPara.CCSelectCondi=SelectedConditions; %If OtherPara.CCSelect=0, the index of conditions used for channel capacity. Select the same as IkB-MainPara.Mutant
OtherPara.CCSelect=TempCCSelect;%Selective; %Default is 0. If 1: select conditions for channel capacity, and 0 is not 
OtherPara.CCSelectCondi=TempCCSelectCondi;
% assignin('base', 'TempCCSelect', TempCCSelect)
% assignin('base', 'TempCCSelect2', TempCCSelectCondi)
%% Run simulation for the picked number of conditions

    
[X, condition, OtherPara,ID]=LoadDataset(OtherPara,MainPara);
%assignin('base', 'X', X([1 5],1));ddd

Y=X;
Nq=length(X);%ddd

if OtherPara.HMM==1 && OtherPara.Dataset==1
HistogramData(OtherPara.Avector,condition,foldername);
close all;
end

%% Vector method
%permutation is inside the function
if OtherPara.HMM==3
    disp('Use vector method.');
    [filenamehmm HMM]=VectorMethod(Y,OtherPara,MainPara);
    disp('Plot channel capacity!');
    PlotCC(filenamehmm,HMM,OtherPara);
    return;
elseif OtherPara.HMM==4
    disp('Use time-point method.');
    [filenamehmm HMM]=VectorMethod(Y,OtherPara,MainPara);
    disp('Plot channel capacity!');
    PlotCC(filenamehmm,HMM,OtherPara);
    return;
       
end

%% Permute data points if required by the user
if permuteData==1
    [index Y]=PermData(Y,permuteMode);
    %Y2=PermDataInputOrder(Y2,permuteMode,index);
    X=PermDataInputOrder(X,permuteMode,index);
    %X2=PermDataInputOrder(X2,permuteMode,index);  
elseif permuteData==2
    [index Y]=PermData2(Y,permuteMode);
    %Y2=PermDataInputOrder2(Y2,permuteMode,index);
    X=PermDataInputOrder2(X,permuteMode,index);
    %X2=PermDataInputOrder2(X2,permuteMode,index);  
end
    

if RunSimulation~=0
%% Run model training for hidden Markov model
%if RunSimulation~=0
if OtherPara.HMM==1
    disp('Use hidden Markov model.');
    %assignin('base', 'Y', Y);ddd
[SampledTraj,DistanceTotal,estTR,estE,estTR_sc,estE_sc,OtherPara,seqs] = HMMTrain(Y,OtherPara); %HidenMarkov
%save(filename,'ID','DistanceTotal', 'SampledTraj','estTR','estE','estTR_sc','estE_sc','OtherPara','seqs');
elseif OtherPara.HMM==2
    disp('Use time-inhomogeneous Markov model.');
    [filename, HMM]=Markov(Dataset,OtherPara,Y,MainPara.Mutant,ID);
    
    
% elseif OtherPara.HMM==3
%     disp('Use vector method.');
%     VectorMethod(Dataset,Y,OtherPara,MainPara)
%     return;
    
end


%% Analyze difference between sampled trajectories and data, and calculate trajectory probability, trajectory entropy and channel capacity

if OtherPara.HMM==1
    
    
if OtherPara.HMMFitToAll==1
    SampledTraj=1;
end

OtherPara.DistanceAsDistribution=1;% compare data and sample as distribution, rather than select nearest neightbor trajectories...

disp('Evaluating the training goodness of model.')
HMM=HMMDistance(estTR,estE,OtherPara,seqs,Y,Nq,foldername,ID,SampledTraj);


disp('Calculating the time-dependent channel capacity.')
[DistanceConditions, HMM]=HMMStatistics(estTR,estE,OtherPara,seqs,Y,Nq,foldername,ID,HMM,SampledTraj);
%assignin('base', 'HMM', HMM);ddd


%% Training goodness: likelihood

%% now here
if OtherPara.PartialTrainTest==1  
    display('Caluclate log-likelihood of test dataset.')
    YTest=cell(Nq,1);
    HMM.TrajProb2=cell(Nq,Nq);
    parfor i=1:Nq     
        YTest{i}=Y{i}(:,end-OtherPara.PartialNum+1:end);
        Y2{i}=Y{i}(:,1:end-OtherPara.PartialNum);    
        [SampledTraj2,DistanceTotal,estTR2,estE2,estTR_sc,estE_sc,OtherPara2,seqs2] = Partial_Train(Y2(i),OtherPara); %HidenMarkov
        [DistanceConditions2, HMM_testdata,seqs2]=Partial_Statistics(estTR2,estE2,OtherPara2,Y2(i),YTest(i),Nq,foldername,ID,HMM,SampledTraj2);
        TrajTemp{i}=HMM_testdata.TrajProb;% now here..
    end
    for i=1:Nq  
    HMM.TrajProb2(i,i)=TrajTemp{i};
    end
end



%% Save the simulation result

save(filename,'DistanceConditions','ID','DistanceTotal', 'SampledTraj','estTR','estE','estTR_sc','estE_sc','OtherPara','seqs','HMM');
% if OtherPara.Jackknife==1
%     continue;
% end

%% Plos results: preliminary figures
HMMPlotting(HMM,seqs,Nq,ID,OtherPara,foldername,SampledTraj);

% if CorrelationToFeatures==1
% disp('Comparing the features of the data and sampled trajectories.')
% HMMCorrelationToFeatures(Y,featureNames2,HMM,Nq,ID,OtherPara,foldername);
% end



end
end



if OtherPara.CCSelect==1
%     if RunSimulation==0
%         disp('Please first train the model for all the conditions, by setting MainPara.CCSelect=0;');
%     end
    if length(OtherPara.CCSelectCondi)>size(HMM.TrajProb,2)  || max(OtherPara.CCSelectCondi)>size(HMM.TrajProb,2) || length(OtherPara.CCSelectCondi)==1
        disp('The seletive conditions index is out of range.');
        return;
    else
    
    OtherPara.Jackknife=0;
    %disp(OtherPara.HMM);dddd
     disp('Calculating the time-dependent channel capacity of seletive conditions.')
     if OtherPara.HMM==1
        [DistanceConditions, HMM]=HMMStatistics(estTR,estE,OtherPara,seqs,Y,Nq,foldername,ID,HMM,SampledTraj);
     elseif OtherPara.HMM==2
              load(filename);%ddddd
         OtherPara.CCSelect=MainPara.CCSelect;%Selective; %Default is 0. If 1: select conditions for channel capacity, and 0 is not 
         OtherPara.CCSelectCondi=MainPara.CCSelectCondi;
         OtherPara.PartialTrainTest=0;
         
         OtherPara.Jackknife=0;
         [HMM]=MMStatisticsSelect(fX_all,OtherPara,Y,seqs,Nq,HMM);
         %[HMM,info1, Q,fX_all,XAverage, TrajEntropy,metrics,aux,PointsDelay,OtherPara] = getCCMarkov(Y,5,OtherPara);
     end
     
     
kk=1;Nowexist=1;
while Nowexist~=0
filename3=[foldername,'\CCofSelectedCondition_trial_',num2str(kk),'.mat'];
if exist(filename3)
    kk=kk+1;
else
    Nowexist=0;
end
end
CCSelected=HMM.MI_Full;%*HMM.IntersectionRatio;%have already scaled by Jackknife
disp(foldername);
save(filename3,'CCSelected','OtherPara');
    end

end


%end

if MainPara.Plotting==1
    disp('Plot channel capacity!');
    PlotCC(filename,HMM,OtherPara);
end

disp('DISC is done!')
close all;
end
