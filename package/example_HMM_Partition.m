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

%This script calculates the average channel capacity between equal-partitioned data 
%within each stimulus versus the number of hidden states for hidden Markov
%model. The ratio between the numbers of hidden and emission states is
%fixed as 2.



%% Load data and experiment information

function filename=example_HMM_Partition(MainPara,OtherPara)

filefolder=pwd;

Dataset=MainPara.Dataset;%1 is our data set and TNFR_BK; 2 is p53 after Mdmx; 3 is Erk dataset; and 4 is synthetic data. 
%MainPara.Mutant=MainPara.Mutant;% Default is 0. For Dataset=1, 1 is WT or 0 is IkB-MainPara.Mutant dataset

if Dataset==0
    LabelName2=[MainPara.LabelName];%'User';
    
    
elseif Dataset==1
    LabelName2='NFkB';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
    if MainPara.Mutant==1
     LabelName2='NFkB_IkBMu';%'New3_HMM_TNFBk_13';%'New3_HMM_IkBMu_6';%: Ade's new full data with IkBMu and TNFBk
    end
elseif Dataset==2
    LabelName2='p53';%'p53_New_dataset_2';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
elseif Dataset==3
    LabelName2='Erk';%'p53_New_dataset_2';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
elseif Dataset==4
    TrajModes=1;LabelName2=['Synthetic_2_',num2str(TrajModes)];
end


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

OtherPara.TotalTimeLength=MainPara.TotalTimeLength;
Jackknife=1;% Default is 0. If 1: do Jackknife analysis. 
CalculateTrajProb=1;%Default is 1. If 1: Calculate trajectory probability etc for each trajectory
Selective=0;%Default is 0. If 1: select conditions for channel capacity, and 0 is not 
SelectedConditions=[12 11 6 17];%[6 17];%[6 ,12 11 17];%If OtherPara.CCSelect=0, the index of conditions used for channel capacity.
OtherPara.permuteData=permuteData;
OtherPara.permuteMode=permuteMode;
OtherPara.TimePointsToUse=TimePointsToUse;
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
OtherPara.CCSelect=Selective; %Default is 0. If 1: select conditions for channel capacity, and 0 is not 
if MainPara.Mutant==1
    SelectedConditions=[3 4];%
end
OtherPara.CCSelectCondi=SelectedConditions; %If OtherPara.CCSelect=0, the index of conditions used for channel capacity. Select the same as IkB-MainPara.Mutant
OtherPara.PointWiseOptim=1;%Default is 1. If 1: Calculate channel capacity by optimizing at each time point



OtherPara.sampleSize=10000;% If fitting each trajectory to a HMM, i.e., OtherPara.HMMFitToAll=0: the number of sampling to get best fit each traj from HMM
OtherPara.estTR_Num=50;%If fitting each trajectory to a HMM, i.e., OtherPara.HMMFitToAll=0: Number of HMM model to save for each condition

CorrelationToFeatures=0;%Default is 0. If 1: compare the features of data and sampled trajectories
OtherPara.Delayembedding=0;%Default is 0. use Delay embedding or not: 0 no, 1 and 2 differnt use


[X, condition, OtherPara,ID]=LoadDataset(OtherPara,MainPara);

OtherPara.PartialTrainTest=0;


LabelName2=[LabelName2,'_HMM_Partition'];
if MainPara.Model==2
LabelName2=[LabelName2,'_Markov'];
end
mkdir(LabelName2);% second number is ratio; default is 1;



for jjj=1:length(ID)%7:12%1:6%1:length(stim_ID)%:-1:1%10:17%1:9%14:20%13:-1:1%1:13%12:-1:1%13:18%19:26 %13:18%10:-1:1%11:17%18:26 %[26 2 8 1]%size(data.nfkb,2):-1:24]
%ID(1,jjj)=data.nfkb(jjj).ID;
%ID(1,jjj)=stim_ID(jjj);
RunSimulation=1;
%LabelName='New3_HMM_IkBMu_6';%'New3_HMM_TNFBk_13';%'New3_HMM_IkBMu_6';%: Ade's new full data with IkBMu and TNFBk
stim_ID2 =ID(1,jjj);
%Ratio=2.^[2 4 5 6 6.59]/32;%2.^[1 2 3 4 5 6]/2.^5;
Ratio=OtherPara.stateScan;%2.^[1:6];%[30 40 50 60]/32;%2.^[1 2 3 4 5 6]/2.^5;
%Ratio=[2/32 4/32];



for m=1:length(Ratio)
LabelName3=[LabelName2,'\Condition_',num2str(stim_ID2)];% 
%1 is bin 10:10:70, 01 is 2:2:8 



% % % OtherPara.binsize=32;%[2.^[1:6]];%10:10:80;%2:2:8 ;%:10:20;%70;
% % % OtherPara.state=round(Ratio(m)*OtherPara.binsize);%[2.^[1:6]];%10:10:80;%2:2:8 ;%:10:20;%70;
% % % %LabelName='New3_HMM_ParitialTrain_1';% is 650 
% % % 
% % % %'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
% % % %'New3_HMM_Condition_5_4_3';%:10TNF+100P3CSK4+control;
% % % %'New3_HMM_Condition_5_6_5';%:10TNF+10cpG+control; 'New3_HMM_Condition_16';%: Ade's new WT data
% % % %LabelName='New2_HMM_Condition_11';%Ade's old WT data
% % % %'New2_HMM_Condition_11';%'p53_Condition_2';%'New_HMM_Condition_11';%'p53_Condition_1';%'HMM_Condition_11_P1';
% % % %'ERK_Condition_1_OnePhase': condition 1 is fig6D, 2 is 6E lower, 3 is 6E upper; 
% % % %'ERK_Condition_1_TwoPhase_1': One phase is not separate time window, Two/Three is to separate before and after drug
% % % TimePointsToUse=143;
% % % ScaleRelation=1;
% % % OtherPara.PartialTrainTest=1;OtherPara.PartialNum=100;
% % % 
% % % 
% % % OtherPara.Delayembedding=0;%use Delay embedding or not: 0 no, 1 and 2 differnt use
% % % OtherPara.Dataset=Dataset; 
% % % OtherPara.HMM=1;%use Hidden Markov
% % % OtherPara.HMMFitToAll=1; % Fit to all traj to get a HMM and use it to sample
% % % OtherPara.HMMAllCondition=0; %Fit to all condition and use it to get mutual information later
% % % OtherPara.PartialTraining=0;   
% % % OtherPara.PartialRatio=1;%1/2;% , Try use half traj to train HMM
% % % OtherPara.ModelBasedStatistic=1;
% % % OtherPara.TrajBasedStatistic=1;
% % % OtherPara.CalculateStatistics=1;
% % % OtherPara.TrajBasedNormalize=0;%normalize
% % % OtherPara.CalculateDistance=1;
% % % OtherPara.EnhancedSample=2;


%CorrelationToFeatures=0;




if OtherPara.HMM==1
%     OtherPara.binsizeRatio=1;%0;%does not matter
%     OtherPara.stateRatio=Ratio.*OtherPara.binsizeRatio;%10-30,30-30,|30-50,60-50, |40-100, 100-100
    OtherPara.binsize=max(2,round(Ratio(m)/2));%0;%does not matter
    OtherPara.state=Ratio(m);%10-30,30-30,|30-50,60-50, |40-100, 100-100
   
    
    
    %disp(OtherPara.stateRatio);disp(OtherPara.binsizeRatio);
    OtherPara.sampleSize=10000;% number of sampling to get best fit traj for each HMM
    OtherPara.estTR_Num=50;%Number of HMM model to save for each condition
    LabelName=[LabelName3,'_States',num2str(OtherPara.state),'_bin_',num2str(OtherPara.binsize)];
    if OtherPara.PartialTraining==1
        LabelName=[LabelName,'Train_',num2str(1/OtherPara.PartialRatio)];
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
   
foldername=[filefolder,'\',LabelName];
filename=[filefolder,'\',LabelName,'.mat'];
figurefolder=[filefolder,'\Figures\Overfitting\',LabelName2];
mkdir(figurefolder);
  %OtherPara.figurenamehmmSample=[foldername,'\HMMTrajs_States_',num2str(OtherPara.stateRatio),'bin_',num2str(OtherPara.binsizeRatio)];

disp(filename);
%%
if exist(filename) 
      disp('load filename');
      load(filename); %OtherPara.TrajBasedNormalize=0;
   RunSimulation=0;
   
   Summary.CC(m,jjj)=mean(HMM.MI_Full);
Summary.HiddenStates(m,jjj)=OtherPara.state;
   
   
else
    %continue;
    RunSimulation=1;
%     if OtherPara.TrajBasedNormalize==1
%         filename2=[filefolder,'\',LabelName,'_normalized.mat'];
%         if exist(filename2) 
%             load(filename2); filename=filename2;
%             RunSimulation=0;
%         end
%     end  
end

%disp(OtherPara.TrajBasedNormalize);

% OtherPara.Dataset=Dataset; %1 is our data set; 2 is p53 after Mdmx; 3 is ERK; 
% OtherPara.ModelBasedStatistic=1;
% OtherPara.TrajBasedStatistic=1;   
% OtherPara.EnhancedSample=2;
% OtherPara.CalculateStatistics=1;
% OtherPara.CalculateDistance=1;
% OtherPara.Jacknife=0;
% OtherPara.CCSelect=0;
% OtherPara.PointWiseOptim=0;%Calculate channel capacity by optimizing at each time point
% if ScaleRelation==1
%        OtherPara.TrajBasedStatistic=0;
%        OtherPara.PointWiseOptim=1;
% end    

else
    foldername=[filefolder,'\',LabelName3];
filename=[filefolder,'\',LabelName3,'.mat'];
 

end

OtherPara.figurenamehmmSample=[foldername,'\HMMTrajs_States',num2str(OtherPara.state),'_bin_',num2str(OtherPara.binsize)];

if ~exist(foldername)
mkdir(foldername);
end
%% Pick experiment IDs
tic;
%stim_ID =[ 546 ,566 ,778  ,780 ,783 , 650];%





%%
if RunSimulation~=0
[X, condition, OtherPara,ID]=LoadDataset(OtherPara,MainPara);
    
Nq=2;%length(OtherPara.binsize);%length(X);% Y{1,j} is  j-th condition: with row time series, column trajs
%Nq=length(Ratio);%length(X);% Y{1,j} is  j-th condition: with row time series, column trajs
%Y=X;

% for i=jjj%1:size(condition,2)
%     OtherPara.condition2(i)=cellstr(condition{i});
% end

%for i=1:1%Nq
   
% % % %      A=X{jjj};
% % % %     
% % % % %     B=X2{jjj};
% % % %     if size(A,1)==1 % Unidimensional case: remove NaN individuals
% % % %         A=A(:,~isnan(A));
% % % %         A=A(:,~isinf(A));
% % % %        % B=B(:,~isnan(B));
% % % %        % B=B(:,~isinf(B));
% % % %     else % Multidimensional case: remove individuals with NaN in any dimension
% % % %         A=A(:,sum(isnan(A))==0);
% % % %         A=A(:,sum(isinf(A))==0); 
% % % %        % B=B(:,sum(isnan(B))==0);
% % % %        % B=B(:,sum(isinf(B))==0); 
% % % %     end
% % % %     A(A<0)=eps; %Make negative to be zero
% % % %    % B(B<0)=eps; %Make negative to be zero
% % % %     
% % % %     if OtherPara.Dataset==1
% % % %      Avector{i} = sort(A(:));
% % % %      AllValues{i}=Avector{i}(1:round(length(Avector{i})*0.99));%90% percentile of cells;
% % % %      MaxPercentile(i)=max(AllValues{i});%OtherPara.YUpLimit(i)=ceil(max(max(A)));
% % % %      
% % % %      A(A>10)=10; % To make MaxValue smaller for better training.
% % % %      OtherPara.YUpLimit(i)=ceil(MaxPercentile(i));
% % % %      OtherPara.YUpLimit(2)=ceil(MaxPercentile(i));
% % % %     end
% % % %     



    YTemp{1}=X{jjj};
%     OtherPara.MaxValue=max(OtherPara.MaxValue,max(max(YTemp{1})));
%     
%     OtherPara.MinValue=min(OtherPara.MinValue,max(min(min(YTemp{1})),0));
%     
    
   
    %if OtherPara.PartialTrainTest==1
    
        PartialNum=round(size(YTemp{1},2)/2);%Equal partition trajectories...
        Y{1}=YTemp{1}(:,end-PartialNum+1:end);
        Y{2}=YTemp{1}(:,1:end-PartialNum);    
        %ddd
    %end
%end
%Y2=Y;
% X=Y';
 %X2=Y';

%%


if OtherPara.HMM==1
   % Y{1}=Y{jjj};
   disp('Use hidden Markov model.');
   [SampledTraj,DistanceTotal,estTR{m},estE{m},estTR_sc{m},estE_sc{m},OtherPara,seqs] = HMMTrain(Y,OtherPara);
%[SampledTraj,DistanceTotal,estTR{m},estE{m},estTR_sc{m},estE_sc{m},OtherPara,seqs] = Partial_Train(Y,OtherPara); %HidenMarkov
%save(filename,'ID','DistanceTotal', 'SampledTraj','estTR','estE','estTR_sc','estE_sc','OtherPara','seqs');

else
%     [fitI,info1,I,Q,fX_all,XAverage,TrajEntropy,OtherPara,metrics,aux,PointsDelay]= jacknifeCC(X,5); 
%     d=size(X{1},1); 
    %save(filename,'ID','XAverage','TrajEntropy','fX_all','d','OtherPara','metrics','aux','PointsDelay');
%     if OtherPara.Delayembedding~=0
%         DelayEmbedding(foldername,PointsDelay, OtherPara);
%     end
    %cd(LabelName2);
   % disp('Use time-inhomogeneous Markov model.');
   % assignin('base', 'Y', Y);
    % now here....
    %[filename, HMM]=Markov(Dataset,OtherPara,Y,MainPara.Mutant,ID);
end
toc;



%%
if OtherPara.HMM==1
    
    
if OtherPara.HMMFitToAll==1
    SampledTraj=1;
else
  HMM=[];%disp(SampledTraj);
  
end

% if CorrelationToFeatures==1
% HMMCorrelationToFeatures(Y2,featureNames2,HMM,Nq,ID,OtherPara,foldername);
% ddd
% end

HMM.AA=0;
OtherPara.DistanceAsDistribution=1;%

% disp('Evaluating the training goodness of model.')
% HMM=HMMDistance(estTR{m},estE{m},OtherPara,seqs,Y,Nq,foldername,ID,HMM,SampledTraj);

disp('Calculating the time-dependent channel capacity.')
[DistanceConditions, HMM]=HMMStatistics(estTR{m},estE{m},OtherPara,seqs,Y,Nq,foldername,ID,HMM,SampledTraj);

%save(filename,'DistanceConditions','ID','DistanceTotal', 'SampledTraj','estTR','estE','estTR_sc','estE_sc','OtherPara','seqs','HMM');
%if ScaleRelation==0
 
%end

HMMPlotting(HMM,seqs,Nq,ID,OtherPara,foldername,SampledTraj);
MISummary=0;
close all;
save(filename,'DistanceConditions','ID','DistanceTotal', 'SampledTraj','estTR','estE','estTR_sc','estE_sc','OtherPara','seqs','HMM','MISummary');
end
end

A=HMM.MI_Full;
B = A(~isnan(A));

Summary.CC(m,jjj)=mean(B);
Summary.HiddenStates(m,jjj)=OtherPara.state;
end
end


Summary.CC=Summary.CC/2;% a rough scaling rescale

h2=figure ('position', [00, 10, 800, 600]);
%plot(Summary.HiddenStates(:,1),Summary.CC(m,jjj),'o','markersize',10); hold on;
%plot(Summary.HiddenStates(:,1),Summary.CC(m,jjj),'o','markersize',10); hold on;
errorbar(Summary.HiddenStates(:,1),mean(Summary.CC,2),std(Summary.CC,0,2),'o','markersize',10,'Color','b','MarkerFaceColor','b'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
cc=linspecer(jjj);
transparence=0.3;
% for kk=1:jjj
%     pp=plot(Summary.HiddenStates(:,1),Summary.CC(:,kk),'-','Color',cc(kk,:),'linewidth',1.5);
%     pp.Color(4) = transparence;
% end
xlim([0 max(Summary.HiddenStates(:,1))]);
ylim([0 3]);
%set(gca, 'xscale', 'log');
assignin('base', 'Summary', Summary);
xlabel('Number of hidden states');
ylabel('Max. mutual information (bits)');
set(gca,'FontSize',22,'linewidth',2);
figurenamehmm2=[figurefolder,'\CCVersusHiddenStates.jpg'];
saveas(gcf,figurenamehmm2);
filename3=[figurefolder,'\CCVersusHiddenStates.mat'];
save(filename3,'Summary');

close all;
end