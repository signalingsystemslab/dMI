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
%within each stimulus versus the number of states for time-inhomogeneous Markov model.


function filename=Markov_partition(MainPara,OtherPara)
%% Load data and experiment information
%clear;
filefolder=pwd;
global ID

%% Parameter input
OtherPara.HMM=MainPara.Model;% Default is 1. If 1: use Hidden Markov
TimePointsToUse=MainPara.TimePointsToUse;%72; %Default is 150. The number of time points used 
permuteData=MainPara.permuteData;% Default is 0. If 1: permute ordering of time points; If 0: no permutation
permuteMode=MainPara.permuteMode;% If 1: descend, 2: ascend, 3: random permutation.
PartialTraining=MainPara.PartialTraining;%Default is 0. If 1: Use part of trajectory in each condition to train the model

Jackknife=1;% Default is 0. If 1: do Jackknife analysis. 
CalculateTrajProb=1;%Default is 1. If 1: Calculate trajectory probability etc for each trajectory
Selective=0;%Default is 0. If 1: select conditions for channel capacity, and 0 is not 
SelectedConditions=[12 11 6 17];%[6 17];%[6 ,12 11 17];%If OtherPara.CCSelect=0, the index of conditions used for channel capacity.
OtherPara.permuteData=permuteData;
OtherPara.permuteMode=permuteMode;
OtherPara.TotalTimeLength=MainPara.TotalTimeLength;
OtherPara.TimePointsToUse=TimePointsToUse;
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


StatesScan=OtherPara.stateScan;%2.^[1:6];
%Dataset=1;%1 is our data set; 2 is p53 after Mdmx;  3 is Erk

%MainPara.Dataset=Dataset;
OtherPara.PartialTrainTest=0;
% % if MainPara.Dataset==0
% %     %LabelName=OtherPara.LabelName;%'User';
% %    filefolder=[filefolder,'\',MainPara.LabelName];
% % elseif MainPara.Dataset==1 
% % filefolder=[filefolder,'\NFkB'];
% % elseif MainPara.Dataset==2 % 2 is p53 after Mdmx;
% %          filefolder=[filefolder,'\p53'];
% % elseif MainPara.Dataset==3 % 3 is Erk
% %        filefolder=[filefolder,'\Erk'];
% % end

Dataset=MainPara.Dataset;
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

 LabelName2=[LabelName2,'_Markov_Partition'];
 filefolder=[filefolder,'\',LabelName2];
     mkdir(filefolder);
     %disp(filefolder);ddd
figurefolder=[pwd,'\Figures\Overfitting\',LabelName2];
mkdir(figurefolder);

for m=1:length(StatesScan)%32;

    
    OtherPara.Dataset=MainPara.Dataset;
    
    MainPara.state=StatesScan(m);%disp(MainPara.state);
    OtherPara.state=StatesScan(m);
    
    %OtherPara.foldername=foldername;
[X, condition, OtherPara,ID]=LoadDataset(OtherPara,MainPara);
Y=X;

    LabelName=['Markov_States',num2str(MainPara.state)];
    
% %% Pre-processing on all the trajectories
% 
% Nq=length(X);% Y{1,j} is  j-th condition: with row time series, column trajs
% disp('number of conditions');disp(Nq);
% Y=cell(Nq,1);%Y2=cell(Nq,1);
% OtherPara.TimeWiseNum=size(X{1},1);% Train point with early time points, time-wise
% if TimePointsToUse<140
%     OtherPara.TimeWiseNum=TimePointsToUse;
% end
% 
% OtherPara.MaxValue=-Inf;OtherPara.MinValue=Inf;
% OtherPara.ID=ID;
% for i=1:size(condition,2)
%     OtherPara.condition2(i)=cellstr(condition{i});
% end
% 
% for i=1:Nq
%    
%      A=X{i};
%     if size(A,1)==1 % Unidimensional case: remove NaN individuals
%         A=A(:,~isnan(A));
%         A=A(:,~isinf(A));
%        % B=B(:,~isnan(A));
%         %B=B(:,~isinf(A));
%     else % Multidimensional case: remove individuals with NaN in any dimension
%         A=A(:,sum(isnan(A))==0);
%         A=A(:,sum(isinf(A))==0); 
%        % B=B(:,sum(isnan(A))==0);
%        % B=B(:,sum(isinf(A))==0); 
%     end
%     A(A<0)=eps; %Make negative to be zero
%     %B(B<0)=eps; %Make negative to be zero
%     
%     if OtherPara.Dataset==1 || OtherPara.Dataset==4
%      Avector{i} = sort(A(:));
%      AllValues{i}=Avector{i}(1:round(length(Avector{i})*0.99));%90% percentile of cells;
%      MaxPercentile(i)=max(AllValues{i});%OtherPara.YUpLimit(i)=ceil(max(max(A)));
%      
%      A(A>10)=10; % To make MaxValue smaller for better training.
%      OtherPara.YUpLimit(i)=ceil(MaxPercentile(i));
%     
%     end
%     
%     Y{i}=A;
%     OtherPara.MaxValue=max(OtherPara.MaxValue,max(max(Y{i})));
%     
%     OtherPara.MinValue=min(OtherPara.MinValue,max(min(min(Y{i})),0));   
% 
%     if size(Y{i},2)<OtherPara.PartialNum
%        OtherPara.PartialNum=size(Y{i},2);     
%     end
% end
% X=Y;


for jjj=1:length(X)
    LabelName2=[LabelName,'condition_',num2str(ID(jjj))];
    filename=[filefolder,'\',LabelName2,'.mat'];
    
    if exist(filename) 
      disp('load filename');
      load(filename); %OtherPara.TrajBasedNormalize=0;
       Summary.CC(m,jjj)=mean(HMM.MI_Full);
Summary.HiddenStates(m,jjj)=OtherPara.state;  
      
    else
        
        
    YTemp{1}=X{jjj};
    OtherPara.PartialNum=round(size(YTemp{1},2)/2);%Equal partition trajectories...
    Y2{1}=YTemp{1}(:,end-OtherPara.PartialNum+1:end);
    Y2{2}=YTemp{1}(:,1:end-OtherPara.PartialNum);  
   % assignin('base', 'Y2', Y2);
[HMM,fitI,info1,I,Q,fX_all,XAverage,TrajEntropy,estTR,estE,OtherPara,metrics,aux]= JackknifeCCMarkov(Y2,5,12,OtherPara); 


 

save(filename,'ID','XAverage','TrajEntropy','fX_all','OtherPara','metrics','aux','HMM');
end



Summary.CC(m,jjj)=mean(HMM.MI_Full);
Summary.HiddenStates(m,jjj)=OtherPara.state;
end
end


Summary.CC=Summary.CC;% a rough scaling rescale

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
