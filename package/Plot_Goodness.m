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

%This script plots the quantification on the training goodness of 
%hidden Markov model, including false nearest neighbor probability, 
%relative KL-divergence, and rescaled likelihood. It was used for our
%NFkB dataset, p53 and Erk datasets. Here, the ratio between the number of
%hidden and emission states is fixed, given by Ratio.

function Plot_Goodness(MainPara,OtherPara)
%% Specify the dataset and name the folder for output

filefolder=pwd;
OtherPara.RatioScan=2;% The ratio between the number of hidden and emission states.
Dataset=MainPara.Dataset;%1 is our data set and TNFR_BK; 2 is p53 after Mdmx; 3 is Erk dataset; and 4 is synthetic data. 

if Dataset==0
    LabelName=MainPara.LabelName;%'User';
    
    
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


% % if OtherPara.HMM==1    
% %     LabelName=[LabelName,'_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize)];
% %     if OtherPara.PartialTraining==1
% %         LabelName=[LabelName,'Partial_',num2str(1/OtherPara.PartialRatio)];
% %     end
% %     
% %     if OtherPara.HMMFitToAll==1
% %     LabelName=[LabelName];
% %     end
% % 
% % 
% % if permuteData==1
% %     LabelName=[LabelName,'_Permute_',num2str(permuteMode)];
% % elseif permuteData==2
% %     LabelName=[LabelName,'_Permute2_',num2str(permuteMode)];
% % end
% % 
% % if TimePointsToUse<140
% %     LabelName=[LabelName,'_Points_',num2str(TimePointsToUse)];
% % end
% % end

OtherPara.HMM=MainPara.Model;% Default is 1. If 1: use Hidden Markov
TimePointsToUse=MainPara.TimePointsToUse;%72; %Default is 150. The number of time points used 
permuteData=MainPara.permuteData;% Default is 0. If 1: permute ordering of time points; If 0: no permutation
permuteMode=MainPara.permuteMode;% If 1: descend, 2: ascend, 3: random permutation.
PartialTraining=MainPara.PartialTraining;%Default is 0. If 1: Use part of trajectory in each condition to train the model
% OtherPara.PartialRatio=1/2; %If OtherPara.PartialTraining=1:, use this ratio of the total number of trajectories to train HMM
% 
% OtherPara.binsize=4;%The number of emission states for hidden Markov model
% OtherPara.state=4;%128;%64;%The number of hidden states for hidden Markov model


Jackknife=1;% Default is 0. If 1: do Jackknife analysis. 
CalculateTrajProb=1;%Default is 1. If 1: Calculate trajectory probability etc for each trajectory
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

if MainPara.Model==1 %HMM
if min(abs(OtherPara.Scancase-[1,2,3]))>0
   disp('Please re-specify the parameter OtherPara.Scancase to 1 or 2 or 3');
   return;
end
 LabelName=[LabelName,'_HMM'];
if OtherPara.Scancase==1
        stateScan=OtherPara.stateScan;
        binsizeScan=max(round(OtherPara.stateScan/OtherPara.RatioScan),2);
        figurefolder=[filefolder,'\Figures\Goodness\',LabelName,'_FixStatesRatio_',num2str(OtherPara.RatioScan)];
        XaxisState=stateScan;
elseif OtherPara.Scancase==2 % HMM vary hidden states
        stateScan=OtherPara.stateScan;
        binsizeScan=OtherPara.binsizeScanFix*ones(1,length(stateScan));
        figurefolder=[filefolder,'\Figures\Goodness\',LabelName,'_ScanHiddenFixEmi_',num2str(OtherPara.binsizeScanFix)];
        XaxisState=stateScan;
elseif OtherPara.Scancase==3 % HMM vary emission states        
        binsizeScan=OtherPara.binsizeScan;
        stateScan=OtherPara.stateScanFix*ones(1,length(binsizeScan));
        XaxisState=binsizeScan;
        figurefolder=[filefolder,'\Figures\Goodness\',LabelName,'_ScanEmiFixHidden_',num2str(OtherPara.stateScanFix)];
end
    
end

if MainPara.Model==2 %MM???
    if min(abs(OtherPara.Scancase-[1]))>0
   disp('To plot goodness, please re-specify the parameter OtherPara.Scancase to 1');
   return;
    end
 
%     if Dataset==1
%     LabelName='Markov_NFkB';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
%     if Mutant==1
%      LabelName='Markov_NFkB_IkBMu';%'New3_HMM_TNFBk_13';%'New3_HMM_IkBMu_6';%: Ade's new full data with IkBMu and TNFBk
%     end
% elseif Dataset==2
%     LabelName='Markov_p53';%'p53_New_dataset_2';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
% elseif Dataset==3
%     LabelName='Markov_Erk';%'p53_New_dataset_2';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
%     end
LabelName=[LabelName,'_Markov'];%LabelName=['Markov_',LabelName];
stateScan=OtherPara.stateScan;
binsizeScan=OtherPara.stateScan;
figurefolder=[filefolder,'\Figures\Goodness\',LabelName];

XaxisState=stateScan;
end
 

if MainPara.Model~=2 && MainPara.Model~=1
    disp('Change MainPara.Model to be 1 or 2.')
    return;
end
       
mkdir(figurefolder)
LabelNameTemp=LabelName;
for jjj=1:length(stateScan)%Scan state size
    
if OtherPara.HMM==1   
    
    OtherPara.sampleSize=10000;% number of sampling to get best fit traj for each HMM
    OtherPara.estTR_Num=50;%Number of HMM model to save for each condition
    LabelName=[LabelNameTemp,'_States_',num2str(stateScan(jjj)),'bin_',num2str(binsizeScan(jjj))];    
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



elseif OtherPara.HMM~=1
    
    disp(stateScan(jjj));
    LabelName=[LabelNameTemp,'_States_',num2str(stateScan(jjj))];
if OtherPara.PartialTraining==1
    LabelName=[LabelName,'Partial_',num2str(1/OtherPara.PartialRatio)];
end

if OtherPara.HMMFitToAll==1
LabelName=[LabelName];
end
if OtherPara.permuteData==1
    LabelName=[LabelName,'_Permute_',num2str(OtherPara.permuteMode)];
elseif OtherPara.permuteData==2
    LabelName=[LabelName,'_Permute2_',num2str(OtherPara.permuteMode)];
end

if OtherPara.TimePointsToUse<140
    LabelName=[LabelName,'_Points_',num2str(OtherPara.TimePointsToUse)];
end


end


foldername=[filefolder,'\',LabelName];
filename=[filefolder,'\',LabelName,'.mat'];
OtherPara.figurenamehmmSample=[foldername,'\HMMTrajs_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize)];


%disp(filename);
%RatioLabel=['HMM_ScanStates_Ratio',num2str(OtherPara.RatioScan)];



disp(filename);%ddddddddd
%%
if exist(filename) 
    %disp(jjj);
    load(filename);

    ConditionIndex=[];
    Nq=size(HMM.TrajProb,2);%length(OtherPara.binsize);
    %binSizeScale=cell(1,Nq);binSizeScale2=cell(1,Nq);
 NumofPoints=100; %sum up first NumofPoints to avoid inifinte of late time points
for Id_index=1:Nq% conditions
    for jj=1:size(HMM.TrajProb{Id_index,Id_index},2) %Time points for rescaled likelihood
       % assignin('base', 'HMM', HMM);
    binSizeScale{Id_index}(:,jj)=log(mean(HMM.TrajProb{Id_index,Id_index}(:,jj),1))+jj*log(binsizeScan(jjj));
    binSizeScale2{Id_index}(:,jj)=jj*log(binsizeScan(jjj));
    end
    if isfield(HMM,'TrajProb2') %exist(HMM.TrajProb2)
    for jj=1:size(HMM.TrajProb2{Id_index,Id_index},2)
    binSizeScale2{Id_index}(:,jj)=log(mean(HMM.TrajProb2{Id_index,Id_index}(:,jj),1))+jj*log(binsizeScan(jjj));
    end
    Likelihood_TestData(jjj,Id_index)=mean(binSizeScale2{Id_index})/size(HMM.TrajProb2{Id_index,Id_index},2);% renoramlize the number back
    AIC(jjj,Id_index)=Likelihood_TestData(jjj,Id_index)-Id_index^2/100;
    end
    
    binsizeSummary(jjj,Id_index)=binsizeScan(jjj); 
    Likelihood(jjj,Id_index)=mean(binSizeScale{Id_index})/size(HMM.TrajProb{Id_index,Id_index},2);% renoramlize the number back
    FalseNeighbor1(jjj,Id_index)=HMM.FalseNeighborProb1(Id_index); 
    FalseNeighbor2(jjj,Id_index)=HMM.FalseNeighborProb12(Id_index);% renoramlize the number back
    KLTemp1=HMM.KLDivergence{Id_index};
    KLTemp1(find(KLTemp1==Inf))=[];
    KLTemp2=HMM.KLEntropy{Id_index};
    KLTemp2(find(KLTemp2==Inf))=[];
    KLDivergence(jjj,Id_index)=mean(KLTemp1);
    KLEntropy(jjj,Id_index)=mean(KLTemp2);
end
  

  ConditionIndex=[ConditionIndex,Id_index];
end
end

%ddd


%%
transparence=0.6;
cc=linspecer(Nq);


 h1=figure ('position', [00, 10, 1700, 1000]);
 Temp1=[];
for Id_index=1:Nq
   % if min(abs(Id_index-ConditionIndex))==0
        disp(Id_index);
pp=plot(XaxisState,Likelihood(:,Id_index),'-*','linewidth',3,'Color',cc(Id_index,:));hold on;
pp.Color(4) = transparence;
% hMarkers = pp.MarkerHandle;
% hMarkers.FaceColorData(4)= transparence;
    Temp1=[Temp1,Likelihood(:,Id_index)];
   % end
     %p1.Color(4) = 0.9;
end
%plot(binsizeSummary(jjj,:),Likelihood(jjj,:),'-*','linewidth',3,'Color',cc(jjj,:));hold on;

b=errorbar(XaxisState,mean(Temp1,2),std(Temp1',1),'o','markersize',20,'linewidth',3,'Color','b','MarkerFaceColor','b'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
b.CapSize = 12;
h=legend(OtherPara.condition2,'Location','bestoutside');
%title(h,'States number');
set(h,'FontSize',26,'fontWeight','normal');
%set(gca,'XLim',[0 max(XaxisState)],'Xtick',2.^[1:6],'Xticklabel',2.^[1:6]) ; 
set(gca,'XLim',[0 max(XaxisState)],'Xtick',[0:round(max(XaxisState)/4):max(XaxisState)],'Xticklabel',[0:round(max(XaxisState)/4):max(XaxisState)]) ; 
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 limsy(2)]);
%set(gca,'xscale','log');
xlabel('Number of hidden states');
ylabel('Rescaled log-likeli. of training set');
set(gca,'FontSize',40,'linewidth',2);
figurenamehmm=[figurefolder,'\Likelihood_TrainingData_Condi_',num2str(Nq),'States_',num2str(stateScan(1)),'bin_',num2str(binsizeScan(1)),'.jpg'];
print(gcf, '-djpeg', '-r400',figurenamehmm)
%[B1,I1] = sort(Likelihood(:,DeleteCondiNumStatesIndex));
% for kk=1:length(I1)
% BadCondition1(kk)=data.nfkb(I1(kk)).ID;
% end



 if isfield(HMM,'TrajProb2')
     h1=figure ('position', [00, 10, 1700, 1000]);
     Temp1=[];
for Id_index=1:Nq
    %if min(abs(Id_index-ConditionIndex))==0
        disp(Id_index);
p1=plot(XaxisState,Likelihood_TestData(:,Id_index),'-*','linewidth',3,'Color',cc(Id_index,:));hold on;
Temp1=[Temp1,Likelihood_TestData(:,Id_index)];
   % end
     %p1.Color(4) = 0.9;
end
assignin('base', 'Temp1', Temp1);
b=errorbar(XaxisState,mean(Temp1,2),std(Temp1',1),'o','markersize',20,'linewidth',3,'Color','b','MarkerFaceColor','b'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
b.CapSize = 12;
h=legend(OtherPara.condition2,'Location','bestoutside');
%title(h,'States number');
set(h,'FontSize',26,'fontWeight','normal');
set(gca,'XLim',[0 max(XaxisState)],'Xtick',[0:round(max(XaxisState)/4):max(XaxisState)],'Xticklabel',[0:round(max(XaxisState)/4):max(XaxisState)]) ; 
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 limsy(2)]);
%set(gca,'xscale','log');
xlabel('Number of hidden states');
ylabel('Rescaled log-likelihood of test set');
set(gca,'FontSize',40,'linewidth',2);
figurenamehmm=[figurefolder,'\Likelihood_TestData_Condi_',num2str(Nq),'States_',num2str(stateScan(1)),'bin_',num2str(binsizeScan(1)),'.jpg'];
print(gcf, '-djpeg', '-r400',figurenamehmm)
 end
 

 
 
% % %   h1=figure ('position', [00, 10, 1700, 1000]);
% % %      Temp1=[];
% % % for Id_index=1:Nq
% % %     %if min(abs(Id_index-ConditionIndex))==0
% % %         disp(Id_index);
% % % p1=plot(XaxisState,AIC(:,Id_index),'-*','linewidth',3,'Color',cc(Id_index,:));hold on;
% % % Temp1=[Temp1,AIC(:,Id_index)];
% % %     %end
% % %      %p1.Color(4) = 0.9;
% % % end
% % % b=errorbar(XaxisState,mean(Temp1,2),std(Temp1',1),'o','markersize',20,'linewidth',3,'Color','b','MarkerFaceColor','b'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
% % % b.CapSize = 12;
% % % h=legend(OtherPara.condition2,'Location','bestoutside');
% % % %title(h,'States number');
% % % set(h,'FontSize',26,'fontWeight','normal');
% % % set(gca,'XLim',[0 max(XaxisState)],'Xtick',[0:round(max(XaxisState)/4):max(XaxisState)],'Xticklabel',[0:round(max(XaxisState)/4):max(XaxisState)]) ; 
% % % limsy=get(gca,'YLim');
% % % set(gca,'Ylim',[0 limsy(2)]);
% % % %set(gca,'xscale','log');
% % % xlabel('Number of hidden states');
% % % ylabel('AIC likelihood of test data');
% % % set(gca,'FontSize',40,'linewidth',2);
% % % figurenamehmm=[figurefolder,'\AIC_sampled_Condi_',num2str(Nq),'States_',num2str(stateScan(1)),'bin_',num2str(binsizeScan(1)),'.jpg'];
% % % print(gcf, '-djpeg', '-r400',figurenamehmm)





 h1=figure ('position', [00, 10, 1700, 1000]);
 Temp1=[];
for Id_index=1:Nq
    %if min(abs(Id_index-ConditionIndex))==0
       % disp(jjj);
        FalseNeighborAverage(:,Id_index)=(FalseNeighbor1(:,Id_index)+FalseNeighbor2(:,Id_index))/2;
pp=plot(XaxisState,FalseNeighborAverage(:,Id_index),'-*','linewidth',3,'Color',cc(Id_index,:));hold on;
pp.Color(4) = transparence;
    Temp1=[Temp1,FalseNeighborAverage(:,Id_index)];
   % end
end
b=errorbar(XaxisState,mean(Temp1,2),std(Temp1',1),'o','markersize',20,'linewidth',3,'Color','b','MarkerFaceColor','b'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
b.CapSize = 12;
h=legend(OtherPara.condition2,'Location','bestoutside');
%title(h,'States number');
set(h,'FontSize',26,'fontWeight','normal');
%set(gca,'XLim',[0 max(XaxisState)],'Xtick',2.^[1:6],'Xticklabel',2.^[1:6]) ; 
set(gca,'XLim',[0 max(XaxisState)],'Xtick',[0:round(max(XaxisState)/4):max(XaxisState)],'Xticklabel',[0:round(max(XaxisState)/4):max(XaxisState)]) ; 
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 0.5]);
%set(gca,'xscale','log');
xlabel('Number of hidden states');
ylabel('False k-nearest neigh. prob.');
set(gca,'FontSize',40,'linewidth',2);
figurenamehmm=[figurefolder,'\FalseNearestNeighbor_k',num2str(HMM.kNNnum),'Condi_',num2str(Nq),'States_',num2str(stateScan(1)),'bin_',num2str(binsizeScan(1)),'.jpg'];
print(gcf, '-djpeg', '-r400',figurenamehmm)
%[B2,I2] = sort(FalseNeighborAverage(:,DeleteCondiNumStatesIndex));
% for kk=1:length(I2)
% BadCondition2(kk)=data.nfkb(I2(kk)).ID;
% end


% % %  h1=figure ('position', [00, 10, 1700, 1000]);
% % %  Temp1=[];
% % % for Id_index=1:Nq
% % %     %if min(abs(Id_index-ConditionIndex))==0
% % %        % disp(jjj);
% % %         FalseNeighborAverage(:,Id_index)=FalseNeighbor1(:,Id_index);
% % % pp=plot(XaxisState,FalseNeighborAverage(:,Id_index),'-*','linewidth',3,'Color',cc(Id_index,:));hold on;
% % % pp.Color(4) = transparence;
% % %     Temp1=[Temp1,FalseNeighborAverage(:,Id_index)];
% % %    % end
% % % end
% % % b=errorbar(XaxisState,mean(Temp1,2),std(Temp1',1),'o','markersize',20,'linewidth',3,'Color','b','MarkerFaceColor','b'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
% % % b.CapSize = 12;
% % % h=legend(OtherPara.condition2,'Location','bestoutside');
% % % %title(h,'States number');
% % % set(h,'FontSize',26,'fontWeight','normal');
% % % %set(gca,'XLim',[0 max(XaxisState)],'Xtick',2.^[1:6],'Xticklabel',2.^[1:6]) ; 
% % % set(gca,'XLim',[0 max(XaxisState)],'Xtick',[0:round(max(XaxisState)/4):max(XaxisState)],'Xticklabel',[0:round(max(XaxisState)/4):max(XaxisState)]) ; 
% % % limsy=get(gca,'YLim');
% % % set(gca,'Ylim',[0 0.5]);
% % % %set(gca,'xscale','log');
% % % xlabel('Number of hidden states');
% % % ylabel('False k-nearest neigh. prob.');
% % % set(gca,'FontSize',40,'linewidth',2);
% % % figurenamehmm=[figurefolder,'\FalseNearestNeighborSToD_k',num2str(HMM.kNNnum),'Condi_',num2str(Nq),'States_',num2str(stateScan(1)),'bin_',num2str(binsizeScan(1)),'.jpg'];
% % % print(gcf, '-djpeg', '-r400',figurenamehmm)
% % % 
% % % 
% % % h1=figure ('position', [00, 10, 1700, 1000]);
% % %  Temp1=[];
% % % for Id_index=1:Nq
% % %    % if min(abs(Id_index-ConditionIndex))==0
% % %        % disp(jjj);
% % %         FalseNeighborAverage(:,Id_index)=FalseNeighbor2(:,Id_index);
% % % pp=plot(XaxisState,FalseNeighborAverage(:,Id_index),'-*','linewidth',3,'Color',cc(Id_index,:));hold on;
% % % pp.Color(4) = transparence;
% % %     Temp1=[Temp1,FalseNeighborAverage(:,Id_index)];
% % %    % end
% % % end
% % % b=errorbar(XaxisState,mean(Temp1,2),std(Temp1',1),'o','markersize',20,'linewidth',3,'Color','b','MarkerFaceColor','b'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
% % % b.CapSize = 12;
% % % h=legend(OtherPara.condition2,'Location','bestoutside');
% % % %title(h,'States number');
% % % set(h,'FontSize',26,'fontWeight','normal');
% % % %set(gca,'XLim',[0 max(XaxisState)],'Xtick',2.^[1:6],'Xticklabel',2.^[1:6]) ; 
% % % set(gca,'XLim',[0 max(XaxisState)],'Xtick',[0:round(max(XaxisState)/4):max(XaxisState)],'Xticklabel',[0:round(max(XaxisState)/4):max(XaxisState)]) ; 
% % % limsy=get(gca,'YLim');
% % % set(gca,'Ylim',[0 0.5]);
% % % %set(gca,'xscale','log');
% % % xlabel('Number of hidden states');
% % % ylabel('False k-nearest neigh. prob.');
% % % set(gca,'FontSize',40,'linewidth',2);
% % % figurenamehmm=[figurefolder,'\FalseNearestNeighborDToS_k',num2str(HMM.kNNnum),'Condi_',num2str(Nq),'States_',num2str(stateScan(1)),'bin_',num2str(binsizeScan(1)),'.jpg'];
% % % print(gcf, '-djpeg', '-r400',figurenamehmm)


filenameFNP=[figurefolder,'\FalseNearestNeighbor_k',num2str(HMM.kNNnum),'Condi_',num2str(Nq),'States_',num2str(stateScan(1)),'bin_',num2str(binsizeScan(1)),'.mat'];
save(filenameFNP,'FalseNeighbor1','FalseNeighbor2','HMM');




Temp1=[];
 h1=figure ('position', [00, 10, 1700, 1000]);
for Id_index=1:Nq
   % if min(abs(Id_index-ConditionIndex))==0
        disp(Id_index);
        KLRatio(:,Id_index)=KLDivergence(:,Id_index)./KLEntropy(:,Id_index);
pp=plot(XaxisState, KLRatio(:,Id_index),'-*','linewidth',3,'Color',cc(Id_index,:));hold on;
        pp.Color(4) = transparence;
    Temp1=[Temp1,KLRatio(:,Id_index)];
  %  end
end
b=errorbar(XaxisState,mean(Temp1,2),std(Temp1',1),'o','markersize',20,'linewidth',3,'Color','b','MarkerFaceColor','b'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
b.CapSize = 12;

h=legend(OtherPara.condition2,'Location','bestoutside');
%title(h,'States number');
set(h,'FontSize',26,'fontWeight','normal');
%set(gca,'XLim',[0 max(XaxisState)],'Xtick',2.^[1:6],'Xticklabel',2.^[1:6]) ; 
set(gca,'XLim',[0 max(XaxisState)],'Xtick',[0:round(max(XaxisState)/4):max(XaxisState)],'Xticklabel',[0:round(max(XaxisState)/4):max(XaxisState)]) ; 
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 limsy(2)]);
%set(gca,'xscale','log');
xlabel('Number of hidden states');
ylabel('KL Divergence / entropy');
set(gca,'FontSize',40,'linewidth',2);
figurenamehmm=[figurefolder,'\MeanKLDivergenceRatio_Condi_',num2str(Nq),'States_',num2str(stateScan(1)),'bin_',num2str(binsizeScan(1)),'.jpg'];
print(gcf, '-djpeg', '-r400',figurenamehmm)
%[B,I] = sort(KLRatio(:,DeleteCondiNumStatesIndex));
% for kk=1:length(I)
% BadCondition3(kk)=data.nfkb(I(kk)).ID;
% end

close all;

return
end
