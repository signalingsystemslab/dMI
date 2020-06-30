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

%This script plots the time course of average channel capacity of all
%stimuli with using hidden Markov model. It compares the WT NFkB data and
%the IkB-mutant dataset.

%% Specify the dataset and name the folder for output


function PlotCC(filename,HMM,OtherPara)
filefolder=pwd;
figurefolder=[filefolder,'\Figures\MaximumMI'];
mkdir(figurefolder)
% % % % % 
% % % % % Dataset=1;%1 is our data set; 2 is p53 after Mdmx; 
% % % % % ScaleRelation=0;
% % % % % ScanPairState=[64];%,120];%Control to Roy's
% % % % % ScanPairBin=[32];%,30];
% % % % % % ScanPairState=[80,100,120,100,100];%,120];%WT
% % % % % % ScanPairBin=[30,30,30,40,20];%,30];
% % % % % % ScanPairState=[80,100,120];%for IkB-mutant
% % % % % % ScanPairBin=[30,30,30];
% % % % % ScaledByJackknife=1;
% % % % % IntersectionRatio=[1.0116, 1]; %WT 64-32, IkB-mutant 96-32
% % % % % 
% % % % % OtherPara.Delayembedding=0;%use Delay embedding or not: 0 no, 1 and 2 differnt use
% % % % % OtherPara.Dataset=Dataset; 
% % % % % Para.HMM=1;%use Hidden Markov
% % % % % OtherPara.HMMFitToAll=1; % Fit to all traj to get a HMM and use it to sample
% % % % % OtherPara.HMMAllCondition=0; %Fit to all condition and use it to get mutual information later
% % % % % OtherPara.PartialTraining=0;   
% % % % % OtherPara.PartialRatio=1;%1/2;% , Try use half traj to train HMM
% % % % % OtherPara.ModelBasedStatistic=1;
% % % % % OtherPara.TrajBasedStatistic=1;
% % % % % OtherPara.CalculateStatistics=1;
% % % % % OtherPara.TrajBasedNormalize=0;%normalize
% % % % % OtherPara.CalculateDistance=1;
% % % % % OtherPara.EnhancedSample=2;
% % % % % CorrelationToFeatures=0;
% % % % % kk=1;jump=0;
% % % % % for DataType=[1 2]%:2% 1 is WT and 2 is IkB-mutant, 3 is Roy's method on WT data...
% % % % % for i=1:size(ScanPairState,2)
% % % % %     if DataType==1
% % % % %         LabelName='New3_HMM_Condition_All_17';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
% % % % %     elseif DataType==2
% % % % %         LabelName='New3_HMM_IkBMu_6';%'New3_HMM_TNFBk_13';%'New3_HMM_IkBMu_6';%: Ade's new full data with IkBMu and TNFBk
% % % % %     elseif DataType==3
% % % % %         LabelName='Roy';
% % % % %         filename=[filefolder,'\',LabelName,'.mat'];
% % % % %         %load(filename);
% % % % %         jump=1;
% % % % %     end
% % % % % if i==1 
% % % % %     if DataType==1
% % % % %     figurefolder=[figurefolder,'\Controls',];
% % % % % mkdir(figurefolder)
% % % % %     end
% % % % % end
% % % % % %'New3_HMM_Condition_5_4_3';%:10TNF+100P3CSK4+control;
% % % % % %'New3_HMM_Condition_5_6_5';%:10TNF+10cpG+control; 'New3_HMM_Condition_16';%: Ade's new WT data
% % % % % %LabelName='New2_HMM_Condition_11';%Ade's old WT data
% % % % % %'New2_HMM_Condition_11';%'p53_Condition_2';%'New_HMM_Condition_11';%'p53_Condition_1';%'HMM_Condition_11_P1';
% % % % % %'ERK_Condition_1_OnePhase': condition 1 is fig6D, 2 is 6E lower, 3 is 6E upper; 
% % % % % %'ERK_Condition_1_TwoPhase_1': One phase is not separate time window, Two/Three is to separate before and after drug
% % % % % 
% % % % % if jump~=1
% % % % % if Para.HMM==1
% % % % %     OtherPara.state=ScanPairState(i);%10-30,30-30,|30-50,60-50, |40-100, 100-100
% % % % %     OtherPara.binsize=ScanPairBin(i);%0;
% % % % %     OtherPara.sampleSize=10000;% number of sampling to get best fit traj for each HMM
% % % % %     OtherPara.estTR_Num=50;%Number of HMM model to save for each condition
% % % % %     LabelName=[LabelName,'_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize)];
% % % % %     if OtherPara.PartialTraining==1
% % % % %         LabelName=[LabelName,'Train_',num2str(1/OtherPara.PartialRatio)];
% % % % %     end
% % % % %     
% % % % %     if OtherPara.HMMFitToAll==1
% % % % %     LabelName=[LabelName,'_FitToAll'];
% % % % %     end
% % % % % end
% % % % % 
% % % % % foldername=[filefolder,'\',LabelName];
% % % % % filename=[filefolder,'\',LabelName,'.mat'];
% % % % % OtherPara.figurenamehmmSample=[foldername,'\HMMTrajs_States_',num2str(OtherPara.state),'bin_',num2str(OtherPara.binsize)];
% % % % % end 
% % % % % %%
% % % % % if DataType==1
% % % % %     %filename=['D:\information theory\New3_HMM_Condition_All_27_States_96bin_32_FitToAll.mat'];
% % % % %     filename=['D:\information theory\New3_HMM_Condition_All_17_States_',num2str(ScanPairState(i)),'bin_32_FitToAll.mat'];
% % % % % elseif DataType==2
% % % % % %         filename=['D:\information theory\New3_HMM_IkBMu_6_States_100bin_30_FitToAll.mat'];
% % % % % filename=['D:\information theory\New3_HMM_IkBMuWiControl_3_States_64bin_32_FitToAll.mat'];
% % % % % end
% % % % %     disp(filename);
% % % % % if exist(filename) 
% % % % % % % %     disp(i);  
% % % % % % % %     load(filename);
% % % % % % % %     if jump~=1
% % % % % % % %     if kk==1
% % % % % % % %         c = categorical(OtherPara.condition2);
% % % % % % % %         Nq=size(HMM.Entropy,2);
% % % % % % % %     end
% % % % % % % %         DataSummary.Entropy{kk}=HMM.Entropy;
% % % % % % % %         DataSummary.EntropyProduction{kk}=HMM.EntropyProduction;
% % % % % % % %         TimeWindow=20;
% % % % % %         for Id_index2=1:Nq
% % % % % %         HMM.MemoryDecayCondition(Id_index2)=mean(mean(HMM.TrajMemory{Id_index2}(:,1:TimeWindow),1));
% % % % % %         end
% % % % % %        DataSummary.MemoryDecayCondition{kk}=HMM.MemoryDecayCondition;
% % % % %         DataSummary.MI_Full{kk}=HMM.MI_Full(1:143);
% % % % % % % %       if ScaledByJackknife==1
% % % % % % % %      DataSummary.MI_Full{kk}=DataSummary.MI_Full{kk}*IntersectionRatio(kk);
% % % % % % % %      end
% % % % %    
% % % % %    
% % % % % AverageNum=3;   
% % % % % if DataType==1
% % % % % Fig=figure ('position', [00, 10, 800, 600]);
% % % % % 
% % % % % %errorbar([1:size(DataSummary.MI_Full{kk},2)]*5/60,movmean(DataSummary.MI_Full{kk},AverageNum),movstd(DataSummary.MI_Full{kk},AverageNum),...
% % % % %     %'.','markersize',20,'LineWidth',1,'color','r'); hold on;
% % % % % plot([1:size(DataSummary.MI_Full{kk},2)]*5/60,DataSummary.MI_Full{kk},...
% % % % %     '.','markersize',20,'color','r'); hold on;
% % % % % 
% % % % % x=[1:size(DataSummary.MI_Full{kk},2)]*5/60;
% % % % % curve1=DataSummary.MI_Full{kk};
% % % % % 
% % % % % end
% % % % % 
% % % % % if DataType==2
% % % % % %     errorbar([1:size(DataSummary.MI_Full{kk},2)]*5/60,movmean(DataSummary.MI_Full{kk},AverageNum),movstd(DataSummary.MI_Full{kk},AverageNum),...
% % % % % %     '.','markersize',20,'LineWidth',1,'color',[0.4940, 0.1840, 0.5560]	); hold on;
% % % % %  plot([1:size(DataSummary.MI_Full{kk},2)]*5/60,DataSummary.MI_Full{kk},...
% % % % %     '*','markersize',7,'MarkerEdgeColor','r','MarkerFaceColor','r'	); hold on;%[0.4940, 0.1840, 0.5560]
% % % % % 
% % % % % curve2=DataSummary.MI_Full{kk};
% % % % % 
% % % % % end
% % % % % 
% % % % % legendlabel{i}=[num2str(ScanPairState(i)),'-',num2str(ScanPairBin(i))];
% % % % %     else
% % % % % %         plot([1:2:35]*5/60,CC([1:2:35]),'.','markersize',20,'LineWidth',1); hold on;
% % % % % %     end
% % % % %     
% % % % %     
% % % % %  kk=kk+1; 
% % % % % end
% % % % % 
% % % % % 
% % % % % end
% % % % % end

% %Uncomment if wanted to show legend of WT and mutant...
% x2 = [x, fliplr(x)];
% inBetween = [curve1, fliplr(curve2)];
% h=fill(x2, inBetween, 'r');%[0.5 0.5 0.5]);
% set(h,'facealpha',.2)
% h=legend('WT','I\kappaB-mutant','Location','best');
% title(h,'Genotype');


% % % % % % xx=0:0.01:0.5;
% % % % % % x3 = [xx, fliplr(xx)];
% % % % % % curveY1=0*ones(1,length(xx));
% % % % % % curveY2=3*ones(1,length(xx));
% % % % % % inBetween = [curveY1, fliplr(curveY2)];
% % % % % % h1=fill(x3, inBetween, [1.000000 0.500000 0.000000]);%[0.5 0.5 0.5]);
% % % % % % set(h1,'facealpha',.2)
% % % % % % xx=0.5:0.01:2;
% % % % % % x3 = [xx, fliplr(xx)];
% % % % % % curveY1=0*ones(1,length(xx));
% % % % % % curveY2=3*ones(1,length(xx));
% % % % % % inBetween = [curveY1, fliplr(curveY2)];
% % % % % % h2=fill(x3, inBetween, 'y');%[0.5 0.5 0.5]);
% % % % % % set(h2,'facealpha',.2)
% % % % % % xx=2:0.01:11.95;
% % % % % % x3 = [xx, fliplr(xx)];
% % % % % % curveY1=0*ones(1,length(xx));
% % % % % % curveY2=3*ones(1,length(xx));
% % % % % % inBetween = [curveY1, fliplr(curveY2)];
% % % % % % h3=fill(x3, inBetween, 'g');%[0.5 0.5 0.5]);
% % % % % % set(h3,'facealpha',.2)
% % % % % % h=legend([h1,h2,h3],'Peak amplitude','Oscillation','Integral and duration','Location','northeast');
% % % % % % %title(h,'Genotype');

% % % % xx=[1:5];
% % % % x3 = [x(xx), fliplr(x(xx))];
% % % % inBetween = [curve1(xx), fliplr(curve2(xx))];
% % % % h1=fill(x3, inBetween, 'y');%[0.5 0.5 0.5]);
% % % % set(h1,'facealpha',.2)
% % % % xx=[5:24];
% % % % x3 = [x(xx), fliplr(x(xx))];
% % % % inBetween = [curve1(xx), fliplr(curve2(xx))];
% % % % h2=fill(x3, inBetween, 'g');%[0.5 0.5 0.5]);
% % % % set(h2,'facealpha',.2)
% % % % xx=[24:143];
% % % % x3 = [x(xx), fliplr(x(xx))];
% % % % inBetween = [curve1(xx), fliplr(curve2(xx))];
% % % % h3=fill(x3, inBetween, 'm');%[0.5 0.5 0.5]);
% % % % set(h3,'facealpha',.2)
% % % % h=legend([h1,h2,h3],'Peak amplitude','Oscillation','Integral, duration','Location','b');
% % % % title(h,'Source of information');
Selected=0;
if OtherPara.CCSelect==1
CCTemp=HMM.MI_Full;
NumOfSelected=length(OtherPara.CCSelectCondi);
CCSelectCondi=OtherPara.CCSelectCondi;
figurefolder=[figurefolder,'\SelectedCondi_',num2str(NumOfSelected)];
mkdir(figurefolder);
Selected=1;
end


load(filename)


if Selected==1
HMM.MI_Full=CCTemp;
OtherPara.CCSelectCondi=CCSelectCondi;
end


Fig=figure ('position', [00, 10, 800, 600]);

%errorbar([1:size(DataSummary.MI_Full{kk},2)]*5/60,movmean(DataSummary.MI_Full{kk},AverageNum),movstd(DataSummary.MI_Full{kk},AverageNum),...
    %'.','markersize',20,'LineWidth',1,'color','r'); hold on;
Time=[1:size(HMM.MI_Full,2)];%*5/60;
DataSummary=HMM.MI_Full;
plot(Time,max(HMM.MI_Full,0),...
    '.','markersize',20,'color','r'); hold on;

%xlim([0 max(Time)]);
xlabel('Time (h)');
if OtherPara.Dataset==2 % 2 is p53 after Mdmx;
         xticks([1:20:122]);xticklabels({'0','10','20','30','40','50','60'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    elseif OtherPara.Dataset==3 % 3 is Erk
         xticks([1:40:130]);xticklabels({'0','5','10','15'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    else  %elseif OtherPara.Dataset==1
        xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
       % xlim([1 max(Time)]);
    
        xlim([1 max(OtherPara.TimePointsToUse)])
        TimeUnit=round(OtherPara.TotalTimeLength/6);
        xTime=0:TimeUnit:5*TimeUnit;
        xTimePoints=round(xTime*OtherPara.TimePointsToUse/OtherPara.TotalTimeLength);    
        xTimePoints(1)=1;
        xticks(xTimePoints);xticklabels(cellstr(string(xTime)));%xticklabels({'0','2','4','6','8','10','12'});
    end
ylim([0 max(2,ceil(max(DataSummary)))]);%grid on;
%ylabel('Channel capacity (bits)');
ylabel('Maximum mutual information (bits)');
% set(h,'FontSize',18,'fontWeight','normal');
set(gca,'FontSize',24,'linewidth',2);
C = strsplit(filename,'\');
LabelName11=C{end};
C2=strsplit(LabelName11,'.');
LabelName12=C2{1};
figurenamehmm=[figurefolder,'\MaxMI_',LabelName12,'_Condi_',num2str(min(length(OtherPara.condition2),length(OtherPara.CCSelectCondi))),'.jpg'];
print(gcf, '-djpeg', '-r300',figurenamehmm)%saveas(gcf,figurenamehmm);
%close all;


 filenamehmm=[figurefolder,'\MaxMI_',LabelName11];
 save(filenamehmm,'DataSummary','OtherPara','Time');
 filenamehmm2=[figurefolder,'\MaxMI_',LabelName12,'.xlsx'];
 
 DataTable2=table(Time',DataSummary');
 writetable(DataTable2,filenamehmm2);

 
close all;
end
