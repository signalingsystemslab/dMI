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

%This script is the main script to train a time-inhomogeneous Markov model, use it to
%sample, evaluate training goodness, and calculate channel-capacity.

function [filename, HMM]=Markov(Dataset,OtherPara,X,Mutant,ID)%1 is our data set and TNFR_BK; 2 is p53 after Mdmx; 3 is Erk dataset;

OtherPara.HMM=0;
%% Load data and experiment information
% clear;
 filefolder=pwd;
% global ID
% fh = load('nfkb_dynamics_ade29-Jul-2019.mat');
% data = fh.dataTbl;
% expts = readtable('stimulus_info.xlsx'); 
% disp(expts)


%% Parameter input

OtherPara.PartialTrainTest=1;
OtherPara.PartialNum=50;

if Dataset==0
    LabelName=[OtherPara.LabelName];%'User';
    
    
elseif Dataset==1
    LabelName='NFkB';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
    if Mutant==1
     LabelName='NFkB_IkBMu';%'New3_HMM_TNFBk_13';%'New3_HMM_IkBMu_6';%: Ade's new full data with IkBMu and TNFBk
    end
elseif Dataset==2
    LabelName='p53';%'p53_New_dataset_2';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
elseif Dataset==3
    LabelName='Erk';%'p53_New_dataset_2';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';%'New3_HMM_Condition_AllNoRepli_15';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_5_4_3';%'New3_HMM_Condition_All_27';%'New3_HMM_Condition_All_33';
end
LabelName=[LabelName,'_Markov'];

if OtherPara.Scancase==4
    LabelName2=[LabelName,'_Partition'];
    mkdir(LabelName2);
else
    LabelName2=LabelName;
end
 
LabelName=[LabelName2,'_States_',num2str(OtherPara.binsize)];
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

filename=[filefolder,'\',LabelName,'.mat'];
foldername=[filefolder,'\',LabelName];
if ~exist(foldername)
    mkdir(foldername);
end


%% For Likelihood of test dataset
for i=1:length(X) %%for likelihood test
if OtherPara.PartialTrainTest==1        
        OtherPara.Ytest{i}=X{i}(:,end-OtherPara.PartialNum+1:end);
        YPartial{i}=X{i}(:,1:end-OtherPara.PartialNum);
end
end


%% Run model training for time-inhomogeneous Markov model. Analyze difference between sampled trajectories and data, and calculate trajectory probability, trajectory entropy and channel capacity
%OtherPara.PartialTrainTest=1;%%for likelihood test

if  exist(filename)
   disp('Have conducted the model training and mutual information calculation.');
   disp(['File name: ',filename]);
     
    load(filename);
    
    
else
    
if OtherPara.PartialTrainTest==1
[HMM,fitI,info1,I,Q,fX_all,XAverage,TrajEntropy,estTR,estE,OtherPara,metrics,aux]= JackknifeCCMarkov(YPartial,5,12,OtherPara); 
for i=1:size(HMM.TrajProb,1)
Temp2(i,i)=HMM.TrajProb(i,i);
end
end

OtherPara.PartialTrainTest=0;
[HMM,fitI,info1,I,Q,fX_all,XAverage,TrajEntropy,estTR,estE,OtherPara,metrics,aux]= JackknifeCCMarkov(X,5,12,OtherPara); 

HMM.TrajProb2=Temp2;%%for likelihood test
 d=size(X{1},1); 
%filename=[filefolder,'\Condition_',num2str(size(ID,2)),'.mat'];
% if exist(filename) 
%     load(filename);
% else
%     save(filename,'ID','XAverage','TrajEntropy','fX_all','estTR','estE','d','OtherPara','metrics','aux'); 
% end




%% Plot transition matrices if needed

% 
% for j=1:size(fX_all,1)
% vidfile = VideoWriter([foldername,'\movie',num2str(ID(j)),'.mp4'],'MPEG-4');
% vidfile.FrameRate = 3;
% open(vidfile);
% for k=1:size(fX_all,2)
%     figure
%     %figure ('position', [00, 10, 800, 600]);
%     %[N,c]=hist3(fX_all{1,k}','Edges',{edges edges},'CDataMode','auto','FaceColor','interp');
%     surfc(fX_all{j,k});
%     h = colorbar;
%     ylabel(h, 'Transition probability')
%     caxis([0 0.5]);
%     colormap jet;
%     view(2)
%     xlabel('NF\kappaB post-value state')
%     title(['Time: ',num2str(round(k/143*12,2)),' h'])
%     ylabel('NF\kappaB pre-value state')
%     set(gca,'FontSize',20);
% F(k) = getframe(gcf); 
% writeVideo(vidfile,F(k));
% close all;
% end
% close(vidfile);
% end
% 
% 
% 
% for j=6:6%1:1%size(fX_all,1)
% 
% for k=1:size(fX_all,2)
%     if  k<20
%      myFig=figure('position', [00, 10, 600, 600]);
%         %figure ('position', [00, 10, 800, 600]);
%     %[N,c]=hist3(fX_all{1,k}','Edges',{edges edges},'CDataMode','auto','FaceColor','interp');
%     surfc(fX_all{j,k});
%     h = colorbar;
%     ylabel(h, 'Probability')
%     caxis([0 1]);
%     set(gca,'xtick',[])
%     set(gca,'xticklabel',[])
%     set(gca,'ytick',[])
%     set(gca,'yticklabel',[])
%     colormap jet;
%     view(2)
%     xlabel('Post-state')
%     title(['Time: ',num2str(round(k/143*12,2)),' h'])
%     ylabel('Pre-state')
%     set(gca,'FontSize',40);
%     set(findall(myFig, 'Type', 'Text'),'FontWeight', 'Normal')
%     figurename2=[foldername,'\TransitionMatrix_',num2str(ID(j)),'_',num2str(k),'.jpg'];
%  saveas(gcf,figurename2); 
% close all;
%     end
% end
% end
% 
% ddddddd
% close all;

%%
% if Para.HMM==1
%     
% Nq=size(estE,2);  
% sampleSize=1000;
% TrajSample=cell(1,Nq);
% TrajSample_Smooth=cell(1,Nq);
% opt= maketicks(size(TrajSample{1,1},2),[0.01, 4],0); opt.Name = '';
% HM= figure ('position', [00, 10, 600, 800]);clf;
% for Id_index=1:Nq
%     for jjj=1:sampleSize
% [seq,states] = hmmgenerate(d,estTR{1,Id_index},estE{1,Id_index});
% TrajSample{1,Id_index}=[TrajSample{1,Id_index}',seq']';
%     end
%     
% 
% for jjj=1:sampleSize
% TrajSample_Smooth{1,Id_index}(jjj,:) = max(smooth(1:d,TrajSample{1,Id_index}(jjj,:),0.1,'loess'),0);
% end
% 
% %[tblB,index]=sortrows(TrajSample_Smooth{1,Id_index},[1:1:120]);
% %s=surf(TrajSample_Smooth{1,Id_index});
% h=colormapStack(TrajSample_Smooth{1,Id_index},[],opt, HM);
% hh = colorbar;
% s.EdgeColor = 'none';view(2);
% xlabel('Time (h)');ylabel('Cell label');
% % zlim([0 4])
% xticks([0:24:150]);xticklabels({'0','2','4','6','8','10','12'});
% set(gca,'FontSize',20);
%      figurename2=[foldername,'\SampleTrajHMM_',num2str(ID(Id_index)),'.jpg'];
%  saveas(gcf,figurename2); 
% end
% end
%% make a movie
% OtherPara.binsize=length(fX_all{1,1});
% OtherPara.MaxValue=10; 
% OtherPara.MinValue=0; 
% edges=linspace(OtherPara.MinValue,OtherPara.MaxValue,OtherPara.binsize);
% for j=1:size(fX_all,1)
% vidfile = VideoWriter([foldername,'\movie',num2str(ID(j)),'.mp4'],'MPEG-4');
% vidfile.FrameRate = 3;
% open(vidfile);
% for k=1:size(fX_all,2)
%     figure
%     %[N,c]=hist3(fX_all{1,k}','Edges',{edges edges},'CDataMode','auto','FaceColor','interp');
%     surfc(fX_all{j,k});
% %      x = 0:1:9;
% %     z = 0:1:9;
% %     [X,Z] = ndgrid(x,z);
% %     mesh(X,Z,fX_all{1,k})
%     %mesh(0:1:9, 0:1:9, fX_all{1,k})
%     colorbar
%     view(2)
%     xlabel('NFkB post-value state')
%     title([num2str(round(k/143*12,2)),' h'])
%     ylabel('NFkB pre-value state')
%     set(gca,'FontSize',20);
%     
% % 
% F(k) = getframe(gcf); 
% writeVideo(vidfile,F(k));
% close all;
% end
% close(vidfile);
% end



%% Sample trajecotrise from time-inhomogeneous Markov model and plot the heatmaps
   
       
    Tosmooth=1;


   TrajSample=cell(1,size(fX_all,1));
  disp('do simulation to sample from Markov model');
  
for j=1:size(fX_all,1)% jth condition
     sampleSize(j)=size(X{j},2);
for ii=1:sampleSize(j)     
chain = zeros(1,size(fX_all,2)+1);
Trajchain = zeros(1,size(fX_all,2)+1);

%fX_all{j,1}
Joint=fX_all{j,1};
Marginal=sum(Joint,2);   

cumulative_distribution = cumsum(Marginal);
r = rand();
chain(1) = find(cumulative_distribution>r,1);
Trajchain(1)=chain(1);%XAverage{j,1}(1,chain(1));

for k=1:size(fX_all,2)
    %Temp=transition_probabilities;
    Temp=cumulative_distribution;
    Joint=fX_all{j,k};
    Marginal=sum(Joint,2);
    transition_probabilities=Joint;
    transition_probabilities(Marginal>0,:)=Joint(Marginal>0,:)./Marginal(Marginal>0);
    
    this_step_distribution =transition_probabilities(chain(k),:);
    cumulative_distribution = cumsum(this_step_distribution);
    r = rand();
    if isempty(find(cumulative_distribution>r,1))
        chain(k+1)=chain(k)-1;
       % Trajchain(k+1)=XAverage{j,1}(1,chain(k));
    else
        chain(k+1) = find(cumulative_distribution>r,1);
    end
    %Trajchain(k+1)=XAverage{j,1}(1,chain(k));
    Trajchain(k+1)=chain(k);
end

TrajSample{1,j}(ii,:)=Trajchain;
end

end


%disp('now herelll')
TrajSample_Smooth=cell(1,size(fX_all,1));
    for ll=1:size(TrajSample,2)
    sampleSize(ll)=size(X{ll},2);
YUpLimit=5;  
if OtherPara.Dataset~=1
    YUpLimit=OtherPara.YUpLimit(ll);
end
YUpLimit=OtherPara.YUpLimit(ll);
%opt= maketicks(size(TrajSample{1,1},2),[0.01, 4],0); opt.Name = '';
opt= maketicks(1:size(TrajSample{1,ll},2),[0, YUpLimit],0); opt.Name = 'ddd';
HM= figure ('position', [00, 10, 600, 800]);%('name','650');  
%Trajectories is M x N matrix here
        Id_index=ll;
        TrajSample{1,Id_index}=TrajSample{1,Id_index}/OtherPara.binsize*OtherPara.MaxValue;
  if Tosmooth==0
    s=surf(TrajSample{1,Id_index});
    s.EdgeColor = 'none';
    view(2);set(gca,'FontSize',20);
  else
    lengthTime=size(TrajSample{1,Id_index},2);
    disp(lengthTime);
    
    for jjj=1:sampleSize(ll)
    TrajSample_Smooth{1,Id_index}(jjj,:) = max(smooth(1:lengthTime,TrajSample{1,Id_index}(jjj,:),0.07,'loess'),0);
    %aa=smooth([1:1:lengthTime]',[TrajSample{1,Id_index}(jjj,:)]',0.1,'loess');
    %TrajSample_Smooth=max(smooth(1:lengthTime,TrajSample{1,Id_index}(jjj,:),0.1,'rloess'),0);
    end
   
    %s=surf(TrajSample_Smooth{1,Id_index});
    
    %disp(size(TrajSample_Smooth{1,Id_index}));ddd
    [metrics,aux] = nfkbmetrics(TrajSample_Smooth{1,Id_index});
    [aaa,index]=sortrows(metrics.peakfreq);
    %disp(size(index));
    tblB=TrajSample_Smooth{1,Id_index}(flipud(index),:);
    %[tblB,index]=sortrows(TrajSample_Smooth{1,Id_index},[1:1:100]);
    
    h=colormapStack(tblB,[],opt, HM);
    
     %s=surf(tblB);
    hh = colorbar;
    
     s.EdgeColor = 'none';view(2);
    xlabel('Time (h)');ylabel('Single cells');
    if OtherPara.Dataset==2 % 2 is p53 after Mdmx;
         xticks([1:20:122]);xticklabels({'0','10','20','30','40','50','60'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    elseif OtherPara.Dataset==3 % 3 is Erk
         xticks([1:40:130]);xticklabels({'0','5','10','15'});
         %xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
      elseif OtherPara.Dataset==1
        xticks([1:24:130]);xticklabels({'0','2','4','6','8','10','12'});
    else
        xlim([1 max(OtherPara.TimePointsToUse)])
        TimeUnit=round(OtherPara.TotalTimeLength/6);
        xTime=0:TimeUnit:5*TimeUnit;
        xTimePoints=round(xTime*OtherPara.TimePointsToUse/OtherPara.TotalTimeLength);    
        xTimePoints(1)=1;
        xticks(xTimePoints);xticklabels(cellstr(string(xTime)));%xticklabels({'0','2','4','6','8','10','12'});
    end
    set(gca,'FontSize',36);
  end
   figurename2=[foldername,'\SampleTraj_',num2str(ID(Id_index)),'States_',num2str(OtherPara.state),'.jpg'];
 saveas(gcf,figurename2); 
    end
    
 %% Quantification on the training goodness of hidden Markov model, including False nearest neighbor probability, Relative KL-divergence, Rescaled likelihood.
 
  seqs=X;
OtherPara.EnhancedSample=1;
OtherPara.binsize=OtherPara.state;
%OtherPara.state=OtherPara.state;
OtherPara.DistanceAsDistribution=1;
Nq=length(X);%length(StatesScan);
HMM.Temp=1;
[HMM]=MarkovDistance(estTR,estE,OtherPara,seqs,X,Nq,foldername,HMM,ID,TrajSample_Smooth);

close all;
  
%save(filename,'ID','XAverage','TrajEntropy','fX_all','TrajSample','estTR','estE','d','OtherPara','metrics','aux','TrajSample','TrajSample_Smooth');


%% Save the simulation result
HMMPlotting(HMM,seqs,Nq,ID,OtherPara,foldername,TrajSample_Smooth);
save(filename,'Nq','seqs','ID','XAverage','TrajEntropy','fX_all','TrajSample','estTR','estE','d','OtherPara','metrics','aux','TrajSample','HMM','TrajSample_Smooth','foldername');

end


close all;
end