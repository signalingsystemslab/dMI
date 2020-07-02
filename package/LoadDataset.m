function [X, condition, OtherPara,ID]=LoadDataset(OtherPara,MainPara)

if OtherPara.Dataset==0
% Load data and experiment information from user;
disp('Load test dataset or import data from user');
%load('data1.mat');
load([pwd,'\',OtherPara.UserDataName]);
condition=OtherPara.UserCondition;
if exist('OtherPara.UserID')
    ID=OtherPara.UserID;%[1,2];
else
    ID=1:1:length(condition);
end

MaxUniform=0;
    for jj=1:length(X)
        MaxUniform=max(max(max(X{jj})),MaxUniform);
    end
    for jj=1:length(X)
        X{jj}=X{jj}/MaxUniform*10;
        OtherPara.YUpLimit(jj)=round(max(max(X{jj})),1);
    end
    
elseif OtherPara.Dataset==1
% Load data and experiment information for dataset 1

filename1='NFkB_WT.mat';
if MainPara.Mutant==1
filename1='NFkB_KO.mat';
end

if exist(filename1)
    disp('Import data');
    load(filename1);
else
    disp('Please put the data from the supplementary of the paper into the code folder');
end
%fh = load('nfkb_dynamics_ade29-Jul-2019.mat');
%data = fh.dataTbl;
%expts = readtable('stimulus_info.xlsx'); 
%disp(expts)
%[IDload, txt2]= xlsread('stimulus_info.xlsx',1,'A1:A54');

% Choose the experimental ID of the selected conditions:
% For our dataset, the index number for the experiment under the various conditions
% 783 control,650 10TNF, 779 100LMW-PIC, 610 100cpG, 720 100P3CSk, 756 100LPS, 566 50PIC, 780 33PIC,
% 778 33LMW-PIC,  546 10P3CSk, 548 1P3CSk
%Remove conditions with <300 trajs: 147 445 446 447 450 757. After removing:
%stim_ID =[ 548 779 650 754 756 664 777 581 720 157  610 546 778 780 755 663 783]; %selected WT: very few replicate
%IkBMu
% if MainPara.Mutant==1
% %stim_ID =[ 751,759,752,760,753,761];% IkBMu with replicate
% stim_ID =[ 759,760,761,783];% IkBMu without replicate
% end

%disp(stim_ID);

%index=ismember(IDload, stim_ID);
%ID_ordered=IDload(index)';
%Stimulus=expts{:,3}(index);
%Concentration=expts{:,4}(index);
%Unit=expts{:,5}(index);

%ctrl_ID = [];
%stim_ID = IDload';%546:548; 
%ID_ordered = [ctrl_ID, ID_ordered];
[ID I]=sort(ID_ordered);
% condition=cell(1,size(Unit,1));
% for i=1:size(Unit,1)
%     str=[num2str(Concentration(i)),Unit(i),Stimulus(i)];%,num2str(ID_ordered(i))];
%     condition{i}=join(str," ");
% end
% Stimulus=Stimulus(I);
% condition=condition(I);

% Get the row names of trajectories of the selected conditions
%idx =  ismember(data.ID, ID);

% Specify the number of time points, where 'txt' are the time points 
%[ss, txt]= xlsread('TimeSeries.xlsx',1,'A1:A143');%A143

txt2=txt;
featureNames = txt';
display(ID);

% Construct a master dataset from the measured time series.
%X= construct_mi_mat(data(idx,:), featureNames); 


elseif OtherPara.Dataset==2
%% load data of p53 signaling response
     
    ID=[1 2];
    X=cell(1,length(ID));condition=cell(1,length(ID));
   % condition={'1','2'};
    condition={'p53 with oscillation','p53 without oscillation'};
    load('Mdmx_siRNA_traces.mat');
    load('UV16cnt_traces.mat');
    %MaxUniform=max(max(max(MdmXsiRNA)),max(max(UV16cnt_traces)));
    X{1}=MdmXsiRNA;%MdmXsiRNA;
    X{2}=UV16cnt_traces;%MdmXsiRNA;
    MaxUniform=0;
    for jj=1:length(X)
        MaxUniform=max(max(max(X{jj})),MaxUniform);
    end
    for jj=1:length(X)
        X{jj}=X{jj}/MaxUniform*10;
        OtherPara.YUpLimit(jj)=round(max(max(X{jj})),1);
    end
    X{2}=X{2}(1:size(X{1},1),:);
    %X2=X;
elseif OtherPara.Dataset==3
%% load data of ERK, p38, JNK signaling response
    
    ID=[1 2 3]; %ERK, p38, JNK
    X=cell(1,length(ID));condition=cell(1,length(ID));
    condition={'ERK', 'p38', 'JNK'};
    load('DataERK.mat');
    TimeWindow1=50;%Visually
    ErkCondition=1;%1,2,3
    X{1}=Data.ERK{1,ErkCondition}';
    X{2}=Data.p38{1,ErkCondition}';
    X{3}=Data.JNK{1,ErkCondition}';
    
     MaxUniform=0;
    for jj=1:length(X)
        MaxUniform=max(max(max(X{jj})),MaxUniform);
    end
    for jj=1:length(X)
        X{jj}=X{jj}/MaxUniform*10;
        OtherPara.YUpLimit(jj)=round(max(max(X{jj})),1);
    end
    
    %X2=X; 
elseif OtherPara.Dataset==4
    
%% load synthetic data
%1 is frequency modulation, and 2 is amplitude modulation
TrajFile=[foldername,'\Trajectory.mat'];
ID=[1 2];

if ~exist(TrajFile)
tt = linspace(0,12,TimePointsToUse);
for ii=1:500
    if TrajModes==1
        X{1}(:,ii)=5*(sin((tt+rand*4)*pi)+1);
        X{2}(:,ii)=5*(sin((tt+rand*4)*pi/2)+1);
        Mode='Frequencey modulation ';
    elseif TrajModes==2
        X{1}(:,ii)=5*rand*(sin(tt*pi)+1);
        X{2}(:,ii)=5*rand*(sin(tt*pi/2)+1);
         Mode='Amplitude modulation ';
    end
end
%X2=X;
for jj=1:length(X)
    condition{jj}=[Mode,num2str(jj)];
end
disp('save trajectory file');
save(TrajFile,'X','condition');

else
  disp('load trajectory file');
  load(TrajFile);  
end


end



%% Pre-processing on all the trajectories

Nq=length(X);% Y{1,j} is  j-th condition: with row time series, column trajs
disp('number of conditions');disp(Nq);
Y=cell(Nq,1);%Y2=cell(Nq,1);
OtherPara.TimeWiseNum=size(X{1},1);% Train point with early time points, time-wise
if OtherPara.TimePointsToUse<140
    OtherPara.TimeWiseNum=OtherPara.TimePointsToUse;
end

OtherPara.MaxValue=-Inf;OtherPara.MinValue=Inf;
OtherPara.ID=ID;
for i=1:size(condition,2)
    OtherPara.condition2(i)=cellstr(condition{i});
end

for i=1:Nq
   
     A=X{i};
    
     %B=X{i};
    if size(A,1)==1 % Unidimensional case: remove NaN individuals
        A=A(:,~isnan(A));
        A=A(:,~isinf(A));
       % B=B(:,~isnan(A));
        %B=B(:,~isinf(A));
    else % Multidimensional case: remove individuals with NaN in any dimension
        A=A(:,sum(isnan(A))==0);
        A=A(:,sum(isinf(A))==0); 
       % B=B(:,sum(isnan(A))==0);
       % B=B(:,sum(isinf(A))==0); 
    end
    A(A<0)=eps; %Make negative to be zero
    %B(B<0)=eps; %Make negative to be zero
    
    if OtherPara.Dataset==1 || OtherPara.Dataset==4
     OtherPara.Avector{i} = sort(A(:));
     AllValues{i}=OtherPara.Avector{i}(1:round(length(OtherPara.Avector{i})*0.99));%90% percentile of cells;
     MaxPercentile(i)=max(AllValues{i});%OtherPara.YUpLimit(i)=ceil(max(max(A)));
     
     A(A>10)=10; % To make MaxValue smaller for better training.
     OtherPara.YUpLimit(i)=ceil(MaxPercentile(i));
    
    end
    
    Y{i}=A;
    OtherPara.MaxValue=max(OtherPara.MaxValue,max(max(Y{i})));
    
    OtherPara.MinValue=min(OtherPara.MinValue,max(min(min(Y{i})),0));   
   % B(isnan(B))=0;
    %Y2{i}=B;
    

    if size(Y{i},2)<OtherPara.PartialNum
       OtherPara.PartialNum=size(Y{i},2);     
    end
end

X=Y;
%Basic statistics: plot histgram of all the data


end
