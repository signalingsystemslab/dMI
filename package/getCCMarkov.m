%Sub-script to train a time-inhomogeneous Markov model by calcuating the transition matrix, and calculate channel-capacity.

function [HMM,info1,Q,fX_all,XAverage,TrajEntropy,metrics,aux,PointsDelay,OtherPara] = getCCMarkov (Y,K,OtherPara,varargin)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% getCC calculates channel capacity for an input cell matrix, Y.
%
% Y         Cell array (size = 1 x n) of sets of individual responses to n different inputs.
%           Each cell is of size [r,c] -> r dimensions x c individuals (r is same across
%           all cells of Y) 
% K         Use Kth nearest neighbor (passed to approxmiate KNN algorithm)
% idisp     (true/false) show verbose output
% varargin  (optional) entropies provided
%
%
% fitI      Extrapolated mutual information for entire set
% info1     Mutual information for full set
% I         Mutual information for all subsets
% Q         "Optimal" weights (of all possible inputs) that maximizes mutual information
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Create new array, X, which holds nan-filtered data 
Nq=length(Y);% Y{1,j} is  j-th condition: with row time series, column trajs
PointsDelay=0;
X=Y;
XAverage=cell(Nq,1);
% % % % X=cell(Nq,1);
% % % % 
% % % % OtherPara.MaxValue=-Inf;OtherPara.MinValue=Inf;
% % % % for i=1:Nq
% % % %     A=Y{i};
% % % %     if size(A,1)==1 % Unidimensional case: remove NaN individuals
% % % %         A=A(:,~isnan(A));
% % % %         A=A(:,~isinf(A));
% % % %     else % Multidimensional case: remove individuals with NaN in any dimension
% % % %         A=A(:,sum(isnan(A))==0);
% % % %         A=A(:,sum(isinf(A))==0); 
% % % %     end
% % % %     A(A<0)=eps; %Make negative to be zero
% % % %     X{i}=A;
% % % %     OtherPara.MaxValue=max(OtherPara.MaxValue,max(max(X{i})));
% % % %     OtherPara.MinValue=min(OtherPara.MinValue,max(min(min(X{i})),0));
% % % % end
% % % % % X  is transponse of Y with removing a few trajs with NaN, typicall <5%
% % % % 
% % % % if OtherPara.Dataset==1
% % % % OtherPara.MaxValue=10; %change smaller: ignore too large
% % % % end
% % % % %OtherPara.MaxValue=round(18.1520826097308); %change to the biggest for all data
% % % % %OtherPara.MaxValue=8;%Mannul put, because max learned has 15, which is over the figures'.
% % % % d=size(X{1},1);
OtherPara.binsize=OtherPara.state;%32;
OtherPara.MemoryLength=20;% number of time points 
OtherPara.kMeans=10;
OtherPara.NonMarkovian=0;

%%
for  i=1:Nq
   % [metrics(i),aux(i)]=[1,1];
[metrics(i),aux(i)] = nfkbmetrics(X{i}');%,varargin);
end



%%
%Delay embedding's mutual information for optimal delay T
OtherPara.Delayembedding=0;
info1=0;Q=0;fX_all=0;%TrajEntropy=0;


% % % if OtherPara.Delayembedding==2
% % %     
% % % OtherPara.DelayLength=2; 
% % % 
% % % OtherPara.DelayDim=round(140/OtherPara.DelayLength)-1;% Two points, i.e. two dimension first...
% % % disp(OtherPara.DelayDim);
% % % PointsDelay=cell(Nq,1);
% % % %OtherPara.DelayLength*OtherPara.DelayDim;
% % % NumOfVector=OtherPara.DelayLength;
% % % for i=1:Nq
% % %     PointsDelayTemp=cell(NumOfVector,size(X{i},2));%,Nq); % row-th delay length for col-th traj
% % %   
% % % for ll=1:size(X{i},2) % Trajs
% % %     for m=1:OtherPara.DelayDim
% % %     
% % %     for k=1:d-OtherPara.DelayLength*(m-1) %Time points
% % %      
% % %     PointsDelayTemp{m,ll}(k,:)=X{i}(k:OtherPara.DelayLength:k+OtherPara.DelayLength*(m-1),ll);
% % %     end
% % %     
% % %     end
% % %     
% % %     
% % %     
% % % end
% % %     
% % %     PointsDelay{i,1}=PointsDelayTemp;
% % % %end
% % % %disp(size(PointsDelay{1,1}));
% % % %disp(size(PointsDelay{2,1}));
% % % 
% % % end
% % % 
% % % 
% % % return;
% % %     
% % % elseif OtherPara.Delayembedding==1 % find data to detect best delay length
% % % 
% % % 
% % % OtherPara.DelayLength=15;  
% % % OtherPara.DelayDim=3;% Two points, i.e. two dimension first...
% % % PointsDelay=cell(Nq,1);
% % % %OtherPara.DelayLength*OtherPara.DelayDim;
% % % NumOfVector=OtherPara.DelayLength;
% % % for i=1:Nq
% % %     PointsDelayTemp=cell(NumOfVector,size(X{i},2));%,Nq); % row-th delay length for col-th traj
% % % for j=1:OtherPara.DelayLength    
% % %     for k=1:d-j*(OtherPara.DelayDim-1) %Time points
% % %     for ll=1:size(X{i},2) % Trajs  
% % %     PointsDelayTemp{j,ll}(k,:)=X{i}(k:j:k+j*(OtherPara.DelayDim-1),ll);
% % %     end
% % %     
% % %     
% % %     
% % %     
% % %     end
% % %     
% % %     PointsDelay{i,1}=PointsDelayTemp;
% % % end
% % % %disp(size(PointsDelay{1,1}));
% % % %disp(size(PointsDelay{2,1}));
% % % 
% % % end
% % % 
% % % 
% % % 
% % % return;
% % % end



%%
%Here: repeat calling conditional_probabilities twice and revise the
%calculation on hRS and hR. Need a new getInfo with taking the two fX and
%HR. So, first a function return two fX, and a function return new hRS, hR
%and hRS-hR
% Calculate conditional probabilities and entropies - if getting infinite entropy, increase added noise
global normalization
timeunit=1;% hour, 5min revise later
normalization=0;
noiselvl = eps;
d=size(X{1},1); %dimension %vector dimension or # of time series 


%%
%Non-Markovian memeory length trajs



% % % if OtherPara.NonMarkovian==1
% % % Xsize=zeros(Nq+1,1);
% % % Xcluster=cell(Nq,1);
% % % for k=1:d-OtherPara.MemoryLength   
% % %     XX=[];
% % %     for i=1:Nq
% % %     XX=[XX,X{i}(k:k+OtherPara.MemoryLength,:)];
% % %     if k==1
% % %     Xsize(i+1)=size(X{i},2);
% % %     end
% % %     end
% % % 
% % %     idx3 = kmeans(XX',OtherPara.kMeans);%,'Distance','cityblock');
% % %     for i=1:Nq
% % %     Xcluster{i}(k,:)=idx3(sum(Xsize(1:i))+1:sum(Xsize(1:i+1)));
% % %     for m=1:OtherPara.kMeans
% % %         [minValue,closestIndex] = min(abs( Xcluster{i}(k,:)-m));
% % %         XAverage{i}(k,m) = mean(X{i}(k+OtherPara.MemoryLength,closestIndex),2);
% % %     end
% % %     %
% % %     end
% % % end
% % % X=Xcluster;
% % % end
% disp(size(X{i}))
% disp(size(Xcluster{i}))
% disp(size(XAverage{i}))
% ddd
%%

%fX_all=cell(Nq,Nq,d-1);
%fX1_all=cell(Nq,Nq,d-1);
fX_all=[];


for k=1:d-1
    XX=cell(Nq,1);
    for i=1:Nq
    XX{i}=X{i}(k:k+1,:);
    end
do_flag = 1;
OtherPara.k=k;
while do_flag
    %[fX1, fX] = conditional_probabilities_kNN(XX,K,noiselvl,varargin);
    fX = conditional_probabilities_bin(XX,OtherPara);
    
    %Need to check how sparse the transtion matrix: whether increase sample
    %size of reduce bin number
    suminf = 0;
    for i = 1:numel(fX)
        suminf = suminf+sum(isinf(fX{i}(:)));
    end
    if suminf==0 
        do_flag = 0;
    else
        noiselvl = noiselvl*10;
    end
%     suminf1 = 0;suminf = 0;
%     
%     for i = 1:numel(fX)
%         suminf1 = suminf1+sum(isinf(fX1{i}(:)));
%         suminf = suminf+sum(isinf(fX{i}(:)));
%     end
%     if suminf1==0 && suminf==0
%         do_flag = 0;
%     else
%         noiselvl = noiselvl*10;
%     end
    
%     %Show transition matirx
%     TransitionMarix(XX,fX1,fX);

end
%fX=cell2mat(fX);
%fX=fX/sum(sum([fX{:,:}]));% normalize the transition probability

fX_all=[fX_all,fX];

end
%fX_all is matrix of cell, fX_all{i,j} is the transition matrix
%for the i-th condition and j-th time point






% % If probabilities are provided, assign into fX
% if nargin>3 
%     if numel(varargin{1})>0    
%         fX=Ftrue;    
%     end
% end


%%

% TrajEntropy=0;
% return;



%%

%fX_all
edges=linspace(OtherPara.MinValue,OtherPara.MaxValue,OtherPara.binsize);

if OtherPara.PartialTrainTest==0
for kk=1:Nq
    NumOfTraj=size(X{kk},2);
    for jj=1:NumOfTraj

    [N,edgesNew,seqs{kk}(jj,:)] = histcounts(X{kk}(:,jj)',edges);
    end
end


elseif OtherPara.PartialTrainTest==1 %%for likelihood test
    for kk=1:Nq
    NumOfTraj=size(OtherPara.Ytest{kk},2);
    for jj=1:NumOfTraj

    [N,edgesNew,seqs{kk}(jj,:)] = histcounts(OtherPara.Ytest{kk}(:,jj)',edges);
    end
    end
end


%   assignin('base', 'edges1', edges);
%   assignin('base', 'seqs', seqs); ddd
%  assignin('base', 'fX_all', fX_all); 
%assignin('base', 'OtherPara2', OtherPara);ddd
[HMM]=MMStatistics(fX_all,OtherPara,X,seqs,Nq);
% assignin('base', 'HMM1', HMM); 




%%
% (STEP 2) Optimize input probabilities for maximum channel capacity
% Channel capacities are calculated as the difference between the sum of conditional entropies H(R|S),
% and the total entropy of the response, H(R).
%
%     - H(R|S) = sum(hRS);
%     - H(R) = doubly-weighted sum (see getInfo.m) 


% % Set options (toggle verbose/non-verbose display)
% if idisp==0
%     %options = optimoptions('fmincon','Algorithm','active-set','Display','off');
%     options = optimoptions('fmincon','Algorithm','interior-point');
% elseif idisp==1
%     %options = optimoptions('fmincon','Algorithm','active-set','Display','iter');
%     options = optimoptions('fmincon','Algorithm','interior-point','Display','iter');
% end
options = optimoptions('fmincon','Algorithm','interior-point');

% Define information calculation function (getInfo.m), initialize bounds for optimization.
%infoF= @(q) -getInfo(q,fX1_all,fX_all,timeunit);
infoF= @(q) -getInfoModelMarkov(q,fX_all);
Aeq=ones(1,Nq);
beq=1;
LB=zeros(Nq,1);
UB=ones(Nq,1);  

% assignin('base', 'PointsDelay', PointsDelay); 
 
if nargin>4 % q's are provided; skip optimization and output.
    if numel(varargin{2})>0
        Qstart=varargin{2};
        info1=-infoF(Qstart);  
    end
else % Optimize q for information content, given probability distributions 
    Qstart=ones(Nq,1);
    Qstart=Qstart/sum(Qstart);    
    qfit=fmincon(infoF,Qstart,[],[],Aeq,beq,LB,UB,[],options);
    Q=qfit;
    %info1=-infoF(qfit);
    %[xx ,TFrate_OP]=getInfo(qfit,fX1_all,fX_all,timeunit);
    [info1,aaa, HMM.TrajEntropyModelBased]=getInfoModelMarkov(qfit,fX_all);
    assignin('base', 'aaa', aaa); 
    TrajEntropy.optimal=aaa;
    [f,bbb ccc]=getInfoModelMarkov(Qstart,fX_all);
    TrajEntropy.uniform=bbb;
end


