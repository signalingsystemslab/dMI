% Train a hidden Markov model for each stimulus condition's path ensemble

function [SampledTraj,DistanceTotal, estTR,estE,estTR_sc,estE_sc,OtherPara,seqs] = HMMTrain (X,OtherPara)


% Create new array, X, which holds nan-filtered data 
Nq=length(X);% Y{1,j} is  j-th condition: with row time series, column trajs


% X=cell(Nq,1);
% OtherPara.MaxValue=-Inf;OtherPara.MinValue=Inf;
% for i=1:Nq
%     A=Y{i};
%     if size(A,1)==1 % Unidimensional case: remove NaN individuals
%         A=A(:,~isnan(A));
%         A=A(:,~isinf(A));
%     else % Multidimensional case: remove individuals with NaN in any dimension
%         A=A(:,sum(isnan(A))==0);
%         A=A(:,sum(isinf(A))==0); 
%     end
%     A(A<0)=eps; %Make negative to be zero
%     X{i}=A;
%     OtherPara.MaxValue=max(OtherPara.MaxValue,max(max(X{i})));
%     OtherPara.MinValue=min(OtherPara.MinValue,max(min(min(X{i})),0));
%     %disp('regularize trajs');
% end




%[SampledTraj,DistanceTotal, estTR,estE,estTR_sc,estE_sc,OtherPara,seqs] = HiddenMarkov (X,OtherPara,Nq);








%function [SampledTraj,DistanceTotal, estTR,estE,estTR_sc,estE_sc,OtherPara,seqs] = HiddenMarkov (X, OtherPara,Nq)


%tic;
% OtherPara.state=20;
% OtherPara.binsize=40;
%OtherPara.MaxValue=4;
disp(OtherPara.MaxValue);
disp(OtherPara.MinValue);


Conversion=(OtherPara.MaxValue-OtherPara.MinValue)/OtherPara.binsize;
%disp('ee');disp(Conversion);ddd
OtherPara.Conversion=Conversion;


trans=rand(OtherPara.state,OtherPara.state);
emis=rand(OtherPara.state,OtherPara.binsize);%ones(OtherPara.state,OtherPara.binsize);


trans=trans./sum(trans,2);
emis=emis./sum(emis,2);



edges=linspace(OtherPara.MinValue,OtherPara.MaxValue,OtherPara.binsize+1);


estTR=cell(1,Nq+1);
estE=cell(1,Nq+1);
SampledTraj=cell(1,Nq);
DistanceTotal=cell(1,Nq);


seqs=cell(1,Nq);
for kk=1:Nq
    %Used # of trajs
    NumOfTraj=size(X{kk},2);
    if OtherPara.PartialTraining==1   
        NumOfTraj=round(NumOfTraj*OtherPara.PartialRatio);
    end 
    disp(['Number of trajectories under the stimulus is: ', num2str(NumOfTraj)]);
     OtherPara.NumTraj(kk)= NumOfTraj;
    %NumOfTraj=9;


OtherPara.estTR_Num=10;estTR_sc=cell(1,OtherPara.estTR_Num);estE_sc=cell(1,OtherPara.estTR_Num);

maxiter=1000;





for jj=1:NumOfTraj


[N,edgesNew,seqs{kk}(jj,:)] = histcounts(X{kk}(:,jj)',edges);

%seqs{kk}(jj,:)=seqs{kk}(jj,:);

warning('off');
if OtherPara.HMMFitToAll==1
    if jj==1
        OtherPara.TrajLength(kk)=length(seqs{kk}(jj,:));        
    end
    continue;
end


[estTRTemp,estETemp] = hmmtrain(seqs{kk}(jj,:),trans,emis,'maxiterations',maxiter,'algorithm','BaumWelch');


if jj<=OtherPara.estTR_Num
        estTR_sc{1,jj}=estTRTemp;estE_sc{1,jj}=estETemp;
end
warning('off');
%OtherPara.sampleSize=10000;%TrajSample=[];
Distance=Inf;
for jjj=1:OtherPara.sampleSize
[seqSampleTemp,states] = hmmgenerate(length(seqs{kk}(jj,:)),estTRTemp,estETemp);
DistTemp=sqrt(sum((seqs{kk}(jj,:)-seqSampleTemp).^2)).*Conversion;
%disp(DistTemp);dd
if DistTemp<Distance
    seqSample=seqSampleTemp.*Conversion;
    Distance=DistTemp;
end
end
SampledTraj{1,kk}(jj,:)=seqSample;%.*Conversion;
DistanceTotal{1,kk}(jj,:)=Distance;

figure('position', [00, 10, 1200, 600])
if jj<=9
subplot(3,3,jj)
%plot(1:length(seqs{kk}(jj,:)),seqSample.*Conversion,'Color','b','linewidth',1);hold on;
%plot(1:length(seqs{kk}(jj,:)),seqs{kk}(jj,:).*Conversion,'r','linewidth',1);hold on;
plot(1:length(seqs{kk}(jj,:)),seqSample,'Color','b','linewidth',1);hold on;
plot(1:length(seqs{kk}(jj,:)),seqs{kk}(jj,:).*Conversion,'r','linewidth',1);hold on;
plot(1:length(seqs{kk}(jj,:)),X{kk}(:,jj)','k','linewidth',1);hold on;
xlabel('Time (h)');
if OtherPara.Dataset==2
    xticks([0:20:120]);xticklabels({'0','10','20','30','40','50','60'});ylabel('p53 (a.u.)');
else
xticks([0:24:150]);xticklabels({'0','2','4','6','8','10','12'});ylabel('NFkB (a.u.)');
end
set(gca,'FontSize',12);
figurename2=[OtherPara.figurenamehmmSample,'_',num2str(OtherPara.ID(kk)),'.jpg'];
 saveas(gcf,figurename2); 
 
 
end


disp(['Number of trajectories X Number of time points:',num2str(size(SampledTraj{1,kk}))]);
end

%toc;
end




if OtherPara.HMMFitToAll==1
    warning('off');
    disp('Use all trajectories of each stimulus to train hidden Markov model and use it to sample');
    %disp(Nq);
   % assignin('base', 'aa', seqs);assignin('base', 'aa1', trans);assignin('base', 'aa2', emis);
parfor kk=1:Nq
    disp(['Training the model for the condition ',num2str(kk)]);
    disp(['Number of trajectories, Number of time points: ',num2str(size(seqs{kk}(:,1:OtherPara.TimeWiseNum)))]);
    %
     %NumOfTraj=size(X{kk},2); %Used # of trajs
    %seqs=cell(1,NumOfTraj);
    %for jj=1:NumOfTraj
    %[N,edgesNew,seqs{kk}(jj,:)] = histcounts(X{kk}(:,jj)',edges);
    %seqs{kk}(jj,:)=seqs{kk}(jj,:);
   % end
    
[estTRTemp,estETemp] = hmmtrain(seqs{kk}(:,1:OtherPara.TimeWiseNum),trans,emis,'maxiterations',maxiter,'algorithm','BaumWelch');
 %disp('Done');
% later make it as that from fitting to all trajectories
estTR{1,kk}=estTRTemp;
estE{1,kk}=estETemp;
warning('off');
end


if OtherPara.HMMAllCondition==1
    countTraj=1;%seqs2test=0;
for kk=1:Nq
    disp('Fit to all condition and use it to get mutual information');
     NumOfTraj=size(X{kk},2); 
    for jj=1:NumOfTraj        
         %disp(size(X{kk}(:,jj)'));
    % [N,edgesNew,seqs2test] = histcounts(X{kk}(:,jj)',edges);
       
    [N,edgesNew,seqs2(countTraj,:)] = histcounts(X{kk}(:,jj)',edges);
    countTraj=countTraj+1;
    %seqs{kk}(jj,:)=seqs{kk}(jj,:);
    end
    disp(size(seqs2));
end
disp('start train HMM by all conditions');
[estTRTemp,estETemp] = hmmtrain(seqs2,trans,emis,'maxiterations',maxiter,'algorithm','BaumWelch');
 disp('Done');
estTR{1,Nq+1}=estTRTemp;
estE{1,Nq+1}=estETemp;

toc;
end



end



