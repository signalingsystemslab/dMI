function [SampledTraj,DistanceTotal, estTR,estE,estTR_sc,estE_sc,OtherPara,seqs] = Partial_Train (X,OtherPara)


% Create new array, X, which holds nan-filtered data 
%Nq=length(X);% Y{1,j} is  j-th condition: with row time series, column trajs

%run paralle with different parameters.... 10:10:70 for now
%OtherPara.binsize=10;%:10:20;%70;
Nq=length(OtherPara.binsize);
%OtherPara.state=round(OtherPara.binsize.*OtherPara.stateRatio/OtherPara.binsizeRatio);

tic;
% OtherPara.state=20;
% OtherPara.binsize=40;
%OtherPara.MaxValue=4;
% disp(OtherPara.MaxValue);
% disp(OtherPara.MinValue);


maxiter=1000;
estTR=cell(1,Nq+1);
estE=cell(1,Nq+1);
SampledTraj=cell(1,Nq);
DistanceTotal=cell(1,Nq);
seqs=cell(1,Nq);edges=cell(1,Nq);


%disp('Fit to all traj to get a HMM and use it to sample');
NumOfTraj=size(X{1},2);
%assignin('base', 'NumOfTraj', NumOfTraj);
%OtherPara.NumTraj= NumOfTraj;
binsize=OtherPara.binsize;    
MaxValue=OtherPara.MaxValue;
MinValue=OtherPara.MinValue;
state=OtherPara.state;
TimeWiseNum=OtherPara.TimeWiseNum;
%OtherPara.TrajLength=length(X{1}(:,1));

warning('off');
% assignin('base', 'state', state);
% assignin('base', 'binsize', binsize);

for kk=1:Nq
Conversion(kk)=(MaxValue-MinValue)/binsize(kk);
%Conversion=Conversion;
edges{kk}=linspace(MinValue,MaxValue,binsize(kk)+1);
%disp( edges{kk});
%disp(state(kk));
%disp('bin');
disp([binsize(kk) state(kk)]);
for jj=1:NumOfTraj
[N,edgesNew,seqs{kk}(jj,:)] = histcounts(X{1}(:,jj)',edges{kk});
end
end


trans=rand(state(kk),state(kk));
emis=rand(state(kk),binsize(kk));%ones(state,binsize);
trans=trans./sum(trans,2);
emis=emis./sum(emis,2);




for kk=1:Nq
    
  %disp(Nq);ddd;
    
[estTRTemp,estETemp] = hmmtrain(seqs{kk}(:,1:TimeWiseNum),trans,emis,'maxiterations',maxiter,'algorithm','BaumWelch');
 disp('Done');
% later make it as that from fitting to all trajectories
estTR{1,kk}=estTRTemp;
estE{1,kk}=estETemp;
warning('off');

end

toc;
OtherPara.Conversion=Conversion;
OtherPara.edges=edges;
estTR_sc=1;estE_sc=1;




