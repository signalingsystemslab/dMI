%Function to calculate channel-capacity

function [f,TrajEntropy,TrajEntropyModelBased]= getInfoModelMarkov(q,fX_all)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% GETINFO calculates mutual information for:
% 
% q         a set of input probabilities (sum(q)=1)
% Fcond     an [nxn] conditional probability matrix
% hRS       a vector of conditional entropies
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
% Get rows of conditional probability matrix - one "query", all "references"
% Sum across this row, individually weighting each point by appropriate
% weight.
NumTime=size(fX_all,2);
Nq=size(fX_all,1);%number of conditions
TrajEntropy=zeros(1,NumTime);


%Hr=0;%hRSTemp=cell(1,Nq);
for kk=1:NumTime %Do path average step by step...
%display(kk);
    
% disp(Fcond2);disp(size(Fcond2));

for s1=1:Nq
Joint=fX_all{s1,kk};
%display(Joint);
Marginal=sum(Joint,2);
assignin('base', 'Joint', Joint);%dddd
%Method 1 and 2:
hRS(s1)=-sum(Joint(Joint>eps).*log2(Joint(Joint>eps)))...
    +sum(Marginal(Marginal>eps).*log2(Marginal(Marginal>eps)));

% %Method 3:
% hRS(s1)=-sum(log2(Joint(Joint>eps)))...
%     +sum(log2(Marginal(Marginal>eps)));


TrajEntropyModelBased{s1}(kk)=hRS(s1);
end
% disp(hRS);ddd
    sz = size(cat(3,fX_all{:,kk}));
    Sum_Joint = reshape(reshape(cat(3,fX_all{:,kk}),[],sz(3))*q,[sz(1:2) 1]); %weighted sum transition matrix for all stimulus

    temp=sum(cat(3,fX_all{:,kk}),2);
    sz = size(temp);
    Sum_Marginal = reshape(reshape(temp,[],sz(3))*q,[sz(1:2) 1]);%weighted sum marginal of transition matrix for all stimulus
 
    % trajector entropy for the sum transition matrix
    %Method 1:
    Hr=-sum(Sum_Joint(Sum_Joint>eps).*log2(Sum_Joint(Sum_Joint>eps)))...
    +sum(Sum_Marginal(Sum_Marginal>eps).*log2(Sum_Marginal(Sum_Marginal>eps)));

assignin('base', 'Sum_Joint', Sum_Joint);

%      %Method 2:
%     Hr=-sum(Joint(Sum_Joint>eps).*log2(Sum_Joint(Sum_Joint>eps)))...
%     +sum(Marginal(Sum_Marginal>eps).*log2(Sum_Marginal(Sum_Marginal>eps)));
% 
%     %Method 3:
%     Hr=-sum(log2(Sum_Joint(Sum_Joint>eps)))...
%     +sum(log2(Sum_Marginal(Sum_Marginal>eps)));




TrajEntropy(kk)=Hr-hRS*q;

end

%%
% Entropy = difference of H(R) and H(R|S), summed over the probabilities of each response. 
f=sum(TrajEntropy);
% f=sum(Hr)-hRS*q;
%f_total=sum(f*timeunit);
