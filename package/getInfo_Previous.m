function [f,Hr] = getInfo(q,Fcond,hRS)
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
Nq=length(Fcond);%number of conditions
F=cell(Nq,1); %Nq x Nq matrix
Wq=cell(Nq,1); %
for k=1:Nq
    F{k}=cat(2,Fcond{k,:})*q; %Fcond{k,:}=conditional probability of a cell's response due to the stimulus given all other responses, q= probability of all conditions
    Wq{k}=ones(size(F{k}))*q(k)/sum(F{k}>realmin);    %probability of a condition/number of responses in that condition          
end

% Combine all probabilities/weights
F_all=cat(1,F{:}); %marginalized probabilty of all responses
Wq=cat(1,Wq{:});%length(Wq) = total number of cells, sum (Wq) = 1

%F_all(F_all<=eps)=nan; 
% Unweight any zero-values probability estimates
Wq(F_all<=realmin) = 0;

% Calculate non-conditional entropy (via eqn 2.13)
Hr = -sum(log2(F_all).*Wq); % negation of the weighted sum of a response given all other responses,

% Entropy = difference of H(R) and H(R|S), summed over the probabilities of each response. 
f=Hr-hRS*q;    
%f=abs(Hr-hRS*q); 
%f=max(f,0);

% assignin('base', 'Hr', Hr);
% assignin('base', 'hRS', hRS);
% assignin('base', 'q', q);


    