function [fX, hRS] = condProb(X,K,noiselvl, varargin)
%  Use K-nearest neighbors to calculate: (d-dimensional) probability densities
%    - probability densities (n x n cell matrix [n= # of inputs])
%    - conditional entropies
%
% Thus, a [n1xn2] conditional probability matrix calculates the probability of getting a response 
% (a query, pulled from response set X{n1}), given another set of reference responses (from X{n2})
% From Selimkhanov, Taylor  2014 conditional_probabilities nested function
% in getCC

% Define number of dimensions (e.g. timepoints) and volume of unit sphere of dimension=d
d=size(X{1},1); %dimension
V=pi^(d/2)/gamma(d/2+1);

Nq=length(X); % # of inputs
fX=cell(Nq,Nq); % Conditional probabilities
hRS=nan(1,Nq); % Conditional entropies       

% Add a small amount of random noise to each condition (avoids divide-by-zero issues later)
for i = 1:Nq
    X{i}=X{i}+noiselvl*rand(size(X{i}));
end

for s1=1:Nq % Outer loop: "query points" for KNN algorithm; rows of hRS
    for s2=1:Nq  % Inner loop: reference points; columns of hRS          
        Nsamp=size(X{s2},2); 
        % Diagonal of fX (same-well comparisons) -> Probability densities + conditional entropies
        if s1==s2           
            [~,DistK]=annquery(X{s2},X{s1},K+1); % Calculate K-nearest distances for all points
            Dk=DistK(K+1,:)'; % Get distance to (k+1)th nearest neighbor for all samples
            fX{s1,s2}=(K)./(Nsamp*V*Dk.^d); % Formula (2.10) for probability distribution             
            
            if nargin>3 % (optional probabilities are provided)
                if numel(varargin{1})>0  % (ensure entropies are actually a vector)
                    Ftrue=varargin{1};
                    hRS(s1)=-sum(log2(Ftrue{s1,s2}(Ftrue{s1,s2}>eps)))/nnz((Ftrue{s1,s2}>eps));
                else 
                    hRS(s1)=-sum(log2(fX{s1,s2}(fX{s1,s2}>eps)))/nnz((fX{s1,s2}>eps));
                end
            else % Calculate conditional entropies using using eqn 2.12
                hRS(s1)=-sum(log2(fX{s1,s2}(fX{s1,s2}>eps)))/nnz((fX{s1,s2}>eps));
            end
        
        % Off-diagonal (different-well comparision) -> Probability densities only.
        else    
            if isequal(X{s2},X{s1})
                [~,DistK]=annquery(X{s2},X{s1},K+1);
                Dk=DistK(K+1,:)';                
                fX{s1,s2}=(K)./(Nsamp*V*Dk.^d); 
            else              
                [~,DistK]=annquery(X{s2},X{s1},K);
                Dk=DistK(K,:)';
                fX{s1,s2}=(K)./((Nsamp+1)*V*Dk.^d); 
            end              
        end
    end
end