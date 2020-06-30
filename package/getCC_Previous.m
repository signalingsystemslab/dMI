function [I,Q,fX] = getCC (Y,K,idisp,varargin)
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
Nq=length(Y);
X=cell(Nq,1);
for i=1:Nq
    A=Y{i};
    if size(A,1)==1 % Unidimensional case: remove NaN individuals
        A=A(:,~isnan(A));
        A=A(:,~isinf(A));
    else % Multidimensional case: remove individuals with NaN in any dimension
        A=A(:,sum(isnan(A))==0);
        A=A(:,sum(isinf(A))==0); 
    end
    X{i}=A;   
end
%assignin('base', 'Z', X);

% Calculate conditional probabilities and entropies - if getting infinite entropy, increase added noise
do_flag = 1;
noiselvl = eps;
while do_flag
    [fX, hRS] = conditional_probabilities(X,K,noiselvl,varargin);
    suminf = 0;
    for i = 1:numel(fX)
        suminf = suminf+sum(isinf(fX{i}(:)));
    end
    if suminf==0
        do_flag = 0;
    else
        noiselvl = noiselvl*10;
    end
end


% If probabilities are provided, assign into fX
if nargin>3 
    if numel(varargin{1})>0    
        fX=Ftrue;    
    end
end
% assignin('base', 'fX', fX);
% assignin('base', 'hRS', hRS);
% assignin('base', 'Xafter', X);
% (STEP 2) Optimize input probabilities for maximum channel capacity
% Channel capacities are calculated as the difference between the sum of conditional entropies H(R|S),
% and the total entropy of the response, H(R).
%
%     - H(R|S) = sum(hRS);
%     - H(R) = doubly-weighted sum (see getInfo.m) 


% Set options (toggle verbose/non-verbose display)
if idisp==0
    options = optimoptions('fmincon','Algorithm','active-set','Display','off');
elseif idisp==1
    options = optimoptions('fmincon','Algorithm','active-set','Display','iter');
end

% Define information calculation function (getInfo.m), initialize bounds for optimization.
infoF = @(q) -getInfo_Previous(q,fX,hRS);
Aeq=ones(1,Nq);
beq=1;
LB=zeros(Nq,1);%+eps;%+realmin;
UB=ones(Nq,1);  

if nargin>4 % q's are provided; skip optimization and output.
    if numel(varargin{2})>0
        Qstart=varargin{2};
        I=-infoF(Qstart);  
    end
else % Optimize q for information content, given probability distributions 
    Qstart=ones(Nq,1);
    Qstart=Qstart/sum(Qstart);    
    qfit=fmincon(infoF,Qstart,[],[],Aeq,beq,LB,UB,[],options);
    Q=qfit;
    I=-infoF(qfit);     
end





function [fX, hRS] = conditional_probabilities(X,K,noiselvl, varargin)
%  Use K-nearest neighbors to calculate: (d-dimensional) probability densities
%    - probability densities (n x n cell matrix [n= # of inputs])
%    - conditional entropies
%
% Thus, a [n1xn2] conditional probability matrix calculates the probability of getting a response 
% (a query, pulled from response set X{n1}), given another set of reference responses (from X{n2})


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
            assignin('base', 'Dk', Dk);
            assignin('base', 'DistK', DistK);
            assignin('base', 'fX', fX);
            %disp(nargin);
            if nargin>5 % (optional probabilities are provided)
                if numel(varargin{1})>0  % (ensure entropies are actually a vector)
                    Ftrue=varargin{1};
                    hRS(s1)=-sum(log2(Ftrue{s1,s2}(Ftrue{s1,s2}>realmin)))/nnz((Ftrue{s1,s2}>realmin));
                else 
                    hRS(s1)=-sum(log2(fX{s1,s2}(fX{s1,s2}>realmin)))/nnz((fX{s1,s2}>realmin));
                end
            else % Calculate conditional entropies using using eqn 2.12
                 if nnz((fX{s1,s2}>realmin))==0
                     hRS(s1)=0;%My revision to the code
                 else
                     hRS(s1)=-sum(log2(fX{s1,s2}(fX{s1,s2}>realmin)))/nnz((fX{s1,s2}>realmin));
                 end
                  %hRS(s1)=-sum(log2(fX{s1,s2}(fX{s1,s2}>eps)))/nnz((fX{s1,s2}>eps));
                assignin('base', 'fX', fX);
                assignin('base', 'hRS', hRS);
%                 dddd
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



