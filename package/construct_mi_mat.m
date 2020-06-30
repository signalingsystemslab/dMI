function X= construct_mi_mat (Data, Names, varargin)
% Constructs the input matrix for jacknifeCC & getCC from data_table
%--------------------------------------------------------------------------
% INPUTS
% REQUIRED
%     Data:    table | rows = cells,  columns = features
%     Names:      cell array of variable names 
% OPTIONAL:
%     GrpVar:   
% OUTPUTS
%   X:            matrix | rows = features, columns = cells

p = inputParser;
addRequired(p, 'Data', @istable);
addRequired(p, 'Names', @(x) iscell(x) ||ischar(x)) ;
addOptional(p, 'GrpVar','ID', @ischar);
parse(p, Data,Names, varargin{:});
G= p.Results.GrpVar;
fun = @(x) {transpose(horzcat(x))};
if iscategorical(Data.(G)); Data.(G) =removecats(Data.(G)); end
expt_idx=grp2idx(Data.(G));
X= splitapply(fun,Data{:,Names},expt_idx)';

end


