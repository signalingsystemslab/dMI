function [nfkb, all_dims, names_1D] = loadnfkb()
%%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [nfkb, all_dims, names_1D] = loadnfkb()
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% LOADNFKB loads the master set of NFkB runs used for (most) BMDM analyses
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% If nfkb.mat doesn't exist, load, measure, and save the individual experiments
ids = [
274 % ctrl

290	% 2015-07-01	0.33ngTNF
291	% 2015-07-01	1ngTNF
292	% 2015-07-01	3.3ngTNF
293	% 2015-07-01	10ngTNF
294	% 2015-07-01	33ngTNF

300	% 2015-07-13	0.33ngLPS
289	% 2015-06-30	1ngLPS
302	% 2015-07-13	3.3ngLPS
303	% 2015-07-13	10ngLPS
304	% 2015-07-13	33ngLPS
311 % 2015-07-23    100ngLPS
310 % 2015-07-23    333ngLPS

305	% 2015-07-14	3.3ug_polyIC 
295	% 2015-07-02	10ug_polyIC
296	% 2015-07-02	33ug_polyIC
306	% 2015-07-14	100ug_polyIC 

281 % 2015-06-29	10nM_CpG
282	% 2015-06-29	33nM_CpG
283	% 2015-06-29	100nM_CpG
284	% 2015-06-29	333nM_CpG
285	% 2015-06-29	1uM_CpG

326	% 2015-09-03	1ngP3CSK
325	% 2015-09-03	3.3ngP3CSK
327 % 2015-09-03	10ngP3CSK
309 % 2015-07-14	33ngP3CSK
298 % 2015-07-02	100ngP3CSK
];

% Define name of nfkb file
P = mfilename('fullpath');
P2 = mfilename;
nfkb_name = [P(1:(length(P)-length(P2))) , 'nfkb.mat'];

% Combine processed data into a single structure
if ~exist(nfkb_name,'file')
    nfkb = struct;
    for i =1:length(ids)
        [nfkb(i).metrics, nfkb(i).aux] = nfkbmetrics(ids(i));
        nfkb(i).id = ids(i);
        nfkb(i).ids = ids;
    end
    save(nfkb_name, 'nfkb')
else
    load(nfkb_name);
end

% % ADD windows of integration
% if ~isfield(nfkb(1).metrics,'intwin1')
%     for i = 1:length(nfkb)
%         if isfield(nfkb(i).metrics,'intwin_1')
%             nfkb(i).metrics = rmfield(nfkb(i).metrics,'intwin_1');
%             nfkb(i).metrics = rmfield(nfkb(i).metrics,'intwin_3');
%         end
%         
%         nfkb(i).metrics = addwindows(nfkb(i).metrics);
%     end
%     save(nfkb_name, 'nfkb');
% 
% end


% If extra output arguments are defined, combine all metrics/measurements into one big matrix 
% (record dimension names as well)
if nargout~=1
    mfields = fieldnames(nfkb(1).metrics);
    names_1D = {};
    nfkb(ids==281) = []; % Drop the "off" CpG set

    for i = 1:length(nfkb)
        % Downsample time series, integrals, and derivatives slightly
        nfkb(i).metrics.time_series = nfkb(i).metrics.time_series(:,[1:36,37:2:180]);
        nfkb(i).metrics.integrals = nfkb(i).metrics.integrals(:,1:2:240);
        nfkb(i).metrics.derivatives = nfkb(i).metrics.derivatives(:,[1:36,37:2:180]);

        % Combine all dimensions into a single gigantic cell thing
        all_dims{i} = [];
        for j = 1:length(mfields);
            all_dims{i} = cat(2,all_dims{i},nfkb(i).metrics.(mfields{j}));
            if i==1
                for k = 1:size(nfkb(i).metrics.(mfields{j}),2)
                    names_1D = cat(1,names_1D,[mfields{j},'_',numseq(k,3)]);
                end
            end
        end
    end

    % Reload full nfkb.mat
    load(nfkb_name)
end

% If called without output, assign everything into base workspace
if nargout<1
    assignin('base','nfkb',nfkb)
    assignin('base','all_dims',all_dims)
    assignin('base','names_1D',names_1D)
end

