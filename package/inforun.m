function info = inforun(all_dims, names, high_dim, info)
%%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% info = inforun(all_dims, names, savename)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% INFORUN performs information optimizing over an arbitrary set of dimensions, specified as columns in 'all_dims', and
% strings in 'names'. Vectors of size up to the value specified by 'high_dim' will be calculated.
%
% (Search depth = 24 candidates)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin<3
    high_dim = 9;
end

if nargin<4
    info = [];
    low_dim = 1;
else
    low_dim = length(info);
end

depth = 16; % Candidates to carry forawrd at each dimension increase
sz = size(all_dims{1},2); % Total number of metrics


% Fit the difference between the top set's increase, and the maximal increase of any set, by dimension
delta_I = (1/sqrt(2)).^(1:6) - (1/2).^(1:6);
delta_I(1) = 0.3;
delta_I = [delta_I, 0.1*ones(1,high_dim-length(delta_I))];

for dim = low_dim:high_dim
    if (low_dim==1) || (dim>low_dim) 
        % 1-D scan: cycle all dimensions, and calculate information capacity 
        if dim==1
            all_I = zeros(sz,1);
            all_Q = zeros(sz,length(all_dims));
            all_names = names;
            parfor dim1 = 1:sz
                tic
                X= cell(size(all_dims));
                for expt =1:length(all_dims)
                    X{expt} = all_dims{expt}(:,dim1)';
                end
                [fitI,I,~,Q] = jacknifeCC(X,9,12);
                if abs(fitI-I)> 0.2 % (ensure that output I isn't an artifact of jacknife process)
                    fitI = I;
                end
                all_I(dim1) = fitI;
                all_Q(dim1,:) = Q;
                t = toc;
                disp(['dim. ',num2str(dim1),': ',names{dim1},'. I = ',num2str(fitI), ' (t = ',num2str(round(t*100)/100),' sec)'])
            end

        % Multidimensional scan: use top dimensions from previous iteration, scan against remaining dimensions    
        else
            % Initialize multidimensional matrix
            all_I = zeros(depth*sz,1);
            all_Q = zeros(depth*sz,1,length(all_dims));
            all_names = cell(depth*sz,1);
            names = info(1).names;

            % Make a list of all combinations of dimensions
            dim_idx = repmat(top_dims,[sz,1]);
            tmp = [];
            for i = 1:sz
                tmp = cat(1,tmp,i*ones(depth,1));
            end
            dim_idx = [dim_idx,tmp];

            parfor j = 1:size(dim_idx,1)
                if numel(dim_idx(j,:))==numel(unique(dim_idx(j,:))) % Skip any vector with redundancies
                    tic
                    X= cell(size(all_dims));
                    for expt =1:length(all_dims)
                        X{expt} = all_dims{expt}(:,dim_idx(j,:))';
                    end
                    [fitI,I,~,Q] = jacknifeCC(X,9,12)
                    if abs(I-fitI) > 0.2 % (ensure that output I isn't an artifact of jacknife process)
                        fitI = I;
                    end
                    all_I(j) = fitI;
                    all_Q(j,1,:) = Q;
                    str = '';
                    for i = 1:length(dim_idx(j,:))
                        str = [str,names{dim_idx(j,i)},'+'];
                    end
                    all_names{j} = str(1:end-1);
                    t = toc;
                    disp([all_names{j},'. I = ',num2str(fitI), ' (t = ',num2str(round(t*100)/100),' sec)'])
                end
            end
        end

        % Combine all data into a single structure, and save
        if dim==1
            info(dim).I = all_I;
            info(dim).Q = all_Q;
            info(dim).names = all_names;
            info(dim).idx{1} = (1:sz)';
        else
            info(dim).I = reshape(all_I,[depth,sz]);
            info(dim).Q = reshape(all_Q,[depth,sz,length(all_dims)]);
            info(dim).names = reshape(all_names,[depth,sz]);
            info(dim).idx{1} = top_dims;
            info(dim).idx{2} = (1:sz)';
        end
            save info.mat info
    end

    % Identify the subset of dimensions that will be able to "catch up" to the top-scoring set
    % Within this subset, identify the most orthogonal - move these forward to the next round.
    info(dim).I(isnan(info(dim).I)) = 0;
    [sorted_I, order] =sort(info(dim).I(:),'descend');
    order = order(sorted_I > (sorted_I(1) - delta_I(dim)));
    [r, c] = ind2sub(size(info(dim).I),order);

    if dim==1
        top_Q = squeeze(info(1).Q(r(1),:))';
        top_dims = info(dim).idx{1}(r(1));
    else
        top_Q = squeeze(info(dim).Q(r(1),c(1),:));
        top_dims = [info(dim).idx{1}(r(1),:),c(1)];

    end
    top_order = order(1);

    for top_ind = 1:(depth-1)
        radius = 0;
        for i = 2:length(r)
            if dim==1
                test_dim = info(dim).idx{1}(r(i));
                test_Q = squeeze(info(1).Q(r(i),:))';
            else
                test_dim = [info(dim).idx{1}(r(i),:),c(i)];
                test_Q = squeeze(info(dim).Q(r(i),c(i),:));
            end
            if ~ismember(order(i),top_order);
                tmp = zeros(1,size(top_Q,2));
                for j = 1:size(top_Q,2)
                    tmp(j) = norm(top_Q(:,j)- test_Q);
                end
                if (min(tmp) > radius)
                    radius = min(tmp);
                    new_Q = test_Q;
                    new_dim = test_dim;
                end
                new_order = order(i);
            end
        end
       if radius > 0
           disp(num2str(radius))
           top_dims = [top_dims;new_dim];
           top_Q = [top_Q, new_Q];
           top_order = [top_order; new_order];
       end
    end
    
  
    disp(['Dimension ', num2str(dim),' complete. Moving ON...'])
end
        