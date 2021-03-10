function[out] = CoxLeeTest(Y1, Y2, X, confLevel, nRandomize)
    length_X = size(X,1);
    p_values = zeros(length_X,1);
    for i = 1:length_X
        [~ , p_values(i)] = ttest2(Y1(i,:)',Y2(i,:)');
    end
    [sorted_p_values, sort_index] = sort(p_values);
    qstar = zeros(length_X,nRandomize);
    for l = 1:nRandomize
        pstar = zeros(length_X,1);
        for j = 1:length_X
            [y1star, y2star] = permute(Y1(j,:)',Y2(j,:)');
            [~, pstar(j)] = ttest2(y1star,y2star);
        end
        sorted_pstar = pstar(sort_index);
        for j = 1:length_X
            qstar(j,l) = min(sorted_pstar(j:end));
        end
    end
    r = zeros(length_X,1);
    for j = 1:length_X
        r(j) = (1/nRandomize)*length(find(qstar(j,:)<sorted_p_values(j)));
    end
    alpha = 1 -confLevel;
    jPlus = find(r > alpha, 1);
    if jPlus > 1
        out.differ = true;
        out.statSigPointIndices = sort(sort_index(1:(jPlus-1)));
    else
        out.differ = false;
        out.statSigPointIndices = [];
    end
    [~,r_order] = sort(sort_index);
    out.corrected_p_values = r(r_order);
end

function[permuted_y1, permuted_y2] = permute(y1,y2)
    if size(y1,2) ~=1 || size(y2,2) ~=1
        error('y1 and y2 must be column vectors');
    end
    length_y1 = length(y1);
    length_y2 = length(y2);
    combined_vector = [y1;y2]; %assuming column vectors
    permuted_vector = combined_vector(randperm(length(combined_vector)));
    permuted_y1 = permuted_vector(1:length_y1);
    permuted_y2 = permuted_vector(length_y1+1:length_y1+length_y2);
end
