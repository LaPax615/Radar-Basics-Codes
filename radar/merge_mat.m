function vector_1 = merge_mat(matrice)
vector_1 = [];

% Merge the rows of fi
    for k = 1:length(matrice(:,1))
        vector_1 = [vector_1, matrice(k, :)];
    end
end