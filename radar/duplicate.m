function vector_f = duplicate(vector_0, K)
    % Initialize fi as a matrix of size K x length(fi_0)
    matrice = convert_to_matrix(vector_0, K);

    % Initialize fi_combined as an empty row vector
    vector_f = merge_mat(matrice);
end
