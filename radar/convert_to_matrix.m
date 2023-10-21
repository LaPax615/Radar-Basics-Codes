function fi = convert_to_matrix(fi_0, K)
    % Initialize fi as a matrix of size K x length(fi_0)
    fi = zeros(K, length(fi_0));

    % Fill each row of fi with fi_0
    for k = 1:K
        fi(k, :) = fi_0;
    end
end
