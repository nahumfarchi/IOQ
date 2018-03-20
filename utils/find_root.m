function root_ = find_root(c0_, n_)
    M_ = zeros(n_, n_);
    for i_ = 1:(n_-1)
        M_(i_+1,i_) = 1;
    end
    M_(1, n_) = c0_;
    [~, D_] = eig(M_);
    D_ = sort(diag(D_));
    root_ = D_(1);
end