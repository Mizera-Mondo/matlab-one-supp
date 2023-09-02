function se = NMSE(X, X_e)
%NMSE Calculate the normalized mean squared error between true values X and
%estimation X_e.
se = norm(X - X_e, "fro")/norm(X, "fro");
se = se.^2;
end

