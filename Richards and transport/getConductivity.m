function K = getConductivity(p, theta, K_s, n)   % i added K_s
    K_pos = K_s.*theta.^(1/2).*(1 - (1 - theta.^(n./(n-1))).^((n-1)./n)).^2;
   
    m = 1-(1/n);
    alpha = 0.0335;
    K_pos = ((1 - (alpha .* abs(p)) .^ (n - 1) ...
              .* (1 + (alpha .* abs(p)) .^ n) .^ (-m)) .^ 2 ...
              ./ (1 + (alpha .* abs(p)) .^ n) .^ (m ./ 2));
    
    K_neg = K_s;
    
    
    neg = p <= 0;
    K = K_neg.*neg + K_pos.*~neg;
   
end