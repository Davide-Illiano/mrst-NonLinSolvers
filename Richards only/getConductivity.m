function K = getConductivity(p, theta, K_s, n)
    K_pos = K_s*theta.^(1/2).*(1 - (1 - theta.^(n./(n-1))).^((n-1)./n)).^2;
    
    K_neg = K_s;
    
    
    neg = p <= 0;
    K = K_neg.*neg + K_pos.*~neg;
   
  
end