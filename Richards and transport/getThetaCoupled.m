function theta = getThetaCoupled(p, c, theta_R, theta_S, alpha, n,a,b)

    theta_neg = theta_R + (theta_S - theta_R).*(1./(1 + abs(alpha*(1./(1-b.*log(c./(a) +1))).^(1).*p)).^n).^((n-1)/n);
    
    %theta_neg = theta_R + (theta_S - theta_R).*(1./(1 + abs(alpha.*abs(p))).^n).^(1-(1/n));
    
    theta_pos = theta_S;
    
    neg = p < 0;
    theta = theta_neg.*neg + theta_pos.*~neg; 
    
end