function theta = getTheta(p,  theta_R, theta_S, alpha, n)
    theta_neg = theta_R + (theta_S - theta_R).*(1./(1 + alpha*abs(p).^n)).^((n-1)/n);
    %theta_neg = theta_R;
    theta_pos = theta_S;
    
    neg = p <= 0;
    theta = theta_neg.*neg + theta_pos.*~neg; 
    
    %% In case of Salinity problem we need a and b as parameters 
    %m = 1-1/n;
    %theta = (theta_S - theta_R)./(1+abs(alpha.*(1-b*log(c./a+1)).^(-1).*p).^n).^m +theta_R;
    %theta = (theta_R - theta_S)./(1+abs(-alpha*(1-b*log(c./a+1)).*p).^n).^m +theta_R;
    
end