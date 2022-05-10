function [theta] = optimal_theta(y,zeta)
 sum_den = 0;
 sum_num = 0;
    for i = 1:length(y(:,1))
        sum_num = (sum_num +  (y(i,1))'*zeta(i,:));
        sum_den = sum_den + zeta(i,:)'*zeta(i,:);
    end
    
    theta = (sum_num / sum_den)';

end

