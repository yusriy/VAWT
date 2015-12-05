%Function to calculate linear regression coefficients a0 and a1

function [a0 a1] = linear_regression(x,y)

n = length(x);
a1 = (n*sum(x.*y) - sum(x)*sum(y))/(n*sum(x.^2) - (sum(x))^2);

a0 = mean(y) - a1*mean(x);

end