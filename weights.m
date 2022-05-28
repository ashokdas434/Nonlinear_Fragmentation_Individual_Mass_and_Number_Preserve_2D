%% Function to formulate the weights
function [w1,w2_b,w2_d] = weights(x1,x2,B)

I1 = length(x1); I2 = length(x2);
w1 = zeros(I1,I2); w2_b = zeros(I1,I2); w2_d = zeros(I1,I2); % Initialization
nu = 2; % number of fragments per breakage

for i=1:I1
    for j=1:I2
        a=1:i; b=1:j;
        B_temp = B(a,b,i,j);
        temp1 = x1(a)* sum(B_temp,2); temp2 = sum(B_temp,1)*x2(b)';

        w1(i,j) = (temp1+temp2) /(x1(i)+x2(j));

        w2_b(i,j) = (x1(i)+x2(j))*(nu-1)/ ( (x1(i)+x2(j))*sum(B_temp(:))-temp1 -temp2 );
        w2_d(i,j) = w2_b(i,j) * w1(i,j);
    end
end
w2_b(1,1) = 0; w2_d(1,1) = 0;

return