N=2;
omega = 0.5;
b0 = 0.25;
v = sqrt(pi/2)/6;
const = 1; % ros/a
myeta = (3/const)*sqrt((sqrt(pi)*erf(v))/(2*v*exp(-(v^2))));

alpha = 10;
beta = 4;
Pm = 0.75;
%myeta = 0.4;
g = 2*b0*(1-Pm);
omegaBar = omega + 2*b0*Pm + 2*sqrt(2*b0*Pm*omega);
%u = 1:0.5:40;
u_db = 1:0.5:40;
u = 10.^(u_db/10);
avgber = zeros(1, length(u));

for i2 = 1:length(u)
    A = ((2*alpha^(alpha/2))/((g^(1+alpha/2))*gamma(alpha)))*(g*beta/(g*beta + omegaBar))^(beta+alpha/2);
    product = 1;
    for i=1:N
        sum = 0;
        for k=1:beta
            
            a_k = nchoosek(beta-1,k-1)*(((g*beta + omegaBar)^(1-k/2))/factorial(k-1))*((omegaBar/g)^(k-1))*((alpha/beta)^(k/2));
            z = ((alpha*beta*k*(g+omegaBar)/(g*beta + omegaBar))^2)/(16*u(i2)); 
            temp1 = ((myeta*myeta*A*a_k*(2^(alpha+k-4)))/pi)*((alpha*beta/(g*beta+omegaBar))^(-alpha/2-k/2));
            %temp2 = meijerG([1 (2+myeta*myeta)/2],[],[myeta*myeta/2 alpha/2 (1+alpha)/2 k/2 (1+k)/2],[], z);
            temp2 = meijerG([1],[(2+myeta*myeta)/2],[myeta*myeta/2, alpha/2 , (1+alpha)/2, k/2 , (1+k)/2],[], z);
            sum = sum + (temp1 * temp2);
            
        end
        product = product * sum;
    end
    avgber(:,i2) = product * 0.5;
end

disp('Done 1');

alpha = 10;
beta = 4;
Pm = 0.75;
%myeta = 0.4;
g = 2*b0*(1-Pm);
omegaBar = omega + 2*b0*Pm + 2*sqrt(2*b0*Pm*omega);
%u = 1:0.5:40;
u_db = 1:0.5:40;
u = 10.^(u_db/10);
avgber1 = zeros(1, length(u));

for i2 = 1:length(u)
    A = ((2*alpha^(alpha/2))/((g^(1+alpha/2))*gamma(alpha)))*(g*beta/(g*beta + omegaBar))^(beta+alpha/2);
    product = 1;
    for i=1:N
        if (i==2)
            e = u(i2)/2;
        else
            e = u(i2);
        end
        sum = 0;
        for k=1:beta
            
            a_k = nchoosek(beta-1,k-1)*(((g*beta + omegaBar)^(1-k/2))/factorial(k-1))*((omegaBar/g)^(k-1))*((alpha/beta)^(k/2));
            z = ((alpha*beta*k*(g+omegaBar)/(g*beta + omegaBar))^2)/(16*e); 
            temp1 = ((myeta*myeta*A*a_k*(2^(alpha+k-4)))/pi)*((alpha*beta/(g*beta+omegaBar))^(-alpha/2-k/2));
            %temp2 = meijerG([1 (2+myeta*myeta)/2],[],[myeta*myeta/2 alpha/2 (1+alpha)/2 k/2 (1+k)/2],[], z);
            temp2 = meijerG([1],[(2+myeta*myeta)/2],[myeta*myeta/2, alpha/2 , (1+alpha)/2, k/2 , (1+k)/2],[], z);
            sum = sum + (temp1 * temp2);
            
        end
        product = product * sum;
    end
    avgber1(:,i2) = product * 0.5;
end

disp('Done 2');
alpha = 10;
beta = 4;
Pm = 0.75;
%myeta = 0.4;
g = 2*b0*(1-Pm);
omegaBar = omega + 2*b0*Pm + 2*sqrt(2*b0*Pm*omega);
%u = 1:0.5:40;
u_db = 1:0.5:40;
u = 10.^(u_db/10);
avgber2 = zeros(1, length(u));

for i2 = 1:length(u)
    A = ((2*alpha^(alpha/2))/((g^(1+alpha/2))*gamma(alpha)))*(g*beta/(g*beta + omegaBar))^(beta+alpha/2);
    product = 1;
    for i=1:N
        if (i==2)
            e = u(i2)/4;
        else
            e = u(i2);
        end
        sum = 0;
        for k=1:beta

            a_k = nchoosek(beta-1,k-1)*(((g*beta + omegaBar)^(1-k/2))/factorial(k-1))*((omegaBar/g)^(k-1))*((alpha/beta)^(k/2));
            z = ((alpha*beta*k*(g+omegaBar)/(g*beta + omegaBar))^2)/(16*e); 
            temp1 = ((myeta*myeta*A*a_k*(2^(alpha+k-4)))/pi)*((alpha*beta/(g*beta+omegaBar))^(-alpha/2-k/2));
            %temp2 = meijerG([1 (2+myeta*myeta)/2],[],[myeta*myeta/2 alpha/2 (1+alpha)/2 k/2 (1+k)/2],[], z);
            temp2 = meijerG([1],[(2+myeta*myeta)/2],[myeta*myeta/2, alpha/2 , (1+alpha)/2, k/2 , (1+k)/2],[], z);
            sum = sum + (temp1 * temp2);
            
        end
        product = product * sum;
    end
    avgber2(:,i2) = product * 0.5;
end
disp('DONE !!!');

semilogy(u_db, avgber, u_db, avgber1, u_db, avgber2);
legend('mu1 = mu2','mu2 = mu1/2','mu2 = mu1/4');
title('Unbalanced SNR');
xlabel('SNR (mu1) in dB');
ylabel('BER');
grid on;