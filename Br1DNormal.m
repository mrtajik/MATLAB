function [xAv,x2Av] = Br1DNormal(k,nTrials,N)

dt = 1;
T = N;
dt = dt/k;

N = k*N;

global t;
t = 0:dt:T; 

global x;
x = zeros(length(t),nTrials);

dx = randn*(1/sqrt(k));
for j = 1:nTrials
    for i = 1:N
        
        if rand < 0.5
            x(i+1,j) = x(i,j) + dx;
            
        else
              x(i+1,j) = x(i,j) - dx;
        end 
       
    end
 
end

xAv = mean(x(end,:));
x2Av = mean(x(end,:).^2);

end