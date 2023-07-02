function [xAv,x2Av,yAv,y2Av,zAv,z2Av] = Br3DUniform(k,nTrials,N)

dt = 1;
T = N;
dt = dt/k;

N = k*N;
t = 0:dt:T; 
global x;
x = zeros(length(t),nTrials);
global y;
y = zeros(length(t),nTrials);
global z;
z = zeros(length(t),nTrials);

dr = 1/sqrt(k);

for j = 1:nTrials
    for i = 1:N
      random = rand;
        if random < 1/2
            x(i+1,j) = x(i,j) + dr;
            y(i+1,j) = y(i,j);
            z(i+1,j) = z(i,j);
        else
           if  random < 2/6
                z(i+1,j) = z(i,j)+dr;
                x(i+1,j) = x(i,j);
                y(i+1,j) = y(i,j);
                
                
            else
                if random < 1/2
                    y(i+1,j) = y(i,j) + dr;
                    x(i+1,j) = x(i,j);
                    z(i+1,j) = z(i,j);
                    
                else
                    if random < 4/6
                        x(i+1,j) = x(i,j) - dr;
                        z(i+1,j) = z(i,j);
                        y(i+1,j) = y(i,j);
                    else
                        if random < 5/6
                            z(i+1,j) = z(i,j) - dr;
                            y(i+1,j) = y(i,j);
                            x(i+1,j) = x(i,j);
                        else
                           
                            y(i+1,j) = y(i,j) - dr;
                            x(i+1,j) = x(i,j);
                             z(i+1,j) = z(i,j);
                            
                        end
                    end
                end
            end
            
         end
    end
end

yAv = mean(y(end,:));
y2Av = mean(y(end,:).^2);

xAv = mean(x(end,:));
x2Av = mean(x(end,:).^2);


zAv = mean(z(end,:));
z2Av = mean(z(end,:).^2);
end