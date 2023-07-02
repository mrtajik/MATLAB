
% Assignment 9

% 1.Random numbers drawn from a uniform distribution:
% a) In MATLAB generate N=10000 random numbers sampled from a uniform distribution between 0 to 1 (use the command rand).
clear;
clc;
X = rand(10000,1);

% b) Plot a histogram for these random samples.
h = histogram(X);
title('Histogram Part b');

% c) Normalize the histogram such that it estimates the probability density function of your samples.
figure
histogram(X,'normalization','pdf')
[f,xi] = ksdensity(X);
hold on
plot(xi,f)
title('Graph of Part c')


% d) Find the mean and variance (use commands mean and var).
M = mean(X);
V = var(X);
fprintf('Mean = %d\n', M);
fprintf('Variance = %d\n', V);
%% 

% 2.Random numbers drawn from a normal distribution:
% a) In MATLAB generate N=10000 random numbers sampled from a normal distribution (use the command randn).
clear;
clc;
Y = randn(10000,1);

% b) Plot a histogram for these random samples.
h = histogram(Y);
title('Histogram Part b');

% c) Normalize the histogram such that it estimates the probability density function of your samples.
figure
histogram(Y,'normalization','pdf')
[f,xi] = ksdensity(Y);
hold on
plot(xi,f)
title('Graph of Part c')


% d) Find the mean and variance (use commands mean and var).
M = mean(Y);
V = var(Y);
fprintf('Mean = %d\n', M);
fprintf('Variance = %d\n', V);

%% 

% 3.Exploring the Central Limit Theorem:
% a) Generate 5 different set(X1,X2,X3,X4,X5) of N=10,000 random numbers from a uniform distribution between -1 to 1.
clc;
clear;
N = 10000;
X1= unifrnd(-1 , 1 ,N,1);
X2= unifrnd(-1 , 1 ,N,1);
X3= unifrnd(-1 , 1 ,N,1);
X4= unifrnd(-1 , 1 ,N,1);
X5= unifrnd(-1 , 1 ,N,1);

% b) Build the sample average X.
SAMPLE_AVARAGE_X = (X1+X2+X3+X4+X5)/5;

% c) Build the normalized sample Z.
Z = normalize(SAMPLE_AVARAGE_X);

% d) Plot a histogram for Z.
histogram(Z)
title('Histogram Part d');


