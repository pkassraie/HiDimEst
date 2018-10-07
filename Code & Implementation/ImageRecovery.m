%% 
% Read Image, Convert to grayscale
X_col = imread('X.jpg');
X = double(rgb2gray(X_col));
Sampling = zeros(size(X));
[d1,d2] = size(X);
r = 50;

i=1;
s = 2;
while i<=d1
    j=1;
    while j<=d2
        Sampling(i,j) = 1;
        j = j+8;
    end
    i = i+ 8;
end
obs = X;
obs(Sampling==1) = 1;
obs(Sampling==0) = 0;
ind = find(Sampling==1);
imagesc(X)

% Projection
[U,S,V] = svd(obs);
P = ones([r,r]);
P = padarray(P,[d1-r d2-r],0,'post');
S = S.*P;
obs = U*S*V';
%obs(ind)=X(ind);  
figure()
imshow(uint8(obs))
snr(j) = 10*(log10(norm(X)/norm(X-obs)));

%%


X_col = imread('X.jpg');
X = double(rgb2gray(X_col));
Sampling = zeros(size(X));
[d1,d2] = size(X);
rank = 50;

i=1;
s = 2;
while i<=d1
    j=1;
    while j<=d2
        Sampling(i,j) = 1;
        j = j+8;
    end
    i = i+ 8;
end
obs = X;
obs = obs.*Sampling;
ind = find(Sampling==1);


% Projection
[U,S,V] = svd(obs);
P = ones([r,r]);
P = padarray(P,[d1-r d2-r],0,'post');
S = S.*P;
obs = U*S*V';

%obs(ind)=X(ind);  

imagesc(obs)

%% 
% Read Image, Convert to grayscale
X_col = imread('X.jpg');
X = double(rgb2gray(X_col));
[d1,d2] = size(X);
r = 50;
snr = [0 0 0 0 0];
M = [10, r*(d1+d2), d1*d2];
% Set Rank

for j=1:3
% Create random observation vector, uniform
m = M(j);
prob = m/(d1*d2);
mask =rand(size(X),'like',X);
obs = X;
obs(mask<prob) = 1;
obs(mask>=prob) = 0;
ind = find(mask==1);


% Projection

[U,S,V] = svd(obs);
P = ones([r,r]);
P = padarray(P,[d1-r d2-r],0,'post');
S = S.*P;
obs = U*S*V';
%obs(ind)=X(ind);  
figure()
imshow(uint8(obs))
snr(j) = 10*(log10(norm(X)/norm(X-obs)));
end





