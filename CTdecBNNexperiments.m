%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pure decomposition using BNN
%
% Author: Shunsuke Ono (ono@sp.ce.titech.ac.jp)
% Last version: Mar 16, 2014
% Article: S. Ono, T. Miyata, and I. Yamada,
% "Cartoon-Texture Image Decomposition Using Blockwise Low-Rank Texture Characterization,"
% IEEE Transactions on Image Processing, vol. 23, no. 3, pp. 1128-1142, 2014.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
addpath subfunctions
addpath images

%%%%%%%%%%%%%%%%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%
imname = 'Barbara256.png';
blocksize = 16; % the block size of BNN
shiftstep = 8; % the shift step number of BNN
%-------shear parameter settings-------
% theta: 0 to 45 degrees
% direction: 'r', 'l', 't', or 'b'
theta{1} = 0;
direction{1} = [];
theta{2} = 45;
direction{2} = 'r';
theta{3} = 45;
direction{3} = 'l';
%---------------------------------
lambda = 0.3; % the weight of TV term
gamma = 0.1; % the step size of ADMM

%===============================================================
% comment out one 'problemtype' and the corresponding parameters
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
problemtype = 'PureDec'; % pure decomposition
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% problemtype = 'Blur+Miss'; % deblurring with missing pixels and noise
% psfsize = 3; % the size of Gaussian blur kernel
% missrate = 0.2; % the rate of missing pixels
% sigma = 0.05; % noise standard deviation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%===============================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_org = double(imread(imname))/255;
n = size(u_org);
N = prod(n);
K = numel(theta); % the number of subtexture components

%% definitions

% difference operators
D = @(z) cat(3, z([2:n(1), n(1)],:) - z, z(:,[2:n(2), n(2)])-z);
Dt = @(z) [-z(1,:,1); - z(2:n(1)-1,:,1) + z(1:n(1)-2,:,1); z(n(1)-1,:,1)] ...
    +[-z(:,1,2), - z(:,2:n(2)-1,2) + z(:,1:n(2)-2,2), z(:,n(2)-1,2)];

% periodically expanding operators
P = @(z) PeriodicExpansion(z, blocksize, shiftstep);
Pt = @(z) PeriodicExpansionTrans(z, shiftstep);

% shear operators
S = cell(K,1);
St = cell(K,1);
for j = 1:K
    if theta{j} == 0
        S{j} = @(z) z;
        St{j} = @(z) z;
    else
        S{j} = @(z) Shear(z, theta{j}, direction{j});
        St{j} = @(z) ShearTrans(z, theta{j}, direction{j});
    end
end

% observation operator
switch problemtype
    case 'PureDec'
        Phi = @(z) z;
        Phit = @(z) z;
        u_obsv = u_org;
        epsilon = 0;
    case 'Blur+Miss'
        blu = zeros(n);
        blu(1:psfsize, 1:psfsize) = fspecial('gaussian', [psfsize psfsize], 1);
        blu = circshift(blu, [-(psfsize-1)/2 -(psfsize-1)/2]);
        h = fft2(blu);
        ht = conj(h);
        h = repmat(h, [1 1]);
        ht = repmat(ht, [1 1]);
        B = @(z) real(ifft2((fft2(z)).*h));
        Bt = @(z) real(ifft2((fft2(z)).*ht));
        
        dr = randperm(N)';
        mesnum =round(N*(1-missrate));
        OM = dr(1:mesnum);
        OMind = zeros(n);
        OMind(OM) = 1;
        R = @(z) z.*OMind;
        Rt = @(z) z.*OMind;
        
        Phi = @(z) R(B(z));
        Phit = @(z) Bt(Rt(z));
        u_obsv = Phi(u_org) + R(sigma*randn(n));
        epsilon = sqrt(mesnum*sigma^2);
end

% variables
c = u_obsv;
t = cell(K,1);
for j = 1:K
    t{j} = zeros(n);
end
s = cell(2*K+4, 1);
for j = 1:K
    s{j} = P(S{j}(t{j}));
    s{j+K} = t{j};
end
s{2*K+1} = D(u_obsv);
s{2*K+2} = u_obsv;
s{2*K+3} = Phi(u_obsv);
d = s;

% for the 1st step
Delta = @(z) Dt(D(z));
Gamma = @(z) z + Phit(Phi(z));
I_K2Gamma = @(z) z + (K/2)*Gamma(z);
A = @(z) reshape(Gamma(reshape(z,n)) + I_K2Gamma(Delta(reshape(z,n))),N,1);

% for the 2nd step
prox = cell(2*K+3,1);
for j = 1:K
    prox{j} = @(z) ProxPreBNN(z, gamma, blocksize);
    prox{j+K} = @(z) ProjAverageConst(z, 0);
end
prox{2*K+1} = @(z) ProxTVnorm(z, lambda*gamma);
prox{2*K+2} = @(z) ProjDynamicRangeConstraint(z, [0,1]);
prox{2*K+3} = @(z) ProjL2ball(z, u_obsv, epsilon);

%% main loop

stopcri = 1e-1; % stopping criterion
maxiter = 2000; % maximum number of iteration
disp('ADMM is running...')
for i = 1:maxiter
    
    v = cellfun(@(z1, z2) z1 - z2, s, d, 'UniformOutput', false);
    temp = zeros(n);
    for j = 1:K
        temp = temp + St{j}(Pt(v{j})) + v{j+K};
    end
    rhs = I_K2Gamma(Dt(v{2*K+1})) + v{2*K+2} + Phit(v{2*K+3})-0.5*Gamma(temp);

    cpre = c;
    [c, flag] = pcg(A, rhs(:));
    c = reshape(c, n);
    error = sqrt(sum(sum(sum((c-cpre).^2))));
    tpre = t;
    sumt = zeros(n);
    for j = 1:K
        t{j} = 0.5*(Delta(c) - Dt(v{2*K+1}) + St{j}(Pt(v{j})) + v{j+K});
        error = error + sqrt(sum(sum(sum((t{j}-tpre{j}).^2))));
        s{j} = prox{j}(P(S{j}(t{j}))+d{j});
        s{j+K} = prox{j+K}(t{j}+d{j+K});
        sumt = sumt + t{j}; 
    end
    u = c + sumt;
    if error < stopcri
        break;
    end
    s{2*K+1} = prox{2*K+1}(D(c)+d{2*K+1});
    s{2*K+2} = prox{2*K+2}(u+d{2*K+2});
    s{2*K+3} = prox{2*K+3}(Phi(u)+d{2*K+3});
    
    for j = 1:K
        d{j} = d{j} + P(S{j}(t{j})) - s{j};
        d{j+K} = d{j+K} + t{j} - s{j+K};
    end
    d{2*K+1} = d{2*K+1} + D(c) - s{2*K+1};
    d{2*K+2} = d{2*K+2} + u - s{2*K+2};
    d{2*K+3} = d{2*K+3} + Phi(u) - s{2*K+3};
    
end

%% result plot

psnr = EvalImgQuality(u, u_org, 'PSNR');
ssim = EvalImgQuality(u, u_org, 'SSIM');
disp(['PSNR = ', num2str(psnr,4), ', SSIM = ', num2str(ssim,4)]);
    
plotsize = [2, 3];
ImgPlot(u_obsv, 'Given image', 1, [plotsize,1]);
ImgPlot(c, 'Cartoon', 1, [plotsize,2]);
ImgPlot(sumt+0.5, 'Texture', 1, [plotsize,3]);
ImgPlot(u, 'Result image', 1, [plotsize,4]);
ImgPlot(u_org, 'Original image', 1, [plotsize,5]);
plotsize = [ceil(K/3), 3];
for j = 1:K
    ImgPlot(t{j}+0.5, ['Subtexture', num2str(j)], 2, [plotsize,j]);
end

