function DoGs = DiffofGauss(img,S_0,sigma_0)
% parameters setting 
    S = 3;          %scale resolution; the number of scale samples. 
                    %It totally needs S+3 blurred images and S+2 Difference 
                    %of Gaussian(DOG) per octave.
    if (nargin >=2)
        S = S_0;
    end    
    sigma0 = 1.6;
    if (nargin >=3)
        sigma0 = sigma_0;   %the scale of keypoint
    end
    sigma_pre = 0.5;    
    k = 2^(1/S);    %constant multiplicative factor
    % The image resolution will decrease 1/4 with incrementation of octave.
    % the minimum resolution 4x4 is required.    
    OCTAVE_NUM = floor(log2(min(size(img))))-2;    
    SCALE_NUM = S+3;
    DOGS_NUM = S+2;

%pre-smoothing    
% assume the origin image has a blur of at least \sigma = 0.5
% double the size of the input image using blinear interpolation
% then the doubled image has \sigma = 1.0 relative to its new spacing.
   % img = imresize(img,2,'bilinear');     
   img = doubleSize(img);
%
% g ** f(sqrt(t1^2+t2^2)) = g ** f(t1) ** f(t2)  % ** means convolution
for i=1:OCTAVE_NUM
    sigma_Iter = (1/k)*sigma0*sqrt(k^2-1);
    for j=1:SCALE_NUM
        % the additional/extra smoothing
        % sigma_cur^2 - sigma_pre^2
        if i==1 && j==1
            sigma_head = sqrt(sigma0^2-(sigma_pre*2)^2);
            [M,N,~] = size(img);
            Octave{1} = zeros(M,N,SCALE_NUM);
            Octave{1}(:,:,1) = imGauss2Dsmooth(img,sigma_head);
        elseif j==1            
            %Each octave has S+3 blurred images, S+2 DoGs and S scales, 
            % which is indexed with 2,3,...,S,S+1. For the continuity of
            % scale space, the first images of the next octave has the same 
            % scale with the (S+1)th images of previous stacks/octave.
            %
            %down-sampled by a factor a 1/2, then it doubled the \sigma
            %automatically. 
            [M,N,~]=size(Octave{i-1});
            Octave{i} = zeros(floor(M/2+0.5),floor(N/2+0.5),SCALE_NUM);
            Octave{i}(:,:,1) = Octave{i-1}(1:2:end,1:2:end,S+1);
        else
            sigma_Iter = k*sigma_Iter;
            Octave{i}(:,:,j) = imGauss2Dsmooth(Octave{i}(:,:,j-1),sigma_Iter);
        end        
    end
end

for i=1:OCTAVE_NUM
    [M,N,~] = size(Octave{i});
    DoGs{i} = zeros(M,N,DOGS_NUM);
    for j=1:DOGS_NUM
        DoGs{i}(:,:,j) =  Octave{i}(:,:,j+1)-Octave{i}(:,:,j);
    end
end

function J = doubleSize(I)
[M,N]=size(I) ;
J = zeros(2*M,2*N) ;
J(1:2:end,1:2:end) = I ;
J(2:2:end-1,2:2:end-1) = ...
  0.25*I(1:end-1,1:end-1) + ...
  0.25*I(2:end,1:end-1) + ...
  0.25*I(1:end-1,2:end) + ...
  0.25*I(2:end,2:end) ;
J(2:2:end-1,1:2:end) = ...
  0.5*I(1:end-1,:) + ...
    0.5*I(2:end,:) ;
J(1:2:end,2:2:end-1) = ...
  0.5*I(:,1:end-1) + ...
    0.5*I(:,2:end) ;