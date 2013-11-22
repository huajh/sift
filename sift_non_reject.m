function [Descriptors,KeyPoints,DoGs,Params] = sift_non_reject(img)
%% Reference:   
%   Lowe, David G. "Distinctive image features from scale-invariant keypoints." 
%   International journal of computer vision 60.2 (2004): 91-110.
%
%   author:         Junhao Hua
%                   Department of Information Science & Electronic Engineering(ISEE),
%                   ZheJiang University
%
%   email:          huajh7 AT gmail DOT com
%   Last update:    2013/10/22
%%

%%convert rgb to gray
if ndims(img) == 3
    img = rgb2gray(im2double(img));
end
% normalized
img = img-min(img(:));
img = img/max(img(:));

%% Calculate the Difference of Gaussian
sigma0 = 1.6;
S = 3;
DoGs = DiffofGauss(img,S,sigma0);
Params.sigma0 = sigma0;
Params.S = S;
OCTAVE_NUM = length(DoGs);
Params.O = OCTAVE_NUM;

threshold = 0.04 / S / 2 ; % 
MAGNIF_SIZE = 3;
SUBREGIONS_BINS = 4;
ORIENTS_BINS = 8;
KeyPoints = [];
Descriptors = [];

for i = 1:OCTAVE_NUM            
    %% Local Extrema Detection
    idx = findlocalmax(DoGs{i},0.8*threshold); % discard week points.
    idx = [idx,findlocalmax(-DoGs{i},0.8*threshold)];
    % convert a index into an array     
    [M,N,O]=size(DoGs{i}); % S == O
    [Row_idx,Col_idx,sca_idx] = ind2sub([M,N,O],idx); 
    
    % resharp index range [1,M-1],[1,M-1],[1,S-2]
    y = Row_idx - 1;
    x = Col_idx - 1;
    sca_idx = sca_idx - 1;
    LocalMax_OCT = [x(:)';y(:)';sca_idx(:)'];
    
    %% Orientation assignment
    KeyPoints_OCT = CalcKeypointsOrient(LocalMax_OCT,DoGs{i},S,sigma0);
    %KeyPoints_OCT = siftormx(KeyPoints_OCT,DoGs{i},S,-1,sigma0);
    
    %% Add to KeyPoints Vector    
    x     = 2^(i-2) * KeyPoints_OCT(1,:); %zoom to origin size
    y     = 2^(i-2) * KeyPoints_OCT(2,:);
    sigma = 2^(i-2) * sigma0 * 2.^(KeyPoints_OCT(3,:)/S);
    KeyPoints = [KeyPoints,[x(:)';y(:)';sigma(:)'; KeyPoints_OCT(4,:)] ];            
    
    %% Local Image descriptor     
    descrp = CalcDescriptor(KeyPoints_OCT,DoGs{i},S,sigma0,MAGNIF_SIZE,SUBREGIONS_BINS,ORIENTS_BINS);    
    Descriptors = [Descriptors, descrp];
    
    
end
