function [ smoothed ] = imGauss2Dsmooth( bw,sigma)
%GUASSSMOOTH Summary of this function goes here
%   bw: gray image range [0,1]
%   sigma: the scale

width=3*sigma;
support = -width:0.1:width;
gauss2D = exp(-(support/sigma).^2/2);
gauss2D = gauss2D / sum(gauss2D);
smoothed = conv2(bw,gauss2D,'same');
end

