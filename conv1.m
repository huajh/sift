im = imread('data/img3.jpg');
figure(1);
image(im);

% linear filter / convolution
bw = double(rgb2gray(im));

figure(2);
image(bw);
gradkernel = [-1 1];
dx = abs(conv2(bw,gradkernel,'same'));
image(dx);
colorbar; colormap gray
[dx,dy] = gradient(bw);
gradmag = sqrt(dx.^2 + dy.^2);
image(gradmag);

% smoothing with gaussian filter
%figure(3);
sigma=3;
width=3*sigma;
support = -width:width;
gauss2D = exp(-(support/sigma).^2/2);
gauss2D = gauss2D / sum(gauss2D);
smoothed = conv2(conv2(bw,gauss2D,'same'),gauss2D','same');
image(smoothed);
colormap(gray(255));
gauss3D = gauss2D'*gauss2D;
smoothed = conv2(bw,gauss3D,'same');
%image(smoothed);

% edge detection with smoothed images
%figure(4);

[dx,dy] = gradient(smoothed);
gradmag = sqrt(dx.^2+dy.^2);
gmax = max(max(gradmag));
imshow(gradmag);
colormap(gray(gmax));

% edge normal
figure(5);
hold on;
image(smoothed);
colormap(gray(255));
[m,n] = size(gradmag);
edges = (gradmag > 0.3*gmax);
indx = find(edges);
[posx,posy] = meshgrid(1:n,1:m);
posx2 = posx(indx);posy2 = posy(indx);
gm2 = gradmag(indx);
sintheta = dx(indx)./gm2;
costheta = -dy(indx)./gm2;

quiver(posx2,posy2,gm2.*sintheta/10, -gm2.*costheta/10, 0);
hold off;



























