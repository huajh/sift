
% load image and normalized the pixel value in the range [0,1]
bw1 = imreadbw('books/box.png');

bw1 = bw1-min(bw1(:));
bw1 = bw1/max(bw1(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sift
[desc1,keyp1,DoGs1,Params1] = sift(bw1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot scale space
figure(1); clf; PlotScaleSpace(DoGs1,Params1); colormap gray;
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot key points
figure(3); clf;
imagesc(bw1) ; colormap gray ;
hold on ;
h = plotKeyPoints( keyp1,'style','arrow' ) ; set(h,'LineWidth',1.5,'Color','g') ;
%h = plotKeyPoints( keyp1,'style','arrow' ) ; set(h,'LineWidth',1,'Color','k') ;




