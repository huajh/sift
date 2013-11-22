
% load image and normalized the pixel value in the range [0,1]
bw1 = imreadbw('books/box.png');
bw1 = bw1-min(bw1(:));
bw1 = bw1/max(bw1(:));
figure(1);
imagesc(bw1) ; colormap gray ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot key points
% sift

[desc1,keyp1,DoGs1,Params1] = sift(bw1);
figure(3); clf;
imagesc(bw1) ; colormap gray ;
hold on ;
h = plotKeyPoints( keyp1,'style','arrow' ) ; set(h,'LineWidth',1.5,'Color','g') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5) ; clf ;
K=9 ;
for k=1:K
  tightsubplot(K,k) ;
  hold on ;
  h=plotsiftdescriptor(d(:,matches(1,k))) ;
 % lh=plotsiftdescriptor(ld(:,matches(2,k))) ;
  set(h,'LineWidth',2) ;
  set(lh,'Color','r') ;
  axis tight;axis square;axis off;
end