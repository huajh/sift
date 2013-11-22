
% load image and normalized the pixel value in the range [0,1]
bw1 = imreadbw('books/box.png');
bw2 = imreadbw('books/box_in_scene.png');

bw1 = bw1-min(bw1(:));
bw1 = bw1/max(bw1(:));
bw2 = bw2-min(bw2(:));
bw2 = bw2/max(bw2(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sift
[desc1,keyp1,DoGs1,Params1] = sift(bw1);
[desc2,keyp2,DoGs2,Params2] = sift(bw2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot scale space
figure(1); clf; PlotScaleSpace(DoGs1,Params1); colormap gray;
%figure(2); clf; PlotScaleSpace(DoGs2,Params2); colormap gray;
drawnow;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot key points
figure(3); clf;
subplot(1,2,1) ; imagesc(bw1) ; colormap gray ;
hold on ;
h = plotKeyPoints( keyp1,'style','arrow' ) ; set(h,'LineWidth',2,'Color','g') ;
h = plotKeyPoints( keyp1,'style','arrow' ) ; set(h,'LineWidth',1,'Color','k') ;

subplot(1,2,2) ; imagesc(bw2) ; colormap gray ;
hold on ;
h=plotKeyPoints( keyp2,'style','arrow') ; set(h,'LineWidth',2,'Color','g') ;
%h=plotKeyPoints( keyp2,'style','arrow') ; set(h,'LineWidth',1,'Color','k') ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Feature/keyPoints matching.

desc1 =uint8(512*desc1);
desc2 =uint8(512*desc2);

match_pairs = siftmatch(desc1,desc2);

figure(4); clf;
plotmatchPairs(bw1,bw2,keyp1(1:2,:),keyp2(1:2,:),match_pairs);
drawnow;





