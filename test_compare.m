% load image and normalized the pixel value in the range [0,1]
bw1 = imreadbw('plane/plane.jpg');
bw1 = bw1-min(bw1(:));
bw1 = bw1/max(bw1(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sift

[~,keyp1,~,~] = sift_non_reject(bw1);
[~,keyp2,~,~] = sift_non_elim_edge_resp(bw1);
[~,keyp3,~,~] = sift(bw1);
key_num1 = size(keyp1,2)
key_num2 = size(keyp2,2)
key_num3 = size(keyp3,2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot key points
figure(1); clf;
subplot(2,2,1) ; imagesc(bw1) ; colormap gray ; xlabel('(a)');
hold on ;
subplot(2,2,2) ; imagesc(bw1) ; colormap gray ; xlabel('(b)');
h = plotKeyPoints( keyp2,'style','arrow' ) ; set(h,'LineWidth',1.5,'Color','g') ;
subplot(2,2,3) ; imagesc(bw1) ; colormap gray ; xlabel('(c)');
h = plotKeyPoints( keyp3,'style','arrow' ) ; set(h,'LineWidth',1.5,'Color','g') ;
subplot(2,2,4) ; imagesc(bw1) ; colormap gray ;xlabel('(d)');
h = plotKeyPoints( keyp1,'style','arrow' ) ; set(h,'LineWidth',1.5,'Color','g') ;