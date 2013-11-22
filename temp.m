%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot key points
figure(3); 

imagesc(bw2) ; colormap gray ;
hold on ;
h=plotKeyPoints( keyp2,'style','arrow') ; set(h,'LineWidth',1.5,'Color','g') ;
%h=plotKeyPoints( keyp2,'style','arrow') ; set(h,'LineWidth',1,'Color','k') ;
