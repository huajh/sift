nexp = 2 ;
switch nexp
  case 1
    I1=imreadbw('data/landscape-a.jpg') ; % I1=I1(1:2:end,:) ;
    I2=imreadbw('data/landscape-b.jpg') ; % I2=I2(1:2:end,:) ;
    I1c=double(imread('data/landscape-a.jpg'))/255.0 ;
    I2c=double(imread('data/landscape-b.jpg'))/255.0 ;

  case 2
    I1=imreadbw('data/vessel-1.pgm') ; I1c = I1 ;
    I2=imreadbw('data/vessel-2.pgm') ; I2c = I2 ;
end

I1=imsmooth(I1,.1) ;
I2=imsmooth(I2,.1) ;

I1=I1-min(I1(:)) ;
I1=I1/max(I1(:)) ;
I2=I2-min(I2(:)) ;
I2=I2/max(I2(:)) ;

S=3 ;

fprintf('Computing frames and descriptors.\n') ;
[frames1,descr1,gss1,dogss1] = sift( I1, 'Verbosity', 1, 'Threshold', ...
                                     0.005, 'NumLevels', S ) ;
[frames2,descr2,gss2,dogss2] = sift( I2, 'Verbosity', 1, 'Threshold', ...
                                     0.005, 'NumLevels', S ) ;

figure(11) ; clf ; plotss(dogss1) ; colormap gray ;
figure(12) ; clf ; plotss(dogss2) ; colormap gray ;
drawnow ;

figure(2) ; clf ;
tightsubplot(1,2,1) ; imagesc(I1) ; colormap gray ; axis image ;
hold on ;
h=plotsiftframe( frames1 ) ; set(h,'LineWidth',2,'Color','g') ;
h=plotsiftframe( frames1 ) ; set(h,'LineWidth',1,'Color','k') ;

tightsubplot(1,2,2) ; imagesc(I2) ; colormap gray ; axis image ;
hold on ;
h=plotsiftframe( frames2 ) ; set(h,'LineWidth',2,'Color','g') ;
h=plotsiftframe( frames2 ) ; set(h,'LineWidth',1,'Color','k') ;

fprintf('Computing matches.\n') ;
% By passing to integers we greatly enhance the matching speed (we use
% the scale factor 512 as Lowe's, but it could be greater without
% overflow)
descr1=uint8(512*descr1) ;
descr2=uint8(512*descr2) ;
tic ;
matches=siftmatch( descr1, descr2 ) ;
fprintf('Matched in %.3f s\n', toc) ;

figure(3) ; clf ;
plotmatches(I1c,I2c,frames1(1:2,:),frames2(1:2,:),matches,...
  'Stacking','v') ;
drawnow ;

% Movie
figure(4) ; set(gcf,'Position',[10 10 1024 512]) ;
figure(4) ; clf ;
tightsubplot(1,1);
imagesc(I1) ; colormap gray ; axis image ; hold on ;
h=plotsiftframe( frames1 ) ; set(h,'LineWidth',1,'Color','g') ;
h=plot(frames1(1,:),frames1(2,:),'r.') ;
MOV(1)=getframe ;

figure(4) ; clf ;
tightsubplot(1,1);
imagesc(I2) ; colormap gray ; axis image ; hold on ;
h=plotsiftframe( frames2 ) ; set(h,'LineWidth',1,'Color','g') ;
h=plot(frames2(1,:),frames2(2,:),'r.') ;
MOV(2)=getframe ;
