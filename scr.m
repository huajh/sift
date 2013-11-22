match_pairs = siftmatch(desc1,desc2);

figure(4); clf;
plotmatchPairs(bw1,bw2,keyp1(1:2,:),keyp2(1:2,:),match_pairs);
drawnow;
