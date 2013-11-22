
surf(fspecial('gaussian',40,4))
surf(fspecial('gaussian',40,8))
surf(fspecial('gaussian',40,8) - fspecial('gaussian',40,4))


im=imread('3063.jpg');
bw = double(im(:,:,1))/256;
 for i =1:10
     gaussD = fspecial('gaussian',40,2*i)-fspecial('gaussian',40,i);
     mesh(gaussD) ; drawnow
     res = abs(conv2(bw,gaussD,'same'));
     res = res/max(max(res));
  %   imshow(res);
     title(['\bf i = ' num2str(i)]); drawnow     
 end