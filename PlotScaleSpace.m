function PlotScaleSpace(DoGs,Params)

if nargin > 2
    error('Too many arguments.') ;
end
sigma0 = Params.sigma0;
S = Params.S;
nlevels = Params.S+2;
OCTAVE_NUM = Params.O;
for oi=1:OCTAVE_NUM
  for si=1:nlevels
    tightsubplot(nlevels,OCTAVE_NUM, nlevels*(oi-1)+si) ;
    s = si-1;
    o = oi-1;
    sigma = sigma0 * 2^(s/S + o-1) ;
    F=squeeze(DoGs{oi}(:,:,si)) ;
    [M,N]=size(F) ;
    imagesc(squeeze(DoGs{oi}(:,:,si))) ; axis image ; axis off ;
    h=text(M/10,N/20,sprintf('(o,s)=(%d,%d), sigma=%.2f',o,s,sigma)) ;
    set(h,'BackgroundColor','w','Color','k') ;
  end
end


