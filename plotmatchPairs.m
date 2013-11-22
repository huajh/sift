function h=plotmatchPairs(I1,I2,P1,P2,matches,varargin)
stack='h' ;
interactive=0 ;
only_interactive=0 ;

for k=1:2:length(varargin)
  switch lower(varargin{k})
   case 'stacking'
    stack=varargin{k+1} ;
   case 'interactive'
    interactive=varargin{k+1};
   otherwise
    error(['[Unknown option ''', varargin{k}, '''.']) ;
  end
end

% --------------------------------------------------------------------
%                                                           Do the job
% --------------------------------------------------------------------

[M1,N1,K1]=size(I1) ;
[M2,N2,K2]=size(I2) ;

switch stack
  case 'h'
    N3=N1+N2 ;
    M3=max(M1,M2) ;
    oj=N1 ;
    oi=0 ;
  case 'v'
    M3=M1+M2 ;
    N3=max(N1,N2) ;
    oj=0 ;
    oi=M1 ;
  case 'd'
    M3=M1+M2 ;
    N3=N1+N2 ;
    oj=N1 ;
    oi=M1 ;
  case 'o'
    M3=max(M1,M2) ;
    N3=max(N1,N2) ;
    oj=0;
    oi=0;
  otherwise
    error(['Unkown stacking type '''], stack, ['''.']) ;
end

% Combine the two images. In most cases just place one image next to
% the other. If the stacking is 'o', however, combine the two images
% linearly.
I=zeros(M3,N3,K1) ;
if stack ~= 'o'
  I(1:M1,1:N1,:) = I1 ;
  I(oi+(1:M2),oj+(1:N2),:) = I2 ;
else
  I(oi+(1:M2),oj+(1:N2),:) = I2 ;
  I(1:M1,1:N1,:) = I(1:M1,1:N1,:) + I1 ;
  I(1:min(M1,M2),1:min(N1,N2),:) = 0.5 * I(1:min(M1,M2),1:min(N1,N2),:) ;
end

axes('Position', [0 0 1 1]) ;
imagesc(I) ; colormap gray ; hold on ; axis image ; axis off ;

K = size(matches, 2) ;
nans = NaN * ones(1,K) ;

x = [ P1(1,matches(1,:)) ; P2(1,matches(2,:))+oj ; nans ] ;
y = [ P1(2,matches(1,:)) ; P2(2,matches(2,:))+oi ; nans ] ;

% if interactive > 1 we do not drive lines, but just points.
if(interactive > 1)
  h = plot(x(:),y(:),'ro') ;
else
  h = line(x(:)', y(:)') ;
end
set(h,'Marker','o','Color','g') ;
plot(x(:),y(:),'gs','LineWidth',1,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','green',...
                'MarkerSize',5)
% --------------------------------------------------------------------
%                                                          Interactive
% --------------------------------------------------------------------

if(~interactive), return ; end

sel1 = unique(matches(1,:)) ;
sel2 = unique(matches(2,:)) ;

K1 = length(sel1) ; %size(P1,2) ;
K2 = length(sel2) ; %size(P2,2) ;
X = [ P1(1,sel1) P2(1,sel2)+oj ;
      P1(2,sel1) P2(2,sel2)+oi ; ] ;

fig = gcf ;
is_hold = ishold ;
hold on ;

% save the handlers for later to restore
dhandler = get(fig,'WindowButtonDownFcn') ;
uhandler = get(fig,'WindowButtonUpFcn') ;
mhandler = get(fig,'WindowButtonMotionFcn') ;
khandler = get(fig,'KeyPressFcn') ;
pointer  = get(fig,'Pointer') ;

set(fig,'KeyPressFcn',        @key_handler) ;
set(fig,'WindowButtonDownFcn',@click_down_handler) ;
set(fig,'WindowButtonUpFcn',  @click_up_handler) ;
set(fig,'Pointer','crosshair') ;

data.exit        = 0 ;   % signal exit to the interactive mode
data.selected    = [] ;  % currently selected feature
data.X           = X ;   % feature anchors

highlighted = [] ;       % currently highlighted feature
hh = [] ;                % hook of the highlight plot

guidata(fig,data) ;
while ~ data.exit
  uiwait(fig) ;
  data = guidata(fig) ;
  if(any(size(highlighted) ~= size(data.selected)) || ...
     any(highlighted ~= data.selected) )

    highlighted = data.selected ;

    % delete previous highlight
    if( ~isempty(hh) )
      delete(hh) ;
    end

    hh=[] ;

    % each selected feature uses its own color
    c=1 ;
    colors=[1.0 0.0 0.0 ;
            0.0 1.0 0.0 ;
            0.0 0.0 1.0 ;
            1.0 1.0 0.0 ;
            0.0 1.0 1.0 ;
            1.0 0.0 1.0 ] ;

    % more than one feature might be seleted at one time...
    for this=highlighted

      % find matches
      if( this <= K1 )
      sel=find(matches(1,:)== sel1(this)) ;
      else
        sel=find(matches(2,:)== sel2(this-K1)) ;
      end
      K=length(sel) ;

      % plot matches
      x = [ P1(1,matches(1,sel)) ; P2(1,matches(2,sel))+oj ; nan*ones(1,K) ] ;
      y = [ P1(2,matches(1,sel)) ; P2(2,matches(2,sel))+oi ; nan*ones(1,K) ] ;

      hh = [hh line(x(:)', y(:)',...
                    'Marker','*',...
                    'Color',colors(c,:),...
                    'LineWidth',3)];

      if( size(P1,1) == 4 )
        f1 = unique(P1(:,matches(1,sel))','rows')' ;
        hp=plotKeyPoints(f1);
        set(hp,'Color',colors(c,:)) ;
        hh=[hh hp] ;
      end

      if( size(P2,1) == 4 )
        f2 = unique(P2(:,matches(2,sel))','rows')' ;
        f2(1,:)=f2(1,:)+oj ;
        f2(2,:)=f2(2,:)+oi ;
        hp=plotKeyPoints(f2);
        set(hp,'Color',colors(c,:)) ;
        hh=[hh hp] ;
      end

      c=c+1 ;
    end

    drawnow ;
  end
end

if( ~isempty(hh) )
  delete(hh) ;
end

if ~is_hold
  hold off ;
end

set(fig,'WindowButtonDownFcn',  dhandler) ;
set(fig,'WindowButtonUpFcn',    uhandler) ;
set(fig,'WindowButtonMotionFcn',mhandler) ;
set(fig,'KeyPressFcn',          khandler) ;
set(fig,'Pointer',              pointer ) ;

% ====================================================================
function data=selection_helper(data)
% --------------------------------------------------------------------
P = get(gca, 'CurrentPoint') ;
P = [P(1,1); P(1,2)] ;

d = (data.X(1,:) - P(1)).^2 + (data.X(2,:) - P(2)).^2 ;
dmin=min(d) ;
idx=find(d==dmin) ;

data.selected = idx ;

% ====================================================================
function click_down_handler(obj,event)
% --------------------------------------------------------------------
% select a feature and change motion handler for dragging

[obj,fig]=gcbo ;
data = guidata(fig) ;
data.mhandler = get(fig,'WindowButtonMotionFcn') ;
set(fig,'WindowButtonMotionFcn',@motion_handler) ;
data = selection_helper(data) ;
guidata(fig,data) ;
uiresume(obj) ;

% ====================================================================
function click_up_handler(obj,event)
% --------------------------------------------------------------------
% stop dragging

[obj,fig]=gcbo ;
data = guidata(fig) ;
set(fig,'WindowButtonMotionFcn',data.mhandler) ;
guidata(fig,data) ;
uiresume(obj) ;

% ====================================================================
function motion_handler(obj,event)
% --------------------------------------------------------------------
% select features while dragging

data = guidata(obj) ;
data = selection_helper(data);
guidata(obj,data) ;
uiresume(obj) ;

% ====================================================================
function key_handler(obj,event)
% --------------------------------------------------------------------
% use keypress to exit

data = guidata(gcbo) ;
data.exit = 1 ;
guidata(obj,data) ;
uiresume(gcbo) ;

