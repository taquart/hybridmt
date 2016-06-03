function [u,v] = drawhudsonnet(varargin)
%DRAWHUDSONNET Draw Hudson network and/or provide MT projection to [u,v] space.
%   Use DRAWHUDSONNET allows either to generate a Hudson plot figure following 
%   the Hudson, J.A., R.G. Pearce, and R.M.Rogers (1989), "Source type plot 
%   for inversion of the moment tensor", J. Geophys. Res., 94, 765–774. or to
%   project (and eventually plot) moment tensor data to [u,v] space of Hudson
%   network.
%
%   Syntax
%
%   DRAWHUDSONNET with no arguments just produce a Hudson plot.
%
%   [U,V] = DRAWHUDSONNET(M) produce a Hudson plot and project moment tensor 
%   data into [U,V] space. Vectors U and V can be than used to plot the points
%   in the generated figure, e.g. by calling: plot(U,V,'ok'); Moment tensor
%   matrix M can be:
%
%   * 3-by-3 matrix containing full moment tensor in a form
%     M = [m11 m12 m13; m21 m22 m23; m31 m32 m33] (note the moment tensor must
%     be symmetric, i.e. m12=m21, m32=m23 and m31=m13). 
%   * n-by-6 matrix where each row corresponds to a different moment tensor 
%     Each row contains 6 independent moment tensors components in 
%     the following format: M(i,:) = [m11 m12 m13 m22 m23 m33]; 
%   * n-by-9 matrix where each row corresponds to a different moment tensor. 
%     Each row contains 9 moment tensor components in the following format:
%     M(i,:) = [m11 m12 m13 m21 m22 m23 m31 m32 m33]; In this case columns
%     [4 7 8] are ignored and matrix is reduced to n-by-6 matrix.
%
%   DRAWHUDSONNET(M) with no output arguments handling produces a Hudson plot 
%   and immediately projects moment tensor data onto Hudson plot.
%
%   EXAMPLES
%
%   drawhudsonnet(eye(3)) generates Hudson plot and projects single moment
%   tensor of the form [1 0 0; 0 1 0; 0 0 0] onto it.
%
%   drawhudsonnet([1 0 0 1 0 1; 0 0 1 0 0 0; 1 0 0 1 0 -2]) generates Hudson
%   plot in projects three moment tensor onto it. The three moment tensors 
%   correspond to pure positive Isotropic, double-couple and negative CLVD
%   moment tensor.

%   Copyright 2015-2016 Grzegorz Kwiatek <kwiatek@gfz-potsdam.de>
%   $Revision: 1.0.2 $  $Date: 2016.06.03 $


if nargin == 1
  M = varargin{1};
if ismatrix(M) && size(M,1) == 3 && size(M,2) == 3
  [k,t] = ae_srctype(M);
  [u,v] = ae_uvt(t,k);
elseif ismatrix(M) && (size(M,2) == 6 || size(M,2) == 9)
  if size(M,2) == 9
    M = M(:,[1 2 3 5 6 9]);
  end
  u = nan(size(M,1),1);
  v = nan(size(M,1),1);
  
  for i=1:size(M,1)
    M_MATRIX = [ M(i,1) M(i,2) M(i,3); ...
      M(i,2) M(i,4) M(i,5); ...
      M(i,3) M(i,5) M(i,6)];
    [k,t] = ae_srctype(M_MATRIX);
    [u(i),v(i)] = ae_uvt(t,k);
  end
else
  error('Inappropriate matrix size. Either 3-by-3, n-by-6 or n-by-9 matrices are allowed.');
end
else
  error('Only one parameter is possible.');
end

set(gcf,'Color','w');
hold on;
fill([-1 0 0],[0 0 1],[0.95 0.95 0.95],'EdgeColor','none');
fill([0 1 0],[0 0 -1],[0.95 0.95 0.95],'EdgeColor','none');
hold off;
axis equal;
set(gca,'XLim',[-4/3 4/3]);
set(gca,'YLim',[-1 1]);
line([0 4/3 0 -4/3 0]',[1 1/3 -1 -1/3 1]','Color','k');
line([-1 1]',[0 0]','Color','k');
line([0 0]',[-1 1]','Color','k');
fs = 10; ha = 'HorizontalAlignment'; va = 'VerticalAlignment';
VV = (-1:0.1:1)';
[x,y] = ae_uvt(VV,0.5*ones(size(VV))); line(x,y,'Color','k','LineStyle','--');
[x,y] = ae_uvt(VV,-0.5*ones(size(VV))); line(x,y,'Color','k','LineStyle','--');
[x,y] = ae_uvt(0.5*ones(size(VV)),VV); line(x,y,'Color','k','LineStyle','--');
[x,y] = ae_uvt(-0.5*ones(size(VV)),VV); line(x,y,'Color','k','LineStyle','--');

offset = 0.02;
[x,y] = ae_uvt(1,1); text(x,y+offset,'Explosion','FontSize',fs,ha,'Center',va,'Bottom');
[x,y] = ae_uvt(1,-1); text(x,y-offset,'Implosion','FontSize',fs,ha,'Center',va,'Top');
[x,y] = ae_uvt(1,0); text(x,y-offset,'CLVD (-)','FontSize',fs,ha,'Center',va,'Top');
[x,y] = ae_uvt(1,-5/9); text(x,y-offset,'Anticrack','FontSize',fs,ha,'Left',va,'Top');
[x,y] = ae_uvt(-1,0); text(x,y-offset,'CLVD (+)','FontSize',fs,ha,'Center',va,'Top');
[x,y] = ae_uvt(-1,5/9); text(x,y+offset,'Tensile Crack','FontSize',fs,ha,'Right',va,'Bottom');
[x,y] = ae_uvt(0,0); text(x,y+offset,'DC','FontSize',fs,ha,'Center',va,'Bottom');
[x,y] = ae_uvt(-1,1/3); text(x,y,'LVD (+)','FontSize',fs,ha,'Right',va,'Bottom');
[x,y] = ae_uvt(1,-1/3); text(x,y,'LVD (-)','FontSize',fs,ha,'Left',va,'Top');

[x,y] = ae_uvt([1 1 1 1 -1 -1 0 -1 1],[1 -1 0 -5/9 0 5/9 0 1/3 -1/3]);
hold on;
plot(x,y,'ok','MarkerSize',4,'MarkerFaceColor','k');
hold off;

grid off;
set(gca,'Visible','off');

if nargout == 0 && ~isempty(u)
hold on;
  plot(u,v,'ok','MarkerFaceColor','r');
hold off;
end

%-------------------------------------------------------------------------
function [k,T] = ae_srctype(M)

% Hudson and Bowers transformation, 1999
[dummy,DD] = eig(M); %#ok<ASGLU>
EGV = DD([1 5 9]);
EGV = sort(EGV,'descend');
MM = mean(EGV);
EGV = EGV - MM; % EGV - deviatoric eigenvalues.
k = MM / (abs(MM) + max(abs(EGV(1)),abs(EGV(3))));

if (abs(MM) + max(abs(EGV(1)),abs(EGV(3)))) == 0
  error('Tensor cannot be translated into [k,T] space.');
end

if max(abs(EGV(1)),abs(EGV(3))) == 0
  T = 0;
else
  T = 2 * EGV(2) / max(abs(EGV(1)),abs(EGV(3)));
end

%-------------------------------------------------------------------------
function [UU,VV] = ae_uvt(TT,KK)

TAU = TT.*(1-abs(KK));
UU = NaN*ones(size(TT));
VV = UU;

% 2nd and 4th quadrants
II = (TAU > 0 & KK < 0) | (TAU < 0 & KK > 0);
UU(II) = TAU(II);
VV(II) = KK(II);

% First quadrant, Region A
II = (TAU < 4*KK) & (TAU >= 0 & KK >= 0);
UU(II) = TAU(II)./(1-TAU(II)/2);
VV(II) = KK(II)./(1-TAU(II)/2);

% First quadrant, Region B
II = (TAU >= 4*KK) & (TAU >= 0 & KK >= 0);
UU(II) = TAU(II)./(1-2*KK(II));
VV(II) = KK(II)./(1-2*KK(II));

% Third quadrant
II = (TAU >= 4*KK) & (TAU <= 0 & KK <= 0);
UU(II) = TAU(II)./(1+TAU(II)/2);
VV(II) = KK(II)./(1+TAU(II)/2);

II = (TAU < 4*KK) & (TAU <= 0 & KK <= 0);
UU(II) = TAU(II)./(1+2*KK(II));
VV(II) = KK(II)./(1+2*KK(II));

if sum(isnan(VV))
  error('Error plotting point in [uu,vv] space.');
end

