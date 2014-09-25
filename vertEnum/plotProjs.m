clear;close all;

dir='/Users/jean/Desktop/Publications/PaperPolytopeRec/Figs/';

poly.d=3;
poly.m=20;
while 1
    
    poly.A=randn(poly.m,poly.d);
    poly.b=rand(poly.m,1);
    poly.P=Polyhedron(poly.A,poly.b);

    if (~poly.P.isEmptySet()) && (poly.P.isFullDim())
        break;
    end
end



% Projections
projs.Pxy=poly.P.projection([1 2],'mplp');
projs.Pxz=poly.P.projection([1 3],'mplp');
projs.Pyz=poly.P.projection([2 3],'mplp');

% Vertices of projections
projs.Vxy=projs.Pxy.V;
projs.Vxz=projs.Pxz.V;
projs.Vyz=projs.Pyz.V;

% Compute convex hull
projs.chxy=convhull(projs.Vxy(:,1),projs.Vxy(:,2));
projs.chxz=convhull(projs.Vxz(:,1),projs.Vxz(:,2));
projs.chyz=convhull(projs.Vyz(:,1),projs.Vyz(:,2));

col_plane=[.9 .9 .9];
col_proj=[.3 .3 .3];
ft_size=50;

ox=10;dx=0.8;
oy=10;dy=0.8;
oz=10;dz=0.6;

% XY projection
fig=figure('Color','w');
patch([ox-dx,ox-dx,ox+dx,ox+dx],[oy-dy,oy+dy,oy+dy,oy-dy],(oz-dz)*ones(1,4),col_plane);hold on;
patch(projs.Vxy(projs.chxy,1)+ox,projs.Vxy(projs.chxy,2)+oy,(oz-dz)*ones(1,length(projs.chxy)),col_proj);hold on;
text(ox+.8*dx,oy+.8*dy,oz-dz,'$\mathcal{H}_1$','Interpreter','LaTex','fontSize',ft_size,'fontName','Times New Roman');
% YZ projection
patch((ox-dx)*ones(1,4),[oy-dy,oy-dy,oy+dy,oy+dy],[oz-dz,oz+dz,oz+dz,oz-dz],col_plane);hold on;
patch((ox-dx)*ones(1,length(projs.chyz)),projs.Vyz(projs.chyz,1)+oy,projs.Vyz(projs.chyz,2)+oz,col_proj);hold on;
text(ox-dx,oy+.8*dy,oz+.8*dz,'$\mathcal{H}_2$','Interpreter','LaTex','fontSize',ft_size,'fontName','Times New Roman');
% XZ projection
patch([ox-dx,ox-dx,ox+dx,ox+dx],(oy-dy)*ones(1,4),[oz-dz,oz+dz,oz+dz,oz-dz],col_plane);hold on;
patch(projs.Vxz(projs.chxz,1)+ox,(oy-dy)*ones(1,length(projs.chxz)),projs.Vxz(projs.chxz,2)+oz,col_proj);hold on;
text(ox+.8*dx,oy-dy,oz+.8*dz,'$\mathcal{H}_3$','Interpreter','LaTex','fontSize',ft_size,'fontName','Times New Roman');
% 3-polytope
plot(Polyhedron([poly.P.V(:,1)+ox,poly.P.V(:,2)+oy,poly.P.V(:,3)+oz]),...
                 struct('shade',0.5,'edgecolor',[.2 .2 .2],'linestyle','--','linewidth',3,'color','w','DisplayName','Fitted polytope'));
text(ox+0.3,oy+0.3,oz-0.2,'$P(A\tilde{R}^{\mathrm{T}},B\tilde{t})$','Interpreter','LaTex','fontSize',ft_size,'fontName','Times New Roman');
axis equal;axis off;grid off;
set(gcf, 'renderer','opengl');
fprintf('Arrange before printing to pdf...\n');
pause;
myaa;
save2pdf([dir,'polyViews.pdf'],gcf,600);









%patch([0 30 30 0],ones(1,4).*0.01,[-o1 -o1 o1 o1],col_plane);
%patch(ones(4,1).*0.01, [0 30 30 0], [-o1 -o1 o1 o1],col_plane);
%patch(projs.Vxy(projs.chxy,2)+o1, ones(length(projs.chxy),1).*0.02,projs.Vxy(projs.chxy,1),col_proj);
%patch(ones(length(projs.chxz),1).*0.02,projs.Vxz(projs.chxz,2)+o1,projs.Vxz(projs.chxz,1),col_proj);
%xlabel('-X-');
%ylabel('-Y-');
%zlabel('-Z-');
%plot3(V1est(chv1est,2)+o1,zeros(length(chv1est),1),V1est(chv1est,1),'color','r','LineWidth',3,'LineStyle','--')
%plot3(zeros(length(chv2est),1),V2est(chv2est,2)+o2,V2est(chv2est,1),'color','r','LineWidth',3,'LineStyle','--')



