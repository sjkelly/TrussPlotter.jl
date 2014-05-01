module TrussPlotter

using PyPlot;


function plot_undeformed(xn,f,Idb,Ucomp,Rcomp,ien,nel,nen,nsd,ndf,nnp,axial)
    display_factor = 0.1;
    eps = 1e-4;

    # characteristic distances
    xmax=maximum(maximum(xn[1,:]));
    xmin=minimum(minimum(xn[1,:]));

    if (nsd >1)
        ymax=maximum(maximum(xn[2,:]));
        ymin=minimum(minimum(xn[2,:]));
        
        Lcar=maximum([xmax-xmin;ymax-ymin]);
    else
        Lcar=xmax-xmin;
    end;

    axis("equal")
    title("Undeformed mesh and BCs");
    plot_mesh_undeformed(nel,ien,xn,nnp,nsd);
    numbers(nel,ien,xn,nnp,nsd);
    #plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xn,nsd);
    #plot_bc_force(type,Lcar,f,display_factor,nnp,xn,nsd);
    #view(nsd);
    #hold off;
end

#function plot_results(type,xn,f,Idb,Ucomp,Rcomp,ien,nel,nen,nsd,ndf,nnp,axial)
#    display_factor=0.1;
#    eps=1e-4;

#    # characteristic distances
#    xmax=max(max(xn(1,:)));
#    xmin=min(min(xn(1,:)));

#    if (nsd >1)
#        ymax=max(max(xn(2,:)));
#        ymin=min(min(xn(2,:)));
#        
#        Lcar=max([xmax-xmin;ymax-ymin]);
#    else
#        Lcar=xmax-xmin;
#    end;

#    figure;
#    axis equal;
#    title('Undeformed and deformed mesh');
#    hold on;
#    set(gcf, 'Color', [1,1,1]); #Background color white
#    plot_mesh_undeformed(nel,ien,xn,nnp,nsd);
#    numbers(nel,ien,xn,nnp,nsd);
#    plot_mesh_deformed(type,xn,Ucomp,Idb,display_factor,Lcar,nel,ien,ndf,nsd,nen,axial);
#    view(nsd);
#    hold off;
#end

#function plot_results(type,xn,f,Idb,Ucomp,Rcomp,ien,nel,nen,nsd,ndf,nnp,axial)
#    display_factor=0.1;
#    eps=1e-4;

#    # characteristic distances
#    xmax=max(max(xn(1,:)));
#    xmin=min(min(xn(1,:)));

#    if (nsd >1)
#        ymax=max(max(xn(2,:)));
#        ymin=min(min(xn(2,:)));
#        
#        Lcar=max([xmax-xmin;ymax-ymin]);
#    else
#        Lcar=xmax-xmin;
#    end;

#    figure;
#    axis equal;
#    title('Reactions');
#    hold on;
#    set(gcf, 'Color', [1,1,1]); #Background color
#    plot_mesh_undeformed(nel,ien,xn,nnp,nsd);
#    plot_reactions(type,Lcar,Rcomp,Idb,display_factor,xn,nnp,ndf,nsd)
#    view(nsd);
#    hold off;
#end;


##############################################################################################
#                                       PLOT FUNCTIONS                                       #
##############################################################################################
############################
# undeformed configuration #
############################
function plot_mesh_undeformed(nel,ien,xn,nnp,nsd)
    for e=1:nel
        node1=ien[1,e];
        node2=ien[2,e];
        if (nsd == 3)
            
        end
        x0=[xn[1,node1];xn[1,node2]];
        if (nsd > 1)
            y0=[xn[2,node1];xn[2,node2]];
            if nsd ==3
                z0=[xn[3,node1];xn[3,node2]];
            end
        else
            y0=[0;0];
        end;
        if nsd < 3
            plot(x0,y0, color="Black", marker=".");
        else
            line(x0,y0,z0);
        end
    end;
end;



###################################################
# numbers the nodes and elements - truss and beam #
###################################################
function numbers(nel,ien,xn,nnp,nsd)
    fontsize = 12;
    # number the nodes in the undeformed configuration
    for n=1:nnp
        if (nsd == 3)
            text(xn[1,n],xn[2,n], xn[3,n], string(n),fontsize=fontsize, color="Blue");
        end
        if (nsd == 2)
            text(xn[1,n],xn[2,n], string(n),fontsize=fontsize, color="Blue");
        end
        if (nsd == 1)
            text(xn[1,n],0, string(n),fontsize=fontsize, color="Blue");
        end
    end
    # number the elements in the undeformed configuration
    for elt=1:nel
        node1=ien[1,elt];
        node2=ien[2,elt];
        #compute the coordinates of the middle point
        xg=0.5*(xn[:,node1]+xn[:,node2]);
        s=string("(",elt,")");
        if (nsd == 3)
            text(xg[1],xg[2], xg[3], s, fontsize=fontsize, color="Blue");
        end
        if (nsd == 2)
            text(xg[1],xg[2], s, fontsize=fontsize, color="Blue");
        end
        if (nsd == 1)
            text(xg[1],0, s,fontsize=fontsize, color="Blue");
        end;
    end
end

###########################
## deformed configuration #
###########################
#function plot_mesh_deformed(type,xn,Ucomp,Idb,display_factor,Lcar,nel,ien,ndf,nsd,nen,axial)

#switch type
#    case 'truss'
#        scale=display_factor*Lcar/max(max(abs(Ucomp))); # scale factor for the displacements
#        legend_comp_flag = 0;
#        legend_ten_flag = 0;
#        for e=1:nel
#            node1=ien(1,e);
#            node2=ien(2,e);
#            xt(:,1)=xn(:,node1);
#            xt(:,2)=xn(:,node2);
#            
#            if (nsd == 1)
#                xt(2,1)=0;
#                xt(2,2)=0;
#            end;
#            
#            for i=1:ndf
#                if Idb(i,node1)~=0
#                    xt(i,1)=xt(i,1)+scale*Ucomp(i,node1);
#                else
#                    xt(i,1)=xt(i,1)+scale*Ucomp(i,node1);
#                end;
#                if Idb(i,node2)~=0
#                    xt(i,2)=xt(i,2)+scale*Ucomp(i,node2);
#                else
#                    xt(i,2)=xt(i,2)+scale*Ucomp(i,node2);
#                end;
#            end;
#            if nsd < 3
#                plot(xt(1,:),xt(2,:)','r-o');
#            end
#            if nsd == 3
#                if (axial(2,e) > 0)
#                   h(1,1) = plot3(xt(1,:),xt(2,:),xt(3,:), 'r-o');
#                   if legend_ten_flag == 0
#                       legend_ten_flag = 1;
#                   end
#                else
#                   h(2,1) = plot3(xt(1,:),xt(2,:),xt(3,:), 'b-o');
#                   if legend_comp_flag == 0
#                       legend_comp_flag = 1;
#                   end
#                end
#            end
#        end
#        if legend_comp_flag && legend_ten_flag
#            legend(h,'Tension', 'Compression')
#        else
#            if legend_comp_flag
#                legend(h(2,1),'Compression')
#            else
#                if legend_ten_flag
#                    legend(h(1,1),'Tension')
#                end
#            end
#        end
#    case 'beam'
#        xi=-1:0.1:1;
#        nxi=size(xi,2);
#        
#        Umax=0;
#        #######################################################################################
#        #                                   SCALE FACTOR                                      #
#        #######################################################################################
#        for e=1:nel
#            node1=ien(1,e);
#            node2=ien(2,e);
#            if nsd == 2 
#            # displacement + rotation
#            Un(1,1)=Ucomp(1,node1);
#            Un(1,2)=Ucomp(2,node1);
#            Un(1,3)=Ucomp(3,node1);
#            Un(2,1)=Ucomp(1,node2);
#            Un(2,2)=Ucomp(2,node2);
#            Un(2,3)=Ucomp(3,node2);
#            elseif nsd ==3
#            Un(1,1)=Ucomp(1,node1);
#            Un(1,2)=Ucomp(2,node1);
#            Un(1,3)=Ucomp(3,node1);
#            Un(1,4)=Ucomp(1,node1);
#            Un(1,5)=Ucomp(2,node1);
#            Un(1,6)=Ucomp(3,node1);
#            Un(2,1)=Ucomp(1,node2);
#            Un(2,2)=Ucomp(2,node2);
#            Un(2,3)=Ucomp(3,node2);
#            Un(2,4)=Ucomp(1,node2);
#            Un(2,5)=Ucomp(2,node2);
#            Un(2,6)=Ucomp(3,node2);
#            end
#            
#            # intial postion
#            x0=zeros(nxi,nsd);
#            for n=1:nxi
#                x0(n,:)=((1-xi(n))/2)*xn(:,node1)'+((1+xi(n))/2)*xn(:,node2)';
#            end
#            
#            Le=sqrt((xn(1,node1)-xn(1,node2))^2+(xn(2,node1)-xn(2,node2))^2);
#            
#            # g(i,j) - vector i coordinate j
#            g=zeros(nsd,nsd);
#            g(1,:)=(1/Le)*(xn(:,node2)-xn(:,node1))';
#            
#            g(2,1)=-g(1,2);
#            g(2,2)=g(1,1);
#            
#            # to make g_3 = z - so that the rotation are the same in both systems of coordinates
#            cross_prod=g(1,1)*g(2,2)-g(1,2)*g(2,1);
#            if cross_prod < 0
#                g(2,:)=-g(2,:);
#            end;
#            
#            #transformation matrices
#            beta=g;
#            betap=inv(beta);
#            
#            # transform displacement coordinates in the global system
#            un=zeros(nen,nsd);
#            
#            for n=1:nen
#                for i=1:nsd
#                    for j=1:nsd
#                        un(n,i)=un(n,i)+Un(n,j)*betap(j,i);
#                    end;
#                end;
#            end;
#            
#            for n=1:nen
#                un(n,3)=Un(n,3);
#            end;
#            
#            # transverse displacement in the local system of coordinates
#            for i=1:nxi
#                N1(i)=(1/4)*(1-xi(i))^2*(2+xi(i));
#                N2(i)=(Le/8)*(1-xi(i)^2)*(1-xi(i));
#                N3(i)=(1/4)*(1+xi(i))^2*(2-xi(i));
#                N4(i)=(Le/8)*(-1+xi(i)^2)*(1+xi(i));
#                u(i,2)=N1(i)*un(1,2)+N2(i)*un(1,3)+N3(i)*un(2,2)+N4(i)*un(2,3);
#            end;
#            
#            # horizontal displacement in the local system of coordinates
#            for i=1:nxi
#                u(i,1)=(1/2)*((1-xi(i))*un(1,1)+(1+xi(i))*un(2,1));
#            end;
#            
#            # transform from local to global system
#            U=zeros(nxi,nsd);
#            for i=1:nsd
#                for n=1:nxi
#                    for j=1:nsd
#                        U(n,i)=U(n,i)+u(n,j)*beta(j,i);
#                    end;
#                end;
#            end;
#            
#            # scale factor
#            if (max(max(abs(U)))>Umax)
#                Umax=max(max(abs(U)));
#            end;
#        end;
#        
#        scale=display_factor*Lcar/Umax;
#        
#        #######################################################################################
#        #                                       PLOT                                          #
#        #######################################################################################
#        for e=1:nel
#            node1=ien(1,e);
#            node2=ien(2,e);
#            
#            # displacement + rotation
#            Un(1,1)=Ucomp(1,node1);
#            Un(1,2)=Ucomp(2,node1);
#            Un(1,3)=Ucomp(3,node1);
#            Un(2,1)=Ucomp(1,node2);
#            Un(2,2)=Ucomp(2,node2);
#            Un(2,3)=Ucomp(3,node2);
#            
#            # intial postion
#            x0=zeros(nxi,nsd);
#            for n=1:nxi
#                x0(n,:)=((1-xi(n))/2)*xn(:,node1)'+((1+xi(n))/2)*xn(:,node2)';
#            end
#            
#            Le=sqrt((xn(1,node1)-xn(1,node2))^2+(xn(2,node1)-xn(2,node2))^2);
#            
#            # g(i,j) - vector i coordinate j
#            g=zeros(nsd,nsd);
#            g(1,:)=(1/Le)*(xn(:,node2)-xn(:,node1))';
#            
#            g(2,1)=-g(1,2);
#            g(2,2)=g(1,1);
#            
#            # to make g_3 = z - so that the rotation are the same in both systems of coordinates
#            cross_prod=g(1,1)*g(2,2)-g(1,2)*g(2,1);
#            if cross_prod < 0
#                g(2,:)=-g(2,:);
#            end;
#            
#            #transformation matrices
#            beta=g;
#            betap=inv(beta);
#            
#            # transform displacement coordinates in the global system
#            un=zeros(nen,nsd);
#            
#            for n=1:nen
#                for i=1:nsd
#                    for j=1:nsd
#                        un(n,i)=un(n,i)+Un(n,j)*betap(j,i);
#                    end;
#                end;
#            end;
#            
#            for n=1:nen
#                un(n,3)=Un(n,3);
#            end;
#            
#            # transverse displacement in the local system of coordinates
#            for i=1:nxi
#                N1(i)=(1/4)*(1-xi(i))^2*(2+xi(i));
#                N2(i)=(Le/8)*(1-xi(i)^2)*(1-xi(i));
#                N3(i)=(1/4)*(1+xi(i))^2*(2-xi(i));
#                N4(i)=(Le/8)*(-1+xi(i)^2)*(1+xi(i));
#                u(i,2)=N1(i)*un(1,2)+N2(i)*un(1,3)+N3(i)*un(2,2)+N4(i)*un(2,3);
#            end;
#            
#            # horizontal displacement in the local system of coordinates
#            for i=1:nxi
#                u(i,1)=(1/2)*((1-xi(i))*un(1,1)+(1+xi(i))*un(2,1));
#            end;
#            
#            # transform from local to global system
#            U=zeros(nxi,nsd);
#            for i=1:nsd
#                for n=1:nxi
#                    for j=1:nsd
#                        U(n,i)=U(n,i)+u(n,j)*beta(j,i);
#                    end;
#                end;
#            end;
#            
#            # deformed configuration in the global system
#            xt=x0+scale*U;
#            
#            plot(xt(:,1),xt(:,2),'b-');
#        end;
#        
#    otherwise
#        disp('type not supported by plot_results');
#end

########################################
# boundary conditions on displacements #
########################################
#function plot_bc_displacements(type,Lcar,display_factor,nnp,Idb,xn,nsd)
#    alpha=display_factor*Lcar; # scale factor  for bc symbols
#    for P=1:nnp
#        if (nsd >1)
#            if ((Idb(1,P) ~= 0) && (Idb(2,P) ~= 0))
#                bc_symbols(xn(:,P),alpha,3,nsd);
#            end;
#            if ((Idb(1,P) ~= 0) && (Idb(2,P) == 0))
#                bc_symbols(xn(:,P),alpha,2,nsd);
#            end;
#            if ((Idb(1,P) == 0) && (Idb(2,P) ~= 0))
#                bc_symbols(xn(:,P),alpha,1,nsd);
#            end;
#        else
#            if (Idb(1,P) ~= 0)
#                bc_symbols([xn(1,P),0],alpha,2,nsd);
#            end;
#        end;
#    end;
#end;


#################################
## boundary conditions on force #
#################################
#function plot_bc_force(type,Lcar,f,display_factor,nnp,xn,nsd)
#delta=display_factor*Lcar/max(max(abs(f))); # scale factor for force b.c.
#alpha=2*display_factor*Lcar; # scale factor for moment
#for N=1:nnp
#    if ( nsd == 3)
#        if ( (f(1,N) ~= 0 ) || (f(2,N) ~= 0) || (f(3,N) ~= 0))
#            quiver3(xn(1,N),xn(2,N),xn(3,N),f(1,N),f(2,N),f(3,N),delta,'r','LineWidth',2);
#        end;
#    end
#    if ( nsd == 2)
#        if ( (f(1,N) ~= 0 ) || (f(2,N) ~= 0))
#            quiver(xn(1,N),xn(2,N),f(1,N),f(2,N),delta,'r');
#        end;
#    end
#    if ( nsd == 1) 
#        if (f(1,N) ~=0)
#            quiver(xn(1,N),0,f(1,N),0,delta,'r');
#        end;
#    end;
#end;

#if (strcmp(type, 'beam'))
#    for N=1:nnp
#        if (f(3,N) > eps)
#            plot_moments([xn(1,N);xn(2,N)],alpha,0);
#        end;
#        if (f(3,N) < -eps)
#            plot_moments([xn(1,N);xn(2,N)],alpha,1);
#        end;
#    end;
#end;


##############
## reactions #
##############
#function plot_reactions(type,Lcar,Rcomp,Idb,display_factor,xn,nnp,ndf,nsd)
#beta=3*display_factor*Lcar/max(max(abs(Rcomp))); # scale factor for reactions
#alpha=2*display_factor*Lcar; # scale factor for moment
#for N=1:nnp
#    RN=zeros(ndf);
#    if (nsd == 3)
#        if ( (Idb(1,N) ~= 0 ) || (Idb(2,N) ~= 0) || (Idb(3,N) ~= 0))
#            h = quiver3(xn(1,N),xn(2,N),xn(3,N),Rcomp(1,N),Rcomp(2,N),Rcomp(3,N),beta,'r-o','LineWidth',2);
#            adjust_quiver_arrowhead_size(h, 0.5/beta);
#        end;
#    end
#    if (nsd == 2)
#        if ((Idb(1,N) ~= 0) || (Idb(2,N) ~= 0))
#            h = quiver(xn(1,N),xn(2,N),Rcomp(1,N),Rcomp(2,N),beta,'r-o', 'LineWidth', 2);
#            adjust_quiver_arrowhead_size(h, 0.5/beta);
#        end;
#        if (strcmp(type, 'beam') && (Idb(3,N) ~= 0))
#            if (Rcomp(3,N) > eps)
#                plot_moments([xn(1,N);xn(2,N)],alpha,0);
#            end;
#            if (Rcomp(3,N) < -eps)
#                plot_moments([xn(1,N);xn(2,N)],alpha,1);
#            end;
#        end;
#    end
#    if (nsd == 1)
#        if (Idb(1,N) ~= 0)
#            h = quiver(xn(1,N),0,Rcomp(1,N),0,beta,'k', 'LineWidth', 2);
#            adjust_quiver_arrowhead_size(h, 0.5/beta);
#        end;
#    end;
#end

###############################################################################################
##                              BOUNDARY CONDITIONS SYMBOLS                                   #
###############################################################################################
#function bc_symbols(xp,alpha,symbol,nsd)
#if (nsd < 3)
#    switch symbol
#        case 1
#            # v fixed
#            x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
#            y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];

#            line(x,y,'Color','k','LineWidth',1.2);

#            for i=0:3,
#                circle([xp(1)-(3/8)*alpha+i*alpha/4; xp(2)-(7/8)*alpha],alpha/8);
#            end;


#        case 2
#            # u fixed
#            x=[xp(1);xp(1)-(3/4)*alpha;xp(1)-(3/4)*alpha;xp(1)];
#            y=[xp(2);xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha;xp(2)];

#            line(x,y,'Color','k','LineWidth',1.2);

#            for i=0:3,
#                circle([xp(1)-(7/8)*alpha;xp(2)-(3/8)*alpha+i*alpha/4],alpha/8);
#            end;

#        case 3
#            # u and v fixed
#            x=[xp(1);xp(1)-alpha/2;xp(1)+alpha/2;xp(1)];
#            y=[xp(2);xp(2)-(3/4)*alpha;xp(2)-(3/4)*alpha;xp(2)];

#            line(x,y,'Color','k','LineWidth',1.2);

#            for i=0:3,
#                line([xp(1)-(alpha/4)+i*alpha/4;xp(1)-(alpha/2)+i*(alpha/4)], ...
#                    [xp(2)-(3/4)*alpha;xp(2)-alpha],'Color','k','LineWidth',1.2);
#            end;

#        case 4
#            # v and theta fixed
#            x=[xp(1)-alpha/2;xp(1)+alpha/2];
#            y=[xp(2);xp(2)];

#            line(x,y,'Color','k','LineWidth',1.2);

#            for i=0:3,
#                circle([xp(1)-(3/8)*alpha+i*alpha/4; xp(2)-(1/8)*alpha],alpha/8);
#            end;

#        case 5
#            # u and theta fixed
#            x=[xp(1);xp(1)];
#            y=[xp(2)+(1/2)*alpha;xp(2)-(1/2)*alpha];

#            line(x,y,'Color','k','LineWidth',1.2);

#            for i=0:3,
#                circle([xp(1)-(1/8)*alpha;xp(2)-(3/8)*alpha+i*alpha/4],alpha/8);
#            end;

#        case 6
#            # u, v and theta fixed
#            line([xp(1)-alpha/2;xp(1)+alpha/2],[xp(2),xp(2)],'Color','k','LineWidth',1.2);
#            for i=0:3,
#                line([xp(1)-alpha/2+(i+1)*alpha/4, xp(1)-alpha/2+i*alpha/4],[xp(2),xp(2)-alpha/4]...
#                    ,'Color','k','LineWidth',1.2);
#            end;
#    end
#end

function circle(x0,r)
    theta=0:0.1:2*pi;
    x=r*cos(theta)+x0(1);
    y=r*sin(theta)+x0(2);

    plot(x,y,color="Black",linewidth=1.2);
end;



###############################################################################################
##                                       MOMENTS                                              #
###############################################################################################
#function plot_moments(xp,alpha,sign)
#switch sign
#    case 0  #positive moment
#        r=alpha/4;
#        theta=0:0.1:3*pi/2;
#        x=r*cos(theta)+xp(1);
#        y=r*sin(theta)+xp(2);
#        
#        plot(x,y,'k','LineWidth',1.2);
#        line([xp(1), xp(1)-alpha/8],[xp(2)-alpha/4, xp(2)-alpha/8],...
#            'Color','k','LineWidth',1.2);
#        line([xp(1), xp(1)-alpha/8],[xp(2)-alpha/4, xp(2)-3*alpha/8],...
#            'Color','k','LineWidth',1.2);
#        
#    case 1  # negative moment
#        r=alpha/4;
#        theta=pi:-0.1:-pi/2;
#        x=r*cos(theta)+xp(1);
#        y=r*sin(theta)+xp(2);
#        
#        plot(x,y,'k','LineWidth',1.2);
#        line([xp(1), xp(1)+alpha/8],[xp(2)-alpha/4, xp(2)-alpha/8],...
#            'Color','k','LineWidth',1.2);
#        line([xp(1), xp(1)+alpha/8],[xp(2)-alpha/4, xp(2)-3*alpha/8],...
#            'Color','k','LineWidth',1.2);
#end;


end # module
