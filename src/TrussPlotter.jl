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
    plot_bc_displacements(Lcar,display_factor,nnp,Idb,xn,nsd);
    plot_bc_force(Lcar,f,display_factor,nnp,xn,nsd);
end

function plot_deformed(xn,f,Idb,Ucomp,Rcomp,ien,nel,nen,nsd,ndf,nnp,axial)
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


    axis("equal");
    title("Undeformed and deformed mesh");
    plot_mesh_undeformed(nel,ien,xn,nnp,nsd);
    plot_mesh_deformed(xn,Ucomp,Idb,display_factor,Lcar,nel,ien,ndf,nsd,nen,axial);
    numbers(nel,ien,xn,nnp,nsd);
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
    for elt=1:nel
        node1=ien[1,elt];
        node2=ien[2,elt];
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
            plot(x0,y0,z0);
        end
    end;
end;



###################################################
# numbers the nodes and elements - truss and beam #
###################################################
function numbers(nel,ien,xn,nnp,nsd)
    fontsize = 10;
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
end;

##########################
# deformed configuration #
##########################
function plot_mesh_deformed(xn,Ucomp,Idb,display_factor,Lcar,nel,ien,ndf,nsd,nen,axial)
    scale=display_factor*Lcar/maximum(maximum(abs(Ucomp))); # scale factor for the displacements
    legend_comp_flag = 0;
    legend_ten_flag = 0;
    h = cell(2);
    xt = zeros(ndf,2);
    for elt = 1:nel
        node1 = ien[1,elt];
        node2 = ien[2,elt];
        xt[:,1]=xn[:,node1];
        xt[:,2]=xn[:,node2];
        
        if (nsd == 1)
            xt[2,1]=0;
            xt[2,2]=0;
        end;
        
        for i = 1:ndf
            if Idb[i,node1] != 0
                xt[i,1] = xt[i,1] + scale*Ucomp[i,node1];
            else
                xt[i,1] = xt[i,1] + scale*Ucomp[i,node1];
            end;
            if Idb[i,node2] != 0
                xt[i,2] = xt[i,2] + scale*Ucomp[i,node2];
            else
                xt[i,2] = xt[i,2] + scale*Ucomp[i,node2];
            end;
        end;
        if nsd < 3
            if (axial[2,elt] > 0.0)
               ten = plot(xt[1,:]',xt[2,:]', color="Red", marker=".", label="Tension");
               if legend_ten_flag == 0
                   legend_ten_flag = 1;
               end
            else
               comp = plot(xt[1,:]',xt[2,:]', color="Blue", marker=".", label="Compression");
               if legend_comp_flag == 0
                   legend_comp_flag = 1;
               end
            end
        end
        if nsd == 3
            if (axial[2,elt] > 0.0)
               plot(xt[1,:]',xt[2,:]',xt[3,:]', color="Red", marker=".");
               if legend_ten_flag == 0
                   legend_ten_flag = 1;
               end
            else
               plot(xt[1,:]',xt[2,:]',xt[3,:]', color="Blue", marker=".");
               if legend_comp_flag == 0
                   legend_comp_flag = 1;
               end
            end
        end
    end
    #legend(numpoints=2);
end;

########################################
# boundary conditions on displacements #
########################################
function plot_bc_displacements(Lcar,display_factor,nnp,Idb,xn,nsd)
    alpha=display_factor*Lcar; # scale factor  for bc symbols
    for P=1:nnp
        if (nsd >1)
            if ((Idb[1,P] != 0) && (Idb[2,P] != 0))
                bc_symbols(xn[:,P],alpha,3,nsd);
            end;
            if ((Idb[1,P] != 0) && (Idb[2,P] == 0))
                bc_symbols(xn[:,P],alpha,2,nsd);
            end;
            if ((Idb[1,P] == 0) && (Idb[2,P] != 0))
                bc_symbols(xn[:,P],alpha,1,nsd);
            end;
        else
            if (Idb[1,P] != 0)
                bc_symbols([xn[1,P],0],alpha,2,nsd);
            end;
        end;
    end;
end;


################################
# boundary conditions on force #
################################
function plot_bc_force(Lcar,f,display_factor,nnp,xn,nsd)
    delta=display_factor*Lcar/maximum(maximum(abs(f))); # scale factor for force b.c.
    alpha=2*display_factor*Lcar; # scale factor for moment
    for N=1:nnp
        if ( nsd == 3)
            if ( (f[1,N] != 0 ) || (f[2,N] != 0) || (f[3,N] != 0))
                #quiver3(xn(1,N),xn(2,N),xn(3,N),f(1,N),f(2,N),f(3,N),delta,'r','LineWidth',2);
            end;
        end
        if ( nsd == 2)
            if ( (f[1,N] != 0 ) || (f[2,N] != 0))
                quiver(xn[1,N],xn[2,N],f[1,N],f[2,N],delta,color="Red");
            end;
        end
        if ( nsd == 1)
            if (f[1,N] !=0)
                quiver(xn[1,N],0,f[1,N],0,delta,color="Red");
            end;
        end;
    end;
end



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

##############################################################################################
#                              BOUNDARY CONDITIONS SYMBOLS                                   #
##############################################################################################
function circle(x0,r)
    theta=0:0.1:2*pi;
    x=r*cos(theta).+x0[1];
    y=r*sin(theta).+x0[2];

    plot(x,y,color="Black",linewidth=1.2);
end;

function bc_symbols(xp,alpha,symbol,nsd)
    if symbol== 1
        # v fixed
        x=[xp[1];xp[1]-alpha/2;xp[1]+alpha/2;xp[1]];
        y=[xp[2];xp[2]-(3/4)*alpha;xp[2]-(3/4)*alpha;xp[2]];

        plot(x,y,color="Black",linewidth=1.2);

        for i=0:3
            circle([xp[1]-(3/8)*alpha+i*alpha/4; xp[2]-(7/8)*alpha],alpha/8);
        end;


    elseif symbol == 2
        # u fixed
        x=[xp[1];xp[1]-(3/4)*alpha;xp[1]-(3/4)*alpha;xp[1]];
        y=[xp[2];xp[2]+(1/2)*alpha;xp[2]-(1/2)*alpha;xp[2]];

        plot(x,y,color="Black",linewidth=1.2);

        for i=0:3
            circle([xp[1]-(7/8)*alpha;xp[2]-(3/8)*alpha+i*alpha/4],alpha/8);
        end;

    elseif symbol == 3
        # u and v fixed
        x=[xp[1];xp[1]-alpha/2;xp[1]+alpha/2;xp[1]];
        y=[xp[2];xp[2]-(3/4)*alpha;xp[2]-(3/4)*alpha;xp[2]];

        plot(x,y,color="Black",linewidth=1.2);

        for i=0:3
            plot([xp[1]-(alpha/4)+i*alpha/4;xp[1]-(alpha/2)+i*(alpha/4)], 
                 [xp[2]-(3/4)*alpha;xp[2]-alpha],color="Black",linewidth=1.2);
        end;

    elseif symbol == 4
        # v and theta fixed
        x=[xp[1]-alpha/2;xp[1]+alpha/2];
        y=[xp[2];xp[2]];

        plot(x,y,color="Black",linewidth=1.2);

        for i=0:3
            circle([xp[1]-(3/8)*alpha+i*alpha/4; xp[2]-(1/8)*alpha],alpha/8);
        end;

    elseif symbol == 5
        # u and theta fixed
        x=[xp[1];xp[1]];
        y=[xp[2]+(1/2)*alpha;xp[2]-(1/2)*alpha];

        plot(x,y,color="Black",linewidth=1.2);

        for i=0:3
            circle([xp[1]-(1/8)*alpha;xp[2]-(3/8)*alpha+i*alpha/4],alpha/8);
        end;

    elseif symbol == 6
        # u, v and theta fixed
        plot([xp[1]-alpha/2;xp[1]+alpha/2],[xp[2],xp[2]], color="Black", linewidth=1.2);
        for i=0:3
            plot([xp[1]-alpha/2+(i+1)*alpha/4, xp[1]-alpha/2+i*alpha/4],[xp[2],xp[2]-alpha/4]
                ,color="Black",linewidth=1.2);
        end;
    end
end;

end # module
