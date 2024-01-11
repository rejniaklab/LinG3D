function tumorGrowth_example1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a companion code for the paper "LinG3D: Visualizing the  %%
%% Spatio-Temporal Dynamics of Clonal Evolution" by A. Hu, A.M.E.   %%
%% Ojwang', K.D. Olumoyin, and K.A. Rejniak                         %%
%% This code generates tumor growth data 'cell_history', 'cellXY_##'%%
%% and 'cellID_##' (where ## is the iteration number) used in all   %%
%% routines for visualization of the 3D linage trees of the tumor.  %%
%%                                                                  %%
%% for the examples discussed in the paper use:                     %%
%% example 1: pathdata='exampleB05';  probmut=0.050;                %%
%% example 2: pathdata='exampleB005'; probmut=0.005;                %%
%%                                                                  %%
%% October 31, 2022                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% example 1
probmut=0.050;            % prob of mutation [0.05 or 0.005]
pathname='exampleB05';

% example 2
%probmut=0.005;           % prob of mutation [0.05 or 0.005]
%pathname='exampleB005';


% switch for saving data: toSave=1 save data; toSave=0 do not save data
toSave=1;
if toSave==1
  dirnamefigs=[pathname,'/figs']; mkdir(dirnamefigs) 
  dirnamedata=[pathname,'/data']; mkdir(dirnamedata)
end


 rng(7); % fixed random number generator for comparison between examples


Ncells=1;                % number of cells
cellXY=[0,0];            % cell coordinates [microns] 
Amat=18*60;              % cell average maturation age [sec]
age=[0.75*Amat];         % cell current age
ageDiv=[Amat];           % cell division age
drugcell=[0];            % drug level absorbed by the cell
mutate=[0];              % cell mutation state

rad=5;                   % cell radius microns
Nneig=5;                 % number of neighbors for overcrowding   
probdie=0.8;             % prob of random dying
probDrugdie=0.85;        % prob of drug-induced dying
Nmutate=0;               % number of mutation 
Cellupt=0.015;           % cell uptake
drugth=4;                % threshold for cell killing level

stif=50;                 % springs stiffness 
nu=250;                  % viscosity of medium
                    
Nvess=5;                 % number of vessel
vess=[-45,-45;-15,14;-65,65;75,35;80,-80];  % coordinate of each vessel 
Infvess=[5,5,5,5,5];     % influx of the vessel 
Radvess=10;              % radius of the vessel


%domain 
xmin=-100; xmax=100; ymin=-100; ymax=100;  % computation al domain
hg=5;                                      % grid width 
xxg=xmin:hg:xmax; yyg=ymin:hg:ymax;        % grid for drawing contours
Ngx=length(xxg); Ngy=length(yyg);          % size of grid in x- & y-

% drug
drug=zeros(Ngx+2,Ngy+2); % drug domain
diff=10;                 % drug diffusion coefficient
step=3;                  % radius of sensing the drug by the cell


dt=0.25;                 % time step [sec]
Niter=100000;            % number of iterations
Ndrugi=60000;            % initial iteration for drug injection
Ndraw=250;               % frequency of drawning and saving


% data required for drawing the 3D trees
cellID=[1];              % unigue cell index
NcellID=1;               % curent maximal cell index number
cell_history=[1,0,0,0,0]; % data format needed for defining the 3D trees
  % [cell ID, clone ID, mother ID, birth iter, div/death iter]
 


%draw the figure 
figure ('position',[250,100,900,750]);


% main loop
for iter=0:Niter
  
   %boundary conditions
   drug(1,:)=drug(2,:);
   drug(Ngx+2,:)=drug(Ngx+1,:);
   drug(:,1)=drug(:,2);
   drug(:,Ngy+2)=drug(:,Ngy+1);
    
   %influx of drug 
   if iter>=Ndrugi
     for k=1:Nvess
       vi=2+floor((vess(k,1)-xmin)/hg); 
       vj=2+floor((vess(k,2)-ymin)/hg);
       for i=-3:3
         for j=-3:3
           if (vi+i+1<Ngx+2)&&(vi+i+1>0)&&(vj+j+1<Ngy+2)&&(vj+j+1>0)
             dx=xmin+hg*(vi+i-1)-vess(k,1); dy=ymin+hg*(vj+j-1)-vess(k,2);
             dxy=sqrt(dx^2+dy^2);
             if dxy<Radvess drug(vi+i+1,vj+j+1)=Infvess(k)*dt; end  
           end % if 
         end % for j
       end % for i
     end % for k
   end %for iter>
    
   %diffusion
   for i=2:Ngx+1
     for j=2:Ngy+1
       drug(i,j)=drug(i,j)+(diff*dt/(hg*hg))*(drug(i-1,j)+drug(i+1,j)+...
                 drug(i,j+1)+drug(i,j-1)-4*drug(i,j));
     end % for j
   end % for i
    
   %drug uptake
   for k=1:Ncells
     ci=2+floor((cellXY(k,1)-xmin)/hg);
     cj=2+floor((cellXY(k,2)-ymin)/hg);
     for i=-3:3
       for j=-3:3
         if (ci+i+1<Ngx+2)&&(ci+i+1>0)&&(cj+j+1<Ngy+2)&&(cj+j+1>0)
           dx=xmin+hg*(ci+i-1)-cellXY(k,1); dy=ymin+hg*(cj+j-1)-cellXY(k,2);
           dxy=sqrt(dx^2+dy^2);
           if dxy<rad
             drugcell(k)=drugcell(k)+dt*Cellupt*drug(ci+i+1,cj+j+1);
             drug(ci+i+1,cj+j+1)=max(0,drug(ci+i+1,cj+j+1)*(1-dt*Cellupt));
           end % if dxy
         end % if
       end % for j
     end % for i 
   end % for k
 
    
   % cell random death
   Ncells2=0;
   for k=1:Ncells
     if (rand<probdie)&&(age(k)>2*ageDiv(k)) 
       cell_history(cellID(k),5)=iter;     
       % [cell ID, clone ID, mother ID, birth iter, div/death iter]
     else
       Ncells2=Ncells2+1;
       cellXY2(Ncells2,1:2)=cellXY(k,1:2); age2(Ncells2)=age(k); 
       ageDiv2(Ncells2)=ageDiv(k); mutate2(Ncells2)=mutate(k); 
       drugcell2(Ncells2)=drugcell(k); cellID2(Ncells2)=cellID(k);        
     end  
   end
   if Ncells2>0
     clear Ncells cellXY age ageDiv mutate drugcell cellID 
     Ncells=Ncells2;
     cellXY=cellXY2; age=age2; ageDiv=ageDiv2;
     mutate=mutate2; drugcell=drugcell2; cellID=cellID2;
     clear Ncells2 cellXY2 age2 ageDiv2 mutate2 drugcell2 cellID2  
   end

   % drug induced death
   Ncells2=0;
   for k=1:Ncells
     if  (mutate(k)==0)&&(drugcell(k)>drugth)&&(rand>probDrugdie)
       cell_history(cellID(k),5)=iter;      
       % [cell ID, clone ID, mother ID, birth iter, div/death iter]
     else  
       Ncells2=Ncells2+1;
       cellXY2(Ncells2,1:2)=cellXY(k,1:2); age2(Ncells2)=age(k);
       ageDiv2(Ncells2)=ageDiv(k); mutate2(Ncells2)=mutate(k);
       drugcell2(Ncells2)=drugcell(k); cellID2(Ncells2)=cellID(k); 
     end 
   end
   if Ncells2>0
     clear Ncells cellXY age ageDiv mutate drugcell cellID  
     Ncells=Ncells2;
     cellXY=cellXY2; age=age2; ageDiv=ageDiv2;
     mutate=mutate2; drugcell=drugcell2; cellID=cellID2;
     clear Ncells2 cellXY2 age2 ageDiv2 mutate2 drugcell2 cellID2 
   end
    
    
   age=age+dt;          % increase cell age
  

   % calculate nearby cell
   neig=zeros(Ncells,1);
   for i=1:Ncells-1
     for j=i+1:Ncells
       dx=cellXY(i,1)-cellXY(j,1); dy=cellXY(i,2)-cellXY(j,2);     
       dxy=sqrt(dx^2+dy^2);
       if (dxy<3*rad)&&(dxy>0) 
         neig(i,1)=neig(i,1)+1; neig(j,1)=neig(j,1)+1; 
       end 
     end % if j
   end % if i
 
  
   % cell division
   for ij=1:Ncells
     if (age(ij)>ageDiv(ij))&&(neig(ij)<Nneig)
       cell_history(cellID(ij),5)=iter-1;    
        % [cell ID, clone ID, mother ID, birth iter, div/death iter]

       Ncells=Ncells+1;   
       alpha=2*pi*rand;
       cellXY(Ncells,1:2)= cellXY(ij,1:2)+0.5*rad*[cos(alpha),sin(alpha)];
       age(Ncells)=0;
       ageDiv(Ncells)=ageDiv(ij)+(rand-0.5)*10;
       neig(Ncells)=0;
       drugcell(Ncells)=0;
       mutate(Ncells)=mutate(ij);
        
       NcellID=NcellID+1;                   
       cell_history(NcellID,1:5)=[NcellID,mutate(ij),cellID(ij),iter,0];  
       cellID(Ncells)=NcellID;               

       age(ij)=0;
       ageDiv(ij)=ageDiv(ij)+(rand-0.5)*10;
       neig(ij)=0;
       drugcell(ij)=0;
          
       NcellID=NcellID+1;                   
       cell_history(NcellID,1:5)=[NcellID,mutate(ij),cellID(ij),iter,0];  
       cellID(ij)=NcellID;                  
        
       if (rand<probmut)&&(iter>=0.8*Ndrugi) %prob of mutation  
         ageDiv(Ncells)=0.5*Amat;
         mutate(Ncells)=Nmutate+1;
         ageDiv(ij)=0.5*Amat;
         mutate(ij)=Nmutate+1;
         Nmutate=Nmutate+1;
          
         cell_history(NcellID-1,2)=Nmutate;  
         cell_history(NcellID  ,2)=Nmutate;  
         % [cell ID, clone ID, mother ID, birth iter, div/death iter]
       end     
     end
   end
  
  
    % cell-cell repulsive forces
    force=zeros(Ncells,2);
    for i=1:Ncells-1
      for j=i+1:Ncells
        dx=cellXY(i,1)-cellXY(j,1);
        dy=cellXY(i,2)-cellXY(j,2);  
        dxy=sqrt(dx^2+dy^2);
        if (dxy<2*rad)&&(dxy>0)
          force(i,1:2)=force(i,1:2)+stif*(2*rad-dxy)*[dx,dy]/dxy;
          force(j,1:2)=force(j,1:2)-stif*(2*rad-dxy)*[dx,dy]/dxy;
        end % if dxy
      end % if j
    end % if i
  
    % vessel-cell repulsive forces
    for i=1:Ncells
      for j=1:Nvess
        dx=cellXY(i,1)-vess(j,1);
        dy=cellXY(i,2)-vess(j,2);  
        dxy=sqrt(dx^2+dy^2);
        if (dxy<(rad+Radvess))&&(dxy>0)
          force(i,1:2)=force(i,1:2)+2*stif*(rad+Radvess-dxy)*[dx,dy]/dxy;
        end % if dxy
      end % if j
    end % if i
  
    
    % cell relocation
    cellXY=cellXY+dt*force/nu;
  
    % remove cells outside the domain
    Ncells2=0;
    for k=1:Ncells
      if (cellXY(k,1)>xmin)&&(cellXY(k,1)<xmax)&&(cellXY(k,2)>ymin)&&(cellXY(k,2)<ymax) 
        Ncells2=Ncells2+1;
        cellXY2(Ncells2,1:2)=cellXY(k,1:2);
        age2(Ncells2)=age(k);
        ageDiv2(Ncells2)=ageDiv(k);
        mutate2(Ncells2)=mutate(k);
        drugcell2(Ncells2)=drugcell(k);        
        cellID2(Ncells2)=cellID(k);       
      else
        cell_history(cellID(k),5)=iter;     
        % [cell ID, clone ID, mother ID, birth iter, div/death iter]
      end
    end    
    if Ncells2>0
      clear Ncells cellXY age ageDiv mutate drugcell cellID  
      Ncells=Ncells2;
      cellXY=cellXY2; age=age2; ageDiv=ageDiv2;
      mutate=mutate2; drugcell=drugcell2; cellID=cellID2;          
      clear Ncells2 cellXY2 age2 ageDiv2 mutate2 drugcell2 cellID2
    end
     
        
  
    % draw and save the data
    if (mod(iter,Ndraw)==0)||(iter==Niter)
       draw_image(xxg,yyg,drug,Ngx,Ngy,Infvess,dt,Nvess,vess,Radvess,...
       cellXY,rad,Ncells,mutate,iter,Niter,probmut,xmin,xmax,ymin,ymax)
       
              
       if toSave==1
         filename=[dirnamefigs,'/fig_',num2str(iter),'.jpg'];
         print(filename,'-djpeg')
        
         filename=[dirnamedata,'/cellXY_',num2str(iter),'.txt'];
         save(filename,'cellXY','-ascii')
         cellID=cellID';
         filename=[dirnamedata,'/cellID_',num2str(iter),'.txt'];
         save(filename,'cellID','-ascii')
       end   
     end  % if mod  
   end  % for iter
    
  
   % save data for the 3D lineage trees
   if toSave==1
     filename=[dirnamedata,'/cell_history.txt'];
     save(filename,'cell_history','-ascii')
    
     filename=[dirnamedata,'/drug.txt'];
     save(filename,'drug','-ascii')
   end    
        
end  % end of "function"

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%-----------------------------------------------------------------------
function draw_image (xxg,yyg,drug,Ngx,Ngy,Infvess,dt,Nvess,vess,Radvess,...
cellXY,rad,Ncells,mutate,iter,Niter,probmut,xmin,xmax,ymin,ymax)

    clf
    
    %draw drug distribution
    if max(max(drug))>0 
      contourf(xxg,yyg,drug(2:Ngx+1,2:Ngy+1)',[0:0.1:max(Infvess)],'edgecolor','none');
    else
      fill([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],'b')
    end
    colormap(jet)
    colorbar
    caxis([0,max(Infvess)*dt])
    hold on
    axis equal
    
 
    % draw vessels
    viscircles(vess,Radvess*ones(Nvess,1),'color','r');
    plot(vess(:,1),vess(:,2),'ro','markersize',45,'markerfacecolor',[204,0,0]/255)
 
    % draw cells
    color=DefineNewColorPalette;
    Ncolor=size(color,1); 
    for jj=0:max(mutate)  
      ind=find(mutate==jj);
      colorj=mod(jj,Ncolor)+1;
      viscircles(cellXY(ind,1:2),rad*ones(length(ind),1),'color',color(colorj,1:3));
    end
    
    axis([xmin,xmax,ymin,ymax])
    
    title(['iteration ',num2str(iter),' of ',num2str(Niter),...
      ';   # of cells: ',num2str(Ncells),';   # of mutations: ',...
      num2str(max(mutate)),';   P^{mut}=',num2str(probmut)],'FontSize',15)
   pause(0.1)
    
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function color=DefineNewColorPalette
    
  color=[255,  0,255;255,  0,  0;  0,255,255;  0,  0,255;  0,255,  0;...
           0,  0,  0;255,191,  0;255,255,  0;191,255,  0;128,128,  0;...
         255,182,193;  0,191,255;  0,128,255;250,235,215;128,  0,255;...
         154,205, 50;255,  0,128;102,  0,  0;102, 77,  0;  0,102,102;...
         204,204,255;255,204,255;153,204,255;255,153,153;  0,153,  0;...
           0,153,153;153,  0, 77;255,228,225;128,  0,  0;102,102,153;...
         153,255,204;218,112,214;255,128,  0;192,192,192;128,128,128;...
          75,  0,130;165, 42, 42;216,191,216;220, 20, 60;245,222,179;...
         255, 99, 71;255,127, 80;205, 92, 92;240,128,128;233,150,122;...
         250,128,114;255,160,122;255, 69,  0;255,140,  0;255,165,  0;...
         255,215,  0;184,134, 11;218,165, 32;  0,100,  0;255,240,245;...
         188,143,143;255,248,220; 50,205, 50;144,238,144;152,251,152;...
         143,188,143;  0,250,154;  0,255,127; 46,139, 87;102,205,170;...
          60,179,113; 32,178,170; 47, 79, 79;  0,128,128;  0,139,139;...
         240,230,140;245,245,220;224,255,255;  0,206,209;255,228,181;...
         255, 20,147;175,238,238;127,255,212;176,224,230; 95,158,160;...
          70,130,180;100,149,237;222,184,135; 30,144,255;238,232,170;...
         189,183,107;107,142, 35;124,252,  0;127,255,  0;173,255,  47;...
         178, 34, 34;221,160,221;255,235,205]/255;          
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%-----------------------------------------------------------------------








