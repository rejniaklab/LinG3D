function linG3DClone(pathData,cloneNum,IsGradient,xmin,xmax,...
ymin,ymax,tmin,tmax,fileStep,toPrint)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a companion code for the paper "LinG3D: Visualizing the  %%
%% Spatio-Temporal Dynamics of Clonal Evolution" by A. Hu, A.M.E.   %%
%% Ojwang', K.D. Olumoyin, and K.A. Rejniak                         %%
%% This code generates the 3D lineage tree of all cells from one    %%
%% clone of number specified in 'cloneNum' with data from the       %%
%% directory 'pathData'.                                            %%
%%                                                                  %%
%% The following parameters need to be specified:                   %%
%%   pathData  -- directory with input data                         %%
%%   cloneNum  -- clone number to be drawn                          %%
%%   IsGradient -- 1 to draw drug in the background, 0 not to draw  %%
%%   xmin,xmax,ymin,ymax -- dimensions of the spacial domain        %%
%%   tmin, tmax          -- dimensions of the temporal domain       %%
%%   fileStep            -- frequency of the sampled data           %% 
%%   toPrint    -- 1 to save the generated figure, 0 not to save    %% 
%% It requires the following data in the pathData/data/ directory:  %%
%%   cell_history.txt -- file with info about each cell             %%
%%   cellID_##.txt    -- cell IDs in a file with index number ##    %% 
%%   cellXY_##.txt    -- cell coordinates in a file with index ##   %%
%%   drug.txt         -- concentration of a drug for background     %%
%% for the examples discussed in the paper use:                     %%
%%   example 1: pathData='exampleB05';  cloneNum between 0 and 9    %%
%%   example 2: pathData='exampleB005'; cloneNum between 0 and 147  %%
%%   example 3: pathData='exampleExp';  numClones=10;               %%
%%                                                                  %%
%% January 10, 2024                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




 dataDirectory='/data/';  % directory with cell and drug data
 % create directory to save individual clone figures  
 if toPrint==1
   pathFigs=[pathData,'/fig_clones']; 
   if ~exist(pathFigs) mkdir(pathFigs); end
 end


 timeStep=(tmax-tmin)/(2.5*(xmax-xmin)); %frequency of sampled temporal data

 % draw 3D figure
 figure
 grid
 axis([xmin,xmax,0,tmax/timeStep,ymin,ymax])
 axis equal
 view(3)
 hold on
 col=DefineColorPalette; Ncol=size(col,1);
 
 % draw background with drug gradient
 if IsGradient==1
   drug=load([pathData,dataDirectory,'drug.txt']);
   DrawBackground(drug,tmax,timeStep,xmin,xmax,ymin,ymax)  
 end
 
 
 % load cell history file
 hist=load([pathData,dataDirectory,'cell_history.txt']); 
   % [cell ID, clone ID, mother ID, birth iter, div/death iter]
 

 disp(['clone=',num2str(cloneNum)]);
 

 % load indices of all cells from cloneNum
 indLast=find((hist(:,2)==cloneNum)&(hist(:,3)<=tmax));

 % define matrix of line segments (3D branches) to draw
 matrix_to_draw=zeros(1,6);  % [x1,k1,y1,x2,k2,y2]
 Nmatrix=0;

  for ii=1:length(indLast)  % for every cell with index in indLast
    if mod(ii,100)==0 disp('... calculating'); end   

    cellNum=hist(indLast(ii),1);  % cell ID
    mothNum=hist(indLast(ii),3);  % mother ID
    strtNum=hist(indLast(ii),4);  % cell birth
    endNum =hist(indLast(ii),5); 
    endNum=max(tmax,min(endNum,tmax)); % cell div/death/tmax 

    % find all appearances of the cellNum  
    kkStart=fileStep*floor(strtNum/fileStep); % initial file number 
    kkEnd  =fileStep*floor(endNum/fileStep);  % final file number
    for kk=kkEnd:-fileStep:kkStart+fileStep   % inspect all files
      % cell ID and cell XY from the first file    
      fileMeID =load([pathData,dataDirectory,'cellID_',num2str(kk),'.txt']);
      fileMeXY =load([pathData,dataDirectory,'cellXY_',num2str(kk),'.txt']);
      indMe =find(fileMeID==cellNum); % find current indices of cellID
      % cell ID and cell XY from the second file  
      fileMe2ID=load([pathData,dataDirectory,'cellID_',num2str(kk-fileStep),'.txt']);      
      fileMe2XY=load([pathData,dataDirectory,'cellXY_',num2str(kk-fileStep),'.txt']);
      indMe2=find(fileMe2ID==cellNum); % find current indices of cellID
      if isempty(indMe)
      elseif isempty(indMe2)
        while kkStart<hist(mothNum,4) % find file with the grand-mother cell
          mothNum=hist(mothNum,3);
        end
        fileMe2ID=load([pathData,dataDirectory,'cellID_',num2str(kkStart),'.txt']);      
        fileMe2XY=load([pathData,dataDirectory,'cellXY_',num2str(kkStart),'.txt']);
        indMe2=find(fileMe2ID==mothNum); % find current indices of mother cellID
        if isempty(indMe2)
        else
          Nmatrix=Nmatrix+1; % save branch to draw [x1,t1,y1,x2,t2,y2]
          matrix_to_draw(Nmatrix,1:6)=[fileMeXY(indMe,1),kkStart+fileStep,...
            fileMeXY(indMe,2),fileMe2XY(indMe2,1),kkStart,fileMe2XY(indMe2,2)];      
        end   
      else  
        Nmatrix=Nmatrix+1; % save branch to draw [x1,t1,y1,x2,t2,y2]
        matrix_to_draw(Nmatrix,1:6)=[fileMe2XY(indMe2,1),kk-fileStep,...
            fileMe2XY(indMe2,2),fileMeXY(indMe,1),kk,fileMeXY(indMe,2)];
      end
    end
  end

  % drawing clones
  NumCol=mod(cloneNum,Ncol)+1;  % clone color 

  title(['cells from clone=',num2str(cloneNum)])
  axis equal
  axis([xmin,xmax,tmin,tmax/timeStep,ymin,ymax])
  view(3)
  ylabel(['iterations/time x',num2str(timeStep)])
  hYLabel = get(gca,'YLabel');
  set(hYLabel,'rotation',-35) %,'VerticalAlignment','middle')
  hold on  
    
  for ii=1:Nmatrix
    plot3([matrix_to_draw(ii,1),matrix_to_draw(ii,4)],...
      [matrix_to_draw(ii,2),matrix_to_draw(ii,5)]/timeStep,...
      [matrix_to_draw(ii,3),matrix_to_draw(ii,6)],'color',...
      col(NumCol,1:3),'linewidth',1.5)
    if mod(ii,100)==0 pause(0.01); end    
  end
  view(3)
  axis equal
  axis([xmin,xmax,tmin,tmax/timeStep,ymin,ymax])
  pause(0.1)

  if toPrint==1
    print('-djpeg100',[pathFigs,'/tree_clone_',num2str(cloneNum),'.jpg'])
  end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  DrawBackground(drug,tmax,timeStep,xmin,xmax,ymin,ymax)  
 
  drugmin=min(min(drug)); drugmax=max(max(drug)); drugstep=(drugmax-drugmin)/4;
  
  kk=tmax/timeStep;
  [Nx,Ny]=size(drug); hgx=(xmax-xmin)/Nx; hgy=(ymax-ymin)/Ny;

  for ii=1:Nx
    for jj=1:Ny
      if (drug(ii,jj)>=drugmin)&&(drug(ii,jj)<drugmin+drugstep)
        patch([xmin+(ii-1)*hgx,xmin+ii*hgx,xmin+ii*hgx,xmin+(ii-1)*hgx,...
        xmin+(ii-1)*hgx],[kk,kk,kk,kk,kk],[ymin+(jj-1)*hgy,ymin+(jj-1)*hgy,...
        ymin+jj*hgy,ymin+jj*hgy,ymin+(jj-1)*hgy],'b','edgecolor','none')
        hold on
      elseif (drug(ii,jj)>=drugmin+drugstep)&&(drug(ii,jj)<drugmin+2*drugstep)
        patch([xmin+(ii-1)*hgx,xmin+ii*hgx,xmin+ii*hgx,xmin+(ii-1)*hgx,...
        xmin+(ii-1)*hgx],[kk,kk,kk,kk,kk],[ymin+(jj-1)*hgy,ymin+(jj-1)*hgy,...
        ymin+jj*hgy,ymin+jj*hgy,ymin+(jj-1)*hgy],'c','edgecolor','none')
        hold on
      elseif (drug(ii,jj)>=drugmin+2*drugstep)&&(drug(ii,jj)<drugmin+3*drugstep)
        patch([xmin+(ii-1)*hgx,xmin+ii*hgx,xmin+ii*hgx,xmin+(ii-1)*hgx,...
        xmin+(ii-1)*hgx],[kk,kk,kk,kk,kk],[ymin+(jj-1)*hgy,ymin+(jj-1)*hgy,...
        ymin+jj*hgy,ymin+jj*hgy,ymin+(jj-1)*hgy],'y','edgecolor','none')
        hold on
      else
        patch([xmin+(ii-1)*hgx,xmin+ii*hgx,xmin+ii*hgx,xmin+(ii-1)*hgx,...
        xmin+(ii-1)*hgx],[kk,kk,kk,kk,kk],[ymin+(jj-1)*hgy,ymin+(jj-1)*hgy,...
        ymin+jj*hgy,ymin+jj*hgy,ymin+(jj-1)*hgy],'r','edgecolor','none')
        hold on
      end
    end
  end
  alpha(0.25)
  hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function col=DefineColorPalette
  col=[255,0,255;255,0,0;0,255,255;0,0,255;0,255,0;0,0,0;255,191,0;...
       255,255,0;191,255,0;128,128,0;255,182,193;0,191,255;0,128,255;...
       250,235,215;128,0,255;154,205,50;255,0,128;102,0,0;102,77,0;...
       0,102,102;204,204,255;255,204,255;153,204,255;255,153,153;0,153,0;...
       0,153,153;153,0,77;255,228,225;128,0,0;102,102,153;153,255,204;...
       218,112,214;255,128,0;192,192,192;128,128,128;75,0,130;165,42,42;...
       216,191,216;220,20,60;245,222,179;255,99,71;255,127,80;205,92,92;...
       240,128,128;233,150,122;250,128,114;255,160,122;255,69,0;...
       255,140,0;255,165,0;255,215,0;184,134,11;218,165,32;0,100,0;...
       255,240,245;188,143,143;255,248,220;50,205,50;144,238,144;...
       152,251,152;143,188,143;0,250,154;0,255,127;46,139,87;102,205,170;...
       60,179,113;32,178,170;47,79,79;0,128,128;0,139,139;240,230,140;...
       245,245,220;224,255,255;0,206,209;255,228,181;255,20,147;...
       175,238,238;127,255,212;176,224,230;95,158,160;70,130,180;...
       100,149,237;222,184,135;30,144,255;238,232,170;189,183,107;...
       107,142,35;124,252,0;127,255,0;173,255,47;...
       178,34,34;221,160,221;255,235,205;]/255;             
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

