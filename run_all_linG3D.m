
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This is a companion code for the paper "LinG3D: Visualizing the  %%
%% Spatio-Temporal Dynamics of Clonal Evolution" by A. Hu, A.M.E.   %%
%% Ojwang', K.D. Olumoyin, and K.A. Rejniak                         %%
%%                                                                  %%
%% This code executes all LinG3D routines for 3 examples discussed  %%
%% in the paper.                                                    %%     
%% The following parameters need to be specified:                   %%
%%   pathData  -- directory with input data                         %%
%%   cloneNum  -- clone number to be drawn                          %%
%%   numClones  -- total number of clones in the data               %%
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


% syntaxt: 
% linG3DAll(pathData,numClones,IsGradient,xmin,xmax,ymin,ymax,tmin,tmax,fileStep,toPrint)
% linG3DAliveAll(pathData,numClones,IsGradient,xmin,xmax,ymin,ymax,tmin,tmax,fileStep,toPrint)
% linG3DClone(pathData,cloneNum,IsGradient,xmin,xmax,ymin,ymax,tmin,tmax,fileStep,toPrint)
% linG3DAliveClone(pathData,cloneNum,IsGradient,xmin,xmax,ymin,ymax,tmin,tmax,fileStep,toPrint)


% Example 1:
linG3DAll('exampleB005',9,1,-100,100,-100,100,0,100000,2000,1)
linG3DAliveAll('exampleB005',9,1,-100,100,-100,100,0,100000,2000,1)
for numClone=0:9 
  linG3DClone('exampleB005',numClone,1,-100,100,-100,100,0,100000,2000,1)
end
for numClone=0:9 
  linG3DAliveClone('exampleB005',numClone,1,-100,100,-100,100,0,100000,2000,1)
end


% Example 2:
linG3DAll('exampleB05',147,1,-100,100,-100,100,0,100000,2000,1)
linG3DAliveAll('exampleB05',147,1,-100,100,-100,100,0,100000,2000,1)
for numClone=0:147 
  linG3DClone('exampleB05',numclone,1,-100,100,-100,100,0,100000,2000,1)
end
for numClone=0:147 
  linG3DAliveClone('exampleB05',numclone,1,-100,100,-100,100,0,100000,2000,1)
end


% Example 3
linG3DAll('exampleExp',10,0,0,1500,0,1000,0,864,4,1)
linG3DAliveAll('exampleExp',10,0,0,1500,0,1000,0,864,4,1)
for numClone=1:10 
  linG3DAliveClone('exampleExp',numClone,0,0,1500,0,1000,0,864,4,1)
end
for numClone=1:10
  linG3DClone('exampleExp',numClone,0,0,1500,0,1000,0,864,4,1)
end



