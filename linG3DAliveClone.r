linG3DAliveClone <- function() {
######################################################################
## This is a companion code for the paper "LinG3D: Visualizing the  ##
## Spatio-Temporal Dynamics of Clonal Evolution" by A. Hu, A.M.E.   ##
## Ojwang', K.D. Olumoyin, and K.A. Rejniak                         ##
## This code generates the 3D lineage tree of one clone of number   ##
## specified in 'cloneNum' taking into account only the cells that  ##
## survived to the end of simulation. It uses data from directory   ##
## 'pathdata'.                                                      ##
##                                                                  ##
## The following parameters need to be specified:                   ##
##   pathdata  -- directory with input data                         ##
##   CloneNum  -- clone number to be drawn                          ##
## It requires the following data in the pathdata/data/ directory:  ##
##   cell_history.txt -- file with info about each cell             ##
##   cellID_##.txt    -- cell IDs in a file with index number ##    ## 
##   cellXY_##.txt    -- cell coordinates in a file with index ##   ##
##   drug.txt         -- concentration of a drug for background     ## 
## The following parameters are project-dependent and should be     ##
## specified for a given project:                                   ##
##   xmin,xmax,ymin,ymax -- dimensions of the spacial domain        ##
##   tmin, tmax          -- dimensions of the temporal domain       ##
##   fileStep            -- frequency of the sampled data           ## 
## for the examples discussed in the paper use:                     ##
##   example 1: pathdata='exampleB05';  cloneNum between 0 and 9    ##
##   example 2: pathdata='exampleB005'; cloneNum between 0 and 147  ##
##                                                                  ##
## October 31, 2022                                                 ##
######################################################################
  
library(readr)  # to read text files
library(rapportools)  # for is.empty
library(rgl) # for 3D 


options(scipen = 999)  # disables printing in scientific notation 
pathdata <- 'exampleB005'
#pathdata <- 'exampleB05'
cloneNum <- 0
toPrint <- 1  # save the final figure
IsGradient <- 1  # draw drug gradient in the background 1-yes; 0-no;
dataDirectory <- '/data/' # directory with cell and drug data

xmin <- -100; xmax <- 100;
ymin <- xmin; ymax <- xmax  # 2D domain boundaries
tmin <- 0; tmax <- 100000  # time/iteration boundaries
timeStep <- (tmax-tmin)/(2.5*(xmax-xmin))
fileStep <- 2000  # frequency of data

 if (toPrint==1){
   pathFigs <- paste0(pathdata,'/fig_clonesAlive')
   dir.create(pathFigs)
 }

#################################################################################
#---------------------prepare 3d view using rgl---------------------------------#

open3d()
axes3d()
view3d()
add_axes <- function(x, y, z, axis.col = NA) {   
  minnum <- function(k){c(-max(abs(k)), max(abs(k))) * 1.05}  
  xminnum <- minnum(x);  zminnum <- minnum(z); 
  rgl.lines(xminnum, c(0, 0), c(0, 0), color = NA)
  rgl.lines(c(0, 0), y, c(0, 0), color = NA)
  rgl.lines(c(0, 0), c(0, 0), zminnum, color = NA) } # needed to properly position the 3D
add_axes(c(xmin,xmax),c(0,tmax/timeStep),c(ymin,ymax))  # set the axes
aspect3d(1, 2, 1)  # required aspect ratio x=1, y=2, z=1
grid3d(c("x+-","y+","z-"), col = "#DCDCDC", n = 6)  # set grid on the visible axes only
goodview <- matrix(c(0.6949471,-0.7133674,0.09030721,0,
                     0.1982032,0.3107649,0.92959183,0,
                     -0.6912048,-0.6281178,0.35735679,0,
                     0,0,0,1), nrow = 4,ncol = 4, byrow = TRUE)  # user view captured from screen
highlevel()
par3d(windowRect = c(0,0, 1200, 1200), userMatrix = goodview, zoom=0.75, cex=2)  # set final view

##################################################################################
#---------------------prepare the color palette----------------------------------#

   col <- c("#FF00FF","#FF0000","#00FFFF","#0000FF","#00FF00","#000000","#FFBF00",
            "#FFFF00","#BFFF00","#808000","#FFB6C1","#00BFFF","#0080FF","#FAEBD7",
            "#8000FF","#9ACD32","#FF0080","#660000","#664D00","#006666","#CCCCFF",
            "#FFCCFF","#99CCFF","#FF9999","#009900","#009999","#99004D","#FFE4E1",
            "#800000","#666699","#99FFCC","#DA70D6","#FF8000","#C0C0C0","#808080",
            "#4B0082","#A52A2A","#D8BFD8","#DC143C","#F5DEB3","#FF6347","#FF7F50",
            "#CD5C5C","#F08080","#E9967A","#FA8072","#FFA07A","#FF4500","#FF8C00",
            "#FFA500","#FFD700","#B8860B","#DAA520","#006400","#FFF0F5","#BC8F8F",
            "#FFF8DC","#32CD32","#90EE90","#98FB98","#8FBC8F","#00FA9A","#00FF7F",
            "#2E8B57","#66CDAA","#3CB371","#20B2AA","#2F4F4F","#008080","#008B8B",
            "#F0E68C","#F5F5DC","#E0FFFF","#00CED1","#FFE4B5","#FF1493","#AFEEEE",
            "#7FFFD4","#B0E0E6","#5F9EA0","#4682B4","#6495ED","#DEB887","#1E90FF",
            "#EEE8AA","#BDB76B","#6B8E23","#7CFC00","#7FFF00","#ADFF2F","#B22222",
            "#DDA0DD","#FFEBCD"	) #}
   Ncol <- length(col)
   
#################################################################################
  
DrawBackground <- function(drug,tmax,timeStep,xmin,xmax,ymin,ymax){
     
     drugmin <- min(min(drug)); drugmax <- max(max(drug)); drugstep <- (drugmax-drugmin)/4
     
     kk <- tmax/timeStep
     Nx <- nrow(drug);  Ny <- ncol(drug); hgx <- (xmax-xmin)/Nx; hgy <- (ymax-ymin)/Ny
     
     for (ii in 1:Nx){
       for (jj in 1:Ny){
         if ((drug[ii,jj]>=drugmin)&&(drug[ii,jj]<drugmin+drugstep)){
           polygon3d(x=c(xmin+(ii-1)*hgx,xmin+ii*hgx,xmin+ii*hgx,xmin+(ii-1)*hgx,xmin+(ii-1)*hgx),
                     y=c(kk,kk,kk,kk,kk),
                     z=c(ymin+(jj-1)*hgy,ymin+(jj-1)*hgy,ymin+jj*hgy,ymin+jj*hgy,ymin+(jj-1)*hgy),
                     col = 'blue',fill=TRUE,add = TRUE,plot = TRUE,coords = c(x=1,z=3),alpha=0.25)
           
         } else if ((drug[ii,jj]>=drugmin+drugstep)&&(drug[ii,jj]<drugmin+2*drugstep)) {
           polygon3d(x=c(xmin+(ii-1)*hgx,xmin+ii*hgx,xmin+ii*hgx,xmin+(ii-1)*hgx,xmin+(ii-1)*hgx),
                     y=c(kk,kk,kk,kk,kk),
                     z=c(ymin+(jj-1)*hgy,ymin+(jj-1)*hgy,ymin+jj*hgy,ymin+jj*hgy,ymin+(jj-1)*hgy),
                     col = "#00FFFF",fill=TRUE,add = TRUE,plot = TRUE,coords = c(x=1,z=3),alpha=0.25)
           
         } else if ((drug[ii,jj]>=drugmin+2*drugstep)&&(drug[ii,jj]<drugmin+3*drugstep)){
           polygon3d(x=c(xmin+(ii-1)*hgx,xmin+ii*hgx,xmin+ii*hgx,xmin+(ii-1)*hgx,xmin+(ii-1)*hgx),
                     y=c(kk,kk,kk,kk,kk),
                     z=c(ymin+(jj-1)*hgy,ymin+(jj-1)*hgy,ymin+jj*hgy,ymin+jj*hgy,ymin+(jj-1)*hgy),
                     col = "yellow",fill=TRUE,add = TRUE,plot = TRUE,coords = c(x=1,z=3),alpha=0.25)
           
         } else {
           polygon3d(x=c(xmin+(ii-1)*hgx,xmin+ii*hgx,xmin+ii*hgx,xmin+(ii-1)*hgx,xmin+(ii-1)*hgx),
                     y=c(kk,kk,kk,kk,kk),
                     z=c(ymin+(jj-1)*hgy,ymin+(jj-1)*hgy,ymin+jj*hgy,ymin+jj*hgy,ymin+(jj-1)*hgy),
                     col="red",fill=TRUE,add = TRUE,plot = TRUE,coords = c(x=1,z=3),alpha=0.25)
         }
       }
     }
   }
  
 # draw background with drug gradient
 if (IsGradient==1){
   drug <- read.table(paste0(pathdata,dataDirectory,'drug.txt'),header = F)
   DrawBackground(drug,tmax,timeStep,xmin,xmax,ymin,ymax)
 }

#################################################################################
   
# load cell history file
hist <- read.table(paste0(pathdata,dataDirectory,'cell_history.txt'),header = F)
# [cell ID, clone ID, mother ID, birth iter, div/death iter]

# load indices of all survived cells
cellID <- as.matrix(read.table(paste0(pathdata,dataDirectory,'cellID_',tmax,'.txt'),header = F))

print(paste0('clone = ',cloneNum))

# identify all survived cells from a given clone
indLast <- c(); Nlast <- 0
  for (ii in 1:length(cellID)) { # all cells in the last file
    if (hist[cellID[ii],2]==cloneNum){
      Nlast <- Nlast+1
      indLast[Nlast] <- cellID[ii]
    }
  }
Numlast <- Nlast
  
if (Numlast>0) {  # will take care of empty clones
  
  for (ii in 1:Numlast){
    moNum <- hist[indLast[ii],3]
    while ((moNum>0)&&(hist[moNum,2]==cloneNum)) {
      Nlast <- Nlast+1
      indLast[Nlast] <- moNum
      moNum <- hist[moNum,3]
    }
  }
  indLast <- unique(indLast)  # remove repeated cells

  Nmatrix <- 0
  # run to find Nmatrix
  for (ii in 1:length(indLast)) { # for every cell with index in indLast
    if (ii%%100==0) { print('... calculating'); }
    
    cellNum <- hist[indLast[ii],1]# cell ID
    mothNum <- hist[indLast[ii],3]# mother ID
    strtNum <- hist[indLast[ii],4]# cell birth
    Num  <- hist[indLast[ii],5]
    Num <- max(tmax,min(Num,tmax))# cell div/death/tmax
    
    # find all appearances of the cellNum
    kkStart <- fileStep*floor(strtNum/fileStep)  # initial file number
    kkEnd   <- fileStep*floor(Num/fileStep)  # final file number
    
    i <- kkEnd; vec <- c(); j <- 1
    while (i>=(kkStart+fileStep)) { vec[j] <- i; i <- i-fileStep; j <- j+1;}
    
    for (kk in vec) {
      #for (kk in kkEnd:-fileStep:kkStart+fileStep){# inspect all files
      # cell ID and cell XY from the first file
      fileMeID  <- as.matrix(read.table(paste0(pathdata,dataDirectory,'cellID_',kk,'.txt'),header = F),dimnames=NULL)
      colnames(fileMeID) <- NULL;
      fileMeXY  <- as.matrix(read.table(paste0(pathdata,dataDirectory,'cellXY_',kk,'.txt'),header = F),dimnames=NULL)
      colnames(fileMeXY) <- NULL;
      indMe = which(fileMeID==cellNum)  # which current indices of cellID
      # cell ID and cell XY from the second file
      fileMe2ID <- as.matrix(read.table(paste0(pathdata,dataDirectory,'cellID_',kk-fileStep,'.txt'),header = F),dimnames=NULL)
      colnames(fileMe2ID) <- NULL;
      fileMe2XY <- as.matrix(read.table(paste0(pathdata,dataDirectory,'cellXY_',kk-fileStep,'.txt'),header = F),dimnames=NULL)
      colnames(fileMe2XY) <- NULL;
      indMe2 = which(fileMe2ID==cellNum)  # which current indices of cellID
      
      if (is.empty(indMe)){
      } else if (is.empty(indMe2)){
        while (kkStart<hist[mothNum,4]){# find file with the grand-mother cell
          mothNum <- hist[mothNum,3]
        }
        fileMe2ID <- as.matrix(read.table(paste0(pathdata,dataDirectory,'cellID_',kkStart,'.txt'),header = F),dimnames=NULL)
        colnames(fileMe2ID) <- NULL;
        fileMe2XY <- as.matrix(read.table(paste0(pathdata,dataDirectory,'cellXY_',kkStart,'.txt'),header = F),dimnames=NULL)
        colnames(fileMe2XY) <- NULL;
        indMe2=which(fileMe2ID==mothNum)# which current indices of mother cellID
        
        if (is.empty(indMe2)){
        } else {
          Nmatrix <- Nmatrix+1# save branch to draw [x1,t1,y1,x2,t2,y2]
        }
      } else {
        Nmatrix <- Nmatrix+1# save branch to draw [x1,t1,y1,x2,t2,y2]
      }
    }
  }
  
if (Nmatrix>0){
  # define matrix of line segments (3D branches) to draw
  matrix_to_draw <- matrix(rep(0,Nmatrix*6), nrow = Nmatrix,ncol = 6) # [x1,k1,y1,x2,k2,y2]
  Nmatrix <- 0

  for (ii in 1:length(indLast)) { # for every cell with index in indLast
    if (ii%%100==0) { print('... calculating'); }

    cellNum <- hist[indLast[ii],1]# cell ID
    mothNum <- hist[indLast[ii],3]# mother ID
    strtNum <- hist[indLast[ii],4]# cell birth
    Num  <- hist[indLast[ii],5]
    Num <- max(tmax,min(Num,tmax))# cell div/death/tmax

    # find all appearances of the cellNum
    kkStart <- fileStep*floor(strtNum/fileStep)  # initial file number
    kkEnd   <- fileStep*floor(Num/fileStep)  # final file number
    
    i <- kkEnd; vec <- c(); j <- 1
    while (i>=(kkStart+fileStep)) { vec[j] <- i; i <- i-fileStep; j <- j+1;}
    
    for (kk in vec) {
      # cell ID and cell XY from the first file
      fileMeID  <- as.matrix(read.table(paste0(pathdata,dataDirectory,'cellID_',kk,'.txt'),header = F),dimnames=NULL)
      colnames(fileMeID) <- NULL;
      fileMeXY  <- as.matrix(read.table(paste0(pathdata,dataDirectory,'cellXY_',kk,'.txt'),header = F),dimnames=NULL)
      colnames(fileMeXY) <- NULL;
      indMe = which(fileMeID==cellNum)  # which current indices of cellID
      # cell ID and cell XY from the second file
      fileMe2ID <- as.matrix(read.table(paste0(pathdata,dataDirectory,'cellID_',kk-fileStep,'.txt'),header = F),dimnames=NULL)
      colnames(fileMe2ID) <- NULL;
      fileMe2XY <- as.matrix(read.table(paste0(pathdata,dataDirectory,'cellXY_',kk-fileStep,'.txt'),header = F),dimnames=NULL)
      colnames(fileMe2XY) <- NULL;
      indMe2 = which(fileMe2ID==cellNum)  # which current indices of cellID
      
      if (is.empty(indMe)){
      } else if (is.empty(indMe2)){
        while (kkStart<hist[mothNum,4]){  # find file with the grand-mother cell
          mothNum <- hist[mothNum,3]
        }
        fileMe2ID <- as.matrix(read.table(paste0(pathdata,dataDirectory,'cellID_',kkStart,'.txt'),header = F),dimnames=NULL)
        colnames(fileMe2ID) <- NULL;
        fileMe2XY <- as.matrix(read.table(paste0(pathdata,dataDirectory,'cellXY_',kkStart,'.txt'),header = F),dimnames=NULL)
        colnames(fileMe2XY) <- NULL;
        indMe2=which(fileMe2ID==mothNum)  # which current indices of mother cellID
        
        if (is.empty(indMe2)){
        } else {
          Nmatrix <- Nmatrix+1  # save branch to draw [x1,t1,y1,x2,t2,y2]
          matrix_to_draw[Nmatrix,1] <- fileMeXY[indMe,1]
          matrix_to_draw[Nmatrix,2] <- kkStart+fileStep
          matrix_to_draw[Nmatrix,3] <- fileMeXY[indMe,2]
          matrix_to_draw[Nmatrix,4] <- fileMe2XY[indMe2,1]
          matrix_to_draw[Nmatrix,5] <- kkStart
          matrix_to_draw[Nmatrix,6] <- fileMe2XY[indMe2,2]
        }
      } else {
        Nmatrix <- Nmatrix+1  # save branch to draw [x1,t1,y1,x2,t2,y2]
        matrix_to_draw[Nmatrix,1] <- fileMe2XY[indMe2,1]
        matrix_to_draw[Nmatrix,2] <- kk-fileStep
        matrix_to_draw[Nmatrix,3] <- fileMe2XY[indMe2,2]
        matrix_to_draw[Nmatrix,4] <- fileMeXY[indMe,1]
        matrix_to_draw[Nmatrix,5] <- kk
        matrix_to_draw[Nmatrix,6] <- fileMeXY[indMe,2]
      }
    }
  }

# drawing clones
NumCol <- (cloneNum%%Ncol)+1 # clone color
} # end if Nmatrix>0

 bgplot3d({ plot.new(); 
   title(paste0(main = 'survived cells from clone = ',cloneNum), line = 0,cex.main=2.5) 
   text(x=0.2,y=0.122,"iterations/time x 200",srt=-41,cex=3) })
 
 for (ii in 1:Nmatrix){
   lines3d(x=c(matrix_to_draw[ii,1],matrix_to_draw[ii,4]),
           y=c(matrix_to_draw[ii,2],matrix_to_draw[ii,5])/timeStep,
           z=c(matrix_to_draw[ii,3],matrix_to_draw[ii,6]),
           col=col[NumCol],lwd = 3,plot = TRUE, add = TRUE,bty="b",
           xlim=c(-100,100),ylim=c(0,tmax/timeStep),zlim=c(-100,100))
   
if (toPrint==1){
     rgl.snapshot(paste0(pathFigs,'/tree_clone_',cloneNum,'.png'), fmt = "png", top = TRUE)}
 }
}  # end if Numlast > 0  
  
print("No survived clones")

}  # end function



