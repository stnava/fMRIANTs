
## ----, echo=TRUE,message=FALSE,cache=FALSE-------------------------------
library(ANTsR)
library(RKRNS)
basedir<-"/Users/stnava/data/fMRIANTs/" # FIXME for your study
setwd(basedir)
motionAcc<-0 # motion accuracy - 0 is for testing, 1 or 2 real studies
gtemplate<-antsImageRead("data2/template2mil.nii.gz",3)
t1<-antsImageRead("data2/t1_2mil.nii.gz",3)
ipttrn<-glob2rx(paste("*boldTask.nii.gz",sep=''))
fns<- paste(basedir,list.files(path=basedir,
  pattern = ipttrn ,recursive=T),sep='/')
designfn<-paste(basedir,"data2/designmat.csv",sep='')
if ( all(dim(gtemplate)==1) | ! file.exists(designfn)
     | all(dim(t1)==1) |  !file.exists(fns[1]) ) 
  stop(paste("Check your working directory",basedir))


## ----, echo=TRUE,cache=FALSE---------------------------------------------
prefix<-paste(basedir,"TEMP",sep='') # paste(tempfile())
figpre<-paste(basedir,"ANTsRfMRI_FAQ_files/figure-html/fig",sep='')
figpre<-paste(basedir,"TEMP_fig",sep='')


## ----, echo=TRUE,message=FALSE,cache=FALSE-------------------------------
if ( ! exists("templateMaps") )
  {
  templateMaps<-antsRegistration(gtemplate,t1,
    typeofTransform = "SyN",outprefix = prefix)
  antsImageWrite(templateMaps$warpedmovout,
    paste(prefix,'t1warped.nii.gz',sep=''))
  }


## ----, echo=TRUE,cache=FALSE---------------------------------------------
fn<-fns[1]
if ( ! file.exists(fn) )
  {
  fn<-file.choose()
  }
bold<-antsImageRead(fn,4)


## ----getavg, echo=TRUE,warning=F,message=F,cache=FALSE-------------------
avg<-getAverageOfTimeSeries(bold)


## ----, echo=F,warning=F,message=F,cache=FALSE----------------------------
onm=paste(figpre,'slices.png',sep='')
plotANTsImage(avg,slices='18x26x4',axis=3,outname=onm)


## ----, echo=TRUE,warning=F,message=F,cache=FALSE-------------------------
N3BiasFieldCorrection(3,avg,avg,2)
N3BiasFieldCorrection(3,avg,avg,2)
mask<-getMask(avg,mean(avg),Inf,2)
mat<-timeseries2matrix(bold, mask)


## ----, echo=TRUE,warning=F,message=F,cache=FALSE-------------------------
dvars<-computeDVARS(mat)
plot(ts(dvars))
# detect outliers?


## ----, echo=TRUE,warning=F,message=F,cache=FALSE-------------------------
  boldarr<-as.array(bold)
  bold2d<-as.antsImage(t(boldarr[20,20,,]))
  onm2da=paste(figpre,'slices2da.png',sep='')
  plotANTsImage(bold2d,outname=onm2da)


## ----, echo=TRUE,warning=F,message=F,cache=FALSE-------------------------
  tr<-antsGetSpacing(bold)[4]
  nbold<-dim(bold)[4]


## ----, eval=TRUE, cache=FALSE--------------------------------------------
  designmat<-read.csv(designfn)


## ----, echo=F,warning=F,message=F,eval=TRUE,cache=FALSE------------------
image(data.matrix(designmat))


## ----, echo=TRUE,warning=F,message=F,eval=TRUE,cache=FALSE---------------
  boldpre<-antsPreprocessfMRI( bold, mask, residualizeMatrix=T,
      spatialSmoothingType='gaussian', useMotionCorrectedImage=T,
      motionCorrectionAccuracyLevel = motionAcc,
      spatialSmoothingParameters=0, numberOfCompCorComponents=6 )


## ----, echo=TRUE,warning=F,message=F,cache=FALSE-------------------------
  boldarr<-as.array(boldpre$cleanBoldImage)
  bold2d<-as.antsImage(t(boldarr[20,20,,]))
  onm2db=paste(figpre,'slices2db.png',sep='')
  plotANTsImage(bold2d,outname=onm2db)


## ----, echo=TRUE,warning=F,message=F,eval=TRUE,cache=FALSE---------------
plot(ts(boldpre$DVARSpost))


## ----, echo=TRUE,warning=F,message=F,eval=TRUE,cache=FALSE---------------
plot(ts(boldpre$FD))


## ----, echo=TRUE,warning=F,message=F,eval=TRUE,cache=FALSE---------------
plot(ts(boldpre$nuisanceVariables[,"compcorr1"]))


## ----, echo=TRUE,warning=F,message=F,eval=TRUE,cache=FALSE---------------
plot(ts(boldpre$nuisanceVariables[,"compcorr2"]))


## ----, echo=TRUE,warning=F,message=F,eval=TRUE,cache=FALSE---------------
  print(colnames(boldpre$nuisanceVariables))
  nuis<-cbind(boldpre$nuisanceVariables,
    DVARS=boldpre$DVARS,FD=boldpre$FD)


## ----, echo=TRUE,warning=F,message=F,eval=TRUE,cache=FALSE---------------
  bold<-antsImageClone( boldpre$cleanBoldImage )
  SmoothImage(4,bold,3.0,bold)
  mat<-timeseries2matrix( bold, mask  )


## ----, echo=TRUE,warning=F,message=F,eval=TRUE,cache=FALSE---------------
hrf<-ts( hemodynamicRF(scans=20, onsets=0, durations=12, rt=tr,cc=0.5,a1=8,a2=9))
plot(hrf)


## ----, echo=FALSE,eval=TRUE,cache=FALSE----------------------------------
 hrfdesignmat <- designmat
 for (i in 1:ncol(hrfdesignmat)) {
   hrfdesignmat[, i] <- conv(hrfdesignmat[, i], hrf)[1:nrow(hrfdesignmat)]
   }
plot( ts(hrfdesignmat[, 1]) )


## ----, echo=TRUE,eval=TRUE,cache=FALSE-----------------------------------
  mdl<-lm( mat ~ A+B+C+D+nuis , data=hrfdesignmat )
  execbetas<-bigLMStats( mdl, 0.0001 )
  print( rownames( execbetas$beta.t ) )
# NOTE - coefficient differences
  betas<-makeImage( mask, 
    execbetas$beta["C",]-execbetas$beta["A",] )


## ----, echo=TRUE,eval=TRUE,cache=FALSE-----------------------------------
  antsImageWrite(betas,paste(prefix,'tbetas.nii.gz',sep=''))
  antsImageWrite(avg,paste(prefix,'avg.nii.gz',sep=''))
  write.csv( nuis, paste(prefix,'nuis.csv',sep=''),row.names=F)


## ----, echo=TRUE,eval=TRUE,cache=FALSE-----------------------------------
  if (!exists("t1brain")) t1brain<-antsImageClone(t1)
  if (!exists("mytx"))
    {
    mytx<-antsRegistration(fixed=avg, moving=t1brain,
      typeofTransform="SyNBold", 
      outprefix=paste(prefix,"B",sep=''))
    }


## ----, echo=TRUE,eval=TRUE,cache=FALSE-----------------------------------
  wnm<-templateMaps$fwdtransforms[[1]]
  anm<-templateMaps$fwdtransforms[[2]]
    fulltx<-c(wnm,anm,mytx$invtransforms[1],
      mytx$invtransforms[2])
    wbeta<-antsApplyTransforms(fixed=gtemplate,
      moving=avg, transformlist=fulltx,
      whichtoinvert=c(F,F,T,F) )
    wbetafn<-paste(prefix,'avgw.nii.gz',sep='')
    antsImageWrite(wbeta,wbetafn)


## ----, echo=TRUE,eval=TRUE,cache=FALSE-----------------------------------
    wbeta<-antsApplyTransforms(fixed=gtemplate,
      moving=betas, transformlist=fulltx,
      whichtoinvert=c(F,F,T,F) )
    wbetafn<-paste(prefix,'tbetasw.nii.gz',sep='')
    antsImageWrite(wbeta,wbetafn)


## ----, echo=TRUE---------------------------------------------------------
# see antsr documentation for how this might be done 
# on the warped beta images


## ----, echo=T,warning=F,message=F,cache=FALSE----------------------------
onm=paste(figpre,'betas.png',sep='')
th<-paste("3.0x",as.character(max(wbeta)),sep='')
plotANTsImage(gtemplate,functional=list(wbeta), thresh=th,
  color='jet',slices='48x62x4',axis=3,outname=onm)


## ----, echo=T------------------------------------------------------------
hist( betas[ mask  == 1  ]   )


