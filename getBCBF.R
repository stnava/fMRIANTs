#!/home/avants/R_builds/R-3.1.1/bin/Rscript
# ! /usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if ( length(args) < 3 )
  {
  cat(" args[[1]] - name of asl file \n")
  cat(" args[[2]] - act t1 directory \n")
  cat(" args[[3]] - output directory & file prefix: e.g asl_study/subject_date_CBF \n")
  cat(" args[[4]] - optional - roi image from template \n")
  q("no")
  }
library(ANTsR)
if ( length(args) > 2 ) {
    outfn<-args[[3]]
  } else {
    outfn<-tempfile()
  }
dirchar<-"/"
outDir<-unlist(strsplit(outfn,dirchar))
myid<-outDir[length(outDir)]
outDir<-paste(outDir[1:(length(outDir)-1) ],
  collapse=dirchar)
# print(paste("Make",outDir))
dir.create(file.path(outDir), showWarnings = FALSE,
  recursive=TRUE )
bcbfn<-paste(outfn,'_bcbf.nii.gz',sep='')
rcbfn<-paste(outfn,'_roicbf.nii.gz',sep='')
rcbfnraw<-paste(outfn,'_roiraw.csv',sep='')
rcbfncsv<-paste(outfn,'_roicbf.csv',sep='')
rvolcsv<-paste(outfn,'_roivol.csv',sep='')
rctcsv<-paste(outfn,'_roict.csv',sep='')
if ( length(args) > 0 ) {
  aslfn<-args[[1]]
} else {
  asldir<-"/Users/stnava/Downloads/PNC_ASL/"
  aslfn<-paste(asldir,
    '005008_ep2d_se_pcasl_PHC_1200ms_SEQ03.nii.gz',sep='')
}
if ( length(args) > 1 ) {
  tdir<-args[[2]]
  tfn<-Sys.glob(paste(tdir,"*BrainSegmentation0N4.nii.gz",sep=''))
  tsfn<-Sys.glob(paste(tdir,"*BrainSegmentation.nii.gz",sep=''))
  templateimg<-antsImageRead(tfn,3)
  templateseg<-antsImageRead(tsfn,3)
  if ( all(dim(templateimg)==1) ) stop(paste("missing",tfn))
  if ( all(dim(templateseg)==1) ) stop(paste("missing",tsfn))
  tlist<-list()
  for ( i in 1:6 )
    {
    popfn<-paste("*BrainSegmentationPosteriors00",i,".nii.gz",sep='')
    pfn<-Sys.glob(paste(tdir,popfn,sep=''))
    pimg<-antsImageRead(pfn,3)
    if ( all(dim(pimg)==1) ) stop(paste("missing",pfn))
    tlist[[i]]<-pimg
    }
  } else {
    tdir<-"/Users/stnava/data/fMRIANTs/data3/"
    tfn<-paste(tdir,'PCASL.nii.gz',sep='')
    templateimg<-antsImageRead(tfn,4)
    templateimg<-getAverageOfTimeSeries(templateimg)
    N3BiasFieldCorrection(3,templateimg,templateimg,2)
    N3BiasFieldCorrection(3,templateimg,templateimg,2)
    maskavg<-getMask(templateimg,mean(templateimg),Inf,2)
    templateimg[maskavg==0]<-0
    tsfn<-paste(tdir,'seg2pcasl.nii.gz',sep='')
    templateseg<-antsImageRead(tsfn,3)
    tlist<-list()
    for ( i in 1:6 )
     {
     pfn<-paste(tdir,"seg2pcaslprob",i,".nii.gz",sep='')
     pimg<-antsImageRead(pfn,3)
     tlist[[i]]<-pimg
     }
  }
# now get subject ASL data
pcasl<-antsImageRead(aslfn,4)
avg<-getAverageOfTimeSeries(pcasl)
N3BiasFieldCorrection(3,avg,avg,2)
N3BiasFieldCorrection(3,avg,avg,2)
maskavg<-getMask(avg,mean(avg),Inf,2)
avg[maskavg==0]<-0
removeLeadAndTrail<-TRUE
if ( removeLeadAndTrail )
  {
  tmat<-timeseries2matrix(pcasl,maskavg)
  tmat<-tmat[5:(nrow(tmat)-4),]
  pcasl<-matrix2timeseries(pcasl,maskavg,tmat)
  }
preg<-antsRegistration( avg, templateimg, typeofTransform='SyNBold',
  outprefix=outfn )
wseg<-antsApplyTransforms( avg, templateseg,
   preg$fwdtransforms,
   "NearestNeighbor")
problist<-list()
for ( i in 1:6 )
  {
  wimg<-antsApplyTransforms( avg,
    tlist[[i]], preg$fwdtransforms )
  problist[[i]]<-wimg
  }
if ( ! file.exists(bcbfn) | TRUE )
  {
  bcbf<-bayesianCBF( pcasl, wseg, problist,
      myPriorStrength=30.0,
      useDataDrivenMask=3,
      denoisingComponents=1:12,
      robustnessvalue=0.925 )
  # reduce robval b/c we are already rejecting a lot via removeLeadAndTrail
  nz<-bcbf$nz
  bcbf<-bcbf$bcbf
  antsImageWrite( bcbf,bcbfn )
  } else {  nz=(-1); bcbf<-antsImageRead(bcbfn,3) }
bcbfsummary<-data.frame( id=myid, nz=nz,
  csfsum=sum(problist[[1]][ problist[[1]] > 0.5 ]),
  gmsum=sum(problist[[2]][ problist[[2]] > 0.5 ]),
  wmsum=sum(problist[[3]][ problist[[3]] > 0.5 ] ),
  dgsum=sum(problist[[4]][ problist[[4]] > 0.5 ] ) ,
  csfcbf=mean(bcbf[ problist[[1]] > 0.5 ]),
  gmcbf=mean(bcbf[ problist[[2]] > 0.5 ]),
  wmcbf=mean(bcbf[ problist[[3]] > 0.5 ] ),
  dgcbf=mean(bcbf[ problist[[4]] > 0.5 ] )
  )
#########
write.csv(bcbfsummary,paste(outfn,'_bcbf.csv',sep=''),row.names=F)
if ( length(args) > 3 ) {
  # map ROIs from the template to subject
  roifn<-args[[4]]
  data("DesikanKillianyTourville",package='ANTsR')
  data("aal",package='ANTsR')
  myRoiDescriptor<-DesikanKillianyTourville
  if ( length(grep("AAL",roifn)) == 1 ) myRoiDescriptor<-aal
  nrois<-nrow(myRoiDescriptor)
  myroidfvol<-matrix(rep(rep(0,nrois),1),nrow=1)
  myroidfvol<-data.frame(myroidfvol)
  myroidfcbf<-data.frame(myroidfvol)
  myroidfct<-data.frame(myroidfvol)
  colnames(myroidfvol)<-myRoiDescriptor$label_name
  colnames(myroidfcbf)<-myRoiDescriptor$label_name
  colnames(myroidfct)<-myRoiDescriptor$label_name
  twfn<-Sys.glob(paste(tdir,'*SubjectToTemplate1Warp.nii.gz',sep=''))
  tafn<-Sys.glob(paste(tdir,'*T1SubjectToTemplate0GenericAffine.mat',sep=''))
  swfn<-Sys.glob(paste(tdir,'*TemplateToSubject0Warp.nii.gz',sep=''))
  safn<-Sys.glob(paste(tdir,'*TemplateToSubject1GenericAffine.mat',sep=''))
  ftxlist<-c( unlist(preg$fwdtransforms),safn,swfn)
  wroi<-antsApplyTransforms( avg, antsImageRead(roifn,3),
    ftxlist,
    "NearestNeighbor")
  # map ROIs to cortex
  wroi[  problist[[2]] < 0.25 ] <- 0
  antsImageWrite( wroi, rcbfn )
  # do some gm roi stuff then use LabelStats to write a csv ....
  # reorg csv in DKT (or AAL) format
  ImageMath(3,rcbfnraw,"LabelStats",wroi,bcbf)
  spc<-antsGetSpacing(wroi)
  volelt<-spc[1]*spc[2]*spc[3]
  roivals<-read.csv(rcbfnraw)
  colct<-1
  for ( anat in colnames(myroidfcbf) )
    {
    num<-as.numeric(myRoiDescriptor$label_num[
      myRoiDescriptor$label_name==anat])
      if ( sum(roivals$label == num ) == 1 )
        {
        # FIXME - someone else should check this!!!
        ctval<-roivals$count[ roivals$label == num  ]
        myvol<-ctval*volelt
        mycbf<-roivals$mass[ roivals$label == num  ]/ctval
        myroidfvol[1,colct ]<-myvol
        myroidfcbf[1, colct]<-mycbf
        myroidfct[1, colct]<-ctval
        }
      colct<-colct+1
    }
  write.csv( myroidfvol , rvolcsv )
  write.csv( myroidfcbf , rcbfncsv )
  write.csv( myroidfct  , rctcsv )
}
