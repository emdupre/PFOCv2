## A script to extract and plot design scores from PLS to R.
## v2.0, emdupre

###############################################
## This needs to be set before running the script.
## Change to your results filename, your group 
## names to the order you entered them in PLS, 
## and list the order in which you want your
## conditions displayed.

resultMat="PFOC_YoungOld_MeanCentering_NoNull_fMRIresult.mat"
groupLevels=c("Young","Old")
condLevels=c("Past", "Future", "Other", "Control")

###############################################

## If error in readMat, execute following command in MATLAB:
## save('your_PLS_filename_result.mat', '-v7')

## First, load relevant packages
library(ggplot2)
library(wesanderson)
library(R.matlab)

## Read in your dset from a result.mat.
dmat<-readMat(resultMat)
Estimate<-as.data.frame(dmat$boot.result[,,1]$orig.usc)
colnames(Estimate)<-paste("Estimate_LV",1:ncol(Estimate),sep="")

UL<-as.data.frame(dmat$boot.result[,,1]$ulusc)
colnames(UL)<-paste("UL_LV",1:ncol(UL),sep="")

LL<-as.data.frame(dmat$boot.result[,,1]$llusc)
colnames(LL)<-paste("LL_LV",1:ncol(LL),sep="")

Condition<-as.factor(t(as.data.frame(dmat$cond.name)))
rownames(Condition)<-NULL

numRepeat<-ncol(Estimate)/ncol(dmat$subj.group)
Group <- as.factor(c(rep(groupLevels[1],numRepeat),rep(groupLevels[2],numRepeat)))

## Create new data structure with this information.
## I've called mine PFOC
PFOC<-data.frame(Estimate,LL,UL,Condition,Group)

## Set factor levels for user preferred display
PFOC$Condition <- factor(PFOC$Condition,levels = condLevels)
PFOC$Group <- factor(PFOC$Group,levels = groupLevels)

## Determine significant LVs for plotting
Significance<-rowSums(as.data.frame(dmat$perm.result[,,1]$s.prob)<0.05)

estimates<-c()
lbs<-c()
ubs<-c()

for(i in 1:length(Significance)){
  if (Significance[i]=="1")
  {estimates<-c(estimates,paste("Estimate_LV",i,sep=""))
    lbs<-c(lbs,paste("LL_LV",i,sep=""))
    ubs<-c(ubs,paste("UL_LV",i,sep=""))}
}

toPlot<-data.frame(cbind(estimates,lbs,ubs))

## Finally, plot significant LVs with ggplot.

## Here, I’ve chosen the Royal1 color scheme. 
## Other options and examples available on the ‘wesanderson’ github.
for(i in 1:nrow(toPlot)){
  print(ggplot(PFOC, aes(x=Condition, y=PFOC[[as.character(toPlot[i,1])]], fill=Group)) + 
          geom_bar(width=.75,position=position_dodge(), stat="identity",
                   size=.2) +
          geom_errorbar(aes(ymin=PFOC[[as.character(toPlot[i,2])]], 
                            ymax=PFOC[[as.character(toPlot[i,3])]]),
                        width=.1,
                        position=position_dodge(.75),
                        colour="black") +
          theme_minimal(base_size = 28, base_family = "Arial") +
          theme(axis.text.y = element_blank()) +
          theme(axis.title.y = element_blank()) +
          theme(axis.title.x = element_text(margin = margin(t= 22))) +
          scale_fill_manual(values=wes_palette("Royal1"))) 
}
