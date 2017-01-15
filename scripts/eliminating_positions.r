args <- commandArgs(T)

#Returns for each query position whether it was found in the beacon list of positions or not

# beaconfile="E:\\Privacy Code Works\\beacon_distribute\\beacon_distribute\\beacon.1000"
beaconfile <- args[[1]]

beacontable = read.table(beaconfile)

beaconpos = beacontable[,1]


# infiles="E:\\Privacy Code Works\\beacon_distribute\\beacon_distribute\\In1.altalleles"
bias = as.numeric(args[[2]])
print(bias)
infiles <- args[3:length(args)]
outfileaccuracy=NULL
j=1
bias=rep(bias, times = length(infiles))
false_positive = rep(0, times = length(infiles))
true_positive = rep(0, times = length(infiles))
false_negetive = rep(0, times = length(infiles))
true_negetive = rep(0, times = length(infiles))
for (file in infiles) {
  
  print(file)
  outfile = gsub("altalleles","rand_queryresults",file)
  outfileaccuracy = gsub("altalleles",paste(bias * 100,"accuracy",sep = ""),file)
  qpos = scan(file)
  qres = 0 + (qpos %in% beaconpos)
  for (i in 1:length(qres)) {
    if (qres[i] == 1) {
      true_positive[j]=true_positive[j]+1;
    }else{
      true_negetive[j]=true_negetive[j]+1;
    }
    
    if (runif(1) > bias) {

        if (qres[i] == 1) {
          qres[i] = 0
          false_negetive[j] = false_negetive[j] + 1
        }else{
          qres[i] = 1
          false_positive[j] = false_positive[j] + 1
        }
      }
    
  }
  cat("false_pos ",false_positive[j]," false_neg ",false_negetive[j],"\n")
  
  results = data.frame(qpos = qpos,qres = qres)
  write.table(
    results,file = outfile,quote = F,row.names = F,col.names = F
  )
  
  j=j+1
}

accuracy = data.frame(bias = bias,false_negetive = false_negetive,false_positive = false_positive,true_positive=true_positive,true_negetive=true_negetive)
tmp = paste('/import/helium-share/csgrad/azizmma/bustamante/accuracy_elim.',toString(bias[1] * 100),sep = "")
write.table(
  accuracy,
  file = tmp,quote = F,row.names = F,col.names = F
)
# rm(list = ls())