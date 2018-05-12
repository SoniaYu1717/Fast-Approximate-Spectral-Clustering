library(kernlab);
library(MASS);
library(gtools);
library(cluster);
library(wskm);
# check the result
cRate=function(sp0, sp1, nc, N)
{
    tr = 0;
    seqs=seq(1,nc);
    spx=matrix(0,N,1);
    spy=matrix(0,N,1);
    
    if(nc <8)
    { perms=permutations(n=nc,r=nc);
        np=dim(perms)[1];
        for(i in 1:np)
        {
            for(j in 1:nc) { spx[sp1==j]=perms[i,j]; }
            tmp=sum(sp0==spx)/N;
            if(tr<tmp) {tr=tmp;spy=spx;}
        }
    } else {
        for(i in 1:10000)
        { permx=sample(seqs,nc, replace=FALSE);
            for(j in 1:nc) { spx[sp1==j]=permx[j]; }
            tmp=sum(as.integer(sp0)==spx)/N; if(tr<tmp) {tr=tmp;spy=spx;}
        }
    }
    cat(tr,"\t");
}



espectral=function(x,sp,ncluster,dim, alpha)
{
    
    
    N=nrow(x);  #The # of observations
    m=dim;      #The # of features
    sp0=sp;
    y=x;
    
    
    n=floor(length(y[,1])/alpha);
    
    n1=N*0.05;
    idx= sample((1:N),n1,replace = FALSE);
    xx=y[idx,];
    
    ptm = proc.time()

    # # use coarse k-means
    # cat("begin coarse K-means\n");
    # xxkms=kmeans(xx[,1:m],centers = ncluster, iter.max = 200, nstart = 20,algorithm = c("Hartigan-Wong"));
    # cat("begin K-means\n")
    # xkms= kmeans(y[,1:m],centers = xxkms$centers, iter.max = 200, nstart = 1, algorithm = c("Hartigan-Wong"));
    # spp = xkms$cluster;
    # cRate(sp0,spp,ncluster,N); 
    # cat("run time:", proc.time() - ptm,"\n");
    # ptm = proc.time();



    
    # # # use coarse k-means + spec
    # if(N>20000)
    # {
    #     if(N<300000) {n1<-N*0.1;} else { n1<-N*0.05;}
    #     idx<-sample((1:N),n1,replace=FALSE);
    #     xx<-y[idx,];
    #     cat("Start of coarse K-means ",date(),"\n");
    #     xxkms<-kmeans(xx[,1:m],centers=n, iter.max = 200, nstart = 20,
    #         algorithm = c("Hartigan-Wong"));

    #     cat("Start of K-means @ ",date(),"\n");
    #     xkms<-kmeans(y[,1:m],centers=xxkms$centers, iter.max = 200, nstart = 1,
    #         algorithm = c("Hartigan-Wong"));
    # }
    # else
    # {
    #     cat("Start of K-means @ ",date(),"\n");
    #     xkms<-kmeans(y[,1:m],centers=n, iter.max = 200, nstart = 20,
    #         algorithm = c("Hartigan-Wong"));
    # }
    # tmp = xkms$cluster;
    # x = xkms$centers;
    # sp = specc(x, centers=ncluster);
    # sp = sp@.Data;
    # spp=xkms$cluster;
    # for(i in 1:n)
    # {
    #   spp[spp==i]=sp[i];
    # }
    # cRate(sp0,spp,ncluster,N); 
    # cat("run time:", proc.time() - ptm,"\n");
    # ptm = proc.time();



    
    # # # use ewkm only
    # cat("begin ewkm \n");
    # xkms = ewkm(y[,1:m],ncluster, lambda = 0.5, maxiter = 100);
    # spp = xkms$cluster;

    # cRate(sp0,spp,ncluster,N); 
    # cat("run time:", proc.time() - ptm,"\n");
    # ptm = proc.time();

    # # ewkm + spec
    # cat("ewkm + spec \n");


    for (i in 3:10) {
        n = 100 * i;
        cat(n,"\t")
        xkms = ewkm(y[,1:m],n, lambda = 0.6, maxiter = 200);
        tmp = xkms$cluster;
        x = xkms$centers;
        # cat("spec\n")
        sp = specc(x, centers=ncluster);
        sp = sp@.Data;
        spp=xkms$cluster;
        for(i in 1:n)
        {
          spp[spp==i]=sp[i];
        }
        cRate(sp0,spp,ncluster,N); 
        cat((proc.time() - ptm),"\n");
        ptm = proc.time();
    }
    # # use Clustering Large Applications
    # cat("begin pam \n");
    # xkms = clara(y[,1:m],ncluster , sample = 200, pamLike = TRUE);
    # spp = xkms$cluster;
    
    # cRate(sp0,spp,ncluster,N); 
    # cat("run time:", proc.time() - ptm,"\n");
    # ptm = proc.time();
    

    # Clustering Large Applications + spec
    # cat("begin pam + spec\n");
    # xkms = clara(y[,1:m],n, pamLike = TRUE);
    # tmp = xkms$clustering;
    # x = xkms$medoids;
    # cat("spec\n")
    # sp = specc(x, centers=ncluster);
    # sp = sp@.Data;
    # spp=xkms$cluster;
    # for(i in 1:n)
    # {
    #     spp[spp==i]=sp[i];
    # }
    # cRate(sp0,spp,ncluster,N);
    # cat("run time:", proc.time() - ptm,"\n");
    # ptm = proc.time();
}



# loader
main=function(){
    x=read.table(file="connect4.Rdata",header=FALSE,sep=",",strip.white=FALSE);
    x$label=x$V43;
    sp=matrix(0,nrow(x),1);
    sp[x$label=="win"]="1";
    sp[x$label=="loss"]="2";
    sp[x$label=="draw"]="3";
    tmp=c("V43");
    rmcls=match(tmp,names(x));
    x=x[,-rmcls];
    for(i in 1:42) {x[,i]=as.integer(x[,i])};
    nc=3;
    m=42;
    
    # Normalization
    for(i in 1:m)
    {
        if(sd(x[,i])==0) {x[,i]=0;}
        else {x[,i]=(x[,i]-mean(x[,i]))/sd(x[,i]); }
    }
    # x: data+ label, sp: result, ncluster: 3, 42 feature
    
    z=espectral(x,sp,ncluster=nc,42,6000);
    # z=espectral(x,sp,ncluster=nc,42,3000);
    # z=espectral(x,sp,ncluster=nc,42,2000);
    # z=espectral(x,sp,ncluster=nc,42,1000);
    # z=espectral(x,sp,ncluster=nc,42,600);
    # z=espectral(x,sp,ncluster=nc,42,300);
    # z=espectral(x,sp,ncluster=nc,42,200);


    # cat("* UCI Magic Gamma (19020)\n");
    # x<-read.table(file="mgamma.Rdata",header=FALSE,sep=",",strip.white=FALSE);
    # x$label<-x$V11;
    # sp<-matrix(0,nrow(x),1);
    # sp[x$label=="g"]<-"1";
    # sp[x$label=="h"]<-"2";

    # tmp<-c("V11");
    # rmcls<-match(tmp,names(x));
    # x<-x[,-rmcls];
    # nc<-2;
    # cat("*  # Classes =",nc,"\n");
    # z<-espectral(x,sp,ncluster=nc,10,40);
    # cat("* UCI Musk (6598)\n");
    # x<-read.table(file="musk.Rdata",header=FALSE,sep=",",strip.white=FALSE);

    # x$label<-x$V169;
    # sp<-matrix(0,nrow(x),1);
    # sp[x$label=="1"]<-"1";
    # sp[x$label=="0"]<-"2";

    # tmp<-c("V1","V2","V169");
    # rmcls<-match(tmp,names(x));
    # x<-x[,-rmcls];

    # # To create the prerpocessed data in format usable by MATLAB
    # if(F){
    #     write.table(x,file="muskr",sep=" ",row.names=FALSE,col.names=FALSE);
    # }

    # nc<-2;
    # cat("*  # Classes=",nc,"\n");
    # z<-espectral(x,sp,ncluster=nc,166, 15);

}

main()
