#OUwie likelihood calculator

#written by Jeremy M. Beaulieu

#Allows the user to calculate the likelihood given a specified set of parameter values.

OUwie.fixed<-function(phy, data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"), simmap.tree=FALSE, root.age=NULL, scaleHeight=FALSE, root.station=FALSE, get.root.theta=FALSE, shift.point=0.5, alpha=NULL, sigma.sq=NULL, theta=NULL, clade=NULL, mserr="none", check.identify=TRUE, quiet=FALSE){
    
    if(model=="BMS" & root.station==TRUE){
        warning("By setting root.station=TRUE, you have specified the group means model of Thomas et al. 2006", call.=FALSE, immediate.=TRUE)
        get.root.theta = FALSE
    }

    if(is.factor(data[,3])==TRUE){
        stop("Check the format of the data column. It's reading as a factor.", .call=FALSE)
    }
    
    if(is.null(root.age)){
        if(any(branching.times(phy)<0)){
            stop("Looks like your tree is producing negative branching times. Must input known root age of tree.", .call=FALSE)
        }
    }
    
    if(model == "OUMV" |  model == "OUMA" |  model == "OUMVA"){
        if(root.station == TRUE){
            stop("Assuming stationarity at the root is no longer allowed under these models. Try OU1 and OUM.", call. = FALSE)
        }
    }

    if(model == "BM1" |  model == "BMS" | model == "TrendyM" | model == "TrendyMS"){
        if(root.station == FALSE){
            get.root.theta = TRUE
        }
    }

    if(check.identify == TRUE){
        identifiable <- check.identify(phy=phy, data=data, simmap.tree=simmap.tree, quiet=TRUE)
        if(identifiable == 0){
            warning("The supplied regime painting is unidentifiable for the regime painting. All regimes form connected subtrees.", call. = FALSE, immediate.=TRUE)
        }
    }

    #Makes sure the data is in the same order as the tip labels
    if(mserr=="none" | mserr=="est"){
        data<-data.frame(data[,2], data[,3], row.names=data[,1])
        data<-data[phy$tip.label,]
    }
    if(mserr=="known"){
        if(!dim(data)[2]==4){
            stop("You specified measurement error should be incorporated, but this information is missing.", call. = FALSE)
        }
        else{
            if(is.factor(data[,4]) == TRUE){
                stop("Check the format of the measurement error column. It's reading as a factor.",  call. = FALSE)
            }else{
                data<-data.frame(data[,2], data[,3], data[,4], row.names=data[,1])
                data<-data[phy$tip.label,]
            }
        }
    }
    
    #Values to be used throughout
    n=max(phy$edge[,1])
    ntips=length(phy$tip.label)
    #Will label the clade of interest if the user so chooses:
    if(is.null(clade)){
        phy=phy
    }
    if(!is.null(clade) & simmap.tree==FALSE){
        node<-mrca(phy)[clade[1],clade[2]]
        int<-c(node,Descendants(phy,node, "all"))
        tips<-int[int<ntips]
        data[,1]<-1
        data[tips,1]<-2
        int<-int[int>ntips]
        phy$node.label<-rep(1,phy$Nnode)
        pp<-int-length(phy$tip.label)
        phy$node.label[pp]<-2
    }
    if (is.character(model)) {
        if (model == "BM1"| model == "OU1"){
            simmap.tree=FALSE
        }
    }
    if(simmap.tree==TRUE){
        k<-length(colnames(phy$mapped.edge))
        tot.states<-factor(colnames(phy$mapped.edge))
        tip.states<-factor(data[,1])
        data[,1]<-as.numeric(tip.states)
        
        #A boolean for whether the root theta should be estimated -- default is that it should be.
        root.station=root.station
        #Obtains the state at the root
        root.state=which(colnames(phy$mapped.edge)==names(phy$maps[[1]][1]))
        ##Begins the construction of the edges matrix -- similar to the ouch format##
        edges=cbind(c(1:(n-1)),phy$edge,MakeAgeTable(phy, root.age=root.age))
        
        if(scaleHeight==TRUE){
            edges[,4:5]<-edges[,4:5]/max(MakeAgeTable(phy, root.age=root.age))
            root.age = 1
        }
        
        edges=edges[sort.list(edges[,3]),]
        
        #Resort the edge matrix so that it looks like the original matrix order
        edges=edges[sort.list(edges[,1]),]
    }
    if(simmap.tree==FALSE){
        #Obtain a a list of all the regime states. This is a solution for instances when tip states and
        #the internal nodes are not of equal length:
        tot.states<-factor(c(phy$node.label,as.character(data[,1])))
        k<-length(levels(tot.states))
        int.states<-factor(phy$node.label)
        phy$node.label=as.numeric(int.states)
        tip.states<-factor(data[,1])
        data[,1]<-as.numeric(tip.states)
 
        #A boolean for whether the root theta should be estimated -- default is that it should be.
        root.station=root.station
        if (is.character(model)) {
            if (model == "BM1"| model == "OU1"){
                ##Begins the construction of the edges matrix -- similar to the ouch format##
                #Makes a vector of absolute times in proportion of the total length of the tree
                k=length(levels(tip.states))
                phy$node.label<-sample(c(1:k),phy$Nnode, replace=T)
                int.states=length(levels(tip.states))
                #Since we only really have one global regime, make up the internal nodes -- this could be improved
                phy$node.label<-as.numeric(int.states)
                #Obtain root state -- for both models assume the root state to be 1 since no other state is used even if provided in the tree
                root.state <- 1
                #New tree matrix to be used for subsetting regimes
                edges=cbind(c(1:(n-1)),phy$edge,MakeAgeTable(phy, root.age=root.age))
                
                if(scaleHeight == TRUE){
                    edges[,4:5] <-  edges[,4:5]/max(MakeAgeTable(phy, root.age=root.age))
                    root.age = 1
                }
                
                edges=edges[sort.list(edges[,3]),]
                
                regime <- matrix(0,nrow=length(edges[,1]),ncol=k)
                regime[,1]<-1
                if (k>1) {
                    regime[,2:k]<-0
                }
                
                edges=cbind(edges,regime)
            }
            else{
                #Obtain root state and internal node labels
                root.state<-phy$node.label[1]
                int.state<-phy$node.label[-1]
                #New tree matrix to be used for subsetting regimes
                edges=cbind(c(1:(n-1)),phy$edge,MakeAgeTable(phy, root.age=root.age))
                if(scaleHeight==TRUE){
                    edges[,4:5]<-edges[,4:5]/max(MakeAgeTable(phy, root.age=root.age))
                    root.age = 1
                }
                edges=edges[sort.list(edges[,3]),]
                
                mm<-c(data[,1],int.state)
                regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
                #Generates an indicator matrix from the regime vector
                for (i in 1:length(mm)) {
                    regime[i,mm[i]] <- 1
                }
                #Finishes the edges matrix
                edges=cbind(edges,regime)
            }
        }
        #Resort the edge matrix so that it looks like the original matrix order
        edges=edges[sort.list(edges[,1]),]
    }
    x<-as.matrix(data[,2])
    
    #Matches the model with the appropriate parameter matrix structure
    if (is.character(model)) {
        Rate.mat <- matrix(1, 2, k)
        
        if (model == "BM1"){
            np <- 1
            index <- matrix(TRUE,2,k)
            Rate.mat[1,1:k] <- 1e-10
            Rate.mat[2,1:k] <- sigma.sq
            param.count <- np+1
            bool <- TRUE
        }
        if (model == "BMS"){
            np <- k
            index<-matrix(TRUE,2,k)
            Rate.mat[1,1:k] <- 1e-10
            Rate.mat[2,1:k] <- sigma.sq
            if(root.station == TRUE){
                param.count <- np+k
            }
            if(root.station == FALSE){
                param.count <- np+1
            }
            bool <- root.station
        }
        if (model == "OU1"){
            np <- 2
            index <- matrix(TRUE,2,k)
            Rate.mat[1,1:k] <- alpha
            Rate.mat[2,1:k] <- sigma.sq
            param.count <- np+1
            bool <- root.station
        }
        if (model == "OUM"){
            np <- 2
            index <- matrix(TRUE,2,k)
            Rate.mat[1,1:k] <- alpha
            Rate.mat[2,1:k] <- sigma.sq
            param.count <- np+k
            bool <- root.station
        }
        if (model == "OUMV") {
            np <- k+1
            index <- matrix(TRUE,2,k)
            Rate.mat[1,1:k] <- alpha
            Rate.mat[2,1:k] <- sigma.sq
            param.count <- np+k
            bool <- root.station
        }
        if (model == "OUMA") {
            np <- k+1
            index <- matrix(TRUE,2,k)
            Rate.mat[1,1:k] <- alpha
            Rate.mat[2,1:k] <- sigma.sq
            param.count <- np+k
            bool <- root.station
        }
        if (model == "OUMVA") {
            np <- k*2
            index <- matrix(TRUE,2,k)
            Rate.mat[1,1:k] <- alpha
            Rate.mat[2,1:k] <- sigma.sq
            param.count <- np+k
            bool <- root.station
        }
    }
    
    obj<-NULL
    
    #Likelihood function for estimating model parameters
    dev.fixed <- function(){
        N <- length(x[,1])
        V <- varcov.ou(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=bool, shift.point=shift.point)
        if(get.root.theta == TRUE){
            W <- weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=FALSE, shift.point=shift.point)
        }else{
            W <- weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=TRUE, shift.point=shift.point)
        }
        if(mserr=="known"){
            diag(V)<-diag(V)+(data[,3]^2)
        }
        if(is.null(theta)){
            theta<-pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)%*%t(W)%*%pseudoinverse(V)%*%x
            se<-sqrt(diag(pseudoinverse(t(W)%*%pseudoinverse(V)%*%W)))
        }
        else{
            theta=theta
            se=rep(NA,length(theta))
        }
        
        theta.est <- cbind(theta,se)
        res <- W%*%theta-x
        #When the model includes alpha, the values of V can get too small, the modulus does not seem correct and the loglik becomes unstable. This is one solution:
        DET <- sum(log(abs(Re(diag(qr(V)$qr)))))
        #However, sometimes this fails (not sure yet why) so I just toggle between this and another approach:
        if(!is.finite(DET)){
            DET <- determinant(V, logarithm=TRUE)
            logl <- -.5*(t(x-W%*%theta)%*%pseudoinverse(V)%*%(x-W%*%theta))-.5*as.numeric(DET$modulus)-.5*(N*log(2*pi))
        }else{
            logl <- -.5*(t(x-W%*%theta)%*%pseudoinverse(V)%*%(x-W%*%theta))-.5*as.numeric(DET)-.5*(N*log(2*pi))
        }
        list(-logl,theta.est,res)
    }
    if(quiet==FALSE){
        cat("Calculating likelihood using fixed parameter values:",c(alpha,sigma.sq,theta), "\n")
    }
    
    fixed.fit <- dev.fixed()
    loglik<- -fixed.fit[[1]]
    
    if(get.root.theta == TRUE){
        W <- weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=FALSE, shift.point=shift.point)
    }else{
        W <- weight.mat(phy, edges, Rate.mat, root.state=root.state, simmap.tree=simmap.tree, root.age=root.age, scaleHeight=scaleHeight, assume.station=TRUE, shift.point=shift.point)
    }
    
    if(get.root.theta == TRUE){
        rates <- cbind(NA, Rate.mat, NA)
        thetas <- cbind(t(fixed.fit[[2]][,1]), NA)
        rates <- rbind(rates, thetas)
        weights <- cbind(W, data[,1])
        regime.weights <- rbind(rates, weights)
        rownames(regime.weights) <- c("alpha", "sigma.sq", "theta", phy$tip.label)
        colnames(regime.weights) <- c("Root", levels(tot.states), "Tip_regime")
    }else{
        rates <- cbind(Rate.mat, NA)
        thetas <- cbind(t(fixed.fit[[2]][,1]), NA)
        rates <- rbind(rates, thetas)
        weights <- cbind(W, data[,1])
        regime.weights <- rbind(rates, weights)
        rownames(regime.weights) <- c("alpha", "sigma.sq", "theta", phy$tip.label)
        colnames(regime.weights) <- c(levels(tot.states), "Tip_regime")
    }
    
    obj = list(loglik = loglik, AIC = -2*loglik+2*param.count, AICc=-2*loglik+(2*param.count*(ntips/(ntips-param.count-1))), BIC=-2*loglik + log(ntips) * param.count, model=model, param.count=param.count, solution=Rate.mat, theta=fixed.fit[[2]], tot.states=tot.states, simmap.tree=simmap.tree, root.age=root.age, shift.point=shift.point, data=data, phy=phy, root.station=root.station, scaleHeight=scaleHeight, res=fixed.fit[[3]], get.root.theta=get.root.theta, regime.weights=regime.weights)
    class(obj)<-"OUwie.fixed"
    return(obj)
}


print.OUwie.fixed<-function(x, ...){
    
    ntips=Ntip(x$phy)
    output<-data.frame(x$loglik,x$AIC,x$AICc,x$BIC,x$model,ntips, row.names="")
    names(output)<-c("lnL","AIC","AICc","BIC","model","ntax")
    cat("\nFit\n")
    print(output)
    cat("\n")
    
    if (is.character(x$model)) {
        if (x$model == "BM1" | x$model == "BMS"){
            param.est <- x$solution
            rownames(param.est)<-c("alpha","sigma.sq")
            if(x$root.station==FALSE){
                theta.mat <- matrix(t(x$theta[1,]), 2, length(levels(x$tot.states)))
            }
            else{
                theta.mat<-matrix(t(x$theta), 2, length(levels(x$tot.states)))
            }
            colnames(theta.mat)<-c("estimate", "se")
            if(x$simmap.tree==FALSE){
                colnames(param.est) <- colnames(theta.mat) <- levels(x$tot.states)
            }
            if(x$simmap.tree==TRUE){
                colnames(param.est) <- colnames(theta.mat) <- c(colnames(x$phy$mapped.edge))
            }
            cat("Rates\n")
            print(param.est)
            cat("\n")
            cat("Optima\n")
            print(theta.mat)
            cat("\n")
        }
        if (x$get.root.theta == FALSE){
            if (x$model == "OU1" | x$model == "OUM"| x$model == "OUMV"| x$model == "OUMA" | x$model == "OUMVA"){
                param.est<- x$solution
                rownames(param.est)<-c("alpha","sigma.sq")
                theta.mat<-matrix(t(x$theta), 2, length(levels(x$tot.states)))
                rownames(theta.mat)<-c("estimate", "se")
                if(x$simmap.tree==FALSE){
                    colnames(param.est) <- colnames(theta.mat)<- levels(x$tot.states)
                }
                if(x$simmap.tree==TRUE){
                    colnames(param.est) <- colnames(theta.mat) <- c(colnames(x$phy$mapped.edge))
                }
                cat("\nRates\n")
                print(param.est)
                cat("\n")
                cat("Optima\n")
                print(theta.mat)
                cat("\n")
            }
        }
        if (x$get.root.theta == TRUE){
            if (x$model == "OU1" | x$model == "OUM"| x$model == "OUMV"| x$model == "OUMA" | x$model == "OUMVA"){
                param.est<- x$solution
                rownames(param.est)<-c("alpha","sigma.sq")
                theta.mat<-matrix(t(x$theta), 2, length(levels(x$tot.states))+1)
                rownames(theta.mat)<-c("estimate", "se")
                if(x$simmap.tree==FALSE){
                    colnames(param.est) <- levels(x$tot.states)
                    colnames(theta.mat)<-c("Root", levels(x$tot.states))
                }
                if(x$simmap.tree==TRUE){
                    colnames(param.est) <- c(colnames(x$phy$mapped.edge))
                    colnames(theta.mat)<-c("Root", colnames(x$phy$mapped.edge))
                }
                cat("\nRates\n")
                print(param.est)
                cat("\n")
                cat("Optima\n")
                print(theta.mat)
                cat("\n")
            }
        }
    }
}

