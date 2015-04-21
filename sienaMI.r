
sienaMI <- function(data, effects, model, behaviorName, noImp = 5,
maxWaveToImpute = NULL, uponly = NULL, downonly = NULL)
{
    depvarsNo <- length(data$depvars)
    behNo <- c(1:depvarsNo)[names(data$depvars) == behaviorName]
    beh <- data$depvars[[behNo]][,1,]
    n <- dim(data$depvars[[behNo]])[1]
    M <- dim(data$depvars[[behNo]])[3]
    imp <- array(rep(NA,n*M*noImp), dim=c(noImp, n, M))

    #### Setting waves to impute
    if(is.null(maxWaveToImpute))
    {
        maxWave <- M
    } else {
        maxWave <- maxWaveToImpute
    }

    wavesToImpute <- c(1:maxWave)

    if (any(!(wavesToImpute %in% c(1:M))) )
    {
        stop("Invalid maxWaveToImpute option.")
    }

    beh <- beh[,wavesToImpute]

    library(arm)

    if (is.null(uponly))
    {
        uponly <- FALSE
    }
    if (is.null(downonly))
    {
        downonly <- FALSE
    }

    maxBeh <- max(beh, na.rm=TRUE)
    minBeh <- min(beh, na.rm=TRUE)
	
	if (maxBeh > 1)
	{
		stop("Not yet implemented for non-binary behavior data.")
	}

    ## Imputing by wave ##
    for (m in 1:noImp)
    {
		cat("\nImputation ");cat(m);cat("\n\n")
        behImp <- beh

        ## Wave 1 ##
        cat("Imputing data for wave 1\n")
        mean.match <- rbinom(sum(is.na(behImp[,1])), maxBeh, prob=mean(behImp[,1],
        na.rm=TRUE)/maxBeh)
        behImp[is.na(behImp[,1]),1] <- mean.match

        ## Waves 2+ ##

        for (w in 2:M)
        {

            if (w <= maxWave)
            {

            cat("Imputing data for wave ");cat(w);cat("\n")

            if (uponly)
            {
                behImp[,w] <- pmax(behImp[,w-1],behImp[,w])
                behImp[behImp[,w-1]==maxBeh,w:maxWave] <- maxBeh
            }

            if (downonly)
            {
                behImp[,w] <- pmin(behImp[,w-1],behImp[,w])
                behImp[behImp[,w-1]==minBeh,w:maxWave] <- minBeh
            }

            ## Non rate effects
            targetsNonRate <- RSiena:::actorTargets(data, effects, model,
            behaviorName, w-1, imputedData = behImp[,w-1])

            targets <- targetsNonRate[[3]]

            if (any(targetsNonRate[[2]] == "rate"))
            {
                ## Rate effects
                targetsRate <- RSiena:::actorTargets(data, effects, model,
                behaviorName, w-1, behImp[,w-1], rate = TRUE)

                targets[,targetsNonRate$effectType=="rate"] <-
                targetsRate[[3]][,targetsNonRate$effectType=="rate"]
            }

            ## At risk means that the actor can make a change

            at.risk <- rep(TRUE, n)

            if (uponly)
            {
                at.risk[behImp[,w-1]==maxBeh] <- FALSE

            } else if (downonly) {

                at.risk[behImp[,w-1]==minBeh] <- FALSE
            }

            at.risk <- at.risk&!is.na(behImp[,w])

                bayes1 <- bayesglm(behImp[at.risk,w] ~ targets[at.risk,],
                family = binomial(link = 'cloglog'))
                post1 <- sim(bayes1, 1)@coef
                targets <- cbind(1, targets)
                p <- 1 - exp(-exp(targets%*%t(post1)))

            behImp[is.na(behImp[,w]),w] <- rbinom(sum(is.na(behImp[,w])),maxBeh,
            p[is.na(behImp[,w])])

            }

        }
        imp[m,,1:maxWave] <- behImp
    }

    effectNames <- targetsNonRate[[1]]
    list(effects=effectNames, imputedData=imp, behaviorName=behaviorName)

}





siena07_MI <- function(imputed_data, model, data, effects, carryForwardTheta =
						TRUE, batch=TRUE, prevAns=NULL,...)
{
behaviorName <- imputed_data[[3]]
imputed_data <- imputed_data[[2]]

n.imp <- dim(imputed_data)[1]
ans <- vector("list", n.imp)

## Checking linear and quadratic effects

linearEffects <- effects[effects$shortName=="linear" & effects$name==
                                                    behaviorName,"include"]
quadEffects <- effects[effects$shortName=="quad" & effects$name==
                                                    behaviorName,"include"]

for (i in 1:n.imp)
{
    beh.data <- imputed_data[i,,]
    imputed.behaviour <- sienaNet(beh.data, type="behavior")
    assign(behaviorName, imputed.behaviour)
    data.objects <- c(names(mydata[[3]]),names(mydata[[4]]),names(mydata[[5]]),
    names(mydata[[6]]),names(mydata[[7]]),names(mydata[[8]]))
    data.arguments.list <- vector("list", length(data.objects))
    names(data.arguments.list) <- data.objects

    for (j in 1:length(data.objects))
    {
        data.arguments.list[[j]] <- eval(parse(text=data.objects[j]))
    }

    newData <- do.call("sienaDataCreate", data.arguments.list)

    ## check linear and quadratic effects

    newEffects <- getEffects(newData)
    newLinearEffects <- newEffects[newEffects$shortName=="linear" &
                                    newEffects$name == behaviorName,"include"]
    newQuadEffects <- newEffects[newEffects$shortName=="quad" &
                                    newEffects$name == behaviorName,"include"]

    if (length(linearEffects)>0 & length(newLinearEffects)>0)
    {
        linearEffects <- linearEffects & newLinearEffects
    } else if (length(linearEffects)>0) {
        linearEffects <- rep(FALSE, 3)
    }
    if (length(quadEffects)>0 & length(newQuadEffects)>0)
    {
        quadEffects <- quadEffects & newQuadEffects
    } else if (length(quadEffects)>0){
        quadEffects <- rep(FALSE, 3)
    }
}

if (length(linearEffects) > 0)
{
	effects[effects$shortName=="linear" &
			effects$name== behaviorName,"include"] <- linearEffects
}
if (length(quadEffects)>0)
{
	effects[effects$shortName=="quad" &
			effects$name== behaviorName,"include"] <- quadEffects
}

## Analysing each of the imputed datasets
analyse <- rep(TRUE, n.imp)
if (!is.null(prevAns))
{
	for (i in 1:n.imp)
	{
		if (max(abs(prevAns[[2]][[i]]$tstat)) < 0.1)
		{
			analyse[i] <- FALSE
		}
	}
}

for (i in c(1:n.imp))
{
	if (analyse[i])
	{
	cat("\n\nAnalysing imputed dataset ");cat(i);cat("\n\n")
	## recreate behaviour
	beh.data <- imputed_data[i,,]
	imputed.behaviour <- sienaNet(beh.data, type="behavior")
	assign(behaviorName, imputed.behaviour)

	## recreate the data object
	data.objects <- c(names(mydata[[3]]),names(mydata[[4]]),names(mydata[[5]]),
				names(mydata[[6]]),names(mydata[[7]]),names(mydata[[8]]))
	data.arguments.list <- vector("list", length(data.objects))
	names(data.arguments.list) <- data.objects

	for (j in 1:length(data.objects))
	{
		data.arguments.list[[j]] <- eval(parse(text=data.objects[j]))
	}
	newData <- do.call("sienaDataCreate", data.arguments.list)

	## recreate the effects
	newEffects <- effects

	## carry forward theta from previous ans
	if (carryForwardTheta & (i>1))
	{
		initialValues <- newEffects$initialValue
		initialValues[newEffects$include] <- ans[[i-1]]$theta
		newEffects$initialValue[newEffects$include & (newEffects$name !=
			behaviorName)] <- initialValues[newEffects$include &
			(newEffects$name != behaviorName)]
	}

	ans[[i]] <- siena07(model, data=newData, effects=newEffects,
					batch = batch,prevAns = prevAns[[2]][[i]],...)
	} else {
	ans[[i]] <- prevAns[[2]][[i]]
	}
}

## Combining the results using Rubin's rules ##

p <- length(ans[[1]]$theta)
thetaMat <- matrix(rep(NA, n.imp*p), ncol = n.imp)
seMat <- matrix(rep(NA, n.imp*p), ncol = n.imp)

notConverged <- 0

for (i in 1:n.imp)
{
	if (max(abs(ans[[i]]$tstat))>0.3)
	{
		notConverged <- notConverged + 1
	}
	thetaMat[,i] <- ans[[i]]$theta
	seMat[,i] <- sqrt(diag(ans[[i]]$cov))
}

if (notConverged > 0)
{
	cat("Problems with convergence for ")
	cat(notConverged)
	cat(" analyses.\n")
}

theta <- apply(thetaMat,1,mean)
se <- (1+1/n.imp)*sqrt(apply(thetaMat,1,var) + apply(seMat^2, 1, mean))

list(combined_results=cbind(effectName=effects$effectName[effects$include],
        Estimate=theta,StandardError=se), ans=ans)
}
