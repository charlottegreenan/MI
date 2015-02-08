
#### Loading the sienaMI and siena07_MI functions

    source("sienaMI.r")

#### Reading data

 	friend.data.w1 <- as.matrix(read.table("s50-network1.dat"))
	friend.data.w2 <- as.matrix(read.table("s50-network2.dat"))
	friend.data.w3 <- as.matrix(read.table("s50-network3.dat"))
	drink <- as.matrix(read.table("s50-alcohol.dat"))
   	smoke <- as.matrix(read.table("s50-smoke.dat"))
	n <- dim(friend.data.w1)[1]

#### Making smoking binary

    smoke[smoke==1] <- 0
	smoke[smoke %in% c(2,3)] <- 1

#### Giving smoking some missing values
	smoke[c(11,29,24,46,43,7,2,14,34,5),1] <- NA
	smoke[c(39,9,12,8,43,45,2,31,33,35),2] <- NA
	smoke[c(26,31,25,30,19,28,43,36,6,23),3] <- NA

#### Loading RSiena: the package 'arm' must be installed

	library(RSiena)

#### Specifying the data, model and effects

	friendship <- sienaNet(
                     array( c( friend.data.w1, friend.data.w2,
						friend.data.w3), dim = c( n, n, 3 ) ) )

	smokeBeh <- sienaNet(smoke, type='behavior')
	drinkCov <- varCovar(drink)

	mydata <- sienaDataCreate(friendship, smokeBeh, drinkCov)
	myeff <- getEffects(mydata)
	myModel <- sienaModelCreate(useStdInits=FALSE, projname='mi_example')

	myeff[myeff$effectName=='transitive triplets' &
 		myeff$type=='eval', 'include'] <- TRUE
	myeff[myeff$effectName=='smokeBeh similarity' &
 		myeff$type=='eval', 'include'] <- TRUE
	myeff[myeff$effectName=='effect drinkCov on rate smokeBeh' &
 		myeff$type=='rate', 'include'] <- TRUE
    myeff[myeff$effectName=='average exposure effect on rate smokeBeh' &
        myeff$type=='rate', 'include'] <- TRUE

#### Getting the imputations; here we want to use diffusion model and impute
#### 'uponly' data, so we set 'uponly = TRUE'

   	imp <- sienaMI(mydata, myeff, myModel, 'smokeBeh', uponly = TRUE)

#### For example, we may decide we wanted to impute using drinking,
#### but we don't want drinking in the siena analysis model, so we needed it in
#### the effects in sienaMI to impute, but not in siena07_MI
#### to analyse our model:

    myeff[myeff$effectName=='effect drinkCov on rate smokeBeh' &
        myeff$type=='rate', 'include'] <- FALSE

#### Running siena07 for each set of imputed data, and combining the results

    mi_ans <- siena07_MI(imp, myModel, mydata, myeff)

#### Returns a list of the results for each imputed dataset

	mi_ans[[2]]

#### Returns the overall estimates

	mi_ans[[1]]

#### Run again to get better convergence

	mi_ans2 <- siena07_MI(imp, myModel, mydata, myeff, prevAns = mi_ans)
