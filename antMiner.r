
getDataRangeMatrix <- function(featureSet, bins=4) {
	dataRangeMatrix <- NULL
	m <- ncol(featureSet)
	for (i in 1:m) {
		dataRangeMatrix <- cbind(dataRangeMatrix, quantile( featureSet[,i], probs = seq(0, 1, 1/bins) ) )
	}
	dataRangeMatrix[nrow(dataRangeMatrix),] <- Inf
	dataRangeMatrix
}

getInfoT <- function(dataSet, dataRangeMatrix) {
	rows <- nrow(dataRangeMatrix) - 1
	cols <- ncol(dataRangeMatrix)
	
	examples <- nrow(dataSet)
	classes <- unique(dataSet[,ncol(dataSet)])
	classNum <- length(classes)
	
	infoT <- matrix(0,rows, cols)
	for (i in 1:rows) {
		for (j in 1:cols) {
			for (k in classes) {
				logicVectorCond <-  (dataRangeMatrix[i,j] <= dataSet[,j] ) & ( dataSet[,j]  < dataRangeMatrix[i+1,j] ) 
				logicVectorClass <- (dataRangeMatrix[i,j] <= dataSet[,j] ) & ( dataSet[,j]  < dataRangeMatrix[i+1,j] ) & ( dataSet[,ncol(dataSet)] == k )
				prob <- length(logicVectorClass[logicVectorClass==TRUE]) / length(logicVectorCond[logicVectorCond==TRUE]) 
				if (prob != 0) 
					infoT[i,j] = infoT[i,j] - prob*log(prob,2)
			}
		}
	}
	infoT
}

getIta <- function(infoT, classNum ) {
	Ita <- ( log(classNum, 2) - infoT ) / sum( log(classNum, 2) - infoT ) 
}

updateTau <- function(Tau , fitnessValue, rules, trailPersistance=0.9) {
	cols <- ncol(rules)-1
	Tau <- Tau - (1-trailPersistance)*Tau
	for (i in 1:cols)
		Tau[ rules[1,i] , rules[2,i] ] <- Tau[ rules[1,i] , rules[2,i] ] + trailPersistance*fitnessValue
	Tau
}

fitness <- function(truth, response, pClass) {
	( length(which( ( response==pClass & truth==pClass ) == TRUE ))  - length(which( ( response==pClass & truth!=pClass ) == TRUE )) ) / length(truth)
}

antMiner <- function(dataSet, bins=4, rowLimit = 10 , n_ants=50, exploitProb = 0.8, t_n = 100) {

	bestRuleList <- NULL
	dataRangeMatrices <- NULL

	for (t in 1:t_n) {
	
		features <- ncol(dataSet) - 1
		examples <- nrow(dataSet)
		classNum <- unique(dataSet[,ncol(dataSet)])
	
		dataRangeMatrix <- getDataRangeMatrix(dataSet[,-ncol(dataSet)], bins)
		infoT <- getInfoT(dataSet, dataRangeMatrix)
		Ita <- infoT / sum(infoT)#getIta(infoT, length(classNum))

		Tau <- matrix(1,features,bins)
		ruleList <- NULL
		fitnessList <- NULL
		for ( ant in 1:n_ants ) {
			probMatrix <- Tau*Ita
			rules <- NULL
			for (i in 1:features) {
				if (sum(probMatrix)==0)
					break
				if (runif(1) < exploitProb)
					rules <- cbind( rules, which( probMatrix == max(probMatrix) , arr.ind=TRUE )[1,] )
				else 
					rules <- cbind( rules, which ( probMatrix == sample(probMatrix,1,prob=probMatrix/sum(probMatrix)) , arr.ind=TRUE )[1,] )
				probMatrix[,rules[2,i]] <- 0
			}
			logicVector <- rep(TRUE, examples)
			for (i in 1:ncol(rules))
				logicVector <- logicVector & ( (dataRangeMatrix[ rules[1,i] , rules[2,i] ] <= dataSet[,rules[2,i]] ) &  (dataSet[,rules[2,i]] < dataRangeMatrix [ rules[1,i] + 1 , rules[2,i] ] ) )
			classRule <- Mode( dataSet[which(logicVector==TRUE),features+1] )
			rules <- cbind(rules, c(classRule,classRule))
			response <- rep( (-classRule)^(-1) , examples )
			response[which(logicVector==TRUE)] <- classRule
			fitnessValue <- fitness(dataSet[,features+1], response, classRule)
			Tau <- updateTau(Tau , fitnessValue, rules, trailPersistance=0.9)
			fitnessList <- c(fitnessList, fitnessValue)
			ruleList[[ant]] <- rules
		}
		bestRuleList[[t]] <- ruleList[[ which(fitnessList == max(fitnessList))[1] ]]
		dataRangeMatrices[[t]] <- dataRangeMatrix
		logicVector <- rep(TRUE, examples)
		for ( i in 1:(ncol(bestRuleList[[t]])-1) )
			logicVector <- logicVector & ( (dataRangeMatrix[ bestRuleList[[t]] [1,i] , bestRuleList[[t]] [2,i] ] <= dataSet[,bestRuleList[[t]] [2,i]] ) &  (dataSet[,bestRuleList[[t]] [2,i]] < dataRangeMatrix [ bestRuleList[[t]] [1,i] + 1 , bestRuleList[[t]] [2,i] ] ) )
		dataSet <- dataSet[-which(logicVector==TRUE),]
		if(nrow(dataSet) < rowLimit)
			break
	}

	for ( i in 1:(t-1) ) {
		cat('Rule ',i,':\n')
		rule <- bestRuleList[[i]]
		for (j in 1:(ncol(rule)-1)) {
			cat('(',rule[2,j],') : ', '[', dataRangeMatrices[[i]] [ rule[1,j] , rule[2,j] ]  ,',',  dataRangeMatrices[[i]] [ rule[1,j] + 1 , rule[2,j] ]  ,')\n')
		}
		cat('class : ',rule[1,ncol(rule)],'\n\n')
	}
	return (c(dataRangeMatrices , bestRuleList))
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

dataSet <- as.matrix(iris[c(1:50,101:150),-5 ])
classes <- c(rep(1,50), rep(-1,50))
dataSet <- cbind(dataSet, classes)
answer <- antMiner (dataSet, bins=4, rowLimit = 10 , n_ants=50, exploitProb = 0.8, t_n = 100) 
