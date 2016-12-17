#' Gini Decomposition Function using Bhattacharya and Mahalanobis (1967) Method
#'
#' This Function allows you to decompose Gini coefficient
#' @param param is the segments
#' @keywords Gini Inequality Decomposition
#' @export
#' @examples Ginidecomp(income, gender, param=2)
#' Ginidecomp()






calcGini<-function(x, y=NULL, group=NULL, param=2)
{
  if (param<=0) return(NULL)
  if (!is.numeric(x)) return(NULL)
  xNA<-sum(as.numeric(is.na(x)))
  if (!is.null(y))
  {
    if(!is.numeric(y)) return(NULL)
    yNA<-sum(as.numeric(is.na(y)))
  }
  weighted<-FALSE
  groupNA<-NULL
  if (is.null(group)) group<-rep(1, length(x))
  else 
  {
    if (!is.numeric(group)) return(NULL)
    weighted<-TRUE
    groupNA<-sum(as.numeric(is.na(group)))
  }
  if (is.null(y)) df<-cbind("x"=x, "group"=group)
  else df<-cbind("x"=x, "y"=y, "group"=group)
  df<-df[complete.cases(df),, drop=FALSE]
  if (nrow(df)==0) return (NULL)
  if (any(df[,"x"]<0)) return(NULL)
  if (sum(df[,"x"])==0) return(NULL)
  index<-0
  names(param)<-"param"
  if (nrow(df)>1)
  {
    if (param != 1)
    {
      df[,"group"]<-df[,"group"]/sum(df[,"group"])
      xMean<-weighted.mean(df[,"x"],df[,"group"])
      if (is.null(y)) df<-df[order(df[,"x"]),]
      else df<-df[order(df[,"y"]),]
      sp<-cumsum(df[,"group"])
      sp[length(sp)]<-1
      sm<-c(0, sp[-length(sp)])
      groupr<-(((1-sm)^param)-((1-sp)^param))/(param*df[,"group"])
      groupvar<-sum(df[,"group"] * ((groupr - weighted.mean(groupr, df[,"group"]))^2))
      reg1<-coef(lm((-param*groupvar*df[,"x"]/xMean)~groupr, weights=df[,"group"]))
      names(reg1)<-NULL
      index<-reg1[2]
    }
  }
  if (is.null(y)) names(index)<-"SGini"
  else names(index)<-"SConc"
  SG<-list(ineq=   list(index=index,
                        parameter=param),
           nas=    NULL)
  if (is.null(y)) SG[["nas"]]<-list(xNA=xNA, groupNA=groupNA,
                                    totalNA=length(x)-nrow(df))
  else SG[["nas"]]<-list(xNA=xNA, yNA=yNA, groupNA=groupNA,
                          totalNA=length(x)-nrow(df))
  class(SG)<-"ICI"
  return(SG)
}






Ginidecomp<-function(x, z, group=NULL, param=2)
{
  if (param<=0) return(NULL)
  if (!is.numeric(x)) return(NULL)
  xNA<-sum(as.numeric(is.na(x)))
  if (!is.factor(z)) return(NULL)
  zNA<-sum(as.numeric(is.na(z)))
  weighted<-FALSE
  groupNA<-NULL
  if (is.null(group)) group<-rep(1, length(x))
  else 
  {
    if (!is.numeric(group)) return(NULL)
    weighted<-TRUE
    groupNA<-sum(as.numeric(is.na(group)))
  }
  df<-data.frame("x"=x, "z"=z, "group"=group)
  df<-df[complete.cases(df),, drop=FALSE]
  if (nrow(df)==0) return (NULL)
  if (any(df[,"x"]<0)) return(NULL)
  if (sum(df[,"x"])==0) return(NULL)
  if (nrow(df)==1) return(NULL)
  names(param)<-"param"
  lx<-length(x)
  df[, "z"]<-factor(df[,"z"], exclude=NULL)
  df[, "group"]<-df[, "group"]/sum(df[, "group"])
  df<-df[order(df[,"x"]),]
  dfSplit<-split(df[,c("x","group")], df[,"z"])
  xMean<-weighted.mean(df[,"x"],df[,"group"])
  xMeanbetween<-sapply(dfSplit,
                  function(df) weighted.mean(df[,"x"],df[,"group"]), simplify=TRUE)
  popweightbetween<-sapply(dfSplit, function(df) sum(df[,"group"]), simplify=TRUE)
  incomesharebetween<-popweightbetween*xMeanbetween/xMean
  if (weighted)
  {
    SGinibetween<-sapply(dfSplit,
          function(df) calcGini(df[,"x"],df[,"group"], param)[["ineq"]][["index"]],
          simplify=FALSE)
    SGini<-calcGini(df[,"x"],df[,"group"], param)[["ineq"]][["index"]]
  }
  else
  {
    SGinibetween<-sapply(dfSplit,
            function(df) calcGini(df[,"x"],NULL, param)[["ineq"]][["index"]],
            simplify=FALSE)
    SGini<-calcGini(df[,"x"],NULL, param)[["ineq"]][["index"]]
  }
  SGinibetween<-unlist(SGinibetween)
  names(SGinibetween)<-names(popweightbetween)
  names(SGini)<-"SGini"
  
  


SGinibetweenshare<-popweightbetween*incomesharebetween*SGinibetween
    SGinibetweenTot<-sum(SGinibetweenshare)
    SGiniwithin<- calcGini(xMeanbetween, popweightbetween, param)[["ineq"]][["index"]]
    SGinioverlap<- SGini - SGinibetweenTot - SGiniwithin
    names(SGiniwithin)<-NULL
    names(SGinioverlap)<-NULL
    SGD<-list(Overall=SGini,
              Within=SGinibetweenTot,
              Between=SGiniwithin,
              Overlap=SGinioverlap,
              GiniGroups=SGinibetween,
              ShareGiniGroups=SGinibetweenshare,
              Popweightbetween=popweightbetween,
              Incomesharebetween=incomesharebetween)

return(SGD)
}