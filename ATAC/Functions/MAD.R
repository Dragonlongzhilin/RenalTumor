#### Unsymmetric Distributions and the Double MAD
#ref: https://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/

DoubleMAD <- function(x, zero.mad.action="warn"){
   # The zero.mad.action determines the action in the event of an MAD of zero.
   # Possible values: "stop", "warn", "na" and "warn and na".
   x         <- x[!is.na(x)]
   m         <- median(x)
   abs.dev   <- abs(x - m)
   left.mad  <- median(abs.dev[x<=m])
   right.mad <- median(abs.dev[x>=m])
   if (left.mad == 0 || right.mad == 0){
      if (zero.mad.action == "stop") stop("MAD is 0")
      if (zero.mad.action %in% c("warn", "warn and na")) warning("MAD is 0")
      if (zero.mad.action %in% c(  "na", "warn and na")){
         if (left.mad  == 0) left.mad  <- NA
         if (right.mad == 0) right.mad <- NA
      }
   }
   return(c(left.mad, right.mad))
}

DoubleMADsFromMedian <- function(x, zero.mad.action="warn"){
   # The zero.mad.action determines the action in the event of an MAD of zero.
   # Possible values: "stop", "warn", "na" and "warn and na".
   two.sided.mad <- DoubleMAD(x, zero.mad.action)
   m <- median(x, na.rm=TRUE)
   x.mad <- rep(two.sided.mad[1], length(x))
   x.mad[x > m] <- two.sided.mad[2]
   mad.distance <- abs(x - m) / x.mad
   mad.distance[x==m] <- 0
   return(mad.distance)
}