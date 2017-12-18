# Rigbe G. Weldatsadik
# 10.12.2017
# rigbe.weldatsadik@helsinki.fi
# Two protugese schools' Students Performance of Math and Portugese language subjects (https://archive.ics.uci.edu/ml/datasets/Student+Performance)

library(dplyr)

# reading the datasets
math <- read.table("data/student-mat.csv", sep = ";" , header=TRUE)
por <- read.table("data/student-por.csv", sep = ";" , header=TRUE)

# exploring
# the math dataset has 395 observation and 33 variables(of factor and int types), while the por dataset has the same 33 variables but 649 observations
str(math)
str(por)
dim(math)
dim(por)

# I will merge the maths and protugese datasets based on the common students that took both classes and based on every variable other than the 3 grades and the variables'studytime','schoolsupyes','famsupyes','paidclass' which might be specific to the two subjects (For instance, I would imagine there would be few instances of extra paid classes for the portugese language than mathematics)

math_por <- merge(math,por,by=c("school","sex","age","address","famsize","Pstatus",
                             "Medu","Fedu","Mjob","Fjob","reason","nursery","internet",
                             "guardian","guardian","traveltime","failures",
                             "activities","higher","romantic",
                             "famrel","freetime","goout","Dalc","Walc","health","absences"),suffixes = c(".math",".por"))


# So after the merging there are 85 students whose chosen info matches in both datasets and the 3 grades and the non-merged variables will be suffixed with math and por to indicate which data set the values have come from
dim(math_por)
str(math_por)
glimpse(math_por)


# saving the merged data
write.table(math_por,file = "data/math_por.csv", sep = "\t")

# I have done some more data wrangling for the LDA part, particularly scaling and creating a binary variable from the continuous grade3 values
library(corrplot)
library(GGally)
library(ggplot2)

# since in the LDA the explanatory variables should be continous, I will select only the variable whose class is not factor
math_por2 <- math_por[sapply(math_por,class)!= "factor"]

summary(math_por2)

# I am checking the correlation between the different numeric variables so that I will be able to remove those collinear variable later when fitting the LDA
cor_matrix <- cor(math_por2) 

corrplot(cor_matrix, method="circle",type="upper",tl.pos = "d",tl.cex = 0.8)

# then i am scaling the data which is necessary before the LDA analysis
mathPor2_scaled <- as.data.frame(scale(math_por2))
str(mathPor2_scaled)

# from the scaled data I am creating a binary variable for the math and portugese final grades such that values above 13 (or in the scaled version above 0.65) are 'pass' and those below are 'fail'
G3_binary.por = ifelse(mathPor2_scaled$G3.por >= 0.65,'pass','fail')
G3_binary.math = ifelse(mathPor2_scaled$G3.math >= 0.65,'pass','fail')

# then i am deleting the original grades and replacing them with the binary variables
mathPor2_scaled <- mathPor2_scaled %>% dplyr::select(-one_of('G3.por','G3.math'))

mathPor2_scaled <- data.frame(mathPor2_scaled, G3_binary.por,G3_binary.math)

# finally i am writing the result in to a table
write.table(mathPor2_scaled,file = "data/mathPor2_scaled.csv", sep = "\t")

# this a is a function that will calculate VIF values iteratively for lm fitting
# vif function (taken from this website https://beckmw.wordpress.com/2013/02/05/collinearity-and-stepwise-vif-selection/)

vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  library(fmsb)
  
  if(any(!'data.frame' %in% class(in_frame))) in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}

# Because we will be using the following dataset in the lm fitting, I am using that to calculate the VIFs to come up with the predictors with low VIFs. These form.in variables i then use in my lm fitting
math2 <- math_por %>% dplyr::select(-one_of('G1.math','G2.math','paid.por','G1.por','G2.por','G3.por','studytime.por','schoolsup.por','famsup.por'))

keep.dat<-vif_func(in_frame=math2,thresh=3,trace=F)
form.in<-paste('y ~',paste(keep.dat,collapse='+'))

por2 <- math_por %>% dplyr::select(-one_of('G1.por','G2.por','paid.math','G1.math','G2.math','G3.math','studytime.math','schoolsup.math','famsup.math'))
keep.dat<-vif_func(in_frame=por2,thresh=3,trace=F)
form.in.por<-paste('y ~',paste(keep.dat,collapse='+'))
