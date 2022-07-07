### In this script, all analysis step from R Markdown documents (part 1- 5) are combined.


### load libraries

library("tidyverse")
library("uwot")
library("factoextra")


### load the data

MS_Table <- read.table("RDeeP_HeLa_Mitosis.csv", header = TRUE, row.names = 1, sep = ";")
head(rownames(MS_Table),12)
head(colnames(MS_Table),12)






### Analysis part 1



## data clean-up

# check whether all data are numeric
sum(apply(MS_Table, 2, is.numeric)) == ncol(MS_Table)
glimpse(MS_Table)

# remove proteins that contain only zeros in at least one replicate
del_row <- c() # save  proteins that contain only zeros in at least one replicate (row number) in del_row

for (i in c("Ctrl_Rep1", "Ctrl_Rep2", "Ctrl_Rep3", "RNase_Rep1", "RNase_Rep2", "RNase_Rep3")) {
  # iterate through all replicates
  
  for (j in 1:nrow(MS_Table)) {
    # iterate through all values
    
    if (sum(MS_Table[j, grep(i, colnames(MS_Table))]) <= 0) {
      del_row <- append(del_row, j)
    }
  }
}
del_row_sort <- sort(unique(del_row))
del_row_sort
MS_Table <- MS_Table[-del_row_sort,] # remove the rows saved in del_row_sort



## normalize replicates

for (i in seq(1, 148, by = 3)) {
  # iterate through the fractions (each fraction has 3 replicates)  
  
  #print(i)
  
  # compute the sums of every replicate
  sums <- c(sum(MS_Table[, i]), sum(MS_Table[, i + 1]), sum(MS_Table[, i + 2]))
  
  # compute the mean value of the two closest sums
  if (abs(sums[1] - sums[2]) < abs(sums[1] - sums[3]) &
      abs(sums[1] - sums[2]) < abs(sums[2] - sums[3])) {
    mean_sum <- mean(c(sums[1], sums[2]))
  } else if (abs(sums[2] - sums[3]) < abs(sums[1] - sums[2]) &
             abs(sums[2] - sums[3]) < abs(sums[1] - sums[3])) {
    mean_sum <- mean(c(sums[2], sums[3]))
  } else {
    mean_sum <- mean(c(sums[1], sums[3]))
  }
  #print(mean_sum)
  
  # compute the normalization factors
  norm_factor <- mean_sum / sums
  #print(norm_factor)
  
  # multiply the normalization factors with their corresponding replicate
  MS_Table[, i] <- MS_Table[, i] * norm_factor[1]
  MS_Table[, i + 1] <- MS_Table[, i + 1] * norm_factor[2]
  MS_Table[, i + 2] <- MS_Table[, i + 2] * norm_factor[3]
}

counter <- 0
# there is the need to use round because the numbers seem to be slightly different behind the comma
print(sum(MS_Table[,7]), digits = 17)
print(sum(MS_Table[,8]), digits = 17)

for (i in seq(1, 148, by = 3)) {
  if (round(sum(MS_Table[, i])) != round(sum(MS_Table[, i + 1])) &
      round(sum(MS_Table[, i])) != round(sum(MS_Table[, i + 2]))) {
    counter <- counter + 1
    print(i)
  }
}
counter == 0



## normalize protein amount per fraction

# create df that contains all rows of control replicate 1
x <- grep("Ctrl_Rep1", colnames(MS_Table))
Ctrl_Rep1 <- MS_Table[, x]

# create df that contains all rows of RNase replicate 1
y <- grep("RNase_Rep1", colnames(MS_Table))
RNase_Rep1 <- MS_Table[, y]

# create df that contains all rows of control replicate 2
z <- grep("Ctrl_Rep2", colnames(MS_Table))
Ctrl_Rep2 <- MS_Table[, z]

# create df that contains all rows of RNase replicate 2
a <- grep("RNase_Rep2", colnames(MS_Table))
RNase_Rep2 <- MS_Table[, a]

# create df that contains all rows of control replicate 3
b <- grep("Ctrl_Rep3", colnames(MS_Table))
Ctrl_Rep3 <- MS_Table[, b]

# create df that contains all rows of RNase replicate 3
c <- grep("RNase_Rep3", colnames(MS_Table))
RNase_Rep3 <- MS_Table[, c]

# normalization step 1 of control replicate 1
for (i in 1:nrow(Ctrl_Rep1)) {
  x <- sum(Ctrl_Rep1[i, ]) / 100
  Ctrl_Rep1[i, ] = Ctrl_Rep1[i, ] / x
}

# normalization step 1 of RNase replicate 1
for (i in 1:nrow(RNase_Rep1)) {
  x <- sum(RNase_Rep1[i, ]) / 100
  RNase_Rep1[i, ] = RNase_Rep1[i, ] / x
}

# normalization step 1 of control replicate 2
for (i in 1:nrow(Ctrl_Rep2)) {
  x <- sum(Ctrl_Rep2[i, ]) / 100
  Ctrl_Rep2[i, ] = Ctrl_Rep2[i, ] / x
}

# normalization step 1 of RNase replicate 2
for (i in 1:nrow(RNase_Rep2)) {
  x <- sum(RNase_Rep2[i, ]) / 100
  RNase_Rep2[i, ] = RNase_Rep2[i, ] / x
}

# normalization step 1 of control replicate 3
for (i in 1:nrow(Ctrl_Rep3)) {
  x <- sum(Ctrl_Rep3[i, ]) / 100
  Ctrl_Rep3[i, ] = Ctrl_Rep3[i, ] / x
}

# normalization step 1 of RNase replicate 3
for (i in 1:nrow(RNase_Rep3)) {
  x <- sum(RNase_Rep3[i, ]) / 100
  RNase_Rep3[i, ] = RNase_Rep3[i, ] / x
}

# test if the sum of each row divided by 100 is the number of rows will tell us whether the sum is 100
sum(Ctrl_Rep1) / 100 == nrow(Ctrl_Rep1)
sum(RNase_Rep1) / 100 == nrow(RNase_Rep1)
sum(Ctrl_Rep2) / 100 == nrow(Ctrl_Rep2)
sum(RNase_Rep2) / 100 == nrow(RNase_Rep2)
sum(Ctrl_Rep3) / 100 == nrow(Ctrl_Rep3)
sum(RNase_Rep3) / 100 == nrow(RNase_Rep3)

df_normalized <- as.data.frame(cbind(Ctrl_Rep1, Ctrl_Rep2, Ctrl_Rep3, RNase_Rep1, RNase_Rep2, RNase_Rep3))
dim(df_normalized) == dim(MS_Table)
sum(df_normalized[1,]) / 6 == 100
df_normalized_copy <- df_normalized



## evaluating reproducibility of our data

tdf_normalized <- as.data.frame(t(df_normalized))

# correlation between ctrl1 and ctrl2 for every protein
c1c2 <- c()
for (i in 1:ncol(tdf_normalized)) {
  c1c2 <- append(c1c2, (cor(tdf_normalized[1:25, i], tdf_normalized[26:50, i])))
}

# correlation between ctrl1 and ctrl3 for every protein
c1c3 <- c()
for (i in 1:ncol(tdf_normalized)) {
  c1c3 <- append(c1c3, (cor(tdf_normalized[1:25, i], tdf_normalized[51:75, i])))
}

# correlation between ctrl2 and ctrl3 for every protein 
c2c3 <- c()
for (i in 1:ncol(tdf_normalized)) {
  c2c3 <- append(c2c3, (cor(tdf_normalized[26:50, i], tdf_normalized[51:75, i])))
}

# correlation between rnase1 and rnase2 for every protein
r1r2 <- c()
for (i in 1:ncol(tdf_normalized)) {
  r1r2 <- append(r1r2, (cor(tdf_normalized[76:100, i], tdf_normalized[101:125, i])))
}

# correlation between rnase1 an rnase3 for every protein
r1r3 <- c()
for (i in 1:ncol(tdf_normalized)) {
  r1r3 <- append(r1r3, (cor(tdf_normalized[76:100, i], tdf_normalized[126:150, i])))
}

# correlation between rnase2 and rnase3 for every protein
r2r3 <- c()
for (i in 1:ncol(tdf_normalized)) {
  r2r3 <- append(r2r3, (cor(tdf_normalized[101:125, i], tdf_normalized[126:150, i])))
}

tres_values <- c()

for (t in seq(0, 0.8, by = 0.05)) {
  tc <- c()
  
  for (i in 1:ncol(tdf_normalized)) {
    if (c1c2[i] < t) {
      tc <- append(tc, colnames(tdf_normalized)[i])
    }
  }
  for (i in 1:ncol(tdf_normalized)) {
    if (c1c3[i] < t) {
      tc <- append(tc, colnames(tdf_normalized)[i])
    }
  }
  for (i in 1:ncol(tdf_normalized)) {
    if (c2c3[i] < t) {
      tc <- append(tc, colnames(tdf_normalized)[i])
    }
  }
  for (i in 1:ncol(tdf_normalized)) {
    if (r1r2[i] < t) {
      tc <- append(tc, colnames(tdf_normalized)[i])
    }
  }
  for (i in 1:ncol(tdf_normalized)) {
    if (r1r3[i] < t) {
      tc <- append(tc, colnames(tdf_normalized)[i])
    }
  }
  for (i in 1:ncol(tdf_normalized)) {
    if (r2r3[i] < t) {
      tc <- append(tc, colnames(tdf_normalized)[i])
    }
  }
  tres_values <- append(tres_values, length(unique(tc)))
}

plot(seq(0, 0.8, 0.05), tres_values, type = "b", xlab = "threshold value", ylab = "number of removed proteins", main = "")

t = 0.6
tc <- c()

for (i in 1:ncol(tdf_normalized)) {
  if (c1c2[i] < t) {
    tc <- append(tc, colnames(tdf_normalized)[i])
  } 
}
for (i in 1:ncol(tdf_normalized)) {
  if (c1c3[i] < t) {
    tc <- append(tc, colnames(tdf_normalized)[i])
  } 
}
for (i in 1:ncol(tdf_normalized)) {
  if (c2c3[i] < t) {
    tc <- append(tc, colnames(tdf_normalized)[i])
  } 
}
for (i in 1:ncol(tdf_normalized)) {
  if (r1r2[i] < t) {
    tc <- append(tc, colnames(tdf_normalized)[i])
  } 
}
for (i in 1:ncol(tdf_normalized)) {
  if (r1r3[i] < t) {
    tc <- append(tc, colnames(tdf_normalized)[i])
  } 
}
for (i in 1:ncol(tdf_normalized)) {
  if (r2r3[i] < t) {
    tc <- append(tc, colnames(tdf_normalized)[i])
  } 
}

length(unique(tc))
utc <- unique(tc)
remove_proteins <- c()

for (i in 1:length(utc)) {
  remove_proteins <- append(remove_proteins, which(rownames(df_normalized) == utc[i]))
}
if (length(remove_proteins) != 0) {
  df_normalized <- df_normalized[-remove_proteins,]
}
dim(df_normalized)



## combining replicates using mean value

df_Ctrl <- df_normalized[, grep("Ctrl", colnames(df_normalized))]
df_RNase <- df_normalized[, grep("RNase", colnames(df_normalized))]

# control condition
df_Ctrl_c <- data.frame()
for (i in 1:nrow(df_Ctrl)) {
  for (j in 1:25) {
    x <- mean(c(df_Ctrl[i, j], df_Ctrl[i, j + 25], df_Ctrl[i, j + 50]))
    df_Ctrl_c[i, j] = x
  }
}

# RNase condition
df_RNase_c <- data.frame()
for (i in 1:nrow(df_RNase)) {
  for (j in 1:25) {
    x <- mean(c(df_RNase[i, j], df_RNase[i, j + 25], df_RNase[i, j + 50]))
    df_RNase_c[i, j] = x
  }
}

# test if mean values in df_Ctrl_c are computed correctly
score1 <- 0
for (i in 1:nrow(df_Ctrl)) {
  for (j in 1:25) {
    if (mean(c(df_Ctrl[i, j], df_Ctrl[i, j + 25], df_Ctrl[i, j + 50])) != df_Ctrl_c[i, j]) {
      score1 <- score1 + 1
    }
  }
}
score1 == 0

# test if mean values in df_RNase_c are computed correctly
score2 <- 0
for (i in 1:nrow(df_RNase)) {
  for (j in 1:25) {
    if (mean(c(df_RNase[i, j], df_RNase[i, j + 25], df_RNase[i, j + 50])) != df_RNase_c[i, j]) {
      score2 <- score2 + 1
    }
  }
}
score2 == 0

df_normalized_combined <- cbind(df_Ctrl_c, df_RNase_c)
rownames(df_normalized_combined) <- rownames(df_normalized)

cnc <- sapply(1:25, function(x) {paste("Fraction", x, "_Ctrl", sep = "")})
cnr <- sapply(1:25, function(x) {paste("Fraction", x, "_RNase", sep = "")})
colnames(df_normalized_combined)[1:25] <- cnc
colnames(df_normalized_combined)[26:50] <- cnr





### Analysis part 2



## global maxima for combined values

d <- dim(df_normalized_combined)
rows <- d[1] 
columns <- d[2]
fractions_ctrl <- 25
fractions_RNase <- 25

df_global_maxima <- data.frame()
proteins <- rownames(df_normalized_combined)
fractionsnames <- colnames(df_normalized_combined)

for (z in 1:rows) {
  proteinname <- proteins[z]
  
  maxfc <- df_normalized_combined[z, 1]
  maxfcc <- 1
  
  for (fc in 2:25) {
    
    if (maxfc < df_normalized_combined[z, fc]) {
      maxfc <- df_normalized_combined[z, fc]
      maxfcc <- fc
    }
  }
  
  maxrna <- df_normalized_combined[z, 26]
  maxfcr <- 26
  
  for (fc in 27:50) {
    
    if (maxrna < df_normalized_combined[z, fc]) {
      maxrna <- df_normalized_combined[z, fc]
      maxfcr <- fc
    }
  }
  
  df_global_maxima <-
    rbind(
      df_global_maxima,
      data.frame(
        #global_maximum_ctrl_fraction1 = fractionsnames[maxfcc],
        ctrl_global_maximum_fraction = maxfcc,
        ctrl_global_maximum_value = maxfc,
        #global_maximum_rnase_fraction1 = fractionsnames[maxfcr],
        rnase_global_maximum_fraction = maxfcr - 25,
        rnase_global_maximum_value = maxrna
      )
    )
}

rownames(df_global_maxima) <- rownames(df_normalized_combined)
#View(df_global_maxima)
dim(df_global_maxima)



## global maxima for every replicate

df_ctrl1 <- df_normalized[,grep("Ctrl_Rep1", colnames(df_normalized))]
df_ctrl2 <- df_normalized[,grep("Ctrl_Rep2", colnames(df_normalized))]
df_ctrl3 <- df_normalized[,grep("Ctrl_Rep3", colnames(df_normalized))]
df_rnase1 <- df_normalized[,grep("RNase_Rep1", colnames(df_normalized))]
df_rnase2 <- df_normalized[,grep("RNase_Rep2", colnames(df_normalized))]
df_rnase3 <- df_normalized[,grep("RNase_Rep3", colnames(df_normalized))]

find_global_maxima <- function(df) {
  
  d <- dim(df)
  rows <- d[1]
  columns <- d[2]
  fractions_ctrl <- 25
  #fractions_RNase <- 25
  
  df_output <- data.frame()
  proteins <- rownames(df)
  fractionsnames <- colnames(df)
  
  for (z in 1:rows) {
    proteinname <- proteins[z]
    
    maxfc <- df[z, 1]
    maxfcc <- 1
    
    for (fc in 2:25) {
      if (maxfc < df[z, fc]) {
        maxfc <- df[z, fc]
        maxfcc <- fc
      }
    }
    
    #maxrna <- df[z, 26]
    #maxfcr <- 26
    
    #for (fc in 27:50) {
    #if (maxrna < df[z, fc]) {
    #maxrna <- df[z, fc]
    #maxfcr <- fc
    #}
    #}
    
    df_output <-
      rbind(
        df_output,
        data.frame(
          #global_maximum_ctrl_fraction1 = fractionsnames[maxfcc],
          ctrl_global_maximum_fraction = maxfcc,
          ctrl_global_maximum_value = maxfc
          #global_maximum_rnase_fraction1 = fractionsnames[maxfcr],
          #rnase_global_maximum_fraction = maxfcr - 25,
          #rnase_global_maximum_value = maxrna
        )
      )
  }
  
  rownames(df_output) <- rownames(df)
  #View(df_global_maxima)
  #dim(df_output)
  
  return(df_output)
}

df_ctrl1_global_maxima <- find_global_maxima(df_ctrl1)
df_ctrl2_global_maxima <- find_global_maxima(df_ctrl2)
df_ctrl3_global_maxima <- find_global_maxima(df_ctrl3)
df_rnase1_global_maxima <- find_global_maxima(df_rnase1)
df_rnase2_global_maxima <- find_global_maxima(df_rnase2)
df_rnase3_global_maxima <- find_global_maxima(df_rnase3)

df_global_maxima_rep <-
  cbind(
    df_ctrl1_global_maxima,
    df_ctrl2_global_maxima,
    df_ctrl3_global_maxima,
    df_rnase1_global_maxima,
    df_rnase2_global_maxima,
    df_rnase3_global_maxima
  )
colnames(df_global_maxima_rep) <-
  c(
    "ctrl1_global_maximum_fraction",
    "ctrl1_global_maximum_value",
    "ctrl2_global_maximum_fraction",
    "ctrl2_global_maximum_value",
    "ctrl3_global_maximum_fraction",
    "ctrl3_global_maximum_value",
    "rnase1_global_maximum_fraction",
    "rnase1_global_maximum_value",
    "rnase2_global_maximum_fraction",
    "rnase2_global_maximum_value",
    "rnase3_global_maximum_fraction",
    "rnase3_global_maximum_value"
  )



## local maxima 

# local maxima ctrl

# create df that will contain the local maxima, because we do not know yet how many local maxima each protein will have we need to make space for more than one per protein and cut the rest of later
df_local_maxima_ctrl <-
  data.frame(
    row.names = rownames(df_normalized_combined),
    ctrl_local_maximum1_fraction = matrix(0, nrow(df_normalized_combined), 1),
    ctrl_local_maximum1_value = matrix(0, nrow(df_normalized_combined), 1),
    ctrl_local_maximum2_fraction = matrix(0, nrow(df_normalized_combined), 1),
    ctrl_local_maximum2_value = matrix(0, nrow(df_normalized_combined), 1),
    ctrl_local_maximum3_fraction = matrix(0, nrow(df_normalized_combined), 1),
    ctrl_local_maximum3_value = matrix(0, nrow(df_normalized_combined), 1),
    ctrl_local_maximum4_fraction = matrix(0, nrow(df_normalized_combined), 1),
    ctrl_local_maximum4_value = matrix(0, nrow(df_normalized_combined), 1),
    ctrl_local_maximum5_fraction = matrix(0, nrow(df_normalized_combined), 1),
    ctrl_local_maximum5_value = matrix(0, nrow(df_normalized_combined), 1)
  )

for (i in 1:nrow(df_normalized_combined)) {
  # first loop iterating through the proteins
  over_4 <- c()
  # store indices of values over 4 for every protein in here (only temporary)
  
  for (j in 1:25) {
    # second loop iterating through all Ctrl fractions
    if (df_normalized_combined[i, j] > 4) {
      # only values > 4 are saved
      over_4 <- append(over_4, j)
    }
  }
  #print(over_4)
  local_maxima <- c()
  # store indices of values over 4 and with no higher neighbours (only temporary)
  
  for (v in over_4) {
    # iterate through all values that are > 4
    #print(df_normalized_combined[i, v])
    
    # values on boarders have to be treated differently
    
    # first value (no value to the left)
    if (v == 1) {
      if (df_normalized_combined[i, v] > df_normalized_combined[i, v + 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v + 2] &
          df_normalized_combined[i, v] != max(df_normalized_combined[i, 1:25])) {
        local_maxima <- append(local_maxima, v)
      }
      
      # second value (only one value to the left)
    } else if (v == 2) {
      if (df_normalized_combined[i, v] > df_normalized_combined[i, v - 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v + 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v + 2] &
          df_normalized_combined[i, v] != max(df_normalized_combined[i, 1:25])) {
        local_maxima <- append(local_maxima, v)
      }
      
      # 24th value (only one value to the right)
    } else if (v == 24) {
      if (df_normalized_combined[i, v] > df_normalized_combined[i, v - 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v - 2] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v + 1] &
          df_normalized_combined[i, v] != max(df_normalized_combined[i, 1:25])) {
        local_maxima <- append(local_maxima, v)
      }
      
      # 25th value (no value to the right)
    } else if (v == 25) {
      if (df_normalized_combined[i, v] > df_normalized_combined[i, v - 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v - 2] &
          df_normalized_combined[i, v] != max(df_normalized_combined[i, 1:25])) {
        local_maxima <- append(local_maxima, v)
      }
      
      # all other values to not have to be treated differently
    } else {
      if (df_normalized_combined[i, v] > df_normalized_combined[i, v - 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v - 2] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v + 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v + 2] &
          df_normalized_combined[i, v] != max(df_normalized_combined[i, 1:25])) {
        local_maxima <- append(local_maxima, v)
      }
    }
  }
  #print(paste("protein", as.character(i)))
  #print(local_maxima)
  
  # now we need to save the protein name as well as the local maximum's position and value
  
  if (length(local_maxima) == 1) {
    # proteins with one local maximum
    
    df_local_maxima_ctrl[i, 1] <- local_maxima[1]
    df_local_maxima_ctrl[i, 2] <- df_normalized_combined[i, local_maxima[1]]
    
  } else if (length(local_maxima) == 2) {
    # proteins with 2 local maxima
    
    df_local_maxima_ctrl[i, 1] <- local_maxima[1]
    df_local_maxima_ctrl[i, 2] <- df_normalized_combined[i, local_maxima[1]]
    df_local_maxima_ctrl[i, 3] <- local_maxima[2]
    df_local_maxima_ctrl[i, 4] <-df_normalized_combined[i, local_maxima[2]]
    
  } else if (length(local_maxima) == 3) {
    # proteins with 3 local maxima
    
    df_local_maxima_ctrl[i, 1] <- local_maxima[1]
    df_local_maxima_ctrl[i, 2] <- df_normalized_combined[i, local_maxima[1]]
    df_local_maxima_ctrl[i, 3] <- local_maxima[2]
    df_local_maxima_ctrl[i, 4] <-df_normalized_combined[i, local_maxima[2]]
    df_local_maxima_ctrl[i, 5] <- local_maxima[3]
    df_local_maxima_ctrl[i, 6] <-df_normalized_combined[i, local_maxima[3]]
    
  } else if (length(local_maxima) == 4) {
    # proteins that have 4 local maxima
    
    df_local_maxima_ctrl[i, 1] <- local_maxima[1]
    df_local_maxima_ctrl[i, 2] <- df_normalized_combined[i, local_maxima[1]]
    df_local_maxima_ctrl[i, 3] <- local_maxima[2]
    df_local_maxima_ctrl[i, 4] <-df_normalized_combined[i, local_maxima[2]]
    df_local_maxima_ctrl[i, 5] <- local_maxima[3]
    df_local_maxima_ctrl[i, 6] <-df_normalized_combined[i, local_maxima[3]]
    df_local_maxima_ctrl[i, 7] <- local_maxima[4]
    df_local_maxima_ctrl[i, 8] <-df_normalized_combined[i, local_maxima[4]]
    
  } else if (length(local_maxima) == 5) {
    # proteins with 5  local maxima
    
    df_local_maxima_ctrl[i, 1] <- local_maxima[1]
    df_local_maxima_ctrl[i, 2] <- df_normalized_combined[i, local_maxima[1]]
    df_local_maxima_ctrl[i, 3] <- local_maxima[2]
    df_local_maxima_ctrl[i, 4] <-df_normalized_combined[i, local_maxima[2]]
    df_local_maxima_ctrl[i, 5] <- local_maxima[3]
    df_local_maxima_ctrl[i, 6] <-df_normalized_combined[i, local_maxima[3]]
    df_local_maxima_ctrl[i, 7] <- local_maxima[4]
    df_local_maxima_ctrl[i, 8] <-df_normalized_combined[i, local_maxima[4]]
    df_local_maxima_ctrl[i, 9] <- local_maxima[5]
    df_local_maxima_ctrl[i, 10] <-df_normalized_combined[i, local_maxima[5]]
    
  }
}

# local maxima rnase

# create df that will contain the local maxima, because we do not know yet how many local maxima each protein will have we need to make space for more than one per protein and cut the rest of later
df_local_maxima_rnase <-
  data.frame(
    row.names = rownames(df_normalized_combined),
    rnase_local_maximum1_fraction = matrix(0, nrow(df_normalized_combined), 1),
    rnase_local_maximum1_value = matrix(0, nrow(df_normalized_combined), 1),
    rnase_local_maximum2_fraction = matrix(0, nrow(df_normalized_combined), 1),
    rnase_local_maximum2_value = matrix(0, nrow(df_normalized_combined), 1),
    rnase_local_maximum3_fraction = matrix(0, nrow(df_normalized_combined), 1),
    rnase_local_maximum3_value = matrix(0, nrow(df_normalized_combined), 1),
    rnase_local_maximum4_fraction = matrix(0, nrow(df_normalized_combined), 1),
    rnase_local_maximum4_value = matrix(0, nrow(df_normalized_combined), 1),
    rnase_local_maximum5_fraction = matrix(0, nrow(df_normalized_combined), 1),
    rnase_local_maximum5_value = matrix(0, nrow(df_normalized_combined), 1)
  )

for (i in 1:nrow(df_normalized_combined)) {
  # first loop iterating through the proteins
  over_4 <- c()
  # store indices of values over 4 for every protein in here (only temporary)
  
  for (j in 26:50) {
    # second loop iterating through all Ctrl fractions
    if (df_normalized_combined[i, j] > 4) {
      # only values > 4 are saved
      over_4 <- append(over_4, j)
    }
  }
  #print(over_4)
  local_maxima <- c()
  # store indices of values over 4 and with no higher neighbours (only temporary)
  
  for (v in over_4) {
    # iterate through all values that are > 4
    #print(df_normalized_combined[i, v])
    
    # values on boarders have to be treated differently
    
    # 26th value (no value to the left)
    if (v == 26) {
      if (df_normalized_combined[i, v] > df_normalized_combined[i, v + 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v + 2] &
          df_normalized_combined[i, v] != max(df_normalized_combined[i, 26:50])) {
        local_maxima <- append(local_maxima, v)
      }
      
      # 27th value (only one value to the left)
    } else if (v == 27) {
      if (df_normalized_combined[i, v] > df_normalized_combined[i, v - 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v + 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v + 2] &
          df_normalized_combined[i, v] != max(df_normalized_combined[i, 26:50])) {
        local_maxima <- append(local_maxima, v)
      }
      
      # 49th value (only one value to the right)
    } else if (v == 49) {
      if (df_normalized_combined[i, v] > df_normalized_combined[i, v - 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v - 2] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v + 1] &
          df_normalized_combined[i, v] != max(df_normalized_combined[i, 26:50])) {
        local_maxima <- append(local_maxima, v)
      }
      
      # 50th value (no value to the right)
    } else if (v == 50) {
      if (df_normalized_combined[i, v] > df_normalized_combined[i, v - 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v - 2] &
          df_normalized_combined[i, v] != max(df_normalized_combined[i, 26:50])) {
        local_maxima <- append(local_maxima, v)
      }
      
      # all other values to not have to be treated differently
    } else {
      if (df_normalized_combined[i, v] > df_normalized_combined[i, v - 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v - 2] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v + 1] &
          df_normalized_combined[i, v] > df_normalized_combined[i, v + 2] &
          df_normalized_combined[i, v] != max(df_normalized_combined[i, 26:50])) {
        local_maxima <- append(local_maxima, v)
      }
    }
  }
  #print(paste("protein", as.character(i)))
  #print(local_maxima)
  
  # now we need to save the protein name as well as the local maximum's position and value
  
  if (length(local_maxima) == 1) {
    # proteins with one local maximum
    
    df_local_maxima_rnase[i, 1] <- local_maxima[1] - 25
    df_local_maxima_rnase[i, 2] <- df_normalized_combined[i, local_maxima[1]]
    
  } else if (length(local_maxima) == 2) {
    # proteins with 2 local maxima
    
    df_local_maxima_rnase[i, 1] <- local_maxima[1] - 25
    df_local_maxima_rnase[i, 2] <- df_normalized_combined[i, local_maxima[1]]
    df_local_maxima_rnase[i, 3] <- local_maxima[2] - 25
    df_local_maxima_rnase[i, 4] <-df_normalized_combined[i, local_maxima[2]]
    
  } else if (length(local_maxima) == 3) {
    # proteins with 3 local maxima
    
    df_local_maxima_rnase[i, 1] <- local_maxima[1] - 25
    df_local_maxima_rnase[i, 2] <- df_normalized_combined[i, local_maxima[1]]
    df_local_maxima_rnase[i, 3] <- local_maxima[2] - 25
    df_local_maxima_rnase[i, 4] <-df_normalized_combined[i, local_maxima[2]]
    df_local_maxima_rnase[i, 5] <- local_maxima[3] - 25
    df_local_maxima_rnase[i, 6] <-df_normalized_combined[i, local_maxima[3]]
    
  } else if (length(local_maxima) == 4) {
    # proteins that have 4 local maxima
    
    df_local_maxima_rnase[i, 1] <- local_maxima[1] - 25
    df_local_maxima_rnase[i, 2] <- df_normalized_combined[i, local_maxima[1]]
    df_local_maxima_rnase[i, 3] <- local_maxima[2] - 25
    df_local_maxima_rnase[i, 4] <-df_normalized_combined[i, local_maxima[2]]
    df_local_maxima_rnase[i, 5] <- local_maxima[3] - 25
    df_local_maxima_rnase[i, 6] <-df_normalized_combined[i, local_maxima[3]]
    df_local_maxima_rnase[i, 7] <- local_maxima[4] - 25
    df_local_maxima_rnase[i, 8] <-df_normalized_combined[i, local_maxima[4]]
    
  } else if (length(local_maxima) == 5) {
    # proteins with 5  local maxima
    
    df_local_maxima_rnase[i, 1] <- local_maxima[1] - 25
    df_local_maxima_rnase[i, 2] <- df_normalized_combined[i, local_maxima[1]]
    df_local_maxima_rnase[i, 3] <- local_maxima[2] - 25
    df_local_maxima_rnase[i, 4] <-df_normalized_combined[i, local_maxima[2]]
    df_local_maxima_rnase[i, 5] <- local_maxima[3] - 25
    df_local_maxima_rnase[i, 6] <-df_normalized_combined[i, local_maxima[3]]
    df_local_maxima_rnase[i, 7] <- local_maxima[4] - 25
    df_local_maxima_rnase[i, 8] <-df_normalized_combined[i, local_maxima[4]]
    df_local_maxima_rnase[i, 9] <- local_maxima[5] - 25
    df_local_maxima_rnase[i, 10] <-df_normalized_combined[i, local_maxima[5]]
    
  }
}

# combine local maxima

df_local_maxima <- as.data.frame(cbind(df_local_maxima_ctrl, df_local_maxima_rnase))
rownames(df_local_maxima) <- rownames(df_normalized_combined)





### Analysis part 3


## shift of global maxima > n fractions

# create vector that will contain the shift in fractions for every protein
shift_global <- c()

# iterate through every protein and calculate how large the shift is
for (i in 1:nrow(df_global_maxima)) {
  
  x <- df_global_maxima[i, 3] - df_global_maxima[i, 1]
  shift_global <- append(shift_global, x)
  
}

# now we can see the shift for every protein
names(shift_global) <- rownames(df_global_maxima)

n = 1 # shift must be greater than n
shift_proteins <- c() # in here, the proteins which are shifting are saved
shift_proteins_index <- c()
counter <- 0

for (i in 1:nrow(df_global_maxima)) {
  
  if (shift_global[i] < -n) {
    shift_proteins <- append(shift_proteins, rownames(df_global_maxima)[i])
    shift_proteins_index <- append(shift_proteins_index, i)
  } else if (shift_global[i] > n) {
    shift_proteins <- append(shift_proteins, rownames(df_global_maxima)[i])
    shift_proteins_index <- append(shift_proteins_index, i)
  } else if (shift_global[i] < -n | shift_global[i] > n) {
    counter <- counter + 1
  } 
}

counter # must be zero
length(shift_proteins) # shift proteins gives us the names of the proteins which shift more than 1 fraction
length(shift_proteins_index) # shift proteins index gives us the indices of the shifting proteins



## significant difference of protein amount at the global maximum (replicates)

pval_rep <- c()

for (i in 1:nrow(df_global_maxima_rep)) {
  
  ctrl_vec <-
    c(df_normalized[i, df_global_maxima_rep[i, 1]],
      df_normalized[i, df_global_maxima_rep[i, 3] + 25],
      df_normalized[i, df_global_maxima_rep[i, 5] + 50])
  rnase_vec <-
    c(df_normalized[i, df_global_maxima_rep[i, 1] + 75],
      df_normalized[i, df_global_maxima_rep[i, 3] + 100],
      df_normalized[i, df_global_maxima_rep[i, 5] + 125])
  
  if (abs(sum(ctrl_vec)) > 299.99 & abs(sum(rnase_vec) > 299.99)) {
    pval_rep <- append(pval_rep, 1)
  } else {
    x <- t.test(ctrl_vec, rnase_vec)
    pval_rep <- append(pval_rep, x$p.value)
  }
}

pval_adjust_rep <- p.adjust(pval_rep, method = "fdr")
pval_low_rep <- which(pval_adjust_rep < 0.05)
sig_dif_proteins_rep <- rownames(df_normalized)[pval_low_rep]
sig_dif_proteins_index_rep <- pval_low_rep

length(sig_dif_proteins_rep) # names of proteins with significant t test (replicate values)
length(sig_dif_proteins_index_rep) # index of proteins with significant t test (replicate values)



## RNA-dependent proteins

rdeep_rep <- c(intersect(sig_dif_proteins_rep, shift_proteins))
rdeep_index_rep <- c(intersect(sig_dif_proteins_index_rep, shift_proteins_index))

length(rdeep_rep)
length(rdeep_index_rep)
df_rdeep <- data.frame(row.names = rdeep_rep, index = rdeep_index_rep)




### Analysis part 4



## Dimension Reduction

df <- df_normalized_combined # we can use this variable to easily switch between combined and not combined

# df_umap1: pca and umap
pca <- prcomp(df)
df_umap1 <- as.data.frame(umap(pca$x))

# df_umap2 : only umap
df_umap2 <- as.data.frame(umap(df))

# df_pca: only pca
df_pca <- as.data.frame(prcomp(df)$x)

ggplot(df_umap1, aes(x = V1, y = V2)) +
  geom_point()

ggplot(df_umap2, aes(x = V1, y = V2)) +
  geom_point()

ggplot(df_pca, aes(x = PC1, y = PC2)) +
  geom_point()

# add shift information in order to colour proteins
n = 1
df_umap <- cbind(df_umap1, matrix(0, nrow(df_normalized_combined), 1))
colnames(df_umap) <- c("V1", "V2", "shift")

for (i in 1:nrow(df_umap)) {
  if (shift_global[i] < -n) {
    df_umap[i, 3] <- "left shift"
  } else if (shift_global[i] > n) {
    df_umap[i, 3] <- "right shift"
  } else {
    df_umap[i, 3] <- "no shift"
  }
}

# add information in order to colour proteins
for (i in 1:nrow(df_umap)) {
  if (rownames(df_umap)[i] %in% rownames(df_rdeep)) {
    df_umap$rdeep[i] <- "rdeep"
  } else {
    df_umap$rdeep[i] <- "non rdeep"
  }
}

alpha <- ifelse(df_umap$shift == "no shift", 0.1, 1)

ggplot(df_umap, aes(x = V1, y = V2, color = shift)) +
  geom_point(alpha = alpha) +
  ggtitle("UMAP of normalized and combined data")

ggplot(df_umap, aes(x = V1, y = V2, color = rdeep)) +
  geom_point() +
  ggtitle("UMAP of normalized and combined data")

x <- grep("left", df_umap$shift)
y <- grep("right", df_umap$shift)
z <- c(x, y)

df_umap_shift <- df_umap[z, ]
df_umap_no_shift <- df_umap[-z, ]

ggplot(df_umap_shift, aes(x = V1, y = V2, color = shift)) +
  geom_point()
ggplot(df_umap_no_shift, aes(x = V1, y = V2, color = shift)) +
  geom_point()

n = 1
df_pca_shift <- cbind(df_pca, matrix(0, nrow(df_normalized_combined), 1))
colnames(df_pca_shift)[51] <- "shift"

for (i in 1:nrow(df_pca_shift)) {
  if (shift_global[i] < -n) {
    df_pca_shift[i, 51] <- "left shift"
  } else if (shift_global[i] > n) {
    df_pca_shift[i, 51] <- "right shift"
  } else {
    df_pca_shift[i, 51] <- "no shift"
  }
  
}

ggplot(df_pca_shift, aes(x = PC1, y = PC2, color = shift)) +
  geom_point()



## group proteins according to their shifting behaviour

n = 1
df_global_maxima_shift <- cbind(df_global_maxima, matrix(0, nrow(df_normalized_combined), 1))
colnames(df_global_maxima_shift)[5] <- "shift"

for (i in 1:nrow(df_global_maxima_shift)) {
  if (shift_global[i] < -n) {
    df_global_maxima_shift[i, 5] <- "left shift"
  } else if (shift_global[i] > n) {
    df_global_maxima_shift[i, 5] <- "right shift"
  } else {
    df_global_maxima_shift[i, 5] <- "no shift"
  }
  
}

ggplot(
  df_global_maxima_shift,
  aes(y = ctrl_global_maximum_fraction, x = rnase_global_maximum_fraction, color = shift)) +
  geom_count() +
  ggtitle("shifting behaviour of proteins") +
  xlab("global maximum RNase condition") + ylab("global maximum Ctrl conditions")

ggplot(
  df_global_maxima_shift,
  aes(y = ctrl_global_maximum_fraction, x = rnase_global_maximum_fraction, color = shift)) +
  geom_jitter() +
  ggtitle("shifting behaviour of proteins") +
  xlab("global maximum RNase condition") + ylab("global maximum Ctrl conditions")



## kmeans with protein amount

df_cluster <- as.data.frame(cbind(shift_global, pval_adjust_rep))
rownames(df_cluster) <- rownames(df_normalized)

for (i in 1:nrow(df_cluster)) {
  if (rownames(df_cluster)[i] %in% rownames(df_rdeep)) {
    df_cluster$rdeep[i] <- "rdeep"
  } else {
    df_cluster$rdeep[i] <- "non rdeep"
  }
}

df_cluster$ctrl_maximum <- df_global_maxima$ctrl_global_maximum_fraction
df_cluster$rnase_maximum <- df_global_maxima$rnase_global_maximum_fraction

wss = sapply(2:10, function(x) { 
  kmeans(df_normalized_combined, centers = x, nstart = 100)$tot.withinss
})
plot(2:10, wss, xlab = "numbers of centers", ylab = "total within-cluster sum of squares", type = "b")

tdf_normalized_combined <- t(df_normalized_combined)
km <- kmeans(df_normalized_combined, centers = 3, nstart = 100)

fviz_cluster(km, data = df_normalized_combined,
             #palette = c("#2E9FDF", "#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

df_cluster$cluster <- as.character(km$cluster)

ggplot(
  df_cluster,
  aes(y = ctrl_maximum, x = rnase_maximum, color = cluster)) +
  geom_count() +
  ggtitle("shifting behaviour of proteins") +
  xlab("global maximum RNase condition") + ylab("global maximum Ctrl conditions")

df_umap$cluster <- as.character(km$cluster)

ggplot(df_umap, aes(x = V1, y = V2, color = cluster)) +
  geom_point() +
  ggtitle("UMAP of normalized and combined data")



## kmeans with shift

km2 <- kmeans(df_cluster[,1], centers = 3)

df_umap$cluster2 <- as.character(km2$cluster)

ggplot(df_umap, aes(x = V1, y = V2, color = cluster2)) +
  geom_point() +
  ggtitle("UMAP of normalized and combined data")






### Analysis part 5



## predict shift with correlation

tdf <- as.data.frame(t(df_normalized_combined))
cor(tdf$EED_HUMAN[1:25], tdf$EED_HUMAN[26:50])  # RDeeP
cor(tdf$`1433B_HUMAN`[1:25], tdf$`1433B_HUMAN`[26:50]) # non-RDeeP

# create df that will contain important information for linear modelling
df_lm <- data.frame(row.names = colnames(tdf))

# add vector containing correlation between rnase and control condition to df
vec_cor <- apply(tdf, 2, function(x) {
  cor(x[1:25], x[26:50])
})
df_lm$correlation <- vec_cor
df_lm$shift <- shift_global

x <- round(nrow(df_lm) * 0.8)
y <- x + 1
df_lm_train <- df_lm[1:x, ]
df_lm_test <- df_lm[y:nrow(df_lm), ]

linear_model <- lm(shift ~ correlation, data = df_lm_train)
summary(linear_model)

# normal distribution of residuals? 
hist(linear_model$residuals, breaks = 20)

qqnorm(linear_model$residuals)
qqline(linear_model$residuals)

## correlation residuals x-values? 
cor(df_lm_train$correlation, linear_model$residuals)
plot(df_lm_train$correlation, linear_model$residuals, pch = 20)

pm <- predict.lm(linear_model, newdata = df_lm_test, se.fit = TRUE, interval = "confidence")

df_lm_predict <- as.data.frame(cbind(df_lm_test$shift, pm$fit))

plot(
  df_lm_predict$V1,
  df_lm_predict$fit,
  pch = 20,
  col = 'blue',
  xlab = 'Real values',
  ylab = 'Predicted values'
); abline(0, 1, col = "blue", lwd = 5)




## same model but absolut shift values

# df for model training
df_lm_train_ab <- df_lm_train
df_lm_train_ab$shift <- abs(df_lm_train_ab$shift)

# df for model testing
df_lm_test_ab <- df_lm_test
df_lm_test_ab$shift <- abs(df_lm_test_ab$shift)

linear_model_ab <- lm(shift ~ correlation, data = df_lm_train_ab)
summary(linear_model_ab)

# normal distribution of residuals? 
hist(linear_model_ab$residuals, breaks = 20)

qqnorm(linear_model_ab$residuals)
qqline(linear_model_ab$residuals)

## correlation residuals x-values? 
cor(df_lm_train_ab$correlation, linear_model_ab$residuals)
plot(df_lm_train_ab$correlation, linear_model_ab$residuals, pch = 20)

pm_ab <- predict.lm(linear_model_ab, newdata = df_lm_test_ab, se.fit = TRUE, interval = "confidence")

df_lm_predict_ab <- as.data.frame(cbind(df_lm_test_ab$shift, pm_ab$fit))

plot(
  df_lm_predict_ab$V1,
  df_lm_predict_ab$fit,
  pch = 20,
  col = 'blue',
  xlab = 'Real values',
  ylab = 'Predicted values'
); abline(0, 1, col = "blue", lwd = 5)











