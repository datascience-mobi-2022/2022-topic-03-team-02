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

sum(apply(MS_Table, 2, is.numeric)) == ncol(MS_Table)
glimpse(MS_Table)

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

##normalize replicates

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








