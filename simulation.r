suppressPackageStartupMessages(library(tidyr, quietly = TRUE))
suppressPackageStartupMessages(library(dplyr, quietly = TRUE))
suppressPackageStartupMessages(library(ggplot2, quietly = TRUE))
suppressPackageStartupMessages(library(magrittr, quietly = TRUE))

theme_set(
    theme_bw() +
        theme(axis.text = element_text(size = 20), legend.position = "none", legend.text = element_text(size = 20), legend.title = element_blank(), axis.title = element_text(size = 30), legend.key = element_blank(), legend.background = element_blank())
)
options(scipen = 999)

# Mother matrix initialization
cat("Enter number of rows in matrix:\n")
row_num <- as.numeric(readline())
if ((row_num %% 2) != 0) {
    cat("\nEnter even number of rows\n")
    quit
}
cat("Enter number of columns in matrix:\n")
col_num <- as.numeric(readline())
mat1 <- matrix(0, nrow = row_num, ncol = col_num)

matrix_initiate <- function(matrix, eu_hist_number, het_hist_number) {
    # takes in matrix(0,0) & ouputs matrix with 1s randomly distributed based on exponential distribution
    num1 <- row_num / 2
    avg_hist_euchrom <- c(eu_hist_number)
    avg_hist_het <- c(het_hist_number)
    num_hist_euchrom <- ceiling(rexp(num1, 1 / avg_hist_euchrom)) # using exponential distribution (with given mean) to distribute histones in first-half of the mother matrix
    num_hist_euchrom <- replace(num_hist_euchrom, num_hist_euchrom > col_num, col_num)
    num_hist_het <- ceiling(rexp(num1, 1 / avg_hist_het)) # using exponential distribution (with given mean) to distribute histones in second-half of the mother matrix
    num_hist_het <- replace(num_hist_het, num_hist_het > col_num, col_num)
    num_hist_vec <- c(num_hist_euchrom, num_hist_het)
    list_sample_eu <- list()
    list_sample_het <- list()

    # sampling where to put histones randomly in the mother matrix
    for (i in 1:nrow(matrix)) {
        if (i > 0 && i <= num1) {
            sampler <- sample(1:col_num, num_hist_vec[i])
            list_sample_eu[[i]] <- c(sampler)
        }

        if (i > num1 && i <= row_num) {
            sampler <- sample(1:col_num, num_hist_vec[i])
            list_sample_het[[i]] <- c(sampler)
        }
    }

    # Putting histones in random positions in mother matrix
    l <- 1
    for (i in 1:nrow(matrix)) {
        if (i > 0 && i <= num1) {
            index <- list_sample_eu[i]
            for (j in 1:length(index[[l]])) {
                k <- index[[l]][j]
                matrix[i, k] <- 1
            }
        }

        if (i > num1 && i <= row_num) {
            index <- list_sample_het[i]
            for (j in 1:length(index[[l]])) {
                k <- index[[l]][j]
                matrix[i, k] <- 1
            }
        }
    }
    return(matrix)
}


cat("\n")
# cat("Enter the Diffusion Constant value:\n")
diff_constant_e <- 1.2    #c(1.2) --> for main text   #c(3.533) --> Escobar fits
diff_constant_l <- 0.07   #c(0.07)                   #c(0.07)
diff_constant_el <-0.0001 #c(0.0001)                #c(0.01)
# print(mat1)

# Daughter matrices initialization
dat1 <- matrix(2, nrow = row_num, ncol = col_num)
dat2 <- dat1


# to calculate all the possible distance dependent rates between 2 matrices
def_distance <- function(matrix, row_number, col_number) {
    xp <- as.data.frame(which(matrix < 3, arr.ind = TRUE))
    xp <- arrange(xp, row)
    uni_row <- unique(xp$row)
    uni_col <- unique(xp$col)
    dist <- c()
    val <- c()
    constant <- 0.5
    for (i in 1:length(uni_col)) {
        for (j in 1:length(uni_row)) {
            val[j] <- sqrt((i - uni_row[j])^2)
            dist <- c(dist, val[j])
        }
    }
    all_rates <- 1 / (constant + dist)
    rate_mat <- matrix(all_rates, nrow = row_number, ncol = col_number, byrow = FALSE)
    dist_mat <- matrix(dist, nrow = row_number, ncol = col_number, byrow = FALSE)
    return(dist_mat)
}

# defining function to select between daughter matrices with probability of 0.5
select_matrix <- function(x, y) {
    rando <- runif(1, min = 0, max = 1)
    if (rando > 0.50 & rando <= 1) {
        out <- x
    } else {
        out <- y
    }
    return(out)
}

# Calculating probabilities using the diffusion constants
def_probability <- function(matrix, row_number, column_number, diff_constant_e, diff_constant_l, diff_constant_el) { # nolint
    dimen <- as.numeric(((row_number)^2) * 2)
    dimen2 <- as.numeric((row_number)^2)
    num1 <- row_number / 2
    num2 <- row_number * 2
    num3 <- dimen2 / 2
    prob_index1 <- seq(1, dimen2, by = row_number)
    prob_index2 <- seq(row_number, dimen2, by = row_number)
    prob_index3 <- seq(1, dimen2, by = num1)
    prob_index4 <- seq(num1, dimen2, by = num1)
    var <- def_distance(matrix, row_number, column_number)
    rates <- c(var[1:row_number, 1:row_number])
    change_var <- seq(1, num2, 2)
    rate <- c()
    for (i in 1:row_number) {
        # for loop to calculate probabilities based on sum of rates
        if (i <= num1) {
            index1 <- change_var[i]
            p <- prob_index3[index1]
            q <- prob_index4[index1]
            rate_val <- c(rates[p:q])
            rate <- c(rate, rate_val + 1)     # Early to Early
            rate <- c(rate, rate_val + 1)     # Early to Late
        }

        if (i > num1) {
            index1 <- change_var[i] + 1
            r <- prob_index3[index1]
            s <- prob_index4[index1]
            rate_val <- c(rates[r:s])
            rate <- c(rate, rate_val + 1)    # Late to Early 
            rate <- c(rate, rate_val)        # Late to Late
        }
    }
    # rate <- replace(rate, rate[1:num3] == 1, 0)
    rates <- rate
    dist_mat <- rate
    time <- c(1)
    pie <- c(3.1415)
    func_part.1 <- 1 / sqrt(4 * pie * diff_constant_e * time)
    func_part.2 <- 4 * diff_constant_e * time
    func_part.3 <- 1 / sqrt(4 * pie * diff_constant_l * time)
    func_part.4 <- 4 * diff_constant_l * time
    func_part.5 <- 1 / sqrt(4 * pie * diff_constant_el * time)
    func_part.6 <- 4 * diff_constant_el * time
    prob <- c()
    for (i in 1:num2) {
        # for loop to calculate probabilities based on sum of rates
        p <- prob_index3[i]
        q <- prob_index4[i]
        value <- sum(rates[p:q])
        rate_val <- c(rates[p:q])
        if (i <= row_number) {
            if ((i %% 2) == 0) {
                str <- as.vector(func_part.5 * exp(-(rate_val)^2 / func_part.6) / 2)
                prob <- c(prob, str)
            } else {
                str <- as.vector(func_part.1 * exp(-(rate_val)^2 / func_part.2) / 2)
                prob <- c(prob, str)
            }
        }

        if (i > row_number) {
            if ((i %% 2) == 0) {
                str <- as.vector(func_part.3 * exp(-(rate_val)^2 / func_part.4) / 2)
                prob <- c(prob, str)
            } else {
                str <- as.vector(func_part.5 * exp(-(rate_val)^2 / func_part.6) / 2)
                prob <- c(prob, str)
            }
        }
    }

    rate <- c()
    probability <- c()
    for (i in 1:row_number) {
        # for loop to make copy of rate and probability for 2 daughter matrices
        f <- prob_index1[i]
        g <- prob_index2[i]
        val <- c(rates[f:g], rates[f:g])
        vari <- c(prob[f:g], prob[f:g])
        rate <- c(rate, val)
        probability <- c(probability, vari)
    }

    # simplifying the distance matrix
    dist_mat <- matrix(dist_mat, nrow = row_number, ncol = row_number)
    dist_vec <- c()
    for (i in 1:nrow(dist_mat)) {
        if (i == nrow(dist_mat)) {
            next
        }
        k <- i + 1
        for (j in k:nrow(dist_mat)) {
            dist_val <- dist_mat[i, j]
            dist_vec <- c(dist_vec, dist_val)
        }
    }

    probability <- split(probability, ceiling(seq_along(probability) / num2))
    return(list(dist_vec, probability))
}

mark_copy_new <- function(matrix, column_number) { # nolint
    var <- 2:(column_number - 1)
    for (i in 1:ncol(matrix)) {
        for (j in 1:nrow(matrix)) {
            if (matrix[j, i] == 0 || matrix[j, i] == 1) {
                next
            }
            if (i == 1 && matrix[j, i] == 2 && matrix[j, i] != matrix[j, i + 1]) {
                matrix[j, i] <- matrix[j, i + 1]
                next
            }
            if (i == 1 && matrix[j, i] == 2 && matrix[j, i] == matrix[j, i + 1]) {
                matrix[j, i] <- sample(0:1, 1)
                next
            }

            if (i == ncol(matrix) && matrix[j, i] == 2 && matrix[j, i] != matrix[j, i - 1]) { # nolint
                matrix[j, i] <- matrix[j, i - 1]
                next
            }
            if (i == ncol(matrix) && matrix[j, i] == 2 && matrix[j, i] == matrix[j, i - 1]) { # nolint
                matrix[j, i] <- sample(0:1, 1)
                next
            }
            if (i %in% var && matrix[j, i] == 2 && matrix[j, i] == matrix[j, i - 1] && matrix[j, i] == matrix[j, i + 1]) { # nolint
                matrix[j, i] <- sample(0:1, 1)
                next
            }
            if (i %in% var && matrix[j, i] == 2 && matrix[j, i - 1] == 0 && matrix[j, i + 1] == 0) { # nolint
                matrix[j, i] <- 0
                next
            }
            if (i %in% var && matrix[j, i] == 2 && matrix[j, i - 1] == 1 && matrix[j, i + 1] == 1) { # nolint
                matrix[j, i] <- 1
                next
            }
            if (i %in% var && matrix[j, i] == 2 && matrix[j, i - 1] == matrix[j, i] && matrix[j, i + 1] == 1) { # nolint
                matrix[j, i] <- 1
                next
            }
            if (i %in% var && matrix[j, i] == 2 && matrix[j, i - 1] == matrix[j, i] && matrix[j, i + 1] == 0) { # nolint
                matrix[j, i] <- 0
                next
            }
            if (i %in% var && matrix[j, i] == 2 && matrix[j, i - 1] == 1 && matrix[j, i + 1] == matrix[j, i]) { # nolint
                matrix[j, i] <- 1
                next
            }
            if (i %in% var && matrix[j, i] == 2 && matrix[j, i - 1] == 0 && matrix[j, i + 1] == matrix[j, i]) { # nolint
                matrix[j, i] <- 0
                next
            }
            if (i %in% var && matrix[j, i] == 2 && matrix[j, i - 1] == 1 && matrix[j, i + 1] == 0) { # nolint
                matrix[j, i] <- sample(0:1, 1)
                next
            }
            if (i %in% var && matrix[j, i] == 2 && matrix[j, i - 1] == 0 && matrix[j, i + 1] == 1) { # nolint
                matrix[j, i] <- sample(0:1, 1)
                next
            }
        }
    }
    return(matrix)
}


Dij_calculate <- function(matrix) {
    dij_vec <- c()
    for (i in 1:nrow(matrix)) {
        if (i == nrow(matrix)) {
            next
        }
        count1 <- sum(matrix[i, ] == 1)
        k <- i + 1
        for (j in k:nrow(matrix)) {
            count2 <- sum(matrix[j, ] == 1)
            dij <- abs(log10(count1 / count2))
            dij_vec <- c(dij_vec, dij)
        }
    }
    return(dij_vec)
}

rep_time <- function(matrix) { # calculating Dij # nolint
    time_vec <- c()
    nrow <- nrow(matrix)
    for (i in 1:nrow(matrix)) {
        if (i == nrow(matrix)) {
            next
        }
        k <- i + 1
        for (j in k:nrow(matrix)) {
            if ((i <= (nrow / 2)) && (j <= (nrow / 2))) {
                temp <- c("Early-Early")
                time_vec <- c(time_vec, temp)
            }

            if ((i > (nrow / 2)) && (j > (nrow / 2))) {
                temp <- c("Late-Late")
                time_vec <- c(time_vec, temp)
            }

            if ((i <= (nrow / 2)) && (j > (nrow / 2))) {
                temp <- c("Early-Late")
                time_vec <- c(time_vec, temp)
            }
        }
    }
    return(time_vec)
}

# creating multiple variables to loop between indices during simulation run
dimen <- ((row_num)^2) * 2
dimen2 <- (row_num)^2
prob_index3 <- seq(1, dimen2, by = row_num)
prob_index4 <- seq(row_num, dimen2, by = row_num)
row_index <- seq(1, (row_num * 2), by = 1)
options(scipen = 999) # turns off scientific notation
probability <- def_probability(mat1, row_num, col_num, diff_constant_e, diff_constant_l, diff_constant_el) # calling function which returns probability values based on distances
dist_vec <- probability[[1]]
probability <- probability[[2]] # calling function which returns probability values based on distances


# to store the list of matrices generated during simulation
check_chr <- c(1:row_num)
check_mat <- c(1:row_num)
v <- row_num + 1
check_mat3 <- c(v:(row_num * 2))
check_mat2 <- c(1:row_num, 1:row_num)
prob_var <- 1 / (row_num)
row_prob <- rep(prob_var, row_num)
cell_cycles <- c(20)
iteration_num <- c(20)     # 100 used for paper
spatial_proximity <- (-dist_vec)
hist_diff <- c()  

# to plot D(i,j) vs spatial proximity from first initialized mother matrix
# This is to compare the first state of histones in chromatin vs the change as multiple cell cycles occur over a period of time
cat("\n")
cat("Plotting the initial state...")
cat("\n")
dij <- Dij_calculate(mat1)


Replication_time <- rep_time(mat1)

df.plot.1 <- data.frame(spatial_proximity, dij, Replication_time)
for (x in iteration_num) {
    mat1 <- matrix(0, nrow = row_num, ncol = col_num)
    mat1 <- matrix_initiate(mat1, 130, 130)
    dij <- Dij_calculate(mat1)
    temp <- data.frame(spatial_proximity, dij, Replication_time)
    df.plot.1 <- rbind(df.plot.1, temp)
}
myplot0 <- ggplot(df.plot.1, aes(x = spatial_proximity, y = dij, fill = Replication_Time)) +
    xlab("-Distance") +
    ylab("Q(i,j)") +
    coord_cartesian(ylim = c(0, 1.5)) +
    theme_bw() +
    theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    geom_boxplot(aes(fill = factor(Replication_time, levels = c("Early-Early", "Late-Late", "Early-Late")), cut_width(spatial_proximity, 3, boundary = 0)), outlier.alpha = 0.1)

pdf(paste0("initial_state.pdf"))
suppressWarnings(print(myplot0))
invisible(dev.off())
cat("Done!")
cat("\n")

mat1 <- matrix(0, nrow = row_num, ncol = col_num)
mat1 <- matrix_initiate(mat1, 130, 130)
print("ANCESTOR MATRIX")
for (i in 1:nrow(mat1)) {
  print(c(length(which(mat1[i,]==0)),length(which(mat1[i,]==1))))
}

cat("\n")
cat("Running Simulation...")
cat("\n")

# simulation start
for (y in 1:iteration_num) {
    for (z in 1:cell_cycles) {
        for (a in 1:ncol(mat1)) {
            var1 <- row_prob
            var2 <- check_chr
            var3 <- row_index
            var4 <- probability
            for (b in 1:nrow(mat1)) {
                chrNo <- sample(var2, 1, prob = var1) # randomly choosing row (chromosome)
                var1[chrNo] <- 0 # reducing probability for chromosome to 0 for next iteration
                variable <- var4[[chrNo]] # storing probability values in variable for the particular chr No./ row No.
                val <- mat1[chrNo, a] # catching the value at the index
                index <- sample(var3, 1, prob = variable) # choosing where to transfer based on probabilities
                var4[names(var4)] <- Map(`[<-`, var4[names(var4)], index, 0) # removing the probability from list of vectors for next iterarion

                if (index %in% check_mat == TRUE) {
                    # transferring the index value/ Histone
                    dat1[index, a] <- val
                }
                if (index %in% check_mat3 == TRUE) {
                    # transferring the index value/ Histone
                    index2 <- check_mat2[index]
                    dat2[index2, a] <- val
                }
            }
        }
        # mark-copying (before selection of daughter matrices)
        dat1 <- mark_copy_new(dat1, col_num)
        dat2 <- mark_copy_new(dat2, col_num)
        mat1 <- select_matrix(dat1, dat2) # randomly choosing between daughter matrices and assigning as mother matrix # nolint
        dat1 <- matrix(2, nrow = row_num, ncol = col_num) # reassigning daughter matrices to all zero values
        dat2 <- dat1
        if (z == cell_cycles) {
            dij <- Dij_calculate(mat1)
        }
    } # finishes loop over cell_cycles
  
    if (y == 1) {
        # Replication_time <- c(rep("Early-Early", 0.25 * dimen2), rep("Early-Late", 0.50 * dimen2), rep("Late-Late", 0.25 * dimen2))
        df.plot <- data.frame(spatial_proximity, dij, Replication_time)
    } else {
        # Replication_time <- c(rep("Early-Early", 0.25 * dimen2), rep("Early-Late", 0.50 * dimen2), rep("Late-Late", 0.25 * dimen2))
        temp_df <- data.frame(spatial_proximity, dij, Replication_time)
        df.plot <- rbind(df.plot, temp_df)
    }
  
    # Compute average histone level in early and late segments and record the difference
    sum_early=0
    sum_late=0
    for (i in 1:row_num) {
      if (i<=(row_num/2)) {sum_early=sum_early+length(which(mat1[i,]==1))} 
      if (i>(row_num/2)) {sum_late=sum_late+length(which(mat1[i,]==1))}
    }
    hist_diff=c(hist_diff, (sum_early/(row_num/2) - sum_late/(row_num/2)))
    
    # Reset matrices before next iteration
    if (y == iteration_num) next
    mat1 <- matrix(0, nrow = row_num, ncol = col_num)
    mat1 <- matrix_initiate(mat1, 130, 130)
    dat1 <- matrix(2, nrow = row_num, ncol = col_num)
    dat2 <- dat1
}   # Finishes loop over iteration_num

cat("\n")
cat("Simulation Run Complete!")
cat("\n")


#Plotting the results
df.plot <- subset(df.plot, Replication_time!="Early-Late")
myplot <- ggplot(df.plot, aes(x = spatial_proximity, y = dij, fill = Replication_Time)) +
    xlab("-Distance") +
    ylab("Q(i,j)") +
    coord_cartesian(ylim = c(0, 0.6)) +
    theme(axis.line = element_line(colour = "black")) + #theme_bw() + #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank() +
    geom_boxplot(aes(fill = factor(Replication_time, levels = c("Early-Early", "Late-Late", "Early-Late")), cut_width(spatial_proximity, 3, boundary = 0)), outlier.alpha = 0.1)
myplot + theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))
pdf(paste0("after_simulation.pdf"))
suppressWarnings(print(myplot))
invisible(dev.off())

# Printing the final matrix and writing out signal values
print("FINAL MATRIX")
sum_early=0
sum_late=0
signal_matrix=data.frame()
for (i in 1:row_num) {
  if (i<=(row_num/2)) {
    sum_early=sum_early+length(which(mat1[i,]==1))
    signal_matrix[i,"Signal"]=length(which(mat1[i,]==1))
    signal_matrix[i,"Replication_Time"] = "Early-Replicating DNA"
    } 
  if (i>(row_num/2)) {
    sum_late=sum_late+length(which(mat1[i,]==1))
    signal_matrix[i,"Signal"]=length(which(mat1[i,]==1))
    signal_matrix[i,"Replication_Time"] = "Late-Replicating DNA"
    }
  print(c(length(which(mat1[i,]==0)),length(which(mat1[i,]==1))))
}
#cat("Average histone mark level in Early is", sum_early/(row_num/2))
#cat("Average histone mark level in Late is", sum_late/(row_num/2))
#write.csv(signal_matrix,"Desktop/GM06990_1Mb_CombinedModificationPlots/signal.csv", 
#            quote = FALSE, row.names = FALSE)


# Extracting the median dij values and performing linear regression.
med_EE <- layer_data(myplot)[4][[1]][c(1,3,5,7)] # Extracts early-early median values, ONLY for rows = 20.
med_LL <- layer_data(myplot)[4][[1]][c(2,4,6,8)] # Extracts late-late median values, ONLY for rows = 20.
xaxis <- c(-10.5,-7.5,-4.5,-1.5)                 # X-axis values, ONLY for rows = 20
reg_df <- data.frame(EE=med_EE, LL=med_LL, Dist=xaxis)
reg_df <- gather(reg_df, key = "Replication Time", value = "Median Dij", 1:2)

regress_plot <- ggplot(reg_df, aes(Dist,`Median Dij`, color=`Replication Time`)) +
          geom_point(size=4.0) +
          geom_smooth(method='lm') +
          theme_bw()+ 
          xlab("-Distance")+ylab("Median Q(i,j)")+
          theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))

pdf(paste0("regression.pdf"))
suppressWarnings(print(regress_plot))
invisible(dev.off())
reg_EE <- lm(med_EE ~ xaxis)
reg_LL <- lm(med_LL ~ xaxis)

# Plot distribution of signal differences 
dist_signal <- data.frame(hist_diff = hist_diff)
sig_diff_plot <- ggplot(dist_signal, aes(x=hist_diff), ) + 
  geom_density(color="black", fill="lightpink", alpha=0.4) +
  xlab("Signal difference")+ylab("Density")+
  theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20))
  
pdf(paste0("signal_difference.pdf"))
suppressWarnings(print(sig_diff_plot))
invisible(dev.off())
cat(paste0("\nAll Plots Saved in the current working directory."))
cat("\n")