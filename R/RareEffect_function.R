# Read Group file and split variants by functional annotations
read_groupfile <- function(groupfile_name, gene_name) {
    groupfile <- file(groupfile_name, "r")
    line <- 0
    var <- NULL
    anno <- NULL

    while (TRUE) {
        line <- line + 1
        marker_group_line <- readLines(groupfile, n = 1)

        if (length(marker_group_line) == 0) {
            break
        }

        marker_group_line_list <- strsplit(marker_group_line, split = c(" +", "\t"))[[1]]

        if (marker_group_line_list[1] == gene_name) {
            if (marker_group_line_list[2] == "var") {
                var <- marker_group_line_list
            } else {
                anno <- marker_group_line_list
            }
        }
    }

    lof_idx <- which(anno == "lof")
    mis_idx <- which(anno == "missense")
    syn_idx <- which(anno == "synonymous")

    lof_var <- var[lof_idx]
    mis_var <- var[mis_idx]
    syn_var <- var[syn_idx]

    out <- list(lof_var, mis_var, syn_var)
    close(groupfile)
    return(out)
}

read_matrix_by_one_marker <- function(objGeno, var_list, sampleID) {
    n_samples <- length(sampleID)
    mat <- Matrix::Matrix(0, nrow = n_samples, ncol = 0, sparse = TRUE)
    if (length(var_list) == 0) {
        print("No variants in the list")
        return(mat)
    }
    for (i in 1:length(var_list)) {
        t_GVec <- rep(0, n_samples)
        idx <- which(var_list[i] == objGeno$markerInfo$ID)
        SAIGE::Unified_getOneMarker(t_genoType = objGeno$genoType,
            t_gIndex_prev = objGeno$markerInfo$genoIndex_prev[idx],
            t_gIndex = objGeno$markerInfo$genoIndex[idx],
            t_ref = "2",
            t_alt = "1",
            t_marker = objGeno$markerInfo$ID[idx],
            t_pd = objGeno$markerInfo$POS[idx],
            t_chr = toString(objGeno$markerInfo$CHROM[idx]),
            t_altFreq = 0,
            t_altCounts = 0,
            t_missingRate = 0,
            t_imputeInfo = 0,
            t_isOutputIndexForMissing = TRUE,
            t_indexForMissing = 0,
            t_isOnlyOutputNonZero = FALSE,
            t_indexForNonZero = 0,
            t_GVec = t_GVec,
            t_isImputation = FALSE
        )

        t_GVec_sp <- as(t_GVec, "sparseVector")
        t_GVec_sp_mat <- as(t_GVec_sp, "Matrix")
        
        mat <- cbind(mat, t_GVec_sp_mat)
    }

    colnames(mat) <- var_list
    rownames(mat) <- sampleID

    return(mat)
}

collapse_matrix <- function(objGeno, var_list, sampleID, modglmm, macThreshold = 10) {
    n_samples <- length(sampleID)
    mat <- Matrix::Matrix(0, nrow = n_samples, ncol = 0, sparse = TRUE)
    if (length(var_list) == 0) {
        print("No variants in the list")
        return(mat)
    }
    mat <- read_matrix_by_one_marker(objGeno, var_list, sampleID)
    mat <- mat[which(rownames(mat) %in% modglmm$sampleID), ]
    MAF <- colSums(mat) / (2 * nrow(mat))
    idx <- which(((MAF < 0.01) | (MAF > 0.99)) & ((MAF < 1) & (MAF > 0)))
    mat <- mat[, idx]
    mat_rare <- mat[, which(colSums(mat) >= macThreshold)]
    mat_UR <- mat[, which(colSums(mat) < macThreshold)]
    UR_rowsum <- rowSums(mat_UR)
    UR_rowsum[which(UR_rowsum > 1)] <- 1
    mat_UR_collapsed <- cbind(mat_rare, UR_rowsum)
    return(mat_UR_collapsed)
}

get_range <- function(v) {
    # input are vectors of variants by functional annotation
    start_pos <- Inf
    end_pos <- 0
    for (i in 1:3) {
        if (length(v[[i]]) > 0) {
            pos <- as.numeric(str_split_fixed(v[[i]], ":", 4)[,2]) # Modified in 0.3
            sub_start_pos <- min(pos)
            sub_end_pos <- max(pos)
            print(sub_start_pos)
            print(sub_end_pos)
            if (start_pos > sub_start_pos) {
                start_pos <- sub_start_pos
            }
            if (end_pos < sub_end_pos) {
                end_pos <- sub_end_pos
            }
        }
    }
    return (list(start_pos, end_pos))
}


# Read phenotype file
read_pheno <- function(pheno_file, pheno_code, iid_col = "f.eid") {
    pheno <- data.table::fread(pheno_file, quote = "")
    out <- subset(pheno, select = c(iid_col, pheno_code))
    return(out)
}


calc_log_lik <- function(delta, S, UtY, Y_UUtY) {
    k <- length(S)
    n <- nrow(Y_UUtY)

    log_lik1 <- 0
    for (i in 1:k) {
        log_lik1 <- log_lik1 + (UtY[i, ])^2 / (S[i] + delta)
    }

    log_lik2 <- 1 / delta * sum((Y_UUtY)^2)

    out <- -0.5 * (n * log(2 * pi) + sum(log(S + delta)) + (n - k) * log(delta)
                   + n + n * log(1 / n * (log_lik1 + log_lik2)))

    return(as.numeric(out))
}

calc_post_beta <- function(K, G, delta, S, UtY, U) {
    K_sparse <- as(K, "dgCMatrix")
    if (length(S) == 1) {
        S <- as.matrix(S)
    }

    out <- K_sparse %*% t(G) %*% U %*% diag(1 / (S + delta)) %*% (UtY)
    return(out)
}

# Run FaST-LMM to obtain posterior beta
fast_lmm <- function(G, Y) {
    # if sum(G) == 0, let effect size = 0
    print("Estimating beta using FaST-LMM")
    print(dim(G))
    G[is.na(G)] <- 0
    print(sum(G))
    if (sum(G, na.rm = T) == 0) {
        return (list(as.matrix(0), 0, 1e6))
    }
    Y <- as.matrix(Y)
    K <- diag(1, nrow = ncol(G))
    L <- chol(K)
    W <- G %*% L
    W_sparse <- as(W, "dgCMatrix")
    svd_mat <- sparsesvd::sparsesvd(W_sparse)

    U <- svd_mat$u
    S <- (svd_mat$d)^2

    UtY <- t(U) %*% Y

    # Y_UUtY = Y - UUtY
    Y_UUtY <- Y - U %*% UtY

    opt <- optim(par = 1, fn = calc_log_lik, S = S, UtY = UtY, Y_UUtY = Y_UUtY,
                 method = c("Brent"), lower = 0, upper = 1e6, control = list(fnscale = -1))
    opt_delta <- opt$par

    tr_GtG <- sum(diag(t(G) %*% G %*% K))

    post_beta <- calc_post_beta(K, G, opt_delta, S, UtY, U)

    return (list(post_beta, tr_GtG, opt_delta))
}


calc_gene_effect_size <- function(G, lof_ncol, post_beta, Y) {
    beta_lof <- post_beta[1:lof_ncol]
    beta_lof <- abs(beta_lof)
    lof_prs <- G[, 1:lof_ncol, drop = F] %*% beta_lof
    lof_prs_norm <- (lof_prs - mean(lof_prs)) / sd(lof_prs)
    m1 <- lm(Y ~ lof_prs_norm[, 1])

    return(m1$coefficients[2])
}

mom_estimator_marginal <- function(G, y) {
    n <- length(y)
    G <- as(G, "dgCMatrix")
    G[is.na(G)] <- 0
    # tr(G Sigma G^T G Sigma G^T) = sum((G Sigma G^T )^2) = sum(G^T G Sigma)^2
    Sigma <- diag(1, nrow = ncol(G))

    system.time({
        t1 <- sum((crossprod(G) %*% Sigma)^2)
    })
    t2 <- sum(diag(crossprod(G) %*% Sigma))
    A <- matrix(c(t1, t2, t2, n), ncol = 2)
    c1 <- as.numeric(t(y) %*% G %*% Sigma %*% t(G) %*% y)
    c2 <- sum(y^2)
    b <- matrix(c(c1, c2), ncol = 1)
    var_comp <- solve(A) %*% b
    print(var_comp)
    # h2_mom_marginal <- var_comp[1, 1] / sum(var_comp)
    return (var_comp)
}

mom_estimator_joint <- function(G1, G2, G3, y) {
    n <- length(y)

    G1 <- as(G1, "dgCMatrix")
    G1[is.na(G1)] <- 0
    G2 <- as(G2, "dgCMatrix")
    G2[is.na(G2)] <- 0
    G3 <- as(G3, "dgCMatrix")
    G3[is.na(G3)] <- 0
    Sigma1 <- diag(1, nrow = ncol(G1)) # L1 %*% t(L1)
    Sigma2 <- diag(1, nrow = ncol(G2)) # L2 %*% t(L2)
    Sigma3 <- diag(1, nrow = ncol(G3)) # L3 %*% t(L3)

    L1 <- chol(Sigma1)
    L2 <- chol(Sigma2)
    L3 <- chol(Sigma3)

    t11 <- sum((crossprod(G1) %*% Sigma1)^2)
    t22 <- sum((crossprod(G2) %*% Sigma2)^2)
    t33 <- sum((crossprod(G3) %*% Sigma3)^2)

    t12 <- sum((t(G1 %*% L1) %*% (G2 %*% L2))^2)
    t13 <- sum((t(G1 %*% L1) %*% (G3 %*% L3))^2)
    t23 <- sum((t(G2 %*% L2) %*% (G3 %*% L3))^2)

    t14 <- sum(diag(t(G1) %*% G1 %*% Sigma1))
    t24 <- sum(diag(t(G2) %*% G2 %*% Sigma2))
    t34 <- sum(diag(t(G3) %*% G3 %*% Sigma3))

    A <- matrix(c(t11, t12, t13, t14, t12, t22, t23, t24, t13, t23, t33, t34, t14, t24, t34, n), ncol = 4)

    c1 <- as.numeric(t(y) %*% G1 %*% Sigma1 %*% t(G1) %*% y)
    c2 <- as.numeric(t(y) %*% G2 %*% Sigma2 %*% t(G2) %*% y)
    c3 <- as.numeric(t(y) %*% G3 %*% Sigma3 %*% t(G3) %*% y)
    c4 <- sum(y^2)
    b <- matrix(c(c1, c2, c3, c4), ncol = 1)

    var_comp <- solve(A) %*% b
    print(var_comp)
    #  h2_mom_joint <- c(var_comp[1, 1], var_comp[2, 1], var_comp[3, 1]) / sum(var_comp)
    return (var_comp)
}

# calculate_single_blup <- function(G, delta, Sigma, y) {
#     # Henderson Mixed Model Equation: beta = (G^T G + delta * Sigma^-1)^-1 * G^T y
#     G <- as(G, "dgCMatrix")
#     beta <- solve(t(G) %*% G + delta * solve(Sigma)) %*% t(G) %*% y
#     return (beta)
# }

calculate_joint_blup <- function(G1, G2, G3, tau1, tau2, tau3, Sigma1, Sigma2, Sigma3, psi, y) {
    G1 <- as(G1, "dgCMatrix")
    G2 <- as(G2, "dgCMatrix")
    G3 <- as(G3, "dgCMatrix")
    G <- cbind(G1, G2, G3)
    G[is.na(G)] <- 0

    Sigma <- bdiag(tau1 * Sigma1, tau2 * Sigma2, tau3 * Sigma3)
    beta <- solve(t(G) %*% G / psi + solve(Sigma)) %*% t(G) %*% y / psi
    return (beta)
}
