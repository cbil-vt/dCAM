# ----
source('./example/config.R')
source('./example/dev_prep.R')

sim_dat_file <- file.path(PROJECT_DIR, "/data_ref/dataset_realistic_simulation_K3_n10_noise0.10.rds")
dataset <- load_rds(sim_dat_file)

# ----
fit <- run_dcam(dataset,
                penalty      = PENALTY,
                K            = K_USE,  K_range = K_RANGE,
                lambda_grid  = LAMBDA_GRID,
                lambda_rule  = LAMBDA_RULE,
                lambda_fixed = LAMBDA_FIXED,
                params       = PENALTY_PARAMS,
                cam_params   = CAM3_PARAMS,   cam_seed = CAM_SEED,
                cam_dup_cor  = CAM_DUP_COR_THRES,
                bic_zero_rule = BIC_ZERO_RULE,
                bic_zero_frac = BIC_ZERO_FRAC,
                bic_zero_abs  = BIC_ZERO_ABS,
                reweight     = REWEIGHT,
                freeze_A     = FREEZE_A,
                max_iter     = MAX_ITER, tol = TOL,
                max_fpr      = MAX_FPR)

# ----
q <- estimation_quality(fit, dataset$truth)
cat(sprintf("\n  agreement with the truth:  signatures=%.3f  proportions=%.3f  change=%.3f\n",
            q$S_cor, q$A_cor, q$delta_cor))

up <- compute_deg_metrics(fit$delta_S, dataset$truth$is_up,   MAX_FPR, "up")
dn <- compute_deg_metrics(fit$delta_S, dataset$truth$is_down, MAX_FPR, "down")
cat(sprintf("  detection (pAUC):          up=%.4f  down=%.4f\n",
            up$pooled_pauc_norm, dn$pooled_pauc_norm))


