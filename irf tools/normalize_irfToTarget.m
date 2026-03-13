function [IRF_final,IRF_final_boot]=normalize_irfToTarget(IRF, IRF_boot,target_bps,shock_idx,rate_idx)

impatto_iniziale_tasso = IRF(rate_idx, shock_idx, 1);

norm_factor_debug="norm_factor"
norm_factor = target_bps / impatto_iniziale_tasso

IRF_final = IRF * norm_factor;
IRF_final_boot = IRF_boot * norm_factor;