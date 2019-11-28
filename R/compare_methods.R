# Helper functions to compare
unlog <- function(x){sign(x)*2^(abs(x))}
logfc <- function(x){sign(x)*log(abs(x), base = 2)}


linear_model <- readRDS("processed/tt_complex.RDS")
select_samples <- readRDS("processed/juanjo_res.RDS")

linear_model <- linear_model
pdf("plots/compare_methods.pdf")
for (i in seq_along(linear_model)) {
    lm_values <- linear_model[[i]]$logFC
    names(lm_values) <- rownames(linear_model[[i]])
    s_samples <- logfc(select_samples[[1]][, i])
    contrasting <- gsub("^fc_", "", colnames(select_samples[[1]]))
    lm_values <- lm_values[names(s_samples)]
    cc <- cor.test(lm_values, s_samples)
    plot(lm_values, s_samples, pch = ".",
         main = paste0("logFC\n", contrasting[i]),
         xlab = "Linear model", ylab = "Selecting just the samples",
        sub = paste("cor r =", round(cc$estimate, 5)))
}
dev.off()

plot(linear_model[[1]]$AveExpr, linear_model[[1]]$P.Value, pch = ".")
