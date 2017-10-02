args = commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("At least two arguments must be supplied (input file).n", call. = FALSE)
}

p_t_df = read.table(args[1], header = FALSE)
p_t = p_t_df[['V1']]

library(qvalue)
l <- lfdr(p_t, pi0.method = "bootstrap")

write.table(data.frame(l), file = args[2], row.names = FALSE)
