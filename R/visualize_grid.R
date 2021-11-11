library(ggplot2)

#' Creates a graphic representation of the input adaptive grid search
#' @name visualize_grid
#' @param grid_search \code{List} object returned from \code{optMLE_grid()} with option \code{return_full_grid = TRUE}.
#' @return
#' @export
visualize_grid <- function(grid_search) {

  opt_long <- grid_search$all_opt[, c("grid", "audit_step", grep("n", colnames(grid_search$all_opt), value = TRUE))]

  num_strat <- length(grep("n", colnames(grid_search$all_opt), value = TRUE))
  opt_long <- data.frame(grid = rep(opt_long$grid, times = num_strat),
                       audit_step = rep(opt_long$audit_step, times = num_strat),
                       stack(x = opt_long[, -c(1:2)]))

  opt_long$grid_lab <- with(opt_long, paste0("Grid ", grid, " (", audit_step, " people)"))
  opt_long$optimal <- TRUE
  opt_long$grid <- factor(unique(opt_long$grid)[order(unique(opt_long$grid))],
                          labels = unique(opt_long$grid_lab)[order(unique(opt_long$grid))])

  library(ggplot2)

  g$full_grid_search %>% dplyr::select(-dplyr::starts_with("pi")) %>%
    dplyr::group_by(grid) %>%
    dplyr::mutate(min_var = min(Vbeta), is_min_des = Vbeta == min_var) %>%
    tidyr::gather("strat", "n", -c(1, 6:8)) %>%
    dplyr::group_by(grid, min_var, strat, n) %>%
    dplyr::summarize(is_min_des = max(is_min_des)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(is_min_des = is_min_des == 1,
                  strat = factor(strat, levels = c("n00", "n01", "n10", "n11"),
                                 labels = c("(Y*=0,\nX*=0)", "(Y*=0,\nX*=1)", "(Y*=1,\nX*=0)", "(Y*=1,\nX*=1)")),
                  grid = factor(grid, labels = c("1 (40 people)", "2 (20 people)", "3 (10 people)", "4 (5 people)", "5 (1 person)"))) %>%
    ggplot(aes(x = strat, y = n, col = factor(grid))) +
    geom_point(shape = 3, size = 2, position = position_dodge(width = 0.75), stroke = 1.1) +
    geom_line(position = position_dodge(width = 0.75), size = 1) +
    geom_point(data = opt_long, aes(x = strat, y = n, group = factor(grid)), color = "gold", size = 4, shape = 17, position = position_dodge(width = 0.75)) +
    theme_bw(base_size = 16) +
    scale_color_manual(values = park_palette("voyageurs")[1:5], name = "Grid (scale):") + #scale_color_discrete(name = "Grid:") +
    theme(legend.position = c(0.25, 0.91), legend.background = element_rect(fill = alpha("white", 0))) + coord_flip() +
    #theme(legend.position = "right") + coord_flip() +
    xlab("Phase I stratum") + ylab("Stratum sample size") + ylim(c(50, 250)) +
    guides(col=guide_legend(nrow=2, byrow=TRUE))
}
