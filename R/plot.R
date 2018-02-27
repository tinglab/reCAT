require(ggplot2)

my_theme <- function()
{
  theme(
    plot.title = element_text(size = 40),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black"),
    axis.title.x = element_text(size = 40),
    axis.title.y = element_text(size = 40),
    axis.text.x = element_text(size = 40, color = "black"),
    axis.text.y = element_text(size = 40, color = "black"),
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 30),
    legend.key = element_rect(fill = "white"))
}


plot_bayes <- function(bayes_score, ordIndex, cls_result = NULL, cls_ord = NULL, colorbar = 0)
{
  if (colorbar == 0)
  {
    cell = bayes_score[ordIndex, ]
    
    data <- data.frame(x = rep(1:length(ordIndex), 3),
                       y = c(cell$G1.score, cell$S.score, cell$G2M.score),
                       z = factor(c(rep("G1 score", length(ordIndex)), rep("S score",length(ordIndex)), rep("G2/M score", length(ordIndex))), levels = c("G1 score", "S score", "G2/M score")))
    
    p <- ggplot(data)
    p + geom_point(aes(x = x, y = y, color = z), size = 2) +
      geom_line(aes(x = x, y = y, color = z, group = z), size = 0.5) +
      scale_color_manual(name = "", values = c("G1 score" = "#6699FF", "S score" = "#66CC00", "G2/M score" = "#FF3366"), 
                         breaks = c("G1 score", "S score", "G2/M score"),
                         guide = guide_legend(keywidth = 2, keyheight = 2)) +
      labs(title = "", x = "Time Series (t)", y = "Bayes-scores") +
      my_theme()
  }
  else
  {
    if (length(unique(cls_result)) == 3)
    {
      G1.id = which(cls_result == 1)
      S.id = which(cls_result == 2)
      G2M.id = which(cls_result == 3)
      col = rep("G1", length(ordIndex))
      col[S.id] = "S"
      col[G2M.id] = "G2/M"
      
      if (is.null(cls_ord))
      {
        cls_ord = c(1:length(ordIndex))
      }
      
      #col = col[cls_ord]
      cell = bayes_score[ordIndex, ][cls_ord, ]
      
      data <- data.frame(x = rep(1:length(ordIndex), 3),
                         y = c(cell$G1.score, cell$S.score, cell$G2M.score),
                         z = factor(c(rep("G1 score", length(ordIndex)), rep("S score",length(ordIndex)), rep("G2/M score", length(ordIndex))), levels = c("G1 score", "S score", "G2/M score")))
      data3 = data.frame(x = c(1:length(ordIndex)), y = max(bayes_score)+50, color = col)
      
      p <- ggplot(data)
      p + geom_point(aes(x = x, y = y, color = z), size = 2) +
        geom_line(aes(x = x, y = y, color = z, group = z), size = 0.5) +
        geom_segment(data = data3, aes(x = x, xend = x, y = y, yend = y+50, color = col), size = 3) +
        scale_color_manual(name = "", values = c("G1 score" = "#6699FF", "S score" = "#66CC00", "G2/M score" = "#FF3366", 
                                                 "G1" = "#6699FF", "S" = "#66CC00", "G2/M" = "#FF3366"), 
                           breaks = c("G1", "S", "G2/M", "G1 score", "S score", "G2/M score"),
                           guide = guide_legend(keywidth = 2, keyheight = 2)) +
        labs(title = "", x = "Time Series (t)", y = "Bayes-scores") +
        my_theme()
    }
    else
    {
      G0.id = which(cls_result == 1)
      G1.id = which(cls_result == 2)
      S.id = which(cls_result == 3)
      G2M.id = which(cls_result == 4)
      col = rep("G0", length(ordIndex))
      col[G1.id] = "G1"
      col[S.id] = "S"
      col[G2M.id] = "G2/M"
      
      if (is.null(cls_ord))
      {
        cls_ord = c(1:length(ordIndex))
      }
      #col = col[cls_ord]
      cell = bayes_score[ordIndex, ][cls_ord, ]
      
      data <- data.frame(x = rep(1:length(ordIndex), 3),
                         y = c(cell$G1.score, cell$S.score, cell$G2M.score),
                         z = factor(c(rep("G1 score", length(ordIndex)), rep("S score",length(ordIndex)), rep("G2/M score", length(ordIndex))), levels = c("G1 score", "S score", "G2/M score")))
      data3 = data.frame(x = c(1:length(ordIndex)), y = max(bayes_score)+50, color = col)
      
      p <- ggplot(data)
      p + geom_point(aes(x = x, y = y, color = z), size = 2) +
        geom_line(aes(x = x, y = y, color = z, group = z), size = 0.5) +
        geom_segment(data = data3, aes(x = x, xend = x, y = y, yend = y+50, color = col), size = 3) +
        scale_color_manual(name = "", values = c("G1 score" = "#6699FF", "S score" = "#66CC00", "G2/M score" = "#FF3366", 
                                                "G0" = "black", "G1" = "#6699FF", "S" = "#66CC00", "G2/M" = "#FF3366"), 
                           breaks = c("G0", "G1", "S", "G2/M", "G1 score", "S score", "G2/M score"),
                           guide = guide_legend(keywidth = 2, keyheight = 2)) +
        labs(title = "", x = "Time Series (t)", y = "Bayes-scores") +
        my_theme()
    }
  }
}


plot_mean <- function(mean_score, ordIndex, cls_result = NULL, cls_ord = NULL, colorbar = 0)
{
  if (colorbar == 0)
  {
    cell = mean_score[ordIndex, ]
    
    data2 <- data.frame(x = rep(1:length(ordIndex),6), 
                        y = c(cell$G1Score, cell$G1SScore, cell$SScore, cell$G2Score, cell$G2MScore, cell$MScore),
                        z = factor(c(rep("G1 score", length(ordIndex)), rep("G1/S score", length(ordIndex)), rep("S score", length(ordIndex)), 
                                     rep("G2 score", length(ordIndex)), rep("G2/M score", length(ordIndex)), rep("M score", length(ordIndex))), 
                                   levels = c("G1 score", "G1/S score", "S score", "G2 score", "G2/M score", "M score")))
    
    p <- ggplot(data2)
    p + geom_point(aes(x = x, y = y, color = z), size = 2) +
      geom_line(aes(x = x, y = y, color = z, group = z), size = 0.5) +
      scale_color_manual(name = "", values = c("G1 score" = "#6699FF","G1/S score" = "black", "S score" = "#66CC00", 
                                               "G2 score" = "orange", "G2/M score" = "#FF3366", "M score" = "yellow"), 
                         breaks = c("G1 score", "G1/S score", "S score", "G2 score", "G2/M score", "M score"),
                         guide = guide_legend(keywidth = 2, keyheight = 2)) +
      labs(title = "", x = "Time Series (t)", y = "Mean-scores") +
      my_theme()
  }
  else
  {
    if (length(unique(cls_result)) == 3)
    {
      G1.id = which(cls_result == 1)
      S.id = which(cls_result == 2)
      G2M.id = which(cls_result == 3)
      col = rep("G1", length(ordIndex))
      col[S.id] = "S"
      col[G2M.id] = "G2/M"
      
      if (is.null(cls_ord))
      {
        cls_ord = c(1:length(ordIndex))
      }
      #col = col[cls_ord]
      cell = mean_score[ordIndex, ][cls_ord, ]
      
      data2 <- data.frame(x = rep(1:length(ordIndex),6), 
                          y = c(cell$G1Score, cell$G1SScore, cell$SScore, cell$G2Score, cell$G2MScore, cell$MScore),
                          z = factor(c(rep("G1 score", length(ordIndex)), rep("G1/S score", length(ordIndex)), rep("S score", length(ordIndex)), 
                                       rep("G2 score", length(ordIndex)), rep("G2/M score", length(ordIndex)), rep("M score", length(ordIndex))), 
                                     levels = c("G1 score", "G1/S score", "S score", "G2 score", "G2/M score", "M score")))
      
      data3 = data.frame(x = c(1:length(ordIndex)), y = max(mean_score)+1, color = col)
      
      p <- ggplot(data2)
      p + geom_point(aes(x = x, y = y, color = z), size = 2) +
        geom_line(aes(x = x, y = y, color = z, group = z), size = 0.5) +
        geom_segment(data = data3, aes(x = x, xend = x, y = y, yend = y+0.5, color = color), size = 3) +
        scale_color_manual(name = "", values = c("G1 score" = "#6699FF","G1/S score" = "black", "S score" = "#66CC00", 
                                                 "G2 score" = "orange", "G2/M score" = "#FF3366", "M score" = "yellow",
                                                 "G1" = "#6699FF", "S" = "#66CC00", "G2/M" = "#FF3366"), 
                           breaks = c("G1", "S", "G2/M", "G1 score", "G1/S score", "S score", "G2 score", "G2/M score", "M score"),
                           guide = guide_legend(keywidth = 2, keyheight = 2)) +
        labs(title = "", x = "Time Series (t)", y = "Mean-scores") +
        my_theme()
    }
    else
    {
      G0.id = which(cls_result == 1)
      G1.id = which(cls_result == 2)
      S.id = which(cls_result == 3)
      G2M.id = which(cls_result == 4)
      col = rep("G0", length(ordIndex))
      col[G1.id] = "G1"
      col[S.id] = "S"
      col[G2M.id] = "G2/M"
      
      
      if (is.null(cls_ord))
      {
        cls_ord = c(1:length(ordIndex))
      }
      #col = col[cls_ord]
      cell = mean_score[ordIndex, ][cls_ord, ]
      
      data2 <- data.frame(x = rep(1:length(ordIndex),6), 
                          y = c(cell$G1Score, cell$G1SScore, cell$SScore, cell$G2Score, cell$G2MScore, cell$MScore),
                          z = factor(c(rep("G1 score", length(ordIndex)), rep("G1/S score", length(ordIndex)), rep("S score", length(ordIndex)), 
                                       rep("G2 score", length(ordIndex)), rep("G2/M score", length(ordIndex)), rep("M score", length(ordIndex))), 
                                     levels = c("G1 score", "G1/S score", "S score", "G2 score", "G2/M score", "M score")))
      
      data3 = data.frame(x = c(1:length(ordIndex)), y = max(mean_score)+1, color = col)
      
      p <- ggplot(data2)
      p + geom_point(aes(x = x, y = y, color = z), size = 2) +
        geom_line(aes(x = x, y = y, color = z, group = z), size = 0.5) +
        geom_segment(data = data3, aes(x = x, xend = x, y = y, yend = y+0.5, color = color), size = 3) +
        scale_color_manual(name = "", values = c("G1 score" = "#6699FF","G1/S score" = "black", "S score" = "#66CC00", 
                                                 "G2 score" = "orange", "G2/M score" = "#FF3366", "M score" = "yellow",
                                                "G0" = "black", "G1" = "#6699FF", "S" = "#66CC00", "G2/M" = "#FF3366"), 
                           breaks = c("G0", "G1", "S", "G2/M", "G1 score", "G1/S score", "S score", "G2 score", "G2/M score", "M score"),
                           guide = guide_legend(keywidth = 2, keyheight = 2)) +
        labs(title = "", x = "Time Series (t)", y = "Mean-scores") +
        my_theme()
    }
  }
}
