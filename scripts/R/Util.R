.gg.plot <- function(...) {
    ggplot(...) +
        theme_classic() +
        theme(axis.title = element_text(size=8)) +
        theme(axis.text = element_text(size=6)) +
        theme(legend.spacing = unit(.1, "lines"),
              legend.key.size = unit(.5, "lines"),
              legend.text = element_text(size=5),
              legend.title = element_text(size=6),
              panel.background = element_rect(fill='transparent'),
              plot.background = element_rect(fill='transparent', color=NA),
              legend.background = element_rect(fill='transparent', size=0.05),
              legend.box.background = element_rect(fill='transparent'))
}