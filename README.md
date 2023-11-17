# MyPackage

## TODO
1. Extract or compute the complete information matrix, this is useful for the confidence interval
   - we can implement vcov.heritMod. `vcov` is a generic function
3. compute the corresponding p-value
4. Parallelize computation for rows when method is used on a matrix (Gene expression data)
5. Think how you might create plot for vpc
      library(ggplot2)

# Example data
data <- data.frame(
  Level = factor(c("Between-Group", "Within-Group"), levels = c("Between-Group", "Within-Group")),
  Proportion = c(0.25, 0.75)
)

# Create a donut-style pie chart
ggplot(data, aes(x = "", y = Proportion, fill = Level)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  scale_fill_manual(values = c("#FFD700", "#C0C0C0")) +
  labs(title = "Variance Partitioning Coefficients",
       subtitle = "Two-Level Hierarchical Model") +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) +
  geom_text(aes(label = scales::percent(Proportion)), position = position_stack(vjust = 0.5)) +
  annotate("text", x = 0, y = 0, label = " ", size = 8)  # Add a blank circle in the middle

# Display the plot
ggsave("VPC_donut_chart.png", width = 8, height = 6)

