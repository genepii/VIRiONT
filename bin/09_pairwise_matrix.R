# Génération de la Heatmap avec échelle ajustée (Visualisation contrastée)
p <- ggplot(plot_data, aes(x = comp1, y = comp2, fill = count)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = count), color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("#2b83ba", "#abdda4", "#fdae61", "#d7191c"), # Bleu -> Vert -> Orange -> Rouge
    trans = "sqrt",                                        # Transformation racine carrée pour faire ressortir les faibles valeurs (0 à 80)
    breaks = c(0, 10, 50, 100, 300, 700),
    name = "Mismatches"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, face = "bold"),
    axis.text.y = element_text(size = 8, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave("matrix_comp.pdf", plot = p, width = 9, height = 9)