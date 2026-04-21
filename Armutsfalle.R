
library(ggplot2)

# 1. Parameter definieren
alpha <- 0.5      # Produktionselastizität des Kapitals
s <- 0.3          # Sparquote
n_delta <- 0.1    # Bevölkerungswachstum (n) + Abschreibungsrate (delta)
c_bar <- 0.5      # Autonomer Konsum (Subsistenzniveau)

# 2. Daten für den Plot generieren
k <- seq(0, 10, by = 0.05) # Kapitalstock pro Kopf (x-Achse)

# Ökonomische Funktionen
y <- k^alpha
investition <- s * (y - c_bar)  # Ersparnis/Investitionskurve verschiebt sich nach unten
abschreibung <- n_delta * k     # Gerade der Kapitalverwässerung

df <- data.frame(k, investition, abschreibung)

# 3. Analytische Berechnung der Steady States (Schnittpunkte)
# Gleichung: n_delta * k = s * k^alpha - s * c_bar
# Substituiert mit x = k^0.5 führt das zu einer quadratischen Gleichung
a <- 1
b <- -(s / n_delta)
c <- (s * c_bar) / n_delta

# pq-Formel / Mitternachtsformel
x1 <- (-b - sqrt(b^2 - 4*a*c)) / (2*a)
x2 <- (-b + sqrt(b^2 - 4*a*c)) / (2*a)

k_crit <- x1^2  # Instabiles Gleichgewicht (Armutsfalle)
k_star <- x2^2  # Stabiles Gleichgewicht

# 4. Plot mit ggplot2 erstellen
plot <- ggplot(df, aes(x = k)) +
  # Linien zeichnen
  geom_line(aes(y = investition, color = "Investition: s(f(k) - c)"), linewidth = 1.2) +
  geom_line(aes(y = abschreibung, color = "Verwässerung: (n+δ)k"), linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  
  # Steady States markieren
  geom_point(aes(x = k_crit, y = n_delta * k_crit), size = 4, color = "red") +
  geom_point(aes(x = k_star, y = n_delta * k_star), size = 4, color = "darkgreen") +
  
  # Vertikale gestrichelte Linien für die Schwellenwerte
  geom_segment(aes(x = k_crit, y = 0, xend = k_crit, yend = n_delta * k_crit), linetype = "dotted", color = "red") +
  geom_segment(aes(x = k_star, y = 0, xend = k_star, yend = n_delta * k_star), linetype = "dotted", color = "darkgreen") +
  
  # Pfeile für die Dynamik (Bewegung des Kapitalstocks)
  # Unterhalb k_crit: Kapital schrumpft (Pfeil nach links)
  annotate("segment", x = 2.5, xend = 1.5, y = 0.05, arrow = arrow(length = unit(0.3, "cm")), color = "red", linewidth = 1) +
  # Zwischen k_crit und k_star: Kapital wächst (Pfeil nach rechts)
  annotate("segment", x = 3.5, xend = 4.5, y = 0.05, arrow = arrow(length = unit(0.3, "cm")), color = "darkgreen", linewidth = 1) +
  # Oberhalb k_star: Kapital schrumpft (Pfeil nach links)
  annotate("segment", x = 8.5, xend = 7.5, y = 0.05, arrow = arrow(length = unit(0.3, "cm")), color = "red", linewidth = 1) +
  
  # Beschriftungen hinzufügen
  annotate("text", x = k_crit, y = -0.05, label = "k_crit\n(Instabil)", color = "red", fontface = "bold") +
  annotate("text", x = k_star, y = -0.05, label = "k*\n(Stabil)", color = "darkgreen", fontface = "bold") +
  annotate("text", x = 1.5, y = 0.1, label = "Armutsfalle", color = "red", fontface = "italic") +
  annotate("text", x = 4, y = 0.1, label = "Wachstum \n(Big Push Bereich)", color = "darkgreen", fontface = "italic") +
  
  # Design und Layout
  scale_color_manual(values = c("Investition: s(f(k) - c)" = "#2980b9",
                                "Verwässerung: (n+δ)k" = "#e67e22")) +
  labs(
    title = "Solow-Modell mit autonomem Konsum (Armutsfalle)",
    subtitle = "Wenn die Investition unterhalb der Verwässerung liegt, sinkt der Kapitalstock.",
    x = "Kapitalstock pro Kopf (k)",
    y = "Investition / Abschreibung",
    color = "Legende:"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

# Plot anzeigen
ggsave(
  filename = "solow_armutsfalle.png", # Name der Datei
  plot = plot,                        # Das Objekt, das gespeichert werden soll
  width = 10,                         # Breite in Zoll (Standard)
  height = 7,                         # Höhe in Zoll
  dpi = 300,                          # Auflösung (300 ist Standard für Druckqualität)
  bg = "white"                        # Hintergrundfarbe (verhindert Transparenzprobleme)
)
print(plot)

plot

