# Falls Pakete fehlen: install.packages(c("shiny", "dplyr", "ggplot2", "tidyr"))
library(shiny)
library(dplyr)
library(ggplot2)
library(tidyr)

# --- 1. User Interface (UI) ---
ui <- fluidPage(
  
  # Titel der App
  titlePanel("Interaktives Solow-Wachstumsmodell"),
  
  # Seitenleisten-Layout: Links die Slider, rechts der Plot
  sidebarLayout(
    sidebarPanel(
      h4("Parameter anpassen"),
      sliderInput("alpha", "Kapitalproduktionselastizität (alpha):", 
                  min = 0.1, max = 0.9, value = 0.4, step = 0.05),
      sliderInput("delta", "Abschreibungsrate (delta):", 
                  min = 0.01, max = 0.3, value = 0.1, step = 0.01),
      sliderInput("s", "Sparquote (s):", 
                  min = 0.05, max = 0.95, value = 0.2, step = 0.05),
      sliderInput("n", "Bevölkerungswachstum (n):", 
                  min = -0.05, max = 0.1, value = 0.0, step = 0.01),
      sliderInput("K1", "Startkapital (K[1]):", 
                  min = 0.1, max = 15, value = 2, step = 0.1),
      hr(),
      helpText("Verändern Sie die Schieberegler, um zu sehen, wie sich das Kapital (k), die Produktion (y) und der Konsum (c) pro Kopf über die Zeit entwickeln.")
    ),
    
    # Hauptbereich für den Plot
    mainPanel(
      plotOutput("solowPlot", height = "700px")
    )
  )
)

# --- 2. Server-Logik ---
server <- function(input, output) {
  
  # Reactive Plot-Generierung
  output$solowPlot <- renderPlot({
    
    # 1. Parameter aus den UI-Slidern auslesen
    alpha <- input$alpha
    delta <- input$delta
    s <- input$s
    n <- input$n
    K_start <- input$K1
    
    L_0 <- 1
    T_steps <- 101
    t <- 0:(T_steps - 1)
    
    # Bevölkerung berechnen
    L <- L_0 * (1 + n)^t
    
    # 2. Kapital-Iteration
    K <- numeric(T_steps)
    K[1] <- K_start 
    
    for (i in 2:T_steps) {
      K[i] <- (1 - delta) * K[i-1] + s * (K[i-1]^alpha) * (L[i-1]^(1-alpha))
    }
    
    # 3. Steady States berechnen
    # Beachten Sie: delta + n im Nenner, um Bevölkerungswachstum zu berücksichtigen
    k_ss <- (s / (delta + n))^(1 / (1 - alpha))
    y_ss <- k_ss^alpha
    c_ss <- (1 - s) * y_ss
    
    # 4. Datenaufbereitung für Facets
    df_facets <- tibble(t = t, K = K, L = L) %>%
      mutate(
        k = K / L,
        y = k^alpha,
        c = (1 - s) * y
      ) %>%
      select(t, k, y, c) %>%
      pivot_longer(cols = -t, names_to = "Variable", values_to = "Wert") %>%
      mutate(Variable = factor(Variable, levels = c("k", "y", "c"), 
                               labels = c("Kapital pro Kopf (k)", "Produktion (y)", "Konsum (c)")))
    
    # Hilfs-Dataframe für Steady-State-Linien und Labels
    df_ss <- tibble(
      Variable = factor(c("Kapital pro Kopf (k)", "Produktion (y)", "Konsum (c)"), 
                        levels = c("Kapital pro Kopf (k)", "Produktion (y)", "Konsum (c)")),
      ss_wert = c(k_ss, y_ss, c_ss),
      label_text = paste("Steady State:", round(c(k_ss, y_ss, c_ss), 2))
    )
    
    # 5. Plot mit facet_wrap zeichnen
    ggplot(df_facets, aes(x = t, y = Wert)) +
      geom_line(color = "steelblue", linewidth = 1.2) +
      geom_hline(data = df_ss, aes(yintercept = ss_wert), linetype = "dashed", color = "firebrick", linewidth = 0.8) +
      geom_text(data = df_ss, aes(x = T_steps * 0.85, y = ss_wert, label = label_text), 
                vjust = 1.5, color = "firebrick", fontface = "bold", size = 4.5) +
      facet_wrap(~ Variable, scales = "free_y", ncol = 1) +
      labs(
        title = "Dynamik im Solow-Wachstumsmodell",
        subtitle = "Beobachten Sie die Konvergenz zu den Steady States in Echtzeit",
        x = "Zeit (t)",
        y = "Wert pro Kopf"
      ) +
      theme_minimal() +
      theme(
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 12)
      )
  })
}

# --- 3. App starten ---
shinyApp(ui = ui, server = server)