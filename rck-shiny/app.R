library(shiny)
library(ggplot2)
library(deSolve)
library(dplyr)
library(splines)

# --- RCK Modell Logik in R ---

# Funktion für dc/dk (Sattelpfad Berechnung)
rck_dcdk <- function(k, c, params) {
  # k wird hier als die "Zeitvariable" in der ODE betrachtet
  # c ist die Zustandsvariable y[1]
  
  rho <- params$rho
  alpha <- params$alpha
  theta <- params$theta
  xi <- params$xi
  delta <- params$delta
  phi <- params$phi
  
  k_safe <- max(k, 1e-6)
  c_safe <- max(c[1], 1e-6)
  
  # dc/dt
  dc_dt <- (c_safe / rho) * (alpha * (k_safe^(alpha - 1)) - theta - (xi + delta) - rho * phi)
  
  # dk/dt
  dk_dt <- k_safe^alpha - c_safe - (phi + xi + delta) * k_safe
  
  # Verhindere Division durch Null nahe des k_dot=0 Lokus
  if (abs(dk_dt) < 1e-9) {
    # Hohe Steigung, um Richtung beizubehalten
    dcdk <- 1e6 * sign(dc_dt) 
  } else {
    dcdk <- dc_dt / dk_dt
  }
  
  return(list(dcdk))
}


# --- UI Definition ---
ui <- fluidPage(
  # Custom CSS für Ästhetik
  tags$head(
    tags$style(HTML("
            @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;700&display=swap');
            body { 
                font-family: 'Inter', sans-serif; 
                background-color: #f7f7f7; 
            }
            .sidebar {
                padding: 20px;
                background-color: #ffffff;
                border-radius = 12px;
                box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
            }
            .main-content {
                padding: 20px;
            }
            h2, h4 {
                color: #004d99;
                font-weight: 700;
                margin-bottom: 15px;
            }
            .plot-container {
                margin-bottom: 30px;
                padding: 15px;
                border: 1px solid #e0e0e0;
                border-radius: 8px;
                background-color: #fff;
                box-shadow: 0 2px 5px rgba(0, 0, 0, 0.05);
            }
            .ss-output {
                background-color: #e6f0ff;
                border: 1px solid #cce0ff;
                border-radius: 8px;
                padding: 15px;
                margin-top: 20px;
            }
        "))
  ),
  
  titlePanel(
    div(
      img(src = "https://placehold.co/40x40/004d99/ffffff?text=RCK", style = "display: inline; margin-right: 10px; border-radius: 5px;"),
      "Interaktives RCK Wachstumsmodell (R-Shiny)"
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      id = "sidebar",
      class = "sidebar",
      h4("Modellparameter pro effektive Einheit"),
      sliderInput("rho", label = HTML("Rel. Risikoaversion (&rho;):"), min = 1.0, max = 5.0, value = 2.0, step = 0.1),
      sliderInput("alpha", label = HTML("Kapitalanteil (&alpha;):"), min = 0.1, max = 0.9, value = 0.3, step = 0.01),
      sliderInput("theta", label = HTML("Zeitpräferenzrate (&theta;):"), min = 0.01, max = 0.1, value = 0.02, step = 0.005),
      sliderInput("xi", label = HTML("Bevölkerungswachstum (&xi;):"), min = 0.0, max = 0.05, value = 0.01, step = 0.005),
      sliderInput("delta", label = HTML("Abschreibungsrate (&delta;):"), min = 0.01, max = 0.15, value = 0.08, step = 0.005),
      sliderInput("phi", label = HTML("Produktivitätswachstum (&phi;):"), min = 0.0, max = 0.05, value = 0.03, step = 0.005),
      
      tags$hr(),
      
      h4("Simulationsparameter"),
      sliderInput("k0", label = HTML("Startkapital (k&₀;):"), min = 0.1, max = 10.0, value = 4.0, step = 0.1),
      sliderInput("t_end", label = "Zeitraum (T):", min = 50, max = 200, value = 100, step = 10)
    ),
    
    mainPanel(
      class = "main-content",
      tabsetPanel(
        tabPanel("Phasendiagramm", 
                 div(class = "plot-container", plotOutput("phasePlot", height = "600px"))
        ),
        tabPanel("Zeitdynamik", 
                 div(class = "plot-container", plotOutput("kPlot", height = "400px")),
                 div(class = "plot-container", plotOutput("cPlot", height = "400px")),
                 div(class = "plot-container", plotOutput("yPlot", height = "400px"))
        ),
        tabPanel("Steady State Werte",
                 uiOutput("ss_output")
        )
      )
    )
  )
)

# --- Server Logik ---
server <- function(input, output, session) {
  
  # Reactive: Berechnet Steady State Werte und Hilfswerte
  model_params_ss <- reactive({
    params <- list(
      rho = input$rho, alpha = input$alpha, theta = input$theta,
      xi = input$xi, delta = input$delta, phi = input$phi
    )
    
    rho <- params$rho
    alpha <- params$alpha
    theta <- params$theta
    xi <- params$xi
    delta <- params$delta
    phi <- params$phi
    
    error_msg <- NULL
    kss <- NaN
    
    # Steady State Kapital (kss)
    if (alpha != 1) {
      denominator <- theta + xi + delta + rho * phi
      if (denominator > 0) {
        kss_base <- alpha / denominator
        if (kss_base > 0) {
          kss <- kss_base^(1 / (1 - alpha))
        }
      }
    }
    
    # Weitere Steady States
    css <- NaN
    yss <- NaN
    
    if (!is.nan(kss) && kss > 0) {
      css <- kss^alpha - (xi + delta + phi) * kss
      yss <- kss^alpha
      if (css < 0) {
        error_msg <- "Steady State Konsum (c̄) ist negativ (c̄ < 0)."
        kss <- NaN; css <- NaN; yss <- NaN
      }
    } else {
      error_msg <- "Steady State Kapital (k̄) ist nicht positiv."
    }
    
    # Maximales Kapital k_max
    maintenance_cost <- phi + xi + delta
    if (maintenance_cost == 0) {
      kmax <- Inf
    } else {
      kmax <- (1 / maintenance_cost)^(1 / (1 - alpha))
    }
    
    return(list(
      kss = kss, css = css, yss = yss, 
      kmax = kmax, params = params, error = error_msg
    ))
  })
  
  # Reactive: Berechnet den Sattelpfad c(k) mit Spline-Interpolation
  saddle_path_func <- reactive({
    ss_data <- model_params_ss()
    if (!is.null(ss_data$error)) return(NULL)
    
    kss <- ss_data$kss
    css <- ss_data$css
    kmax <- ss_data$kmax
    params <- ss_data$params
    
    eps <- 1e-4 # Epsilon Störung
    
    # Wrapper für die ODE-Funktion für deSolve::ode
    ode_func_wrapper <- function(k, c, parms) {
      rck_dcdk(k, c, parms)
    }
    
    # 1. Pfad unterhalb kss (k von kss nach 0, Rückwärtsintegration)
    # Mehr Punkte für höhere Auflösung
    k_seq_below <- seq(kss, 0.01, length.out = 150) 
    sol_below <- ode(
      y = c(c = css - eps), 
      times = k_seq_below,  
      func = ode_func_wrapper, 
      parms = params
    )
    
    # 2. Pfad oberhalb kss (k von kss nach kmax, Vorwärtsintegration)
    # Mehr Punkte für höhere Auflösung
    k_seq_above <- seq(kss, kmax * 0.99, length.out = 150)
    sol_above <- ode(
      y = c(c = css + eps), 
      times = k_seq_above,  
      func = ode_func_wrapper, 
      parms = params
    )
    
    # Kombiniere und bereinige Daten
    df_below <- as.data.frame(sol_below) %>% 
      rename(k = time, c = c) %>%
      filter(k > 1e-4, c > 1e-4) # Filtern nach positiven Werten
    
    df_above <- as.data.frame(sol_above) %>% 
      rename(k = time, c = c) %>%
      filter(k > 1e-4, c > 1e-4)
    
    # Kombiniere, sortiere, entferne Duplikate und NAs
    df_combined <- bind_rows(df_below, df_above) %>%
      arrange(k) %>%
      distinct(k, .keep_all = TRUE) %>%
      na.omit() 
    
    if (nrow(df_combined) < 10) return(NULL)
    
    # Spline erstellen (Methode "hyman" ist oft robuster und formtreuer)
    sp <- tryCatch({
      splinefun(df_combined$k, df_combined$c, method = "hyman")
    }, error = function(e) {
      # Bei Spline-Fehler NULL zurückgeben
      NULL 
    })
    
    return(sp)
  })
  
  # Reactive: Berechnet die Zeit-Dynamik k(t), c(t), y(t)
  dynamics_data <- reactive({
    sp <- saddle_path_func()
    ss_data <- model_params_ss()
    
    if (is.null(sp) || !is.null(ss_data$error)) return(NULL)
    
    k0 <- input$k0
    t_end <- input$t_end
    params <- ss_data$params
    
    # Dynamische ODE-Funktion für dk/dt = f(t, k)
    k_dot_func <- function(t, k, parms) {
      k_val <- k[1]
      if (k_val <= 1e-6) return(list(0))
      
      # Hole c vom Sattelpfad
      c_val <- tryCatch({
        sp(k_val)
      }, error = function(e) {
        # Wenn k_val außerhalb des spline-Bereichs ist
        return(0)
      })
      
      if (c_val < 0) return(list(0))
      
      # dk/dt = k^alpha - c - (n + g + delta)k
      dk <- k_val^parms$alpha - c_val - (parms$phi + parms$xi + parms$delta) * k_val
      
      return(list(dk))
    }
    
    # Löse k(t)
    t_seq <- seq(0, t_end, length.out = 300)
    sol <- ode(
      y = c(k = k0), 
      times = t_seq, 
      func = k_dot_func, 
      parms = params
    )
    
    df <- as.data.frame(sol) %>% rename(k = k)
    
    # Berechne c(t) und y(t)
    df <- df %>%
      mutate(
        c = sp(k),
        y = k^params$alpha
      ) %>%
      filter(c >= 0, k > 0) # Konsum und Kapital müssen positiv sein
    
    return(df)
  })
  
  # OUTPUT: Steady State Werte
  output$ss_output <- renderUI({
    ss_data <- model_params_ss()
    if (!is.null(ss_data$error)) {
      return(tags$p(class = "text-danger", paste("FEHLER:", ss_data$error)))
    }
    
    kss <- ss_data$kss
    css <- ss_data$css
    yss <- ss_data$yss
    kmax <- ss_data$kmax
    
    div(
      class = "ss-output",
      h4("Stabile Steady-State Werte (pro effektive Einheit):"),
      tags$ul(
        tags$li(HTML(sprintf("<b>Kapital (k̄):</b> %.4f", kss))),
        tags$li(HTML(sprintf("<b>Konsum (c̄):</b> %.4f", css))),
        tags$li(HTML(sprintf("<b>Output (ȳ):</b> %.4f", yss)))
      ),
      tags$p(HTML(sprintf("Maximales Kapital $k_{max}$: %.4f", kmax)))
    )
  })
  
  # OUTPUT: Phasendiagramm
  output$phasePlot <- renderPlot({
    ss_data <- model_params_ss()
    dyn_data <- dynamics_data()
    
    if (!is.null(ss_data$error)) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Warten auf gültige Parameter oder Fehler"))
    }
    
    kss <- ss_data$kss
    css <- ss_data$css
    kmax <- ss_data$kmax
    params <- ss_data$params
    
    k_grid <- seq(0.01, kmax * 1.05, length.out = 300)
    # k_dot=0 Lokus
    c_k0 <- k_grid^params$alpha - (params$phi + params$xi + params$delta) * k_grid
    
    plot_df <- data.frame(k = k_grid, c_k0 = c_k0)
    
    # Datenrahmen für Steady State Punkt
    ss_point_df <- data.frame(k = kss, c = css)
    
    p <- ggplot(plot_df, aes(x = k)) +
      # k_dot=0 Lokus
      geom_line(aes(y = c_k0, color = "k_dot=0"), linewidth = 1.5) +
      # c_dot=0 Lokus
      geom_vline(xintercept = kss, linetype = "dashed", color = "red", linewidth = 1) +
      
      # Steady State (Explizite Daten verwenden, um Warnung zu vermeiden)
      geom_point(data = ss_point_df, aes(x = k, y = c, color = "SS"), size = 4) +
      annotate("text", x = kss, y = css, label = "SS", hjust = -0.1, vjust = 0.5)
    
    # Sattelpfad
    sp <- saddle_path_func()
    if (!is.null(sp)) {
      sp_df <- data.frame(k = k_grid, c = sp(k_grid)) %>%
        filter(c >= 0, c < max(c_k0) * 1.5)
      p <- p + geom_line(data = sp_df, aes(x = k, y = c, color = "Sattelpfad"), linewidth = 2)
    }
    
    # Dynamik Trajektorie
    if (!is.null(dyn_data) && nrow(dyn_data) > 1) {
      # Datenrahmen für den Startpunkt der Trajektorie
      start_point_df <- dyn_data %>% slice(1)
      
      p <- p + 
        geom_path(data = dyn_data, aes(x = k, y = c, color = "Trajektorie"), linetype = "dotted", linewidth = 1) +
        geom_point(data = start_point_df, aes(x = k, y = c, color = "Start"), size = 4)
    }
    
    # Plot-Einstellungen
    p <- p +
      scale_color_manual(
        name = "",
        values = c("k_dot=0" = "blue", "SS" = "black", "Sattelpfad" = "darkgreen", "Trajektorie" = "orange", "Start" = "orange"),
        labels = c("k_dot=0" = expression(dot(k)==0), "SS" = "Steady State", "Sattelpfad" = "Sattelpfad", "Trajektorie" = "Trajektorie", "Start" = "Startpunkt")
      ) +
      labs(
        title = "Phasendiagramm und Sattelpfad",
        x = expression("Kapital pro effektive Einheit " (k)),
        y = expression("Konsum pro effektive Einheit " (c))
      ) +
      coord_cartesian(xlim = c(0, kmax * 1.05), ylim = c(0, max(c_k0) * 1.2)) +
      theme_minimal() +
      theme(legend.position = "bottom", 
            plot.title = element_text(hjust = 0.5, color = "#004d99", face = "bold"))
    
    return(p)
  })
  
  # Funktion zur Generierung von Zeitreihen-Plots
  plot_time_series <- function(var_name, ss_label, title, color) {
    dyn_data <- dynamics_data()
    ss_data <- model_params_ss()
    
    if (is.null(dyn_data) || !is.null(ss_data$error)) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste0(title, ": Warten auf Daten")))
    }
    
    ss_val <- ss_data[[paste0(var_name, "ss")]]
    
    p <- ggplot(dyn_data, aes(x = time, y = .data[[var_name]])) +
      geom_line(linewidth = 1.5, color = color) +
      geom_hline(yintercept = ss_val, linetype = "dashed", color = "black", linewidth = 1, show.legend = TRUE) +
      labs(
        title = title,
        x = "Zeit (t)",
        y = paste0(ss_label, " pro effektive Einheit")
      ) +
      annotate("text", x = max(dyn_data$time) * 0.95, y = ss_val, label = ss_label, hjust = 1, vjust = -0.5) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    return(p)
  }
  
  # OUTPUT: Kapital (k) Plot
  output$kPlot <- renderPlot({
    plot_time_series("k", expression(bar(k)), "Kapitalentwicklung k(t)", "#004d99")
  })
  
  # OUTPUT: Konsum (c) Plot
  output$cPlot <- renderPlot({
    plot_time_series("c", expression(bar(c)), "Konsumentwicklung c(t)", "#28b463")
  })
  
  # OUTPUT: Output (y) Plot
  output$yPlot <- renderPlot({
    plot_time_series("y", expression(bar(y)), "Outputentwicklung y(t)", "#e67e22")
  })
  
}

shinyApp(ui = ui, server = server)