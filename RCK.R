# RCK

## Pakete laden

library(R6)
library(deSolve)
#library(pracma)
library(ggplot2)
library(dplyr)

##Code

RCKmod <- R6Class("RCKmod",
                  # ----------------------------------------------------
                  # Öffentliche Felder (Parameter und Konsumfunktion)
                  # ----------------------------------------------------
                  public = list(
                    rho = NULL,
                    alpha = NULL,
                    theta = NULL,
                    xi = NULL,
                    delta = NULL,
                    phi = NULL,
                    cFunc = NULL,
                    kmax = NULL,
                    kss = NULL,
                    css = NULL,
                    
                    # ----------------------------------------------------
                    # Konstruktor (__init__)
                    # ----------------------------------------------------
                    initialize = function(rho, alpha, theta, xi, delta, phi) {
                      # Parameter zuweisen
                      self$rho <- rho
                      self$alpha <- alpha
                      self$theta <- theta
                      self$xi <- xi
                      self$delta <- delta
                      self$phi <- phi
                      
                      # Maximale Kapitalintensität (kmax)
                      self$kmax <- (1/(self$phi + self$xi + self$delta))^(1/(1-self$alpha))
                      
                      # Steady State Kapital (kss)
                      self$kss <- (self$alpha/(self$theta + self$xi + self$delta + self$rho * self$phi))^(1/(1-self$alpha))
                      # Steady State Konsum (css)
                      self$css <- self$kss^self$alpha - (self$xi + self$delta + self$phi) * self$kss
                      
                      # Modell lösen
                      self$solve()
                    },
                    
                    # ----------------------------------------------------
                    # Cobb-Douglas Produktionsfunktion (output)
                    # ----------------------------------------------------
                    output = function(k) {
                      return(k^self$alpha)
                    },
                    
                    # ----------------------------------------------------
                    # Konsum-Differentialgleichung (dcdt)
                    # ----------------------------------------------------
                    dcdt = function(c, k) {
                      # dc/dt = c/rho * [ alpha*k^(alpha-1) - theta - (xi + delta) - rho*phi ]
                      dc <- c / self$rho * (self$alpha * k^(self$alpha - 1) - self$theta - (self$xi + self$delta) - self$rho * self$phi)
                      return(dc)
                    },
                    
                    # ----------------------------------------------------
                    # Kapital-Differentialgleichung (dkdt)
                    # ----------------------------------------------------
                    dkdt = function(c, k) {
                      # dk/dt = output(k) - c - (phi + xi + delta)*k
                      dk <- self$output(k) - c - (self$phi + self$xi + self$delta) * k
                      return(dk)
                    },
                    
                    # ----------------------------------------------------
                    # Differentialgleichung für die Zeiteliminierungsmethode (dcdk)
                    # ----------------------------------------------------
                    dcdk = function(k, y, parms) {
                      # y[1] ist c, da y die Vektorvariable ist
                      c <- y[1] 
                      
                      # Zur Vermeidung von 0/0 am Steady State:
                      # Wir verwenden einen sehr kleinen Wert für dk, falls dieser Null ist.
                      dk <- self$dkdt(c, k)
                      if (abs(dk) < 1e-10) {
                        dk <- 1e-10 # Kleinen Wert setzen, um Division durch Null zu vermeiden
                      }
                      
                      dcdk_val <- self$dcdt(c, k) / dk
                      
                      # Die deSolve::ode-Funktion benötigt einen Vektor als Rückgabe
                      return(list(dcdk_val)) 
                    },
                    
                    # ----------------------------------------------------
                    # Löser (solve)
                    # ----------------------------------------------------
                    solve = function(eps = 10^(-8), npoints = 400) {
                      
                      # K-Bereiche definieren
                      k_below <- seq(self$kss, 0.0001, length.out = npoints)
                      k_above <- seq(self$kss, self$kmax, length.out = npoints)
                      
                      # Für die Integration muss die Sequenz von k aufsteigend/absteigend sein
                      # odeint: Lösen entlang der Integrationsvariable (k)
                      # In deSolve::ode: func ist f(x, y, parms), wobei x die Integrationsvariable ist (hier k)
                      
                      # k_below: Integration von kss nach 0.0001 (Richtung umkehren für die Integration)
                      c_below_sol <- ode(
                        y = self$css - eps, 
                        times = rev(k_below), # Muss in Integrationsrichtung sein (von kss nach 0.0001)
                        func = self$dcdk, 
                        parms = NULL
                      )
                      # Die Ergebnisse in der richtigen Reihenfolge speichern (von klein nach groß)
                      c_below <- c_below_sol[ , 2] %>% rev()
                      k_below_sorted <- c_below_sol[ , 1] %>% rev()
                      
                      # k_above: Integration von kss nach kmax
                      c_above_sol <- ode(
                        y = self$css + eps, 
                        times = k_above, # Integrieren von kss nach kmax
                        func = self$dcdk, 
                        parms = NULL
                      )
                      # Ergebnisse speichern
                      c_above <- c_above_sol[ , 2]
                      k_above_sorted <- c_above_sol[ , 1]
                      
                      # Gesamte Datenpunkte
                      k <- c(k_below_sorted, k_above_sorted[-1]) # kss nicht duplizieren
                      c <- c(c_below, c_above[-1])
                      
                      # Konsumfunktion als Interpolation (pracma::interp1)
                      self$cFunc <- function(k_input) {
                        # interp1(x, y, xi, method)
                        return(interp1(k, c, k_input, method = "linear")) 
                      }
                    },
                    
                    # ----------------------------------------------------
                    # Differentialgleichung für k mit optimalem c (dkdt_opt)
                    # ----------------------------------------------------
                    dkdt_opt = function(t, y, parms) {
                      # y[1] ist k, da y die Vektorvariable ist
                      k <- y[1]
                      
                      # Hole c vom Interpolator
                      c_opt <- self$cFunc(k)
                      
                      # Berechne dk/dt
                      dk_val <- self$dkdt(c_opt, k)
                      
                      return(list(dk_val))
                    },
                    
                    # ----------------------------------------------------
                    # Kapital-Dynamik-Simulation (k_dynamics)
                    # ----------------------------------------------------
                    k_dynamics = function(k0, t) {
                      # Löse die ODE für k über die Zeit
                      k_sol <- ode(
                        y = c(k = k0), # Initialwert
                        times = t,     # Zeitpunkte
                        func = self$dkdt_opt, 
                        parms = NULL   # Keine zusätzlichen Parameter
                      )
                      
                      # Rückgabe des Vektors der k-Werte
                      return(k_sol[, "k"])
                    },
                    
                    # ----------------------------------------------------
                    # k0-Locus (k0locus)
                    # ----------------------------------------------------
                    k0locus = function(k) {
                      # Konsum, der dk/dt = 0 ergibt: c = output(k) - (phi + xi + delta)*k
                      return(self$output(k) - (self$phi + self$xi + self$delta) * k)
                    },
                    
                    # ----------------------------------------------------
                    # Phasendiagramm-Plot (phase_diagram)
                    # ----------------------------------------------------
                    phase_diagram = function(npoints = 200, arrows = FALSE, n_arrows = 5) {
                      
                      k_vals <- seq(0.01, self$kmax, length.out = npoints)
                      
                      # Berechne die Werte für die Pfade
                      c_k0locus <- self$k0locus(k_vals)
                      c_saddle <- self$cFunc(k_vals)
                      
                      # Erstelle einen Datenrahmen für ggplot2
                      plot_data <- data.frame(k = k_vals, c_k0 = c_k0locus, c_sp = c_saddle)
                      
                      p <- ggplot(plot_data, aes(x = k)) +
                        # Plot k0 locus
                        geom_line(aes(y = c_k0, color = "kdot=0 locus"), linewidth = 1) +
                        # Plot c0 locus
                        geom_vline(xintercept = self$kss, linetype = "dashed", color = "blue", linewidth = 1) +
                        geom_text(aes(x = self$kss, y = max(c_k0, c_saddle) * 0.9, label = "cdot=0 locus"), 
                                  color = "blue", angle = 90, hjust = 1.1) +
                        # Plot saddle path
                        geom_line(aes(y = c_sp, color = "Saddle path"), linewidth = 1) +
                        # Plot steady state
                        geom_point(aes(x = self$kss, y = self$css, color = "Steady state"), size = 4, shape = 8) +
                        
                        # Titel und Beschriftungen
                        labs(
                          title = "Phase Diagramm und Konsumregel\n(durch Effizienzeinheiten normalisiert)",
                          x = "k",
                          y = "c",
                          color = "Legende"
                        ) +
                        scale_color_manual(values = c("kdot=0 locus" = "darkgreen", 
                                                      "Saddle path" = "red", 
                                                      "Steady state" = "red")) +
                        theme_minimal()
                      
                      # Füge Vektorfeld-Pfeile hinzu (Quiver Plot)
                      if (arrows) {
                        k_grid <- seq(k_vals[1], k_vals[length(k_vals)], length.out = n_arrows)
                        # Erstelle y-Gitter basierend auf dem Bereich der Sattelpfad-Werte
                        c_range <- range(c_saddle, na.rm = TRUE)
                        c_grid <- seq(c_range[1], c_range[2], length.out = n_arrows)
                        
                        # Erstelle ein Gitter von k und c
                        grid_data <- expand.grid(k = k_grid, c = c_grid)
                        
                        # Berechne dc/dt und dk/dt an jedem Gitterpunkt
                        dc <- self$dcdt(grid_data$c, grid_data$k)
                        dk <- self$dkdt(grid_data$c, grid_data$k)
                        
                        # Normalisieren der Vektoren
                        M <- sqrt(dk^2 + dc^2)
                        M[M == 0] <- 1
                        dk_norm <- dk / M
                        dc_norm <- dc / M
                        
                        # Addiere die Quiver-Schicht zu ggplot (mit geom_segment)
                        p <- p + geom_segment(
                          data = grid_data,
                          aes(x = k, y = c, xend = k + dk_norm * 0.05, yend = c + dc_norm * 0.05), # Skalierung für bessere Sichtbarkeit
                          arrow = arrow(length = unit(0.01, "npc")),
                          alpha = 0.5
                        )
                      }
                      
                      print(p)
                    }
                  )
)



## Funktionen abfragen

# Parameter definieren (wie im Python-Code)
rho_val <- 2
alpha_val <- 0.3
theta_val <- 0.03
xi_val <- 0
delta_val <- 0.02
phi_val <- 0

# Modell erstellen und lösen
RCKmodExample <- RCKmod$new(
  rho = rho_val, 
  alpha = alpha_val, 
  theta = theta_val, 
  xi = xi_val, 
  delta = delta_val, 
  phi = phi_val
)

# Testen der Konsumregel
k_test <- RCKmodExample$kss / 2
c_test <- RCKmodExample$cFunc(k_test)
cat(sprintf("Konsum bei k = %.2f ist c = %.2f\n", k_test, c_test))

# Optional: Phasendiagramm plotten
# RCKmodExample$phase_diagram(arrows = TRUE, n_arrows = 12)

# ----------------------------------------------------
# Simulation der Kapital- und Konsumdynamik
# ----------------------------------------------------

# Gitter der Zeitpunkte erstellen
t <- seq(0, 100, length.out = 100)

# Startkapital
k0 <- 4

# Kapitaldynamik finden
k_sim <- RCKmodExample$k_dynamics(k0, t)

# Datenrahmen für das Plotten erstellen
sim_data <- data.frame(Time = t, k = k_sim)

## Plot 1: Kapital (k) Dynamik
p_k <- ggplot(sim_data, aes(x = Time, y = k)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_hline(yintercept = RCKmodExample$kss, linetype = "dashed", color = "black") +
  annotate("text", x = max(t) * 0.9, y = RCKmodExample$kss, label = expression(bar(k)), 
           vjust = -0.5, size = 5) +
  labs(title = "Kapitalintensität (k) über die Zeit",
       x = "Zeit",
       y = "k") +
  theme_minimal()

print(p_k)

# Konsum finden
c_sim <- RCKmodExample$cFunc(k_sim)

# Datenrahmen für das Konsum-Plot erstellen
sim_data_c <- data.frame(Time = t, c = c_sim)

## Plot 2: Konsum (c) Dynamik
p_c <- ggplot(sim_data_c, aes(x = Time, y = c)) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_hline(yintercept = RCKmodExample$css, linetype = "dashed", color = "black") +
  annotate("text", x = max(t) * 0.9, y = RCKmodExample$css, label = expression(bar(c)), 
           vjust = -0.5, size = 5) +
  labs(title = "Konsum (c) über die Zeit",
       x = "Zeit",
       y = "c") +
  theme_minimal()

print(p_c) # ist das so richtig?


## Phasendiagramm

RCKmodExample$phase_diagram(arrows = TRUE, n_arrows = 12)
