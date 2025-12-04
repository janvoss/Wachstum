import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.interpolate import make_interp_spline
from shiny import App, render, ui, reactive
import matplotlib.pyplot as plt
from typing import List

# --- RCK Modell Logik in Python ---

def rck_dcdk(k: float, c: float, params: dict) -> float:
    """
    Berechnet dc/dk = (dc/dt) / (dk/dt) für die Zeiteliminierungsmethode.
    """
    rho = params['rho']
    alpha = params['alpha']
    theta = params['theta']
    xi = params['xi']
    delta = params['delta']
    phi = params['phi']
    
    k_safe = max(k, 1e-6)
    c_safe = max(c, 1e-6)
    
    # dc/dt
    dc_dt = (c_safe / rho) * (alpha * (k_safe**(alpha - 1)) - theta - (xi + delta) - rho * phi)
    
    # dk/dt
    dk_dt = k_safe**alpha - c_safe - (phi + xi + delta) * k_safe
    
    # Verhindere Division durch Null nahe des k_dot=0 Lokus
    if abs(dk_dt) < 1e-9:
        # Vertikale Tangente oder Steady State
        # Wir geben einen großen Wert zurück, um die Steigung anzunähern, 
        # oder 0, wenn wir genau im SS wären (wird durch Solver meist übersprungen)
        return 1e6 * np.sign(dc_dt)

    return dc_dt / dk_dt

# --- SHINY UI (Benutzeroberfläche) ---

app_ui = ui.page_sidebar(
    
    # 1. SIDEBAR
    ui.sidebar(
        ui.h4("Modellparameter pro effektive Einheit"),
        
        ui.input_slider("rho", "Rel. Risikoaversion (ρ):", min=1.0, max=5.0, value=2.0, step=0.1),
        ui.input_slider("alpha", "Kapitalanteil (α):", min=0.1, max=0.9, value=0.3, step=0.01),
        ui.input_slider("theta", "Zeitpräferenzrate (θ):", min=0.01, max=0.1, value=0.02, step=0.005),
        ui.input_slider("xi", "Bevölkerungswachstum (ξ):", min=0.0, max=0.05, value=0.01, step=0.005),
        ui.input_slider("delta", "Abschreibungsrate (δ):", min=0.01, max=0.15, value=0.08, step=0.005),
        ui.input_slider("phi", "Produktivitätswachstum (φ):", min=0.0, max=0.05, value=0.03, step=0.005),
        
        ui.tags.hr(),
        
        ui.h4("Simulationsparameter"),
        ui.input_slider("k0", "Startkapital (k₀):", min=0.1, max=10.0, value=4.0, step=0.1),
        ui.input_slider("t_end", "Zeitraum (T):", min=50, max=200, value=100, step=10),
        
        id="sidebar"
    ),

    # 2. HAUPTINHALT
    ui.tags.head(
        ui.tags.style(
            """
            body { font-family: 'Inter', sans-serif; background-color: #f7f7f7; }
            .sidebar { padding: 20px; background-color: #ffffff; border-right: 1px solid #e0e0e0; border-radius: 12px; box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1); }
            .main-panel { padding: 20px; }
            h2, h4 { color: #004d99; font-weight: 700; margin-bottom: 15px; }
            .plot-container { margin-bottom: 30px; padding: 15px; border: 1px solid #e0e0e0; border-radius: 8px; background-color: #fff; box-shadow: 0 2px 5px rgba(0, 0, 0, 0.05); }
            .ss-output { background-color: #e6f0ff; border: 1px solid #cce0ff; border-radius: 8px; padding: 15px; margin-top: 20px; }
            """
        )
    ),

    ui.div(
        {"class": "main-panel"},
        ui.navset_tab(
            ui.nav_panel("Phasendiagramm", 
                         ui.div({"class": "plot-container"}, ui.output_plot("phasePlot", height="600px"))
            ),
            ui.nav_panel("Zeitdynamik", 
                         ui.div({"class": "plot-container"}, ui.output_plot("kPlot", height="400px")),
                         ui.div({"class": "plot-container"}, ui.output_plot("cPlot", height="400px")),
                         ui.div({"class": "plot-container"}, ui.output_plot("yPlot", height="400px"))
            ),
            ui.nav_panel("Steady State Werte",
                         ui.output_ui("ss_output")
            )
        )
    ),
    
    title=ui.tags.div(
        ui.tags.img(src="https://placehold.co/40x40/004d99/ffffff?text=RCK", style="display: inline; margin-right: 10px; border-radius: 5px;"),
        "Interaktives RCK Wachstumsmodell (Python)"
    )
)


# --- SHINY SERVER ---
def server(input, output, session):
    
    @reactive.Calc
    def model_params_ss():
        """Berechnet Steady State Werte und Hilfswerte."""
        params = {
            'rho': input.rho(), 'alpha': input.alpha(), 'theta': input.theta(),
            'xi': input.xi(), 'delta': input.delta(), 'phi': input.phi()
        }
        
        rho, alpha, theta, xi, delta, phi = params.values()

        # Steady State Kapital (kss)
        if alpha == 1:
            kss = np.nan
        else:
            denominator = theta + xi + delta + rho * phi
            if denominator <= 0:
                kss = np.nan
            else:
                kss_base = alpha / denominator
                if kss_base <= 0:
                    kss = np.nan
                else:
                    kss = kss_base**(1 / (1 - alpha))

        # Weitere Steady States
        css = np.nan
        yss = np.nan
        error_msg = None
        
        if not np.isnan(kss) and kss > 0:
            css = kss**alpha - (xi + delta + phi) * kss
            yss = kss**alpha
            if css < 0:
                error_msg = "Steady State Konsum (c̄) ist negativ (c̄ < 0)."
        else:
            error_msg = "Steady State Kapital (k̄) ist nicht positiv."

        # Maximales Kapital k_max
        maintenance_cost = phi + xi + delta
        if maintenance_cost == 0:
            kmax = np.inf
        else:
            kmax = (1 / maintenance_cost)**(1 / (1 - alpha))
        
        return {
            'kss': kss, 'css': css, 'yss': yss, 
            'kmax': kmax, 'params': params, 'error': error_msg
        }

    @reactive.Calc
    def saddle_path_func():
        """
        Berechnet den Sattelpfad c(k) mit der Zeit-Eliminierungs-Methode (Time Elimination Method).
        Wir lösen dc/dk statt dc/dt und dk/dt.
        """
        ss_data = model_params_ss()
        if ss_data['error']:
            return None
        
        kss, css, kmax, params = ss_data['kss'], ss_data['css'], ss_data['kmax'], ss_data['params']
        
        # Epsilon Störung, um dc/dk = 0/0 im Steady State zu vermeiden
        eps = 1e-4 # Kleine Störung des Konsums
        
        # Wrapper für solve_ivp: y = c, t = k
        # solve_ivp erwartet fun(t, y), also fun(k, c)
        fun_ode = lambda k, c: [rck_dcdk(k, c[0], params)]

        # 1. Pfad unterhalb kss (k von kss nach 0)
        # Wir starten leicht unterhalb von css (Sattelpfadsteigung ist positiv)
        # Wenn wir k verringern, muss c auch verringert werden, also starten wir bei css - eps.
        k_span_below = (kss, 0.01) # Rückwärtsintegration
        init_below = [css - eps]
        
        sol_below = solve_ivp(fun_ode, k_span_below, init_below, dense_output=True, 
                              method='RK45', rtol=1e-6, atol=1e-8)
        
        # 2. Pfad oberhalb kss (k von kss nach kmax)
        # Wir starten leicht oberhalb von css
        k_span_above = (kss, kmax * 0.99) # Vorwärtsintegration
        init_above = [css + eps]
        
        sol_above = solve_ivp(fun_ode, k_span_above, init_above, dense_output=True, 
                              method='RK45', rtol=1e-6, atol=1e-8)
        
        # Kombiniere Daten
        # sol.t ist k, sol.y[0] ist c
        
        # Below (umkehren, damit k aufsteigend ist)
        k_below = sol_below.t[::-1]
        c_below = sol_below.y[0][::-1]
        
        # Above
        k_above = sol_above.t
        c_above = sol_above.y[0]
        
        k_combined = np.concatenate([k_below, k_above])
        c_combined = np.concatenate([c_below, c_above])
        
        # Datenbereinigung
        valid_idx = (k_combined > 0) & (c_combined > 0) & np.isfinite(k_combined) & np.isfinite(c_combined)
        k_clean = k_combined[valid_idx]
        c_clean = c_combined[valid_idx]
        
        # Sortieren (sollte eigentlich schon sortiert sein, aber sicherheitshalber)
        sort_idx = np.argsort(k_clean)
        k_clean = k_clean[sort_idx]
        c_clean = c_clean[sort_idx]
        
        # Duplikate entfernen (insb. nahe kss)
        k_clean, unique_idx = np.unique(k_clean, return_index=True)
        c_clean = c_clean[unique_idx]

        if len(k_clean) < 10:
            return None

        # Spline erstellen
        # k=1 (linear) ist oft robuster bei der Zeit-Eliminierung nahe SS
        try:
            sp = make_interp_spline(k_clean, c_clean, k=1) 
        except Exception:
            return None
            
        return sp

    @reactive.Calc
    def dynamics_data():
        """Berechnet die Zeit-Dynamik, indem k entlang des berechneten Sattelpfads bewegt wird."""
        sp = saddle_path_func()
        ss_data = model_params_ss()
        
        if sp is None or ss_data['error']:
            return None
        
        k0 = input.k0()
        t_end = input.t_end()
        params = ss_data['params']
        kmax = ss_data['kmax']

        # Nur simulieren, wenn k0 im gültigen Bereich des Splines liegt
        # Wir prüfen den Definitionsbereich des Splines (ungefähr 0 bis kmax)
        if k0 <= 0 or k0 > kmax:
            return None
            
        # Differentialgleichung für k(t), wobei c(t) durch den Sattelpfad c(k) bestimmt wird
        def k_dot_func(t, y):
            k_val = y[0]
            if k_val <= 0: return 0
            
            # Hole c vom Sattelpfad
            try:
                c_val = sp(k_val).item()
            except ValueError:
                return 0 # Außerhalb des Bereichs
                
            if c_val < 0: return 0
            
            # dk/dt = k^alpha - c - (n + g + delta)k
            dk = k_val**params['alpha'] - c_val - (params['phi'] + params['xi'] + params['delta']) * k_val
            return dk

        # Löse k(t)
        t_span = (0, t_end)
        sol = solve_ivp(k_dot_func, t_span, [k0], dense_output=True, 
                        method='RK45', t_eval=np.linspace(0, t_end, 300))
        
        if not sol.success:
            return None
            
        k_t = sol.y[0]
        c_t = sp(k_t)
        y_t = k_t**params['alpha']
        
        df = pd.DataFrame({
            't': sol.t,
            'k': k_t,
            'c': c_t,
            'y': y_t
        })
        
        return df

    @output
    @render.ui
    def ss_output():
        """Zeigt Steady State Werte an."""
        ss_data = model_params_ss()
        if ss_data['error']:
            return ui.tags.p(f"FEHLER: {ss_data['error']}")
        
        return ui.div(
            {"class": "ss-output"},
            ui.h4("Stabile Steady-State Werte (pro effektive Einheit):"),
            ui.tags.ul(
                ui.tags.li(ui.HTML(f"<b>Kapital (k̄):</b> {ss_data['kss']:.4f}")),
                ui.tags.li(ui.HTML(f"<b>Konsum (c̄):</b> {ss_data['css']:.4f}")),
                ui.tags.li(ui.HTML(f"<b>Output (ȳ):</b> {ss_data['yss']:.4f}"))
            ),
            ui.tags.p(ui.HTML(f"Maximales Kapital $k_{{max}}$: {ss_data['kmax']:.4f}"))
        )

    @output
    @render.plot
    def phasePlot():
        ss_data = model_params_ss()
        dyn_data = dynamics_data()
        
        if ss_data['error']: return None

        kss, css, kmax, params = ss_data['kss'], ss_data['css'], ss_data['kmax'], ss_data['params']

        fig, ax = plt.subplots(figsize=(9, 9))
        
        # Lokus-Linien
        k_grid = np.linspace(0.01, kmax * 1.05, 300)
        c_k0 = k_grid**params['alpha'] - (params['phi'] + params['xi'] + params['delta']) * k_grid
        
        ax.plot(k_grid, c_k0, label=r'$\dot{k}=0$', color='blue')
        ax.axvline(x=kss, label=r'$\dot{c}=0$', color='red', linestyle='--')
        
        # Sattelpfad
        sp = saddle_path_func()
        if sp is not None:
            c_sp = sp(k_grid)
            # Filter für Plot (nur positive c)
            valid = (c_sp >= 0) & (c_sp < c_k0.max()*1.5)
            ax.plot(k_grid[valid], c_sp[valid], label='Sattelpfad', color='green', linewidth=2)

        # Steady State
        ax.plot(kss, css, 'ko', markersize=8, label='SS')

        # Dynamik Startpunkt
        if dyn_data is not None and len(dyn_data) > 0:
            k0 = dyn_data['k'].iloc[0]
            c0 = dyn_data['c'].iloc[0]
            ax.plot(k0, c0, 'o', color='orange', markersize=8, label='Start')
            # Plotte Verlauf im Phasendiagramm
            ax.plot(dyn_data['k'], dyn_data['c'], color='orange', linestyle=':')

        ax.set_title("Phasendiagramm")
        ax.set_xlabel("Kapital k")
        ax.set_ylabel("Konsum c")
        ax.set_ylim(0, c_k0.max() * 1.2)
        ax.set_xlim(0, kmax * 1.05)
        ax.legend()
        ax.grid(True, alpha=0.3)
        return fig

    @output
    @render.plot
    def kPlot():
        return plot_time_series('k', r"Kapital $k(t)$", r'$\bar{k}$')

    @output
    @render.plot
    def cPlot():
        return plot_time_series('c', r"Konsum $c(t)$", r'$\bar{c}$')

    @output
    @render.plot
    def yPlot():
        return plot_time_series('y', r"Output $y(t)$", r'$\bar{y}$')

    def plot_time_series(var_name, title, ss_label):
        dyn_data = dynamics_data()
        ss_data = model_params_ss()
        
        if dyn_data is None or ss_data['error']: return None
        
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(dyn_data['t'], dyn_data[var_name], linewidth=2.5, color='#004d99')
        
        # SS Linie
        ss_val = ss_data[f'{var_name}ss']
        ax.axhline(y=ss_val, color='black', linestyle='--', label=ss_label)
        
        ax.set_title(title)
        ax.set_xlabel("Zeit t")
        ax.legend()
        ax.grid(True, alpha=0.3)
        return fig

app = App(app_ui, server)