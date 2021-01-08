# Load libraries
library(shiny)
library(shinydashboard)
library(rhandsontable)
library(patchwork)
library(tidyverse)

# Body of app (Layout + UI)
body <-
  dashboardBody(
    tabItems(
      # Magma Density Tab Layout
      tabItem('density',
              # Row
              fluidRow(
                # Data Box
                box(width = 8,
                    title = 'User Input',
                    status = 'primary',
                    solidHeader = TRUE,
                    # Data table
                    rHandsontableOutput(outputId = 'table', height = '315px'),
                    # Toggle box
                    checkboxInput('toggle', 'Toggle All', value = FALSE)),
                # UI Box
                box(width = 4,
                    title = 'Parameters',
                    status = 'primary',
                    solidHeader = TRUE,
                    numericInput('Tlow1', 'T Range (Celcius)', value = 800, min = 0, max = 2000, step = 25),
                    numericInput('Thigh1', NULL, value = 1500, min = 0, max = 2000, step = 25),
                    numericInput('Plow1', 'P Range (bar)', value = 1000, min = 0, max = 20000, step = 100),
                    numericInput('Phigh1', NULL, value = 10000, min = 0, max = 20000, step = 100),
                    # Select plot
                    selectInput('type1', label = 'Plot Type', choices = c('Summary', 'Summary Density', 'Density T', 'Density P'), selected = 'Summary'),
                    # Calculate density button
                    actionButton('calc1', 'Calculate Density')),
                # Plot Box
                box(width = 7,
                  title = 'Visualization',
                  status = 'primary',
                  solidHeader = TRUE,
                  plotOutput(outputId = 'p1')),
                # Results Box
                box(width = 5,
                    title = 'Density Results',
                    status = 'primary',
                    solidHeader = TRUE,
                    # Results Table
                    rHandsontableOutput(outputId = 'table2', height = '400px')))),
      # Magma viscosity tab layout
      tabItem('viscosity',
              # Row
              fluidRow(
                # Results Box
                box(width = 6,
                  title = 'Density Results',
                  status = 'primary',
                  solidHeader = TRUE,
                  # Results Table
                  rHandsontableOutput(outputId = 'table3', height = '400px')),
                # Plot box
                box(width = 6,
                    title = 'Visualization',
                    status = 'primary',
                    solidHeader = TRUE,
                    plotOutput(outputId = 'p2')),
                # UI box
                box(width = 2,
                    title = 'Parameters',
                    status = 'primary',
                    solidHeader = TRUE,
                    # Crystal fraction
                    numericInput('phi', 'Crystal Fraction', value = 0, min = 0, max = 1, step = 0.05),
                    # Plot type
                    selectInput('type2', label = 'Plot Type', choices = c('Summary', 'Viscosity T', 'Viscosity H2O'), selected = 'Viscosity T'),
                    # Calculate viscosity button
                    actionButton('calc2', 'Calculate Viscosity')),
                # Results Box
                box(width = 10,
                  title = 'Viscosity Results',
                  status = 'primary',
                  solidHeader = TRUE,
                  # Results Table
                  rHandsontableOutput(outputId = 'table4', height = '360px')))),
      # Stokes velocity tab layout
      tabItem('stokes',
              # Row
              fluidRow(
                # Plot box
                box(width = 12,
                    title = 'Visualization',
                    status = 'primary',
                    solidHeader = TRUE,
                    plotOutput(outputId = 'p3')),
                # Results Box
                box(width = 10,
                  title = 'Viscosity Results',
                  status = 'primary',
                  solidHeader = TRUE,
                  # Results Table
                  rHandsontableOutput(outputId = 'table5', height = '360px')),
                # UI box
                box(width = 2,
                    title = 'Parameters',
                    status = 'primary',
                    solidHeader = TRUE,
                    # Crystal Parameters
                    'Crystal Parameters',
                    br(),
                    numericInput('rad', 'Radius (m)', value = 0.001, min = 0, max = 0.1, step = 0.0001),
                    numericInput('rho', 'Density (kg/m^3)', value = 3000, min = 1000, max = 5000, step = 25),
                    # Calculate viscosity button
                    actionButton('calc3', 'Calculate Velocity')))),
      # Options tab
      tabItem('options',
              # Row
              fluidRow(
                # First box (download configurations)
                box(width = 2,
                    title = 'Configure Download',
                    status = 'primary',
                    solidHeader = TRUE,
                    selectInput('p.select', 'Select Plot', choices = c('Magma Density', 'Magma Viscosity', 'Stokes Velocity')),
                    textInput('filename', 'Filename (no extension)'),
                    numericInput('pdf_width', 'Plot width (inches)', min = 5.5, max = 17, step = 0.25, value = 5.5),
                    numericInput('pdf_height', 'Plot height (inches)', min = 4.25, max = 11, step = 0.25, value = 4.25),
                    downloadButton('dlplot', 'Get Plot')),
                # Plot box
                box(width = 10,
                    title = 'Stokes Velocity',
                    status = 'primary',
                    solidHeader = TRUE,
                    plotOutput(outputId = 'p6', height = '340px')),
                # Plot box
                box(
                    title = 'Magma Viscosity',
                    status = 'primary',
                    solidHeader = TRUE,
                    plotOutput(outputId = 'p5')),
                # Plot box
                box(
                    title = 'Magma Density',
                    status = 'primary',
                    solidHeader = TRUE,
                    plotOutput(outputId = 'p4')))),
    # Final tab layout
    tabItem('ack',
            # Row
            fluidRow(
              # Text output for acknowledgements
              box(width = 12,
                  status = 'primary',
                  solidHeader = TRUE,
                  verbatimTextOutput('acknowledgements'))))))
# Magma density function
magma_density <- function(SiO2, TiO2, Al2O3, Fe2O3, FeO, MgO, CaO, Na2O, K2O, H2O, P, T, IDs = NULL){
  # Save inputs & tidy
  d <- tibble(
    SiO2 = SiO2,
    TiO2 = TiO2,
    Al2O3 = Al2O3,
    Fe2O3 = Fe2O3,
    FeO = FeO,
    MgO = MgO,
    CaO = CaO,
    Na2O = Na2O,
    K2O = K2O,
    H2O = H2O
  ) %>%
    replace(is.na(.), 0) %>%
    mutate(across(where(is.logical), as.numeric))
  
  # Molecular Weights
  mw.SiO2 <- 60.0855
  mw.TiO2 <- 79.88
  mw.Al2O3 <- 101.96
  mw.Fe2O3 <- 159.69
  mw.FeO <- 71.85
  mw.MgO <- 40.3
  mw.CaO <- 56.08
  mw.Na2O <- 61.98
  mw.K2O <- 94.2
  mw.H2O <- 18.02
  
  # Partial Molar Volumes
  # Volumes for SiO2, Al2O3, MgO, CaO, Na2O, K2O at Tref=1773 K (Lange, 1997; CMP)
  # Volume for H2O at Tref=1273 K (Ochs and Lange, 1999)
  # Volume for FeO at Tref=1723 K (Guo et al., 2014)
  # Volume for Fe2O3 at Tref=1723 K (Liu and Lange, 2006)
  # Volume for TiO2 at Tref=1773 K (Lange and Carmichael, 1987)
  mv.SiO2 <- 26.86
  mv.TiO2 <- 28.32
  mv.Al2O3 <- 37.42
  mv.Fe2O3 <- 41.50
  mv.FeO <- 12.68
  mv.MgO <- 12.02
  mv.CaO <- 16.90
  mv.Na2O <- 29.65
  mv.K2O <- 47.28
  mv.H2O <- 22.9
  
  # Partial Molar Volume uncertainties
  # value = 0 if not reported
  sig.mv.SiO2 <- 0.03
  sig.mv.TiO2 <- 0
  sig.mv.Al2O3 <- 0.09
  sig.mv.Fe2O3 <- 0
  sig.mv.FeO <- 0
  sig.mv.MgO <- 0.07
  sig.mv.CaO <- 0.06
  sig.mv.Na2O <- 0.07
  sig.mv.K2O <- 0.10
  sig.mv.H2O <- 0.60
  
  # dV/dT values
  # MgO, CaO, Na2O, K2O Table 4 (Lange, 1997)
  # SiO2, TiO2, Al2O3 Table 9 (Lange and Carmichael, 1987)
  # H2O from Ochs & Lange (1999)
  # Fe2O3 from Liu & Lange (2006)
  # FeO from Guo et al (2014)
  dvdt.SiO2 <- 0.0
  dvdt.TiO2 <- 0.00724
  dvdt.Al2O3 <- 0.00262
  dvdt.Fe2O3 <- 0.0
  dvdt.FeO <- 0.00369
  dvdt.MgO <- 0.00327
  dvdt.CaO <- 0.00374
  dvdt.Na2O <- 0.00768
  dvdt.K2O <- 0.01208
  dvdt.H2O <- 0.0095
  
  # dV/dT uncertainties
  # value <- 0 if not reported
  sig.dvdt.SiO2 <- 0
  sig.dvdt.TiO2 <- 0
  sig.dvdt.Al2O3 <- 0
  sig.dvdt.Fe2O3 <- 0
  sig.dvdt.FeO <- 0
  sig.dvdt.MgO <- 0
  sig.dvdt.CaO <- 0
  sig.dvdt.Na2O <- 0
  sig.dvdt.K2O <- 0
  sig.dvdt.H2O <- 0.00080
  
  # dV/dP values
  # Anhydrous component data from Kess and Carmichael (1991)
  # H2O data from Ochs & Lange (1999)
  dvdp.SiO2 <- -0.000189
  dvdp.TiO2 <- -0.000231
  dvdp.Al2O3 <- -0.000226
  dvdp.Fe2O3 <- -0.000253
  dvdp.FeO <- -0.000045
  dvdp.MgO <- 0.000027
  dvdp.CaO <- 0.000034
  dvdp.Na2O <- -0.00024
  dvdp.K2O <- -0.000675
  dvdp.H2O <- -0.00032
  
  # dV/dP uncertainties
  sig.dvdp.SiO2 <- 0.000002
  sig.dvdp.TiO2 <- 0.000006
  sig.dvdp.Al2O3 <- 0.000009
  sig.dvdp.Fe2O3 <- 0.000009
  sig.dvdp.FeO <- 0.000003
  sig.dvdp.MgO <- 0.000007
  sig.dvdp.CaO <- 0.000005
  sig.dvdp.Na2O <- 0.000005
  sig.dvdp.K2O <- 0.000014
  sig.dvdp.H2O <- 0.000060
  
  # Tref values
  T.ref.SiO2 <- 1773
  T.ref.TiO2 <- 1773
  T.ref.Al2O3 <- 1773
  T.ref.Fe2O3 <- 1723
  T.ref.FeO <- 1723
  T.ref.MgO <- 1773
  T.ref.CaO <- 1773
  T.ref.Na2O <- 1773
  T.ref.K2O <- 1773
  T.ref.H2O <- 1273
  
  # Sum weights
  d <- d %>% mutate(total = rowSums(., na.rm = T))
  
  # Normalize
  d.norm <- d %>% 
    transmute(norm.SiO2 = SiO2/total,
              norm.TiO2 = TiO2/total,
              norm.Al2O3 = Al2O3/total,
              norm.Fe2O3 = Fe2O3/total,
              norm.FeO = FeO/total,
              norm.MgO = MgO/total,
              norm.CaO = CaO/total,
              norm.Na2O = Na2O/total,
              norm.K2O = K2O/total,
              norm.H2O = H2O/total) %>%
    mutate(total = rowSums(., na.rm = T))
  
  # Calculate molar proportions
  d.mol.p <- d.norm %>%
    transmute(mol.p.SiO2 = norm.SiO2/mw.SiO2,
              mol.p.TiO2 = norm.TiO2/mw.TiO2,
              mol.p.Al2O3 = norm.Al2O3/mw.Al2O3,
              mol.p.Fe2O3 = norm.Fe2O3/mw.Fe2O3,
              mol.p.FeO = norm.FeO/mw.FeO,
              mol.p.MgO = norm.MgO/mw.MgO,
              mol.p.CaO = norm.CaO/mw.CaO,
              mol.p.Na2O = norm.Na2O/mw.Na2O,
              mol.p.K2O = norm.K2O/mw.K2O,
              mol.p.H2O = norm.H2O/mw.H2O) %>%
    mutate(total = rowSums(., na.rm = T))
  
  # Normalize
  d.mol.frac <- d.mol.p %>%
    transmute(mol.frac.SiO2 = mol.p.SiO2/total,
              mol.frac.TiO2 = mol.p.TiO2/total,
              mol.frac.Al2O3 = mol.p.Al2O3/total,
              mol.frac.Fe2O3 = mol.p.Fe2O3/total,
              mol.frac.FeO = mol.p.FeO/total,
              mol.frac.MgO = mol.p.MgO/total,
              mol.frac.CaO = mol.p.CaO/total,
              mol.frac.Na2O = mol.p.Na2O/total,
              mol.frac.K2O = mol.p.K2O/total,
              mol.frac.H2O = mol.p.H2O/total) %>%
    mutate(total = rowSums(., na.rm = T))
  
  # Convert T to K
  TK <- T + 273.15
  
  # Calculate densities m/V = mol.frac * mw / mv + (dVdT * (TK - Tref)) + (dVdP * (P - 1))
  
  d.dens <- d.mol.frac %>%
    transmute(rho.SiO2 = (mol.frac.SiO2 * mw.SiO2) / (mv.SiO2 + (dvdt.SiO2 * (TK - T.ref.SiO2)) + (dvdp.SiO2 * (P - 1))),
              rho.TiO2 = (mol.frac.TiO2 * mw.TiO2) / (mv.TiO2 + (dvdt.TiO2 * (TK - T.ref.TiO2)) + (dvdp.TiO2 * (P - 1))),
              rho.Al2O3 = (mol.frac.Al2O3 * mw.Al2O3) / (mv.Al2O3 + (dvdt.Al2O3 * (TK - T.ref.Al2O3)) + (dvdp.Al2O3 * (P - 1))),
              rho.Fe2O3 = (mol.frac.Fe2O3 * mw.Fe2O3) / (mv.Fe2O3 + (dvdt.Fe2O3 * (TK - T.ref.Fe2O3)) + (dvdp.Fe2O3 * (P - 1))),
              rho.FeO = (mol.frac.FeO * mw.FeO) / (mv.FeO + (dvdt.FeO * (TK - T.ref.FeO)) + (dvdp.FeO * (P - 1))),
              rho.MgO = (mol.frac.MgO * mw.MgO) / (mv.MgO + (dvdt.MgO * (TK - T.ref.MgO)) + (dvdp.MgO * (P - 1))),
              rho.CaO = (mol.frac.CaO * mw.CaO) / (mv.CaO + (dvdt.CaO * (TK - T.ref.CaO)) + (dvdp.CaO * (P - 1))),
              rho.Na2O = (mol.frac.Na2O * mw.Na2O) / (mv.Na2O + (dvdt.Na2O * (TK - T.ref.Na2O)) + (dvdp.Na2O * (P - 1))),
              rho.K2O = (mol.frac.K2O * mw.K2O) / (mv.K2O + (dvdt.K2O * (TK - T.ref.K2O)) + (dvdp.K2O * (P - 1))),
              rho.H2O = (mol.frac.H2O * mw.H2O) / (mv.H2O + (dvdt.H2O * (TK - T.ref.H2O)) + (dvdp.H2O * (P - 1))))
  # Calculate liquid volume = (mv + (dVdT * (TK - Tref)) + (dVdP * (P - 1))) * mol.frac
  d.vol.liq <- d.mol.frac %>%
    transmute(vol.liq.SiO2 = (mv.SiO2 + (dvdt.SiO2 * (TK - T.ref.SiO2)) + (dvdp.SiO2 * (P - 1))) * mol.frac.SiO2,
              vol.liq.TiO2 = (mv.TiO2 + (dvdt.TiO2 * (TK - T.ref.TiO2)) + (dvdp.TiO2 * (P - 1))) * mol.frac.TiO2,
              vol.liq.Al2O3 = (mv.Al2O3 + (dvdt.Al2O3 * (TK - T.ref.Al2O3)) + (dvdp.Al2O3 * (P - 1))) * mol.frac.Al2O3,
              vol.liq.Fe2O3 = (mv.Fe2O3 + (dvdt.Fe2O3 * (TK - T.ref.Fe2O3)) + (dvdp.Fe2O3 * (P - 1))) * mol.frac.Fe2O3,
              vol.liq.FeO = (mv.FeO + (dvdt.FeO * (TK - T.ref.FeO)) + (dvdp.FeO * (P - 1))) * mol.frac.FeO,
              vol.liq.MgO = (mv.MgO + (dvdt.MgO * (TK - T.ref.MgO)) + (dvdp.MgO * (P - 1))) * mol.frac.MgO,
              vol.liq.CaO = (mv.CaO + (dvdt.CaO * (TK - T.ref.CaO)) + (dvdp.CaO * (P - 1))) * mol.frac.CaO,
              vol.liq.Na2O = (mv.Na2O + (dvdt.Na2O * (TK - T.ref.Na2O)) + (dvdp.Na2O * (P - 1))) * mol.frac.Na2O,
              vol.liq.K2O = (mv.K2O + (dvdt.K2O * (TK - T.ref.K2O)) + (dvdp.K2O * (P - 1))) * mol.frac.K2O,
              vol.liq.H2O = (mv.H2O + (dvdt.H2O * (TK - T.ref.H2O)) + (dvdp.H2O * (P - 1))) * mol.frac.H2O) %>%
    mutate(total = rowSums(., na.rm = T))
  
  # Calculate weights
  d.w <- d.mol.frac %>%
    transmute(w.SiO2 = mw.SiO2 * mol.frac.SiO2,
              w.TiO2 = mw.TiO2 * mol.frac.TiO2,
              w.Al2O3 = mw.Al2O3 * mol.frac.Al2O3,
              w.Fe2O3 = mw.Fe2O3 * mol.frac.Fe2O3,
              w.FeO = mw.FeO * mol.frac.FeO,
              w.MgO = mw.MgO * mol.frac.MgO,
              w.CaO = mw.CaO * mol.frac.CaO,
              w.Na2O = mw.Na2O * mol.frac.Na2O,
              w.K2O = mw.K2O * mol.frac.K2O,
              w.H2O = mw.H2O * mol.frac.H2O) %>%
    mutate(total = rowSums(., na.rm = T))
  
  # Calculate density of melt [gcm^-3]
  d.rho.liq <- d.w$total / d.vol.liq$total
  
  # Calculate uncertainties
  d.sig.mv <- tibble(SiO2 = sig.mv.SiO2 / mv.SiO2,
                     TiO2 = sig.mv.TiO2 / mv.TiO2,
                     Al2O3 = sig.mv.Al2O3 / mv.Al2O3,
                     Fe2O3 = sig.mv.Fe2O3 / mv.Fe2O3,
                     FeO = sig.mv.FeO / mv.FeO,
                     MgO = sig.mv.MgO / mv.MgO,
                     CaO = sig.mv.CaO / mv.CaO,
                     Na2O = sig.mv.Na2O / mv.Na2O,
                     K2O = sig.mv.K2O / mv.K2O,
                     H2O = sig.mv.H2O / mv.H2O)
  d.sig.dvdt <- tibble(SiO2 = sig.dvdt.SiO2 / 1,
                       TiO2 = sig.dvdt.TiO2 / dvdt.TiO2,
                       Al2O3 = sig.dvdt.Al2O3 / dvdt.Al2O3,
                       Fe2O3 = 0,
                       FeO = sig.dvdt.FeO / dvdt.FeO,
                       MgO = sig.dvdt.MgO / dvdt.MgO,
                       CaO = sig.dvdt.CaO / dvdt.CaO,
                       Na2O = sig.dvdt.Na2O / dvdt.Na2O,
                       K2O = sig.dvdt.K2O / dvdt.K2O,
                       H2O = sig.dvdt.H2O / dvdt.H2O)
  d.sig.dvdp <- tibble(SiO2 = sig.dvdp.SiO2 / dvdp.SiO2,
                       TiO2 = sig.dvdp.TiO2 / dvdp.TiO2,
                       Al2O3 = sig.dvdp.Al2O3 / dvdp.Al2O3,
                       Fe2O3 = sig.dvdp.Fe2O3 / dvdp.Fe2O3,
                       FeO = sig.dvdp.FeO / dvdp.FeO,
                       MgO = sig.dvdp.MgO / dvdp.MgO,
                       CaO = sig.dvdp.CaO / dvdp.CaO,
                       Na2O = sig.dvdp.Na2O / dvdp.Na2O,
                       K2O = sig.dvdp.K2O / dvdp.K2O,
                       H2O = sig.dvdp.H2O / dvdp.H2O)
  # Calculate percent error for vol liquid = sqrt(sig.mv^2 + sig.dvdt^2 + sig.dvdp^2)
  perc.error.vol.liq <- purrr::pmap_df(list(d.sig.mv, d.sig.dvdt, d.sig.dvdp), ~sqrt(..1^2 + ..2^2 + ..3^2))
  
  d.sig.vol.liq <- d.vol.liq %>%
    transmute(sig.SiO2 = vol.liq.SiO2 * perc.error.vol.liq$SiO2,
              sig.TiO2 = vol.liq.TiO2 * perc.error.vol.liq$TiO2,
              sig.Al2O3 = vol.liq.Al2O3 * perc.error.vol.liq$Al2O3,
              sig.Fe2O3 = vol.liq.Fe2O3 * perc.error.vol.liq$Fe2O3,
              sig.FeO = vol.liq.FeO * perc.error.vol.liq$FeO,
              sig.MgO = vol.liq.MgO * perc.error.vol.liq$MgO,
              sig.CaO = vol.liq.CaO * perc.error.vol.liq$CaO,
              sig.Na2O = vol.liq.Na2O * perc.error.vol.liq$Na2O,
              sig.K2O = vol.liq.K2O * perc.error.vol.liq$K2O,
              sig.H2O = vol.liq.H2O * perc.error.vol.liq$H2O) %>%
    mutate(total = rowSums(., na.rm = T))
  
  # Calculate uncertainty for density
  d.sig.rho <- d.sig.vol.liq$total / d.vol.liq$total
  
  # Summary data
  return(tibble(ID = IDs,
                SiO2 = d.norm$norm.SiO2*100,
                H2O = d.norm$norm.H2O*100,
                T = T,
                P = P/1000,
                rho.melt = d.rho.liq,
                sig.rho = d.sig.rho))
}
