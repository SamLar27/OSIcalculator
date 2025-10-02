OSI_lookup <- function(FEV1_PCT, Tif_PCT, table = NULL, load_path = NULL) {
  # --- Load table if needed ---
  if (is.null(table)) {
    if (is.null(load_path)) {
      user_name <- Sys.getenv("USER")
      load_path <- file.path(
        "/Users", user_name,
        "Library/CloudStorage/OneDrive-USherbrooke/Recherche/Projets/Pneumologie/ORACLE/ORACLE_OSI/Analysis"
      )
    }
    rdata_file <- file.path(load_path, "Table_FEV1_Tif_OSI.RData")
    rds_file   <- file.path(load_path, "Table_FEV1_Tif_OSI.rds")
    if (file.exists(rdata_file)) {
      load(rdata_file)  # should create Table_FEV1_Tif_OSI
      if (exists("Table_FEV1_Tif_OSI", inherits = FALSE)) table <- Table_FEV1_Tif_OSI
    } else if (file.exists(rds_file)) {
      table <- readRDS(rds_file)
    } else {
      stop("Lookup table not provided and no local file found.")
    }
  }
  
  # --- Normalize & checks ---
  table <- as.data.frame(table)
  if (!all(c("FEV1_PCT","Tif_PCT","OSI") %in% names(table))) {
    stop("`table` must contain columns: FEV1_PCT, Tif_PCT, OSI")
  }
  if (!is.numeric(FEV1_PCT) || !is.numeric(Tif_PCT)) {
    stop("FEV1_PCT and Tif_PCT must be numeric")
  }
  if (length(FEV1_PCT) != length(Tif_PCT)) {
    stop("FEV1_PCT and Tif_PCT must have the same length.")
  }
  
  # Ensure numeric keys in the table
  table$FEV1_PCT <- as.numeric(table$FEV1_PCT)
  table$Tif_PCT  <- as.numeric(table$Tif_PCT)
  
  # --- Snap to nearest grid present in the table ---
  grid_fev1 <- sort(unique(table$FEV1_PCT))
  grid_tif  <- sort(unique(table$Tif_PCT))
  
  snap_to_grid <- function(x, grid) {
    idx  <- findInterval(x, grid, all.inside = TRUE)
    idx2 <- pmin(idx + 1L, length(grid))
    left  <- grid[idx]
    right <- grid[idx2]
    pick_right <- abs(x - right) < abs(x - left)
    ifelse(pick_right, right, left)
  }
  
  n <- length(FEV1_PCT)
  out <- rep(NA_real_, n)
  na_mask <- is.na(FEV1_PCT) | is.na(Tif_PCT)
  if (all(na_mask)) return(out)
  
  fev1_snap <- snap_to_grid(FEV1_PCT[!na_mask], grid_fev1)
  tif_snap  <- snap_to_grid(Tif_PCT[!na_mask],  grid_tif)
  
  # --- Key-based lookup (no data.table [,] semantics) ---
  make_key <- function(a,b) paste0(a, "||", b)
  table$.__key__ <- make_key(table$FEV1_PCT, table$Tif_PCT)
  qkey <- make_key(fev1_snap, tif_snap)
  
  idx <- match(qkey, table$.__key__)
  out[!na_mask] <- table$OSI[idx]
  
  if (any(is.na(out[!na_mask]))) {
    warning(sprintf("OSI not found for %d snapped pair(s) (nearest-grid lookup).",
                    sum(is.na(out[!na_mask]))))
  }
  out
}
