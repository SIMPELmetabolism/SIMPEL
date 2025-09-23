# .github/scripts/fetch_traffic.R
library(httr)
library(jsonlite)

# Repo + token (injected by GitHub Actions)
owner <- Sys.getenv("GITHUB_OWNER")
repo  <- Sys.getenv("GITHUB_REPO")
token <- Sys.getenv("GH_PAT")

headers <- add_headers(Authorization = paste("token", token))

# --- Helper function to fetch traffic ---
fetch_traffic <- function(url) {
  res <- content(GET(url, headers), as = "parsed", type = "application/json")

  # Daily data: if empty, create empty data frame with proper columns
  daily <- if (!is.null(res[[basename(url)]])) {
    df <- as.data.frame(res[[basename(url)]])
    df$fetched_at <- Sys.Date()
    df
  } else {
    data.frame(timestamp=character(),
               count=integer(),
               uniques=integer(),
               fetched_at=as.Date(character()))
  }
  
  # Totals: always create 1-row data frame
  totals <- data.frame(
    metric     = basename(url),
    total      = ifelse(is.null(res$count), 0, res$count),
    uniques    = ifelse(is.null(res$uniques), 0, res$uniques),
    fetched_at = Sys.Date()
  )
  
  list(totals = totals, daily = daily)
}

# --- API URLs ---
url_clones <- paste0("https://api.github.com/repos/", owner, "/", repo, "/traffic/clones")
url_views  <- paste0("https://api.github.com/repos/", owner, "/", repo, "/traffic/views")

# --- Fetch both metrics ---
clones <- fetch_traffic(url_clones)
views  <- fetch_traffic(url_views)

# --- Save directory ---
dir.create("traffic_logs", showWarnings = FALSE)

# --- Append helper ---
append_or_write <- function(df, file) {
  if (file.exists(file)) {
    old <- read.csv(file)
    new <- rbind(old, df)
    # remove exact duplicates (in case of overlap)
    new <- unique(new)
    write.csv(new, file, row.names = FALSE)
  } else {
    write.csv(df, file, row.names = FALSE)
  }
}

# --- Save daily stats ---
append_or_write(clones$daily, "traffic_logs/clones_daily.csv")
append_or_write(views$daily,  "traffic_logs/views_daily.csv")

# --- Save totals ---
append_or_write(rbind(clones$totals, views$totals),
                "traffic_logs/traffic_totals.csv")
