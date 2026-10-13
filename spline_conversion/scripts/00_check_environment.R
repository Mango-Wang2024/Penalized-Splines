#!/usr/bin/env Rscript
# Read-only environment check. This script installs nothing.
cat("R:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n\n")
missing <- character()
for (pkg in c("mgcv", "nlme")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, ":", as.character(packageVersion(pkg)), "\n")
  } else {
    missing <- c(missing,pkg)
    cat(pkg, ": MISSING\n")
  }
}
if (length(missing)) {
  cat("\nIn R, install the missing packages with:\n")
  cat(sprintf('install.packages(c(%s), repos="https://cloud.r-project.org")\n',
              paste(sprintf('"%s"',missing),collapse=", ")))
  quit(status=1)
}
cat("\nReady for the reference export. Versions will be recorded, not silently pinned.\n")
