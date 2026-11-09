#!/usr/bin/env Rscript
# Generate every non-handwritten file of the withinr R package.
#
# Run from the repository root:
#
#   Rscript scripts/prepare-r-package.R            # glue + docs (online cargo build)
#   Rscript scripts/prepare-r-package.R --vendor   # ... plus offline vendored deps (CRAN)
#   Rscript scripts/prepare-r-package.R --no-document   # glue only, no compile
#   Rscript scripts/prepare-r-package.R --clean    # delete everything generated
#
# Handwritten sources are `withinr/R/*.R` (except extendr-wrappers.R),
# `withinr/src/rust/{Cargo.toml,Cargo.lock,src/}`, DESCRIPTION, tests, and docs.
# Everything listed in `generated` below is produced here and ignored by Git.
if (!file.exists("Cargo.toml") || !dir.exists("withinr")) {
  stop("Run this script from the within repository root.")
}
args <- commandArgs(trailingOnly = TRUE)
if (any(!args %in% c("--vendor", "--no-document", "--clean"))) {
  stop("Unknown argument; use --vendor, --no-document, or --clean.")
}
root <- normalizePath(".", winslash = "/")
pkg <- file.path(root, "withinr")

# extendr/rextendr scaffolding: FFI entrypoint, build templates, configure
# scripts, the Rust `document` binary, and the R `.Call` wrappers it writes.
scaffolding <- c(
  "src/entrypoint.c", "src/Makevars.in", "src/Makevars.win.in", "src/withinr-win.def",
  "src/rust/document.rs", "tools/msrv.R", "tools/config.R",
  "configure", "configure.win", "cleanup", "cleanup.win"
)
vendored <- c("src/rust/vendor.tar.xz", "src/rust/vendor-config.toml", "inst/COPYRIGHTS")
documentation <- c("NAMESPACE", "man")
# Build leftovers; a stale extracted `src/vendor` would shadow the tarball.
leftovers <- c(
  "src/vendor", "src/rust/vendor", "src/.cargo", "src/Makevars", "src/Makevars.win",
  file.path("src", list.files(file.path(pkg, "src"), "\\.(o|so|dll)$"))
)
generated <- c(scaffolding, "R/extendr-wrappers.R", "src/rust/within-core", leftovers)

if ("--clean" %in% args) {
  unlink(file.path(pkg, c(generated, vendored, documentation, "src/rust/target")),
         recursive = TRUE)
  if (dir.exists(file.path(pkg, "tools")) && !length(dir(file.path(pkg, "tools")))) {
    unlink(file.path(pkg, "tools"), recursive = TRUE)
  }
  message("Removed generated withinr files.")
  quit(save = "no")
}

for (tool in c("rextendr", "roxygen2", "pkgload")) {
  if (!requireNamespace(tool, quietly = TRUE)) stop("Install the R package '", tool, "' first.")
}
if (packageVersion("rextendr") < "0.5.0") stop("rextendr >= 0.5.0 is required.")
if (!nzchar(Sys.which("cargo"))) stop("cargo is not on PATH; install Rust (>= 1.85).")

# Start from a clean slate so no stale artifact can leak into the build (the
# Cargo target directory is kept for incremental compiles). Without --vendor a
# stale vendor archive is removed too, and cargo resolves dependencies online.
unlink(file.path(pkg, c(generated, if (!"--no-document" %in% args) documentation,
                        if (!"--vendor" %in% args) vendored)), recursive = TRUE)

# 1. Snapshot the workspace's Rust core into the package, so the source
#    package builds without this repository or a crates.io release of within.
#    The small nested workspace keeps inherited manifest fields.
snapshot <- file.path(pkg, "src/rust/within-core")
dir.create(file.path(snapshot, "crates"), recursive = TRUE)
manifest <- readLines(file.path(root, "Cargo.toml"))
manifest <- manifest[!grepl('"crates/within-py"', manifest, fixed = TRUE)]
writeLines(manifest, file.path(snapshot, "Cargo.toml"))
stopifnot(file.copy(file.path(root, "LICENSE"), snapshot, overwrite = TRUE))
for (crate in c("within", "schwarz-precond")) {
  from <- file.path(root, "crates", crate)
  to <- file.path(snapshot, "crates", crate)
  dir.create(to)
  # Runtime source only: drop test/benchmark manifest sections.
  manifest <- readLines(file.path(from, "Cargo.toml"))
  manifest <- sub("^\\[package\\]$", '[package]\nworkspace = "../.."', manifest)
  section <- cumsum(grepl("^\\[", manifest))
  headers <- manifest[grepl("^\\[", manifest)]
  drop <- which(grepl("^\\[dev-dependencies|^\\[\\[bench", headers))
  writeLines(manifest[!section %in% drop], file.path(to, "Cargo.toml"))
  stopifnot(file.copy(file.path(from, "src"), to, recursive = TRUE))
  stopifnot(file.copy(file.path(from, "README.md"), to))
}

# 2. Generate the extendr scaffolding in a disposable package, so
#    use_extendr(overwrite = TRUE) cannot touch the handwritten bridge.
scaffold <- tempfile("withinr-scaffold-")
dir.create(scaffold)
stopifnot(file.copy(file.path(pkg, "DESCRIPTION"), scaffold))
rextendr::use_extendr(scaffold, crate_name = "withinr", lib_name = "withinr",
                      overwrite = TRUE, quiet = TRUE)
for (name in c(scaffolding, "R/extendr-wrappers.R")) {
  dest <- file.path(pkg, name)
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file.path(scaffold, name), dest, overwrite = TRUE))
}
Sys.chmod(file.path(pkg, c("configure", "configure.win", "cleanup", "cleanup.win")), "0755")
unlink(scaffold, recursive = TRUE)

# Policy additions to rextendr's build templates: locked dependency resolution
# and at most two compiler jobs (CRAN), even in development builds.
for (name in c("src/Makevars.in", "src/Makevars.win.in")) {
  dest <- file.path(pkg, name)
  lines <- readLines(dest)
  patched <- gsub("cargo (build|run) @CRAN_FLAGS@",
                  "CARGO_BUILD_JOBS=2 cargo \\1 --locked @CRAN_FLAGS@", lines)
  if (sum(patched != lines) != 2L) {
    stop("The rextendr template in ", name, " changed; update the patch in this script.")
  }
  writeLines(patched, dest)
}

# 3. Optionally vendor the locked crates.io dependencies for offline builds.
if ("--vendor" %in% args) {
  rextendr::vendor_crates(pkg, clean = FALSE)
  # Some WASI-only transitive dependencies ship precompiled WebAssembly/object
  # fixtures. CRAN source packages cannot contain executable binaries, and
  # native R builds never use them: drop them and update Cargo's checksums.
  vendor <- file.path(pkg, "src/rust/vendor")
  for (crate in list.dirs(vendor, recursive = FALSE)) {
    files <- list.files(crate, recursive = TRUE, all.files = TRUE)
    binaries <- files[grepl("\\.(wasm|o|a|so|dll|dylib|exe)$", files)]
    if (length(binaries)) {
      checksum_path <- file.path(crate, ".cargo-checksum.json")
      checksum <- jsonlite::read_json(checksum_path)
      for (name in binaries) {
        unlink(file.path(crate, name))
        checksum$files[[name]] <- NULL
      }
      jsonlite::write_json(checksum, checksum_path, auto_unbox = TRUE)
    }
  }
  processx::run("tar", c("-cJ", "--no-xattrs", "-f", "vendor.tar.xz", "vendor"),
                wd = file.path(pkg, "src/rust"))
  unlink(vendor, recursive = TRUE)
  # Attribution for every bundled crate; original notices stay in the archive.
  metadata <- processx::run("cargo", c("metadata", "--locked", "--offline",
    "--format-version=1", "--manifest-path", file.path(pkg, "src/rust/Cargo.toml")))
  crates <- jsonlite::fromJSON(metadata$stdout, simplifyVector = FALSE)$packages
  crates <- crates[order(vapply(crates, function(x) x$name, character(1)))]
  lines <- c("Bundled Rust source attribution (generated; do not edit).",
    "Original license and copyright notices are retained in src/rust/vendor.tar.xz",
    "and src/rust/within-core. The source manifest of each crate is authoritative.",
    "Precompiled WASI objects/fixtures are omitted. This package targets native R.", "")
  for (crate in crates) {
    lines <- c(lines, paste(crate$name, crate$version),
      paste("Authors:", paste(unlist(crate$authors), collapse = "; ")),
      paste("License:", crate$license), paste("Repository:", crate$repository), "")
  }
  dir.create(file.path(pkg, "inst"), showWarnings = FALSE)
  writeLines(lines, file.path(pkg, "inst/COPYRIGHTS"), useBytes = TRUE)
}

# 4. Compile, write the real `.Call` wrappers, and document.
if (!"--no-document" %in% args) {
  # Bootstrap NAMESPACE (useDynLib) from source without loading a DLL.
  roxygen2::roxygenise(pkg, load_code = roxygen2::load_source)
  # pkgload compiles the crate and runs its `document` binary, which writes
  # R/extendr-wrappers.R from extendr's metadata; roxygen then writes the
  # final NAMESPACE and man/ pages.
  roxygen2::roxygenise(pkg)
}
message("withinr sources prepared in ", pkg)
