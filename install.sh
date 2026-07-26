#!/usr/bin/env bash
# FlexAppeal developer installer -- installs pixi if needed, then solves the
# development environment declared in pixi.toml (OpenMM, MDTraj, Flask, pytest).
#
# This is for working ON FlexAppeal. If you just want to RUN a simulation, you
# do not need any of this: download a bundle from the Prepare tab and run it.
#
# Usage: ./install.sh
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

BLUE='\033[38;2;30;115;190m'
GREEN='\033[38;2;0;208;132m'
AMBER='\033[38;2;252;185;0m'
RESET='\033[0m'
BOLD='\033[1m'

info()  { printf "${BLUE}ℹ${RESET} %s\n" "$1"; }
ok()    { printf "${GREEN}✓${RESET} %s\n" "$1"; }
step()  { printf "${BLUE}→${RESET} %s\n" "$1"; }
warn()  { printf "${AMBER}⚠${RESET} %s\n" "$1"; }

printf '%b' "\n${BOLD}${BLUE}"
cat <<'BANNER'
 ___ _   _____   __ __  ___ ___ ___  __  _
| __| | | __\ \_/ //  \| _,\ _,\ __|/  \| |
| _|| |_| _| > , <| /\ | v_/ v_/ _|| /\ | |_
|_| |___|___/_/ \_\_||_|_| |_| |___|_||_|___|
BANNER
printf '%s' "${RESET}"
printf "\n"
printf "  ${BLUE}molecular dynamics, prepared and analysed${RESET}\n\n"

case "$(uname -s)" in
    Darwin) ;;
    Linux)  ;;
    *)
        warn "unrecognised platform '$(uname -s)' -- FlexAppeal is developed on macOS (Apple Silicon) and deployed on Linux (x86-64)."
        ;;
esac

if ! command -v pixi >/dev/null 2>&1; then
    step "pixi not found -- installing it (https://pixi.sh)"
    curl -fsSL https://pixi.sh/install.sh | sh
    # The installer adds pixi to shell rc files for future sessions; make it
    # available in this one without requiring a shell restart.
    export PATH="$HOME/.pixi/bin:$PATH"
    if ! command -v pixi >/dev/null 2>&1; then
        warn "pixi installed but not yet on PATH in this shell."
        info "Run: source ~/.zshrc  (or restart your terminal), then re-run ./install.sh"
        exit 1
    fi
else
    ok "pixi already installed ($(pixi --version))"
fi

step "solving and installing the FlexAppeal development environment"
info "this pulls OpenMM, AmberTools and the OpenFF toolkit -- several GB, one-time"
pixi install

echo
ok "install complete."
info "Run the tests:        pixi run test"
info "Serve the app:        pixi run serve   (then open http://127.0.0.1:8004)"
info "Open a shell in it:   pixi shell"
