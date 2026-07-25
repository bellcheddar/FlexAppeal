#!/usr/bin/env bash
# Push FlexAppeal from your Mac to the droplet and restart the web service.
# Run from the repo root:  bash deploy/deploy.sh
#
# Reads DROPLET_SSH / DROPLET_PATH from .env (see deploy/.env.example).
# Idempotent, and excludes the venv, scratch and secrets so the server's own
# state is never clobbered.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if [[ -f .env ]]; then set -a; source .env; set +a; fi
DROPLET_SSH="${DROPLET_SSH:-}"
DROPLET_PATH="${DROPLET_PATH:-/opt/flexappeal}"
SSH_KEY="${SSH_KEY:-}"

if [[ -z "$DROPLET_SSH" ]]; then
    echo "DROPLET_SSH is not set. Copy deploy/.env.example to .env and fill it in."
    exit 1
fi

SSH_OPTS=()
[[ -n "$SSH_KEY" ]] && SSH_OPTS=(-e "ssh -i ${SSH_KEY/#\~/$HOME}")

echo "==> Syncing to ${DROPLET_SSH}:${DROPLET_PATH}"
# The vendored JavaScript is ~10 MB and rsync will only send it once, but it is
# genuinely needed on the server -- nginx serves it from disk. Not excluded.
#
# .pixi/ is the local development environment (several GB of OpenMM); the droplet
# neither has nor needs it. web_scratch/ holds live user sessions.
#
# ${arr[@]+"${arr[@]}"} expands to nothing when empty without tripping `set -u`,
# which macOS's bash 3.2 requires.
rsync -az --delete ${SSH_OPTS[@]+"${SSH_OPTS[@]}"} \
    --exclude '.venv/' --exclude '.venv-dev/' --exclude '.pixi/' \
    --exclude '__pycache__/' --exclude '*.pyc' --exclude '.git/' \
    --exclude '.env' --exclude 'web_scratch/' --exclude '.pytest_cache/' \
    ./ "${DROPLET_SSH}:${DROPLET_PATH}/"

echo "==> Installing dependencies and restarting on the droplet"
SSH_CMD=(ssh)
[[ -n "$SSH_KEY" ]] && SSH_CMD=(ssh -i "${SSH_KEY/#\~/$HOME}")
"${SSH_CMD[@]}" "$DROPLET_SSH" bash -s <<REMOTE
set -euo pipefail
cd "${DROPLET_PATH}"

if [[ ! -x .venv/bin/python ]]; then
    echo "No virtual environment yet -- run deploy/provision.sh as root first."
    exit 1
fi

sudo -u flexappeal env PIP_NO_CACHE_DIR=1 ./.venv/bin/pip install --quiet -r requirements.txt

# rsync runs as root, so every file it wrote is root-owned and mode 0600 --
# rsync preserves the Mac's permissions. The service user then hits
# PermissionError at RUNTIME rather than at deploy time, which is a confusing
# way to find out. Chown on every deploy, not only at provisioning.
#
# The venv and live scratch sessions are pruned: nothing in this sync touched
# them, and a blanket chown over a session an analysis has open is asking for
# trouble.
sudo find "${DROPLET_PATH}" \
    -path "${DROPLET_PATH}/.venv" -prune -o \
    -path "${DROPLET_PATH}/web_scratch" -prune -o \
    -exec chown flexappeal:flexappeal {} +

sudo systemctl restart flexappeal-web.service
sleep 1
sudo systemctl --no-pager --lines=3 status flexappeal-web.service || true
REMOTE

echo "==> Verifying"
SERVER_NAME="${SERVER_NAME:-flexappeal.mdeller.com}"
CODE=$(curl -s -o /dev/null -w '%{http_code}' "https://${SERVER_NAME}/healthz" || echo "000")
echo "    https://${SERVER_NAME}/healthz -> ${CODE}"
[[ "${CODE}" == "200" ]] && echo "==> Deployed." || echo "!! Not healthy; check: journalctl -u flexappeal-web -n 50"
