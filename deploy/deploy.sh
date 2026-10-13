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

# Measured, not claimed. These three are what the explorer reads, and a deploy
# that is suddenly much bigger than the last one is worth seeing before it runs
# rather than wondering about afterwards.
SERVED=$(du -ch web_runs data/refstates.db 2>/dev/null | tail -1 | cut -f1)
echo "==> Run library and explorer data to send: ${SERVED:-none found}"

echo "==> Syncing to ${DROPLET_SSH}:${DROPLET_PATH}"
# The vendored JavaScript is ~10 MB and rsync will only send it once, but it is
# genuinely needed on the server -- nginx serves it from disk. Not excluded.
#
# .pixi/ is the local development environment (several GB of OpenMM); the droplet
# neither has nor needs it. web_scratch/ holds live user sessions.
#
# ${arr[@]+"${arr[@]}"} expands to nothing when empty without tripping `set -u`,
# which macOS's bash 3.2 requires.
#
# --chmod, because -a preserves the source's permissions and a repo checked out
# on a Mac can easily be 0700. That lands /opt/flexappeal itself as drwx------
# owned by the syncing user, and the service user then cannot traverse its own
# directory: the first symptom is pip "Permission denied" two lines later, which
# names neither the directory nor the mode.
# `protect web_runs/**`, because that directory has two writers. The four
# reference-state runs are packed here and travel with the code, and anything a
# visitor keeps from the Analysis tab is written there by the server. With a
# bare --delete the next deploy would silently remove every run a user had kept,
# and the directory was empty on both sides until now, so nothing would have
# noticed. Protect stops the receiver deleting what the sender does not have; it
# does not stop the sender replacing a file it does have.
rsync -az --delete --filter='protect web_runs/**' --chmod=D755,F644 \
    ${SSH_OPTS[@]+"${SSH_OPTS[@]}"} \
    --exclude '.venv/' --exclude '.venv-dev/' --exclude '.pixi/' \
    --exclude '__pycache__/' --exclude '*.pyc' --exclude '.git/' \
    --exclude '.env' --exclude 'web_scratch/' --exclude '.pytest_cache/' \
    --exclude 'data/md/' --exclude 'data/cache/' --exclude 'data/refs/' \
    --exclude 'examples/*/run/' --exclude 'examples/*/terminal.log' \
    --exclude 'examples/*/analyse.log' \
    ./ "${DROPLET_SSH}:${DROPLET_PATH}/"

echo "==> Installing dependencies and restarting on the droplet"
SSH_CMD=(ssh)
[[ -n "$SSH_KEY" ]] && SSH_CMD=(ssh -i "${SSH_KEY/#\~/$HOME}")
# `<<'REMOTE'`, quoted, and the path passed as an argument instead of
# interpolated. Unquoted, the local shell expands this heredoc before sending
# it, so every backtick in the comments below was command substitution run on
# the Mac: the first real run of this script printed `data/md/: is a directory`
# and `data/refstates.db: Permission denied`, sent a remote script with those
# substitutions blanked out, and died without installing, chowning or
# restarting anything. The rsync above had already succeeded, so the deploy
# looked like it had worked while the service went on running the old code.
"${SSH_CMD[@]}" "$DROPLET_SSH" bash -s -- "$DROPLET_PATH" <<'REMOTE'
set -euo pipefail
DROPLET_PATH="$1"
cd "${DROPLET_PATH}"

if [[ ! -x .venv/bin/python ]]; then
    echo "No virtual environment yet -- run deploy/provision.sh as root first."
    exit 1
fi

# The working data stays here. `data/md/` is the raw simulation output and it
# is already 565 MB with one row a sixth of the way through its 50 ns; four
# finished rows are several gigabytes of trajectory, checkpoint and system XML.
# `data/cache/` and `data/refs/` are downloaded structures, which exist so a
# rebuild does not hit KLIFS and GPCRdb again. The web app reads none of the
# three: it reads `data/refstates.db`, `data/analysis/` and `data/systems/`,
# which are not excluded and whose size is printed below rather than asserted
# here. It used to say "under 5 MB"; it passed 5 MB while the second row was
# still a third of the way through, because `data/analysis/` carries a viewer
# bundle per row and those grow with the trajectory. Section 6.3 caps a
# per-system web bundle at 150 MB, so there is room, but the number in a
# comment goes stale and the number from `du` cannot. The example run's
# untracked output is the same story, 183 MB of trajectory and terminal log
# that nothing serves.
#
# Explicit excludes rather than `--filter=':- .gitignore'`, which would be
# tidier and would be wrong: `data/refstates.db` is itself gitignored, so the
# filter would quietly ship a droplet with no explorer database. The app treats
# an absent database as an empty section rather than an error, by design, so
# the page would come up looking fine with the reference-states tab showing
# nothing at all.

# Ownership before pip, not after it. rsync runs as root, so every file it wrote
# is root-owned, and the service user then hits PermissionError at RUNTIME rather
# than at deploy time -- a confusing way to find out. Chown on every deploy, not
# only at provisioning.
#
# The ordering is the part that matters. With this after the install, a deploy
# that lands the tree unreadable fails in pip and exits before the step that
# would have fixed it, so every subsequent deploy fails the same way until
# someone chowns by hand. The install needs a traversable directory as much as
# the running service does.
#
# The venv and live scratch sessions are pruned: nothing in this sync touched
# them, and a blanket chown over a session an analysis has open is asking for
# trouble.
sudo find "${DROPLET_PATH}" \
    -path "${DROPLET_PATH}/.venv" -prune -o \
    -path "${DROPLET_PATH}/web_scratch" -prune -o \
    -exec chown flexappeal:flexappeal {} +
sudo chmod 755 "${DROPLET_PATH}"

sudo -u flexappeal env PIP_NO_CACHE_DIR=1 ./.venv/bin/pip install --quiet -r requirements.txt

sudo systemctl restart flexappeal-web.service
sleep 1
sudo systemctl --no-pager --lines=3 status flexappeal-web.service || true
REMOTE

# Verified on more than /healthz, because the old code is healthy too. The
# first real run of this script synced the files and then died before the
# restart, so the service went on serving a build from July: /healthz answered
# 200 the whole time while /refstates/ answered 404. A check that passes
# whether or not the deploy worked is not a check. Each route below is one a
# different part of this app owns, so a page tree that half loaded says so.
echo "==> Verifying"
SERVER_NAME="${SERVER_NAME:-flexappeal.mdeller.com}"
FAILED=0
for ROUTE in /healthz / /analysis /example /runs /runs/abl1_apo; do
    CODE=$(curl -s -o /dev/null -w '%{http_code}' \
           "https://${SERVER_NAME}${ROUTE}" || echo "000")
    printf '    %-22s -> %s\n' "${ROUTE}" "${CODE}"
    [[ "${CODE}" == "200" ]] || FAILED=1
done
# The library's own pages, fetched until they are warm.
#
# A run with a ligand measures its binding site the first time anybody opens
# it: 12 seconds for a 500-frame receptor, against 0.6 once the answer is
# cached. The cache is on disk and keyed on the results file, so that cost is
# paid once ever rather than once per worker per restart, and this loop is what
# pays it: a new or changed results file is measured here, by the deploy,
# instead of by whoever opens the page next.
#
# Fetched until two in a row come back quickly, because a request lands on
# whichever worker is free. The times are printed rather than hidden: a page
# that is still slow after the cache has been written is something to look at,
# and so is one that is slow on every deploy, because that would mean the cache
# is not being kept.
echo "==> Warming the run pages"
for RUN in $(curl -s "https://${SERVER_NAME}/runs" \
             | grep -oE 'href="/runs/[a-zA-Z0-9_-]+"' \
             | cut -d'"' -f2 | sort -u); do
    TIMES=""
    QUICK=0
    for _ in 1 2 3 4 5 6 7 8; do
        TAKEN=$(curl -s -o /dev/null -w '%{time_starttransfer}' \
                "https://${SERVER_NAME}${RUN}")
        TIMES="${TIMES} ${TAKEN}s"
        # Two consecutive fetches under two seconds means every worker that
        # answered has the measurement cached.
        if [[ "${TAKEN%%.*}" -lt 2 ]]; then
            QUICK=$((QUICK + 1))
            [[ "${QUICK}" -ge 2 ]] && break
        else
            QUICK=0
        fi
    done
    printf '    %-26s%s\n' "${RUN}" "${TIMES}"
done

# And one that must be gone. The fifth tab was removed on 13 October and a
# droplet still answering it is one running code from before that, which is
# the failure this whole section exists to catch: every route above answers
# 200 on the old build too.
GONE=$(curl -s -o /dev/null -w '%{http_code}' \
       "https://${SERVER_NAME}/refstates/" || echo "000")
printf '    %-22s -> %s (404 expected)\n' "/refstates/" "${GONE}"
[[ "${GONE}" == "404" ]] || FAILED=1
if [[ "${FAILED}" == "0" ]]; then
    echo "==> Deployed."
else
    echo "!! Not healthy; check: journalctl -u flexappeal-web -n 50"
    exit 1
fi
