#!/usr/bin/env bash
# One-time (and idempotent) provisioning of FlexAppeal on the droplet.
# Run as root ON THE DROPLET:  bash /opt/flexappeal/deploy/provision.sh
#
# Safe to re-run: every step either creates or updates in place. Re-running is
# in fact the intended way to apply a change to the nginx site or a systemd unit.
#
# Before the first run, flexappeal.mdeller.com must resolve to this droplet.
# certbot proves domain control over HTTP, so without the DNS record the TLS
# step fails and the site is left serving plain HTTP.
set -euo pipefail

APP=flexappeal
APP_PATH=/opt/${APP}
SERVER_NAME="${SERVER_NAME:-flexappeal.mdeller.com}"
BIND_ADDR="${BIND_ADDR:-127.0.0.1:8004}"
CERTBOT_EMAIL="${CERTBOT_EMAIL:-marc@marcdeller.com}"

if [[ $EUID -ne 0 ]]; then echo "Run this as root."; exit 1; fi

echo "==> Packages"
apt-get update -qq
apt-get install -y -qq python3-venv python3-pip python3-dev build-essential \
    nginx certbot python3-certbot-nginx rsync

echo "==> Service user"
if ! id -u "${APP}" >/dev/null 2>&1; then
    useradd --system --shell /usr/sbin/nologin --home "${APP_PATH}" "${APP}"
fi

echo "==> Ownership"
# Before the venv, not after. deploy.sh rsyncs as root, so everything under
# /opt/flexappeal arrives root-owned; building the virtual environment as the
# service user then fails with a bare
#     Error: [Errno 13] Permission denied: '/opt/flexappeal/.venv'
# which says nothing about ownership being the cause. Chowned again at the end
# to catch anything the steps below create as root.
mkdir -p "${APP_PATH}"
chown -R "${APP}:${APP}" "${APP_PATH}"

echo "==> Virtual environment"
if [[ ! -x "${APP_PATH}/.venv/bin/python" ]]; then
    sudo -u "${APP}" python3 -m venv "${APP_PATH}/.venv"
fi
sudo -u "${APP}" env PIP_NO_CACHE_DIR=1 "${APP_PATH}/.venv/bin/pip" install --quiet --upgrade pip
sudo -u "${APP}" env PIP_NO_CACHE_DIR=1 "${APP_PATH}/.venv/bin/pip" install --quiet -r "${APP_PATH}/requirements.txt"

echo "==> Scratch directory"
mkdir -p "${APP_PATH}/web_scratch"
chown -R "${APP}:${APP}" "${APP_PATH}/web_scratch"

echo "==> systemd units"
install -m 644 "${APP_PATH}/deploy/${APP}-web.service" /etc/systemd/system/
install -m 644 "${APP_PATH}/deploy/${APP}-scratch-clean.service" /etc/systemd/system/
install -m 644 "${APP_PATH}/deploy/${APP}-scratch-clean.timer" /etc/systemd/system/
systemctl daemon-reload
systemctl enable --now "${APP}-web.service"
systemctl enable --now "${APP}-scratch-clean.timer"

echo "==> nginx"
# The rate-limit zone has to be declared at http scope, which a sites-available
# file is not -- hence the separate conf.d snippet.
install -m 644 "${APP_PATH}/deploy/nginx-${APP}-limits.conf" /etc/nginx/conf.d/
sed -e "s|__SERVER_NAME__|${SERVER_NAME}|g" \
    -e "s|__BIND_ADDR__|${BIND_ADDR}|g" \
    "${APP_PATH}/deploy/nginx-${APP}.conf" > "/etc/nginx/sites-available/${APP}"
ln -sf "/etc/nginx/sites-available/${APP}" "/etc/nginx/sites-enabled/${APP}"
nginx -t
systemctl reload nginx

echo "==> TLS"
# certbot runs unconditionally on every provision, never skipped because a
# certificate already exists. The vhost above was just re-templated from the
# plain-HTTP source, so its TLS block is gone; skipping certbot here leaves the
# site with no SSL server block at all, and nginx then answers HTTPS for this
# name from whichever other vhost happens to hold a certificate. That is not
# hypothetical -- it is how boltzmaker.mdeller.com once served AlphaFraud's.
#
# Retried because certbot's own renewal timer can be holding the lock.
if ! host "${SERVER_NAME}" >/dev/null 2>&1; then
    echo "!! ${SERVER_NAME} does not resolve. Add an A record pointing at this"
    echo "!! droplet, then re-run this script. Skipping TLS; the site is HTTP-only."
else
    certbot_ok=0
    for attempt in 1 2 3; do
        if certbot --nginx -d "${SERVER_NAME}" --non-interactive --agree-tos \
                -m "${CERTBOT_EMAIL}" --redirect; then
            certbot_ok=1; break
        fi
        echo "    certbot attempt ${attempt}/3 failed, retrying in 10s..."
        sleep 10
    done
    [[ $certbot_ok -eq 1 ]] || echo "!! certbot failed three times; the site is HTTP-only."

    # certbot does not enable HTTP/2 on nginx 1.24, so patch it in. Idempotent.
    #
    # Both listen directives have to be patched, not just one. In nginx 1.24
    # `ssl` and `http2` are protocol options on the listening socket rather than
    # per-server settings, so if one vhost declares 0.0.0.0:443 differently from
    # the others nginx warns "protocol options redefined for 0.0.0.0:443" and
    # honours whichever server block was parsed first. The first version of this
    # regex put the space inside the optional group:
    #
    #     listen(\s+\[::\]:)?443 ssl;
    #
    # which requires "listen443" when the group is absent, so the IPv6 line was
    # patched and the IPv4 line silently was not. nginx 1.25.1 replaces all of
    # this with a separate `http2 on;` directive; this box is on 1.24.
    python3 - "$SERVER_NAME" <<'PATCH'
import re, sys, pathlib
path = pathlib.Path("/etc/nginx/sites-available/flexappeal")
text = path.read_text()
patched = re.sub(r"listen(\s+\[::\]:)?\s*443 ssl;",
                 lambda m: m.group(0)[:-1] + " http2;", text)
if patched != text:
    path.write_text(patched)
    print("    HTTP/2 enabled")
PATCH
    nginx -t && systemctl reload nginx
fi

echo "==> Ownership"
chown -R "${APP}:${APP}" "${APP_PATH}"

echo "==> Status"
systemctl --no-pager --lines=3 status "${APP}-web.service" || true
echo
echo "Provisioned. Check:  curl -sI https://${SERVER_NAME}/healthz"
