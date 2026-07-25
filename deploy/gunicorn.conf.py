"""Gunicorn config for the FlexAppeal web service.

Referenced by deploy/flexappeal-web.service.
"""
import os

# 8004 on the shared droplet: 8000 AlphaFraud, 8001 chem_sage, 8002 chatPDB,
# 8003 BoltzMaker. Parameterised rather than hardcoded so the port lives in one
# place (.env) instead of being duplicated here and in the nginx proxy_pass.
bind = os.environ.get("BIND_ADDR", "127.0.0.1:8004")

# Two, not three. Each worker can hold a parsed .fxa plus its numpy arrays, and
# this box has 3.8 GB shared with four other applications.
workers = int(os.environ.get("WEB_WORKERS", "2"))
worker_class = "sync"

# MUST stay in sync with `proxy_read_timeout` in nginx-flexappeal.conf.
# Uploading and parsing a 250 MB results file is the slow path; the re-analysis
# route does not block a worker at all, because it forks.
timeout = 300
graceful_timeout = 30
keepalive = 5

# Log to stdout/stderr so journald captures everything (journalctl -u flexappeal-web).
accesslog = "-"
errorlog = "-"
loglevel = os.environ.get("LOG_LEVEL", "info")
proc_name = "flexappeal-web"
