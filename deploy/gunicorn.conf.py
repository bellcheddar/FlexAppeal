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

# Import the application once in the master and fork, rather than importing it
# again in every worker on its first request.
#
# This was added while chasing a 12.7 second first request, on the theory that
# it was numpy and MDTraj being imported inside whichever request arrived
# first. **That theory was wrong**: the first request was still 12.8 seconds
# with preloading on, and the cost turned out to be the binding-site
# measurement, which is cached on disk now (4bv). Timed on the box afterwards,
# importing the whole application takes 0.21 s.
#
# Kept anyway, for the reason that survives: on a box with 3.8 GB shared
# between five applications, preloading puts the modules in pages the workers
# share after the fork rather than in each worker's own heap. It is a memory
# measure, not a latency one, and the comment says so rather than leaving a
# disproved explanation in the file for the next reader to believe.
#
# The usual caution about preloading is code that opens a file or a socket at
# import time and then has it shared across forks. This app opens the database
# per request and `web_scratch` per session, so there is nothing to share.
preload_app = True

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
