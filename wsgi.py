"""Gunicorn entrypoint:  gunicorn --config deploy/gunicorn.conf.py wsgi:app

See deploy/flexappeal-web.service.
"""

from flexappeal.webapp import create_app

app = create_app()
