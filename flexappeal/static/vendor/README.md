# Vendored libraries

Committed rather than fetched from a CDN, so the Analysis tab works with no
third-party network access at all and cannot break when someone else's CDN does.
Both are MIT licensed; their notices are reproduced in the repository `LICENSE`.

| File | Version | Licence | Upstream |
|---|---|---|---|
| `plotly.min.js` | 2.35.2 | MIT | https://github.com/plotly/plotly.js |
| `molstar.js`, `molstar.css` | 5.11.0 | MIT | https://github.com/molstar/molstar |

`molstar.version` records the exact version, so a refresh is a deliberate act
with a diff rather than a silent drift.

To update:

```bash
V=$(curl -s https://registry.npmjs.org/molstar/latest | python3 -c "import json,sys; print(json.load(sys.stdin)['version'])")
curl -sS -o molstar.js  "https://cdn.jsdelivr.net/npm/molstar@${V}/build/viewer/molstar.js"
curl -sS -o molstar.css "https://cdn.jsdelivr.net/npm/molstar@${V}/build/viewer/molstar.css"
echo "$V" > molstar.version
```
