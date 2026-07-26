/* Analysis tab: draw the Plotly panels, drive the Mol* viewer, and link them.
 *
 * Vanilla JS, no framework. Both libraries are vendored, never fetched from a
 * CDN, so this page works with no third-party network access at all.
 */

(function () {
  'use strict';

  /* One pass through the trajectory takes this long, however many frames it
     holds. Mol*'s own parameter is durationInS and it clamps to 1-120, so this
     is in seconds, not milliseconds. Fifty is slow enough to follow a loop
     opening and closing and short enough that the whole run has been seen
     before anyone scrolls past. */
  const PLAYBACK_SECONDS = 50;

  /* ------------------------------------------------------- drop zone */

  function initDropZone() {
    const zone = document.getElementById('drop-zone');
    const input = document.getElementById('fxa-input');
    const hint = document.getElementById('drop-hint');
    if (!zone || !input) return;

    input.addEventListener('change', function () {
      if (input.files.length) hint.textContent = input.files[0].name;
    });

    ['dragenter', 'dragover'].forEach(function (event) {
      zone.addEventListener(event, function (e) {
        e.preventDefault();
        zone.classList.add('is-over');
      });
    });
    ['dragleave', 'drop'].forEach(function (event) {
      zone.addEventListener(event, function (e) {
        e.preventDefault();
        zone.classList.remove('is-over');
      });
    });
    zone.addEventListener('drop', function (e) {
      if (!e.dataTransfer || !e.dataTransfer.files.length) return;
      input.files = e.dataTransfer.files;
      hint.textContent = e.dataTransfer.files[0].name;
    });
  }

  /* ---------------------------------------------------------- plots */

  const PLOT_CONFIG = {
    responsive: true,
    displaylogo: false,
    // lasso and box-select do nothing useful on a time series and only crowd
    // the mode bar; the rest stays for zoom, pan and PNG export.
    modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d'],
    toImageButtonOptions: { format: 'png', scale: 2 }
  };

  const drawn = [];

  function initPlots() {
    if (typeof Plotly === 'undefined') return;

    document.querySelectorAll('.md-plot[data-figure]').forEach(function (node) {
      let figure;
      try {
        figure = JSON.parse(node.dataset.figure);
      } catch (err) {
        node.innerHTML = '<p class="md-plot-error">This panel could not be drawn.</p>';
        return;
      }
      Plotly.newPlot(node, figure.data, figure.layout, PLOT_CONFIG);
      drawn.push(node);
    });

    window.addEventListener('resize', function () {
      drawn.forEach(function (node) { Plotly.Plots.resize(node); });
    });
  }

  /* --------------------------------------------------------- viewer */

  /* Frame linking: the viewer and the time-series plots share one x axis
   * (simulated time), so moving in either should move the other. Kept
   * deliberately simple -- a vertical line on each time plot, updated as the
   * viewer's frame changes. */

  function timePlots() {
    return drawn.filter(function (node) {
      const layout = node.layout;
      return layout && layout.xaxis && layout.xaxis.title &&
        String(layout.xaxis.title.text || '').indexOf('Time') === 0;
    });
  }

  function markFrame(ns) {
    timePlots().forEach(function (node) {
      Plotly.relayout(node, {
        shapes: [{
          type: 'line', x0: ns, x1: ns, y0: 0, y1: 1, yref: 'paper',
          line: { color: '#6b7c93', width: 1, dash: 'dot' }
        }]
      });
    });
  }

  /* Start the trajectory looping as soon as it is loaded.

     Mol* registers its animations by name and does not export the model-index
     one on the global bundle, so it is looked up rather than imported -- which
     also means a version that renames or drops it degrades to a static
     structure with working manual controls, instead of throwing inside the
     load chain and losing the viewer entirely.

     Skipped for anyone who has asked their system for reduced motion. A
     3D scene that starts moving on its own is exactly what that setting is
     about, and the playback controls are still there for them. */
  function autoplay(viewer) {
    if (window.matchMedia
        && window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
      return false;
    }
    try {
      const manager = viewer.plugin.managers.animation;
      const animation = manager.animations.find(function (a) {
        return a.name === 'built-in.animate-model-index';
      });
      if (!animation) return false;
      manager.play(animation, {
        mode: { name: 'loop', params: { direction: 'forward' } },
        duration: { name: 'fixed', params: { durationInS: PLAYBACK_SECONDS } }
      });
      return true;
    } catch (err) {
      return false;      /* a static structure is a fine outcome; a broken page is not */
    }
  }

  function initViewer() {
    const host = document.getElementById('molstar-viewer');
    const status = document.getElementById('viewer-status');
    if (!host || typeof molstar === 'undefined') return;

    const structureUrl = host.dataset.structure;
    const trajectoryUrl = host.dataset.trajectory;
    const duration = parseFloat(host.dataset.duration || '0');
    const frames = parseInt(host.dataset.frames || '0', 10);

    function say(message, isError) {
      if (!status) return;
      status.textContent = message;
      status.classList.toggle('is-error', Boolean(isError));
    }

    molstar.Viewer.create(host, {
      layoutIsExpanded: false,
      layoutShowControls: false,
      layoutShowSequence: true,
      layoutShowLog: false,
      layoutShowLeftPanel: false,
      viewportShowExpand: true,
      viewportShowSelectionMode: false,
      viewportShowAnimation: true,
      pdbProvider: 'rcsb',
      emdbProvider: 'rcsb'
    }).then(function (viewer) {
      /* Topology plus coordinates: Mol* reads the PDB for connectivity and the
       * XTC for the frames, which is exactly how the .fxa stores them. */
      return viewer.loadTrajectory({
        model: { kind: 'model-url', url: structureUrl, format: 'pdb' },
        coordinates: { kind: 'coordinates-url', url: trajectoryUrl,
                       format: 'xtc', isBinary: true },
        preset: 'default'
      }).then(function () {
        const playing = frames > 1 && autoplay(viewer);
        say(frames
          ? frames + ' frames'
              + (playing
                  ? ', playing on a ' + PLAYBACK_SECONDS + ' second loop — the controls under the viewport pause it.'
                  : ' loaded — use the playback controls under the viewport.')
          : 'Structure loaded.');

        /* Poll the animation frame rather than subscribing to Mol*'s internal
         * state: the public surface for "which frame is showing" is not stable
         * across versions, and a 250 ms poll is imperceptible for this. */
        if (duration > 0 && frames > 1) {
          let last = -1;
          setInterval(function () {
            try {
              const cell = viewer.plugin.state.data.select(
                molstar.StateSelection.Generators.rootsOfType(
                  molstar.PluginStateObject.Molecule.Model))[0];
              const index = cell && cell.obj && cell.obj.data
                ? (cell.obj.data.modelNum - 1) : -1;
              if (index >= 0 && index !== last) {
                last = index;
                markFrame((index / (frames - 1)) * duration);
              }
            } catch (err) { /* viewer not ready, or a version without this shape */ }
          }, 250);
        }
      });
    }).catch(function (err) {
      say('The 3D viewer could not load this structure (' + err + '). '
        + 'The plots below are unaffected.', true);
    });
  }

  /* ----------------------------------------------------------- boot */

  document.addEventListener('DOMContentLoaded', function () {
    initDropZone();
    initPlots();
    initViewer();
  });
})();

/* ------------------------------------------------------- re-analysis -----
 * Fire a bounded server-side job, then poll. The work happens in a detached
 * subprocess, so this is a status poll rather than a long-held request.
 */

(function () {
  'use strict';

  const POLL_MS = 1500;
  const PLOT_CONFIG = {
    responsive: true, displaylogo: false,
    modeBarButtonsToRemove: ['lasso2d', 'select2d', 'autoScale2d']
  };

  function init() {
    const card = document.getElementById('reanalyse-card');
    if (!card) return;

    const button = document.getElementById('re-run');
    const status = document.getElementById('re-status');
    const output = document.getElementById('re-results');
    const selection = document.getElementById('re-selection');
    const selectionB = document.getElementById('re-selection-b');
    const selectionBField = document.getElementById('re-selection-b-field');

    function chosen() {
      return Array.from(document.querySelectorAll('input[name="re-metric"]:checked'))
        .map(function (el) { return el.value; });
    }

    /* The second selection only means anything for a distance. */
    function syncDistanceField() {
      selectionBField.hidden = chosen().indexOf('distance') === -1;
    }
    document.querySelectorAll('input[name="re-metric"]').forEach(function (el) {
      el.addEventListener('change', syncDistanceField);
    });
    syncDistanceField();

    function say(message, kind) {
      status.textContent = message;
      status.className = 'md-reanalyse-status' + (kind ? ' is-' + kind : '');
    }

    function draw(panels) {
      output.innerHTML = '';
      if (!panels || !panels.length) {
        say('That produced nothing to plot.', 'error');
        return;
      }
      panels.forEach(function (panel) {
        const wrap = document.createElement('div');
        wrap.className = 'md-panel md-reanalyse-panel';
        const heading = document.createElement('h4');
        heading.textContent = panel.title;
        const target = document.createElement('div');
        target.className = 'md-plot';
        target.id = 'plot-' + panel.id;
        wrap.appendChild(heading);
        wrap.appendChild(target);
        output.appendChild(wrap);
        Plotly.newPlot(target, panel.figure.data, panel.figure.layout, PLOT_CONFIG);
      });
    }

    let timer = null;

    function poll() {
      fetch(card.dataset.status)
        .then(function (r) { return r.json(); })
        .then(function (data) {
          if (data.status === 'running') {
            say('Working' + (data.seconds ? ' (' + data.seconds + 's)' : '') + '…');
            timer = setTimeout(poll, POLL_MS);
            return;
          }
          button.disabled = false;
          if (data.status === 'ready') {
            const m = data.metrics || {};
            say('Done — ' + m.n_atoms + ' atoms, ' + m.n_frames + ' frames.', 'ok');
            draw(data.figures);
          } else {
            say(data.message || 'That did not work.', 'error');
          }
        })
        .catch(function () {
          button.disabled = false;
          say('Lost contact with the server.', 'error');
        });
    }

    button.addEventListener('click', function () {
      const metrics = chosen();
      if (!metrics.length) {
        say('Choose at least one metric.', 'error');
        return;
      }
      clearTimeout(timer);
      button.disabled = true;
      say('Starting…');
      output.innerHTML = '';

      fetch(card.dataset.endpoint, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          metrics: metrics,
          selection: selection.value,
          selection_b: selectionB ? selectionB.value : ''
        })
      })
        .then(function (r) { return r.json().then(function (d) { return [r.status, d]; }); })
        .then(function (pair) {
          const code = pair[0], data = pair[1];
          if (code === 202) { timer = setTimeout(poll, POLL_MS); return; }
          button.disabled = false;
          say(data.message || 'That request was refused.', 'error');
        })
        .catch(function () {
          button.disabled = false;
          say('Could not reach the server.', 'error');
        });
    });
  }

  document.addEventListener('DOMContentLoaded', init);
})();
