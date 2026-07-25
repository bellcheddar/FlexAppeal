/* Analysis tab: draw the Plotly panels, drive the Mol* viewer, and link them.
 *
 * Vanilla JS, no framework. Both libraries are vendored, never fetched from a
 * CDN, so this page works with no third-party network access at all.
 */

(function () {
  'use strict';

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
        say(frames
          ? frames + ' frames loaded — use the playback controls under the viewport.'
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
