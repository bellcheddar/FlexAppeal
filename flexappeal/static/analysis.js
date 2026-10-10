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

  /* Neon presentation, matching DeepFlex's pocket viewer so the same complex
   * looks like the same complex in both places.
   *
   * The constants are DeepFlex's, which got them from BINMAN, and every one of
   * them is tuning somebody already paid for:
   *
   *   - carbons are pinned to one vivid colour rather than coloured by chain. A
   *     chain-coloured ligand can come out the same hue as its own shell and
   *     vanish into the glow; the colour a ligand gets should not depend on
   *     which chain it happens to sit beside. Oxygen, nitrogen and sulfur keep
   *     their element colours, so the molecule still reads chemically.
   *   - sticks at 0.34 rather than the default. Thinner and the glow sits over
   *     them, and the ligand reads as a cloud of unconnected dots: the bonds are
   *     there and simply cannot be seen.
   *   - the halo is a shell, not a fog. The probe radius inflates the surface
   *     clear of the atoms, and its emissive is below the sticks' own so the
   *     bloom pass cannot wash them out.
   *   - bloom runs in `emissive` mode, so only geometry carrying emissive > 0
   *     blooms and the protein behind the ligand stays readable.
   *
   * All of it is guarded. A viewer that cannot apply the styling shows the
   * default presentation, which is a worse-looking panel rather than no panel.
   */
  const NEON = {
    background: 0x07080F,
    neon: 0x00F5D4,
    carbon: 0xFF2D9E,
    stick: { sizeFactor: 0.34, aspectRatio: 1, emissive: 0.6 },
    cloud: { alpha: 0.09, emissive: 0.7, probeRadius: 2.4, resolution: 0.5 },
    bloom: { strength: 1.3, radius: 0.6, threshold: 0,
             mode: 'emissive', transparency: true }
  };

  function ligandComponents(plugin) {
    const out = [];
    try {
      const current = plugin.managers.structure.hierarchy.current;
      (current.structures || []).forEach(function (ref) {
        (ref.components || []).forEach(function (component) {
          const key = String(component.key || '').toLowerCase();
          let label = '';
          try {
            label = String(component.cell.obj.label || '').toLowerCase();
          } catch (err) { /* unlabelled is still testable by key */ }
          if (key.indexOf('ligand') >= 0 || label.indexOf('ligand') >= 0) {
            out.push(component);
          }
        });
      });
    } catch (err) { /* no hierarchy yet */ }
    return out;
  }

  function applyNeon(viewer) {
    const plugin = viewer.plugin;
    if (!plugin || !plugin.canvas3d) return;

    try {
      plugin.canvas3d.setProps({
        renderer: { backgroundColor: NEON.background },
        postprocessing: { bloom: { name: 'on', params: NEON.bloom } }
      });
    } catch (err) { /* a light viewport is not worth failing the mount over */ }

    const components = ligandComponents(plugin);
    if (!components.length) return;

    /* Retune the ligand's existing sticks in place. Only a ball-and-stick is
       touched: renaming a type while keeping another type's params is how a
       representation ends up with nonsense settings. */
    try {
      const update = plugin.build();
      components.forEach(function (component) {
        (component.representations || []).forEach(function (rep) {
          if (!rep.cell) return;
          update.to(rep.cell).update(function (old) {
            if (!old.type || old.type.name !== 'ball-and-stick') return;
            old.type.params = old.type.params || {};
            old.type.params.sizeFactor = NEON.stick.sizeFactor;
            old.type.params.aspectRatio = NEON.stick.aspectRatio;
            old.type.params.emissive = NEON.stick.emissive;
            /* No hydrogens on the ligand. They treble the stick count on a
               small molecule and bury the heavy-atom skeleton that the pose,
               the contacts and the 2D diagram beside it are all about. The
               diagram strips them too (measure.depiction uses RemoveHs), so
               the two views of the ligand now show the same atoms. */
            old.type.params.ignoreHydrogens = true;
            old.colorTheme = {
              name: 'element-symbol',
              params: {
                carbonColor: { name: 'uniform', params: { value: NEON.carbon } }
              }
            };
          });
        });
      });
      update.commit();
    } catch (err) { /* default sticks are a fine outcome */ }

    /* The shell. Added as a second representation rather than replacing the
       sticks, because it is the glow around them. */
    try {
      const builder = plugin.builders.structure.representation;
      components.forEach(function (component) {
        builder.addRepresentation(component.cell, {
          type: 'molecular-surface',
          typeParams: {
            alpha: NEON.cloud.alpha, emissive: NEON.cloud.emissive,
            probeRadius: NEON.cloud.probeRadius,
            resolution: NEON.cloud.resolution, ignoreLight: true
          },
          color: 'uniform',
          colorParams: { value: NEON.neon }
        }).catch(function () {
          /* Surface meshing fails on a very small or disordered ligand.
             Inflated spheres are the same idea with no mesh to build. */
          return builder.addRepresentation(component.cell, {
            type: 'spacefill',
            typeParams: {
              alpha: NEON.cloud.alpha, emissive: NEON.cloud.emissive,
              sizeFactor: 1.7, ignoreLight: true
            },
            color: 'uniform',
            colorParams: { value: NEON.neon }
          }).catch(function () { /* no halo; the sticks still glow */ });
        });
      });
    } catch (err) { /* no halo; the sticks still glow */ }
  }

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

  /* The cursor is appended to whatever shapes a panel already has, not
     substituted for them. Several binding panels carry their own: the shaded
     equilibration head, the 5 A drift gate, the quarter-retained contact
     threshold. Replacing the list wiped all of those the first time the slider
     moved, so a panel lost its own annotations as soon as anybody scrubbed it. */
  const baseShapes = new WeakMap();

  function markFrame(ns) {
    timePlots().forEach(function (node) {
      if (!baseShapes.has(node)) {
        const existing = (node.layout && node.layout.shapes) || [];
        baseShapes.set(node, existing.map(function (s) { return s; }));
      }
      const shapes = baseShapes.get(node).concat([{
        type: 'line', x0: ns, x1: ns, y0: 0, y1: 1, yref: 'paper',
        line: { color: '#1a1a2e', width: 1.5 }
      }]);
      Plotly.relayout(node, { shapes: shapes });
    });
  }

  /* ------------------------------------------------------- frame slider */

  /* Step a multi-model structure WITHOUT reloading it.

     The frame is a parameter on whichever state cell carries a modelIndex, found
     by that rather than by transformer name: the name differs between the
     trajectory and the model-from-trajectory transforms, and this only needs the
     thing with the knob on it. Updating it costs a parameter commit and never
     touches the camera, so zoom and orientation stay where the user left them.
     Loading a frame as a structure instead would make Mol* tear the scene down
     and re-frame it, which is what makes a scrub flash and jump. */
  function setFrame(viewer, index) {
    try {
      const plugin = viewer.plugin;
      if (!plugin || !plugin.state || !plugin.build) return false;
      const cells = [];
      plugin.state.data.cells.forEach(function (cell) {
        const params = cell.transform && cell.transform.params;
        if (params && typeof params.modelIndex === 'number') cells.push(cell);
      });
      if (!cells.length) return false;
      const update = plugin.build();
      cells.forEach(function (cell) {
        update.to(cell.transform.ref).update(function (old) {
          old.modelIndex = index;
        });
      });
      update.commit();
      return true;
    } catch (err) {
      return false;
    }
  }

  /* Wire the slider, the readout, the play button and the plot cursor together.

     Returns true when the slider is driving playback, so the caller knows not
     to also start Mol*'s own animation loop: two things stepping the same model
     index fight, and the slider's position then disagrees with what is drawn. */
  function wireSlider(viewer, frames, duration) {
    const slider = document.getElementById('scrub');
    if (!slider || frames < 2) return false;
    const readout = document.getElementById('scrub-readout');
    const play = document.getElementById('scrub-play');

    function show() {
      const index = parseInt(slider.value, 10);
      const ns = duration > 0 ? (duration * index) / (frames - 1) : index;
      if (readout) {
        readout.textContent = duration > 0
          ? ns.toFixed(duration < 10 ? 3 : 1) + ' ns'
          : 'frame ' + index;
      }
      markFrame(ns);
      /* No debounce: a model switch is a parameter update rather than a load,
         so dragging can drive it directly. */
      setFrame(viewer, index);
    }
    slider.addEventListener('input', show);

    let timer = null;
    /* Set only by the person, never by the page. Without it, anything that
       restarts playback would override a deliberate pause. */
    let paused = window.matchMedia
      && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

    function stop(byUser) {
      if (byUser) paused = true;
      if (timer) { clearInterval(timer); timer = null; }
      if (play) {
        play.textContent = 'Play';
        play.setAttribute('aria-label', 'Play the trajectory');
      }
    }
    function start(byUser) {
      if (byUser) paused = false;
      if (timer || paused) return;
      if (play) {
        play.textContent = 'Pause';
        play.setAttribute('aria-label', 'Pause the trajectory');
      }
      const every = Math.max(40, Math.round((PLAYBACK_SECONDS * 1000) / frames));
      timer = setInterval(function () {
        slider.value = String((parseInt(slider.value, 10) + 1) % frames);
        show();
      }, every);
    }
    if (play) {
      play.addEventListener('click', function () {
        if (timer) { stop(true); } else { start(true); }
      });
    }
    /* Dragging pauses: the person has taken the wheel, and a loop that keeps
       advancing underneath them is fighting the hand on the control. */
    slider.addEventListener('pointerdown', function () { stop(true); });

    show();
    start(false);
    return true;
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

  /* Mol* injects its own logo as a link with no text and no label, which
   * axe-core reports as a serious link-name violation. It is the library's
   * markup rather than ours, so it is labelled after the viewer builds rather
   * than fixed at the source. Named rather than hidden: it is a real link to
   * the project and crediting it is the point of it being there. */
  function nameTheMolstarLogo(host) {
    /* Retried, because Mol* builds its viewport controls after the promise
     * that loads the structure resolves: labelling it once on that resolution
     * found no logo and left the violation in place. Bounded rather than a
     * MutationObserver left running for the life of the page. */
    let attempts = 0;
    (function attempt() {
      const logo = host.querySelector('.msp-logo');
      if (logo) {
        if (!logo.getAttribute('aria-label')) {
          logo.setAttribute('aria-label', 'Mol* molecular viewer, opens molstar.org');
        }
        if (!logo.getAttribute('title')) {
          logo.setAttribute('title', 'Mol* molecular viewer');
        }
        return;
      }
      if (++attempts < 40) setTimeout(attempt, 100);
    })();
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
        /* The slider is the one thing stepping the model, where it exists.
         *
         * It replaces both halves of what was here: Mol*'s own model-index
         * animation, and a 250 ms poll that watched the frame change so the plot
         * cursor could follow. Two things driving one model index fight, and the
         * slider's position then disagrees with what is drawn; and a poll could
         * only ever read the frame, never set it, so there was no way to go to a
         * frame a plot made you curious about. Driving it from one control makes
         * the structure and the cursor the same fact.
         *
         * Mol*'s animation stays as the fallback for a page without the slider,
         * which is any single-frame structure. */
        nameTheMolstarLogo(host);
        applyNeon(viewer);
        const scrubbing = wireSlider(viewer, frames, duration);
        const playing = scrubbing || (frames > 1 && autoplay(viewer));
        say(frames
          ? frames + ' frames'
              + (scrubbing
                  ? ', on a ' + PLAYBACK_SECONDS + ' second loop. Drag the slider to a frame and every time plot marks it.'
                  : (playing
                      ? ', playing on a ' + PLAYBACK_SECONDS + ' second loop — the controls under the viewport pause it.'
                      : ' loaded — use the playback controls under the viewport.'))
          : 'Structure loaded.');
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
